#include "cimport_cache.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define CIMPORT_CACHE_SCRATCH ((size_t)1024 * 1024)
#define CIMPORT_CACHE_FLOAT_MISSING UINT32_C(0x7fc00001)

int cimport_cache_init(CImportCacheBuilder *b, int columns)
{
    memset(b, 0, sizeof(*b));
    if (columns <= 0)
        return -1;
    b->columns = columns;
    b->capacity = CIMPORT_CACHE_SCRATCH / (16 * (size_t)columns);
    if (b->capacity > 4096)
        b->capacity = 4096;
    if (!b->capacity)
        b->capacity = 1;
    b->values = ctools_safe_malloc3(b->capacity, (size_t)columns, sizeof(double));
    b->fields = ctools_safe_malloc3(b->capacity, (size_t)columns, sizeof(uint64_t));
    b->plan = ctools_safe_malloc2((size_t)columns, sizeof(*b->plan));
    b->numeric = malloc((size_t)columns);
    if (!b->values || !b->fields || !b->plan || !b->numeric)
    {
        cimport_cache_destroy(b);
        return -1;
    }
    return 0;
}

void cimport_cache_destroy(CImportCacheBuilder *b)
{
    free(b->values);
    free(b->fields);
    free(b->plan);
    free(b->numeric);
    memset(b, 0, sizeof(*b));
}

void cimport_cache_append(CImportCacheBuilder *b, const CImportFieldRef *fields, int count,
                          const double *values, uint64_t start, uint64_t end)
{
    if (!b->rows)
    {
        b->file_start = start;
        b->wide = false;
    }
    b->file_end = end;
    for (int col = 0; col < b->columns; col++)
    {
        size_t index = (size_t)col * b->capacity + b->rows;
        b->values[index] = values[col];
        uint64_t offset = col < count ? fields[col].offset - b->file_start : 0;
        if (offset > UINT32_MAX)
            b->wide = true;
        b->fields[index] = (offset << 32) | (col < count ? fields[col].length : 0);
    }
    b->rows++;
}

static size_t cimport_cache_plan(CImportCacheBuilder *b, double missing, bool as_float,
                                 bool as_double, bool inline_text)
{
    size_t bytes = sizeof(CImportCachedBlock) + (size_t)b->columns * sizeof(*b->plan);
    for (int col = 0; col < b->columns; col++)
    {
        CImportCachedColumn *p = &b->plan[col];
        memset(p, 0, sizeof(*p));
        size_t width = 0;
        if (!b->numeric[col])
        {
            const uint64_t *fields = b->fields + (size_t)col * b->capacity;
            p->length = (uint32_t)fields[0];
            bool fixed = true;
            for (size_t r = 1; r < b->rows; r++)
                if ((uint32_t)fields[r] != p->length)
                {
                    fixed = false;
                    break;
                }
            p->kind = fixed ? CIMPORT_CACHE_FIXED_TEXT : CIMPORT_CACHE_TEXT;
            width = fixed ? 4 : 8;
            if (inline_text)
            {
                uint32_t longest = 0;
                for (size_t r = 0; r < b->rows; r++)
                {
                    uint32_t length = (uint32_t)fields[r] & 0x7FFFFFFFu;
                    if (length > longest)
                        longest = length;
                }
                /* Raw width bounds the decoded width. Avoid copying large text
                 * sparsely filling a block; source references are cheaper there. */
                if (longest <= 64)
                {
                    p->kind = CIMPORT_CACHE_INLINE_TEXT;
                    width = p->length = longest ? longest + 1 : 0;
                }
            }
        }
        else
        {
            const double *values = b->values + (size_t)col * b->capacity;
            double low = DBL_MAX, high = -DBL_MAX;
            bool integer = true, finite = true, seen = false, any_missing = false;
            for (size_t r = 0; r < b->rows; r++)
            {
                double value = values[r];
                if (value == missing)
                {
                    any_missing = true;
                    continue;
                }
                seen = true;
                if (!isfinite(value) || value > missing)
                    finite = false;
                if (floor(value) != value || (value == 0 && signbit(value)))
                    integer = false;
                if (value < low)
                    low = value;
                if (value > high)
                    high = value;
            }
            double range = high - low;
            if (!seen)
                p->kind = CIMPORT_CACHE_MISSING;
            else if (finite && integer && range == 0 && !any_missing)
            {
                p->kind = CIMPORT_CACHE_CONSTANT;
                p->base = low;
            }
            else if (finite && integer && range <= (double)UINT32_MAX - 1)
            {
                p->base = low;
                if (range <= UINT8_MAX - 1)
                {
                    p->kind = CIMPORT_CACHE_U8;
                    width = 1;
                }
                else if (range <= UINT16_MAX - 1)
                {
                    p->kind = CIMPORT_CACHE_U16;
                    width = 2;
                }
                else
                {
                    p->kind = CIMPORT_CACHE_U32;
                    width = 4;
                }
            }
            else if (finite && !integer && !as_double && low >= -FLT_MAX && high <= FLT_MAX &&
                     (as_float || (low >= -2147483647.0 && high <= 2147483620.0)))
            {
                /* Provisional float storage is consumed only if the final Stata
                 * column is float. A later double/string promotion reparses the
                 * affected block from its original bytes, preserving precision. */
                p->kind = CIMPORT_CACHE_FLOAT;
                width = 4;
            }
            else
            {
                p->kind = CIMPORT_CACHE_DOUBLE;
                width = 8;
            }
        }
        bytes = (bytes + 7) & ~(size_t)7;
        p->offset = bytes;
        bytes += width * b->rows;
    }
    return bytes;
}

CImportCachedBlock *cimport_cache_pack(CImportCacheBuilder *b, CImportCacheBudget *budget,
                                       double missing, bool as_float, bool as_double,
                                       const char *file_data, char quote)
{
    if (!b->rows || b->wide)
        return NULL;
    size_t bytes = cimport_cache_plan(b, missing, as_float, as_double, file_data != NULL);
    uint_fast64_t used = atomic_load(&budget->used);
    do
    {
        if (used > budget->limit || bytes > budget->limit - used)
            return NULL;
    } while (!atomic_compare_exchange_weak(&budget->used, &used, used + bytes));
    CImportCachedBlock *block = malloc(bytes);
    if (!block)
    {
        atomic_fetch_sub(&budget->used, bytes);
        return NULL;
    }
    block->next = NULL;
    block->file_start = b->file_start;
    block->file_end = b->file_end;
    block->rows = b->rows;
    block->bytes = bytes;
    memcpy(block->columns, b->plan, (size_t)b->columns * sizeof(*b->plan));
    for (int col = 0; col < b->columns; col++)
    {
        const CImportCachedColumn *p = &b->plan[col];
        const double *values = b->values + (size_t)col * b->capacity;
        const uint64_t *fields = b->fields + (size_t)col * b->capacity;
        void *data = (char *)block + p->offset;
        switch (p->kind)
        {
        case CIMPORT_CACHE_MISSING:
        case CIMPORT_CACHE_CONSTANT:
            break;
        case CIMPORT_CACHE_U8:
            for (size_t r = 0; r < b->rows; r++)
                ((uint8_t *)data)[r] =
                    values[r] == missing ? UINT8_MAX : (uint8_t)(values[r] - p->base);
            break;
        case CIMPORT_CACHE_U16:
            for (size_t r = 0; r < b->rows; r++)
                ((uint16_t *)data)[r] =
                    values[r] == missing ? UINT16_MAX : (uint16_t)(values[r] - p->base);
            break;
        case CIMPORT_CACHE_U32:
            for (size_t r = 0; r < b->rows; r++)
                ((uint32_t *)data)[r] =
                    values[r] == missing ? UINT32_MAX : (uint32_t)(values[r] - p->base);
            break;
        case CIMPORT_CACHE_FLOAT:
            for (size_t r = 0; r < b->rows; r++)
            {
                float value = (float)values[r];
                uint32_t bits;
                memcpy(&bits, &value, sizeof(bits));
                ((uint32_t *)data)[r] = values[r] == missing ? CIMPORT_CACHE_FLOAT_MISSING : bits;
            }
            break;
        case CIMPORT_CACHE_DOUBLE:
            memcpy(data, values, b->rows * sizeof(double));
            break;
        case CIMPORT_CACHE_TEXT:
            memcpy(data, fields, b->rows * sizeof(uint64_t));
            break;
        case CIMPORT_CACHE_FIXED_TEXT:
            for (size_t r = 0; r < b->rows; r++)
                ((uint32_t *)data)[r] = (uint32_t)(fields[r] >> 32);
            break;
        case CIMPORT_CACHE_INLINE_TEXT:
            if (p->length)
            {
                for (size_t r = 0; r < b->rows; r++)
                {
                    CImportFieldRef field = {b->file_start + (fields[r] >> 32),
                                             (uint32_t)fields[r]};
                    char *text = (char *)data + r * p->length;
                    if (!quote || !(field.length & CIMPORT_FIELD_QUOTED_FLAG))
                    {
                        const char *src = file_data + field.offset;
                        size_t length = CIMPORT_FIELD_LENGTH(field);
                        while (length && (src[length - 1] == '\r' || src[length - 1] == '\n'))
                            length--;
                        memcpy(text, src, length);
                        text[length] = '\0';
                    }
                    else
                        cimport_extract_field_fast(file_data, &field, text, (int)p->length, quote);
                }
            }
            break;
        }
    }
    b->rows = 0;
    return block;
}

void cimport_cache_free(CImportCachedBlock *block, CImportCacheBudget *budget)
{
    while (block)
    {
        CImportCachedBlock *next = block->next;
        atomic_fetch_sub(&budget->used, block->bytes);
        free(block);
        block = next;
    }
}

double cimport_cache_number(const CImportCachedBlock *block, int col, size_t row, double missing)
{
    const CImportCachedColumn *p = &block->columns[col];
    const void *data = (const char *)block + p->offset;
    switch (p->kind)
    {
    case CIMPORT_CACHE_CONSTANT:
        return p->base;
    case CIMPORT_CACHE_U8:
    {
        uint8_t value = ((const uint8_t *)data)[row];
        return value == UINT8_MAX ? missing : p->base + value;
    }
    case CIMPORT_CACHE_U16:
    {
        uint16_t value = ((const uint16_t *)data)[row];
        return value == UINT16_MAX ? missing : p->base + value;
    }
    case CIMPORT_CACHE_U32:
    {
        uint32_t value = ((const uint32_t *)data)[row];
        return value == UINT32_MAX ? missing : p->base + value;
    }
    case CIMPORT_CACHE_FLOAT:
    {
        uint32_t bits = ((const uint32_t *)data)[row];
        float value;
        memcpy(&value, &bits, sizeof(value));
        return bits == CIMPORT_CACHE_FLOAT_MISSING ? missing : (double)value;
    }
    case CIMPORT_CACHE_DOUBLE:
        return ((const double *)data)[row];
    default:
        return missing;
    }
}

CImportFieldRef cimport_cache_field(const CImportCachedBlock *block, int col, size_t row)
{
    const CImportCachedColumn *p = &block->columns[col];
    const void *data = (const char *)block + p->offset;
    if (p->kind == CIMPORT_CACHE_FIXED_TEXT)
        return (CImportFieldRef){block->file_start + ((const uint32_t *)data)[row], p->length};
    uint64_t field = ((const uint64_t *)data)[row];
    return (CImportFieldRef){block->file_start + (field >> 32), (uint32_t)field};
}

const char *cimport_cache_text(const CImportCachedBlock *block, int col, size_t row)
{
    const CImportCachedColumn *p = &block->columns[col];
    return p->length ? (const char *)block + p->offset + row * p->length : "";
}
