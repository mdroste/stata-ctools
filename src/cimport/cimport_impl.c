/*
 * cimport_impl.c - High-Performance CSV Import for Stata
 *
 * Part of the ctools suite. Replaces "import delimited" with a
 * C-accelerated implementation featuring:
 *   - Multi-threaded parallel parsing
 *   - SIMD-accelerated newline/delimiter scanning
 *   - Custom fast float parser
 *   - Arena allocator for memory efficiency
 *   - Column-major caching for fast SPI loading
 *   - Optional direct DTA file writing (fast mode)
 *
 * Based on fastimport by Claude (Anthropic)
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdatomic.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <pthread.h>

/* Platform-specific includes for memory-mapped files */
#ifdef _WIN32
    #include <windows.h>
    #include <io.h>
    #define CIMPORT_WINDOWS 1
#else
    #include <sys/mman.h>
    #include <sys/stat.h>
    #include <fcntl.h>
    #include <unistd.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "cimport_impl.h"

/* High-resolution timing - use shared timer module */
#define cimport_get_time_ms() ctools_timer_ms()

/* SIMD headers */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define CIMPORT_USE_NEON 1
#elif defined(__SSE2__)
    #include <emmintrin.h>
    #define CIMPORT_USE_SSE2 1
#endif

/* ============================================================================
 * Configuration Constants
 * Use shared constants from ctools_config.h where applicable.
 * ============================================================================ */

/* Shared constants (from ctools_config.h):
   - CTOOLS_MAX_COLUMNS: max CSV columns (32767)
   - CTOOLS_MAX_VARNAME_LEN: max variable name length (32)
   - CTOOLS_MAX_STRING_LEN: max string variable length (2045)
   - CTOOLS_IMPORT_CHUNK_SIZE: bytes per parsing chunk (8MB)
   - CTOOLS_IO_MAX_THREADS: max threads for I/O (16)
   - CTOOLS_ARENA_BLOCK_SIZE: arena allocator block size (1MB)
*/

/* Local: cimport may use more threads for parsing than general I/O */
#define CIMPORT_MAX_THREADS      32

/* ============================================================================
 * Type Definitions
 * ============================================================================ */

typedef enum {
    CIMPORT_COL_UNKNOWN = 0,
    CIMPORT_COL_NUMERIC,
    CIMPORT_COL_STRING
} CImportColumnType;

typedef enum {
    CIMPORT_BINDQUOTES_LOOSE = 0,   /* Default: each line is a row, ignore quotes for row boundaries */
    CIMPORT_BINDQUOTES_STRICT = 1   /* Respect quotes: quoted fields can span multiple lines */
} CImportBindQuotesMode;

typedef enum {
    CIMPORT_NUM_DOUBLE = 0,
    CIMPORT_NUM_FLOAT,
    CIMPORT_NUM_LONG,
    CIMPORT_NUM_INT,
    CIMPORT_NUM_BYTE
} CImportNumericSubtype;

typedef struct {
    char name[CTOOLS_MAX_VARNAME_LEN];
    CImportColumnType type;
    int max_strlen;
    CImportNumericSubtype num_subtype;
    bool is_integer;
    double min_value;
    double max_value;
    int max_decimal_digits;
} CImportColumnInfo;

typedef struct {
    uint64_t offset;
    uint32_t length;
} CImportFieldRef;

typedef struct {
    uint16_t num_fields;
    CImportFieldRef fields[];
} CImportParsedRow;

typedef struct CImportArenaBlock {
    struct CImportArenaBlock *next;
    size_t used;
    size_t capacity;
    char data[];
} CImportArenaBlock;

typedef struct {
    CImportArenaBlock *first;
    CImportArenaBlock *current;
    size_t total_allocated;
} CImportArena;

typedef struct {
    bool seen_string;
    bool seen_non_empty;
    int max_field_len;
} CImportColumnParseStats;

typedef struct {
    CImportParsedRow **rows;
    size_t num_rows;
    size_t capacity;
    CImportArena arena;
    CImportColumnParseStats *col_stats;
    int num_col_stats;
    int max_fields_in_chunk;
} CImportParsedChunk;

typedef struct {
    double *numeric_data;
    char **string_data;
    CImportArena string_arena;
    size_t count;
} CImportColumnCache;

typedef struct {
    char *file_data;
    size_t file_size;
    char *filename;
    CImportColumnInfo *columns;
    int num_columns;
    size_t total_rows;
    char delimiter;
    char quote_char;
    bool has_header;
    CImportBindQuotesMode bindquotes;
    CImportParsedChunk *chunks;
    int num_chunks;
    size_t *row_offsets;
    CImportColumnCache *col_cache;
    char *string_arena;
    size_t string_arena_size;
    size_t string_arena_used;
    int num_threads;
    atomic_int error_code;
    char error_message[256];
    bool is_loaded;
    bool cache_ready;
    bool verbose;
    int max_fields_seen;
    double time_mmap;
    double time_parse;
    double time_type_infer;
    double time_cache;
#ifdef CIMPORT_WINDOWS
    HANDLE file_handle;
    HANDLE mapping_handle;
#endif
} CImportContext;

/* Global cached context */
static CImportContext *g_cimport_ctx = NULL;

/* ============================================================================
 * Arena Allocator
 * ============================================================================ */

static void cimport_arena_init(CImportArena *arena) {
    arena->first = NULL;
    arena->current = NULL;
    arena->total_allocated = 0;
}

static void *cimport_arena_alloc(CImportArena *arena, size_t size) {
    size = (size + 7) & ~7;

    if (!arena->current || arena->current->used + size > arena->current->capacity) {
        size_t block_size = CTOOLS_ARENA_BLOCK_SIZE;
        if (size > block_size - sizeof(CImportArenaBlock)) {
            block_size = size + sizeof(CImportArenaBlock);
        }

        CImportArenaBlock *block = malloc(block_size);
        if (!block) return NULL;

        block->next = NULL;
        block->used = 0;
        block->capacity = block_size - sizeof(CImportArenaBlock);

        if (arena->current) {
            arena->current->next = block;
        } else {
            arena->first = block;
        }
        arena->current = block;
        arena->total_allocated += block_size;
    }

    void *ptr = arena->current->data + arena->current->used;
    arena->current->used += size;
    return ptr;
}

static void cimport_arena_free(CImportArena *arena) {
    CImportArenaBlock *block = arena->first;
    while (block) {
        CImportArenaBlock *next = block->next;
        free(block);
        block = next;
    }
    arena->first = NULL;
    arena->current = NULL;
    arena->total_allocated = 0;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

static void cimport_display_msg(const char *msg) {
    SF_display((char *)msg);
}

static void cimport_display_error(const char *msg) {
    SF_error((char *)msg);
}

static inline bool cimport_is_whitespace(char c) {
    return c == ' ' || c == '\t' || c == '\r';
}

/* ============================================================================
 * Fast Float Parser
 * NOTE: Fast double parsing is now provided by ctools_types.h:
 *   - ctools_parse_double_fast(str, len, result, missval)
 *   - ctools_pow10_table[]
 * ============================================================================ */

/* ============================================================================
 * SIMD-Accelerated Scanning
 * ============================================================================ */

#if CIMPORT_USE_NEON
static inline const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim) {
    const uint8x16_t v_newline = vdupq_n_u8('\n');
    const uint8x16_t v_delim = vdupq_n_u8(delim);
    const uint8x16_t v_quote = vdupq_n_u8('"');

    while (ptr + 16 <= end) {
        uint8x16_t chunk = vld1q_u8((const uint8_t *)ptr);
        uint8x16_t cmp_nl = vceqq_u8(chunk, v_newline);
        uint8x16_t cmp_delim = vceqq_u8(chunk, v_delim);
        uint8x16_t cmp_quote = vceqq_u8(chunk, v_quote);
        uint8x16_t cmp = vorrq_u8(vorrq_u8(cmp_nl, cmp_delim), cmp_quote);

        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        if (vgetq_lane_u64(cmp64, 0) || vgetq_lane_u64(cmp64, 1)) {
            for (int i = 0; i < 16 && ptr + i < end; i++) {
                char c = ptr[i];
                if (c == '\n' || c == delim || c == '"') return ptr + i;
            }
        }
        ptr += 16;
    }

    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}

#elif CIMPORT_USE_SSE2
static inline const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim) {
    const __m128i v_newline = _mm_set1_epi8('\n');
    const __m128i v_delim = _mm_set1_epi8(delim);
    const __m128i v_quote = _mm_set1_epi8('"');

    while (ptr + 16 <= end) {
        __m128i chunk = _mm_loadu_si128((const __m128i *)ptr);
        __m128i cmp_nl = _mm_cmpeq_epi8(chunk, v_newline);
        __m128i cmp_delim = _mm_cmpeq_epi8(chunk, v_delim);
        __m128i cmp_quote = _mm_cmpeq_epi8(chunk, v_quote);
        __m128i cmp = _mm_or_si128(_mm_or_si128(cmp_nl, cmp_delim), cmp_quote);
        int mask = _mm_movemask_epi8(cmp);

        if (mask) {
            int pos = __builtin_ctz(mask);
            return ptr + pos;
        }
        ptr += 16;
    }

    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}

#else
static inline const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim) {
    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}
#endif

/* ============================================================================
 * CSV Parsing Engine
 * ============================================================================ */

/* Find next row boundary - LOOSE mode: each line is a row (ignore quotes) */
static inline const char *cimport_find_next_row_loose(const char *ptr, const char *end) {
    while (ptr < end) {
        if (*ptr == '\n') {
            return ptr + 1;
        }
        ptr++;
    }
    return end;
}

/* Find next row boundary - STRICT mode: respect quotes (fields can span lines) */
static inline const char *cimport_find_next_row_strict(const char *ptr, const char *end, char quote) {
    bool in_quotes = false;

    while (ptr < end) {
        char c = *ptr;

        if (c == quote) {
            in_quotes = !in_quotes;
        } else if (!in_quotes && c == '\n') {
            return ptr + 1;
        }
        ptr++;
    }
    return end;
}

/* Wrapper that chooses based on bindquotes mode */
static inline const char *cimport_find_next_row(const char *ptr, const char *end, char quote, CImportBindQuotesMode bindquotes) {
    if (bindquotes == CIMPORT_BINDQUOTES_STRICT) {
        return cimport_find_next_row_strict(ptr, end, quote);
    } else {
        return cimport_find_next_row_loose(ptr, end);
    }
}

static int cimport_parse_row_fast(const char *start, const char *end, char delim, char quote,
                                   CImportFieldRef *fields, int max_fields, const char *file_base) {
    int field_count = 0;
    const char *ptr = start;
    const char *field_start = start;
    bool in_quotes = false;

    while (ptr < end && field_count < max_fields) {
        if (!in_quotes) {
            const char *found = cimport_find_delim_or_newline_simd(ptr, end, delim);

            if (found < end) {
                char c = *found;

                if (c == '"') {
                    in_quotes = true;
                    ptr = found + 1;
                    continue;
                }

                fields[field_count].offset = (uint64_t)(field_start - file_base);
                fields[field_count].length = (uint32_t)(found - field_start);
                field_count++;

                if (c == '\r' || c == '\n') {
                    return field_count;
                }

                field_start = found + 1;
                ptr = field_start;
                continue;
            }
            ptr = end;
        } else {
            char c = *ptr;
            if (c == quote) {
                if (ptr + 1 < end && *(ptr + 1) == quote) {
                    ptr += 2;
                    continue;
                }
                in_quotes = false;
            }
            ptr++;
        }
    }

    if (field_start < end && field_count < max_fields) {
        const char *field_end = end;
        while (field_end > field_start && (field_end[-1] == '\r' || field_end[-1] == '\n')) {
            field_end--;
        }
        if (field_end > field_start) {
            fields[field_count].offset = (uint64_t)(field_start - file_base);
            fields[field_count].length = (uint32_t)(field_end - field_start);
            field_count++;
        }
    }

    return field_count;
}

static inline bool cimport_field_contains_quote(const char *src, int len, char quote) {
    while (len >= 8) {
        if (src[0] == quote || src[1] == quote || src[2] == quote || src[3] == quote ||
            src[4] == quote || src[5] == quote || src[6] == quote || src[7] == quote) {
            return true;
        }
        src += 8;
        len -= 8;
    }
    while (len > 0) {
        if (*src == quote) return true;
        src++;
        len--;
    }
    return false;
}

static int cimport_extract_field_fast(const char *file_base, CImportFieldRef *field, char *output, int max_len, char quote) {
    const char *src = file_base + field->offset;
    int src_len = field->length;
    int out_len = 0;

    while (src_len > 0 && (*src == ' ' || *src == '\t')) {
        src++;
        src_len--;
    }

    while (src_len > 0 && (src[src_len-1] == ' ' || src[src_len-1] == '\t' ||
                           src[src_len-1] == '\r' || src[src_len-1] == '\n')) {
        src_len--;
    }

    bool is_quoted = (src_len >= 2 && src[0] == quote && src[src_len-1] == quote);

    if (is_quoted) {
        src++;
        src_len -= 2;

        for (int i = 0; i < src_len && out_len < max_len - 1; i++) {
            if (src[i] == quote && i + 1 < src_len && src[i + 1] == quote) {
                output[out_len++] = quote;
                i++;
            } else {
                output[out_len++] = src[i];
            }
        }
    } else {
        int copy_len = (src_len < max_len - 1) ? src_len : max_len - 1;
        memcpy(output, src, copy_len);
        out_len = copy_len;
    }

    output[out_len] = '\0';
    return out_len;
}

static inline int cimport_extract_field_unquoted(const char *file_base, CImportFieldRef *field, char *output, int max_len) {
    const char *src = file_base + field->offset;
    int src_len = field->length;

    while (src_len > 0 && (*src == ' ' || *src == '\t')) {
        src++;
        src_len--;
    }

    while (src_len > 0 && (src[src_len-1] == ' ' || src[src_len-1] == '\t' ||
                           src[src_len-1] == '\r' || src[src_len-1] == '\n')) {
        src_len--;
    }

    int copy_len = (src_len < max_len - 1) ? src_len : max_len - 1;
    memcpy(output, src, copy_len);
    output[copy_len] = '\0';
    return copy_len;
}

static inline bool cimport_field_looks_numeric(const char *src, int len) {
    while (len > 0 && (*src == ' ' || *src == '\t')) { src++; len--; }
    if (len == 0) return true;
    if (*src == '"' && len >= 2) { src++; len -= 2; }
    if (len == 0) return true;

    char c = *src;
    if (c >= '0' && c <= '9') return true;
    if (c == '-' || c == '+' || c == '.') return true;

    if (len >= 2 && (c == 'N' || c == 'n')) {
        char c2 = src[1];
        if (c2 == 'A' || c2 == 'a') return true;
        if (len >= 3 && (c2 == 'a' || c2 == 'A')) {
            char c3 = src[2];
            if (c3 == 'N' || c3 == 'n') return true;
        }
    }

    return false;
}

static inline bool cimport_analyze_numeric_fast(const char *file_base, CImportFieldRef *field, char quote,
                                                 double *out_value, bool *out_is_integer) {
    const char *src = file_base + field->offset;
    int len = field->length;

    while (len > 0 && (*src == ' ' || *src == '\t' || *src == quote)) { src++; len--; }
    while (len > 0 && (src[len-1] == ' ' || src[len-1] == '\t' || src[len-1] == quote ||
                       src[len-1] == '\r' || src[len-1] == '\n')) { len--; }

    if (len == 0) return false;

    if (len == 1 && *src == '.') return false;
    if (len == 2 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'A' || src[1] == 'a')) return false;
    if (len == 3 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'a' || src[1] == 'A') &&
        (src[2] == 'N' || src[2] == 'n')) return false;

    double val;
    if (!ctools_parse_double_fast(src, len, &val, SV_missval)) return false;

    *out_value = val;

    *out_is_integer = true;
    for (int i = 0; i < len; i++) {
        if (src[i] == '.') {
            *out_is_integer = false;
            break;
        }
    }

    return true;
}

static CImportNumericSubtype cimport_determine_numeric_subtype(CImportColumnInfo *col) {
    if (col->min_value < -2147483647.0 || col->max_value > 2147483620.0) {
        return CIMPORT_NUM_DOUBLE;
    }

    if (!col->is_integer) {
        return CIMPORT_NUM_FLOAT;
    }

    if (col->min_value >= -127 && col->max_value <= 100) {
        return CIMPORT_NUM_BYTE;
    }
    if (col->min_value >= -32767 && col->max_value <= 32740) {
        return CIMPORT_NUM_INT;
    }

    return CIMPORT_NUM_LONG;
}

/* ============================================================================
 * Memory-Mapped File I/O (Cross-Platform)
 * ============================================================================ */

#ifdef CIMPORT_WINDOWS

static int cimport_mmap_file(CImportContext *ctx, const char *filename) {
    /* Open file */
    ctx->file_handle = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ,
                                    NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (ctx->file_handle == INVALID_HANDLE_VALUE) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot open file: error %lu", GetLastError());
        return -1;
    }

    /* Get file size */
    LARGE_INTEGER file_size;
    if (!GetFileSizeEx(ctx->file_handle, &file_size)) {
        CloseHandle(ctx->file_handle);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot get file size: error %lu", GetLastError());
        return -1;
    }

    ctx->file_size = (size_t)file_size.QuadPart;
    if (ctx->file_size == 0) {
        CloseHandle(ctx->file_handle);
        strcpy(ctx->error_message, "File is empty");
        return -1;
    }

    /* Create file mapping */
    ctx->mapping_handle = CreateFileMappingA(ctx->file_handle, NULL, PAGE_READONLY,
                                              0, 0, NULL);
    if (ctx->mapping_handle == NULL) {
        CloseHandle(ctx->file_handle);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot create file mapping: error %lu", GetLastError());
        return -1;
    }

    /* Map view of file */
    ctx->file_data = (char *)MapViewOfFile(ctx->mapping_handle, FILE_MAP_READ,
                                            0, 0, ctx->file_size);
    if (ctx->file_data == NULL) {
        CloseHandle(ctx->mapping_handle);
        CloseHandle(ctx->file_handle);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot map file: error %lu", GetLastError());
        return -1;
    }

    return 0;
}

static void cimport_munmap_file(CImportContext *ctx) {
    if (ctx->file_data) {
        UnmapViewOfFile(ctx->file_data);
        ctx->file_data = NULL;
    }
    if (ctx->mapping_handle) {
        CloseHandle(ctx->mapping_handle);
        ctx->mapping_handle = NULL;
    }
    if (ctx->file_handle) {
        CloseHandle(ctx->file_handle);
        ctx->file_handle = NULL;
    }
}

#else /* Unix/POSIX */

static int cimport_mmap_file(CImportContext *ctx, const char *filename) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot open file: %s", strerror(errno));
        return -1;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot stat file: %s", strerror(errno));
        return -1;
    }

    ctx->file_size = st.st_size;
    if (ctx->file_size == 0) {
        close(fd);
        strcpy(ctx->error_message, "File is empty");
        return -1;
    }

    ctx->file_data = mmap(NULL, ctx->file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);

    if (ctx->file_data == MAP_FAILED) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot mmap file: %s", strerror(errno));
        ctx->file_data = NULL;
        return -1;
    }

#ifdef __linux__
    madvise(ctx->file_data, ctx->file_size, MADV_SEQUENTIAL | MADV_WILLNEED);
#elif defined(__APPLE__)
    madvise(ctx->file_data, ctx->file_size, MADV_SEQUENTIAL);
#endif

    return 0;
}

static void cimport_munmap_file(CImportContext *ctx) {
    if (ctx->file_data) {
        munmap(ctx->file_data, ctx->file_size);
        ctx->file_data = NULL;
    }
}

#endif /* CIMPORT_WINDOWS */

/* ============================================================================
 * Variable Name Sanitization
 * ============================================================================ */

static void cimport_sanitize_varname(const char *input, int len, char *output) {
    int j = 0;

    while (len > 0 && (cimport_is_whitespace(*input) || *input == '"')) {
        input++;
        len--;
    }
    while (len > 0 && (cimport_is_whitespace(input[len-1]) || input[len-1] == '"')) {
        len--;
    }

    for (int i = 0; i < len && j < CTOOLS_MAX_VARNAME_LEN - 1; i++) {
        char c = input[i];
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_') {
            output[j++] = c;
        } else if (c >= '0' && c <= '9' && j > 0) {
            output[j++] = c;
        }
    }

    if (j == 0) {
        strcpy(output, "v");
        j = 1;
    }

    output[j] = '\0';
}

/* ============================================================================
 * Header and Type Inference
 * ============================================================================ */

static int cimport_parse_header(CImportContext *ctx) {
    CImportFieldRef fields[CTOOLS_MAX_COLUMNS];
    const char *end = cimport_find_next_row(ctx->file_data, ctx->file_data + ctx->file_size, ctx->quote_char, ctx->bindquotes);

    ctx->num_columns = cimport_parse_row_fast(ctx->file_data, end, ctx->delimiter, ctx->quote_char,
                                               fields, CTOOLS_MAX_COLUMNS, ctx->file_data);

    if (ctx->num_columns <= 0) {
        strcpy(ctx->error_message, "No columns found in header");
        return -1;
    }

    ctx->columns = calloc(ctx->num_columns, sizeof(CImportColumnInfo));
    if (!ctx->columns) return -1;

    char name_buf[CTOOLS_MAX_STRING_LEN];
    for (int i = 0; i < ctx->num_columns; i++) {
        if (ctx->has_header) {
            cimport_extract_field_fast(ctx->file_data, &fields[i], name_buf, sizeof(name_buf), ctx->quote_char);
            cimport_sanitize_varname(name_buf, strlen(name_buf), ctx->columns[i].name);
        } else {
            snprintf(ctx->columns[i].name, CTOOLS_MAX_VARNAME_LEN, "v%d", i + 1);
        }
        ctx->columns[i].type = CIMPORT_COL_UNKNOWN;
        ctx->columns[i].max_strlen = 0;
    }

    return 0;
}

static void cimport_infer_column_types(CImportContext *ctx) {
    for (int i = 0; i < ctx->num_columns; i++) {
        ctx->columns[i].is_integer = true;
        ctx->columns[i].min_value = DBL_MAX;
        ctx->columns[i].max_value = -DBL_MAX;
        ctx->columns[i].max_decimal_digits = 0;
        ctx->columns[i].num_subtype = CIMPORT_NUM_DOUBLE;
    }

    for (int col = 0; col < ctx->num_columns; col++) {
        CImportColumnInfo *colinfo = &ctx->columns[col];
        bool is_string = false;
        int max_len = 0;

        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            if (chunk->col_stats && col < chunk->num_col_stats) {
                CImportColumnParseStats *stats = &chunk->col_stats[col];
                if (stats->seen_string) {
                    is_string = true;
                }
                if (stats->max_field_len > max_len) {
                    max_len = stats->max_field_len;
                }
            }
        }

        if (is_string) {
            colinfo->type = CIMPORT_COL_STRING;
            colinfo->max_strlen = max_len > 0 ? max_len : 1;
            if (colinfo->max_strlen > CTOOLS_MAX_STRING_LEN) {
                colinfo->max_strlen = CTOOLS_MAX_STRING_LEN;
            }
        } else {
            colinfo->type = CIMPORT_COL_NUMERIC;
        }
    }

    size_t total_data_rows = 0;
    for (int c = 0; c < ctx->num_chunks; c++) {
        size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;
        total_data_rows += ctx->chunks[c].num_rows - start_row;
    }

    for (int col_idx = 0; col_idx < ctx->num_columns; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        if (col->type != CIMPORT_COL_NUMERIC) continue;

        size_t sample_count = 0;
        size_t max_samples = 1000;
        size_t sample_step = (total_data_rows > max_samples) ? (total_data_rows / max_samples) : 1;
        size_t row_idx = 0;

        for (int c = 0; c < ctx->num_chunks && sample_count < max_samples; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

            for (size_t r = start_row; r < chunk->num_rows && sample_count < max_samples; r++) {
                if (row_idx % sample_step == 0) {
                    CImportParsedRow *row = chunk->rows[r];
                    if (col_idx < row->num_fields) {
                        CImportFieldRef *field = &row->fields[col_idx];
                        double val;
                        bool is_int;

                        if (cimport_analyze_numeric_fast(ctx->file_data, field, ctx->quote_char, &val, &is_int)) {
                            if (!is_int) {
                                col->is_integer = false;
                            }
                            if (val < col->min_value) {
                                col->min_value = val;
                            }
                            if (val > col->max_value) {
                                col->max_value = val;
                            }
                            sample_count++;
                        }
                    }
                }
                row_idx++;
            }
        }
    }

    for (int i = 0; i < ctx->num_columns; i++) {
        if (ctx->columns[i].type == CIMPORT_COL_UNKNOWN) {
            ctx->columns[i].type = CIMPORT_COL_NUMERIC;
        }
        if (ctx->columns[i].type == CIMPORT_COL_STRING && ctx->columns[i].max_strlen == 0) {
            ctx->columns[i].max_strlen = 1;
        }
        if (ctx->columns[i].max_strlen > CTOOLS_MAX_STRING_LEN) {
            ctx->columns[i].max_strlen = CTOOLS_MAX_STRING_LEN;
        }

        if (ctx->columns[i].type == CIMPORT_COL_NUMERIC) {
            if (ctx->columns[i].min_value == DBL_MAX) {
                ctx->columns[i].num_subtype = CIMPORT_NUM_DOUBLE;
            } else {
                ctx->columns[i].num_subtype = cimport_determine_numeric_subtype(&ctx->columns[i]);
            }
        }
    }
}

/* ============================================================================
 * Context Management
 * ============================================================================ */

static void cimport_free_context(CImportContext *ctx) {
    if (!ctx) return;

    if (ctx->chunks) {
        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            free(chunk->rows);
            free(chunk->col_stats);
            cimport_arena_free(&chunk->arena);
        }
        free(ctx->chunks);
    }

    if (ctx->col_cache) {
        for (int i = 0; i < ctx->num_columns; i++) {
            free(ctx->col_cache[i].numeric_data);
            free(ctx->col_cache[i].string_data);
            cimport_arena_free(&ctx->col_cache[i].string_arena);
        }
        free(ctx->col_cache);
    }
    free(ctx->string_arena);

    free(ctx->columns);
    free(ctx->filename);
    free(ctx->row_offsets);
    cimport_munmap_file(ctx);
    free(ctx);
}

static void cimport_clear_cached_context(void) {
    if (g_cimport_ctx) {
        cimport_free_context(g_cimport_ctx);
        g_cimport_ctx = NULL;
    }
}

/* Forward declarations */
static CImportContext *cimport_parse_csv(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes);
static void cimport_build_column_cache(CImportContext *ctx);
static ST_retcode cimport_do_scan(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes);
static ST_retcode cimport_do_load(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes);

/* ============================================================================
 * Parallel Chunk Parsing
 * ============================================================================ */

typedef struct {
    CImportContext *ctx;
    int chunk_id;
    const char *start;
    const char *end;
    CImportParsedChunk *chunk;
} CImportChunkParseTask;

static size_t cimport_find_safe_boundary(const char *data, size_t file_size, size_t target, char quote, CImportBindQuotesMode bindquotes) {
    if (target >= file_size) return file_size;

    const char *ptr = data + target;
    const char *start = data;

    /* Find the nearest newline before target */
    while (ptr > start && ptr[-1] != '\n') {
        ptr--;
    }

    /* In loose mode, just use the nearest newline - no quote checking needed */
    if (bindquotes == CIMPORT_BINDQUOTES_LOOSE) {
        return ptr - data;
    }

    /* In strict mode, check if we're inside a quoted field */
    size_t scan_limit = 10000;
    const char *scan_start = ptr;
    if ((size_t)(ptr - start) > scan_limit) {
        scan_start = ptr - scan_limit;
        while (scan_start > start && scan_start[-1] != '\n') scan_start--;
    } else {
        scan_start = start;
    }

    int quote_count = 0;
    for (const char *p = scan_start; p < ptr; p++) {
        if (*p == quote) quote_count++;
    }

    if (quote_count % 2 == 1) {
        ptr = cimport_find_next_row_strict(ptr, data + file_size, quote);
    }

    return ptr - data;
}

static void *cimport_parse_chunk_parallel(void *arg) {
    CImportChunkParseTask *task = (CImportChunkParseTask *)arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = task->chunk;

    cimport_arena_init(&chunk->arena);
    chunk->capacity = 16384;
    chunk->rows = malloc(sizeof(CImportParsedRow *) * chunk->capacity);
    chunk->num_rows = 0;

    chunk->num_col_stats = ctx->num_columns;
    chunk->col_stats = calloc(ctx->num_columns, sizeof(CImportColumnParseStats));

    if (!chunk->rows || !chunk->col_stats) {
        atomic_store(&ctx->error_code, 1);
        return NULL;
    }

    CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
    const char *ptr = task->start;
    const char *end = task->end;

    while (ptr < end && atomic_load(&ctx->error_code) == 0) {
        const char *row_end = cimport_find_next_row(ptr, end, ctx->quote_char, ctx->bindquotes);
        if (row_end > end) row_end = end;

        int num_fields = cimport_parse_row_fast(ptr, row_end, ctx->delimiter, ctx->quote_char,
                                                 field_buf, CTOOLS_MAX_COLUMNS, ctx->file_data);

        if (num_fields == 0 || (num_fields == 1 && field_buf[0].length == 0)) {
            ptr = row_end;
            continue;
        }

        bool is_header_row = (ctx->has_header && task->chunk_id == 0 && chunk->num_rows == 0);
        if (!is_header_row) {
            for (int f = 0; f < num_fields && f < chunk->num_col_stats; f++) {
                CImportColumnParseStats *stats = &chunk->col_stats[f];
                CImportFieldRef *field = &field_buf[f];

                if ((int)field->length > stats->max_field_len) {
                    stats->max_field_len = field->length;
                }

                if (!stats->seen_string && field->length > 0) {
                    stats->seen_non_empty = true;
                    const char *src = ctx->file_data + field->offset;
                    if (!cimport_field_looks_numeric(src, field->length)) {
                        stats->seen_string = true;
                    }
                }
            }
        }

        if (chunk->num_rows >= chunk->capacity) {
            chunk->capacity *= 2;
            CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * chunk->capacity);
            if (!new_rows) {
                atomic_store(&ctx->error_code, 1);
                break;
            }
            chunk->rows = new_rows;
        }

        size_t row_size = sizeof(CImportParsedRow) + sizeof(CImportFieldRef) * num_fields;
        CImportParsedRow *row = cimport_arena_alloc(&chunk->arena, row_size);
        if (!row) {
            atomic_store(&ctx->error_code, 1);
            break;
        }

        row->num_fields = num_fields;
        memcpy(row->fields, field_buf, sizeof(CImportFieldRef) * num_fields);

        chunk->rows[chunk->num_rows++] = row;

        if (num_fields > chunk->max_fields_in_chunk) {
            chunk->max_fields_in_chunk = num_fields;
        }

        ptr = row_end;
    }

    return NULL;
}

/* ============================================================================
 * Main Parse Function
 * ============================================================================ */

static CImportContext *cimport_parse_csv(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes) {
    CImportContext *ctx = calloc(1, sizeof(CImportContext));
    if (!ctx) return NULL;

    double t_start, t_end;

    ctx->delimiter = delimiter;
    ctx->quote_char = '"';
    ctx->has_header = has_header;
    ctx->bindquotes = bindquotes;
    ctx->filename = strdup(filename);
    ctx->verbose = verbose;
    atomic_init(&ctx->error_code, 0);

    t_start = cimport_get_time_ms();
    if (cimport_mmap_file(ctx, filename) != 0) {
        cimport_display_error(ctx->error_message);
        cimport_free_context(ctx);
        return NULL;
    }
    t_end = cimport_get_time_ms();
    ctx->time_mmap = t_end - t_start;

    if (cimport_parse_header(ctx) != 0) {
        cimport_display_error("Failed to parse header");
        cimport_free_context(ctx);
        return NULL;
    }

#ifdef CIMPORT_WINDOWS
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    int cpu_count = (int)sysinfo.dwNumberOfProcessors;
#else
    int cpu_count = (int)sysconf(_SC_NPROCESSORS_ONLN);
#endif
    if (cpu_count <= 0) cpu_count = 4;
    if (cpu_count > CIMPORT_MAX_THREADS) cpu_count = CIMPORT_MAX_THREADS;
    ctx->num_threads = cpu_count;

    t_start = cimport_get_time_ms();

    size_t min_chunk_size = 4 * 1024 * 1024;
    int num_chunks = ctx->num_threads;

    if (ctx->file_size < min_chunk_size * 2) {
        num_chunks = 1;
    }

    ctx->num_chunks = num_chunks;
    ctx->chunks = calloc(ctx->num_chunks, sizeof(CImportParsedChunk));
    if (!ctx->chunks) {
        cimport_free_context(ctx);
        return NULL;
    }

    if (num_chunks == 1) {
        CImportParsedChunk *chunk = &ctx->chunks[0];
        cimport_arena_init(&chunk->arena);
        chunk->capacity = 65536;
        chunk->rows = malloc(sizeof(CImportParsedRow *) * chunk->capacity);
        chunk->num_rows = 0;

        chunk->num_col_stats = ctx->num_columns;
        chunk->col_stats = calloc(ctx->num_columns, sizeof(CImportColumnParseStats));

        if (!chunk->rows || !chunk->col_stats) {
            cimport_free_context(ctx);
            return NULL;
        }

        CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
        const char *ptr = ctx->file_data;
        const char *end = ctx->file_data + ctx->file_size;

        while (ptr < end) {
            const char *row_end = cimport_find_next_row(ptr, end, ctx->quote_char, ctx->bindquotes);

            int num_fields = cimport_parse_row_fast(ptr, row_end, ctx->delimiter, ctx->quote_char,
                                                     field_buf, CTOOLS_MAX_COLUMNS, ctx->file_data);

            if (num_fields == 0 || (num_fields == 1 && field_buf[0].length == 0)) {
                ptr = row_end;
                continue;
            }

            bool is_header_row = (ctx->has_header && chunk->num_rows == 0);
            if (!is_header_row) {
                for (int f = 0; f < num_fields && f < chunk->num_col_stats; f++) {
                    CImportColumnParseStats *stats = &chunk->col_stats[f];
                    CImportFieldRef *field = &field_buf[f];

                    if ((int)field->length > stats->max_field_len) {
                        stats->max_field_len = field->length;
                    }

                    if (!stats->seen_string && field->length > 0) {
                        stats->seen_non_empty = true;
                        const char *src = ctx->file_data + field->offset;
                        if (!cimport_field_looks_numeric(src, field->length)) {
                            stats->seen_string = true;
                        }
                    }
                }
            }

            if (chunk->num_rows >= chunk->capacity) {
                chunk->capacity *= 2;
                CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * chunk->capacity);
                if (!new_rows) {
                    cimport_free_context(ctx);
                    return NULL;
                }
                chunk->rows = new_rows;
            }

            size_t row_size = sizeof(CImportParsedRow) + sizeof(CImportFieldRef) * num_fields;
            CImportParsedRow *row = cimport_arena_alloc(&chunk->arena, row_size);
            if (!row) {
                cimport_free_context(ctx);
                return NULL;
            }

            row->num_fields = num_fields;
            memcpy(row->fields, field_buf, sizeof(CImportFieldRef) * num_fields);
            chunk->rows[chunk->num_rows++] = row;

            if (num_fields > chunk->max_fields_in_chunk) {
                chunk->max_fields_in_chunk = num_fields;
            }

            ptr = row_end;
        }
    } else {
        size_t chunk_size = ctx->file_size / num_chunks;

        size_t *boundaries = malloc((num_chunks + 1) * sizeof(size_t));
        if (!boundaries) {
            cimport_free_context(ctx);
            return NULL;
        }

        boundaries[0] = 0;
        for (int i = 1; i < num_chunks; i++) {
            size_t target = i * chunk_size;
            boundaries[i] = cimport_find_safe_boundary(ctx->file_data, ctx->file_size, target, ctx->quote_char, ctx->bindquotes);
        }
        boundaries[num_chunks] = ctx->file_size;

        pthread_t threads[CIMPORT_MAX_THREADS];
        CImportChunkParseTask tasks[CIMPORT_MAX_THREADS];

        for (int i = 0; i < num_chunks; i++) {
            tasks[i].ctx = ctx;
            tasks[i].chunk_id = i;
            tasks[i].start = ctx->file_data + boundaries[i];
            tasks[i].end = ctx->file_data + boundaries[i + 1];
            tasks[i].chunk = &ctx->chunks[i];

            pthread_create(&threads[i], NULL, cimport_parse_chunk_parallel, &tasks[i]);
        }

        for (int i = 0; i < num_chunks; i++) {
            pthread_join(threads[i], NULL);
        }

        free(boundaries);

        if (atomic_load(&ctx->error_code) != 0) {
            cimport_free_context(ctx);
            return NULL;
        }
    }

    t_end = cimport_get_time_ms();
    ctx->time_parse = t_end - t_start;

    ctx->total_rows = 0;
    for (int c = 0; c < ctx->num_chunks; c++) {
        ctx->total_rows += ctx->chunks[c].num_rows;
    }
    if (ctx->has_header && ctx->total_rows > 0) {
        ctx->total_rows--;
    }

    ctx->max_fields_seen = ctx->num_columns;
    for (int c = 0; c < ctx->num_chunks; c++) {
        if (ctx->chunks[c].max_fields_in_chunk > ctx->max_fields_seen) {
            ctx->max_fields_seen = ctx->chunks[c].max_fields_in_chunk;
        }
    }

    if (ctx->max_fields_seen > ctx->num_columns) {
        int old_num = ctx->num_columns;
        int new_num = ctx->max_fields_seen;

        CImportColumnInfo *new_columns = realloc(ctx->columns, new_num * sizeof(CImportColumnInfo));
        if (!new_columns) {
            cimport_free_context(ctx);
            return NULL;
        }
        ctx->columns = new_columns;

        for (int i = old_num; i < new_num; i++) {
            memset(&ctx->columns[i], 0, sizeof(CImportColumnInfo));
            snprintf(ctx->columns[i].name, CTOOLS_MAX_VARNAME_LEN, "v%d", i + 1);
            ctx->columns[i].type = CIMPORT_COL_UNKNOWN;
            ctx->columns[i].max_strlen = 0;
        }

        ctx->num_columns = new_num;
    }

    t_start = cimport_get_time_ms();
    cimport_infer_column_types(ctx);
    t_end = cimport_get_time_ms();
    ctx->time_type_infer = t_end - t_start;

    ctx->is_loaded = true;
    ctx->cache_ready = false;

    return ctx;
}

/* ============================================================================
 * Column Cache Builder
 * ============================================================================ */

typedef struct {
    CImportContext *ctx;
    int col_start;
    int col_end;
    int thread_id;
} CImportCacheBuildTask;

static void *cimport_build_cache_worker(void *arg) {
    CImportCacheBuildTask *task = (CImportCacheBuildTask *)arg;
    CImportContext *ctx = task->ctx;
    char field_buf[CTOOLS_MAX_STRING_LEN + 1];

    for (int col_idx = task->col_start; col_idx < task->col_end; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        CImportColumnCache *cache = &ctx->col_cache[col_idx];
        size_t row_idx = 0;

        if (col->type == CIMPORT_COL_STRING) {
            cimport_arena_init(&cache->string_arena);

            for (int c = 0; c < ctx->num_chunks; c++) {
                CImportParsedChunk *chunk = &ctx->chunks[c];
                size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

                for (size_t r = start_row; r < chunk->num_rows; r++) {
                    CImportParsedRow *row = chunk->rows[r];
                    char *str;

                    if (col_idx < row->num_fields) {
                        CImportFieldRef *field = &row->fields[col_idx];
                        const char *src = ctx->file_data + field->offset;

                        int len;
                        if (!cimport_field_contains_quote(src, field->length, ctx->quote_char)) {
                            len = cimport_extract_field_unquoted(ctx->file_data, field, field_buf, CTOOLS_MAX_STRING_LEN);
                        } else {
                            len = cimport_extract_field_fast(ctx->file_data, field, field_buf, CTOOLS_MAX_STRING_LEN, ctx->quote_char);
                        }
                        field_buf[len] = '\0';

                        str = cimport_arena_alloc(&cache->string_arena, len + 1);
                        if (str) {
                            memcpy(str, field_buf, len + 1);
                        }
                    } else {
                        str = cimport_arena_alloc(&cache->string_arena, 1);
                        if (str) str[0] = '\0';
                    }

                    cache->string_data[row_idx++] = str;
                }
            }
        } else {
            double missing = SV_missval;

            for (int c = 0; c < ctx->num_chunks; c++) {
                CImportParsedChunk *chunk = &ctx->chunks[c];
                size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

                for (size_t r = start_row; r < chunk->num_rows; r++) {
                    CImportParsedRow *row = chunk->rows[r];
                    double val;

                    if (col_idx < row->num_fields) {
                        CImportFieldRef *field = &row->fields[col_idx];
                        const char *src = ctx->file_data + field->offset;

                        if (!ctools_parse_double_fast(src, field->length, &val, missing)) {
                            val = missing;
                        }
                    } else {
                        val = missing;
                    }

                    cache->numeric_data[row_idx++] = val;
                }
            }
        }

        cache->count = row_idx;
    }

    return NULL;
}

static void cimport_build_column_cache(CImportContext *ctx) {
    if (ctx->cache_ready) return;

    double t_start = cimport_get_time_ms();

    ctx->col_cache = calloc(ctx->num_columns, sizeof(CImportColumnCache));
    if (!ctx->col_cache) return;

    for (int i = 0; i < ctx->num_columns; i++) {
        CImportColumnInfo *col = &ctx->columns[i];
        CImportColumnCache *cache = &ctx->col_cache[i];

        if (col->type == CIMPORT_COL_STRING) {
            cache->string_data = calloc(ctx->total_rows, sizeof(char *));
            cache->numeric_data = NULL;
        } else {
            cache->numeric_data = malloc(ctx->total_rows * sizeof(double));
            cache->string_data = NULL;
        }
        cache->count = 0;
    }

    int num_threads = ctx->num_threads;
    if (num_threads > ctx->num_columns) num_threads = ctx->num_columns;
    if (num_threads < 1) num_threads = 1;

    if (num_threads == 1 || ctx->num_columns < 4) {
        CImportCacheBuildTask task = {ctx, 0, ctx->num_columns, 0};
        cimport_build_cache_worker(&task);
    } else {
        pthread_t threads[CIMPORT_MAX_THREADS];
        CImportCacheBuildTask tasks[CIMPORT_MAX_THREADS];

        int cols_per_thread = (ctx->num_columns + num_threads - 1) / num_threads;

        for (int t = 0; t < num_threads; t++) {
            tasks[t].ctx = ctx;
            tasks[t].col_start = t * cols_per_thread;
            tasks[t].col_end = (t + 1) * cols_per_thread;
            if (tasks[t].col_end > ctx->num_columns) tasks[t].col_end = ctx->num_columns;
            tasks[t].thread_id = t;

            pthread_create(&threads[t], NULL, cimport_build_cache_worker, &tasks[t]);
        }

        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }
    }

    ctx->cache_ready = true;

    double t_end = cimport_get_time_ms();
    ctx->time_cache = t_end - t_start;
}

/* ============================================================================
 * SCAN Mode
 * ============================================================================ */

static ST_retcode cimport_do_scan(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes) {
    char msg[512];
    char macro_val[65536];

    cimport_clear_cached_context();

    g_cimport_ctx = cimport_parse_csv(filename, delimiter, has_header, verbose, bindquotes);
    if (!g_cimport_ctx) {
        return 601;
    }

    CImportContext *ctx = g_cimport_ctx;

    snprintf(macro_val, sizeof(macro_val), "%zu", ctx->total_rows);
    SF_macro_save("_cimport_nobs", macro_val);

    snprintf(macro_val, sizeof(macro_val), "%d", ctx->num_columns);
    SF_macro_save("_cimport_nvar", macro_val);

    char *p = macro_val;
    for (int i = 0; i < ctx->num_columns; i++) {
        if (i > 0) *p++ = ' ';
        strcpy(p, ctx->columns[i].name);
        p += strlen(ctx->columns[i].name);
    }
    *p = '\0';
    SF_macro_save("_cimport_varnames", macro_val);

    p = macro_val;
    for (int i = 0; i < ctx->num_columns; i++) {
        if (i > 0) *p++ = ' ';
        *p++ = (ctx->columns[i].type == CIMPORT_COL_STRING) ? '1' : '0';
    }
    *p = '\0';
    SF_macro_save("_cimport_vartypes", macro_val);

    p = macro_val;
    for (int i = 0; i < ctx->num_columns; i++) {
        if (i > 0) *p++ = ' ';
        if (ctx->columns[i].type == CIMPORT_COL_STRING) {
            *p++ = '0';
        } else {
            *p++ = '0' + (int)ctx->columns[i].num_subtype;
        }
    }
    *p = '\0';
    SF_macro_save("_cimport_numtypes", macro_val);

    p = macro_val;
    for (int i = 0; i < ctx->num_columns; i++) {
        if (i > 0) *p++ = ' ';
        int len = ctx->columns[i].type == CIMPORT_COL_STRING ? ctx->columns[i].max_strlen : 0;
        p += sprintf(p, "%d", len);
    }
    *p = '\0';
    SF_macro_save("_cimport_strlens", macro_val);

    snprintf(msg, sizeof(msg), "Scanned %zu rows, %d columns\n", ctx->total_rows, ctx->num_columns);
    cimport_display_msg(msg);

    if (verbose) {
        snprintf(msg, sizeof(msg), "  [C timing] mmap: %.2f ms, parse: %.2f ms, type_infer: %.2f ms\n",
                 ctx->time_mmap, ctx->time_parse, ctx->time_type_infer);
        cimport_display_msg(msg);
    }

    return 0;
}

/* ============================================================================
 * LOAD Mode
 * ============================================================================ */

static ST_retcode cimport_do_load(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes) {
    char msg[512];

    CImportContext *ctx = g_cimport_ctx;
    bool used_cache = (ctx != NULL && ctx->is_loaded && ctx->filename && strcmp(ctx->filename, filename) == 0);

    if (!used_cache) {
        cimport_clear_cached_context();
        g_cimport_ctx = cimport_parse_csv(filename, delimiter, has_header, verbose, bindquotes);
        ctx = g_cimport_ctx;
        if (!ctx) {
            return 601;
        }
    }

    if (!ctx->cache_ready) {
        cimport_build_column_cache(ctx);
        if (!ctx->cache_ready) {
            cimport_display_error("Failed to build column cache\n");
            cimport_clear_cached_context();
            return 459;
        }
    }

    double t_load_start = cimport_get_time_ms();
    long total_fields = 0;

    ST_int stata_nvars = SF_nvars();
    ST_int stata_nobs = SF_nobs();

    if (verbose) {
        snprintf(msg, sizeof(msg), "  [debug] SF_nvars()=%d, SF_nobs()=%d, ctx->num_columns=%d, ctx->total_rows=%zu\n",
                 stata_nvars, stata_nobs, ctx->num_columns, ctx->total_rows);
        cimport_display_msg(msg);
    }

    if (stata_nvars < ctx->num_columns) {
        snprintf(msg, sizeof(msg), "Error: Stata has %d variables but CSV has %d columns\n",
                 stata_nvars, ctx->num_columns);
        cimport_display_error(msg);
        cimport_clear_cached_context();
        return 198;
    }

    if ((size_t)stata_nobs < ctx->total_rows) {
        snprintf(msg, sizeof(msg), "Error: Stata has %d observations but CSV has %zu rows\n",
                 stata_nobs, ctx->total_rows);
        cimport_display_error(msg);
        cimport_clear_cached_context();
        return 198;
    }

    for (int col_idx = 0; col_idx < ctx->num_columns; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        CImportColumnCache *cache = &ctx->col_cache[col_idx];
        ST_int var = col_idx + 1;
        size_t nrows = cache->count;

        if (col->type == CIMPORT_COL_STRING) {
            char **strings = cache->string_data;

            for (size_t p = 0; p < 16 && p < nrows; p++) {
                __builtin_prefetch(strings[p], 0, 0);
            }

            size_t i = 0;
            size_t nrows_aligned = nrows & ~(size_t)3;

            for (; i < nrows_aligned; i += 4) {
                if (i + 16 < nrows) {
                    __builtin_prefetch(&strings[i + 16], 0, 0);
                    __builtin_prefetch(strings[i + 12], 0, 0);
                    __builtin_prefetch(strings[i + 13], 0, 0);
                    __builtin_prefetch(strings[i + 14], 0, 0);
                    __builtin_prefetch(strings[i + 15], 0, 0);
                }

                SF_sstore(var, (ST_int)(i + 1), strings[i]);
                SF_sstore(var, (ST_int)(i + 2), strings[i + 1]);
                SF_sstore(var, (ST_int)(i + 3), strings[i + 2]);
                SF_sstore(var, (ST_int)(i + 4), strings[i + 3]);
            }

            for (; i < nrows; i++) {
                SF_sstore(var, (ST_int)(i + 1), strings[i]);
            }
            total_fields += nrows;
        } else {
            double *values = cache->numeric_data;

            __builtin_prefetch(values, 0, 0);
            __builtin_prefetch(values + 8, 0, 0);
            __builtin_prefetch(values + 16, 0, 0);
            __builtin_prefetch(values + 24, 0, 0);

            size_t i = 0;
            size_t nrows_aligned = nrows & ~(size_t)7;

            for (; i < nrows_aligned; i += 8) {
                __builtin_prefetch(values + i + 32, 0, 0);
                __builtin_prefetch(values + i + 40, 0, 0);

                SF_vstore(var, (ST_int)(i + 1), values[i]);
                SF_vstore(var, (ST_int)(i + 2), values[i + 1]);
                SF_vstore(var, (ST_int)(i + 3), values[i + 2]);
                SF_vstore(var, (ST_int)(i + 4), values[i + 3]);
                SF_vstore(var, (ST_int)(i + 5), values[i + 4]);
                SF_vstore(var, (ST_int)(i + 6), values[i + 5]);
                SF_vstore(var, (ST_int)(i + 7), values[i + 6]);
                SF_vstore(var, (ST_int)(i + 8), values[i + 7]);
            }

            for (; i < nrows; i++) {
                SF_vstore(var, (ST_int)(i + 1), values[i]);
            }
            total_fields += nrows;
        }
    }

    double t_load_end = cimport_get_time_ms();
    double time_spi = t_load_end - t_load_start;

    snprintf(msg, sizeof(msg), "Loaded %zu observations\n", ctx->total_rows);
    cimport_display_msg(msg);

    if (verbose) {
        snprintf(msg, sizeof(msg), "  [C timing] SPI store: %.2f ms (%ld fields, %.0f fields/sec)\n",
                 time_spi, total_fields, total_fields / (time_spi / 1000.0));
        cimport_display_msg(msg);
    }

    cimport_clear_cached_context();

    return 0;
}

/* ============================================================================
 * Plugin Entry Point
 * ============================================================================ */

ST_retcode cimport_main(const char *args) {
    if (args == NULL || strlen(args) == 0) {
        cimport_display_error("cimport: no arguments specified\n");
        return 198;
    }

    /* Parse arguments: mode filename [delimiter] [noheader] [verbose] [bindquotes=loose|strict] */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *mode = NULL;
    char *filename = NULL;
    char delimiter = ',';
    bool has_header = true;
    bool verbose = false;
    CImportBindQuotesMode bindquotes = CIMPORT_BINDQUOTES_LOOSE;  /* Default: loose (matches Stata) */

    char *token = strtok(args_copy, " ");
    int arg_idx = 0;

    while (token != NULL) {
        if (arg_idx == 0) {
            mode = token;
        } else if (arg_idx == 1) {
            filename = token;
        } else {
            if (strcmp(token, "noheader") == 0) {
                has_header = false;
            } else if (strcmp(token, "verbose") == 0) {
                verbose = true;
            } else if (strcmp(token, "tab") == 0) {
                delimiter = '\t';
            } else if (strcmp(token, "bindquotes=strict") == 0) {
                bindquotes = CIMPORT_BINDQUOTES_STRICT;
            } else if (strcmp(token, "bindquotes=loose") == 0) {
                bindquotes = CIMPORT_BINDQUOTES_LOOSE;
            } else if (strlen(token) == 1) {
                delimiter = token[0];
            } else if (strlen(token) == 3 && token[0] == '"' && token[2] == '"') {
                delimiter = token[1];
            }
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    if (mode == NULL || filename == NULL) {
        cimport_display_error("cimport: mode and filename required\n");
        return 198;
    }

    if (strcmp(mode, "scan") == 0) {
        return cimport_do_scan(filename, delimiter, has_header, verbose, bindquotes);
    } else if (strcmp(mode, "load") == 0) {
        return cimport_do_load(filename, delimiter, has_header, verbose, bindquotes);
    } else {
        cimport_display_error("cimport: invalid mode. Use 'scan' or 'load'\n");
        return 198;
    }
}
