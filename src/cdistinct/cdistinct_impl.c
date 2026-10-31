/* Exact tuple counts. Loaded rows and permutations use zero-based indices.
 * All allocations are local to a call. */
#include "cdistinct_impl.h"
#include "ctools_types.h"
#include "ctools_order.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include <stdio.h>
#include <stdlib.h>

static int missing_key(const stata_data *d, size_t row, int nkeys)
{
    for (int k = 0; k < nkeys; k++)
    {
        const stata_variable *v = &d->vars[k];
        if (v->type == STATA_TYPE_DOUBLE ? SF_is_missing(v->data.dbl[row])
                                         : !v->data.str[row] || !v->data.str[row][0])
            return 1;
    }
    return 0;
}

/* The table borrows rows from the loaded columns. Zero denotes an empty slot;
 * occupied slots hold row + 1 and a hash used only to reject unequal keys.
 * Hash matches always receive an exact comparison, including missing codes. */
typedef struct
{
    size_t row;
    uint64_t hash;
} cdistinct_entry;

static uint64_t cdistinct_mix(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static int cdistinct_row_hash(const stata_data *d, size_t row, int nk, int missing, double missval,
                              uint64_t *hash)
{
    uint64_t h = 0;
    for (int k = 0; k < nk; k++)
    {
        const stata_variable *v = &d->vars[k];
        uint64_t bits = 0;
        if (v->type == STATA_TYPE_DOUBLE)
        {
            double x = v->data.dbl[row];
            if (!missing && x >= missval)
                return 0;
            /* Stata treats negative and positive zero as the same value. */
            if (x != 0)
                memcpy(&bits, &x, sizeof(bits));
        }
        else
        {
            const unsigned char *s = (const unsigned char *)v->data.str[row];
            if (!missing && !*s)
                return 0;
            bits = UINT64_C(14695981039346656037);
            for (; *s; s++)
                bits = (bits ^ *s) * UINT64_C(1099511628211);
        }
        h = cdistinct_mix(h ^ bits ^ UINT64_C(0x9e3779b97f4a7c15));
    }
    *hash = h;
    return 1;
}

static ST_retcode cdistinct_grow(cdistinct_entry **table, size_t *capacity)
{
    if (*capacity > SIZE_MAX / 2 / sizeof(**table))
        return 920;
    size_t next = *capacity ? *capacity * 2 : 64;
    cdistinct_entry *entries = calloc(next, sizeof(*entries));
    if (!entries)
        return 920;
    for (size_t j = 0; j < *capacity; j++)
    {
        if (!(*table)[j].row)
            continue;
        size_t slot = (*table)[j].hash & (next - 1);
        while (entries[slot].row)
            slot = (slot + 1) & (next - 1);
        entries[slot] = (*table)[j];
    }
    free(*table);
    *table = entries;
    *capacity = next;
    return 0;
}

static ST_retcode cdistinct_hash_count(const stata_data *d, const int *keys, int nk, int missing,
                                       size_t *total, size_t *distinct, int *sort)
{
    cdistinct_entry *table = NULL;
    size_t capacity = 0;
    ST_retcode rc = 0;
    double missval = SV_missval;
    *total = *distinct = 0;
    *sort = 0;
    for (size_t i = 0; i < d->nobs; i++)
    {
        uint64_t hash;
        if (!cdistinct_row_hash(d, i, nk, missing, missval, &hash))
            continue;
        (*total)++;
        if (!capacity && (rc = cdistinct_grow(&table, &capacity)))
            break;
        size_t slot = hash & (capacity - 1);
        while (table[slot].row)
        {
            if (table[slot].hash == hash &&
                !ctools_compare_rows(d, table[slot].row - 1, keys, d, i, keys, nk))
                break;
            slot = (slot + 1) & (capacity - 1);
        }
        if (table[slot].row)
            continue;
        if (*distinct >= capacity - capacity / 4)
        {
            /* Bound the trial table at 64 KiB. High-cardinality inputs then
             * use the shared radix/merge sorter instead of a growing table. */
            if (capacity == 4096)
            {
                *sort = 1;
                break;
            }
            if ((rc = cdistinct_grow(&table, &capacity)))
                break;
            slot = hash & (capacity - 1);
            while (table[slot].row)
                slot = (slot + 1) & (capacity - 1);
        }
        table[slot] = (cdistinct_entry){i + 1, hash};
        (*distinct)++;
    }
    free(table);
    return rc;
}

static ST_retcode cdistinct_fast_count(stata_data *d, const int *keys, int nk, int missing,
                                       size_t *total, size_t *distinct)
{
    size_t n = d->nobs, i = 0;
    *total = *distinct = 0;
    /* An ordered input needs neither a hash table nor a permutation. */
    for (; i < n; i++)
    {
        int cmp = i ? ctools_compare_rows(d, i - 1, keys, d, i, keys, nk) : -1;
        if (cmp > 0)
            break;
        if (!missing && missing_key(d, i, nk))
            continue;
        if (cmp || !*total)
            (*distinct)++;
        (*total)++;
    }
    if (i == n)
        return 0;
    if (i < n)
    {
        int sort = 0;
        ST_retcode rc = cdistinct_hash_count(d, keys, nk, missing, total, distinct, &sort);
        if (rc || !sort)
            return rc;
        if (n > SIZE_MAX / sizeof(*d->sort_order))
            return 920;
        d->sort_order = malloc(n * sizeof(*d->sort_order));
        if (!d->sort_order)
            return 920;
        for (size_t j = 0; j < n; j++)
            d->sort_order[j] = j;
        rc = ctools_stata_rc(ctools_order_stable(d, keys, nk));
        if (rc)
            return rc;
    }
    *total = *distinct = 0;
    size_t previous = 0;
    for (size_t p = 0; p < n; p++)
    {
        size_t row = d->sort_order ? d->sort_order[p] : p;
        if (!missing && missing_key(d, row, nk))
            continue;
        if (!*total || ctools_compare_rows(d, previous, keys, d, row, keys, nk))
            (*distinct)++;
        (*total)++;
        previous = row;
    }
    return 0;
}

ST_retcode cdistinct_main(const char *args)
{
    int nk = 0;
    if (ctools_parse_int_option(args, "nkeys", &nk) != CTOOLS_PARSE_OK || nk < 1 ||
        nk >= SF_nvars() || SF_nvars() != nk + 1 || SF_var_is_string(nk + 1))
        return 198;
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    int *cols = malloc((size_t)nk * sizeof(*cols));
    ST_retcode rc = 920;
    if (!cols)
        goto done;
    for (int k = 0; k < nk; k++)
        cols[k] = k + 1;
    double start = ctools_timer_seconds();
    rc = ctools_stata_rc(ctools_data_load(&fd, cols, nk, 1, SF_nobs(),
                                          CTOOLS_LOAD_MARKER | CTOOLS_LOAD_NO_SORT_ORDER));
    if (rc)
        goto done;
    double loaded = ctools_timer_seconds();
    for (int k = 0; k < nk; k++)
        cols[k] = k;
    size_t total, distinct;
    rc = cdistinct_fast_count(&fd.data, cols, nk, ctools_parse_bool_option(args, "missing"), &total,
                              &distinct);
    if (rc)
        goto done;
    char value[32];
    snprintf(value, sizeof(value), "%zu", total);
    rc = SF_macro_save("_ctools_count_N", value);
    if (rc)
        goto done;
    snprintf(value, sizeof(value), "%zu", distinct);
    rc = SF_macro_save("_ctools_count_distinct", value);
    if (ctools_parse_bool_option(args, "verbose"))
        ctools_msg("cdistinct", "Load %.3fs; count %.3fs", loaded - start,
                   ctools_timer_seconds() - loaded);
done:
    free(cols);
    ctools_filtered_data_free(&fd);
    return rc;
}
