/* LSD sort of mixed numeric and string keys against a stable reference merge sort.
   Strings are allocated at their exact length, so AddressSanitizer reports any read
   past a terminator. Keys cover empty strings, every length up to the width, shared
   prefixes of 1-15 bytes or a two-valued first byte, bytes above 127, NULL strings,
   -0, extended missing values, heavy duplicates, and identity or shuffled incoming
   orders; row counts cross the single-thread, cache-split and parallel thresholds. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_sort.h"
#include "ctools_convert.h"

static int threads = 1;
int ctools_get_max_threads(void) { return threads; }

static uint64_t rng = 88172645463325252ULL;
static uint64_t next(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static int uni(int n) { return (int)(next() % (uint64_t)n); }

static const stata_variable *ref_vars;
static const int *ref_keys;
static int ref_nkeys;

static int compare(perm_idx_t a, perm_idx_t b)
{
    for (int k = 0; k < ref_nkeys; k++) {
        const stata_variable *v = &ref_vars[ref_keys[k] - 1];
        if (v->type == STATA_TYPE_DOUBLE) {
            uint64_t x = ctools_double_to_sortable(v->data.dbl[a]);
            uint64_t y = ctools_double_to_sortable(v->data.dbl[b]);
            if (x != y) return x < y ? -1 : 1;
        } else {
            /* strcmp compares bytes as unsigned char */
            int r = strcmp(v->data.str[a] ? v->data.str[a] : "", v->data.str[b] ? v->data.str[b] : "");
            if (r) return r < 0 ? -1 : 1;
        }
    }
    return 0;
}

static void merge_sort(perm_idx_t *a, perm_idx_t *t, size_t n)
{
    if (n < 2) return;
    size_t m = n / 2, i = 0, j = m, k = 0;
    merge_sort(a, t, m);
    merge_sort(a + m, t, n - m);
    while (i < m && j < n) t[k++] = compare(a[j], a[i]) < 0 ? a[j++] : a[i++];
    while (i < m) t[k++] = a[i++];
    while (j < n) t[k++] = a[j++];
    memcpy(a, t, n * sizeof(perm_idx_t));
}

static const unsigned char alphabet[] = {' ', 'a', 'b', 'z', '0', '~', 0x80, 0xC8, 0xFF};

/* Length 0..width; the first `prefix` bytes are shared by every value of the column,
   or with prefix < 0 the first byte is 'a' or 'b', which leaves two buckets of
   many rows that still differ in later bytes. */
static char *random_string(int width, int prefix)
{
    int len = uni(width + 1);
    char *s = malloc((size_t)len + 1);
    assert(s);
    for (int i = 0; i < len; i++) {
        s[i] = (char)(i < prefix ? 'P' : alphabet[uni(sizeof alphabet)]);
        if (i == 0 && prefix < 0) s[i] = uni(2) ? 'a' : 'b';
    }
    s[len] = 0;
    return s;
}

static void fill_column(stata_variable *v, size_t n)
{
    static const int widths[] = {1, 3, 7, 8, 9, 15, 16, 24};
    v->nobs = n;
    if (uni(3) == 0) {
        int kind = uni(3);
        v->type = STATA_TYPE_DOUBLE;
        v->data.dbl = malloc(n * sizeof(double));
        assert(v->data.dbl);
        for (size_t i = 0; i < n; i++) {
            double x = kind == 0 ? (double)uni(4) : kind == 1 ? (double)uni(2000001) / 1000 - 1000
                                                              : (uni(9) ? (double)uni(3) - 1 : -0.0);
            if (uni(25) == 0) {
                /* ., .a, ..., .z */
                uint64_t bits = 0x7FE0000000000000ULL + ((uint64_t)uni(27) << 40);
                memcpy(&x, &bits, sizeof x);
            }
            v->data.dbl[i] = x;
        }
        return;
    }
    int width = widths[uni(8)], prefix = uni(3) == 0 ? -1 : uni(2) ? uni(width) : 0;
    int npool = 1 + uni(6);
    char *pool[6];
    for (int p = 0; p < npool; p++) pool[p] = random_string(width, prefix);
    int few = uni(3) == 0;
    v->type = STATA_TYPE_STRING;
    v->str_maxlen = (size_t)width;
    v->data.str = malloc(n * sizeof(char *));
    assert(v->data.str);
    for (size_t i = 0; i < n; i++) {
        if (uni(1000) == 0) {
            v->data.str[i] = NULL;
        } else if (few) {
            const char *s = pool[uni(npool)];
            v->data.str[i] = malloc(strlen(s) + 1);
            assert(v->data.str[i]);
            strcpy(v->data.str[i], s);
        } else {
            v->data.str[i] = random_string(width, prefix);
        }
    }
    for (int p = 0; p < npool; p++) free(pool[p]);
}

static void run(size_t n, int nvars, int team, int shuffle)
{
    stata_variable vars[3];
    stata_data data = {n, (size_t)nvars, vars, malloc(n * sizeof(perm_idx_t))};
    perm_idx_t *expect = malloc(n * sizeof(perm_idx_t)), *scratch = malloc(n * sizeof(perm_idx_t));
    int keys[3];
    assert(data.sort_order && expect && scratch);
    memset(vars, 0, sizeof vars);
    for (int c = 0; c < nvars; c++) {
        fill_column(&vars[c], n);
        keys[c] = nvars - c;  /* key order differs from storage order */
    }
    for (size_t i = 0; i < n; i++) data.sort_order[i] = (perm_idx_t)i;
    for (size_t i = n; shuffle && i > 1; i--) {
        size_t j = next() % i;
        perm_idx_t t = data.sort_order[i - 1];
        data.sort_order[i - 1] = data.sort_order[j];
        data.sort_order[j] = t;
    }
    memcpy(expect, data.sort_order, n * sizeof(perm_idx_t));
    ref_vars = vars;
    ref_keys = keys;
    ref_nkeys = nvars;
    merge_sort(expect, scratch, n);

    threads = team;
    if (ctools_sort_radix_lsd_order_only(&data, keys, (size_t)nvars) != STATA_OK ||
        memcmp(expect, data.sort_order, n * sizeof(perm_idx_t)) != 0) {
        fprintf(stderr, "mismatch: n=%zu keys=%d threads=%d shuffled=%d\n", n, nvars, team, shuffle);
        abort();
    }
    for (int c = 0; c < nvars; c++) {
        if (vars[c].type == STATA_TYPE_DOUBLE) {
            free(vars[c].data.dbl);
        } else {
            for (size_t i = 0; i < n; i++) free(vars[c].data.str[i]);
            free(vars[c].data.str);
        }
    }
    free(data.sort_order);
    free(expect);
    free(scratch);
}

int main(void)
{
    /* Parallel runs of tied records at least 200,000 long, and runs split for the cache. */
    for (int iteration = 0; iteration < 120; iteration++) {
        int band = uni(10);
        size_t n = band < 5 ? 1 + (size_t)uni(200)
                 : band < 8 ? 1000 + (size_t)uni(20000)
                 : band < 9 ? 140000 + (size_t)uni(60000)
                            : 420000 + (size_t)uni(200000);
        run(n, 1 + uni(3), 1 + uni(4), uni(2));
    }
    return 0;
}
