/*
    cmerge_impl.c
    Optimized C-Accelerated Merge for Stata

    Plugin calls, in order:
    - load_using: load keys and selected payload columns, sort, and cache them.
    - plan: in the untouched master, load the master keys alone, sort and join
      them. Join errors therefore leave the master unchanged.
    - execute: write the output. An identity output reuses the master already
      in Stata; otherwise whole master observations are staged row by row and
      written once in output order, and the staged rows are kept as an undo
      copy. Datasets whose strings are too wide to stage load, assemble and
      keep full columns instead.
    - rollback or commit: restore the master from the undo copy, or release it.
    - On the column path, output string arrays borrow contents from the loaded
      inputs, which remain alive until output storage and cleanup finish.

    Merge types: 1:1, m:1, 1:m, m:m, and 1:1 by observation number.

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <pthread.h>
#include <stdatomic.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/* C11 atomics for platforms without GCC/Clang/MSVC intrinsics */
#if !defined(_WIN32) && !defined(__GNUC__) && !defined(__clang__)
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_ATOMICS__)
#define CMERGE_USE_C11_ATOMICS 1
#include <stdatomic.h>
#endif
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_threads.h"
#include "ctools_spi.h"
#include "ctools_config.h"
#include "ctools_rowstage.h"
#include "ctools_arena.h"
#include "ctools_release.h"
#include "ctools_parse.h"
#include "ctools_file.h"
#include "cmerge_impl.h"

/* Helper modules */
#include "ctools_sort_pairs.h"
#include "cmerge_keys.h"
#include "cmerge_group_search.h"
#include "cmerge_join.h"
#include "cmerge_io.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CMERGE_MAX_KEYVARS 32
/* Longest Stata str# value */
#define CMERGE_STR_MAXLEN 2045

/* Per-variable index lists arrive through locals (the token names a local
 * the ado set), because a list in the plugin argument exceeds Stata's
 * plugin-call argument length on wide data. The buffers grow to the widest
 * merge seen and are released by cmerge_cleanup_cache. */
static int *g_cm_ku_indices = NULL;     /* load_using: keepusing positions */
static size_t g_cm_ku_cap = 0;
static int *g_cm_placeholders = NULL;   /* execute: keepusing destinations */
static size_t g_cm_placeholders_cap = 0;
static int *g_cm_shared = NULL;         /* execute: 1 if a destination is a master variable */
static size_t g_cm_shared_cap = 0;

static int *cmerge_grow(int **buf, size_t *cap, size_t n)
{
    if (n == 0) n = 1;
    if (n > *cap || !*buf) {
        int *grown = (int *)realloc(*buf, n * sizeof(int));
        if (!grown) return NULL;
        *buf = grown;
        *cap = n;
    }
    return *buf;
}

/* Read n integers from local macro `name` (the ado's local, prefixed by "_").
 * Fails on truncation, missing values or trailing tokens. */
static int cmerge_macro_ints(const char *local, int *out, size_t n)
{
    char name[40];
    if (snprintf(name, sizeof(name), "_%s", local) >= (int)sizeof(name)) return 0;
    if (n > (size_t)INT_MAX / 16) return 0;
    size_t size = n * 12 + 16;
    char *buf = (char *)malloc(size);
    if (!buf) return 0;
    buf[0] = '\0';
    int ok = SF_macro_use(name, buf, (int)size) == 0 && strlen(buf) < size - 1;
    const char *p = buf;
    ok = ok && ctools_parse_int_array(out, n, &p) == 0;
    while (ok && (*p == ' ' || *p == '\t')) p++;
    ok = ok && *p == '\0';
    free(buf);
    return ok;
}

/* Timer: use ctools_timer_ms() directly from ctools_timer.h */

/* ============================================================================
 * Global using data cache (used by cmerge_io.c)
 * ============================================================================ */

cmerge_using_cache_t g_using_cache = {0};

/* Release the keepusing columns and, with strl, their whole strL values. */
static void cmerge_free_keepusing(void)
{
    if (g_using_cache.ku_spans || g_using_cache.ku_arena) {
        for (int v = 0; v < g_using_cache.n_keepusing; v++) {
            if (g_using_cache.ku_spans) free(g_using_cache.ku_spans[v]);
            if (g_using_cache.ku_arena) free(g_using_cache.ku_arena[v]);
        }
    }
    free(g_using_cache.ku_spans);
    free(g_using_cache.ku_arena);
    g_using_cache.ku_spans = NULL;
    g_using_cache.ku_arena = NULL;
    g_using_cache.has_strl = 0;
    free(g_using_cache.perm);
    g_using_cache.perm = NULL;
    ctools_filtered_data_free(&g_using_cache.keepusing);
}

/* Put the key columns of keys in the order given (row i becomes old row
 * order[i]). A string key is copied into one block of fixed slots in the new
 * order: after a sort its pointers would otherwise point all over the old
 * order, and the join, which reads the keys in sequence, would wait on memory
 * at every comparison. On failure the columns are unchanged. */
static stata_retcode cmerge_order_keys(stata_data *keys, int nkeys, const perm_idx_t *order)
{
    const size_t n = keys->nobs;
    for (int k = 0; k < nkeys; k++) {
        stata_variable *var = &keys->vars[k];
        stata_variable moved;
        memset(&moved, 0, sizeof(moved));
        moved.nobs = n;
        moved.type = var->type;
        if (var->type == STATA_TYPE_DOUBLE) {
            double *values = ctools_safe_cacheline_alloc2(n ? n : 1, sizeof(double));
            if (!values) return STATA_ERR_MEMORY;
            const double *old = var->data.dbl;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if (n >= MIN_OBS_PER_THREAD * 2)
            #endif
            for (size_t i = 0; i < n; i++) values[i] = old[order[i]];
            moved.data.dbl = values;
        } else {
            char *const *old = var->data.str;
            size_t width = var->str_maxlen;
            if (width == 0 || width >= CMERGE_STR_MAXLEN) {
                /* no declared width: the longest value decides the slot */
                size_t longest = 0;
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static) reduction(max:longest) if (n >= MIN_OBS_PER_THREAD * 2)
                #endif
                for (size_t i = 0; i < n; i++) {
                    size_t len = old[i] ? strnlen(old[i], CMERGE_STR_MAXLEN) : 0;
                    if (len > longest) longest = len;
                }
                width = longest;
            }
            const size_t slot = width + 1;
            size_t bytes;
            if (ctools_safe_mul_size(n ? n : 1, slot, &bytes) != 0) return STATA_ERR_MEMORY;
            char **ptrs = ctools_safe_cacheline_alloc2(n ? n : 1, sizeof(char *));
            char *flat = (char *)malloc(bytes);
            ctools_string_arena *arena = (ctools_string_arena *)calloc(1, sizeof(*arena));
            if (!ptrs || !flat || !arena) {
                ctools_aligned_free(ptrs);
                free(flat);
                free(arena);
                return STATA_ERR_MEMORY;
            }
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if (n >= MIN_OBS_PER_THREAD * 2)
            #endif
            for (size_t i = 0; i < n; i++) {
                const char *s = old[order[i]];
                char *dst = flat + i * slot;
                size_t len = s ? strnlen(s, width) : 0;
                memcpy(dst, s ? s : "", len);
                dst[len] = '\0';
                ptrs[i] = dst;
            }
            arena->base = flat;
            arena->capacity = arena->used = bytes;
            arena->mode = CTOOLS_STRING_ARENA_NO_FALLBACK;
            arena->has_fallback = 0;
            moved.data.str = ptrs;
            moved.str_maxlen = width;
            moved._arena = arena;
        }
        stata_variable_free(var);
        *var = moved;
    }
    if (keys->sort_order) {
        for (size_t i = 0; i < n; i++) keys->sort_order[i] = (perm_idx_t)i;
    }
    return STATA_OK;
}

/* The longest value of each string key, as its str_maxlen: the loader may
 * report the 2045-byte maximum, and the sort dispatcher chooses its engine by
 * this length. The keys are only read. */
static void cmerge_tighten_string_keys(stata_data *keys, int nkeys)
{
    const size_t n = keys->nobs;
    for (int k = 0; k < nkeys; k++) {
        stata_variable *var = &keys->vars[k];
        if (var->type != STATA_TYPE_STRING || var->str_maxlen < CMERGE_STR_MAXLEN) continue;
        char *const *s = var->data.str;
        size_t longest = 0;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(max:longest) if (n >= MIN_OBS_PER_THREAD * 2)
        #endif
        for (size_t i = 0; i < n; i++) {
            size_t len = s[i] ? strnlen(s[i], CMERGE_STR_MAXLEN) : 0;
            if (len > longest) longest = len;
        }
        var->str_maxlen = longest;
    }
}

/* A copy of keys' sort order, or NULL when it is the identity (or on
 * allocation failure, with *failed set). */
static perm_idx_t *cmerge_take_order(const stata_data *keys, int *failed)
{
    const size_t n = keys->nobs;
    *failed = 0;
    int identity = 1;
    for (size_t i = 0; i < n && identity; i++) identity = keys->sort_order[i] == (perm_idx_t)i;
    if (identity) return NULL;
    perm_idx_t *order = (perm_idx_t *)ctools_safe_malloc2(n, sizeof(perm_idx_t));
    if (!order) {
        *failed = 1;
        return NULL;
    }
    memcpy(order, keys->sort_order, n * sizeof(perm_idx_t));
    return order;
}

/* Observations obs1 .. obs1 + nobs - 1 of strL variable var_idx (one-based
 * plugin index) as spans into one block: any length, binary values marked.
 * The SPI reads strLs one value at a time, and only on the calling thread:
 * reading rows of one strL variable from several threads corrupted Stata's
 * heap. Each value is read into a probe-sized slot first; only a value that
 * may be longer is measured and read again, so a short value costs one read
 * and the binary flag. */
#define CMERGE_STRL_PROBE 256
static stata_retcode cmerge_read_strl(int var_idx, size_t obs1, size_t nobs,
                                      ctools_dta_span **spans_out, char **arena_out)
{
    *spans_out = NULL;
    *arena_out = NULL;
    ctools_dta_span *spans = ctools_safe_calloc2(nobs ? nobs : 1, sizeof(*spans));
    size_t *at = ctools_safe_malloc2(nobs ? nobs : 1, sizeof(size_t));
    size_t cap = nobs * 32 + CMERGE_STRL_PROBE + 1, used = 0;
    char *arena = spans && at ? malloc(cap) : NULL;
    if (!spans || !at || !arena) {
        free(spans);
        free(at);
        free(arena);
        return STATA_ERR_MEMORY;
    }
    stata_retcode rc = STATA_OK;
    for (size_t i = 0; i < nobs && rc == STATA_OK; i++) {
        ST_int obs = (ST_int)(obs1 + i);
        size_t want = CMERGE_STRL_PROBE;
        for (int pass = 0; pass < 2; pass++) {
            if (cap - used < want + 1) {
                size_t grown = cap;
                while (grown - used < want + 1) grown = grown > SIZE_MAX / 2 ? used + want + 1 : grown * 2;
                char *p = realloc(arena, grown);
                if (!p) {
                    rc = STATA_ERR_MEMORY;
                    break;
                }
                arena = p;
                cap = grown;
            }
            int got = SF_strldata(var_idx, obs, arena + used, (int)want + 1);
            if (got < 0) {
                rc = STATA_ERR_STATA_READ;
                break;
            }
            if ((size_t)got < want) { /* the whole value */
                spans[i].len = (size_t)got;
                break;
            }
            int len = SF_sdatalen(var_idx, obs);
            if (len < 0) {
                rc = STATA_ERR_STATA_READ;
                break;
            }
            if ((size_t)len <= want) { /* it fit the slot */
                if ((size_t)got != (size_t)len) rc = STATA_ERR_STATA_READ;
                spans[i].len = (size_t)len;
                break;
            }
            if (pass) {
                rc = STATA_ERR_STATA_READ;
                break;
            }
            want = (size_t)len; /* read it again whole */
        }
        if (rc != STATA_OK) break;
        at[i] = used;
        arena[used + spans[i].len] = '\0';
        used += spans[i].len + 1;
        spans[i].binary = spans[i].len && SF_var_is_binary(var_idx, obs) ? 1 : 0;
    }
    if (rc == STATA_OK)
        for (size_t i = 0; i < nobs; i++) spans[i].ptr = arena + at[i];
    free(at);
    if (rc != STATA_OK) {
        free(spans);
        free(arena);
        return rc;
    }
    *spans_out = spans;
    *arena_out = arena;
    return STATA_OK;
}

/* Load the keepusing columns (one-based plugin indices). With strl (option
 * tempdta) a strL column is read whole into g_using_cache.ku_spans and its
 * keepusing entry is an empty string column; the others load as usual. */
static stata_retcode cmerge_load_keepusing(int *indices, int n, int strl)
{
    int any_strl = 0;
    for (int v = 0; strl && v < n; v++) any_strl |= SF_var_is_strl(indices[v]) != 0;
    if (!any_strl)
        return ctools_data_load(&g_using_cache.keepusing, indices, (size_t)n, 0, 0,
                                CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_ROW_TILES);

    const size_t nobs = (size_t)SF_nobs();
    int *plain = ctools_safe_malloc2((size_t)n, sizeof(int));
    int *where = ctools_safe_malloc2((size_t)n, sizeof(int));
    ctools_filtered_data tmp;
    ctools_filtered_data_init(&tmp);
    ctools_filtered_data_init(&g_using_cache.keepusing);
    g_using_cache.n_keepusing = n;
    g_using_cache.ku_spans = calloc((size_t)n, sizeof(*g_using_cache.ku_spans));
    g_using_cache.ku_arena = calloc((size_t)n, sizeof(*g_using_cache.ku_arena));
    stata_data *data = &g_using_cache.keepusing.data;
    data->vars = ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, (size_t)n, sizeof(stata_variable));
    stata_retcode rc = (plain && where && g_using_cache.ku_spans && g_using_cache.ku_arena && data->vars)
                           ? STATA_OK : STATA_ERR_MEMORY;
    if (rc == STATA_OK) {
        memset(data->vars, 0, (size_t)n * sizeof(stata_variable));
        data->nvars = (size_t)n;
        data->nobs = nobs;
        g_using_cache.has_strl = 1;
    }
    int nplain = 0;
    for (int v = 0; v < n && rc == STATA_OK; v++) {
        if (!SF_var_is_strl(indices[v])) {
            where[nplain] = v;
            plain[nplain++] = indices[v];
            continue;
        }
        data->vars[v].type = STATA_TYPE_STRING;
        data->vars[v].nobs = nobs;
        rc = cmerge_read_strl(indices[v], 1, nobs, &g_using_cache.ku_spans[v], &g_using_cache.ku_arena[v]);
    }
    if (rc == STATA_OK && nplain) {
        rc = ctools_data_load(&tmp, plain, (size_t)nplain, 0, 0, CTOOLS_LOAD_SKIP_IF);
        if (rc == STATA_OK && (!tmp.data.vars || tmp.data.nobs != nobs)) rc = STATA_ERR_STATA_READ;
        for (int p = 0; p < nplain && rc == STATA_OK; p++) {
            data->vars[where[p]] = tmp.data.vars[p]; /* the cache owns the column now */
            memset(&tmp.data.vars[p], 0, sizeof(stata_variable));
        }
    }
    ctools_filtered_data_free(&tmp);
    free(plain);
    free(where);
    if (rc != STATA_OK) cmerge_free_keepusing();
    return rc;
}

/* Spinlock for protecting cache access during concurrent operations.
 * Prevents race conditions if Stata interrupts/restarts during execution. */
static volatile int g_cache_lock = 0;
static volatile int g_cache_in_use = 0;  /* Set while execute is running */

/* Platform-specific spinlock acquire/release */
static void cache_lock_acquire(void)
{
#if defined(_WIN32)
    while (InterlockedCompareExchange((volatile long *)&g_cache_lock, 1, 0) != 0) {
        /* Spin - brief contention only */
    }
#elif defined(__GNUC__) || defined(__clang__)
    while (__sync_lock_test_and_set(&g_cache_lock, 1)) {
        /* Spin */
    }
#elif defined(CMERGE_USE_C11_ATOMICS)
    /* C11 atomics fallback */
    int expected = 0;
    while (!atomic_compare_exchange_weak((atomic_int *)&g_cache_lock, &expected, 1)) {
        expected = 0;
        /* Spin */
    }
#else
    /* Last resort fallback: volatile spinlock. Not perfectly atomic but provides
     * basic protection. The race window is small (between read and write). */
    while (g_cache_lock) {
        /* Spin */
    }
    g_cache_lock = 1;
#endif
}

static void cache_lock_release(void)
{
#if defined(_WIN32)
    InterlockedExchange((volatile long *)&g_cache_lock, 0);
#elif defined(__GNUC__) || defined(__clang__)
    __sync_lock_release(&g_cache_lock);
#elif defined(CMERGE_USE_C11_ATOMICS)
    /* C11 atomics fallback */
    atomic_store((atomic_int *)&g_cache_lock, 0);
#else
    /* Last resort fallback */
    g_cache_lock = 0;
#endif
}

/* Helper to cleanup cache and clear in-use flag on error paths.
 * This frees the cached using data since the merge operation failed
 * and the cache is no longer useful. Prevents memory leaks in OOM conditions. */
/* The join assumes key order, so a "sorted" hint (the option, or a Stata sort
 * marker) is checked in one pass: 0 if some row sorts before its predecessor.
 * key_cols are 0-based positions in data->vars. Numeric order matches the
 * join (missing codes are ordered doubles above all numbers); strings compare
 * bytewise as in cmerge_compare_keys. */
static int cmerge_keys_in_order(const stata_data *data, const int *key_cols, int nkeys)
{
    for (size_t i = 1; i < data->nobs; i++) {
        int c = 0;
        for (int k = 0; k < nkeys && c == 0; k++) {
            const stata_variable *v = &data->vars[key_cols[k]];
            if (v->type == STATA_TYPE_DOUBLE) {
                double a = v->data.dbl[i - 1], b = v->data.dbl[i];
                c = (a < b) ? -1 : (a > b);
            } else {
                const char *a = v->data.str[i - 1] ? v->data.str[i - 1] : "";
                const char *b = v->data.str[i] ? v->data.str[i] : "";
                c = strcmp(a, b);
            }
        }
        if (c > 0) return 0;
    }
    return 1;
}

static void cache_cleanup_on_error(void)
{
    cache_lock_acquire();
    if (g_using_cache.loaded) {
        ctools_filtered_data_free(&g_using_cache.keys);
        cmerge_free_keepusing();
        g_using_cache.loaded = 0;
    }
    g_cache_in_use = 0;
    cache_lock_release();
}

/* ============================================================================
 * Phase 1: Load using dataset (keys + keepusing only)
 * ============================================================================ */

static ST_retcode cmerge_load_using(const char *args)
{
    ST_retcode rc;

    /* Parse: "load_using <nkeys> <key_idx1>...<key_idxN> n_keepusing <count>
              [keepusing_local <name>] [verbose]"; local <name> holds the
              keepusing positions <idx1>...<idxM>. */
    char args_copy[4096];
    if (strlen(args) >= sizeof(args_copy)) return 198;
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    int nkeys = 0;
    int key_indices[CMERGE_MAX_KEYVARS];
    int n_keepusing = 0;
    int *keepusing_indices = NULL;
    char keepusing_local[32] = {0};
    int verbose = 0;
    int sorted = 0;
    int sorted_hint = 0;   /* Stata sort marker: trusted only if the check passes */
    int merge_by_n = 0;
    int strl = 0;          /* option tempdta: read strL keepusing columns whole */

    char *token = strtok(args_copy, " ");
    int arg_idx = 0;
    int in_keepusing_indices = 0;
    int keepusing_idx = 0;

    while (token != NULL) {
        if (arg_idx == 0) {
            if (!ctools_safe_atoi(token, &nkeys)) {
                SF_error("cmerge: invalid nkeys value\n");
                return 198;
            }
            if (nkeys > CMERGE_MAX_KEYVARS) {
                SF_error("cmerge: too many key variables (max 32)\n");
                return 198;
            }
        }
        else if (arg_idx <= nkeys) {
            if (!ctools_safe_atoi(token, &key_indices[arg_idx - 1])) {
                SF_error("cmerge: invalid key index value\n");
                return 198;
            }
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            if (token) {
                if (!ctools_safe_atoi(token, &n_keepusing)) {
                    SF_error("cmerge: invalid n_keepusing value\n");
                    return 198;
                }
                if (n_keepusing < 0 || n_keepusing > SF_nvars()) {
                    SF_error("cmerge: invalid n_keepusing value\n");
                    return 198;
                }
                keepusing_indices = cmerge_grow(&g_cm_ku_indices, &g_cm_ku_cap,
                                                (size_t)n_keepusing);
                if (!keepusing_indices) {
                    SF_error("cmerge: out of memory\n");
                    return 920;
                }
            }
        }
        else if (strcmp(token, "keepusing_local") == 0) {
            token = strtok(NULL, " ");
            if (!token || strlen(token) >= sizeof(keepusing_local)) {
                SF_error("cmerge: invalid keepusing list\n");
                return 198;
            }
            strcpy(keepusing_local, token);
        }
        else if (strcmp(token, "keepusing_indices") == 0) {
            in_keepusing_indices = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_indices && keepusing_indices && keepusing_idx < n_keepusing) {
            if (!ctools_safe_atoi(token, &keepusing_indices[keepusing_idx])) {
                SF_error("cmerge: invalid keepusing index value\n");
                return 198;
            }
            keepusing_idx++;
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "sorted_hint") == 0) {
            sorted_hint = 1;
        }
        else if (strcmp(token, "merge_by_n") == 0) {
            merge_by_n = 1;
        }
        else if (strcmp(token, "strl") == 0) {
            strl = 1;
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    (void)verbose;  /* Timing output moved to ado file */
    if (n_keepusing > 0 && keepusing_local[0] &&
        !cmerge_macro_ints(keepusing_local, keepusing_indices, (size_t)n_keepusing)) {
        SF_error("cmerge: invalid keepusing list\n");
        return 198;
    }
    if (n_keepusing > 0 && !keepusing_local[0] && keepusing_idx != n_keepusing) {
        SF_error("cmerge: invalid keepusing list\n");
        return 198;
    }

    /* Acquire lock before modifying global cache */
    cache_lock_acquire();

    /* Check if cache is currently in use by execute */
    if (g_cache_in_use) {
        cache_lock_release();
        SF_error("cmerge: Cannot load using data while merge is in progress\n");
        return 459;
    }

    /* Free any previous cache */
    if (g_using_cache.loaded) {
        ctools_filtered_data_free(&g_using_cache.keys);
        cmerge_free_keepusing();
        g_using_cache.loaded = 0;
    }

    cache_lock_release();

    double t_phase1_start = ctools_timer_ms();

    double t_load_keys = 0;
    double t_load_keepusing = 0;

    /* For merge_by_n, we don't need key variables - just keepusing */
    if (merge_by_n) {
        /* Initialize empty keys structure */
        ctools_filtered_data_init(&g_using_cache.keys);
        g_using_cache.keys.data.nobs = SF_nobs();  /* Get nobs from Stata */

        /* Load keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = ctools_timer_ms();
            rc = cmerge_load_keepusing(keepusing_indices, n_keepusing, strl);
            if (rc != STATA_OK) {
                ctools_filtered_data_free(&g_using_cache.keys);
                cmerge_free_keepusing();
                SF_error("cmerge: Failed to load keepusing variables\n");
                return ctools_stata_rc(rc);
            }
            if (g_using_cache.keepusing.data.vars == NULL) {
                ctools_filtered_data_free(&g_using_cache.keys);
                cmerge_free_keepusing();
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = ctools_timer_ms() - t_load_keepusing_start;
            g_using_cache.keys.data.nobs = g_using_cache.keepusing.data.nobs;
        } else {
            ctools_filtered_data_init(&g_using_cache.keepusing);
        }
    } else {
        /* Standard merge: Load ONLY key variables */
        double t_load_keys_start = ctools_timer_ms();
        rc = ctools_data_load(&g_using_cache.keys, key_indices, nkeys, 0, 0,
                              CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_ROW_TILES);
        if (rc != STATA_OK) {
            ctools_filtered_data_free(&g_using_cache.keys);
            SF_error("cmerge: Failed to load using keys\n");
            return ctools_stata_rc(rc);
        }
        /* Critical null check - prevents crash if load failed silently */
        if (g_using_cache.keys.data.vars == NULL) {
            ctools_filtered_data_free(&g_using_cache.keys);
            SF_error("cmerge: FATAL - keys.vars is NULL after load\n");
            return 920;
        }
        t_load_keys = ctools_timer_ms() - t_load_keys_start;

        /* Load ONLY keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = ctools_timer_ms();
            rc = cmerge_load_keepusing(keepusing_indices, n_keepusing, strl);
            if (rc != STATA_OK) {
                ctools_filtered_data_free(&g_using_cache.keys);
                cmerge_free_keepusing();
                SF_error("cmerge: Failed to load keepusing variables\n");
                return ctools_stata_rc(rc);
            }
            /* Critical null check */
            if (g_using_cache.keepusing.data.vars == NULL) {
                ctools_filtered_data_free(&g_using_cache.keys);
                cmerge_free_keepusing();
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = ctools_timer_ms() - t_load_keepusing_start;
        } else {
            ctools_filtered_data_init(&g_using_cache.keepusing);
        }
    }

    g_using_cache.nobs = g_using_cache.keys.data.nobs;
    g_using_cache.nkeys = nkeys;
    g_using_cache.n_keepusing = n_keepusing;
    g_using_cache.merge_by_n = merge_by_n;

    double t_sort = 0;
    double t_apply_perm = 0;

    /* merge errors when data promised sorted (the option) are not; a Stata
     * sort marker is trusted only after the same check, else the data are
     * sorted below */
    if ((sorted || sorted_hint) && !merge_by_n && nkeys > 0) {
        int key_cols[CMERGE_MAX_KEYVARS];
        int ok = (nkeys <= CMERGE_MAX_KEYVARS);
        for (int k = 0; ok && k < nkeys; k++) key_cols[k] = k;
        ok = ok && cmerge_keys_in_order(&g_using_cache.keys.data, key_cols, nkeys);
        if (!ok && sorted) {
            ctools_filtered_data_free(&g_using_cache.keys);
            cmerge_free_keepusing();
            SF_error("cmerge: using data not sorted\n");
            return 5;
        }
        sorted = ok;
    }

    /* Sort the using keys (skip if sorted option specified or merge_by_n).
     * The dispatcher's algorithms are all stable, so rows with equal keys
     * keep the using data's order, which 1:m and m:m pair by. Only the keys
     * are put in key order; the keepusing columns stay where they are and
     * are read through g_using_cache.perm. */
    if (!sorted && !merge_by_n) {
        double t_sort_start = ctools_timer_ms();
        int sort_vars[CMERGE_MAX_KEYVARS];
        for (int i = 0; i < nkeys; i++) sort_vars[i] = i + 1;
        cmerge_tighten_string_keys(&g_using_cache.keys.data, nkeys);
        stata_retcode sort_rc = ctools_sort_dispatch(&g_using_cache.keys.data, sort_vars,
                                                     (size_t)nkeys, SORT_ALG_AUTO);
        int failed = 0;
        perm_idx_t *order = NULL;
        if (sort_rc == STATA_OK) {
            order = cmerge_take_order(&g_using_cache.keys.data, &failed);
            if (failed) sort_rc = STATA_ERR_MEMORY;
        }
        if (sort_rc == STATA_OK && order)
            sort_rc = cmerge_order_keys(&g_using_cache.keys.data, nkeys, order);
        if (sort_rc != STATA_OK) {
            free(order);
            ctools_filtered_data_free(&g_using_cache.keys);
            cmerge_free_keepusing();
            SF_error("cmerge: Failed to sort using data\n");
            return ctools_stata_rc(sort_rc);
        }
        g_using_cache.perm = order;
        t_sort = ctools_timer_ms() - t_sort_start;
    }

    g_using_cache.loaded = 1;

    double t_phase1_total = ctools_timer_ms() - t_phase1_start;

    /* Save Phase 1 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_using_nobs", (double)g_using_cache.nobs);
    SF_scal_save("_cmerge_p1_load_keys", t_load_keys / 1000.0);
    SF_scal_save("_cmerge_p1_load_keepusing", t_load_keepusing / 1000.0);
    SF_scal_save("_cmerge_p1_sort", t_sort / 1000.0);
    SF_scal_save("_cmerge_p1_apply_perm", t_apply_perm / 1000.0);
    SF_scal_save("_cmerge_p1_total", t_phase1_total / 1000.0);
    CTOOLS_SAVE_THREAD_INFO("_cmerge");

    return 0;
}

/* ============================================================================
 * Phase 2: plan, execute, then commit or roll back
 *
 * plan:     runs in the untouched master. Loads the master keys alone, checks a
 *           promised sort order, sorts and joins the keys, and keeps the output
 *           plan, so key, sortedness and uniqueness errors leave Stata unchanged.
 * execute:  runs after the ado has added the new variables and observations.
 *           Before the first write to an original master cell it copies every
 *           original master variable: staged rows, or full columns when the
 *           strings are too wide to stage. The output is written from that
 *           copy, and a failed write restores the original observations from it
 *           before returning the error. An identity output reuses the master
 *           already in Stata. With nocopy (the ado's snapshot or fast modes) no
 *           copy is kept.
 * rollback: after a successful execute (Break, or a failure in the ado's
 *           remaining steps), rewrites the original observations from the copy.
 * commit:   releases the undo copy once the ado has finished.
 * ============================================================================ */

/* The join, computed by plan and consumed by execute. */
typedef struct {
    int ready;
    cmerge_type_t merge_type;
    int merge_by_n;
    int nkeys;
    int key_idx[CMERGE_MAX_KEYVARS]; /* one-based plugin indices of the master keys */
    size_t master_nobs;
    size_t orig_nvars;               /* master variables when planned (the plan call lists all) */
    size_t output_nobs;
    cmerge_output_spec_t *specs;     /* owned [output_nobs] */
    int32_t *orig_rows;              /* owned [output_nobs]: zero-based master row, -1 using-only */
    size_t n_master_only, n_using_only, n_matched;
    int identity;                    /* the output rows are the master rows, in order */
    int hold_undo;                   /* keep an undo copy even if no original cell is written */
    int nocopy;                      /* keep no undo copy (the ado's snapshot or fast mode) */
    double t_load, t_sort, t_join, t_reorder, t_total; /* milliseconds */
} cmerge_plan_t;

/* Every original master variable, observations 1..master_nobs in their original
 * order, held from before the first write to an original master cell until
 * the ado commits or rolls back: staged rows, or loaded columns when the
 * master cannot be staged. */
enum { CMERGE_UNDO_NONE = 0, CMERGE_UNDO_ROWS, CMERGE_UNDO_COLUMNS };
typedef struct {
    int held;                  /* CMERGE_UNDO_* */
    size_t master_nobs;
    ctools_rowstage rows;
    ctools_filtered_data cols;
    int *col_vars;             /* owned [n_cols]: one-based plugin index of each column */
    size_t n_cols;
} cmerge_undo_t;

static cmerge_plan_t g_plan;
static cmerge_undo_t g_undo;

static void cmerge_plan_free(void)
{
    free(g_plan.specs);
    free(g_plan.orig_rows);
    memset(&g_plan, 0, sizeof(g_plan));
}

static void cmerge_undo_free(void)
{
    ctools_rowstage_free(&g_undo.rows);
    ctools_filtered_data_free(&g_undo.cols);
    free(g_undo.col_vars);
    g_undo.col_vars = NULL;
    g_undo.n_cols = 0;
    g_undo.held = CMERGE_UNDO_NONE;
    g_undo.master_nobs = 0;
}

/* Write the original observations back from the undo copy. The in-order row
 * store allocates nothing. */
static stata_retcode cmerge_undo_restore(void)
{
    if ((size_t)SF_nobs() < g_undo.master_nobs) return STATA_ERR_INVALID_INPUT;
    if (g_undo.held == CMERGE_UNDO_ROWS)
        return ctools_rowstage_store(&g_undo.rows, NULL, g_undo.master_nobs, 1);
    if (g_undo.held == CMERGE_UNDO_COLUMNS)
        return ctools_data_store_selective(&g_undo.cols.data, g_undo.col_vars, g_undo.n_cols, 1);
    return STATA_OK;
}

/* Everything one merge holds between plugin calls. */
static void cmerge_release_all(void)
{
    cmerge_undo_free();
    cmerge_plan_free();
    cache_cleanup_on_error();
}

/* Arguments: "<merge_type> <nkeys> <key_idx>... master_nobs <N>
 * preserve_order <0|1> [sorted] [sorted_hint] [hold_undo] [nocopy]
 * [merge_by_n] [verbose]".
 * The plugin varlist lists every master variable in its current order; key
 * indices are one-based within it, and execute uses the same indices.
 * hold_undo: the ado will drop rows (keep()), so an undo copy is needed even
 * when no original cell is written. nocopy: keep no undo copy at all. */
static ST_retcode cmerge_plan(const char *args)
{
    double t_total_start = ctools_timer_ms();
    cmerge_plan_free();
    cmerge_undo_free();

    char args_copy[4096];
    if (strlen(args) >= sizeof(args_copy)) {
        cache_cleanup_on_error();
        return 198;
    }
    memcpy(args_copy, args, strlen(args) + 1);

    cmerge_plan_t plan;
    memset(&plan, 0, sizeof(plan));
    int preserve_order = 0, sorted = 0, sorted_hint = 0, value = 0, bad = 0, arg_idx = 0;
    int keep_mask = 0, update_split = 0; /* keep(): bit c set when code c is kept */
    for (char *token = strtok(args_copy, " "); token && !bad;
         token = strtok(NULL, " "), arg_idx++) {
        if (arg_idx == 0) {
            bad = !ctools_safe_atoi(token, &value) || value < MERGE_1_1 || value > MERGE_M_M;
            plan.merge_type = (cmerge_type_t)value;
        } else if (arg_idx == 1) {
            bad = !ctools_safe_atoi(token, &plan.nkeys) || plan.nkeys < 0 ||
                  plan.nkeys > CMERGE_MAX_KEYVARS;
        } else if (arg_idx < 2 + plan.nkeys) {
            int *key = &plan.key_idx[arg_idx - 2];
            bad = !ctools_safe_atoi(token, key) || *key < 1 || *key > SF_nvars();
        } else if (strcmp(token, "master_nobs") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atozu(token, &plan.master_nobs);
        } else if (strcmp(token, "preserve_order") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &preserve_order);
        } else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        } else if (strcmp(token, "sorted_hint") == 0) {
            sorted_hint = 1;
        } else if (strcmp(token, "hold_undo") == 0) {
            plan.hold_undo = 1;
        } else if (strcmp(token, "keep_mask") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &keep_mask) || keep_mask < 0 || keep_mask > 62;
        } else if (strcmp(token, "update_split") == 0) {
            update_split = 1;
        } else if (strcmp(token, "nocopy") == 0) {
            plan.nocopy = 1;
        }
        /* merge_by_n is taken from the using cache; verbose needs nothing here */
    }

    cache_lock_acquire();
    int loaded = g_using_cache.loaded && !g_cache_in_use &&
                 (g_using_cache.merge_by_n || g_using_cache.keys.data.vars != NULL);
    cache_lock_release();
    if (!loaded) {
        cache_cleanup_on_error();
        SF_error("cmerge: Using data not loaded. Call load_using first.\n");
        return 459;
    }
    plan.merge_by_n = g_using_cache.merge_by_n;
    if (bad || plan.master_nobs == 0 || plan.master_nobs > (size_t)SF_nobs() ||
        (!plan.merge_by_n && (plan.nkeys < 1 || plan.nkeys != g_using_cache.nkeys))) {
        cache_cleanup_on_error();
        SF_error("cmerge: invalid plan arguments\n");
        return 198;
    }

    const size_t master_nobs = plan.master_nobs;
    ST_retcode rc = 0;
    ctools_filtered_data keys;
    ctools_filtered_data_init(&keys);
    perm_idx_t *sort_perm = NULL;
    cmerge_output_spec_t *specs = NULL;
    int32_t *orig_rows = NULL;
    size_t output_nobs = 0;
    double t_start = ctools_timer_ms();

    if (plan.merge_by_n) {
        /* Observation i of the master meets observation i of the using data. */
        size_t using_nobs = g_using_cache.nobs;
        output_nobs = master_nobs > using_nobs ? master_nobs : using_nobs;
        if (output_nobs > INT32_MAX) {
            SF_error("cmerge: Merge result exceeds Stata observation limit\n");
            rc = 920;
            goto done;
        }
        specs = ctools_safe_malloc2(output_nobs, sizeof(cmerge_output_spec_t));
        if (!specs) {
            rc = 920;
            goto done;
        }
        for (size_t i = 0; i < output_nobs; i++) {
            int in_master = i < master_nobs, in_using = i < using_nobs;
            specs[i].master_sorted_row = in_master ? (int32_t)i : -1;
            specs[i].using_sorted_row = in_using ? (int32_t)i : -1;
            specs[i].merge_result = in_master && in_using ? MERGE_RESULT_BOTH
                                    : in_master          ? MERGE_RESULT_MASTER_ONLY
                                                         : MERGE_RESULT_USING_ONLY;
        }
        plan.t_join = ctools_timer_ms() - t_start;
    } else {
        stata_retcode load_rc = ctools_data_load(&keys, plan.key_idx, (size_t)plan.nkeys, 1,
                                                 master_nobs, CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_ROW_TILES);
        if (load_rc != STATA_OK || keys.data.vars == NULL) {
            SF_error("cmerge: Failed to load master data\n");
            rc = load_rc != STATA_OK ? ctools_stata_rc(load_rc) : 920;
            goto done;
        }
        plan.t_load = ctools_timer_ms() - t_start;

        /* merge errors when data promised sorted (the option) are not; a Stata
         * sort marker is trusted only after the same check, else the data are
         * sorted below */
        int key_cols[CMERGE_MAX_KEYVARS];
        for (int k = 0; k < plan.nkeys; k++) key_cols[k] = k;
        if (sorted || sorted_hint) {
            int ok = cmerge_keys_in_order(&keys.data, key_cols, plan.nkeys);
            if (!ok && sorted) {
                SF_error("cmerge: master data not sorted\n");
                rc = 5;
                goto done;
            }
            sorted = ok;
        }

        if (!sorted) {
            t_start = ctools_timer_ms();
            /* Every algorithm here is stable, so the dispatcher's pick for the
             * key types does not change the order. One-based key positions. */
            int sort_vars[CMERGE_MAX_KEYVARS];
            for (int k = 0; k < plan.nkeys; k++) sort_vars[k] = k + 1;
            cmerge_tighten_string_keys(&keys.data, plan.nkeys);
            stata_retcode sort_rc =
                ctools_sort_dispatch(&keys.data, sort_vars, (size_t)plan.nkeys, SORT_ALG_AUTO);
            if (sort_rc != STATA_OK) {
                rc = ctools_stata_rc(sort_rc);
                goto done;
            }
            /* sort_perm[i]: the master row that is sorted row i */
            int failed = 0;
            sort_perm = cmerge_take_order(&keys.data, &failed);
            if (failed) {
                rc = 920;
                goto done;
            }
            /* Put the keys in sorted order for the join. */
            if (sort_perm && cmerge_order_keys(&keys.data, plan.nkeys, sort_perm) != STATA_OK) {
                rc = 920;
                goto done;
            }
            plan.t_sort = ctools_timer_ms() - t_start;
        }

        t_start = ctools_timer_ms();
        int64_t joined = cmerge_sorted_join(&keys.data, &g_using_cache.keys.data, plan.nkeys,
                                            plan.merge_type, &specs);
        if (joined < 0) {
            specs = NULL;
            if (joined == -4) {
                SF_error("cmerge: key variable types must match\n");
                rc = 106;
            } else if (joined == -2 || joined == -3) {
                SF_error(joined == -2 ?
                    "cmerge: keys do not uniquely identify master observations\n" :
                    "cmerge: keys do not uniquely identify using observations\n");
                rc = 459;
            } else {
                SF_error("cmerge: Merge join failed\n");
                rc = 920;
            }
            goto done;
        }
        /* Stata max obs is 2^31-1; reject if exceeded */
        if (joined > INT32_MAX) {
            SF_error("cmerge: Merge result exceeds Stata observation limit\n");
            rc = 920;
            goto done;
        }
        output_nobs = (size_t)joined;
        plan.t_join = ctools_timer_ms() - t_start;
    }

    /* Original (zero-based) master row of each output row. */
    orig_rows = ctools_safe_malloc2(output_nobs ? output_nobs : 1, sizeof(int32_t));
    if (!orig_rows) {
        rc = 920;
        goto done;
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if (output_nobs >= MIN_OBS_PER_THREAD * 2)
    #endif
    for (size_t i = 0; i < output_nobs; i++) {
        int32_t sorted_row = specs[i].master_sorted_row;
        orig_rows[i] = sorted_row < 0 ? -1
                       : sort_perm   ? (int32_t)sort_perm[sorted_row]
                                     : sorted_row;
    }

    /* merge's order: each master row's first output row, in key order, then
     * every other row (using-only rows and the later 1:m or m:m pairings of a
     * master row) in join order. A stable partition of the join; skipped when
     * every row is a first master row. */
    if (!preserve_order && !plan.merge_by_n && output_nobs > 0) {
        uint8_t *seen = (uint8_t *)calloc(master_nobs ? master_nobs : 1, 1);
        uint8_t *head = (uint8_t *)malloc(output_nobs);
        if (!seen || !head) {
            free(seen);
            free(head);
            rc = 920;
            goto done;
        }
        size_t n_tail = 0;
        for (size_t i = 0; i < output_nobs; i++) {
            int32_t m = specs[i].master_sorted_row;
            head[i] = m >= 0 && !seen[m];
            if (head[i]) seen[m] = 1;
            n_tail += !head[i];
        }
        free(seen);
        /* Rows already in place: the tail is the trailing block */
        size_t first_tail = output_nobs;
        for (size_t i = 0; i < output_nobs; i++) {
            if (!head[i]) { first_tail = i; break; }
        }
        if (n_tail > 0 && output_nobs - first_tail != n_tail) {
            t_start = ctools_timer_ms();
            cmerge_output_spec_t *new_specs =
                ctools_safe_malloc2(output_nobs, sizeof(cmerge_output_spec_t));
            int32_t *new_orig_rows = ctools_safe_malloc2(output_nobs, sizeof(int32_t));
            if (!new_specs || !new_orig_rows) {
                free(head);
                free(new_specs);
                free(new_orig_rows);
                rc = 920;
                goto done;
            }
            size_t h = 0, t = output_nobs - n_tail;
            for (size_t i = 0; i < output_nobs; i++) {
                size_t dest = head[i] ? h++ : t++;
                new_specs[dest] = specs[i];
                new_orig_rows[dest] = orig_rows[i];
            }
            free(specs);
            free(orig_rows);
            specs = new_specs;
            orig_rows = new_orig_rows;
            plan.t_reorder = ctools_timer_ms() - t_start;
        }
        free(head);
    }

    /* preserve_order: master rows in their original order, using-only rows after. */
    if (preserve_order && output_nobs > 0) {
        t_start = ctools_timer_ms();
        ctools_order_pair_t *pairs = ctools_safe_malloc2(output_nobs, sizeof(ctools_order_pair_t));
        cmerge_output_spec_t *new_specs =
            ctools_safe_malloc2(output_nobs, sizeof(cmerge_output_spec_t));
        int32_t *new_orig_rows = ctools_safe_malloc2(output_nobs, sizeof(int32_t));
        if (!pairs || !new_specs || !new_orig_rows) {
            free(pairs);
            free(new_specs);
            free(new_orig_rows);
            rc = 920;
            goto done;
        }
        for (size_t i = 0; i < output_nobs; i++) {
            pairs[i].orig_idx = i;
            pairs[i].order_key = orig_rows[i] >= 0
                ? (size_t)orig_rows[i]
                : master_nobs + (size_t)specs[i].using_sorted_row;
        }
        ctools_sort_order_pairs(pairs, output_nobs);
        for (size_t i = 0; i < output_nobs; i++) {
            new_specs[i] = specs[pairs[i].orig_idx];
            new_orig_rows[i] = orig_rows[pairs[i].orig_idx];
        }
        free(pairs);
        free(specs);
        free(orig_rows);
        specs = new_specs;
        orig_rows = new_orig_rows;
        plan.t_reorder = ctools_timer_ms() - t_start;
    }

    /* Count the results; identity: every output row is master row i, in order. */
    {
        size_t n1 = 0, n2 = 0, n3 = 0;
        int identity = 1;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:n1, n2, n3) reduction(&:identity) if (output_nobs >= MIN_OBS_PER_THREAD * 2)
        #endif
        for (size_t i = 0; i < output_nobs; i++) {
            int8_t m = specs[i].merge_result;
            n1 += (m == MERGE_RESULT_MASTER_ONLY);
            n2 += (m == MERGE_RESULT_USING_ONLY);
            n3 += (m == MERGE_RESULT_BOTH);
            identity &= specs[i].master_sorted_row == (int32_t)i && orig_rows[i] == (int32_t)i;
        }
        plan.n_master_only = n1;
        plan.n_using_only = n2;
        plan.n_matched = n3;
        plan.identity = identity;
    }


done:
    free(sort_perm);
    ctools_filtered_data_free(&keys);
    if (rc) {
        free(specs);
        free(orig_rows);
        cache_cleanup_on_error();
        return rc;
    }
    plan.output_nobs = output_nobs;
    plan.specs = specs;
    plan.orig_rows = orig_rows;
    plan.orig_nvars = (size_t)SF_nvars();
    /* keep() that keeps every row's code drops nothing, so it needs no undo
     * copy beyond what the written cells need. With update, matched rows
     * become 3, 4 or 5 only during execute: then all three must be kept. */
    if (plan.hold_undo && keep_mask) {
        const int matched_bits = update_split ? (1 << 3) | (1 << 4) | (1 << 5) : (1 << 3);
        int drops = (plan.n_master_only && !(keep_mask & (1 << 1))) ||
                    (plan.n_using_only && !(keep_mask & (1 << 2))) ||
                    (plan.n_matched && (keep_mask & matched_bits) != matched_bits);
        if (!drops) plan.hold_undo = 0;
    }
    plan.t_total = ctools_timer_ms() - t_total_start;
    plan.ready = 1;
    g_plan = plan;
    SF_scal_save("_cmerge_plan_N", (double)output_nobs);
    return 0;
}

/* Load every original master variable (one-based plugin indices) as full
 * columns in their original order, for data that cannot be staged. */
static stata_retcode cmerge_load_columns(ctools_filtered_data *master, int *vars, size_t nvars,
                                         size_t master_nobs)
{
    ctools_filtered_data_init(master);
    stata_retcode rc = ctools_data_load(master, vars, nvars, 1, master_nobs,
                                        CTOOLS_LOAD_SKIP_IF);
    if (rc == STATA_OK && master->data.vars == NULL) rc = STATA_ERR_MEMORY;
    if (rc != STATA_OK) ctools_filtered_data_free(master);
    return rc;
}

/* Column path for master data that cannot be staged: gather the output rows
 * from the loaded columns of every original master variable (vars, one-based
 * plugin indices, in the order of master's columns) and store them. */
static stata_retcode cmerge_store_columns(const ctools_filtered_data *master, int *vars,
                                          size_t nvars, double *t_permute)
{
    const size_t output_nobs = g_plan.output_nobs;
    const int32_t *orig_rows = g_plan.orig_rows;
    const cmerge_output_spec_t *specs = g_plan.specs;
    const size_t using_nobs = g_using_cache.nobs;
    stata_retcode rc;

    double t_start = ctools_timer_ms();
    stata_data output;
    stata_data_init(&output);
    size_t vars_bytes;
    if (ctools_safe_mul_size(nvars, sizeof(stata_variable), &vars_bytes) == 0) {
        output.vars = (stata_variable *)ctools_cacheline_alloc(vars_bytes);
    }
    if (!output.vars) {
        return STATA_ERR_MEMORY;
    }
    memset(output.vars, 0, vars_bytes);
    output.nobs = output_nobs;
    output.nvars = nvars;

    /* Output string arrays borrow contents from the loaded master columns and
     * the using keys, which both outlive the store. */
    int alloc_failed = 0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (size_t j = 0; j < nvars; j++) {
        int failed;
        #ifdef _OPENMP
        #pragma omp atomic read
        #endif
        failed = alloc_failed;
        if (failed) continue;

        const stata_variable *src = &master->data.vars[j];
        stata_variable *dst = &output.vars[j];
        const stata_variable *using_key = NULL;
        for (int k = 0; k < g_plan.nkeys; k++) {
            if (g_plan.key_idx[k] == vars[j]) using_key = &g_using_cache.keys.data.vars[k];
        }
        dst->nobs = output_nobs;
        dst->type = src->type;
        if (src->type == STATA_TYPE_DOUBLE) {
            double *values = (double *)ctools_safe_cacheline_alloc2(output_nobs ? output_nobs : 1,
                                                                    sizeof(double));
            if (!values) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }
            for (size_t i = 0; i < output_nobs; i++) {
                int32_t row = orig_rows[i], using_row = specs[i].using_sorted_row;
                values[i] = row >= 0 ? src->data.dbl[row]
                            : (using_key && using_row >= 0 && (size_t)using_row < using_nobs)
                                ? using_key->data.dbl[using_row]
                                : SV_missval;
            }
            dst->data.dbl = values;
        } else {
            dst->buffer_storage = CTOOLS_BUFFER_HEAP;
            dst->string_storage = CTOOLS_STRINGS_BORROWED;
            char **values = (char **)ctools_safe_calloc2(output_nobs ? output_nobs : 1,
                                                         sizeof(char *));
            if (!values) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }
            for (size_t i = 0; i < output_nobs; i++) {
                int32_t row = orig_rows[i], using_row = specs[i].using_sorted_row;
                char *value = row >= 0 ? src->data.str[row]
                              : (using_key && using_row >= 0 && (size_t)using_row < using_nobs)
                                  ? using_key->data.str[using_row]
                                  : NULL;
                values[i] = value ? value : "";
            }
            dst->data.str = values;
        }
    }
    *t_permute = ctools_timer_ms() - t_start;

    rc = alloc_failed ? STATA_ERR_MEMORY
                      : ctools_data_store_selective(&output, vars, nvars, 1);
    /* Output buffers only borrow strings; the loaded columns stay with the caller. */
    stata_data_free(&output);
    return rc;
}

/* An execute that fails after writing restores the original observations from
 * the undo copy before reporting the error: the runtime releases this
 * command's state after any failed call, so a later rollback would find no
 * copy. The ado then undoes its own changes. */
static ST_retcode cmerge_execute_failed(stata_retcode rc)
{
    if (g_undo.held && cmerge_undo_restore() != STATA_OK) {
        SF_error("cmerge: could not restore the master observations\n");
    }
    cmerge_release_all();
    return ctools_stata_rc(rc);
}

/* Everything execute writes, in one pass over the output rows: each worker
 * owns a run of rows and moves it in tiles, writing every variable of a tile
 * before the next, so each Stata observation is touched while it is in cache
 * instead of once per variable. Per tile: the master variables from the
 * staged rows (when rows != NULL), the using keys of using-only rows, the
 * keepusing variables, then _merge. update and replace read the master value
 * a shared variable has in the output, written earlier in the same tile. */
typedef struct {
    const ctools_rowstage *rows;        /* staged master rows to write, or NULL */
    const int32_t *orig_rows;           /* [output_nobs] zero-based master row, -1 none */
    const cmerge_output_spec_t *specs;  /* [output_nobs] */
    size_t output_nobs;
    int nkeys;
    const int *key_idx;                 /* one-based master key variables */
    int n_keepusing;
    const int *ku_dest;                 /* one-based destination of each keepusing column */
    const int *ku_shared;               /* destination is a master variable */
    int update_mode, replace_mode;
    int merge_var_idx;                  /* one-based, 0 none */
} cmerge_writer_t;

#define CMERGE_WRITE_TILE_BYTES ((size_t)64 * 1024)
#define CMERGE_WRITE_TILE_MAX 256
#define CMERGE_WRITE_CELLS_PER_THREAD 65536
#define CMERGE_WRITE_PREFETCH_ROWS 8

/* Keep the first nonzero code */
static void cmerge_record_error(atomic_int *error, int rc)
{
    int expected = 0;
    if (rc) atomic_compare_exchange_strong(error, &expected, rc);
}

static stata_retcode cmerge_write_output(const cmerge_writer_t *w, size_t *n_updated_out,
                                         size_t *n_conflict_out, int *values_changed_out)
{
    const size_t N = w->output_nobs;
    const ctools_rowstage *rs = w->rows;
    const size_t stride = rs ? rs->stride : 0;
    *n_updated_out = *n_conflict_out = 0;
    *values_changed_out = 0;
    if (N == 0) return STATA_OK;

    size_t tile = CMERGE_WRITE_TILE_MAX;
    if (rs && stride && CMERGE_WRITE_TILE_BYTES / stride < tile) tile = CMERGE_WRITE_TILE_BYTES / stride;
    if (tile < 16) tile = 16;

    size_t row_cells = (rs ? rs->nvars : 0) + (size_t)w->n_keepusing + (w->merge_var_idx > 0) + 1;
    int threads = 1;
#ifdef _OPENMP
    {
        size_t cells;
        threads = ctools_get_max_threads();
        if (threads < 1) threads = 1;
        if (ctools_safe_mul_size(N, row_cells, &cells) == 0 &&
            cells / CMERGE_WRITE_CELLS_PER_THREAD < (size_t)threads)
            threads = (int)(cells / CMERGE_WRITE_CELLS_PER_THREAD);
        if (threads < 1) threads = 1;
    }
#endif

    /* A using-only row has no master row: missing numbers, empty strings. */
    char *blank = NULL;
    if (rs) {
        blank = (char *)calloc(1, stride);
        if (!blank) return STATA_ERR_MEMORY;
        for (size_t j = 0; j < rs->nvars; j++)
            if (rs->width[j] == 0) *(double *)(void *)(blank + rs->offset[j]) = SV_missval;
    }

    const stata_variable *ku_cols = g_using_cache.keepusing.data.vars;
    const stata_variable *ukeys = g_using_cache.keys.data.vars;
    const int32_t using_nobs = (int32_t)g_using_cache.nobs;
    const int checking = w->update_mode || w->replace_mode;
    const int counting = w->update_mode && w->n_keepusing > 0;
    atomic_int error = 0;
    size_t n_updated = 0, n_conflict = 0;
    int values_changed = 0;
    (void)threads;

    #ifdef _OPENMP
    #pragma omp parallel num_threads(threads) if (threads > 1) reduction(+:n_updated, n_conflict) reduction(|:values_changed)
    #endif
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        size_t base = N / team, extra = N % team;
        size_t start = tid * base + (tid < extra ? tid : extra);
        size_t end = start + base + (tid < extra);
        char *buf = rs ? (char *)ctools_safe_malloc2(tile, stride) : NULL;
        uint8_t flags[CMERGE_WRITE_TILE_MAX];
        int32_t ku_row[CMERGE_WRITE_TILE_MAX];
        char str_buf[2049];
        if (rs && !buf) cmerge_record_error(&error, STATA_ERR_MEMORY);

        for (size_t t0 = start; t0 < end && (!rs || buf) && !atomic_load(&error); t0 += tile) {
            const size_t cnt = end - t0 < tile ? end - t0 : tile;
            const ST_int obs0 = (ST_int)(t0 + 1);
            int io_rc = 0;

            /* master variables from the staged rows */
            if (rs) {
                for (size_t i = 0; i < cnt; i++) {
                    int32_t row = w->orig_rows[t0 + i];
#if defined(__GNUC__) || defined(__clang__)
                    if (t0 + i + CMERGE_WRITE_PREFETCH_ROWS < end) {
                        int32_t ahead = w->orig_rows[t0 + i + CMERGE_WRITE_PREFETCH_ROWS];
                        if (ahead >= 0)
                            for (size_t b = 0; b < stride && b < 1024; b += CACHE_LINE_SIZE)
                                __builtin_prefetch(rs->rows + (size_t)ahead * stride + b, 0, 0);
                    }
#endif
                    memcpy(buf + i * stride, row < 0 ? blank : rs->rows + (size_t)row * stride, stride);
                }
                for (size_t j = 0; j < rs->nvars && !io_rc; j++) {
                    const ST_int idx = (ST_int)rs->var_idx[j];
                    const char *slot = buf + rs->offset[j];
                    if (rs->width[j] == 0) {
                        for (size_t i = 0; i < cnt && !io_rc; i++, slot += stride)
                            io_rc = (_stata_)->safestore(idx, obs0 + (ST_int)i,
                                                         *(const double *)(const void *)slot);
                    } else {
                        for (size_t i = 0; i < cnt && !io_rc; i++, slot += stride)
                            io_rc = (_stata_)->sstore(idx, obs0 + (ST_int)i, (char *)slot);
                    }
                }
                /* using-only rows take their key values from the using data */
                for (size_t i = 0; i < cnt && !io_rc; i++) {
                    if (w->orig_rows[t0 + i] >= 0) continue;
                    int32_t u = w->specs[t0 + i].using_sorted_row;
                    if (u < 0 || u >= using_nobs) continue;
                    for (int k = 0; k < w->nkeys && !io_rc; k++) {
                        const stata_variable *key = &ukeys[k];
                        ST_int kidx = (ST_int)w->key_idx[k];
                        if (key->type == STATA_TYPE_DOUBLE) {
                            io_rc = (_stata_)->safestore(kidx, obs0 + (ST_int)i, key->data.dbl[u]);
                        } else {
                            char *value = key->data.str[u];
                            io_rc = SF_sstore(kidx, obs0 + (ST_int)i, value ? value : "");
                        }
                    }
                }
            }

            /* keepusing variables */
            if (w->n_keepusing > 0) {
                for (size_t i = 0; i < cnt; i++) {
                    int32_t u = w->specs[t0 + i].using_sorted_row;
                    ku_row[i] = u >= 0 && u < using_nobs ? cmerge_using_row(u) : -1;
                    flags[i] = 0;
                }
            }
            for (int kv = 0; kv < w->n_keepusing && !io_rc; kv++) {
                const stata_variable *src = &ku_cols[kv];
                const ST_int dest = (ST_int)w->ku_dest[kv];
                const int shared = w->ku_shared[kv];
#if defined(__GNUC__) || defined(__clang__)
                /* unsorted using rows are read at random: fetch ahead (a
                 * string's pointer first, its bytes once the pointer is in) */
                for (size_t i = 0; i < cnt && i < 16; i++)
                    if (ku_row[i] >= 0)
                        __builtin_prefetch(src->type == STATA_TYPE_DOUBLE ? (const void *)&src->data.dbl[ku_row[i]]
                                                                          : (const void *)&src->data.str[ku_row[i]], 0, 0);
#endif
                for (size_t i = 0; i < cnt && !io_rc; i++) {
                    const int8_t result = w->specs[t0 + i].merge_result;
                    const int32_t r = ku_row[i];
                    const ST_int obs = obs0 + (ST_int)i;
#if defined(__GNUC__) || defined(__clang__)
                    if (i + 16 < cnt && ku_row[i + 16] >= 0)
                        __builtin_prefetch(src->type == STATA_TYPE_DOUBLE ? (const void *)&src->data.dbl[ku_row[i + 16]]
                                                                          : (const void *)&src->data.str[ku_row[i + 16]], 0, 0);
                    if (src->type != STATA_TYPE_DOUBLE && i + 8 < cnt && ku_row[i + 8] >= 0 &&
                        src->data.str[ku_row[i + 8]])
                        __builtin_prefetch(src->data.str[ku_row[i + 8]], 0, 0);
#endif
                    if (r < 0 || result == MERGE_RESULT_MASTER_ONLY) continue;
                    int write = result == MERGE_RESULT_USING_ONLY || !shared;
                    if (!write && result == MERGE_RESULT_BOTH && checking) {
                        /* update fills missing master values; replace also
                         * overwrites with nonmissing using values. merge counts
                         * a missing master value as updated whenever the using
                         * value differs, extended missing included. */
                        if (src->type == STATA_TYPE_DOUBLE) {
                            double cur;
                            if ((io_rc = SF_vdata(dest, obs, &cur))) break;
                            double uval = src->data.dbl[r];
                            int mmiss = SF_is_missing(cur), umiss = SF_is_missing(uval);
                            write = w->replace_mode ? (!umiss || mmiss) : mmiss;
                            if (mmiss) {
                                if (uval != cur) flags[i] |= CMERGE_ROW_UPDATED;
                            } else if (!umiss && cur != uval) {
                                flags[i] |= CMERGE_ROW_CONFLICT;
                            }
                        } else {
                            if ((io_rc = SF_sdata(dest, obs, str_buf))) break;
                            const char *uval = src->data.str[r] ? src->data.str[r] : "";
                            int mmiss = str_buf[0] == '\0', umiss = uval[0] == '\0';
                            write = w->replace_mode ? (!umiss || mmiss) : mmiss;
                            if (mmiss && !umiss) flags[i] |= CMERGE_ROW_UPDATED;
                            else if (!mmiss && !umiss && strcmp(str_buf, uval) != 0)
                                flags[i] |= CMERGE_ROW_CONFLICT;
                        }
                    }
                    if (!write) continue;
                    if (src->type == STATA_TYPE_DOUBLE) {
                        io_rc = (_stata_)->safestore(dest, obs, src->data.dbl[r]);
                    } else {
                        const char *val = src->data.str[r] ? src->data.str[r] : "";
                        io_rc = SF_sstore(dest, obs, (char *)val);
                    }
                }
            }

            /* _merge, with update's codes 4 and 5 */
            for (size_t i = 0; i < cnt && !io_rc; i++) {
                const int8_t result = w->specs[t0 + i].merge_result;
                double code = (double)result;
                if (counting && result == MERGE_RESULT_BOTH) {
                    uint8_t f = flags[i];
                    if (f & CMERGE_ROW_CONFLICT) {
                        n_conflict++;
                        code = 5.0;
                    } else if (f & CMERGE_ROW_UPDATED) {
                        n_updated++;
                        code = 4.0;
                    }
                    if ((f & CMERGE_ROW_UPDATED) || (w->replace_mode && (f & CMERGE_ROW_CONFLICT)))
                        values_changed = 1;
                }
                if (w->merge_var_idx > 0)
                    io_rc = (_stata_)->safestore((ST_int)w->merge_var_idx, obs0 + (ST_int)i, code);
            }
            if (io_rc) cmerge_record_error(&error, ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
        }
        free(buf);
    }
    free(blank);
    ctools_memory_barrier();
    stata_retcode rc = (stata_retcode)atomic_load(&error);
    *n_updated_out = n_updated;
    *n_conflict_out = n_conflict;
    *values_changed_out = values_changed;
    return rc;
}

/* Arguments: "<merge_type> <nkeys> <key_idx>... master_nobs <N> master_nvars <K>
 * n_keepusing <n> [keepusing_placeholders <idx>...] [merge_var_idx <idx>]
 * preserve_order <0|1> [update] [replace] [shared_flags <0|1>...] [fault_store]
 * [placeholders_local <name>] [shared_local <name>]", where the two locals hold
 * the placeholder and shared-flag lists in place of the inline lists
 * [sorted|sorted_hint] [merge_by_n] [verbose]".
 * The varlist lists every current variable; the original master variables come
 * first, in the order plan saw them. fault_store is a validation hook: the
 * store fails after writing half of the output rows. */

/* ============================================================================
 * Option tempdta: the output as a .dta
 *
 * The SPI cannot write strL values, so with strL variables in the master or
 * kept from the using data the ado has execute write the whole output to the
 * .dta named by local __cmerge_dta_path and loads it with use. The plugin
 * varlist is every current variable (the master's, then the appended new
 * variables and _merge); local __cmerge_dta_types holds their storage types
 * and __cmerge_dta_names their names, in order. Master columns are gathered
 * through the output's master rows and new keepusing columns through its
 * using rows; keys, shared variables and _merge are assembled first, by the
 * rules the SPI writers apply. Stata is not modified, so no undo copy is
 * needed: a failure leaves the master as it was. The ado restores the
 * characteristics and value labels after use.
 * ============================================================================ */

enum { CM_DTA_MASTER = 0, CM_DTA_KEY, CM_DTA_SHARED, CM_DTA_NEW, CM_DTA_MERGE };

/* A string cell of a loaded column: char ** or whole strL spans. */
static inline void cm_dta_str(const stata_variable *var, const ctools_dta_span *spans, int64_t row,
                              const char **p, size_t *n, int *binary)
{
    *binary = 0;
    if (row < 0) {
        *p = "";
        *n = 0;
    } else if (spans) {
        *p = spans[row].ptr ? spans[row].ptr : "";
        *n = spans[row].ptr ? spans[row].len : 0;
        *binary = spans[row].binary;
    } else {
        const char *s = var->data.str[row];
        *p = s ? s : "";
        *n = s ? strlen(s) : 0;
    }
}

static ctools_dta_type cm_dta_parse_type(const char *t, int *width)
{
    *width = 0;
    if (!strcmp(t, "byte")) return CTOOLS_DTA_BYTE;
    if (!strcmp(t, "int")) return CTOOLS_DTA_INT;
    if (!strcmp(t, "long")) return CTOOLS_DTA_LONG;
    if (!strcmp(t, "float")) return CTOOLS_DTA_FLOAT;
    if (!strcmp(t, "double")) return CTOOLS_DTA_DOUBLE;
    if (!strcmp(t, "strL")) return CTOOLS_DTA_STRL;
    if (!strncmp(t, "str", 3) && ctools_safe_atoi(t + 3, width) && *width >= 1 && *width <= 2045)
        return CTOOLS_DTA_STR;
    return (ctools_dta_type)0;
}

/* A local into a malloc'd buffer (empty when undefined). */
static char *cm_dta_local(const char *name)
{
    char macro[64];
    snprintf(macro, sizeof(macro), "_%s", name);
    for (size_t cap = (size_t)1 << 16; cap <= ((size_t)1 << 26); cap *= 4) {
        char *buf = calloc(cap, 1);
        if (!buf) return NULL;
        if (SF_macro_use(macro, buf, (int)cap) != 0) buf[0] = 0;
        if (strlen(buf) < cap - 1) return buf;
        free(buf);
    }
    return NULL;
}

static size_t cm_dta_tokens(char *s, char **tok, size_t n)
{
    size_t k = 0;
    for (char *p = s; *p;) {
        while (*p == ' ') p++;
        if (!*p) break;
        char *e = p;
        while (*e && *e != ' ') e++;
        if (k < n) tok[k] = p;
        k++;
        if (*e) *e++ = 0;
        p = e;
    }
    return k;
}

static ST_retcode cmerge_execute_dta(size_t K, const int *master_key_indices, int nkeys,
                                     int n_keepusing, const int *placeholders, const int *shared_flags,
                                     int merge_var_idx, int update_mode, int replace_mode, int fault_store,
                                     double t_total_start)
{
    const size_t N = g_plan.output_nobs, master_nobs = g_plan.master_nobs, orig_nvars = g_plan.orig_nvars;
    const cmerge_output_spec_t *specs = g_plan.specs;
    const int32_t *orig_rows = g_plan.orig_rows;
    const size_t using_nobs = g_using_cache.nobs;
    const double miss = SV_missval;
    ST_retcode rc = 0;
    double t_start = ctools_timer_ms(), t_load = 0, t_write = 0;

    char *path = cm_dta_local("__cmerge_dta_path"), *types = cm_dta_local("__cmerge_dta_types"),
         *names = cm_dta_local("__cmerge_dta_names"), *sort = cm_dta_local("__cmerge_dta_sort");
    char **ttok = calloc(K + 1, sizeof(char *)), **ntok = calloc(K + 1, sizeof(char *)),
         **stok = calloc(K + 1, sizeof(char *));
    int *role = calloc(K + 1, sizeof(int)), *which = calloc(K + 1, sizeof(int));
    int *sortlist = calloc(K + 1, sizeof(int));
    ctools_dta_column *cols = calloc(K + 1, sizeof(*cols));
    const stata_variable **mvar = calloc(K + 1, sizeof(*mvar));
    ctools_dta_span **mspans = calloc(K + 1, sizeof(*mspans));
    char **marena = calloc(K + 1, sizeof(*marena));
    void **built = calloc(K + 1, sizeof(void *)); /* columns assembled here */
    int *plain = calloc(K + 1, sizeof(int)), *plain_at = calloc(K + 1, sizeof(int));
    int64_t *mrow = ctools_safe_malloc2(N ? N : 1, sizeof(int64_t));
    int64_t *urow = ctools_safe_malloc2(N ? N : 1, sizeof(int64_t));
    uint8_t *row_flags = update_mode && N ? calloc(N, 1) : NULL;
    ctools_filtered_data mcols;
    ctools_filtered_data_init(&mcols);
    if (!path || !types || !names || !sort || !ttok || !ntok || !stok || !role || !which || !sortlist ||
        !cols || !mvar || !mspans || !marena || !built || !plain || !plain_at || !mrow || !urow ||
        (update_mode && N && !row_flags)) {
        rc = 920;
        goto done;
    }
    if (!*path || cm_dta_tokens(types, ttok, K) != K || cm_dta_tokens(names, ntok, K) != K) {
        SF_error("cmerge: invalid tempdta arguments\n");
        rc = 198;
        goto done;
    }
    size_t nsort = cm_dta_tokens(sort, stok, K);

    /* What each variable is in the output */
    for (int kv = 0; kv < n_keepusing && !rc; kv++) {
        int d = placeholders[kv];
        if (d < 1 || (size_t)d > K) rc = 198;
        else {
            role[d - 1] = shared_flags[kv] ? CM_DTA_SHARED : CM_DTA_NEW;
            which[d - 1] = kv;
        }
    }
    if (merge_var_idx > 0) role[merge_var_idx - 1] = CM_DTA_MERGE;
    for (int k = 0; k < nkeys && !rc; k++) {
        int d = master_key_indices[k];
        if (d < 1 || (size_t)d > K) rc = 198;
        else if (role[d - 1] == CM_DTA_MASTER) {
            role[d - 1] = CM_DTA_KEY;
            which[d - 1] = k;
        }
    }
    for (size_t j = 0; j < K && !rc; j++) {
        int from_master = role[j] == CM_DTA_MASTER || role[j] == CM_DTA_KEY || role[j] == CM_DTA_SHARED;
        if (from_master != (j < orig_nvars)) rc = 198;
    }
    if (rc) {
        SF_error("cmerge: invalid tempdta variable roles\n");
        goto done;
    }

    /* The master columns in their original order: strL whole, the rest
     * through the parallel loader */
    size_t nplain = 0;
    for (size_t j = 0; j < orig_nvars && !rc; j++) {
        if (!SF_var_is_strl((int)j + 1)) {
            plain_at[nplain] = (int)j;
            plain[nplain++] = (int)j + 1;
        } else if (role[j] == CM_DTA_KEY) {
            rc = 109; /* merge refuses strL keys before this */
        } else if (master_nobs) {
            stata_retcode lrc = cmerge_read_strl((int)j + 1, 1, master_nobs, &mspans[j], &marena[j]);
            if (lrc != STATA_OK) rc = ctools_stata_rc(lrc);
        }
    }
    if (!rc && nplain && master_nobs) {
        stata_retcode lrc = ctools_data_load(&mcols, plain, nplain, 1, master_nobs, CTOOLS_LOAD_SKIP_IF);
        if (lrc != STATA_OK || !mcols.data.vars) rc = lrc != STATA_OK ? ctools_stata_rc(lrc) : 920;
        for (size_t p = 0; p < nplain && !rc; p++) mvar[plain_at[p]] = &mcols.data.vars[p];
    }
    t_load = ctools_timer_ms() - t_start;
    if (rc) goto done;

    for (size_t i = 0; i < N; i++) {
        mrow[i] = orig_rows[i];
        int32_t u = specs[i].using_sorted_row;
        /* keepusing columns are in the using data's own order */
        urow[i] = specs[i].merge_result != MERGE_RESULT_MASTER_ONLY && u >= 0 && (size_t)u < using_nobs
                      ? cmerge_using_row(u) : -1;
    }

    /* Columns assembled from master and using rows */
    size_t values_changed = 0;
    for (size_t j = 0; j < K && !rc; j++) {
        int width;
        ctools_dta_type type = cm_dta_parse_type(ttok[j], &width);
        int string = type == CTOOLS_DTA_STR || type == CTOOLS_DTA_STRL;
        if (!type) {
            rc = 198;
            break;
        }
        const stata_variable *mv = j < orig_nvars ? mvar[j] : NULL;
        const ctools_dta_span *ms = j < orig_nvars ? mspans[j] : NULL;
        if (j < orig_nvars && master_nobs && !mv && !ms) {
            rc = 459;
            break;
        }
        if (role[j] == CM_DTA_KEY) {
            const stata_variable *uk = &g_using_cache.keys.data.vars[which[j]];
            if (!string) {
                double *v = ctools_safe_malloc2(N ? N : 1, sizeof(double));
                if (!(built[j] = v)) { rc = 920; break; }
                for (size_t i = 0; i < N; i++) {
                    int32_t u = specs[i].using_sorted_row;
                    v[i] = mrow[i] >= 0 ? mv->data.dbl[mrow[i]]
                           : u >= 0 && (size_t)u < using_nobs ? uk->data.dbl[u] : miss;
                }
            } else {
                const char **v = ctools_safe_malloc2(N ? N : 1, sizeof(char *));
                if (!(built[j] = (void *)v)) { rc = 920; break; }
                for (size_t i = 0; i < N; i++) {
                    int32_t u = specs[i].using_sorted_row;
                    const char *s = mrow[i] >= 0 ? mv->data.str[mrow[i]]
                                    : u >= 0 && (size_t)u < using_nobs ? uk->data.str[u] : "";
                    v[i] = s ? s : "";
                }
            }
        } else if (role[j] == CM_DTA_SHARED) {
            const int kv = which[j];
            const stata_variable *uv = &g_using_cache.keepusing.data.vars[kv];
            const ctools_dta_span *us = g_using_cache.ku_spans ? g_using_cache.ku_spans[kv] : NULL;
            if (!string) {
                double *v = ctools_safe_malloc2(N ? N : 1, sizeof(double));
                if (!(built[j] = v)) { rc = 920; break; }
                for (size_t i = 0; i < N; i++) {
                    double cur = mrow[i] >= 0 ? mv->data.dbl[mrow[i]] : miss;
                    int write = specs[i].merge_result == MERGE_RESULT_USING_ONLY;
                    if (specs[i].merge_result == MERGE_RESULT_BOTH && (update_mode || replace_mode)) {
                        double uval = urow[i] >= 0 ? uv->data.dbl[urow[i]] : miss;
                        int mmiss = !(cur < miss), umiss = !(uval < miss);
                        write = replace_mode ? (!umiss || mmiss) : mmiss;
                        uint8_t flag = 0;
                        if (mmiss) {
                            if (uval != cur) flag = CMERGE_ROW_UPDATED;
                        } else if (!umiss && cur != uval) {
                            flag = CMERGE_ROW_CONFLICT;
                        }
                        if (row_flags) row_flags[i] |= flag;
                    }
                    v[i] = write && urow[i] >= 0 ? uv->data.dbl[urow[i]] : cur;
                }
            } else {
                ctools_dta_span *v = ctools_safe_malloc2(N ? N : 1, sizeof(ctools_dta_span));
                if (!(built[j] = v)) { rc = 920; break; }
                for (size_t i = 0; i < N; i++) {
                    const char *mp, *up;
                    size_t mn, un;
                    int mb, ub;
                    cm_dta_str(mv, ms, mrow[i], &mp, &mn, &mb);
                    int write = specs[i].merge_result == MERGE_RESULT_USING_ONLY;
                    if (specs[i].merge_result == MERGE_RESULT_BOTH && (update_mode || replace_mode)) {
                        cm_dta_str(uv, us, urow[i], &up, &un, &ub);
                        int mmiss = mn == 0, umiss = un == 0;
                        write = replace_mode ? (!umiss || mmiss) : mmiss;
                        uint8_t flag = 0;
                        if (mmiss && !umiss) flag = CMERGE_ROW_UPDATED;
                        else if (!mmiss && !umiss && (mn != un || memcmp(mp, up, mn))) flag = CMERGE_ROW_CONFLICT;
                        if (row_flags) row_flags[i] |= flag;
                    }
                    if (write && urow[i] >= 0) {
                        cm_dta_str(uv, us, urow[i], &up, &un, &ub);
                        v[i] = (ctools_dta_span){up, un, ub};
                    } else {
                        v[i] = (ctools_dta_span){mp, mn, mb};
                    }
                }
            }
        }
    }
    if (rc) goto done;

    /* _merge, with update's codes 4 and 5 */
    size_t n_updated = 0, n_conflict = 0;
    if (row_flags) {
        for (size_t i = 0; i < N; i++) {
            if (specs[i].merge_result != MERGE_RESULT_BOTH) continue;
            if (row_flags[i] & CMERGE_ROW_CONFLICT) n_conflict++;
            else if (row_flags[i] & CMERGE_ROW_UPDATED) n_updated++;
            if ((row_flags[i] & CMERGE_ROW_UPDATED) || (replace_mode && (row_flags[i] & CMERGE_ROW_CONFLICT)))
                values_changed = 1;
        }
    }
    for (size_t j = 0; j < K && !rc; j++) {
        if (role[j] != CM_DTA_MERGE) continue;
        double *v = ctools_safe_malloc2(N ? N : 1, sizeof(double));
        if (!(built[j] = v)) { rc = 920; break; }
        for (size_t i = 0; i < N; i++) {
            double code = (double)specs[i].merge_result;
            if (row_flags && specs[i].merge_result == MERGE_RESULT_BOTH) {
                if (row_flags[i] & CMERGE_ROW_CONFLICT) code = 5.0;
                else if (row_flags[i] & CMERGE_ROW_UPDATED) code = 4.0;
            }
            v[i] = code;
        }
    }
    if (rc) goto done;

    /* The file's columns */
    for (size_t j = 0; j < K && !rc; j++) {
        int width;
        ctools_dta_type type = cm_dta_parse_type(ttok[j], &width);
        int string = type == CTOOLS_DTA_STR || type == CTOOLS_DTA_STRL;
        ctools_dta_column *c = &cols[j];
        c->name = ntok[j];
        c->type = type;
        c->width = width;
        if (built[j]) {
            c->data = built[j];
            c->source = !string ? CTOOLS_DTA_SRC_F64
                        : role[j] == CM_DTA_KEY ? CTOOLS_DTA_SRC_CSTR : CTOOLS_DTA_SRC_SPAN;
        } else if (role[j] == CM_DTA_NEW) {
            const int kv = which[j];
            const stata_variable *uv = &g_using_cache.keepusing.data.vars[kv];
            const ctools_dta_span *us = g_using_cache.ku_spans ? g_using_cache.ku_spans[kv] : NULL;
            c->rows = urow;
            c->flags = CTOOLS_DTA_COL_SHARE; /* m:1 repeats using rows */
            if (us) c->source = CTOOLS_DTA_SRC_SPAN, c->data = us;
            else if (uv->type == STATA_TYPE_DOUBLE) c->source = CTOOLS_DTA_SRC_F64, c->data = uv->data.dbl;
            else c->source = CTOOLS_DTA_SRC_CSTR, c->data = uv->data.str;
        } else { /* a master column */
            c->rows = mrow;
            c->flags = CTOOLS_DTA_COL_SHARE; /* 1:m repeats master rows */
            if (mspans[j]) c->source = CTOOLS_DTA_SRC_SPAN, c->data = mspans[j];
            else if (mvar[j] && mvar[j]->type == STATA_TYPE_DOUBLE) c->source = CTOOLS_DTA_SRC_F64, c->data = mvar[j]->data.dbl;
            else if (mvar[j]) c->source = CTOOLS_DTA_SRC_CSTR, c->data = mvar[j]->data.str;
            else c->source = string ? CTOOLS_DTA_SRC_CSTR : CTOOLS_DTA_SRC_F64, c->data = NULL;
        }
        if (string != (c->source != CTOOLS_DTA_SRC_F64)) {
            SF_error("cmerge: storage types do not match the merged data\n");
            rc = 198;
        }
        if (!c->data && N) rc = 459;
    }
    if (rc) goto done;

    int same_order = g_plan.identity && g_plan.n_using_only == 0 && N == master_nobs && !values_changed;
    size_t ns = 0;
    for (size_t s = 0; same_order && s < nsort && s < K; s++) {
        int p;
        if (!ctools_safe_atoi(stok[s], &p) || p < 1 || (size_t)p > K) break;
        sortlist[ns++] = p;
    }
    ctools_dta_spec spec;
    ctools_dta_spec_init(&spec);
    spec.nobs = N;
    spec.ncols = K;
    spec.cols = cols;
    spec.sortlist = ns ? sortlist : NULL;
    spec.nsort = ns;
    spec.threads = ctools_get_max_threads();
    spec.temporary = 1;
    ctools_dta_result res;
    double t_w = ctools_timer_ms();
    rc = ctools_dta_write(path, &spec, &res);
    t_write = ctools_timer_ms() - t_w;
    if (rc && res.message[0]) {
        char line[320];
        snprintf(line, sizeof(line), "cmerge: %s\n", res.message);
        SF_error(line);
    }
    if (!rc && fault_store) { /* validation hook: the output could not be written */
        ctools_remove(path);
        rc = 920;
    }
    if (rc) goto done;

    double t_total = g_plan.t_total + ctools_timer_ms() - t_total_start;
    SF_scal_save("_cmerge_N", (double)N);
    SF_scal_save("_cmerge_N_1", (double)g_plan.n_master_only);
    SF_scal_save("_cmerge_N_2", (double)g_plan.n_using_only);
    SF_scal_save("_cmerge_N_3", (double)(g_plan.n_matched - n_updated - n_conflict));
    SF_scal_save("_cmerge_N_4", (double)n_updated);
    SF_scal_save("_cmerge_N_5", (double)n_conflict);
    SF_scal_save("_cmerge_same_order", (double)same_order);
    SF_scal_save("_cmerge_p2_load_master", (g_plan.t_load + t_load) / 1000.0);
    SF_scal_save("_cmerge_p2_sort_master", g_plan.t_sort / 1000.0);
    SF_scal_save("_cmerge_p2_merge_join", g_plan.t_join / 1000.0);
    SF_scal_save("_cmerge_p2_reorder", g_plan.t_reorder / 1000.0);
    SF_scal_save("_cmerge_p2_permute", 0.0);
    SF_scal_save("_cmerge_p2_stream", 0.0);
    SF_scal_save("_cmerge_p2_store", t_write / 1000.0);
    SF_scal_save("_cmerge_p2_write_meta", 0.0);
    SF_scal_save("_cmerge_p2_cleanup", 0.0);
    SF_scal_save("_cmerge_p2_total", t_total / 1000.0);
    SF_scal_save("_cmerge_n_output_vars", (double)K);
    CTOOLS_SAVE_THREAD_INFO("_cmerge");

done:
    for (size_t j = 0; built && j < K; j++) free(built[j]);
    for (size_t j = 0; mspans && j < K; j++) free(mspans[j]);
    for (size_t j = 0; marena && j < K; j++) free(marena[j]);
    ctools_filtered_data_free(&mcols);
    free(path), free(types), free(names), free(sort), free(ttok), free(ntok), free(stok);
    free(role), free(which), free(sortlist), free(cols), free(mvar), free(mspans), free(marena);
    free(built), free(plain), free(plain_at), free(mrow), free(urow), free(row_flags);
    /* The output exists (or failed); nothing is held for an undo. */
    cmerge_release_all();
    return rc;
}

static ST_retcode cmerge_execute(const char *args)
{
    double t_start, t_total_start;

    t_total_start = ctools_timer_ms();

    /* Acquire lock and mark cache as in-use to prevent concurrent clear */
    cache_lock_acquire();

    if (!g_using_cache.loaded || !g_plan.ready) {
        cache_lock_release();
        cmerge_release_all();
        SF_error("cmerge: Using data not loaded or join not planned.\n");
        return 459;
    }

    /* Mark cache as in-use before releasing lock */
    g_cache_in_use = 1;
    cache_lock_release();

    /* Parse arguments - use static buffers to avoid stack overflow.
     * Thread safety: Stata plugins are single-threaded, so these statics
     * cannot race. They're reinitialized on each call. */
    static char args_copy[8192];
    if (strlen(args) >= sizeof(args_copy)) {
        cmerge_release_all();
        return 198;
    }
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    cmerge_type_t merge_type = MERGE_1_1;
    int nkeys = 0;
    static int master_key_indices[CMERGE_MAX_KEYVARS];
    size_t master_nobs = 0;
    size_t master_nvars = 0;
    int n_keepusing = 0;
    int *keepusing_placeholder_indices = NULL;
    int merge_var_idx = 0;
    int update_mode = 0;
    int replace_mode = 0;
    int fault_store = 0;
    int dta = 0;           /* option tempdta: write the output as a .dta */
    int *shared_flags = NULL;
    char placeholders_local[32] = {0};
    char shared_local[32] = {0};

    int arg_idx = 0;
    int in_keepusing_placeholders = 0;
    int keepusing_idx = 0;
    int bad = 0;

    char *token = strtok(args_copy, " ");
    int parse_tmp;  /* Temporary for safe parsing */
    while (token != NULL && !bad) {
        if (arg_idx == 0) {
            bad = !ctools_safe_atoi(token, &parse_tmp);
            merge_type = (cmerge_type_t)parse_tmp;
        }
        else if (arg_idx == 1) {
            bad = !ctools_safe_atoi(token, &nkeys) || nkeys < 0 || nkeys > CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx >= 2 && arg_idx < 2 + nkeys) {
            bad = !ctools_safe_atoi(token, &master_key_indices[arg_idx - 2]);
        }
        else if (strcmp(token, "master_nobs") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atozu(token, &master_nobs);
        }
        else if (strcmp(token, "master_nvars") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atozu(token, &master_nvars) ||
                  master_nvars > (size_t)SF_nvars();
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &n_keepusing) || n_keepusing < 0 ||
                  n_keepusing > SF_nvars();
            if (!bad) {
                keepusing_placeholder_indices = cmerge_grow(
                    &g_cm_placeholders, &g_cm_placeholders_cap, (size_t)n_keepusing);
                shared_flags = cmerge_grow(&g_cm_shared, &g_cm_shared_cap, (size_t)n_keepusing);
                bad = !keepusing_placeholder_indices || !shared_flags;
                if (!bad) memset(shared_flags, 0, (size_t)n_keepusing * sizeof(int));
            }
        }
        else if (strcmp(token, "placeholders_local") == 0 || strcmp(token, "shared_local") == 0) {
            char *dest = token[0] == 'p' ? placeholders_local : shared_local;
            token = strtok(NULL, " ");
            bad = !token || strlen(token) >= sizeof(placeholders_local);
            if (!bad) strcpy(dest, token);
        }
        else if (strcmp(token, "keepusing_placeholders") == 0) {
            in_keepusing_placeholders = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_placeholders && keepusing_placeholder_indices &&
                 keepusing_idx < n_keepusing) {
            bad = !ctools_safe_atoi(token, &keepusing_placeholder_indices[keepusing_idx]);
            keepusing_idx++;
        }
        else if (strcmp(token, "merge_var_idx") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &merge_var_idx);
            in_keepusing_placeholders = 0;
        }
        else if (strcmp(token, "preserve_order") == 0) {
            /* Applied by plan */
            token = strtok(NULL, " ");
            bad = !token;
        }
        else if (strcmp(token, "update") == 0) {
            update_mode = 1;
        }
        else if (strcmp(token, "replace") == 0) {
            replace_mode = 1;
        }
        else if (strcmp(token, "fault_store") == 0) {
            fault_store = 1;
        }
        else if (strcmp(token, "dta") == 0) {
            dta = 1;
        }
        else if (strcmp(token, "shared_flags") == 0) {
            bad = !shared_flags && n_keepusing > 0;
            for (int i = 0; i < n_keepusing && !bad; i++) {
                token = strtok(NULL, " ");
                if (token) {
                    bad = !ctools_safe_atoi(token, &shared_flags[i]);
                } else {
                    shared_flags[i] = 0;
                }
            }
        }
        /* sorted, sorted_hint, merge_by_n and verbose were applied by plan */
        arg_idx++;
        token = strtok(NULL, " ");
    }

    if (!bad && n_keepusing > 0 && placeholders_local[0]) {
        bad = !cmerge_macro_ints(placeholders_local, keepusing_placeholder_indices,
                                 (size_t)n_keepusing);
    }
    if (!bad && n_keepusing > 0 && shared_local[0]) {
        bad = !cmerge_macro_ints(shared_local, shared_flags, (size_t)n_keepusing);
    }
    if (!bad && n_keepusing > 0 && !placeholders_local[0] && keepusing_idx != n_keepusing) {
        bad = 1;
    }

    /* The plan must describe this master. */
    int matches = !bad && merge_type == g_plan.merge_type && nkeys == g_plan.nkeys &&
                  master_nobs == g_plan.master_nobs && master_nobs <= (size_t)SF_nobs() &&
                  (dta || g_plan.output_nobs <= (size_t)SF_nobs()) &&
                  (dta || !g_using_cache.has_strl);
    for (int k = 0; matches && k < nkeys; k++) {
        matches = master_key_indices[k] == g_plan.key_idx[k];
    }
    /* _merge must be one of the variables appended after the master's own,
     * or the merge codes would overwrite a master variable. */
    if (merge_var_idx != 0 &&
        (merge_var_idx < 0 || (size_t)merge_var_idx <= g_plan.orig_nvars ||
         merge_var_idx > SF_nvars() || SF_var_is_string(merge_var_idx))) {
        matches = 0;
    }
    if (!matches) {
        cmerge_release_all();
        SF_error("cmerge: invalid execute arguments\n");
        return 198;
    }

    if (dta) {
        return cmerge_execute_dta(master_nvars, master_key_indices, nkeys, n_keepusing,
                                  keepusing_placeholder_indices, shared_flags, merge_var_idx,
                                  update_mode, replace_mode, fault_store, t_total_start);
    }

    const size_t output_nobs = g_plan.output_nobs;
    cmerge_output_spec_t *output_specs = g_plan.specs;
    const int32_t *master_orig_rows = g_plan.orig_rows;
    const size_t n_using_only = g_plan.n_using_only;

    /* Original master variables: everything except the new keepusing
     * placeholders and _merge, which are all missing here and are filled by
     * the writers below. One-based plugin variable indices. */
    int *stage_vars = NULL;
    size_t n_stage_vars = 0;
    if (master_nvars > 0) {
        stage_vars = ctools_safe_malloc2(master_nvars, sizeof(int));
        if (!stage_vars) {
            cmerge_release_all();
            SF_error("cmerge: Failed to allocate var indices\n");
            return 920;
        }
        for (size_t v = 0; v < master_nvars; v++) {
            int stata_idx = (int)(v + 1);
            int placeholder = (stata_idx == merge_var_idx);
            for (int kv = 0; kv < n_keepusing && !placeholder; kv++) {
                placeholder = (keepusing_placeholder_indices[kv] == stata_idx && !shared_flags[kv]);
            }
            for (int k = 0; k < nkeys; k++) {
                if (master_key_indices[k] == stata_idx) placeholder = 0;
            }
            if (!placeholder) stage_vars[n_stage_vars++] = stata_idx;
        }
    }
    int stage_ok = n_stage_vars > 0 &&
                   ctools_rowstage_supported(stage_vars, n_stage_vars, NULL, master_nobs);

    /* Only a true identity output can reuse the master already in Stata:
     * using-only keys can insert rows anywhere in the output. */
    int skip_master_store = g_plan.identity && n_using_only == 0 && output_nobs == master_nobs;
    /* update and replace rewrite matched master values of shared variables. */
    int shared_rewritten = 0;
    for (int kv = 0; kv < n_keepusing && (update_mode || replace_mode); kv++) {
        shared_rewritten |= shared_flags[kv] != 0;
    }
    int writes_master = n_stage_vars > 0 && (!skip_master_store || shared_rewritten);
    int need_copy = !g_plan.nocopy && n_stage_vars > 0 && (writes_master || g_plan.hold_undo);
    int write_rows = n_stage_vars > 0 && !skip_master_store;

    double t_load = g_plan.t_load;
    double t_permute = 0.0;
    double t_store = 0.0;
    /* The master data the output is written from: whole master observations in
     * their original order, or full columns when they cannot be staged. Kept
     * as the undo copy unless nocopy. Stata has not been modified by this call
     * yet, so a failed load leaves it intact. */
    ctools_rowstage rows;
    ctools_filtered_data cols;
    ctools_rowstage_init(&rows);
    ctools_filtered_data_init(&cols);
    int have = CMERGE_UNDO_NONE;
    if (need_copy || write_rows) {
        t_start = ctools_timer_ms();
        stata_retcode load_rc = stage_ok
            ? ctools_rowstage_load(&rows, stage_vars, n_stage_vars, NULL, 1, master_nobs)
            : cmerge_load_columns(&cols, stage_vars, n_stage_vars, master_nobs);
        t_load += ctools_timer_ms() - t_start;
        if (load_rc != STATA_OK) {
            free(stage_vars);
            cmerge_release_all();
            SF_error("cmerge: Failed to load master data\n");
            return ctools_stata_rc(load_rc);
        }
        have = stage_ok ? CMERGE_UNDO_ROWS : CMERGE_UNDO_COLUMNS;
    }
    if (need_copy) {
        g_undo.held = have;
        g_undo.master_nobs = master_nobs;
        g_undo.rows = rows;
        g_undo.cols = cols;
        g_undo.col_vars = have == CMERGE_UNDO_COLUMNS ? stage_vars : NULL;
        g_undo.n_cols = have == CMERGE_UNDO_COLUMNS ? n_stage_vars : 0;
        if (g_undo.col_vars) stage_vars = NULL;   /* now owned by the undo copy */
        ctools_rowstage_init(&rows);
        ctools_filtered_data_init(&cols);
    }
    const ctools_rowstage *src_rows = need_copy ? &g_undo.rows : &rows;
    const ctools_filtered_data *src_cols = need_copy ? &g_undo.cols : &cols;
    int *col_vars = need_copy && g_undo.col_vars ? g_undo.col_vars : stage_vars;

    /* One row-tiled pass writes the master rows, the keys of using-only rows,
     * the keepusing variables and _merge. A master that could not be staged
     * (column copy) and the fault_store validation hook keep the writers
     * that move one variable at a time. */
    size_t n_updated = 0, n_conflict = 0;
    int values_changed = 0;   /* update/replace rewrote some master values */
    double t_write_meta = 0.0;
    const int fused = !fault_store && (have == CMERGE_UNDO_ROWS || !write_rows);
    if (fused) {
        t_start = ctools_timer_ms();
        cmerge_writer_t wr;
        memset(&wr, 0, sizeof(wr));
        wr.rows = write_rows ? src_rows : NULL;
        wr.orig_rows = master_orig_rows;
        wr.specs = output_specs;
        wr.output_nobs = output_nobs;
        wr.nkeys = nkeys;
        wr.key_idx = master_key_indices;
        wr.n_keepusing = n_keepusing;
        wr.ku_dest = keepusing_placeholder_indices;
        wr.ku_shared = shared_flags;
        wr.update_mode = update_mode;
        wr.replace_mode = replace_mode;
        wr.merge_var_idx = merge_var_idx;
        stata_retcode write_rc = cmerge_write_output(&wr, &n_updated, &n_conflict, &values_changed);
        t_store = ctools_timer_ms() - t_start;
        /* Without an undo copy, the source of the output goes now (on the
         * release helper once the merge is written). */
        if (write_rc == STATA_OK) {
            ctools_rowstage_release(&rows);
            ctools_filtered_data_release(&cols);
        } else {
            ctools_rowstage_free(&rows);
            ctools_filtered_data_free(&cols);
        }
        free(stage_vars);
        if (write_rc != STATA_OK) {
            SF_error("cmerge: Failed to store output data\n");
            return cmerge_execute_failed(write_rc);
        }
    } else {
        if (write_rows) {
            stata_retcode store_rc;
            t_start = ctools_timer_ms();
            if (have == CMERGE_UNDO_ROWS) {
                /* Write each output row once from the staged rows. master_orig_rows
                 * holds zero-based master rows, with -1 (CTOOLS_ROW_NONE as
                 * perm_idx_t) for using-only rows. */
                size_t n_store = fault_store ? output_nobs / 2 : output_nobs;
                store_rc = ctools_rowstage_store(src_rows, (const perm_idx_t *)master_orig_rows,
                                                 n_store, 1);

                /* Using-only rows take their key values from the using data. */
                for (size_t i = 0; i < n_store && store_rc == STATA_OK && n_using_only > 0; i++) {
                    if (master_orig_rows[i] >= 0) continue;
                    int32_t using_row = output_specs[i].using_sorted_row;
                    if (using_row < 0 || (size_t)using_row >= g_using_cache.nobs) continue;
                    for (int k = 0; k < nkeys && store_rc == STATA_OK; k++) {
                        const stata_variable *key = &g_using_cache.keys.data.vars[k];
                        ST_int key_idx = (ST_int)master_key_indices[k];
                        int io_rc;
                        if (key->type == STATA_TYPE_DOUBLE) {
                            io_rc = (_stata_)->safestore(key_idx, (ST_int)(i + 1),
                                                         key->data.dbl[using_row]);
                        } else {
                            char *value = key->data.str[using_row];
                            io_rc = SF_sstore(key_idx, (ST_int)(i + 1), value ? value : "");
                        }
                        store_rc = ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                    }
                }
            } else {
                store_rc = cmerge_store_columns(src_cols, col_vars, n_stage_vars, &t_permute);
            }
            if (store_rc == STATA_OK && fault_store) store_rc = STATA_ERR_MEMORY;
            t_store = ctools_timer_ms() - t_start;

            if (store_rc != STATA_OK) {
                free(stage_vars);
                ctools_rowstage_free(&rows);
                ctools_filtered_data_free(&cols);
                SF_error("cmerge: Failed to store output data\n");
                return cmerge_execute_failed(store_rc);
            }
        }
        /* Without an undo copy, the source of the output goes now. */
        ctools_rowstage_free(&rows);
        ctools_filtered_data_free(&cols);
        free(stage_vars);

        t_start = ctools_timer_ms();

        /* update: per-row "missing updated" / "nonmissing conflict" flags, filled
         * by the keepusing writers and folded into _merge codes 4 and 5 below. */
        uint8_t *row_flags = NULL;
        if (update_mode && n_keepusing > 0 && output_nobs > 0) {
            row_flags = calloc(output_nobs, 1);
            if (!row_flags) {
                SF_error("cmerge: memory allocation failed\n");
                return cmerge_execute_failed(STATA_ERR_MEMORY);
            }
        }


        /* ===================================================================
         * Step 8: Write keepusing variables from cache (parallel)
         * =================================================================== */

        if (n_keepusing > 0) {
            /* One task per variable and row range: every task owns distinct cells,
             * and a few variables still occupy every worker. */
            size_t ku_chunks = 1;
            if (output_nobs >= (size_t)MIN_OBS_PER_THREAD * 2) {
                int workers = ctools_get_max_threads();
                ku_chunks = output_nobs / MIN_OBS_PER_THREAD;
                if (workers > 0 && ku_chunks > (size_t)workers) ku_chunks = (size_t)workers;
            }
            size_t ku_tasks = (size_t)n_keepusing * ku_chunks;
            cmerge_keepusing_write_args_t *ku_args_array = malloc(ku_tasks * sizeof(cmerge_keepusing_write_args_t));
            if (!ku_args_array) {
                /* Fallback to sequential on allocation failure */
                for (int kv = 0; kv < n_keepusing; kv++) {
                    cmerge_keepusing_write_args_t ku_args;
                    ku_args.keepusing_idx = kv;
                    ku_args.dest_idx = (ST_int)keepusing_placeholder_indices[kv];
                    ku_args.is_shared = shared_flags[kv];
                    ku_args.update_mode = update_mode;
                    ku_args.replace_mode = replace_mode;
                    ku_args.specs = output_specs;
                    ku_args.row_start = 0;
                    ku_args.output_nobs = output_nobs;
                    ku_args.row_flags = row_flags;
                    ku_args.success = 0;
                    cmerge_write_keepusing_var_thread(&ku_args);
                    if (!ku_args.success) {
                        free(row_flags);
                        return cmerge_execute_failed(ctools_spi_status(ku_args.error, STATA_ERR_STATA_WRITE));
                    }
                }
            } else {
                /* Prepare all args */
                for (size_t task = 0; task < ku_tasks; task++) {
                    int kv = (int)(task / ku_chunks);
                    size_t chunk = task % ku_chunks;
                    size_t share = output_nobs / ku_chunks, extra = output_nobs % ku_chunks;
                    ku_args_array[task].keepusing_idx = kv;
                    ku_args_array[task].dest_idx = (ST_int)keepusing_placeholder_indices[kv];
                    ku_args_array[task].is_shared = shared_flags[kv];
                    ku_args_array[task].update_mode = update_mode;
                    ku_args_array[task].replace_mode = replace_mode;
                    ku_args_array[task].specs = output_specs;
                    ku_args_array[task].row_start = chunk * share + (chunk < extra ? chunk : extra);
                    ku_args_array[task].output_nobs = ku_args_array[task].row_start + share + (chunk < extra);
                    ku_args_array[task].row_flags = row_flags;
                    ku_args_array[task].success = 0;
                }

                /* Submit all keepusing writes to thread pool */
                ctools_persistent_pool *pool = ctools_get_global_pool();
                int pool_ok = 0;
                if (pool != NULL && ku_tasks >= 1) {
                    if (ctools_persistent_pool_submit_batch(pool, cmerge_write_keepusing_var_thread,
                                                             ku_args_array, ku_tasks,
                                                             sizeof(cmerge_keepusing_write_args_t)) == 0) {
                        ctools_persistent_pool_wait(pool);
                        pool_ok = 1;
                    }
                }
                if (!pool_ok) {
                    /* Sequential fallback if no pool or submission failed */
                    for (size_t task = 0; task < ku_tasks; task++) {
                        cmerge_write_keepusing_var_thread(&ku_args_array[task]);
                    }
                }
                stata_retcode keepusing_rc = STATA_OK;
                for (size_t task = 0; task < ku_tasks; task++) {
                    if (!ku_args_array[task].success && keepusing_rc == STATA_OK)
                        keepusing_rc = ctools_spi_status(ku_args_array[task].error, STATA_ERR_STATA_WRITE);
                }
                free(ku_args_array);
                if (keepusing_rc != STATA_OK) {
                    free(row_flags);
                    return cmerge_execute_failed(keepusing_rc);
                }
            }
        }

        /* ===================================================================
         * Step 9: Write _merge variable (after the update flags are known)
         * Note: SF_vstore is the only Stata API for writing values.
         * Workers write disjoint, contiguous row ranges of the one variable.
         * With update, matched rows become 4 (missing updated) or 5 (nonmissing
         * conflict, which takes precedence), as with merge.
         * =================================================================== */

        if (row_flags) {
            for (size_t i = 0; i < output_nobs; i++) {
                if (output_specs[i].merge_result != MERGE_RESULT_BOTH) continue;
                if (row_flags[i] & CMERGE_ROW_CONFLICT) n_conflict++;
                else if (row_flags[i] & CMERGE_ROW_UPDATED) n_updated++;
                if ((row_flags[i] & CMERGE_ROW_UPDATED) ||
                    (replace_mode && (row_flags[i] & CMERGE_ROW_CONFLICT)))
                    values_changed = 1;
            }
        }
        if (merge_var_idx > 0) {
            ST_int mvar = (ST_int)merge_var_idx;
            int marker_error = 0;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if(output_nobs >= MIN_OBS_PER_THREAD * 2) reduction(max:marker_error)
            #endif
            for (size_t i = 0; i < output_nobs; i++) {
                double code = (double)output_specs[i].merge_result;
                if (row_flags && output_specs[i].merge_result == MERGE_RESULT_BOTH) {
                    if (row_flags[i] & CMERGE_ROW_CONFLICT) code = 5.0;
                    else if (row_flags[i] & CMERGE_ROW_UPDATED) code = 4.0;
                }
                int error = (_stata_)->safestore(mvar, (ST_int)(i + 1), code);
                if (error > marker_error) marker_error = error;
            }
            if (marker_error) {
                free(row_flags);
                return cmerge_execute_failed(ctools_spi_status(marker_error, STATA_ERR_STATA_WRITE));
            }
        }
        free(row_flags);

        t_write_meta = ctools_timer_ms() - t_start;
    }

    /* ===================================================================
     * Cleanup: the undo copy stays until commit or rollback
     * =================================================================== */

    double t_cleanup_start = ctools_timer_ms();
    const size_t n_master_only = g_plan.n_master_only;
    const size_t n_matched = g_plan.n_matched;
    const int is_identity_permutation = g_plan.identity;
    const double t_sort = g_plan.t_sort, t_merge = g_plan.t_join;
    const double t_reorder = g_plan.t_reorder, t_plan = g_plan.t_total;
    ctools_release_later(g_plan.specs);
    ctools_release_later(g_plan.orig_rows);
    g_plan.specs = NULL;
    g_plan.orig_rows = NULL;
    cmerge_plan_free();

    /* Free using cache (protected by lock); its columns go on the release
     * helper, the strL spans of the .dta path at once */
    cache_lock_acquire();
    ctools_filtered_data_release(&g_using_cache.keys);
    if (!g_using_cache.ku_spans && !g_using_cache.ku_arena)
        ctools_filtered_data_release(&g_using_cache.keepusing);
    ctools_release_later(g_using_cache.perm);
    g_using_cache.perm = NULL;
    cmerge_free_keepusing();
    g_using_cache.loaded = 0;
    g_cache_in_use = 0;
    cache_lock_release();

    double t_cleanup = ctools_timer_ms() - t_cleanup_start;

    /* Calculate total time (plan and execute) */
    double t_total = t_plan + ctools_timer_ms() - t_total_start;

    /* Save merge results to Stata scalars */
    SF_scal_save("_cmerge_N", (double)output_nobs);
    SF_scal_save("_cmerge_N_1", (double)n_master_only);
    SF_scal_save("_cmerge_N_2", (double)n_using_only);
    SF_scal_save("_cmerge_N_3", (double)(n_matched - n_updated - n_conflict));
    SF_scal_save("_cmerge_N_4", (double)n_updated);
    SF_scal_save("_cmerge_N_5", (double)n_conflict);
    /* 1 when the output rows are the master rows in their original order
     * and no master value was rewritten (a sort variable may have changed) */
    SF_scal_save("_cmerge_same_order",
                 (double)(is_identity_permutation && n_using_only == 0 &&
                          output_nobs == master_nobs && !values_changed));

    /* Save Phase 2 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_p2_load_master", t_load / 1000.0);
    SF_scal_save("_cmerge_p2_sort_master", t_sort / 1000.0);
    SF_scal_save("_cmerge_p2_merge_join", t_merge / 1000.0);
    SF_scal_save("_cmerge_p2_reorder", t_reorder / 1000.0);
    SF_scal_save("_cmerge_p2_permute", t_permute / 1000.0);
    SF_scal_save("_cmerge_p2_stream", 0.0);
    SF_scal_save("_cmerge_p2_store", t_store / 1000.0);
    SF_scal_save("_cmerge_p2_write_meta", t_write_meta / 1000.0);
    SF_scal_save("_cmerge_p2_cleanup", t_cleanup / 1000.0);
    SF_scal_save("_cmerge_p2_total", t_total / 1000.0);
    SF_scal_save("_cmerge_n_output_vars", (double)master_nvars);
    CTOOLS_SAVE_THREAD_INFO("_cmerge");

    return STATA_OK;
}

/* Write the original observations back from the undo copy, then release
 * everything the merge held. The varlist must list the original master
 * variables first, in their original order and with the storage types they
 * had during execute; the dataset must have at least the original number of
 * observations. Without an undo copy only the release happens. */
static ST_retcode cmerge_rollback(void)
{
    stata_retcode rc = g_undo.held ? cmerge_undo_restore() : STATA_OK;
    cmerge_release_all();
    if (rc != STATA_OK) {
        SF_error("cmerge: could not restore the master observations\n");
        return ctools_stata_rc(rc);
    }
    return 0;
}

/* ============================================================================
 * Clear cached using data
 * ============================================================================ */

static ST_retcode cmerge_clear(void)
{
    cache_lock_acquire();

    /* Check if cache is currently in use by execute */
    if (g_cache_in_use) {
        cache_lock_release();
        SF_error("cmerge: Cannot clear cache while merge is in progress\n");
        return 459;
    }

    if (g_using_cache.loaded) {
        ctools_filtered_data_free(&g_using_cache.keys);
        cmerge_free_keepusing();
        g_using_cache.loaded = 0;
    }

    cache_lock_release();
    cmerge_plan_free();
    cmerge_undo_free();
    return 0;
}

/* ============================================================================
 * Plugin entry point
 * ============================================================================ */

ST_retcode cmerge_main(const char *args)
{

    if (args == NULL || strlen(args) == 0) {
        SF_error("cmerge: no arguments specified\n");
        return 198;
    }

    /* Use static buffer to avoid potential stack issues */
    static char args_copy[8192];
    if (strlen(args) >= sizeof(args_copy)) return 198;
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *phase = strtok(args_copy, " ");
    const char *rest = args + strlen(phase);
    while (*rest == ' ') rest++;

    if (strcmp(phase, "load_using") == 0) {
        return cmerge_load_using(rest);
    }
    else if (strcmp(phase, "plan") == 0) {
        return cmerge_plan(rest);
    }
    else if (strcmp(phase, "execute") == 0) {
        return cmerge_execute(rest);
    }
    else if (strcmp(phase, "rollback") == 0) {
        return cmerge_rollback();
    }
    else if (strcmp(phase, "commit") == 0) {
        /* The merge stands: the undo copy goes on the release helper */
        ctools_rowstage_release(&g_undo.rows);
        ctools_filtered_data_release(&g_undo.cols);
        cmerge_release_all();
        return 0;
    }
    else if (strcmp(phase, "clear") == 0) {
        return cmerge_clear();
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "cmerge: unknown phase '%s'\n", phase);
        SF_error(msg);
        return 198;
    }
}

/* ============================================================================
 * Cleanup function for ctools_cleanup system
 * ============================================================================ */

void cmerge_cleanup_cache(void)
{
    /* The dispatcher calls this only after workers have joined. A previous
     * interrupted execute may have left the in-use marker set. */
    cache_lock_acquire();
    g_cache_in_use = 0;
    ctools_filtered_data_free(&g_using_cache.keys);
    cmerge_free_keepusing();
    memset(&g_using_cache, 0, sizeof(g_using_cache));
    cache_lock_release();
    cmerge_plan_free();
    cmerge_undo_free();
    free(g_cm_ku_indices);
    free(g_cm_placeholders);
    free(g_cm_shared);
    g_cm_ku_indices = g_cm_placeholders = g_cm_shared = NULL;
    g_cm_ku_cap = g_cm_placeholders_cap = g_cm_shared_cap = 0;
}
