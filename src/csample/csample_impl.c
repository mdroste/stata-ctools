/*
    csample_impl.c
    High-performance random sampling without replacement for Stata

    Performance approach:
    - Thread-local PRNG (xoshiro256**) for parallel random number generation
    - Fisher-Yates partial shuffle for O(k) sampling of k from n
    - Parallel group processing for by() operations
    - Single-pass marking of observations to keep/drop

    Author: ctools project
    License: MIT
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "ctools_error.h"
#include "csample_impl.h"

#define CSAMPLE_MODULE "csample"

/* ===========================================================================
   High-quality PRNG: xoshiro256** (Blackman & Vigna)
   Much better statistical properties than rand()
   =========================================================================== */

typedef struct {
    uint64_t s[4];
} xoshiro256_state;

static inline uint64_t rotl64(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t xoshiro256_next(xoshiro256_state *state) {
    const uint64_t result = rotl64(state->s[1] * 5, 7) * 9;
    const uint64_t t = state->s[1] << 17;

    state->s[2] ^= state->s[0];
    state->s[3] ^= state->s[1];
    state->s[1] ^= state->s[2];
    state->s[0] ^= state->s[3];

    state->s[2] ^= t;
    state->s[3] = rotl64(state->s[3], 45);

    return result;
}

/* SplitMix64 for seeding */
static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void xoshiro256_seed(xoshiro256_state *state, uint64_t seed) {
    uint64_t sm_state = seed;
    state->s[0] = splitmix64(&sm_state);
    state->s[1] = splitmix64(&sm_state);
    state->s[2] = splitmix64(&sm_state);
    state->s[3] = splitmix64(&sm_state);
}

/* Generate uniform random in [0, n) */
static inline size_t xoshiro256_uniform(xoshiro256_state *state, size_t n) {
    if (n == 0) return 0;
    /* Rejection sampling to avoid modulo bias */
    uint64_t threshold = (UINT64_MAX - n + 1) % n;
    uint64_t r;
    do {
        r = xoshiro256_next(state);
    } while (r < threshold);
    return (size_t)(r % n);
}

/* ===========================================================================
   Group Detection (same pattern as cwinsor)
   =========================================================================== */

typedef struct {
    size_t start;
    size_t count;
} group_info;

typedef struct {
    double **by_data;
    size_t nby;
} sort_context;

static sort_context g_sort_ctx;

static int compare_by_indices(const void *a, const void *b) {
    size_t ia = *(const size_t *)a;
    size_t ib = *(const size_t *)b;
    const double miss = SV_missval;

    for (size_t k = 0; k < g_sort_ctx.nby; k++) {
        double va = g_sort_ctx.by_data[k][ia];
        double vb = g_sort_ctx.by_data[k][ib];

        int ma = (va >= miss);
        int mb = (vb >= miss);

        if (ma && mb) continue;
        if (ma) return 1;
        if (mb) return -1;

        if (va < vb) return -1;
        if (va > vb) return 1;
    }
    return (ia < ib) ? -1 : (ia > ib) ? 1 : 0;
}

static void argsort_by_vars(double **by_data, size_t nobs, size_t nby, size_t *indices) {
    for (size_t i = 0; i < nobs; i++) {
        indices[i] = i;
    }
    g_sort_ctx.by_data = by_data;
    g_sort_ctx.nby = nby;
    qsort(indices, nobs, sizeof(size_t), compare_by_indices);
}

static void detect_groups(double **by_data, size_t nobs, size_t nby,
                          group_info *groups, size_t *ngroups) {
    if (nobs == 0) {
        *ngroups = 0;
        return;
    }

    const double miss = SV_missval;
    size_t g = 0;
    groups[0].start = 0;
    groups[0].count = 1;

    for (size_t i = 1; i < nobs; i++) {
        int same_group = 1;

        for (size_t b = 0; b < nby && same_group; b++) {
            double prev = by_data[b][i - 1];
            double curr = by_data[b][i];

            int prev_miss = (prev >= miss);
            int curr_miss = (curr >= miss);

            if (prev_miss && curr_miss) continue;
            if (prev_miss || curr_miss || prev != curr) {
                same_group = 0;
            }
        }

        if (same_group) {
            groups[g].count++;
        } else {
            g++;
            groups[g].start = i;
            groups[g].count = 1;
        }
    }

    *ngroups = g + 1;
}

/* ===========================================================================
   Argument Parsing
   =========================================================================== */

typedef struct {
    double percent;     /* Percentage to keep (0-100), or -1 if using count */
    size_t count;       /* Number to keep, or 0 if using percent */
    int verbose;
    uint64_t seed;      /* Random seed (0 = use time-based) */
} csample_config;

static int parse_option(const char *args, const char *name) {
    if (!args || !name) return 0;
    const char *p = strstr(args, name);
    if (!p) return 0;
    if (p > args && *(p-1) != ' ' && *(p-1) != '\t') return 0;
    char after = *(p + strlen(name));
    if (after != '\0' && after != ' ' && after != '\t' && after != '=') return 0;
    return 1;
}

static double parse_double_option(const char *args, const char *key, double def) {
    if (!args || !key) return def;
    char pat[64];
    snprintf(pat, sizeof(pat), "%s=", key);
    const char *p = strstr(args, pat);
    if (!p) return def;
    char *endptr;
    double result = strtod(p + strlen(pat), &endptr);
    if (endptr == p + strlen(pat)) return def;  /* No conversion */
    return result;
}

static size_t parse_size_option(const char *args, const char *key, size_t def) {
    if (!args || !key) return def;
    char pat[64];
    snprintf(pat, sizeof(pat), "%s=", key);
    const char *p = strstr(args, pat);
    if (!p) return def;
    char *endptr;
    unsigned long long val = strtoull(p + strlen(pat), &endptr, 10);
    if (endptr == p + strlen(pat)) return def;  /* No conversion */
    return (size_t)val;
}

static int parse_int_array(int *arr, size_t count, const char **start) {
    const char *p = *start;
    for (size_t i = 0; i < count; i++) {
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') return -1;
        char *end;
        arr[i] = (int)strtol(p, &end, 10);
        if (end == p) return -1;
        p = end;
    }
    *start = p;
    return 0;
}

/* ===========================================================================
   Main Entry Point
   =========================================================================== */

ST_retcode csample_main(const char *args) {
    double t_start, t_load, t_sort, t_sample, t_store;
    int *by_indices = NULL;
    double **by_data = NULL;
    uint8_t *keep_flags = NULL;
    group_info *groups = NULL;
    size_t *sort_indices = NULL;
    size_t *work_indices = NULL;
    int num_threads = 1;
    xoshiro256_state *thread_rngs = NULL;

    t_start = ctools_timer_seconds();
    t_sort = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CSAMPLE_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    csample_config config;
    config.percent = parse_double_option(args, "percent", -1.0);
    config.count = parse_size_option(args, "count", 0);
    config.verbose = parse_option(args, "verbose");
    config.seed = (uint64_t)parse_size_option(args, "seed", 0);

    /* Validate - must have either percent or count */
    if (config.percent < 0 && config.count == 0) {
        ctools_error(CSAMPLE_MODULE, "must specify percent= or count=");
        return 198;
    }
    if (config.percent >= 0 && (config.percent < 0 || config.percent > 100)) {
        ctools_error(CSAMPLE_MODULE, "percent must be between 0 and 100");
        return 198;
    }

    /* Parse: keep_idx nby by_indices... */
    const char *p = args;
    while (*p == ' ' || *p == '\t') p++;

    char *end;
    int keep_idx = (int)strtol(p, &end, 10);
    if (end == p || keep_idx < 1) {
        ctools_error(CSAMPLE_MODULE, "invalid keep_idx");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    size_t nby = (size_t)strtol(p, &end, 10);
    if (end == p) {
        ctools_error(CSAMPLE_MODULE, "invalid nby");
        return 198;
    }
    p = end;

    if (nby > 0) {
        by_indices = (int *)malloc(nby * sizeof(int));
        if (!by_indices) {
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }
        if (parse_int_array(by_indices, nby, &p) != 0) {
            free(by_indices);
            ctools_error(CSAMPLE_MODULE, "failed to parse by indices");
            return 198;
        }
    }

    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    size_t nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        if (by_indices) free(by_indices);
        ctools_error(CSAMPLE_MODULE, "no observations");
        return 2000;
    }

    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif

    /* Use time-based seed if not specified */
    if (config.seed == 0) {
        config.seed = (uint64_t)time(NULL) ^ ((uint64_t)clock() << 32);
    }

    /* === Load Phase === */
    double load_start = ctools_timer_seconds();

    /* Allocate keep flags */
    keep_flags = (uint8_t *)calloc(nobs, sizeof(uint8_t));
    if (!keep_flags) {
        if (by_indices) free(by_indices);
        ctools_error_alloc(CSAMPLE_MODULE);
        return 920;
    }

    /* Load by-variables if specified */
    if (nby > 0) {
        by_data = (double **)malloc(nby * sizeof(double *));
        if (!by_data) {
            free(keep_flags);
            free(by_indices);
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }

        for (size_t b = 0; b < nby; b++) {
            by_data[b] = (double *)malloc(nobs * sizeof(double));
            if (!by_data[b]) {
                for (size_t j = 0; j < b; j++) free(by_data[j]);
                free(by_data);
                free(keep_flags);
                free(by_indices);
                ctools_error_alloc(CSAMPLE_MODULE);
                return 920;
            }

            int idx = by_indices[b];
            for (size_t i = 0; i < nobs; i++) {
                double val;
                if (SF_vdata(idx, (ST_int)(obs1 + i), &val) != 0) val = SV_missval;
                by_data[b][i] = val;
            }
        }
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort Phase (if by-variables) === */
    double sort_start = ctools_timer_seconds();

    size_t ngroups = 1;
    groups = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    if (!groups) {
        if (by_data) {
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
        }
        free(keep_flags);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CSAMPLE_MODULE);
        return 920;
    }

    if (nby > 0) {
        sort_indices = (size_t *)malloc(nobs * sizeof(size_t));
        if (!sort_indices) {
            free(groups);
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
            free(keep_flags);
            free(by_indices);
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }

        argsort_by_vars(by_data, nobs, nby, sort_indices);

        /* Reorder by_data arrays by sort order */
        double *temp = (double *)malloc(nobs * sizeof(double));
        if (!temp) {
            free(sort_indices);
            free(groups);
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
            free(keep_flags);
            free(by_indices);
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }

        for (size_t b = 0; b < nby; b++) {
            for (size_t i = 0; i < nobs; i++) {
                temp[i] = by_data[b][sort_indices[i]];
            }
            memcpy(by_data[b], temp, nobs * sizeof(double));
        }
        free(temp);

        detect_groups(by_data, nobs, nby, groups, &ngroups);
    } else {
        groups[0].start = 0;
        groups[0].count = nobs;
    }

    t_sort = ctools_timer_seconds() - sort_start;

    /* === Sample Phase === */
    double sample_start = ctools_timer_seconds();

    /* Allocate thread-local RNGs */
    thread_rngs = (xoshiro256_state *)malloc(num_threads * sizeof(xoshiro256_state));
    if (!thread_rngs) {
        if (sort_indices) free(sort_indices);
        free(groups);
        if (by_data) {
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
        }
        free(keep_flags);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CSAMPLE_MODULE);
        return 920;
    }

    /* Seed each thread's RNG uniquely */
    for (int t = 0; t < num_threads; t++) {
        xoshiro256_seed(&thread_rngs[t], config.seed + (uint64_t)t * 0x9e3779b97f4a7c15ULL);
    }

    /* Allocate work indices (max group size needed) */
    size_t max_group_size = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (groups[g].count > max_group_size)
            max_group_size = groups[g].count;
    }

    work_indices = (size_t *)malloc(num_threads * max_group_size * sizeof(size_t));
    if (!work_indices) {
        free(thread_rngs);
        if (sort_indices) free(sort_indices);
        free(groups);
        if (by_data) {
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
        }
        free(keep_flags);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CSAMPLE_MODULE);
        return 920;
    }

    /* Process each group */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for (size_t g = 0; g < ngroups; g++) {
        #ifdef _OPENMP
        int tid = omp_get_thread_num();
        #else
        int tid = 0;
        #endif

        size_t start = groups[g].start;
        size_t count = groups[g].count;

        if (count == 0) continue;

        /* Calculate how many to keep in this group */
        size_t keep_count;
        if (config.count > 0) {
            /* Fixed count per group */
            keep_count = (config.count < count) ? config.count : count;
        } else {
            /* Percentage */
            keep_count = (size_t)round((config.percent / 100.0) * (double)count);
            if (keep_count > count) keep_count = count;
        }

        if (keep_count == count) {
            /* Keep all - mark all as selected */
            for (size_t i = 0; i < count; i++) {
                size_t obs_idx = (nby > 0) ? sort_indices[start + i] : (start + i);
                keep_flags[obs_idx] = 1;
            }
        } else if (keep_count == 0) {
            /* Keep none - flags already 0 */
        } else {
            /* Use Fisher-Yates to select keep_count from count */
            size_t *indices = &work_indices[tid * max_group_size];

            /* Initialize indices for this group */
            for (size_t i = 0; i < count; i++) {
                indices[i] = i;
            }

            /* Fisher-Yates partial shuffle */
            xoshiro256_state *rng = &thread_rngs[tid];
            for (size_t i = 0; i < keep_count; i++) {
                size_t j = i + xoshiro256_uniform(rng, count - i);
                size_t tmp = indices[i];
                indices[i] = indices[j];
                indices[j] = tmp;
            }

            /* Mark selected observations */
            for (size_t i = 0; i < keep_count; i++) {
                size_t local_idx = indices[i];
                size_t obs_idx = (nby > 0) ? sort_indices[start + local_idx] : (start + local_idx);
                keep_flags[obs_idx] = 1;
            }
        }
    }

    t_sample = ctools_timer_seconds() - sample_start;

    /* === Store Phase === */
    double store_start = ctools_timer_seconds();

    /* Write keep flags to Stata */
    size_t obs_idx = 0;  /* Index into keep_flags array */
    for (ST_int obs = obs1; obs <= obs2; obs++) {
        /* Only process observations that satisfy the if condition */
        if (!SF_ifobs(obs)) continue;

        double val = keep_flags[obs_idx] ? 1.0 : 0.0;
        SF_vstore(keep_idx, obs, val);
        obs_idx++;
    }

    t_store = ctools_timer_seconds() - store_start;

    /* Count kept/dropped */
    size_t kept = 0;
    for (size_t i = 0; i < nobs; i++) {
        if (keep_flags[i]) kept++;
    }
    size_t dropped = nobs - kept;

    /* === Cleanup === */
    free(work_indices);
    free(thread_rngs);
    if (sort_indices) free(sort_indices);
    free(groups);
    if (by_data) {
        for (size_t b = 0; b < nby; b++) free(by_data[b]);
        free(by_data);
    }
    free(keep_flags);
    if (by_indices) free(by_indices);

    /* Store timing scalars */
    double total = ctools_timer_seconds() - t_start;
    SF_scal_save("_csample_time_load", t_load);
    SF_scal_save("_csample_time_sort", t_sort);
    SF_scal_save("_csample_time_sample", t_sample);
    SF_scal_save("_csample_time_store", t_store);
    SF_scal_save("_csample_time_total", total);
    SF_scal_save("_csample_nobs", (double)nobs);
    SF_scal_save("_csample_kept", (double)kept);
    SF_scal_save("_csample_dropped", (double)dropped);
    SF_scal_save("_csample_ngroups", (double)ngroups);

    #ifdef _OPENMP
    SF_scal_save("_csample_threads", (double)num_threads);
    #else
    SF_scal_save("_csample_threads", 1.0);
    #endif

    if (config.verbose) {
        ctools_msg(CSAMPLE_MODULE, "%zu obs, %zu groups, kept %zu, dropped %zu, %d threads",
                   nobs, ngroups, kept, dropped, num_threads);
        ctools_msg(CSAMPLE_MODULE, "load=%.4f sort=%.4f sample=%.4f store=%.4f total=%.4f",
                   t_load, t_sort, t_sample, t_store, total);
    }

    return 0;
}
