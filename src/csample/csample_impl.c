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
#include <errno.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Checked SPI writes resolve positions in the plugin varlist, even after
 * estimation has introduced internal Stata variables. */
#define SD_SAFEMODE
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_parse.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "csample_impl.h"

#define CSAMPLE_MODULE "csample"

#include "ctools_rng.h"
#include "ctools_group.h"

/* ===========================================================================
   Argument Parsing
   =========================================================================== */

typedef struct {
    double percent;     /* Percentage to keep (0-100), or -1 if using count */
    size_t count;       /* Number to keep, or 0 if using percent */
    int verbose;
    uint64_t seed;      /* Random seed (0 = use time-based) */
} csample_config;





/* Parse key=<unsigned decimal integer> as a whole token. Returns 1 when the
 * option is present and valid, 0 when absent, and -1 when malformed: a partial
 * parse (for example "2.30584300921e+18" read as 2) must not be accepted. */


/* The wrapper draws the seed from Stata's RNG as two exact 31-bit halves
 * (a single ~62-bit value loses digits in a Stata macro). Returns 1 when a
 * seed was supplied, 0 when not, -1 when malformed. */




/* ===========================================================================
   Main Entry Point
   =========================================================================== */

ST_retcode csample_main(const char *args) {
    double t_start, t_load, t_sort, t_sample, t_store;
    int *by_indices = NULL;
    int *sort_vars = NULL;
    uint8_t *keep_flags = NULL;
    ctools_group_range *groups = NULL;
    perm_idx_t *sort_perm = NULL;  /* Saved permutation for mapping back */
    size_t *work_indices = NULL;
    ctools_filtered_data by_filtered;
    perm_idx_t *obs_map = NULL;
    perm_idx_t *if_map = NULL;  /* obs_map when built here by SF_ifobs */
    ST_int obs_start = 0;  /* For null obs_map: first observation number */
    stata_retcode rc = STATA_OK;
    int num_threads = 1;
    ctools_rng *thread_rngs = NULL;

    ctools_filtered_data_init(&by_filtered);
    t_start = ctools_timer_seconds();
    t_sort = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CSAMPLE_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    csample_config config;
    config.percent = ctools_parse_double_option(args, "percent", -1.0);
    uint64_t count_val = 0;
    int count_rc = ctools_parse_u64_option(args, "count", &count_val);
    config.count = (size_t)count_val;
    config.verbose = ctools_parse_bool_option(args, "verbose");
    config.seed = 0;
    int seed_rc = ctools_parse_seed_option(args, &config.seed);
    if (count_rc < 0 || seed_rc < 0 || count_val > SIZE_MAX || !isfinite(config.percent)) {
        ctools_error(CSAMPLE_MODULE, "invalid count= or seed option");
        return 198;
    }
    int skip_if = ctools_parse_bool_option(args, "skipif");
    /* touse: the last plugin variable is the sample marker, in place of the
     * if condition; the shared loader reads it in parallel. */
    const int selection = skip_if ? CTOOLS_LOAD_SKIP_IF
                        : ctools_parse_bool_option(args, "touse") ? CTOOLS_LOAD_MARKER
                                                                  : CTOOLS_LOAD_CHECK_IF;

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

    int keep_idx;
    if (ctools_parse_next_int(&p, &keep_idx) || keep_idx < 1) {
        ctools_error(CSAMPLE_MODULE, "invalid keep_idx");
        return 198;
    }

    while (*p == ' ' || *p == '\t') p++;
    int nby_value;
    if (ctools_parse_next_int(&p, &nby_value) || nby_value < 0 ||
        nby_value > SF_nvars()) {
        ctools_error(CSAMPLE_MODULE, "invalid nby");
        return 198;
    }
    size_t nby = (size_t)nby_value;

    if (nby > 0) {
        by_indices = (int *)malloc(nby * sizeof(int));
        if (!by_indices) {
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }
        if (ctools_parse_int_array(by_indices, nby, &p) != 0) {
            free(by_indices);
            ctools_error(CSAMPLE_MODULE, "failed to parse by indices");
            return 198;
        }
    }

    #ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    #endif

    /* Use time-based seed if not specified */
    if (seed_rc == 0) {
        config.seed = (uint64_t)time(NULL) ^ ((uint64_t)clock() << 32);
    }

    /* === Load Phase using ctools_data_load === */
    double load_start = ctools_timer_seconds();

    size_t nobs;
    size_t ngroups = 1;

    /* Load by-variables if specified using filtered loading */
    if (nby > 0) {
        rc = ctools_data_load(&by_filtered, by_indices, nby, 0, 0, selection);
        if (rc != STATA_OK) {
            ctools_filtered_data_free(&by_filtered);
            free(by_indices);
            ctools_error(CSAMPLE_MODULE, "failed to load data");
            return ctools_stata_rc(rc);
        }
        nobs = by_filtered.data.nobs;
        obs_map = by_filtered.obs_map;
    } else {
        /* No by-variables */
        ST_int in1 = SF_in1();
        ST_int in2 = SF_in2();

        if (skip_if) {
            /* No if/in and no by-variables: skip data load entirely */
            nobs = (size_t)(in2 - in1 + 1);
            obs_map = NULL;
            obs_start = in1;
        } else if (selection == CTOOLS_LOAD_MARKER) {
            /* The sample marker, loaded as the only column, selects the rows */
            int marker = SF_nvars();
            rc = ctools_data_load(&by_filtered, &marker, 1, 0, 0,
                                  CTOOLS_LOAD_MARKER | CTOOLS_LOAD_NO_SORT_ORDER);
            if (rc != STATA_OK) {
                ctools_filtered_data_free(&by_filtered);
                ctools_error(CSAMPLE_MODULE, "failed to load data");
                return ctools_stata_rc(rc);
            }
            nobs = by_filtered.data.nobs;
            obs_map = by_filtered.obs_map;
        } else {
            /* Has if/in: must check SF_ifobs */
            nobs = 0;
            for (ST_int obs = in1; obs <= in2; obs++) {
                if (SF_ifobs(obs)) nobs++;
            }

            if (nobs == 0) {
                if (by_indices) free(by_indices);
                ctools_error(CSAMPLE_MODULE, "no observations");
                return 2000;
            }

            obs_map = if_map = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
            if (!obs_map) {
                if (by_indices) free(by_indices);
                ctools_error_alloc(CSAMPLE_MODULE);
                return 920;
            }

            size_t idx = 0;
            for (ST_int obs = in1; obs <= in2; obs++) {
                if (SF_ifobs(obs)) {
                    obs_map[idx++] = (perm_idx_t)obs;
                }
            }
        }
    }

    if (nobs == 0) {
        ctools_filtered_data_free(&by_filtered);
        free(if_map);
        if (by_indices) free(by_indices);
        ctools_error(CSAMPLE_MODULE, "no observations");
        return 2000;
    }

    /* Allocate keep flags */
    keep_flags = (uint8_t *)calloc(nobs, sizeof(uint8_t));
    if (!keep_flags) {
        ctools_filtered_data_free(&by_filtered);
        free(if_map);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CSAMPLE_MODULE);
        return 920;
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort Phase using ctools_sort_dispatch (if by-variables) === */
    double sort_start = ctools_timer_seconds();

    groups = (ctools_group_range *)malloc((nobs + 1) * sizeof(ctools_group_range));
    if (!groups) {
        ctools_filtered_data_free(&by_filtered);
        free(if_map);
        free(keep_flags);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CSAMPLE_MODULE);
        return 920;
    }

    if (nby > 0) {
        /* Build sort_vars array (1-based indices into loaded data) */
        sort_vars = (int *)malloc(nby * sizeof(int));
        if (!sort_vars) {
            free(groups);
            ctools_filtered_data_free(&by_filtered);
            free(keep_flags);
            free(by_indices);
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }
        for (size_t i = 0; i < nby; i++) {
            sort_vars[i] = (int)(i + 1);  /* 1-based */
        }

        /* Sort using counting sort for integer group variables (optimal for ties),
           fall back to LSD radix if data isn't suitable for counting sort */
        rc = ctools_sort_dispatch(&by_filtered.data, sort_vars, nby, SORT_ALG_COUNTING);
        if (rc != STATA_OK) {
            free(sort_vars);
            free(groups);
            ctools_filtered_data_free(&by_filtered);
            free(keep_flags);
            free(by_indices);
            ctools_error(CSAMPLE_MODULE, "sort failed");
            return ctools_stata_rc(rc);
        }

        /* Save the sort permutation before applying (for mapping back to original indices) */
        sort_perm = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
        if (!sort_perm) {
            free(sort_vars);
            free(groups);
            ctools_filtered_data_free(&by_filtered);
            free(keep_flags);
            free(by_indices);
            ctools_error_alloc(CSAMPLE_MODULE);
            return 920;
        }
        memcpy(sort_perm, by_filtered.data.sort_order, nobs * sizeof(perm_idx_t));

        /* Apply permutation to physically reorder the data */
        rc = ctools_apply_permutation(&by_filtered.data);
        if (rc != STATA_OK) {
            free(sort_perm);
            free(sort_vars);
            free(groups);
            ctools_filtered_data_free(&by_filtered);
            free(keep_flags);
            free(by_indices);
            ctools_error(CSAMPLE_MODULE, "permutation failed");
            return ctools_stata_rc(rc);
        }

        free(sort_vars);

        /* Detect groups from sorted data */
        ctools_group_ranges(&by_filtered.data, NULL, NULL, nby, groups, &ngroups);
    } else {
        groups[0].start = 0;
        groups[0].count = nobs;
    }

    t_sort = ctools_timer_seconds() - sort_start;

    /* === Sample Phase === */
    double sample_start = ctools_timer_seconds();

    /* Only active groups need workers; bound aggregate shuffle scratch. */
    size_t max_group_size = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (groups[g].count > max_group_size) max_group_size = groups[g].count;
    }
    num_threads = ctools_workspace_threads(num_threads, ngroups,
                                           max_group_size * sizeof(size_t));
    thread_rngs = (ctools_rng *)ctools_safe_malloc2((size_t)num_threads, sizeof(ctools_rng));
    work_indices = (size_t *)ctools_safe_malloc3((size_t)num_threads, max_group_size, sizeof(size_t));
    if (!thread_rngs || !work_indices) {
        ctools_error_alloc(CSAMPLE_MODULE);
        rc = 920;
        goto cleanup;
    }

    /* Process each group */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
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
                size_t obs_idx = (nby > 0) ? sort_perm[start + i] : (start + i);
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

            /* Seed RNG based on group index (not thread) for reproducibility */
            ctools_rng_seed(&thread_rngs[tid], config.seed + (uint64_t)g * 0x9e3779b97f4a7c15ULL);
            ctools_rng *rng = &thread_rngs[tid];

            /* Fisher-Yates partial shuffle */
            for (size_t i = 0; i < keep_count; i++) {
                size_t j = i + ctools_rng_uniform(rng, count - i);
                size_t tmp = indices[i];
                indices[i] = indices[j];
                indices[j] = tmp;
            }

            /* Mark selected observations */
            for (size_t i = 0; i < keep_count; i++) {
                size_t local_idx = indices[i];
                size_t obs_idx = (nby > 0) ? sort_perm[start + local_idx] : (start + local_idx);
                keep_flags[obs_idx] = 1;
            }
        }
    }

    t_sample = ctools_timer_seconds() - sample_start;

    /* === Store Phase using obs_map === */
    double store_start = ctools_timer_seconds();

    /* Write keep flags to Stata using obs_map (or arithmetic if NULL) */
    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs_map ? (ST_int)obs_map[i] : (obs_start + (ST_int)i);
        double val = keep_flags[i] ? 1.0 : 0.0;
        rc = SF_vstore(keep_idx, obs, val);
        if (rc) goto cleanup;
    }

    t_store = ctools_timer_seconds() - store_start;

    /* Count kept/dropped with parallel reduction */
    size_t kept = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:kept) schedule(static)
    #endif
    for (size_t i = 0; i < nobs; i++) {
        if (keep_flags[i]) kept++;
    }
    size_t dropped = nobs - kept;

cleanup:
    /* === Cleanup === */
    free(work_indices);
    free(thread_rngs);
    if (sort_perm) free(sort_perm);
    free(groups);
    ctools_filtered_data_free(&by_filtered);
    free(if_map);
    free(keep_flags);
    if (by_indices) free(by_indices);

    if (rc) return rc;

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
