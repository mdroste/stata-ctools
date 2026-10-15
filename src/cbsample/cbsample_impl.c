/*
    cbsample_impl.c
    High-performance bootstrap sampling with replacement for Stata

    Performance approach:
    - Thread-local PRNG (xoshiro256**) for parallel random number generation
    - O(n) sampling with replacement
    - Parallel strata processing
    - Cluster-aware bootstrap (resample clusters, include all within)

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

/* Resolve output indices in the plugin varlist, including hidden Stata variables. */
#define SD_SAFEMODE
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_parse.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_spi.h"
#include "cbsample_impl.h"

#define CBSAMPLE_MODULE "cbsample"

#include "ctools_rng.h"
#include "ctools_group.h"

/* ===========================================================================
   Argument Parsing
   =========================================================================== */

typedef struct {
    size_t n;           /* Number of observations in bootstrap sample (0 = same as original) */
    int verbose;
    uint64_t seed;
} cbsample_config;



/* Parse key=<unsigned decimal integer> as a whole token. Returns 1 when the
 * option is present and valid, 0 when absent, and -1 when malformed: a partial
 * parse (for example "2.30584300921e+18" read as 2) must not be accepted. */


/* The wrapper draws the seed from Stata's RNG as two exact 31-bit halves
 * (a single ~62-bit value loses digits in a Stata macro). Returns 1 when a
 * seed was supplied, 0 when not, -1 when malformed. */




/* ===========================================================================
   Main Entry Point
   =========================================================================== */

ST_retcode cbsample_main(const char *args) {
    double t_start, t_load, t_sort_cluster, t_sample, t_store;
    int *cluster_indices = NULL;
    int *strata_indices = NULL;
    int *combined_indices = NULL;
    int *sort_vars = NULL;
    double *freq_weights = NULL;
    ctools_group_range *clusters = NULL;
    ctools_group_range *strata = NULL;
    perm_idx_t *sort_perm = NULL;  /* Saved permutation for mapping back */
    ctools_filtered_data by_filtered;
    perm_idx_t *obs_map = NULL;
    ST_int obs_start = 0;  /* For null obs_map: first observation number */
    stata_retcode rc = STATA_OK;
    int num_threads = 1;
    ctools_rng *thread_rngs = NULL;

    ctools_filtered_data_init(&by_filtered);
    t_start = ctools_timer_seconds();
    t_sort_cluster = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CBSAMPLE_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    cbsample_config config;
    uint64_t n_val = 0;
    int n_rc = ctools_parse_u64_option(args, "n", &n_val);
    config.n = (size_t)n_val;
    config.verbose = ctools_parse_bool_option(args, "verbose");
    config.seed = 0;
    int seed_rc = ctools_parse_seed_option(args, &config.seed);
    if (n_rc < 0 || seed_rc < 0 || n_val > SIZE_MAX) {
        ctools_error(CBSAMPLE_MODULE, "invalid n= or seed option");
        return 198;
    }
    int skip_if = ctools_parse_bool_option(args, "skipif");

    /* Parse: freq_idx ncluster nstrata cluster_indices... strata_indices... */
    const char *p = args;
    while (*p == ' ' || *p == '\t') p++;

    int freq_idx;
    if (ctools_parse_next_int(&p, &freq_idx) || freq_idx < 1) {
        ctools_error(CBSAMPLE_MODULE, "invalid freq_idx");
        return 198;
    }

    while (*p == ' ' || *p == '\t') p++;
    int ncluster_value;
    if (ctools_parse_next_int(&p, &ncluster_value) || ncluster_value < 0 ||
        ncluster_value > SF_nvars()) {
        ctools_error(CBSAMPLE_MODULE, "invalid ncluster");
        return 198;
    }
    size_t ncluster = (size_t)ncluster_value;

    while (*p == ' ' || *p == '\t') p++;
    int nstrata_value;
    if (ctools_parse_next_int(&p, &nstrata_value) || nstrata_value < 0 ||
        nstrata_value > SF_nvars()) {
        ctools_error(CBSAMPLE_MODULE, "invalid nstrata");
        return 198;
    }
    size_t nstrata = (size_t)nstrata_value;

    if (ncluster > 0) {
        cluster_indices = (int *)malloc(ncluster * sizeof(int));
        if (!cluster_indices) {
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }
        if (ctools_parse_int_array(cluster_indices, ncluster, &p) != 0) {
            free(cluster_indices);
            ctools_error(CBSAMPLE_MODULE, "failed to parse cluster indices");
            return 198;
        }
    }

    if (nstrata > 0) {
        strata_indices = (int *)malloc(nstrata * sizeof(int));
        if (!strata_indices) {
            if (cluster_indices) free(cluster_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }
        if (ctools_parse_int_array(strata_indices, nstrata, &p) != 0) {
            if (cluster_indices) free(cluster_indices);
            free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "failed to parse strata indices");
            return 198;
        }
    }

    #ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    #endif

    if (seed_rc == 0) {
        config.seed = (uint64_t)time(NULL) ^ ((uint64_t)clock() << 32);
    }

    /* === Load Phase using ctools_data_load === */
    double load_start = ctools_timer_seconds();

    size_t total_by = nstrata + ncluster;
    size_t nclusters = 1;
    size_t nstrata_groups = 1;
    size_t nobs;

    if (total_by > 0) {
        /* Build combined indices array: strata first, then cluster */
        combined_indices = (int *)malloc(total_by * sizeof(int));
        if (!combined_indices) {
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }

        for (size_t s = 0; s < nstrata; s++) {
            combined_indices[s] = strata_indices[s];
        }
        for (size_t c = 0; c < ncluster; c++) {
            combined_indices[nstrata + c] = cluster_indices[c];
        }

        /* Load data using filtered loading */
        rc = ctools_data_load(&by_filtered, combined_indices, total_by, 0, 0, 0);
        if (rc != STATA_OK) {
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "failed to load data");
            return ctools_stata_rc(rc);
        }
        nobs = by_filtered.data.nobs;
        obs_map = by_filtered.obs_map;
    } else {
        /* No strata/cluster */
        ST_int in1 = SF_in1();
        ST_int in2 = SF_in2();

        if (skip_if) {
            /* No if/in and no by-variables: skip data load entirely */
            nobs = (size_t)(in2 - in1 + 1);
            obs_map = NULL;
            obs_start = in1;
        } else {
            /* Has if/in: must check SF_ifobs */
            nobs = 0;
            for (ST_int obs = in1; obs <= in2; obs++) {
                if (SF_ifobs(obs)) nobs++;
            }

            if (nobs == 0) {
                if (cluster_indices) free(cluster_indices);
                if (strata_indices) free(strata_indices);
                ctools_error(CBSAMPLE_MODULE, "no observations");
                return 2000;
            }

            obs_map = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
            if (!obs_map) {
                if (cluster_indices) free(cluster_indices);
                if (strata_indices) free(strata_indices);
                ctools_error_alloc(CBSAMPLE_MODULE);
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
        if (total_by > 0) ctools_filtered_data_free(&by_filtered);
        else if (obs_map) free(obs_map);
        if (combined_indices) free(combined_indices);
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error(CBSAMPLE_MODULE, "no observations");
        return 2000;
    }

    /* n=0 means draw each stratum's original number of sampling units. */

    freq_weights = (double *)calloc(nobs, sizeof(double));
    if (!freq_weights) {
        if (total_by > 0) ctools_filtered_data_free(&by_filtered);
        else if (obs_map) free(obs_map);
        if (combined_indices) free(combined_indices);
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort and Detect Groups Phase using ctools_sort_dispatch === */
    double sort_start = ctools_timer_seconds();

    /* Allocate groups */
    clusters = (ctools_group_range *)malloc((nobs + 1) * sizeof(ctools_group_range));
    strata = (ctools_group_range *)malloc((nobs + 1) * sizeof(ctools_group_range));
    if (!clusters || !strata) {
        if (clusters) free(clusters);
        if (strata) free(strata);
        if (total_by > 0) {
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
        } else {
            free(obs_map);
        }
        free(freq_weights);
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    if (total_by > 0) {
        /* Build sort_vars array (1-based indices into loaded data) */
        sort_vars = (int *)malloc(total_by * sizeof(int));
        if (!sort_vars) {
            free(clusters);
            free(strata);
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }
        for (size_t i = 0; i < total_by; i++) {
            sort_vars[i] = (int)(i + 1);  /* 1-based */
        }

        /* Sort using counting sort for integer group variables (optimal for ties),
           fall back to LSD radix if data isn't suitable for counting sort */
        rc = ctools_sort_dispatch(&by_filtered.data, sort_vars, total_by, SORT_ALG_COUNTING);
        if (rc != STATA_OK) {
            free(sort_vars);
            free(clusters);
            free(strata);
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "sort failed");
            return ctools_stata_rc(rc);
        }

        /* Save the sort permutation before applying (for mapping back to original indices) */
        sort_perm = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
        if (!sort_perm) {
            free(sort_vars);
            free(clusters);
            free(strata);
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }
        memcpy(sort_perm, by_filtered.data.sort_order, nobs * sizeof(perm_idx_t));

        /* Apply permutation to physically reorder the data */
        rc = ctools_apply_permutation(&by_filtered.data);
        if (rc != STATA_OK) {
            free(sort_perm);
            free(sort_vars);
            free(clusters);
            free(strata);
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "permutation failed");
            return ctools_stata_rc(rc);
        }

        free(sort_vars);
    }

    /* Detect strata groups (using first nstrata variables) */
    if (nstrata > 0) {
        /* Create a temporary view with just strata vars for group detection */
        stata_data strata_view;
        strata_view.nobs = by_filtered.data.nobs;
        strata_view.nvars = nstrata;
        strata_view.vars = by_filtered.data.vars;  /* First nstrata vars are strata */
        strata_view.sort_order = NULL;
        ctools_group_ranges(&strata_view, NULL, NULL, nstrata, strata, &nstrata_groups);
    } else {
        strata[0].start = 0;
        strata[0].count = nobs;
        nstrata_groups = 1;
    }

    /* Detect clusters (using cluster variables at offset nstrata) */
    if (ncluster > 0) {
        /* Create a temporary view with just cluster vars */
        stata_data cluster_view;
        cluster_view.nobs = by_filtered.data.nobs;
        cluster_view.nvars = total_by;
        cluster_view.vars = by_filtered.data.vars;  /* Include stratum boundaries. */
        cluster_view.sort_order = NULL;
        ctools_group_ranges(&cluster_view, NULL, NULL, total_by, clusters, &nclusters);
    } else {
        /* No clustering: each observation is its own "cluster" */
        for (size_t i = 0; i < nobs; i++) {
            clusters[i].start = i;
            clusters[i].count = 1;
        }
        nclusters = nobs;
    }

    t_sort_cluster = ctools_timer_seconds() - sort_start;

    /* === Bootstrap Sample Phase === */
    double sample_start = ctools_timer_seconds();

    /* Allocate thread-local RNGs */
    /* Strata are processed serially for reproducible draws. */
    num_threads = 1;
    thread_rngs = (ctools_rng *)malloc(sizeof(ctools_rng));
    if (!thread_rngs) {
        if (sort_perm) free(sort_perm);
        free(clusters);
        free(strata);
        if (total_by > 0) {
            ctools_filtered_data_free(&by_filtered);
            free(combined_indices);
        } else {
            free(obs_map);
        }
        free(freq_weights);
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    /* Process each stratum - use stratum-based seeding for reproducibility */
    if (ncluster > 0) {
        /* Cluster bootstrap: resample clusters with replacement within each stratum */
        /* For each stratum, identify which clusters are in it and sample from those */

        /* Strata and clusters are runs of the rows sorted by (strata,
         * cluster), both in row order, and cluster runs also break at every
         * stratum boundary, so each stratum owns a contiguous run of
         * clusters: one cursor walks them once across all strata. */
        size_t cl_cursor = 0;
        for (size_t st = 0; st < nstrata_groups; st++) {
            size_t st_start = strata[st].start;
            size_t st_end = st_start + strata[st].count;

            /* Find clusters within this stratum */
            while (cl_cursor < nclusters && clusters[cl_cursor].start < st_start) cl_cursor++;
            size_t first_cluster = cl_cursor;
            while (cl_cursor < nclusters && clusters[cl_cursor].start < st_end) cl_cursor++;
            size_t last_cluster = cl_cursor;

            size_t n_clusters_in_stratum = last_cluster - first_cluster;
            if (n_clusters_in_stratum == 0) continue;

            size_t target_n = config.n ? config.n : n_clusters_in_stratum;
            if (target_n > n_clusters_in_stratum) { rc = 498; goto cleanup; }

            /* Seed RNG based on stratum index for reproducibility */
            ctools_rng_seed(&thread_rngs[0], config.seed + (uint64_t)st * 0x9e3779b97f4a7c15ULL);
            ctools_rng *rng = &thread_rngs[0];

            /* Explicit n counts clusters within each stratum. */
            for (size_t draw = 0; draw < target_n; draw++) {
                size_t cluster_idx = first_cluster + ctools_rng_uniform(rng, n_clusters_in_stratum);
                size_t cl_start = clusters[cluster_idx].start;
                size_t cl_count = clusters[cluster_idx].count;

                /* Increment frequency weight for all obs in this cluster */
                for (size_t j = 0; j < cl_count; j++) {
                    size_t obs_idx = sort_perm ? sort_perm[cl_start + j] : (cl_start + j);
                    freq_weights[obs_idx] += 1.0;
                }
            }
        }
    } else {
        /* No clustering: simple bootstrap of observations */
        /* For reproducibility, use sequential processing with stratum-based seeding */
        for (size_t st = 0; st < nstrata_groups; st++) {
            size_t st_start = strata[st].start;
            size_t st_count = strata[st].count;

            size_t target_n = config.n ? config.n : st_count;
            if (target_n > st_count) { rc = 498; goto cleanup; }

            /* Seed RNG based on stratum index for reproducibility */
            ctools_rng_seed(&thread_rngs[0], config.seed + (uint64_t)st * 0x9e3779b97f4a7c15ULL);
            ctools_rng *rng = &thread_rngs[0];

            for (size_t draw = 0; draw < target_n; draw++) {
                size_t local_idx = ctools_rng_uniform(rng, st_count);
                size_t obs_idx = sort_perm ? sort_perm[st_start + local_idx] : (st_start + local_idx);
                freq_weights[obs_idx] += 1.0;
            }
        }
    }

    t_sample = ctools_timer_seconds() - sample_start;

    /* === Store Phase using obs_map === */
    double store_start = ctools_timer_seconds();

    /* Write frequency weights to Stata using obs_map (or arithmetic if NULL) */
    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs_map ? (ST_int)obs_map[i] : (obs_start + (ST_int)i);
        rc = SF_vstore(freq_idx, obs, freq_weights[i]);
        if (rc) goto cleanup;
    }

    t_store = ctools_timer_seconds() - store_start;

    /* Count statistics */
    size_t n_selected = 0;
    double total_weight = 0.0;
    for (size_t i = 0; i < nobs; i++) {
        if (freq_weights[i] > 0) n_selected++;
        total_weight += freq_weights[i];
    }

cleanup:
    /* No Stata output is written until all strata have valid draw counts. */
    free(thread_rngs);
    if (sort_perm) free(sort_perm);
    free(clusters);
    free(strata);
    if (total_by > 0) {
        ctools_filtered_data_free(&by_filtered);
        free(combined_indices);
    } else if (obs_map) {
        free(obs_map);
    }
    free(freq_weights);
    if (cluster_indices) free(cluster_indices);
    if (strata_indices) free(strata_indices);

    if (rc) {
        if (rc == 498) SF_error("cbsample: resample size exceeds sampling units in a stratum\n");
        return rc;
    }

    /* Store timing scalars */
    double total = ctools_timer_seconds() - t_start;
    SF_scal_save("_cbsample_time_load", t_load);
    SF_scal_save("_cbsample_time_sort", t_sort_cluster);
    SF_scal_save("_cbsample_time_sample", t_sample);
    SF_scal_save("_cbsample_time_store", t_store);
    SF_scal_save("_cbsample_time_total", total);
    SF_scal_save("_cbsample_nobs", (double)nobs);
    SF_scal_save("_cbsample_n_selected", (double)n_selected);
    SF_scal_save("_cbsample_total_weight", total_weight);
    SF_scal_save("_cbsample_nclusters", (double)nclusters);
    SF_scal_save("_cbsample_nstrata", (double)nstrata_groups);

    #ifdef _OPENMP
    SF_scal_save("_cbsample_threads", (double)num_threads);
    #else
    SF_scal_save("_cbsample_threads", 1.0);
    #endif

    if (config.verbose) {
        ctools_msg(CBSAMPLE_MODULE, "%zu obs, %zu clusters, %zu strata, selected %zu, total weight %.0f",
                   nobs, nclusters, nstrata_groups, n_selected, total_weight);
        ctools_msg(CBSAMPLE_MODULE, "load=%.4f sort=%.4f sample=%.4f store=%.4f total=%.4f",
                   t_load, t_sort_cluster, t_sample, t_store, total);
    }

    return 0;
}
