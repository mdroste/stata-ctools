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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "ctools_error.h"
#include "ctools_spi.h"
#include "cbsample_impl.h"

#define CBSAMPLE_MODULE "cbsample"

/* ===========================================================================
   High-quality PRNG: xoshiro256** (Blackman & Vigna)
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

static inline size_t xoshiro256_uniform(xoshiro256_state *state, size_t n) {
    if (n == 0) return 0;
    uint64_t threshold = (UINT64_MAX - n + 1) % n;
    uint64_t r;
    do {
        r = xoshiro256_next(state);
    } while (r < threshold);
    return (size_t)(r % n);
}

/* ===========================================================================
   Group Detection
   =========================================================================== */

typedef struct {
    size_t start;
    size_t count;
} group_info;

/* Check if two adjacent observations are in the same group */
static inline int same_group_check(stata_data *data, size_t nvars, size_t i, double miss) {
    for (size_t b = 0; b < nvars; b++) {
        double prev = data->vars[b].data.dbl[i - 1];
        double curr = data->vars[b].data.dbl[i];

        int prev_miss = (prev >= miss);
        int curr_miss = (curr >= miss);

        if (prev_miss && curr_miss) continue;
        if (prev_miss || curr_miss || prev != curr) {
            return 0;
        }
    }
    return 1;
}

/* Detect groups from sorted stata_data (data must already be sorted/permuted)
   Parallelized for large datasets using boundary detection */
static void detect_groups_from_data(stata_data *data, size_t nvars,
                                     group_info *groups, size_t *ngroups) {
    size_t nobs = data->nobs;
    if (nobs == 0) {
        *ngroups = 0;
        return;
    }

    const double miss = SV_missval;

    #ifdef _OPENMP
    int num_threads = omp_get_max_threads();

    /* For small datasets or few threads, use sequential version */
    if (nobs < 10000 || num_threads <= 1) {
    #endif
        /* Sequential version */
        size_t g = 0;
        groups[0].start = 0;
        groups[0].count = 1;

        for (size_t i = 1; i < nobs; i++) {
            if (same_group_check(data, nvars, i, miss)) {
                groups[g].count++;
            } else {
                g++;
                groups[g].start = i;
                groups[g].count = 1;
            }
        }
        *ngroups = g + 1;
        return;
    #ifdef _OPENMP
    }

    /* Parallel version: first find all group boundaries */
    /* Use a bitmap to mark boundaries (where new group starts) */
    uint8_t *boundaries = (uint8_t *)calloc(nobs, sizeof(uint8_t));
    if (!boundaries) {
        /* Fallback to sequential */
        size_t g = 0;
        groups[0].start = 0;
        groups[0].count = 1;
        for (size_t i = 1; i < nobs; i++) {
            if (same_group_check(data, nvars, i, miss)) {
                groups[g].count++;
            } else {
                g++;
                groups[g].start = i;
                groups[g].count = 1;
            }
        }
        *ngroups = g + 1;
        return;
    }

    boundaries[0] = 1;  /* First obs is always a boundary */

    /* Phase 1: Parallel boundary detection */
    #pragma omp parallel for schedule(static)
    for (size_t i = 1; i < nobs; i++) {
        if (!same_group_check(data, nvars, i, miss)) {
            boundaries[i] = 1;
        }
    }

    /* Phase 2: Sequential group construction (must be sequential for ordering) */
    size_t g = 0;
    for (size_t i = 0; i < nobs; i++) {
        if (boundaries[i]) {
            if (g > 0) {
                /* Finalize previous group's count */
                groups[g - 1].count = i - groups[g - 1].start;
            }
            groups[g].start = i;
            g++;
        }
    }
    /* Finalize last group */
    if (g > 0) {
        groups[g - 1].count = nobs - groups[g - 1].start;
    }

    free(boundaries);
    *ngroups = g;
    #endif
}

/* ===========================================================================
   Argument Parsing
   =========================================================================== */

typedef struct {
    size_t n;           /* Number of observations in bootstrap sample (0 = same as original) */
    int verbose;
    uint64_t seed;
} cbsample_config;

static int parse_option(const char *args, const char *name) {
    if (!args || !name) return 0;
    const char *p = strstr(args, name);
    if (!p) return 0;
    if (p > args && *(p-1) != ' ' && *(p-1) != '\t') return 0;
    char after = *(p + strlen(name));
    if (after != '\0' && after != ' ' && after != '\t' && after != '=') return 0;
    return 1;
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

ST_retcode cbsample_main(const char *args) {
    double t_start, t_load, t_sort_cluster, t_sample, t_store;
    int *cluster_indices = NULL;
    int *strata_indices = NULL;
    int *combined_indices = NULL;
    int *sort_vars = NULL;
    double *freq_weights = NULL;
    group_info *clusters = NULL;
    group_info *strata = NULL;
    perm_idx_t *sort_perm = NULL;  /* Saved permutation for mapping back */
    stata_data by_data;
    stata_retcode rc;
    int num_threads = 1;
    xoshiro256_state *thread_rngs = NULL;

    stata_data_init(&by_data);
    t_start = ctools_timer_seconds();
    t_sort_cluster = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CBSAMPLE_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    cbsample_config config;
    config.n = parse_size_option(args, "n", 0);
    config.verbose = parse_option(args, "verbose");
    config.seed = (uint64_t)parse_size_option(args, "seed", 0);

    /* Parse: freq_idx ncluster nstrata cluster_indices... strata_indices... */
    const char *p = args;
    while (*p == ' ' || *p == '\t') p++;

    char *end;
    int freq_idx = (int)strtol(p, &end, 10);
    if (end == p || freq_idx < 1) {
        ctools_error(CBSAMPLE_MODULE, "invalid freq_idx");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    size_t ncluster = (size_t)strtol(p, &end, 10);
    if (end == p) {
        ctools_error(CBSAMPLE_MODULE, "invalid ncluster");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    size_t nstrata = (size_t)strtol(p, &end, 10);
    if (end == p) {
        ctools_error(CBSAMPLE_MODULE, "invalid nstrata");
        return 198;
    }
    p = end;

    if (ncluster > 0) {
        cluster_indices = (int *)malloc(ncluster * sizeof(int));
        if (!cluster_indices) {
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }
        if (parse_int_array(cluster_indices, ncluster, &p) != 0) {
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
        if (parse_int_array(strata_indices, nstrata, &p) != 0) {
            if (cluster_indices) free(cluster_indices);
            free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "failed to parse strata indices");
            return 198;
        }
    }

    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    size_t nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error(CBSAMPLE_MODULE, "no observations");
        return 2000;
    }

    /* Default n to nobs if not specified */
    if (config.n == 0) {
        config.n = nobs;
    }

    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif

    if (config.seed == 0) {
        config.seed = (uint64_t)time(NULL) ^ ((uint64_t)clock() << 32);
    }

    /* === Load Phase using ctools_data_load_selective === */
    double load_start = ctools_timer_seconds();

    freq_weights = (double *)calloc(nobs, sizeof(double));
    if (!freq_weights) {
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    size_t total_by = nstrata + ncluster;
    size_t nclusters = 1;
    size_t nstrata_groups = 1;

    if (total_by > 0) {
        /* Build combined indices array: strata first, then cluster */
        combined_indices = (int *)malloc(total_by * sizeof(int));
        if (!combined_indices) {
            free(freq_weights);
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

        /* Load data using standard ctools function */
        rc = ctools_data_load_selective(&by_data, combined_indices, total_by, 0, 0);
        if (rc != STATA_OK) {
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "failed to load data");
            return 920;
        }
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort and Detect Groups Phase using ctools_sort_dispatch === */
    double sort_start = ctools_timer_seconds();

    /* Allocate groups */
    clusters = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    strata = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    if (!clusters || !strata) {
        if (clusters) free(clusters);
        if (strata) free(strata);
        if (total_by > 0) {
            stata_data_free(&by_data);
            free(combined_indices);
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
            stata_data_free(&by_data);
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
        rc = ctools_sort_dispatch(&by_data, sort_vars, total_by, SORT_ALG_COUNTING);
        if (rc != STATA_OK) {
            free(sort_vars);
            free(clusters);
            free(strata);
            stata_data_free(&by_data);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "sort failed");
            return 920;
        }

        /* Save the sort permutation before applying (for mapping back to original indices) */
        sort_perm = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
        if (!sort_perm) {
            free(sort_vars);
            free(clusters);
            free(strata);
            stata_data_free(&by_data);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }
        memcpy(sort_perm, by_data.sort_order, nobs * sizeof(perm_idx_t));

        /* Apply permutation to physically reorder the data */
        rc = ctools_apply_permutation(&by_data);
        if (rc != STATA_OK) {
            free(sort_perm);
            free(sort_vars);
            free(clusters);
            free(strata);
            stata_data_free(&by_data);
            free(combined_indices);
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error(CBSAMPLE_MODULE, "permutation failed");
            return 920;
        }

        free(sort_vars);
    }

    /* Detect strata groups (using first nstrata variables) */
    if (nstrata > 0) {
        /* Create a temporary view with just strata vars for group detection */
        stata_data strata_view;
        strata_view.nobs = by_data.nobs;
        strata_view.nvars = nstrata;
        strata_view.vars = by_data.vars;  /* First nstrata vars are strata */
        strata_view.sort_order = NULL;
        detect_groups_from_data(&strata_view, nstrata, strata, &nstrata_groups);
    } else {
        strata[0].start = 0;
        strata[0].count = nobs;
        nstrata_groups = 1;
    }

    /* Detect clusters (using cluster variables at offset nstrata) */
    if (ncluster > 0) {
        /* Create a temporary view with just cluster vars */
        stata_data cluster_view;
        cluster_view.nobs = by_data.nobs;
        cluster_view.nvars = ncluster;
        cluster_view.vars = &by_data.vars[nstrata];  /* Cluster vars start after strata */
        cluster_view.sort_order = NULL;
        detect_groups_from_data(&cluster_view, ncluster, clusters, &nclusters);
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
    thread_rngs = (xoshiro256_state *)malloc(num_threads * sizeof(xoshiro256_state));
    if (!thread_rngs) {
        if (sort_perm) free(sort_perm);
        free(clusters);
        free(strata);
        if (total_by > 0) {
            stata_data_free(&by_data);
            free(combined_indices);
        }
        free(freq_weights);
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    for (int t = 0; t < num_threads; t++) {
        xoshiro256_seed(&thread_rngs[t], config.seed + (uint64_t)t * 0x9e3779b97f4a7c15ULL);
    }

    /* Process each stratum */
    if (ncluster > 0) {
        /* Cluster bootstrap: resample clusters with replacement within each stratum */
        /* For each stratum, identify which clusters are in it and sample from those */

        for (size_t st = 0; st < nstrata_groups; st++) {
            size_t st_start = strata[st].start;
            size_t st_count = strata[st].count;

            /* Find clusters within this stratum */
            /* Clusters are sorted by (strata, cluster), so find range */
            size_t first_cluster = 0;
            size_t last_cluster = nclusters;

            /* Binary search for first cluster in stratum */
            for (size_t cl = 0; cl < nclusters; cl++) {
                if (clusters[cl].start >= st_start &&
                    clusters[cl].start < st_start + st_count) {
                    first_cluster = cl;
                    break;
                }
            }
            /* Find last cluster in stratum */
            for (size_t cl = first_cluster; cl < nclusters; cl++) {
                if (clusters[cl].start >= st_start + st_count) {
                    last_cluster = cl;
                    break;
                }
            }

            size_t n_clusters_in_stratum = last_cluster - first_cluster;
            if (n_clusters_in_stratum == 0) continue;

            /* Calculate target sample size for this stratum (proportional) */
            size_t target_n = (size_t)round((double)config.n * (double)st_count / (double)nobs);
            if (target_n == 0) target_n = 1;

            /* Sample clusters with replacement */
            xoshiro256_state *rng = &thread_rngs[0];  /* Use single thread for now */

            /* We need to sample enough clusters to get target_n observations */
            /* Simple approach: sample n_clusters_in_stratum clusters */
            for (size_t draw = 0; draw < n_clusters_in_stratum; draw++) {
                size_t cluster_idx = first_cluster + xoshiro256_uniform(rng, n_clusters_in_stratum);
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
        #ifdef _OPENMP
        if (nstrata_groups == 1 && num_threads > 1 && config.n >= 1000) {
            /* Single stratum with many draws: use thread-local arrays to avoid atomics */
            /* Each thread accumulates into its own array, then we merge */
            double **thread_weights = (double **)malloc(num_threads * sizeof(double *));
            if (!thread_weights) {
                /* Fallback to atomic version */
                goto atomic_fallback;
            }

            int alloc_ok = 1;
            for (int t = 0; t < num_threads; t++) {
                thread_weights[t] = (double *)calloc(nobs, sizeof(double));
                if (!thread_weights[t]) {
                    for (int j = 0; j < t; j++) free(thread_weights[j]);
                    free(thread_weights);
                    alloc_ok = 0;
                    break;
                }
            }

            if (alloc_ok) {
                size_t target_n = config.n;
                size_t draws_per_thread = target_n / num_threads;
                size_t remainder = target_n % num_threads;

                #pragma omp parallel
                {
                    int tid = omp_get_thread_num();
                    xoshiro256_state *rng = &thread_rngs[tid];
                    double *local_weights = thread_weights[tid];
                    size_t my_draws = draws_per_thread + (tid < (int)remainder ? 1 : 0);

                    for (size_t draw = 0; draw < my_draws; draw++) {
                        size_t local_idx = xoshiro256_uniform(rng, nobs);
                        size_t obs_idx = sort_perm ? sort_perm[local_idx] : local_idx;
                        local_weights[obs_idx] += 1.0;
                    }
                }

                /* Merge thread-local arrays into freq_weights */
                #pragma omp parallel for
                for (size_t i = 0; i < nobs; i++) {
                    double sum = 0.0;
                    for (int t = 0; t < num_threads; t++) {
                        sum += thread_weights[t][i];
                    }
                    freq_weights[i] = sum;
                }

                for (int t = 0; t < num_threads; t++) {
                    free(thread_weights[t]);
                }
                free(thread_weights);
            } else {
                goto atomic_fallback;
            }
        } else {
            atomic_fallback:
            /* Multiple strata or single thread: process strata in parallel with atomics */
            #pragma omp parallel for schedule(dynamic, 1)
            for (size_t st = 0; st < nstrata_groups; st++) {
                int tid = omp_get_thread_num();
                size_t st_start = strata[st].start;
                size_t st_count = strata[st].count;

                size_t target_n = (size_t)round((double)config.n * (double)st_count / (double)nobs);
                if (target_n == 0) target_n = 1;

                xoshiro256_state *rng = &thread_rngs[tid];

                for (size_t draw = 0; draw < target_n; draw++) {
                    size_t local_idx = xoshiro256_uniform(rng, st_count);
                    size_t obs_idx = sort_perm ? sort_perm[st_start + local_idx] : (st_start + local_idx);

                    #pragma omp atomic
                    freq_weights[obs_idx] += 1.0;
                }
            }
        }
        #else
        /* No OpenMP: simple sequential sampling */
        for (size_t st = 0; st < nstrata_groups; st++) {
            size_t st_start = strata[st].start;
            size_t st_count = strata[st].count;

            size_t target_n = (size_t)round((double)config.n * (double)st_count / (double)nobs);
            if (target_n == 0) target_n = 1;

            xoshiro256_state *rng = &thread_rngs[0];

            for (size_t draw = 0; draw < target_n; draw++) {
                size_t local_idx = xoshiro256_uniform(rng, st_count);
                size_t obs_idx = sort_perm ? sort_perm[st_start + local_idx] : (st_start + local_idx);
                freq_weights[obs_idx] += 1.0;
            }
        }
        #endif
    }

    t_sample = ctools_timer_seconds() - sample_start;

    /* === Store Phase === */
    double store_start = ctools_timer_seconds();

    /* Write frequency weights to Stata - use batched writes for speed */
    /* Check if we have contiguous observations (no if condition filtering) */
    /* If nobs == range size, then no observations were filtered by if condition */
    int contiguous = (nobs == (size_t)(obs2 - obs1 + 1));

    if (contiguous) {
        /* Fast path: no if condition, use batched writes */
        ST_int stata_var = (ST_int)freq_idx;
        size_t i = 0;
        size_t i_end_16 = nobs - (nobs % 16);

        /* Batched writes of 16 values */
        for (; i < i_end_16; i += 16) {
            SF_VSTORE_BATCH16(stata_var, obs1 + i, &freq_weights[i]);
        }

        /* Handle remaining elements */
        for (; i < nobs; i++) {
            SF_vstore(stata_var, (ST_int)(obs1 + i), freq_weights[i]);
        }
    } else {
        /* Slow path: if condition present, check each observation */
        size_t obs_idx = 0;
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (!SF_ifobs(obs)) continue;
            SF_vstore(freq_idx, obs, freq_weights[obs_idx]);
            obs_idx++;
        }
    }

    t_store = ctools_timer_seconds() - store_start;

    /* Count statistics */
    size_t n_selected = 0;
    double total_weight = 0.0;
    for (size_t i = 0; i < nobs; i++) {
        if (freq_weights[i] > 0) n_selected++;
        total_weight += freq_weights[i];
    }

    /* === Cleanup === */
    free(thread_rngs);
    if (sort_perm) free(sort_perm);
    free(clusters);
    free(strata);
    if (total_by > 0) {
        stata_data_free(&by_data);
        free(combined_indices);
    }
    free(freq_weights);
    if (cluster_indices) free(cluster_indices);
    if (strata_indices) free(strata_indices);

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
