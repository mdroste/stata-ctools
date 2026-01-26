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
    double **cluster_data = NULL;
    double **strata_data = NULL;
    double *freq_weights = NULL;
    group_info *clusters = NULL;
    group_info *strata = NULL;
    size_t *sort_indices = NULL;
    int num_threads = 1;
    xoshiro256_state *thread_rngs = NULL;

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

    /* === Load Phase === */
    double load_start = ctools_timer_seconds();

    freq_weights = (double *)calloc(nobs, sizeof(double));
    if (!freq_weights) {
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    /* Load cluster variables */
    if (ncluster > 0) {
        cluster_data = (double **)malloc(ncluster * sizeof(double *));
        if (!cluster_data) {
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }

        for (size_t c = 0; c < ncluster; c++) {
            cluster_data[c] = (double *)malloc(nobs * sizeof(double));
            if (!cluster_data[c]) {
                for (size_t j = 0; j < c; j++) free(cluster_data[j]);
                free(cluster_data);
                free(freq_weights);
                if (cluster_indices) free(cluster_indices);
                if (strata_indices) free(strata_indices);
                ctools_error_alloc(CBSAMPLE_MODULE);
                return 920;
            }

            int idx = cluster_indices[c];
            for (size_t i = 0; i < nobs; i++) {
                double val;
                if (SF_vdata(idx, (ST_int)(obs1 + i), &val) != 0) val = SV_missval;
                cluster_data[c][i] = val;
            }
        }
    }

    /* Load strata variables */
    if (nstrata > 0) {
        strata_data = (double **)malloc(nstrata * sizeof(double *));
        if (!strata_data) {
            if (cluster_data) {
                for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
                free(cluster_data);
            }
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }

        for (size_t s = 0; s < nstrata; s++) {
            strata_data[s] = (double *)malloc(nobs * sizeof(double));
            if (!strata_data[s]) {
                for (size_t j = 0; j < s; j++) free(strata_data[j]);
                free(strata_data);
                if (cluster_data) {
                    for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
                    free(cluster_data);
                }
                free(freq_weights);
                if (cluster_indices) free(cluster_indices);
                if (strata_indices) free(strata_indices);
                ctools_error_alloc(CBSAMPLE_MODULE);
                return 920;
            }

            int idx = strata_indices[s];
            for (size_t i = 0; i < nobs; i++) {
                double val;
                if (SF_vdata(idx, (ST_int)(obs1 + i), &val) != 0) val = SV_missval;
                strata_data[s][i] = val;
            }
        }
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort and Detect Groups Phase === */
    double sort_start = ctools_timer_seconds();

    size_t nclusters = 1;
    size_t nstrata_groups = 1;

    /* Allocate groups */
    clusters = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    strata = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    if (!clusters || !strata) {
        if (clusters) free(clusters);
        if (strata) free(strata);
        if (strata_data) {
            for (size_t s = 0; s < nstrata; s++) free(strata_data[s]);
            free(strata_data);
        }
        if (cluster_data) {
            for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
            free(cluster_data);
        }
        free(freq_weights);
        if (cluster_indices) free(cluster_indices);
        if (strata_indices) free(strata_indices);
        ctools_error_alloc(CBSAMPLE_MODULE);
        return 920;
    }

    /* Sort by combined (strata, cluster) to detect nested groups */
    /* First, sort by strata + cluster together */
    size_t total_by = nstrata + ncluster;
    double **combined_data = NULL;

    if (total_by > 0) {
        combined_data = (double **)malloc(total_by * sizeof(double *));
        if (!combined_data) {
            free(clusters);
            free(strata);
            if (strata_data) {
                for (size_t s = 0; s < nstrata; s++) free(strata_data[s]);
                free(strata_data);
            }
            if (cluster_data) {
                for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
                free(cluster_data);
            }
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }

        /* Strata first, then cluster */
        for (size_t s = 0; s < nstrata; s++) {
            combined_data[s] = strata_data[s];
        }
        for (size_t c = 0; c < ncluster; c++) {
            combined_data[nstrata + c] = cluster_data[c];
        }

        sort_indices = (size_t *)malloc(nobs * sizeof(size_t));
        if (!sort_indices) {
            free(combined_data);
            free(clusters);
            free(strata);
            if (strata_data) {
                for (size_t s = 0; s < nstrata; s++) free(strata_data[s]);
                free(strata_data);
            }
            if (cluster_data) {
                for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
                free(cluster_data);
            }
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }

        argsort_by_vars(combined_data, nobs, total_by, sort_indices);

        /* Reorder all by_data arrays */
        double *temp = (double *)malloc(nobs * sizeof(double));
        if (!temp) {
            free(sort_indices);
            free(combined_data);
            free(clusters);
            free(strata);
            if (strata_data) {
                for (size_t s = 0; s < nstrata; s++) free(strata_data[s]);
                free(strata_data);
            }
            if (cluster_data) {
                for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
                free(cluster_data);
            }
            free(freq_weights);
            if (cluster_indices) free(cluster_indices);
            if (strata_indices) free(strata_indices);
            ctools_error_alloc(CBSAMPLE_MODULE);
            return 920;
        }

        for (size_t s = 0; s < nstrata; s++) {
            for (size_t i = 0; i < nobs; i++) {
                temp[i] = strata_data[s][sort_indices[i]];
            }
            memcpy(strata_data[s], temp, nobs * sizeof(double));
        }
        for (size_t c = 0; c < ncluster; c++) {
            for (size_t i = 0; i < nobs; i++) {
                temp[i] = cluster_data[c][sort_indices[i]];
            }
            memcpy(cluster_data[c], temp, nobs * sizeof(double));
        }
        free(temp);
        free(combined_data);
    }

    /* Detect strata groups (using only strata variables) */
    if (nstrata > 0) {
        detect_groups(strata_data, nobs, nstrata, strata, &nstrata_groups);
    } else {
        strata[0].start = 0;
        strata[0].count = nobs;
        nstrata_groups = 1;
    }

    /* Detect clusters within each stratum (using cluster variables) */
    if (ncluster > 0) {
        /* Since data is sorted by (strata, cluster), we detect clusters globally */
        detect_groups(cluster_data, nobs, ncluster, clusters, &nclusters);
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
        if (sort_indices) free(sort_indices);
        free(clusters);
        free(strata);
        if (strata_data) {
            for (size_t s = 0; s < nstrata; s++) free(strata_data[s]);
            free(strata_data);
        }
        if (cluster_data) {
            for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
            free(cluster_data);
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
                    size_t obs_idx = sort_indices ? sort_indices[cl_start + j] : (cl_start + j);
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
                        size_t obs_idx = sort_indices ? sort_indices[local_idx] : local_idx;
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
                    size_t obs_idx = sort_indices ? sort_indices[st_start + local_idx] : (st_start + local_idx);

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
                size_t obs_idx = sort_indices ? sort_indices[st_start + local_idx] : (st_start + local_idx);
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
    if (sort_indices) free(sort_indices);
    free(clusters);
    free(strata);
    if (strata_data) {
        for (size_t s = 0; s < nstrata; s++) free(strata_data[s]);
        free(strata_data);
    }
    if (cluster_data) {
        for (size_t c = 0; c < ncluster; c++) free(cluster_data[c]);
        free(cluster_data);
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
