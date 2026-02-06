/*
 * cqreg_hdfe.c
 *
 * HDFE integration for quantile regression.
 * Uses creghdfe's CG solver for partialling out fixed effects.
 * Part of the ctools suite.
 */

#include "cqreg_hdfe.h"
#include "cqreg_linalg.h"
#include "../creghdfe/creghdfe_types.h"
#include "../creghdfe/creghdfe_solver.h"
#include "../creghdfe/creghdfe_hdfe.h"
#include "../ctools_runtime.h"
#include "../ctools_config.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * Hash Table for Factor Creation
 * ============================================================================ */

#define CQREG_HASH_EMPTY -1
#define CQREG_HASH_LOAD 0.7

typedef struct {
    ST_int *keys;
    ST_int *values;
    ST_int capacity;
    ST_int size;
} cqreg_hash_table;

static cqreg_hash_table *hash_create(ST_int capacity)
{
    cqreg_hash_table *ht = (cqreg_hash_table *)malloc(sizeof(cqreg_hash_table));
    if (ht == NULL) return NULL;

    ht->capacity = capacity;
    ht->size = 0;
    ht->keys = (ST_int *)malloc(capacity * sizeof(ST_int));
    ht->values = (ST_int *)malloc(capacity * sizeof(ST_int));

    if (ht->keys == NULL || ht->values == NULL) {
        free(ht->keys);
        free(ht->values);
        free(ht);
        return NULL;
    }

    for (ST_int i = 0; i < capacity; i++) {
        ht->keys[i] = CQREG_HASH_EMPTY;
    }

    return ht;
}

static void hash_free(cqreg_hash_table *ht)
{
    if (ht == NULL) return;
    free(ht->keys);
    free(ht->values);
    free(ht);
}

static ST_int hash_function(ST_int key, ST_int capacity)
{
    /* Multiplicative hash */
    unsigned int h = (unsigned int)key * 2654435761u;
    return (ST_int)(h % (unsigned int)capacity);
}

static ST_int hash_insert(cqreg_hash_table *ht, ST_int key)
{
    ST_int idx = hash_function(key, ht->capacity);

    while (ht->keys[idx] != CQREG_HASH_EMPTY) {
        if (ht->keys[idx] == key) {
            return ht->values[idx];  /* Already exists */
        }
        idx = (idx + 1) % ht->capacity;
    }

    /* Insert new entry */
    ht->keys[idx] = key;
    ht->values[idx] = ht->size + 1;  /* 1-indexed levels */
    ht->size++;

    return ht->values[idx];
}

/* ============================================================================
 * Factor Creation
 * ============================================================================ */

ST_int cqreg_hdfe_create_factor(ST_int var_idx,
                                ST_int in1, ST_int in2,
                                ST_int N_filtered,
                                ST_int *levels,
                                ST_double *counts,
                                ST_int *num_levels)
{
    ST_int idx, obs;

    /* Create hash table with room for all unique values */
    ST_int capacity = (ST_int)(N_filtered / CQREG_HASH_LOAD) + 1;
    cqreg_hash_table *ht = hash_create(capacity);
    if (ht == NULL) return -1;

    /* Read variable and build level assignments, filtering by if-condition */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        /* Skip observations that don't pass the if-condition */
        if (!SF_ifobs(obs)) continue;

        ST_double val;
        if (SF_vdata(var_idx, obs, &val) != 0) {
            hash_free(ht);
            return -2;
        }

        /* Convert to integer (FE variables should be integer-valued) */
        ST_int int_val = (ST_int)val;
        levels[idx] = hash_insert(ht, int_val);
        idx++;
    }

    /* Verify we got the expected number of observations */
    if (idx != N_filtered) {
        hash_free(ht);
        return -3;  /* Observation count mismatch */
    }

    *num_levels = ht->size;

    /* Count observations per level */
    memset(counts, 0, ht->size * sizeof(ST_double));
    for (idx = 0; idx < N_filtered; idx++) {
        counts[levels[idx] - 1] += 1.0;  /* levels are 1-indexed */
    }

    hash_free(ht);
    return 0;
}

/* ============================================================================
 * Singleton Detection
 * ============================================================================ */

ST_int cqreg_hdfe_detect_singletons(cqreg_state *state, ST_int *mask)
{
    HDFE_State *hdfe = (HDFE_State *)state->hdfe_state;
    if (hdfe == NULL) return 0;

    ST_int N = hdfe->N;
    ST_int G = hdfe->G;
    ST_int total_dropped = 0;
    ST_int iter = 0;
    ST_int max_iter = 100;

    /* Initialize mask to keep all */
    for (ST_int i = 0; i < N; i++) {
        mask[i] = 1;
    }

    /* Iteratively find and drop singletons */
    while (iter < max_iter) {
        ST_int found = 0;

        for (ST_int g = 0; g < G; g++) {
            FE_Factor *f = &hdfe->factors[g];

            /* Reset counts */
            memset(f->counts, 0, f->num_levels * sizeof(ST_double));

            /* Count non-dropped observations per level */
            for (ST_int i = 0; i < N; i++) {
                if (mask[i]) {
                    ST_int level = f->levels[i] - 1;  /* 0-indexed */
                    f->counts[level] += 1.0;
                }
            }

            /* Mark singletons for dropping */
            for (ST_int i = 0; i < N; i++) {
                if (mask[i]) {
                    ST_int level = f->levels[i] - 1;
                    if (f->counts[level] <= 1.0) {
                        mask[i] = 0;
                        found++;
                        total_dropped++;
                    }
                }
            }
        }

        if (found == 0) break;
        iter++;
    }

    return total_dropped;
}

/* ============================================================================
 * HDFE Initialization
 * ============================================================================ */

ST_int cqreg_hdfe_init(cqreg_state *state,
                       const ST_int *fe_vars,
                       ST_int G,
                       ST_int N_filtered,
                       ST_int in1, ST_int in2,
                       ST_int maxiter,
                       ST_double tolerance)
{
    ST_int g, t;
    HDFE_State *hdfe = NULL;

    if (G <= 0 || fe_vars == NULL || N_filtered <= 0) {
        return -1;
    }

    /* Allocate HDFE state */
    hdfe = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    if (hdfe == NULL) return -1;

    hdfe->G = G;
    hdfe->N = N_filtered;  /* Use filtered observation count */
    hdfe->in1 = in1;
    hdfe->in2 = in2;
    hdfe->has_weights = 0;
    hdfe->weights = NULL;
    hdfe->maxiter = (maxiter > 0) ? maxiter : 10000;
    hdfe->tolerance = (tolerance > 0) ? tolerance : 1e-8;
    hdfe->verbose = 0;

    /* Allocate factors */
    hdfe->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    if (hdfe->factors == NULL) {
        free(hdfe);
        return -1;
    }

    /* Create each factor */
    for (g = 0; g < G; g++) {
        FE_Factor *f = &hdfe->factors[g];

        /* Allocate level assignments for filtered observations */
        f->levels = (ST_int *)malloc(N_filtered * sizeof(ST_int));
        if (f->levels == NULL) {
            goto cleanup;
        }

        /* Temporary counts (will be resized) */
        ST_double *tmp_counts = (ST_double *)calloc(N_filtered, sizeof(ST_double));
        if (tmp_counts == NULL) {
            goto cleanup;
        }

        /* Create factor from Stata variable, filtering by if-condition */
        ST_int num_levels;
        if (cqreg_hdfe_create_factor(fe_vars[g], in1, in2, N_filtered,
                                     f->levels, tmp_counts, &num_levels) != 0) {
            free(tmp_counts);
            goto cleanup;
        }

        f->num_levels = num_levels;
        f->has_intercept = 1;
        f->num_slopes = 0;

        /* Allocate proper-sized counts */
        f->counts = (ST_double *)malloc(num_levels * sizeof(ST_double));
        f->means = (ST_double *)calloc(num_levels, sizeof(ST_double));
        if (f->counts == NULL || f->means == NULL) {
            free(tmp_counts);
            goto cleanup;
        }
        memcpy(f->counts, tmp_counts, num_levels * sizeof(ST_double));
        free(tmp_counts);

        /* Compute inverse counts for fast division in projection */
        f->inv_counts = (ST_double *)malloc(num_levels * sizeof(ST_double));
        if (f->inv_counts) {
            for (ST_int lev = 0; lev < num_levels; lev++) {
                f->inv_counts[lev] = (f->counts[lev] > 0) ? 1.0 / f->counts[lev] : 0.0;
            }
        }
        f->inv_weighted_counts = NULL;  /* cqreg doesn't use weights */
        f->weighted_counts = NULL;
    }

    /* Determine thread count using unified policy */
    hdfe->num_threads = ctools_get_max_threads();

    /* Allocate per-thread CG buffers */
    hdfe->thread_cg_r = (ST_double **)calloc(hdfe->num_threads, sizeof(ST_double *));
    hdfe->thread_cg_u = (ST_double **)calloc(hdfe->num_threads, sizeof(ST_double *));
    hdfe->thread_cg_v = (ST_double **)calloc(hdfe->num_threads, sizeof(ST_double *));
    hdfe->thread_proj = (ST_double **)calloc(hdfe->num_threads, sizeof(ST_double *));
    hdfe->thread_fe_means = (ST_double **)calloc(hdfe->num_threads * G, sizeof(ST_double *));

    if (hdfe->thread_cg_r == NULL || hdfe->thread_cg_u == NULL ||
        hdfe->thread_cg_v == NULL || hdfe->thread_proj == NULL ||
        hdfe->thread_fe_means == NULL) {
        goto cleanup;
    }

    for (t = 0; t < hdfe->num_threads; t++) {
        hdfe->thread_cg_r[t] = (ST_double *)ctools_cacheline_alloc(N_filtered * sizeof(ST_double));
        hdfe->thread_cg_u[t] = (ST_double *)ctools_cacheline_alloc(N_filtered * sizeof(ST_double));
        hdfe->thread_cg_v[t] = (ST_double *)ctools_cacheline_alloc(N_filtered * sizeof(ST_double));
        hdfe->thread_proj[t] = (ST_double *)ctools_cacheline_alloc(N_filtered * sizeof(ST_double));

        if (hdfe->thread_cg_r[t] == NULL || hdfe->thread_cg_u[t] == NULL ||
            hdfe->thread_cg_v[t] == NULL || hdfe->thread_proj[t] == NULL) {
            goto cleanup;
        }

        for (g = 0; g < G; g++) {
            ST_int nl = hdfe->factors[g].num_levels;
            hdfe->thread_fe_means[t * G + g] = (ST_double *)calloc(nl, sizeof(ST_double));
            if (hdfe->thread_fe_means[t * G + g] == NULL) {
                goto cleanup;
            }
        }
    }

    /* Compute degrees of freedom absorbed */
    hdfe->df_a = 0;
    for (g = 0; g < G; g++) {
        hdfe->df_a += hdfe->factors[g].num_levels - 1;
    }
    hdfe->mobility_groups = 1;  /* Will be computed later if G >= 2 */

    hdfe->factors_initialized = 1;

    /* Store in state */
    state->hdfe_state = hdfe;
    state->has_hdfe = 1;
    state->G = G;

    return 0;

cleanup:
    /* Clean up on error */
    if (hdfe != NULL) {
        if (hdfe->factors != NULL) {
            for (g = 0; g < G; g++) {
                free(hdfe->factors[g].levels);
                free(hdfe->factors[g].counts);
                if (hdfe->factors[g].inv_counts) free(hdfe->factors[g].inv_counts);
                free(hdfe->factors[g].means);
            }
            free(hdfe->factors);
        }

        if (hdfe->thread_cg_r != NULL) {
            for (t = 0; t < hdfe->num_threads; t++) {
                ctools_aligned_free(hdfe->thread_cg_r[t]);
            }
            free(hdfe->thread_cg_r);
        }
        if (hdfe->thread_cg_u != NULL) {
            for (t = 0; t < hdfe->num_threads; t++) {
                ctools_aligned_free(hdfe->thread_cg_u[t]);
            }
            free(hdfe->thread_cg_u);
        }
        if (hdfe->thread_cg_v != NULL) {
            for (t = 0; t < hdfe->num_threads; t++) {
                ctools_aligned_free(hdfe->thread_cg_v[t]);
            }
            free(hdfe->thread_cg_v);
        }
        if (hdfe->thread_proj != NULL) {
            for (t = 0; t < hdfe->num_threads; t++) {
                ctools_aligned_free(hdfe->thread_proj[t]);
            }
            free(hdfe->thread_proj);
        }
        if (hdfe->thread_fe_means != NULL) {
            for (t = 0; t < hdfe->num_threads * G; t++) {
                free(hdfe->thread_fe_means[t]);
            }
            free(hdfe->thread_fe_means);
        }

        free(hdfe);
    }

    return -1;
}

/* ============================================================================
 * Partial Out
 * ============================================================================ */

ST_int cqreg_hdfe_partial_out_column(cqreg_state *state,
                                     ST_double *col,
                                     ST_int N,
                                     ST_int thread_id)
{
    (void)N;  /* Unused - HDFE state knows N internally */
    HDFE_State *hdfe = (HDFE_State *)state->hdfe_state;
    if (hdfe == NULL || !hdfe->factors_initialized) {
        return -1;
    }

    /* Use the creghdfe CG solver */
    return cg_solve_column_threaded(hdfe, col, thread_id);
}

ST_int cqreg_hdfe_partial_out(cqreg_state *state,
                              ST_double *y,
                              ST_double *X,
                              ST_int N, ST_int K)
{
    HDFE_State *hdfe = (HDFE_State *)state->hdfe_state;
    if (hdfe == NULL || !hdfe->factors_initialized) {
        return -1;
    }

    ST_int k;
    ST_int rc = 0;

    /* Partial out y first (single-threaded) */
    ST_int iter_y = cg_solve_column_threaded(hdfe, y, 0);
    if (iter_y < 0) {
        ctools_error("cqreg", "HDFE: CG solver failed for y");
        return -1;
    }

    /* Partial out X columns in parallel */
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(K > 1)
#endif
    for (k = 0; k < K; k++) {
        ST_int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        if (tid >= hdfe->num_threads) tid = hdfe->num_threads - 1;
#endif
        ST_int iter = cg_solve_column_threaded(hdfe, &X[k * N], tid);
        if (iter < 0) {
            rc = -1;
        }
    }

    return rc;
}

/* ============================================================================
 * Accessors
 * ============================================================================ */

ST_int cqreg_hdfe_get_df_absorbed(const cqreg_state *state)
{
    if (state == NULL || state->hdfe_state == NULL) {
        return 0;
    }

    HDFE_State *hdfe = (HDFE_State *)state->hdfe_state;

    /* df_a = sum of (num_levels - 1) - (mobility_groups - 1) if G >= 2 */
    ST_int df = hdfe->df_a;
    if (hdfe->G >= 2 && hdfe->mobility_groups > 1) {
        df -= (hdfe->mobility_groups - 1);
    }

    return df;
}

ST_int cqreg_hdfe_get_num_singletons(const cqreg_state *state)
{
    if (state == NULL) return 0;
    return state->N_original - state->N;
}

/* ============================================================================
 * Cleanup
 * ============================================================================ */

void cqreg_hdfe_cleanup(cqreg_state *state)
{
    if (state == NULL || state->hdfe_state == NULL) {
        return;
    }

    HDFE_State *hdfe = (HDFE_State *)state->hdfe_state;
    ST_int g, t;

    /* Validate state before cleanup to prevent crashes from corrupted data */
    if (hdfe->G < 0 || hdfe->G > 1000000 ||
        hdfe->num_threads < 0 || hdfe->num_threads > 10000) {
        /* State appears corrupted - just NULL out and return */
        state->hdfe_state = NULL;
        state->has_hdfe = 0;
        return;
    }

    /* Free factors */
    if (hdfe->factors != NULL) {
        for (g = 0; g < hdfe->G; g++) {
            if (hdfe->factors[g].levels != NULL) free(hdfe->factors[g].levels);
            if (hdfe->factors[g].counts != NULL) free(hdfe->factors[g].counts);
            if (hdfe->factors[g].inv_counts != NULL) free(hdfe->factors[g].inv_counts);
            if (hdfe->factors[g].means != NULL) free(hdfe->factors[g].means);
            /* Clear pointers after free to detect double-free */
            hdfe->factors[g].levels = NULL;
            hdfe->factors[g].counts = NULL;
            hdfe->factors[g].inv_counts = NULL;
            hdfe->factors[g].means = NULL;
        }
        free(hdfe->factors);
        hdfe->factors = NULL;
    }

    /* Free per-thread buffers */
    if (hdfe->thread_cg_r != NULL) {
        for (t = 0; t < hdfe->num_threads; t++) {
            if (hdfe->thread_cg_r[t] != NULL) {
                ctools_aligned_free(hdfe->thread_cg_r[t]);
                hdfe->thread_cg_r[t] = NULL;
            }
        }
        free(hdfe->thread_cg_r);
        hdfe->thread_cg_r = NULL;
    }
    if (hdfe->thread_cg_u != NULL) {
        for (t = 0; t < hdfe->num_threads; t++) {
            if (hdfe->thread_cg_u[t] != NULL) {
                ctools_aligned_free(hdfe->thread_cg_u[t]);
                hdfe->thread_cg_u[t] = NULL;
            }
        }
        free(hdfe->thread_cg_u);
        hdfe->thread_cg_u = NULL;
    }
    if (hdfe->thread_cg_v != NULL) {
        for (t = 0; t < hdfe->num_threads; t++) {
            if (hdfe->thread_cg_v[t] != NULL) {
                ctools_aligned_free(hdfe->thread_cg_v[t]);
                hdfe->thread_cg_v[t] = NULL;
            }
        }
        free(hdfe->thread_cg_v);
        hdfe->thread_cg_v = NULL;
    }
    if (hdfe->thread_proj != NULL) {
        for (t = 0; t < hdfe->num_threads; t++) {
            if (hdfe->thread_proj[t] != NULL) {
                ctools_aligned_free(hdfe->thread_proj[t]);
                hdfe->thread_proj[t] = NULL;
            }
        }
        free(hdfe->thread_proj);
        hdfe->thread_proj = NULL;
    }
    if (hdfe->thread_fe_means != NULL) {
        ST_int total_fe_means = hdfe->num_threads * hdfe->G;
        for (t = 0; t < total_fe_means; t++) {
            if (hdfe->thread_fe_means[t] != NULL) {
                free(hdfe->thread_fe_means[t]);
                hdfe->thread_fe_means[t] = NULL;
            }
        }
        free(hdfe->thread_fe_means);
        hdfe->thread_fe_means = NULL;
    }

    if (hdfe->weights != NULL) {
        free(hdfe->weights);
        hdfe->weights = NULL;
    }
    free(hdfe);

    state->hdfe_state = NULL;
    state->has_hdfe = 0;
}
