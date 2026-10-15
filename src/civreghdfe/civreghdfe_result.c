#include "civreghdfe_result.h"
#include "ctools_spi.h"
#include "ctools_config.h"
#include <stdio.h>

/* Expand omitted columns and reorder [exog,endog] to Stata's [endog,exog].
 * A failed map allocation must not fall back to indexing the compact V with
 * full-model indices. Checked SPI writes also protect undersized host matrices. */
static ST_retcode store_iv_matrices(const ST_double *beta, const ST_double *V, ST_int K_exog_orig,
                                    ST_int K_endog_orig, ST_int K_total, const ST_int *is_collinear)
{
    ST_int full_K = K_exog_orig + K_endog_orig;
    ST_int *map = (ST_int *)malloc((size_t)full_K * sizeof(ST_int));
    if (!map)
        return 920;
    ST_int compact = 0;
    for (ST_int k = 0; k < full_K; k++)
        map[k] = is_collinear[k] ? -1 : compact++;
    if (compact != K_total)
    {
        free(map);
        return 503;
    }

    ST_retcode rc = 0;
    for (ST_int i = 0; i < full_K && !rc; i++)
    {
        ST_int fi = i < K_endog_orig ? K_exog_orig + i : i - K_endog_orig;
        ST_int ci = map[fi];
        rc = ctools_mat_store("__civreghdfe_b", 1, i + 1, ci < 0 ? 0.0 : beta[ci]);
        for (ST_int j = 0; j < full_K && !rc; j++)
        {
            ST_int fj = j < K_endog_orig ? K_exog_orig + j : j - K_endog_orig;
            ST_int cj = map[fj];
            ST_double value = ci < 0 || cj < 0 ? 0.0 : V[(size_t)cj * K_total + ci];
            rc = ctools_mat_store("__civreghdfe_V", i + 1, j + 1, value);
        }
    }
    free(map);
    return rc;
}
ST_retcode civreghdfe_post_results(const civreghdfe_result *result)
{
    ST_retcode output_rc = 0;
    /* Store results to Stata using checked wrappers */
    /* Scalars - use error-checking wrappers to prevent silent failures */
    if ((output_rc = ctools_scal_save("__civreghdfe_N", result->N_eff)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__civreghdfe_df_r", (ST_double)result->df_r_val)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__civreghdfe_df_a", (ST_double)result->df_a)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__civreghdfe_df_a_for_vce", (ST_double)result->df_a_for_vce)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__civreghdfe_K", (ST_double)(result->K_total + result->n_abs_x))))
        return output_rc;
    if ((output_rc = ctools_scal_save("__civreghdfe_rss", result->rss)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__civreghdfe_tss", result->tss)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__civreghdfe_r2", result->r2)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__civreghdfe_rmse", result->rmse)))
        return output_rc;

    if (result->has_cluster2)
    {
        /* Two-way clustering: N_clust = min(G1, G2), store both counts separately */
        ST_int min_clust = (result->num_clusters < result->num_clusters2) ? result->num_clusters
                                                                          : result->num_clusters2;
        if ((output_rc = ctools_scal_save("__civreghdfe_N_clust", (ST_double)min_clust)))
            return output_rc;
        if ((output_rc =
                 ctools_scal_save("__civreghdfe_N_clust1", (ST_double)result->num_clusters)))
            return output_rc;
    }
    else if (result->has_cluster)
    {
        if ((output_rc = ctools_scal_save("__civreghdfe_N_clust", (ST_double)result->num_clusters)))
            return output_rc;
    }

    if (result->has_cluster2)
    {
        if ((output_rc =
                 ctools_scal_save("__civreghdfe_N_clust2", (ST_double)result->num_clusters2)))
            return output_rc;
    }

    /* Store result->lambda for LIML/Fuller */
    if (result->est_method == 1 || result->est_method == 2)
    {
        if ((output_rc = ctools_scal_save("__civreghdfe_lambda", result->lambda)))
            return output_rc;
    }

    ST_retcode matrix_rc =
        store_iv_matrices(result->beta, result->V, result->K_exog_orig, result->K_endog_orig,
                          result->K_total, result->is_collinear_x);
    if (result->residual_rc)
        return result->residual_rc;
    if (matrix_rc)
        return matrix_rc;

    /* First stage F-stats */
    for (ST_int e = 0; e < result->K_endog; e++)
    {
        char name[64];
        snprintf(name, sizeof(name), "__civreghdfe_F1_%d", (int)(e + 1));
        if ((output_rc = ctools_scal_save(name, result->first_stage_F[e])))
            return output_rc;
    }

    return 0;
}
