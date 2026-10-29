#include "cppmlhdfe_result.h"
#include "ctools_spi.h"
#include "ctools_config.h"
#include <stdio.h>

ST_retcode cppmlhdfe_post_results(const cppmlhdfe_result *result)
{
    ST_int i, j, k, g;
    char scalar_name[64];
    ST_retcode output_rc = 0;
    /* result->N reporting: for fweights, report sum of weights; otherwise report obs count */
    {
        ST_double N_report = (ST_double)result->N;
        if (result->weight_type == 2 && result->has_weights && result->w_user_compact)
        {
            N_report = 0.0;
            for (i = 0; i < result->N; i++)
                N_report += result->w_user_compact[i];
        }
        if ((output_rc = ctools_scal_save("__cppmlhdfe_N", N_report)))
            return output_rc;
    }

    if ((output_rc =
             ctools_scal_save("__cppmlhdfe_num_singletons", (ST_double)result->num_singletons)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__cppmlhdfe_num_separated", (ST_double)result->num_separated)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_K_keep", (ST_double)result->K_keep)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_df_a", (ST_double)result->df_a)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__cppmlhdfe_mobility_groups", (ST_double)result->mobility_groups)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_deviance", result->deviance * result->stdev_y *
                                                                  result->user_weight_scale)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_ll", result->ll)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_ll_0", result->ll_0)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_irls_iterations", (ST_double)result->irls_iter)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__cppmlhdfe_irls_converged", (ST_double)result->irls_converged)))
        return output_rc;

    for (g = 0; g < result->G; g++)
    {
        snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_num_levels_%d", g + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)result->factors[g].num_levels)))
            return output_rc;
    }

    /* Store betas */
    for (k = 0; k < result->K_keep; k++)
    {
        snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_beta_%d", k + 1);
        if ((output_rc =
                 ctools_scal_save(scalar_name, result->beta_final ? result->beta_final[k] : 0.0)))
            return output_rc;
    }

    /* Store VCE matrix */
    if (result->V_keep)
    {
        for (i = 0; i < result->K_keep; i++)
        {
            for (j = 0; j < result->K_keep; j++)
            {
                if ((output_rc = ctools_mat_store("__cppmlhdfe_V", i + 1, j + 1,
                                                  result->V_keep[i * result->K_keep + j])))
                    return output_rc;
            }
        }
    }

    if ((output_rc = ctools_scal_save("__cppmlhdfe_cluster_terms_done", result->cluster_terms)))
        return output_rc;

    /* Store cluster count */
    if (result->vcetype == 2 && result->num_clusters > 0)
    {
        if ((output_rc = ctools_scal_save("__cppmlhdfe_N_clust", (ST_double)result->num_clusters)))
            return output_rc;
    }

    return 0;
}
