#include "creghdfe_result.h"
#include "ctools_spi.h"
#include "ctools_config.h"
#include <stdio.h>

ST_retcode creghdfe_post_results(const creghdfe_result *result)
{
    ST_int i, j, k, g;
    char scalar_name[64];
    ST_retcode output_rc = 0;
    /* HDFE init results */
    /* For fweight, report result->N as sum of weights (like reghdfe) */
    if (result->weight_type == 2 && result->has_weights)
    {
        if ((output_rc = ctools_scal_save("__creghdfe_N", result->g_state->sum_weights)))
            return output_rc;
    }
    else
    {
        if ((output_rc = ctools_scal_save("__creghdfe_N", (ST_double)result->N)))
            return output_rc;
    }
    if ((output_rc =
             ctools_scal_save("__creghdfe_num_singletons", (ST_double)result->num_singletons)))
        return output_rc;
    for (g = 0; g < result->G; g++)
    {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_num_levels_%d", g + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)result->factors[g].num_levels)))
            return output_rc;
    }
    if ((output_rc = ctools_scal_save("__creghdfe_df_a", (ST_double)result->df_a)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__creghdfe_df_a_nested_computed", (ST_double)result->df_a_nested)))
        return output_rc;
    if ((output_rc =
             ctools_scal_save("__creghdfe_mobility_groups", (ST_double)result->mobility_groups)))
        return output_rc;

    /* OLS results */
    /* For fweight, report result->N as sum of weights (like reghdfe) */
    if (result->weight_type == 2 && result->has_weights)
    {
        if ((output_rc = ctools_scal_save("__creghdfe_ols_N", result->g_state->sum_weights)))
            return output_rc;
    }
    else
    {
        if ((output_rc = ctools_scal_save("__creghdfe_ols_N", (ST_double)result->N)))
            return output_rc;
    }
    if ((output_rc = ctools_scal_save("__creghdfe_K_keep", (ST_double)result->K_keep)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__creghdfe_has_cons", 1.0)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__creghdfe_rss", result->rss)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__creghdfe_tss_within", result->tss_within)))
        return output_rc;
    if ((output_rc = ctools_scal_save("__creghdfe_tss", result->tss[0])))
        return output_rc;

    /* Store betas */
    for (k = 0; k < result->K_keep; k++)
    {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_beta_%d", k + 1);
        if ((output_rc = ctools_scal_save(scalar_name, result->beta_keep[k])))
            return output_rc;
    }
    if ((output_rc = ctools_scal_save("__creghdfe_cons", result->beta_keep[result->K_keep])))
        return output_rc;

    /* Store VCE matrix directly */
    for (i = 0; i < result->K_with_cons; i++)
    {
        for (j = 0; j < result->K_with_cons; j++)
        {
            if ((output_rc = ctools_mat_store("__creghdfe_V", i + 1, j + 1,
                                              result->V_keep[i * result->K_with_cons + j])))
                return output_rc;
        }
    }

    /* Store number of clusters */
    if (result->vcetype == 2 && result->num_clusters > 0)
    {
        if ((output_rc = ctools_scal_save("__creghdfe_N_clust", (ST_double)result->num_clusters)))
            return output_rc;
    }

    return 0;
}
