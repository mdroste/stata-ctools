#include "cppmlhdfe_internal.h"

void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N)
{
    ST_int g, i, lev;
    ST_double sum_w = 0.0;

    for (i = 0; i < N; i++)
        sum_w += irls_w[i];
    S->sum_weights = sum_w;

    /* each factor's counts are an independent pass over the rows */
#pragma omp parallel for schedule(static) private(i, lev) if (S->G > 1 && N > 100000)
    for (g = 0; g < S->G; g++)
    {
        FE_Factor *f = &S->factors[g];
        if (f->slope_center)
        {
            memset(f->weighted_counts, 0, f->num_levels * sizeof(ST_double));
            memset(f->slope_center, 0, f->num_levels * sizeof(ST_double));
            for (i = 0; i < N; i++)
            {
                lev = f->levels[i] - 1;
                f->weighted_counts[lev] += irls_w[i];
                f->slope_center[lev] += irls_w[i] * f->slope[i];
            }
            for (i = 0; i < f->num_levels; i++)
                if (f->weighted_counts[i] > 0)
                    f->slope_center[i] /= f->weighted_counts[i];
        }
        memset(f->weighted_counts, 0, f->num_levels * sizeof(ST_double));
        for (i = 0; i < N; i++)
        {
            lev = f->levels[i] - 1;
            ST_double loading =
                f->slope ? f->slope[i] - (f->slope_center ? f->slope_center[lev] : 0.0) : 1.0;
            f->weighted_counts[lev] += irls_w[i] * loading * loading;
        }
        if (f->inv_weighted_counts)
        {
            for (i = 0; i < f->num_levels; i++)
            {
                f->inv_weighted_counts[i] =
                    (f->weighted_counts[i] > 0.0) ? 1.0 / f->weighted_counts[i] : 0.0;
            }
        }
    }
}

ST_double ppml_compute_deviance(const ST_double *y, const ST_double *mu, const ST_double *w_user,
                                ST_int N)
{
    ST_int i;
    ST_double dev = 0.0;

    for (i = 0; i < N; i++)
    {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double d;
        if (y[i] > 0.0)
        {
            d = y[i] * log(y[i] / mu[i]) - (y[i] - mu[i]);
        }
        else
        {
            d = mu[i]; /* 0*log(0/mu) - (0 - mu) = mu */
        }
        dev += w * d;
    }

    return 2.0 * dev;
}

ST_double ppml_compute_loglik(const ST_double *y, const ST_double *mu, const ST_double *w_user,
                              ST_int N)
{
    ST_int i;
    ST_double ll = 0.0;

    for (i = 0; i < N; i++)
    {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double contrib;
        if (y[i] > 0.0)
        {
            contrib = y[i] * log(mu[i]) - mu[i] - lgamma(y[i] + 1.0);
        }
        else
        {
            contrib = -mu[i];
        }
        ll += w * contrib;
    }

    return ll;
}
