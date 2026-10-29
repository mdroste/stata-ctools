#include <math.h>
#include <string.h>
#include "ctools_ols.h"
#include "cpsmatch_stats.h"

typedef struct
{
    size_t n;
    double anchor, mean;
    dd_real shifted_sum, squares;
} cpsmatch_moment;

static int cpsmatch_stat_nonmissing(double x)
{
    return isfinite(x) && x < SV_missval;
}

static double cpsmatch_stat_value(double x)
{
    return cpsmatch_stat_nonmissing(x) ? x : SV_missval;
}

static void cpsmatch_moment_add(cpsmatch_moment *moment, double x)
{
    if (moment->n == 0)
        moment->anchor = x;
    moment->shifted_sum = dd_add_d(moment->shifted_sum, x - moment->anchor);
    moment->n++;
}

static void cpsmatch_moment_deviation(cpsmatch_moment *moment, double x)
{
    double d = x - moment->mean;
    dd_real square = two_prod(d, d);
    moment->squares = dd_add_d(moment->squares, square.hi);
    moment->squares = dd_add_d(moment->squares, square.lo);
}

static double cpsmatch_moment_variance(const cpsmatch_moment *moment)
{
    if (moment->n < 2)
        return SV_missval;
    /* Match summarize's squared deviations from its double-precision mean.
     * Re-centering that sum changes results for large-offset outcomes. */
    double variance = (moment->squares.hi + moment->squares.lo) / (moment->n - 1);
    return cpsmatch_stat_value(variance < 0.0 ? 0.0 : variance);
}

void cpsmatch_compute_statistics(size_t n, const double *treatment, const double *outcome,
                                 const double *weight, const double *support,
                                 cpsmatch_statistics *result)
{
    /* Full control/treated arms, supported treated, and used controls. The
     * control variance is unweighted; the matched control mean is weighted. */
    cpsmatch_moment moments[4] = {0};
    dd_real weighted_sum = {0}, sum_weight = {0}, sum_weight_squared = {0};
    double weighted_anchor = 0.0;
    size_t n_weighted = 0;
    memset(result, 0, sizeof(*result));

    for (size_t i = 0; i < n; i++)
    {
        double t = treatment[i];
        if (t != 0.0 && t != 1.0)
            continue;
        double s = support ? support[i] : 1.0;
        double w = weight ? weight[i] : 0.0;
        if (t == 0.0)
        {
            result->n0_on += s == 1.0;
            result->n0_off += s == 0.0;
        }
        else
        {
            result->n1_on += s == 1.0;
            result->n1_off += s == 0.0;
        }
        if (!outcome)
            continue;
        if (t == 0.0 && cpsmatch_stat_nonmissing(w))
        {
            /* gen double weight_squared = _weight^2 rounds each product to
             * double and turns overflow into missing before summarize. */
            double square = w * w;
            if (cpsmatch_stat_nonmissing(square))
                sum_weight_squared = dd_add_d(sum_weight_squared, square);
        }
        double y = outcome[i];
        if (!cpsmatch_stat_nonmissing(y))
            continue;
        cpsmatch_moment_add(&moments[(int)t], y);
        if (t == 1.0 && s == 1.0)
            cpsmatch_moment_add(&moments[2], y);
        if (t == 0.0 && w > 0.0 && cpsmatch_stat_nonmissing(w))
        {
            cpsmatch_moment_add(&moments[3], y);
            if (s == 1.0)
            {
                if (n_weighted == 0)
                    weighted_anchor = y;
                weighted_sum = dd_add_d(weighted_sum, w * (y - weighted_anchor));
                sum_weight = dd_add_d(sum_weight, w);
                n_weighted++;
            }
        }
    }
    for (size_t j = 0; j < 4; j++)
    {
        cpsmatch_moment *moment = &moments[j];
        moment->mean =
            moment->n
                ? cpsmatch_stat_value(moment->anchor +
                                      (moment->shifted_sum.hi + moment->shifted_sum.lo) / moment->n)
                : SV_missval;
    }
    if (outcome)
    {
        for (size_t i = 0; i < n; i++)
        {
            double t = treatment[i], y = outcome[i];
            if ((t != 0.0 && t != 1.0) || !cpsmatch_stat_nonmissing(y))
                continue;
            double s = support ? support[i] : 1.0;
            double w = weight ? weight[i] : 0.0;
            cpsmatch_moment_deviation(&moments[(int)t], y);
            if (t == 1.0 && s == 1.0)
                cpsmatch_moment_deviation(&moments[2], y);
            if (t == 0.0 && w > 0.0 && cpsmatch_stat_nonmissing(w))
                cpsmatch_moment_deviation(&moments[3], y);
        }
    }

    result->y1 = moments[2].mean;
    result->y0 = n_weighted
                     ? cpsmatch_stat_value(weighted_anchor + (weighted_sum.hi + weighted_sum.lo) /
                                                                 (sum_weight.hi + sum_weight.lo))
                     : SV_missval;
    result->n1 = (double)moments[2].n;
    result->n0 = (double)n_weighted;
    result->var1 = cpsmatch_moment_variance(&moments[2]);
    result->var0 = cpsmatch_moment_variance(&moments[3]);
    result->weight_squared = cpsmatch_stat_value(sum_weight_squared.hi + sum_weight_squared.lo);
    result->y1_u = moments[1].mean;
    result->y0_u = moments[0].mean;
    result->n1_u = (double)moments[1].n;
    result->n0_u = (double)moments[0].n;
    result->var1_u = cpsmatch_moment_variance(&moments[1]);
    result->var0_u = cpsmatch_moment_variance(&moments[0]);
}

ST_retcode cpsmatch_save_statistics(const cpsmatch_statistics *result)
{
    const char *names[] = {"_cpsmatch_y1",     "_cpsmatch_y0",    "_cpsmatch_n1",
                           "_cpsmatch_n0",     "_cpsmatch_var1",  "_cpsmatch_var0",
                           "_cpsmatch_w2",     "_cpsmatch_y1_u",  "_cpsmatch_y0_u",
                           "_cpsmatch_n1_u",   "_cpsmatch_n0_u",  "_cpsmatch_var1_u",
                           "_cpsmatch_var0_u", "_cpsmatch_n0_on", "_cpsmatch_n0_off",
                           "_cpsmatch_n1_on",  "_cpsmatch_n1_off"};
    const double values[] = {result->y1,
                             result->y0,
                             result->n1,
                             result->n0,
                             result->var1,
                             result->var0,
                             result->weight_squared,
                             result->y1_u,
                             result->y0_u,
                             result->n1_u,
                             result->n0_u,
                             result->var1_u,
                             result->var0_u,
                             result->n0_on,
                             result->n0_off,
                             result->n1_on,
                             result->n1_off};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++)
    {
        ST_retcode rc = SF_scal_save((char *)names[i], values[i]);
        if (rc)
            return rc;
    }
    return 0;
}
