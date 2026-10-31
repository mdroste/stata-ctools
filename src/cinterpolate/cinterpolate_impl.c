/* Linear or local cubic interpolation after one stable index sort. */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "cinterpolate_impl.h"
#include "ctools_types.h"
#include "ctools_order.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"

static double line(double x, double x0, double y0, double x1, double y1)
{
    double m = (y1 - y0) / (x1 - x0);
    double b = y0 - m * x0;
    double v = m * x + b;
    return isfinite(v) && v < SV_missval && v > -SV_missval ? v : SV_missval;
}

/* Four distinct anchors in sorted-position space. Ratio products avoid cubing
 * large x differences; extended precision reduces cancellation where available.
 * This is local Lagrange interpolation, not a cubic spline. */
static double cubic_value(double target, const double *x, const double *out, const perm_idx_t *p,
                          const size_t anchors[4])
{
    long double value = 0;
    for (int i = 0; i < 4; i++)
    {
        long double term = out[p[anchors[i]]];
        for (int j = 0; j < 4; j++)
        {
            if (i != j)
                term *= ((long double)target - x[p[anchors[j]]]) /
                        ((long double)x[p[anchors[i]]] - x[p[anchors[j]]]);
        }
        value += term;
    }
    return isfinite(value) && value < SV_missval && value > -SV_missval ? (double)value
                                                                        : SV_missval;
}

ST_retcode cinterpolate_main(const char *args)
{
    int nby, epolate, verbose;
    if (ctools_parse_next_int(&args, &nby) || nby < 0 || ctools_parse_next_int(&args, &epolate) ||
        ctools_parse_next_int(&args, &verbose) || SF_nvars() != nby + 3)
        return 198;
    int cubic = ctools_parse_bool_option(args, "cubic");
    if (cubic && epolate)
        return 198;
    double t0 = ctools_timer_seconds();
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    int rc = 0;
    int *vars = malloc((nby + 2) * sizeof(*vars));
    int *keys = malloc((nby + 1) * sizeof(*keys));
    size_t *groups = NULL;
    double *out = NULL;
    size_t nmissing = (size_t)SF_nobs();
    if (!vars || !keys)
    {
        rc = 920;
        goto done;
    }
    for (int k = 0; k < nby + 2; k++)
        vars[k] = k + 1;
    for (int k = 0; k < nby; k++)
        keys[k] = k + 2;
    keys[nby] = 1;
    int flags =
        ctools_parse_bool_option(args, "skipif") ? CTOOLS_LOAD_SKIP_IF : CTOOLS_LOAD_CHECK_IF;
    if ((rc = ctools_stata_rc(ctools_data_load(&fd, vars, nby + 2, 0, 0, flags))))
    {
        goto done;
    }
    size_t n = fd.data.nobs;
    if (!n)
        goto done;
    if (fd.data.vars[0].type != STATA_TYPE_DOUBLE || fd.data.vars[1].type != STATA_TYPE_DOUBLE)
    {
        rc = 109;
        goto done;
    }
    double tload = ctools_timer_seconds();
    if ((rc = ctools_stata_rc(ctools_order_stable(&fd.data, keys, nby + 1))))
    {
        goto done;
    }
    double tsort = ctools_timer_seconds();
    groups = malloc((n + 1) * sizeof(*groups));
    out = malloc(n * sizeof(*out));
    if (!groups || !out)
    {
        rc = 920;
        goto done;
    }
    const double *x = fd.data.vars[1].data.dbl, *y = fd.data.vars[0].data.dbl;
    const perm_idx_t *p = fd.data.sort_order;
/* Group discovery and initialization are bandwidth-sized parallel scans.
 * Compact flags in place afterward; writes never overtake unread flags. */
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if (n > 10000)
    for (size_t i = 0; i < n; i++)
    {
        out[i] = SV_missval;
        groups[i] =
            i && nby && ctools_compare_rows(&fd.data, p[i - 1], keys, &fd.data, p[i], keys, nby);
    }
    size_t ng = 1;
    groups[0] = 0;
    if (nby)
        for (size_t i = 1; i < n; i++)
            if (groups[i])
                groups[ng++] = i;
    groups[ng] = n;
#pragma omp parallel for schedule(dynamic, 16)                                                     \
    num_threads(ctools_get_max_threads()) if (n > 10000 && ng > 1)
    for (size_t g = 0; g < ng; g++)
    {
        size_t begin = groups[g], end = groups[g + 1];
        size_t first = end, last = end, second = end;
        size_t anchors[4], nanchors = 0;
        for (size_t a = begin; a < end;)
        {
            size_t b = a + 1, count = 0;
            long double sum = 0;
            while (b < end && x[p[b]] == x[p[a]])
                b++;
            if (x[p[a]] < SV_missval)
            {
                for (size_t j = a; j < b; j++)
                    if (y[p[j]] < SV_missval)
                    {
                        sum += y[p[j]];
                        count++;
                    }
                if (count)
                {
                    double mean = (double)(sum / count);
                    for (size_t j = a; j < b; j++)
                        out[p[j]] = mean;
                    if (cubic)
                    {
                        /* Fill only the middle interval once its second right
                         * anchor is known. Each gap is visited at most once. */
                        if (nanchors == 4)
                        {
                            for (int k = 0; k < 3; k++)
                                anchors[k] = anchors[k + 1];
                            nanchors = 3;
                        }
                        anchors[nanchors++] = a;
                        if (nanchors == 4)
                            for (size_t j = anchors[1] + 1; j < anchors[2]; j++)
                                if (out[p[j]] >= SV_missval)
                                    out[p[j]] = cubic_value(x[p[j]], x, out, p, anchors);
                        a = b;
                        continue;
                    }
                    if (first == end)
                        first = a;
                    else
                    {
                        if (second == end)
                            second = a;
                        for (size_t j = last + 1; j < a; j++)
                            if (out[p[j]] >= SV_missval)
                                out[p[j]] = line(x[p[j]], x[p[last]], out[p[last]], x[p[a]], mean);
                    }
                    last = a;
                }
            }
            a = b;
        }
        if (epolate && second != end)
        {
            /* Native ipolate obtains endpoint slopes from adjacent rows after
             * interpolation. Reusing original anchors changes rounding when
             * the neighboring row was interpolated. Skip duplicate x values. */
            size_t left = end, right = end;
            for (size_t j = first; j + 1 < end; j++)
            {
                if (x[p[j + 1]] >= SV_missval || x[p[j]] == x[p[j + 1]] ||
                    out[p[j]] >= SV_missval || out[p[j + 1]] >= SV_missval)
                    continue;
                if (left == end)
                    left = j;
                right = j;
            }
            if (left != end)
            {
                for (size_t j = begin; j < first; j++)
                    out[p[j]] =
                        line(x[p[j]], x[p[left]], out[p[left]], x[p[left + 1]], out[p[left + 1]]);
                for (size_t j = last + 1; j < end; j++)
                    if (x[p[j]] < SV_missval && out[p[j]] >= SV_missval)
                        out[p[j]] = line(x[p[j]], x[p[right]], out[p[right]], x[p[right + 1]],
                                         out[p[right + 1]]);
            }
        }
    }
    double tcompute = ctools_timer_seconds();
    /* Include excluded observations: their initialized destination stays missing. */
    for (size_t i = 0; i < n; i++)
        if (out[i] < SV_missval)
            nmissing--;
    if ((rc = ctools_stata_rc(ctools_store_filtered(out, n, nby + 3, fd.obs_map))))
    {
        goto done;
    }
    ctools_verbose("cinterpolate", verbose,
                   "load %.4fs; sort %.4fs; interpolate %.4fs; store %.4fs", tload - t0,
                   tsort - tload, tcompute - tsort, ctools_timer_seconds() - tcompute);
done:
    if (!rc)
    {
        char count[32];
        snprintf(count, sizeof(count), "%zu", nmissing);
        rc = SF_macro_save("_cinterpolate_nmissing", count);
    }
    free(vars);
    free(keys);
    free(groups);
    free(out);
    ctools_filtered_data_free(&fd);
    return rc;
}
