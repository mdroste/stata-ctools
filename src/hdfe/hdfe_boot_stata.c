/*
 * hdfe_boot_stata.c
 *
 * The Stata side of the pairs bootstrap shared by the HDFE estimators (see
 * hdfe_boot_stata.h).
 * Part of the ctools Stata plugin suite
 */

#define SD_SAFEMODE
#include "hdfe_boot_stata.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../ctools_spi.h"

/* Value q (0-based) of those the ado packed into the call's variables
 * (_ctools_hdfe_boot.ado): variable q / nobs + 1, observation q % nobs + 1 */
static ST_retcode packed_value(ST_int nobs, size_t q, ST_double *x)
{
    return SF_vdata((ST_int)(q / (size_t)nobs) + 1, (ST_int)(q % (size_t)nobs) + 1, x);
}

static void draws_error(const char *cmd, const char *what)
{
    char msg[256];
    snprintf(msg, sizeof(msg), "%s: %s\n", cmd, what);
    SF_error(msg);
}

static ST_retcode read_uniforms(const hdfe_boot_design *d, int cluster_mode, ST_int count,
                                const char *cmd, ST_int *mult)
{
    const ST_int Ud = hdfe_boot_draw_units(d);
    if (Ud < 1)
    {
        draws_error(cmd, "draws(stata) draws at most 2^31 - 1 observations (the fweights stand "
                         "for more); use draws(fast)");
        return 198;
    }
    const ST_int nobs = SF_nobs();
    const size_t n = 2 * (size_t)count * (size_t)Ud;
    if (nobs < 1 || Ud < 1 || (size_t)SF_nvars() != (n + (size_t)nobs - 1) / (size_t)nobs)
    {
        draws_error(cmd, "boot_chunk got the wrong number of draw variables");
        return 198;
    }
    ST_double *ud = malloc((size_t)Ud * sizeof(ST_double));
    ST_double *ut = malloc((size_t)Ud * sizeof(ST_double));
    ST_retcode rc = (!ud || !ut) ? 920 : 0;
    for (ST_int r = 0; r < count && !rc; r++)
    {
        const size_t base = 2 * (size_t)r * (size_t)Ud;
        for (ST_int u = 0; u < Ud && !rc; u++)
        {
            rc = packed_value(nobs, base + (size_t)u, &ud[u]);
            if (!rc)
                rc = packed_value(nobs, base + (size_t)Ud + (size_t)u, &ut[u]);
        }
        if (!rc)
            rc = hdfe_boot_map_stata(d, ud, ut, cluster_mode, mult + (size_t)r * (size_t)d->U);
    }
    free(ud);
    free(ut);
    return rc;
}

static ST_retcode read_supplied(const hdfe_boot_design *d, ST_int count, ST_int N,
                                const ST_int *unit_of_row, const ST_int *obs_of_row,
                                const char *cmd, ST_int *mult)
{
    if (SF_nvars() != count)
    {
        draws_error(cmd, "boot_chunk got the wrong number of _bsdraws() variables");
        return 198;
    }
    const ST_int U = d->U;
    for (ST_int r = 0; r < count; r++)
    {
        ST_int *m = mult + (size_t)r * (size_t)U;
        for (ST_int u = 0; u < U; u++)
            m[u] = -1;
        for (ST_int i = 0; i < N; i++)
        {
            ST_double x = 0.0;
            ST_retcode rc = SF_vdata(r + 1, obs_of_row[i], &x);
            if (rc)
                return rc;
            if (SF_is_missing(x) || !(x >= 0.0) || x != floor(x) || x > 2147483647.0)
            {
                draws_error(cmd, "_bsdraws() takes nonnegative integer multiplicities");
                return 198;
            }
            ST_int u = unit_of_row[i];
            if (m[u] >= 0 && m[u] != (ST_int)x)
            {
                draws_error(cmd, "_bsdraws(): the rows of a resampling unit have different "
                                 "multiplicities");
                return 198;
            }
            m[u] = (ST_int)x;
        }
        for (ST_int u = 0; u < U; u++)
        {
            if (m[u] < 0)
                m[u] = 0;
        }
    }
    return 0;
}

ST_retcode hdfe_boot_read_draws(const hdfe_boot_design *d, int cluster_mode, int supplied,
                                ST_int count, ST_int N, const ST_int *unit_of_row,
                                const ST_int *obs_of_row, const char *cmd, ST_int **out)
{
    *out = NULL;
    if (count < 1 || d->U < 1)
        return 198;
    ST_int *mult = malloc((size_t)count * (size_t)d->U * sizeof(ST_int));
    if (!mult)
        return 920;
    ST_retcode rc = supplied
                        ? read_supplied(d, count, N, unit_of_row, obs_of_row, cmd, mult)
                        : read_uniforms(d, cluster_mode, count, cmd, mult);
    if (rc)
    {
        free(mult);
        return rc;
    }
    *out = mult;
    return 0;
}

ST_retcode hdfe_boot_save_dup_strata(const hdfe_boot_design *d, const char *cmd)
{
    double off[HDFE_BOOT_DUP_MAX], size[HDFE_BOOT_DUP_MAX];
    ST_int K = 0;
    double pos = 0.0;
    for (ST_int s = 0; s < d->S; s++)
    {
        const ST_int lo = d->stratum_start ? d->stratum_start[s] : 0;
        const ST_int hi = d->stratum_start ? d->stratum_start[s + 1] : d->U;
        double n = 0.0;
        if (d->unit_freq)
        {
            for (ST_int u = lo; u < hi; u++)
                n += d->unit_freq[u];
        }
        else
            n = (double)(hi - lo);
        if (n >= 2.0)
        {
            if (K < HDFE_BOOT_DUP_MAX)
            {
                off[K] = pos;
                size[K] = n;
            }
            K++;
        }
        pos += n;
    }
    char name[64];
    snprintf(name, sizeof(name), "__%s_bs_dupK", cmd);
    ST_retcode rc = ctools_scal_save(name, K <= HDFE_BOOT_DUP_MAX ? (ST_double)K : SV_missval);
    for (ST_int j = 0; j < K && K <= HDFE_BOOT_DUP_MAX && !rc; j++)
    {
        snprintf(name, sizeof(name), "__%s_bs_dupo%d", cmd, (int)(j + 1));
        rc = ctools_scal_save(name, off[j]);
        if (!rc)
        {
            snprintf(name, sizeof(name), "__%s_bs_dupn%d", cmd, (int)(j + 1));
            rc = ctools_scal_save(name, size[j]);
        }
    }
    return rc;
}

ST_retcode hdfe_boot_store_rows(const hdfe_boot_acc *acc, ST_int first, ST_int count,
                                const char *name)
{
    ST_retcode rc = 0;
    for (ST_int r = 0; r < count && !rc; r++)
    {
        const ST_double *row = acc->b + (size_t)(first + r) * (size_t)acc->K_b;
        for (ST_int k = 0; k < acc->K_b && !rc; k++)
            rc = ctools_mat_store(name, r + 1, k + 1, isfinite(row[k]) ? row[k] : SV_missval);
    }
    return rc;
}
