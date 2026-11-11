/*
 * cppmlhdfe_boot.c
 *
 * vce(bootstrap) for cppmlhdfe (see cppmlhdfe_boot.h).
 * Part of the ctools Stata plugin suite
 */

/* Checked SPI access: indices resolve in the plugin varlist */
#define SD_SAFEMODE
#include "cppmlhdfe_boot.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cppmlhdfe_fit.h"
#include "hdfe/hdfe_boot.h"
#include "hdfe/hdfe_boot_stata.h"
#include "hdfe/hdfe_sample.h"
#include "hdfe/hdfe_solver.h"
#include "hdfe/hdfe_wild.h"
#include "../ctools_config.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_matrix.h"
#include "../ctools_ols.h"
#include "../ctools_parse.h"
#include "../ctools_platform.h"
#include "../ctools_spi.h"

typedef struct
{
    ST_int nobs; /* SF_nobs() when retained: the data may not change between calls */
    ST_int N, K_x, G;
    ST_double *y, *X, *offset, *w_user, *eta, *stdev_x;
    ST_double stdev_y;
    ST_int weight_type;
    ST_int **levels;
    ST_int *num_levels;
    /* Heterogeneous slopes: factor g's slope variable on the compact rows
     * (NULL for an intercept, slope NULL without slopes), and whether it is
     * centered within its levels */
    ST_double **slope;
    ST_int slope_centered[10];
    /* Centered slope factor g: the intercept factor on the same levels, or
     * -1; slope_est (L_total, slope nodes only): the full fit estimates the
     * level's slope (its variable varies within the level) */
    ST_int slope_intercept[10];
    unsigned char *slope_est;
    PPML_Options options;
    ST_double tolerance, irls_tol, inner_tol_override;
    ST_int maxiter, irls_maxiter, keep_singletons, check_fe, check_relu;
    ST_int *is_collinear;
    ST_int K_keep;
    ST_double *b_hat;
    ST_int *obs_of_row;
    /* Resampling frame */
    ST_int *unit_of_row, *stratum_start;
    ST_double *unit_freq; /* fweighted rows resampled: the fweights, or NULL */
    hdfe_boot_design design;
    int cluster_mode;
    hdfe_boot_acc acc;
    /* fevce(bootstrap): levels of all effects as nodes offset[g] + level - 1,
     * and each node's full-sample connected set over all effects jointly
     * (the minimum-norm normalization is global) */
    ST_int *node_off, L_total;
    ST_int *joint_comp, n_joint;
    int fe_on, dse_on;
    /* wildbootstrap(): clusters (N dense ids, or NULL: each row), the tested
     * columns of X, their analytic standard errors, the replicates' weights */
    ST_int *wb_vce, *wb_boot, wb_G, wb_C;
    ST_int *wb_test, wb_ntest;
    ST_double *wb_scale;
    hdfe_wild_score_model wb_score;
} ppml_boot_model;

static ppml_boot_model *g_pboot = NULL;

static void wild_free(ppml_boot_model *m);
static ST_retcode save_parallel(const ppml_boot_model *m);
static ST_retcode slope_levels(ppml_boot_model *m);

static void model_free(ppml_boot_model *m)
{
    if (!m)
        return;
    wild_free(m);
    free(m->y);
    free(m->X);
    free(m->offset);
    free(m->w_user);
    free(m->eta);
    free(m->stdev_x);
    for (ST_int g = 0; m->levels && g < m->G; g++)
        free(m->levels[g]);
    free(m->levels);
    for (ST_int g = 0; m->slope && g < m->G; g++)
        free(m->slope[g]);
    free(m->slope);
    free(m->slope_est);
    free(m->num_levels);
    free(m->is_collinear);
    free(m->b_hat);
    free(m->obs_of_row);
    free(m->unit_of_row);
    free(m->stratum_start);
    free(m->unit_freq);
    hdfe_boot_acc_free(&m->acc);
    free(m->node_off);
    free(m->joint_comp);
    free(m);
}

void cppmlhdfe_boot_release(void)
{
    model_free(g_pboot);
    g_pboot = NULL;
}

void cppmlhdfe_boot_install(void)
{
    cppmlhdfe_boot_retain_hook = cppmlhdfe_boot_retain;
    cppmlhdfe_boot_release_hook = cppmlhdfe_boot_release;
}

/* ========================================================================
 * Retain
 * ======================================================================== */

typedef struct
{
    double stratum, cluster;
    ST_int row;
} unit_key;

static int unit_key_cmp(const void *a, const void *b)
{
    const unit_key *x = a, *y = b;
    if (x->stratum != y->stratum)
        return x->stratum < y->stratum ? -1 : 1;
    if (x->cluster != y->cluster)
        return x->cluster < y->cluster ? -1 : 1;
    return (x->row > y->row) - (x->row < y->row);
}

/* Units in bsample's order: strata ascending; within a stratum, clusters in
 * ascending value, or rows in data order. A cluster that spans strata is a
 * unit within each stratum, as bsample numbers clusters within strata. */
static ST_retcode build_frame(ppml_boot_model *m, ST_int unit_var, ST_int strata_var)
{
    const ST_int N = m->N;
    unit_key *keys = malloc((size_t)N * sizeof(*keys));
    m->unit_of_row = malloc((size_t)N * sizeof(ST_int));
    if (!keys || !m->unit_of_row)
    {
        free(keys);
        return 920;
    }
    for (ST_int i = 0; i < N; i++)
    {
        double s = 0.0, c = 0.0;
        ST_retcode rc = 0;
        if (strata_var > 0)
            rc = SF_vdata(strata_var, m->obs_of_row[i], &s);
        if (!rc && unit_var > 0)
            rc = SF_vdata(unit_var, m->obs_of_row[i], &c);
        if (rc)
        {
            free(keys);
            return rc;
        }
        keys[i].stratum = s;
        keys[i].cluster = c;
        keys[i].row = i;
    }
    m->cluster_mode = unit_var > 0;
    qsort(keys, (size_t)N, sizeof(*keys), unit_key_cmp);
    ST_int U = 0, S = 0;
    for (ST_int i = 0; i < N; i++)
    {
        int new_stratum = i == 0 || keys[i].stratum != keys[i - 1].stratum;
        int new_unit = !m->cluster_mode || new_stratum || keys[i].cluster != keys[i - 1].cluster;
        S += new_stratum;
        U += new_unit;
        m->unit_of_row[keys[i].row] = U - 1;
    }
    m->stratum_start = malloc((size_t)(S + 1) * sizeof(ST_int));
    if (!m->stratum_start)
    {
        free(keys);
        return 920;
    }
    S = 0;
    for (ST_int i = 0; i < N; i++)
    {
        if (i == 0 || keys[i].stratum != keys[i - 1].stratum)
            m->stratum_start[S++] = m->unit_of_row[keys[i].row];
    }
    m->stratum_start[S] = U;
    /* Frequency-weighted rows stand for fw observations each: resample those
     * (plan 5.3; observation-level resampling only) */
    if (m->w_user && m->weight_type == 2 && !m->cluster_mode)
    {
        m->unit_freq = malloc((size_t)U * sizeof(ST_double));
        if (!m->unit_freq)
        {
            free(keys);
            return 920;
        }
        for (ST_int i = 0; i < N; i++)
            m->unit_freq[m->unit_of_row[i]] = m->w_user[i];
    }
    free(keys);
    m->design.U = U;
    m->design.S = S;
    m->design.stratum_start = m->stratum_start;
    m->design.unit_freq = m->unit_freq;
    m->design.seed = 0;
    return hdfe_boot_design_check(&m->design);
}

/* Connected sets of the effects' levels over rows rows[0..n) (all rows when
 * rows is NULL): comp[node] numbers each node's set densely */
static ST_retcode joint_sets(const ppml_boot_model *m, const ST_int *rows, ST_int n,
                             ST_int *comp, ST_int *n_comp)
{
    ctools_UnionFind *uf = ctools_uf_create(m->L_total);
    ST_int *dense = malloc((size_t)m->L_total * sizeof(ST_int));
    if (!uf || !dense)
    {
        if (uf)
            ctools_uf_destroy(uf);
        free(dense);
        return 920;
    }
    /* The intercepts only: a slope is identified within its level, and its
     * nodes stay sets of their own */
    ST_int g0 = 0;
    while (g0 < m->G && m->slope && m->slope[g0])
        g0++;
    for (ST_int r = 0; r < n && g0 < m->G; r++)
    {
        const ST_int i = rows ? rows[r] : r;
        const ST_int a = m->node_off[g0] + m->levels[g0][i] - 1;
        for (ST_int g = g0 + 1; g < m->G; g++)
        {
            if (m->slope && m->slope[g])
                continue;
            ctools_uf_union(uf, a, m->node_off[g] + m->levels[g][i] - 1);
        }
    }
    for (ST_int v = 0; v < m->L_total; v++)
        dense[v] = -1;
    ST_int k = 0;
    for (ST_int v = 0; v < m->L_total; v++)
    {
        const ST_int root = ctools_uf_find(uf, v);
        if (dense[root] < 0)
            dense[root] = k++;
        comp[v] = dense[root];
    }
    *n_comp = k;
    free(dense);
    ctools_uf_destroy(uf);
    return 0;
}

static ST_double *copy_doubles(const ST_double *src, size_t n)
{
    ST_double *dst = malloc((n ? n : 1) * sizeof(*dst));
    if (dst && n)
        memcpy(dst, src, n * sizeof(*dst));
    return dst;
}

ST_retcode cppmlhdfe_boot_retain(const cppmlhdfe_boot_source *src)
{
    cppmlhdfe_boot_release();
    ppml_boot_model *m = calloc(1, sizeof(*m));
    if (!m)
        return 920;
    const ST_int N = src->N, K_x = src->K_x, G = src->G;
    m->nobs = SF_nobs();
    m->N = N;
    m->K_x = K_x;
    m->G = G;
    m->stdev_y = src->stdev_y;
    m->weight_type = src->weight_type;
    m->options = *src->options;
    m->tolerance = src->tolerance;
    m->irls_tol = src->irls_tol;
    m->inner_tol_override = src->inner_tol_override;
    m->maxiter = src->maxiter;
    m->irls_maxiter = src->irls_maxiter;
    m->keep_singletons = src->keep_singletons;
    m->check_fe = src->check_fe;
    m->check_relu = src->check_relu;
    m->K_keep = src->K_keep;
    m->y = copy_doubles(src->y, (size_t)N);
    m->X = copy_doubles(src->X, (size_t)N * (size_t)K_x);
    m->offset = src->offset ? copy_doubles(src->offset, (size_t)N) : NULL;
    m->w_user = src->w_user ? copy_doubles(src->w_user, (size_t)N) : NULL;
    m->eta = copy_doubles(src->eta, (size_t)N);
    m->stdev_x = copy_doubles(src->stdev_x, (size_t)K_x);
    m->b_hat = copy_doubles(src->b_hat, (size_t)src->K_keep + 1);
    m->is_collinear = malloc((size_t)K_x * sizeof(ST_int));
    m->levels = calloc((size_t)G, sizeof(ST_int *));
    m->num_levels = malloc((size_t)G * sizeof(ST_int));
    m->obs_of_row = malloc((size_t)N * sizeof(ST_int));
    ST_retcode rc = 0;
    if (!m->y || !m->X || (src->offset && !m->offset) || (src->w_user && !m->w_user) ||
        !m->eta || !m->stdev_x || !m->b_hat || !m->is_collinear || !m->levels ||
        !m->num_levels || !m->obs_of_row)
        rc = 920;
    if (!rc)
    {
        memcpy(m->is_collinear, src->is_collinear, (size_t)K_x * sizeof(ST_int));
        for (ST_int i = 0; i < N; i++)
            m->obs_of_row[i] = (ST_int)src->obs_map[i];
    }
    for (ST_int g = 0; g < G && !rc; g++)
    {
        m->num_levels[g] = src->num_levels[g];
        if (!(m->levels[g] = malloc((size_t)N * sizeof(ST_int))))
            rc = 920;
        else
            memcpy(m->levels[g], src->levels[g], (size_t)N * sizeof(ST_int));
    }
    for (ST_int g = 0; g < G && !rc && src->slope; g++)
    {
        if (!src->slope[g])
            continue;
        if (!m->slope && !(m->slope = calloc((size_t)G, sizeof(ST_double *))))
            rc = 920;
        else if (!(m->slope[g] = copy_doubles(src->slope[g], (size_t)N)))
            rc = 920;
        else
            m->slope_centered[g] = src->slope_centered ? src->slope_centered[g] : 0;
    }
    if (!rc)
        rc = build_frame(m, (ST_int)ppml_option("bs_unit_idx", 0),
                         (ST_int)ppml_option("bs_strata_idx", 0));
    if (!rc && !(m->node_off = malloc(((size_t)G + 1) * sizeof(ST_int))))
        rc = 920;
    if (!rc)
    {
        m->node_off[0] = 0;
        for (ST_int g = 0; g < G; g++)
            m->node_off[g + 1] = m->node_off[g] + m->num_levels[g];
        m->L_total = m->node_off[G];
        if (G >= 2)
        {
            if (!(m->joint_comp = malloc((size_t)m->L_total * sizeof(ST_int))))
                rc = 920;
            else
                rc = joint_sets(m, NULL, N, m->joint_comp, &m->n_joint);
        }
    }
    if (!rc)
        rc = slope_levels(m);
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_bs_U", hdfe_boot_draw_total(&m->design));
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_bs_S", (ST_double)m->design.S);
    if (!rc)
        rc = hdfe_boot_save_dup_strata(&m->design, "cppmlhdfe");
    if (!rc)
        rc = save_parallel(m);
    if (rc)
    {
        model_free(m);
        return rc;
    }
    g_pboot = m;
    return 0;
}

/* ========================================================================
 * One replicate
 * ======================================================================== */

/* Sample screening on the drawn rows (alive[r]), as ppml_select_fe_sample_ex
 * screens the full sample with each row counted mult times: a level whose
 * rows total one draw is a singleton, a level without a positive outcome is
 * separated; repeated until nothing changes. Returns the rows left. */
static ST_int screen(const ppml_boot_model *m, const ST_int *rows, const ST_int *mult,
                     ST_int N1, unsigned char *alive, ST_int **count, unsigned char **positive)
{
    const int check = m->check_fe && !m->keep_singletons;
    for (;;)
    {
        for (ST_int g = 0; g < m->G; g++)
        {
            memset(count[g], 0, (size_t)m->num_levels[g] * sizeof(ST_int));
            memset(positive[g], 0, (size_t)m->num_levels[g]);
            for (ST_int r = 0; r < N1; r++)
            {
                if (!alive[r])
                    continue;
                ST_int l = m->levels[g][rows[r]] - 1;
                count[g][l] += mult[r];
                if (m->y[rows[r]] > 0)
                    positive[g][l] = 1;
            }
        }
        ST_int removed = 0;
        for (ST_int r = 0; r < N1; r++)
        {
            if (!alive[r])
                continue;
            int drop = 0;
            for (ST_int g = 0; g < m->G && !drop; g++)
            {
                ST_int l = m->levels[g][rows[r]] - 1;
                drop = (!m->keep_singletons && count[g][l] == 1) ||
                       (check && !(m->slope && m->slope[g]) && !positive[g][l]);
            }
            if (drop)
            {
                alive[r] = 0;
                removed++;
            }
        }
        if (!removed)
            break;
    }
    ST_int n = 0;
    for (ST_int r = 0; r < N1; r++)
        n += alive[r];
    return n;
}

/* flat[l] for slope factor g over rows idx[0 .. n) of the model (idx NULL:
 * rows 0 .. n): 1 where the slope variable does not vary within level l
 * (centered slope) or is zero there (uncentered), 0 where it does, -1 where
 * the level has none of the rows. */
static ST_retcode slope_flat(const ppml_boot_model *m, ST_int g, const ST_int *idx, ST_int n,
                             signed char *flat)
{
    const ST_int L = m->num_levels[g];
    ST_double *lo = malloc((size_t)(L > 0 ? L : 1) * sizeof(ST_double));
    ST_double *hi = malloc((size_t)(L > 0 ? L : 1) * sizeof(ST_double));
    if (!lo || !hi)
    {
        free(lo);
        free(hi);
        return 920;
    }
    for (ST_int l = 0; l < L; l++)
    {
        lo[l] = INFINITY;
        hi[l] = -INFINITY;
    }
    for (ST_int k = 0; k < n; k++)
    {
        const ST_int i = idx ? idx[k] : k, l = m->levels[g][i] - 1;
        lo[l] = fmin(lo[l], m->slope[g][i]);
        hi[l] = fmax(hi[l], m->slope[g][i]);
    }
    for (ST_int l = 0; l < L; l++)
    {
        const ST_double scale = fmax(fabs(lo[l]), fabs(hi[l]));
        const int is_flat = m->slope_centered[g] ? !(hi[l] - lo[l] > 1e-12 * scale)
                                                 : !(scale > 0.0);
        flat[l] = !(hi[l] >= lo[l]) ? -1 : (signed char)is_flat;
    }
    free(lo);
    free(hi);
    return 0;
}

/* For fevce(bootstrap) with slopes: each centered slope's intercept factor
 * (the one on the same levels, as ppml_fe_minnorm pairs them) and the
 * levels whose slope the full fit estimates */
static ST_retcode slope_levels(ppml_boot_model *m)
{
    for (ST_int g = 0; g < m->G; g++)
        m->slope_intercept[g] = -1;
    if (!m->slope)
        return 0;
    if (!(m->slope_est = calloc((size_t)(m->L_total > 0 ? m->L_total : 1), 1)))
        return 920;
    for (ST_int g = 0; g < m->G; g++)
    {
        if (!m->slope[g])
            continue;
        for (ST_int h = 0; m->slope_centered[g] && h < m->G; h++)
        {
            if (!m->slope[h] && m->num_levels[h] == m->num_levels[g] &&
                !memcmp(m->levels[h], m->levels[g], (size_t)m->N * sizeof(ST_int)))
            {
                m->slope_intercept[g] = h;
                break;
            }
        }
        signed char *flat = malloc((size_t)(m->num_levels[g] > 0 ? m->num_levels[g] : 1));
        if (!flat)
            return 920;
        ST_retcode rc = slope_flat(m, g, NULL, m->N, flat);
        for (ST_int l = 0; !rc && l < m->num_levels[g]; l++)
            m->slope_est[m->node_off[g] + l] = flat[l] == 0;
        free(flat);
        if (rc)
            return rc;
    }
    return 0;
}

/* fevce(bootstrap): the replicate's fixed effects as the full fit saves
 * them. d = z - resid - Xb - mean(d) on the replicate's rows (IRLS
 * weights), its minimum-norm coefficients on the effects, mapped to the
 * full-sample levels. A level is estimable when drawn and its full-sample
 * connected set did not split in the replicate (plan 8.2); the others stay
 * NaN. d is the sum of the effects on every full-sample row whose levels
 * are all estimable. A recovery that does not converge leaves them NaN. */
static ST_retcode replicate_effects(const ppml_boot_model *m, HDFE_State *S,
                                    const ST_int *rows, const ST_int *row2, ST_int N2,
                                    const ppml_irls_io *io, hdfe_boot_rep *out)
{
    const ST_int G = m->G, K_keep = io->K_keep;
    const ST_double *aug = io->aug_copy, *z = io->z_orig, *X = io->X, *w = io->irls_w;
    ST_retcode rc = 0;
    ST_double *cert = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *fe = malloc((size_t)m->L_total * sizeof(ST_double));
    unsigned char *present = calloc((size_t)m->L_total, 1);
    ST_double *alpha[10] = {0};
    ST_int *kept = NULL, *comp = NULL, *first = NULL;
    unsigned char *split = NULL;
    if (!cert || !fe || !present)
    {
        rc = 920;
        goto done;
    }
    {
        dd_real sw = {0.0, 0.0}, sd = {0.0, 0.0};
        for (ST_int k = 0; k < N2; k++)
        {
            ST_double resid = aug[k];
            for (ST_int c = 0; c < K_keep; c++)
                resid -= aug[(size_t)(io->keep_idx[c] + 1) * N2 + k] * io->beta_final[c];
            ST_double d = z[k] - resid;
            for (ST_int c = 0; c < K_keep; c++)
                d -= X[(size_t)io->keep_idx[c] * N2 + k] * io->beta_final[c];
            sw = dd_add_d(sw, w[k]);
            sd = dd_add_d(sd, w[k] * d);
            cert[k] = resid;
        }
        const ST_double mean_d = (sd.hi + sd.lo) / (sw.hi + sw.lo);
        for (ST_int k = 0; k < N2; k++)
        {
            ST_double d = z[k] - cert[k] - mean_d;
            for (ST_int c = 0; c < K_keep; c++)
                d -= X[(size_t)io->keep_idx[c] * N2 + k] * io->beta_final[c];
            cert[k] = d;
        }
    }
    for (ST_int g = 0; g < G; g++)
    {
        if (!(alpha[g] = calloc((size_t)S->factors[g].num_levels, sizeof(ST_double))))
        {
            rc = 920;
            goto done;
        }
    }
    rc = ppml_fe_minnorm(S, cert, alpha, 0);
    if (rc == 430)
    {
        rc = 0;
        goto done;
    }
    if (rc)
        goto done;
    for (ST_int v = 0; v < m->L_total; v++)
        fe[v] = NAN;
    for (ST_int k = 0; k < N2; k++)
    {
        const ST_int i = rows[row2[k]];
        for (ST_int g = 0; g < G; g++)
        {
            const ST_int node = m->node_off[g] + m->levels[g][i] - 1;
            fe[node] = alpha[g][S->factors[g].levels[k] - 1];
            present[node] = 1;
        }
    }
    if (!(kept = malloc((size_t)(N2 > 0 ? N2 : 1) * sizeof(ST_int))))
    {
        rc = 920;
        goto done;
    }
    for (ST_int k = 0; k < N2; k++)
        kept[k] = rows[row2[k]];
    /* A slope is not estimable in a level where the replicate's slope
     * variable does not vary (centered) or is zero (uncentered). When the
     * full fit estimates that slope, the level's intercept is not estimable
     * either: the replicate identifies only intercept + slope x value, and
     * the minimum-norm solution sets the slope to 0. */
    for (ST_int g = 0; m->slope && g < G; g++)
    {
        if (!m->slope[g])
            continue;
        signed char *flat = malloc((size_t)(m->num_levels[g] > 0 ? m->num_levels[g] : 1));
        if (!flat)
        {
            rc = 920;
            goto done;
        }
        rc = slope_flat(m, g, kept, N2, flat);
        const ST_int h = m->slope_intercept[g];
        for (ST_int l = 0; !rc && l < m->num_levels[g]; l++)
        {
            if (flat[l] != 1)
                continue;
            fe[m->node_off[g] + l] = NAN;
            if (h >= 0 && m->slope_est && m->slope_est[m->node_off[g] + l])
                fe[m->node_off[h] + l] = NAN;
        }
        free(flat);
        if (rc)
            goto done;
    }
    if (G >= 2)
    {
        ST_int n_comp = 0;
        comp = malloc((size_t)m->L_total * sizeof(ST_int));
        first = malloc((size_t)m->n_joint * sizeof(ST_int));
        split = calloc((size_t)m->n_joint, 1);
        if (!comp || !first || !split)
        {
            rc = 920;
            goto done;
        }
        if ((rc = joint_sets(m, kept, N2, comp, &n_comp)))
            goto done;
        for (ST_int j = 0; j < m->n_joint; j++)
            first[j] = -1;
        for (ST_int v = 0; v < m->L_total; v++)
        {
            if (!present[v])
                continue;
            const ST_int j = m->joint_comp[v];
            if (first[j] < 0)
                first[j] = comp[v];
            else if (first[j] != comp[v])
                split[j] = 1;
        }
        for (ST_int g = 0; g < G; g++)
        {
            if (m->slope && m->slope[g])
                continue;
            for (ST_int v = m->node_off[g]; v < m->node_off[g + 1]; v++)
            {
                if (split[m->joint_comp[v]])
                {
                    fe[v] = NAN;
                    out->fe_split = 1;
                }
            }
        }
    }
    if (out->fe)
        memcpy(out->fe, fe, (size_t)m->L_total * sizeof(ST_double));
    if (out->d)
    {
        for (ST_int i = 0; i < m->N; i++)
        {
            ST_double sum = 0.0;
            for (ST_int g = 0; g < G; g++)
            {
                const ST_double a = fe[m->node_off[g] + m->levels[g][i] - 1];
                sum += (m->slope && m->slope[g]) ? a * m->slope[g][i] : a;
            }
            out->d[i] = sum;
        }
    }

done:
    for (ST_int g = 0; g < G; g++)
        free(alpha[g]);
    free(cert);
    free(fe);
    free(present);
    free(kept);
    free(comp);
    free(first);
    free(split);
    return rc;
}

/* Slope factors of S, built on n compact rows (the k-th is row
 * rows[row2[k]] of the model, or k when rows is NULL): as the fit sets
 * them up, centered within the level when an intercept on the same levels
 * is absorbed (ppml_update_fe_weights computes the centers). */
static ST_retcode slope_factors(const ppml_boot_model *m, HDFE_State *S, const ST_int *rows,
                                const ST_int *row2, ST_int n)
{
    for (ST_int g = 0; m->slope && g < m->G; g++)
    {
        if (!m->slope[g])
            continue;
        FE_Factor *f = &S->factors[g];
        if (!(f->slope = malloc((size_t)(n > 0 ? n : 1) * sizeof(ST_double))))
            return 920;
        for (ST_int k = 0; k < n; k++)
            f->slope[k] = m->slope[g][rows ? rows[row2[k]] : k];
        f->has_intercept = 0;
        f->num_slopes = 1;
        if (m->slope_centered[g] &&
            !(f->slope_center = calloc((size_t)f->num_levels, sizeof(ST_double))))
            return 920;
    }
    return 0;
}

/* The fit on the rows alive: 0 with out->b and out->status set, *dropped
 * the rows a separation check removed (the caller screens again and refits),
 * or an error that aborts the bootstrap. */
static ST_retcode fit_rows(const ppml_boot_model *m, const ST_int *rows, const ST_int *mult,
                           ST_int N1, unsigned char *alive, ST_int N2, hdfe_boot_rep *out,
                           ST_int *dropped)
{
    const ST_int K_x = m->K_x, K_aug = K_x + 1, G = m->G;
    /* Column parallelism with the threads the engine left over; each
     * column is projected by one thread, so the result does not change */
    const ST_int threads = out->threads > 1 ? out->threads : 1;
    PPML_Options options = m->options;
    ST_retcode rc = 0;
    *dropped = 0;
    HDFE_State S;
    memset(&S, 0, sizeof(S));
    ppml_irls_io io;
    memset(&io, 0, sizeof(io));

    ST_int *mask = malloc((size_t)N1 * sizeof(ST_int));
    ST_int **lv = calloc((size_t)G, sizeof(ST_int *));
    ST_int *level_max = malloc((size_t)G * sizeof(ST_int));
    ST_int *row2 = malloc((size_t)N2 * sizeof(ST_int));
    ST_double *y = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *X = malloc((size_t)N2 * (size_t)K_x * sizeof(ST_double));
    ST_double *off = m->offset ? malloc((size_t)N2 * sizeof(ST_double)) : NULL;
    ST_double *w = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *mu = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *eta = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *irls_w = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *z_orig = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *aug = malloc((size_t)N2 * (size_t)K_aug * sizeof(ST_double));
    ST_double *previous = calloc((size_t)N2, sizeof(ST_double));
    ST_double *cert = calloc((size_t)N2, sizeof(ST_double));
    ST_double *xtx = malloc((size_t)K_x * (size_t)K_x * sizeof(ST_double));
    ST_double *xty = malloc((size_t)K_x * sizeof(ST_double));
    ST_double *beta = malloc((size_t)K_x * sizeof(ST_double));
    ST_double *inv_xx = malloc((size_t)K_x * (size_t)K_x * sizeof(ST_double));
    ST_double *xb_work = malloc((size_t)N2 * sizeof(ST_double));
    ST_double *fe_work = malloc((size_t)N2 * sizeof(ST_double));
    ST_int *sep_m = calloc((size_t)N2, sizeof(ST_int));
    /* the simplex check's positive-row weights: the copies, times the
     * fweight where the full fit uses it (aw and pw count each row once) */
    ST_double *w_pos = malloc((size_t)(N2 > 0 ? N2 : 1) * sizeof(ST_double));
    ST_int *keep_irls = NULL;
    ST_double *xtx_k = NULL, *xty_k = NULL, *beta_k = NULL, *inv_xx_k = NULL;
    if (!mask || !lv || !level_max || !row2 || !y || !X || (m->offset && !off) || !w || !mu ||
        !eta || !irls_w || !z_orig || !aug || !previous || !cert || !xtx || !xty || !beta ||
        !inv_xx || !xb_work || !fe_work || !sep_m || !w_pos)
    {
        rc = 920;
        goto done;
    }

    /* The rows alive, compact */
    for (ST_int r = 0, k = 0; r < N1; r++)
    {
        mask[r] = alive[r];
        if (!alive[r])
            continue;
        const ST_int i = rows[r];
        row2[k] = r;
        y[k] = m->y[i];
        for (ST_int c = 0; c < K_x; c++)
            X[(size_t)c * N2 + k] = m->X[(size_t)c * m->N + i];
        if (off)
            off[k] = m->offset[i];
        /* multiplicity times the user's weight, as in the expanded data;
         * when fweighted observations are drawn (unit_freq), the
         * multiplicity is the count of them drawn */
        w[k] = (double)mult[r] * ((m->w_user && !m->unit_freq) ? m->w_user[i] : 1.0);
        w_pos[k] = (double)mult[r] *
                   ((m->w_user && m->weight_type == 2 && !m->unit_freq) ? m->w_user[i] : 1.0);
        k++;
    }
    for (ST_int g = 0; g < G; g++)
    {
        level_max[g] = m->num_levels[g];
        if (!(lv[g] = malloc((size_t)N1 * sizeof(ST_int))))
        {
            rc = 920;
            goto done;
        }
        for (ST_int r = 0; r < N1; r++)
            lv[g][r] = m->levels[g][rows[r]];
    }
    if ((rc = hdfe_build_factors(&S, G, lv, level_max, N1, mask, N2, NULL, 1)))
        goto done;
    if ((rc = slope_factors(m, &S, rows, row2, N2)))
        goto done;
    S.K = K_aug;
    S.maxiter = m->maxiter;
    S.tolerance = m->tolerance;
    S.verbose = 0;
    S.num_threads = threads;
    S.factors_initialized = 1;
    S.has_weights = 1;
    S.weight_type = 1;
    S.weights = irls_w;
    S.sum_weights = 0.0;
    if (ctools_hdfe_alloc_buffers(&S, 0, K_aug) != 0)
    {
        rc = 920;
        goto done;
    }

    /* Separation checks the full fit ran: rows they find leave the
     * replicate, which is then screened again, as the full fit refits */
    if (options.check_simplex)
    {
        int nsep = ppml_simplex(&options, &S, y, X, w, 2, w_pos, N2, K_x, sep_m);
        if (nsep < 0)
        {
            /* The check did not converge (430): this replicate fails, as a
             * replicate whose IRLS does not converge; only memory aborts */
            if (nsep == -920)
                rc = 920;
            else
                out->status = HDFE_BOOT_NOCONV;
            goto done;
        }
        if (nsep)
        {
            for (ST_int k = 0; k < N2; k++)
                if (sep_m[k])
                    alive[row2[k]] = 0;
            *dropped = nsep;
            goto done;
        }
    }
    if (m->check_relu)
    {
        ST_int nsep = ppml_relu_ex(&options, &S, y, X, N2, K_x, sep_m, 1, cert);
        if (nsep < 0)
        {
            /* A projection that did not converge (430) or relu_strict()'s
             * iteration limit (9010): this replicate fails; only memory
             * aborts */
            if (nsep == -920)
                rc = 920;
            else
                out->status = HDFE_BOOT_NOCONV;
            goto done;
        }
        if (nsep)
        {
            for (ST_int k = 0; k < N2; k++)
                if (sep_m[k])
                    alive[row2[k]] = 0;
            *dropped = nsep;
            goto done;
        }
    }

    /* Warm start from the full fit's linear predictor */
    ST_double user_weight_scale = 0.0;
    for (ST_int k = 0; k < N2; k++)
        user_weight_scale += w[k] / N2;
    for (ST_int k = 0; k < N2; k++)
    {
        const ST_int i = rows[row2[k]];
        mu[k] = exp(m->eta[i] + (off ? off[k] : 0.0));
        if (!(mu[k] >= PPML_MU_FLOOR))
            mu[k] = PPML_MU_FLOOR;
        eta[k] = log(mu[k]) - (off ? off[k] : 0.0);
        irls_w[k] = mu[k] * (w[k] / user_weight_scale);
        if (irls_w[k] < DBL_MIN)
            irls_w[k] = DBL_MIN;
    }
    ppml_update_fe_weights(&S, irls_w, N2);

    /* Regressors collinear with the effects leave the IRLS, as in the fit */
    ST_int K_irls = K_x;
    if (options.remove_collinear)
    {
        for (ST_int k = 0; k < N2; k++)
            aug[k] = 0.0;
        memcpy(aug + N2, X, (size_t)N2 * (size_t)K_x * sizeof(ST_double));
        HDFE_SolveResult projection = ppml_partial(&options, &S, aug, N2, K_aug, S.num_threads);
        if (projection.status == 430)
        {
            out->status = HDFE_BOOT_NOCONV;
            goto done;
        }
        if (projection.status)
        {
            rc = projection.status;
            goto done;
        }
        ST_double *xtx_pre = calloc((size_t)K_x * (size_t)K_x, sizeof(ST_double));
        ST_int *pre = calloc((size_t)K_x, sizeof(ST_int));
        if (!xtx_pre || !pre)
        {
            free(xtx_pre);
            free(pre);
            rc = 920;
            goto done;
        }
        ctools_matmul_ata(aug + N2, N2, K_x, xtx_pre);
        detect_collinearity(xtx_pre, K_x, pre, 0);
        ST_int n_coll = 0;
        for (ST_int c = 0; c < K_x; c++)
            n_coll += pre[c] != 0;
        K_irls = K_x - n_coll;
        if (n_coll > 0 && K_irls > 0)
        {
            keep_irls = malloc((size_t)K_irls * sizeof(ST_int));
            xtx_k = malloc((size_t)K_irls * (size_t)K_irls * sizeof(ST_double));
            xty_k = malloc((size_t)K_irls * sizeof(ST_double));
            beta_k = malloc((size_t)K_irls * sizeof(ST_double));
            inv_xx_k = malloc((size_t)K_irls * (size_t)K_irls * sizeof(ST_double));
            if (!keep_irls || !xtx_k || !xty_k || !beta_k || !inv_xx_k)
                rc = 920;
            for (ST_int c = 0, j = 0; c < K_x && !rc; c++)
                if (!pre[c])
                    keep_irls[j++] = c;
        }
        free(xtx_pre);
        free(pre);
        if (rc)
            goto done;
    }

    io.report = 0;
    io.iteration_log = 0;
    io.options = &options;
    io.S = &S;
    io.y = y;
    io.X = X;
    io.offset = off;
    io.w_user = w;
    io.user_weight_scale = user_weight_scale;
    io.tolerance = m->tolerance;
    io.irls_tol = m->irls_tol;
    io.stdev_y = m->stdev_y;
    io.inner_tol_override = m->inner_tol_override;
    io.N = N2;
    io.K_x = K_x;
    io.K_aug = K_aug;
    io.G = G;
    io.irls_maxiter = m->irls_maxiter;
    io.num_threads = S.num_threads;
    io.verbose = 0;
    io.has_weights = 1;
    io.weight_type = 2;
    io.K_irls = K_irls;
    io.irls_keep_idx = keep_irls;
    io.mu = mu;
    io.eta = eta;
    io.irls_w = irls_w;
    io.z_orig = z_orig;
    io.aug_copy = aug;
    io.previous_working = previous;
    io.xtx = xtx;
    io.xty = xty;
    io.irls_beta = beta;
    io.inv_xx = inv_xx;
    io.xb_work = xb_work;
    io.fe_work = fe_work;
    io.xtx_k = xtx_k;
    io.xty_k = xty_k;
    io.beta_k = beta_k;
    io.inv_xx_k = inv_xx_k;
    io.sep_m = sep_m;
    memset(sep_m, 0, (size_t)N2 * sizeof(ST_int));
    rc = ppml_irls_loop(&io);
    if (rc == 430 || (!rc && !io.converged))
    {
        rc = 0;
        out->status = HDFE_BOOT_NOCONV;
        goto done;
    }
    if (rc)
        goto done;
    out->iterations = io.iterations;

    /* Rows the mu check separated leave the replicate, which is refit */
    if (options.check_mu)
    {
        ST_int nsep = 0;
        for (ST_int k = 0; k < N2; k++)
        {
            if (sep_m[k])
            {
                alive[row2[k]] = 0;
                nsep++;
            }
        }
        if (nsep)
        {
            *dropped = nsep;
            goto done;
        }
    }

    rc = ppml_final_collinear(&io);
    if (rc == 430)
    {
        rc = 0;
        out->status = HDFE_BOOT_NOCONV;
        goto done;
    }
    if (rc)
        goto done;
    for (ST_int c = 0; c < K_x; c++)
    {
        if (io.is_collinear[c] != m->is_collinear[c])
        {
            /* As the bootstrap prefix: a replicate whose omitted
             * coefficients differ from the full sample's is missing */
            out->status = HDFE_BOOT_COLLINEAR;
            goto done;
        }
    }
    rc = ppml_final_solve(&io);
    if (rc == 498)
    {
        rc = 0;
        out->status = HDFE_BOOT_COLLINEAR;
        goto done;
    }
    if (rc == 430)
    {
        rc = 0;
        out->status = HDFE_BOOT_NOCONV;
        goto done;
    }
    if (rc)
        goto done;

    /* b in the units of the data and the constant, as the full fit posts
     * them: mean(log mu) - xbar'b - mean(offset) + log(sd(y)) */
    {
        const ST_int K_keep = io.K_keep;
        dd_real sw = {0.0, 0.0}, slog = {0.0, 0.0}, soff = {0.0, 0.0};
        ST_double sum_w = 0.0, constant;
        for (ST_int k = 0; k < N2; k++)
            sum_w += irls_w[k];
        for (ST_int k = 0; k < N2; k++)
        {
            sw = dd_add_d(sw, irls_w[k]);
            slog = dd_add_d(slog, irls_w[k] * log(mu[k]));
            if (off)
                soff = dd_add_d(soff, irls_w[k] * off[k]);
        }
        const ST_double wsum = sw.hi + sw.lo;
        constant = (slog.hi + slog.lo) / wsum;
        for (ST_int c = 0; c < K_keep; c++)
        {
            const ST_double *xc = X + (size_t)io.keep_idx[c] * N2;
            ST_double wmean = 0.0;
            for (ST_int k = 0; k < N2; k++)
                wmean += irls_w[k] * xc[k];
            constant -= (wmean / sum_w) * io.beta_final[c];
            out->b[c] = io.beta_final[c] / m->stdev_x[io.keep_idx[c]];
        }
        if (off)
            constant -= (soff.hi + soff.lo) / wsum;
        constant += log(m->stdev_y);
        out->b[K_keep] = constant;
        out->status = HDFE_BOOT_OK;
    }
    if (out->fe || out->d)
        rc = replicate_effects(m, &S, rows, row2, N2, &io, out);

done:
    free(io.beta_prev);
    free(io.eta_refit);
    free(io.eta_refit_prev);
    free(io.is_collinear);
    free(io.keep_idx);
    free(io.beta_final);
    free(io.w_reg);
    free(io.xtx_keep);
    free(io.xty_keep);
    free(io.inv_xx_final);
    free(io.V_keep);
    free(io.xty_dd);
    S.weights = NULL;
    if (S.factors)
        ctools_hdfe_state_cleanup(&S);
    for (ST_int g = 0; lv && g < G; g++)
        free(lv[g]);
    free(lv);
    free(mask);
    free(level_max);
    free(row2);
    free(y);
    free(X);
    free(off);
    free(w);
    free(mu);
    free(eta);
    free(irls_w);
    free(z_orig);
    free(aug);
    free(previous);
    free(cert);
    free(xtx);
    free(xty);
    free(beta);
    free(inv_xx);
    free(xb_work);
    free(fe_work);
    free(sep_m);
    free(w_pos);
    free(keep_irls);
    free(xtx_k);
    free(xty_k);
    free(beta_k);
    free(inv_xx_k);
    return rc;
}

/* Screen-and-fit passes per replicate run until no separation check
 * removes rows, as the full fit refits: each pass that continues removes at
 * least one of the N1 drawn rows, so there are at most N1 + 1 */

static ST_retcode fit_replicate(const void *mv, void *ws, const ST_int *unit_mult,
                                hdfe_boot_rep *out)
{
    const ppml_boot_model *m = mv;
    const ST_int N = m->N, G = m->G;
    (void)ws;
    ST_retcode rc = 0;
    ST_int *rows = malloc((size_t)N * sizeof(ST_int));
    ST_int *mult = malloc((size_t)N * sizeof(ST_int));
    unsigned char *alive = malloc((size_t)N);
    ST_int **count = calloc((size_t)G, sizeof(ST_int *));
    unsigned char **positive = calloc((size_t)G, sizeof(unsigned char *));
    if (!rows || !mult || !alive || !count || !positive)
        rc = 920;
    for (ST_int g = 0; g < G && !rc; g++)
    {
        count[g] = malloc((size_t)m->num_levels[g] * sizeof(ST_int));
        positive[g] = malloc((size_t)m->num_levels[g]);
        if (!count[g] || !positive[g])
            rc = 920;
    }
    ST_int N1 = 0;
    for (ST_int i = 0; i < N && !rc; i++)
    {
        ST_int k = unit_mult[m->unit_of_row[i]];
        if (k <= 0)
            continue;
        rows[N1] = i;
        mult[N1] = k;
        alive[N1] = 1;
        N1++;
    }
    for (ST_int v = 0; out->fe && v < m->L_total; v++)
        out->fe[v] = NAN;
    for (ST_int i = 0; out->d && i < N; i++)
        out->d[i] = NAN;
    out->status = HDFE_BOOT_OTHER;
    for (ST_int pass = 0; pass <= N1 && !rc; pass++)
    {
        ST_int N2 = N1 ? screen(m, rows, mult, N1, alive, count, positive) : 0;
        if (N2 == 0)
        {
            out->status = HDFE_BOOT_EMPTY;
            break;
        }
        ST_int dropped = 0;
        rc = fit_rows(m, rows, mult, N1, alive, N2, out, &dropped);
        if (rc || !dropped)
            break;
        out->status = HDFE_BOOT_OTHER;
    }
    for (ST_int g = 0; count && g < G; g++)
        free(count[g]);
    for (ST_int g = 0; positive && g < G; g++)
        free(positive[g]);
    free(count);
    free(positive);
    free(rows);
    free(mult);
    free(alive);
    return rc;
}

static size_t fit_bytes(const void *mv)
{
    const ppml_boot_model *m = mv;
    /* X, its partialled copy and the augmented columns, with the simplex
     * check's two more copies of X, and about 24 row vectors */
    return (size_t)m->N * ((size_t)m->K_x * 4 + 24) * sizeof(ST_double) +
           (size_t)m->N * (size_t)m->G * 3 * sizeof(ST_int);
}

static void *fit_new(const void *mv)
{
    (void)mv;
    return malloc(1);
}

static void fit_free(void *ws)
{
    free(ws);
}

/* Each further thread of a replicate's projections: three CG buffers and
 * the per-level means */
static size_t fit_thread_bytes(const void *mv)
{
    const ppml_boot_model *m = mv;
    return (size_t)m->N * 3 * sizeof(ST_double) + (size_t)m->L_total * sizeof(ST_double);
}

static const hdfe_boot_model_ops fit_ops = {fit_bytes, fit_new, fit_free, fit_replicate,
                                             fit_thread_bytes};

/* Memory for replicates in flight: half of physical memory */
static size_t boot_budget(void)
{
    uint64_t physical = ctools_physical_memory_bytes();
    return physical ? (size_t)(physical / 2) : (size_t)1 << 31;
}

/* Scalar __cppmlhdfe_bs_par: the replicates boot_chunk runs at once, so that
 * the ado sizes its chunks as multiples of it */
static ST_retcode save_parallel(const ppml_boot_model *m)
{
    ST_int par = hdfe_boot_parallel(&m->design, &fit_ops, m, m->K_keep + 1, 0, 0,
                                    ctools_get_max_threads(), boot_budget());
    return ctools_scal_save("__cppmlhdfe_bs_par", (ST_double)(par > 0 ? par : 1));
}

/* ========================================================================
 * Plugin calls
 * ======================================================================== */

static ST_retcode check_retained(void)
{
    if (!g_pboot)
        return 198;
    if (SF_nobs() != g_pboot->nobs)
    {
        SF_error("cppmlhdfe: the data changed between the fit and its bootstrap\n");
        return 459;
    }
    return 0;
}

/* Value q (0-based) of those the ado packed into the call's variables
 * (_ctools_hdfe_boot.ado): variable q / nobs + 1, observation q % nobs + 1 */
static ST_retcode packed_value(ST_int nobs, size_t q, ST_double *x)
{
    return SF_vdata((ST_int)(q / (size_t)nobs) + 1, (ST_int)(q % (size_t)nobs) + 1, x);
}

/* The call passed exactly the variables that n packed values fill */
static int packed_vars_ok(size_t n)
{
    const ST_int nobs = SF_nobs();
    return nobs > 0 && n > 0 && (size_t)SF_nvars() == (n + (size_t)nobs - 1) / (size_t)nobs;
}

ST_retcode cppmlhdfe_boot_chunk(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    ppml_boot_model *m = g_pboot;
    int first = -1, count = -1, total = -1;
    if (ctools_parse_int_option(args, "first", &first) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "count", &count) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "total", &total) != CTOOLS_PARSE_OK || first < 0 ||
        count < 1 || total < 2 || first + count > total)
    {
        SF_error("cppmlhdfe: invalid boot_chunk arguments\n");
        return 198;
    }
    const int stata_draws = ctools_parse_bool_option(args, "stata");
    const int supplied_draws = ctools_parse_bool_option(args, "supplied");
    if (first == 0)
    {
        hdfe_boot_acc_free(&m->acc);
        m->fe_on = ctools_parse_bool_option(args, "fe");
        m->dse_on = ctools_parse_bool_option(args, "dse");
        if ((rc = hdfe_boot_acc_init(&m->acc, m->K_keep + 1, m->fe_on ? m->L_total : 0,
                                     m->dse_on ? m->N : 0, total)))
            return rc;
        if (m->fe_on && ctools_parse_bool_option(args, "ci"))
        {
            rc = hdfe_boot_acc_keep_fe(&m->acc, boot_budget());
            if (rc == 909)
                SF_error("cppmlhdfe: fevce(bootstrap, ci) keeps every replicate's fixed effects, "
                         "and they do not fit in memory; use fewer reps() or omit ci\n");
            if (rc)
                return rc;
        }
        if (!stata_draws && !supplied_draws)
        {
            uint64_t seed = 0;
            if (ctools_parse_seed_option(args, &seed) != CTOOLS_PARSE_OK)
            {
                SF_error("cppmlhdfe: boot_chunk needs seedhi= and seedlo=\n");
                return 198;
            }
            m->design.seed = seed;
        }
    }
    if (m->acc.R_cap != total || m->acc.R != first)
        return 198;

    ST_int *supplied = NULL;
    if (stata_draws || supplied_draws)
    {
        /* draws(stata): bsample's uniforms, mapped to unit multiplicities
         * (on the fweight-expanded frame when fweighted rows are resampled);
         * _bsdraws(): the multiplicities themselves */
        if ((rc = hdfe_boot_read_draws(&m->design, m->cluster_mode, supplied_draws, count, m->N,
                                       m->unit_of_row, m->obs_of_row, "cppmlhdfe", &supplied)))
            return rc;
    }

    rc = hdfe_boot_run(&m->design, supplied, &fit_ops, m, first, count, ctools_get_max_threads(),
                       boot_budget(), &m->acc);
    free(supplied);
    if (rc == 909)
        SF_error("cppmlhdfe: not enough memory for one bootstrap replicate\n");
    if (rc)
        return rc;
    if (ctools_parse_bool_option(args, "rows"))
        rc = hdfe_boot_store_rows(&m->acc, first, count, "__cppmlhdfe_bs_rows");
    for (int r = 0; r < count && !rc; r++)
        rc = ctools_mat_store("__cppmlhdfe_bs_status", 1, r + 1,
                              (ST_double)m->acc.status[first + r]);
    return rc;
}

static ST_retcode store_vector(const char *name, const ST_double *v, ST_int n, ST_int rows)
{
    ST_retcode rc = 0;
    for (ST_int j = 0; j < n && !rc; j++)
    {
        for (ST_int r = 0; r < rows && !rc; r++)
        {
            ST_double x = v[(size_t)r * n + j];
            rc = ctools_mat_store(name, r + 1, j + 1, isfinite(x) ? x : SV_missval);
        }
    }
    return rc;
}

static int cmp_int(const void *a, const void *b)
{
    ST_int x = *(const ST_int *)a, y = *(const ST_int *)b;
    return (x > y) - (x < y);
}

/* fevce(bootstrap): the SE of each saved effect (the effects flagged in
 * matrix __cppmlhdfe_bs_fe_targets, 1 x G) in variables 1.. of the call, then
 * their valid-replicate counts with nreps, then d's SE with dse, on the
 * estimation sample; matrix __cppmlhdfe_bs_fe_table (G x 5: levels, never
 * estimable, always estimable, min and median valid replicates) and scalars
 * __cppmlhdfe_bs_fe_split_reps, __cppmlhdfe_bs_fe_N_reps */
static ST_retcode finish_effects(const ppml_boot_model *m, const char *args)
{
    const ST_int G = m->G, N = m->N;
    const int nreps = ctools_parse_bool_option(args, "nreps");
    char tmat[] = "__cppmlhdfe_bs_fe_targets";
    int target[10] = {0}, T = 0;
    for (ST_int g = 0; g < G && m->fe_on; g++)
    {
        ST_double x = 0;
        if (SF_mat_el(tmat, 1, g + 1, &x))
            return 198;
        target[g] = x != 0;
        T += target[g];
    }
    const int ci = ctools_parse_bool_option(args, "ci");
    const ST_int nd = T * (1 + nreps + 2 * ci); /* d's SE follows the per-effect outputs */
    if (SF_nvars() != nd + m->dse_on || (ci && !m->acc.fe_draws))
        return 198;
    ST_retcode rc = 0;
    const size_t Lc = (size_t)(m->L_total > 0 ? m->L_total : 1);
    ST_double *fe_se = malloc((size_t)(m->L_total > 0 ? m->L_total : 1) * sizeof(ST_double));
    ST_double *lb = ci ? malloc(Lc * sizeof(ST_double)) : NULL;
    ST_double *ub = ci ? malloc(Lc * sizeof(ST_double)) : NULL;
    if (ci && (!lb || !ub))
        rc = 920;
    if (!rc && ci)
        rc = hdfe_boot_fe_ci(&m->acc, ctools_parse_double_option(args, "level", 95.0), lb, ub);
    ST_double *d_se = m->dse_on ? malloc((size_t)N * sizeof(ST_double)) : NULL;
    ST_int *sorted = malloc((size_t)(m->L_total > 0 ? m->L_total : 1) * sizeof(ST_int));
    if (!rc && (!fe_se || (m->dse_on && !d_se) || !sorted))
        rc = 920;
    if (!rc)
        hdfe_boot_fe_se(&m->acc, m->fe_on ? fe_se : NULL, d_se);
    for (ST_int i = 0; i < N && !rc; i++)
    {
        ST_int col = 0;
        for (ST_int g = 0; g < G && !rc && m->fe_on; g++)
        {
            if (!target[g])
                continue;
            const ST_int node = m->node_off[g] + m->levels[g][i] - 1;
            const ST_double se = fe_se[node];
            rc = SF_vstore(col + 1, m->obs_of_row[i], isfinite(se) ? se : SV_missval);
            if (!rc && nreps)
                rc = SF_vstore(T + col + 1, m->obs_of_row[i], (ST_double)m->acc.fe_n[node]);
            if (!rc && ci)
                rc = SF_vstore(T * (1 + nreps) + col + 1, m->obs_of_row[i],
                               isfinite(lb[node]) ? lb[node] : SV_missval);
            if (!rc && ci)
                rc = SF_vstore(T * (2 + nreps) + col + 1, m->obs_of_row[i],
                               isfinite(ub[node]) ? ub[node] : SV_missval);
            col++;
        }
        if (!rc && m->dse_on)
            rc = SF_vstore(nd + 1, m->obs_of_row[i], isfinite(d_se[i]) ? d_se[i] : SV_missval);
    }
    for (ST_int g = 0; g < G && !rc && m->fe_on; g++)
    {
        const ST_int L = m->num_levels[g];
        ST_int never = 0, always = 0;
        for (ST_int l = 0; l < L; l++)
        {
            const ST_int n = m->acc.fe_n[m->node_off[g] + l];
            sorted[l] = n;
            never += n == 0;
            always += n == m->acc.R;
        }
        qsort(sorted, (size_t)L, sizeof(ST_int), cmp_int);
        const ST_double median = L % 2 ? sorted[L / 2] : 0.5 * (sorted[L / 2 - 1] + sorted[L / 2]);
        const ST_double row[5] = {(ST_double)L, (ST_double)never, (ST_double)always,
                                  (ST_double)sorted[0], median};
        for (int c = 0; c < 5 && !rc; c++)
            rc = ctools_mat_store("__cppmlhdfe_bs_fe_table", g + 1, c + 1, row[c]);
    }
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_bs_fe_split_reps", (ST_double)m->acc.fe_split_reps);
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_bs_fe_N_reps", (ST_double)m->acc.R);
    free(fe_se);
    free(lb);
    free(ub);
    free(d_se);
    free(sorted);
    return rc;
}

ST_retcode cppmlhdfe_boot_finish(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    ppml_boot_model *m = g_pboot;
    if ((m->fe_on || m->dse_on) && (rc = finish_effects(m, args)))
        return rc;
    if (!ctools_parse_bool_option(args, "coef"))
        return 0;
    double level = ctools_parse_double_option(args, "level", 95.0);
    int mse = ctools_parse_bool_option(args, "mse");
    hdfe_boot_summary s;
    rc = hdfe_boot_summarize(&m->acc, m->b_hat, mse, level, &s);
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_bs_N_reps", (ST_double)s.N_reps);
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_bs_N_misreps", (ST_double)s.N_misreps);
    if (rc == 2000)
    {
        SF_error("cppmlhdfe: insufficient observations to compute bootstrap standard errors\n");
        return rc;
    }
    const ST_int Kb = m->K_keep + 1;
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_V", s.V, Kb, Kb);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_b", s.b_bs, Kb, 1);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_se", s.se, Kb, 1);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_bias", s.bias, Kb, 1);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_z0", s.z0, Kb, 1);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_ci_normal", s.ci_normal, Kb, 2);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_ci_percentile", s.ci_percentile, Kb, 2);
    if (!rc)
        rc = store_vector("__cppmlhdfe_bs_ci_bc", s.ci_bc, Kb, 2);
    for (int k = 0; k < HDFE_BOOT_NSTATUS && !rc; k++)
        rc = ctools_mat_store("__cppmlhdfe_bs_counts", 1, k + 1,
                              (ST_double)m->acc.status_count[k]);
    hdfe_boot_summary_free(&s);
    return rc;
}

/* ========================================================================
 * wildbootstrap(): the score bootstrap of the coefficients (plan 5.9, 8.4)
 * ======================================================================== */

typedef struct
{
    double v;
    ST_int row;
} value_row;

static int value_row_cmp(const void *a, const void *b)
{
    const value_row *x = a, *y = b;
    if (x->v != y->v)
        return x->v < y->v ? -1 : 1;
    return (x->row > y->row) - (x->row < y->row);
}

/* 0-based ids of plugin variable var on the rows, in ascending order of its
 * values (as egen group() and boottest number clusters) */
static ST_retcode dense_ids(const ppml_boot_model *m, ST_int var, ST_int *ids, ST_int *n)
{
    const ST_int N = m->N;
    value_row *vr = malloc((size_t)N * sizeof(*vr));
    if (!vr)
        return 920;
    for (ST_int i = 0; i < N; i++)
    {
        double x;
        ST_retcode rc = SF_vdata(var, m->obs_of_row[i], &x);
        if (rc || SF_is_missing(x))
        {
            free(vr);
            return rc ? rc : 198;
        }
        vr[i].v = x;
        vr[i].row = i;
    }
    qsort(vr, (size_t)N, sizeof(*vr), value_row_cmp);
    ST_int k = -1;
    for (ST_int i = 0; i < N; i++)
    {
        if (i == 0 || vr[i].v != vr[i - 1].v)
            k++;
        ids[vr[i].row] = k;
    }
    *n = k + 1;
    free(vr);
    return 0;
}

static void wild_free(ppml_boot_model *m)
{
    free(m->wb_vce);
    free(m->wb_boot);
    free(m->wb_test);
    free(m->wb_scale);
    m->wb_vce = m->wb_boot = m->wb_test = NULL;
    m->wb_scale = NULL;
    hdfe_wild_score_free(&m->wb_score);
}

/* Solve the n x n symmetric positive definite system M x = r (column-major)
 * by Cholesky; 198 when M is not positive definite */
static ST_retcode spd_solve(ST_int n, const ST_double *M, const ST_double *r, ST_double *x)
{
    ST_double *L = calloc((size_t)n * (size_t)n, sizeof(ST_double));
    if (!L)
        return 920;
    ST_retcode rc = 0;
#define IX(i, j) ((size_t)(i) + (size_t)(j) * (size_t)n)
    for (ST_int j = 0; j < n && !rc; j++)
    {
        ST_double s = M[IX(j, j)];
        for (ST_int k = 0; k < j; k++)
            s -= L[IX(j, k)] * L[IX(j, k)];
        if (!(s > 1e-14 * M[IX(j, j)]) || !isfinite(s))
        {
            rc = 198;
            break;
        }
        L[IX(j, j)] = sqrt(s);
        for (ST_int i = j + 1; i < n; i++)
        {
            ST_double t = M[IX(i, j)];
            for (ST_int k = 0; k < j; k++)
                t -= L[IX(i, k)] * L[IX(j, k)];
            L[IX(i, j)] = t / L[IX(j, j)];
        }
    }
    for (ST_int i = 0; i < n && !rc; i++)
    {
        ST_double t = r[i];
        for (ST_int k = 0; k < i; k++)
            t -= L[IX(i, k)] * x[k];
        x[i] = t / L[IX(i, i)];
    }
    for (ST_int i = n - 1; i >= 0 && !rc; i--)
    {
        ST_double t = x[i];
        for (ST_int k = i + 1; k < n; k++)
            t -= L[IX(k, i)] * x[k];
        x[i] = t / L[IX(i, i)];
    }
#undef IX
    free(L);
    return rc;
}

/*
 * Restricted fits for the score bootstrap. An evaluator holds the full
 * sample's factors and the IRLS arrays, built once and reused for every
 * fit of the tests and of their interval searches; the search runs several
 * evaluators side by side.
 */
typedef struct
{
    HDFE_State S;
    ST_int cap;   /* threads its projection buffers serve */
    ST_int K_aug; /* the tested column, then the other kept regressors */
    ST_double *Xr, *off, *y, *mu, *eta, *irls_w, *z_orig, *aug, *previous;
    ST_double *xtx, *xty, *beta, *inv_xx, *xb_work, *fe_work, *gam, *M, *rhs, *e, *t;
    ST_int *sep_m;
} ppml_reval;

static void reval_free(ppml_reval *ev)
{
    if (!ev)
        return;
    ev->S.weights = NULL;
    if (ev->S.factors)
        ctools_hdfe_state_cleanup(&ev->S);
    free(ev->Xr);
    free(ev->off);
    free(ev->y);
    free(ev->mu);
    free(ev->eta);
    free(ev->irls_w);
    free(ev->z_orig);
    free(ev->aug);
    free(ev->previous);
    free(ev->xtx);
    free(ev->xty);
    free(ev->beta);
    free(ev->inv_xx);
    free(ev->xb_work);
    free(ev->fe_work);
    free(ev->gam);
    free(ev->M);
    free(ev->rhs);
    free(ev->e);
    free(ev->t);
    free(ev->sep_m);
    free(ev);
}

/* An evaluator's memory with projection buffers for `threads` threads, B
 * replicates */
static size_t reval_bytes(const ppml_boot_model *m, ST_int threads, ST_int B)
{
    const size_t N = (size_t)m->N, Kr = m->K_keep > 1 ? (size_t)m->K_keep - 1 : 1;
    const size_t cap = (size_t)(threads < 1 ? 1 : (threads < m->K_keep ? threads : m->K_keep));
    return N * (2 * Kr + 11) * sizeof(ST_double) + N * (1 + (size_t)m->G) * sizeof(ST_int) +
           cap * 3 * N * sizeof(ST_double) +
           (size_t)m->L_total * (4 + cap) * sizeof(ST_double) + ((size_t)B + 1) * sizeof(ST_double);
}

static ppml_reval *reval_new(const ppml_boot_model *m, ST_int threads, ST_int B)
{
    const ST_int N = m->N, G = m->G, Kr = m->K_keep - 1;
    const ST_int K_aug = Kr + 1, KK = Kr > 0 ? Kr : 1;
    ppml_reval *ev = calloc(1, sizeof(*ev));
    if (!ev || Kr < 0)
    {
        free(ev);
        return NULL;
    }
    ev->K_aug = K_aug;
    ev->Xr = malloc((size_t)N * (size_t)KK * sizeof(ST_double));
    ev->off = malloc((size_t)N * sizeof(ST_double));
    ev->y = malloc((size_t)N * sizeof(ST_double));
    ev->mu = malloc((size_t)N * sizeof(ST_double));
    ev->eta = malloc((size_t)N * sizeof(ST_double));
    ev->irls_w = malloc((size_t)N * sizeof(ST_double));
    ev->z_orig = malloc((size_t)N * sizeof(ST_double));
    ev->aug = malloc((size_t)N * (size_t)K_aug * sizeof(ST_double));
    ev->previous = calloc((size_t)N, sizeof(ST_double));
    ev->xtx = malloc((size_t)KK * (size_t)KK * sizeof(ST_double));
    ev->xty = malloc((size_t)KK * sizeof(ST_double));
    ev->beta = malloc((size_t)KK * sizeof(ST_double));
    ev->inv_xx = malloc((size_t)KK * (size_t)KK * sizeof(ST_double));
    ev->xb_work = malloc((size_t)N * sizeof(ST_double));
    ev->fe_work = malloc((size_t)N * sizeof(ST_double));
    ev->gam = malloc((size_t)KK * sizeof(ST_double));
    ev->M = malloc((size_t)KK * (size_t)KK * sizeof(ST_double));
    ev->rhs = malloc((size_t)KK * sizeof(ST_double));
    ev->e = malloc((size_t)N * sizeof(ST_double));
    ev->t = malloc(((size_t)B + 1) * sizeof(ST_double));
    ev->sep_m = calloc((size_t)N, sizeof(ST_int));
    ST_int *mask = malloc((size_t)N * sizeof(ST_int));
    ST_int *level_max = malloc((size_t)G * sizeof(ST_int));
    if (!ev->Xr || !ev->off || !ev->y || !ev->mu || !ev->eta || !ev->irls_w || !ev->z_orig ||
        !ev->aug || !ev->previous || !ev->xtx || !ev->xty || !ev->beta || !ev->inv_xx ||
        !ev->xb_work || !ev->fe_work || !ev->gam || !ev->M || !ev->rhs || !ev->e || !ev->t ||
        !ev->sep_m || !mask || !level_max)
    {
        free(mask);
        free(level_max);
        reval_free(ev);
        return NULL;
    }
    for (ST_int i = 0; i < N; i++)
    {
        mask[i] = 1;
        ev->y[i] = m->y[i];
    }
    for (ST_int g = 0; g < G; g++)
        level_max[g] = m->num_levels[g];
    ST_retcode rc = hdfe_build_factors(&ev->S, G, m->levels, level_max, N, mask, N, NULL, 1);
    free(mask);
    free(level_max);
    if (!rc)
        rc = slope_factors(m, &ev->S, NULL, NULL, N);
    if (rc)
    {
        reval_free(ev);
        return NULL;
    }
    HDFE_State *S = &ev->S;
    S->K = K_aug;
    S->maxiter = m->maxiter;
    S->tolerance = m->tolerance;
    S->verbose = 0;
    S->num_threads = threads > 1 ? threads : 1;
    S->factors_initialized = 1;
    S->has_weights = 1;
    S->weight_type = 1;
    S->weights = ev->irls_w;
    if (ctools_hdfe_alloc_buffers(S, 0, K_aug) != 0)
    {
        reval_free(ev);
        return NULL;
    }
    ev->cap = S->num_threads;
    return ev;
}

/*
 * The fit under H0: b_col = b0 (data units) on the full fit's sample: x_col
 * leaves the regressors and enters the offset as x_col b0, the IRLS starts
 * from ev->mu, which receives the restricted means. Then ev->e, each row's
 * contribution to the efficient score of b_col: w (y - mu) r, r the residual
 * of x_col on the effects and the other kept regressors with weights w mu
 * (the Hessian of the unrestricted model at the restricted estimate, as
 * boottest takes it). *ok is 0 when the fit or the projection fails. The
 * projections use up to `threads` threads (one per column, so the result
 * does not depend on them).
 */
static ST_retcode restricted_fit(const ppml_boot_model *m, ppml_reval *ev, ST_int col, double b0,
                                 ST_int threads, int *ok)
{
    const ST_int N = m->N, G = m->G;
    *ok = 0;
    ST_int Kr = 0;
    for (ST_int c = 0; c < m->K_x; c++)
        Kr += !m->is_collinear[c] && c != col;
    if (Kr + 1 != ev->K_aug)
        return 198;
    const ST_int K_aug = ev->K_aug;
    PPML_Options options = m->options;
    options.check_mu = 0; /* the sample is the full fit's */
    HDFE_State *S = &ev->S;
    S->num_threads = threads < 1 ? 1 : (threads < ev->cap ? threads : ev->cap);
    S->tolerance = m->tolerance;
    ST_double *Xr = ev->Xr, *off = ev->off, *y = ev->y, *mu = ev->mu, *eta = ev->eta;
    ST_double *irls_w = ev->irls_w, *aug = ev->aug, *e = ev->e;
    ST_retcode rc = 0;
    ppml_irls_io io;
    memset(&io, 0, sizeof(io));
    const double bstd = b0 * m->stdev_x[col];
    for (ST_int i = 0; i < N; i++)
        off[i] = (m->offset ? m->offset[i] : 0.0) + m->X[(size_t)col * N + i] * bstd;
    for (ST_int c = 0, k = 0; c < m->K_x; c++)
    {
        if (m->is_collinear[c] || c == col)
            continue;
        memcpy(Xr + (size_t)k * N, m->X + (size_t)c * N, (size_t)N * sizeof(ST_double));
        k++;
    }
    const ST_double *w = m->w_user;
    ST_double scale = 1.0;
    if (w)
    {
        scale = 0.0;
        for (ST_int i = 0; i < N; i++)
            scale += w[i] / N;
    }
    for (ST_int i = 0; i < N; i++)
    {
        if (!(mu[i] >= PPML_MU_FLOOR))
            mu[i] = PPML_MU_FLOOR;
        eta[i] = log(mu[i]) - off[i];
        irls_w[i] = mu[i] * ((w ? w[i] : 1.0) / scale);
        if (irls_w[i] < DBL_MIN)
            irls_w[i] = DBL_MIN;
    }
    ppml_update_fe_weights(S, irls_w, N);

    io.report = 0;
    io.options = &options;
    io.S = S;
    io.y = y;
    io.X = Xr;
    io.offset = off;
    io.w_user = w;
    io.user_weight_scale = scale;
    io.tolerance = m->tolerance;
    io.irls_tol = m->irls_tol;
    io.stdev_y = m->stdev_y;
    io.inner_tol_override = m->inner_tol_override;
    io.N = N;
    io.K_x = Kr;
    io.K_aug = K_aug;
    io.G = G;
    io.irls_maxiter = m->irls_maxiter;
    io.num_threads = S->num_threads;
    io.has_weights = w != NULL;
    io.weight_type = m->weight_type;
    io.K_irls = Kr;
    io.mu = mu;
    io.eta = eta;
    io.irls_w = irls_w;
    io.z_orig = ev->z_orig;
    io.aug_copy = aug;
    io.previous_working = ev->previous;
    io.xtx = ev->xtx;
    io.xty = ev->xty;
    io.irls_beta = ev->beta;
    io.inv_xx = ev->inv_xx;
    io.xb_work = ev->xb_work;
    io.fe_work = ev->fe_work;
    io.sep_m = ev->sep_m;
    rc = ppml_irls_loop(&io);
    if (rc == 430 || (!rc && !io.converged))
    {
        rc = 0;
        goto done;
    }
    if (rc)
        goto done;

    /* x_col partialled on the effects and the other regressors, with the
     * weights at the restricted means */
    for (ST_int i = 0; i < N; i++)
    {
        irls_w[i] = mu[i] * ((w ? w[i] : 1.0) / scale);
        if (irls_w[i] < DBL_MIN)
            irls_w[i] = DBL_MIN;
    }
    ppml_update_fe_weights(S, irls_w, N);
    S->tolerance = io.target_inner_tol;
    memcpy(aug, m->X + (size_t)col * N, (size_t)N * sizeof(ST_double));
    memcpy(aug + N, Xr, (size_t)N * (size_t)Kr * sizeof(ST_double));
    HDFE_SolveResult pr = ppml_partial(&options, S, aug, N, K_aug, S->num_threads);
    if (pr.status == 430)
        goto done;
    if (pr.status)
    {
        rc = pr.status;
        goto done;
    }
    if (Kr > 0)
    {
        ST_double *M = ev->M, *rhs = ev->rhs, *gam = ev->gam;
        for (ST_int a = 0; a < Kr; a++)
        {
            const ST_double *xa = aug + (size_t)(a + 1) * N;
            ST_double s = 0.0;
            for (ST_int i = 0; i < N; i++)
                s += irls_w[i] * xa[i] * aug[i];
            rhs[a] = s;
            for (ST_int b = 0; b <= a; b++)
            {
                const ST_double *xb = aug + (size_t)(b + 1) * N;
                ST_double t = 0.0;
                for (ST_int i = 0; i < N; i++)
                    t += irls_w[i] * xa[i] * xb[i];
                M[a + (size_t)b * Kr] = M[b + (size_t)a * Kr] = t;
            }
        }
        ST_retcode src = spd_solve(Kr, M, rhs, gam);
        if (src == 198)
            goto done;
        if (src)
        {
            rc = src;
            goto done;
        }
        for (ST_int a = 0; a < Kr; a++)
        {
            const ST_double *xa = aug + (size_t)(a + 1) * N;
            for (ST_int i = 0; i < N; i++)
                aug[i] -= gam[a] * xa[i];
        }
    }
    for (ST_int i = 0; i < N; i++)
        e[i] = (w ? w[i] : 1.0) * (y[i] - mu[i]) * aug[i];
    *ok = 1;

done:
    free(io.beta_prev);
    free(io.eta_refit);
    free(io.eta_refit_prev);
    return rc;
}

ST_retcode cppmlhdfe_wild_setup(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    ppml_boot_model *m = g_pboot;
    wild_free(m);
    const ST_int N = m->N, K_x = m->K_x;
    const int has_vce = ctools_parse_bool_option(args, "vce");
    const int has_boot = ctools_parse_bool_option(args, "boot");
    if (SF_nvars() != has_vce + has_boot || m->K_keep < 1)
        return 198;
    char tmat[] = "__cppmlhdfe_wb_test", smat[] = "__cppmlhdfe_wb_scale";
    const ST_int T = SF_col(tmat);
    if (T < 1 || SF_col(smat) != T)
        return 198;
    m->wb_test = malloc((size_t)T * sizeof(ST_int));
    m->wb_scale = malloc((size_t)T * sizeof(ST_double));
    if (!m->wb_test || !m->wb_scale)
        return 920;
    for (ST_int j = 0; j < T && !rc; j++)
    {
        ST_double x, s;
        if ((rc = SF_mat_el(tmat, 1, j + 1, &x)) || (rc = SF_mat_el(smat, 1, j + 1, &s)))
            break;
        int p = (int)x;
        if (p < 1 || p > K_x || m->is_collinear[p - 1])
        {
            SF_error("cppmlhdfe: wildbootstrap() cannot test an omitted coefficient\n");
            rc = 198;
            break;
        }
        m->wb_test[j] = p - 1;
        m->wb_scale[j] = s;
    }
    m->wb_ntest = T;
    m->wb_G = m->wb_C = N;
    if (!rc && has_vce)
    {
        if (!(m->wb_vce = malloc((size_t)N * sizeof(ST_int))))
            rc = 920;
        else
            rc = dense_ids(m, 1, m->wb_vce, &m->wb_G);
    }
    if (!rc && (has_boot || has_vce))
    {
        if (!(m->wb_boot = malloc((size_t)N * sizeof(ST_int))))
            rc = 920;
        else if (has_boot)
            rc = dense_ids(m, 1 + has_vce, m->wb_boot, &m->wb_C);
        else
        {
            memcpy(m->wb_boot, m->wb_vce, (size_t)N * sizeof(ST_int));
            m->wb_C = m->wb_G;
        }
    }
    if (rc)
    {
        wild_free(m);
        return rc;
    }
    return ctools_scal_save("__cppmlhdfe_wb_C", (ST_double)m->wb_C);
}

ST_retcode cppmlhdfe_wild_chunk(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    ppml_boot_model *m = g_pboot;
    int first = -1, count = -1, total = -1;
    if (ctools_parse_int_option(args, "first", &first) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "count", &count) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "total", &total) != CTOOLS_PARSE_OK || first < 0 ||
        count < 1 || total < 1 || first + count > total || !m->wb_test)
    {
        SF_error("cppmlhdfe: invalid wild_chunk arguments\n");
        return 198;
    }
    const int enumerate = ctools_parse_bool_option(args, "enum");
    if (first == 0)
    {
        hdfe_wild_score_free(&m->wb_score);
        if ((rc = hdfe_wild_score_setup(m->N, m->wb_vce, m->wb_G, m->wb_boot, m->wb_C,
                                        m->w_user, total, enumerate, &m->wb_score)))
            return rc;
    }
    if (m->wb_score.B != total)
        return 198;
    if (enumerate)
        return 0;
    /* The ado packs the chunk's weights replicate after replicate */
    const size_t values = (size_t)count * (size_t)m->wb_C;
    if (!packed_vars_ok(values))
        return 198;
    const ST_int nobs = SF_nobs();
    ST_double *v = malloc(values * sizeof(ST_double));
    if (!v)
        return 920;
    for (size_t q = 0; q < values && !rc; q++)
        rc = packed_value(nobs, q, &v[q]);
    if (!rc)
        rc = hdfe_wild_score_weights(&m->wb_score, v, first, count);
    free(v);
    return rc;
}

/* The interval search's grid: 2 half + 1 points, evaluated at once */
#define PPML_WILD_HALF 4
#define PPML_WILD_NPT (2 * PPML_WILD_HALF + 1)

typedef struct
{
    const ppml_boot_model *m;
    ST_int col;
    int ptype;
    ppml_reval **ev; /* P evaluators */
    ST_int P, inner;
    const ST_double *mu0; /* the full fit's means: where every chain starts */
    ST_double **chain;    /* each chain's last restricted means, or NULL:
                             every fit starts from mu0 */
    int n_chain;
    ST_retcode rc;
} wild_pctx;

/*
 * p(b0) at a batch of points (hdfe_wild_invert_batch): a restricted fit
 * started from its chain's last means, then the replicates' statistics;
 * NaN when the fit fails. The points run side by side, an evaluator each,
 * with the threads left over for each fit's own parallelism. A point's p
 * depends only on b0 and its chain, not on the threads.
 */
static void wild_p_batch(void *ctx, const double *b0, const int *from, const int *to, int n,
                         double *pv)
{
    wild_pctx *c = ctx;
    const ST_int N = c->m->N;
    ST_retcode rcs[PPML_WILD_NPT];
    for (int k = 0; k < n; k++)
    {
        pv[k] = NAN;
        rcs[k] = 0;
    }
    if (c->rc || n > PPML_WILD_NPT)
    {
        if (!c->rc)
            c->rc = 198;
        return;
    }
    const ST_int nthr = n < c->P ? n : c->P;
    (void)nthr;
#ifdef _OPENMP
    const int levels = omp_get_max_active_levels();
    if (nthr > 1 && c->inner > 1 && levels < 2)
        omp_set_max_active_levels(2);
#pragma omp parallel num_threads(nthr) if (nthr > 1)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        omp_set_num_threads(c->inner);
#endif
        ppml_reval *ev = c->ev[tid];
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int k = 0; k < n; k++)
        {
            const ST_double *start =
                (c->chain && from[k] >= 0 && from[k] < c->n_chain) ? c->chain[from[k]] : c->mu0;
            memcpy(ev->mu, start, (size_t)N * sizeof(ST_double));
            int ok = 0;
            rcs[k] = restricted_fit(c->m, ev, c->col, b0[k], c->inner, &ok);
            if (!rcs[k] && ok)
            {
                rcs[k] = hdfe_wild_score_t(&c->m->wb_score, ev->e, c->inner, ev->t);
                if (!rcs[k])
                    pv[k] = hdfe_wild_pvalue_t(ev->t, c->m->wb_score.B, c->ptype, NULL);
                if (!rcs[k] && c->chain && to[k] >= 0 && to[k] < c->n_chain)
                    memcpy(c->chain[to[k]], ev->mu, (size_t)N * sizeof(ST_double));
            }
        }
    }
#ifdef _OPENMP
    omp_set_max_active_levels(levels);
#endif
    for (int k = 0; k < n && !c->rc; k++)
        c->rc = rcs[k];
}

static ST_double stata_value(double x)
{
    return isfinite(x) ? x : SV_missval;
}

ST_retcode cppmlhdfe_wild_finish(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    ppml_boot_model *m = g_pboot;
    const hdfe_wild_score_model *sm = &m->wb_score;
    if (!m->wb_test || sm->B < 1)
        return 198;
    if (sm->got != sm->B)
    {
        SF_error("cppmlhdfe: wild_finish before every replicate's weights arrived\n");
        return 198;
    }
    int ptype = 0;
    if (ctools_parse_int_option(args, "ptype", &ptype) != CTOOLS_PARSE_OK || ptype < 0 ||
        ptype > 3)
        ptype = 0;
    const double level = ctools_parse_double_option(args, "level", 95.0);
    const int ci = ctools_parse_bool_option(args, "ci");
    const ST_int N = m->N;
    const ST_int all = ctools_get_max_threads() > 1 ? ctools_get_max_threads() : 1;

    /* Evaluators: one for the tests' fits (all threads), P for the
     * interval search's points, side by side within the memory budget,
     * and each chain's last means when they fit too */
    const size_t budget = boot_budget();
    const size_t eval = reval_bytes(m, all, sm->B);
    const int n_chain = PPML_WILD_NPT + 2;
    const size_t chain_bytes = (size_t)n_chain * (size_t)N * sizeof(ST_double);
    ST_int P = 1;
    int chains = 0;
    if (ci)
    {
        chains = budget >= chain_bytes + eval;
        const size_t room = budget - (chains ? chain_bytes : 0);
        P = all < PPML_WILD_NPT ? all : PPML_WILD_NPT;
        if (eval > 0 && room / eval < (size_t)P)
            P = (ST_int)(room / eval);
        if (P < 1)
            P = 1;
    }
    const ST_int inner = all / P > 1 ? all / P : 1;
    ppml_reval **ev = calloc((size_t)P, sizeof(*ev));
    ST_double *mu0 = malloc((size_t)N * sizeof(ST_double));
    ST_double *chain_block = chains ? malloc(chain_bytes) : NULL;
    ST_double *chain[PPML_WILD_NPT + 2] = {0};
    if (!ev || !mu0 || (chains && !chain_block))
        rc = 920;
    for (ST_int k = 0; k < P && !rc; k++)
    {
        if (!(ev[k] = reval_new(m, k == 0 ? all : inner, sm->B)))
        {
            /* fewer points side by side; not even one fit is an error */
            if (k == 0)
                rc = 920;
            P = k;
        }
    }
    for (ST_int i = 0; i < N && !rc; i++)
        mu0[i] = exp(m->eta[i] + (m->offset ? m->offset[i] : 0.0));
    for (int q = 0; q < n_chain && chain_block && !rc; q++)
    {
        chain[q] = chain_block + (size_t)q * (size_t)N;
        memcpy(chain[q], mu0, (size_t)N * sizeof(ST_double));
    }

    ST_int min_feasible = sm->B;
    for (ST_int j = 0; j < m->wb_ntest && !rc; j++)
    {
        const ST_int col = m->wb_test[j];
        memcpy(ev[0]->mu, mu0, (size_t)N * sizeof(ST_double));
        int ok = 0;
        rc = restricted_fit(m, ev[0], col, 0.0, all, &ok);
        if (rc)
            break;
        double p = NAN;
        if (ok)
        {
            if ((rc = hdfe_wild_score_t(sm, ev[0]->e, all, ev[0]->t)))
                break;
            ST_int F = 0;
            p = hdfe_wild_pvalue_t(ev[0]->t, sm->B, ptype, &F);
            if (F < min_feasible)
                min_feasible = F;
        }
        rc = ctools_mat_store("__cppmlhdfe_wb_p", 1, j + 1, stata_value(p));
        /* the observed score statistic, which the replicates' are compared to */
        if (!rc)
            rc = ctools_mat_store("__cppmlhdfe_wb_score", 1, j + 1,
                                  ok ? stata_value(ev[0]->t[0]) : SV_missval);
        if (!rc)
            rc = ctools_mat_store("__cppmlhdfe_wb_fail", 1, j + 1, ok ? 0.0 : 1.0);
        if (rc || !ci || !ok)
            continue;
        /* the coefficient in the data's units */
        ST_int kept = 0;
        for (ST_int c = 0; c < col; c++)
            kept += !m->is_collinear[c];
        for (int q = 0; q < n_chain && chain_block; q++)
            memcpy(chain[q], mu0, (size_t)N * sizeof(ST_double));
        wild_pctx c = {m, col, ptype, ev, P, inner, mu0, chain_block ? chain : NULL, n_chain, 0};
        hdfe_wild_ci set;
        hdfe_wild_invert_batch(wild_p_batch, &c, m->b_hat[kept], m->wb_scale[j],
                               1.0 - level / 100.0, PPML_WILD_HALF, 1e-4, &set);
        if ((rc = c.rc))
            break;
        const double row[7] = {set.lo,       set.hi,       set.open_lo,         set.open_hi,
                               set.disjoint, set.empty,    set.failed > 0};
        for (int k = 0; k < 2 && !rc; k++)
            rc = ctools_mat_store("__cppmlhdfe_wb_ci", j + 1, k + 1, stata_value(row[k]));
        for (int k = 0; k < 5 && !rc; k++)
            rc = ctools_mat_store("__cppmlhdfe_wb_cistat", j + 1, k + 1, row[2 + k]);
    }
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_wb_reps", (ST_double)sm->B);
    if (!rc)
        rc = ctools_scal_save("__cppmlhdfe_wb_misreps", (ST_double)(sm->B - min_feasible));
    for (ST_int k = 0; ev && k < P; k++)
        reval_free(ev[k]);
    free(ev);
    free(mu0);
    free(chain_block);
    return rc;
}
