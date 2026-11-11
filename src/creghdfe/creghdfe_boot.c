/*
 * creghdfe_boot.c
 *
 * vce(bootstrap) for creghdfe (see creghdfe_boot.h). A replicate is the
 * full-sample model re-fit with each row weighted by its draw multiplicity:
 * singletons judged on the multiplicities, the partialled full-sample
 * columns projected again under the replicate weights (they give the same
 * result as the raw columns and start closer to it), then creghdfe's own
 * collinearity rules and OLS (creghdfe_core.c).
 * Part of the ctools Stata plugin suite
 */

/* Checked SPI access: indices resolve in the plugin varlist (as in
 * creghdfe_regress.c) */
#define SD_SAFEMODE
#include "creghdfe_boot.h"

#include <limits.h>

#include "creghdfe_core.h"
#include "creghdfe_hdfe.h"
#include "hdfe/hdfe_boot.h"
#include "hdfe/hdfe_boot_stata.h"
#include "hdfe/hdfe_fe_recover.h"
#include "hdfe/hdfe_sample.h"
#include "hdfe/hdfe_solver.h"
#include "hdfe/hdfe_wild.h"
#include "../ctools_config.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_parse.h"
#include "../ctools_platform.h"
#include "../ctools_spi.h"

typedef struct
{
    ST_int nobs; /* SF_nobs() when retained: the data may not change between calls */
    ST_int N, K, G;
    ST_double *raw, *tilde, *stdevs, *weights;
    ST_int weight_type;
    ST_int **levels; /* G x N, 1-based full-sample levels */
    ST_int *num_levels;
    ST_int *is_collinear;
    ST_double *b_hat;
    ST_int K_keep;
    ST_int maxiter, drop_singletons, use_quad;
    ST_double tolerance, collinear_tol;
    /* Resampling frame */
    ST_int *unit_of_row, *stratum_start;
    ST_double *unit_freq;
    hdfe_boot_design design;
    int cluster_mode;
    hdfe_boot_acc acc;
    /* Heterogeneous slopes: the fit's slope factors (slope_idx, fe_absvar,
     * num_slopes of slope_in) and each slope factor's variable on the
     * compact rows (NULL for an intercept); slope NULL without slopes */
    creghdfe_input slope_in;
    ST_double **slope;
    ST_int is_intercept[10]; /* factor g is an intercept */
    ST_int n_intercepts;
    /* Slope factor g: the intercept factor of its absorbed term, or -1;
     * slope_est (L_total, slope nodes only): the full fit estimates the
     * level's slope (its variable varies within the level) */
    ST_int slope_intercept[10];
    unsigned char *slope_est;
    /* fevce(bootstrap): the saved fixed effects (plan 5.7) */
    ST_int *obs_of_row;   /* N: 1-based Stata observation of each compact row */
    ST_int *offset;       /* G + 1: first node of each effect among L_total levels */
    ST_int L_total;
    ST_int *joint_comp;   /* L_total: full-sample connected set of each level, all effects */
    ST_int n_joint;
    ST_int **pair_comp;   /* G: for g >= 1, the set of each level of g - 1 and g */
    ST_int *n_pair;
    int fe_target, dse_target;
    /* wildbootstrap() (plan 7.5): the partialled kept regressors in the
     * units of the data, residuals, A = (X'WX)^-1, cluster ids */
    ST_int *wb_vce, *wb_boot; /* N dense ids, or NULL: each row */
    ST_int wb_G, wb_C;
    ST_int *wb_test;          /* tested kept columns */
    ST_int wb_ntest;
    int wb_null, wb_projected;
    ST_double *wb_X, *wb_u, *wb_A, *wb_beta;
    hdfe_wild_model wild;
    HDFE_State wb_S;          /* the projected path's factors */
    int wb_S_ready;
} boot_model;

static boot_model *g_boot = NULL;

static void wild_free(boot_model *m);
static ST_retcode save_parallel(const boot_model *m);

static void model_free(boot_model *m)
{
    if (!m)
        return;
    wild_free(m);
    free(m->raw);
    free(m->tilde);
    free(m->stdevs);
    free(m->weights);
    for (ST_int g = 0; m->levels && g < m->G; g++)
        free(m->levels[g]);
    free(m->levels);
    free(m->num_levels);
    free(m->is_collinear);
    free(m->b_hat);
    free(m->unit_of_row);
    free(m->stratum_start);
    free(m->unit_freq);
    for (ST_int g = 0; m->slope && g < m->G; g++)
        free(m->slope[g]);
    free(m->slope);
    free(m->slope_est);
    hdfe_boot_acc_free(&m->acc);
    free(m->obs_of_row);
    free(m->offset);
    free(m->joint_comp);
    for (ST_int g = 0; m->pair_comp && g < m->G; g++)
        free(m->pair_comp[g]);
    free(m->pair_comp);
    free(m->n_pair);
    free(m);
}

void creghdfe_boot_release(void)
{
    model_free(g_boot);
    g_boot = NULL;
}

void creghdfe_boot_install(void)
{
    creghdfe_boot_retain_hook = creghdfe_boot_retain;
    creghdfe_boot_release_hook = creghdfe_boot_release;
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
 * ascending value, or rows in data order (the compact rows follow the
 * Stata observation order). */
static ST_retcode build_frame(boot_model *m, const creghdfe_boot_source *src)
{
    ST_int N = m->N;
    unit_key *keys = malloc((size_t)N * sizeof(*keys));
    m->unit_of_row = malloc((size_t)N * sizeof(ST_int));
    if (!keys || !m->unit_of_row)
    {
        free(keys);
        return 920;
    }
    ST_int idx = 0;
    for (ST_int row = 0; row < src->N_orig; row++)
    {
        if (!src->mask[row])
            continue;
        double s = 0.0, c = 0.0;
        ST_retcode rc = 0;
        if (src->strata_var > 0)
            rc = SF_vdata(src->strata_var, (ST_int)src->obs_map[row], &s);
        if (!rc && src->unit_var > 0)
            rc = SF_vdata(src->unit_var, (ST_int)src->obs_map[row], &c);
        if (rc)
        {
            free(keys);
            return rc;
        }
        keys[idx].stratum = s;
        keys[idx].cluster = src->unit_var > 0 ? c : 0.0;
        keys[idx].row = idx;
        idx++;
    }
    /* A cluster that spans strata is a separate unit in each stratum, as
     * bsample numbers clusters within strata */
    m->cluster_mode = src->unit_var > 0;
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
     * (observation-level resampling only) */
    if (m->weights && m->weight_type == 2 && !m->cluster_mode)
    {
        m->unit_freq = malloc((size_t)U * sizeof(ST_double));
        if (!m->unit_freq)
        {
            free(keys);
            return 920;
        }
        for (ST_int i = 0; i < N; i++)
            m->unit_freq[m->unit_of_row[i]] = m->weights[i];
    }
    free(keys);
    m->design.U = U;
    m->design.S = S;
    m->design.stratum_start = m->stratum_start;
    m->design.unit_freq = m->unit_freq;
    m->design.seed = 0;
    return hdfe_boot_design_check(&m->design);
}

/* Connected sets of levels: for all effects together (joint) and for each
 * consecutive pair g - 1, g, over the rows kept[0..n), levels in m->levels
 * (full-sample numbering). Nodes: effect g's levels at offset[g]. With
 * comp != NULL, writes each node's dense set number; with first != NULL,
 * instead marks in split[] the joint sets whose present nodes fall in more
 * than one set here (first: per-set scratch, -1 filled by the caller). */
static ST_retcode level_sets(const boot_model *m, const ST_int *kept, ST_int n,
                             const ST_int *facs, ST_int nf, const ST_int *node_off, ST_int *comp,
                             ST_int *n_comp)
{
    /* Nodes: the listed factors' levels, one after another (comp in that
     * order) */
    ST_int off[11];
    off[0] = 0;
    for (ST_int j = 0; j < nf; j++)
        off[j + 1] = off[j] + m->num_levels[facs[j]];
    (void)node_off;
    ST_int nodes = off[nf];
    ctools_UnionFind *uf = ctools_uf_create(nodes);
    if (!uf)
        return 920;
    for (ST_int r = 0; r < n; r++)
    {
        ST_int i = kept[r];
        ST_int a = m->levels[facs[0]][i] - 1;
        for (ST_int j = 1; j < nf; j++)
            ctools_uf_union(uf, a, off[j] + m->levels[facs[j]][i] - 1);
    }
    ST_int *dense = malloc((size_t)nodes * sizeof(ST_int));
    if (!dense)
    {
        ctools_uf_destroy(uf);
        return 920;
    }
    for (ST_int v = 0; v < nodes; v++)
        dense[v] = -1;
    ST_int k = 0;
    for (ST_int v = 0; v < nodes; v++)
    {
        ST_int root = ctools_uf_find(uf, v);
        if (dense[root] < 0)
            dense[root] = k++;
        comp[v] = dense[root];
    }
    *n_comp = k;
    free(dense);
    ctools_uf_destroy(uf);
    return 0;
}

/* The intercept factors, in order (connected sets and the normalization
 * concern the intercepts: a slope is identified within its level) */
static ST_int intercept_list(const boot_model *m, ST_int *facs)
{
    ST_int n = 0;
    for (ST_int g = 0; g < m->G; g++)
    {
        if (m->is_intercept[g])
            facs[n++] = g;
    }
    return n;
}

static ST_retcode full_sample_sets(boot_model *m)
{
    const ST_int G = m->G, N = m->N;
    m->offset = malloc((size_t)(G + 1) * sizeof(ST_int));
    if (!m->offset)
        return 920;
    m->offset[0] = 0;
    for (ST_int g = 0; g < G; g++)
        m->offset[g + 1] = m->offset[g] + m->num_levels[g];
    m->L_total = m->offset[G];
    ST_int facs[10];
    const ST_int nf = intercept_list(m, facs);
    m->n_intercepts = nf;
    if (nf < 2)
        return 0;
    ST_int *all = malloc((size_t)N * sizeof(ST_int));
    ST_int *comp = malloc((size_t)m->L_total * sizeof(ST_int));
    m->joint_comp = malloc((size_t)m->L_total * sizeof(ST_int));
    m->pair_comp = calloc((size_t)G, sizeof(ST_int *));
    m->n_pair = calloc((size_t)G, sizeof(ST_int));
    if (!all || !comp || !m->joint_comp || !m->pair_comp || !m->n_pair)
    {
        free(all);
        free(comp);
        return 920;
    }
    for (ST_int i = 0; i < N; i++)
        all[i] = i;
    /* joint_comp by node; slope nodes are in no set (-1) */
    for (ST_int v = 0; v < m->L_total; v++)
        m->joint_comp[v] = -1;
    ST_retcode rc = level_sets(m, all, N, facs, nf, m->offset, comp, &m->n_joint);
    for (ST_int j = 0, q = 0; j < nf && !rc; j++)
    {
        for (ST_int l = 0; l < m->num_levels[facs[j]]; l++)
            m->joint_comp[m->offset[facs[j]] + l] = comp[q++];
    }
    /* pairs of consecutive intercepts p - 1, p (in their own numbering) */
    for (ST_int p = 1; p < nf && !rc; p++)
    {
        m->pair_comp[p] = malloc((size_t)(m->num_levels[facs[p - 1]] + m->num_levels[facs[p]]) *
                                 sizeof(ST_int));
        if (!m->pair_comp[p])
            rc = 920;
        else
            rc = level_sets(m, all, N, facs + p - 1, 2, m->offset, m->pair_comp[p],
                            &m->n_pair[p]);
    }
    free(all);
    free(comp);
    return rc;
}

/* Marks split[j] for every full-sample joint set j that the replicate's rows
 * (kept[0..n)) leave in more than one piece, for all effects together or for
 * any consecutive pair (plan 5.7). present[node] flags the drawn levels. */
static ST_retcode replicate_splits(const boot_model *m, const ST_int *kept, ST_int n,
                                   const unsigned char *present, unsigned char *split)
{
    ST_int facs[10];
    const ST_int nf = intercept_list(m, facs);
    ST_retcode rc = 0;
    ST_int *comp = malloc((size_t)m->L_total * sizeof(ST_int));
    ST_int *first = malloc((size_t)m->L_total * sizeof(ST_int));
    ST_int *full_joint = malloc((size_t)m->L_total * sizeof(ST_int));
    if (!comp || !first || !full_joint)
        rc = 920;
    /* the full-sample joint sets in the listed intercepts' node order */
    for (ST_int j = 0, q = 0; j < nf && !rc; j++)
    {
        for (ST_int l = 0; l < m->num_levels[facs[j]]; l++)
            full_joint[q++] = m->joint_comp[m->offset[facs[j]] + l];
    }
    for (ST_int p = 0; p < nf && !rc; p++)
    {
        /* pass 0: all intercepts; then each pair p - 1, p */
        const ST_int *list = p == 0 ? facs : facs + p - 1;
        const ST_int nl = p == 0 ? nf : 2;
        ST_int k = 0;
        if ((rc = level_sets(m, kept, n, list, nl, m->offset, comp, &k)))
            break;
        const ST_int *full = p == 0 ? full_joint : m->pair_comp[p];
        ST_int n_full = p == 0 ? m->n_joint : m->n_pair[p];
        for (ST_int c = 0; c < n_full; c++)
            first[c] = -1;
        for (ST_int j = 0, v = 0; j < nl; j++)
        {
            for (ST_int l = 0; l < m->num_levels[list[j]]; l++, v++)
            {
                ST_int node = m->offset[list[j]] + l;
                if (!present[node])
                    continue;
                if (first[full[v]] < 0)
                    first[full[v]] = comp[v];
                else if (first[full[v]] != comp[v])
                    split[m->joint_comp[node]] = 1;
            }
        }
    }
    free(comp);
    free(first);
    free(full_joint);
    return rc;
}

static ST_double *copy_doubles(const ST_double *src, size_t n)
{
    ST_double *dst = malloc((n ? n : 1) * sizeof(*dst));
    if (dst && n)
        memcpy(dst, src, n * sizeof(*dst));
    return dst;
}

ST_retcode creghdfe_boot_retain(const creghdfe_boot_source *src)
{
    creghdfe_boot_release();
    boot_model *m = calloc(1, sizeof(*m));
    if (!m)
        return 920;
    const size_t NK = (size_t)src->N * (size_t)src->K;
    const ST_int K_x = src->K - 1;
    m->nobs = SF_nobs();
    m->N = src->N;
    m->K = src->K;
    m->G = src->G;
    m->weight_type = src->weight_type;
    m->K_keep = src->K_keep;
    m->maxiter = src->maxiter;
    m->drop_singletons = src->drop_singletons;
    m->use_quad = src->use_quad;
    m->tolerance = src->tolerance;
    m->collinear_tol = src->collinear_tol;
    m->raw = copy_doubles(src->raw, NK);
    m->tilde = copy_doubles(src->tilde, NK);
    m->stdevs = copy_doubles(src->stdevs, (size_t)src->K);
    m->weights = src->weights ? copy_doubles(src->weights, (size_t)src->N) : NULL;
    m->b_hat = copy_doubles(src->b_hat, (size_t)src->K_keep + 1);
    m->is_collinear = malloc((size_t)(K_x > 0 ? K_x : 1) * sizeof(ST_int));
    m->levels = calloc((size_t)src->G, sizeof(ST_int *));
    m->num_levels = malloc((size_t)src->G * sizeof(ST_int));
    ST_retcode rc = 0;
    if (!m->raw || !m->tilde || !m->stdevs || (src->weights && !m->weights) || !m->b_hat ||
        !m->is_collinear || !m->levels || !m->num_levels)
    {
        rc = 920;
        goto fail;
    }
    if (K_x > 0)
        memcpy(m->is_collinear, src->is_collinear, (size_t)K_x * sizeof(ST_int));
    for (ST_int g = 0; g < src->G; g++)
    {
        m->num_levels[g] = src->factors[g].num_levels;
        m->levels[g] = malloc((size_t)src->N * sizeof(ST_int));
        if (!m->levels[g])
        {
            rc = 920;
            goto fail;
        }
        memcpy(m->levels[g], src->factors[g].levels, (size_t)src->N * sizeof(ST_int));
    }
    for (ST_int g = 0; g < src->G; g++)
        m->is_intercept[g] = 1;
    if (src->in && src->in->num_slopes > 0)
    {
        m->slope_in = *src->in;
        if (!(m->slope = calloc((size_t)src->G, sizeof(ST_double *))))
        {
            rc = 920;
            goto fail;
        }
        for (ST_int g = 0; g < src->G; g++)
        {
            if (!src->in->slope_idx[g])
                continue;
            m->is_intercept[g] = 0;
            if (!src->factors[g].slope ||
                !(m->slope[g] = copy_doubles(src->factors[g].slope, (size_t)src->N)))
            {
                rc = src->factors[g].slope ? 920 : 198;
                goto fail;
            }
        }
    }
    if ((rc = build_frame(m, src)))
        goto fail;
    m->obs_of_row = malloc((size_t)src->N * sizeof(ST_int));
    if (!m->obs_of_row)
    {
        rc = 920;
        goto fail;
    }
    for (ST_int row = 0, i = 0; row < src->N_orig; row++)
    {
        if (src->mask[row])
            m->obs_of_row[i++] = (ST_int)src->obs_map[row];
    }
    if ((rc = full_sample_sets(m)))
        goto fail;
    if (m->slope)
    {
        /* The levels whose slope the full fit estimates (creghdfe_setup_slopes
         * zeroes counts where reghdfe drops the slope) */
        if (!(m->slope_est = calloc((size_t)(m->L_total > 0 ? m->L_total : 1), 1)))
        {
            rc = 920;
            goto fail;
        }
        for (ST_int g = 0; g < src->G; g++)
        {
            m->slope_intercept[g] = -1;
            if (m->is_intercept[g])
                continue;
            for (ST_int h = 0; h < src->G; h++)
            {
                if (m->is_intercept[h] && src->in->fe_absvar[h] == src->in->fe_absvar[g])
                {
                    m->slope_intercept[g] = h;
                    break;
                }
            }
            for (ST_int l = 0; l < m->num_levels[g]; l++)
                m->slope_est[m->offset[g] + l] = src->factors[g].counts[l] > 0.0;
        }
    }
    if ((rc = ctools_scal_save("__creghdfe_bs_U", hdfe_boot_draw_total(&m->design))))
        goto fail;
    if ((rc = ctools_scal_save("__creghdfe_bs_S", (ST_double)m->design.S)))
        goto fail;
    if ((rc = hdfe_boot_save_dup_strata(&m->design, "creghdfe")))
        goto fail;
    if ((rc = save_parallel(m)))
        goto fail;
    g_boot = m;
    return 0;
fail:
    model_free(m);
    return rc;
}

/* ========================================================================
 * One replicate
 * ======================================================================== */

static size_t fit_bytes(const void *mv)
{
    const boot_model *m = mv;
    /* columns, raw norms, weights, counts, levels, CG buffers */
    return (size_t)m->N * ((size_t)m->K * 2 + 10) * sizeof(ST_double) +
           (size_t)m->N * (size_t)m->G * 2 * sizeof(ST_int) +
           (size_t)m->L_total * (2 * sizeof(ST_double) + 3 * sizeof(ST_int));
}

/* Each further thread of a replicate's projection: three CG buffers and
 * the per-level means */
static size_t fit_thread_bytes(const void *mv)
{
    const boot_model *m = mv;
    return (size_t)m->N * 3 * sizeof(ST_double) + (size_t)m->L_total * sizeof(ST_double);
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

/* The replicate's fixed effects, as creghdfe's savefe recovers them, in
 * full-sample level order; levels not drawn, or in a connected set that the
 * replicate splits, stay missing. Also d per compact row. */
static ST_retcode replicate_effects(const boot_model *m, HDFE_State *S, const ST_int *rows,
                                    const ST_int *mask, ST_int N1, ST_int N2,
                                    const creghdfe_ols_out *ols, hdfe_boot_rep *out)
{
    const ST_int G = m->G, N = m->N;
    ST_retcode rc = 0;
    ST_int *kept = malloc((size_t)N2 * sizeof(ST_int));
    ST_double *error = malloc((size_t)N2 * sizeof(ST_double));
    ST_double **alpha = calloc((size_t)G, sizeof(ST_double *));
    unsigned char *present = calloc((size_t)m->L_total, 1);
    unsigned char *split = calloc((size_t)(m->n_joint > 0 ? m->n_joint : 1), 1);
    if (!kept || !error || !alpha || !present || !split)
    {
        rc = 920;
        goto done;
    }
    for (ST_int r = 0, r2 = 0; r < N1; r++)
    {
        if (mask[r])
            kept[r2++] = rows[r];
    }
    /* d = y - _cons - Xb - residual on the kept rows (savefe's formula) */
    const ST_int Kk = ols->K_keep;
    for (ST_int r2 = 0; r2 < N2; r2++)
    {
        ST_int i = kept[r2];
        ST_double fit_s = 0.0;
        for (ST_int k = 0; k < Kk; k++)
            fit_s += ols->data_keep[(size_t)(k + 1) * N2 + r2] * ols->beta_keep[k];
        ST_double resid = (ols->data_keep[r2] - fit_s) * m->stdevs[0];
        ST_double e = m->raw[i] - out->b[Kk] - resid;
        for (ST_int k = 0; k < Kk; k++)
            e -= out->b[k] * m->raw[(size_t)ols->keep_idx[k] * N + i];
        error[r2] = e;
    }
    for (ST_int g = 0; g < G; g++)
    {
        if (!(alpha[g] = calloc((size_t)S->factors[g].num_levels, sizeof(ST_double))))
        {
            rc = 920;
            goto done;
        }
    }
    rc = m->slope ? creghdfe_fe_recover_slopes(S, &m->slope_in, error, alpha)
                  : hdfe_fe_recover(S, error, alpha);
    if (rc == 430)
    {
        rc = 0; /* not converged: this replicate's effects stay missing */
        goto done;
    }
    if (rc)
        goto done;
    for (ST_int r2 = 0; r2 < N2; r2++)
    {
        for (ST_int g = 0; g < G; g++)
            present[m->offset[g] + m->levels[g][kept[r2]] - 1] = 1;
    }
    if (m->n_intercepts >= 2)
    {
        if ((rc = replicate_splits(m, kept, N2, present, split)))
            goto done;
        for (ST_int j = 0; j < m->n_joint; j++)
            out->fe_split |= split[j];
    }
    /* Replicate level q of effect g is the full-sample level of its rows */
    for (ST_int r2 = 0; r2 < N2; r2++)
    {
        for (ST_int g = 0; g < G; g++)
        {
            ST_int node = m->offset[g] + m->levels[g][kept[r2]] - 1;
            /* a slope is identified within its level (and its coefficient
             * is zero where reghdfe drops it) */
            if (m->n_intercepts >= 2 && m->joint_comp[node] >= 0 &&
                split[m->joint_comp[node]])
                continue;
            if (!m->is_intercept[g] && !(S->factors[g].counts[S->factors[g].levels[r2] - 1] > 0.0))
                continue;
            out->fe[node] = alpha[g][S->factors[g].levels[r2] - 1];
        }
    }
    /* A level whose slope the full fit estimates but this replicate cannot
     * (its drawn rows share one value of the slope variable) identifies only
     * intercept + slope x value: the recovery sets the slope to 0, and the
     * intercept is then not the full fit's quantity either */
    for (ST_int r2 = 0; m->slope_est && r2 < N2; r2++)
    {
        for (ST_int g = 0; g < G; g++)
        {
            const ST_int h = m->is_intercept[g] ? -1 : m->slope_intercept[g];
            if (h < 0 || !m->slope_est[m->offset[g] + m->levels[g][kept[r2]] - 1] ||
                S->factors[g].counts[S->factors[g].levels[r2] - 1] > 0.0)
                continue;
            out->fe[m->offset[h] + m->levels[h][kept[r2]] - 1] = NAN;
        }
    }
    /* d of every compact row whose levels are all estimable: intercepts,
     * plus slopes times their variable */
    for (ST_int i = 0; out->d && i < N; i++)
    {
        ST_double d = 0.0;
        for (ST_int g = 0; g < G; g++)
        {
            ST_double a = out->fe[m->offset[g] + m->levels[g][i] - 1];
            d += m->is_intercept[g] ? a : a * m->slope[g][i];
        }
        out->d[i] = d;
    }
done:
    for (ST_int g = 0; alpha && g < G; g++)
        free(alpha[g]);
    free(alpha);
    free(kept);
    free(error);
    free(present);
    free(split);
    return rc;
}

static ST_retcode fit_replicate(const void *mv, void *ws, const ST_int *unit_mult,
                                hdfe_boot_rep *out)
{
    const boot_model *m = mv;
    const ST_int N = m->N, K = m->K, G = m->G, K_x = K - 1;
    (void)ws;
    ST_retcode rc = 0;
    HDFE_State S;
    memset(&S, 0, sizeof(S));
    creghdfe_ols_out ols;
    memset(&ols, 0, sizeof(ols));

    /* Rows drawn at least once, and their count weights for the singleton
     * rule: a level whose rows total one copy is a singleton, as in the
     * expanded data (with fweights, copies times fweight; when the fweighted
     * observations themselves are drawn, unit_freq, the count is the number
     * of them drawn, and so is the weight) */
    ST_int *rows = malloc((size_t)N * sizeof(ST_int));
    ST_double *count = malloc((size_t)N * sizeof(ST_double));
    ST_int **lv = calloc((size_t)G, sizeof(ST_int *));
    ST_int *mask = malloc((size_t)N * sizeof(ST_int));
    ST_double *omega = malloc((size_t)N * sizeof(ST_double));
    ST_double *data = NULL, *norm2 = NULL, *count_kept = NULL, *means = NULL;
    ST_double *xtx = NULL, *xty = NULL, *scratch = NULL;
    ST_int *fe_flagged = NULL, *is_coll = NULL;
    /* Fixed effects that this replicate cannot estimate stay missing */
    for (ST_int v = 0; out->fe && v < m->L_total; v++)
        out->fe[v] = NAN;
    for (ST_int i = 0; out->d && i < N; i++)
        out->d[i] = NAN;
    if (!rows || !count || !lv || !mask || !omega)
    {
        rc = 920;
        goto done;
    }
    ST_int N1 = 0;
    for (ST_int i = 0; i < N; i++)
    {
        ST_int mult = unit_mult[m->unit_of_row[i]];
        if (mult <= 0)
            continue;
        rows[N1] = i;
        count[N1] = (double)mult *
                    ((m->weights && m->weight_type == 2 && !m->unit_freq) ? m->weights[i] : 1.0);
        N1++;
    }
    if (N1 == 0)
    {
        out->status = HDFE_BOOT_EMPTY;
        goto done;
    }
    for (ST_int g = 0; g < G; g++)
    {
        if (!(lv[g] = malloc((size_t)N1 * sizeof(ST_int))))
        {
            rc = 920;
            goto done;
        }
        for (ST_int r = 0; r < N1; r++)
            lv[g][r] = m->levels[g][rows[r]];
    }
    ST_int N2 = N1;
    if (m->drop_singletons)
    {
        if (ctools_remove_singletons(lv, G, N1, mask, count, 0) < 0)
        {
            rc = 920;
            goto done;
        }
        N2 = 0;
        for (ST_int r = 0; r < N1; r++)
            N2 += mask[r] != 0;
    }
    else
    {
        for (ST_int r = 0; r < N1; r++)
            mask[r] = 1;
    }
    if (N2 == 0)
    {
        out->status = HDFE_BOOT_EMPTY;
        goto done;
    }

    /* Replicate weights: multiplicity times the fit's weight (the drawn
     * observations alone when the fweighted observations are drawn) */
    for (ST_int r = 0; r < N1; r++)
    {
        ST_int i = rows[r];
        omega[r] = (double)unit_mult[m->unit_of_row[i]] *
                   ((m->weights && !m->unit_freq) ? m->weights[i] : 1.0);
    }
    ST_int level_max[10];
    for (ST_int g = 0; g < G; g++)
        level_max[g] = m->num_levels[g];
    if ((rc = hdfe_build_factors(&S, G, lv, level_max, N1, mask, N2, omega, 0)))
        goto done;
    S.has_weights = 1;
    S.weight_type = 2;
    S.maxiter = m->maxiter;
    S.tolerance = m->tolerance;
    /* Column parallelism with the threads the engine left over; each
     * column is projected by one thread, so the result does not change */
    S.num_threads = out->threads > 1 ? out->threads : 1;
    S.factors_initialized = 1;
    if (m->slope)
    {
        /* Slope factors on the drawn rows, under the replicate weights */
        const ST_double *cols[10] = {0};
        ST_double *sbuf = malloc((size_t)G * (size_t)N1 * sizeof(ST_double));
        if (!sbuf)
        {
            rc = 920;
            goto done;
        }
        for (ST_int g = 0; g < G; g++)
        {
            if (!m->slope[g])
                continue;
            ST_double *col = sbuf + (size_t)g * N1;
            for (ST_int r = 0; r < N1; r++)
                col[r] = m->slope[g][rows[r]];
            cols[g] = col;
        }
        rc = creghdfe_setup_slopes(&S, &m->slope_in, cols, mask, N1, omega);
        free(sbuf);
        if (rc)
            goto done;
    }
    S.weights = malloc((size_t)N2 * sizeof(ST_double));
    data = malloc((size_t)N2 * (size_t)K * sizeof(ST_double));
    norm2 = calloc((size_t)K, sizeof(ST_double));
    count_kept = malloc((size_t)N2 * sizeof(ST_double));
    means = calloc((size_t)K, sizeof(ST_double));
    size_t KX = (size_t)(K_x > 0 ? K_x : 1);
    xtx = malloc(KX * KX * sizeof(ST_double));
    xty = malloc(KX * sizeof(ST_double));
    scratch = malloc(KX * KX * sizeof(ST_double));
    fe_flagged = malloc(KX * sizeof(ST_int));
    is_coll = malloc(KX * sizeof(ST_int));
    if (!S.weights || !data || !norm2 || !count_kept || !means || !xtx || !xty || !scratch ||
        !fe_flagged || !is_coll)
    {
        rc = 920;
        goto done;
    }
    /* The kept rows: the partialled full-sample columns (projecting them
     * again under the replicate weights gives what the raw columns would),
     * the raw norms and the weighted means of the raw columns */
    ST_double sum_omega = 0.0;
    ST_int r2 = 0;
    for (ST_int r = 0; r < N1; r++)
    {
        if (!mask[r])
            continue;
        ST_int i = rows[r];
        S.weights[r2] = omega[r];
        count_kept[r2] = (double)unit_mult[m->unit_of_row[i]];
        sum_omega += omega[r];
        for (ST_int k = 0; k < K; k++)
        {
            double x = m->raw[(size_t)k * N + i];
            double xs = x / m->stdevs[k];
            data[(size_t)k * N2 + r2] = m->tilde[(size_t)k * N + i];
            norm2[k] += count_kept[r2] * xs * xs;
            means[k] += omega[r] * x;
        }
        r2++;
    }
    for (ST_int k = 0; k < K; k++)
        means[k] /= sum_omega;
    S.sum_weights = sum_omega;
    if (ctools_hdfe_alloc_buffers(&S, 0, K) != 0)
    {
        rc = 920;
        goto done;
    }
    HDFE_SolveResult proj = partial_out_columns(&S, data, N2, K, 1);
    out->iterations = proj.iterations;
    if (proj.status == 430)
    {
        out->status = HDFE_BOOT_NOCONV;
        goto done;
    }
    if (proj.status)
    {
        rc = proj.status;
        goto done;
    }

    ST_int depvar_explained = 0;
    if ((rc = creghdfe_collinear(data, N2, K, norm2, count_kept, S.weights, 2, m->collinear_tol,
                                 fe_flagged, is_coll, &depvar_explained, xtx, xty, scratch)))
        goto done;
    for (ST_int k = 0; k < K_x; k++)
    {
        if (is_coll[k] != m->is_collinear[k])
        {
            /* As the bootstrap prefix: a replicate whose omitted
             * coefficients differ from the full sample's is missing */
            out->status = HDFE_BOOT_COLLINEAR;
            goto done;
        }
    }

    {
    creghdfe_ols_in in = {
        .data = data, .N = N2, .K = K, .is_collinear = is_coll, .xtx = xtx, .xty_all = xty,
        .weights = S.weights, .weight_type = 2, .sum_weights = sum_omega,
        .means = means, .stdevs = m->stdevs, .standardize = 1, .use_quad = m->use_quad,
    };
    ST_retcode ols_rc = creghdfe_ols(&in, &ols);
    if (ols_rc == 198)
    {
        out->status = HDFE_BOOT_COLLINEAR;
        goto done;
    }
    if (ols_rc)
    {
        rc = ols_rc;
        goto done;
    }
    }
    /* Back to the units of the data: b_k * sd(y) / sd(x_k); the constant
     * from the replicate's weighted means. Without an absorbed intercept
     * (slopes alone) the model has no constant: the full fit posts 0, and
     * savefe's recovery, which has no level mean to absorb one, must not
     * subtract one either. */
    ST_double xb = 0.0;
    for (ST_int k = 0; k < ols.K_keep; k++)
    {
        ST_int col = ols.keep_idx[k];
        out->b[k] = ols.beta_keep[k] * m->stdevs[0] / m->stdevs[col];
        xb += means[col] * out->b[k];
    }
    out->b[ols.K_keep] = m->n_intercepts > 0 ? means[0] - xb : 0.0;
    out->status = HDFE_BOOT_OK;
    if (out->fe)
    {
        rc = replicate_effects(m, &S, rows, mask, N1, N2, &ols, out);
        if (rc)
            goto done;
    }

done:
    if (ols.data_keep != data)
        free(ols.data_keep);
    free(ols.keep_idx);
    free(ols.xtx_keep);
    free(ols.xty_keep);
    free(ols.beta_keep);
    free(ols.inv_xx_keep);
    free(ols.means_x);
    free(S.weights);
    S.weights = NULL;
    if (S.factors)
        ctools_hdfe_state_cleanup(&S);
    for (ST_int g = 0; lv && g < G; g++)
        free(lv[g]);
    free(lv);
    free(rows);
    free(count);
    free(mask);
    free(omega);
    free(data);
    free(norm2);
    free(count_kept);
    free(means);
    free(xtx);
    free(xty);
    free(scratch);
    free(fe_flagged);
    free(is_coll);
    return rc;
}

static const hdfe_boot_model_ops fit_ops = {fit_bytes, fit_new, fit_free, fit_replicate,
                                             fit_thread_bytes};

/* Memory for replicates in flight: half of physical memory */
static size_t boot_budget(void)
{
    uint64_t physical = ctools_physical_memory_bytes();
    return physical ? (size_t)(physical / 2) : (size_t)1 << 31;
}

/* Scalar __creghdfe_bs_par: the replicates boot_chunk runs at once, so that
 * the ado sizes its chunks as multiples of it */
static ST_retcode save_parallel(const boot_model *m)
{
    ST_int par = hdfe_boot_parallel(&m->design, &fit_ops, m, m->K_keep + 1, 0, 0,
                                    ctools_get_max_threads(), boot_budget());
    return ctools_scal_save("__creghdfe_bs_par", (ST_double)(par > 0 ? par : 1));
}

/* ========================================================================
 * Plugin calls
 * ======================================================================== */

static ST_retcode check_retained(void)
{
    if (!g_boot)
        return 198;
    if (SF_nobs() != g_boot->nobs)
    {
        SF_error("creghdfe: the data changed between the fit and its bootstrap\n");
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

ST_retcode creghdfe_boot_chunk(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    boot_model *m = g_boot;
    int first = -1, count = -1, total = -1;
    if (ctools_parse_int_option(args, "first", &first) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "count", &count) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "total", &total) != CTOOLS_PARSE_OK || first < 0 ||
        count < 1 || total < 2 || first + count > total)
    {
        SF_error("creghdfe: invalid boot_chunk arguments\n");
        return 198;
    }
    const int stata_draws = ctools_parse_bool_option(args, "stata");
    const int supplied_draws = ctools_parse_bool_option(args, "supplied");
    if (first == 0)
    {
        hdfe_boot_acc_free(&m->acc);
        m->fe_target = ctools_parse_bool_option(args, "fe");
        m->dse_target = m->fe_target && ctools_parse_bool_option(args, "dse");
        if ((rc = hdfe_boot_acc_init(&m->acc, m->K_keep + 1, m->fe_target ? m->L_total : 0,
                                     m->dse_target ? m->N : 0, total)))
            return rc;
        if (m->fe_target && ctools_parse_bool_option(args, "ci"))
        {
            rc = hdfe_boot_acc_keep_fe(&m->acc, boot_budget());
            if (rc == 909)
                SF_error("creghdfe: fevce(bootstrap, ci) keeps every replicate's fixed effects, "
                         "and they do not fit in memory; use fewer reps() or omit ci\n");
            if (rc)
                return rc;
        }
        if (!stata_draws && !supplied_draws)
        {
            uint64_t seed = 0;
            if (ctools_parse_seed_option(args, &seed) != CTOOLS_PARSE_OK)
            {
                SF_error("creghdfe: boot_chunk needs seedhi= and seedlo=\n");
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
                                       m->unit_of_row, m->obs_of_row, "creghdfe", &supplied)))
            return rc;
    }

    rc = hdfe_boot_run(&m->design, supplied, &fit_ops, m, first, count, ctools_get_max_threads(),
                       boot_budget(), &m->acc);
    free(supplied);
    if (rc == 909)
        SF_error("creghdfe: not enough memory for one bootstrap replicate\n");
    if (rc)
        return rc;
    if (ctools_parse_bool_option(args, "rows"))
        rc = hdfe_boot_store_rows(&m->acc, first, count, "__creghdfe_bs_rows");
    for (int r = 0; r < count && !rc; r++)
        rc = ctools_mat_store("__creghdfe_bs_status", 1, r + 1,
                              (ST_double)m->acc.status[first + r]);
    return rc;
}

static ST_retcode store_vector(const char *name, const ST_double *v, ST_int n, ST_int rows)
{
    ST_retcode rc = 0;
    for (ST_int r = 0; r < rows && !rc; r++)
    {
        for (ST_int k = 0; k < n && !rc; k++)
        {
            ST_double x = v[(size_t)r * n + k];
            rc = ctools_mat_store(name, r + 1, k + 1, isfinite(x) ? x : SV_missval);
        }
    }
    return rc;
}

static int cmp_int(const void *a, const void *b)
{
    ST_int x = *(const ST_int *)a, y = *(const ST_int *)b;
    return (x > y) - (x < y);
}

/* fevce(bootstrap): each saved effect's bootstrap SE in variables 1..G of the
 * call (G more with the valid-replicate counts when nreps, then d's SE when
 * dse), on the estimation sample; matrix __creghdfe_bs_fe_table (G x 5:
 * levels, never estimable, always estimable, min and median valid
 * replicates) and scalars __creghdfe_bs_fe_split_reps, __creghdfe_bs_fe_N_reps */
static ST_retcode finish_effects(const boot_model *m, const char *args)
{
    const ST_int G = m->G, N = m->N;
    const int nreps = ctools_parse_bool_option(args, "nreps");
    const int ci = ctools_parse_bool_option(args, "ci");
    const ST_int nd = G * (1 + nreps + 2 * ci); /* d's SE follows the per-effect outputs */
    if (SF_nvars() != nd + m->dse_target || (ci && !m->acc.fe_draws))
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
    ST_double *d_se = m->dse_target ? malloc((size_t)N * sizeof(ST_double)) : NULL;
    ST_int *sorted = malloc((size_t)(m->L_total > 0 ? m->L_total : 1) * sizeof(ST_int));
    if (!rc && (!fe_se || (m->dse_target && !d_se) || !sorted))
        rc = 920;
    if (!rc)
        hdfe_boot_fe_se(&m->acc, fe_se, d_se);
    for (ST_int i = 0; i < N && !rc; i++)
    {
        for (ST_int g = 0; g < G && !rc; g++)
        {
            ST_int node = m->offset[g] + m->levels[g][i] - 1;
            ST_double se = fe_se[node];
            rc = SF_vstore(g + 1, m->obs_of_row[i], isfinite(se) ? se : SV_missval);
            if (!rc && nreps)
                rc = SF_vstore(G + g + 1, m->obs_of_row[i], (ST_double)m->acc.fe_n[node]);
            if (!rc && ci)
                rc = SF_vstore(G * (1 + nreps) + g + 1, m->obs_of_row[i],
                               isfinite(lb[node]) ? lb[node] : SV_missval);
            if (!rc && ci)
                rc = SF_vstore(G * (2 + nreps) + g + 1, m->obs_of_row[i],
                               isfinite(ub[node]) ? ub[node] : SV_missval);
        }
        if (!rc && m->dse_target)
            rc = SF_vstore(nd + 1, m->obs_of_row[i], isfinite(d_se[i]) ? d_se[i] : SV_missval);
    }
    for (ST_int g = 0; g < G && !rc; g++)
    {
        ST_int L = m->num_levels[g], never = 0, always = 0;
        for (ST_int l = 0; l < L; l++)
        {
            ST_int n = m->acc.fe_n[m->offset[g] + l];
            sorted[l] = n;
            never += n == 0;
            always += n == m->acc.R;
        }
        qsort(sorted, (size_t)L, sizeof(ST_int), cmp_int);
        ST_double median = L % 2 ? sorted[L / 2] : 0.5 * (sorted[L / 2 - 1] + sorted[L / 2]);
        const ST_double row[5] = {(ST_double)L, (ST_double)never, (ST_double)always,
                                  (ST_double)sorted[0], median};
        for (int c = 0; c < 5 && !rc; c++)
            rc = ctools_mat_store("__creghdfe_bs_fe_table", g + 1, c + 1, row[c]);
    }
    if (!rc)
        rc = ctools_scal_save("__creghdfe_bs_fe_split_reps", (ST_double)m->acc.fe_split_reps);
    if (!rc)
        rc = ctools_scal_save("__creghdfe_bs_fe_N_reps", (ST_double)m->acc.R);
    free(fe_se);
    free(lb);
    free(ub);
    free(d_se);
    free(sorted);
    return rc;
}

ST_retcode creghdfe_boot_finish(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    boot_model *m = g_boot;
    if (m->fe_target && (rc = finish_effects(m, args)))
        return rc;
    if (!ctools_parse_bool_option(args, "coef"))
        return 0;
    double level = ctools_parse_double_option(args, "level", 95.0);
    int mse = ctools_parse_bool_option(args, "mse");
    hdfe_boot_summary s;
    rc = hdfe_boot_summarize(&m->acc, m->b_hat, mse, level, &s);
    if (!rc)
        rc = ctools_scal_save("__creghdfe_bs_N_reps", (ST_double)s.N_reps);
    if (!rc)
        rc = ctools_scal_save("__creghdfe_bs_N_misreps", (ST_double)s.N_misreps);
    if (rc == 2000)
    {
        SF_error("creghdfe: insufficient observations to compute bootstrap standard errors\n");
        return rc;
    }
    const ST_int Kb = m->K_keep + 1;
    if (!rc)
        rc = store_vector("__creghdfe_bs_V", s.V, Kb, Kb);
    if (!rc)
        rc = store_vector("__creghdfe_bs_b", s.b_bs, Kb, 1);
    if (!rc)
        rc = store_vector("__creghdfe_bs_se", s.se, Kb, 1);
    if (!rc)
        rc = store_vector("__creghdfe_bs_bias", s.bias, Kb, 1);
    if (!rc)
        rc = store_vector("__creghdfe_bs_z0", s.z0, Kb, 1);
    if (!rc)
        rc = store_vector("__creghdfe_bs_ci_normal", s.ci_normal, Kb, 2);
    if (!rc)
        rc = store_vector("__creghdfe_bs_ci_percentile", s.ci_percentile, Kb, 2);
    if (!rc)
        rc = store_vector("__creghdfe_bs_ci_bc", s.ci_bc, Kb, 2);
    for (int k = 0; k < HDFE_BOOT_NSTATUS && !rc; k++)
        rc = ctools_mat_store("__creghdfe_bs_counts", 1, k + 1,
                              (ST_double)m->acc.status_count[k]);
    hdfe_boot_summary_free(&s);
    return rc;
}

/* ========================================================================
 * wildbootstrap(): the wild bootstrap of the coefficients (plan 5.9, 7.5)
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

/* 0-based ids of plugin variable var on the compact rows, numbered in
 * ascending order of its values (as egen group() and boottest number
 * clusters) */
static ST_retcode dense_ids(const boot_model *m, ST_int var, ST_int *ids, ST_int *n)
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

/* Inverse of an n x n symmetric positive definite matrix (column-major) by
 * Cholesky; 198 when it is not positive definite */
static ST_retcode spd_inverse(ST_int n, const ST_double *M, ST_double *inv)
{
    ST_double *L = calloc((size_t)n * (size_t)n, sizeof(ST_double));
    ST_double *Li = calloc((size_t)n * (size_t)n, sizeof(ST_double));
    ST_retcode rc = (!L || !Li) ? 920 : 0;
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
    for (ST_int j = 0; j < n && !rc; j++)
    {
        Li[IX(j, j)] = 1.0 / L[IX(j, j)];
        for (ST_int i = j + 1; i < n; i++)
        {
            ST_double t = 0;
            for (ST_int k = j; k < i; k++)
                t -= L[IX(i, k)] * Li[IX(k, j)];
            Li[IX(i, j)] = t / L[IX(i, i)];
        }
    }
    for (ST_int a = 0; a < n && !rc; a++)
    {
        for (ST_int b = 0; b < n; b++)
        {
            ST_double t = 0;
            for (ST_int k = a > b ? a : b; k < n; k++)
                t += Li[IX(k, a)] * Li[IX(k, b)];
            inv[IX(a, b)] = t;
        }
    }
#undef IX
    free(L);
    free(Li);
    return rc;
}

static void wild_free(boot_model *m)
{
    free(m->wb_vce);
    free(m->wb_boot);
    free(m->wb_test);
    free(m->wb_X);
    free(m->wb_u);
    free(m->wb_A);
    free(m->wb_beta);
    m->wb_vce = m->wb_boot = m->wb_test = NULL;
    m->wb_X = m->wb_u = m->wb_A = m->wb_beta = NULL;
    hdfe_wild_free(&m->wild);
    if (m->wb_S_ready)
    {
        free(m->wb_S.weights);
        m->wb_S.weights = NULL;
        if (m->wb_S.factors)
            ctools_hdfe_state_cleanup(&m->wb_S);
        m->wb_S_ready = 0;
    }
}

ST_retcode creghdfe_wild_setup(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    boot_model *m = g_boot;
    wild_free(m);
    const ST_int N = m->N, K_x = m->K - 1, Kw = m->K_keep;
    const int has_vce = ctools_parse_bool_option(args, "vce");
    const int has_boot = ctools_parse_bool_option(args, "boot");
    const int partialled = ctools_parse_bool_option(args, "partialled");
    const int force_projected = ctools_parse_bool_option(args, "projected");
    m->wb_null = !ctools_parse_bool_option(args, "nonull");
    if (SF_nvars() != has_vce + has_boot || Kw < 1)
        return 198;

    /* Tested coefficients: regressor positions (1-based) in matrix
     * __creghdfe_wb_test, mapped to the kept columns */
    char tmat[] = "__creghdfe_wb_test";
    const ST_int T = SF_col(tmat);
    if (T < 1)
        return 198;
    ST_int *kept_of = malloc((size_t)(K_x > 0 ? K_x : 1) * sizeof(ST_int));
    m->wb_test = malloc((size_t)T * sizeof(ST_int));
    if (!kept_of || !m->wb_test)
    {
        free(kept_of);
        return 920;
    }
    ST_int nk = 0;
    for (ST_int k = 0; k < K_x; k++)
        kept_of[k] = m->is_collinear[k] ? -1 : nk++;
    if (nk != Kw)
    {
        free(kept_of);
        return 198;
    }
    for (ST_int j = 0; j < T && !rc; j++)
    {
        ST_double x;
        if ((rc = SF_mat_el(tmat, 1, j + 1, &x)))
            break;
        int p = (int)x;
        if (p < 1 || p > K_x || kept_of[p - 1] < 0)
        {
            SF_error("creghdfe: wildbootstrap() cannot test an omitted coefficient\n");
            rc = 198;
            break;
        }
        m->wb_test[j] = kept_of[p - 1];
    }
    m->wb_ntest = T;

    /* The partialled kept regressors and y in the units of the data, the
     * residuals, and A = (X'WX)^-1 */
    m->wb_X = malloc((size_t)N * (size_t)Kw * sizeof(ST_double));
    m->wb_u = malloc((size_t)N * sizeof(ST_double));
    m->wb_A = malloc((size_t)Kw * (size_t)Kw * sizeof(ST_double));
    m->wb_beta = malloc((size_t)Kw * sizeof(ST_double));
    ST_double *M = calloc((size_t)Kw * (size_t)Kw, sizeof(ST_double));
    if (!rc && (!m->wb_X || !m->wb_u || !m->wb_A || !m->wb_beta || !M))
        rc = 920;
    for (ST_int k = 0; k < K_x && !rc; k++)
    {
        if (kept_of[k] < 0)
            continue;
        ST_double *x = m->wb_X + (size_t)kept_of[k] * N;
        const ST_double *t = m->tilde + (size_t)(k + 1) * N;
        for (ST_int i = 0; i < N; i++)
            x[i] = t[i] * m->stdevs[k + 1];
    }
    free(kept_of);
    if (!rc)
    {
        for (ST_int a = 0; a < Kw; a++)
            m->wb_beta[a] = m->b_hat[a];
        for (ST_int i = 0; i < N; i++)
        {
            ST_double u = m->tilde[i] * m->stdevs[0];
            for (ST_int a = 0; a < Kw; a++)
                u -= m->wb_X[(size_t)a * N + i] * m->wb_beta[a];
            m->wb_u[i] = u;
        }
        for (ST_int a = 0; a < Kw; a++)
            for (ST_int b = 0; b <= a; b++)
            {
                const ST_double *xa = m->wb_X + (size_t)a * N, *xb = m->wb_X + (size_t)b * N;
                ST_double s = 0;
                for (ST_int i = 0; i < N; i++)
                    s += (m->weights ? m->weights[i] : 1.0) * xa[i] * xb[i];
                M[a + (size_t)b * Kw] = M[b + (size_t)a * Kw] = s;
            }
        rc = spd_inverse(Kw, M, m->wb_A);
    }
    free(M);

    /* Clusters: the vce() clusters and the bootstrap clusters (the vce()
     * clusters unless bootcluster()); none: each row */
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

    /* Every absorbed effect nested in the bootstrap clusters: the aggregate
     * path is exact */
    int nested = 1;
    for (ST_int g = 0; g < m->G && !rc && nested; g++)
    {
        ST_int *seen = malloc((size_t)m->num_levels[g] * sizeof(ST_int));
        if (!seen)
        {
            rc = 920;
            break;
        }
        for (ST_int l = 0; l < m->num_levels[g]; l++)
            seen[l] = -1;
        for (ST_int i = 0; i < N && nested; i++)
        {
            ST_int l = m->levels[g][i] - 1, c = m->wb_boot ? m->wb_boot[i] : i;
            if (seen[l] < 0)
                seen[l] = c;
            else if (seen[l] != c)
                nested = 0;
        }
        free(seen);
    }
    if (rc)
    {
        wild_free(m);
        return rc;
    }
    m->wb_projected = force_projected || (!nested && !partialled);
    const int path = m->wb_projected ? (nested ? 3 : 1) : (nested ? 0 : 2);
    rc = ctools_scal_save("__creghdfe_wb_C", (ST_double)m->wb_C);
    if (!rc)
        rc = ctools_scal_save("__creghdfe_wb_G", (ST_double)m->wb_G);
    if (!rc)
        rc = ctools_scal_save("__creghdfe_wb_path", (ST_double)path);
    return rc;
}

/* The projected path's partialling: the full sample's factors with the
 * fit's weights, built once */
static ST_retcode wild_projection(boot_model *m)
{
    const ST_int N = m->N, G = m->G;
    HDFE_State *S = &m->wb_S;
    memset(S, 0, sizeof(*S));
    ST_int *mask = malloc((size_t)N * sizeof(ST_int));
    ST_int *level_max = malloc((size_t)G * sizeof(ST_int));
    ST_double *omega = malloc((size_t)N * sizeof(ST_double));
    ST_retcode rc = (!mask || !level_max || !omega) ? 920 : 0;
    ST_double sum = 0;
    for (ST_int i = 0; i < N && !rc; i++)
    {
        mask[i] = 1;
        omega[i] = m->weights ? m->weights[i] : 1.0;
        sum += omega[i];
    }
    for (ST_int g = 0; g < G && !rc; g++)
        level_max[g] = m->num_levels[g];
    if (!rc)
        rc = hdfe_build_factors(S, G, m->levels, level_max, N, mask, N, omega, 0);
    free(level_max);
    if (rc)
    {
        free(mask);
        free(omega);
        if (S->factors)
            ctools_hdfe_state_cleanup(S);
        return rc;
    }
    S->has_weights = 1;
    S->weight_type = 2;
    S->maxiter = m->maxiter;
    S->tolerance = m->tolerance;
    S->num_threads = ctools_get_max_threads();
    S->factors_initialized = 1;
    if (m->slope)
        rc = creghdfe_setup_slopes(S, &m->slope_in, (const ST_double *const *)m->slope, mask, N,
                                   omega);
    free(mask);
    if (rc)
    {
        free(omega);
        ctools_hdfe_state_cleanup(S);
        return rc;
    }
    S->weights = omega;
    S->sum_weights = sum;
    m->wb_S_ready = 1;
    if (ctools_hdfe_alloc_buffers(S, 0, 1 << 20) != 0)
        return 920;
    return 0;
}

static ST_retcode wild_project(void *ctx, ST_double *cols, ST_int ncols)
{
    boot_model *m = ctx;
    HDFE_SolveResult r = partial_out_columns(&m->wb_S, cols, m->N, ncols, m->wb_S.num_threads);
    if (r.status == 430)
        SF_error("creghdfe: the wild bootstrap's partialling did not converge; increase iterate()\n");
    return r.status;
}

ST_retcode creghdfe_wild_chunk(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    boot_model *m = g_boot;
    int first = -1, count = -1, total = -1;
    if (ctools_parse_int_option(args, "first", &first) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "count", &count) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "total", &total) != CTOOLS_PARSE_OK || first < 0 ||
        count < 1 || total < 1 || first + count > total || !m->wb_X)
    {
        SF_error("creghdfe: invalid wild_chunk arguments\n");
        return 198;
    }
    const int enumerate = ctools_parse_bool_option(args, "enum");
    if (first == 0)
    {
        hdfe_wild_free(&m->wild);
        hdfe_wild_input in;
        memset(&in, 0, sizeof(in));
        in.N = m->N;
        in.K = m->K_keep;
        in.X = m->wb_X;
        in.w = m->weights;
        in.resid = m->wb_u;
        in.A = m->wb_A;
        in.beta = m->wb_beta;
        in.vce_id = m->wb_vce;
        in.G = m->wb_G;
        in.boot_id = m->wb_boot;
        in.C = m->wb_C;
        in.n_test = m->wb_ntest;
        in.test = m->wb_test;
        in.null_imposed = m->wb_null;
        in.projected = m->wb_projected;
        if ((rc = hdfe_wild_setup(&in, total, &m->wild)))
            return rc;
        if (m->wb_projected && !m->wb_S_ready && (rc = wild_projection(m)))
            return rc;
    }
    if (m->wild.B != total || m->wild.done != first)
        return 198;

    ST_double *v = NULL;
    const ST_int C = m->wb_C;
    if (!enumerate)
    {
        /* The ado packs the chunk's weights replicate after replicate */
        const size_t values = (size_t)count * (size_t)C;
        if (!packed_vars_ok(values))
            return 198;
        const ST_int nobs = SF_nobs();
        v = malloc(values * sizeof(ST_double));
        if (!v)
            return 920;
        for (size_t q = 0; q < values && !rc; q++)
            rc = packed_value(nobs, q, &v[q]);
        if (rc)
        {
            free(v);
            return rc;
        }
    }
    /* A quarter of physical memory: hdfe_wild_run batches its columns at
     * about 2 GB and uses more only for one replicate that needs it */
    uint64_t physical = ctools_physical_memory_bytes();
    size_t budget = physical ? (size_t)(physical / 4) : (size_t)1 << 30;
    rc = hdfe_wild_run(&m->wild, v, first, count, ctools_get_max_threads(),
                       m->wb_projected ? wild_project : NULL, m, budget);
    free(v);
    if (rc == 909)
        SF_error("creghdfe: not enough memory for one wild bootstrap replicate\n");
    return rc;
}

typedef struct
{
    const hdfe_wild_model *w;
    ST_int j;
    int ptype;
} wild_pctx;

static double wild_p(void *ctx, double b0)
{
    const wild_pctx *c = ctx;
    return hdfe_wild_pvalue(c->w, c->j, b0, c->ptype, NULL);
}

static ST_double stata_value(double x)
{
    return isfinite(x) ? x : SV_missval;
}

ST_retcode creghdfe_wild_finish(const char *args)
{
    ST_retcode rc = check_retained();
    if (rc)
        return rc;
    boot_model *m = g_boot;
    const hdfe_wild_model *w = &m->wild;
    if (!w->five || w->done != w->B)
        return 198;
    int ptype = 0;
    if (ctools_parse_int_option(args, "ptype", &ptype) != CTOOLS_PARSE_OK || ptype < 0 ||
        ptype > 3)
        ptype = 0;
    const double level = ctools_parse_double_option(args, "level", 95.0);
    const int ci = ctools_parse_bool_option(args, "ci");
    ST_int min_feasible = w->B;
    for (ST_int j = 0; j < w->n_test && !rc; j++)
    {
        ST_int F = 0;
        const double p = hdfe_wild_pvalue(w, j, 0.0, ptype, &F);
        if (F < min_feasible)
            min_feasible = F;
        rc = ctools_mat_store("__creghdfe_wb_p", 1, j + 1, stata_value(p));
        if (rc || !ci)
            continue;
        wild_pctx c = {w, j, ptype};
        hdfe_wild_ci set;
        hdfe_wild_invert(wild_p, &c, w->b_hat[j], w->se[j], 1.0 - level / 100.0, &set);
        const double row[7] = {set.lo,       set.hi,       set.open_lo,         set.open_hi,
                               set.disjoint, set.empty,    set.failed > 0};
        for (int k = 0; k < 2 && !rc; k++)
            rc = ctools_mat_store("__creghdfe_wb_ci", j + 1, k + 1, stata_value(row[k]));
        for (int k = 0; k < 5 && !rc; k++)
            rc = ctools_mat_store("__creghdfe_wb_cistat", j + 1, k + 1, row[2 + k]);
    }
    if (!rc)
        rc = ctools_scal_save("__creghdfe_wb_reps", (ST_double)w->B);
    if (!rc)
        rc = ctools_scal_save("__creghdfe_wb_misreps", (ST_double)(w->B - min_feasible));
    return rc;
}
