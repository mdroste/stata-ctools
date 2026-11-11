/*
 * civreghdfe_boot.c
 *
 * vce(bootstrap) and fevce(bootstrap) for civreghdfe (see civreghdfe_boot.h).
 * A replicate is the full-sample model re-fit with each row weighted by its
 * draw multiplicity, as the expanded data would be: singletons judged on the
 * copies, the partialled full-sample columns projected again under the
 * replicate weights (they give what the raw columns would, and start closer
 * to it), partial(), the fit's own collinearity rules, then the SPI-free
 * coefficient core (ivest_fit_core) with the copy-aware moment covariance
 * for GMM2S and CUE (docs/development/BOOTSTRAP_PLAN.md, 9.2 and 9.3).
 * Part of the ctools Stata plugin suite
 */

/* Checked SPI access: indices resolve in the plugin varlist */
#define SD_SAFEMODE
#include "civreghdfe_boot.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "civreghdfe_impl.h"
#include "civreghdfe_estimate.h"
#include "civreghdfe_omega.h"
#include "civreghdfe_wre.h"
#include "hdfe/hdfe_boot.h"
#include "hdfe/hdfe_boot_stata.h"
#include "hdfe/hdfe_wild.h"
#include "hdfe/hdfe_fe_recover.h"
#include "hdfe/hdfe_sample.h"
#include "hdfe/hdfe_solver.h"
#include "../ctools_config.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_matrix.h"
#include "../ctools_parse.h"
#include "../ctools_platform.h"
#include "../ctools_spi.h"

typedef struct
{
    ST_int nobs; /* SF_nobs() when retained: the data may not change between calls */
    ST_int N, G, Kx, Ke, Kz, ncols;
    ST_double *y, *X_endog, *X_exog, *Z_excl, *tilde, *weights;
    ST_int weight_type;
    ST_int **levels; /* G x N, 1-based full-sample levels */
    ST_int *num_levels;
    ST_int *partial; /* Kx, or NULL */
    ST_int K_x;      /* exogenous regressors after partial() */
    ST_int *collin_x, *collin_z;
    ST_int endo_in_z;
    ST_int K_b;
    ST_double *b_hat;
    ST_int est_method;
    ST_double kclass_user, fuller_alpha;
    ST_int center, dofminus, maxiter, keep_singletons;
    ST_double tolerance;
    /* The GMM moment covariance of the replicates */
    ST_int s_kind, s_units;
    ST_int *cluster_ids; /* N, 1-based (s_kind cluster on another variable), or NULL */
    /* Resampling frame */
    ST_int *obs_of_row;
    ST_int *unit_of_row, *stratum_start;
    ST_double *unit_freq;
    hdfe_boot_design design;
    int cluster_mode, fw_rows;
    hdfe_boot_acc acc;
    /* fevce(bootstrap): levels as nodes offset[g] + level - 1, their
     * full-sample connected sets (all effects; each consecutive pair) */
    ST_int *offset, L_total;
    ST_int *joint_comp, n_joint;
    ST_int **pair_comp, *n_pair;
    ST_int *fe_targets; /* G flags: the saved effects */
    int fe_target, dse_target;
    /* wildbootstrap() (civreghdfe_wre.h) */
    ST_double N_eff;
    ST_int K_iv_fit;
    civ_wre_input wb_in;
    ST_double *wb_y, *wb_Y2, *wb_X1, *wb_X2;
    ST_int *wb_test, *wb_vce, *wb_boot;
    int wb_ready;
    civ_wre_model *wre;
    HDFE_State wb_S;
    int wb_S_ready;
} civ_boot_model;

static civ_boot_model *g_civboot = NULL;

static ST_retcode civb_save_parallel(const civ_boot_model *m);

static void civb_wild_free(civ_boot_model *m)
{
    civ_wre_free(m->wre);
    m->wre = NULL;
    free(m->wb_y);
    free(m->wb_Y2);
    free(m->wb_X1);
    free(m->wb_X2);
    free(m->wb_test);
    free(m->wb_vce);
    free(m->wb_boot);
    m->wb_y = m->wb_Y2 = m->wb_X1 = m->wb_X2 = NULL;
    m->wb_test = m->wb_vce = m->wb_boot = NULL;
    m->wb_ready = 0;
    if (m->wb_S_ready)
    {
        free(m->wb_S.weights);
        m->wb_S.weights = NULL;
        if (m->wb_S.factors)
            ctools_hdfe_state_cleanup(&m->wb_S);
        m->wb_S_ready = 0;
    }
}

static void civb_model_free(civ_boot_model *m)
{
    if (!m)
        return;
    civb_wild_free(m);
    free(m->y);
    free(m->X_endog);
    free(m->X_exog);
    free(m->Z_excl);
    free(m->tilde);
    free(m->weights);
    for (ST_int g = 0; m->levels && g < m->G; g++)
        free(m->levels[g]);
    free(m->levels);
    free(m->num_levels);
    free(m->partial);
    free(m->collin_x);
    free(m->collin_z);
    free(m->b_hat);
    free(m->cluster_ids);
    free(m->obs_of_row);
    free(m->unit_of_row);
    free(m->stratum_start);
    free(m->unit_freq);
    hdfe_boot_acc_free(&m->acc);
    free(m->offset);
    free(m->joint_comp);
    for (ST_int g = 0; m->pair_comp && g < m->G; g++)
        free(m->pair_comp[g]);
    free(m->pair_comp);
    free(m->n_pair);
    free(m->fe_targets);
    free(m);
}

void civreghdfe_boot_release(void)
{
    civb_model_free(g_civboot);
    g_civboot = NULL;
}

/* ========================================================================
 * Retain
 * ======================================================================== */

typedef struct
{
    double stratum, cluster;
    ST_int row;
} civb_unit_key;

static int civb_unit_key_cmp(const void *a, const void *b)
{
    const civb_unit_key *x = a, *y = b;
    if (x->stratum != y->stratum)
        return x->stratum < y->stratum ? -1 : 1;
    if (x->cluster != y->cluster)
        return x->cluster < y->cluster ? -1 : 1;
    return (x->row > y->row) - (x->row < y->row);
}

/* Units in bsample's order: strata ascending; within a stratum, clusters in
 * ascending value, or rows in data order. A cluster that spans strata is a
 * unit within each stratum, as bsample numbers clusters within strata.
 * Frequency-weighted rows resampled one by one stand for fw observations
 * each: those are resampled (plan 5.3). */
static ST_retcode civb_build_frame(civ_boot_model *m, ST_int unit_var, ST_int strata_var)
{
    const ST_int N = m->N;
    civb_unit_key *keys = malloc((size_t)N * sizeof(*keys));
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
    qsort(keys, (size_t)N, sizeof(*keys), civb_unit_key_cmp);
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
    free(keys);
    m->fw_rows = m->weights && m->weight_type == 2 && !m->cluster_mode;
    if (m->fw_rows)
    {
        if (!(m->unit_freq = malloc((size_t)U * sizeof(ST_double))))
            return 920;
        for (ST_int i = 0; i < N; i++)
            m->unit_freq[m->unit_of_row[i]] = m->weights[i];
    }
    m->design.U = U;
    m->design.S = S;
    m->design.stratum_start = m->stratum_start;
    m->design.unit_freq = m->unit_freq;
    m->design.seed = 0;
    return hdfe_boot_design_check(&m->design);
}

/* Connected sets of levels over the rows kept[0..n), for effects g0..g1
 * together (as creghdfe_boot.c): comp numbers each node's set densely */
static ST_retcode civb_level_sets(const civ_boot_model *m, const ST_int *kept, ST_int n, ST_int g0,
                                  ST_int g1, ST_int *comp, ST_int *n_comp)
{
    const ST_int *off = m->offset;
    ST_int nodes = off[g1 + 1] - off[g0];
    ctools_UnionFind *uf = ctools_uf_create(nodes);
    if (!uf)
        return 920;
    for (ST_int r = 0; r < n; r++)
    {
        ST_int i = kept[r];
        ST_int a = m->levels[g0][i] - 1;
        for (ST_int g = g0 + 1; g <= g1; g++)
            ctools_uf_union(uf, a, off[g] - off[g0] + m->levels[g][i] - 1);
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

static ST_retcode civb_full_sample_sets(civ_boot_model *m)
{
    const ST_int G = m->G, N = m->N;
    m->offset = malloc((size_t)(G + 1) * sizeof(ST_int));
    if (!m->offset)
        return 920;
    m->offset[0] = 0;
    for (ST_int g = 0; g < G; g++)
        m->offset[g + 1] = m->offset[g] + m->num_levels[g];
    m->L_total = m->offset[G];
    if (G < 2)
        return 0;
    ST_int *all = malloc((size_t)N * sizeof(ST_int));
    m->joint_comp = malloc((size_t)m->L_total * sizeof(ST_int));
    m->pair_comp = calloc((size_t)G, sizeof(ST_int *));
    m->n_pair = calloc((size_t)G, sizeof(ST_int));
    if (!all || !m->joint_comp || !m->pair_comp || !m->n_pair)
    {
        free(all);
        return 920;
    }
    for (ST_int i = 0; i < N; i++)
        all[i] = i;
    ST_retcode rc = civb_level_sets(m, all, N, 0, G - 1, m->joint_comp, &m->n_joint);
    for (ST_int g = 1; g < G && !rc; g++)
    {
        m->pair_comp[g] = malloc((size_t)(m->num_levels[g - 1] + m->num_levels[g]) * sizeof(ST_int));
        if (!m->pair_comp[g])
            rc = 920;
        else
            rc = civb_level_sets(m, all, N, g - 1, g, m->pair_comp[g], &m->n_pair[g]);
    }
    free(all);
    return rc;
}

/* Marks split[j] for every full-sample joint set j that the replicate's rows
 * leave in more than one piece, for all effects together or for any
 * consecutive pair (plan 5.7); present[node] flags the drawn levels */
static ST_retcode civb_replicate_splits(const civ_boot_model *m, const ST_int *kept, ST_int n,
                                        const unsigned char *present, unsigned char *split)
{
    const ST_int G = m->G;
    ST_retcode rc = 0;
    ST_int *comp = malloc((size_t)m->L_total * sizeof(ST_int));
    ST_int *first = malloc((size_t)m->L_total * sizeof(ST_int));
    if (!comp || !first)
        rc = 920;
    for (ST_int g1 = 0; g1 < G && !rc; g1++)
    {
        ST_int lo = g1 == 0 ? 0 : g1 - 1, hi = g1 == 0 ? G - 1 : g1, k = 0;
        if ((rc = civb_level_sets(m, kept, n, lo, hi, comp, &k)))
            break;
        const ST_int *full = g1 == 0 ? m->joint_comp : m->pair_comp[g1];
        ST_int n_full = g1 == 0 ? m->n_joint : m->n_pair[g1];
        for (ST_int c = 0; c < n_full; c++)
            first[c] = -1;
        for (ST_int v = 0; v < m->offset[hi + 1] - m->offset[lo]; v++)
        {
            ST_int node = m->offset[lo] + v;
            if (!present[node])
                continue;
            if (first[full[v]] < 0)
                first[full[v]] = comp[v];
            else if (first[full[v]] != comp[v])
                split[m->joint_comp[node]] = 1;
        }
    }
    free(comp);
    free(first);
    return rc;
}

static ST_double *civb_copy_doubles(const ST_double *src, size_t n)
{
    ST_double *dst = malloc((n ? n : 1) * sizeof(*dst));
    if (dst && n)
        memcpy(dst, src, n * sizeof(*dst));
    return dst;
}

static ST_int *civb_copy_ints(const ST_int *src, size_t n)
{
    ST_int *dst = malloc((n ? n : 1) * sizeof(*dst));
    if (dst && n)
        memcpy(dst, src, n * sizeof(*dst));
    return dst;
}

static ST_int civb_scalar(const char *name, ST_int dflt)
{
    ST_double v = 0.0;
    char buf[64];
    snprintf(buf, sizeof(buf), "%s", name);
    return SF_scal_use(buf, &v) == 0 && v >= 0 && v < (ST_double)INT_MAX ? (ST_int)v : dflt;
}

ST_retcode civreghdfe_boot_retain(const civreghdfe_boot_source *src)
{
    civreghdfe_boot_release();
    civ_boot_model *m = calloc(1, sizeof(*m));
    if (!m)
        return 920;
    const ST_int N = src->N, G = src->G;
    m->nobs = SF_nobs();
    m->N = N;
    m->G = G;
    m->Kx = src->K_exog;
    m->Ke = src->K_endog;
    m->Kz = src->K_excl;
    m->ncols = src->K_excl + 1 + src->K_endog + src->K_exog;
    m->weight_type = src->weight_type;
    m->endo_in_z = src->endo_in_z;
    m->K_b = src->K_b;
    m->est_method = src->est_method;
    m->kclass_user = src->kclass_user;
    m->fuller_alpha = src->fuller_alpha;
    m->center = src->center;
    m->dofminus = src->dofminus;
    m->maxiter = src->maxiter;
    m->keep_singletons = src->keep_singletons;
    m->tolerance = src->tolerance;
    m->N_eff = src->N_eff;
    m->K_iv_fit = src->K_iv;
    ST_int n_partial = 0;
    for (ST_int k = 0; src->partial && k < src->K_exog; k++)
        n_partial += src->partial[k] != 0;
    m->K_x = src->K_exog - n_partial;
    const ST_int K_coll = m->K_x + src->K_endog;
    ST_retcode rc = 0;

    /* The replicates' GMM moment covariance: the fit's (robust or one-way
     * cluster); the other estimators do not need one for b */
    m->s_kind = CIV_OMEGA_IID;
    if (src->est_method == 4 || src->est_method == 5)
    {
        if (src->kernel_type > 0 || src->kiefer || src->vce_type > 2)
        {
            SF_error("civreghdfe: the bootstrap of GMM2S and CUE supports the robust and one-way "
                     "cluster moment covariance\n");
            civb_model_free(m);
            return 198;
        }
        m->s_kind = src->vce_type == 2 ? CIV_OMEGA_CLUSTER
                    : src->vce_type == 1 ? CIV_OMEGA_ROBUST
                                         : CIV_OMEGA_IID;
    }
    m->s_units = m->s_kind == CIV_OMEGA_CLUSTER && civb_scalar("__civreghdfe_bs_s_units", 0) == 1;

    m->y = civb_copy_doubles(src->y, (size_t)N);
    m->X_endog = civb_copy_doubles(src->X_endog, (size_t)N * (size_t)src->K_endog);
    m->X_exog = civb_copy_doubles(src->X_exog, (size_t)N * (size_t)src->K_exog);
    m->Z_excl = civb_copy_doubles(src->Z_excl, (size_t)N * (size_t)src->K_excl);
    m->tilde = civb_copy_doubles(src->tilde, (size_t)N * (size_t)m->ncols);
    m->weights = src->weights ? civb_copy_doubles(src->weights, (size_t)N) : NULL;
    m->partial = src->partial ? civb_copy_ints(src->partial, (size_t)src->K_exog) : NULL;
    m->collin_x = civb_copy_ints(src->collin_x, (size_t)K_coll);
    m->collin_z = civb_copy_ints(src->collin_z, (size_t)src->K_excl);
    m->b_hat = civb_copy_doubles(src->b_hat, (size_t)src->K_b);
    if (m->s_kind == CIV_OMEGA_CLUSTER && !m->s_units)
        m->cluster_ids = src->cluster_ids ? civb_copy_ints(src->cluster_ids, (size_t)N) : NULL;
    m->levels = calloc((size_t)(G > 0 ? G : 1), sizeof(ST_int *));
    m->num_levels = malloc((size_t)(G > 0 ? G : 1) * sizeof(ST_int));
    m->obs_of_row = malloc((size_t)N * sizeof(ST_int));
    if (!m->y || (src->K_endog && !m->X_endog) || (src->K_exog && !m->X_exog) ||
        (src->K_excl && !m->Z_excl) || !m->tilde || (src->weights && !m->weights) ||
        (src->partial && !m->partial) || !m->collin_x || !m->collin_z || !m->b_hat ||
        (m->s_kind == CIV_OMEGA_CLUSTER && !m->s_units && !m->cluster_ids) || !m->levels ||
        !m->num_levels || !m->obs_of_row)
        rc = 920;
    for (ST_int i = 0; !rc && i < N; i++)
        m->obs_of_row[i] = (ST_int)src->obs_map[i];
    for (ST_int g = 0; g < G && !rc; g++)
    {
        m->num_levels[g] = src->factors[g].num_levels;
        if (!(m->levels[g] = civb_copy_ints(src->factors[g].levels, (size_t)N)))
            rc = 920;
    }
    if (!rc)
        rc = civb_build_frame(m, civb_scalar("__civreghdfe_bs_unit_idx", 0),
                              civb_scalar("__civreghdfe_bs_strata_idx", 0));
    if (!rc)
        rc = civb_full_sample_sets(m);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_bs_U", hdfe_boot_draw_total(&m->design));
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_bs_S", (ST_double)m->design.S);
    if (!rc)
        rc = hdfe_boot_save_dup_strata(&m->design, "civreghdfe");
    if (!rc)
        rc = civb_save_parallel(m);
    if (rc)
    {
        civb_model_free(m);
        return rc;
    }
    g_civboot = m;
    return 0;
}

/* ========================================================================
 * One replicate
 * ======================================================================== */

/* The raw values of partialled column c ([excluded, y, endogenous, exogenous]) */
static const ST_double *civb_raw(const civ_boot_model *m, ST_int c)
{
    const size_t N = (size_t)m->N;
    if (c < m->Kz)
        return m->Z_excl + (size_t)c * N;
    c -= m->Kz;
    if (c == 0)
        return m->y;
    c -= 1;
    if (c < m->Ke)
        return m->X_endog + (size_t)c * N;
    return m->X_exog + (size_t)(c - m->Ke) * N;
}

static size_t civb_fit_bytes(const void *mv)
{
    const civ_boot_model *m = mv;
    const size_t total = (size_t)m->Kx + (size_t)m->ncols;
    return (size_t)m->N * (total + 16) * sizeof(ST_double) +
           (size_t)m->N * ((size_t)m->G * 2 + 6) * sizeof(ST_int) +
           (size_t)m->L_total * (2 * sizeof(ST_double) + 3 * sizeof(ST_int)) +
           total * total * 8 * sizeof(ST_double);
}

static size_t civb_fit_thread_bytes(const void *mv)
{
    const civ_boot_model *m = mv;
    return (size_t)m->N * 3 * sizeof(ST_double) + (size_t)m->L_total * sizeof(ST_double);
}

static void *civb_fit_new(const void *mv)
{
    (void)mv;
    return malloc(1);
}

static void civb_fit_free(void *ws)
{
    free(ws);
}

static void civb_fit_io_free(civ_fit_io *io)
{
    ST_double *arrays[] = {io->ZtZ,   io->ZtZ_inv,   io->ZtX,   io->Zty,         io->XtPzX,
                           io->XtPzy, io->temp1,     io->X_all, io->resid,       io->X_hat,
                           io->ZtZ_inv_Zty, io->XtX, io->Xty,   io->XkX,         io->Xky,
                           io->XkX_copy, io->beta_temp, io->S_main, io->gmm_hessian_inv,
                           io->W_bread, io->W_zx};
    for (size_t a = 0; a < sizeof(arrays) / sizeof(arrays[0]); a++)
        free(arrays[a]);
    civ_omega_free(&io->om);
}

/* Compact the columns of X (N x K, column-major) whose flag is 0 */
static ST_int civb_compact(ST_double *X, ST_int N, ST_int K, const ST_int *drop)
{
    ST_int k2 = 0;
    for (ST_int k = 0; k < K; k++)
    {
        if (drop[k])
            continue;
        if (k2 != k)
            memmove(X + (size_t)k2 * N, X + (size_t)k * N, (size_t)N * sizeof(ST_double));
        k2++;
    }
    return k2;
}

typedef struct
{
    ST_int *rows, *mask, *krow, **lv, *coll_x, *coll_z, *isp, *drop, *cl;
    ST_double *cnt, *om1, *w2, *mult2, *A, *ss_pre, *ss_unc, *buf, *M, *M2, *beta, *err;
    ST_double **alpha;
    unsigned char *present, *split;
    /* partial() with saved effects: the partialled columns' cross products
     * before the FWL step, Kp x (Kp + 1 + Kx + Ke) ([W, y, X_exog, X_endog],
     * each column partialled by the fixed effects) */
    ST_double *pw;
    /* their FWL solve: WW (Kp x Kp), rhs and gamma (Kp each); dropped (Kp) */
    ST_double *pgw;
    ST_int *pdrop;
} civb_rep_buf;

static void civb_rep_buf_free(civb_rep_buf *b, ST_int G)
{
    for (ST_int g = 0; b->lv && g < G; g++)
        free(b->lv[g]);
    for (ST_int g = 0; b->alpha && g < G; g++)
        free(b->alpha[g]);
    free(b->lv);
    free(b->alpha);
    void *p[] = {b->rows, b->mask, b->krow, b->coll_x, b->coll_z, b->isp, b->drop, b->cl,
                 b->cnt,  b->om1,  b->w2,   b->mult2,  b->A,      b->ss_pre, b->ss_unc, b->buf,
                 b->M,    b->M2,   b->beta, b->err,    b->present, b->split, b->pw,
                 b->pgw,  b->pdrop};
    for (size_t a = 0; a < sizeof(p) / sizeof(p[0]); a++)
        free(p[a]);
}

/* The replicate's fixed effects as the fit's savefe recovers them (y - Xb
 * absorbed with creghdfe's normalization, the constant added to the first
 * effect), in full-sample level order; levels not drawn, or in a connected
 * set the replicate splits, stay missing. Also d per compact row. */
static ST_retcode civb_effects(const civ_boot_model *m, HDFE_State *S, civb_rep_buf *b, ST_int N2,
                               const ST_double *y_t, const ST_double *Xx_t, const ST_double *Xe_t,
                               ST_int K_x, ST_int K_e, hdfe_boot_rep *out)
{
    const ST_int G = m->G, N = m->N;
    ST_retcode rc = 0;
    if (!(b->err = malloc((size_t)N2 * sizeof(ST_double))) ||
        !(b->alpha = calloc((size_t)G, sizeof(ST_double *))) ||
        !(b->present = calloc((size_t)(m->L_total > 0 ? m->L_total : 1), 1)) ||
        !(b->split = calloc((size_t)(m->n_joint > 0 ? m->n_joint : 1), 1)))
        return 920;
    /* partial(): the partialled regressors' coefficients, by FWL: those of
     * y - Xb on them, all partialled by the fixed effects (the structural
     * residual is orthogonal to them) */
    ST_int Kp = 0;
    for (ST_int k = 0; m->partial && k < m->Kx; k++)
        Kp += m->partial[k] != 0;
    const ST_double *gamma = NULL;
    if (Kp > 0)
    {
        if (!b->pw)
            return 198;
        const ST_int ld = Kp + 1 + m->Kx + m->Ke;
        if (!(b->pgw = malloc(((size_t)Kp * (size_t)Kp + 2 * (size_t)Kp) * sizeof(ST_double))) ||
            !(b->pdrop = calloc((size_t)Kp, sizeof(ST_int))))
            return 920;
        ST_double *WW = b->pgw, *rhs = WW + (size_t)Kp * Kp, *gam = rhs + Kp;
        ST_int *dropped = b->pdrop;
        gamma = gam;
        for (ST_int p = 0; p < Kp; p++)
        {
            const ST_double *row = b->pw + (size_t)p * ld;
            for (ST_int q = 0; q < Kp; q++)
                WW[(size_t)p * Kp + q] = row[q];
            ST_double v = row[Kp];
            /* the kept exogenous regressors (post-partial order), then the
             * endogenous ones */
            for (ST_int k = 0, j = 0, c = 0; k < m->Kx; k++)
            {
                if (m->partial[k])
                    continue;
                if (!b->coll_x[j])
                    v -= row[Kp + 1 + k] * b->beta[c++];
                j++;
            }
            for (ST_int k = 0, c = 0; k < m->Ke; k++)
            {
                if (b->coll_x[m->K_x + k])
                    continue;
                v -= row[Kp + 1 + m->Kx + k] * b->beta[K_x + c];
                c++;
            }
            rhs[p] = v;
        }
        /* Cholesky with redundant columns dropped (their coefficient is
         * zero). The replicate's FWL step decided which columns it dropped
         * (absorbed by the fixed effects, or collinear with earlier ones:
         * isp == 2); drop exactly those, since an absorbed column is left
         * as rounding noise whose own diagonal no relative test catches.
         * The relative test below only guards what that step kept. */
        for (ST_int k = 0, p = 0; b->isp && k < m->Kx; k++)
        {
            if (m->partial[k])
                dropped[p++] = b->isp[k] == 2;
        }
        for (ST_int j = 0; j < Kp; j++)
        {
            ST_double dj = WW[(size_t)j * Kp + j];
            const ST_double scale = dj;
            for (ST_int k = 0; k < j; k++)
                dj -= WW[(size_t)j * Kp + k] * WW[(size_t)j * Kp + k];
            if (dropped[j] || !(dj > 1e-12 * scale) || !(scale > 0.0))
            {
                dropped[j] = 1;
                for (ST_int i2 = j; i2 < Kp; i2++)
                    WW[(size_t)i2 * Kp + j] = 0.0;
                continue;
            }
            const ST_double ljj = sqrt(dj);
            WW[(size_t)j * Kp + j] = ljj;
            for (ST_int i2 = j + 1; i2 < Kp; i2++)
            {
                ST_double s = WW[(size_t)i2 * Kp + j];
                for (ST_int k = 0; k < j; k++)
                    s -= WW[(size_t)i2 * Kp + k] * WW[(size_t)j * Kp + k];
                WW[(size_t)i2 * Kp + j] = s / ljj;
            }
        }
        for (ST_int j = 0; j < Kp; j++)
        {
            if (dropped[j])
            {
                gam[j] = 0.0;
                continue;
            }
            ST_double s = rhs[j];
            for (ST_int k = 0; k < j; k++)
                s -= WW[(size_t)j * Kp + k] * gam[k];
            gam[j] = s / WW[(size_t)j * Kp + j];
        }
        for (ST_int j = Kp - 1; j >= 0; j--)
        {
            if (dropped[j])
                continue;
            ST_double s = gam[j];
            for (ST_int k = j + 1; k < Kp; k++)
                s -= WW[(size_t)k * Kp + j] * gam[k];
            gam[j] = s / WW[(size_t)j * Kp + j];
        }
    }
    /* y - Xb on the raw columns less the structural residual on the
     * partialled ones */
    ST_double sw = 0.0, sd = 0.0;
    for (ST_int r = 0; r < N2; r++)
    {
        const ST_int i = b->krow[r];
        ST_double e = y_t[r], d = m->y[i];
        for (ST_int k = 0, j = 0, c = 0, p = 0; k < m->Kx; k++)
        {
            if (m->partial && m->partial[k])
            {
                d -= m->X_exog[(size_t)k * N + i] * gamma[p++];
                continue;
            }
            if (b->coll_x[j++])
                continue;
            e -= Xx_t[(size_t)c * N2 + r] * b->beta[c];
            d -= m->X_exog[(size_t)k * N + i] * b->beta[c];
            c++;
        }
        for (ST_int k = 0, c = 0; k < m->Ke; k++)
        {
            if (b->coll_x[m->K_x + k])
                continue;
            e -= Xe_t[(size_t)c * N2 + r] * b->beta[K_x + c];
            d -= m->X_endog[(size_t)k * N + i] * b->beta[K_x + c];
            c++;
        }
        b->err[r] = d - e;
        sw += b->w2[r];
        sd += b->w2[r] * b->err[r];
    }
    (void)K_e;
    const ST_double cons = sd / sw;
    for (ST_int r = 0; r < N2; r++)
        b->err[r] -= cons;
    for (ST_int g = 0; g < G; g++)
    {
        if (!(b->alpha[g] = calloc((size_t)S->factors[g].num_levels, sizeof(ST_double))))
            return 920;
    }
    rc = hdfe_fe_recover(S, b->err, b->alpha);
    if (rc == 430)
        return 0; /* not converged: this replicate's effects stay missing */
    if (rc)
        return rc;
    for (ST_int r = 0; r < N2; r++)
    {
        for (ST_int g = 0; g < G; g++)
            b->present[m->offset[g] + m->levels[g][b->krow[r]] - 1] = 1;
    }
    if (G >= 2)
    {
        if ((rc = civb_replicate_splits(m, b->krow, N2, b->present, b->split)))
            return rc;
        for (ST_int j = 0; j < m->n_joint; j++)
            out->fe_split |= b->split[j];
    }
    for (ST_int r = 0; r < N2; r++)
    {
        for (ST_int g = 0; g < G; g++)
        {
            ST_int node = m->offset[g] + m->levels[g][b->krow[r]] - 1;
            if (G >= 2 && b->split[m->joint_comp[node]])
                continue;
            ST_double a = b->alpha[g][S->factors[g].levels[r] - 1];
            out->fe[node] = g == 0 ? a + cons : a;
        }
    }
    for (ST_int i = 0; out->d && i < N; i++)
    {
        ST_double d = 0.0;
        for (ST_int g = 0; g < G; g++)
            d += out->fe[m->offset[g] + m->levels[g][i] - 1];
        out->d[i] = d;
    }
    return 0;
}

static ST_retcode civb_fit_replicate(const void *mv, void *ws, const ST_int *unit_mult,
                                     hdfe_boot_rep *out)
{
    const civ_boot_model *m = mv;
    const ST_int N = m->N, G = m->G, Kx = m->Kx, Ke = m->Ke, Kz = m->Kz, ncols = m->ncols;
    const ST_int total_cols = Kx + ncols;
    const ST_int col_excl = Kx, col_y = Kx + Kz, col_endog = col_y + 1, col_exog = col_endog + Ke;
    (void)ws;
    ST_retcode rc = 0;
    HDFE_State S;
    memset(&S, 0, sizeof(S));
    civb_rep_buf b;
    memset(&b, 0, sizeof(b));
    civ_fit_io io;
    memset(&io, 0, sizeof(io));
    int io_used = 0;

    for (ST_int v = 0; out->fe && v < m->L_total; v++)
        out->fe[v] = NAN;
    for (ST_int i = 0; out->d && i < N; i++)
        out->d[i] = NAN;

    /* Rows drawn at least once, and their counts in the expanded data (the
     * singleton rule and N): copies, times the fweight when clusters or
     * weightless rows are drawn; drawn fweighted observations otherwise */
    b.rows = malloc((size_t)N * sizeof(ST_int));
    b.cnt = malloc((size_t)N * sizeof(ST_double));
    b.om1 = malloc((size_t)N * sizeof(ST_double));
    b.mask = malloc((size_t)N * sizeof(ST_int));
    b.lv = calloc((size_t)(G > 0 ? G : 1), sizeof(ST_int *));
    if (!b.rows || !b.cnt || !b.om1 || !b.mask || !b.lv)
    {
        rc = 920;
        goto done;
    }
    ST_int N1 = 0;
    for (ST_int i = 0; i < N; i++)
    {
        const ST_int mult = unit_mult[m->unit_of_row[i]];
        if (mult <= 0)
            continue;
        const ST_double w = m->weights ? m->weights[i] : 1.0;
        b.rows[N1] = i;
        b.cnt[N1] = (double)mult * ((m->weights && m->weight_type == 2 && !m->fw_rows) ? w : 1.0);
        b.om1[N1] = (double)mult * (m->fw_rows ? 1.0 : w);
        N1++;
    }
    if (N1 == 0)
    {
        out->status = HDFE_BOOT_EMPTY;
        goto done;
    }
    for (ST_int g = 0; g < G; g++)
    {
        if (!(b.lv[g] = malloc((size_t)N1 * sizeof(ST_int))))
        {
            rc = 920;
            goto done;
        }
        for (ST_int r = 0; r < N1; r++)
            b.lv[g][r] = m->levels[g][b.rows[r]];
    }
    ST_int N2 = N1;
    if (G > 0 && !m->keep_singletons)
    {
        if (ctools_remove_singletons(b.lv, G, N1, b.mask, b.cnt, 0) < 0)
        {
            rc = 920;
            goto done;
        }
        N2 = 0;
        for (ST_int r = 0; r < N1; r++)
            N2 += b.mask[r] != 0;
    }
    else
    {
        for (ST_int r = 0; r < N1; r++)
            b.mask[r] = 1;
    }
    if (N2 == 0)
    {
        out->status = HDFE_BOOT_EMPTY;
        goto done;
    }

    /* The kept rows: weights (aweights and pweights normalized to the
     * expanded count, as the fit normalizes them to N), multiplicities */
    b.krow = malloc((size_t)N2 * sizeof(ST_int));
    b.w2 = malloc((size_t)N2 * sizeof(ST_double));
    b.mult2 = malloc((size_t)N2 * sizeof(ST_double));
    b.A = malloc((size_t)N2 * (size_t)total_cols * sizeof(ST_double));
    b.ss_pre = calloc((size_t)total_cols, sizeof(ST_double));
    b.ss_unc = calloc((size_t)total_cols, sizeof(ST_double));
    b.buf = malloc((size_t)N2 * sizeof(ST_double));
    if (!b.krow || !b.w2 || !b.mult2 || !b.A || !b.ss_pre || !b.ss_unc || !b.buf)
    {
        rc = 920;
        goto done;
    }
    ST_double N_eff = 0.0, sum_w = 0.0, sum_m = 0.0;
    for (ST_int r = 0, r2 = 0; r < N1; r++)
    {
        if (!b.mask[r])
            continue;
        const ST_int i = b.rows[r];
        b.krow[r2] = i;
        b.w2[r2] = b.om1[r];
        b.mult2[r2] = (double)unit_mult[m->unit_of_row[i]];
        N_eff += b.cnt[r];
        sum_w += b.w2[r2];
        sum_m += b.mult2[r2];
        r2++;
    }
    /* (the factors' weighted counts below must use the same weights as the
     * projection, so om1 is normalized too) */
    if (m->weights && m->weight_type != 2 && sum_w > 0.0)
    {
        const ST_double scale = sum_m / sum_w;
        for (ST_int r2 = 0; r2 < N2; r2++)
            b.w2[r2] *= scale;
        for (ST_int r = 0; r < N1; r++)
            b.om1[r] *= scale;
    }
    const ST_int wtype = (m->weights && m->weight_type != 2) ? m->weight_type : 2;
    for (ST_int c = 0; c < ncols; c++)
    {
        const ST_double *t = m->tilde + (size_t)c * N;
        ST_double *dst = b.A + (size_t)(col_excl + c) * N2;
        for (ST_int r2 = 0; r2 < N2; r2++)
            dst[r2] = t[b.krow[r2]];
    }

    /* The fixed effects, partialled out again under the replicate weights,
     * each column at unit scale as the fit projects it */
    if (G > 0)
    {
        ST_int level_max[64];
        if (G > 64)
        {
            rc = 198;
            goto done;
        }
        for (ST_int g = 0; g < G; g++)
            level_max[g] = m->num_levels[g];
        if ((rc = hdfe_build_factors(&S, G, b.lv, level_max, N1, b.mask, N2, b.om1, 0)))
            goto done;
        S.has_weights = 1;
        S.weight_type = 2;
        S.maxiter = m->maxiter;
        S.tolerance = m->tolerance;
        S.num_threads = out->threads > 1 ? out->threads : 1;
        S.factors_initialized = 1;
        if (!(S.weights = malloc((size_t)N2 * sizeof(ST_double))))
        {
            rc = 920;
            goto done;
        }
        memcpy(S.weights, b.w2, (size_t)N2 * sizeof(ST_double));
        S.sum_weights = sum_m;
        if (ctools_hdfe_alloc_buffers(&S, 0, ncols) != 0)
        {
            rc = 920;
            goto done;
        }
        ST_double scale[1024];
        ST_double *scales = ncols <= 1024 ? scale : malloc((size_t)ncols * sizeof(ST_double));
        if (!scales)
        {
            rc = 920;
            goto done;
        }
        for (ST_int c = 0; c < ncols; c++)
        {
            const ST_double *raw = civb_raw(m, c);
            ST_double largest = 0.0;
            for (ST_int r2 = 0; r2 < N2; r2++)
                largest = fmax(largest, fabs(raw[b.krow[r2]]));
            int exponent = 0;
            scales[c] = 1.0;
            if (largest > 0.0 && isfinite(largest))
            {
                frexp(largest, &exponent);
                scales[c] = ldexp(1.0, exponent - 1);
                ST_double *x = b.A + (size_t)(col_excl + c) * N2;
                for (ST_int r2 = 0; r2 < N2; r2++)
                    x[r2] /= scales[c];
            }
        }
        HDFE_SolveResult proj = partial_out_columns(&S, b.A + (size_t)col_excl * N2, N2, ncols, 1);
        for (ST_int c = 0; c < ncols; c++)
        {
            if (scales[c] == 1.0)
                continue;
            ST_double *x = b.A + (size_t)(col_excl + c) * N2;
            for (ST_int r2 = 0; r2 < N2; r2++)
                x[r2] *= scales[c];
        }
        if (scales != scale)
            free(scales);
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
    }

    /* Pre-partialling sums of squares of the raw columns (the fit's rules
     * for absorbed columns and redundant partial() columns) */
    for (ST_int c = 0; c < ncols; c++)
    {
        const ST_double *raw = civb_raw(m, c);
        for (ST_int r2 = 0; r2 < N2; r2++)
            b.buf[r2] = raw[b.krow[r2]];
        b.ss_pre[col_excl + c] = civ_wss_dev(b.buf, b.w2, N2);
        b.ss_unc[col_excl + c] = civ_wss(b.buf, b.w2, N2);
    }
    if (Kx > 0)
        memcpy(b.A, b.A + (size_t)col_exog * N2, (size_t)N2 * (size_t)Kx * sizeof(ST_double));
    ST_double *ss_exog = b.ss_pre + col_exog, *ss_unc_exog = b.ss_unc + col_exog;
    const ST_double *ss_unc_endog = b.ss_unc + col_endog, *ss_unc_excl = b.ss_unc + col_excl;

    /* partial(): with saved effects, the partialled columns' cross
     * products first (civb_effects recovers their coefficients) */
    if (m->partial && out->fe && G > 0)
    {
        ST_int Kp = 0;
        for (ST_int k = 0; k < Kx; k++)
            Kp += m->partial[k] != 0;
        const ST_int ld = Kp + 1 + Kx + Ke;
        if (!(b.pw = malloc((size_t)(Kp > 0 ? Kp : 1) * (size_t)ld * sizeof(ST_double))))
        {
            rc = 920;
            goto done;
        }
        for (ST_int k = 0, p = 0; k < Kx; k++)
        {
            if (!m->partial[k])
                continue;
            const ST_double *wk = b.A + (size_t)(col_exog + k) * N2;
            ST_double *row = b.pw + (size_t)p * ld;
            for (ST_int c = 0; c < ld; c++)
            {
                const ST_double *x;
                if (c < Kp)
                {
                    ST_int q = -1, kk = 0;
                    for (; kk < Kx; kk++)
                    {
                        if (m->partial[kk] && ++q == c)
                            break;
                    }
                    x = b.A + (size_t)(col_exog + kk) * N2;
                }
                else if (c == Kp)
                    x = b.A + (size_t)col_y * N2;
                else if (c < Kp + 1 + Kx)
                    x = b.A + (size_t)(col_exog + c - Kp - 1) * N2;
                else
                    x = b.A + (size_t)(col_endog + c - Kp - 1 - Kx) * N2;
                ST_double s = 0.0;
                for (ST_int r2 = 0; r2 < N2; r2++)
                    s += b.w2[r2] * wk[r2] * x[r2];
                row[c] = s;
            }
            p++;
        }
    }

    /* partial() */
    ST_int K_x = Kx, K_iv = Kx + Kz;
    if (m->partial)
    {
        const ST_double absorb_tol = (m->tolerance / 10.0 < 1e-6) ? m->tolerance / 10.0 : 1e-6;
        if (!(b.isp = civb_copy_ints(m->partial, (size_t)Kx)))
        {
            rc = 920;
            goto done;
        }
        ST_retcode fwl_rc = civ_fwl_partial_out(b.A, N2, total_cols, col_exog, 0, Kx, b.isp, b.w2,
                                                ss_exog, absorb_tol, 0);
        if (fwl_rc == 920)
        {
            rc = 920;
            goto done;
        }
        if (fwl_rc)
        {
            out->status = HDFE_BOOT_COLLINEAR;
            goto done;
        }
        ST_int n = 0;
        for (ST_int k = 0; k < Kx; k++)
        {
            if (b.isp[k])
                continue;
            if (n != k)
            {
                memmove(b.A + (size_t)(col_exog + n) * N2, b.A + (size_t)(col_exog + k) * N2,
                        (size_t)N2 * sizeof(ST_double));
                ss_exog[n] = ss_exog[k];
                ss_unc_exog[n] = ss_unc_exog[k];
            }
            n++;
        }
        /* Z: the exogenous columns kept, then the excluded instruments */
        ST_int z = 0;
        for (ST_int k = 0; k < Kx; k++)
        {
            if (b.isp[k])
                continue;
            if (z != k)
                memmove(b.A + (size_t)z * N2, b.A + (size_t)k * N2, (size_t)N2 * sizeof(ST_double));
            z++;
        }
        for (ST_int k = 0; k < Kz; k++, z++)
            memmove(b.A + (size_t)z * N2, b.A + (size_t)(Kx + k) * N2, (size_t)N2 * sizeof(ST_double));
        K_x = n;
        K_iv = n + Kz;
    }
    if (K_x != m->K_x)
    {
        out->status = HDFE_BOOT_COLLINEAR;
        goto done;
    }
    ST_double *y_t = b.A + (size_t)col_y * N2;
    ST_double *Xe_t = b.A + (size_t)col_endog * N2;
    ST_double *Xx_t = b.A + (size_t)col_exog * N2;
    ST_double *Z_t = b.A;

    /* The fit's collinearity rules: columns absorbed by the fixed effects,
     * then invsym() of X'WX in [endogenous, exogenous] order. A replicate
     * whose omitted columns differ from the fit's is missing, as in the
     * bootstrap prefix. */
    const ST_int K_total = K_x + Ke;
    b.coll_x = calloc((size_t)(K_total > 0 ? K_total : 1), sizeof(ST_int));
    b.drop = calloc((size_t)(K_total + K_iv + 1), sizeof(ST_int));
    b.M = malloc((size_t)(K_total + K_iv + 1) * (size_t)(K_total + K_iv + 1) * sizeof(ST_double));
    b.M2 = malloc((size_t)(K_total + K_iv + 1) * (size_t)(K_total + K_iv + 1) * sizeof(ST_double));
    b.coll_z = calloc((size_t)(K_iv > 0 ? K_iv : 1), sizeof(ST_int));
    if (!b.coll_x || !b.drop || !b.M || !b.M2 || !b.coll_z)
    {
        rc = 920;
        goto done;
    }
    for (ST_int k = 0; k < K_x; k++)
        b.coll_x[k] = !(civ_wss(Xx_t + (size_t)k * N2, b.w2, N2) > CIV_FE_NOISE * ss_unc_exog[k]);
    for (ST_int k = 0; k < Ke; k++)
        b.coll_x[K_x + k] =
            !(civ_wss(Xe_t + (size_t)k * N2, b.w2, N2) > CIV_FE_NOISE * ss_unc_endog[k]);
    if (K_total > 0)
    {
        const ST_double *X_ne = Ke > 0 ? Xe_t : Xx_t;
        ctools_matmul_atdb(X_ne, X_ne, b.w2, N2, K_total, K_total, b.M);
        for (ST_int k = 0; k < Ke; k++)
            b.drop[k] = b.coll_x[K_x + k];
        for (ST_int k = 0; k < K_x; k++)
            b.drop[Ke + k] = b.coll_x[k];
        if (civ_invsym_drops(b.M, K_total, 0, b.drop) < 0)
        {
            rc = 920;
            goto done;
        }
        for (ST_int k = 0; k < Ke; k++)
            b.coll_x[K_x + k] = b.drop[k];
        for (ST_int k = 0; k < K_x; k++)
            b.coll_x[k] = b.drop[Ke + k];
    }
    for (ST_int k = 0; k < K_total; k++)
    {
        if ((b.coll_x[k] != 0) != (m->collin_x[k] != 0))
        {
            out->status = HDFE_BOOT_COLLINEAR;
            goto done;
        }
    }
    const ST_int K_xk = civb_compact(Xx_t, N2, K_x, b.coll_x);
    const ST_int K_ek = civb_compact(Xe_t, N2, Ke, b.coll_x + K_x);
    /* Z's exogenous block follows X's; the excluded instruments move down */
    {
        ST_int z = civb_compact(Z_t, N2, K_x, b.coll_x);
        for (ST_int k = 0; k < Kz; k++, z++)
        {
            if (z != K_x + k)
                memmove(Z_t + (size_t)z * N2, Z_t + (size_t)(K_x + k) * N2,
                        (size_t)N2 * sizeof(ST_double));
        }
        K_iv = z;
    }
    /* Excluded instruments absorbed by the fixed effects, then invsym() of
     * Z'WZ with the exogenous regressors swept first */
    for (ST_int j = 0; j < Kz; j++)
        b.coll_z[K_xk + j] =
            !(civ_wss(Z_t + (size_t)(K_xk + j) * N2, b.w2, N2) > CIV_FE_NOISE * ss_unc_excl[j]);
    if (K_iv > 0)
    {
        ctools_matmul_atdb(Z_t, Z_t, b.w2, N2, K_iv, K_iv, b.M);
        if (civ_invsym_drops(b.M, K_iv, K_xk, b.coll_z) < 0)
        {
            rc = 920;
            goto done;
        }
    }
    for (ST_int j = 0; j < Kz; j++)
    {
        if ((b.coll_z[K_xk + j] != 0) != (m->collin_z[j] != 0))
        {
            out->status = HDFE_BOOT_COLLINEAR;
            goto done;
        }
    }
    K_iv = civb_compact(Z_t, N2, K_iv, b.coll_z);
    /* Endogenous regressors collinear with the instruments */
    if (K_ek > 0 && K_iv > 0)
    {
        const ST_int KJ = K_iv + K_ek;
        ST_double *J = b.M, *zx = b.M2;
        ctools_matmul_atdb(Z_t, Z_t, b.w2, N2, K_iv, K_iv, b.M2);
        for (ST_int c = 0; c < K_iv; c++)
            for (ST_int r = 0; r < K_iv; r++)
                J[(size_t)c * KJ + r] = b.M2[(size_t)c * K_iv + r];
        ctools_matmul_atdb(Z_t, Xe_t, b.w2, N2, K_iv, K_ek, zx);
        for (ST_int c = 0; c < K_ek; c++)
            for (ST_int r = 0; r < K_iv; r++)
            {
                J[(size_t)(K_iv + c) * KJ + r] = zx[(size_t)c * K_iv + r];
                J[(size_t)r * KJ + K_iv + c] = zx[(size_t)c * K_iv + r];
            }
        ctools_matmul_atdb(Xe_t, Xe_t, b.w2, N2, K_ek, K_ek, zx);
        for (ST_int c = 0; c < K_ek; c++)
            for (ST_int r = 0; r < K_ek; r++)
                J[(size_t)(K_iv + c) * KJ + K_iv + r] = zx[(size_t)c * K_ek + r];
        for (ST_int k = 0; k < KJ; k++)
            b.drop[k] = 0;
        ST_int j_drops = civ_invsym_drops(J, KJ, KJ, b.drop);
        if (j_drops < 0)
        {
            rc = 920;
            goto done;
        }
        if (j_drops != m->endo_in_z)
        {
            out->status = HDFE_BOOT_COLLINEAR;
            goto done;
        }
    }
    const ST_int K_kept = K_xk + K_ek;
    if (K_iv < K_kept || K_kept == 0)
    {
        out->status = HDFE_BOOT_COLLINEAR;
        goto done;
    }

    /* The GMM clusters: the vce() clusters, or the drawn units, each copy
     * its own cluster */
    ST_int n_cl = 0;
    if (m->s_kind == CIV_OMEGA_CLUSTER)
    {
        if (!(b.cl = malloc((size_t)N2 * sizeof(ST_int))))
        {
            rc = 920;
            goto done;
        }
        for (ST_int r2 = 0; r2 < N2; r2++)
            b.cl[r2] = m->s_units ? m->unit_of_row[b.krow[r2]] + 1 : m->cluster_ids[b.krow[r2]];
        if (ctools_remap_cluster_ids(b.cl, N2, &n_cl) != 0)
        {
            rc = 920;
            goto done;
        }
    }

    if (!(b.beta = calloc((size_t)K_kept, sizeof(ST_double))))
    {
        rc = 920;
        goto done;
    }
    ST_double lambda = 1.0;
    io.report = 0;
    io.y = y_t;
    io.X_exog = K_xk > 0 ? Xx_t : NULL;
    io.X_endog = K_ek > 0 ? Xe_t : NULL;
    io.Z = Z_t;
    io.weights = b.w2;
    io.weight_type = wtype;
    io.N = N2;
    io.N_eff = N_eff;
    io.K_exog = K_xk;
    io.K_endog = K_ek;
    io.K_iv = K_iv;
    io.K_total = K_kept;
    io.beta = b.beta;
    io.lambda_out = &lambda;
    io.est_method = m->est_method;
    io.kclass_user = m->kclass_user;
    io.fuller_alpha = m->fuller_alpha;
    io.s_kind = m->s_kind;
    io.s_noniid = m->s_kind != CIV_OMEGA_IID;
    io.center = m->center;
    io.cluster_ids = b.cl;
    io.num_clusters = n_cl;
    io.dofminus = m->dofminus;
    io.mult = b.mult2;
    io.mult_groups = m->s_units;
    io_used = 1;
    ST_retcode fit_rc = ivest_fit_core(&io);
    if (fit_rc == 920)
    {
        rc = 920;
        goto done;
    }
    if (fit_rc)
    {
        out->status = HDFE_BOOT_OTHER;
        goto done;
    }
    for (ST_int c = 0; c < K_kept; c++)
    {
        if (!isfinite(b.beta[c]))
        {
            out->status = HDFE_BOOT_OTHER;
            goto done;
        }
    }

    /* e(b)'s layout: [endogenous, exogenous], zeros where omitted */
    for (ST_int i = 0, ce = 0, cx = 0; i < K_total; i++)
    {
        if (i < K_x)
            out->b[Ke + i] = b.coll_x[i] ? 0.0 : b.beta[cx++];
        else
            out->b[i - K_x] = b.coll_x[i] ? 0.0 : b.beta[K_xk + ce++];
    }
    out->status = HDFE_BOOT_OK;
    if (out->fe && G > 0)
        rc = civb_effects(m, &S, &b, N2, y_t, Xx_t, Xe_t, K_xk, K_ek, out);

done:
    if (io_used)
        civb_fit_io_free(&io);
    free(S.weights);
    S.weights = NULL;
    if (S.factors)
        ctools_hdfe_state_cleanup(&S);
    civb_rep_buf_free(&b, G);
    return rc;
}

static const hdfe_boot_model_ops civb_fit_ops = {civb_fit_bytes, civb_fit_new, civb_fit_free,
                                                 civb_fit_replicate, civb_fit_thread_bytes};

/* Memory for replicates in flight: half of physical memory */
static size_t civb_budget(void)
{
    uint64_t physical = ctools_physical_memory_bytes();
    return physical ? (size_t)(physical / 2) : (size_t)1 << 31;
}

/* Scalar __civreghdfe_bs_par: the replicates boot_chunk runs at once */
static ST_retcode civb_save_parallel(const civ_boot_model *m)
{
    ST_int par = hdfe_boot_parallel(&m->design, &civb_fit_ops, m, m->K_b, 0, 0,
                                    ctools_get_max_threads(), civb_budget());
    return ctools_scal_save("__civreghdfe_bs_par", (ST_double)(par > 0 ? par : 1));
}

/* ========================================================================
 * Plugin calls
 * ======================================================================== */

static ST_retcode civb_check_retained(void)
{
    if (!g_civboot)
        return 198;
    if (SF_nobs() != g_civboot->nobs)
    {
        SF_error("civreghdfe: the data changed between the fit and its bootstrap\n");
        return 459;
    }
    return 0;
}

static ST_retcode civb_packed_value(ST_int nobs, size_t q, ST_double *x)
{
    return SF_vdata((ST_int)(q / (size_t)nobs) + 1, (ST_int)(q % (size_t)nobs) + 1, x);
}

static int civb_packed_vars_ok(size_t n)
{
    const ST_int nobs = SF_nobs();
    return nobs > 0 && n > 0 && (size_t)SF_nvars() == (n + (size_t)nobs - 1) / (size_t)nobs;
}

ST_retcode civreghdfe_boot_chunk(const char *args)
{
    ST_retcode rc = civb_check_retained();
    if (rc)
        return rc;
    civ_boot_model *m = g_civboot;
    int first = -1, count = -1, total = -1;
    if (ctools_parse_int_option(args, "first", &first) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "count", &count) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "total", &total) != CTOOLS_PARSE_OK || first < 0 ||
        count < 1 || total < 2 || first + count > total)
    {
        SF_error("civreghdfe: invalid boot_chunk arguments\n");
        return 198;
    }
    const int stata_draws = ctools_parse_bool_option(args, "stata");
    const int supplied_draws = ctools_parse_bool_option(args, "supplied");
    if (first == 0)
    {
        hdfe_boot_acc_free(&m->acc);
        m->fe_target = ctools_parse_bool_option(args, "fe") && m->G > 0;
        m->dse_target = m->fe_target && ctools_parse_bool_option(args, "dse");
        if (m->fe_target)
        {
            free(m->fe_targets);
            if (!(m->fe_targets = calloc((size_t)m->G, sizeof(ST_int))))
                return 920;
            for (ST_int g = 0; g < m->G; g++)
            {
                ST_double v = 0.0;
                if (SF_mat_el("__civreghdfe_bs_fe_targets", 1, g + 1, &v))
                    return 198;
                m->fe_targets[g] = v != 0.0;
            }
        }
        if ((rc = hdfe_boot_acc_init(&m->acc, m->K_b, m->fe_target ? m->L_total : 0,
                                     m->dse_target ? m->N : 0, total)))
            return rc;
        if (m->fe_target && ctools_parse_bool_option(args, "ci"))
        {
            rc = hdfe_boot_acc_keep_fe(&m->acc, civb_budget());
            if (rc == 909)
                SF_error("civreghdfe: fevce(bootstrap, ci) keeps every replicate's fixed effects, "
                         "and they do not fit in memory; use fewer reps() or omit ci\n");
            if (rc)
                return rc;
        }
        if (!stata_draws && !supplied_draws)
        {
            uint64_t seed = 0;
            if (ctools_parse_seed_option(args, &seed) != CTOOLS_PARSE_OK)
            {
                SF_error("civreghdfe: boot_chunk needs seedhi= and seedlo=\n");
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
                                       m->unit_of_row, m->obs_of_row, "civreghdfe", &supplied)))
            return rc;
    }

    rc = hdfe_boot_run(&m->design, supplied, &civb_fit_ops, m, first, count,
                       ctools_get_max_threads(), civb_budget(), &m->acc);
    free(supplied);
    if (rc == 909)
        SF_error("civreghdfe: not enough memory for one bootstrap replicate\n");
    if (rc)
        return rc;
    if (ctools_parse_bool_option(args, "rows"))
        rc = hdfe_boot_store_rows(&m->acc, first, count, "__civreghdfe_bs_rows");
    for (int r = 0; r < count && !rc; r++)
        rc = ctools_mat_store("__civreghdfe_bs_status", 1, r + 1,
                              (ST_double)m->acc.status[first + r]);
    return rc;
}

static ST_retcode civb_store_vector(const char *name, const ST_double *v, ST_int n, ST_int rows)
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

static int civb_cmp_int(const void *a, const void *b)
{
    ST_int x = *(const ST_int *)a, y = *(const ST_int *)b;
    return (x > y) - (x < y);
}

/* fevce(bootstrap): each saved effect's bootstrap SE in the call's
 * variables (then the valid-replicate counts with nreps, then d's SE with
 * dse), on the estimation sample; matrix __civreghdfe_bs_fe_table (saved
 * effects x 5: levels, never estimable, always estimable, min and median
 * valid replicates) and scalars __civreghdfe_bs_fe_split_reps,
 * __civreghdfe_bs_fe_N_reps */
static ST_retcode civb_finish_effects(const civ_boot_model *m, const char *args)
{
    const ST_int G = m->G, N = m->N;
    const int nreps = ctools_parse_bool_option(args, "nreps");
    ST_int n_t = 0;
    for (ST_int g = 0; g < G; g++)
        n_t += m->fe_targets[g];
    const int ci = ctools_parse_bool_option(args, "ci");
    const ST_int nd = n_t * (1 + nreps + 2 * ci); /* d's SE follows the per-effect outputs */
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
        for (ST_int g = 0, t = 0; g < G && !rc; g++)
        {
            if (!m->fe_targets[g])
                continue;
            ST_int node = m->offset[g] + m->levels[g][i] - 1;
            ST_double se = fe_se[node];
            rc = SF_vstore(t + 1, m->obs_of_row[i], isfinite(se) ? se : SV_missval);
            if (!rc && nreps)
                rc = SF_vstore(n_t + t + 1, m->obs_of_row[i], (ST_double)m->acc.fe_n[node]);
            if (!rc && ci)
                rc = SF_vstore(n_t * (1 + nreps) + t + 1, m->obs_of_row[i],
                               isfinite(lb[node]) ? lb[node] : SV_missval);
            if (!rc && ci)
                rc = SF_vstore(n_t * (2 + nreps) + t + 1, m->obs_of_row[i],
                               isfinite(ub[node]) ? ub[node] : SV_missval);
            t++;
        }
        if (!rc && m->dse_target)
            rc = SF_vstore(nd + 1, m->obs_of_row[i], isfinite(d_se[i]) ? d_se[i] : SV_missval);
    }
    for (ST_int g = 0, t = 0; g < G && !rc; g++)
    {
        if (!m->fe_targets[g])
            continue;
        ST_int L = m->num_levels[g], never = 0, always = 0;
        for (ST_int l = 0; l < L; l++)
        {
            ST_int n = m->acc.fe_n[m->offset[g] + l];
            sorted[l] = n;
            never += n == 0;
            always += n == m->acc.R;
        }
        qsort(sorted, (size_t)L, sizeof(ST_int), civb_cmp_int);
        ST_double median = L % 2 ? sorted[L / 2] : 0.5 * (sorted[L / 2 - 1] + sorted[L / 2]);
        const ST_double row[5] = {(ST_double)L, (ST_double)never, (ST_double)always,
                                  (ST_double)sorted[0], median};
        t++;
        for (int c = 0; c < 5 && !rc; c++)
            rc = ctools_mat_store("__civreghdfe_bs_fe_table", t, c + 1, row[c]);
    }
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_bs_fe_split_reps", (ST_double)m->acc.fe_split_reps);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_bs_fe_N_reps", (ST_double)m->acc.R);
    free(fe_se);
    free(lb);
    free(ub);
    free(d_se);
    free(sorted);
    return rc;
}

ST_retcode civreghdfe_boot_finish(const char *args)
{
    ST_retcode rc = civb_check_retained();
    if (rc)
        return rc;
    civ_boot_model *m = g_civboot;
    if (m->fe_target && (rc = civb_finish_effects(m, args)))
        return rc;
    if (!ctools_parse_bool_option(args, "coef"))
        return 0;
    double level = ctools_parse_double_option(args, "level", 95.0);
    int mse = ctools_parse_bool_option(args, "mse");
    hdfe_boot_summary s;
    rc = hdfe_boot_summarize(&m->acc, m->b_hat, mse, level, &s);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_bs_N_reps", (ST_double)s.N_reps);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_bs_N_misreps", (ST_double)s.N_misreps);
    if (rc == 2000)
    {
        SF_error("civreghdfe: insufficient observations to compute bootstrap standard errors\n");
        return rc;
    }
    const ST_int Kb = m->K_b;
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_V", s.V, Kb, Kb);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_b", s.b_bs, Kb, 1);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_se", s.se, Kb, 1);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_bias", s.bias, Kb, 1);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_z0", s.z0, Kb, 1);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_ci_normal", s.ci_normal, Kb, 2);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_ci_percentile", s.ci_percentile, Kb, 2);
    if (!rc)
        rc = civb_store_vector("__civreghdfe_bs_ci_bc", s.ci_bc, Kb, 2);
    for (int k = 0; k < HDFE_BOOT_NSTATUS && !rc; k++)
        rc = ctools_mat_store("__civreghdfe_bs_counts", 1, k + 1,
                              (ST_double)m->acc.status_count[k]);
    hdfe_boot_summary_free(&s);
    return rc;
}

/* ========================================================================
 * wildbootstrap(): WRE (civreghdfe_wre.h)
 * ======================================================================== */

typedef struct
{
    double v;
    ST_int row;
} civb_value_row;

static int civb_value_row_cmp(const void *a, const void *b)
{
    const civb_value_row *x = a, *y = b;
    if (x->v != y->v)
        return x->v < y->v ? -1 : 1;
    return (x->row > y->row) - (x->row < y->row);
}

/* 0-based ids of plugin variable var on the compact rows, numbered in
 * ascending order of its values (as boottest numbers clusters) */
static ST_retcode civb_dense_ids(const civ_boot_model *m, ST_int var, ST_int *ids, ST_int *n)
{
    const ST_int N = m->N;
    civb_value_row *vr = malloc((size_t)N * sizeof(*vr));
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
    qsort(vr, (size_t)N, sizeof(*vr), civb_value_row_cmp);
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

static size_t civb_wild_budget(void)
{
    uint64_t physical = ctools_physical_memory_bytes();
    /* A quarter of physical memory: the WRE sizes its caches and batches at
     * CIV_WRE_SOFT_BYTES and uses more only for one replicate that needs it */
    size_t budget = physical ? (size_t)(physical / 4) : (size_t)1 << 30;
    return budget;
}

ST_retcode civreghdfe_wild_setup(const char *args)
{
    ST_retcode rc = civb_check_retained();
    if (rc)
        return rc;
    civ_boot_model *m = g_civboot;
    civb_wild_free(m);
    const ST_int N = m->N, Kx = m->Kx, Ke = m->Ke, Kz = m->Kz, K_x = m->K_x;
    const int has_vce = ctools_parse_bool_option(args, "vce");
    const int has_boot = ctools_parse_bool_option(args, "boot");
    const int partialled = ctools_parse_bool_option(args, "partialled");
    const int force_projected = ctools_parse_bool_option(args, "projected");
    int vcetype = -1;
    if (ctools_parse_int_option(args, "vcetype", &vcetype) != CTOOLS_PARSE_OK || vcetype < 0 ||
        vcetype > 2 || (vcetype == 2) != has_vce || SF_nvars() != has_vce + has_boot)
    {
        SF_error("civreghdfe: invalid wild_setup arguments\n");
        return 198;
    }
    if (m->est_method > 3)
    {
        SF_error("civreghdfe: wildbootstrap() covers 2SLS, LIML, Fuller and k-class\n");
        return 198;
    }
    if (m->endo_in_z)
    {
        SF_error("civreghdfe: wildbootstrap() is not available when an endogenous regressor is "
                 "collinear with the instruments\n");
        return 198;
    }

    /* The model's columns, the fixed effects partialled out: kept
     * endogenous regressors; kept exogenous regressors, then partial()'s
     * (regressors here, the redundant ones left out); kept excluded
     * instruments */
    const ST_double *t_excl = m->tilde, *t_y = m->tilde + (size_t)Kz * N;
    const ST_double *t_endog = t_y + N, *t_exog = t_endog + (size_t)Ke * N;
    ST_int K2 = 0, K1 = 0, Kz2 = 0, Kp = 0;
    for (ST_int e = 0; e < Ke; e++)
        K2 += !m->collin_x[K_x + e];
    for (ST_int k = 0, ip = 0; k < Kx; k++)
    {
        if (m->partial && m->partial[k])
            Kp++;
        else
            K1 += !m->collin_x[ip++];
    }
    for (ST_int j = 0; j < Kz; j++)
        Kz2 += !m->collin_z[j];
    if (K2 < 1)
    {
        SF_error("civreghdfe: wildbootstrap() needs an endogenous regressor\n");
        return 198;
    }
    m->wb_y = malloc((size_t)N * sizeof(ST_double));
    m->wb_Y2 = malloc((size_t)N * K2 * sizeof(ST_double));
    m->wb_X1 = malloc((size_t)N * (size_t)(K1 + Kp > 0 ? K1 + Kp : 1) * sizeof(ST_double));
    m->wb_X2 = malloc((size_t)N * (size_t)(Kz2 > 0 ? Kz2 : 1) * sizeof(ST_double));
    ST_double *gs = malloc((size_t)N * (size_t)(Kp > 0 ? Kp : 1) * sizeof(ST_double));
    if (!m->wb_y || !m->wb_Y2 || !m->wb_X1 || !m->wb_X2 || !gs)
    {
        free(gs);
        civb_wild_free(m);
        return 920;
    }
    memcpy(m->wb_y, t_y, (size_t)N * sizeof(ST_double));
    for (ST_int e = 0, q = 0; e < Ke; e++)
        if (!m->collin_x[K_x + e])
            memcpy(m->wb_Y2 + (size_t)(q++) * N, t_endog + (size_t)e * N, (size_t)N * sizeof(ST_double));
    ST_int n1 = 0;
    for (ST_int k = 0, ip = 0; k < Kx; k++)
    {
        if (m->partial && m->partial[k])
            continue;
        if (!m->collin_x[ip++])
            memcpy(m->wb_X1 + (size_t)(n1++) * N, t_exog + (size_t)k * N, (size_t)N * sizeof(ST_double));
    }
    /* partial(): a column is kept when it is not absorbed by the fixed
     * effects and not spanned by the partial() columns kept before it
     * (weighted Gram-Schmidt on copies) */
    ST_int ng = 0;
    for (ST_int k = 0; k < Kx; k++)
    {
        if (!(m->partial && m->partial[k]))
            continue;
        const ST_double *x = t_exog + (size_t)k * N;
        ST_double *r = gs + (size_t)ng * N;
        memcpy(r, x, (size_t)N * sizeof(ST_double));
        const ST_double ss_raw = civ_wss(m->X_exog + (size_t)k * N, m->weights, N);
        const ST_double ss = civ_wss(x, m->weights, N);
        if (!(ss > CIV_FE_NOISE * ss_raw) || !(ss > 0.0))
            continue;
        for (ST_int q = 0; q < ng; q++)
        {
            const ST_double *g = gs + (size_t)q * N;
            ST_double num = 0, den = 0;
            for (ST_int i = 0; i < N; i++)
            {
                const ST_double w = m->weights ? m->weights[i] : 1.0;
                num += w * g[i] * r[i];
                den += w * g[i] * g[i];
            }
            const ST_double f = num / den;
            for (ST_int i = 0; i < N; i++)
                r[i] -= f * g[i];
        }
        if (!(civ_wss(r, m->weights, N) > 1e-9 * ss))
            continue;
        memcpy(m->wb_X1 + (size_t)(n1++) * N, x, (size_t)N * sizeof(ST_double));
        ng++;
    }
    free(gs);
    for (ST_int j = 0, q = 0; j < Kz; j++)
        if (!m->collin_z[j])
            memcpy(m->wb_X2 + (size_t)(q++) * N, t_excl + (size_t)j * N, (size_t)N * sizeof(ST_double));

    /* Tested coefficients: positions in e(b) = [endogenous, exogenous] */
    char tmat[] = "__civreghdfe_wb_test";
    const ST_int T = SF_col(tmat);
    if (T < 1 || !(m->wb_test = malloc((size_t)T * sizeof(ST_int))))
    {
        civb_wild_free(m);
        return T < 1 ? 198 : 920;
    }
    for (ST_int j = 0; j < T && !rc; j++)
    {
        ST_double x;
        if ((rc = SF_mat_el(tmat, 1, j + 1, &x)))
            break;
        const int p = (int)x - 1;
        int omitted = p < 0 || p >= Ke + K_x;
        ST_int idx = 0;
        if (!omitted && p < Ke)
        {
            omitted = m->collin_x[K_x + p] != 0;
            for (ST_int e = 0; e < p; e++)
                idx += !m->collin_x[K_x + e];
        }
        else if (!omitted)
        {
            const ST_int ip = p - Ke;
            omitted = m->collin_x[ip] != 0;
            idx = K2;
            for (ST_int q = 0; q < ip; q++)
                idx += !m->collin_x[q];
        }
        if (omitted)
        {
            SF_error("civreghdfe: wildbootstrap() cannot test an omitted coefficient\n");
            rc = 198;
            break;
        }
        m->wb_test[j] = idx;
    }

    /* Clusters */
    ST_int G = 0, C = 0;
    if (!rc && has_vce)
    {
        if (!(m->wb_vce = malloc((size_t)N * sizeof(ST_int))))
            rc = 920;
        else
            rc = civb_dense_ids(m, 1, m->wb_vce, &G);
    }
    if (!rc && (has_boot || has_vce))
    {
        if (!(m->wb_boot = malloc((size_t)N * sizeof(ST_int))))
            rc = 920;
        else if (has_boot)
            rc = civb_dense_ids(m, 1 + has_vce, m->wb_boot, &C);
        else
        {
            memcpy(m->wb_boot, m->wb_vce, (size_t)N * sizeof(ST_int));
            C = G;
        }
    }
    if (!m->wb_boot)
        C = N;

    /* Every absorbed effect nested in the bootstrap clusters: the
     * perturbed basis needs no projection */
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
            const ST_int l = m->levels[g][i] - 1, c = m->wb_boot ? m->wb_boot[i] : i;
            if (seen[l] < 0)
                seen[l] = c;
            else if (seen[l] != c)
                nested = 0;
        }
        free(seen);
    }
    if (rc)
    {
        civb_wild_free(m);
        return rc;
    }
    const int projected = force_projected || (!nested && !partialled);
    const int path = projected ? (nested ? 3 : 1) : (nested ? 0 : 2);

    civ_wre_input *in = &m->wb_in;
    memset(in, 0, sizeof(*in));
    in->N = N;
    in->K2 = K2;
    in->K1 = n1;
    in->Kz = Kz2;
    in->y = m->wb_y;
    in->Y2 = m->wb_Y2;
    in->X1 = m->wb_X1;
    in->X2 = m->wb_X2;
    in->w = m->weights;
    in->est = m->est_method;
    in->kclass = m->kclass_user;
    in->fuller_adj =
        m->est_method == 2 ? m->fuller_alpha / ((ST_double)m->N_eff - (ST_double)m->K_iv_fit) : 0.0;
    in->vce = vcetype;
    in->vce_id = m->wb_vce;
    in->G = G;
    in->boot_id = m->wb_boot;
    in->C = C;
    in->n_test = T;
    in->test = m->wb_test;
    in->null_imposed = !ctools_parse_bool_option(args, "nonull");
    in->projected = projected;
    m->wb_ready = 1;
    rc = ctools_scal_save("__civreghdfe_wb_C", (ST_double)C);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_wb_path", (ST_double)path);
    return rc;
}

/* The projected path's partialling: the full sample's factors with the
 * fit's weights, built once */
static ST_retcode civb_wild_projection(civ_boot_model *m)
{
    const ST_int N = m->N, G = m->G;
    HDFE_State *S = &m->wb_S;
    memset(S, 0, sizeof(*S));
    ST_int *mask = malloc((size_t)N * sizeof(ST_int));
    ST_int *level_max = malloc((size_t)(G > 0 ? G : 1) * sizeof(ST_int));
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
    free(mask);
    free(level_max);
    if (rc)
    {
        free(omega);
        if (S->factors)
            ctools_hdfe_state_cleanup(S);
        memset(S, 0, sizeof(*S));
        return rc;
    }
    S->has_weights = 1;
    S->weight_type = 2;
    S->maxiter = m->maxiter;
    S->tolerance = m->tolerance;
    S->num_threads = ctools_get_max_threads();
    S->factors_initialized = 1;
    S->weights = omega;
    S->sum_weights = sum;
    m->wb_S_ready = 1;
    if (ctools_hdfe_alloc_buffers(S, 0, 1 << 20) != 0)
        return 920;
    return 0;
}

static ST_retcode civb_wild_project(void *ctx, ST_double *cols, ST_int ncols)
{
    civ_boot_model *m = ctx;
    HDFE_SolveResult r = partial_out_columns(&m->wb_S, cols, m->N, ncols, m->wb_S.num_threads);
    if (r.status == 430)
        SF_error("civreghdfe: the wild bootstrap's partialling did not converge; increase "
                 "maxiter()\n");
    return r.status;
}

ST_retcode civreghdfe_wild_chunk(const char *args)
{
    ST_retcode rc = civb_check_retained();
    if (rc)
        return rc;
    civ_boot_model *m = g_civboot;
    int first = -1, count = -1, total = -1;
    if (ctools_parse_int_option(args, "first", &first) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "count", &count) != CTOOLS_PARSE_OK ||
        ctools_parse_int_option(args, "total", &total) != CTOOLS_PARSE_OK || first < 0 ||
        count < 1 || total < 1 || first + count > total || !m->wb_ready)
    {
        SF_error("civreghdfe: invalid wild_chunk arguments\n");
        return 198;
    }
    const int enumerate = ctools_parse_bool_option(args, "enum");
    if (first == 0)
    {
        civ_wre_free(m->wre);
        m->wre = NULL;
        rc = civ_wre_setup(&m->wb_in, total, &m->wre);
        if (rc == 198)
            SF_error("civreghdfe: wildbootstrap() could not set up the model (an underidentified "
                     "or singular instrument set)\n");
        if (rc)
            return rc;
    }
    if (!m->wre || civ_wre_B(m->wre) != total)
        return 198;
    ST_double *v = NULL;
    const ST_int C = civ_wre_C(m->wre);
    if (!enumerate)
    {
        /* The ado packs the chunk's weights replicate after replicate */
        const size_t values = (size_t)count * (size_t)C;
        if (!civb_packed_vars_ok(values))
            return 198;
        const ST_int nobs = SF_nobs();
        if (!(v = malloc(values * sizeof(ST_double))))
            return 920;
        for (size_t q = 0; q < values && !rc; q++)
            rc = civb_packed_value(nobs, q, &v[q]);
        if (rc)
        {
            free(v);
            return rc;
        }
    }
    rc = civ_wre_weights(m->wre, v, first, count, civb_wild_budget());
    free(v);
    if (rc == 909)
        SF_error("civreghdfe: not enough memory for the wild bootstrap's weights\n");
    return rc;
}

/* p(b0) at a batch of points for hdfe_wild_invert_batch */
typedef struct
{
    civ_boot_model *m;
    ST_int j;
    int ptype, null_imposed;
    ST_int R;
    double *bk, *var, *t; /* npts x R, npts x R, R */
    int cap;
    double *bk0, *var0;   /* nonull: the statistic's pieces, which do not depend on b0 */
    ST_retcode rc;
} civb_wpctx;

static ST_retcode civb_wild_eval(civb_wpctx *c, int n, const double *b0)
{
    if (n > c->cap)
    {
        free(c->bk);
        free(c->var);
        c->bk = malloc((size_t)n * c->R * sizeof(double));
        c->var = malloc((size_t)n * c->R * sizeof(double));
        if (!c->bk || !c->var)
        {
            c->cap = 0;
            return 920;
        }
        c->cap = n;
    }
    return civ_wre_eval(c->m->wre, c->j, n, b0, c->bk, c->var, ctools_get_max_threads(),
                        c->m->wb_in.projected ? civb_wild_project : NULL, c->m,
                        civb_wild_budget());
}

/* t[0 .. R) at b0 from point k of the last evaluation (or the nonull pieces) */
static void civb_wild_tvec(const civb_wpctx *c, int k, double b0)
{
    const ST_int R = c->R;
    const double *bk = c->null_imposed ? c->bk + (size_t)k * R : c->bk0;
    const double *var = c->null_imposed ? c->var + (size_t)k * R : c->var0;
    for (ST_int r = 0; r < R; r++)
    {
        const double centre = (c->null_imposed || r == 0) ? b0 : bk[0];
        c->t[r] = (var[r] > 0.0 && isfinite(var[r]) && isfinite(bk[r]))
                      ? (bk[r] - centre) / sqrt(var[r])
                      : NAN;
    }
}

static void civb_wild_pbatch(void *ctx, const double *b0, const int *from, const int *to, int n,
                             double *pv)
{
    civb_wpctx *c = ctx;
    (void)from;
    (void)to;
    for (int k = 0; k < n; k++)
        pv[k] = NAN;
    if (c->rc)
        return;
    if (c->null_imposed && (c->rc = civb_wild_eval(c, n, b0)))
        return;
    for (int k = 0; k < n; k++)
    {
        civb_wild_tvec(c, k, b0[k]);
        pv[k] = hdfe_wild_pvalue_t(c->t, c->R - 1, c->ptype, NULL);
    }
}

static ST_double civb_stata_value(double x)
{
    return isfinite(x) ? x : SV_missval;
}

ST_retcode civreghdfe_wild_finish(const char *args)
{
    ST_retcode rc = civb_check_retained();
    if (rc)
        return rc;
    civ_boot_model *m = g_civboot;
    if (!m->wre || !m->wb_ready)
        return 198;
    int ptype = 0;
    if (ctools_parse_int_option(args, "ptype", &ptype) != CTOOLS_PARSE_OK || ptype < 0 ||
        ptype > 3)
        ptype = 0;
    const double level = ctools_parse_double_option(args, "level", 95.0);
    const int ci = ctools_parse_bool_option(args, "ci");
    if (m->wb_in.projected && !m->wb_S_ready && (rc = civb_wild_projection(m)))
        return rc;
    rc = civ_wre_prepare(m->wre, ctools_get_max_threads(),
                         m->wb_in.projected ? civb_wild_project : NULL, m, civb_wild_budget());
    if (rc == 909)
        SF_error("civreghdfe: not enough memory for one wild bootstrap replicate\n");
    if (rc)
        return rc;
    const ST_int B = civ_wre_B(m->wre), R = B + 1;
    ST_int min_feasible = B;
    civb_wpctx c;
    memset(&c, 0, sizeof(c));
    c.m = m;
    c.ptype = ptype;
    c.null_imposed = m->wb_in.null_imposed;
    c.R = R;
    if (!(c.t = malloc((size_t)R * sizeof(double))))
        return 920;
    for (ST_int j = 0; j < m->wb_in.n_test && !rc; j++)
    {
        c.j = j;
        c.rc = 0;
        const double zero = 0.0;
        if ((rc = civb_wild_eval(&c, 1, &zero)))
            break;
        if (!c.null_imposed)
        {
            free(c.bk0);
            free(c.var0);
            c.bk0 = malloc((size_t)R * sizeof(double));
            c.var0 = malloc((size_t)R * sizeof(double));
            if (!c.bk0 || !c.var0)
            {
                rc = 920;
                break;
            }
            memcpy(c.bk0, c.bk, (size_t)R * sizeof(double));
            memcpy(c.var0, c.var, (size_t)R * sizeof(double));
        }
        civb_wild_tvec(&c, 0, 0.0);
        ST_int F = 0;
        const double p = hdfe_wild_pvalue_t(c.t, B, ptype, &F);
        if (F < min_feasible)
            min_feasible = F;
        rc = ctools_mat_store("__civreghdfe_wb_p", 1, j + 1, civb_stata_value(p));
        if (rc || !ci)
            continue;
        /* the interval around the estimate, in steps of its analytic
         * standard error */
        ST_double centre = NAN, scale = NAN;
        ST_double pos;
        if ((rc = SF_mat_el("__civreghdfe_wb_test", 1, j + 1, &pos)) ||
            (rc = SF_mat_el("__civreghdfe_wb_scale", 1, j + 1, &scale)))
            break;
        const int pi = (int)pos - 1;
        if (pi >= 0 && pi < m->K_b)
            centre = m->b_hat[pi];
        if (!isfinite(centre) || !(scale > 0.0))
        {
            /* no analytic scale to search on: no bounds, marked empty */
            rc = ctools_mat_store("__civreghdfe_wb_cistat", j + 1, 4, 1.0);
            continue;
        }
        hdfe_wild_ci set;
        hdfe_wild_invert_batch(civb_wild_pbatch, &c, centre, scale, 1.0 - level / 100.0, 16,
                               1e-6, &set);
        if ((rc = c.rc))
            break;
        const double row[7] = {set.lo,       set.hi,       set.open_lo,         set.open_hi,
                               set.disjoint, set.empty,    set.failed > 0};
        for (int k = 0; k < 2 && !rc; k++)
            rc = ctools_mat_store("__civreghdfe_wb_ci", j + 1, k + 1, civb_stata_value(row[k]));
        for (int k = 0; k < 5 && !rc; k++)
            rc = ctools_mat_store("__civreghdfe_wb_cistat", j + 1, k + 1, row[2 + k]);
    }
    free(c.t);
    free(c.bk);
    free(c.var);
    free(c.bk0);
    free(c.var0);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_wb_reps", (ST_double)B);
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_wb_misreps", (ST_double)(B - min_feasible));
    if (!rc)
        rc = ctools_scal_save("__civreghdfe_wb_cached", (ST_double)civ_wre_cached(m->wre));
    return rc;
}

/* ========================================================================
 * Entry point
 * ======================================================================== */

ST_retcode civreghdfe_main(const char *args)
{
    if (args == NULL || strlen(args) == 0)
    {
        SF_error("civreghdfe: No subcommand specified\n");
        return 198;
    }
    char op[32] = {0};
    sscanf(args, "%31s", op);
    /* vce(bootstrap) is linked here, not into the fit (civreghdfe_boot.h) */
    civreghdfe_boot_retain_hook = civreghdfe_boot_retain;
    civreghdfe_boot_release_hook = civreghdfe_boot_release;
    if (!strcmp(op, "iv_regression"))
    {
        /* retain keeps the model for the bootstrap (ctools_runtime.c holds
         * the command in phase 1 until boot_finish) */
        civreghdfe_boot_retain_requested = ctools_parse_bool_option(args, "retain");
        ST_retcode rc = civreghdfe_fit();
        civreghdfe_boot_retain_requested = 0;
        return rc;
    }
    if (!strcmp(op, "boot_chunk"))
        return civreghdfe_boot_chunk(args);
    if (!strcmp(op, "boot_finish"))
        return civreghdfe_boot_finish(args);
    if (!strcmp(op, "wild_setup"))
        return civreghdfe_wild_setup(args);
    if (!strcmp(op, "wild_chunk"))
        return civreghdfe_wild_chunk(args);
    if (!strcmp(op, "wild_finish"))
        return civreghdfe_wild_finish(args);
    if (!strcmp(op, "boot_release"))
    {
        /* The phase reset at the start of this call freed the model */
        civreghdfe_boot_release();
        return 0;
    }
    SF_error("civreghdfe: Unknown subcommand\n");
    return 198;
}
