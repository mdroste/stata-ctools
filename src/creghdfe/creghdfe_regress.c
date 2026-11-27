/*
 * creghdfe_regress.c
 *
 * Full regression command for creghdfe
 * Orchestrates HDFE init, partial out, and OLS
 * Part of the ctools Stata plugin suite
 */

/* Resolve output indices in the plugin varlist, including hidden Stata variables. */
#define SD_SAFEMODE
#include "creghdfe_regress.h"
#include "creghdfe_result.h"
#include "creghdfe_input.h"
#include "hdfe/hdfe_remap.h"
#include "hdfe/hdfe_sample.h"
#include "hdfe/hdfe_missing.h"
#include "creghdfe_hdfe.h"
#include "hdfe/hdfe_solver.h"
#include "hdfe/hdfe_fe_recover.h"
#include "hdfe/hdfe_indiv.h"
#include "../ctools_ols.h"
#include "../ctools_matrix.h"
#include "creghdfe_vce.h"
#include "creghdfe_core.h"
#include "creghdfe_boot.h" /* types only; called through creghdfe_boot_retain_hook */
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */
#include "../ctools_release.h"
#include "../ctools_threads.h"

/* The keep flags of the loaded rows are NULL when every row is kept. A step
 * that takes them as an array (slopes, the bootstrap) gets one of ones. */
static ST_retcode creghdfe_full_mask(ST_int **mask, ST_int N_orig)
{
    if (*mask) return 0;
    ST_int *all = (ST_int *)ctools_safe_malloc2((size_t)(N_orig > 0 ? N_orig : 1), sizeof(ST_int));
    if (!all) return 920;
    ST_int i;
    #pragma omp parallel for schedule(static) if(N_orig > 200000)
    for (i = 0; i < N_orig; i++) all[i] = 1;
    *mask = all;
    return 0;
}

/* One serial pass of the per-level counts (hdfe_sample.h), for a helper
 * thread: the rows of a factor, or its weights in row order */
typedef struct {
    HDFE_State *S;
    ST_int g, N_in, weights_pass;
    const ST_int *keep;
    const ST_double *weights;
} creghdfe_count_job;

static void *creghdfe_count_job_run(void *arg)
{
    const creghdfe_count_job *job = (const creghdfe_count_job *)arg;
    if (job->weights_pass) hdfe_count_weights(job->S, job->g, job->N_in, job->keep, job->weights);
    else hdfe_count_rows_serial(job->S, job->g);
    return NULL;
}

/* Values of the regression's matrix (rows times columns) from which the
 * columns are built on the helper threads while the projection is prepared
 * (STEP 5-7). Below it the two take well under a millisecond together, less
 * than handing work to the helper threads can cost when they have to be
 * started for it (the loader starts them only for larger data). */
#define CREGHDFE_AHEAD_MIN_VALUES ((size_t)1000000)

/* One column of the regression's matrix (creghdfe_core.h) as a task: its
 * mean first, from its sum in row order, when the caller has not taken the
 * sums (take_mean), then its sum of squares, its scale and the column
 * itself. Every pass is the column's own, in row order on one thread, so
 * the tasks give the same numbers on helper threads, as iterations of an
 * OpenMP loop or one after another. */
typedef struct {
    const creghdfe_columns *columns;
    const creghdfe_columns_opts *opts;
    ST_int k, take_mean;
    ST_double mean_divisor;
    ST_double *means; /* the array behind opts->means */
    ST_double *data, *raw, *tss, *stdevs, *norm2;
} creghdfe_column_job;

static void *creghdfe_column_job_run(void *arg)
{
    const creghdfe_column_job *job = (const creghdfe_column_job *)arg;
    if (job->take_mean)
        job->means[job->k] = creghdfe_column_sum(job->columns, job->k) / job->mean_divisor;
    creghdfe_column_build(job->columns, job->opts, job->k, job->data, job->raw, job->tss,
                          job->stdevs, job->norm2);
    return NULL;
}

/* The same column when a single effect is absorbed and the projection on it
 * is done as the column is written (creghdfe_core.h): the mean and the scale
 * first, which do not need the effect's counts ... */
static void creghdfe_column_job_scale(const creghdfe_column_job *job)
{
    if (job->take_mean)
        job->means[job->k] = creghdfe_column_sum(job->columns, job->k) / job->mean_divisor;
    creghdfe_column_scale(job->columns, job->opts, job->k, job->tss, job->stdevs);
}

/* ... then, with the counts in effect, the column itself, less its level
 * means. sums: scratch of the calling thread for the level sums. */
static void creghdfe_column_job_absorb(const creghdfe_column_job *job, const creghdfe_absorbed *effect,
                                       ST_double *sums)
{
    creghdfe_column_absorb(job->columns, job->opts, job->k, effect, job->stdevs, job->norm2, sums,
                           job->data);
}

/* Memory returned on a helper thread while the command works on: a loaded
 * variable that is no longer read, and what the projection had prepared.
 * Returning a block of 80 MB takes about 2 ms on Windows. Unlike
 * ctools_release_later() these do not keep the memory until the command
 * returns: the task starts at once. */
static void *creghdfe_variable_free_run(void *arg)
{
    stata_variable_free((stata_variable *)arg);
    return NULL;
}

static void *creghdfe_prepared_free_run(void *arg)
{
    partial_out_prepared_free((struct hdfe_coef_prepared *)arg);
    return NULL;
}

/* Recover the fixed effects from the FE component d = y - Xb - _cons -
 * residual (not the already orthogonal regression residual) with
 * hdfe_fe_recover (reghdfe's store_alphas normalization), reading y and X
 * back from Stata and storing each observation's level effect in variables
 * first_var .. first_var + G - 1. */
static ST_retcode store_fixed_effects(HDFE_State *S, const perm_idx_t *obs_map,
    const ST_int *mask, ST_int N_orig, ST_int K_keep, const ST_int *keep_idx,
    const ST_double *beta, const ST_double *resid, ST_int first_var,
    const creghdfe_input *in)
{
    ST_int N = S->N, G = S->G, idx = 0;
    ST_retcode rc = 0;
    ST_double *error = malloc((size_t)N * sizeof(*error));
    ST_double **alpha = calloc((size_t)G, sizeof(*alpha));
    if (!error || !alpha) { rc = 920; goto cleanup; }
    for (ST_int g = 0; g < G; g++) {
        alpha[g] = calloc((size_t)S->factors[g].num_levels, sizeof(**alpha));
        if (!alpha[g]) { rc = 920; goto cleanup; }
    }
    for (ST_int row = 0; row < N_orig; row++) {
        if (mask && !mask[row]) continue;
        ST_double y, x;
        rc = SF_vdata(1, (ST_int)obs_map[row], &y);
        if (rc) goto cleanup;
        error[idx] = y - beta[K_keep] - resid[idx];
        for (ST_int k = 0; k < K_keep; k++) {
            rc = SF_vdata(keep_idx[k] + 1, (ST_int)obs_map[row], &x);
            if (rc) goto cleanup;
            error[idx] -= beta[k] * x;
        }
        idx++;
    }
    /* With slopes, each absorbed effect's intercept and slopes are fitted
     * together (creghdfe_hdfe.c); a slope factor stores its coefficient */
    rc = in->num_slopes ? creghdfe_fe_recover_slopes(S, in, error, alpha)
                        : hdfe_fe_recover(S, error, alpha);
    if (rc == 430)
        SF_error("creghdfe: the fixed effects for savefe did not converge; try a larger iterate()\n");
    if (rc) goto cleanup;
    for (ST_int g = 0; g < G; g++) {
        idx = 0;
        for (ST_int row = 0; row < N_orig; row++) {
            if (mask && !mask[row]) continue;
            rc = SF_vstore(first_var + g, (ST_int)obs_map[row], alpha[g][S->factors[g].levels[idx++] - 1]);
            if (rc) goto cleanup;
        }
    }
cleanup:
    for (ST_int g = 0; alpha && g < G; g++) free(alpha[g]);
    free(alpha); free(error);
    return rc;
}

/* A mobility group's smallest value of the first absorbed effect */
typedef struct {
    ST_double value;
    ST_int root;
} creghdfe_group_min;

static int creghdfe_group_min_cmp(const void *a, const void *b)
{
    ST_double va = ((const creghdfe_group_min *)a)->value;
    ST_double vb = ((const creghdfe_group_min *)b)->value;
    return (va > vb) - (va < vb);
}

/* groupvar(): the mobility groups of effects ga and gb (0-based), numbered
 * as reghdfe numbers them. Union-find over the levels of both effects;
 * reghdfe's zigzag search starts each new subgraph at the next unvisited
 * level of the first effect in sorted order, so groups are numbered by the
 * smallest raw value of the first effect they contain (the compact level
 * numbering follows first appearance instead). fea_values holds that
 * effect's raw values for the N_orig loaded rows. */
static ST_retcode build_groupvar(const HDFE_State *S, ST_int ga, ST_int gb,
    const ST_double *fea_values, const ST_int *mask, ST_int N_orig, ST_int N,
    ST_int *group_assignments)
{
    const ST_int *fe1_compact = S->factors[ga].levels;
    const ST_int *fe2_compact = S->factors[gb].levels;
    const ST_int n1 = S->factors[ga].num_levels;
    if (n1 > INT_MAX - S->factors[gb].num_levels) return 920;
    ST_int total_nodes = n1 + S->factors[gb].num_levels;
    ctools_UnionFind *uf = ctools_uf_create(total_nodes);
    ST_int *root_to_group = (ST_int *)calloc(total_nodes, sizeof(ST_int));
    creghdfe_group_min *group_min = (creghdfe_group_min *)malloc(
        (size_t)total_nodes * sizeof(creghdfe_group_min));
    if (!uf || !root_to_group || !group_min) {
        ctools_uf_destroy(uf);
        free(root_to_group);
        free(group_min);
        return 920;
    }
    for (ST_int idx = 0; idx < N; idx++)
        ctools_uf_union(uf, fe1_compact[idx] - 1, n1 + fe2_compact[idx] - 1);
    /* root_to_group temporarily holds 1 + the root's group_min slot */
    ST_int num_roots = 0, idx = 0;
    for (ST_int row = 0; row < N_orig; row++) {
        if (mask && !mask[row]) continue;
        ST_int root = ctools_uf_find(uf, fe1_compact[idx++] - 1);
        if (root_to_group[root] == 0) {
            group_min[num_roots].value = fea_values[row];
            group_min[num_roots].root = root;
            root_to_group[root] = ++num_roots;
        } else if (fea_values[row] < group_min[root_to_group[root] - 1].value) {
            group_min[root_to_group[root] - 1].value = fea_values[row];
        }
    }
    qsort(group_min, (size_t)num_roots, sizeof(creghdfe_group_min), creghdfe_group_min_cmp);
    for (ST_int a = 0; a < num_roots; a++) root_to_group[group_min[a].root] = a + 1;
    for (idx = 0; idx < N; idx++)
        group_assignments[idx] = root_to_group[ctools_uf_find(uf, fe1_compact[idx] - 1)];
    ctools_uf_destroy(uf);
    free(root_to_group);
    free(group_min);
    return 0;
}

/*
 * FULLY COMBINED: HDFE init + Partial out + OLS in one shot
 * This reads ALL variables exactly once, eliminating data transfer overhead
 *
 * Expected varlist: depvar indepvars... fe1_var fe2_var ... feG_var [cluster_var]
 * Expected scalars:
 *   __creghdfe_K: number of data variables (depvar + indepvars)
 *   __creghdfe_G: number of FE groups
 *   __creghdfe_drop_singletons: whether to drop singletons
 *   __creghdfe_maxiter, __creghdfe_tolerance: CG solver parameters
 *   __creghdfe_verbose, __creghdfe_standardize
 *   __creghdfe_vce_type (0=unadjusted, 1=robust, 2=cluster)
 *   __creghdfe_compute_dof: whether to compute DOF/mobility groups
 *
 * Returns via scalars:
 *   __creghdfe_N: final number of observations
 *   __creghdfe_num_singletons: total singletons dropped
 *   __creghdfe_num_levels_1, __creghdfe_num_levels_2, ...: levels per FE
 *   __creghdfe_df_a: degrees of freedom absorbed
 *   __creghdfe_fe_redundant_1, ...: redundant coefficients per FE
 *   __creghdfe_mobility_groups: mobility groups (if G>=2)
 *   __creghdfe_ols_N, __creghdfe_K_keep, etc.: OLS results
 */
ST_retcode creghdfe_full_regression(int retain)
{
    ST_int K, G, N_orig, N, in1, in2;
    ST_int k, g, i, j, idx;
    ST_double val;
    ST_int drop_singletons, verbose, compute_dof;
    ST_int num_singletons;
    ST_int *mask = NULL;  /* 1 = keep, 0 = drop */
    FactorData *factors = NULL;
    /* group() and individual(): the individual effect of each row's group */
    HDFE_Indiv indiv;
    memset(&indiv, 0, sizeof(indiv));
    /* What the load and singleton removal learn about the factors, for the
     * steps after them (hdfe_sample.h) */
    HDFE_FactorOrder fe_order[10];
    HDFE_Singletons fe_sing;
    memset(fe_order, 0, sizeof(fe_order));
    memset(&fe_sing, 0, sizeof(fe_sing));
    /* The serial passes of the per-level counts, on helper threads from the
     * time the factors are built until count_pool is waited for (it is
     * non-NULL while any may be running) */
    creghdfe_count_job count_jobs[20];
    ST_int n_count_jobs = 0, buffers_ready = 0;
    ctools_persistent_pool *count_pool = NULL;
    /* The columns' tasks of STEP 5-7, and what the projection had prepared
     * for it while they ran */
    creghdfe_column_job *column_jobs = NULL;
    struct hdfe_coef_prepared *prepared = NULL;
    /* non-NULL once a task that returns memory may be running on it */
    ctools_persistent_pool *release_pool = NULL;
    /* 1 when the columns were written with a single effect projected out;
     * the threads' level sums while they are */
    ST_int absorbed = 0;
    ST_double *level_sums = NULL;
    char scalar_name[64];
    double t_start, t_load, t_copy, t_remap, t_singleton, t_dof, t_partial, t_ols, t_vce;
    ST_int mobility_groups = 1;
    ST_double df_a = 0;
    ST_int compute_resid = 0;
    ST_int resid_var_idx = 0;
    ST_retcode output_rc = 0;

    /* Data arrays */
    ST_double *data = NULL;  /* N x K matrix (column-major) */
    ST_double *xtx = NULL, *xtx_keep = NULL, *xty_keep = NULL;
    ST_double *beta_keep = NULL, *inv_xx_keep = NULL, *V_keep = NULL;
    ST_double *V_psd = NULL;    /* reghdfe's corrected V made PSD, when it is indefinite */
    ST_int psd_made = 0;
    ST_double *data_keep = NULL;
    ST_int *is_collinear = NULL, *keep_idx = NULL;
    ST_int num_collinear, K_keep, K_x, vcetype;
    ST_double df_a_nested;
    ST_double df_r;
    ST_double tss_within, rss;
    ST_double *means = NULL, *stdevs = NULL, *tss = NULL, *norm2 = NULL;
    ST_double *xty_all = NULL, *inv_rmcoll = NULL, *fe_raw = NULL;
    ST_int maxiter, standardize;
    ST_int num_threads;
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    perm_idx_t *obs_map = NULL;  /* Maps filtered index to 1-based Stata obs */

    ST_int *root_to_group = NULL;
    /* the loaded weight column (filtered.data owns it), until STEP 5 */
    const ST_double *weights = NULL;
    double *cluster_raw_values = NULL;
    ST_int *group_assignments = NULL;
    ST_double *data_compact = NULL;
    ST_double *means_compact = NULL;
    ST_int *compact_map = NULL;
    ST_double *weights_compact = NULL;
    const ST_double **loaded_cols = NULL; /* the loaded column of each regression column */
    ST_double *fv_level = NULL;           /* factor-variable level of each column, or -1 */
    ST_double *means_x = NULL;
    ST_double *inv_xx_x = NULL;
    ST_double *data_with_cons = NULL;
    double *cluster_values_compact = NULL;
    ST_double *resid = NULL;
    ST_double *bs_raw = NULL, *bs_tilde = NULL, *bs_ones = NULL; /* vce(bootstrap) */
    ST_int *cluster_ids_all = NULL; /* multi-way clustering: one column per dimension */
    ST_double *means_zero = NULL;   /* slopes only: the regression has no constant */
    ST_int num_clusters_dim[10] = {0};
    ST_int exact[10] = {0};


    t_start = get_time_sec();

    in1 = SF_in1();
    in2 = SF_in2();

    creghdfe_input input;
    ST_retcode input_rc = creghdfe_read_input(&input);
    if (input_rc) return input_rc;
    input.bs_retain = retain;
    K = input.K;
    G = input.G;
    drop_singletons = input.drop_singletons;
    verbose = input.verbose;
    maxiter = input.maxiter;
    standardize = input.standardize;
    vcetype = input.vcetype;
    ST_int has_weights = input.has_weights;
    ST_int weight_type = input.weight_type;
    compute_dof = input.compute_dof;
    ST_int dof_adjust_type = input.dof_adjust_type;
    ST_int dof_clusters = input.dof_clusters;
    ST_int savefe = input.savefe;
    ST_int savefe_var_idx = input.savefe_var_idx;
    ST_int compute_groupvar = input.compute_groupvar;
    ST_int groupvar_var_idx = input.groupvar_var_idx;
    compute_resid = input.compute_resid;
    resid_var_idx = input.resid_var_idx;
    ST_int use_quad = input.use_quad;
    ST_double tolerance = input.tolerance;
    df_a_nested = input.df_a_nested;

    /* With an individual effect the other effects may be none */
    if (G < (input.indiv ? 0 : 1) || G > 10) {
        SF_error("creghdfe: invalid number of FE groups (must be 1-10)\n");
        return 198;
    }
    /* Cluster variables after the FE variables: the dimensions of
     * vce(cluster) (vce type 2), or the time variable of vce(dkraay) (3) */
    const ST_int ncl = vcetype >= 2 ? input.num_cluster : 0;
    /* Every absorbed factor a slope (absorb(id#c.x) alone): nothing absorbs
     * the constant, and reghdfe then fits without one */
    const ST_int no_intercept = input.num_slopes == G && !input.indiv;
    if (ncl > 10 || (vcetype == 3 && ncl != 1)) {
        SF_error("creghdfe: invalid number of cluster variables\n");
        return 198;
    }
    if (K < 1) {
        SF_error("creghdfe: need a dependent variable (K < 1)\n");
        return 198;
    }
    if (K > 1000) {
        SF_error("creghdfe: too many variables (K > 1000)\n");
        return 198;
    }

    /* Every variable position below follows the wrapper's layout: K data
     * vars, G FE vars, [cluster], [weight], [resid], [groupvar], [G savefe
     * vars], final-sample flag, [sample marker]. Refuse any other length
     * rather than read or write the wrong variables. */
    {
        ST_int expected_nvars = K + G + ncl + (has_weights != 0) +
                                (compute_resid != 0) + (compute_groupvar != 0) +
                                (savefe ? G : 0) + 1 + input.touse +
                                (input.bs_unit_idx > 0) + (input.bs_strata_idx > 0) +
                                input.num_slopes + (input.indiv ? 3 : 0);
        if (SF_nvars() != expected_nvars) {
            char msg[160];
            snprintf(msg, sizeof(msg),
                     "creghdfe: plugin received %d variables, expected %d\n",
                     (int)SF_nvars(), (int)expected_nvars);
            SF_error(msg);
            return 198;
        }
    }

    /* Determine threads */
#ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    if (num_threads > K) num_threads = K;
    if (num_threads < 1) num_threads = 1;
#else
    num_threads = 1;
#endif

    /* ================================================================
     * STEP 1: PARALLEL data loading using ctools_data_load
     * This handles if/in filtering at load time, loading only filtered observations.
     * ================================================================ */

    /* Calculate total variables to load:
     * K numeric (depvar + indepvars) + G FE vars + ncl cluster vars + (1 weight if has_weights),
     * then the slope variables of slope factors */
    ST_int total_vars = K + G + ncl + (has_weights ? 1 : 0) + input.num_slopes + (input.indiv ? 1 : 0);

    /* Build variable indices array (1-based for Stata) */
    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("creghdfe: memory allocation failed\n");
        return 920;
    }

    /* Variables 1..K are numeric data, K+1..K+G are FE vars */
    for (i = 0; i < total_vars - input.num_slopes - (input.indiv ? 1 : 0); i++) {
        var_indices[i] = i + 1;
    }
    for (g = 0; g < G; g++) {
        if (input.slope_idx[g]) var_indices[i++] = input.slope_idx[g];
    }
    /* the group number of each row, last */
    if (input.indiv) var_indices[i++] = input.indiv_gidx;

    /* Load all variables in parallel with if/in filtering */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    int load_flags = CTOOLS_LOAD_NO_SORT_ORDER | (input.touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, load_flags);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("creghdfe: parallel data load failed\n");
        return ctools_stata_rc(load_rc);
    }

    ST_double screen_requested = 0;
    const int screen_missing = input.touse &&
        SF_scal_use("__creghdfe_screen", &screen_requested) == 0 && screen_requested == 1;
    double t_screen = 0;
    if (screen_missing) {
        double start = get_time_sec();
        ST_retcode screen_rc = hdfe_screen_missing(&filtered, SF_nvars());
        t_screen = get_time_sec() - start;
        if (screen_rc) { ctools_filtered_data_free(&filtered); return screen_rc; }
    }

    /* Get filtered observation count and obs_map */
    N_orig = (ST_int)filtered.data.nobs;
    obs_map = filtered.obs_map;

    if (N_orig <= 0) {
        ctools_filtered_data_free(&filtered);
        SF_error("creghdfe: no observations\n");
        return screen_missing ? 2000 : 198;
    }

    /* Verify all variables were loaded successfully (defensive check) */
    for (i = 0; i < total_vars; i++) {
        if (filtered.data.vars[i].nobs > 0) {
            if (filtered.data.vars[i].type == STATA_TYPE_DOUBLE && filtered.data.vars[i].data.dbl == NULL) {
                char errmsg[128];
                snprintf(errmsg, sizeof(errmsg), "creghdfe: variable %d failed to load (NULL data pointer, nobs=%zu)\n",
                         i + 1, filtered.data.vars[i].nobs);
                SF_error(errmsg);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
            if (filtered.data.vars[i].type == STATA_TYPE_STRING && filtered.data.vars[i].data.str == NULL) {
                char errmsg[128];
                snprintf(errmsg, sizeof(errmsg), "creghdfe: string variable %d failed to load (NULL data pointer, nobs=%zu)\n",
                         i + 1, filtered.data.vars[i].nobs);
                SF_error(errmsg);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
        }
    }

    t_load = get_time_sec();

    filtered.obs_map = NULL;  /* obs_map is now owned by this command. */
    /* The individuals of each row's group, read over the whole dataset */
    if (input.indiv) {
        output_rc = hdfe_indiv_load(&indiv, filtered.data.vars[total_vars - 1].data.dbl, N_orig,
                                    input.indiv_gidx, input.indiv_iidx, input.indiv_tidx,
                                    input.indiv_mean);
        if (output_rc) {
            SF_error(output_rc == 198 ? "creghdfe: group() numbers the groups wrongly\n"
                                      : "creghdfe: memory allocation failed for the individual effect\n");
            goto cleanup;
        }
    }
    /* Allocate factor structures (no hash tables needed with sort-based remapping) */
    factors = (FactorData *)calloc(G, sizeof(FactorData));

    /* The data matrix (column-major) is built in STEP 5 from the loaded
     * columns, for the retained rows only */
    means = (ST_double *)malloc(K * sizeof(ST_double));
    stdevs = (ST_double *)malloc(K * sizeof(ST_double));
    tss = (ST_double *)malloc(K * sizeof(ST_double));

    if (!factors || !means || !stdevs || !tss) {

        SF_error("creghdfe: memory allocation failed\n");
        output_rc = 920;
        goto cleanup;
    }

    /* mask, the keep flags of the loaded rows, stays NULL while every row
     * is kept: singleton removal makes it when it drops a row */
    /* The weights are read from the loaded column until STEP 5 compacts
     * them, before the loaded data are released */
    if (has_weights) weights = filtered.data.vars[K + G + ncl].data.dbl;

    t_copy = get_time_sec();

    /* FE variables to 1-based levels in order of first appearance, each in
     * parallel row blocks; the load also notes the rows alone in a level */
    ST_retcode remap_rc;
    {
        const ST_double *fe_values[10];
        for (g = 0; g < G; g++) fe_values[g] = filtered.data.vars[K + g].data.dbl;
        remap_rc = hdfe_load_factors_ordered(factors, fe_order, G, fe_values, N_orig);
    }

    /* Save cluster variable raw values BEFORE freeing filtered data (if clustering).
     * Also check if cluster variable matches any FE variable (common optimization). */

    ST_int cluster_matches_fe = -1;  /* -1 = no match, 0..G-1 = which FE it matches */

    if (vcetype >= 2) {
        /* Cluster variables are at positions K+G.. in filtered.data */
        cluster_raw_values = (double *)ctools_safe_malloc3((size_t)N_orig, (size_t)ncl, sizeof(double));
        if (cluster_raw_values) {
            /* Extract cluster values - direct copy since data is filtered */
            for (ST_int d = 0; d < ncl; d++)
                memcpy(cluster_raw_values + (size_t)d * N_orig, filtered.data.vars[K + G + d].data.dbl,
                       N_orig * sizeof(double));

            /* Check if cluster values match any FE variable (common case: vce(cluster i) with absorb(i t)).
             * One-way clustering only: dkraay numbers its periods in time order. */
            double *cluster_data = filtered.data.vars[K + G].data.dbl;
            for (g = 0; g < G && vcetype == 2 && ncl == 1; g++) {
                if (input.slope_idx[g]) continue;
                double *fe_data = filtered.data.vars[K + g].data.dbl;

                /* Quick rejection: check first, middle, and last values */
                if (cluster_data[0] != fe_data[0] ||
                    cluster_data[N_orig/2] != fe_data[N_orig/2] ||
                    cluster_data[N_orig-1] != fe_data[N_orig-1]) {
                    continue;
                }

                /* Full comparison - data is already filtered, compare directly */
                ST_int all_match = 1;
                for (idx = 0; idx < N_orig && all_match; idx++) {
                    if (cluster_data[idx] != fe_data[idx]) {
                        all_match = 0;
                    }
                }

                if (all_match) {
                    cluster_matches_fe = g;
                    break;
                }
            }
        }
    }

    /* filtered.data stays loaded until STEP 5 has compacted its columns;
     * obs_map is kept for write-back of residuals/groupvar/savefe */

    if (remap_rc) {
        SF_error("creghdfe: FE remapping failed\n");
        output_rc = remap_rc;
        goto cleanup;
    }

    t_remap = get_time_sec();

    /* ================================================================
     * STEP 2: Iteratively drop singletons using shared utility
     * ================================================================ */
    num_singletons = 0;
    N = N_orig;

    if (input.indiv) {
        /* With an individual effect, also the groups of an individual found
         * in no other group (fweights are refused there): its own removal,
         * on flags for every row */
        if ((output_rc = creghdfe_full_mask(&mask, N_orig))) goto cleanup;
        if (drop_singletons) {
            ST_int *fe_levels[10];
            for (g = 0; g < G; g++) fe_levels[g] = factors[g].levels;
            num_singletons = hdfe_indiv_singletons(fe_levels, G, &indiv, mask);
            /* Count remaining observations (a failure leaves N = 0, handled below) */
            N = 0;
            if (num_singletons >= 0) {
                for (i = 0; i < N_orig; i++) {
                    if (mask[i]) N++;
                }
            }
        }
    } else if (drop_singletons) {
        /* As in reghdfe, fweights drop only levels whose total weight is 1;
         * aweights and pweights use the unweighted rule. The removal starts
         * from the rows the load found alone in a level; the flags of the
         * rows stay with fe_sing until the factors are built */
        ST_int *fe_levels[10], fe_level_max[10];
        for (g = 0; g < G; g++) {
            fe_levels[g] = factors[g].levels;
            fe_level_max[g] = factors[g].num_levels;
        }
        if (hdfe_drop_singletons(&fe_sing, fe_levels, fe_level_max, fe_order, G, N_orig,
                                 (has_weights && weight_type == 2) ? weights : NULL, (verbose >= 1))) {
            num_singletons = -1;  /* a failure leaves N = 0, handled below */
            N = 0;
        } else {
            num_singletons = fe_sing.num_dropped;
            N = N_orig - num_singletons;
        }
    }

    /* Check if all observations were dropped as singletons */
    if (N == 0) {
        if (num_singletons >= 0) {
            if ((output_rc = ctools_scal_save("__creghdfe_N", 0.0))) goto cleanup;
            if ((output_rc = ctools_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons))) goto cleanup;
            if ((output_rc = ctools_scal_save("__creghdfe_K_keep", 0.0))) goto cleanup;
            SF_error("creghdfe: all observations are singletons\n");
        } else {
            SF_error("creghdfe: singleton removal failed (out of memory)\n");
        }

        output_rc = num_singletons >= 0 ? 2001 : 920;
        goto cleanup;  /* 2001: all singletons */
    }

    /* The individual effect on the rows kept */
    if (input.indiv) {
        if ((output_rc = hdfe_indiv_compact(&indiv, mask))) goto cleanup;
        if ((output_rc = ctools_scal_save("__creghdfe_indiv_levels", (ST_double)indiv.NI))) goto cleanup;
    }

    t_singleton = get_time_sec();

    if (compute_groupvar) {
        group_assignments = (ST_int *)calloc(N, sizeof(ST_int));
        if (!group_assignments) { output_rc = 920; goto cleanup; }
    }

    /* ================================================================
     * STEP 3-4: Set up global state for CG solver on the retained rows.
     * Redundant coefficients and mobility groups are counted below from the
     * compacted factors.
     * ================================================================ */
    cleanup_state();
    g_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));

    if (!g_state) {
        output_rc = 920;
        goto cleanup;
    }

    g_state->K = K;
    g_state->in1 = in1;
    g_state->in2 = in2;
    g_state->maxiter = maxiter;
    g_state->tolerance = tolerance;
    g_state->verbose = verbose;
    g_state->num_threads = num_threads;
    g_state->factors_initialized = 1;
    g_state->has_weights = has_weights;
    g_state->weight_type = weight_type;
    g_state->weights = NULL;  /* Will be set after compaction */
    g_state->sum_weights = 0.0;

    /* Factors: surviving levels renumbered 1..num_levels, with counts and the
     * (not yet normalized) weighted counts. The loaded levels are already in
     * order of first appearance: with no row dropped they are the factors'
     * as they are, with a few dropped they are shifted in place */
    ST_int counts_pending = 0;
    if (input.indiv) {
        ST_int *level_ptrs[10];
        ST_int level_max[10];
        for (g = 0; g < G; g++) {
            level_ptrs[g] = factors[g].levels;
            level_max[g] = factors[g].num_levels;
        }
        ST_retcode factor_rc = hdfe_build_factors(g_state, G, level_ptrs, level_max, N_orig, mask, N,
                                                  has_weights ? weights : NULL, 0);
        if (factor_rc) {
            output_rc = factor_rc;
            goto cleanup;
        }
    } else {
        ST_retcode factor_rc = hdfe_build_factors_ordered(g_state, G, factors, fe_order,
                                                          drop_singletons ? &fe_sing : NULL, N_orig,
                                                          has_weights ? weights : NULL, 0,
                                                          &counts_pending);
        if (factor_rc) {
            output_rc = factor_rc;
            goto cleanup;
        }
        /* the keep flags are this function's from here on (NULL when no
         * row was dropped) */
        mask = fe_sing.keep;
        fe_sing.keep = NULL;
    }

    /* Rows (and weights) of every level. The passes that are serial (a
     * factor with many levels in rows of no order, weights in row order)
     * run on helper threads while the redundant coefficients are counted
     * below and, without weights, while the columns are built in STEP 5-7;
     * the others take the row blocks here */
    if (counts_pending) {
        ST_int in_blocks[10];
        for (g = 0; g < G; g++) {
            in_blocks[g] = hdfe_count_rows_parallel(g_state, g, fe_order[g].rows.local);
            if (!in_blocks[g])
                count_jobs[n_count_jobs++] = (creghdfe_count_job){g_state, g, N_orig, 0, mask, NULL};
        }
        for (g = 0; g < G && has_weights; g++)
            count_jobs[n_count_jobs++] = (creghdfe_count_job){g_state, g, N_orig, 1, mask, weights};
        if (n_count_jobs > 0) {
            count_pool = ctools_get_global_pool();
            if (!count_pool ||
                ctools_persistent_pool_submit_batch(count_pool, creghdfe_count_job_run, count_jobs,
                                                    (size_t)n_count_jobs, sizeof(count_jobs[0])) != 0) {
                count_pool = NULL;
                for (ST_int job = 0; job < n_count_jobs; job++) creghdfe_count_job_run(&count_jobs[job]);
            }
        }
        for (g = 0; g < G; g++)
            if (in_blocks[g]) hdfe_count_rows(g_state, g, fe_order[g].rows.local);
    }

    /* Heterogeneous slopes: their variables follow the weight in the
     * loaded columns. They replace the counts of their factors */
    if (input.num_slopes) {
        if (count_pool) {
            ctools_persistent_pool_wait(count_pool);
            count_pool = NULL;
        }
        const ST_double *slope_cols[10] = {0};
        ST_int pos = K + G + ncl + (has_weights ? 1 : 0);
        for (g = 0; g < G; g++)
            if (input.slope_idx[g]) slope_cols[g] = filtered.data.vars[pos++].data.dbl;
        if ((output_rc = creghdfe_full_mask(&mask, N_orig))) goto cleanup;
        if ((output_rc = creghdfe_setup_slopes(g_state, &input, slope_cols, mask, N_orig,
                                               has_weights ? weights : NULL))) goto cleanup;
    }

    /* Levels left after singleton removal; the individual effect has one
     * coefficient per individual, none counted as redundant (reghdfe's
     * dof_update_individual_fe() counts none) */
    df_a = input.indiv ? indiv.NI : 0;
    for (g = 0; g < G; g++) {
        factors[g].num_levels = g_state->factors[g].num_levels;
        df_a += factors[g].num_levels;
    }

    /* Redundant coefficients of each effect, counted as reghdfe counts them
     * (dofadjustments: 0=all, 1=none, 2=firstpair, 3=pairwise). Cluster
     * nesting is not known yet; STEP 9 recounts when it finds a nested effect. */
    ST_int redundant[10] = {0};
    ST_int dof_pairs = dof_adjust_type == 1 ? 0 :
        dof_adjust_type == 2 ? CTOOLS_DOF_FIRSTPAIR : CTOOLS_DOF_PAIRWISE;
    if (compute_dof || !dof_pairs) {
        if ((output_rc = creghdfe_redundant(g_state, &input, NULL, dof_pairs, redundant, exact)))
            goto cleanup;
        for (g = 0; g < G; g++) df_a -= redundant[g];
    }

    /* Mobility groups of the first two effects, for e(N_hdfe_extended) and
     * groupvar(); one more is reported for each later effect, as before. */
    if (dof_pairs && compute_dof && G >= 2) {
        mobility_groups = redundant[1];
        if (G > 2 && dof_pairs == CTOOLS_DOF_PAIRWISE) mobility_groups += G - 2;

        /* groupvar() is filled in STEP 9, once cluster nesting is known */
    } else if (!dof_pairs || G == 1) {
        if (dof_pairs) mobility_groups = 0;
    }
    /* groupvar() numbers the groups by the raw values of an effect, which
     * are released with the loaded data in STEP 5 */
    if (group_assignments && G >= 2) {
        fe_raw = (ST_double *)ctools_safe_malloc3((size_t)N_orig, (size_t)G, sizeof(ST_double));
        if (!fe_raw) { output_rc = 920; goto cleanup; }
        for (g = 0; g < G; g++)
            memcpy(fe_raw + (size_t)g * N_orig, filtered.data.vars[K + g].data.dbl,
                   (size_t)N_orig * sizeof(ST_double));
    }
    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_redundant_%d", g + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)redundant[g]))) goto cleanup;
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_exact_%d", g + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)exact[g]))) goto cleanup;
    }

    g_state->df_a = df_a;
    g_state->mobility_groups = mobility_groups;

    /* With weights STEP 5-7 rescale the weighted counts: all counts are
     * complete before it. Without, the passes still running are waited for
     * just before the projection, which is the first to read the counts */
    if (count_pool && has_weights) {
        ctools_persistent_pool_wait(count_pool);
        count_pool = NULL;
    }

    t_dof = get_time_sec();

    /* Allocate inv_counts, inv_weighted_counts, and thread buffers (from
     * the counts: once they are complete) */
    g_state->num_threads = num_threads;
    if (!count_pool) {
        if (ctools_hdfe_alloc_buffers(g_state, 0, K) != 0) {
            output_rc = 920;
            goto cleanup;
        }
        buffers_ready = 1;
    }

    /* ================================================================
     * STEP 5-7: The regression's columns from the loaded variables: the
     * kept rows of each, its mean, total sum of squares and scale, and the
     * standardized matrix, written once (creghdfe_core.h)
     * ================================================================ */
    data_compact = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
    means_compact = (ST_double *)ctools_safe_malloc2((size_t)K, sizeof(ST_double));
    norm2 = (ST_double *)malloc((size_t)K * sizeof(ST_double));
    loaded_cols = (const ST_double **)malloc((size_t)K * sizeof(*loaded_cols));

    if (!data_compact || !means_compact || !norm2 || !loaded_cols) {
        output_rc = 920;
        goto cleanup;
    }
    for (k = 0; k < K; k++) loaded_cols[k] = filtered.data.vars[k].data.dbl;

    /* The loaded row of every kept row; not needed when every row is kept */
    if (N < N_orig) {
        compact_map = (ST_int *)malloc(N * sizeof(ST_int));
        if (!compact_map) {
            output_rc = 920;
            goto cleanup;
        }
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                compact_map[idx] = i;
                idx++;
            }
        }
    }

    /* The weights of the kept rows (rows in parallel) */
    ST_double sum_w = 0.0;
    if (has_weights) {
        weights_compact = (ST_double *)malloc((size_t)N * sizeof(ST_double));
        if (!weights_compact) {
            output_rc = 920;
            goto cleanup;
        }
        #pragma omp parallel for schedule(static) if(N > 200000)
        for (idx = 0; idx < N; idx++) {
            weights_compact[idx] = weights[compact_map ? compact_map[idx] : idx];
        }
    }

    /* Factor-variable indicators built here instead of by fvrevar (the
     * ado's _ctools_fvnative): with __creghdfe_fvspec, regressor k (column
     * k of the loaded data) holds a factor variable when entry k of
     * __creghdfe_fvlev is a level (>= 0), and becomes 1 where the variable
     * equals it, 0 elsewhere: the doubles fvrevar's indicator holds */
    {
        ST_double use_fv = 0;
        if (SF_scal_use("__creghdfe_fvspec", &use_fv) == 0 && use_fv == 1 && K > 1) {
            fv_level = (ST_double *)malloc((size_t)K * sizeof(ST_double));
            if (!fv_level) { output_rc = 920; goto cleanup; }
            fv_level[0] = -1;
            for (k = 1; k < K; k++) {
                if (SF_mat_el("__creghdfe_fvlev", 1, k, &fv_level[k]) != 0 ||
                    !(fv_level[k] >= -1 && fv_level[k] == floor(fv_level[k]))) {
                    SF_error("creghdfe: invalid factor-variable specification\n");
                    output_rc = 198; goto cleanup;
                }
            }
            if ((output_rc = ctools_scal_save("__creghdfe_fvdone", 1.0))) goto cleanup;
        }
    }

    const creghdfe_columns columns = {
        .cols = loaded_cols, .level = fv_level, .rows = compact_map, .N = N, .K = K,
        .weights = has_weights ? weights_compact : NULL,
    };

    /* Means (weighted if using weights): each column's sum in row order.
     * With weights the sums are taken here, with the sum of the weights
     * that their normalization below needs; without, every column's task
     * takes its own (creghdfe_column_job) */
    if (has_weights) creghdfe_column_sums(&columns, means_compact, &sum_w);

    /* Divide by sum of weights (or N if unweighted) */
    ST_double mean_divisor = (has_weights) ? sum_w : (ST_double)N;
    if (!isfinite(mean_divisor) || mean_divisor <= 0.0) {
        output_rc = 498; goto cleanup;
    }
    for (k = 0; k < K && has_weights; k++) {
        means_compact[k] /= mean_divisor;
    }
    free(means);
    means = means_compact; means_compact = NULL;

    /* The compacted weights (above) to g_state */
    if (has_weights) {
        /* For aweight/pweight: normalize weights so sum(w) = N
         * This is what reghdfe does (reghdfe.mata line 3598)
         * For fweight: keep raw weights (no normalization) */
        if (weight_type == 1 || weight_type == 3) {
            const ST_double scale = (ST_double)N / sum_w;
            #pragma omp parallel for schedule(static) if(N > 200000)
            for (idx = 0; idx < N; idx++) {
                weights_compact[idx] *= scale;
            }
            /* Also normalize the weighted_counts in factors and recompute inv_weighted_counts */
            for (g = 0; g < G; g++) {
                FE_Factor *f = &g_state->factors[g];
                if (f->weighted_counts) {
                    #pragma omp parallel for schedule(static) if(f->num_levels > 200000)
                    for (i = 0; i < f->num_levels; i++) {
                        f->weighted_counts[i] *= scale;
                        /* Recompute inv_weighted_counts after scaling */
                        if (f->inv_weighted_counts)
                            f->inv_weighted_counts[i] =
                                (f->weighted_counts[i] > 0) ? 1.0 / f->weighted_counts[i] : 0.0;
                    }
                }
            }
            /* After normalization, sum_w should be N for aw/pw */
            /* But we keep the original sum_w for fweight calculations */
        }

        /* columns.weights is this array: from here on the weights the
         * regression uses */
        g_state->weights = weights_compact; weights_compact = NULL;
        g_state->sum_weights = sum_w;

    }

    /* vce(bootstrap) refits from the unstandardized, unpartialled columns */
    if (input.bs_retain) {
        bs_raw = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
        if (!bs_raw) { output_rc = 920; goto cleanup; }
    }

    /* Each column's total sum of squares about its mean (compensated, or in
     * double-double arithmetic with the quad option, to match Mata's
     * quadcross()) and standard deviation; then the column divided by it,
     * written to the matrix together with reghdfe's norm2, its uncentered
     * unweighted sum of squares in the units that are partialled, for the
     * collinear-with-FE rule. For aweight/pweight the weights have been
     * normalized so sum(w) = N. */
    {
        const creghdfe_columns_opts column_opts = {
            .means = means, .weight_type = weight_type,
            .sum_weights = has_weights ? g_state->sum_weights : 0.0,
            .standardize = standardize, .use_quad = use_quad, .uncentered = no_intercept,
        };
        column_jobs = (creghdfe_column_job *)malloc((size_t)K * sizeof(*column_jobs));
        if (!column_jobs) { output_rc = 920; goto cleanup; }
        for (k = 0; k < K; k++) {
            column_jobs[k] = (creghdfe_column_job){
                .columns = &columns, .opts = &column_opts, .k = k, .take_mean = !has_weights,
                .mean_divisor = mean_divisor, .means = means, .data = data_compact, .raw = bs_raw,
                .tss = tss, .stdevs = stdevs, .norm2 = norm2,
            };
        }
        /* A column is one thread's work, so the columns keep K threads busy
         * and leave the others idle. When the projection will be solved for
         * the effects' coefficients, what that solve can do before it sees a
         * column (hdfe_solver.h: about half of it) is done meanwhile: the
         * columns run on the helper threads, and this thread prepares the
         * solve with the OpenMP team. Neither reads what the other writes. */
        const int ahead = G >= 2 && !input.indiv && !input.num_slopes && !input.map &&
                          input.transform == 0 && input.acceleration == 1 &&
                          (size_t)N * (size_t)K >= CREGHDFE_AHEAD_MIN_VALUES;
        ctools_persistent_pool *column_pool = ahead ? ctools_get_global_pool() : NULL;
        if (column_pool &&
            ctools_persistent_pool_submit_batch(column_pool, creghdfe_column_job_run, column_jobs,
                                                (size_t)K, sizeof(column_jobs[0])) == 0) {
            /* The tasks read the K regression columns. The other loaded
             * variables (effects, clusters, weight, slopes) were copied or
             * turned into levels above and go back now, each as a task of
             * its own, so that the preparation's memory does not come on
             * top of theirs. If they cannot be queued they go back with
             * the rest below. */
            weights = NULL;
            if (filtered.data.nvars > (size_t)K)
                (void)ctools_persistent_pool_submit_batch(column_pool, creghdfe_variable_free_run,
                                                          &filtered.data.vars[K],
                                                          filtered.data.nvars - (size_t)K,
                                                          sizeof(filtered.data.vars[0]));
            prepared = partial_out_prepare(g_state, K);
            /* the columns, and with them the count passes still running */
            ctools_persistent_pool_wait(column_pool);
            count_pool = NULL;
        } else if (G == 1 && !input.indiv && !input.num_slopes && !input.map &&
                   input.transform == 0 && input.acceleration == 1 && !use_quad && !bs_raw &&
                   maxiter >= 1 && isfinite(tolerance) && tolerance > 0 &&
                   !fe_order[0].rows.local &&
                   !hdfe_levels_grouped(g_state->factors[0].levels, N)) {
            /* A single absorbed effect, projected by partial_out_columns()
             * as it stands, in rows that are not grouped by it: the
             * projection is each value less the mean of its level, and
             * the columns are written with it done (creghdfe_core.h),
             * which saves the matrix a reading. (The load noted whether
             * the levels keep to their own stretch of rows, as in a
             * panel sorted by its unit: such rows, and any others that
             * the scan finds grouped, go to partial_out_columns().)
             * First what does not need the effect's counts ... */
            #pragma omp parallel for schedule(static) if(K > 1 && (size_t)N * (size_t)K > 200000)
            for (k = 0; k < K; k++) creghdfe_column_job_scale(&column_jobs[k]);
            /* ... then, with the counts complete and their inverses made,
             * the columns */
            if (count_pool) {
                ctools_persistent_pool_wait(count_pool);
                count_pool = NULL;
            }
            if (!buffers_ready) {
                if (ctools_hdfe_alloc_buffers(g_state, 0, K) != 0) {
                    output_rc = 920;
                    goto cleanup;
                }
                buffers_ready = 1;
            }
            {
                /* the counts and the weights that the projection on the
                 * rows uses (hdfe_solver.c) */
                const FE_Factor *f = &g_state->factors[0];
                const int weighted = g_state->weights != NULL;
                const creghdfe_absorbed effect = {
                    .levels = f->levels, .num_levels = f->num_levels,
                    .inv = (weighted && f->inv_weighted_counts) ? f->inv_weighted_counts
                                                                : f->inv_counts,
                    .counts = (weighted && f->weighted_counts) ? f->weighted_counts : f->counts,
                };
                /* The level sums of a column are its thread's scratch.
                 * The whole team runs the loop, a column per thread: a
                 * team of K threads would cost the next region its other
                 * threads. */
                const int team = omp_get_max_threads();
                const size_t lanes = (size_t)(K < team ? K : team);
                level_sums = (ST_double *)ctools_safe_malloc3(lanes, (size_t)f->num_levels,
                                                              sizeof(ST_double));
                if (!level_sums) { output_rc = 920; goto cleanup; }
                #pragma omp parallel for schedule(static) if(K > 1 && (size_t)N * (size_t)K > 200000)
                for (k = 0; k < K; k++) {
                    creghdfe_column_job_absorb(&column_jobs[k], &effect,
                        level_sums + (size_t)omp_get_thread_num() * (size_t)f->num_levels);
                }
                free(level_sums); level_sums = NULL;
            }
            absorbed = 1;
        } else {
            #pragma omp parallel for schedule(static) if(K > 1 && (size_t)N * (size_t)K > 200000)
            for (k = 0; k < K; k++) creghdfe_column_job_run(&column_jobs[k]);
        }
        free(column_jobs); column_jobs = NULL;
    }

    /* The loaded columns are no longer needed; use the matrix */
    free(fv_level); fv_level = NULL;
    free(compact_map); compact_map = NULL;
    free((void *)loaded_cols); loaded_cols = NULL;
    weights = NULL;
    stata_data_free_threaded(&filtered.data);
    data = data_compact; data_compact = NULL;

    /* The counts of every level are complete from here on (see STEP 3-4) */
    if (count_pool) {
        ctools_persistent_pool_wait(count_pool);
        count_pool = NULL;
    }
    if (!buffers_ready) {
        if (ctools_hdfe_alloc_buffers(g_state, 0, K) != 0) {
            output_rc = 920;
            goto cleanup;
        }
        buffers_ready = 1;
    }

    /* Partial out via CG solver using shared helper */
    HDFE_SolveResult projection = {0, 1, 0};
    if (absorbed) {
        /* done as the columns were written: what partial_out_columns()
         * reports for a single effect */
        projection = ctools_cancel_requested() ? (HDFE_SolveResult){1, 0, 1}
                                               : (HDFE_SolveResult){0, 1, 1};
    } else if (input.indiv) {
        /* reghdfe projects on the individual effect with LSMR; solved here
         * a hundred times tighter than tolerance(), its atol and btol */
        projection.status = hdfe_indiv_partial_out(g_state, &indiv, data, N, K, tolerance / 100.0,
                                                   maxiter, &projection.iterations);
    } else if (input.map || input.transform != 0 || input.acceleration != 1) {
        /* reghdfe's MAP iteration, asked for by option (technique(),
         * transform(), acceleration() or preconditioner()) */
        projection = partial_out_columns_options(g_state, data, N, K, input.transform,
                                                 input.acceleration);
    } else {
        /* with what was prepared while the columns were built, if anything */
        projection = partial_out_columns_prepared(g_state, data, N, K, num_threads, prepared);
    }
    if (prepared) {
        /* its memory goes back while the regression is fitted */
        ctools_persistent_pool *pool = ctools_get_global_pool();
        if (pool && ctools_persistent_pool_submit(pool, creghdfe_prepared_free_run, prepared) == 0)
            release_pool = pool;
        else
            partial_out_prepared_free(prepared);
        prepared = NULL;
    }
    if (projection.status == 920) {
        SF_error("creghdfe: memory allocation failed in the fixed-effect projection\n");
        output_rc = 920;
        goto cleanup;
    }
    if (projection.status == 498) {
        SF_error("creghdfe: numerical breakdown in fixed-effect projection; coefficients are not reliable\n");
        SF_error("check fixed-effect connectivity and scaling; abort(0) cannot override numerical breakdown\n");
        output_rc = 498;
        goto cleanup;
    }
    /* abort(0): as reghdfe, the projection after iterate() iterations stands
     * and e(converged) is 0 */
    ST_int projection_converged = 1;
    if (projection.status == 430 && !input.abort_nonconv) {
        projection.status = 0;
        projection_converged = 0;
    }
    if ((output_rc = ctools_scal_save("__creghdfe_converged", (ST_double)projection_converged)))
        goto cleanup;
    if (projection.status) {
        if (projection.status != 1) SF_error("creghdfe: fixed-effect projection did not converge\n");

        output_rc = projection.status;
        goto cleanup;
    }
    ST_int max_iters = projection.iterations;
    /* ... and starts its projections from the partialled columns */
    if (input.bs_retain) {
        bs_tilde = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
        if (!bs_tilde) { output_rc = 920; goto cleanup; }
        memcpy(bs_tilde, data, (size_t)N * K * sizeof(ST_double));
    }

    /* Save iteration count to Stata scalar */
    if ((output_rc = ctools_scal_save("__creghdfe_iterations", (ST_double)max_iters))) goto cleanup;

    t_partial = get_time_sec();

    /* ================================================================
     * STEP 8: OLS with collinearity detection
     * ================================================================ */
    K_x = K - 1;  /* Excluding y */
    /* For fweight, use sum(weights) as effective N in df calculation (reghdfe.mata line 3594) */
    ST_double N_eff_df = (weight_type == 2) ? g_state->sum_weights : N;
    df_r = N_eff_df - df_a;

    /* Allocate for collinearity detection - cast to size_t to prevent 32-bit overflow */
    const size_t K_x_alloc = (K_x > 0) ? (size_t)K_x : 1;
    xtx = (ST_double *)malloc(K_x_alloc * K_x_alloc * sizeof(ST_double));
    xty_all = (ST_double *)malloc(K_x_alloc * sizeof(ST_double));
    inv_rmcoll = (ST_double *)malloc(K_x_alloc * K_x_alloc * sizeof(ST_double));
    is_collinear = (ST_int *)malloc(K_x_alloc * sizeof(ST_int));

    if (!xtx || !xty_all || !inv_rmcoll || !is_collinear) {
        output_rc = 920;
        goto cleanup;
    }

    memset(is_collinear, 0, K_x_alloc * sizeof(ST_int));

    /* Collinearity: reghdfe's rule for regressors collinear with the fixed
     * effects (a depvar they explain is zeroed), then reghdfe_rmcoll's
     * invsym rule (creghdfe_core.c) */
    ST_int depvar_explained = 0;
    ST_int *fe_flagged = (ST_int *)calloc(K_x_alloc, sizeof(ST_int));
    if (!fe_flagged) {
        output_rc = 920;
        goto cleanup;
    }
    output_rc = creghdfe_collinear(data, N, K, norm2, NULL, has_weights ? g_state->weights : NULL,
                                   weight_type, input.collinear_tol, fe_flagged, is_collinear,
                                   &depvar_explained, xtx, xty_all, inv_rmcoll);
    for (k = 1; k < K && !output_rc; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_collinear_%d", k);
        output_rc = ctools_scal_save(scalar_name, (ST_double)fe_flagged[k - 1]);
    }
    free(fe_flagged);
    if (output_rc) goto cleanup;
    if ((output_rc = ctools_scal_save("__creghdfe_depvar_explained", (ST_double)depvar_explained))) goto cleanup;

    /* Recount collinear */
    num_collinear = 0;
    for (k = 0; k < K_x; k++) {
        if (is_collinear[k]) num_collinear++;
    }

    K_keep = K_x - num_collinear;

    /* Store collinearity flags */
    if ((output_rc = ctools_scal_save("__creghdfe_num_collinear", (ST_double)num_collinear))) goto cleanup;
    for (k = 0; k < K_x; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_collinear_%d", k + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)is_collinear[k]))) goto cleanup;
    }

    ST_int sample_var_idx = 0;
    if (SF_scal_use("__creghdfe_sample_idx", &val) == 0)
        sample_var_idx = (ST_int)val;
    /* A wrapper that passes its sample marker as the output (it is 1 on every
     * loaded row) sets __creghdfe_sample_prefilled: only the dropped rows are
     * then written, one store per dropped row, where a separate sample
     * variable takes one per kept row and a fill in Stata.
     * __creghdfe_marked tells the wrapper that the marker was updated. */
    ST_int sample_prefilled = 0;
    if (SF_scal_use("__creghdfe_sample_prefilled", &val) == 0)
        sample_prefilled = (val == 1.0);
    if (sample_var_idx > 0 && sample_prefilled) {
        for (i = 0; i < N_orig && !output_rc && N < N_orig; i++) {
            if (!mask[i]) output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 0.0);
        }
        if (!output_rc) output_rc = ctools_scal_save("__creghdfe_marked", 1.0);
    } else if (sample_var_idx > 0) {
        for (i = 0; i < N_orig && !output_rc; i++) {
            if (!mask || mask[i]) output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 1.0);
        }
    }

    if (output_rc) goto cleanup;

    /* OLS of y on the kept regressors and the constant (creghdfe_core.c);
     * its arrays become this function's, freed at cleanup */
    {
        /* Without a constant, the regression through the origin: with zero
         * means creghdfe_ols's constant is zero and its inverse is the
         * regressors' own, bordered by the constant (cleared below) */
        if (no_intercept) {
            means_zero = (ST_double *)calloc((size_t)K, sizeof(ST_double));
            if (!means_zero) { output_rc = 920; goto cleanup; }
        }
        creghdfe_ols_in ols_in = {
            .data = data, .N = N, .K = K, .is_collinear = is_collinear,
            .xtx = xtx, .xty_all = xty_all,
            .weights = has_weights ? g_state->weights : NULL, .weight_type = weight_type,
            .sum_weights = has_weights ? g_state->sum_weights : 0.0,
            .means = no_intercept ? means_zero : means, .stdevs = stdevs,
            .standardize = standardize, .use_quad = use_quad,
        };
        creghdfe_ols_out ols;
        ST_retcode ols_rc = creghdfe_ols(&ols_in, &ols);
        keep_idx = ols.keep_idx;
        data_keep = ols.data_keep;
        xtx_keep = ols.xtx_keep;
        xty_keep = ols.xty_keep;
        beta_keep = ols.beta_keep;
        inv_xx_keep = ols.inv_xx_keep;
        means_x = ols.means_x;
        tss_within = ols.tss_within;
        rss = ols.rss;
        if (ols_rc == 198) SF_error("creghdfe: X'X not positive definite\n");
        if (ols_rc) {
            output_rc = ols_rc;
            goto cleanup;
        }
    }
    t_ols = get_time_sec();

    ST_int K_with_cons = K_keep + 1;
    /* No constant: its row and column of the bread are zero, so every VCE
     * leaves the regressors' block as without the constant and the
     * constant's variance zero (its column of X is then irrelevant) */
    if (no_intercept) {
        for (i = 0; i < K_with_cons; i++) {
            inv_xx_keep[(size_t)i * K_with_cons + K_keep] = 0.0;
            inv_xx_keep[(size_t)K_keep * K_with_cons + i] = 0.0;
        }
        beta_keep[K_keep] = 0.0;
    }
    V_keep = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    V_psd = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    if (!V_keep || !V_psd) {
        output_rc = 920;
        goto cleanup;
    }
    ST_double xb_mean = 0.0;

    /* Build data with constant for the robust, cluster and dkraay VCEs (the
     * unadjusted VCE does not use it) - cast to size_t to prevent 32-bit overflow */
    if (vcetype >= 1) {
        data_with_cons = (ST_double *)malloc((size_t)N * (K_with_cons + 1) * sizeof(ST_double));
        if (!data_with_cons) {
            output_rc = 920;
            goto cleanup;
        }

        /* y (partialled), X with the means added back, and the constant column */
        #pragma omp parallel for schedule(static) if((size_t)N * (size_t)K_with_cons > 200000)
        for (k = 0; k <= K_with_cons; k++) {
            ST_double *dst = data_with_cons + (size_t)k * N;
            if (k == 0) {
                memcpy(dst, data_keep, (size_t)N * sizeof(ST_double));
            } else if (k == K_with_cons) {
                for (ST_int row = 0; row < N; row++) dst[row] = 1.0;
            } else {
                const ST_double *src = data_keep + (size_t)k * N;
                const ST_double mean_k = means_x[k - 1];
                for (ST_int row = 0; row < N; row++) dst[row] = src[row] + mean_k;
            }
        }
    }

    /* ================================================================
     * STEP 9: Compute VCE
     * ================================================================ */

    /* Read cluster variable if clustering */
    ST_double df_a_nested_computed = 0;  /* Track FEs nested within cluster */
    ST_int nested[10] = {0};

    if (vcetype >= 2 && ncl > 1) {
        /* Multi-way clustering: every dimension numbered in sorted order;
         * an effect is nested when it is nested in any of them */
        cluster_ids_all = (ST_int *)ctools_safe_malloc3((size_t)N, (size_t)ncl, sizeof(ST_int));
        cluster_values_compact = (double *)malloc(N * sizeof(double));
        if (!cluster_ids_all || !cluster_values_compact || !cluster_raw_values) {
            output_rc = 920;
            goto cleanup;
        }
        num_clusters = 0;
        for (ST_int d = 0; d < ncl; d++) {
            idx = 0;
            for (i = 0; i < N_orig; i++)
                if (!mask || mask[i]) cluster_values_compact[idx++] = cluster_raw_values[(size_t)d * N_orig + i];
            if (ctools_numeric_to_cluster_ids(cluster_values_compact, N, cluster_ids_all + (size_t)d * N,
                                              &num_clusters_dim[d]) != 0) {
                output_rc = 920;
                goto cleanup;
            }
            if (d == 0 || num_clusters_dim[d] < num_clusters) num_clusters = num_clusters_dim[d];
            for (g = 0; g < G && dof_clusters; g++) {
                if (input.slope_idx[g] || nested[g]) continue;
                ST_int is_nested = ctools_fe_nested_in_cluster(
                    g_state->factors[g].levels, g_state->factors[g].num_levels,
                    cluster_ids_all + (size_t)d * N, N);
                if (is_nested < 0) { output_rc = 920; goto cleanup; }
                nested[g] = is_nested;
            }
        }
        free(cluster_values_compact); cluster_values_compact = NULL;
        for (g = 0; g < G; g++) {
            if (nested[g]) df_a_nested_computed += g_state->factors[g].num_levels;
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
            if ((output_rc = ctools_scal_save(scalar_name, (ST_double)nested[g]))) goto cleanup;
        }
        free(cluster_raw_values);
        cluster_raw_values = NULL;
    } else if (vcetype >= 2) {
        cluster_ids = (ST_int *)malloc(N * sizeof(ST_int));
        if (!cluster_ids) {
            output_rc = 920;
            goto cleanup;
        }

        /* FAST PATH: If cluster variable matches an FE variable (e.g., vce(cluster i) with absorb(i t)),
         * we can reuse the existing FE levels as cluster IDs - no hashing or SF_vdata needed! */
        if (cluster_matches_fe >= 0) {
            /* Use existing FE levels (converted to 0-based) as cluster_ids */
            ST_int *fe_levels_for_cluster = g_state->factors[cluster_matches_fe].levels;
            num_clusters = g_state->factors[cluster_matches_fe].num_levels;

            for (idx = 0; idx < N; idx++) {
                cluster_ids[idx] = fe_levels_for_cluster[idx] - 1;  /* Convert 1-based to 0-based */
            }

            /* Check FE nesting within cluster */
            for (g = 0; g < G; g++) {
                ST_int is_nested;
                if (!dof_clusters || input.slope_idx[g]) {
                    is_nested = 0;
                } else if (g == cluster_matches_fe) {
                    is_nested = 1;  /* FE is the cluster variable - trivially nested */
                } else {
                    is_nested = ctools_fe_nested_in_cluster(
                        g_state->factors[g].levels, g_state->factors[g].num_levels,
                        cluster_ids, N);
                    if (is_nested < 0) { output_rc = 920; goto cleanup; }
                }
                nested[g] = is_nested;
                if (is_nested) df_a_nested_computed += g_state->factors[g].num_levels;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)is_nested))) goto cleanup;
            }
        } else {
            /* Cluster variable is different from all FE variables.
             * Use saved cluster_raw_values and remap with sort (still much faster than SF_vdata + hash). */
            if (!cluster_raw_values) {
                output_rc = 920;
                goto cleanup;
            }

            /* Extract cluster values for non-singleton observations and remap using sort */
            cluster_values_compact = (double *)malloc(N * sizeof(double));

            if (!cluster_values_compact) {
                output_rc = 920;
                goto cleanup;
            }

            /* Compact cluster values (remove singletons) */
            idx = 0;
            for (i = 0; i < N_orig; i++) {
                if (!mask || mask[i]) {
                    cluster_values_compact[idx] = cluster_raw_values[i];
                    idx++;
                }
            }

            /* Dense 0-based ids in sorted value order */
            if (ctools_numeric_to_cluster_ids(cluster_values_compact, N, cluster_ids, &num_clusters) != 0) {
                output_rc = 920;
                goto cleanup;
            }

            free(cluster_values_compact); cluster_values_compact = NULL;

            /* Check if any FE is nested within cluster */
            for (g = 0; g < G; g++) {
                ST_int is_nested = (!dof_clusters || input.slope_idx[g]) ? 0 : ctools_fe_nested_in_cluster(
                    g_state->factors[g].levels, g_state->factors[g].num_levels,
                    cluster_ids, N);
                if (is_nested < 0) { output_rc = 920; goto cleanup; }
                nested[g] = is_nested;
                if (is_nested) df_a_nested_computed += g_state->factors[g].num_levels;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)is_nested))) goto cleanup;
            }
        }

        /* Free saved cluster raw values */
        free(cluster_raw_values);
        cluster_raw_values = NULL;
        num_clusters_dim[0] = num_clusters;
    } else {
        /* No clustering - save 0 for all FE nested status */
        for (g = 0; g < G; g++) {
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
            if ((output_rc = ctools_scal_save(scalar_name, 0.0))) goto cleanup;
        }
    }

    /* Use computed df_a_nested if we found nested FEs, otherwise use passed value */
    if (df_a_nested_computed > 0) {
        df_a_nested = df_a_nested_computed;

        /* reghdfe sets nested effects aside before it counts the redundant
         * coefficients of the others, so the count depends on the nesting.
         * df_a keeps the nested levels; they are removed through df_a_nested. */
        if (compute_dof || !dof_pairs) {
            if ((output_rc = creghdfe_redundant(g_state, &input, nested, dof_pairs, redundant, exact)))
                goto cleanup;
            df_a = input.indiv ? indiv.NI : 0;
            for (g = 0; g < G; g++) {
                df_a += g_state->factors[g].num_levels - (nested[g] ? 0 : redundant[g]);
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_redundant_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)redundant[g]))) goto cleanup;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_exact_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)exact[g]))) goto cleanup;
            }
            df_r = N_eff_df - df_a;
            g_state->df_a = df_a;
        }
    }

    df_r -= K_keep;

    /* groupvar(), as reghdfe creates it: the mobility groups of the first two
     * effects that are not nested within the clusters, when degrees of
     * freedom are adjusted for mobility groups (not dofadjustments(none) or
     * clusters alone). With one effect every observation is in group 1, as
     * documented. Otherwise no variable is created (the wrapper drops it). */
    if (group_assignments) {
        ST_int ga = -1, gb = -1, created = 0;
        if (G == 1) {
            for (idx = 0; idx < N; idx++) group_assignments[idx] = 1;
            created = 1;
            ga = 0;
        } else if (dof_pairs && compute_dof) {
            for (g = 0; g < G && gb < 0; g++) {
                if (nested[g] || input.slope_idx[g]) continue;
                if (ga < 0) ga = g;
                else gb = g;
            }
            if (gb >= 0) {
                if ((output_rc = build_groupvar(g_state, ga, gb, fe_raw + (size_t)ga * N_orig,
                                                mask, N_orig, N, group_assignments))) goto cleanup;
                created = 1;
            }
        }
        if ((output_rc = ctools_scal_save("__creghdfe_groupvar_created", (ST_double)created))) goto cleanup;
        if ((output_rc = ctools_scal_save("__creghdfe_groupvar_fe1", (ST_double)(ga + 1)))) goto cleanup;
        if ((output_rc = ctools_scal_save("__creghdfe_groupvar_fe2", (ST_double)(gb + 1)))) goto cleanup;
        if (!created) { free(group_assignments); group_assignments = NULL; }
    }
    free(fe_raw); fe_raw = NULL;

    /* Compute residuals (needed for robust/cluster VCE and optionally stored) */
    resid = (ST_double *)malloc(N * sizeof(ST_double));
    if (!resid) {
        SF_error("creghdfe: memory allocation failed for residuals\n");
        output_rc = 920;
        free(data_with_cons); data_with_cons = NULL;
        goto cleanup;
    }
    if (resid) {
        /* y - X*beta: each fitted value accumulates the products over the
         * regressors in order, in plain double arithmetic (no FMA, as Mata's
         * matrix multiply) */
        memset(resid, 0, (size_t)N * sizeof(ST_double));
        ctools_matvec_add(data_keep + N, N, K_keep, beta_keep, resid);
        #pragma omp parallel for schedule(static) if(N > 200000)
        for (idx = 0; idx < N; idx++) {
            resid[idx] = data_keep[idx] - resid[idx];
        }

        /* Recompute RSS from residuals: use quad-precision or Kahan. */
        if (use_quad) {
            if (has_weights) {
                rss = dd_sum_sq_weighted(resid, g_state->weights, N);
            } else {
                rss = dd_sum_sq(resid, N);
            }
        } else {
            rss = 0.0;
            if (has_weights) {
                for (idx = 0; idx < N; idx++) {
                    volatile ST_double sq = resid[idx] * resid[idx];
                    rss += g_state->weights[idx] * sq;
                }
            } else {
                for (idx = 0; idx < N; idx++) {
                    volatile ST_double sq = resid[idx] * resid[idx];
                    rss += sq;
                }
            }
        }
    }

    if (vcetype == 0) {
        output_rc = compute_vce_unadjusted(inv_xx_keep, rss, df_r, K_with_cons, V_keep);
    } else if (vcetype >= 1) {
        if (resid) {
            /* Pass weights for weighted VCE (NULL if no weights) */
            ST_double *vce_weights = has_weights ? g_state->weights : NULL;
            /* N_eff = sum(weights) for fweight, else N */
            ST_double N_eff = (weight_type == 2) ? g_state->sum_weights : N;
            if (vcetype == 1) {
                output_rc = compute_vce_robust(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type, N, N_eff, K_with_cons, df_a, V_keep);
            } else if (vcetype == 2 && ncl == 1) {
                ST_int df_m_cluster = K_keep;
                output_rc = compute_vce_cluster(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type, cluster_ids, N, N_eff, K_with_cons, num_clusters, V_keep, df_m_cluster, df_a, df_a_nested);
            } else if (vcetype == 2) {
                ST_int psd_fixed = 0;
                output_rc = compute_vce_multicluster(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type,
                    cluster_ids_all, num_clusters_dim, ncl, N, N_eff, K_with_cons, V_keep, K_keep, df_a,
                    df_a_nested, &psd_fixed, V_psd, &psd_made);
                if (!output_rc) output_rc = ctools_scal_save("__creghdfe_psd_fixed", (ST_double)psd_fixed);
            } else {
                ST_int bandwidth_used = 0;
                output_rc = compute_vce_dkraay(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type,
                    cluster_ids, num_clusters, input.dkraay_bw, N, N_eff, K_with_cons, V_keep, K_keep,
                    df_a - df_a_nested, &bandwidth_used, V_psd, &psd_made);
                if (!output_rc) output_rc = ctools_scal_save("__creghdfe_dkraay_bw", (ST_double)bandwidth_used);
            }
            /* The clusters of each dimension, and the fewest of them */
            for (ST_int d = 0; d < ncl && vcetype >= 2 && !output_rc; d++) {
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_N_clust_%d", d + 1);
                output_rc = ctools_scal_save(scalar_name, (ST_double)num_clusters_dim[d]);
            }
            if (vcetype >= 2 && !output_rc)
                output_rc = ctools_scal_save("__creghdfe_N_clust", (ST_double)num_clusters);
        }
    }

    free(data_with_cons); data_with_cons = NULL;
    if (output_rc) {
        SF_error("creghdfe: covariance calculation failed\n");
        free(resid); resid = NULL;
        goto cleanup;
    }

    /* Wald F of the kept regressors, as reghdfe computes it for every VCE:
     * b' invsym(V) b / df_m in the units of the standardized regression, and
     * missing when invsym() finds V rank deficient (too few clusters, a
     * perfect fit, a saturated model). Without standardize, b and V are put
     * in those units first. */
    {
        ST_double F_wald = SV_missval;
        if (K_keep > 0) {
            ST_double *Vs = (ST_double *)malloc((size_t)K_keep * K_keep * sizeof(ST_double));
            ST_double *Vinv = (ST_double *)malloc((size_t)K_keep * K_keep * sizeof(ST_double));
            ST_double *bs = (ST_double *)malloc((size_t)K_keep * sizeof(ST_double));
            if (!Vs || !Vinv || !bs) {
                free(Vs); free(Vinv); free(bs);
                output_rc = 920;
                goto cleanup;
            }
            for (i = 0; i < K_keep; i++) {
                ST_double si = standardize ? 1.0 : stdevs[keep_idx[i]] / stdevs[0];
                bs[i] = beta_keep[i] * si;
                for (j = 0; j < K_keep; j++) {
                    ST_double sj = standardize ? 1.0 : stdevs[keep_idx[j]] / stdevs[0];
                    Vs[(size_t)i * K_keep + j] = V_keep[(size_t)i * K_with_cons + j] * si * sj;
                }
            }
            ST_int rank_V = 0;
            if (ctools_invsym(Vs, K_keep, Vinv, &rank_V) == 0 && rank_V == K_keep) {
                ST_double w = 0.0;
                for (i = 0; i < K_keep; i++)
                    for (j = 0; j < K_keep; j++)
                        w += bs[i] * Vinv[(size_t)i * K_keep + j] * bs[j];
                if (isfinite(w)) F_wald = w / K_keep;
            }
            free(Vs); free(Vinv); free(bs);
        }
        if ((output_rc = ctools_scal_save("__creghdfe_F", F_wald))) {
            free(resid); resid = NULL;
            goto cleanup;
        }
    }

    /* Destandardize if needed */
    if (standardize) {
        ST_double stdev_y = stdevs[0];

        /* Destandardize residuals first */
        if (resid) {
            for (idx = 0; idx < N; idx++) {
                resid[idx] *= stdev_y;
            }
        }

        /* Match reghdfe: multiply standardized RSS/TSS_within by stdev_y^2.
         * reghdfe computes RSS from standardized residuals via quadcross(),
         * then does sol.rss = sol.rss * stdev_y ^ 2 (reghdfe.mata:3717). */
        rss = rss * stdev_y * stdev_y;
        tss_within = tss_within * stdev_y * stdev_y;

        for (k = 0; k < K_keep; k++) {
            ST_int orig_idx = keep_idx[k];
            ST_double stdev_x = stdevs[orig_idx];
            beta_keep[k] = beta_keep[k] * stdev_y / stdev_x;
        }

        for (i = 0; i < K_with_cons; i++) {
            ST_double stdev_x_i = (i < K_keep) ? (stdevs[keep_idx[i]] / stdev_y) : (1.0 / stdev_y);
            for (j = 0; j < K_with_cons; j++) {
                ST_double stdev_x_j = (j < K_keep) ? (stdevs[keep_idx[j]] / stdev_y) : (1.0 / stdev_y);
                V_keep[i * K_with_cons + j] = V_keep[i * K_with_cons + j] / (stdev_x_i * stdev_x_j);
                if (psd_made) V_psd[i * K_with_cons + j] = V_psd[i * K_with_cons + j] / (stdev_x_i * stdev_x_j);
            }
        }

        /* Re-compute constant */
        xb_mean = 0.0;
        for (k = 0; k < K_keep; k++) {
            xb_mean += means[keep_idx[k]] * beta_keep[k];
        }
        beta_keep[K_keep] = no_intercept ? 0.0 : means[0] - xb_mean;
    }

    /* Store residuals back to Stata if requested.
     * Uses obs_map to write to correct Stata observations. */
    if (compute_resid && resid_var_idx > 0 && resid) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (!mask || mask[i]) {
                if (!output_rc) output_rc = SF_vstore(resid_var_idx, (ST_int)obs_map[i], resid[idx]);
                idx++;
            }
        }
    }

    /* Store groupvar (mobility group assignments) back to Stata if requested */
    if (compute_groupvar && groupvar_var_idx > 0 && group_assignments) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (!mask || mask[i]) {
                if (!output_rc) output_rc = SF_vstore(groupvar_var_idx, (ST_int)obs_map[i], (ST_double)group_assignments[idx]);
                idx++;
            }
        }
    }

    if (savefe && savefe_var_idx > 0 && resid && !output_rc) {
        output_rc = store_fixed_effects(g_state, obs_map, mask, N_orig,
            K_keep, keep_idx, beta_keep, resid, savefe_var_idx, &input);
    }

    /* Free residuals and group assignments */
    free(resid);
    resid = NULL;
    free(group_assignments);
    group_assignments = NULL;

    if (output_rc) goto cleanup;
    t_vce = get_time_sec();

    /* ================================================================
     * STEP 10: Store all results
     * ================================================================ */

    const creghdfe_result result = {
        .weight_type = weight_type,
        .has_weights = has_weights,
        .N = N,
        .num_singletons = num_singletons,
        .G = G,
        .mobility_groups = mobility_groups,
        .K_keep = K_keep,
        .K_with_cons = K_with_cons,
        .vcetype = vcetype,
        .num_clusters = num_clusters,
        .df_a = df_a,
        .df_a_nested = df_a_nested,
        .rss = rss,
        .tss_within = tss_within,
        .tss = tss,
        .beta_keep = beta_keep,
        .V_keep = V_keep,
        .factors = factors,
        .g_state = g_state,
    };
    output_rc = creghdfe_post_results(&result);
    if (output_rc) goto cleanup;
    /* The PSD version of e(V), which the wrapper posts when Stata's ereturn
     * post would set reghdfe's indefinite one to zero (HDFE-27) */
    if ((output_rc = ctools_scal_save("__creghdfe_psd_made", (ST_double)psd_made))) goto cleanup;
    for (i = 0; i < K_with_cons && psd_made; i++)
        for (j = 0; j < K_with_cons; j++)
            if ((output_rc = ctools_mat_store("__creghdfe_V_psd", i + 1, j + 1,
                                              V_psd[(size_t)i * K_with_cons + j]))) goto cleanup;
    /* e(sumweights): the sum of the weights as given (N unweighted) */
    if ((output_rc = ctools_scal_save("__creghdfe_sumweights",
                                      has_weights ? g_state->sum_weights : (ST_double)N))) goto cleanup;

    /* vce(bootstrap): keep what the replicates need (creghdfe_boot.c) */
    if (input.bs_retain) {
        bs_ones = (ST_double *)malloc((size_t)K * sizeof(ST_double));
        if (!bs_ones) { output_rc = 920; goto cleanup; }
        for (k = 0; k < K; k++) bs_ones[k] = 1.0;
        /* the bootstrap reads the keep flags of every loaded row */
        if ((output_rc = creghdfe_full_mask(&mask, N_orig))) goto cleanup;
        const creghdfe_boot_source source = {
            .N = N, .K = K, .G = G, .raw = bs_raw, .tilde = bs_tilde,
            .stdevs = standardize ? stdevs : bs_ones,
            .weights = has_weights ? g_state->weights : NULL,
            .weight_type = has_weights ? weight_type : 0,
            .factors = g_state->factors, .is_collinear = is_collinear,
            .b_hat = beta_keep, .K_keep = K_keep,
            .maxiter = maxiter, .drop_singletons = drop_singletons, .use_quad = use_quad,
            .tolerance = tolerance, .collinear_tol = input.collinear_tol,
            .mask = mask, .obs_map = obs_map, .N_orig = N_orig,
            .unit_var = input.bs_unit_idx, .strata_var = input.bs_strata_idx,
            .in = &input,
        };
        if (!creghdfe_boot_retain_hook) { output_rc = 198; goto cleanup; }
        if ((output_rc = creghdfe_boot_retain_hook(&source))) goto cleanup;
    }

    /* Save timing results to Stata scalars */
    if ((output_rc = ctools_scal_save("_creghdfe_time_load", t_load - t_start - t_screen))) goto cleanup;
    if (screen_missing && (output_rc = ctools_scal_save("_creghdfe_time_screen", t_screen))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_copy", t_copy - t_load))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_remap", t_remap - t_copy))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_singleton", t_singleton - t_remap))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_dof", t_dof - t_singleton))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_partial", t_partial - t_dof))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_ols", t_ols - t_partial))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_vce", t_vce - t_ols))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_total", t_vce - t_start))) goto cleanup;
    CTOOLS_SAVE_THREAD_INFO("_creghdfe");

    /* ================================================================
     * Cleanup
     * ================================================================ */
cleanup:
    /* No count pass may outlive the factors it fills */
    if (count_pool) {
        ctools_persistent_pool_wait(count_pool);
        count_pool = NULL;
    }
    /* ... and no task that returns memory may outlive the command */
    if (release_pool) {
        ctools_persistent_pool_wait(release_pool);
        release_pool = NULL;
    }
    free(column_jobs);
    free(level_sums);
    partial_out_prepared_free(prepared);
    /* Everything below is returned by a helper thread after the plugin
     * call (ctools_release.h): at 10M rows that is 1-2 GB, 40 ms or more
     * on Windows */
    ctools_filtered_data_release(&filtered);
    ctools_release_later_aligned(obs_map);
    for (g = 0; g < G; g++) {
        if (factors) {
            ctools_release_later(factors[g].levels);
            ctools_release_later(factors[g].counts);
        }
    }
    ctools_release_later(factors); ctools_release_later(mask);
    for (g = 0; g < 10; g++) {
        ctools_release_later(fe_order[g].rows.first_row);
        ctools_release_later(fe_order[g].rows.single_rows);
    }
    ctools_release_later(fe_sing.keep); ctools_release_later(fe_sing.dropped);
    for (g = 0; fe_sing.affected && g < fe_sing.G; g++) {
        ctools_release_later(fe_sing.affected[g].level);
        ctools_release_later(fe_sing.affected[g].left);
    }
    ctools_release_later(fe_sing.affected);
    for (g = 0; fe_sing.rows_left && g < fe_sing.G; g++) ctools_release_later(fe_sing.rows_left[g]);
    ctools_release_later(fe_sing.rows_left);
    if (data_keep == data) data_keep = NULL;  /* shared when nothing was omitted */
    ctools_release_later(data); ctools_release_later(means);
    ctools_release_later(stdevs); ctools_release_later(tss);
    ctools_release_later(data_keep); ctools_release_later(xtx_keep);
    ctools_release_later(xty_keep); ctools_release_later(beta_keep);
    ctools_release_later(inv_xx_keep); ctools_release_later(V_keep); ctools_release_later(V_psd);
    ctools_release_later(keep_idx); ctools_release_later(xtx); ctools_release_later(is_collinear);
    ctools_release_later(norm2); ctools_release_later(xty_all);
    ctools_release_later(inv_rmcoll); ctools_release_later(fe_raw);
    ctools_release_later(cluster_ids);
    ctools_release_later(cluster_ids_all);
    ctools_release_later(means_zero);
    /* weights pointed into the loaded data, released with them */
    ctools_release_later(cluster_raw_values);
    ctools_release_later(group_assignments);
    ctools_release_later(data_compact);
    ctools_release_later(means_compact);
    ctools_release_later(compact_map);
    ctools_release_later(loaded_cols);
    ctools_release_later(fv_level);
    ctools_release_later(weights_compact);
    ctools_release_later(means_x);
    ctools_release_later(inv_xx_x);
    ctools_release_later(data_with_cons);
    ctools_release_later(cluster_values_compact);
    ctools_release_later(root_to_group);
    ctools_release_later(resid);
    ctools_release_later(bs_raw); ctools_release_later(bs_tilde); ctools_release_later(bs_ones);
    hdfe_indiv_free(&indiv);
    cleanup_state_later();
    return output_rc;
}

ST_retcode do_full_regression(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    return creghdfe_full_regression(0);
}
