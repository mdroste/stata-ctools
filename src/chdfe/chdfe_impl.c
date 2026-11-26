/* Standalone within transformation. The shared HDFE state owns its factors
 * and solver buffers; this adapter owns the loaded data and compact weights. */
#include "chdfe_impl.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include "hdfe/hdfe_sample.h"
#include "hdfe/hdfe_remap.h"
#include "hdfe/hdfe_solver.h"
#include "hdfe/hdfe_missing.h"
#include "ctools_release.h"
#include <limits.h>
#include <float.h>

static int save_number(const char *name, double x)
{
    char value[64];
    snprintf(value, sizeof(value), "%.17g", x);
    return SF_macro_save((char *)name, value);
}

static int save_factor(const char *field, int g, int value)
{
    char name[64];
    snprintf(name, sizeof(name), "_chdfe_%s%d", field, g + 1);
    return save_number(name, value);
}

/* A power-of-two rescaling preserves mantissas while keeping squared norms
 * in range. This is undone before publication and does not change projections. */
static double chdfe_unit_scale(double *x, int n)
{
    double largest = 0;
    for (int i = 0; i < n; i++)
        largest = fmax(largest, fabs(x[i]));
    int exponent = 0;
    if (largest == 0)
        return 1;
    frexp(largest, &exponent);
    double scale = ldexp(1, exponent - 1);
    for (int i = 0; i < n; i++) {
        double normalized = x[i] / scale;
        /* The solver forms squared norms. Do not erase nonzero inputs or
         * entire groups when their dynamic range exceeds those norms. */
        if (x[i] != 0 && fabs(normalized) < sqrt(DBL_MIN)) return NAN;
        x[i] = normalized;
    }
    return scale;
}

static void chdfe_release_state(void *arg)
{
    ctools_hdfe_state_cleanup((HDFE_State *)arg);
    free(arg);
}

ST_retcode chdfe_main(const char *args)
{
    int k, g, c, weighted, keepsingletons, maxiter, flags, pool, transform, accel, groupout;
    const char *cursor = args;
    if (!args || ctools_parse_next_int(&cursor, &k) || ctools_parse_next_int(&cursor, &g) ||
        ctools_parse_next_int(&cursor, &c) || ctools_parse_next_int(&cursor, &weighted) ||
        ctools_parse_next_int(&cursor, &keepsingletons) ||
        ctools_parse_next_int(&cursor, &maxiter) || ctools_parse_next_int(&cursor, &flags) ||
        ctools_parse_next_int(&cursor, &pool) || ctools_parse_next_int(&cursor, &transform) ||
        ctools_parse_next_int(&cursor, &accel) || ctools_parse_next_int(&cursor, &groupout))
        return 198;
    if (k < 1 || g < 1 || c < 0 || k > SF_nvars() || g > SF_nvars() || c > SF_nvars() ||
        weighted < 0 || weighted > 3 || keepsingletons < 0 || keepsingletons > 1 || maxiter < 1 ||
        pool < 1 || flags < 0 || flags > 15 || transform < 0 || transform > 3 || accel < 0 ||
        accel > 4 || groupout < 0 || groupout > 1)
        return 198;
    double tolerance = ctools_parse_double_option(cursor, "tol", 1e-8);
    if (!isfinite(tolerance) || tolerance < 1e-16 || tolerance >= 1)
        return 198;
    const int rowspace = ctools_parse_bool_option(cursor, "rowspace");
    /* fillall: every output (and sample()/groupvar()) is written on every
     * observation, missing (or 0 for sample()) outside the sample, so the
     * wrapper creates them without filling; _chdfe_filled acknowledges it.
     * screen: rows missing in any loaded variable are dropped here (the
     * marker cleared) instead of by the wrapper's markouts. */
    const int fillall = ctools_parse_bool_option(cursor, "fillall");
    const int screen = ctools_parse_bool_option(cursor, "screen");
    ST_retcode rc = 920;
    HDFE_State state = {0};
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    HDFE_LoadedFactor *raw = calloc(g, sizeof(*raw));
    /* What the load and singleton removal learn about the factors, for the
     * steps after them (hdfe_sample.h) */
    HDFE_FactorOrder *order = calloc(g, sizeof(*order));
    HDFE_Singletons sing = {0};
    int *slopes = calloc(g, sizeof(*slopes)), **levels = calloc(g, sizeof(*levels));
    int *sizes = calloc(g, sizeof(*sizes)), *keep = NULL, *columns = NULL;
    int *nested = calloc(g, sizeof(*nested)), *redundant = calloc(g, sizeof(*redundant));
    double **values = calloc(g, sizeof(*values)), *weights = NULL, *data = NULL, *scales = NULL;
    perm_idx_t *obs = NULL;
    FE_Factor *intercepts = calloc(g, sizeof(*intercepts));
    int *imap = calloc(g, sizeof(*imap)), *inested = calloc(g, sizeof(*inested));
    int *ired = calloc(g, sizeof(*ired)), *clusters = NULL, *mobility = NULL;
    const double **store_cols = NULL;
    int *store_vars = NULL;
    double *store_mobility = NULL;
    ctools_UnionFind *uf = NULL;
    if (!raw || !order || !slopes || !levels || !sizes || !nested || !redundant || !values ||
        !intercepts || !imap || !inested || !ired)
        goto cleanup;
    int ns = 0;
    for (int j = 0; j < g; j++)
    {
        if (ctools_parse_next_int(&cursor, &slopes[j]) || slopes[j] < 0 || slopes[j] > 2)
        {
            rc = 198;
            goto cleanup;
        }
        ns += slopes[j] != 0; /* 2 means center slope within its category */
    }
    if ((size_t)k * 2 + g + ns + c + (weighted != 0) + groupout + 2 != (size_t)SF_nvars())
    {
        rc = 198;
        goto cleanup;
    }
    int inputs = k + g + ns + c + (weighted != 0);
    int sampleout = inputs + k + 1;
    if (SF_nvars() != sampleout + groupout + 1)
    {
        rc = 198;
        goto cleanup;
    }
    for (int j = 1; j <= SF_nvars(); j++)
    {
        if (SF_var_is_string(j))
        {
            rc = 109;
            goto cleanup;
        }
    }
    columns = malloc((size_t)inputs * sizeof(*columns));
    if (!columns)
        goto cleanup;
    for (int j = 0; j < inputs; j++)
        columns[j] = j + 1;
    double start = ctools_timer_seconds();
    rc = ctools_stata_rc(ctools_data_load(&fd, columns, inputs, 0, 0,
                                          CTOOLS_LOAD_MARKER | CTOOLS_LOAD_NO_SORT_ORDER));
    if (rc)
        goto cleanup;
    if (screen && (rc = hdfe_screen_missing(&fd, SF_nvars())))
        goto cleanup;
    int nin = (int)fd.data.nobs;
    if (!nin)
    {
        rc = 2000;
        goto cleanup;
    }
    /* Without screen the ado marked missing values out: validate the private
     * protocol before remapping (Stata holds no infinities or NaNs, and every
     * value at or above missval is missing). */
    if (!screen)
    {
        const double missval = SV_missval;
        int invalid = 0;
#pragma omp parallel for schedule(static) reduction(| : invalid) if (inputs > 1 && (size_t)nin * inputs > 200000)
        for (int j = 0; j < inputs; j++)
        {
            const double *column = fd.data.vars[j].data.dbl;
            int bad = 0;
            for (int i = 0; i < nin; i++)
                bad |= !(column[i] < missval && column[i] > -INFINITY);
            invalid |= bad;
        }
        if (invalid)
        {
            rc = 198;
            goto cleanup;
        }
    }
    const double *win = weighted ? fd.data.vars[inputs - 1].data.dbl : NULL;
    if (win)
        for (int i = 0; i < nin; i++)
        {
            if (win[i] <= 0 || (weighted == 2 && win[i] != floor(win[i])))
            {
                rc = 401;
                goto cleanup;
            }
        }
    if (weighted && !isfinite(chdfe_unit_scale(fd.data.vars[inputs - 1].data.dbl, nin))) {
        ctools_error("chdfe", "weight range is too wide to preserve the projection; rescale or restrict the input");
        rc = 498;
        goto cleanup;
    }
    for (int j = 0; j < g; j++)
        values[j] = fd.data.vars[k + j].data.dbl;
    /* Levels in order of first appearance, with the rows alone in a level
     * noted; singleton removal starts from those rows, and the loaded levels
     * are the factors' as they are when no row is dropped (shifted in place
     * when a few are) */
    rc = hdfe_load_factors_ordered(raw, order, g, (const double *const *)values, nin);
    if (rc)
        goto cleanup;
    for (int j = 0; j < g; j++)
    {
        levels[j] = raw[j].levels;
        sizes[j] = raw[j].num_levels;
    }
    if (!keepsingletons && hdfe_drop_singletons(&sing, levels, sizes, order, g, nin, NULL, 0))
    {
        rc = 920;
        goto cleanup;
    }
    const int dropped = sing.num_dropped;
    int n = nin - dropped;
    if (!n)
    {
        rc = 2000;
        goto cleanup;
    }
    rc = hdfe_build_factors_ordered(&state, g, raw, order, &sing, nin, win, weighted != 0, NULL);
    if (rc)
        goto cleanup;
    /* The keep flags of the loaded rows: NULL when every row is kept */
    keep = sing.keep;
    sing.keep = NULL;
    weights = weighted ? malloc((size_t)n * sizeof(*weights)) : NULL;
    obs = malloc((size_t)n * sizeof(*obs));
    if (!obs || (weighted && !weights))
    {
        rc = 920;
        goto cleanup;
    }
    if (n == nin)
    {
        memcpy(obs, fd.obs_map, (size_t)n * sizeof(*obs));
        if (weights)
            memcpy(weights, win, (size_t)n * sizeof(*weights));
    }
    else
        for (int i = 0, q = 0; i < nin; i++)
            if (keep[i])
            {
                obs[q] = fd.obs_map[i];
                if (weights)
                    weights[q] = win[i];
                q++;
            }
    state.weights = weights;
    state.has_weights = weighted != 0;
    state.tolerance = tolerance;
    state.maxiter = maxiter;
    state.num_threads = ctools_get_max_threads();
    int slope_col = k + g;
    for (int j = 0; j < g; j++)
        if (slopes[j])
        {
            FE_Factor *f = &state.factors[j];
            int l = f->num_levels;
            f->has_intercept = 0;
            f->num_slopes = 1;
            f->slope = malloc((size_t)n * sizeof(*f->slope));
            if (!f->slope)
            {
                rc = 920;
                goto cleanup;
            }
            for (int i = 0, q = 0; i < nin; i++)
                if (!keep || keep[i])
                    f->slope[q++] = fd.data.vars[slope_col].data.dbl[i];
            slope_col++;
            double slope_scale = chdfe_unit_scale(f->slope, n);
            if (!isfinite(slope_scale)) {
                ctools_error("chdfe", "slope range is too wide to preserve the projection; rescale or restrict the input");
                rc = 498;
                goto cleanup;
            }
            if (slopes[j] == 2)
            {
                f->slope_center = calloc(l, sizeof(*f->slope_center));
                double *deviation = calloc(l, sizeof(*deviation));
                if (!f->slope_center || !deviation)
                {
                    free(deviation);
                    rc = 920;
                    goto cleanup;
                }
                /* Center relative to an observed value in each group. A
                 * constant slope then centers to exactly zero even at extreme
                 * scales, rather than amplifying summation roundoff. */
                for (int i = 0; i < n; i++)
                    f->slope_center[f->levels[i] - 1] = f->slope[i];
                for (int i = 0; i < n; i++) {
                    int h = f->levels[i] - 1;
                    deviation[h] += (f->slope[i] - f->slope_center[h]) * (weights ? weights[i] : 1);
                }
                for (int h = 0; h < l; h++)
                    f->slope_center[h] += deviation[h] / (weights ? f->weighted_counts[h] : f->counts[h]);
                free(deviation);
            }
            memset(f->counts, 0, (size_t)l * sizeof(*f->counts));
            if (weights)
                memset(f->weighted_counts, 0, (size_t)l * sizeof(*f->weighted_counts));
            double *lo = malloc((size_t)l * sizeof(*lo)), *hi = malloc((size_t)l * sizeof(*hi));
            if (!lo || !hi)
            {
                free(lo);
                free(hi);
                rc = 920;
                goto cleanup;
            }
            for (int h = 0; h < l; h++)
            {
                lo[h] = INFINITY;
                hi[h] = -INFINITY;
            }
            for (int i = 0; i < n; i++)
            {
                int h = f->levels[i] - 1;
                double x = f->slope[i] - (f->slope_center ? f->slope_center[h] : 0);
                f->counts[h] += x * x;
                if (weights)
                    f->weighted_counts[h] += weights[i] * x * x;
                lo[h] = fmin(lo[h], x);
                hi[h] = fmax(hi[h], x);
            }
            for (int h = 0; weights && h < l; h++) {
                if (f->counts[h] > 0 && !(f->weighted_counts[h] > 0)) {
                    free(lo); free(hi);
                    ctools_error("chdfe", "weighted slope norm underflows; rescale or restrict the input");
                    rc = 498;
                    goto cleanup;
                }
            }
            if (flags & 8)
                for (int h = 0; h < l; h++)
                    redundant[j] += fabs(lo[h]) * slope_scale < sqrt(2.220446049250313e-16) &&
                                    fabs(hi[h]) * slope_scale < sqrt(2.220446049250313e-16);
            free(lo);
            free(hi);
        }
    /* Cluster nesting and pairwise connected components use the same C
     * helpers as creghdfe. Slopes do not enter the intercept graph. */
    clusters = malloc((size_t)n * sizeof(*clusters));
    double *cluster_values = malloc((size_t)n * sizeof(*cluster_values));
    if (!clusters || !cluster_values)
    {
        free(cluster_values);
        rc = 920;
        goto cleanup;
    }
    if (flags & 4)
        for (int h = 0; h < c; h++)
        {
            for (int i = 0, q = 0; i < nin; i++)
                if (!keep || keep[i])
                    cluster_values[q++] = fd.data.vars[k + g + ns + h].data.dbl[i];
            int nc;
            if (ctools_numeric_to_cluster_ids(cluster_values, n, clusters, &nc))
            {
                rc = 920;
                break;
            }
            for (int j = 0; j < g; j++)
                if (!slopes[j])
                {
                    int nest = ctools_fe_nested_in_cluster(
                        state.factors[j].levels, state.factors[j].num_levels, clusters, n);
                    if (nest < 0)
                    {
                        rc = 920;
                        break;
                    }
                    nested[j] |= nest;
                }
            if (rc)
                break;
        }
    free(cluster_values);
    if (rc)
        goto cleanup;
    int ni = 0;
    for (int j = 0; j < g; j++)
        if (!slopes[j])
        {
            intercepts[ni] = state.factors[j];
            imap[ni] = j;
            inested[ni++] = nested[j];
        }
    if (ni && ctools_hdfe_redundant_intercepts(intercepts, ni, n, inested, flags & 3, ired))
    {
        rc = 920;
        goto cleanup;
    }
    /* hdfe 3 starts with one redundant intercept after the first, even when
     * the first is nested. Unlike newer reghdfe it does not subtract an
     * additional intercept from the first nonnested effect. */
    for (int j = 0; j < ni; j++)
        redundant[imap[j]] = ired[j];
    int first = -1, second = -1;
    for (int j = 0; j < ni; j++)
        if (!inested[j])
        {
            if (first < 0)
                first = j;
            else if (second < 0)
                second = j;
        }
    if (first == 0)
        redundant[imap[first]] = 0;
    if (groupout && second >= 0 && (flags & 3))
    {
        FE_Factor *a = &intercepts[first], *b = &intercepts[second];
        if (a->num_levels > INT_MAX - b->num_levels)
        {
            rc = 920;
            goto cleanup;
        }
        uf = ctools_uf_create(a->num_levels + b->num_levels);
        mobility = calloc(a->num_levels, sizeof(*mobility));
        if (!uf || !mobility)
        {
            rc = 920;
            goto cleanup;
        }
        for (int i = 0; i < n; i++)
            ctools_uf_union(uf, a->levels[i] - 1, a->num_levels + b->levels[i] - 1);
        /* Components are numbered by their smallest original category code,
         * not by the solver's first-appearance numbering. */
        size_t roots = (size_t)a->num_levels + b->num_levels;
        double *root_min = malloc(roots * sizeof(*root_min));
        double *level_min = malloc((size_t)a->num_levels * sizeof(*level_min));
        if (!root_min || !level_min)
        {
            free(root_min);
            free(level_min);
            rc = 920;
            goto cleanup;
        }
        for (size_t h = 0; h < roots; h++)
            root_min[h] = INFINITY;
        for (int i = 0, q = 0; i < nin; i++)
            if (!keep || keep[i])
            {
                int root = ctools_uf_find(uf, a->levels[q++] - 1);
                root_min[root] = fmin(root_min[root], fd.data.vars[k + imap[first]].data.dbl[i]);
            }
        for (int l = 0; l < a->num_levels; l++)
            level_min[l] = root_min[ctools_uf_find(uf, l)];
        int ng = 0;
        int remap_rc = remap_values_sorted(level_min, a->num_levels, mobility, &ng);
        free(root_min);
        free(level_min);
        if (remap_rc)
        {
            rc = 920;
            goto cleanup;
        }
    }
    int batch = pool < k ? pool : k;
    if (ctools_hdfe_alloc_buffers(&state, 0, batch))
    {
        rc = 920;
        goto cleanup;
    }
    data = malloc((size_t)n * batch * sizeof(*data));
    scales = malloc((size_t)batch * sizeof(*scales));
    store_cols = malloc((size_t)(batch > 2 ? batch : 2) * sizeof(*store_cols));
    store_vars = malloc((size_t)(batch > 2 ? batch : 2) * sizeof(*store_vars));
    if (!data || !scales || !store_cols || !store_vars)
    {
        rc = 920;
        goto cleanup;
    }
    double prepared = ctools_timer_seconds(), solve_time = 0, write_time = 0;
    int iterations = 0;
    for (int base = 0; base < k; base += batch)
    {
        int width = k - base < batch ? k - base : batch;
        /* The columns are independent: compact, rescale and check each on its
         * own thread; the first failure in column order is reported as before */
        int column_error = 0;
#pragma omp parallel for schedule(static) reduction(max : column_error) if (width > 1 && (size_t)n * width > 200000)
        for (int j = 0; j < width; j++) {
            double *column = data + (size_t)j * n;
            const double *loaded = fd.data.vars[base + j].data.dbl;
            if (n == nin)
                memcpy(column, loaded, (size_t)n * sizeof(*column));
            else
                for (int i = 0, q = 0; i < nin; i++)
                    if (keep[i])
                        column[q++] = loaded[i];
            scales[j] = chdfe_unit_scale(column, n);
            int status = 0;
            if (!isfinite(scales[j]))
                status = 2 * (width - j) + 1;
            for (int i = 0; !status && weights && i < n; i++) {
                double value = column[i];
                if (value != 0 && weights[i] * value * value == 0)
                    status = 2 * (width - j);
            }
            if (status > column_error) column_error = status;
        }
        if (column_error) {
            ctools_error("chdfe", (column_error & 1)
                ? "variable range is too wide to preserve the projection; rescale or restrict the input"
                : "weighted variable norm underflows; rescale or restrict the input");
            rc = 498;
            goto cleanup;
        }
        double tick = ctools_timer_seconds();
        /* transform() or acceleration() given: that iteration, on the rows */
        HDFE_SolveResult result =
            rowspace || transform != 0 || accel != 1
                ? partial_out_columns_options(&state, data, n, width, transform, accel)
                : partial_out_columns(&state, data, n, width, state.num_threads);
        solve_time += ctools_timer_seconds() - tick;
        rc = result.status;
        if (rc)
        {
            if (rc != 1) ctools_error("chdfe", "fixed effects did not converge");
            goto cleanup;
        }
        if (result.iterations > iterations)
            iterations = result.iterations;
        tick = ctools_timer_seconds();
        /* A retained singleton intercept spans its row exactly. Publishing
         * zero avoids reporting solver roundoff as a nonzero residual (zero
         * stays zero when the column's scale is applied in the store). */
        if (keepsingletons)
            for (int h = 0; h < g; h++)
                if (!slopes[h])
                {
                    FE_Factor *f = &state.factors[h];
                    for (int i = 0; i < n; i++)
                        if (f->counts[f->levels[i] - 1] == 1)
                            for (int j = 0; j < width; j++)
                                data[(size_t)j * n + i] = 0;
                }
        /* One pass over the records for the whole batch: the outputs are
         * adjacent in each record (one st_addvar), the power-of-two scale is
         * undone and every value checked as it is stored. With fillall the
         * rows outside the sample get missing, as a filled variable has. */
        {
            int unstorable = 0;
            for (int j = 0; j < width; j++)
            {
                store_cols[j] = data + (size_t)j * n;
                store_vars[j] = inputs + base + j + 1;
            }
            rc = ctools_stata_rc(ctools_store_filtered_cols(
                store_cols, scales, NULL, store_vars, (size_t)width, (size_t)n, obs,
                fillall ? (size_t)SF_nobs() : 0, &unstorable));
            if (rc)
                goto cleanup;
            if (unstorable)
            {
                rc = 498;
                goto cleanup;
            }
        }
        write_time += ctools_timer_seconds() - tick;
    }
    /* With "prefilled" the wrapper passes its sample marker (1 on every
     * loaded row) as the output: only the dropped rows are written, one store
     * per dropped row, and _chdfe_marked tells the wrapper that it was done. */
    const int prefilled = ctools_parse_bool_option(cursor, "prefilled");
    if (prefilled)
    {
        for (int i = 0; i < nin && n < nin; i++)
            if (!keep[i])
            {
                rc = (_stata_)->safestore(sampleout, fd.obs_map[i], 0);
                if (rc)
                    goto cleanup;
            }
        if ((rc = save_number("_chdfe_marked", 1)))
            goto cleanup;
    }
    if (fillall)
    {
        /* sample() (1 in the sample, 0 elsewhere) and groupvar() (missing
         * elsewhere) together, on every observation */
        double scale_one[2] = {1.0, 1.0}, other[2];
        size_t nout = 0;
        if (!prefilled)
        {
            store_cols[nout] = NULL;
            store_vars[nout] = sampleout;
            other[nout++] = 0.0;
        }
        if (mobility)
        {
            store_mobility = malloc((size_t)n * sizeof(*store_mobility));
            if (!store_mobility)
            {
                rc = 920;
                goto cleanup;
            }
            const int *codes = intercepts[first].levels;
#pragma omp parallel for schedule(static) if (n > 200000)
            for (int i = 0; i < n; i++)
                store_mobility[i] = mobility[codes[i] - 1];
            store_cols[nout] = store_mobility;
            store_vars[nout] = sampleout + 1;
            other[nout++] = SV_missval;
        }
        if (nout && (rc = ctools_stata_rc(ctools_store_filtered_cols(
                         store_cols, scale_one, other, store_vars, nout, (size_t)n, obs,
                         (size_t)SF_nobs(), NULL))))
            goto cleanup;
        if ((rc = save_number("_chdfe_filled", 1)))
            goto cleanup;
    }
    else
    {
        for (int i = 0; !prefilled && i < n; i++)
        {
            rc = (_stata_)->safestore(sampleout, obs[i], 1);
            if (rc)
                goto cleanup;
        }
        for (int i = 0; mobility && i < n; i++)
        {
            rc = (_stata_)->safestore(sampleout + 1, obs[i],
                                      mobility[intercepts[first].levels[i] - 1]);
            if (rc)
                goto cleanup;
        }
    }
    if (ctools_parse_bool_option(cursor, "storeids"))
    {
        double *codes = malloc((size_t)n * sizeof(*codes));
        if (!codes)
        {
            rc = 920;
            goto cleanup;
        }
        for (int j = 0; j < g + c; j++)
        {
            int source = j < g ? k + j : k + g + ns + j - g;
            for (int i = 0, q = 0; i < nin; i++)
                if (!keep || keep[i])
                    codes[q++] = fd.data.vars[source].data.dbl[i];
            int ng;
            if (remap_values_sorted(codes, n, clusters, &ng))
            {
                rc = 920;
                break;
            }
            for (int i = 0; i < n; i++)
            {
                rc = (_stata_)->safestore(source + 1, obs[i], clusters[i]);
                if (rc)
                    break;
            }
            if (rc)
                break;
        }
        free(codes);
        if (rc)
            goto cleanup;
    }
    int df = 0, m = 0, mnested = 0, nonnested = 0;
    for (int j = 0; j < g; j++)
    {
        df += state.factors[j].num_levels - redundant[j];
        m += redundant[j];
        if (nested[j])
            mnested += redundant[j];
        int exact = nested[j];
        if (!slopes[j] && !nested[j])
        {
            nonnested++;
            exact = nonnested == 1 || (nonnested == 2 && (flags & 3));
        }
        if ((rc = save_factor("K", j, state.factors[j].num_levels)) ||
            (rc = save_factor("M", j, redundant[j])) ||
            (rc = save_factor("nested", j, nested[j])) || (rc = save_factor("exact", j, exact)))
            goto cleanup;
    }
    if ((rc = save_number("_chdfe_N", n)) || (rc = save_number("_chdfe_df", df)) ||
        (rc = save_number("_chdfe_mobility", m)) || (rc = save_number("_chdfe_nested", mnested)) ||
        (rc = save_number("_chdfe_singletons", dropped)) ||
        (rc = save_number("_chdfe_iterations", iterations)) ||
        (rc = save_number("_chdfe_group_saved", mobility != NULL)))
        goto cleanup;
    if (ctools_parse_bool_option(cursor, "verbose"))
        ctools_msg("chdfe", "Load and prepare %.3fs; project %.3fs; store %.3fs; %d iterations",
                   prepared - start, solve_time, write_time, iterations);
cleanup:
    /* The loaded columns, factors, solver buffers and row arrays (over 1 GB
     * at 10M rows) are returned by the helper thread after the call
     * (ctools_release.h); small blocks are freed here */
    {
        HDFE_State *held = ctools_release_active() ? malloc(sizeof(*held)) : NULL;
        if (held)
        {
            *held = state;
            ctools_release_later_call(chdfe_release_state, held);
        }
        else
            ctools_hdfe_state_cleanup(&state);
    }
    ctools_filtered_data_release(&fd);
    if (raw)
        for (int j = 0; j < g; j++)
        {
            ctools_release_later(raw[j].levels);
            ctools_release_later(raw[j].counts);
        }
    free(raw);
    for (int j = 0; order && j < g; j++)
    {
        ctools_release_later(order[j].rows.first_row);
        ctools_release_later(order[j].rows.single_rows);
    }
    free(order);
    ctools_release_later(sing.keep);
    ctools_release_later(sing.dropped);
    for (int j = 0; sing.affected && j < sing.G; j++)
    {
        ctools_release_later(sing.affected[j].level);
        ctools_release_later(sing.affected[j].left);
    }
    free(sing.affected);
    free(slopes);
    free(levels);
    free(sizes);
    ctools_release_later(keep);
    free(columns);
    free(nested);
    free(redundant);
    free(values);
    ctools_release_later(weights);
    ctools_release_later(data);
    free(scales);
    ctools_release_later(obs);
    free(intercepts);
    free(imap);
    free(inested);
    free(ired);
    ctools_release_later(clusters);
    free(mobility);
    free(store_cols);
    free(store_vars);
    ctools_release_later(store_mobility);
    ctools_uf_destroy(uf);
    return rc;
}
