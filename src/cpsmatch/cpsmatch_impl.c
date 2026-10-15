#include "cpsmatch_internal.h"
#include "cpsmatch_impl.h"

/* Nearest-neighbor matching with ties(), with replacement, step by step as
 * psmatch2's match_pscore() does it (its weights are the reference):
 * - controls and treated are in psmatch2's order (pscore, ascending or
 *   descending, then observation), and treated are visited in that order;
 * - the nearest control is found by scanning on from the previous one in the
 *   current direction while the distance does not increase (the first of
 *   equal distances wins); once it reaches the last control, the direction
 *   turns backward and the scan no longer moves;
 * - ties() adds the controls that follow it in that direction at no greater
 *   distance, and neighbor(k) then adds the k-1 closest around it from both
 *   sides, so a control can enter a treated's list twice;
 * - every list entry gets 1/(list length), and a match needs a distance below
 *   caliper() (strict, as psmatch2). */
static int cpsmatch_match_psmatch2_ties(const sorted_control_t *controls, size_t n_controls,
                                        const sorted_control_t *treated, size_t n_treated,
                                        const size_t *treated_idx, const size_t *control_idx,
                                        const perm_idx_t *obs_map,
                                        int n_neighbors, double caliper,
                                        double *out_weight, double *out_match, double *out_nn,
                                        double *out_support, int *n_matched, int *n_off)
{
    if (n_controls == 0 || n_treated == 0) return 0;
    size_t *list = malloc((n_controls + (size_t)n_neighbors + 1) * sizeof(size_t));
    if (!list) return 920;
    const ssize_t j1 = (ssize_t)n_controls - 1;
    const double cal = caliper > 0.0 ? caliper : INFINITY;
    ssize_t jmatch = 0;
    int forward = 1;
    for (size_t t = 0; t < n_treated; t++) {
        double p = treated[t].pscore;
        size_t obs_idx = treated_idx[treated[t].orig_idx];
        ssize_t j = jmatch;
        double dif0, dif1 = fabs(p - controls[j].pscore);
        while (j >= 0 && j < j1) {
            j += forward ? 1 : -1;
            if (j < 0 || j > j1) break;
            dif0 = dif1;
            dif1 = fabs(p - controls[j].pscore);
            if (dif1 > dif0) j = j1;
            if (dif1 < dif0) jmatch = j;
        }
        if (fabs(p - controls[jmatch].pscore) < cal) {
            size_t n = 0;
            list[n++] = (size_t)jmatch;
            /* ties(): onward in the scan direction */
            double tie = fabs(p - controls[jmatch].pscore);
            for (ssize_t k = jmatch + (forward ? 1 : -1);
                 k >= 0 && k <= j1 && fabs(p - controls[k].pscore) <= tie;
                 k += forward ? 1 : -1) {
                list[n++] = (size_t)k;
            }
            /* neighbor(k): the k-1 closest around the nearest, left on equal
             * distance */
            ssize_t pos0 = jmatch - 1, pos1 = jmatch + 1, jm = jmatch;
            double d0 = 0.0, d1 = 0.0;
            for (int k = 1; k < n_neighbors; ) {
                if (pos0 >= 0) d0 = fabs(p - controls[pos0].pscore);
                if (pos1 <= j1) d1 = fabs(p - controls[pos1].pscore);
                if (pos0 >= 0 && (d0 <= d1 || pos1 > j1)) jm = pos0--;
                else if (pos1 <= j1) jm = pos1++;
                else break;
                if (!(fabs(p - controls[jm].pscore) < cal)) break;
                list[n++] = (size_t)jm;
                k++;
            }
            for (size_t m = 0; m < n; m++) {
                out_weight[control_idx[controls[list[m]].orig_idx]] += 1.0 / (double)n;
            }
            out_match[obs_idx] = (double)obs_map[control_idx[controls[jmatch].orig_idx]];
            out_nn[obs_idx] = (double)n;
            out_weight[obs_idx] = 1.0;
            (*n_matched)++;
        } else {
            out_support[obs_idx] = 0.0;
            (*n_off)++;
        }
        if (jmatch == j1 && t + 1 < n_treated) forward = 0;
    }
    free(list);
    return 0;
}

ST_retcode cpsmatch_main(const char *args)
{
    cpsmatch_options opts;
    double t_start, t_load, t_setup, t_sort, t_match, t_store, t_total;
    ST_retcode rc = 0;

    t_start = ctools_timer_seconds();
    t_setup = 0.0;
    t_sort = 0.0;

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cpsmatch: no arguments specified\n");
        return 198;
    }

    if (cpsmatch_parse_options(args, &opts) != 0) {
        SF_error("cpsmatch: failed to parse arguments\n");
        return 198;
    }

    opts.out_nn_idx = opts.outcome_idx > 0 ? 7 : 6;

    if (!opts.with_replace && (opts.method != MATCH_NEAREST || opts.n_neighbors != 1)) {
        SF_error("cpsmatch: noreplacement requires nearest-neighbor matching with neighbor(1)\n");
        return 198;
    }

    if (opts.nobs == 0) {
        SF_error("cpsmatch: no observations\n");
        return 2000;
    }

    /* Get observation range bounds (used for ctools_data_load) */
    (void)SF_in1();  /* obs1/obs2 handled by filtered load */
    (void)SF_in2();

    /* ========================================================================
     * PHASE 1: Load data with if/in filtering using ctools_data_load()
     * This loads ONLY observations that pass the if-condition, reducing memory
     * and eliminating the need for double iteration with SF_ifobs().
     * ======================================================================== */
    t_load = ctools_timer_seconds();

    /* Build array of variable indices to load */
    int nvars_to_load = (opts.outcome_idx > 0) ? 3 : 2;
    int *var_indices = malloc(nvars_to_load * sizeof(int));
    if (!var_indices) {
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }
    var_indices[0] = opts.treat_idx;
    var_indices[1] = opts.pscore_idx;
    if (opts.outcome_idx > 0) {
        var_indices[2] = opts.outcome_idx;
    }

    /* Load filtered data - only observations passing if-condition are loaded */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices,
                                                       nvars_to_load, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER);
    free(var_indices);

    if (load_rc != STATA_OK) {
        SF_error("cpsmatch: failed to load data\n");
        ctools_filtered_data_free(&filtered);
        return ctools_stata_rc(load_rc);
    }

    /* Get pointers to the filtered data arrays */
    size_t n_filtered = filtered.data.nobs;
    double *treatment = filtered.data.vars[0].data.dbl;
    double *pscore = filtered.data.vars[1].data.dbl;
    perm_idx_t *obs_map = filtered.obs_map;  /* Maps filtered index -> Stata obs */
    /* outcome is loaded for potential ATT/ATE computation but not currently used */
    double *outcome __attribute__((unused)) = (opts.outcome_idx > 0) ? filtered.data.vars[2].data.dbl : NULL;

    t_load = ctools_timer_seconds() - t_load;

    /* ========================================================================
     * PHASE 2: Separate treated and control groups (single pass, no SF_ifobs)
     * Since data is already filtered, we just need to check treatment values.
     * ======================================================================== */
    t_setup = ctools_timer_seconds();

    /* Count treated and controls - no SF_ifobs needed, already filtered */
    size_t n_treated = 0, n_controls = 0;
    for (size_t i = 0; i < n_filtered; i++) {
        /* psmatch2 excludes observations without a score from matching */
        if (SF_is_missing(treatment[i]) || SF_is_missing(pscore[i])) continue;
        if (treatment[i] != 0.0) n_treated++;
        else n_controls++;
    }

    if (n_treated == 0 || n_controls == 0) {
        ctools_filtered_data_free(&filtered);
        SF_error("cpsmatch: need both treated and control observations\n");
        return 198;
    }

    /* Build index arrays - indices now refer to position in filtered array */
    size_t *treated_idx = malloc(n_treated * sizeof(size_t));
    size_t *control_idx = malloc(n_controls * sizeof(size_t));
    double *treated_pscore = malloc(n_treated * sizeof(double));
    double *control_pscore = malloc(n_controls * sizeof(double));

    if (!treated_idx || !control_idx || !treated_pscore || !control_pscore) {
        ctools_filtered_data_free(&filtered);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    size_t ti = 0, ci = 0;
    for (size_t i = 0; i < n_filtered; i++) {
        if (SF_is_missing(treatment[i]) || SF_is_missing(pscore[i])) continue;
        if (treatment[i] != 0.0) {
            treated_idx[ti] = i;  /* Index in filtered array */
            treated_pscore[ti] = pscore[i];
            ti++;
        } else {
            control_idx[ci] = i;  /* Index in filtered array */
            control_pscore[ci] = pscore[i];
            ci++;
        }
    }

    /* ========================================================================
     * PHASE 3: Compute common support region
     * ======================================================================== */

    double ps_min_treated = DBL_MAX, ps_max_treated = -DBL_MAX;
    double ps_min_control = DBL_MAX, ps_max_control = -DBL_MAX;

    for (size_t i = 0; i < n_treated; i++) {
        if (!SF_is_missing(treated_pscore[i])) {
            if (treated_pscore[i] < ps_min_treated) ps_min_treated = treated_pscore[i];
            if (treated_pscore[i] > ps_max_treated) ps_max_treated = treated_pscore[i];
        }
    }

    for (size_t i = 0; i < n_controls; i++) {
        if (!SF_is_missing(control_pscore[i])) {
            if (control_pscore[i] < ps_min_control) ps_min_control = control_pscore[i];
            if (control_pscore[i] > ps_max_control) ps_max_control = control_pscore[i];
        }
    }

    /* Common support: psmatch2 restricts TREATED to be within control pscore range.
     * Controls are not restricted - they're always on support for ATT estimation. */
    double common_min = ps_min_control;
    double common_max = ps_max_control;

    /* ========================================================================
     * PHASE 4: Perform matching
     * ======================================================================== */

    /* Allocate output arrays - sized by n_filtered (not nobs_range) for memory efficiency.
     * These are indexed by position in the filtered array, and use obs_map for write-back. */
    double *out_weight = calloc(n_filtered, sizeof(double));
    double *out_match = malloc(n_filtered * sizeof(double));
    double *out_support = malloc(n_filtered * sizeof(double));
    double *out_nn = malloc(n_filtered * sizeof(double));

    if (!out_weight || !out_match || !out_support || !out_nn) {
        ctools_filtered_data_free(&filtered);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        free(out_weight);
        free(out_match);
        free(out_support); free(out_nn);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    /* Initialize to missing */
    for (size_t i = 0; i < n_filtered; i++) {
        out_match[i] = SV_missval;
        out_support[i] = SV_missval;
        out_nn[i] = SV_missval;
    }

    for (size_t t = 0; t < n_treated; t++) out_nn[treated_idx[t]] = 0.0;

    /* Caliper is always in absolute units (propensity score scale) for psmatch2 compatibility */
    double caliper = opts.caliper;

    int n_matched_treated = 0;
    int n_off_support = 0;

    /* End setup phase, start sort/index phase */
    t_setup = ctools_timer_seconds() - t_setup;
    t_sort = ctools_timer_seconds();

    /* For matching without replacement, we need to:
     * 1. Sort treated by propensity score (descending by default, like psmatch2)
     * 2. Process sequentially to avoid race conditions
     * 3. Mark excess treated as off-support when controls run out
     */

    /* Create sorted order for treated (by propensity score) */
    size_t *treated_order = malloc(n_treated * sizeof(size_t));
    if (!treated_order) {
        ctools_filtered_data_free(&filtered);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        free(out_weight);
        free(out_match);
        free(out_support); free(out_nn);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    /* For noreplacement, we need to sort treated by propensity score.
     * Use qsort with sorted_control_t structure for O(n log n) instead of O(n²) insertion sort */
    sorted_control_t *sorted_treated = NULL;
    if (!opts.with_replace) {
        sorted_treated = malloc(n_treated * sizeof(sorted_control_t));
        if (!sorted_treated) {
            ctools_filtered_data_free(&filtered);
            free(treated_idx);
            free(control_idx);
            free(treated_pscore);
            free(control_pscore);
            free(out_weight);
            free(out_match);
            free(out_support); free(out_nn);
            free(treated_order);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }

        /* Populate and sort */
        for (size_t i = 0; i < n_treated; i++) {
            sorted_treated[i].orig_idx = i;
            sorted_treated[i].obs_idx = treated_idx[i];
            sorted_treated[i].pscore = treated_pscore[i];
        }

        /* Sort: ascending (default) or descending */
        if (opts.descending) {
            qsort(sorted_treated, n_treated, sizeof(sorted_control_t), cpsmatch_compare_by_pscore_desc);
        } else {
            qsort(sorted_treated, n_treated, sizeof(sorted_control_t), cpsmatch_compare_by_pscore);
        }

        /* Copy sorted indices to treated_order */
        for (size_t i = 0; i < n_treated; i++) {
            treated_order[i] = sorted_treated[i].orig_idx;
        }
        free(sorted_treated);
    } else {
        /* For with-replacement, order doesn't matter - use original order */
        for (size_t i = 0; i < n_treated; i++) {
            treated_order[i] = i;
        }
    }

    /* Count available controls for matching without replacement */
    int n_controls_available = (int)n_controls;

    /* Matching loop - sequential for noreplacement, parallel otherwise */
    if (!opts.with_replace) {
        /* ================================================================
         * FAST sequential matching for noreplacement using sorted index
         * ================================================================ */

        /* Build sorted control index for O(log n) lookups */
        sorted_control_t *sorted_controls = malloc(n_controls * sizeof(sorted_control_t));
        int *sorted_available = malloc(n_controls * sizeof(int));  /* Track available in sorted order */

        if (!sorted_controls || !sorted_available) {
            free(sorted_controls);
            free(sorted_available);
            ctools_filtered_data_free(&filtered);
            free(treated_idx);
            free(control_idx);
            free(treated_pscore);
            free(control_pscore);
            free(out_weight);
            free(out_match);
            free(out_support); free(out_nn);
            free(treated_order);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }

        /* Populate and sort control array */
        for (size_t i = 0; i < n_controls; i++) {
            sorted_controls[i].orig_idx = i;
            sorted_controls[i].obs_idx = obs_map[control_idx[i]];  /* Stata observation number for stable sort */
            sorted_controls[i].pscore = control_pscore[i];
            sorted_available[i] = 1;
        }
        qsort(sorted_controls, n_controls, sizeof(sorted_control_t), cpsmatch_compare_by_pscore);

        /* End sort phase, start matching phase */
        t_sort = ctools_timer_seconds() - t_sort;
        t_match = ctools_timer_seconds();

        /* Sequential matching using binary search */
        for (size_t t_idx = 0; t_idx < n_treated; t_idx++) {
            size_t t = treated_order[t_idx];
            size_t obs_idx = treated_idx[t];
            double ps = treated_pscore[t];

            /* Check common support */
            if (opts.common_support) {
                if (ps < common_min || ps > common_max) {
                    out_support[obs_idx] = 0.0;
                    n_off_support++;
                    continue;
                }
            }

            /* Check if any controls are available */
            if (n_controls_available <= 0) {
                /* No more controls - treated is off support */
                out_support[obs_idx] = 0.0;
                n_off_support++;
                continue;
            }

            out_support[obs_idx] = 1.0;

            /* Binary search to find insertion point */
            size_t insert_pos = cpsmatch_binary_search_left(sorted_controls, n_controls, ps);

            /* Compare the nearest available control on each side of insertion.
             * Equal distances choose the left candidate deterministically. */
            ssize_t left = (ssize_t)insert_pos - 1;
            size_t right = insert_pos;
            size_t best_idx = SIZE_MAX;
            double best_dist = DBL_MAX;
            while (left >= 0 && !sorted_available[left]) left--;
            while (right < n_controls && !sorted_available[right]) right++;
            if (left >= 0) {
                double dist = fabs(ps - sorted_controls[left].pscore);
                if (caliper <= 0.0 || dist <= caliper) {
                    best_idx = (size_t)left;
                    best_dist = dist;
                }
            }
            if (right < n_controls) {
                double dist = fabs(ps - sorted_controls[right].pscore);
                if (dist < best_dist && (caliper <= 0.0 || dist <= caliper)) {
                    best_idx = right;
                    best_dist = dist;
                }
            }

            if (best_idx != SIZE_MAX) {
                n_matched_treated++;
                out_nn[obs_idx] = 1.0;
                n_controls_available--;

                /* Mark control as used */
                sorted_available[best_idx] = 0;

                /* Get original control index in filtered array */
                size_t orig_ctrl_idx = sorted_controls[best_idx].orig_idx;
                size_t match_filtered_idx = control_idx[orig_ctrl_idx];

                /* Store match ID - use obs_map to get Stata observation number */
                out_match[obs_idx] = (double)obs_map[match_filtered_idx];

                /* Set weights - index by filtered position */
                out_weight[match_filtered_idx] = 1.0;
                out_weight[obs_idx] = 1.0;
            } else {
                /* No match found (e.g., caliper constraint) - treated is off support */
                out_support[obs_idx] = 0.0;
                n_off_support++;
            }
        }

        free(sorted_controls);
        free(sorted_available);
    } else {
        /* ================================================================
         * FAST matching for with-replacement using sorted control index
         * ================================================================ */

        /* Build sorted control index for O(log n) lookups */
        sorted_control_t *sorted_controls = malloc(n_controls * sizeof(sorted_control_t));
        if (!sorted_controls) {
            ctools_filtered_data_free(&filtered);
            free(treated_idx);
            free(control_idx);
            free(treated_pscore);
            free(control_pscore);
            free(out_weight);
            free(out_match);
            free(out_support); free(out_nn);
            free(treated_order);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }

        /* Populate sorted control array */
        for (size_t i = 0; i < n_controls; i++) {
            sorted_controls[i].orig_idx = i;
            sorted_controls[i].obs_idx = obs_map[control_idx[i]];  /* Stata observation number for stable sort */
            sorted_controls[i].pscore = control_pscore[i];
        }

        /* Sort controls by propensity score */
        qsort(sorted_controls, n_controls, sizeof(sorted_control_t), cpsmatch_compare_by_pscore);

        /* End sort phase, start matching phase */
        t_sort = ctools_timer_seconds() - t_sort;
        t_match = ctools_timer_seconds();

        if (opts.method == MATCH_NEAREST && opts.ties) {
            /* psmatch2's sequential ties() algorithm, in its sort order */
            sorted_control_t *ctrl = malloc((n_controls ? n_controls : 1) * sizeof(sorted_control_t));
            sorted_control_t *trt = malloc((n_treated ? n_treated : 1) * sizeof(sorted_control_t));
            if (!ctrl || !trt) {
                free(ctrl); free(trt);
                rc = 920;
            } else {
                int (*cmp)(const void *, const void *) = opts.descending ?
                    cpsmatch_compare_by_pscore_desc : cpsmatch_compare_by_pscore;
                memcpy(ctrl, sorted_controls, n_controls * sizeof(sorted_control_t));
                if (opts.descending) qsort(ctrl, n_controls, sizeof(sorted_control_t), cmp);
                size_t n_on = 0;
                for (size_t t = 0; t < n_treated; t++) {
                    size_t obs_idx = treated_idx[t];
                    double ps = treated_pscore[t];
                    if (opts.common_support && (ps < common_min || ps > common_max)) {
                        out_support[obs_idx] = 0.0;
                        n_off_support++;
                        continue;
                    }
                    out_support[obs_idx] = 1.0;
                    trt[n_on].orig_idx = t;
                    trt[n_on].obs_idx = obs_map[obs_idx];
                    trt[n_on].pscore = ps;
                    n_on++;
                }
                qsort(trt, n_on, sizeof(sorted_control_t), cmp);
                rc = cpsmatch_match_psmatch2_ties(ctrl, n_controls, trt, n_on, treated_idx,
                                                  control_idx, obs_map, opts.n_neighbors, caliper,
                                                  out_weight, out_match, out_nn, out_support,
                                                  &n_matched_treated, &n_off_support);
                free(ctrl); free(trt);
            }
            if (rc) {
                free(sorted_controls);
                ctools_filtered_data_free(&filtered);
                free(treated_idx);
                free(control_idx);
                free(treated_pscore);
                free(control_pscore);
                free(out_weight);
                free(out_match);
                free(out_support); free(out_nn);
                free(treated_order);
                SF_error("cpsmatch: memory allocation failed\n");
                return 920;
            }
        } else {
        /* Optimized parallel matching using fast O(log n) algorithms */

        /* Get number of threads */
        int n_threads = 1;
        #ifdef _OPENMP
        n_threads = ctools_get_max_threads();
        #endif

        /* Allocate thread-local weight buffers (shared by all methods) */
        double **thread_weights = malloc(n_threads * sizeof(double *));
        if (thread_weights) {
            for (int i = 0; i < n_threads; i++) {
                thread_weights[i] = calloc(n_filtered, sizeof(double));
            }
        }

        /* Allocate per-thread workspaces (for NN/kernel only) */
        match_workspace *thread_ws = NULL;
        if (opts.method != MATCH_RADIUS) {
            thread_ws = malloc(n_threads * sizeof(match_workspace));
            if (!thread_ws) {
                if (thread_weights) {
                    for (int i = 0; i < n_threads; i++) free(thread_weights[i]);
                    free(thread_weights);
                }
                free(sorted_controls);
                ctools_filtered_data_free(&filtered);
                free(treated_idx);
                free(control_idx);
                free(treated_pscore);
                free(control_pscore);
                free(out_weight);
                free(out_match);
                free(out_support); free(out_nn);
                free(treated_order);
                SF_error("cpsmatch: memory allocation failed\n");
                return 920;
            }
            for (int i = 0; i < n_threads; i++) {
                cpsmatch_match_workspace_init(&thread_ws[i], MATCH_WORKSPACE_INITIAL_CAP);
            }
        }

        if (opts.method == MATCH_RADIUS) {
            double radius = (caliper > 0.0) ? caliper : CPSMATCH_DEFAULT_CALIPER;

            #pragma omp parallel reduction(+:n_matched_treated, n_off_support)
            {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif

                double *local_weights = (thread_weights && thread_weights[tid]) ? thread_weights[tid] : out_weight;
                int use_local = (thread_weights && thread_weights[tid]) ? 1 : 0;

                #pragma omp for schedule(dynamic, 64)
                for (size_t t = 0; t < n_treated; t++) {
                    size_t obs_idx = treated_idx[t];
                    double ps = treated_pscore[t];

                    /* Check common support */
                    if (opts.common_support) {
                        if (ps < common_min || ps > common_max) {
                            out_support[obs_idx] = 0.0;
                            n_off_support++;
                            continue;
                        }
                    }
                    out_support[obs_idx] = 1.0;

                    /* Fast range query - no allocation */
                    radius_match_info info;
                    cpsmatch_find_radius_range_fast(ps, sorted_controls, n_controls, radius, &info);

                    if (info.n_matches > 0) {
                        n_matched_treated++;
                        out_nn[obs_idx] = (double)info.n_matches;

                        /* Store first match ID - use obs_map for Stata observation number */
                        size_t first_match_filtered = control_idx[sorted_controls[info.start].orig_idx];
                        out_match[obs_idx] = (double)obs_map[first_match_filtered];

                        /* Accumulate weights - use local buffer if available */
                        double weight = 1.0 / info.n_matches;
                        for (size_t i = info.start; i < info.end; i++) {
                            size_t match_idx = control_idx[sorted_controls[i].orig_idx];
                            if (use_local) {
                                local_weights[match_idx] += weight;
                            } else {
                                #pragma omp atomic
                                out_weight[match_idx] += weight;
                            }
                        }

                        /* Weight for treated observation = 1 */
                        if (use_local) {
                            local_weights[obs_idx] = 1.0;
                        } else {
                            out_weight[obs_idx] = 1.0;
                        }
                    } else {
                        /* No match found - treated is off support */
                        out_support[obs_idx] = 0.0;
                        n_off_support++;
                    }
                }
            }
        } else {
            /* NN and kernel matching with thread-local weights and workspaces */
            #pragma omp parallel reduction(+:n_matched_treated, n_off_support)
            {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif

                double *local_weights = (thread_weights && thread_weights[tid]) ? thread_weights[tid] : out_weight;
                int use_local = (thread_weights && thread_weights[tid]) ? 1 : 0;
                match_workspace *ws = (thread_ws) ? &thread_ws[tid] : NULL;

                #pragma omp for schedule(dynamic, 64)
                for (size_t t = 0; t < n_treated; t++) {
                    size_t obs_idx = treated_idx[t];
                    double ps = treated_pscore[t];

                    /* Check common support */
                    if (opts.common_support) {
                        if (ps < common_min || ps > common_max) {
                            out_support[obs_idx] = 0.0;
                            n_off_support++;
                            continue;
                        }
                    }
                    out_support[obs_idx] = 1.0;

                    match_result result;
                    memset(&result, 0, sizeof(result));

                    int match_rc = 0;

                    switch (opts.method) {
                        case MATCH_NEAREST:
                            match_rc = cpsmatch_find_nearest_neighbors_fast(
                                ps, sorted_controls, control_idx, n_controls,
                                opts.n_neighbors, caliper, opts.ties, ws, &result);
                            break;

                        case MATCH_KERNEL:
                            match_rc = cpsmatch_find_kernel_matches_fast(
                                ps, sorted_controls, control_idx, n_controls,
                                opts.bandwidth, opts.kernel, caliper, ws, &result);
                            break;

                        default:
                            match_rc = cpsmatch_find_nearest_neighbors_fast(
                                ps, sorted_controls, control_idx, n_controls,
                                opts.n_neighbors, caliper, opts.ties, ws, &result);
                            break;
                    }

                    if (match_rc == 0 && result.n_matches > 0) {
                        n_matched_treated++;
                        out_nn[obs_idx] = (double)result.n_matches;

                        /* Store first match ID - use obs_map for Stata observation number */
                        out_match[obs_idx] = (double)obs_map[result.match_ids[0]];

                        /* Accumulate weights using thread-local buffer (no atomics) */
                        for (int m = 0; m < result.n_matches; m++) {
                            size_t match_filtered_idx = result.match_ids[m];
                            if (use_local) {
                                local_weights[match_filtered_idx] += result.match_weights[m];
                            } else {
                                #pragma omp atomic
                                out_weight[match_filtered_idx] += result.match_weights[m];
                            }
                        }

                        /* Weight for treated observation = 1 */
                        if (use_local) {
                            local_weights[obs_idx] = 1.0;
                        } else {
                            out_weight[obs_idx] = 1.0;
                        }
                    } else {
                        /* No match found (e.g., caliper/bandwidth constraint) */
                        out_support[obs_idx] = 0.0;
                        n_off_support++;
                    }
                    /* No free needed - workspace owns the buffers */
                }
            }
        }

        /* Shared merge phase: reduce thread_weights into out_weight */
        if (thread_weights) {
            #pragma omp parallel for schedule(static)
            for (size_t j = 0; j < n_filtered; j++) {
                double sum = 0.0;
                for (int i = 0; i < n_threads; i++) {
                    if (thread_weights[i]) {
                        sum += thread_weights[i][j];
                    }
                }
                out_weight[j] += sum;
            }
            for (int i = 0; i < n_threads; i++) {
                free(thread_weights[i]);
            }
            free(thread_weights);
        }

        /* Free per-thread workspaces */
        if (thread_ws) {
            for (int i = 0; i < n_threads; i++) {
                cpsmatch_match_workspace_free(&thread_ws[i]);
            }
            free(thread_ws);
        }
        } /* end parallel matching */

        free(sorted_controls);
    }

    free(treated_order);

    /* Mark controls on support - for psmatch2 compatibility, controls are always
     * on support for ATT estimation. The common support restriction only applies
     * to treated observations. */
    for (size_t i = 0; i < n_controls; i++) {
        size_t obs_idx = control_idx[i];
        out_support[obs_idx] = 1.0;
    }

    t_match = ctools_timer_seconds() - t_match;

    /* ========================================================================
     * PHASE 5: Store results back to Stata using obs_map
     * ======================================================================== */
    t_store = ctools_timer_seconds();

    /* Store results - use obs_map to write to correct Stata observations.
     * Output arrays are indexed by filtered position (0 to n_filtered-1).
     * Variable positions from the ado are already one-based. */
    for (size_t i = 0; i < n_filtered; i++) {
        ST_int obs = (ST_int)obs_map[i];

        rc = SF_vstore(opts.out_nn_idx, obs, out_nn[i]);
        if (rc) goto cleanup;

        if (opts.out_weight_idx > 0) {
            rc = SF_vstore(opts.out_weight_idx, obs, out_weight[i]);
            if (rc) goto cleanup;
        }

        if (opts.out_match_idx > 0) {
            rc = SF_vstore(opts.out_match_idx, obs, out_match[i]);
            if (rc) goto cleanup;
        }

        if (opts.out_support_idx > 0) {
            rc = SF_vstore(opts.out_support_idx, obs, out_support[i]);
            if (rc) goto cleanup;
        }
    }

    t_store = ctools_timer_seconds() - t_store;

    /* ========================================================================
     * PHASE 6: Store timing and diagnostics
     * ======================================================================== */

    t_total = ctools_timer_seconds() - t_start;

    SF_scal_save("_cpsmatch_time_load", t_load);
    SF_scal_save("_cpsmatch_time_setup", t_setup);
    SF_scal_save("_cpsmatch_time_sort", t_sort);
    SF_scal_save("_cpsmatch_time_match", t_match);
    SF_scal_save("_cpsmatch_time_store", t_store);
    SF_scal_save("_cpsmatch_time_total", t_total);

    SF_scal_save("_cpsmatch_n_treated", (double)n_treated);
    SF_scal_save("_cpsmatch_n_controls", (double)n_controls);
    SF_scal_save("_cpsmatch_n_matched", (double)n_matched_treated);
    SF_scal_save("_cpsmatch_n_off_support", (double)n_off_support);
    SF_scal_save("_cpsmatch_common_min", common_min);
    SF_scal_save("_cpsmatch_common_max", common_max);

    /* Thread diagnostics */
    CTOOLS_SAVE_THREAD_INFO("_cpsmatch");

cleanup:
    /* Cleanup - use stata_data_free for input data (frees treatment, pscore, outcome) */
    ctools_filtered_data_free(&filtered);
    free(treated_idx);
    free(control_idx);
    free(treated_pscore);
    free(control_pscore);
    free(out_weight);
    free(out_match);
    free(out_support); free(out_nn);

    return rc;
}
