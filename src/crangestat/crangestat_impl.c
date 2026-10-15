#include "crangestat_internal.h"
#include "crangestat_impl.h"

ST_retcode crangestat_main(const char *args)
{
    double t_start, t_load, t_sort, t_groups, t_compute, t_store;
    int *source_indices = NULL;
    int *result_indices = NULL;
    int *by_indices = NULL;
    int *load_indices = NULL;
    double **source_data = NULL;
    double **result_data = NULL;
    double *key_data = NULL;
    double **work_arrays = NULL;
    stat_spec *specs = NULL;
    group_info *groups = NULL;
    size_t *inv_perm = NULL;
    perm_idx_t *fwd_perm = NULL;
    ctools_filtered_data filtered;
    perm_idx_t *obs_map = NULL;
    int num_threads = 1;

    t_start = ctools_timer_seconds();

    if (!args || strlen(args) == 0) {
        ctools_error(CRANGESTAT_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    crangestat_config config;
    config.exclude_self = ctools_parse_bool_option(args, "excludeself");
    config.verbose = ctools_parse_bool_option(args, "verbose");
    config.interval_low = ctools_parse_double_option(args, "low", 0.0);
    config.interval_high = ctools_parse_double_option(args, "high", 0.0);
    if (!isfinite(config.interval_low) || !isfinite(config.interval_high)) return 198;
    config.low_is_missing = ctools_parse_bool_option(args, "low_miss");
    config.high_is_missing = ctools_parse_bool_option(args, "high_miss");

    /*
       Argument format:
       "nstats nsource nby key_idx source_indices... result_indices... by_indices... stat_types... [options]"

       stat_types are space-separated integers corresponding to stat_type enum
    */

    const char *p = args;
    while (*p == ' ' || *p == '\t') p++;

    int nstats_value;
    if (ctools_parse_next_int(&p, &nstats_value) || nstats_value < 0 ||
        nstats_value > SF_nvars() || nstats_value == 0) {
        ctools_error(CRANGESTAT_MODULE, "invalid nstats");
        return 198;
    }
    size_t nstats = (size_t)nstats_value;

    while (*p == ' ' || *p == '\t') p++;
    int nsource_value;
    if (ctools_parse_next_int(&p, &nsource_value) || nsource_value < 0 ||
        nsource_value > SF_nvars() || nsource_value == 0) {
        ctools_error(CRANGESTAT_MODULE, "invalid nsource");
        return 198;
    }
    size_t nsource = (size_t)nsource_value;

    while (*p == ' ' || *p == '\t') p++;
    int nby_value;
    if (ctools_parse_next_int(&p, &nby_value) || nby_value < 0 ||
        nby_value > SF_nvars()) {
        ctools_error(CRANGESTAT_MODULE, "invalid nby");
        return 198;
    }
    size_t nby = (size_t)nby_value;

    while (*p == ' ' || *p == '\t') p++;
    int key_idx;
    if (ctools_parse_next_int(&p, &key_idx) || key_idx < 1 || key_idx > SF_nvars()) {
        ctools_error(CRANGESTAT_MODULE, "invalid key_idx");
        return 198;
    }

    /* Allocate arrays */
    source_indices = (int *)malloc(nsource * sizeof(int));
    result_indices = (int *)malloc(nstats * sizeof(int));
    specs = (stat_spec *)malloc(nstats * sizeof(stat_spec));

    if (!source_indices || !result_indices || !specs) {
        if (source_indices) free(source_indices);
        if (result_indices) free(result_indices);
        if (specs) free(specs);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    /* Parse source variable indices */
    for (size_t i = 0; i < nsource; i++) {
        while (*p == ' ' || *p == '\t') p++;
        if (ctools_parse_next_int(&p, &source_indices[i]) || source_indices[i] < 1 || source_indices[i] > SF_nvars()) {
            free(source_indices);
            free(result_indices);
            free(specs);
            ctools_error(CRANGESTAT_MODULE, "failed to parse source indices");
            return 198;
        }
    }

    /* Parse result variable indices */
    for (size_t i = 0; i < nstats; i++) {
        while (*p == ' ' || *p == '\t') p++;
        if (ctools_parse_next_int(&p, &result_indices[i]) || result_indices[i] < 1 || result_indices[i] > SF_nvars()) {
            free(source_indices);
            free(result_indices);
            free(specs);
            ctools_error(CRANGESTAT_MODULE, "failed to parse result indices");
            return 198;
        }
    }

    /* Parse by-variable indices if present */
    if (nby > 0) {
        by_indices = (int *)malloc(nby * sizeof(int));
        if (!by_indices) {
            free(source_indices);
            free(result_indices);
            free(specs);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }

        for (size_t i = 0; i < nby; i++) {
            while (*p == ' ' || *p == '\t') p++;
            if (ctools_parse_next_int(&p, &by_indices[i]) || by_indices[i] < 1 || by_indices[i] > SF_nvars()) {
                free(source_indices);
                free(result_indices);
                free(specs);
                free(by_indices);
                ctools_error(CRANGESTAT_MODULE, "failed to parse by indices");
                return 198;
            }
        }
    }

    /* Parse stat types (as a group) */
    for (size_t i = 0; i < nstats; i++) {
        while (*p == ' ' || *p == '\t') p++;
        int stat_int;
        if (ctools_parse_next_int(&p, &stat_int) || stat_int < 0 || stat_int >= STAT_COUNT_ALL) {
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error(CRANGESTAT_MODULE, "failed to parse stat types");
            return 198;
        }
        specs[i].type = (stat_type)stat_int;
        specs[i].result_var_idx = result_indices[i];
    }

    /* Parse source var mappings (as a group) */
    for (size_t i = 0; i < nstats; i++) {
        while (*p == ' ' || *p == '\t') p++;
        int src_idx;
        if (ctools_parse_next_int(&p, &src_idx) || src_idx < 0 || src_idx >= (int)nsource) {
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error(CRANGESTAT_MODULE, "failed to parse source var for stat");
            return 198;
        }
        specs[i].source_var_idx = src_idx;
    }

    /* Counts that disagree with the index lists shift every later field, so
     * a result variable must be distinct from the key, source and by
     * variables and from every other result variable */
    {
        unsigned char *role = (unsigned char *)calloc((size_t)SF_nvars() + 1, 1);
        if (!role) {
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }
        role[key_idx] = 1;
        for (size_t i = 0; i < nsource; i++) role[source_indices[i]] = 1;
        for (size_t i = 0; i < nby; i++) role[by_indices[i]] = 1;
        int overlap = 0;
        for (size_t i = 0; i < nstats && !overlap; i++) {
            if (role[result_indices[i]]) overlap = 1;
            role[result_indices[i]] = 1;
        }
        free(role);
        if (overlap) {
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error(CRANGESTAT_MODULE, "result variable overlaps an input or another result");
            return 198;
        }
    }

    #ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    #endif

    /*
     * === Load Phase ===
     * Use ctools_data_load for optimized parallel data loading with if/in filtering.
     *
     * Variable layout in loaded data:
     *   [0..nby-1]           = by-variables (for grouping)
     *   [nby]                = key variable (for window bounds)
     *   [nby+1..nby+nsource] = source variables (for stat computation)
     */
    double load_start = ctools_timer_seconds();

    size_t nvars_to_load = nby + 1 + nsource;
    load_indices = (int *)malloc(nvars_to_load * sizeof(int));
    if (!load_indices) {
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    /* Build load index array: [by_vars..., key, source_vars...] */
    size_t idx = 0;
    for (size_t b = 0; b < nby; b++) {
        load_indices[idx++] = by_indices[b];
    }
    load_indices[idx++] = key_idx;
    for (size_t v = 0; v < nsource; v++) {
        load_indices[idx++] = source_indices[v];
    }

    /* Load data using filtered loading (handles if/in at load time) */
    ctools_filtered_data_init(&filtered);

    stata_retcode rc = ctools_data_load(&filtered, load_indices, nvars_to_load, 0, 0, 0);
    if (rc != STATA_OK) {
            ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "failed to load data");
        return ctools_stata_rc(rc);
    }

    /* Get filtered observation count and obs_map */
    size_t nobs = filtered.data.nobs;
    obs_map = filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "no observations");
        return 2000;
    }

    t_load = ctools_timer_seconds() - load_start;

    /*
     * === Sort Phase ===
     * Sort by (by_vars, key) using optimized algorithm selection.
     *
     * PERFORMANCE: Use SORT_ALG_COUNTING which:
     *   1. Tries counting sort first (O(n+k), optimal for small-range keys like dates)
     *   2. Falls back to LSD radix if range is too large (still O(n*d), faster than comparison sorts)
     *
     * This provides massive speedups for typical rangestat workloads where the key
     * variable (e.g., date, time period) has a small number of unique values.
     * Benchmark: 20M obs with 20 unique dates: IPS4o=3.1s, Counting=0.01s (280x faster)
     */
    double sort_start = ctools_timer_seconds();

    /* Build sort variable indices: [1, 2, ..., nby, nby+1] (1-based into loaded data) */
    size_t nsort = nby + 1;
    int *sort_vars = (int *)malloc(nsort * sizeof(int));
    if (!sort_vars) {
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    for (size_t i = 0; i < nsort; i++) {
        sort_vars[i] = (int)(i + 1);  /* 1-based indices */
    }

    /* Sort using counting sort (with LSD radix fallback for large ranges) */
    rc = ctools_sort_dispatch(&filtered.data, sort_vars, nsort, SORT_ALG_COUNTING);
    if (rc != STATA_OK) {
        free(sort_vars);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "sort failed");
        return ctools_stata_rc(rc);
    }

    /*
     * Compute inverse permutation before applying.
     * sort_order[sorted_idx] = original_idx
     * inv_perm[original_idx] = sorted_idx
     *
     * This lets us map results (in sorted order) back to original observation order.
     */
    inv_perm = (size_t *)malloc(nobs * sizeof(size_t));
    if (!inv_perm) {
        free(sort_vars);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    /* Compute inverse: inv_perm[original_idx] = sorted_idx */
    for (size_t i = 0; i < nobs; i++) {
        inv_perm[filtered.data.sort_order[i]] = i;
    }

    /* Save forward permutation for sequential store-phase reads.
     * sort_order[sorted_idx] = original_idx, so iterating sorted_idx gives
     * sequential reads from data[] and scattered writes via obs_map[]. */
    fwd_perm = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
    if (fwd_perm) {
        memcpy(fwd_perm, filtered.data.sort_order, nobs * sizeof(perm_idx_t));
    }
    /* Non-fatal: falls back to inv_perm-based store if alloc fails */

    /* Apply permutation to sort data in place */
    rc = ctools_apply_permutation(&filtered.data);
    if (rc != STATA_OK) {
        free(fwd_perm);
        free(inv_perm);
        free(sort_vars);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "permutation failed");
        return ctools_stata_rc(rc);
    }

    free(sort_vars);

    t_sort = ctools_timer_seconds() - sort_start;

    /*
     * Set up pointers to sorted data for computation.
     * Data layout after loading:
     *   filtered.data.vars[0..nby-1] = by-variables
     *   filtered.data.vars[nby] = key variable
     *   filtered.data.vars[nby+1..nby+nsource] = source variables
     */
    key_data = filtered.data.vars[nby].data.dbl;

    /* Set up source_data pointers (don't allocate, just point into stata_data) */
    source_data = (double **)malloc(nsource * sizeof(double *));
    if (!source_data) {
        free(fwd_perm);
        free(inv_perm);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }
    for (size_t v = 0; v < nsource; v++) {
        source_data[v] = filtered.data.vars[nby + 1 + v].data.dbl;
    }

    /* Allocate result arrays (these are separate, not in stata_data) */
    result_data = (double **)malloc(nstats * sizeof(double *));
    if (!result_data) {
        free(source_data);
        free(fwd_perm);
        free(inv_perm);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    for (size_t s = 0; s < nstats; s++) {
        result_data[s] = (double *)malloc(nobs * sizeof(double));
        if (!result_data[s]) {
            for (size_t j = 0; j < s; j++) free(result_data[j]);
            free(result_data);
            free(source_data);
            free(fwd_perm);
            free(inv_perm);
            ctools_filtered_data_free(&filtered);
            free(load_indices);
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }
        /* Initialize to missing */
        for (size_t i = 0; i < nobs; i++) {
            result_data[s][i] = SV_missval;
        }
    }

    /* === Group Detection Phase === */
    double groups_start = ctools_timer_seconds();

    size_t ngroups = 1;
    groups = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    size_t *group_starts = (size_t *)malloc((nobs + 1) * sizeof(size_t));
    int *group_keys = (int *)malloc((nby > 0 ? nby : 1) * sizeof(int));
    if (!groups || !group_starts || !group_keys) {
        free(groups);
        free(group_starts);
        free(group_keys);
        free(fwd_perm);
        free(inv_perm);
        for (size_t s = 0; s < nstats; s++) free(result_data[s]);
        free(result_data);
        free(source_data);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    if (nby > 0) {
        /* The loaded columns are physically sorted by (by-vars, key), so
         * adjacent rows give Stata by-groups: strings compare byte-wise and
         * each missing code forms its own group. */
        for (size_t b = 0; b < nby; b++) group_keys[b] = (int)b;
        ngroups = ctools_group_starts(&filtered.data, NULL, group_keys, nby, group_starts);
        for (size_t g = 0; g < ngroups; g++) {
            groups[g].start = group_starts[g];
            groups[g].count = group_starts[g + 1] - group_starts[g];
        }
    } else {
        groups[0].start = 0;
        groups[0].count = nobs;
        ngroups = 1;
    }

    free(group_starts);
    free(group_keys);

    /* Shrink groups array from (nobs+1) to actual ngroups.
     * For 10M rows with 1000 groups: 160MB -> 16KB. */
    if (ngroups < nobs + 1) {
        group_info *shrunk = (group_info *)realloc(groups, ngroups * sizeof(group_info));
        if (shrunk) groups = shrunk;  /* realloc failure is non-fatal */
    }

    t_groups = ctools_timer_seconds() - groups_start;

    /* === Computation Phase === */
    double compute_start = ctools_timer_seconds();

    /* Compute max group size for work array allocation */
    size_t max_group_size = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (groups[g].count > max_group_size)
            max_group_size = groups[g].count;
    }

    size_t avg_group_size = nobs / ngroups;
    /* Groups below the within-group parallel threshold cannot amortize a
     * separate team per group. Keep them in one cross-group team. Retain the
     * existing small-average policy for highly unbalanced groups. The extra
     * cases fit the cross-group path's 16-source sparse-table buffer. */
    int use_cross_group_parallel = (ngroups >= CRANGESTAT_PARALLEL_GROUPS_MIN &&
                                    (avg_group_size < CRANGESTAT_SMALL_GROUP_THRESHOLD ||
                                     (max_group_size < CRANGESTAT_PARALLEL_THRESHOLD &&
                                      nsource <= 16)));
    int needs_work = 0;
    for (size_t s = 0; s < nstats; s++)
        needs_work |= crangestat_stat_needs_work_array(specs[s].type);
    num_threads = ctools_workspace_threads(num_threads,
        use_cross_group_parallel ? ngroups : max_group_size,
        needs_work ? max_group_size * sizeof(double) : 0);

    /* Allocate work arrays for each thread (sized to max_group_size, not nobs) */
    work_arrays = (double **)calloc((size_t)num_threads, sizeof(double *));
    if (!work_arrays) {
        free(groups);
        free(fwd_perm);
        free(inv_perm);
        for (size_t s = 0; s < nstats; s++) free(result_data[s]);
        free(result_data);
        free(source_data);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    for (int t = 0; needs_work && t < num_threads; t++) {
        work_arrays[t] = (double *)malloc(max_group_size * sizeof(double));
        if (!work_arrays[t]) {
            for (int j = 0; j < t; j++) free(work_arrays[j]);
            free(work_arrays);
            free(groups);
            free(fwd_perm);
            free(inv_perm);
            for (size_t s = 0; s < nstats; s++) free(result_data[s]);
            free(result_data);
            free(source_data);
            ctools_filtered_data_free(&filtered);
            free(load_indices);
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }
    }

    /*
     * OPTIMIZATION: Cross-group parallelism for many small groups
     *
     * When there are many by-groups with small average size, parallelizing
     * across groups is more efficient than parallelizing within each group.
     * This avoids OpenMP region setup/teardown overhead per group.
     */


    if (use_cross_group_parallel) {
        /*
         * Cross-group parallel processing: single parallel region over all groups.
         * Each thread processes complete groups independently.
         *
         * Pre-allocate per-thread O(1) buffers (prefix arrays reused across groups,
         * sparse tables built/freed per group since groups are small).
         */

        /* Determine which source variables need O(1) structures */
        int *a_needs_prefix = (int *)calloc(nsource, sizeof(int));
        int *a_needs_sparse = (int *)calloc(nsource, sizeof(int));
        int a_has_o1 = (a_needs_prefix && a_needs_sparse &&
                        max_group_size >= CRANGESTAT_PREFIX_SUM_THRESHOLD);

        if (a_has_o1) {
            for (size_t s = 0; s < nstats; s++) {
                int src_idx = specs[s].source_var_idx;
                if (src_idx < 0 || src_idx >= (int)nsource) continue;
                if (crangestat_stat_uses_prefix_sums(specs[s].type)) a_needs_prefix[src_idx] = 1;
                if (crangestat_stat_uses_sparse_table(specs[s].type)) a_needs_sparse[src_idx] = 1;
            }
        }

        /* Per-thread prefix arrays: thread t, source var v → a_prefix[t * nsource + v] */
        prefix_arrays *a_prefix = NULL;
        if (a_has_o1) {
            a_prefix = (prefix_arrays *)calloc((size_t)num_threads * nsource, sizeof(prefix_arrays));
            if (a_prefix) {
                for (int t = 0; t < num_threads; t++) {
                    for (size_t v = 0; v < nsource; v++) {
                        if (!a_needs_prefix[v]) continue;
                        size_t idx = (size_t)t * nsource + v;
                        a_prefix[idx].prefix_sum = (double *)malloc((max_group_size + 1) * sizeof(double));
                        a_prefix[idx].prefix_sum2 = (double *)malloc((max_group_size + 1) * sizeof(double));
                        a_prefix[idx].prefix_count = (size_t *)malloc((max_group_size + 1) * sizeof(size_t));
                        if (!a_prefix[idx].prefix_sum || !a_prefix[idx].prefix_sum2 ||
                            !a_prefix[idx].prefix_count) {
                            /* Disable prefix for this var on alloc failure */
                            a_needs_prefix[v] = 0;
                            if (a_prefix[idx].prefix_sum) { free(a_prefix[idx].prefix_sum); a_prefix[idx].prefix_sum = NULL; }
                            if (a_prefix[idx].prefix_sum2) { free(a_prefix[idx].prefix_sum2); a_prefix[idx].prefix_sum2 = NULL; }
                            if (a_prefix[idx].prefix_count) { free(a_prefix[idx].prefix_count); a_prefix[idx].prefix_count = NULL; }
                        }
                    }
                }
            } else {
                a_has_o1 = 0;
            }
        }

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        #endif
        for (size_t g = 0; g < ngroups; g++) {
            #ifdef _OPENMP
            int tid = omp_get_thread_num();
            #else
            int tid = 0;
            #endif

            double *work = work_arrays[tid];
            size_t g_start = groups[g].start;
            size_t g_count = groups[g].count;
            size_t g_end = g_start + g_count;
            /* Stata missing codes are finite; an unbounded window must still
             * exclude missing interval keys in every optimized path. */
            while (g_end > g_start && key_data[g_end - 1] >= SV_missval) g_end--;
            g_count = g_end - g_start;
            if (g_count == 0) continue;

            /* Build O(1) structures for this group */
            int group_has_o1 = a_has_o1 && (g_count >= CRANGESTAT_PREFIX_SUM_THRESHOLD);
            sparse_table group_sparse[16]; /* Stack-allocated for small nsource */
            sparse_table *group_sparse_ptr = (nsource <= 16) ? group_sparse : NULL;
            int sparse_built = 0;

            if (group_has_o1) {
                /* Build prefix arrays (reuse per-thread buffers) */
                for (size_t v = 0; v < nsource; v++) {
                    if (a_needs_prefix[v]) {
                        size_t idx = (size_t)tid * nsource + v;
                        a_prefix[idx].nobs = g_count;
                        crangestat_build_prefix_arrays(source_data[v], g_start, g_count, &a_prefix[idx]);
                    }
                }

                /* Build sparse tables per group (groups are small, so alloc overhead is minimal) */
                if (group_sparse_ptr) {
                    memset(group_sparse_ptr, 0, nsource * sizeof(sparse_table));
                    for (size_t v = 0; v < nsource; v++) {
                        if (a_needs_sparse[v]) {
                            if (crangestat_build_sparse_table(source_data[v], g_start, g_count,
                                                   &group_sparse_ptr[v], SV_missval) != 0) {
                                /* Release completed sources before using the fallback. */
                                for (size_t u = 0; u <= v; u++)
                                    crangestat_free_sparse_table(&group_sparse_ptr[u]);
                                group_sparse_ptr = NULL;
                                break;
                            }
                        }
                    }
                    if (group_sparse_ptr) sparse_built = 1;
                }
            }

            /* Process each observation in this group with sliding window.
             * Since data is sorted by key, window bounds only advance forward. */
            size_t win_lo = g_start;
            size_t win_hi = g_start;
            for (size_t i = g_start; i < g_end; i++) {
                double key_val = key_data[i];

                /* Missing keys sort to end — stop processing */
                if (key_val >= SV_missval) break;

                /* Without excludeself, equal keys have identical windows and
                 * results, including order statistics and first/last values.
                 * Reuse only within this group; missing keys stay untouched. */
                if (!config.exclude_self && i > g_start && key_val == key_data[i - 1]) {
                    for (size_t s = 0; s < nstats; s++)
                        result_data[s][i] = result_data[s][i - 1];
                    continue;
                }

                /* Compute window bounds */
                double low_bound = config.low_is_missing ? -DBL_MAX : key_val + config.interval_low;
                double high_bound = config.high_is_missing ? DBL_MAX : key_val + config.interval_high;

                /* Slide window pointers forward (O(1) amortized) */
                while (win_lo < g_end && key_data[win_lo] < low_bound) win_lo++;
                while (win_hi < g_end && key_data[win_hi] <= high_bound) win_hi++;

                size_t win_start = win_lo;
                size_t win_end = win_hi;

                /* Compute each statistic for this observation */
                for (size_t s = 0; s < nstats; s++) {
                    int src_idx = specs[s].source_var_idx;
                    if (src_idx < 0 || src_idx >= (int)nsource) continue;

                    double val;
                    stat_type stype = specs[s].type;

                    if (crangestat_stat_needs_work_array(stype)) {
                        val = crangestat_compute_window_stat_complex(
                            stype,
                            source_data[src_idx],
                            win_start,
                            win_end,
                            i,
                            config.exclude_self,
                            work
                        );
                    } else if (group_has_o1 && a_needs_prefix[src_idx] &&
                               crangestat_stat_uses_prefix_sums(stype)) {
                        size_t pidx = (size_t)tid * nsource + (size_t)src_idx;
                        val = crangestat_compute_stat_with_prefix(
                            stype,
                            &a_prefix[pidx],
                            source_data[src_idx],
                            g_start,
                            win_start,
                            win_end,
                            i,
                            config.exclude_self
                        );
                    } else if (group_has_o1 && sparse_built && a_needs_sparse[src_idx] &&
                               crangestat_stat_uses_sparse_table(stype)) {
                        val = crangestat_compute_minmax_with_sparse(
                            stype,
                            &group_sparse_ptr[src_idx],
                            source_data[src_idx],
                            g_start,
                            win_start,
                            win_end,
                            i,
                            config.exclude_self
                        );
                    } else {
                        val = crangestat_compute_window_stat_simple(
                            stype,
                            source_data[src_idx],
                            win_start,
                            win_end,
                            i,
                            config.exclude_self
                        );
                    }

                    result_data[s][i] = val;
                }
            }

            /* Free per-group sparse tables */
            if (sparse_built) {
                for (size_t v = 0; v < nsource; v++) {
                    if (a_needs_sparse[v]) {
                        crangestat_free_sparse_table(&group_sparse_ptr[v]);
                    }
                }
            }
        }

        /* Cleanup per-thread prefix arrays */
        if (a_prefix) {
            for (int t = 0; t < num_threads; t++) {
                for (size_t v = 0; v < nsource; v++) {
                    size_t idx = (size_t)t * nsource + v;
                    if (a_prefix[idx].prefix_sum) free(a_prefix[idx].prefix_sum);
                    if (a_prefix[idx].prefix_sum2) free(a_prefix[idx].prefix_sum2);
                    if (a_prefix[idx].prefix_count) free(a_prefix[idx].prefix_count);
                }
            }
            free(a_prefix);
        }
        if (a_needs_prefix) free(a_needs_prefix);
        if (a_needs_sparse) free(a_needs_sparse);
    } else {
        /* Standard path: process groups sequentially with internal parallelism */

        /*
         * OPTIMIZATION: O(1) range queries using prefix sums and sparse tables
         *
         * Pre-allocate O(1) data structures once at max_group_size, then
         * rebuild per group to avoid per-group alloc/free overhead.
         */
        int use_o1_queries = (max_group_size >= CRANGESTAT_PREFIX_SUM_THRESHOLD);

        /* Determine which source variables need which data structures */
        int *needs_prefix = NULL;
        int *needs_sparse = NULL;
        prefix_arrays *prefix_arr = NULL;
        sparse_table *sparse_tbl = NULL;

        if (use_o1_queries) {
            needs_prefix = (int *)calloc(nsource, sizeof(int));
            needs_sparse = (int *)calloc(nsource, sizeof(int));
            prefix_arr = (prefix_arrays *)calloc(nsource, sizeof(prefix_arrays));
            sparse_tbl = (sparse_table *)calloc(nsource, sizeof(sparse_table));

            if (!needs_prefix || !needs_sparse || !prefix_arr || !sparse_tbl) {
                use_o1_queries = 0;
                if (needs_prefix) free(needs_prefix);
                if (needs_sparse) free(needs_sparse);
                if (prefix_arr) free(prefix_arr);
                if (sparse_tbl) free(sparse_tbl);
                needs_prefix = NULL;
                needs_sparse = NULL;
                prefix_arr = NULL;
                sparse_tbl = NULL;
            }
        }

        if (use_o1_queries) {
            for (size_t s = 0; s < nstats; s++) {
                int src_idx = specs[s].source_var_idx;
                if (src_idx < 0 || src_idx >= (int)nsource) continue;

                if (crangestat_stat_uses_prefix_sums(specs[s].type)) {
                    needs_prefix[src_idx] = 1;
                }
                if (crangestat_stat_uses_sparse_table(specs[s].type)) {
                    needs_sparse[src_idx] = 1;
                }
            }

            /* Pre-allocate buffers at max_group_size */
            for (size_t v = 0; v < nsource; v++) {
                if (needs_prefix[v]) {
                    prefix_arr[v].prefix_sum = (double *)malloc((max_group_size + 1) * sizeof(double));
                    prefix_arr[v].prefix_sum2 = (double *)malloc((max_group_size + 1) * sizeof(double));
                    prefix_arr[v].prefix_count = (size_t *)malloc((max_group_size + 1) * sizeof(size_t));
                    prefix_arr[v].nobs = max_group_size;

                    if (!prefix_arr[v].prefix_sum || !prefix_arr[v].prefix_sum2 ||
                        !prefix_arr[v].prefix_count) {
                        needs_prefix[v] = 0;
                        if (prefix_arr[v].prefix_sum) free(prefix_arr[v].prefix_sum);
                        if (prefix_arr[v].prefix_sum2) free(prefix_arr[v].prefix_sum2);
                        if (prefix_arr[v].prefix_count) free(prefix_arr[v].prefix_count);
                        prefix_arr[v].prefix_sum = NULL;
                        prefix_arr[v].prefix_sum2 = NULL;
                        prefix_arr[v].prefix_count = NULL;
                    }
                }

                if (needs_sparse[v]) {
                    if (crangestat_alloc_sparse_table(&sparse_tbl[v], max_group_size) != 0) {
                        needs_sparse[v] = 0;
                    }
                }
            }
        }

        for (size_t g = 0; g < ngroups; g++) {
            size_t g_start = groups[g].start;
            size_t g_count = groups[g].count;
            size_t g_end = g_start + g_count;
            /* Stata missing codes are finite; an unbounded window must still
             * exclude missing interval keys in every optimized path. */
            while (g_end > g_start && key_data[g_end - 1] >= SV_missval) g_end--;
            g_count = g_end - g_start;
            if (g_count == 0) continue;

            /*
             * OPTIMIZATION: When excludeself is false, observations with the same
             * key value get the same window and same result. We process in "runs"
             * of same-key observations, computing once per unique key.
             *
             * This can provide massive speedups when there are many duplicate keys.
             * For example, panel data with 100K IDs and 20 time periods:
             * - Without optimization: 2M window computations
             * - With optimization: 20 window computations
             *
             * For small groups (< 1000 obs), skip run detection as overhead exceeds benefit.
             */
            if (!config.exclude_self && g_count >= CRANGESTAT_RUN_OPT_THRESHOLD) {
                /* First pass: identify run boundaries (starts of each unique key) */
                size_t *run_starts = (size_t *)malloc((g_count + 1) * sizeof(size_t));
                if (!run_starts) {
                    /* Fall back to standard path on allocation failure */
                    goto standard_obs_loop;
                }

                size_t num_runs = 0;
                run_starts[num_runs++] = g_start;

                for (size_t i = g_start + 1; i < g_end; i++) {
                    double prev_key = key_data[i - 1];
                    double curr_key = key_data[i];

                    /* Check for run boundary (different key or missing) */
                    int prev_miss = (prev_key >= SV_missval);
                    int curr_miss = (curr_key >= SV_missval);

                    if (prev_miss != curr_miss || (!prev_miss && prev_key != curr_key)) {
                        run_starts[num_runs++] = i;
                    }
                }
                run_starts[num_runs] = g_end;  /* Sentinel */

                /* Mostly distinct keys benefit from the linear sliding scan,
                 * rather than two binary searches and one scheduled task for
                 * almost every row. No arithmetic or window definition changes. */
                if (num_runs > g_count / 2) {
                    free(run_starts);
                    goto standard_obs_loop;
                }

                /* Build O(1) structures for this group (reuse pre-allocated buffers) */
                int group_has_o1 = use_o1_queries && (g_count >= CRANGESTAT_PREFIX_SUM_THRESHOLD);
                if (group_has_o1) {
                    for (size_t v = 0; v < nsource; v++) {
                        if (needs_prefix[v]) {
                            prefix_arr[v].nobs = g_count;
                            crangestat_build_prefix_arrays(source_data[v], g_start, g_count, &prefix_arr[v]);
                        }
                        if (needs_sparse[v]) {
                            crangestat_rebuild_sparse_table(&sparse_tbl[v], source_data[v], g_start,
                                                 g_count, SV_missval);
                        }
                    }
                }

                /* Process runs in parallel */
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
                #endif
                for (size_t r = 0; r < num_runs; r++) {
                    #ifdef _OPENMP
                    int tid = omp_get_thread_num();
                    #else
                    int tid = 0;
                    #endif

                    double *work = work_arrays[tid];
                    size_t run_start = run_starts[r];
                    size_t run_end = run_starts[r + 1];
                    double key_val = key_data[run_start];

                    /* Skip missing keys */
                    if (key_val >= SV_missval) {
                        continue;
                    }

                    /* Compute window bounds */
                    double low_bound = config.low_is_missing ? -DBL_MAX : key_val + config.interval_low;
                    double high_bound = config.high_is_missing ? DBL_MAX : key_val + config.interval_high;

                    /* Find window using binary search within group */
                    size_t win_start = crangestat_lower_bound(key_data, g_start, g_end, low_bound);
                    size_t win_end = crangestat_upper_bound(key_data, g_start, g_end, high_bound);

                    /* Compute each statistic ONCE for this run */
                    for (size_t s = 0; s < nstats; s++) {
                        int src_idx = specs[s].source_var_idx;
                        if (src_idx < 0 || src_idx >= (int)nsource) continue;

                        double val;
                        stat_type stype = specs[s].type;

                        if (crangestat_stat_needs_work_array(stype)) {
                            val = crangestat_compute_window_stat_complex(
                                stype,
                                source_data[src_idx],
                                win_start,
                                win_end,
                                0,  /* self_idx unused when exclude_self=0 */
                                0,  /* exclude_self=0 */
                                work
                            );
                        } else if (group_has_o1 && needs_prefix[src_idx] &&
                                   crangestat_stat_uses_prefix_sums(stype)) {
                            val = crangestat_compute_stat_with_prefix(
                                stype,
                                &prefix_arr[src_idx],
                                source_data[src_idx],
                                g_start,
                                win_start,
                                win_end,
                                0,
                                0
                            );
                        } else if (group_has_o1 && needs_sparse[src_idx] &&
                                   crangestat_stat_uses_sparse_table(stype)) {
                            val = crangestat_compute_minmax_with_sparse(
                                stype,
                                &sparse_tbl[src_idx],
                                source_data[src_idx],
                                g_start,
                                win_start,
                                win_end,
                                0,
                                0
                            );
                        } else {
                            val = crangestat_compute_window_stat_simple(
                                stype,
                                source_data[src_idx],
                                win_start,
                                win_end,
                                0,
                                0
                            );
                        }

                        /* Apply result to ALL observations in this run */
                        for (size_t i = run_start; i < run_end; i++) {
                            result_data[s][i] = val;
                        }
                    }
                }

                free(run_starts);
            } else {
standard_obs_loop:;
                /* Rebuild O(1) structures for this group (reuse pre-allocated buffers) */
                int group_has_o1 = use_o1_queries && (g_count >= CRANGESTAT_PREFIX_SUM_THRESHOLD);
                if (group_has_o1) {
                    for (size_t v = 0; v < nsource; v++) {
                        if (needs_prefix[v]) {
                            prefix_arr[v].nobs = g_count;
                            crangestat_build_prefix_arrays(source_data[v], g_start, g_count, &prefix_arr[v]);
                        }
                        if (needs_sparse[v]) {
                            crangestat_rebuild_sparse_table(&sparse_tbl[v], source_data[v], g_start,
                                                 g_count, SV_missval);
                        }
                    }
                }

                /*
                 * OPTIMIZATION: Chunked parallel sliding window
                 *
                 * Each thread gets a contiguous chunk of non-missing observations.
                 * Each thread initializes its window with ONE binary search for its
                 * first observation, then slides the pointers forward (O(1) amortized).
                 * Total: O(P * log N + N) instead of O(N * log N).
                 */

                /* Find boundary before missing keys (they sort to end) */
                size_t g_end_nonmiss = g_end;
                while (g_end_nonmiss > g_start && key_data[g_end_nonmiss - 1] >= SV_missval)
                    g_end_nonmiss--;

                size_t nonmiss_count = g_end_nonmiss - g_start;

                #ifdef _OPENMP
                #pragma omp parallel num_threads(num_threads) if(nonmiss_count > CRANGESTAT_PARALLEL_THRESHOLD)
                #endif
                {
                    #ifdef _OPENMP
                    int tid = omp_get_thread_num();
                    int nthreads_actual = omp_get_num_threads();
                    #else
                    int tid = 0;
                    int nthreads_actual = 1;
                    #endif

                    double *work = work_arrays[tid];

                    /* Each thread gets a contiguous chunk */
                    size_t chunk = (nonmiss_count + (size_t)nthreads_actual - 1) / (size_t)nthreads_actual;
                    size_t my_start = g_start + (size_t)tid * chunk;
                    size_t my_end = my_start + chunk;
                    if (my_start > g_end_nonmiss) my_start = g_end_nonmiss;
                    if (my_end > g_end_nonmiss) my_end = g_end_nonmiss;

                    if (my_start < my_end) {
                        /* Initialize window with ONE binary search for first obs */
                        double first_key = key_data[my_start];
                        double first_low = config.low_is_missing ? -DBL_MAX : first_key + config.interval_low;
                        double first_high = config.high_is_missing ? DBL_MAX : first_key + config.interval_high;

                        size_t win_lo = crangestat_lower_bound(key_data, g_start, g_end_nonmiss, first_low);
                        size_t win_hi = crangestat_upper_bound(key_data, g_start, g_end_nonmiss, first_high);

                        for (size_t i = my_start; i < my_end; i++) {
                            double key_val = key_data[i];
                            double low_bound = config.low_is_missing ? -DBL_MAX : key_val + config.interval_low;
                            double high_bound = config.high_is_missing ? DBL_MAX : key_val + config.interval_high;

                            /* Slide window pointers forward */
                            while (win_lo < g_end_nonmiss && key_data[win_lo] < low_bound) win_lo++;
                            while (win_hi < g_end_nonmiss && key_data[win_hi] <= high_bound) win_hi++;

                            size_t win_start = win_lo;
                            size_t win_end = win_hi;

                            /* Compute each statistic for this observation */
                            for (size_t s = 0; s < nstats; s++) {
                                int src_idx = specs[s].source_var_idx;
                                if (src_idx < 0 || src_idx >= (int)nsource) continue;

                                double val;
                                stat_type stype = specs[s].type;

                                if (crangestat_stat_needs_work_array(stype)) {
                                    val = crangestat_compute_window_stat_complex(
                                        stype,
                                        source_data[src_idx],
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self,
                                        work
                                    );
                                } else if (group_has_o1 && needs_prefix[src_idx] &&
                                           crangestat_stat_uses_prefix_sums(stype)) {
                                    val = crangestat_compute_stat_with_prefix(
                                        stype,
                                        &prefix_arr[src_idx],
                                        source_data[src_idx],
                                        g_start,
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self
                                    );
                                } else if (group_has_o1 && needs_sparse[src_idx] &&
                                           crangestat_stat_uses_sparse_table(stype)) {
                                    val = crangestat_compute_minmax_with_sparse(
                                        stype,
                                        &sparse_tbl[src_idx],
                                        source_data[src_idx],
                                        g_start,
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self
                                    );
                                } else {
                                    val = crangestat_compute_window_stat_simple(
                                        stype,
                                        source_data[src_idx],
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self
                                    );
                                }

                                result_data[s][i] = val;
                            }
                        }
                    }
                }
            }
        }

        /* Cleanup O(1) data structures (once, after all groups) */
        if (use_o1_queries) {
            for (size_t v = 0; v < nsource; v++) {
                if (needs_prefix[v]) {
                    free(prefix_arr[v].prefix_sum);
                    free(prefix_arr[v].prefix_sum2);
                    free(prefix_arr[v].prefix_count);
                }
                if (needs_sparse[v]) {
                    crangestat_free_sparse_table(&sparse_tbl[v]);
                }
            }
            free(needs_prefix);
            free(needs_sparse);
            free(prefix_arr);
            free(sparse_tbl);
        }
    }

    t_compute = ctools_timer_seconds() - compute_start;

    /* === Store Phase === */
    double store_start = ctools_timer_seconds();

    /*
     * OPTIMIZATION: Parallel result storage
     * Store multiple result variables in parallel - each variable is independent.
     *
     * When fwd_perm is available, iterate in sorted order for sequential data[] reads.
     * fwd_perm[sorted_idx] = original_idx, so data[sorted_idx] is sequential.
     */
    int store_rc = 0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads) reduction(max:store_rc) if(nstats > 1)
    #endif
    for (size_t s = 0; s < nstats; s++) {
        ST_retcode local_rc = 0;
        int result_idx = specs[s].result_var_idx;
        double *data = result_data[s];

        if (fwd_perm) {
            /* Sequential reads from data[sorted_idx], scattered writes via obs_map */
            for (size_t sorted_idx = 0; sorted_idx < nobs; sorted_idx++) {
                perm_idx_t orig_idx = fwd_perm[sorted_idx];
                local_rc = SF_vstore(result_idx, (ST_int)obs_map[orig_idx], data[sorted_idx]);
                if (local_rc) break;
            }
        } else {
            for (size_t i = 0; i < nobs; i++) {
                local_rc = SF_vstore(result_idx, (ST_int)obs_map[i], data[inv_perm[i]]);
                if (local_rc) break;
            }
        }
        if (local_rc > store_rc) store_rc = local_rc;
    }

    t_store = ctools_timer_seconds() - store_start;

    /* === Cleanup === */
    for (int t = 0; t < num_threads; t++) {
        free(work_arrays[t]);
    }
    free(work_arrays);
    free(groups);
    free(inv_perm);
    if (fwd_perm) free(fwd_perm);

    /* Free result arrays (allocated separately from stata_data) */
    for (size_t s = 0; s < nstats; s++) free(result_data[s]);
    free(result_data);

    /* Free pointer arrays (not the underlying data, which is in stata_data) */
    free(source_data);

    /* Free the main data structure (includes key_data, source vars, by vars) */
    ctools_filtered_data_free(&filtered);
    free(load_indices);

    /* Free remaining allocations */
    free(source_indices);
    free(result_indices);
    free(specs);
    if (by_indices) free(by_indices);

    if (store_rc) return store_rc;

    /* Store timing scalars */
    double total = ctools_timer_seconds() - t_start;
    SF_scal_save("_crangestat_time_load", t_load);
    SF_scal_save("_crangestat_time_sort", t_sort);
    SF_scal_save("_crangestat_time_groups", t_groups);
    SF_scal_save("_crangestat_time_compute", t_compute);
    SF_scal_save("_crangestat_time_store", t_store);
    SF_scal_save("_crangestat_time_total", total);
    SF_scal_save("_crangestat_nobs", (double)nobs);
    SF_scal_save("_crangestat_nstats", (double)nstats);
    SF_scal_save("_crangestat_ngroups", (double)ngroups);

    CTOOLS_SAVE_THREAD_INFO("_crangestat");

    /* Verbose output handled by .ado file using stored scalars */

    return 0;
}
