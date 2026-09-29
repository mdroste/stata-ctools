stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered,
                                            int var_idx, perm_idx_t *obs_map)
{
    if (!n_filtered) return STATA_OK;
    if (!values) return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = validate_variable(var_idx, 0);
    if (rc) return rc;
    if (!obs_map) return STATA_ERR_INVALID_INPUT;
    size_t available = (size_t)SF_nobs();
    size_t min_obs = available, max_obs = 0;
    int increasing = 1, decreasing = 1;
    for (size_t i = 0; i < n_filtered; i++) {
        size_t obs = (size_t)obs_map[i];
        if (obs < 1 || obs > available) return STATA_ERR_INVALID_INPUT;
        if (i && obs_map[i] <= obs_map[i - 1]) increasing = 0;
        if (i && obs_map[i] >= obs_map[i - 1]) decreasing = 0;
        if (obs < min_obs) min_obs = obs;
        if (obs > max_obs) max_obs = obs;
    }
    int distinct = increasing || decreasing;
    #ifdef _OPENMP
    if (!distinct && n_filtered >= MIN_OBS_PER_THREAD * 2 &&
        ctools_get_max_threads() > 1) {
        /* by-group commands supply unique permutations. Prove uniqueness
         * with bounded scratch rather than serializing their writes. Sparse
         * ranges or allocation failure safely retain sequential behavior. */
        size_t bytes = (max_obs - min_obs) / 8 + 1;
        if (bytes / sizeof(perm_idx_t) <= n_filtered) {
            unsigned char *seen = (unsigned char *)calloc(bytes, 1);
            if (seen) {
                distinct = 1;
                for (size_t i = 0; i < n_filtered; i++) {
                    size_t bit = (size_t)obs_map[i] - min_obs;
                    unsigned char mask = (unsigned char)(1u << (bit & 7));
                    if (seen[bit >> 3] & mask) { distinct = 0; break; }
                    seen[bit >> 3] |= mask;
                }
                free(seen);
            }
        }
    }
    #endif
    ST_IIID store_fn = IO_VSTORE_FN;
    /* Repeated destinations require input-order last-write-wins. */
    if (!distinct) {
        for (size_t i = 0; i < n_filtered; i++)
            if (store_fn(var_idx, (ST_int)obs_map[i], values[i]))
                return store_status(STATA_ERR_STATA_WRITE);
        return STATA_OK;
    }
    atomic_int error = 0;
    #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < n_filtered; i++) {
        if (store_fn(var_idx, (ST_int)obs_map[i], values[i])) record_io_error(&error, STATA_ERR_STATA_WRITE);
    }
    return store_status(atomic_load(&error));
}
