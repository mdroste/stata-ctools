#include "ctools_transport_internal.h"

stata_retcode ctools_io_build_obs_map(size_t obs_start, size_t obs_end, int flags,
                                      perm_idx_t **obs_map_out, size_t *n_filtered_out,
                                      size_t *n_range_out, int *was_filtered_out)
{
    ST_int obs1, obs2;
    size_t n_range, n_filtered;
    perm_idx_t *obs_map = NULL;
    int is_identity = 0;
    int skip_if_check = (flags & CTOOLS_LOAD_SKIP_IF) != 0;

    if (obs_start > INT_MAX || obs_end > INT_MAX || (obs_start && obs_end && obs_end < obs_start))
        return STATA_ERR_INVALID_INPUT;

    /* Resolve observation range */
    if (obs_start == 0 || obs_end == 0)
    {
        obs1 = SF_in1();
        obs2 = SF_in2();
        if (obs1 < 1 || obs2 < obs1)
        {
            *obs_map_out = NULL;
            *n_filtered_out = 0;
            *n_range_out = 0;
            *was_filtered_out = 0;
            return STATA_OK; /* Empty but valid */
        }
    }
    else
    {
        obs1 = (ST_int)obs_start;
        obs2 = (ST_int)obs_end;
    }

    if (ctools_io_validate_range((size_t)obs1, (size_t)(obs2 - obs1 + 1)))
        return STATA_ERR_INVALID_INPUT;
    n_range = (size_t)(obs2 - obs1 + 1);

    /* Fast path: skip SF_ifobs checks entirely */
    if (skip_if_check)
    {
        is_identity = 1;
        n_filtered = n_range;
        obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
        if (obs_map == NULL)
            return STATA_ERR_MEMORY;
        for (size_t i = 0; i < n_filtered; i++)
        {
            obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
        }
        *obs_map_out = obs_map;
        *n_filtered_out = n_filtered;
        *n_range_out = n_range;
        *was_filtered_out = 0;
        return STATA_OK;
    }

    /* Quick identity detection for small ranges */
    if (n_range <= 32)
    {
        is_identity = ctools_io_probe_for_identity(obs1, obs2, n_range);
        if (is_identity)
        {
            n_filtered = n_range;
            obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
            if (obs_map == NULL)
                return STATA_ERR_MEMORY;
            for (size_t i = 0; i < n_filtered; i++)
            {
                obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
            }
            *obs_map_out = obs_map;
            *n_filtered_out = n_filtered;
            *n_range_out = n_range;
            *was_filtered_out = 0;
            return STATA_OK;
        }
    }

    /* Scan the all-selected prefix without allocating scratch. Once a row is
     * excluded, remember the remaining selection bits so the expensive SPI
     * predicate is called only once per observation. At most N/8 extra bytes
     * are live, even for a very sparse sample. Allocation failure retains the
     * two-pass path, rather than making scratch memory a new requirement. */
    ST_BI selected_fn = (_stata_)->selobs;
    size_t prefix = 0;
    while (prefix < n_range && selected_fn(obs1 + (ST_int)prefix))
        prefix++;
    n_filtered = prefix;
    uint64_t *selected = NULL;
    size_t words = 0;
    if (prefix < n_range)
    {
        words = (n_range - prefix + 63) / 64;
        selected = (uint64_t *)calloc(words, sizeof(uint64_t));
        for (size_t i = prefix + 1; i < n_range; i++)
        {
            if (selected_fn(obs1 + (ST_int)i))
            {
                n_filtered++;
                if (selected)
                    selected[(i - prefix) / 64] |= UINT64_C(1) << ((i - prefix) % 64);
            }
        }
    }

    if (n_filtered == 0)
    {
        free(selected);
        *obs_map_out = NULL;
        *n_filtered_out = 0;
        *n_range_out = n_range;
        *was_filtered_out = 1;
        return STATA_OK;
    }

    is_identity = (n_filtered == n_range);
    obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
    if (obs_map == NULL)
    {
        free(selected);
        return STATA_ERR_MEMORY;
    }

    for (size_t i = 0; i < prefix; i++)
        obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
    size_t idx = prefix;
    if (selected)
    {
        for (size_t w = 0; w < words; w++)
        {
            uint64_t bits = selected[w];
            while (bits)
            {
                unsigned bit = 0;
#if defined(__GNUC__) || defined(__clang__)
                bit = (unsigned)__builtin_ctzll(bits);
#else
                while (!(bits & (UINT64_C(1) << bit)))
                    bit++;
#endif
                obs_map[idx++] = (perm_idx_t)((size_t)obs1 + prefix + w * 64 + bit);
                bits &= bits - 1;
            }
        }
        free(selected);
    }
    else if (!is_identity)
    {
        for (size_t i = prefix + 1; i < n_range; i++)
            if (selected_fn(obs1 + (ST_int)i))
                obs_map[idx++] = (perm_idx_t)(obs1 + (ST_int)i);
    }

    *obs_map_out = obs_map;
    *n_filtered_out = n_filtered;
    *n_range_out = n_range;
    *was_filtered_out = !is_identity;
    return STATA_OK;
}

int ctools_io_probe_for_identity(ST_int obs1, ST_int obs2, size_t n_range)
{
    /* Probe size: 16 at each end, or full range if small */
    size_t probe_size = 16;
    if (n_range <= 32)
    {
        /* Small range: check all observations */
        for (ST_int obs = obs1; obs <= obs2; obs++)
        {
            if (!SF_ifobs(obs))
            {
                return 0; /* Filtering detected */
            }
        }
        return 1; /* Identity confirmed */
    }

    /* Probe first 16 observations */
    for (size_t i = 0; i < probe_size; i++)
    {
        if (!SF_ifobs(obs1 + (ST_int)i))
        {
            return 0; /* Filtering detected */
        }
    }

    /* Probe last 16 observations */
    for (size_t i = 0; i < probe_size; i++)
    {
        if (!SF_ifobs(obs2 - (ST_int)i))
        {
            return 0; /* Filtering detected */
        }
    }

    return 1; /* Likely identity (will verify during full pass if needed) */
}

stata_retcode ctools_data_select(ctools_filtered_data *result, size_t obs_start, size_t obs_end,
                                 int flags)
{
    if (result == NULL)
        return STATA_ERR_INVALID_INPUT;
    ctools_filtered_data_init(result);
    size_t n_filtered = 0, n_range = 0;
    int was_filtered = 0;
    perm_idx_t *obs_map = NULL;
    stata_retcode rc = ctools_io_build_obs_map(obs_start, obs_end, flags, &obs_map, &n_filtered,
                                               &n_range, &was_filtered);
    if (rc != STATA_OK)
        return ctools_io_load_status(rc);
    result->obs_map = obs_map;
    result->n_range = n_range;
    result->was_filtered = was_filtered;
    result->data.nobs = n_filtered;
    return STATA_OK;
}

stata_retcode ctools_io_plan_load(ctools_load_plan *plan, int *var_indices, size_t nvars,
                                  size_t obs_start, size_t obs_end, int flags,
                                  const int *str_widths)
{
    memset(plan, 0, sizeof(*plan));
    size_t n_range = 0, n_filtered = 0;
    perm_idx_t *obs_map = NULL;
    int was_filtered = 0;
    int *auto_indices = NULL, *auto_str_widths = NULL;
    /* "Load all variables" mode */
    if (var_indices == NULL)
    {
        nvars = (size_t)SF_nvars();
        if (nvars == 0)
        {
            return STATA_OK;
        }
        auto_indices = (int *)malloc(nvars * sizeof(int));
        if (auto_indices == NULL)
        {
            return STATA_ERR_MEMORY;
        }
        for (size_t i = 0; i < nvars; i++)
        {
            auto_indices[i] = (int)(i + 1);
        }
        var_indices = auto_indices;
    }
    else if (nvars == 0)
    {
        return STATA_ERR_INVALID_INPUT;
    }

    for (size_t j = 0; j < nvars; j++)
    {
        stata_retcode valid = ctools_io_validate_variable(var_indices[j], -1);
        if (valid)
        {
            free(auto_indices);
            return valid;
        }
    }

    /* Auto-detect string widths from Stata local */
    if (str_widths == NULL)
    {
        auto_str_widths = ctools_io_ctools_read_strw_from_stata((size_t)SF_nvars());
        str_widths = auto_str_widths;
    }

    /* Build observation map (shared filtering logic) */
    stata_retcode rc = ctools_io_build_obs_map(obs_start, obs_end, flags, &obs_map, &n_filtered,
                                               &n_range, &was_filtered);
    if (rc != STATA_OK)
    {
        free(auto_str_widths);
        free(auto_indices);
        return ctools_io_load_status(rc);
    }

    *plan = (ctools_load_plan){var_indices, auto_indices, str_widths, auto_str_widths, nvars,
                               obs_map,     n_range,      n_filtered, was_filtered};
    return STATA_OK;
}
