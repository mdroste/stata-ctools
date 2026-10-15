#include "cpsmatch_internal.h"

double cpsmatch_kernel_epan(double u)
{
    if (fabs(u) > 1.0)
        return 0.0;
    return 0.75 * (1.0 - u * u);
}

double cpsmatch_kernel_normal(double u)
{
    return exp(-0.5 * u * u) / sqrt(2.0 * M_PI);
}

double cpsmatch_kernel_biweight(double u)
{
    if (fabs(u) > 1.0)
        return 0.0;
    double v = 1.0 - u * u;
    return 0.9375 * v * v;
}

double cpsmatch_kernel_uniform(double u)
{
    if (fabs(u) > 1.0)
        return 0.0;
    return 0.5;
}

double cpsmatch_kernel_tricube(double u)
{
    if (fabs(u) > 1.0)
        return 0.0;
    double v = 1.0 - fabs(u * u * u);
    return 70.0 / 81.0 * v * v * v;
}

double cpsmatch_apply_kernel(double u, kernel_type_t kernel)
{
    switch (kernel)
    {
    case KERNEL_EPAN:
        return cpsmatch_kernel_epan(u);
    case KERNEL_NORMAL:
        return cpsmatch_kernel_normal(u);
    case KERNEL_BIWEIGHT:
        return cpsmatch_kernel_biweight(u);
    case KERNEL_UNIFORM:
        return cpsmatch_kernel_uniform(u);
    case KERNEL_TRICUBE:
        return cpsmatch_kernel_tricube(u);
    default:
        return cpsmatch_kernel_epan(u);
    }
}

void cpsmatch_match_workspace_init(match_workspace *ws, size_t initial_cap)
{
    ws->candidates = malloc(initial_cap * sizeof(idx_dist_pair));
    ws->candidates_cap = ws->candidates ? initial_cap : 0;
    ws->match_ids = malloc(initial_cap * sizeof(size_t));
    ws->match_weights = malloc(initial_cap * sizeof(double));
    ws->match_cap = (ws->match_ids && ws->match_weights) ? initial_cap : 0;
}

void cpsmatch_match_workspace_free(match_workspace *ws)
{
    free(ws->candidates);
    free(ws->match_ids);
    free(ws->match_weights);
    ws->candidates = NULL;
    ws->match_ids = NULL;
    ws->match_weights = NULL;
    ws->candidates_cap = 0;
    ws->match_cap = 0;
}

int cpsmatch_match_workspace_ensure_candidates(match_workspace *ws, size_t needed)
{
    if (needed <= ws->candidates_cap)
        return 0;
    size_t new_cap = ws->candidates_cap ? ws->candidates_cap : MATCH_WORKSPACE_INITIAL_CAP;
    while (new_cap < needed)
        new_cap *= 2;
    idx_dist_pair *p = realloc(ws->candidates, new_cap * sizeof(idx_dist_pair));
    if (!p)
        return -1;
    ws->candidates = p;
    ws->candidates_cap = new_cap;
    return 0;
}

int cpsmatch_match_workspace_ensure_match(match_workspace *ws, size_t needed)
{
    if (needed <= ws->match_cap)
        return 0;
    size_t new_cap = ws->match_cap ? ws->match_cap : MATCH_WORKSPACE_INITIAL_CAP;
    while (new_cap < needed)
        new_cap *= 2;
    size_t *p1 = realloc(ws->match_ids, new_cap * sizeof(size_t));
    double *p2 = realloc(ws->match_weights, new_cap * sizeof(double));
    if (!p1 || !p2)
    {
        if (p1)
            ws->match_ids = p1;
        if (p2)
            ws->match_weights = p2;
        return -1;
    }
    ws->match_ids = p1;
    ws->match_weights = p2;
    ws->match_cap = new_cap;
    return 0;
}

int cpsmatch_compare_by_pscore(const void *a, const void *b)
{
    const sorted_control_t *pa = (const sorted_control_t *)a;
    const sorted_control_t *pb = (const sorted_control_t *)b;
    if (pa->pscore < pb->pscore)
        return -1;
    if (pa->pscore > pb->pscore)
        return 1;
    /* Tie-breaker: sort by Stata observation index (ascending) */
    if (pa->obs_idx < pb->obs_idx)
        return -1;
    if (pa->obs_idx > pb->obs_idx)
        return 1;
    return 0;
}

int cpsmatch_compare_by_pscore_desc(const void *a, const void *b)
{
    const sorted_control_t *pa = (const sorted_control_t *)a;
    const sorted_control_t *pb = (const sorted_control_t *)b;
    if (pa->pscore > pb->pscore)
        return -1;
    if (pa->pscore < pb->pscore)
        return 1;
    /* Tie-breaker: sort by Stata observation index (ascending) */
    if (pa->obs_idx < pb->obs_idx)
        return -1;
    if (pa->obs_idx > pb->obs_idx)
        return 1;
    return 0;
}

size_t cpsmatch_binary_search_left(const sorted_control_t *sorted, size_t n, double target)
{
    size_t lo = 0, hi = n;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (sorted[mid].pscore < target)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid;
        }
    }
    return lo;
}

size_t cpsmatch_binary_search_right(const sorted_control_t *sorted, size_t n, double target)
{
    size_t lo = 0, hi = n;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (sorted[mid].pscore <= target)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid;
        }
    }
    return lo;
}

int cpsmatch_find_nearest_neighbors_fast(double pscore_treated,
                                         const sorted_control_t *sorted_controls,
                                         size_t *control_indices, size_t n_controls,
                                         int n_neighbors, double caliper, int handle_ties,
                                         match_workspace *ws, match_result *result)
{
    if (n_controls == 0)
    {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Binary search to find insertion point */
    size_t insert_pos = cpsmatch_binary_search_left(sorted_controls, n_controls, pscore_treated);

    /* For k-NN, we need to find k nearest neighbors around insert_pos
     * Use two pointers expanding outward from insert_pos */

    /* Ensure workspace has room for initial candidates */
    size_t initial_cap = (size_t)n_neighbors * 2 + 10;
    if (initial_cap > n_controls)
        initial_cap = n_controls;
    if (cpsmatch_match_workspace_ensure_candidates(ws, initial_cap) != 0)
        return -1;

    int n_candidates = 0;

    /* Two-pointer expansion from insert_pos */
    ssize_t left = (ssize_t)insert_pos - 1;
    size_t right = insert_pos;

    /* First, find the k nearest neighbors */
    while (n_candidates < n_neighbors && (left >= 0 || right < n_controls))
    {
        double dist_left =
            (left >= 0) ? fabs(pscore_treated - sorted_controls[left].pscore) : DBL_MAX;
        double dist_right =
            (right < n_controls) ? fabs(pscore_treated - sorted_controls[right].pscore) : DBL_MAX;
        /* Choose by which sides remain, not by the DBL_MAX sentinel, which a
         * real distance can equal or exceed */
        int take_left = left >= 0 && (right >= n_controls || dist_left <= dist_right);

        if (take_left)
        {
            if (caliper > 0.0 && dist_left > caliper)
            {
                left = -1; /* No more valid on left */
            }
            else
            {
                ws->candidates[n_candidates].idx = sorted_controls[left].orig_idx;
                ws->candidates[n_candidates].dist = dist_left;
                n_candidates++;
                left--;
            }
        }
        else
        {
            if (caliper > 0.0 && dist_right > caliper)
            {
                right = n_controls; /* No more valid on right */
            }
            else
            {
                ws->candidates[n_candidates].idx = sorted_controls[right].orig_idx;
                ws->candidates[n_candidates].dist = dist_right;
                n_candidates++;
                right++;
            }
        }
    }

    if (n_candidates == 0)
    {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Handle ties: include all observations at the same distance as the k-th */
    if (handle_ties && n_candidates >= n_neighbors)
    {
        double last_dist = ws->candidates[n_candidates - 1].dist;

        /* Continue expanding to find ties */
        while (left >= 0 || right < n_controls)
        {
            double dist_left =
                (left >= 0) ? fabs(pscore_treated - sorted_controls[left].pscore) : DBL_MAX;
            double dist_right = (right < n_controls)
                                    ? fabs(pscore_treated - sorted_controls[right].pscore)
                                    : DBL_MAX;

            double min_dist = (dist_left < dist_right) ? dist_left : dist_right;
            if (min_dist > last_dist + 1e-12)
                break; /* No more ties */

            /* Grow workspace if needed */
            if ((size_t)n_candidates >= ws->candidates_cap)
            {
                if (cpsmatch_match_workspace_ensure_candidates(ws, (size_t)n_candidates + 1) != 0)
                    return -1;
            }

            int take_left = left >= 0 && (right >= n_controls || dist_left <= dist_right);
            if (take_left && dist_left <= last_dist + 1e-12)
            {
                ws->candidates[n_candidates].idx = sorted_controls[left].orig_idx;
                ws->candidates[n_candidates].dist = dist_left;
                n_candidates++;
                left--;
            }
            else if (right < n_controls && dist_right <= last_dist + 1e-12)
            {
                ws->candidates[n_candidates].idx = sorted_controls[right].orig_idx;
                ws->candidates[n_candidates].dist = dist_right;
                n_candidates++;
                right++;
            }
            else
            {
                break;
            }
        }
    }

    /* Ensure workspace match buffers are large enough */
    if (cpsmatch_match_workspace_ensure_match(ws, (size_t)n_candidates) != 0)
        return -1;

    double total_weight = (double)n_candidates;
    for (int i = 0; i < n_candidates; i++)
    {
        ws->match_ids[i] = control_indices[ws->candidates[i].idx];
        ws->match_weights[i] = 1.0 / total_weight;
    }

    /* Point result at workspace buffers (no ownership transfer) */
    result->match_ids = ws->match_ids;
    result->match_weights = ws->match_weights;
    result->n_matches = n_candidates;
    result->on_support = 1;

    return 0;
}

void cpsmatch_find_radius_range_fast(double pscore_treated, const sorted_control_t *sorted_controls,
                                     size_t n_controls, double radius, radius_match_info *info)
{
    if (n_controls == 0 || radius <= 0.0)
    {
        info->n_matches = 0;
        info->on_support = 0;
        return;
    }

    /* Binary search to find range [low_target, high_target] */
    double low_target = pscore_treated - radius;
    double high_target = pscore_treated + radius;

    info->start = cpsmatch_binary_search_left(sorted_controls, n_controls, low_target);
    info->end = cpsmatch_binary_search_right(sorted_controls, n_controls, high_target);
    info->n_matches = (int)(info->end - info->start);
    info->on_support = (info->n_matches > 0) ? 1 : 0;
}

int cpsmatch_find_kernel_matches_fast(double pscore_treated,
                                      const sorted_control_t *sorted_controls,
                                      size_t *control_indices, size_t n_controls, double bandwidth,
                                      kernel_type_t kernel, double caliper, match_workspace *ws,
                                      match_result *result)
{
    if (n_controls == 0 || bandwidth <= 0.0)
    {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Determine search range based on kernel type and caliper */
    double search_radius = bandwidth; /* Most kernels have support [-1, 1] */
    if (kernel == KERNEL_NORMAL)
    {
        search_radius = 6.0 * bandwidth; /* Normal has infinite support, use 6 SD for precision */
    }
    if (caliper > 0.0 && caliper < search_radius)
    {
        search_radius = caliper;
    }

    /* Binary search to find range */
    double low_target = pscore_treated - search_radius;
    double high_target = pscore_treated + search_radius;

    size_t start = cpsmatch_binary_search_left(sorted_controls, n_controls, low_target);
    size_t end = cpsmatch_binary_search_right(sorted_controls, n_controls, high_target);

    if (start >= end)
    {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Ensure workspace can hold max possible matches */
    size_t range_size = end - start;
    if (cpsmatch_match_workspace_ensure_match(ws, range_size) != 0)
        return -1;

    /* Single pass: fill match buffers and accumulate total weight */
    int n_matches = 0;
    double total_weight = 0.0;

    for (size_t i = start; i < end; i++)
    {
        double dist = fabs(pscore_treated - sorted_controls[i].pscore);
        if (caliper > 0.0 && dist > caliper)
            continue;

        double u = dist / bandwidth;
        double w = cpsmatch_apply_kernel(u, kernel);
        if (w > 0.0)
        {
            ws->match_ids[n_matches] = control_indices[sorted_controls[i].orig_idx];
            ws->match_weights[n_matches] = w;
            total_weight += w;
            n_matches++;
        }
    }

    if (n_matches == 0)
    {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Normalize weights in-place */
    for (int i = 0; i < n_matches; i++)
    {
        ws->match_weights[i] /= total_weight;
    }

    /* Point result at workspace buffers (no ownership transfer) */
    result->match_ids = ws->match_ids;
    result->match_weights = ws->match_weights;
    result->n_matches = n_matches;
    result->on_support = 1;

    return 0;
}
