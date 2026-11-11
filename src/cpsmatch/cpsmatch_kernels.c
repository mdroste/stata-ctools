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

/* The controls within distance r of p, as psmatch2 tests it: |c - p| <= r.
 * Over controls sorted by score the rounded distance falls, then rises, so
 * the matches are one run, found by two binary searches. */
void cpsmatch_distance_window(const sorted_control_t *sorted, size_t n, double p,
                                     double r, size_t *start, size_t *end)
{
    size_t lo = 0, hi = n;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        double c = sorted[mid].pscore;
        if (c < p && !(fabs(c - p) <= r))
            lo = mid + 1;
        else
            hi = mid;
    }
    *start = lo;
    hi = n;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        double c = sorted[mid].pscore;
        if (c <= p || fabs(c - p) <= r)
            lo = mid + 1;
        else
            hi = mid;
    }
    *end = lo;
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

    cpsmatch_distance_window(sorted_controls, n_controls, pscore_treated, radius, &info->start,
                             &info->end);
    info->n_matches = (int)(info->end - info->start);
    info->on_support = (info->n_matches > 0) ? 1 : 0;
}
