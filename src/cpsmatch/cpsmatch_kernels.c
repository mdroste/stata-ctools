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

/* ========================================================================
 * Kernel matching through moment sums
 *
 * The Epanechnikov, biweight, uniform and tricube weights are polynomials
 * in u = q - x on either side of the query q (|u|^3 for tricube), so the
 * sums a kernel pass needs over a window of sorted items, sum_x K(|q - x|)
 * times 1, an outcome or an item weight, follow from the window's sums of
 * (z - x)^i for an anchor z near q: K(|q - x|) = sum_k c_k (d + (z - x))^k
 * with d = q - z. The windows (the items with a positive weight, found with
 * the very test the direct pass makes) slide with the sorted queries; the
 * sums are double-double and are formed again from the window whenever the
 * query has moved a bandwidth from the anchor, so every term stays within a
 * few bandwidths and the cancellation stays bounded. Results are the exact
 * sums correctly rounded but for some units in the 30th digit, where the
 * direct pass carries the rounding of its sequential double sums.
 * ======================================================================== */

#include "ctools_ols.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define CPSMATCH_KMAX_DEG 9

static inline dd_real kdd_quick(double a, double b)
{
    double s = a + b;
    return (dd_real){s, b - (s - a)};
}
static inline dd_real kdd_add(dd_real a, dd_real b)
{
    dd_real s = two_sum(a.hi, b.hi), t = two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = kdd_quick(s.hi, s.lo);
    s.lo += t.lo;
    return kdd_quick(s.hi, s.lo);
}
static inline dd_real kdd_neg(dd_real a)
{
    return (dd_real){-a.hi, -a.lo};
}
static inline dd_real kdd_mul(dd_real a, dd_real b)
{
    dd_real p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return kdd_quick(p.hi, p.lo);
}
static inline dd_real kdd_d(double x)
{
    return (dd_real){x, 0.0};
}
static inline dd_real kdd_div(dd_real a, dd_real b)
{
    double q1 = a.hi / b.hi;
    dd_real r = kdd_add(a, kdd_neg(kdd_mul(kdd_d(q1), b)));
    double q2 = r.hi / b.hi;
    r = kdd_add(r, kdd_neg(kdd_mul(kdd_d(q2), b)));
    double q3 = r.hi / b.hi;
    return kdd_add(kdd_quick(q1, q2), kdd_d(q3));
}

/* The direct pass's test: does an item at distance dist get a positive
 * weight? (|dif| <= bw, then the kernel) */
static inline int kmom_positive(double q, double x, double bw, kernel_type_t kernel)
{
    double dist = fabs(q - x);
    if (!(dist <= bw)) return 0;
    return cpsmatch_apply_kernel(dist / bw, kernel) > 0.0;
}

/* Polynomial coefficients in u for the items left of the query (u >= 0)
 * and right of it (u < 0); constants that normalization cancels dropped */
typedef struct
{
    int deg;
    dd_real left[CPSMATCH_KMAX_DEG + 1], right[CPSMATCH_KMAX_DEG + 1];
} kmom_poly;

static void kmom_poly_make(kmom_poly *p, kernel_type_t kernel, double bw)
{
    memset(p, 0, sizeof(*p));
    dd_real ib = kdd_div(kdd_d(1.0), kdd_d(bw)), pw[10];
    pw[0] = kdd_d(1.0);
    for (int i = 1; i < 10; i++) pw[i] = kdd_mul(pw[i - 1], ib);
    switch (kernel) {
    case KERNEL_EPAN: /* 1 - u^2/b^2 */
        p->deg = 2;
        p->left[0] = kdd_d(1.0);
        p->left[2] = kdd_neg(pw[2]);
        break;
    case KERNEL_BIWEIGHT: /* (1 - u^2/b^2)^2 */
        p->deg = 4;
        p->left[0] = kdd_d(1.0);
        p->left[2] = kdd_mul(kdd_d(-2.0), pw[2]);
        p->left[4] = pw[4];
        break;
    case KERNEL_UNIFORM:
        p->deg = 0;
        p->left[0] = kdd_d(1.0);
        break;
    default: /* tricube: (1 - |u|^3/b^3)^3, |u|^3 = u^3 on the left, -u^3 on the right */
        p->deg = 9;
        p->left[0] = p->right[0] = kdd_d(1.0);
        p->left[3] = kdd_mul(kdd_d(-3.0), pw[3]);
        p->right[3] = kdd_mul(kdd_d(3.0), pw[3]);
        p->left[6] = p->right[6] = kdd_mul(kdd_d(3.0), pw[6]);
        p->left[9] = kdd_neg(pw[9]);
        p->right[9] = pw[9];
        return;
    }
    for (int i = 0; i <= p->deg; i++) p->right[i] = p->left[i];
}

/* A window of items [l, r) and its sums of wt * (z - x)^i (and of
 * wt * (z - x)^i * y_k with the count of missing y_k) */
typedef struct
{
    size_t l, r;
    dd_real s[CPSMATCH_KMAX_DEG + 1];
    dd_real *y;    /* ny * (deg + 1) */
    size_t *miss;  /* ny */
} kmom_win;

typedef struct
{
    const double *x;          /* item positions, ascending */
    const dd_real *wt;        /* item weights (NULL: one) */
    const double *const *yv;  /* ny item values (in item order) */
    int ny, deg;
} kmom_items;

static void kmom_win_reset(kmom_win *w, const kmom_items *it, size_t at)
{
    w->l = w->r = at;
    memset(w->s, 0, sizeof(w->s));
    if (it->ny) {
        memset(w->y, 0, (size_t)it->ny * (size_t)(it->deg + 1) * sizeof(dd_real));
        memset(w->miss, 0, (size_t)it->ny * sizeof(size_t));
    }
}

static void kmom_win_item(kmom_win *w, const kmom_items *it, size_t j, double z, int sign)
{
    dd_real zeta = two_sum(z, -it->x[j]); /* exact */
    dd_real pw[CPSMATCH_KMAX_DEG + 1];
    pw[0] = it->wt ? it->wt[j] : kdd_d(1.0);
    for (int i = 1; i <= it->deg; i++) pw[i] = kdd_mul(pw[i - 1], zeta);
    for (int i = 0; i <= it->deg; i++) w->s[i] = kdd_add(w->s[i], sign > 0 ? pw[i] : kdd_neg(pw[i]));
    for (int k = 0; k < it->ny; k++) {
        double v = it->yv[k][j];
        if (!(v < SV_missval)) {
            w->miss[k] = sign > 0 ? w->miss[k] + 1 : w->miss[k] - 1;
            continue;
        }
        dd_real *yk = w->y + (size_t)k * (size_t)(it->deg + 1);
        for (int i = 0; i <= it->deg; i++) {
            dd_real t = kdd_mul(pw[i], kdd_d(v));
            yk[i] = kdd_add(yk[i], sign > 0 ? t : kdd_neg(t));
        }
    }
}

/* Move the window to [l2, r2) (l2 >= l, r2 >= r) */
static void kmom_win_move(kmom_win *w, const kmom_items *it, size_t l2, size_t r2, double z)
{
    if (l2 >= w->r) {
        kmom_win_reset(w, it, l2);
    } else {
        for (size_t j = w->l; j < l2; j++) kmom_win_item(w, it, j, z, -1);
        w->l = l2;
    }
    for (size_t j = w->r; j < r2; j++) kmom_win_item(w, it, j, z, +1);
    w->r = r2;
}

/* g_i = sum_k c_k binom(k, i) d^(k - i): the window sum is sum_i g_i s_i */
static void kmom_coefs(const dd_real *c, int deg, dd_real d, dd_real *g)
{
    static const double binom[10][10] = {
        {1}, {1, 1}, {1, 2, 1}, {1, 3, 3, 1}, {1, 4, 6, 4, 1}, {1, 5, 10, 10, 5, 1},
        {1, 6, 15, 20, 15, 6, 1}, {1, 7, 21, 35, 35, 21, 7, 1}, {1, 8, 28, 56, 70, 56, 28, 8, 1},
        {1, 9, 36, 84, 126, 126, 84, 36, 9, 1}};
    dd_real dp[CPSMATCH_KMAX_DEG + 1];
    dp[0] = kdd_d(1.0);
    for (int i = 1; i <= deg; i++) dp[i] = kdd_mul(dp[i - 1], d);
    for (int i = 0; i <= deg; i++) {
        dd_real sum = kdd_d(0.0);
        for (int k = i; k <= deg; k++) {
            if (c[k].hi == 0.0 && c[k].lo == 0.0) continue;
            sum = kdd_add(sum, kdd_mul(kdd_mul(c[k], kdd_d(binom[k][i])), dp[k - i]));
        }
        g[i] = sum;
    }
}

/* For each query q[t] (ascending): the positive window [s, e) of the
 * items, total = sum K, and with values ysum[t*ny + k] = sum K y_k and
 * ymiss[t*ny + k] = missing y_k in the window. work: per thread, 2 * ny *
 * (deg + 1) dd and 2 * ny counts. */
/* The weight before the constant normalization cancels, straight from the
 * distance (the kernel's own formula without its constant) */
static inline double kmom_free(double q, double x, double bw, kernel_type_t kernel)
{
    const double u = fabs(q - x) / bw;
    switch (kernel) {
    case KERNEL_EPAN:
        return 1.0 - u * u;
    case KERNEL_BIWEIGHT: {
        const double v = 1.0 - u * u;
        return v * v;
    }
    case KERNEL_UNIFORM:
        return 1.0;
    default: {
        const double v = 1.0 - u * u * u;
        return v * v * v;
    }
    }
}

/* A window's sums straight from its items: for the windows whose moment
 * sums cancel too much (all weights tiny, at the bandwidth's edge) */
static void kmom_direct(double q, const kmom_items *it, size_t s, size_t e, double bw,
                        kernel_type_t kernel, dd_real *total, dd_real *ysum)
{
    dd_real tot = kdd_d(0.0);
    for (int k = 0; k < it->ny; k++) ysum[k] = kdd_d(0.0);
    for (size_t j = s; j < e; j++) {
        dd_real w = kdd_d(kmom_free(q, it->x[j], bw, kernel));
        if (it->wt) w = kdd_mul(w, it->wt[j]);
        tot = kdd_add(tot, w);
        for (int k = 0; k < it->ny; k++) {
            const double v = it->yv[k][j];
            if (v < SV_missval) ysum[k] = kdd_add(ysum[k], kdd_mul(w, kdd_d(v)));
        }
    }
    *total = tot;
}

/* Sums whose terms exceed the result by more than this are formed again
 * from the window's items: the double-double sums then still carry some
 * 20 digits */
#define KMOM_CANCEL 1e8

static void kmom_pass(const double *q, size_t nq, const kmom_items *it, size_t nx,
                      const kmom_poly *poly, kernel_type_t kernel, double bw, int threads,
                      dd_real *work_y, size_t *work_miss, size_t *cnt, dd_real *total,
                      dd_real *ysum, size_t *ymiss)
{
    const int ny = it->ny, D = it->deg;
    const size_t chunks = nq < 4096 ? 1 : (size_t)threads * 8;
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if (threads > 1)
    for (size_t c = 0; c < chunks; c++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        const size_t t0 = nq * c / chunks, t1 = nq * (c + 1) / chunks;
        if (t0 >= t1) continue;
        kmom_win L, R;
        L.y = work_y + (size_t)tid * 2 * (size_t)ny * (size_t)(D + 1);
        R.y = L.y + (size_t)ny * (size_t)(D + 1);
        L.miss = work_miss + (size_t)tid * 2 * (size_t)ny;
        R.miss = L.miss + ny;
        /* the first query's window, by binary search */
        const double q0 = q[t0];
        size_t lo = 0, hi = nx;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (it->x[mid] < q0 && !kmom_positive(q0, it->x[mid], bw, kernel)) lo = mid + 1;
            else hi = mid;
        }
        size_t s = lo, e, m;
        lo = s, hi = nx;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (it->x[mid] <= q0 || kmom_positive(q0, it->x[mid], bw, kernel)) lo = mid + 1;
            else hi = mid;
        }
        e = lo;
        lo = 0, hi = nx;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (it->x[mid] <= q0) lo = mid + 1;
            else hi = mid;
        }
        m = lo;
        double z = q0;
        int fresh = 1;
        for (size_t t = t0; t < t1; t++) {
            const double qt = q[t];
            if (t > t0) {
                while (s < nx && it->x[s] < qt && !kmom_positive(qt, it->x[s], bw, kernel)) s++;
                if (e < s) e = s;
                while (e < nx && (it->x[e] <= qt || kmom_positive(qt, it->x[e], bw, kernel))) e++;
                while (m < nx && it->x[m] <= qt) m++;
            }
            size_t mc = m < s ? s : (m > e ? e : m);
            if (fresh || qt - z > bw) {
                /* the window's sums about a new anchor at this query */
                z = qt;
                kmom_win_reset(&L, it, s);
                kmom_win_reset(&R, it, mc);
                fresh = 0;
            }
            kmom_win_move(&L, it, s, mc, z);
            kmom_win_move(&R, it, mc, e, z);
            cnt[t] = e - s;
            dd_real d = two_sum(qt, -z), gl[CPSMATCH_KMAX_DEG + 1], gr[CPSMATCH_KMAX_DEG + 1];
            kmom_coefs(poly->left, D, d, gl);
            kmom_coefs(poly->right, D, d, gr);
            /* Each item is within two bandwidths of the anchor, so a term
             * g_i s_i is at most |g_i| (sum of item weights) (2 bw)^i: the
             * sum's cancellation is that bound over the sum. The outcome
             * sums share it, relative to the outcomes' size. */
            dd_real tot = kdd_d(0.0);
            double scale = 0.0, span = 1.0;
            const double wl = fabs(L.s[0].hi), wr = fabs(R.s[0].hi);
            for (int i = 0; i <= D; i++) {
                tot = kdd_add(tot, kdd_add(kdd_mul(gl[i], L.s[i]), kdd_mul(gr[i], R.s[i])));
                scale += (fabs(gl[i].hi) * wl + fabs(gr[i].hi) * wr) * span;
                span *= 2.0 * bw;
            }
            for (int k = 0; k < ny; k++) {
                const dd_real *yl = L.y + (size_t)k * (size_t)(D + 1);
                const dd_real *yr = R.y + (size_t)k * (size_t)(D + 1);
                dd_real sum = kdd_d(0.0);
                for (int i = 0; i <= D; i++)
                    sum = kdd_add(sum, kdd_add(kdd_mul(gl[i], yl[i]), kdd_mul(gr[i], yr[i])));
                ysum[t * (size_t)ny + k] = sum;
            }
            if (e > s && !(fabs(tot.hi) * KMOM_CANCEL > scale))
                kmom_direct(qt, it, s, e, bw, kernel, &tot, ny ? ysum + t * (size_t)ny : NULL);
            total[t] = tot;
            for (int k = 0; k < ny; k++) ymiss[t * (size_t)ny + k] = L.miss[k] + R.miss[k];
        }
    }
}

int cpsmatch_kernel_moments(kernel_type_t kernel, double bw, const double *tps, size_t nt,
                            const double *dps, size_t nd, const double *const *dy, int ny,
                            int threads, double *donor_weight, double *nn, int *on,
                            double *const *matched)
{
    if (threads < 1) threads = 1;
    kmom_poly poly;
    kmom_poly_make(&poly, kernel, bw);
    const int D = poly.deg;
    const size_t nyk = ny > 0 ? (size_t)ny : 1;
    size_t *cnt = malloc((nt ? nt : 1) * sizeof(size_t));
    dd_real *total = malloc((nt ? nt : 1) * sizeof(dd_real));
    dd_real *ysum = malloc((nt ? nt : 1) * nyk * sizeof(dd_real));
    size_t *ymiss = malloc((nt ? nt : 1) * nyk * sizeof(size_t));
    dd_real *inv = malloc((nt ? nt : 1) * sizeof(dd_real));
    dd_real *wsum = malloc((nd ? nd : 1) * sizeof(dd_real));
    size_t *dcnt = malloc((nd ? nd : 1) * sizeof(size_t));
    dd_real *work_y = malloc((size_t)threads * 2 * nyk * (size_t)(D + 1) * sizeof(dd_real));
    size_t *work_miss = malloc((size_t)threads * 2 * nyk * sizeof(size_t));
    int rc = 0;
    if (!cnt || !total || !ysum || !ymiss || !inv || !wsum || !dcnt || !work_y || !work_miss) {
        rc = 920;
        goto done;
    }
    /* the targets' windows of donors: totals and outcome sums */
    kmom_items donors = {dps, NULL, dy, ny, D};
    kmom_pass(tps, nt, &donors, nd, &poly, kernel, bw, threads, work_y, work_miss, cnt, total,
              ysum, ymiss);
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t t = 0; t < nt; t++) {
        on[t] = cnt[t] > 0;
        inv[t] = kdd_d(0.0);
        if (!on[t]) continue;
        nn[t] = (double)cnt[t];
        inv[t] = kdd_div(kdd_d(1.0), total[t]);
        for (int k = 0; k < ny; k++) {
            if (!matched || !matched[k]) continue;
            if (ymiss[t * (size_t)ny + k]) {
                matched[k][t] = SV_missval;
                continue;
            }
            dd_real mean = kdd_div(ysum[t * (size_t)ny + k], total[t]);
            matched[k][t] = mean.hi + mean.lo;
        }
    }
    /* the donors' weights: sum over the targets that weigh a donor of
     * K / (the target's total) */
    kmom_items targets = {tps, inv, NULL, 0, D};
    kmom_pass(dps, nd, &targets, nt, &poly, kernel, bw, threads, work_y, work_miss, dcnt, wsum,
              NULL, NULL);
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t j = 0; j < nd; j++) donor_weight[j] = dcnt[j] ? wsum[j].hi + wsum[j].lo : 0.0;
done:
    free(cnt);
    free(total);
    free(ysum);
    free(ymiss);
    free(inv);
    free(wsum);
    free(dcnt);
    free(work_y);
    free(work_miss);
    return rc;
}
