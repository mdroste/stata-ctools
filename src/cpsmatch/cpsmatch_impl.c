#include "cpsmatch_internal.h"
#include "cpsmatch_impl.h"
#include "cpsmatch_stats.h"
#include "cpsmatch_index.h"
#include "ctools_ols.h"

/* Propensity-score matching as psmatch2 4.1.0 performs it.
 *
 * Nearest-neighbor matching works on psmatch2's row order: the sample sorted
 * by support (on support first), _treated, pscore (negated with descending)
 * and observation. A row's position in that order is psmatch2's _id; the
 * supported controls and the supported treated are two blocks of it. A pass
 * matches each row of a target block, in order, to rows of a donor block:
 * - the treated to the controls (ATT); with ate first the controls to the
 *   treated; with ai(#) each block to itself, skipping the row itself;
 * - the nearest donor is found by scanning on from the previous one in the
 *   current direction, over donors not yet used when noreplace, while the
 *   distance does not increase (the first of equal distances wins); once a
 *   match is the last donor, the direction turns backward;
 * - with replacement, ties includes every distinct donor through the kth
 *   distance boundary. Without replacement (neighbor(1)), it retains the
 *   reference's forward scan over unused donors at the nearest distance;
 * - every list entry gets 1/(list length), and a match needs a distance below
 *   caliper() (strict); a target left without one is off support;
 * - without replacement the matched donors are used up and the scan resumes
 *   at the next unused donor; once it runs off the donors the remaining
 *   targets are off support. Only the nearest donor and its ties are taken
 *   that way: in the ai(#) pass the other # - 1 matches are the positions
 *   around the nearest, used or not, as in psmatch2 (the scalar
 *   _cpsmatch_ai_pass tells the ado file that the plugin does this).
 * With replacement and without ties the scan finds what a binary search
 * finds (the closer neighbor, the lower on equal distance, the first of a run
 * of equal scores), so those passes run in parallel.
 *
 * One call matches once for every outcome (the match depends on the score
 * alone) and, for the ai() correction, for every score-model column. Donor
 * weights add up in the targets' order, as a serial pass adds them; the
 * other outputs belong to one row each and are formed in parallel.
 *
 * Kernel and radius matching weigh every donor of the other arm within the
 * bandwidth; with ate the controls are matched first, and a control left
 * without a match is then no longer a donor for the treated, as in psmatch2. */

/* Stata's missing values, tested without a host callback per element */
static inline int psm_missing(double x, double missval)
{
    return !(x < missval);
}

/* Threads for a pass over n items */
static int psm_threads(size_t n)
{
    int threads = 1;
#ifdef _OPENMP
    threads = ctools_get_max_threads();
    if (threads < 1 || n < 16384)
        threads = 1;
#else
    (void)n;
#endif
    return threads;
}

/* gen byte: truncation toward zero, missing outside -127..100 */
static inline double psm_byte(double x)
{
    if (psm_missing(x, SV_missval))
        return SV_missval;
    double t = trunc(x);
    if (t < -127.0 || t > 100.0)
        return SV_missval;
    return t + 0.0;
}

/* One in-sample row in psmatch2's matching order */
typedef struct
{
    double arm;      /* stored _treated */
    double key;      /* pscore, negated with descending */
    perm_idx_t obs;  /* Stata observation */
    perm_idx_t row;  /* filtered row */
    int support;     /* support before matching */
} psm_order_t;

static int psm_compare_order(const void *a, const void *b)
{
    const psm_order_t *pa = (const psm_order_t *)a, *pb = (const psm_order_t *)b;
    if (pa->support != pb->support) return pa->support > pb->support ? -1 : 1;
    if (pa->arm != pb->arm) return pa->arm < pb->arm ? -1 : 1;
    if (pa->key != pb->key) return pa->key < pb->key ? -1 : 1;
    return (pa->obs > pb->obs) - (pa->obs < pb->obs);
}

/* A stable sort of a block by score: an LSD radix sort of order-preserving
 * 64-bit images of the keys, so equal scores keep the block's observation
 * order. */
typedef struct
{
    uint64_t key;
    uint32_t at;
} psm_radix_t;

static uint64_t psm_radix_key(double x)
{
    uint64_t bits;
    if (x == 0.0) x = 0.0; /* -0 sorts with +0, as they compare */
    memcpy(&bits, &x, sizeof(bits));
    return (bits >> 63) ? ~bits : bits | ((uint64_t)1 << 63);
}

enum { PSM_RADIX_BITS = 11, PSM_RADIX_BUCKETS = 1 << PSM_RADIX_BITS };

/* Sort keys[0..n) stably by key; returns the array holding the result
 * (keys or spare). Each thread counts the digits of a contiguous chunk and
 * scatters it in chunk order, so the result does not depend on the thread
 * count. hist holds threads * PSM_RADIX_BUCKETS entries. */
static psm_radix_t *psm_radix_sort(psm_radix_t *keys, psm_radix_t *spare, size_t n, size_t *hist,
                                   int threads)
{
    if (n < 2) return keys;
    if (threads < 1) threads = 1;
    psm_radix_t *src = keys, *dst = spare;
    int skip = 0, team_size = 1;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
#pragma omp single
        team_size = (int)team;
        const size_t lo = n * tid / team, hi = n * (tid + 1) / team;
        for (int shift = 0; shift < 64; shift += PSM_RADIX_BITS)
        {
            size_t *h = hist + tid * PSM_RADIX_BUCKETS;
            memset(h, 0, PSM_RADIX_BUCKETS * sizeof(size_t));
            for (size_t i = lo; i < hi; i++)
                h[(src[i].key >> shift) & (PSM_RADIX_BUCKETS - 1)]++;
#pragma omp barrier
#pragma omp single
            {
                skip = 0;
                size_t sum = 0;
                for (size_t b = 0; b < PSM_RADIX_BUCKETS && !skip; b++)
                {
                    size_t total = 0;
                    for (int t = 0; t < team_size; t++)
                        total += hist[(size_t)t * PSM_RADIX_BUCKETS + b];
                    if (total == n) skip = 1;
                }
                if (!skip)
                {
                    for (size_t b = 0; b < PSM_RADIX_BUCKETS; b++)
                        for (int t = 0; t < team_size; t++)
                        {
                            size_t c = hist[(size_t)t * PSM_RADIX_BUCKETS + b];
                            hist[(size_t)t * PSM_RADIX_BUCKETS + b] = sum;
                            sum += c;
                        }
                }
            }
            if (!skip)
            {
                for (size_t i = lo; i < hi; i++)
                    dst[h[(src[i].key >> shift) & (PSM_RADIX_BUCKETS - 1)]++] = src[i];
            }
#pragma omp barrier
#pragma omp single
            {
                if (!skip)
                {
                    psm_radix_t *swap = src;
                    src = dst;
                    dst = swap;
                }
            }
        }
    }
    return src;
}

static void psm_sort_block(psm_order_t *block, size_t n, psm_order_t *tmp, psm_radix_t *keys,
                           psm_radix_t *spare, size_t *hist, int threads)
{
    if (n < 2) return;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < n; i++)
    {
        keys[i].key = psm_radix_key(block[i].key);
        keys[i].at = (uint32_t)i;
    }
    keys = psm_radix_sort(keys, spare, n, hist, threads);
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < n; i++) tmp[i] = block[keys[i].at];
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < n; i++) block[i] = tmp[i];
}

/* The bins of psmatch2's order: supported rows coded below 0, the
 * supported controls, supported rows coded in (0, 1), the supported
 * treated, then the rest. */
static int psm_bin(const psm_order_t *r)
{
    if (!r->support) return 4;
    if (r->arm < 0.0) return 0;
    if (r->arm == 0.0) return 1;
    if (r->arm < 1.0) return 2;
    return r->arm == 1.0 ? 3 : 4;
}

/* Stable counting sort of order[0..n) into out by bin; bin_start gets the
 * six bin boundaries. hist holds threads * 5 entries. */
static void psm_bin_rows(const psm_order_t *order, size_t n, psm_order_t *out, size_t *bin_start,
                         size_t *hist, int threads)
{
    int team_size = 1;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
#pragma omp single
        team_size = (int)team;
        const size_t lo = n * tid / team, hi = n * (tid + 1) / team;
        size_t *h = hist + tid * 5;
        memset(h, 0, 5 * sizeof(size_t));
        for (size_t p = lo; p < hi; p++) h[psm_bin(&order[p])]++;
#pragma omp barrier
#pragma omp single
        {
            size_t sum = 0;
            for (int b = 0; b < 5; b++)
            {
                bin_start[b] = sum;
                for (int t = 0; t < team_size; t++)
                {
                    size_t c = hist[(size_t)t * 5 + b];
                    hist[(size_t)t * 5 + b] = sum;
                    sum += c;
                }
            }
            bin_start[5] = sum;
        }
        for (size_t p = lo; p < hi; p++) out[h[psm_bin(&order[p])]++] = order[p];
    }
}

/* Match lists of a pass: target t's donors (positions) are
 * item[start[t] .. start[t] + count[t]). */
typedef struct
{
    size_t *start;
    uint32_t *count;
    perm_idx_t *item;
    size_t len, cap;
} psm_lists;

static void psm_lists_free(psm_lists *l)
{
    free(l->start);
    free(l->count);
    free(l->item);
    memset(l, 0, sizeof(*l));
}

static int psm_lists_init(psm_lists *l, size_t n_targets, size_t capacity)
{
    memset(l, 0, sizeof(*l));
    l->start = ctools_safe_malloc2(n_targets ? n_targets : 1, sizeof(size_t));
    l->count = ctools_safe_calloc2(n_targets ? n_targets : 1, sizeof(uint32_t));
    l->item = ctools_safe_malloc2(capacity ? capacity : 1, sizeof(perm_idx_t));
    l->cap = capacity ? capacity : 1;
    if (!l->start || !l->count || !l->item) {
        psm_lists_free(l);
        return -1;
    }
    return 0;
}

static int psm_lists_push(psm_lists *l, ssize_t position)
{
    if (l->len == l->cap) {
        size_t cap = l->cap * 2;
        perm_idx_t *grown = realloc(l->item, cap * sizeof(perm_idx_t));
        if (!grown) return -1;
        l->item = grown;
        l->cap = cap;
    }
    l->item[l->len++] = (perm_idx_t)position;
    return 0;
}

/* psmatch2's next_unmatched: the next position in the scan direction, past
 * the target itself and, without replacement, past used donors. Positions
 * outside the sample are returned for the caller to reject. */
static ssize_t psm_step(ssize_t j, int forward, ssize_t self, int noreplace,
                        const unsigned char *used, ssize_t n_pos)
{
    for (;;) {
        j += forward ? 1 : -1;
        if (j < 0 || j >= n_pos) return j;
        if (j == self) continue;
        if (!noreplace || !used[j]) return j;
    }
}

/* The sequential pass over targets [i0, i1] and donors [j0, j1] (inclusive
 * positions; a range is empty when its end precedes its start). self: the
 * two ranges are one block, matched to itself. */
static int psm_scan(const double *ps, ssize_t n_pos, ssize_t i0, ssize_t i1, ssize_t j0,
                    ssize_t j1, int self, int neighbors, double cal, int ties, int noreplace,
                    unsigned char *used, psm_lists *out)
{
    memset(used, 0, (size_t)n_pos);
    ssize_t jm = j0, i = i0;
    int forward = 1;
    for (; i <= i1 && jm >= j0 && jm <= j1; i++) {
        size_t t = (size_t)(i - i0);
        ssize_t me = self ? i : -1;
        out->start[t] = out->len;
        out->count[t] = 0;
        if (self) {
            if (i == j0) jm++;
            if (i == jm) jm--;
            /* a block of one row has no other donor */
            if (jm > j1) break;
        }
        ssize_t j = jm;
        double d0, d1 = fabs(ps[i] - ps[j]);
        while (j >= j0 && j < j1) {
            j = psm_step(j, forward, me, noreplace, used, n_pos);
            if (j < j0 || j > j1) break;
            d0 = d1;
            d1 = fabs(ps[i] - ps[j]);
            if (d1 > d0) break;
            if (d1 < d0) jm = j;
        }
        uint32_t n = 0;
        double nearest = fabs(ps[i] - ps[jm]);
        if (nearest < cal) {
            used[jm] = 1;
            if (psm_lists_push(out, jm)) return -1;
            n++;
            /* Without replacement, ties are the unused donors that follow in
             * the scan direction at no greater distance (psmatch2's
             * match_ties); they are used up too. */
            if (ties && noreplace) {
                for (ssize_t k = psm_step(jm, forward, me, noreplace, used, n_pos);
                     k >= j0 && k <= j1 && fabs(ps[i] - ps[k]) <= nearest;
                     k = psm_step(k, forward, me, noreplace, used, n_pos)) {
                    if (psm_lists_push(out, k)) return -1;
                    n++;
                    used[k] = 1;
                }
            }
            ssize_t p0 = jm - 1, p1 = jm + 1;
            if (p0 == me) p0--;
            if (p1 == me) p1++;
            double e0 = 0.0, e1 = 0.0;
            /* The other neighbors, walking outward from the nearest.
             * With replacement: each donor through the kth-distance boundary
             * exactly once, even when the nearest run already has k.
             * Without replacement (psmatch2's match2; neighbors is above 1
             * only in the ai() pass): the next neighbors - 1 positions, used
             * or not, and none of them is marked as used. After ties a donor
             * can so be listed twice, as in psmatch2. */
            double boundary = nearest;
            uint32_t around = 0;
            for (; !noreplace || around + 1 < (uint32_t)neighbors;) {
                if (p0 >= j0) e0 = fabs(ps[i] - ps[p0]);
                if (p1 <= j1) e1 = fabs(ps[i] - ps[p1]);
                ssize_t pick;
                if (p0 >= j0 && (p1 > j1 || e0 <= e1)) {
                    pick = p0--;
                    if (p0 == me) p0--;
                } else if (p1 <= j1) {
                    pick = p1++;
                    if (p1 == me) p1++;
                } else {
                    break;
                }
                double distance = fabs(ps[i] - ps[pick]);
                if (!(distance < cal)) break;
                if (!noreplace && n >= (uint32_t)neighbors && (!ties || distance > boundary))
                    break;
                if (psm_lists_push(out, pick)) return -1;
                n++;
                around++;
                boundary = distance;
            }
        }
        out->count[t] = n;
        if (jm == j1 && i < i1) forward = 0;
        if (noreplace && n > 0) jm = psm_step(jm, forward, me, noreplace, used, n_pos);
    }
    /* the donors ran out: the remaining targets have no match */
    for (; i <= i1; i++) {
        out->start[i - i0] = out->len;
        out->count[i - i0] = 0;
    }
    return 0;
}

static size_t psm_lower_bound(const double *x, size_t n, double target)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (x[mid] < target) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/* The same pass with replacement and without ties, one target at a time:
 * lists have at most neighbors entries, at stride neighbors. */
static void psm_nearest(const double *ps, ssize_t i0, ssize_t i1, ssize_t j0, ssize_t j1,
                        int neighbors, double cal, psm_lists *out, int threads)
{
    const double *donor = ps + j0;
    const size_t n_donors = (size_t)(j1 - j0 + 1);
    const ssize_t n_targets = i1 - i0 + 1;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (ssize_t t = 0; t < n_targets; t++) {
        const double p = ps[i0 + t];
        perm_idx_t *list = out->item + (size_t)t * (size_t)neighbors;
        out->start[t] = (size_t)t * (size_t)neighbors;
        out->count[t] = 0;
        size_t at = psm_lower_bound(donor, n_donors, p), jm = at;
        if (at == n_donors || (at > 0 && fabs(p - donor[at - 1]) <= fabs(p - donor[at])))
            jm = psm_lower_bound(donor, n_donors, donor[at - 1]);
        if (!(fabs(p - donor[jm]) < cal)) continue;
        uint32_t n = 0;
        list[n++] = (perm_idx_t)(j0 + (ssize_t)jm);
        ssize_t p0 = (ssize_t)jm - 1;
        size_t p1 = jm + 1;
        while ((int)n < neighbors) {
            size_t pick;
            if (p0 >= 0 && (p1 >= n_donors || fabs(p - donor[p0]) <= fabs(p - donor[p1])))
                pick = (size_t)p0--;
            else if (p1 < n_donors)
                pick = p1++;
            else
                break;
            if (!(fabs(p - donor[pick]) < cal)) break;
            list[n++] = (perm_idx_t)(j0 + (ssize_t)pick);
        }
        out->count[t] = n;
    }
    out->len = (size_t)n_targets * (size_t)neighbors;
}

/* Outputs, by filtered row */
typedef struct
{
    double *weight, *support, *nn, *id, *ids, *pdif;
    double **matched; /* ny outcomes */
    double **self;    /* nself: the outcomes' same-arm values, then the score columns' */
    int neighbors, ny, nself;
} psm_outputs;

/* Record a matching pass between blocks (not ai()): psmatch2's _n#, _nn,
 * _weight, _pdif and matched outcomes, and support for unmatched targets.
 * Donor weights add up serially in the targets' order in position space
 * (wpos), as a serial pass adds them by row; the other outputs belong to
 * one target each. ypos[k] holds outcome k in position order. */
static void psm_apply(const psm_lists *l, ssize_t i0, ssize_t i1, const psm_order_t *order,
                      const double *ps, double *wpos, double *const *ypos, psm_outputs *o,
                      int threads)
{
    const ssize_t n_t = i1 - i0 + 1;
    for (ssize_t t = 0; t < n_t; t++) {
        const uint32_t n = l->count[t];
        const perm_idx_t *list = l->item + l->start[t];
        for (uint32_t m = 0; m < n; m++) wpos[list[m]] += 1.0 / (double)n;
    }
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (ssize_t t = 0; t < n_t; t++) {
        const size_t row = order[i0 + t].row;
        const uint32_t n = l->count[t];
        if (n == 0) {
            o->support[row] = 0.0;
            continue;
        }
        const perm_idx_t *list = l->item + l->start[t];
        o->nn[row] = (double)n;
        if (o->ids) {
            for (uint32_t m = 0; m < n && (int)m < o->neighbors; m++)
                o->ids[row * (size_t)o->neighbors + m] = (double)list[m] + 1.0;
        }
        if (o->pdif) o->pdif[row] = fabs(ps[i0 + t] - ps[list[0]]);
        for (int k = 0; k < o->ny; k++) {
            if (!o->matched || !o->matched[k]) continue;
            const double *y = ypos[k];
            double matched = 0.0;
            int missing = 0;
            for (uint32_t m = 0; m < n; m++) {
                const double v = y[list[m]];
                if (psm_missing(v, SV_missval)) missing = 1;
                else matched = matched + v / (double)n;
            }
            o->matched[k][row] = missing ? SV_missval : matched;
        }
    }
}

/* Record an ai() pass for one column (src in position order): the mean of
 * the same-arm matches, or with altvariance the conditional variance of
 * Abadie et al. (2004, p. 303) */
static void psm_apply_self(const psm_lists *l, ssize_t i0, ssize_t i1, const psm_order_t *order,
                           const double *src, int altvariance, double *self, int threads)
{
    const ssize_t n_t = i1 - i0 + 1;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (ssize_t t = 0; t < n_t; t++) {
        const ssize_t i = i0 + t;
        const size_t row = order[i].row;
        const uint32_t n = l->count[t];
        if (n == 0 || psm_missing(self[row], SV_missval)) continue;
        const perm_idx_t *list = l->item + l->start[t];
        double mean = self[row];
        const double own = src[i];
        int missing = psm_missing(own, SV_missval);
        for (uint32_t m = 0; m < n; m++) {
            const double v = src[list[m]];
            if (psm_missing(v, SV_missval)) missing = 1;
            else mean = mean + v / (double)n;
        }
        if (missing) {
            self[row] = SV_missval;
            continue;
        }
        if (altvariance) {
            const double center = ((double)n * mean + own) / ((double)n + 1.0);
            double sum = (own - center) * (own - center);
            for (uint32_t m = 0; m < n; m++) {
                const double d = src[list[m]] - center;
                sum = sum + d * d;
            }
            mean = sum / (double)n;
        }
        self[row] = mean;
    }
}

/* Kernel matching pass over every donor pair: targets are filtered rows;
 * donor_sorted (donor_rows by orig_idx) is ordered by score. A target
 * weighs the donors with |dif| <= bandwidth (the normal kernel: all, within
 * the window where its weight underflows to zero), normalized to sum to
 * one; the weights add into weight, and targets get _nn (positive weights),
 * support and the weighted mean of each outcome. Donor scores, outcomes and
 * weight sums are kept in sorted order, so a target reads one contiguous
 * run. */
static int psm_kernel_pass(const cpsmatch_options *opts, const double *ps,
                           const double *const *ycols, int ny, const size_t *targets,
                           size_t n_targets, const sorted_control_t *donor_sorted,
                           const size_t *donor_rows, size_t n_donors, double *weight, double *nn,
                           double *support, double *const *matched)
{
    const double missval = SV_missval;
    /* psmatch2 weighs by normalden(dif/bw), which is even in bw; the other
     * kernels need |dif| <= bw and so match nothing when bw <= 0 */
    double bw = opts->bandwidth;
    if (opts->kernel == KERNEL_NORMAL) bw = fabs(bw);
    if (n_targets == 0) return 0;
    if (n_donors == 0 || !(bw > 0.0)) {
        for (size_t t = 0; t < n_targets; t++) support[targets[t]] = 0.0;
        return 0;
    }
    if (opts->kernel != KERNEL_NORMAL && !opts->kernel_direct) {
        /* the polynomial kernels: moment sums over sorted windows (targets
         * come in score order) */
        const int threads = psm_threads(n_targets + n_donors);
        const size_t nyk = ny > 0 ? (size_t)ny : 1;
        double *tps = ctools_safe_malloc2(n_targets, sizeof(double));
        double *dps = ctools_safe_malloc2(n_donors, sizeof(double));
        double *dyv = ctools_safe_malloc2(n_donors * nyk, sizeof(double));
        double *tnn = ctools_safe_malloc2(n_targets, sizeof(double));
        int *on = ctools_safe_malloc2(n_targets, sizeof(int));
        double *tmatched = ctools_safe_malloc2(n_targets * nyk, sizeof(double));
        double *dweight = ctools_safe_malloc2(n_donors, sizeof(double));
        const double **dy = ctools_safe_calloc2(nyk, sizeof(double *));
        double **tm = ctools_safe_calloc2(nyk, sizeof(double *));
        int rc = 0;
        if (!tps || !dps || !dyv || !tnn || !on || !tmatched || !dweight || !dy || !tm) rc = 920;
        if (!rc) {
            for (size_t t = 0; t < n_targets; t++) tps[t] = ps[targets[t]];
            for (size_t j = 0; j < n_donors; j++) dps[j] = donor_sorted[j].pscore;
            for (int k = 0; k < ny; k++) {
                double *col = dyv + (size_t)k * n_donors;
                for (size_t j = 0; j < n_donors; j++)
                    col[j] = ycols[k][donor_rows[donor_sorted[j].orig_idx]];
                dy[k] = col;
                tm[k] = (matched && matched[k]) ? tmatched + (size_t)k * n_targets : NULL;
            }
            rc = cpsmatch_kernel_moments(opts->kernel, bw, tps, n_targets, dps, n_donors, dy, ny,
                                         threads, dweight, tnn, on, tm);
        }
        if (!rc) {
            for (size_t t = 0; t < n_targets; t++) {
                const size_t row = targets[t];
                if (!on[t]) {
                    support[row] = 0.0;
                    continue;
                }
                nn[row] = tnn[t];
                for (int k = 0; k < ny; k++)
                    if (tm[k]) matched[k][row] = tm[k][t];
            }
            for (size_t j = 0; j < n_donors; j++)
                weight[donor_rows[donor_sorted[j].orig_idx]] += dweight[j];
        }
        free(tps);
        free(dps);
        free(dyv);
        free(tnn);
        free(on);
        free(tmatched);
        free(dweight);
        free((void *)dy);
        free(tm);
        return rc;
    }
    /* exp(-u^2/2) underflows to exactly 0 beyond |u| = 38.6 */
    const double window = opts->kernel == KERNEL_NORMAL ? 38.7 * bw : bw;
    int n_threads = psm_threads(n_targets);
    const size_t nyk = ny > 0 ? (size_t)ny : 1;
    int failed = 0;
    double *dy = ctools_safe_malloc2(n_donors * nyk, sizeof(double));
    double *acc = ctools_safe_malloc2((size_t)n_threads * nyk * 2, sizeof(double));
    double **local = ctools_safe_calloc2((size_t)n_threads, sizeof(double *));
    /* a target's kernel weights, in its thread's workspace */
    match_workspace *ws = ctools_safe_calloc2((size_t)n_threads, sizeof(match_workspace));
    if (!dy || !acc || !local || !ws) failed = 1;
    for (int k = 0; k < n_threads && !failed; k++) {
        local[k] = ctools_safe_calloc2(n_donors, sizeof(double));
        cpsmatch_match_workspace_init(&ws[k], MATCH_WORKSPACE_INITIAL_CAP);
        if (!local[k] || !ws[k].match_weights) failed = 1;
    }
    if (!failed) {
        for (int k = 0; k < ny; k++)
            for (size_t j = 0; j < n_donors; j++)
                dy[(size_t)k * n_donors + j] = ycols[k][donor_rows[donor_sorted[j].orig_idx]];
#pragma omp parallel num_threads(n_threads) reduction(| : failed)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            double *sum_wy = acc + (size_t)tid * nyk * 2, *miss = sum_wy + nyk;
#pragma omp for schedule(dynamic, 64)
            for (size_t t = 0; t < n_targets; t++) {
                size_t row = targets[t], start, end;
                double p = ps[row];
                cpsmatch_distance_window(donor_sorted, n_donors, p, window, &start, &end);
                if (cpsmatch_match_workspace_ensure_match(&ws[tid], end - start)) {
                    failed = 1;
                    continue;
                }
                double *w = ws[tid].match_weights, total = 0.0;
                for (size_t j = start; j < end; j++) {
                    double dist = fabs(p - donor_sorted[j].pscore);
                    w[j - start] = (opts->kernel != KERNEL_NORMAL && !(dist <= bw))
                                       ? 0.0 : cpsmatch_apply_kernel(dist / bw, opts->kernel);
                    total += w[j - start];
                }
                if (!(total > 0.0)) {
                    support[row] = 0.0;
                    continue;
                }
                double sum_w = 0.0;
                int n = 0;
                for (int k = 0; k < ny; k++) sum_wy[k] = miss[k] = 0.0;
                for (size_t j = start; j < end; j++) {
                    if (!(w[j - start] > 0.0)) continue;
                    double wn = w[j - start] / total;
                    local[tid][j] += wn;
                    for (int k = 0; k < ny; k++) {
                        const double v = dy[(size_t)k * n_donors + j];
                        if (psm_missing(v, missval)) miss[k] = 1.0;
                        sum_wy[k] += wn * v;
                    }
                    sum_w += wn;
                    n++;
                }
                nn[row] = (double)n;
                for (int k = 0; k < ny; k++)
                    if (matched && matched[k])
                        matched[k][row] = miss[k] != 0.0 ? missval : sum_wy[k] / sum_w;
            }
        }
    }
    if (!failed) {
        for (size_t j = 0; j < n_donors; j++) {
            double sum = 0.0;
            for (int k = 0; k < n_threads; k++) sum += local[k][j];
            weight[donor_rows[donor_sorted[j].orig_idx]] += sum;
        }
    }
    for (int k = 0; local && k < n_threads; k++) free(local[k]);
    for (int k = 0; ws && k < n_threads; k++) cpsmatch_match_workspace_free(&ws[k]);
    free(dy);
    free(acc);
    free(local);
    free(ws);
    return failed ? 920 : 0;
}

/* Radius matching pass: each donor within the radius gets 1/count, through
 * the exact endpoint sums of cpsmatch_match_radius; the matched outcome is a
 * range mean over double-double prefix sums. */
static int psm_radius_pass(size_t n_rows, const double *ps, const double *const *ycols, int ny,
                           const size_t *targets, size_t n_targets,
                           const sorted_control_t *donor_sorted, const size_t *donor_rows,
                           size_t n_donors, const perm_idx_t *obs_map, double radius,
                           double *weight, double *nn, double *support, double *const *matched)
{
    double *target_ps = ctools_safe_malloc2(n_targets ? n_targets : 1, sizeof(double));
    double *scratch_w = ctools_safe_calloc2(n_rows ? n_rows : 1, sizeof(double));
    double *scratch_m = ctools_safe_malloc2(n_rows ? n_rows : 1, sizeof(double));
    dd_real *prefix = ny ? ctools_safe_calloc2(n_donors + 1, sizeof(dd_real)) : NULL;
    size_t *missing = ny ? ctools_safe_calloc2(n_donors + 1, sizeof(size_t)) : NULL;
    perm_idx_t *starts = ctools_safe_malloc2(n_targets ? n_targets : 1, sizeof(perm_idx_t));
    int rc = 0;
    if (!target_ps || !scratch_w || !scratch_m || !starts || (ny && (!prefix || !missing))) {
        rc = 920;
        goto done;
    }
    for (size_t t = 0; t < n_targets; t++) target_ps[t] = ps[targets[t]];
    int n_matched = 0, n_off = 0;
    rc = cpsmatch_match_radius(donor_sorted, n_donors, target_ps, targets, n_targets, donor_rows,
                               obs_map, radius, 0, 0.0, 0.0, scratch_w, scratch_m, nn, support,
                               &n_matched, &n_off, starts);
    if (rc) goto done;
    for (size_t j = 0; j < n_donors; j++) {
        size_t row = donor_rows[donor_sorted[j].orig_idx];
        weight[row] += scratch_w[row];
    }
    for (int k = 0; k < ny; k++) {
        const double *y = ycols[k];
        double *out = matched ? matched[k] : NULL;
        if (!out) continue;
        for (size_t j = 0; j < n_donors; j++) {
            double v = y[donor_rows[donor_sorted[j].orig_idx]];
            int absent = psm_missing(v, SV_missval);
            prefix[j + 1] = dd_add_d(prefix[j], absent ? 0.0 : v);
            missing[j + 1] = missing[j] + (size_t)absent;
        }
#pragma omp parallel for schedule(static) num_threads(psm_threads(n_targets))
        for (size_t t = 0; t < n_targets; t++) {
            size_t row = targets[t];
            if (nn[row] <= 0.0) continue;
            size_t start = starts[t], end = start + (size_t)nn[row];
            if (missing[end] != missing[start]) {
                out[row] = SV_missval;
                continue;
            }
            double sum = (prefix[end].hi - prefix[start].hi) + (prefix[end].lo - prefix[start].lo);
            double scale = fmax(fabs(prefix[end].hi), fabs(prefix[start].hi));
            if (!isfinite(sum) || fabs(sum) <= 1e-8 * scale) {
                /* Prefix subtraction can erase a small local result after
                 * enormous unrelated outcomes, even with double-double sums.
                 * Recompute only this interval; average before summing to
                 * avoid overflowing a finite local mean. */
                dd_real mean = {0.0, 0.0};
                for (size_t j = start; j < end; j++)
                    mean = dd_add_d(mean, y[donor_rows[donor_sorted[j].orig_idx]] / nn[row]);
                out[row] = mean.hi + mean.lo;
            } else {
                out[row] = sum / nn[row];
            }
        }
    }
done:
    free(target_ps);
    free(scratch_w);
    free(scratch_m);
    free(prefix);
    free(missing);
    free(starts);
    return rc;
}

/* Several columns written in one pass over Stata's records: workers own
 * contiguous ranges and write every column of a tile of rows before the
 * next tile. With fill_rows every observation 1..fill_rows is written, the
 * ones outside obs_map (increasing) as missing. */
typedef struct
{
    const double *col;
    int var;
} psm_store_col;

enum { PSM_STORE_TILE = 2048 };

static ST_retcode psm_store(const psm_store_col *cols, size_t ncols, size_t n,
                            const perm_idx_t *obs_map, size_t fill_rows)
{
    const size_t rows = fill_rows ? fill_rows : n;
    if (!ncols || !rows) return 0;
    const int threads = psm_threads(rows);
    ST_retcode first_rc = 0;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
        const size_t lo = rows * tid / team, hi = rows * (tid + 1) / team;
        size_t q = lo;
        if (fill_rows) {
            /* the first filtered row at or after observation lo + 1 */
            size_t a = 0, b = n;
            while (a < b) {
                size_t mid = a + (b - a) / 2;
                if ((size_t)obs_map[mid] < lo + 1) a = mid + 1;
                else b = mid;
            }
            q = a;
        }
        ST_retcode rc = 0;
        for (size_t tile = lo; tile < hi && !rc; tile += PSM_STORE_TILE) {
            const size_t end = hi - tile > PSM_STORE_TILE ? tile + PSM_STORE_TILE : hi;
            size_t next = q;
            for (size_t c = 0; c < ncols && !rc; c++) {
                const double *col = cols[c].col;
                const int var = cols[c].var;
                size_t p = q;
                if (fill_rows) {
                    for (size_t o = tile + 1; o <= end && !rc; o++) {
                        double v = SV_missval;
                        if (p < n && (size_t)obs_map[p] == o) v = col[p++];
                        rc = SF_vstore(var, (ST_int)o, v);
                    }
                } else {
                    for (; p < end && !rc; p++) rc = SF_vstore(var, (ST_int)obs_map[p], col[p]);
                }
                next = p;
            }
            q = next;
        }
        if (rc) {
#pragma omp critical(psm_store_rc)
            if (!first_rc) first_rc = rc;
        }
    }
    return first_rc;
}

/* Donors of a kernel pass, ordered by score (their block is) */
static void psm_kernel_donors(const double *ps, const perm_idx_t *obs_map, const size_t *rows,
                              size_t n, sorted_control_t *sorted)
{
    for (size_t k = 0; k < n; k++) {
        sorted[k].orig_idx = k;
        sorted[k].obs_idx = obs_map[rows[k]];
        sorted[k].pscore = ps[rows[k]];
    }
}

/* Stata's arithmetic on possibly missing values (the ai() helpers) */
static inline double psm_sub(double a, double b)
{
    return (psm_missing(a, SV_missval) || psm_missing(b, SV_missval)) ? SV_missval : a - b;
}
static inline double psm_add(double a, double b)
{
    return (psm_missing(a, SV_missval) || psm_missing(b, SV_missval)) ? SV_missval : a + b;
}
static inline double psm_mul(double a, double b)
{
    return (psm_missing(a, SV_missval) || psm_missing(b, SV_missval)) ? SV_missval : a * b;
}
static inline double psm_div(double a, double b)
{
    return (psm_missing(a, SV_missval) || psm_missing(b, SV_missval) || b == 0.0) ? SV_missval
                                                                                  : a / b;
}

/* ------------------------------------------------------------------------
 * The Abadie and Imbens (2016) correction terms, as psmatch2 4.1.0 (and
 * cpsmatch's former Mata function _cpsmatch_ai2016) computes them, over the
 * observations on support (the plugin call's if) in observation order: with
 * c a row vector over the score model's coefficients and V their variance,
 * each term is c V c': (ATE, ATT minus, ATT plus, ATU minus, ATU plus).
 * Sums run down the rows in order and skip missing terms, as Mata's colsum.
 * Variables: _pscore, the score's derivative, the outcome, _self_ and _
 * of the outcome, _treated, _id, _n1.._nK, the k model columns and their
 * same-arm means. Arguments: k K ai hascons n0 n1 att atu matrix keep...
 * ------------------------------------------------------------------------ */
static ST_retcode cpsmatch_ai2016_main(const char *args)
{
    int k = 0, K = 0, ai = 0, hascons = 0;
    double n0 = 0, n1 = 0, att = 0, atu = 0;
    char vname[64] = {0};
    const char *p = args;
    char *end;
    k = (int)strtol(p, &end, 10); p = end;
    K = (int)strtol(p, &end, 10); p = end;
    ai = (int)strtol(p, &end, 10); p = end;
    hascons = (int)strtol(p, &end, 10); p = end;
    double *targets[4] = {&n0, &n1, &att, &atu};
    for (int j = 0; j < 4; j++) {
        while (*p == ' ') p++;
        if (*p == '.' && (p[1] == ' ' || p[1] == '\0')) {
            *targets[j] = SV_missval;
            p++;
        } else {
            *targets[j] = strtod(p, &end);
            if (end == p) return 198;
            p = end;
        }
    }
    while (*p == ' ') p++;
    size_t len = strcspn(p, " ");
    if (len == 0 || len >= sizeof(vname) || k < 0 || K < 1 || ai < 1) return 198;
    memcpy(vname, p, len);
    p += len;
    const int kc = k + (hascons ? 1 : 0);
    int keep[128];
    if (kc > 128) return 198;
    for (int j = 0; j < kc; j++) {
        keep[j] = (int)strtol(p, &end, 10);
        if (end == p) return 198;
        p = end;
    }
    /* clear: then this outcome's _<outcome> is cleared for the controls;
     * final: and _n#, _nn and _weight as psmatch2 leaves them when it
     * matched the controls for the correction alone (without ate) */
    int clear = 0, final = 0;
    while (*p) {
        while (*p == ' ') p++;
        size_t w = strcspn(p, " ");
        if (w == 5 && !strncmp(p, "clear", 5)) clear = 1;
        else if (w == 5 && !strncmp(p, "final", 5)) final = 1;
        else if (w) return 198;
        p += w;
    }
    ST_retcode rc = 0;
    /* ... then _support, _weight, _nn; the sample marker last */
    const int nvars = 7 + K + 2 * k + 3;
    int *vars = ctools_safe_malloc2((size_t)nvars, sizeof(int));
    double *V = ctools_safe_malloc2((size_t)(kc ? kc * kc : 1), sizeof(double));
    if (!vars || !V) {
        free(vars);
        free(V);
        return 920;
    }
    for (int j = 0; j < nvars; j++) vars[j] = j + 1;
    for (int a = 0; a < kc && !rc; a++)
        for (int b = 0; b < kc && !rc; b++)
            rc = SF_mat_el(vname, keep[a], keep[b], &V[a * kc + b]);
    if (rc) {
        free(vars);
        free(V);
        return rc;
    }
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, vars, (size_t)nvars, 0, 0,
                                             CTOOLS_LOAD_NO_SORT_ORDER | CTOOLS_LOAD_MARKER);
    free(vars);
    if (load_rc != STATA_OK) {
        free(V);
        ctools_filtered_data_free(&filtered);
        return ctools_stata_rc(load_rc);
    }
    /* the rows on support, in observation order */
    const size_t n_all = filtered.data.nobs;
    const double *S_all = filtered.data.vars[7 + K + 2 * k].data.dbl;
    size_t *sel = ctools_safe_malloc2(n_all ? n_all : 1, sizeof(size_t));
    if (!sel) {
        free(V);
        ctools_filtered_data_free(&filtered);
        return 920;
    }
    size_t n = 0;
    for (size_t r = 0; r < n_all; r++)
        if (S_all[r] == 1.0) sel[n++] = r;
#define AI_COL(c, i) (filtered.data.vars[(c)].data.dbl[sel[(i)]])
    double q[5];
    for (int j = 0; j < 5; j++) q[j] = SV_missval;
    double *own = NULL, *other = NULL, *count = NULL, *c1 = NULL, *c0 = NULL;
    double *sums = NULL;
    perm_idx_t *rowof = NULL;
    if (n > 0) {
        const size_t nobs = (size_t)SF_nobs();
        const int threads = psm_threads(n);
        for (size_t i = 0; i < n && !rc; i++)
            if (!(AI_COL(0, i) > 0.0) || !(AI_COL(0, i) < 1.0)) {
                SF_error("cpsmatch: propensity score outside (0,1) in the ai() correction\n");
                rc = 3498;
            }
        if (!rc) {
            own = ctools_safe_malloc2(n * (size_t)(k ? k : 1), sizeof(double));
            other = ctools_safe_calloc2(n * (size_t)(k ? k : 1), sizeof(double));
            count = ctools_safe_calloc2(n, sizeof(double));
            c1 = ctools_safe_malloc2(n * (size_t)(k ? k : 1), sizeof(double));
            c0 = ctools_safe_malloc2(n * (size_t)(k ? k : 1), sizeof(double));
            sums = ctools_safe_calloc2((size_t)kc * 5, sizeof(double));
            rowof = ctools_safe_malloc2(nobs + 1, sizeof(perm_idx_t));
            if (!own || !other || !count || !c1 || !c0 || !sums || !rowof) rc = 920;
        }
        if (!rc) {
            const double frac = (double)ai / ((double)ai + 1.0);
            for (int j = 0; j < k; j++) {
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
                for (size_t i = 0; i < n; i++)
                    own[(size_t)j * n + i] =
                        psm_mul(psm_mul(frac, psm_sub(AI_COL(7 + K + j, i), AI_COL(7 + K + k + j, i))),
                                psm_sub(AI_COL(2, i), AI_COL(3, i)));
            }
            /* rowof[id] = 1::n (the ids are distinct) */
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t o = 0; o <= nobs; o++) rowof[o] = PERM_IDX_MAX;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t i = 0; i < n; i++) {
                const double id = AI_COL(6, i);
                if (!psm_missing(id, SV_missval) && id >= 1 && id <= (double)nobs)
                    rowof[(size_t)id] = (perm_idx_t)i;
            }
            /* each row's sums over its matches, in the matches' order */
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t i = 0; i < n; i++) {
                for (int m = 0; m < K; m++) {
                    const double donor = AI_COL(7 + m, i);
                    if (psm_missing(donor, SV_missval)) continue;
                    count[i] += 1.0;
                    size_t d = (size_t)donor;
                    if (donor < 1 || d > nobs || rowof[d] == PERM_IDX_MAX) continue;
                    const size_t r = rowof[d];
                    for (int j = 0; j < k; j++)
                        other[(size_t)j * n + i] =
                            psm_add(other[(size_t)j * n + i], own[(size_t)j * n + r]);
                }
            }
            for (size_t i = 0; i < n && !rc; i++)
                if (count[i] == 0.0) {
                    SF_error("cpsmatch: an observation on support has no match in the ai() correction\n");
                    rc = 3498;
                }
        }
        if (!rc) {
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t i = 0; i < n; i++) {
                const double t = AI_COL(5, i) == 1.0 ? 1.0 : 0.0;
                for (int j = 0; j < k; j++) {
                    const size_t at = (size_t)j * n + i;
                    other[at] = psm_div(other[at], count[i]);
                    c1[at] = psm_add(other[at], psm_mul(t, psm_sub(own[at], other[at])));
                    c0[at] = psm_add(own[at], psm_mul(t, psm_sub(other[at], own[at])));
                }
            }
            /* sums[e * kc + j]: column j of cA, cT, dT, cU, dU before
             * scaling, each summed down the rows in order (Mata's colsum);
             * the 5 * kc sums are independent */
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if (threads > 1)
            for (int task = 0; task < 5 * kc; task++) {
                const int e = task / kc, j = task % kc;
                double sum = 0.0;
                for (size_t i = 0; i < n; i++) {
                    const double t = AI_COL(5, i) == 1.0 ? 1.0 : 0.0;
                    const double pi = AI_COL(0, i), qi = 1.0 - pi;
                    double x, a1, a0;
                    if (j < k) {
                        x = AI_COL(7 + K + j, i);
                        a1 = c1[(size_t)j * n + i];
                        a0 = c0[(size_t)j * n + i];
                    } else {
                        x = 1.0;
                        a1 = a0 = 0.0;
                    }
                    double term;
                    if (e == 0) {
                        term = psm_add(psm_div(a1, pi), psm_div(a0, qi));
                    } else {
                        const double delta =
                            psm_mul(psm_sub(2.0 * t, 1.0), psm_sub(AI_COL(3, i), AI_COL(4, i)));
                        const double aT = psm_sub(delta, att), aU = psm_sub(delta, atu);
                        /* Mata binds - before + (and / before *):
                         * a + b - c is a + (b - c) */
                        if (e == 1)
                            term = psm_add(psm_add(psm_mul(x, aT), a1), psm_mul(psm_div(pi, qi), a0));
                        else if (e == 2)
                            term = psm_add(psm_mul(x, aT), psm_sub(a1, a0));
                        else if (e == 3)
                            term = psm_sub(psm_sub(psm_mul(x, aU), a0), psm_mul(psm_div(qi, pi), a1));
                        else
                            term = psm_add(psm_mul(x, aU), psm_sub(a1, a0));
                    }
                    const double v = psm_mul(AI_COL(1, i), term);
                    if (!psm_missing(v, SV_missval)) sum = sum + v;
                }
                sums[e * kc + j] = sum;
            }
            const double div[5] = {n0 + n1, n1, n1, n0, n0};
            for (int e = 0; e < 5; e++) {
                double c[128], r[128];
                for (int j = 0; j < kc; j++) c[j] = psm_div(sums[e * kc + j], div[e]);
                /* (c * V) * c' */
                for (int b = 0; b < kc; b++) {
                    double s = 0.0;
                    for (int a = 0; a < kc; a++) s = psm_add(s, psm_mul(c[a], V[a * kc + b]));
                    r[b] = s;
                }
                double s = 0.0;
                for (int b = 0; b < kc; b++) s = psm_add(s, psm_mul(r[b], c[b]));
                q[e] = s;
            }
        }
    }
#undef AI_COL
    free(own);
    free(other);
    free(count);
    free(c1);
    free(c0);
    free(sums);
    free(rowof);
    free(V);
    free(sel);
    /* psmatch2's clean-up of the controls' results (see the arguments) */
    if (!rc && (clear || final) && n_all) {
        const double *tr = filtered.data.vars[5].data.dbl;
        double *buffer = ctools_safe_malloc2(n_all * (size_t)(final ? K + 3 : 1), sizeof(double));
        psm_store_col *cols = ctools_safe_malloc2((size_t)K + 3, sizeof(psm_store_col));
        if (!buffer || !cols) rc = 920;
        size_t nc = 0;
        double *at = buffer;
        if (!rc && clear) {
            const double *my = filtered.data.vars[4].data.dbl;
            for (size_t r = 0; r < n_all; r++) at[r] = tr[r] == 0.0 ? SV_missval : my[r];
            cols[nc++] = (psm_store_col){at, 5};
            at += n_all;
        }
        if (!rc && final) {
            for (int m = 0; m < K; m++) {
                const double *v = filtered.data.vars[7 + m].data.dbl;
                for (size_t r = 0; r < n_all; r++) at[r] = tr[r] == 0.0 ? SV_missval : v[r];
                cols[nc++] = (psm_store_col){at, 8 + m};
                at += n_all;
            }
            const double *S = filtered.data.vars[7 + K + 2 * k].data.dbl;
            const double *W = filtered.data.vars[8 + K + 2 * k].data.dbl;
            const double *NN = filtered.data.vars[9 + K + 2 * k].data.dbl;
            for (size_t r = 0; r < n_all; r++) {
                double w = (tr[r] == 1.0 && S[r] == 1.0) ? 1.0 : W[r];
                at[r] = (S[r] == 0.0 || w == 0.0) ? SV_missval : w;
            }
            cols[nc++] = (psm_store_col){at, 9 + K + 2 * k};
            at += n_all;
            for (size_t r = 0; r < n_all; r++) at[r] = tr[r] == 0.0 ? SV_missval : NN[r];
            cols[nc++] = (psm_store_col){at, 10 + K + 2 * k};
        }
        if (!rc) rc = psm_store(cols, nc, n_all, filtered.obs_map, 0);
        free(buffer);
        free(cols);
    }
    ctools_filtered_data_free(&filtered);
    if (rc) return rc;
    static const char *names[5] = {"_cpsmatch_psq1", "_cpsmatch_psq2", "_cpsmatch_psq3",
                                   "_cpsmatch_psq4", "_cpsmatch_psq5"};
    for (int e = 0; e < 5 && !rc; e++) rc = SF_scal_save((char *)names[e], q[e]);
    return rc;
}

enum { PSM_FIT_MAXK = 64 };

/* ------------------------------------------------------------------------
 * The Mata engine's matching loops (Mahalanobis metric, and local linear
 * regression with kernels other than the Epanechnikov one), over the data
 * as _cpsmatch_mata sorted it: a position is an observation. Arithmetic
 * follows the Mata code it replaces, operation for operation: the
 * quadratic form row by row, (delta * W) :* delta summed across, Mata's
 * quadsum as an exact sum rounded once, missing values propagating.
 * ------------------------------------------------------------------------ */
typedef struct
{
    double d;
    double pos; /* 1-based position */
} psm_dpair;

/* psmatch2's _select: the k smallest distances by Hoare's selection, in the
 * order it leaves them (a[0..n), k >= 1) */
static void psm_hoare_select(psm_dpair *a, size_t n, size_t k)
{
    /* 1-based indices, as in the Mata code */
    ptrdiff_t lo = 1, hi = (ptrdiff_t)n, kk = (ptrdiff_t)k;
    while (lo < hi) {
        const double pivot = a[kk - 1].d;
        ptrdiff_t i = lo, j = hi;
        do {
            while (a[i - 1].d < pivot) i++;
            while (pivot < a[j - 1].d) j--;
            if (i <= j) {
                psm_dpair s = a[i - 1];
                a[i - 1] = a[j - 1];
                a[j - 1] = s;
                i++;
                j--;
            }
        } while (i <= j);
        if (j < kk) lo = i;
        if (kk < i) hi = j;
    }
}

/* Mata's quadratic form of the row difference: rowsum((delta * W) :* delta) */
static inline double psm_qform(const double *const *X, int nx, size_t r, size_t i, const double *W)
{
    double delta[PSM_FIT_MAXK], s = 0.0;
    for (int l = 0; l < nx; l++) delta[l] = psm_sub(X[l][r], X[l][i]);
    for (int m = 0; m < nx; m++) {
        double t = 0.0;
        for (int l = 0; l < nx; l++) t = psm_add(t, psm_mul(delta[l], W[l * nx + m]));
        const double v = psm_mul(t, delta[m]);
        if (!psm_missing(v, SV_missval)) s = s + v; /* rowsum skips missing */
    }
    return s;
}

static int psm_parse_dot(const char **p, double *v)
{
    char *end;
    while (**p == ' ') (*p)++;
    if (**p == '.' && ((*p)[1] == ' ' || (*p)[1] == '\0')) {
        *v = SV_missval;
        (*p)++;
        return 0;
    }
    *v = strtod(*p, &end);
    if (end == *p) return -1;
    *p = end;
    return 0;
}

/* _cpsmatch_nn. Arguments: t0 t1 c0 c1 neighbors caliper W self ate altvar
 * nx nout. Variables: the nx metric columns, _weight, _support, [_nn
 * _n1.._nK] (not self), the nout outcomes and their matched columns. */
static ST_retcode cpsmatch_mnn_main(const char *args)
{
    double v[6];
    const char *p = args;
    for (int j = 0; j < 6; j++)
        if (psm_parse_dot(&p, &v[j])) return 198;
    char wname[64];
    int self = 0, ate = 0, altvar = 0, nx = 0, nout = 0;
    if (sscanf(p, " %63s %d %d %d %d %d", wname, &self, &ate, &altvar, &nx, &nout) != 6 ||
        nx < 1 || nx > PSM_FIT_MAXK || nout < 0)
        return 198;
    const ptrdiff_t t0 = (ptrdiff_t)v[0], t1 = (ptrdiff_t)v[1], c0 = (ptrdiff_t)v[2],
                    c1 = (ptrdiff_t)v[3];
    const int K = (int)v[4];
    const double caliper = v[5];
    if (K < 1) return 198;
    const int nids = self ? 0 : K + 1; /* _nn and _n1.._nK */
    const int nvars = nx + 2 + nids + 2 * nout;
    ST_retcode rc = 0;
    double W[PSM_FIT_MAXK * PSM_FIT_MAXK];
    for (int a = 0; a < nx && !rc; a++)
        for (int b = 0; b < nx && !rc; b++) rc = SF_mat_el(wname, a + 1, b + 1, &W[a * nx + b]);
    if (rc) return rc;
    int *vars = ctools_safe_malloc2((size_t)nvars, sizeof(int));
    if (!vars) return 920;
    for (int j = 0; j < nvars; j++) vars[j] = j + 1;
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    stata_retcode lrc =
        ctools_data_load(&fd, vars, (size_t)nvars, 0, 0, CTOOLS_LOAD_NO_SORT_ORDER | CTOOLS_LOAD_CHECK_IF);
    free(vars);
    if (lrc != STATA_OK) {
        ctools_filtered_data_free(&fd);
        return ctools_stata_rc(lrc);
    }
    const size_t n = fd.data.nobs;
    const double *X[PSM_FIT_MAXK];
    for (int l = 0; l < nx; l++) X[l] = fd.data.vars[l].data.dbl;
    double *weight = fd.data.vars[nx].data.dbl, *support = fd.data.vars[nx + 1].data.dbl;
    double *nn = self ? NULL : fd.data.vars[nx + 2].data.dbl;
    double *ids0 = self ? NULL : fd.data.vars[nx + 3].data.dbl;
    (void)ids0;
    const int ybase = nx + 2 + nids;
    if (t1 < t0 || c1 < c0 || t0 < 1 || (size_t)t1 > n || (c1 >= c0 && (c0 < 1 || (size_t)c1 > n))) {
        if (!self && t1 >= t0 && t0 >= 1 && (size_t)t1 <= n) {
            for (ptrdiff_t i = t0; i <= t1; i++) support[i - 1] = 0.0;
            psm_store_col col = {support, nx + 2};
            rc = psm_store(&col, 1, n, fd.obs_map, 0);
        }
        ctools_filtered_data_free(&fd);
        return rc;
    }
    const size_t n_t = (size_t)(t1 - t0 + 1), n_d = (size_t)(c1 - c0 + 1);
    const int threads = psm_threads(n_t * n_d / 64 + 1);
    psm_dpair *work = ctools_safe_malloc2((size_t)threads * n_d, sizeof(psm_dpair));
    double *picks = ctools_safe_malloc2(n_t * (size_t)K, sizeof(double));
    int *npick = ctools_safe_calloc2(n_t, sizeof(int));
    unsigned *uses = self ? NULL : ctools_safe_calloc2(n, sizeof(unsigned));
    if (!work || !picks || !npick || (!self && !uses)) {
        free(work);
        free(picks);
        free(npick);
        free(uses);
        ctools_filtered_data_free(&fd);
        return 920;
    }
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads) if (threads > 1)
    for (ptrdiff_t i = t0; i <= t1; i++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        const size_t t = (size_t)(i - t0);
        npick[t] = -1; /* skipped */
        if (self && !ate && weight[i - 1] == 0.0) continue;
        psm_dpair *a = work + (size_t)tid * n_d;
        size_t m = 0;
        for (ptrdiff_t r = c0; r <= c1; r++) {
            if (self && r == i) continue;
            a[m].d = psm_qform(X, nx, (size_t)(r - 1), (size_t)(i - 1), W);
            a[m].pos = (double)r;
            m++;
        }
        const size_t k = (size_t)K < m ? (size_t)K : m;
        int count = 0;
        if (k >= 1) {
            psm_hoare_select(a, m, k);
            for (size_t q = 0; q < k; q++)
                if (a[q].d <= caliper || (psm_missing(caliper, SV_missval) && !psm_missing(a[q].d, SV_missval)))
                    picks[t * (size_t)K + (size_t)count++] = a[q].pos;
        }
        npick[t] = count;
    }
    /* outputs, target by target (the targets' rows are their own) */
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (ptrdiff_t i = t0; i <= t1; i++) {
        const size_t t = (size_t)(i - t0), row = (size_t)(i - 1);
        const int cnt = npick[t];
        if (cnt < 0) continue;
        const double *pk = picks + t * (size_t)K;
        if (!self) {
            nn[row] = (double)cnt;
            for (int q = 0; q < cnt; q++) fd.data.vars[nx + 3 + q].data.dbl[row] = pk[q];
        }
        for (int o = 0; o < nout; o++) {
            const double *Y = fd.data.vars[ybase + o].data.dbl;
            double *M = fd.data.vars[ybase + nout + o].data.dbl;
            double mi = M[row];
            for (int q = 0; q < cnt; q++) mi = psm_add(mi, psm_div(Y[(size_t)pk[q] - 1], (double)cnt));
            if (altvar && self) {
                const double mean = psm_div(psm_add(psm_mul((double)cnt, mi), Y[row]), (double)cnt + 1.0);
                double d0 = psm_sub(Y[row], mean), acc = psm_mul(d0, d0);
                for (int q = 0; q < cnt; q++) {
                    const double dq = psm_sub(Y[(size_t)pk[q] - 1], mean);
                    acc = psm_add(acc, psm_mul(dq, dq));
                }
                mi = psm_div(acc, (double)cnt);
            }
            M[row] = mi;
        }
        if (cnt < 1 && !self) support[row] = 0.0;
    }
    /* weight[pick] += 1 / neighbors, as many times as a donor is picked
     * (the same addend, so the order of the additions does not matter) */
    if (!self) {
        for (size_t t = 0; t < n_t; t++)
            for (int q = 0; q < npick[t]; q++) uses[(size_t)picks[t * (size_t)K + (size_t)q] - 1]++;
        const double add = 1.0 / (double)K;
        for (size_t r = 0; r < n; r++)
            for (unsigned u = 0; u < uses[r]; u++) weight[r] = psm_add(weight[r], add);
    }
    psm_store_col cols[4 + 2 * PSM_FIT_MAXK];
    size_t nc = 0;
    cols[nc++] = (psm_store_col){weight, nx + 1};
    cols[nc++] = (psm_store_col){support, nx + 2};
    for (int q = 0; q < nids; q++) cols[nc++] = (psm_store_col){fd.data.vars[nx + 2 + q].data.dbl, nx + 3 + q};
    for (int o = 0; o < nout && nc < sizeof(cols) / sizeof(cols[0]); o++)
        cols[nc++] = (psm_store_col){fd.data.vars[ybase + nout + o].data.dbl, ybase + nout + o + 1};
    rc = psm_store(cols, nc, n, fd.obs_map, 0);
    free(work);
    free(picks);
    free(npick);
    free(uses);
    ctools_filtered_data_free(&fd);
    return rc;
}

/* _cpsmatch_kernel's weight before normalization (missing outside the
 * kernel's support, normal kernel aside) */
static inline double psm_mata_kernel(double d, kernel_type_t kernel, double bw)
{
    const double u = psm_div(d, bw);
    const int inside = !psm_missing(d, SV_missval) && fabs(d) <= bw;
    if (kernel == KERNEL_NORMAL)
        return psm_missing(u, SV_missval) ? SV_missval : 0.39894228040143267794 * exp(-0.5 * u * u);
    if (psm_missing(u, SV_missval) || !inside) return SV_missval; /* w :/ 0 */
    double w;
    if (kernel == KERNEL_EPAN) w = 1.0 - u * u;
    else if (kernel == KERNEL_BIWEIGHT) {
        const double t = 1.0 - u * u;
        w = t * t;
    } else if (kernel == KERNEL_UNIFORM) w = 1.0;
    else {
        const double a = fabs(u), t = 1.0 - a * a * a;
        w = t * t * t;
    }
    return w;
}

/* _cpsmatch_smooth. Arguments: c0 c1 t0 t1 llr mahal W kernel bw ate nx
 * nout. Variables: the nx metric columns, _weight, _support, _treated,
 * the nout outcomes and their matched columns. */
static ST_retcode cpsmatch_msmooth_main(const char *args)
{
    double v[4], bw;
    const char *p = args;
    for (int j = 0; j < 4; j++)
        if (psm_parse_dot(&p, &v[j])) return 198;
    int llr = 0, mahal = 0, kernel = 0, ate = 0, nx = 0, nout = 0;
    char wname[64];
    if (sscanf(p, " %d %d %63s %d %lf %d %d %d", &llr, &mahal, wname, &kernel, &bw, &ate, &nx,
               &nout) != 8 || nx < 1 || nx > PSM_FIT_MAXK || nout < 0)
        return 198;
    const ptrdiff_t c0 = (ptrdiff_t)v[0], c1 = (ptrdiff_t)v[1], t0 = (ptrdiff_t)v[2],
                    t1 = (ptrdiff_t)v[3];
    ST_retcode rc = 0;
    double W[PSM_FIT_MAXK * PSM_FIT_MAXK];
    if (mahal)
        for (int a = 0; a < nx && !rc; a++)
            for (int b = 0; b < nx && !rc; b++) rc = SF_mat_el(wname, a + 1, b + 1, &W[a * nx + b]);
    if (rc) return rc;
    const int nvars = nx + 3 + 2 * nout;
    int *vars = ctools_safe_malloc2((size_t)nvars, sizeof(int));
    if (!vars) return 920;
    for (int j = 0; j < nvars; j++) vars[j] = j + 1;
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    stata_retcode lrc =
        ctools_data_load(&fd, vars, (size_t)nvars, 0, 0, CTOOLS_LOAD_NO_SORT_ORDER | CTOOLS_LOAD_CHECK_IF);
    free(vars);
    if (lrc != STATA_OK) {
        ctools_filtered_data_free(&fd);
        return ctools_stata_rc(lrc);
    }
    const size_t n = fd.data.nobs;
    const double *X[PSM_FIT_MAXK];
    for (int l = 0; l < nx; l++) X[l] = fd.data.vars[l].data.dbl;
    double *weight = fd.data.vars[nx].data.dbl, *support = fd.data.vars[nx + 1].data.dbl;
    const double *treated = fd.data.vars[nx + 2].data.dbl;
    const int ybase = nx + 3;
    size_t *donors = ctools_safe_malloc2(n ? n : 1, sizeof(size_t));
    double *Q = ctools_safe_malloc2(n ? n : 1, sizeof(double));  /* per target: quadsum(w) */
    double *B = ctools_safe_malloc2(n ? n : 1, sizeof(double));  /* llr: b */
    double *S = ctools_safe_malloc2(n ? n : 1, sizeof(double));  /* llr: s */
    unsigned char *has = ctools_safe_calloc2(n ? n : 1, 1);      /* the target has weights */
    if (!donors || !Q || !B || !S || !has) {
        free(donors); free(Q); free(B); free(S); free(has);
        ctools_filtered_data_free(&fd);
        return 920;
    }
    const kernel_type_t kt = (kernel_type_t)kernel;
#define MS_D(j, i) (mahal ? psm_qform(X, nx, (j), (i), W) : psm_sub(X[0][(j)], X[0][(i)]))
    for (int phase = ate ? 0 : 1; phase < 2 && !rc; phase++) {
        const ptrdiff_t a0 = phase ? t0 : c0, a1 = phase ? t1 : c1;
        if (a1 < a0) continue;
        const double arm = phase ? 0.0 : 1.0; /* the donors' arm */
        size_t nd = 0;
        for (size_t r = 0; r < n; r++)
            if (support[r] == 1.0 && treated[r] == arm) donors[nd++] = r;
        const size_t n_t = (size_t)(a1 - a0 + 1);
        const int threads = psm_threads(n_t * (nd ? nd : 1) / 64 + 1);
        if (nd == 0) {
            for (ptrdiff_t i = a0; i <= a1; i++)
                if (!llr || nout) support[i - 1] = 0.0;
            continue;
        }
        /* per target: normalizer, llr moments, outcomes */
#pragma omp parallel for schedule(dynamic, 8) num_threads(threads) if (threads > 1)
        for (ptrdiff_t i = a0; i <= a1; i++) {
            const size_t row = (size_t)(i - 1);
            cpsmatch_xsum q = {0}, sw = {0}, ss = {0}, sb = {0};
            size_t ok = 0, pos = 0;
            for (size_t j = 0; j < nd; j++) {
                const double w = psm_mata_kernel(MS_D(donors[j], row), kt, bw);
                if (!psm_missing(w, SV_missval)) cpsmatch_xsum_add(&q, w);
            }
            const double qsum = cpsmatch_xsum_value(&q);
            Q[row] = qsum;
            for (size_t j = 0; j < nd; j++) {
                const double d = MS_D(donors[j], row);
                const double w = psm_div(psm_mata_kernel(d, kt, bw), qsum);
                if (psm_missing(w, SV_missval)) continue;
                ok++;
                if (llr && w > 0.0) {
                    pos++;
                    cpsmatch_xsum_add(&sw, w);
                    const double wd = psm_mul(w, d);
                    if (!psm_missing(wd, SV_missval)) cpsmatch_xsum_add(&ss, wd);
                }
            }
            has[row] = ok > 0;
            if (!llr) {
                if (!ok) {
                    support[row] = 0.0;
                    continue;
                }
                for (int o = 0; o < nout; o++) {
                    const double *Y = fd.data.vars[ybase + o].data.dbl;
                    double *M = fd.data.vars[ybase + nout + o].data.dbl;
                    cpsmatch_xsum a = {0}, b = {0};
                    size_t cnt = 0;
                    for (size_t j = 0; j < nd; j++) {
                        const double y = Y[donors[j]];
                        if (psm_missing(y, SV_missval)) continue;
                        const double w = psm_div(psm_mata_kernel(MS_D(donors[j], row), kt, bw), qsum);
                        if (psm_missing(w, SV_missval)) continue;
                        const double wy = psm_mul(w, y);
                        if (!psm_missing(wy, SV_missval)) cpsmatch_xsum_add(&a, wy);
                        cpsmatch_xsum_add(&b, w);
                        cnt++;
                    }
                    if (cnt) {
                        const double mean = psm_div(cpsmatch_xsum_value(&a), cpsmatch_xsum_value(&b));
                        if (!psm_missing(mean, SV_missval)) M[row] = mean;
                    }
                }
                continue;
            }
            /* llr: the local linear weights and, per outcome, the intercept */
            const double swv = cpsmatch_xsum_value(&sw), s = cpsmatch_xsum_value(&ss);
            const double mean = psm_div(s, swv);
            double bval = SV_missval;
            if (pos > 1) {
                for (size_t j = 0; j < nd; j++) {
                    const double d = MS_D(donors[j], row);
                    const double w = psm_div(psm_mata_kernel(d, kt, bw), qsum);
                    if (psm_missing(w, SV_missval) || !(w > 0.0)) continue;
                    const double dm = psm_sub(d, mean), t = psm_mul(w, psm_mul(dm, dm));
                    if (!psm_missing(t, SV_missval)) cpsmatch_xsum_add(&sb, t);
                }
                bval = psm_div(cpsmatch_xsum_value(&sb), swv);
            }
            B[row] = bval;
            S[row] = s;
            for (int o = 0; o < nout; o++) {
                const double *Y = fd.data.vars[ybase + o].data.dbl;
                double *M = fd.data.vars[ybase + nout + o].data.dbl;
                cpsmatch_xsum a = {0}, b = {0}, c = {0};
                size_t cnt = 0;
                for (size_t j = 0; j < nd; j++) {
                    const double y = Y[donors[j]];
                    const double d = MS_D(donors[j], row);
                    const double w = psm_div(psm_mata_kernel(d, kt, bw), qsum);
                    if (psm_missing(y, SV_missval) || psm_missing(w, SV_missval) || !(w > 0.0)) continue;
                    cnt++;
                    cpsmatch_xsum_add(&a, w);
                    const double wd = psm_mul(w, d), wy = psm_mul(w, y);
                    if (!psm_missing(wd, SV_missval)) cpsmatch_xsum_add(&b, wd);
                    if (!psm_missing(wy, SV_missval)) cpsmatch_xsum_add(&c, wy);
                }
                if (cnt < 2) {
                    support[row] = 0.0;
                    continue;
                }
                const double swo = cpsmatch_xsum_value(&a);
                const double xbar = psm_div(cpsmatch_xsum_value(&b), swo);
                const double ybar = psm_div(cpsmatch_xsum_value(&c), swo);
                cpsmatch_xsum sxx = {0}, sxy = {0};
                for (size_t j = 0; j < nd; j++) {
                    const double y = Y[donors[j]];
                    const double d = MS_D(donors[j], row);
                    const double w = psm_div(psm_mata_kernel(d, kt, bw), qsum);
                    if (psm_missing(y, SV_missval) || psm_missing(w, SV_missval) || !(w > 0.0)) continue;
                    const double dx = psm_sub(d, xbar), dy = psm_sub(y, ybar);
                    const double t1v = psm_mul(w, psm_mul(dx, dx));
                    const double t2v = psm_mul(psm_mul(w, dx), dy);
                    if (!psm_missing(t1v, SV_missval)) cpsmatch_xsum_add(&sxx, t1v);
                    if (!psm_missing(t2v, SV_missval)) cpsmatch_xsum_add(&sxy, t2v);
                }
                const double vxx = cpsmatch_xsum_value(&sxx), vxy = cpsmatch_xsum_value(&sxy);
                M[row] = vxx > 0.0 ? psm_sub(ybar, psm_mul(psm_div(vxy, vxx), xbar)) : ybar;
            }
        }
        /* donor weights: each donor's additions in the targets' order */
#pragma omp parallel for schedule(dynamic, 8) num_threads(threads) if (threads > 1)
        for (size_t j = 0; j < nd; j++) {
            const size_t dj = donors[j];
            double wj = weight[dj];
            for (ptrdiff_t i = a0; i <= a1; i++) {
                const size_t row = (size_t)(i - 1);
                if (!has[row]) continue;
                const double d = MS_D(dj, row);
                const double w = psm_div(psm_mata_kernel(d, kt, bw), Q[row]);
                if (psm_missing(w, SV_missval)) continue;
                if (!llr) wj = psm_add(wj, w);
                else {
                    /* Mata binds / before * and - before +: w * (b +
                     * s^2 - s*d) / b is w * ((b + (s^2 - s*d)) / b) */
                    const double b = B[row], s = S[row];
                    wj = psm_add(wj, psm_mul(w, psm_div(psm_add(b, psm_sub(psm_mul(s, s), psm_mul(s, d))), b)));
                }
            }
            weight[dj] = wj;
        }
    }
#undef MS_D
    psm_store_col cols[3 + PSM_FIT_MAXK];
    size_t nc = 0;
    cols[nc++] = (psm_store_col){weight, nx + 1};
    cols[nc++] = (psm_store_col){support, nx + 2};
    for (int o = 0; o < nout && nc < sizeof(cols) / sizeof(cols[0]); o++)
        cols[nc++] = (psm_store_col){fd.data.vars[ybase + nout + o].data.dbl, ybase + nout + o + 1};
    rc = psm_store(cols, nc, n, fd.obs_map, 0);
    free(donors);
    free(Q);
    free(B);
    free(S);
    free(has);
    ctools_filtered_data_free(&fd);
    return rc;
}

/* Rows outside the sample, for _id after the sample's rows: by score (its
 * negative with descending; missing scores last), then by observation */
typedef struct
{
    perm_idx_t *obs; /* outside observations, increasing */
    double *score;   /* their _pscore */
    double *id;      /* their _id, in obs order */
    size_t n;
} psm_outside;

static int psm_outside_build(psm_outside *o, const perm_idx_t *obs_map, size_t n_rows,
                             size_t nobs, int pscore_var, int descending, size_t first,
                             psm_radix_t *keys, psm_radix_t *spare, size_t *hist)
{
    size_t q = 0, k = 0;
    for (size_t obs = 1; obs <= nobs; obs++) {
        if (q < n_rows && (size_t)obs_map[q] == obs) {
            q++;
            continue;
        }
        o->obs[k++] = (perm_idx_t)obs;
    }
    o->n = k;
    for (size_t j = 0; j < k; j++) {
        double v = SV_missval;
        int rc = SF_vdata(pscore_var, (ST_int)o->obs[j], &v);
        if (rc) return rc;
        o->score[j] = v;
        keys[j].key = psm_radix_key(descending && !psm_missing(v, SV_missval) ? -v : v);
        keys[j].at = (uint32_t)j;
    }
    psm_radix_t *sorted = psm_radix_sort(keys, spare, k, hist, psm_threads(k));
    for (size_t r = 0; r < k; r++) o->id[sorted[r].at] = (double)(first + r + 1);
    return 0;
}

ST_retcode cpsmatch_main(const char *args)
{
    cpsmatch_options opts;
    ST_retcode rc = 0;
    double t_start = ctools_timer_seconds(), t_load, t_setup, t_sort = 0.0, t_match = 0.0;
    double t_store = 0.0, t_stats = 0.0;

    if (args == NULL || *args == '\0') {
        SF_error("cpsmatch: no arguments specified\n");
        return 198;
    }
    {
        const char *a = args;
        while (*a == ' ') a++;
        if (!strncmp(a, "ai2016", 6) && (a[6] == ' ' || a[6] == '\0'))
            return cpsmatch_ai2016_main(a + 6);
        if (!strncmp(a, "mnn", 3) && (a[3] == ' ' || a[3] == '\0'))
            return cpsmatch_mnn_main(a + 3);
        if (!strncmp(a, "msmooth", 7) && (a[7] == ' ' || a[7] == '\0'))
            return cpsmatch_msmooth_main(a + 7);
    }
    if (cpsmatch_parse_options(args, &opts) != 0) {
        SF_error("cpsmatch: failed to parse arguments\n");
        return 198;
    }
    if (opts.out_nn_idx == 0) opts.out_nn_idx = opts.outcome_idx > 0 ? 7 : 6;
    if (opts.nobs == 0) {
        SF_error("cpsmatch: no observations\n");
        return 2000;
    }
    const int nearest = opts.method == MATCH_NEAREST;
    const int K = opts.n_neighbors;
    const int ny = opts.outcome_idx > 0 ? (opts.ny > 0 ? opts.ny : 1) : 0;
    const int with_self = nearest && opts.ai > 0 && opts.self_idx > 0 && ny > 0;
    const int nx = (nearest && opts.ai > 0 && opts.xcol_idx > 0 && opts.xself_idx > 0 && opts.nx > 0)
                       ? opts.nx : 0;
    const int nself = (with_self ? ny : 0) + nx;
    const int with_matched = ny > 0 && opts.matched_idx > 0;

    /* Load treatment, score, then the optional outcomes, support, arm,
     * support score and score-model columns */
    t_load = ctools_timer_seconds();
    int nvars = 0, col_y = -1, col_yeff = -1, col_presup = -1, col_arm = -1, col_sup = -1;
    int col_x = -1;
    int *var_indices = ctools_safe_malloc2((size_t)(6 + 2 * ny + nx), sizeof(int));
    if (!var_indices) {
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }
    var_indices[nvars++] = opts.treat_idx;
    var_indices[nvars++] = opts.pscore_idx;
    if (ny) {
        col_y = nvars;
        for (int k = 0; k < ny; k++) var_indices[nvars++] = opts.outcome_idx + k;
        if (opts.yeff_idx > 0) {
            col_yeff = nvars;
            for (int k = 0; k < ny; k++) var_indices[nvars++] = opts.yeff_idx + k;
        }
    }
    if (opts.presupport_idx > 0) { col_presup = nvars; var_indices[nvars++] = opts.presupport_idx; }
    if (opts.arm_idx > 0 && !opts.make_arm) { col_arm = nvars; var_indices[nvars++] = opts.arm_idx; }
    if (opts.supscore_idx > 0) { col_sup = nvars; var_indices[nvars++] = opts.supscore_idx; }
    if (nx) {
        col_x = nvars;
        for (int j = 0; j < nx; j++) var_indices[nvars++] = opts.xcol_idx + j;
    }

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    int load_flags =
        CTOOLS_LOAD_NO_SORT_ORDER | (opts.touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc =
        ctools_data_load(&filtered, var_indices, (size_t)nvars, 0, 0, load_flags);
    free(var_indices);
    if (load_rc != STATA_OK) {
        SF_error("cpsmatch: failed to load data\n");
        ctools_filtered_data_free(&filtered);
        return ctools_stata_rc(load_rc);
    }
    const size_t n_rows = filtered.data.nobs;
    const double *treatment = filtered.data.vars[0].data.dbl;
    const double *ps = filtered.data.vars[1].data.dbl;
    const double *presup = col_presup >= 0 ? filtered.data.vars[col_presup].data.dbl : NULL;
    const double *supscore = col_sup >= 0 ? filtered.data.vars[col_sup].data.dbl : ps;
    const perm_idx_t *obs_map = filtered.obs_map;
    const int threads = psm_threads(n_rows);
    t_load = ctools_timer_seconds() - t_load;
    t_setup = ctools_timer_seconds();

    /* Every buffer is allocated before the first store */
    const size_t nr = n_rows ? n_rows : 1;
    psm_outputs out = {0};
    out.neighbors = K;
    out.ny = ny;
    out.nself = nself;
    const double **ycol = ctools_safe_calloc2((size_t)(ny ? ny : 1), sizeof(double *));
    const double **yeff = ctools_safe_calloc2((size_t)(ny ? ny : 1), sizeof(double *));
    const double **xcol = ctools_safe_calloc2((size_t)(nx ? nx : 1), sizeof(double *));
    double **ypos = ctools_safe_calloc2((size_t)(ny + nx ? ny + nx : 1), sizeof(double *));
    out.matched = ctools_safe_calloc2((size_t)(ny ? ny : 1), sizeof(double *));
    out.self = ctools_safe_calloc2((size_t)(nself ? nself : 1), sizeof(double *));
    psm_order_t *order = ctools_safe_malloc2(nr, sizeof(psm_order_t));
    psm_order_t *order_tmp = ctools_safe_malloc2(nr, sizeof(psm_order_t));
    psm_radix_t *radix = ctools_safe_malloc2(nr, sizeof(psm_radix_t));
    psm_radix_t *radix_spare = ctools_safe_malloc2(nr, sizeof(psm_radix_t));
    size_t *hist = ctools_safe_malloc2((size_t)(threads > 1 ? threads : 1) * PSM_RADIX_BUCKETS,
                                       sizeof(size_t));
    size_t *chunk = ctools_safe_calloc2((size_t)threads + 1, sizeof(size_t));
    double *ps_pos = ctools_safe_malloc2(nr, sizeof(double));
    double *wpos = nearest ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    unsigned char *used = ctools_safe_malloc2(nr, 1);
    size_t *target_rows = ctools_safe_malloc2(nr, sizeof(size_t));
    size_t *donor_rows = ctools_safe_malloc2(nr, sizeof(size_t));
    sorted_control_t *donor_sorted = ctools_safe_malloc2(nr, sizeof(sorted_control_t));
    double *arm_made = opts.make_arm ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    cpsmatch_xsum *xwork =
        opts.effects ? ctools_safe_malloc2((size_t)threads * CPSMATCH_EFFECT_SUMS,
                                           sizeof(cpsmatch_xsum)) : NULL;
    psm_store_col *store_cols = ctools_safe_malloc2((size_t)(12 + K + 2 * ny + nself),
                                                    sizeof(psm_store_col));
    out.weight = ctools_safe_malloc2(nr, sizeof(double));
    out.support = ctools_safe_malloc2(nr, sizeof(double));
    out.nn = ctools_safe_malloc2(nr, sizeof(double));
    out.id = opts.out_match_idx > 0 ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    out.ids =
        opts.neighbor_ids_idx > 0 ? ctools_safe_malloc2(nr, (size_t)K * sizeof(double)) : NULL;
    out.pdif = opts.pdif_idx > 0 ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    int failed = !ycol || !yeff || !xcol || !ypos || !out.matched || !out.self || !order ||
                 !order_tmp || !radix || !radix_spare || !hist || !chunk || !ps_pos ||
                 (nearest && !wpos) || !used || !target_rows || !donor_rows || !donor_sorted ||
                 (opts.make_arm && !arm_made) || (opts.effects && !xwork) || !store_cols ||
                 !out.weight || !out.support || !out.nn || (opts.out_match_idx > 0 && !out.id) ||
                 (opts.neighbor_ids_idx > 0 && !out.ids) || (opts.pdif_idx > 0 && !out.pdif);
    for (int k = 0; k < ny && !failed; k++) {
        ycol[k] = filtered.data.vars[col_y + k].data.dbl;
        yeff[k] = col_yeff >= 0 ? filtered.data.vars[col_yeff + k].data.dbl : ycol[k];
        if (with_matched && !(out.matched[k] = ctools_safe_malloc2(nr, sizeof(double)))) failed = 1;
    }
    for (int j = 0; j < nx && !failed; j++) xcol[j] = filtered.data.vars[col_x + j].data.dbl;
    for (int j = 0; j < nself && !failed; j++)
        if (!(out.self[j] = ctools_safe_malloc2(nr, sizeof(double)))) failed = 1;
    if (nearest)
        for (int j = 0; j < ny + nx && !failed; j++)
            if (!(ypos[j] = ctools_safe_malloc2(nr, sizeof(double)))) failed = 1;
    /* rows outside the sample, for _id (and the one-row ai() arm) */
    psm_outside outside = {0};
    const int need_outside = nearest && (opts.outside_ids || (opts.single && opts.ai > 0));
    size_t nobs_all = 0;
    if (!failed && need_outside) {
        nobs_all = (size_t)SF_nobs();
        size_t n_out = nobs_all > n_rows ? nobs_all - n_rows : 0;
        size_t m = n_out ? n_out : 1;
        outside.obs = ctools_safe_malloc2(m, sizeof(perm_idx_t));
        outside.score = ctools_safe_malloc2(m, sizeof(double));
        outside.id = ctools_safe_malloc2(m, sizeof(double));
        psm_radix_t *ok = ctools_safe_malloc2(m, sizeof(psm_radix_t));
        psm_radix_t *os = ctools_safe_malloc2(m, sizeof(psm_radix_t));
        size_t *oh = ctools_safe_malloc2((size_t)psm_threads(n_out) * PSM_RADIX_BUCKETS, sizeof(size_t));
        if (!outside.obs || !outside.score || !outside.id || !ok || !os || !oh) failed = 1;
        else if ((rc = psm_outside_build(&outside, obs_map, n_rows, nobs_all, opts.pscore_idx,
                                         opts.descending, n_rows, ok, os, oh)))
            failed = 2;
        free(ok);
        free(os);
        free(oh);
    }
    psm_lists lists = {0};
    if (failed) {
        if (failed == 1) rc = 920;
        goto cleanup;
    }
    const double *arm = opts.make_arm ? arm_made
                        : col_arm >= 0 ? filtered.data.vars[col_arm].data.dbl : treatment;

    /* Arms, the control score range (r(common_min), r(common_max)) and the
     * support score's range in each arm */
    size_t n_treated = 0, n_controls = 0;
    double cmin = INFINITY, cmax = -INFINITY, smin0 = INFINITY, smax0 = -INFINITY;
    double smin1 = INFINITY, smax1 = -INFINITY;
    double *armw = arm_made;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)                  \
    reduction(+ : n_treated, n_controls) reduction(min : cmin, smin0, smin1)                    \
    reduction(max : cmax, smax0, smax1)
    for (size_t i = 0; i < n_rows; i++) {
        if (armw) armw[i] = psm_byte(treatment[i]);
        if (psm_missing(treatment[i], SV_missval) || psm_missing(ps[i], SV_missval)) continue;
        const double a = arm[i], s = supscore[i];
        if (a == 1.0) {
            n_treated++;
            if (!psm_missing(s, SV_missval)) {
                if (s < smin1) smin1 = s;
                if (s > smax1) smax1 = s;
            }
        } else if (a == 0.0) {
            n_controls++;
            if (ps[i] < cmin) cmin = ps[i];
            if (ps[i] > cmax) cmax = ps[i];
            if (!psm_missing(s, SV_missval)) {
                if (s < smin0) smin0 = s;
                if (s > smax0) smax0 = s;
            }
        }
    }
    if (!n_controls) cmin = cmax = SV_missval;
    if (smin0 > smax0) smin0 = smax0 = SV_missval;
    if (smin1 > smax1) smin1 = smax1 = SV_missval;

    /* psmatch2's variables before matching: _weight = _treated on support
     * (zero for the treated with ate), _nn zero on support, the matched
     * outcome zero for the rows to be matched. Common support (summarize
     * the score in the other arm, then replace _support = 0 outside its
     * range) compares as Stata does, missing above every number. */
    size_t n_valid = 0;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
        const size_t lo = n_rows * tid / team, hi = n_rows * (tid + 1) / team;
        size_t valid = 0;
        for (size_t i = lo; i < hi; i++)
            valid += !psm_missing(treatment[i], SV_missval) && !psm_missing(ps[i], SV_missval);
        chunk[tid + 1] = valid;
#pragma omp barrier
#pragma omp single
        {
            chunk[0] = 0;
            for (size_t t = 0; t < team; t++) chunk[t + 1] += chunk[t];
            n_valid = chunk[team];
        }
        size_t at = chunk[tid];
        for (size_t i = lo; i < hi; i++) {
            out.weight[i] = out.support[i] = out.nn[i] = SV_missval;
            if (out.id) out.id[i] = SV_missval;
            if (out.ids) for (int m = 0; m < K; m++) out.ids[i * (size_t)K + m] = SV_missval;
            if (out.pdif) out.pdif[i] = SV_missval;
            for (int k = 0; k < ny; k++) if (out.matched[k]) out.matched[k][i] = SV_missval;
            for (int j = 0; j < nself; j++) out.self[j][i] = SV_missval;
            if (psm_missing(treatment[i], SV_missval) || psm_missing(ps[i], SV_missval)) continue;
            const double a = arm[i], sc = supscore[i];
            int s = 1;
            if (presup) s = presup[i] == 1.0;
            else if (opts.common_support && a == 1.0) s = !(sc < smin0 || sc > smax0);
            else if (opts.common_support == 2 && a == 0.0) s = !(sc < smin1 || sc > smax1);
            order[at].arm = a;
            order[at].key = opts.descending ? -ps[i] : ps[i];
            order[at].obs = obs_map[i];
            order[at].row = (perm_idx_t)i;
            order[at].support = s;
            at++;
            out.support[i] = s;
            if (!s) continue;
            out.weight[i] = (opts.ate && a == 1.0) ? 0.0 : a;
            out.nn[i] = 0.0;
            if (opts.ate || a == 1.0)
                for (int k = 0; k < ny; k++) if (out.matched[k]) out.matched[k][i] = 0.0;
        }
    }
    if (opts.common_support == 0 && !presup) {
        /* nothing changes support before matching */
    }
    t_setup = ctools_timer_seconds() - t_setup;

    t_sort = ctools_timer_seconds();
    /* Bin the rows (stably), then sort each bin: the supported controls and
     * treated are contiguous blocks of the order. Kernel and radius
     * matching sort only the blocks. */
    size_t bin_start[6] = {0};
    psm_bin_rows(order, n_valid, order_tmp, bin_start, hist, psm_threads(n_valid));
    {
        psm_order_t *swap = order;
        order = order_tmp;
        order_tmp = swap;
    }
    for (int b = 0; b < 5; b++) {
        size_t n = bin_start[b + 1] - bin_start[b];
        if (b == 1 || b == 3)
            psm_sort_block(order + bin_start[b], n, order_tmp, radix, radix_spare, hist,
                           psm_threads(n));
        else if (nearest)
            qsort(order + bin_start[b], n, sizeof(psm_order_t), psm_compare_order);
    }
    const ssize_t c0 = (ssize_t)bin_start[1], c1 = (ssize_t)bin_start[2] - 1;
    const ssize_t t0 = (ssize_t)bin_start[3], t1 = (ssize_t)bin_start[4] - 1;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t p = 0; p < n_valid; p++) {
        const size_t row = order[p].row;
        ps_pos[p] = ps[row];
        if (out.id) out.id[row] = (double)p + 1.0;
        if (wpos) wpos[p] = out.weight[row];
        for (int k = 0; k < ny && nearest; k++) ypos[k][p] = ycol[k][row];
        for (int j = 0; j < nx; j++) ypos[ny + j][p] = xcol[j][row];
    }
    t_sort = ctools_timer_seconds() - t_sort;

    t_match = ctools_timer_seconds();
    const double cal = opts.caliper != 0.0 ? opts.caliper : INFINITY;
    int single_error = 0;
    if (nearest) {
        /* passes: with ate the controls to the treated first, then the
         * treated to the controls */
        for (int pass = opts.ate ? 0 : 1; pass < 2 && !rc; pass++) {
            ssize_t i0 = pass ? t0 : c0, i1 = pass ? t1 : c1;
            ssize_t j0 = pass ? c0 : t0, j1 = pass ? c1 : t1;
            size_t n_t = i1 >= i0 ? (size_t)(i1 - i0 + 1) : 0;
            int fast = opts.with_replace && !opts.ties && j1 >= j0;
            if (psm_lists_init(&lists, n_t, fast ? n_t * (size_t)K : n_t + 16)) { rc = 920; break; }
            if (fast && n_t) psm_nearest(ps_pos, i0, i1, j0, j1, K, cal, &lists, psm_threads(n_t));
            else if (n_t && psm_scan(ps_pos, (ssize_t)n_valid, i0, i1, j0, j1, 0, K, cal,
                                     opts.ties, !opts.with_replace, used, &lists)) {
                rc = 920;
                break;
            }
            if (n_t) psm_apply(&lists, i0, i1, order, ps_pos, wpos, ypos, &out, psm_threads(n_t));
            psm_lists_free(&lists);
        }
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
        for (size_t p = 0; p < n_valid; p++) out.weight[order[p].row] = wpos[p];
        /* ai(#): each arm matched to itself, controls then treated, for
         * every outcome and score-model column */
        if (!rc && nself) {
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t i = 0; i < n_rows; i++) {
                const double v = out.support[i] == 1.0 ? 0.0 : SV_missval;
                for (int j = 0; j < nself; j++) out.self[j][i] = v;
            }
            const int ny_self = with_self ? ny : 0;
            for (int pass = 0; pass < 2 && !rc; pass++) {
                ssize_t i0 = pass ? t0 : c0, i1 = pass ? t1 : c1;
                size_t n_t = i1 >= i0 ? (size_t)(i1 - i0 + 1) : 0;
                if (psm_lists_init(&lists, n_t, n_t + 16)) { rc = 920; break; }
                if (n_t && psm_scan(ps_pos, (ssize_t)n_valid, i0, i1, i0, i1, 1, opts.ai, cal,
                                    opts.ties, !opts.with_replace, used, &lists)) {
                    rc = 920;
                    break;
                }
                for (int j = 0; j < nself && n_t; j++) {
                    const double *src = j < ny_self ? ypos[j] : ypos[ny + (j - ny_self)];
                    psm_apply_self(&lists, i0, i1, order, src, j < ny_self && opts.altvariance,
                                   out.self[j], psm_threads(n_t));
                }
                psm_lists_free(&lists);
            }
        }
        if (!rc && opts.ai > 0 && opts.single) {
            const int ny_self = with_self ? ny : 0;
            /* psmatch2's same-arm pass on an arm with one row on support
             * before matching reads the row after it in the matching order
             * (a row of the other arm, off support or outside the sample)
             * and takes it as the match if its score is within the caliper;
             * with no row after it psmatch2 stops (r(3301)) and cpsmatch with
             * r(2001), after storing. */
            for (int a = 0; a < 2 && !rc && opts.single && !single_error; a++) {
                size_t count = 0, at = 0;
                for (size_t p = 0; p < n_valid; p++)
                    if (order[p].support && order[p].arm == (double)a) {
                        count++;
                        at = p;
                    }
                if (count != 1) continue;
                const size_t row = order[at].row;
                int next_in = 0;
                size_t next_row = 0, next_out = 0;
                if (at + 1 < n_valid) {
                    next_in = 1;
                    next_row = order[at + 1].row;
                } else {
                    size_t found = outside.n;
                    for (size_t j = 0; j < outside.n; j++)
                        if (outside.id[j] == (double)(n_valid + 1)) found = j;
                    if (found == outside.n) {
                        single_error = a ? 2 : 1;
                        break;
                    }
                    next_out = found;
                }
                if (out.support[row] != 1.0) continue;
                const double next_ps = next_in ? ps[next_row] : outside.score[next_out];
                const double dist = fabs(psm_sub(ps[row], next_ps));
                if (psm_missing(next_ps, SV_missval) ||
                    !(opts.caliper == 0.0 ? 1 : dist < opts.caliper))
                    continue;
                for (int j = 0; j < nself; j++) {
                    double value = SV_missval, own;
                    const int is_y = j < ny_self;
                    const double *colv = is_y ? ycol[j] : xcol[j - ny_self];
                    const int var = is_y ? opts.outcome_idx + j : opts.xcol_idx + (j - ny_self);
                    if (next_in) value = colv[next_row];
                    else if ((rc = SF_vdata(var, (ST_int)outside.obs[next_out], &value))) break;
                    if (opts.altvariance && is_y) {
                        own = colv[row];
                        const double center = psm_div(psm_add(value, own), 2.0);
                        const double d1 = psm_sub(own, center), d2 = psm_sub(value, center);
                        value = psm_add(psm_mul(d1, d1), psm_mul(d2, d2));
                    }
                    out.self[j][row] = value;
                }
            }
        }
    } else {
        /* kernel and radius: controls first with ate; their donors are the
         * supported treated, and the treated's donors are then the controls
         * still on support */
        for (int pass = opts.ate ? 0 : 1; pass < 2 && !rc; pass++) {
            ssize_t i0 = pass ? t0 : c0, i1 = pass ? t1 : c1;
            ssize_t j0 = pass ? c0 : t0, j1 = pass ? c1 : t1;
            size_t n_t = 0, n_d = 0;
            for (ssize_t p = i0; p <= i1; p++) target_rows[n_t++] = order[p].row;
            for (ssize_t p = j0; p <= j1; p++)
                if (out.support[order[p].row] == 1.0) donor_rows[n_d++] = order[p].row;
            psm_kernel_donors(ps, obs_map, donor_rows, n_d, donor_sorted);
            if (opts.method == MATCH_RADIUS) {
                /* psmatch2's radius is |dif| <= caliper: a negative caliper
                 * matches nothing, so only 0 (none given) takes the default */
                double radius = opts.caliper != 0.0 ? opts.caliper : CPSMATCH_DEFAULT_CALIPER;
                rc = psm_radius_pass(n_rows, ps, ycol, with_matched ? ny : 0, target_rows, n_t,
                                     donor_sorted, donor_rows, n_d, obs_map, radius, out.weight,
                                     out.nn, out.support, out.matched);
            } else {
                rc = psm_kernel_pass(&opts, ps, ycol, with_matched ? ny : 0, target_rows, n_t,
                                     donor_sorted, donor_rows, n_d, out.weight, out.nn,
                                     out.support, out.matched);
            }
        }
    }
    if (rc) goto cleanup;

    /* psmatch2: replace _weight = . if _weight==0 | _support==0, and no
     * matched outcome off support */
    size_t n_matched = 0, n_off = 0, n_controls_on = 0, n_controls_off = 0;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)                  \
    reduction(+ : n_matched, n_off, n_controls_on, n_controls_off)
    for (size_t i = 0; i < n_rows; i++) {
        if (out.weight[i] == 0.0 || out.support[i] == 0.0) out.weight[i] = SV_missval;
        if (out.support[i] == 0.0) {
            for (int k = 0; k < ny; k++) if (out.matched[k]) out.matched[k][i] = SV_missval;
            for (int j = 0; j < nself; j++) out.self[j][i] = SV_missval;
        }
        if (!psm_missing(treatment[i], SV_missval) && !psm_missing(ps[i], SV_missval)) {
            if (arm[i] == 1.0) {
                if (out.support[i] == 1.0) n_matched++;
                else n_off++;
            } else if (arm[i] == 0.0) {
                if (out.support[i] == 1.0) n_controls_on++;
                else n_controls_off++;
            }
        }
    }
    /* psmatch2's radius without a caliper is uniform-kernel matching with a
     * missing bandwidth: every target spreads its weight over all _N
     * observations, so on support every _weight gains (number of targets)
     * / _N (the ado file passes radall with an unreachable radius) */
    if (opts.radius_all) {
        size_t n_targets = 0;
        for (size_t i = 0; i < n_rows; i++)
            if (out.support[i] == 1.0 && (arm[i] == 1.0 || opts.ate)) n_targets++;
        const double share = (double)n_targets / (double)SF_nobs();
        for (size_t i = 0; i < n_rows; i++) {
            if (out.support[i] != 1.0) continue;
            double w = psm_add(opts.ate ? 0.0 : arm[i], share);
            out.weight[i] = w == 0.0 ? SV_missval : w;
        }
    }
    t_match = ctools_timer_seconds() - t_match;

    {
        double started = ctools_timer_seconds();
        int skip_last = 0;
        if (opts.statistics) {
            cpsmatch_statistics statistics;
            const double *y = ny ? ycol[ny - 1] : NULL;
            cpsmatch_compute_statistics(n_rows, arm, y, out.weight, out.support, &statistics);
            if (y && opts.float_squares) {
                /* psmatch2 squares the weights with gen w2 = _weight^2, in the
                 * default storage type: with float, each square is rounded */
                dd_real sum = {0.0, 0.0};
                for (size_t i = 0; i < n_rows; i++) {
                    double w = out.weight[i];
                    if (arm[i] != 0.0 || psm_missing(w, SV_missval)) continue;
                    float square = (float)(w * w);
                    if (isfinite(square)) sum = dd_add_d(sum, (double)square);
                }
                statistics.weight_squared = sum.hi + sum.lo;
            }
            /* the ado file takes these for the last outcome when no control
             * is off support */
            skip_last = statistics.n0_off == 0.0;
            if ((rc = cpsmatch_save_statistics(&statistics))) goto cleanup;
        }
        if (opts.effects && ny) {
            cpsmatch_effect_weights ew;
            cpsmatch_compute_effect_weights(n_rows, arm, out.support, out.weight, threads, xwork,
                                            &ew);
            const char *wn[4] = {"w2c", "w2t", "w1c", "w1t"};
            const double wv[4] = {ew.w2c, ew.w2t, ew.w1c, ew.w1t};
            char name[40];
            for (int j = 0; j < 4 && !rc; j++) {
                snprintf(name, sizeof(name), "_cpsmatch_e_%s", wn[j]);
                rc = SF_scal_save(name, wv[j]);
            }
            for (int k = 0; k < ny && !rc; k++) {
                if (skip_last && k == ny - 1) break;
                cpsmatch_effects e;
                cpsmatch_compute_effects(n_rows, arm, out.support, out.weight, yeff[k],
                                         out.matched[k], with_self ? out.self[k] : NULL, opts.ai,
                                         opts.ate, opts.samplevar, opts.altvariance, K, threads,
                                         xwork, &e);
                const char *names[17] = {"u1u", "n1u", "var1u", "u0u", "n0u", "var0u",
                                         "m1t", "n1",  "m0t",   "m1u", "m0u", "n0",
                                         "var1", "var0", "sumvt", "sumvu", "sumve"};
                const double values[17] = {e.u1u, e.n1u, e.var1u, e.u0u, e.n0u, e.var0u,
                                           e.m1t, e.n1,  e.m0t,   e.m1u, e.m0u, e.n0,
                                           e.var1, e.var0, e.sumvt, e.sumvu, e.sumve};
                for (int j = 0; j < 17 && !rc; j++) {
                    snprintf(name, sizeof(name), "_cpsmatch_e%d_%s", k + 1, names[j]);
                    rc = SF_scal_save(name, values[j]);
                }
            }
            if (rc) goto cleanup;
        }
        t_stats = ctools_timer_seconds() - started;
    }

    /* Ties in the score among the controls (with ate, also the treated)
     * where one of them is on support (warnings): psmatch2's check */
    if (opts.warn && nearest) {
        int dup = 0;
        for (int a = 0; a < (opts.ate ? 2 : 1) && !dup; a++) {
            size_t m = 0;
            for (size_t p = 0; p < n_valid; p++)
                if (order[p].arm == (double)a) {
                    radix[m].key = psm_radix_key(ps_pos[p]);
                    radix[m].at = (uint32_t)p;
                    m++;
                }
            psm_radix_t *sorted = psm_radix_sort(radix, radix_spare, m, hist, psm_threads(m));
            for (size_t r = 0; r < m && !dup;) {
                size_t e = r + 1;
                int on = order[sorted[r].at].support;
                while (e < m && ps_pos[sorted[e].at] == ps_pos[sorted[r].at]) {
                    on |= order[sorted[e].at].support;
                    e++;
                }
                if (e - r > 1 && on) dup = 1;
                r = e;
            }
        }
        if ((rc = SF_scal_save("_cpsmatch_dup", (double)dup))) goto cleanup;
    }

    t_store = ctools_timer_seconds();
    {
        size_t nc = 0;
        store_cols[nc++] = (psm_store_col){out.nn, opts.out_nn_idx};
        if (opts.out_weight_idx > 0) store_cols[nc++] = (psm_store_col){out.weight, opts.out_weight_idx};
        if (opts.out_support_idx > 0) store_cols[nc++] = (psm_store_col){out.support, opts.out_support_idx};
        if (out.id) store_cols[nc++] = (psm_store_col){out.id, opts.out_match_idx};
        if (out.pdif) store_cols[nc++] = (psm_store_col){out.pdif, opts.pdif_idx};
        for (int k = 0; k < ny; k++)
            if (out.matched[k]) store_cols[nc++] = (psm_store_col){out.matched[k], opts.matched_idx + k};
        for (int j = 0; j < nself; j++) {
            const int ny_self = with_self ? ny : 0;
            const int var = j < ny_self ? opts.self_idx + j : opts.xself_idx + (j - ny_self);
            store_cols[nc++] = (psm_store_col){out.self[j], var};
        }
        if (opts.make_arm && opts.arm_idx > 0) store_cols[nc++] = (psm_store_col){arm_made, opts.arm_idx};
        const size_t fill_rows = opts.fill ? (size_t)SF_nobs() : 0;
        rc = psm_store(store_cols, nc, n_rows, obs_map, fill_rows);
        /* _n1.._nK: column m of the row-major ids, through a strided copy */
        for (int m = 0; out.ids && m < K && !rc; m++) {
            double *col = ps_pos; /* free by now: reuse as the column */
            for (size_t i = 0; i < n_rows; i++) col[i] = out.ids[i * (size_t)K + m];
            psm_store_col one = {col, opts.neighbor_ids_idx + m};
            rc = psm_store(&one, 1, n_rows, obs_map, fill_rows);
        }
        if (!rc && out.id && opts.outside_ids && outside.n) {
            psm_store_col one = {outside.id, opts.out_match_idx};
            rc = psm_store(&one, 1, outside.n, outside.obs, 0);
        }
    }
    if (rc) goto cleanup;
    t_store = ctools_timer_seconds() - t_store;

    if ((rc = SF_scal_save("_cpsmatch_time_load", t_load))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_time_setup", t_setup))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_time_sort", t_sort))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_time_match", t_match))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_time_store", t_store))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_time_stats", t_stats))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_time_total", ctools_timer_seconds() - t_start))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_n_treated", (double)n_treated))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_n_controls", (double)n_controls))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_n_matched", (double)n_matched))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_n_off_support", (double)n_off))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_n_controls_on", (double)n_controls_on))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_n_controls_off", (double)n_controls_off))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_common_min", cmin))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_common_max", cmax))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_threads_max", (double)ctools_get_max_threads()))) goto cleanup;
    if ((rc = SF_scal_save("_cpsmatch_openmp_enabled", (double)CTOOLS_OPENMP_ENABLED))) goto cleanup;
    /* Revision of the ai() same-arm pass: 2 since noreplacement takes only
     * the first match without replacement. The ado file refuses the results
     * of ai(#) above 1 with noreplacement from a plugin without it. */
    if ((rc = SF_scal_save("_cpsmatch_ai_pass", 2.0))) goto cleanup;
    if (single_error) {
        SF_error(single_error == 2
                     ? "cpsmatch: ai() has no same-arm match for the only treated observation on support\n"
                     : "cpsmatch: ai() has no same-arm match for the only untreated observation on support\n");
        rc = 2001;
    }

cleanup:
    if (rc == 920) SF_error("cpsmatch: memory allocation failed\n");
    psm_lists_free(&lists);
    free(order);
    free(order_tmp);
    free(radix);
    free(radix_spare);
    free(hist);
    free(chunk);
    free(ps_pos);
    free(wpos);
    free(used);
    free(target_rows);
    free(donor_rows);
    free(donor_sorted);
    free(arm_made);
    free(xwork);
    free(store_cols);
    free(out.weight);
    free(out.support);
    free(out.nn);
    free(out.id);
    free(out.ids);
    free(out.pdif);
    for (int k = 0; out.matched && k < ny; k++) free(out.matched[k]);
    for (int j = 0; out.self && j < nself; j++) free(out.self[j]);
    for (int j = 0; ypos && j < ny + nx; j++) free(ypos[j]);
    free(out.matched);
    free(out.self);
    free(ypos);
    free((void *)ycol);
    free((void *)yeff);
    free((void *)xcol);
    free(outside.obs);
    free(outside.score);
    free(outside.id);
    ctools_filtered_data_free(&filtered);
    return rc;
}
