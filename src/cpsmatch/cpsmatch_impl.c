#include "cpsmatch_internal.h"
#include "cpsmatch_impl.h"
#include "cpsmatch_stats.h"
#include "cpsmatch_index.h"
#include "ctools_ols.h"

/* Propensity-score matching as psmatch2 4.0.12 performs it.
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
 * - ties adds the (unused) donors that follow it in that direction at no
 *   greater distance, and neighbor(k) then adds the k-1 closest around it
 *   from both sides, the left on equal distance, so a donor can enter a
 *   list twice;
 * - every list entry gets 1/(list length), and a match needs a distance below
 *   caliper() (strict); a target left without one is off support;
 * - without replacement the matched donors are used up and the scan resumes
 *   at the next unused donor; once it runs off the donors the remaining
 *   targets are off support.
 * With replacement and without ties the scan finds what a binary search
 * finds (the closer neighbor, the lower on equal distance, the first of a run
 * of equal scores), so those passes run in parallel.
 *
 * Kernel and radius matching weigh every donor of the other arm within the
 * bandwidth; with ate the controls are matched first, and a control left
 * without a match is then no longer a donor for the treated, as in psmatch2. */

/* Stata's missing values, tested without a host callback per element */
static inline int psm_missing(double x, double missval)
{
    return !(x < missval);
}

/* One in-sample row in psmatch2's matching order */
typedef struct
{
    double arm;      /* stored _treated */
    double key;      /* pscore, negated with descending; zero for kernel matching */
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
 * order. tmp holds n rows; keys and spare hold n entries each. */
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

/* Sort keys[0..n) stably by key; returns the array holding the result
 * (keys or spare) */
static psm_radix_t *psm_radix_sort(psm_radix_t *keys, psm_radix_t *spare, size_t n)
{
    enum { BITS = 11, BUCKETS = 1 << BITS };
    size_t count[BUCKETS];
    if (n < 2) return keys;
    for (int shift = 0; shift < 64; shift += BITS) {
        uint64_t first = (keys[0].key >> shift) & (BUCKETS - 1);
        int same = 1;
        memset(count, 0, sizeof(count));
        for (size_t i = 0; i < n; i++) {
            uint64_t digit = (keys[i].key >> shift) & (BUCKETS - 1);
            same &= digit == first;
            count[digit]++;
        }
        if (same) continue;
        size_t sum = 0;
        for (size_t b = 0; b < BUCKETS; b++) {
            size_t c = count[b];
            count[b] = sum;
            sum += c;
        }
        for (size_t i = 0; i < n; i++)
            spare[count[(keys[i].key >> shift) & (BUCKETS - 1)]++] = keys[i];
        psm_radix_t *swap = keys;
        keys = spare;
        spare = swap;
    }
    return keys;
}

static void psm_sort_block(psm_order_t *block, size_t n, psm_order_t *tmp, psm_radix_t *keys,
                           psm_radix_t *spare)
{
    if (n < 2) return;
    for (size_t i = 0; i < n; i++) {
        keys[i].key = psm_radix_key(block[i].key);
        keys[i].at = (uint32_t)i;
    }
    keys = psm_radix_sort(keys, spare, n);
    for (size_t i = 0; i < n; i++) tmp[i] = block[keys[i].at];
    memcpy(block, tmp, n * sizeof(psm_order_t));
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
            if (ties) {
                for (ssize_t k = psm_step(jm, forward, me, noreplace, used, n_pos);
                     k >= j0 && k <= j1 && fabs(ps[i] - ps[k]) <= nearest;
                     k = psm_step(k, forward, me, noreplace, used, n_pos)) {
                    if (psm_lists_push(out, k)) return -1;
                    n++;
                    if (noreplace) used[k] = 1;
                }
            }
            ssize_t p0 = jm - 1, p1 = jm + 1;
            if (p0 == me) p0--;
            if (p1 == me) p1++;
            double e0 = 0.0, e1 = 0.0;
            for (int k = 1; k < neighbors;) {
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
                if (!(fabs(ps[i] - ps[pick]) < cal)) break;
                if (psm_lists_push(out, pick)) return -1;
                n++;
                k++;
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
                        int neighbors, double cal, psm_lists *out)
{
    const double *donor = ps + j0;
    const size_t n_donors = (size_t)(j1 - j0 + 1);
    const ssize_t n_targets = i1 - i0 + 1;
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads())
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
    double *weight, *support, *nn, *id, *ids, *pdif, *matched, *self;
    int neighbors;
} psm_outputs;

/* Record a matching pass between blocks (not ai()): psmatch2's _n#, _nn,
 * _weight, _pdif and matched outcome, and support for unmatched targets. */
static void psm_apply(const psm_lists *l, ssize_t i0, ssize_t i1, const psm_order_t *order,
                      const double *ps, const double *y, psm_outputs *o)
{
    for (ssize_t i = i0; i <= i1; i++) {
        size_t t = (size_t)(i - i0), row = order[i].row;
        uint32_t n = l->count[t];
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
        if (o->pdif) o->pdif[row] = fabs(ps[i] - ps[list[0]]);
        double matched = 0.0;
        int missing = 0;
        for (uint32_t m = 0; m < n; m++) {
            size_t donor = order[list[m]].row;
            o->weight[donor] += 1.0 / (double)n;
            if (y) {
                if (psm_missing(y[donor], SV_missval)) missing = 1;
                else matched = matched + y[donor] / (double)n;
            }
        }
        if (o->matched) o->matched[row] = missing ? SV_missval : matched;
    }
}

/* Record an ai() pass: the mean outcome of the same-arm matches, or with
 * altvariance the conditional variance of Abadie et al. (2004, p. 303) */
static void psm_apply_self(const psm_lists *l, ssize_t i0, ssize_t i1, const psm_order_t *order,
                           const double *y, int altvariance, double *self)
{
    for (ssize_t i = i0; i <= i1; i++) {
        size_t t = (size_t)(i - i0), row = order[i].row;
        uint32_t n = l->count[t];
        if (n == 0 || psm_missing(self[row], SV_missval)) continue;
        const perm_idx_t *list = l->item + l->start[t];
        double mean = self[row];
        int missing = psm_missing(y[row], SV_missval);
        for (uint32_t m = 0; m < n; m++) {
            double v = y[order[list[m]].row];
            if (psm_missing(v, SV_missval)) missing = 1;
            else mean = mean + v / (double)n;
        }
        if (missing) {
            self[row] = SV_missval;
            continue;
        }
        if (altvariance) {
            double center = ((double)n * mean + y[row]) / ((double)n + 1.0);
            double sum = (y[row] - center) * (y[row] - center);
            for (uint32_t m = 0; m < n; m++) {
                double d = y[order[list[m]].row] - center;
                sum = sum + d * d;
            }
            mean = sum / (double)n;
        }
        self[row] = mean;
    }
}

/* Kernel matching pass: targets are filtered rows; donor_sorted (donor_rows
 * by orig_idx) is ordered by score. A target weighs the donors with
 * |dif| <= bandwidth (the normal kernel: all, within the window where its
 * weight underflows to zero), normalized to sum to one; the weights add into
 * weight, and targets get _nn (positive weights), support and the weighted
 * mean outcome. Donor scores, outcomes and weight sums are kept in sorted
 * order, so a target reads one contiguous run. */
static int psm_kernel_pass(const cpsmatch_options *opts, const double *ps, const double *y,
                           const size_t *targets, size_t n_targets,
                           const sorted_control_t *donor_sorted, const size_t *donor_rows,
                           size_t n_donors, double *weight, double *nn, double *support,
                           double *matched)
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
    /* exp(-u^2/2) underflows to exactly 0 beyond |u| = 38.6 */
    const double window = opts->kernel == KERNEL_NORMAL ? 38.7 * bw : bw;
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = ctools_get_max_threads();
#endif
    int failed = 0;
    double *dy = ctools_safe_malloc2(n_donors, sizeof(double));
    double **local = ctools_safe_calloc2((size_t)n_threads, sizeof(double *));
    /* a target's kernel weights, in its thread's workspace */
    match_workspace *ws = ctools_safe_calloc2((size_t)n_threads, sizeof(match_workspace));
    if (!dy || !local || !ws) failed = 1;
    for (int k = 0; k < n_threads && !failed; k++) {
        local[k] = ctools_safe_calloc2(n_donors, sizeof(double));
        cpsmatch_match_workspace_init(&ws[k], MATCH_WORKSPACE_INITIAL_CAP);
        if (!local[k] || !ws[k].match_weights) failed = 1;
    }
    if (!failed) {
        for (size_t j = 0; j < n_donors; j++)
            dy[j] = y ? y[donor_rows[donor_sorted[j].orig_idx]] : 0.0;
#pragma omp parallel num_threads(n_threads) reduction(| : failed)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
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
                double sum_w = 0.0, sum_wy = 0.0;
                int n = 0, missing = 0;
                for (size_t j = start; j < end; j++) {
                    if (!(w[j - start] > 0.0)) continue;
                    double wn = w[j - start] / total;
                    local[tid][j] += wn;
                    if (psm_missing(dy[j], missval)) missing = 1;
                    sum_wy += wn * dy[j];
                    sum_w += wn;
                    n++;
                }
                nn[row] = (double)n;
                if (matched) matched[row] = missing ? missval : sum_wy / sum_w;
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
    free(local);
    free(ws);
    return failed ? 920 : 0;
}

/* Radius matching pass: each donor within the radius gets 1/count, through
 * the exact endpoint sums of cpsmatch_match_radius; the matched outcome is a
 * range mean over double-double prefix sums. */
static int psm_radius_pass(size_t n_rows, const double *ps, const double *y, const size_t *targets,
                           size_t n_targets, const sorted_control_t *donor_sorted,
                           const size_t *donor_rows, size_t n_donors, const perm_idx_t *obs_map,
                           double radius, double *weight, double *nn, double *support,
                           double *matched)
{
    double *target_ps = ctools_safe_malloc2(n_targets ? n_targets : 1, sizeof(double));
    double *scratch_w = ctools_safe_calloc2(n_rows ? n_rows : 1, sizeof(double));
    double *scratch_m = ctools_safe_malloc2(n_rows ? n_rows : 1, sizeof(double));
    dd_real *prefix = y ? ctools_safe_calloc2(n_donors + 1, sizeof(dd_real)) : NULL;
    size_t *missing = y ? ctools_safe_calloc2(n_donors + 1, sizeof(size_t)) : NULL;
    perm_idx_t *starts = ctools_safe_malloc2(n_targets ? n_targets : 1, sizeof(perm_idx_t));
    int rc = 0;
    if (!target_ps || !scratch_w || !scratch_m || !starts || (y && (!prefix || !missing))) {
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
    if (y) {
        for (size_t j = 0; j < n_donors; j++) {
            double v = y[donor_rows[donor_sorted[j].orig_idx]];
            int absent = psm_missing(v, SV_missval);
            prefix[j + 1] = dd_add_d(prefix[j], absent ? 0.0 : v);
            missing[j + 1] = missing[j] + (size_t)absent;
        }
        for (size_t t = 0; t < n_targets; t++) {
            size_t row = targets[t];
            if (nn[row] <= 0.0) continue;
            size_t start = starts[t], end = start + (size_t)nn[row];
            if (missing[end] != missing[start]) {
                matched[row] = SV_missval;
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
                matched[row] = mean.hi + mean.lo;
            } else {
                matched[row] = sum / nn[row];
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

/* Donors of a kernel pass, ordered by score; rows are in observation
 * order, which the stable sort keeps among equal scores */
static void psm_kernel_donors(const double *ps, const perm_idx_t *obs_map, const size_t *rows,
                              size_t n, sorted_control_t *sorted, psm_radix_t *keys,
                              psm_radix_t *spare)
{
    for (size_t j = 0; j < n; j++) {
        keys[j].key = psm_radix_key(ps[rows[j]]);
        keys[j].at = (uint32_t)j;
    }
    keys = psm_radix_sort(keys, spare, n);
    for (size_t k = 0; k < n; k++) {
        size_t j = keys[k].at;
        sorted[k].orig_idx = j;
        sorted[k].obs_idx = obs_map[rows[j]];
        sorted[k].pscore = ps[rows[j]];
    }
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

    /* Load treatment, score, then the optional outcome, support and arm */
    t_load = ctools_timer_seconds();
    int var_indices[5], nvars = 0, col_outcome = -1, col_presup = -1, col_arm = -1;
    var_indices[nvars++] = opts.treat_idx;
    var_indices[nvars++] = opts.pscore_idx;
    if (opts.outcome_idx > 0) { col_outcome = nvars; var_indices[nvars++] = opts.outcome_idx; }
    if (opts.presupport_idx > 0) { col_presup = nvars; var_indices[nvars++] = opts.presupport_idx; }
    if (opts.arm_idx > 0) { col_arm = nvars; var_indices[nvars++] = opts.arm_idx; }

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    int load_flags =
        CTOOLS_LOAD_NO_SORT_ORDER | (opts.touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc =
        ctools_data_load(&filtered, var_indices, (size_t)nvars, 0, 0, load_flags);
    if (load_rc != STATA_OK) {
        SF_error("cpsmatch: failed to load data\n");
        ctools_filtered_data_free(&filtered);
        return ctools_stata_rc(load_rc);
    }
    const size_t n_rows = filtered.data.nobs;
    const double *treatment = filtered.data.vars[0].data.dbl;
    const double *ps = filtered.data.vars[1].data.dbl;
    const double *y = col_outcome >= 0 ? filtered.data.vars[col_outcome].data.dbl : NULL;
    const double *presup = col_presup >= 0 ? filtered.data.vars[col_presup].data.dbl : NULL;
    /* psmatch2 matches on _treated (gen byte _treated = treatvar): treated
     * are _treated == 1 and controls _treated == 0. Rows with another code
     * stay in the sample, on support, with _weight = _treated. */
    const double *arm = col_arm >= 0 ? filtered.data.vars[col_arm].data.dbl : treatment;
    const perm_idx_t *obs_map = filtered.obs_map;
    t_load = ctools_timer_seconds() - t_load;
    t_setup = ctools_timer_seconds();

    /* Every buffer is allocated before the first store */
    const size_t nr = n_rows ? n_rows : 1;
    psm_outputs out = {0};
    out.neighbors = K;
    psm_order_t *order = ctools_safe_malloc2(nr, sizeof(psm_order_t));
    psm_order_t *order_tmp = ctools_safe_malloc2(nr, sizeof(psm_order_t));
    psm_radix_t *radix = ctools_safe_malloc2(nr, sizeof(psm_radix_t));
    psm_radix_t *radix_spare = ctools_safe_malloc2(nr, sizeof(psm_radix_t));
    double *ps_pos = ctools_safe_malloc2(nr, sizeof(double));
    unsigned char *used = ctools_safe_malloc2(nr, 1);
    size_t *target_rows = ctools_safe_malloc2(nr, sizeof(size_t));
    size_t *donor_rows = ctools_safe_malloc2(nr, sizeof(size_t));
    sorted_control_t *donor_sorted = ctools_safe_malloc2(nr, sizeof(sorted_control_t));
    out.weight = ctools_safe_malloc2(nr, sizeof(double));
    out.support = ctools_safe_malloc2(nr, sizeof(double));
    out.nn = ctools_safe_malloc2(nr, sizeof(double));
    out.id = opts.out_match_idx > 0 ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    out.ids =
        opts.neighbor_ids_idx > 0 ? ctools_safe_malloc2(nr, (size_t)K * sizeof(double)) : NULL;
    out.pdif = opts.pdif_idx > 0 ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    out.matched = y && opts.matched_idx > 0 ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    out.self =
        y && opts.self_idx > 0 && opts.ai > 0 ? ctools_safe_malloc2(nr, sizeof(double)) : NULL;
    psm_lists lists = {0};
    if (!order || !order_tmp || !radix || !radix_spare || !ps_pos || !used || !target_rows ||
        !donor_rows || !donor_sorted || !out.weight || !out.support || !out.nn ||
        (opts.out_match_idx > 0 && !out.id) ||
        (opts.neighbor_ids_idx > 0 && !out.ids) || (opts.pdif_idx > 0 && !out.pdif) ||
        (y && opts.matched_idx > 0 && !out.matched) ||
        (y && opts.self_idx > 0 && opts.ai > 0 && !out.self)) {
        rc = 920;
        goto cleanup;
    }

    /* Arms and the control score range (r(common_min), r(common_max)) */
    size_t n_treated = 0, n_controls = 0;
    double cmin = SV_missval, cmax = SV_missval;
    for (size_t i = 0; i < n_rows; i++) {
        if (psm_missing(treatment[i], SV_missval) || psm_missing(ps[i], SV_missval)) continue;
        if (arm[i] == 1.0) n_treated++;
        else if (arm[i] == 0.0) {
            n_controls++;
            if (psm_missing(cmin, SV_missval) || ps[i] < cmin) cmin = ps[i];
            if (psm_missing(cmax, SV_missval) || ps[i] > cmax) cmax = ps[i];
        }
    }

    /* psmatch2's variables before matching: _weight = _treated on support
     * (zero for the treated with ate), _nn zero on support, the matched
     * outcome zero for the rows to be matched. */
    size_t n_valid = 0;
    for (size_t i = 0; i < n_rows; i++) {
        out.weight[i] = out.support[i] = out.nn[i] = SV_missval;
        if (out.id) out.id[i] = SV_missval;
        if (out.ids) for (int m = 0; m < K; m++) out.ids[i * (size_t)K + m] = SV_missval;
        if (out.pdif) out.pdif[i] = SV_missval;
        if (out.matched) out.matched[i] = SV_missval;
        if (out.self) out.self[i] = SV_missval;
        if (psm_missing(treatment[i], SV_missval) || psm_missing(ps[i], SV_missval)) continue;
        int s = 1;
        if (presup) s = presup[i] == 1.0;
        else if (opts.common_support && arm[i] == 1.0) s = !(ps[i] < cmin || ps[i] > cmax);
        order[n_valid].arm = arm[i];
        order[n_valid].key = nearest ? (opts.descending ? -ps[i] : ps[i]) : 0.0;
        order[n_valid].obs = obs_map[i];
        order[n_valid].row = (perm_idx_t)i;
        order[n_valid].support = s;
        n_valid++;
        out.support[i] = s;
        if (!s) continue;
        out.weight[i] = (opts.ate && arm[i] == 1.0) ? 0.0 : arm[i];
        out.nn[i] = 0.0;
        if (out.matched && (opts.ate || arm[i] == 1.0)) out.matched[i] = 0.0;
    }
    t_setup = ctools_timer_seconds() - t_setup;

    t_sort = ctools_timer_seconds();
    /* Bin the rows (stably), then sort each bin: the supported controls and
     * treated are contiguous blocks of the order. Kernel and radius
     * matching need only the blocks. */
    size_t bin_start[6] = {0};
    for (size_t p = 0; p < n_valid; p++) bin_start[psm_bin(&order[p]) + 1]++;
    for (int b = 0; b < 5; b++) bin_start[b + 1] += bin_start[b];
    {
        size_t fill[5];
        memcpy(fill, bin_start, sizeof(fill));
        for (size_t p = 0; p < n_valid; p++) order_tmp[fill[psm_bin(&order[p])]++] = order[p];
        psm_order_t *swap = order;
        order = order_tmp;
        order_tmp = swap;
    }
    if (nearest) {
        for (int b = 0; b < 5; b++) {
            size_t n = bin_start[b + 1] - bin_start[b];
            if (b == 1 || b == 3)
                psm_sort_block(order + bin_start[b], n, order_tmp, radix, radix_spare);
            else
                qsort(order + bin_start[b], n, sizeof(psm_order_t), psm_compare_order);
        }
    }
    const ssize_t c0 = (ssize_t)bin_start[1], c1 = (ssize_t)bin_start[2] - 1;
    const ssize_t t0 = (ssize_t)bin_start[3], t1 = (ssize_t)bin_start[4] - 1;
    for (size_t p = 0; p < n_valid; p++) {
        ps_pos[p] = ps[order[p].row];
        if (out.id) out.id[order[p].row] = (double)p + 1.0;
    }
    t_sort = ctools_timer_seconds() - t_sort;

    t_match = ctools_timer_seconds();
    const double cal = opts.caliper != 0.0 ? opts.caliper : INFINITY;
    if (nearest) {
        /* passes: with ate the controls to the treated first, then the
         * treated to the controls */
        for (int pass = opts.ate ? 0 : 1; pass < 2 && !rc; pass++) {
            ssize_t i0 = pass ? t0 : c0, i1 = pass ? t1 : c1;
            ssize_t j0 = pass ? c0 : t0, j1 = pass ? c1 : t1;
            size_t n_t = i1 >= i0 ? (size_t)(i1 - i0 + 1) : 0;
            int fast = opts.with_replace && !opts.ties && j1 >= j0;
            if (psm_lists_init(&lists, n_t, fast ? n_t * (size_t)K : n_t + 16)) { rc = 920; break; }
            if (fast && n_t) psm_nearest(ps_pos, i0, i1, j0, j1, K, cal, &lists);
            else if (n_t && psm_scan(ps_pos, (ssize_t)n_valid, i0, i1, j0, j1, 0, K, cal,
                                     opts.ties, !opts.with_replace, used, &lists)) {
                rc = 920;
                break;
            }
            psm_apply(&lists, i0, i1, order, ps_pos, y, &out);
            psm_lists_free(&lists);
        }
        /* ai(#): each arm matched to itself, controls then treated */
        if (!rc && out.self) {
            for (size_t i = 0; i < n_rows; i++)
                out.self[i] = out.support[i] == 1.0 ? 0.0 : SV_missval;
            for (int pass = 0; pass < 2 && !rc; pass++) {
                ssize_t i0 = pass ? t0 : c0, i1 = pass ? t1 : c1;
                size_t n_t = i1 >= i0 ? (size_t)(i1 - i0 + 1) : 0;
                if (psm_lists_init(&lists, n_t, n_t + 16)) { rc = 920; break; }
                if (n_t && psm_scan(ps_pos, (ssize_t)n_valid, i0, i1, i0, i1, 1, opts.ai, cal,
                                    opts.ties, !opts.with_replace, used, &lists)) {
                    rc = 920;
                    break;
                }
                psm_apply_self(&lists, i0, i1, order, y, opts.altvariance, out.self);
                psm_lists_free(&lists);
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
            psm_kernel_donors(ps, obs_map, donor_rows, n_d, donor_sorted, radix, radix_spare);
            if (opts.method == MATCH_RADIUS) {
                /* psmatch2's radius is |dif| <= caliper: a negative caliper
                 * matches nothing, so only 0 (none given) takes the default */
                double radius = opts.caliper != 0.0 ? opts.caliper : CPSMATCH_DEFAULT_CALIPER;
                rc = psm_radius_pass(n_rows, ps, y, target_rows, n_t, donor_sorted, donor_rows,
                                     n_d, obs_map, radius, out.weight, out.nn, out.support,
                                     out.matched);
            } else {
                rc = psm_kernel_pass(&opts, ps, y, target_rows, n_t, donor_sorted, donor_rows,
                                     n_d, out.weight, out.nn, out.support, out.matched);
            }
        }
    }
    if (rc) goto cleanup;

    /* psmatch2: replace _weight = . if _weight==0 | _support==0, and no
     * matched outcome off support */
    size_t n_matched = 0, n_off = 0, n_controls_on = 0, n_controls_off = 0;
    for (size_t i = 0; i < n_rows; i++) {
        if (out.weight[i] == 0.0 || out.support[i] == 0.0) out.weight[i] = SV_missval;
        if (out.support[i] == 0.0) {
            if (out.matched) out.matched[i] = SV_missval;
            if (out.self) out.self[i] = SV_missval;
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
    t_match = ctools_timer_seconds() - t_match;

    if (opts.statistics) {
        cpsmatch_statistics statistics;
        double started = ctools_timer_seconds();
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
        rc = cpsmatch_save_statistics(&statistics);
        t_stats = ctools_timer_seconds() - started;
        if (rc) goto cleanup;
    }

    t_store = ctools_timer_seconds();
    for (size_t i = 0; i < n_rows && !rc; i++) {
        ST_int obs = (ST_int)obs_map[i];
        rc = SF_vstore(opts.out_nn_idx, obs, out.nn[i]);
        if (!rc && opts.out_weight_idx > 0) rc = SF_vstore(opts.out_weight_idx, obs, out.weight[i]);
        if (!rc && opts.out_support_idx > 0)
            rc = SF_vstore(opts.out_support_idx, obs, out.support[i]);
        if (!rc && out.id) rc = SF_vstore(opts.out_match_idx, obs, out.id[i]);
        for (int m = 0; out.ids && m < K && !rc; m++)
            rc = SF_vstore(opts.neighbor_ids_idx + m, obs, out.ids[i * (size_t)K + m]);
        if (!rc && out.pdif) rc = SF_vstore(opts.pdif_idx, obs, out.pdif[i]);
        if (!rc && out.matched) rc = SF_vstore(opts.matched_idx, obs, out.matched[i]);
        if (!rc && out.self) rc = SF_vstore(opts.self_idx, obs, out.self[i]);
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

cleanup:
    if (rc == 920) SF_error("cpsmatch: memory allocation failed\n");
    psm_lists_free(&lists);
    free(order);
    free(order_tmp);
    free(radix);
    free(radix_spare);
    free(ps_pos);
    free(used);
    free(target_rows);
    free(donor_rows);
    free(donor_sorted);
    free(out.weight);
    free(out.support);
    free(out.nn);
    free(out.id);
    free(out.ids);
    free(out.pdif);
    free(out.matched);
    free(out.self);
    ctools_filtered_data_free(&filtered);
    return rc;
}
