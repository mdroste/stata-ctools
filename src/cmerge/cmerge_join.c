/*
 * cmerge_join.c
 * Sorted merge join algorithm for cmerge
 *
 * Both key sets arrive sorted. The join runs in chunks: the master is cut at
 * the first row of a key group, the using data at the first row whose key is
 * not below that group's key, so every key falls in one chunk on both sides
 * and the chunks' outputs, concatenated, are the serial join's output.
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "cmerge_join.h"
#include "cmerge_keys.h"
#include "cmerge_group_search.h"
#include "../ctools_config.h"  /* For ctools_safe_mul_size */
#ifdef _OPENMP
#include <omp.h>
#endif

#define CM_JOIN_MAX_KEYS 32
/* Rows (both sides) a chunk should have before another worker is worth it */
#ifndef CM_JOIN_ROWS_PER_CHUNK
#define CM_JOIN_ROWS_PER_CHUNK 262144
#endif
#define CM_JOIN_MAX_CHUNKS 64

/* The key columns of one side. Numeric keys compare as doubles: every Stata
 * missing value is a double above all nonmissing values and the missing
 * values ascend . < .a < ... < .z, so plain comparison is merge's order;
 * strings compare bytewise (strcmp, unsigned), NULL as "". */
typedef struct {
    const double *num[CM_JOIN_MAX_KEYS];
    char *const *str[CM_JOIN_MAX_KEYS];
} cm_side_t;

typedef struct {
    int nkeys;
    int single_num;            /* one numeric key */
    uint8_t is_str[CM_JOIN_MAX_KEYS];
    cm_side_t m, u;
    cmerge_type_t type;
} cm_join_t;

static inline int cm_cmp(const cm_join_t *J, const cm_side_t *A, size_t a,
                         const cm_side_t *B, size_t b)
{
    if (J->single_num) {
        double x = A->num[0][a], y = B->num[0][b];
        return (x > y) - (x < y);
    }
    for (int k = 0; k < J->nkeys; k++) {
        if (!J->is_str[k]) {
            double x = A->num[k][a], y = B->num[k][b];
            if (x < y) return -1;
            if (x > y) return 1;
        } else {
            const char *x = A->str[k][a], *y = B->str[k][b];
            int c = strcmp(x ? x : "", y ? y : "");
            if (c) return c < 0 ? -1 : 1;
        }
    }
    return 0;
}

/* 1 when adjacent rows i-1 and i of one side share a key, for some i in [start, end) */
static int cm_has_duplicate(const cm_join_t *J, const cm_side_t *S, size_t start, size_t end)
{
    if (J->single_num) {
        const double *v = S->num[0];
        for (size_t i = start; i < end; i++)
            if (v[i - 1] == v[i]) return 1;
        return 0;
    }
    for (size_t i = start; i < end; i++)
        if (cm_cmp(J, S, i - 1, S, i) == 0) return 1;
    return 0;
}

static int cm_side_unique(const cm_join_t *J, const cm_side_t *S, size_t n, int chunks)
{
    if (n < 2) return 1;
    int dup = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(chunks) reduction(|:dup) if (chunks > 1)
#endif
    for (int c = 0; c < chunks; c++) {
        size_t lo = 1 + (n - 1) * (size_t)c / (size_t)chunks;
        size_t hi = 1 + (n - 1) * (size_t)(c + 1) / (size_t)chunks;
        dup |= cm_has_duplicate(J, S, lo, hi);
    }
    (void)chunks;
    return !dup;
}

/* First using row in [lo, hi) whose key is not below master row m's key */
static size_t cm_using_lower_bound(const cm_join_t *J, size_t m, size_t lo, size_t hi)
{
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (cm_cmp(J, &J->u, mid, &J->m, m) < 0) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static inline void cm_emit(cmerge_output_spec_t *o, int32_t m, int32_t u, int8_t r)
{
    o->master_sorted_row = m;
    o->using_sorted_row = u;
    o->merge_result = r;
}

/* The serial join of master rows [m0, m1) and using rows [u0, u1). Writes at
 * most (m1 - m0) + (u1 - u0) rows to out and returns their number. */
static size_t cm_join_chunk(const cm_join_t *J, size_t m0, size_t m1, size_t u0, size_t u1,
                            cmerge_output_spec_t *out)
{
    size_t n = 0, mi = m0, ui = u0;
    /* the current master/using row has been written as matched */
    int m_matched = 0, u_matched = 0;
    while (mi < m1 && ui < u1) {
        int c = cm_cmp(J, &J->m, mi, &J->u, ui);
        if (c < 0) {
            if (!m_matched) cm_emit(&out[n++], (int32_t)mi, -1, MERGE_RESULT_MASTER_ONLY);
            mi++;
            m_matched = 0;
        } else if (c > 0) {
            if (!u_matched) cm_emit(&out[n++], -1, (int32_t)ui, MERGE_RESULT_USING_ONLY);
            ui++;
            u_matched = 0;
        } else if (J->type == MERGE_1_1) {
            cm_emit(&out[n++], (int32_t)mi++, (int32_t)ui++, MERGE_RESULT_BOTH);
        } else if (J->type == MERGE_M_1) {
            /* every master row of the group meets the one using row */
            cm_emit(&out[n++], (int32_t)mi++, (int32_t)ui, MERGE_RESULT_BOTH);
            u_matched = 1;
        } else if (J->type == MERGE_1_M) {
            cm_emit(&out[n++], (int32_t)mi, (int32_t)ui++, MERGE_RESULT_BOTH);
            m_matched = 1;
        } else {
            /* m:m pairs the group's rows in order; the shorter side repeats
             * its last row, as merge does */
            size_t me = mi + 1, ue = ui + 1;
            while (me < m1 && cm_cmp(J, &J->m, mi, &J->m, me) == 0) me++;
            while (ue < u1 && cm_cmp(J, &J->u, ui, &J->u, ue) == 0) ue++;
            size_t mc = me - mi, uc = ue - ui, pairs = mc > uc ? mc : uc;
            for (size_t i = 0; i < pairs; i++) {
                size_t a = i < mc ? i : mc - 1, b = i < uc ? i : uc - 1;
                cm_emit(&out[n++], (int32_t)(mi + a), (int32_t)(ui + b), MERGE_RESULT_BOTH);
            }
            mi = me;
            ui = ue;
        }
    }
    if (mi < m1 && m_matched) mi++;
    for (; mi < m1; mi++) cm_emit(&out[n++], (int32_t)mi, -1, MERGE_RESULT_MASTER_ONLY);
    if (ui < u1 && u_matched) ui++;
    for (; ui < u1; ui++) cm_emit(&out[n++], -1, (int32_t)ui, MERGE_RESULT_USING_ONLY);
    return n;
}

int64_t cmerge_sorted_join(
    stata_data *master_keys, stata_data *using_keys,
    int nkeys, cmerge_type_t merge_type,
    cmerge_output_spec_t **output_specs_out)
{
    const size_t m_nobs = master_keys->nobs;
    const size_t u_nobs = using_keys->nobs;

    *output_specs_out = NULL;
    for (int k = 0; k < nkeys; k++) {
        if (master_keys->vars[k].type != using_keys->vars[k].type) return -4;
    }
    if (nkeys < 1 || nkeys > CM_JOIN_MAX_KEYS) return -1;
    if (m_nobs > INT32_MAX || u_nobs > INT32_MAX) return -1;

    cm_join_t J;
    memset(&J, 0, sizeof(J));
    J.nkeys = nkeys;
    J.type = merge_type;
    for (int k = 0; k < nkeys; k++) {
        J.is_str[k] = master_keys->vars[k].type == STATA_TYPE_STRING;
        if (J.is_str[k]) {
            J.m.str[k] = master_keys->vars[k].data.str;
            J.u.str[k] = using_keys->vars[k].data.str;
        } else {
            J.m.num[k] = master_keys->vars[k].data.dbl;
            J.u.num[k] = using_keys->vars[k].data.dbl;
        }
    }
    J.single_num = nkeys == 1 && !J.is_str[0];

    size_t total = m_nobs + u_nobs;
    int chunks = 1;
#ifdef _OPENMP
    {
        int threads = ctools_get_max_threads();
        size_t want = total / CM_JOIN_ROWS_PER_CHUNK;
        if (want > (size_t)threads) want = (size_t)threads;
        if (want > CM_JOIN_MAX_CHUNKS) want = CM_JOIN_MAX_CHUNKS;
        if (want > 1) chunks = (int)want;
    }
#endif

    /* Validate every group, including unmatched keys, before emitting rows.
     * -2/-3 identify cardinality errors separately from allocation failure;
     * a nonunique master is reported before a nonunique using file. */
    if ((merge_type == MERGE_1_1 || merge_type == MERGE_1_M) && !cm_side_unique(&J, &J.m, m_nobs, chunks))
        return -2;
    if ((merge_type == MERGE_1_1 || merge_type == MERGE_M_1) && !cm_side_unique(&J, &J.u, u_nobs, chunks))
        return -3;

    /* Chunk boundaries: master rows that start a key group, and the using
     * rows where their keys begin. */
    size_t mb[CM_JOIN_MAX_CHUNKS + 1], ub[CM_JOIN_MAX_CHUNKS + 1];
    mb[0] = 0;
    ub[0] = 0;
    for (int c = 1; c < chunks; c++) {
        size_t b = m_nobs * (size_t)c / (size_t)chunks;
        if (b < mb[c - 1]) b = mb[c - 1];
        while (b > mb[c - 1] && b < m_nobs && cm_cmp(&J, &J.m, b - 1, &J.m, b) == 0) b++;
        mb[c] = b;
        ub[c] = b < m_nobs ? cm_using_lower_bound(&J, b, ub[c - 1], u_nobs) : u_nobs;
    }
    mb[chunks] = m_nobs;
    ub[chunks] = u_nobs;

    /* Each chunk writes at most its own rows on both sides; with one chunk
     * that region is the whole output. */
    size_t alloc_size;
    if (ctools_safe_mul_size(total ? total : 1, sizeof(cmerge_output_spec_t), &alloc_size) != 0)
        return -1;
    cmerge_output_spec_t *specs = malloc(alloc_size);
    if (!specs) return -1;
    if (chunks == 1) {
        size_t n = cm_join_chunk(&J, 0, m_nobs, 0, u_nobs, specs);
        *output_specs_out = specs;
        return (int64_t)n;
    }

    size_t count[CM_JOIN_MAX_CHUNKS], at[CM_JOIN_MAX_CHUNKS + 1];
    cmerge_output_spec_t *scratch = malloc(alloc_size);
    if (!scratch) {
        free(specs);
        return -1;
    }
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(chunks)
#endif
    for (int c = 0; c < chunks; c++) {
        count[c] = cm_join_chunk(&J, mb[c], mb[c + 1], ub[c], ub[c + 1],
                                 scratch + mb[c] + ub[c]);
    }
    at[0] = 0;
    for (int c = 0; c < chunks; c++) at[c + 1] = at[c] + count[c];
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(chunks)
#endif
    for (int c = 0; c < chunks; c++) {
        if (count[c])
            memcpy(specs + at[c], scratch + mb[c] + ub[c], count[c] * sizeof(cmerge_output_spec_t));
    }
    free(scratch);
    *output_specs_out = specs;
    return (int64_t)at[chunks];
}
