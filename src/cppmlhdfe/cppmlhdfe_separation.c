/*
 * cppmlhdfe_separation.c
 *
 * FE-based sample selection for PPML estimation: iteratively drops
 * singletons and FE groups whose outcome is all zero until stable.
 * Part of the ctools Stata plugin suite
 */

#include "cppmlhdfe_separation.h"
#include <string.h>

ST_int ppml_select_fe_sample(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated)
{
    return ppml_select_fe_sample_ex(y, levels, num_levels, G, N, mask,
        num_singletons, num_separated, NULL, 0, 1);
}

/* Rows a block removed in the current round, in row order */
typedef struct {
    ST_int *rows;
    size_t len, cap;
} ppml_removed;

static int ppml_removed_push(ppml_removed *list, ST_int row)
{
    if (list->len == list->cap) {
        size_t cap = list->cap ? 2 * list->cap : 1024;
        ST_int *grown = (ST_int *)malloc(cap * sizeof(ST_int));
        if (!grown) return -1;
        if (list->len) memcpy(grown, list->rows, list->len * sizeof(ST_int));
        free(list->rows);
        list->rows = grown;
        list->cap = cap;
    }
    list->rows[list->len++] = row;
    return 0;
}

/*
 * Each round removes, at once, every kept row that is alone in a level (a
 * singleton) or in a level with no positive outcome (separated), judged on
 * the rows kept when the round starts; rounds repeat until none is removed.
 * The per-level counts of kept rows and of kept positive outcomes are built
 * once, one factor per thread, and each round subtracts the rows it removed
 * instead of recounting; the rows are judged in parallel row blocks. The
 * rounds, and so the removals and their classification, are those of a
 * recount every round.
 */
ST_int ppml_select_fe_sample_ex(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated,
    const ST_int *slope_only, ST_int keep_singletons, ST_int check_separation)
{
    ST_int **counts = calloc((size_t)(G > 0 ? G : 1), sizeof(*counts));
    ST_int **positive = calloc((size_t)(G > 0 ? G : 1), sizeof(*positive));
    const int B = ctools_hdfe_row_blocks(N);
    ppml_removed *removed = calloc((size_t)B, sizeof(*removed));
    ST_int rc = -1, g;
    int bad = 0, b;
    *num_singletons = *num_separated = 0;
    if (!counts || !positive || !removed) goto cleanup;
    for (g = 0; g < G; g++) {
        if (num_levels[g] < 1) goto cleanup;
        counts[g] = calloc((size_t)num_levels[g], sizeof(**counts));
        positive[g] = calloc((size_t)num_levels[g], sizeof(**positive));
        if (!counts[g] || !positive[g]) goto cleanup;
    }

    /* Every row's levels must be valid; the kept rows are counted */
    #pragma omp parallel for schedule(dynamic) reduction(| : bad) if (B > 1)
    for (g = 0; g < G; g++) {
        const ST_int *lev = levels[g];
        const ST_int L = num_levels[g];
        ST_int *c = counts[g], *p = positive[g];
        for (ST_int i = 0; i < N; i++) {
            if (lev[i] < 1 || lev[i] > L) {
                bad = 1;
                continue;
            }
            if (!mask[i]) continue;
            c[lev[i] - 1]++;
            if (y[i] > 0) p[lev[i] - 1]++;
        }
    }
    if (bad) goto cleanup;

    for (;;) {
        ST_int round_singletons = 0, round_separated = 0;
        int failed = 0;
        #pragma omp parallel for schedule(static) if (B > 1) \
            reduction(+ : round_singletons, round_separated) reduction(| : failed)
        for (b = 0; b < B; b++) {
            ST_int i0, i1;
            ctools_hdfe_block_range(N, B, b, &i0, &i1);
            ppml_removed *list = &removed[b];
            list->len = 0;
            for (ST_int i = i0; i < i1; i++) {
                if (!mask[i]) continue;
                ST_int singleton = 0, separated = 0;
                for (ST_int h = 0; h < G; h++) {
                    ST_int level = levels[h][i] - 1;
                    singleton |= !keep_singletons && counts[h][level] == 1;
                    separated |= check_separation && !(slope_only && slope_only[h]) &&
                                 positive[h][level] == 0;
                }
                if (singleton || separated) {
                    if (ppml_removed_push(list, i)) {
                        failed = 1;
                        break;
                    }
                    mask[i] = 0;
                    if (singleton) round_singletons++;
                    else round_separated++;
                }
            }
        }
        if (failed) goto cleanup;
        *num_singletons += round_singletons;
        *num_separated += round_separated;
        if (!round_singletons && !round_separated) break;

        /* The next round's counts: the removed rows leave their levels */
        #pragma omp parallel for schedule(dynamic) if (B > 1)
        for (g = 0; g < G; g++) {
            const ST_int *lev = levels[g];
            ST_int *c = counts[g], *p = positive[g];
            for (int t = 0; t < B; t++) {
                const ppml_removed *list = &removed[t];
                for (size_t k = 0; k < list->len; k++) {
                    const ST_int i = list->rows[k];
                    c[lev[i] - 1]--;
                    if (y[i] > 0) p[lev[i] - 1]--;
                }
            }
        }
    }
    rc = 0;
cleanup:
    for (g = 0; g < G; g++) {
        if (counts) free(counts[g]);
        if (positive) free(positive[g]);
    }
    for (b = 0; removed && b < B; b++) free(removed[b].rows);
    free(counts); free(positive); free(removed);
    return rc;
}
