#define SD_SAFEMODE
#include "hdfe_missing.h"
#include "../ctools_config.h"
#include <stdlib.h>
#include <string.h>

/* A column slice fits in L1/L2. Workers own disjoint row blocks, so the
 * screening uses all available threads even for a model with few columns. */
#define HDFE_SCREEN_BLOCK 4096

ST_retcode hdfe_markout_main(const char *args)
{
    (void)args;
    const ST_int marker = SF_nvars();
    if (marker < 2) return 198;
    for (ST_int j = 1; j <= marker; j++)
        if (SF_var_is_string(j)) return 109;
    const ST_int first = SF_in1(), last = SF_in2();
    if (last < first) return 0;
    if (first < 1 || last > SF_nobs()) return 198;
    const size_t n = (size_t)last - first + 1;
    const size_t blocks = 1 + (n - 1) / HDFE_SCREEN_BLOCK;
    int threads = ctools_get_max_threads();
    if (threads < 1) threads = 1;
    if ((size_t)threads > blocks) threads = (int)blocks;
    const int parallel = n >= 32768 && threads > 1;
#ifndef _OPENMP
    (void)parallel;
#endif
    /* Match the numeric loader: validated indices/ranges allow unchecked
     * reads. Writes always use the checked, private-variable-safe callback. */
    ST_IIIDp read_value = (_stata_)->vdata;
    ST_IIID store_marker = (_stata_)->safestore;
    const double missing = SV_missval;
    int error = 0;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel) reduction(max:error)
    for (size_t b = 0; b < blocks; b++) {
        double values[HDFE_SCREEN_BLOCK];
        unsigned char keep[HDFE_SCREEN_BLOCK], selected[HDFE_SCREEN_BLOCK];
        const size_t offset = b * HDFE_SCREEN_BLOCK;
        const size_t count = n - offset < HDFE_SCREEN_BLOCK ? n - offset : HDFE_SCREEN_BLOCK;
        const ST_int row = first + (ST_int)offset;
        int rc = 0, active = 0;
        for (size_t i = 0; i < count; i++) {
            double v;
            rc = read_value(marker, row + (ST_int)i, &v);
            if (rc) break;
            selected[i] = keep[i] = (v != 0 && v < missing);
            active |= keep[i];
        }
        for (ST_int j = 1; j < marker && !rc && active; j++) {
            for (size_t i = 0; i < count; i++) {
                values[i] = 0;
                if (keep[i]) {
                    rc = read_value(j, row + (ST_int)i, values + i);
                    if (rc) break;
                }
            }
            if (rc) break;
            active = 0;
#pragma omp simd reduction(|:active)
            for (size_t i = 0; i < count; i++) {
                keep[i] &= values[i] < missing;
                active |= keep[i];
            }
        }
        if (!rc) {
            for (size_t i = 0; i < count; i++) if (selected[i] && !keep[i]) {
                int write_rc = store_marker(marker, row + (ST_int)i, 0.0);
                if (write_rc > rc) rc = write_rc;
            }
        }
        if (rc > error) error = rc;
    }
    return error;
}

ST_retcode hdfe_screen_missing(ctools_filtered_data *fd, ST_int marker)
{
    if (!fd || marker < 1 || fd->data.sort_order) return 198;
    const size_t n = fd->data.nobs, k = fd->data.nvars;
    if (!n) return 0;
    if (!k || !fd->data.vars || !fd->obs_map) return 198;
    for (size_t j = 0; j < k; j++)
        if (fd->data.vars[j].type != STATA_TYPE_DOUBLE ||
            fd->data.vars[j].nobs != n || !fd->data.vars[j].data.dbl) return 198;

    /* The load compared every value it read with the missing values: a
     * column it found without one is not read again, and when that is every
     * column there is nothing to do. Columns it did not note are read. */
    size_t n_check = 0;
    for (size_t j = 0; j < k; j++)
        n_check += fd->data.vars[j].missing_state != CTOOLS_MISSING_NONE;
    if (!n_check) return 0;

    const size_t blocks = 1 + (n - 1) / HDFE_SCREEN_BLOCK;
    unsigned char *bad = calloc(blocks, 1);
    if (!bad) return 920;
    int threads = ctools_get_max_threads();
    if (threads < 1) threads = 1;
    if ((size_t)threads > blocks) threads = (int)blocks;
    const int parallel = n >= 32768 && threads > 1;
#ifndef _OPENMP
    (void)parallel;
#endif
    const double missing = SV_missval;
    size_t bad_blocks = 0;

    /* Common case: just vectorized reads and one flag per block. No shared
     * row-mask writes, atomics, per-value SPI calls, or data copies. */
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel) reduction(+:bad_blocks)
    for (size_t b = 0; b < blocks; b++) {
        const size_t first = b * HDFE_SCREEN_BLOCK;
        const size_t last = n - first < HDFE_SCREEN_BLOCK ? n : first + HDFE_SCREEN_BLOCK;
        int any = 0;
        for (size_t j = 0; j < k && !any; j++) {
            if (fd->data.vars[j].missing_state == CTOOLS_MISSING_NONE) continue;
            const double *v = fd->data.vars[j].data.dbl;
#pragma omp simd reduction(|:any)
            for (size_t i = first; i < last; i++) any |= !(v[i] < missing);
        }
        bad[b] = (unsigned char)any;
        bad_blocks += any;
    }
    if (!bad_blocks) { free(bad); return 0; }

    /* Only blocks with missing values need row flags. Good blocks are
     * compacted in bulk, so sparse missingness does not add K branchy scans. */
    unsigned char *keep = malloc(n);
    if (!keep) { free(bad); return 920; }
    size_t dropped = 0;
    int write_rc = 0;
    /* Resolve plugin-visible indices, including private Stata variables.
     * Cache the callback outside the workers and independently of callers'
     * SD_FASTMODE preprocessor state. */
    ST_IIID store_marker = (_stata_)->safestore;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel) reduction(+:dropped) reduction(max:write_rc)
    for (size_t b = 0; b < blocks; b++) {
        if (!bad[b]) continue;
        const size_t first = b * HDFE_SCREEN_BLOCK;
        const size_t last = n - first < HDFE_SCREEN_BLOCK ? n : first + HDFE_SCREEN_BLOCK;
        memset(keep + first, 1, last - first);
        for (size_t j = 0; j < k; j++) {
            if (fd->data.vars[j].missing_state == CTOOLS_MISSING_NONE) continue;
            const double *v = fd->data.vars[j].data.dbl;
#pragma omp simd
            for (size_t i = first; i < last; i++) keep[i] &= (v[i] < missing);
        }
        for (size_t i = first; i < last; i++) {
            if (!keep[i]) {
                dropped++;
                int rc = store_marker(marker, (ST_int)fd->obs_map[i], 0.0);
                if (rc > write_rc) write_rc = rc;
            }
        }
    }
    if (write_rc) { free(keep); free(bad); return write_rc; }

    /* In-place compaction is parallel across columns, never across row
     * blocks: a later block must not overwrite an earlier worker's input.
     * The observation map is an independent final task. */
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
    for (size_t j = 0; j <= k; j++) {
        size_t out = 0;
        for (size_t b = 0; b < blocks; b++) {
            const size_t first = b * HDFE_SCREEN_BLOCK;
            const size_t last = n - first < HDFE_SCREEN_BLOCK ? n : first + HDFE_SCREEN_BLOCK;
            if (!bad[b]) {
                if (out != first) {
                    if (j == k) memmove(fd->obs_map + out, fd->obs_map + first, (last - first) * sizeof(*fd->obs_map));
                    else memmove(fd->data.vars[j].data.dbl + out, fd->data.vars[j].data.dbl + first, (last - first) * sizeof(double));
                }
                out += last - first;
            } else {
                for (size_t i = first; i < last; i++) if (keep[i]) {
                    if (j == k) fd->obs_map[out] = fd->obs_map[i];
                    else fd->data.vars[j].data.dbl[out] = fd->data.vars[j].data.dbl[i];
                    out++;
                }
            }
        }
        if (j < k) fd->data.vars[j].nobs = n - dropped;
    }
    fd->data.nobs = n - dropped;
    fd->was_filtered = 1;
    free(keep);
    free(bad);
    return 0;
}
