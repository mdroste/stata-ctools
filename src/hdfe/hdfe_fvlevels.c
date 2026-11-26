/*
 * hdfe_fvlevels.c - the levels of factor variables on a sample, for the
 * wrappers that build i.f indicators in the plugin instead of fvrevar
 *
 * plugin call ctools_plugin f1 ... fk touse, "fvlevels"
 *
 * For each fj, the distinct values on the rows where touse is nonzero, in
 * ascending order, as the local macro ctools_fvl_<j> ("1 2 5"), and
 * ctools_fvl_ok_<j>: 1 when every value there is an integer from 0 to
 * CTOOLS_FVL_MAX (missing values are skipped, as fvexpand skips them), 0
 * otherwise; the wrapper then takes Stata's own fvrevar/fvexpand path,
 * which also gives Stata's messages for noninteger or negative values.
 */
#define SD_SAFEMODE
#include "hdfe_missing.h"
#include "../ctools_config.h"
#include "../ctools_threads.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define CTOOLS_FVL_MAX 65535
#define CTOOLS_FVL_BLOCK 4096

ST_retcode hdfe_fvlevels_main(const char *args)
{
    (void)args;
    const ST_int marker = SF_nvars();
    if (marker < 2) return 198;
    for (ST_int j = 1; j <= marker; j++)
        if (SF_var_is_string(j)) return 109;
    const ST_int first = SF_in1(), last = SF_in2();
    if (first < 1 || last > SF_nobs()) return 198;
    const size_t n = last >= first ? (size_t)last - first + 1 : 0;
    const size_t blocks = n ? 1 + (n - 1) / CTOOLS_FVL_BLOCK : 0;
    const ST_int k = marker - 1;
    const size_t span = (size_t)CTOOLS_FVL_MAX + 1;
    int threads = ctools_get_max_threads();
    if (threads < 1) threads = 1;
    if (blocks && (size_t)threads > blocks) threads = (int)blocks;
    if (!blocks) threads = 1;
    /* one presence table per thread and variable (64 KB each), OR'ed after */
    unsigned char *seen = (unsigned char *)calloc((size_t)threads * (size_t)k, span);
    int *bad = (int *)calloc((size_t)k, sizeof(int));
    if (!seen || !bad) {
        free(seen); free(bad);
        return 920;
    }
    ST_IIIDp read_value = (_stata_)->vdata;
    const double missing = SV_missval;
    int error = 0;
#pragma omp parallel num_threads(threads) if(threads > 1 && n >= 32768) reduction(max:error)
    {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        unsigned char keep[CTOOLS_FVL_BLOCK];
        int my_bad[64] = {0};   /* k is checked below */
        int *badv = k <= 64 ? my_bad : NULL;
#pragma omp for schedule(static)
        for (size_t b = 0; b < blocks; b++) {
            const size_t offset = b * CTOOLS_FVL_BLOCK;
            const size_t count = n - offset < CTOOLS_FVL_BLOCK ? n - offset : CTOOLS_FVL_BLOCK;
            const ST_int row = first + (ST_int)offset;
            int rc = 0, active = 0;
            for (size_t i = 0; i < count && !rc; i++) {
                double v;
                rc = read_value(marker, row + (ST_int)i, &v);
                keep[i] = !rc && v != 0 && v < missing;
                active |= keep[i];
            }
            for (ST_int j = 0; j < k && !rc && active; j++) {
                unsigned char *tab = seen + ((size_t)t * (size_t)k + (size_t)j) * span;
                for (size_t i = 0; i < count; i++) {
                    if (!keep[i]) continue;
                    double v;
                    rc = read_value(j + 1, row + (ST_int)i, &v);
                    if (rc) break;
                    if (!(v < missing)) continue;
                    if (v >= 0 && v <= CTOOLS_FVL_MAX && v == floor(v)) tab[(size_t)v] = 1;
                    else if (badv) badv[j] = 1;
                    else {
#pragma omp atomic write
                        bad[j] = 1;
                    }
                }
            }
            if (rc > error) error = rc;
        }
        if (badv)
            for (ST_int j = 0; j < k; j++)
                if (badv[j]) {
#pragma omp atomic write
                    bad[j] = 1;
                }
    }
    ST_retcode rc = error;
    char name[64];
    for (ST_int j = 0; j < k && !rc; j++) {
        /* OR the thread tables into the first */
        unsigned char *tab0 = seen + (size_t)j * span;
        for (int t = 1; t < threads; t++) {
            const unsigned char *tab = seen + ((size_t)t * (size_t)k + (size_t)j) * span;
            for (size_t v = 0; v < span; v++) tab0[v] |= tab[v];
        }
        size_t nlev = 0;
        for (size_t v = 0; v < span; v++) nlev += tab0[v];
        char *list = (char *)malloc(nlev * 7 + 1);
        if (!list) { rc = 920; break; }
        size_t len = 0;
        list[0] = '\0';
        for (size_t v = 0; v < span; v++)
            if (tab0[v]) len += (size_t)snprintf(list + len, 8, len ? " %u" : "%u", (unsigned)v);
        snprintf(name, sizeof(name), "_ctools_fvl_%d", (int)j + 1);
        rc = SF_macro_save(name, list);
        free(list);
        if (rc) break;
        snprintf(name, sizeof(name), "_ctools_fvl_ok_%d", (int)j + 1);
        rc = SF_macro_save(name, (bad[j] || !nlev) ? "0" : "1");
    }
    free(seen);
    free(bad);
    return rc;
}
