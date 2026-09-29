/* Inclusive interval join. Sort the using index once; two binary searches per
 * master row. Output size is checked before Stata's dataset is resized. */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "crangejoin_impl.h"
#include "ctools_types.h"
#include "ctools_order.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"

static struct {
    ctools_filtered_data using, master;
    int *ukeys, *mkeys;
    size_t *lo, *hi, *offset, *groups;
    double *sorted_key;
    size_t ngroups;
    int nby, key, nm, verbose, ready;
    size_t nout;
} state;

void crangejoin_cleanup_cache(void)
{
    ctools_filtered_data_free(&state.using);
    ctools_filtered_data_free(&state.master);
    free(state.ukeys); free(state.mkeys);
    free(state.lo); free(state.hi); free(state.offset);
    free(state.groups); free(state.sorted_key);
    memset(&state, 0, sizeof(state));
}

static int compare_group(size_t m, size_t u)
{
    return ctools_compare_rows(&state.master.data, m, state.mkeys,
                              &state.using.data, u, state.ukeys, state.nby);
}

static size_t bound(double value, int upper, size_t l, size_t r)
{
    while (l < r) {
        size_t mid = l + (r-l)/2;
        double x = state.sorted_key[mid];
        if (upper ? x <= value : x < value) l = mid+1;
        else r = mid;
    }
    return l;
}

static size_t find_group(size_t m)
{
    size_t l = 0, r = state.ngroups;
    while (l < r) {
        size_t mid = l + (r-l)/2;
        if (compare_group(m, state.using.data.sort_order[state.groups[mid]]) > 0) l = mid+1;
        else r = mid;
    }
    if (l < state.ngroups && !compare_group(m, state.using.data.sort_order[state.groups[l]])) return l;
    return state.ngroups;
}

ST_retcode crangejoin_main(const char *args)
{
    int rc = 0;
    if (!strcmp(args, "clear")) { crangejoin_cleanup_cache(); return 0; }
    if (!strncmp(args, "using ", 6)) {
        crangejoin_cleanup_cache(); args += 6;
        if (ctools_parse_next_int(&args, &state.nby) || state.nby < 0 || state.nby >= SF_nvars() ||
            ctools_parse_next_int(&args, &state.key) || state.key <= state.nby || state.key > SF_nvars() ||
            ctools_parse_next_int(&args, &state.verbose)) { rc = 198; goto fail; }
        state.key--;
        state.ukeys = malloc((state.nby+1) * sizeof(*state.ukeys));
        if (!state.ukeys) { rc = 920; goto fail; }
        for (int k = 0; k < state.nby; k++) state.ukeys[k] = k;
        state.ukeys[state.nby] = state.key;
        double start = ctools_timer_seconds();
        if (ctools_data_load(&state.using, NULL, 0, 0, 0, CTOOLS_LOAD_SKIP_IF)) { rc = 920; goto fail; }
        if (!state.using.data.nobs) { rc = 2000; goto fail; }
        if (state.using.data.vars[state.key].type != STATA_TYPE_DOUBLE) { rc = 109; goto fail; }
        for (size_t i = 0; i < state.using.data.nobs; i++)
            if (state.using.data.vars[state.key].data.dbl[i] >= SV_missval) { rc = 198; goto fail; }
        double loaded = ctools_timer_seconds();
        if (ctools_order_stable(&state.using.data, state.ukeys, state.nby+1)) { rc = 920; goto fail; }
        double sorted = ctools_timer_seconds();
        size_t n = state.using.data.nobs;
        state.sorted_key = malloc(n * sizeof(*state.sorted_key));
        state.groups = malloc((state.nby ? n+1 : 2) * sizeof(*state.groups));
        if (!state.sorted_key || !state.groups) { rc = 920; goto fail; }
        const perm_idx_t *p = state.using.data.sort_order;
        const double *x = state.using.data.vars[state.key].data.dbl;
        #pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if(n > 10000)
        for (size_t i = 0; i < n; i++) state.sorted_key[i] = x[p[i]];
        state.groups[0] = 0; state.ngroups = 1;
        if (state.nby) for (size_t i = 1; i < n; i++)
            if (ctools_compare_rows(&state.using.data, p[i-1], state.ukeys,
                                    &state.using.data, p[i], state.ukeys, state.nby))
                state.groups[state.ngroups++] = i;
        state.groups[state.ngroups] = n;
        state.ready = 1;
        ctools_verbose("crangejoin", state.verbose, "using load %.4fs; index sort %.4fs; group/key index %.4fs",
            loaded-start, sorted-loaded, ctools_timer_seconds()-sorted);
        return 0;
    }
    if (!strncmp(args, "prepare ", 8)) {
        args += 8;
        if (state.ready != 1 || ctools_parse_next_int(&args, &state.nm) || state.nm < 1 ||
            SF_nvars() != state.nm + 2 + (int)state.using.data.nvars - state.nby) { rc = 198; goto fail; }
        state.mkeys = malloc((state.nby+1) * sizeof(*state.mkeys));
        int *vars = malloc((state.nm+2) * sizeof(*vars));
        if (!state.mkeys || !vars) { free(vars); rc = 920; goto fail; }
        for (int k = 0; k < state.nby; k++) {
            if (ctools_parse_next_int(&args, &state.mkeys[k]) || state.mkeys[k] < 1 || state.mkeys[k] > state.nm) {
                free(vars); rc = 198; goto fail;
            }
            state.mkeys[k]--;
        }
        for (int k = 0; k < state.nm+2; k++) vars[k] = k+1;
        double start = ctools_timer_seconds();
        rc = ctools_data_load(&state.master, vars, state.nm+2, 0, 0, CTOOLS_LOAD_SKIP_IF);
        free(vars);
        if (rc) { rc = 920; goto fail; }
        for (int k = 0; k < state.nby; k++)
            if (state.master.data.vars[state.mkeys[k]].type != state.using.data.vars[k].type) { rc = 106; goto fail; }
        if (state.master.data.vars[state.nm].type != STATA_TYPE_DOUBLE || state.master.data.vars[state.nm+1].type != STATA_TYPE_DOUBLE) { rc = 109; goto fail; }
        size_t n = state.master.data.nobs;
        if (!n) { rc = 2000; goto fail; }
        state.lo = malloc(n * sizeof(*state.lo)); state.hi = malloc(n * sizeof(*state.hi));
        state.offset = malloc((n+1) * sizeof(*state.offset));
        if (!state.lo || !state.hi || !state.offset) { rc = 920; goto fail; }
        double loaded = ctools_timer_seconds();
        const double *low = state.master.data.vars[state.nm].data.dbl;
        const double *high = state.master.data.vars[state.nm+1].data.dbl;
        #pragma omp parallel num_threads(ctools_get_max_threads()) if(n > 10000)
        {
            size_t previous = n, group = state.ngroups;
            #pragma omp for schedule(static)
            for (size_t i = 0; i < n; i++) {
                if (low[i] > high[i]) { state.lo[i] = state.hi[i] = 0; continue; }
                if (!state.nby) group = 0;
                else if (previous == n || ctools_compare_rows(&state.master.data, previous, state.mkeys,
                                                              &state.master.data, i, state.mkeys, state.nby))
                    group = find_group(i);
                previous = i;
                if (group == state.ngroups) state.lo[i] = state.hi[i] = 0;
                else {
                    size_t end = state.groups[group+1];
                    state.lo[i] = bound(low[i], 0, state.groups[group], end);
                    state.hi[i] = bound(high[i], 1, state.lo[i], end);
                }
            }
        }
        state.offset[0] = 0;
        for (size_t i = 0; i < n; i++) {
            size_t count = state.hi[i] - state.lo[i];
            if (!count) count = 1;
            if (count > INT_MAX - state.offset[i]) {
                ctools_error("crangejoin", "joined output exceeds the Stata plugin observation limit");
                rc = 920; goto fail;
            }
            state.offset[i+1] = state.offset[i] + count;
        }
        state.nout = state.offset[n];
        char result[32]; snprintf(result, sizeof(result), "%zu", state.nout);
        rc = SF_macro_save("_crangejoin_nout", result);
        if (rc) goto fail;
        state.ready = 2;
        ctools_verbose("crangejoin", state.verbose, "master load %.4fs; match/count %.4fs; %zu output rows",
            loaded-start, ctools_timer_seconds()-loaded, state.nout);
        return 0;
    }
    if (!strcmp(args, "write")) {
        int nu = (int)state.using.data.nvars - state.nby;
        if (state.ready != 2 || (size_t)SF_nobs() != state.nout || SF_nvars() != state.nm+2+nu) { rc = 198; goto fail; }
        for (int col = 0; col < state.nm+nu; col++) {
            int ismaster = col < state.nm;
            int dest = ismaster ? col+1 : col+3;
            const stata_variable *v = ismaster ? state.master.data.vars+col : state.using.data.vars+state.nby+col-state.nm;
            if (SF_var_is_string(dest) != (v->type == STATA_TYPE_STRING) || SF_var_is_strl(dest)) { rc = 109; goto fail; }
        }
        double start = ctools_timer_seconds();
        /* Tile by output rows, not columns: narrow results and a single
         * master row with millions of matches can both use every worker.
         * Workers own disjoint row ranges, avoiding cross-column false sharing. */
        const size_t tile = 4096, nm = state.master.data.nobs;
        #pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if(state.nout > 10000) reduction(max:rc)
        for (size_t begin = 0; begin < state.nout; begin += tile) {
            size_t end = begin + tile < state.nout ? begin + tile : state.nout;
            size_t l = 0, r = nm;
            while (l < r) {
                size_t mid = l + (r-l)/2;
                if (state.offset[mid+1] <= begin) l = mid+1;
                else r = mid;
            }
            for (int col = 0; col < state.nm+nu; col++) {
                int ismaster = col < state.nm;
                int dest = ismaster ? col+1 : col+3;
                const stata_variable *v = ismaster ? state.master.data.vars+col : state.using.data.vars+state.nby+col-state.nm;
                for (size_t i = l; i < nm && state.offset[i] < end; i++) {
                    int matched = state.hi[i] > state.lo[i];
                    size_t first = state.offset[i] > begin ? state.offset[i] : begin;
                    size_t stop = state.offset[i+1] < end ? state.offset[i+1] : end;
                    for (size_t j = first; j < stop; j++) {
                        size_t source = ismaster ? i : (matched ? state.using.data.sort_order[state.lo[i]+j-state.offset[i]] : 0);
                        int err;
                        if (v->type == STATA_TYPE_STRING)
                            err = SF_sstore(dest, (ST_int)(j+1), (ismaster || matched) ? v->data.str[source] : "");
                        else
                            err = SF_vstore(dest, (ST_int)(j+1), (ismaster || matched) ? v->data.dbl[source] : SV_missval);
                        if (err > rc) rc = err;
                    }
                }
            }
        }
        ctools_verbose("crangejoin", state.verbose, "write %.4fs", ctools_timer_seconds()-start);
        crangejoin_cleanup_cache();
        return rc;
    }
    rc = 198;
fail:
    crangejoin_cleanup_cache();
    return rc;
}
