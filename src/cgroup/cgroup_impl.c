/*
    cgroup_impl.c
    Group identifiers in C, replacing egen group() in the ado files.

    Keys are loaded through the shared transport, ordered by the stable sorts
    csort uses, and numbered in one pass over that order, so the numbering is
    egen group()'s: tuples in sort order, with observations that have any
    missing key left missing.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_parse.h"
#include "cgroup_impl.h"

/* egen group(..., label) defines at most this many label values. */
#define CGROUP_MAX_LABELS 65536

static int cgroup_key_missing(const stata_data *data, size_t nkeys, size_t row)
{
    for (size_t k = 0; k < nkeys; k++) {
        const stata_variable *v = &data->vars[k];
        if (v->type == STATA_TYPE_DOUBLE) {
            if (SF_is_missing(v->data.dbl[row])) return 1;
        } else if (!v->data.str[row] || v->data.str[row][0] == '\0') {
            return 1;
        }
    }
    return 0;
}

static int cgroup_same_keys(const stata_data *data, size_t nkeys, size_t a, size_t b)
{
    for (size_t k = 0; k < nkeys; k++) {
        const stata_variable *v = &data->vars[k];
        if (v->type == STATA_TYPE_DOUBLE) {
            if (v->data.dbl[a] != v->data.dbl[b]) return 0;
        } else if (strcmp(v->data.str[a], v->data.str[b]) != 0) {
            return 0;
        }
    }
    return 1;
}

ST_retcode cgroup_main(const char *args)
{
    int nkeys = 0;
    if (!args || ctools_parse_int_option(args, "nkeys", &nkeys) != CTOOLS_PARSE_OK || nkeys < 1) {
        SF_error("group: invalid arguments\n");
        return 198;
    }
    const int out = ctools_parse_bool_option(args, "out");
    const int want_reps = ctools_parse_bool_option(args, "reps");
    /* Without an if condition (skipif), or with a sample marker after the
     * other variables in its place (touse), rows are selected without one
     * SF_ifobs() call per observation. */
    const int touse = ctools_parse_bool_option(args, "touse");
    const int load_flags = ctools_parse_bool_option(args, "skipif") ? CTOOLS_LOAD_SKIP_IF
                         : touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF;
    const ST_int dest = (ST_int)nkeys + 1;
    if (SF_nvars() != nkeys + (out ? 1 : 0) + (touse ? 1 : 0) || (out && SF_var_is_string(dest))) {
        SF_error("group: expected the key variables and a numeric destination\n");
        return 198;
    }

    int *key_vars = (int *)malloc((size_t)nkeys * sizeof(int));
    if (!key_vars) return 920;
    for (int k = 0; k < nkeys; k++) key_vars[k] = k + 1;

    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    stata_retcode load_rc = ctools_data_load(&fd, key_vars, (size_t)nkeys, 0, 0, load_flags);
    if (load_rc != STATA_OK) {
        free(key_vars);
        ctools_filtered_data_free(&fd);
        if (load_rc == STATA_ERR_UNSUPPORTED_TYPE) return 109;
        return ctools_stata_rc(load_rc);
    }

    const size_t nobs = fd.data.nobs;
    stata_retcode rc = STATA_OK;
    if (nobs > 1) rc = ctools_sort_dispatch(&fd.data, key_vars, (size_t)nkeys, SORT_ALG_AUTO);
    free(key_vars);

    /* ids[i]: group of loaded row i, 0 when a key is missing. */
    uint32_t *ids = nobs ? (uint32_t *)calloc(nobs, sizeof(uint32_t)) : NULL;
    size_t reps_cap = 0, groups = 0, previous = 0;
    perm_idx_t *reps = NULL;
    int have_previous = 0;
    if (rc == STATA_OK && nobs && !ids) rc = STATA_ERR_MEMORY;
    for (size_t pos = 0; rc == STATA_OK && pos < nobs; pos++) {
        size_t row = nobs > 1 ? fd.data.sort_order[pos] : pos;
        if (cgroup_key_missing(&fd.data, (size_t)nkeys, row)) continue;
        if (!have_previous || !cgroup_same_keys(&fd.data, (size_t)nkeys, previous, row)) {
            groups++;
            if (want_reps) {
                if (groups > CGROUP_MAX_LABELS) {
                    SF_error("too many values for label define\n");
                    ctools_filtered_data_free(&fd);
                    free(ids);
                    free(reps);
                    return 134;
                }
                if (groups > reps_cap) {
                    size_t cap = reps_cap ? reps_cap * 2 : 64;
                    perm_idx_t *grown = (perm_idx_t *)realloc(reps, cap * sizeof(perm_idx_t));
                    if (!grown) {
                        rc = STATA_ERR_MEMORY;
                        break;
                    }
                    reps = grown;
                    reps_cap = cap;
                }
                reps[groups - 1] = fd.obs_map[row];
            }
        }
        ids[row] = (uint32_t)groups;
        previous = row;
        have_previous = 1;
    }

    /* Every observation of the in-range: its group, or missing. */
    if (rc == STATA_OK && out) {
        ST_int obs1 = SF_in1(), obs2 = SF_in2();
        size_t nrange = obs2 >= obs1 ? (size_t)(obs2 - obs1) + 1 : 0;
        double *values = nrange ? (double *)malloc(nrange * sizeof(double)) : NULL;
        if (nrange && !values) rc = STATA_ERR_MEMORY;
        if (rc == STATA_OK) {
            for (size_t i = 0; i < nrange; i++) values[i] = SV_missval;
            for (size_t i = 0; i < nobs; i++)
                if (ids[i]) values[fd.obs_map[i] - (perm_idx_t)obs1] = (double)ids[i];
            int store_error = 0;
            const ST_IIID store_fn = (_stata_)->safestore;
            #pragma omp parallel for schedule(static) if (nrange >= MIN_OBS_PER_THREAD * 2)
            for (size_t i = 0; i < nrange; i++) {
                int io_rc = store_fn(dest, obs1 + (ST_int)i, values[i]);
                if (io_rc) {
                    #pragma omp atomic write
                    store_error = io_rc;
                }
            }
            if (store_error) rc = store_error == 1 ? STATA_ERR_CANCELLED : STATA_ERR_STATA_WRITE;
        }
        free(values);
    }

    if (rc == STATA_OK && want_reps) {
        /* Up to 11 characters per observation number. */
        char *list = (char *)malloc(groups * 12 + 1);
        if (!list) rc = STATA_ERR_MEMORY;
        else {
            size_t used = 0;
            list[0] = '\0';
            for (size_t g = 0; g < groups; g++)
                used += (size_t)snprintf(list + used, 13, g ? " %u" : "%u", (unsigned)reps[g]);
            if (SF_macro_save("_ctools_group_reps", list) != 0) rc = STATA_ERR_MEMORY;
            free(list);
        }
    }

    ctools_filtered_data_free(&fd);
    free(ids);
    free(reps);
    if (rc != STATA_OK) {
        if (rc == STATA_ERR_STATA_WRITE) SF_error("group: failed to write Stata data\n");
        return ctools_stata_rc(rc);
    }
    SF_scal_save("_ctools_group_n", (double)groups);
    return 0;
}
