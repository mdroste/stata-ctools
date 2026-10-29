/* Isolated REAL-SPI transport benchmark. Build separately from the dispatcher.
 * Arguments: label threads mode verify; mode=identity|sorted|scattered|filtered|read.
 * Timed regions exclude verification, permutation, and input restoration.
 * Verification rereads every cell; Stata also checks the restored datasignature.
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_threads.h"

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static int check(const stata_data *data, const perm_idx_t *map,
                 const perm_idx_t *order, size_t first)
{
    for (size_t j = 0; j < data->nvars; j++) {
        const stata_variable *v = data->vars + j;
        for (size_t i = 0; i < data->nobs; i++) {
            size_t row = order ? order[i] : i;
            ST_int obs = map ? (ST_int)map[i] : (ST_int)(first + i);
            if (v->type == STATA_TYPE_DOUBLE) {
                double value;
                if ((_stata_)->safevdata((ST_int)j + 1, obs, &value) ||
                    memcmp(&value, v->data.dbl + row, sizeof(value))) return 9;
            } else {
                char value[2046];
                if (SF_var_is_strl((ST_int)j + 1)) {
                    int len = SF_sdatalen((ST_int)j + 1, obs);
                    if (len < 0 || len > 2045 || SF_strldata((ST_int)j + 1, obs, value, len + 1) != len) return 9;
                    value[len] = 0;
                } else if (SF_sdata((ST_int)j + 1, obs, value)) return 9;
                if (strcmp(value, v->data.str[row] ? v->data.str[row] : "")) return 9;
            }
        }
    }
    return 0;
}

STDLL stata_call(int argc, char *argv[])
{
    if (argc != 4) return 198;
    /* Match ctools_plugin.c's initialization when using its static runtime. */
    setenv("KMP_DUPLICATE_LIB_OK", "TRUE", 0);
    ctools_set_max_threads(atoi(argv[1]));
    int filtered = !strcmp(argv[2], "filtered");
    int filtered_read = !strcmp(argv[2], "filteredread");
    int check_if = filtered || filtered_read || !strcmp(argv[2], "checked");
    int sorted = !strcmp(argv[2], "sorted");
    int scattered = !strcmp(argv[2], "scattered");
    int read_only = filtered_read || !strcmp(argv[2], "read");
    int verify = atoi(argv[3]);
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    double start = now();
    int rc = ctools_data_load(&fd, NULL, 0, 0, 0,
                             check_if ? CTOOLS_LOAD_CHECK_IF : CTOOLS_LOAD_SKIP_IF);
    double load = now() - start;
    if (rc) return rc;
    stata_data *data = &fd.data;
    size_t n = data->nobs, first = n ? fd.obs_map[0] : 1;
    void **original = NULL;
    if (verify) {
        size_t selected = 0;
        for (ST_int obs = SF_in1(); obs <= SF_in2(); obs++) {
            if (!check_if || SF_ifobs(obs)) {
                if (selected >= n || fd.obs_map[selected++] != (perm_idx_t)obs) {
                    rc = 9; goto done;
                }
            }
        }
        if (selected != n) { rc = 9; goto done; }
    }
    if (verify && (rc = check(data, fd.obs_map, NULL, first))) goto done;
    if (sorted || scattered) {
        uint64_t rng = 0x123456789abcdefULL;
        for (size_t i = n; i > 1; i--) {
            rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
            size_t k = rng % i;
            perm_idx_t tmp = data->sort_order[i - 1];
            data->sort_order[i - 1] = data->sort_order[k]; data->sort_order[k] = tmp;
        }
    }
    if (scattered) {
        original = calloc(data->nvars, sizeof(void *));
        if (!original) { rc = 920; goto done; }
        for (size_t j = 0; j < data->nvars; j++) {
            stata_variable *v = data->vars + j;
            void *copy = malloc(n * (v->type == STATA_TYPE_DOUBLE ? sizeof(double) : sizeof(char *)));
            if (!copy) { rc = 920; goto done; }
            if (v->type == STATA_TYPE_DOUBLE) {
                for (size_t i = 0; i < n; i++) ((double *)copy)[i] = v->data.dbl[data->sort_order[i]];
                original[j] = v->data.dbl; v->data.dbl = copy;
            } else {
                for (size_t i = 0; i < n; i++) ((char **)copy)[i] = v->data.str[data->sort_order[i]];
                original[j] = v->data.str; v->data.str = copy;
            }
        }
    }
    start = now();
    if (filtered) {
        for (size_t j = 0; j < data->nvars && !rc; j++) {
            stata_variable *v = data->vars + j;
            rc = v->type == STATA_TYPE_DOUBLE
                ? ctools_store_filtered_rowpar(v->data.dbl, n, (int)j + 1, fd.obs_map)
                : ctools_store_filtered_str(v->data.str, n, (int)j + 1, fd.obs_map);
        }
    } else if (!read_only) {
        rc = sorted ? ctools_data_store_sorted(data, first) : ctools_data_store(data, first);
    }
    double store = now() - start;
    if (!rc && verify && !read_only)
        rc = check(data, filtered ? fd.obs_map : NULL, sorted ? data->sort_order : NULL, first);
    if (original) {
        for (size_t j = 0; j < data->nvars; j++) {
            stata_variable *v = data->vars + j;
            if (original[j]) {
                if (v->type == STATA_TYPE_DOUBLE) { free(v->data.dbl); v->data.dbl = original[j]; }
                else { free(v->data.str); v->data.str = original[j]; }
            }
        }
        free(original); original = NULL;
    }
    if (!rc && (sorted || scattered)) rc = ctools_data_store(data, first);
    if (!rc) {
        char line[512];
        snprintf(line, sizeof(line), "TRANSPORT,%s,%s,%zu,%zu,%d,%.9f,%.9f,%d\n",
                 argv[0], argv[2], n, data->nvars, atoi(argv[1]), load, store, verify);
        SF_display(line);
    }
done:
    if (original) {
        for (size_t j = 0; j < data->nvars; j++) if (original[j]) {
            stata_variable *v = data->vars + j;
            if (v->type == STATA_TYPE_DOUBLE) { free(v->data.dbl); v->data.dbl = original[j]; }
            else { free(v->data.str); v->data.str = original[j]; }
        }
        free(original);
    }
    ctools_filtered_data_free(&fd);
    return rc;
}
