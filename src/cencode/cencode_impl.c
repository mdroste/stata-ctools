/*
 * cencode_impl.c
 *
 * High-performance string encoding for Stata
 * Part of the ctools suite
 *
 * Algorithm (fixed-width string sources):
 * 1. Each worker reads its own run of observations straight from Stata and
 *    records, per observation, its string's code in a worker-local hash table
 * 2. Merge the workers' unique strings, sort them bytewise, and number them
 *    1, 2, 3, ...
 * 3. Each worker writes the final codes for its observations
 * With an existing value label, strings it already labels keep their (smallest)
 * value and the others are numbered after its largest value, as encode does.
 * strL sources load the column first (ctools_data_load) and encode it serially.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"
#include "ctools_runtime.h"
#include "ctools_config.h"
#include "ctools_hash.h"
#include "ctools_parse.h"
#include "cencode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CENCODE_HASH_INIT_SIZE  1024

/* ============================================================================
 * Sorted String List
 * ============================================================================ */

typedef struct {
    char *str;
    int original_code;
    int final_code;
} cencode_string_entry;

static int cencode_string_compare(const void *a, const void *b)
{
    const cencode_string_entry *ea = (const cencode_string_entry *)a;
    const cencode_string_entry *eb = (const cencode_string_entry *)b;
    return strcmp(ea->str, eb->str);
}

/* ============================================================================
 * Streaming encoder for fixed-width string sources
 * ============================================================================ */

/* Longest fixed-width string value, plus its terminator. */
#define CENCODE_STR_BYTES 2046

typedef struct {
    ctools_str_hash_table table; /* this chunk's unique strings -> local codes 1.. */
    int *final;                  /* owned [table.count + 1]: local code -> final code */
    size_t start, end;           /* zero-based rows of the range */
    int status;                  /* 0, or a Stata return code */
} cencode_chunk;

/* Largest code a long variable can hold. */
#define CENCODE_MAX_CODE 2147483620.0

typedef struct {
    int64_t value;
    const char *text;
} cencode_label_entry;

static int cencode_label_entry_compare(const void *a, const void *b)
{
    int64_t x = ((const cencode_label_entry *)a)->value;
    int64_t y = ((const cencode_label_entry *)b)->value;
    return (x > y) - (x < y);
}

/*
    Existing label handed over by cencode.ado in local __ce_map (maplen bytes):
    "<count> " then "<value> <bytes> <text> " per mapping, then "|". Values are
    integers. Fills table (text -> smallest value with that text), *count and
    *max_value. *map receives the owned macro copy that the keys point into.
    Returns 0, or a Stata return code.
*/
static int cencode_label_read(size_t maplen, char **map, ctools_str_hash_table *table,
                              size_t *count, int64_t *max_value)
{
    char *buf = (char *)calloc(maplen + 2, 1), *p, *end, *stop;
    cencode_label_entry *entries = NULL;
    int rc = 0;
    *map = buf;
    *count = 0;
    *max_value = 0;
    if (!buf) return 920;
    if (SF_macro_use("___ce_map", buf, (ST_int)maplen + 1) != 0 || strlen(buf) != maplen)
        return 198;
    p = buf;
    end = buf + maplen;
    long declared = strtol(p, &stop, 10);
    if (stop == p || declared < 0 || declared > CTOOLS_MAX_LABELS || stop >= end || *stop != ' ')
        return 198;
    p = stop + 1;
    entries = (cencode_label_entry *)calloc((size_t)declared + 1, sizeof(*entries));
    if (!entries) return 920;
    for (long e = 0; e < declared && !rc; e++) {
        long long value = strtoll(p, &stop, 10);
        if (stop == p || stop >= end || *stop != ' ' || value < -2147483647LL ||
            value > 2147483620LL) {
            rc = 198;
            break;
        }
        p = stop + 1;
        long bytes = strtol(p, &stop, 10);
        if (stop == p || bytes < 0 || bytes > CENCODE_STR_BYTES || stop >= end || *stop != ' ') {
            rc = 198;
            break;
        }
        p = stop + 1;
        if ((size_t)(end - p) < (size_t)bytes + 1 || p[bytes] != ' ') {
            rc = 198;
            break;
        }
        p[bytes] = '\0';
        entries[e].value = (int64_t)value;
        entries[e].text = p;
        p += bytes + 1;
    }
    if (!rc) {
        /* Ascending values: the first insert of a text keeps the smallest value. */
        qsort(entries, (size_t)declared, sizeof(*entries), cencode_label_entry_compare);
        for (long e = 0; e < declared && !rc; e++) {
            if (entries[e].text[0] != '\0' &&
                ctools_str_hash_insert_value(table, entries[e].text, (int)entries[e].value) != 0)
                rc = 920;
        }
        if (declared > 0) *max_value = entries[declared - 1].value;
        *count = (size_t)declared;
    }
    free(entries);
    return rc;
}

static void cencode_save_times(double parse, double load, double collect, double sort,
                               double encode, double labels, double total)
{
    SF_scal_save("_cencode_time_parse", parse);
    SF_scal_save("_cencode_time_load", load);
    SF_scal_save("_cencode_time_collect", collect);
    SF_scal_save("_cencode_time_sort", sort);
    SF_scal_save("_cencode_time_encode", encode);
    SF_scal_save("_cencode_time_labels", labels);
    SF_scal_save("_cencode_time_total", total);
}

/*
    Encode the fixed-width string variable var_idx into numeric gen_idx for
    observations SF_in1()..SF_in2(). has_if evaluates the if condition once per
    observation on the calling thread. fill_all writes every observation of
    the range, missing where there is no code, for a destination created
    without initial values; otherwise only coded observations are written.

    label_len > 0 extends the existing label in local __ce_map (see
    cencode_label_read); the label file then lists only the new mappings. When
    encode itself must decide (noextend with a string the label lacks, codes
    past the long range, or too many mappings), nothing is written and
    scalar _cencode_fallback is 1.
*/
static ST_retcode cencode_stream(int var_idx, int gen_idx, const char *labelfile, int fill_all,
                                 int has_if, int touse, size_t label_len, int noextend,
                                 double t_start)
{
    double t_parse = ctools_timer_seconds();
    ST_int obs1 = SF_in1(), obs2 = SF_in2();
    size_t nobs = obs2 >= obs1 ? (size_t)(obs2 - obs1) + 1 : 0;
    int rc = 0;

    size_t nchunks = 1;
    int workers = ctools_get_max_threads();
    if (workers > 1 && nobs >= (size_t)MIN_OBS_PER_THREAD * 2) {
        nchunks = nobs / MIN_OBS_PER_THREAD;
        if (nchunks > (size_t)workers) nchunks = (size_t)workers;
    }

    uint8_t *selected = NULL;
    uint32_t *codes = nobs ? (uint32_t *)ctools_safe_malloc2(nobs, sizeof(uint32_t)) : NULL;
    cencode_chunk *chunks = (cencode_chunk *)calloc(nchunks, sizeof(cencode_chunk));
    if ((nobs && !codes) || !chunks) rc = 920;
    if (!rc && has_if && nobs) {
        selected = (uint8_t *)calloc(nobs, 1);
        if (!selected) rc = 920;
        if (!rc && touse) {
            /* A sample marker follows the other variables in place of the if
             * condition; the shared loader reads it in parallel. */
            ctools_filtered_data marked;
            stata_retcode select_rc = ctools_data_select(&marked, 0, 0, CTOOLS_LOAD_MARKER);
            if (select_rc != STATA_OK) rc = ctools_stata_rc(select_rc);
            for (size_t k = 0; !rc && k < marked.data.nobs; k++)
                selected[marked.obs_map[k] - (perm_idx_t)obs1] = 1;
            ctools_filtered_data_free(&marked);
        }
        for (size_t i = 0; !rc && !touse && i < nobs; i++) selected[i] = SF_ifobs(obs1 + (ST_int)i) != 0;
    }
    for (size_t c = 0; !rc && c < nchunks; c++) {
        size_t base = nobs / nchunks, extra = nobs % nchunks;
        chunks[c].start = c * base + (c < extra ? c : extra);
        chunks[c].end = chunks[c].start + base + (c < extra);
        if (ctools_str_hash_init(&chunks[c].table, CENCODE_HASH_INIT_SIZE) != 0) rc = 920;
    }

    /* Pass 1: every worker reads its observations and codes them locally. */
    if (!rc) {
        const ST_IIIS read_fn = (_stata_)->sdata;
        #pragma omp parallel for schedule(dynamic, 1) if (nchunks > 1)
        for (size_t c = 0; c < nchunks; c++) {
            cencode_chunk *chunk = &chunks[c];
            char value[CENCODE_STR_BYTES + 2];
            for (size_t i = chunk->start; i < chunk->end; i++) {
                codes[i] = 0;
                if (selected && !selected[i]) continue;
                int io_rc = read_fn(var_idx, obs1 + (ST_int)i, value);
                if (io_rc) {
                    chunk->status = io_rc == 1 ? 1 : 459;
                    break;
                }
                if (value[0] == '\0') continue;
                int code = ctools_str_hash_insert(&chunk->table, value,
                                                  ctools_str_hash_compute(value));
                if (code < 0) {
                    chunk->status = 920;
                    break;
                }
                if (chunk->table.count > CTOOLS_MAX_LABELS) {
                    chunk->status = 134;
                    break;
                }
                codes[i] = (uint32_t)code;
            }
        }
        for (size_t c = 0; c < nchunks && !rc; c++) rc = chunks[c].status;
    }
    double t_load = ctools_timer_seconds();

    /* Merge the chunks' strings: chunk local code -> merged code. */
    ctools_str_hash_table all, existing;
    int all_ready = 0, existing_ready = 0, fallback = 0;
    char *label_map = NULL;
    size_t n_unique = 0, max_len = 0, n_existing = 0, n_new = 0;
    int64_t max_existing = 0;
    if (!rc) {
        if (ctools_str_hash_init(&all, CENCODE_HASH_INIT_SIZE) != 0) rc = 920;
        else all_ready = 1;
    }
    if (!rc && label_len > 0) {
        if (ctools_str_hash_init(&existing, CENCODE_HASH_INIT_SIZE) != 0) rc = 920;
        else {
            existing_ready = 1;
            rc = cencode_label_read(label_len, &label_map, &existing, &n_existing,
                                    &max_existing);
            if (rc == 198) SF_error("cencode: malformed value-label mapping\n");
        }
    }
    for (size_t c = 0; !rc && c < nchunks; c++) {
        ctools_str_hash_table *t = &chunks[c].table;
        chunks[c].final = (int *)calloc(t->count + 1, sizeof(int));
        if (!chunks[c].final) {
            rc = 920;
            break;
        }
        for (size_t e = 0; e < t->capacity; e++) {
            if (!t->entries[e].key) continue;
            int merged = ctools_str_hash_insert(&all, t->entries[e].key, t->entries[e].hash);
            if (merged < 0) {
                rc = 920;
                break;
            }
            chunks[c].final[t->entries[e].value] = merged;
        }
        if (!rc && all.count > CTOOLS_MAX_LABELS) rc = 134;
    }
    double t_collect = ctools_timer_seconds();

    /* Number the merged strings in bytewise order: strings the existing label
     * already holds keep its value; the others follow its largest value. The
     * first n_new entries of sorted are the new label mappings. */
    cencode_string_entry *sorted = NULL;
    int *merged_final = NULL;
    if (!rc) {
        n_unique = all.count;
        sorted = (cencode_string_entry *)malloc((n_unique + 1) * sizeof(*sorted));
        merged_final = (int *)calloc(n_unique + 1, sizeof(int));
        if (!sorted || !merged_final) rc = 920;
    }
    if (!rc) {
        size_t k = 0, kept = n_unique;
        for (size_t e = 0; e < all.capacity && k + (n_unique - kept) < n_unique; e++) {
            if (!all.entries[e].key) continue;
            const char *key = all.entries[e].key;
            size_t len = strlen(key);
            if (len > max_len) max_len = len;
            int value;
            if (existing_ready && ctools_str_hash_lookup(&existing, key, &value)) {
                merged_final[all.entries[e].value] = value;
                kept--;
                continue;
            }
            sorted[k].str = all.entries[e].key;
            sorted[k].original_code = all.entries[e].value;
            k++;
        }
        n_new = k;
        qsort(sorted, n_new, sizeof(*sorted), cencode_string_compare);
        if (existing_ready && n_new > 0 &&
            (noextend || (double)max_existing + (double)n_new > CENCODE_MAX_CODE ||
             n_existing + n_new > CTOOLS_MAX_LABELS))
            fallback = 1;
        for (size_t k2 = 0; !fallback && k2 < n_new; k2++) {
            sorted[k2].final_code = (int)(max_existing + (int64_t)k2 + 1);
            merged_final[sorted[k2].original_code] = sorted[k2].final_code;
        }
        for (size_t c = 0; !fallback && c < nchunks; c++)
            for (size_t l = 1; l <= chunks[c].table.count; l++)
                chunks[c].final[l] = merged_final[chunks[c].final[l]];
    }
    double t_sort = ctools_timer_seconds();

    /* Pass 2: each chunk writes its observations' final codes. */
    if (!rc && !fallback) {
        const ST_IIID store_fn = (_stata_)->safestore;
        #pragma omp parallel for schedule(dynamic, 1) if (nchunks > 1)
        for (size_t c = 0; c < nchunks; c++) {
            cencode_chunk *chunk = &chunks[c];
            for (size_t i = chunk->start; i < chunk->end; i++) {
                uint32_t local = codes[i];
                if (!local && !fill_all) continue;
                double value = local ? (double)chunk->final[local] : SV_missval;
                int io_rc = store_fn(gen_idx, obs1 + (ST_int)i, value);
                if (io_rc) {
                    chunk->status = io_rc == 1 ? 1 : 459;
                    break;
                }
            }
        }
        for (size_t c = 0; c < nchunks && !rc; c++) rc = chunks[c].status;
    }
    double t_encode = ctools_timer_seconds();

    /* New label mappings for the ado to install in one Mata call. */
    if (!rc && !fallback && n_new > 0 && labelfile[0] != '\0') {
        const char **label_strings = (const char **)malloc(n_new * sizeof(char *));
        int *label_codes = (int *)malloc(n_new * sizeof(int));
        if (!label_strings || !label_codes) rc = 920;
        for (size_t k = 0; !rc && k < n_new; k++) {
            label_strings[k] = sorted[k].str;
            label_codes[k] = sorted[k].final_code;
        }
        if (!rc && ctools_label_write_file(label_strings, label_codes, n_new, labelfile) != 0) {
            SF_error("cencode: failed to write label file\n");
            rc = 920;
        }
        free(label_strings);
        free(label_codes);
    }
    double t_labels = ctools_timer_seconds();

    free(sorted);
    free(merged_final);
    if (all_ready) ctools_str_hash_free(&all);
    if (existing_ready) ctools_str_hash_free(&existing);
    free(label_map);
    for (size_t c = 0; chunks && c < nchunks; c++) {
        free(chunks[c].final);
        ctools_str_hash_free(&chunks[c].table);
    }
    free(chunks);
    free(codes);
    free(selected);

    if (rc) {
        if (rc == 134)
            SF_error("cencode: too many unique values (max 65536)\n");
        else if (rc == 920)
            SF_error("cencode: memory allocation failed\n");
        else if (rc == 459)
            SF_error("cencode: failed to read or write Stata data\n");
        return rc == 134 ? 198 : rc;
    }
    SF_scal_save("_cencode_fallback", fallback ? 1.0 : 0.0);
    SF_scal_save("_cencode_n_unique", (double)n_unique);
    SF_scal_save("_cencode_n_new", (double)n_new);
    SF_scal_save("_cencode_max_len", (double)max_len);
    SF_scal_save("_cencode_protocol", 2.0);
    cencode_save_times(t_parse - t_start, t_load - t_parse, t_collect - t_load, t_sort - t_collect,
                       t_encode - t_sort, t_labels - t_encode, ctools_timer_seconds() - t_start);
    CTOOLS_SAVE_THREAD_INFO("_cencode");
    return 0;
}

/* ============================================================================
 * Main Implementation
 * ============================================================================ */

ST_retcode cencode_main(const char *args)
{
    double t_start, t_parse, t_collect, t_sort, t_encode, t_labels, t_total;
    int var_idx = 0, gen_idx = 0;
    size_t nobs;

    /*
     * cencode.ado passes the source and the destination:
     *   plugin call ctools_plugin srcvar destvar, "cencode 1 2 ..."
     *
     * var_idx = source string variable's position in the plugin varlist
     * gen_idx = destination numeric variable's position in the plugin varlist
     */

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cencode: no arguments specified\n");
        return 198;
    }

    /* First two tokens are positional: var_idx gen_idx */
    const char *cursor = args;
    if (ctools_parse_next_int(&cursor, &var_idx) != 0) {
        SF_error("cencode: invalid var_idx value\n");
        return 198;
    }
    if (ctools_parse_next_int(&cursor, &gen_idx) != 0) {
        SF_error("cencode: invalid gen_idx value\n");
        return 198;
    }

    /* Named options. labelfile=<bytes> announces the label file path in local
     * __ce_labelfile: a temporary path can contain spaces, which the
     * argument tokenizer would split. A missing or cut path is an error. */
    char labelfile[4096] = "";
    int labelfile_len = 0;
    ctools_parse_status lf_status = ctools_parse_int_option(args, "labelfile", &labelfile_len);
    if (lf_status == CTOOLS_PARSE_INVALID ||
        (lf_status == CTOOLS_PARSE_OK &&
         (labelfile_len <= 0 || (size_t)labelfile_len >= sizeof(labelfile) - 1 ||
          SF_macro_use("___ce_labelfile", labelfile, (ST_int)sizeof(labelfile) - 1) != 0 ||
          strlen(labelfile) != (size_t)labelfile_len))) {
        SF_error("cencode: invalid label file path\n");
        return 198;
    }

    if (var_idx < 1 || gen_idx < 1) {
        SF_error("cencode: invalid variable indices\n");
        return 198;
    }

    /* Check source is a string variable */
    if (!SF_var_is_string(var_idx)) {
        SF_error("cencode: source variable must be a string variable\n");
        return 198;
    }

    /* Fixed-width sources stream from Stata; strL sources load the column. */
    if (ctools_parse_bool_option(args, "stream") && !SF_var_is_strl(var_idx)) {
        if (var_idx > SF_nvars() || gen_idx > SF_nvars() || SF_var_is_string(gen_idx)) {
            SF_error("cencode: invalid variable indices\n");
            return 198;
        }
        int label_len = 0;
        if (ctools_parse_int_option(args, "labelmap", &label_len) < 0 || label_len < 0) {
            SF_error("cencode: invalid labelmap length\n");
            return 198;
        }
        return cencode_stream(var_idx, gen_idx, labelfile,
                              ctools_parse_bool_option(args, "fillall"),
                              ctools_parse_bool_option(args, "ifcond"),
                              ctools_parse_bool_option(args, "touse"), (size_t)label_len,
                              ctools_parse_bool_option(args, "noextend"), t_start);
    }

    t_parse = ctools_timer_seconds();
    double time_parse = t_parse - t_start;

    /* ========================================================================
     * Bulk load source string variable with if/in filtering
     * Uses ctools_data_load() to load only filtered observations,
     * eliminating the need for a separate if_mask array.
     * ======================================================================== */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    int load_var_idx = var_idx;
    /* Without an if condition (no ifcond) every row of the in-range is
     * selected; touse passes a sample marker last in place of the if */
    int load_flags = CTOOLS_LOAD_NO_SORT_ORDER |
                     (!ctools_parse_bool_option(args, "ifcond") ? CTOOLS_LOAD_SKIP_IF
                      : ctools_parse_bool_option(args, "touse") ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc = ctools_data_load(&filtered, &load_var_idx, 1, 0, 0, load_flags);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("cencode: failed to load source string data\n");
        return ctools_stata_rc(load_rc);
    }

    nobs = filtered.data.nobs;  /* Number of filtered observations */
    perm_idx_t *obs_map = filtered.obs_map;  /* Maps filtered index -> Stata obs */

    if (nobs == 0) {
        /* No observations match if/in - match Stata's encode behavior (return success) */
        ctools_filtered_data_free(&filtered);
        SF_scal_save("_cencode_n_unique", 0);
        SF_scal_save("_cencode_max_len", 0);
        SF_scal_save("_cencode_protocol", 2.0);
        /* (no label file or macros needed for this early-return path) */
        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", ctools_timer_seconds() - t_parse);
        SF_scal_save("_cencode_time_collect", 0);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", 0);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", ctools_timer_seconds() - t_start);
        return 0;
    }

    char **src_strings = filtered.data.vars[0].data.str;

    double time_load = ctools_timer_seconds() - t_parse;

    /* ========================================================================
     * Process loaded strings: collect unique values AND store codes
     * ======================================================================== */

    ctools_str_hash_table ht;
    if (ctools_str_hash_init(&ht, CENCODE_HASH_INIT_SIZE) != 0) {
        ctools_filtered_data_free(&filtered);
        SF_error("cencode: memory allocation failed for hash table\n");
        return 920;
    }

    /* Allocate array to store codes for each filtered observation */
    int *obs_codes = malloc(nobs * sizeof(int));
    if (!obs_codes) {
        ctools_str_hash_free(&ht);
        ctools_filtered_data_free(&filtered);
        SF_error("cencode: memory allocation failed for observation codes\n");
        return 920;
    }
    int collect_error = 0;

    /* Process filtered observations - no if_mask check needed since data is pre-filtered */
    for (size_t i = 0; i < nobs && !collect_error; i++) {
        char *str_ptr = src_strings[i];

        /* Empty or NULL string -> missing */
        if (!str_ptr || str_ptr[0] == '\0') {
            obs_codes[i] = 0;  /* missing */
            continue;
        }

        uint32_t hash = ctools_str_hash_compute(str_ptr);
        int code = ctools_str_hash_insert(&ht, str_ptr, hash);
        if (code < 0) {
            collect_error = 1;
        } else {
            obs_codes[i] = code;  /* Store the original code */
        }
    }

    if (collect_error) {
        free(obs_codes);
        ctools_str_hash_free(&ht);
        ctools_filtered_data_free(&filtered);
        SF_error("cencode: memory allocation failed during collection\n");
        return 920;
    }

    t_collect = ctools_timer_seconds();
    double time_collect = t_collect - t_parse - time_load;

    /* ========================================================================
     * Normal path - sort and create new labels
     * ======================================================================== */

    size_t n_unique = ht.count;

    if (n_unique == 0) {
        /* All values were empty/missing - write all missing */
        double *miss_values = malloc(nobs * sizeof(double));
        stata_retcode store_rc = STATA_OK;
        if (miss_values) {
            for (size_t i = 0; i < nobs; i++) miss_values[i] = SV_missval;
            store_rc = ctools_store_filtered_rowpar(miss_values, nobs, gen_idx, obs_map);
            free(miss_values);
        } else {
            store_rc = STATA_ERR_MEMORY;
        }
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        if (store_rc) return ctools_stata_rc(store_rc);
        SF_scal_save("_cencode_n_unique", 0);
        SF_scal_save("_cencode_max_len", 0);
        SF_scal_save("_cencode_protocol", 2.0);

        t_total = ctools_timer_seconds();
        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", time_load);
        SF_scal_save("_cencode_time_collect", time_collect);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", 0);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", t_total - t_start);
        return 0;
    }

    if (n_unique > CTOOLS_MAX_LABELS) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cencode: too many unique values (%zu, max %d)\n",
                 n_unique, CTOOLS_MAX_LABELS);
        SF_error(msg);
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        return 198;
    }

    /* ========================================================================
     * Sort unique strings and build code mapping
     * ======================================================================== */

    cencode_string_entry *sorted_strings = malloc(n_unique * sizeof(cencode_string_entry));
    if (!sorted_strings) {
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }

    size_t idx = 0;
    for (size_t i = 0; i < ht.capacity && idx < n_unique; i++) {
        if (ht.entries[i].key) {
            sorted_strings[idx].str = ht.entries[i].key;
            sorted_strings[idx].original_code = ht.entries[i].value;
            idx++;
        }
    }

    qsort(sorted_strings, n_unique, sizeof(cencode_string_entry), cencode_string_compare);

    /* Check for overflow in (n_unique + 1) */
    if (n_unique >= SIZE_MAX) {
        free(sorted_strings);
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        SF_error("cencode: too many unique values\n");
        return 920;
    }

    /* Use overflow-safe allocation */
    int *code_map = ctools_safe_malloc2(n_unique + 1, sizeof(int));
    if (!code_map) {
        free(sorted_strings);
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }
    memset(code_map, 0, (n_unique + 1) * sizeof(int));

    for (size_t i = 0; i < n_unique; i++) {
        sorted_strings[i].final_code = (int)(i + 1);
        code_map[sorted_strings[i].original_code] = sorted_strings[i].final_code;
    }

    t_sort = ctools_timer_seconds();
    double time_sort = t_sort - t_collect;

    /* ========================================================================
     * Pass 2: Apply code_map to stored codes, write encoded values
     * via ctools_store_filtered
     * ======================================================================== */

    /* Build double array with mapped codes and use row-parallel store */
    {
        double *enc_values = malloc(nobs * sizeof(double));
        if (!enc_values) {
            free(obs_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            SF_error("cencode: memory allocation failed\n");
            return 920;
        }
        for (size_t i = 0; i < nobs; i++) {
            int orig_code = obs_codes[i];
            enc_values[i] = (orig_code > 0) ? (double)code_map[orig_code] : SV_missval;
        }
        stata_retcode store_rc = ctools_store_filtered_rowpar(enc_values, nobs, gen_idx, obs_map);
        free(enc_values);
        if (store_rc) {
            free(obs_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            return ctools_stata_rc(store_rc);
        }
    }

    free(obs_codes);

    t_encode = ctools_timer_seconds();
    double time_encode = t_encode - t_sort;

    /* ========================================================================
     * Write label data for the .ado to install in one Mata call
     * ======================================================================== */

    size_t max_len = 0;
    for (size_t i = 0; i < n_unique; i++) {
        size_t len = strlen(sorted_strings[i].str);
        if (len > max_len) max_len = len;
    }
    SF_scal_save("_cencode_n_unique", (double)n_unique);
    SF_scal_save("_cencode_max_len", (double)max_len);
    SF_scal_save("_cencode_protocol", 2.0);

    if (labelfile[0] != '\0') {
        /* Build arrays for file writing */
        const char **label_strings = malloc(n_unique * sizeof(const char *));
        int *label_codes = malloc(n_unique * sizeof(int));
        if (!label_strings || !label_codes) {
            free(label_strings);
            free(label_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            SF_error("cencode: memory allocation failed\n");
            return 920;
        }

        for (size_t i = 0; i < n_unique; i++) {
            label_strings[i] = sorted_strings[i].str;
            label_codes[i] = sorted_strings[i].final_code;
        }

        if (ctools_label_write_file(label_strings, label_codes, n_unique, labelfile) != 0) {
            free(label_strings);
            free(label_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            SF_error("cencode: failed to write label file\n");
            return 920;
        }

        free(label_strings);
        free(label_codes);
    }

    t_labels = ctools_timer_seconds();
    double time_labels = t_labels - t_encode;

    t_total = ctools_timer_seconds();
    double time_total = t_total - t_start;

    free(code_map);
    free(sorted_strings);
    ctools_filtered_data_free(&filtered);
    ctools_str_hash_free(&ht);

    /* Store timing */
    SF_scal_save("_cencode_time_parse", time_parse);
    SF_scal_save("_cencode_time_load", time_load);
    SF_scal_save("_cencode_time_collect", time_collect);
    SF_scal_save("_cencode_time_sort", time_sort);
    SF_scal_save("_cencode_time_encode", time_encode);
    SF_scal_save("_cencode_time_labels", time_labels);
    SF_scal_save("_cencode_time_total", time_total);

    CTOOLS_SAVE_THREAD_INFO("_cencode");

    return 0;
}
