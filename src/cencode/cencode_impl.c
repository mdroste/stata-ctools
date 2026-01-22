/*
 * cencode_impl.c
 *
 * High-performance string encoding for Stata
 * Part of the ctools suite
 *
 * Optimized streaming algorithm:
 * 1. Pass 1: Stream through strings, building hash table of ONLY unique values
 *    - No need to store all N strings, just K unique ones (K << N typically)
 * 2. Sort unique strings alphabetically, assign codes 1, 2, 3, ...
 * 3. Pass 2: Re-read strings from Stata, look up codes, write to output
 *
 * Note: Stata's Plugin Interface serializes data access internally.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_config.h"
#include "ctools_arena.h"
#include "cencode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CENCODE_HASH_INIT_SIZE  1024
#define CENCODE_HASH_LOAD_FACTOR 0.75
#define CENCODE_MAX_LABELS      65536
#define CENCODE_STR_BUF_SIZE    2048

/* ============================================================================
 * Hash Table for String -> Integer Mapping
 * ============================================================================ */

typedef struct {
    char *key;
    int value;
    uint32_t hash;
} cencode_hash_entry;

typedef struct {
    cencode_hash_entry *entries;
    size_t capacity;
    size_t count;
    ctools_string_arena *arena;
} cencode_hash_table;

static inline uint32_t cencode_hash_string(const char *s)
{
    uint32_t hash = 2166136261u;
    while (*s) {
        hash ^= (uint8_t)*s++;
        hash *= 16777619u;
    }
    return hash;
}

static int cencode_hash_init(cencode_hash_table *ht, size_t initial_capacity)
{
    ht->capacity = initial_capacity;
    ht->count = 0;
    ht->entries = calloc(initial_capacity, sizeof(cencode_hash_entry));
    if (!ht->entries) return -1;

    ht->arena = ctools_string_arena_create(initial_capacity * 64,
                                            CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
    if (!ht->arena) {
        free(ht->entries);
        ht->entries = NULL;
        return -1;
    }
    return 0;
}

static void cencode_hash_free(cencode_hash_table *ht)
{
    if (ht->arena) {
        ctools_string_arena_free(ht->arena);
        ht->arena = NULL;
    }
    if (ht->entries) {
        free(ht->entries);
        ht->entries = NULL;
    }
    ht->count = 0;
    ht->capacity = 0;
}

static int cencode_hash_resize(cencode_hash_table *ht, size_t new_capacity)
{
    cencode_hash_entry *old_entries = ht->entries;
    size_t old_capacity = ht->capacity;

    ht->entries = calloc(new_capacity, sizeof(cencode_hash_entry));
    if (!ht->entries) {
        ht->entries = old_entries;
        return -1;
    }
    ht->capacity = new_capacity;

    for (size_t i = 0; i < old_capacity; i++) {
        if (old_entries[i].key) {
            uint32_t hash = old_entries[i].hash;
            size_t idx = hash % new_capacity;
            while (ht->entries[idx].key) {
                idx = (idx + 1) % new_capacity;
            }
            ht->entries[idx] = old_entries[i];
        }
    }

    free(old_entries);
    return 0;
}

static int cencode_hash_insert(cencode_hash_table *ht, const char *key, uint32_t precomputed_hash)
{
    if ((double)ht->count / ht->capacity >= CENCODE_HASH_LOAD_FACTOR) {
        if (cencode_hash_resize(ht, ht->capacity * 2) != 0) {
            return -1;
        }
    }

    size_t idx = precomputed_hash % ht->capacity;

    while (ht->entries[idx].key) {
        if (ht->entries[idx].hash == precomputed_hash &&
            strcmp(ht->entries[idx].key, key) == 0) {
            return ht->entries[idx].value;
        }
        idx = (idx + 1) % ht->capacity;
    }

    char *key_copy = ctools_string_arena_strdup(ht->arena, key);
    if (!key_copy) return -1;

    ht->entries[idx].key = key_copy;
    ht->entries[idx].hash = precomputed_hash;
    ht->entries[idx].value = (int)(ht->count + 1);
    ht->count++;

    return ht->entries[idx].value;
}

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
 * Main Implementation
 * ============================================================================ */

ST_retcode cencode_main(const char *args)
{
    double t_start, t_parse, t_collect, t_sort, t_encode, t_labels, t_total;
    int var_idx = 0, gen_idx = 0;
    char label_name[256] = "";
    int noextend = 0;
    size_t nobs;
    ST_retcode rc = 0;

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cencode: no arguments specified\n");
        return 198;
    }

    char args_copy[1024];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *token = strtok(args_copy, " ");
    if (!token) {
        SF_error("cencode: missing source variable index\n");
        return 198;
    }
    var_idx = atoi(token);

    token = strtok(NULL, " ");
    if (!token) {
        SF_error("cencode: missing generate variable index\n");
        return 198;
    }
    gen_idx = atoi(token);

    while ((token = strtok(NULL, " ")) != NULL) {
        if (strncmp(token, "label=", 6) == 0) {
            strncpy(label_name, token + 6, sizeof(label_name) - 1);
            label_name[sizeof(label_name) - 1] = '\0';
        } else if (strcmp(token, "noextend") == 0) {
            noextend = 1;
        }
    }
    (void)noextend;

    if (var_idx < 1 || gen_idx < 1) {
        SF_error("cencode: invalid variable indices\n");
        return 198;
    }

    if (!SF_var_is_string(var_idx)) {
        SF_error("cencode: source variable must be a string variable\n");
        return 198;
    }

    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        SF_error("cencode: no observations in sample\n");
        return 2000;
    }

    t_parse = ctools_timer_seconds();
    double time_parse = t_parse - t_start;

    int is_strl = SF_var_is_strl(var_idx);

    /* ========================================================================
     * Pass 1: Stream through strings, collect unique values AND store codes
     * We store the original code for each observation to avoid re-reading
     * strings in Pass 2.
     * ======================================================================== */

    cencode_hash_table ht;
    if (cencode_hash_init(&ht, CENCODE_HASH_INIT_SIZE) != 0) {
        SF_error("cencode: memory allocation failed for hash table\n");
        return 920;
    }

    /* Allocate array to store original codes for each observation */
    int *obs_codes = malloc(nobs * sizeof(int));
    if (!obs_codes) {
        cencode_hash_free(&ht);
        SF_error("cencode: memory allocation failed for observation codes\n");
        return 920;
    }
    memset(obs_codes, 0, nobs * sizeof(int));  /* 0 = missing/empty */

    char *buf = malloc(CENCODE_STR_BUF_SIZE + 1);
    char *large_buf = NULL;
    size_t large_buf_size = 0;

    if (!buf) {
        free(obs_codes);
        cencode_hash_free(&ht);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }

    int collect_error = 0;

    for (size_t i = 0; i < nobs && !collect_error; i++) {
        ST_int obs = obs1 + (ST_int)i;
        int slen = SF_sdatalen(var_idx, obs);

        if (slen <= 0) {
            obs_codes[i] = 0;  /* missing */
            continue;
        }

        char *str_ptr;

        if ((size_t)slen <= CENCODE_STR_BUF_SIZE) {
            if (is_strl) {
                SF_strldata(var_idx, obs, buf, slen + 1);
            } else {
                SF_sdata(var_idx, obs, buf);
            }
            buf[slen] = '\0';
            str_ptr = buf;
        } else {
            if ((size_t)slen > large_buf_size) {
                free(large_buf);
                large_buf_size = (size_t)slen + 256;
                large_buf = malloc(large_buf_size + 1);
                if (!large_buf) {
                    collect_error = 1;
                    continue;
                }
            }
            if (is_strl) {
                SF_strldata(var_idx, obs, large_buf, slen + 1);
            } else {
                SF_sdata(var_idx, obs, large_buf);
            }
            large_buf[slen] = '\0';
            str_ptr = large_buf;
        }

        uint32_t hash = cencode_hash_string(str_ptr);
        int code = cencode_hash_insert(&ht, str_ptr, hash);
        if (code < 0) {
            collect_error = 1;
        } else {
            obs_codes[i] = code;  /* Store the original code */
        }
    }

    free(large_buf);
    large_buf = NULL;
    large_buf_size = 0;
    free(buf);
    buf = NULL;

    if (collect_error) {
        free(obs_codes);
        cencode_hash_free(&ht);
        SF_error("cencode: memory allocation failed during collection\n");
        return 920;
    }

    t_collect = ctools_timer_seconds();
    double time_collect = t_collect - t_parse;

    size_t n_unique = ht.count;

    if (n_unique == 0) {
        for (size_t i = 0; i < nobs; i++) {
            SF_vstore(gen_idx, obs1 + (ST_int)i, SV_missval);
        }
        SF_scal_save("_cencode_n_unique", 0);
        SF_scal_save("_cencode_n_chunks", 1);
        SF_macro_save("_cencode_labels_0", "");

        free(obs_codes);
        cencode_hash_free(&ht);

        t_total = ctools_timer_seconds();
        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", time_collect);
        SF_scal_save("_cencode_time_collect", 0);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", 0);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", t_total - t_start);
        return 0;
    }

    if (n_unique > CENCODE_MAX_LABELS) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cencode: too many unique values (%zu, max %d)\n",
                 n_unique, CENCODE_MAX_LABELS);
        SF_error(msg);
        free(obs_codes);
        cencode_hash_free(&ht);
        return 198;
    }

    /* ========================================================================
     * Sort unique strings and build code mapping
     * ======================================================================== */

    cencode_string_entry *sorted_strings = malloc(n_unique * sizeof(cencode_string_entry));
    if (!sorted_strings) {
        free(obs_codes);
        cencode_hash_free(&ht);
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

    int *code_map = malloc((n_unique + 1) * sizeof(int));
    if (!code_map) {
        free(sorted_strings);
        free(obs_codes);
        cencode_hash_free(&ht);
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
     * No need to re-read strings - we stored the original codes in Pass 1
     * ======================================================================== */

    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs1 + (ST_int)i;
        int orig_code = obs_codes[i];

        if (orig_code > 0) {
            SF_vstore(gen_idx, obs, (double)code_map[orig_code]);
        } else {
            SF_vstore(gen_idx, obs, SV_missval);
        }
    }

    free(obs_codes);

    t_encode = ctools_timer_seconds();
    double time_encode = t_encode - t_sort;

    /* ========================================================================
     * Store label mapping in Stata macros
     * ======================================================================== */

    SF_scal_save("_cencode_n_unique", (double)n_unique);

    size_t chunk = 0;
    size_t labels_in_chunk = 0;
    char macro_name[64];
    char *macro_buf = malloc(32000);
    size_t macro_pos = 0;

    if (!macro_buf) {
        free(code_map);
        free(sorted_strings);
        cencode_hash_free(&ht);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }

    macro_buf[0] = '\0';

    for (size_t i = 0; i < n_unique; i++) {
        int code = sorted_strings[i].final_code;
        const char *str = sorted_strings[i].str;

        char escaped[2048];
        size_t esc_pos = 0;
        for (const char *p = str; *p && esc_pos < sizeof(escaped) - 2; p++) {
            if (*p == '"' || *p == '\\' || *p == '|') {
                escaped[esc_pos++] = '\\';
            }
            escaped[esc_pos++] = *p;
        }
        escaped[esc_pos] = '\0';

        char entry[2200];
        int entry_len = snprintf(entry, sizeof(entry), "%s%d|%s",
                                 (labels_in_chunk > 0) ? "||" : "",
                                 code, escaped);

        if (macro_pos + entry_len >= 31000) {
            snprintf(macro_name, sizeof(macro_name), "_cencode_labels_%zu", chunk);
            SF_macro_save(macro_name, macro_buf);

            chunk++;
            labels_in_chunk = 0;
            macro_pos = 0;
            macro_buf[0] = '\0';

            entry_len = snprintf(entry, sizeof(entry), "%d|%s", code, escaped);
        }

        memcpy(macro_buf + macro_pos, entry, entry_len + 1);
        macro_pos += entry_len;
        labels_in_chunk++;
    }

    if (labels_in_chunk > 0) {
        snprintf(macro_name, sizeof(macro_name), "_cencode_labels_%zu", chunk);
        SF_macro_save(macro_name, macro_buf);
    }

    SF_scal_save("_cencode_n_chunks", (double)(chunk + 1));

    if (label_name[0] != '\0') {
        SF_macro_save("_cencode_label_name", label_name);
    }

    free(macro_buf);

    t_labels = ctools_timer_seconds();
    double time_labels = t_labels - t_encode;

    t_total = ctools_timer_seconds();
    double time_total = t_total - t_start;

    free(code_map);
    free(sorted_strings);
    cencode_hash_free(&ht);

    /* Store timing */
    SF_scal_save("_cencode_time_parse", time_parse);
    SF_scal_save("_cencode_time_load", time_collect);
    SF_scal_save("_cencode_time_collect", 0);
    SF_scal_save("_cencode_time_sort", time_sort);
    SF_scal_save("_cencode_time_encode", time_encode);
    SF_scal_save("_cencode_time_labels", time_labels);
    SF_scal_save("_cencode_time_total", time_total);

#ifdef _OPENMP
    SF_scal_save("_cencode_openmp_enabled", 1.0);
    SF_scal_save("_cencode_threads_max", (double)omp_get_max_threads());
#else
    SF_scal_save("_cencode_openmp_enabled", 0.0);
    SF_scal_save("_cencode_threads_max", 1.0);
#endif

    return rc;
}
