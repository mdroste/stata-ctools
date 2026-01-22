/*
 * cdecode_impl.c
 *
 * High-performance numeric-to-string decoding for Stata
 * Part of the ctools suite
 *
 * Algorithm:
 * 1. Parse value-label mappings passed from Stata
 * 2. Build hash table mapping integer values to string labels
 * 3. For each observation:
 *    a. Read numeric value
 *    b. Look up corresponding label in hash table
 *    c. Write label string to destination variable
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_config.h"
#include "ctools_arena.h"
#include "cdecode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CDECODE_HASH_INIT_SIZE  256
#define CDECODE_HASH_LOAD_FACTOR 0.75
#define CDECODE_MAX_LABELS      65536
#define CDECODE_MAX_LABEL_LEN   2045

/* ============================================================================
 * Hash Table for Integer -> String Mapping
 * ============================================================================ */

typedef struct {
    int key;           /* Numeric value */
    char *value;       /* Label string */
    int occupied;      /* 1 if slot is used */
} cdecode_hash_entry;

typedef struct {
    cdecode_hash_entry *entries;
    size_t capacity;
    size_t count;
    ctools_string_arena *arena;
} cdecode_hash_table;

static int cdecode_hash_init(cdecode_hash_table *ht, size_t initial_capacity)
{
    ht->capacity = initial_capacity;
    ht->count = 0;
    ht->entries = calloc(initial_capacity, sizeof(cdecode_hash_entry));
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

static void cdecode_hash_free(cdecode_hash_table *ht)
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

static int cdecode_hash_resize(cdecode_hash_table *ht, size_t new_capacity)
{
    cdecode_hash_entry *old_entries = ht->entries;
    size_t old_capacity = ht->capacity;

    ht->entries = calloc(new_capacity, sizeof(cdecode_hash_entry));
    if (!ht->entries) {
        ht->entries = old_entries;
        return -1;
    }
    ht->capacity = new_capacity;

    for (size_t i = 0; i < old_capacity; i++) {
        if (old_entries[i].occupied) {
            /* Simple hash function for integers */
            unsigned int hash = (unsigned int)old_entries[i].key;
            hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
            hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
            hash = (hash >> 16) ^ hash;

            size_t idx = hash % new_capacity;
            while (ht->entries[idx].occupied) {
                idx = (idx + 1) % new_capacity;
            }
            ht->entries[idx] = old_entries[i];
        }
    }

    free(old_entries);
    return 0;
}

static int cdecode_hash_insert(cdecode_hash_table *ht, int key, const char *value)
{
    if ((double)ht->count / ht->capacity >= CDECODE_HASH_LOAD_FACTOR) {
        if (cdecode_hash_resize(ht, ht->capacity * 2) != 0) {
            return -1;
        }
    }

    /* Simple hash function for integers */
    unsigned int hash = (unsigned int)key;
    hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
    hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
    hash = (hash >> 16) ^ hash;

    size_t idx = hash % ht->capacity;

    while (ht->entries[idx].occupied) {
        if (ht->entries[idx].key == key) {
            /* Update existing entry */
            return 0;
        }
        idx = (idx + 1) % ht->capacity;
    }

    char *value_copy = ctools_string_arena_strdup(ht->arena, value);
    if (!value_copy) return -1;

    ht->entries[idx].key = key;
    ht->entries[idx].value = value_copy;
    ht->entries[idx].occupied = 1;
    ht->count++;

    return 0;
}

static const char *cdecode_hash_lookup(cdecode_hash_table *ht, int key)
{
    unsigned int hash = (unsigned int)key;
    hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
    hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
    hash = (hash >> 16) ^ hash;

    size_t idx = hash % ht->capacity;
    size_t start = idx;

    while (ht->entries[idx].occupied) {
        if (ht->entries[idx].key == key) {
            return ht->entries[idx].value;
        }
        idx = (idx + 1) % ht->capacity;
        if (idx == start) break; /* Wrapped around */
    }

    return NULL;
}

/* ============================================================================
 * Label Parsing
 * ============================================================================ */

/*
 * Parse value-label pairs from encoded string.
 * Format: "value|label||value|label||..."
 * Returns 0 on success, -1 on error.
 */
static int parse_labels(const char *labels_str, cdecode_hash_table *ht, int *max_len)
{
    if (!labels_str || *labels_str == '\0') {
        return 0;  /* No labels is OK */
    }

    char *labels_copy = strdup(labels_str);
    if (!labels_copy) return -1;

    *max_len = 0;
    char *ptr = labels_copy;

    while (*ptr) {
        /* Find value */
        char *pipe = strchr(ptr, '|');
        if (!pipe) break;

        *pipe = '\0';
        int value = atoi(ptr);
        ptr = pipe + 1;

        /* Find end of label (either || or end of string) */
        char *end = strstr(ptr, "||");
        char *label;

        if (end) {
            *end = '\0';
            label = ptr;
            ptr = end + 2;
        } else {
            label = ptr;
            ptr = ptr + strlen(ptr);
        }

        /* Unescape special characters */
        char unescaped[CDECODE_MAX_LABEL_LEN + 1];
        size_t j = 0;
        for (size_t i = 0; label[i] && j < CDECODE_MAX_LABEL_LEN; i++) {
            if (label[i] == '\\' && label[i+1]) {
                i++;
                if (label[i] == '|') {
                    unescaped[j++] = '|';
                } else if (label[i] == '\\') {
                    unescaped[j++] = '\\';
                } else if (label[i] == '"') {
                    unescaped[j++] = '"';
                } else {
                    unescaped[j++] = label[i];
                }
            } else {
                unescaped[j++] = label[i];
            }
        }
        unescaped[j] = '\0';

        /* Track max length */
        int len = (int)strlen(unescaped);
        if (len > *max_len) {
            *max_len = len;
        }

        /* Insert into hash table */
        if (cdecode_hash_insert(ht, value, unescaped) != 0) {
            free(labels_copy);
            return -1;
        }
    }

    free(labels_copy);
    return 0;
}

/* ============================================================================
 * Main Implementation
 * ============================================================================ */

ST_retcode cdecode_main(const char *args)
{
    double t_start, t_parse, t_decode, t_total;
    int src_idx = 0, dst_idx = 0;
    int maxlen = 0;
    const char *labels_str = NULL;
    size_t nobs;
    ST_retcode rc = 0;

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cdecode: no arguments specified\n");
        return 198;
    }

    /* Make a copy for parsing */
    size_t args_len = strlen(args);
    char *args_copy = malloc(args_len + 1);
    if (!args_copy) {
        SF_error("cdecode: memory allocation failed\n");
        return 920;
    }
    strcpy(args_copy, args);

    /* Find labels= parameter first (it may contain spaces in quoted strings) */
    char *labels_ptr = strstr(args_copy, "labels=");
    if (labels_ptr) {
        labels_str = labels_ptr + 7;
        *labels_ptr = '\0'; /* Terminate args before labels= */
    }

    /* Parse other arguments */
    char *token = strtok(args_copy, " \t");
    if (!token) {
        free(args_copy);
        SF_error("cdecode: missing source variable index\n");
        return 198;
    }
    src_idx = atoi(token);

    token = strtok(NULL, " \t");
    if (!token) {
        free(args_copy);
        SF_error("cdecode: missing destination variable index\n");
        return 198;
    }
    dst_idx = atoi(token);

    /* Parse optional maxlen */
    while ((token = strtok(NULL, " \t")) != NULL) {
        if (strncmp(token, "maxlen=", 7) == 0) {
            maxlen = atoi(token + 7);
        }
    }

    if (src_idx < 1 || dst_idx < 1) {
        free(args_copy);
        SF_error("cdecode: invalid variable indices\n");
        return 198;
    }

    /* Source should be numeric */
    if (SF_var_is_string(src_idx)) {
        free(args_copy);
        SF_error("cdecode: source variable must be numeric\n");
        return 198;
    }

    /* Get observation range */
    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        free(args_copy);
        SF_scal_save("_cdecode_n_decoded", 0);
        SF_scal_save("_cdecode_n_missing", 0);
        SF_scal_save("_cdecode_n_unlabeled", 0);
        return 0;
    }

    t_parse = ctools_timer_seconds();

    /* Build hash table from labels */
    cdecode_hash_table ht;
    if (cdecode_hash_init(&ht, CDECODE_HASH_INIT_SIZE) != 0) {
        free(args_copy);
        SF_error("cdecode: memory allocation failed for hash table\n");
        return 920;
    }

    int detected_maxlen = 0;
    if (labels_str && *labels_str) {
        if (parse_labels(labels_str, &ht, &detected_maxlen) != 0) {
            cdecode_hash_free(&ht);
            free(args_copy);
            SF_error("cdecode: failed to parse labels\n");
            return 920;
        }
    }

    free(args_copy);
    args_copy = NULL;

    /* Use detected maxlen if not specified */
    if (maxlen == 0) {
        maxlen = detected_maxlen > 0 ? detected_maxlen : 1;
    }

    /* Counters */
    int n_decoded = 0;
    int n_missing = 0;
    int n_unlabeled = 0;

    /* Process observations */
    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs1 + (ST_int)i;

        /* Check if observation is in the sample (respects if/in conditions) */
        if (!SF_ifobs(obs)) {
            /* Not in sample - store empty string */
            SF_sstore(dst_idx, obs, "");
            continue;
        }

        double dval;
        SF_vdata(src_idx, obs, &dval);

        if (SF_is_missing(dval)) {
            /* Missing value -> empty string */
            SF_sstore(dst_idx, obs, "");
            n_missing++;
            continue;
        }

        /* Convert to integer for label lookup */
        int ival = (int)round(dval);

        /* Look up label */
        const char *label = cdecode_hash_lookup(&ht, ival);

        if (label) {
            SF_sstore(dst_idx, obs, (char *)label);
            n_decoded++;
        } else {
            /* No label found - store empty string */
            SF_sstore(dst_idx, obs, "");
            n_unlabeled++;
        }
    }

    t_decode = ctools_timer_seconds();
    t_total = t_decode - t_start;

    cdecode_hash_free(&ht);

    /* Store results */
    SF_scal_save("_cdecode_n_decoded", (double)n_decoded);
    SF_scal_save("_cdecode_n_missing", (double)n_missing);
    SF_scal_save("_cdecode_n_unlabeled", (double)n_unlabeled);
    SF_scal_save("_cdecode_n_labels", (double)ht.count);
    SF_scal_save("_cdecode_max_len", (double)detected_maxlen);
    SF_scal_save("_cdecode_time_parse", t_parse - t_start);
    SF_scal_save("_cdecode_time_decode", t_decode - t_parse);
    SF_scal_save("_cdecode_time_total", t_total);

#ifdef _OPENMP
    SF_scal_save("_cdecode_openmp_enabled", 1.0);
    SF_scal_save("_cdecode_threads_max", (double)omp_get_max_threads());
#else
    SF_scal_save("_cdecode_openmp_enabled", 0.0);
    SF_scal_save("_cdecode_threads_max", 1.0);
#endif

    return rc;
}
