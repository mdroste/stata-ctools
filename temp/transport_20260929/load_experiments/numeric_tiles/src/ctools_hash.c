/*
 * ctools_hash.c
 * String -> integer hash table and value-label output for ctools
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_hash.h"

/* ============================================================================
 * String -> Integer Hash Table Implementation
 * ============================================================================ */

/*
 * FNV-1a hash function.
 * Fast and provides good distribution for typical string data.
 */
uint32_t ctools_str_hash_compute(const char *s)
{
    uint32_t hash = 2166136261u;  /* FNV offset basis */
    while (*s) {
        hash ^= (uint8_t)*s++;
        hash *= 16777619u;        /* FNV prime */
    }
    return hash;
}

int ctools_str_hash_init(ctools_str_hash_table *ht, size_t initial_capacity)
{
    if (!ht) return -1;

    ht->capacity = initial_capacity;
    ht->count = 0;
    ht->entries = calloc(initial_capacity, sizeof(ctools_str_hash_entry));
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

void ctools_str_hash_free(ctools_str_hash_table *ht)
{
    if (!ht) return;

    /* Free strdup'd fallback keys that aren't owned by the arena */
    if (ht->arena && ht->arena->has_fallback && ht->entries) {
        for (size_t i = 0; i < ht->capacity; i++) {
            if (ht->entries[i].key &&
                !ctools_string_arena_owns(ht->arena, ht->entries[i].key)) {
                free(ht->entries[i].key);
            }
        }
    }

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

static int ctools_str_hash_resize(ctools_str_hash_table *ht, size_t new_capacity)
{
    ctools_str_hash_entry *old_entries = ht->entries;
    size_t old_capacity = ht->capacity;

    ht->entries = calloc(new_capacity, sizeof(ctools_str_hash_entry));
    if (!ht->entries) {
        ht->entries = old_entries;
        return -1;
    }
    ht->capacity = new_capacity;

    /* Rehash all existing entries */
    for (size_t i = 0; i < old_capacity; i++) {
        if (old_entries[i].key) {
            uint32_t hash = old_entries[i].hash;
            size_t idx = hash % new_capacity;
            size_t probes = 0;

            while (ht->entries[idx].key && probes < new_capacity) {
                idx = (idx + 1) % new_capacity;
                probes++;
            }

            if (probes >= new_capacity) {
                /* Table full - should never happen with proper load factor */
                free(ht->entries);
                ht->entries = old_entries;
                ht->capacity = old_capacity;
                return -1;
            }

            ht->entries[idx] = old_entries[i];
        }
    }

    free(old_entries);
    return 0;
}

int ctools_str_hash_insert(ctools_str_hash_table *ht, const char *key,
                            uint32_t precomputed_hash)
{
    if (!ht || !key) return -1;

    /* Resize if load factor exceeded */
    if ((double)ht->count / ht->capacity >= CTOOLS_HASH_LOAD_FACTOR) {
        if (ctools_str_hash_resize(ht, ht->capacity * 2) != 0) {
            return -1;
        }
    }

    size_t idx = precomputed_hash % ht->capacity;
    size_t probes = 0;

    /* Linear probing */
    while (ht->entries[idx].key && probes < ht->capacity) {
        if (ht->entries[idx].hash == precomputed_hash &&
            strcmp(ht->entries[idx].key, key) == 0) {
            /* Key already exists - return existing value */
            return ht->entries[idx].value;
        }
        idx = (idx + 1) % ht->capacity;
        probes++;
    }

    if (probes >= ht->capacity) {
        return -1;  /* Table full */
    }

    /* Insert new entry */
    char *key_copy = ctools_string_arena_strdup(ht->arena, key);
    if (!key_copy) return -1;

    ht->entries[idx].key = key_copy;
    ht->entries[idx].hash = precomputed_hash;
    ht->entries[idx].value = (int)(ht->count + 1);  /* 1-based indexing */
    ht->count++;

    return ht->entries[idx].value;
}

int ctools_str_hash_insert_value(ctools_str_hash_table *ht, const char *key,
                                  int value)
{
    if (!ht || !key) return -1;

    /* Resize if load factor exceeded */
    if ((double)ht->count / ht->capacity >= CTOOLS_HASH_LOAD_FACTOR) {
        if (ctools_str_hash_resize(ht, ht->capacity * 2) != 0) {
            return -1;
        }
    }

    uint32_t hash = ctools_str_hash_compute(key);
    size_t idx = hash % ht->capacity;
    size_t probes = 0;

    /* Linear probing */
    while (ht->entries[idx].key && probes < ht->capacity) {
        if (ht->entries[idx].hash == hash &&
            strcmp(ht->entries[idx].key, key) == 0) {
            /* Existing mapping is retained; status is independent of its code. */
            return 0;
        }
        idx = (idx + 1) % ht->capacity;
        probes++;
    }

    if (probes >= ht->capacity) {
        return -1;  /* Table full */
    }

    /* Insert new entry with specified value */
    char *key_copy = ctools_string_arena_strdup(ht->arena, key);
    if (!key_copy) return -1;

    ht->entries[idx].key = key_copy;
    ht->entries[idx].hash = hash;
    ht->entries[idx].value = value;
    ht->count++;

    return 0;
}

int ctools_str_hash_lookup(ctools_str_hash_table *ht, const char *key, int *value)
{
    if (!ht || !key || !value || ht->count == 0) return 0;

    uint32_t hash = ctools_str_hash_compute(key);
    size_t idx = hash % ht->capacity;
    size_t probes = 0;

    while (ht->entries[idx].key && probes < ht->capacity) {
        if (ht->entries[idx].hash == hash &&
            strcmp(ht->entries[idx].key, key) == 0) {
            *value = ht->entries[idx].value;
            return 1;
        }
        idx = (idx + 1) % ht->capacity;
        probes++;
    }

    return 0;  /* Not found */
}

/* ============================================================================
 * Value Label Output (cencode)
 * ============================================================================ */

/* Emit only ASCII numeric byte literals as Mata source. User text never becomes
 * Stata macro syntax, quoting delimiters, or executable Mata code. */
int ctools_label_write_stata_file(const char **strings, const int *codes,
                                   size_t n_labels, const char *label_name,
                                   const char *filepath)
{
    FILE *fp = fopen(filepath, "w");
    if (!fp) return -1;
    fprintf(fp, "tempname __ctools_label_text\n");
    for (size_t i = 0; i < n_labels; i++) {
        const unsigned char *p = (const unsigned char *)(strings[i] ? strings[i] : "");
        size_t length = strlen((const char *)p);
        /* Bound each expression below Mata's token limit, including str2045. */
        fprintf(fp, "mata: `__ctools_label_text' = \"\"\n");
        for (size_t offset = 0; offset < length; offset += 128) {
            size_t end = offset + 128 < length ? offset + 128 : length;
            fprintf(fp, "mata: `__ctools_label_text' = `__ctools_label_text' + invtokens(char((");
            for (size_t j = offset; j < end; j++)
                fprintf(fp, "%s%u", j == offset ? "" : ",", (unsigned)p[j]);
            fprintf(fp, ")), \"\")\n");
        }
        fprintf(fp, "mata: st_vlmodify(\"%s\", %d, `__ctools_label_text')\n", label_name, codes[i]);
    }
    if (n_labels) fprintf(fp, "mata: mata drop `__ctools_label_text'\n");
    int failed = ferror(fp);
    if (fclose(fp) != 0) failed = 1;
    return failed ? -1 : 0;
}
