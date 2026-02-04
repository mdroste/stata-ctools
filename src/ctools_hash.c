/*
 * ctools_hash.c
 * Unified hash table implementations for ctools
 *
 * Extracted from cencode_impl.c and cdecode_impl.c for code reuse.
 */

#include <stdlib.h>
#include <string.h>
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
            /* Key already exists - return existing value */
            return ht->entries[idx].value;
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

    return value;
}

int ctools_str_hash_lookup(ctools_str_hash_table *ht, const char *key)
{
    if (!ht || !key || ht->count == 0) return 0;

    uint32_t hash = ctools_str_hash_compute(key);
    size_t idx = hash % ht->capacity;
    size_t probes = 0;

    while (ht->entries[idx].key && probes < ht->capacity) {
        if (ht->entries[idx].hash == hash &&
            strcmp(ht->entries[idx].key, key) == 0) {
            return ht->entries[idx].value;
        }
        idx = (idx + 1) % ht->capacity;
        probes++;
    }

    return 0;  /* Not found */
}

/* ============================================================================
 * Integer -> String Hash Table Implementation
 * ============================================================================ */

/*
 * Integer hash function using bit mixing.
 * Provides good distribution for sequential and sparse integer keys.
 */
uint32_t ctools_int_hash_compute(int key)
{
    uint32_t hash = (uint32_t)key;
    hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
    hash = ((hash >> 16) ^ hash) * 0x45d9f3b;
    hash = (hash >> 16) ^ hash;
    return hash;
}

int ctools_int_hash_init(ctools_int_hash_table *ht, size_t initial_capacity)
{
    if (!ht) return -1;

    ht->capacity = initial_capacity;
    ht->count = 0;
    ht->entries = calloc(initial_capacity, sizeof(ctools_int_hash_entry));
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

void ctools_int_hash_free(ctools_int_hash_table *ht)
{
    if (!ht) return;

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

static int ctools_int_hash_resize(ctools_int_hash_table *ht, size_t new_capacity)
{
    ctools_int_hash_entry *old_entries = ht->entries;
    size_t old_capacity = ht->capacity;

    ht->entries = calloc(new_capacity, sizeof(ctools_int_hash_entry));
    if (!ht->entries) {
        ht->entries = old_entries;
        return -1;
    }
    ht->capacity = new_capacity;

    /* Rehash all existing entries */
    for (size_t i = 0; i < old_capacity; i++) {
        if (old_entries[i].occupied) {
            uint32_t hash = ctools_int_hash_compute(old_entries[i].key);
            size_t idx = hash % new_capacity;
            size_t probes = 0;

            while (ht->entries[idx].occupied && probes < new_capacity) {
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

int ctools_int_hash_insert(ctools_int_hash_table *ht, int key, const char *value)
{
    if (!ht || !value) return -1;

    /* Resize if load factor exceeded */
    if ((double)ht->count / ht->capacity >= CTOOLS_HASH_LOAD_FACTOR) {
        if (ctools_int_hash_resize(ht, ht->capacity * 2) != 0) {
            return -1;
        }
    }

    uint32_t hash = ctools_int_hash_compute(key);
    size_t idx = hash % ht->capacity;
    size_t probes = 0;

    /* Linear probing */
    while (ht->entries[idx].occupied && probes < ht->capacity) {
        if (ht->entries[idx].key == key) {
            /* Key already exists - do not update */
            return 0;
        }
        idx = (idx + 1) % ht->capacity;
        probes++;
    }

    if (probes >= ht->capacity) {
        return -1;  /* Table full */
    }

    /* Insert new entry */
    char *value_copy = ctools_string_arena_strdup(ht->arena, value);
    if (!value_copy) return -1;

    ht->entries[idx].key = key;
    ht->entries[idx].value = value_copy;
    ht->entries[idx].occupied = 1;
    ht->count++;

    return 0;
}

const char *ctools_int_hash_lookup(ctools_int_hash_table *ht, int key)
{
    if (!ht || ht->count == 0) return NULL;

    uint32_t hash = ctools_int_hash_compute(key);
    size_t idx = hash % ht->capacity;
    size_t start = idx;

    while (ht->entries[idx].occupied) {
        if (ht->entries[idx].key == key) {
            return ht->entries[idx].value;
        }
        idx = (idx + 1) % ht->capacity;
        if (idx == start) break;  /* Wrapped around */
    }

    return NULL;  /* Not found */
}
