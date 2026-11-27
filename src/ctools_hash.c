/*
 * ctools_hash.c
 * String -> integer hash table and value-label output for ctools
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ctools_file.h"

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_hash.h"

/* ============================================================================
 * String -> Integer Hash Table Implementation
 * ============================================================================ */

/*
 * Multiply-xorshift over 8-byte words (the tail zero-padded), seeded with
 * the length. Several times faster than a byte-at-a-time hash on the
 * strings cencode sees, and well mixed in the low bits the table uses.
 */
uint32_t ctools_str_hash_compute(const char *s)
{
    size_t len = strlen(s);
    uint64_t h = UINT64_C(0x9E3779B97F4A7C15) ^ (uint64_t)len;
    while (len >= 8) {
        uint64_t w;
        memcpy(&w, s, 8);
        h = (h ^ w) * UINT64_C(0xFF51AFD7ED558CCD);
        h ^= h >> 32;
        s += 8;
        len -= 8;
    }
    if (len) {
        uint64_t w = 0;
        memcpy(&w, s, len);
        h = (h ^ w) * UINT64_C(0xC4CEB9FE1A85EC53);
        h ^= h >> 29;
    }
    h *= UINT64_C(0xFF51AFD7ED558CCD);
    return (uint32_t)(h ^ (h >> 32));
}

/* Capacities are powers of two, so a slot is hash & (capacity - 1). */
static size_t ctools_str_hash_capacity(size_t wanted)
{
    size_t capacity = 16;
    while (capacity < wanted && capacity <= SIZE_MAX / 4) capacity <<= 1;
    return capacity;
}

int ctools_str_hash_init(ctools_str_hash_table *ht, size_t initial_capacity)
{
    if (!ht) return -1;

    /* capacity stays 0 until the allocation succeeds: callers test it to
       decide whether the table is initialised. */
    size_t capacity = ctools_str_hash_capacity(initial_capacity);
    ht->capacity = 0;
    ht->count = 0;
    /* Keys are copied into blocks of the same size as the initial table's
       key space; the blocks are allocated as they are needed. */
    ctools_arena_init(&ht->keys, capacity * 64);
    ht->entries = calloc(capacity, sizeof(ctools_str_hash_entry));
    if (!ht->entries) return -1;

    ht->capacity = capacity;
    return 0;
}

void ctools_str_hash_free(ctools_str_hash_table *ht)
{
    if (!ht) return;

    if (ht->capacity || ht->entries) ctools_arena_free(&ht->keys);
    free(ht->entries);
    ht->entries = NULL;
    ht->count = 0;
    ht->capacity = 0;
}

/* Whether one more key keeps the load factor below 3/4. */
static int ctools_str_hash_room(const ctools_str_hash_table *ht)
{
    return ht->count * 4 < ht->capacity * 3;
}

static int ctools_str_hash_resize(ctools_str_hash_table *ht, size_t new_capacity)
{
    ctools_str_hash_entry *old_entries = ht->entries;
    size_t old_capacity = ht->capacity;

    ctools_str_hash_entry *entries = calloc(new_capacity, sizeof(ctools_str_hash_entry));
    if (!entries) return -1;

    /* Rehash all existing entries; the new table has room for every one. */
    const size_t mask = new_capacity - 1;
    for (size_t i = 0; i < old_capacity; i++) {
        if (old_entries[i].key) {
            size_t idx = old_entries[i].hash & mask;
            while (entries[idx].key) idx = (idx + 1) & mask;
            entries[idx] = old_entries[i];
        }
    }

    ht->entries = entries;
    ht->capacity = new_capacity;
    free(old_entries);
    return 0;
}

/*
 * The slot that holds key, or the empty slot where it belongs. The table
 * always has an empty slot (load factor below 3/4).
 */
static size_t ctools_str_hash_slot(const ctools_str_hash_table *ht, const char *key,
                                   uint32_t hash)
{
    const size_t mask = ht->capacity - 1;
    size_t idx = hash & mask;
    while (ht->entries[idx].key) {
        if (ht->entries[idx].hash == hash && strcmp(ht->entries[idx].key, key) == 0)
            break;
        idx = (idx + 1) & mask;
    }
    return idx;
}

/* Copy key into slot idx with value; returns 0, or -1 when out of memory. */
static int ctools_str_hash_fill(ctools_str_hash_table *ht, size_t idx, const char *key,
                                uint32_t hash, int value)
{
    char *key_copy = ctools_arena_strdup(&ht->keys, key);
    if (!key_copy) return -1;

    ht->entries[idx].key = key_copy;
    ht->entries[idx].hash = hash;
    ht->entries[idx].value = value;
    ht->count++;
    return 0;
}

int ctools_str_hash_insert(ctools_str_hash_table *ht, const char *key,
                            uint32_t precomputed_hash)
{
    if (!ht || !key || !ht->capacity) return -1;

    /* Grow before the load factor reaches 3/4 */
    if (!ctools_str_hash_room(ht) && ctools_str_hash_resize(ht, ht->capacity * 2) != 0)
        return -1;

    size_t idx = ctools_str_hash_slot(ht, key, precomputed_hash);
    if (ht->entries[idx].key) return ht->entries[idx].value; /* existing key */

    /* New key: value = count + 1 (1-based indexing) */
    if (ctools_str_hash_fill(ht, idx, key, precomputed_hash, (int)(ht->count + 1)) != 0)
        return -1;
    return ht->entries[idx].value;
}

int ctools_str_hash_insert_value(ctools_str_hash_table *ht, const char *key,
                                  int value)
{
    if (!ht || !key || !ht->capacity) return -1;

    /* Grow before the load factor reaches 3/4 */
    if (!ctools_str_hash_room(ht) && ctools_str_hash_resize(ht, ht->capacity * 2) != 0)
        return -1;

    uint32_t hash = ctools_str_hash_compute(key);
    size_t idx = ctools_str_hash_slot(ht, key, hash);
    /* An existing mapping is retained; status is independent of its code. */
    if (ht->entries[idx].key) return 0;
    return ctools_str_hash_fill(ht, idx, key, hash, value);
}

int ctools_str_hash_lookup(ctools_str_hash_table *ht, const char *key, int *value)
{
    if (!ht || !key || !value || ht->count == 0) return 0;

    size_t idx = ctools_str_hash_slot(ht, key, ctools_str_hash_compute(key));
    if (!ht->entries[idx].key) return 0; /* Not found */
    *value = ht->entries[idx].value;
    return 1;
}

/* ============================================================================
 * Value Label Output (cencode)
 * ============================================================================ */

/* Transfer data, never Stata/Mata source. Each ASCII code/byte-count header is
 * followed by exactly that many raw bytes, so quotes, macros and newlines need
 * no escaping. Binary mode keeps this framing identical on every platform. */
int ctools_label_write_file(const char **strings, const int *codes, size_t n_labels,
                            const char *filepath)
{
    FILE *fp = ctools_fopen(filepath, "wb");
    if (!fp)
        return -1;
    fprintf(fp, "CTOOLS_LABELS_1\n");
    for (size_t i = 0; i < n_labels; i++)
    {
        const char *text = strings[i] ? strings[i] : "";
        size_t length = strlen(text);
        fprintf(fp, "%d %zu\n", codes[i], length);
        if (fwrite(text, 1, length, fp) != length)
        {
            fclose(fp);
            return -1;
        }
    }
    int failed = ferror(fp);
    if (fclose(fp) != 0)
        failed = 1;
    return failed ? -1 : 0;
}

/* ============================================================================
 * Value Label Input (cencode, cdecode)
 * ============================================================================ */

int ctools_label_map_read_file(const char *path, size_t len, char *buf)
{
    if (!path || !path[0])
        return 198;
    FILE *fp = ctools_fopen(path, "rb");
    if (!fp)
        return 198;
    /* One byte more than announced shows a longer file. */
    size_t got = fread(buf, 1, len + 1, fp);
    int failed = ferror(fp);
    fclose(fp);
    if (failed || got != len)
        return 198;
    buf[len] = '\0';
    return strlen(buf) == len ? 0 : 198;
}
