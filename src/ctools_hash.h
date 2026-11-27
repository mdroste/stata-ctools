/*
 * ctools_hash.h
 * String -> integer hash table and value-label output for ctools
 *
 * The hash table uses:
 * - Open addressing with linear probing over a power-of-two capacity
 * - Automatic resizing at 75% load factor
 * - A growing block arena for the keys
 */

#ifndef CTOOLS_HASH_H
#define CTOOLS_HASH_H

#include <stddef.h>
#include <stdint.h>
#include "ctools_arena.h"

/* Default configuration */
#define CTOOLS_HASH_LOAD_FACTOR  0.75

/* ============================================================================
 * String -> Integer Hash Table (for cencode)
 *
 * Maps unique string values to integer codes.
 * ============================================================================ */

typedef struct {
    char *key;          /* String key (owned by the table's arena) */
    int value;          /* Integer value */
    uint32_t hash;      /* Cached hash for faster probing */
} ctools_str_hash_entry;

typedef struct {
    ctools_str_hash_entry *entries;
    size_t capacity;    /* a power of two; 0 until initialised */
    size_t count;
    ctools_arena keys;  /* copies of the keys */
} ctools_str_hash_table;

/*
 * Hash of a string (word-at-a-time multiply-xorshift).
 */
uint32_t ctools_str_hash_compute(const char *s);

/*
 * Initialize a string->int hash table with room for about
 * initial_capacity * 3/4 keys before it grows (the capacity is rounded up
 * to a power of two).
 * Returns 0 on success, -1 on allocation failure.
 */
int ctools_str_hash_init(ctools_str_hash_table *ht, size_t initial_capacity);

/*
 * Free all memory associated with a string->int hash table.
 * Safe to call on an uninitialized or already-freed table.
 */
void ctools_str_hash_free(ctools_str_hash_table *ht);

/*
 * Insert a string key with auto-assigned value.
 * If key exists, returns the existing value.
 * If key is new, assigns value = count + 1 (1-based indexing).
 *
 * @param ht              Hash table
 * @param key             String key to insert
 * @param precomputed_hash  Hash value from ctools_str_hash_compute()
 * @return                Assigned value (positive), or -1 on error
 */
int ctools_str_hash_insert(ctools_str_hash_table *ht, const char *key,
                            uint32_t precomputed_hash);

/*
 * Insert a string key with a specific value.
 * If key exists, keeps its existing value (does not update).
 * Useful for loading existing label mappings.
 *
 * @param ht     Hash table
 * @param key    String key to insert
 * @param value  Value to associate with the key
 * @return       0 on success, -1 on allocation failure
 */
int ctools_str_hash_insert_value(ctools_str_hash_table *ht, const char *key,
                                  int value);

/*
 * Look up a string key.
 *
 * @param ht   Hash table
 * @param key  String key to find
 * @param value Output: associated signed code when found
 * @return     1 if found, 0 if absent
 */
int ctools_str_hash_lookup(ctools_str_hash_table *ht, const char *key, int *value);

/* ============================================================================
 * Value Label Output (for cencode)
 * ============================================================================ */

#define CTOOLS_MAX_LABELS       65536

/*
 * Write CTOOLS_LABELS_1 data for cencode's Mata label reader: an ASCII magic
 * line, then "code byte_count\n" followed by exactly byte_count raw bytes for
 * each entry. The strings are borrowed; the file contains no executable code.
 * Returns 0 on success, -1 on error.
 */
int ctools_label_write_file(const char **strings, const int *codes, size_t n_labels,
                            const char *filepath);

/* ============================================================================
 * Value Label Input (cencode, cdecode)
 * ============================================================================ */

/*
 * Read a label mapping too long for a macro: exactly len bytes from the
 * file at path into buf (at least len + 2 bytes). Returns 0, or 198 when the
 * file is missing, of another length or holds a NUL byte.
 */
int ctools_label_map_read_file(const char *path, size_t len, char *buf);

#endif /* CTOOLS_HASH_H */
