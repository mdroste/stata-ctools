/*
 * ctools_hash.h
 * String -> integer hash table and value-label output for ctools
 *
 * The hash table uses:
 * - Open addressing with linear probing
 * - Automatic resizing at 75% load factor
 * - String arena for efficient memory management
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
    char *key;          /* String key (owned by arena) */
    int value;          /* Integer value */
    uint32_t hash;      /* Cached hash for faster probing */
} ctools_str_hash_entry;

typedef struct {
    ctools_str_hash_entry *entries;
    size_t capacity;
    size_t count;
    ctools_string_arena *arena;
} ctools_str_hash_table;

/*
 * Compute FNV-1a hash of a string.
 * Fast, good distribution for typical string data.
 */
uint32_t ctools_str_hash_compute(const char *s);

/*
 * Initialize a string->int hash table.
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
 * Write label entries as a .do file that rebuilds each label with Mata
 * st_vlmodify() from byte literals. Used by cencode to pass new labels back
 * to the .ado file.
 * Returns 0 on success, -1 on error.
 */
int ctools_label_write_stata_file(const char **strings, const int *codes,
                                   size_t n_labels, const char *label_name,
                                   const char *filepath);

#endif /* CTOOLS_HASH_H */
