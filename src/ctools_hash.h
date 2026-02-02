/*
 * ctools_hash.h
 * Unified hash table implementations for ctools
 *
 * Provides two specialized hash table types:
 * 1. ctools_str_hash_table - String keys to integer values (for cencode)
 * 2. ctools_int_hash_table - Integer keys to string values (for cdecode)
 *
 * Both implementations use:
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
#define CTOOLS_HASH_INIT_SIZE    1024
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
 * If key exists, returns the existing value (does not update).
 * Useful for loading existing label mappings.
 *
 * @param ht     Hash table
 * @param key    String key to insert
 * @param value  Value to associate with the key
 * @return       The value (new or existing), or -1 on error
 */
int ctools_str_hash_insert_value(ctools_str_hash_table *ht, const char *key,
                                  int value);

/*
 * Look up a string key.
 *
 * @param ht   Hash table
 * @param key  String key to find
 * @return     Associated value if found, 0 if not found
 */
int ctools_str_hash_lookup(ctools_str_hash_table *ht, const char *key);

/* ============================================================================
 * Integer -> String Hash Table (for cdecode)
 *
 * Maps integer values to string labels.
 * ============================================================================ */

typedef struct {
    int key;            /* Integer key */
    char *value;        /* String value (owned by arena) */
    int occupied;       /* 1 if slot is used, 0 otherwise */
} ctools_int_hash_entry;

typedef struct {
    ctools_int_hash_entry *entries;
    size_t capacity;
    size_t count;
    ctools_string_arena *arena;
} ctools_int_hash_table;

/*
 * Compute hash of an integer key.
 * Uses bit mixing for good distribution.
 */
uint32_t ctools_int_hash_compute(int key);

/*
 * Initialize an int->string hash table.
 * Returns 0 on success, -1 on allocation failure.
 */
int ctools_int_hash_init(ctools_int_hash_table *ht, size_t initial_capacity);

/*
 * Free all memory associated with an int->string hash table.
 * Safe to call on an uninitialized or already-freed table.
 */
void ctools_int_hash_free(ctools_int_hash_table *ht);

/*
 * Insert an integer key with a string value.
 * If key exists, does nothing (keeps original value).
 *
 * @param ht     Hash table
 * @param key    Integer key
 * @param value  String value to associate
 * @return       0 on success, -1 on error
 */
int ctools_int_hash_insert(ctools_int_hash_table *ht, int key, const char *value);

/*
 * Look up an integer key.
 *
 * @param ht   Hash table
 * @param key  Integer key to find
 * @return     Associated string value if found, NULL if not found
 */
const char *ctools_int_hash_lookup(ctools_int_hash_table *ht, int key);

#endif /* CTOOLS_HASH_H */
