/*
 * ctools_hash.h
 * Hash tables and label utilities for ctools
 *
 * Provides:
 * 1. ctools_str_hash_table - String keys to integer values (for cencode)
 * 2. ctools_int_hash_table - Integer keys to string values (for cdecode)
 * 3. Label escape/unescape/parse/serialize utilities (for cdecode/cencode)
 *
 * Both hash table implementations use:
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

/* ============================================================================
 * Label Utilities
 *
 * Shared label escape/unescape/parse/serialize for cdecode and cencode.
 *
 * Escape convention (for serialization):
 *   " -> \"    \ -> \\    | -> \|
 *
 * Unescape convention (for parsing):
 *   \" -> "    \\ -> \    \| -> |
 * ============================================================================ */

#define CTOOLS_MAX_LABEL_LEN    2045
#define CTOOLS_MAX_LABELS       65536
#define CTOOLS_LABEL_BUF_SIZE   2048
#define CTOOLS_MACRO_CHUNK_SIZE 31000

/*
 * Unescape a label string.
 * Converts \| -> |, \\ -> \, \" -> "
 * Writes at most max_len characters (excluding null terminator) to dst.
 */
void ctools_label_unescape(const char *src, char *dst, size_t max_len);

/*
 * Escape a label string for serialization.
 * Escapes: " -> \", \ -> \\, | -> \|
 * Writes at most dst_size-1 characters to dst (always null-terminated).
 * Returns number of bytes written (excluding null terminator).
 */
size_t ctools_label_escape(const char *src, char *dst, size_t dst_size);

/*
 * Parse value-label pairs from a file (one "value|label" per line).
 * Used by cdecode. Calls ctools_label_unescape internally.
 * If max_label_len is non-NULL, tracks the maximum unescaped label length.
 * Returns 0 on success, -1 on error.
 */
int ctools_label_parse_file(const char *filepath,
                            ctools_int_hash_table *ht,
                            int *max_label_len);

/*
 * Parse value-label pairs from a "value|label||value|label||..." string.
 * Used by cencode (noextend path). Calls ctools_label_unescape internally.
 * Returns 0 on success, -1 on error.
 */
int ctools_label_parse_string(const char *str,
                              ctools_str_hash_table *ht);

/*
 * Serialize sorted label entries into chunked Stata macros.
 * Used by cencode to pass new labels back to the .ado file.
 * Calls ctools_label_escape internally.
 *
 * Format per entry: "code|escaped_string", entries separated by "||"
 * Chunks are saved as Stata macros named "<prefix>_0", "<prefix>_1", etc.
 * Each chunk is at most CTOOLS_MACRO_CHUNK_SIZE bytes.
 *
 * prefix: macro name prefix (e.g. "_cencode_labels")
 * Returns number of chunks written, or -1 on error.
 */
int ctools_label_serialize_macros(const char **strings, const int *codes,
                                  size_t n_labels, const char *prefix);

/*
 * Parse a Stata `label save` format file directly.
 * Lines have format: label define <name> <value> `"text"', modify
 *
 * These functions parse label save files, eliminating the need
 * for intermediate file formats or Mata preprocessing.
 */

/*
 * Parse into int->string hash table (for cdecode).
 * If max_label_len is non-NULL, tracks the maximum label length.
 * Returns 0 on success, -1 on error.
 */
int ctools_label_parse_stata_file_int(const char *filepath,
                                       ctools_int_hash_table *ht,
                                       int *max_label_len);

/*
 * Parse into string->int hash table (for cencode noextend).
 * Returns 0 on success, -1 on error.
 */
int ctools_label_parse_stata_file_str(const char *filepath,
                                       ctools_str_hash_table *ht);

/*
 * Scan for max label length only, without building a hash table.
 * Used by cdecode_scan for pre-flight maxlen detection.
 * Returns 0 on success, -1 on error.
 */
int ctools_label_scan_stata_file(const char *filepath, int *max_label_len);

/*
 * Write label entries as a Stata .do file with `label define` commands.
 * Used by cencode to pass new labels back to the .ado file.
 * Format: label define <name> <code> `"text"'[, add]
 * Returns 0 on success, -1 on error.
 */
int ctools_label_write_stata_file(const char **strings, const int *codes,
                                   size_t n_labels, const char *label_name,
                                   const char *filepath);

#endif /* CTOOLS_HASH_H */
