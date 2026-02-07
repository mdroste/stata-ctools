/*
 * ctools_hash.c
 * Hash tables and label utilities for ctools
 *
 * Extracted from cencode_impl.c and cdecode_impl.c for code reuse.
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

/* ============================================================================
 * Label Escape / Unescape
 * ============================================================================ */

void ctools_label_unescape(const char *src, char *dst, size_t max_len)
{
    size_t j = 0;
    for (size_t i = 0; src[i] && j < max_len; i++) {
        if (src[i] == '\\' && src[i + 1]) {
            i++;
            if (src[i] == '|') {
                dst[j++] = '|';
            } else if (src[i] == '\\') {
                dst[j++] = '\\';
            } else if (src[i] == '"') {
                dst[j++] = '"';
            } else {
                /* Unknown escape — pass through the escaped char */
                dst[j++] = src[i];
            }
        } else {
            dst[j++] = src[i];
        }
    }
    dst[j] = '\0';
}

size_t ctools_label_escape(const char *src, char *dst, size_t dst_size)
{
    size_t pos = 0;
    for (const char *p = src; *p && pos < dst_size - 2; p++) {
        if (*p == '"' || *p == '\\' || *p == '|') {
            dst[pos++] = '\\';
        }
        dst[pos++] = *p;
    }
    dst[pos] = '\0';
    return pos;
}

/* ============================================================================
 * Label File Parsing (cdecode)
 * ============================================================================ */

int ctools_label_parse_file(const char *filepath,
                            ctools_int_hash_table *ht,
                            int *max_label_len)
{
    if (!filepath || *filepath == '\0') {
        return 0;  /* No file is OK */
    }

    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        return -1;
    }

    if (max_label_len) {
        *max_label_len = 0;
    }

    char line[CTOOLS_MAX_LABEL_LEN + 32];

    while (fgets(line, sizeof(line), fp)) {
        /* Remove trailing newline */
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }

        if (len == 0) continue;

        /* Find the pipe separator between value and label */
        char *pipe = strchr(line, '|');
        if (!pipe) continue;

        *pipe = '\0';
        int value;
        if (!ctools_safe_atoi(line, &value)) {
            fclose(fp);
            return -1;  /* Invalid value in labels file */
        }
        char *label = pipe + 1;

        /* Unescape the label */
        char unescaped[CTOOLS_MAX_LABEL_LEN + 1];
        ctools_label_unescape(label, unescaped, CTOOLS_MAX_LABEL_LEN);

        /* Track max length */
        if (max_label_len) {
            int lbl_len = (int)strlen(unescaped);
            if (lbl_len > *max_label_len) {
                *max_label_len = lbl_len;
            }
        }

        /* Insert into hash table */
        if (ctools_int_hash_insert(ht, value, unescaped) != 0) {
            fclose(fp);
            return -1;
        }
    }

    fclose(fp);
    return 0;
}

/* ============================================================================
 * Label String Parsing (cencode noextend path)
 * ============================================================================ */

int ctools_label_parse_string(const char *str,
                              ctools_str_hash_table *ht)
{
    if (!str || *str == '\0') return 0;

    size_t len = strlen(str);
    char *buf = malloc(len + 1);
    if (!buf) return -1;
    memcpy(buf, str, len + 1);

    char *pos = buf;
    while (*pos) {
        /* Parse value */
        char *pipe = strchr(pos, '|');
        if (!pipe) break;
        *pipe = '\0';
        int value;
        if (!ctools_safe_atoi(pos, &value)) {
            free(buf);
            return -1;
        }
        pos = pipe + 1;

        /* Find end of string (|| or end of buffer) */
        char *end = strstr(pos, "||");
        if (end) {
            *end = '\0';
        }

        /* Unescape the string in-place */
        char unescaped[CTOOLS_MAX_LABEL_LEN + 1];
        ctools_label_unescape(pos, unescaped, CTOOLS_MAX_LABEL_LEN);

        /* Insert into hash table */
        if (ctools_str_hash_insert_value(ht, unescaped, value) < 0) {
            free(buf);
            return -1;
        }

        if (end) {
            pos = end + 2;  /* Skip || */
        } else {
            break;
        }
    }

    free(buf);
    return 0;
}

/* ============================================================================
 * Label Serialization (cencode)
 * ============================================================================ */

int ctools_label_serialize_macros(const char **strings, const int *codes,
                                  size_t n_labels, const char *prefix)
{
    char *macro_buf = malloc(CTOOLS_MACRO_CHUNK_SIZE + CTOOLS_LABEL_BUF_SIZE);
    if (!macro_buf) return -1;

    macro_buf[0] = '\0';
    size_t macro_pos = 0;
    size_t chunk = 0;
    size_t labels_in_chunk = 0;
    char macro_name[64];

    for (size_t i = 0; i < n_labels; i++) {
        int code = codes[i];
        const char *str = strings[i];

        char escaped[CTOOLS_LABEL_BUF_SIZE];
        ctools_label_escape(str, escaped, sizeof(escaped));

        char entry[CTOOLS_LABEL_BUF_SIZE + 32];
        int entry_len = snprintf(entry, sizeof(entry), "%s%d|%s",
                                 (labels_in_chunk > 0) ? "||" : "",
                                 code, escaped);

        if (macro_pos + (size_t)entry_len >= CTOOLS_MACRO_CHUNK_SIZE) {
            snprintf(macro_name, sizeof(macro_name), "%s_%zu", prefix, chunk);
            SF_macro_save(macro_name, macro_buf);

            chunk++;
            labels_in_chunk = 0;
            macro_pos = 0;
            macro_buf[0] = '\0';

            entry_len = snprintf(entry, sizeof(entry), "%d|%s", code, escaped);
        }

        memcpy(macro_buf + macro_pos, entry, (size_t)entry_len + 1);
        macro_pos += (size_t)entry_len;
        labels_in_chunk++;
    }

    if (labels_in_chunk > 0) {
        snprintf(macro_name, sizeof(macro_name), "%s_%zu", prefix, chunk);
        SF_macro_save(macro_name, macro_buf);
    }

    free(macro_buf);
    return (int)(chunk + 1);
}

/* ============================================================================
 * Stata `label save` File Parsing
 *
 * Parses files produced by Stata's `label save` command, which have lines like:
 *   label define lblname 1 `"Label Text"', modify
 *   label define lblname 2 "Simple", modify
 * ============================================================================ */

/*
 * Parse a single line from a Stata `label save` format file.
 * Returns 1 if successfully parsed, 0 if line doesn't match.
 */
static int ctools_parse_label_save_line(const char *line, int *out_value,
                                         char *out_label, int *out_label_len)
{
    /* Must start with "label define " */
    if (strncmp(line, "label define ", 13) != 0) return 0;
    const char *p = line + 13;

    /* Skip label name */
    while (*p && *p != ' ') p++;
    if (!*p) return 0;
    while (*p == ' ') p++;

    /* Parse integer value */
    char *endptr;
    long val = strtol(p, &endptr, 10);
    if (endptr == p) return 0;
    *out_value = (int)val;

    p = endptr;
    while (*p == ' ') p++;

    /* Try compound quotes: `"..."' */
    if (p[0] == '`' && p[1] == '"') {
        p += 2;
        const char *close = strstr(p, "\"'");
        if (!close) return 0;

        size_t len = (size_t)(close - p);
        if (len > CTOOLS_MAX_LABEL_LEN) len = CTOOLS_MAX_LABEL_LEN;
        memcpy(out_label, p, len);
        out_label[len] = '\0';
        *out_label_len = (int)len;
        return 1;
    }

    /* Try simple quotes: "..." */
    if (*p == '"') {
        p++;
        const char *close = strchr(p, '"');
        if (!close) return 0;

        size_t len = (size_t)(close - p);
        if (len > CTOOLS_MAX_LABEL_LEN) len = CTOOLS_MAX_LABEL_LEN;
        memcpy(out_label, p, len);
        out_label[len] = '\0';
        *out_label_len = (int)len;
        return 1;
    }

    return 0;
}

int ctools_label_parse_stata_file_int(const char *filepath,
                                       ctools_int_hash_table *ht,
                                       int *max_label_len)
{
    if (!filepath || *filepath == '\0') return 0;

    FILE *fp = fopen(filepath, "r");
    if (!fp) return -1;

    if (max_label_len) *max_label_len = 0;

    char line[CTOOLS_MAX_LABEL_LEN + 256];
    char label[CTOOLS_MAX_LABEL_LEN + 1];
    int value, label_len;

    while (fgets(line, sizeof(line), fp)) {
        /* Strip trailing newline */
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';

        if (ctools_parse_label_save_line(line, &value, label, &label_len)) {
            if (ctools_int_hash_insert(ht, value, label) != 0) {
                fclose(fp);
                return -1;
            }
            if (max_label_len && label_len > *max_label_len)
                *max_label_len = label_len;
        }
    }

    fclose(fp);
    return 0;
}

int ctools_label_parse_stata_file_str(const char *filepath,
                                       ctools_str_hash_table *ht)
{
    if (!filepath || *filepath == '\0') return 0;

    FILE *fp = fopen(filepath, "r");
    if (!fp) return -1;

    char line[CTOOLS_MAX_LABEL_LEN + 256];
    char label[CTOOLS_MAX_LABEL_LEN + 1];
    int value, label_len;

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';

        if (ctools_parse_label_save_line(line, &value, label, &label_len)) {
            if (ctools_str_hash_insert_value(ht, label, value) < 0) {
                fclose(fp);
                return -1;
            }
        }
    }

    fclose(fp);
    return 0;
}

int ctools_label_scan_stata_file(const char *filepath, int *max_label_len)
{
    if (!filepath || *filepath == '\0') return 0;

    FILE *fp = fopen(filepath, "r");
    if (!fp) return -1;

    if (max_label_len) *max_label_len = 0;

    char line[CTOOLS_MAX_LABEL_LEN + 256];
    char label[CTOOLS_MAX_LABEL_LEN + 1];
    int value, label_len;

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';

        if (ctools_parse_label_save_line(line, &value, label, &label_len)) {
            if (max_label_len && label_len > *max_label_len)
                *max_label_len = label_len;
        }
    }

    fclose(fp);
    return 0;
}

/* Check if a string contains control characters that Stata's .do parser
 * can't handle (tab expands to spaces, newline breaks the line) */
static int label_string_needs_escape(const char *s)
{
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (*p < 32) return 1;
    }
    return 0;
}

/* Write a string expression using char() for control characters.
 * Output: `"seg1"' + char(9) + `"seg2"' + char(10) + `"seg3"'
 * For use in: local __lbl_tmp = <this expression> */
static void write_escaped_string_expr(FILE *fp, const char *s)
{
    const unsigned char *p = (const unsigned char *)s;
    const unsigned char *seg_start = p;
    int first = 1;

    for (;; p++) {
        if (*p == '\0' || *p < 32) {
            /* Write accumulated safe segment */
            if (p > seg_start) {
                if (!first) fprintf(fp, " + ");
                fprintf(fp, "`\"");
                fwrite(seg_start, 1, (size_t)(p - seg_start), fp);
                fprintf(fp, "\"'");
                first = 0;
            }
            if (*p == '\0') break;
            /* Write char(N) for the control character */
            if (!first) fprintf(fp, " + ");
            fprintf(fp, "char(%d)", (int)*p);
            first = 0;
            seg_start = p + 1;
        }
    }
    if (first) {
        /* Empty string */
        fprintf(fp, "`\"\"'");
    }
}

int ctools_label_write_stata_file(const char **strings, const int *codes,
                                   size_t n_labels, const char *label_name,
                                   const char *filepath)
{
    FILE *fp = fopen(filepath, "w");
    if (!fp) return -1;

    for (size_t i = 0; i < n_labels; i++) {
        const char *add = (i == 0) ? "" : ", add";

        if (!label_string_needs_escape(strings[i])) {
            /* Simple case: no control characters */
            fprintf(fp, "label define %s %d `\"%s\"'%s\n",
                    label_name, codes[i], strings[i], add);
        } else {
            /* Escaped case: use local + char() to preserve control chars */
            fprintf(fp, "local __lbl_tmp = ");
            write_escaped_string_expr(fp, strings[i]);
            fprintf(fp, "\n");
            fprintf(fp, "label define %s %d `\"`__lbl_tmp'\"'%s\n",
                    label_name, codes[i], add);
        }
    }

    fclose(fp);
    return 0;
}
