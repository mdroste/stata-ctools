/*
 * cmerge_memory.h
 * Memory allocation utilities for cmerge
 *
 * Provides:
 * - Cache-line aligned memory allocation
 * - String arena allocator for efficient bulk string operations
 */

#ifndef CMERGE_MEMORY_H
#define CMERGE_MEMORY_H

#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Cache-Line Aligned Memory Allocation
 * ============================================================================ */

/* Allocate memory aligned to cache line boundary for optimal memory access */
void *cmerge_aligned_alloc(size_t size);

/* Free cache-line aligned memory */
void cmerge_aligned_free(void *ptr);

/* ============================================================================
 * String Arena Allocator
 *
 * Allocates strings from a contiguous block to reduce per-string malloc
 * overhead and enable O(1) bulk free instead of O(n) individual frees.
 * ============================================================================ */

typedef struct {
    char *base;         /* Base pointer to arena memory */
    size_t capacity;    /* Total arena capacity in bytes */
    size_t used;        /* Currently used bytes */
    int alloc_failed;   /* Set to 1 if any allocation failed */
} cmerge_string_arena;

/* Create a string arena with given capacity */
cmerge_string_arena *cmerge_arena_create(size_t capacity);

/* Allocate a string copy from the arena. Falls back to strdup if full.
 * Sets arena->alloc_failed on allocation failure.
 * ALWAYS returns a valid pointer (never NULL) - falls back to static empty string. */
char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s);

/* Check if a pointer was allocated from the arena (or is the static fallback) */
int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr);

/* Free the arena. Note: strings from fallback strdup must be freed separately */
void cmerge_arena_free(cmerge_string_arena *arena);

#endif /* CMERGE_MEMORY_H */
