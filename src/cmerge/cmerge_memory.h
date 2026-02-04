/*
 * cmerge_memory.h
 * Memory allocation utilities for cmerge
 *
 * Provides:
 * - String arena allocator for efficient bulk string operations
 *
 * Note: For cache-line aligned memory allocation, use the shared functions
 * from ctools_config.h: ctools_cacheline_alloc() and ctools_aligned_free()
 *
 * The cmerge string arena is now a typedef to the shared ctools_string_arena.
 * The cmerge_arena_* functions are thin wrappers that configure the arena
 * with STRDUP_FALLBACK mode (falls back to strdup when arena is full).
 */

#ifndef CMERGE_MEMORY_H
#define CMERGE_MEMORY_H

#include <stdlib.h>
#include <string.h>
#include "../ctools_config.h"
#include "../ctools_arena.h"

/* ============================================================================
 * String Arena Allocator
 *
 * Allocates strings from a contiguous block to reduce per-string malloc
 * overhead and enable O(1) bulk free instead of O(n) individual frees.
 *
 * The cmerge arena uses STATIC_FALLBACK mode, meaning it ALWAYS returns a
 * valid pointer (never NULL) - falls back to static empty string.
 * ============================================================================ */

/* cmerge_string_arena is now a typedef to the shared implementation */
typedef ctools_string_arena cmerge_string_arena;

/* Create a string arena with given capacity.
 * Uses STATIC_FALLBACK mode - strdup never returns NULL. */
cmerge_string_arena *cmerge_arena_create(size_t capacity);

/* Allocate a string copy from the arena.
 * Falls back to static empty string when arena is full.
 * ALWAYS returns a valid pointer (never NULL). */
char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s);

/* Check if a pointer was allocated from the arena (or is the static fallback) */
int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr);

/* Free the arena */
void cmerge_arena_free(cmerge_string_arena *arena);

#endif /* CMERGE_MEMORY_H */
