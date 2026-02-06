/*
 * cmerge_memory.h
 * Memory allocation utilities for cmerge
 *
 * Provides:
 * - Growing arena allocator for efficient bulk string operations
 *
 * Note: For cache-line aligned memory allocation, use the shared functions
 * from ctools_config.h: ctools_cacheline_alloc() and ctools_aligned_free()
 *
 * The cmerge string arena uses the growing arena (ctools_arena) which
 * allocates from a chain of blocks. When a block fills up, a new one is
 * allocated automatically. This handles strings of any length up to
 * Stata's str2045 limit without fallback to per-string malloc.
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
 * Uses the growing arena (ctools_arena) which allocates from a chain of
 * fixed-size blocks. Bump-pointer allocation within each block gives O(1)
 * per-string cost. When a block is exhausted, a new block is added.
 * All blocks are freed together in O(blocks) time.
 *
 * Thread safety: Each string variable gets its own arena, so parallel
 * OpenMP loops have zero contention.
 * ============================================================================ */

/* cmerge_string_arena wraps the growing arena */
typedef ctools_arena cmerge_string_arena;

/* Create a growing string arena with given initial block size.
 * Returns NULL on allocation failure. */
cmerge_string_arena *cmerge_arena_create(size_t block_size);

/* Allocate a string copy from the arena.
 * The arena grows automatically if the current block is full.
 * Returns NULL only on malloc failure (new block allocation). */
char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s);

/* Check if a pointer was allocated from this arena */
int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr);

/* Free the arena and all its blocks */
void cmerge_arena_free(cmerge_string_arena *arena);

#endif /* CMERGE_MEMORY_H */
