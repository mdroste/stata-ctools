/*
 * ctools_arena.h
 * Unified arena allocator for ctools
 *
 * Provides two arena types:
 * 1. Growing Arena (ctools_arena) - Multi-block, for unknown sizes
 * 2. String Arena (ctools_string_arena) - Single block with optional fallback
 *
 * Usage:
 * - Use ctools_arena for bulk allocations where size is unknown upfront
 *   (e.g., parsing CSV data with unknown row counts)
 * - Use ctools_string_arena for string pooling with known/estimated counts
 *   (e.g., loading Stata string variables)
 */

#ifndef CTOOLS_ARENA_H
#define CTOOLS_ARENA_H

#include <stdlib.h>
#include <stddef.h>
#include <string.h>

/* Default block size for growing arena (1MB) */
#define CTOOLS_ARENA_DEFAULT_BLOCK_SIZE (1024 * 1024)

/* ============================================================================
 * Type 1: Growing Arena (multi-block, for unknown sizes)
 *
 * Allocates from a chain of fixed-size blocks. When the current block is
 * exhausted, a new block is added. All blocks are freed together in O(n)
 * where n is the number of blocks (typically small).
 * ============================================================================ */

typedef struct CToolsArenaBlock {
    struct CToolsArenaBlock *next;  /* Next block in chain */
    size_t used;                    /* Bytes used in this block */
    size_t capacity;                /* Total capacity of this block */
    char data[];                    /* Flexible array member for data */
} CToolsArenaBlock;

typedef struct {
    CToolsArenaBlock *first;        /* First block in chain */
    CToolsArenaBlock *current;      /* Current block for allocations */
    size_t block_size;              /* Size of new blocks */
    size_t total_allocated;         /* Total bytes allocated (for stats) */
} ctools_arena;

/* Initialize an arena with specified block size */
void ctools_arena_init(ctools_arena *arena, size_t block_size);

/* Allocate memory from the arena (8-byte aligned)
 * Returns NULL on allocation failure */
void *ctools_arena_alloc(ctools_arena *arena, size_t size);

/* Duplicate a string into the arena */
char *ctools_arena_strdup(ctools_arena *arena, const char *s);

/* Reset arena (keep blocks but mark as empty) */
void ctools_arena_reset(ctools_arena *arena);

/* Free all memory in the arena */
void ctools_arena_free(ctools_arena *arena);

/* ============================================================================
 * Type 2: String Arena (single block with optional fallback)
 *
 * Optimized for string pooling where an upper bound on total size is known
 * or can be estimated. Supports three fallback modes when the arena is full:
 *
 * - NO_FALLBACK: Return NULL when arena is full
 * - STRDUP_FALLBACK: Fall back to strdup() (caller must track and free)
 * - STATIC_FALLBACK: Return a static empty string (never NULL, never freed)
 *
 * This is more memory-efficient than the growing arena for string-heavy
 * workloads because it uses a single contiguous block.
 * ============================================================================ */

typedef enum {
    CTOOLS_STRING_ARENA_NO_FALLBACK,    /* Return NULL when full */
    CTOOLS_STRING_ARENA_STRDUP_FALLBACK, /* Fall back to strdup */
    CTOOLS_STRING_ARENA_STATIC_FALLBACK  /* Fall back to static empty string */
} ctools_string_arena_mode;

typedef struct {
    char *base;                         /* Base pointer to arena memory */
    size_t capacity;                    /* Total arena capacity in bytes */
    size_t used;                        /* Currently used bytes */
    ctools_string_arena_mode mode;      /* Fallback mode */
    int has_fallback;                   /* Non-zero if any fallback allocations were made */
} ctools_string_arena;

/* Static empty string for STATIC_FALLBACK mode - defined in ctools_arena.c */
extern char ctools_static_empty_string[];

/* Create a string arena with given capacity and fallback mode
 * Returns NULL on allocation failure */
ctools_string_arena *ctools_string_arena_create(size_t capacity, ctools_string_arena_mode mode);

/* Duplicate a string into the arena
 * Behavior when full depends on mode:
 * - NO_FALLBACK: Returns NULL
 * - STRDUP_FALLBACK: Uses strdup (sets has_fallback flag)
 * - STATIC_FALLBACK: Returns static empty string (never NULL) */
char *ctools_string_arena_strdup(ctools_string_arena *arena, const char *s);

/* Check if a pointer was allocated from the arena (or is the static fallback)
 * Returns non-zero if owned by arena or is the static fallback */
int ctools_string_arena_owns(ctools_string_arena *arena, const char *ptr);

/* Free the arena
 * Note: Strings from STRDUP_FALLBACK must be freed separately by caller */
void ctools_string_arena_free(ctools_string_arena *arena);

#endif /* CTOOLS_ARENA_H */
