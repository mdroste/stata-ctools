/*
 * ctools_arena.c
 * Unified arena allocator implementation for ctools
 *
 * See ctools_arena.h for documentation.
 */

#include <stdlib.h>
#include <string.h>
#include "ctools_arena.h"

/* ============================================================================
 * Static empty string for STATIC_FALLBACK mode
 *
 * Using volatile to prevent compiler from optimizing away pointer comparisons.
 * This is used as a last-resort fallback that is never freed.
 * ============================================================================ */

char ctools_static_empty_string[] = "";

/* ============================================================================
 * Growing Arena Implementation
 * ============================================================================ */

/* Create a new block with given capacity */
static CToolsArenaBlock *ctools_arena_new_block(size_t capacity)
{
    CToolsArenaBlock *block = (CToolsArenaBlock *)malloc(
        sizeof(CToolsArenaBlock) + capacity);
    if (block == NULL) {
        return NULL;
    }
    block->next = NULL;
    block->used = 0;
    block->capacity = capacity;
    return block;
}

void ctools_arena_init(ctools_arena *arena, size_t block_size)
{
    if (arena == NULL) return;

    if (block_size == 0) {
        block_size = CTOOLS_ARENA_DEFAULT_BLOCK_SIZE;
    }

    arena->first = NULL;
    arena->current = NULL;
    arena->block_size = block_size;
    arena->total_allocated = 0;
}

void *ctools_arena_alloc(ctools_arena *arena, size_t size)
{
    if (arena == NULL || size == 0) {
        return NULL;
    }

    /* Align to 8 bytes */
    size = (size + 7) & ~(size_t)7;

    /* Try current block first */
    if (arena->current != NULL &&
        arena->current->used + size <= arena->current->capacity) {
        void *ptr = arena->current->data + arena->current->used;
        arena->current->used += size;
        arena->total_allocated += size;
        return ptr;
    }

    /* Need a new block - use larger of block_size and requested size */
    size_t new_capacity = arena->block_size;
    if (size > new_capacity) {
        new_capacity = size;
    }

    CToolsArenaBlock *new_block = ctools_arena_new_block(new_capacity);
    if (new_block == NULL) {
        return NULL;
    }

    /* Link to chain */
    if (arena->first == NULL) {
        arena->first = new_block;
    } else if (arena->current != NULL) {
        arena->current->next = new_block;
    }
    arena->current = new_block;

    /* Allocate from new block */
    void *ptr = new_block->data;
    new_block->used = size;
    arena->total_allocated += size;
    return ptr;
}

char *ctools_arena_strdup(ctools_arena *arena, const char *s)
{
    if (s == NULL) return NULL;

    size_t len = strlen(s) + 1;
    char *copy = (char *)ctools_arena_alloc(arena, len);
    if (copy != NULL) {
        memcpy(copy, s, len);
    }
    return copy;
}

void ctools_arena_reset(ctools_arena *arena)
{
    if (arena == NULL) return;

    /* Reset all blocks to empty but keep them allocated */
    CToolsArenaBlock *block = arena->first;
    while (block != NULL) {
        block->used = 0;
        block = block->next;
    }
    arena->current = arena->first;
    arena->total_allocated = 0;
}

void ctools_arena_free(ctools_arena *arena)
{
    if (arena == NULL) return;

    CToolsArenaBlock *block = arena->first;
    while (block != NULL) {
        CToolsArenaBlock *next = block->next;
        free(block);
        block = next;
    }

    arena->first = NULL;
    arena->current = NULL;
    arena->total_allocated = 0;
}

/* ============================================================================
 * String Arena Implementation
 * ============================================================================ */

ctools_string_arena *ctools_string_arena_create(size_t capacity, ctools_string_arena_mode mode)
{
    ctools_string_arena *arena = (ctools_string_arena *)malloc(sizeof(ctools_string_arena));
    if (arena == NULL) {
        return NULL;
    }

    arena->base = (char *)malloc(capacity);
    if (arena->base == NULL) {
        free(arena);
        return NULL;
    }

    arena->capacity = capacity;
    arena->used = 0;
    arena->mode = mode;
    arena->has_fallback = 0;

    return arena;
}

char *ctools_string_arena_strdup(ctools_string_arena *arena, const char *s)
{
    if (s == NULL) {
        return NULL;
    }

    size_t len = strlen(s) + 1;

    /* Try to allocate from arena (with overflow protection) */
    if (arena != NULL &&
        arena->used <= arena->capacity - len &&
        len <= arena->capacity) {
        char *ptr = arena->base + arena->used;
        memcpy(ptr, s, len);
        arena->used += len;
        return ptr;
    }

    /* Arena full - handle based on mode */
    if (arena == NULL) {
        /* No arena - just use strdup */
        return strdup(s);
    }

    switch (arena->mode) {
        case CTOOLS_STRING_ARENA_NO_FALLBACK:
            return NULL;

        case CTOOLS_STRING_ARENA_STRDUP_FALLBACK: {
            char *result = strdup(s);
            if (result != NULL) {
                arena->has_fallback = 1;
            }
            return result;
        }

        case CTOOLS_STRING_ARENA_STATIC_FALLBACK:
            arena->has_fallback = 1;
            /* Try to allocate empty string in arena first */
            if (arena->used < arena->capacity) {
                char *ptr = arena->base + arena->used;
                ptr[0] = '\0';
                arena->used += 1;
                return ptr;
            }
            /* Last resort: return static empty string */
            return ctools_static_empty_string;
    }

    /* Should not reach here */
    return NULL;
}

int ctools_string_arena_owns(ctools_string_arena *arena, const char *ptr)
{
    if (ptr == NULL) {
        return 0;
    }

    /* Check if pointer is the static empty string */
    if (ptr == ctools_static_empty_string) {
        return 1;  /* Treat as owned (should not be freed) */
    }

    /* Check if pointer is within arena bounds */
    if (arena != NULL) {
        return (ptr >= arena->base && ptr < arena->base + arena->capacity);
    }

    return 0;
}

void ctools_string_arena_free(ctools_string_arena *arena)
{
    if (arena != NULL) {
        free(arena->base);
        free(arena);
    }
}
