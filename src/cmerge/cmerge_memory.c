/*
 * cmerge_memory.c
 * Memory allocation utilities for cmerge
 *
 * Wrappers around ctools_arena (growing arena) for merge string output.
 * The growing arena allocates from blocks and grows automatically,
 * handling strings of any length without fallback to per-string malloc.
 */

#include <stdlib.h>
#include <string.h>
#include "cmerge_memory.h"

/* ============================================================================
 * String Arena Allocator - Wrappers around ctools_arena (growing arena)
 * ============================================================================ */

cmerge_string_arena *cmerge_arena_create(size_t block_size)
{
    cmerge_string_arena *arena = (cmerge_string_arena *)malloc(sizeof(cmerge_string_arena));
    if (!arena) return NULL;

    /* Minimum 4KB blocks to avoid excessive block allocation overhead */
    if (block_size < 4096) block_size = 4096;

    ctools_arena_init(arena, block_size);
    return arena;
}

char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s)
{
    if (!arena) {
        /* Fallback if arena creation failed */
        return s ? strdup(s) : NULL;
    }
    return ctools_arena_strdup(arena, s);
}

int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr)
{
    if (!arena || !ptr) return 0;

    /* Check all blocks in the chain */
    CToolsArenaBlock *block = arena->first;
    while (block) {
        if (ptr >= block->data && ptr < block->data + block->capacity) return 1;
        block = block->next;
    }
    return 0;
}

void cmerge_arena_free(cmerge_string_arena *arena)
{
    if (arena) {
        ctools_arena_free(arena);  /* Free all blocks */
        free(arena);               /* Free the struct itself */
    }
}
