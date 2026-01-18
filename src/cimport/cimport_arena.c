/*
 * cimport_arena.c
 * Arena allocator for cimport CSV parser
 */

#include <stdlib.h>
#include <limits.h>
#include "cimport_arena.h"

void cimport_arena_init(CImportArena *arena)
{
    arena->first = NULL;
    arena->current = NULL;
    arena->total_allocated = 0;
}

void *cimport_arena_alloc(CImportArena *arena, size_t size)
{
    /* Align to 8 bytes */
    size = (size + 7) & ~7;

    /* Check for overflow in used + size calculation */
    size_t needed = 0;
    if (arena->current) {
        if (size > SIZE_MAX - arena->current->used) return NULL;  /* overflow */
        needed = arena->current->used + size;
    }

    if (!arena->current || needed > arena->current->capacity) {
        size_t block_size = CTOOLS_ARENA_BLOCK_SIZE;
        size_t min_data_size = block_size - sizeof(CImportArenaBlock);
        if (size > min_data_size) {
            /* Check for overflow in size + sizeof addition */
            if (size > SIZE_MAX - sizeof(CImportArenaBlock)) return NULL;
            block_size = size + sizeof(CImportArenaBlock);
        }

        CImportArenaBlock *block = malloc(block_size);
        if (!block) return NULL;

        block->next = NULL;
        block->used = 0;
        block->capacity = block_size - sizeof(CImportArenaBlock);

        if (arena->current) {
            arena->current->next = block;
        } else {
            arena->first = block;
        }
        arena->current = block;
        arena->total_allocated += block_size;
    }

    void *ptr = arena->current->data + arena->current->used;
    arena->current->used += size;
    return ptr;
}

void cimport_arena_free(CImportArena *arena)
{
    CImportArenaBlock *block = arena->first;
    while (block) {
        CImportArenaBlock *next = block->next;
        free(block);
        block = next;
    }
    arena->first = NULL;
    arena->current = NULL;
    arena->total_allocated = 0;
}
