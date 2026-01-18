/*
 * cimport_arena.h
 * Arena allocator for cimport CSV parser
 *
 * Provides efficient bulk allocation with O(1) deallocation.
 * Used for parsed rows and string data during CSV import.
 */

#ifndef CIMPORT_ARENA_H
#define CIMPORT_ARENA_H

#include <stdlib.h>
#include <stddef.h>
#include "../ctools_config.h"

/* Arena block structure - linked list of allocation blocks */
typedef struct CImportArenaBlock {
    struct CImportArenaBlock *next;
    size_t used;
    size_t capacity;
    char data[];  /* Flexible array member */
} CImportArenaBlock;

/* Arena allocator - manages multiple blocks */
typedef struct {
    CImportArenaBlock *first;
    CImportArenaBlock *current;
    size_t total_allocated;
} CImportArena;

/* Initialize an arena (must be called before use) */
void cimport_arena_init(CImportArena *arena);

/* Allocate memory from the arena (8-byte aligned)
 * Returns NULL on allocation failure */
void *cimport_arena_alloc(CImportArena *arena, size_t size);

/* Free all memory in the arena (O(1) for all allocations) */
void cimport_arena_free(CImportArena *arena);

#endif /* CIMPORT_ARENA_H */
