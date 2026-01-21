/*
 * cmerge_memory.c
 * Memory allocation utilities for cmerge
 *
 * Thin wrappers around ctools_arena functions, configured with
 * STATIC_FALLBACK mode for cmerge's never-NULL semantics.
 */

#include <stdlib.h>
#include <string.h>
#include "cmerge_memory.h"

/* ============================================================================
 * String Arena Allocator - Wrappers around ctools_string_arena
 * ============================================================================ */

cmerge_string_arena *cmerge_arena_create(size_t capacity)
{
    /* Create arena with STATIC_FALLBACK mode - never returns NULL from strdup */
    return ctools_string_arena_create(capacity, CTOOLS_STRING_ARENA_STATIC_FALLBACK);
}

char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s)
{
    return ctools_string_arena_strdup(arena, s);
}

int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr)
{
    return ctools_string_arena_owns(arena, ptr);
}

void cmerge_arena_free(cmerge_string_arena *arena)
{
    ctools_string_arena_free(arena);
}
