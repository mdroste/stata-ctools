/*
 * cmerge_memory.c
 * Memory allocation utilities for cmerge
 */

#include <stdlib.h>
#include <string.h>
#include "cmerge_memory.h"
#include "../ctools_config.h"

#ifdef _WIN32
#include <malloc.h>
#endif

/* ============================================================================
 * Cache-Line Aligned Memory Allocation
 * ============================================================================ */

void *cmerge_aligned_alloc(size_t size)
{
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, CACHE_LINE_SIZE);
#else
    if (posix_memalign(&ptr, CACHE_LINE_SIZE, size) != 0) {
        return NULL;
    }
#endif
    return ptr;
}

void cmerge_aligned_free(void *ptr)
{
    if (ptr == NULL) return;
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/* ============================================================================
 * String Arena Allocator
 * ============================================================================ */

/* Static empty string for fallback when all allocations fail.
 * Use volatile to prevent compiler from optimizing away pointer comparisons.
 * Store both the array and a pointer to it so we can reliably detect it later. */
static char cmerge_empty_string_storage[1] = "";
static char * volatile cmerge_empty_string = cmerge_empty_string_storage;

/* Check if a pointer is the static fallback string (must not be freed) */
static int cmerge_is_static_string(const char *ptr)
{
    /* Compare against both the storage address and the pointer variable
     * to handle any compiler optimization scenarios */
    return (ptr == cmerge_empty_string_storage ||
            ptr == cmerge_empty_string ||
            ptr == (char *)&cmerge_empty_string_storage[0]);
}

cmerge_string_arena *cmerge_arena_create(size_t capacity)
{
    cmerge_string_arena *arena = (cmerge_string_arena *)malloc(sizeof(cmerge_string_arena));
    if (arena == NULL) return NULL;

    arena->base = (char *)malloc(capacity);
    if (arena->base == NULL) {
        free(arena);
        return NULL;
    }

    arena->capacity = capacity;
    arena->used = 0;
    arena->alloc_failed = 0;
    return arena;
}

char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s)
{
    size_t len = strlen(s) + 1;

    if (arena != NULL && arena->used + len <= arena->capacity) {
        char *ptr = arena->base + arena->used;
        memcpy(ptr, s, len);
        arena->used += len;
        return ptr;
    }

    /* Fallback to regular strdup */
    char *result = strdup(s);
    if (result == NULL) {
        if (arena != NULL) {
            arena->alloc_failed = 1;
            /* Return pointer to empty string in arena if possible */
            if (arena->used < arena->capacity) {
                char *ptr = arena->base + arena->used;
                ptr[0] = '\0';
                arena->used += 1;
                return ptr;
            }
        }
        /* Last resort: return static empty string (never NULL)
         * Using the storage directly ensures cmerge_is_static_string() will detect it */
        return cmerge_empty_string_storage;
    }
    return result;
}

int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr)
{
    if (ptr == NULL) return 0;
    /* Static empty string should not be freed - use robust check */
    if (cmerge_is_static_string(ptr)) return 1;
    if (arena == NULL) return 0;
    return (ptr >= arena->base && ptr < arena->base + arena->capacity);
}

void cmerge_arena_free(cmerge_string_arena *arena)
{
    if (arena != NULL) {
        free(arena->base);
        free(arena);
    }
}
