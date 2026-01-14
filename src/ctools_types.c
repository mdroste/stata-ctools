/*
    ctools_types.c
    Implementation of common utility functions for Stata-C data structures
*/

#include <stdlib.h>
#include <string.h>
#include "ctools_types.h"
#include "ctools_config.h"

/*
    NOTE: Use ctools_aligned_free() from ctools_config.h for freeing aligned memory.
    This ensures consistent behavior across all platforms including Windows.
*/
#define aligned_free_internal ctools_aligned_free

void stata_data_init(stata_data *data)
{
    if (data == NULL) return;

    data->nobs = 0;
    data->nvars = 0;
    data->vars = NULL;
    data->sort_order = NULL;
}

/*
    Internal: Free a string arena (matches arena_create in ctools_data_io.c).
    This is duplicated here to avoid circular dependencies.
*/
typedef struct {
    char *base;
    size_t capacity;
    size_t used;
} string_arena_internal;

static void arena_free_internal(void *arena_ptr)
{
    if (arena_ptr != NULL) {
        string_arena_internal *arena = (string_arena_internal *)arena_ptr;
        free(arena->base);
        free(arena);
    }
}

void stata_data_free(stata_data *data)
{
    size_t i;

    if (data == NULL) return;

    /* Free each variable's data (allocated with aligned_alloc_cacheline) */
    if (data->vars != NULL) {
        for (i = 0; i < data->nvars; i++) {
            if (data->vars[i].type == STATA_TYPE_DOUBLE) {
                /* Numeric data uses aligned allocation */
                aligned_free_internal(data->vars[i].data.dbl);
            } else if (data->vars[i].type == STATA_TYPE_STRING) {
                if (data->vars[i].data.str != NULL) {
                    /* Check if strings were allocated via arena (fast path)
                       or individually via strdup (slow path) */
                    if (data->vars[i]._arena != NULL) {
                        /* Arena allocation: O(1) bulk free - just free the arena.
                           All strings are inside the arena's contiguous memory block. */
                        arena_free_internal(data->vars[i]._arena);
                        data->vars[i]._arena = NULL;
                    } else {
                        /* Legacy path: individual strdup allocations - O(n) frees.
                           This path is only used if arena allocation failed. */
                        for (size_t j = 0; j < data->vars[i].nobs; j++) {
                            free(data->vars[i].data.str[j]);
                        }
                    }
                    /* Free the pointer array (aligned allocation) */
                    aligned_free_internal(data->vars[i].data.str);
                }
            }
        }
        /* Free vars array (aligned) */
        aligned_free_internal(data->vars);
    }

    /* Free sort order array (aligned) */
    aligned_free_internal(data->sort_order);

    /* Reset the structure */
    stata_data_init(data);
}
