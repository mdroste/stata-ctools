#include "ctools_group.h"
#include "ctools_config.h"
#include <stdlib.h>
#include <string.h>

/* Check if two adjacent observations are in the same group */
static int same_group_check(const stata_data *data, const perm_idx_t *order, const int *keys,
                            size_t nvars, size_t i)
{
    size_t previous = order ? order[i - 1] : i - 1;
    size_t current = order ? order[i] : i;
    for (size_t b = 0; b < nvars; b++)
    {
        const stata_variable *var = &data->vars[keys ? (size_t)keys[b] : b];
        if (var->type == STATA_TYPE_STRING)
        {
            const char *prev = var->data.str[previous];
            const char *curr = var->data.str[current];
            if (strcmp(prev ? prev : "", curr ? curr : "") != 0)
                return 0;
        }
        else if (var->data.dbl[previous] != var->data.dbl[current])
        {
            return 0;
        }
    }
    return 1;
}

/* Detect groups from sorted stata_data (data must already be sorted/permuted)
   Parallelized for large datasets using boundary detection */
void ctools_group_ranges(const stata_data *data, const perm_idx_t *order, const int *keys,
                         size_t nvars, ctools_group_range *groups, size_t *ngroups)
{
    size_t nobs = data->nobs;
    if (nobs == 0)
    {
        *ngroups = 0;
        return;
    }

#ifdef _OPENMP
    int num_threads = ctools_get_max_threads();

    /* For small datasets or few threads, use sequential version */
    if (nobs < 10000 || num_threads <= 1)
    {
#endif
        /* Sequential version */
        size_t g = 0;
        groups[0].start = 0;
        groups[0].count = 1;

        for (size_t i = 1; i < nobs; i++)
        {
            if (same_group_check(data, order, keys, nvars, i))
            {
                groups[g].count++;
            }
            else
            {
                g++;
                groups[g].start = i;
                groups[g].count = 1;
            }
        }
        *ngroups = g + 1;
        return;
#ifdef _OPENMP
    }

    /* Parallel version: first find all group boundaries */
    uint8_t *boundaries = (uint8_t *)calloc(nobs, sizeof(uint8_t));
    if (!boundaries)
    {
        /* Fallback to sequential */
        size_t g = 0;
        groups[0].start = 0;
        groups[0].count = 1;
        for (size_t i = 1; i < nobs; i++)
        {
            if (same_group_check(data, order, keys, nvars, i))
            {
                groups[g].count++;
            }
            else
            {
                g++;
                groups[g].start = i;
                groups[g].count = 1;
            }
        }
        *ngroups = g + 1;
        return;
    }

    boundaries[0] = 1; /* First obs is always a boundary */

/* Phase 1: Parallel boundary detection */
#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (size_t i = 1; i < nobs; i++)
    {
        if (!same_group_check(data, order, keys, nvars, i))
        {
            boundaries[i] = 1;
        }
    }

    /* Phase 2: Sequential group construction */
    size_t g = 0;
    for (size_t i = 0; i < nobs; i++)
    {
        if (boundaries[i])
        {
            if (g > 0)
            {
                groups[g - 1].count = i - groups[g - 1].start;
            }
            groups[g].start = i;
            g++;
        }
    }
    if (g > 0)
    {
        groups[g - 1].count = nobs - groups[g - 1].start;
    }

    free(boundaries);
    *ngroups = g;
#endif
}
