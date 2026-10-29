#ifndef CPSMATCH_INDEX_H
#define CPSMATCH_INDEX_H

#include "ctools_data.h"

/* Path-compressed predecessor/successor indexes over sorted control positions.
 * The caller owns n+1 initialized parent entries in each direction. Successor
 * coordinates are positions 0..n (n is absent); predecessor coordinates are
 * positions plus one, with zero absent. Deletion links to the next live root. */
static inline perm_idx_t cpsmatch_available_find(perm_idx_t *parent, perm_idx_t position)
{
    while (parent[position] != position)
    {
        parent[position] = parent[parent[position]];
        position = parent[position];
    }
    return position;
}

#endif
