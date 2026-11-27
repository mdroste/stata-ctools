/*
    cfuzzyjoin_lev.h
    Bit-parallel Levenshtein distance for patterns of at most 64 bytes
*/

#ifndef CFUZZYJOIN_LEV_H
#define CFUZZYJOIN_LEV_H

#include <stdint.h>

typedef struct
{
    uint64_t peq[256];
    uint32_t m;
} cfj_myers;

/* Pattern a, m <= 64 bytes */
void cfj_myers_init(cfj_myers *p, const unsigned char *a, uint32_t m);
/* The distance between the pattern and b when it is at most k; otherwise
 * some value above k */
uint32_t cfj_myers_dist(const cfj_myers *p, const unsigned char *b, uint32_t n, uint32_t k);

#endif /* CFUZZYJOIN_LEV_H */
