/*
    cfuzzyjoin_lev.c
    Levenshtein distance with unit costs over bytes, as strgroup computes it

    Patterns of at most 64 bytes use the bit-parallel algorithm of Myers
    (1999) in Hyyrö's form for the distance between whole strings: one word
    operation sequence per byte of the other string. Longer strings use the
    dynamic program restricted to the band |i - j| <= k (Ukkonen 1985). Both
    stop as soon as the distance is known to exceed the bound k.
*/

#include <string.h>

#include "cfuzzyjoin_internal.h"
#include "cfuzzyjoin_lev.h"

void cfj_myers_init(cfj_myers *p, const unsigned char *a, uint32_t m)
{
    memset(p->peq, 0, sizeof(p->peq));
    p->m = m;
    for (uint32_t i = 0; i < m; i++)
        p->peq[a[i]] |= (uint64_t)1 << i;
}

uint32_t cfj_myers_dist(const cfj_myers *p, const unsigned char *b, uint32_t n, uint32_t k)
{
    const uint32_t m = p->m;
    if (m == 0)
        return n;
    uint64_t pv = ~(uint64_t)0, mv = 0;
    const uint64_t last = (uint64_t)1 << (m - 1);
    uint32_t score = m;
    for (uint32_t j = 0; j < n; j++)
    {
        const uint64_t eq = p->peq[b[j]];
        const uint64_t xv = eq | mv;
        const uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
        uint64_t ph = mv | ~(xh | pv);
        uint64_t mh = pv & xh;
        if (ph & last)
            score++;
        else if (mh & last)
            score--;
        /* the top boundary row is D[0][j] = j: a +1 enters at bit 0 */
        ph = (ph << 1) | 1;
        mh <<= 1;
        pv = mh | ~(xv | ph);
        mv = ph & xv;
        /* each remaining byte lowers the final distance by at most one */
        if (score > k && score - k > n - 1 - j)
            return k + 1;
    }
    return score;
}

/* Banded dynamic program; la <= lb, lb - la <= k. work: 2 * (lb + 2). */
static uint32_t banded(const unsigned char *a, uint32_t la, const unsigned char *b, uint32_t lb,
                       uint32_t k, uint32_t *work)
{
    const uint32_t inf = k + 1;
    uint32_t *prev = work, *cur = work + (lb + 2);
    uint32_t top = lb < k ? lb : k;
    for (uint32_t j = 0; j <= top; j++)
        prev[j] = j;
    if (top + 1 <= lb)
        prev[top + 1] = inf;
    for (uint32_t i = 1; i <= la; i++)
    {
        uint32_t lo = i > k ? i - k : 1;
        uint32_t hi = i + k < lb ? i + k : lb;
        cur[lo - 1] = (lo == 1 && i <= k) ? i : inf;
        uint32_t rowmin = inf;
        const unsigned char ca = a[i - 1];
        for (uint32_t j = lo; j <= hi; j++)
        {
            uint32_t v = prev[j - 1] + (ca != b[j - 1]);
            uint32_t up = prev[j] + 1;
            uint32_t left = cur[j - 1] + 1;
            if (up < v)
                v = up;
            if (left < v)
                v = left;
            if (v > inf)
                v = inf;
            cur[j] = v;
            if (v < rowmin)
                rowmin = v;
        }
        if (hi < lb)
            cur[hi + 1] = inf;
        if (rowmin > k)
            return inf;
        uint32_t *t = prev;
        prev = cur;
        cur = t;
    }
    return prev[lb];
}

uint32_t cfj_lev_bounded(const unsigned char *a, uint32_t la, const unsigned char *b, uint32_t lb,
                         uint32_t k, uint32_t *work)
{
    while (la && lb && *a == *b)
        a++, b++, la--, lb--;
    while (la && lb && a[la - 1] == b[lb - 1])
        la--, lb--;
    if (la > lb)
    {
        const unsigned char *t = a;
        a = b;
        b = t;
        uint32_t tl = la;
        la = lb;
        lb = tl;
    }
    if (la == 0)
        return lb;
    if (lb - la > k)
        return k + 1;
    if (la <= 64)
    {
        cfj_myers p;
        cfj_myers_init(&p, a, la);
        return cfj_myers_dist(&p, b, lb, k);
    }
    return banded(a, la, b, lb, k, work);
}

uint32_t cfj_lev_exact(const unsigned char *a, uint32_t la, const unsigned char *b, uint32_t lb,
                       uint32_t *work)
{
    return cfj_lev_bounded(a, la, b, lb, la > lb ? la : lb, work);
}
