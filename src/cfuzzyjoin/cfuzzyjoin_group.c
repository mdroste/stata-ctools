/*
    cfuzzyjoin_group.c
    group: strgroup's grouping of similar strings

    Two strings match when their Levenshtein distance divided by the length
    of the shorter string (normalize(shorter), the default), of the longer
    one (longer) or by one (none) is at most the threshold; groups are the
    connected components of that relation, numbered in order of their first
    observation. With by-groups the groups never cross a by-group and are
    numbered by-group after by-group; with first, strings are also split by
    their first character (other, then a to z), as strgroup does it.

    The plugin call (one per command):

        group nby=<k> normalize=<0|1|2> threshold(<text>) [first] [clean]
              [verbose]

    varlist: the string variable, the k by-variables, the new group variable
    and the sample marker (last). The data are sorted by the by-variables.

    Within a part (a by-group, or one first character of it) equal strings
    are one node. Nodes are scanned in order of length: a node is compared
    with the longer ones until the length difference alone exceeds the
    allowed distance, pairs already in one group are skipped, and a bound
    from the characters only one string has saves most distance
    computations. Matching pairs are joined in a lock-free union-find, so
    the parts' nodes are scanned in parallel; the groups do not depend on
    the order in which pairs are found.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdatomic.h>

#include "cfuzzyjoin_internal.h"
#include "cfuzzyjoin_lev.h"
#include "ctools_config.h"
#include "ctools_cancel.h"
#include "ctools_runtime.h"
#include "ctools_spi.h"
#include "ctools_parse.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Parts at least this large are scanned by all threads together; smaller
 * ones are handed out whole, one thread each. When there are fewer parts than
 * threads (no by-groups, say), parts down to CFJ_GROUP_SHARED strings are
 * also scanned by all threads. */
#define CFJ_GROUP_BIG 2048
#define CFJ_GROUP_SHARED 128

typedef struct
{
    const unsigned char *s;
    uint32_t len;
    uint64_t mask; /* byte classes present (byte & 63) */
} node;

typedef struct
{
    double t;
    int normalize; /* 0 shorter, 1 longer, 2 none */
} rule;

/* The largest distance d <= cap with d / norm <= t, as strgroup compares
 * (an unsigned distance divided by a double length) */
static uint32_t kmax_for(const rule *r, uint32_t norm, uint32_t cap)
{
    double n = (double)norm;
    double x = r->t * n;
    uint32_t d = x >= (double)cap ? cap : (x > 0 ? (uint32_t)x : 0);
    while (d < cap && (double)(d + 1) / n <= r->t)
        d++;
    while (d > 0 && (double)d / n > r->t)
        d--;
    return d;
}

static uint32_t uf_find(_Atomic uint32_t *par, uint32_t x)
{
    for (;;)
    {
        uint32_t p = atomic_load_explicit(&par[x], memory_order_relaxed);
        if (p == x)
            return x;
        uint32_t g = atomic_load_explicit(&par[p], memory_order_relaxed);
        if (g == p)
            return p;
        atomic_compare_exchange_weak_explicit(&par[x], &p, g, memory_order_relaxed,
                                              memory_order_relaxed);
        x = g;
    }
}

static void uf_union(_Atomic uint32_t *par, uint32_t a, uint32_t b)
{
    for (;;)
    {
        a = uf_find(par, a);
        b = uf_find(par, b);
        if (a == b)
            return;
        /* the root with the larger index goes under the smaller: no cycles */
        if (a < b)
        {
            uint32_t t = a;
            a = b;
            b = t;
        }
        uint32_t expect = a;
        if (atomic_compare_exchange_strong_explicit(&par[a], &expect, b, memory_order_relaxed,
                                                    memory_order_relaxed))
            return;
    }
}

static int popcount64(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0;
    while (x)
        x &= x - 1, c++;
    return c;
#endif
}

/* Compare node order[p] with the longer nodes after it */
static void scan_one(const node *nd, const uint32_t *order, uint32_t nn, uint32_t p, const rule *r,
                     uint32_t maxlen, _Atomic uint32_t *par, uint32_t *work, cfj_myers *my,
                     uint64_t *pairs)
{
    const uint32_t i = order[p];
    const node *a = &nd[i];
    uint32_t kfix = 0;
    if (r->normalize == 0)
        kfix = kmax_for(r, a->len, maxlen);
    else if (r->normalize == 2)
        kfix = kmax_for(r, 1, maxlen);
    int prepared = 0;
    for (uint32_t q = p + 1; q < nn; q++)
    {
        const uint32_t j = order[q];
        const node *b = &nd[j];
        const uint32_t diff = b->len - a->len;
        const uint32_t k = r->normalize == 1 ? kmax_for(r, b->len, b->len) : kfix;
        if (diff > k)
            break;
        if (uf_find(par, i) == uf_find(par, j))
            continue;
        int ao = popcount64(a->mask & ~b->mask), bo = popcount64(b->mask & ~a->mask);
        if ((uint32_t)(ao > bo ? ao : bo) > k)
            continue;
        uint32_t dist;
        (*pairs)++;
        if (a->len <= 64)
        {
            if (!prepared)
            {
                cfj_myers_init(my, a->s, a->len);
                prepared = 1;
            }
            dist = cfj_myers_dist(my, b->s, b->len, k);
        }
        else
            dist = cfj_lev_bounded(a->s, a->len, b->s, b->len, k, work);
        if (dist <= k)
            uf_union(par, i, j);
    }
}

typedef struct
{
    size_t first, count; /* rows of the part, in rows[] */
    uint32_t nnodes;
    uint32_t ngroups;
} part;

typedef struct
{
    uint64_t *hash;
    uint32_t *slot;
    size_t mask;
} strtab;

static int strtab_reserve(strtab *t, size_t n)
{
    size_t c = 64;
    while (c < n * 2)
        c *= 2;
    if (c <= t->mask + 1 && t->slot)
    {
        memset(t->slot, 0, (t->mask + 1) * sizeof(uint32_t));
        return 0;
    }
    free(t->hash);
    free(t->slot);
    t->hash = (uint64_t *)malloc(c * sizeof(uint64_t));
    t->slot = (uint32_t *)calloc(c, sizeof(uint32_t));
    t->mask = c - 1;
    return t->hash && t->slot ? 0 : -1;
}

/* Number the nodes of one part and its groups. node_of[r] is the part's node
 * of row r (relative to the part); group_of_row[r] its group within the part
 * (1-based). Returns 0 or an error. big: scan with every thread. */
static int run_part(part *pt, const uint32_t *rows, const unsigned char *const *str,
                    const uint32_t *len, const rule *r, strtab *tab, uint32_t *node_of,
                    uint32_t *group_of_row, int threads, int big, uint64_t *pairs)
{
    const size_t n = pt->count;
    if (strtab_reserve(tab, n))
        return STATA_ERR_MEMORY;
    node *nd = (node *)malloc((n ? n : 1) * sizeof(node));
    if (!nd)
        return STATA_ERR_MEMORY;
    uint32_t nn = 0, maxlen = 0;
    for (size_t k = 0; k < n; k++)
    {
        const uint32_t row = rows[pt->first + k];
        const unsigned char *s = str[row];
        const uint32_t l = len[row];
        uint64_t h = cfj_hash_bytes((const char *)s, l);
        size_t at = (size_t)h & tab->mask;
        uint32_t found = 0;
        while (tab->slot[at])
        {
            uint32_t id = tab->slot[at] - 1;
            if (tab->hash[at] == h && nd[id].len == l && !memcmp(nd[id].s, s, l))
            {
                found = id + 1;
                break;
            }
            at = (at + 1) & tab->mask;
        }
        if (!found)
        {
            nd[nn].s = s;
            nd[nn].len = l;
            uint64_t m = 0;
            for (uint32_t c = 0; c < l; c++)
                m |= (uint64_t)1 << (s[c] & 63);
            nd[nn].mask = m;
            if (l > maxlen)
                maxlen = l;
            tab->slot[at] = nn + 1;
            tab->hash[at] = h;
            found = ++nn;
        }
        node_of[k] = found - 1;
    }
    pt->nnodes = nn;

    _Atomic uint32_t *par = (_Atomic uint32_t *)malloc((nn ? nn : 1) * sizeof(*par));
    uint32_t *order = (uint32_t *)malloc((nn ? nn : 1) * sizeof(uint32_t));
    if (!par || !order)
    {
        free(nd);
        free(par);
        free(order);
        return STATA_ERR_MEMORY;
    }
    for (uint32_t i = 0; i < nn; i++)
        atomic_init(&par[i], i);
    /* counting sort of the nodes by length (stable) */
    {
        uint32_t *cnt = (uint32_t *)calloc((size_t)maxlen + 2, sizeof(uint32_t));
        if (!cnt)
        {
            free(nd);
            free(par);
            free(order);
            return STATA_ERR_MEMORY;
        }
        for (uint32_t i = 0; i < nn; i++)
            cnt[nd[i].len + 1]++;
        for (uint32_t l = 0; l <= maxlen; l++)
            cnt[l + 1] += cnt[l];
        for (uint32_t i = 0; i < nn; i++)
            order[cnt[nd[i].len]++] = i;
        free(cnt);
    }

    int rc = STATA_OK;
    if (nn > 1)
    {
        if (big && threads > 1)
        {
            atomic_int failed = 0;
            uint64_t total = 0;
#pragma omp parallel num_threads(threads) reduction(+ : total)
            {
                uint32_t *work = (uint32_t *)malloc(2 * ((size_t)maxlen + 2) * sizeof(uint32_t));
                cfj_myers *my = (cfj_myers *)malloc(sizeof(cfj_myers));
                if (!work || !my)
                    atomic_store(&failed, STATA_ERR_MEMORY);
#pragma omp for schedule(dynamic, 16)
                for (uint32_t p = 0; p < nn; p++)
                {
                    if (atomic_load_explicit(&failed, memory_order_relaxed))
                        continue;
                    if ((p & 255) == 0 && ctools_cancel_requested())
                    {
                        atomic_store(&failed, STATA_ERR_CANCELLED);
                        continue;
                    }
                    scan_one(nd, order, nn, p, r, maxlen, par, work, my, &total);
                }
                free(work);
                free(my);
            }
            *pairs += total;
            rc = atomic_load(&failed);
        }
        else
        {
            uint32_t *work = (uint32_t *)malloc(2 * ((size_t)maxlen + 2) * sizeof(uint32_t));
            cfj_myers *my = (cfj_myers *)malloc(sizeof(cfj_myers));
            if (!work || !my)
                rc = STATA_ERR_MEMORY;
            for (uint32_t p = 0; p < nn && !rc; p++)
            {
                if ((p & 1023) == 1023 && ctools_cancel_requested())
                    rc = STATA_ERR_CANCELLED;
                else
                    scan_one(nd, order, nn, p, r, maxlen, par, work, my, pairs);
            }
            free(work);
            free(my);
        }
    }

    /* groups in order of first observation: nodes are in that order */
    if (!rc)
    {
        uint32_t *gid = order; /* reused: group of each root */
        for (uint32_t i = 0; i < nn; i++)
            gid[i] = 0;
        uint32_t ng = 0;
        for (size_t k = 0; k < n; k++)
        {
            uint32_t root = uf_find(par, node_of[k]);
            if (!gid[root])
                gid[root] = ++ng;
            group_of_row[k] = gid[root];
        }
        pt->ngroups = ng;
    }
    free(nd);
    free((void *)par);
    free(order);
    return rc;
}

static int same_by(const ctools_filtered_data *fd, int nby, size_t a, size_t b)
{
    for (int v = 1; v <= nby; v++)
    {
        const stata_variable *var = &fd->data.vars[v];
        if (var->type == STATA_TYPE_STRING)
        {
            const char *x = var->data.str[a] ? var->data.str[a] : "";
            const char *y = var->data.str[b] ? var->data.str[b] : "";
            if (strcmp(x, y))
                return 0;
        }
        else
        {
            double x = var->data.dbl[a], y = var->data.dbl[b];
            if (!(x == y || (x != x && y != y)))
                return 0;
        }
    }
    return 1;
}

static int bucket_of(unsigned char c)
{
    if (c >= 'a' && c <= 'z')
        return 1 + (c - 'a');
    if (c >= 'A' && c <= 'Z')
        return 1 + (c - 'A');
    return 0;
}

ST_retcode cfj_group(const char *args)
{
    int nby = 0, normalize = 0;
    char ttext[64];
    if (cfj_opt_int(args, "nby", &nby) || cfj_opt_int(args, "normalize", &normalize) ||
        !cfj_opt_text(args, "threshold", ttext, sizeof(ttext)) || nby < 0 || normalize < 0 ||
        normalize > 2 || SF_nvars() != nby + 3)
        return 198;
    rule r;
    char *end;
    r.t = strtod(ttext, &end);
    r.normalize = normalize;
    if (end == ttext || !(r.t > 0))
        return 198;
    /* strgroup's plugin matches as if the threshold had seven decimal places
     * (found by sweeping thresholds around exact ratios such as 1/7, 1/8 and
     * 1/20): use it rounded the same way */
    {
        char rounded[64];
        snprintf(rounded, sizeof(rounded), "%.7f", r.t);
        /* below 0.00000005 that is 0: only equal strings match */
        r.t = strtod(rounded, NULL);
    }
    const int first = cfj_flag(args, "first");
    const int clean = cfj_flag(args, "clean");
    const int verbose = cfj_flag(args, "verbose");
    const int threads = ctools_get_max_threads();
    const double t0 = ctools_timer_seconds();

    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    stata_retcode lrc = cfj_load(&fd, 1, 1 + nby, CTOOLS_LOAD_MARKER);
    if (lrc)
    {
        ctools_filtered_data_free(&fd);
        if (lrc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return ctools_stata_rc(lrc);
    }
    if (fd.data.vars[0].type != STATA_TYPE_STRING)
    {
        ctools_filtered_data_free(&fd);
        return 109;
    }
    const size_t n = fd.data.nobs;
    const double t1 = ctools_timer_seconds();

    stata_retcode rc = STATA_OK;
    const unsigned char **str = (const unsigned char **)malloc((n ? n : 1) * sizeof(char *));
    uint32_t *len = (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    uint32_t *key = (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    double *out = (double *)malloc((n ? n : 1) * sizeof(double));
    if (!str || !len || !key || !out)
        rc = STATA_ERR_MEMORY;

    /* Cleaned strings, parts (by-group * 27 + first character) */
    size_t nparts_max = 0, used = 0;
    uint32_t byid = 0;
    if (!rc)
    {
        size_t prev = (size_t)-1;
        for (size_t i = 0; i < n; i++)
        {
            const char *s = fd.data.vars[0].data.str[i];
            if (!s)
                s = "";
            size_t a = 0, b = strlen(s);
            if (clean)
            {
                while (a < b && s[a] == ' ')
                    a++;
                while (b > a && s[b - 1] == ' ')
                    b--;
            }
            str[i] = (const unsigned char *)s + a;
            len[i] = (uint32_t)(b - a);
            out[i] = SV_missval;
            if (b == a)
            {
                key[i] = UINT32_MAX;
                continue;
            }
            if (nby && prev != (size_t)-1 && !same_by(&fd, nby, prev, i))
                byid++;
            prev = i;
            key[i] = byid * 27u + (first ? (uint32_t)bucket_of(str[i][0]) : 0u);
            used++;
        }
        nparts_max = ((size_t)byid + 1) * 27;
    }

    /* Rows of each part, in observation order (counting sort by part) */
    uint32_t *rows = NULL;
    size_t *start = NULL;
    part *parts = NULL;
    size_t nparts = 0;
    if (!rc)
    {
        start = (size_t *)calloc(nparts_max + 1, sizeof(size_t));
        rows = (uint32_t *)malloc((used ? used : 1) * sizeof(uint32_t));
        if (!start || !rows)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc)
    {
        for (size_t i = 0; i < n; i++)
            if (key[i] != UINT32_MAX)
                start[key[i] + 1]++;
        for (size_t p = 0; p < nparts_max; p++)
        {
            if (start[p + 1])
                nparts++;
            start[p + 1] += start[p];
        }
        parts = (part *)calloc(nparts ? nparts : 1, sizeof(part));
        if (!parts)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc)
    {
        size_t *fill = (size_t *)malloc((nparts_max + 1) * sizeof(size_t));
        if (!fill)
            rc = STATA_ERR_MEMORY;
        else
        {
            memcpy(fill, start, (nparts_max + 1) * sizeof(size_t));
            for (size_t i = 0; i < n; i++)
                if (key[i] != UINT32_MAX)
                    rows[fill[key[i]]++] = (uint32_t)i;
            free(fill);
            size_t k = 0;
            for (size_t p = 0; p < nparts_max; p++)
            {
                if (start[p + 1] > start[p])
                {
                    parts[k].first = start[p];
                    parts[k].count = start[p + 1] - start[p];
                    k++;
                }
            }
        }
    }

    /* Group each part */
    uint32_t *node_of = NULL, *gpart = NULL;
    uint64_t pairs = 0;
    if (!rc)
    {
        node_of = (uint32_t *)malloc((used ? used : 1) * sizeof(uint32_t));
        gpart = (uint32_t *)malloc((used ? used : 1) * sizeof(uint32_t));
        if (!node_of || !gpart)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc)
    {
        /* big parts one after another, all threads on each */
        const size_t big = nparts < (size_t)threads ? CFJ_GROUP_SHARED : CFJ_GROUP_BIG;
        strtab tab = {0};
        for (size_t p = 0; p < nparts && !rc; p++)
        {
            if (parts[p].count < big)
                continue;
            rc = run_part(&parts[p], rows, str, len, &r, &tab, node_of + parts[p].first,
                          gpart + parts[p].first, threads, 1, &pairs);
        }
        free(tab.hash);
        free(tab.slot);
        /* the small parts handed out whole */
        if (!rc)
        {
            atomic_int failed = 0;
            uint64_t total = 0;
            int team = threads > 1 && nparts > 1 ? threads : 1;
#pragma omp parallel num_threads(team) if (team > 1) reduction(+ : total)
            {
                strtab local = {0};
#pragma omp for schedule(dynamic, 8)
                for (size_t p = 0; p < nparts; p++)
                {
                    if (parts[p].count >= big ||
                        atomic_load_explicit(&failed, memory_order_relaxed))
                        continue;
                    if ((p & 63) == 0 && ctools_cancel_requested())
                    {
                        atomic_store(&failed, STATA_ERR_CANCELLED);
                        continue;
                    }
                    int prc = run_part(&parts[p], rows, str, len, &r, &local,
                                       node_of + parts[p].first, gpart + parts[p].first, 1, 0,
                                       &total);
                    if (prc)
                        atomic_store(&failed, prc);
                }
                free(local.hash);
                free(local.slot);
            }
            pairs += total;
            rc = atomic_load(&failed);
        }
    }

    /* Number the groups part after part */
    double ngroups = 0;
    size_t nnodes = 0;
    if (!rc)
    {
        uint64_t base = 0;
        for (size_t p = 0; p < nparts; p++)
        {
            for (size_t k = 0; k < parts[p].count; k++)
                out[rows[parts[p].first + k]] = (double)(base + gpart[parts[p].first + k]);
            base += parts[p].ngroups;
            nnodes += parts[p].nnodes;
        }
        ngroups = (double)base;
    }
    const double t2 = ctools_timer_seconds();
    if (!rc)
        rc = ctools_store_filtered(out, n, nby + 2, fd.obs_map);
    const double t3 = ctools_timer_seconds();

    free(str);
    free(len);
    free(key);
    free(out);
    free(rows);
    free(start);
    free(parts);
    free(node_of);
    free(gpart);
    ctools_filtered_data_free(&fd);
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return ctools_stata_rc(rc);
    }
    ST_retcode src = 0;
    src = src ? src : ctools_scal_save("_cfuzzyjoin_N", (double)used);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_ngroups", ngroups);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_nunique", (double)nnodes);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_nparts", (double)nparts);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_pairs", (double)pairs);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_load", t1 - t0);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_compute", t2 - t1);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_store", t3 - t2);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_threads", (double)threads);
    if (src)
        return src;
    ctools_verbose(CFJ_MODULE, verbose,
                   "group: %zu strings, %zu distinct within %zu parts, %.0f distance computations, "
                   "%.0f groups; load %.4fs, group %.4fs, store %.4fs",
                   used, nnodes, nparts, (double)pairs, ngroups, t1 - t0, t2 - t1, t3 - t2);
    return 0;
}
