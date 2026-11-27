/*
    cfuzzyjoin_vocab.c
    Gram vocabulary and per-string gram profiles

    The strings of every set (master, then using) are cut into contiguous
    chunks, one per task. A task makes each string's grams, numbers them in
    a dictionary of its own and keeps each string's distinct local ids with
    their counts. The task dictionaries are then merged in chunk order, so a
    gram's global id is the place of its first appearance over all strings:
    the same for any thread count. Each profile is finally renumbered and
    sorted by global id.
*/

#include <stdlib.h>
#include <string.h>

#include "cfuzzyjoin_internal.h"
#include "ctools_config.h"
#include "ctools_cancel.h"

#ifdef _OPENMP
#include <omp.h>
#endif

uint64_t cfj_hash_bytes(const char *s, size_t len)
{
    uint64_t h = 0x9E3779B97F4A7C15ull ^ ((uint64_t)len * 0xff51afd7ed558ccdull);
    while (len >= 8)
    {
        uint64_t k;
        memcpy(&k, s, 8);
        h ^= k * 0xc4ceb9fe1a85ec53ull;
        h = (h << 31) | (h >> 33);
        h *= 0x9E3779B97F4A7C15ull;
        s += 8;
        len -= 8;
    }
    uint64_t k = 0;
    for (size_t i = 0; i < len; i++)
        k |= (uint64_t)(unsigned char)s[i] << (8 * i);
    h ^= k * 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdull;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33;
    return h;
}

/* ----------------------------------------------------------------------------
 * Dictionary: gram bytes -> id in insertion order
 * -------------------------------------------------------------------------- */

typedef struct
{
    uint64_t *hash;
    uint32_t *slot; /* id + 1; 0 empty */
    size_t mask, n;
    char *bytes;
    size_t used, cap;
    uint64_t *off;
    uint32_t *len;
    uint64_t *ghash;
    size_t gcap;
} dict;

static void dict_free(dict *d)
{
    free(d->hash);
    free(d->slot);
    free(d->bytes);
    free(d->off);
    free(d->len);
    free(d->ghash);
    memset(d, 0, sizeof(*d));
}

static int dict_init(dict *d, size_t cap)
{
    memset(d, 0, sizeof(*d));
    size_t c = 64;
    while (c < cap * 2)
        c *= 2;
    d->hash = (uint64_t *)malloc(c * sizeof(uint64_t));
    d->slot = (uint32_t *)calloc(c, sizeof(uint32_t));
    if (!d->hash || !d->slot)
        return -1;
    d->mask = c - 1;
    return 0;
}

static int dict_rehash(dict *d)
{
    size_t c = (d->mask + 1) * 2;
    uint64_t *h = (uint64_t *)malloc(c * sizeof(uint64_t));
    uint32_t *s = (uint32_t *)calloc(c, sizeof(uint32_t));
    if (!h || !s)
    {
        free(h);
        free(s);
        return -1;
    }
    for (size_t id = 0; id < d->n; id++)
    {
        size_t i = (size_t)d->ghash[id] & (c - 1);
        while (s[i])
            i = (i + 1) & (c - 1);
        s[i] = (uint32_t)id + 1;
        h[i] = d->ghash[id];
    }
    free(d->hash);
    free(d->slot);
    d->hash = h;
    d->slot = s;
    d->mask = c - 1;
    return 0;
}

/* The id of s, inserting it when new. -1 when out of memory. */
static long long dict_put(dict *d, const char *s, size_t len, uint64_t h)
{
    size_t i = (size_t)h & d->mask;
    while (d->slot[i])
    {
        if (d->hash[i] == h)
        {
            uint32_t id = d->slot[i] - 1;
            if (d->len[id] == len && !memcmp(d->bytes + d->off[id], s, len))
                return id;
        }
        i = (i + 1) & d->mask;
    }
    if (d->n >= 0xFFFFFFF0u)
        return -1;
    if (d->n + 1 > d->gcap)
    {
        size_t c = d->gcap ? d->gcap * 2 : 256;
        uint64_t *off = (uint64_t *)realloc(d->off, c * sizeof(uint64_t));
        if (!off)
            return -1;
        d->off = off;
        uint32_t *l = (uint32_t *)realloc(d->len, c * sizeof(uint32_t));
        if (!l)
            return -1;
        d->len = l;
        uint64_t *gh = (uint64_t *)realloc(d->ghash, c * sizeof(uint64_t));
        if (!gh)
            return -1;
        d->ghash = gh;
        d->gcap = c;
    }
    if (d->used + len + 1 > d->cap)
    {
        size_t c = d->cap ? d->cap : 4096;
        while (c < d->used + len + 1)
            c *= 2;
        char *b = (char *)realloc(d->bytes, c);
        if (!b)
            return -1;
        d->bytes = b;
        d->cap = c;
    }
    memcpy(d->bytes + d->used, s, len);
    uint32_t id = (uint32_t)d->n++;
    d->off[id] = d->used;
    d->len[id] = (uint32_t)len;
    d->ghash[id] = h;
    d->used += len;
    d->slot[i] = id + 1;
    d->hash[i] = h;
    if (d->n * 2 > d->mask + 1 && dict_rehash(d))
        return -1;
    return id;
}

/* ----------------------------------------------------------------------------
 * Vocabulary
 * -------------------------------------------------------------------------- */

void cfj_vocab_free(cfj_vocab *v)
{
    free(v->bytes);
    free(v->off);
    memset(v, 0, sizeof(*v));
}

void cfj_profile_free(cfj_profile *p)
{
    free(p->off);
    free(p->id);
    free(p->cnt);
    memset(p, 0, sizeof(*p));
}

void cfj_vocab_index_free(cfj_vocab_index *ix)
{
    free(ix->slot);
    memset(ix, 0, sizeof(*ix));
}

stata_retcode cfj_vocab_index_build(const cfj_vocab *v, cfj_vocab_index *ix)
{
    size_t c = 64;
    while (c < v->n * 2 + 2)
        c *= 2;
    ix->slot = (uint32_t *)calloc(c, sizeof(uint32_t));
    if (!ix->slot)
        return STATA_ERR_MEMORY;
    ix->mask = c - 1;
    for (size_t id = 0; id < v->n; id++)
    {
        const char *s = v->bytes + v->off[id];
        size_t len = (size_t)(v->off[id + 1] - v->off[id]);
        size_t i = (size_t)cfj_hash_bytes(s, len) & ix->mask;
        while (ix->slot[i])
            i = (i + 1) & ix->mask;
        ix->slot[i] = (uint32_t)id + 1;
    }
    return STATA_OK;
}

long long cfj_vocab_find(const cfj_vocab *v, const cfj_vocab_index *ix, const char *s, size_t len)
{
    size_t i = (size_t)cfj_hash_bytes(s, len) & ix->mask;
    while (ix->slot[i])
    {
        uint32_t id = ix->slot[i] - 1;
        size_t l = (size_t)(v->off[id + 1] - v->off[id]);
        if (l == len && !memcmp(v->bytes + v->off[id], s, len))
            return id;
        i = (i + 1) & ix->mask;
    }
    return -1;
}

/* ----------------------------------------------------------------------------
 * Profiles
 * -------------------------------------------------------------------------- */

typedef struct
{
    dict d;
    size_t first, last; /* strings [first, last) of the concatenation */
    uint32_t *nd;       /* distinct grams of each string */
    uint32_t *id;       /* local ids, then global ids */
    uint32_t *cnt;
    size_t used, cap;
    uint32_t *map; /* local -> global */
    int rc;
} task;

static int cmp_u32(const void *a, const void *b)
{
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return x < y ? -1 : x > y;
}

static void sort_u32(uint32_t *a, size_t n)
{
    if (n < 24)
    {
        for (size_t i = 1; i < n; i++)
        {
            uint32_t v = a[i];
            size_t j = i;
            while (j && a[j - 1] > v)
            {
                a[j] = a[j - 1];
                j--;
            }
            a[j] = v;
        }
        return;
    }
    qsort(a, n, sizeof(uint32_t), cmp_u32);
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

/* (id, count) pairs sorted by id */
static void sort_pairs(uint32_t *id, uint32_t *cnt, size_t n)
{
    if (n > 32)
    {
        uint64_t *w = (uint64_t *)malloc(n * sizeof(uint64_t));
        if (w)
        {
            for (size_t i = 0; i < n; i++)
                w[i] = ((uint64_t)id[i] << 32) | cnt[i];
            qsort(w, n, sizeof(uint64_t), cmp_u64);
            for (size_t i = 0; i < n; i++)
            {
                id[i] = (uint32_t)(w[i] >> 32);
                cnt[i] = (uint32_t)w[i];
            }
            free(w);
            return;
        }
    }
    for (size_t i = 1; i < n; i++)
    {
        uint32_t v = id[i], c = cnt[i];
        size_t j = i;
        while (j && id[j - 1] > v)
        {
            id[j] = id[j - 1];
            cnt[j] = cnt[j - 1];
            j--;
        }
        id[j] = v;
        cnt[j] = c;
    }
}

static const char *string_at(const cfj_strset *sets, int nsets, size_t k, size_t *len)
{
    for (int s = 0; s < nsets; s++)
    {
        if (k < sets[s].n)
        {
            const char *p = sets[s].s[k];
            if (!p)
                p = "";
            *len = strlen(p);
            return p;
        }
        k -= sets[s].n;
    }
    *len = 0;
    return "";
}

static void run_task(task *t, const cfj_method *m, const cfj_strset *sets, int nsets)
{
    cfj_grams g;
    cfj_grams_init(&g);
    uint32_t *tmp = NULL;
    size_t tcap = 0;
    size_t n = t->last - t->first;
    t->nd = (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    if (!t->nd || dict_init(&t->d, 1024))
    {
        t->rc = STATA_ERR_MEMORY;
        goto done;
    }
    for (size_t k = t->first; k < t->last; k++)
    {
        if (((k - t->first) & 1023) == 1023 && ctools_cancel_requested())
        {
            t->rc = STATA_ERR_CANCELLED;
            goto done;
        }
        size_t len;
        const char *s = string_at(sets, nsets, k, &len);
        if (cfj_extract(m, s, len, &g))
        {
            t->rc = STATA_ERR_MEMORY;
            goto done;
        }
        if (g.n > tcap)
        {
            size_t c = tcap ? tcap : 64;
            while (c < g.n)
                c *= 2;
            uint32_t *p = (uint32_t *)realloc(tmp, c * sizeof(uint32_t));
            if (!p)
            {
                t->rc = STATA_ERR_MEMORY;
                goto done;
            }
            tmp = p;
            tcap = c;
        }
        for (size_t j = 0; j < g.n; j++)
        {
            const char *gs = g.buf + g.off[j];
            long long id = dict_put(&t->d, gs, g.len[j], cfj_hash_bytes(gs, g.len[j]));
            if (id < 0)
            {
                t->rc = STATA_ERR_MEMORY;
                goto done;
            }
            tmp[j] = (uint32_t)id;
        }
        sort_u32(tmp, g.n);
        size_t nd = 0;
        for (size_t j = 0; j < g.n; j++)
            if (!j || tmp[j] != tmp[j - 1])
                nd++;
        if (t->used + nd > t->cap)
        {
            size_t c = t->cap ? t->cap : 1024;
            while (c < t->used + nd)
                c *= 2;
            uint32_t *a = (uint32_t *)realloc(t->id, c * sizeof(uint32_t));
            if (!a)
            {
                t->rc = STATA_ERR_MEMORY;
                goto done;
            }
            t->id = a;
            uint32_t *b = (uint32_t *)realloc(t->cnt, c * sizeof(uint32_t));
            if (!b)
            {
                t->rc = STATA_ERR_MEMORY;
                goto done;
            }
            t->cnt = b;
            t->cap = c;
        }
        for (size_t j = 0; j < g.n; j++)
        {
            if (!j || tmp[j] != tmp[j - 1])
            {
                t->id[t->used] = tmp[j];
                t->cnt[t->used] = 1;
                t->used++;
            }
            else if (!m->unit)
                t->cnt[t->used - 1]++;
        }
        t->nd[k - t->first] = (uint32_t)nd;
    }
done:
    free(tmp);
    cfj_grams_free(&g);
}

stata_retcode cfj_profiles_build(const cfj_method *m, const cfj_strset *sets, int nsets,
                                 cfj_vocab *vocab, cfj_profile *profiles, int threads)
{
    size_t total = 0;
    for (int s = 0; s < nsets; s++)
        total += sets[s].n;
    memset(vocab, 0, sizeof(*vocab));
    for (int s = 0; s < nsets; s++)
        memset(&profiles[s], 0, sizeof(profiles[s]));

    int ntasks = threads < 1 ? 1 : threads;
    if ((size_t)ntasks > total / 512 + 1)
        ntasks = (int)(total / 512 + 1);
    task *tasks = (task *)calloc((size_t)ntasks, sizeof(task));
    if (!tasks)
        return STATA_ERR_MEMORY;
    for (int t = 0; t < ntasks; t++)
    {
        tasks[t].first = total * (size_t)t / (size_t)ntasks;
        tasks[t].last = total * (size_t)(t + 1) / (size_t)ntasks;
    }
    int team = ntasks;
    (void)team; /* used by the OpenMP pragmas only */
#pragma omp parallel for schedule(static, 1) num_threads(team) if (team > 1)
    for (int t = 0; t < ntasks; t++)
        run_task(&tasks[t], m, sets, nsets);

    stata_retcode rc = STATA_OK;
    for (int t = 0; t < ntasks && !rc; t++)
        rc = (stata_retcode)tasks[t].rc;

    /* Merge the task dictionaries in chunk order */
    dict global;
    memset(&global, 0, sizeof(global));
    if (!rc && dict_init(&global, 4096))
        rc = STATA_ERR_MEMORY;
    for (int t = 0; t < ntasks && !rc; t++)
    {
        dict *d = &tasks[t].d;
        tasks[t].map = (uint32_t *)malloc((d->n ? d->n : 1) * sizeof(uint32_t));
        if (!tasks[t].map)
        {
            rc = STATA_ERR_MEMORY;
            break;
        }
        for (size_t id = 0; id < d->n; id++)
        {
            long long gid = dict_put(&global, d->bytes + d->off[id], d->len[id], d->ghash[id]);
            if (gid < 0)
            {
                rc = STATA_ERR_MEMORY;
                break;
            }
            tasks[t].map[id] = (uint32_t)gid;
        }
        free(d->hash);
        d->hash = NULL;
        free(d->slot);
        d->slot = NULL;
    }

    /* Profiles: offsets per set, then renumbered pairs */
    for (int s = 0; s < nsets && !rc; s++)
    {
        profiles[s].n = sets[s].n;
        profiles[s].off = (uint64_t *)calloc(sets[s].n + 1, sizeof(uint64_t));
        if (!profiles[s].off)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc)
    {
        /* distinct counts into the set offsets */
        for (int t = 0; t < ntasks; t++)
        {
            for (size_t i = tasks[t].first; i < tasks[t].last; i++)
            {
                size_t r = i;
                int s = 0;
                while (r >= sets[s].n)
                    r -= sets[s].n, s++;
                profiles[s].off[r + 1] = tasks[t].nd[i - tasks[t].first];
            }
        }
        for (int s = 0; s < nsets && !rc; s++)
        {
            for (size_t i = 0; i < sets[s].n; i++)
                profiles[s].off[i + 1] += profiles[s].off[i];
            size_t np = (size_t)profiles[s].off[sets[s].n];
            profiles[s].id = (uint32_t *)malloc((np ? np : 1) * sizeof(uint32_t));
            profiles[s].cnt = (uint32_t *)malloc((np ? np : 1) * sizeof(uint32_t));
            if (!profiles[s].id || !profiles[s].cnt)
                rc = STATA_ERR_MEMORY;
        }
    }
    if (!rc)
    {
#pragma omp parallel for schedule(static, 1) num_threads(team) if (team > 1)
        for (int t = 0; t < ntasks; t++)
        {
            task *tk = &tasks[t];
            size_t at = 0;
            for (size_t i = tk->first; i < tk->last; i++)
            {
                size_t r = i;
                int s = 0;
                while (r >= sets[s].n)
                    r -= sets[s].n, s++;
                size_t nd = tk->nd[i - tk->first];
                uint32_t *id = profiles[s].id + profiles[s].off[r];
                uint32_t *cnt = profiles[s].cnt + profiles[s].off[r];
                for (size_t j = 0; j < nd; j++)
                {
                    id[j] = tk->map[tk->id[at + j]];
                    cnt[j] = tk->cnt[at + j];
                }
                if (t > 0)
                    sort_pairs(id, cnt, nd);
                at += nd;
            }
        }
    }
    if (!rc)
    {
        vocab->n = global.n;
        vocab->off = (uint64_t *)malloc((global.n + 1) * sizeof(uint64_t));
        if (!vocab->off)
            rc = STATA_ERR_MEMORY;
        else
        {
            for (size_t id = 0; id < global.n; id++)
                vocab->off[id] = global.off[id];
            vocab->off[global.n] = global.used;
            vocab->bytes = global.bytes;
            global.bytes = NULL;
        }
    }
    dict_free(&global);
    for (int t = 0; t < ntasks; t++)
    {
        dict_free(&tasks[t].d);
        free(tasks[t].nd);
        free(tasks[t].id);
        free(tasks[t].cnt);
        free(tasks[t].map);
    }
    free(tasks);
    if (rc)
    {
        cfj_vocab_free(vocab);
        for (int s = 0; s < nsets; s++)
            cfj_profile_free(&profiles[s]);
    }
    return rc;
}
