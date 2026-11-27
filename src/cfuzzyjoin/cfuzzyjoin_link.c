/*
    cfuzzyjoin_link.c
    link: record linkage with reclink's matching scores

    reclink (M. Blasnik, version 1.7) compares each master record with the
    distinct value combinations of the using matching variables. Records that
    equal a combination on every variable are perfect matches (score 1) and
    take that combination out of the search. The others are scored against
    every remaining combination that passes the blocking rules: required()
    variables equal, and at least one or-blocking variable equal. The
    combinations with the highest score, when it exceeds minscore(), are the
    matches (all of them when they tie), and each brings every using record
    that has it.

    The score of a pair follows reclink's Stata statements exactly, including
    the rounding of its float variables after each statement and Stata's
    reading of a*b/c as a*(b/c):

      bigram (string variables): the share of the master string's bigrams
        (by position) that occur anywhere in the using string, scaled by
        2/(slen+mlen-2); strings shorter than three bytes score
        min/(mlen+slen-1) when one contains the other, else 0; Winkler's
        adjustment adds wink*((1-b)/10) for each of the first one to four
        bytes on which the strings agree; equal strings score 1, an empty
        using string 0.
      M += wm*(.5*b^2 + (b > minbigram ? biscale*(b-minbigram)^(1/3) : 0))
      NM += wnm*(b < minbigram-.2 ? 1 : b < minbigram ? 1-b^2 : (1-b)^2)
            when the using value is not missing
      NM += .3*wnm when exactly one of the master value (as text) and the
            using value is missing
      score = M/(M+NM), 0 when undefined.

    Plugin calls:

        link_master            varlist: idmaster, then the k matching
                               variables. Checks that idmaster identifies the
                               observations.
        link_exclude           varlist: the exclude() file's idmaster and
                               idusing variables.
        link_using k=<k> p=<p> varlist: idusing, the k using matching
                               variables, then p payload variables. Locals
                               __cfj_kind (b bigram, e exact string, n
                               numeric), __cfj_req and __cfj_orb (0/1),
                               __cfj_wm and __cfj_wnm hold one entry per
                               matching variable; minscore(), minbigram() and
                               biscale() the numbers as Stata wrote them.
                               Finds the matches and lays out the output.
        link_layout            varlist: a new variable that receives the
                               number of output rows of each master record.
        link_store             varlist: the master record of each row
                               (after expand), the score, idusing, the
                               payload variables, _merge, the output order,
                               then (flag urow) the using observation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdatomic.h>

#include "cfuzzyjoin_internal.h"
#include "ctools_config.h"
#include "ctools_cancel.h"
#include "ctools_runtime.h"
#include "ctools_spi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define KIND_BIGRAM 'b'
#define KIND_EXACT 'e'
#define KIND_NUMERIC 'n'

typedef struct
{
    uint32_t m;  /* master row (0-based) */
    uint32_t u;  /* using row (0-based); UINT32_MAX unmatched */
    float score;
} entry;

static struct
{
    int phase; /* 1 master loaded, 2 using matched, 3 laid out */
    int k, p;
    ctools_filtered_data mfd, ufd;
    /* exclude() ids */
    ctools_filtered_data xfd;
    int have_exclude;
    /* output */
    entry *out;
    size_t nout, nmatched_rows;
    uint32_t *first; /* first entry of each master row, nm + 1 */
    size_t nm;
    double perfect, matched_masters, using_left, master_left, using_n, master_n;
    double candidates;
} lk;

void cfj_link_cleanup(void)
{
    ctools_filtered_data_free(&lk.mfd);
    ctools_filtered_data_free(&lk.ufd);
    ctools_filtered_data_free(&lk.xfd);
    free(lk.out);
    free(lk.first);
    memset(&lk, 0, sizeof(lk));
    ctools_filtered_data_init(&lk.mfd);
    ctools_filtered_data_init(&lk.ufd);
    ctools_filtered_data_init(&lk.xfd);
}

static ST_retcode fail(ST_retcode rc)
{
    cfj_link_cleanup();
    return rc;
}

/* ----------------------------------------------------------------------------
 * Values
 * -------------------------------------------------------------------------- */

static int is_str(const stata_variable *v)
{
    return v->type == STATA_TYPE_STRING;
}

static const char *sval(const stata_variable *v, size_t row)
{
    const char *s = v->data.str[row];
    return s ? s : "";
}

static uint64_t hash_value(const stata_variable *v, size_t row)
{
    if (is_str(v))
    {
        const char *s = sval(v, row);
        return cfj_hash_bytes(s, strlen(s));
    }
    double x = v->data.dbl[row];
    if (x == 0)
        x = 0; /* -0 equals 0 */
    uint64_t b;
    memcpy(&b, &x, 8);
    b ^= b >> 31;
    b *= 0x9E3779B97F4A7C15ull;
    return b ^ (b >> 29);
}

static int same_value(const stata_variable *a, size_t ra, const stata_variable *b, size_t rb)
{
    if (is_str(a))
        return !strcmp(sval(a, ra), sval(b, rb));
    return a->data.dbl[ra] == b->data.dbl[rb];
}

/* Stata's missing(): "" or any missing code */
static int value_missing(const stata_variable *v, size_t row)
{
    if (is_str(v))
        return sval(v, row)[0] == '\0';
    return cfj_missing(v->data.dbl[row]);
}

/* ----------------------------------------------------------------------------
 * Index of combinations by the values of some variables
 * -------------------------------------------------------------------------- */

typedef struct
{
    int nv;
    int vars[64];
    size_t mask;
    uint32_t *slot; /* bucket + 1 */
    uint64_t *shash;
    uint32_t *rep;   /* representative using row of each bucket */
    uint32_t *start; /* nb + 1 */
    uint32_t *list;  /* combinations, ascending */
    size_t nb;
} vindex;

static void vindex_free(vindex *x)
{
    free(x->slot);
    free(x->shash);
    free(x->rep);
    free(x->start);
    free(x->list);
    memset(x, 0, sizeof(*x));
}

static uint64_t hash_row(const ctools_filtered_data *fd, const int *vars, int nv, size_t row)
{
    uint64_t h = 0x84222325cbf29ce4ull;
    for (int i = 0; i < nv; i++)
    {
        h ^= hash_value(&fd->data.vars[1 + vars[i]], row) + 0x9E3779B97F4A7C15ull + (h << 6) +
             (h >> 2);
    }
    return h;
}

static int same_row(const ctools_filtered_data *a, size_t ra, const ctools_filtered_data *b,
                    size_t rb, const int *vars, int nv)
{
    for (int i = 0; i < nv; i++)
        if (!same_value(&a->data.vars[1 + vars[i]], ra, &b->data.vars[1 + vars[i]], rb))
            return 0;
    return 1;
}

/* Group the given rows (combination representatives, or all rows) by vars.
 * bucket_of[i] receives the bucket of item i; the buckets are numbered in
 * order of first appearance. */
static int vindex_build(vindex *x, const ctools_filtered_data *fd, const uint32_t *rows, size_t n,
                        uint32_t *bucket_of)
{
    size_t c = 64;
    while (c < n * 2 + 2)
        c *= 2;
    x->mask = c - 1;
    x->slot = (uint32_t *)calloc(c, sizeof(uint32_t));
    x->shash = (uint64_t *)malloc(c * sizeof(uint64_t));
    x->rep = (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    uint32_t *own = bucket_of ? NULL : (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    if (!bucket_of)
        bucket_of = own;
    if (!x->slot || !x->shash || !x->rep || !bucket_of)
    {
        free(own);
        return -1;
    }
    x->nb = 0;
    for (size_t i = 0; i < n; i++)
    {
        size_t row = rows ? rows[i] : i;
        uint64_t h = hash_row(fd, x->vars, x->nv, row);
        size_t at = (size_t)h & x->mask;
        uint32_t b = 0;
        while (x->slot[at])
        {
            uint32_t cand = x->slot[at] - 1;
            if (x->shash[at] == h && same_row(fd, x->rep[cand], fd, row, x->vars, x->nv))
            {
                b = cand + 1;
                break;
            }
            at = (at + 1) & x->mask;
        }
        if (!b)
        {
            x->rep[x->nb] = (uint32_t)row;
            x->slot[at] = (uint32_t)x->nb + 1;
            x->shash[at] = h;
            b = (uint32_t)++x->nb;
        }
        bucket_of[i] = b - 1;
    }
    x->start = (uint32_t *)calloc(x->nb + 1, sizeof(uint32_t));
    x->list = (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    if (!x->start || !x->list)
    {
        free(own);
        return -1;
    }
    for (size_t i = 0; i < n; i++)
        x->start[bucket_of[i] + 1]++;
    for (size_t b = 0; b < x->nb; b++)
        x->start[b + 1] += x->start[b];
    uint32_t *fill = (uint32_t *)malloc((x->nb + 1) * sizeof(uint32_t));
    if (!fill)
    {
        free(own);
        return -1;
    }
    memcpy(fill, x->start, (x->nb + 1) * sizeof(uint32_t));
    for (size_t i = 0; i < n; i++)
        x->list[fill[bucket_of[i]]++] = (uint32_t)i;
    free(fill);
    free(own);
    return 0;
}

/* The bucket whose values equal master row mrow's, or -1 */
static long long vindex_find(const vindex *x, const ctools_filtered_data *ufd,
                             const ctools_filtered_data *mfd, size_t mrow)
{
    uint64_t h = hash_row(mfd, x->vars, x->nv, mrow);
    size_t at = (size_t)h & x->mask;
    while (x->slot[at])
    {
        uint32_t b = x->slot[at] - 1;
        if (x->shash[at] == h && same_row(ufd, x->rep[b], mfd, mrow, x->vars, x->nv))
            return b;
        at = (at + 1) & x->mask;
    }
    return -1;
}

/* ----------------------------------------------------------------------------
 * Scores
 * -------------------------------------------------------------------------- */

typedef struct
{
    char kind;
    int req, orb;
    double wm, wnm;
} varspec;

typedef struct
{
    double minscore, minbigram, biscale;
} params;

/* The using side of a bigram variable: each combination's distinct bigrams */
typedef struct
{
    uint32_t *start; /* ncomb + 1 */
    uint16_t *code;
} ubigrams;

static inline uint16_t bigram_code(const unsigned char *s)
{
    return (uint16_t)((s[0] << 8) | s[1]);
}

/* strpos(hay, needle) > 0 (needle may be empty: true) */
static int contains(const char *hay, size_t lh, const char *needle, size_t ln)
{
    if (ln == 0)
        return 1;
    if (ln > lh)
        return 0;
    for (size_t i = 0; i + ln <= lh; i++)
        if (hay[i] == needle[0] && !memcmp(hay + i, needle, ln))
            return 1;
    return 0;
}

/* bigram1: the master string ms (mlen bytes, its bigram counts in mcnt)
 * against the using string us */
static float bigram_score(const char *ms, size_t mlen, const uint16_t *mcnt, const char *us,
                          size_t slen, const uint16_t *ucode, uint32_t nucode)
{
    if (slen == 0)
        return 0.0f;
    if (slen == mlen && !memcmp(ms, us, mlen))
        return 1.0f;
    float gen = 0.0f;
    if (mlen > 2)
    {
        uint32_t c = 0;
        for (uint32_t j = 0; j < nucode; j++)
            c += mcnt[ucode[j]];
        gen = (float)c;
        /* gen*2/(slen+mlen-2), read by Stata as gen*(2/(slen+mlen-2)) */
        gen = (float)((double)gen * (2.0 / ((double)slen + (double)mlen - 2.0)));
    }
    if (mlen < 3 || slen < 3)
    {
        if (!contains(ms, mlen, us, slen) && !contains(us, slen, ms, mlen))
            gen = 0.0f;
        else
            gen = (float)((double)(mlen < slen ? mlen : slen) / ((double)mlen + (double)slen - 1.0));
    }
    int wink = 0;
    for (size_t i = 1; i <= 4; i++)
    {
        size_t a = i < slen ? i : slen, b = i < mlen ? i : mlen;
        if (a == b && !memcmp(us, ms, a))
            wink++;
    }
    /* gen + wink*(1-gen)/10, read as gen + wink*((1-gen)/10) */
    gen = (float)((double)gen + (double)wink * ((1.0 - (double)gen) / 10.0));
    return gen;
}

/* The text reclink's bigram routine receives for a master value: it travels
 * through the string option match() of a syntax command, which trims its
 * spaces, joins its words with single spaces (tabs stay as they are) and drops
 * double quotes. out holds at least strlen(s) + 1 bytes. */
static size_t syntax_text(const char *s, char *out)
{
    size_t n = 0;
    int pending = 0;
    for (const char *p = s; *p; p++)
    {
        if (*p == '"')
            continue;
        if (*p == ' ')
        {
            pending = n > 0;
            continue;
        }
        if (pending)
            out[n++] = ' ';
        pending = 0;
        out[n++] = *p;
    }
    out[n] = '\0';
    return n;
}

/* The master record, prepared once */
typedef struct
{
    char **mstr;     /* per bigram variable: the master text as bigram1 sees it */
    size_t *mlen;
    uint16_t **mcnt; /* per bigram variable: 65536 counts */
    uint16_t *touched;
    size_t ntouched, cap;
} mprep;

typedef struct
{
    const ctools_filtered_data *mfd, *ufd;
    const uint32_t *comb_rep;
    const varspec *spec;
    const params *par;
    int k;
    const ubigrams *ub; /* per variable (bigram ones) */
    const int *bslot;   /* bigram table slot of each variable, -1 */
} ctx;

static float pair_score(const ctx *c, const mprep *mp, size_t mrow, uint32_t comb)
{
    const size_t urow = c->comb_rep[comb];
    float M = 0.0f, NM = 0.0f;
    for (int i = 0; i < c->k; i++)
    {
        const stata_variable *mv = &c->mfd->data.vars[1 + i];
        const stata_variable *uv = &c->ufd->data.vars[1 + i];
        const varspec *sp = &c->spec[i];
        /* the master value as reclink's macro holds it: "" for an empty
         * string, "." (not empty) for a missing number */
        const int mempty = is_str(mv) && sval(mv, mrow)[0] == '\0';
        const int umiss = value_missing(uv, urow);
        if (!mempty)
        {
            if (!sp->req)
            {
                float bi;
                if (sp->kind == KIND_BIGRAM)
                {
                    const int b = c->bslot[i];
                    const char *us = sval(uv, urow);
                    const ubigrams *ub = &c->ub[i];
                    bi = bigram_score(mp->mstr[b], mp->mlen[b], mp->mcnt[b], us, strlen(us),
                                      ub->code + ub->start[comb],
                                      ub->start[comb + 1] - ub->start[comb]);
                }
                else
                    bi = same_value(mv, mrow, uv, urow) ? 1.0f : 0.0f;
                const double b = (double)bi;
                const double mb = c->par->minbigram;
                double add = .5 * (b * b) + (b > mb ? c->par->biscale * pow(b - mb, 1.0 / 3.0) : 0.0);
                M = (float)((double)M + sp->wm * add);
                if (!umiss)
                {
                    double nm = b < (mb - .2) ? 1.0 : (b < mb ? 1.0 - b * b : (1.0 - b) * (1.0 - b));
                    NM = (float)((double)NM + sp->wnm * nm);
                }
            }
            else
                M = (float)((double)M + sp->wm);
        }
        if (mempty + umiss == 1)
            NM = (float)((double)NM + sp->wnm * 0.3);
    }
    double s = (double)M / ((double)M + (double)NM);
    float score = (float)s;
    if (score != score)
        score = 0.0f;
    return score;
}

/* ----------------------------------------------------------------------------
 * Plugin calls
 * -------------------------------------------------------------------------- */

static int check_ids(const stata_variable *v, size_t n, const char *what)
{
    /* missing values, then duplicates, as isid reports them */
    for (size_t i = 0; i < n; i++)
    {
        if (value_missing(v, i))
        {
            cfj_say("{err}variable %s should never be missing", what);
            return 459;
        }
    }
    vindex x;
    memset(&x, 0, sizeof(x));
    ctools_filtered_data tmp;
    memset(&tmp, 0, sizeof(tmp));
    /* index on the id alone: wrap it as a one-variable dataset view */
    stata_variable vars[2];
    memset(vars, 0, sizeof(vars));
    vars[1] = *v;
    tmp.data.vars = vars;
    tmp.data.nvars = 2;
    tmp.data.nobs = n;
    x.nv = 1;
    x.vars[0] = 0;
    if (vindex_build(&x, &tmp, NULL, n, NULL))
    {
        vindex_free(&x);
        return 920;
    }
    int dup = x.nb != n;
    vindex_free(&x);
    if (dup)
    {
        cfj_say("{err}variable %s does not uniquely identify the observations", what);
        return 459;
    }
    return 0;
}

ST_retcode cfj_link_master(const char *args)
{
    char idname[64] = "idmaster";
    cfj_link_cleanup();
    cfj_opt_text(args, "id", idname, sizeof(idname));
    if (SF_nvars() < 2)
        return fail(198);
    stata_retcode rc = cfj_load(&lk.mfd, 1, SF_nvars(), CTOOLS_LOAD_SKIP_IF);
    if (rc)
        return fail(ctools_stata_rc(rc));
    lk.nm = lk.mfd.data.nobs;
    int irc = check_ids(&lk.mfd.data.vars[0], lk.nm, idname);
    if (irc)
        return fail(irc);
    lk.k = SF_nvars() - 1;
    lk.phase = 1;
    return 0;
}

ST_retcode cfj_link_exclude(const char *args)
{
    (void)args;
    if (lk.phase != 1 || SF_nvars() != 2)
        return fail(198);
    ctools_filtered_data_free(&lk.xfd);
    stata_retcode rc = cfj_load(&lk.xfd, 1, 2, CTOOLS_LOAD_SKIP_IF);
    if (rc)
        return fail(ctools_stata_rc(rc));
    if (is_str(&lk.xfd.data.vars[0]) != is_str(&lk.mfd.data.vars[0]))
    {
        cfj_say("{err}type mismatch: idmaster in the exclude() file and in the master data");
        return fail(106);
    }
    lk.have_exclude = 1;
    return 0;
}

/* Is the id at row r of v among the exclude() ids of column col (rows where
 * both ids are present)? Uses an index built once. */
typedef struct
{
    vindex x;
    ctools_filtered_data view;
    stata_variable vars[2];
    uint32_t *rows;
} idset;

static int idset_build(idset *s, const ctools_filtered_data *xfd, int col)
{
    memset(s, 0, sizeof(*s));
    size_t n = xfd->data.nobs, m = 0;
    s->rows = (uint32_t *)malloc((n ? n : 1) * sizeof(uint32_t));
    if (!s->rows)
        return -1;
    for (size_t i = 0; i < n; i++)
        if (!value_missing(&xfd->data.vars[0], i) && !value_missing(&xfd->data.vars[1], i))
            s->rows[m++] = (uint32_t)i;
    s->vars[1] = xfd->data.vars[col];
    s->view.data.vars = s->vars;
    s->view.data.nvars = 2;
    s->view.data.nobs = n;
    s->x.nv = 1;
    s->x.vars[0] = 0;
    return vindex_build(&s->x, &s->view, s->rows, m, NULL);
}

static int idset_has(const idset *s, const stata_variable *v, size_t row)
{
    stata_variable vars[2];
    memset(vars, 0, sizeof(vars));
    vars[1] = *v;
    ctools_filtered_data view;
    memset(&view, 0, sizeof(view));
    view.data.vars = vars;
    view.data.nvars = 2;
    view.data.nobs = row + 1;
    return vindex_find(&s->x, &s->view, &view, row) >= 0;
}

static void idset_free(idset *s)
{
    vindex_free(&s->x);
    free(s->rows);
    memset(s, 0, sizeof(*s));
}

/* Read n tokens of local `name` */
static int read_list(const char *name, char *buf, size_t size, char **tok, int n)
{
    char local[64];
    snprintf(local, sizeof(local), "_%s", name);
    buf[0] = '\0';
    if (SF_macro_use(local, buf, (int)size))
        return -1;
    int k = 0;
    char *p = buf;
    while (*p && k < n)
    {
        while (*p == ' ')
            p++;
        if (!*p)
            break;
        tok[k++] = p;
        while (*p && *p != ' ')
            p++;
        if (*p)
            *p++ = '\0';
    }
    return k == n ? 0 : -1;
}

static int cmp_ctx_numeric;
static const stata_variable *cmp_ctx_id;

static int cmp_entries(const void *pa, const void *pb)
{
    const entry *a = (const entry *)pa, *b = (const entry *)pb;
    if (cmp_ctx_numeric)
    {
        double x = cmp_ctx_id->data.dbl[a->u], y = cmp_ctx_id->data.dbl[b->u];
        if (x < y)
            return -1;
        if (x > y)
            return 1;
    }
    else
    {
        int c = strcmp(sval(cmp_ctx_id, a->u), sval(cmp_ctx_id, b->u));
        if (c)
            return c;
    }
    return a->m < b->m ? -1 : a->m > b->m;
}

ST_retcode cfj_link_using(const char *args)
{
    int k = 0, p = 0;
    if (lk.phase != 1 || cfj_opt_int(args, "k", &k) || cfj_opt_int(args, "p", &p) || k != lk.k ||
        k < 1 || k > 64 || p < 0 || SF_nvars() != 1 + k + p)
        return fail(198);
    const int verbose = cfj_flag(args, "verbose");
    char tx1[64], tx2[64], tx3[64];
    params par;
    if (!cfj_opt_text(args, "minscore", tx1, sizeof(tx1)) ||
        !cfj_opt_text(args, "minbigram", tx2, sizeof(tx2)) ||
        !cfj_opt_text(args, "biscale", tx3, sizeof(tx3)))
        return fail(198);
    par.minscore = strtod(tx1, NULL);
    par.minbigram = strtod(tx2, NULL);
    par.biscale = strtod(tx3, NULL);

    /* per-variable specification */
    varspec spec[64];
    {
        size_t size = (size_t)k * 40 + 64;
        char *buf = (char *)malloc(size);
        char *tok[64];
        if (!buf)
            return fail(920);
        int ok = 1;
        if (read_list("__cfj_kind", buf, size, tok, k))
            ok = 0;
        for (int i = 0; ok && i < k; i++)
            spec[i].kind = tok[i][0];
        if (ok && read_list("__cfj_req", buf, size, tok, k))
            ok = 0;
        for (int i = 0; ok && i < k; i++)
            spec[i].req = tok[i][0] == '1';
        if (ok && read_list("__cfj_orb", buf, size, tok, k))
            ok = 0;
        for (int i = 0; ok && i < k; i++)
            spec[i].orb = tok[i][0] == '1';
        if (ok && read_list("__cfj_wm", buf, size, tok, k))
            ok = 0;
        for (int i = 0; ok && i < k; i++)
            spec[i].wm = strtod(tok[i], NULL);
        if (ok && read_list("__cfj_wnm", buf, size, tok, k))
            ok = 0;
        for (int i = 0; ok && i < k; i++)
            spec[i].wnm = strtod(tok[i], NULL);
        free(buf);
        if (!ok)
            return fail(198);
    }

    const double t0 = ctools_timer_seconds();
    stata_retcode rc = cfj_load(&lk.ufd, 1, 1 + k + p, CTOOLS_LOAD_SKIP_IF);
    if (rc)
        return fail(ctools_stata_rc(rc));
    const size_t nu = lk.ufd.data.nobs, nm = lk.nm;
    lk.p = p;
    for (int i = 0; i < k; i++)
    {
        if (is_str(&lk.mfd.data.vars[1 + i]) != is_str(&lk.ufd.data.vars[1 + i]))
        {
            cfj_say("{err}type mismatch: matching variable %d is a string in one dataset and "
                    "numeric in the other",
                    i + 1);
            return fail(106);
        }
        if (spec[i].kind == KIND_NUMERIC && is_str(&lk.mfd.data.vars[1 + i]))
            spec[i].kind = KIND_EXACT;
    }
    if (lk.have_exclude && is_str(&lk.xfd.data.vars[1]) != is_str(&lk.ufd.data.vars[0]))
    {
        cfj_say("{err}type mismatch: idusing in the exclude() file and in the using data");
        return fail(106);
    }

    /* Exclusions */
    uint8_t *mskip = (uint8_t *)calloc(nm ? nm : 1, 1);
    uint32_t *urows = (uint32_t *)malloc((nu ? nu : 1) * sizeof(uint32_t));
    if (!mskip || !urows)
    {
        free(mskip);
        free(urows);
        return fail(920);
    }
    size_t nuse = 0, mleft = nm;
    if (lk.have_exclude)
    {
        idset xs_m, xs_u;
        if (idset_build(&xs_m, &lk.xfd, 0) || idset_build(&xs_u, &lk.xfd, 1))
        {
            idset_free(&xs_m);
            idset_free(&xs_u);
            free(mskip);
            free(urows);
            return fail(920);
        }
        for (size_t i = 0; i < nm; i++)
            if (idset_has(&xs_m, &lk.mfd.data.vars[0], i))
                mskip[i] = 1, mleft--;
        for (size_t i = 0; i < nu; i++)
            if (!idset_has(&xs_u, &lk.ufd.data.vars[0], i))
                urows[nuse++] = (uint32_t)i;
        idset_free(&xs_m);
        idset_free(&xs_u);
    }
    else
    {
        for (size_t i = 0; i < nu; i++)
            urows[nuse++] = (uint32_t)i;
    }
    lk.master_n = (double)nm;
    lk.using_n = (double)nu;
    lk.master_left = (double)mleft;
    lk.using_left = (double)nuse;

    /* idusing identifies the remaining using observations (isid) */
    {
        stata_variable vars[2];
        memset(vars, 0, sizeof(vars));
        vars[1] = lk.ufd.data.vars[0];
        ctools_filtered_data view;
        memset(&view, 0, sizeof(view));
        view.data.vars = vars;
        view.data.nvars = 2;
        view.data.nobs = nu;
        char idname[64] = "idusing";
        cfj_opt_text(args, "id", idname, sizeof(idname));
        for (size_t i = 0; i < nuse; i++)
        {
            if (value_missing(&lk.ufd.data.vars[0], urows[i]))
            {
                cfj_say("{err}variable %s should never be missing", idname);
                free(mskip);
                free(urows);
                return fail(459);
            }
        }
        vindex x;
        memset(&x, 0, sizeof(x));
        x.nv = 1;
        x.vars[0] = 0;
        int brc = vindex_build(&x, &view, urows, nuse, NULL);
        int dup = !brc && x.nb != nuse;
        vindex_free(&x);
        if (brc || dup)
        {
            if (dup)
                cfj_say("{err}variable %s does not uniquely identify the observations", idname);
            free(mskip);
            free(urows);
            return fail(brc ? 920 : 459);
        }
    }

    /* Combinations of the matching variables */
    vindex comb;
    memset(&comb, 0, sizeof(comb));
    comb.nv = k;
    for (int i = 0; i < k; i++)
        comb.vars[i] = i;
    uint32_t *comb_of = (uint32_t *)malloc((nuse ? nuse : 1) * sizeof(uint32_t));
    if (!comb_of || vindex_build(&comb, &lk.ufd, urows, nuse, comb_of))
    {
        free(comb_of);
        vindex_free(&comb);
        free(mskip);
        free(urows);
        return fail(920);
    }
    const size_t ncomb = comb.nb;
    /* comb.rep holds using rows; comb.list holds positions in urows */

    /* Perfect matches */
    uint32_t *perfect = (uint32_t *)malloc((nm ? nm : 1) * sizeof(uint32_t));
    uint8_t *gone = (uint8_t *)calloc(ncomb ? ncomb : 1, 1);
    rc = (!perfect || !gone) ? STATA_ERR_MEMORY : STATA_OK;
    size_t nperfect = 0;
    if (!rc)
    {
#pragma omp parallel for schedule(static) reduction(+ : nperfect) if (nm > 20000)
        for (size_t i = 0; i < nm; i++)
        {
            perfect[i] = UINT32_MAX;
            if (mskip[i])
                continue;
            long long b = vindex_find(&comb, &lk.ufd, &lk.mfd, i);
            if (b >= 0)
            {
                perfect[i] = (uint32_t)b;
                nperfect++;
            }
        }
        for (size_t i = 0; i < nm; i++)
            if (perfect[i] != UINT32_MAX)
                gone[perfect[i]] = 1;
    }
    lk.perfect = (double)nperfect;

    /* Blocking indexes over the combinations */
    int orbvars[64], reqvars[64], norb = 0, nreq = 0;
    for (int i = 0; i < k; i++)
    {
        if (spec[i].req)
            reqvars[nreq++] = i;
        else if (spec[i].orb)
            orbvars[norb++] = i;
    }
    vindex req;
    memset(&req, 0, sizeof(req));
    vindex orb[64];
    memset(orb, 0, sizeof(orb));
    uint32_t *crow = NULL; /* representative using row of each combination */
    if (!rc)
    {
        crow = (uint32_t *)malloc((ncomb ? ncomb : 1) * sizeof(uint32_t));
        if (!crow)
            rc = STATA_ERR_MEMORY;
        else
            memcpy(crow, comb.rep, ncomb * sizeof(uint32_t));
    }
    if (!rc && nreq)
    {
        req.nv = nreq;
        memcpy(req.vars, reqvars, sizeof(int) * (size_t)nreq);
        if (vindex_build(&req, &lk.ufd, crow, ncomb, NULL))
            rc = STATA_ERR_MEMORY;
    }
    for (int o = 0; o < norb && !rc; o++)
    {
        orb[o].nv = 1;
        orb[o].vars[0] = orbvars[o];
        if (vindex_build(&orb[o], &lk.ufd, crow, ncomb, NULL))
            rc = STATA_ERR_MEMORY;
    }

    /* Bigrams of the using strings, per combination */
    ubigrams ub[64];
    memset(ub, 0, sizeof(ub));
    int bslot[64], nbig = 0;
    for (int i = 0; i < k; i++)
        bslot[i] = spec[i].kind == KIND_BIGRAM ? nbig++ : -1;
    for (int i = 0; i < k && !rc; i++)
    {
        if (bslot[i] < 0)
            continue;
        const stata_variable *uv = &lk.ufd.data.vars[1 + i];
        ub[i].start = (uint32_t *)calloc(ncomb + 1, sizeof(uint32_t));
        if (!ub[i].start)
        {
            rc = STATA_ERR_MEMORY;
            break;
        }
        size_t total = 0;
        for (size_t c = 0; c < ncomb; c++)
        {
            size_t l = strlen(sval(uv, crow[c]));
            total += l > 1 ? l - 1 : 0;
        }
        ub[i].code = (uint16_t *)malloc((total ? total : 1) * sizeof(uint16_t));
        if (!ub[i].code)
        {
            rc = STATA_ERR_MEMORY;
            break;
        }
        size_t at = 0;
        for (size_t c = 0; c < ncomb; c++)
        {
            const unsigned char *s = (const unsigned char *)sval(uv, crow[c]);
            size_t l = strlen((const char *)s);
            size_t from = at;
            for (size_t j = 0; j + 1 < l; j++)
                ub[i].code[at++] = bigram_code(s + j);
            /* distinct codes: sort the run, drop repeats */
            size_t nrun = at - from;
            uint16_t *run = ub[i].code + from;
            for (size_t a = 1; a < nrun; a++)
            {
                uint16_t v = run[a];
                size_t b = a;
                while (b && run[b - 1] > v)
                {
                    run[b] = run[b - 1];
                    b--;
                }
                run[b] = v;
            }
            size_t w = 0;
            for (size_t a = 0; a < nrun; a++)
                if (!a || run[a] != run[w - 1])
                    run[w++] = run[a];
            at = from + w;
            ub[i].start[c + 1] = (uint32_t)at;
        }
    }
    const double t1 = ctools_timer_seconds();

    /* Fuzzy matches */
    ctx cx;
    cx.mfd = &lk.mfd;
    cx.ufd = &lk.ufd;
    cx.comb_rep = crow;
    cx.spec = spec;
    cx.par = &par;
    cx.k = k;
    cx.ub = ub;
    cx.bslot = bslot;

    const int threads = ctools_get_max_threads();
    int team = threads;
    if ((size_t)team > nm / 16 + 1)
        team = (int)(nm / 16 + 1);
    /* results: per master, the best combinations */
    typedef struct
    {
        uint32_t m, c;
        float s;
    } hit;
    hit **hits = NULL;
    size_t *nhits = NULL;
    double cand_total = 0;
    if (!rc)
    {
        hits = (hit **)calloc((size_t)team, sizeof(hit *));
        nhits = (size_t *)calloc((size_t)team, sizeof(size_t));
        if (!hits || !nhits)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc)
    {
        atomic_int failed = 0;
#pragma omp parallel num_threads(team) if (team > 1) reduction(+ : cand_total)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t cap = 256, n = 0;
            hit *mine = (hit *)malloc(cap * sizeof(hit));
            mprep mp;
            memset(&mp, 0, sizeof(mp));
            mp.mcnt = (uint16_t **)calloc((size_t)(nbig ? nbig : 1), sizeof(uint16_t *));
            mp.mstr = (char **)calloc((size_t)(nbig ? nbig : 1), sizeof(char *));
            mp.mlen = (size_t *)calloc((size_t)(nbig ? nbig : 1), sizeof(size_t));
            uint32_t *stamp = (uint32_t *)calloc(ncomb ? ncomb : 1, sizeof(uint32_t));
            uint32_t *cand = (uint32_t *)malloc((ncomb ? ncomb : 1) * sizeof(uint32_t));
            uint32_t *best = (uint32_t *)malloc((ncomb ? ncomb : 1) * sizeof(uint32_t));
            int ok = mine && mp.mcnt && mp.mstr && mp.mlen && stamp && cand && best;
            for (int v = 0; ok && v < k; v++)
            {
                if (bslot[v] < 0)
                    continue;
                const int b = bslot[v];
                mp.mcnt[b] = (uint16_t *)calloc(65536, sizeof(uint16_t));
                mp.mstr[b] = (char *)malloc(lk.mfd.data.vars[1 + v].str_maxlen + 2048);
                ok = mp.mcnt[b] != NULL && mp.mstr[b] != NULL;
            }
            if (!ok)
                atomic_store(&failed, STATA_ERR_MEMORY);
#pragma omp for schedule(dynamic, 16)
            for (size_t i = 0; i < nm; i++)
            {
                if (!ok || atomic_load_explicit(&failed, memory_order_relaxed))
                    continue;
                if ((i & 255) == 0 && ctools_cancel_requested())
                {
                    atomic_store(&failed, STATA_ERR_CANCELLED);
                    continue;
                }
                if (mskip[i] || perfect[i] != UINT32_MAX)
                    continue;
                /* candidates */
                size_t nc = 0;
                const uint32_t mark = (uint32_t)i + 1;
                if (nreq)
                {
                    long long b = vindex_find(&req, &lk.ufd, &lk.mfd, i);
                    if (b < 0)
                        continue;
                    for (uint32_t q = req.start[b]; q < req.start[b + 1]; q++)
                    {
                        uint32_t c = req.list[q];
                        if (gone[c])
                            continue;
                        if (norb)
                        {
                            int any = 0;
                            for (int o = 0; o < norb && !any; o++)
                            {
                                const stata_variable *mv = &lk.mfd.data.vars[1 + orbvars[o]];
                                if (is_str(mv) && !sval(mv, i)[0])
                                    continue;
                                any = same_value(mv, i, &lk.ufd.data.vars[1 + orbvars[o]], crow[c]);
                            }
                            if (!any)
                                continue;
                        }
                        cand[nc++] = c;
                    }
                }
                else if (norb)
                {
                    for (int o = 0; o < norb; o++)
                    {
                        const stata_variable *mv = &lk.mfd.data.vars[1 + orbvars[o]];
                        if (is_str(mv) && !sval(mv, i)[0])
                            continue;
                        long long b = vindex_find(&orb[o], &lk.ufd, &lk.mfd, i);
                        if (b < 0)
                            continue;
                        for (uint32_t q = orb[o].start[b]; q < orb[o].start[b + 1]; q++)
                        {
                            uint32_t c = orb[o].list[q];
                            if (gone[c] || stamp[c] == mark)
                                continue;
                            stamp[c] = mark;
                            cand[nc++] = c;
                        }
                    }
                }
                else
                {
                    for (uint32_t c = 0; c < ncomb; c++)
                        if (!gone[c])
                            cand[nc++] = c;
                }
                if (!nc)
                    continue;
                cand_total += (double)nc;
                /* the master's bigram counts */
                mp.ntouched = 0;
                for (int v = 0; v < k; v++)
                {
                    if (bslot[v] < 0)
                        continue;
                    const int b = bslot[v];
                    mp.mlen[b] = syntax_text(sval(&lk.mfd.data.vars[1 + v], i), mp.mstr[b]);
                    const unsigned char *s = (const unsigned char *)mp.mstr[b];
                    size_t l = mp.mlen[b];
                    uint16_t *t = mp.mcnt[b];
                    for (size_t j = 0; j + 1 < l; j++)
                        t[bigram_code(s + j)]++;
                }
                float top = -1.0f;
                size_t nbest = 0;
                for (size_t q = 0; q < nc; q++)
                {
                    float s = pair_score(&cx, &mp, i, cand[q]);
                    if (s > top)
                    {
                        top = s;
                        nbest = 0;
                    }
                    if (s == top)
                        best[nbest++] = cand[q];
                }
                for (int v = 0; v < k; v++)
                {
                    if (bslot[v] < 0)
                        continue;
                    const unsigned char *s = (const unsigned char *)mp.mstr[bslot[v]];
                    size_t l = mp.mlen[bslot[v]];
                    uint16_t *t = mp.mcnt[bslot[v]];
                    for (size_t j = 0; j + 1 < l; j++)
                        t[bigram_code(s + j)] = 0;
                }
                if (!((double)top > par.minscore))
                    continue;
                for (size_t q = 0; q < nbest; q++)
                {
                    if (n == cap)
                    {
                        hit *more = (hit *)realloc(mine, 2 * cap * sizeof(hit));
                        if (!more)
                        {
                            atomic_store(&failed, STATA_ERR_MEMORY);
                            break;
                        }
                        mine = more;
                        cap *= 2;
                    }
                    mine[n].m = (uint32_t)i;
                    mine[n].c = best[q];
                    mine[n].s = top;
                    n++;
                }
            }
            hits[tid] = mine;
            nhits[tid] = n;
            for (int b = 0; mp.mcnt && b < nbig; b++)
                free(mp.mcnt[b]);
            for (int b = 0; mp.mstr && b < nbig; b++)
                free(mp.mstr[b]);
            free(mp.mcnt);
            free(mp.mstr);
            free(mp.mlen);
            free(stamp);
            free(cand);
            free(best);
        }
        rc = atomic_load(&failed);
    }
    lk.candidates = cand_total;
    const double t2 = ctools_timer_seconds();

    /* Entries: every using record of each matched combination */
    if (!rc)
    {
        uint32_t *nent = (uint32_t *)calloc(nm + 1, sizeof(uint32_t));
        if (!nent)
            rc = STATA_ERR_MEMORY;
        /* combination -> its using rows, ascending */
        uint32_t *cstart = comb.start, *clist = comb.list;
        size_t nmatched = 0;
        for (size_t i = 0; !rc && i < nm; i++)
            if (perfect[i] != UINT32_MAX)
                nent[i + 1] += cstart[perfect[i] + 1] - cstart[perfect[i]];
        for (int t = 0; !rc && t < team; t++)
            for (size_t h = 0; h < nhits[t]; h++)
                nent[hits[t][h].m + 1] += cstart[hits[t][h].c + 1] - cstart[hits[t][h].c];
        if (!rc)
        {
            lk.first = (uint32_t *)malloc((nm + 1) * sizeof(uint32_t));
            size_t total = 0;
            for (size_t i = 0; i < nm; i++)
            {
                size_t e = nent[i + 1];
                if (e)
                    nmatched++;
                total += e ? e : 1;
            }
            lk.out = (entry *)malloc((total ? total : 1) * sizeof(entry));
            if (!lk.first || !lk.out)
                rc = STATA_ERR_MEMORY;
            else
            {
                /* master order first; sorted into the output order below */
                size_t at = 0;
                for (size_t i = 0; i < nm; i++)
                {
                    lk.first[i] = (uint32_t)at;
                    at += nent[i + 1] ? nent[i + 1] : 1;
                }
                lk.first[nm] = (uint32_t)at;
                uint32_t *fill = nent; /* reuse: entries written so far */
                memset(fill, 0, (nm + 1) * sizeof(uint32_t));
                for (size_t i = 0; i < nm; i++)
                {
                    if (perfect[i] == UINT32_MAX)
                        continue;
                    uint32_t c = perfect[i];
                    for (uint32_t q = cstart[c]; q < cstart[c + 1]; q++)
                    {
                        entry *e = &lk.out[lk.first[i] + fill[i]++];
                        e->m = (uint32_t)i;
                        e->u = urows[clist[q]];
                        e->score = 1.0f;
                    }
                }
                for (int t = 0; t < team; t++)
                {
                    for (size_t h = 0; h < nhits[t]; h++)
                    {
                        const uint32_t i = hits[t][h].m, c = hits[t][h].c;
                        for (uint32_t q = cstart[c]; q < cstart[c + 1]; q++)
                        {
                            entry *e = &lk.out[lk.first[i] + fill[i]++];
                            e->m = i;
                            e->u = urows[clist[q]];
                            e->score = hits[t][h].s;
                        }
                    }
                }
                for (size_t i = 0; i < nm; i++)
                {
                    if (fill[i])
                        continue;
                    entry *e = &lk.out[lk.first[i]];
                    e->m = (uint32_t)i;
                    e->u = UINT32_MAX;
                    e->score = 0.0f;
                }
                lk.nout = total;
                lk.nmatched_rows = 0;
                for (size_t e = 0; e < total; e++)
                    if (lk.out[e].u != UINT32_MAX)
                        lk.nmatched_rows++;
                /* each master's entries by idusing, so copies of a master
                 * row take its matches in that order */
                cmp_ctx_numeric = !is_str(&lk.ufd.data.vars[0]);
                cmp_ctx_id = &lk.ufd.data.vars[0];
                for (size_t i = 0; i < nm; i++)
                {
                    size_t a = lk.first[i], b = lk.first[i + 1];
                    if (b - a > 1)
                        qsort(lk.out + a, b - a, sizeof(entry), cmp_entries);
                }
            }
            lk.matched_masters = (double)nmatched;
        }
        free(nent);
    }

    for (int t = 0; hits && t < team; t++)
        free(hits[t]);
    free(hits);
    free(nhits);
    for (int i = 0; i < k; i++)
    {
        free(ub[i].start);
        free(ub[i].code);
    }
    for (int o = 0; o < norb; o++)
        vindex_free(&orb[o]);
    vindex_free(&req);
    vindex_free(&comb);
    free(comb_of);
    free(crow);
    free(perfect);
    free(gone);
    free(mskip);
    free(urows);
    const double t3 = ctools_timer_seconds();
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return fail(ctools_stata_rc(rc));
    }
    lk.phase = 2;
    ST_retcode src = 0;
    src = src ? src : ctools_scal_save("_cfuzzyjoin_perfect", lk.perfect);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_matched", lk.matched_masters);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_master_left", lk.master_left);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_using_left", lk.using_left);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_using_N", lk.using_n);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_ncomb", (double)ncomb);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_tomatch", (double)(mleft - nperfect));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_nout", (double)lk.nout);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_expand", (double)(lk.nout > nm));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_pairs", lk.candidates);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_index", t1 - t0);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_compute", t2 - t1);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_layout", t3 - t2);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_threads", (double)threads);
    if (src)
        return fail(src);
    ctools_verbose(CFJ_MODULE, verbose,
                   "link: %zu master, %zu using obs, %zu combinations, %.0f perfect, %.0f pairs "
                   "scored; index %.4fs, score %.4fs, layout %.4fs",
                   nm, nu, ncomb, lk.perfect, lk.candidates, t1 - t0, t2 - t1, t3 - t2);
    return 0;
}

ST_retcode cfj_link_layout(const char *args)
{
    (void)args;
    if (lk.phase != 2 || SF_nvars() != 1 || (size_t)SF_nobs() != lk.nm)
        return fail(198);
    double *cnt = (double *)malloc((lk.nm ? lk.nm : 1) * sizeof(double));
    if (!cnt)
        return fail(920);
    for (size_t i = 0; i < lk.nm; i++)
        cnt[i] = (double)(lk.first[i + 1] - lk.first[i]);
    stata_retcode rc = STATA_OK;
    for (size_t i = 0; i < lk.nm && !rc; i++)
        if (CFJ_VSTORE(1, (ST_int)(i + 1), cnt[i]))
            rc = STATA_ERR_STATA_WRITE;
    free(cnt);
    if (rc)
        return fail(ctools_stata_rc(rc));
    lk.phase = 3;
    return 0;
}

ST_retcode cfj_link_store(const char *args)
{
    const int want_urow = cfj_flag(args, "urow");
    const int nv = 1 + 1 + 1 + lk.p + 1 + 1 + (want_urow ? 1 : 0);
    if (lk.phase != 3 || SF_nvars() != nv || (size_t)SF_nobs() != lk.nout)
    {
        cfj_say("{err}cfuzzyjoin link: store expects %d variables and %zu observations "
                "(got %d and %d, phase %d)",
                nv, lk.nout, (int)SF_nvars(), (int)SF_nobs(), lk.phase);
        return fail(198);
    }
    const size_t nrow = lk.nout;
    /* the entry of each row: originals are rows 1..nm, copies follow */
    uint32_t *ent = (uint32_t *)malloc((nrow ? nrow : 1) * sizeof(uint32_t));
    uint32_t *seen = (uint32_t *)calloc(lk.nm ? lk.nm : 1, sizeof(uint32_t));
    uint32_t *pos = (uint32_t *)malloc((nrow ? nrow : 1) * sizeof(uint32_t));
    if (!ent || !seen || !pos)
    {
        free(ent);
        free(seen);
        free(pos);
        return fail(920);
    }
    stata_retcode rc = STATA_OK;
    for (size_t r = 0; r < nrow && !rc; r++)
    {
        size_t m = r;
        if (r >= lk.nm)
        {
            double v;
            if (CFJ_VDATA(1, (ST_int)(r + 1), &v) || cfj_missing(v) || v < 1 || v > (double)lk.nm)
            {
                rc = STATA_ERR_STATA_READ;
                break;
            }
            m = (size_t)v - 1;
        }
        uint32_t e = lk.first[m] + seen[m]++;
        if (e >= lk.first[m + 1])
        {
            cfj_say("{err}cfuzzyjoin link: observation %zu repeats master record %zu too often",
                    r + 1, m + 1);
            rc = STATA_ERR_INVALID_INPUT;
            break;
        }
        ent[r] = e;
    }
    /* output order: matched entries by idusing (then master row), then the
     * unmatched master rows in master order */
    if (!rc)
    {
        uint32_t *order = (uint32_t *)malloc((nrow ? nrow : 1) * sizeof(uint32_t));
        if (!order)
            rc = STATA_ERR_MEMORY;
        else
        {
            size_t nmat = 0;
            for (size_t e = 0; e < nrow; e++)
                if (lk.out[e].u != UINT32_MAX)
                    order[nmat++] = (uint32_t)e;
            /* sort matched entry indices by (idusing, master) */
            entry *tmp = (entry *)malloc((nmat ? nmat : 1) * sizeof(entry));
            if (!tmp)
                rc = STATA_ERR_MEMORY;
            else
            {
                for (size_t q = 0; q < nmat; q++)
                {
                    tmp[q] = lk.out[order[q]];
                    tmp[q].m = order[q]; /* sort key: entry index for ties */
                }
                cmp_ctx_numeric = !is_str(&lk.ufd.data.vars[0]);
                cmp_ctx_id = &lk.ufd.data.vars[0];
                qsort(tmp, nmat, sizeof(entry), cmp_entries);
                for (size_t q = 0; q < nmat; q++)
                    pos[tmp[q].m] = (uint32_t)q + 1;
                size_t at = nmat;
                for (size_t e = 0; e < nrow; e++)
                    if (lk.out[e].u == UINT32_MAX)
                        pos[e] = (uint32_t)++at;
            }
            free(tmp);
            free(order);
        }
    }
    /* write */
    if (!rc)
    {
        const int vscore = 2, vid = 3, vpay = 4, vmerge = 4 + lk.p, vord = 5 + lk.p,
                  vurow = 6 + lk.p;
        const stata_variable *uid = &lk.ufd.data.vars[0];
        atomic_int err = 0;
        int team = ctools_get_max_threads();
        if ((size_t)team > nrow / 4096 + 1)
            team = (int)(nrow / 4096 + 1);
#pragma omp parallel for schedule(static) num_threads(team) if (team > 1)
        for (size_t r = 0; r < nrow; r++)
        {
            const entry *e = &lk.out[ent[r]];
            const ST_int obs = (ST_int)(r + 1);
            int w = 0;
            w |= CFJ_VSTORE(vord, obs, (double)pos[ent[r]]);
            if (e->u == UINT32_MAX)
            {
                w |= CFJ_VSTORE(vmerge, obs, 1.0);
            }
            else
            {
                w |= CFJ_VSTORE(vscore, obs, (double)e->score);
                if (is_str(uid))
                    w |= SF_sstore(vid, obs, (char *)sval(uid, e->u));
                else
                    w |= CFJ_VSTORE(vid, obs, uid->data.dbl[e->u]);
                for (int j = 0; j < lk.p; j++)
                {
                    const stata_variable *pv = &lk.ufd.data.vars[1 + lk.k + j];
                    if (is_str(pv))
                        w |= SF_sstore(vpay + j, obs, (char *)sval(pv, e->u));
                    else
                        w |= CFJ_VSTORE(vpay + j, obs, pv->data.dbl[e->u]);
                }
                w |= CFJ_VSTORE(vmerge, obs, 3.0);
                if (want_urow)
                    w |= CFJ_VSTORE(vurow, obs, (double)e->u + 1.0);
            }
            if (w)
                atomic_store(&err, 1);
        }
        if (atomic_load(&err))
            rc = STATA_ERR_STATA_WRITE;
    }
    free(ent);
    free(seen);
    free(pos);
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return fail(ctools_stata_rc(rc));
    }
    cfj_link_cleanup();
    return 0;
}
