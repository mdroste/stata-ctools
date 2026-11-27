/*
    cfuzzyjoin_match.c
    match: matchit's similarity join of two datasets, and of two columns

    Each text becomes a vector of gram counts (cfuzzyjoin_grams.c). With
    weights, gram g counts w(g) times as much, where w is 1/f, 1/sqrt(f) or
    1/(log(f)+1) of the gram's frequency f: its count over the master and
    using texts, the frequency in wgtfile() (1 for grams the file lacks), or
    both added when stopwordsauto is also given, as matchit does it. Grams
    more frequent than swthreshold() times the number of observations of
    both datasets are stopwords and take no part. With N the weighted dot
    product of two vectors and D1, D2 their squared norms, the scores are

        jaccard    N / sqrt(D1*D2)
        simple     2*N / (D1+D2)
        minsimple  N / min(D1, D2), at most 1

    and 0 when the denominator is not positive. Only pairs that share a gram
    are scored (matchit's index), and those scoring at least threshold() are
    kept. Sums run over grams in vocabulary order, so the results do not
    depend on the thread count; without weights they are exact integers and
    every score is matchit's to the last bit.

    Plugin calls (local __cfj_sim holds the similarity function):

        match_master            varlist: idmaster txtmaster
        match_mweights          varlist: grams freq (the wgtfile() data),
                                after match_master
        match_weights           the same, before match_columns
        match_using <options>   varlist: idusing txtusing. Scores the pairs;
                                saves the pair count and the storage types
                                compress would leave the result with.
        match_store             varlist: the five result variables
        match_columns <options> varlist: text1 text2 score; the score of
                                every observation

    options: score(jaccard|simple|minsimple) weights(none|simple|root|log)
             threshold(#) [wgtfile] [stopwords swthreshold(#)] [diagnose]
             [verbose]
*/

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

enum
{
    SCORE_JACCARD = 0,
    SCORE_SIMPLE,
    SCORE_MINSIMPLE
};
enum
{
    W_NONE = 0,
    W_SIMPLE,
    W_ROOT,
    W_LOG
};

typedef struct
{
    int phase; /* 10 weights file only, 1 master, 2 pairs */
    ctools_filtered_data mfd, ufd, wfd;
    int have_wfile;
    /* pairs */
    uint32_t *pm, *pu;
    double *ps;
    size_t npairs;
} match_state;

static match_state mt;

void cfj_match_cleanup(void)
{
    ctools_filtered_data_free(&mt.mfd);
    ctools_filtered_data_free(&mt.ufd);
    ctools_filtered_data_free(&mt.wfd);
    free(mt.pm);
    free(mt.pu);
    free(mt.ps);
    memset(&mt, 0, sizeof(mt));
    ctools_filtered_data_init(&mt.mfd);
    ctools_filtered_data_init(&mt.ufd);
    ctools_filtered_data_init(&mt.wfd);
}

static ST_retcode fail(ST_retcode rc)
{
    cfj_match_cleanup();
    return rc;
}

typedef struct
{
    cfj_method method;
    int score, weights, stopwords, wfile, diagnose, verbose;
    double threshold, swthreshold;
} options;

static int read_method(cfj_method *m)
{
    char buf[512];
    buf[0] = '\0';
    if (SF_macro_use("___cfj_sim", buf, (int)sizeof(buf)))
        buf[0] = '\0';
    const char *what = NULL;
    if (cfj_method_parse(buf, m, &what))
    {
        cfj_say("{err}similmethod(%s): %s", buf, what ? what : "invalid");
        return 198;
    }
    return 0;
}

static int read_options(const char *args, options *o, int columns)
{
    char text[64];
    memset(o, 0, sizeof(*o));
    if (read_method(&o->method))
        return 198;
    if (!cfj_opt_text(args, "score", text, sizeof(text)) || !strcmp(text, "jaccard"))
        o->score = SCORE_JACCARD;
    else if (!strcmp(text, "simple"))
        o->score = SCORE_SIMPLE;
    else if (!strcmp(text, "minsimple"))
        o->score = SCORE_MINSIMPLE;
    else
        return 198;
    if (!cfj_opt_text(args, "weights", text, sizeof(text)) || !strcmp(text, "none"))
        o->weights = W_NONE;
    else if (!strcmp(text, "simple"))
        o->weights = W_SIMPLE;
    else if (!strcmp(text, "root"))
        o->weights = W_ROOT;
    else if (!strcmp(text, "log"))
        o->weights = W_LOG;
    else
        return 198;
    o->threshold = 0.5;
    if (!columns)
    {
        if (cfj_opt_text(args, "threshold", text, sizeof(text)))
            o->threshold = strtod(text, NULL);
        o->stopwords = cfj_flag(args, "stopwords");
        o->swthreshold = 0.2;
        if (cfj_opt_text(args, "swthreshold", text, sizeof(text)))
            o->swthreshold = strtod(text, NULL);
        o->diagnose = cfj_flag(args, "diagnose");
    }
    o->wfile = cfj_flag(args, "wgtfile");
    o->verbose = cfj_flag(args, "verbose");
    return 0;
}

static double weight_of(int kind, double f)
{
    switch (kind)
    {
    case W_SIMPLE:
        return 1 / f;
    case W_ROOT:
        return 1 / sqrt(f);
    case W_LOG:
        return 1 / (log(f) + 1);
    default:
        return 1;
    }
}

static double score_of(int kind, double num, double d1, double d2)
{
    double denom;
    switch (kind)
    {
    case SCORE_SIMPLE:
        denom = d1 + d2;
        if (denom <= 0)
            return 0;
        /* 2*numerator/denom, read as 2*(numerator/denom) */
        return 2 * (num / denom);
    case SCORE_MINSIMPLE:
    {
        denom = d1 * d2;
        double mn = d1 < d2 ? d1 : d2;
        if (denom <= 0)
            return 0;
        if (num > mn)
            return 1;
        return num / mn;
    }
    default:
        denom = d1 * d2;
        if (denom <= 0)
            return 0;
        return num / sqrt(denom);
    }
}

/* The gram frequencies of the weighting and the stopwords: freq[g] (data
 * counts and/or the wgtfile() frequency), present[g]. data: add the counts
 * of the profiles. */
static stata_retcode frequencies(const cfj_vocab *v, const cfj_profile *pr, int npr, int data,
                                 const ctools_filtered_data *wfd, double **freq_out,
                                 uint8_t **present_out)
{
    double *freq = (double *)calloc(v->n ? v->n : 1, sizeof(double));
    uint8_t *present = (uint8_t *)calloc(v->n ? v->n : 1, 1);
    if (!freq || !present)
    {
        free(freq);
        free(present);
        return STATA_ERR_MEMORY;
    }
    if (wfd)
    {
        cfj_vocab_index ix;
        memset(&ix, 0, sizeof(ix));
        if (cfj_vocab_index_build(v, &ix))
        {
            free(freq);
            free(present);
            return STATA_ERR_MEMORY;
        }
        const stata_variable *g = &wfd->data.vars[0], *f = &wfd->data.vars[1];
        for (size_t i = 0; i < wfd->data.nobs; i++)
        {
            const char *s = g->data.str[i] ? g->data.str[i] : "";
            long long id = cfj_vocab_find(v, &ix, s, strlen(s));
            /* load_weights_to_array: a later row of the same gram wins */
            if (id >= 0)
            {
                freq[id] = f->data.dbl[i];
                present[id] = 1;
            }
        }
        cfj_vocab_index_free(&ix);
    }
    if (data)
    {
        for (int s = 0; s < npr; s++)
        {
            const cfj_profile *p = &pr[s];
            size_t np = (size_t)p->off[p->n];
            for (size_t j = 0; j < np; j++)
            {
                freq[p->id[j]] += (double)p->cnt[j];
                present[p->id[j]] = 1;
            }
        }
    }
    *freq_out = freq;
    *present_out = present;
    return STATA_OK;
}

/* Squared norm of one profile, stopwords left out */
static double norm2(const cfj_profile *p, size_t i, const double *w, const uint8_t *stop)
{
    double d = 0;
    for (uint64_t j = p->off[i]; j < p->off[i + 1]; j++)
    {
        uint32_t g = p->id[j];
        if (stop && stop[g])
            continue;
        if (w)
        {
            double a = (double)p->cnt[j] * w[g];
            d = d + a * a;
        }
        else
        {
            double c = (double)p->cnt[j];
            d = d + c * c;
        }
    }
    return d;
}

/* ----------------------------------------------------------------------------
 * Types that compress leaves
 * -------------------------------------------------------------------------- */

typedef struct
{
    int nonint, any;
    double lo, hi;
} numtype;

static void numtype_add(numtype *t, double x)
{
    if (cfj_missing(x))
        return;
    if (x != floor(x))
        t->nonint = 1;
    if (!t->any || x < t->lo)
        t->lo = x;
    if (!t->any || x > t->hi)
        t->hi = x;
    t->any = 1;
}

/* 1 byte, 2 int, 3 long, 5 double (as compress demotes a double) */
static int numtype_code(const numtype *t)
{
    if (!t->any)
        return 1;
    if (t->nonint)
        return 5;
    if (t->lo >= -127 && t->hi <= 100)
        return 1;
    if (t->lo >= -32767 && t->hi <= 32740)
        return 2;
    if (t->lo >= -2147483647.0 && t->hi <= 2147483620.0)
        return 3;
    return 5;
}

/* ----------------------------------------------------------------------------
 * Diagnose
 * -------------------------------------------------------------------------- */

typedef struct
{
    uint32_t g;
    double key;
} ranked;

static int cmp_ranked(const void *pa, const void *pb)
{
    const ranked *a = (const ranked *)pa, *b = (const ranked *)pb;
    if (a->key > b->key)
        return -1;
    if (a->key < b->key)
        return 1;
    return a->g < b->g ? -1 : a->g > b->g;
}

static void gram_text(const cfj_vocab *v, uint32_t g, char *out, size_t size)
{
    size_t l = (size_t)(v->off[g + 1] - v->off[g]);
    if (l > size - 1)
        l = size - 1;
    memcpy(out, v->bytes + v->off[g], l);
    out[l] = '\0';
    for (size_t i = 0; i < l; i++)
        if (out[i] == '{' || out[i] == '}' || out[i] == '\n' || out[i] == '\r')
            out[i] = '?';
}

static void diagnose(const cfj_vocab *v, const cfj_profile *pm, const cfj_profile *pu)
{
    size_t ng = v->n;
    double *fm = (double *)calloc(ng ? ng : 1, sizeof(double));
    double *fu = (double *)calloc(ng ? ng : 1, sizeof(double));
    ranked *r = (ranked *)malloc((ng ? ng : 1) * sizeof(ranked));
    if (!fm || !fu || !r)
    {
        free(fm);
        free(fu);
        free(r);
        return;
    }
    for (size_t j = 0; j < (size_t)pm->off[pm->n]; j++)
        fm[pm->id[j]] += pm->cnt[j];
    for (size_t j = 0; j < (size_t)pu->off[pu->n]; j++)
        fu[pu->id[j]] += pu->cnt[j];
    const double mN = (double)pm->n, uN = (double)pu->n;
    char gt[64];
    cfj_say("{txt} ");
    cfj_say("{txt}Performing preliminary diagnosis");
    cfj_say("{txt}--------------------------------");
    for (int side = 0; side < 2; side++)
    {
        const double *f = side ? fu : fm;
        const double N = side ? uN : mN;
        size_t n = 0;
        for (size_t g = 0; g < ng; g++)
            if (f[g] > 0)
                r[n].g = (uint32_t)g, r[n].key = f[g], n++;
        qsort(r, n, sizeof(ranked), cmp_ranked);
        cfj_say("{txt} ");
        cfj_say("{txt}Analyzing %s file", side ? "Using" : "Master");
        cfj_say("{txt}List of most frequent grams in %s file:", side ? "Using" : "Master");
        cfj_say("{txt}  %-32s %12s %14s", "grams", "freq", "grams_per_obs");
        for (size_t i = 0; i < n && i < 20; i++)
        {
            gram_text(v, r[i].g, gt, sizeof(gt));
            cfj_say("{res}  %-32s %12.0f %14.4f", gt, f[r[i].g], N > 0 ? f[r[i].g] / N : 0.0);
        }
    }
    double maxspace = 0;
    size_t n = 0;
    for (size_t g = 0; g < ng; g++)
    {
        double space = 0;
        if (fm[g] > 0 && fu[g] > 0)
        {
            space = (fm[g] / mN) * (fu[g] / uN) * 100;
            if (space > 100)
                space = 100;
            r[n].g = (uint32_t)g, r[n].key = fm[g] * fu[g], n++;
        }
        if (space > maxspace)
            maxspace = space;
    }
    double psmax = floor((100 - maxspace) / .01 + .5) * .01;
    cfj_say("{txt} ");
    cfj_say("{txt}Overall diagnosis");
    cfj_say("{txt}Pairs being compared: Master(%.0f) x Using(%.0f) = %.0f", mN, uN, mN * uN);
    cfj_say("{txt}Estimated maximum reduction by indexation (%%):%.2f", psmax);
    cfj_say("{txt}(note: this is an indication, final results may differ)");
    cfj_say("{txt} ");
    cfj_say("{txt}List of grams with greater negative impact to indexation:");
    cfj_say("{txt}(note: values are estimated, final results may differ)");
    qsort(r, n, sizeof(ranked), cmp_ranked);
    cfj_say("{txt}  %-32s %14s %17s %14s", "grams", "crosspairs", "max_common_space",
            "grams_per_obs");
    for (size_t i = 0; i < n && i < 20; i++)
    {
        uint32_t g = r[i].g;
        double space = (fm[g] / mN) * (fu[g] / uN) * 100;
        if (space > 100)
            space = 100;
        gram_text(v, g, gt, sizeof(gt));
        cfj_say("{res}  %-32s %14.0f %17.2f %14.4f", gt, fm[g] * fu[g], space,
                (fm[g] + fu[g]) / (mN + uN));
    }
    cfj_say("{txt} ");
    free(fm);
    free(fu);
    free(r);
}

/* ----------------------------------------------------------------------------
 * Plugin calls
 * -------------------------------------------------------------------------- */

/* A call that starts a match keeps a weights file loaded just before */
static void start_fresh(void)
{
    if (mt.phase == 10 && mt.have_wfile)
    {
        ctools_filtered_data keep = mt.wfd;
        ctools_filtered_data_init(&mt.wfd);
        cfj_match_cleanup();
        mt.wfd = keep;
        mt.have_wfile = 1;
    }
    else
        cfj_match_cleanup();
}

ST_retcode cfj_match_master(const char *args)
{
    (void)args;
    start_fresh();
    if (SF_nvars() != 2)
        return fail(198);
    stata_retcode rc = cfj_load(&mt.mfd, 1, 2, CTOOLS_LOAD_SKIP_IF);
    if (rc)
        return fail(ctools_stata_rc(rc));
    if (mt.mfd.data.vars[0].type == STATA_TYPE_STRING || mt.mfd.data.vars[1].type != STATA_TYPE_STRING)
        return fail(109);
    mt.phase = 1;
    return 0;
}

ST_retcode cfj_match_weights(const char *args, int after_master)
{
    (void)args;
    if (after_master)
    {
        if (mt.phase != 1)
            return fail(198);
        ctools_filtered_data_free(&mt.wfd);
    }
    else
        cfj_match_cleanup();
    if (SF_nvars() != 2)
        return fail(198);
    stata_retcode rc = cfj_load(&mt.wfd, 1, 2, CTOOLS_LOAD_SKIP_IF);
    if (rc)
        return fail(ctools_stata_rc(rc));
    if (mt.wfd.data.vars[0].type != STATA_TYPE_STRING || mt.wfd.data.vars[1].type == STATA_TYPE_STRING)
        return fail(109);
    mt.have_wfile = 1;
    mt.phase = after_master ? 1 : 10;
    return 0;
}

typedef struct
{
    uint32_t *m, *u;
    double *s;
    size_t n, cap;
} block;

static int block_push(block *b, uint32_t m, uint32_t u, double s)
{
    if (b->n == b->cap)
    {
        size_t c = b->cap ? b->cap * 2 : 256;
        uint32_t *a = (uint32_t *)realloc(b->m, c * sizeof(uint32_t));
        if (!a)
            return -1;
        b->m = a;
        a = (uint32_t *)realloc(b->u, c * sizeof(uint32_t));
        if (!a)
            return -1;
        b->u = a;
        double *d = (double *)realloc(b->s, c * sizeof(double));
        if (!d)
            return -1;
        b->s = d;
        b->cap = c;
    }
    b->m[b->n] = m;
    b->u[b->n] = u;
    b->s[b->n] = s;
    b->n++;
    return 0;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

static void sort_u64(uint64_t *a, size_t n)
{
    if (n > 24)
    {
        qsort(a, n, sizeof(uint64_t), cmp_u64);
        return;
    }
    for (size_t i = 1; i < n; i++)
    {
        uint64_t v = a[i];
        size_t j = i;
        while (j && a[j - 1] > v)
        {
            a[j] = a[j - 1];
            j--;
        }
        a[j] = v;
    }
}

ST_retcode cfj_match_using(const char *args)
{
    options o;
    if (mt.phase != 1 || SF_nvars() != 2)
        return fail(198);
    if (read_options(args, &o, 0))
        return fail(198);
    if (o.wfile && !mt.have_wfile)
        return fail(198);
    const double t0 = ctools_timer_seconds();
    stata_retcode rc = cfj_load(&mt.ufd, 1, 2, CTOOLS_LOAD_SKIP_IF);
    if (rc)
        return fail(ctools_stata_rc(rc));
    if (mt.ufd.data.vars[0].type == STATA_TYPE_STRING || mt.ufd.data.vars[1].type != STATA_TYPE_STRING)
        return fail(109);
    const size_t nm = mt.mfd.data.nobs, nu = mt.ufd.data.nobs;
    const int threads = ctools_get_max_threads();

    cfj_strset sets[2] = {{(const char *const *)mt.mfd.data.vars[1].data.str, nm},
                          {(const char *const *)mt.ufd.data.vars[1].data.str, nu}};
    cfj_vocab voc;
    cfj_profile prof[2];
    rc = cfj_profiles_build(&o.method, sets, 2, &voc, prof, threads);
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return fail(ctools_stata_rc(rc));
    }
    const double t1 = ctools_timer_seconds();
    if (o.diagnose)
        diagnose(&voc, &prof[0], &prof[1]);

    /* weights and stopwords */
    double *freq = NULL, *w = NULL;
    uint8_t *present = NULL, *stop = NULL;
    const int need_freq = o.weights != W_NONE || o.stopwords;
    if (need_freq)
    {
        const int data = (o.weights != W_NONE && !o.wfile) || o.stopwords;
        rc = frequencies(&voc, prof, 2, data, o.wfile ? &mt.wfd : NULL, &freq, &present);
    }
    if (!rc && o.weights != W_NONE)
    {
        w = (double *)malloc((voc.n ? voc.n : 1) * sizeof(double));
        if (!w)
            rc = STATA_ERR_MEMORY;
        else
            for (size_t g = 0; g < voc.n; g++)
                w[g] = weight_of(o.weights, present[g] ? freq[g] : 1.0);
    }
    size_t nstop = 0;
    if (!rc && o.stopwords)
    {
        stop = (uint8_t *)calloc(voc.n ? voc.n : 1, 1);
        if (!stop)
            rc = STATA_ERR_MEMORY;
        else
        {
            const double limit = o.swthreshold * (double)(nm + nu);
            for (size_t g = 0; g < voc.n; g++)
                if (present[g] && freq[g] > limit)
                    stop[g] = 1, nstop++;
        }
    }

    /* index of the using profiles: postings by gram */
    uint64_t *pstart = NULL;
    uint32_t *prow = NULL, *pcnt = NULL;
    double *d2 = NULL;
    if (!rc)
    {
        pstart = (uint64_t *)calloc(voc.n + 1, sizeof(uint64_t));
        d2 = (double *)malloc((nu ? nu : 1) * sizeof(double));
        if (!pstart || !d2)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc)
    {
        const cfj_profile *pu = &prof[1];
        size_t np = (size_t)pu->off[nu];
        for (size_t j = 0; j < np; j++)
            if (!stop || !stop[pu->id[j]])
                pstart[pu->id[j] + 1]++;
        for (size_t g = 0; g < voc.n; g++)
            pstart[g + 1] += pstart[g];
        size_t total = (size_t)pstart[voc.n];
        prow = (uint32_t *)malloc((total ? total : 1) * sizeof(uint32_t));
        pcnt = (uint32_t *)malloc((total ? total : 1) * sizeof(uint32_t));
        uint64_t *fill = (uint64_t *)malloc((voc.n + 1) * sizeof(uint64_t));
        if (!prow || !pcnt || !fill)
            rc = STATA_ERR_MEMORY;
        else
        {
            memcpy(fill, pstart, (voc.n + 1) * sizeof(uint64_t));
            for (size_t u = 0; u < nu; u++)
            {
                for (uint64_t j = pu->off[u]; j < pu->off[u + 1]; j++)
                {
                    uint32_t g = pu->id[j];
                    if (stop && stop[g])
                        continue;
                    prow[fill[g]] = (uint32_t)u;
                    pcnt[fill[g]] = pu->cnt[j];
                    fill[g]++;
                }
            }
        }
        free(fill);
        if (!rc)
        {
#pragma omp parallel for schedule(static) if (nu > 10000)
            for (size_t u = 0; u < nu; u++)
                d2[u] = norm2(pu, u, w, stop);
        }
    }
    const double t2 = ctools_timer_seconds();

    /* score: blocks of master rows */
    const size_t bsize = 256;
    const size_t nblocks = (nm + bsize - 1) / bsize;
    block *blocks = NULL;
    double space = 0;
    if (!rc)
    {
        blocks = (block *)calloc(nblocks ? nblocks : 1, sizeof(block));
        if (!blocks)
            rc = STATA_ERR_MEMORY;
    }
    if (!rc && nm && nu)
    {
        int team = ctools_workspace_threads(threads, nblocks,
                                            nu * (sizeof(double) + 2 * sizeof(uint32_t)));
        atomic_int failed = 0;
        const cfj_profile *pm = &prof[0];
#pragma omp parallel num_threads(team) if (team > 1) reduction(+ : space)
        {
            double *acc = (double *)malloc(nu * sizeof(double));
            uint32_t *mark = (uint32_t *)calloc(nu, sizeof(uint32_t));
            uint32_t *touched = (uint32_t *)malloc(nu * sizeof(uint32_t));
            /* the pairs kept for one master row, as (using row << 32 | slot) */
            uint64_t *kept = (uint64_t *)malloc(nu * sizeof(uint64_t));
            double *kscore = (double *)malloc(nu * sizeof(double));
            if (!acc || !mark || !touched || !kept || !kscore)
                atomic_store(&failed, STATA_ERR_MEMORY);
#pragma omp for schedule(dynamic, 1)
            for (size_t b = 0; b < nblocks; b++)
            {
                if (atomic_load_explicit(&failed, memory_order_relaxed))
                    continue;
                if (ctools_cancel_requested())
                {
                    atomic_store(&failed, STATA_ERR_CANCELLED);
                    continue;
                }
                block *bl = &blocks[b];
                const size_t lo = b * bsize, hi = lo + bsize < nm ? lo + bsize : nm;
                for (size_t i = lo; i < hi; i++)
                {
                    const uint32_t stamp = (uint32_t)i + 1;
                    size_t nt = 0;
                    for (uint64_t j = pm->off[i]; j < pm->off[i + 1]; j++)
                    {
                        const uint32_t g = pm->id[j];
                        if (stop && stop[g])
                            continue;
                        const double cm = (double)pm->cnt[j];
                        const double wsq = w ? w[g] * w[g] : 1.0;
                        for (uint64_t q = pstart[g]; q < pstart[g + 1]; q++)
                        {
                            const uint32_t u = prow[q];
                            const double add = w ? (cm * (double)pcnt[q]) * wsq : cm * (double)pcnt[q];
                            if (mark[u] != stamp)
                            {
                                mark[u] = stamp;
                                acc[u] = add;
                                touched[nt++] = u;
                            }
                            else
                                acc[u] = acc[u] + add;
                        }
                    }
                    space += (double)nt;
                    if (!nt)
                        continue;
                    const double d1 = norm2(pm, i, w, stop);
                    /* score in the order the rows were met; only the kept
                     * pairs are put in using order */
                    size_t nk = 0;
                    for (size_t q = 0; q < nt; q++)
                    {
                        const uint32_t u = touched[q];
                        const double s = score_of(o.score, acc[u], d1, d2[u]);
                        if (s >= o.threshold)
                        {
                            kscore[nk] = s;
                            kept[nk] = ((uint64_t)u << 32) | (uint64_t)nk;
                            nk++;
                        }
                    }
                    sort_u64(kept, nk);
                    for (size_t q = 0; q < nk; q++)
                    {
                        if (block_push(bl, (uint32_t)i, (uint32_t)(kept[q] >> 32),
                                       kscore[(uint32_t)kept[q]]))
                        {
                            atomic_store(&failed, STATA_ERR_MEMORY);
                            break;
                        }
                    }
                }
            }
            free(acc);
            free(mark);
            free(touched);
            free(kept);
            free(kscore);
        }
        rc = atomic_load(&failed);
    }
    const double t3 = ctools_timer_seconds();

    /* concatenate the blocks */
    size_t npairs = 0;
    if (!rc)
    {
        for (size_t b = 0; b < nblocks; b++)
            npairs += blocks[b].n;
        mt.pm = (uint32_t *)malloc((npairs ? npairs : 1) * sizeof(uint32_t));
        mt.pu = (uint32_t *)malloc((npairs ? npairs : 1) * sizeof(uint32_t));
        mt.ps = (double *)malloc((npairs ? npairs : 1) * sizeof(double));
        if (!mt.pm || !mt.pu || !mt.ps)
            rc = STATA_ERR_MEMORY;
        else
        {
            size_t at = 0;
            for (size_t b = 0; b < nblocks; b++)
            {
                memcpy(mt.pm + at, blocks[b].m, blocks[b].n * sizeof(uint32_t));
                memcpy(mt.pu + at, blocks[b].u, blocks[b].n * sizeof(uint32_t));
                memcpy(mt.ps + at, blocks[b].s, blocks[b].n * sizeof(double));
                at += blocks[b].n;
            }
            mt.npairs = npairs;
        }
    }
    for (size_t b = 0; blocks && b < nblocks; b++)
    {
        free(blocks[b].m);
        free(blocks[b].u);
        free(blocks[b].s);
    }
    free(blocks);
    free(pstart);
    free(prow);
    free(pcnt);
    free(d2);
    free(freq);
    free(present);
    free(w);
    free(stop);
    const size_t nvocab = voc.n;
    cfj_vocab_free(&voc);
    cfj_profile_free(&prof[0]);
    cfj_profile_free(&prof[1]);
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return fail(ctools_stata_rc(rc));
    }

    /* the types compress would choose */
    numtype tim = {0}, tiu = {0}, ts = {0};
    size_t lm = 0, lu = 0;
    {
        const stata_variable *im = &mt.mfd.data.vars[0], *xm = &mt.mfd.data.vars[1];
        const stata_variable *iu = &mt.ufd.data.vars[0], *xu = &mt.ufd.data.vars[1];
        uint8_t *seen_m = (uint8_t *)calloc(nm ? nm : 1, 1);
        uint8_t *seen_u = (uint8_t *)calloc(nu ? nu : 1, 1);
        if (!seen_m || !seen_u)
        {
            free(seen_m);
            free(seen_u);
            return fail(920);
        }
        for (size_t q = 0; q < npairs; q++)
        {
            numtype_add(&ts, mt.ps[q]);
            seen_m[mt.pm[q]] = 1;
            seen_u[mt.pu[q]] = 1;
        }
        for (size_t i = 0; i < nm; i++)
        {
            if (!seen_m[i])
                continue;
            numtype_add(&tim, im->data.dbl[i]);
            size_t l = xm->data.str[i] ? strlen(xm->data.str[i]) : 0;
            if (l > lm)
                lm = l;
        }
        for (size_t i = 0; i < nu; i++)
        {
            if (!seen_u[i])
                continue;
            numtype_add(&tiu, iu->data.dbl[i]);
            size_t l = xu->data.str[i] ? strlen(xu->data.str[i]) : 0;
            if (l > lu)
                lu = l;
        }
        free(seen_m);
        free(seen_u);
    }
    mt.phase = 2;
    ST_retcode src = 0;
    const double saved = (nm && nu) ? 100 - trunc(space * (100.0 / (double)nu / (double)nm)) : 0;
    src = src ? src : ctools_scal_save("_cfuzzyjoin_npairs", (double)npairs);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_type_idm", (double)numtype_code(&tim));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_type_idu", (double)numtype_code(&tiu));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_type_score", (double)numtype_code(&ts));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_len_txtm", (double)(lm ? lm : 1));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_len_txtu", (double)(lu ? lu : 1));
    src = src ? src : ctools_scal_save("_cfuzzyjoin_saved", saved);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_pairs", space);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_ngrams", (double)nvocab);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_nstop", (double)nstop);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_grams", t1 - t0);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_index", t2 - t1);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_compute", t3 - t2);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_threads", (double)threads);
    if (src)
        return fail(src);
    ctools_verbose(CFJ_MODULE, o.verbose,
                   "match: %zu x %zu obs, %zu grams (%zu stopwords), %.0f pairs scored, %zu kept; "
                   "grams %.4fs, index %.4fs, score %.4fs",
                   nm, nu, nvocab, nstop, space, npairs, t1 - t0, t2 - t1, t3 - t2);
    return 0;
}

ST_retcode cfj_match_store(const char *args)
{
    (void)args;
    if (mt.phase != 2 || SF_nvars() != 5 || (size_t)SF_nobs() != mt.npairs)
        return fail(198);
    const stata_variable *im = &mt.mfd.data.vars[0], *xm = &mt.mfd.data.vars[1];
    const stata_variable *iu = &mt.ufd.data.vars[0], *xu = &mt.ufd.data.vars[1];
    atomic_int err = 0;
    const size_t n = mt.npairs;
    int team = ctools_get_max_threads();
    if ((size_t)team > n / 4096 + 1)
        team = (int)(n / 4096 + 1);
#pragma omp parallel for schedule(static) num_threads(team) if (team > 1)
    for (size_t q = 0; q < n; q++)
    {
        const ST_int obs = (ST_int)(q + 1);
        const uint32_t m = mt.pm[q], u = mt.pu[q];
        int w = CFJ_VSTORE(1, obs, im->data.dbl[m]);
        w |= SF_sstore(2, obs, xm->data.str[m] ? xm->data.str[m] : (char *)"");
        w |= CFJ_VSTORE(3, obs, iu->data.dbl[u]);
        w |= SF_sstore(4, obs, xu->data.str[u] ? xu->data.str[u] : (char *)"");
        w |= CFJ_VSTORE(5, obs, mt.ps[q]);
        if (w)
            atomic_store(&err, 1);
    }
    int bad = atomic_load(&err);
    cfj_match_cleanup();
    return bad ? 459 : 0;
}

ST_retcode cfj_match_columns(const char *args)
{
    options o;
    start_fresh();
    if (SF_nvars() != 3 || read_options(args, &o, 1) || (o.wfile && !mt.have_wfile))
        return fail(198);
    const double t0 = ctools_timer_seconds();
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    stata_retcode rc = cfj_load(&fd, 1, 2, CTOOLS_LOAD_SKIP_IF);
    if (rc)
    {
        ctools_filtered_data_free(&fd);
        return fail(ctools_stata_rc(rc));
    }
    if (fd.data.vars[0].type != STATA_TYPE_STRING || fd.data.vars[1].type != STATA_TYPE_STRING)
    {
        ctools_filtered_data_free(&fd);
        return fail(109);
    }
    const size_t n = fd.data.nobs;
    const int threads = ctools_get_max_threads();
    cfj_strset sets[2] = {{(const char *const *)fd.data.vars[0].data.str, n},
                          {(const char *const *)fd.data.vars[1].data.str, n}};
    cfj_vocab voc;
    cfj_profile prof[2];
    rc = cfj_profiles_build(&o.method, sets, 2, &voc, prof, threads);
    double *freq = NULL, *w = NULL, *out = NULL;
    uint8_t *present = NULL;
    if (!rc && o.weights != W_NONE)
    {
        rc = frequencies(&voc, prof, 2, !o.wfile, o.wfile ? &mt.wfd : NULL, &freq, &present);
        if (!rc)
        {
            w = (double *)malloc((voc.n ? voc.n : 1) * sizeof(double));
            if (!w)
                rc = STATA_ERR_MEMORY;
            else
                for (size_t g = 0; g < voc.n; g++)
                    w[g] = weight_of(o.weights, present[g] ? freq[g] : 1.0);
        }
    }
    if (!rc)
    {
        out = (double *)malloc((n ? n : 1) * sizeof(double));
        if (!out)
            rc = STATA_ERR_MEMORY;
    }
    const double t1 = ctools_timer_seconds();
    if (!rc)
    {
        const cfj_profile *a = &prof[0], *b = &prof[1];
#pragma omp parallel for schedule(static) if (n > 20000)
        for (size_t i = 0; i < n; i++)
        {
            double num = 0;
            uint64_t p = a->off[i], q = b->off[i];
            const uint64_t pe = a->off[i + 1], qe = b->off[i + 1];
            while (p < pe && q < qe)
            {
                if (a->id[p] < b->id[q])
                    p++;
                else if (a->id[p] > b->id[q])
                    q++;
                else
                {
                    uint32_t g = a->id[p];
                    double c = (double)a->cnt[p] * (double)b->cnt[q];
                    num = num + (w ? c * (w[g] * w[g]) : c);
                    p++, q++;
                }
            }
            out[i] = score_of(o.score, num, norm2(a, i, w, NULL), norm2(b, i, w, NULL));
        }
        rc = ctools_store_filtered(out, n, 3, fd.obs_map);
    }
    const double t2 = ctools_timer_seconds();
    free(freq);
    free(present);
    free(w);
    free(out);
    cfj_vocab_free(&voc);
    cfj_profile_free(&prof[0]);
    cfj_profile_free(&prof[1]);
    ctools_filtered_data_free(&fd);
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CFJ_MODULE);
        return fail(ctools_stata_rc(rc));
    }
    ST_retcode src = ctools_scal_save("_cfuzzyjoin_t_grams", t1 - t0);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_t_compute", t2 - t1);
    src = src ? src : ctools_scal_save("_cfuzzyjoin_threads", (double)threads);
    if (src)
        return fail(src);
    ctools_verbose(CFJ_MODULE, o.verbose, "match columns: %zu obs; grams %.4fs, score %.4fs", n,
                   t1 - t0, t2 - t1);
    cfj_match_cleanup();
    return 0;
}
