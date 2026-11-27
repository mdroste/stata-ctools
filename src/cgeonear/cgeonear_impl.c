/* cgeonear: geonear's nearest neighbors from a k-d tree.
 *
 * Plugin calls, in order (the runtime refuses any other sequence):
 *   base     master data, variables id lat lon: the base locations and the
 *            location of each observation (kept)
 *   using    the neighbor data, variables id lat lon: the neighbor
 *            locations (kept)
 *   compute  wide: variables nid1 km1 ... nidK kmK, written for every
 *            observation, and everything is freed. long: no variables;
 *            the pairs are kept and __cgeonear_nout says how many
 *   store    long only, variables baseid nborid dist of a dataset of
 *            __cgeonear_nout observations: writes the pairs
 *   clear    frees everything
 *   dist     geodist's form (cgeonear_dist.c), which keeps nothing
 *
 * Results are those of comparing every base with every neighbor: for each
 * base location, the neighbors ranked by distance and then by identifier
 * (geonear's documented order of ties), the nearest nearcount() of them and,
 * in long form, every one within(). geonear reaches the same pairs by
 * splitting the base locations into regions; it never drops a neighbor that
 * could be nearer than the ones it keeps. The search is in
 * cgeonear_search.c; this file reads and writes Stata's data. */
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cgeonear_impl.h"
#include "cgeonear_internal.h"
#include "cgeonear_search.h"
#include "cgeo_math.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include "ctools_types.h"

typedef struct
{
    int stage; /* 1 base read, 2 neighbors read, 3 long pairs kept */
    cgeo_locs base, nbor;
    uint32_t *obs_loc; /* base observation -> base location, or CGEO_NONE */
    size_t nobs;
    cgeo_result pairs; /* long form, between compute and store */
} cgeonear_state;

static cgeonear_state S;

void cgeonear_cleanup_state(void)
{
    cgeo_locs_free(&S.base);
    cgeo_locs_free(&S.nbor);
    free(S.obs_loc);
    cgeo_result_free(&S.pairs);
    memset(&S, 0, sizeof(S));
}

static void save_fail(cgeo_fail fail)
{
    static const char *names[] = {"", "notconst", "noobs", "string", "lat", "lon"};
    SF_macro_save("_cgeonear_fail", (char *)names[fail]);
}

static int read_cfg(const char *args, cgeo_cfg *c, int *timing)
{
    memset(c, 0, sizeof(*c));
    c->wide = ctools_parse_bool_option(args, "wide");
    c->sphere = ctools_parse_bool_option(args, "sphere");
    c->miles = ctools_parse_bool_option(args, "miles");
    c->missing = SV_missval;
    *timing = ctools_parse_bool_option(args, "verbose");
    int k = 0, limit = -1;
    if (ctools_parse_int_option(args, "nearcount", &k) != CTOOLS_PARSE_OK || k < 0)
        return 198;
    ctools_parse_status ls = ctools_parse_int_option(args, "limit", &limit);
    if (ls == CTOOLS_PARSE_INVALID || (ls == CTOOLS_PARSE_OK && limit < 0))
        return 198;
    c->k = (size_t)k;
    c->limit = limit < 0 ? SIZE_MAX : (size_t)limit;
    if (c->wide && !c->k)
        return 198;
    double r = 0, a = 0, f = 0, w = 0;
    if (c->sphere)
    {
        if (SF_scal_use("__cgeonear_R", &r) || !(r > 0))
            return 198;
        c->r = r;
    }
    else
    {
        if (SF_scal_use("__cgeonear_a", &a) || SF_scal_use("__cgeonear_f", &f) || !(a > 0) ||
            !(f > 0) || !(f < 1))
            return 198;
        cgeo_ellipsoid_set(&c->e, a, f);
    }
    if (!c->wide && SF_scal_use("__cgeonear_within", &w))
        return 198;
    c->within = w;
    return 0;
}

static ST_retcode phase_read(int base)
{
    cgeo_fail fail;
    cgeo_locs *l = base ? &S.base : &S.nbor;
    size_t nobs = 0;
    int in_order = 0;
    int rc = cgeo_locs_read(l, base ? &S.obs_loc : NULL, &nobs, &in_order, &fail);
    save_fail(fail);
    if (rc)
        return rc;
    if (base)
    {
        S.nobs = nobs;
        S.stage = 1;
        SF_scal_save("__cgeonear_nbase", (double)l->n);
        /* wide form: the ado reorders the data by baseid only when needed */
        SF_scal_save("__cgeonear_inorder", (double)in_order);
    }
    else
    {
        S.stage = 2;
        SF_scal_save("__cgeonear_nnbor", (double)l->n);
    }
    return 0;
}

/* Base location -> neighbor location with the same identifier (ignoreself);
 * both tables are in identifier order */
static uint32_t *self_matches(void)
{
    uint32_t *self = malloc((S.base.n ? S.base.n : 1) * sizeof(*self));
    if (!self)
        return NULL;
    size_t j = 0;
    for (size_t i = 0; i < S.base.n; i++)
    {
        while (j < S.nbor.n && cgeo_cmp_id(&S.nbor, j, &S.base, i) < 0)
            j++;
        self[i] = j < S.nbor.n && cgeo_same_id(&S.nbor, j, &S.base, i) ? (uint32_t)j : CGEO_NONE;
    }
    return self;
}

/* Wide form: every observation of the master data, its base location's
 * neighbors (missing for observations without one) */
static int write_wide(const cgeo_cfg *c, const cgeo_result *res)
{
    const size_t k = c->k;
    if ((size_t)SF_nvars() != 2 * k || (size_t)SF_nobs() != S.nobs)
        return 198;
    int str = S.nbor.id_is_str;
    for (size_t j = 0; j < k; j++)
        if ((SF_var_is_string((int)(2 * j + 1)) != 0) != str || SF_var_is_string((int)(2 * j + 2)))
            return 109;
    ST_IIID vstore = (_stata_)->safestore;
    ST_IIIS sstore = (_stata_)->sstore;
    atomic_int failed = 0;
    const size_t n = S.nobs;
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if (n > 20000)
    for (size_t i = 0; i < n; i++)
    {
        uint32_t g = S.obs_loc[i];
        for (size_t j = 0; j < k; j++)
        {
            int rc;
            ST_int var = (ST_int)(2 * j + 1), obs = (ST_int)(i + 1);
            uint32_t nb = g == CGEO_NONE ? CGEO_NONE : res->loc[(size_t)g * k + j];
            if (nb == CGEO_NONE)
            {
                rc = str ? sstore(var, obs, "") : vstore(var, obs, SV_missval);
                if (!rc)
                    rc = vstore(var + 1, obs, SV_missval);
            }
            else
            {
                rc =
                    str ? sstore(var, obs, S.nbor.id_str[nb]) : vstore(var, obs, S.nbor.id_num[nb]);
                if (!rc)
                    rc = vstore(var + 1, obs, res->val[(size_t)g * k + j]);
            }
            if (rc)
                atomic_store(&failed, 1);
        }
    }
    return atomic_load(&failed) ? 459 : 0;
}

static ST_retcode phase_compute(const char *args)
{
    cgeo_cfg c;
    int timing;
    int rc = read_cfg(args, &c, &timing);
    if (rc)
        return rc;
    if (S.stage != 2)
        return 198;
    int ignore = ctools_parse_bool_option(args, "ignoreself");
    if (S.nbor.n < c.k + (size_t)ignore)
        return 2001;
    /* geonear compares the identifiers: string with numeric is r(109) */
    if (ignore && S.base.id_is_str != S.nbor.id_is_str)
        return 109;
    cgeo_math_init();
    int threads = ctools_get_max_threads();
    double t0 = ctools_timer_seconds();
    cgeo_index index;
    cgeo_index_init(&index);
    cgeo_result res;
    cgeo_result_init(&res);
    uint32_t *self = NULL;
    if (ignore && !(self = self_matches()))
    {
        rc = 920;
        goto done;
    }
    rc = cgeo_index_build(&index, &c, S.nbor.lat, S.nbor.lon, S.nbor.n, threads);
    if (rc)
        goto done;
    double t1 = ctools_timer_seconds();
    rc = cgeo_search_all(&c, &index, S.base.lat, S.base.lon, S.base.n, self, threads, &res);
    if (rc)
        goto done;
    double t2 = ctools_timer_seconds();
    SF_scal_save("__cgeonear_ndist", (double)res.evals);
    SF_scal_save("__cgeonear_threads", (double)threads);
    if (timing)
    {
        SF_scal_save("__cgeonear_time_index", t1 - t0);
        SF_scal_save("__cgeonear_time_search", t2 - t1);
    }
    if (c.wide)
    {
        rc = write_wide(&c, &res);
        if (timing)
            SF_scal_save("__cgeonear_time_store", ctools_timer_seconds() - t2);
        goto done;
    }
    if (res.nout >= (size_t)2147483647)
    {
        rc = 920;
        goto done;
    }
    S.pairs = res;
    cgeo_result_init(&res);
    S.stage = 3;
    SF_scal_save("__cgeonear_nout", (double)S.pairs.nout);
done:
    cgeo_result_free(&res);
    cgeo_index_free(&index);
    free(self);
    return rc;
}

/* Long form: the pairs, into a dataset of exactly that many observations */
static ST_retcode phase_store(void)
{
    if (S.stage != 3)
        return 198;
    if (SF_nvars() != 3 || (size_t)SF_nobs() != S.pairs.nout)
        return 198;
    int bstr = S.base.id_is_str, nstr = S.nbor.id_is_str;
    if ((SF_var_is_string(1) != 0) != bstr || (SF_var_is_string(2) != 0) != nstr ||
        SF_var_is_string(3))
        return 109;
    ST_IIID vstore = (_stata_)->safestore;
    ST_IIIS sstore = (_stata_)->sstore;
    atomic_int failed = 0;
    const size_t nb = S.base.n;
    const cgeo_result *p = &S.pairs;
#pragma omp parallel for schedule(dynamic, 64)                                                     \
    num_threads(ctools_get_max_threads()) if (p->nout > 20000)
    for (size_t b = 0; b < nb; b++)
    {
        for (size_t r = p->row_start[b]; r < p->row_start[b + 1] && !atomic_load(&failed); r++)
        {
            ST_int obs = (ST_int)(r + 1);
            uint32_t nbr = p->nbor[r];
            int rc = bstr ? sstore(1, obs, S.base.id_str[b]) : vstore(1, obs, S.base.id_num[b]);
            if (!rc)
                rc = nstr ? sstore(2, obs, S.nbor.id_str[nbr]) : vstore(2, obs, S.nbor.id_num[nbr]);
            if (!rc)
                rc = vstore(3, obs, p->dist[r]);
            if (rc)
                atomic_store(&failed, 1);
        }
    }
    return atomic_load(&failed) ? 459 : 0;
}

ST_retcode cgeonear_main(const char *args)
{
    char op[16] = {0};
    if (args)
        sscanf(args, "%15s", op);
    const char *rest = args ? args + strlen(op) : "";
    if (!strcmp(op, "base"))
    {
        cgeonear_cleanup_state();
        return phase_read(1);
    }
    if (!strcmp(op, "using"))
        return phase_read(0);
    if (!strcmp(op, "compute"))
        return phase_compute(rest);
    if (!strcmp(op, "store"))
        return phase_store();
    /* geodist's form: no state; the runtime frees any held by the
     * nearest-neighbor form before it runs */
    if (!strcmp(op, "dist"))
        return cgeonear_dist_main(rest);
    if (!strcmp(op, "clear"))
    {
        cgeonear_cleanup_state();
        return 0;
    }
    return 198;
}
