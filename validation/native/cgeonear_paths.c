/* cgeonear's k-d tree search against a comparison of every base with every
 * neighbor, and the distance kernels against geodist's results.
 *
 * The search prunes on straight-line distances and ranks on the exact ones;
 * here every case is also solved by brute force with the same exact
 * distances, and the neighbors, their order (ties by location number) and
 * the reported values must agree exactly: wide and long form, the sphere and
 * the ellipsoid, nearcount(), within() (also 0, which keeps neighbors at the
 * same place), limit(), ignoreself, co-located points, grids with many equal
 * distances, the poles and the date line, trees deep enough to be built in
 * parallel, and one or several threads. */
#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ctools_config.h"
#include "ctools_threads.h"

static int threads = 1;
int ctools_get_max_threads(void)
{
    return threads;
}
/* ctools_sort_order_pairs falls back to one thread without a pool */
ctools_persistent_pool *ctools_get_global_pool(void)
{
    return NULL;
}
int ctools_persistent_pool_submit_batch(ctools_persistent_pool *pool, ctools_thread_func fn,
                                        void *args, size_t count, size_t size)
{
    (void)pool;
    (void)fn;
    (void)args;
    (void)count;
    (void)size;
    return -1;
}
int ctools_persistent_pool_wait(ctools_persistent_pool *pool)
{
    (void)pool;
    return -1;
}

#include "cgeonear/cgeo_math.c"
#include "cgeonear/cgeonear_tree.c"
#include "cgeonear/cgeonear_search.c"
#include "ctools_sort_pairs.c"

static const double MISSING = 8.98846567431158e307;

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static double uniform(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) / 9007199254740992.0;
}
static double normal(void)
{
    double u = uniform(), v = uniform();
    if (u < 1e-300)
        u = 1e-300;
    return sqrt(-2 * log(u)) * cos(2 * CGEO_PI * v);
}

enum
{
    GLOBAL,
    CLUSTER,
    DUPLICATES,
    GRID
};

static void make_points(int kind, size_t n, double *lat, double *lon)
{
    for (size_t i = 0; i < n; i++)
    {
        switch (kind)
        {
        case GLOBAL:
            lat[i] = asin(2 * uniform() - 1) * 180 / CGEO_PI;
            lon[i] = 360 * uniform() - 180;
            break;
        case CLUSTER:
            lat[i] = 40 + 2 * normal();
            lon[i] = -100 + 2 * normal();
            break;
        case DUPLICATES:
            /* few places, each shared by several locations */
            lat[i] = 37 + floor(uniform() * 6) * 0.5;
            lon[i] = -109 + floor(uniform() * 6) * 0.7;
            break;
        default:
            /* integer grid: many neighbors at the same distance */
            lat[i] = floor(uniform() * 9) - 4;
            lon[i] = floor(uniform() * 9) - 4;
            break;
        }
    }
    if (kind == GLOBAL && n > 8)
    {
        lat[0] = 90;
        lat[1] = -90;
        lon[2] = 180;
        lon[3] = -180;
        lat[4] = 0;
        lat[5] = 0;
    }
}

/* What the search must find for base b: every neighbor ranked, then
 * geonear's selection */
static size_t brute(const cgeo_cfg *c, double blat, double blon, const double *nlat,
                    const double *nlon, size_t nn, uint32_t self, cgeo_cand *all)
{
    size_t m = 0;
    for (size_t j = 0; j < nn; j++)
    {
        if ((uint32_t)j == self)
            continue;
        all[m].key = cgeo_rank_distance(c, blat, blon, nlat[j], nlon[j]);
        all[m].loc = (uint32_t)j;
        m++;
    }
    qsort(all, m, sizeof(*all), cand_cmp);
    if (c->wide)
        return m < c->k ? m : c->k;
    size_t keep = c->k < m ? c->k : m, nw = 0;
    while (nw < m && all[nw].key <= c->within)
        nw++;
    if (nw > keep)
        keep = nw;
    if (keep > c->limit)
        keep = c->limit;
    return keep;
}

static int same_value(double a, double b)
{
    return a == b || (isnan(a) && isnan(b));
}

static size_t cases_run = 0;

static void check_case(const cgeo_cfg *c, int bkind, size_t nb, int nkind, size_t nn, int self_set)
{
    double *blat = malloc(nb * sizeof(double)), *blon = malloc(nb * sizeof(double));
    double *nlat = malloc(nn * sizeof(double)), *nlon = malloc(nn * sizeof(double));
    cgeo_cand *all = malloc((nn + 1) * sizeof(*all));
    uint32_t *self = NULL;
    assert(blat && blon && nlat && nlon && all);
    make_points(nkind, nn, nlat, nlon);
    if (self_set)
    {
        /* ignoreself on one set of locations: base b is neighbor b */
        assert(nb <= nn);
        memcpy(blat, nlat, nb * sizeof(double));
        memcpy(blon, nlon, nb * sizeof(double));
        self = malloc(nb * sizeof(*self));
        assert(self);
        for (size_t b = 0; b < nb; b++)
            self[b] = (uint32_t)b;
    }
    else
        make_points(bkind, nb, blat, blon);
    cgeo_index x;
    cgeo_index_init(&x);
    assert(!cgeo_index_build(&x, c, nlat, nlon, nn, threads));
    cgeo_result r[2];
    int saved = threads;
    for (int pass = 0; pass < 2; pass++)
    {
        threads = pass ? 4 : 1;
        cgeo_result_init(&r[pass]);
        assert(!cgeo_search_all(c, &x, blat, blon, nb, self, threads, &r[pass]));
    }
    threads = saved;
    for (size_t b = 0; b < nb; b++)
    {
        size_t m = brute(c, blat[b], blon[b], nlat, nlon, nn, self ? self[b] : CGEO_NONE, all);
        for (int pass = 0; pass < 2; pass++)
        {
            const cgeo_result *res = &r[pass];
            if (c->wide)
            {
                for (size_t j = 0; j < c->k; j++)
                {
                    uint32_t want = j < m ? all[j].loc : CGEO_NONE;
                    double value = j < m ? cgeo_report_value(c, all[j].key) : c->missing;
                    if (res->loc[b * c->k + j] != want ||
                        !same_value(res->val[b * c->k + j], value))
                    {
                        fprintf(stderr, "wide mismatch: base %zu rank %zu got %u want %u\n", b, j,
                                res->loc[b * c->k + j], want);
                        exit(1);
                    }
                }
                continue;
            }
            size_t first = res->row_start[b], count = res->row_start[b + 1] - first;
            if (count != m)
            {
                fprintf(stderr, "long mismatch: base %zu has %zu pairs, want %zu\n", b, count, m);
                exit(1);
            }
            for (size_t j = 0; j < m; j++)
                if (res->nbor[first + j] != all[j].loc ||
                    !same_value(res->dist[first + j], cgeo_report_value(c, all[j].key)))
                {
                    fprintf(stderr, "long mismatch: base %zu pair %zu\n", b, j);
                    exit(1);
                }
        }
    }
    assert(r[0].evals == r[1].evals);
    for (int pass = 0; pass < 2; pass++)
        cgeo_result_free(&r[pass]);
    cgeo_index_free(&x);
    free(self);
    free(all);
    free(blat);
    free(blon);
    free(nlat);
    free(nlon);
    cases_run++;
}

static cgeo_cfg config(int wide, int sphere, size_t k, double within, size_t limit, int miles)
{
    cgeo_cfg c;
    memset(&c, 0, sizeof(c));
    c.wide = wide;
    c.sphere = sphere;
    c.miles = miles;
    c.k = k;
    c.limit = limit;
    c.r = 6371;
    cgeo_ellipsoid_set(&c.e, 6378137, 1 / 298.257223563);
    c.within = within;
    c.missing = MISSING;
    return c;
}

/* geodist's results in StataNow 19.5 for Windows (r(distance)) */
static void check_kernels(void)
{
    const double d2r = CGEO_PI / 180;
    cgeo_ellipsoid e;
    cgeo_ellipsoid_set(&e, 6378137, 1 / 298.257223563);
    double su1, cu1, su2, cu2, ba, smd;
    int passes;
    cgeo_reduced(&e, d2r * 42.265837, &su1, &cu1);
    cgeo_reduced(&e, d2r * 42.207667, &su2, &cu2);
    int st = cgeo_vincenty(&e, su1, cu1, su2, cu2, d2r * (-83.356022 - -83.748696), 24, &ba, &smd,
                           &passes);
    assert(st == CGEO_VINCENTY_CONVERGED && passes == 4);
    double miles = ba * (((smd) / 1000) / 1.609344);
    assert(fabs(miles / 20.53644088672285 - 1) < 1e-13);
    double sphere =
        cgeo_haversine_deg(42.265837, -83.748696, 42.207667, -83.356022, d2r, 6371 / 1.609344);
    assert(fabs(sphere / 20.48530988426987 - 1) < 1e-13);
    /* coincident points; antipodal points do not converge */
    cgeo_reduced(&e, d2r * 10, &su1, &cu1);
    assert(cgeo_vincenty(&e, su1, cu1, su1, cu1, 0, 24, &ba, &smd, &passes) ==
           CGEO_VINCENTY_COINCIDENT);
    cgeo_reduced(&e, d2r * -10, &su2, &cu2);
    assert(cgeo_vincenty(&e, su1, cu1, su2, cu2, d2r * (-160 - 20), 24, &ba, &smd, &passes) ==
           CGEO_VINCENTY_FAILED);
    /* Stata's asin and tan */
    assert(cgeo_asin(1) == atan(INFINITY));
    assert(cgeo_tan(0.5) == cgeo_sin(0.5) / cgeo_cos(0.5));
    assert(cgeo_haversine_rad(0.3, cgeo_cos(0.3), 1.0, 0.3, cgeo_cos(0.3), 1.0) == 0);
}

int main(void)
{
    cgeo_math_init();
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    check_kernels();
    const size_t none = SIZE_MAX;
    /* wide form */
    for (int sphere = 1; sphere >= 0; sphere--)
    {
        size_t ks[] = {1, 3, 7};
        for (int i = 0; i < 3; i++)
        {
            cgeo_cfg c = config(1, sphere, ks[i], 0, none, i == 1);
            check_case(&c, CLUSTER, 400, CLUSTER, 900, 0);
            check_case(&c, GLOBAL, 300, GLOBAL, 1200, 0);
            check_case(&c, DUPLICATES, 200, DUPLICATES, 150, 0);
            check_case(&c, GRID, 200, GRID, 120, 0);
            check_case(&c, CLUSTER, 100, CLUSTER, 100, 1);
            check_case(&c, GRID, 60, GRID, 60, 1);
            check_case(&c, GLOBAL, ks[i] + 1, GLOBAL, ks[i] + 1, 1);
        }
        /* a tree built by parallel tasks */
        cgeo_cfg c = config(1, sphere, 2, 0, none, 0);
        check_case(&c, CLUSTER, sphere ? 40 : 12, CLUSTER, 140000, 0);
        check_case(&c, GLOBAL, sphere ? 40 : 12, GLOBAL, 140000, 0);
    }
    /* long form: nearcount(), within() (0 keeps co-located neighbors), limit() */
    for (int sphere = 1; sphere >= 0; sphere--)
    {
        struct
        {
            size_t k;
            double within;
            size_t limit;
        } runs[] = {{1, 0, none}, {3, 0, none}, {0, 80, none}, {2, 80, none},
                    {2, 80, 4},   {0, 5000, 3}, {4, 0, 4},     {0, 0, none}};
        for (size_t i = 0; i < sizeof(runs) / sizeof(runs[0]); i++)
        {
            cgeo_cfg c = config(0, sphere, runs[i].k, runs[i].within, runs[i].limit, i == 3);
            check_case(&c, CLUSTER, 300, CLUSTER, 700, 0);
            check_case(&c, DUPLICATES, 150, DUPLICATES, 150, 0);
            check_case(&c, GRID, 100, GRID, 100, 1);
            check_case(&c, GLOBAL, 100, GLOBAL, 600, 0);
        }
        cgeo_cfg c = config(0, sphere, 1, 25, none, 0);
        check_case(&c, CLUSTER, sphere ? 30 : 10, CLUSTER, 100000, 0);
    }
    /* every pair: within() beyond any distance */
    cgeo_cfg all = config(0, 1, 0, 25000, none, 0);
    check_case(&all, GLOBAL, 40, GLOBAL, 60, 0);
    /* the radius is part of the long-form ranking distance */
    cgeo_cfg radius = config(0, 1, 2, 30, none, 0);
    radius.r = 1000;
    check_case(&radius, CLUSTER, 200, CLUSTER, 300, 0);
    printf("cgeonear search: %zu cases agree with brute force\n", cases_run);
    return 0;
}
