/* cgeonear's distance form: geodist's distances, one per observation or for
 * one pair. Reached as the plugin command "cgeonear dist ...":
 *
 * "cgeonear dist vars kinds=VVCC [sphere] [miles] bodies=# check=# [out]
 *  [marker] [allmissing]"
 *   Plugin variables: the coordinates that are variables, in the order
 *   lat1 lon1 lat2 lon2; then the new double variable (out); then the
 *   if/in marker (marker). kinds says which of the four are variables (V)
 *   and which constants (C, in scalars __cgeonear_dc1 .. __cgeonear_dc4).
 *   geodist marks out the sample one coordinate at a time and checks the
 *   range of each variable on the observations still in it; the first
 *   coordinate (among the first `check`) whose check fails is returned in
 *   local cgeonear_dviol. With out, and no failed check, every observation
 *   of the new variable is written: the distance in the sample, missing
 *   elsewhere, .a where Vincenty's iteration did not converge (counted in
 *   local cgeonear_dnfail). __cgeonear_dnsample is the size of the sample.
 * "cgeonear dist scalar [sphere] [miles] bodies=#"
 *   One distance from the four constants: __cgeonear_dd, and geodist's
 *   iteration count in __cgeonear_diter (local cgeonear_dfailed is 1 when
 *   geodist would report a failure to converge).
 * Model: sphere uses __cgeonear_dfac (radius, or radius / 1.609344 for
 * miles); otherwise __cgeonear_da and __cgeonear_df (semi-major axis in
 * meters, flattening). bodies is maxiter() - 1, geodist's number of
 * passes. */
#include <math.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cgeonear_internal.h"
#include "cgeo_math.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include "ctools_types.h"

#define CGEO_PI 3.14159265358979323846

typedef struct
{
    int sphere, miles, bodies;
    double fac; /* sphere */
    cgeo_ellipsoid e;
    double d2r;
} geodist_cfg;

/* Stata's .a: the next missing code after . */
static double missing_a(void)
{
    return SV_missval + ldexp(1.0, 1011);
}

static int read_model(const char *args, geodist_cfg *c)
{
    memset(c, 0, sizeof(*c));
    c->sphere = ctools_parse_bool_option(args, "sphere");
    c->miles = ctools_parse_bool_option(args, "miles");
    if (ctools_parse_int_option(args, "bodies", &c->bodies) != CTOOLS_PARSE_OK)
        return 198;
    c->d2r = CGEO_PI / 180;
    if (c->sphere)
        return SF_scal_use("__cgeonear_dfac", &c->fac) ? 198 : 0;
    double a, f;
    if (SF_scal_use("__cgeonear_da", &a) || SF_scal_use("__cgeonear_df", &f))
        return 198;
    cgeo_ellipsoid_set(&c->e, a, f);
    return 0;
}

/* geodist's ellipsoid distance from precomputed reduced latitudes */
static double vincenty_value(const geodist_cfg *c, double su1, double cu1, double su2, double cu2,
                             double lon1, double lon2, int *status, int *passes)
{
    double ba, smd;
    *status = cgeo_vincenty(&c->e, su1, cu1, su2, cu2, c->d2r * (lon2 - lon1), c->bodies, &ba, &smd,
                            passes);
    if (*status == CGEO_VINCENTY_COINCIDENT)
        return 0;
    if (*status == CGEO_VINCENTY_FAILED)
        return missing_a();
    return c->miles ? ba * (((smd) / 1000) / 1.609344) : ba * ((smd) / 1000);
}

static void reduced_deg(const geodist_cfg *c, double lat, double *su, double *cu)
{
    cgeo_reduced(&c->e, c->d2r * lat, su, cu);
}

static ST_retcode geodist_scalar(const char *args)
{
    geodist_cfg c;
    int rc = read_model(args, &c);
    if (rc)
        return rc;
    double v[4];
    char name[32];
    for (int j = 0; j < 4; j++)
    {
        snprintf(name, sizeof(name), "__cgeonear_dc%d", j + 1);
        if (SF_scal_use(name, &v[j]))
            return 198;
    }
    cgeo_math_init();
    double d;
    int iter = 0, failed = 0;
    if (c.sphere)
        d = cgeo_haversine_deg(v[0], v[1], v[2], v[3], c.d2r, c.fac);
    else
    {
        /* geodist: while ++iter < maxiter & more { ... }; failure when the
         * loop ends with iter == maxiter, also after converging on its
         * last pass */
        if (c.bodies < 1)
            return 111;
        double su1, cu1, su2, cu2;
        int status, passes;
        reduced_deg(&c, v[0], &su1, &cu1);
        reduced_deg(&c, v[2], &su2, &cu2);
        d = vincenty_value(&c, su1, cu1, su2, cu2, v[1], v[3], &status, &passes);
        /* coincident points leave the loop with continue, break; otherwise
         * the loop's test has incremented iter once more */
        iter = status == CGEO_VINCENTY_COINCIDENT ? passes : passes + 1;
        if (iter == c.bodies + 1)
        {
            failed = 1;
            d = missing_a();
        }
    }
    SF_scal_save("__cgeonear_dd", d);
    SF_scal_save("__cgeonear_diter", (double)iter);
    SF_macro_save("_cgeonear_dfailed", failed ? "1" : "0");
    return 0;
}

static ST_retcode geodist_vars(const char *args)
{
    geodist_cfg c;
    int rc = read_model(args, &c);
    if (rc)
        return rc;
    char kinds[8] = {0};
    int check = 4;
    if (ctools_parse_string_option(args, "kinds", kinds, sizeof(kinds)) != CTOOLS_PARSE_OK ||
        strlen(kinds) != 4 || ctools_parse_int_option(args, "check", &check) != CTOOLS_PARSE_OK)
        return 198;
    int out = ctools_parse_bool_option(args, "out");
    int marker = ctools_parse_bool_option(args, "marker");
    int allmissing = ctools_parse_bool_option(args, "allmissing");
    int nv = 0, col[4];
    double konst[4] = {0, 0, 0, 0};
    char name[32];
    for (int j = 0; j < 4; j++)
    {
        if (kinds[j] == 'V')
            col[j] = nv++;
        else if (kinds[j] == 'C')
        {
            col[j] = -1;
            snprintf(name, sizeof(name), "__cgeonear_dc%d", j + 1);
            if (SF_scal_use(name, &konst[j]))
                return 198;
        }
        else
            return 198;
    }
    if (!nv || SF_nvars() != nv + out + marker)
        return 198;
    int outvar = nv + 1;
    if (out && SF_var_is_string(outvar))
        return 109;
    cgeo_math_init();
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    double *res = NULL;
    int vars[5];
    for (int j = 0; j < nv; j++)
        vars[j] = j + 1;
    if (marker)
        vars[nv] = SF_nvars();
    int flags = CTOOLS_LOAD_NO_SORT_ORDER | (marker ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_SKIP_IF);
    rc = ctools_stata_rc(ctools_data_load(&fd, vars, (size_t)nv, 0, 0, flags));
    if (rc)
        goto done;
    for (int j = 0; j < nv; j++)
        if (fd.data.vars[j].type != STATA_TYPE_DOUBLE)
        {
            rc = 109;
            goto done;
        }
    const size_t n = fd.data.nobs;
    const double *colv[4] = {NULL, NULL, NULL, NULL};
    for (int j = 0; j < 4; j++)
        if (col[j] >= 0)
            colv[j] = fd.data.vars[col[j]].data.dbl;
    /* geodist's checks: coordinate j's range on the observations whose
     * variables up to j are all nonmissing */
    atomic_int first_bad = 5;
    size_t nsample = 0;
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads())                    \
    reduction(+ : nsample) if (n > 50000)
    for (size_t i = 0; i < n; i++)
    {
        int in = 1;
        for (int j = 0; j < 4 && in; j++)
        {
            if (!colv[j])
                continue;
            double x = colv[j][i];
            if (x >= SV_missval)
            {
                in = 0;
                break;
            }
            double lim = (j & 1) ? 180 : 90;
            if (j < check && (x > lim || x < -lim))
            {
                int cur = atomic_load(&first_bad);
                while (j + 1 < cur && !atomic_compare_exchange_weak(&first_bad, &cur, j + 1))
                    ;
            }
        }
        nsample += (size_t)in;
    }
    int viol = atomic_load(&first_bad);
    char buf[32];
    snprintf(buf, sizeof(buf), "%d", viol == 5 ? 0 : viol);
    SF_macro_save("_cgeonear_dviol", buf);
    SF_scal_save("__cgeonear_dnsample", (double)nsample);
    if (viol != 5 || check < 4 || !out)
        goto done;
    if (!nsample)
        allmissing = 1;
    size_t nobs = (size_t)SF_nobs();
    res = malloc((nobs ? nobs : 1) * sizeof(*res));
    if (!res)
    {
        rc = 920;
        goto done;
    }
    for (size_t i = 0; i < nobs; i++)
        res[i] = SV_missval;
    size_t nfail = 0;
    if (!allmissing)
    {
        /* reduced latitudes of constant coordinates, computed once */
        double su_k[2] = {0, 0}, cu_k[2] = {0, 0};
        if (!c.sphere)
            for (int s = 0; s < 2; s++)
                if (!colv[2 * s])
                    reduced_deg(&c, konst[2 * s], &su_k[s], &cu_k[s]);
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads())                    \
    reduction(+ : nfail) if (n > 20000)
        for (size_t i = 0; i < n; i++)
        {
            double v[4];
            int in = 1;
            for (int j = 0; j < 4; j++)
            {
                v[j] = colv[j] ? colv[j][i] : konst[j];
                if (colv[j] && v[j] >= SV_missval)
                    in = 0;
            }
            if (!in)
                continue;
            double d;
            if (c.sphere)
                d = cgeo_haversine_deg(v[0], v[1], v[2], v[3], c.d2r, c.fac);
            else
            {
                double su1 = su_k[0], cu1 = cu_k[0], su2 = su_k[1], cu2 = cu_k[1];
                if (colv[0])
                    reduced_deg(&c, v[0], &su1, &cu1);
                if (colv[2])
                    reduced_deg(&c, v[2], &su2, &cu2);
                int status, passes;
                d = vincenty_value(&c, su1, cu1, su2, cu2, v[1], v[3], &status, &passes);
                nfail += status == CGEO_VINCENTY_FAILED;
            }
            res[fd.obs_map[i] - 1] = d;
        }
    }
    snprintf(buf, sizeof(buf), "%zu", nfail);
    SF_macro_save("_cgeonear_dnfail", buf);
    ST_IIID vstore = (_stata_)->safestore;
    atomic_int failed = 0;
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if (nobs > 20000)
    for (size_t i = 0; i < nobs; i++)
        if (vstore(outvar, (ST_int)(i + 1), res[i]))
            atomic_store(&failed, 1);
    if (atomic_load(&failed))
        rc = 459;
done:
    free(res);
    ctools_filtered_data_free(&fd);
    return rc;
}

ST_retcode cgeonear_dist_main(const char *args)
{
    char op[16] = {0};
    if (!args)
        args = "";
    while (*args == ' ')
        args++;
    sscanf(args, "%15s", op);
    const char *rest = args + strlen(op);
    if (!strcmp(op, "vars"))
        return geodist_vars(rest);
    if (!strcmp(op, "scalar"))
        return geodist_scalar(rest);
    return 198;
}
