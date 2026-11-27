/* Geodetic distances computed as geonear and geodist compute them in Stata.
 *
 * Each kernel repeats the arithmetic of the Stata expression it replaces in
 * Stata's order of evaluation, so that the results agree to the last bit:
 *   - a*b/c is a*(b/c): Stata gives division precedence over multiplication
 *     (a*b*c/d is (a*b)*(c/d), a*b/c/d is a*((b/c)/d));
 *   - x^2 is x*x;
 *   - tan(x) is sin(x)/cos(x) and asin(y) is atan(y/sqrt(1-y*y));
 *   - on Windows, sin, cos and atan2 are those of the Universal C Runtime
 *     (ucrtbase.dll), which Stata calls; MinGW's own differ in the last bit
 *     for a few percent of arguments. atan is the C library's (MinGW's
 *     matches Stata's). cgeo_math_init() resolves them; until it has run,
 *     and where ucrtbase.dll is missing, the C library's functions are used.
 * Measured on StataNow 19.5 (Windows) over 100,000 random arguments: sin,
 * cos, tan, atan, asin and sqrt identical in every case, atan2 in all but 7.
 *
 * Missing values: a failed Vincenty iteration is reported through the
 * status, never as a NaN the caller must recognise. */
#ifndef CGEO_MATH_H
#define CGEO_MATH_H

#include <math.h>

typedef double (*cgeo_fn1)(double);
typedef double (*cgeo_fn2)(double, double);

extern cgeo_fn1 cgeo_sin_fn;
extern cgeo_fn1 cgeo_cos_fn;
extern cgeo_fn2 cgeo_atan2_fn;

/* Resolve the runtime's functions. Call on the Stata thread before any
 * parallel region; idempotent. Returns 1 when Stata's runtime was found. */
int cgeo_math_init(void);

static inline double cgeo_sin(double x)
{
    return cgeo_sin_fn(x);
}
static inline double cgeo_cos(double x)
{
    return cgeo_cos_fn(x);
}
static inline double cgeo_atan2(double y, double x)
{
    return cgeo_atan2_fn(y, x);
}
static inline double cgeo_tan(double x)
{
    return cgeo_sin(x) / cgeo_cos(x);
}
/* Stata's asin(); y = 1 gives atan(+inf) = pi/2 */
static inline double cgeo_asin(double y)
{
    return atan(y / sqrt(1.0 - y * y));
}

/* Stata's min(1, x) for the nonmissing x of the haversine formula */
static inline double cgeo_min1(double x)
{
    return x > 1.0 ? 1.0 : x;
}

/* geonear's great-circle angle (radians) between two points in radians:
 *   2 * asin(min(1, sqrt(sin((lat2 - lat1) / 2)^2 +
 *       cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2)))
 * with cos1 = cos(lat1), cos2 = cos(lat2). Point 1 is the base location. */
static inline double cgeo_haversine_rad(double lat1, double cos1, double lon1, double lat2,
                                        double cos2, double lon2)
{
    double s1 = cgeo_sin((lat2 - lat1) / 2);
    double s2 = cgeo_sin((lon2 - lon1) / 2);
    double h = s1 * s1 + (cos1 * cos2) * (s2 * s2);
    return 2 * cgeo_asin(cgeo_min1(sqrt(h)));
}

/* geodist's great-circle formula on coordinates in degrees:
 *   2 * asin(min(1, sqrt(sin((lat2 - lat1) * d2r / 2)^2 +
 *       cos(lat1 * d2r) * cos(lat2 * d2r) * sin((lon2 - lon1) * d2r / 2)^2)))
 *   * radius [/ 1.609344]
 * factor is radius, or radius / 1.609344 for miles (Stata's precedence). */
static inline double cgeo_haversine_deg(double lat1, double lon1, double lat2, double lon2,
                                        double d2r, double factor)
{
    double s1 = cgeo_sin((lat2 - lat1) * (d2r / 2));
    double s2 = cgeo_sin((lon2 - lon1) * (d2r / 2));
    double h = s1 * s1 + (cgeo_cos(lat1 * d2r) * cgeo_cos(lat2 * d2r)) * (s2 * s2);
    return (2 * cgeo_asin(cgeo_min1(sqrt(h)))) * factor;
}

/* A reference ellipsoid as geonear and geodist hold it: a (meters), f
 * (flattening, 1 / the reciprocal the user gives), b = a - a * f and the
 * constant (a^2 - b^2) / (b^2) of u^2. */
typedef struct
{
    double a, f, b, one_minus_f, usq;
} cgeo_ellipsoid;

static inline void cgeo_ellipsoid_set(cgeo_ellipsoid *e, double a, double f)
{
    e->a = a;
    e->f = f;
    e->b = a - a * f;
    e->one_minus_f = 1 - f;
    e->usq = (a * a - e->b * e->b) / (e->b * e->b);
}

/* Reduced latitude of a latitude in radians: U = atan((1 - f) * tan(lat)) */
static inline void cgeo_reduced(const cgeo_ellipsoid *e, double lat, double *sin_u, double *cos_u)
{
    double u = atan(e->one_minus_f * cgeo_tan(lat));
    *sin_u = cgeo_sin(u);
    *cos_u = cgeo_cos(u);
}

enum
{
    CGEO_VINCENTY_CONVERGED = 0,
    CGEO_VINCENTY_COINCIDENT = 1, /* sin(sigma) == 0: the distance is 0 */
    CGEO_VINCENTY_FAILED = 2      /* no convergence within the iterations */
};

/* Vincenty's (1975) inverse solution as geonear_ellipsoid and geodist run it:
 * at most `bodies` passes of the lambda iteration (24 for geonear and for
 * geodist's default maxiter(25)), converged when |lambda - lambda_old| <=
 * 1e-12. On convergence *ba = b * A and *smd = sigma - delta_sigma, so that
 * geonear's distance in km is ba * (smd / 1000). *passes counts the passes
 * run (the pass that converged, or found coincident points, included). */
static inline int cgeo_vincenty(const cgeo_ellipsoid *e, double sin_u1, double cos_u1,
                                double sin_u2, double cos_u2, double l, int bodies, double *ba,
                                double *smd, int *passes)
{
    const double f = e->f;
    const double cu1cu2 = cos_u1 * cos_u2, su1su2 = sin_u1 * sin_u2;
    const double cu1su2 = cos_u1 * sin_u2, su1cu2 = sin_u1 * cos_u2;
    const double two_su1 = 2 * sin_u1;
    double lambda = l;
    double sin_sigma = 0, cos_sigma = 0, sigma = 0, cos_sq_alpha = 0, c2m = 0;
    int converged = 0, pass = 0;
    while (pass < bodies)
    {
        pass++;
        double sl = cgeo_sin(lambda), cl = cgeo_cos(lambda);
        double p = cos_u2 * sl;
        double q = cu1su2 - su1cu2 * cl;
        sin_sigma = sqrt(p * p + q * q);
        if (sin_sigma == 0)
        {
            *passes = pass;
            *ba = 0;
            *smd = 0;
            return CGEO_VINCENTY_COINCIDENT;
        }
        cos_sigma = su1su2 + cu1cu2 * cl;
        sigma = cgeo_atan2(sin_sigma, cos_sigma);
        double sin_alpha = cu1cu2 * (sl / sin_sigma);
        cos_sq_alpha = 1 - sin_alpha * sin_alpha;
        /* both points on the equator: Stata's division by zero is missing,
         * which the next statement replaces by 0 */
        c2m = cos_sq_alpha == 0 ? 0 : cos_sigma - two_su1 * (sin_u2 / cos_sq_alpha);
        double lambda_old = lambda;
        double c = ((f / 16) * cos_sq_alpha) * (4 + f * (4 - 3 * cos_sq_alpha));
        lambda =
            l + (((1 - c) * f) * sin_alpha) *
                    (sigma + (c * sin_sigma) * (c2m + (c * cos_sigma) * (-1 + 2 * (c2m * c2m))));
        if (fabs(lambda - lambda_old) <= 1e-12)
        {
            converged = 1;
            break;
        }
    }
    *passes = pass;
    if (!converged)
        return CGEO_VINCENTY_FAILED;
    double u_sq = cos_sq_alpha * e->usq;
    double big_a = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));
    double big_b = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)));
    double delta_sigma =
        (big_b * sin_sigma) *
        (c2m + (big_b / 4) * (cos_sigma * (-1 + 2 * (c2m * c2m)) -
                              (((big_b / 6) * c2m) * (-3 + 4 * (sin_sigma * sin_sigma))) *
                                  (-3 + 4 * (c2m * c2m))));
    *ba = e->b * big_a;
    *smd = sigma - delta_sigma;
    return CGEO_VINCENTY_CONVERGED;
}

#endif
