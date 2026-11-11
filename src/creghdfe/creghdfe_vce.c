/*
 * creghdfe_vce.c
 *
 * Variance-Covariance Estimation: unadjusted, robust (HC1), clustered
 * Implements reghdfe-compatible VCE formulas using shared ctools_vce engine.
 * Part of the ctools Stata plugin suite
 */

#include <string.h>
#include "creghdfe_vce.h"
#include "../ctools_matrix.h"

/* A model without residual degrees of freedom (a saturated model) has no
 * VCE; reghdfe then posts zeros (its V is missing and set to zero). */
static ST_retcode zero_vce(ST_int K, ST_double *V)
{
    memset(V, 0, (size_t)K * K * sizeof(ST_double));
    return 0;
}

/* The sandwich D*M*D is symmetric only up to rounding; as reghdfe's
 * _makesymmetric, the lower triangle is copied into the upper one. */
static void make_symmetric(ST_int K, ST_double *V)
{
    for (ST_int i = 0; i < K; i++)
        for (ST_int j = i + 1; j < K; j++)
            V[(size_t)i * K + j] = V[(size_t)j * K + i];
}

/* ========================================================================
 * Compute unadjusted VCE: V = sigma^2 * (X'X)^(-1)
 * where sigma^2 = RSS / df_r
 * ======================================================================== */

ST_retcode compute_vce_unadjusted(
    const ST_double *inv_xx,  /* (X'X)^(-1), K x K */
    ST_double rss,
    ST_double df_r,
    ST_int K,
    ST_double *V              /* Output: K x K VCE matrix */
)
{
    ST_int i;
    if (!isfinite(df_r)) return 498;
    if (df_r <= 0.0) return zero_vce(K, V);
    ST_double sigma2 = rss / df_r;

    for (i = 0; i < K * K; i++) {
        V[i] = sigma2 * inv_xx[i];
    }
    return 0;
}

/* ========================================================================
 * Compute robust VCE (HC1): V = D * M * D * dof_adj
 * where D = (X'X)^(-1), M = X'WX with W = diag(e^2), and X includes constant
 * This matches reghdfe.mata lines 3929-3930
 * ======================================================================== */

ST_retcode compute_vce_robust(
    const ST_double *data,    /* N x (K_keep+2) matrix: y, X1...X_K_keep, constant */
    const ST_double *resid,   /* N x 1 pre-computed residuals from partialled X */
    const ST_double *inv_xx,  /* K_with_cons x K_with_cons inverse (X vars + constant) */
    const ST_double *weights, /* N x 1 weights (NULL for unweighted) */
    ST_int weight_type,       /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_int N,
    ST_double N_eff,             /* Effective N (= sum(weights) for fweight, else N) */
    ST_int K_with_cons,       /* Number of X vars + constant (excludes y) */
    ST_double df_a,              /* Degrees of freedom absorbed by FEs */
    ST_double *V              /* Output: K_with_cons x K_with_cons */
)
{
    /* df_m is K_with_cons - 1 (X vars only, not constant) for consistency with reghdfe */
    ST_int df_m = K_with_cons - 1;
    /* HC1 adjustment: N_eff / (N_eff - df_m - df_a) - matches reghdfe.mata line 3924 */
    if (!isfinite(N_eff)) return 498;
    if (N_eff <= df_m + df_a) return zero_vce(K_with_cons, V);
    ST_double dof_adj = N_eff / (N_eff - df_m - df_a);

    /* X_eff = data + N (skip y column, point to X1...X_K_keep, constant) */
    ctools_vce_data d;
    d.X_eff = data + N;
    d.D = inv_xx;
    d.resid = resid;
    d.weights = weights;
    d.weight_type = weight_type;
    d.N = N;
    d.K = K_with_cons;
    d.normalize_weights = 1;

    ST_retcode rc = ctools_vce_robust(&d, dof_adj, V);
    if (!rc) make_symmetric(K_with_cons, V);
    return rc;
}

/* ========================================================================
 * Compute clustered VCE using reghdfe's formula:
 * V = D * M * D * dof_adj
 * where D = (X'X)^(-1), M = sum_c (X_c'e_c)(e_c'X_c), X includes constant
 * dof_adj = (N-1)/(N - nested_adj - df_m - df_a) * M/(M-1)
 * where nested_adj = 1 if df_a_nested > 0, else 0
 * This matches reghdfe.mata lines 4021, 4026, 4031
 * ======================================================================== */

ST_retcode compute_vce_cluster(
    const ST_double *data,      /* N x (K_keep+2) matrix: y, X1...X_K_keep, constant */
    const ST_double *resid,     /* N x 1 pre-computed residuals from partialled X */
    const ST_double *inv_xx,    /* K_with_cons x K_with_cons inverse (X vars + constant) */
    const ST_double *weights,   /* N x 1 weights (NULL for unweighted) */
    ST_int weight_type,         /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    const ST_int *cluster_ids,  /* N x 1 cluster IDs (0-indexed) */
    ST_int N,
    ST_double N_eff,               /* Effective N (= sum(weights) for fweight, else N) */
    ST_int K_with_cons,         /* Number of X vars + constant (excludes y) */
    ST_int num_clusters,
    ST_double *V,               /* Output: K_with_cons x K_with_cons */
    ST_int df_m,                /* Degrees of freedom model (number of X vars, excluding constant) */
    ST_double df_a,                /* Degrees of freedom absorbed by FEs */
    ST_double df_a_nested          /* Degrees of freedom nested within cluster (for adjustment) */
)
{
    /* Handle degenerate case: single cluster (VCE undefined, like reghdfe) */
    if (num_clusters <= 1) {
        return 498;
    }

    /* reghdfe's DOF adjustment: (N-1)/(N - nested_adj - df_m - S.df_a) * M/(M-1)
     * Use N_eff for fweight (sum of weights) */
    ST_double df_a_adjusted = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;
    if (!isfinite(N_eff)) return 498;
    if (N_eff <= nested_adj + df_m + df_a_adjusted) return zero_vce(K_with_cons, V);
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (N_eff - nested_adj - df_m - df_a_adjusted)) * ((ST_double)num_clusters / (num_clusters - 1));

    /* X_eff = data + N (skip y column, point to X1...X_K_keep, constant) */
    ctools_vce_data d;
    d.X_eff = data + N;
    d.D = inv_xx;
    d.resid = resid;
    d.weights = weights;
    d.weight_type = weight_type;
    d.N = N;
    d.K = K_with_cons;
    d.normalize_weights = 1;

    ST_retcode rc = ctools_vce_cluster(&d, cluster_ids, num_clusters, dof_adj, V);
    if (!rc) make_symmetric(K_with_cons, V);
    return rc;
}

/* ========================================================================
 * Multi-way clustering and Driscoll-Kraay (reghdfe_vce_cluster with more
 * than one cluster variable, reghdfe_vce_dkraay)
 * ======================================================================== */

/* Eigen decomposition of a symmetric n x n matrix by cyclic Jacobi
 * rotations: lambda (n) and the eigenvectors in the columns of U (n x n,
 * U[i*n+j] = component i of vector j). Returns 0 or 920. */
static ST_retcode sym_eigen(const ST_double *A, ST_int n, ST_double *lambda, ST_double *U)
{
    ST_double *a = (ST_double *)malloc((size_t)n * n * sizeof(*a));
    if (!a) return 920;
    memcpy(a, A, (size_t)n * n * sizeof(*a));
    for (ST_int i = 0; i < n; i++)
        for (ST_int j = 0; j < n; j++) U[(size_t)i * n + j] = i == j;
    ST_double total = 0.0;
    for (size_t i = 0; i < (size_t)n * n; i++) total += a[i] * a[i];
    for (int sweep = 0; sweep < 100; sweep++) {
        ST_double off = 0.0;
        for (ST_int p = 0; p < n; p++)
            for (ST_int q = p + 1; q < n; q++) off += a[(size_t)p * n + q] * a[(size_t)p * n + q];
        if (!(off > 1e-32 * total)) break;
        for (ST_int p = 0; p < n; p++) {
            for (ST_int q = p + 1; q < n; q++) {
                ST_double apq = a[(size_t)p * n + q];
                if (apq == 0.0) continue;
                ST_double theta = (a[(size_t)q * n + q] - a[(size_t)p * n + p]) / (2.0 * apq);
                ST_double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                ST_double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (ST_int k = 0; k < n; k++) {
                    ST_double akp = a[(size_t)k * n + p], akq = a[(size_t)k * n + q];
                    a[(size_t)k * n + p] = c * akp - s * akq;
                    a[(size_t)k * n + q] = s * akp + c * akq;
                }
                for (ST_int k = 0; k < n; k++) {
                    ST_double apk = a[(size_t)p * n + k], aqk = a[(size_t)q * n + k];
                    a[(size_t)p * n + k] = c * apk - s * aqk;
                    a[(size_t)q * n + k] = s * apk + c * aqk;
                }
                for (ST_int k = 0; k < n; k++) {
                    ST_double ukp = U[(size_t)k * n + p], ukq = U[(size_t)k * n + q];
                    U[(size_t)k * n + p] = c * ukp - s * ukq;
                    U[(size_t)k * n + q] = s * ukp + c * ukq;
                }
            }
        }
    }
    for (ST_int i = 0; i < n; i++) lambda[i] = a[(size_t)i * n + i];
    free(a);
    return 0;
}

/* Negative eigenvalues of the n x n matrix V set to zero (V = U diag(max(l,
 * 0)) U'); *fixed tells whether any was negative */
static ST_retcode psd_clip(ST_double *V, ST_int n, ST_int *fixed)
{
    *fixed = 0;
    if (n < 1) return 0;
    ST_double *lambda = (ST_double *)malloc((size_t)n * sizeof(*lambda));
    ST_double *U = (ST_double *)malloc((size_t)n * n * sizeof(*U));
    ST_retcode rc = (lambda && U) ? sym_eigen(V, n, lambda, U) : 920;
    if (!rc) {
        for (ST_int i = 0; i < n; i++) if (lambda[i] < 0.0) *fixed = 1;
        if (*fixed) {
            for (ST_int i = 0; i < n; i++) {
                for (ST_int j = 0; j < n; j++) {
                    ST_double v = 0.0;
                    for (ST_int k = 0; k < n; k++)
                        if (lambda[k] > 0.0) v += U[(size_t)i * n + k] * lambda[k] * U[(size_t)j * n + k];
                    V[(size_t)i * n + j] = v;
                }
            }
        }
    }
    free(lambda); free(U);
    return rc;
}

/* The smallest eigenvalue of V (n square) scaled to unit diagonal (rows with
 * a zero diagonal unscaled), in *min_eig */
static ST_retcode scaled_min_eigen(const ST_double *V, ST_int n, ST_double *min_eig)
{
    ST_double *S = (ST_double *)malloc((size_t)n * n * sizeof(*S));
    ST_double *s = (ST_double *)malloc((size_t)n * sizeof(*s));
    ST_double *lambda = (ST_double *)malloc((size_t)n * sizeof(*lambda));
    ST_double *U = (ST_double *)malloc((size_t)n * n * sizeof(*U));
    ST_retcode rc = (S && s && lambda && U) ? 0 : 920;
    if (!rc) {
        for (ST_int i = 0; i < n; i++) {
            ST_double d = V[(size_t)i * n + i];
            s[i] = d > 0.0 ? sqrt(d) : 1.0;
        }
        for (ST_int i = 0; i < n; i++)
            for (ST_int j = 0; j < n; j++) S[(size_t)i * n + j] = V[(size_t)i * n + j] / (s[i] * s[j]);
        rc = sym_eigen(S, n, lambda, U);
    }
    if (!rc) {
        *min_eig = lambda[0];
        for (ST_int i = 1; i < n; i++) if (lambda[i] < *min_eig) *min_eig = lambda[i];
    }
    free(S); free(s); free(lambda); free(U);
    return rc;
}

/* reghdfe's construction can be indefinite: the constant's row and column
 * come from the clipped full matrix, the coefficient block B from clipping
 * the original block on its own. Stata's ereturn post may then set e(V) to
 * zero, depending on rounding (with reghdfe, on the order of the cluster
 * variables). The repair keeps B, replaces the constant's covariances c by
 * their part in the column space of B, and raises its variance d to the
 * least value that makes V positive semidefinite, c' B^+ c, when it is below
 * that. */
static ST_retcode repair_constant(ST_double *V, ST_int K_with_cons)
{
    const ST_int K = K_with_cons - 1, n = K_with_cons;
    ST_double *B = (ST_double *)malloc((size_t)K * K * sizeof(*B));
    ST_double *lambda = (ST_double *)malloc((size_t)K * sizeof(*lambda));
    ST_double *U = (ST_double *)malloc((size_t)K * K * sizeof(*U));
    ST_double *c = (ST_double *)malloc((size_t)K * sizeof(*c));
    ST_double *v = (ST_double *)calloc((size_t)K, sizeof(*v));
    ST_retcode rc = (B && lambda && U && c && v) ? 0 : 920;
    if (!rc) {
        for (ST_int i = 0; i < K; i++) {
            c[i] = V[(size_t)i * n + K];
            for (ST_int j = 0; j < K; j++) B[(size_t)i * K + j] = V[(size_t)i * n + j];
        }
        rc = sym_eigen(B, K, lambda, U);
    }
    if (!rc) {
        ST_double lmax = 0.0;
        for (ST_int k = 0; k < K; k++) if (lambda[k] > lmax) lmax = lambda[k];
        const ST_double cut = lmax * K * 1e-13;
        ST_double q = 0.0;
        for (ST_int k = 0; k < K; k++) {
            if (!(lambda[k] > cut)) continue;
            ST_double proj = 0.0;
            for (ST_int i = 0; i < K; i++) proj += U[(size_t)i * K + k] * c[i];
            for (ST_int i = 0; i < K; i++) v[i] += U[(size_t)i * K + k] * proj;
            q += proj * proj / lambda[k];
        }
        for (ST_int i = 0; i < K; i++) V[(size_t)i * n + K] = V[(size_t)K * n + i] = v[i];
        if (V[(size_t)K * n + K] < q) V[(size_t)K * n + K] = q;
    }
    free(B); free(lambda); free(U); free(c); free(v);
    return rc;
}

/* reghdfe_fix_psd() with the constant last: when V (K_with_cons square) has
 * a negative eigenvalue its eigenvalues are clipped, and the coefficient
 * block is then the original block, clipped in turn if it needs it. When
 * that leaves V indefinite beyond rounding and V_psd is not NULL, V_psd gets
 * V with the constant's row and column repaired (repair_constant) and
 * *psd_made is 1; V itself stays reghdfe's. */
static ST_retcode fix_psd(ST_double *V, ST_int K_with_cons, ST_int *fixed,
                          ST_double *V_psd, ST_int *psd_made)
{
    if (psd_made) *psd_made = 0;
    const ST_int K = K_with_cons - 1;
    ST_double *block = (ST_double *)malloc((size_t)(K > 0 ? K : 1) * (K > 0 ? K : 1) * sizeof(*block));
    if (!block) return 920;
    for (ST_int i = 0; i < K; i++)
        for (ST_int j = 0; j < K; j++) block[(size_t)i * K + j] = V[(size_t)i * K_with_cons + j];
    ST_retcode rc = psd_clip(V, K_with_cons, fixed);
    if (!rc && *fixed && K > 0) {
        ST_int block_fixed = 0;
        rc = psd_clip(block, K, &block_fixed);
        for (ST_int i = 0; i < K && !rc; i++)
            for (ST_int j = 0; j < K; j++) V[(size_t)i * K_with_cons + j] = block[(size_t)i * K + j];
        ST_double min_eig = 0.0;
        if (!rc && V_psd) rc = scaled_min_eigen(V, K_with_cons, &min_eig);
        if (!rc && V_psd && min_eig < -1e-8) {
            memcpy(V_psd, V, (size_t)K_with_cons * K_with_cons * sizeof(*V_psd));
            rc = repair_constant(V_psd, K_with_cons);
            if (!rc && psd_made) *psd_made = 1;
        }
    }
    free(block);
    return rc;
}

ST_retcode compute_vce_multicluster(
    const ST_double *data, const ST_double *resid, const ST_double *inv_xx,
    const ST_double *weights, ST_int weight_type, const ST_int *cluster_ids,
    const ST_int *num_clusters, ST_int Q, ST_int N, ST_double N_eff, ST_int K_with_cons,
    ST_double *V, ST_int df_m, ST_double df_a, ST_double df_a_nested, ST_int *psd_fixed,
    ST_double *V_psd, ST_int *psd_made)
{
    *psd_fixed = 0;
    if (psd_made) *psd_made = 0;
    if (Q < 2 || Q > 10) return 198;
    ST_int N_clust = num_clusters[0];
    for (ST_int d = 1; d < Q; d++) if (num_clusters[d] < N_clust) N_clust = num_clusters[d];
    if (N_clust <= 1) return 498;

    /* One small-sample factor for every term, with the fewest clusters */
    ST_double df_a_adjusted = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;
    if (!isfinite(N_eff)) return 498;
    if (N_eff <= nested_adj + df_m + df_a_adjusted) return zero_vce(K_with_cons, V);
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (N_eff - nested_adj - df_m - df_a_adjusted)) *
                        ((ST_double)N_clust / (N_clust - 1));

    ctools_vce_data d;
    d.X_eff = data + N;
    d.D = inv_xx;
    d.resid = resid;
    d.weights = weights;
    d.weight_type = weight_type;
    d.N = N;
    d.K = K_with_cons;
    d.normalize_weights = 1;

    const size_t KK = (size_t)K_with_cons * K_with_cons;
    ST_int *ids = (ST_int *)malloc((size_t)N * sizeof(*ids));
    ST_int *joined = (ST_int *)malloc((size_t)N * sizeof(*joined));
    ST_double *term = (ST_double *)malloc(KK * sizeof(*term));
    ST_retcode rc = (ids && joined && term) ? 0 : 920;
    memset(V, 0, KK * sizeof(*V));
    /* Cameron, Gelbach and Miller: every intersection of the cluster
     * variables, added with an odd number of them and subtracted with an
     * even number */
    for (ST_int subset = 1; subset < (1 << Q) && !rc; subset++) {
        ST_int members = 0, groups = 0;
        for (ST_int dim = 0; dim < Q && !rc; dim++) {
            if (!(subset & (1 << dim))) continue;
            const ST_int *dim_ids = cluster_ids + (size_t)dim * N;
            if (!members++) {
                memcpy(ids, dim_ids, (size_t)N * sizeof(*ids));
                groups = num_clusters[dim];
            } else {
                if (ctools_intersect_cluster_ids(ids, dim_ids, N, joined, &groups)) rc = 920;
                else memcpy(ids, joined, (size_t)N * sizeof(*ids));
            }
        }
        if (rc) break;
        rc = ctools_vce_cluster(&d, ids, groups, dof_adj, term);
        ST_double sign = (members % 2) ? 1.0 : -1.0;
        for (size_t i = 0; i < KK && !rc; i++) V[i] += sign * term[i];
    }
    free(ids); free(joined); free(term);
    if (rc) return rc;
    make_symmetric(K_with_cons, V);
    return fix_psd(V, K_with_cons, psd_fixed, V_psd, psd_made);
}

ST_retcode compute_vce_dkraay(
    const ST_double *data, const ST_double *resid, const ST_double *inv_xx,
    const ST_double *weights, ST_int weight_type, const ST_int *time_ids, ST_int T,
    ST_int bandwidth, ST_int N, ST_double N_eff, ST_int K_with_cons, ST_double *V,
    ST_int df_m, ST_double df_a, ST_int *bandwidth_used, ST_double *V_psd, ST_int *psd_made)
{
    if (psd_made) *psd_made = 0;
    const ST_int K = K_with_cons;
    if (T <= 1) return 498;
    /* Bandwidth = lags + 1; reghdfe's default is Newey and West's (1994)
     * floor(4 (T/100)^(2/9)) lags, and at most T - 1 lags are used */
    ST_int lags = bandwidth > 0 ? bandwidth - 1 : (ST_int)floor(4.0 * pow(T / 100.0, 2.0 / 9.0));
    if (lags < 0) lags = 0;
    if (lags >= T) lags = T - 1;
    *bandwidth_used = lags + 1;
    if (!isfinite(N_eff)) return 498;
    if (N_eff <= df_m + df_a) return zero_vce(K, V);
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (N_eff - df_m - df_a)) * ((ST_double)T / (T - 1));

    /* Weighted residuals as for the clustered VCE: fweights as they are,
     * aweights and pweights normalized to sum to N */
    ST_double wscale = 1.0;
    if (weights && (weight_type == 1 || weight_type == 3)) {
        ST_double sum_w = 0.0;
        for (ST_int i = 0; i < N; i++) sum_w += weights[i];
        if (!isfinite(sum_w) || sum_w <= 0) return 498;
        wscale = (ST_double)N / sum_w;
    }
    long double *h = (long double *)calloc((size_t)T * K, sizeof(*h));
    long double *M = (long double *)calloc((size_t)K * K, sizeof(*M));
    ST_double *DM = (ST_double *)malloc((size_t)K * K * sizeof(*DM));
    if (!h || !M || !DM) { free(h); free(M); free(DM); return 920; }
    /* h_t: the moments of period t summed over its observations */
    for (ST_int i = 0; i < N; i++) {
        ST_int t = time_ids[i];
        if (t < 0 || t >= T) { free(h); free(M); free(DM); return 498; }
        ST_double e = resid[i] * (weights ? weights[i] * wscale : 1.0);
        for (ST_int k = 0; k < K; k++) h[(size_t)t * K + k] += (long double)e * data[(size_t)(k + 1) * N + i];
    }
    /* Meat: sum_t h_t' h_t plus the Bartlett-weighted autocovariances */
    for (ST_int t = 0; t < T; t++)
        for (ST_int a = 0; a < K; a++)
            for (ST_int b = 0; b < K; b++) M[(size_t)a * K + b] += h[(size_t)t * K + a] * h[(size_t)t * K + b];
    for (ST_int j = 1; j <= lags; j++) {
        long double wj = 1.0L - (long double)j / (lags + 1);
        for (ST_int a = 0; a < K; a++) {
            for (ST_int b = 0; b < K; b++) {
                long double om = 0.0L, om_t = 0.0L;
                for (ST_int t = j; t < T; t++) {
                    om += h[(size_t)t * K + a] * h[(size_t)(t - j) * K + b];
                    om_t += h[(size_t)t * K + b] * h[(size_t)(t - j) * K + a];
                }
                M[(size_t)a * K + b] += wj * (om + om_t);
            }
        }
    }
    /* V = D M D * dof_adj */
    for (ST_int a = 0; a < K; a++)
        for (ST_int b = 0; b < K; b++) {
            long double s = 0.0L;
            for (ST_int c = 0; c < K; c++) s += (long double)inv_xx[(size_t)a * K + c] * M[(size_t)c * K + b];
            DM[(size_t)a * K + b] = (ST_double)s;
        }
    for (ST_int a = 0; a < K; a++)
        for (ST_int b = 0; b < K; b++) {
            long double s = 0.0L;
            for (ST_int c = 0; c < K; c++) s += (long double)DM[(size_t)a * K + c] * inv_xx[(size_t)c * K + b];
            V[(size_t)a * K + b] = (ST_double)(s * dof_adj);
        }
    free(h); free(M); free(DM);
    make_symmetric(K, V);
    ST_int fixed = 0;
    return fix_psd(V, K, &fixed, V_psd, psd_made);
}
