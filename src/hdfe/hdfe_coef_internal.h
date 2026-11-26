/*
 * hdfe_coef_internal.h
 *
 * Private interfaces of the coefficient-space projection (hdfe_coef.h)
 * Part of the ctools Stata plugin suite
 *
 * Conventions: levels are 1-based as in FE_Factor; coefficient vectors are
 * 0-based and hold the K columns of a batch interleaved, [level * K + k].
 * Every function is reentrant. Allocations use malloc/calloc/free only (the
 * native harnesses interpose exactly those).
 */

#ifndef CTOOLS_HDFE_COEF_INTERNAL_H
#define CTOOLS_HDFE_COEF_INTERNAL_H

#include <stdint.h>

#include "hdfe_coef.h"
#include "../ctools_cancel.h"

/* Columns solved in one batch: bounds the interleaved vectors and the
 * kernels' stack arrays. More columns run in successive batches on the same
 * structure. */
#define HDFE_COEF_KMAX 16

/* Private per-chunk buffers of level sums may use this much memory. */
#define HDFE_COEF_BUFFER_BYTES ((size_t)32 << 20)

/* Loops shorter than this many element operations run on one thread. */
#define HDFE_COEF_PARALLEL_MIN 100000

#if defined(__GNUC__) || defined(__clang__)
#define HDFE_COEF_KERNEL static inline __attribute__((always_inline))
#else
#define HDFE_COEF_KERNEL static inline
#endif

/* Run CALL with KK a compile-time constant equal to K, so that the kernels'
 * loops over the columns of a batch have fixed trip counts. */
#define HDFE_COEF_K_CASE(n, CALL)                                                                  \
    case n:                                                                                        \
    {                                                                                              \
        enum                                                                                       \
        {                                                                                          \
            KK = n                                                                                 \
        };                                                                                         \
        CALL;                                                                                      \
    }                                                                                              \
    break;
#define HDFE_COEF_K_DISPATCH(K, CALL)                                                              \
    switch (K)                                                                                     \
    {                                                                                              \
        HDFE_COEF_K_CASE(1, CALL)                                                                  \
        HDFE_COEF_K_CASE(2, CALL)                                                                  \
        HDFE_COEF_K_CASE(3, CALL)                                                                  \
        HDFE_COEF_K_CASE(4, CALL)                                                                  \
        HDFE_COEF_K_CASE(5, CALL)                                                                  \
        HDFE_COEF_K_CASE(6, CALL)                                                                  \
        HDFE_COEF_K_CASE(7, CALL)                                                                  \
        HDFE_COEF_K_CASE(8, CALL)                                                                  \
        HDFE_COEF_K_CASE(9, CALL)                                                                  \
        HDFE_COEF_K_CASE(10, CALL)                                                                 \
        HDFE_COEF_K_CASE(11, CALL)                                                                 \
        HDFE_COEF_K_CASE(12, CALL)                                                                 \
        HDFE_COEF_K_CASE(13, CALL)                                                                 \
        HDFE_COEF_K_CASE(14, CALL)                                                                 \
        HDFE_COEF_K_CASE(15, CALL)                                                                 \
    default:                                                                                       \
    {                                                                                              \
        enum                                                                                       \
        {                                                                                          \
            KK = HDFE_COEF_KMAX                                                                    \
        };                                                                                         \
        CALL;                                                                                      \
    }                                                                                              \
    break;                                                                                         \
    }

/* Threads a parallel region opened here would get */
static inline int hdfe_coef_threads(void)
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/* ------------------------------------------------------------------------
 * Rows grouped by one effect
 * ---------------------------------------------------------------------- */

/*
 * rowptr (nlev + 1, owned by the caller afterwards): the rows of level l
 * occupy positions rowptr[l-1] .. rowptr[l]-1 of the grouped order, in their
 * original order. When the rows already arrive grouped in ascending level
 * order, *sorted is 1 and *dest is NULL; otherwise dest[i] (N, owned by the
 * caller) is the grouped position of row i.
 *
 * With second != NULL the same pass checks 1 <= second[i] <= nsecond and
 * sets *increasing when, within every level of an already grouped input, the
 * second effect's levels strictly increase (so no cell repeats).
 *
 * Returns 0, 198 (a level out of range) or 920.
 */
ST_retcode hdfe_coef_group_rows(ST_int N, const ST_int *levels, ST_int nlev, const ST_int *second,
                                ST_int nsecond, ST_int **rowptr, ST_int **dest, int *sorted,
                                int *increasing);

/* ------------------------------------------------------------------------
 * Cross-tabulation of two effects: one entry per distinct cell
 * ---------------------------------------------------------------------- */

typedef struct
{
    ST_int nrow, ncol;    /* levels of the row effect and of the column effect */
    ST_int nnz;           /* distinct cells */
    ST_int *rowptr;       /* nrow + 1 */
    const ST_int *col;    /* 1-based level of the column effect of each cell */
    const ST_double *val; /* cell weights; NULL when every cell weighs one */
    ST_int *col_own;      /* storage behind col, NULL when it borrows the level array */
    ST_double *val_own;   /* likewise for val */
    int sorted;           /* the rows arrived grouped by the row effect */
} hdfe_crosstab;

/* Cells of (row_levels, col_levels) with weights w (NULL: counts). The
 * result borrows col_levels and w when the rows arrive grouped by the row
 * effect and no cell repeats: both must outlive it. Returns 0, 198 or 920;
 * C is freed on failure. */
ST_retcode hdfe_crosstab_build(hdfe_crosstab *C, ST_int N, const ST_int *row_levels, ST_int nrow,
                               const ST_int *col_levels, ST_int ncol, const ST_double *w);
void hdfe_crosstab_free(hdfe_crosstab *C);

/* Fixed boundaries (nchunks + 1 entries, owned by the caller) that split
 * rows[0..n) (or 0..n-1 when rows == NULL) into chunks with about equal
 * numbers of cells of rowptr. NULL on allocation failure. */
ST_int *hdfe_coef_chunks(const ST_int *rowptr, const ST_int *rows, ST_int n, int nchunks);

/* ------------------------------------------------------------------------
 * A cross-tabulation kept while the weights change (hdfe_coef_plan)
 *
 * The cells of two effects depend on their levels alone; only the cell
 * weights follow the row weights. hdfe_crossplan_build() finds the cells
 * once and records where every row goes, and hdfe_crossplan_values() then
 * gives, for any row weights, the cross-tabulation hdfe_crosstab_build()
 * would return for them: the same cells in the same order, every cell's
 * weight the sum of its rows' in the same order.
 * ---------------------------------------------------------------------- */

typedef struct
{
    hdfe_crosstab C;   /* a view: its arrays are the ones below or the caller's */
    ST_int N;
    ST_int *grouprow;  /* nrow + 1: the rows of each level of the row effect, grouped */
    ST_int *cellptr;   /* nrow + 1 when cells repeat (then C.rowptr), else NULL */
    ST_int *order;     /* N: the row at each grouped position; NULL: rows arrive grouped */
    ST_int *cell;      /* N: the cell of each grouped position; NULL: one row per cell */
    ST_int *col;       /* storage behind C.col, NULL when it borrows the level array */
    ST_double *val;    /* storage behind C.val, NULL when it borrows the weights */
    ST_int *chunk;     /* nchunk + 1: levels of the row effect, about equal rows each */
    int nchunk;
} hdfe_crossplan;

/* Returns 0, 198 or 920; P is freed on failure. col_levels must outlive P. */
ST_retcode hdfe_crossplan_build(hdfe_crossplan *P, ST_int N, const ST_int *row_levels,
                                ST_int nrow, const ST_int *col_levels, ST_int ncol);
/* Sets P->C.val for the row weights w (N, not NULL; borrowed when the rows
 * arrive grouped and no cell repeats, so w must outlive the use of P->C) */
void hdfe_crossplan_values(hdfe_crossplan *P, const ST_double *w);
void hdfe_crossplan_free(hdfe_crossplan *P);

/* A buffer that grows to the largest size asked of it; contents are not kept */
typedef struct
{
    void *p;
    size_t bytes;
} hdfe_coef_buf;
static inline void *hdfe_coef_buf_get(hdfe_coef_buf *b, size_t bytes)
{
    if (!b->p || bytes > b->bytes)
    {
        free(b->p);
        b->p = malloc(bytes ? bytes : 1);
        b->bytes = b->p ? bytes : 0;
    }
    return b->p;
}
static inline void hdfe_coef_buf_free(hdfe_coef_buf *b)
{
    free(b->p);
    b->p = NULL;
    b->bytes = 0;
}

/* What hdfe_coef_partial_out_plan() keeps between calls: the factors it was
 * built for, and the structure of the two-effect or the multi-effect solver
 * with its release function. */
struct hdfe_coef_plan
{
    ST_int N;
    int G;
    const ST_int *levels[HDFE_COEF_MAX_EFFECTS];
    ST_int nlev[HDFE_COEF_MAX_EFFECTS];
    void *kept;
    void (*kept_free)(void *kept);
};
/* Records S's factors as the ones plan->kept describes */
void hdfe_coef_plan_bind(hdfe_coef_plan *plan, const HDFE_State *S, void *kept,
                         void (*kept_free)(void *kept));

/* ------------------------------------------------------------------------
 * Passes over the data rows
 *
 * B is the eliminated effect; the coefficients of the nk kept effects are
 * stacked, effect h starting at off[h], m = off[nk] in all.
 * ---------------------------------------------------------------------- */

typedef struct
{
    ST_int N, nB, m;
    int nk;
    int sorted;                               /* rows arrive grouped by B */
    const ST_int *lvB;                        /* 1-based level of B of every row */
    const ST_int *lvk[HDFE_COEF_MAX_EFFECTS]; /* 1-based level of each kept effect */
    ST_int off[HDFE_COEF_MAX_EFFECTS + 1];
    const ST_double *w;    /* row weights, NULL when unweighted */
    const ST_double *invB; /* nB: 1 / weighted count of each level */
} hdfe_coef_rows;

/* Chunks of rows that accumulate level sums into private buffers */
#define HDFE_COEF_ROW_CHUNKS 64
#define HDFE_COEF_ROW_CHUNK_ROWS 16384

/* How the accumulating passes run for batches of up to K columns: over
 * fixed chunks of rows, each with a private buffer of level sums, while
 * those buffers fit HDFE_COEF_BUFFER_BYTES; otherwise one column per
 * thread. Either way every sum is added in row order within a chunk (or
 * column) and the chunks are combined in order, so the results do not depend
 * on the thread count. */
typedef struct
{
    int nrc; /* chunks of rows; 0: one column per thread */
    ST_int bound[HDFE_COEF_ROW_CHUNKS + 1];
    int private_b;    /* the level sums of B need a buffer per chunk (rows not grouped) */
    ST_double *tbuf;  /* nrc * nB * K when private_b */
    ST_double *gbuf;  /* nrc * m * K, or (nB + m) * K column-major when nrc == 0 */
    ST_double *ssbuf; /* nrc * 2 * K */
    size_t tbytes, gbytes, ssbytes; /* what the three buffers hold */
} hdfe_coef_rowplan;

/* R is all zeros before its first use. Initializing it again (for another K,
 * on the same rows) keeps the buffers that are large enough; on failure
 * (920) R is freed. */
ST_retcode hdfe_coef_rowplan_init(hdfe_coef_rowplan *R, const hdfe_coef_rows *D, int K);
void hdfe_coef_rowplan_free(hdfe_coef_rowplan *R);

/* For K columns X (column-major, N rows): a0 (nB x K) = level means of B,
 * g (m x K) = Dk' W (x - D_B a0), f0[k] = (x - D_B a0)' W (x - D_B a0) and
 * scale[k] = x' W x. Working on the B-demeaned column keeps level means much
 * larger than the within variation from costing digits. */
void hdfe_coef_rows_rhs(const hdfe_coef_rows *D, const hdfe_coef_rowplan *R, const ST_double *X,
                        int K, ST_double *a0, ST_double *g, ST_double *f0, ST_double *scale);

/* X -= D_B a + Dk beta */
void hdfe_coef_rows_residual(const hdfe_coef_rows *D, ST_double *X, int K, const ST_double *a,
                             const ST_double *beta);

/* ------------------------------------------------------------------------
 * Preconditioned conjugate gradients for several right-hand sides
 * ---------------------------------------------------------------------- */

/* out = A in; both n x K interleaved */
typedef void (*hdfe_coef_apply_fn)(void *ctx, const ST_double *in, ST_double *out);

typedef struct
{
    ST_int n;                 /* unknowns per column */
    int K;                    /* columns */
    hdfe_coef_apply_fn apply; /* out = S in */
    void *apply_ctx;
    hdfe_coef_apply_fn precond; /* out = M^-1 in; NULL: the diagonal minv */
    void *precond_ctx;
    const ST_double *minv; /* n, used when precond == NULL */
    ST_double tolerance;
    ST_int maxiter;
    int warm; /* x holds a starting point; otherwise it is zeroed */
    /* Iterations allowed beyond maxiter while every unfinished column is
     * within HDFE_COEF_NEAR of its stopping rule (a trial run that is about
     * to finish is not cut short) */
    ST_int extend;
} hdfe_coef_pcg_problem;

/* A column is about to converge when its last step's decrease is at most
 * this multiple of what the stopping rule accepts: two orders of magnitude
 * in the error, a dozen steps of an iteration that needs sixty in all. */
#define HDFE_COEF_NEAR 1e4

/* Doubles of scratch hdfe_coef_pcg() needs */
#define HDFE_COEF_PCG_BLOCK 8192
#define HDFE_COEF_PCG_WORK(n, K)                                                                   \
    ((size_t)4 * (size_t)(n) * (size_t)(K) +                                                       \
     (((size_t)(n) + HDFE_COEF_PCG_BLOCK - 1) / HDFE_COEF_PCG_BLOCK) * (size_t)(K))

/*
 * How the solve of one column of a batch ended. A solve in coefficient
 * space carries rounding relative to the column it starts from, while the
 * tolerance is relative to what is left once the effects are removed. A
 * column that the effects explain almost entirely (level means far larger
 * than the within variation, say) therefore reaches its rounding first
 * (ROUNDING); solving once more from the residuals just written, a much
 * smaller column, recovers the digits. The drivers do that for at most
 * HDFE_COEF_ROUNDS rounds.
 */
/* not finished */
#define HDFE_COEF_COLUMN_ACTIVE 0
/* to the tolerance, or to rounding where no further round can help */
#define HDFE_COEF_COLUMN_CONVERGED 1
/* to rounding; another round would gain digits */
#define HDFE_COEF_COLUMN_ROUNDING 2
#define HDFE_COEF_ROUNDS 3

/*
 * Solves S x = g for the K columns. f0[k] is the residual sum of squares of
 * column k at x = 0, the reference of the stopping rule (hdfe_coef.h), and
 * scale[k] the sum of squares of the column the right-hand side was formed
 * from, the reference of its rounding.
 * state (K) holds one of the values above per column: a column that is not
 * ACTIVE on entry is left as it is (its x too); on return every column
 * that finished has its outcome. work is HDFE_COEF_PCG_WORK(n, K) doubles
 * of scratch. Polls for cancellation once per iteration. Returns 0, 430
 * (maxiter reached with columns still ACTIVE; x holds the last iterate) or
 * 1; *iterations is the largest count over the columns.
 */
ST_retcode hdfe_coef_pcg(const hdfe_coef_pcg_problem *P, const ST_double *g, ST_double *x,
                         const ST_double *f0, const ST_double *scale, ST_double *work,
                         unsigned char *state, ST_int *iterations);

/* The outcome of a direct solve b of S b = g, for each of K columns:
 * weakest is the factor's smallest pivot relative to its diagonal. */
void hdfe_coef_direct_state(ST_int n, int K, const ST_double *g, const ST_double *b,
                            const ST_double *f0, const ST_double *scale, ST_double tolerance,
                            ST_double weakest, unsigned char *state);

/* The contiguous range of columns that covers every ROUNDING column of
 * state[0..K): on return *first is advanced to its start and *count is its
 * length. Returns 0 when there is none. */
int hdfe_coef_refine_span(const unsigned char *state, int K, int *first, int *count);

/* ------------------------------------------------------------------------
 * Dense Cholesky of a positive semidefinite matrix
 * ---------------------------------------------------------------------- */

/*
 * A: row-major n x n, lower triangle read and overwritten with the factor.
 * pin (n, or NULL): rows marked 1 are fixed at zero (the caller's choice of
 * one level per connected component). A row whose pivot falls to rounding
 * level is fixed at zero as well. dead[j] (n) receives 1 for every fixed row.
 * *weakest receives the smallest ratio of a pivot to its row's diagonal
 * before elimination (1 when no row is left free). Returns 0 or 920.
 */
ST_retcode hdfe_coef_chol_factor(ST_double *A, ST_int n, const unsigned char *pin,
                                 unsigned char *dead, ST_double *weakest);

/* Solves L L' x = b in place for K interleaved columns, with fixed rows zero. */
void hdfe_coef_chol_solve(const ST_double *L, ST_int n, const unsigned char *dead, ST_double *x,
                          int K);

/* ------------------------------------------------------------------------
 * Approximate Cholesky factor of a graph Laplacian
 * ---------------------------------------------------------------------- */

typedef struct hdfe_lap hdfe_lap;

/* Factor of the Laplacian with edges (eu[e], ev[e]) of weight ew[e] > 0 on
 * vertices 0..n-1. The edge arrays are scratch (they are rewritten).
 * Deterministic: the sampling uses a fixed seed. NULL when memory runs out. */
hdfe_lap *hdfe_lap_factor(ST_int n, ST_int m, ST_int *eu, ST_int *ev, ST_double *ew);
void hdfe_lap_free(hdfe_lap *L);

/* y = (approximate) pseudo-inverse of the Laplacian times y, in place, for K
 * interleaved columns. One vertex of each component is held at zero. */
void hdfe_lap_solve(const hdfe_lap *L, ST_double *y, int K);

/* The sampling rule of the factorization for a star: writes the len - 1
 * edges that replace the star with leaves cs[].row and weights cs[].val
 * (sorted ascending by weight on return) into eu/ev/ew and returns their
 * number. state is the generator's state. */
typedef struct
{
    ST_int entry;
    ST_int row;
    ST_double val;
} hdfe_lap_leaf;
ST_int hdfe_lap_sample_star(hdfe_lap_leaf *cs, ST_double *cum, ST_int len, uint64_t *state,
                            ST_int *eu, ST_int *ev, ST_double *ew);

/* ------------------------------------------------------------------------
 * Three or more effects (hdfe_coef_multi.c)
 * ---------------------------------------------------------------------- */

/* plan: NULL, or the plan to keep the structure in (S has weights then) */
HDFE_SolveResult hdfe_coef_multi_partial_out(const HDFE_State *S, hdfe_coef_plan *plan,
                                             ST_double *data, ST_int N, ST_int K,
                                             hdfe_coef_method method);

#endif /* CTOOLS_HDFE_COEF_INTERNAL_H */
