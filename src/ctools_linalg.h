#ifndef CTOOLS_LINALG_H
#define CTOOLS_LINALG_H
/* Row-major dense lower factor L, leading dimension n; diagonal must be nonzero.
 * No allocation, no refactorization, ordinary double arithmetic in index order.
 * b and x may alias exactly; L cannot alias either. n==0 touches no pointers.
 * These kernels do not translate statuses or apply singularity thresholds. */
static inline void ctools_solve_lower_factor(const double *L, const double *b, int n, double *x)
{
    for (int i = 0; i < n; i++)
    {
        double sum = b[i];
        for (int j = 0; j < i; j++)
            sum -= L[i * n + j] * x[j];
        x[i] = sum / L[i * n + i];
    }
}
static inline void ctools_solve_lower_transpose(const double *L, const double *b, int n, double *x)
{
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
            sum -= L[j * n + i] * x[j];
        x[i] = sum / L[i * n + i];
    }
}
#endif
