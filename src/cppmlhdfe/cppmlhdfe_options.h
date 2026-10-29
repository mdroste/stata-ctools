#ifndef CPPMLHDFE_OPTIONS_H
#define CPPMLHDFE_OPTIONS_H
typedef struct
{
    int standardize, min_ok, exact_solver, exact_partial, heuristic, remove_collinear;
    int step_halving, max_halving, guess, guess_idx, check_mu, check_simplex;
    int relu_maxiter, relu_strict, relu_accelerate, simplex_maxiter;
    int tag_idx, certificate_idx;
    double start_tol, halving_memory, mu_tol, relu_tol, relu_zero_tol, simplex_tol;
    int acceleration, transform, pool;
    double btol, conlim;
    int ref_exact_partial; /* exact partialling explicitly requested */
} PPML_Options;

#define PPML_OPTIONS_DEFAULT                                                                       \
    {1,    1, 0, 1,    1,  1,    0,    2,    0,     0, 0, 0, 100,  0,  0,                          \
     1000, 0, 0, 1e-4, .9, 1e-6, 1e-4, 1e-8, 1e-12, 0, 0, 0, 1e-8, 1e8, 0}
void ppml_read_options(PPML_Options *options);
double ppml_option(const char *name, double fallback);
#endif
