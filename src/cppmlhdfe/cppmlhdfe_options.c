#include "cppmlhdfe_options.h"
#include "stplugin.h"
#include <stdio.h>

double ppml_option(const char *name, double fallback)
{
    char key[64];
    double value;
    snprintf(key, sizeof(key), "__cppmlhdfe_%s", name);
    return SF_scal_use(key, &value) == 0 ? value : fallback;
}

void ppml_read_options(PPML_Options *options)
{
#define OPT(field, name, value) options->field = ppml_option(name, value)
    OPT(standardize, "standardize_data", 1);
    OPT(min_ok, "min_ok", 1);
    OPT(exact_solver, "use_exact_solver", 0);
    OPT(exact_partial, "use_exact_partial", 1);
    OPT(heuristic, "use_heuristic_tol", 1);
    OPT(remove_collinear, "remove_collinear", 1);
    OPT(step_halving, "use_step_halving", 0);
    OPT(max_halving, "max_step_halving", 2);
    OPT(guess, "guess", 0);
    OPT(guess_idx, "guess_idx", 0);
    OPT(check_mu, "check_mu", 0);
    OPT(check_simplex, "check_simplex", 0);
    OPT(relu_maxiter, "relu_maxiter", 100);
    OPT(relu_strict, "relu_strict", 0);
    OPT(relu_accelerate, "relu_accelerate", 0);
    OPT(simplex_maxiter, "simplex_maxiter", 1000);
    OPT(tag_idx, "tag_idx", 0);
    OPT(certificate_idx, "certificate_idx", 0);
    OPT(start_tol, "start_inner_tol", 1e-4);
    OPT(halving_memory, "step_halving_memory", .9);
    OPT(mu_tol, "mu_tol", 1e-6);
    OPT(relu_tol, "relu_tol", 1e-4);
    OPT(relu_zero_tol, "relu_zero_tol", 1e-8);
    OPT(simplex_tol, "simplex_tol", 1e-12);
    OPT(acceleration, "acceleration", 0);
    OPT(transform, "transform", 0);
    OPT(pool, "pool", 0);
    OPT(btol, "btol", 1e-8);
    OPT(conlim, "conlim", 1e8);
    OPT(ref_exact_partial, "ref_exact_partial", 0);
#undef OPT
}
