/*
    civreghdfe_ffirst.h
    ivreg2's first-stage and weak-instrument statistics (first, ffirst,
    sfirst): the system of reduced-form and first-stage regressions, the
    e(first) matrix (partial and Shea R2, F, Sanderson-Windmeijer and
    Angrist-Pischke statistics), Anderson-Rubin and Stock-Wright tests.
*/

#ifndef CIVREGHDFE_FFIRST_H
#define CIVREGHDFE_FFIRST_H

#include "../stplugin.h"
#include "civreghdfe_omega.h"

/* Rows of e(first), in ivreg2's order */
#define CIV_FFIRST_ROWS 21

typedef struct {
    ST_double N;          /* ivreg2's N: sum of fweights, else #obs */
    ST_int dofminus;      /* large-sample dofminus() */
    ST_int sdofminus;     /* small-sample adjustment (absorbed effects, partialled-out regressors) */
    ST_int L_dof;         /* instrument count in the dof (absorbed columns included) */
    ST_int cluster;       /* clustered S (one- or two-way, Driscoll-Kraay) */
    ST_int nclust;        /* clusters (two-way: the smaller dimension) */
    ST_int robust;        /* non-iid S: Wald statistics use the robust VCE */
} civ_ffirst_dof;

/*
    Compute the statistics and store them for the ado:
      matrix __civreghdfe_firstmat (21 x K1), __civreghdfe_sfb (1 x (K1+1)L),
      __civreghdfe_sfV and __civreghdfe_sfS ((K1+1)L square), with the
      equations [y, X1...] and, within each, the instruments in ivreg2's
      order (excluded, then included); scalars __civreghdfe_ar*,
      __civreghdfe_sstat*, __civreghdfe_rmse_rf, __civreghdfe_ff_K1 and
      __civreghdfe_ff_L. The matrices must exist with at least these
      dimensions (the ado creates them).
    y, X1 (N x K1), Z (N x L, the K2 included instruments first, then the
    excluded ones); X2 is the first K2 columns of Z. Returns 0, 920 or a
    store error; a statistic that cannot be computed is missing.
*/
ST_retcode civ_ffirst(const ST_double *y, const ST_double *X1, ST_int K1,
                      const ST_double *Z, ST_int L, ST_int K2,
                      const ST_double *w, ST_int N, const civ_omega *om,
                      const civ_ffirst_dof *dof);

#endif /* CIVREGHDFE_FFIRST_H */
