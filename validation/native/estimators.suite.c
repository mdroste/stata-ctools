/* Estimators tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_cppmlhdfe_native_source)
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "hdfe/hdfe_solver.c"
#include "cppmlhdfe/cppmlhdfe_irls.c"
#include "cppmlhdfe/cppmlhdfe_result.c"
#include "cppmlhdfe/cppmlhdfe_input.c"
#include "cppmlhdfe/cppmlhdfe_options.c"
#include "cppmlhdfe/cppmlhdfe_solvers.c"
#include "cppmlhdfe/cppmlhdfe_simplex.c"
#include "cppmlhdfe/cppmlhdfe_relu.c"
#include "cppmlhdfe/cppmlhdfe_moments.c"
static PPML_Options ppml_options = PPML_OPTIONS_DEFAULT;
static int message(char *s) { (void)s; return 0; }
static void init(HDFE_State *s, int n, int groups) {
 memset(s,0,sizeof(*s));s->N=n;s->G=groups;s->num_threads=1;
 s->tolerance=1e-12;s->maxiter=10000;s->has_weights=1;
 s->weights=calloc(n,sizeof(double));s->factors=calloc(groups,sizeof(FE_Factor));
 for(int g=0;g<groups;g++) {
  FE_Factor *f=&s->factors[g];f->num_levels=2;f->has_intercept=1;
  f->levels=calloc(n,sizeof(int));f->counts=calloc(2,sizeof(double));
  f->weighted_counts=calloc(2,sizeof(double));
  for(int i=0;i<n;i++){f->levels[i]=1+i%2;f->counts[i%2]++;}
 }
 assert(!ctools_hdfe_alloc_buffers(s,0,4));
}
static void destroy(HDFE_State *s) {free(s->weights);ctools_hdfe_state_cleanup(s);}
int main(void) {
 mock.spouterr=message;mock.spoutsml=message;
 for(int run=0;run<30;run++) {
  HDFE_State s;init(&s,60,1);
  double y[60],x[120];int sep[60];
  for(int i=0;i<60;i++){y[i]=i<15?0:1;x[i]=sin(i);x[60+i]=x[i]+(i<15);}
  assert(ppml_relu(&ppml_options, &s,y,x,60,2,sep,1)==15);
  for(int i=0;i<60;i++)assert(sep[i]==(i<15));
  /* A tiny positive outcome must impose the equality constraint too. */
  y[0]=1e-200;assert(ppml_relu(&ppml_options, &s,y,x,60,2,sep,1)==0);
  /* Zeros distributed over a full-rank interior are not separation. */
  for(int i=0;i<60;i++){y[i]=(i%5)?1:0;x[60+i]=cos(i);}
  assert(ppml_relu(&ppml_options, &s,y,x,60,2,sep,1)==0);
  destroy(&s);
 }
 HDFE_State s;init(&s,60,2);
 double y[60],x[60];int sep[60];
 for(int i=0;i<60;i++) {
  s.factors[0].levels[i]=1+(i>=20);
  s.factors[1].levels[i]=1+(i>=40);
  y[i]=(i>=20 && i<40)?0:1;x[i]=sin(i);
 }
 assert(ppml_relu(&ppml_options, &s,y,x,60,1,sep,1)==20);
 for(int i=0;i<60;i++)assert(sep[i]==(i>=20&&i<40));
 destroy(&s);
 init(&s,60,2);
 s.factors[1].slope=calloc(60,sizeof(double));
 s.factors[1].slope_center=calloc(2,sizeof(double));
 for(int i=0;i<60;i++) {
  s.factors[1].slope[i]=i/60.0;
  s.weights[i]=1+(i%3);
  y[i]=(i%2?3:2)+(i%2?-2:4)*s.factors[1].slope[i];
 }
 ppml_update_fe_weights(&s,s.weights,60);
 HDFE_SolveResult r=partial_out_columns(&s,y,60,1,1);
 assert(!r.status);
 for(int i=0;i<60;i++)assert(fabs(y[i])<1e-9);
 destroy(&s);
 {
  double a[]={1,0, 0,1, 1,1};int mask[4]={0};
  assert(ppml_simplex_tableau(a,3,2,100,mask)==3);
  double b[]={1,0, -1,0, 0,1, 0,-1};
  assert(ppml_simplex_tableau(b,4,2,100,mask)==0);
  double c[]={1,1, 1,-1, 0,1, 0,-1};
  assert(ppml_simplex_tableau(c,4,2,100,mask)==2);
  assert(mask[0] && mask[1] && !mask[2] && !mask[3]);
  double duplicate[]={1,1, 2,2, -1,-1};
  assert(ppml_simplex_tableau(duplicate,3,2,100,mask)==0);
 }
 init(&s,60,2);
 for(int i=0;i<60;i++) {
  s.factors[0].levels[i]=1+(i%2);
  s.factors[1].levels[i]=1+(i>=30);
  s.weights[i]=1+(i%3);y[i]=3*(i%2)+7*(i>=30)+sin(i);
 }
 ppml_update_fe_weights(&s,s.weights,60);
 double baseline[60];memcpy(baseline,y,sizeof(y));
 assert(!partial_out_columns(&s,baseline,60,1,1).status);
 ppml_options.btol=1e-12;ppml_options.conlim=1e8;
 assert(!ppml_lsmr(&ppml_options, &s,y,60,1).status);
 for(int i=0;i<60;i++)assert(fabs(y[i]-baseline[i])<1e-9);
 for(int method=0;method<5;method++)for(int transform=0;transform<3;transform++) {
  if(transform==1 && (method==0 || method==4))continue;
  for(int i=0;i<60;i++)y[i]=3*(i%2)+7*(i>=30)+sin(i);
  ppml_options.acceleration=method;ppml_options.transform=transform;
  assert(!ppml_map_options(&ppml_options, &s,y,60,1).status);
  for(int i=0;i<60;i++)assert(fabs(y[i]-baseline[i])<1e-8);
 }
 destroy(&s);
 return 0;
}
#endif

#if defined(CTOOLS_CASE_cqreg_native)
/* Bounded cqreg checks and timings; built by test_cqreg_native.py. */
#include <assert.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ctools_config.h"
#include "cqreg_types.h"
#include "cqreg_linalg.h"
#include "cqreg_sparsity.h"
#include "cqreg_blas.h"
/* Include the solver to exercise the numerical QR retry directly as well. */
#include "cqreg_fn.c"

ST_plugin *_stata_;
static size_t allocated_bytes;
int cqreg_record_alloc(void **p, size_t alignment, size_t bytes) {
    int rc = posix_memalign(p, alignment, bytes);
    if (rc == 0) allocated_bytes += bytes;
    return rc;
}
int ctools_get_max_threads(void) { return 4; }
void stata_data_free(stata_data *data) { (void)data; }
double ctools_timer_seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
void ctools_msg(const char *module, const char *fmt, ...) {
    (void)module;
    va_list args; va_start(args, fmt); vfprintf(stderr, fmt, args); va_end(args);
    fputc('\n', stderr);
}
void ctools_error(const char *module, const char *fmt, ...) {
    (void)module;
    va_list args; va_start(args, fmt); vfprintf(stderr, fmt, args); va_end(args);
    fputc('\n', stderr);
}
static uint64_t seed = 42;
static double uniform(void) {
    seed = seed * UINT64_C(6364136223846793005) + 1;
    return (seed >> 11) * 0x1p-53;
}
static int compare(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static void selection_checks(void) {
    const int sizes[] = {1, 2, 3, 9, 63, 64, 65, 257, 8193};
    const double quantiles[] = {0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999};
    for (size_t s = 0; s < sizeof(sizes)/sizeof(*sizes); s++) {
        int n = sizes[s];
        double *y = malloc(n*sizeof(*y)), *sorted = malloc(n*sizeof(*sorted));
        double *copy = malloc(n*sizeof(*copy));
        assert(y && sorted && copy);
        for (int pattern = 0; pattern < 6; pattern++) {
            for (int i = 0; i < n; i++) {
                y[i] = pattern == 0 ? uniform() : pattern == 1 ? i :
                    pattern == 2 ? n-i : pattern == 3 ? 7 :
                    pattern == 4 ? i%5 : (i < n/2 ? i : n-i);
            }
            memcpy(copy,y,n*sizeof(*y)); memcpy(sorted,y,n*sizeof(*y));
            qsort(sorted,n,sizeof(*sorted),compare);
            for (size_t q = 0; q < sizeof(quantiles)/sizeof(*quantiles); q++) {
                int k = (int)floor((n+1.0)*quantiles[q]);
                if (k < 1) k = 1;
                if (k > n) k = n;
                assert(cqreg_compute_quantile(y,n,quantiles[q]) == sorted[k-1]);
                assert(memcmp(y,copy,n*sizeof(*y)) == 0);
#ifndef CQREG_BASELINE
                memcpy(copy,y,n*sizeof(*y));
                assert(cqreg_select(copy,n,k-1) == sorted[k-1]);
                for (int i=0;i<k-1;i++) assert(copy[i]<=copy[k-1]);
                for (int i=k;i<n;i++) assert(copy[i]>=copy[k-1]);
                memcpy(copy,y,n*sizeof(*y));
#endif
            }
        }
        free(y); free(sorted); free(copy);
    }
}
static void gemv_checks(void) {
    double a[] = {1,2,3,4,5,6}, x[] = {2,-1,3}, y[] = {NAN,NAN,NAN};
    blas_dgemv(0,3,2,1,a,3,x,0,y);
    assert(y[0] == -2 && y[1] == -1 && y[2] == 0);
    y[0] = y[1] = NAN;
    blas_dgemv(1,3,2,1,a,3,x,0,y);
    assert(y[0] == 9 && y[1] == 21);
    blas_dgemv(1,3,2,2,a,3,x,.5,y);
    assert(y[0] == 22.5 && y[1] == 52.5);
}
static void solver_checks(void) {
    enum {N=303,K=3};
    double x[N*K], y[N], b[K];
    for (int group = 0; group < N/3; group++) {
        double a = 2*uniform()-1, c = 2*uniform()-1;
        for (int e = 0; e < 3; e++) {
            int i=3*group+e;
            x[i]=a; x[N+i]=c; x[2*N+i]=1;
            y[i]=2*a-3*c+4+e-1;
        }
    }
    for (int qr=0; qr<2; qr++) for (int aux=0; aux<2; aux++) {
        cqreg_ipm_config c; cqreg_ipm_config_init(&c);
        c.tol_gap=1e-12; c.maxiter=100; c.skip_crossover=aux;
        cqreg_ipm_state *s=cqreg_ipm_create_lite(N,K,&c); assert(s);
        for (int q=1;q<=9;q+=4) {
            int rc=cqreg_fn_solve_core(s,y,x,q/10.0,b,qr);
            assert(rc>0 && s->converged);
            assert(fabs(b[0]-2)<1e-8 && fabs(b[1]+3)<1e-8);
            assert(fabs(b[2]-(q==1?3:q==5?4:5))<1e-8);
            if (!aux) for (int i=0;i<N;i++) {
                double r=y[i]-x[i]*b[0]-x[N+i]*b[1]-b[2];
                assert(fabs((s->u[i]-s->v[i])-r)<1e-10);
            }
        }
        s->config.maxiter=1;
        assert(cqreg_fn_solve(s,y,x,.5,b)<0 && !s->converged);
        cqreg_ipm_free(s);
    }
}
static void benchmark(int n, int k) {
    double *x=malloc((size_t)n*k*sizeof(*x)), *y=malloc(n*sizeof(*y));
    double *b=malloc(k*sizeof(*b)); assert(x && y && b);
    seed=42;
    for (int i=0;i<n;i++) {
        y[i]=2*uniform()-1;
        for (int j=0;j<k;j++) {
            x[(size_t)j*n+i]=j==k-1?1:2*uniform()-1;
            y[i]+=(j+1)*x[(size_t)j*n+i];
        }
    }
    double start=ctools_timer_seconds();
    double q=cqreg_compute_quantile(y,n,.5);
    printf("quantile,%d,%d,%.9f,%.17g\n",n,k,ctools_timer_seconds()-start,q);
    cqreg_ipm_config c; cqreg_ipm_config_init(&c);
    c.tol_gap=1e-12; c.maxiter=100;
    allocated_bytes=0;
#ifdef CQREG_BASELINE
    cqreg_ipm_state *s=cqreg_ipm_create(n,k,&c);
#else
    cqreg_ipm_state *s=cqreg_ipm_create_lite(n,k,&c);
#endif
    assert(s);
    size_t main_bytes=allocated_bytes;
    c.skip_crossover=1;
    allocated_bytes=0;
    cqreg_ipm_state *aux=cqreg_ipm_create_lite(n,k,&c); assert(aux);
    printf("workspace,%d,%d,0,%zu\n",n,k,main_bytes+allocated_bytes);
    cqreg_ipm_free(aux);
    start=ctools_timer_seconds();
    int rc=cqreg_fn_solve(s,y,x,.5,b); assert(rc>0 && s->converged);
    double elapsed=ctools_timer_seconds()-start, obj=0;
    for(int i=0;i<n;i++) obj+=.5*fabs(s->r_primal[i]);
    printf("solve,%d,%d,%.9f,%.17g,%d",n,k,elapsed,obj,rc);
    for(int j=0;j<k;j++) printf(",%.17g",b[j]);
    puts(""); cqreg_ipm_free(s);
    /* A modest tied-residual case exposes quadratic selection without a big fit. */
    int tied_n=n<12000?n:12000;
    for(int i=0;i<tied_n;i++) y[i]=i%3-1;
    cqreg_sparsity_state *sp=cqreg_sparsity_create(tied_n,.5,CQREG_BW_HSHEATHER);
    assert(sp); start=ctools_timer_seconds();
    double sparsity=cqreg_estimate_sparsity(sp,y);
    printf("tied_sparsity,%d,%d,%.9f,%.17g\n",tied_n,k,ctools_timer_seconds()-start,sparsity);
    ctools_aligned_free(sp->sorted_resid); free(sp); free(x); free(y); free(b);
}
int main(int argc, char **argv) {
    if (argc>1) { benchmark(atoi(argv[1]),atoi(argv[2])); return 0; }
    selection_checks(); gemv_checks(); solver_checks(); puts("cqreg native checks passed");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_regress)
/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls, scalar_at, scalar_calls, matrix_at, matrix_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int clustered = 0, weighted = 0, mode = 0;
static ST_retcode message(char *s) { return 0; }
static ST_int one(void) { return 1; }
static ST_int obs(void) { return 12; }
static ST_int nvars(void) { return 4 + (mode == 2) + clustered + !!weighted + (mode == 2) + 2*(mode == 3); }
static ST_retcode scalar(char *k, double *v) {
    *v = 0;
    if (!strcmp(k, "__creghdfe_K"))
        *v = 2;
    else if (!strcmp(k, "__creghdfe_G"))
        *v = mode == 2 ? 2 : 1;
    else if (!strcmp(k, "__creghdfe_maxiter"))
        *v = 10000;
    else if (!strcmp(k, "__creghdfe_tolerance"))
        *v = 1e-10;
    else if (!strcmp(k, "__creghdfe_standardize"))
        *v = 1;
    else if (!strcmp(k, "__creghdfe_vce_type"))
        *v = clustered ? 2 : 0;
    else if (!strcmp(k, "__creghdfe_has_weights"))
        *v = !!weighted;
    else if (!strcmp(k, "__creghdfe_weight_type"))
        *v = weighted == 2 ? 2 : 1;
    else if (!strcmp(k, "__creghdfe_sample_idx"))
        *v = nvars();
    else if (!strcmp(k, "__creghdfe_compute_dof") || !strcmp(k, "__creghdfe_compute_groupvar"))
        *v = mode == 2;
    else if (!strcmp(k, "__creghdfe_groupvar_idx"))
        *v = mode == 2 ? nvars() - 1 : 0;
    else if (!strcmp(k, "__creghdfe_savefe") || !strcmp(k, "__creghdfe_compute_resid"))
        *v = mode==3;
    else if (!strcmp(k, "__creghdfe_savefe_idx"))
        *v = mode==3?nvars()-1:0;
    else if (!strcmp(k, "__creghdfe_resid_var_idx"))
        *v = mode==3?nvars()-2:0;
    return 0;
}
static ST_retcode read_value(ST_int v, ST_int row, double *x) {
    int i=row-1;*x=v==1?3+2*i+sin(i):i;return 0;
}
static ST_retcode save(char *k, double v) { return ++scalar_calls == scalar_at ? 459 : 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    return ++store_calls == store_at ? 459 : 0;
}
static ST_retcode mat(char *k, ST_int i, ST_int j, double x) {
    return ++matrix_calls == matrix_at ? 459 : 0;
}
static ST_boolean missing(double x) { return x >= 8.98846567431158e307; }
int ctools_get_max_threads(void) { return 1; }
int ctools_get_max_threads_used(void) { return 1; }
int ctools_get_num_procs(void) { return 1; }
double ctools_timer_seconds(void) { return 0; }
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void stata_data_free(stata_data *d) {
    for (size_t j = 0; j < d->nvars; j++)
        free(d->vars[j].data.dbl);
    free(d->vars);
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d) {
    stata_data_free(&d->data);
    free(d->obs_map);
    d->obs_map = NULL;
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = 12;
    d->obs_map = malloc(12 * sizeof(perm_idx_t));
    if (!d->obs_map)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < 12; i++)
        d->obs_map[i] = i + 1;
    for (size_t j = 0; j < n; j++) {
        d->data.vars[j].nobs = 12;
        d->data.vars[j].type = STATA_TYPE_DOUBLE;
        double *p = malloc(12 * sizeof(double));
        d->data.vars[j].data.dbl = p;
        if (!p)
            return STATA_ERR_MEMORY;
        for (int i = 0; i < 12; i++)
            p[i] = j == 0                              ? (mode == 1 ? 7 : 3 + 2 * i + sin(i))
                   : j == 1                            ? i
                   : j == 2                            ? i % 3
                   : j == 3 && mode == 2               ? i % 2
                   : j == 3 + (mode == 2) && clustered ? i % 4
                                                       : (weighted == 2 ? 300000000 : 1);
    }
    return STATA_OK;
}
#include "creghdfe/creghdfe_hdfe.c"
#include "creghdfe/creghdfe_regress.c"
#include "creghdfe/creghdfe_result.c"
#include "creghdfe/creghdfe_input.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "creghdfe/creghdfe_vce.c"
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "ctools_spi.c"
int main(int argc, char **argv) {
    store_at = argc > 5 ? atoi(argv[5]) : 0;
    scalar_at = argc > 6 ? atoi(argv[6]) : 0;
    matrix_at = argc > 7 ? atoi(argv[7]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    clustered = argc > 2 ? atoi(argv[2]) : 1;
    weighted = argc > 3 ? atoi(argv[3]) : 1;
    mode = argc > 4 ? atoi(argv[4]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.nvars = nvars;
    mock.scalaruse = scalar;
    mock.scalsave = save;
    mock.spouterr = message;
    mock.spoutsml = message;
    mock.safevdata = read_value;
    mock.vdata = read_value;
    mock.safestore = store;
    mock.matstore = mat;
    mock.safematstore = mat;
    mock.ismissing = missing;
    mock.missval = 8.98846567431158e307;
    int rc = do_full_regression(0, NULL);
    cleanup_state();
    int first_rc = rc, first_calls = calls;
    assert(alive == 0);
    fail_at = 0;
    store_at = scalar_at = matrix_at = 0;
    calls = 0;
    assert(do_full_regression(0, NULL) == 0);
    assert(alive == 0);
    rc = first_rc;
    calls = first_calls;
    printf("rc=%d allocations=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_ppml)
/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls, scalar_at, scalar_calls, matrix_at, matrix_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int clustered = 0, weighted = 0, mode = 0;
static ST_retcode message(char *s) { return 0; }
static ST_int one(void) { return 1; }
static ST_int obs(void) { return 12; }
static ST_int nvars(void) { return 4 + clustered + weighted; }
static ST_retcode scalar(char *k, double *v) {
    *v = 0;
    if (!strcmp(k, "__cppmlhdfe_K"))
        *v = 2;
    else if (!strcmp(k, "__cppmlhdfe_G"))
        *v = 1;
    else if (!strcmp(k, "__cppmlhdfe_maxiter"))
        *v = 10000;
    else if (!strcmp(k, "__cppmlhdfe_tolerance"))
        *v = 1e-10;
    else if (!strcmp(k, "__cppmlhdfe_standardize"))
        *v = 1;
    else if (!strcmp(k, "__cppmlhdfe_vce_type"))
        *v = clustered ? 2 : 0;
    else if (!strcmp(k, "__cppmlhdfe_cluster_terms"))
        *v = clustered;
    else if (!strcmp(k, "__cppmlhdfe_irls_maxiter"))
        *v = 100;
    else if (!strcmp(k, "__cppmlhdfe_irls_tol"))
        *v = 1e-8;
    else if (!strcmp(k, "__cppmlhdfe_has_weights"))
        *v = weighted;
    else if (!strcmp(k, "__cppmlhdfe_weight_type"))
        *v = 1;
    return 0;
}
static ST_retcode save(char *k, double v) { return ++scalar_calls == scalar_at ? 459 : 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    return ++store_calls == store_at ? 459 : 0;
}
static ST_retcode mat(char *k, ST_int i, ST_int j, double x) {
    return ++matrix_calls == matrix_at ? 459 : 0;
}
static ST_boolean missing(double x) { return x >= 8.98846567431158e307; }
int ctools_get_max_threads(void) { return 1; }
int ctools_get_max_threads_used(void) { return 1; }
int ctools_get_num_procs(void) { return 1; }
double ctools_timer_seconds(void) { return 0; }
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void stata_data_free(stata_data *d) {
    for (size_t j = 0; j < d->nvars; j++)
        free(d->vars[j].data.dbl);
    free(d->vars);
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d) {
    stata_data_free(&d->data);
    free(d->obs_map);
    d->obs_map = NULL;
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = 12;
    d->obs_map = malloc(12 * sizeof(perm_idx_t));
    if (!d->obs_map)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < 12; i++)
        d->obs_map[i] = i + 1;
    for (size_t j = 0; j < n; j++) {
        d->data.vars[j].nobs = 12;
        d->data.vars[j].type = STATA_TYPE_DOUBLE;
        double *p = malloc(12 * sizeof(double));
        d->data.vars[j].data.dbl = p;
        if (!p)
            return STATA_ERR_MEMORY;
        for (int i = 0; i < 12; i++)
            p[i] = j == 0                ? (mode == 1 ? 7 : 3 + 2 * i + sin(i))
                   : j == 1              ? i
                   : j == 2              ? i % 3
                   : j == 3 && clustered ? i % 4
                                         : 1;
    }
    return STATA_OK;
}
#include "cppmlhdfe/cppmlhdfe_irls.c"
#include "cppmlhdfe/cppmlhdfe_result.c"
#include "cppmlhdfe/cppmlhdfe_input.c"
#include "cppmlhdfe/cppmlhdfe_options.c"
#include "cppmlhdfe/cppmlhdfe_solvers.c"
#include "cppmlhdfe/cppmlhdfe_simplex.c"
#include "cppmlhdfe/cppmlhdfe_relu.c"
#include "cppmlhdfe/cppmlhdfe_moments.c"
#include "cppmlhdfe/cppmlhdfe_separation.c"
#include "creghdfe/creghdfe_hdfe.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "creghdfe/creghdfe_vce.c"
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "ctools_spi.c"
int main(int argc, char **argv) {
    store_at = argc > 5 ? atoi(argv[5]) : 0;
    scalar_at = argc > 6 ? atoi(argv[6]) : 0;
    matrix_at = argc > 7 ? atoi(argv[7]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    clustered = argc > 2 ? atoi(argv[2]) : 1;
    weighted = argc > 3 ? atoi(argv[3]) : 1;
    mode = argc > 4 ? atoi(argv[4]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.nvars = nvars;
    mock.scalaruse = scalar;
    mock.scalsave = save;
    mock.spouterr = message;
    mock.spoutsml = message;
    mock.safestore = store;
    mock.matstore = mat;
    mock.safematstore = mat;
    mock.ismissing = missing;
    mock.missval = 8.98846567431158e307;
    int rc = do_ppml_regression(0, NULL);
    cleanup_ppml_state();
    int first_rc = rc, first_calls = calls;
    assert(alive == 0);
    fail_at = 0;
    store_at = scalar_at = matrix_at = 0;
    calls = 0;
    assert(do_ppml_regression(0, NULL) == 0);
    assert(alive == 0);
    rc = first_rc;
    calls = first_calls;
    printf("rc=%d allocations=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_iv)
/* Full IV command: injected allocations must fail cleanly or preserve estimates. */
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <math.h>

static ST_plugin host;
ST_plugin *_stata_ = &host;
enum
{
    OBS = 40
};
static int weighted, two_fe, singleton;
static int save_at, save_calls; /* nonzero: the save_at-th scalar save fails */
static const char *missing_scalar; /* reported as not set, its destination untouched */
static double coefficient, variance, sample_size;
static ST_retcode message(char *s)
{
    (void)s;
    return 0;
}
static ST_int one(void)
{
    return 1;
}
static ST_int obs(void)
{
    return OBS;
}
static ST_retcode scalar(char *key, double *value)
{
    if (missing_scalar && !strcmp(key, missing_scalar))
        return 111;
    *value = 0;
    if (!strcmp(key, "__civreghdfe_K_endog") || !strcmp(key, "__civreghdfe_K_iv"))
        *value = 1;
    else if (!strcmp(key, "__civreghdfe_G"))
        *value = 1 + two_fe;
    else if (!strcmp(key, "__civreghdfe_has_weights"))
        *value = weighted;
    else if (!strcmp(key, "__civreghdfe_weight_type"))
        *value = 1;
    else if (!strcmp(key, "__civreghdfe_maxiter"))
        *value = 10000;
    else if (!strcmp(key, "__civreghdfe_tolerance"))
        *value = 1e-10;
    return 0;
}
static ST_retcode save(char *key, double value)
{
    if (save_at && ++save_calls == save_at)
        return 459;
    if (!strcmp(key, "__civreghdfe_N"))
        sample_size = value;
    return 0;
}
static ST_retcode matrix(char *key, ST_int row, ST_int col, double value)
{
    (void)row;
    (void)col;
    if (!strcmp(key, "__civreghdfe_b"))
        coefficient = value;
    if (!strcmp(key, "__civreghdfe_V"))
        variance = value;
    return 0;
}
static ST_boolean missing(double x)
{
    return x >= 8.98846567431158e307;
}
int ctools_get_max_threads(void)
{
    return 1;
}
int ctools_get_max_threads_used(void)
{
    return 1;
}
int ctools_get_num_procs(void)
{
    return 1;
}
double ctools_timer_seconds(void)
{
    return 0;
}
void ctools_filtered_data_init(ctools_filtered_data *data)
{
    memset(data, 0, sizeof(*data));
}
void ctools_filtered_data_free(ctools_filtered_data *data)
{
    for (size_t j = 0; j < data->data.nvars; j++)
        free(data->data.vars[j].data.dbl);
    free(data->data.vars);
    free(data->obs_map);
    memset(data, 0, sizeof(*data));
}
stata_retcode ctools_data_load(ctools_filtered_data *data, int *indices, size_t n, size_t start,
                               size_t end, int flags)
{
    (void)indices;
    (void)start;
    (void)end;
    (void)flags;
    data->data.vars = calloc(n, sizeof(stata_variable));
    if (!data->data.vars)
        return STATA_ERR_MEMORY;
    data->data.nvars = n;
    data->data.nobs = OBS;
    data->obs_map = malloc(OBS * sizeof(perm_idx_t));
    if (!data->obs_map)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < OBS; i++)
        data->obs_map[i] = i + 1;
    for (size_t j = 0; j < n; j++)
    {
        double *p = malloc(OBS * sizeof(double));
        data->data.vars[j].data.dbl = p;
        if (!p)
            return STATA_ERR_MEMORY;
        data->data.vars[j].nobs = OBS;
        data->data.vars[j].type = STATA_TYPE_DOUBLE;
        for (int i = 0; i < OBS; i++)
        {
            double z = sin(0.7 * i) + 0.1 * i;
            double x = z + cos(1.3 * i);
            p[i] = j == 0             ? 2 * x + 3 * (i % 4) + sin(2.1 * i)
                   : j == 1           ? x
                   : j == 2           ? z
                   : j == 3           ? (singleton && (i == 3 || i == 18) ? 100 + i : i % 4)
                   : j == 4 && two_fe ? (i / 3) % 3
                                      : 1 + (i % 3);
        }
    }
    return STATA_OK;
}

#include "civreghdfe/civreghdfe_impl.c"
#include "civreghdfe/civreghdfe_result.c"
#include "civreghdfe/civreghdfe_input.c"
#include "civreghdfe/civreghdfe_estimate.c"
#include "civreghdfe/civreghdfe_matrix.c"
#include "civreghdfe/civreghdfe_omega.c"
#include "civreghdfe/civreghdfe_tests.c"
#include "civreghdfe/civreghdfe_vce.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "ctools_spi.c"

int main(int argc, char **argv)
{
    int failure = argc > 1 ? atoi(argv[1]) : 0;
    weighted = argc > 2 ? atoi(argv[2]) : 0;
    two_fe = argc > 3 ? atoi(argv[3]) : 0;
    singleton = argc > 4 ? atoi(argv[4]) : 0;
    int save_failure = argc > 5 ? atoi(argv[5]) : 0;
    const char *missing = argc > 6 ? argv[6] : NULL;
    host.nobs1 = one;
    host.nobs2 = obs;
    host.nobs = obs;
    host.scalaruse = scalar;
    host.scalsave = save;
    host.spouterr = message;
    host.spoutsml = message;
    host.matstore = matrix;
    host.safematstore = matrix;
    host.ismissing = missing;
    host.missval = 8.98846567431158e307;
    assert(do_iv_regression() == 0);
    assert(alive == 0 && iv_state == NULL);
    double expected_b = coefficient, expected_v = variance, expected_n = sample_size;
    int allocations = calls;
    fail_at = failure;
    calls = 0;
    save_at = save_failure;
    save_calls = 0;
    missing_scalar = missing;
    coefficient = variance = sample_size = NAN;
    int rc = do_iv_regression();
    printf("failure=%d rc=%d allocations=%d alive=%d b=%.17g expected=%.17g\n", failure, rc,
           allocations, alive, coefficient, expected_b);
    leftovers();
    fflush(stdout);
    assert(alive == 0 && iv_state == NULL);
    /* A failed scalar save is the command's status, never ignored; a
       configuration scalar that is not set is an error, not the value of the
       scalar read before it */
    if (missing_scalar)
        assert(rc == 198);
    else if (save_at)
        assert(rc == 459);
    else if (rc == 0)
    {
        assert(fabs(coefficient - expected_b) < 1e-10);
        assert(fabs(variance - expected_v) < 1e-10);
        assert(sample_size == expected_n);
    }
    else
        assert(rc == 920 || rc == 430);
    /* A failed command must not contaminate the next command in the same host. */
    fail_at = 0;
    save_at = 0;
    missing_scalar = NULL;
    assert(do_iv_regression() == 0);
    assert(alive == 0 && iv_state == NULL);
    assert(fabs(coefficient - expected_b) < 1e-10);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_bins_hdfe)
/* Compare batch projections with an independent dummy-variable least-squares solve. */
#include "alloc.h"
#include "stplugin.h"
#include <math.h>
#include "cbinscatter/cbinscatter_resid.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "ctools_hdfe_utils.c"
#include "ctools_ols.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
static ST_retcode message(char *s)
{
    (void)s;
    return 0;
}
enum
{
    OBS = 48,
    COLS = 5,
    DUMMIES = 6
};

static void reference(const double *input, const int *fe, const double *weights, double *out)
{
    double gram[DUMMIES][DUMMIES + 1] = {{0}};
    for (int i = 0; i < OBS; i++)
    {
        double row[DUMMIES] = {
            1, fe[i] == 2, fe[i] == 3, fe[OBS + i] == 2, fe[OBS + i] == 3, fe[OBS + i] == 4};
        double w = weights ? weights[i] : 1;
        for (int j = 0; j < DUMMIES; j++)
        {
            for (int k = 0; k < DUMMIES; k++)
                gram[j][k] += w * row[j] * row[k];
            gram[j][DUMMIES] += w * row[j] * input[i];
        }
    }
    for (int j = 0; j < DUMMIES; j++)
    {
        double pivot = gram[j][j];
        assert(fabs(pivot) > 1e-6);
        for (int k = j; k <= DUMMIES; k++)
            gram[j][k] /= pivot;
        for (int l = 0; l < DUMMIES; l++)
        {
            if (l == j)
                continue;
            double factor = gram[l][j];
            for (int k = j; k <= DUMMIES; k++)
                gram[l][k] -= factor * gram[j][k];
        }
    }
    for (int i = 0; i < OBS; i++)
    {
        out[i] = input[i] - gram[0][DUMMIES];
        if (fe[i] > 1)
            out[i] -= gram[fe[i] - 1][DUMMIES];
        if (fe[OBS + i] > 1)
            out[i] -= gram[1 + fe[OBS + i]][DUMMIES];
    }
}

int main(int argc, char **argv)
{
    host.spoutsml = host.spouterr = message;
    int failure = argc > 1 ? atoi(argv[1]) : 0;
    int weighted = argc > 2 ? atoi(argv[2]) : 0;
    int first = argc > 3 ? atoi(argv[3]) : 0;
    int maxiter = argc > 4 ? atoi(argv[4]) : 10000;
    int fe[2 * OBS];
    double weights[OBS], values[COLS][OBS], expected[COLS][OBS], *vars[COLS];
    for (int i = 0; i < OBS; i++)
    {
        fe[i] = 1 + (i / 16) % 3;
        fe[OBS + i] = 1 + (i / 12) % 4;
        weights[i] = 1 + i % 5;
        values[0][i] = 0;         /* The old stopping rule immediately accepts this column. */
        values[1][i] = i % 7 < 3; /* Bin indicator with nonzero FE group means. */
        values[2][i] = sin(0.7 * i) + i * 0.05;
        values[3][i] = cos(1.1 * i);
        values[4][i] = 2 * fe[i] - 3 * fe[OBS + i]; /* Fully absorbed column. */
    }
    const double *w = weighted ? weights : NULL;
    /* Also test a nonzero first column already orthogonal to every FE. */
    if (first == 1)
    {
        reference(values[3], fe, w, expected[3]);
        memcpy(values[0], expected[3], sizeof(values[0]));
    }
    for (int v = 0; v < COLS; v++)
    {
        reference(values[v], fe, w, expected[v]);
        vars[v] = values[(v + (first == 2)) % COLS];
    }
    fail_at = failure;
    int dropped = -1;
    int rc =
        hdfe_residualize_batch(vars, COLS, fe, OBS, 2, w, weighted, maxiter, 1e-12, 0, &dropped);
    double error = 0, group_mean = 0;
    for (int v = 0; v < COLS; v++)
    {
        for (int i = 0; i < OBS; i++)
            error = fmax(error, fabs(values[v][i] - expected[v][i]));
        for (int g = 0; g < 2; g++)
        {
            for (int level = 1; level <= 4; level++)
            {
                double sum = 0, count = 0;
                for (int i = 0; i < OBS; i++)
                {
                    if (fe[g * OBS + i] != level)
                        continue;
                    double wi = w ? w[i] : 1;
                    sum += wi * values[v][i];
                    count += wi;
                }
                if (count)
                    group_mean = fmax(group_mean, fabs(sum / count));
            }
        }
    }
    printf("failure=%d rc=%d allocations=%d alive=%d error=%.9g group_mean=%.9g\n", failure, rc,
           calls, alive, error, group_mean);
    leftovers();
    fflush(stdout);
    assert(alive == 0);
    if (failure)
        assert(rc == CBINSCATTER_ERR_MEMORY);
    else if (maxiter <= 1)
        assert(rc == CBINSCATTER_ERR_HDFE);
    else
    {
        assert(rc == 0 && dropped == 0);
        assert(error < 1e-9 && group_mean < 1e-10);
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_5)
#include "ctools_matrix.c"
int main(void) {
 double X[]={1,2,3,4},D[]={1.0/30},e[]={1,-1,2,1},w[]={1,2,1,2},V[1]; int ids[]={0,0,1,1};
 ctools_vce_data d={.X_eff=X,.D=D,.resid=e,.weights=w,.N=4,.K=1,.weight_type=1,.normalize_weights=1};
 for(int cluster=0;cluster<2;cluster++) {
  calls=0;fail_at=0;V[0]=123;
  assert((cluster?ctools_vce_cluster(&d,ids,2,1,V):ctools_vce_robust(&d,1,V))==0);
  int allocations=calls;assert(V[0]>0 && V[0]!=123 && live_count==0);
  for(int f=1;f<=allocations;f++) {calls=0;fail_at=f;V[0]=123;
   assert((cluster?ctools_vce_cluster(&d,ids,2,1,V):ctools_vce_robust(&d,1,V))==920);
   assert(V[0]==123 && live_count==0);
  }
 }
 fail_at=0; memset(e,0,sizeof(e));assert(ctools_vce_robust(&d,1,V)==0 && V[0]==0);
 V[0]=123;assert(ctools_vce_cluster(&d,ids,1,1,V)==498 && V[0]==123);
 e[0]=1e308;V[0]=123;assert(ctools_vce_robust(&d,1,V)==498 && V[0]==123);
 assert(live_count==0);return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_9)
#include "cbinscatter/cbinscatter_resid.c"
int main(void) {
 double gram[9]={0},rhs[3]={0},beta[3];
 for(int i=0;i<5;i++){double x[]={1,1e-9*i,2e-9*i},y=3+2*i;
  for(int j=0;j<3;j++){rhs[j]+=x[j]*y;for(int k=0;k<3;k++)gram[j*3+k]+=x[j]*x[k];}}
 assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==0);
 int allocations=calls;assert(live_count==0);
 for(int i=0;i<5;i++)assert(fabs(beta[0]+1e-9*i*beta[1]+2e-9*i*beta[2]-(3+2*i))<1e-10);
 for(int f=1;f<=allocations;f++) {calls=0;fail_at=f;beta[0]=beta[1]=beta[2]=123;
  assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==920);
  assert(beta[0]==123 && beta[1]==123 && beta[2]==123 && live_count==0);}
 fail_at=0;gram[0]=-1;
 assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==199);
 assert(live_count==0);return 0;
}
#endif
