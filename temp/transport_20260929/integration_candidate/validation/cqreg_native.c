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
