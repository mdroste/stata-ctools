/* Deterministic Break injection: no SPI on workers, prompt cancellation of
 * every projection family, and a fresh successful command after cancellation. */
#undef NDEBUG
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include "ctools_cancel.h"
#include "ctools_hdfe_utils.c"
#include "ctools_ols.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_indiv.c"
#include "hdfe/hdfe_fe_recover.c"
#include "hdfe/hdfe_boot.c"
#include "cppmlhdfe/cppmlhdfe_solvers.c"
#include "cppmlhdfe/cppmlhdfe_moments.c"
#include "cppmlhdfe/cppmlhdfe_simplex.c"
#include "cppmlhdfe/cppmlhdfe_relu.c"

/* Loading is outside this kernel fixture; accidental SPI use must fail. */
ST_plugin *_stata_;
void ctools_filtered_data_init(ctools_filtered_data *fd) { (void)fd; abort(); }
void ctools_filtered_data_free(ctools_filtered_data *fd) { (void)fd; abort(); }
stata_retcode ctools_data_load(ctools_filtered_data *fd, int *v, size_t k,
                              size_t first, size_t last, int flags)
{
    (void)fd; (void)v; (void)k; (void)first; (void)last; (void)flags;
    abort();
}

static pthread_t host;
static int polls, cancel_after;
static int poll_host(void)
{
    assert(pthread_equal(pthread_self(), host));
    return ++polls >= cancel_after;
}
static void begin(int after)
{
    host = pthread_self(); polls = 0; cancel_after = after;
    ctools_cancel_begin(poll_host);
}
static void reset_data(double *data, int n, int k)
{
    for (int j = 0; j < k; j++)
        for (int i = 0; i < n; i++)
            data[(size_t)j*n+i] = sin(i*.37+j) + .02*(i%17);
}
static void setup(HDFE_State *s, int weighted, int threads)
{
    memset(s, 0, sizeof(*s));
    s->N = 512; s->G = 2; s->num_threads = threads;
    s->maxiter = INT_MAX-1; s->tolerance = 1e-30;
    s->has_weights = weighted;
    if (weighted) {
        s->weights = calloc(s->N, sizeof(double));
        for (int i=0; i<s->N; i++) s->weights[i] = 1 + i%5;
    }
    s->factors = calloc(s->G, sizeof(FE_Factor));
    for (int g=0; g<s->G; g++) {
        FE_Factor *f = &s->factors[g];
        f->num_levels = 32; f->has_intercept = 1;
        f->levels = calloc(s->N, sizeof(int));
        f->counts = calloc(32, sizeof(double));
        f->weighted_counts = calloc(32, sizeof(double));
        for (int i=0; i<s->N; i++) {
            int l = (i/16 + (g && i%16==0))%32;
            f->levels[i]=1+l; f->counts[l]++;
            f->weighted_counts[l] += weighted ? s->weights[i] : 1;
        }
    }
    assert(!ctools_hdfe_alloc_buffers(s, 1, 8));
}
static void cancelled_result(HDFE_SolveResult r)
{
    assert(r.status == 1 && !r.converged);
    assert(polls == cancel_after);
    assert(ctools_cancel_requested());
    ctools_cancel_end();
    assert(!ctools_cancel_requested());
}

static int live_workspaces;
static size_t cancel_fit_bytes(const void *m) { (void)m; return 8; }
static void *cancel_fit_new(const void *m)
{
    (void)m; live_workspaces++; return malloc(8);
}
static void cancel_fit_free(void *p) { free(p); live_workspaces--; }
static ST_retcode cancel_fit(const void *m, void *ws, const int *mult, hdfe_boot_rep *out)
{
    (void)m; (void)ws; (void)mult; (void)out;
    while (!ctools_cancel_requested()) { }
    return 1;
}
static void bootstrap_cancel(int threads)
{
    hdfe_boot_design d={4,1,NULL,NULL,17};
    hdfe_boot_model_ops ops={cancel_fit_bytes,cancel_fit_new,cancel_fit_free,cancel_fit,NULL};
    hdfe_boot_acc acc;
    assert(!hdfe_boot_acc_init(&acc,1,0,0,8));
    begin(5);
    assert(hdfe_boot_run(&d,NULL,&ops,NULL,0,8,threads,1<<20,&acc)==1);
    assert(!live_workspaces && acc.R==0 && polls==5);
    ctools_cancel_end();
    hdfe_boot_acc_free(&acc);
}

#ifdef _OPENMP
/* The host has no work; a worker cannot finish until the host polls while
 * waiting. A barrier without polling deadlocks (the runner has a timeout). */
static void waiting_host(void)
{
    begin(4);
    ctools_cancel_team team = {0};
    int size = 0;
#pragma omp parallel num_threads(2) shared(size)
    {
        if (omp_get_thread_num() == 0) size = omp_get_num_threads();
        else while (!ctools_cancel_requested()) { }
        ctools_cancel_wait(&team);
    }
    assert(size == 2 && polls == 4);
    ctools_cancel_end();
    /* Nested teams must neither poll from another team's master nor reuse
     * its arrival counter. */
    omp_set_max_active_levels(2);
    begin(4);
    ctools_cancel_team outer = {0};
#pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() != 0) {
            ctools_cancel_team inner = {0};
#pragma omp parallel num_threads(2)
            {
                while (!ctools_cancel_requested()) { }
                ctools_cancel_wait(&inner);
            }
        }
        ctools_cancel_wait(&outer);
    }
    assert(polls == 4);
    ctools_cancel_end();
}
#endif

int main(void)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
    waiting_host();
#endif
    for (int threads=1; threads<=4; threads*=4)
    for (int weighted=0; weighted<=1; weighted++) {
        HDFE_State s; setup(&s, weighted, threads);
        double data[512*8];
        for (int transform=0; transform<4; transform++)
        for (int accel=0; accel<5; accel++) {
            reset_data(data,s.N,8); begin(4);
            cancelled_result(partial_out_columns_options(&s,data,s.N,8,transform,accel));
        }
        for (int g=0; g<=1; g++) {
            s.G=g; reset_data(data,s.N,8); begin(1);
            cancelled_result(partial_out_columns(&s,data,s.N,8,threads));
        }
        s.G=2;
        /* Individual/group LSMR and fixed-effect recovery use the same flag. */
        int ptr[513], ind[512];
        for (int i=0; i<=512; i++) ptr[i]=i;
        for (int i=0; i<512; i++) ind[i]=(i/8)%32;
        HDFE_Indiv indiv={1,512,32,ptr,ind};
        reset_data(data,s.N,8); begin(4);
        assert(hdfe_indiv_partial_out(&s,&indiv,data,s.N,8,1e-30,INT_MAX-1,NULL)==1);
        assert(polls==4); ctools_cancel_end();
        double a0[32]={0},a1[32]={0}; double *alpha[]={a0,a1};
        reset_data(data,s.N,1); begin(2);
        assert(hdfe_fe_recover(&s,data,alpha)==1);
        ctools_cancel_end();
        if (weighted) {
            PPML_Options options=PPML_OPTIONS_DEFAULT;
            reset_data(data,s.N,8); begin(4);
            cancelled_result(ppml_reference_partial(&s,data,s.N,8,threads));
            for (int method=0; method<5; method++) {
                options.acceleration=method;
                reset_data(data,s.N,8); begin(4);
                cancelled_result(ppml_map_options(&options,&s,data,s.N,8));
            }
            options.btol=1e-30;
            reset_data(data,s.N,8); begin(4);
            cancelled_result(ppml_lsmr(&options,&s,data,s.N,8));
        }
        s.maxiter=10000; s.tolerance=1e-8;
        reset_data(data,s.N,8); begin(INT_MAX);
        HDFE_SolveResult ok=partial_out_columns(&s,data,s.N,8,threads);
        assert(!ok.status && ok.converged && polls>0);
        ctools_cancel_end();
        free(s.weights); s.weights=NULL; ctools_hdfe_state_cleanup(&s);
    }
    bootstrap_cancel(1);
    bootstrap_cancel(4);
    /* Analytic covariance transform with cancellation in the constant's
     * variance: all entries and symmetry, with no reference estimation. */
    double V[]={1,1,1, 1,1e16,1e16, 1,1e16,1e16+4};
    const double means[]={1,1};
    ppml_uncenter_vce(V,means,2);
    assert(V[2]==-1 && V[5]==-1 && V[6]==-1 && V[7]==-1 && V[8]==3);
    puts("PASS HDFE cancellation families, reset, and PPML covariance transform");
    return 0;
}
