/*
 * cpplmhdfe_irls.c
 *
 * IRLS (Iteratively Reweighted Least Squares) loop for PPML estimation
 * with high-dimensional fixed effects.
 *
 * Algorithm:
 *   1. Load data, set up FE, remove singletons
 *   2. Detect separation (FE-level: all y=0 in a group)
 *   3. Initialize beta from OLS of log(max(y,1)) on partialled X
 *   4. IRLS loop:
 *      a. Compute IRLS weights: w_irls = mu * w_user
 *      b. Compute working depvar: z = eta + (y - mu)/mu
 *      c. Partial out [z, X] via CG solver with IRLS weights
 *      d. Solve WLS: beta = (X~'X~)^{-1} X~'z~
 *      e. Update eta, mu; check deviance convergence
 *   5. VCE: sandwich with score residuals (y - mu)
 *   6. Store results
 *
 * Part of the ctools Stata plugin suite
 */

/* Resolve output indices in the plugin varlist, including hidden Stata variables. */
#define SD_SAFEMODE
#include <float.h>
#include <limits.h>
#include "cpplmhdfe_irls.h"
#include "cpplmhdfe_separation.h"
#include "../creghdfe/creghdfe_utils.h"
#include "../creghdfe/creghdfe_solver.h"
#include "../ctools_matrix.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_types.h"
#include "../ctools_spi.h"
#include "../ctools_ols.h"

/* Forward declarations for helpers defined at bottom of file */
static void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N);
static ST_double ppml_compute_deviance(const ST_double *y, const ST_double *mu,
                                        const ST_double *w_user, ST_int N);
static ST_double ppml_compute_loglik(const ST_double *y, const ST_double *mu,
                                      const ST_double *w_user, ST_int N);

typedef struct {
    int standardize, min_ok, exact_solver, exact_partial, heuristic, remove_collinear;
    int step_halving, max_halving, guess, guess_idx, check_mu, check_simplex;
    int relu_maxiter, relu_strict, relu_accelerate, simplex_maxiter;
    int tag_idx, certificate_idx;
    double start_tol, halving_memory, mu_tol, relu_tol, relu_zero_tol, simplex_tol;
    int acceleration,transform,pool;
    double btol,conlim;
} PPML_Options;
static PPML_Options ppml_options = {1,1,0,1,1,1,0,2,0,0,0,0,100,0,0,1000,0,0,
                                    1e-4,.9,1e-6,1e-4,1e-8,1e-12,0,0,0,1e-8,1e8};
static double ppml_option(const char *name,double fallback)
{
    char key[64];double value;
    snprintf(key,sizeof(key),"__cpplmhdfe_%s",name);
    return SF_scal_use(key,&value)==0 ? value : fallback;
}
static void ppml_read_options(void)
{
#define OPT(field,name,value) ppml_options.field=ppml_option(name,value)
    OPT(standardize,"standardize_data",1); OPT(min_ok,"min_ok",1);
    OPT(exact_solver,"use_exact_solver",0); OPT(exact_partial,"use_exact_partial",1);
    OPT(heuristic,"use_heuristic_tol",1); OPT(remove_collinear,"remove_collinear",1);
    OPT(step_halving,"use_step_halving",0); OPT(max_halving,"max_step_halving",2);
    OPT(guess,"guess",0); OPT(guess_idx,"guess_idx",0); OPT(check_mu,"check_mu",0);
    OPT(check_simplex,"check_simplex",0); OPT(relu_maxiter,"relu_maxiter",100);
    OPT(relu_strict,"relu_strict",0); OPT(relu_accelerate,"relu_accelerate",0);
    OPT(simplex_maxiter,"simplex_maxiter",1000); OPT(tag_idx,"tag_idx",0);
    OPT(certificate_idx,"certificate_idx",0); OPT(start_tol,"start_inner_tol",1e-4);
    OPT(halving_memory,"step_halving_memory",.9); OPT(mu_tol,"mu_tol",1e-6);
    OPT(relu_tol,"relu_tol",1e-4); OPT(relu_zero_tol,"relu_zero_tol",1e-8);
    OPT(simplex_tol,"simplex_tol",1e-12);
    OPT(acceleration,"acceleration",0); OPT(transform,"transform",0); OPT(pool,"pool",0);
    OPT(btol,"btol",1e-8); OPT(conlim,"conlim",1e8);
#undef OPT
}

/* Reference CG stopping controls are essential in the separation problem:
 * a near-singular projection can otherwise change the retained sample. */
static double ppml_quad_dot(const double *x, const double *y, const double *w, int n)
{
    dd_real sum = {0,0};
    for (int i=0;i<n;i++) {
        double xy=x[i]*y[i];
        dd_real prod=two_prod(xy,w ? w[i] : 1.0);
        sum=dd_add_d(dd_add_d(sum,prod.hi),prod.lo);
    }
    return sum.hi+sum.lo;
}

/* Rank-revealing weighted QR, with a second orthogonalization pass. */
static int ppml_qrsolve(const double *X,const double *y,const double *w,int N,int K,double *b,double *resid)
{
    double *q=malloc((size_t)N*(K?K:1)*sizeof(double));
    double *R=calloc((size_t)(K?K*K:1),sizeof(double));
    double *qy=calloc(K?K:1,sizeof(double));
    int *order=malloc((K?K:1)*sizeof(int));
    if(!q||!R||!qy||!order){free(q);free(R);free(qy);free(order);return 920;}
    double maximum=0;int rank=0;
    for(int k=0;k<K;k++) {
        order[k]=k;if(b)b[k]=0;
        for(int i=0;i<N;i++)q[(size_t)k*N+i]=X[(size_t)k*N+i]*sqrt(w?w[i]:1);
        maximum=fmax(maximum,ppml_quad_dot(q+(size_t)k*N,q+(size_t)k*N,NULL,N));
    }
    for(int k=0;k<K;k++) {
        int pivot=k;double norm=0;
        for(int j=k;j<K;j++) {
            double z=ppml_quad_dot(q+(size_t)j*N,q+(size_t)j*N,NULL,N);
            if(z>norm){norm=z;pivot=j;}
        }
        if(norm<=fmax(DBL_MIN,maximum*1e-24))break;
        int index=order[k];order[k]=order[pivot];order[pivot]=index;
        for(int i=0;i<N;i++){double z=q[(size_t)k*N+i];q[(size_t)k*N+i]=q[(size_t)pivot*N+i];q[(size_t)pivot*N+i]=z;}
        for(int i=0;i<k;i++){double z=R[i*K+k];R[i*K+k]=R[i*K+pivot];R[i*K+pivot]=z;}
        R[k*K+k]=sqrt(norm);
        for(int i=0;i<N;i++)q[(size_t)k*N+i]/=R[k*K+k];
        dd_real sum={0,0};
        for(int i=0;i<N;i++)sum=dd_add_d(sum,q[(size_t)k*N+i]*sqrt(w?w[i]:1)*y[i]);
        qy[k]=sum.hi+sum.lo;
        for(int j=k+1;j<K;j++)for(int pass=0;pass<2;pass++) {
            double dot=ppml_quad_dot(q+(size_t)k*N,q+(size_t)j*N,NULL,N);
            R[k*K+j]+=dot;
            for(int i=0;i<N;i++)q[(size_t)j*N+i]-=dot*q[(size_t)k*N+i];
        }
        rank++;
    }
    if(resid)for(int i=0;i<N;i++) {
        double z=0;for(int k=0;k<rank;k++)z+=q[(size_t)k*N+i]*qy[k];
        resid[i]=y[i]-z/sqrt(w?w[i]:1);
    }
    if(b)for(int k=rank-1;k>=0;k--) {
        double z=qy[k];for(int j=k+1;j<rank;j++)z-=R[k*K+j]*b[order[j]];
        b[order[k]]=z/R[k*K+k];
    }
    free(q);free(R);free(qy);free(order);return 0;
}

static int ppml_joint_slope(const HDFE_State *S,int g)
{
    if(S->factors[g].slope || g+1>=S->G)return 0;
    const FE_Factor *f=&S->factors[g], *h=&S->factors[g+1];
    return h->slope && h->slope_center && f->num_levels==h->num_levels &&
        !memcmp(f->levels,h->levels,(size_t)S->N*sizeof(int)) &&
        !(g+2<S->G && S->factors[g+2].slope &&
          S->factors[g+2].num_levels==f->num_levels &&
          !memcmp(f->levels,S->factors[g+2].levels,(size_t)S->N*sizeof(int)));
}
static void ppml_project(HDFE_State *S, double *y, int g)
{
    FE_Factor *f=&S->factors[g];
    double *a=S->thread_fe_means[g];
    int joint=ppml_joint_slope(S,g);
    FE_Factor *h=joint ? &S->factors[g+1] : NULL;
    double *b=joint ? S->thread_fe_means[g+1] : NULL;
    memset(a,0,f->num_levels*sizeof(double));
    for(int i=0;i<S->N;i++) {
        int l=f->levels[i]-1;
        double x=f->slope ? f->slope[i]-(f->slope_center?f->slope_center[l]:0) : 1;
        a[l]+=y[i]*x*S->weights[i];
    }
    for(int l=0;l<f->num_levels;l++)a[l]*=f->inv_weighted_counts[l];
    if(joint) {
        memset(b,0,f->num_levels*sizeof(double));
        for(int i=0;i<S->N;i++) {
            int l=f->levels[i]-1;
            b[l]+=(h->slope[i]-h->slope_center[l])*S->weights[i]*(y[i]-a[l]);
        }
        for(int l=0;l<f->num_levels;l++)b[l]*=h->inv_weighted_counts[l];
    }
    for(int i=0;i<S->N;i++) {
        int l=f->levels[i]-1;
        double x=f->slope ? f->slope[i]-(f->slope_center?f->slope_center[l]:0) : 1;
        double fitted=x*a[l];
        if(joint)fitted+=(h->slope[i]-h->slope_center[l])*b[l];
        y[i]-=fitted;
    }
}
/* Recover raw per-level intercept/slope coefficients from the saved sum.
 * All factors participate even when only some outputs were requested. */
static int ppml_save_effects(HDFE_State *S,double *d,const perm_idx_t *obs_map)
{
    int targets[10],any=0,rc=0;
    double *alpha[10]={0};
    for(int g=0;g<S->G;g++) {
        char key[32];snprintf(key,sizeof(key),"save_fe_%d",g+1);
        targets[g]=(int)ppml_option(key,0);any|=targets[g]>0;
    }
    if(!any)return 0;
    for(int g=0;g<S->G;g++) {
        alpha[g]=calloc((size_t)S->factors[g].num_levels,sizeof(double));
        if(!alpha[g]){rc=920;goto done;}
    }
    double scale=1;
    for(int i=0;i<S->N;i++)scale=fmax(scale,fabs(d[i]));
    int iteration;
    for(iteration=0;iteration<S->maxiter;iteration++) {
        for(int g=0;g<S->G;g++) {
            FE_Factor *f=&S->factors[g];int joint=ppml_joint_slope(S,g);
            ppml_project(S,d,g);
            double *a=S->thread_fe_means[g];
            for(int l=0;l<f->num_levels;l++)alpha[g][l]+=a[l];
            if(joint) {
                FE_Factor *h=&S->factors[g+1];double *b=S->thread_fe_means[g+1];
                for(int l=0;l<f->num_levels;l++){alpha[g][l]-=h->slope_center[l]*b[l];alpha[g+1][l]+=b[l];}
                g++;
            }
            else if(f->slope_center)for(int h=0;h<S->G;h++) {
                FE_Factor *intercept=&S->factors[h];
                if(!intercept->slope && intercept->num_levels==f->num_levels &&
                    !memcmp(intercept->levels,f->levels,(size_t)S->N*sizeof(int))) {
                    for(int l=0;l<f->num_levels;l++)alpha[h][l]-=f->slope_center[l]*a[l];
                    break;
                }
            }
        }
        double error=0;for(int i=0;i<S->N;i++)error=fmax(error,fabs(d[i]));
        if(error<=fmax(S->tolerance,1e-10)*scale)break;
    }
    if(iteration==S->maxiter){rc=430;goto done;}
    for(int g=0;g<S->G && !rc;g++)if(targets[g]) {
        for(int i=0;i<S->N && !rc;i++)rc=SF_vstore(targets[g],(int)obs_map[i],alpha[g][S->factors[g].levels[i]-1]);
    }
done:
    for(int g=0;g<S->G;g++)free(alpha[g]);
    return rc;
}

static void ppml_transform(HDFE_State *S,const double *y,double *ans)
{
    memcpy(ans,y,(size_t)S->N*sizeof(double));
    int groups[10],n=0;
    for(int g=0;g<S->G;g++) { groups[n++]=g;if(ppml_joint_slope(S,g))g++; }
    for(int g=0;g<n;g++)ppml_project(S,ans,groups[g]);
    for(int g=n-2;g>=0;g--)ppml_project(S,ans,groups[g]);
    for(int i=0;i<S->N;i++)ans[i]=y[i]-ans[i];
}

static HDFE_SolveResult ppml_reference_partial(HDFE_State *S, double *data, int N, int K, int threads)
{
    (void)threads;
    if (S->G==1) return partial_out_columns(S,data,N,K,1);
    double *r=S->thread_cg_r[0], *u=S->thread_cg_u[0], *v=S->thread_cg_v[0];
    int maxiter=0;
    for (int k=0;k<K;k++) {
        double *y=data+(size_t)k*N;
        double initial=ppml_quad_dot(y,y,NULL,N);
        double potential=ppml_quad_dot(y,y,S->weights,N);
        ppml_transform(S,y,r);
        double ssr=ppml_quad_dot(r,r,S->weights,N);
        memcpy(u,r,N*sizeof(double));
        int iter;
        for (iter=1;iter<=S->maxiter;iter++) {
            ppml_transform(S,u,v);
            double alpha=ssr/fmax(ppml_quad_dot(u,v,S->weights,N),DBL_EPSILON);
            double recent=alpha*ssr;
            potential-=recent;
            for(int i=0;i<N;i++) { y[i]-=alpha*u[i];r[i]-=alpha*v[i]; }
            double old=ssr;
            ssr=ppml_quad_dot(r,r,S->weights,N);
            double beta=ssr/fmax(old,DBL_EPSILON);
            for(int i=0;i<N;i++)u[i]=r[i]+beta*u[i];
            double error=fabs(recent)<1e-15 ? 0 : recent/fmax(potential,DBL_EPSILON);
            if(error<=S->tolerance*S->tolerance)break;
        }
        if(iter>S->maxiter)return (HDFE_SolveResult){430,0,iter};
        if(iter>maxiter)maxiter=iter;
        if(ppml_quad_dot(y,y,NULL,N)<=initial*fmin(1e-6,S->tolerance/10))memset(y,0,N*sizeof(double));
    }
    return (HDFE_SolveResult){0,1,maxiter};
}

#include "cpplmhdfe_solvers.h"
#include "cpplmhdfe_simplex.h"

/* Iterated rectifier separation, using the equality-constrained weighted
 * least-squares construction in Correia, Guimaraes and Zylkin.  Positive-y
 * observations receive the Van Loan penalty ceil(eps^(-1/2)).  Partialling
 * out FEs keeps memory O(N K), independent of the number of FE categories.
 * A pivoted, twice-orthogonalized QR avoids squaring the condition number of
 * the constrained regression.  No fitted-mean cutoff is used as evidence of
 * separation: arbitrarily small positive fitted means can be valid. */
static int ppml_relu_basis(HDFE_State *S,const double *X,int N,int K,double *q,int nthreads)
{
    if (!K) return 0;
    int rank=0;

        memcpy(q, X, (size_t)N * K * sizeof(*q));
        HDFE_SolveResult projection = ppml_reference_partial(S, q, N, K, nthreads);
        if (projection.status) return -projection.status;
        for (ST_int k = 0; k < K; k++) {
            ST_double norm = 0, original_norm = 0;
            for (ST_int i = 0; i < N; i++) {
                q[(size_t)k*N+i] *= sqrt(S->weights[i]);
                norm = hypot(norm, q[(size_t)k*N+i]);
                original_norm = hypot(original_norm, X[(size_t)k*N+i]*sqrt(S->weights[i]));
            }
            /* Do not amplify projection noise in an FE-collinear regressor
             * into an artificial separating direction. */
            if (norm < fmax(1e-9, 1e-8*original_norm)) memset(q + (size_t)k*N, 0, (size_t)N*sizeof(*q));
            else for (ST_int i = 0; i < N; i++) q[(size_t)k*N+i] /= norm;
        }
        for (ST_int k = 0; k < K; k++) {
            ST_int pivot = k;
            ST_double best = 0;
            for (ST_int j = k; j < K; j++) {
                ST_double norm = fast_dot(q+(size_t)j*N, q+(size_t)j*N, N);
                if (norm > best) { best = norm; pivot = j; }
            }
            if (best < 1e-20) break;
            for (ST_int i = 0; i < N; i++) {
                ST_double tmp = q[(size_t)k*N+i];
                q[(size_t)k*N+i] = q[(size_t)pivot*N+i];
                q[(size_t)pivot*N+i] = tmp;
                q[(size_t)k*N+i] /= sqrt(best);
            }
            rank++;
            for (ST_int j = k+1; j < K; j++) {
                for (ST_int pass = 0; pass < 2; pass++) {
                    ST_double dot = fast_dot(q+(size_t)k*N, q+(size_t)j*N, N);
                    for (ST_int i = 0; i < N; i++) q[(size_t)j*N+i] -= dot*q[(size_t)k*N+i];
                }
            }
        }
        return rank;
}

static ST_int ppml_relu_ex(HDFE_State *S, const ST_double *y, const ST_double *X,
    ST_int N, ST_int K, ST_int *separated, ST_int nthreads, ST_double *certificate)
{
    ST_int zeros = 0, rank = 0, rc = 0;
    ST_double saved_tol = S->tolerance;
    const ST_double penalty = 67108864.0;
    memset(separated, 0, (size_t)N * sizeof(*separated));
    for (ST_int i = 0; i < N; i++) zeros += y[i] == 0;
    if (!zeros) return 0;
    ST_double *q = malloc((size_t)N * (K + 1) * sizeof(*q));
    ST_double *u = malloc((size_t)N * sizeof(*u));
    ST_double *r = malloc((size_t)N * sizeof(*r));
    ST_double *fit = malloc((size_t)N * sizeof(*fit));
    ST_double *last = calloc((size_t)N*4,sizeof(*last));
    if (!q || !u || !r || !fit || !last) { rc = -920; goto done; }
    for (ST_int i = 0; i < N; i++) {
        S->weights[i] = y[i] == 0 ? 1.0 : penalty;
        u[i] = y[i] == 0;
    }
    ppml_update_fe_weights(S, S->weights, N);
    S->tolerance = fmax(1e-13,ppml_options.relu_tol*ppml_options.relu_tol);
    ST_double norm_zero = 0.0, norm_positive = 0.0;
    for (ST_int k = 0; k < K; k++) {
        ST_double zero_norm = 0.0, positive_norm = 0.0;
        for (ST_int i = 0; i < N; i++) {
            if (y[i] == 0) zero_norm += fabs(X[(size_t)k*N+i]);
            else positive_norm += fabs(X[(size_t)k*N+i]);
        }
        norm_zero = fmax(norm_zero, zero_norm);
        norm_positive = fmax(norm_positive, positive_norm);
    }
    ST_double ratio = norm_positive > 0 ? fmax(norm_zero / norm_positive, 1.0) : 1.0;
    rank=ppml_relu_basis(S,X,N,K,q,nthreads);
    if(rank<0){rc=rank;goto done;}
    double ee_cumulative=0,ee_boundary=zeros,progress1=0,progress2=0,acceleration_value=1;
    int candidates1=-1,candidates2=-1,stuck=0;
    for (ST_int iter = 0; iter < ppml_options.relu_maxiter; iter++) {
        for(int i=0;i<N;i++) { r[i]=u[i]+last[i]-last[(size_t)N+i];last[(size_t)N+i]=u[i]; }
        HDFE_SolveResult projection = ppml_reference_partial(S, r, N, 1, nthreads);
        memcpy(last,r,(size_t)N*sizeof(*r));
        if (projection.status) { rc = -projection.status; goto done; }
        for (ST_int i = 0; i < N; i++) r[i] *= sqrt(S->weights[i]);
        for (ST_int pass = 0; pass < 2; pass++) {
            for (ST_int k = 0; k < rank; k++) {
                ST_double dot = fast_dot(q+(size_t)k*N, r, N);
                for (ST_int i = 0; i < N; i++) r[i] -= dot*q[(size_t)k*N+i];
            }
        }
        ST_double delta = fmax(fast_dot(u,u,N)/penalty*ratio*ratio, 1e-8) + ppml_options.relu_tol;
        ST_double min_fit = 0.0, min_resid = 0.0;
        for (ST_int i = 0; i < N; i++) {
            r[i] /= sqrt(S->weights[i]);
            fit[i] = y[i] == 0 ? u[i] - r[i] : 0.0;
            if (fit[i] >= -0.1*delta && fit[i] <= delta) fit[i] = 0;
            if (y[i] == 0) {
                min_fit = fmin(min_fit, fit[i]);
                min_resid = fmin(min_resid, fabs(r[i]) < ppml_options.relu_zero_tol ? 0 : r[i]);
            }
        }
        int candidates=0;
        for(int i=0;i<N;i++)candidates+=y[i]==0 && fit[i]>delta;
        double ee=ppml_quad_dot(r,r,NULL,N);
        ee_cumulative+=iter?2*ee:ee;
        double progress=100*ee_cumulative/ee_boundary;
        if(ppml_options.relu_accelerate && iter>=3 && progress-progress2<1 && candidates==candidates2)stuck=1;
        progress2=progress1;progress1=progress;candidates2=candidates1;candidates1=candidates;
        if (min_fit >= 0 || min_resid >= 0) {
            for (ST_int i = 0; i < N; i++) {
                separated[i] = y[i] == 0 && fit[i] > 0 && (min_fit >= 0 || r[i] <= delta);
                if(certificate)certificate[i]=separated[i] ? fit[i] : 0;
                rc += separated[i];
            }
            goto done;
        }
        int changed_weights=0,accelerated=0;
        if(stuck)for(int i=0;i<N;i++)accelerated+=y[i]==0 && last[(size_t)3*N+i]<1.01*last[(size_t)2*N+i] && last[(size_t)2*N+i]<1.01*fit[i] && fit[i]<-.1*delta;
        acceleration_value=accelerated?fmin(256,4*acceleration_value):1;
        for(int i=0;i<N;i++) {
            double w=y[i]>0?penalty:1;
            if(stuck && y[i]==0 && last[(size_t)3*N+i]<1.01*last[(size_t)2*N+i] && last[(size_t)2*N+i]<1.01*fit[i] && fit[i]<-.1*delta)w=acceleration_value;
            changed_weights|=S->weights[i]!=w;S->weights[i]=w;
            last[(size_t)3*N+i]=last[(size_t)2*N+i];last[(size_t)2*N+i]=fit[i];
            u[i]=fmax(fit[i],0);
        }
        if(changed_weights) {
            ppml_update_fe_weights(S,S->weights,N);
            rank=ppml_relu_basis(S,X,N,K,q,nthreads);
            if(rank<0){rc=rank;goto done;}
        }
    }
    /* As in ppmlhdfe's default non-strict ReLU check, a limit does not
     * certify any observation as separated. IRLS still must converge. */
    if(ppml_options.relu_strict)rc=-9010;
done:
    S->tolerance = saved_tol;
    free(q); free(u); free(r); free(fit); free(last);
    return rc;
}

/* Kept as a small entry point for validation/test_cpplmhdfe_native.py. */
__attribute__((unused))
static ST_int ppml_relu(HDFE_State *S,const double *y,const double *X,int N,int K,int *mask,int threads)
{ return ppml_relu_ex(S,y,X,N,K,mask,threads,NULL); }

/* Global state pointer for CG solver (same pattern as creghdfe) */
static HDFE_State *g_ppml_state = NULL;

static void cleanup_ppml_state(void)
{
    if (g_ppml_state) {
        free(g_ppml_state->weights);
        g_ppml_state->weights = NULL;
        ctools_hdfe_state_cleanup(g_ppml_state);
        free(g_ppml_state);
        g_ppml_state = NULL;
    }
}

/* ========================================================================
 * Main PPML regression
 * ======================================================================== */

ST_retcode do_ppml_regression(int argc, char *argv[])
{
    ST_int K, G, N_orig, N, in1, in2;
    ST_retcode output_rc = 0;
    ST_int sample_var_idx = 0;
    ST_int k, g, i, j, idx;
    ST_double val;
    ST_int verbose;
    ST_int num_singletons;
    ST_int *mask = NULL;
    PPMLFactorData *factors = NULL;
    char scalar_name[64];
    double t_start, t_load, t_remap, t_singleton, t_dof, t_irls, t_vce;
    ST_int mobility_groups = 1;
    ST_int df_a = 0;
    ST_int num_threads;

    /* Data arrays */
    ST_double *data = NULL;     /* N x K matrix (column-major) */
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    perm_idx_t *obs_map = NULL;

    (void)argc;
    (void)argv;

    ppml_read_options();
    ctools_scal_save("__cpplmhdfe_api", 3);
    t_start = get_time_sec();

    in1 = SF_in1();
    in2 = SF_in2();

    /* Read parameters */
    SF_scal_use("__cpplmhdfe_K", &val); K = (ST_int)val;
    SF_scal_use("__cpplmhdfe_G", &val); G = (ST_int)val;
    SF_scal_use("__cpplmhdfe_verbose", &val); verbose = (ST_int)val;

    ST_int maxiter;
    ST_double tolerance;
    SF_scal_use("__cpplmhdfe_maxiter", &val); maxiter = (ST_int)val;
    SF_scal_use("__cpplmhdfe_tolerance", &val); tolerance = val;

    ST_int vcetype;
    SF_scal_use("__cpplmhdfe_vce_type", &val); vcetype = (ST_int)val;

    ST_int cluster_terms = vcetype == 2 ? 1 : 0;
    if (SF_scal_use("__cpplmhdfe_cluster_terms", &val) == 0)
        cluster_terms = (ST_int)val;
    if (vcetype == 2 && (cluster_terms < 1 || cluster_terms > 1023)) return 198;

    ST_int has_weights = 0, weight_type = 0;
    if (SF_scal_use("__cpplmhdfe_has_weights", &val) == 0)
        has_weights = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_weight_type", &val) == 0)
        weight_type = (ST_int)val;

    ST_int has_offset = 0;
    if (SF_scal_use("__cpplmhdfe_has_offset", &val) == 0)
        has_offset = (ST_int)val;

    /* IRLS parameters */
    ST_int irls_maxiter = 10000;
    ST_double irls_tol = 1e-8;
    if (SF_scal_use("__cpplmhdfe_irls_maxiter", &val) == 0)
        irls_maxiter = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_irls_tol", &val) == 0)
        irls_tol = val;

    /* DOF adjustment */
    ST_int dof_adjust_type = 0;
    if (SF_scal_use("__cpplmhdfe_dof_adjust_type", &val) == 0)
        dof_adjust_type = (ST_int)val;

    ST_int dof_flags = 13; /* pairwise=1, firstpair=2, clusters=4, continuous=8 */
    if (SF_scal_use("__cpplmhdfe_dof_flags", &val) == 0) dof_flags = (ST_int)val;

    ST_int slope_only[10] = {0}, slope_idx[10] = {0};
    ST_int keep_singletons = 0;
    if (SF_scal_use("__cpplmhdfe_keepsingletons", &val) == 0) keep_singletons = (ST_int)val;
    for (g = 0; g < G && g < 10; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_slope_%d", g+1);
        if (SF_scal_use(scalar_name, &val) == 0) slope_idx[g] = (ST_int)val;
        slope_only[g] = slope_idx[g] > 0;
    }

    ST_int check_fe = 1, check_relu = 1;
    if (SF_scal_use("__cpplmhdfe_check_fe", &val) == 0) check_fe = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_check_relu", &val) == 0) check_relu = (ST_int)val;

    /* Validation */
    if (G < 1 || G > 10) {
        SF_error("cpplmhdfe: invalid number of FE groups (must be 1-10)\n");
        return 198;
    }
    if (K < 2) {
        SF_error("cpplmhdfe: need at least depvar and one indepvar (K < 2)\n");
        return 198;
    }

    /* Determine threads */
#ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    if (num_threads > K) num_threads = K;
    if (num_threads < 1) num_threads = 1;
#else
    num_threads = 1;
#endif

    /* ================================================================
     * STEP 1: Data loading
     * varlist: depvar indepvars... fe1..feG [cluster] [weight] [offset]
     * ================================================================ */
    ST_int total_vars = K + G + cluster_terms + (has_weights ? 1 : 0) + (has_offset ? 1 : 0);

    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("cpplmhdfe: memory allocation failed\n");
        return 920;
    }
    for (i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;
    }

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, 0 | CTOOLS_LOAD_NO_SORT_ORDER);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("cpplmhdfe: parallel data load failed\n");
        return 920;
    }

    N_orig = (ST_int)filtered.data.nobs;
    obs_map = filtered.obs_map;

    if (N_orig <= 0) {
        ctools_filtered_data_free(&filtered);
        SF_error("cpplmhdfe: no observations\n");
        return 2000;
    }

    t_load = get_time_sec();

    /* Allocate working arrays */
    factors = (PPMLFactorData *)calloc(G, sizeof(PPMLFactorData));
    mask = (ST_int *)malloc(N_orig * sizeof(ST_int));
    data = (ST_double *)malloc((size_t)N_orig * K * sizeof(ST_double));

    ST_double *weights = NULL;
    if (has_weights) {
        weights = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    ST_double *offset_arr = NULL;
    if (has_offset) {
        offset_arr = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    if (!factors || !mask || !data || (has_weights && !weights) || (has_offset && !offset_arr)) {
        if (factors) free(factors);
        if (mask) free(mask);
        if (data) free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        ctools_filtered_data_free(&filtered);
        SF_error("cpplmhdfe: memory allocation failed\n");
        return 920;
    }

    for (i = 0; i < N_orig; i++) mask[i] = 1;

    /* Per-factor level arrays */
    for (g = 0; g < G; g++) {
        factors[g].levels = (ST_int *)malloc(N_orig * sizeof(ST_int));
        factors[g].num_obs = N_orig;
        factors[g].counts = NULL;
        factors[g].num_levels = 0;
        if (!factors[g].levels) {
            for (i = 0; i < g; i++) free(factors[i].levels);
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            ctools_filtered_data_free(&filtered);
            return 920;
        }
    }

    /* Copy numeric data to column-major */
    for (k = 0; k < K; k++) {
        double *src = filtered.data.vars[k].data.dbl;
        double *dst = &data[(size_t)k * N_orig];
        memcpy(dst, src, N_orig * sizeof(double));
    }

    /* Copy weights */
    if (has_weights) {
        ST_int weight_pos = K + G + cluster_terms;
        memcpy(weights, filtered.data.vars[weight_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Copy offset */
    if (has_offset) {
        ST_int offset_pos = K + G + cluster_terms + (has_weights ? 1 : 0);
        memcpy(offset_arr, filtered.data.vars[offset_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Validate y >= 0 */
    {
        ST_double *y = data;  /* First column */
        for (i = 0; i < N_orig; i++) {
            if (y[i] < 0.0) {
                SF_error("cpplmhdfe: dependent variable must be non-negative\n");
                for (g = 0; g < G; g++) free(factors[g].levels);
                free(factors); free(mask); free(data);
                if (weights) free(weights);
                if (offset_arr) free(offset_arr);
                ctools_filtered_data_free(&filtered);
                return 198;
            }
        }
    }

    /* FE remap + counting */
    int remap_failed = 0;
    for (g = 0; g < G; g++) {
        ST_int *counts_g = NULL;
        if (remap_and_count(filtered.data.vars[K + g].data.dbl, N_orig,
                            factors[g].levels, &factors[g].num_levels,
                            &counts_g, NULL, NULL) == 0) {
            factors[g].counts = counts_g;
        } else {
            remap_failed = 1;
        }
    }

    /* Save cluster raw values before freeing filtered data */
    double *cluster_raw_values = NULL;
    if (vcetype == 2) {
        cluster_raw_values = malloc((size_t)N_orig * cluster_terms * sizeof(double));
        if (cluster_raw_values) {
            for (ST_int term = 0; term < cluster_terms; term++)
                memcpy(cluster_raw_values + (size_t)term * N_orig,
                       filtered.data.vars[K + G + term].data.dbl,
                       N_orig * sizeof(double));
        }
    }

    /* Free filtered data, keep obs_map */
    stata_data_free(&filtered.data);

    /* Verify remapping */
    for (g = 0; g < G; g++) {
        if (remap_failed || factors[g].num_levels == 0) {
            for (i = 0; i < G; i++) {
                free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            if (cluster_raw_values) free(cluster_raw_values);
            ctools_aligned_free(obs_map);
            SF_error("cpplmhdfe: FE remapping failed\n");
            return 920;
        }
    }

    t_remap = get_time_sec();

    /* ================================================================
     * STEP 2: Singleton removal
     * ================================================================ */
    ST_int num_separated = 0;
    ST_int *selection_levels[10];
    ST_int selection_nlevels[10];
    for (g = 0; g < G; g++) {
        selection_levels[g] = factors[g].levels;
        selection_nlevels[g] = factors[g].num_levels;
    }
    ST_int selection_rc = ppml_select_fe_sample_ex(data, selection_levels,
        selection_nlevels, G, N_orig, mask, &num_singletons, &num_separated,
        slope_only, keep_singletons, check_fe && !keep_singletons);
    /* ppmlhdfe counts the initial FE screen as singleton removal. Explicit
     * regressor/rectifier separation is reported separately. */
    num_singletons += num_separated;
    num_separated = 0;
    N = 0;
    for (i = 0; i < N_orig; i++) if (mask[i]) N++;
    if (selection_rc || N == 0 || (vcetype == 2 && !cluster_raw_values)) {
        for (g = 0; g < G; g++) {
            free(factors[g].levels); free(factors[g].counts);
        }
        free(factors); free(mask); free(data); free(weights); free(offset_arr);
        free(cluster_raw_values); ctools_aligned_free(obs_map);
        if (selection_rc || N > 0) {
            SF_error("cpplmhdfe: sample selection allocation failed\n");
            return 920;
        }
        SF_error("cpplmhdfe: no observations remain after singleton/separation removal\n");
        return 2001;
    }
    if (verbose && (num_singletons || num_separated)) {
        char msg[160];
        snprintf(msg, sizeof(msg), "cpplmhdfe: dropped %d singletons and %d separated observations\n",
                 num_singletons, num_separated);
        SF_display(msg);
    }

    t_singleton = get_time_sec();

    /* Recount levels after singleton removal */
    ST_int orig_num_levels[10];
    for (g = 0; g < G; g++) {
        orig_num_levels[g] = factors[g].num_levels;
        memset(factors[g].counts, 0, orig_num_levels[g] * sizeof(ST_int));
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                ST_int level = factors[g].levels[i] - 1;
                if (level >= 0 && level < orig_num_levels[g])
                    factors[g].counts[level]++;
            }
        }
        ST_int nlev = 0;
        for (i = 0; i < orig_num_levels[g]; i++)
            if (factors[g].counts[i] > 0) nlev++;
        factors[g].num_levels = nlev;
    }

    /* ================================================================
     * STEP 3: DOF computation
     * ================================================================ */
    df_a = 0;
    for (g = 0; g < G; g++)
        df_a += factors[g].num_levels;

    if (dof_adjust_type == 1) {
        if (G >= 2) {
            mobility_groups = 1;
            df_a -= 1;
        }
    } else if (G >= 2) {
        ST_int *fe1_c = (ST_int *)malloc(N * sizeof(ST_int));
        ST_int *fe2_c = (ST_int *)malloc(N * sizeof(ST_int));
        if (fe1_c && fe2_c) {
            ST_int *remap1 = (ST_int *)calloc(orig_num_levels[0] + 1, sizeof(ST_int));
            ST_int *remap2 = (ST_int *)calloc(orig_num_levels[1] + 1, sizeof(ST_int));
            if (remap1 && remap2) {
                ST_int next1 = 1, next2 = 1;
                idx = 0;
                for (i = 0; i < N_orig; i++) {
                    if (mask[i]) {
                        ST_int lev1 = factors[0].levels[i];
                        ST_int lev2 = factors[1].levels[i];
                        if (remap1[lev1] == 0) remap1[lev1] = next1++;
                        if (remap2[lev2] == 0) remap2[lev2] = next2++;
                        fe1_c[idx] = remap1[lev1];
                        fe2_c[idx] = remap2[lev2];
                        idx++;
                    }
                }
                mobility_groups = count_connected_components(
                    fe1_c, fe2_c, N, factors[0].num_levels, factors[1].num_levels);
                if (mobility_groups < 0) mobility_groups = 1;
                free(remap1); free(remap2);
            }
        }
        if (fe1_c) free(fe1_c);
        if (fe2_c) free(fe2_c);
        df_a -= mobility_groups;
        if (G > 2 && (dof_adjust_type == 0 || dof_adjust_type == 3)) {
            ST_int extra = G - 2;
            df_a -= extra;
            mobility_groups += extra;
        }
    } else {
        mobility_groups = 0;
    }

    t_dof = get_time_sec();

    /* ================================================================
     * STEP 4: Build HDFE state (compacted)
     * ================================================================ */
    cleanup_ppml_state();
    g_ppml_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    if (!g_ppml_state) {
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) ctools_aligned_free(obs_map);
        return 920;
    }

    g_ppml_state->G = G;
    g_ppml_state->N = N;
    g_ppml_state->K = K;
    g_ppml_state->in1 = in1;
    g_ppml_state->in2 = in2;
    g_ppml_state->maxiter = maxiter;
    g_ppml_state->tolerance = tolerance;
    g_ppml_state->verbose = verbose;
    g_ppml_state->num_threads = num_threads;
    g_ppml_state->factors_initialized = 1;
    g_ppml_state->df_a = df_a;
    g_ppml_state->mobility_groups = mobility_groups;
    g_ppml_state->has_weights = 1;  /* IRLS always uses weights */
    g_ppml_state->weight_type = 1;  /* Treat as aweight for CG solver */
    g_ppml_state->weights = NULL;
    g_ppml_state->sum_weights = 0.0;

    /* Allocate and remap factors */
    g_ppml_state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    if (!g_ppml_state->factors) {
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) ctools_aligned_free(obs_map);
        return 920;
    }

    for (g = 0; g < G; g++) {
        g_ppml_state->factors[g].num_levels = factors[g].num_levels;
        g_ppml_state->factors[g].max_level = factors[g].num_levels - 1;
        g_ppml_state->factors[g].has_intercept = 1;
        g_ppml_state->factors[g].levels = (ST_int *)malloc(N * sizeof(ST_int));
        g_ppml_state->factors[g].counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        g_ppml_state->factors[g].weighted_counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        g_ppml_state->factors[g].means = NULL;

        if (!g_ppml_state->factors[g].levels || !g_ppml_state->factors[g].counts ||
            !g_ppml_state->factors[g].weighted_counts) {
            cleanup_ppml_state();
            for (i = 0; i < G; i++) {
                free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            if (cluster_raw_values) free(cluster_raw_values);
            if (obs_map) ctools_aligned_free(obs_map);
            return 920;
        }

        /* Remap to contiguous levels */
        ST_int *remap = (ST_int *)calloc(orig_num_levels[g] + 1, sizeof(ST_int));
        if (!remap) {
            cleanup_ppml_state();
            for (i = 0; i < G; i++) {
                free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            if (cluster_raw_values) free(cluster_raw_values);
            if (obs_map) ctools_aligned_free(obs_map);
            return 920;
        }
        ST_int next_level = 1;
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                ST_int old_level = factors[g].levels[i];
                if (remap[old_level] == 0)
                    remap[old_level] = next_level++;
                g_ppml_state->factors[g].levels[idx] = remap[old_level];
                g_ppml_state->factors[g].counts[remap[old_level] - 1] += 1.0;
                idx++;
            }
        }
        free(remap);
    }

    /* ================================================================
     * STEP 5: Compact data (remove singleton rows)
     * ================================================================ */
    ST_double *data_compact = (ST_double *)malloc((size_t)N * K * sizeof(ST_double));
    ST_double *w_user_compact = NULL;
    ST_double *offset_compact = NULL;

    if (has_weights) w_user_compact = (ST_double *)malloc(N * sizeof(ST_double));
    if (has_offset) offset_compact = (ST_double *)malloc(N * sizeof(ST_double));

    if (!data_compact || (has_weights && !w_user_compact) || (has_offset && !offset_compact)) {
        if (data_compact) free(data_compact);
        if (w_user_compact) free(w_user_compact);
        if (offset_compact) free(offset_compact);
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) ctools_aligned_free(obs_map);
        return 920;
    }

    /* Compact data columns */
    for (k = 0; k < K; k++) {
        const ST_double *src_col = data + (size_t)k * N_orig;
        ST_double *dst_col = data_compact + (size_t)k * N;
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                dst_col[idx++] = src_col[i];
            }
        }
    }

    /* Compact weights */
    if (has_weights) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) w_user_compact[idx++] = weights[i];
        }
    }

    /* Compact offset */
    if (has_offset) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) offset_compact[idx++] = offset_arr[i];
        }
    }

    free(data);
    data = data_compact;
    if (weights) free(weights);
    weights = NULL;
    if (offset_arr) free(offset_arr);
    offset_arr = NULL;

    /* Compact all row identities with exactly the same original-row mask.
     * After this point every array, including clusters and obs_map, has N rows. */
    idx = 0;
    for (i = 0; i < N_orig; i++) {
        if (mask[i]) {
            obs_map[idx] = obs_map[i];
            if (cluster_raw_values) {
                for (ST_int term = 0; term < cluster_terms; term++)
                    cluster_raw_values[(size_t)term * N_orig + idx] =
                        cluster_raw_values[(size_t)term * N_orig + i];
            }
            idx++;
        }
    }
    free(mask);
    mask = NULL;

    /* ================================================================
     * STEP 7: Allocate IRLS arrays and HDFE buffers
     * ================================================================ */
    ST_double *mu = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *eta = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *irls_w = (ST_double *)malloc(N * sizeof(ST_double));
    /* z_orig: saves pre-partialled working depvar for FE reconstruction.
     * aug_copy: [z, X1, ..., Xk] for CG partialling (modified in place). */
    ST_int K_x = K - 1;
    ST_int K_aug = K_x + 1;  /* working depvar + X vars */
    ST_double *z_orig = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *aug_copy = (ST_double *)malloc((size_t)N * K_aug * sizeof(ST_double));

    /* Pre-allocate IRLS working buffers (reused across iterations) */
    ST_double *xtx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    ST_double *xty = (ST_double *)malloc(K_x * sizeof(ST_double));
    ST_double *irls_beta = (ST_double *)malloc(K_x * sizeof(ST_double));
    ST_double *inv_xx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    ST_int *sep_m = (ST_int *)malloc(N * sizeof(ST_int));

    if (!mu || !eta || !irls_w || !z_orig || !aug_copy ||
        !xtx || !xty || !irls_beta || !inv_xx || !sep_m) {
        if (mu) free(mu);
        if (eta) free(eta);
        if (irls_w) free(irls_w);
        if (z_orig) free(z_orig);
        if (aug_copy) free(aug_copy);
        if (xtx) free(xtx);
        if (xty) free(xty);
        if (irls_beta) free(irls_beta);
        if (inv_xx) free(inv_xx);
        if (sep_m) free(sep_m);
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); if (mask) free(mask); free(data);
        if (w_user_compact) free(w_user_compact);
        if (offset_compact) free(offset_compact);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) ctools_aligned_free(obs_map);
        return 920;
    }

    /* Continuous FE loadings follow the same compact sample as y and X. */
    for (g = 0; g < G; g++) {
        if (!slope_only[g]) continue;
        FE_Factor *f = &g_ppml_state->factors[g];
        f->slope = malloc((size_t)N*sizeof(*f->slope));
        if (!f->slope) {
            /* Allocation failure uses the common path once IRLS owns weights. */
            output_rc = 920;
            break;
        }
        f->has_intercept = 0;
        f->num_slopes = 1;
        for (j = 0; j < G; j++) {
            FE_Factor *intercept = &g_ppml_state->factors[j];
            if (!slope_only[j] && f->num_levels == intercept->num_levels &&
                !memcmp(f->levels, intercept->levels, (size_t)N*sizeof(ST_int))) {
                f->slope_center = calloc((size_t)f->num_levels, sizeof(ST_double));
                if (!f->slope_center) output_rc = 920;
                break;
            }
        }
        if (output_rc) break;
        for (i = 0; i < N; i++) {
            ST_retcode read_rc = SF_vdata(slope_idx[g], (ST_int)obs_map[i], &f->slope[i]);
            if (read_rc) { output_rc = read_rc; break; }
        }
        if (output_rc) break;
    }

    /* Set IRLS weights for HDFE state */
    g_ppml_state->weights = irls_w;

    /* Allocate HDFE buffers (alloc_proj=0, max_columns=K_aug) */
    if (ctools_hdfe_alloc_buffers(g_ppml_state, 0, K_aug) != 0) {
        free(mu); free(eta);
        free(z_orig); free(aug_copy);
        free(xtx); free(xty); free(irls_beta); free(inv_xx); free(sep_m);
        /* irls_w is owned by g_ppml_state->weights, freed by cleanup */
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); if (mask) free(mask); free(data);
        if (w_user_compact) free(w_user_compact);
        if (offset_compact) free(offset_compact);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) ctools_aligned_free(obs_map);
        return 920;
    }

    ST_double *stdev_x = NULL;
    ST_int *irls_keep_idx = NULL;
    ST_double *xtx_k = NULL;
    ST_double *xty_k = NULL;
    ST_double *beta_k = NULL;
    ST_double *inv_xx_k = NULL;
    ST_int *is_collinear = NULL;
    ST_int *keep_idx = NULL;
    ST_double *beta_final = NULL;
    ST_double *inv_xx_final = NULL;
    ST_double *V_keep = NULL;
    ST_double *xtx_keep = NULL;
    ST_double *xty_keep = NULL;
    ST_double *w_reg = NULL;
    ST_double *means_x = NULL;
    ST_double *inv_xx_ext = NULL;
    ST_double *X_eff = NULL;
    ST_double *vce_resid = NULL;
    ST_double *V_ext = NULL;
    ST_double *separation_certificate = NULL;
    ST_double *previous_working = NULL;
    if (output_rc) goto ppml_cleanup;
    /* ================================================================
     * STEP 8: Initialize IRLS
     * Initialize eta from OLS of log(max(y,1)) on X with unit weights
     * ================================================================ */
    ST_double *y = data;  /* First column of data */
    ST_double *X = data + N;  /* Columns 1..K-1 */

    /* Standardize by sample stdev before IRLS. Retain the reference's 1e-3
     * regressor scale floor, but keep response scaling invariant even for
     * very small outcomes. */
    ST_double stdev_y = 1.0;
    stdev_x = (ST_double *)malloc(K_x * sizeof(ST_double));
    if (!stdev_x) { output_rc = 920; goto ppml_cleanup; }
    {
        /* Center before squaring: E[x^2]-E[x]^2 loses all precision when
         * the mean is large relative to the within-sample variation. */
        for (k = -1; k < K_x; k++) {
            ST_double *column = k < 0 ? y : X + (size_t)k * N;
            dd_real total = {0.0, 0.0};
            for (i = 0; i < N; i++) total = dd_add_d(total, column[i] / N);
            ST_double mean = total.hi + total.lo;
            ST_double max_delta = 0.0;
            for (i = 0; i < N; i++) max_delta = fmax(max_delta, fabs(column[i] - mean));
            dd_real squares = {0.0, 0.0};
            if (max_delta > 0.0) {
                for (i = 0; i < N; i++) {
                    ST_double delta = (column[i] - mean) / max_delta;
                    squares = dd_add_d(squares, delta * delta);
                }
            }
            ST_double sd = max_delta * sqrt((squares.hi + squares.lo) / (N - 1));
            /* The response scale must not change the fit or trigger spurious
             * separation for small positive outcomes. */
            if (k < 0) {
                stdev_y = sd > 0.0 ? sd : fmax(fabs(mean), 1e-3);
                sd = stdev_y;
            } else {
                sd = fmax(sd, 1e-3);
                stdev_x[k] = sd;
            }
            if(!ppml_options.standardize) {
                sd=1;
                if(k<0)stdev_y=1;else stdev_x[k]=1;
            }
            ST_double fixed_scales = 0;
            SF_scal_use("__cpplmhdfe_scale_fixed", &fixed_scales);
            if (fixed_scales) {
                SF_mat_el("__cpplmhdfe_input_scales", 1, k+2, &sd);
                if (k < 0) stdev_y = sd; else stdev_x[k] = sd;
            } else ctools_mat_store("__cpplmhdfe_input_scales", 1, k+2, sd);
            for (i = 0; i < N; i++) column[i] /= sd;
        }
    }

    separation_certificate=calloc((size_t)N,sizeof(double));
    previous_working=calloc((size_t)N,sizeof(double));
    if(!separation_certificate || !previous_working) {output_rc=920;goto ppml_cleanup;}
    if(ppml_options.check_simplex) {
        int nsep=ppml_simplex(g_ppml_state,y,X,w_user_compact,weight_type,N,K_x,sep_m);
        if(nsep<0){output_rc=-nsep;goto ppml_cleanup;}
        if(nsep) {
            sample_var_idx=(int)ppml_option("sample_idx",0);
            for(i=0;i<N && !output_rc;i++)if(!sep_m[i])output_rc=SF_vstore(sample_var_idx,(int)obs_map[i],1);
            ctools_scal_save("__cpplmhdfe_refit",1);
            ctools_scal_save("__cpplmhdfe_num_singletons",num_singletons);
            ctools_scal_save("__cpplmhdfe_num_separated",num_separated+nsep);
            goto ppml_cleanup;
        }
    }
    ST_int relu_separated = check_relu ? ppml_relu_ex(g_ppml_state, y, X, N, K_x, sep_m, num_threads,separation_certificate) : 0;
    if (relu_separated < 0) { output_rc = -relu_separated; goto ppml_cleanup; }
    /* ReLU diagnostics refer to its input sample, excluding rows already
     * removed by FE screening, singleton removal, or simplex. */
    if(check_relu) for(i=0;i<N && !output_rc;i++) {
        if(ppml_options.tag_idx)output_rc=SF_vstore(ppml_options.tag_idx,(int)obs_map[i],sep_m[i]);
        if(ppml_options.certificate_idx)output_rc=SF_vstore(ppml_options.certificate_idx,(int)obs_map[i],separation_certificate[i]);
    }
    if(output_rc || ppml_option("diagnostic",0))goto ppml_cleanup;
    if (relu_separated) {
        if (verbose) {
            char message[128];
            snprintf(message, sizeof(message), "cpplmhdfe: ReLU identified %d separated observations\n", relu_separated);
            SF_display(message);
        }
        if (SF_scal_use("__cpplmhdfe_sample_idx", &val) == 0) sample_var_idx = (ST_int)val;
        for (i = 0; i < N && !output_rc; i++)
            if (!sep_m[i]) output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 1.0);
        ctools_scal_save("__cpplmhdfe_refit", 2);
        ctools_scal_save("__cpplmhdfe_num_singletons", num_singletons);
        ctools_scal_save("__cpplmhdfe_num_separated", num_separated + relu_separated);
        goto ppml_cleanup;
    }

    /* Initial means: simple, OLS, or a user-supplied mean variable. */
    {
        double ym=0,ws=0;
        for(i=0;i<N;i++){double w=w_user_compact?w_user_compact[i]:1;ym+=w*y[i];ws+=w;irls_w[i]=w;}
        ym/=ws;
        if(ppml_options.guess==1) {
            for(i=0;i<N;i++)aug_copy[i]=log(y[i]+ym/100);
            memcpy(aug_copy+N,X,(size_t)N*K_x*sizeof(double));
            ppml_update_fe_weights(g_ppml_state,irls_w,N);
            g_ppml_state->tolerance=fmax(ppml_options.start_tol,tolerance);
            HDFE_SolveResult pr=ppml_partial(g_ppml_state,aug_copy,N,K_aug,num_threads);
            if(pr.status){output_rc=pr.status;goto ppml_cleanup;}
            output_rc=ppml_qrsolve(aug_copy+N,aug_copy,irls_w,N,K_x,NULL,mu);
            if(output_rc)goto ppml_cleanup;
            for(i=0;i<N;i++)mu[i]=exp(log1p(y[i])-mu[i]);
        }
        else if(ppml_options.guess==2) {
            double average=0,sd=0;
            for(i=0;i<N;i++) {
                output_rc=SF_vdata(ppml_options.guess_idx,(int)obs_map[i],&mu[i]);
                if(output_rc || !isfinite(mu[i]) || mu[i]>=SV_missval){output_rc=198;goto ppml_cleanup;}
                average+=mu[i]/N;
            }
            if(ppml_options.standardize) {
                for(i=0;i<N;i++)sd+=(mu[i]-average)*(mu[i]-average);
                sd=fmax(sqrt(sd/fmax(N-1,1)),1e-3);
                for(i=0;i<N;i++)mu[i]/=sd;
            }
        }
        else for(i=0;i<N;i++)mu[i]=.5*(y[i]+ym);
        for(i=0;i<N;i++) {
            mu[i]=fmax(mu[i],fmax(.05*y[i],1e-3));
            eta[i]=log(mu[i])-(offset_compact?offset_compact[i]:0);
        }
    }

    ST_double user_weight_scale = 1.0;
    if (w_user_compact) {
        user_weight_scale = 0.0;
        for (i = 0; i < N; i++) user_weight_scale += w_user_compact[i] / N;
    }

    /* A common scale cancels from weighted projections and the sandwich. */
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        irls_w[i] = mu[i] * (w_u / user_weight_scale);
        if (irls_w[i] < DBL_MIN) irls_w[i] = DBL_MIN;  /* Floor to prevent zero weights */
    }

    /* Set weighted_counts for initial weights */
    ppml_update_fe_weights(g_ppml_state, irls_w, N);

    /* ================================================================
     * STEP 8b: Pre-IRLS collinearity detection
     * Partial out FE from X with initial weights, then detect collinearity
     * in X'X to identify variables collinear with the absorbed FEs.
     * ================================================================ */

    ST_int K_irls = K_x;  /* number of non-collinear X vars for IRLS */
    if(ppml_options.remove_collinear) {
        /* Set up [z_dummy, X] for partialling - z doesn't matter, just X */
        for (i = 0; i < N; i++)
            aug_copy[i] = 0.0;  /* dummy z */
        for (k = 0; k < K_x; k++)
            memcpy(&aug_copy[(size_t)(k + 1) * N], &X[(size_t)k * N], N * sizeof(ST_double));

        HDFE_SolveResult projection = ppml_partial(g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (projection.status) {
            SF_error("cpplmhdfe: fixed-effect projection did not converge\n");
            output_rc = projection.status;
            goto ppml_cleanup;
        }

        /* Build X'X from partialled X columns */
        ST_double *xtx_pre = (ST_double *)calloc(K_x * K_x, sizeof(ST_double));
        if (!xtx_pre) { output_rc = 920; goto ppml_cleanup; }
        if (xtx_pre) {
            for (i = 0; i < K_x; i++) {
                for (j = 0; j <= i; j++) {
                    ST_double d = fast_dot(&aug_copy[(size_t)(i+1)*N], &aug_copy[(size_t)(j+1)*N], N);
                    xtx_pre[i * K_x + j] = d;
                    xtx_pre[j * K_x + i] = d;
                }
            }
            ST_int *pre_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));
            if (!pre_collinear) { free(xtx_pre); output_rc = 920; goto ppml_cleanup; }
            if (pre_collinear) {
                detect_collinearity(xtx_pre, K_x, pre_collinear, verbose);
                ST_int n_coll = 0;
                for (k = 0; k < K_x; k++)
                    if (pre_collinear[k]) n_coll++;
                K_irls = K_x - n_coll;
                if (n_coll > 0 && K_irls > 0) {
                    irls_keep_idx = (ST_int *)malloc(K_irls * sizeof(ST_int));
                    if (!irls_keep_idx) { free(pre_collinear); free(xtx_pre); output_rc = 920; goto ppml_cleanup; }
                    if (irls_keep_idx) {
                        idx = 0;
                        for (k = 0; k < K_x; k++)
                            if (!pre_collinear[k])
                                irls_keep_idx[idx++] = k;
                    }
                    if (verbose) {
                        char msg[128];
                        snprintf(msg, sizeof(msg),
                            "(cpplmhdfe: %d variable(s) collinear with FE, excluded from IRLS)\n", n_coll);
                        SF_display(msg);
                    }
                }
                free(pre_collinear);
            }
            free(xtx_pre);
        }
    }

    /* Allocate compact IRLS arrays if collinearity detected */




    if (irls_keep_idx && K_irls < K_x) {
        xtx_k = (ST_double *)malloc(K_irls * K_irls * sizeof(ST_double));
        xty_k = (ST_double *)malloc(K_irls * sizeof(ST_double));
        beta_k = (ST_double *)malloc(K_irls * sizeof(ST_double));
        inv_xx_k = (ST_double *)malloc(K_irls * K_irls * sizeof(ST_double));
        if (!xtx_k || !xty_k || !beta_k || !inv_xx_k) { output_rc = 920; goto ppml_cleanup; }
    }

    /* ================================================================
     * STEP 9: IRLS LOOP
     * ================================================================ */
    ST_double deviance = 1e30, deviance_old;
    ST_int irls_iter;
    ST_int irls_converged = 0, ok_iterations=0, halving_count=0;
    memset(sep_m,0,(size_t)N*sizeof(int));
    ST_double log_eps_history[3] = {0.0, 0.0, 0.0};
    ST_int eps_history_count = 0;

    /* Adaptive CG tolerance matching ppmlhdfe:
     * Start at 1e-4 (fast early iterations), tighten toward 1e-9
     * (precise final iterations) based on relative deviance change. */
    ST_double target_inner_tol = fmax(1e-12, fmin(1e-9, 0.1 * tolerance));
    if (SF_scal_use("__cpplmhdfe_inner_tol", &val) == 0 && val > 0) target_inner_tol = val;
    ST_double start_inner_tol = ppml_options.start_tol;
    g_ppml_state->tolerance = ppml_options.heuristic ? fmax(start_inner_tol,target_inner_tol) : target_inner_tol;

    for (irls_iter = 0; irls_iter < irls_maxiter; irls_iter++) {
        deviance_old = deviance;

        /* (a) Compute working depvar: z = eta + (y - mu)/mu
         * Save in z_orig for FE reconstruction; copy into aug_copy for CG. */
        for (i = 0; i < N; i++) {
            ST_double z = eta[i] + (y[i] - mu[i]) / mu[i];
            z_orig[i] = z;
            aug_copy[i] = ppml_options.exact_partial || !irls_iter ? z : aug_copy[i]+z-previous_working[i];
            previous_working[i]=z;
        }

        /* Reuse the previous residuals when inexact partialling is requested. */
        if (ppml_options.exact_partial || !irls_iter) for (k = 0; k < K_x; k++) {
            memcpy(&aug_copy[(size_t)(k + 1) * N], &X[(size_t)k * N], N * sizeof(ST_double));
        }

        /* (d) Partial out via CG solver using current IRLS weights */
        HDFE_SolveResult projection = ppml_partial(g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (projection.status) {
            SF_error("cpplmhdfe: fixed-effect projection did not converge\n");
            output_rc = projection.status;
            goto ppml_cleanup;
        }

        /* (e) Solve WLS: beta = (X~'W X~)^{-1} X~'W z~ using partialled data */
        memset(irls_beta, 0, K_x * sizeof(ST_double));

        if (ppml_options.exact_solver || !ppml_options.remove_collinear) {
            output_rc=ppml_qrsolve(aug_copy+N,aug_copy,irls_w,N,K_x,irls_beta,NULL);
            if(output_rc)goto ppml_cleanup;
        } else if (irls_keep_idx && K_irls < K_x && xtx_k && xty_k && beta_k && inv_xx_k) {
            /* Compact solve: only non-collinear columns.
             * Use dd_real (double-double) precision for X'WX and X'Wz to match
             * ppmlhdfe/reghdfe's quadcross(). Small precision differences in
             * each IRLS iteration accumulate in the mu trajectory, which then
             * gets amplified through cond(X'WX) in the final VCE. */
            ST_int ki, kj;
            memset(xtx_k, 0, K_irls * K_irls * sizeof(ST_double));
            memset(xty_k, 0, K_irls * sizeof(ST_double));

            for (ki = 0; ki < K_irls; ki++) {
                ST_int ci = irls_keep_idx[ki];
                const ST_double *xi = &aug_copy[(size_t)(ci + 1) * N];
                const ST_double *z_par = &aug_copy[0];
                {
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xi[i];
                        dd_real prod = two_prod(wxi, z_par[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xty_k[ki] = acc.hi + acc.lo;
                }
                for (kj = 0; kj <= ki; kj++) {
                    ST_int cj = irls_keep_idx[kj];
                    const ST_double *xj = &aug_copy[(size_t)(cj + 1) * N];
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xi[i];
                        dd_real prod = two_prod(wxi, xj[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xtx_k[ki * K_irls + kj] = acc.hi + acc.lo;
                    xtx_k[kj * K_irls + ki] = acc.hi + acc.lo;
                }
            }

            memcpy(inv_xx_k, xtx_k, K_irls * K_irls * sizeof(ST_double));
            if (ctools_cholesky(inv_xx_k, K_irls) != 0) {
                break;
            }
            ctools_invert_from_cholesky(inv_xx_k, K_irls, inv_xx_k);

            for (ki = 0; ki < K_irls; ki++) {
                beta_k[ki] = 0.0;
                for (kj = 0; kj < K_irls; kj++)
                    beta_k[ki] += inv_xx_k[ki * K_irls + kj] * xty_k[kj];
            }
            /* Expand to full beta */
            for (ki = 0; ki < K_irls; ki++)
                irls_beta[irls_keep_idx[ki]] = beta_k[ki];
        } else if (K_irls > 0) {
            /* Full solve: all columns (no collinearity detected).
             * Use dd_real for X'WX and X'Wz to match quadcross() precision. */
            memset(xtx, 0, K_x * K_x * sizeof(ST_double));
            memset(xty, 0, K_x * sizeof(ST_double));
            for (j = 0; j < K_x; j++) {
                const ST_double *xj_ptr = &aug_copy[(size_t)(j + 1) * N];
                const ST_double *z_par = &aug_copy[0];
                {
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xj_ptr[i];
                        dd_real prod = two_prod(wxi, z_par[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xty[j] = acc.hi + acc.lo;
                }
                for (k = 0; k <= j; k++) {
                    const ST_double *xk_ptr = &aug_copy[(size_t)(k + 1) * N];
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xj_ptr[i];
                        dd_real prod = two_prod(wxi, xk_ptr[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xtx[j * K_x + k] = acc.hi + acc.lo;
                    xtx[k * K_x + j] = acc.hi + acc.lo;
                }
            }

            memcpy(inv_xx, xtx, K_x * K_x * sizeof(ST_double));
            if (ctools_cholesky(inv_xx, K_x) != 0) break;
            ctools_invert_from_cholesky(inv_xx, K_x, inv_xx);

            for (i = 0; i < K_x; i++) {
                irls_beta[i] = 0.0;
                for (j = 0; j < K_x; j++)
                    irls_beta[i] += inv_xx[i * K_x + j] * xty[j];
            }
        }

        /* (f) Update eta and mu.
         * z_pred[i] = X_tilde * beta + FE_component
         * where FE_component = z_orig - z_partialled (original z minus CG result)
         */
        ST_double max_eta_change = 0.0;
        for (i = 0; i < N; i++) {
            ST_double xb = 0.0;
            for (k = 0; k < K_x; k++) {
                xb += aug_copy[(size_t)(k + 1) * N + i] * irls_beta[k];
            }
            /* FE component: original z minus partialled z */
            ST_double fe_part = z_orig[i] - aug_copy[i];
            ST_double z_hat = xb + fe_part;
            ST_double off = (offset_compact != NULL) ? offset_compact[i] : 0.0;
            /* eta = X*beta + FE (without offset); z didn't include offset */
            if (!isfinite(z_hat)) { output_rc = 430; break; }
            /* Zero outcomes with negligible means may have no finite FE
             * coefficient after a tolerance-limited separation check. They
             * must not prevent convergence of the identified parameters. */
            if (y[i] > 0 || mu[i] > 1e-12)
                max_eta_change = fmax(max_eta_change, fabs(z_hat - eta[i]));
            eta[i] = z_hat;
            mu[i] = exp(eta[i] + off);
            /* Clamp mu to prevent overflow/underflow */
            if (mu[i] > 1e18) mu[i] = 1e18;
            if (mu[i] < 1e-18) mu[i] = 1e-18;
        }

        if (output_rc) break;
        if(ppml_options.check_mu) {
            double min_eta=DBL_MAX;
            for(i=0;i<N;i++)if(y[i]>0)min_eta=fmin(min_eta,eta[i]+(offset_compact?offset_compact[i]:0));
            double cutoff=log(ppml_options.mu_tol)+fmin(min_eta+5,0);
            for(i=0;i<N;i++) {
                sep_m[i]|=y[i]==0 && eta[i]+(offset_compact?offset_compact[i]:0)<=cutoff;
                if(sep_m[i])mu[i]=1.4210854715202004e-14;
            }
        }

        /* (g) Compute deviance */
        deviance = ppml_compute_deviance(y, mu, w_user_compact, N) / user_weight_scale;

        /* (h) Check convergence (matching ppmlhdfe's criterion) */
        if (irls_iter > 0) {
            ST_double denom_eps = deviance < deviance_old ? deviance : deviance_old;
            if (denom_eps < 0.1) denom_eps = 0.1;
            ST_double rel_change = fabs(deviance - deviance_old) / denom_eps;
            if (verbose) {
                char msg[256];
                snprintf(msg, sizeof(msg),
                    "  iter %d: deviance=%.10e eps=%.4e\n",
                    irls_iter + 1, deviance, rel_change);
                SF_display(msg);
            }
            /* Match the reference's decision to finish with an exact solve:
             * extrapolate log(eps) with a cubic through iterations 0..3,
             * where the intercept is zero. One-way demeaning is exact; for
             * multiple FE dimensions also require the target inner tolerance. */
            ST_double predicted_eps = DBL_MAX;
            if (eps_history_count == 3) {
                predicted_eps = exp(4.0 * log_eps_history[0] -
                                    6.0 * log_eps_history[1] +
                                    4.0 * log_eps_history[2]);
            }
            int final_solve = ppml_options.exact_solver || g_ppml_state->tolerance <= 11.0 * irls_tol ||
                              g_ppml_state->tolerance <= target_inner_tol ||
                              predicted_eps <= irls_tol;
            /* Deviance is locally quadratic; covariance depends on weights
             * to first order. Also require stable fitted log-means so an
             * early deviance plateau cannot leave inaccurate standard errors
             * in explicitly tight fits. Preserve default ppmlhdfe stopping
             * behavior; it can intentionally return less accurate weights. */
            if (rel_change < irls_tol && final_solve &&
                (irls_tol >= 1e-10 || max_eta_change <= 1e-8) &&
                (G == 1 || g_ppml_state->tolerance <= target_inner_tol)) {
                ok_iterations++;
                if(ok_iterations>=ppml_options.min_ok) {
                    irls_converged = 1;
                    irls_iter++;
                    break;
                }
            }

            if(ppml_options.step_halving && deviance>deviance_old && halving_count<ppml_options.max_halving) {
                for(i=0;i<N;i++) {
                    double wu=w_user_compact?w_user_compact[i]:1;
                    double old=log(irls_w[i]*user_weight_scale/wu)-(offset_compact?offset_compact[i]:0);
                    eta[i]=ppml_options.halving_memory*old+(1-ppml_options.halving_memory)*eta[i];
                    if(halving_count && eta[i]<-10)eta[i]=-10;
                    mu[i]=exp(eta[i]+(offset_compact?offset_compact[i]:0));
                }
                halving_count++;ok_iterations=0;deviance=deviance_old;
            } else halving_count=0;

            if (rel_change > 0.0 && isfinite(rel_change)) {
                if (eps_history_count < 3) log_eps_history[eps_history_count++] = log(rel_change);
                else {
                    log_eps_history[0] = log_eps_history[1];
                    log_eps_history[1] = log_eps_history[2];
                    log_eps_history[2] = log(rel_change);
                }
            }

            /* Adaptive CG tolerance: tighten as IRLS converges.
             * Matches ppmlhdfe's scheme: new_tol = max(target, eps * 0.01) */
            if (ppml_options.heuristic && rel_change < 0.1) {
                ST_double new_tol = fmax(target_inner_tol, rel_change * 0.01);
                if (new_tol < g_ppml_state->tolerance) {
                    g_ppml_state->tolerance = new_tol;
                }
            }
        } else if (verbose) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                "  iter %d: deviance=%.10e eps=.\n",
                irls_iter + 1, deviance);
            SF_display(msg);
        }

        /* (i) Update IRLS weights for next iteration */
        for (i = 0; i < N; i++) {
            ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
            irls_w[i] = mu[i] * (w_u / user_weight_scale);
            if (irls_w[i] < DBL_MIN) irls_w[i] = DBL_MIN;
        }
        ppml_update_fe_weights(g_ppml_state, irls_w, N);


    }

    if (!irls_converged && !output_rc) {
        SF_error("cpplmhdfe: IRLS did not converge\n");
        output_rc = 430;
    }
    if (output_rc) {
        free(mu); free(eta); free(z_orig); free(aug_copy);
        free(xtx); free(xty); free(irls_beta); free(inv_xx); free(sep_m);
        free(irls_keep_idx); free(xtx_k); free(xty_k); free(beta_k); free(inv_xx_k);
        free(stdev_x); free(w_user_compact); free(offset_compact);
        free(separation_certificate); free(previous_working);
        free(cluster_raw_values); ctools_aligned_free(obs_map); free(data);
        for (g = 0; g < G; g++) {
            free(factors[g].levels); free(factors[g].counts);
        }
        free(factors);
        cleanup_ppml_state(); /* also owns irls_w */
        return output_rc;
    }

    if(ppml_options.check_mu) {
        int nsep=0;for(i=0;i<N;i++)nsep+=sep_m[i];
        if(nsep) {
            sample_var_idx=(int)ppml_option("sample_idx",0);
            for(i=0;i<N && !output_rc;i++) {
                if(!sep_m[i])output_rc=SF_vstore(sample_var_idx,(int)obs_map[i],1);
            }
            ctools_scal_save("__cpplmhdfe_refit",3);
            ctools_scal_save("__cpplmhdfe_num_singletons",num_singletons);
            ctools_scal_save("__cpplmhdfe_num_separated",num_separated+nsep);
            goto ppml_cleanup;
        }
    }

    if (verbose) {
        char msg[256];
        if (irls_converged) {
            snprintf(msg, sizeof(msg), "(IRLS converged in %d iterations, deviance = %.8g)\n",
                     irls_iter, deviance);
        } else {
            snprintf(msg, sizeof(msg), "(IRLS did NOT converge in %d iterations, deviance = %.8g)\n",
                     irls_maxiter, deviance);
        }
        SF_display(msg);
    }

    t_irls = get_time_sec();

    /* ================================================================
     * STEP 10: Final OLS on converged data for beta, VCE
     * ================================================================ */

    /* Ensure final partialled data uses tightest CG tolerance.
     * If the adaptive scheme hasn't tightened to target yet (e.g. fast
     * convergence in few iterations), do one extra pass at target_inner_tol. */
    if (target_inner_tol < g_ppml_state->tolerance) {
        g_ppml_state->tolerance = target_inner_tol;
        memcpy(aug_copy, z_orig, N * sizeof(ST_double));
        if(ppml_options.exact_partial || !irls_iter) {
            for (k = 0; k < K_x; k++)
                memcpy(&aug_copy[(size_t)(k + 1) * N], &X[(size_t)k * N], N * sizeof(ST_double));
        }
        HDFE_SolveResult projection = ppml_partial(g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (projection.status) {
            SF_error("cpplmhdfe: fixed-effect projection did not converge\n");
            output_rc = projection.status;
            goto ppml_cleanup;
        }
    }

    /* Collinearity detection */
    is_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));
    if (!is_collinear) { output_rc = 920; goto ppml_cleanup; }
    ST_int num_collinear = 0;
    ST_int K_keep;

    if (is_collinear) {
        ST_double *xtx_final = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
        if (!xtx_final) { output_rc = 920; goto ppml_cleanup; }
        if (xtx_final) {
            for (i = 0; i < K_x; i++) {
                for (j = 0; j < K_x; j++) {
                    xtx_final[i * K_x + j] = fast_dot(&aug_copy[(size_t)(i+1)*N], &aug_copy[(size_t)(j+1)*N], N);
                }
            }
            detect_collinearity(xtx_final, K_x, is_collinear, verbose);
            free(xtx_final);
        }
        for (k = 0; k < K_x; k++)
            if (is_collinear[k]) num_collinear++;
    }
    K_keep = K_x - num_collinear;

    /* Store collinearity flags */
    ctools_scal_save("__cpplmhdfe_num_collinear", (ST_double)num_collinear);
    for (k = 0; k < K_x; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_collinear_%d", k + 1);
        ctools_scal_save(scalar_name, (ST_double)(is_collinear ? is_collinear[k] : 0));
    }

    /* Build non-collinear index */
    keep_idx = (ST_int *)malloc(K_keep * sizeof(ST_int));
    if (!keep_idx && K_keep > 0) { output_rc = 920; goto ppml_cleanup; }
    if (keep_idx && is_collinear) {
        idx = 0;
        for (k = 0; k < K_x; k++) {
            if (!is_collinear[k])
                keep_idx[idx++] = k;
        }
    }

    /* Final solve: beta and inv(X'WX)
     * Use converged IRLS weights for the final normal equations so stored
     * coefficients match the terminal PPML optimum. */
    beta_final = (ST_double *)calloc(K_keep, sizeof(ST_double));
    if (!beta_final && K_keep > 0) { output_rc = 920; goto ppml_cleanup; }




    ST_int final_use_weights = 1;
    ST_double N_ref = (ST_double)N;
    if (has_weights && weight_type == 2 && w_user_compact) {
        N_ref = 0.0;
        for (i = 0; i < N; i++) N_ref += w_user_compact[i];
    }
    ST_double sum_irls_w = 0.0;
    for (i = 0; i < N; i++) sum_irls_w += irls_w[i];
    ST_double w_scale = (sum_irls_w > 0.0) ? (N_ref / sum_irls_w) : 1.0;
    w_reg = (ST_double *)malloc(N * sizeof(ST_double));
    if (!w_reg) { output_rc = 920; goto ppml_cleanup; }
    if (w_reg) {
        for (i = 0; i < N; i++) w_reg[i] = irls_w[i] * w_scale;
    }

    if (K_keep > 0 && keep_idx && beta_final) {
        xtx_keep = (ST_double *)calloc(K_keep * K_keep, sizeof(ST_double));
        xty_keep = (ST_double *)calloc(K_keep, sizeof(ST_double));
        inv_xx_final = (ST_double *)malloc(K_keep * K_keep * sizeof(ST_double));
        V_keep = (ST_double *)calloc(K_keep * K_keep, sizeof(ST_double));

        if (!xtx_keep || !xty_keep || !inv_xx_final || !V_keep) { output_rc = 920; goto ppml_cleanup; }
        if (xtx_keep && xty_keep && inv_xx_final && V_keep) {

            /* Build X'X (or X'WX for fweight) using dd_real precision. */
            for (i = 0; i < K_keep; i++) {
                for (j = 0; j <= i; j++) {
                    dd_real acc = {0.0, 0.0};
                    const ST_double *xi = &aug_copy[(size_t)(keep_idx[i]+1)*N];
                    const ST_double *xj = &aug_copy[(size_t)(keep_idx[j]+1)*N];
                    for (idx = 0; idx < N; idx++) {
                        ST_double wi = final_use_weights ? (w_reg ? w_reg[idx] : irls_w[idx]) : 1.0;
                        ST_double val = wi * xi[idx];
                        dd_real prod = two_prod(val, xj[idx]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xtx_keep[i * K_keep + j] = acc.hi + acc.lo;
                    xtx_keep[j * K_keep + i] = acc.hi + acc.lo;
                }
                {
                    dd_real acc = {0.0, 0.0};
                    const ST_double *xi = &aug_copy[(size_t)(keep_idx[i]+1)*N];
                    const ST_double *z_part = aug_copy;
                    for (idx = 0; idx < N; idx++) {
                        ST_double wi = final_use_weights ? (w_reg ? w_reg[idx] : irls_w[idx]) : 1.0;
                        ST_double val = wi * xi[idx];
                        dd_real prod = two_prod(val, z_part[idx]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xty_keep[i] = acc.hi + acc.lo;
                }
            }

            /* Solve with iterative refinement for maximum beta precision. */
            memcpy(inv_xx_final, xtx_keep, K_keep * K_keep * sizeof(ST_double));
            if (ctools_cholesky(inv_xx_final, K_keep) != 0) { output_rc = 498; goto ppml_cleanup; }
            {
                ctools_invert_from_cholesky(inv_xx_final, K_keep, inv_xx_final);
                for (i = 0; i < K_keep; i++) {
                    beta_final[i] = 0.0;
                    for (j = 0; j < K_keep; j++)
                        beta_final[i] += inv_xx_final[i * K_keep + j] * xty_keep[j];
                }
                /* Iterative refinement */
                for (ST_int refine = 0; refine < 2; refine++) {
                    ST_double r_keep[64];
                    if (K_keep > 64) break;
                    for (i = 0; i < K_keep; i++) {
                        const ST_double *xi = &aug_copy[(size_t)(keep_idx[i]+1)*N];
                        const ST_double *z_part = aug_copy;
                        dd_real acc = {0.0, 0.0};
                        for (idx = 0; idx < N; idx++) {
                            ST_double wi = final_use_weights ? (w_reg ? w_reg[idx] : irls_w[idx]) : 1.0;
                            ST_double val = wi * xi[idx];
                            dd_real prod = two_prod(val, z_part[idx]);
                            acc = dd_add_d(acc, prod.hi);
                            acc = dd_add_d(acc, prod.lo);
                        }
                        for (j = 0; j < K_keep; j++) {
                            dd_real prod = two_prod(xtx_keep[i * K_keep + j], beta_final[j]);
                            acc = dd_add_d(acc, -prod.hi);
                            acc = dd_add_d(acc, -prod.lo);
                        }
                        r_keep[i] = acc.hi + acc.lo;
                    }
                    for (i = 0; i < K_keep; i++) {
                        ST_double delta = 0.0;
                        for (j = 0; j < K_keep; j++)
                            delta += inv_xx_final[i * K_keep + j] * r_keep[j];
                        beta_final[i] += delta;
                    }
                }
            }
        }
    }

    /* ================================================================
     * STEP 11: VCE computation
     * Matching creghdfe/reghdfe: add weighted means back to X, include
     * constant column, extend inv_xx via block partition formula.
     * This accounts for the intercept in the sandwich VCE.
     * ================================================================ */

    /* Compute IRLS-weighted means of original X columns (before partialling).
     * These means are added back to the partialled X for VCE computation.
     * Uses IRLS weights = mean(x, HDFE.weight) matching ppmlhdfe. */
    ST_int K_with_cons = K_keep + 1;
    means_x = (ST_double *)calloc(K_keep, sizeof(ST_double));
    if (!means_x && K_keep > 0) { output_rc = 920; goto ppml_cleanup; }
    if (means_x && K_keep > 0) {
        ST_double sum_w_eff = 0.0;
        for (i = 0; i < N; i++) sum_w_eff += irls_w[i];
        for (k = 0; k < K_keep; k++) {
            ST_double wmean = 0.0;
            const ST_double *xk = &X[(size_t)keep_idx[k] * N];
            for (i = 0; i < N; i++) {
                wmean += irls_w[i] * xk[i];
            }
            means_x[k] = wmean / sum_w_eff;
        }
    }

    /* Extend inv_xx_final to (K_keep+1)x(K_keep+1) using block partition formula.
     * The (K+1)th element is the constant/intercept.
     * inv([A b; b' c]) with A=X'WX, b=X'W1, c=1'W1.
     * Using: side = -inv(A) * mean_x, corner = 1/c - mean_x' * side */
    inv_xx_ext = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    if (!inv_xx_ext) { output_rc = 920; goto ppml_cleanup; }
    if (inv_xx_ext && inv_xx_final && means_x && K_keep > 0) {
        ST_double *side = (ST_double *)calloc(K_keep, sizeof(ST_double));
        if (!side) { output_rc = 920; goto ppml_cleanup; }
        if (side) {
            for (j = 0; j < K_keep; j++) {
                for (i = 0; i < K_keep; i++) {
                    side[j] -= means_x[i] * inv_xx_final[i * K_keep + j];
                }
            }

            /* N_corner = 1'W1 where W is the converged IRLS weight matrix. */
            ST_double N_corner = N_ref;
            ST_double corner = 1.0 / N_corner;
            for (i = 0; i < K_keep; i++) {
                corner -= means_x[i] * side[i];
            }

            /* Build extended matrix */
            for (i = 0; i < K_keep; i++) {
                for (j = 0; j < K_keep; j++) {
                    inv_xx_ext[i * K_with_cons + j] = inv_xx_final[i * K_keep + j];
                }
                inv_xx_ext[i * K_with_cons + K_keep] = side[i];
                inv_xx_ext[K_keep * K_with_cons + i] = side[i];
            }
            inv_xx_ext[K_keep * K_with_cons + K_keep] = corner;

            free(side);
        }
    }

    if (!K_keep) inv_xx_ext[0] = 1.0 / N_ref;

    /* Build X_eff with means added back + constant column (K_with_cons columns) */

    {
        X_eff = (ST_double *)malloc((size_t)N * K_with_cons * sizeof(ST_double));
        if (!X_eff) { output_rc = 920; goto ppml_cleanup; }
        if (X_eff) {
            for (k = 0; k < K_keep; k++) {
                const ST_double *xp = &aug_copy[(size_t)(keep_idx[k] + 1) * N];
                ST_double mk = means_x ? means_x[k] : 0.0;
                for (i = 0; i < N; i++) {
                    X_eff[(size_t)k * N + i] = xp[i] + mk;
                }
            }
            /* Constant column */
            for (i = 0; i < N; i++) {
                X_eff[(size_t)K_keep * N + i] = 1.0;
            }
        }
    }

    /* Compute final OLS residual from partialled data. */
    vce_resid = (ST_double *)malloc(N * sizeof(ST_double));
    if (!vce_resid) { output_rc = 920; goto ppml_cleanup; }
    if (vce_resid) {
        if (beta_final && K_keep > 0) {
            for (i = 0; i < N; i++) {
                ST_double xb = 0.0;
                for (k = 0; k < K_keep; k++) {
                    xb += aug_copy[(size_t)(keep_idx[k] + 1) * N + i] * beta_final[k];
                }
                vce_resid[i] = aug_copy[i] - xb;
            }
        }
        else {
            for (i = 0; i < N; i++)
                vce_resid[i] = aug_copy[i];
        }
    }

    /* Count each marginal cluster and detect FE nesting on the retained sample.
     * Subset indices are bitmasks; powers of two are marginal dimensions. */
    ST_int nested[10] = {0};
    if (vcetype == 2) {
        cluster_ids = malloc(N * sizeof(ST_int));
        if (!cluster_ids) { output_rc = 920; goto ppml_cleanup; }
        num_clusters = N;
        ST_int dim = 0;
        for (ST_int subset = 1; subset <= cluster_terms; subset <<= 1) {
            ST_int ncl;
            if (ctools_numeric_to_cluster_ids(
                    cluster_raw_values + (size_t)(subset - 1) * N_orig,
                    N, cluster_ids, &ncl)) {
                output_rc = 920; goto ppml_cleanup;
            }
            if (ncl < num_clusters) num_clusters = ncl;
            snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_N_clust_%d", ++dim);
            ctools_scal_save(scalar_name, ncl);
            for (g = 0; g < G; g++) {
                if (slope_only[g] || !(dof_flags & 4)) continue;
                ST_int is_nested = ctools_fe_nested_in_cluster(
                    g_ppml_state->factors[g].levels,
                    g_ppml_state->factors[g].num_levels, cluster_ids, N);
                if (is_nested < 0) { output_rc = 920; goto ppml_cleanup; }
                nested[g] |= is_nested;
            }
        }
        if (num_clusters < 2) {
            SF_error("cpplmhdfe: each cluster dimension requires at least two retained clusters\n");
            output_rc = 459; goto ppml_cleanup;
        }
    }
    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_fe_nested_%d", g + 1);
        ctools_scal_save(scalar_name, nested[g]);
    }

    /* Match reghdfe's pairwise mobility and cluster-nesting accounting.
     * These degrees of freedom describe the model; PPML's asymptotic
     * sandwich uses only N/(N-1) or G/(G-1), not residual-DF scaling. */
    ST_int redundant[10] = {0}, exact[10] = {0};
    ST_int df_initial = 0, df_nested = 0, df_redundant = 0, first_intercept = -1, pairs = 0;
    for (g = 0; g < G; g++) if (nested[g]) df_nested += factors[g].num_levels;
    for (g = 0; g < G; g++) {
        FE_Factor *fg = &g_ppml_state->factors[g];
        df_initial += fg->num_levels;
        if (nested[g]) { redundant[g] = fg->num_levels; exact[g] = 1; continue; }
        if (slope_only[g]) {
            if (dof_flags & 8) {
                ST_double *lo = malloc((size_t)fg->num_levels*sizeof(*lo));
                ST_double *hi = malloc((size_t)fg->num_levels*sizeof(*hi));
                if (!lo || !hi) { free(lo); free(hi); output_rc = 920; goto ppml_cleanup; }
                for (i = 0; i < fg->num_levels; i++) { lo[i] = DBL_MAX; hi[i] = -DBL_MAX; }
                for (i = 0; i < N; i++) {
                    ST_int l = fg->levels[i]-1;
                    lo[l] = fmin(lo[l], fg->slope[i]); hi[l] = fmax(hi[l], fg->slope[i]);
                }
                ST_int has_int = 0;
                for (j = 0; j < G; j++) {
                    FE_Factor *fj = &g_ppml_state->factors[j];
                    if (!slope_only[j] && fg->num_levels == fj->num_levels &&
                        !memcmp(fg->levels, fj->levels, (size_t)N*sizeof(ST_int))) has_int = 1;
                }
                for (i = 0; i < fg->num_levels; i++)
                    redundant[g] += has_int ? hi[i]-lo[i] <= 1e-6 : fabs(lo[i])+fabs(hi[i]) <= 1e-6;
                free(lo); free(hi);
            }
            continue;
        }
        redundant[g] = first_intercept >= 0 || df_nested > 0;
        if (first_intercept < 0) { first_intercept = g; exact[g] = 1; }
        if (!(dof_flags & 3)) continue;
        for (j = 0; j < g; j++) {
            if (slope_only[j] || nested[j]) continue;
            if ((dof_flags & 2) && pairs) break;
            ST_int groups = count_connected_components(g_ppml_state->factors[j].levels,
                fg->levels, N, g_ppml_state->factors[j].num_levels, fg->num_levels);
            if (groups < 0) { output_rc = 920; goto ppml_cleanup; }
            if (!pairs) exact[g] = 1;
            pairs++;
            if (groups > redundant[g]) redundant[g] = groups;
        }
    }
    for (g = 0; g < G; g++) {
        df_redundant += redundant[g];
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_fe_redundant_%d", g+1);
        ctools_scal_save(scalar_name, redundant[g]);
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_fe_exact_%d", g+1);
        ctools_scal_save(scalar_name, exact[g]);
    }
    df_a = df_initial - df_redundant;
    ctools_scal_save("__cpplmhdfe_df_initial", df_initial);
    ctools_scal_save("__cpplmhdfe_df_nested", df_nested);
    ctools_scal_save("__cpplmhdfe_df_redundant", df_redundant);

    /* Compute VCE: PPML sandwich estimator matching creghdfe/reghdfe.
     * Uses extended (K_keep+1)x(K_keep+1) system with constant column,
     * then extracts K_keep×K_keep submatrix for reporting.
     *
     * For fweight: D = extended H^{-1} (raw), resid = (y-mu), weights = fw
     * For others: D extended & normalized, resid = OLS_resid, weights = irls_w as aweight */
  /* (K_with_cons)x(K_with_cons) full VCE */
    if (vce_resid && X_eff && inv_xx_ext) {
        V_ext = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
        if (!V_ext) { output_rc = 920; goto ppml_cleanup; }

        if (vcetype == 0) {
            /* Unadjusted VCE: V = H^{-1} (Poisson dispersion = 1) */
            for (i = 0; i < K_keep; i++)
                for (j = 0; j < K_keep; j++)
                    V_keep[i * K_keep + j] = inv_xx_final[i * K_keep + j];

        } else if (V_ext && ((vcetype == 1) || (vcetype == 2 && cluster_ids))) {

            ctools_vce_data d;
            /* For a single sandwich, the centered slope system is equivalent
             * and avoids cancellation from large regressor means. Multiway
             * PSD correction needs the full system including the constant. */
            ST_int vce_K = K_with_cons;
            if (cluster_terms <= 1) {
                for (k = 0; k < K_keep; k++)
                    memcpy(X_eff + (size_t)k * N,
                           aug_copy + (size_t)(keep_idx[k] + 1) * N,
                           N * sizeof(ST_double));
            }
            if (cluster_terms <= 1) {
                for (i = 0; i < K_keep; i++) {
                    inv_xx_ext[i*K_with_cons+K_keep] = 0;
                    inv_xx_ext[K_keep*K_with_cons+i] = 0;
                }
                inv_xx_ext[K_keep*K_with_cons+K_keep] = 1.0/N_ref;
            }
            d.X_eff = X_eff;
            d.D = inv_xx_ext;
            ST_double *vce_resid_work = vce_resid;
            if (has_weights && weight_type == 2 && w_user_compact) {
                vce_resid_work = (ST_double *)malloc(N * sizeof(ST_double));
                if (!vce_resid_work) { output_rc = 920; goto ppml_cleanup; }
                memcpy(vce_resid_work, vce_resid, N * sizeof(ST_double));
                for (i = 0; i < N; i++) {
                    ST_double fw = w_user_compact[i];
                    ST_double wr = w_reg ? w_reg[i] : irls_w[i];
                    if (fw > 0.0 && vce_resid_work) {
                        vce_resid_work[i] *= wr / fw;
                    }
                }
                d.weights = w_user_compact;
                d.weight_type = 2;
                d.normalize_weights = 0;
            }
            else {
                d.weights = w_reg ? w_reg : irls_w;
                d.weight_type = 1;
                d.normalize_weights = 1;
            }
            d.resid = vce_resid_work ? vce_resid_work : vce_resid;
            d.N = N;
            d.K = vce_K;

            if (vcetype == 1) {
                ST_double dof_adj = N_ref / (N_ref - 1.0);
                output_rc = ctools_vce_robust(&d, dof_adj, V_ext);
            } else {
                ST_double dof_adj = (ST_double)num_clusters / (num_clusters - 1);
                ST_double *V_term = calloc((size_t)K_with_cons * K_with_cons, sizeof(ST_double));
                if (!V_term) output_rc = 920;
                for (ST_int subset = 1; subset <= cluster_terms && !output_rc; subset++) {
                    ST_int ncl, bits = 0;
                    for (ST_int bit = subset; bit; bit >>= 1) bits += bit & 1;
                    if (ctools_numeric_to_cluster_ids(
                            cluster_raw_values + (size_t)(subset - 1) * N_orig,
                            N, cluster_ids, &ncl)) {
                        output_rc = 920; break;
                    }
                    output_rc = ctools_vce_cluster(&d, cluster_ids, ncl, dof_adj, V_term);
                    ST_double sign = bits % 2 ? 1.0 : -1.0;
                    for (i = 0; i < K_with_cons * K_with_cons && !output_rc; i++)
                        V_ext[i] += sign * V_term[i];
                }
                free(V_term);
            }
            if (vce_resid_work && vce_resid_work != vce_resid) free(vce_resid_work);

            if (output_rc) {
                SF_error("cpplmhdfe: covariance calculation failed\n");
                goto ppml_cleanup;
            }

            /* Centered-system sandwich is numerically stable even when a
             * regressor has a huge mean. Transform only the constant row and
             * column back to the reported parameterization. */
            if (cluster_terms <= 1) {
                ST_double corner = V_ext[K_keep*K_with_cons+K_keep];
                for (i = 0; i < K_keep; i++) {
                    corner -= 2*means_x[i]*V_ext[i*K_with_cons+K_keep];
                    for (j = 0; j < K_keep; j++)
                        corner += means_x[i]*means_x[j]*V_ext[i*K_with_cons+j];
                }
                for (i = 0; i < K_keep; i++) {
                    ST_double side = V_ext[i*K_with_cons+K_keep];
                    for (j = 0; j < K_keep; j++) side -= V_ext[i*K_with_cons+j]*means_x[j];
                    V_ext[i*K_with_cons+K_keep] = V_ext[K_keep*K_with_cons+i] = side;
                }
                V_ext[K_keep*K_with_cons+K_keep] = corner;
            }
            for (i = 0; i < K_with_cons; i++) {
                ctools_mat_store("__cpplmhdfe_scales", 1, i+1,
                                 i < K_keep ? stdev_x[keep_idx[i]] : 1.0);
                for (j = 0; j < K_with_cons; j++)
                    ctools_mat_store("__cpplmhdfe_V_full", i+1, j+1, V_ext[i*K_with_cons+j]);
            }

            /* Extract K_keep×K_keep submatrix (dropping constant row/col) */
            for (i = 0; i < K_keep; i++)
                for (j = 0; j < K_keep; j++)
                    V_keep[i * K_keep + j] = V_ext[i * vce_K + j];
        }
    }
    free(V_ext); V_ext = NULL;

    t_vce = get_time_sec();

    /* ================================================================
     * STEP 12: Un-standardize and store results
     * ================================================================ */

    ST_double constant = log(stdev_y), mean_d = 0;
    for (i = 0; i < N; i++) {
        ST_double off = offset_compact ? offset_compact[i] : 0;
        constant += irls_w[i] / sum_irls_w * (log(mu[i]) - off);
        ST_double d = z_orig[i] - vce_resid[i];
        for (k = 0; k < K_keep; k++) d -= X[(size_t)keep_idx[k]*N+i]*beta_final[k];
        mean_d += irls_w[i] / sum_irls_w * d;
    }
    for (k = 0; k < K_keep; k++) constant -= means_x[k]*beta_final[k];
    ctools_scal_save("__cpplmhdfe_cons", constant);
    ST_int d_var_idx = 0;
    if (SF_scal_use("__cpplmhdfe_d_idx", &val) == 0) d_var_idx = (ST_int)val;
    if (d_var_idx > 0) {
        for (i = 0; i < N && !output_rc; i++) {
            ST_double d = z_orig[i] - vce_resid[i] - mean_d;
            for (k = 0; k < K_keep; k++) d -= X[(size_t)keep_idx[k]*N+i]*beta_final[k];
            output_rc = SF_vstore(d_var_idx, (ST_int)obs_map[i], d);
        }
    }

    for(i=0;i<N;i++) {
        separation_certificate[i]=z_orig[i]-vce_resid[i]-mean_d;
        for(k=0;k<K_keep;k++)separation_certificate[i]-=X[(size_t)keep_idx[k]*N+i]*beta_final[k];
    }
    output_rc=ppml_save_effects(g_ppml_state,separation_certificate,obs_map);
    if(output_rc)goto ppml_cleanup;

    /* Un-standardize y and mu for log-likelihood computation.
     * y_actual = y_std * stdev_y, mu_actual = mu_std * stdev_y */
    for (i = 0; i < N; i++) {
        y[i] *= stdev_y;
        mu[i] *= stdev_y;
    }

    /* Un-standardize beta: beta_actual = beta_std / stdev_x */
    if (stdev_x && keep_idx && beta_final) {
        for (k = 0; k < K_keep; k++)
            beta_final[k] /= stdev_x[keep_idx[k]];
    }

    /* Un-standardize V: V_actual[i,j] = V_std[i,j] / (stdev_x[i] * stdev_x[j]) */
    if (stdev_x && keep_idx && V_keep) {
        for (i = 0; i < K_keep; i++)
            for (j = 0; j < K_keep; j++)
                V_keep[i * K_keep + j] /= stdev_x[keep_idx[i]] * stdev_x[keep_idx[j]];
    }

    ST_double ll = ppml_compute_loglik(y, mu, w_user_compact, N);

    /* Compute log-likelihood of null model (intercept only): ll_0 */
    ST_double y_mean = 0.0, w_sum = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        y_mean += w_u * y[i];
        w_sum += w_u;
    }
    y_mean /= w_sum;
    if (y_mean < 1e-18) y_mean = 1e-18;
    ST_double ll_0 = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        if (y[i] > 0.0) {
            ll_0 += w_u * (y[i] * log(y_mean) - y_mean - lgamma(y[i] + 1.0));
        } else {
            ll_0 += w_u * (-y_mean);
        }
    }

    if (SF_scal_use("__cpplmhdfe_sample_idx", &val) == 0)
        sample_var_idx = (ST_int)val;
    if (sample_var_idx > 0) {
        for (i = 0; i < N && !output_rc; i++)
            output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 1.0);
    }

    /* N reporting: for fweights, report sum of weights; otherwise report obs count */
    {
        ST_double N_report = (ST_double)N;
        if (weight_type == 2 && has_weights && w_user_compact) {
            N_report = 0.0;
            for (i = 0; i < N; i++) N_report += w_user_compact[i];
        }
        ctools_scal_save("__cpplmhdfe_N", N_report);
    }

    ctools_scal_save("__cpplmhdfe_num_singletons", (ST_double)num_singletons);
    ctools_scal_save("__cpplmhdfe_num_separated", (ST_double)num_separated);
    ctools_scal_save("__cpplmhdfe_K_keep", (ST_double)K_keep);
    ctools_scal_save("__cpplmhdfe_df_a", (ST_double)df_a);
    ctools_scal_save("__cpplmhdfe_mobility_groups", (ST_double)mobility_groups);
    ctools_scal_save("__cpplmhdfe_deviance", deviance * stdev_y * user_weight_scale);
    ctools_scal_save("__cpplmhdfe_ll", ll);
    ctools_scal_save("__cpplmhdfe_ll_0", ll_0);
    ctools_scal_save("__cpplmhdfe_irls_iterations", (ST_double)irls_iter);
    ctools_scal_save("__cpplmhdfe_irls_converged", (ST_double)irls_converged);

    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_num_levels_%d", g + 1);
        ctools_scal_save(scalar_name, (ST_double)factors[g].num_levels);
    }

    /* Store betas */
    for (k = 0; k < K_keep; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_beta_%d", k + 1);
        ctools_scal_save(scalar_name, beta_final ? beta_final[k] : 0.0);
    }

    /* Store VCE matrix */
    if (V_keep) {
        for (i = 0; i < K_keep; i++) {
            for (j = 0; j < K_keep; j++) {
                ctools_mat_store("__cpplmhdfe_V", i + 1, j + 1, V_keep[i * K_keep + j]);
            }
        }
    }


    ctools_scal_save("__cpplmhdfe_cluster_terms_done", cluster_terms);

    /* Store cluster count */
    if (vcetype == 2 && num_clusters > 0) {
        ctools_scal_save("__cpplmhdfe_N_clust", (ST_double)num_clusters);
    }

    /* Timing */
    ctools_scal_save("_cpplmhdfe_time_load", t_load - t_start);
    ctools_scal_save("_cpplmhdfe_time_remap", t_remap - t_load);
    ctools_scal_save("_cpplmhdfe_time_singleton", t_singleton - t_remap);
    ctools_scal_save("_cpplmhdfe_time_dof", t_dof - t_singleton);
    ctools_scal_save("_cpplmhdfe_time_irls", t_irls - t_dof);
    ctools_scal_save("_cpplmhdfe_time_vce", t_vce - t_irls);
    ctools_scal_save("_cpplmhdfe_time_total", t_vce - t_start);
    CTOOLS_SAVE_THREAD_INFO("_cpplmhdfe");

    /* ================================================================
     * Cleanup
     * ================================================================ */
ppml_cleanup:
    free(V_ext);
    free(separation_certificate);
    free(previous_working);
    free(cluster_raw_values);
    free(mu); free(eta);
    free(z_orig); free(aug_copy);
    free(xtx); free(xty); free(irls_beta); free(inv_xx); free(sep_m);
    if (irls_keep_idx) free(irls_keep_idx);
    if (xtx_k) free(xtx_k);
    if (xty_k) free(xty_k);
    if (beta_k) free(beta_k);
    if (inv_xx_k) free(inv_xx_k);
    if (vce_resid) free(vce_resid);
    if (beta_final) free(beta_final);
    if (inv_xx_final) free(inv_xx_final);
    if (V_keep) free(V_keep);
    if (xtx_keep) free(xtx_keep);
    if (xty_keep) free(xty_keep);
    if (w_reg) free(w_reg);
    if (is_collinear) free(is_collinear);
    if (keep_idx) free(keep_idx);
    if (X_eff) free(X_eff);
    if (means_x) free(means_x);
    if (inv_xx_ext) free(inv_xx_ext);
    if (cluster_ids) free(cluster_ids);
    if (stdev_x) free(stdev_x);
    if (w_user_compact) free(w_user_compact);
    if (offset_compact) free(offset_compact);

    for (g = 0; g < G; g++) {
        free(factors[g].levels);
        if (factors[g].counts) free(factors[g].counts);
    }
    free(factors);
    if (mask) free(mask);
    free(data);
    if (obs_map) ctools_aligned_free(obs_map);

    /* State owns irls_w; cleanup releases it together with the FE buffers. */
    cleanup_ppml_state();

    return output_rc;
}

/* ========================================================================
 * Helper: Update FE weighted counts from IRLS weights
 * ======================================================================== */
static void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N)
{
    ST_int g, i, lev;
    ST_double sum_w = 0.0;

    for (i = 0; i < N; i++)
        sum_w += irls_w[i];
    S->sum_weights = sum_w;

    for (g = 0; g < S->G; g++) {
        FE_Factor *f = &S->factors[g];
        if (f->slope_center) {
            memset(f->weighted_counts, 0, f->num_levels*sizeof(ST_double));
            memset(f->slope_center, 0, f->num_levels*sizeof(ST_double));
            for (i = 0; i < N; i++) {
                lev = f->levels[i]-1;
                f->weighted_counts[lev] += irls_w[i];
                f->slope_center[lev] += irls_w[i]*f->slope[i];
            }
            for (i = 0; i < f->num_levels; i++)
                if (f->weighted_counts[i] > 0) f->slope_center[i] /= f->weighted_counts[i];
        }
        memset(f->weighted_counts, 0, f->num_levels * sizeof(ST_double));
        for (i = 0; i < N; i++) {
            lev = f->levels[i] - 1;
            ST_double loading = f->slope ? f->slope[i] - (f->slope_center ? f->slope_center[lev] : 0.0) : 1.0;
            f->weighted_counts[lev] += irls_w[i]*loading*loading;
        }
        if (f->inv_weighted_counts) {
            for (i = 0; i < f->num_levels; i++) {
                f->inv_weighted_counts[i] = (f->weighted_counts[i] > 0.0) ?
                    1.0 / f->weighted_counts[i] : 0.0;
            }
        }
    }
}

/* ========================================================================
 * Helper: Compute Poisson deviance
 * D = 2 * sum[y*log(y/mu) - (y - mu)]  (y*log(y/mu) = 0 when y = 0)
 * ======================================================================== */
static ST_double ppml_compute_deviance(const ST_double *y, const ST_double *mu,
                                        const ST_double *w_user, ST_int N)
{
    ST_int i;
    ST_double dev = 0.0;

    for (i = 0; i < N; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double d;
        if (y[i] > 0.0) {
            d = y[i] * log(y[i] / mu[i]) - (y[i] - mu[i]);
        } else {
            d = mu[i];  /* 0*log(0/mu) - (0 - mu) = mu */
        }
        dev += w * d;
    }

    return 2.0 * dev;
}

/* ========================================================================
 * Helper: Compute Poisson log-likelihood
 * ll = sum[y*log(mu) - mu - log(y!)]
 * ======================================================================== */
static ST_double ppml_compute_loglik(const ST_double *y, const ST_double *mu,
                                      const ST_double *w_user, ST_int N)
{
    ST_int i;
    ST_double ll = 0.0;

    for (i = 0; i < N; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double contrib;
        if (y[i] > 0.0) {
            contrib = y[i] * log(mu[i]) - mu[i] - lgamma(y[i] + 1.0);
        } else {
            contrib = -mu[i];
        }
        ll += w * contrib;
    }

    return ll;
}
