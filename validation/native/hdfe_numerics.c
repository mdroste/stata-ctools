/* Regression for a weak weighted FE bridge and a healthy tiny CG step.
 * The old absolute denominator floor incorrectly accepts the former. */
#include <assert.h>
#include "hdfe/hdfe_solver.c"

static void run_case(double scale, double weight_scale, int weak, int weighted, int threads)
{
    int a[6]={1,1,2,2,1,1}, b[6]={2,2,2,2,1,1};
    double w[6], y[6], expected[6], counts[2][2]={{0}}, inv[2][2];
    double means[4][2][2], r[4][6], u[4][6], v[4][6];
    double *rp[4], *up[4], *vp[4], *mp[8];
    FE_Factor factors[2]={{0}};
    HDFE_State S={0};
    for (int i=0;i<6;i++) {
        double sign=(i%2)?1:-1;
        w[i]=(weak && a[i]!=b[i]?1e-14:2)*weight_scale;
        expected[i]=sign*scale;
        y[i]=(sign+(weak?1e6*(a[i]-b[i]):3*a[i]-2*b[i]))*scale;
        counts[0][a[i]-1]+=weighted?w[i]:1;
        counts[1][b[i]-1]+=weighted?w[i]:1;
    }
    for (int j=0;j<2;j++) {
        factors[j].num_levels=2; factors[j].levels=j?b:a;
        factors[j].counts=counts[j]; factors[j].inv_counts=inv[j];
        factors[j].weighted_counts=counts[j]; factors[j].inv_weighted_counts=inv[j];
        for (int k=0;k<2;k++) inv[j][k]=1/counts[j][k];
    }
    for (int j=0;j<threads;j++) {
        rp[j]=r[j]; up[j]=u[j]; vp[j]=v[j];
        mp[2*j]=means[j][0]; mp[2*j+1]=means[j][1];
    }
    S.G=2; S.N=6; S.K=1; S.has_weights=weighted;
    S.weights=weighted?w:NULL; S.factors=factors;
    S.maxiter=1000; S.tolerance=1e-12; S.num_threads=threads;
    S.thread_cg_r=rp; S.thread_cg_u=up; S.thread_cg_v=vp; S.thread_fe_means=mp;
    HDFE_SolveResult result=partial_out_columns(&S,y,6,1,threads);
    if (weak) {
        if (result.status==0) {
            assert(result.converged);
            for (int i=0;i<6;i++) assert(fabs((y[i]-expected[i])/scale)<1e-8);
        } else assert(result.status==498 && !result.converged);
    } else {
        assert(result.status==0 && result.converged);
        for (int i=0;i<6;i++) assert(fabs((y[i]-expected[i])/scale)<1e-10);
        /* Exact zero residual and idempotent input must not be rejected. */
        result=partial_out_columns(&S,y,6,1,threads);
        assert(result.status==0 && result.converged);
        for (int i=0;i<6;i++) y[i]=(3*a[i]-2*b[i])*scale;
        result=partial_out_columns(&S,y,6,1,threads);
        assert(result.status==0 && result.converged);
        for (int i=0;i<6;i++) assert(fabs(y[i]/scale)<1e-12);
        memset(y,0,sizeof(y));
        result=partial_out_columns(&S,y,6,1,threads);
        assert(result.status==0 && result.converged);
    }
}

int main(void)
{
    for (int threads=1;threads<=4;threads*=4)
        for (int exp=-40;exp<=40;exp+=40)
            for (int wexp=-10;wexp<=10;wexp+=10) {
                run_case(ldexp(1,exp),ldexp(1,wexp),1,1,threads);
                run_case(ldexp(1,exp),ldexp(1,wexp),0,1,threads);
                run_case(ldexp(1,exp),1,0,0,threads);
            }
    puts("PASS HDFE numerical breakdown and scale invariance");
    return 0;
}
