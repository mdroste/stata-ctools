/* Independent binscatter transformation oracles and allocation sweeps. */
#include "newcommands_mock.h"
#include "cbinscatter/cbinscatter_bins.c"
#include "cbinscatter/cbinscatter_fit.c"
#include "cbinscatter/cbinscatter_resid.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "ctools_hdfe_utils.c"
#include "ctools_ols.c"
static ST_boolean ismissing(double x){return x>=mock.missval;}
static ST_retcode message(char *s){(void)s;return 0;}
static int last_count;
static void reset(int failure)
{
    assert(!live);allocations=0;fail_alloc=failure;
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(threads);
#endif
}
static void finish(ByGroupResult *g)
{
    last_count=allocations;free(g->bins);free(g->fit_coefs);assert(!live);
}
static int keys_case(int weighted,int medians,int failure)
{
    reset(failure);ByGroupResult g={0};
    const double labels[]={1,2147483647.,2147483648.,2147483649.,0x1p40,0x1p52};
    double key[25],x[25],y[25],w[25];int ids[25];
    /* Interleaved bins: equal half-weights at the central gap yield medians
     * x=5 and y=10; weighted means remain independently computable. */
    for(int i=0;i<24;i++) {
        int b=i%6,r=i/6;key[i]=labels[b];x[i]=b+(double[]){0,2,8,20}[r];
        y[i]=2*x[i];w[i]=(double[]){.1,.2,.1,.2}[r];
    }
    key[24]=mock.missval;x[24]=y[24]=999;w[24]=1;
    int rc=compute_bins_by_key(y,x,key,weighted?w:NULL,25,1,6,&g,ids);
    if(!rc) {
        assert(g.num_bins==6&&ids[24]==0);
        set_external_bin_ids(key,ids,25,&g);
        if(medians)rc=apply_bin_medians(y,x,weighted?w:NULL,25,ids,1,&g);
        if(!rc)for(int b=0;b<6;b++) {
            assert(g.bins[b].bin_id==b+1&&g.bins[b].external_id==labels[b]);
            assert(g.bins[b].n_obs==4);
            double expected=b+(medians?5:weighted?26./3:7.5);
            if(fabs(g.bins[b].y_mean-2*expected)>=1e-12)
                fprintf(stderr,"bin %d weighted=%d medians=%d allocation=%d y=%.17g expected=%.17g\n",b,weighted,medians,failure,g.bins[b].y_mean,2*expected);
            assert(fabs(g.bins[b].x_mean-expected)<1e-12);
            assert(fabs(g.bins[b].y_mean-2*expected)<1e-12);
            assert(fabs(g.bins[b].sum_weights-(weighted?.6:4))<1e-12);
        }
        if(!rc)for(int i=0;i<24;i++)assert(ids[i]==i%6+1);
    }
    finish(&g);return rc;
}
static int sampled_case(int weighted,int failure)
{
    reset(failure);
    double x[]={8,100,0,200,2,300,20,400},w[]={.1,1,.1,1,.2,1,.2,1};
    unsigned char sel[]={1,0,1,0,1,0,1,0};double *cut=NULL;
    int rc=compute_sampled_cutpoints(x,weighted?w:NULL,sel,8,4,&cut);
    if(!rc) {
        const double unweighted[]={0,1,5,14,20},weighted_cut[]={0,2,5,20,20};
        for(int i=0;i<5;i++)assert(cut[i]==(weighted?weighted_cut:unweighted)[i]);
    }
    last_count=allocations;free(cut);assert(!live);return rc;
}
static int fit_case(int weighted,int order,int controls,int absorb,int many,int failure)
{
    reset(0);int n=many?196608:48;
    /* Allocate input outside the production allocation-failure sequence. */
    double *x=malloc(n*sizeof(double)),*y=malloc(n*sizeof(double)),*w=malloc(n*sizeof(double));
    double *c=malloc(2*n*sizeof(double));int *fe=malloc(n*sizeof(int));
    assert(x&&y&&w&&c&&fe);int input_live=live;allocations=0;fail_alloc=failure;
    double cmean=0,fmean=0,sw=0;
    for(int i=0;i<n;i++) {
        x[i]=(i%4)-1;w[i]=weighted?1+i%3:1;c[i]=(i/4)%3;c[n+i]=2*c[i];
        fe[i]=1+(i/12)%2;
        y[i]=3+2*x[i]+(order>1?.5*x[i]*x[i]:0)+(controls?4*c[i]:0)+(absorb?5*fe[i]:0);
        sw+=w[i];cmean+=w[i]*c[i];fmean+=w[i]*fe[i];
    }
    cmean/=sw;fmean/=sw;
    BinscatterConfig config={.linetype=order,.maxiter=10000,.tolerance=1e-12,.weight_type=weighted};
    ByGroupResult g={0};
    int rc=fit_group_polynomial(y,x,weighted?w:NULL,n,controls?c:NULL,controls?2:0,
                                absorb?fe:NULL,absorb,&config,&g);
    if(!rc) {
        for(int j=0;j<4;j++) {
            double xx=j-1,got=0;for(int k=g.fit_order;k>=0;k--)got=got*xx+g.fit_coefs[k];
            double expected=3+2*xx+(order>1?.5*xx*xx:0)+(controls?4*cmean:0)+(absorb?5*fmean:0);
            assert(fabs(got-expected)<1e-8);
        }
        assert(g.fit_n==n);
    }
    last_count=allocations;free(g.fit_coefs);assert(live==input_live);
    free(x);free(y);free(w);free(c);free(fe);assert(!live);return rc;
}
static int few_case(int order,int weighted,int failure)
{
    reset(failure);double x[]={1,2,1,2,1,2},y[]={1,10,3,12,5,14},w[]={1,1,2,2,1,1};
    BinscatterConfig config={.linetype=order};ByGroupResult g={0};
    int rc=fit_group_polynomial(y,x,weighted?w:NULL,6,NULL,0,NULL,0,&config,&g);
    if(!rc) {
        assert(g.fit_order==order&&g.fit_n==6);
        for(int j=0;j<2;j++) { double xx=j+1,value=0;
            for(int k=order;k>=0;k--)value=value*xx+g.fit_coefs[k];
            assert(fabs(value-(j?12:3))<1e-12);
        }
        for(int k=1;k<order;k++)assert(g.fit_coefs[k]==0);
    }
    finish(&g);return rc;
}
int main(void)
{
    mock.missval=8.98846567431158e307;mock.ismissing=ismissing;mock.spouterr=mock.spoutsml=message;
    for(threads=1;threads<=4;threads*=4)for(int weighted=0;weighted<2;weighted++) {
        for(int medians=0;medians<2;medians++) {
            assert(!keys_case(weighted,medians,0));int n=last_count;
            for(int f=1;f<=n;f++)assert(keys_case(weighted,medians,f)==920);
        }
        assert(!sampled_case(weighted,0));int n=last_count;
        for(int f=1;f<=n;f++)assert(sampled_case(weighted,f)==920);
        for(int order=2;order<=3;order++) {
            assert(!few_case(order,weighted,0));n=last_count;
            for(int f=1;f<=n;f++)assert(few_case(order,weighted,f)==920);
        }
        for(int ctrl=0;ctrl<2;ctrl++)for(int absorb=0;absorb<2;absorb++) {
            assert(!fit_case(weighted,2,ctrl,absorb,0,0));n=last_count;
            for(int f=1;f<=n;f++)assert(fit_case(weighted,2,ctrl,absorb,0,f)==920);
            assert(!fit_case(weighted,2,ctrl,absorb,1,0));
        }
    }
    reset(0);double x[]={1,2},*cut=NULL;unsigned char sel[]={0,0};
    for(int nq=1;nq<=3;nq++)assert(compute_sampled_cutpoints(x,NULL,sel,2,nq,&cut)==198&&!cut&&!live);
    sel[0]=1;assert(!compute_sampled_cutpoints(x,NULL,sel,2,2,&cut));
    assert(cut[0]==1&&cut[1]==1&&cut[2]==1);free(cut);assert(!live);
    puts("PASS binscatter keys, medians, subsamples, polynomial fits, controls, absorption and faults");return 0;
}
