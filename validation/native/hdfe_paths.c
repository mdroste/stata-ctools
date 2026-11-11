/* Analytic paired residuals are orthogonal to every supplied slope/category. */
#include "newcommands_mock.h"
void ctools_msg(const char *m,const char *f,...){(void)m;(void)f;}
#include "ctools_parse.c"
#include "ctools_order.c"
#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "hdfe/hdfe_solver.c"
#include "chdfe/chdfe_impl.c"
static ST_boolean ismissing(double x){return x>=mock.missval;}
static ST_retcode message(char *s){(void)s;return 0;}
static const char *totals[]={"N","df","mobility","nested","singletons","iterations","group_saved"};
static const char *fields[]={"K","M","nested","exact"};
static double meta[7],factor[4][3];
static unsigned meta_seen,factor_seen[4];
static ST_retcode capture_macro(char *name,char *value)
{
    assert(!strncmp(name,"_chdfe_",7));name+=7;
    for(int j=0;j<7;j++)if(!strcmp(name,totals[j])) {
        assert(!(meta_seen&(1u<<j)));meta_seen|=1u<<j;meta[j]=atof(value);return 0;
    }
    for(int j=0;j<4;j++)for(int g=0;g<3;g++) {
        char key[32];snprintf(key,sizeof(key),"%s%d",fields[j],g+1);
        if(!strcmp(name,key)) { assert(!(factor_seen[j]&(1u<<g)));
            factor_seen[j]|=1u<<g;factor[j][g]=atof(value);return 0; }
    }
    abort();
}
static void case_run(int mode,int weighted,int cluster,int pool,double scale)
{
    int singleton=mode==1&&scale==1;
    assert(!live);clear_fixture();resize_fixture(96+singleton);
    int g=mode==0||mode==4?1:mode==2?3:2;
    int slope[3]={0};
    if(mode==0||mode==4)slope[0]=1;
    if(mode>=1&&mode<=3)slope[1]=2;
    if(mode==2)slope[2]=2;
    int ns=0;for(int j=0;j<g;j++)ns+=slope[j]!=0;
    int inputs=2+g+ns+cluster+weighted;nv=inputs+4;
    double expected[96],design[3][96],weights[96];
    for(int i=0;i<96;i++) {
        int pair=i/2,a=pair%3,b=(pair/3)%2,component=pair/24;
        double t=(pair/6)%4-1,u=component,fit=0;
        weights[i]=weighted?(i%2?2:1):1;
        expected[i]=(i%2?-1./weights[i]:1);
        if(mode==0)fit=(2+a)*t;
        if(mode>=1&&mode<=3)fit=3*a+(mode==3?0:(2+a)*t)+(mode==2?(4-a)*u:0);
        if(mode==5)fit=3*a+2*(2*a+b);
        if(mode==6)fit=3*(a+3*component)-2*(b+2*component);
        numbers[0][i]=(fit+expected[i])*scale;
        numbers[1][i]=(2*fit-3*expected[i])*scale;
        for(int j=0;j<g;j++)numbers[2+j][i]=a;
        if(mode==5)numbers[3][i]=2*a+b;
        if(mode==6) { numbers[2][i]=a+3*component;numbers[3][i]=b+2*component; }
        int col=2+g;
        for(int j=0;j<g;j++) {
            double s=mode==4?0:mode==3?1+a:j==2?u:t;
            design[j][i]=slope[j]?s:1;
            if(slope[j])numbers[col++][i]=s*scale;
        }
        if(cluster)numbers[col++][i]=numbers[2][i];
        if(weighted)numbers[col++][i]=weights[i]*scale;
        assert(col==inputs);numbers[nv-1][i]=1;
    }
    if(singleton) {
        for(int col=0;col<nv;col++)numbers[col][96]=col<inputs?999:mock.missval;
        numbers[nv-1][96]=1;
    }
    char args[256];int n=snprintf(args,sizeof(args),"2 %d %d %d 0 10000 13 %d 0 1 0",g,cluster,weighted,pool);
    for(int j=0;j<g;j++)n+=snprintf(args+n,sizeof(args)-n," %d",slope[j]);
    snprintf(args+n,sizeof(args)-n," tol=1e-12");
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(threads);
#endif
    for(int repeat=0;repeat<2;repeat++) {
        meta_seen=0;memset(factor_seen,0,sizeof(factor_seen));
        assert(!chdfe_main(args)&&!live);
        assert(meta_seen==127&&meta[0]==96&&meta[4]==singleton&&meta[5]>=0&&meta[6]==0);
        if(singleton)assert(numbers[inputs+2][96]==mock.missval&&numbers[inputs][96]==mock.missval);
        int df=0,mobility=0,nested=0;
        for(int j=0;j<g;j++) {
            int k=mode==5&&j==1?6:mode==6?(j==0?6:4):3;
            int nest=cluster&&!slope[j]&&(mode!=6||j==0);
            /* In mode 5 each finer B level is also nested in cluster A. */
            int red=nest?k:mode==3&&j==1?3:mode==4?3:
                mode==5&&j==1?3:mode==6&&j==1?2:0;
            if(mode==6&&cluster&&j==1)red=1;
            if(slope[j]&&scale<1e-100)red=k;
            int exact=nest||!slope[j];
            if(factor[0][j]!=k||factor[1][j]!=red||factor[2][j]!=nest||factor[3][j]!=exact)
                fprintf(stderr,"mode=%d weighted=%d cluster=%d scale=%g factor=%d K/M/nest/exact=%g/%g/%g/%g expected=%d/%d/%d/%d\n",mode,weighted,cluster,scale,j,factor[0][j],factor[1][j],factor[2][j],factor[3][j],k,red,nest,exact);
            assert(factor[0][j]==k&&factor[1][j]==red);
            assert(factor[2][j]==nest&&factor[3][j]==exact);
            df+=k-red;mobility+=red;if(nest)nested+=red;
        }
        for(int j=0;j<4;j++)assert(factor_seen[j]==(1u<<g)-1);
        assert(meta[1]==df&&meta[2]==mobility&&meta[3]==nested);
        for(int v=0;v<2;v++)for(int i=0;i<96;i++) {
            double residual=numbers[inputs+v][i]/scale;
            assert(fabs(residual-(v?-3:1)*expected[i])<1e-9);
            assert(numbers[inputs+2][i]==1);
        }
        /* Explicit weighted orthogonality for every category/slope column. */
        for(int j=0;j<g;j++)for(int level=0;level<6;level++) {
            double inner=0;
            for(int i=0;i<96;i++)if(numbers[2+j][i]==level)
                inner+=weights[i]*design[j][i]*numbers[inputs][i]/scale;
            assert(fabs(inner)<1e-8);
        }
        /* Re-project the already transformed columns: idempotence. */
        for(int v=0;v<2;v++)memcpy(numbers[v],numbers[inputs+v],96*sizeof(double));
    }
    clear_fixture();
}
int main(void)
{
    mock.nobs=nobs;mock.nvars=nvars;mock.nvar=nvars;mock.isstr=isstr;
    mock.ismissing=ismissing;mock.safestore=putnum;mock.macresave=capture_macro;
    mock.spouterr=mock.spoutsml=message;mock.missval=8.98846567431158e307;
    double scales[]={1,1e-200,1e200};
    for(threads=1;threads<=4;threads*=4)for(int mode=0;mode<=6;mode++)
        for(int weighted=0;weighted<2;weighted++)for(int cluster=0;cluster<2;cluster++)
            for(int pool=1;pool<=2;pool++)for(int s=0;s<3;s++)
                case_run(mode,weighted,cluster,pool,scales[s]);
    puts("PASS HDFE slope residuals, orthogonality, idempotence, metadata, nesting and scaling");return 0;
}
