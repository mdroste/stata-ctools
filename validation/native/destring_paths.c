/* Production conversion pipeline, with independent expected values and faults. */
#include "newcommands_mock.h"
#include "ctools_runtime.h"
static int scalar_calls,fail_scalar;
static double scalars[9];
static char failed_counts[256];
static const char *scalar_names[]={"_cdestring_n_converted","_cdestring_n_failed",
    "_cdestring_time_parse","_cdestring_time_load","_cdestring_time_convert",
    "_cdestring_time_store","_cdestring_time_total","_cdestring_threads_max",
    "_cdestring_openmp_enabled"};
static ST_retcode save_scalar(char *name,double value)
{
    if(++scalar_calls==fail_scalar)return 459;
    for(int i=0;i<9;i++)if(!strcmp(name,scalar_names[i])) { scalars[i]=value;return 0; }
    abort();
}
static ST_retcode save_macro(char *name,char *value)
{
    assert(!strcmp(name,"_cdestring_failed"));if(fail_macro)return 459;
    assert(strlen(value)<sizeof(failed_counts));strcpy(failed_counts,value);return 0;
}
static ST_retcode error_message(char *message){(void)message;return 0;}
stata_retcode ctools_store_filtered_rowpar(double *x,size_t n,int col,perm_idx_t *map)
{ return ctools_store_filtered(x,n,col,map); }
#include "ctools_types.c"
#include "cdestring/cdestring_impl.c"

static void init(int rows,int vars)
{
    assert(!live);clear_fixture();resize_fixture(rows);nv=2*vars;
    for(int v=0;v<vars;v++)types[v]=1;
    for(int v=0;v<vars;v++)for(int r=0;r<rows;r++)numbers[vars+v][r]=-12345;
    scalar_calls=fail_scalar=0;failed_counts[0]=0;
    for(int i=0;i<9;i++)scalars[i]=-999;
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(threads);
#endif
}
static int run(const char *args)
{ int rc=cdestring_main(args);assert(!live);return rc; }
static double missing(int letter)
{
    uint64_t bits;double value=mock.missval;memcpy(&bits,&value,8);
    bits+=(uint64_t)letter<<40;memcpy(&value,&bits,8);return value;
}
static void value_case(const char *text,const char *option,int valid,double expected)
{
    init(1,1);text_cell(0,0,text);char args[4096];
    snprintf(args,sizeof(args),"nvars=1 seqpairs force %s",option);
    assert(!run(args));
    if(scalars[0]!=valid||numbers[1][0]!=(valid?expected:mock.missval))
        fprintf(stderr,"conversion: [%s] [%s], valid=%g value=%.17g expected=%.17g\n",text,option,scalars[0],numbers[1][0],expected);
    assert(scalars[0]==valid&&scalars[1]==!valid);
    assert(numbers[1][0]==(valid?expected:mock.missval));
}
static void values(void)
{
    value_case("  +12.5\t","",1,12.5);
    value_case("-1.25d+2","",1,-125);
    value_case(".5D-2","",1,.005);
    value_case("+1.8000000000000X+001","",1,3);
    value_case("1.0000000000001X-3ff","",1,0x0.0000000000001p-1022);
    value_case("-1.0X+3ff","",1,-0x1p1023);
    const char *bad[]={"1d","1e+","NaN","NA","Inf","1e999","1 2",
        "1.0X+400","1.0X-400","1.12345678901234X+0","0.1X-3ff","1,2"};
    for(size_t i=0;i<sizeof(bad)/sizeof(*bad);i++)value_case(bad[i],"",0,0);
    value_case("","",1,missing(0));value_case(".","",1,missing(0));
    for(int i=1;i<=26;i++) { char s[]={'.',(char)('a'+i-1),0};value_case(s,"",1,missing(i)); }
    value_case("1.8e38","float",0,0);value_case("1e38","float",1,1e38);
    value_case("-1.8e38","float",0,0);
    value_case("1,25","dpcomma",1,1.25);value_case("1.25","dpcomma",0,0);
    value_case("1,2,3","dpcomma",0,0);value_case(",","dpcomma",0,0);
    value_case("12%","percent",1,.12);
    value_case("\xc2\xa0" "12" "\xe3\x80\x80","",1,12);
    value_case("é12é","ignore=é",1,12);
    value_case("ê12","ignore=é",0,0);
    value_case("\xc3" "12","ignore=é asbytes",1,12);
    value_case("\xc3" "12","ignore=é",0,0);
    value_case("\xff" "12" "\xed\xa0\x80","illegal",1,12);
    value_case("\xe3\x80\x80" "12" "\xc2\xa0","illegal",1,12);
    value_case("$ 1\t2\\","ignore=$\\s\\t\\\\",1,12);
    char options[2048];strcpy(options,"ignore=");memset(options+7,'x',1000);
    strcpy(options+1007,"é$");value_case("é$12",options,1,12);
    char long_value[2046];memset(long_value,'$',2043);strcpy(long_value+2043,"12");
    value_case(long_value,"ignore=$",1,12);
    /* Hexadecimal missing bit patterns are retained under force while the
     * source is counted as nonnumeric, matching real()/destring. */
    init(2,1);text_cell(0,0,"1.0010000000000X+3ff");text_cell(0,1,"1.fffffffffffffX+3ff");
    assert(!run("nvars=1 seqpairs force"));
    assert(scalars[0]==0&&scalars[1]==2);
    assert(numbers[1][0]==missing(1)&&numbers[1][1]==0x1.fffffffffffffp1023);
}
static int last_allocations,last_stores;
static int pipeline(int rows,int vars,int force,int percent,int filter,int af,int sf,int mf,int vf)
{
    init(rows,vars);int fail[4]={0},converted=0,total_failed=0;
    for(int v=0;v<vars;v++)for(int r=0;r<rows;r++) {
        selected[r]=!filter||r%3!=1;
        const char *s=v==0&&r%7==0?"bad":v==1&&r%5==0?"50%":v==2&&r%11==0?".a":"12.5";
        text_cell(v,r,s);
        if(selected[r]) {
            if((v==0&&r%7==0)||(!percent&&v==1&&r%5==0))fail[v]++;
            else converted++;
        }
    }
    for(int v=0;v<vars;v++)total_failed+=fail[v];
    fail_alloc=af;fail_store=sf;fail_macro=mf;fail_scalar=vf;
    char args[128];snprintf(args,sizeof(args),"nvars=%d seqpairs %s %s %s",vars,
        force?"force":"",percent?"percent":"",filter?"":"skipif");
    int rc=run(args);last_allocations=allocations;last_stores=stores;
    if(!rc) {
        assert(scalars[0]==converted&&scalars[1]==total_failed);
        char counts[256];int pos=0;
        for(int v=0;v<vars;v++)pos+=snprintf(counts+pos,sizeof(counts)-pos,"%s%d",v?" ":"",fail[v]);
        assert(!strcmp(counts,failed_counts));
        for(int v=0;v<vars;v++)for(int r=0;r<rows;r++) {
            double expected=-12345;
            if(selected[r]&&(force||!fail[v])) {
                expected=v==0&&r%7==0?missing(0):v==1&&r%5==0?(percent?50:missing(0)):
                    v==2&&r%11==0?missing(1):12.5;
                if(percent&&v==1)expected=expected>=mock.missval?mock.missval:expected/100;
            }
            assert(numbers[vars+v][r]==expected);
        }
#ifdef _OPENMP
        if(threads>1&&vars>1&&force)assert(worker_store);
#endif
    }
    clear_fixture();return rc;
}
int main(void)
{
    mock.nobs=nobs;mock.nvars=nvars;mock.nvar=nvars;mock.isstr=isstr;mock.isstrl=isstrl;
    mock.store=putnum;mock.safestore=putnum;mock.scalsave=save_scalar;mock.macresave=save_macro;
    mock.spouterr=error_message;mock.missval=8.98846567431158e307;
    values();
    for(threads=1;threads<=4;threads*=4) {
        for(int vars=1;vars<=3;vars+=2)for(int force=0;force<2;force++)
            for(int percent=0;percent<2;percent++)for(int filter=0;filter<2;filter++)
                assert(!pipeline(12001,vars,force,percent,filter,0,0,0,0));
        assert(!pipeline(12001,3,1,1,1,0,0,0,0));int n=last_stores;
        int points[]={1,n/2,n};
        for(int i=0;i<3;i++)assert(pipeline(12001,3,1,1,1,0,points[i],0,0)==459);
        assert(pipeline(20,3,1,1,0,0,0,1,0)==459);
        for(int i=1;i<=9;i++)assert(pipeline(20,3,1,1,0,0,0,0,i)==459);
        assert(!pipeline(4,3,1,1,0,0,0,0,0));n=last_allocations;
        for(int i=1;i<=n;i++)assert(pipeline(4,3,1,1,0,i,0,0,0)==920);
        for(int fault=0;fault<=8;fault++) {
            init(3,2);memset(selected,0,3);fail_macro=fault==8;fail_scalar=fault;
            assert(run("nvars=2 seqpairs verbose")== (fault?459:0));
            if(!fault) { assert(!strcmp(failed_counts,"0 0"));for(int i=0;i<7;i++)assert(scalars[i]==0); }
        }
        assert(!pipeline(20,3,1,1,0,0,0,0,0));
    }
    clear_fixture();puts("PASS destring values, filters, parallelism, allocation and publication faults");return 0;
}
