/* Production conversion pipeline, with independent expected values and faults. */
#include "newcommands_mock.h"
#include "ctools_runtime.h"
static int scalar_calls,fail_scalar;
static double scalars[10];
static char failed_counts[256],type_words[256],missing_words[256],seen_words[256];
static const char *raw_ignore = "";
static int ignore_read;
static const char *scalar_names[]={"_cdestring_n_converted","_cdestring_n_failed",
    "_cdestring_time_parse","_cdestring_time_load","_cdestring_time_convert",
    "_cdestring_time_store","_cdestring_time_total","_cdestring_threads_max",
    "_cdestring_openmp_enabled","_cdestring_nobs"};
static ST_retcode save_scalar(char *name,double value)
{
    if(++scalar_calls==fail_scalar)return 459;
    for(int i=0;i<10;i++)if(!strcmp(name,scalar_names[i])) { scalars[i]=value;return 0; }
    abort();
}
static ST_retcode save_macro(char *name,char *value)
{
    if(!strcmp(name,"_cdestring_ignore_read")) {
        assert(!strcmp(value,"1"));ignore_read=1;return fail_macro?459:0;
    }
    if(fail_macro)return 459;
    char *into=!strcmp(name,"_cdestring_failed")?failed_counts:
        !strcmp(name,"_cdestring_types")?type_words:
        !strcmp(name,"_cdestring_missing")?missing_words:
        !strcmp(name,"_cdestring_seen")?seen_words:NULL;
    assert(into&&strlen(value)<sizeof(failed_counts));strcpy(into,value);return 0;
}
static ST_retcode read_macro(char *name,char *value,ST_int size)
{
    assert(!strcmp(name,"___cds_ignore"));
    assert(strlen(raw_ignore)<(size_t)size);strcpy(value,raw_ignore);return 0;
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
    scalar_calls=fail_scalar=ignore_read=0;
    failed_counts[0]=type_words[0]=missing_words[0]=seen_words[0]=0;
    for(int i=0;i<10;i++)scalars[i]=-999;
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
static void release_regressions(void)
{
    int counts[]={INT_MAX,INT_MAX,123};
    assert(cdestring_total(counts,3)==UINT64_C(4294967417));
    raw_ignore="\" \\t";
    init(1,1);text_cell(0,0,"\"12\" \\t");
    assert(!run("nvars=1 seqpairs force ignorelen=4"));
    assert(ignore_read&&numbers[1][0]==12&&scalars[0]==1&&scalars[1]==0);
    raw_ignore="";
    /* A valid decimal longer than the parser's stack buffer must never
     * become missing under force when any allocation fails. */
    for(int exponent=0;exponent<2;exponent++) {
        char value[160];strcpy(value,"1.");memset(value+2,'1',140);
        strcpy(value+142,exponent?"d+0":"");
        int count=0;
        for(int fault=0;fault<=count;fault++) {
            init(1,1);text_cell(0,0,value);fail_alloc=fault;
            int rc=run("nvars=1 seqpairs force skipif");
            if(!fault) {
                assert(!rc&&numbers[1][0]>1.11&&numbers[1][0]<1.12);
                count=allocations;
            } else {
                assert(rc==920&&stores==0&&numbers[1][0]==-12345);
            }
        }
    }
    clear_fixture();
}
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
/* parse keeps the numbers and describes them; the destinations exist only
 * for store, one per variable parse reported a type for. */
static void staged_rows(void)
{
    const char *a[]={"1","2","100",".a","","-127"};
    const char *b[]={"1.5","x","3","4","5","6"};
    const char *c[]={"$1,000","2","3","40000","5","6"};
    for(int r=0;r<6;r++) { text_cell(0,r,a[r]);text_cell(1,r,b[r]);text_cell(2,r,c[r]); }
}
static void staged(void)
{
    const double a[]={1,2,100,missing(1),missing(0),-127};
    const double c[]={1000,2,3,40000,5,6};
    raw_ignore="$,";
    /* force: every variable is stored, the nonnumeric value as missing */
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif force ignorelen=2")&&live>0&&!stores);
    assert(!strcmp(failed_counts,"0 1 0")&&!strcmp(type_words,"byte double long"));
    assert(!strcmp(missing_words,"2 1 0")&&!strcmp(seen_words,". . 0,1"));
    assert(scalars[9]==6&&scalars[0]==17&&scalars[1]==1);
    assert(!cdestring_main("store")&&!live);
    for(int r=0;r<6;r++) {
        assert(numbers[0][r]==a[r]&&numbers[2][r]==c[r]);
        assert(numbers[1][r]==(r==0?1.5:r==1?mock.missval:r+1));
    }
    /* without force the variable with a nonnumeric value has no destination */
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif ignorelen=2")&&live>0);
    assert(!strcmp(type_words,"byte . long")&&!strcmp(missing_words,"2 . 0"));
    nv=3;assert(cdestring_main("store")==198&&!live&&!stores);
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif ignorelen=2"));
    nv=2;assert(!cdestring_main("store")&&!live);
    for(int r=0;r<6;r++)assert(numbers[0][r]==a[r]&&numbers[1][r]==c[r]);
    /* discard, a store with nothing parsed, and another call free the numbers */
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif force ignorelen=2")&&live>0);
    assert(!cdestring_main("discard")&&!live&&!stores);
    assert(cdestring_main("store")==198&&!live);
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif force ignorelen=2")&&live>0);
    nv=6;assert(!run("nvars=3 seqpairs force ignorelen=2"));
    /* the data must not change between the calls */
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif force ignorelen=2"));
    nr=5;assert(cdestring_main("store")==459&&!live&&!stores);nr=6;
    raw_ignore="";
    /* float: values are rounded to float before the type is chosen */
    init(3,3);nv=3;
    text_cell(0,0,"0.1");text_cell(0,1,"16777217");text_cell(0,2,"3");
    text_cell(1,0,"1");text_cell(1,1,"2");text_cell(1,2,"40000");
    text_cell(2,0,"1");text_cell(2,1,"32740");text_cell(2,2,"-32767");
    assert(!cdestring_main("parse nvars=3 skipif float"));
    assert(!strcmp(type_words,"float float int"));
    assert(!cdestring_main("store")&&!live);
    assert(numbers[0][0]==(double)0.1f&&numbers[0][1]==16777216&&numbers[2][1]==32740);
    /* no observation selected: byte destinations that stay missing */
    init(3,1);nv=1;memset(selected,0,3);
    assert(!cdestring_main("parse nvars=1")&&scalars[9]==0);
    assert(!strcmp(type_words,"byte")&&!strcmp(missing_words,"0")&&!strcmp(seen_words,"."));
    assert(!cdestring_main("store")&&!live&&!stores);
    /* an allocation failure in parse keeps nothing */
    raw_ignore="$,";
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif force ignorelen=2"));
    int count=allocations;
    assert(!cdestring_main("discard")&&!live);
    for(int fault=1;fault<=count;fault++) {
        init(6,3);staged_rows();nv=3;fail_alloc=fault;
        assert(cdestring_main("parse nvars=3 skipif force ignorelen=2")==920&&!live&&!stores);
    }
    /* a failed write frees the numbers */
    init(6,3);staged_rows();nv=3;
    assert(!cdestring_main("parse nvars=3 skipif force ignorelen=2"));
    fail_store=1;assert(cdestring_main("store")==459&&!live);
    raw_ignore="";clear_fixture();
}
int main(void)
{
    mock.nobs=nobs;mock.nvars=nvars;mock.nvar=nvars;mock.isstr=isstr;mock.isstrl=isstrl;
    mock.store=putnum;mock.safestore=putnum;mock.scalsave=save_scalar;mock.macresave=save_macro;
    mock.macuse=read_macro;
    mock.spouterr=error_message;mock.missval=8.98846567431158e307;
    values();
    for(threads=1;threads<=4;threads*=4) {
        release_regressions();
        staged();
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
    clear_fixture();puts("PASS destring values, filters, staged parse/store, parallelism, allocation and publication faults");return 0;
}
