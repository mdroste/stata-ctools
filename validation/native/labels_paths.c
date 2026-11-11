/* Production encode/decode protocols, mapping oracles and injected failures. */
#include "newcommands_mock.h"
static atomic_int reads;
static int fail_read,fail_scalar,scalar_calls,fail_labels,label_count;
static double decoded,unique,fallback;
static const char *wire;
static int label_codes[16];static char label_text[16][64];
static ST_int first(void){return 1;}
static ST_boolean ifobs(ST_int r){
#ifdef _OPENMP
    assert(omp_get_thread_num()==0);
#endif
    return selected[r-1];
}
static ST_retcode message(char *s){(void)s;return 0;}
static ST_boolean ismissing(double x){return x>=mock.missval;}
static ST_retcode get_number(ST_int v,ST_int r,double *x)
{
    if(atomic_fetch_add(&reads,1)+1==fail_read)return 459;
    assert(v>=1&&v<=nv&&r>=1&&r<=nr);*x=numbers[v-1][r-1];return 0;
}
static ST_retcode get_string(ST_int v,ST_int r,char *s)
{
    if(atomic_fetch_add(&reads,1)+1==fail_read)return 459;
    assert(v==1&&r>=1&&r<=nr);strcpy(s,string_cell(0,r-1));return 0;
}
static ST_retcode get_macro(char *key,char *out,ST_int size)
{
    if(fail_macro)return 459;
    const char *s=!strcmp(key,"___ce_labelfile")?"labels":wire;
    assert(s&&strlen(s)<(size_t)size);strcpy(out,s);return 0;
}
static ST_retcode put_scalar(char *key,double x)
{
    if(++scalar_calls==fail_scalar)return 459;
    if(!strcmp(key,"_cdecode_n_decoded"))decoded=x;
    if(!strcmp(key,"_cencode_n_unique"))unique=x;
    if(!strcmp(key,"_cencode_fallback"))fallback=x;
    return 0;
}
stata_retcode ctools_data_select(ctools_filtered_data *fd,size_t start,size_t end,int flags)
{
    (void)start;(void)end;(void)flags;ctools_filtered_data_init(fd);
    fd->obs_map=malloc((size_t)nr*sizeof(perm_idx_t));if(!fd->obs_map)return STATA_ERR_MEMORY;
    for(int i=0;i<nr;i++)if(selected[i])fd->obs_map[fd->data.nobs++]=i+1;
    return STATA_OK;
}
stata_retcode ctools_store_filtered_rowpar(double *x,size_t n,int col,perm_idx_t *map)
{ return ctools_store_filtered(x,n,col,map); }
static char *tracked_strdup(const char *s)
{ size_t n=strlen(s)+1;char *p=malloc(n);if(p)memcpy(p,s,n);return p; }
#define strdup tracked_strdup
#include "ctools_arena.c"
/* The binary label writer has its own file-framing tests. Here intercept its
 * boundary to check the encoder's independent mapping and writer-error path. */
#define ctools_label_write_file production_label_write_file
#include "ctools_hash.c"
#undef ctools_label_write_file
int ctools_label_write_file(const char **text,const int *codes,size_t n,const char *path)
{
    assert(!strcmp(path,"labels")&&n<=16);if(fail_labels)return -1;label_count=(int)n;
    for(size_t i=0;i<n;i++){label_codes[i]=codes[i];assert(strlen(text[i])<64);strcpy(label_text[i],text[i]);}
    return 0;
}
#include "ctools_parse.c"
#include "cencode/cencode_impl.c"
#include "cdecode/cdecode_impl.c"
static void init(int rows)
{
    assert(!live);clear_fixture();resize_fixture(rows);nv=2;
    reads=0;fail_read=fail_scalar=scalar_calls=fail_labels=label_count=0;
    decoded=unique=fallback=-999;
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(threads);
#endif
}
static double extmissing(int n)
{ uint64_t bits=CDECODE_SYSMISS_BITS+((uint64_t)n<<40);double d;memcpy(&d,&bits,8);return d; }
static int last_alloc,last_store,last_read,last_scalar;
static void done(void)
{
    assert(!live);last_alloc=allocations;last_store=stores;last_read=reads;last_scalar=scalar_calls;
    clear_fixture();
}
static int decode_case(int rows,int sparse,int filter,int af,int rf,int sf,int vf)
{
    init(rows);types[1]=1;nv=filter==2?3:2;
    wire=sparse?"5 -7 5 alpha 1000000000 1 b 3 5 alpha .a 2 é .z 0  |":
        "5 -7 5 alpha 2 1 b 3 5 alpha .a 2 é .z 0  |";
    double values[]={-7,sparse?1e9:2,3,extmissing(1),extmissing(26),0,.5,extmissing(0),NAN,INFINITY,-INFINITY};
    const char *texts[]={"alpha","b","alpha","é","","","","","","",""};
    int count=0;
    for(int i=0;i<rows;i++) { numbers[0][i]=values[i%11];selected[i]=filter==3?0:!filter||i%3!=1;
        numbers[2][i]=selected[i];text_cell(1,i,"");if(selected[i]&&*texts[i%11])count++; }
    fail_alloc=af;fail_read=rf;fail_store=sf;fail_scalar=vf;
    char args[200];snprintf(args,sizeof(args),"v2 maplen=%zu width=5 %s",strlen(wire),
        filter==2?"ifcond touse":filter?"ifcond":"fillall");
    int rc=cdecode_main(args);
    if(!rc) {
        assert(decoded==count&&scalar_calls==4);
        for(int i=0;i<rows;i++)assert(!strcmp(string_cell(1,i),selected[i]?texts[i%11]:""));
#ifdef _OPENMP
        if(threads>1&&rows>=2*MIN_OBS_PER_THREAD)assert(worker_store);
#endif
    }
    done();return rc;
}
static int encode_case(int rows,int stream,int existing,int filter,int af,int rf,int sf,int vf,int lf)
{
    init(rows);types[0]=1;nv=filter==2?3:2;
    wire="4 -7 5 alpha 3 5 alpha 9 2 é 1000 1 b |";
    const char *text[]={"b","alpha","é","","zz"};
    const double fresh[]={2,1,4,0,3},mapped[]={1000,-7,9,0,1001};
    for(int i=0;i<rows;i++) { text_cell(0,i,text[i%5]);numbers[1][i]=mock.missval;
        selected[i]=filter==3?0:!filter||i%3!=1;numbers[2][i]=selected[i]; }
    fail_alloc=af;fail_read=rf;fail_store=sf;fail_scalar=vf;fail_labels=lf;
    char args[256];snprintf(args,sizeof(args),"1 2 %s labelfile=6 labelmap=%zu %s",stream?"stream":"",
        existing?strlen(wire):0,filter==2?"ifcond touse":filter?"ifcond":"fillall");
    int rc=cencode_main(args);
    if(!rc) {
        for(int i=0;i<rows;i++) {
            double expected=selected[i]&&i%5!=3?(existing?mapped:fresh)[i%5]:mock.missval;
            assert(numbers[1][i]==expected);
        }
        if(filter==3)assert(unique==0&&label_count==0);
        if(rows>=5&&filter!=3) {
            if(unique!=4||label_count!=(existing?1:4))
                fprintf(stderr,"encode rows=%d stream=%d existing=%d filter=%d allocation=%d unique=%g labels=%d\n",rows,stream,existing,filter,af,unique,label_count);
            assert(unique==4&&label_count==(existing?1:4));
            if(existing)assert(label_codes[0]==1001&&!strcmp(label_text[0],"zz"));
            else { const char *sorted[]={"alpha","b","zz","é"};
                for(int i=0;i<4;i++)assert(label_codes[i]==i+1&&!strcmp(label_text[i],sorted[i])); }
        }
#ifdef _OPENMP
        if(stream&&threads>1&&rows>=2*MIN_OBS_PER_THREAD)assert(worker_store);
#endif
    }
    done();return rc;
}
static void malformed(void)
{
    const char *bad[]={"0 x","0 |x","1 1 1 a ","1 1 2 a |","1 1 -1 a |",
        "1 999999999999999999999999999999 1 a |","2 1 1 a 1 1 b |","1 1 99999999999 a |"};
    for(size_t i=0;i<sizeof(bad)/sizeof(*bad);i++)for(int encode=0;encode<2;encode++) {
        init(1);wire=bad[i];char args[128];
        if(encode) { types[0]=1;text_cell(0,0,"a");
            snprintf(args,sizeof(args),"1 2 stream labelmap=%zu",strlen(wire));assert(cencode_main(args)==198); }
        else { types[1]=1;snprintf(args,sizeof(args),"v2 maplen=%zu width=5",strlen(wire));assert(cdecode_main(args)==198); }
        assert(!stores);done();
    }
    init(1);types[1]=1;wire="2 .a 1 a .a 1 b |";char args[128];
    snprintf(args,sizeof(args),"v2 maplen=%zu width=5",strlen(wire));assert(cdecode_main(args)==198);done();
}
int main(void)
{
    mock.nobs=nobs;mock.nobs1=first;mock.nobs2=nobs;mock.nvars=nvars;mock.nvar=nvars;
    mock.isstr=isstr;mock.isstrl=isstrl;mock.ismissing=ismissing;mock.selobs=ifobs;
    mock.vdata=get_number;mock.sdata=get_string;mock.safestore=putnum;mock.store=putnum;mock.sstore=putstr;
    mock.macuse=get_macro;mock.scalsave=put_scalar;mock.spouterr=mock.spoutsml=message;
    mock.missval=8.98846567431158e307;
    malformed();
    for(threads=1;threads<=4;threads*=4) {
        for(int sparse=0;sparse<2;sparse++)for(int filter=0;filter<3;filter++) {
            assert(!decode_case(12001,sparse,filter,0,0,0,0));int ns=last_store,nr_=last_read;
            int sp[]={1,ns/2,ns},rp[]={1,nr_/2,nr_};
            for(int j=0;j<3;j++) { assert(decode_case(12001,sparse,filter,0,0,sp[j],0)==459);
                assert(decode_case(12001,sparse,filter,0,rp[j],0,0)==459); }
            assert(!decode_case(10,sparse,filter,0,0,0,0));int n=last_alloc;
            for(int f=1;f<=n;f++)assert(decode_case(10,sparse,filter,f,0,0,0)==920);
            for(int f=1;f<=4;f++)assert(decode_case(10,sparse,filter,0,0,0,f)==459);
        }
        for(int stream=0;stream<2;stream++)for(int existing=0;existing<=stream;existing++)for(int filter=0;filter<3;filter++) {
            assert(!encode_case(12001,stream,existing,filter,0,0,0,0,0));int ns=last_store,nr_=last_read;
            int sp[]={1,ns/2,ns},rp[]={1,nr_/2,nr_};
            for(int j=0;j<3;j++) { assert(encode_case(12001,stream,existing,filter,0,0,sp[j],0,0)==459);
                if(stream)assert(encode_case(12001,stream,existing,filter,0,rp[j],0,0,0)==459); }
            assert(!encode_case(10,stream,existing,filter,0,0,0,0,0));int n=last_alloc,sc=last_scalar;
            for(int f=1;f<=n;f++)assert(encode_case(10,stream,existing,filter,f,0,0,0,0)==920);
            for(int f=1;f<=sc;f++)assert(encode_case(10,stream,existing,filter,0,0,0,f,0)==459);
            assert(encode_case(10,stream,existing,filter,0,0,0,0,1)==920);
            assert(!encode_case(0,stream,existing,filter,0,0,0,0,0));
        }
        assert(!decode_case(0,0,1,0,0,0,0));
        assert(!decode_case(10,0,3,0,0,0,0));
        for(int stream=0;stream<2;stream++)assert(!encode_case(10,stream,0,3,0,0,0,0,0));
        assert(!decode_case(200001,1,2,0,0,0,0));
        int dn=last_store;
        assert(!encode_case(200001,1,1,2,0,0,0,0,0));int en=last_store;
        int dp[]={1,dn/2,dn},ep[]={1,en/2,en};
        for(int j=0;j<3;j++) {
            assert(decode_case(200001,1,2,0,0,dp[j],0)==459);
            assert(encode_case(200001,1,1,2,0,0,ep[j],0,0)==459);
        }
    }
    init(2);types[0]=1;text_cell(0,0,"a");text_cell(0,1,"new");numbers[1][0]=numbers[1][1]=-99;
    wire="1 -7 1 a |";char args[128];snprintf(args,sizeof(args),"1 2 stream noextend labelmap=%zu",strlen(wire));
    assert(!cencode_main(args)&&fallback==1&&!stores&&numbers[1][0]==-99&&numbers[1][1]==-99);done();
    for(int encode=0;encode<2;encode++) {
        init(1);wire="1 1 1 a |";fail_macro=1;
        if(encode){types[0]=1;text_cell(0,0,"a");assert(cencode_main("1 2 stream labelmap=9")==198);}
        else {types[1]=1;assert(cdecode_main("v2 maplen=9 width=5")==198);}
        assert(!stores);done();
    }
    puts("PASS encode/decode mapping, protocols, transfer/publication faults, allocation and parallel chunks");return 0;
}
