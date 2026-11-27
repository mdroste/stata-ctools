/* Production encode/decode protocols, mapping oracles and injected failures. */
#include "newcommands_mock.h"
/* File paths are converted with the C library's allocator (Windows), so
 * allocation failures stay those of the code under test on every system. */
#undef malloc
#undef calloc
#undef free
#include "ctools_file.h"
#define malloc tracked_malloc
#define calloc tracked_calloc
#define free tracked_free
static atomic_int reads;
static int fail_read,fail_scalar,scalar_calls,fail_labels,label_count;
static double decoded,unique,fallback;
static const char *wire;
static int label_codes[16];static char label_text[16][64];
/* Protocol 3: the in range, a strL source, mapping and label-dataset files. */
static int range_lo=1,range_hi,strl_src;
static double protocol,unsupported,route;
static char map_path[1024],dta_path[1024];
static ST_int in1(void){return range_lo;}
static ST_int in2(void){return range_hi?range_hi:nr;}
static ST_boolean isstrl_v3(ST_int v){return strl_src&&v==1;}
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
/* strLs are read on the calling thread only; "<binary>" holds a NUL byte. */
static ST_int get_strl(ST_int v,ST_int r,char *s,ST_int len)
{
#ifdef _OPENMP
    assert(omp_get_thread_num()==0);
#endif
    if(atomic_fetch_add(&reads,1)+1==fail_read)return -1;
    assert(v==1&&r>=1&&r<=nr&&len>0);
    const char *t=string_cell(0,r-1);
    if(!strcmp(t,"<binary>")){t="a\0b";int n=len<3?len:3;memcpy(s,t,(size_t)n);return n;}
    int n=(int)strlen(t);if(n>len)n=len;memcpy(s,t,(size_t)n);return n;
}
static ST_retcode get_macro(char *key,char *out,ST_int size)
{
    if(fail_macro)return 459;
    const char *s=!strcmp(key,"___ce_labelfile")?"labels":
        !strcmp(key,"___ce_mapfile")||!strcmp(key,"___cd_mapfile")?map_path:
        !strcmp(key,"___ce_labeldta")?dta_path:!strcmp(key,"___ce_dtavar")?"dvar":
        !strcmp(key,"___ce_dtalab")?"dlab":wire;
    assert(s&&strlen(s)<(size_t)size);strcpy(out,s);return 0;
}
static ST_retcode put_scalar(char *key,double x)
{
    if(++scalar_calls==fail_scalar)return 459;
    if(!strcmp(key,"_cdecode_n_decoded"))decoded=x;
    if(!strcmp(key,"_cencode_n_unique"))unique=x;
    if(!strcmp(key,"_cencode_fallback"))fallback=x;
    if(!strcmp(key,"_cencode_protocol"))protocol=x;
    if(!strcmp(key,"_cencode_unsupported"))unsupported=x;
    if(!strcmp(key,"_cencode_label_route"))route=x;
    return 0;
}
/* The marker is read within the in range, as the loader does. */
stata_retcode ctools_data_select(ctools_filtered_data *fd,size_t start,size_t end,int flags)
{
    (void)start;(void)end;(void)flags;ctools_filtered_data_init(fd);
    fd->obs_map=malloc((size_t)nr*sizeof(perm_idx_t)+1);if(!fd->obs_map)return STATA_ERR_MEMORY;
    for(int i=in1()-1;i<in2();i++)if(selected[i])fd->obs_map[fd->data.nobs++]=i+1;
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
/* The .dta schema module grows its buffers with realloc, which allocation
 * tracking does not cover; it keeps the C library's allocator. */
#undef malloc
#undef calloc
#undef free
#include "formats/ctools_dta_schema.c"
#define malloc tracked_malloc
#define calloc tracked_calloc
#define free tracked_free
#include "cencode/cencode_impl.c"
#include "cdecode/cdecode_impl.c"
static void init(int rows)
{
    assert(!live);clear_fixture();resize_fixture(rows);nv=2;
    reads=0;fail_read=fail_scalar=scalar_calls=fail_labels=label_count=0;
    decoded=unique=fallback=protocol=unsupported=route=-999;
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
/* ---- Protocol 3 ---- */
#ifdef _WIN32
#include <process.h>
#define lp_getpid _getpid
#else
#include <unistd.h>
#define lp_getpid getpid
#endif
/* A file of this test in the system's temporary directory. */
static void temp_path(char *out,size_t cap,const char *name)
{
    const char *dir=getenv("TMPDIR");
    if(!dir||!*dir)dir=getenv("TEMP");
    if(!dir||!*dir)dir=getenv("TMP");
#ifdef _WIN32
    if(!dir||!*dir)dir=".";
#else
    if(!dir||!*dir)dir="/tmp";
#endif
    snprintf(out,cap,"%s/labels_paths_%d_%s",dir,(int)lp_getpid(),name);
}
static void write_text(const char *path,const char *text)
{
    FILE *f=fopen(path,"wb");assert(f);
    assert(fwrite(text,1,strlen(text),f)==strlen(text));assert(!fclose(f));
}
enum { DTA_NEW=3000 };
/* -1 alpha, -2 "", else k for s<k> */
static int dta_kind(int i)
{ return i>=10&&i<10+DTA_NEW?i-10:i%5==0?-1:i%5==1?-2:i%DTA_NEW; }
/* Rows 10..3009 hold s0000..s2999 once each, other rows alpha, "" or an s
 * value: more new mappings than st_vlmodify() installs quickly, so the
 * label comes back as a dataset (dtacost=0), read back here. The in range
 * leaves out rows at both ends, which fillrows sets to missing. */
static int encode_dta_case(int existing,int from_file,int af,int sf,int vf)
{
    const int rows=12001;init(rows);types[0]=1;range_lo=3;range_hi=rows-2;
    wire="4 -7 5 alpha 3 5 alpha 9 2 é 1000 1 b |";
    if(from_file)write_text(map_path,wire);
    char text[32];
    for(int i=0;i<rows;i++) {
        int kind=dta_kind(i);
        if(kind==-1)text_cell(0,i,"alpha");
        else if(kind==-2)text_cell(0,i,"");
        else {snprintf(text,sizeof(text),"s%04d",kind);text_cell(0,i,text);}
        numbers[1][i]=-99;
    }
    remove(dta_path);
    fail_alloc=af;fail_store=sf;fail_scalar=vf;
    char args[600];
    snprintf(args,sizeof(args),"v3 1 2 labelfile=6 labeldta=%zu dtacost=0 stream fillrows labelmap=%zu %s",
        strlen(dta_path),existing?strlen(wire):0,from_file?"mapfile":"");
    int rc=cencode_main(args);
    range_lo=1;range_hi=0;
    if(!rc) {
        assert(protocol==3&&unsupported==0&&route==2&&label_count==0);
        for(int i=0;i<rows;i++) {
            double want=mock.missval;
            int kind=dta_kind(i);
            if(i>=2&&i<rows-2&&kind!=-2)want=kind==-1?(existing?-7:1):(existing?1001:2)+kind;
            assert(numbers[1][i]==want);
        }
        ctools_dta_schema s;char msg[256];
        assert(!ctools_dta_schema_read(dta_path,&s,msg,sizeof(msg)));
        assert(s.nvars==1&&!strcmp(s.vars[0].name,"dvar")&&s.vars[0].type==CTOOLS_DTA_CODE_BYTE);
        long at=ctools_dta_schema_label_find(&s,"dlab");assert(at>=0&&s.nlabels==1);
        size_t count;assert(!ctools_dta_schema_label_count(&s,(size_t)at,&count));
        assert(count==(size_t)DTA_NEW+(existing?4:1));
        const int old_values[]={-7,3,9,1000};const char *old_texts[]={"alpha","alpha","é","b"};
        for(size_t k=0;k<count;k++) {
            int32_t v;const char *t;assert(!ctools_dta_schema_label_entry(&s,(size_t)at,k,&v,&t));
            size_t fresh=existing?k-4:k-1;
            if(existing&&k<4)assert(v==old_values[k]&&!strcmp(t,old_texts[k]));
            else if(!existing&&k==0)assert(v==1&&!strcmp(t,"alpha"));
            else {snprintf(text,sizeof(text),"s%04zu",fresh);
                assert(v==(existing?1001:2)+(int32_t)fresh&&!strcmp(t,text));}
        }
        ctools_dta_schema_free(&s);
    }
    done();return rc;
}
/* A strL source: read on the calling thread; bad=1 puts a value over 2,045
 * bytes in, bad=2 a binary one, which encode must take. */
static int encode_strl_case(int rows,int bad,int af,int rf,int sf)
{
    init(rows);types[0]=1;strl_src=1;
    for(int i=0;i<rows;i++) {text_cell(0,i,i%3?"b":i%5?"a":"");numbers[1][i]=-99;}
    if(bad==1) {static char big[2101];memset(big,'x',2100);big[2100]=0;text_cell(0,rows/2,big);}
    if(bad==2)text_cell(0,rows/2,"<binary>");
    fail_alloc=af;fail_read=rf;fail_store=sf;
    int rc=cencode_main("v3 1 2 labelfile=6 stream fillrows");
    strl_src=0;
    if(!rc) {
        assert(protocol==3&&unsupported==(bad?1:0));
        if(bad)assert(!stores&&numbers[1][0]==-99);
        else {
            assert(route==1&&label_count==2&&!strcmp(label_text[0],"a")&&!strcmp(label_text[1],"b"));
            for(int i=0;i<rows;i++)assert(numbers[1][i]==(i%3?2:i%5?1:mock.missval));
        }
    }
    done();return rc;
}
/* cdecode v3: the mapping from a file and fillrows with if and in. */
static int decode_v3_case(int rows,int filter,int af,int sf)
{
    init(rows);types[1]=1;nv=filter==2?3:2;range_lo=2;range_hi=rows>4?rows-3:rows;
    wire="5 -7 5 alpha 2 1 b 3 5 alpha .a 2 é .z 0  |";
    write_text(map_path,wire);
    double values[]={-7,2,3,extmissing(1),extmissing(26),0,.5};
    const char *texts[]={"alpha","b","alpha","é","","",""};
    for(int i=0;i<rows;i++) { numbers[0][i]=values[i%7];selected[i]=!filter||i%3!=1;
        numbers[2][i]=selected[i];text_cell(1,i,"zz"); }
    fail_alloc=af;fail_store=sf;
    char args[200];snprintf(args,sizeof(args),"v3 maplen=%zu width=5 mapfile fillrows %s",strlen(wire),
        filter==2?"ifcond touse":filter?"ifcond":"");
    int rc=cdecode_main(args);
    int lo=range_lo,hi=range_hi;range_lo=1;range_hi=0;
    if(!rc) for(int i=0;i<rows;i++) {
        int in=i+1>=lo&&i+1<=hi&&selected[i];
        assert(!strcmp(string_cell(1,i),in?texts[i%7]:""));
    }
    done();return rc;
}
static void protocol3(void)
{
    /* The question alone */
    init(1);types[0]=1;assert(!cencode_main("v3")&&protocol==3&&scalar_calls==1&&!stores);done();
    init(1);types[0]=1;fail_scalar=1;assert(cencode_main(" v3 ")==459);done();
    temp_path(map_path,sizeof(map_path),"map.tmp");
    temp_path(dta_path,sizeof(dta_path),"label.dta");
    for(int existing=0;existing<2;existing++)for(int from_file=0;from_file<=existing;from_file++) {
        assert(!encode_dta_case(existing,from_file,0,0,0));int n=last_alloc,ns=last_store,sc=last_scalar;
        for(int f=1;f<=n;f++)assert(encode_dta_case(existing,from_file,f,0,0)==920);
        int sp[]={1,ns/2,ns};for(int j=0;j<3;j++)assert(encode_dta_case(existing,from_file,0,sp[j],0)==459);
        for(int f=1;f<=sc;f++)assert(encode_dta_case(existing,from_file,0,0,f)==459);
    }
    /* A label dataset that cannot be written */
    temp_path(dta_path,sizeof(dta_path),"no_such_directory/label.dta");
    assert(encode_dta_case(0,0,0,0,0)==693);
    temp_path(dta_path,sizeof(dta_path),"label.dta");
    for(int bad=0;bad<3;bad++) {
        /* Reads stop at a value encode must take, so every injected one fails. */
        assert(!encode_strl_case(12001,bad,0,0,0));int nr_=last_read,ns=last_store;
        int rp[]={1,nr_/2,nr_};for(int j=0;j<3;j++)assert(encode_strl_case(12001,bad,0,rp[j],0)==459);
        if(!bad){int sp[]={1,ns/2,ns};for(int j=0;j<3;j++)assert(encode_strl_case(12001,bad,0,0,sp[j])==459);}
        assert(!encode_strl_case(10,bad,0,0,0));int n=last_alloc;
        if(!bad)for(int f=1;f<=n;f++)assert(encode_strl_case(10,bad,f,0,0)==920);
    }
    for(int filter=0;filter<3;filter++) {
        assert(!decode_v3_case(12001,filter,0,0));int ns=last_store;
        int sp[]={1,ns/2,ns};for(int j=0;j<3;j++)assert(decode_v3_case(12001,filter,0,sp[j])==459);
        assert(!decode_v3_case(10,filter,0,0));int n=last_alloc;
        for(int f=1;f<=n;f++)assert(decode_v3_case(10,filter,f,0)==920);
    }
    /* A v2 call cannot ask for a file: it reads the macro */
    init(1);types[1]=1;wire="1 1 1 a |";write_text(map_path,"garbage!!");
    assert(cdecode_main("v2 maplen=9 width=5 mapfile fillall")==0&&decoded==0);done();
    /* A missing or shorter mapping file */
    init(1);types[1]=1;remove(map_path);assert(cdecode_main("v3 maplen=9 width=5 mapfile")==198);done();
    init(1);types[0]=1;text_cell(0,0,"a");write_text(map_path,"1 1 1 a ");
    assert(cencode_main("v3 1 2 stream mapfile labelmap=9")==198);done();
    remove(map_path);remove(dta_path);
}
int main(void)
{
    mock.nobs=nobs;mock.nobs1=in1;mock.nobs2=in2;mock.nvars=nvars;mock.nvar=nvars;
    mock.isstr=isstr;mock.isstrl=isstrl_v3;mock.strldata=get_strl;mock.ismissing=ismissing;mock.selobs=ifobs;
    mock.vdata=get_number;mock.sdata=get_string;mock.safestore=putnum;mock.store=putnum;mock.sstore=putstr;
    mock.macuse=get_macro;mock.scalsave=put_scalar;mock.spouterr=mock.spoutsml=message;
    mock.missval=8.98846567431158e307;
    malformed();
    for(threads=1;threads<=4;threads*=4)protocol3();
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
