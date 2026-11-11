/* Focused regressions for confirmed November 6 audit defects. */
#ifdef NDEBUG
#error Native audit regressions require assertions
#endif
#if defined(CTOOLS_CASE_astra_encode)
#include <assert.h>
#include <string.h>
#include "ctools_types.h"
#include "ctools_parse.c"
#include "ctools_arena.c"
#include "ctools_hash.c"
#include "cencode/cencode_impl.c"
#include "cgroup/cgroup_impl.c"
static ST_plugin host;
ST_plugin *_stata_=&host;
static int rows=3, fail_scalar, scalar_calls;
static char *input[]={"a","b","c"};
static double output[3];
static perm_idx_t map[]={1,2,3}, order[]={0,1,2};
static stata_variable col;
static ST_int one(void){return 1;}
static ST_int last(void){return rows;}
static ST_int vars(void){return 2;}
static ST_boolean isstring(ST_int v){return v==1;}
static ST_boolean islong(ST_int v){(void)v;return 0;}
static ST_boolean ismissing(ST_double x){return x>=host.missval;}
static ST_retcode message(char *s){(void)s;return 0;}
static ST_retcode scalar(char *s, double x){(void)s;(void)x;return ++scalar_calls==fail_scalar?459:0;}
static ST_retcode readstring(ST_int v,ST_int i,char *s){assert(v==1);strcpy(s,input[i-1]);return 0;}
static ST_retcode writevalue(ST_int v,ST_int i,double x){assert(v==2&&i>=1&&i<=rows);output[i-1]=x;return 0;}
static ST_retcode getmacro(char *n,char *s,ST_int size){(void)n;assert(size>9);strcpy(s,"1 7 1 b |");return 0;}
int ctools_get_max_threads(void){return 1;}
int ctools_io_strl_readable(int v,int i){(void)v;(void)i;abort();}
double ctools_timer_seconds(void){return 0;}
void ctools_filtered_data_init(ctools_filtered_data *d){memset(d,0,sizeof(*d));}
void ctools_filtered_data_free(ctools_filtered_data *d){(void)d;}
stata_retcode ctools_data_select(ctools_filtered_data *d,size_t a,size_t b,int f){(void)d;(void)a;(void)b;(void)f;abort();}
stata_retcode ctools_data_load(ctools_filtered_data *d,int *v,size_t n,size_t a,size_t b,int f){
    (void)v;(void)n;(void)a;(void)b;(void)f;
    col.type=STATA_TYPE_STRING;col.nobs=rows;col.data.str=input;
    d->data.vars=&col;d->data.nvars=1;d->data.nobs=rows;d->data.sort_order=order;d->obs_map=map;
    return STATA_OK;
}
stata_retcode ctools_sort_dispatch(stata_data *d,int *v,size_t n,sort_algorithm_t a){(void)d;(void)v;(void)n;(void)a;return STATA_OK;}
stata_retcode ctools_store_filtered_rowpar(double *v,size_t n,int j,perm_idx_t *m){
    for(size_t i=0;i<n;i++)writevalue(j,m[i],v[i]);return STATA_OK;
}
int main(void){
    host.nobs1=one;host.nobs2=last;host.nobs=last;host.nvars=vars;
    host.isstr=isstring;host.isstrl=islong;host.ismissing=ismissing;
    host.sdata=readstring;host.safestore=writevalue;host.store=writevalue;
    host.scalsave=scalar;host.spouterr=message;host.macuse=getmacro;host.missval=8.98846567431158e307;
    const char *args[]={"1 2", "1 2 stream", "1 2 stream labelmap=9", "1 2 stream labelmap=9 noextend"};
    for(int mode=0;mode<4;mode++)for(int shape=0;shape<3;shape++){
        rows=shape==0?0:3;input[0]=shape==1?"":"a";input[1]=shape==1?"":"b";input[2]=shape==1?"":"c";
        fail_scalar=scalar_calls=0;assert(cencode_main(args[mode])==0);
        int total=scalar_calls;assert(total>=10);
        for(int f=1;f<=total;f++){
            fail_scalar=f;scalar_calls=0;assert(cencode_main(args[mode])==459);
        }
        fail_scalar=scalar_calls=0;assert(cencode_main(args[mode])==0);
    }
    rows=3;input[0]="a";input[1]="b";input[2]="c";
    fail_scalar=1;scalar_calls=0;assert(cgroup_main("nkeys=1 out skipif")==459);
    fail_scalar=0;scalar_calls=0;assert(cgroup_main("nkeys=1 out skipif")==0);
    assert(output[0]==1&&output[1]==2&&output[2]==3);
    return 0;
}
#endif
#if defined(CTOOLS_CASE_astra_parquet)
#include <assert.h>
#include "io/parquet_decode.c"
/* The page decoder test must never enter file/page I/O or decompression. */
int pq_page_header(const uint8_t *p, size_t n, pq_page *h, size_t *used) {
    (void)p; (void)n; (void)h; (void)used; abort();
}
const char *pq_decompress(int codec, const uint8_t *src, size_t n, uint8_t *dst, size_t out) {
    (void)codec; (void)src; (void)n; (void)dst; (void)out; abort();
}
int main(void) {
    /* Independently encoded DELTA_BINARY_PACKED headers: block=128,
       miniblocks=4, count, first=0; no delta payload for count <= 1. */
    for (int width = 0; width <= 1; width++) {
        for (int count = 0; count <= 3; count++) {
            uint8_t bad[] = {0x80, 1, 4, (uint8_t)count, 0};
            pq_vals v = {0};
            if (count == 2) continue;
            assert(delta_strings(bad, sizeof(bad), &v, 2, 1, width) != 0);
            vals_free(&v);
        }
        /* One value, zero prefix, one-byte suffix 'a'. */
        const uint8_t good[] = {0x80,1,4,1,0, 0x80,1,4,1,2, 'a'};
        pq_vals v = {0};
        assert(delta_strings(good, sizeof(good), &v, 1, 1, width) == 0);
        assert(v.count == 1);
        assert(width ? v.fixed[0] == 'a' : v.strs[0].n == 1 && v.strs[0].p[0] == 'a');
        vals_free(&v);
    }
    pq_values out = {0};
    pq_field fd = {0};
    pq_vals absent = {0};
    fd.arrow = PA_UTF8; fd.physical = PQ_BYTE_ARRAY;
    assert(emit(&out, &fd, SIZE_MAX, &absent, 0) != 0);
    absent.count = 1;
    assert(emit(&out, &fd, SIZE_MAX, &absent, 0) != 0);
    fd.physical = PQ_FLBA;
    assert(emit(&out, &fd, SIZE_MAX, &absent, 0) != 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_astra_merge)
#include <assert.h>
#include <string.h>
#include "cmerge/cmerge_io.c"
cmerge_using_cache_t g_using_cache;
static ST_plugin spi;
ST_plugin *_stata_ = &spi;
static int fail_read, fail_store, calls;
static ST_retcode readnum(ST_int v, ST_int i, ST_double *out) {
    (void)v; calls++;
    if (i == fail_read) return 459;
    *out = 0; return 0;
}
static ST_retcode store(ST_int v, ST_int i, ST_double x) {
    (void)v; (void)x; calls++; return i == fail_store ? 459 : 0;
}
static ST_retcode readstr(ST_int v, ST_int i, char *out) {
    double d; int rc = readnum(v,i,&d); if (!rc) strcpy(out, ""); return rc;
}
static ST_retcode strstore(ST_int v, ST_int i, char *s) { (void)s; return store(v,i,0); }
static ST_boolean missing(ST_double x) { return x >= spi.missval; }
int main(void) {
    spi.safevdata = readnum; spi.vdata = readnum; spi.safestore = store;
    spi.sdata = readstr; spi.sstore = strstore; spi.missval = 8.98846567431158e307;
    spi.ismissing = missing;
    double nums[] = {1,2,3}; char *strs[] = {"a","b","c"};
    stata_variable column = {0}; column.nobs = 3;
    g_using_cache.keepusing.data.vars = &column;
    cmerge_output_spec_t specs[3] = {0};
    for (int i=0;i<3;i++) { specs[i].using_sorted_row=i; specs[i].merge_result=MERGE_RESULT_BOTH; }
    for (int string=0;string<2;string++) {
        column.type=string ? STATA_TYPE_STRING : STATA_TYPE_DOUBLE;
        if (string) column.data.str=strs; else column.data.dbl=nums;
        for (int read=0;read<2;read++) for (int row=1;row<=3;row++) {
            cmerge_keepusing_write_args_t a = {0};
            a.dest_idx=1; a.output_nobs=3; a.specs=specs;
            a.is_shared=read; a.replace_mode=1;
            fail_read=read ? row : 0; fail_store=read ? 0 : row; calls=0;
            assert(cmerge_write_keepusing_var_thread(&a) != NULL);
            assert(!a.success && a.error == 459);
            fail_read=fail_store=0;
            assert(cmerge_write_keepusing_var_thread(&a) == NULL && a.success && !a.error);
        }
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_astra_stream)
#include <assert.h>
#include <stdatomic.h>
#include "csort/csort_stream.c"
static ST_plugin spi;
ST_plugin *_stata_ = &spi;
static double numbers[4][65];
static char strings[2][65][16];
static int fail_var, fail_row, fail_write, status=459;
static atomic_int read_count, store_count;
int ctools_get_max_threads(void) { return 2; }
double ctools_timer_seconds(void) { return 0; }
/* The transfer entry point under test must not call the key-loader/sorter. */
void stata_data_init(stata_data *d) { (void)d; abort(); }
void ctools_filtered_data_init(ctools_filtered_data *d) { (void)d; abort(); }
void ctools_filtered_data_free(ctools_filtered_data *d) { (void)d; abort(); }
stata_retcode ctools_data_load(ctools_filtered_data *d,int *v,size_t n,size_t a,size_t b,int f) {
    (void)d;(void)v;(void)n;(void)a;(void)b;(void)f;abort();
}
stata_retcode ctools_sort_dispatch(stata_data *d,int *v,size_t n,sort_algorithm_t a) {
    (void)d;(void)v;(void)n;(void)a;abort();
}
stata_retcode ctools_apply_permutation(stata_data *d) { (void)d;abort(); }
stata_retcode ctools_store_filtered(double *v,size_t n,int j,perm_idx_t *m) {
    (void)v;(void)n;(void)j;(void)m;abort();
}
stata_retcode ctools_store_filtered_str(char **v,size_t n,int j,perm_idx_t *m) {
    (void)v;(void)n;(void)j;(void)m;abort();
}
static ST_boolean isstr(ST_int v) { return v>=5; }
static ST_retcode getnum(ST_int v,ST_int i,double *out) {
    atomic_fetch_add(&read_count,1);
    if (!fail_write && v==fail_var && i==fail_row) return status;
    *out=numbers[v-1][i-1]; return 0;
}
static ST_retcode putnum(ST_int v,ST_int i,double x) {
    atomic_fetch_add(&store_count,1);
    if (fail_write && v==fail_var && i==fail_row) return status;
    numbers[v-1][i-1]=x; return 0;
}
static ST_retcode getstr(ST_int v,ST_int i,char *out) {
    atomic_fetch_add(&read_count,1);
    if (!fail_write && v==fail_var && i==fail_row) return status;
    strcpy(out,strings[v-5][i-1]); return 0;
}
static ST_retcode putstr(ST_int v,ST_int i,char *s) {
    atomic_fetch_add(&store_count,1);
    if (fail_write && v==fail_var && i==fail_row) return status;
    strcpy(strings[v-5][i-1],s); return 0;
}
static void reset_data(void) {
    for(int v=0;v<4;v++) for(int i=0;i<65;i++) numbers[v][i]=v*100+i+1;
    for(int v=0;v<2;v++) for(int i=0;i<65;i++) snprintf(strings[v][i],16,"%d:%d",v,i+1);
    atomic_store(&read_count,0); atomic_store(&store_count,0);
}
int main(void) {
    spi.isstr=isstr; spi.safevdata=getnum; spi.safestore=putnum;
    spi.sdata=getstr; spi.sstore=putstr;
    int vars[]={1,2,3,4,5,6}, widths[]={0,0,0,0,15,15};
    perm_idx_t perm[33],map[33];
    for(int i=0;i<33;i++) {perm[i]=32-i;map[i]=2*i+1;}
    const int rows[]={1,33,65};
    for(int batch=1;batch<=2;batch++) {
        for(int v=1;v<=6;v++) for(int w=0;w<2;w++) for(int r=0;r<3;r++) {
            reset_data(); fail_var=v;fail_row=rows[r];fail_write=w;
            stata_retcode rc=csort_stream_apply_permutation(perm,33,vars,6,map,0,batch,widths,NULL);
            assert(rc == (w ? STATA_ERR_STATA_WRITE : STATA_ERR_STATA_READ));
        }
        reset_data();fail_var=0;
        assert(csort_stream_apply_permutation(perm,33,vars,6,map,0,batch,widths,NULL)==STATA_OK);
        for(int v=0;v<4;v++) for(int i=0;i<65;i++)
            assert(numbers[v][i] == v*100 + (i%2 ? i+1 : 65-i));
        for(int v=0;v<2;v++) for(int i=0;i<65;i++) {
            char expected[16];snprintf(expected,16,"%d:%d",v,i%2 ? i+1 : 65-i);
            assert(!strcmp(strings[v][i],expected));
        }
    }
    reset_data();fail_var=5;fail_row=33;fail_write=0;
    assert(csort_stream_apply_permutation(perm,33,vars,6,map,0,1,NULL,NULL)!=STATA_OK);
    assert(atomic_load(&store_count)==0); /* Width preflight failed before writing. */
    reset_data();fail_var=1;fail_row=1;fail_write=0;status=1;
    assert(csort_stream_apply_permutation(perm,33,vars,6,map,0,1,widths,NULL)==STATA_ERR_CANCELLED);
    return 0;
}
#endif
