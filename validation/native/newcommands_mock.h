/* Dynamic SPI fixtures with thread-safe allocation and publication injection. */
#undef NDEBUG
#include <assert.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "stplugin.h"
#include "ctools_types.h"

static ST_plugin mock;
ST_plugin *_stata_=&mock;
enum { COLS=16 };
static int nr,nv,threads=1,types[COLS];
static double *numbers[COLS];
static char **strings[COLS];
static unsigned char *selected;
static size_t capacity;
static atomic_int allocations,live,stores,worker_store;
static int fail_alloc,fail_store,fail_macro;
static char widths[16000],nout[64],nmissing[64];
static const char *delimiters[]={"::",":"};

static void text_cell(int col,int row,const char *s)
{
    size_t n=strlen(s)+1;
    char *p=realloc(strings[col][row],n);assert(p);
    strings[col][row]=p;memcpy(p,s,n);
}
static const char *string_cell(int col,int row)
{ return strings[col][row] ? strings[col][row] : ""; }
static void resize_fixture(int rows)
{
    if ((size_t)rows>capacity) {
        for(int c=0;c<COLS;c++) {
            numbers[c]=realloc(numbers[c],(size_t)rows*sizeof(double));
            strings[c]=realloc(strings[c],(size_t)rows*sizeof(char*));
            assert(numbers[c] && strings[c]);
            for(size_t r=capacity;r<(size_t)rows;r++) { numbers[c][r]=mock.missval;strings[c][r]=NULL; }
        }
        selected=realloc(selected,(size_t)rows);assert(selected);
        memset(selected+capacity,1,(size_t)rows-capacity);
        capacity=rows;
    }
    nr=rows;
}
static void clear_data(void)
{
    for(int c=0;c<COLS;c++) {
        for(size_t r=0;r<capacity;r++) free(strings[c][r]);
        free(strings[c]);free(numbers[c]);strings[c]=NULL;numbers[c]=NULL;
    }
    free(selected);selected=NULL;capacity=0;nr=nv=0;
    memset(types,0,sizeof(types));
}
static void clear_fixture(void)
{
    clear_data();allocations=stores=worker_store=0;fail_alloc=fail_store=fail_macro=0;
    widths[0]=nout[0]=nmissing[0]=0;
}
typedef union { max_align_t align; struct { uint64_t magic; } info; } tracked_header;
static void *tracked_malloc(size_t n)
{
    if(atomic_fetch_add(&allocations,1)+1==fail_alloc) return NULL;
    if(n>SIZE_MAX-sizeof(tracked_header)) return NULL;
    tracked_header *h=malloc(sizeof(*h)+n);if(!h)return NULL;
    h->info.magic=UINT64_C(0x415354524136);atomic_fetch_add(&live,1);return h+1;
}
static void *tracked_calloc(size_t n,size_t s)
{
    if(s && n>SIZE_MAX/s)return NULL;
    void *p=tracked_malloc(n*s);if(p)memset(p,0,n*s);return p;
}
static void tracked_free(void *p)
{
    if(!p)return;
    tracked_header *h=(tracked_header*)p-1;
    assert(h->info.magic==UINT64_C(0x415354524136));h->info.magic=0;
    atomic_fetch_sub(&live,1);free(h);
}
static int store_failure(void)
{
#ifdef _OPENMP
    if(omp_get_thread_num()>0)atomic_store(&worker_store,1);
#endif
    return atomic_fetch_add(&stores,1)+1==fail_store;
}
static ST_int nobs(void){return nr;}
static ST_int nvars(void){return nv;}
static ST_boolean isstr(ST_int v){assert(v>0&&v<=nv);return types[v-1];}
static ST_boolean isstrl(ST_int v){(void)v;return 0;}
static ST_retcode putnum(ST_int v,ST_int r,double x)
{
    assert(v>0&&v<=nv&&r>0&&r<=nr);if(store_failure())return 459;
    numbers[v-1][r-1]=x;return 0;
}
static ST_retcode putstr(ST_int v,ST_int r,char *x)
{
    assert(v>0&&v<=nv&&r>0&&r<=nr&&strlen(x)<=2045);
    if(store_failure())return 459;text_cell(v-1,r-1,x);return 0;
}
static ST_retcode macro_save(char *key,char *value)
{
    if(fail_macro)return 459;
    if(!strcmp(key,"_csplit_widths"))strcpy(widths,value);
    else if(!strcmp(key,"_crangejoin_nout"))strcpy(nout,value);
    else if(!strcmp(key,"_cinterpolate_nmissing"))strcpy(nmissing,value);
    else abort();return 0;
}
static ST_retcode macro_get(char *key,char *value,ST_int size)
{
    int d=0;assert(sscanf(key,"_csplit_delim%d",&d)==1&&d>=1&&d<=2);
    assert(size>(int)strlen(delimiters[d-1]));strcpy(value,delimiters[d-1]);return 0;
}
int ctools_get_max_threads(void){return threads;}
double ctools_timer_seconds(void){return 0;}
void ctools_verbose(const char *m,int v,const char *f,...){(void)m;(void)v;(void)f;}
void ctools_error(const char *m,const char *f,...){(void)m;(void)f;}

#define malloc tracked_malloc
#define calloc tracked_calloc
#define free tracked_free
void ctools_filtered_data_init(ctools_filtered_data *fd){memset(fd,0,sizeof(*fd));}
void ctools_filtered_data_free(ctools_filtered_data *fd)
{
    for(size_t c=0;fd->data.vars&&c<fd->data.nvars;c++) {
        stata_variable *v=fd->data.vars+c;
        if(v->type==STATA_TYPE_STRING) {
            for(size_t r=0;v->data.str&&r<v->nobs;r++)free(v->data.str[r]);
            free(v->data.str);
        } else free(v->data.dbl);
    }
    free(fd->data.vars);free(fd->data.sort_order);free(fd->obs_map);memset(fd,0,sizeof(*fd));
}
stata_retcode ctools_data_load(ctools_filtered_data *fd,int *vars,size_t n,size_t first,size_t last,int flags)
{
    (void)first;(void)last;if(!vars)n=nv;
    size_t rows=0;
    for(int r=0;r<nr;r++)if((flags&CTOOLS_LOAD_SKIP_IF)||selected[r])rows++;
    fd->data.vars=calloc(n,sizeof(stata_variable));if(!fd->data.vars)return STATA_ERR_MEMORY;
    fd->data.nvars=n;fd->data.nobs=rows;
    fd->data.sort_order=malloc(rows*sizeof(perm_idx_t));fd->obs_map=malloc(rows*sizeof(perm_idx_t));
    if(!fd->data.sort_order||!fd->obs_map)return STATA_ERR_MEMORY;
    for(int r=0,i=0;r<nr;r++)if((flags&CTOOLS_LOAD_SKIP_IF)||selected[r]) {
        fd->data.sort_order[i]=i;fd->obs_map[i++]=r+1;
    }
    for(size_t c=0;c<n;c++) {
        int from=vars?vars[c]-1:(int)c;stata_variable *v=fd->data.vars+c;
        v->nobs=rows;v->type=types[from]?STATA_TYPE_STRING:STATA_TYPE_DOUBLE;
        if(types[from]) {
            v->data.str=calloc(rows,sizeof(char*));if(!v->data.str)return STATA_ERR_MEMORY;
            for(size_t i=0;i<rows;i++) {
                const char *s=string_cell(from,fd->obs_map[i]-1);size_t len=strlen(s)+1;
                v->data.str[i]=malloc(len);if(!v->data.str[i])return STATA_ERR_MEMORY;
                memcpy(v->data.str[i],s,len);
            }
        } else {
            v->data.dbl=malloc(rows*sizeof(double));if(!v->data.dbl)return STATA_ERR_MEMORY;
            for(size_t i=0;i<rows;i++)v->data.dbl[i]=numbers[from][fd->obs_map[i]-1];
        }
    }
    return STATA_OK;
}
stata_retcode ctools_store_filtered(double *x,size_t n,int col,perm_idx_t *map)
{
    for(size_t i=0;i<n;i++)if(putnum(col,map[i],x[i]))return STATA_ERR_STATA_WRITE;
    return STATA_OK;
}
