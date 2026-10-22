#include "cpsmatch/cpsmatch_internal.h"
static int fail_growth;
static void *probe_realloc(void *p,size_t n) {if(fail_growth) return NULL; return realloc(p,n);}
#define realloc probe_realloc
#include "cpsmatch/cpsmatch_kernels.c"
#undef realloc
#include "cpsmatch/cpsmatch_options.c"
#include "cpsmatch/cpsmatch_impl.c"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static double tr[301], ps[301], outputs[7][302], matched;
static perm_idx_t map[301];
static stata_variable columns[2];
static ST_int one(void){return 1;} static ST_int count(void){return 301;}
static ST_retcode store(ST_int v,ST_int i,ST_double x){outputs[v][i]=x;return 0;}
static ST_retcode scal(char *name,double x){if(!strcmp(name,"_cpsmatch_n_matched")) matched=x;return 0;}
static ST_boolean missing(double x){return x>=8.98846567431158e307;}
static ST_retcode message(char *s){fputs(s,stderr);return 0;}
double ctools_timer_seconds(void){return 0;}
int ctools_get_max_threads(void){return 1;}
void ctools_filtered_data_init(ctools_filtered_data *fd){memset(fd,0,sizeof(*fd));}
void ctools_filtered_data_free(ctools_filtered_data *fd){(void)fd;}
stata_retcode ctools_data_load(ctools_filtered_data *fd,int *v,size_t nv,size_t s,size_t e,int flags){
 (void)v;(void)nv;(void)s;(void)e;(void)flags;
 columns[0].data.dbl=tr;columns[1].data.dbl=ps;
 fd->data.nobs=301;fd->data.nvars=2;fd->data.vars=columns;fd->obs_map=map;return STATA_OK;
}
int main(void){
 mock.ismissing=missing;mock.nobs1=one;mock.nobs2=count;mock.safestore=store;mock.store=store;mock.scalsave=scal;mock.spouterr=message;mock.missval=8.98846567431158e307;
 tr[0]=1;for(int i=0;i<301;i++){ps[i]=0.5;map[i]=i+1;}
 const char *args="1 2 0 3 4 5 301 2 1 0 0 0 1 1 0 0";
 for(fail_growth=0;fail_growth<2;fail_growth++){
 int rc=cpsmatch_main(args);printf("fail_growth=%d rc=%d matched=%.0f treated_support=%.0f\n",fail_growth,rc,matched,outputs[5][1]);
 }
}
