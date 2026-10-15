
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_types.h"
static ST_plugin mock;ST_plugin *_stata_=&mock;
static ST_boolean missing(ST_double x){return x>=8.98846567431158e307;}
int ctools_get_max_threads(void) {return 4;}
#include "ctools_sort_ALGORITHM.c"
int main(void) {
 mock.ismissing=missing;mock.missval=8.98846567431158e307;
 const size_t N=400001;int key=1;
 stata_variable var={.type=STATA_TYPE_STRING,.nobs=N};
 char **s=malloc(N*sizeof(char*));double *x=malloc(N*sizeof(double));
 perm_idx_t *order=malloc(N*sizeof(perm_idx_t));unsigned char *seen=malloc(N);
 stata_data d={.nobs=N,.nvars=1,.vars=&var,.sort_order=order};
 for(int string=0;string<2;string++) {
  for(size_t i=0;i<N;i++){s[i]=i%7==0?"":(i%3==0?"a":(i%2?"az":"ba"));x[i]=(double)((N-i)%2345);order[i]=i;}
  var.type=string?STATA_TYPE_STRING:STATA_TYPE_DOUBLE;if(string)var.data.str=s;else var.data.dbl=x;
  assert(ctools_sort_ALGORITHM_order_only(&d,&key,1)==0);memset(seen,0,N);
  for(size_t i=0;i<N;i++){assert(order[i]<N && !seen[order[i]]);seen[order[i]]=1;if(i){if(string)assert(strcmp(s[order[i-1]],s[order[i]])<=0);else assert(x[order[i-1]]<=x[order[i]]);}}
 }
 free(s);free(x);free(order);free(seen);return 0;
}
