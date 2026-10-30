/* Sort contracts tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_order_parallel_source)
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_order.h"
#include "ctools_config.h"
static int threads = 1, keys[] = {0,1};
static stata_data data;
static size_t *rank_in;
int ctools_get_max_threads(void) { return threads; }
static int reference(const void *pa, const void *pb) {
    size_t a=*(const perm_idx_t*)pa, b=*(const perm_idx_t*)pb;
    for(int k=0;k<2;k++) {
        int c;
        if(data.vars[k].type==STATA_TYPE_STRING) c=strcmp(data.vars[k].data.str[a],data.vars[k].data.str[b]);
        else {double x=data.vars[k].data.dbl[a],y=data.vars[k].data.dbl[b];c=(x>y)-(x<y);}
        if(c)return c;
    }
    return (rank_in[a]>rank_in[b])-(rank_in[a]<rank_in[b]);
}
static unsigned long long seed=923842;
static unsigned rnd(void) { seed=seed*6364136223846793005ULL+1;return seed>>32; }
static void check(size_t n,int strings) {
    data.nobs=n;data.nvars=2;data.vars=calloc(2,sizeof(*data.vars));
    data.sort_order=malloc((n+1)*sizeof(perm_idx_t));
    rank_in=malloc((n+1)*sizeof(*rank_in));
    perm_idx_t *expected=malloc((n+1)*sizeof(*expected)), *input=malloc((n+1)*sizeof(*input));
    char (*text)[24]=malloc((n+1)*24);
    for(int k=0;k<2;k++) {data.vars[k].type=STATA_TYPE_DOUBLE;data.vars[k].data.dbl=malloc((n+1)*sizeof(double));}
    if(strings){free(data.vars[0].data.dbl);data.vars[0].type=STATA_TYPE_STRING;data.vars[0].data.str=malloc((n+1)*sizeof(char*));}
    for(size_t i=0;i<n;i++) {
        double value=(int)(rnd()%1024)-512;
        if(i%17==0)value=0;
        if(i%19==0)value=-0.0;
        if(i%23==0) {value=8.98846567431158e307;for(unsigned j=0;j<i%27;j++)value=nextafter(value,INFINITY);}
        if(strings){snprintf(text[i],24,"group%u",rnd()%101);data.vars[0].data.str[i]=text[i];}
        else data.vars[0].data.dbl[i]=value;
        data.vars[1].data.dbl[i]=(int)(rnd()%13)-6;
        input[i]=(perm_idx_t)i;
    }
    for(size_t i=n;i>1;i--){size_t j=rnd()%i;perm_idx_t t=input[i-1];input[i-1]=input[j];input[j]=t;}
    for(size_t i=0;i<n;i++)rank_in[input[i]]=i;
    memcpy(expected,input,n*sizeof(*expected));qsort(expected,n,sizeof(*expected),reference);
    int counts[]={1,2,12};
    for(int t=0;t<3;t++) {
        threads=counts[t];memcpy(data.sort_order,input,n*sizeof(*input));
        assert(ctools_order_stable(&data,keys,2)==STATA_OK);
        assert(!memcmp(data.sort_order,expected,n*sizeof(*expected)));
        assert(ctools_order_stable(&data,keys,2)==STATA_OK);
        assert(!memcmp(data.sort_order,expected,n*sizeof(*expected)));
    }
    free(text);free(expected);free(input);free(rank_in);free(data.sort_order);
    free(data.vars[0].data.dbl);free(data.vars[1].data.dbl);free(data.vars);
}
int main(void) {
    size_t sizes[]={0,1,257,32767,32768,32769,200003};
    for(size_t s=0;s<sizeof(sizes)/sizeof(*sizes);s++)for(int strings=0;strings<2;strings++)check(sizes[s],strings);
    puts("PASS stable parallel order: numeric/string keys, signed zeros, extended missings, partial tiles, incoming tie order, threads 1/2/12");
}
#endif

#if defined(CTOOLS_CASE_p1_native_2)
#include <assert.h>
#include "ctools_ols.h"
__attribute__((noinline)) static double sum(const double *x, int n) {
    dd_real total = {0, 0};
    for (int i = 0; i < n; i++) total = dd_add_d(total, x[i]);
    return total.hi + total.lo;
}
int main(void) {
    volatile double a = 1e16, b = 1;
    dd_real pair = two_sum(a, b);
    assert(pair.lo == 1);
    double x[] = {1e16, 1, -1e16};
    assert(sum(x, 3) == 1);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_p1_native_3)
#include "ctools_types.c"
#include "ctools_platform.c"
int ctools_get_max_threads(void) { return 4; }
int main(void) {
    stata_variable vars[4] = {0};
    perm_idx_t perm[2] = {1, 0};
    stata_data data = {0};
    data.nvars = 4; data.nobs = 2; data.vars = vars; data.sort_order = perm;
    /* Two pointer arrays and two aligned buffers for each of four workers. */
    for (int failure = 1; failure <= 10; failure++) {
        calls = 0; fail_at = failure;
        assert(ctools_apply_permutation(&data) == STATA_ERR_MEMORY);
        assert(live == 0);
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_7)
#include "ctools_types.h"
int ctools_get_max_threads(void) {return 4;}
#include SORT_IMPLEMENTATION
int main(void) {
 perm_idx_t order[128];uint64_t keys[128];char *strings[128];
 for(int i=0;i<128;i++){keys[i]=128-i;strings[i]=i%2?"az":"ba";}
 for(int string=0;string<2;string++) {
  calls=0;fail_at=0;for(int i=0;i<128;i++)order[i]=i;
  assert((string?STRING_CALL:NUMERIC_CALL)==0);int allocations=calls;assert(live_count==0);
  for(int f=1;f<=allocations;f++) {
   calls=0;fail_at=f;for(int i=0;i<128;i++)order[i]=i;
   assert((string?STRING_CALL:NUMERIC_CALL)!=0);assert(live_count==0);
  }
 }
 return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_sort)
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_types.h"
static ST_plugin mock;ST_plugin *_stata_=&mock;
static ST_boolean missing(ST_double x){return x>=8.98846567431158e307;}
int ctools_get_max_threads(void) {return 4;}
#include SORT_IMPLEMENTATION
int main(void) {
 mock.ismissing=missing;mock.missval=8.98846567431158e307;
 const size_t N=400001;int key=1;
 stata_variable var={.type=STATA_TYPE_STRING,.nobs=N};
 char **s=malloc(N*sizeof(char*));double *x=malloc(N*sizeof(double));
 perm_idx_t *order=malloc(N*sizeof(perm_idx_t));unsigned char *seen=malloc(N);
 stata_data d={.nobs=N,.nvars=1,.vars=&var,.sort_order=order};
 for(int string=0;string<STRING_VARIANTS;string++) {
  for(size_t i=0;i<N;i++){s[i]=i%7==0?"":(i%3==0?"a":(i%2?"az":"ba"));x[i]=(double)((N-i)%2345);order[i]=i;}
  var.type=string?STATA_TYPE_STRING:STATA_TYPE_DOUBLE;if(string)var.data.str=s;else var.data.dbl=x;
  assert(SORT_ORDER(&d,&key,1)==0);memset(seen,0,N);
  for(size_t i=0;i<N;i++){assert(order[i]<N && !seen[order[i]]);seen[order[i]]=1;if(i){if(string)assert(strcmp(s[order[i-1]],s[order[i]])<=0);else assert(x[order[i-1]]<=x[order[i]]);}}
 }
 free(s);free(x);free(order);free(seen);return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_timsort)
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
int ctools_get_max_threads(void) { return 1; }
#include "ctools_sort_timsort.c"
int main(void) {
    double x[] = {3, 1, 2, 4};
    char *str[] = {"b", "a", "b", "a"};
    stata_variable vars[3] = {{.type = STATA_TYPE_DOUBLE, .nobs = 4},
                              {.type = STATA_TYPE_STRING, .nobs = 4},
                              {.type = STATA_TYPE_DOUBLE, .nobs = 4}};
    vars[0].data.dbl = x;
    vars[1].data.str = str;
    vars[2].data.dbl = x;
    perm_idx_t order[4];
    int keys[] = {1, 2, 3}, allocations = 0;
    stata_data data = {.nobs = 4, .nvars = 3, .vars = vars, .sort_order = order};
    for (int f = 0; f <= allocations; f++) {
        calls = 0;
        fail_at = f;
        for (int i = 0; i < 4; i++)
            order[i] = i;
        int rc = timsort_by_composite_keys(&data, keys, 3);
        if (!f) {
            assert(!rc);
            allocations = calls;
        } else
            assert(rc == STATA_ERR_MEMORY);
        assert(alive == 0);
    }
    return 0;
}
#endif
