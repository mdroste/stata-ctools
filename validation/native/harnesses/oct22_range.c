#include <stdio.h>
#include <assert.h>
#include "crangestat/crangestat_kernels.c"
ST_plugin mock;
ST_plugin *_stata_ = &mock;
int main(void) {
 mock.missval = 8.98846567431158e307;
 double data[] = {1e16,1,2,3};
 double sum[5],sq[5]; size_t count[5];
 prefix_arrays pa = {.prefix_sum=sum,.prefix_sum2=sq,.prefix_count=count};
 crangestat_build_prefix_arrays(data,0,4,&pa);
 for(size_t i=0;i<4;i++) {
  assert(crangestat_compute_stat_with_prefix(STAT_MEAN,&pa,data,0,i,i+1,i,0)==data[i]);
  assert(crangestat_compute_stat_with_prefix(STAT_SUM,&pa,data,0,i,i+1,i,0)==data[i]);
 }
 assert(crangestat_compute_stat_with_prefix(STAT_SUM,&pa,data,0,1,4,2,0)==6);
 assert(crangestat_compute_stat_with_prefix(STAT_MEAN,&pa,data,0,1,4,2,1)==2);
 puts("oct22 prefix cancellation passed");
}
