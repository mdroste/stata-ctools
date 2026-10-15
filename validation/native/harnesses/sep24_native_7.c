
#include "ctools_types.h"
int ctools_get_max_threads(void) {return 4;}
#include "ctools_sort_ALGORITHM.c"
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
