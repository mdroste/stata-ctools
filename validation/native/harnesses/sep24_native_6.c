
#include "ctools_arena.c"
int main(void) {
 ctools_arena a;ctools_arena_init(&a,32);
 size_t capacity=0;
 for(int run=0;run<50;run++) {
  assert(ctools_arena_alloc(&a,24));assert(ctools_arena_alloc(&a,64));
  void *aligned=ctools_arena_alloc_aligned(&a,53,64);assert(aligned && (uintptr_t)aligned%64==0);
  assert(!ctools_arena_alloc(&a,SIZE_MAX));assert(!ctools_arena_alloc_aligned(&a,SIZE_MAX,64));
  if(run==0)capacity=live_count;else assert(capacity==(size_t)live_count);
  ctools_arena_reset(&a);
 }
 ctools_arena_free(&a);assert(live_count==0);return 0;
}
