
#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include "ctools_threads.h"
static int calls, fail_at;
static int create(pthread_t *t,const pthread_attr_t *a,void *(*f)(void*),void *p) {
 if(++calls==fail_at)return EAGAIN;return pthread_create(t,a,f,p);
}
#define pthread_create create
#include "ctools_threads.c"
int main(void) {
 for(int run=0;run<100;run++)for(int f=1;f<=4;f++) {
  ctools_persistent_pool pool;calls=0;fail_at=f;
  assert(ctools_persistent_pool_init(&pool,4)!=0);
  assert(!pool.initialized && !pool.workers && pool.num_workers==0);
  fail_at=0;assert(ctools_persistent_pool_init(&pool,2)==0);
  ctools_persistent_pool_destroy(&pool);
 }
 return 0;
}
