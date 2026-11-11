/* Bounded, fault-injectable allocations for adversarial file tests. */
#undef NDEBUG
#include <assert.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
static atomic_int pq_allocations,pq_live;
static atomic_size_t pq_live_bytes;
static int pq_fail_at,pq_denied;
typedef union { max_align_t align; struct { size_t n; uint64_t magic; } info; } pq_alloc_header;
static void *pq_test_malloc(size_t n)
{
    if(atomic_fetch_add(&pq_allocations,1)+1==pq_fail_at)return NULL;
    if(n>16*1024*1024 || atomic_load(&pq_live_bytes)>64*1024*1024-n) { pq_denied=1;return NULL; }
    pq_alloc_header *p=malloc(sizeof(*p)+n+64);if(!p)return NULL;
    p->info.n=n;p->info.magic=UINT64_C(0x50415251554554);
    memset((unsigned char *)(p+1)+n,0xa7,64);
    atomic_fetch_add(&pq_live,1);atomic_fetch_add(&pq_live_bytes,n);return p+1;
}
static void pq_test_free(void *p)
{
    if(!p)return;pq_alloc_header *h=(pq_alloc_header*)p-1;
    assert(h->info.magic==UINT64_C(0x50415251554554));
    for(int i=0;i<64;i++)assert(((unsigned char *)p)[h->info.n+i]==0xa7);
    h->info.magic=0;
    atomic_fetch_sub(&pq_live,1);atomic_fetch_sub(&pq_live_bytes,h->info.n);free(h);
}
static void *pq_test_calloc(size_t n,size_t s)
{
    if(s&&n>SIZE_MAX/s){pq_denied=1;return NULL;}
    void *p=pq_test_malloc(n*s);if(p)memset(p,0,n*s);return p;
}
static void *pq_test_realloc(void *p,size_t n)
{
    void *out=pq_test_malloc(n);if(!out)return NULL;
    if(p) { pq_alloc_header *h=(pq_alloc_header*)p-1;
        memcpy(out,p,n<h->info.n?n:h->info.n);pq_test_free(p); }
    return out;
}
#define malloc pq_test_malloc
#define calloc pq_test_calloc
#define realloc pq_test_realloc
#define free pq_test_free
