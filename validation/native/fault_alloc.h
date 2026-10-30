
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static int calls, fail_at, live_count;
static void *live[10000];
static void record(void *p) { if(p) { for(int i=0;i<10000;i++) if(!live[i]) {live[i]=p;live_count++;return;} abort(); } }
static void *a_malloc(size_t n) { if(++calls==fail_at)return NULL;void *p=malloc(n);if(p)memset(p,0xa5,n);record(p);return p; }
static void *a_calloc(size_t n,size_t s) {if(++calls==fail_at)return NULL;void *p=calloc(n,s);record(p);return p;}
static int a_align(void **p,size_t a,size_t n) {if(++calls==fail_at)return ENOMEM;int rc=posix_memalign(p,a,n);if(!rc)record(*p);return rc;}
static void a_free(void *p) {if(!p)return;for(int i=0;i<10000;i++)if(live[i]==p){live[i]=NULL;live_count--;free(p);return;}assert(!"invalid free");}
#define malloc a_malloc
#define calloc a_calloc
#define posix_memalign a_align
#define free a_free
