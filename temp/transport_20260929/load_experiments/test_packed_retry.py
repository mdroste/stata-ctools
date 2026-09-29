from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
# Reuse the independent value model and compilation driver, inserting precise
# arena-allocation failure injection before source inclusion and final checks.
base=Path(__file__).with_name('test_load_edges.py').read_text()
base=base.replace('#include "ctools_data_io.c"',r'''#include "ctools_arena.h"
static atomic_int fail_large_arena, failed_large_seen, retry_seen;
static size_t expected_retry;
static ctools_string_arena *test_arena_create(size_t capacity,ctools_string_arena_mode mode) {
    if(capacity==expected_retry) atomic_fetch_add(&retry_seen,1);
    if(capacity>1000000 && atomic_exchange(&fail_large_arena,0)) {
        atomic_fetch_add(&failed_large_seen,1);return NULL;
    }
    return ctools_string_arena_create(capacity,mode);
}
#define ctools_string_arena_create test_arena_create
#include "ctools_data_io.c"
#undef ctools_string_arena_create''',1)
base=base.replace('    ctools_destroy_global_pool();\n    return 0;\n}',r'''    ctools_destroy_global_pool();
    shape=2;
    for(ncols=1;ncols<=4;ncols+=3) for(filter=0;filter<=1;filter++) {
      size_t selected_count=0;
      for(int row=3;row<=nrows-2;row++) selected_count+=selected(row)!=0;
      expected_retry=selected_count*64;
      fail_large_arena=1;failed_large_seen=0;retry_seen=0;
      assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_CHECK_IF)==STATA_OK);
      assert(failed_large_seen==1 && retry_seen==1);
      size_t out=0;
      for(int row=3;row<=nrows-2;row++) if(selected(row)) {
        char value[2046];
        for(int j=1;j<=ncols;j++) {
          wide_text(j,row,value);
          assert(!strcmp(fd.data.vars[j-1].data.str[out],value));
        }
        out++;
      }
      assert(out==fd.data.nobs);ctools_filtered_data_free(&fd);
    }
    ctools_destroy_global_pool();
    return 0;
}''',1)
exec(compile(base,str(Path(__file__).with_name('test_load_edges.py')),'exec'))
