#include "ctools_threads.c"
#include <assert.h>
#include <limits.h>
int main(void) {
    ctools_set_max_threads(INT_MAX);
    assert(ctools_get_max_threads() == CTOOLS_MAX_THREADS);
    assert(ctools_workspace_threads(64, 1, 80000000) == 1);
    assert(ctools_workspace_threads(64, 100, 80000000) == 3);
    assert(ctools_workspace_threads(64, 100, 0) == 64);
    ctools_set_max_threads(2);
    assert(ctools_get_max_threads() == 2);
    ctools_reset_max_threads();
    assert(ctools_get_max_threads() > 0);
    return 0;
}
