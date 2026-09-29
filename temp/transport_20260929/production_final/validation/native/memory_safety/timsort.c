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
