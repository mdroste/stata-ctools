#include "alloc.h"
#include "stplugin.h"
#include <math.h>
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
#include "creghdfe/creghdfe_utils.c"
int main(void) {
    double cases[][6] = {{1, 1, 2, 2, -1, -1},
                         {1.5, 1.5, 2.5, 2.5, -1, -1},
                         {1e30, 1e30, -1e30, -1e30, 0, 0},
                         {-0x1p63, 0, 0x1p62, 1, -0x1p63, 1},
                         {0x1p63, 0x1p63, -0x1p64, 2, 2, -0x1p64}};
    double weights[] = {1, 2, 3, 4, 5, 6};
    for (size_t c = 0; c < sizeof(cases) / sizeof(*cases); c++) {
        for (int weighted = 0; weighted < 2; weighted++) {
            int allocations = 0;
            for (int f = 0; f <= allocations; f++) {
                ST_int levels[6], nlevels = 123, *counts = (void *)1;
                double *wcounts = (void *)1;
                calls = 0;
                fail_at = f;
                int rc = remap_and_count(cases[c], 6, levels, &nlevels, &counts,
                                         weighted ? weights : NULL, &wcounts);
                if (!f)
                    allocations = calls;
                if (f) {
                    assert(rc && nlevels == 0 && !counts && !wcounts);
                } else {
                    assert(!rc && nlevels > 0 && counts);
                    int total = 0;
                    double sum = 0;
                    for (int j = 0; j < nlevels; j++) {
                        total += counts[j];
                        if (weighted)
                            sum += wcounts[j];
                    }
                    assert(total == 6 && (!weighted || sum == 21));
                    for (int i = 0; i < 6; i++)
                        for (int j = 0; j < 6; j++)
                            assert((levels[i] == levels[j]) == (cases[c][i] == cases[c][j]));
                    free(counts);
                    free(wcounts);
                }
                assert(alive == 0);
            }
        }
    }
    return 0;
}
