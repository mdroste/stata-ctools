
#include <assert.h>
#include "ctools_ols.h"
__attribute__((noinline)) static double sum(const double *x, int n) {
    dd_real total = {0, 0};
    for (int i = 0; i < n; i++) total = dd_add_d(total, x[i]);
    return total.hi + total.lo;
}
int main(void) {
    volatile double a = 1e16, b = 1;
    dd_real pair = two_sum(a, b);
    assert(pair.lo == 1);
    double x[] = {1e16, 1, -1e16};
    assert(sum(x, 3) == 1);
    return 0;
}
