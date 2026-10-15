
#include "ctools_types.c"
int ctools_get_max_threads(void) { return 4; }
int main(void) {
    stata_variable vars[4] = {0};
    perm_idx_t perm[2] = {1, 0};
    stata_data data = {0};
    data.nvars = 4; data.nobs = 2; data.vars = vars; data.sort_order = perm;
    /* Two pointer arrays and two aligned buffers for each of four workers. */
    for (int failure = 1; failure <= 10; failure++) {
        calls = 0; fail_at = failure;
        assert(ctools_apply_permutation(&data) == STATA_ERR_MEMORY);
        assert(live == 0);
    }
    return 0;
}
