
#include "cbinscatter/cbinscatter_impl.c"
static ST_retcode scalars(char *s, ST_double d) { (void)s; (void)d; return 0; }
static ST_retcode matrix(char *s, ST_int row, ST_int col, ST_double d) {
    (void)s; (void)row; (void)d; return col == 4 ? 503 : 0;
}
int main(void) {
    setup(); mock.scalsave = scalars; mock.safematstore = matrix;
    BinStats bin = {.n_obs = 3};
    ByGroupResult group = {.num_bins = 1, .bins = &bin};
    BinscatterResults results = {.num_by_groups = 1, .groups = &group};
    BinscatterConfig config = {0};
    assert(store_results(&results, &config) == 503);
    return 0;
}
