
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
#include "ctools_data_io.c"
int main(void) {
    setup();
    int indices[] = {1, 2}, widths[] = {8, 8};
    ctools_filtered_data fd;
    for (int pool = 0; pool <= 1; pool++) {
        ctools_set_max_threads(pool ? 4 : 1);
        for (int columns = 1; columns <= 2; columns++)
        for (int strings = 0; strings <= 1; strings++)
        for (int hints = 0; hints <= 1; hints++)
        for (filter = 0; filter <= 1; filter++) {
            string_var = strings ? 1 : 0;
            bad_row = 3;
            assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, hints ? widths : NULL) == STATA_ERR_STATA_READ);
            assert(fd.data.vars == NULL && fd.obs_map == NULL);
            bad_row = 0;
            assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, hints ? widths : NULL) == 0);
            bad_row = 3;
            assert(ctools_data_store_ex(&fd.data, indices, columns, 1) == STATA_ERR_STATA_WRITE);
            if (!strings) {
                assert(ctools_store_filtered(fd.data.vars[0].data.dbl, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
                assert(ctools_store_filtered_rowpar(fd.data.vars[0].data.dbl, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
            } else {
                assert(ctools_store_filtered_str(fd.data.vars[0].data.str, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
            }
            ctools_filtered_data_free(&fd);
        }
        ctools_destroy_global_pool();
    }
    bad_row = filter = 0; string_var = 1; strl = 1; reads = 0;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == 0 && reads == 10);
    ctools_filtered_data_free(&fd);
    reads = 0; strl_length = 2046;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0);
    strl_length = 4; binary_value = 1;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0);
    binary_value = 0; strl = 0; widths[0] = 1;
    assert(ctools_data_load_ex(&fd, indices, 1, 1, 10, 0, widths) == STATA_ERR_STATA_READ);
    indices[0] = 3;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_INVALID_INPUT);
    indices[0] = 1;
    assert(ctools_data_load(&fd, indices, 1, 10, 1, 0) == STATA_ERR_INVALID_INPUT);
    assert(ctools_data_load(&fd, indices, 1, 1, nrows+1, 0) == STATA_ERR_INVALID_INPUT);
    double values[] = {1, 2}; perm_idx_t map[] = {1, 6001}; writes = 0; string_var = 0;
    assert(ctools_store_filtered_rowpar(values, 2, 1, map) == STATA_ERR_INVALID_INPUT && writes == 0);
    int64_t rows[] = {0, 2, 1}; bad_row = 3; writes = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_STATA_READ && writes == 0);
    string_var = 1; writes = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_STATA_READ && writes == 0);
    strl = 1; reads = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0 && writes == 0);
    ctools_destroy_global_pool();
    return 0;
}
