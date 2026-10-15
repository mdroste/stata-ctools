
#include "cimport/cimport_impl.c"
int main(void) {
    setup(); nrows = 3; bad_row = 2;
    /* One column of three rows ("1", "2", "3"): a store tile must stop at the
     * failing observation and report its code, for string and numeric columns. */
    CImportParsedRow *rows[3];
    for (int i = 0; i < 3; i++) {
        rows[i] = calloc(1, cimport_row_size(1));
        CImportFieldRef field = {.offset = (uint64_t)i, .length = 1};
        assert(cimport_fill_row(rows[i], &field, 1));
        double value = i + 1;
        memcpy(&cimport_row_values(rows[i])[0], &value, sizeof(value));
    }
    CImportParsedChunk chunk = {.rows = rows, .num_rows = 3};
    CImportContext ctx = {.chunks = &chunk, .num_chunks = 1, .total_rows = 3, .num_columns = 1,
                          .file_data = "123", .decimal_separator = '.', .quote_char = '"'};
    unsigned char kind[1] = {CIMPORT_STORE_STRING};
    string_var = 1; writes = 0;
    CImportStoreTile tile = {&ctx, kind, 0, 0, 3, 1, 0};
    assert(cimport_store_tile_worker(&tile) != NULL);
    assert(tile.error == 459 && writes == 2);
    string_var = 0; kind[0] = CIMPORT_STORE_NUMERIC; writes = 0;
    tile = (CImportStoreTile){&ctx, kind, 0, 0, 3, 1, 0};
    assert(cimport_store_tile_worker(&tile) != NULL);
    assert(tile.error == 459 && writes == 2);
    /* strL columns are transferred by the ado: the tile writes nothing */
    kind[0] = CIMPORT_STORE_STRL; writes = 0;
    tile = (CImportStoreTile){&ctx, kind, 0, 0, 3, 1, 0};
    assert(cimport_store_tile_worker(&tile) == NULL && tile.error == 0 && writes == 0);
    for (int i = 0; i < 3; i++) free(rows[i]);
    return 0;
}
