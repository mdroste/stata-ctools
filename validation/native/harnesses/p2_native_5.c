
#include "cimport/cimport_xlsx.c"
int main(void) {
    setup(); nrows = 3; string_var = 1; bad_row = 2;
    char *strings[] = {"a", "b", "c"}; double values[] = {1, 2, 3};
    CImportColumnInfo col = {.type = CIMPORT_COL_STRING};
    CImportColumnCache cache = {.count = 3, .string_data = strings, .numeric_data = values};
    xlsx_store_task task = {.col = &col, .cache = &cache, .var = 1, .stata_nobs = 3};
    assert(xlsx_store_worker(&task) != NULL && task.error == 459);
    string_var = 0; col.type = CIMPORT_COL_NUMERIC;
    assert(xlsx_store_worker(&task) != NULL && task.error == 459);
    task.stata_nobs = 2; writes = 0;
    assert(xlsx_store_worker(&task) != NULL && task.error == 198 && writes == 0);
    return 0;
}
