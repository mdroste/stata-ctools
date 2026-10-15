#ifndef CTOOLS_EXCEL_ADDRESS_H
#define CTOOLS_EXCEL_ADDRESS_H
#include <stdbool.h>
/* Outputs are one based. If end is NULL, require end-of-string; otherwise stop
 * at the first suffix (range adapters consume ':' themselves). On failure no
 * output is changed. Dollar signs are allowed only when absolute is true. */
bool ctools_excel_address(const char *text, int max_col, int max_row, bool absolute, int *col,
                          int *row, const char **end);
#endif
