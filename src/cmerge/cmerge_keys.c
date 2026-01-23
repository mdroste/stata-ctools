/*
 * cmerge_keys.c
 * Key comparison functions for cmerge
 */

#include <string.h>
#include "cmerge_keys.h"
#include "stplugin.h"

/* Check if all keys are numeric (for fast path selection) */
int cmerge_all_keys_numeric(stata_data *data, int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        if (data->vars[k].type != STATA_TYPE_DOUBLE) {
            return 0;
        }
    }
    return 1;
}

/* Fast path: compare all-numeric keys without type checking in loop */
int cmerge_compare_keys_numeric(stata_data *data_a, size_t row_a,
                                 stata_data *data_b, size_t row_b,
                                 int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        double val_a = data_a->vars[k].data.dbl[row_a];
        double val_b = data_b->vars[k].data.dbl[row_b];

        /* Missing value handling - must distinguish extended missing (.a, .b, etc.) */
        int miss_a = SF_is_missing(val_a);
        int miss_b = SF_is_missing(val_b);

        if (miss_a | miss_b) {
            if (miss_a && miss_b) {
                /* Both missing - compare actual values to distinguish .a, .b, .z, etc. */
                if (val_a < val_b) return -1;
                if (val_a > val_b) return 1;
                continue;  /* Same missing type */
            }
            return miss_a ? 1 : -1;
        }

        /* Use subtraction for branchless comparison when possible */
        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
    }
    return 0;
}

/* Fast path: compare same-dataset numeric keys */
int cmerge_compare_keys_numeric_same(stata_data *data, size_t row_a,
                                      size_t row_b, int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        double *col = data->vars[k].data.dbl;
        double val_a = col[row_a];
        double val_b = col[row_b];

        int miss_a = SF_is_missing(val_a);
        int miss_b = SF_is_missing(val_b);

        if (miss_a | miss_b) {
            if (miss_a && miss_b) {
                /* Both missing - compare actual values to distinguish .a, .b, .z, etc. */
                if (val_a < val_b) return -1;
                if (val_a > val_b) return 1;
                continue;  /* Same missing type */
            }
            return miss_a ? 1 : -1;
        }

        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
    }
    return 0;
}

/* General path: handles mixed string/numeric keys */
int cmerge_compare_keys(stata_data *data_a, size_t row_a,
                         stata_data *data_b, size_t row_b,
                         int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        stata_variable *var_a = &data_a->vars[k];
        stata_variable *var_b = &data_b->vars[k];

        if (var_a->type == STATA_TYPE_DOUBLE) {
            double val_a = var_a->data.dbl[row_a];
            double val_b = var_b->data.dbl[row_b];

            int miss_a = SF_is_missing(val_a);
            int miss_b = SF_is_missing(val_b);

            if (miss_a && miss_b) {
                /* Both missing - compare actual values to distinguish .a, .b, .z, etc. */
                if (val_a < val_b) return -1;
                if (val_a > val_b) return 1;
                continue;  /* Same missing type */
            }
            if (miss_a) return 1;
            if (miss_b) return -1;

            if (val_a < val_b) return -1;
            if (val_a > val_b) return 1;
        } else {
            char *str_a = var_a->data.str[row_a];
            char *str_b = var_b->data.str[row_b];

            if (str_a == NULL) str_a = "";
            if (str_b == NULL) str_b = "";

            int cmp = strcmp(str_a, str_b);
            if (cmp != 0) return (cmp < 0) ? -1 : 1;
        }
    }
    return 0;
}

/* General path: same-dataset comparison */
int cmerge_compare_keys_same(stata_data *data, size_t row_a, size_t row_b, int nkeys)
{
    return cmerge_compare_keys(data, row_a, data, row_b, nkeys);
}
