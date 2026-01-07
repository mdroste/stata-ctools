/*
    ctools_types.c
    Implementation of common utility functions for Stata-C data structures
*/

#include <stdlib.h>
#include <string.h>
#include "ctools_types.h"

void stata_data_init(stata_data *data)
{
    if (data == NULL) return;

    data->nobs = 0;
    data->nvars = 0;
    data->vars = NULL;
    data->sort_order = NULL;
}

void stata_data_free(stata_data *data)
{
    size_t i, j;

    if (data == NULL) return;

    /* Free each variable's data */
    if (data->vars != NULL) {
        for (i = 0; i < data->nvars; i++) {
            if (data->vars[i].type == STATA_TYPE_DOUBLE) {
                free(data->vars[i].data.dbl);
            } else if (data->vars[i].type == STATA_TYPE_STRING) {
                if (data->vars[i].data.str != NULL) {
                    for (j = 0; j < data->vars[i].nobs; j++) {
                        free(data->vars[i].data.str[j]);
                    }
                    free(data->vars[i].data.str);
                }
            }
        }
        free(data->vars);
    }

    /* Free sort order array */
    free(data->sort_order);

    /* Reset the structure */
    stata_data_init(data);
}
