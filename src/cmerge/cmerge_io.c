/*
 * cmerge_io.c
 * Data streaming and I/O functions for cmerge
 */

#include <string.h>
#include "cmerge_io.h"
#include "stplugin.h"

void *cmerge_write_keepusing_var_thread(void *arg)
{
    cmerge_keepusing_write_args_t *a = (cmerge_keepusing_write_args_t *)arg;
    size_t output_nobs = a->output_nobs;
    ST_int dest_idx = a->dest_idx;
    int is_shared = a->is_shared;
    int update_mode = a->update_mode;
    int replace_mode = a->replace_mode;
    cmerge_output_spec_t *specs = a->specs;
    stata_variable *src = &g_using_cache.keepusing.data.vars[a->keepusing_idx];

    a->success = 0;

    if (src->type == STATA_TYPE_DOUBLE) {
        /* Numeric keepusing variable */
        for (size_t i = 0; i < output_nobs; i++) {
            int32_t using_row = specs[i].using_sorted_row;
            int8_t merge_result = specs[i].merge_result;

            /* Determine if we should write using value */
            int should_write = 0;

            if (merge_result == MERGE_RESULT_USING_ONLY) {
                /* Always write for using-only rows */
                should_write = 1;
            }
            else if (merge_result == MERGE_RESULT_BOTH) {
                /* Matched row - depends on shared status and update/replace */
                if (!is_shared) {
                    /* Non-shared var: always write using value for matched rows */
                    should_write = 1;
                }
                else if (replace_mode) {
                    /* Shared var with replace: overwrite only if using value is non-missing */
                    should_write = (using_row >= 0 && using_row < (int32_t)src->nobs &&
                                    !SF_is_missing(src->data.dbl[using_row]));
                    /* If using is missing, still fill missing master values */
                    if (!should_write) {
                        double current_val;
                        SF_vdata(dest_idx, (ST_int)(i + 1), &current_val);
                        should_write = SF_is_missing(current_val);
                    }
                }
                else if (update_mode) {
                    /* Shared var with update: only if master value is missing */
                    double current_val;
                    SF_vdata(dest_idx, (ST_int)(i + 1), &current_val);
                    should_write = SF_is_missing(current_val);
                }
                /* else: shared var without update/replace - don't overwrite master */
            }

            if (should_write && using_row >= 0 && using_row < (int32_t)src->nobs) {
                SF_vstore(dest_idx, (ST_int)(i + 1), src->data.dbl[using_row]);
            }
        }
    } else {
        /* String keepusing variable */
        char str_buf[2049];
        for (size_t i = 0; i < output_nobs; i++) {
            int32_t using_row = specs[i].using_sorted_row;
            int8_t merge_result = specs[i].merge_result;

            /* Determine if we should write using value */
            int should_write = 0;

            if (merge_result == MERGE_RESULT_USING_ONLY) {
                /* Always write for using-only rows */
                should_write = 1;
            }
            else if (merge_result == MERGE_RESULT_BOTH) {
                /* Matched row - depends on shared status and update/replace */
                if (!is_shared) {
                    /* Non-shared var: always write using value for matched rows */
                    should_write = 1;
                }
                else if (replace_mode) {
                    /* Shared var with replace: overwrite only if using value is non-empty */
                    int using_nonempty = (using_row >= 0 && using_row < (int32_t)src->nobs &&
                                          src->data.str[using_row] != NULL &&
                                          src->data.str[using_row][0] != '\0');
                    if (using_nonempty) {
                        should_write = 1;
                    }
                    else {
                        /* If using is empty, still fill empty master values */
                        SF_sdata(dest_idx, (ST_int)(i + 1), str_buf);
                        should_write = (str_buf[0] == '\0');
                    }
                }
                else if (update_mode) {
                    /* Shared var with update: only if master value is missing (empty string) */
                    SF_sdata(dest_idx, (ST_int)(i + 1), str_buf);
                    should_write = (str_buf[0] == '\0');
                }
                /* else: shared var without update/replace - don't overwrite master */
            }

            if (should_write && using_row >= 0 && using_row < (int32_t)src->nobs) {
                const char *val = src->data.str[using_row] ? src->data.str[using_row] : "";
                SF_sstore(dest_idx, (ST_int)(i + 1), (char *)val);
            }
        }
    }

    a->success = 1;
    return NULL;
}
