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
                else if (update_mode || replace_mode) {
                    /* Shared var: update fills missing master values; replace
                     * also overwrites with nonmissing using values. */
                    double current_val;
                    SF_vdata(dest_idx, (ST_int)(i + 1), &current_val);
                    double using_val = (using_row >= 0 && using_row < (int32_t)src->nobs)
                                       ? src->data.dbl[using_row] : SV_missval;
                    int master_missing = SF_is_missing(current_val);
                    int using_missing = SF_is_missing(using_val);
                    should_write = replace_mode ? (!using_missing || master_missing)
                                                : master_missing;
                    /* merge counts a missing master value as updated whenever
                     * the using value differs, extended missing included
                     * (. <- .a and .a <- . are code 4) */
                    uint8_t flag = 0;
                    if (master_missing) {
                        if (using_val != current_val) flag = CMERGE_ROW_UPDATED;
                    } else if (!using_missing && current_val != using_val) {
                        flag = CMERGE_ROW_CONFLICT;
                    }
                    if (flag && a->row_flags) {
                        __atomic_fetch_or(&a->row_flags[i], flag, __ATOMIC_RELAXED);
                    }
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
                else if (update_mode || replace_mode) {
                    /* Shared var: "" is missing; update fills empty master
                     * values; replace also overwrites with nonempty using values. */
                    SF_sdata(dest_idx, (ST_int)(i + 1), str_buf);
                    int master_missing = (str_buf[0] == '\0');
                    const char *using_val = (using_row >= 0 && using_row < (int32_t)src->nobs &&
                                             src->data.str[using_row] != NULL)
                                            ? src->data.str[using_row] : "";
                    int using_missing = (using_val[0] == '\0');
                    should_write = replace_mode ? (!using_missing || master_missing)
                                                : master_missing;
                    uint8_t flag = 0;
                    if (master_missing && !using_missing) flag = CMERGE_ROW_UPDATED;
                    else if (!master_missing && !using_missing &&
                             strcmp(str_buf, using_val) != 0) flag = CMERGE_ROW_CONFLICT;
                    if (flag && a->row_flags) {
                        __atomic_fetch_or(&a->row_flags[i], flag, __ATOMIC_RELAXED);
                    }
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
