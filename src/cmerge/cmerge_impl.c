/*
    cmerge_impl.c
    C-Accelerated Merge for Stata - Full Native Implementation

    High-performance drop-in replacement for Stata's merge command.
    Performs the entire merge operation in C without falling back to Stata.

    Architecture (Two-Phase Plugin Call):
    Phase 1 - "load_using":
        - Called when using dataset is current in Stata
        - Loads using data into C memory
        - Sorts using data on key variables (if needed)
        - Stores in static buffer for phase 2

    Phase 2 - "execute":
        - Called when master dataset is current in Stata
        - Loads master data into C memory
        - Sorts master data on key variables (if needed)
        - Performs sorted merge join in C
        - Stores merged result back to Stata

    Merge Types:
    - 1:1: Each key value appears at most once in each dataset
    - m:1: Key may appear multiple times in master, once in using
    - 1:m: Key appears once in master, may appear multiple times in using
    - m:m: Key may appear multiple times in both (not recommended)

    Thread Safety:
    - Parallel data loading per variable
    - Parallel merge output construction by variable
    - Single-threaded merge join (inherently sequential)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "cmerge_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CMERGE_MAX_KEYVARS 32
#define CMERGE_MAX_VARNAME 128

/* Merge type codes */
typedef enum {
    MERGE_1_1 = 0,
    MERGE_M_1 = 1,
    MERGE_1_M = 2,
    MERGE_M_M = 3
} cmerge_type_t;

/* Merge result codes (matches Stata's _merge values) */
typedef enum {
    MERGE_RESULT_MASTER_ONLY = 1,   /* master only */
    MERGE_RESULT_USING_ONLY = 2,    /* using only */
    MERGE_RESULT_BOTH = 3           /* matched */
} cmerge_result_t;

/* Timing uses ctools_timer.h */
#define cmerge_get_time_ms() ctools_timer_ms()

/* ============================================================================
 * Static buffer for using dataset (persists between plugin calls)
 * ============================================================================ */

static stata_data g_using_data;
static int g_using_loaded = 0;
static int g_using_nkeys = 0;
static int g_using_key_indices[CMERGE_MAX_KEYVARS];
static char g_using_varnames[256][CMERGE_MAX_VARNAME];  /* Variable names */
static size_t g_using_nvars_saved = 0;

/* ============================================================================
 * Key comparison functions
 * ============================================================================ */

/*
    Compare two rows by their key variables.
    Returns: -1 if row_a < row_b, 0 if equal, 1 if row_a > row_b
*/
static int cmerge_compare_keys(stata_data *data_a, size_t row_a,
                                stata_data *data_b, size_t row_b,
                                int *key_indices_a, int *key_indices_b,
                                int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        int idx_a = key_indices_a[k];
        int idx_b = key_indices_b[k];
        stata_variable *var_a = &data_a->vars[idx_a];
        stata_variable *var_b = &data_b->vars[idx_b];

        if (var_a->type == STATA_TYPE_DOUBLE) {
            double val_a = var_a->data.dbl[row_a];
            double val_b = var_b->data.dbl[row_b];

            /* Handle missing values - they sort high */
            int miss_a = SF_is_missing(val_a);
            int miss_b = SF_is_missing(val_b);

            if (miss_a && miss_b) continue;  /* Both missing = equal */
            if (miss_a) return 1;            /* Missing sorts high */
            if (miss_b) return -1;

            if (val_a < val_b) return -1;
            if (val_a > val_b) return 1;
        } else {
            /* String comparison */
            char *str_a = var_a->data.str[row_a];
            char *str_b = var_b->data.str[row_b];

            if (str_a == NULL) str_a = "";
            if (str_b == NULL) str_b = "";

            int cmp = strcmp(str_a, str_b);
            if (cmp != 0) return (cmp < 0) ? -1 : 1;
        }
    }
    return 0;  /* All keys equal */
}

/*
    Compare a row's keys within single dataset
*/
static int cmerge_compare_keys_internal(stata_data *data,
                                         size_t row_a, size_t row_b,
                                         int *key_indices, int nkeys)
{
    return cmerge_compare_keys(data, row_a, data, row_b,
                               key_indices, key_indices, nkeys);
}

/* ============================================================================
 * Merge output row tracking
 * ============================================================================ */

typedef struct {
    size_t master_row;     /* Index in master, or SIZE_MAX if none */
    size_t using_row;      /* Index in using, or SIZE_MAX if none */
    int merge_result;      /* _merge value */
} cmerge_output_row_t;

/* ============================================================================
 * Sorted merge join algorithm
 * ============================================================================ */

/*
    Perform sorted merge join.
    Both datasets must be sorted on their key variables.

    Returns: Number of output rows, or -1 on error
*/
static int64_t cmerge_sorted_merge_join(
    stata_data *master, int *master_key_indices,
    stata_data *using_data, int *using_key_indices,
    int nkeys, cmerge_type_t merge_type,
    cmerge_output_row_t **output_rows)
{
    size_t m_nobs = master->nobs;
    size_t u_nobs = using_data->nobs;

    /* Estimate output size - for m:m could be larger */
    size_t est_output = m_nobs + u_nobs;
    size_t output_capacity = est_output;
    size_t output_count = 0;

    cmerge_output_row_t *out = malloc(output_capacity * sizeof(cmerge_output_row_t));
    if (!out) return -1;

    size_t m_idx = 0;
    size_t u_idx = 0;

    while (m_idx < m_nobs || u_idx < u_nobs) {
        /* Grow output array if needed */
        if (output_count >= output_capacity) {
            output_capacity *= 2;
            cmerge_output_row_t *new_out = realloc(out, output_capacity * sizeof(cmerge_output_row_t));
            if (!new_out) {
                free(out);
                return -1;
            }
            out = new_out;
        }

        if (m_idx >= m_nobs) {
            /* Master exhausted, output remaining using rows */
            out[output_count].master_row = SIZE_MAX;
            out[output_count].using_row = u_idx;
            out[output_count].merge_result = MERGE_RESULT_USING_ONLY;
            output_count++;
            u_idx++;
        }
        else if (u_idx >= u_nobs) {
            /* Using exhausted, output remaining master rows */
            out[output_count].master_row = m_idx;
            out[output_count].using_row = SIZE_MAX;
            out[output_count].merge_result = MERGE_RESULT_MASTER_ONLY;
            output_count++;
            m_idx++;
        }
        else {
            /* Compare keys */
            int cmp = cmerge_compare_keys(master, m_idx,
                                          using_data, u_idx,
                                          master_key_indices,
                                          using_key_indices,
                                          nkeys);

            if (cmp < 0) {
                /* Master key < Using key: master only */
                out[output_count].master_row = m_idx;
                out[output_count].using_row = SIZE_MAX;
                out[output_count].merge_result = MERGE_RESULT_MASTER_ONLY;
                output_count++;
                m_idx++;
            }
            else if (cmp > 0) {
                /* Master key > Using key: using only */
                out[output_count].master_row = SIZE_MAX;
                out[output_count].using_row = u_idx;
                out[output_count].merge_result = MERGE_RESULT_USING_ONLY;
                output_count++;
                u_idx++;
            }
            else {
                /* Keys match - handle based on merge type */
                size_t m_start = m_idx;
                size_t u_start = u_idx;

                /* Count matching rows in master */
                size_t m_count = 1;
                while (m_idx + m_count < m_nobs &&
                       cmerge_compare_keys_internal(master, m_start,
                                                    m_idx + m_count,
                                                    master_key_indices,
                                                    nkeys) == 0) {
                    m_count++;
                }

                /* Count matching rows in using */
                size_t u_count = 1;
                while (u_idx + u_count < u_nobs &&
                       cmerge_compare_keys_internal(using_data, u_start,
                                                    u_idx + u_count,
                                                    using_key_indices,
                                                    nkeys) == 0) {
                    u_count++;
                }

                /* Generate output based on merge type */
                switch (merge_type) {
                    case MERGE_1_1:
                        /* 1:1 - should have exactly one match each */
                        out[output_count].master_row = m_start;
                        out[output_count].using_row = u_start;
                        out[output_count].merge_result = MERGE_RESULT_BOTH;
                        output_count++;
                        break;

                    case MERGE_M_1:
                        /* m:1 - multiple master rows map to one using row */
                        for (size_t i = 0; i < m_count; i++) {
                            if (output_count >= output_capacity) {
                                output_capacity *= 2;
                                cmerge_output_row_t *new_out = realloc(out, output_capacity * sizeof(cmerge_output_row_t));
                                if (!new_out) { free(out); return -1; }
                                out = new_out;
                            }
                            out[output_count].master_row = m_start + i;
                            out[output_count].using_row = u_start;
                            out[output_count].merge_result = MERGE_RESULT_BOTH;
                            output_count++;
                        }
                        break;

                    case MERGE_1_M:
                        /* 1:m - one master row maps to multiple using rows */
                        for (size_t i = 0; i < u_count; i++) {
                            if (output_count >= output_capacity) {
                                output_capacity *= 2;
                                cmerge_output_row_t *new_out = realloc(out, output_capacity * sizeof(cmerge_output_row_t));
                                if (!new_out) { free(out); return -1; }
                                out = new_out;
                            }
                            out[output_count].master_row = m_start;
                            out[output_count].using_row = u_start + i;
                            out[output_count].merge_result = MERGE_RESULT_BOTH;
                            output_count++;
                        }
                        break;

                    case MERGE_M_M:
                        /* m:m - full Cartesian product of matching rows */
                        for (size_t i = 0; i < m_count; i++) {
                            for (size_t j = 0; j < u_count; j++) {
                                if (output_count >= output_capacity) {
                                    output_capacity *= 2;
                                    cmerge_output_row_t *new_out = realloc(out, output_capacity * sizeof(cmerge_output_row_t));
                                    if (!new_out) { free(out); return -1; }
                                    out = new_out;
                                }
                                out[output_count].master_row = m_start + i;
                                out[output_count].using_row = u_start + j;
                                out[output_count].merge_result = MERGE_RESULT_BOTH;
                                output_count++;
                            }
                        }
                        break;
                }

                m_idx = m_start + m_count;
                u_idx = u_start + u_count;
            }
        }
    }

    *output_rows = out;
    return (int64_t)output_count;
}

/* ============================================================================
 * Parallel output variable construction
 * ============================================================================ */

typedef struct {
    cmerge_output_row_t *output_rows;
    size_t output_nobs;
    stata_data *source;       /* Master or using data */
    int source_var_idx;       /* Variable index in source */
    int is_master;            /* 1 if from master, 0 if from using */
    stata_variable *out_var;  /* Output variable to populate */
} construct_var_args_t;

static void *construct_var_thread(void *arg)
{
    construct_var_args_t *args = (construct_var_args_t *)arg;
    cmerge_output_row_t *rows = args->output_rows;
    size_t nobs = args->output_nobs;
    stata_variable *src_var = &args->source->vars[args->source_var_idx];
    stata_variable *out = args->out_var;

    out->type = src_var->type;
    out->nobs = nobs;
    out->str_maxlen = src_var->str_maxlen;

    if (src_var->type == STATA_TYPE_DOUBLE) {
        out->data.dbl = malloc(nobs * sizeof(double));
        if (!out->data.dbl) return NULL;

        double missing = SV_missval;
        double *src_data = src_var->data.dbl;
        double *dst_data = out->data.dbl;
        int is_master = args->is_master;

        /* 4x loop unrolling for better performance */
        size_t i = 0;
        size_t nobs_4 = nobs - (nobs % 4);

        for (; i < nobs_4; i += 4) {
            size_t r0 = is_master ? rows[i].master_row : rows[i].using_row;
            size_t r1 = is_master ? rows[i+1].master_row : rows[i+1].using_row;
            size_t r2 = is_master ? rows[i+2].master_row : rows[i+2].using_row;
            size_t r3 = is_master ? rows[i+3].master_row : rows[i+3].using_row;

            dst_data[i]   = (r0 == SIZE_MAX) ? missing : src_data[r0];
            dst_data[i+1] = (r1 == SIZE_MAX) ? missing : src_data[r1];
            dst_data[i+2] = (r2 == SIZE_MAX) ? missing : src_data[r2];
            dst_data[i+3] = (r3 == SIZE_MAX) ? missing : src_data[r3];
        }

        /* Handle remaining elements */
        for (; i < nobs; i++) {
            size_t src_row = is_master ? rows[i].master_row : rows[i].using_row;
            dst_data[i] = (src_row == SIZE_MAX) ? missing : src_data[src_row];
        }
    } else {
        out->data.str = malloc(nobs * sizeof(char *));
        if (!out->data.str) return NULL;

        char **src_data = src_var->data.str;
        char **dst_data = out->data.str;
        int is_master = args->is_master;

        for (size_t i = 0; i < nobs; i++) {
            size_t src_row = is_master ? rows[i].master_row : rows[i].using_row;

            if (src_row == SIZE_MAX) {
                dst_data[i] = strdup("");
            } else {
                char *s = src_data[src_row];
                dst_data[i] = strdup(s ? s : "");
            }
        }
    }

    return arg;  /* Success */
}

/* ============================================================================
 * Phase 1: Load using dataset
 * ============================================================================ */

static ST_retcode cmerge_load_using(const char *args)
{
    char msg[256];
    ST_retcode rc;
    double t_start, t_end;

    /* Parse args: "nkeys key1 key2 ... [sorted] [verbose]" */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    int nkeys = 0;
    int sorted = 0;
    int verbose = 0;
    int arg_idx = 0;

    char *token = strtok(args_copy, " ");
    while (token != NULL) {
        if (arg_idx == 0) {
            nkeys = atoi(token);
            if (nkeys > CMERGE_MAX_KEYVARS) nkeys = CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx <= nkeys) {
            /* Key variable index (1-based from Stata, convert to 0-based) */
            g_using_key_indices[arg_idx - 1] = atoi(token) - 1;
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        /* Ignore other options (mergetype, keep_excludes_using) - no longer used */
        arg_idx++;
        token = strtok(NULL, " ");
    }

    g_using_nkeys = nkeys;

    /* Free any previously loaded using data */
    if (g_using_loaded) {
        stata_data_free(&g_using_data);
        g_using_loaded = 0;
    }

    /* Load using data from Stata */
    if (verbose) {
        SF_display("cmerge: Loading using dataset into C...\n");
    }
    t_start = cmerge_get_time_ms();

    size_t nvars = SF_nvars();
    g_using_nvars_saved = nvars;

    /* Save variable names for later use */
    for (size_t v = 0; v < nvars && v < 256; v++) {
        /* Get variable name via macro - variable names are 1-indexed */
        /* We'll use a simple approach: store the index, .ado passes names */
        g_using_varnames[v][0] = '\0';
    }

    stata_data_init(&g_using_data);
    rc = ctools_data_load(&g_using_data, nvars);

    if (rc != STATA_OK) {
        SF_error("cmerge: Failed to load using dataset\n");
        return 459;
    }

    t_end = cmerge_get_time_ms();

    if (verbose) {
        snprintf(msg, sizeof(msg), "  Loaded %zu obs, %zu vars in %.1f ms\n",
                 g_using_data.nobs, g_using_data.nvars, t_end - t_start);
        SF_display(msg);
    }

    /* Sort if needed for sorted merge join */
    if (!sorted && g_using_data.nobs > 1) {
        if (verbose) {
            SF_display("cmerge: Sorting using dataset...\n");
        }
        t_start = cmerge_get_time_ms();

        /* Convert 0-based to 1-based for sort function */
        int *sort_vars = malloc(nkeys * sizeof(int));
        for (int i = 0; i < nkeys; i++) {
            sort_vars[i] = g_using_key_indices[i] + 1;
        }

        rc = ctools_sort_radix_lsd(&g_using_data, sort_vars, nkeys);
        free(sort_vars);

        if (rc != STATA_OK) {
            SF_error("cmerge: Failed to sort using dataset\n");
            stata_data_free(&g_using_data);
            return 459;
        }

        t_end = cmerge_get_time_ms();
        if (verbose) {
            snprintf(msg, sizeof(msg), "  Sorted in %.1f ms\n", t_end - t_start);
            SF_display(msg);
        }
    }

    g_using_loaded = 1;

    /* Return using nobs and nvars via scalars */
    SF_scal_save("_cmerge_using_nobs", (double)g_using_data.nobs);
    SF_scal_save("_cmerge_using_nvars", (double)g_using_data.nvars);

    return 0;
}

/* ============================================================================
 * Phase 2: Execute merge
 * ============================================================================ */

static ST_retcode cmerge_execute(const char *args)
{
    char msg[512];
    ST_retcode rc;
    double t_start, t_end, t_total_start;

    t_total_start = cmerge_get_time_ms();

    /* Check that using data is loaded */
    if (!g_using_loaded) {
        SF_error("cmerge: Using data not loaded. Call load_using first.\n");
        return 459;
    }

    /* Parse args: "merge_type nkeys key1 key2 ... [sorted] [verbose] [nogenerate] [master_nobs N]" */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    cmerge_type_t merge_type = MERGE_1_1;
    int nkeys = 0;
    int master_key_indices[CMERGE_MAX_KEYVARS];
    int sorted = 0;
    int verbose = 0;
    int nogenerate = 0;
    size_t actual_master_nobs = 0;  /* If set, use this instead of SF_nobs() for master */
    int arg_idx = 0;

    char *token = strtok(args_copy, " ");
    while (token != NULL) {
        if (arg_idx == 0) {
            merge_type = (cmerge_type_t)atoi(token);
        }
        else if (arg_idx == 1) {
            nkeys = atoi(token);
            if (nkeys > CMERGE_MAX_KEYVARS) nkeys = CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx >= 2 && arg_idx < 2 + nkeys) {
            master_key_indices[arg_idx - 2] = atoi(token) - 1;  /* 0-based */
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        else if (strcmp(token, "nogenerate") == 0) {
            nogenerate = 1;
        }
        else if (strcmp(token, "master_nobs") == 0) {
            /* Next token is the actual master nobs */
            token = strtok(NULL, " ");
            if (token) {
                actual_master_nobs = (size_t)atol(token);
            }
        }
        /* Ignore other options (keep_excludes_using) - no longer used */
        arg_idx++;
        token = strtok(NULL, " ");
    }

    /* Load master data */
    if (verbose) {
        SF_display("cmerge: Loading master dataset...\n");
    }
    t_start = cmerge_get_time_ms();

    stata_data master;
    stata_data_init(&master);
    size_t master_nvars = SF_nvars();

    rc = ctools_data_load(&master, master_nvars);
    if (rc != STATA_OK) {
        SF_error("cmerge: Failed to load master dataset\n");
        return 459;
    }

    /* If actual_master_nobs is set, we've expanded the dataset for output
       and need to truncate to the actual master data */
    if (actual_master_nobs > 0 && actual_master_nobs < master.nobs) {
        master.nobs = actual_master_nobs;
        /* Also need to truncate each variable's data */
        for (size_t v = 0; v < master.nvars; v++) {
            master.vars[v].nobs = actual_master_nobs;
        }
    }

    t_end = cmerge_get_time_ms();
    double time_load_master = t_end - t_start;

    if (verbose) {
        snprintf(msg, sizeof(msg), "  Loaded %zu obs, %zu vars in %.1f ms\n",
                 master.nobs, master.nvars, time_load_master);
        SF_display(msg);
    }

    double time_sort_master = 0;
    double time_merge = 0;
    cmerge_output_row_t *output_rows = NULL;
    int64_t output_nobs = 0;
    size_t n_master_only = 0, n_using_only = 0, n_matched = 0;

    /* Sort master if needed for sorted merge join */
    if (!sorted && master.nobs > 1) {
        if (verbose) {
            SF_display("cmerge: Sorting master dataset...\n");
        }
        t_start = cmerge_get_time_ms();

        int *sort_vars = malloc(nkeys * sizeof(int));
        for (int i = 0; i < nkeys; i++) {
            sort_vars[i] = master_key_indices[i] + 1;
        }

        rc = ctools_sort_radix_lsd(&master, sort_vars, nkeys);
        free(sort_vars);

        if (rc != STATA_OK) {
            SF_error("cmerge: Failed to sort master dataset\n");
            stata_data_free(&master);
            return 459;
        }

        t_end = cmerge_get_time_ms();
        time_sort_master = t_end - t_start;

        if (verbose) {
            snprintf(msg, sizeof(msg), "  Sorted in %.1f ms\n", time_sort_master);
            SF_display(msg);
        }
    }

    /* Perform sorted merge join */
    if (verbose) {
        SF_display("cmerge: Performing sorted merge join...\n");
    }
    t_start = cmerge_get_time_ms();

    output_nobs = cmerge_sorted_merge_join(
        &master, master_key_indices,
        &g_using_data, g_using_key_indices,
        nkeys, merge_type,
        &output_rows);

    if (output_nobs < 0) {
        SF_error("cmerge: Merge join failed (out of memory)\n");
        stata_data_free(&master);
        return 920;
    }

    t_end = cmerge_get_time_ms();
    time_merge = t_end - t_start;

    if (verbose) {
        snprintf(msg, sizeof(msg), "  Merge produced %lld rows in %.1f ms\n",
                 (long long)output_nobs, time_merge);
        SF_display(msg);
    }

    /* Count merge results for reporting */
    for (int64_t i = 0; i < output_nobs; i++) {
        switch (output_rows[i].merge_result) {
            case MERGE_RESULT_MASTER_ONLY: n_master_only++; break;
            case MERGE_RESULT_USING_ONLY: n_using_only++; break;
            case MERGE_RESULT_BOTH: n_matched++; break;
        }
    }

    /* Build output dataset */
    if (verbose) {
        SF_display("cmerge: Building output dataset...\n");
    }
    t_start = cmerge_get_time_ms();

    /* Determine which using variables to add (non-key) */
    int *using_nonkey_vars = malloc(g_using_data.nvars * sizeof(int));
    size_t n_using_nonkey = 0;

    for (size_t v = 0; v < g_using_data.nvars; v++) {
        int is_key = 0;
        for (int k = 0; k < g_using_nkeys; k++) {
            if ((size_t)g_using_key_indices[k] == v) {
                is_key = 1;
                break;
            }
        }
        if (!is_key) {
            using_nonkey_vars[n_using_nonkey++] = (int)v;
        }
    }

    /* Total output variables = master vars + using non-key vars + _merge */
    size_t output_nvars = master_nvars + n_using_nonkey;
    if (!nogenerate) {
        output_nvars++;  /* For _merge */
    }

    /* Allocate output data structure */
    stata_data output;
    stata_data_init(&output);
    output.nobs = (size_t)output_nobs;
    output.nvars = output_nvars;
    output.vars = calloc(output_nvars, sizeof(stata_variable));

    if (!output.vars) {
        free(output_rows);
        free(using_nonkey_vars);
        stata_data_free(&master);
        return 920;
    }

    /* Parallel construction of output variables */
    size_t n_threads = master_nvars + n_using_nonkey;
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    construct_var_args_t *thread_args = malloc(n_threads * sizeof(construct_var_args_t));

    if (!threads || !thread_args) {
        free(threads);
        free(thread_args);
        free(output_rows);
        free(using_nonkey_vars);
        stata_data_free(&master);
        stata_data_free(&output);
        return 920;
    }

    size_t thread_idx = 0;

    /* Master variables */
    for (size_t v = 0; v < master_nvars; v++) {
        thread_args[thread_idx].output_rows = output_rows;
        thread_args[thread_idx].output_nobs = (size_t)output_nobs;
        thread_args[thread_idx].source = &master;
        thread_args[thread_idx].source_var_idx = (int)v;
        thread_args[thread_idx].is_master = 1;
        thread_args[thread_idx].out_var = &output.vars[thread_idx];

        pthread_create(&threads[thread_idx], NULL, construct_var_thread, &thread_args[thread_idx]);
        thread_idx++;
    }

    /* Using non-key variables */
    for (size_t i = 0; i < n_using_nonkey; i++) {
        int v = using_nonkey_vars[i];
        thread_args[thread_idx].output_rows = output_rows;
        thread_args[thread_idx].output_nobs = (size_t)output_nobs;
        thread_args[thread_idx].source = &g_using_data;
        thread_args[thread_idx].source_var_idx = v;
        thread_args[thread_idx].is_master = 0;
        thread_args[thread_idx].out_var = &output.vars[thread_idx];

        pthread_create(&threads[thread_idx], NULL, construct_var_thread, &thread_args[thread_idx]);
        thread_idx++;
    }

    /* Wait for all threads */
    for (size_t t = 0; t < n_threads; t++) {
        pthread_join(threads[t], NULL);
    }

    free(threads);
    free(thread_args);

    /* Construct _merge variable if requested */
    if (!nogenerate) {
        stata_variable *merge_var = &output.vars[output_nvars - 1];
        merge_var->type = STATA_TYPE_DOUBLE;
        merge_var->nobs = (size_t)output_nobs;
        merge_var->data.dbl = malloc((size_t)output_nobs * sizeof(double));

        if (merge_var->data.dbl) {
            for (int64_t i = 0; i < output_nobs; i++) {
                merge_var->data.dbl[i] = (double)output_rows[i].merge_result;
            }
        }
    }

    t_end = cmerge_get_time_ms();
    double time_build = t_end - t_start;

    if (verbose) {
        snprintf(msg, sizeof(msg), "  Built output in %.1f ms\n", time_build);
        SF_display(msg);
    }

    /* Store output back to Stata
       NOTE: The .ado has already expanded the dataset to accommodate
       the maximum possible output size (master + using). We write all
       output rows, and the .ado will trim excess afterward. */
    if (verbose) {
        SF_display("cmerge: Storing results to Stata...\n");
    }
    t_start = cmerge_get_time_ms();

    rc = ctools_data_store(&output, 1);

    t_end = cmerge_get_time_ms();
    double time_store = t_end - t_start;

    if (verbose) {
        snprintf(msg, sizeof(msg), "  Stored %zu obs in %.1f ms\n",
                 output.nobs, time_store);
        SF_display(msg);
    }

    /* Save merge counts to scalars */
    SF_scal_save("_cmerge_N_1", (double)n_master_only);
    SF_scal_save("_cmerge_N_2", (double)n_using_only);
    SF_scal_save("_cmerge_N_3", (double)n_matched);
    SF_scal_save("_cmerge_N", (double)output_nobs);

    double total_time = cmerge_get_time_ms() - t_total_start;
    SF_scal_save("_cmerge_time", total_time / 1000.0);

    if (verbose) {
        snprintf(msg, sizeof(msg), "\nTotal C time: %.1f ms\n", total_time);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "  Load master: %.1f ms, Sort master: %.1f ms\n",
                 time_load_master, time_sort_master);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "  Merge: %.1f ms, Build: %.1f ms, Store: %.1f ms\n",
                 time_merge, time_build, time_store);
        SF_display(msg);
    }

    /* Cleanup */
    free(output_rows);
    free(using_nonkey_vars);
    stata_data_free(&master);
    stata_data_free(&output);

    /* Free using data after merge */
    stata_data_free(&g_using_data);
    g_using_loaded = 0;

    return rc;
}

/* ============================================================================
 * Clear cached using data (for cleanup/error handling)
 * ============================================================================ */

static ST_retcode cmerge_clear(void)
{
    if (g_using_loaded) {
        stata_data_free(&g_using_data);
        g_using_loaded = 0;
    }
    return 0;
}

/* ============================================================================
 * Plugin entry point
 * ============================================================================ */

ST_retcode cmerge_main(const char *args)
{
    if (args == NULL || strlen(args) == 0) {
        SF_error("cmerge: no arguments specified\n");
        SF_error("Usage: plugin call ctools_plugin, \"cmerge <phase> <args>\"\n");
        SF_error("  Phases: load_using, execute, clear\n");
        return 198;
    }

    /* Parse phase from first argument */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *phase = strtok(args_copy, " ");
    const char *rest = args + strlen(phase);
    while (*rest == ' ') rest++;  /* Skip spaces */

    if (strcmp(phase, "load_using") == 0) {
        return cmerge_load_using(rest);
    }
    else if (strcmp(phase, "execute") == 0) {
        return cmerge_execute(rest);
    }
    else if (strcmp(phase, "clear") == 0) {
        return cmerge_clear();
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "cmerge: unknown phase '%s'\n", phase);
        SF_error(msg);
        SF_error("Valid phases: load_using, execute, clear\n");
        return 198;
    }
}
