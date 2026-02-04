/*
 * test_civreghdfe_asan.c
 * Standalone test harness for civreghdfe with AddressSanitizer
 *
 * Build with:
 *   clang -fsanitize=address -g -O1 -DTEST_MOCK_SPI \
 *     -I. -Isrc -Isrc/civreghdfe -Isrc/creghdfe \
 *     src/test_civreghdfe_asan.c \
 *     src/civreghdfe/*.c src/creghdfe/*.c \
 *     src/ctools_*.c \
 *     -lm -lomp -o test_civreghdfe_asan
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Include the mock SPI header instead of real stplugin.h */
#include "test_stplugin_mock.h"

/* Include ctools_types for STATA_OK and other definitions */
#include "ctools_types.h"

/* ============================================================================
 * Mock SPI Implementation
 * ============================================================================ */

/* Global test state */
static ST_int g_in1 = 1;
static ST_int g_in2 = 200;
static ST_int g_nvar = 10;

/* Scalar storage */
#define MAX_SCALARS 200
static struct {
    char name[64];
    ST_double value;
    int used;
} g_scalars[MAX_SCALARS];
static int g_num_scalars = 0;

/* Test data arrays - column-major */
static ST_double *g_test_data = NULL;
static int g_test_nobs = 0;
static int g_test_nvar = 0;

/* Mock SPI functions */
ST_int SF_in1(void) { return g_in1; }
ST_int SF_in2(void) { return g_in2; }
ST_int SF_nvars(void) { return g_nvar; }
ST_int SF_nobs(void) { return g_in2; }
ST_int SF_ifobs(ST_int obs) { (void)obs; return 1; }
ST_int SF_var_is_string(ST_int var) { (void)var; return 0; }
ST_int SF_var_is_strl(ST_int var) { (void)var; return 0; }

ST_int SF_scal_use(const char *name, ST_double *value) {
    for (int i = 0; i < g_num_scalars; i++) {
        if (g_scalars[i].used && strcmp(g_scalars[i].name, name) == 0) {
            *value = g_scalars[i].value;
            return 0;
        }
    }
    *value = 0;
    return 0;
}

ST_int SF_scal_save(const char *name, ST_double value) {
    /* Find existing or empty slot */
    for (int i = 0; i < g_num_scalars; i++) {
        if (g_scalars[i].used && strcmp(g_scalars[i].name, name) == 0) {
            g_scalars[i].value = value;
            return 0;
        }
    }
    if (g_num_scalars < MAX_SCALARS) {
        strncpy(g_scalars[g_num_scalars].name, name, 63);
        g_scalars[g_num_scalars].name[63] = '\0';
        g_scalars[g_num_scalars].value = value;
        g_scalars[g_num_scalars].used = 1;
        g_num_scalars++;
    }
    return 0;
}

ST_int SF_vdata(ST_int var, ST_int obs, ST_double *value) {
    if (g_test_data && var >= 1 && var <= g_test_nvar && obs >= 1 && obs <= g_test_nobs) {
        *value = g_test_data[(var - 1) * (size_t)g_test_nobs + (obs - 1)];
        return 0;
    }
    *value = 8.988e307 + 1;  /* SV_missval */
    return 0;
}

ST_int SF_vstore(ST_int var, ST_int obs, ST_double value) {
    if (g_test_data && var >= 1 && var <= g_test_nvar && obs >= 1 && obs <= g_test_nobs) {
        g_test_data[(var - 1) * (size_t)g_test_nobs + (obs - 1)] = value;
    }
    return 0;
}

ST_int SF_sdata(ST_int var, ST_int obs, char *value) {
    (void)var; (void)obs;
    value[0] = '\0';
    return 0;
}

ST_int SF_sstore(ST_int var, ST_int obs, char *value) {
    (void)var; (void)obs; (void)value;
    return 0;
}

ST_int SF_mat_store(const char *name, ST_int row, ST_int col, ST_double value) {
    (void)name; (void)row; (void)col; (void)value;
    return 0;
}

ST_int SF_mat_el(const char *name, ST_int row, ST_int col, ST_double *value) {
    (void)name; (void)row; (void)col;
    *value = 0;
    return 0;
}

void SF_error(const char *msg) { fprintf(stderr, "ERROR: %s", msg); }
void SF_display(const char *msg) { printf("%s", msg); }

/* ============================================================================
 * Test Setup
 * ============================================================================ */

static void reset_scalars(void) {
    for (int i = 0; i < g_num_scalars; i++) {
        g_scalars[i].used = 0;
    }
    g_num_scalars = 0;
}

static void setup_test_data(int nobs, int seed, int test_variant) {
    srand(seed);

    g_test_nobs = nobs;
    g_test_nvar = 12;  /* y, x_endog, x_exog, z1, z2, id, time, weight, cluster, cluster2, extra1, extra2 */
    g_in1 = 1;
    g_in2 = nobs;
    g_nvar = g_test_nvar;

    if (g_test_data) {
        free(g_test_data);
    }
    g_test_data = (ST_double *)calloc((size_t)g_test_nobs * g_test_nvar, sizeof(ST_double));
    if (!g_test_data) {
        fprintf(stderr, "Failed to allocate test data\n");
        exit(1);
    }

    /* Generate random data */
    int num_ids = 20 + (test_variant % 10);
    for (int i = 0; i < nobs; i++) {
        int id = (i % num_ids) + 1;
        int time = (i / num_ids) + 1;

        double z1 = (rand() / (double)RAND_MAX - 0.5) * 2;
        double z2 = (rand() / (double)RAND_MAX - 0.5) * 2;
        double x_exog = (rand() / (double)RAND_MAX - 0.5) * 2;
        double x_endog = 0.5 * z1 + 0.3 * z2 + (rand() / (double)RAND_MAX - 0.5);
        double y = 1 + 0.5 * x_endog + 0.3 * x_exog + (rand() / (double)RAND_MAX - 0.5);
        double weight = 0.5 + rand() / (double)RAND_MAX;

        /* Store in column-major order (var index - 1) * nobs + (obs index - 1) */
        g_test_data[0 * nobs + i] = y;
        g_test_data[1 * nobs + i] = x_endog;
        g_test_data[2 * nobs + i] = x_exog;
        g_test_data[3 * nobs + i] = z1;
        g_test_data[4 * nobs + i] = z2;
        g_test_data[5 * nobs + i] = id;
        g_test_data[6 * nobs + i] = time;
        g_test_data[7 * nobs + i] = weight;
        g_test_data[8 * nobs + i] = id;    /* cluster = id */
        g_test_data[9 * nobs + i] = time;  /* cluster2 = time */
        g_test_data[10 * nobs + i] = x_exog * 2;  /* extra1 */
        g_test_data[11 * nobs + i] = z1 + z2;     /* extra2 */
    }

    /* Reset and setup scalars for civreghdfe */
    reset_scalars();

    SF_scal_save("__civreghdfe_K_endog", 1);
    SF_scal_save("__civreghdfe_K_exog", 1);
    SF_scal_save("__civreghdfe_K_iv", 2);
    SF_scal_save("__civreghdfe_G", 1);  /* One FE (id) */

    /* Vary VCE options based on test_variant */
    int vce_type = test_variant % 3;  /* 0=unadjusted, 1=robust, 2=cluster */
    int has_cluster = (vce_type == 2) ? 1 : 0;

    SF_scal_save("__civreghdfe_has_cluster", has_cluster);
    SF_scal_save("__civreghdfe_has_cluster2", 0);
    SF_scal_save("__civreghdfe_has_weights", (test_variant % 5 == 0) ? 1 : 0);
    SF_scal_save("__civreghdfe_weight_type", 1);  /* fweight */
    SF_scal_save("__civreghdfe_vce_type", vce_type);
    SF_scal_save("__civreghdfe_maxiter", 16000);
    SF_scal_save("__civreghdfe_tolerance", 1e-8);
    SF_scal_save("__civreghdfe_verbose", 0);
    SF_scal_save("__civreghdfe_est_method", test_variant % 3);  /* 0=2SLS, 1=LIML, 2=Fuller */
    SF_scal_save("__civreghdfe_kclass", 1.0);
    SF_scal_save("__civreghdfe_fuller_alpha", 1.0);
    SF_scal_save("__civreghdfe_kernel_type", 0);
    SF_scal_save("__civreghdfe_bw", 0);
    SF_scal_save("__civreghdfe_kiefer", 0);
    SF_scal_save("__civreghdfe_dofminus", 0);
    SF_scal_save("__civreghdfe_sdofminus", 0);
    SF_scal_save("__civreghdfe_nopartialsmall", 0);
    SF_scal_save("__civreghdfe_center", 0);
    SF_scal_save("__civreghdfe_partial_nrows", 0);
}

/* Forward declaration */
extern ST_retcode civreghdfe_main(const char *args);

/* Stub functions for missing symbols */
void cmerge_cleanup_cache(void) {}
void cimport_cleanup_cache(void) {}
void cexport_cleanup_state(void) {}

/* Stub sort functions - not needed for civreghdfe testing */
stata_retcode ctools_sort_radix_lsd_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}
stata_retcode ctools_sort_radix_msd_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}
stata_retcode ctools_sort_timsort_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}
stata_retcode ctools_sort_merge_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}
stata_retcode ctools_sort_ips4o_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}
stata_retcode ctools_sort_counting_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}
stata_retcode ctools_sort_sample_order_only(stata_data *data, int *sort_vars, size_t nsort) {
    (void)data; (void)sort_vars; (void)nsort; return STATA_OK;
}

/* ============================================================================
 * Main Test Driver
 * ============================================================================ */

int main(int argc, char **argv) {
    int num_iterations = 150;  /* More than the 126 tests in validation */

    if (argc > 1) {
        num_iterations = atoi(argv[1]);
    }

    printf("===========================================\n");
    printf("civreghdfe ASan Test Harness\n");
    printf("Running %d iterations\n", num_iterations);
    printf("===========================================\n");

#ifdef _OPENMP
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#else
    printf("OpenMP: disabled\n");
#endif

    for (int i = 0; i < num_iterations; i++) {
        /* Vary test parameters */
        int nobs = 200 + (i % 100);
        setup_test_data(nobs, 12345 + i, i);

        ST_retcode rc = civreghdfe_main("iv_regression");

        if (i % 25 == 0 || rc != STATA_OK) {
            printf("Iteration %3d: nobs=%d, rc=%d\n", i + 1, nobs, rc);
        }
    }

    if (g_test_data) {
        free(g_test_data);
        g_test_data = NULL;
    }

    printf("===========================================\n");
    printf("All %d iterations complete!\n", num_iterations);
    printf("===========================================\n");

    return 0;
}
