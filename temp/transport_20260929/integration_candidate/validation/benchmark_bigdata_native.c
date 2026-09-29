/* Key-kernel scaling only: excludes Stata I/O and wide payload columns.
 * Arguments: rows, threads, repetitions, sort_algorithm_t code, workload.
 * Workloads: 0 = 10,000 groups; 1 = large integers; 2 = doubles; 3 = m:1 join.
 */
#include "cmerge/cmerge_join.h"
#include "ctools_config.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double seconds(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1e-9 * t.tv_nsec;
}
ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_boolean missing(ST_double x) { return x >= 8.98846567431158e307; }
int main(int argc, char **argv) {
  assert(argc == 6);
  size_t n = strtoull(argv[1], 0, 10);
  int threads = atoi(argv[2]);
  int repeats = atoi(argv[3]);
  int algorithm = atoi(argv[4]);
  int kind = atoi(argv[5]);
  mock.missval = 8.98846567431158e307;
  mock.ismissing = missing;
  ctools_set_max_threads(threads);
  if (kind == 3) {
    size_t nu = (n + 9) / 10;
    double *m = malloc(n * sizeof(double)), *u = malloc(nu * sizeof(double));
    assert(m && u);
    for (size_t i = 0; i < n; i++)
      m[i] = (double)(i / 10);
    for (size_t i = 0; i < nu; i++)
      u[i] = (double)i;
    stata_variable mv = {.type = STATA_TYPE_DOUBLE, .nobs = n, .data.dbl = m};
    stata_variable uv = {.type = STATA_TYPE_DOUBLE, .nobs = nu, .data.dbl = u};
    stata_data md = {.nobs = n, .nvars = 1, .vars = &mv},
               ud = {.nobs = nu, .nvars = 1, .vars = &uv};
    for (int rep = 0; rep < repeats; rep++) {
      cmerge_output_spec_t *out = NULL;
      double start = seconds();
      int64_t count = cmerge_sorted_join(&md, &ud, 1, MERGE_M_1, &out);
      double elapsed = seconds() - start;
      assert(count == (int64_t)n);
      for (size_t i = 0; i < n; i++)
        assert(out[i].master_sorted_row == (int32_t)i &&
               out[i].using_sorted_row == (int32_t)(i / 10) &&
               out[i].merge_result == 3);
      printf("NATIVE,%zu,%d,%d,%d,%d,%.9f\n", n, threads, rep, algorithm, kind,
             elapsed);
      fflush(stdout);
      free(out);
    }
    free(m);
    free(u);
    return 0;
  }
  double *values = malloc(n * sizeof(double));
  perm_idx_t *order = malloc(n * sizeof(perm_idx_t));
  unsigned char *seen = malloc(n);
  assert(values && order && seen);
  stata_variable var = {
      .nobs = n, .type = STATA_TYPE_DOUBLE, .data.dbl = values};
  stata_data data = {.nobs = n, .nvars = 1, .vars = &var, .sort_order = order};
  int key = 1;
  uint64_t state = 24092026;
  for (size_t i = 0; i < n; i++) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    values[i] = kind == 0 ? (double)(state % 10000)
                          : (kind == 1 ? (double)(state % (n * 4))
                                       : (double)(state >> 11) * 0x1p-53);
  }
  for (int rep = 0; rep < repeats; rep++) {
    for (size_t i = 0; i < n; i++)
      order[i] = i;
    double start = seconds();
    assert(ctools_sort_dispatch(&data, &key, 1, (sort_algorithm_t)algorithm) ==
           0);
    double elapsed = seconds() - start;
    memset(seen, 0, n);
    for (size_t i = 0; i < n; i++) {
      assert(order[i] < n && !seen[order[i]]);
      seen[order[i]] = 1;
      if (i)
        assert(values[order[i - 1]] <= values[order[i]]);
    }
    printf("NATIVE,%zu,%d,%d,%d,%d,%.9f\n", n, threads, rep, algorithm, kind,
           elapsed);
    fflush(stdout);
  }
  free(values);
  free(order);
  free(seen);
  return 0;
}
