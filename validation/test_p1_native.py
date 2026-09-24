"""Standalone P1 arithmetic and allocation-failure checks (no Stata process)."""
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CC = os.environ.get("CC", "clang")
COMMON = r'''
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
static int calls, fail_at, live;
static void *fi_malloc(size_t n) {
    if (++calls == fail_at) return NULL;
    void *p = malloc(n);
    if (p) { memset(p, 0xa5, n); live++; }
    return p;
}
static void *fi_calloc(size_t n, size_t size) {
    if (++calls == fail_at) return NULL;
    void *p = calloc(n, size);
    if (p) live++;
    return p;
}
static int fi_aligned(void **p, size_t align, size_t n) {
    if (++calls == fail_at) return ENOMEM;
    int rc = posix_memalign(p, align, n);
    if (!rc) live++;
    return rc;
}
static void fi_free(void *p) { if (p) live--; free(p); }
#define malloc fi_malloc
#define calloc fi_calloc
#define posix_memalign fi_aligned
#define free fi_free
'''


def run_test(directory, name, source, *, openmp=False, sanitize=True):
    cfile = directory / f"{name}.c"
    binary = directory / name
    cfile.write_text(source)
    flags = ["-std=c11", "-O3",
             "-fno-fast-math", "-ffp-contract=off", "-I", str(ROOT / "src")]
    if sanitize:
        flags += ["-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
                  "-fno-omit-frame-pointer"]
    libraries = []
    if platform.system() == "Darwin":
        flags += ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE", "-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"]
        if openmp:
            prefix = os.environ.get("LIBOMP_PREFIX")
            if not prefix:
                prefix = subprocess.check_output(["brew", "--prefix", "libomp"], text=True).strip()
            flags += ["-Xpreprocessor", "-fopenmp", "-I", prefix + "/include",
                      "-L", prefix + "/lib", "-lomp"]
    else:
        flags += ["-DSYSTEM=STUNIX", "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections"]
        libraries = ["-lm"]
        if openmp:
            flags += ["-fopenmp"]
    subprocess.run([CC, *flags, str(cfile), *libraries, "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True, timeout=60)
    print(f"PASS {name}", flush=True)


with tempfile.TemporaryDirectory(prefix="ctools-p1-native-") as tmp:
    directory = Path(tmp)
    run_test(directory, "compensated_sum", r'''
#include <assert.h>
#include "ctools_ols.h"
__attribute__((noinline)) static double sum(const double *x, int n) {
    dd_real total = {0, 0};
    for (int i = 0; i < n; i++) total = dd_add_d(total, x[i]);
    return total.hi + total.lo;
}
int main(void) {
    volatile double a = 1e16, b = 1;
    dd_real pair = two_sum(a, b);
    assert(pair.lo == 1);
    double x[] = {1e16, 1, -1e16};
    assert(sum(x, 3) == 1);
    return 0;
}
''')
    for module, allocators in [
        ("ctools_sort_radix_lsd.c", ["radix_context", "string_context"]),
        ("ctools_sort_radix_msd.c", ["msd_context"]),
    ]:
        body = COMMON + f'\n#include "{module}"\nint main(void) {{\n'
        for allocator in allocators:
            body += f'''
    calls = 0; fail_at = 0;
    {allocator}_t *ctx_{allocator} = {allocator}_alloc(4);
    assert(ctx_{allocator});
    int count_{allocator} = calls;
    {allocator}_free(ctx_{allocator});
    assert(live == 0);
    for (int failure = 1; failure <= count_{allocator}; failure++) {{
        calls = 0; fail_at = failure;
        ctx_{allocator} = {allocator}_alloc(4);
        assert(!ctx_{allocator});
        assert(live == 0);
    }}
'''.replace(f"{allocator}_t", {
                "radix_context": "radix_sort_context_t",
                "string_context": "string_sort_context_t",
                "msd_context": "msd_sort_context_t",
            }[allocator])
        run_test(directory, "test_" + module[:-2], body + "return 0; }\n")
    run_test(directory, "permutation_allocations", COMMON + r'''
#include "ctools_types.c"
int ctools_get_max_threads(void) { return 4; }
int main(void) {
    stata_variable vars[4] = {0};
    perm_idx_t perm[2] = {1, 0};
    stata_data data = {0};
    data.nvars = 4; data.nobs = 2; data.vars = vars; data.sort_order = perm;
    /* Two pointer arrays and two aligned buffers for each of four workers. */
    for (int failure = 1; failure <= 10; failure++) {
        calls = 0; fail_at = failure;
        assert(ctools_apply_permutation(&data) == STATA_ERR_MEMORY);
        assert(live == 0);
    }
    return 0;
}
''', openmp=True)
