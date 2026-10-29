"""Standalone P1 arithmetic and allocation-failure checks (no Stata process)."""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness
import os
from pathlib import Path
import platform
import subprocess
import tempfile

from validation import ROOT
CC = os.environ.get("CC", "clang")
COMMON = harness("p1_native_common.c")


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
    compile_native([CC, *flags, str(cfile), *libraries, "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True, timeout=60)
    print(f"PASS {name}", flush=True)


with tempfile.TemporaryDirectory(prefix="ctools-p1-native-") as tmp:
    directory = Path(tmp)
    run_test(directory, "compensated_sum", harness("p1_native_2.c"))
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
    run_test(directory, "permutation_allocations", COMMON + harness("p1_native_3.c"), openmp=True)
