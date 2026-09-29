"""Small native cqreg regressions and optional before/after timings.

Default checks use at most 8,193 values and 303-observation fits. Timings default
to 20,000 rows, capped at 200,000; no Stata installation is required.

  python3 validation/test_cqreg_native.py
  python3 validation/test_cqreg_native.py --benchmark --rows 200000 --columns 12
  python3 validation/test_cqreg_native.py --benchmark --baseline-dir /tmp/base/cqreg

The baseline directory must contain a saved copy of src/cqreg from before edits.
Set CTOOLS_TEST_OPENMP_PREFIX to test with OpenMP. Set CTOOLS_SANITIZERS to
address,undefined to check memory access as well as numerical results.
"""
import argparse
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def build(directory, output, *, baseline=False, fallback=False):
    flags = ["-O3", "-fno-fast-math", "-ffp-contract=off", "-funroll-loops",
             "-Wno-deprecated-declarations", "-I" + str(directory),
             "-I" + str(ROOT / "src"), "-I" + str(ROOT / "src/cqreg")]
    link = ["-lm"]
    if platform.system() == "Darwin":
        flags += ["-DSYSTEM=APPLEMAC"]
        link += ["-framework", "Accelerate"]
    else:
        flags += ["-DSYSTEM=STUNIX"]
    if fallback:
        flags += ["-DUSE_BLAS=0"]
    if baseline:
        flags += ["-DCQREG_BASELINE"]
    if os.environ.get("CTOOLS_SANITIZERS"):
        flags += ["-fsanitize=" + os.environ["CTOOLS_SANITIZERS"],
                  "-fno-sanitize-recover=all", "-g"]
    prefix = os.environ.get("CTOOLS_TEST_OPENMP_PREFIX")
    if prefix:
        flags += ["-Xpreprocessor", "-fopenmp", "-I" + prefix + "/include"]
        link += [prefix + "/lib/libomp.a"]
    elif platform.system() != "Darwin":
        flags += ["-fopenmp"]
    sources = [ROOT / "validation/cqreg_native.c"]
    # Count requested aligned workspace bytes without allocating huge datasets.
    counted_types = output.with_suffix(".c")
    counted_types.write_text(
        '#include <stdlib.h>\n'
        'extern int cqreg_record_alloc(void **, size_t, size_t);\n'
        '#define posix_memalign cqreg_record_alloc\n'
        f'#include "{directory / "cqreg_types.c"}"\n')
    sources += [counted_types]
    sources += [directory / ("cqreg_" + s + ".c")
                for s in ("linalg", "blas", "sparsity")]
    subprocess.run([os.environ.get("CC", "cc"), *flags, *map(str, sources),
                    *link, "-o", str(output)], check=True, timeout=60)


def run(binary, args=()):
    env = dict(os.environ, OMP_NUM_THREADS="4")
    result = subprocess.run([str(binary), *map(str, args)], env=env,
                            check=True, capture_output=True, text=True, timeout=30)
    print(result.stdout, end="", flush=True)
    return [line.split(",") for line in result.stdout.splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--rows", type=int, default=20000)
    parser.add_argument("--columns", type=int, default=6)
    parser.add_argument("--baseline-dir", type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if not 100 <= args.rows <= 200000 or not 2 <= args.columns <= 32:
        parser.error("use 100..200000 rows and 2..32 columns to keep runs bounded")
    if not 1 <= args.repetitions <= 5:
        parser.error("use 1..5 repetitions")
    with tempfile.TemporaryDirectory(prefix="cqreg-native-") as tmp:
        tmp = Path(tmp)
        current = tmp / "current"
        build(ROOT / "src/cqreg", current)
        run(current)
        fallback = tmp / "fallback"
        build(ROOT / "src/cqreg", fallback, fallback=True)
        run(fallback)
        if args.benchmark:
            if args.baseline_dir:
                baseline = tmp / "baseline"
                build(args.baseline_dir.resolve(), baseline, baseline=True)
            new_times, old_times = {}, {}
            for _ in range(args.repetitions):
                current_results = run(current, (args.rows, args.columns))
                if args.baseline_dir:
                    old_results = run(baseline, (args.rows, args.columns))
                    assert len(current_results) == len(old_results)
                    for new, old in zip(current_results, old_results):
                        assert new[:3] == old[:3]
                        if new[0] == "workspace":
                            assert int(new[4]) < int(old[4])
                            print(f"workspace allocation saved: {int(old[4])-int(new[4]):,} bytes")
                            continue
                        assert new[4:] == old[4:], (new, old)
                        new_times.setdefault(new[0], []).append(float(new[3]))
                        old_times.setdefault(old[0], []).append(float(old[3]))
            for name in new_times:
                before = statistics.median(old_times[name])
                after = statistics.median(new_times[name])
                print(f"{name}: {before:.6f}s -> {after:.6f}s ({before/after:.2f}x)")


if __name__ == "__main__":
    main()
