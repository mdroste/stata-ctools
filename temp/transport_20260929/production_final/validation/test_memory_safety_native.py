"""Command-level leak, allocation failure, store failure and overflow regressions.

Compile actual production sources with a mocked Stata host. A tracked allocator
checks local/global ownership, while ASan/UBSan check invalid access/arithmetic.
Each regression and PPML failure is followed by a successful call in-process.
"""
import os
from pathlib import Path
import platform
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "validation/native/memory_safety"
SYSTEM = "APPLEMAC" if platform.system() == "Darwin" else "OPUNIX"
CC = os.environ.get("CC", "/usr/bin/clang" if SYSTEM == "APPLEMAC" else "cc")
SAN = os.environ.get("CTOOLS_SANITIZERS", "address,undefined,float-cast-overflow")
ENV = dict(os.environ, UBSAN_OPTIONS="halt_on_error=1", ASAN_OPTIONS="detect_leaks=0", OMP_NUM_THREADS="4", OMP_THREAD_LIMIT="4")

def main():
    with tempfile.TemporaryDirectory(prefix="ctools-memory-safety-") as tmp:
        def build(name, source=None, flags=()):
            exe = Path(tmp) / name
            omp = os.environ.get("CTOOLS_TEST_OPENMP_PREFIX")
            parallel = (["-Xpreprocessor", "-fopenmp", "-I"+omp+"/include", omp+"/lib/libomp.a"]
                        if omp else (["-fopenmp"] if SYSTEM != "APPLEMAC" else []))
            command = [CC, "-std=gnu11", "-O1", "-g", f"-DSYSTEM={SYSTEM}",
                       "-I", str(ROOT / "src"), "-I", str(FIXTURES),
                       f"-fsanitize={SAN}", "-ffunction-sections", "-fdata-sections",
                       *flags, *parallel, str(FIXTURES / f"{source or name}.c"), "-o", str(exe),
                       "-Wl,-dead_strip" if SYSTEM == "APPLEMAC" else "-Wl,--gc-sections",
                       "-lm", "-pthread"]
            subprocess.run(command, check=True, cwd=ROOT, capture_output=True, text=True)
            return exe

        def run(exe, *args):
            result = subprocess.run([str(exe), *map(str, args)], cwd=ROOT, env=ENV,
                                    check=True, capture_output=True, text=True, timeout=30)
            m = re.search(r"rc=(\d+) (?:allocations|calls)=(\d+) alive=(\d+)", result.stdout)
            if m:
                rc, calls, alive = map(int, m.groups())
                assert alive == 0, result.stdout
                return rc, calls
            return None

        cases = 0
        for name in ("regress", "ppml", "sampling", "bootstrap", "rangestat"):
            exe = build(name, "sampling" if name == "bootstrap" else name,
                        ["-DBOOTSTRAP"] if name == "bootstrap" else [])
            variants = [(0,0,0), (1,1,0), (1,2,0), (1,1,2), (0,0,1), (1,1,3)] if name == "regress" else (
                       [(0,0,0), (1,1,0)] if name == "ppml" else [()])
            for variant in variants:
                rc, count = run(exe, 0, *variant)
                assert rc == 0, (name, variant, rc)
                for failure in range(1, count + 1):
                    rc, _ = run(exe, failure, *variant)
                    # Optional accelerator allocations may fall back successfully.
                    assert rc in (0, 920, 430), (name, variant, failure, rc)
                    cases += 1
                if name != "ppml":
                    for store in (1, 6, 12):
                        rc, _ = run(exe, 0, *variant, store)
                        assert rc == 459, (name, variant, store, rc)
                        cases += 1
            if name == "regress":
                for scalar in (1, 8, 20):
                    assert run(exe, 0, 1, 1, 0, 0, scalar)[0] == 459
                    cases += 1
                for matrix in (1, 2, 4):
                    assert run(exe, 0, 1, 1, 0, 0, 0, matrix)[0] == 459
                    cases += 1
            print(f"{name}: allocation and output failures passed", flush=True)
        for name in ("remap", "timsort", "threads"):
            run(build(name))
            print(f"{name}: passed", flush=True)
        exe = build("styles")
        for failure in range(5):
            run(exe, failure)
        print(f"styles: all parser/growth failures passed; {cases} command injections", flush=True)

if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(exc.stdout or "")
        print(exc.stderr or "")
        raise
