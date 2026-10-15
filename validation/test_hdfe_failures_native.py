"""IV allocation ownership and batch binscatter convergence regressions."""

import os
import platform
import re
import subprocess
import tempfile
from pathlib import Path

from native.support.compiler import compile_native

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "validation/native/memory_safety"


def main():
    apple = platform.system() == "Darwin"
    env = dict(os.environ, ASAN_OPTIONS="detect_leaks=0", UBSAN_OPTIONS="halt_on_error=1")
    with tempfile.TemporaryDirectory(prefix="ctools-hdfe-failures-") as tmp:

        def build(name):
            binary = Path(tmp) / name
            flags = (
                ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip"]
                if apple
                else [
                    "-DSYSTEM=OPUNIX",
                    "-ffunction-sections",
                    "-fdata-sections",
                    "-Wl,--gc-sections",
                ]
            )
            compile_native(
                [
                    os.environ.get("CC", "/usr/bin/clang" if apple else "cc"),
                    "-std=gnu11",
                    "-O1",
                    "-g",
                    "-fno-fast-math",
                    "-ffp-contract=off",
                    "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "address,undefined"),
                    "-I",
                    str(ROOT / "src"),
                    "-I",
                    str(FIXTURES),
                    *flags,
                    str(FIXTURES / (name + ".c")),
                    "-lm",
                    "-pthread",
                    "-o",
                    str(binary),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=ROOT,
            )
            return binary

        def run(binary, *args):
            result = subprocess.run(
                [str(binary), *map(str, args)],
                check=True,
                capture_output=True,
                text=True,
                env=env,
                timeout=30,
            )
            return int(re.search(r"allocations=(\d+)", result.stdout)[1])

        iv = build("iv")
        cases = 0
        for variant in [(0, 0, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]:
            count = run(iv, 0, *variant)
            for failure in range(1, count + 1):
                run(iv, failure, *variant)
                cases += 1
            print(
                f"PASS IV {variant}: {count} allocation failures, estimates and recovery",
                flush=True,
            )
        bins = build("bins_hdfe")
        for weighted in (0, 1):
            for first in (0, 1, 2):
                run(bins, 0, weighted, first)
            count = run(bins, 0, weighted)
            for failure in range(1, count + 1):
                run(bins, failure, weighted)
                cases += 1
            for maxiter in (0, 1):
                run(bins, 0, weighted, 0, maxiter)
        print(
            f"PASS batch HDFE: dummy-OLS parity, weighted orthogonality, column order, iteration limits; {cases} fault injections",
            flush=True,
        )


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(exc.stdout or "")
        print(exc.stderr or "")
        raise
