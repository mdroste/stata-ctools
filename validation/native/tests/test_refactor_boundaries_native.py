"""Public header isolation and shared utility contracts introduced by the refactor."""

from validation.native.support.runner import CC, ROOT, harness, run_test
from pathlib import Path
import subprocess
import tempfile
import platform


def main():
    subprocess.run(["python3", str(ROOT / "scripts/sync_architecture.py"), "--check"], check=True)
    headers = [
        "ctools_data.h",
        "ctools_data_io.h",
        "ctools_sort.h",
        "ctools_convert.h",
        "ctools_alloc.h",
        "ctools_platform.h",
        "ctools_rng.h",
        "ctools_group.h",
        "ctools_linalg.h",
        "ctools_text.h",
        "formats/ctools_output.h",
        "formats/excel_address.h",
        "hdfe/hdfe_types.h",
        "hdfe/hdfe_solver.h",
        "hdfe/hdfe_remap.h",
    ]
    with tempfile.TemporaryDirectory(prefix="ctools-boundaries-") as tmp:
        directory = Path(tmp)
        for header in headers:
            source = directory / "header.c"
            source.write_text(f'#include "{header}"\nint main(void) {{ return 0; }}\n')
            system = "APPLEMAC" if platform.system() == "Darwin" else "STUNIX"
            subprocess.run(
                [
                    CC,
                    "-std=c11",
                    "-fsyntax-only",
                    "-DSYSTEM=" + system,
                    "-I",
                    str(ROOT / "src"),
                    str(source),
                ],
                check=True,
            )
        print("PASS standalone public headers", flush=True)
        # Header extraction must not preempt a translation unit's SPI mode.
        # Inspect actual production preprocessing under the release-mode flags.
        for name in ("creghdfe/creghdfe_regress.c", "cppmlhdfe/cppmlhdfe_irls.c"):
            macros = subprocess.check_output(
                [
                    CC,
                    "-std=c11",
                    "-E",
                    "-dM",
                    "-DSD_FASTMODE",
                    "-DSYSTEM=" + system,
                    "-I",
                    str(ROOT / "src"),
                    str(ROOT / "src" / name),
                ],
                text=True,
            )
            for macro, callback in (("SF_vdata", "safevdata"), ("SF_vstore", "safestore")):
                definition = next(
                    line
                    for line in macros.splitlines()
                    if line.startswith("#define " + macro + "(")
                )
                assert callback in definition, (name, definition)
        print("PASS estimator checked-SPI include order", flush=True)
        run_test(
            directory,
            "refactor_boundaries",
            harness("refactor_boundaries.c"),
            ["ctools_group.c", "ctools_platform.c", "ctools_runtime.c", "formats/excel_address.c"],
        )


if __name__ == "__main__":
    main()
