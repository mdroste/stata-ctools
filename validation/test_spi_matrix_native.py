"""Release-mode matrix writes must stay inside host-owned buffers."""
from native.support.compiler import compile_native
from native.support.runner import harness
from pathlib import Path
import os
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = harness("spi_matrix_native_source.c")


PUBLICATION = harness("spi_matrix_native_publication.c")


def main():
    system = "APPLEMAC" if platform.system() == "Darwin" else "OPUNIX"
    cc = os.environ.get("CC", "/usr/bin/clang" if system == "APPLEMAC" else "cc")
    with tempfile.TemporaryDirectory(prefix="ctools-spi-matrix-") as tmp:
        source = Path(tmp) / "matrix.c"
        source.write_text(SOURCE)
        for mode in ("SD_FASTMODE", "SD_SAFEMODE"):
            binary = Path(tmp) / mode
            compile_native([cc, "-std=c11", "-g", "-O1", "-fsanitize=address,undefined",
                            "-fno-omit-frame-pointer", *(["-DSD_FASTMODE"] if mode == "SD_FASTMODE" else []), f"-DSYSTEM={system}",
                            "-I", str(ROOT / "src"), str(source),
                            str(ROOT / "src/ctools_spi.c"), "-o", str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=30)
            print(f"PASS {mode}: invalid matrix writes, injected error, and recovery")
        source.write_text(PUBLICATION)
        binary = Path(tmp) / "publication"
        link = ["-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"] if system == "APPLEMAC" else ["-Wl,--gc-sections", "-lm"]
        compile_native([cc, "-std=gnu11", "-g", "-O1", "-fsanitize=address,undefined",
                        "-ffunction-sections", "-fdata-sections", "-DSD_FASTMODE", f"-DSYSTEM={system}",
                        "-I", str(ROOT / "src"), str(source), str(ROOT / "src/ctools_spi.c"),
                        *link, "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=30)
        print("PASS IV/quantile publication: allocation failure, 16 write failures, shapes, recovery")


if __name__ == "__main__":
    main()
