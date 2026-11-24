"""Required cqreg result writes must propagate the original host error."""
import os
import platform
import subprocess
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import compile_native


def test_cqreg_publication_native():
    with tempfile.TemporaryDirectory(prefix="cqreg-publication-") as directory:
        output = Path(directory) / "publication"
        system = platform.system()
        flags = ["-std=gnu11", "-O2", "-flto", "-DSD_FASTMODE", "-I" + str(ROOT / "src")]
        if system == "Darwin":
            flags += ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"]
        else:
            flags += ["-DSYSTEM=" + ("STWIN32" if system == "Windows" else "STUNIX"),
                      "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections"]
        compile_native([os.environ.get("CC", "cc"), *flags,
                        str(ROOT / "validation/native/cqreg_publication.c"), "-o", str(output)], check=True)
        subprocess.run([str(output)], check=True, timeout=60)
