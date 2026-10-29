"""Native Windows runtime checks, invoked by the MSYS2 build job."""
from pathlib import Path
import os
import subprocess
import tempfile

from validation import ROOT
def main():
    if os.name != 'nt':
        raise SystemExit('This regression requires a real Windows worker')
    with tempfile.TemporaryDirectory(prefix='ctools-win-paths-') as tmp:
        binary = Path(tmp)/'paths.exe'
        subprocess.run([os.environ.get('CC', 'gcc'), '-std=c11', '-O1', '-DSYSTEM=STWIN32', '-DSD_FASTMODE',
                        '-I', str(ROOT/'src'), str(ROOT/'validation/native/harnesses/oct21_windows.c'),
                        str(ROOT/'src/formats/ctools_output.c'), '-o', str(binary)], check=True)
        subprocess.run([str(binary), tmp], check=True)
        build, source = Path(tmp)/'build', Path(tmp)/'src'
        build.mkdir(); source.mkdir()
        for name in ('ctools.plugin', 'ctools_windows.plugin', 'stale.o'): (build/name).touch()
        (source/'stale.o').touch()
        subprocess.run(['make', 'clean', 'DETECTED_OS=Windows', f'BUILD_DIR={build.as_posix()}',
                        f'SRC_DIR={source.as_posix()}'], cwd=ROOT, check=True)
        assert not list(build.iterdir()) and not list(source.iterdir())


if __name__ == '__main__':
    main()
