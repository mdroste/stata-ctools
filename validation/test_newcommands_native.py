"""UBSan and allocation-failure checks for the new command kernels.

The mock loader owns real heap copies; every allocation in each command is
failed in turn. This checks cleanup of staged caches and recovery on later calls.
Stata reference tests separately exercise the real SPI, syntax, and metadata.
Set CTOOLS_SANITIZERS=address,undefined to also run ASan (enabled in Linux CI).
"""
from native.support.compiler import compile_native
from native.support.runner import harness
from pathlib import Path
import os
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = harness("newcommands_native_source.c")


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-new-native-') as tmp:
        cfile = Path(tmp) / 'test.c'
        binary = Path(tmp) / 'test'
        cfile.write_text(SOURCE)
        flags = ['-std=c11', '-O1', '-g', '-fsanitize=' + os.environ.get('CTOOLS_SANITIZERS', 'undefined'),
                 '-fno-omit-frame-pointer', '-fno-fast-math', '-ffp-contract=off',
                 '-Wno-unknown-pragmas', '-I', str(ROOT / 'src')]
        flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE'] if platform.system() == 'Darwin' else ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE']
        compile_native([os.environ.get('CC','clang'), *flags, str(cfile),
                        str(ROOT/'src/ctools_parse.c'), '-lm', '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60)

if __name__ == '__main__':
    main()
