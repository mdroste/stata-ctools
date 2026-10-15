"""ASan/UBSan checks for XLSX's dedicated length-aware strL loader."""
from native.support.compiler import compile_native
from native.support.runner import harness
import os
import platform
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = harness("xlsx_strl_native_source.c")


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-xlsx-strl-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        binary = directory / 'test'
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fsanitize=address,undefined',
                 '-fno-omit-frame-pointer', '-pthread', '-I', str(ROOT / 'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-Wl,-dead_strip',
                      '-Wl,-undefined,dynamic_lookup']
        else:
            flags += ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                      '-fdata-sections', '-Wl,--gc-sections', '-lm']
        modules = ['ctools_data_io.c', 'ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']
        compile_native([os.environ.get('CC', 'cc'), *flags, str(source),
                        *[str(ROOT / 'src' / name) for name in modules], '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})
        print('PASS XLSX long text, filtering, SPI failures, binary rejection, UTF-16 limit (ASan/UBSan)')


if __name__ == '__main__':
    main()
