"""ASan/UBSan regression for BIFF continuation, formulas and OLE stream edits."""
from native.support.compiler import compile_native
from native.support.runner import harness
import os
import platform
import subprocess
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SOURCE = harness("xls_edit_native_source.c")

def main():
    vendor = ROOT / 'src/vendor'
    with tempfile.TemporaryDirectory(prefix='ctools-xls-edit-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'; source.write_text(SOURCE)
        binary = directory / 'test'
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                 '-I', str(ROOT/'src'), '-I', str(vendor/'libxls/include'),
                 '-I', str(vendor/'libxls/include/libxls'), '-D_DARWIN_C_SOURCE', '-D_GNU_SOURCE']
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup', '-liconv']
        else:
            flags += ['-DSYSTEM=STUNIX', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', '-lm']
        compile_native([os.environ.get('CC','cc'), *flags, str(source),
                        *map(str,sorted((vendor/'libxls/src').glob('*.c'))), '-o', str(binary)],check=True)
        subprocess.run([str(binary),str(directory)],check=True,timeout=60,
                       env={**os.environ,'UBSAN_OPTIONS':'halt_on_error=1'})
        print('PASS XLS long SST continuations, formulas, preserved OLE streams, record bounds (ASan/UBSan)')
if __name__ == '__main__':
    main()
