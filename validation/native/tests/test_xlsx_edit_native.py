"""ASan/UBSan checks for lossless worksheet/style merge operations."""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness
import os
import platform
import subprocess
import tempfile
from pathlib import Path

from validation import ROOT
SOURCE = harness("xlsx_edit_native_source.c")

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-xlsx-edit-') as tmp:
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
        compile_native([os.environ.get('CC', 'cc'), *flags, str(source),
                        str(ROOT / 'src/formats/xlsx_xml.c'), '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})
        print('PASS XLSX edits preserve styles/formulas/XML, namespaces, truncated tags (ASan/UBSan)')

if __name__ == '__main__':
    main()
