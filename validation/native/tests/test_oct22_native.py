"""Regression tests for the October 22 audit."""
from pathlib import Path
import os
import platform
import subprocess
import tempfile
from validation.native.support.compiler import compile_native

from validation import ROOT
def main():
    with tempfile.TemporaryDirectory(prefix='ctools-oct22-') as tmp:
        directory = Path(tmp)
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-D_GNU_SOURCE', '-D_DARWIN_C_SOURCE',
                 '-pthread', '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-I', str(ROOT/'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip']
            libraries = []
        else:
            flags += ['-DSYSTEM=STUNIX', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections']
            libraries = ['-lm']
        for name in ['xml', 'export', 'match_oom', 'range']:
            sources = [ROOT/('validation/native/harnesses/oct22_'+name+'.c')]
            if name == 'xml':
                sources += [ROOT/'src'/p for p in ['formats/xlsx_xml.c','ctools_types.c','ctools_arena.c','ctools_threads.c']]
                sources += sorted((ROOT/'src/vendor/miniz').glob('*.c')) + sorted((ROOT/'src/vendor/libdeflate').rglob('*.c'))
            if name == 'export':
                sources += [ROOT/'src/ctools_hash.c', ROOT/'src/ctools_arena.c']
            binary = directory/name
            compile_native([os.environ.get('CC','cc'), *flags, *map(str,sources), *libraries, '-o',str(binary)],check=True)
            subprocess.run([str(binary)],check=True,env={**os.environ,'UBSAN_OPTIONS':'halt_on_error=1'})

if __name__=='__main__': main()
