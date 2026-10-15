"""ASan/UBSan tests for native ReLU separation and heterogeneous projections."""
from native.support.compiler import compile_native
from native.support.runner import harness
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = harness("cpplmhdfe_native_source.c")

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-ppml-native-') as tmp:
        src=Path(tmp)/'test.c'; binary=Path(tmp)/'test'
        src.write_text(SOURCE)
        flags=['-std=c11','-O1','-g','-fno-fast-math','-fno-omit-frame-pointer',
               '-fsanitize=address,undefined','-I',str(ROOT/'src')]
        if platform.system()=='Darwin':
            flags += ['-DSYSTEM=APPLEMAC','-Wl,-dead_strip']
        else:
            flags += ['-DSYSTEM=STUNIX','-ffunction-sections','-fdata-sections','-Wl,--gc-sections']
        compile_native([os.environ.get('CC','/usr/bin/clang' if platform.system()=='Darwin' else 'clang'),*flags,str(src),'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=60)
        print('PASS: native regressor/FE separation, overlap, tiny positive outcomes, heterogeneous slopes, simplex, LSMR (ASan/UBSan)')

if __name__=='__main__': main()
