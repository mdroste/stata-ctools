"""Exercise stable numeric radix and tiled merges against an independent qsort.

Uses the real OpenMP runtime (including teams smaller than the requested count).
Set CTOOLS_LIBOMP_PREFIX on macOS. --source-root can test a frozen baseline.
"""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness, openmp_prefix
from pathlib import Path
import argparse
import os
import platform
import subprocess
import tempfile

from validation import ROOT
SOURCE = harness("order_parallel_source.c")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-root',type=Path,default=ROOT)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='ctools-order-') as tmp:
        cfile=Path(tmp)/'test.c'; binary=Path(tmp)/'test'
        cfile.write_text(SOURCE)
        flags=['-std=c11','-O2','-g','-fsanitize=undefined','-fno-fast-math','-ffp-contract=off',
               '-Wall','-Wextra','-Werror','-I',str(args.source_root/'src')]
        if platform.system()=='Darwin':
            prefix=openmp_prefix()
            flags+=['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Xpreprocessor','-fopenmp','-I',str(prefix/'include')]
            libs=[str(prefix/'lib/libomp.a')]
        else:
            flags+=['-DSYSTEM=STUNIX','-D_GNU_SOURCE','-fopenmp'];libs=[]
        compile_native([os.environ.get('CC','clang'),*flags,str(cfile),str(args.source_root/'src/ctools_order.c'),*libs,'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=120)
        subprocess.run([str(binary)],env={**os.environ,'OMP_THREAD_LIMIT':'2'},check=True,timeout=120)

if __name__=='__main__':
    main()
