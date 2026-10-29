"""Regression contracts for parsing, command sessions, ownership, and statuses."""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness
from pathlib import Path
import os
import platform
import subprocess
import tempfile

from validation import ROOT
SOURCE = harness("contracts_native_source.c")


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-contracts-') as tmp:
        c=Path(tmp)/'test.c';exe=Path(tmp)/'test';c.write_text(SOURCE)
        flags=['-std=c11','-O1','-g','-fsanitize='+os.environ.get('CTOOLS_SANITIZERS','undefined'),'-I'+str(ROOT/'src'),'-pthread','-ffunction-sections','-fdata-sections']
        flags += ['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE'] if platform.system()=='Darwin' else ['-DSYSTEM=STUNIX','-D_GNU_SOURCE']
        flags += ['-Wl,-dead_strip'] if platform.system()=='Darwin' else ['-Wl,--gc-sections']
        compile_native([os.environ.get('CC','clang'),*flags,str(c),*[str(ROOT/'src'/p) for p in ['ctools_parse.c','ctools_ols.c','ctools_runtime.c','ctools_threads.c','ctools_arena.c']],'-lm','-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)


if __name__=='__main__':main()
