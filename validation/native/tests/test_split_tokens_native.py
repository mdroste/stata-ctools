"""Sanitizer coverage for csplit's in-place single-delimiter storage path.

Reuse the existing SPI/allocator harness without changing its assertions.
"""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness
import os
from pathlib import Path
import platform
import subprocess
import tempfile
from validation.native.support.runner import ROOT, harness
SOURCE = harness("newcommands_native_source.c")

EXTRA = harness("split_tokens_native_extra.c")

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-split-tokens-') as tmp:
        cfile=Path(tmp)/'test.c';binary=Path(tmp)/'test'
        cfile.write_text('#define CTOOLS_HARNESS_MAIN original_main\n' + SOURCE + '\n#undef CTOOLS_HARNESS_MAIN\n' + EXTRA)
        flags=['-std=c11','-O1','-g','-fsanitize='+os.environ.get('CTOOLS_SANITIZERS','undefined'),
               '-fno-omit-frame-pointer','-fno-fast-math','-ffp-contract=off','-Wno-unknown-pragmas','-I',str(ROOT/'src')]
        flags+=['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE'] if platform.system()=='Darwin' else ['-DSYSTEM=STUNIX','-D_GNU_SOURCE']
        compile_native([os.environ.get('CC','clang'),*flags,str(cfile),str(ROOT/'src/ctools_parse.c'),'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=60)

if __name__=='__main__':
    main()
