"""ASan/UBSan checks of sample selection, cluster identity, and merge contracts."""
from native.support.compiler import compile_native
from native.support.runner import harness
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

SOURCE = harness("sep22_native_source.c")

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-sep22-native-') as tmp:
        cfile=Path(tmp)/'test.c'; binary=Path(tmp)/'test'
        cfile.write_text(SOURCE)
        flags=['-std=c11','-O1','-g','-fno-fast-math','-fno-omit-frame-pointer',
               '-fsanitize=address,undefined','-I',str(ROOT/'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup']
        else:
            flags += ['-DSYSTEM=STUNIX','-ffunction-sections','-fdata-sections','-Wl,--gc-sections']
        compile_native([os.environ.get('CC', '/usr/bin/clang' if platform.system() == 'Darwin' else 'clang'),*flags,str(cfile),'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=60)
        print('PASS: full-double clusters; original-row selection; singleton chains; every selection allocation failure (ASan/UBSan)')

if __name__ == '__main__': main()
