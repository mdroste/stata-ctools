from pathlib import Path
import sys,subprocess
ROOT=Path('/Users/mike/Documents/GitHub/stata-ctools')
sys.path.insert(0,str(ROOT/'validation'))
from native.support.compiler import compile_native
flags=['-std=c11','-O1','-g','-DSD_FASTMODE','-D_GNU_SOURCE','-D_DARWIN_C_SOURCE','-pthread','-fsanitize=address,undefined','-fno-omit-frame-pointer','-I',str(ROOT/'src'),'-DSYSTEM=APPLEMAC','-Wl,-dead_strip']
sources=[Path('/tmp/ctools-io-audit/probe.c'),*[ROOT/'src'/p for p in ['formats/xlsx_xml.c','ctools_types.c','ctools_arena.c','ctools_threads.c']],*sorted((ROOT/'src/vendor/miniz').glob('*.c')),*sorted((ROOT/'src/vendor/libdeflate').rglob('*.c'))]
compile_native(['cc',*flags,*map(str,sources),'-o','/tmp/ctools-io-audit/probe'],check=True)
subprocess.run(['/tmp/ctools-io-audit/probe'],check=True)
