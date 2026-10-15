"""Native regressions for SEP24: allocation failure, sorting, pool shutdown, file publication.

Optional CTOOLS_TEST_OPENMP_PREFIX builds the same sort tests with OpenMP and
constrained teams. All binaries run under UBSan with bounded timeouts; set CTOOLS_SANITIZERS=address,undefined to also enable ASan.
"""
from native.support.compiler import compile_native
from native.support.runner import harness
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
ALLOC=harness("sep24_native_alloc.c")

VCE=ALLOC+harness("sep24_native_5.c")

ARENA=ALLOC+harness("sep24_native_6.c")

POOL=harness("sep24_native_pool.c")

SORT_FAULT=ALLOC+harness("sep24_native_7.c")
SORT=harness("sep24_native_sort.c")

ATOMIC=harness("sep24_native_atomic.c")

KEY_GUARD=ALLOC+harness("sep24_native_8.c")

RANK_SOLVE=ALLOC+harness("sep24_native_9.c")

def compile_run(tmp,name,source,flags=(),env=None,args=()):
    cfile=tmp/(name+'.c');binary=tmp/name;cfile.write_text(source)
    base=['-std=c11','-O1','-g','-fno-fast-math','-fno-omit-frame-pointer','-fsanitize='+os.environ.get('CTOOLS_SANITIZERS','undefined'),'-I',str(ROOT/'src')]
    if platform.system()=='Darwin':base+=['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup']
    else:base+=['-DSYSTEM=STUNIX','-D_POSIX_C_SOURCE=200809L','-ffunction-sections','-Wl,--gc-sections']
    compile_native([os.environ.get('CC','/usr/bin/clang' if platform.system()=='Darwin' else 'clang'),*base,str(cfile),*flags,'-pthread','-lm','-o',str(binary)],check=True,capture_output=True,text=True)
    run_env=dict(os.environ,ASAN_OPTIONS='detect_leaks=0:symbolize=0:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:symbolize=0',**(env or {}))
    subprocess.run([str(binary),*map(str,args)],check=True,env=run_env,timeout=45,capture_output=True,text=True)
    print('PASS',name,flush=True)

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-sep24-native-') as directory:
        tmp=Path(directory)
        for name,source in [('key_type_guard',KEY_GUARD),('rank_solve',RANK_SOLVE),('vce_faults',VCE),('arena_reset',ARENA),('pool_create_failure',POOL),('atomic_output',ATOMIC)]:
            compile_run(tmp,name,source,args=[tmp/'literal threads(2) path'] if name=='atomic_output' else [])
        for alg,num,string in [('sample','sample_sort_numeric_impl(order,keys,128,4)','sample_sort_string_impl(order,strings,128,4,2)'),('merge','parallel_merge_sort_numeric(order,keys,128,4)','parallel_merge_sort_string(order,strings,128,4,2)')]:
            compile_run(tmp,alg+'_faults',SORT_FAULT.replace('ALGORITHM',alg).replace('NUMERIC_CALL',num).replace('STRING_CALL',string))
        for alg in ['sample','merge','radix_lsd','counting','ips4o']:
            # Public LSD name omits "radix_".
            source=SORT.replace('ALGORITHM',alg)
            if alg == 'counting': source=source.replace('string<2','string<1')
            extra=[str(ROOT/'src/ctools_arena.c')] if alg in ('counting','ips4o') else []
            compile_run(tmp,alg+'_noomp',source,extra)
            prefix=os.environ.get('CTOOLS_TEST_OPENMP_PREFIX')
            if prefix:
                flags=['-Xpreprocessor','-fopenmp','-I',str(Path(prefix)/'include'),str(Path(prefix)/'lib/libomp.a')]
                for limit in ['1','2']:
                    compile_run(tmp,alg+'_team'+limit,source,flags+extra,{'OMP_THREAD_LIMIT':limit,'OMP_DYNAMIC':'TRUE'})

if __name__=='__main__':
    try:main()
    except subprocess.CalledProcessError as e:
        print(e.stdout or '',e.stderr or '',flush=True);raise
