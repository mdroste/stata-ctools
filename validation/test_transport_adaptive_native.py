"""Adaptive transport regressions against an independent deterministic SPI model.

Requires the adaptive-load changes: opt-in CTOOLS_LOAD_NO_SORT_ORDER, numeric
row tiles for wide data, and packed maximum-width string arenas with retry.
An older source snapshot is expected to fail rather than silently skip them.

CTOOLS_TEST_SOURCE selects an isolated source tree. Serial and OpenMP builds
run under UBSan by default; CTOOLS_SANITIZERS=address,undefined also enables
ASan on systems with a working runtime. No Stata executable is invoked.
"""
from native.support.compiler import compile_native
from native.support.runner import harness
from pathlib import Path
import os
import subprocess
import tempfile
from native.support.runner import SPI, ROOT

SOURCE = SPI + harness("transport_adaptive_native_1.c")


def main():
    source_dir = Path(os.environ.get('CTOOLS_TEST_SOURCE', ROOT / 'src'))
    sanitizers = os.environ.get('CTOOLS_SANITIZERS', 'undefined')
    with tempfile.TemporaryDirectory(prefix='ctools-adaptive-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        for omp in (False, True):
            binary = directory / ('test-' + str(omp))
            flags = ['-std=c11', '-O2', '-DSD_FASTMODE', '-fno-fast-math',
                     '-fsanitize=' + sanitizers, '-fno-omit-frame-pointer',
                     '-pthread', '-I' + str(source_dir)]
            if os.uname().sysname == 'Darwin':
                flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE',
                          '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup']
                if omp:
                    prefix = os.environ.get('LIBOMP_PREFIX', '/opt/homebrew/opt/libomp')
                    flags += ['-Xpreprocessor', '-fopenmp', '-I' + prefix + '/include',
                              prefix + '/lib/libomp.a']
            else:
                flags += ['-DSYSTEM=OPUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                          '-fdata-sections', '-Wl,--gc-sections', '-lm']
                if omp:
                    flags += ['-fopenmp']
            compile_native([os.environ.get('CC', 'clang'), *flags, str(source),
                            *[str(source_dir / name) for name in
                              ['ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']],
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=120,
                           env=dict(os.environ, UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',
                                    KMP_BLOCKTIME='0', OMP_WAIT_POLICY='PASSIVE'))
            print(f'PASS adaptive transport sanitizers={sanitizers} '
                  f'openmp={omp} source={source_dir}', flush=True)


if __name__ == '__main__':
    main()
