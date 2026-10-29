"""Row-major staging (ctools_rowstage) against an independent SPI model.

Covers eligibility, staged contents and extracted columns, in-order and
reordered stores (repeated and absent sources, more output than staged rows),
range and index rejection before any write, read/write failures, cancellation,
stale width hints, and allocation failure, in serial and OpenMP builds.
Set CTOOLS_TEST_SOURCE to test an isolated source snapshot. Runs without Stata.
"""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness
from pathlib import Path
import os
import subprocess
import tempfile
from validation.native.support.runner import SPI, ROOT

SOURCE = SPI + harness("rowstage_native.c")


def main():
    source_dir = Path(os.environ.get('CTOOLS_TEST_SOURCE', ROOT / 'src'))
    with tempfile.TemporaryDirectory(prefix='ctools-rowstage-native-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        for omp in (False, True):
            binary = directory / ('test-' + str(omp))
            flags = ['-std=c11', '-O2', '-DSD_FASTMODE', '-fno-fast-math',
                     '-fsanitize=' + os.environ.get('CTOOLS_SANITIZERS', 'undefined'),
                     '-fno-omit-frame-pointer', '-pthread', '-I' + str(source_dir)]
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
                            *[str(source_dir / f) for f in
                              ['ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']],
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=120,
                           env=dict(os.environ, UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1'))
            print(f'PASS rowstage openmp={omp} source={source_dir}', flush=True)


if __name__ == '__main__':
    main()
