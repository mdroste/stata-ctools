"""LSD record sort of numeric and string keys against a stable reference sort.

Covers string keys of every length up to their width (empty, shared prefixes,
bytes above 127, NULL), extended missing values and -0, duplicates, mixed key
lists in both orders, and shuffled incoming orders, below and above the
parallel thresholds, in serial and OpenMP builds. Strings are allocated at their
exact length: run with CTOOLS_SANITIZERS=address,undefined to catch any read
past a terminator. Runs without Stata.
"""
from validation.native.support.compiler import compile_native
from validation.native.support.runner import harness, openmp_prefix, ROOT
import os
from pathlib import Path
import platform
import subprocess
import tempfile


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-sort-keys-native-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(harness('sort_keys_native.c'))
        for omp in (False, True):
            binary = directory / ('test-' + str(omp))
            flags = ['-std=c11', '-O1', '-g', '-fno-fast-math', '-fno-omit-frame-pointer',
                     '-fsanitize=' + os.environ.get('CTOOLS_SANITIZERS', 'undefined'),
                     '-pthread', '-I', str(ROOT / 'src')]
            if platform.system() == 'Darwin':
                flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-Wl,-dead_strip',
                          '-Wl,-undefined,dynamic_lookup']
                if omp:
                    prefix = str(openmp_prefix())
                    flags += ['-Xpreprocessor', '-fopenmp', '-I', prefix + '/include',
                              prefix + '/lib/libomp.a']
            else:
                flags += ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                          '-fdata-sections', '-Wl,--gc-sections', '-lm']
                if omp:
                    flags += ['-fopenmp']
            compile_native([os.environ.get('CC', 'clang'), *flags, str(source),
                            str(ROOT / 'src/ctools_sort_radix_lsd.c'), '-o', str(binary)],
                           check=True)
            subprocess.run([str(binary)], check=True, timeout=300,
                           env=dict(os.environ, ASAN_OPTIONS='detect_leaks=0:abort_on_error=1',
                                    UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1'))
            print(f'PASS sort_keys openmp={omp}', flush=True)


if __name__ == '__main__':
    main()
