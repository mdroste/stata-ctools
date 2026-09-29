"""Build and run isolated sort/join scaling measurements, without launching Stata.

Example: python3 validation/benchmark_bigdata_native.py --rows 100000000
The Stata end-to-end harness is benchmark_bigdata.do; these are key kernels only.
"""
import argparse
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, default=1000000)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--algorithms', nargs='+', type=int, default=[7])
    parser.add_argument('--workloads', nargs='+', type=int, default=[0, 1, 2, 3])
    args = parser.parse_args()
    if not 0 < args.rows <= 2147483647 or args.threads < 1 or args.repetitions < 1:
        parser.error('positive rows, threads, and repetitions are required; rows must fit SPI')
    modules = ['ctools_types.c', 'ctools_threads.c', 'ctools_arena.c']
    modules += [p.name for p in sorted((ROOT / 'src').glob('ctools_sort_*.c'))]
    modules += ['cmerge/cmerge_join.c', 'cmerge/cmerge_keys.c', 'cmerge/cmerge_group_search.c']
    flags = ['-O3', '-fno-fast-math', '-ffp-contract=off', '-flto', '-DSD_FASTMODE', '-pthread', '-I'+str(ROOT/'src')]
    if platform.system() == 'Darwin':
        prefix = os.environ.get('LIBOMP_PREFIX') or subprocess.check_output(['brew', '--prefix', 'libomp'], text=True).strip()
        flags += (['-mcpu=apple-m1'] if platform.machine() == 'arm64' else ['-march=x86-64', '-mtune=generic'])
        flags += ['-DSYSTEM=APPLEMAC', '-Xpreprocessor', '-fopenmp', '-I'+prefix+'/include', '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup']
        libraries = [prefix+'/lib/libomp.a']
    else:
        flags += ['-DSYSTEM=STUNIX', '-fopenmp', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections']
        libraries = ['-lm']
    with tempfile.TemporaryDirectory(prefix='ctools-kernel-benchmark-') as directory:
        binary = Path(directory)/'benchmark'
        subprocess.run([os.environ.get('CC', 'clang' if platform.system()=='Darwin' else 'gcc'), *flags, str(ROOT/'validation/benchmark_bigdata_native.c'), *[str(ROOT/'src'/m) for m in modules], *libraries, '-o', str(binary)], check=True)
        for kind in args.workloads:
            for algorithm in (args.algorithms[:1] if kind == 3 else args.algorithms):
                subprocess.run([str(binary), str(args.rows), str(args.threads), str(args.repetitions), str(algorithm), str(kind)], check=True, timeout=600)


if __name__ == '__main__':
    main()
