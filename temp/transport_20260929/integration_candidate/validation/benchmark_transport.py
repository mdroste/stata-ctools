"""Prepare a real-Stata transport benchmark (does not launch Stata).

Builds temporary plugins from NAME=SOURCE_DIRECTORY variants and writes run.do.
Run that file via the machine's stata shell alias. Cases alternate variant
order on each repetition, warm each plugin, verify every cell on an extra pass,
and compare Stata datasignatures before/after. Output: TRANSPORT CSV log records.
macOS toolchain defaults match local benchmark builds; no distribution files change.
"""
import argparse
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant', action='append', required=True, help='NAME=src-directory')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--rows', type=int, default=1000000)
    p.add_argument('--columns', type=int, default=20)
    p.add_argument('--repetitions', type=int, default=5)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--shapes', nargs='+', default=['numeric', 'mixed', 'strings'])
    p.add_argument('--widths', nargs='+', type=int, default=[8, 32, 244])
    p.add_argument('--modes', nargs='+', default=['identity', 'scattered', 'sorted', 'filtered'])
    p.add_argument('--hints', nargs='+', type=int, default=[1, 0])
    p.add_argument('--predicate', default='mod(_n,3)!=0')
    p.add_argument('--prepare-only', action='store_true', help='reuse existing compiled variants')
    args = p.parse_args()
    if min(args.rows, args.columns, args.repetitions, args.threads) < 1:
        p.error('sizes, repetitions, and threads must be positive')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    prefix = Path(os.environ.get('LIBOMP_PREFIX', '/opt/homebrew/opt/libomp'))
    cc = os.environ.get('CC', '/opt/homebrew/opt/llvm/bin/clang')
    variants = []
    commands = []
    for spec in args.variant:
        name, source = spec.split('=', 1)
        source = Path(source).resolve()
        binary = out / (name + '.plugin')
        cmd = [cc, '-O3', '-g', '-fPIC', '-DSYSTEM=APPLEMAC', '-DSD_FASTMODE',
               '-fno-fast-math', '-ffp-contract=off', '-funroll-loops', '-ftree-vectorize',
               '-flto', '-fno-strict-aliasing', '-arch', 'arm64', '-mcpu=apple-m1',
               '-mmacosx-version-min=26.0', '-Xpreprocessor', '-fopenmp',
               '-I'+str(prefix/'include'), '-I'+str(source), '-bundle',
               '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup',
               '-Wl,-exported_symbol,_pginit', '-Wl,-exported_symbol,_stata_call',
               str(ROOT/'validation/benchmark_transport.c'),
               *[str(source/f) for f in ['stplugin.c', 'ctools_data_io.c',
                                         'ctools_types.c', 'ctools_threads.c', 'ctools_arena.c']],
               str(prefix/'lib/libomp.a'), '-o', str(binary)]
        if not args.prepare_only:
            subprocess.run(cmd, check=True)
        commands.append(cmd)
        variants.append((name, binary))
    import json
    (out/'build_commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    lines = ['clear all', 'set more off', 'set linesize 255',
             f'log using "{out}/transport.log", text replace',
             f'adopath ++ "{ROOT}/build"']
    for i, (_, binary) in enumerate(variants):
        lines.append(f'program io{i}, plugin using("{binary}")')
    lines += ['program run_benchmark', 'version 16']
    for shape in args.shapes:
        for width in ([8] if shape == 'numeric' else args.widths):
            lines += ['clear', f'quietly set obs {args.rows}']
            # Keep wide fixtures below Stata's compiled-program size limit.
            # Expanding the padding once per variable can exceed that limit
            # before any benchmark calls execute.
            if shape != 'numeric':
                lines.append(f'local transport_pad "{"x" * width}"')
            for j in range(1, args.columns+1):
                if shape in ('strings', 'long') or (shape == 'mixed' and j % 2 == 0):
                    kind = 'strL' if shape == 'long' else f'str{width}'
                    # Full-width, varying values plus empty strings; UTF-8 tested separately.
                    lines.append(f'quietly gen {kind} v{j} = cond(mod(_n,17)==0,"",substr(string(mod(_n+{j},1000000),"%06.0f")+"`transport_pad\'",1,{width}))')
                else:
                    kind = ['byte', 'int', 'long', 'float', 'double'][(j-1) % 5]
                    expr = f'mod(_n+{j},101)-50' if kind in ('byte', 'int') else f'_n+{j}'
                    if kind in ('float', 'double'):
                        expr = f'({expr})/7'
                    lines += [f'quietly gen {kind} v{j} = {expr}',
                              f'quietly replace v{j} = . if mod(_n,101)==0',
                              f'quietly replace v{j} = .z if mod(_n,103)==0']
            lines += ['quietly datasignature', 'local signature "`r(datasignature)\'"']
            for hints in args.hints:
                lines.append('_ctools_strw v*' if hints else 'local __ctools_strw ""')
                for mode in (['read'] if shape == 'long' else args.modes):
                    case = f'{shape}_w{width}_h{hints}'
                    # Warmup and the extra verification pass are excluded from
                    # reported medians: verification can change cache residency.
                    for rep in range(args.repetitions+2):
                        indices = range(len(variants)) if rep % 2 == 0 else reversed(range(len(variants)))
                        for i in indices:
                            label = f'{variants[i][0]}_{case}_r{rep}'
                            select = ' if ' + args.predicate if mode == 'filtered' else ''
                            lines.append(f'plugin call io{i} v*{select}, "{label}" "{args.threads}" "{mode}" "{int(rep == args.repetitions+1)}"')
                    lines += ['quietly datasignature', 'assert "`r(datasignature)\'" == "`signature\'"']
    lines += ['end', 'capture noisily run_benchmark', 'local rc = _rc',
              'di "TRANSPORT_COMPLETE RC=`rc\'"', 'log close', 'exit, clear']
    draft = out/'run.do.tmp'
    draft.write_text('\n'.join(lines)+'\n')
    draft.replace(out/'run.do')
    print(out/'run.do')


if __name__ == '__main__':
    main()
