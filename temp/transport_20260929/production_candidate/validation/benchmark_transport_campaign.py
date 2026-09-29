"""Prepare and summarize isolated real-Stata transport comparisons.

This program NEVER launches Stata. It writes frozen benchmark plugins, do-files,
a manifest, and a run-plan.tsv of exact stata commands for sequential execution.
Each process loads one plugin only. Blocks alternate variant order (AB, BA), and
summaries first reduce repetitions to process medians to avoid pseudoreplication.

Example:
  python3 validation/benchmark_transport_campaign.py prepare \
    --variant baseline=/tmp/baseline/src --variant candidate=/tmp/candidate/src \
    --cases /tmp/cases.json --output /tmp/campaign --blocks 4 --repetitions 7
  python3 validation/benchmark_transport_campaign.py summarize /tmp/campaign

Cases are JSON objects with name, rows, columns, and optional shape (numeric,
mixed, strings, long), width, actual_width, storage (cycle/byte/int/long/float/
double), threads, hints, mode, predicate, host_columns, layout (front/spread/
reverse), and obs_range ([first,last]). Use --write-example to inspect a pilot.
Peak RSS from run-plan resource files includes Stata, fixtures, and verification;
use --isolate-cases to attribute it to one case. It is not plugin-only RSS.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[1]
NUMERIC = ['byte', 'int', 'long', 'float', 'double']
DEFAULTS = dict(shape='numeric', width=32, storage='cycle', threads=8, hints=1,
                mode='identity', predicate='mod(_n,3)!=0', layout='front')
EXAMPLE = [
    dict(name='tiny_numeric', rows=32, columns=20),
    dict(name='wide_numeric', rows=100000, columns=256),
    dict(name='wide_mixed', rows=100000, columns=256, shape='mixed'),
    dict(name='subset_spread', rows=250000, columns=8, host_columns=256,
         layout='spread', shape='numeric'),
    dict(name='sparse_mixed', rows=1000000, columns=20, shape='mixed',
         mode='filtered', predicate='mod(_n,100)==0'),
    dict(name='empty_mixed', rows=100000, columns=20, shape='mixed',
         mode='filtered', predicate='0'),
    dict(name='strl', rows=200000, columns=2, shape='long', mode='read'),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_name(value):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', value):
        raise ValueError(f'Use letters, digits, and underscores in name: {value!r}')
    return value


def validate_case(case):
    c = DEFAULTS | case
    safe_name(c['name'])
    for key in ('rows', 'columns', 'threads', 'width'):
        if not isinstance(c[key], int) or c[key] < 1:
            raise ValueError(f'{c["name"]}: {key} must be positive integer')
    if c['rows'] > 2147483647 or c['width'] > 2045:
        raise ValueError('Rows and strings must fit the supported SPI range')
    c.setdefault('host_columns', c['columns'])
    c.setdefault('actual_width', c['width'])
    if c['host_columns'] < c['columns'] or not 0 <= c['actual_width'] <= c['width']:
        raise ValueError('host_columns >= columns; 0 <= actual_width <= width required')
    if c['shape'] not in ('numeric', 'mixed', 'strings', 'long'):
        raise ValueError(f'Unknown shape {c["shape"]}')
    if c['storage'] not in ['cycle', *NUMERIC] or c['layout'] not in ('front', 'spread', 'reverse'):
        raise ValueError('Invalid storage or layout')
    if c['mode'] not in ('identity', 'checked', 'sorted', 'scattered', 'filtered', 'filteredread', 'read'):
        raise ValueError(f'Unknown mode {c["mode"]}')
    if c['shape'] == 'long' and c['mode'] not in ('read', 'filteredread'):
        raise ValueError('strL cases must use read or filteredread')
    if c['hints'] not in (0, 1) or any(ch in c['predicate'] for ch in '\n\r'):
        raise ValueError('Invalid hints or predicate')
    if 'obs_range' in c:
        first, last = c['obs_range']
        if not 1 <= first <= last <= c['rows']:
            raise ValueError('obs_range must be inside the dataset')
    return c


def selected_indices(c):
    k, host = c['columns'], c['host_columns']
    if c['layout'] == 'spread':
        return [1 + j * (host - 1) // max(k - 1, 1) for j in range(k)]
    indices = list(range(1, k + 1))
    return indices[::-1] if c['layout'] == 'reverse' else indices


def case_driver(c, variant, repetitions):
    indices = selected_indices(c)
    selected = {idx: j for j, idx in enumerate(indices, 1)}
    varlist = ' '.join(f'v{idx}' for idx in indices)
    lines = ['clear', f'quietly set obs {c["rows"]}',
             f'local transport_pad "{"x" * c["actual_width"]}"']
    for idx in range(1, c['host_columns'] + 1):
        # Shape describes selected columns; unselected host columns are doubles.
        j = selected.get(idx, 0)
        string = j and (c['shape'] in ('strings', 'long') or (c['shape'] == 'mixed' and j % 2 == 0))
        if string:
            kind = 'strL' if c['shape'] == 'long' else f'str{c["width"]}'
            expr = f'cond(mod(_n,17)==0,"",substr(string(mod(_n+{j},1000000),"%06.0f")+"`transport_pad\'",1,{c["actual_width"]}))'
        else:
            kind = (NUMERIC[(j - 1) % 5] if c['storage'] == 'cycle' else c['storage']) if j else 'double'
            expr = f'mod(_n+{idx},101)-50' if kind in ('byte', 'int') else f'_n+{idx}'
            if kind in ('float', 'double'):
                expr = f'({expr})/7'
        lines.append(f'quietly gen {kind} v{idx} = {expr}')
        if not string and j:
            lines += [f'quietly replace v{idx} = . if mod(_n,101)==0',
                      f'quietly replace v{idx} = .z if mod(_n,103)==0']
    lines += ['quietly datasignature', 'local signature "`r(datasignature)\'"',
              f'_ctools_strw {varlist}' if c['hints'] else 'local __ctools_strw ""']
    qualifier = (' if ' + c['predicate']) if c['mode'] in ('filtered', 'filteredread') else ''
    if 'obs_range' in c:
        qualifier += f' in {c["obs_range"][0]}/{c["obs_range"][1]}'
    for rep in range(repetitions + 2):
        label = f'{variant}_{c["name"]}_r{rep}'
        lines.append(f'plugin call transport_io {varlist}{qualifier}, "{label}" "{c["threads"]}" "{c["mode"]}" "{int(rep == repetitions + 1)}"')
    lines += ['quietly datasignature', 'assert "`r(datasignature)\'" == "`signature\'"']
    return '\n'.join(lines) + '\n'


def compile_plugin(source, binary, harness, flags):
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        raise ValueError('This campaign builder currently requires macOS arm64; port the toolchain explicitly')
    prefix = Path(os.environ.get('LIBOMP_PREFIX', '/opt/homebrew/opt/libomp'))
    cmd = [os.environ.get('CC', '/opt/homebrew/opt/llvm/bin/clang'), '-O3', '-g',
           '-fPIC', '-DSYSTEM=APPLEMAC', '-DSD_FASTMODE', '-fno-fast-math',
           '-ffp-contract=off', '-funroll-loops', '-ftree-vectorize', '-flto',
           '-fno-strict-aliasing', '-arch', 'arm64', '-mcpu=apple-m1',
           '-mmacosx-version-min='+os.environ.get('MACOSX_DEPLOYMENT_TARGET', '26.0'),
           '-Xpreprocessor', '-fopenmp',
           '-I'+str(prefix/'include'), '-I'+str(source), '-bundle',
           '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup',
           '-Wl,-exported_symbol,_pginit', '-Wl,-exported_symbol,_stata_call',
           *flags, str(harness),
           *[str(source/f) for f in ('stplugin.c', 'ctools_data_io.c', 'ctools_types.c',
                                     'ctools_threads.c', 'ctools_arena.c')],
           str(prefix/'lib/libomp.a'), '-o', str(binary)]
    subprocess.run(cmd, check=True)
    return cmd


def prepare(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'manifest.json').exists():
        raise ValueError('Campaign already exists; choose a fresh output to preserve its provenance')
    if min(args.blocks, args.repetitions) < 1:
        raise ValueError('blocks and repetitions must be positive')
    cases = [validate_case(c) for c in json.loads(args.cases.read_text())]
    if args.only:
        wanted = set(args.only)
        cases = [c for c in cases if c['name'] in wanted]
        if wanted != {c['name'] for c in cases}:
            raise ValueError('Requested case missing from --cases')
    if not cases or len({c['name'] for c in cases}) != len(cases):
        raise ValueError('Case names must be nonempty and unique')
    if args.no_sort_order and any(c['mode'] in ('sorted', 'scattered') for c in cases):
        raise ValueError('--no-sort-order is incompatible with sorted/scattered cases')
    # Freeze the benchmark harness as well as production sources. Cleanup is a
    # separate timing because the original harness excludes free/pool teardown.
    harness = out/'benchmark_transport_campaign.c'
    original = (ROOT/'validation/benchmark_transport.c').read_text()
    if args.no_sort_order:
        load_flags = 'check_if ? CTOOLS_LOAD_CHECK_IF : CTOOLS_LOAD_SKIP_IF);'
        if original.count(load_flags) != 1:
            raise ValueError('Benchmark load flags changed; review opt-out instrumentation')
        original = original.replace(load_flags, '(check_if ? CTOOLS_LOAD_CHECK_IF : CTOOLS_LOAD_SKIP_IF) | CTOOLS_LOAD_NO_SORT_ORDER);')
        original = original.replace('#include "ctools_threads.h"', '#include "ctools_threads.h"\n#ifndef CTOOLS_LOAD_NO_SORT_ORDER\n#define CTOOLS_LOAD_NO_SORT_ORDER 0x02\n#endif')
        original = original.replace('        SF_display(line);', '''        SF_display(line);
        snprintf(line, sizeof(line), "TRANSPORT_ORDER_BYTES,%s,%zu\\n", argv[0],
                 data->sort_order ? data->nobs * sizeof(perm_idx_t) : (size_t)0);
        SF_display(line);''')
    needle = '    ctools_filtered_data_free(&fd);\n    return rc;'
    if original.count(needle) != 1:
        raise ValueError('Benchmark harness cleanup changed; review instrumentation')
    replacement = '''    double cleanup_start = now();
    ctools_filtered_data_free(&fd);
    POOL_CLEANUP
    double cleanup = now() - cleanup_start;
    if (!rc) {
        char cleanup_line[512];
        snprintf(cleanup_line, sizeof(cleanup_line), "TRANSPORT_CLEANUP,%s,%.9f\\n", argv[0], cleanup);
        SF_display(cleanup_line);
    }
    return rc;'''.replace('POOL_CLEANUP', 'ctools_destroy_global_pool();' if args.lifecycle == 'production' else '/* Keep the pool for kernel-only measurements. */')
    harness.write_text(original.replace(needle, replacement))
    flags = defaultdict(list)
    environment = {}
    for spec in args.env:
        key, value = spec.split('=', 1)
        if not re.fullmatch(r'(?:OMP|KMP)_[A-Z0-9_]+', key):
            raise ValueError('--env is limited to explicit OMP_* or KMP_* runtime settings')
        environment[key] = value
    for spec in args.cflag:
        name, flag = spec.split('=', 1)
        flags[name].append(flag)
    variants = []
    for spec in args.variant:
        name, source = spec.split('=', 1)
        safe_name(name)
        if name in [v['name'] for v in variants]:
            raise ValueError('Variant names must be unique')
        source = Path(source).resolve()
        binary = out/(name+'.plugin')
        before = {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*')) if p.is_file()}
        command = compile_plugin(source, binary, harness, flags[name])
        after = {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*')) if p.is_file()}
        if before != after:
            raise ValueError(f'{name}: source changed during build; freeze the snapshot first')
        variants.append(dict(name=name, source=str(source), source_sha256=before,
                             plugin=str(binary), plugin_sha256=sha(binary), command=command))
    if set(flags) - {v['name'] for v in variants}:
        raise ValueError('--cflag references an unknown variant')
    helper = ROOT/'build/_ctools_strw.ado'
    (out/'_ctools_strw.ado').write_bytes(helper.read_bytes())
    runs = []
    for block in range(args.blocks):
        ordered = variants if block % 2 == 0 else variants[::-1]
        groups = [[c] for c in cases] if args.isolate_cases else [cases]
        for group in groups:
            for v in ordered:
                name = f'b{block}_{v["name"]}' + ('_'+group[0]['name'] if args.isolate_cases else '')
                directory = out/name
                directory.mkdir()
                lines = ['clear all', 'set more off', 'set linesize 255',
                         f'log using "{directory}/transport.log", text replace',
                         f'adopath ++ "{out}"', f'program transport_io, plugin using("{v["plugin"]}")',
                         'program run_campaign', 'version 16']
                for c in group:
                    casefile = directory/(c['name']+'.do')
                    casefile.write_text(case_driver(c, v['name'], args.repetitions))
                    lines.append(f'do "{casefile}"')
                lines += ['end', 'capture noisily run_campaign', 'local rc = _rc',
                          'di "TRANSPORT_COMPLETE RC=`rc\'"', 'log close', 'exit, clear']
                driver = directory/'run.do'
                driver.write_text('\n'.join(lines)+'\n')
                shell = 'stata -q -b do '+shlex.quote(str(driver))
                resource = directory/'resources.txt'
                prefix = ['/usr/bin/time', '-l']
                if environment:
                    prefix += ['env', *[key+'='+value for key, value in environment.items()]]
                command = shlex.join([*prefix, '/bin/zsh', '-lic', shell])+' 2>'+shlex.quote(str(resource))
                runs.append(dict(name=name, block=block, variant=v['name'], cases=[c['name'] for c in group],
                                 driver=str(driver), log=str(directory/'transport.log'),
                                 resources=str(resource), command=command))
    manifest = dict(schema=1, platform=platform.platform(), lifecycle=args.lifecycle,
                    no_sort_order=args.no_sort_order,
                    blocks=args.blocks, repetitions=args.repetitions, isolate_cases=args.isolate_cases,
                    preparation_environment={k: os.environ.get(k) for k in ('OMP_WAIT_POLICY', 'KMP_BLOCKTIME', 'OMP_NUM_THREADS')},
                    requested_runtime_environment=environment,
                    harness_sha256=sha(harness), width_helper_sha256=sha(helper),
                    cases=cases, variants=variants, runs=runs)
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    with (out/'run-plan.tsv').open('w') as f:
        f.write('run\tcommand\n')
        for r in runs:
            f.write(r['name']+'\t'+r['command']+'\n')
    print(f'{len(cases)} cases, {len(variants)} variants, {len(runs)} isolated Stata processes')
    print(out/'run-plan.tsv')


def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(args):
    out = args.directory.resolve()
    m = json.loads((out/'manifest.json').read_text())
    raw, process, resources, incomplete = [], [], [], []
    for run in m['runs']:
        path = Path(run['log'])
        text = path.read_text() if path.exists() else ''
        if '\nTRANSPORT_COMPLETE RC=0\n' not in text:
            incomplete.append(run['name'])
            continue
        timings, cleanups, order_bytes = {}, {}, {}
        for line in text.splitlines():
            if line.startswith('TRANSPORT,'):
                _, label, mode, n, k, threads, load, store, verify = line.split(',')
                if label in timings:
                    raise ValueError(f'Duplicate observation {run["name"]}: {label}')
                timings[label] = dict(mode=mode, selected_rows=int(n), columns=int(k), threads=int(threads),
                                      load=float(load), store=float(store), verify=int(verify))
            elif line.startswith('TRANSPORT_CLEANUP,'):
                _, label, cleanup = line.split(',')
                if label in cleanups:
                    raise ValueError(f'Duplicate cleanup observation {run["name"]}: {label}')
                cleanups[label] = float(cleanup)
            elif line.startswith('TRANSPORT_ORDER_BYTES,'):
                _, label, size = line.split(',')
                if label in order_bytes:
                    raise ValueError(f'Duplicate allocation observation {run["name"]}: {label}')
                order_bytes[label] = int(size)
        expected = {f'{run["variant"]}_{name}_r{rep}' for name in run['cases'] for rep in range(m['repetitions']+2)}
        if set(timings) != expected or set(cleanups) != expected:
            raise ValueError(f'{run["name"]}: missing or unexpected records')
        if m.get('no_sort_order') and set(order_bytes) != expected:
            raise ValueError(f'{run["name"]}: missing sort-order allocation records')
        resource_path = Path(run['resources'])
        resource_text = resource_path.read_text() if resource_path.exists() else ''
        resource = dict(run=run['name'], block=run['block'], variant=run['variant'])
        for key, label in [('max_rss_bytes', 'maximum resident set size'), ('swaps', 'swaps'),
                           ('page_reclaims', 'page reclaims'), ('page_faults', 'page faults')]:
            found = re.search(r'^\s*(\d+)\s+'+label+r'\s*$', resource_text, flags=re.MULTILINE)
            resource[key] = int(found.group(1)) if found else ''
        resources.append(resource)
        for name in run['cases']:
            values = []
            for rep in range(m['repetitions'] + 2):
                label = f'{run["variant"]}_{name}_r{rep}'
                item = timings[label]
                if item['verify'] != int(rep == m['repetitions']+1):
                    raise ValueError(f'{label}: inconsistent verification marker')
                row = dict(run=run['name'], block=run['block'], variant=run['variant'], case=name,
                           repetition=rep, **item, cleanup=cleanups[label],
                           sort_order_bytes=order_bytes.get(label, ''))
                row['transport'] = row['load'] + row['store'] + row['cleanup']
                raw.append(row)
                if rep and not item['verify']:
                    values.append(row)
            process.append(dict(run=run['name'], block=run['block'], variant=run['variant'], case=name,
                                **{phase: statistics.median(v[phase] for v in values)
                                   for phase in ('load', 'store', 'cleanup', 'transport')}))
    if incomplete and not args.allow_incomplete:
        raise ValueError('Incomplete runs (use --allow-incomplete to inspect partial evidence): '+', '.join(incomplete))
    write_csv(out/'raw.csv', raw)
    write_csv(out/'process.csv', process)
    write_csv(out/'resources.csv', resources)
    summary = []
    baseline = args.baseline or m['variants'][0]['name']
    grouped = defaultdict(dict)
    for row in process:
        grouped[(row['case'], row['variant'])][row['block']] = row
    for c in m['cases']:
        bases = grouped.get((c['name'], baseline), {})
        for variant in [v['name'] for v in m['variants'] if v['name'] != baseline]:
            candidates = grouped.get((c['name'], variant), {})
            blocks = sorted(bases.keys() & candidates.keys())
            if not blocks:
                continue
            for phase in ('load', 'store', 'cleanup', 'transport'):
                b = [bases[i][phase] for i in blocks]
                n = [candidates[i][phase] for i in blocks]
                ratios = [x/y for x, y in zip(b, n) if y > 0]
                summary.append(dict(case=c['name'], baseline=baseline, candidate=variant, phase=phase,
                                    process_pairs=len(blocks), baseline_median=statistics.median(b),
                                    candidate_median=statistics.median(n),
                                    paired_speedup_median=statistics.median(ratios) if ratios else '',
                                    paired_speedup_min=min(ratios) if ratios else '',
                                    paired_speedup_max=max(ratios) if ratios else '',
                                    baseline_min=min(b), baseline_max=max(b), candidate_min=min(n), candidate_max=max(n)))
    write_csv(out/'summary.csv', summary)
    print(f'{len(raw)} records; {len(process)} process medians; {len(summary)} paired phase summaries')
    if incomplete:
        print('Incomplete: '+', '.join(incomplete))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-example', type=Path)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    prep.add_argument('--variant', action='append', required=True)
    prep.add_argument('--cflag', action='append', default=[], help='NAME=compiler-flag, repeatable')
    prep.add_argument('--env', action='append', default=[], help='OMP_*/KMP_*=value runtime override, repeatable')
    prep.add_argument('--cases', type=Path, required=True)
    prep.add_argument('--output', type=Path, required=True)
    prep.add_argument('--blocks', type=int, default=4)
    prep.add_argument('--repetitions', type=int, default=7)
    prep.add_argument('--lifecycle', choices=['production', 'persistent'], default='production')
    prep.add_argument('--no-sort-order', action='store_true', help='Pass the optional read-only sort-order allocation opt-out flag')
    prep.add_argument('--isolate-cases', action='store_true')
    prep.add_argument('--only', nargs='+')
    summary = sub.add_parser('summarize')
    summary.add_argument('directory', type=Path)
    summary.add_argument('--baseline')
    summary.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    if args.write_example:
        args.write_example.write_text(json.dumps(EXAMPLE, indent=2)+'\n')
    elif args.action == 'prepare':
        prepare(args)
    elif args.action == 'summarize':
        summarize(args)
    else:
        parser.error('Choose prepare, summarize, or --write-example')


if __name__ == '__main__':
    main()
