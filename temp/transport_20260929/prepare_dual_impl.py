"""Prepare, build, or summarize a one-plugin transport A/B control.

Preparation never compiles or invokes Stata. Build is a separate explicit mode.
Both data_io translation units use the final shared headers and runtime; every
externally defined data_io function is renamed with a per-TU -D definition.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SHARED_C = ('stplugin.c', 'ctools_types.c', 'ctools_threads.c', 'ctools_arena.c')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def external_functions(source):
    # All definitions in these frozen files have explicit void/stata_retcode
    # return types at column zero; static definitions are deliberately excluded.
    pattern = r'^(?:void|stata_retcode)\s+(ctools_\w+)\s*\([^;{}]*\)\s*\{'
    names = re.findall(pattern, source.read_text(), re.M)
    if len(names) != len(set(names)) or len(names) != 13:
        raise ValueError(f'Review public definition extraction: {source}: {names}')
    return names


def harness_source():
    text = (ROOT/'validation/benchmark_transport.c').read_text()
    text = text.replace('/* Isolated REAL-SPI transport benchmark.',
                        '/* One-plugin, dual-implementation REAL-SPI transport benchmark.')
    text = text.replace('Arguments: label threads mode verify;',
                        'Arguments: label threads mode verify implementation(base|cand);')
    declarations = '''
typedef stata_retcode (*load_fn_t)(ctools_filtered_data *, int *, size_t,
                                  size_t, size_t, int);
typedef stata_retcode (*store_fn_t)(stata_data *, size_t);
typedef void (*lifetime_fn_t)(ctools_filtered_data *);
typedef stata_retcode (*numeric_map_fn_t)(double *, size_t, int, perm_idx_t *);
typedef stata_retcode (*string_map_fn_t)(char **, size_t, int, perm_idx_t *);
#define DECLARE_IO(prefix) \\
extern stata_retcode prefix##ctools_data_load(ctools_filtered_data *, int *, size_t, size_t, size_t, int); \\
extern stata_retcode prefix##ctools_data_store(stata_data *, size_t); \\
extern stata_retcode prefix##ctools_data_store_sorted(stata_data *, size_t); \\
extern void prefix##ctools_filtered_data_init(ctools_filtered_data *); \\
extern void prefix##ctools_filtered_data_free(ctools_filtered_data *); \\
extern stata_retcode prefix##ctools_store_filtered_rowpar(double *, size_t, int, perm_idx_t *); \\
extern stata_retcode prefix##ctools_store_filtered_str(char **, size_t, int, perm_idx_t *);
DECLARE_IO(base_)
DECLARE_IO(cand_)
'''
    text = text.replace('#include "ctools_threads.h"', '#include "ctools_threads.h"\n'+declarations)
    selection = '''    if (argc != 5) return 198;
    int candidate = !strcmp(argv[4], "cand");
    if (!candidate && strcmp(argv[4], "base")) return 198;
    /* All implementation selection occurs before either timed region. The
     * shared thread pool/OpenMP runtime is linked and initialized only once. */
    load_fn_t load_fn = candidate ? cand_ctools_data_load : base_ctools_data_load;
    store_fn_t identity_fn = candidate ? cand_ctools_data_store : base_ctools_data_store;
    store_fn_t sorted_fn = candidate ? cand_ctools_data_store_sorted : base_ctools_data_store_sorted;
    lifetime_fn_t init_fn = candidate ? cand_ctools_filtered_data_init : base_ctools_filtered_data_init;
    lifetime_fn_t free_fn = candidate ? cand_ctools_filtered_data_free : base_ctools_filtered_data_free;
    numeric_map_fn_t numeric_map_fn = candidate ? cand_ctools_store_filtered_rowpar : base_ctools_store_filtered_rowpar;
    string_map_fn_t string_map_fn = candidate ? cand_ctools_store_filtered_str : base_ctools_store_filtered_str;'''
    text = text.replace('    if (argc != 4) return 198;', selection)
    start = text.index('STDLL stata_call')
    prefix, body = text[:start], text[start:]
    replacements = {
        'ctools_filtered_data_init(&fd)': 'init_fn(&fd)',
        'ctools_data_load(&fd': 'load_fn(&fd',
        'ctools_store_filtered_rowpar(v->': 'numeric_map_fn(v->',
        'ctools_store_filtered_str(v->': 'string_map_fn(v->',
        'ctools_data_store_sorted(data, first)': 'sorted_fn(data, first)',
        'ctools_data_store(data, first)': 'identity_fn(data, first)',
    }
    for old, new in replacements.items():
        if old not in body:
            raise ValueError(f'Harness changed: {old}')
        body = body.replace(old, new)
    body = body.replace('    if (rc) return rc;\n    stata_data *data', '''    if (rc) {
        free_fn(&fd);
        ctools_destroy_global_pool();
        return rc;
    }
    stata_data *data''')
    needle = '    ctools_filtered_data_free(&fd);\n    return rc;'
    if body.count(needle) != 1:
        raise ValueError('Harness cleanup changed')
    body = body.replace(needle, '''    double cleanup_start = now();
    free_fn(&fd);
    ctools_destroy_global_pool();
    double cleanup = now() - cleanup_start;
    if (!rc) {
        char line[512];
        snprintf(line, sizeof(line), "TRANSPORT_CLEANUP,%s,%.9f\\n", argv[0], cleanup);
        SF_display(line);
    }
    return rc;''')
    return prefix+body


def campaign_module():
    spec = importlib.util.spec_from_file_location('transport_campaign', ROOT/'validation/benchmark_transport_campaign.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(args):
    out = args.directory.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'manifest.json').exists():
        raise ValueError('Choose a fresh output; existing provenance is immutable')
    sources = {'base': args.baseline.resolve(), 'cand': args.candidate.resolve()}
    shared = out/'shared'
    shared.mkdir()
    # Only the four common C bodies are linked. Freeze all top-level headers so
    # both implementations consume exactly the same ABI and inline helpers.
    for source in sorted(sources['cand'].glob('*.h')):
        shutil.copy2(source, shared/source.name)
    for name in SHARED_C:
        shutil.copy2(sources['cand']/name, shared/name)
    variants = {}
    for name, source in sources.items():
        dest = out/(name+'_ctools_data_io.c')
        shutil.copy2(source/'ctools_data_io.c', dest)
        names = external_functions(dest)
        variants[name] = dict(source=str(source/'ctools_data_io.c'), frozen_source=str(dest),
                              sha256=sha(dest), renamed_symbols={symbol: name+'_'+symbol for symbol in names})
    if set(variants['base']['renamed_symbols']) != set(variants['cand']['renamed_symbols']):
        raise ValueError('Public data_io definitions differ; review namespacing')
    harness = out/'benchmark_dual.c'
    harness.write_text(harness_source())
    shutil.copy2(ROOT/'build/_ctools_strw.ado', out/'_ctools_strw.ado')
    campaign = campaign_module()
    cases = [campaign.validate_case(c) for c in [
        dict(name='mixed_w8', rows=500000, columns=20, shape='mixed', width=8),
        dict(name='mixed_w32', rows=500000, columns=20, shape='mixed', width=32),
        dict(name='mixed_w244', rows=500000, columns=20, shape='mixed', width=244),
        dict(name='numeric_control', rows=1000000, columns=20),
        dict(name='large_sorted_control', rows=1000000, columns=20, mode='sorted'),
        dict(name='wide_byte', rows=200000, columns=64, storage='byte'),
    ]]
    expected = []
    runs = []
    for block in range(args.blocks):
        directory = out/f'b{block}'
        directory.mkdir()
        case_order = cases if block % 2 == 0 else cases[::-1]
        lines = ['clear all', 'set more off', 'set linesize 255',
                 f'log using "{directory}/transport.log", text replace',
                 f'adopath ++ "{out}"', f'program transport_io, plugin using("{out}/dual.plugin")',
                 'program run_campaign', 'version 16']
        for case in case_order:
            # Reuse exactly the fixtures/hints from the established campaign.
            original = campaign.case_driver(case, 'unused', 0).splitlines()
            first_call = next(i for i, line in enumerate(original) if line.startswith('plugin call '))
            content = original[:first_call]
            variables = ' '.join(f'v{i}' for i in range(1, case['columns']+1))
            for rep in range(args.repetitions+2):
                order = ['base', 'cand'] if (rep+block) % 2 == 0 else ['cand', 'base']
                for position, name in enumerate(order):
                    verify = int(rep == args.repetitions+1)
                    label = f'{name}_{case["name"]}_b{block}_r{rep}'
                    content.append(f'plugin call transport_io {variables}, "{label}" "{case["threads"]}" "{case["mode"]}" "{verify}" "{name}"')
                    # Sorted writeback is restored by the C harness outside
                    # timing. Check it before the next implementation sees it.
                    if case['mode'] == 'sorted' or verify:
                        content += ['quietly datasignature', 'assert "`r(datasignature)\'" == "`signature\'"']
                    expected.append(dict(label=label, block=block, case=case['name'], variant=name,
                                         rep=rep, position=position, verify=verify,
                                         expected_rows=case['rows'], expected_columns=case['columns'],
                                         expected_mode=case['mode'], expected_threads=case['threads'],
                                         timed=1 <= rep <= args.repetitions))
            content += ['quietly datasignature', 'assert "`r(datasignature)\'" == "`signature\'"']
            driver = directory/(case['name']+'.do')
            driver.write_text('\n'.join(content)+'\n')
            lines.append(f'do "{driver}"')
        lines += ['end', 'capture noisily run_campaign', 'local rc = _rc',
                  'di "DUAL_TRANSPORT_COMPLETE RC=`rc\'"', 'log close', 'exit, clear']
        driver = directory/'run.do'
        driver.write_text('\n'.join(lines)+'\n')
        command = shlex.join(['/bin/zsh', '-lic', 'oldstata -q -b do '+shlex.quote(str(driver))])
        runs.append(dict(block=block, driver=str(driver), log=str(directory/'transport.log'), command=command))
    prefix = Path(os.environ.get('LIBOMP_PREFIX', '/private/tmp/ctools-memory-audit/libomp'))
    cc = os.environ.get('CC', 'clang')
    flags = [cc, '-O3', '-g', '-fPIC', '-DSYSTEM=APPLEMAC', '-DSD_FASTMODE', '-fno-fast-math',
             '-ffp-contract=off', '-funroll-loops', '-ftree-vectorize', '-flto',
             '-fno-strict-aliasing', '-arch', 'arm64', '-mcpu=apple-m1',
             '-mmacosx-version-min='+os.environ.get('MACOSX_DEPLOYMENT_TARGET', '26.0'),
             '-Xpreprocessor', '-fopenmp', '-I'+str(prefix/'include'), '-I'+str(shared)]
    commands = []
    objects = []
    for name, info in variants.items():
        obj = out/(name+'_ctools_data_io.o')
        defines = ['-D'+symbol+'='+renamed for symbol, renamed in info['renamed_symbols'].items()]
        commands.append(flags+defines+['-c', info['frozen_source'], '-o', str(obj)])
        objects.append(str(obj))
    for source in [harness, *[shared/name for name in SHARED_C]]:
        obj = out/(source.stem+'.o')
        commands.append(flags+['-c', str(source), '-o', str(obj)])
        objects.append(str(obj))
    binary = out/'dual.plugin'
    commands.append(flags+['-bundle', '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup',
                           '-Wl,-exported_symbol,_pginit', '-Wl,-exported_symbol,_stata_call',
                           *objects, str(prefix/'lib/libomp.a'), '-o', str(binary)])
    # Force shared headers/runtimes and source snapshots to remain fixed between
    # prepare and build; no shell substitutions occur in any compiler command.
    manifest = dict(schema=1, lifecycle='production', repetitions=args.repetitions, blocks=args.blocks,
                    variants=variants, cases=cases, runs=runs, expected=expected,
                    harness_sha256=sha(harness), preparation_script_sha256=sha(Path(__file__)),
                    shared_sha256={p.name: sha(p) for p in sorted(shared.iterdir())},
                    helper_sha256=sha(out/'_ctools_strw.ado'), commands=commands,
                    compiler=cc, libomp=str(prefix/'lib/libomp.a'),
                    runtime_environment={k: os.environ.get(k) for k in ('OMP_WAIT_POLICY', 'KMP_BLOCKTIME', 'OMP_NUM_THREADS')})
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (out/'run-plan.tsv').write_text('block\tcommand\n'+''.join(f'{r["block"]}\t{r["command"]}\n' for r in runs))
    (out/'README.md').write_text('''# Same-plugin A/B causal control

Preparation freezes baseline and candidate data_io bodies as separate translation
units. All 13 externally defined functions receive base_ or cand_ prefixes.
Both use the final headers, types, arena, thread pool, SPI wrapper, and one OpenMP
runtime. The runtime selector chooses function pointers before load/store timing.

Each block uses one Stata process and the same plugin. Pairs alternate AB/BA;
block 2 reverses both the initial order and case order. Each fixture has one
warmup pair, 11 measured pairs, and one all-cell verification pair. Sorted stores
restore the original fixture outside timing, and Stata checks its datasignature
before the next call. All other cases check signatures after each verification
call and at case completion. Cleanup/free/pool destruction is timed separately.

No load-time environment override is introduced. Source namespacing changes code
placement, so this is not an identical-body flag-toggle control, and shared final
inline helpers mean this isolates the data_io bodies rather than every baseline
source change. Existing process-isolated comparisons remain separate evidence.

Compile only when other timing has stopped:
`python3 temp/transport_20260929/prepare_dual_impl.py build DIRECTORY`

Run the exact oldstata commands in run-plan.tsv sequentially, then:
`python3 temp/transport_20260929/prepare_dual_impl.py summarize DIRECTORY`
''')
    print(f'Prepared {len(cases)} cases, {args.blocks} blocks, {len(expected)} calls; no compilation or Stata execution.')
    print(out/'manifest.json')


def build(args):
    out = args.directory.resolve()
    manifest = json.loads((out/'manifest.json').read_text())
    for info in manifest['variants'].values():
        assert sha(Path(info['frozen_source'])) == info['sha256']
    for name, expected in manifest['shared_sha256'].items():
        assert sha(out/'shared'/name) == expected
    assert sha(out/'benchmark_dual.c') == manifest['harness_sha256']
    with (out/'build.log').open('w') as log:
        for command in manifest['commands']:
            log.write(shlex.join(command)+'\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    # A full Mach-O audit also catches forgotten public data_io definitions or
    # an accidental unresolved unprefixed reference after final LTO/linking.
    result = subprocess.run(['/usr/bin/nm', '-m', str(out/'dual.plugin')], capture_output=True, text=True, check=True)
    (out/'symbols.txt').write_text(result.stdout)
    for name in manifest['variants']['base']['renamed_symbols']:
        if re.search(r'\s_'+re.escape(name)+r'(?:\s|$)', result.stdout):
            raise ValueError(f'Unexpected unprefixed data_io symbol: {name}')
    (out/'build_result.json').write_text(json.dumps(dict(exit_code=0, binary_sha256=sha(out/'dual.plugin'),
        libomp_sha256=sha(Path(manifest['libomp'])), symbol_audit='no unprefixed data_io symbols'), indent=2)+'\n')
    print(out/'dual.plugin')


def summarize(args):
    out = args.directory.resolve()
    manifest = json.loads((out/'manifest.json').read_text())
    observations = {}
    cleanup = {}
    for run in manifest['runs']:
        text = Path(run['log']).read_text()
        if '\nDUAL_TRANSPORT_COMPLETE RC=0\n' not in text:
            raise ValueError(f'Incomplete or failing run: {run["log"]}')
        for line in text.splitlines():
            if line.startswith('TRANSPORT,'):
                _, label, mode, n, k, threads, load, store, verify = line.split(',')
                if label in observations:
                    raise ValueError(f'Duplicate observation: {label}')
                observations[label] = dict(load=float(load), store=float(store), transport=float(load)+float(store),
                                           rows=int(n), columns=int(k), threads=int(threads), mode=mode, verify=int(verify))
            elif line.startswith('TRANSPORT_CLEANUP,'):
                _, label, value = line.split(',')
                if label in cleanup:
                    raise ValueError(f'Duplicate cleanup: {label}')
                cleanup[label] = float(value)
    labels = {item['label'] for item in manifest['expected']}
    if labels != set(observations) or labels != set(cleanup):
        raise ValueError(f'Unexpected/missing rows: {labels ^ set(observations)}, cleanup: {labels ^ set(cleanup)}')
    rows = []
    for expected in manifest['expected']:
        row = observations[expected['label']]
        if row['verify'] != expected['verify']:
            raise ValueError(f'Wrong verify flag: {expected["label"]}')
        for dimension in ('rows', 'columns', 'mode', 'threads'):
            if row[dimension] != expected['expected_'+dimension]:
                raise ValueError(f'Wrong {dimension}: {expected["label"]}')
        rows.append(expected | row | dict(cleanup=cleanup[expected['label']],
                     total=row['transport']+cleanup[expected['label']]))
    with (out/'raw.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    summaries = []
    for case in manifest['cases']:
        for metric in ('load', 'store', 'transport', 'cleanup', 'total'):
            selected = [row for row in rows if row['case'] == case['name'] and row['timed']]
            pairs = []
            block_ratios = []
            for block in range(manifest['blocks']):
                group = [row for row in selected if row['block'] == block]
                values = {name: [row[metric] for row in group if row['variant'] == name] for name in ('base', 'cand')}
                block_ratios.append(statistics.median(values['base'])/statistics.median(values['cand']))
                for rep in range(1, manifest['repetitions']+1):
                    pair = {row['variant']: row[metric] for row in group if row['rep'] == rep}
                    pairs.append(pair['base']/pair['cand'])
            summaries.append(dict(case=case['name'], metric=metric,
                baseline_median=statistics.median(row[metric] for row in selected if row['variant']=='base'),
                candidate_median=statistics.median(row[metric] for row in selected if row['variant']=='cand'),
                paired_speedup_median=statistics.median(pairs), paired_speedup_min=min(pairs), paired_speedup_max=max(pairs),
                block_speedups=';'.join(f'{value:.6f}' for value in block_ratios)))
    with (out/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader(); writer.writerows(summaries)
    print(out/'summary.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'build', 'summarize'))
    parser.add_argument('directory', type=Path)
    parser.add_argument('--blocks', type=int, default=2)
    parser.add_argument('--repetitions', type=int, default=11)
    parser.add_argument('--baseline', type=Path, default=HERE/'production_baseline/src')
    parser.add_argument('--candidate', type=Path, default=HERE/'production_final/src')
    args = parser.parse_args()
    if min(args.blocks, args.repetitions) < 1:
        parser.error('blocks and repetitions must be positive')
    globals()[args.action](args)
