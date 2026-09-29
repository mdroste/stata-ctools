"""Prepare or build a single-plugin numeric-scheduler comparison; never launch Stata."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import sys

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO/'validation'))
spec = importlib.util.spec_from_file_location('campaign', REPO/'validation/benchmark_transport_campaign.py')
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)
OUT = ROOT/'numeric_scheduler_toggle'
BASE_SOURCE = ROOT/'integration_candidate/src'
SOURCE = OUT/'src'
REPETITIONS = 11


def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'manifest.json').exists():
        raise ValueError('Experiment already prepared; preserve its provenance')
    shutil.copytree(BASE_SOURCE, SOURCE)
    io = SOURCE/'ctools_data_io.c'
    code = io.read_text()
    gate1 = '(nvars >= IO_NUMERIC_TILE_MIN_VARS && n_filtered >= IO_NUMERIC_TILE_MIN_ROWS &&'
    gate2 = '(string_bytes == 0 && nvars >= IO_NUMERIC_TILE_MIN_VARS && n_filtered >= IO_NUMERIC_TILE_MIN_ROWS &&'
    assert code.count(gate1) == code.count(gate2) == 1
    code = 'extern int ctools_bench_numeric_tiles;\n'+code
    code = code.replace(gate1, '(ctools_bench_numeric_tiles && nvars >= IO_NUMERIC_TILE_MIN_VARS && n_filtered >= IO_NUMERIC_TILE_MIN_ROWS &&')
    code = code.replace(gate2, '(ctools_bench_numeric_tiles && string_bytes == 0 && nvars >= IO_NUMERIC_TILE_MIN_VARS && n_filtered >= IO_NUMERIC_TILE_MIN_ROWS &&')
    io.write_text(code)
    harness = (REPO/'validation/benchmark_transport.c').read_text()
    assert harness.count('if (argc != 4) return 198;') == 1
    harness = harness.replace('#include \"ctools_threads.h\"', '#include \"ctools_threads.h\"\nint ctools_bench_numeric_tiles = 0;')
    harness = harness.replace('if (argc != 4) return 198;', 'if (argc != 4 && argc != 5) return 198;\n    ctools_bench_numeric_tiles = argc == 5 && atoi(argv[4]);')
    harness = harness.replace('        SF_display(line);', '''        SF_display(line);
        snprintf(line, sizeof(line), "TRANSPORT_ORDER_BYTES,%s,%zu\\n", argv[0],
                 data->sort_order ? data->nobs * sizeof(perm_idx_t) : (size_t)0);
        SF_display(line);''')
    needle = '    ctools_filtered_data_free(&fd);\n    return rc;'
    assert harness.count(needle) == 1
    harness = harness.replace(needle, '''    double cleanup_start = now();
    ctools_filtered_data_free(&fd);
    ctools_destroy_global_pool();
    double cleanup = now() - cleanup_start;
    if (!rc) {
        char line[512];
        snprintf(line, sizeof(line), "TRANSPORT_CLEANUP,%s,%.9f\\n", argv[0], cleanup);
        SF_display(line);
    }
    return rc;''')
    (OUT/'benchmark.c').write_text(harness)
    (OUT/'_ctools_strw.ado').write_bytes((REPO/'build/_ctools_strw.ado').read_bytes())
    cases = [dict(name=f'n200000_k{k}_{storage}', rows=200000, columns=k, storage=storage, mode='read')
             for k in (64,128,256) for storage in ('byte','int','long','float','double')]
    cases += [dict(name='n50000_k128_cycle', rows=50000, columns=128, mode='read'),
              dict(name='n10000_k2000_cycle', rows=10000, columns=2000, mode='read')]
    cases = [campaign.validate_case(c) for c in cases]
    runs = []
    for block in range(2):
        directory = OUT/f'b{block}'
        directory.mkdir(exist_ok=True)
        lines = ['clear all', 'set more off', 'set linesize 255',
                 f'log using "{directory}/transport.log", text replace', f'adopath ++ "{OUT}"',
                 f'program transport_io, plugin using("{OUT}/toggle.plugin")',
                 'program run_campaign', 'version 16']
        for case in cases:
            generated = campaign.case_driver(case, 'columns', REPETITIONS).splitlines()
            expanded = []
            for line in generated:
                if line.startswith('plugin call '):
                    # Exactly the same fixture/verification as case_driver. The
                    # two modes share one loaded binary and alternate call order.
                    label = line.split('"')[1]
                    rep = int(label.rsplit('_r', 1)[1])
                    flags = (0, 1) if (rep + block) % 2 == 0 else (1, 0)
                    for flag in flags:
                        call = line if not flag else line.replace('"columns_', '"tiles_', 1)
                        expanded.append(call+f' "{flag}"')
                else:
                    expanded.append(line)
            path = directory/(case['name']+'.do')
            path.write_text('\n'.join(expanded)+'\n')
            lines.append(f'do "{path}"')
        lines += ['end', 'capture noisily run_campaign', 'local rc = _rc',
                  'di "TRANSPORT_COMPLETE RC=`rc\'"', 'log close', 'exit, clear']
        driver = directory/'run.do'
        driver.write_text('\n'.join(lines)+'\n')
        resource = directory/'resources.txt'
        shell = 'oldstata -q -b do '+shlex.quote(str(driver))
        command = shlex.join(['/usr/bin/time', '-l', '/bin/zsh', '-lic', shell])+' 2>'+shlex.quote(str(resource))
        runs.append(dict(name=f'b{block}', block=block, variant='runtime_toggle', cases=[c['name'] for c in cases],
                         driver=str(driver), log=str(directory/'transport.log'), resources=str(resource), command=command))
    hashes = {str(p.relative_to(SOURCE)): campaign.sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file()}
    manifest = dict(schema='single_plugin_numeric_scheduler_v1', base_source=str(BASE_SOURCE), base_io_sha256=campaign.sha(BASE_SOURCE/'ctools_data_io.c'), source=str(SOURCE), source_sha256=hashes,
                    repetitions=REPETITIONS, blocks=2, lifecycle='production', cases=cases, runs=runs,
                    harness_sha256=campaign.sha(OUT/'benchmark.c'), plugin=str(OUT/'toggle.plugin'),
                    requested_runtime_environment={}, compiler='clang', libomp_prefix='/private/tmp/ctools-memory-audit/libomp',
                    compiled=False, summarize_script=str(ROOT/'summarize_numeric_scheduler_toggle.py'))
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (OUT/'run-plan.tsv').write_text('run\tcommand\n'+''.join(r['name']+'\t'+r['command']+'\n' for r in runs))
    print(f'Prepared {len(cases)} cases, one binary, 2 processes, 11 timed paired calls per case/process; NOT compiled')


def build():
    manifest = json.loads((OUT/'manifest.json').read_text())
    current = {str(p.relative_to(SOURCE)): campaign.sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file()}
    if current != manifest['source_sha256']:
        raise ValueError('Frozen source changed after preparation; review and prepare again')
    os.environ['CC'] = manifest['compiler']
    os.environ['LIBOMP_PREFIX'] = manifest['libomp_prefix']
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = '26.0'
    command = campaign.compile_plugin(SOURCE, OUT/'toggle.plugin', OUT/'benchmark.c', [])
    after = {str(p.relative_to(SOURCE)): campaign.sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file()}
    if after != current:
        raise ValueError('Source changed during build')
    manifest.update(compiled=True, command=command, plugin_sha256=campaign.sha(OUT/'toggle.plugin'))
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(OUT/'run-plan.tsv')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', action='store_true')
    args = p.parse_args()
    if args.build:
        build()
    else:
        prepare()
