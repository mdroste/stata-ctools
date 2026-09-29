"""Plan or archive completed transport evidence. Never run builds, Git, or Stata.

Default is plan only. After root declares the experiments/acceptance complete:
  python3 temp/transport_20260929/archive_transport_evidence.py --copy --completed \
    --campaign final_revision --extra path/to/acceptance.log
Frozen manifests are copied verbatim. An additional path map describes relocation.
"""
import argparse
import csv
import difflib
import hashlib
import gzip
import io
import json
from pathlib import Path
import re
import shutil
import tarfile

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
CAMPAIGNS = ['tiny_campaign','tiny_weighted_screen','numeric_screen','str2045_screen',
             'numeric_tile_sweep_full','store_screen','runtime_focus','no_sort_order_focused',
             'final_comparison','final_allocation','no_order_toggle','numeric_scheduler_toggle','map_comparison',
             'numeric_small_byte_toggle','numeric_sixm_toggle']
PROFILES = ['profile_final_baseline','profile_final_candidate']
MODULES = ['stplugin.c','ctools_data_io.c','ctools_types.c','ctools_threads.c','ctools_arena.c']
SUFFIXES = {'.csv','.json','.c','.h','.do','.ado','.log','.txt','.md','.py','.patch','.diff','.tsv'}
ROOT_FILES = ['campaign_index.csv','campaign_raw.csv','campaign_phase_summary.csv','campaign_resources.csv',
              'component_selected.csv','final_selected.csv','profile_leaf_counts.csv','report_draft.md',
              'starting_hashes.json','final_evidence_index.csv','collect_campaign_tables.py','profile_transport.py',
              'prepare_no_order_toggle.py','summarize_no_order_toggle.py',
              'prepare_numeric_scheduler_toggle.py','summarize_numeric_scheduler_toggle.py',
              'prepare_numeric_small_byte_toggle.py','summarize_numeric_small_byte_toggle.py',
              'prepare_numeric_sixm_toggle.py','summarize_numeric_sixm_toggle.py','render_transport_report.py',
              'archive_transport_evidence.py','run_campaign.py','benchmark_adaptive_commands.do',
              'build_checks.json','build_checks_final.json','final_toolchain.json','final_full_builds.json','final_build_Makefile.txt','production_baseline_build.log','production_candidate_build.log',
              'production_final_build.log','production_final_baseline_build.log','prepare_dual_impl.py','run_dual.py','install_final.json',
              'run_commands.py','run_commands_final.py','run_final_timings.py','run_final_acceptance_timings.py']


def sha(data):
    return hashlib.sha256(data).hexdigest()


def transport_sources(source):
    """Five compiled modules plus their recursively quoted local includes."""
    files, pending = {}, [source / name for name in MODULES]
    while pending:
        p = pending.pop().resolve()
        relative = str(p.relative_to(source))
        if relative in files:
            continue
        data = p.read_bytes()
        files[relative] = data
        for name in re.findall(rb'^\s*#\s*include\s*"([^"]+)"', data, re.MULTILINE):
            name = name.decode()
            possibilities = [p.parent/name, source/name]
            match = next((q.resolve() for q in possibilities if q.is_file()), None)
            if match is None:
                raise ValueError(f'Unresolved local include {name} in {p}')
            match.relative_to(source)
            pending.append(match)
    return files


def copy_file(source, target, inventory):
    data = source.read_bytes()
    compress = (source.suffix in {'.log','.txt'} and len(data) > 128*1024) or (source.suffix == '.do' and len(data) > 8*1024)
    archived = gzip.compress(data, mtime=0) if compress else data
    if compress:
        target = target.with_name(target.name+'.gz')
    if target.exists():
        raise ValueError(f'Archive path collision: {target}')
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(archived)
    inventory.append(dict(original=str(source), archived=str(target.relative_to(DEST)),
                          source_sha256=sha(data), archived_sha256=sha(archived),
                          source_bytes=len(data), archived_bytes=len(archived),
                          compression='gzip' if compress else 'none'))


def copy_text_tree(source, dest, inventory):
    for p in sorted(source.rglob('*')):
        if not p.is_file() or p.suffix not in SUFFIXES:
            continue
        # Batch run.log repeats the explicit transport log; retain the latter.
        if p.name == 'run.log' and p.with_name('transport.log').exists():
            continue
        relative = p.relative_to(source)
        if 'src' in relative.parts or any(x.endswith('.dSYM') for x in relative.parts):
            continue
        copy_file(p,dest/relative,inventory)


p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--destination', type=Path, default=REPO/'docs/benchmarks/transport_adaptive_20260929')
p.add_argument('--campaign', action='append', default=[], help='Additional completed campaign directory name')
p.add_argument('--extra', action='append', type=Path, default=[], help='Additional acceptance text/source file or tree')
p.add_argument('--copy', action='store_true')
p.add_argument('--completed', action='store_true', help='Root has declared all selected work complete')
a = p.parse_args()
DEST = a.destination.resolve()
names = list(dict.fromkeys(CAMPAIGNS+a.campaign))
if not a.copy:
    print(json.dumps(dict(destination=str(DEST), campaigns=names, profiles=PROFILES,
                          source_policy='Deduplicated gzip archives of five transport modules plus recursive local headers; diffs against frozen component baseline',
                          large_logs='Lossless .gz for .log/.txt files above 128 KiB; inventory records original-content and archived-byte SHA-256',
                          excluded='Compiled plugins, dSYM bundles, unrelated source trees, unused campaigns',
                          extra=[str(x) for x in a.extra], action='PLAN ONLY; no archive files copied'),indent=2))
    raise SystemExit()
if not a.completed:
    raise SystemExit('Refusing to copy before root declares completion (--completed)')
if DEST.exists():
    raise SystemExit('Destination already exists; choose a fresh destination to preserve the archive')

# Validate everything before creating the durable destination.
manifests = {}
for name in names:
    m = json.loads((ROOT/name/'manifest.json').read_text())
    for run in m['runs']:
        if '\nTRANSPORT_COMPLETE RC=0\n' not in Path(run['log']).read_text():
            raise ValueError(f'Campaign incomplete: {name}, {run["name"]}')
    if not (ROOT/name/'raw.csv').exists() or not (ROOT/name/'summary.csv').exists():
        raise ValueError(f'Campaign has not been summarized: {name}')
    manifests[name] = m
for name in PROFILES:
    m = json.loads((ROOT/name/'profile_metadata.json').read_text())
    if m['exit_code'] != 0 or not m['profiled_timings_excluded_from_benchmark_evidence']:
        raise ValueError(f'Profile not accepted: {name}')

sources = {}
for name,m in manifests.items():
    variants = m.get('variants', [m])
    for v in variants:
        source = Path(v['source']).resolve()
        files = transport_sources(source)
        recorded = v['source_sha256']
        for relative,data in files.items():
            if recorded.get(relative) != sha(data):
                raise ValueError(f'Frozen source no longer matches recorded hash: {source/relative}')
        identity = sha(json.dumps({n:sha(d) for n,d in sorted(files.items())},sort_keys=True).encode())[:20]
        sources[str(source)] = dict(id=identity,files=files,recorded=recorded)
base_source = str((ROOT/'baseline/src').resolve())
base_files = sources[base_source]['files']
# Preserve both the previous evidence candidate and accepted final caller patch.
# Full plugins still require the other repository/build dependencies.
production_pairs = {}
for target_side in ('candidate', 'final'):
    production_payload, production_patch, production_hashes = {}, [], {}
    production_maps = {side: json.loads((ROOT/('production_'+side)/'source_hashes.json').read_text())
                       for side in ('baseline', target_side)}
    changed = [name for name in sorted(production_maps['baseline'].keys() | production_maps[target_side].keys())
               if Path(name).suffix in {'.c','.h','.inc'} and
               production_maps['baseline'].get(name) != production_maps[target_side].get(name)]
    for name in changed:
        pair = {}
        for side in ('baseline', target_side):
            source = ROOT/('production_'+side)/'src'/name
            data = source.read_bytes() if source.exists() else b''
            recorded = production_maps[side].get(name)
            if recorded is not None and recorded != sha(data):
                raise ValueError(f'Production source no longer matches recorded hash: {source}')
            if recorded is not None:
                production_payload[side+'/src/'+name] = data
                production_hashes[side+'/src/'+name] = recorded
            pair[side] = data.decode().splitlines(keepends=True)
        production_patch.extend(difflib.unified_diff(pair['baseline'],pair[target_side],
                                fromfile='a/src/'+name,tofile='b/src/'+name))
    production_pairs[target_side] = dict(payload=production_payload, patch=production_patch,
                                        hashes=production_hashes, changed=changed)
DEST.mkdir(parents=True)
inventory, source_map = [], []
for name in names+PROFILES:
    copy_text_tree(ROOT/name,DEST/name,inventory)
for name in ROOT_FILES:
    source = ROOT/name
    if source.exists():
        copy_file(source,DEST/name,inventory)
for source in [REPO/'validation/benchmark_transport_campaign.py', REPO/'validation/benchmark_transport.c',
               REPO/'validation/benchmark_transport_commands.do', REPO/'validation/benchmark_clock.c',
               *[REPO/'validation'/name for name in ('test_transport_native.py','test_transport_scheduling.py',
                    'test_transport_adaptive_native.py','test_transport_store_native.py')]]:
    copy_file(source,DEST/'tools'/source.name,inventory)
for source in a.extra:
    source=source.resolve()
    try:
        relative_extra = source.relative_to(ROOT)
    except ValueError:
        relative_extra = Path(source.name)
    target=DEST/'acceptance'/relative_extra
    if source.is_dir():
        copy_text_tree(source,target,inventory)
    else:
        if source.suffix not in SUFFIXES:
            raise ValueError(f'Only text/source acceptance files are supported: {source}')
        copy_file(source,target,inventory)
written=set()
for origin,entry in sorted(sources.items()):
    identity=entry['id'];archive=DEST/'sources'/(identity+'.tar.gz')
    archive.parent.mkdir(exist_ok=True)
    if identity not in written:
        with tarfile.open(archive,'w:gz') as tar:
            for name,data in sorted(entry['files'].items()):
                info=tarfile.TarInfo('src/'+name);info.size=len(data);info.mode=0o644;info.mtime=0
                tar.addfile(info,io.BytesIO(data))
        patch=[]
        for name in sorted(base_files.keys() | entry['files'].keys()):
            before=base_files.get(name,b'').decode().splitlines(keepends=True)
            after=entry['files'].get(name,b'').decode().splitlines(keepends=True)
            patch.extend(difflib.unified_diff(before,after,fromfile='a/'+name,tofile='b/'+name))
        (archive.parent/(identity+'.patch')).write_text(''.join(patch))
        (archive.parent/(identity+'.files.json')).write_text(json.dumps({n:sha(d) for n,d in sorted(entry['files'].items())},indent=2)+'\n')
        written.add(identity)
    source_map.append(dict(original_source=origin,bundle=str(archive.relative_to(DEST)),
                           bundle_sha256=sha(archive.read_bytes()),
                           base_bundle='sources/'+sources[base_source]['id']+'.tar.gz',
                           patch='sources/'+identity+'.patch',
                           files_sha256={n:sha(d) for n,d in sorted(entry['files'].items())}))
(DEST/'source_archive_map.json').write_text(json.dumps(source_map,indent=2)+'\n')
for target_side, pair in production_pairs.items():
    stem = 'production_'+target_side+'_changed_sources'
    production_archive = DEST/'sources'/(stem+'.tar.gz')
    with tarfile.open(production_archive,'w:gz') as tar:
        for name,data in sorted(pair['payload'].items()):
            info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644;info.mtime=0
            tar.addfile(info,io.BytesIO(data))
    production_diff = DEST/'sources'/('production_baseline_to_'+target_side+'.patch')
    production_diff.write_text(''.join(pair['patch']))
    (DEST/'sources'/(stem+'.json')).write_text(json.dumps(dict(
        changed_files=pair['changed'], files_sha256=pair['hashes'],
        bundle=str(production_archive.relative_to(DEST)), bundle_sha256=sha(production_archive.read_bytes()),
        patch=str(production_diff.relative_to(DEST)), patch_sha256=sha(production_diff.read_bytes()),
        note='Both sides of changed .c/.h/.inc files; excludes noncompiled .orig backup files. Other source/build dependencies are identified by full source maps.'
    ),indent=2)+'\n')
# Only the human-readable draft needs relative links adjusted for compression;
# all frozen manifests, harnesses, and drivers are copied byte for byte.
draft = DEST/'report_draft.md'
if draft.exists():
    content = draft.read_text()
    for item in inventory:
        try:
            relative = str(Path(item['original']).relative_to(ROOT))
        except ValueError:
            continue
        content = content.replace('('+relative+')','('+item['archived']+')')
    draft.write_text(content)
    record=next(i for i in inventory if i['archived']=='report_draft.md')
    record['archived_sha256']=sha(draft.read_bytes())
    record['archived_bytes']=draft.stat().st_size
    record['transformation']='relative links adjusted to archived artifact paths'
(DEST/'artifact_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
(DEST/'ARCHIVE.md').write_text('''# Frozen transport evidence

Original campaign manifests, source hashes, and recorded compiler commands are unchanged. Logs and profile text larger than 128 KiB are compressed losslessly as `.gz`; Generated do-files larger than 8 KiB are also compressed losslessly; CSV files, manifests, short drivers, and small logs remain directly readable. Decompress archived do-files into the scratch run directory before relocating their recorded paths. Original absolute paths in manifests, run plans, and drivers remain unchanged as historical provenance; rewrite paths only in a scratch reproduction copy after decompression. `artifact_inventory.json` maps original artifacts to archive paths and records both original-content and archived-byte SHA-256 hashes. The human-readable appendix adjusts relative links to relocated or compressed artifacts; its transformation and archived hash are recorded. `source_archive_map.json` maps each original source directory to a deduplicated gzip source archive and a patch against the frozen component baseline.

The hashed transport source archives each contain the five C modules linked by the standalone transport benchmark and their recursively included local headers. Extract into a scratch directory, verify the `.files.json` hashes, and substitute that directory for the original source prefix in the recorded compiler command. Substitute the archived harness path and a scratch output plugin path as well. Compiler and libomp paths are environment-specific and are part of the recorded provenance; comparisons require matching them between variants. The exact runtime-toggle harnesses are archived inside their campaign directories. Every Stata invocation must use the machine's `oldstata` wrapper through interactive login zsh; the archived run plans record it.

These are sufficient sources for the standalone transport plugins, not complete repository snapshots for rebuilding every command. Full recorded source-hash maps remain in the original manifests. Additional merged acceptance artifacts are stored separately under `acceptance/`. `sources/production_candidate_changed_sources.tar.gz` and `sources/production_final_changed_sources.tar.gz` preserve both sides of all changed production source files, including the 11 caller flags, with separate baseline-to-candidate/final patches and per-file hashes. The frozen 64-column gate benchmark and subsequent corrected-gate evidence remain separate; do not overwrite either provenance record.

Compiled plugins and dSYM bundles are intentionally omitted. CPU profile text, including its module-address information, is preserved for interpreting stack samples. Profiled timing records are excluded from transport benchmark summaries. Peak process RSS includes Stata and fixture construction; sample counts are not CPU percentages.
''')
print(f'Archived {len(names)} campaigns, {len(PROFILES)} profiles, {len(written)} distinct transport source bundles to {DEST}')
