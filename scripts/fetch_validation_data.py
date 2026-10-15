"""Prepare pinned Stata example data for offline validation (never run by Stata)."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / 'validation/fixtures'
# examples: official Stata datasets required by the offline gate.
# hdfe: Sergio Correia's HDFE benchmark collection (https://www.scorreia.com/data/hdfe/),
#       about 227 MB, used only by the opt-in validation/run_hdfe_datasets.py.
COLLECTIONS = {
    'examples': {
        'manifest': FIXTURES / 'manifest.json', 'directory': FIXTURES / 'data',
        'base': 'https://www.stata-press.com/data/r18/', 'max_time': '90',
        'names': '''airline auto2 bdesop cancer countxmpl grunfeld lutkepohl2 klein
sysdsn1 lbw nhanes2 nhanes2f nlsw88 nlswork pig bplong ships lifeexp union
uslifeexp2 wpi1 census lifeexp educ99gdp bplong census voter ships'''.split()},
    'hdfe': {
        'manifest': FIXTURES / 'hdfe_manifest.json', 'directory': FIXTURES / 'data/hdfe',
        'base': 'https://f001.backblazeb2.com/file/sergio-public-data/', 'max_time': '900',
        'names': '''credit2 credit soccer synthetic-complete synthetic-uniform-easy
synthetic-uniform-hard synthetic-uniform-harder synthetic-assortative synthetic-zigzag
enron github patents workers schools directors'''.split()},
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collection', choices=sorted(COLLECTIONS), default='examples')
    parser.add_argument('--directory', type=Path, help='cache directory (default: the collection\'s)')
    parser.add_argument('--check-only', action='store_true', help='verify cache without downloading')
    parser.add_argument('--record', action='store_true', help='maintainer: record source checksums')
    args = parser.parse_args()
    collection = COLLECTIONS[args.collection]
    MANIFEST, BASE = collection['manifest'], collection['base']
    if args.directory is None: args.directory = collection['directory']
    args.directory.mkdir(parents=True, exist_ok=True)
    manifest = {} if args.record else json.loads(MANIFEST.read_text())

    def fetch(name):
        target = args.directory / (name + '.dta')
        entry = manifest.get(name, {'url': BASE + name + '.dta'})
        if not target.exists():
            if args.check_only: raise RuntimeError(f'missing fixture: {name}; run scripts/fetch_validation_data.py --collection {args.collection}')
            partial = target.with_suffix('.part')
            result = subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                                     '--retry', '2', '--max-time', collection['max_time'], entry['url'], '-o', str(partial)])
            if result.returncode:
                partial.unlink(missing_ok=True)
                raise RuntimeError(f'download failed: {name}')
            partial.replace(target)
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        if not args.record and digest != entry['sha256']:
            raise RuntimeError(f'checksum mismatch: {name}; remove cached file and retry')
        return name, {'url': entry['url'], 'sha256': digest}

    names = sorted(set(collection['names'])) if args.record else sorted(manifest)
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = {}
        errors = []
        futures = [pool.submit(fetch, name) for name in names]
        for future in as_completed(futures):
            try:
                name, entry = future.result()
                records[name] = entry
            except RuntimeError as error:
                errors.append(str(error))
        if errors:
            raise SystemExit("\n".join(sorted(errors)))
    if args.record:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(json.dumps(records, indent=2, sort_keys=True) + '\n')
    print(f'Verified {len(records)} offline datasets in {args.directory}')

if __name__ == '__main__':
    main()
