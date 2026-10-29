"""Create a versioned, flat ZIP from a complete staged all-platform package."""
import argparse
import json
from pathlib import Path
import re
import sys
import tempfile
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.package_contract import check_package
from scripts.stage_package import PLATFORMS


def archive_package(source, output, release_file=ROOT / 'config/release.json'):
    version = json.loads(release_file.read_text())['version']
    if not isinstance(version, str) or not re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+', version):
        raise ValueError('config/release.json must contain a major.minor.patch version')
    errors = check_package(source, integrity=True)
    if errors:
        raise ValueError('\n'.join(errors))
    info = json.loads((source / 'BUILD_INFO.json').read_text())
    if set(info['platforms']) != set(PLATFORMS):
        raise ValueError('release archives require all four platform plugins')
    if f'd Version: {version}' not in (source / 'ctools.pkg').read_text().splitlines():
        raise ValueError('staged package version differs from config/release.json')
    # Use the verified manifest, excluding caches, logs and other staging debris.
    entries = sorted([*info['sha256'], 'BUILD_INFO.json'])
    output.mkdir(parents=True, exist_ok=True)
    destination = output / f'ctools-{version}.zip'
    with tempfile.TemporaryDirectory(prefix='.ctools-archive-', dir=output) as tmp:
        candidate = Path(tmp) / destination.name
        with ZipFile(candidate, 'w', compression=ZIP_DEFLATED) as archive:
            for entry in entries:
                archive.write(source / entry, arcname=entry)
        candidate.replace(destination)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='directory produced by stage_package.py')
    parser.add_argument('--output', type=Path, default=ROOT / 'dist')
    args = parser.parse_args()
    try:
        result = archive_package(args.source, args.output)
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise SystemExit(str(error))
    print(result)


if __name__ == '__main__':
    main()
