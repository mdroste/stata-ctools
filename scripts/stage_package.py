"""Stage a complete Stata net-install directory for one platform or all platforms."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'validation'))
from check_package import check_package

PLATFORMS = ('mac_arm', 'mac_x86', 'windows', 'linux')

def native_platform():
    system, machine = platform.system(), platform.machine().lower()
    if system == 'Darwin': return 'mac_arm' if machine in ('arm64', 'aarch64') else 'mac_x86'
    if system == 'Windows': return 'windows'
    if system == 'Linux' and machine in ('x86_64', 'amd64'): return 'linux'
    raise ValueError(f'unsupported platform: {system} {machine}')


def stage(source, output, selected, revision):
    source, output = source.resolve(), output.resolve()
    if output == source or source in output.parents:
        raise ValueError('stage outside the source build directory')
    if output.exists():
        raise ValueError(f'output already exists; choose a fresh staging directory: {output}')
    lines = (source / 'ctools.pkg').read_text().splitlines()
    keep = {f'ctools_{p}.plugin' for p in selected}
    staged_lines = [line for line in lines if not line.startswith('d Platforms: ')
                    and not (line.startswith('f ') and line.endswith('.plugin') and line[2:] not in keep)]
    staged_lines.insert(1, 'd Platforms: ' + ' '.join(selected))
    entries = [line[2:] for line in staged_lines if line.startswith('f ')]
    for entry in entries + ['stata.toc']:
        if Path(entry).name != entry or not (source / entry).is_file():
            raise ValueError(f'missing or invalid package entry: {entry}')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.ctools-package-', dir=output.parent) as tmp:
        stage_dir = Path(tmp) / 'package'
        stage_dir.mkdir()
        for entry in entries + ['stata.toc']:
            shutil.copy2(source / entry, stage_dir / entry)
        (stage_dir / 'ctools.pkg').write_text('\n'.join(staged_lines) + '\n')
        errors = check_package(stage_dir)
        if errors: raise ValueError('\n'.join(errors))
        hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(stage_dir.iterdir())}
        (stage_dir / 'BUILD_INFO.json').write_text(json.dumps({
            'revision': revision, 'platforms': selected, 'sha256': hashes}, indent=2) + '\n')
        stage_dir.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT/'build')
    parser.add_argument('--output', type=Path, default=ROOT/'dist/ctools')
    parser.add_argument('--platform', choices=(*PLATFORMS, 'all'), default=native_platform())
    parser.add_argument('--revision', required=True, help='revision embedded when compiling the selected plugins')
    args = parser.parse_args()
    selected = list(PLATFORMS) if args.platform == 'all' else [args.platform]
    try: result = stage(args.source,args.output,selected,args.revision)
    except ValueError as error: raise SystemExit(str(error))
    print(f'Complete {", ".join(selected)} package: {result}')

if __name__ == '__main__': main()
