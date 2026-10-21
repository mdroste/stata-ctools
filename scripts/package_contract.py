"""Shared package invariants for staging and validation."""

from pathlib import Path
import hashlib
import json
import re


def check_package(directory, metadata_only=False, integrity=False):
    entries = [
        line.split(maxsplit=1)[1]
        for line in (directory / "ctools.pkg").read_text().splitlines()
        if line.startswith("f ")
    ]
    errors = []
    platforms = {"mac_arm", "mac_x86", "windows", "linux"}
    scope = [
        line.removeprefix("d Platforms: ").split()
        for line in (directory / "ctools.pkg").read_text().splitlines()
        if line.startswith("d Platforms: ")
    ]
    if scope:
        if len(scope) != 1 or not scope[0] or not set(scope[0]) <= platforms:
            errors.append("invalid explicit platform scope")
        else:
            platforms = set(scope[0])
    required_plugins = {f"ctools_{name}.plugin" for name in platforms}
    listed_plugins = {name for name in entries if name.endswith(".plugin")}
    for name in sorted(listed_plugins - required_plugins):
        errors.append(f"binary outside declared platform scope: {name}")
    for name in sorted(required_plugins - set(entries)):
        errors.append(f"platform binary not in manifest: {name}")
    for name in ("LICENSE", "THIRD_PARTY_NOTICES"):
        if name not in entries:
            errors.append(f"license file not in manifest: {name}")
    if len(entries) != len(set(entries)):
        errors.append("duplicate manifest entries")
    inventory = {p.name for p in directory.iterdir() if p.suffix in {".ado", ".sthlp"}}
    for name in sorted(inventory - set(entries)):
        errors.append(f"not in manifest: {name}")
    for name in entries + ["stata.toc"]:
        if metadata_only and name.endswith(".plugin"):
            continue
        if not (directory / name).is_file():
            errors.append(f"missing package file: {name}")
    if integrity:
        if metadata_only:
            errors.append('integrity validation cannot skip binary payloads')
        try:
            info = json.loads((directory / 'BUILD_INFO.json').read_text())
            hashes = info['sha256']
            expected = set(entries) | {'ctools.pkg', 'stata.toc'}
            if not isinstance(hashes, dict) or set(hashes) != expected:
                errors.append('checksum inventory does not match package manifest')
            if set(info['platforms']) != platforms or len(info['platforms']) != len(platforms):
                errors.append('checksum platform scope does not match package manifest')
            if not isinstance(info['revision'], str) or not info['revision']:
                errors.append('missing build revision')
            if isinstance(hashes, dict):
                for name, digest in hashes.items():
                    if Path(name).name != name or not isinstance(digest, str) or not re.fullmatch(r'[0-9a-f]{64}', digest):
                        errors.append(f'invalid checksum entry: {name}')
                    elif not (directory / name).is_file():
                        errors.append(f'missing checksummed file: {name}')
                    else:
                        with (directory / name).open('rb') as payload:
                            hasher = hashlib.sha256()
                            for chunk in iter(lambda: payload.read(1024 * 1024), b''):
                                hasher.update(chunk)
                            actual = hasher.hexdigest()
                        if actual != digest:
                            errors.append(f'checksum mismatch: {name}')
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(f'invalid or missing BUILD_INFO.json: {error}')
    return errors
