"""Check the Stata install manifest against a staged release directory."""
import argparse
from pathlib import Path


def check_package(directory, metadata_only=False):
    entries = [line.split(maxsplit=1)[1] for line in
               (directory / "ctools.pkg").read_text().splitlines()
               if line.startswith("f ")]
    errors = []
    platforms = {"mac_arm", "mac_x86", "windows", "linux"}
    scope = [line.removeprefix("d Platforms: ").split() for line in
             (directory / "ctools.pkg").read_text().splitlines() if line.startswith("d Platforms: ")]
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
        if name not in entries: errors.append(f"license file not in manifest: {name}")
    if len(entries) != len(set(entries)):
        errors.append("duplicate manifest entries")
    inventory = {p.name for p in directory.iterdir()
                 if p.suffix in {".ado", ".sthlp"}}
    for name in sorted(inventory - set(entries)):
        errors.append(f"not in manifest: {name}")
    for name in entries + ["stata.toc"]:
        if metadata_only and name.endswith(".plugin"):
            continue
        if not (directory / name).is_file():
            errors.append(f"missing package file: {name}")
    return errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--metadata-only", action="store_true",
                        help="check source metadata before platform builds")
    args = parser.parse_args()
    errors = check_package(args.directory, args.metadata_only)
    if errors:
        raise SystemExit("\n".join(errors))
    print("Package metadata checked" if args.metadata_only else "Package complete")
