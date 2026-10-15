"""Stage editable package sources into the published build/ install tree."""

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def generated():
    for folder in ("ado", "help", "metadata"):
        for source in sorted((ROOT / "package" / folder).iterdir()):
            if source.is_file():
                yield ROOT / "build" / source.name, source.read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    stale = []
    for output, content in generated():
        if not output.exists() or output.read_bytes() != content:
            if args.check:
                stale.append(output.name)
            else:
                output.parent.mkdir(exist_ok=True)
                temporary = output.with_suffix(output.suffix + ".tmp")
                temporary.write_bytes(content)
                temporary.replace(output)
    if stale:
        raise SystemExit("Stale package copies: " + ", ".join(stale))
    print("Package sources checked" if args.check else "Package sources staged")


if __name__ == "__main__":
    main()
