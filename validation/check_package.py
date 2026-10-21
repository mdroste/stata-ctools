"""Check a Stata install manifest against a staged directory."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.package_contract import check_package


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--metadata-only", action="store_true",
                        help="check source metadata before platform builds")
    parser.add_argument("--integrity", action="store_true", help="require and verify the staged release checksum manifest")
    args = parser.parse_args()
    errors = check_package(args.directory, args.metadata_only, args.integrity)
    if errors:
        raise SystemExit("\n".join(errors))
    print("Package metadata checked" if args.metadata_only else "Package complete")
