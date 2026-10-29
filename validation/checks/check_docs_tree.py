"""Reject dated snapshots and run artifacts under docs/.

docs/ holds current documentation. docs/benchmarks/ holds only the evidence that
scripts regenerate in place; dated experiments belong in the gitignored reports/.
"""
import re
import subprocess
import sys
from pathlib import Path

from validation import ROOT
BENCHMARKS = ("docs/benchmarks/hdfe_datasets_", "docs/benchmarks/readme_profile/")
DATED = re.compile(r"(?<!\d)(19|20)\d{2}-?[01]\d-?[0-3]\d(?!\d)")
ARTIFACTS = {".do", ".py", ".plugin", ".log", ".patch", ".dta", ".o"}


def main():
    files = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "docs"],
        cwd=ROOT, check=True, capture_output=True, text=True).stdout.split()
    files = [name for name in files if (ROOT / name).exists()]
    errors = []
    for name in files:
        if DATED.search(name):
            errors.append(f"{name}: dated path; move snapshots to reports/")
        elif Path(name).suffix in ARTIFACTS:
            errors.append(f"{name}: run artifact; keep drivers in validation/ and output in temp/")
        elif name.startswith("docs/benchmarks/") and not name.startswith(BENCHMARKS):
            errors.append(f"{name}: not regenerated evidence; see docs/development/README.md")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"docs tree checked ({len(files)} files)")


if __name__ == "__main__":
    sys.exit(main())
