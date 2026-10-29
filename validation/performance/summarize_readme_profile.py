"""Print one-shot speedup ranges; optionally fill the main README table."""

import argparse
import csv
import math
import re
from pathlib import Path

COMMANDS = {
    "import": "import",
    "export": "export",
    "sort": "sort",
    "merge": "merge",
    "sample": "sample",
    "bsample": "bsample",
    "encode": "encode",
    "decode": "decode",
    "destring": "destring",
    "ipolate": "ipolate",
    "split": "split",
    "rangejoin": "rangejoin",
    "qreg": "qreg",
    "winsor": "gstats winsor",
    "rangestat": "rangestat",
    "psmatch2": "psmatch2",
    "binscatter": "binscatter",
    "reghdfe": "reghdfe",
    "ivreghdfe": "ivreghdfe",
    "ppmlhdfe": "ppmlhdfe",
}
SIZES = (1_000_000, 10_000_000, 25_000_000)


def rounded(ratio, *, upper):
    """Round range endpoints outward: tenths below 2x, units below 10x, then tens."""
    rounding = math.ceil if upper else math.floor
    if ratio < 2:
        return f"{rounding(ratio * 10) / 10:.1f}"
    step = 1 if ratio < 10 else 10
    return str(rounding(ratio / step) * step)


def summarize(path, sizes):
    rows = {}
    with path.open(newline="") as source:
        for row in csv.DictReader(source):
            key = (row["command"], int(row["observations"]), row["method"])
            if key in rows:
                raise ValueError(f"duplicate timing: {key}")
            rows[key] = row
    expected = {(cmd, n, method) for cmd in COMMANDS for n in sizes for method in ("ref", "c")}
    if rows.keys() != expected:
        raise ValueError(f"unexpected/missing timings: {rows.keys() ^ expected}")
    output = [
        "| Command | Observations | Reference (s) | ctools (s) | Speedup |",
        "|---|---:|---:|---:|---:|",
    ]
    ranges = {}
    for cmd, label in COMMANDS.items():
        ratios = []
        for n in sizes:
            pair = [rows[cmd, n, method] for method in ("ref", "c")]
            times = [float(row["seconds"]) for row in pair]
            valid = all(
                row["rc"] == "0" and row["check_rc"] == "0" and math.isfinite(t) and t > 0
                for row, t in zip(pair, times)
            )
            if valid:
                ratio = times[0] / times[1]
                ratios.append(ratio)
                result = f"{ratio:.2f}×"
            else:
                result = "— (failed check/command or timer resolution)"
            output.append(f"| `{label}` | {n:,} | {times[0]:.3f} | {times[1]:.3f} | {result} |")
        if len(ratios) == len(sizes):
            low = rounded(min(ratios), upper=False)
            high = rounded(max(ratios), upper=True)
            ranges[label] = f"{low}×" if low == high else f"{low}–{high}×"
        else:
            ranges[label] = "—"
    return "\n".join(output) + "\n", ranges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, help="directory containing profile.log and timings.csv"
    )
    parser.add_argument(
        "--sizes", nargs="+", type=int, default=SIZES, help="override only for smoke tests"
    )
    parser.add_argument("--update-readme", type=Path)
    args = parser.parse_args()
    log = (args.directory / "profile.log").read_text()
    if not re.search(r"^PROFILE_COMPLETE$", log, re.MULTILINE) or not re.search(
        r"^PROFILE_EXIT rc=0$", log, re.MULTILINE
    ):
        parser.error("profile did not finish successfully; see profile.log")
    if args.update_readme and tuple(args.sizes) != SIZES:
        parser.error("README publication requires the full 1M/10M/25M run")
    try:
        details, ranges = summarize(args.directory / "timings.csv", args.sizes)
    except (ValueError, KeyError) as error:
        parser.error(str(error))
    (args.directory / "results.md").write_text(details)
    for command, speedup in ranges.items():
        print(f"{command:14} {speedup}")
    if args.update_readme:
        readme = args.update_readme.read_text()
        for command, speedup in ranges.items():
            pattern = rf"^(\| `{re.escape(command)}` \|.*\|) [^|]+ \|$"
            readme, count = re.subn(
                pattern, lambda m, value=speedup: f"{m[1]} {value} |", readme, flags=re.MULTILINE
            )
            if count != 1:
                parser.error(f"expected exactly one README row for {command}; no changes written")
        args.update_readme.write_text(readme)


if __name__ == "__main__":
    main()
