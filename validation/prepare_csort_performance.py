"""Prepare or run serial ABBA csort measurements against two frozen builds.

Each sample retains monotonic wall time and C stage timings; exact comparisons
against native stable sorting run outside timing. No observations are trimmed.
Use --labels b1 for baseline screening, --rows CASE=N for threshold probes, or
--summarize to re-read an existing campaign without launching Stata.
Optional adversarial cases: nearly_sorted, nearly_sorted_middle, ties,
constant_ties, long_prefix.
"""

import argparse
import csv
import hashlib
import json
import math
import platform
import shlex
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import median

from stata_runner import ROOT, benchmark_data, launch_driver, stata_command

CASES = {
    "continuous": 1_000_000,
    "integer": 1_000_000,
    "panel": 1_000_000,
    "mixed": 500_000,
    "long_string": 100_000,
    "long_string_large": 1_000_000,
    "sorted": 1_000_000,
    "reverse": 1_000_000,
    "wide": 200_000,
    "verywide": 10_000,
    "small": 1_000,
}
DEFAULT_CASES = tuple(CASES)
CASES.update(
    {
        "nearly_sorted": 1_000_000,
        "nearly_sorted_middle": 1_000_000,
        "ties": 1_000_000,
        "constant_ties": 1_000_000,
        "long_prefix": 100_000,
    }
)
STAGES = ("wall", "load", "sort", "permute", "store", "stream", "cleanup", "total")
FIELDS = (
    "run",
    "case",
    "rows",
    "variables",
    "threads",
    "rep",
    *STAGES,
    "stream_mode",
    "threads_max",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def case_columns(case):
    return {"wide": 32, "verywide": 256, "ties": 5, "constant_ties": 5}.get(case, 4)


def quote(value):
    value = str(value)
    if any(c in value for c in ('"', "`", "\n", "\r")):
        raise ValueError("unsupported character in Stata argument: " + value)
    return '"' + value + '"'


def summarize(out, manifest):
    rows = []
    groups = defaultdict(list)
    cases = manifest["cases"]
    teams = manifest["threads"]
    reps = manifest["repetitions"]
    for label in manifest["labels"]:
        lines = (out / f"{label}.log").read_text().splitlines()
        if (
            f"CSORT_PERFORMANCE_COMPLETE label={label}" not in lines
            or "CSORT_DRIVER_RC=0" not in lines
        ):
            raise ValueError("incomplete or failed session: " + label)
        checks = {line for line in lines if line.startswith("CSORT_CHECK,")}
        expected_checks = {
            f"CSORT_CHECK,{label},{case},{team},nosortedby,PASS" for case in cases for team in teams
        }
        if checks != expected_checks:
            raise ValueError("missing or unexpected plugin-only checks: " + label)
        seen = set()
        for line in lines:
            if not line.startswith("CSORT_PERF,"):
                continue
            values = next(csv.reader([line]))[1:]
            if len(values) != len(FIELDS):
                raise ValueError("malformed sample: " + line)
            row = dict(zip(FIELDS, values))
            for field in ("rows", "variables", "threads", "rep", "stream_mode", "threads_max"):
                row[field] = int(row[field])
            for stage in STAGES:
                row[stage] = float(row[stage])
            case, team, rep = row["case"], row["threads"], row["rep"]
            key = (case, team, rep)
            columns = manifest.get("columns", {}).get(case, case_columns(case))
            if (
                row["run"] != label
                or case not in cases
                or team not in teams
                or rep not in range(1, reps + 1)
                or key in seen
                or row["rows"] != cases[case]
                or row["variables"] != columns
                or row["stream_mode"] not in (0, 1)
                or row["threads_max"] < 1
                or row["wall"] <= 0
                or any(not math.isfinite(row[s]) or row[s] < 0 for s in STAGES)
            ):
                raise ValueError("invalid or duplicate sample: " + line)
            seen.add(key)
            rows.append(row)
            groups[label[0], case, team].append(row)
        expected = {
            (case, team, rep) for case in cases for team in teams for rep in range(1, reps + 1)
        }
        if seen != expected:
            raise ValueError(f"missing samples in {label}: {expected - seen}")
    summaries = []
    for (variant, case, team), samples in groups.items():
        summaries.append(
            {
                "variant": variant,
                "case": case,
                "threads": team,
                "samples": len(samples),
                **{stage: median(r[stage] for r in samples) for stage in STAGES},
            }
        )
    for name, fields, records in (
        ("timings.csv", FIELDS, rows),
        ("summary.csv", ("variant", "case", "threads", "samples", *STAGES), summaries),
    ):
        with (out / name).open("w", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
    lookup = {(r["variant"], r["case"], r["threads"]): r for r in summaries}
    print("| Case | Threads | Baseline (s) | Candidate (s) | Speedup |")
    print("|---|---:|---:|---:|---:|")
    for case in cases:
        for team in teams:
            before = lookup.get(("b", case, team))
            after = lookup.get(("c", case, team))
            b = f"{before['wall']:.6f}" if before else "—"
            c = f"{after['wall']:.6f}" if after else "—"
            ratio = f"{before['wall'] / after['wall']:.3f}×" if before and after else "—"
            print(f"| {case} | {team} | {b} | {c} | {ratio} |")
    print(f"Retained {len(rows)} measured samples; warmups excluded, no outliers removed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=DEFAULT_CASES)
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 4, 8])
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--rows", action="append", default=[], metavar="CASE=N")
    parser.add_argument(
        "--labels", nargs="+", choices=("b1", "c1", "c2", "b2"), default=["b1", "c1", "c2", "b2"]
    )
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--keep-data", action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    if args.summarize:
        summarize(out, json.loads((out / "manifest.json").read_text()))
        return
    for values in (args.cases, args.threads, args.labels):
        if len(set(values)) != len(values):
            parser.error("duplicate cases, thread counts, or labels")
    if args.reps < 1 or min(args.threads) < 1:
        parser.error("repetitions and thread counts must be positive")
    cases = {case: CASES[case] for case in args.cases}
    for override in args.rows:
        try:
            case, raw = override.split("=")
            n = int(raw)
            if case not in cases or not 1 <= n <= 2_147_483_647:
                raise ValueError
            cases[case] = n
        except ValueError:
            parser.error("--rows requires a selected CASE and N in 1..2147483647")
    builds = {"b": args.baseline, "c": args.candidate}
    for label in args.labels:
        if builds[label[0]] is None:
            parser.error("selected labels require --baseline and/or --candidate")
    out.mkdir(parents=True, exist_ok=True)
    mac = platform.system() == "Darwin"
    if platform.system() == "Windows":
        parser.error("benchmark clock compilation currently supports macOS and Linux")
    plugin = (
        ("ctools_mac_arm.plugin" if platform.machine() == "arm64" else "ctools_mac_x86.plugin")
        if mac
        else "ctools_linux.plugin"
    )
    manifest = {
        "cases": cases,
        "columns": {case: case_columns(case) for case in cases},
        "threads": args.threads,
        "repetitions": args.reps,
        "labels": args.labels,
        "seed": 95183,
        "machine": platform.platform(),
        "warmups_per_case_thread": 1,
        "plugin_only_checks_per_case_thread": 1,
        "builds": {},
    }
    for tag in dict.fromkeys(label[0] for label in args.labels):
        build = builds[tag].resolve()
        files = [build / plugin, *sorted(build.glob("*.ado"))]
        if not (build / "csort.ado").is_file():
            parser.error("frozen build lacks csort.ado: " + str(build))
        manifest["builds"][tag] = {
            "directory": str(build),
            "files": {str(f): digest(f) for f in files},
        }
    source = ROOT / "validation/benchmark_csort_performance.do"
    shutil.copy2(source, out / "benchmark.do")
    flags = ["-std=c11", "-O2", "-I", str(ROOT / "src")]
    flags += (
        ["-bundle", "-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
        if mac
        else ["-shared", "-fPIC", "-DSYSTEM=STUNIX", "-D_POSIX_C_SOURCE=200809L"]
    )
    command = [
        "cc",
        *flags,
        str(ROOT / "validation/benchmark_clock.c"),
        str(ROOT / "src/stplugin.c"),
        "-o",
        str(out / "clock.plugin"),
    ]
    subprocess.run(command, check=True)
    manifest["clock_compile_command"] = command
    manifest["harness_files"] = {
        str(p): digest(p) for p in (source, Path(__file__).resolve(), out / "clock.plugin")
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Separate processes avoid loading two static OpenMP runtimes into Stata.
    with benchmark_data(out, keep=not args.run or args.keep_data) as data:
        for label in args.labels:
            build = manifest["builds"][label[0]]
            call_args = [
                build["directory"],
                label,
                out,
                " ".join(f"{c}:{n}" for c, n in cases.items()),
                " ".join(map(str, args.threads)),
                args.reps,
                data,
                int(args.keep_data),
            ]
            driver = out / f"{label}.do"
            driver.write_text(f"""clear all
set more off
log using {quote(out / (label + ".log"))}, text replace
capture noisily do {quote(out / "benchmark.do")} {" ".join(map(quote, call_args))}
local rc = _rc
di "CSORT_DRIVER_RC=`rc'"
log close
exit, clear
""")
            if args.run:
                for path, expected in build["files"].items():
                    if digest(Path(path)) != expected:
                        raise ValueError("frozen build changed: " + path)
                print("Running " + label, flush=True)
                launch_driver(driver, out / f"{label}.log")
            else:
                print(shlex.join(stata_command(driver)))
        if args.run:
            summarize(out, manifest)


if __name__ == "__main__":
    main()
