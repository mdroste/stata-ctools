"""Bench paths: shared execution with explicit workload-specific checks."""

import argparse
import csv
import hashlib
import json
import math
import platform
import shlex
import subprocess
import threading
import time
from collections import defaultdict
from pathlib import Path
from statistics import median

from validation.stata import ROOT, benchmark_data, launch_driver, stata_command, write_workload

sort_CASES = {
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
sort_DEFAULT_CASES = tuple(sort_CASES)
sort_CASES.update(
    {
        "nearly_sorted": 1_000_000,
        "nearly_sorted_middle": 1_000_000,
        "ties": 1_000_000,
        "constant_ties": 1_000_000,
        "long_prefix": 100_000,
    }
)
sort_STAGES = ("wall", "load", "sort", "permute", "store", "stream", "cleanup", "total")
sort_FIELDS = (
    "run",
    "case",
    "rows",
    "variables",
    "threads",
    "rep",
    *sort_STAGES,
    "stream_mode",
    "threads_max",
)


def sort_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sort_case_columns(case):
    return {"wide": 32, "verywide": 256, "ties": 5, "constant_ties": 5}.get(case, 4)


def sort_quote(value):
    value = str(value)
    if any(c in value for c in ('"', "`", "\n", "\r")):
        raise ValueError("unsupported character in Stata argument: " + value)
    return '"' + value + '"'


def sort_summarize(out, manifest):
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
            if len(values) != len(sort_FIELDS):
                raise ValueError("malformed sample: " + line)
            row = dict(zip(sort_FIELDS, values))
            for field in ("rows", "variables", "threads", "rep", "stream_mode", "threads_max"):
                row[field] = int(row[field])
            for stage in sort_STAGES:
                row[stage] = float(row[stage])
            case, team, rep = row["case"], row["threads"], row["rep"]
            key = (case, team, rep)
            columns = manifest.get("columns", {}).get(case, sort_case_columns(case))
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
                or any(not math.isfinite(row[s]) or row[s] < 0 for s in sort_STAGES)
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
                **{stage: median(r[stage] for r in samples) for stage in sort_STAGES},
            }
        )
    for name, fields, records in (
        ("timings.csv", sort_FIELDS, rows),
        ("summary.csv", ("variant", "case", "threads", "samples", *sort_STAGES), summaries),
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


def sort_main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cases", nargs="+", choices=sort_CASES, default=sort_DEFAULT_CASES)
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 4, 8])
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--rows", action="append", default=[], metavar="CASE=N")
    parser.add_argument(
        "--labels", nargs="+", choices=("b1", "c1", "c2", "b2"), default=["b1", "c1", "c2", "b2"]
    )
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--keep-data", action="store_true")
    args = parser.parse_args(argv)
    out = args.output.resolve()
    if args.summarize:
        sort_summarize(out, json.loads((out / "manifest.json").read_text()))
        return
    for values in (args.cases, args.threads, args.labels):
        if len(set(values)) != len(values):
            parser.error("duplicate cases, thread counts, or labels")
    if args.reps < 1 or min(args.threads) < 1:
        parser.error("repetitions and thread counts must be positive")
    cases = {case: sort_CASES[case] for case in args.cases}
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
        "columns": {case: sort_case_columns(case) for case in cases},
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
            "files": {str(f): sort_digest(f) for f in files},
        }
    source = write_workload("csort_performance", out / "benchmark.do")
    flags = ["-std=c11", "-O2", "-I", str(ROOT / "src")]
    flags += (
        ["-bundle", "-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
        if mac
        else ["-shared", "-fPIC", "-DSYSTEM=STUNIX", "-D_POSIX_C_SOURCE=200809L"]
    )
    command = [
        "cc",
        *flags,
        str(ROOT / "validation/performance/benchmark_clock.c"),
        str(ROOT / "src/stplugin.c"),
        "-o",
        str(out / "clock.plugin"),
    ]
    subprocess.run(command, check=True)
    manifest["clock_compile_command"] = command
    manifest["harness_files"] = {
        str(p): sort_digest(p) for p in (source, Path(__file__).resolve(), out / "clock.plugin")
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
            driver.write_text(
                f"""clear all
set more off
log using {sort_quote(out / (label + ".log"))}, text replace
capture noisily do {sort_quote(out / "benchmark.do")} {" ".join(map(sort_quote, call_args))}
local rc = _rc
di "CSORT_DRIVER_RC=`rc'"
log close
exit, clear
"""
            )
            if args.run:
                for path, expected in build["files"].items():
                    if sort_digest(Path(path)) != expected:
                        raise ValueError("frozen build changed: " + path)
                print("Running " + label, flush=True)
                launch_driver(driver, out / f"{label}.log")
            else:
                print(shlex.join(stata_command(driver)))
        if args.run:
            sort_summarize(out, manifest)


matching_CASES = {
    "nr": 160000,
    "nr_equal": 160000,
    "nr_desc": 160000,
    "nr_common": 160000,
    "radius": 160000,
    "radius_narrow": 160000,
    "radius_default": 80000,
    "nearest": 1000000,
    "kernel": 20000,
}
matching_STAGES = ["wall", "load", "setup", "sort", "match", "stats", "store", "total"]


def matching_digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def matching_quote(value):
    value = str(value)
    if any(c in value for c in ['"', "`", "\n"]):
        raise ValueError("unsupported quote")
    return '"' + value + '"'


def matching_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--variant",
        action="append",
        required=True,
        help="label=build directory; use b1 first with --parity",
    )
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--cases", nargs="+")
    p.add_argument("--teams", nargs="+", type=int, default=[1, 12])
    p.add_argument("--reps", type=int, default=5)
    p.add_argument("--parity", action="store_true")
    p.add_argument(
        "--allow-busy",
        action="store_true",
        help="Record competing Stata activity instead of rejecting timing runs",
    )
    a = p.parse_args(argv)
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    variants = {k: Path(v).resolve() for k, v in (x.split("=", 1) for x in a.variant)}
    cases = {k: int(v) for k, v in (x.split(":") for x in a.cases)} if a.cases else matching_CASES
    files = {
        k: {p.name: matching_digest(p) for p in [*v.glob("*.ado"), v / "ctools_mac_arm.plugin"]}
        for k, v in variants.items()
    }
    manifest = {
        "variants": {k: str(v) for k, v in variants.items()},
        "files": files,
        "cases": cases,
        "teams": a.teams,
        "reps": a.reps,
        "parity": a.parity,
        "allow_busy": a.allow_busy,
        "seed": 95183,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    frozen = write_workload("cpsmatch_paths", out / "benchmark.do")
    clock = out / "clock.plugin"
    subprocess.run(
        [
            "clang",
            "-O2",
            "-std=c11",
            "-bundle",
            "-DSYSTEM=APPLEMAC",
            "-D_DARWIN_C_SOURCE",
            "-I",
            str(ROOT / "src"),
            str(ROOT / "validation/performance/benchmark_clock.c"),
            str(ROOT / "src/stplugin.c"),
            "-o",
            str(clock),
        ],
        check=True,
    )
    reference = out / "reference"
    reference.mkdir(exist_ok=True)
    samples = []
    results = []
    for case, n in cases.items():
        for label, build in variants.items():
            ps = subprocess.check_output(["ps", "-axo", "args="], text=True)
            if not (a.parity or a.allow_busy) and any(
                "/stata-mp " in line for line in ps.splitlines()
            ):
                raise RuntimeError("competing Stata")
            assert all(matching_digest(build / name) == h for name, h in files[label].items())
            log = out / f"{case}_{label}.txt"
            driver = out / f"{case}_{label}.do"
            args = [
                build,
                clock,
                case,
                n,
                " ".join(map(str, a.teams)),
                a.reps,
                label,
                95183,
                reference if a.parity else "",
            ]
            driver.write_text(
                "clear all\nset more off\nset linesize 255\nlog using "
                + matching_quote(log)
                + ", text replace\ncapture noisily do "
                + matching_quote(frozen)
                + " "
                + " ".join(map(matching_quote, args))
                + '\nlocal rc=_rc\ndi "DRIVER_RC=`rc\'"\nlog close\nexit, clear\n'
            )
            finished = threading.Event()
            activity = []

            def monitor(finished=finished, activity=activity):
                while not finished.is_set():
                    ps = subprocess.check_output(["ps", "-axo", "pid=,pcpu=,args="], text=True)
                    activity.append(
                        {
                            "time": time.time(),
                            "jobs": [s.strip() for s in ps.splitlines() if "/stata-mp " in s],
                        }
                    )
                    finished.wait(2)

            worker = threading.Thread(target=monitor)
            worker.start()
            print("RUN", case, label, flush=True)
            try:
                text = launch_driver(driver, log, cwd=out)
            finally:
                finished.set()
                worker.join()
                (out / f"{case}_{label}_activity.json").write_text(
                    json.dumps(activity, indent=2) + "\n"
                )
            assert (
                a.parity
                or a.allow_busy
                or not any("/ctools-stata-" not in job for r in activity for job in r["jobs"])
            )
            assert f"\nPATH_COMPLETE,{label},{case}\n" in text
            seen = set()
            for line in text.splitlines():
                if line.startswith("RESULT,"):
                    results.append(line)
                if not line.startswith("PATH,"):
                    continue
                v = next(csv.reader([line]))[1:]
                assert len(v) == 13
                tag, c, rows, team, rep = v[:5]
                team, rep = int(team), int(rep)
                values = list(map(float, v[5:]))
                assert all(math.isfinite(x) and x >= 0 for x in values)
                assert (team, rep) not in seen and tag == label and c == case and int(rows) == n
                seen.add((team, rep))
                samples.append(
                    dict(
                        label=label,
                        case=case,
                        rows=n,
                        threads=team,
                        rep=rep,
                        **dict(zip(matching_STAGES, values)),
                    )
                )
            assert seen == {(t, r) for t in a.teams for r in range(1, a.reps + 1)}
            if a.parity and label != "b1":
                assert sum(s.startswith("PARITY,") for s in text.splitlines()) == len(a.teams)
    with (out / "timings.tsv").open("w") as f:
        w = csv.DictWriter(f, samples[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(samples)
    (out / "results.txt").write_text("\n".join(results) + "\n")
    groups = defaultdict(list)
    for row in samples:
        groups[row["label"], row["case"], row["threads"]].append(row)
    summary = [
        dict(
            label=k[0],
            case=k[1],
            threads=k[2],
            samples=len(v),
            **{s: median(r[s] for r in v) for s in matching_STAGES},
        )
        for k, v in groups.items()
    ]
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("CPSMATCH_PATH_CAMPAIGN_COMPLETE samples=" + str(len(samples)), flush=True)
