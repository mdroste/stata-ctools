"""Run sequential B-C-C-B benchmarks through the project's Stata launcher.

Benchmark cencode against native encode. Use frozen ado/plugin directories. Every session checks all rows against the
reference command and baseline. Warmups, phases, and individual observations
are retained; only positive repetition numbers enter the timing summaries.
Generated datasets are removed on success or failure unless --keep-data is set.
"""

import argparse
import csv
import hashlib
import json
import shutil
import statistics
import subprocess
from pathlib import Path

from validation.support.stata_runner import benchmark_data, launch_driver

from validation import ROOT
LABELS = ("b1", "c1", "c2", "b2")


def summarize(out):
    records, phases, groups = [], [], {}
    for label in LABELS:
        text = (out / f"{label}.log").read_text()
        if text.splitlines().count("ENCODE_BENCHMARK_RC=0") != 1:
            raise ValueError("Incomplete run: " + label)
        for line in text.splitlines():
            if line.startswith("PERF,"):
                _, tag, case, op, rows, threads, rep, seconds = line.split(",")
                row = (tag, case, op, int(rows), int(threads), int(rep), float(seconds))
                records.append(row)
                if int(rep):
                    key = (case, op, int(rows), int(threads))
                    groups.setdefault(key, {}).setdefault(tag[0], []).append(float(seconds))
            elif line.startswith("PHASE,"):
                phases.append(line.split(",")[1:])
    for filename, header, values in (
        ("timings.csv", "run,case,operation,rows,threads,rep,seconds", records),
        ("phases.csv", "run,case,threads,rep,phase,seconds", phases),
    ):
        with (out / filename).open("w") as f:
            w = csv.writer(f)
            w.writerow(header.split(","))
            w.writerows(values)
    summary = []
    for (case, op, rows, threads), samples in groups.items():
        if op != "ctools":
            continue
        b, c = samples["b"], samples["c"]
        if len(b) != len(c):
            raise ValueError(f"Unbalanced baseline/candidate samples: {case}, {threads}")
        refs = groups[case, "reference", rows, 0]
        ref = refs["b"] + refs["c"]
        bm, cm, rm = map(statistics.median, (b, c, ref))
        summary.append(
            {
                "case": case,
                "rows": rows,
                "threads": threads,
                "baseline_seconds": bm,
                "candidate_seconds": cm,
                "reference_seconds": rm,
                "improvement": bm / cm,
                "versus_reference": rm / cm,
                "samples_per_build": len(b),
                "baseline_min": min(b),
                "baseline_max": max(b),
                "candidate_min": min(c),
                "candidate_max": max(c),
            }
        )
    with (out / "summary.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=summary[0])
        w.writeheader()
        w.writerows(summary)
    for r in summary:
        print(
            f"{r['case']:12} t={r['threads']} {r['baseline_seconds']:.6f} -> "
            f"{r['candidate_seconds']:.6f} ({r['improvement']:.2f}x); "
            f"reference {r['reference_seconds']:.6f} ({r['versus_reference']:.2f}x)"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rows", type=int, default=1000000)
    parser.add_argument("--reps", type=int, default=4)
    parser.add_argument("--threads", default="1 8")
    parser.add_argument("--cases", default="")
    parser.add_argument("--refreps", type=int)
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--keep-data", action="store_true",
                        help="retain generated comparison datasets for debugging")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.summarize_only:
        summarize(out)
        return
    if args.baseline is None or args.candidate is None:
        parser.error("provide --baseline and --candidate, or use --summarize-only")
    if args.rows < 1 or args.reps < 1 or (args.refreps is not None and args.refreps < 1):
        parser.error("rows, repetitions, and reference repetitions must be positive")
    builds = {"b": args.baseline.resolve(), "c": args.candidate.resolve()}
    hashes = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for build in builds.values()
        for p in [*build.glob("*.plugin"), *build.glob("*.ado")]
    }
    (out / "manifest.json").write_text(json.dumps(hashes, indent=2) + "\n")
    subprocess.run(
        [
            "clang",
            "-std=c11",
            "-O2",
            "-I" + str(ROOT / "src"),
            "-bundle",
            "-DSYSTEM=APPLEMAC",
            "-D_DARWIN_C_SOURCE",
            str(ROOT / "validation/performance/benchmark_clock.c"),
            str(ROOT / "src/stplugin.c"),
            "-o",
            str(out / "clock.plugin"),
        ],
        check=True,
    )
    shutil.copy2(ROOT / "validation/performance/benchmark_encode_performance.do", out / "benchmark.do")
    with benchmark_data(out, keep=args.keep_data) as data:
        for label in LABELS:
            if any(hashlib.sha256(Path(p).read_bytes()).hexdigest() != h for p, h in hashes.items()):
                raise ValueError("Frozen build changed")
            driver = out / f"{label}.do"
            driver.write_text(f'''clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{builds[label[0]]}" {label} "{out}" "{args.cases}" {args.reps} "{args.threads}" {args.rows} "{args.refreps if args.refreps is not None else args.reps}" "{data}"
local rc=_rc
di "ENCODE_BENCHMARK_RC=`rc'"
log close
exit, clear
''')
            print("Running " + label, flush=True)
            launch_driver(driver, out / f"{label}.log")
        summarize(out)


if __name__ == "__main__":
    main()
