"""Bench: shared execution with explicit workload-specific checks."""

if __package__ in (None, ""):
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import argparse
import csv
import hashlib
import json
import math
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from statistics import median

from validation import ROOT
from validation.stata import benchmark_data, launch_driver, run_do, stata_command, write_workload

command_CASES = (
    "ip_sorted",
    "ip_shuffled",
    "ip_single",
    "ip_skew",
    "ip_string",
    "sp_short",
    "sp_long",
    "sp_multi",
    "rj_sorted",
    "rj_string",
    "rj_dense",
    "rj_wide",
)


def command_summarize(directory, output_csv, cases=command_CASES):
    rows = []
    timings = defaultdict(list)
    for label in ("b1", "c1", "c2", "b2"):
        text = (directory / (label + ".log")).read_text()
        if (
            f"\nPERFORMANCE_COMPLETE label={label}\n" not in text
            or "\nPERFORMANCE_DRIVER_RC=0\n" not in text
        ):
            raise SystemExit(f"Incomplete or failed run: {label}")
        seen = set()
        for line in text.splitlines():
            if not line.startswith("PERF,"):
                continue
            _, run, case, threads, rep, seconds, n = line.split(",")
            threads, rep, n = int(threads), int(rep), int(n)
            seconds = float(seconds)
            if (
                run != label
                or case not in cases
                or threads not in (1, 4, 12)
                or rep not in range(1, 7)
                or seconds <= 0
            ):
                raise SystemExit(f"Invalid benchmark record: {line}")
            key = (case, threads, rep)
            if key in seen:
                raise SystemExit(f"Duplicate benchmark record: {line}")
            seen.add(key)
            rows.append((run, case, threads, rep, seconds, n))
            timings[label[0], case, threads].append(seconds)
        if len(seen) != len(cases) * 3 * 6:
            raise SystemExit(f"Missing records: {label}: {len(seen)}")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("run", "case", "threads", "rep", "seconds", "output_rows"))
        writer.writerows(rows)
    print(
        "| Case | Before, 12 threads | After, 12 threads | Speedup | After, 1 thread | After, 4 threads |"
    )
    print("|---|---:|---:|---:|---:|---:|")
    for case in cases:
        b = median(timings["b", case, 12])
        c = median(timings["c", case, 12])
        one = median(timings["c", case, 1])
        four = median(timings["c", case, 4])
        print(f"| {case} | {b:.4f} | {c:.4f} | {b/c:.2f}× | {one:.4f} | {four:.4f} |")
    print(f"\nRetained all {len(rows)} observations; no timing outliers removed.")


def command_main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", default=command_CASES, choices=command_CASES)
    args = parser.parse_args(argv)
    command_summarize(args.directory, args.csv, args.cases)


def commands_main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--baseline", required=True, type=Path)
    p.add_argument("--candidate", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--cases", nargs="+", default=[])
    p.add_argument("--run", action="store_true")
    p.add_argument(
        "--keep-data", action="store_true", help="retain generated comparison datasets after --run"
    )
    a = p.parse_args(argv)
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    mac = platform.system() == "Darwin"
    plugin = (
        "ctools_mac_arm.plugin"
        if mac and platform.machine() == "arm64"
        else "ctools_mac_x86.plugin" if mac else "ctools_linux.plugin"
    )
    manifest = {}
    for tag, build in [("b", a.baseline.resolve()), ("c", a.candidate.resolve())]:
        files = [build / plugin] + [
            build / (command + ".ado") for command in ("cipolate", "csplit", "crangejoin")
        ]
        manifest[tag] = {str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    flags = ["-std=c11", "-O2", "-I", str(ROOT / "src")]
    flags += (
        ["-bundle", "-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
        if mac
        else ["-shared", "-fPIC", "-DSYSTEM=STUNIX", "-D_POSIX_C_SOURCE=200809L"]
    )
    subprocess.run(
        [
            "cc",
            *flags,
            str(ROOT / "validation/performance/benchmark_clock.c"),
            str(ROOT / "src/stplugin.c"),
            "-o",
            str(out / "clock.plugin"),
        ],
        check=True,
    )
    write_workload("command_performance", out / "benchmark.do")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Use separate processes: loading two plugins containing static OpenMP
    # runtimes into one Stata session is not a supported comparison method.
    with benchmark_data(out, keep=not a.run or a.keep_data) as data:
        for label in ("b1", "c1", "c2", "b2"):
            build = a.baseline.resolve() if label.startswith("b") else a.candidate.resolve()
            driver = out / (label + ".do")
            driver.write_text(
                f"""clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{build}" {label} "{out}" "{' '.join(a.cases)}" "{data}"
local rc=_rc
di "PERFORMANCE_DRIVER_RC=`rc'"
log close
exit, clear
"""
            )
            if a.run:
                for path, digest in manifest[label[0]].items():
                    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                        raise ValueError("Frozen build changed: " + path)
                print("Running " + label, flush=True)
                launch_driver(driver, out / (label + ".log"))
            else:
                print(shlex.join(stata_command(driver)))
        if a.run:
            command_summarize(out, out / "timings.csv", a.cases or command_CASES)


io_LABELS = ("b1", "c1", "c2", "b2")


def io_summarize(out, data=None):
    data = out if data is None else data
    (out / "comparisons.txt").unlink(missing_ok=True)
    measurements = {}
    records = []
    for label in io_LABELS:
        text = (out / (label + ".log")).read_text()
        if text.splitlines().count("IO_BENCHMARK_RC=0") != 1:
            raise ValueError("incomplete benchmark: " + label)
        for line in text.splitlines():
            if not line.startswith("PERF,"):
                continue
            _, tag, case, op, threads, rep, seconds = line.split(",")
            if tag != label:
                raise ValueError("unexpected benchmark label")
            record = (label, case, op, int(threads), int(rep), float(seconds))
            records.append(record)
            if int(rep):
                key = (case, op, int(threads))
                measurements.setdefault(key, {}).setdefault(label[0], []).append(float(seconds))
    with (out / "timings.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(("run", "case", "operation", "threads", "rep", "seconds"))
        writer.writerows(records)
    summary = []
    for (case, op, threads), samples in measurements.items():
        if not samples.get("b") or len(samples["b"]) != len(samples.get("c", [])):
            raise ValueError("unbalanced samples")
        b, c = (statistics.median(samples[tag]) for tag in ("b", "c"))
        summary.append(
            dict(
                case=case,
                operation=op,
                threads=threads,
                baseline_seconds=b,
                candidate_seconds=c,
                speedup=b / c,
                baseline_min=min(samples["b"]),
                baseline_max=max(samples["b"]),
                candidate_min=min(samples["c"]),
                candidate_max=max(samples["c"]),
                samples_per_build=len(samples["b"]),
            )
        )
    if not summary:
        raise ValueError("no measurements")
    with (out / "summary.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    checked = 0
    for case, op, threads in measurements:
        if not op.startswith("export_"):
            continue
        ext = "csv" if op == "export_csv" else "xlsx"
        expected = data / f"b1_{case}_{threads}.{ext}"
        for label in ("c1", "c2", "b2"):
            actual = data / f"{label}_{case}_{threads}.{ext}"
            if ext == "csv":
                if (
                    hashlib.sha256(expected.read_bytes()).digest()
                    != hashlib.sha256(actual.read_bytes()).digest()
                ):
                    raise ValueError("CSV bytes differ: " + str(actual))
            else:
                with zipfile.ZipFile(expected) as b, zipfile.ZipFile(actual) as c:
                    if b.namelist() != c.namelist() or c.testzip() is not None:
                        raise ValueError("invalid workbook: " + str(actual))
                    for name in b.namelist():
                        if b.read(name) != c.read(name):
                            raise ValueError(f"XLSX part differs: {actual}: {name}")
            checked += 1
    result = f"PASS {checked} exact export comparisons; imported values checked by Stata cf"
    (out / "comparisons.txt").write_text(result + "\n")
    print(result)
    for row in summary:
        print(
            f"{row['case']:8} {row['operation']:12} t={row['threads']:2} "
            f"{row['baseline_seconds']:.6f} -> {row['candidate_seconds']:.6f} s "
            f"({row['speedup']:.2f}x)"
        )
    return summary


def io_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", type=Path)
    p.add_argument("--candidate", type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument(
        "--cases",
        nargs="+",
        choices=["numeric", "mixed", "wide", "strings", "small", "large"],
        default=["numeric", "mixed", "wide", "strings"],
    )
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--run", action="store_true")
    p.add_argument("--summarize-only", action="store_true")
    p.add_argument(
        "--keep-data",
        action="store_true",
        help="retain generated data after --run for debugging or --summarize-only",
    )
    a = p.parse_args(argv)
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if a.summarize_only:
        io_summarize(out)
        return
    if a.baseline is None or a.candidate is None or a.reps < 1:
        p.error("provide --baseline, --candidate, and positive --reps")
    (out / "comparisons.txt").unlink(missing_ok=True)
    mac = platform.system() == "Darwin"
    plugin = (
        ("ctools_mac_arm.plugin" if platform.machine() == "arm64" else "ctools_mac_x86.plugin")
        if mac
        else "ctools_linux.plugin"
    )
    builds = {"b": a.baseline.resolve(), "c": a.candidate.resolve()}
    manifest = {}
    for tag, build in builds.items():
        manifest[tag] = {
            str(f): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in [build / plugin, *sorted(build.glob("*.ado"))]
        }
    flags = ["-std=c11", "-O2", "-I", str(ROOT / "src")]
    flags += (
        ["-bundle", "-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
        if mac
        else ["-shared", "-fPIC", "-DSYSTEM=STUNIX", "-D_POSIX_C_SOURCE=200809L"]
    )
    subprocess.run(
        [
            "cc",
            *flags,
            str(ROOT / "validation/performance/benchmark_clock.c"),
            str(ROOT / "src/stplugin.c"),
            "-o",
            str(out / "clock.plugin"),
        ],
        check=True,
    )
    write_workload("io_performance", out / "benchmark.do")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Prepared-only drivers keep stable paths; executed runs own a scratch directory.
    with benchmark_data(out, keep=not a.run or a.keep_data) as data:
        for label in io_LABELS:
            driver = out / (label + ".do")
            driver.write_text(
                f"""clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{builds[label[0]]}" {label} "{out}" "{' '.join(a.cases)}" {a.reps} "{data}"
local rc=_rc
di "IO_BENCHMARK_RC=`rc'"
log close
exit, clear
"""
            )
            if a.run:
                for path, digest in manifest[label[0]].items():
                    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                        raise ValueError("frozen build changed: " + path)
                print("Running " + label, flush=True)
                launch_driver(driver, out / (label + ".log"))
            else:
                print(shlex.join(stata_command(driver)))
        if a.run:
            io_summarize(out, data)


abba_LABELS = ("b1", "c1", "c2", "b2")


def abba_summarize(out, kind="encode"):
    records, phases, groups = [], [], {}
    for label in abba_LABELS:
        text = (out / f"{label}.log").read_text()
        if text.splitlines().count(f"{kind.upper()}_BENCHMARK_RC=0") != 1:
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
        (
            "phases.csv",
            (
                "run,case,threads,rep,phase,seconds"
                if kind == "encode"
                else "run,case,threads,rep,load,sort,compute,store"
            ),
            phases,
        ),
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


def abba_main(argv=None, *, kind="encode"):
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
    parser.add_argument(
        "--keep-data",
        action="store_true",
        help="retain generated comparison datasets for debugging",
    )
    args = parser.parse_args(argv)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.summarize_only:
        abba_summarize(out, kind)
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
    write_workload(kind + "_performance", out / "benchmark.do")
    with benchmark_data(out, keep=args.keep_data) as data:
        for label in abba_LABELS:
            if any(
                hashlib.sha256(Path(p).read_bytes()).hexdigest() != h for p, h in hashes.items()
            ):
                raise ValueError("Frozen build changed")
            driver = out / f"{label}.do"
            driver.write_text(
                f"""clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{builds[label[0]]}" {label} "{out}" "{args.cases}" {args.reps} "{args.threads}" {args.rows} "{args.refreps if args.refreps is not None else args.reps}" "{data}"
local rc=_rc
di "{kind.upper()}_BENCHMARK_RC=`rc'"
log close
exit, clear
"""
            )
            print("Running " + label, flush=True)
            launch_driver(driver, out / f"{label}.log")
        abba_summarize(out, kind)


audit_CASES = (
    "export import sort merge sample bsample encode decode destring ipolate split "
    "rangejoin winsor rangestat binscatter reghdfe ivreghdfe ppmlhdfe psmatch2 qreg "
    "group_numeric group_mixed group_count group_label winsor_by reghdfe_strcluster "
    "sample_by bsample_by"
).split()


def audit_summarize(out):
    records = []
    checks = {}
    for manifest_path in sorted(out.glob("*.manifest.json")):
        plan = json.loads(manifest_path.read_text())
        label = plan["label"]
        log = (out / (label + ".log")).read_text()
        if f"\nAUDIT_COMPLETE label={label}\n" not in log or not re.search(
            r"^CTOOLS_LAUNCH=[0-9a-f]+ RC=0$", log, re.M
        ):
            raise ValueError("incomplete run: " + label)
        seen = set()
        checked = set()
        for line in log.splitlines():
            if line.startswith("AUDIT,"):
                _, lab, n, case, threads, rep, seconds, rc, rows, columns = line.split(",")
                key = (int(n), case, int(rep))
                if lab != label or key in seen or int(threads) != plan["threads"]:
                    raise ValueError("invalid/duplicate record: " + line)
                seen.add(key)
                elapsed = float(seconds)
                if int(rc) or not math.isfinite(elapsed) or elapsed <= 0:
                    raise ValueError("failed measurement: " + line)
                records.append(
                    dict(
                        label=lab,
                        observations=int(n),
                        command=case,
                        threads=int(threads),
                        repetition=int(rep),
                        seconds=elapsed,
                        output_rows=int(rows),
                        output_columns=int(columns),
                    )
                )
            elif line.startswith("AUDIT_CHECK,"):
                _, lab, n, case, threads, signature, groups = line.split(",")
                key = (int(n), case, int(threads))
                if lab != label or int(threads) != plan["threads"] or key in checked:
                    raise ValueError("invalid/duplicate output check: " + line)
                checked.add(key)
                value = (signature, groups)
                if key in checks and checks[key] != value:
                    raise ValueError(f"output mismatch {key}: {checks[key]} != {value}")
                checks[key] = value
        expected = {
            (n, c, r) for n in plan["sizes"] for c in plan["cases"] for r in range(plan["reps"] + 1)
        }
        if seen != expected:
            raise ValueError(f"incomplete grid {label}: {seen ^ expected}")
        expected_checks = {(n, c, plan["threads"]) for n in plan["sizes"] for c in plan["cases"]}
        if checked != expected_checks:
            raise ValueError(f"incomplete output checks {label}: {checked ^ expected_checks}")
    if not records:
        raise ValueError("no measurements")
    with (out / "raw.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    groups = {}
    for row in records:
        if row["repetition"]:
            key = (row["label"][0], row["observations"], row["command"], row["threads"])
            groups.setdefault(key, []).append(row["seconds"])
    summary = [
        dict(
            build=k[0],
            observations=k[1],
            command=k[2],
            threads=k[3],
            trials=len(v),
            median=statistics.median(v),
            minimum=min(v),
            maximum=max(v),
        )
        for k, v in sorted(groups.items())
    ]
    with (out / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print(
            f"{row['build']} {row['observations']:>8} {row['command']:20} "
            f"t={row['threads']:2} median={row['median']:.6f} "
            f"range=[{row['minimum']:.6f}, {row['maximum']:.6f}]"
        )
    print(f"PASS: {len(records)} calls, {len(checks)} matching output signatures")


def audit_main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--label", default="b1")
    parser.add_argument("--sizes", type=int, nargs="+", default=[2000000, 5000000])
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--cases", nargs="+", choices=audit_CASES, default=audit_CASES)
    parser.add_argument("--summarize", action="store_true")
    args = parser.parse_args(argv)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.summarize:
        audit_summarize(out)
        return
    if (
        args.build is None
        or args.reps < 1
        or args.threads < 1
        or any(n < 1000 or n % 200 for n in args.sizes)
        or not re.fullmatch(r"[a-z][a-z0-9_]*", args.label)
    ):
        parser.error(
            "provide a build, positive threads/repetitions and row counts divisible by 200"
        )
    manifest_path = out / (args.label + ".manifest.json")
    if manifest_path.exists():
        parser.error("label already exists; preserve prior measurements and use a new label")
    mac = platform.system() == "Darwin"
    plugin = (
        ("ctools_mac_arm.plugin" if platform.machine() == "arm64" else "ctools_mac_x86.plugin")
        if mac
        else "ctools_linux.plugin"
    )
    frozen = out / (args.label + "_build")
    frozen.mkdir()
    for path in [args.build / plugin, *sorted(args.build.glob("*.ado"))]:
        shutil.copy2(path, frozen / path.name)
    harness = out / (args.label + ".do")
    write_workload("audit", harness)
    clock = out / "clock.plugin"
    if not clock.exists():
        flags = (
            ["-bundle", "-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
            if mac
            else ["-shared", "-fPIC", "-DSYSTEM=STUNIX", "-D_POSIX_C_SOURCE=200809L"]
        )
        subprocess.run(
            [
                "cc",
                "-std=c11",
                "-O2",
                "-I",
                str(ROOT / "src"),
                *flags,
                str(ROOT / "validation/performance/benchmark_clock.c"),
                str(ROOT / "src/stplugin.c"),
                "-o",
                str(clock),
            ],
            check=True,
        )
    paths = [*sorted(frozen.iterdir()), harness, clock]
    paths += [ROOT / "validation/bench.py", ROOT / "Makefile", ROOT / "scripts/build.mk"]
    manifest = dict(
        label=args.label,
        sizes=args.sizes,
        threads=args.threads,
        reps=args.reps,
        cases=args.cases,
        platform=platform.platform(),
        machine=platform.machine(),
        seed=95183,
        columns=40,
        files={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    run_do(
        harness,
        out / (args.label + ".log"),
        [
            frozen,
            out,
            args.label,
            " ".join(map(str, args.sizes)),
            args.threads,
            args.reps,
            " ".join(args.cases),
        ],
        markers=[f"AUDIT_COMPLETE label={args.label}"],
    )
    print("PASS " + args.label)


def transport_fixture(out):
    from validation.bench_transport import transport_main

    destination = out / "transport-plugin"
    transport_main(
        [
            "--variant",
            f'current={ROOT / "src"}',
            "--output",
            str(destination),
            "--rows",
            "1",
            "--columns",
            "1",
            "--repetitions",
            "1",
        ]
    )
    return destination / "current.plugin"


def main(argv=None):
    from functools import partial

    from validation import bench_paths, bench_profiles, bench_transport

    actions = {
        "commands": commands_main,
        "io": io_main,
        "encode": partial(abba_main, kind="encode"),
        "group": partial(abba_main, kind="group"),
        "audit": audit_main,
        "sort": bench_paths.sort_main,
        "matching": bench_paths.matching_main,
        "transport": bench_transport.transport_main,
        "transport-campaign": bench_transport.campaign_main,
        "bigdata": bench_transport.bigdata_main,
        "readme": bench_profiles.readme_main,
        "hdfe": bench_profiles.hdfe_main,
        "summarize-commands": command_main,
        "summarize-readme": bench_profiles.readme_summary_main,
        "summarize-profile": bench_profiles.profile_main,
    }
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workload", choices=[*actions, "do", "summarize-transport"])
    values = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(values[:1])
    rest = values[1:]
    if args.workload == "summarize-transport":
        if len(rest) != 2:
            parser.error("summarize-transport needs input directory and output prefix")
        return bench_transport.transport_summary_summarize(*map(Path, rest))
    if args.workload == "do":
        child = argparse.ArgumentParser(
            description="Run an opt-in Stata workload with its original arguments"
        )
        child.add_argument("name")
        child.add_argument("--output", type=Path, required=True)
        child.add_argument("arguments", nargs="*")
        request = child.parse_intermixed_args(rest)
        request.output = request.output.resolve()
        request.output.mkdir(parents=True, exist_ok=True)
        source = write_workload(request.name, request.output / "workload.do")
        return run_do(source, request.output / "workload.log", request.arguments)
    return actions[args.workload](rest)


if __name__ == "__main__":
    main()
