"""Bench profiles: shared execution with explicit workload-specific checks."""

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

from validation import ROOT
from validation.stata import run_do, stata_command, write_workload


def readme_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snapshot", required=True, type=Path)
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument(
        "--reference-dependency",
        action="append",
        type=Path,
        default=[],
        help="Additional reference source to hash, e.g. rangestat.ado (repeatable)",
    )
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--sizes", nargs="+", type=int, default=[1000000, 5000000, 20000000])
    p.add_argument("--threads", type=int, default=12)
    p.add_argument("--reps", type=int, default=4)
    a = p.parse_args(argv)
    if a.threads < 1 or a.reps < 2 or a.reps % 2 or any(n < 1000 or n % 1000 for n in a.sizes):
        p.error("Use positive threads, an even repeat count >=2, and sizes divisible by 1000.")
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    snapshot = a.snapshot.resolve()
    build = snapshot / "build"
    refs = a.reference.resolve()
    mac = platform.system() == "Darwin"
    plugin = (
        ("ctools_mac_arm.plugin" if platform.machine() == "arm64" else "ctools_mac_x86.plugin")
        if mac
        else "ctools_linux.plugin"
    )
    if not (build / plugin).is_file():
        p.error("Build the frozen snapshot first.")
    flags = ["-std=c11", "-O2", "-I", str(snapshot / "src")]
    flags += (
        ["-bundle", "-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE", "-mmacosx-version-min=11.0"]
        if mac
        else ["-shared", "-fPIC", "-DSYSTEM=STUNIX", "-D_POSIX_C_SOURCE=200809L"]
    )
    subprocess.run(
        [
            "cc",
            *flags,
            str(ROOT / "validation/performance/benchmark_clock.c"),
            str(snapshot / "src/stplugin.c"),
            "-o",
            str(out / "clock.plugin"),
        ],
        check=True,
    )
    write_workload("readme_commands", out / "benchmark.do")
    files = sorted(f for f in (snapshot / "src").rglob("*") if f.is_file())
    files += sorted(build.glob("*.ado")) + [
        build / plugin,
        snapshot / "Makefile",
        snapshot / "scripts/build.mk",
    ]
    files += sorted(refs.glob("*.ado")) + [out / "clock.plugin", out / "benchmark.do"]
    files += [f.resolve() for f in a.reference_dependency]
    files += [
        ROOT / "validation/performance/benchmark_clock.c",
        ROOT / "validation/bench_profiles.py",
    ]
    for path in [
        "/Applications/Stata/ado/base/i/ipolate.ado",
        "/Applications/Stata/ado/base/s/split.ado",
    ]:
        if Path(path).exists():
            files.append(Path(path))
    hashes = {str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    source_identity = "\n".join(
        f"{f.relative_to(snapshot)} {hashes[str(f)]}"
        for f in sorted((snapshot / "src").rglob("*"))
        if f.is_file()
    )
    manifest = dict(
        platform=platform.platform(),
        machine=platform.machine(),
        sizes=a.sizes,
        threads=a.threads,
        reps=a.reps,
        seed=95183,
        source_tree_sha256=hashlib.sha256(source_identity.encode()).hexdigest(),
        files=hashes,
    )
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for n in a.sizes:
        driver = out / f"run_{n}.do"
        driver.write_text(
            f"""clear all
set more off
log using "{out}/run_{n}.log", text replace
capture noisily do "{out}/benchmark.do" "{build}" "{refs}" "{out}" {n} {a.threads} {a.reps}
local rc=_rc
di "README_DRIVER_RC=`rc'"
log close
exit, clear
"""
        )
        print(shlex.join(stata_command(driver)))


def readme_summary_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--csv", required=True, type=Path)
    p.add_argument(
        "--evidence", type=Path, help="Compact checks and timings, without license details"
    )
    a = p.parse_args(argv)
    manifest = json.loads((a.directory / "manifest.json").read_text())
    cases = ["ipolate", "split", "rangejoin"]
    rows, evidence = [], []
    pattern = re.compile(
        r"^README_BENCH,(\w+),(\d+),(reference|ctools),(\d+),(\d+),\s*([\d.]+),(\d+)$", re.M
    )
    for n in manifest["sizes"]:
        log = (a.directory / f"run_{n}.log").read_text()
        assert re.search(r"^README_DRIVER_RC=0$", log, re.M), f"{n}: driver failed or unfinished"
        assert re.search(rf"^README_BENCHMARK_COMPLETE N={n}$", log, re.M), f"{n}: incomplete run"
        checks = re.findall(r"^README_CHECK,(\w+),(\d+),PASS,(\d+)$", log, re.M)
        outputs = {case: 2 * n - n // 500 if case == "rangejoin" else n for case in cases}
        assert sorted(checks) == sorted(
            (case, str(n), str(outputs[case])) for case in cases
        ), f"{n}: missing exact comparisons"
        seen = set()
        for case, size, method, threads, rep, seconds, output in pattern.findall(log):
            size, threads, rep, output = map(int, (size, threads, rep, output))
            seconds = float(seconds)
            assert case in cases and size == n and threads == manifest["threads"]
            assert 1 <= rep <= manifest["reps"] and output == outputs[case]
            assert math.isfinite(seconds) and seconds > 0
            key = (case, method, rep)
            assert key not in seen, f"duplicate {n}: {key}"
            seen.add(key)
            rows.append((case, n, method, threads, rep, seconds, output))
        expected = {
            (case, method, rep)
            for case in cases
            for method in ["reference", "ctools"]
            for rep in range(1, manifest["reps"] + 1)
        }
        assert seen == expected, f"{n}: missing timing records"
        evidence.extend(
            line
            for line in log.splitlines()
            if line.startswith(
                (
                    "README_CHECK,",
                    "README_BENCH,",
                    "README_BENCHMARK_COMPLETE ",
                    "README_DRIVER_RC=",
                )
            )
        )
    with a.csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "command",
                "input_rows_per_dataset",
                "implementation",
                "threads",
                "repetition",
                "seconds",
                "output_rows",
            ]
        )
        writer.writerows(rows)
    if a.evidence:
        a.evidence.write_text("\n".join(evidence) + "\n")
    print("| Command | Input rows per dataset | Reference (s) | ctools (s) | Speedup |")
    print("|---|---:|---:|---:|---:|")
    for case in cases:
        for n in manifest["sizes"]:
            times = {
                method: statistics.median(row[5] for row in rows if row[:3] == (case, n, method))
                for method in ["reference", "ctools"]
            }
            print(
                f"| `{case}` | {n:,} | {times['reference']:.3f} | {times['ctools']:.3f} | {times['reference'] / times['ctools']:.1f}× |"
            )
    print(
        f'Validated {len(rows)} timings and {len(cases) * len(manifest["sizes"])} exact full-output comparisons.'
    )


profile_COMMANDS = {
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
profile_SIZES = (1_000_000, 10_000_000, 25_000_000)


def profile_rounded(ratio, *, upper):
    """Round range endpoints outward: tenths below 2x, units below 10x, then tens."""
    rounding = math.ceil if upper else math.floor
    if ratio < 2:
        return f"{rounding(ratio * 10) / 10:.1f}"
    step = 1 if ratio < 10 else 10
    return str(rounding(ratio / step) * step)


def profile_summarize(path, sizes):
    rows = {}
    with path.open(newline="") as source:
        for row in csv.DictReader(source):
            key = (row["command"], int(row["observations"]), row["method"])
            if key in rows:
                raise ValueError(f"duplicate timing: {key}")
            rows[key] = row
    expected = {
        (cmd, n, method) for cmd in profile_COMMANDS for n in sizes for method in ("ref", "c")
    }
    if rows.keys() != expected:
        raise ValueError(f"unexpected/missing timings: {rows.keys() ^ expected}")
    output = [
        "| Command | Observations | Reference (s) | ctools (s) | Speedup |",
        "|---|---:|---:|---:|---:|",
    ]
    ranges = {}
    for cmd, label in profile_COMMANDS.items():
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
            low = profile_rounded(min(ratios), upper=False)
            high = profile_rounded(max(ratios), upper=True)
            ranges[label] = f"{low}×" if low == high else f"{low}–{high}×"
        else:
            ranges[label] = "—"
    return "\n".join(output) + "\n", ranges


def profile_main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, help="directory containing profile.log and timings.csv"
    )
    parser.add_argument(
        "--sizes", nargs="+", type=int, default=profile_SIZES, help="override only for smoke tests"
    )
    parser.add_argument("--update-readme", type=Path)
    args = parser.parse_args(argv)
    log = (args.directory / "profile.log").read_text()
    if not re.search(r"^PROFILE_COMPLETE$", log, re.MULTILINE) or not re.search(
        r"^PROFILE_EXIT rc=0$", log, re.MULTILINE
    ):
        parser.error("profile did not finish successfully; see profile.log")
    if args.update_readme and tuple(args.sizes) != profile_SIZES:
        parser.error("README publication requires the full 1M/10M/25M run")
    try:
        details, ranges = profile_summarize(args.directory / "timings.csv", args.sizes)
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


hdfe_OUT = ROOT / "temp/hdfe_datasets"
hdfe_DOC = ROOT / "docs/performance/HDFE_DATASETS.md"
hdfe_EVIDENCE = ROOT / "docs/benchmarks"
hdfe_MARKER = "HDFE_DATASETS_COMPLETE"
# Order of https://www.scorreia.com/data/hdfe/index.html.
hdfe_DATASETS = [
    "credit2",
    "credit",
    "soccer",
    "synthetic-complete",
    "synthetic-uniform-easy",
    "synthetic-uniform-hard",
    "synthetic-uniform-harder",
    "synthetic-assortative",
    "synthetic-zigzag",
    "enron",
    "github",
    "patents",
    "workers",
    "schools",
    "directors",
]
hdfe_LABELS = {
    "specification": "Specification",
    "settings": "Tested settings",
    "n_full": "Rows",
    "n": "N used",
    "singletons": "Singletons dropped",
    "groups": "Mobility groups",
    "df_a": "Absorbed df",
    "df_a_reghdfe": "Absorbed df, reghdfe",
    "df_a_creghdfe": "Absorbed df, creghdfe",
    "iterations": "reghdfe iterations",
    "sf_b": "Digits, b",
    "sf_v": "Digits, V",
    "sf_stats": "Digits, fit",
    "sf_v_default": "Digits, V (default settings)",
    "t_reghdfe": "reghdfe (s)",
    "t_creghdfe": "creghdfe (s)",
    "result": "Result",
}
hdfe_TEXT = ("specification", "settings", "result")


def hdfe_read(name):
    with open(hdfe_OUT / (name + ".csv"), newline="") as handle:
        rows = list(csv.DictReader(handle))
    return sorted(
        rows, key=lambda row: hdfe_DATASETS.index(row["dataset"])
    )  # stable: keeps run order


def hdfe_cases(results):
    """One entry per specification: its default-settings run and the run the test asserts."""
    found = {}
    for row in results:
        key = (row["section"], row["dataset"], row["varlist"], row["absorb"], row["vce"])
        case = found.setdefault(key, {"default": row})
        case["tested"] = row  # the converged run, when there is one, follows the default run
    return list(found.values())


def hdfe_number(value, decimals=None):
    if value in ("", "."):
        return ""
    return f"{round(float(value)):,}" if decimals is None else f"{float(value):.{decimals}f}"


def hdfe_result(row):
    if row["status"] != "fail":
        return row["status"]
    return "fail: absorbed df" if "df_a:" in row["differences"] else "fail"


def hdfe_cell(case, column):
    default, tested = case["default"], case["tested"]
    if column == "specification":
        regressors = "" if tested["varlist"] == "y x1 x2" else tested["varlist"] + ", "
        vce = f" vce({tested['vce']})" if tested["vce"] else ""
        return f"`{regressors}absorb({tested['absorb']}){vce}`"
    if column == "settings":
        return tested["settings"]
    if column == "result":
        return hdfe_result(tested)
    if column == "sf_v_default":
        return hdfe_number(default["sf_v"], 1)
    if column.startswith("sf_"):
        return hdfe_number(tested[column], 1)
    if column.startswith("t_"):  # timings always describe the default-settings run
        return hdfe_number(default[column], 2)
    return hdfe_number(tested["df_a" if column == "df_a_reghdfe" else column])


def hdfe_table(header, rows, text=()):
    align = ["---" if index == 0 or index in text else "---:" for index in range(len(header))]
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join(align) + "|"]
    return "\n".join(lines + ["| " + " | ".join(row) + " |" for row in rows])


def hdfe_section(found, name, columns):
    rows = [
        [case["tested"]["dataset"]] + [hdfe_cell(case, column) for column in columns]
        for case in found
        if case["tested"]["section"] == name
    ]
    text = [index + 1 for index, column in enumerate(columns) if column in hdfe_TEXT]
    return hdfe_table(["Dataset"] + [hdfe_LABELS[column] for column in columns], rows, text)


def hdfe_environment(log):
    stata = dict(re.findall(r"^HDFE_ENV (\w+)=\s*(.*)$", log, re.MULTILINE))
    mac = platform.system() == "Darwin"
    # The plugin Stata actually loaded, as the suite reported it.
    plugin = ROOT / stata.get("plugin", "build/missing")
    published = plugin.resolve().parent == (ROOT / "build").resolve()
    source = f"`build/{plugin.name}`" if published else f"`{plugin.name}` built outside `build/`"

    def shell(*command):
        try:
            return subprocess.run(
                command, cwd=ROOT, capture_output=True, text=True, check=True
            ).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return "unknown"

    cpu = shell("sysctl", "-n", "machdep.cpu.brand_string") if mac else platform.processor()
    revision = shell("git", "rev-parse", "--short", "HEAD")
    if revision != "unknown" and shell("git", "status", "--porcelain", "--untracked-files=no"):
        revision += " with uncommitted changes"
    digest = hashlib.sha256(plugin.read_bytes()).hexdigest() if plugin.exists() else "unknown"
    version = lambda name: stata.get(name, "unknown").removeprefix("*! ")
    return "\n".join(
        [
            f"- Run date: {stata.get('date', 'unknown')}",
            (
                f"- Machine: {cpu}, {os.cpu_count()} logical cores, {platform.system()} {platform.release()} "
                f"({platform.machine()})"
            ),
            f"- Stata: {stata.get('stata', 'unknown')}, licensed for {stata.get('processors', '?')} processors",
            f"- Reference: reghdfe `{version('reghdfe')}`, ftools `{version('ftools')}`",
            f"- ctools: creghdfe `{version('creghdfe')}`, source revision {revision}",
            f"- Plugin: {source}, SHA-256 `{digest}`, default thread count",
        ]
    )


def hdfe_blocks(log):
    results, oracle = hdfe_read("results"), hdfe_read("oracle")
    found = hdfe_cases(results)
    tested = [case["tested"] for case in found]
    summary = re.search(r"^\s+Passed: (\d+)\n\s+Failed: (\d+)", log, re.MULTILINE)
    passed, failed = summary.groups() if summary else ("?", "?")
    by_status = lambda status: sum(row["status"] == status for row in tested)
    absorbed = sum(hdfe_result(row) == "fail: absorbed df" for row in tested)
    default = [case["default"] for case in found]
    seconds = lambda column: sum(
        float(row[column]) for row in default if row[column] not in ("", ".")
    )
    least = lambda rows, column: min(
        float(row[column]) for row in rows if row[column] not in ("", ".")
    )
    passing = [row for row in tested if row["status"] == "pass"]
    on_absorbed = f" ({absorbed} of them on absorbed degrees of freedom)" if absorbed else ""
    exact = hdfe_table(
        [
            "Dataset",
            "Specification",
            "Run",
            "Residual df",
            "Exact residual df",
            "Digits, b",
            "Digits, V",
            "Digits, RSS",
            "Result",
        ],
        [
            [
                row["dataset"],
                f"`{row['varlist']}, absorb({row['absorb']})`",
                f"`{row['run'].strip()}`",
                hdfe_number(row["df_r"]),
                hdfe_number(row["df_r_exact"]),
                hdfe_number(row["sf_b"], 1),
                hdfe_number(row["sf_v"], 1),
                hdfe_number(row["sf_rss"], 1),
                row["status"],
            ]
            for row in oracle
        ],
        text=(1, 2, 8),
    )
    return {
        "environment": hdfe_environment(log),
        "totals": "\n".join(
            [
                f"- Checks: {passed} passed, {failed} failed.",
                (
                    f"- Regression comparisons: {len(tested)} specifications, {by_status('pass')} pass, "
                    f"{by_status('fail')} fail{on_absorbed}, {by_status('error')} errors."
                ),
                (
                    f"- Smallest agreement: coefficients {least(default, 'sf_b'):.1f} digits over every "
                    f"default-settings run; among the passing comparisons, covariance "
                    f"{least(passing, 'sf_v'):.1f} digits and fit statistics {least(passing, 'sf_stats'):.1f} digits."
                ),
                (
                    f"- Default-settings wall clock over those {len(default)} specifications, one run each: "
                    f"reghdfe {seconds('t_reghdfe'):.1f} s, creghdfe {seconds('t_creghdfe'):.1f} s."
                ),
            ]
        ),
        "oneway": hdfe_section(
            found,
            "oneway",
            [
                "specification",
                "n",
                "singletons",
                "df_a",
                "sf_b",
                "sf_v",
                "sf_stats",
                "t_reghdfe",
                "t_creghdfe",
                "result",
            ],
        ),
        "twoway": hdfe_section(
            found,
            "twoway",
            [
                "n_full",
                "singletons",
                "n",
                "df_a",
                "groups",
                "iterations",
                "sf_b",
                "sf_v",
                "sf_stats",
                "t_reghdfe",
                "t_creghdfe",
                "result",
            ],
        ),
        "vce": hdfe_section(
            found,
            "vce",
            [
                "specification",
                "df_a_reghdfe",
                "df_a_creghdfe",
                "sf_v_default",
                "sf_b",
                "sf_v",
                "t_reghdfe",
                "t_creghdfe",
                "result",
            ],
        ),
        "regroup": hdfe_section(
            found,
            "regroup",
            [
                "specification",
                "n",
                "singletons",
                "df_a",
                "groups",
                "sf_b",
                "sf_v",
                "sf_stats",
                "t_reghdfe",
                "t_creghdfe",
                "result",
            ],
        ),
        "multiway": hdfe_section(
            found,
            "multiway",
            [
                "specification",
                "settings",
                "n",
                "df_a_reghdfe",
                "df_a_creghdfe",
                "sf_b",
                "sf_v",
                "t_reghdfe",
                "t_creghdfe",
                "result",
            ],
        ),
        "oracle": exact,
    }


def hdfe_write_doc(log):
    text = hdfe_DOC.read_text()
    for name, body in hdfe_blocks(log).items():
        begin, end = f"<!-- BEGIN GENERATED: {name} -->\n", f"<!-- END GENERATED: {name} -->"
        if text.count(begin) != 1 or text.count(end) != 1:
            raise ValueError(f"{hdfe_DOC.name}: missing generated block {name}")
        head, rest = text.split(begin)
        text = head + begin + body + "\n" + end + rest.split(end)[1]
    hdfe_DOC.write_text(text)
    for name in ("results", "oracle"):
        shutil.copyfile(hdfe_OUT / (name + ".csv"), hdfe_EVIDENCE / f"hdfe_datasets_{name}.csv")


def hdfe_main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-doc", action="store_true", help="regenerate the tables in " + hdfe_DOC.name
    )
    parser.add_argument(
        "--render-only",
        action="store_true",
        help="regenerate the tables from the last run without starting Stata",
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        help="test the plugin and ado files in this directory instead of build/",
    )
    args = parser.parse_args(argv)
    log = hdfe_OUT / "validation.log"
    if args.render_only:
        hdfe_write_doc(log.read_text(errors="replace"))
        return
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/fetch_validation_data.py"),
            "--collection",
            "hdfe",
            "--check-only",
        ],
        cwd=ROOT,
        check=True,
    )
    hdfe_OUT.mkdir(parents=True, exist_ok=True)
    for name in ("results", "oracle"):  # tables are written only by a run that reaches the end
        (hdfe_OUT / (name + ".csv")).unlink(missing_ok=True)
    failure = None
    try:
        build = [args.build_dir.resolve()] if args.build_dir else []
        run_do(
            write_workload("hdfe_datasets", hdfe_OUT / "workload.do"),
            log,
            [hdfe_OUT, *build],
            [hdfe_MARKER],
        )
    except (ValueError, subprocess.CalledProcessError) as error:
        failure = error
    # Failed comparisons still leave complete tables, and those are what the document must show.
    if args.write_doc and (hdfe_OUT / "results.csv").exists():
        hdfe_write_doc(log.read_text(errors="replace"))
    if failure:
        raise SystemExit(f"HDFE dataset validation failed: {failure}; see {log}")
    print(f"HDFE dataset validation passed. Log: {log}")
