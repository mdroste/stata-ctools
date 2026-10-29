"""creghdfe vs reghdfe on the HDFE benchmark collection; opt-in, not part of the offline gate.

Verifies the pinned data cache, runs validate_hdfe_datasets.do through the stata
shell alias, and with --write-doc copies the result tables into docs/benchmarks
and regenerates the marked tables of docs/performance/HDFE_DATASETS.md.
"""

import argparse
import csv
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

from validation.support.stata_runner import ROOT, run_do

OUT = ROOT / "temp/hdfe_datasets"
DOC = ROOT / "docs/performance/HDFE_DATASETS.md"
EVIDENCE = ROOT / "docs/benchmarks"
MARKER = "HDFE_DATASETS_COMPLETE"
# Order of https://www.scorreia.com/data/hdfe/index.html.
DATASETS = [
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
LABELS = {
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
TEXT = ("specification", "settings", "result")


def read(name):
    with open(OUT / (name + ".csv"), newline="") as handle:
        rows = list(csv.DictReader(handle))
    return sorted(rows, key=lambda row: DATASETS.index(row["dataset"]))  # stable: keeps run order


def cases(results):
    """One entry per specification: its default-settings run and the run the test asserts."""
    found = {}
    for row in results:
        key = (row["section"], row["dataset"], row["varlist"], row["absorb"], row["vce"])
        case = found.setdefault(key, {"default": row})
        case["tested"] = row  # the converged run, when there is one, follows the default run
    return list(found.values())


def number(value, decimals=None):
    if value in ("", "."):
        return ""
    return f"{round(float(value)):,}" if decimals is None else f"{float(value):.{decimals}f}"


def result(row):
    if row["status"] != "fail":
        return row["status"]
    return "fail: absorbed df" if "df_a:" in row["differences"] else "fail"


def cell(case, column):
    default, tested = case["default"], case["tested"]
    if column == "specification":
        regressors = "" if tested["varlist"] == "y x1 x2" else tested["varlist"] + ", "
        vce = f" vce({tested['vce']})" if tested["vce"] else ""
        return f"`{regressors}absorb({tested['absorb']}){vce}`"
    if column == "settings":
        return tested["settings"]
    if column == "result":
        return result(tested)
    if column == "sf_v_default":
        return number(default["sf_v"], 1)
    if column.startswith("sf_"):
        return number(tested[column], 1)
    if column.startswith("t_"):  # timings always describe the default-settings run
        return number(default[column], 2)
    return number(tested["df_a" if column == "df_a_reghdfe" else column])


def table(header, rows, text=()):
    align = ["---" if index == 0 or index in text else "---:" for index in range(len(header))]
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join(align) + "|"]
    return "\n".join(lines + ["| " + " | ".join(row) + " |" for row in rows])


def section(found, name, columns):
    rows = [
        [case["tested"]["dataset"]] + [cell(case, column) for column in columns]
        for case in found
        if case["tested"]["section"] == name
    ]
    text = [index + 1 for index, column in enumerate(columns) if column in TEXT]
    return table(["Dataset"] + [LABELS[column] for column in columns], rows, text)


def environment(log):
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


def blocks(log):
    results, oracle = read("results"), read("oracle")
    found = cases(results)
    tested = [case["tested"] for case in found]
    summary = re.search(r"^\s+Passed: (\d+)\n\s+Failed: (\d+)", log, re.MULTILINE)
    passed, failed = summary.groups() if summary else ("?", "?")
    by_status = lambda status: sum(row["status"] == status for row in tested)
    absorbed = sum(result(row) == "fail: absorbed df" for row in tested)
    default = [case["default"] for case in found]
    seconds = lambda column: sum(
        float(row[column]) for row in default if row[column] not in ("", ".")
    )
    least = lambda rows, column: min(
        float(row[column]) for row in rows if row[column] not in ("", ".")
    )
    passing = [row for row in tested if row["status"] == "pass"]
    on_absorbed = f" ({absorbed} of them on absorbed degrees of freedom)" if absorbed else ""
    exact = table(
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
                number(row["df_r"]),
                number(row["df_r_exact"]),
                number(row["sf_b"], 1),
                number(row["sf_v"], 1),
                number(row["sf_rss"], 1),
                row["status"],
            ]
            for row in oracle
        ],
        text=(1, 2, 8),
    )
    return {
        "environment": environment(log),
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
        "oneway": section(
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
        "twoway": section(
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
        "vce": section(
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
        "regroup": section(
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
        "multiway": section(
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


def write_doc(log):
    text = DOC.read_text()
    for name, body in blocks(log).items():
        begin, end = f"<!-- BEGIN GENERATED: {name} -->\n", f"<!-- END GENERATED: {name} -->"
        if text.count(begin) != 1 or text.count(end) != 1:
            raise ValueError(f"{DOC.name}: missing generated block {name}")
        head, rest = text.split(begin)
        text = head + begin + body + "\n" + end + rest.split(end)[1]
    DOC.write_text(text)
    for name in ("results", "oracle"):
        shutil.copyfile(OUT / (name + ".csv"), EVIDENCE / f"hdfe_datasets_{name}.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-doc", action="store_true", help="regenerate the tables in " + DOC.name
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
    args = parser.parse_args()
    log = OUT / "validation.log"
    if args.render_only:
        write_doc(log.read_text(errors="replace"))
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
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ("results", "oracle"):  # tables are written only by a run that reaches the end
        (OUT / (name + ".csv")).unlink(missing_ok=True)
    failure = None
    try:
        build = [args.build_dir.resolve()] if args.build_dir else []
        run_do(ROOT / "validation/performance/validate_hdfe_datasets.do", log, [OUT, *build], [MARKER])
    except (ValueError, subprocess.CalledProcessError) as error:
        failure = error
    # Failed comparisons still leave complete tables, and those are what the document must show.
    if args.write_doc and (OUT / "results.csv").exists():
        write_doc(log.read_text(errors="replace"))
    if failure:
        raise SystemExit(f"HDFE dataset validation failed: {failure}; see {log}")
    print(f"HDFE dataset validation passed. Log: {log}")


if __name__ == "__main__":
    main()
