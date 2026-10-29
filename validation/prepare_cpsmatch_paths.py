"""Serial fresh-process path benchmarks and separate frozen-output comparisons.

Build directories must be frozen; variants are supplied in execution order.
All observations are retained. Use --parity for untimed output comparisons;
never compare its timings because save/merge affects Stata allocation history.
"""

import argparse
import csv
import hashlib
import json
import math
import subprocess
import threading
import time
from collections import defaultdict
from pathlib import Path
from statistics import median

from stata_runner import ROOT, launch_driver

CASES = {
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
STAGES = ["wall", "load", "setup", "sort", "match", "stats", "store", "total"]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def quote(value):
    value = str(value)
    if any(c in value for c in ['"', "`", "\n"]):
        raise ValueError("unsupported quote")
    return '"' + value + '"'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant", action="append", required=True, help="label=build directory")
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
    a = p.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    variants = {k: Path(v).resolve() for k, v in (x.split("=", 1) for x in a.variant)}
    cases = {k: int(v) for k, v in (x.split(":") for x in a.cases)} if a.cases else CASES
    files = {
        k: {p.name: digest(p) for p in [*v.glob("*.ado"), v / "ctools_mac_arm.plugin"]}
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
    source = ROOT / "validation/benchmark_cpsmatch_paths.do"
    frozen = out / "benchmark.do"
    frozen.write_bytes(source.read_bytes())
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
            str(ROOT / "validation/benchmark_clock.c"),
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
            assert all(digest(build / name) == h for name, h in files[label].items())
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
                + quote(log)
                + ", text replace\ncapture noisily do "
                + quote(frozen)
                + " "
                + " ".join(map(quote, args))
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
                        **dict(zip(STAGES, values)),
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
            **{s: median(r[s] for r in v) for s in STAGES},
        )
        for k, v in groups.items()
    ]
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("CPSMATCH_PATH_CAMPAIGN_COMPLETE samples=" + str(len(samples)), flush=True)


if __name__ == "__main__":
    main()
