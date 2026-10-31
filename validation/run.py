"""Run correctness checks. Names select individual commands or native checks."""

if __package__ in (None, ""):
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import argparse
import importlib
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

from validation import ROOT

# Defaults previously exercised by the Linux correctness gate. Individual tests
# retain additional compiler modes and runtime variants inside their functions.
ASAN = {
    "test_newcommands_native",
    "test_sep24_native",
    "test_counting_sort_native",
    "test_csort_ordered_native",
    "test_msd_sort_native",
    "test_sort_keys_native",
    "test_split_tokens_native",
    "test_transport_native",
    "test_transport_scheduling",
    "test_contracts_native",
    "test_thread_policy_native",
}


def native_cases():
    found = {}
    for path in sorted((ROOT / "validation").glob("native_*.py")):
        module = importlib.import_module("validation." + path.stem)
        for name, value in vars(module).items():
            if name.startswith("test_") and callable(value):
                if name in found:
                    raise ValueError("Duplicate native check: " + name)
                found[name] = module.__name__
    return found


def run_native(names, out):
    available = native_cases()
    failures = []
    for name in names:
        env = os.environ.copy()
        if platform.system() == "Linux":
            if name in ASAN:
                env.setdefault("CTOOLS_SANITIZERS", "address,undefined")
            elif name == "test_memory_safety_native":
                env.setdefault("CTOOLS_SANITIZERS", "address,undefined,float-cast-overflow")
        log = out / (name + ".log")
        with log.open("w") as handle:
            code = f"from {available[name]} import {name}; {name}()"
            result = subprocess.run(
                [sys.executable, "-c", code],
                cwd=ROOT,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            failures.append(name)
            print(f"FAIL {name}: {log}", flush=True)
            print(log.read_text(errors="replace")[-6000:])
        else:
            print(f"PASS {name}", flush=True)
    return failures


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("names", nargs="*", help="command, Stata case, or native test name")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--native", action="store_true", help="native tests and repository checks only"
    )
    mode.add_argument("--stata", action="store_true", help="Stata tests only")
    mode.add_argument("--checks", action="store_true", help="repository/harness checks only")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--output", type=Path)
    args, extra = parser.parse_known_args(argv)
    if args.checks and args.names:
        from validation import checks

        return checks.main([*args.names, *extra])
    if extra:
        parser.error("Unknown arguments: " + " ".join(extra))
    from validation import stata

    available, nested = stata.cases()
    native = native_cases()
    if args.list:
        print("Stata: " + " ".join(sorted((set(available) - nested) | {"option_coverage"})))
        print("Native: " + " ".join(native))
        return
    unknown = set(args.names) - available.keys() - native.keys() - {"option_coverage"}
    if unknown:
        parser.error("Unknown tests: " + ", ".join(sorted(unknown)))
    if args.names and args.native and not set(args.names) <= native.keys():
        parser.error("--native accepts native test names only")
    if args.names and args.stata and set(args.names) & native.keys():
        parser.error("--stata accepts Stata case or command names only")
    if args.output:
        out = args.output.resolve()
        out.mkdir(parents=True, exist_ok=True)
    else:
        parent = ROOT / "temp/validation"
        parent.mkdir(parents=True, exist_ok=True)
        out = Path(tempfile.mkdtemp(prefix="run-", dir=parent))
    print("Logs: " + str(out), flush=True)
    failures = []
    if not args.stata and not args.checks:
        selected = [n for n in native if not args.names or n in args.names]
        failures += run_native(selected, out)
    if not args.stata and not args.names:
        from validation import checks

        failures += checks.run(out)
    if not args.native and not args.checks and (not args.names or set(args.names) - native.keys()):
        failures += run_stata(args.names, out)
    if failures:
        raise SystemExit("Failed: " + ", ".join(failures))
    print("All selected checks passed. Logs: " + str(out))


def run_stata(selected, out):
    from validation import stata

    available, nested = stata.cases()
    chosen = {
        name
        for name, (path, _) in available.items()
        if name in selected or (name not in nested and (not selected or path.stem in selected))
    }
    order = [name for name in STATA_ORDER if name in chosen]
    order += sorted(chosen - set(order))
    failures = []
    # Check offline inputs before launching the full reference suite.
    subprocess.run(
        [sys.executable, "scripts/fetch_validation_data.py", "--check-only"], cwd=ROOT, check=True
    )
    masters = [name for name in order if available[name][1]]
    if masters:
        try:
            stata.run_groups(masters, out)
        except (ValueError, subprocess.CalledProcessError) as error:
            print(error, flush=True)
            failures.append("Stata command/contract suites")
    for name in order:
        if name in masters:
            continue
        print("Running " + name, flush=True)
        try:
            if name in (
                "io_parity",
                "io_formats",
                "io_delimited_options",
                "io_excel_options",
                "oct21",
                "oct22",
            ):
                stata.run_fixture(name, out)
            elif name == "transport_strl":
                from validation.bench import transport_fixture

                stata.run_groups([name], out, fixture_args={name: [transport_fixture(out)]})
            else:
                stata.run_groups([name], out)
        except (ValueError, subprocess.CalledProcessError) as error:
            print(error, flush=True)
            failures.append(name)
    commands = set(selected) & {p.stem for p in (ROOT / "build").glob("*.sthlp")}
    if not selected or "option_coverage" in selected or commands:
        try:
            from validation.options import run

            run(
                out / "options",
                commands if selected and "option_coverage" not in selected else None,
            )
        except (ValueError, subprocess.CalledProcessError, SystemExit) as error:
            print(error, flush=True)
            failures.append("option_coverage")
    return failures


STATA_ORDER = (
    "csort",
    "csort_identity",
    "cmerge",
    "cimport",
    "cexport",
    "creghdfe",
    "cqreg",
    "civreghdfe",
    "cdecode",
    "cencode",
    "cdestring",
    "csample",
    "cbsample",
    "cbinscatter",
    "cpsmatch",
    "crangestat",
    "cwinsor",
    "cppmlhdfe",
    "cinterpolate",
    "csplit",
    "crangejoin",
    "audit_p1",
    "audit_p2",
    "sep22",
    "sep24",
    "oct27",
    "memory_safety",
    "command_optimizations",
    "contracts",
    "threads",
    "group_windows",
    "option_contracts",
    "encode_labels",
    "winsor_transport",
    "io_parity",
    "io_formats",
    "io_delimited_options",
    "io_excel_options",
    "cqreg_performance",
    "bigdata",
    "transport_strl",
    "data_io",
    "refactor_timing",
    "coverage_helpers",
    "option_coverage",
    "oct21",
    "performance_paths",
    "oct22",
    "oct22_reference",
    "winsor_rounding",
    "cpsmatch_statistics",
)

if __name__ == "__main__":
    main()
