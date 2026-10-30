"""The single Stata launch boundary: configured executable or local zsh alias.

Stata batch exit status alone does not report do-file failure. Each driver must
save _rc before logging/cleanup; a nonce rejects stale or interrupted output.
"""

import os
import re
import shlex
import subprocess
import tempfile
import uuid
from contextlib import contextmanager
from pathlib import Path

from validation import ROOT


@contextmanager
def benchmark_data(out, *, keep=False):
    """Share generated data across benchmark sessions, then remove it even on failure."""
    if keep:
        yield Path(out)
    else:
        with tempfile.TemporaryDirectory(prefix="data-", dir=out) as tmp:
            yield Path(tmp)


@contextmanager
def io_workspace():
    """Give legacy relative-path I/O tests disposable data, retaining external logs."""
    with tempfile.TemporaryDirectory(prefix="ctools-io-") as tmp:
        workspace = Path(tmp)
        for name in ("build", "validation"):
            (workspace / name).symlink_to(ROOT / name, target_is_directory=True)
        # cio_export_test shares these output paths across the option suites.
        for name in ("native", "result"):
            (workspace / "temp/io_allformats" / name).mkdir(parents=True)
        yield workspace


def stata_command(driver):
    executable = os.environ.get("CTOOLS_STATA_EXECUTABLE")
    driver = str(Path(driver).resolve())
    if executable:
        # A path, not a shell command: spaces need no shell quoting.
        batch = ["/e", "do"] if os.name == "nt" else ["-q", "-b", "do"]
        return [executable, *batch, driver]
    if os.name == "nt":
        raise ValueError("Set CTOOLS_STATA_EXECUTABLE to the licensed Stata executable")
    # Retain the existing local alias workflow; CI configures an executable.
    return ["/bin/zsh", "-lic", "stata -q -b do " + shlex.quote(driver)]


def validate_completion(text, nonce):
    matches = re.findall(r"^CTOOLS_LAUNCH=" + re.escape(nonce) + r" RC=(\d+)$", text, re.M)
    if matches != ["0"]:
        raise ValueError("Stata failed or did not finish this driver: " + repr(matches))


def quote(value):
    """Quote a Stata argument, including paths with spaces or ordinary quotes."""
    value = str(value)
    if "\n" in value or "\r" in value or "`" in value:
        raise ValueError("unsupported Stata path characters")
    return '`"' + value + "\"'"


def launch_driver(driver, log, *, rc_local="rc", cwd=ROOT):
    source = Path(driver).read_text()
    if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", rc_local):
        raise ValueError("invalid return-code local")
    if not source.rstrip().endswith("exit, clear"):
        raise ValueError("batch driver must end with exit, clear")
    if "\nlog close\n" not in source:
        raise ValueError("batch driver must close its log after reporting _rc")
    nonce = uuid.uuid4().hex
    source = source.replace(
        "\nlog close\n", f'\ndi "CTOOLS_LAUNCH={nonce} RC=`{rc_local}\'"\nlog close\n'
    )
    log = Path(log)
    log.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix="ctools-stata-") as tmp:
        executable = Path(tmp) / "driver.do"
        # Stata opens driver.log before reading the do-file. Start in scratch
        # space, then restore the requested working directory for suite paths.
        executable.write_text(f"cd {quote(Path(cwd).resolve())}\n" + source)
        subprocess.run(stata_command(executable), cwd=tmp, check=True)
    text = log.read_text(errors="replace") if log.exists() else ""
    try:
        validate_completion(text, nonce)
    except ValueError:
        print(text[-20000:])
        raise
    return text


def run_do(path, log, args=(), markers=(), *, cwd=ROOT):
    source = f"""clear all
set more off
set linesize 255
capture log close _all
log using {quote(log)}, text replace
capture noisily do {quote(Path(path).resolve())} {' '.join(quote(a) for a in args)}
local rc = _rc
log close
exit, clear
"""
    with tempfile.TemporaryDirectory(prefix="ctools-do-") as tmp:
        driver = Path(tmp) / "driver.do"
        driver.write_text(source)
        text = launch_driver(driver, log, cwd=cwd)
    for marker in markers:
        if marker not in text.splitlines():
            raise ValueError(f"missing completion marker: {marker}")
    return text


def workload_source(name):
    """Freeze ordinary Stata source with the selected workload as its first argument."""
    matches = [
        p
        for p in (ROOT / "validation/performance").glob("*.do")
        if f'if "`__workload\'" == "{name}" {{' in p.read_text()
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one benchmark workload {name}, found {matches}")
    return 'local __workload "' + name + '"\n' + matches[0].read_text()


def write_workload(name, destination):
    destination.write_text(workload_source(name))
    return destination


def cases():
    """Discover named do-file sections and their completion/counter contracts."""
    found = {}
    called = set()
    pattern = re.compile(r'^if "`__case\'" == "([a-z0-9_]+)" \{\n', re.M)
    for path in sorted((ROOT / "validation/stata").glob("*.do")):
        source = path.read_text()
        sections = list(pattern.finditer(source))
        for i, match in enumerate(sections):
            name = match[1]
            body = source[
                match.end() : sections[i + 1].start() if i + 1 < len(sections) else len(source)
            ]
            if name in found:
                raise ValueError("Duplicate Stata case: " + name)
            if f'global CTOOLS_CASE_COMPLETE "{name}"' not in body:
                raise ValueError("Missing completion assignment: " + name)
            found[name] = (path, f'global CTOOLS_COMPONENT_COMPLETE "{name}"' in body)
        called.update(re.findall(r'do "?validation/stata/\w+\.do"? (\w+)', source))
    public = {p.stem for p in (ROOT / "build").glob("*.sthlp")} - {"ctools"}
    missing = public - found.keys()
    if missing:
        raise ValueError(f"Missing public-command tests: {sorted(missing)}")
    return found, called


def validate_results(text, expected):
    found = {}
    for name, passed, failed, skipped, rc, complete in re.findall(
        r"^TEST_RESULT (\w+) PASSED=(\d+) FAILED=(\d+) SKIPPED=(\d+) RC=(\d+) COMPLETE=(\w*)$",
        text,
        re.M,
    ):
        if name in found:
            raise ValueError("Duplicate result: " + name)
        found[name] = (int(passed), int(failed), int(skipped), int(rc), complete)
    if found.keys() != set(expected):
        raise ValueError(
            f"Incomplete results: missing={set(expected)-found.keys()}, extra={found.keys()-set(expected)}"
        )
    for name, counted in expected.items():
        passed, failed, skipped, rc, complete = found[name]
        if rc or complete != name or (counted and (not passed or failed or skipped)):
            raise ValueError(f"Failed or incomplete {name}: {found[name]}")
    return found


def run_groups(names, out, *, fixture_args=None, cwd=ROOT):
    available, _ = cases()
    expected = {name: available[name][1] for name in names}
    fixture_args = fixture_args or {}
    log = out / ("stata.log" if len(names) > 1 else names[0] + ".log")
    lines = [
        "clear all",
        "set more off",
        "set linesize 255",
        "capture log close _all",
        f"log using {quote(log)}, text replace",
        "local rc = 0",
    ]
    revision = os.environ.get("CTOOLS_EXPECTED_REVISION", "")
    expected_platform = os.environ.get("CTOOLS_EXPECTED_PLATFORM", "")
    if expected_platform:
        from scripts.stage_package import native_platform

        if native_platform() != expected_platform or not os.environ.get("CTOOLS_STATA_EXECUTABLE"):
            raise ValueError(
                "Licensed runner platform/executable does not match the candidate artifact"
            )
    if revision:
        if not re.fullmatch("[0-9a-f]{40}", revision):
            raise ValueError("Expected revision must be a full commit SHA")
        lines += [
            "adopath ++ build",
            "ctools, version",
            f'assert "`r(plugin_revision)\'" == "{revision}"',
        ]
    if any(expected.values()):
        lines += [
            "foreach cmd in reghdfe ivreghdfe ppmlhdfe psmatch2 binscatter rangestat rangejoin winsor2 gstats {",
            "    capture which `cmd'",
            "    if _rc {",
            '        di as error "Required validation reference missing: `cmd\'"',
            "        local rc = 9",
            "    }",
            "}",
        ]
    lines += ["if !`rc' {"]
    for name in names:
        path, counted = available[name]
        args = [name, *fixture_args.get(name, ())]
        lines += [
            "global TESTS_PASSED = 0",
            "global TESTS_FAILED = 0",
            "global TESTS_SKIPPED = 0",
            'global CTOOLS_CASE_COMPLETE ""',
            f"capture noisily do {quote(path)} " + " ".join(map(quote, args)),
            "local case_rc = _rc",
            f'di "TEST_RESULT {name} PASSED=$TESTS_PASSED FAILED=$TESTS_FAILED SKIPPED=$TESTS_SKIPPED RC=`case_rc\' COMPLETE=$CTOOLS_CASE_COMPLETE"',
            f'if `case_rc\' | "$CTOOLS_CASE_COMPLETE" != "{name}" local rc = 9',
        ]
        if counted:
            lines += ["if $TESTS_PASSED == 0 | $TESTS_FAILED | $TESTS_SKIPPED local rc = 9"]
    lines += ["}", "log close", "exit, clear"]
    with tempfile.TemporaryDirectory(prefix="driver-", dir=out) as temporary:
        driver = Path(temporary) / "run.do"
        driver.write_text("\n".join(lines) + "\n")
        text = launch_driver(driver, log, cwd=cwd)
    results = validate_results(text, expected)
    for name, values in results.items():
        print(f"PASS {name}" + (f": {values[0]} checks" if expected[name] else ""), flush=True)
    return results


def run_fixture(name, out):
    from validation import fixtures

    with io_workspace() as workspace:
        folder = workspace / "temp" / {"io_formats": "io_allformats"}.get(name, name)
        folder.mkdir(parents=True, exist_ok=True)
        plugin = os.environ.get("CTOOLS_TEST_PLUGIN_DIR", "")
        if name == "io_parity":
            fixtures.parity_fixtures(folder)
            args = [plugin]
        elif name == "io_formats":
            for sub in ("native", "result"):
                (folder / sub).mkdir(exist_ok=True)
            fixtures.shapefiles(folder)
            fixtures.codec_fixtures(folder)
            fixtures.failure_fixtures(folder)
            args = [plugin]
        elif name == "io_delimited_options":
            fixtures.delimited_fixtures(folder)
            args = [folder]
        elif name == "io_excel_options":
            fixtures.excel_fixtures(folder)
            fixtures.formula_fixtures(folder)
            for ext in ("xls", "xlsx"):
                (folder / f"{ext}_invalid.{ext}").write_bytes(b"invalid file header\n")
                (folder / f"{ext}_truncated.{ext}").write_bytes(
                    (folder / f"display_formats.{ext}").read_bytes()[:16]
                )
            args = [folder]
        elif name == "oct21":
            fixtures.integrity_fixtures(folder)
            args = [folder, plugin]
        elif name == "oct22":
            fixtures.regression_fixtures(folder)
            args = [folder]
        else:
            raise ValueError("Unknown fixture case: " + name)
        result = run_groups([name], out, fixture_args={name: args}, cwd=workspace)
        if name == "io_excel_options":
            print(
                f"PASS {fixtures.compare_formula_outputs(folder)} formula/calculation-chain comparisons"
            )
        return result
