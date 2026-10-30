"""Check runner: executable checks and their negative controls."""

import contextlib
import io
import json
import platform
import re
import shlex
import shutil
import tempfile
import unittest
from functools import partial
from pathlib import Path
from unittest.mock import patch

from validation import bench as benchmark
from validation import options as coverage
from validation import stata as launcher
from validation.bench import audit_summarize as summarize
from validation.options import CASES, source, validate_cases


class performance_AuditEvidence(unittest.TestCase):
    def fixture(self, directory, label="b1", *, seconds="0.2", signature="sig"):
        out = Path(directory)
        (out / f"{label}.manifest.json").write_text(
            json.dumps(dict(label=label, sizes=[2000], threads=1, reps=1, cases=["split"]))
        )
        log = (
            f"AUDIT,{label},2000,split,1,0,0.3,0,2000,43\n"
            f"AUDIT,{label},2000,split,1,1,{seconds},0,2000,43\n"
            f"AUDIT_CHECK,{label},2000,split,1,{signature},.\n"
            f"AUDIT_COMPLETE label={label}\nCTOOLS_LAUNCH=abcd RC=0\n"
        )
        (out / f"{label}.log").write_text(log)
        return out, log

    def reject(self, transform):
        with tempfile.TemporaryDirectory() as directory:
            out, log = self.fixture(directory)
            (out / "b1.log").write_text(transform(log))
            with self.assertRaises(ValueError):
                summarize(out)
            self.assertFalse((out / "summary.csv").exists())

    def test_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            out, _ = self.fixture(directory)
            self.fixture(directory, "c1")
            with contextlib.redirect_stdout(io.StringIO()):
                summarize(out)
            self.assertIn("b,2000,split,1,1,0.2", (out / "summary.csv").read_text())
            self.assertEqual(len((out / "raw.csv").read_text().splitlines()), 5)

    def test_missing_trial(self):
        self.reject(lambda log: log.replace("AUDIT,b1,2000,split,1,1,0.2,0,2000,43\n", ""))

    def test_missing_check(self):
        self.reject(lambda log: log.replace("AUDIT_CHECK,b1,2000,split,1,sig,.\n", ""))

    def test_duplicate_trial(self):
        self.reject(lambda log: log + "AUDIT,b1,2000,split,1,1,0.2,0,2000,43\n")

    def test_duplicate_check(self):
        self.reject(lambda log: log + "AUDIT_CHECK,b1,2000,split,1,sig,.\n")

    def test_missing_completion(self):
        self.reject(lambda log: log.replace("CTOOLS_LAUNCH=abcd RC=0", "CTOOLS_LAUNCH=abcd RC=1"))

    def test_nonfinite_or_nonpositive(self):
        for value in ["nan", "inf", "0", "-1"]:
            with self.subTest(value=value):
                self.reject(lambda log: log.replace(",1,0.2,", f",1,{value},"))

    def test_different_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            out, _ = self.fixture(directory)
            self.fixture(directory, "c1", signature="different")
            with self.assertRaisesRegex(ValueError, "output mismatch"):
                summarize(out)
            self.assertFalse((out / "summary.csv").exists())


def options_reject(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError("invalid coverage accepted")


def options_main(argv=None):
    rows = coverage.records()
    assert rows and not [r for r in rows if r["status"] == "gap"]
    # Transaction wrappers must not hide their implementation's public syntax.
    assert {"cencode", "cdecode", "cmerge", "cwinsor", "crangestat", "creghdfe"} <= {
        r["command"] for r in rows
    }
    assert any(r["command"] == "cencode" and r["option"] == "replace" for r in rows)
    assert (coverage.ROOT / "docs/OPTION_COVERAGE.md").read_text() == coverage.report(rows)
    cases = [(0, CASES[0]), (1, CASES[1])]
    complete = "\n".join(
        f"OPTION_CASE={i} VARIANT={v} FAILED=0 RC=0" for i, _ in cases for v in (1, 2, 3)
    )
    assert len(validate_cases(complete, cases)) == 6
    options_reject(lambda: validate_cases(complete.split("\n", 1)[1], cases))
    options_reject(lambda: validate_cases(complete + "\n" + complete, cases))
    options_reject(
        lambda: validate_cases(complete + "\nOPTION_CASE=99 VARIANT=1 FAILED=0 RC=0", cases)
    )
    options_reject(lambda: validate_cases(complete.replace("FAILED=0", "FAILED=2", 1), cases))
    failed = validate_cases(complete.replace("FAILED=0 RC=0", "FAILED=1 RC=198", 1), cases)
    assert failed[(0, 1)] == dict(failed=1, rc=198)
    driver = source(cases, Path("/tmp/option coverage.log"))
    assert driver.rstrip().endswith("exit, clear") and "OPTION_COVERAGE_COMPLETE" in driver
    assert coverage.exercised(
        dict(options=["generate"], candidate="cencode s, gen(out)"), "generate", "GENerate"
    )
    assert not coverage.exercised(
        dict(options=["generate"], candidate="cencode s, agent(out)"), "generate", "GENerate"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "build").mkdir()
        for ado in (coverage.ROOT / "build").glob("*.ado"):
            shutil.copy2(ado, root / "build" / ado.name)
        shutil.copytree(
            coverage.ROOT / "validation",
            root / "validation",
            ignore=shutil.ignore_patterns("*.log", "__pycache__"),
        )
        path = root / "build/csort.ado"
        path.write_text(
            path.read_text().replace("syntax ", "syntax ", 1)
            + "\nprogram define csort\nsyntax, NEWcoverageoption\nend\n"
        )
        with patch.object(coverage, "ROOT", root):
            gaps = [r for r in coverage.records() if r["status"] == "gap"]
            assert [(r["command"], r["option"]) for r in gaps] == [("csort", "newcoverageoption")]
            evidence = root / "validation/option_evidence.json"
            data = json.loads(evidence.read_text())
            data["ctools"]["version"]["evidence"][0]["anchor"] = "ABSENT_TEST_ASSERTION_1234"
            evidence.write_text(json.dumps(data))
            options_reject(coverage.records)
    print(
        "PASS option inventory mutations, stale evidence, incomplete/duplicate results, and failed-case retention"
    )


def reject(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError("invalid input accepted")


def check_io_cleanup(root):
    for fail in (False, True):
        try:
            with launcher.io_workspace() as workspace:
                assert (workspace / "validation/support/io_test_helpers.do").is_file()
                assert (workspace / "build").is_dir()
                (workspace / "temp/io_allformats/native/output.csv").write_text("fixture")
                if fail:
                    raise ValueError("injected test failure")
        except ValueError:
            assert fail
        assert not workspace.exists()

    build = root / "build"
    build.mkdir()
    plugin = (
        ("ctools_mac_arm.plugin" if platform.machine() == "arm64" else "ctools_mac_x86.plugin")
        if platform.system() == "Darwin"
        else "ctools_linux.plugin"
    )
    (build / plugin).write_bytes(b"fake frozen build")
    for mode in ("success", "launch_failure", "mismatch", "keep"):
        out = root / mode
        out.mkdir()
        (out / "unrelated.csv").write_text("preserve me")
        (out / "comparisons.txt").write_text("stale pass")
        workspaces = []

        def fake_launch(driver, log):
            command = next(
                line
                for line in driver.read_text().splitlines()
                if line.startswith("capture noisily do ")
            )
            data = Path(shlex.split(command)[-1])
            workspaces.append(data)
            label = driver.stem
            log.write_text(f"IO_BENCHMARK_RC=0\nPERF,{label},small,export_csv,1,1,0.1\n")
            (data / f"{label}_small_1.csv").write_text(
                "different" if mode == "mismatch" and label == "c1" else "same bytes"
            )
            if mode == "launch_failure":
                raise ValueError("injected Stata failure")

        argv = [
            "prepare_io_performance.py",
            "--baseline",
            str(build),
            "--candidate",
            str(build),
            "--output",
            str(out),
            "--cases",
            "small",
            "--reps",
            "1",
            "--run",
        ]
        if mode == "keep":
            argv.append("--keep-data")
        with patch("sys.argv", argv), patch.object(benchmark.subprocess, "run"), patch.object(
            benchmark, "launch_driver", fake_launch
        ):
            if mode in ("launch_failure", "mismatch"):
                reject(benchmark.io_main)
                assert not (out / "comparisons.txt").exists()
            else:
                benchmark.io_main()
                assert (out / "comparisons.txt").read_text().startswith("PASS 3 ")
                assert (out / "timings.csv").is_file() and (out / "summary.csv").is_file()
        assert (out / "unrelated.csv").read_text() == "preserve me"
        assert (out / "b1.log").is_file()
        assert workspaces
        if mode == "keep":
            assert (out / "c1_small_1.csv").is_file()
            benchmark.io_summarize(out)
        else:
            assert all(not path.exists() for path in workspaces)


def check_dataset_cleanup(root):
    build = root / "build"
    for command in ("cipolate", "csplit", "crangejoin"):
        (build / (command + ".ado")).write_text("fake frozen ado")
    for name in ("command", "encode", "group"):
        runner = benchmark
        run = (
            benchmark.commands_main
            if name == "command"
            else partial(benchmark.abba_main, kind=name)
        )
        summary_name = "command_summarize" if name == "command" else "abba_summarize"
        for mode in ("success", "launch_failure", "summary_failure", "keep", "prepare"):
            if mode == "prepare" and name != "command":
                continue
            out = root / f"{name}_{mode}"
            out.mkdir()
            (out / "unrelated.dta").write_bytes(b"preserve me")
            workspaces = []

            def fake_launch(driver, log):
                command = next(
                    line
                    for line in driver.read_text().splitlines()
                    if line.startswith("capture noisily do ")
                )
                data = Path(shlex.split(command)[-1])
                workspaces.append(data)
                (data / "expected_fixture.dta").write_bytes(b"generated comparison data")
                label = driver.stem
                if name == "command":
                    lines = [f"PERFORMANCE_COMPLETE label={label}", "PERFORMANCE_DRIVER_RC=0"]
                    lines += [
                        f"PERF,{label},ip_sorted,{threads},{rep},0.1,1000"
                        for threads in (1, 4, 12)
                        for rep in range(1, 7)
                    ]
                else:
                    lines = [
                        f"{name.upper()}_BENCHMARK_RC=0",
                        f"PERF,{label},fixture,ctools,1000,1,1,0.1",
                        f"PERF,{label},fixture,reference,1000,0,1,0.2",
                    ]
                log.write_text("\n" + "\n".join(lines) + "\n")
                if mode == "launch_failure":
                    raise ValueError("injected Stata failure")

            argv = [
                runner.__file__,
                "--baseline",
                str(build),
                "--candidate",
                str(build),
                "--output",
                str(out),
            ]
            if name == "command":
                argv += ["--cases", "ip_sorted"]
                if mode != "prepare":
                    argv.append("--run")
            if mode == "keep":
                argv.append("--keep-data")
            with patch("sys.argv", argv), patch.object(runner.subprocess, "run"), patch.object(
                runner, "launch_driver", fake_launch
            ):
                if mode == "summary_failure":
                    with patch.object(runner, summary_name, side_effect=ValueError("bad summary")):
                        reject(run)
                elif mode == "launch_failure":
                    reject(run)
                else:
                    run()
            assert (out / "unrelated.dta").read_bytes() == b"preserve me"
            if mode == "prepare":
                assert not workspaces and (out / "b1.do").is_file()
                assert not list(out.glob("data-*"))
                continue
            assert workspaces and (out / "b1.log").is_file()
            if mode == "keep":
                assert (out / "expected_fixture.dta").is_file()
            else:
                assert all(not path.exists() for path in workspaces)
            if mode in ("success", "keep"):
                assert (out / "timings.csv").is_file()
                if name != "command":
                    assert (out / "summary.csv").is_file()


def test_performance():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(performance_AuditEvidence)
    result = unittest.TextTestRunner().run(suite)
    if not result.wasSuccessful():
        raise AssertionError("Benchmark evidence checks failed")


def test_runner():
    available, _ = launcher.cases()
    assert "memory_safety" in available
    for marker in (
        "",
        "\nCTOOLS_LAUNCH=stale RC=0\n",
        "\nCTOOLS_LAUNCH=now RC=198\n",
        "\nCTOOLS_LAUNCH=now RC=0\nCTOOLS_LAUNCH=now RC=0\n",
    ):
        reject(lambda: launcher.validate_completion(marker, "now"))
    launcher.validate_completion("\nCTOOLS_LAUNCH=now RC=0\n", "now")
    expected = {"csort": True, "cmerge": True}
    log = "\n".join(
        f"TEST_RESULT {n} PASSED=100 FAILED=0 SKIPPED=0 RC=0 COMPLETE={n}" for n in expected
    )
    assert len(launcher.validate_results(log, expected)) == 2
    for bad in (
        log.replace("COMPLETE=cmerge", "COMPLETE="),
        log.replace("RC=0", "RC=679"),
        log.replace("FAILED=0", "FAILED=1"),
        log.replace("SKIPPED=0", "SKIPPED=1"),
        log.replace("PASSED=100", "PASSED=0"),
        log.replace("TEST_RESULT", "MISSING_RESULT"),
        log.split("\n", 1)[0],
        log + "\n" + log,
    ):
        reject(lambda: launcher.validate_results(bad, expected))
    with tempfile.TemporaryDirectory(prefix="ctools runner spaces ") as tmp:
        root = Path(tmp)
        folder = root / "validation/stata"
        folder.mkdir(parents=True)
        (root / "build").mkdir()
        (root / "build/csort.sthlp").touch()
        with patch.object(launcher, "ROOT", root):
            reject(launcher.cases)  # A public command without any tests.
            path = folder / "csort.do"
            path.write_text('if "`__case\'" == "csort" {\n}\n')
            reject(launcher.cases)  # A selected case that can never finish.
            path.write_text('if "`__case\'" == "csort" {\nglobal CTOOLS_CASE_COMPLETE "csort"\n}\n')
            assert set(launcher.cases()[0]) == {"csort"}
            (folder / "duplicate.do").write_text(path.read_text())
            reject(launcher.cases)
        driver = root / "driver with spaces.do"
        logpath = root / "out.log"
        driver.write_text(
            f'log using "{logpath}", text replace\nlocal rc = 0\nlog close\nexit, clear\n'
        )
        launched = []

        def fake_run(argv, **kwargs):
            assert argv[:2] == ["/bin/zsh", "-lic"]
            command = shlex.split(argv[2])
            assert command[:4] == ["stata", "-q", "-b", "do"]
            path = Path(command[4])
            launched.append(path)
            text = path.read_text()
            assert text.rstrip().endswith("exit, clear")
            assert Path(kwargs["cwd"]).resolve() == path.parent
            assert text.startswith("cd " + launcher.quote(root.resolve()) + "\n")
            nonce = re.search("CTOOLS_LAUNCH=([a-f0-9]+)", text)[1]
            logpath.write_text(f"\nCTOOLS_LAUNCH={nonce} RC=0\n")

        with patch.object(launcher.subprocess, "run", fake_run):
            launcher.launch_driver(driver, logpath, cwd=root)
        assert launched and not launched[0].exists()
        driver.write_text("local rc = 0\n")
        reject(lambda: launcher.launch_driver(driver, logpath))
        # The cleanup checks create their own build fixture.
        (root / "build/csort.sthlp").unlink()
        (root / "build").rmdir()
        check_io_cleanup(root)
        check_dataset_cleanup(root)
    print(
        "PASS discovery, stale/missing/duplicate/error/empty/skipped results, path quoting, driver cleanup"
    )
    print(
        "PASS I/O and benchmark cleanup on success/failure, retained logs, keep-data, and preparation"
    )
