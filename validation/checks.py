"""Checks: executable checks and their negative controls."""

import fnmatch
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from scripts import stage_package as packaging
from scripts.package_contract import check_package
from scripts.stage_package import stage
from validation import ROOT, stata
from validation.native import unit
from validation.stata import stata_command


def release_workflow_contract():
    source = (ROOT / ".github/workflows/build.yml").read_text()
    for job in ("commit-plugins", "release"):
        block = source.split("\n  " + job + ":\n")[-1]
        condition = re.search(r"    if: (.*?)\n    needs:", block, re.S).group(1)
        condition = condition.removeprefix(">-").strip()
        assert "needs.package.result" in condition
        for event in ("push", "pull_request", "workflow_dispatch", "release"):
            for ref in ("refs/heads/main", "refs/heads/dev", "refs/tags/v1.0.0"):
                for requested in (False, True):
                    for status in ("success", "failure", "skipped", "cancelled"):
                        for cancelled in (False, True):
                            expression = condition.replace("always()", "True")
                            expression = expression.replace("!cancelled()", repr(not cancelled))
                            for key, value in {
                                "needs.package.result": status,
                                "github.event_name": event,
                                "github.ref": ref,
                                "inputs.publish_plugins": requested,
                            }.items():
                                expression = expression.replace(key, repr(value))
                            expression = " ".join(
                                expression.replace("&&", "and").replace("||", "or").split()
                            )
                            publish = (
                                event == "release"
                                if job == "release"
                                else event == "workflow_dispatch"
                                and requested
                                and ref == "refs/heads/main"
                            )
                            allowed = publish and status == "success" and not cancelled
                            assert eval(expression, {"__builtins__": {}}, {}) == allowed, (
                                job,
                                event,
                                ref,
                                requested,
                                status,
                                cancelled,
                            )
    dispatch = source.split("  workflow_dispatch:\n")[1].split("  release:\n")[0]
    assert "publish_plugins:" in dispatch
    assert "type: boolean" in dispatch and "default: false" in dispatch
    paths = re.findall(r"    paths:\n((?:      - .*\n)+)", source)
    assert len(paths) == 2
    for block in paths:
        patterns = re.findall(r"'([^']+)'", block)
        for name in ("build/new.ado", "build/cimport.sthlp", "build/ctools.pkg"):
            assert any(fnmatch.fnmatchcase(name, pattern) for pattern in patterns)
    assert "stata-validation:" not in source
    assert "CTOOLS_STATA_RUNNERS" not in source
    for platform in ("linux", "mac_arm", "mac_x86", "windows"):
        assert f"name: ctools_{platform}" in source
    assert "needs.sign-windows.result == 'success'" in source
    with patch.dict(os.environ, {"CTOOLS_STATA_EXECUTABLE": "/test path/Stata"}):
        assert stata_command("driver.do")[0] == "/test path/Stata"
        assert stata_command("driver.do")[-1] == str(Path("driver.do").resolve())
    with patch(
        "validation.stata.os",
        SimpleNamespace(
            name="nt", environ={"CTOOLS_STATA_EXECUTABLE": "C:/Program Files/Stata/StataMP.exe"}
        ),
    ):
        command = stata_command("driver.do")
        assert command[0] == "C:/Program Files/Stata/StataMP.exe"
        assert command[1:3] == ["/e", "do"]
    print("PASS publication requires explicit action and successful checks; pushes only build")


def release_packages():
    with tempfile.TemporaryDirectory(prefix="ctools-release-integrity-") as tmp:
        directory = stage(ROOT / "build", Path(tmp) / "package", ["mac_arm"], "test-only")
        assert not check_package(directory, integrity=True)
        info_path = directory / "BUILD_INFO.json"
        original = info_path.read_bytes()
        for name in ("ctools_mac_arm.plugin", "cppmlhdfe_p.ado"):
            payload_path = directory / name
            payload = payload_path.read_bytes()
            payload_path.write_bytes(b"corrupt or truncated payload")
            assert any("checksum mismatch" in e for e in check_package(directory, integrity=True))
            payload_path.write_bytes(payload)
        for mutation in ("missing", "scope", "inventory", "digest", "manifest", "invalid"):
            info = json.loads(original)
            if mutation == "missing":
                info_path.unlink()
            elif mutation == "invalid":
                info_path.write_text("null")
            else:
                if mutation == "scope":
                    info["platforms"] = ["linux"]
                if mutation == "inventory":
                    info["sha256"].pop("cppmlhdfe_p.ado")
                if mutation == "digest":
                    info["sha256"]["ctools.pkg"] = "0" * 64
                if mutation == "manifest":
                    (directory / "ctools.pkg").write_text(
                        (directory / "ctools.pkg").read_text() + "d changed\n"
                    )
                info_path.write_text(json.dumps(info))
            assert check_package(directory, integrity=True), mutation
            if mutation == "manifest":
                p = directory / "ctools.pkg"
                p.write_text(p.read_text().removesuffix("d changed\n"))
            info_path.write_bytes(original)
        assert not check_package(directory, integrity=True)
    print("PASS checksum corruption, inventory, scope, manifest, and missing build-info checks")


def release_cleanup():
    with tempfile.TemporaryDirectory(prefix="ctools-clean-") as tmp:
        directory = Path(tmp)
        build, source = directory / "build", directory / "src"
        build.mkdir()
        source.mkdir()
        (build / "obj").mkdir()
        artifacts = [
            build / "ctools.plugin",
            build / "ctools_windows.plugin",
            build / "x.o",
            source / "x.o",
            build / "obj/x.o",
        ]
        for artifact in artifacts:
            artifact.touch()
        result = subprocess.run(
            [
                "make",
                "-s",
                "-f",
                str(ROOT / "Makefile"),
                "clean",
                "DETECTED_OS=Windows",
                f"BUILD_DIR={build}",
                f"SRC_DIR={source}",
                "LIBOMP_PREFIX=/nonexistent",
            ],
            cwd=ROOT,
            check=True,
        )
        assert all(not p.exists() for p in artifacts)
        # A matching directory cannot be unlinked: genuine cleanup errors propagate.
        (build / "bad.o").mkdir()
        result = subprocess.run(
            ["python3", "scripts/clean_build.py", str(build), str(source)],
            cwd=ROOT,
            capture_output=True,
        )
        assert result.returncode != 0
    print("PASS MSYS/POSIX recipe cleanup, generic plugin, and error propagation")


def package_must_reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("invalid release input was accepted")


def package_packages():
    with tempfile.TemporaryDirectory(prefix="ctools-package-test-") as tmp:
        directory = Path(tmp)
        source = directory / "source"
        source.mkdir()
        # Real source manifests and helpers; dummy binaries test staging only.
        for entry in (ROOT / "build").iterdir():
            if entry.is_file() and (
                entry.suffix in (".ado", ".sthlp", ".pkg")
                or entry.name in ("stata.toc", "LICENSE", "THIRD_PARTY_NOTICES")
            ):
                shutil.copy2(entry, source / entry.name)
        # Every other non-plugin file the manifest lists (e.g. codec licences)
        for line in (ROOT / "build" / "ctools.pkg").read_text().splitlines():
            name = line[2:]
            if line.startswith("f ") and not name.endswith(".plugin"):
                shutil.copy2(ROOT / "build" / name, source / name)
        for target in packaging.PLATFORMS:
            (source / f"ctools_{target}.plugin").write_bytes(b"test fixture, not a plugin")
        for targets in (["mac_arm"], list(packaging.PLATFORMS)):
            staged = packaging.stage(
                source,
                directory / ("stage-" + targets[0] + str(len(targets))),
                targets,
                "test-only",
            )
            assert not packaging.check_package(staged)
            assert (staged / "cppmlhdfe_p.ado").is_file()
            assert {p.name for p in staged.glob("*.plugin")} == {
                f"ctools_{p}.plugin" for p in targets
            }
            (staged / "cppmlhdfe_p.ado").unlink()
            assert packaging.check_package(staged)
        (source / "ctools_linux.plugin").unlink()
        package_must_reject(
            lambda: packaging.stage(
                source, directory / "bad", list(packaging.PLATFORMS), "test-only"
            )
        )
        assert not (directory / "bad").exists()
        (source / "cppmlhdfe_p.ado").unlink()
        package_must_reject(
            lambda: packaging.stage(source, directory / "bad-helper", ["mac_arm"], "test-only")
        )
    print("PASS single/all-platform manifests, helper completeness, and atomic staging")


def package_inventory():
    public = {p.stem for p in (ROOT / "build").glob("*.sthlp")} - {"ctools"}
    registered = {p.stem for p in (ROOT / "validation/stata").glob("c*.do")} - {"contracts"}
    available, _ = stata.cases()
    assert public <= available.keys()
    assert public == registered, (public - registered, registered - public)
    manifest = (ROOT / "build/ctools.pkg").read_text().splitlines()
    for name in public:
        assert (ROOT / "src" / name).is_dir(), name
        for suffix in (".ado", ".sthlp"):
            assert (ROOT / "build" / (name + suffix)).is_file(), name
            assert f"f {name}{suffix}" in manifest, name
    print(f"PASS {len(public)} public-command help/source/package/test inventory")


docs_BENCHMARKS = ("docs/benchmarks/hdfe_datasets_", "docs/benchmarks/readme_profile/")
docs_DATED = re.compile(r"(?<!\d)(19|20)\d{2}-?[01]\d-?[0-3]\d(?!\d)")
docs_ARTIFACTS = {".do", ".py", ".plugin", ".log", ".patch", ".dta", ".o"}


def docs_main(argv=None):
    files = subprocess.run(
        (
            ["/bin/zsh", "-lic", "git ls-files --cached --others --exclude-standard docs"]
            if sys.platform == "darwin"
            else ["git", "ls-files", "--cached", "--others", "--exclude-standard", "docs"]
        ),
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    files = [name for name in files if (ROOT / name).exists()]
    errors = []
    for name in files:
        if docs_DATED.search(name):
            errors.append(f"{name}: dated path; move snapshots to reports/")
        elif Path(name).suffix in docs_ARTIFACTS:
            errors.append(f"{name}: run artifact; keep drivers in validation/ and output in temp/")
        elif name.startswith("docs/benchmarks/") and not name.startswith(docs_BENCHMARKS):
            errors.append(f"{name}: not regenerated evidence; see docs/development/README.md")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"docs tree checked ({len(files)} files)")


def windows_main(argv=None):
    if os.name != "nt":
        raise SystemExit("This regression requires a real Windows worker")
    with tempfile.TemporaryDirectory(prefix="ctools-win-paths-") as tmp:
        binary = Path(tmp) / "paths.exe"
        subprocess.run(
            [
                os.environ.get("CC", "gcc"),
                "-std=c11",
                "-O1",
                "-DSYSTEM=STWIN32",
                "-DSD_FASTMODE",
                "-I",
                str(ROOT / "src"),
                str(unit(Path(tmp), "oct21_windows")),
                str(ROOT / "src/formats/ctools_output.c"),
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary), tmp], check=True)
        build, source = Path(tmp) / "build", Path(tmp) / "src"
        build.mkdir()
        source.mkdir()
        for name in ("ctools.plugin", "ctools_windows.plugin", "stale.o"):
            (build / name).touch()
        (source / "stale.o").touch()
        subprocess.run(
            [
                "make",
                "clean",
                "DETECTED_OS=Windows",
                f"BUILD_DIR={build.as_posix()}",
                f"SRC_DIR={source.as_posix()}",
            ],
            cwd=ROOT,
            check=True,
        )
        assert not list(build.iterdir()) and not list(source.iterdir())


def run(out):
    """Run independently so deliberate mocks cannot contaminate later checks."""
    import subprocess
    import sys

    targets = [
        ("build.compiler_probes", "validation.check_build", "build_compiler_probes"),
        ("build.dependencies", "validation.check_build", "build_dependency_metadata"),
        ("build.fail_fast", "validation.check_build", "build_fail_fast_builds"),
        ("build.incremental", "validation.check_build", "build_incremental_builds"),
        ("release.archive", "validation.check_build", "archive_main"),
        ("release.workflow", "validation.checks", "release_workflow_contract"),
        ("release.integrity", "validation.checks", "release_packages"),
        ("release.cleanup", "validation.checks", "release_cleanup"),
        ("package.staging", "validation.checks", "package_packages"),
        ("package.inventory", "validation.checks", "package_inventory"),
        ("runner", "validation.check_runner", "test_runner"),
        ("options", "validation.check_runner", "options_main"),
        ("benchmark.evidence", "validation.check_runner", "test_performance"),
        ("docs", "validation.checks", "docs_main"),
    ]
    failures = []
    for name, module, function in targets:
        log = out / (name + ".log")
        with log.open("w") as handle:
            result = subprocess.run(
                [sys.executable, "-c", f"from {module} import {function}; {function}()"],
                cwd=ROOT,
                stdout=handle,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            failures.append(name)
            print(f"FAIL {name}: {log}\n" + log.read_text(errors="replace")[-6000:], flush=True)
        else:
            print("PASS " + name, flush=True)
    return failures


def main(argv):
    import argparse

    from validation import check_build

    if argv[0] == "dependencies":
        return check_build.dependencies_main(argv[1:])
    if argv[0] == "windows":
        return windows_main(argv[1:])
    if argv[0] == "docs":
        return docs_main(argv[1:])
    if argv[0] == "options":
        from validation import options

        parser = argparse.ArgumentParser(description="Check the public option inventory")
        parser.add_argument("--write", action="store_true")
        args = parser.parse_args(argv[1:])
        rows = options.records()
        gaps = [f"{r['command']}.{r['option']}" for r in rows if r["status"] == "gap"]
        if gaps:
            raise SystemExit("Options without coverage: " + ", ".join(gaps))
        document = ROOT / "docs/OPTION_COVERAGE.md"
        expected = options.report(rows)
        if args.write:
            document.write_text(expected)
        elif document.read_text() != expected:
            raise SystemExit("Option inventory is stale; run --checks options --write")
        print(f"{len(rows)} public options checked; no coverage gaps")
        return
    if argv[0] == "package":
        from scripts.package_contract import check_package

        parser = argparse.ArgumentParser()
        parser.add_argument("directory", type=Path)
        parser.add_argument("--metadata-only", action="store_true")
        parser.add_argument("--integrity", action="store_true")
        args = parser.parse_args(argv[1:])
        errors = check_package(args.directory, args.metadata_only, args.integrity)
        if errors:
            raise SystemExit("\n".join(errors))
        print("Package checked")
        return
    raise SystemExit("Unknown repository check: " + argv[0])
