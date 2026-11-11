"""Check build: executable checks and their negative controls."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from itertools import product
from pathlib import Path
from zipfile import ZipFile

from scripts.archive_package import archive_package
from scripts.package_contract import check_package
from scripts.stage_package import PLATFORMS, stage
from validation import ROOT
from validation import check_build as deps


def dependencies_output(*args):
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT)


def dependencies_version(value):
    return tuple(int(p) for p in value.split(".")) + (0,) * (3 - len(value.split(".")))


def dependencies_check_macos(path, arch=None, archive=False):
    if arch and arch not in dependencies_output("lipo", "-archs", str(path)).split():
        raise ValueError(f"{path}: missing {arch} architecture")
    loads = dependencies_output("otool", "-l", str(path))
    mins = re.findall(r"\bminos\s+([\d.]+)", loads)
    mins += re.findall(r"cmd LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version ([\d.]+)", loads)
    if not mins or any(dependencies_version(v) > dependencies_version("11.0") for v in mins):
        raise ValueError(
            f"{path}: deployment target must be macOS 11.0 or earlier; found {sorted(set(mins))}"
        )
    if not archive:
        dependencies = re.findall(
            r"^\s+(\S+) \(compatibility", dependencies_output("otool", "-L", str(path)), re.M
        )
        unexpected = [
            d for d in dependencies if not d.startswith(("/usr/lib/", "/System/Library/"))
        ]
        if unexpected:
            raise ValueError(f"{path}: non-system dynamic dependencies: {unexpected}")


# The Linux floor is manylinux2014: CentOS 7's glibc 2.17. The OpenMP runtime
# is LLVM's libomp linked statically, so no OpenMP symbol version may be needed.
dependencies_LINUX_FLOOR = {"GLIBC": "2.17"}
dependencies_LINUX_EXPORTS = {"stata_call", "pginit"}


def dependencies_check_linux(path):
    needed = set(
        re.findall(r"Shared library: \[(.*?)\]", dependencies_output("readelf", "-d", str(path)))
    )
    allowed = {
        "libc.so.6",
        "libm.so.6",
        "libpthread.so.0",
        "librt.so.1",
        "libdl.so.2",
        "libgcc_s.so.1",
        "ld-linux-x86-64.so.2",
    }
    if needed - allowed:
        raise ValueError(f"{path}: unsupported dependencies: {needed - allowed}")
    if "libm.so.6" not in needed:
        raise ValueError(f"{path}: libm.so.6 is not linked; math symbols would rely on the host")
    versions = dependencies_output("readelf", "-W", "--version-info", str(path))
    for prefix, floor in dependencies_LINUX_FLOOR.items():
        newer = sorted(
            {
                v
                for v in re.findall(rf"\b{prefix}_(\d[\d.]*)", versions)
                if dependencies_version(v.rstrip(".")) > dependencies_version(floor)
            }
        )
        if newer:
            raise ValueError(f"{path}: requires {prefix} {newer}, newer than {prefix}_{floor}")
    openmp = sorted(set(re.findall(r"\b((?:GOMP|OMP|GOMP_PLUGIN)_\d[\d.]*)", versions)))
    if openmp:
        raise ValueError(f"{path}: needs a system OpenMP runtime ({openmp}); link the static libomp")
    symbols = dependencies_output("readelf", "-W", "--dyn-syms", str(path))
    # Only the plugin entry points may be exported; a leaked runtime symbol
    # could bind to another OpenMP runtime loaded in the same Stata process.
    exported = {
        name.split("@")[0]
        for name in re.findall(
            r"^\s*\d+:\s+\S+\s+\S+\s+(?:FUNC|OBJECT|TLS)\s+(?:GLOBAL|WEAK)\s+\S+\s+(?!UND)\S+\s+(\S+)", symbols, re.M
        )
    }
    if exported - dependencies_LINUX_EXPORTS:
        raise ValueError(f"{path}: exports symbols besides the plugin entry points: {sorted(exported - dependencies_LINUX_EXPORTS)[:8]}")
    unbound = [
        name
        for name in re.findall(r"^\s*\d+:\s+\S+\s+\S+\s+\S+\s+GLOBAL\s+\S+\s+UND\s+(\S+)", symbols, re.M)
        if "@" not in name
    ]
    if unbound:
        raise ValueError(f"{path}: undefined symbols bound to no library version: {sorted(unbound)[:8]}")
    if "not found" in dependencies_output("ldd", str(path)):
        raise ValueError(f"{path}: unresolved dynamic dependency")


def dependencies_check_windows(path):
    dependencies = {
        d.lower()
        for d in re.findall(r"DLL Name:\s*(\S+)", dependencies_output("objdump", "-p", str(path)))
    }
    allowed = {
        "kernel32.dll",
        "msvcrt.dll",
        "ucrtbase.dll",
        "advapi32.dll",
        "user32.dll",
        "shell32.dll",
        "ole32.dll",
    }
    unexpected = {
        d for d in dependencies if d not in allowed and not d.startswith("api-ms-win-crt-")
    }
    if unexpected:
        raise ValueError(f"{path}: non-system DLL dependencies: {unexpected}")


def dependencies_main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("platform", choices=["macos", "linux", "windows"])
    parser.add_argument("path", type=Path)
    parser.add_argument("--arch", choices=["arm64", "x86_64"])
    parser.add_argument("--archive", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.platform == "macos":
            dependencies_check_macos(args.path, args.arch, args.archive)
        elif args.platform == "linux":
            dependencies_check_linux(args.path)
        else:
            dependencies_check_windows(args.path)
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))
    print(f"Dependency contract checked: {args.path}")


def build_must_reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("incompatible dependency metadata was accepted")


def build_compiler_probes():
    make = shutil.which("make")
    with tempfile.TemporaryDirectory(prefix="ctools-compiler-probes-") as tmp:
        directory = Path(tmp)
        # A controlled PATH, independent of compilers installed on the host.
        for command in ("uname", "find", "sort", "grep", "file", "lipo"):
            path = shutil.which(command)
            if path:
                (directory / command).symlink_to(path)
        env = dict(os.environ, PATH=str(directory))
        rule = "$(info WIN=$(CC_WIN))\n$(info LINUX=$(CC_LINUX))\n$(info FLAGS=$(CFLAGS_WIN) $(LINUX_BASE_FLAGS))\np2-probe: ; @:"
        probe = directory / "probe.mk"
        probe.write_text("include Makefile\n" + rule + "\n")
        command = [
            make,
            "--no-print-directory",
            "-s",
            "-f",
            str(probe),
            "p2-probe",
            "DETECTED_OS=macOS",
            "LIBOMP_PREFIX=/nonexistent",
            "LLVM_MINGW_PREFIX=/nonexistent",
        ]
        absent = subprocess.check_output(command, cwd=ROOT, env=env, text=True)
        assert "WIN=\n" in absent and "LINUX=\n" in absent, absent
        for name in ("x86_64-w64-mingw32-gcc", "x86_64-linux-gnu-gcc"):
            path = directory / name
            path.write_text("#!/bin/sh\nexit 0\n")
            path.chmod(0o755)
        present = subprocess.check_output(command, cwd=ROOT, env=env, text=True)
        assert "WIN=x86_64-w64-mingw32-gcc\n" in present, present
        assert "LINUX=x86_64-linux-gnu-gcc\n" in present, present
        assert "-march=x86-64" in present and "-march=haswell" not in present
    print("PASS compiler availability and x86 baseline")


def build_dependency_metadata():
    mac = {
        "lipo": "arm64\n",
        "otool-l": "cmd LC_BUILD_VERSION\n minos 11.0\n",
        "otool-L": "plugin:\n /usr/lib/libSystem.B.dylib (compatibility version 1.0.0)\n",
    }

    def mac_output(*args):
        return mac[args[0] if args[0] == "lipo" else args[0] + args[1]]

    deps.dependencies_output = mac_output
    deps.dependencies_check_macos(Path("plugin"), "arm64")
    mac["otool-l"] = "cmd LC_BUILD_VERSION\n minos 26.0\n"
    build_must_reject(
        lambda: deps.dependencies_check_macos(Path("runtime.a"), "arm64", archive=True)
    )
    mac["otool-l"] = "cmd LC_BUILD_VERSION\n minos 11.0\n"
    mac["otool-L"] += " /opt/homebrew/lib/libomp.dylib (compatibility version 5.0.0)\n"
    build_must_reject(lambda: deps.dependencies_check_macos(Path("plugin"), "arm64"))
    linux = {
        "-d": "Shared library: [libc.so.6]\nShared library: [libm.so.6]\nShared library: [libpthread.so.0]",
        "--version-info": "Name: GLIBC_2.2.5\nName: GLIBC_2.17",
        "--dyn-syms": (
            "  1: 0000000000000000 0 FUNC GLOBAL DEFAULT UND exp@GLIBC_2.2.5 (2)\n"
            "  2: 0000000000000000 0 NOTYPE WEAK DEFAULT UND __gmon_start__\n"
            "  3: 00000000000a1b20 120 FUNC GLOBAL DEFAULT 12 stata_call\n"
            "  4: 00000000000a1ba0 40 FUNC GLOBAL DEFAULT 12 pginit\n"
        ),
        "ldd": "libm.so.6 => /lib64/libm.so.6",
    }
    deps.dependencies_output = lambda *args: linux["ldd" if args[0] == "ldd" else args[-2]]
    deps.dependencies_check_linux(Path("plugin"))
    accepted = dict(linux)
    for key, value in (
        ("--version-info", "Name: GLIBC_2.18"),
        ("--version-info", "Name: GLIBC_2.34"),
        ("--version-info", "Name: GOMP_4.0"),
        ("--version-info", "Name: OMP_1.0"),
        ("-d", "Shared library: [libc.so.6]\nShared library: [libpthread.so.0]"),
        ("-d", accepted["-d"] + "\nShared library: [libgomp.so.1]"),
        ("-d", accepted["-d"] + "\nShared library: [libstdc++.so.6]"),
        ("--dyn-syms", accepted["--dyn-syms"] + "  5: 0000000000000000 0 NOTYPE GLOBAL DEFAULT UND exp\n"),
        ("--dyn-syms", accepted["--dyn-syms"] + "  5: 00000000000b0000 64 FUNC GLOBAL DEFAULT 12 omp_get_thread_num\n"),
        ("--dyn-syms", accepted["--dyn-syms"] + "  5: 00000000000c0000 8 OBJECT WEAK DEFAULT 23 __kmp_global\n"),
        ("ldd", "libm.so.6 => not found"),
    ):
        linux.update(accepted)
        linux[key] = accepted[key] + "\n" + value if key == "--version-info" else value
        build_must_reject(lambda: deps.dependencies_check_linux(Path("plugin")))
    linux.update(accepted)
    for dll in ("libwinpthread-1.dll", "libgomp-1.dll", "libomp.dll"):
        deps.dependencies_output = lambda *args: "DLL Name: KERNEL32.dll\nDLL Name: " + dll
        build_must_reject(lambda: deps.dependencies_check_windows(Path("plugin")))
    print("PASS dependency rejection fixtures")


# A compiler that fails on the first dependency catches masked loop failures
# even when stale objects and a previous plugin are present.
def build_fail_fast_builds():
    with tempfile.TemporaryDirectory(prefix="ctools-build-failure-") as tmp:
        folder = Path(tmp)
        compiler = folder / "compiler"
        compiler.write_text(
            """#!/bin/sh
case " $* " in
  *" -c "*) echo compile >> "$CTOOLS_COMPILE_LOG"; exit 31 ;;
  *) echo link >> "$CTOOLS_COMPILE_LOG"; exit 0 ;;
esac
"""
        )
        compiler.chmod(0o755)
        for target, plugin, cc, flags, omp in (
            ("macos-arm", "ctools_mac_arm.plugin", "CC_MAC", "CFLAGS_MAC_ARM", "MAC_ARM_HAS_OMP"),
            ("macos-intel", "ctools_mac_x86.plugin", "CC_MAC", "CFLAGS_MAC_X86", "MAC_X86_HAS_OMP"),
            ("windows", "ctools_windows.plugin", "CC_WIN", "CFLAGS_WIN", "WIN_HAS_OMP"),
            ("linux", "ctools_linux.plugin", "CC_LINUX", "CFLAGS_LINUX", "LINUX_HAS_OMP"),
        ):
            build = folder / target
            build.mkdir()
            (build / plugin).write_bytes(b"previous plugin")
            (build / "ld_stale_arm.o").write_bytes(b"stale")
            log = folder / (target + ".log")
            env = dict(os.environ, CTOOLS_COMPILE_LOG=str(log))
            result = subprocess.run(
                [
                    "make",
                    "--no-print-directory",
                    target,
                    f"BUILD_DIR={build}",
                    f"{cc}={compiler}",
                    f"{flags}=",
                    f"{omp}=no",
                    "SOURCES=",
                    "HEADERS=",
                    "LIBDEFLATE_SRCS=first.c second.c",
                ],
                cwd=ROOT,
                env=env,
                capture_output=True,
                text=True,
            )
            assert result.returncode != 0, (target, result.stdout, result.stderr)
            # Ignore compiler feature probes performed while reading Makefile.
            lines = log.read_text().splitlines()
            assert lines.count("compile") == 1, (target, lines)
            assert lines[-1] == "compile", (target, lines)
            assert (build / plugin).read_bytes() == b"previous plugin"
            assert not list(build.glob(".ctools-*"))
    print("PASS four platform builds abort on first compile failure and preserve old plugin")


def build_depfile_paths():
    """Drive-letter targets must not hide changed compiler prerequisites."""
    from scripts.build_plugin import dependencies, dirty

    with tempfile.TemporaryDirectory(prefix="ctools-depfile-") as tmp:
        folder = Path(tmp)
        depfile = folder / "unit.d"
        source = folder / "unit.c"
        fragment = folder / "commands.def"
        obj = folder / "unit.o"
        for path in (source, fragment, obj):
            path.write_text("fixture")
            os.utime(path, ns=(1_000_000_000, 1_000_000_000))
        # GCC/Clang emit a drive colon in the -MT target on Windows. Test it
        # on every host, together with continuations and escaped spaces.
        for target in ("unit.o", r"G:\build\unit.o"):
            paths = [source, fragment, folder / "header with spaces.h"]
            paths[-1].write_text("fixture")
            os.utime(paths[-1], ns=(1_000_000_000, 1_000_000_000))
            words = [p.as_posix().replace(" ", r"\ ") for p in paths]
            depfile.write_text(target + ": " + " \\\n ".join(words) + "\n")
            assert dependencies(depfile) == paths
            assert not dirty(obj, depfile, source)
            os.utime(fragment, ns=(2_000_000_000, 2_000_000_000))
            assert dirty(obj, depfile, source)
            os.utime(fragment, ns=(1_000_000_000, 1_000_000_000))
        fragment.unlink()
        assert dirty(obj, depfile, source)
        depfile.write_text("malformed depfile")
        assert dirty(obj, depfile, source)
    print("PASS depfile drive letters, spaces, continuations, changed/missing .def inputs")


def build_incremental_builds():
    """Exercise real compiler depfiles, inventory changes, and atomic publication."""
    compiler = shutil.which(os.environ.get("CC", "clang"))
    assert compiler, "a C compiler is required"
    with tempfile.TemporaryDirectory(prefix="ctools-incremental-") as tmp:
        folder = Path(tmp)
        header, fragment = folder / "value.h", folder / "value.inc"
        main, extra = folder / "main.c", folder / "extra.c"
        header.write_text("#define VALUE 0\n")
        fragment.write_text("return VALUE;\n")
        main.write_text('#include "value.h"\nint main(void) {\n#include "value.inc"\n}\n')
        extra.write_text("int extra(void) { return 1; }\n")
        output = folder / "program"
        sources = [main, extra]

        def invoke(flags="", ldflags="", success=True):
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts/build_plugin.py"),
                    "--platform=test",
                    "--cc=" + compiler,
                    "--cflags=" + flags,
                    "--ldflags=" + ldflags,
                    "--sources=" + " ".join(map(str, sources)),
                    "--output",
                    str(output),
                    "--objects",
                    str(folder / "objects"),
                ],
                cwd=folder,
                text=True,
                capture_output=True,
            )
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result.stdout

        def touch(path):
            path.touch()
            # Guarantee strictly newer timestamps even on coarse filesystems.
            newest = max(p.stat().st_mtime_ns for p in (folder / "objects").rglob("*.o"))
            os.utime(path, ns=(newest + 1, newest + 1))

        assert invoke().count("Compile ") == 2
        assert "Up to date:" in invoke()
        for dependency in (header, fragment):
            touch(dependency)
            text = invoke()
            assert text.count("Compile ") == 1 and "main.c" in text and "Link " in text
        sources.remove(extra)
        text = invoke()
        assert "Compile " not in text and "Link " in text
        assert "Up to date:" in invoke()
        output.write_bytes(b"external replacement")
        assert "Link " in invoke()
        assert output.read_bytes() != b"external replacement"
        assert invoke("-DCHANGED=1").count("Compile ") == 1
        valid = output.read_bytes()
        invoke("-DCHANGED=1", "-lctools_missing_test_library", success=False)
        assert output.read_bytes() == valid
        main.write_text("this cannot compile\n")
        invoke("-DCHANGED=1", success=False)
        assert output.read_bytes() == valid
        assert not list(folder.glob(".ctools-*"))
        assert not list((folder / "objects").rglob(".compile-*"))
    print(
        "PASS incremental headers/fragments, no-op, source removal, flags, output integrity, and failed publication"
    )


def archive_reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("invalid release archive was accepted")


def archive_main(argv=None):
    workflow = (ROOT / ".github/workflows/build.yml").read_text()
    package_block = workflow.split("\n  package:\n")[1]
    package_condition = re.search(r"    if: >-\n(.*?)\n    runs-on:", package_block, re.S).group(1)
    build_jobs = ("build-linux", "build-macos-arm", "build-macos-intel", "build-windows")
    for event in ("push", "release", "workflow_dispatch", "pull_request"):
        for signing, mac_signing in product(("success", "failure", "skipped", "cancelled"), repeat=2):
            for failed_build in (None, *build_jobs):
                expression = package_condition.replace("always()", "True").replace(
                    "!cancelled()", "True"
                )
                for job in build_jobs:
                    expression = expression.replace(
                        f"needs.{job}.result", repr("failure" if job == failed_build else "success")
                    )
                expression = expression.replace("needs.sign-windows.result", repr(signing))
                expression = expression.replace("needs.sign-macos.result", repr(mac_signing))
                expression = expression.replace("github.event_name", repr(event))
                expression = " ".join(expression.replace("&&", "and").replace("||", "or").split())
                allowed = failed_build is None and (
                    signing == "success" or (event == "pull_request" and signing == "skipped")
                ) and (
                    mac_signing == "success" or (event == "pull_request" and mac_signing == "skipped")
                )
                assert eval(expression, {"__builtins__": {}}, {}) == allowed
    block = workflow.split("\n  release:\n")[-1]
    condition = re.search(r"    if: >-\n(.*?)\n    needs:", block, re.S).group(1)
    for event, ref, publish in (
        ("push", "refs/heads/main", False),
        ("push", "refs/heads/dev", False),
        ("pull_request", "refs/heads/main", False),
        ("workflow_dispatch", "refs/heads/main", False),
        ("release", "refs/tags/v1.0.0", True),
    ):
        for status in ("success", "failure", "skipped", "cancelled"):
            expression = condition.replace("always()", "True").replace("!cancelled()", "True")
            for key, value in {
                "needs.package.result": status,
                "github.event_name": event,
                "github.ref": ref,
            }.items():
                expression = expression.replace(key, repr(value))
            expression = " ".join(expression.replace("&&", "and").replace("||", "or").split())
            assert eval(expression, {"__builtins__": {}}, {}) == (publish and status == "success")
    assert "needs: [package]" in block
    assert "tag_name:" not in block
    assert "Publish main build release" not in block
    with tempfile.TemporaryDirectory(prefix="ctools-zip-") as tmp:
        root = Path(tmp)
        binaries = root / "binaries"
        binaries.mkdir()
        # Packaging fixtures deliberately do not pretend to be executable plugins.
        for platform in PLATFORMS:
            (binaries / f"ctools_{platform}.plugin").write_bytes(platform.encode())
        package = stage(
            ROOT / "build", root / "package", list(PLATFORMS), "archive-test-revision", binaries
        )
        (package / "scratch.log").write_text("must not ship")
        archive = archive_package(package, root / "dist")
        version = json.loads((ROOT / "config/release.json").read_text())["version"]
        assert archive.name == f"ctools-{version}.zip"
        info = json.loads((package / "BUILD_INFO.json").read_text())
        assert info["version"] == version
        assert (
            info["release_date"] == json.loads((ROOT / "config/release.json").read_text())["date"]
        )
        with ZipFile(archive) as zipped:
            assert set(zipped.namelist()) == set(info["sha256"]) | {"BUILD_INFO.json"}
            assert all("/" not in name for name in zipped.namelist())
            assert zipped.testzip() is None
            zipped.extractall(root / "installed")
        assert not check_package(root / "installed", integrity=True)
        assert (
            json.loads((root / "installed/BUILD_INFO.json").read_text())["revision"]
            == "archive-test-revision"
        )
        assert {p.name for p in (root / "installed").glob("*.plugin")} == {
            f"ctools_{platform}.plugin" for platform in PLATFORMS
        }

        single = stage(ROOT / "build", root / "single", ["linux"], "test", binaries)
        archive_reject(lambda: archive_package(single, root / "bad"))
        wrong_release = root / "release.json"
        wrong_release.write_text(json.dumps({"version": "999.0.0"}))
        archive_reject(lambda: archive_package(package, root / "bad", wrong_release))
        wrong_release.write_text(json.dumps({"version": "../escape"}))
        archive_reject(lambda: archive_package(package, root / "bad", wrong_release))
        original_archive = archive.read_bytes()
        (package / "ctools_windows.plugin").write_bytes(b"corrupt")
        archive_reject(lambda: archive_package(package, root / "dist"))
        assert archive.read_bytes() == original_archive
        (package / "ctools_linux.plugin").unlink()
        archive_reject(lambda: archive_package(package, root / "bad"))
        assert not (root / "bad").exists()
    print("PASS flat all-platform ZIP, extracted integrity, version, scope, and corruption checks")
