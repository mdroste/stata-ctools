"""Exercise compiler probes and reject incompatible release dependency metadata."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import sys

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "dependencies", ROOT / "validation/check_dependencies.py"
)
deps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(deps)


def must_reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("incompatible dependency metadata was accepted")


def compiler_probes():
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


def dependency_metadata():
    mac = {
        "lipo": "arm64\n",
        "otool-l": "cmd LC_BUILD_VERSION\n minos 11.0\n",
        "otool-L": "plugin:\n /usr/lib/libSystem.B.dylib (compatibility version 1.0.0)\n",
    }

    def mac_output(*args):
        return mac[args[0] if args[0] == "lipo" else args[0] + args[1]]

    deps.output = mac_output
    deps.check_macos(Path("plugin"), "arm64")
    mac["otool-l"] = "cmd LC_BUILD_VERSION\n minos 26.0\n"
    must_reject(lambda: deps.check_macos(Path("runtime.a"), "arm64", archive=True))
    mac["otool-l"] = "cmd LC_BUILD_VERSION\n minos 11.0\n"
    mac["otool-L"] += " /opt/homebrew/lib/libomp.dylib (compatibility version 5.0.0)\n"
    must_reject(lambda: deps.check_macos(Path("plugin"), "arm64"))
    linux = {
        "-d": "Shared library: [libc.so.6]\nShared library: [libgomp.so.1]",
        "--version-info": "GLIBC_2.35",
        "ldd": "libgomp.so.1 => /usr/lib/libgomp.so.1",
    }
    deps.output = lambda *args: linux["ldd" if args[0] == "ldd" else args[1]]
    deps.check_linux(Path("plugin"))
    linux["--version-info"] = "GLIBC_2.38"
    must_reject(lambda: deps.check_linux(Path("plugin")))
    linux["--version-info"] = "GLIBC_2.35"
    linux["ldd"] = "libgomp.so.1 => not found"
    must_reject(lambda: deps.check_linux(Path("plugin")))
    for dll in ("libwinpthread-1.dll", "libgomp-1.dll", "libomp.dll"):
        deps.output = lambda *args: "DLL Name: KERNEL32.dll\nDLL Name: " + dll
        must_reject(lambda: deps.check_windows(Path("plugin")))
    print("PASS dependency rejection fixtures")


# A compiler that fails on the first dependency catches masked loop failures
# even when stale objects and a previous plugin are present.
def fail_fast_builds():
    with tempfile.TemporaryDirectory(prefix="ctools-build-failure-") as tmp:
        folder = Path(tmp)
        compiler = folder / "compiler"
        compiler.write_text("""#!/bin/sh
case " $* " in
  *" -c "*) echo compile >> "$CTOOLS_COMPILE_LOG"; exit 31 ;;
  *) echo link >> "$CTOOLS_COMPILE_LOG"; exit 0 ;;
esac
""")
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


def incremental_builds():
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
        original = output.read_bytes()
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


if __name__ == "__main__":
    compiler_probes()
    dependency_metadata()
    fail_fast_builds()
    incremental_builds()
