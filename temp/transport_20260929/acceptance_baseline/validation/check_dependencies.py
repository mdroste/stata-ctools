"""Fail release builds whose runtime dependencies exceed the platform contract."""
import argparse
from pathlib import Path
import re
import subprocess


def output(*args):
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT)


def version(value):
    return tuple(int(p) for p in value.split(".")) + (0,) * (3-len(value.split(".")))


def check_macos(path, arch=None, archive=False):
    if arch and arch not in output("lipo", "-archs", str(path)).split():
        raise ValueError(f"{path}: missing {arch} architecture")
    loads = output("otool", "-l", str(path))
    mins = re.findall(r"\bminos\s+([\d.]+)", loads)
    mins += re.findall(r"cmd LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version ([\d.]+)", loads)
    if not mins or any(version(v) > version("11.0") for v in mins):
        raise ValueError(f"{path}: deployment target must be macOS 11.0 or earlier; found {sorted(set(mins))}")
    if not archive:
        dependencies = re.findall(r"^\s+(\S+) \(compatibility", output("otool", "-L", str(path)), re.M)
        unexpected = [d for d in dependencies if not d.startswith(("/usr/lib/", "/System/Library/"))]
        if unexpected:
            raise ValueError(f"{path}: non-system dynamic dependencies: {unexpected}")


def check_linux(path):
    needed = set(re.findall(r"Shared library: \[(.*?)\]", output("readelf", "-d", str(path))))
    allowed = {"libc.so.6", "libm.so.6", "libgomp.so.1", "libpthread.so.0", "librt.so.1", "libdl.so.2", "libgcc_s.so.1", "ld-linux-x86-64.so.2"}
    if needed - allowed:
        raise ValueError(f"{path}: unsupported dependencies: {needed - allowed}")
    versions = re.findall(r"GLIBC_([\d.]+)", output("readelf", "--version-info", str(path)))
    if any(version(v) > version("2.35") for v in versions):
        raise ValueError(f"{path}: requires glibc newer than 2.35")
    if "not found" in output("ldd", str(path)):
        raise ValueError(f"{path}: unresolved dynamic dependency")


def check_windows(path):
    dependencies = {d.lower() for d in re.findall(r"DLL Name:\s*(\S+)", output("objdump", "-p", str(path)))}
    allowed = {"kernel32.dll", "msvcrt.dll", "ucrtbase.dll", "advapi32.dll", "user32.dll", "shell32.dll", "ole32.dll"}
    unexpected = {d for d in dependencies if d not in allowed and not d.startswith("api-ms-win-crt-")}
    if unexpected:
        raise ValueError(f"{path}: non-system DLL dependencies: {unexpected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("platform", choices=["macos", "linux", "windows"])
    parser.add_argument("path", type=Path)
    parser.add_argument("--arch", choices=["arm64", "x86_64"])
    parser.add_argument("--archive", action="store_true")
    args = parser.parse_args()
    try:
        if args.platform == "macos":
            check_macos(args.path, args.arch, args.archive)
        elif args.platform == "linux":
            check_linux(args.path)
        else:
            check_windows(args.path)
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))
    print(f"Dependency contract checked: {args.path}")


if __name__ == "__main__":
    main()
