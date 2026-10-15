"""Incremental compiler driver shared by the four Makefile platform targets.

Objects and dependency files are private to a platform and compiler/flag digest.
Only a successful link replaces the published plugin. Source inventory changes
invalidate linking, including deletions; compiler depfiles track handwritten .inc
files as well as headers. No generated source or network access is required.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile


def dependencies(path):
    text = path.read_text().replace("\\\n", " ")
    return [Path(p) for p in shlex.split(text.split(":", 1)[1])]


def dirty(obj, depfile, source):
    if not obj.exists() or not depfile.exists():
        return True
    try:
        return any(
            not p.exists() or p.stat().st_mtime_ns > obj.stat().st_mtime_ns
            for p in [source, *dependencies(depfile)]
        )
    except (OSError, ValueError):
        return True


def package_flags(source):
    """Only a dependency's own compilation sees its compatibility headers."""
    name = source.as_posix()
    flags = []
    if "/vendor/" in name:
        flags.append("-w")
    if "/vendor/libxls/" in name:
        flags += ["-Isrc/vendor/libxls", "-Isrc/vendor/libxls/include"]
    if "/vendor/readstat/" in name:
        flags += ["-DHAVE_ZLIB=1", "-Isrc/vendor/zlib"]
    return flags


def build(args):
    compiler = shlex.split(args.cc)
    if not compiler:
        raise ValueError("No compiler available for " + args.platform)
    executable = shutil.which(compiler[0]) or compiler[0]
    compiler_state = os.stat(executable)
    flags = shlex.split(args.cflags)
    config = [
        compiler,
        executable,
        compiler_state.st_mtime_ns,
        compiler_state.st_size,
        flags,
        Path(__file__).read_text(),
    ]
    digest = hashlib.sha256(json.dumps(config).encode()).hexdigest()[:16]
    directory = args.objects / args.platform / digest
    directory.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    objects = []
    for source in map(Path, shlex.split(args.sources)):
        # Preserve paths in object names; equal basenames never collide.
        key = hashlib.sha256(str(source).encode()).hexdigest()[:12]
        obj = directory / (source.stem + "-" + key + ".o")
        depfile = obj.with_suffix(".d")
        if dirty(obj, depfile, source):
            print("Compile " + str(source), flush=True)
            with tempfile.TemporaryDirectory(prefix=".compile-", dir=directory) as tmp:
                temp_obj, temp_dep = Path(tmp) / "object.o", Path(tmp) / "object.d"
                subprocess.run(
                    [
                        *compiler,
                        *flags,
                        *package_flags(source),
                        "-MMD",
                        "-MF",
                        str(temp_dep),
                        "-MT",
                        str(obj),
                        "-c",
                        str(source),
                        "-o",
                        str(temp_obj),
                    ],
                    check=True,
                )
                temp_dep.replace(depfile)
                temp_obj.replace(obj)
        objects.append(obj)
    link_flags = shlex.split(args.ldflags)
    # Explicit archives and linker scripts are link inputs, even when no C changes.
    extra_inputs = [Path(x) for x in link_flags if Path(x).is_file()]
    for x in link_flags:
        if "--version-script=" in x:
            extra_inputs.append(Path(x.split("--version-script=", 1)[1]))
    signature = {
        "objects": [[str(p), p.stat().st_mtime_ns] for p in objects],
        "flags": link_flags,
        "inputs": [[str(p), p.stat().st_mtime_ns] for p in extra_inputs],
    }
    state = args.output.with_suffix(args.output.suffix + ".build.json")
    if args.output.exists() and state.exists():
        previous = json.loads(state.read_text())
        if (
            previous.get("inputs") == signature
            and previous.get("output_sha256")
            == hashlib.sha256(args.output.read_bytes()).hexdigest()
        ):
            print("Up to date: " + str(args.output))
            return
    print("Link " + str(args.output), flush=True)
    # Same filesystem as destination: final publication is an atomic rename.
    with tempfile.TemporaryDirectory(prefix=".ctools-link-", dir=args.output.parent) as tmp:
        output = Path(tmp) / args.output.name
        subprocess.run(
            [*compiler, *flags, "-o", str(output), *map(str, objects), *link_flags], check=True
        )
        output.replace(args.output)
    # JSON round trip normalizes tuples so a no-op invocation compares equal.
    temporary = state.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(
            {
                "inputs": signature,
                "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    temporary.replace(state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("platform", "cc", "cflags", "ldflags", "sources"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--objects", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))


if __name__ == "__main__":
    main()
