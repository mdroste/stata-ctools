"""Build the pinned static OpenMP runtime for the Linux plugin.

Run inside the manylinux2014 container (`make linux-openmp`). The archive is
position-independent so it can be linked into the plugin, and GCC's GOMP_*
calls resolve to libomp's GOMP compatibility layer by their plain names, so
the plugin needs no libgomp.so.1 at run time.
"""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import urllib.request

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_macos_openmp import SHA256, VERSION  # noqa: E402  (one pinned runtime for both)

CMAKE_FALLBACK = "3.31.6"


def cmake():
    found = shutil.which("cmake")
    if found:
        return found
    subprocess.run([sys.executable, "-m", "pip", "install", "--user", f"cmake=={CMAKE_FALLBACK}"], check=True)
    found = shutil.which("cmake") or str(Path.home() / ".local/bin/cmake")
    if not Path(found).exists():
        raise SystemExit("cmake is not available")
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    for project, digest in SHA256.items():
        name = f"{project}-{VERSION}.src.tar.xz"
        archive = work / name
        if not archive.exists():
            urllib.request.urlretrieve(f"https://github.com/llvm/llvm-project/releases/download/llvmorg-{VERSION}/{name}", archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != digest:
            raise SystemExit(f"Checksum mismatch: {archive}")
        with tarfile.open(archive) as tar:
            tar.extractall(work, filter="data")
    cmake_link = work / "cmake"
    if not cmake_link.exists():
        cmake_link.symlink_to(f"cmake-{VERSION}.src", target_is_directory=True)
    build = work / "build-linux-x86_64"
    tool = cmake()
    subprocess.run([tool, "-S", str(work / f"openmp-{VERSION}.src"), "-B", str(build),
                    "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_C_COMPILER=gcc", "-DCMAKE_CXX_COMPILER=g++",
                    f"-DCMAKE_INSTALL_PREFIX={args.prefix.resolve()}",
                    "-DCMAKE_POSITION_INDEPENDENT_CODE=ON", "-DLIBOMP_ENABLE_SHARED=OFF",
                    "-DLIBOMP_USE_VERSION_SYMBOLS=OFF", "-DLIBOMP_INSTALL_ALIASES=OFF",
                    "-DOPENMP_ENABLE_LIBOMPTARGET=OFF", "-DLIBOMP_OMPT_SUPPORT=OFF",
                    "-DLIBOMP_USE_HWLOC=OFF"], check=True)
    subprocess.run([tool, "--build", str(build), "--parallel", "4"], check=True)
    subprocess.run([tool, "--install", str(build)], check=True)
    if not (args.prefix / "lib/libomp.a").exists():
        raise SystemExit(f"{args.prefix}/lib/libomp.a was not installed")


if __name__ == "__main__":
    main()
