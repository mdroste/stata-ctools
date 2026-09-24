"""Build the pinned static OpenMP runtime against the ctools macOS baseline."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import tarfile
import urllib.request

VERSION = "21.1.8"
SHA256 = {
    "openmp": "856b023748b41ac7b2c83fd8e9f765ff48a4df2fe6777d2811ef7c7ed8f2f977",
    "cmake": "85735f20fd8c81ecb0a09abb0c267018475420e93b65050cc5b7634eab744de9",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arch", choices=["arm64", "x86_64"])
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
    build = work / f"build-{args.arch}"
    subprocess.run(["cmake", "-S", str(work / f"openmp-{VERSION}.src"), "-B", str(build),
                    "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_OSX_ARCHITECTURES={args.arch}",
                    "-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0", f"-DCMAKE_INSTALL_PREFIX={args.prefix.resolve()}",
                    "-DLIBOMP_ENABLE_SHARED=OFF", "-DOPENMP_ENABLE_LIBOMPTARGET=OFF",
                    "-DLIBOMP_OMPT_SUPPORT=OFF"], check=True)
    subprocess.run(["cmake", "--build", str(build), "--parallel", "4"], check=True)
    subprocess.run(["cmake", "--install", str(build)], check=True)


if __name__ == "__main__":
    main()
