# Platform contract

The distribution plugins use Stata's 3.0 plugin interface and require Stata 14.1 or later. Stata's plugin API limits observations to 2^31−1. The minimum Stata version is an interface requirement; a successful run on a newer Stata version does not certify every older version.

| Plugin | CPU baseline | OS/runtime baseline | External dependencies |
|---|---|---|---|
| macOS ARM | Apple Silicon (M1) | macOS 11 | System Accelerate; static LLVM OpenMP bundled |
| macOS Intel | x86-64, without AVX2/FMA requirements | macOS 11 | System Accelerate; static LLVM OpenMP bundled |
| Linux | x86-64, without AVX2/FMA requirements | glibc 2.17 (manylinux2014, as on RHEL/CentOS 7) or newer | glibc only; static LLVM OpenMP bundled (no `libgomp`) |
| Windows | x86-64, without AVX2/FMA requirements | Windows 10 or newer | System DLLs; GCC/OpenMP/winpthreads linked statically |

These are build and packaging baselines. CI runs on its listed hosted runners and checks imported libraries and deployment versions; it does not emulate an old CPU or certify execution on every minimum OS. Test release candidates on those minimum systems before claiming runtime coverage. Custom builds using OpenBLAS or a dynamic llvm-mingw OpenMP runtime fall outside this distribution contract and must supply their dependencies.

## macOS runtime build

An installed Homebrew OpenMP archive may target a newer macOS release than the plugin. Build the pinned LLVM runtime with the same deployment target (Python 3.12+, CMake, and a C/C++ compiler are needed):

```bash
python3 scripts/build_macos_openmp.py arm64 /tmp/ctools-libomp --work /tmp/ctools-llvm
make macos-arm LIBOMP_PREFIX=/tmp/ctools-libomp
python3 validation/run.py --checks dependencies macos build/ctools_mac_arm.plugin --arch arm64
```

For Intel, substitute `x86_64` and `macos-intel`. The helper verifies source archive hashes and builds OpenMP 21.1.8 for macOS 11. Both Make targets check the archive's architecture and deployment version before linking. With no static archive, a local build uses pthreads without OpenMP; it never silently adds a Homebrew dylib dependency. Distribution CI builds the static runtime explicitly.

Linux distribution CI compiles inside the `quay.io/pypa/manylinux2014_x86_64` container (`make linux-docker`, CentOS 7's glibc 2.17 and devtoolset GCC) and disables optional OpenBLAS. `make linux-openmp` first builds the pinned LLVM OpenMP 21.1.8 runtime there as a static, position-independent archive (`scripts/build_linux_openmp.py`, cached by CI), and the plugin links it in place of `libgomp`: GCC's `GOMP_*` calls resolve to libomp's GOMP layer inside the plugin, and the version script keeps the runtime's symbols local. `python3 validation/run.py --checks dependencies linux build/ctools_linux.plugin` then rejects a plugin that needs a symbol version newer than `GLIBC_2.17`, needs any `GOMP`/`OMP` symbol version or `libgomp.so.1`, does not link `libm.so.6`, exports anything besides `stata_call` and `pginit`, or has an undefined symbol bound to no library version. `make linux-docker-load` loads the plugin with every symbol bound (`RTLD_NOW`) and a clean environment in manylinux2014, plain `centos:7` and `rockylinux:8`, and fails if `ldd` resolves an OpenMP runtime. A plugin built directly on a newer distribution (`make linux`, without `LIBOMP_LINUX`) links that host's `libgomp.so.1` and fails the dependency check. `USE_OPENBLAS=yes` enables detection for custom builds. Windows distribution CI uses MSYS2 GCC; dependency validation rejects external compiler DLLs.

macOS branch and release builds require Developer ID Application signing and
Apple notarization before packaging; see [macOS signing](MACOS_SIGNING.md) for
Apple credentials, GitHub environment configuration, and online ticket lookup.

Windows branch and release builds also require Azure Authenticode signing before
validation and packaging; see [Windows signing](WINDOWS_SIGNING.md) for the Azure
and GitHub environment configuration. Pull-request candidates remain unsigned.

## Worker and scratch limits

`threads(0)` uses the configured default; positive requests accept 1–256 workers and are clamped to the available CPUs and runtime limits. Without an override, the default cap is 16. See the [thread policy](README.md#thread-count-management). Sampling limits workers to available groups; range statistics allocate sorting scratch only for percentiles, the median, and the interquartile range. Both commands reduce concurrent sorting/shuffle work to a 256 MiB scratch target. One worker remains available when a single group needs more than that target; this is not a limit on total command memory.
