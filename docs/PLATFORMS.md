# Platform contract

The distribution plugins use Stata's 3.0 plugin interface and require Stata 14.1 or later. Stata's plugin API limits observations to 2^31−1. The minimum Stata version is an interface requirement; a successful run on a newer Stata version does not certify every older version.

| Plugin | CPU baseline | OS/runtime baseline | External dependencies |
|---|---|---|---|
| macOS ARM | Apple Silicon (M1) | macOS 11 | System Accelerate; static LLVM OpenMP bundled |
| macOS Intel | x86-64, without AVX2/FMA requirements | macOS 11 | System Accelerate; static LLVM OpenMP bundled |
| Linux | x86-64, without AVX2/FMA requirements | glibc 2.35 (Ubuntu 22.04) or newer | `libgomp.so.1` (Ubuntu package `libgomp1`) |
| Windows | x86-64, without AVX2/FMA requirements | Windows 10 or newer | System DLLs; GCC/OpenMP/winpthreads linked statically |

These are build and packaging baselines. CI runs on its listed hosted runners and checks imported libraries and deployment versions; it does not emulate an old CPU or certify execution on every minimum OS. Test release candidates on those minimum systems before claiming runtime coverage. Custom builds using OpenBLAS or a dynamic llvm-mingw OpenMP runtime fall outside this distribution contract and must supply their dependencies.

## macOS runtime build

An installed Homebrew OpenMP archive may target a newer macOS release than the plugin. Build the pinned LLVM runtime with the same deployment target (Python 3.12+, CMake, and a C/C++ compiler are needed):

```bash
python3 scripts/build_macos_openmp.py arm64 /tmp/ctools-libomp --work /tmp/ctools-llvm
make macos-arm LIBOMP_PREFIX=/tmp/ctools-libomp
python3 validation/check_dependencies.py macos build/ctools_mac_arm.plugin --arch arm64
```

For Intel, substitute `x86_64` and `macos-intel`. The helper verifies source archive hashes and builds OpenMP 21.1.8 for macOS 11. Both Make targets check the archive's architecture and deployment version before linking. With no static archive, a local build uses pthreads without OpenMP; it never silently adds a Homebrew dylib dependency. Distribution CI builds the static runtime explicitly.

Linux distribution CI uses Ubuntu 22.04 and disables optional OpenBLAS. `USE_OPENBLAS=yes` enables detection for custom builds. Windows distribution CI uses MSYS2 GCC; dependency validation rejects external compiler DLLs.
