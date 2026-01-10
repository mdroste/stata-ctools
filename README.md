# ctools

[![Build Plugins](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

High-performance C-accelerated drop-in replacements for Stata commands.

## Overview

**ctools** provides drop-in replacements for common Stata commands that run substantially faster on large datasets. All commands are implemented in optimized C with OpenMP parallelization and compile to native plugins for Windows, macOS (Intel and Apple Silicon), and Linux.

| Stata command | ctools command | Purpose | Typical speedup |
| --- | --- | --- | ---: |
| `sort` | `csort` | Sort dataset | **2-4x** |
| `import delimited` | `cimport delimited` | Import CSV/TSV files | **50x** |
| `export delimited` | `cexport delimited` | Export CSV/TSV files | **25x** |
| `merge` | `cmerge` | Merge (join) datasets | **2-5x** |
| `reghdfe` | `creghdfe` | High-dimensional FE regression | **10x** |
| `qreg` | `cqreg` | Quantile regression | **4x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **10x+** |
| `ivreghdfe` | `civreghdfe` | IV regression with HDFE | **10x** |

## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build")
```

### Manual Installation

1. Download or clone this repository
2. Copy the contents of the `build/` directory to a location on Stata's adopath (e.g., your personal ado directory)

The package automatically detects your operating system and architecture, loading the appropriate precompiled plugin.


## Building from Source

You almost certainly do not need to compile this program from the source files, as pre-compiled plugins are automatically built by GitHub and installeed when you follow the installation instructions above. If you want to, then it is highly recommended that your compiler work with OpenMP. The Makefile should automatically detect your platforms and compile with OpenMP when available.

### Requirements

| Platform | Requirements |
| --- | --- |
| macOS | Xcode Command Line Tools, Homebrew with `libomp` |
| Windows | MSVC/MinGW, or cross-compile with `llvm-mingw` |
| Linux | GCC with OpenMP (`libgomp`) |

### Build Commands

```bash
make              # Build for current platform
make all          # Build all platform plugins
make check        # Check build dependencies
make clean        # Remove compiled files
```

### Output Files

| File | Platform |
| --- | --- |
| `build/ctools_mac_arm.plugin` | macOS Apple Silicon |
| `build/ctools_mac_x86.plugin` | macOS Intel |
| `build/ctools_windows.plugin` | Windows x64 |
| `build/ctools_linux.plugin` | Linux x64 |

## Project Structure

```
ctools/
├── src/                        # C source files
│   ├── stplugin.c/h            # Stata Plugin Interface
│   ├── ctools_plugin.c         # Main dispatcher
│   ├── ctools_data_io.c        # Parallel data load/store
│   ├── ctools_sort_radix_lsd.c # Parallel LSD radix sort
│   ├── csort/                  # Sort command
│   ├── cmerge/                 # Merge command
│   ├── cimport/                # CSV import
│   ├── cexport/                # CSV export
│   ├── creghdfe/               # HDFE regression
│   ├── cqreg/                  # Quantile regression
│   ├── cbinscatter/            # Binned scatter plots
│   └── civreghdfe/             # IV regression with HDFE
├── build/                      # Compiled plugins, .ado, .sthlp files
├── benchmark/                  # Performance test .do files
├── Makefile                    # Build system
└── README.md
```

## Performance Notes

- All commands use OpenMP for parallel computation when available
- Data loading uses 8-way loop unrolling for cache efficiency
- Sorting uses a parallel LSD radix sort optimized for Stata's data types
- `creghdfe` uses a conjugate gradient solver with Kaczmarz-style sweeping
- `cqreg` uses a primal-dual Interior Point Method with Mehrotra predictor-corrector

## Compatibility

- Stata 14.0 or later
- Windows (x64), macOS (Intel and Apple Silicon), Linux (x64)
- OpenMP support recommended but not required

## License

MIT License

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.

## Citation

If you use ctools in your research, please cite:

```bibtex
@software{ctools,
  title = {ctools: High-performance C-accelerated commands for Stata},
  url = {https://github.com/mdroste/stata-ctools}
}
```
