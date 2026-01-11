# ctools

[![PluginsSome ](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

High-performance C-accelerated drop-in replacements for Stata commands.

## Overview

**ctools** is a set of drop-in replacements for Stata programs. These programs inherit the syntax and functionality of each program that they replace, but are typically much faster for large datasets.

**ctools** includes replacements for the following Stata programs:

| Stata command | ctools command | Purpose | Typical speedup |
| --- | --- | --- | ---: |
| `sort` | `csort` | Sort dataset | **1-5x** |
| `import delimited` | `cimport delimited` | Import CSV/TSV files | **50x** |
| `export delimited` | `cexport delimited` | Export CSV/TSV files | **25x** |
| `merge` | `cmerge` | Merge (join) datasets | **2-5x** |
| `reghdfe` | `creghdfe` | High-dimensional FE regression | **10-20x** |
| `ivreghdfe` | `civreghdfe` | High-dimensional FE IV (2SLS) regression | **10-20x** |
| `qreg` | `cqreg` | Quantile regression | **4x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **30x+** |
| `ivreghdfe` | `civreghdfe` | IV regression with HDFE | **10-20x** |

Each of these programs is comprised of a wrapper ado file and plugin written in C. 


## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build")
```

### Manual Installation

1. Download or clone this repository
2. Copy the contents of the `build/` directory to a location on Stata's adopath (e.g., your personal ado directory)

The package automatically detects your operating system and architecture, loading the appropriate precompiled plugin.


## Building from Source Files

You almost certainly do not need to compile this program from the source files, as pre-compiled plugins are automatically built by GitHub and are included with the installation methods described above. If you want to build from source, then make sure your compiler can use OpenMP.

```bash
make              # Build for current platform
make all          # Build all platform plugins
make check        # Check build dependencies
make clean        # Remove compiled files
```

## Authorship of this Program

Ctools contains >30k lines of code C code and ~4k lines of Stata code; 99.99% of this was written by Claude Code with Opus 4.5. 


## Performance Notes

- All commands use OpenMP for parallel computation when available.
- All commands attempt to use the same few tricks to speed up computation where possible: 8-way or 16-way loop unrolling, SIMD vectorisation and OpenMP parallelization where appropriate. 
- ctools works by transporting a copy of your data into C (and, depending on the program, possibly from C back to Stata). This means that you need a comfortable amount of memory (RAM) to run these programs. 
- If your CPU has an exceptionally large L3 cache, you might try to mess with ctools_config.h to see if more aggressive prefetching buys you speed. 

## Compatibility

- Stata 14.0 or later
- Windows (x64), macOS (Intel and Apple Silicon), Linux (x64)

## License

This program is MIT-licensed. 

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.
