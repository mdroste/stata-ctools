# ctools

[![PluginsSome ](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

High-performance C-accelerated drop-in replacements for Stata commands.

**This is an initial release. This software has not been widely tested. Please report any bugs you encounter on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** is a set of drop-in replacements for a variety of Stata programs. **ctools** programs inherit the syntax and functionality of the programs they replace, but are usually much faster for large datasets.

| Stata command | ctools command | Purpose | Typical speedup |
| --- | --- | --- | ---: |
| `import delimited` | `cimport delimited` | Import text-delimited files | **40x** |
| `export delimited` | `cexport delimited` | Export text-delimited files | **40x** |
| `sort` | `csort` | Sort dataset | **1-6x** |
| `merge` | `cmerge` | Merge (join) datasets | **1-8x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **2-30x+** |
| `reghdfe` | `creghdfe` | OLS with high-dimensional fixed effects | **10-20x** |
| `ivreghdfe` | `civreghdfe` | IV, GMM, etc. with high-dimensional fixed effects | **10-20x** |
| `qreg` | `cqreg` | Quantile regression | **4x** |


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


You probably do not need to compile the ctools C plugin yourself; GitHub automatically compiles plugins for Windows, Mac (Intel and ARM/M-series), and Linux platforms, and installation will automatically download and select the appropriate plugin for you. If you want to build from source, then make sure your compiler can use OpenMP.

```bash
make              # Build for current platform
make all          # Build all platform plugins
make check        # Check build dependencies
make clean        # Remove compiled files
```

A lot of ctools operations appear to be memory bandwidth bound. If your CPU has a very large L3 cache, you might explore playing around with the settings in ./src/ctools_config.h.

You are welcome to use ctools however you like. For instance, ./src/ctools_data_io.c implements what I think are pretty efficient methods to transport data from Stata to C and vice versa; if you want to write a Stata C plugin, you do not need to reinvent the wheel, and the "common" ctools functions are intended to be modular and portable to other projects.

See DEVELOPERS.md for additional information on ctools' architecture and core logic.

## Threading Model (quick note)

ctools uses OpenMP for tight, regular loops (sorting, linear algebra, HDFE work) and a global persistent thread pool for per-variable I/O and other repeated batch work (to avoid thread creation overhead). See DEVELOPERS.md for a detailed list of where each is used and tuning suggestions.

## Usage Notes and Limitations

- All ctools programs work by copying Stata data in memory, operating on it, and returning results. Because these programs do not operate on Stata datasets in place, -ctools- requires more memory (RAM) than their Stata/Mata-coded counterpart commands. In addition, some of these commands will run faster when they involve fewer variables, or when you have fewer variables in memory. As an extreme case, -csort- requires loading the entire dataset into memory (not just the 'sort key' variables). This means that -csort-'s runtime is heavily dependent on the number of variables in your current dataset (frame), and likewise with -cmerge-. In contrast, other commands like civreghdfe, creghdfe, cbinscatter, and cqreg will only load those variables that are necessary, and do not generally write very much data back to Stata, and this I/O concern is minimal.
- ctools only works with datasets smaller than 2^31-1 observations (~2.147 billion obs). This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of the Stata function interface for C plugins and can only be addressed if Stata updates their API for interfacing with C plugins. 

## Compatibility

- Stata 16.0+ for full compatibility; Stata 14.0+ probably works (untested).

## Authorship

99.9% of the code for ctools repository was written using Claude Opus 4.5 and Claude Code.

## Thanks

Thanks to [Sergio Correia](https://github.com/sergiocorreia) for creating [ftools](https://github.com/sergiocorreia/ftools) and [reghdfe](https://github.com/sergiocorreia/reghdfe), and thanks to [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for creating [gtools](https://github.com/mcaceresb/gtools). Thanks also to Christopher (Kit) Baum for creating ivreg2, from which both ivreghdfe and civreghdfe are based. Lastly, thanks to Sascha Witt for the "In Place Parallel Scalar Supersort (ips4o) sorting algorithm that works really well. 


## License

This program is MIT-licensed. 

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.
