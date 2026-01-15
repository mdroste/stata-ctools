# ctools

[![PluginsSome ](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

High-performance C-accelerated drop-in replacements for Stata commands.

** This is an initial release. This software has not been widely tested. Please report any bugs you encounter on Github Issues.**

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

If your computing environment has a CPU with an exceptionally large L3 cache (e.g. AMD's X3D chips or a Threadripper), you might explore playing around with the settings in ./src/ctools_config.h.

## Usage Notes

- Some of these programs have I/O overhead involved in transporting data from Stata to C. This means -ctools- commands require more memory (RAM) than their Stata/Mata counterparts, as C needs to copy at least some of your data into arrays visible to C. This also means that some commands will run faster when they involve fewer variables (all else equal), or when you have fewer variables in memory. Only the variables required to run a command are pulled into C. For instance, -cbinscatter y x- will only load the variables y and x into memory, and the remainder of your dataset is irrelevant. On the other hand, -csort- and -cmerge- generally require pulling a lot of data into memory (and back); fixing N, the IO overhead of csort and cmerge scale linearly with the width of your data.
- ctools only works with datasets smaller than 2^31-1 observations (about 2.147 billion). This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of the Stata function interface for C plugins and can only be addressed if Stata updates this interface.

## Compatibility

- Stata 14.0+

## Authorship

99.9% of the code for ctools repository was written using Claude Opus 4.5 and Claude Code.

## Thanks

Thanks to [Sergio Correa](https://github.com/sergiocorrea) for creating [ftools](https://github.com/sergiocorrea/ftools)/ and [reghdfe](https://github.com/sergiocorreia/reghdfe), and thanks to [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for creating [gtools](https://github.com/mcaceresb/gtools). Thanks also to Christopher (Kit) Baum for creating ivreg2, from which both ivreghdfe and civreghdfe are based. 


## License

This program is MIT-licensed. 

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.
