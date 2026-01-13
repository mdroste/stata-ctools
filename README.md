# ctools

[![PluginsSome ](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

High-performance C-accelerated drop-in replacements for Stata commands.

** This is an initial release. This software has not been widely tested. Please report any bugs you encounter on Github Issues.**

## Overview

**ctools** is a set of drop-in replacements for a variety of Stata programs. **ctools** programs inherit the syntax and functionality of the programs they replace, but are usually much faster for large datasets.

**ctools** includes replacements for the following Stata programs:

| Stata command | ctools command | Purpose | Typical speedup |
| --- | --- | --- | ---: |
| `import delimited` | `cimport delimited` | Import text-delimited files | **40x** |
| `export delimited` | `cexport delimited` | Export text-delimited files | **25x** |
| `sort` | `csort` | Sort dataset | **1-6x** |
| `merge` | `cmerge` | Merge (join) datasets | **2-8x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **30x+** |
| `reghdfe` | `creghdfe` | OLS with high-dimensional fixed effects | **10-20x** |
| `ivreghdfe` | `civreghdfe` | 2SLS with high-dimensional fixed effects | **10-20x** |
| `qreg` | `cqreg` | Quantile regression | **4x** |

Each of these programs is implemented with its own Stata .ado file (and internal .sthlp documentation).


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

## Usage Notes

- Some of these programs have I/O overhead involved in transporting data from Stata to C. This matters for you in two ways. First, it means that using --ctools-- requires more memory (RAM) than built-in Stata command, particularly for -csort- and -cmerge- which tend to require loading all variables / many variables into memory (cqreg, creghdfe, binscatter etc. only pull in what is needed). In addition, this IO overhead means that -csort- and -cmerge- in particular will run faster when there are fewer variables in memory. 
- The Stata function interface is limited to working with datasets with fewer than 2^32-1 observations (about 2.147 billion). This is a [known limitation]() of the Stata function interface for C plugins. 
- csort includes an additional option (algorithm) allowing you to choose a handful of sorting algorithms, as described in the internal help.

## Compatibility

- Stata 14.0+

## Authorship

99.9% of the code for ctools repository was written using Claude Opus 4.5 and Claude Code. by Claude Opus 4.5 using Claude Code.

## Thanks

Thanks to [Sergio Correa](https://github.com/sergiocorrea) for creating -[ftools](https://github.com/sergiocorrea/ftools)- and [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for creating -[gtools](https://github.com/mcaceresb/gtools)-. 


## License

This program is MIT-licensed. 

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.
