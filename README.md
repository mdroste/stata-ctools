  ```
   █████╗██████╗ █████╗  █████╗ ██╗   ██████╗
  ██╔═══╝╚═██╔═╝██╔══██╗██╔══██╗██║   ██╔═══╝
  ██║      ██║  ██║  ██║██║  ██║██║   ██████╗
  ██║      ██║  ██║  ██║██║  ██║██║   ╚═══██║
  ╚█████╗  ██║  ╚█████╔╝╚█████╔╝█████╗██████║
   ╚════╝  ╚═╝   ╚════╝  ╚════╝ ╚════╝╚═════╝
   some really fast Stata programs
  ```

[![Build](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)
![Version](https://img.shields.io/badge/version-1.0.1_(February_2026)-blue)

**This is an initial release. Please report problems/suggestions on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** programs inherit the syntax and functionality of the programs they replace, but are usually much faster for large datasets.

| Stata command | Replaced with | Description | Typical speedup |
| --- | --- | --- | ---: |
| `import` | `cimport` | Import text-delimited and Excel data | **10-30x** |
| `export` | `cexport` | Export text-delimited and Excel data | **10-30x** |
| `sort` | `csort` | Sort dataset | **1-5x** |
| `merge` | `cmerge` | Merge (join) datasets | **1-5x** |
| `sample` | `csample` | Resampling without replacement | **3-4x** |
| `bsample` | `cbsample` | Resampling with replacement | **3-4x** |
| `encode` | `cencode` | Recast string as labeled numeric | **2-10x** |
| `decode` | `cdecode` | Recast labeled numeric as string | **2-10x** |
| `destring` | `cdestring` | Recast string as numeric type | **5-20x** |
| `qreg` | `cqreg` | Quantile regression | **2-4x** |
| `gstats winsor` | `cwinsor` | Winsorize variables | **2-10x** |
| `rangestat` | `crangestat` | Range statistics of variables | **20-100x** |
| `psmatch2` | `cpsmatch` | Propensity score matching | **10-30x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **10-30x** |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects | **10-30x** |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects | **10-30x** |

Some ctools programs have extended functionality. For instance, `cbinscatter` supports multi-way fixed effects and the procedure to control for covariates characterized by [Cattaneo et al. (2024)](https://www.aeaweb.org/articles?id=10.1257/aer.20221576). See [FEATURES.md](./FEATURES.md) for a brief description of the new features implemented for each command above. Each command also has an associated internal help file (e.g. `help cbinscatter`).

Most ctools commands (e.g. `creghdfe`, `cbinscatter`) will be much faster on most (or maybe all) datasets. On the other hand, `csort` and `cmerge` involve a lot of overhead (needing to read entire datasets from Stata to the C plugin and back), and themselves replace internal compiled Stata routines (`sort`, `merge`) that are only inefficient when datasets become fairly long (millions of obs). As a result, it is possible for `csort` or `cmerge` to be *slower* than  `sort`/`merge` if your dataset is sufficiently wide (many variables) or not very long. If your dataset has at most a few dozen variables and many millions of observations, `csort` and `cmerge` will probably be significantly faster.


## Compatibility and Requirements

ctools is compatible with Stata 14.1+ (plugin interface version 3.0). It does not require any dependencies.

**Note:** ctools does not support datasets exceeding 2^31 (~2.147 billion) observations. This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of Stata's plugin interface for C, and can only be addressed if this interface is updated. Most commands also do not support operating on `strL` variable types. 



## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

### Manual Installation

1. Download or clone this repository
2. Copy the contents of this repository's `build/` folder somewhere  Stata's `adopath` (e.g., your personal ado directory)

These installations will download all of the Stata program and documentation files (.ado and .sthlp) and compiled plugins for Linux, Windows, and Mac. The appropriate (operating system and architecture-dependent) compiled plugin is automatically invoked by ctools.

An installation path via `ssc install` will be provided soon.


## Building from Source

You probably do not need to compile the ctools plugin yourself. GitHub automatically builds plugins for Windows, Mac, and Linux when this repository is updated, and these pre-built plugins are included with standard installation. If you want to build these plugins yourself, make sure OpenMP is available on your system. 

```bash
make              # Build for current platform
make all          # Build all platform plugins
make check        # Check build dependencies
make clean        # Remove compiled files
```

If you have a workstation/server CPU with lots of cache, you might want to try playing around with the settings in [src/ctools_config.h](src/ctools_config.h).

See [DEVELOPERS.md](./DEVELOPERS.md) for additional information on ctools' architecture and core logic.


## Usage Notes

- ctools programs follow a simple structure: (1) copy data from Stata to C; (2) operate on that data (3) return data from C to Stata. Because ctools copies data, it generally requires more memory to run ctools programs than the programs they replace.
- The commands `csort` and `cmerge` involve a lot of overhead when the number of variables is large relative to the number of observations because data transfer overhead is very significant for these two programs. They can be slower Stata's built-in `sort`/`merge` if datasets are very wide (e.g. 100+ variables) and/or not very long (e.g. <1M obs). They're probably faster if your dataset has many millions of observations and fewer 2-3 dozen variables. Your mileage may vary. 


## Authorship and Development

99.9% of the code in this repository was written by Claude Opus 4.5 (through Claude Code), with some debugging and refactoring assistance from OpenAI GPT 5.2 (through Codex).


## Thanks

- [Sergio Correia](https://github.com/sergiocorreia) for [ftools](https://github.com/sergiocorreia/ftools), [reghdfe](https://github.com/sergiocorreia/reghdfe), and [ivreghdfe](https://github.com/sergiocorreia/ivreghdfe) (with [Lars Vilhuber](https://www.ilr.cornell.edu/people/lars-vilhuber))
- [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for [gtools](https://github.com/mcaceresb/gtools) and invaluable contributions to [cowsay](https://github.com/mdroste/stata-cowsay)
- [Christopher (Kit) Baum](https://www.bc.edu/bc-web/schools/morrissey/departments/economics/people/faculty-directory/christopher-baum.html), [Mark E Schaffer](https://www.hw.ac.uk/profiles/uk/school/ebs/faculty/mark-schaffer), and [Steven Stillman](https://www.unibz.it/en/faculties/economics-management/academic-staff/person/36390-steven-stillman) for [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html)
- [Robert Picard](https://ideas.repec.org/f/ppi320.html), [Nicholas J. Cox](https://ideas.repec.org/e/pco34.html), and Roberto Ferrer for [rangestat](https://ideas.repec.org/c/boc/bocode/s458161.html)
- [Sascha Witt](https://github.com/SaschaWitt) for the [In-place Parallel Super Scalar Samplesort (IPS⁴o)](https://github.com/SaschaWitt/ips4o) sorting algorithm.
- [Libdeflate](https://github.com/ebiggers/libdeflate) for really good DEFLATE-based compression/decompression
- Claude Code


## License

This program is MIT-licensed. 


## Contributing

Contributions are welcome. Please open an [Issue](https://github.com/mdroste/stata-ctools/issues) or submit a pull request.
