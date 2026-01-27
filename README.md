# ctools

[![PluginsSome ](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

Extremely fast drop-in replacements for a variety of Stata commands.

**This project is pre-release. It has not been exhaustively tested. Please report any bugs you encounter on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** programs inherit the syntax and functionality of the programs they replace, but are usually much faster for large datasets.

| Stata command | Replaced with | Description | Typical speedup |
| --- | --- | --- | ---: |
| `import` | `cimport` | Import text-delimited and Excel data | **30x** |
| `export` | `cexport` | Export text-delimited and Excel data | **25x** |
| `sort` | `csort` | Sort dataset | **1-5x** |
| `merge` | `cmerge` | Merge (join) datasets | **1-5x** |
| `sample` | `csample` | Resampling without replacement | **2-5x** |
| `bsample` | `cbsample` | Resampling with replacement | **2-5x** |
| `encode` | `cencode` | Recast string as labeled numeric | **10-20x** |
| `decode` | `cdecode` | Recast labeled numeric as string | **10-20x** |
| `destring` | `cdestring` | Recast string as numeric type | **10-20x** |
| `gstats winsor` | `cwinsor` | Winsorize variables | **2-10x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **10-40x** |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects | **10-30x** |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects | **10-30x** |
| `qreg` | `cqreg` | Quantile regression | **2-4x** |

Some ctools programs have extended functionality. For instance, `cbinscatter` supports multi-way (high-dimensional) fixed effects and the alternative procedure to control for covariates described in [Cattaneo et al. (2024)](https://www.aeaweb.org/articles?id=10.1257/aer.20221576); `csort` allows the user to select one of several different implemented parallelized sorting algorithms. Internal help files for each ctools program provide complete documentation.

Speedups will vary depending on a lot of factors. Most commands will run appreciably on pretty much any dataset. On the other hand, `csort` and `cmerge` involve a lot of overhead (needing to read entire datasets from Stata to the C plugin and back), and these commands can actualy be *slower* than built-in `sort`/`merge` if the datset is sufficiently wide (many variables). If your dataset has at most a few dozen variables and many millions of observations, `csort` and `cmerge` will probably speed things up.



## Compatibility and Requirements

ctools is compatible with Stata 16.0+. It does not require any dependencies.


## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

### Manual Installation

1. Download or clone this repository
2. Copy the contents of the `build/` directory to a location on Stata's adopath (e.g., your personal ado directory)

These installations will download all of the Stata program and documentation files (.ado and .sthlp) and compiled plugins for Linux, Windows, and Mac. The appropriate (operating system and architecture-dependent) compiled plugin is automatically invoked by ctools.


## Building from Source Files


You probably do not need to compile the ctools plugin syourself; GitHub automatically builds plugins for Windows, Mac (Intel and ARM/M-series), and Linux. If you want to build plugins from source, make sure OpenMP is available on your system.

```bash
make              # Build for current platform
make all          # Build all platform plugins
make check        # Check build dependencies
make clean        # Remove compiled files
```

If you have a workstation/server CPU with lots of cache memory, you might want to try playing around with the settings in [src/ctools_config.h](src/ctools_config.h).

See [DEVELOPERS.md](DEVELOPERS.md) for additional information on ctools' architecture and core logic.


## Usage Notes and Limitations

- All ctools programs follow a basic structure: (1) copy Stata's data in memory into C structures; (2) operate on those structures; (3) return data from Stata to C. *This means that ctools programs require more memory than the programs they replace.* In addition, some of these commands will run faster when they involve fewer variables, or when you have fewer variables in memory. For instance, the runtime of `csort` is heavily dependent on the number of variables in memory, and can be slower Stata's built-in `sort` if the dataset is relatively wide (e.g. 100+ variables) due to this data transfer overhead.
- ctools does not work with datasets exceeding 2^31 observations (~2.147 billion obs). This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of the Stata function interface for C plugins and can only be addressed with an update to the Stata Plugin Interface. 


## Authorship

99.9% of the code in this repository was written by Claude Opus 4.5 (through Claude Code), with some debugging and refactoring assistance from OpenAI GPT 5.2 (through Codex).


## Thanks

- [Sergio Correia](https://github.com/sergiocorreia) for [ftools](https://github.com/sergiocorreia/ftools), [reghdfe](https://github.com/sergiocorreia/reghdfe), and [ivreghdfe](https://github.com/sergiocorreia/ivreghdfe) (with [Lars Vilhuber](https://www.ilr.cornell.edu/people/lars-vilhuber))
- [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for [gtools](https://github.com/mcaceresb/gtools) and invaluable contributions to [cowsay](https://github.com/mdroste/stata-cowsay) 
- [Christopher (Kit) Baum](https://www.bc.edu/bc-web/schools/morrissey/departments/economics/people/faculty-directory/christopher-baum.html), [Mark E Schaffer](https://www.hw.ac.uk/profiles/uk/school/ebs/faculty/mark-schaffer), and [Steven Stillman](https://www.unibz.it/en/faculties/economics-management/academic-staff/person/36390-steven-stillman) for [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html)
- [Sascha Witt](https://github.com/SaschaWitt) for the [In-place Parallel Super Scalar Samplesort (IPS⁴o)](https://github.com/SaschaWitt/ips4o) sorting algorithm.
- Claude Code


## License

This program is MIT-licensed. 


## Contributing

Contributions are welcome. Please open an [Issue](https://github.com/mdroste/stata-ctools/issues) or submit a pull request.
