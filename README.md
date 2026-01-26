# ctools

[![PluginsSome ](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)

Fast drop-in replacements for a variety of Stata commands.

**This is an initial release. This software has not been widely tested. Please report any bugs you encounter on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** programs inherit the syntax and functionality of the programs they replace, but are usually much faster for large datasets.

| Stata command | Replaced with | Description | Typical speedup |
| --- | --- | --- | ---: |
| `import delimited` | `cimport delimited` | Import text-delimited files | **30x** |
| `export delimited` | `cexport delimited` | Export text-delimited files | **25x** |
| `sort` | `csort` | Sort dataset | **1-5x** |
| `merge` | `cmerge` | Merge (join) datasets | **1-5x** |
| `encode`/`decode` | `cencode`/`cdecode` | Recast encoded string as labeled numeric type | **10-20x** |
| `destring` | `cdestring` | Recast string as numeric type | **10-20x** |
| `gstats winsor` | `cwinsor` | Winsorization | **2-10x** |
| `bsample`/`sample` | `cbsample`/`csample` | Resampling with/without replacement | **2-5x** |
| `binscatter` | `cbinscatter` | Binned scatter plots | **10-40x** |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects | **10-20x** |
| `ivreghdfe` | `civreghdfe` | IV, GMM, etc. with multi-way fixed effects | **10-20x** |
| `qreg` | `cqreg` | Quantile regression | **2-4x** |

Some ctools programs have extended functionality. For instance, `cbinscatter` supports multi-way (high-dimensional) fixed effects and the alternative procedure to control for covariates described in [Cattaneo et al. (2024)](https://www.aeaweb.org/articles?id=10.1257/aer.20221576). Internal help files for each ctools program provide complete documentation.

Speedups will vary depending on a lot of factors. The regression commands (*hdfe) and binscatter will uniformly run faster. On the other hand, csort can be slower than sort if the number of variables in memory is sufficiently large relative to the number of observations due to the overhead involved with copying data between Stata memory and C memory (and back).



## Compatibility and Requirements

ctools is compatible with Stata 16.0+. It does not require any dependencies. Precompiled plugins for Windows, Mac OSX, and Linux are automatically included with installation, and ctools will automatically load the correct platform-specific plugin for your machine.


## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

### Manual Installation

1. Download or clone this repository
2. Copy the contents of the `build/` directory to a location on Stata's adopath (e.g., your personal ado directory)

The package automatically detects your operating system and architecture, loading the appropriate precompiled plugin.


## Building from Source Files


You probably do not need to compile the ctools C plugin yourself; GitHub automatically builds plugins for Windows, Mac (Intel and ARM/M-series), and Linux platforms. The installation instructions above will automatically include these plugins, and ctools chooses the appropriate plugin for your operating system/architecture automatically. If you want to build from source, you should make sure OpenMP is available on your system.

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
- ctools only works with datasets smaller than 2^31-1 observations (~2.147 billion obs). This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of the Stata function interface for C plugins and can only be addressed with an update to the Stata Plugin Interface. 


## Compatibility

- Stata 16.0+ for full compatibility; Stata 14.0+ probably works (untested).


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

Contributions are welcome. Please open an issue or submit a pull request.
