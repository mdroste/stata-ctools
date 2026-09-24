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
![Version](https://img.shields.io/badge/version-1.0.2-blue)

**This is an initial release. Please report problems/suggestions on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** provides C implementations of common Stata data operations and estimators. Supported syntax and statistical behavior vary by command; see the [compatibility table](docs/COMPATIBILITY.md).

| Stata command | Replaced with | Description |
| --- | --- | --- |
| `import` | `cimport` | Import text-delimited and Excel data |
| `export` | `cexport` | Export text-delimited and Excel data |
| `sort` | `csort` | Sort dataset |
| `merge` | `cmerge` | Merge (join) datasets |
| `sample` | `csample` | Resampling without replacement |
| `bsample` | `cbsample` | Resampling with replacement |
| `encode` | `cencode` | Recast string as labeled numeric |
| `decode` | `cdecode` | Recast labeled numeric as string |
| `destring` | `cdestring` | Recast string as numeric type |
| `qreg` | `cqreg` | Quantile regression |
| `gstats winsor` | `cwinsor` | Winsorize variables |
| `rangestat` | `crangestat` | Range statistics of variables |
| `psmatch2` | `cpsmatch` | Propensity score matching |
| `binscatter` | `cbinscatter` | Binned scatter plots |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects |
| `ppmlhdfe` | `cpplmhdfe` | PPML with multi-way fixed effects |

Some ctools programs have extended functionality. For instance, `cbinscatter` supports multi-way fixed effects and the procedure to control for covariates characterized by [Cattaneo et al. (2024)](https://www.aeaweb.org/articles?id=10.1257/aer.20221576). See [FEATURES.md](FEATURES.md) for a brief description of the new features implemented for each command above. Each command also has an associated internal help file (e.g. `help cbinscatter`).

Performance depends on dataset shape, options, and hardware. Data transfer overhead can make sorting or merging wide datasets slower than native Stata. See the [benchmark reporting requirements](docs/COMPATIBILITY.md#performance-evidence) before interpreting speed comparisons.


## Compatibility and Requirements

ctools is compatible with Stata 14.1+ (plugin interface version 3.0). The [platform contract](docs/PLATFORMS.md) specifies CPU/OS baselines and runtime dependencies, including Linux libgomp.

**Note:** ctools does not support datasets exceeding 2^31 (~2.147 billion) observations. This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of Stata's plugin API for C and can only be addressed with an internal Stata update. ctools will gracefully exit with an error if your dataset exceeds this limit.


## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

### Manual Installation

1. Download a complete release archive and validate its extracted directory with `python3 validation/check_package.py /path/to/release/build`.
2. In Stata, run `net install ctools, from("/absolute/path/to/release/build") replace`.

A source checkout's `build/` is not necessarily an installable distribution: it
may contain only the locally compiled binary. After building from source, use
`make package BUILD_REVISION=<revision> PACKAGE_DIR=dist/ctools-<revision>` to
produce a complete package for the host platform. Install from that staged
directory. Its manifest explicitly identifies the platform and `BUILD_INFO.json`
records the revision and file checksums. Full release packages contain all four
platform plugins; publication requires the licensed correctness gate.


## Building from Source Files

Release identity is defined in `release.json`. Run `python3 scripts/sync_release.py` after changing it; CI checks the generated copies. `ctools, version` reports the installed ado and plugin identities.


You probably do not need to compile the ctools plugin yourself; GitHub automatically builds plugins for Windows, Mac, and Linux. If you want or need to build plugins from source, follow the [platform build instructions](docs/PLATFORMS.md), including the deployment-compatible macOS OpenMP runtime.

```bash
make              # Build for current platform
make all          # Build all platform plugins (requires every toolchain)
make package      # Build and stage an installable host-platform package
make check        # Check build dependencies
make clean        # Remove compiled files
```

If you have a workstation/server CPU with lots of cache, you might want to try playing around with the settings in [src/ctools_config.h](src/ctools_config.h).

See [DEVELOPERS.md](DEVELOPERS.md) for additional information on ctools' architecture and core logic.


## Usage Notes

- All ctools programs follow a basic structure: (1) copy data from Stata to C; (2) operate on that data (3) return data from C to Stata. *This means that ctools programs require more memory than the programs they replace.* In addition, some commands will run faster when they involve fewer variables, or when you have fewer variables in memory. For instance,  `csort`'s runtime is heavily dependent on the number of variables in memory, and can be slower Stata's built-in `sort` if the dataset is relatively wide (e.g. 100+ variables) due to this data transfer overhead.
- The default options for `csort` and `cmerge` require a lot of memory. For `csort`, you probably require 2-3 times as much memory as your dataset. For `csort`, you can reduce this memory overhead with the optional argument `stream(#)`, which reads variables a handful at a time rather than all at once (at a modest cost to runtime).

## Issues
- [ ] Precision of accumulated scalar statistics (e.g. total/model/residual sum of squares) associated with regression output (cqreg, creghdfe, civredhfe): only matches replacement to ~7 significant digits.

## Authorship

99.9% of the code in this repository was written by Claude Opus 4.5 (through Claude Code), with some debugging and refactoring assistance from OpenAI GPT 5.2 (through Codex).


## Thanks

- [Sergio Correia](https://github.com/sergiocorreia) for [ftools](https://github.com/sergiocorreia/ftools), [reghdfe](https://github.com/sergiocorreia/reghdfe), and [ivreghdfe](https://github.com/sergiocorreia/ivreghdfe) (with [Lars Vilhuber](https://www.ilr.cornell.edu/people/lars-vilhuber))
- [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for [gtools](https://github.com/mcaceresb/gtools) and invaluable contributions to [cowsay](https://github.com/mdroste/stata-cowsay)
- [Christopher (Kit) Baum](https://www.bc.edu/bc-web/schools/morrissey/departments/economics/people/faculty-directory/christopher-baum.html), [Mark E Schaffer](https://www.hw.ac.uk/profiles/uk/school/ebs/faculty/mark-schaffer), and [Steven Stillman](https://www.unibz.it/en/faculties/economics-management/academic-staff/person/36390-steven-stillman) for [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html)
- [Robert Picard](https://ideas.repec.org/f/ppi320.html), [Nicholas J. Cox](https://ideas.repec.org/e/pco34.html), and Roberto Ferrer for [rangestat](https://ideas.repec.org/c/boc/bocode/s458161.html)
- [Sascha Witt](https://github.com/SaschaWitt) for the [In-place Parallel Super Scalar Samplesort (IPS⁴o)](https://github.com/SaschaWitt/ips4o) sorting algorithm.
- Claude Code


## License

Project code is [MIT-licensed](LICENSE). The distribution includes [third-party notices](THIRD_PARTY_NOTICES).


## Contributing

Contributions are welcome. Please open an [Issue](https://github.com/mdroste/stata-ctools/issues) or submit a pull request.
