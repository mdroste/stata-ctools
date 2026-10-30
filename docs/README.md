# Documentation

## Install and use ctools

- [Installation](INSTALLATION.md): install ctools or build it from source.
- [Usage](USAGE.md): command examples.
- [Features](FEATURES.md): additional features provided by ctools.
- [Compatibility](COMPATIBILITY.md): supported syntax and differences from the commands ctools replaces.
- [Import and export](IO_PARITY.md): supported file formats and Excel workbook behavior.
- [Tests of command options](OPTION_COVERAGE.md): which options are tested and where to find the tests.

In Stata, run `help <command>` for syntax and options. These guides give more
detail about individual commands:

| Data and file operations | Estimation and statistics |
| --- | --- |
| [csort](commands/csort.md) | [cbinscatter](commands/cbinscatter.md) |
| [cmerge](commands/cmerge.md) | [cqreg](commands/cqreg.md) |
| [cencode](commands/cencode.md) | [crangestat](commands/crangestat.md) |
| [cimport](commands/cimport.md) | [creghdfe](commands/creghdfe.md) |
| [cexport](commands/cexport.md) | [civreghdfe](commands/civreghdfe.md) |
| | [cppmlhdfe](commands/cppmlhdfe.md) |

## Build, test and contribute

- [Developer guide](development/README.md): how the code works and how to add commands.
- [Build requirements](development/PLATFORMS.md): compilers, OpenMP and required libraries for each platform.
- [Code style](development/STYLE.md): formatting conventions.
- [Command implementations](development/COMMAND_INVENTORY.md): generated list of commands, C functions and cleanup functions.
- [Testing guide](../validation/README.md): tests, test data and benchmarks.
- [Generated import tables](../tools/generation/README.md): check or regenerate character encoding and locale tables.

## Benchmarks

- [Speed comparisons in the main README](performance/README.md): how they are measured and how to repeat them.
- [HDFE comparisons](performance/HDFE_DATASETS.md): compare fixed-effects regression on example datasets.

Scripts write the published tables and measurements to `benchmarks/`. Save
investigation reports in `reports/` and temporary files in `temp/` at the
repository root; both directories are ignored by Git.
