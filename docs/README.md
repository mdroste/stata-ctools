# Documentation

## Install and use ctools

- [Installation](INSTALLATION.md): install ctools or build it from source.
- [Usage](USAGE.md): command examples.
- [Interpolation](USAGE.md#interpolation-cinterpolate), [string splitting](USAGE.md#string-fields-csplit), and [interval joins](USAGE.md#interval-joins-crangejoin).
- [Features](FEATURES.md): additional features provided by ctools.
- [Compatibility](COMPATIBILITY.md): supported syntax and differences from the commands ctools replaces.
- [Import and export](IO_PARITY.md): supported file formats and Excel workbook behavior.
- [Check command options](OPTION_COVERAGE.md).

In Stata, run `help <command>` for syntax and options. These guides give more
detail about individual commands:

| Data and file operations | Estimation and statistics |
| --- | --- |
| [csort](commands/csort.md) | [cbinscatter](commands/cbinscatter.md) |
| [cmerge](commands/cmerge.md) | [cqreg](commands/cqreg.md) |
| [cencode](commands/cencode.md) | [crangestat](commands/crangestat.md) |
| [cimport](commands/cimport.md) | [creghdfe](commands/creghdfe.md) |
| [cexport](commands/cexport.md) | [civreghdfe](commands/civreghdfe.md) |
| [cdistinct](commands/cdistinct.md) | |
| [chdfe](commands/chdfe.md) | [cppmlhdfe](commands/cppmlhdfe.md) |

## Build, test and contribute

- [Developer guide](development/README.md): how the code works and how to add commands.
- [Build requirements](development/PLATFORMS.md): compilers, OpenMP and required libraries for each platform.
- [Code style](development/STYLE.md): formatting conventions.
- [Command implementations](development/COMMAND_INVENTORY.md): generated list of commands, C functions and cleanup functions.
- [Run the tests](../validation/README.md).
- [Generated import tables](../tools/generation/README.md): check or regenerate character encoding and locale tables.

## Benchmarks

- [Run timing comparisons](performance/README.md).
- [HDFE comparisons](performance/HDFE_DATASETS.md): compare fixed-effects regression on example datasets.
