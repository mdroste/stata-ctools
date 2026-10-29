# Documentation

## Install and use ctools

- [Installation](INSTALLATION.md): Stata installation, release archives and source builds.
- [Usage](USAGE.md): examples for the command families.
- [Features](FEATURES.md): capabilities added by ctools.
- [Compatibility](COMPATIBILITY.md): supported syntax and differences from reference commands.
- [Import/export parity](IO_PARITY.md): file formats and workbook behavior.
- [Option coverage](OPTION_COVERAGE.md): generated inventory and validation evidence.

Stata's `help <command>` contains the installed command reference. The following
command guides provide additional implementation and usage detail:

| Data and file operations | Estimation and statistics |
| --- | --- |
| [csort](commands/csort.md) | [cbinscatter](commands/cbinscatter.md) |
| [cmerge](commands/cmerge.md) | [cqreg](commands/cqreg.md) |
| [cencode](commands/cencode.md) | [crangestat](commands/crangestat.md) |
| [cimport](commands/cimport.md) | [creghdfe](commands/creghdfe.md) |
| [cexport](commands/cexport.md) | [civreghdfe](commands/civreghdfe.md) |
| | [cppmlhdfe](commands/cppmlhdfe.md) |

## Develop and validate

- [Developer guide](development/README.md): architecture, numerical contracts and adding commands.
- [Platform contracts](development/PLATFORMS.md): compilers, OpenMP and binary dependencies.
- [Code style](development/STYLE.md): formatting scope and conventions.
- [Command inventory](development/COMMAND_INVENTORY.md): generated plugin dispatch inventory.
- [Validation guide](../validation/README.md): correctness gates, fixtures and performance tools.
- [Offline generators](../tools/generation/README.md): pinned inputs and generated import tables.

## Performance evidence

- [README benchmark methodology](performance/README.md).
- [HDFE dataset comparisons](performance/HDFE_DATASETS.md).

Reproducible tables and measurements live at fixed paths in `benchmarks/`.
Investigation reports and disposable output belong in the ignored `reports/`
and `temp/` directories at the repository root.
