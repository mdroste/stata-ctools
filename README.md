```text
 ████╗█████╗ ████╗  ████╗ ██╗  █████╗
██╔══╝╚═██╔╝██╔═██╗██╔═██╗██║  ██╔══╝
██║     ██║ ██║ ██║██║ ██║██║  █████╗
██║     ██║ ██║ ██║██║ ██║██║  ╚══██║
╚████╗  ██║ ╚████╔╝╚████╔╝████╗█████║
 ╚═══╝  ╚═╝  ╚═══╝  ╚═══╝ ╚═══╝╚════╝
 some really fast Stata programs
```

[![Build](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml/badge.svg)](https://github.com/mdroste/stata-ctools/actions/workflows/build.yml)
![Version](https://img.shields.io/badge/version-1.2.3-blue)

**This is an initial release. Please report problems/suggestions on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** provides really fast implementations of common Stata routines. The syntax is largely unchanged from the commands they replace, sometimes with additional features; see the [compatibility table](docs/COMPATIBILITY.md).

### Stata built-in commands

| Stata command | Replaced with | Description | Approx. speedup¹ |
| --- | --- | --- | ---: |
| `import` | `cimport` | Import text, Excel, dBase, SAS, SPSS, and shapefiles | 40–60× |
| `export` | `cexport` | Export text, Excel, dBase, SAS, SPSS, and shapefiles | 50–70× |
| `sort` | `csort` | Sort dataset | 1.3–4× |
| `merge` | `cmerge` | Merge (join) datasets | 2–7× |
| `sample` | `csample` | Resampling without replacement | 5–9× |
| `bsample` | `cbsample` | Resampling with replacement | 3–5× |
| `encode` | `cencode` | Recast string as labeled numeric | 15–70× |
| `decode` | `cdecode` | Recast labeled numeric as string | 2–20× |
| `destring` | `cdestring` | Recast string as numeric type | 40–50× |
| `ipolate` | [`cinterpolate`](docs/USAGE.md#interpolation-cinterpolate) | Linear interpolation and extrapolation | 50–80× |
| `split` | [`csplit`](docs/USAGE.md#string-fields-csplit) | Split strings into fields | 70–130× |
| `qreg` | `cqreg` | Quantile regression | 4–6× |

### User-written commands

| Stata command | Replaced with | Description | Approx. speedup¹ |
| --- | --- | --- | ---: |
| `distinct` | [`cdistinct`](docs/commands/cdistinct.md) | Report distinct values | 10–30× |
| `hdfe` | [`chdfe`](docs/commands/chdfe.md) | Remove fixed effects from variables | — |
| `cipolate` | [`cinterpolate, cubic`](docs/USAGE.md#interpolation-cinterpolate) | Local cubic interpolation | — |
| `rangejoin` | [`crangejoin`](docs/USAGE.md#interval-joins-crangejoin) | Join observations within inclusive intervals | 10–20× |
| `winsor2`, `gstats winsor` | `cwinsor` | Winsorize variables | 2–5× vs. `gstats winsor` |
| `rangestat` | `crangestat` | Range statistics of variables | 80–110× |
| `psmatch2` | `cpsmatch` | Propensity score matching | 2–4× |
| `binscatter`, `binsreg` | `cbinscatter` | Binned scatter plots | 8–60× vs. `binscatter` |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects | 10–20× |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects | 10–20× |
| `ppmlhdfe` | `cppmlhdfe` | PPML with multi-way fixed effects | 5–6× |

¹ These are approximate speedups from a small number of benchmarks on data with **1M, 10M, and 25M observations** and 20 mixed-type variables. Values below 1× mean ctools was slower. Your machine's mileage may vary. See [benchmark details and results](docs/performance/README.md).

## Installation

Requires Stata 14.1 or later. In Stata:

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

For a drag-and-drop installation, download `ctools-<version>.zip` from
[GitHub Releases](https://github.com/mdroste/stata-ctools/releases), unzip it,
and copy its contents into a directory on your Stata
[ADOPATH](https://www.stata.com/support/faqs/programming/search-path-for-ado-files/)
(for example, PERSONAL, whose location is shown by `sysdir`). The ZIP includes
all ado and help files and all four platform plugins; Stata selects the right
plugin automatically. Restart Stata after replacing an installed plugin.

## Documentation

See the [documentation index](docs/README.md) for command guides, development and validation.

- [Command features](docs/FEATURES.md), [usage examples](docs/USAGE.md), and Stata's `help <command>`.
- [Compatibility and limitations](docs/COMPATIBILITY.md), including supported options and numerical differences.
- [Installation and building from source](docs/INSTALLATION.md) and the [developer guide](docs/development/README.md).

## Authorship

99.9% of this code was written by large language models. The early stuff was written in January 2026 by Claude Opus 4.5 and ChatGPT 5.2, back when it was still surprising for programs to be written like this. More recent bugfixes and improvements were written by Claude Fable 5.1, Claude Opus 5.5, and GPT 6.0/6.1 Astra/Sol. 


## Thanks

- Sergio Correia for [ftools](https://github.com/sergiocorreia/ftools), [reghdfe](https://github.com/sergiocorreia/reghdfe), and [ivreghdfe](https://github.com/sergiocorreia/ivreghdfe) (with Lars Vilhuber), as well as tons of helpful comments and suggestions
- Sergio Correia, Paulo Guimarães, and Thomas Zylkin for [ppmlhdfe](https://github.com/sergiocorreia/ppmlhdfe)
- Mauricio Caceres Bravo for [gtools](https://github.com/mcaceresb/gtools) and invaluable contributions to [cowsay](https://github.com/mdroste/stata-cowsay)
- Christopher (Kit) Baum, Mark E Schaffer, and Steven Stillman for [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html)
- Robert Picard, Nicholas J. Cox, and Roberto Ferrer for [rangestat](https://ideas.repec.org/c/boc/bocode/s458161.html)
- Nicholas J. Cox for [cipolate](https://ideas.repec.org/c/boc/bocode/s426202.html)
- Robert Picard for [rangejoin](https://ideas.repec.org/c/boc/bocode/s458162.html)
- Edwin Leuven and Barbara Sianesi for [psmatch2](https://ideas.repec.org/c/boc/bocode/s432001.html)
- Yujun Lian for [winsor2](https://ideas.repec.org/c/boc/bocode/s457765.html)
- Michael Stepner and Jessica Laird for [binscatter](https://michaelstepner.com/binscatter/)
- Matias D. Cattaneo, Richard K. Crump, Max H. Farrell, and Yingjie Feng for [binsreg](https://github.com/nppackages/binsreg)
- Michael Axtmann, Sascha Witt, Daniel Ferizovic, and Peter Sanders for the [In-place Parallel Super Scalar Samplesort (IPS⁴o)](https://github.com/SaschaWitt/ips4o) sorting algorithm.


## License

Project code is [MIT-licensed](LICENSE). The distribution includes [third-party notices](docs/licenses/THIRD_PARTY_NOTICES).


## Contributing

Contributions are welcome. Please open an [Issue](https://github.com/mdroste/stata-ctools/issues) or submit a pull request.
