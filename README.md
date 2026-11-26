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
![Version](https://img.shields.io/badge/version-1.3.2-blue)

**This is an initial release. Please report problems/suggestions on [Issues](https://github.com/mdroste/stata-ctools/issues).**

## Overview

**ctools** provides really fast drop-in replacements for common Stata programs. The syntax of ctools is largely unchanged from the commands they replace, sometimes with additional features; see the [compatibility table](docs/COMPATIBILITY.md).

Each ctools command reimplements an existing Stata command and is pinned to one version of it. [Attribution](docs/ATTRIBUTION.md) lists each original command, its authors, and the pinned version. The same information appears in each command's help file.

### Stata built-in commands

| Stata command | Replaced with | Description | Approx. speedup¹ |
| --- | --- | --- | ---: |
| `import` | `cimport` | Import text, Excel, dBase, SAS, SPSS, and shapefiles | 10–20× |
| `export` | `cexport` | Export text, Excel, dBase, SAS, SPSS, and shapefiles | 50–80× |
| `sort` | `csort` | Sort dataset | 1.7–4×² |
| `merge` | `cmerge` | Merge (join) datasets | 2–6×² |
| `joinby` | [`cjoinby`](docs/commands/cjoinby.md) | Form all pairwise combinations within groups | 20–70× |
| `reshape` | [`creshape`](docs/commands/creshape.md) | Convert data from wide to long form and vice versa | 20–40×⁷ |
| `collapse` | [`ccollapse`](docs/commands/ccollapse.md) | Make dataset of summary statistics | 1.1–20×⁴ |
| `egen` | [`cegen`](docs/commands/cegen.md) | Extensions to generate | 3–60×⁵ |
| `sample` | `csample` | Resampling without replacement | 6–10× |
| `bsample` | `cbsample` | Resampling with replacement | 2–4× |
| `encode` | `cencode` | Recast string as labeled numeric | 30–160× |
| `decode` | `cdecode` | Recast labeled numeric as string | 2–5× |
| `destring` | `cdestring` | Recast string as numeric type | 20–40× |
| `ipolate` | [`cinterpolate`](docs/USAGE.md#interpolation-cinterpolate) | Linear interpolation and extrapolation | 50–90× |
| `split` | [`csplit`](docs/USAGE.md#string-fields-csplit) | Split strings into fields | 140–270× |
| `qreg` | `cqreg` | Quantile regression | 8–10×³ |

### User-written commands

See [Attribution](docs/ATTRIBUTION.md#community-contributed-commands) for the authors of these commands and the versions ctools is pinned to.

| Stata command | Replaced with | Description | Approx. speedup¹ |
| --- | --- | --- | ---: |
| `distinct` | [`cdistinct`](docs/commands/cdistinct.md) | Report distinct values | 10–40× |
| `cipolate` | [`cinterpolate, cubic`](docs/USAGE.md#interpolation-cinterpolate) | Local cubic interpolation | 40–80× |
| `rangejoin` | [`crangejoin`](docs/USAGE.md#interval-joins-crangejoin) | Join observations within inclusive intervals | 10–30× |
| `winsor2`<br>`gstats winsor` | `cwinsor` | Winsorize variables | 20–60×<br>1.9–4× |
| `rangestat` | `crangestat` | Range statistics of variables | 30–70× |
| `binscatter`<br>`binscatter2`<br>`binsreg` | `cbinscatter` | Binned scatter plots | 7–40×<br>3–20×<br>3–30× |
| `psmatch2` | `cpsmatch` | Propensity score matching | 3–4× |
| `reghdfe`<br>`fereg` | `creghdfe` | OLS with multi-way fixed effects | 8–20×<br>0.3–0.5× |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects | 10–20× |
| `ppmlhdfe`<br>`feppml` | `cppmlhdfe` | PPML with multi-way fixed effects | 6–9×³<br>0.3–0.5× |
| `hdfe` | [`chdfe`](docs/commands/chdfe.md) | Remove fixed effects from variables | 5–20× |

¹ These are approximate speedups from a small number of benchmarks on data with **1M, 10M, and 25M observations** and 20 mixed-type variables. Where a ctools command replaces several commands, each has its own speedup, in the same order. Values below 1× mean ctools was slower. Your machine's mileage may vary. See [benchmark details and results](docs/performance/README.md).

² `csort` and `cmerge` work better on data with many observations relative to variables; their speedups shrink as the number of variables grows.

³ Measured at 1M and 10M observations only: `qreg` and `ppmlhdfe` take minutes per run at 25M.

⁴ Measured at 10M observations only, on nine combinations of statistics and `by()` variables, and not yet part of the benchmark profile above; the smallest speedup is without `by()`. Against `gcollapse` the range was 4–20×; see the [command page](docs/commands/ccollapse.md#performance).

⁵ Measured at 2M and 10M observations on 33 combinations of `egen` functions and `by()` variables, and not yet part of the benchmark profile above. The smallest speedups are for sums without `by()` and across variables, where `egen` is already fast; the largest are for `mode()`, `tag()` and `seq()` within groups. Against `gegen` the range at 2M was 1.3–6× for the functions `gegen` implements; see the [command page](docs/commands/cegen.md#performance).

⁷ Measured on data of 1M and 10M observations in long form (ten values of `j`, three reshaped variables and three others), in both directions, against `reshape`'s default `favor(speed)`, and not yet part of the benchmark profile above. See the [`creshape` guide](docs/commands/creshape.md#performance).

## Installation

ctools requires Stata 14.1+, macOS 11 or later (Apple Silicon or Intel), Windows 10 or
later (x86-64), or x86-64 Linux with glibc 2.17 or later. See [PLATFORMS.md](docs/development/PLATFORMS.md).

You can install the latest version of ctools from this GitHub repository within Stata:
```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

Once ctools is installed, `ctools, update` checks GitHub for the latest release
and installs it if it is newer than your version.

Alternatively, you can install a recent stable version from SSC:
```stata
ssc install ctools, replace
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
- [Attribution](docs/ATTRIBUTION.md): the original command behind each ctools command, its authors, and the pinned version.
- [Installation and building from source](docs/INSTALLATION.md) and the [developer guide](docs/development/README.md).

## Authorship

99.9% of this code was written by AI. Most of the package was written in January 2026 by Claude Opus 4.5 and ChatGPT 5.2, back when it was still surprising for programs to be written like this. More recent bugfixes and improvements were written by Claude Fable 5.1, Claude Opus 5.5, and GPT 6.0/6.1 Astra/Sol.


## Thanks

ctools reimplements commands written by the people below. [Attribution](docs/ATTRIBUTION.md) lists the version of each command that ctools is pinned to, along with references to cite.

- Sergio Correia for [ftools](https://github.com/sergiocorreia/ftools), [reghdfe](https://github.com/sergiocorreia/reghdfe) (with Noah Constantine), [ivreghdfe](https://github.com/sergiocorreia/ivreghdfe) (with Lars Vilhuber), and [hdfe](https://ideas.repec.org/c/boc/bocode/s457985.html), as well as tons of helpful comments and suggestions
- Sergio Correia, Paulo Guimarães, and Thomas Zylkin for [ppmlhdfe](https://github.com/sergiocorreia/ppmlhdfe)
- Mauricio Caceres Bravo for [gtools](https://github.com/mcaceresb/gtools) and invaluable contributions to [cowsay](https://github.com/mdroste/stata-cowsay)
- Christopher (Kit) Baum, Mark E Schaffer, and Steven Stillman for [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html)
- Frank Kleibergen, Mark E Schaffer, and Frank Windmeijer for [ranktest](https://ideas.repec.org/c/boc/bocode/s456865.html)
- David Roodman for [boottest](https://ideas.repec.org/c/boc/bocode/s458121.html)
- Gary Longton and Nicholas J. Cox for [distinct](https://ideas.repec.org/c/boc/bocode/s424201.html)
- Robert Picard, Nicholas J. Cox, and Roberto Ferrer for [rangestat](https://ideas.repec.org/c/boc/bocode/s458161.html)
- Nicholas J. Cox for [cipolate](https://ideas.repec.org/c/boc/bocode/s426202.html), and for Stata's official `split` and `destring`, which the Stata manual credits to him
- Robert Picard for [rangejoin](https://ideas.repec.org/c/boc/bocode/s458162.html)
- Edwin Leuven and Barbara Sianesi for [psmatch2](https://ideas.repec.org/c/boc/bocode/s432001.html)
- Yujun Lian for [winsor2](https://ideas.repec.org/c/boc/bocode/s457765.html)
- Michael Stepner and Jessica Laird for [binscatter](https://michaelstepner.com/binscatter/)
- Matias D. Cattaneo, Richard K. Crump, Max H. Farrell, and Yingjie Feng for [binsreg](https://github.com/nppackages/binsreg)
- Michael Axtmann, Sascha Witt, Daniel Ferizovic, and Peter Sanders for the [In-place Parallel Super Scalar Samplesort (IPS⁴o)](https://github.com/SaschaWitt/ips4o) sorting algorithm.


## License

Project code is [MIT-licensed](LICENSE). The distribution includes [third-party notices](docs/licenses/THIRD_PARTY_NOTICES) for the libraries compiled into the plugins; after `net install`, `help ctools_licenses` shows both texts.


## Contributing

Contributions are welcome. Please open an [Issue](https://github.com/mdroste/stata-ctools/issues) or submit a pull request.
