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

**ctools** provides really fast implementations of common Stata routines. The syntax is largely unchanged from the commands they replace, sometimes with additional features; see the [compatibility table](docs/COMPATIBILITY.md).

| Stata command | Replaced with | Description | Approx. speedup¹ |
| --- | --- | --- | ---: |
| `import` | `cimport` | Import text, Excel, dBase, SAS, SPSS, and shapefiles | 30–50× |
| `export` | `cexport` | Export text, Excel, dBase, SAS, SPSS, and shapefiles | 40–60× |
| `sort` | `csort` | Sort dataset | 1–3× |
| `merge` | `cmerge` | Merge (join) datasets | 1–7× |
| `sample` | `csample` | Resampling without replacement | 5–8× |
| `bsample` | `cbsample` | Resampling with replacement | 3–4× |
| `encode` | `cencode` | Recast string as labeled numeric | 15–60× |
| `decode` | `cdecode` | Recast labeled numeric as string | 1–10× |
| `destring` | `cdestring` | Recast string as numeric type | 40× |
| `ipolate` | [`cipolate`](docs/USAGE.md#linear-interpolation-cipolate) | Linear interpolation and extrapolation | 10–12× |
| `split` | [`csplit`](docs/USAGE.md#string-fields-csplit) | Split strings into fields | 15–18× |
| `rangejoin` | [`crangejoin`](docs/USAGE.md#interval-joins-crangejoin) | Join observations within inclusive intervals | 10-15× |
| `qreg` | `cqreg` | Quantile regression | 5–6× |
| `gstats winsor` | `cwinsor` | Winsorize variables | 1.5–2× |
| `rangestat` | `crangestat` | Range statistics of variables | 80–110× |
| `psmatch2` | `cpsmatch` | Propensity score matching | 1.5–2.5× |
| `binscatter` | `cbinscatter` | Binned scatter plots | 5–10× |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects | 9.0–13.0× |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects | 14–16× |
| `ppmlhdfe` | `cpplmhdfe` | PPML with multi-way fixed effects | 5–6× |

¹ Approximate speedups from three benchmarks with **1M, 10M, and 25M observations** and 20 mixed-type variables. Import/export timings use CSV. Values below 1× mean ctools was slower. Your mileage may vary. See [benchmark details and results](docs/README_BENCHMARK.md).

## Installation

Requires Stata 14.1 or later. In Stata:

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

## Documentation

- [Command features](FEATURES.md), [usage examples](docs/USAGE.md), and Stata's `help <command>`.
- [Compatibility and limitations](docs/COMPATIBILITY.md), including supported options and numerical differences.
- [Installation and building from source](docs/INSTALLATION.md) and the [developer guide](DEVELOPERS.md).

## Authorship

99.9% of this code was written by LLMs. The early stuff was written in January 2026 by Claude Opus 4.5 and ChatGPT 5.2; recent bugfixes and improvements were written by Claude Fable 5.1, Claude Opus 5.5, and GPT 6.0 Astra/Sol. 


## Thanks

- [Sergio Correia](https://github.com/sergiocorreia) for [ftools](https://github.com/sergiocorreia/ftools), [reghdfe](https://github.com/sergiocorreia/reghdfe), and [ivreghdfe](https://github.com/sergiocorreia/ivreghdfe) (with [Lars Vilhuber](https://www.ilr.cornell.edu/people/lars-vilhuber)), as well as tons of helpful comments and suggestions
- [Mauricio Caceres Bravo](https://mcaceresb.github.io/) for [gtools](https://github.com/mcaceresb/gtools) and invaluable contributions to [cowsay](https://github.com/mdroste/stata-cowsay)
- [Christopher (Kit) Baum](https://www.bc.edu/bc-web/schools/morrissey/departments/economics/people/faculty-directory/christopher-baum.html), [Mark E Schaffer](https://www.hw.ac.uk/profiles/uk/school/ebs/faculty/mark-schaffer), and [Steven Stillman](https://www.unibz.it/en/faculties/economics-management/academic-staff/person/36390-steven-stillman) for [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html)
- [Robert Picard](https://ideas.repec.org/f/ppi320.html), [Nicholas J. Cox](https://ideas.repec.org/e/pco34.html), and Roberto Ferrer for [rangestat](https://ideas.repec.org/c/boc/bocode/s458161.html)
- [Sascha Witt](https://github.com/SaschaWitt) for the [In-place Parallel Super Scalar Samplesort (IPS⁴o)](https://github.com/SaschaWitt/ips4o) sorting algorithm.
- Claude Code


## License

Project code is [MIT-licensed](LICENSE). The distribution includes [third-party notices](THIRD_PARTY_NOTICES).


## Contributing

Contributions are welcome. Please open an [Issue](https://github.com/mdroste/stata-ctools/issues) or submit a pull request.
