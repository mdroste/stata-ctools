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

| Stata command | Replaced with | Description | Approx. speedup¹ |
| --- | --- | --- | ---: |
| `import` | `cimport` | Import text, Excel, dBase, SAS, SPSS, and shapefiles | — |
| `export` | `cexport` | Export text, Excel, dBase, SAS, SPSS, and shapefiles | — |
| `sort` | `csort` | Sort dataset | — |
| `merge` | `cmerge` | Merge (join) datasets | — |
| `sample` | `csample` | Resampling without replacement | — |
| `bsample` | `cbsample` | Resampling with replacement | — |
| `encode` | `cencode` | Recast string as labeled numeric | — |
| `decode` | `cdecode` | Recast labeled numeric as string | — |
| `destring` | `cdestring` | Recast string as numeric type | — |
| `ipolate` | [`cipolate`](#linear-interpolation-cipolate) | Linear interpolation and extrapolation | 16–19× |
| `split` | [`csplit`](#string-fields-csplit) | Split strings into fields | 20–23× |
| `rangejoin` | [`crangejoin`](#interval-joins-crangejoin) | Join observations within inclusive intervals | 15–20× |
| `qreg` | `cqreg` | Quantile regression | — |
| `gstats winsor` | `cwinsor` | Winsorize variables | — |
| `rangestat` | `crangestat` | Range statistics of variables | — |
| `psmatch2` | `cpsmatch` | Propensity score matching | — |
| `binscatter` | `cbinscatter` | Binned scatter plots | — |
| `reghdfe` | `creghdfe` | OLS with multi-way fixed effects | — |
| `ivreghdfe` | `civreghdfe` | 2SLS/GMM with multi-way fixed effects | — |
| `ppmlhdfe` | `cpplmhdfe` | PPML with multi-way fixed effects | — |

¹ Measured against the command in the first column on plausible synthetic datasets with **1M, 5M, and 20M observations**, using Stata/MP 18.5 with 12 processors and `threads(12)` on an Apple M4 Pro. Figures are rounded ranges of median speedups across those sizes. Join sizes refer to **each input**; output reaches 39.96M rows. All nine full-output comparisons matched exactly. See the [benchmark report, timings, and reproducible code](docs/BENCHMARKS_README_COMMANDS.md) for dataset details and limitations. A dash means that command was not benchmarked in this comparison.

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

### Linear interpolation: `cipolate`

`cipolate` implements Stata's **linear** `ipolate`, including grouped interpolation and optional extrapolation. It is distinct from the SSC cubic-interpolation command also named `cipolate`; check `which cipolate` if both packages are installed.

```stata
cipolate yvar xvar [if] [in], generate(newvar) ///
    [by(varlist) epolate threads(#) verbose]
```

- `generate(newvar)` is required and creates a double-precision result. The source variables and observation order are preserved.
- `by(varlist)` interpolates separately within numeric or string groups. The `by:` prefix also works; use either the prefix or `by()`. Missing group values participate, with `.`, `.a`, etc. treated as distinct groups. Group strings must fit within 2,045 bytes.
- At duplicate `xvar` values, nonmissing `yvar` values are averaged. Missing values between known points are filled by linear interpolation. Missing `xvar` values and observations excluded by `if`/`in` receive missing results.
- Without `epolate`, results outside the known range are missing. `epolate` extends the line through the nearest two distinct known x points at each end; groups with fewer than two known points cannot extrapolate.
- `threads(#)` accepts a nonnegative integer and caps C worker threads; `threads(0)`, the default, uses the automatic setting. `verbose` reports phase timings.

```stata
* Fill missing monthly sales separately for each firm.
cipolate sales month, generate(sales_i) by(firm)

* Also extrapolate at the ends of each firm's observed series.
cipolate sales month, generate(sales_e) by(firm) epolate threads(8)
```

The output variable is the result; no command-specific stored results are promised. See `help cipolate` or the [full help file](build/cipolate.sthlp).

### String fields: `csplit`

`csplit` implements `split` for fixed-width strings, including multiple delimiters, field limits, and optional numeric conversion. It preserves the source variable and observation order.

```stata
csplit strvar [if] [in] [, generate(stub) parse(parse_strings) notrim ///
    limit(#) destring force float ignore(strings) percent threads(#) verbose]
```

- `generate(stub)` names outputs `stub1`, `stub2`, etc.; the default stub is the source variable's name. Output names must be new. The number of fields is the maximum needed across selected observations; unselected observations have empty string outputs (missing after numeric conversion).
- With no `parse()`, spaces separate fields and successive remainders are trimmed. `parse("|" ";" "::")` supplies literal delimiter strings, which may contain multiple characters. If delimiters match at the same position, the last listed one takes precedence. This is literal splitting, without CSV-style quoted-field handling.
- By default, outer spaces are trimmed before parsing. With explicit delimiters, spaces within the remaining text are preserved. `notrim` retains outer spaces too and requires explicit `parse()`. Adjacent delimiters preserve empty fields; a trailing delimiter does not create another final field.
- `limit(#)` must be positive and retains at most that many fields; text after those fields is discarded.
- `destring` applies Stata's native numeric conversion to the fields. Its options `force` (unconvertible values become missing), `float` (float storage), `ignore(strings)` (remove specified characters), and `percent` (convert percentage values) require `destring`. Without `force`, a field containing nonnumeric values can remain a string, as with `destring`.
- `threads(#)` accepts a nonnegative integer and caps C workers, with zero selecting the default; `verbose` reports phase timings. The C path supports fixed strings up to 2,045 bytes. `strL` input falls back to native `split`, so it does not receive the C speedup. String bytes are preserved without Unicode normalization.

```stata
* Example code: NY|retail|SKU000042|12345678|USD|2024-01-31
csplit code, generate(part) parse("|")

* Convert comma-separated numeric fields, ignoring dollar signs.
csplit amounts, generate(amount) parse(",") destring ignore("$")
```

The returned macros `r(varlist)` and `r(nvars)` contain the output names and count; scalar `r(k_new)` also contains the count. An empty selected sample returns error 2000. See `help csplit` or the [full help file](build/csplit.sthlp).

### Interval joins: `crangejoin`

`crangejoin` implements the matching interface of Robert Picard's SSC `rangejoin` 1.1.3. Each observation in memory is matched to every using observation whose numeric key falls within its **inclusive** interval and whose `by()` values agree. The result replaces the dataset in memory. Neither `rangejoin` nor `rangestat` needs to be installed to use it.

```stata
crangejoin keyvar low high using filename [, by(varlist) keepusing(varlist) ///
    prefix(string) suffix(string) all threads(#) verbose]
```

- `keyvar` must be numeric in the using dataset. Numeric master variables supplied as `low` and `high` give **absolute bounds**. Literal numbers are offsets from the same-named master key when it exists, and absolute bounds otherwise. For example, `date -1 1` matches dates within one day of the master date.
- A missing bound (`.` or a missing value in a bound variable) leaves that side unbounded. Missing using keys never match. If the master contains `keyvar`, missing master keys never match either, even with explicit bound variables. Reversed intervals do not match. If there are no valid master intervals or no nonmissing using keys, the command returns error 2000.
- `by(varlist)` restricts matches to equal numeric or string groups. Missing group values participate, and extended numeric missings remain distinct. The key cannot also appear in `by()`.
- `keepusing(varlist)` selects using payload variables; key and group variables are included automatically. Using variables that share master names receive `prefix()` and/or `suffix()`; the default is `suffix(_U)` when neither is specified. `all` applies renaming to every using variable except the group variables. Conflicting output names return an error.
- Unmatched master rows remain, with missing using values. Master observation order is preserved; matches within each master row are ordered by using key, then by original using order for ties. The wrapper preserves variable metadata and restores the master dataset on failure.
- `threads(#)` accepts a nonnegative integer and caps C workers, with zero selecting the default; `verbose` reports phase timings. Numeric and fixed-width strings are supported. Master variables and selected using variables cannot be `strL`; recast them to fixed strings first only when lossless. Inputs and output must each have at most 2,147,483,647 rows. Memory use and runtime grow with the number of matches, which can greatly exceed input size.

```stata
* Match firm offers within one day of each inquiry's date.
crangejoin date -1 1 using offers.dta, by(firm) keepusing(price venue)

* Match events to observation-specific absolute date windows.
crangejoin date window_start window_end using events.dta, by(firm)

* Attach all history within each person; watch the resulting match count.
crangejoin date . . using history.dta, by(person) suffix(_history)
```

There is no `if`/`in` syntax: filter the master dataset beforehand. The joined dataset is the result; no command-specific stored results are promised. See `help crangejoin` or the [full help file](build/crangejoin.sthlp).

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
