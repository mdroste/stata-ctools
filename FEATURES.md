# ctools New Features

This document describes features available **only** in ctools commands—functionality that does not exist in the standard Stata commands or popular user-written alternatives they replace.

Most ctools commands use C implementations with parallel processing. `cdecode`
uses Stata's native decoding engine to preserve literal and long value labels.
This document describes added capabilities; see `docs/COMPATIBILITY.md` for the
supported scope and differences from reference commands.

---

## Command overview

| Command | Replaces | New features? |
|---------|----------|---------------|
| `csort` | `sort` | Yes |
| `cmerge` | `merge` | — |
| `cimport` | `import delimited`, `import excel` | Yes |
| `cexport` | `export delimited`, `export excel` | Yes |
| `cencode` | `encode` | Yes |
| `cdecode` | `decode` | Yes |
| `cdestring` | `destring` | — |
| `cwinsor` | `winsor2` | — |
| `csample` | `sample` | Yes |
| `cbsample` | `bsample` | — |
| `crangestat` | `rangestat` | — |
| `creghdfe` | `reghdfe` | — |
| `civreghdfe` | `ivreghdfe` | — |
| `cqreg` | `qreg` | Yes |
| `cpsmatch` | `psmatch2` | — |
| `cbinscatter` | `binscatter` | Yes |

---

## cencode / cdecode

Stata's native `encode` and `decode` commands only process one variable at a time and require creating new variables.

**ctools adds:**

- **Varlist support**: Encode or decode multiple variables in a single command
- **Replace option**: Modify variables in-place instead of creating new ones

```stata
* Encode three variables at once, replacing in-place
cencode var1 var2 var3, replace
```

---

## csort

Stata's native `sort` command provides no user control over memory usage or algorithm selection.

**ctools adds:**

- **Algorithm selection**: `algorithm()` option to choose from `auto`, `counting`, `lsd`, `msd`, `ips4o`, `timsort`, `sample`, or `merge`
- **Streaming mode**: `stream()` option reduces memory usage for wide datasets by permuting non-key variables sequentially

```stata
* Use counting sort for integer data
csort year, algorithm(counting)

* Streaming mode for memory-constrained sorting
csort id, stream(4)
```

---

## cimport / cexport

Stata's `import delimited` and `import excel` (and their `export` counterparts) are separate commands with different syntax. ctools unifies each pair under a single command.

**ctools adds:**

- **Unified interface**: `cimport delimited` and `cimport excel` under one command (same for `cexport`), providing a consistent syntax for both CSV and Excel workflows
- **Locale-aware number parsing** (cimport): `locale()` and `parselocale` options for automatic handling of locale-specific decimal and grouping separators (e.g., `locale(de_DE)` for German number formatting)

```stata
* Import a CSV file
cimport delimited using data.csv, clear

* Import an Excel file with the same command
cimport excel using data.xlsx, clear firstrow

* Import with German number formatting (comma = decimal, period = thousands)
cimport delimited using data_de.csv, clear locale(de_DE) parselocale

* Export to CSV and Excel with the same command family
cexport delimited using output.csv, replace
cexport excel using output.xlsx, replace
```

---

## csample

Both native `sample` and `csample` support `by()` for within-group sampling.

**ctools adds:**

- **Parallel implementation**: sample separately within groups; `csample` uses `count(#)` where native `sample` uses `# , count`

```stata
* Draw a 50% sample within each foreign category
csample 50, by(foreign)

* Draw exactly 5 observations from each group
csample, count(5) by(foreign)
```

---

## cqreg

`absorb()` is currently unavailable: least-squares partialling does not estimate
fixed-effect quantile regression. Use explicit indicators, such as
`cqreg ln_wage age tenure i.idcode i.year`, when the number of groups permits it.

---

## cbinscatter

The original `binscatter` package uses the "classic" method for covariate control, which residualizes both x and y before creating bins.

**ctools adds:**

- **Binsreg controls**: `method(binsreg)` implements the method to control for covariates described by Cattaneo et al. (AER 2024). Cattaneo et al. (2024) make the point that if covariates enter nonlinearly into a given conditional expectation function, then the Frisch-Waugh method of partialing out Y and X with respect to these covariates and running the binned scatterplot on the residuals is not easily interpretable. They propose a simpler method, which is basically to construct bins of X based on the observed (non-residualized) distribution of X, regress Y on X and Z, and evaluate predictions holding fixed Z (e.g. evaluate the model with Z = 0, X = means within bin).

```stata
* Use the binsreg covariate control method
cbinscatter y x, controls(z) method(binsreg)
```
