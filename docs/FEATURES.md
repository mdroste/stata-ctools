# ctools New Features

ctools adds the following features to the commands it replaces. See [Compatibility](COMPATIBILITY.md) for supported options and differences.

Most ctools commands use C implementations with parallel processing. `cdecode`
reads value labels directly from Stata and writes the text in C. It calls
Stata's `decode` when the output requires `strL` or when `maxlength()` cuts
through a multibyte character.

---

## Command overview

| Command | Replaces | New features? |
|---------|----------|---------------|
| `csort` | `sort` | Yes |
| `cmerge` | `merge` | — |
| `cimport` | `import` file-format subcommands | Yes |
| `cexport` | `export` file-format subcommands | Yes |
| `cencode` | `encode` | Yes |
| `cdecode` | `decode` | Yes |
| `cdestring` | `destring` | — |
| `cdistinct` | `distinct` | Counts each variable or joint combinations |
| `chdfe` | `hdfe` | Removes categorical effects and group-specific slopes in C |
| `cinterpolate` | `ipolate`, `cipolate` | Preserves input row order |
| `csplit` | `split` | — |
| `crangejoin` | `rangejoin` | — |
| `cwinsor` | `winsor2` | — |
| `csample` | `sample` | Yes |
| `cbsample` | `bsample` | — |
| `crangestat` | `rangestat` | — |
| `creghdfe` | `reghdfe` | — |
| `civreghdfe` | `ivreghdfe` | — |
| `cppmlhdfe` | `ppmlhdfe` | — |
| `cqreg` | `qreg` | Yes |
| `cpsmatch` | `psmatch2` | — |
| `cbinscatter` | `binscatter` | Yes |

`cinterpolate`, `csplit`, and `crangejoin` support `threads()` and `verbose` for
setting the maximum number of threads and reporting time spent in each computation
stage. Their full syntax, examples, and limits are in the
[usage notes](USAGE.md#interpolation-cinterpolate).
`cinterpolate` defaults to linear interpolation; `cubic` selects four-point local
cubic interpolation as in SSC `cipolate`. Cubic mode does not extrapolate.

By default, each command uses at most one thread per available CPU, up to 16
threads. `ctools, threads(#)` changes this for the session and
`ctools, threads(#) permanently` saves it for later sessions;
`CTOOLS_NUM_THREADS` and `OMP_NUM_THREADS` are also honored. On Linux the usable
CPUs follow taskset/Slurm affinity and container CPU limits. See `help ctools`.

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
- **Streaming mode**: `stream()` option reduces memory usage for wide datasets by reordering variables outside the sort key in batches

```stata
* Use counting sort for integer data
csort year, algorithm(counting)

* Streaming mode for memory-constrained sorting
csort id, stream(4)
```

---

## cimport / cexport

`cimport` and `cexport` use C code to read and write delimited text, XLS/XLSX,
dBase, SAS7BDAT, SAS XPORT5/8, SPSS and ESRI shapefiles. Catalog and ZSAV import
are supported; FRED, Haver, and Parquet are not. See
[import and export formats](IO_PARITY.md) for details.

**ctools adds:**

- **Workbook editing**: XLSX and BIFF8 XLS support adding worksheets, `sheet(..., modify/replace)` and `keepcellfmt`.
- **Long text and encodings**: Delimited fields use strL when needed; UTF-8/16/32 and named code pages are supported.
- **Locale-aware number parsing** (cimport delimited): `parselocale()` supports 1,016 locale names and their conventions for digits, decimal and grouping separators, and signs.

```stata
* Import a CSV file
cimport delimited using data.csv, clear

* Import an Excel file with the same command
cimport excel using data.xlsx, clear firstrow

* Import with German number formatting (comma = decimal, period = thousands)
cimport delimited using data_de.csv, clear parselocale(de_DE)

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

`cqreg` adds `vce(cluster clustvar)` for clustered standard errors. See the
[command guide](commands/cqreg.md) for density and bandwidth controls.

`absorb()` is currently unavailable: least-squares partialling does not estimate
fixed-effect quantile regression. Use explicit indicators, such as
`cqreg ln_wage age tenure i.idcode i.year`, when the number of groups permits it.

---

## cbinscatter

The original `binscatter` package uses the "classic" method for covariate control, which residualizes both x and y before creating bins.

**ctools adds:**

- **Binsreg controls**: `method(binsreg)` forms bins from the raw x values and estimates bin indicators jointly with controls and absorbed effects. It evaluates the adjusted bin means at the controls' sample means. This follows the covariate-adjustment approach of Cattaneo et al. (2024); confidence bands and tests are not implemented. See the [command guide](commands/cbinscatter.md).

```stata
* Use the binsreg covariate control method
cbinscatter y x, controls(z) method(binsreg)
```
