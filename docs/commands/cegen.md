# cegen

Create a variable from a function of the data, with the syntax, functions,
options, storage types, labels, and messages of Stata's `egen`.

`cegen` is pinned to `egen` as it ships with Stata (`egen.ado` 3.4.2, April 13, 2020).
See [Attribution](../ATTRIBUTION.md#cegen) for the authors of the functions.

## Syntax

```stata
cegen [type] newvar = fcn(arguments) [if] [in] [, by(varlist) fcn_options threads(#) verbose]
by varlist: cegen [type] newvar = fcn(arguments) [if] [in] [, fcn_options]
```

`by(varlist)` needs no sort and does not sort the data. The `by` prefix needs
data sorted by `varlist`, as every `by` prefix does. The order of the
observations never changes.

| Option | Meaning |
|---|---|
| `by(varlist)` | Compute the function within groups. Missing by-values form groups, each missing code its own. |
| `threads(#)` | Thread limit (ctools addition). Results do not depend on it. |
| `verbose` | Report the sample, the number of groups, the storage type, and phase times (ctools addition). |

## Functions

Functions marked (by) accept `by()` and the `by` prefix.

| Family | Functions |
|---|---|
| Statistics of an expression (by) | `count()`, `mean()`, `total()` [`missing`], `sd()`, `min()` and `max()` [`missing`], `median()`, `pctile()` [`p(#)`], `iqr()`, `mad()`, `mdev()`, `skew()`, `kurt()` |
| One value per observation (by) | `std()` [`mean(#) sd(#)`], `pc()` [`prop`], `rank()` [`field`, `track`, `unique`], `mode()` [`minmode`, `maxmode`, `nummode(#)`, `missing`], `seq()` [`from(#) to(#) block(#)`] |
| Across variables | `rowtotal()` [`missing`], `rowmean()`, `rowsd()`, `rowmin()`, `rowmax()`, `rowmedian()`, `rowpctile()` [`p(#)`], `rowfirst()`, `rowlast()`, `rowmiss()`, `rownonmiss()` [`strok`], `anycount()`, `anymatch()`, `anyvalue()` [`values(numlist)`], `diff()` |
| Groups, bins, and sequences | `group()` [`missing`, `autotype`, `label`, `label(name, replace truncate(#))`], `tag()` [`missing`], `cut()` [`at(numlist)`, `group(#)`, `icodes`, `label`], `fill()`, `ma()` [`t(#) nomiss`] |
| Strings | `concat()` [`punct() decode maxlength(#) format(%fmt)`], `ends()` [`punct() head last tail trim`] |

The names `egen` kept from before Stata 9 are accepted: `sum()`, `any()`,
`eqany()`, `neqany()`, `rfirst()`, `rlast()`, `rmax()`, `rmean()`, `rmin()`,
`rmiss()`, `robs()`, `rsd()`, and `rsum()`.

An argument of `exp` is any Stata expression. A variable name is read as it
is. Any other expression costs one `generate` over the data before the plugin
runs, within by-groups when it refers to `_n`, `_N`, a subscript, or `sum()`.

## How it works

Each call makes two plugin calls. The first loads the sample (sources,
by-variables, and the sample marker in one parallel pass), computes the
function, and keeps the result. The ado-file then creates the new variable
with its final storage type and without initial values, and the second call
writes every observation in parallel. No temporary variables, sorts, or
`replace` passes run in Stata.

By-groups are found three ways, cheapest first:

1. Data already in by-variable order (the usual panel layout): one parallel
   pass over neighbouring rows.
2. Numeric by-variables in any order: the groups are numbered through the
   parallel remap the HDFE commands use (a direct table for dense integer
   codes, hash partitions for anything else). The moment functions then add
   each value to its group's sum in one pass, with no sort. Order statistics,
   `rank()`, and `mode()` gather each group with one counting pass and sort
   or select within groups in parallel.
3. String by-variables: the stable parallel sorts of `csort`.

In every case the rows of a group are taken in observation order, so a sum is
the same on every run and for any number of threads.

## Results compared with egen

`cegen` follows `egen`'s arithmetic as well as its definitions.

- **By-groups.** `egen` adds up each group in the order of the observations.
  `cegen` does the same, so on data sorted by the by-variables the two give
  identical values. On unsorted data `egen` sorts first, and Stata's sort
  orders ties at random, so `egen`'s own sums can differ in the last digits
  from run to run; `cegen` gives the same value every time.
- **No by-groups.** `egen` takes means and sums from `summarize`. `cegen`
  sums fixed blocks of rows in parallel with compensated summation. The two
  agree to rounding error; on the validation data, and on a million normal
  draws, every moment function gave identical values.
- **Percentiles.** `egen` uses `_pctile` without by-groups and a formula on
  the sorted values with them. The two rules differ for a few combinations of
  group size and percentile (the 33.3rd percentile of 3,000 values is one).
  `cegen` keeps both rules.
- **Storage types.** The new variable has the type `egen` gives it, including
  the cases where `egen` builds the variable with `replace` and Stata widens
  an integer type: `egen byte x = rowtotal(a b)` becomes `float` when a sum
  has a fraction, and its earlier steps were truncated. `cegen` follows the
  same steps.
- **Labels and messages.** Variable labels, value labels of `group()` and
  `cut()`, notes, display formats, and the "missing values generated" message
  are `egen`'s. The warnings of `mode()` say the same in other words.

## Differences from egen

- `rank(), unique`: `egen` breaks ties in the order its sort leaves them,
  which changes between runs. `cegen` breaks them in observation order.
- `seq()` with `by()` and `if`: `egen` numbers a group's observations in the
  order its sort leaves them. `cegen` numbers them in observation order,
  which is what `egen` gives on data sorted by the by-variables.
- `mtr()`, functions from other packages (such as `egenmore`), `rowfirst()`
  and `rowlast()` of string variables, and results or keys that need a `strL`
  longer than 2,045 bytes are run by `egen` itself, with the same arguments.
- In `concat()`, numbers other than integers of up to seven digits, and any
  number with `format()` or `decode`, are converted to text by Stata's
  `string()` and `decode` before the plugin joins them.
- The undocumented `suffix()` option of `rank()` is not supported.
- `cegen` does not sort the data, so it never changes the sort order or the
  sort marker. `egen` restores the order after its own sorts.

## Performance

Measured on one Windows machine (16 cores) with other work running, best of
two runs, 20-observation groups unless noted. These timings are not part of
the README benchmark profile.

| Call | `egen` 10M | `cegen` 10M | Speedup | vs `gegen` at 2M |
|---|---:|---:|---:|---:|
| `mean(x)` | 0.31 s | 0.08 s | 4× | 3× |
| `median(x)` | 7.5 s | 0.17 s | 44× | 3× |
| `mean(x), by(g)`, unsorted | 5.0 s | 0.14 s | 36× | 5× |
| `sd(x), by(g)`, unsorted | 5.1 s | 0.16 s | 32× | 5× |
| `median(x), by(g)`, unsorted | 5.8 s | 0.24 s | 24× | 5× |
| `rank(x), by(g)`, unsorted | 6.3 s | 0.25 s | 25× | 5× |
| `mode(k), minmode by(g)` | 13.5 s | 0.23 s | 59× | 52× |
| `mean(x), by(g)`, sorted | 0.53 s | 0.12 s | 4× | 3× |
| `median(x), by(g)`, sorted | 2.0 s | 0.12 s | 17× | 4× |
| `group(g h)` (8.2M groups) | 7.8 s | 0.72 s | 11× | 1.3× |
| `tag(g)` | 7.2 s | 0.11 s | 63× | 4× |
| `rowtotal(a1 a2 a3 a4)` | 0.45 s | 0.17 s | 3× | 3× |
| `rowmedian(a1 a2 a3 a4)` | 9.6 s | 0.20 s | 48× | 41× |
| `cut(x), group(10)` | 7.6 s | 0.28 s | 27× | 17× |
| `concat(k h)` | 6.4 s | 0.15 s | 42× | 27× |

`gegen` passes the functions it does not implement to `egen`; those rows show
`cegen` against `egen`'s own speed.

Most of what remains is moving the data in and out of Stata: at 10M
observations a column takes about 25 ms to read and 30 ms to write, and
`mean(x), by(g)` spends 35 ms computing.

## Examples

```stata
cegen mean_price = mean(price), by(foreign)
cegen double share = pc(price * 1.1) if mpg > 15, by(foreign rep78)
cegen p90 = pctile(price), p(90) by(rep78)
cegen r = rank(mpg), field by(foreign)
cegen long id = group(foreign rep78), label
cegen total = rowtotal(price weight length)
cegen q = cut(price), group(5)
cegen key = concat(foreign rep78), punct("-")
bysort foreign (price): cegen run = total(price[_n-1])
cegen sd_price = sd(price), by(foreign) threads(4) verbose
```

## Validation

`validation/stata/cegen.do` compares `cegen` with `egen` for every function
and option: values, storage type, display format, variable label, and value
label, or the return code when `egen` refuses the call. It covers data with
and without by-groups, `if` and `in`, sorted and unsorted, extended missing
values, string and numeric keys, and 250,000 observations for the threaded
loops, where results from one and four threads must be identical.

```sh
python3 validation/run.py cegen
```
