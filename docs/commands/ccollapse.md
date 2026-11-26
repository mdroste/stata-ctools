# ccollapse

`ccollapse` replaces the dataset in memory with a dataset of summary statistics, one observation per `by()` group. It takes the syntax of Stata's `collapse` and produces the dataset `collapse` produces. The statistics are computed in parallel in the ctools plugin.

`ccollapse` reimplements `collapse` as shipped with StataNow 19.5 (`collapse.ado` 6.2.2); see [Attribution](../ATTRIBUTION.md#ccollapse).

## Syntax

```stata
ccollapse clist [if] [in] [weight] [, by(varlist) cw fast clabel frame(name) threads(#) verbose]
```

`clist` is any combination of

```stata
[(stat)] varlist
[(stat)] target_var = varname
```

with `stat` one of `mean` (the default), `median`, `p1` to `p99`, `sd`, `semean`, `sebinomial`, `sepoisson`, `sum`, `rawsum`, `count`, `percent`, `max`, `min`, `iqr`, `first`, `last`, `firstnm` and `lastnm`. Variables in `clist` may carry time-series operators. `aweight`s, `fweight`s, `iweight`s and `pweight`s are allowed, with `collapse`'s restrictions on the standard-error statistics.

| Option | Meaning |
| --- | --- |
| `by(varlist)` | Groups over which the statistics are calculated; numeric or string |
| `cw` | Casewise deletion: leave out observations with a missing value in any source variable |
| `fast` | Accepted; no effect (see below) |
| `clabel` | Label targets "stat of varname" (an undocumented option of `collapse`) |
| `frame(name)` | Accepted and unused, as in `collapse` |
| `threads(#)` | Maximum number of threads, 1 to 256 |
| `verbose` | Display the time of each step, the observations, the groups, and how the groups were found |

## Saved results

`ccollapse` leaves these global scalars:

| Scalar | Content |
| --- | --- |
| `_ccollapse_nobs` | Observations used |
| `_ccollapse_ngroups` | Groups |
| `_ccollapse_threads` | Threads used |
| `_ccollapse_method` | How the groups were found: 1 one group, 2 observations already in order, 3 whole-number keys in a small range, 4 hash table, 5 sort |
| `_ccollapse_time_load`, `_time_sample`, `_time_group`, `_time_stats`, `_time_store`, `_time_total` | Seconds per step in the plugin |

## How it works

1. **Read once.** The `by()`, source and weight variables of the selected observations are read in one tiled, parallel pass. An `if` condition is evaluated by Stata into a marker that the plugin reads in parallel; an `in` range and a weight that is a single variable are read as they are.
2. **Find the groups without sorting the data, where possible.**
   - Observations already in `by()` order (panel data, the result of an earlier sort): one parallel pass of adjacent comparisons.
   - Whole-number keys in a small range (years, states, identifiers): each observation's cell in the box of key ranges is its group.
   - Few distinct keys of any type: per-thread hash tables, merged; only the distinct keys are sorted.
   - Otherwise a stable parallel sort of the keys.
3. **Compute.** Sums, moments, extremes and first/last values are accumulated per group in the order of the observations. With the observations in group order every thread takes a run of whole groups; otherwise every thread takes a source variable. Percentiles are selected, not sorted for: threads take groups and run a quickselect narrowed from one requested rank to the next. A group of a million observations or more is selected from by all threads at once, by bracketing each rank with a sorted sample and collecting only the values inside the brackets.
4. **Replace the dataset.** Only now is the dataset changed: it is dropped, `_ccollapse_ngroups` observations are created, the `by()` variables are made anew with their storage type, format, labels, value labels, characteristics and notes, and the plugin writes the keys and the statistics. Dropping the data and starting from the collapsed size is far cheaper than cutting a large dataset down.

Because nothing changes before step 4, an error or Break during the computation leaves the data as they were, and no `preserve` copy is made.

## Differences from collapse

`ccollapse` returns `collapse`'s observations, order, variable names, storage types, formats, labels, sort marker, messages and return codes. The validation suite compares all of these on more than 300 command lines.

- **Exact agreement.** `count`, `percent`, `min`, `max`, `first`, `last`, `firstnm`, `lastnm`, the percentiles and `iqr` are identical. `sum`, `rawsum`, `mean`, `sebinomial` and `sepoisson` are identical on data sorted by the `by()` variables.
- **Order of summation.** On data that are not sorted by the `by()` variables, `collapse` arranges the observations of a group in an order of its own, and that order changes from run to run: two runs of `collapse` on the same unsorted data give sums of inexact numbers that differ in the last digits. `ccollapse` always adds a group's values in the order of the observations, so its result is the same every time and for any number of threads.
- **Standard deviations.** Stata/MP forms the sum of squared deviations of a group in as many pieces as it has processors, so `sd` and `semean` from `collapse` change in the last digit with `set processors`. With one processor `ccollapse` matches `collapse` exactly for unweighted, `fweight`ed and `iweight`ed data; with `aweight`s the two can differ in the last digit.
- **`fast`** has no effect: there is no copy of the data to skip.
- **`sepoisson`** is computed for groups whose total reaches 2^31, where `collapse` stores a value that is not a number.
- **A `by()` variable may also be a source** (under another target name): `ccollapse (count) n = id, by(id)`. `collapse` stops with "failed to setup collapse problem".
- **strL variables.** The plugin cannot write strL values. With a strL `by()` or source variable `ccollapse` calls `collapse`.
- **Old versions.** Under `version 4` or earlier `ccollapse` calls `collapse`, which runs its old command.
- **`graphbar()`**, the undocumented option `graph bar` passes to `collapse`, is not available.

## Performance

Seconds on 10 million observations (Windows, Ryzen 9950X, 16 threads, StataNow 19.5 MP with 4 processors; the smallest of three runs for `ccollapse` and `gcollapse`, one run for `collapse`; `collapse` and `gcollapse` with `fast`). The machine was running other jobs, so read the ratios as approximate.

| Workload | `collapse` | `gcollapse` | `ccollapse` | vs `collapse` | vs `gcollapse` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(mean) x1 x2 x3, by(id)`: 1M groups, data sorted by `id` | 0.40 | 0.95 | 0.12 | 3.2× | 7.7× |
| `(mean) x1 x2 x3, by(year)`: 30 groups | 0.89 | 0.68 | 0.08 | 11× | 8.7× |
| `(mean) (sd) (count) (sum), by(g1 g2)`: 50,000 groups | 0.94 | 0.86 | 0.10 | 9.6× | 8.8× |
| `(mean) x1 x2 (max) x3, by(code)`: 1,000 string keys | 0.71 | 1.41 | 0.10 | 7.1× | 14× |
| `(median) x1 (p90) x2, by(g1)`: 1,000 groups | 0.90 | 0.87 | 0.10 | 9.2× | 8.8× |
| `(mean) x1 x2 (sd) x1 (p50) x2`: no `by()` | 0.12 | 0.43 | 0.11 | 1.1× | 3.9× |
| `(sum) x1 x2, by(big)`: 235,000 sparse integer keys | 1.10 | 0.69 | 0.13 | 8.5× | 5.3× |
| `(mean) x1 x2 (sd) x3 [aw=w], by(g1)` | 2.02 | 1.69 | 0.10 | 21× | 17× |
| Ten statistics, `by(year g2)`: 1,500 groups | 1.19 | 2.24 | 0.11 | 11× | 20× |

A typical call at this size spends about 40 ms reading the data, 10-20 ms finding the groups, 10-30 ms on the statistics and 10-30 ms replacing the dataset.

## Examples

```stata
sysuse auto, clear
ccollapse price mpg, by(foreign)

sysuse auto, clear
ccollapse (mean) price (sd) sd_price=price (median) med_price=price (count) n=price, by(rep78)

webuse census, clear
ccollapse (p25) medage (p75) p75=medage (iqr) iqr=medage [fweight=pop], by(region)

sysuse auto, clear
ccollapse (first) make (last) last_make=make (max) price, by(foreign)
```
