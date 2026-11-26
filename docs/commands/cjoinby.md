# cjoinby

Form all pairwise combinations within groups, with the syntax, options and
results of Stata's `joinby`.

`cjoinby` is pinned to `joinby.ado` 2.1.0 (December 22, 2022). See
[Attribution](../ATTRIBUTION.md#cjoinby).

## Overview

`cjoinby` joins each observation of the dataset in memory (the master data)
with each observation of the using file that has the same values of the key
variables. A group with *m* master and *u* using observations gives *m* × *u*
observations. The result replaces the data in memory.

## Syntax

```stata
cjoinby [varlist] using filename [, unmatched(none|both|master|using) update replace ///
    _merge(varname) nolabel threads(#) verbose snapshot fast]
```

Without a varlist, the join is on the variables the two datasets have in
common, which `cjoinby` names. Missing key values form groups like any other
value, and each extended missing value (`.a`, `.b`, …) is a value of its own.
Numeric keys compare as doubles, strings by their bytes. Key variables cannot
be `strL`.

## Options

| Option | Description |
|--------|-------------|
| `unmatched(none)` | Leave out observations found in only one dataset (the default) |
| `unmatched(both)` | Keep them from both datasets |
| `unmatched(master)` | Keep them from the data in memory |
| `unmatched(using)` | Keep them from the using file |
| `update` | A missing master value of a variable both datasets have takes the using value |
| `replace` | With `update`: a nonmissing using value also replaces a nonmissing master value |
| `_merge(varname)` | Name of the variable marking each observation's source (default `_merge`) |
| `nolabel` | Do not copy value-label definitions from the using file |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |
| `verbose` | Display where the time went |
| `snapshot` | Restore the master after a failure from a `preserve` snapshot (written to disk above `c(max_preservemem)`) |
| `fast` | Keep no copy of the master; a failure or Break can leave the data partly joined |

The words of `unmatched()` may be abbreviated. As in `joinby`, the `_merge`
variable is created only when `unmatched()` or `_merge()` is specified.

## Examples

```stata
* Every child with every parent of the same family
use child, clear
cjoinby family_id using parent.dta

* Keep the observations without a match and mark the source
cjoinby family_id using parent.dta, unmatched(both)
tabulate _merge

* Join on two keys; fill in missing values of the other shared variables
cjoinby firm year using deals.dta, update

* Join on every variable the datasets have in common
cjoinby using deals.dta

* Two keys, four threads, with a timing table
cjoinby firm year using deals.dta, threads(4) verbose
```

## What the result looks like

**Observations.** Each master observation is followed by its pairings, in the
order of the observations in the using file, and the master observations keep
their order. Observations kept from the using file only come after them, in
key order. When the master is sorted, the result is sorted the same way, as
`joinby` (a `sortpreserve` program) leaves it: the using-only observations are
then placed by their values of the sort variables, after the master
observations they tie with. A using-only observation holds its key values and
missing values in the master's other variables.

`joinby` leaves tied observations in an order that depends on
`set sortseed`. `cjoinby`'s order is the same on every run and for every
thread count.

**Variables.** The master's variables come first, then `_merge` if there is
one, then the variables only the using file has, with their formats, labels,
notes and characteristics.

| Variable | In the result |
|----------|---------------|
| Key variables | The wider of the two storage types; a string's display format widens with it |
| Other variables both datasets have | The wider storage type; the master's values, unless `update` |
| … of using-only observations | Missing, unless `update` (this is `joinby`'s behavior) |
| … string in one dataset, numeric in the other | The master's variable, untouched |
| Variables only the using file has | Added after the master's |

`long` with `float` gives `double`. Key and shared variables gain the using
variable's notes and the characteristics the master's variable lacks; with
`update` they also take its variable label and value label where the master's
variable has none.

**Value labels.** Definitions of the using file are copied unless the master
already defines the name. With `nolabel` only the labels attached to the using
file's key variables are copied, as in `joinby`. `joinby` also defines the
value label `__MERGE` on every non-empty join, with or without a `_merge`
variable, and so does `cjoinby`.

## `_merge` values

| Value | Meaning |
|-------|---------|
| 1 | Observation from the master data only |
| 2 | Observation from the using data only |
| 3 | Pairing of a master and a using observation |

`update` does not add codes 4 and 5 (only `joinby`'s value label mentions
them).

## Stored results

`joinby` stores nothing; `cjoinby` stores the following in `r()`:

| Result | Description |
|--------|-------------|
| `r(N)` | Number of observations in the result |
| `r(N_1)` | Observations from the master data only |
| `r(N_2)` | Observations from the using data only |
| `r(N_3)` | Pairings |
| `r(time)` | Elapsed seconds with `verbose`; missing otherwise |
| `r(using)` | Name of the using file |
| `r(keyvars)` | Key variable names |

## Differences from joinby

- The order of observations within ties is fixed (see above).
- Without a varlist, the common variables are found by their exact names;
  `joinby` accepts a master variable name that abbreviates a using variable.
- A using file without observations gives the result of a join without
  matches. `joinby` stops with r(111).
- A key that is a frame alias in the master data is refused with r(106); run
  `frunalias` first. Alias variables that are not keys follow their link
  variable, as in `joinby`.
- An error or Break leaves the master as it was. `joinby` can leave a
  half-joined dataset (for example after `replace` without `update`).
- `e(sample)` stays with the first pairing of each master observation;
  `joinby` leaves it with one of them.
- `cjoinby` stores `r()` results.

## strL variables

`strL` variables are supported everywhere except as keys. The plugin cannot
write `strL` values, so it writes the master and using observation number of
each output row and the values are copied in Stata afterwards. That step runs
at Stata's speed; joins without `strL` variables do not pay for it.

## Performance

A join is three plugin calls:

1. **Index the using data.** The using file is read into a frame of its own
   (Stata 16 or later) and staged row-major in one pass over its
   observations. When the keys are integers with a small combined range, each
   key tuple is a position in a table of groups: one radix sort of 8-byte
   (position, row) words gives the key order, the groups and the table, and
   keys already in order are only checked. Other keys (non-integer numbers,
   sparse identifiers, strings) are put in stable order by csort's radix
   sorts and indexed by a hash table built in parallel.
2. **Plan the join in the untouched master.** Only the master's key variables
   are read. Every master observation is looked up in parallel; runs of one
   key look up once. The counts give the size and layout of the output, so
   errors stop here with nothing to undo.
3. **Write the output.** The master observations that move are staged
   row-major (observations before the first one that moves stay where they
   are and are never read). Workers then write tiles of output rows, every
   variable of a tile before the next tile, so each Stata observation is
   touched once.

Memory: the staged using data and the staged master observations that move,
at 8 bytes per numeric variable and width + 1 per string, plus 8 bytes per
master observation for the layout. `fast` releases the staged master rows
right after the write instead of keeping them for a rollback.

The speedup in the README comes from the benchmark profile, which joins the
profile's 20-variable data at 1M, 10M and 25M observations with a small using
file: see the [published results](../benchmarks/readme_profile/results.md) and
[how the comparison is run](../performance/README.md#compare-pairwise-joins).

The table below is a separate set of timings on the same machine (Ryzen 9
9950X, StataNow 19.5 MP, 16 threads) for joins of two large files: master and
using data with five variables each, best of three runs for `cjoinby`:

| Join | Master × using → output | `joinby` | `cjoinby` | Speedup |
|------|------------------------|---------:|----------:|--------:|
| Integer key, many-to-many | 1M × 1M → 3.0M | 4.0 s | 0.09 s | 47× |
| … `unmatched(both)` | 1M × 1M → 3.3M | 4.9 s | 0.09 s | 57× |
| … sorted master | 1M × 1M → 3.0M | 3.8 s | 0.09 s | 45× |
| … sorted master, `unmatched(both)` | 1M × 1M → 3.3M | 4.8 s | 0.08 s | 63× |
| `str12` key, many-to-many | 1M × 1M → 3.0M | 5.6 s | 0.10 s | 57× |
| Non-integer `double` key | 1M × 1M → 3.0M | 4.2 s | 0.10 s | 43× |
| Lookup (unique using keys), `unmatched(master)` | 1M × 1M → 1.0M | 2.2 s | 0.05 s | 47× |
| Lookup, matches only | 1M × 1M → 0.8M | 2.0 s | 0.06 s | 36× |
| Integer key, many-to-many | 10M × 10M → 30.0M | 41 s | 0.82 s | 50× |
| … `unmatched(both)` | 10M × 10M → 33.1M | 51 s | 0.79 s | 64× |
| … sorted master | 10M × 10M → 30.0M | 46 s | 0.66 s | 70× |
| … sorted master, `unmatched(both)` | 10M × 10M → 33.1M | 51 s | 0.75 s | 68× |
| `str12` key, many-to-many | 10M × 10M → 30.0M | 85 s | 1.2 s | 72× |
| Non-integer `double` key | 10M × 10M → 30.0M | 53 s | 0.99 s | 53× |
| Lookup (unique using keys), `unmatched(master)` | 10M × 10M → 10.0M | 30 s | 0.65 s | 46× |
| Lookup, matches only | 10M × 10M → 8.3M | 28 s | 0.61 s | 46× |

`verbose` prints the time of every step. For a join this size, reading the
using file and adding the new observations in Stata take about as long as the
plugin's own work.
