# creshape

Convert data from wide to long form and back, with the syntax, messages and
results of Stata's `reshape`.

`creshape` is pinned to `reshape.ado` 4.7.1 (October 4, 2023). See
[Attribution](../ATTRIBUTION.md#creshape).

## Overview

```
 long                          wide
+----------------+            +-------------------+
| i  j   stub    |            | i   stub1  stub2  |
|----------------|  creshape  |-------------------|
| 1  1   4.1     | <--------> | 1   4.1    4.5    |
| 1  2   4.5     |            | 2   3.3    3.0    |
| 2  1   3.3     |            +-------------------+
| 2  2   3.0     |
+----------------+
```

`creshape long` turns one observation per `i` into one per `i` and `j`;
`creshape wide` does the reverse. The result is what `reshape` gives: the same
observations in the same order, the same variables with the same storage
types, display formats, labels, notes and characteristics, and the same
settings recorded with the data. `reshape` and `creshape` can therefore be
mixed on one dataset: `creshape long ...` followed by `reshape wide` works, and
so does the reverse.

## Syntax

```stata
creshape long stubnames, i(varlist) [j(varname [values]) string atwl(chars) ///
    favor(speed|memory) verbose threads(#)]
creshape wide stubnames, i(varlist) [j(varname [values]) string atwl(chars) ///
    favor(speed|memory) verbose threads(#)]

creshape long [, verbose threads(#)]     // back to long, settings from the data
creshape wide [, verbose threads(#)]     // back to wide, settings from the data

creshape i varlist                       // advanced syntax, as in reshape
creshape j varname [values] [, string]
creshape xij fvarnames [, atwl(chars)]
creshape xi [varlist]
creshape [query]
creshape clear
creshape error
```

A stub may contain `@` where the value of `j` appears in the names of the wide
form: `inc@r` stands for `inc80r`, `inc81r`, … and for `incr` in long form.

## Options

| Option | Description |
|--------|-------------|
| `i(varlist)` | The variables that identify a logical observation (required) |
| `j(varname [values])` | The variable holding the sub-observation: new in `creshape long`, existing in `creshape wide`. `values` restricts the values used; without `j()` the variable is `_j` |
| `string` | `j` holds strings |
| `atwl(chars)` | Text that replaces `@` in the names of the long form. `reshape` takes it only in `reshape xij`; `creshape` also takes it in the basic syntax |
| `favor(speed\|memory)` | Accepted so that `reshape` commands run unchanged. It is recorded in `_dta[ReS_favor]` for a later `reshape`; `creshape` has one method |
| `verbose` | Display where the time went |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |

## Examples

```stata
* Wide to long, and back
webuse reshape1, clear
creshape long inc ue, i(id) j(year)
creshape wide

* Long to wide
webuse reshape6, clear
creshape wide inc ue, i(id) j(year)

* String values of j
webuse reshape4, clear
creshape long inc, i(id) j(sex) string

* j in the middle of the variable names
webuse reshape3, clear
creshape long inc@r ue, i(id) j(year)

* Timings, on four threads
creshape wide, verbose threads(4)
```

## What the result looks like

**Observations.** The result is sorted by `i` (and by `j` in long form) and
marked as sorted. `creshape long` gives every `i` one observation for every
value of `j`, with missing values where a stub has no variable for that value.
`creshape wide` gives one observation per `i`, with missing values where an
`i` has no observation for a value of `j`.

**Variables.** In order: the `i` variables; `j` and then the stubs (long), or
the stubs for each value of `j` in turn (wide); then every other variable, in
the order of the dataset.

| Variable | In the result |
|----------|---------------|
| `j` (long) | The smallest storage type that holds its values, with that type's default format |
| A stub's variable (long) | The widest storage type of its wide variables (`byte` < `int` < `long` < `float` < `double`; the longest `str#`); the format and value label of the last of them in order of `j`; no variable label; all their characteristics (a later variable's wins) and their notes, each once |
| A stub's variables (wide) | The long variable's type, format, value label, notes and characteristics; the label is the value of `j` and the name, such as `80 inc` |
| Other variables | Unchanged |

**Labels across a round trip.** `creshape wide` sets aside the variable label
and value label of `j` and the variable labels of the long variables in
characteristics of the dataset (`_dta[__JVarLab]`, `_dta[__JValLabName]`,
`_dta[__JValLab]`, `_dta[__XijVarLab…]`), and `creshape long` puts them back,
as `reshape` does.

**Settings.** `_dta[ReS_i]`, `_dta[ReS_j]`, `_dta[ReS_jv]`, `_dta[ReS_Xij]`,
`_dta[ReS_str]`, `_dta[ReS_atwl]`, `_dta[Res_Xi]` and the summary
characteristics `_dta[ReS_Xij_n]`, `_dta[ReS_Xij_wide#]`, `_dta[ReS_Xij_long#]`
are those of `reshape`. The subcommands that only declare, show or clear them
(`creshape i`, `j`, `xij`, `xi`, `query`, `clear`) and `creshape error` are
carried out by `reshape` itself.

## What is checked

A request that would lose information is refused with `reshape`'s message and
return code, before the data are changed:

| Condition | Return code |
|-----------|-------------|
| Wide form: `i` does not identify the observations | 9 |
| Long form: a value of `j` repeats within `i` | 9 |
| Long form: a variable that is not reshaped varies within `i` | 9 |
| `j` has missing values, or a string `j` has leading or trailing blanks | 498 |
| A string `j` holds a value with a space | 111 |
| No variable carries a stub (`creshape long`) | 111 |
| A stub's wide variables mix strings and numbers | 198 |
| A new variable would take an existing name | 110 |
| A value of `j` makes an invalid variable name (a fraction, a negative number) | 198 |
| `j` takes on more values than `c(maxvar)` (`creshape wide`) | 134 |
| The result would have more variables than `c(maxvar)` | 900 |

`creshape error` then lists the problem observations.

## How it works

1. **Plan.** The plugin copies every observation into its own memory, a whole
   observation at a time, and works out the request from the variable names:
   which variables carry which value of `j`, each stub's variables, and the
   variables that stay. It finds the order of the `i` variables (data already
   in order are not sorted), and for long data the values of `j`, the groups
   of `i` and the constant variables. Every check above happens here.
2. **Description.** Stata drops the observations and saves the empty dataset
   to a temporary file. The plugin reads that description (storage types,
   formats, labels, notes, characteristics, value labels), rewrites it for
   the new form and writes it back; Stata loads it with `use` and sets the
   number of observations. Nothing is done in Stata one variable at a time,
   so a result with thousands of variables costs no more than one with ten.
3. **Store.** The plugin writes the values on every thread, each thread a
   block of observations at a time across all variables.

The description names the sort order, and the plugin writes the observations
in that order, so Stata does not sort or check the result again.

**Memory.** While it runs, `creshape` holds the original data in the plugin
(8 bytes per numeric value, the declared width plus one per string) next to
the reshaped data in Stata.

**If something fails.** From the plan to the end of the store the plugin holds
the only copy of the observations, and `creshape` does not respond to Break in
that time. If Stata cannot build the new form (for lack of memory, for
example), the original data are put back exactly: values, order, sort order,
labels, characteristics, `c(changed)` and `c(filename)`.

## Differences from reshape

- **`strL` and very wide strings.** A plugin cannot write `strL` values. A
  dataset with a `strL` variable, a `str2045` variable, or a string variable
  whose values would take more than 2 GB is reshaped by `reshape` itself; the
  result is the same. `verbose` says when this happens.
- **Large values of `j`.** `creshape long` stores `j` in the smallest type
  that holds its values. `reshape` with `favor(speed)` stores values beyond
  2,147,483,620 as missing.
- **Listed values without observations.** With values listed in `j()`,
  `reshape wide` with `favor(speed)` reports "no data" for every value listed
  after one that has no observations, and leaves out their variables.
  `creshape` leaves out only the values that have no observations.
- **Refused requests.** `reshape` leaves the data sorted by `i`; `creshape`
  leaves them as they were.
- **Data that are already wide.** `creshape wide` on data that are already
  wide displays "(data are already wide)" and returns 0, as `reshape` with
  `favor(memory)` does; `reshape` with `favor(speed)`, the default, stops with
  r(111).
- **`if`, `in` and weights.** `creshape` refuses them with r(101) before any
  setting is recorded. `reshape` records them in `_dta[ReS_Xij]` as if they
  were stubs and then stops with r(110) or r(101).
- **Unattached value labels** are kept; `reshape` drops them.
- **No observations.** `creshape` returns the empty dataset in the new form;
  `reshape` with `favor(speed)` stops with a Mata error.
- **Old syntax.** The syntax of Stata 10 and earlier (`reshape groups`,
  `vars`, `cons`) is passed to `reshape`, as are datasets whose settings were
  declared in it.

## Stored results

`creshape long` and `creshape wide` store in `r()`:

| Result | Description |
|--------|-------------|
| `r(N)`, `r(k)` | Observations and variables of the result |
| `r(N_before)`, `r(k_before)` | Observations and variables before |
| `r(time_load)` | Seconds copying the data into the plugin |
| `r(time_plan)` | Seconds ordering and checking the data |
| `r(time_schema)` | Seconds describing the result |
| `r(time_store)` | Seconds writing the result |
| `r(threads_max)` | Maximum number of threads |
| `r(cmd)` | `creshape` |

## Performance

Seconds for wide data of 100,000 and 1,000,000 observations with ten values
of `j` (1M and 10M observations in long form): one `long` id, three stubs
(`double`, `long`, `str8`) and three other variables, in random order.
StataNow 19.5 MP on a Ryzen 9 9950X under Windows 11, measured November 26,
2026 while other jobs were running; `creshape` is the best of three runs,
`reshape` one run.

| Long observations | Command | `reshape`, `favor(speed)` | `reshape`, `favor(memory)` | `creshape` | Speedup |
|---:|---|---:|---:|---:|---:|
| 1M | `long` | 0.50 | 0.88 | 0.023 | 22× / 38× |
| 1M | `wide` | 0.50 | 1.41 | 0.024 | 21× / 59× |
| 10M | `long` | 5.77 | 9.72 | 0.135 | 43× / 72× |
| 10M | `wide` | 4.71 | 12.96 | 0.126 | 37× / 103× |

`favor(speed)` is `reshape`'s default in StataNow 19.5. With `threads(1)` the
10M cases took 0.80 and 0.92 seconds, with `threads(4)` 0.26 and 0.29.

## Validation

`python3 validation/run.py creshape` compares `creshape` with `reshape` on
labelled, random and large datasets in both directions (values in order, and
everything in the description), checks each refused request's return code and
that the data are untouched, fails the command at each stage after the data
have started to change, and covers `strL` data, empty data, 800 variables,
names outside ASCII and the truth of the sort marker.
