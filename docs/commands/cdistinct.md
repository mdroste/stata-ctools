# cdistinct

Report the number of observations and distinct values for each variable,
using the syntax, table, messages, and saved results of `distinct` by Gary
Longton and Nicholas J. Cox.

`cdistinct` is pinned to `distinct` 1.4.3 (September 3, 2023), published in the
*Stata Journal* as package dm0042_5. SSC distributes an older version, 1.2.1
(March 1, 2012), whose table and defaults are different. See
[Attribution](../ATTRIBUTION.md#cdistinct) for details and references.

## Syntax

```stata
cdistinct [varlist] [if] [in] [, missing abbrev(#) joint minimum(#) maximum(#)
    sort(specification) width(#) verbose threads(#)]
```

Omitting the variable list selects all variables. Missing values are excluded
separately for each variable. `missing` includes empty strings and numeric
missing values, with extended missing codes counted separately. `joint`
counts combinations of the full variable list and reports them in one row
labeled `(jointly)`; without `missing`, a missing value in any key excludes
the row. Numeric and string keys can be combined. `by:` is supported.

| Option | Meaning |
|---|---|
| `abbrev(#)` | Abbreviate variable names to at most `#` characters. The default is 32 (names in full); values below 5 are treated as 5. No effect with `joint`. |
| `minimum(#)`, `maximum(#)` | Show only rows whose distinct count is at least, or at most, `#`. The defaults are 1 and no upper bound, so a variable with no nonmissing values is hidden unless `minimum(0)` is given. Reversed limits are exchanged, with a message. If no row is left, a note replaces the table. |
| `sort(specification)` | Sort the rows by `alpha` (name as displayed; A before Z before _ before a before z), `distinct` (distinct count), or `total` (observations counted); add `descending` to reverse. Words may be abbreviated, except that `d` alone is refused as ambiguous. `sort` alone means `sort(distinct)`. No effect with `joint`. |
| `width(#)` | Width of each count column. The minimum and default is 11; above ten million observations the default grows with the number of observations. |
| `threads(#)` | Thread limit (ctools addition). |
| `verbose` | Report load and count times above the table (ctools addition). |

`width()` and `sort` without an argument are accepted by `distinct` 1.4.3 but
are not in its help file.

`cdistinct` counts sorted values in one pass. For other data, it records
distinct values in a table and switches to sorting if that table fills.
Counts are exact, and values and observation order are preserved.
strL values of any length, and binary ones, are counted exactly: values the plugin cannot read (over 2,045 bytes, or binary) are counted through their ranks in Stata's sort order.

## Saved results

| Result | Meaning |
|---|---|
| `r(N)` | Observations counted for the last variable, or for the joint list |
| `r(ndistinct)` | Distinct values for that variable or list |

The results describe the last variable whether or not `minimum()` and
`maximum()` let its row be shown. Both are 0 when no row at all is shown,
which includes a sample with no observations. `joint` with no observations to
count stops with error r(2000) instead.

## Differences from distinct 1.4.3

- `in` selects the same original rows for every variable. In `distinct`,
  internal sorting can change the selection after the first variable. Use a
  fixed sample indicator when comparing that reference command with `cdistinct`.
- With `sort()`, rows that tie on the sort key stay in the order of the
  variable list. `distinct` displays tied rows in an order that can change
  from one run to the next.
- With no observations in memory, `cdistinct` does what both commands do when
  `if` excludes every observation: it displays the note (or a table of zeros
  with `minimum(0)`) and returns 0. `distinct` stops with error r(111) in that
  case unless `joint` is specified.
- `distinct` leaves its table behind in Mata as the global matrices `mout`,
  `ndistinct`, and `ntotal`. `cdistinct` creates no Mata objects.
- The package of `distinct` also contains `distinctgen`, which generates a
  variable holding the number of distinct values. ctools does not reimplement
  it; use `distinctgen` itself.

## Examples

```stata
cdistinct firm year
cdistinct firm year, joint missing
cdistinct, sort(distinct descending)
cdistinct, maximum(10) sort(alpha)
bysort industry: cdistinct firm
```

See `help cdistinct` in Stata.
