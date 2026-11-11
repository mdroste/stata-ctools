# cdistinct

Report the number of observations and distinct values for each variable,
using the syntax and saved results of the user-written `distinct` command.

## Syntax

```stata
cdistinct [varlist] [if] [in] [, missing joint abbrev(#) minimum(#) maximum(#) verbose threads(#)]
```

Omitting the variable list selects all variables. Missing values are excluded
separately for each variable. `missing` includes empty strings and numeric
missing values, with extended missing codes counted separately. `joint`
counts combinations of the full variable list; without `missing`, a missing
value in any key excludes the row. Numeric and string keys can be combined.

`minimum()` and `maximum()` limit displayed rows, not saved results. Reversed
limits are exchanged. `abbrev()` controls the displayed name width, with a
minimum of five. `by:` is supported.

`cdistinct` counts sorted values in one pass. For other data, it records
distinct values in a table and switches to sorting if that table fills.
Counts are exact, and values and observation order are preserved.
strL values of any length, and binary ones, are counted exactly: values the plugin cannot read (over 2,045 bytes, or binary) are counted through their ranks in Stata's sort order.
`threads()` sets the thread limit. `verbose` reports load and count times above
the table.

`in` selects the same original rows for every variable. In `distinct` 1.2.1,
internal sorting can change the selection after the first variable. Use a
fixed sample indicator when comparing that reference command with `cdistinct`.

## Saved results

| Result | Meaning |
|---|---|
| `r(N)` | Eligible observations for the last variable, or joint list |
| `r(ndistinct)` | Distinct values for that variable or list |

An empty sample returns zeros. Display filters do not affect these results.

## Examples

```stata
cdistinct firm year
cdistinct firm year, joint missing
bysort industry: cdistinct firm
```

See `help cdistinct` in Stata.
