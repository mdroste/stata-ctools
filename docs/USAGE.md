# Usage notes

- All ctools programs follow a basic structure: (1) copy data from Stata to C; (2) operate on that data (3) return data from C to Stata. *This means that ctools programs require more memory than the programs they replace.* In addition, some commands will run faster when they involve fewer variables, or when you have fewer variables in memory. For instance,  `csort`'s runtime is heavily dependent on the number of variables in memory, and can be slower than Stata's built-in `sort` if the dataset is relatively wide (e.g. 100+ variables) due to this data transfer overhead.
- The default options for `csort` and `cmerge` require a lot of memory. For `csort`, you probably require 2-3 times as much memory as your dataset. For `csort`, you can reduce this memory overhead with the optional argument `stream(#)`, which reads variables a handful at a time rather than all at once (at a modest cost to runtime).

### Interpolation: `cinterpolate`

`cinterpolate` defaults to Stata's linear `ipolate`, including grouped interpolation and optional extrapolation. The `cubic` option implements local cubic interpolation as in Nicholas J. Cox's SSC [`cipolate`](https://www.stata.com/statalist/archive/2002-07/msg00088.html). Neither reference command is required to run `cinterpolate`.

```stata
cinterpolate yvar xvar [if] [in], generate(newvar) ///
    [by(varlist) epolate cubic threads(#) verbose]
```

- `generate(newvar)` is required and creates a double-precision result. The source variables and observation order are preserved.
- `by(varlist)` interpolates separately within numeric or string groups. The `by:` prefix also works; use either the prefix or `by()`. Missing group values participate, with `.`, `.a`, etc. treated as distinct groups. Group strings must fit within 2,045 bytes.
- At duplicate `xvar` values, nonmissing `yvar` values are averaged. By default, missing values between known points are filled by linear interpolation. Missing `xvar` values and observations excluded by `if`/`in` receive missing results.
- `cubic` fits a polynomial through the nearest two distinct known x points on each side of each gap. It leaves gaps missing when either side has fewer than two anchors, including gaps in the first and last known intervals. This is local Lagrange interpolation, not a spline; it can overshoot. Duplicate x values are averaged before choosing distinct anchors. Unlike legacy SSC `cipolate`, duplicate rows do not create repeated anchors, and extended missing values are treated as missing. `cubic` cannot be combined with `epolate`.
- Without `epolate`, results outside the known range are missing. `epolate` extends the line through the nearest two distinct known x points at each end; groups with fewer than two known points cannot extrapolate.
- `threads(#)` accepts a nonnegative integer and caps C worker threads; `threads(0)`, the default, uses the automatic setting. `verbose` reports phase timings.

```stata
* Fill missing monthly sales separately for each firm.
cinterpolate sales month, generate(sales_i) by(firm)

* Also extrapolate at the ends of each firm's observed series.
cinterpolate sales month, generate(sales_e) by(firm) epolate threads(8)

* Use four neighboring known points for local cubic interpolation.
cinterpolate sales month, generate(sales_c) by(firm) cubic
```

Use `cinterpolate` in scripts that previously called the ctools linear command `cipolate`. ctools no longer ships `cipolate.ado` or `cipolate.sthlp`, leaving that name to the SSC package. When updating an existing installation, remove stale ctools copies of those two files from its ado directory; use `which cipolate` to locate them and check the file header before removal.

The output variable is the result; no command-specific stored results are promised. See `help cinterpolate` or the [full help file](../build/cinterpolate.sthlp).

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

The returned macros `r(varlist)` and `r(nvars)` contain the output names and count; scalar `r(k_new)` also contains the count. An empty selected sample returns error 2000. See `help csplit` or the [full help file](../build/csplit.sthlp).

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

There is no `if`/`in` syntax: filter the master dataset beforehand. The joined dataset is the result; no command-specific stored results are promised. See `help crangejoin` or the [full help file](../build/crangejoin.sthlp).

## Numerical precision

Accumulated scalar statistics (such as total, model, and residual sums of squares)
in `cqreg`, `creghdfe`, and `civreghdfe` may match their reference commands to only
about seven significant digits. See [compatibility](COMPATIBILITY.md) for
command-specific numerical differences.

## Distinct values (cdistinct)

```stata
cdistinct firm year
cdistinct firm year, joint missing
```

`cdistinct` excludes missing values by default and examines each variable
separately unless `joint` is specified. Add `missing` to include missing values.
See [cdistinct](commands/cdistinct.md) for output and saved-result details.

## Removing fixed effects (chdfe)

```stata
chdfe sales capital, absorb(firm year) generate(within_) sample(used)
chdfe sales capital [aw=weight], absorb(firm##c.year) clear keepvars(firm year)
```

Use `generate()` to keep the original data or `clear` to retain only the
transformed sample and requested variables. The command reports absorbed
degrees of freedom in `e(df_a)`. See [chdfe](commands/chdfe.md) for weights,
slopes, solver options, and saved results.
