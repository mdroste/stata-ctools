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
- `by(varlist)` interpolates separately within numeric or string groups. The `by:` prefix also works; use either the prefix or `by()`. Missing group values participate, with `.`, `.a`, etc. treated as distinct groups. strL groups may hold text of any length and binary values.
- At duplicate `xvar` values, nonmissing `yvar` values are averaged. By default, missing values between known points are filled by linear interpolation. Missing `xvar` values and observations excluded by `if`/`in` receive missing results.
- `cubic` fits a polynomial through the nearest two distinct observed x values on each side of a gap. It averages nonmissing y values at repeated x values before selecting those four points. Gaps remain missing when either side has fewer than two points, including gaps in the first and last known intervals. This is local Lagrange interpolation, not a spline; it can overshoot. Extended missing values are treated as missing. `cubic` cannot be combined with `epolate`.
- Without `epolate`, results outside the known range are missing. `epolate` extends the line through the nearest two distinct known x points at each end; groups with fewer than two known points cannot extrapolate.
- `threads(#)` accepts a nonnegative integer and sets the maximum number of threads; `threads(0)`, the default, uses the automatic setting. `verbose` reports time spent in each computation stage.

```stata
* Fill missing monthly sales separately for each firm.
cinterpolate sales month, generate(sales_i) by(firm)

* Also extrapolate at the ends of each firm's observed series.
cinterpolate sales month, generate(sales_e) by(firm) epolate threads(8)

* Use four neighboring known points for local cubic interpolation.
cinterpolate sales month, generate(sales_c) by(firm) cubic
```

Use `cinterpolate` in scripts that previously called the ctools linear command `cipolate`. ctools no longer ships `cipolate.ado` or `cipolate.sthlp`, leaving that name to the SSC package. When updating an existing installation, remove stale ctools copies of those two files from its ado directory; use `which cipolate` to locate them and check the file header before removal.

The result is saved in the output variable. The command does not store command-specific results in `r()` or `e()`. See `help cinterpolate` or the [full help file](../build/cinterpolate.sthlp).

### String fields: `csplit`

`csplit` implements `split` for fixed-width strings, including multiple delimiters, field limits, and optional numeric conversion. It preserves the source variable and observation order.

```stata
csplit strvar [if] [in] [, generate(stub) parse(parse_strings) notrim ///
    limit(#) destring force float ignore(strings) percent threads(#) verbose]
```

- `generate(stub)` names outputs `stub1`, `stub2`, etc.; the default stub is the source variable's name. Output names must be new. The number of fields is the maximum needed across selected observations; unselected observations have empty string outputs (missing after numeric conversion).
- With no `parse()`, spaces separate fields and outer spaces are removed from the remaining text after each field is split off. `parse("|" ";" "::")` supplies literal delimiter strings, which may contain multiple characters. If delimiters match at the same position, the last listed one takes precedence. This is literal splitting, without CSV-style quoted-field handling.
- By default, outer spaces are trimmed before parsing. With explicit delimiters, spaces within the remaining text are preserved. `notrim` retains outer spaces too and requires explicit `parse()`. Adjacent delimiters preserve empty fields; a trailing delimiter does not create another final field.
- `limit(#)` must be positive and retains at most that many fields; text after those fields is discarded.
- `destring` applies Stata's native numeric conversion to the fields. Its options `force` (unconvertible values become missing), `float` (float storage), `ignore(strings)` (remove specified characters), and `percent` (convert percentage values) require `destring`. Without `force`, a field containing nonnumeric values can remain a string, as with `destring`.
- `threads(#)` accepts a nonnegative integer and sets the maximum number of threads, with zero selecting the default; `verbose` reports time spent in each computation stage. The C implementation handles fixed-length strings up to 2,045 bytes. For `strL` input, `csplit` calls Stata's `split`. String bytes are preserved without Unicode normalization.

```stata
* Example code: NY|retail|SKU000042|12345678|USD|2024-01-31
csplit code, generate(part) parse("|")

* Convert comma-separated numeric fields, ignoring dollar signs.
csplit amounts, generate(amount) parse(",") destring ignore("$")
```

The returned macros `r(varlist)` and `r(nvars)` contain the output names and count; scalar `r(k_new)` also contains the count. An empty selected sample returns error 2000. See `help csplit` or the [full help file](../build/csplit.sthlp).

### Interval joins: `crangejoin`

`crangejoin` matches observations by numeric intervals using the syntax of Robert Picard's SSC `rangejoin` 1.1.3. Each observation in memory is matched to every using observation whose numeric key falls within its **inclusive** interval and whose `by()` values agree. The result replaces the dataset in memory. Neither `rangejoin` nor `rangestat` needs to be installed to use it.

```stata
crangejoin keyvar low high using filename [, by(varlist) keepusing(varlist) ///
    prefix(string) suffix(string) all threads(#) verbose]
```

- `keyvar` must be numeric in the using dataset. Numeric master variables supplied as `low` and `high` give **absolute bounds**. Literal numbers are offsets from the same-named master key when it exists, and absolute bounds otherwise. For example, `date -1 1` matches dates within one day of the master date.
- A missing bound (`.` or a missing value in a bound variable) leaves that side unbounded. Missing using keys never match. If the master contains `keyvar`, missing master keys never match either, even with explicit bound variables. Reversed intervals do not match. If there are no valid master intervals or no nonmissing using keys, the command returns error 2000.
- `by(varlist)` restricts matches to equal numeric or string groups. Missing group values participate, and extended numeric missings remain distinct. The key cannot also appear in `by()`.
- `keepusing(varlist)` selects variables to retain from the using dataset; key and group variables are included automatically. Using variables that share master names receive `prefix()` and/or `suffix()`; the default is `suffix(_U)` when neither is specified. `all` applies renaming to every using variable except the group variables. Conflicting output names return an error.
- Unmatched master rows remain, with missing using values. Master observation order is preserved; matches within each master row are ordered by using key, then by original using order for ties. The wrapper preserves variable metadata and restores the master dataset on failure.
- `threads(#)` accepts a nonnegative integer and sets the maximum number of threads, with zero selecting the default; `verbose` reports time spent in each computation stage. Numeric and fixed-width strings are supported. `strL` master and using variables are carried through the join; as in `rangejoin`, `by()` variables cannot be `strL`. Inputs and output must each have at most 2,147,483,647 rows. Memory use and runtime grow with the number of matches, which can greatly exceed input size.

```stata
* Match firm offers within one day of each inquiry's date.
crangejoin date -1 1 using offers.dta, by(firm) keepusing(price venue)

* Match events to observation-specific absolute date windows.
crangejoin date window_start window_end using events.dta, by(firm)

* Attach all history within each person; watch the resulting match count.
crangejoin date . . using history.dta, by(person) suffix(_history)
```

There is no `if`/`in` syntax: filter the master dataset beforehand. The result replaces the dataset in memory. The command does not store command-specific results in `r()` or `e()`. See `help crangejoin` or the [full help file](../build/crangejoin.sthlp).

### Nearest neighbors and distances: `cgeonear`

`cgeonear` finds the nearest neighbors of each location by geodetic distance, with the syntax of Robert Picard's SSC `geonear` 2.0.3, and computes the distance between two points, with the syntax of his `geodist` 1.1.0: a `using` file selects the first form, four coordinates the second. Coordinates are signed decimal degrees. Neither `geonear` nor `geodist` needs to be installed.

```stata
cgeonear baseid baselat baselon using nborfile, neighbors(nborid nborlat nborlon) ///
    [wide long nearcount(#) ignoreself genstub(str) within(#) limit(#) ///
    miles radius(#) ellipsoid a(#) f(#) threads(#) verbose]
cgeonear lat1 lon1 lat2 lon2 [if] [in], generate(newvar) [miles sphere radius(#) ellipsoid(#1,#2)]
```

- Wide form (the default) adds the nearest neighbor's identifier and distance (`nid`, `km_to_nid`, or `nid1 km_to_nid1 ...` with `nearcount()`) to every observation of the data in memory. Long form (`long`) replaces the data with one observation per base location and neighbor; `within()` adds every neighbor within a distance and `limit()` caps the count.
- The nearest-neighbor form measures great-circle (haversine) distances by default and on the WGS 1984 ellipsoid (Vincenty) with `ellipsoid`; the distance form defaults to the ellipsoid, as `geodist` does. Distances agree with `geonear`'s and `geodist`'s to the last digit, except where noted in the [guide](commands/cgeonear.md#accuracy).
- Neighbors at the same distance are ranked by ascending identifier. Each identifier is one location; observations with missing coordinates of a located identifier get its neighbors.

```stata
* The nearest other cell tower of each tower, and every tower within 50 km
cgeonear ctid ctlat ctlon using cell.dta, neighbors(ctid ctlat ctlon) ignoreself
cgeonear ctid0 ctlat ctlon using cell.dta, n(ctid ctlat ctlon) ignoreself long within(50) near(0)

* Distances from each park to a stadium, on the ellipsoid and on a sphere
cgeonear 42.265837 -83.748696 lat lon, generate(d)
cgeonear 42.265837 -83.748696 lat lon, generate(d_sphere) sphere
```

See `help cgeonear` or the [guide](commands/cgeonear.md).

### Pairwise joins: `cjoinby`

`cjoinby` forms all pairwise combinations within groups, with the syntax of Stata's `joinby`. Each observation in memory is paired with every observation of the using file that has the same values of the key variables.

```stata
cjoinby [varlist] using filename [, unmatched(none|both|master|using) update replace ///
    _merge(varname) nolabel threads(#) verbose snapshot fast]
```

Without a varlist the join is on the variables the two datasets share. Observations found in only one dataset are left out unless `unmatched()` keeps them; `_merge` (created with `unmatched()` or `_merge()`) marks the source of each observation. Variables that both datasets have keep the master's values unless `update` (fill in missing values) or `update replace` (also overwrite with nonmissing values) is specified.

```stata
* Every child with every parent of the same family.
use child, clear
cjoinby family_id using parent.dta

* Keep unmatched observations from both datasets and mark them.
cjoinby family_id using parent.dta, unmatched(both)

* Join on two keys; fill in missing values of the other shared variables.
cjoinby firm year using deals.dta, update

* Join on every variable the datasets have in common.
cjoinby using deals.dta
```

Each master observation is followed by its pairings in the using file's order, and a sorted master stays sorted, as after `joinby`. Where `joinby`'s order within ties depends on `set sortseed`, `cjoinby`'s does not vary. The result replaces the dataset in memory; `r(N)`, `r(N_1)`, `r(N_2)` and `r(N_3)` hold the observation counts. See `help cjoinby` or the [guide](commands/cjoinby.md).

### Fuzzy joins: `cfuzzyjoin`

`cfuzzyjoin` joins records whose keys are similar rather than equal. Its three subcommands follow three community commands: `link` is `reclink` (record linkage from several variables), `match` is `matchit` (every pair of similar strings of two datasets, or the similarity of two columns) and `group` is `strgroup` (groups of similar strings by Levenshtein distance).

```stata
cfuzzyjoin link varlist using filename, idmaster(varname) idusing(varname) gen(newvar) [options]
cfuzzyjoin match idmaster txtmaster using filename, idusing(varname) txtusing(varname) [options]
cfuzzyjoin match varname1 varname2 [, options]
[by varlist:] cfuzzyjoin group varname [if] [in], generate(newvar) threshold(#) [options]
```

```stata
* The best match of each firm on name, address and ZIP code.
use firms_a, clear
cfuzzyjoin link name address zip using firms_b, idmaster(id) idusing(rid) gen(score)

* Every pair of names with a bigram similarity of at least 0.7.
use names_a, clear
cfuzzyjoin match id name using names_b, idusing(rid) txtusing(rname) threshold(0.7)

* Group the spellings of firm names within each state.
bysort state: cfuzzyjoin group name, generate(firm) threshold(0.2)
```

`link` adds the score, the matched identifier, the using variables and `_merge` to the data in memory; `match` with `using` replaces the data with the pairs; `group` adds a group number. See `help cfuzzyjoin` or the [guide](commands/cfuzzyjoin.md).

### Reshaping: `creshape`

`creshape` converts data from wide to long form and back, with the syntax of Stata's `reshape`. In wide form each `i` has one observation and one variable per value of `j` (`inc80`, `inc81`, …); in long form each `i` and `j` has one observation.

```stata
creshape long stubnames, i(varlist) [j(varname [values]) string atwl(chars) ///
    favor(speed|memory) verbose threads(#)]
creshape wide stubnames, i(varlist) [j(varname [values]) string atwl(chars) ///
    favor(speed|memory) verbose threads(#)]
creshape long | wide [, verbose threads(#)]
```

```stata
* Wide to long: inc80 inc81 inc82 become inc, with the new variable year
webuse reshape1, clear
creshape long inc ue, i(id) j(year)

* And back, from the settings stored with the data
creshape wide

* Long to wide, string values of j, j inside the names
creshape wide inc ue, i(id) j(year)
creshape long inc, i(id) j(sex) string
creshape long inc@r ue, i(id) j(year)
```

The result is `reshape`'s: the same observations, sorted by `i` (and `j`), and the same variables with their storage types, formats, labels, notes and characteristics. The settings are the dataset characteristics `reshape` uses, so the two commands can be mixed, and `creshape i`, `j`, `xij`, `xi`, `query`, `clear` and `error` work as in `reshape`. A request that would lose information (a repeated `j` within `i`, a variable that varies within `i`) is refused before the data change. Data with `strL` variables are reshaped by `reshape` itself. See `help creshape` or the [guide](commands/creshape.md).

### SQL queries: `csql`

`csql` runs SQL in DuckDB, inside Stata, on the data in memory (table `stata`), on Stata datasets read in place, and on Parquet, CSV and JSON files. It filters, joins and aggregates files larger than memory and loads (`clear`), saves (`saving()`) or lists the result. Install DuckDB once with `csql, install`.

```stata
csql using "transactions.parquet", query("SELECT firm_id, year, SUM(sales) AS sales FROM data GROUP BY ALL") clear
csql, query("SELECT s.*, f.industry FROM stata s JOIN 'firms.dta' f USING (firm_id)") clear
csql, query("SELECT * FROM read_csv('raw/*.csv')") saving("all.parquet", replace)
```

Results take the narrowest Stata types and keep the labels and formats of the variables they pass through. See `help csql` or the [guide](commands/csql.md).

## Numerical precision

Accumulated scalar statistics (such as total, model, and residual sums of squares)
in `cqreg`, `creghdfe`, and `civreghdfe` may match their reference commands to only
about seven significant digits. See [compatibility](COMPATIBILITY.md) for
command-specific numerical differences.

## Extensions to generate (cegen)

```stata
cegen mean_sales = mean(sales), by(firm)
cegen double share = pc(sales) if year > 2000, by(firm year)
cegen long id = group(firm year), label
cegen total = rowtotal(q1 q2 q3 q4)
```

`cegen` takes the syntax and functions of `egen`. `by()` does not need sorted
data and does not sort it. An argument that is a variable name is read
directly; any other expression is evaluated by `generate` first. See
[cegen](commands/cegen.md) for the functions and the differences from `egen`.

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
