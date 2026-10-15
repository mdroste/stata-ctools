# README benchmark

The [one-shot profiler](../validation/profile_readme.do) compares every command
in the README with its native Stata or community-package reference. Each pair
runs **once** on each of **1,000,000, 10,000,000, and 25,000,000 observations**.
The README reports the minimum and maximum reference/ctools time ratios across
these three sizes, rounded to one decimal place. Values below 1 mean ctools was
slower. These are rough workload-specific measurements, not confidence intervals.

## Running it

From the repository root, with a compiled plugin and ado files in `build/`:

```sh
/bin/zsh -lic 'stata -q -b do "/absolute/path/to/stata-ctools/validation/profile_readme.do"'
python3 validation/summarize_readme_profile.py docs/benchmarks/readme_profile --update-readme README.md
```

On this machine Stata must be invoked through the `stata` alias in an interactive
login zsh. Confirm the invocation on other machines. The driver exits with
`exit, clear`, including after a captured error. The summarizer checks completion
markers, exactly one timing per command/method/size, return codes, basic output
checks, and positive finite elapsed times before updating the table. A failed
pair at any size leaves that command's speedup as a dash. It never publishes a
range based on only the successful sizes.

Reference packages must already be installed: `rangejoin`, `rangestat`, `gtools`
(`gstats`), `psmatch2`, `binscatter`, `reghdfe`, `ivreghdfe`, `ppmlhdfe`, and their
dependencies (including `ftools`, `ivreg2`, and `ranktest`). The profiler does not
install or update packages. Missing commands produce explicit nonzero return
codes in the CSV and log; other commands continue.

Optional positional arguments are the output directory, a quoted list of sizes,
and the Stata processor count. Defaults are `docs/benchmarks/readme_profile`,
`"1000000 10000000 25000000"`, and the current Stata processor setting. The
processor count applies to Stata (`set processors`) and so to the reference
commands. ctools commands run without a `threads()` option and use the ctools
default: every logical CPU. Sizes must be multiples of 200 and at least 1,000.
For a quick check:

```sh
/bin/zsh -lic 'stata -q -b do "/absolute/path/to/stata-ctools/validation/profile_readme.do" "/private/tmp/ctools-profile-smoke" "10000"'
python3 validation/summarize_readme_profile.py /private/tmp/ctools-profile-smoke --sizes 10000
```

The profiler replaces files in its output directory. Save previous results
elsewhere before rerunning. It writes `profile.log` (complete Stata output) and
`timings.csv` (each command's time, return code, and check return code). The
summarizer adds `results.md`; without `--update-readme`, it only writes the results
and prints ranges. Nonstandard sizes cannot update the README. Large temporary
datasets and CSVs are deleted after each size. Allow several GB of RAM and tens
of GB of scratch space; native regression and matching can take substantially
longer than the data operations.

## Data and timed work

Every command starts from the same saved, unsorted **20-variable** dataset at
each size: 16 numeric variables spanning byte/int/long/float/double and four
fixed-width strings (`str8`, `str8`, `str16`, `str40`). It contains a 200-month
firm panel, a second categorical effect, continuous regressors and an instrument,
a continuous outcome, a Poisson outcome, treatment, labeled categories, numeric
text, and pipe-delimited text. The seed and sort seed are 95183. The
[do-file](../validation/profile_readme.do) specifies the full generator and exact
commands; no data download is needed.

| Pair | Workload |
|---|---|
| import/export | All 20 columns in comma-delimited text, numeric category codes. Both imports read the native export file; CSV only, not other supported formats. |
| sort | Stable sort on a continuous regressor and unique row ID. |
| merge | Many-to-one merge into a three-column firm lookup with N/200 rows. |
| sample/bsample | Retain 10% without replacement / bootstrap N observations. |
| encode/decode/destring | Four-category strings / four labeled codes / numeric strings. |
| ipolate | Interpolate 30% randomly missing outcomes within firm, with extrapolation. |
| split | Split three pipe-delimited fields. |
| rangejoin | N-row, 20-column inputs, ±0.5 month within firm; using times offset by 0.25, giving N output rows. Select two using payload variables. |
| gstats winsor | Winsorize three regressors at the 1st/99th percentiles in place. |
| rangestat | Mean outcome over a ±2-month window within firm. |
| binscatter | 20 bins; both commands construct the graph with display disabled. |
| reghdfe/ivreghdfe | Three regressors, firm and group effects; IV instruments x1 with z. |
| ppmlhdfe | Poisson outcome, three regressors, firm and group effects. |
| psmatch2 | Estimate logit propensity scores, then one-neighbor ATT matching with replacement. |
| qreg | Median regression of y on three regressors. |

Stata `timer` surrounds the complete command, including ado work and plugin
loading, reads/writes intrinsic to that command, computation, and output.
Dataset generation, saving/reloading starting data, and basic result checks are
outside the timer. There are no warmups or repetitions. The reference runs first,
then ctools, with input and RNG state reset before each. Filesystem caches are
not flushed. Package initialization and execution order can affect these
one-shot measurements, especially at the smallest size.

The driver checks output observation counts, imported column counts, and
estimator sample sizes; it records command failures. It is a timing profiler,
not a full numerical-equivalence test. Broader correctness checks live in
`validation/`. Sampling realizations can differ under the same seed, and
estimators can use different algorithms; see [compatibility](COMPATIBILITY.md).
Twenty variables remain in memory even when a command only uses a few of them.
Memory pressure, string widths, group counts, options, and join expansion can
change the results substantially.

## Results and provenance

Measured October 14, 2026 on an Apple M4 Pro (12 logical CPUs, 48 GiB RAM),
macOS 26.7.1, StataNow/MP 19.5 (revision August 12, 2026). Stata ran with
**four processors**, the current Stata license limit, so the reference commands
used at most four cores. ctools ran at its **default thread count: all 12
logical CPUs**. The speedups therefore include the extra cores ctools can use,
not only per-core efficiency. The previous run (October 12) capped ctools at four
threads as well. The plugin and `build/` wrappers were rebuilt from the
working tree with `make macos-arm` immediately before the run; the manifest
records their hashes. The working tree contains uncommitted work, so its HEAD
is not asserted to identify the compiled source.

All 120 timed calls completed with return code zero; all observation-count,
imported-column-count, and estimator-sample-size checks passed. These checks
do not establish full numerical equivalence.

- [Times at each size](benchmarks/readme_profile/results.md)
- [Raw timing and status records](benchmarks/readme_profile/timings.csv)
- [Machine and file identities](benchmarks/readme_profile/manifest.json)

The full log is retained locally at `docs/benchmarks/readme_profile/profile.log`
and is ignored by Git because Stata's `about` output includes license details.
The [earlier three-command report](BENCHMARKS_README_COMMANDS.md) uses different
sizes, dataset widths, and repeated measurements; it is retained as historical
evidence and does not supply the current README ranges.
