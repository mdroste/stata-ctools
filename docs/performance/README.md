# Measuring ctools performance

The [benchmark script](../../validation/performance/profile_readme.do) measures
the speed comparisons in the main README. It runs each ctools command and the
Stata or community command it replaces once on datasets with **1 million,
10 million and 25 million observations**.

The reported speedup is the comparison command's time divided by ctools' time.
The README shows the smallest and largest ratios across the three dataset sizes,
rounded outward: tenths below 2×, whole numbers below 10×, and multiples of ten
thereafter. For example, 48–57× becomes 40–60×. A value below 1× means ctools
was slower. These ranges describe the three measurements; they are not
confidence intervals.

## Run the benchmark

Build the plugin and keep it with the ado files in `build/`. Install the commands
used for comparison: `rangejoin`, `rangestat`, `gtools` (`gstats`), `psmatch2`,
`binscatter`, `reghdfe`, `ivreghdfe` and `ppmlhdfe`, along with their dependencies,
including `ftools`, `ivreg2` and `ranktest`.

In a fresh Stata session, set the working directory to the repository root and
run:

```stata
do "validation/performance/profile_readme.do"
```

The script closes Stata when it finishes, including after an error. It writes
results to `docs/benchmarks/readme_profile`, replacing files at that location.
To summarize them, run this command from the repository root in a terminal:

```sh
python3 -m validation.performance.summarize_readme_profile docs/benchmarks/readme_profile
```

Add `--update-readme README.md` to update the speedup table. The summarizer
requires a completed run, one timing for each command and size, successful
return codes and output checks, and positive finite times. A failed comparison
at any size produces a dash for that command.

## Run a smaller comparison

The do-file accepts four arguments, in order:

| Argument | Default |
| --- | --- |
| Output directory | `docs/benchmarks/readme_profile` |
| Quoted list of observation counts | `"1000000 10000000 25000000"` |
| Number of Stata processors | Current Stata setting |
| Quoted list of command pairs | All twenty pairs |

For a quick check, run this in a fresh Stata session from the repository root:

```stata
capture mkdir "temp"
do "validation/performance/profile_readme.do" "temp/profile-small" "10000"
```

Then summarize the results in a terminal:

```sh
python3 -m validation.performance.summarize_readme_profile temp/profile-small --sizes 10000
```

To select commands, use the names in `timings.csv`, such as `"reghdfe ivreghdfe
ppmlhdfe"`, as the fourth argument. Include `export` whenever selecting `import`:
those comparisons read the CSV written by Stata's export command. Observation
counts must be multiples of 200 and at least 1,000. Only the three standard sizes
can update the main README.

The processor argument sets Stata's `set processors` value. ctools uses its
default thread count, normally limited to 16 available CPUs, subject to session,
environment and runtime settings. Run `ctools, environment_check` to see the
settings before comparing results.

`profile.log` contains Stata's output. `timings.csv` records each command's time,
return code and output-check result. The summarizer adds `results.md` and prints
the speedup ranges. Missing comparison commands are recorded as errors while the
other comparisons continue. The script removes large temporary datasets and CSVs
after each size. The full run needs several GB of RAM and tens of GB of temporary
disk space; regression and matching can take much longer than data operations.

## Data and operations

Each comparison starts from the same saved, unsorted dataset with **20 variables**:
16 numeric variables spanning byte, int, long, float and double types, and four
strings of widths 8, 8, 16 and 40. The data form a 200-month firm panel and include
another categorical effect, regressors, an instrument, continuous and Poisson
outcomes, treatment, labeled categories, numeric text and pipe-delimited text.
Both random seeds are 95183. The do-file generates all the data locally.

| Command pair | Operation |
| --- | --- |
| import/export | All 20 columns as CSV, with numeric category codes. Both imports read Stata's exported file. |
| sort | Stable sort on a continuous regressor and unique row ID. |
| merge | Many-to-one merge with a three-column firm lookup containing N/200 rows. |
| sample/bsample | Retain 10% without replacement / bootstrap N observations. |
| encode/decode/destring | Convert four-category strings / four labeled codes / numeric strings. |
| ipolate | Interpolate 30% randomly missing outcomes within firm, with extrapolation. |
| split | Split three pipe-delimited fields. |
| rangejoin | Join two N-row, 20-column datasets within firm and ±0.5 month. Using-data times are offset by 0.25, producing N rows. Keep two variables from the using data. |
| gstats winsor | Winsorize three regressors at the 1st/99th percentiles in place. |
| rangestat | Mean outcome over a ±2-month window within firm. |
| binscatter | Construct a 20-bin graph with display disabled. |
| reghdfe/ivreghdfe | Two regressors (`x1`, `x2`), firm and group effects; IV instruments `x1` with `z`. |
| ppmlhdfe | Poisson outcome, two regressors (`x1`, `x2`), firm and group effects. |
| psmatch2 | Estimate logit propensity scores, then one-neighbor ATT matching with replacement. |
| qreg | Median regression of y on three regressors. |

Stata's timer measures the complete command: ado code, plugin loading, data
transfers, computation and output. Preparing and reloading data and checking
results happen outside the timer. The comparison command runs first, followed
by ctools, with the data and random-number state reset before each command.
There are no warmup runs or repetitions, and filesystem caches are left in place.
Loading packages and the order of execution can affect the times, especially
for the smallest dataset.

The script checks observation counts, imported column counts and estimation
sample sizes. Numerical comparisons are in `validation/`. Sampling results can
differ under the same seed, and estimators can use different algorithms; see
[compatibility](../COMPATIBILITY.md). All twenty variables stay in memory even
when a command uses only a few. Available memory, string widths, group counts,
options and the number of rows produced by a join can all affect speed.

## Published results

The published measurements use an Apple M4 Pro with 12 logical CPUs and 48 GiB
RAM, macOS 26.7.1, and StataNow/MP 19.5 (revision August 12, 2026). Stata used
**four processors**, its license limit; ctools used **12 threads**, its default
on that machine. The speedups therefore reflect both the implementation and the
additional cores available to ctools.

Each command pair and dataset size ran in a separate Stata session, with
`export` and `import` together. A session began after other processes used less
than 80% of one CPU in total for a minute; sessions exceeding 150% were rerun.
Package initialization and ado-file loading therefore fall inside every timing.
The [measurement record](../benchmarks/readme_profile/manifest.json) contains
these settings, per-session CPU measurements, and checksums identifying the
plugin and ado files used.

All 120 timed calls and their output-count checks passed. The fastest calls,
such as `decode` at 1 million observations, take only a few milliseconds, so
Stata's timer resolution limits the precision of some speedup estimates.

- [Times and speedups](../benchmarks/readme_profile/results.md)
- [Individual timings and return codes](../benchmarks/readme_profile/timings.csv)
- [Hardware, software and file checksums](../benchmarks/readme_profile/manifest.json)

The full `profile.log` is ignored by Git because Stata's `about` output includes
license details.
