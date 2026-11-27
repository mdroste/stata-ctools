# Run timing comparisons

Complete the [test setup](../../validation/README.md#before-you-start) first.

## Run from Stata

To run the README comparisons and print the timings and speedups, run this in
Stata from the ctools folder. Start with the quick 10,000-observation run:

```stata
do validation/performance/readme_benchmarks.do quick
do validation/performance/readme_benchmarks.do
```

The second command runs the README sizes: 1M, 10M and 25M observations.
Results go to `temp/readme_benchmarks/`. Comparisons whose reference command
(a user-written package) is not installed are skipped and listed, with how to
install it. The top of the do-file lists its options (`sizes()`, `pairs()`,
`processors()`, `output()`).

A comparison is one Stata command and the ctools command that replaces it,
named by the Stata command. A ctools command with several Stata counterparts
has one comparison for each, and each times the ctools command again:

| ctools command | Comparisons |
|---|---|
| `cwinsor` | `winsor` (`gstats winsor`), `winsor2` |
| `cbinscatter` | `binscatter`, `binscatter2`, `binsreg` |
| `creghdfe` | `reghdfe`, `fereg` |
| `cppmlhdfe` | `ppmlhdfe`, `feppml` |

`binsreg` chooses its number of bins from the data unless told, so its
comparison asks for the 20 quantile bins the others use (`nbins(20)`).
`binsreg`, [`binscatter2`](https://github.com/mdroste/stata-binscatter2) and
[`fereg`](https://github.com/johannes-schmieder/fereg) (which includes
`feppml`) are not on SSC. `binscatter2` needs `gtools`.

`qreg`, `ppmlhdfe` and `geonear` take minutes per run at 25M observations, so
those comparisons run at 1M and 10M only. The limits are in `profiles.do`,
`readme_benchmarks.do` and `profile_LIMITS` in `validation/bench_profiles.py`.

Each comparison times only the command. Making the test data, reloading it
before each command and checking each result are not timed but still take
time. To show where a run spends its time, the profiler writes the seconds of
every step to `steps.csv` and prints a `PROFILE_STEP` line after each step.
The report lists every step and ranks the comparisons by total run time. Press
Break to stop a run; the report then covers the steps that finished.

An older ctools release installed its own linear `cipolate.ado`. The profiler
stops if `cipolate` is that copy, and `readme_benchmarks.do` skips the
comparison. Delete the old file and run `ssc install cipolate`.

## Run from a terminal

Run the following commands in a terminal from the ctools folder.

## Start with a small dataset

Run the comparisons with 10,000 observations, then summarize the results:

```sh
python3 validation/bench.py do profile_readme --output temp/timings temp/timings "10000"
python3 validation/bench.py summarize-profile temp/timings --sizes 10000
```

Open `temp/timings/results.md` to read the results. If a comparison fails,
check `temp/timings/profile.log`.

## Run the full comparison

Allow several GB of memory and tens of GB of free disk space. Run:

```sh
python3 validation/bench.py do profile_readme --output temp/timings-full temp/timings-full
python3 validation/bench.py summarize-profile temp/timings-full
```

Open `temp/timings-full/results.md` to read the results.

## Compare distinct counts

The README comparison runs `distinct firm group region sector` against
`cdistinct` on the same 20-variable datasets at 1M, 10M, and 25M rows. It
counts each column separately: firm identifiers, numeric groups, regions,
and string sector codes. Every column's count and the original row order are
checked outside the timed region. The comparison uses one timed call per
command and size, as in the other README profiles.

To run only this pair, with four Stata/MP processors and the default ctools
thread setting:

```sh
python3 validation/bench.py do profile_readme --output temp/timings-distinct temp/timings-distinct "1000000 10000000 25000000" 4 distinct
python3 validation/bench.py summarize-profile temp/timings-distinct --pairs distinct --update-readme README.md
```

`--pairs` limits the summary and README update to the selected comparisons.
The summarizer requires every size a comparison runs at (all three, or 1M and
10M for `qreg`, `ppmlhdfe` and `geonear`) before updating the README. In the README, the
row of a ctools command lists its Stata commands one per line of the first
cell and their speedups in the same order in the last cell. The update changes
only the line of the named command and keeps a footnote marker that follows
the speedup.

## Compare pairwise joins

The README comparison runs `joinby firm` against `cjoinby` on the same
20-variable datasets at 1M, 10M, and 25M rows. The using file holds two
observations for each of the first half of the firms, so the join pairs every
observation of those firms with both and leaves out the other firms; the
result has as many observations as the master and three more variables. The
observation count, the 23 variables and the pairing of every observation with
an event of its firm are checked outside the timed region.

To run only this pair:

```sh
python3 validation/bench.py do profile_readme --output temp/timings-joinby temp/timings-joinby "1000000 10000000 25000000" 4 joinby
python3 validation/bench.py summarize-profile temp/timings-joinby --pairs joinby --update-readme README.md
```

## Compare nearest neighbors and distances

The `geonear` comparison finds, for every observation of the 20-variable
datasets, the nearest of `N`/200 places: the base coordinates are the
profile's normal variables `z` and `noise` (degrees, within a few of 0), and
the places are one per firm, at the coordinates of its first observation. It
runs at 1M and 10M rows (`geonear` takes several minutes at 10M). The
`geodist` comparison computes one ellipsoid distance per observation,
`geodist z noise x2 x3, gen(result)` against `cgeonear`'s distance form with
the same arguments, at all three sizes. The checks outside
the timed region are the observation count, the 22 (or 21) variables and a
nonmissing neighbor and distance (or distance) in every observation.

```sh
python3 validation/bench.py do profile_readme --output temp/timings-geo temp/timings-geo "1000000 10000000 25000000" 4 "geonear geodist"
python3 validation/bench.py summarize-profile temp/timings-geo --pairs geonear geodist --update-readme README.md
```

## Compare fuzzy joins

`reclink`, `matchit` and `strgroup` compare every pair of records or strings in
Stata or Mata (`strgroup` also keeps a matrix of all pairs), so they cannot run
at the sizes of the README profile. `validation/performance/cfuzzyjoin.do`
times each against its `cfuzzyjoin` subcommand on generated firm names with
typos (with a city and a ZIP code for the record linkage), at 1,000, 2,000 and
5,000 observations in each dataset, and `cfuzzyjoin` alone at a larger size.
Outside the timed region it checks that both commands give the same result:
reclink's pairs and scores, matchit's pairs and scores, strgroup's group
numbers. The original commands run once, `cfuzzyjoin` best of ten (three
at the larger size).

```sh
python3 validation/bench.py do cfuzzyjoin_performance --output temp/timings-fuzzy "1000 2000 5000" "50000"
```

The results are on the [cfuzzyjoin page](../commands/cfuzzyjoin.md#performance).

## Compare fixed-effects regressions

Download the additional data, then run:

```sh
python3 scripts/fetch_validation_data.py --collection hdfe
python3 validation/bench.py hdfe
```

Check `temp/hdfe_datasets/validation.log` for failures.

## Compare HDFE ado overhead

Freeze the baseline and candidate ado files in separate directories, with an
identical, protocol-compatible plugin binary in both. From the repository root:

```stata
do validation/support/hdfe_ado_sample.do "path/to/baseline" "path/to/candidate"
do validation/support/hdfe_ado_buffers.do "path/to/baseline" "path/to/candidate"
do validation/support/hdfe_native_screen.do "path/to/baseline" "path/to/candidate"
do validation/support/hdfe_native_expanded.do "path/to/baseline" "path/to/candidate"
```

The sample check compares coefficients, covariance matrices, coefficient names,
observation counts, estimation samples and error codes, including missing values,
weights, singleton groups and expressions using previous estimates. It also
checks PPML's rejection of negative outcomes before singleton removal.
The buffer check covers models without absorbed effects, saved residuals and
fixed effects, excluded rows, output storage types and formats, `by:`, temporary
variable cleanup, and cold/warm clock calls preserving `r()`.
The IV comparisons also check `e(N_full)`: unweighted fits obtain this count
from native row and singleton counts when the ado skips its pre-fit count.

The benchmark generates a 10-million-row balanced panel, warms each command
including its display routines, and prints five verbose timing reports for each
estimator. Run these in separate Stata processes with identical settings:

```stata
do validation/performance/hdfe_ado_overhead.do "path/to/baseline" baseline
do validation/performance/hdfe_ado_overhead.do "path/to/candidate" candidate
```

Optional arguments select the number of individuals (ten rows each), repetitions
and commands, for example `1000000 5 "creghdfe"`. Repeat the process pair in
reverse order to check for changes in machine load.

Keep the log under `temp/`. Compare medians of the tagged `ADO_BENCH_BEGIN` /
`ADO_BENCH_END` blocks. The reported OLS/IV "Plugin call overhead" includes native
result publication outside the plugin's internal timer; it is not entirely ado
execution. The three estimators hand their working memory to a helper thread
(`src/ctools_release.h`), which returns it while Stata posts the results; the
next ctools command waits for it before it starts. Use a separate Stata process for each directory so that
plugin and program state cannot carry over from the other build.

Plain numeric, unweighted models now defer missing-value screening to the
plugin. OLS and PPML scan the loaded numeric inputs before FE remapping; IV uses
its existing parallel screening. OLS reports `Missing screening` separately
from SPI loading. Large samples use the command's default thread budget,
including when the model has fewer variables than available threads;
`threads(1)` selects serial execution.

Weighted and factor-variable models use an earlier native streaming pass so
the marked sample is ready for factor expansion and weight-expression staging.
Stata still evaluates and validates weights and excludes zero/missing weights;
native code screens the model inputs. This preserves invalid-weight errors,
weight expressions using previous estimates, and the observations used to
choose factor levels and base categories. IV retains its original preparation
boundaries and its explicit row count for weighted `e(N_full)`.

The early pass uses bounded per-thread buffers and updates only the private
sample marker. It does not load a second dataset copy or allocate a full row
mask. Numeric and extended missing values are screened in parallel, including
in base/omitted factor inputs. String/interacted/sloped absorbed effects and
PPML offsets/exposures keep their own preparation checks after the model scan.
The early screening cost is included in the reported Stata setup time.

The OLS/PPML screen uses vectorized scans of 4,096-row blocks. An all-valid
sample needs only one flag per block and no compaction. When rows are missing,
the native code clears their private sample markers and compacts columns and
the observation map in stable row order. Column compaction is parallel across
columns: parallel in-place writes across row blocks could overwrite another
worker's input. Numeric missing codes, extended missing values and NaNs are
excluded before singleton/separation decisions.

Time-series terms, bootstrap, individual/group models and auxiliary IV options
that require earlier sample preparation retain the existing screening path.
Separate plugin capabilities preserve compatibility with older plugins, including
ones that support only the scan after loading. Rebuild the plugin and restart
Stata to enable the extension. The native sample checks compare automatic and
one-thread runs, missingness, metadata, weight validation, factor coefficient
names and empty-sample errors against old ado files with the **same new plugin**.

For setup profiling, separate sample marking from option parsing. `profiler
on`/`off` and `profiler report` identify ado helper costs; temporary
`_ctools_lap` calls around `marksample` separate its cost from other built-ins.

Output initialization also matters at large N. The wrappers use `mark` for
constant byte variables and Mata's `st_addvar()` for selected missing-valued
outputs. The latter retains its default missing initialization: skipping that
fill would leave excluded rows undefined. Clock laps reuse the registered clock
plugin; normal estimator loading still checks plugin identity and protocol.

## View published results

Published timings describe the builds and settings recorded with each run.
They are not a fresh measurement of the current checkout. The README profile
measures linear `cinterpolate` against `ipolate`, `cinterpolate, cubic`
against SSC `cipolate`, and `chdfe` against SSC `hdfe`. A published time is
the smallest over the passes the manifest records (`repetitions`); a
comparison it lists under `comparisons_not_measured` has no row. A comparison
refreshed on its own after the shared run (see `command_runs` in the manifest)
was measured on the build recorded there. `qreg` and `ppmlhdfe` were timed at
1M and 10M observations only. `csort` and `cmerge` do best on data with many
observations relative to variables; their speedups shrink as the number of
variables grows.

- [Times and speedups](../benchmarks/readme_profile/results.md)
- [Individual timings](../benchmarks/readme_profile/timings.csv)
- [Computer and software details](../benchmarks/readme_profile/manifest.json)
- [Fixed-effects regression comparisons](HDFE_DATASETS.md)

To refresh the published timing results, collect and summarize a complete run
under `temp/`, then copy only `timings.csv` and `results.md` into
`docs/benchmarks/readme_profile/`. Update `manifest.json` from the actual run's
machine, Stata and reference-package versions, thread settings, build versions,
and source revisions. The profiler writes the CSV and the summarizer writes the
table; neither creates that
run-information file. Keep logs, driver do-files, generated data, and exploratory
results under `reports/performance/` or `temp/`.

For a single-command refresh, replace that command's rows in the published
CSV and regenerate the combined table. Record its machine, settings, checks,
and build/reference hashes in the manifest's `command_runs` entry; those values
override the shared manifest fields for that command. Preserve the measured
rows and run information for commands that were not rerun.
