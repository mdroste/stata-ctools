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
Results go to `temp/readme_benchmarks/`. Comparisons whose SSC reference
command is not installed are skipped and listed. The top of the do-file lists
its options (`sizes()`, `pairs()`, `processors()`, `output()`).

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
The summarizer requires all three sizes before updating the README.

## Compare fixed-effects regressions

Download the additional data, then run:

```sh
python3 scripts/fetch_validation_data.py --collection hdfe
python3 validation/bench.py hdfe
```

Check `temp/hdfe_datasets/validation.log` for failures.

## View published results

Published timings describe the builds and settings recorded with each run.
They are not a fresh measurement of the current checkout. The README profile
measures linear `cinterpolate` against `ipolate`, `cinterpolate, cubic`
against SSC `cipolate`, and `chdfe` against SSC `hdfe`. The published results
do not yet include the `cipolate` and `hdfe` comparisons.

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

For a single-command refresh, replace that command's six rows in the published
CSV and regenerate the combined table. Record its machine, settings, checks,
and build/reference hashes in the manifest's `command_runs` entry; those values
override the shared manifest fields for that command. Preserve the measured
rows and run information for commands that were not rerun.
