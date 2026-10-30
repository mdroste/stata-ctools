# Run timing comparisons

Complete the [test setup](../../validation/README.md#before-you-start) first.
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

## Compare fixed-effects regressions

Download the additional data, then run:

```sh
python3 scripts/fetch_validation_data.py --collection hdfe
python3 validation/bench.py hdfe
```

Check `temp/hdfe_datasets/validation.log` for failures.

## View published results

- [Times and speedups](../benchmarks/readme_profile/results.md)
- [Individual timings](../benchmarks/readme_profile/timings.csv)
- [Computer and software details](../benchmarks/readme_profile/manifest.json)
- [Fixed-effects regression comparisons](HDFE_DATASETS.md)
