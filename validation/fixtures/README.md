# Test data

This folder contains small test files, scripts that create test data, and lists
of datasets to download. Downloaded files go in `validation/fixtures/data/`,
which is ignored by Git. The download lists record each file's URL and SHA-256
checksum so the tests can check that its contents have not changed.

## Stata example datasets

`manifest.json` lists 24 official Stata 18 example datasets. Download them once
from the repository root:

```sh
python3 scripts/fetch_validation_data.py
```

Check the downloaded files without fetching anything:

```sh
python3 scripts/fetch_validation_data.py --check-only
```

The full Stata test script checks these files before running tests. The
`ctools_fixture` helper in `validation/support/validate_setup.do` loads them;
Stata's `sysuse` command loads its bundled examples. Missing or changed files
cause the checks to fail. When intentionally updating the datasets, review the
source changes and use the fetch script's `--record` option to save new checksums.

See [running the tests](../README.md#run-the-tests) for Stata setup and required
comparison commands.

## HDFE benchmark datasets

`hdfe_manifest.json` lists the 15 Stata files in Sergio Correia's HDFE benchmark
collection, about 227 MB in total. Download and check them with:

```sh
python3 scripts/fetch_validation_data.py --collection hdfe
python3 scripts/fetch_validation_data.py --collection hdfe --check-only
```

These files go in `validation/fixtures/data/hdfe/`. Download them when running
`validation/performance/run_hdfe_datasets.py`, which compares fixed-effects
regression times. See the [HDFE benchmark guide](../../docs/performance/HDFE_DATASETS.md)
for datasets and options.

## File format tests

Python and C programs in this folder create small files for import and export
tests, including Excel formats and formulas, malformed files and other cases
specific to each format. The Excel workbooks stored here cover the 1900 and
1904 date systems. Save generated files in `temp/`.
