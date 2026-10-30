# Test data

This folder holds dataset manifests, small test inputs and builders for file
format tests. Downloaded datasets are cached under `validation/fixtures/data/`,
which is ignored by Git. Each manifest records source URLs and SHA-256 checksums
so runs use the same data.

## Stata example datasets

`manifest.json` lists 24 official Stata 18 example datasets. Download them once
from the repository root:

```sh
python3 scripts/fetch_validation_data.py
```

Check the cache without downloading:

```sh
python3 scripts/fetch_validation_data.py --check-only
```

The full Stata audit checks this cache before running tests. The `ctools_fixture`
helper in `validation/support/validate_setup.do` loads the cached examples;
Stata's `sysuse` command loads its bundled examples. Missing files or checksum
mismatches fail validation. To update a manifest, review the source changes and
use the fetch script's `--record` option to record the new checksums.

See [running the tests](../README.md#run-the-tests) for Stata setup, reference
packages and the full audit command.

## HDFE benchmark datasets

`hdfe_manifest.json` lists the 15 Stata files in Sergio Correia's HDFE benchmark
collection, about 227 MB in total. Download and check them with:

```sh
python3 scripts/fetch_validation_data.py --collection hdfe
python3 scripts/fetch_validation_data.py --collection hdfe --check-only
```

These files are cached in `validation/fixtures/data/hdfe/` and used by
`validation/performance/run_hdfe_datasets.py`. Download this collection when
running that benchmark. See the [HDFE benchmark guide](../../docs/performance/HDFE_DATASETS.md)
for workloads and options.

## File format inputs

The Python and C builders create small files for import and export tests,
including Excel formats and formulas, malformed files and other format-specific
cases. The checked-in Excel workbooks cover the 1900 and 1904 date systems.
Keep generated files in `temp/`; the builders and small fixed inputs live here.
