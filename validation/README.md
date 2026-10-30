# Run the tests

Run these commands in a terminal from the ctools folder.

## Before you start

1. Install Python 3 and Stata, then [build ctools](../docs/development/PLATFORMS.md).
2. In Stata, install `reghdfe`, `ivreghdfe`, `ppmlhdfe`, `psmatch2`,
   `binscatter`, `rangestat`, `rangejoin`, `winsor2` and `gtools`, along with
   any packages they require.
3. Make sure the `stata` command works in your terminal. Otherwise, set
   `CTOOLS_STATA_EXECUTABLE` to the location of your Stata application.
4. Download the test data:

```sh
python3 scripts/fetch_validation_data.py
```

## Run all tests

```sh
python3 validation/run.py
```

To run only the Stata tests:

```sh
python3 validation/run.py --stata
```

To run checks that do not use Stata:

```sh
python3 validation/run.py --native
```

## Check one command

Replace `cencode` with the command you want to check:

```sh
python3 validation/run.py cencode
```

To list the available tests:

```sh
python3 validation/run.py --list
```

Read the pass/fail messages in the terminal. If a test fails, open the log
at the printed path under `temp/validation/`.

For timing comparisons, follow the [benchmark instructions](../docs/performance/README.md).
