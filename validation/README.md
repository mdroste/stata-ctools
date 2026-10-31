# Run the tests

Run these commands in a terminal from the ctools folder.

## Before you start

1. Install Python 3 and Stata, then [build ctools](../docs/development/PLATFORMS.md).
2. In Stata, install `reghdfe`, `ivreghdfe`, `ppmlhdfe`, `psmatch2`,
   `binscatter`, `binsreg`, `rangestat`, `rangejoin`, `winsor2`, `cipolate`, `distinct` and `gtools`, along with
   any packages they require.
3. On this Mac, the runner invokes the `stata` alias from `~/.zshrc` through
   `/bin/zsh -lic`. Check it with `/bin/zsh -lic 'type stata'`. On other machines,
   configure the launcher for that environment; `CTOOLS_STATA_EXECUTABLE` accepts
   an executable path (required on Windows). Batch drivers end with `exit, clear`.
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

## Check documentation and generated inventories

```sh
python3 scripts/sync_release.py --check
python3 scripts/sync_architecture.py --check
python3 validation/run.py --checks options
python3 validation/run.py --checks docs
```

The docs-tree check invokes Git. On this Mac, run it through zsh outside the
sandbox, as required for other Git operations.

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
