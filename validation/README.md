# Tests and benchmarks

This folder contains tests for ctools and scripts for measuring its speed.
Stata tests compare ctools with the commands it replaces. C tests check individual
functions. Other checks cover installation files, release versions, documentation
and the scripts that run the tests.

Run the commands below from the repository root.

## Run the tests

List the tests, then run the C tests and repository checks:

```sh
python3 -m validation.suite_registry --list
python3 -m validation.suite_registry --native
```

For Stata tests, build ctools with `make` and set `CTOOLS_STATA_EXECUTABLE` to
the path of your licensed Stata executable. Install the commands used for
comparison: `reghdfe`, `ivreghdfe`, `ppmlhdfe`, `psmatch2`, `binscatter`,
`rangestat`, `rangejoin`, `winsor2` and `gstats`, along with their dependencies.

Download the test datasets once, then run all Stata correctness tests:

```sh
python3 scripts/fetch_validation_data.py
python3 -m validation.run_stata_audit
```

The script checks the downloaded datasets and runs the tests listed in
`suites.json`. A passing run finishes every test group with zero failed or
skipped checks. Read `temp/validation/audit_ci.log` for the main results.
Additional do-file logs are in `temp/registry_logs/`; the log paths for Python
scripts are recorded in `suites.json`.

To check the downloaded files separately, run
`python3 scripts/fetch_validation_data.py --check-only`.
See [test data](fixtures/README.md) for details.

## Find a test or tool

| Location | Contents |
| --- | --- |
| `suites.json` | Test names, paths, required packages, how to run each test and how to check that it finished |
| `suite_registry.py` | Reads and checks that list, and runs groups of tests |
| `validate_all.do`, `suites_master.do` | Main Stata test script and its generated list of tests |
| `run_stata_audit.py` | Runs all Stata correctness tests and checks the logs |
| `stata/` | Stata tests and Python scripts that run them |
| `native/tests/` | Python scripts that compile and run C tests |
| `native/harnesses/`, `native/memory_safety/` | C test programs, including tests of memory handling and allocation failures |
| `native/support/` | Compiler helpers and a test substitute for the Stata Plugin Interface |
| `checks/` | Checks of builds, installation files, release versions, documentation and test scripts |
| `performance/` | Benchmarks and scripts to prepare data and summarize timings |
| `run_performance_audit.py` | Compares builds in separate Stata processes |
| `support/` | Shared test functions, option tests and records of which options are tested |
| `fixtures/` | Test data, download lists and scripts that create test files |
| `probes/` | Scripts for investigating specific problems manually |

Run individual Python checks as modules:

```sh
python3 -m validation.native.tests.test_sort_keys_native
python3 -m validation.checks.check_package build --metadata-only
python3 -m validation.checks.check_docs_tree
```

Benchmarks run separately. These commands describe the available options:

```sh
python3 -m validation.performance.prepare_command_performance --help
python3 -m validation.run_performance_audit --help
```

See [measuring performance](../docs/performance/README.md) for workloads and
how to interpret the results. Save generated data and logs in `temp/`, and
investigation reports in `reports/`. Keep published benchmark results in
`docs/benchmarks/`, with a script that updates them at the same paths.

## Add or move a test

Put the test in the appropriate directory and add it to `suites.json`.
Then regenerate the list used by Stata and check that every test is listed:

```sh
python3 -m validation.suite_registry --write-master
python3 -m validation.checks.test_validation_registry
```

The checker looks directly in the validation root, `stata/`, `native/tests/`,
`checks/`, `performance/` and `support/` for `validate_*.do`, `test_*.do`,
`test_*.py`, `benchmark_*.do` and `benchmark_*.py`. Each matching file needs
an entry. For helpers and manual benchmarks, explain in the entry why they run
separately. Update the checker when adding another test directory. It does not
search `fixtures/`, `native/support/` or `probes/`.

When moving tests, preserve their names, checks, numerical tolerances, simulated
allocation failures, and runs with and without OpenMP. Keep the messages that
confirm each test finished. Put shared functions in support modules so tests
can run independently. Use `support/stata_runner.py` to launch Stata; batch
scripts close their logs and end with `exit, clear`. See
[C tests](native/README.md) for more detail.
