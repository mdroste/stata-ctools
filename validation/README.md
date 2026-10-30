# Validation

This folder contains correctness tests, performance benchmarks and diagnostic
tools for ctools. Stata tests compare commands with their reference implementations;
native tests exercise the C code directly. Repository checks cover packaging,
release metadata, documentation and test infrastructure.

Run the commands below from the repository root. Invoke Python drivers with
`python3 -m validation.<directory>.<module>`.

## Run the tests

List the registered suites and run the native tests and repository checks:

```sh
python3 -m validation.suite_registry --list
python3 -m validation.suite_registry --native
```

For Stata tests, build ctools with `make` and set `CTOOLS_STATA_EXECUTABLE` to
the path of your licensed Stata executable. Install the reference packages
`reghdfe`, `ivreghdfe`, `ppmlhdfe`, `psmatch2`, `binscatter`, `rangestat`,
`rangejoin`, `winsor2` and `gstats` in Stata.

Download the test datasets once, then run the full Stata test suite:

```sh
python3 scripts/fetch_validation_data.py
python3 -m validation.run_stata_audit
```

The audit checks the cached datasets and runs both the master suite and the
additional correctness suites in the registry. A passing run completes every
suite with zero failed or skipped checks. The main log is
`temp/validation/audit_ci.log`; additional do-file suite logs are in
`temp/registry_logs/`. Python drivers use the log paths recorded in `suites.json`.

To check the dataset cache separately, run
`python3 scripts/fetch_validation_data.py --check-only`.
See [test data](fixtures/README.md) for download and cache details.

## Find a test or tool

| Location | Purpose |
| --- | --- |
| `suites.json` | Records each suite's name, path, runner, dependencies and completion marker |
| `suite_registry.py` | Checks the registry, lists suites and runs groups of tests |
| `validate_all.do`, `suites_master.do` | Run the master Stata suite using a generated list of names and paths |
| `run_stata_audit.py` | Runs the full Stata correctness suite and checks its results |
| `stata/` | Stata correctness tests and their Python drivers |
| `native/tests/` | Python drivers that compile and run C tests |
| `native/harnesses/`, `native/memory_safety/` | C tests, including allocation failure and memory safety cases |
| `native/support/` | Compiler helpers and a mock Stata Plugin Interface |
| `checks/` | Checks for builds, packages, releases, documentation and test infrastructure |
| `performance/` | Benchmarks, data preparation, profiling and result summaries |
| `run_performance_audit.py` | Compares performance across builds in separate Stata processes |
| `support/` | Shared Stata and Python helpers, option cases and coverage records |
| `fixtures/` | Dataset manifests with checksums, small test inputs and data builders |
| `probes/` | Manual diagnostic scripts |

For a single check or benchmark, invoke its module directly:

```sh
python3 -m validation.native.tests.test_sort_keys_native
python3 -m validation.checks.check_package build --metadata-only
python3 -m validation.checks.check_docs_tree
python3 -m validation.performance.prepare_command_performance --help
python3 -m validation.run_performance_audit --help
```

Run benchmarks separately from correctness tests. See the
[benchmark methodology](../docs/performance/README.md) for choosing workloads
and comparing results. Save generated data and logs in `temp/`, and investigation
reports in `reports/`. Published benchmark evidence belongs at fixed paths in
`docs/benchmarks/`, with a script that regenerates it.

## Add or move a test

Place suites directly in the appropriate directory and register them in
`suites.json`. Then regenerate the Stata master list and check the registry:

```sh
python3 -m validation.suite_registry --write-master
python3 -m validation.checks.test_validation_registry
```

Discovery covers the validation root, `stata/`, `native/tests/`, `checks/`,
`performance/` and `support/`. It matches `validate_*.do`, `test_*.do`,
`test_*.py`, `benchmark_*.do` and `benchmark_*.py`. Each matching file needs a
registry entry; entries for helpers and manual benchmarks include a reason for
their exclusion from automatic runs. Extend discovery when adding a new suite
directory. Fixture builders, native support files and probes sit outside discovery.

Keep suite names, completion markers, assertions, numerical tolerances,
allocation failure cases and serial/OpenMP variants intact when moving tests.
Put shared code in support modules so each test can run independently. Use
`support/stata_runner.py` to launch Stata; batch drivers close their logs and end
with `exit, clear`. See [native tests](native/README.md) for C test conventions.
