# Validation

Run commands from the repository root. Python drivers use package imports and
are invoked with `python3 -m validation.<directory>.<module>`; they need no
installation or custom `PYTHONPATH`. The three root Python entry points also
accept their direct-script invocation for existing local workflows.

## Correctness gates

```sh
python3 -m validation.suite_registry --list
python3 -m validation.suite_registry --native
python3 scripts/fetch_validation_data.py --check-only
python3 -m validation.run_stata_audit
```

The native gate includes C harnesses and repository/runner contract tests. The
complete Stata gate runs the master suites and the extra registered correctness
suites, and rejects incomplete logs, failed assertions and unexpected skips.
Its main log is `temp/validation/audit_ci.log`; extra-suite logs are under
`temp/registry_logs/`. Stata's automatic batch log stays in temporary launch
space and is removed after the process exits.

Stata drivers use the shared launcher in `support/stata_runner.py`. On the
maintainer's macOS machine it invokes the `stata` alias through `/bin/zsh -lic`.
Licensed CI configures `CTOOLS_STATA_EXECUTABLE` for its platform. Batch drivers
must close their logs and end with `exit, clear`.

## Directory map

| Location | Contents |
| --- | --- |
| `suites.json` | Authoritative suite identities, paths, runners, categories, dependencies and completion markers |
| `suite_registry.py`, `suites_master.do` | Registry checks and generated Stata names/paths |
| `validate_all.do`, `run_stata_audit.py` | Master Stata suite and complete correctness gate |
| `run_performance_audit.py` | Entry point for isolated performance comparisons |
| `stata/` | Stata correctness suites, cases and Python orchestration |
| `native/tests/` | Python drivers that compile and run C tests |
| `native/harnesses/`, `native/memory_safety/` | C fixtures, allocation failure and memory safety cases |
| `native/support/` | Shared compiler, harness and SPI mock support |
| `checks/` | Build, package, release, documentation, coverage and runner contracts |
| `performance/` | Benchmarks, data preparation, profiling and summaries |
| `support/` | Shared Stata/Python helpers and option cases/evidence |
| `fixtures/` | Checksum-pinned manifests, small checked-in inputs and fixture builders |
| `probes/` | Manual diagnostics excluded from correctness gates |

See [native test organization](native/README.md) for the white-box harness
contract and [fixtures](fixtures/README.md) for downloads and offline checks.

## Add or move a suite

Register each suite in `suites.json`, then regenerate and check the master list:

```sh
python3 -m validation.suite_registry --write-master
python3 -m validation.checks.test_validation_registry
```

The registry discovers suite filenames directly in `stata/`, `native/tests/`,
`checks/`, `performance/` and `support/`, as well as the root master runner.
Files with `validate_`, `test_` or `benchmark_` prefixes and registered driver
extensions must have an entry. Helpers and manual benchmarks retain explicit
exclusion reasons. Fixture builders, native support and probes are outside
suite discovery. Put new suites directly in the matching directory; extend
discovery explicitly before adding another suite directory.

The generated Stata master list contains both names and paths. Names and
completion markers identify suites independently of filenames. Keep assertions,
numerical tolerances, fault-injection bounds and serial/OpenMP variants intact
when moving tests. Shared code belongs in support modules; tests must not import
executable tests from one another.

## Individual checks and benchmarks

```sh
python3 -m validation.checks.check_package build --metadata-only
python3 -m validation.checks.check_docs_tree
python3 -m validation.native.tests.test_sort_keys_native
python3 -m validation.performance.prepare_command_performance --help
python3 -m validation.run_performance_audit --help
```

Benchmarks are opt-in and are not part of the correctness gate. Keep generated
data, logs and other scratch output under `temp/`. Store investigation reports
under `reports/`, and publish only script-regenerated evidence at fixed paths
under `docs/benchmarks/`. See the [benchmark methodology](../docs/performance/README.md).
