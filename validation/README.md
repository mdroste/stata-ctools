# Run the tests

Run these commands in a terminal from the ctools folder.

## Before you start

1. Install Python 3 and Stata, then [build ctools](../docs/development/PLATFORMS.md).
2. In Stata, install `reghdfe`, `ivreghdfe`, `ppmlhdfe`, `psmatch2`,
   `binscatter`, `binsreg`, `rangestat`, `rangejoin`, `winsor2`, `cipolate`, `hdfe`, `distinct` and `gtools`, along with
   any packages they require.
3. On this Mac, the runner invokes the `stata` alias from `~/.zshrc` through
   `/bin/zsh -lic`. Check it with `/bin/zsh -lic 'type stata'`. On other machines,
   configure the launcher for that environment; `CTOOLS_STATA_EXECUTABLE` accepts
   an executable path (required on Windows). Batch drivers end with `exit, clear`.
4. Download the test data:

```sh
python3 scripts/fetch_validation_data.py
```

Install `binsreg` from its upstream Stata package (it is not on SSC):

```stata
net install binsreg, from(https://raw.githubusercontent.com/nppackages/binsreg/main/stata) replace
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
python3 validation/run.py --checks seeds
```

The runner executes many named cases in one Stata session, so each case (and
each benchmark workload) must `set seed` before its first random draw; a seed
set by another case does not count. The seeds check enforces this.

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

## Check that the tests reject incorrect results

```sh
python3 validation/run.py --stata comparison_guards coverage_helpers
python3 validation/run.py --native test_harness_guards_native
python3 -c "from validation.check_runner import test_runner; test_runner()"
```

These negative controls deliberately inject incorrect values, missing codes,
dropped columns, reordered rows, stale estimates and graphs, and differences
after a long identical file prefix. Correct controls must pass and each mutation
must be rejected. Matching nonzero native/ctools return codes remain valid
passing error-compatibility tests; different codes or success versus error fail.

The runner rejects incomplete results, recorded failures and skipped checks.
Python optimization is refused because it disables assertions. Native test
compilation keeps C assertions enabled and makes sanitizer findings fatal,
including when caller flags or environment settings would disable those checks.

For timing comparisons, follow the [benchmark instructions](../docs/performance/README.md).

## HDFE cancellation

`python3 validation/run.py --native test_hdfe_cancel_native` runs deterministic
host-poll cancellation tests with and without OpenMP. It checks weighted and
unweighted solver families, nested teams, a host that finishes before its workers,
bootstrap abort/cleanup, and successful reuse after cancellation. Only the calling
thread may invoke the simulated host callback.

For licensed rollback checks, generate a test entry point with
`python3 validation/support/hdfe_cancel_plugin.py --output temp/cancel/ctools_plugin.c`.
Build a separate plugin using that source instead of `src/ctools_plugin.c`, keeping
the normal compiler flags and all other sources. Set `CTOOLS_TEST_PLUGIN_DIR` to
its directory and `CTOOLS_TEST_CANCEL_PLUGIN=1`, then run
`validation/support/hdfe_cancel.do`. This exercises the real dispatch, cleanup and
ado wrappers with deterministic Break responses; it does not measure interactive
keyboard/UI Break latency. Never put the instrumented plugin in `build/` or a
release package.

`python3 validation/run.py --native test_newcommands_paths_native` checks grouped
linear/cubic interpolation, split token bytes, ordered range joins, allocation
and publication failures, and the join output-size boundary. It runs serial and
OpenMP variants with independent expected results. For the corresponding licensed
rollback fixture, build a separate **serial** plugin with
`validation/support/publication_fault.h` included after the normal headers in
`cinterpolate_impl.c`, `csplit_impl.c` and `crangejoin_impl.c`; set
`CTOOLS_TEST_PUBLICATION_PLUGIN=1` and point `CTOOLS_TEST_PLUGIN_DIR` at that
plugin before running `validation/support/newcommands_publication.do`. This
instrumentation is test-only and must not be included in release binaries.

`python3 validation/run.py --native test_destring_paths_native` exercises the
production converter and numeric parser with independent expected values,
serial/OpenMP execution, and allocation/store/result-publication failures.
For licensed staged-output checks, include the same test-only publication header
in `cdestring_impl.c` in a separate serial plugin and run
`validation/support/destring_publication.do` with the two environment variables
above. It verifies every scalar publication in generate/replace and empty-sample
paths, preserves original data/labels, and checks subsequent successful reuse.

`test_bins_paths_native` covers external bin labels beyond the signed-int range,
weighted median ties, sampled quantiles, few-value polynomials, collinear controls,
absorbed effects, allocation failures, and parallel block reductions.
`test_hdfe_paths_native` checks analytically constructed weighted/unweighted slope
residuals, orthogonality, idempotence, complete native metadata, nesting,
disconnected groups, singleton removal, extreme scaling, and pool/thread counts.
Both run through `python3 validation/run.py --native <test_name>` with serial and
OpenMP sanitizer variants. The licensed `chdfe` suite additionally compares every
per-effect DoF/exact/nested/owner field with `hdfe`.

`test_labels_paths_native` exercises the production encoder and decoder, including
200,001-row worker paths, existing-label mappings, malformed protocols and
allocation/transfer/publication failures. The label-file writer itself retains
its separate framing tests; this fixture checks the encoder's writer inputs and
failure return. For licensed wrapper checks, include `publication_fault.h` in
`cencode_impl.c` and `cdecode_impl.c` in a separate serial plugin and run
`validation/support/labels_publication.do` with the fault-plugin environment
variables above. It checks encoding rollback and decoding's complete Stata
fallback after each result-save failure.

`test_parquet_paths_native` reads complete Parquet files from the repository
writer and checked-in independent producers. It checks exact cells, row ranges,
codecs, v1/v2 pages, dictionaries/nulls, malformed metadata and page bodies,
allocation faults, cleanup and a bounded byte-mutation/truncation corpus.
The Linux runner selects ASan/UBSan; `CTOOLS_SANITIZERS` can select an explicitly
supported sanitizer combination. A compiler accepting a flag alone is not
evidence that its sanitizer runtime is active.

For licensed malformed-file rollback checks, generate a persistent corpus:

```sh
python3 -c "from pathlib import Path; from validation.support.parquet_native_fixtures import corpus; corpus(Path('temp/parquet-corpus'))"
```

Set `CTOOLS_PARQUET_CORPUS` to that directory's absolute path and
`CTOOLS_TEST_PLUGIN_DIR` to the candidate plugin directory, then run
`validation/support/parquet_malformed.do`. It checks exact error codes and
unchanged data in default/snapshot modes; no fault-instrumented plugin is needed.
