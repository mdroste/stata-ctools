# Run the tests

Run these commands in a terminal from the ctools folder.

## Before you start

1. Install Python 3 and Stata, then [build ctools](../docs/development/PLATFORMS.md).
2. In Stata, install `reghdfe`, `ivreghdfe`, `ppmlhdfe`, `psmatch2`,
   `binscatter`, `binsreg`, `rangestat`, `rangejoin`, `winsor2`, `cipolate`, `hdfe`, `distinct` and `gtools`, along with
   any packages they require. [docs/ATTRIBUTION.md](../docs/ATTRIBUTION.md)
   lists the versions ctools is pinned to; `which <command>` shows the version
   you have. If you update ctools to a newer version of an original command,
   record the new version in `config/attribution.json` and run
   `python3 scripts/sync_attribution.py`.
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

Install `distinct` from the *Stata Journal* archive. SSC distributes version
1.2.1; `cdistinct` follows 1.4.3, and its suite stops if another version comes
first on the adopath. If the SSC copy is installed, remove it first with
`ado uninstall distinct`.

```stata
net install dm0042_5, from(http://www.stata-journal.com/software/sj23-4) replace
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
python3 scripts/sync_attribution.py --check
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

A helper that compares what two commands print (`_ct_same_output` and similar)
reads the output back from a log. Stata writes nothing to a log while output
is quiet, and `noisily` inside a called program does not undo a `quietly` in
its caller, so a helper called from a suite's `quietly {` block would compare
two empty texts and pass. These helpers call `_ct_loud` (in
`validation/support/validate_setup.do`) first, which stops a quiet call with
r(198). Prefix the call with `noisily` at the level of the `quietly` block, and
call `_ct_loud <name>` at the top of any new helper that captures output.

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

`python3 validation/run.py --native test_hdfe_rounding_native` checks the
row-space conjugate gradient at the rounding of its input: columns with a level
of 40, 2000 and 1e6 against their centered copy at tolerances from 1e-12 to
1e-16, a column the effects explain, a column that is already projected, a
random and a balanced design, with and without weights. It also requires the
plain Hestenes-Stiefel iteration, bit for bit, where the rounding is not
reached, the same bits for one and four threads, and r(430) when one iteration
is not enough.

`python3 validation/run.py --native test_hdfe_grouped_native` checks the
one-pass projections on rows that arrive grouped by a fixed effect, bit for
bit against the general passes written out in the test: one effect whose
levels never decrease (short runs, runs longer than a block, more columns than
are walked together, weights, a lone decrease after the probed rows) against
the scatter and gather, and the fused right-hand side of the coefficient-space
solve (levels without rows, a level of 30,000 rows, one to three kept effects,
every arrangement of the row plan) against the two separate passes.

`python3 validation/run.py --native test_creghdfe_columns_native` checks
creghdfe's column preparation (means, total sums of squares, standard
deviations, norms and the standardized matrix, written once from the loaded
variables) bit for bit against the copy, sum, standardize and dot-product
passes it replaces: dropped rows, factor-variable indicators, each weight type,
the double-double sums, the uncentered total, and matrices on either side of
the size from which the stores bypass the cache, at both alignments.

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

`python3 validation/run.py --native test_six_release_native` checks malformed
range-bound source indices, every range-statistic scalar publication, and
bootstrap id-pass allocation/read/store/Break failures in serial and OpenMP
builds. The six command suites also run
`validation/support/six_release_regressions.do`: quoted ignore characters,
scalar-name collisions, excluded outliers, exact and overflowing interval
bounds, cluster-sorted size expressions, and bootstrap rollback after a real
maxvar failure or injected drop/expand/id-pass/Break errors. The destring native
fixture additionally checks long-decimal allocation failure and aggregate
counts above `INT_MAX` without allocating billion-row datasets.

`test_bins_paths_native` covers external bin labels beyond the signed-int range,
weighted median ties, sampled quantiles, few-value polynomials, collinear controls,
absorbed effects, allocation failures, and parallel block reductions.
`test_hdfe_paths_native` checks analytically constructed weighted/unweighted slope
residuals, orthogonality, idempotence, complete native metadata, nesting,
disconnected groups, singleton removal, extreme scaling, and pool/thread counts.
`test_hdfe_coef_native` holds every method of the coefficient-space projection
(dense, conjugate gradients, Laplacian preconditioner; two to four effects) to a
Gram-Schmidt orthogonalization of the dummies on small designs: repeated cells,
disconnected and nested effects, zero weights, more columns than one batch, and
level means a million times the noise. It also checks the dispatch of
`partial_out_columns()`, the iteration limit, Break at every poll, every
allocation failure, and, with OpenMP, that 300,000-row designs give the same
bits on one thread and on four.
`test_cross_products_native` holds the estimators' cross products
(`ctools_matmul_ata/atb/atdb`), the compensated ones (`ctools_kahan_cross`,
`kahan_dot`, the double-double sums) and the robust and clustered meats to
sequential references: bit for bit wherever the definition is the sequential
sum, within rounding of a double-double sum otherwise. It runs each kernel
again with its scratch allocation refused (the other loop arrangement must
give the same bits), on one thread against four, across block, chunk and wave
edges, with every weight type and with clusters random, sorted and empty, and
builds the clustered meat both ways (`CTOOLS_VCE_CLUSTER_SUMS`).
All three run through `python3 validation/run.py --native <test_name>` with serial and
OpenMP sanitizer variants. The licensed `chdfe` suite additionally compares every
per-effect DoF/exact/nested/owner field with `hdfe`.

`test_hdfe_sample_parallel_native` runs the HDFE sample construction in parallel
row blocks and serially, on dense, negative, signed-zero, non-integer, sparse,
wide-range, single-level and unique codes: the loaded levels must define the same
partition, singleton removal (with and without fweights) must keep the rows a
round-by-round recount keeps, PPML sample selection must remove and classify the
rows a recount every round does, the renumbered factors (levels, counts and
weighted counts) must match bit for bit with and without dropped rows, connected
components must match a rank-based union-find on disconnected designs, and the
nesting check must agree. It also checks that every path runs serially inside a
parallel region.

`test_hdfe_sample_ordered_native` runs the sample construction whose steps share
their work (`hdfe_load_factors_ordered`, `hdfe_drop_singletons`,
`hdfe_build_factors_ordered`, the per-level counts all at once or job by job)
against an independent reference and against the separate steps, on one, three
and eight threads: the kept rows, the levels in order of first appearance among
them, the rows per level and the weights per level summed in row order must match
bit for bit, and the load's first rows and single-row levels must be right. Its
thresholds are small, so a few hundred rows reach the row blocks, the singleton
step's counting pass and its fallback to counting every level, the in-place shift
of the levels over dropped rows, and both parallel counts; a second build uses
other thresholds. It also drops arbitrary rows through keep flags
(`hdfe_singletons_from_keep`, listed and not) and compares the factors with
`hdfe_build_factors` on the same flags, and checks that
`ppml_fe_sample_is_whole` answers yes exactly when the PPML sample selection
removes no row.

`test_hdfe_missing_native` also screens columns that carry the load's note of
missing values (`stata_variable.missing_state`): the result is the same, and data
noted free of them costs no allocation.

`test_hdfe_numerics_native` checks weak-connection breakdown, healthy small CG
steps, fully absorbed columns, and invariance to column/weight units. The
`creghdfe` suite also runs `support/creghdfe_numerics.do`: six-row weak-FE and
eight-row near-collinearity/scaling regressions with analytical coefficient
oracles, and explicit rejection of unsupported solver requests. Performance
comparisons for the affected paths are in `performance/creghdfe_numerics.do`;
select each build with `CTOOLS_TEST_PLUGIN_DIR` in a fresh Stata process.

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
