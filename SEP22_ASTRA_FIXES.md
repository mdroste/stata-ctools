# September 22 audit repairs

This document records the implementation work requested after
[SEP22_ASTRA_AUDIT.md](/Users/Mike/Documents/GitHub/stata-ctools/SEP22_ASTRA_AUDIT.md).
All 16 findings are addressed in the working tree. The original audit is retained
as historical evidence; its pre-repair assessment does not describe this build.
Existing uncommitted work was preserved, and `src/stplugin.c` and `src/stplugin.h`
were not changed. No commit, remote push, or release publication was performed.

The complete Stata gate passes **2,866 checks, with zero failures and zero
unexpected skips**, across all 20 components. It separately records 48 matching
standard-error comparisons as documented method differences, not successful
comparisons. The repaired package has been installed and exercised locally on
Apple Silicon. The other three platform builds and licensed CI publication have
not been executed from this session.

## Resolution by finding

| ID | Original priority | Resolution |
|---|---|---|
| S01 | P1 | Enforce merge cardinality before mutation and restore the master on failure. |
| S02 | P1 | Use plugin-local winsor output indices; roll back rejected operations. |
| S03 | P1 | Preserve full numeric cluster identities in OLS and PPML. |
| S04 | P1 | Map IV numeric clusters before integer conversion; reject insufficient retained clusters. |
| S05 | P1 | Post the actual retained OLS/PPML estimation sample. |
| S06 | P1 | Select the PPML sample before compaction and apply one original-row mask to every row identity. |
| S07 | P1 | Decode through Stata's native label engine with transactional multi-variable output. |
| S08 | P1 | Honor observation/cluster draw counts separately in each stratum. |
| S09 | P2 | Register an explicit PPML prediction handler and reject undefined predictions. |
| S10 | P2 | Evaluate numeric weight expressions once before calling the estimators. |
| S11 | P2 | Store actual matched-neighbor counts. |
| S12 | P2 | Apply empty-side merge filters even with `nogenerate`. |
| S13 | P2 | Stage a complete package with an explicit platform scope and file hashes. |
| S14 | P2 | Make the full gate offline, require component completion, and block publication without licensed validation. |
| S15 | P2 | Correct the installed syntax, output, and unsupported-option documentation. |
| S16 | P3 | Update build/developer instructions and check command inventories. |

### S01 and S12: merge cardinality, rollback, and empty inputs

[cmerge.ado](/Users/Mike/Documents/GitHub/stata-ctools/build/cmerge.ado) checks
uniqueness on each side declared unique. This includes unmatched keys and missing
key values. Invalid `1:1`, `m:1`, or `1:m` requests return `r(459)` instead of
discarding duplicate rows. The join kernel independently checks these contracts
in [cmerge_join.c](/Users/Mike/Documents/GitHub/stata-ctools/src/cmerge/cmerge_join.c).
The wrapper preserves the master and restores it after errors, including a
failure while processing the using frame.

Empty-input branches now filter using their known result code. Thus omitting the
visible merge indicator cannot disable `keep(match)` or another `keep()` choice.
The regression cases verify rejection with unchanged data signatures, duplicate
unmatched keys, missing string keys, and empty-side filtering. The complete merge
component passes 110 checks.

### S02: winsor output and sample boundaries

[cwinsor.ado](/Users/Mike/Documents/GitHub/stata-ctools/build/cwinsor.ado) uses
the positions in the actual plugin varlist, preserving the checked store API.
The surrounding transaction removes incomplete generated variables and restores
inputs after an error. Generated results are missing outside `if/in`; `replace`
leaves excluded input observations unchanged.

Checks cover targets moved between dataset columns, generated and replacement
outputs, grouping, qualifiers, and a later invalid output name after an earlier
valid one. All 55 component checks pass. The destination-index correction already
present at the start of the repair task was preserved and verified.

### S03 and S04: cluster identity and valid covariance inputs

The OLS/PPML shortcut now compares original doubles with FE values instead of
comparing integer casts. General IV cluster mapping uses the new
`ctools_numeric_to_cluster_ids()` helper in
[ctools_hdfe_utils.c](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_hdfe_utils.c).
It sorts full double values before creating dense integer IDs, treating missing
values separately from legitimate negative labels. Both IV cluster dimensions
use this path.

The regression cases check covariance invariance under fractional, negative,
large, and string relabelings, including two-way IV clustering. The original
OLS/PPML design now retains all 40 clusters. IV rejects fewer than two retained
clusters in either requested cluster dimension with `r(459)`; PPML also rejects
an insufficient cluster count before computing its covariance. PPML cluster
allocation/remapping failures return an error rather than using incomplete IDs.

### S05 and S06: estimation samples and PPML separation

OLS and PPML now receive a separate sample-output variable in their plugin
varlists. Their ado wrappers post that returned indicator as `e(sample)`. OLS
keeps an original-observation map through singleton removal. PPML's
`ppml_select_fe_sample()` removes singletons and all-zero FE groups to a fixed
point using an original-length mask, before data compaction. Data, weights,
offsets, cluster values, and observation identities then use that same selection.
The shortened-mask/original-length read is removed.

Direct SPI writes in these paths use `SD_SAFEMODE`, so plugin-local output
indices cannot be interpreted as raw dataset indices. Residual and group-output
store errors are checked. The PPML IRLS-weight ownership was also corrected to
release that allocation during cleanup.

Checks verify sample membership, saved residuals, stored/restored estimates,
singleton exclusion, qualifiers, and equivalent fits after manually deleting
separated groups at the beginning, middle, and end of the data. The PPML cases
combine separation with a singleton, fractional clusters, analytic weights, and
an offset. AddressSanitizer/UndefinedBehaviorSanitizer checks cover the selection
and numeric-cluster helpers, exact-sized masks, singleton chains, and every
selection-allocation failure point. The complete Stata plugin was tested normally,
not inside an ASan-instrumented Stata process.

General separation beyond FE removal remains unsupported. The old IRLS branch
repeatedly counted retained observations as dropped and changed their weights to
tiny values. It now returns `r(430)` with a specific diagnostic when its
`septolerance()` guard triggers. It does not post a misleading successful fit.
The installed help states this boundary. A regression case checks both rejection
and a successful subsequent fit. With frequency weights, `e(N)` is the weight sum;
it need not equal the unweighted count of `e(sample)`.

### S07: literal and long value labels

[cdecode.ado](/Users/Mike/Documents/GitHub/stata-ctools/build/cdecode.ado) stages
native `decode` results for every input before renaming or replacing variables.
This avoids parsing executable `label save` text and preserves quoted text,
literal macro punctuation, tabs/newlines, Unicode, long strings, and labeled
extended missing values. Native `maxlength()` behavior is retained. The public
`r(N_vars)` result is preserved.

The legacy C entry points now reject outdated wrappers with a useful diagnostic.
`threads()` remains accepted for compatibility and does not control native
decoding; `verbose` identifies the native engine. These semantics are documented.
All 219 existing decode checks pass, and the new regression case verifies exact
label text, a 3,000-character label, native truncation, destructive replacement,
and rollback when a later input lacks a value label.

### S08: bootstrap sampling units

[cbsample_impl.c](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c)
distinguishes an omitted count from an explicit positional count. The default
draws each stratum's original number of sampling units; an explicit count draws
that many observations or clusters within each stratum. Stratum identity is part
of cluster identity, so repeated cluster labels in different strata remain
separate units. Oversized requests fail before any frequency output is written.

The wrapper validates an existing numeric `weight()` destination and excludes
missing strata consistently with native `bsample`. Tests check frequency totals,
expanded row counts, cluster draws, unequal strata, reused cluster labels, and
invalid sizes. The full bootstrap component passes 29 checks.

### S09 and S10: prediction and weight expressions

The new installed
[cpplmhdfe_p.ado](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe_p.ado)
is registered in `e(predict)` and included in the package manifest. Explicit
`predict ..., xb` returns the slope index, excluding absorbed effects and the
offset. Default prediction and unsupported fitted means return `r(198)` rather
than silently returning a quantity with the wrong meaning. Stored/restored
estimates retain the handler.

OLS, IV, and PPML evaluate a supplied numeric weight expression once into a
temporary double variable and pass its name to the plugin. Tests compare an
expression such as `[aw=2*w]` with an explicitly generated weight variable.

### S11: actual matching counts

[cpsmatch_impl.c](/Users/Mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_impl.c)
returns actual neighbor counts through an explicit output column. Nearest-neighbor
counts include selected ties, radius counts include eligible controls, kernel
counts include positive-weight controls, and the supported no-replacement case
reports one selected control. An eligible unmatched treated observation receives
zero; controls and observations outside the estimation sample are missing.
The wrapper no longer overwrites these counts with the requested `neighbor()`.
Targeted cases cover ties, too few controls, radius matching, and kernel matching.

### S13: complete, scoped installation packages

[stage_package.py](/Users/Mike/Documents/GitHub/stata-ctools/scripts/stage_package.py)
stages ado/help files, helpers, licenses, and the selected platform binaries into
a fresh directory. Single-platform manifests declare their scope explicitly;
unscoped release manifests still require all four binaries. The checker rejects
missing entries, duplicate entries, and binaries outside the declared scope.
Staging validates before publishing the directory and writes `BUILD_INFO.json`
with the build revision and SHA-256 hashes.

The locally validated package is
[dist/ctools-sep22-validated](/Users/Mike/Documents/GitHub/stata-ctools/dist/ctools-sep22-validated).
It targets Apple Silicon and includes the new prediction helper. A temporary
`net install` verified the helper's installed path, plugin identity, winsorization,
literal decoding and its stored result, bootstrap counts, PPML `xb`, and rejection
of default PPML prediction. User-installed reference packages were not replaced.
The source `build/` remains a development directory with a four-platform manifest;
use the staged package for a local installation.

### S14: complete offline validation and publication policy

The full suite uses 24 official Stata example datasets pinned by URL and SHA-256
in [manifest.json](/Users/Mike/Documents/GitHub/stata-ctools/validation/fixtures/manifest.json).
Preparation downloads and verifies them once; Stata tests only read the local
cache. Missing or changed fixtures fail validation. Generic import round trips
that referenced unavailable dataset names now use available official examples.

[run_stata_audit.py](/Users/Mike/Documents/GitHub/stata-ctools/validation/run_stata_audit.py)
checks reference dependencies, removes the previous log, runs all 20 components,
and requires both per-component completion markers and a unique successful driver
marker. A failed assertion, aborted component, missing summary, empty run, or
unexplained skip cannot pass. The driver captures test failures and still reaches
`exit, clear`, allowing `oldstata` to restore time.

The 48 ATT-SE comparisons are explicitly excluded because `cpsmatch` implements
a different variance estimator from `psmatch2`; ATT comparisons continue to run.
Publication and automatic binary commits require a successful licensed full
Stata job. A missing `CTOOLS_STATA_RUNNER` setting cannot satisfy that requirement.
Pull requests may produce build artifacts without running untrusted code on the
persistent licensed runner, but cannot publish a release through that exception.

### S15 and S16: documentation and maintenance contracts

Installed help now uses the implemented PPML spellings `irlstolerance()`,
`irlsmaxiter()`, and `septolerance()`; documents the prediction/separation limits;
and gives valid bootstrap count and pre-created `weight()` examples. The accepted
but ineffective `cbinscatter, genxq()` option now fails immediately with `r(198)`.
Decode, matching-count, and winsorization help describes the corrected outputs.

The developer guides describe recursive source discovery, dispatcher registration,
cleanup ownership, package helpers, and full validation. They use the required
machine-specific `oldstata`/zsh workflow and retire the obsolete `runstata` and
manual `*_SRCS` build instructions. The public-command inventory check requires
agreement between all 17 commands' help, source directories, package entries, and
full-suite registration. Local Markdown documentation links resolve.

## Additional defects exposed by completing the suite

The offline run reached tests that had previously aborted before execution:

- **Saved OLS fixed effects:** the old implementation stored group means of the
  final residual. It now recovers the additive component `y - Xb - residual`
  and solves for the FE contributions by weighted backfitting. Output writes are
  checked. Saved-FE fits use a tighter internal projection tolerance.
- **Quantile numerical breakdown:** normal equations became unstable when IPM
  weights concentrated on active observations. The solver retries numerical
  breakdown with a scaled weighted Householder QR solve, preserving the original
  convergence criterion. Collinearity detection now accounts for the intercept,
  including constant regressors. Exhausting the requested iteration limit still
  fails rather than reporting convergence.
- **IV time-series syntax:** parsing the first closing parenthesis truncated an
  instrument expression containing `L(1/2).z`. The parser now tracks nested
  parentheses; the regression test executes instead of skipping that case.

Validation corrections were made to repair S14, not to reduce precision. The
seven-significant-figure comparison policy is unchanged. The saved-FE reference
fits now request `reghdfe, tolerance(1e-12)` because its default fit did not agree
with its own tightly converged fit to the precision being tested. Bootstrap
fixtures that requested more draws than available units now assert the native
error contract; valid size fixtures still check exact frequency/row totals.

## Verification and reproducibility

| Check | Result |
|---|---|
| Complete Stata suite | 2,866 passed; 0 failed; 48 documented exclusions; 0 unexpected skips; all 20 components complete |
| New September 22 component | 12 regression groups passed |
| Existing native P1 tests | 4 passed |
| Existing native P2 tests | 6 passed |
| New ASan/UBSan harness | Cluster identity, original-row selection, singleton chains, and allocation-failure cleanup passed |
| Release gate tests | Aborts, stale logs, failed/empty runs, and unexpected skips rejected |
| Packaging/inventory tests | Single/all-platform completeness and all 17 public-command inventories passed |
| Compiler/dependency tests | Compiler availability, CPU baseline, and incompatible dependency rejection passed |
| macOS ARM dependency contract | Passed for rebuilt plugin and compatible static OpenMP |
| Release metadata | Synchronized copies pass `sync_release.py --check` |
| Clean installation | Package installs; expected plugin/helper and command smoke checks pass |
| Working-tree whitespace | `git diff --check` passed |

The complete log is
[validation/audit_ci.log](/Users/Mike/Documents/GitHub/stata-ctools/validation/audit_ci.log).
Detailed temporary build, targeted, installation, and comparison evidence is in
`/tmp/ctools-sep22-fixes/`. The pre-repair snapshot there distinguishes these
repairs from changes already present when this task began. Temporary files may be
removed by the operating system; the substantive results are recorded here.

Runtime validation used StataNow/MP 18.5 on Apple Silicon through the required
`oldstata` wrapper. Its July 1 log timestamps reflect the wrapper's temporary
clock setting. The plugin reports version `1.0.2`, build revision `sep22-fixes`.
The local platform plugin and generic development copy were rebuilt together.
Compilation uses the macOS 11-compatible static OpenMP runtime at
`/tmp/ctools-p2-fixes/libomp-arm21`. Remaining compiler warnings concern unused
functions/parameters/variables; this task did not remove unrelated dormant code.

Native sanitizer tests use Apple Clang on this host. Homebrew Clang 21's ASan
runtime hung before `main()` in an independent startup probe, so it was not used
as evidence of a passing test. Apple Clang's ASan/UBSan run completed normally.

To repeat the full gate after preparing the cache, run from the repository root
outside the sandbox on this machine:

```sh
/bin/zsh -lic 'python3 validation/run_stata_audit.py'
```

The Python runner itself invokes Stata only through `oldstata`. See
[the fixture instructions](/Users/Mike/Documents/GitHub/stata-ctools/validation/fixtures/README.md)
for cache preparation and required reference commands. Linux, Windows, Intel
macOS, minimum Stata versions, and minimum supported operating systems still
require their own execution before making platform-wide validation claims.
