# September 24 audit repairs

This checklist tracks all 21 findings in `SEP24_ASTRA_AUDIT.md`. The original
report and reproduction archive remain unchanged. Repairs were made on top of
the existing working tree; no commit, publication, or pull request was created.

## Itemized changes

| Audit item | Repair | Evidence |
|---|---|---|
| **1 · P1: merge storage loss** | Build the shared-variable schema before writing. Widen numeric types and string widths on the master side, including using-only rows and empty-side paths. Mixed long/float storage promotes to double. | Stata A01/A08; existing merge suite. |
| **2 · P1: destructive output names** | Validate all output names for duplicates and existing variables before creating any output. OLS, range statistics, and encoding restore the dataset when the command fails. Residuals, saved effects, and group IDs require new names. | Stata A02 and A04; OLS/range suites. |
| **3 · P1: lost binsreg adjustment** | Replace the fallback to unadjusted means with a normalized, rank-aware modified Cholesky solve. Omit dependent columns, preserve the retained control span, and return explicit allocation/numerical errors. Use this path with and without absorbed effects. | Stata A03/A16 tests duplicate controls and rescaling by 1e-9; native rank/allocator tests. |
| **4 · P1: ignored solver failure** | Return an explicit status, convergence flag, and iteration count from shared HDFE projection. Propagate failures through OLS, IV, and PPML. Reject nonpositive limits/tolerances and return 430 when projection or PPML IRLS exhausts its limit. | Stata A04; full estimator comparisons. |
| **5 · P1: Excel dates and datetimes** | Correct the 1900 epoch offset; honor the workbook's 1904 epoch. Import the nonexistent 1900 serial 60 as missing. Preserve datetime milliseconds through fractional-day export and datetime import; convert leap-second timestamps with `cofC()`. Add a datetime style distinct from a daily-date style. | Stata A05/A06; explicit 1900/1904 workbook fixtures; daily/datetime round trips. |
| **6 · P1: filename parsing and partial files** | Transfer literal filename, sheet, and missing-value text separately from option tokens, with explicit length validation. Restrict dispatcher thread parsing to its leading token. Write to a reserved sibling temporary file, check write/close results, then publish atomically. Non-replace publication refuses an existing destination. | Stata A05/A06; native temporary-file cancellation, replacement, and no-replace tests. |
| **7 · P1: missing interval keys** | Trim missing keys from each group's searchable population before selecting the optimized or ordinary range-statistics path. Apply the same population rule to unbounded windows and excludeself. | Stata A07 at 999, 1000, 1001, and 2001 rows; full range suite. |
| **8 · P1: incompatible merge keys** | Reject string/numeric key mismatches with 106 regardless of force. Check key types again at the C join entry point before accessing the data union. Force only controls incompatible non-key using values. | Stata A01/A08; native composite-key tests in both directions. |
| **9 · P1: OpenMP capacity assumptions** | Separate usable OpenMP capacity from hardware/pthread capacity. Use logical work-sharing partitions in sample, LSD, counting, and IPS4o sort paths so reduced teams execute every partition. Without OpenMP, choose a valid serial path. | Native 400,001-row numeric/string permutation checks without OpenMP and with thread limits of one and two; full sort suite. Counting sort's native kernel is numeric-only. |
| **10 · P1: invalid cleanup after sort allocation failure** | Zero-initialize owned pointer tables and validate all bucket scratch buffers before publishing the permutation. Initialize short final sampling partitions completely. | Native numeric/string merge and sample tests inject each allocation failure, including aligned allocations, and check resource accounting under sanitizers. |
| **11 · P1: successful zero VCE on failure** | Make shared robust/cluster VCE APIs return status. Propagate failures through OLS, IV, and PPML and guard PPML inference workspaces. Preserve the caller's matrix until every entry is finite. Allocation failures return 920; invalid/nonfinite shared covariance calculations return 498. | Native fault injection covers covariance buffers and cluster-sort allocation, unchanged outputs on failure, invalid clusters, numerical overflow, and valid zero residual variance. Full estimator suites pass. |
| **12 · P2: string merge-sort ordering** | Apply stable character passes from the last byte to the first. Use the appropriate lower/upper bound for left/right merge pivots to preserve equal-key order. | Stata A12 at block and parallel thresholds; native numeric/string order and permutation checks. |
| **13 · P2: truncated export metadata** | Transfer names, storage types, and date styles per column; reject incomplete metadata. Allocate CSV header space from the actual names and escaping requirements. | Stata A13 checks every name and value across 1,100 long-named columns in CSV and XLSX. |
| **14 · P2: omitted-column PPML Wald test** | Select retained coefficient and covariance indices when constructing the joint test, rather than taking the first K entries; use the rank of those restrictions for the test degrees of freedom. | Stata A14 compares explicit joint tests with omitted regressors at the beginning, middle, and end, under robust, clustered, and two-cluster VCE. |
| **15 · P2: partial/empty encoding** | Validate source types, target counts, and new names before the empty-data branch. Treat a multi-variable operation as one dataset transaction, restoring variables and value labels on failure. | Stata A15 exercises failure after an earlier label extension, unchanged source data, and empty-data validation. |
| **16 · P2: binscatter weights and factors** | Materialize expression weights once through a shared helper, exclude missing weights, validate positivity and integral fweights, and materialize factor-variable controls with `fvrevar`. | Stata A03/A16 compares compound weights and factor controls with explicitly generated equivalents. |
| **17 · P2: stale dependency objects** | Compile each platform in a fresh temporary object directory with fail-fast shell behavior. Link an explicit object list, then rename the completed plugin into place. Use a temporary filename with an extension for MinGW compatibility. | Four-platform controlled-compiler tests verify first-failure termination, cleanup, and preservation of a previous plugin. Actual platform builds exercise linking. |
| **18 · P2: lost thread-pool wakeup** | Set the shutdown predicate and broadcast while holding the queue mutex, then join successfully created workers and clear ownership fields. | Native startup failures at worker positions 1–4, repeated 100 times, followed by successful reinitialization under a timeout. |
| **19 · P3: orphaned arena blocks** | Reuse retained successor blocks after reset, for ordinary and aligned allocation. Check rounding, padding, and allocation-accounting overflow. | Native repeated allocate/reset cycles verify bounded capacity, alignment, overflow rejection, and complete freeing. |
| **20 · P3: wrapper infrastructure** | Centralize platform selection and build identity in `_ctools_load`; keep the small caller-local registration required by Stata. Add shared output-name and weight-staging helpers. Retain native `qreg_p` as the active prediction route and remove the unused `cqreg_p.ado` from the package. | Shared-loader repeat checks, all public command suites, and package/helper inventories. |
| **21 · P3: stale documentation** | Update command help, READMEs, developer documentation, and the compatibility matrix. Correct the public auto-sort default, binscatter projection algorithm, implemented date formatting, unsupported timeit examples, IV overidentification results, merge force semantics, failure contracts, and unsubstantiated timing claims. Describe keepcellfmt's limited style-table behavior explicitly. | Metadata, public-command inventory, and package checks; corresponding executable regression examples. |

Round-trip testing also found and repaired discarded Excel inline-string headers
and values. Their arena ownership now survives parallel parsing until cache
creation. Final review added missing permutation frees in range-statistics
allocation-failure paths and the Windows declarations needed by the arena's atomic-store
macro.

## Validation

- **2,877 Stata checks passed, zero failures**, across all 21 suite components.
  The same **48 documented matching-SE method comparisons** remain excluded.
  No numerical comparison tolerance was relaxed.
- The 11 new Stata regression groups are registered in `validate_all.do` and the
  offline release gate. Existing overwrite tests now assert preservation; the
  error-comparison helper restores identical inputs before each command.
- All 23 new native regression groups pass with the system Clang AddressSanitizer and
  UndefinedBehaviorSanitizer. Coverage includes allocation failures, numerical
  overflow, rank normalization, key types, serial sorting, reduced OpenMP teams,
  pool startup failure, arena reuse, and atomic output publication.
- Existing P1, P2, and September 22 native suites pass. Release identity,
  dependency contracts, package completeness, and release-gate tests pass.
- Clang static analysis ran on the changed C translation units. This is supporting
  review evidence, not a claim that every analyzer warning has been eliminated.

The Stata runs used the required `oldstata` wrapper. Its July 1 timestamps reflect
the wrapper's temporary clock adjustment. The package version remains 1.0.2;
rebuilt binaries identify their source revision as `sep24-fixes`.

The companion `SEP24_FIX_EVIDENCE.zip` contains the final test/build logs,
regression sources, this checklist, and SHA-256 hashes of the repaired source
and package files. It is validation evidence, not a cross-platform release.

## Build and compatibility limits

Validation ran on the local macOS ARM Stata runtime. macOS Intel and Windows
cross-builds are checked for their binary dependency contracts; they have not
been executed inside Intel/Windows Stata here. Linux's source changes are ready
for CI, but its binary could not be rebuilt locally because Docker is not running
and no Linux cross-compiler is installed. No Linux plugin is included in the rebuilt artifacts. Do not publish the mixed-platform build
directory as a complete release before the Linux build and licensed platform
gates finish.

`keepcellfmt` retains its limited existing behavior: copying a workbook style
table does not preserve per-cell formatting, other sheets, or arbitrary date
style mappings. Shared string output remains limited to `str2045`; Excel cannot
represent leap seconds. These boundaries are documented rather than presented as
full native-command compatibility.
