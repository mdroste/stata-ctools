# Memory safety audit and repairs — September 28, 2026

All eleven findings below have been addressed in the working tree. Existing unrelated changes were preserved. The original audit evidence and proposed acceptance checks are retained below; the repair summary describes the implemented changes and completed verification.

The strongest evidence is concentrated in failure paths: tracked allocations demonstrate leaks, and native AddressSanitizer/UndefinedBehaviorSanitizer harnesses demonstrate crashes and invalid conversions. Successful baseline runs of the same small fixtures released all tracked allocations. This does **not** establish that every successful command is leak-free.

Priority: **P1** = repair before release because the defect can crash the host, corrupt memory, or silently publish incorrect results; **P2** = repair in the same cleanup work, with less immediate or more conditional impact. “Source-confirmed” means the defective path is present, but I did not reproduce its full runtime consequence.

## Implemented repairs

| Finding | Implemented change |
|---|---|
| 1 | FE remapping initializes empty outputs, publishes counts only after success, and both regression callers check failure before using factors. |
| 2 | Counting remapping checks finite, integral, bounded double codes before conversion; extreme and fractional labels use sorting. |
| 3 | cbinscatter matrix offsets multiply in `size_t` throughout loading, compaction, grouping, binsreg, and residualization; multidimensional malloc calls use checked products. |
| 4 | Effective sample sizes and absorbed/residual degrees of freedom remain doubles through covariance calculations. Invalid covariance denominators return an error. |
| 5 | creghdfe owns initialized buffers through one cleanup path, clearing transferred/released pointers. DOF workspace failures also clean up and return errors. PPML releases its observation map on remap failure; cqreg routes the two leaking allocation failures through cleanup. |
| 6 | Sampling allocation failures release ungrouped observation maps. Grouped loader failures safely release any partial load. |
| 7 | Failed sparse-table builders clear released pointers; cross-group rollback releases every completed table before fallback. |
| 8 | Sample/result stores return errors after cleanup. Regression also propagates scalar/matrix write errors instead of continuing result publication. The current cbsample store-error handling was retained and tested. |
| 9 | XLSX style-parser creation and callback allocation failures return 920. Failure to extract a present styles part returns 610; an absent optional part remains valid. |
| 10 | Timsort cleanup consults initialized numeric pointers instead of uninitialized key-type flags. |
| 11 | Simple range statistics allocate no sorting buffers. Sampling uses only workers with groups to process; bootstrap uses one RNG for its serial draw loop. Scratch concurrency targets 256 MiB, with at least one worker. Explicit thread requests are limited to 1–256 before reaching OpenMP. |

The scratch target applies to sorting/shuffle buffers, not total command memory. A single worker may exceed it for a large group. See [platform documentation](/Users/Mike/Documents/GitHub/stata-ctools/docs/PLATFORMS.md).

## Repair verification

- Passed **707 command allocation/output failure injections** under ASan/UBSan with OpenMP, with zero tracked leaks and successful in-process regression/PPML retries. Added [command failure-injection tests](/Users/Mike/Documents/GitHub/stata-ctools/validation/test_memory_safety_native.py) with tracked allocations and ASan/UBSan, plus fixtures for extreme FE labels, styles, mixed-key timsort, and thread/workspace limits. Regression and PPML retry successfully in the same process after injected failures. Coverage includes weights, clusters, two-way FE/group output, saved effects/residuals, sample flags, scalar/matrix output, and sparse-table fallback. The native suite is registered in Linux CI.
- Added [Stata integration checks](/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_memory_safety.do) and registered them in the complete validation runner. All four integration groups passed through `oldstata`: large/fractional FE labels and **36 billion** total frequency weight; range statistics and thread-limit rejection; filtered/grouped sampling; and cbinscatter controls with/without absorption.
- Seven existing native suites passed: SEP24, new commands, transport (serial and OpenMP), XLSX/strL, PPML, quantile regression, and binary I/O formats. Release-gate/package tests also passed.
- The macOS ARM plugin was rebuilt with Apple Clang and the repository's checksum-pinned LLVM OpenMP 21.1.8 runtime targeting macOS 11. Dependency-contract validation passed. The machine's Homebrew archive was rejected for targeting macOS 26; it was not used for the rebuilt plugin. Runtime sources and temporary build tools remain under `/tmp/ctools-memory-audit`.
- Follow-up static analysis completed for ten affected translation units. Remaining warnings include conservative loop/data-flow paths; warning counts are not treated as confirmed defects. The previous timsort cleanup warning and XLSX parser warnings are absent.

## Original findings and repair plan

The descriptions below record the defects **before repair**. Source links point to the repaired files; quoted line numbers and allocation positions refer to the audit snapshot. Acceptance tests below describe the desired breadth of coverage and are not a claim that every platform or maximum-size case was exercised.


### 1. [P1] Failed FE remapping can crash both regression commands

- **Location:** [remap_counting_impl](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_utils.c), [remap_sort_impl](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_utils.c), [creghdfe caller](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c), [PPML caller](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c).
- **Defect:** The remapper publishes a positive `num_levels` before allocating the corresponding counts. If that allocation fails, it returns an error with `counts_out == NULL`. Both callers ignore the error and subsequently test only whether `num_levels == 0`. The later recount writes through the null counts pointer.
- **Evidence:** With `{1,1,2,2}`, failing the counts allocation produces `rc=-1, num_levels=2, counts=NULL`. Full command harnesses crashed under ASan: two injected failures in weighted `creghdfe`, and one in PPML.
- **Proposed fix:** Initialize every output to an empty state; construct results in local variables and publish them only after all allocations succeed. Record and check every remapper return code, including in the parallel loop, before accessing any factor. Return the memory error through a shared cleanup path.
- **Acceptance test:** Fail each remapper allocation in the counting and sorting implementations, with and without weights; both commands must return a memory error, leak nothing, and succeed on the next invocation.

### 2. [P1] FE codes can overflow conversions and range arithmetic

- **Location:** [remap_and_count](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_utils.c).
- **Defect:** The routine converts a double to `int64_t` before checking whether it is representable. It also computes `vmax - vmin + 1` in signed 64-bit arithmetic without checking the span. Both operations can have undefined behavior even when the original doubles are valid nonmissing Stata values.
- **Evidence:** UBSan reports an invalid conversion for FE values `1e30`. Values including `-2^63` and `2^62` trigger signed overflow in the range subtraction. These are current-source reproductions, independent of the older counting-sort bug.
- **Proposed fix:** Check finiteness, integrality, representability, and the usable range in double precision before casting. Use the sort-based remapper whenever the counting representation is unsuitable. FE values are labels; large or fractional codes need not be rejected.
- **Acceptance test:** Negative and fractional codes, values on either side of the signed 64-bit boundaries, and widely separated codes must produce the same partition as sorting the original doubles.

### 3. [P1] cbinscatter still calculates large matrix offsets in 32 bits

- **Location:** [initial controls copy](/Users/Mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_impl.c), [missing-value scan](/Users/Mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_impl.c), [group copies](/Users/Mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_impl.c), [binsreg adjustment](/Users/Mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_binsreg.c), [residualization](/Users/Mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_resid.c).
- **Defect:** Allocation sizes are often checked in `size_t`, but expressions such as `j * N`, `k * N + i`, and `k * n_group + j` still multiply signed `ST_int` values. A correctly sized allocation therefore does not prevent a wrapped index and out-of-bounds access.
- **Evidence:** Source-confirmed. For example, 100 controls and 22 million observations give `99 * 22,000,000 = 2,178,000,000`, exceeding `INT_MAX`. I did not allocate the multi-gigabyte fixture on this machine.
- **Proposed fix:** Convert an operand to `size_t` before every matrix-offset multiplication; keep row counts distinct from byte counts and offsets. Audit the whole cbinscatter pipeline, including compacted groups and residualization.
- **Acceptance test:** Check index arithmetic around `INT_MAX`, then run a large-data sanitizer case on an appropriately provisioned host. Small matrices should remain numerically identical.

### 4. [P1] Frequency-weight totals overflow the effective sample size

- **Location:** [regression degrees of freedom](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c), [VCE effective N](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c).
- **Defect:** `sum_weights` is converted to `ST_int`. A frequency-weighted sample can represent more than 2^31−1 observations while containing only a few physical rows. The invalid conversion then contaminates degrees of freedom and covariance calculations.
- **Evidence:** A 12-row native fixture with frequency weights of 300,000,000 produces `sum_weights = 3.6e9`; UBSan stops at the conversion on line 1078.
- **Proposed fix:** Carry effective sample sizes and degrees of freedom in a representation that supports large weight sums throughout the VCE interface. Validate before any narrower conversion; an explicit unsupported-size error is preferable to undefined arithmetic if a downstream interface cannot be widened.
- **Acceptance test:** Compare estimates and covariance scaling below and above 2^31−1 total weight, without creating billions of physical rows.

### 5. [P2] Regression commands have incomplete cleanup on early returns

- **Location:** [creghdfe remap failure](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c), [compaction failures](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c), [PPML remap failure](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c), [cqreg solver allocation failure](/Users/Mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_regress.c), [cqreg sparsity allocation failure](/Users/Mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_regress.c).
- **Defect:** Once the loaded columns are freed, `obs_map` remains locally owned. Numerous creghdfe exits omit it and `cluster_raw_values`; some also omit weights, weighted counts, or replacement means. PPML loses `obs_map` on remap failure. cqreg allocates `is_collinear` locally but omits it on the two later allocation-failure returns above.
- **Evidence:** Of 70 injected allocation positions in a small weighted, clustered creghdfe fixture, **41 left tracked allocations live**, even after global-state cleanup; two others crashed as described in finding 1. A representative failure leaked 48 bytes of observation mapping plus 96 bytes of cluster values. Other failures retained four buffers. PPML remap failure retained its 48-byte observation map. cqreg's two leaks are source-confirmed.
- **Proposed fix:** Give each command one initialized ownership structure and one cleanup routine. Route exits through it; clear pointers when ownership moves into global state. Preserve the required aligned deallocator for production observation maps. Normalize allocation failures to the intended memory return code.
- **Acceptance test:** An allocation-failure sweep of the complete command, not only its numerical kernels, must leave zero tracked local allocations after cleanup. Cover weights, clusters, sample flags, saved effects, singleton removal, and solver failures.

### 6. [P2] Filtered sampling leaks its observation map on allocation failure

- **Location:** [csample RNG allocation](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c), [csample workspace allocation](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c), [cbsample group allocation](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c), [cbsample RNG allocation](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c).
- **Defect:** Without by/cluster/strata variables, the filtered observation map is a standalone allocation. These cleanup branches free a filtered-data object only when grouping is present and omit the standalone map otherwise.
- **Evidence:** Complete native command harnesses with 12 selected observations: csample failures 4 and 5 each retain 48 bytes; cbsample failures 3, 4, and 5 do likewise. Successful runs return to zero allocations. The leak is four bytes per selected observation in the production index type.
- **Proposed fix:** Centralize cleanup and explicitly distinguish an owned standalone map from an alias into `by_filtered`.
- **Acceptance test:** Fail every allocation with filtered and unfiltered input, with and without grouping; verify both no leaks and no double frees.

### 7. [P2] crangestat loses completed sparse tables when a later table fails

- **Location:** [cross-group sparse-table construction](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c), [conditional cleanup](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c).
- **Defect:** When the second or later source's table fails to build, the code sets `group_sparse_ptr = NULL` and never sets `sparse_built`. Tables already built for earlier sources are then abandoned. The slower fallback completes and the command reports success.
- **Evidence:** With eight groups of 64 rows and two source variables, an injected failure in the second table returns **rc=0 while leaking four allocations totaling 7,540 bytes**. The normal case leaks nothing.
- **Proposed fix:** Retain the table array, track which entries are initialized, and free all completed entries before fallback. Make failed builders clear pointers they have freed so common cleanup is safe.
- **Acceptance test:** Fail each allocation for every source-table position, including repeated failures across groups and invocations. Check output equivalence with the fallback and zero retained allocations.

### 8. [P1] Output failures are ignored or overwritten

- **Location:** [creghdfe sample flag](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c), [VCE assignment to the same status](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c), [crangestat stores](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c), [csample stores](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c), [cbsample stores](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c).
- **Defect:** creghdfe records a failed `e(sample)` flag write, then replaces that error with the successful robust/cluster covariance status. The sampling commands and crangestat do not check their numeric-store return codes. A command can report success with incomplete or absent output.
- **Evidence:** Making every sample-flag store return 459 still gives rc=0 from the clustered creghdfe fixture; the unadjusted fixture preserves 459. Making every crangestat result store return 459 also gives rc=0. Sampling store omissions are source-confirmed.
- **Proposed fix:** Preserve the first error; stop dependent processing and enter cleanup. Parallel store loops need a shared error accumulator and must propagate it after the workers finish. Do not post successful estimates or timing/status metadata as if results were committed.
- **Acceptance test:** Fail the first, middle, and last numeric/string/matrix store; confirm the command returns an error and the ado wrapper handles any partially staged outputs appropriately.

### 9. [P1] XLSX allocation failures silently discard date styles

- **Location:** [styles callback](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c), [styles parser wrapper](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c), [XML callback-stop semantics](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx_xml.c).
- **Defect:** Failure to create the styles parser returns success. Failure to expand style arrays stops the callback, but the XML parser treats callback-requested stopping as success, and the styles wrapper has no failure flag. A workbook with styles can therefore be imported with missing style information and incorrect date/datetime interpretation.
- **Evidence:** A native fixture containing General and date styles normally gives rc=0 and a date style. Failing any of the four allocations instead gives **rc=0 with zero styles**.
- **Proposed fix:** Distinguish an absent optional styles part from failure to read a present part. Return an allocation error on parser creation failure; carry an explicit callback failure status, as the shared-strings parser already does. Abort the import before publishing values whose interpretation depends on the lost styles.
- **Acceptance test:** Fault-inject parser creation and both array-growth allocations; imports must either retain all styles and correct dates or return an error.

### 10. [P2] Multi-key timsort reads uninitialized cleanup flags

- **Location:** [context allocation and initialization](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_timsort.c), [cleanup loop](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_timsort.c).
- **Defect:** `is_numeric` is allocated without initialization. If an early numeric-key allocation fails, cleanup reads flags for all later, unvisited keys. Their corresponding numeric pointers are initialized to NULL, but the left operand of `&&` is still an uninitialized read.
- **Evidence:** Clang's path-sensitive analyzer identifies this allocation-failure path; direct source inspection confirms it. I did not demonstrate a resulting invalid free or crash.
- **Proposed fix:** Initialize all flags, or free only the initialized numeric-pointer array entries without consulting the flags.
- **Acceptance test:** Fail conversion of the first and middle keys in mixed numeric/string sorts. Use MemorySanitizer on a supported host in addition to allocation tracking; ASan alone does not detect uninitialized integer reads.

### 11. [P2] Excess scratch allocation magnifies memory pressure

- **Location:** [crangestat workspace](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c), [csample workspace](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c), [thread option validation](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_plugin.c), [thread setter](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_threads.c).
- **Defect:** crangestat allocates a maximum-group-sized double workspace for every configured thread even for count/mean operations that do not need it. csample similarly allocates one maximum-group-sized index workspace per configured thread even when there is only one group to process. The dispatcher accepts any positive representable thread count, and the setter forwards it to OpenMP without a practical bound.
- **Impact:** These allocations are normally freed; this is excessive peak memory, not a leak. A 10-million-row group and 64 configured threads require 5.12 GB of such scratch space alone. Allocation failures then expose the cleanup defects above.
- **Proposed fix:** Allocate scratch only for statistics/algorithms that use it, limit workers to available independent tasks, and use checked products and an explicit workspace budget. Validate unreasonable thread requests before changing the OpenMP runtime.
- **Acceptance test:** Track peak live bytes for one large group, many small groups, and simple versus percentile statistics. Verify that simple crangestat operations allocate no sorting workspace and a single sampling group does not allocate one buffer per idle thread.

## Implementation order

1. **Repair shared remapping first (findings 1–2).** Both regression commands inherit its failure contract; fixing their local cleanup without this contract still leaves crashes.
2. **Unify command ownership and cleanup (findings 5–7, 10).** Add allocation-failure sweeps covering complete commands, including retained/global state. Require a successful retry after every injected failure.
3. **Eliminate remaining arithmetic overflow (findings 3–4).** Widen offsets before multiplication and effective counts before their downstream use. Keep unsupported-size errors explicit.
4. **Preserve errors through result publication (findings 8–9).** Exercise output callbacks and parser callbacks separately from allocation cleanup.
5. **Reduce peak workspace and make checks a release requirement (finding 11).** Add sanitizer, leak-accounting, and failure-injection cases for the repaired paths. Repeat on Linux and Windows, including Windows aligned-allocation behavior.

## Original audit coverage

- Clang static analysis completed for all **83 first-party C translation units**, including their referenced headers and `.inc` files, using the macOS configuration without OpenMP. It emitted 88 warnings; these are candidates, **not 88 confirmed defects**. Many zero-size and inconsistent-loop-path warnings did not survive inspection and are not reported as findings.
- Seven existing native test scripts passed with ASan/UBSan enabled where supported by their runners: `test_sep24_native.py`, `test_newcommands_native.py`, `test_transport_native.py`, `test_xlsx_strl_native.py`, `test_cpplmhdfe_native.py`, `test_cqreg_native.py`, and `test_io_formats_native.py`. Transport passed both serial and OpenMP variants. Initial runs of two suites timed out with the PATH-selected toolchain; rerunning unchanged tests with `/usr/bin/clang` passed. This is not evidence of a product deadlock.
- New temporary harnesses compiled the actual command/kernel sources against a mocked Stata interface, with tracked allocations and selected allocation/store failures. Reproduction sources, scripts, and logs are in [/tmp/ctools-memory-audit](/tmp/ctools-memory-audit). Useful entries are `regress.c`, `faults.py`, `regress-faults.json`, `ppml.c`, `sampling.c`, `rangestat.c`, `rangestat-faults.txt`, `remap.c`, `regress_weight.c`, `regress_store.c`, and `xlsx_styles.c`.
- Earlier findings were rechecked where relevant. The current all-singletons cleanup releases its observation map, the current counting-sort conversion guard handles extreme keys, and the current IPS4o string routine bounds its recursive depth. They are not presented here as open bugs.

## Limits

This work addresses the eleven audited issues; it is not an exhaustive statistical or format-parity review. The native harnesses mock host callbacks and loading. Targeted Stata integration passed, but live Stata resident memory was not profiled. Linux/Windows runtime checks, MemorySanitizer, and a cbinscatter dataset whose physical matrix exceeds 2^31 cells remain untested locally. The widened offsets were reviewed and the ordinary-size cbinscatter integration checks passed. Third-party codecs were exercised by existing tests but not independently reviewed line by line.
