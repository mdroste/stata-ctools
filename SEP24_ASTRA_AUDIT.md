# ctools audit — September 24, 2026

The audit identifies **21 items: 11 P1, 7 P2, and 3 P3**. The most consequential problems are silent data loss, incorrect statistical results, and unsafe native-code failure paths. Each item below includes the evidence, the cause, a proposed fix, and acceptance tests. No implementation changes were made.

**Scope and verification.** This is an audit of the current working tree, including its existing uncommitted repairs, rather than of the last Git commit. I reviewed the 17 public command wrappers and their C implementations, shared infrastructure, validation and packaging scripts, build configuration, and documentation. I built a fresh Apple Silicon plugin in `/tmp/ctools-sep24-audit`, ran the complete Stata validation suite, ran the existing native/build/release checks, analyzed all 80 first-party C translation units with Clang, and added targeted reproductions in that isolated directory. Hashes of 420 input files still matched the working tree at the end of testing.

The existing Stata suite reports **2,866 passed, zero failed, and 48 expected exclusions** for matching standard-error methods. All 20 validation components completed. The existing P1, P2, September 22, release-gate, build-contract, release-metadata, and package-metadata checks also passed. The findings below concern cases that those checks do not cover. Warnings from static analysis were investigated rather than counted as bugs by themselves.

Execution was on macOS Apple Silicon. Linux, Windows, and Intel plugins were not executed. Vendored compression/XML dependencies were reviewed at their integration boundaries, not subjected to an independent exhaustive security audit. A passing baseline does not establish correctness outside the cases it exercises.

**Priorities.** P1 means fix before the next release because the defect can corrupt data or results, or crash the host process. P2 means fix in the next correctness/reliability cycle; the trigger is narrower or the impact more contained. P3 means maintenance or a currently unused API defect. Items are ordered approximately by practical urgency within these groups.

The companion [evidence archive](/Users/Mike/Documents/GitHub/stata-ctools/SEP24_AUDIT_EVIDENCE.zip) contains the reproduction do-files, native harnesses, logs, build log, and source hashes. These are diagnostic reproductions, including expected failures, not additions to the release test suite. Stata log timestamps show July 1 because the required `oldstata` wrapper temporarily changes the machine clock. The audit date is September 24.

**1. P1 — `cmerge` silently loses values when shared variables need wider storage types.**

**Evidence:** Stata reproduction, `repro3.do`. The master has `byte value` containing missing and `str1 text` containing an empty string. The using file has `double value=1000.125` and `str20 text="abcdefghijklmnopqrst"`. After `cmerge 1:1 id using ..., update`, both the matched row and a using-only row contain missing numeric values and the string `"a"`. The command returns success.

**Cause:** [cmerge.ado:838](/Users/Mike/Documents/GitHub/stata-ctools/build/cmerge.ado:838) builds storage definitions for new variables but excludes shared variables. Their master storage types remain in place when the plugin writes using values. A correct C-side double or string cannot survive an incompatible Stata destination type.

**Proposed fix:** Build a complete output schema before any dataset mutation. For each shared numeric variable, select a lossless common storage type based on both inputs; for strings, use a width sufficient for both, with an explicit supported policy for `strL`. Apply the same reconciliation to keys that receive using-only values. Recast destination variables inside the existing rollback boundary before writing. Reject incompatible numeric/string data types according to an explicit `force` policy; do not reinterpret their representations. Preserve formats and labels deliberately when choosing the output schema.

**Acceptance tests:** Compare against native `merge` for byte/long/float/double combinations, fractional and out-of-range values, narrow/wide strings, and using-only rows. Cover plain merge, `update`, `replace`, and each supported merge cardinality. Verify that errors leave the original data intact.

**2. P1 — Output names can overwrite inputs before computation.**

**Evidence:** Stata reproduction, `repro1.do`. `creghdfe y x, residuals(x)` returns `r(504)` after replacing all 600 values of `x` with missing. `crangestat (mean) v=v, interval(key -1 1)` returns success after replacing all 999 source values with missing.

**Cause:** [creghdfe.ado:313](/Users/Mike/Documents/GitHub/stata-ctools/build/creghdfe.ado:313) drops the requested residual destination before loading the regression inputs. Its `groupvar()` and saved-FE destinations need the same review. [crangestat.ado:204](/Users/Mike/Documents/GitHub/stata-ctools/build/crangestat.ado:204) clears existing numeric output variables while parsing statistics, allowing an output to alias a source, interval key, or grouping variable.

**Proposed fix:** Parse and validate the complete input/output name graph first. Reject input/output aliases unless the command explicitly supports replacement and can stage it safely. Require new names where appropriate. Compute into temporary variables and rename or copy only after all calculations and stores succeed. Include dataset mutation, output labels, and estimation-result posting in the success/rollback design. A failure must not destroy existing variables.

**Acceptance tests:** Alias each output with the dependent variable, regressor, weight, FE, range key, source, and `by()` variable. Include duplicate output names and failures after the first successful store. Compare values, storage types, labels, observation order, and variable order before and after rejected calls.

**3. P1 — Collinear `binsreg` controls silently remove all adjustment.**

**Evidence:** Stata reproduction, `repro5.do`. With `z2=2*z`, changing `controls(z)` to `controls(z z2)` under `method(binsreg)` changes bin data by as much as **8.7091208**. The result with both controls is exactly equal to the result with no controls. The command reports no error.

**Cause:** [cbinscatter_binsreg.c:181](/Users/Mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_binsreg.c:181) sets every control coefficient to zero when Cholesky fails. This changes the requested estimand from adjusted to raw bin means.

**Proposed fix:** Use a rank-revealing solve that drops only redundant control columns and preserves the fitted control component. Expose or consolidate the collinearity-aware solve already present in the residualization code, with consistent tolerances and an explicit retained-column map. Distinguish rank deficiency from allocation and numerical failures. If an admissible adjusted solution cannot be obtained, return an error rather than raw means.

**Acceptance tests:** Adding a duplicate, rescaled duplicate, zero control, or linear combination must leave adjusted bin estimates invariant within tolerance. Repeat with weights, group splits, and fixed effects. Test near-collinearity separately from exact redundancy.

**4. P1 — HDFE solver failure is ignored, and OLS reports convergence incorrectly.**

**Evidence:** Stata reproduction, `repro3.do`, seed 742913. On the same two-way-FE sample, `iterate(1) tolerance(1e-14)` returns success, announces convergence after one iteration, and reports a coefficient of **0.609696879573123**. Allowing convergence gives **0.678911076385502** in 12 iterations.

**Cause:** [creghdfe_solver.c:358](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_solver.c:358) returns a negative iteration count on exhaustion. [creghdfe_regress.c:1028](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c:1028) takes its absolute value and proceeds; [creghdfe.ado:778](/Users/Mike/Documents/GitHub/stata-ctools/build/creghdfe.ado:778) prints convergence. [civreghdfe_impl.c:750](/Users/Mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:750) and its later projection call ignore the return entirely. [cpplmhdfe_irls.c:869](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:869) similarly removes the failure sign; other projection calls ignore it. These IV/PPML paths are source-confirmed, rather than separate end-to-end reproductions of coefficient error.

**Proposed fix:** Replace the signed-count convention with a result containing status, convergence, and iterations. Propagate it through every projection caller. Validate iteration limits and tolerances. By default, return a convergence error such as `r(430)` without posting a successful fit or committing outputs. If unconverged estimates are intentionally supported, make that an explicit documented policy and expose the status consistently. PPML already distinguishes outer IRLS nonconvergence in its messages and `e(converged)`; retain that distinction and add checks for the inner projections.

**Acceptance tests:** Force nonconvergence with restrictive limits in OLS, IV, and PPML, including the PPML final VCE projection. Verify status, messages, return code, and output preservation. Successful results must remain invariant to increasing an already sufficient iteration limit.

**5. P1 — Excel date import shifts dates, and datetime export discards the time of day.**

**Evidence:** `repro2.do` exports Stata daily dates 0, 1, and 2 using native `export excel`; `cimport excel` reads them as **-1, 0, and 1**. In `repro5.do`, exporting `%tc` value `24sep2026 12:34:56` writes the integer Excel serial **46289**, with no fractional day. Native Stata import consequently recovers a daily date without the time.

**Cause:** [cimport_xlsx.c:2389](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:2389) combines the 21916-day offset with an additional post-February-1900 decrement. The epoch already accounts for the modern-date offset. The reader also has no handling of workbook `date1904` metadata. [cexport.ado:606](/Users/Mike/Documents/GitHub/stata-ctools/build/cexport.ado:606) applies `dofc()`/`dofC()` to datetimes, losing their sub-day component before writing.

**Proposed fix:** Centralize Excel/Stata date conversion. Handle the 1900 system piecewise around Excel's fictitious leap day, explicitly define serial 60 behavior, and read the workbook's 1904-system flag. For datetimes, export the fractional day and a datetime number format; specify how `%tC` leap-second-aware values map to Excel's time convention. Keep numeric date cells separate from requests to export formatted strings. Apply the inverse conversion consistently on import.

**Acceptance tests:** Native-Stata and XML-level fixtures for serials 1, 59, 60, 61, January 1960, modern dates, 1904 workbooks, fractional days, midnight boundaries, and millisecond precision. Test both import and export independently so matching mistakes cannot cancel in a round trip.

**6. P1 — Filenames containing spaces are parsed incorrectly and written to the wrong path.**

**Evidence:** `repro3.do` calls both CSV and XLSX export with a path ending in `space name.*`. Both return success; neither requested file exists. The writers use the path ending in `space` instead.

**Cause:** [cexport_parse.c:96](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_parse.c:96) and [cexport_xlsx.c:470](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_xlsx.c:470) tokenize a combined filename/options string on whitespace. Stata's filename validation therefore checks a different path from the one C opens. Separately, [ctools_plugin.c:98](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_plugin.c:98) searches the whole string for `threads(` and removes it, so option-shaped data is not safely separated from options. The space-path failure is reproduced; the substring issue is source-confirmed.

**Proposed fix:** Pass filenames and sheet names as distinct plugin arguments or dedicated length-checked metadata, not embedded tokens. Parse options only within their own fields. Ensure existence checks, opening, writing, and reported paths refer to precisely the same filename. Write to a temporary sibling file and commit after successful completion, preserving no-replace semantics at the final operation. Review all commands sharing the same dispatcher contract.

**Acceptance tests:** Spaces, tabs where allowed, Unicode, parentheses, and names containing `threads(2)` or other option text. Include an existing file at the truncated path and verify that it is unchanged. Assert the actual requested artifact exists and is readable after every reported success.

**7. P1 — Optimized range statistics include rows with missing interval keys.**

**Evidence:** `repro1.do` creates 1,000 observations, with one missing key and an extreme source value on that row. Native `rangestat` returns mean **1** for `interval(key . .)`; `crangestat` returns **2**. Removing one ordinary observation changes the ctools result back to **1**, showing dependence on the optimization threshold.

**Cause:** The optimized paths at [crangestat_impl.c:1888](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:1888) and [crangestat_impl.c:2136](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:2136) search to the physical group end. For an unbounded upper window, Stata's finite numeric missing codes fall below `DBL_MAX` and enter the source window. The ordinary path computes an end excluding missing keys.

**Proposed fix:** Compute a nonmissing-key boundary once per group and use it in every window search, prefix-sum, sparse-table, and run-based implementation. Keep missing-key target rows out of computation as well. Factor out the boundary logic so optimized branches cannot acquire different sample definitions.

**Acceptance tests:** Samples immediately below, at, and above every dispatch threshold; one and many groups; ordinary and extended missing keys; one-sided and fully unbounded intervals; `excludeself`; all supported statistics. Results must be invariant to the selected execution path.

**8. P1 — `cmerge, force` permits incompatible key types that C interprets through the wrong union member.**

**Evidence:** A native comparator harness intercepts `strcmp` without dereferencing its arguments. Comparing a string master key with numeric using key `1.0` supplies **0x3ff0000000000000** as the second string pointer. This is the bit representation of the double, not a string address. No live Stata crash was needed to establish the invalid access.

**Cause:** [cmerge.ado:460](/Users/Mike/Documents/GitHub/stata-ctools/build/cmerge.ado:460) allows numeric/string key mismatches under `force`. [cmerge_keys.c:105](/Users/Mike/Documents/GitHub/stata-ctools/src/cmerge/cmerge_keys.c:105) chooses both union accesses using only the master variable's type. [cmerge_join.c:50](/Users/Mike/Documents/GitHub/stata-ctools/src/cmerge/cmerge_join.c:50) routes such a pair to that general comparator.

**Proposed fix:** Require matching numeric-versus-string key types regardless of `force`, unless an explicit, separately implemented conversion is requested. Validate again at the C entry point so alternate callers cannot bypass the contract. Reserve `force` for a clearly defined policy on non-key variables. Do not try to handle mismatched representations inside an ordinary comparison routine.

**Acceptance tests:** Both directions of key-type mismatch, with and without `force`, including composite and unmatched keys. Assert a deterministic type error before data mutation. Add native type-contract tests with sanitizers.

**9. P1 — The supported no-OpenMP build can crash in sample sort.**

**Evidence:** `noomp_sample.c`, compiled without OpenMP and with four reported CPUs, invokes the public order-only sample-sort path on 400,001 strings. AddressSanitizer reports a read of an uninitialized pointer in `strcmp`, reached from the string sample-sort implementation.

**Cause:** [ctools_threads.c:82](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_threads.c:82) can report several available CPUs without OpenMP. [ctools_sort_sample.c:855](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_sample.c:855) uses that count to select the parallel implementation. When OpenMP pragmas are ignored, only the code for thread zero executes, while sampling and bucket metadata still assume several participants. The Makefile explicitly supports a pthread-only fallback.

**Proposed fix:** Separate hardware/pthread capacity from the number of usable OpenMP workers. Force these algorithms to a correct serial implementation when OpenMP is absent. In OpenMP builds, either obtain the actual team size before partitioning or distribute all logical partitions with work-sharing constructs that remain correct when the runtime supplies fewer workers. Review analogous manually partitioned regions in the other sort engines.

**Acceptance tests:** Build with OpenMP disabled and run numeric, string, and mixed-key sorts above all parallel thresholds under ASan/UBSan. Also exercise OpenMP dynamic teams and thread limits. Check permutation completeness, key order, stability where promised, and payload integrity.

**10. P1 — Sort allocation failures lead to invalid frees.**

**Evidence:** Independent merge-sort and sample-sort fault harnesses fail the first ordinary allocation and fill later allocations with a recognizable pattern. Each records **four invalid free attempts** during cleanup, instead of a clean memory-error return.

**Cause:** [ctools_sort_merge.c:462](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_merge.c:462) allocates `thread_temps` with `malloc`, then can jump to cleanup before its elements are initialized. Cleanup frees every element. [ctools_sort_sample.c:407](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_sample.c:407) has the same ownership ordering. Review both numeric and string variants.

**Proposed fix:** Zero-initialize pointer tables at allocation, or initialize them immediately before any possible jump to cleanup. Track the number of successfully constructed resources where zero initialization alone is insufficient. Give each routine a single explicit ownership/cleanup structure, and avoid mutating a caller's permutation until prerequisites are available.

**Acceptance tests:** Inject failure at each allocation position, including aligned allocation and reallocation, across all sort implementations. Require a memory-error status, no invalid frees or leaks, and a documented input-preservation guarantee. Run these tests under sanitizers rather than relying only on allocator interception.

**11. P1 — Covariance allocation failures become successful zero standard errors.**

**Evidence:** `fault_vce.c` injects an allocation failure in the shared robust covariance routine. It returns with **V=0**, and its `void` API provides no error signal to the estimator.

**Cause:** [ctools_matrix.c:176](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_matrix.c:176) and [ctools_matrix.c:273](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_matrix.c:273) zero the output matrix and return silently on allocation or cluster-sort failure. Their callers cannot distinguish a computed covariance from a failed calculation. This shared code is used by estimation routines, including OLS and PPML paths.

**Proposed fix:** Return a status code from both covariance functions and propagate it through every wrapper and estimator. Compute into owned temporary storage and publish `e(V)` only on success. Distinguish allocation failure, inadmissible cluster/degree-of-freedom configurations, and numerical failure. A zero matrix must represent a mathematical result, not an error sentinel.

**Acceptance tests:** Inject each covariance-buffer allocation failure and cluster-sort failure. Verify a nonzero command status and absence of newly posted successful estimation results. Keep separate numerical tests for valid zero residual variance so legitimate zeros remain supported.

**12. P2 — The merge-sort string kernel does not implement lexicographic order.**

**Evidence:** `repro3.do` tests all seven explicit algorithms on alternating `"az"`/`"ba"` keys. `csort key, algorithm(merge) nosortedby threads(2)` produces an inversion at N=32, 100, and 30,000. N=31 takes a working small-block path. The other tested algorithms pass these cases.

**Cause:** [ctools_sort_merge.c:193](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_merge.c:193) applies stable character passes from the first character to the last across the whole block. Later characters therefore dominate earlier characters. The comment calls this MSD sorting, but there is no recursive prefix partition. The wrapper's usual final native sort can mask the defect; `nosortedby` exposes it.

**Proposed fix:** Use stable LSD passes from the last character to the first, or implement genuine MSD recursion within equal-prefix buckets. Prefer sharing an already tested string kernel over maintaining a third version. Preserve the intended ordering for empty strings and unequal lengths, and verify stable merging of ties.

**Acceptance tests:** Test the plugin permutation directly, without a final Stata sort, across algorithm and block thresholds. Include empty strings, prefixes, UTF-8 byte sequences, identical keys, multiple keys, and payload columns. Check both sortedness and stability where the API promises it.

**13. P2 — Wide exports silently rename columns when metadata buffers truncate.**

**Evidence:** `repro2.do` exports 1,100 variables named `long_variable_name_for_test_1` through `_1100`. The first damaged CSV header is column **1027**, named `long_variable_na`; subsequent columns fall back to names such as `v1098`, `v1099`, and `v1100`. The export succeeds.

**Cause:** [cexport_parse.c:169](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_parse.c:169) reads all names into a 32,768-byte buffer and accepts incomplete metadata, synthesizing the remaining names. On this Stata runtime, the macro read can truncate while returning success. [cexport_xlsx.c:542](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_xlsx.c:542) has similar fixed buffers for names and types; XLSX corruption was not separately reproduced.

**Proposed fix:** Transfer metadata per variable or in length-checked chunks, or allocate based on a reliable explicit length protocol. Require exactly the expected number of complete names and types. Reject missing metadata instead of inventing replacement names. Size the output header dynamically as well, and perform all metadata validation before opening the destination.

**Acceptance tests:** Cross each metadata-byte boundary using long and short names, with CSV and XLSX and enough columns to stress type metadata. Assert every output name and storage interpretation, not just row/column counts. Include a final token cut in the middle of a name.

**14. P2 — PPML's reported joint Wald statistic depends on the position of an omitted regressor.**

**Evidence:** `repro2.do`, seed 92815. With `zero=0`, `cpplmhdfe y zero x z, absorb(g)` reports **e(F)=70.62785**, while `test x z` gives **80.18204**. Reordering the same regressors to `x z zero` makes the displayed statistic 80.18204. The estimable coefficients are unchanged.

**Cause:** [cpplmhdfe.ado:461](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe.ado:461) takes the first `K_keep` entries of matrices that retain omitted regressors in their original positions. That is not the same as selecting retained regressors.

**Proposed fix:** Construct the joint-test coefficient and covariance matrices using the retained-column index map, excluding omitted terms and any intercept as required by the model-test definition. Use the rank of the tested restrictions for the numerator degrees of freedom. Alternatively, delegate the restriction calculation to Stata's tested postestimation machinery after posting consistent coefficient stripes.

**Acceptance tests:** Omitted variables at the beginning, middle, and end; factor-variable base levels; duplicate regressors; robust and clustered covariance. Compare the reported statistic with an explicit joint test and require permutation invariance.

**15. P2 — `cencode` partially commits failed multi-variable operations and bypasses validation on empty data.**

**Evidence:** `repro5.do` runs `cencode a b, replace` where `a` is a short string and `b` is a 4,096-byte `strL`. The command fails with `r(920)` on `b`, but `a` has already become numeric. In `repro1.do`, an empty-data call with an existing generated target succeeds and changes its storage type, while native `encode` returns `r(110)`.

**Cause:** [cencode.ado:175](/Users/Mike/Documents/GitHub/stata-ctools/build/cencode.ado:175) processes and commits variables one at a time; label changes are also immediate. The early `_N==0` branch at [cencode.ado:32](/Users/Mike/Documents/GitHub/stata-ctools/build/cencode.ado:32) runs before ordinary source-type, output-count, and existing-target checks. Its comment explicitly accommodates leftover test variables, allowing test scaffolding to alter production semantics.

**Proposed fix:** Validate the entire request before treating the empty-data case. Stage every encoded variable and value-label definition, then commit the whole request only after all variables succeed. Preserve or restore existing label definitions when extending them fails. Move cleanup of leftover test variables into the tests themselves. Keep unsupported long strings a clear early error unless the advertised support is expanded.

**Acceptance tests:** Failure on the second or last input must leave all earlier source variables and labels unchanged. Exercise empty data with numeric sources, existing targets, duplicate targets, incorrect output counts, `replace`, `label()`, and `noextend`.

**16. P2 — `cbinscatter` accepts syntax that its plugin call cannot consume.**

**Evidence:** `repro5.do`: `[aw=1+abs(z)]` fails with `r(198)` and “1 invalid name”; `controls(i.category)` on a valid nonnegative integer category fails with `r(101)` and “factor-variable and time-series operators not allowed.” The wrapper declares weights and factor-variable controls.

**Cause:** [cbinscatter.ado:132](/Users/Mike/Documents/GitHub/stata-ctools/build/cbinscatter.ado:132) strips `=` from the weight expression and treats the result as a variable name. Controls are passed to the plugin without materializing the factor-variable expansion.

**Proposed fix:** Evaluate the weight expression into a temporary numeric variable on the candidate sample, validate the selected weight type, and include weight missingness in sample construction. Materialize factor-variable controls with Stata's expansion tools while retaining a map to user-facing terms. Follow the repaired estimator wrappers' staging approach rather than creating another parser. If a subset of this syntax is intentionally unsupported, reject it explicitly at parsing and narrow the documentation; do not advertise `fv` and fail only at the plugin call.

**Acceptance tests:** Constant and compound weights, missing/invalid weights, factor controls, interactions, base levels, and `if`/`in`. Compare with the equivalent call using manually generated weight and indicator variables.

**17. P2 — Dependency compilation loops do not stop on the first compiler failure.**

**Evidence:** Source review of the platform build recipes, beginning at [Makefile:417](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:417). This is a build-reliability finding, not a reproduced bad release artifact.

**Cause:** The libdeflate shell loops issue one compile per source without `set -e` or checking each command's result. The loop's exit status follows its final command. A failed earlier compile can therefore be hidden by a later success. Linking uses wildcard object lists; an object left from an earlier failed build can satisfy the link and silently carry old code into a new plugin. A clean directory may instead fail at link time, so successful stale output depends on leftovers.

**Proposed fix:** Use normal Make object targets with explicit source/header dependencies and an explicit object list, separated by platform and build configuration. At minimum, stop each loop immediately on a nonzero compiler result, remove stale destination objects before compilation, and do not link through a broad wildcard. Publish the final plugin atomically after a successful link and dependency check.

**Acceptance tests:** A mock compiler fails the first and a middle dependency compile while later compilations would succeed. Seed old objects beforehand. Require the build to fail, the linker not to run, and the previous release artifact to remain unchanged. Repeat for each platform recipe.

**18. P2 — Thread-pool initialization can lose its shutdown wakeup after `pthread_create` fails.**

**Evidence:** Source review of [ctools_threads.c:230](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_threads.c:230). A hang was not induced in the running Stata process.

**Cause:** The initialization-failure branch writes `pool->shutdown` and broadcasts `work_available` without taking `queue_mutex`. Workers read the predicate under that mutex before waiting. A worker can observe false, miss the unsynchronized broadcast, and then sleep while the initializer waits forever in `pthread_join`. The unsynchronized predicate access is also a data race.

**Proposed fix:** Acquire `queue_mutex`, set shutdown, broadcast, and unlock before joining successfully created workers. Track their count explicitly and clear pool ownership fields after cleanup. Review all shutdown/error branches for the same predicate/condition-variable discipline.

**Acceptance tests:** Inject `pthread_create` failure after one or several successful creations. Use barriers to force the worker's check/wait interleaving, require bounded completion, and run a ThreadSanitizer harness where supported. Check repeat initialization after failure.

**19. P3 — Resetting the growing arena can orphan retained blocks.**

**Evidence:** `arena_reset.c` creates two blocks, resets the arena, and allocates enough to need the second block again. The original second block is no longer in the chain. No current production caller of `ctools_arena_reset()` was found, which limits present exposure.

**Cause:** [ctools_arena.c:195](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_arena.c:195) clears usage and sets `current=first`, retaining the chain. The next overflow allocation at [ctools_arena.c:104](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_arena.c:104) assigns a newly allocated block to `current->next`, overwriting the retained suffix. The aligned allocation path has the same structure. The orphaned blocks are never reached by `ctools_arena_free()`.

**Proposed fix:** Choose and document one reset contract. To retain capacity, walk/reuse existing successor blocks before appending. For a simpler implementation, free all blocks after the first during reset. Audit alignment-size arithmetic for overflow while changing this shared allocator. If reset is not intended to be supported, remove the unused public entry point rather than leaving an unsafe contract.

**Acceptance tests:** Repeated allocate/reset/free cycles with several block sizes and alignments, including a request larger than the normal block size. Track every allocation and free and verify bounded retained capacity.

**20. P3 — Consolidate wrapper infrastructure and remove ambiguous dormant routes.**

**Evidence:** Platform detection and plugin loading are copied across most command wrappers, and twice within import/export. Weight evaluation, sample staging, output-index mapping, and cleanup also vary across wrappers; items 2, 15, and 16 demonstrate consequences of that divergence. [cqreg.ado:543](/Users/Mike/Documents/GitHub/stata-ctools/build/cqreg.ado:543) selects native `qreg_p`, while a separate `cqreg_p.ado` remains in the tree. Its presence can mislead reviews about which prediction implementation is active; the active route passed the audit's factor-variable `stdp` comparison.

**Proposed fix:** Refactor in small steps after the correctness repairs. First introduce one plugin-loader helper that returns the selected binary and build information with a consistent failure policy. Next extract well-defined helpers for weight materialization and input/output index metadata, leaving command-specific sample rules in their wrappers. Standardize native status propagation and resource cleanup contracts. Audit the actual dispatch and `e(predict)` graph, then either test and wire dormant helpers or remove them from the distributable package. Avoid a wholesale rewrite of the statistical engines.

**Acceptance tests:** Run one shared loader test matrix for OS/architecture, missing binaries, incompatible binaries, and repeat loading. Verify every public command still has the intended entry point. Preserve command-level reference tests throughout each small refactor and check the package manifest whenever a helper is added or removed.

**21. P3 — Documentation overstates or misidentifies implemented behavior.**

**Evidence:** Concrete examples in the current tree:

- [README_cexport.md:29](/Users/Mike/Documents/GitHub/stata-ctools/docs/README_cexport.md:29) calls `datafmt` “not yet implemented,” but the wrapper does implement date/time formatting. That does not establish support for every numeric display format; the documentation should specify the actual subset.
- [README_cbinscatter.md:110](/Users/Mike/Documents/GitHub/stata-ctools/docs/README_cbinscatter.md:110) describes use of the same CG solver as `creghdfe`. The binscatter residualizer implements iterative projection sweeps rather than that shared CG implementation.
- [DEVELOPERS.md:242](/Users/Mike/Documents/GitHub/stata-ctools/DEVELOPERS.md:242) labels LSD the default, while the public `csort` wrapper defaults to `auto` at [csort.ado:34](/Users/Mike/Documents/GitHub/stata-ctools/build/csort.ado:34). Distinguish public defaults from any lower-level API convention.
- Broad replacement/compatibility language needs to reflect the limitations exposed above, especially merge key types, export metadata and dates, and binscatter control syntax. The existing compatibility matrix is a useful place to maintain those boundaries.

**Proposed fix:** Reconcile command help, per-command READMEs, developer documentation, and the compatibility matrix against actual parser options, dispatch routes, and executable examples. Mark support as implemented, intentionally unsupported, or known-broken with an issue reference; avoid describing partial support as complete. Generate mechanical option/default inventories where feasible and retain hand-written explanations for statistical semantics. Update documentation in the same change as each behavior repair.

**Acceptance tests:** Execute small documented examples offline using cached fixtures. Check public option/default inventories against wrappers and keep explicit tests for documented rejection behavior. Require each performance/algorithm description to identify the code path or benchmark supporting it.

**Suggested repair sequence.** Address merge schema/type safety and output transactions first (1, 2, 8, 15), then statistical correctness and failure propagation (3, 4, 7, 11, 14), followed by file fidelity (5, 6, 13), native sorting/thread reliability (9, 10, 12, 18), and wrapper/build contracts (16, 17). Complete the smaller API, refactoring, and documentation work after those changes have stable regression coverage. Promote the reproductions into the appropriate validation suites as the defects are fixed; in particular, test execution thresholds, alternate algorithms without wrapper fallbacks, no-OpenMP builds, allocation failures, and failure-time data preservation.
