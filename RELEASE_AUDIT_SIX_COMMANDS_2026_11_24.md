**Release audit of cdecode cencode crangestat csample cbsample and cdestring**

November 24, 2026. Repository base: `0388e3075907dc080e29849ffcf639d971b9ced8`.

Audit duration: approximately 25 minutes, within the requested 40-minute limit.

**Fix follow-up — November 24, 2026:** all 11 findings below have been addressed in the source and ado wrappers. The original findings and audit evidence are retained below as a historical record.

| Finding | Implemented fix | Regression coverage |
| --- | --- | --- |
| 1 | Transfer raw ignore bytes through a length-checked local macro, with a plugin acknowledgment before committing results. | Literal quotes in generate/replace/force; raw quote, whitespace and backslash bytes in the native fixture. |
| 2 | Exclude the current observation before SIMD summation; use a stable retained-window mean for higher moments. | Excluded large values at both ends of groups of 7, 8, 10, 63, 64 and 1,000 rows. |
| 3 | Sort eligible rows by cluster, or by strata and cluster, before evaluating the sample-size expression. | Subscripted expressions with and without strata and filtering. |
| 4 | Serialize scalar interval offsets with 17 significant digits. | Exact endpoints through both C count and Mata covariance paths. |
| 5 | Qualify plugin scalar reads with scalar() throughout all six wrappers. | Numeric and string variables named after plugin result/diagnostic scalars; encoding label round trips. |
| 6 | Give the shared parser a distinct memory-error status and abort destring before stores on conversion allocation failure. | Every allocation in long decimal and D-exponent conversions; serial/OpenMP builds. |
| 7 | Normalize positive missing-range values and IEEE overflow to unbounded endpoints; reset sliding searches when bounds move backward. Retain valid finite negative values, including values below -SV_missval. | Scalar bounds versus bounds generated in Stata, both signs and endpoints, small/large and single/multiple groups. |
| 8 | Parse and validate both bound-source indices before allocating or loading data. | Missing, negative, fractional, out-of-range and unrepresentable indices return 198 without allocation. |
| 9 | Preserve unweighted modern bootstrap calls, preallocate id/key variables, and reject expanded id-pass data above the plugin row limit. Restore on error or Break. | Real maxvar exhaustion; rollback after drop/expand/idpass/Break; every native id-pass allocation/read/store failure. Legacy version-7 behavior remains native. |
| 10 | Check all eleven range-statistic scalar publications, including thread diagnostics. | Each publication returns its injected error or Break status and frees allocations. |
| 11 | Aggregate converted/failed cell counts in uint64_t. | Synthetic per-column counts totaling 4,294,967,417, without large physical allocations. |

**Final validation:** a fresh Windows plugin built from the working tree passed **1,232 Stata checks, zero failures, zero skips**: cdecode 230, cencode 251, crangestat 338, csample 71, cbsample 109, cdestring 233. This adds 51 Stata checks to the original 1,181. Seven native harnesses passed in serial and OpenMP builds with undefined-behavior and float-cast-overflow instrumentation: labels, destring, group moments, counting sort, parallel ordering, range contracts, and bootstrap id pass. The final scoped diff check was clean.

The tested plugin is an isolated candidate at [ctools_windows.plugin](G:/GitHub/stata-ctools/temp/fix-six-20261124/candidate/ctools_windows.plugin); distribution binaries were not replaced. All platform binaries must be rebuilt with these wrappers before publishing. Linux/macOS execution and physical billion-cell workloads were not performed. Bootstrap rollback now incurs Stata preserve storage for unweighted calls; the help file documents this behavior.

Evidence: [six-command results](G:/GitHub/stata-ctools/temp/fix-six-20261124/suites-run.log), [full Stata log](G:/GitHub/stata-ctools/temp/fix-six-20261124/stata.log), [new Stata regressions](G:/GitHub/stata-ctools/validation/support/six_release_regressions.do), and [native regression entry point](G:/GitHub/stata-ctools/validation/native_six_release.py).

---

**Release recommendation:** fix findings 1–6 before publishing. They allow valid commands to produce incorrect results, lose values or labels, or mistake a resource failure for invalid data. Findings 7–10 address extreme numeric inputs and failure handling. Finding 11 is a large-data counter issue. There are **11 findings: 10 reproduced and one established by source inspection**. No production fixes were applied by this audit.

P1 denotes high-priority correctness or data loss; P2 denotes conditional correctness or failure-handling defects; P3 denotes a lower-priority diagnostic defect.

The audit covered the six ado wrappers, their C and Mata implementations, and shared parsing, label/hash, RNG, grouping, sorting, permutation, selection, data-transfer, and result-publication paths. Runtime checks used Windows and StataNow 19.5. Both the packaged Windows plugin and a fresh build from the working tree passed the existing six-command suite. Linux/macOS execution and physical billion-cell workloads were outside this timebox.

| Command | Existing checks passed per build | Findings |
| --- | ---: | --- |
| cdecode | 229 | No new command-specific correctness defect established |
| cencode | 248 | 5 |
| crangestat | 306 | 2, 4, 7, 8, 10 |
| csample | 70 | 5 |
| cbsample | 99 | 3, 9 |
| cdestring | 229 | 1, 6, 11 |
| Total | **1,181** | |

The passing suite does not cover the reproductions below. Native label/conversion, group/moment, counting/LSD-sort, and stable-order harnesses also passed in serial and OpenMP builds with undefined-behavior instrumentation. Targeted sampling allocation, store, and scalar-publication fault checks passed. The range and decimal-parser fault probes exposed findings 6, 8, and 10.

1. **[P1] cdestring loses valid quoted numbers when ignore() contains a literal double quote.**

   **Location:** [ignore argument construction](G:/GitHub/stata-ctools/build/cdestring.ado:183), [plugin call](G:/GitHub/stata-ctools/build/cdestring.ado:271), and the shared argument parser.

   **Reproduction and impact:** set a string variable to `char(34)+"12"+char(34)` and supply that quote character to `ignore()`. Native `destring` produces 12. `cdestring` reports nonnumeric characters and creates no output without `force`; with `replace force`, it returns success and replaces the valid value with missing. Its message even says the quote was removed, making the failure particularly misleading.

   **Cause:** the wrapper escapes backslashes and whitespace but interpolates a literal quote into the quoted plugin command argument. The ignore data is not transported intact across the Stata/plugin argument boundary.

   **Proposed fix:** pass the ignore bytes through a length-delimited local macro, as the label commands already do, or use a complete encoding that cannot be interpreted as argument syntax. Keep transport decoding separate from the actual ignore-character semantics.

   **Regression checks:** literal quotes alone and combined with other ignore characters; `generate()` and `replace`; with and without `force`; `aschars` and `asbytes`; fixed strings and strL fallback parity. See [reproduction and results](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes2.log).

2. **[P1] crangestat subtracts an excluded outlier after summation and returns incorrect sum, skewness, and kurtosis.**

   **Location:** [SIMD accumulation and subtraction](G:/GitHub/stata-ctools/src/crangestat/crangestat_kernels.c:1001), [skewness](G:/GitHub/stata-ctools/src/crangestat/crangestat_kernels.c:1137), and [kurtosis](G:/GitHub/stata-ctools/src/crangestat/crangestat_kernels.c:1166).

   **Reproduction and impact:** use ten observations with `x = 1e16` in the first row and `x = 1` in the other nine, with an unbounded interval and `excludeself`. For the first row, `rangestat` returns sum 9; `crangestat` returns 0. The remaining nine values are constant, so skewness and kurtosis should be missing, but `crangestat` returns 1 for both. The moment error also reproduced with groups of 64 and 10,000 observations.

   **Cause:** the simple kernel includes the excluded value in floating-point aggregates and subtracts it afterward. The small contributions have already been rounded away. Skewness and kurtosis then use the incorrect aggregate to compute their mean. The existing stable mean/variance paths do not protect these statistics.

   **Proposed fix:** omit self during accumulation, including the SIMD path, and use a stable mean for higher moments. If retaining the subtraction optimization, add a reliable cancellation test and recompute the window without self before using the result.

   **Regression checks:** sizes 7/8 and 63/64, larger groups, both signs of the excluded outlier, nonconstant remaining values, by-groups, and different thread counts. Verify both numeric values and missingness. See [sum results](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes2.log) and [moment results](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes3.log).

3. **[P1] cbsample evaluates clustered sample-size expressions in the wrong observation order.**

   **Location:** [_cbsample_size](G:/GitHub/stata-ctools/build/cbsample.ado:351), especially the sorts at lines 364 and 380.

   **Reproduction and impact:** create six rows with cluster values `6, 5, 4, 3, 2, 1`, then run `cbsample c[1], cluster(c) weight(w)`. Native `bsample` draws one cluster, while `cbsample` draws six: the sum of weights is 1 versus 6. This reproduces both without strata and with a constant `strata(s)` variable.

   **Cause:** native `bsample` sorts by sample membership, strata, and cluster before evaluating the expression. `_cbsample_size` receives no cluster argument and sorts only eligibility or strata. Expressions using subscripts therefore read a different row. This changes the requested number of draws, independently of the intentionally different random-number sequence.

   **Proposed fix:** pass the cluster variables into the helper and evaluate the expression after the same stable ordering as `bsample`. Preserve the correct unstratified versus within-stratum interpretation of subscripts and `_N`, then restore caller order as intended.

   **Regression checks:** `c[1]`, another variable's first/last value, multiple cluster variables, numeric and string clusters, `if/in`, and differing expressions across strata. See [probe and results](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes.log).

4. **[P2] crangestat rounds scalar interval endpoints before using them.**

   **Location:** [low endpoint conversion](G:/GitHub/stata-ctools/build/crangestat.ado:144), [high endpoint conversion](G:/GitHub/stata-ctools/build/crangestat.ado:157), and [argument construction](G:/GitHub/stata-ctools/build/crangestat.ado:552).

   **Reproduction and impact:** create two double keys, 0 and `1.0000000000000002`, and request a count with both interval offsets equal to `1.0000000000000002`. `rangestat` counts one match for the first row; `crangestat` returns missing. The same discrepancy reproduced with `.12345678901234567` and `123456789012345.67`.

   **Cause:** `local int_low_val = ...` and its high-bound equivalent use Stata's default numeric-to-macro formatting. The first example becomes the text `1`. Both C statistics and the Mata flexible-statistic path subsequently use the rounded macros.

   **Proposed fix:** preserve validated numeric tokens where possible, or serialize evaluated values with a round-trip-safe format such as `%24.17g`. Apply the same representation to both C and Mata paths.

   **Regression checks:** exact endpoint inclusion for values needing 17 significant digits, both interval sides, negative offsets, and mixed scalar/variable bounds; include a flexible statistic. See [boundary results](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes3.log) and [formatting demonstration](G:/GitHub/stata-ctools/temp/audit-six-20261124/sampling_extra.log).

5. **[P2] data variables can shadow plugin result scalars, losing cencode labels and breaking csample after deletion.**

   **Location:** [cencode category count](G:/GitHub/stata-ctools/build/cencode.ado:265), [new-label count](G:/GitHub/stata-ctools/build/cencode.ado:280), and [csample summary](G:/GitHub/stata-ctools/build/csample.ado:279).

   **Reproduction and impact:** encode strings `a` and `b` with an existing numeric variable `_cencode_n_unique = 0`. The command succeeds and writes codes 1 and 2, but the value-label definition is absent. With `replace`, the source strings disappear and a subsequent `decode` fails. Other shadow values cause r(610), and `_cencode_n_new = 0` prevents new labels from being added to an existing definition.

   In `csample`, `_csample_dropped = 999` makes a ten-row, 50-percent sample report 999 deletions. Making that variable a string containing `x` produces r(111) after the dataset has already shrunk to five rows.

   **Cause:** unqualified expressions resolve dataset variables before same-named scalars. Some neighboring reads already use `scalar(...)`, but these do not.

   **Proposed fix:** qualify every plugin scalar read with `scalar(...)`, including diagnostic expressions across the six wrappers. Prefer uniquely named temporary scalars or local results for future interfaces.

   **Regression checks:** exact-name and abbreviation collisions, numeric and string variables, `generate()`/`replace`, and label extension. Assert label contents and post-command return codes, not just numeric codes. See [label probes](G:/GitHub/stata-ctools/temp/audit-six-20261124/collisions.log) and [replacement and sampling probes](G:/GitHub/stata-ctools/temp/audit-six-20261124/collisions2.log).

6. **[P2] the shared decimal parser turns allocation failure into a nonnumeric value.**

   **Location:** [long-input allocation](G:/GitHub/stata-ctools/src/ctools_types.c:678) and [cdestring parser consumer](G:/GitHub/stata-ctools/src/cdestring/cdestring_impl.c:288).

   **Reproduction and impact:** a valid 142-character decimal, `1.` followed by 140 ones, normally converts to approximately 1.1111111111111112. In an isolated native harness, failing the parser's temporary-buffer allocation makes `cdestring_main(... force ...)` return 0, store missing, and report one nonnumeric value. Other tested allocation failures return 920. Native `destring` and the normal ctools execution agree on the valid number.

   **Cause:** `ctools_parse_double_fast()` returns the same false status for malformed text and failed `malloc()`. The converter has no way to distinguish the failure and applies ordinary `force` semantics.

   **Proposed fix:** use separate success, invalid-input, and allocation-error statuses and propagate 920 before committing output. Alternatively, eliminate this allocation with a suitably bounded reusable buffer for the conversion path. Review other users of the shared parser when changing its contract.

   **Regression checks:** inject failure at the long/ambiguous decimal fallback, with and without `force`, and ensure replacement preserves the original variable on resource failure. See [fault results](G:/GitHub/stata-ctools/temp/audit-six-20261124/destring_oom.log) and [reproducible harness runner](G:/GitHub/stata-ctools/temp/audit-six-20261124/faults.py).

7. **[P2] crangestat handles overflowing scalar interval bounds differently from rangestat and its own variable-bound path.**

   **Location:** [cross-group bounds](G:/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:709), [run bounds](G:/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:972), and [sliding-loop bounds](G:/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:1112).

   **Reproduction and impact:** use keys `-8e307, 0, 8e307` and count with `interval(key 1e307 .)`. On the last row, `rangestat` returns 3 while `crangestat` returns missing. Creating a double variable `low = key + 1e307` and supplying `interval(key low .)` makes `crangestat` return 3, demonstrating inconsistent equivalent paths.

   **Cause:** Stata converts the overflowing addition into missing, and a missing interval bound means unbounded. The C path instead retains a finite value above Stata's usable range, or an IEEE infinity, and binary-searches with it.

   **Proposed fix:** centralize bound construction and normalize overflow according to Stata arithmetic before applying missing-bound semantics. Apply this consistently to all optimized loops, including their initial sliding pointers.

   **Regression checks:** overflow in either direction, each endpoint separately, repeated-key runs, by-groups, and equivalent precomputed variable bounds. See [comparison](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes3.log).

8. **[P2] malformed crangestat variable-bound indices can crash the native routine.**

   **Location:** [unchecked bound-index parsing](G:/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:40) and [source-array indexing](G:/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:704).

   **Reproduction and impact:** an isolated native harness using an otherwise valid request fails with an access violation or undefined-behavior trap when given `varbounds` without indices, `varbounds lowsrc=100 highsrc=0`, or `varbounds lowsrc=1e100 highsrc=0`. These are malformed internal plugin requests; ordinary current ado syntax supplies valid indices.

   **Cause:** the indices are parsed as doubles, cast to `int` without range/integrality checks, and used without checking `0 <= index < nsource`. Missing indices default to -1. This leaves the plugin vulnerable to a bad direct call or a wrapper/protocol mismatch.

   **Proposed fix:** parse with the checked integer parser, require both indices under `varbounds`, validate against `nsource`, and return 198 before loading or indexing data.

   **Regression checks:** omitted, negative, equal-to-nsource, fractional, huge, NaN, and duplicate options. Run these in isolated native tests; each should return 198 without writes or leaked allocations. See [fault log](G:/GitHub/stata-ctools/temp/audit-six-20261124/range_faults.log).

9. **[P2] cbsample can fail after destructive sampling has already changed the dataset.**

   **Location:** [drop and expansion](G:/GitHub/stata-ctools/build/cbsample.ado:276) followed by [output creation and idcluster pass](G:/GitHub/stata-ctools/build/cbsample.ado:281).

   **Reproduction and impact:** with `set maxvar 2048`, create 12 rows and 2,042 variables including a cluster identifier. `cbsample 6, cluster(c) idcluster(newid)` returns r(900), but leaves six rows and a new public `newid` variable. The missing capacity is discovered only after dropping and expanding observations. Later id-pass allocation or SPI errors also occur after the destructive phase.

   Native `bsample` also changes the dataset before failing on this fixture. This finding is a failure-safety improvement, rather than a compatibility regression.

   **Proposed fix:** allocate required output and key variables before dropping rows. Protect the destructive phase with rollback through successful id creation and final ordering; preflight the expanded observation count against the plugin limit before entering the id pass.

   **Regression checks:** maxvar exhaustion, expansion failure, id-pass allocation/store errors, and interruption. Assert original observation count, values, order, and public variables after failure. See [native and ctools comparison](G:/GitHub/stata-ctools/temp/audit-six-20261124/probes2.log).

10. **[P2] crangestat ignores errors while publishing scalar results.**

    **Location:** [scalar publication](G:/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:1274) and the following `CTOOLS_SAVE_THREAD_INFO` call.

    **Reproduction and impact:** make the mocked Stata scalar-save callback return 459. `crangestat_main()` still returns 0 after storing the result columns. Failed publication can therefore leave absent or stale diagnostics while reporting success, and the ado wrapper does not receive the failure needed to initiate its normal rollback.

    **Cause:** nine `SF_scal_save` return values, plus thread-information publication, are discarded. The sibling encode, destring, and sampling commands check these failures.

    **Proposed fix:** check and propagate every publication return code, including thread diagnostics. Keep cleanup independent of publication success so the wrapper can roll back on the reported failure.

    **Regression checks:** fail each published scalar individually, including a Break-style return code, and verify propagation and wrapper rollback in normal non-fast mode. See [fault log](G:/GitHub/stata-ctools/temp/audit-six-20261124/range_faults.log).

11. **[P3] cdestring aggregate conversion counters overflow above 2,147,483,647 cells.**

    **Location:** [32-bit total declarations](G:/GitHub/stata-ctools/src/cdestring/cdestring_impl.c:469) and [cross-variable aggregation](G:/GitHub/stata-ctools/src/cdestring/cdestring_impl.c:671).

    **Evidence and impact:** source inspection establishes that `total_converted` and `total_failed` are signed `int`, although they sum counts across all variables. The per-variable observation limit does not bound this total: 100 variables with 22 million converted rows produce 2.2 billion conversions. Signed overflow is undefined behavior and can produce invalid published totals. No billion-cell dataset was allocated during this audit; the identified direct effect is on counts and diagnostics.

    **Proposed fix:** use `uint64_t` or `size_t` for aggregate counters and convert to double only at the scalar-publication boundary. Keep intermediate additions in the wider type.

    **Regression checks:** test aggregation with synthetic per-variable counts immediately below and above INT_MAX under undefined-behavior instrumentation; a giant Stata dataset is unnecessary.

**Evidence and follow-through**

The [audit evidence directory](G:/GitHub/stata-ctools/temp/audit-six-20261124) contains Stata reproductions, logs, native fault harnesses, and the isolated candidate build. The [suite log](G:/GitHub/stata-ctools/temp/audit-six-20261124/stata.log) records the passing fresh-build checks. All targeted correctness reproductions described above were confirmed with the fresh candidate, so they are not explained by the packaged binary being older than the checkout.

Promote the focused reproductions into the corresponding Stata/native suites when implementing fixes, then rerun the six-command suite against the rebuilt release artifacts. Retain exact assertions for counts, endpoint membership, labels, and missingness; broad relative-error tolerances would miss several of these failures. Recheck the numerical fixes on Linux and macOS as well as Windows before distributing their binaries.
