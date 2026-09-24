# ctools audit — September 22, 2026

**Repair status:** All 16 findings below have been addressed in the working tree.
See [SEP22_ASTRA_FIXES.md](/Users/Mike/Documents/GitHub/stata-ctools/SEP22_ASTRA_FIXES.md)
for the changes, explicit compatibility boundaries, and the passing 2,866-check
full-suite result. The assessment and reproductions below describe the original
pre-repair snapshot and are retained as historical evidence.

## Assessment

**Do not release the examined tree as a validated replacement for the reference commands.** The most urgent findings are silent row loss on invalid merge cardinalities, incorrect destination indices in winsorization, loss of cluster identity in three estimators, incorrect estimation-sample markers, and an out-of-bounds read in PPML after separation. These are independent of the limitations already described in the new compatibility documentation.

This report concerns the **working tree**, including its existing uncommitted repairs, rather than only commit `2518a9cdf7e4e9b91ae0f7ffaabee43286325da2` on `dev`. It is a fresh audit, not a restatement of `CTOOLS_AUDIT_ASTRA_SEP20.md`. Findings below distinguish runtime reproductions, source-confirmed defects, documentation discrepancies, and unverified risks. No implementation fixes were made to the repository for this audit.

The tree was changing while the audit ran. Compatibility/platform documents, a compiler/dependency test, a corrected `noextend` fixture, and a repaired environment check appeared during the investigation. I rechecked them and **do not list the earlier missing-file, fixture, or environment-check failures as outstanding findings**. The final C source snapshot and file hashes are preserved with the temporary audit evidence.

### Priorities

- **P1 — fix before release:** data corruption, memory safety, incorrect statistical results, or failure of a primary advertised operation.
- **P2 — fix before claiming compatibility:** narrower behavioral errors, incomplete distribution, missing validation coverage, or misleading executable documentation.
- **P3 — maintenance:** development instructions and structural improvements that make further defects more likely.

No P0 issue was established. Within each priority, the ordering below reflects impact and breadth.

| ID | Priority | Finding | Evidence |
|---|---|---|---|
| S01 | P1 | `cmerge` accepts duplicate keys and silently discards rows | Reproduced after the late wrapper fix |
| S02 | P1 | `cwinsor` supplies dataset indices to a plugin-local checked store API | Reproduced |
| S03 | P1 | OLS and PPML conflate distinct clusters through integer casts | Reproduced; standard errors change |
| S04 | P1 | IV regression also truncates numeric cluster IDs | Reproduced; valid relabelings fail |
| S05 | P1 | OLS/PPML `e(sample)` includes dropped singleton/separated observations | Reproduced against references |
| S06 | P1 | PPML reads beyond a shortened mask and misaligns cluster data after separation | Source-confirmed |
| S07 | P1 | `cdecode` silently loses or changes valid label text | Reproduced |
| S08 | P1 | `cbsample` ignores requested cluster draws and changes stratum sample sizes | Reproduced against native `bsample` |
| S09 | P2 | PPML falls through to generic prediction without defining its meaning | Reproduced; partial linear predictor returned |
| S10 | P2 | OLS/PPML reject valid weight expressions | Reproduced |
| S11 | P2 | Matching reports requested, rather than actual, neighbor counts | Reproduced with ties |
| S12 | P2 | Empty-dataset merges ignore `keep()` when `nogenerate` is used | Reproduced |
| S13 | P2 | The local package is not a complete installable distribution | Manifest check fails |
| S14 | P2 | The validation gate leaves substantial behavior untested | Full-suite evidence and CI inspection |
| S15 | P2 | Installed help advertises invalid or unimplemented options | Syntax/source checks; PPML options reproduced |
| S16 | P3 | Developer instructions describe an obsolete build workflow | Documentation/source comparison |

**Finalization status:** the verified full-suite snapshot predates a few late C edits. The new bounded-text strL read path and import macro-buffer adjustments were not runtime-verified here. The late merge-wrapper correction was separately tested and its remaining cardinality defect is S01. This report identifies the evidence boundary instead of treating a moving checkout as one immutable build.

## Scope, method, and validation

I inspected the ado interfaces, plugin dispatcher, shared I/O and grouping infrastructure, C command implementations, package manifest, build scripts, CI workflow, command help, README files, and the earlier audit. The repository contains roughly 123,000 lines across the inspected C/header/ado/do-file inventory; this was a risk-directed audit, not a claim that every line was exhaustively verified.

Runtime checks used StataNow/MP 18.5, revision 26 February 2025, on Apple Silicon. Every Stata invocation went through `/bin/zsh -lic 'oldstata ...'`; the sandboxed attempt failed in Shortcuts, and subsequent exact audit commands ran outside the sandbox. Stata exited cleanly with `exit, clear`, allowing the wrapper to restore network time. **The July 1 timestamps in Stata logs come from that wrapper; they are not the date of this audit.** All Git commands used the required external zsh workflow.

A fresh plugin was built in `/tmp/ctools-sep22-audit/build`, with a final build from a frozen copy of the C sources. The build used the existing deployment-compatible static OpenMP 21.1.8 runtime at `/tmp/ctools-p2-fixes/libomp-arm21`. Both source builds completed; the final plugin passed the macOS architecture/dependency contract check. The build emitted **61 warnings**, largely unused functions/variables. The plugin identified itself as version `1.0.2`, build `sep22-audit`. Repository binaries were not overwritten.

The isolated audit directory contains copies of ado/help/validation files. Synthetic reproductions used local, disposable data. No reference packages were installed or updated, test tolerances were not relaxed, and repository validation tests were not rewritten by this audit.

### Results

- Native P1 checks: **4 passed** — compensated arithmetic and allocation-failure cleanup.
- Final native P2 checks: **6 passed** — signed label lookup, encoding BOMs, shared SPI failures, CSV/XLSX write failures, and matrix-store failures.
- Compiler/dependency-contract checks: **2 groups passed**.
- Release identity synchronization: **passed**.
- Package metadata-only validation: **passed**. Complete package validation: **failed**, with three missing binaries; see S13.
- Updated targeted P2 Stata suite: **14 passed, 0 failed**.
- Late merge-wrapper check: basic merge passed; targeted P1 suite **10 passed, 0 failed**. Duplicate-key and empty-side filtering defects remained reproducible.
- Initial full Stata suite: **1,204 passed, 29 failed, 70 skipped**. A later full run against the frozen final source and corrected fixture is recorded below.

The full rerun against the frozen source snapshot completed with **1,206 passed, 28 failed, 70 skipped** and return code 1. The updated P2 suite accounted for the improvement. Component totals were:

| Component | Passed | Failed | Skipped / limitation |
|---|---:|---:|---|
| csort | 88 | 1 | Script aborted, r(679) |
| cmerge | 0 | 16 | Ordinary merge error 2 in audited wrapper |
| cimport | 313 | 0 | 22 skipped |
| cexport | 16 | 1 | Script aborted, r(679) |
| creghdfe | 17 | 1 | Script aborted, r(679) |
| cqreg | 33 | 1 | Script aborted, r(679) |
| civreghdfe | 16 | 1 | Script aborted, r(679) |
| cdecode | 74 | 1 | Script aborted, r(679) |
| cencode | 74 | 1 | Script aborted, r(679) |
| cdestring | 166 | 1 | strL comparison in audited source |
| csample | 27 | 0 | |
| cbsample | 29 | 0 | |
| cbinscatter | 84 | 1 | Script aborted, r(679) |
| cpsmatch | 67 | 0 | 48 skipped |
| crangestat | 91 | 1 | Script aborted, r(679) |
| cwinsor | 0 | 1 | Store error 459 |
| cpplmhdfe | 88 | 0 | |
| audit_p1 | 9 | 1 | Merge error 2 in audited wrapper |
| audit_p2 | 14 | 0 | Corrected fixture |

**Snapshot boundary:** additional working-tree edits arrived after this full run, including the missing merge-varlist correction, bounded text-strL reads, and import macro-buffer bounds. The separate late-merge verification described under S01 covers the wrapper fix and new cardinality checks; the new strL implementation and these later C changes were not included in the full-suite totals. Those totals must not be represented as a run of the subsequently changing checkout.

The initial suite's 29 failures were not 29 independent implementation bugs: 16 came from the subsequently repaired missing merge varlist discussed under S01; the P1 merge regression was another; one was the subsequently corrected `noextend` fixture; one was the documented `strL` limitation; one was S02; and nine component scripts aborted with `r(679)`. A direct `webuse nlswork, clear` returned `web error 679` on this machine. Do not count the unexecuted portions of those component scripts as passing.

PPML's existing suite passed all **88 checks** in the initial full run. That is positive evidence for the earlier covariance repairs, but does not establish correct sample markers or safe clustering after separation. The current bootstrap suite also passed all **29 checks**, despite S08; its coverage does not establish the requested draw-count semantics.

## Detailed findings

### S01 — P1: `cmerge` does not enforce key cardinality and silently loses duplicate rows

**Location:** [src/cmerge/cmerge_join.c:211](/Users/Mike/Documents/GitHub/stata-ctools/src/cmerge/cmerge_join.c:211), particularly lines 217–225 and 235–257; [build/cmerge.ado:910](/Users/Mike/Documents/GitHub/stata-ctools/build/cmerge.ado:910).

**Current remaining defect:** the join kernel calculates duplicate-group sizes but does not validate them against the requested merge type. `MERGE_1_1` emits only the first master/using pair; `MERGE_M_1` uses only the first using row; `MERGE_1_M` uses only the first master row. Violating uniqueness on a declared “1” side is therefore accepted and can discard data.

**Reproduction:** save a using dataset with IDs 1 and 2. In the master, create two observations with ID 1 and distinct row markers 1 and 2. Native `merge 1:1` rejects the duplicate master key with `r(459)`. The late-fixed `cmerge 1:1` returns success, emits the master row with marker 1 plus the using-only ID 2 row, and silently discards the master row with marker 2. The exact native/ctools comparison and retained markers are recorded in `cardinality.log`.

**Proposed fix:** validate uniqueness on every side declared unique, before modifying the master. This includes duplicate groups that do not match the other dataset. Return the native-compatible error with a useful side/key diagnostic. Do not silently convert a `1:1`, `m:1`, or `1:m` request into a different relation. Keep all input/output staging transactional.

**Required regression checks:** duplicates in master only, using only, both, and unmatched groups for every merge type; distinct row markers to detect loss; missing-value and string keys; and unchanged observation counts/data after rejection. A successful valid `1:1` test cannot establish this contract.

**Repair observed during this audit:** the earlier wrapper used an undefined local `current_varlist` for its phase-two plugin call. With the new shared index checks, an ordinary valid two-row merge returned `r(2)` after expanding the master to four observations. The full-suite snapshot recorded 16 failures from that problem. A late edit changed the call to the actual `__all_phase2_vars` list. I verified that the same valid merge now succeeds with two observations, and the targeted P1 suite now passes **10/10**. The undefined-varlist defect is therefore closed, not the outstanding basis for this finding. The failure branch still contains no general restoration of an expanded master; error rollback should be covered in the cardinality repair.

### S02 — P1: `cwinsor` uses the wrong index space for storing results

**Location:** [build/cwinsor.ado:163](/Users/Mike/Documents/GitHub/stata-ctools/build/cwinsor.ado:163), lines 163–200; [src/cwinsor/cwinsor_impl.c:675](/Users/Mike/Documents/GitHub/stata-ctools/src/cwinsor/cwinsor_impl.c:675); [src/ctools_data_io.c:1597](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_data_io.c:1597).

The ado computes destination positions in the entire Stata dataset, then calls the plugin with only `target_varlist by`. Its comment says that stores use global dataset positions. The current shared store functions validate and write through checked callbacks in the **plugin-local** variable space. These contracts conflict.

**Reproduced:**

```stata
sysuse auto, clear
capture cwinsor price, suffix(_w)   // r(459)
capture cwinsor price, replace      // r(459)
```

With a dataset containing only `x` as its first variable, `cwinsor x, replace` succeeds. Thus dataset column placement changes whether the same operation works. The full winsor validation aborts at its first plugin smoke check.

**Proposed fix:** use destination indices `1..nvars_target` for the current selective plugin varlist, or explicitly pass a complete list and calculate every index from it. Keep the checked callback boundary. Ensure a failed generated-output operation removes incomplete outputs or restores the previous dataset.

**Required regression checks:** move an identical target among first/middle/last columns; test generated outputs, multiple targets, mixed string/numeric grouping, `if/in`, and `replace`. Assert that read-only grouping and unrelated columns remain unchanged.

### S03 — P1: OLS and PPML incorrectly recognize different cluster variables as fixed effects

**Location:** [src/creghdfe/creghdfe_regress.c:342](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c:342), lines 354–370 and 1336 onward; [src/cpplmhdfe/cpplmhdfe_irls.c:286](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:286), lines 296–305 and 1469 onward.

An optimization checks whether cluster values equal an absorbed FE variable, casts both values to `ST_int`, and, on a match, substitutes the FE's group IDs. Different numeric clusters with the same integer part therefore become the same cluster. Large numeric labels also make these casts unsafe outside the integer range.

**Reproduction design:** 400 observations, `g = ceil(_n/20)`, and `cl = g + cond(mod(_n,2),.1,.2)`. There are 20 fixed-effect groups and 40 clusters. `egen cid = group(cl)` preserves the cluster partition while giving it integer labels. Fit the same model twice, using `vce(cluster cl)` and `vce(cluster cid)`.

| Estimator | Clusters using `cl` | Clusters using `cid` | SE of `x`, `cl` | SE of `x`, `cid` |
|---|---:|---:|---:|---:|
| `creghdfe` | 20 | 40 | 0.04892324 | 0.05085301 |
| `cpplmhdfe` | 20 | 40 | 0.04028746 | 0.04398411 |

The outcomes were generated with seed 45678 as recorded in `repro.do`. Relabeling groups must not change inference.

**Proposed fix:** compare original doubles without conversion for an exact-value shortcut, or compare canonical group partitions. Use one tested, type-preserving remapper for general cluster IDs. Derive FE nesting from partitions, not numeric-code resemblance.

**Required regression checks:** nonconsecutive integers, fractional labels, negative labels, labels above 32-bit range, string labels, and clusters nested in/crossing FEs. Compare complete covariance matrices, cluster counts, and degrees-of-freedom adjustments under bijective relabelings.

### S04 — P1: IV regression truncates numeric clusters before remapping them

**Location:** [src/civreghdfe/civreghdfe_impl.c:201](/Users/Mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:201), especially lines 216 and 233; subsequent cluster remapping at lines 917–945.

Both numeric cluster dimensions are cast directly from `double` to `ST_int`. Later remapping cannot recover groups that were already combined by truncation. String clusters take a different, identity-preserving path.

**Reproduced:** generate 20 groups as `g = ceil(_n/20)` in a 400-observation IV model, then relabel them as `fraction = g/100` and `large = 1e12+g`. The model with `vce(cluster g)` succeeds with 20 clusters and an `x` standard error of 0.04430205. Both relabelings fail with `r(504)` and `invsym(): matrix has missing values`. The fractional version has been reduced to a single integer code before remapping.

**Proposed fix:** remap numeric doubles directly to dense IDs before storing anything in an integer array. Keep missing-value status separate from a valid negative group code; apply the same rules to both cluster dimensions.

**Required regression checks:** one-way and two-way clustering under equivalent string, fractional, negative, and large-number labels. Reject genuinely insufficient cluster counts with an explicit diagnostic rather than allowing an invalid covariance matrix to reach Mata.

### S05 — P1: `e(sample)` does not describe the observations used by OLS or PPML

**Location:** [build/creghdfe.ado:350](/Users/Mike/Documents/GitHub/stata-ctools/build/creghdfe.ado:350) and line 716; [build/cpplmhdfe.ado:314](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe.ado:314) and line 540.

The main regression plugin calls do not update the ado's original `touse` variable to reflect observations dropped internally. The ado posts that original marker as `e(sample)` while reporting the reduced `N_final` as `e(N)`.

**Reproduced:**

| Case | Reference `e(N)` / sample count | ctools `e(N)` / sample count |
|---|---:|---:|
| OLS, one singleton among 401 observations | 400 / 400 | 400 / 401 |
| PPML, one singleton among 401 observations | 400 / 400 | 400 / 401 |
| PPML, one all-zero FE group among 400 observations | 380 / 380 | 380 / 400 |

Native references were installed `reghdfe` and `ppmlhdfe`. The ctools commands in this reproduction return success. This affects any subsequent summary, sample comparison, prediction restriction, or estimation that uses `if e(sample)`.

**Proposed fix:** maintain an original-observation map through every removal/compaction step. Return the final estimation indicator through an explicitly passed output variable and post that indicator. The recent IV residual/sample repair provides a relevant pattern, but OLS and PPML need their own removal-path tests.

**Required regression checks:** `count if e(sample) == e(N)`; singleton chains across multiple FEs; separation; missing regressors/weights; `if/in`; and consistency among sample markers, saved residuals, and stored/restored estimates. An exploratory OLS call with saved residuals also returned `esample() invalid` (`r(471)`); the confirmed marker mismatch above does not depend on that additional behavior.

### S06 — P1: PPML separation shortens `mask`, but clustered VCE still indexes its original length

**Location:** [src/cpplmhdfe/cpplmhdfe_irls.c:692](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:692), especially lines 715–721; [same file:1488](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:1488).

**Source-confirmed memory defect; no sanitizer-confirmed Stata crash is claimed.** After FE separation drops rows, the code frees the original mask, allocates `N_new` entries, and fills them with ones. For a cluster variable that does not match an FE, the later VCE path loops to **`N_orig`** and reads `mask[i]`. Because separation makes `N_new < N_orig`, that loop reads past the allocation. The `idx < N` guard protects a destination write, not the preceding out-of-bounds mask read.

There is also an alignment error before the out-of-bounds portion: the new mask no longer identifies the retained original rows, but the code uses it to select from `cluster_raw_values` in the original row order. If removed rows occur near the start or middle, retained outcomes can receive other rows' cluster labels.

**Trigger:** PPML with at least one FE group whose outcomes are all zero, plus `vce(cluster cl)` where `cl` is not recognized as an absorbed FE. Filtering or singleton removal can make the mapping still more complicated.

**Proposed fix:** compact cluster arrays alongside outcomes, regressors, weights, offsets, and FE IDs using a single retained-row map. Keep separate, explicitly named original-sample and compact-sample indices. Do not reuse an all-ones compact mask as an original-row selector. Treat allocation/remapping failure as an error before any partially compacted state is used.

**Required regression checks:** AddressSanitizer coverage for separation plus independent clustering; removed blocks at the start/middle/end; `if/in` and prior singleton removal; equivalence to manually dropping separated observations before estimation. Existing unclustered or FE-clustered PPML parity tests do not exercise the faulty branch.

### S07 — P1: `cdecode` can return success with missing or altered label text

**Location:** [src/ctools_hash.c:621](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_hash.c:621), the `label save` parser and fixed line buffers through line 755; [src/ctools_hash.h:169](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_hash.h:169); [src/cdecode/cdecode_impl.c:207](/Users/Mike/Documents/GitHub/stata-ctools/src/cdecode/cdecode_impl.c:207); [build/cdecode.ado:193](/Users/Mike/Documents/GitHub/stata-ctools/build/cdecode.ado:193).

Three independent cases were reproduced:

1. **Labeled extended missing values:** `label define lab 1 "one" .a "refused"` followed by native `decode` returns `"refused"` for `.a`. `cdecode` returns an empty string and `_rc == 0`. The C loop skips missing values before attempting label lookup, and its integer-key parser does not represent `.a`.
2. **Literal text:** labels containing a backtick/local-macro delimiter, compound quotes, or embedded tab/newline/carriage-return characters do not round-trip. In a seven-row comparison, rows 2, 3, and 4 failed exact equality; the control-character row became empty. Dollar-sign text, Unicode, ordinary quotes, and a 2,045-byte simple label passed in this fixture.
3. **Long labels:** a 3,000-character label created with `st_vlmodify()` decodes to 3,000 characters with native `decode`, but to **zero characters** with `cdecode`, again with `_rc == 0`. The fixed-width line reader can fail to find the closing quote; the wrapper silently caps its destination at 2,045 bytes.

The decode loop also ignores both `SF_sstore` return values and increments the decoded count after an attempted store. Thus repairing shared I/O alone does not make this command report its own write failures.

**Proposed fix:** read value-label data through a representation that preserves literal bytes and extended missing codes, preferably a structured Stata/Mata interface instead of reparsing executable `label save` text. Support long-string output where possible; otherwise reject unsupported label lengths explicitly before replacing anything. Implement `maxlength()` by deliberately truncating text according to its documented contract, rather than relying on a failed destination write. Check every store and retain the original numeric variable on error.

**Required regression checks:** all extended missings, signed numeric codes, nested quotes, backslashes, macro delimiters, multiline labels, Unicode byte boundaries, 2,045/2,046/3,000-byte labels, explicit `maxlength`, and output-store fault injection. Test both `generate()` and `replace`.

### S08 — P1: bootstrap sample-size semantics differ from native `bsample`

**Location:** [src/cbsample/cbsample_impl.c:628](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c:628), lines 631–642 and 656–668; [build/cbsample.ado:59](/Users/Mike/Documents/GitHub/stata-ctools/build/cbsample.ado:59).

The cluster branch calculates `target_n` but ignores it, drawing `n_clusters_in_stratum` clusters regardless of the requested sample size. The observation branch allocates the requested count proportionally across strata and forces a minimum of one. Native `bsample` interprets an explicit count as the number of observations, or clusters, **within each stratum**.

**Reproduced on 40 observations in four equal groups of ten:**

| Command | Resulting observations |
|---|---:|
| `bsample 1, cluster(cl)` | 10 |
| `cbsample 1, cluster(cl)` | 40 |
| `bsample 2, strata(cl)` | 8 |
| `cbsample 2, strata(cl)` | 4 |

These differences change the bootstrap experiment, not merely the random seed or selected IDs. Existing passing tests do not cover these cases.

**Proposed fix:** preserve whether a positional count was supplied. Use native per-stratum draw-count semantics, with cluster counts rather than row counts when clustering. Define defaults from each stratum's available sampling units; validate oversized requests consistently with native behavior. Remove the unused calculation and minimum-one workaround after the contract is corrected.

**Required regression checks:** explicit/default counts; unequal strata; unequal cluster sizes; cluster plus strata; one requested draw; frequency-weight output; repeated cluster labels across strata; and invalid counts. Compare numbers of draws and expanded row multiplicities, not identical random selections.

### S09 — P2: PPML has no predictor and silently falls back to generic linear prediction

**Location:** [build/cpplmhdfe.ado:540](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe.ado:540), ereturn block through line 577; [docs/COMPATIBILITY.md](/Users/Mike/Documents/GitHub/stata-ctools/docs/COMPATIBILITY.md).

The new compatibility table correctly says that no dedicated PPML predictor is provided. However, the estimator does not register a predictor that enforces this limitation. `predict double mu` succeeds using generic prediction, returning only the available slope index. In the reproduction it produced values as low as **−0.64638136**. This is not a Poisson conditional mean and omits absorbed effects and any offset/exposure contribution.

This is a postestimation-interface issue, not evidence that the estimated PPML coefficients are wrong. The explicit documentation caveat reduces, but does not remove, the risk of a plausible command silently producing the wrong intended object.

**Proposed fix:** register a dedicated predictor. Until means/FE reconstruction are implemented, reject default prediction and unsupported statistics with a clear error, and optionally support explicit `xb` with a precise definition. A full predictor must distinguish slope index, index including FEs/offset, and `exp(index)`.

**Required regression checks:** default prediction, explicit supported statistics, offsets/exposure, saved/restored estimates, and nonnegative fitted means compared with `ppmlhdfe` when means are supported.

### S10 — P2: advertised weight support is restricted to variable names

**Location:** [build/creghdfe.ado:77](/Users/Mike/Documents/GitHub/stata-ctools/build/creghdfe.ado:77); [build/cpplmhdfe.ado:72](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe.ado:72). IV's handling at [build/civreghdfe.ado:467](/Users/Mike/Documents/GitHub/stata-ctools/build/civreghdfe.ado:467) should be included in the repair review.

The OLS and PPML wrappers strip the equals sign from the parsed weight expression and treat the remaining text as a variable name. Thus `[aw=2*w]` reaches `markout` as the nonexistent variable `2*w`.

**Reproduced:** both `creghdfe y x [aw=2*w], absorb(g)` and `cpplmhdfe y x [aw=2*w], absorb(g)` return `r(111)` and `variable 2*w not found` on otherwise valid synthetic data. The syntax/help advertise Stata weight types without this expression restriction.

**Proposed fix:** evaluate the parsed expression once into a temporary double variable over the marked sample. Apply missing/positivity/frequency-integrality checks to its values. Pass only that variable to C. Alternatively, explicitly document and validate a narrower contract, though that would retain a compatibility gap.

**Required regression checks:** expressions versus precomputed equivalent variables, constants, arithmetic with parentheses, missing/zero/negative weights, nonintegral fweights, and all advertised weight types.

### S11 — P2: matching's `_nn` is not the number of neighbors actually used

**Location:** [build/cpsmatch.ado:272](/Users/Mike/Documents/GitHub/stata-ctools/build/cpsmatch.ado:272); [build/cpsmatch.sthlp:88](/Users/Mike/Documents/GitHub/stata-ctools/build/cpsmatch.sthlp:88).

After matching, the ado sets `_nn = neighbor` for every supported treated observation. That is the requested count. Ties can increase the actual count; calipers or a small control pool can reduce it. Radius and kernel methods have their own actual support sizes.

**Reproduced:** two controls and one treated observation all have propensity score `.5`. `cpsmatch treat, pscore(ps) outcome(y) ties` gives the two controls weights `.5` each, showing both are used, but `_nn` for the treated observation is **1**. Help defines `_nn` as the number of neighbors used.

**Proposed fix:** return each treated observation's actual match count from C and store it directly. Document `_nn` for radius/kernel matching and unmatched/control rows. Do not infer it from aggregate matching weights.

**Required regression checks:** ties at the cutoff, insufficient eligible controls, calipers, radius/kernel support, and without-replacement exhaustion. ATT may remain correct while `_nn` is wrong; test the output variables separately.

### S12 — P2: empty-dataset merges ignore `keep()` with `nogenerate`

**Location:** [build/cmerge.ado:588](/Users/Mike/Documents/GitHub/stata-ctools/build/cmerge.ado:588), especially lines 636–666.

The empty-side path remains faulty after the late phase-two wrapper repair. It applies `keep()` only when `nogenerate` is absent. With an empty using dataset, `cmerge 1:1 id using ..., keep(match) nogenerate` returned success and retained both master-only observations in the two-row reproduction, even though there were no matches. The correct retained count is zero. The ordinary branch has a temporary indicator for filtering; the empty branch needs equivalent logic.

**Proposed fix:** derive filtering directly from the known empty-side result code, or use a private temporary indicator regardless of `nogenerate`. Omitting a visible merge indicator must not change row selection.

**Required regression checks:** empty master/using/both, every `keep`/`assert` combination with/without `nogenerate`, and unchanged storage types and metadata. A suspected parallel type-list issue was also tested: two empty-using `str10` variables both retained `str10`; that suspicion is not reported as a defect.

### S13 — P2: the local build directory is not a complete installable package

**Location:** [build/ctools.pkg](/Users/Mike/Documents/GitHub/stata-ctools/build/ctools.pkg); [README.md:53](/Users/Mike/Documents/GitHub/stata-ctools/README.md:53); [validation/check_package.py](/Users/Mike/Documents/GitHub/stata-ctools/validation/check_package.py).

`python3 validation/check_package.py build` fails because these manifest entries are absent:

```text
ctools_mac_x86.plugin
ctools_windows.plugin
ctools_linux.plugin
```

The two local binaries present are the ARM platform plugin and generic `ctools.plugin`. Metadata-only validation intentionally skips platform-file existence and is not evidence of a complete installation directory.

**Scope qualification:** this is a confirmed defect in treating the examined local `build/` as a complete distribution. The live `main/build` installation endpoint and remote release artifacts were not inspected, so this report does **not** claim that every current remote installation is broken. The revised CI packaging step now assembles all four platform binaries and validates completeness, which is the right direction.

**Proposed fix:** publish/install only a fully staged directory that passes the complete checker, or produce platform-specific manifests intentionally. Clarify the source-checkout/manual-install distinction. Build ado/help and binaries from the same revision and record that revision with the package.

**Required regression checks:** clean `net install` from the staged directory; all manifest entries exist; every public ado's helpers are present; each supported architecture loads its matching plugin; no fallback silently masks a missing or stale platform artifact.

### S14 — P2: validation does not yet provide a comprehensive correctness gate

**Location:** [validation/validate_all.do:30](/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_all.do:30); [validation/run_stata_audit.py](/Users/Mike/Documents/GitHub/stata-ctools/validation/run_stata_audit.py); [.github/workflows/build.yml:195](/Users/Mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:195).

The framework now reports nonzero script returns and rejects missing/dimension-mismatched comparisons. Those earlier false-success problems are improved. Remaining coverage gaps are material:

- Nine components aborted with `r(679)` in the initial full run. Their remaining tests never executed. The direct web-data probe also failed with that code.
- The full run reported 70 skips: 22 import checks and 48 matching checks. Skips must remain visible in release evidence.
- The licensed CI job is optional; package publication explicitly accepts that job being skipped. When configured, it runs the two targeted audit suites, not the complete command suite.
- Passing bootstrap and PPML suites did not catch S08 and S05/S06 respectively.
- In the full-run snapshot, the suite compared a `strL` destring case against native success while the shared-I/O contract rejected `strL`. A bounded-text strL read implementation arrived afterward. Verify that change against the original case and add binary/oversized/boundary cases; the recorded failure is not evidence about the later implementation.

**Proposed fix:** make essential tests deterministic and offline, with vendored or generated fixtures and explicit reference dependency checks. Add invariants for sample membership, group relabeling, sampling-unit counts, output state, and rollback. Define a release policy requiring either a successful licensed correctness run or a separately recorded equivalent run; a skipped optional job is not statistical validation. Align intentional limitations with clearly named tests rather than weakening numerical thresholds.

**Required regression checks:** intentionally broken helper assertions fail CI; all component scripts reach a completion marker; unexpected skips/aborts fail release validation; the minimal cases in this report run on every release candidate.

### S15 — P2: installed help contains executable inaccuracies

**Locations and exact corrections:**

| Documentation | What it says | Actual behavior / proposed correction |
|---|---|---|
| [cpplmhdfe.sthlp:39](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe.sthlp:39), also options prose around lines 108–118 | `irls_tolerance()`, `irls_maxiter()`, `separation_tolerance()` | All three return `r(198)`. The ado accepts `irlstolerance()`, `irlsmaxiter()`, and `septolerance()`. Correct help or implement tested aliases. |
| [cbsample.sthlp:67](/Users/Mike/Documents/GitHub/stata-ctools/build/cbsample.sthlp:67) | `n(#)`, default `n(_N)` | The wrapper accepts a positional count, not an `n()` option. Show `cbsample #` and document cluster/stratum units after S08 is fixed. |
| [cbsample.sthlp:81](/Users/Mike/Documents/GitHub/stata-ctools/build/cbsample.sthlp:81) | `weight(newvar)` creates a new variable | The ado requires the variable to exist and replaces it. Document `weight(varname)` with a pre-created numeric destination, or change the interface deliberately. |
| [cbsample.sthlp:46](/Users/Mike/Documents/GitHub/stata-ctools/build/cbsample.sthlp:46) | Duplicates are implemented by keeping one copy and dropping unselected observations | The ado explicitly calls `expand`; describe actual replication/frequency-weight output. |
| [cbinscatter.sthlp:40](/Users/Mike/Documents/GitHub/stata-ctools/build/cbinscatter.sthlp:40), line 157 | `genxq()` generates bin assignments | The ado at lines 433–436 only prints “not yet implemented” and creates no variable. Its README admits this limitation. Implement it, or reject the option and make help consistent. |
| [cdecode.sthlp:93](/Users/Mike/Documents/GitHub/stata-ctools/build/cdecode.sthlp:93) | Automatic width follows the longest label | The wrapper caps width at 2045, and longer/escaped labels can become empty; see S07. State the actual boundary until fixed. |

**Proposed fix:** maintain one option/result inventory per command and execute help examples against deterministic fixtures. Explicitly test that advertised output variables are created. A message saying an accepted option is unimplemented is weaker than rejecting the request before doing expensive work.

The compatibility and platform pages that were missing at the start now exist, and the relative Markdown links checked at the final documentation pass resolved. Those earlier broken links are closed observations, not findings in this table.

### S16 — P3: developer instructions describe obsolete source registration and Stata launch procedures

**Location:** [DEVELOPERS.md:959](/Users/Mike/Documents/GitHub/stata-ctools/DEVELOPERS.md:959); [CLAUDE.md:5](/Users/Mike/Documents/GitHub/stata-ctools/CLAUDE.md:5) and its adding-command section; [Makefile:59](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:59).

The developer guide instructs contributors to add `NEWCMD_SRCS`/`NEWCMD_HEADERS` and append to `SRCS`. The current Makefile discovers C and header files automatically and uses `SOURCES`; that recipe is obsolete. The repository's older `runstata` instruction also conflicts with the user's current machine-wide `oldstata` requirement. The command inventories in these development documents omit newer modules such as PPML.

**Proposed fix:** update source-discovery, dispatcher/cleanup registration, manifest, and test-registration steps together. Replace stale local launch instructions with the applicable machine workflow, or clearly separate portable project instructions from machine-specific setup. Include an inventory check so adding a command requires an explicit decision about help, packaging, cleanup, and validation.

## Findings from September 20 that should not simply be repeated

The current tree contains substantial repairs. Evidence supporting their status includes:

- Strict floating-point flags and native compensation tests now pass.
- The targeted allocation-failure checks pass for the repaired sort/permutation paths. This does not cover all PPML allocation paths or S06.
- The audited shared I/O validates indices, rejects `strL`, and propagates tested read/write failures. A later bounded-text strL implementation arrived after the full source snapshot and needs its own verification. The closed missing-varlist defect under S01 and remaining S02 are caller-contract errors exposed by that improvement; do not remove validation to conceal them.
- The targeted tests pass for sampling qualifier preservation, typed sampling groups, existing label preservation, binscatter groups, import width rejection, IV residual/sample behavior, stable range variance, and quantile convergence/absorb rejection.
- The current `cqreg` cluster-relabeling tests pass. Its repaired mapping does not imply that IV/OLS/PPML use the same safe path.
- The older `cqreg_p.ado` has a questionable factor-variable `stdp` implementation, but the public estimator registers native `qreg_p`. Actual `predict, stdp` matched `_predict, stdp` in the audit. I therefore do **not** report the unused helper as a current public prediction bug; remove or clearly mark unused code during maintenance.
- PPML's 88-check suite passes. The September 20 broad covariance discrepancy is not established as still open. S03, S05, and S06 identify narrower remaining defects.
- The updated `noextend` fixture now separately checks native/ctools rejection of an undefined `D` label, then compares signed-label results on `A/B/C`. Its updated targeted suite passes all 14 checks.
- `ctools, environment_check` now makes a real plugin identity call and reports the audit build successfully. The earlier use of `program list` returned `r(111)` even for a loaded plugin; that observed false failure was fixed during this audit.
- Release metadata is synchronized; compiler-selection and dependency rejection tests pass; CI now includes broader path triggers and complete staged packaging.

These statements concern the tests actually executed, not certification of all command options or platforms.

## Proposed repair order and release criteria

1. **Restore primary data operations:** S01 and S02, with rollback checks. Address S07 before allowing destructive label replacement on unsupported inputs.
2. **Correct inference and memory safety:** S03–S06. Share typed group-remapping and retained-row mapping where practical, but verify each estimator's call boundary independently.
3. **Correct resampling and postestimation contracts:** S08–S12. Verify draw counts, actual matched-neighbor counts, predictions, and weight-expression equivalence.
4. **Make verification complete enough to trust:** S14, with offline fixtures and a licensed release gate. Run the new cases in addition to existing coefficient/VCE comparisons.
5. **Finish the distribution and documentation:** S13, S15, and S16; validate a clean staged installation and the documented supported platforms.

Before release, every P1 above should have a minimal regression test. Require sample-count consistency, cluster-relabeling invariance, correct resampling units, exact supported-label round trips, and unchanged user data after a rejected operation. Run sanitizer tests on the PPML separation/cluster path. A complete passing run must not hide an aborted component or unexplained skip.

Useful follow-on refactoring includes a single plugin loader with an enforced interface/version handshake, a common explicit index-space contract for checked SPI calls, typed grouping helpers shared across estimators, and command contexts with one ownership/cleanup path. Keep these changes tied to the regression cases; a broad rewrite is not necessary to repair the concrete defects.

## Evidence and limitations

Temporary evidence is under `/tmp/ctools-sep22-audit/`:

- `build.log`, `build-final.log`, `build/ctools_mac_arm.plugin`: fresh compilation and final frozen-source build.
- `audit_suites.do/.log`: initial targeted P1/P2 runs.
- `full_suite.do/.log`, `full_suite_final.do/.log`: whole-suite runs and machine-readable completion/return markers.
- `repro.do/.log`: bootstrap counts, cluster aliasing, ordinary merge failure, long labels, weight expressions, and the prediction control case.
- `repro2.do/.log`: estimator sample markers, IV cluster relabeling, winsorization, PPML prediction/help, and the web-data failure probe.
- `final_checks.do/.log`: corrected fixture/environment recheck, labeled extended missing values, literal-label round trips, and matching counts.
- `late_merge.do/.log`, `cardinality.do/.log`: late wrapper fix, updated P1 run, duplicate-key loss, and empty-input merge behavior.
- `native-p2-final.log`, `src-final/`, `file_hashes.txt`: final native checks and retained source identity.
- Initial native logs are `/tmp/sep22-native-p1.log` and `/tmp/sep22-native-p2.log`.

The reproductions record seeds and exact commands. A few exploratory cases were rejected by native Stata itself—most notably an attempt to attach a value label to ordinary `.`—and were corrected before drawing conclusions. Native unlabeled numeric values also decoded to empty strings in the control test; that behavior is **not** reported as a ctools defect. Filenames containing `threads(...)` worked in the tested import/export paths.

Temporary evidence may be removed by the operating system; the main inputs, results, causes, and required repairs are therefore included in this report. No Windows/Linux/Intel runtime execution, minimum-version Stata certification, minimum-OS execution, remote release inspection, full XLSX fuzzing, or large-memory stress campaign was performed. Native checks used UndefinedBehaviorSanitizer where configured; the full Stata plugin was not run under AddressSanitizer. The source-confirmed PPML out-of-bounds finding should be verified under that tool before accepting its repair.
