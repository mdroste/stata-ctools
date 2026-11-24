# Release audit: cqreg, cpsmatch, chdfe, and crangejoin

**Date:** November 24, 2026. **Scope:** the four public commands, their Ado/Mata helpers, native implementations, and the shared code on their execution paths. **Original audit disposition:** audit and proposed fixes; production code was not changed during the audit. The subsequent requested implementation is recorded in [Implementation follow-up](#implementation-follow-up).

**Timing:** approximately 22 minutes, within the requested 40-minute limit. Work began at 10:07 a.m. Pacific; final evidence/link verification completed at 10:28 a.m.

**Release recommendation:** fix F01, F04, and F08 before release. They can silently produce seriously incorrect inference, destroy an input variable, or join the wrong file. Resolve the remaining five findings through a fix or an explicit, enforced limitation. F05 and F06 reproduce behavior in the installed `psmatch2`; fixing them requires a deliberate compatibility decision, not merely another parity test.

All selected existing Stata checks passed against both the supplied Windows plugin and a fresh source build. The findings below come from additional reproductions, mathematical invariants, and fault injection. Passing the existing suite therefore does not establish release readiness.

## Itemized findings

| ID | Priority | Component | Finding |
| --- | --- | --- | --- |
| F01 | P1 | cqreg | Regressor location/scale can corrupt coefficients and covariance while reporting convergence |
| F02 | P2 | cqreg | Multiplying all importance weights by a small constant changes the fitted solution |
| F03 | P2 | cqreg | Required scalar-publication failures are discarded |
| F04 | P1 | cpsmatch | Generated output names can collide with and destroy inputs |
| F05 | P2 | cpsmatch | `ties neighbor(k)` can duplicate donors and omit boundary ties |
| F06 | P2 | cpsmatch | Translating a Mahalanobis covariate can change all matches |
| F07 | P2 | chdfe | Global normalization can erase entire small-weight or small-slope groups |
| F08 | P1 | crangejoin | Literal macro characters in the using filename can select a different file |

P1 means a high-priority release blocker in this assessment. P2 means a confirmed issue to fix or explicitly constrain; it is not a statement that the affected output is safe to use.

### F01 — cqreg accepts inaccurate fits and covariance for uncentered or extreme-scale regressors

**Locations:** [collinearity screening and design preparation](G:/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:1120), [solver initialization](G:/GitHub/stata-ctools/src/cqreg/cqreg_fn.c:1133), [convergence and crossover](G:/GitHub/stata-ctools/src/cqreg/cqreg_fn.c:1499), and [covariance inversion](G:/GitHub/stata-ctools/src/cqreg/cqreg_vce.c:109).

**Reproduction:**

```stata
clear
set seed 71234
set obs 120
generate double x0 = rnormal()
generate double y = 3 + 2*x0 + rnormal()
generate double x = x0
cqreg y x
replace x = x0 + 1e8
cqreg y x
replace x = x0 + 1e12
cqreg y x
```

An intercept makes these translations the same statistical model. In the audit, the centered fit had slope **1.79763999416**, slope SE **0.119821743599**, and check-loss objective **48.994090134624**. At an offset of `1e8`, the slope remained approximately correct but its SE fell to **0.0013761695315**, about **87 times too small**. At `1e12`, the slope was **0.0972086271946** and the objective was **101.66312775089**, over twice the available centered-fit objective. Both calls returned success and `e(convcode)==0`.

The broader scale probe also found success with a zero rescaled slope at `x=x0*1e160`, and erroneous omission at `x=x0*1e-200`. These are additional manifestations of unsafe regressor scaling; they are much more extreme than the location example.

**Cause:** response normalization and centered collinearity screening do not normalize the design used throughout the main solve and VCE. Raw cross-products lose the small variation around a large regressor level, and can overflow/underflow at extreme units. The VCE uses normal-equation inversion with a regularization fallback. A finite result and the existing solver-gap test do not certify accurate coefficients or covariance; crossover compares against the current solver objective, which itself may be wrong.

**Proposed fix:** center nonconstant design columns and scale their spreads safely before every main/auxiliary solve and covariance calculation. Apply the exact coefficient/intercept transformation and full covariance transformation when publishing. Use a rank-revealing factorization or explicit conditioning checks for the covariance, and validate feasibility/optimality independently before declaring convergence. Return a clear numerical error if a safe result cannot be represented; do not silently substitute a regularized covariance for valid inference.

**Acceptance checks:** translated and rescaled versions must preserve fitted values, objective, slopes after conversion, and slope SEs to tolerances justified by input rounding. Cover IID, robust, cluster, and the supported density estimators, plus collinear columns and factor-variable expansion. The raw shifted `qreg` call is **not** an adequate oracle: in this example the installed `qreg` omits the shifted regressor at large offsets. Compare against the well-conditioned equivalent model and its transformed covariance.

**Evidence:** [numeric reproduction and reference comparison](G:/GitHub/stata-ctools/temp/audit-four-20261124/numeric_candidate.log), [broader scale probe](G:/GitHub/stata-ctools/temp/audit-four-20261124/probes_candidate.log).

### F02 — cqreg's convergence test is not invariant to common importance-weight scaling

**Locations:** [weighted rows](G:/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:1245), [absolute gap tolerance](G:/GitHub/stata-ctools/src/cqreg/cqreg_fn.c:1243), and [crossover objective tolerance](G:/GitHub/stata-ctools/src/cqreg/cqreg_fn.c:1550).

**Reproduction:** use the centered data from F01, generate a constant double weight, and compare `cqreg y x [iw=w]` with `w=1`, `1e-10`, and `1e-20`.

| Common weight | Slope | Objective divided by common weight | Convergence code |
| --- | ---: | ---: | ---: |
| 1 | 1.797639994160 | 48.994090134624 | 0 |
| 1e-10 | 1.781257080845 | 48.994161226400 | 0 |
| 1e-20 | 1.938867061447 | 49.538431871400 | 0 |

Multiplying every positive importance weight by the same constant does not change the minimizer. The installed `qreg` accepts `1e-10` and returns the original slope and normalized objective. It rejects some smaller-weight cases, so the mathematical invariance, rather than reference output, is the oracle there.

**Cause:** the implementation multiplies both response and design rows by the original weights, then uses an absolute duality-gap tolerance. The crossover allowance also contains an absolute `1e-10` term. At small common scales, these tolerances can accept a materially inferior solution relative to the actual objective.

**Proposed fix:** normalize weights by a safe common factor for numerical solving and evaluate scale-aware primal/dual feasibility and objective-gap criteria. Preserve original weight semantics separately: frequency-weight counts, bandwidth calculations, reported sums, and VCE must not be accidentally changed by solver normalization. Ensure crossover cannot undo an otherwise valid accuracy guarantee.

**Acceptance checks:** constant rescalings and heterogeneous positive importance weights must preserve coefficients and normalized objectives. Include pweights as routed by the wrapper, and separately test frequency-weight accounting. Retain explicit errors for cases whose output cannot be represented.

**Evidence:** [successful numeric probe](G:/GitHub/stata-ctools/temp/audit-four-20261124/numeric_candidate.log).

### F03 — cqreg ignores failures when publishing required scalars

**Locations:** [store_scalar](G:/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:98), [store_results](G:/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:768), [collinearity publication](G:/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:1171), and [basis/nonuniqueness publication](G:/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:1331).

**Reproduction:** a native harness includes the actual production implementation, supplies a two-coefficient state, makes the host's scalar-save callback return error 459, and allows matrix stores to succeed. It prints:

```text
AUDIT_CQREG_PUBLICATION rc=0 rejected_scalars=11 attempted=11
```

Thus `store_results()` reports success after all eleven required scalar writes in this configuration fail. Its matrix writes do propagate errors. The wrapper subsequently consumes these scalars to assemble estimates.

**Impact:** a failure in the plugin-to-Stata result handoff is hidden, leaving an incomplete result protocol and potentially stale scratch values. This audit proves the discarded callback errors; it does **not** claim to have reproduced an incorrect posted `e(b)` under ordinary host memory pressure. The wrapper may instead fail later with a less useful error.

**Proposed fix:** return and check the host status for every required scalar write, including collinearity and basis metadata. Abort publication with the original error, clean scratch results, and preserve the command's intended estimation-state behavior. Treat optional timing values separately from required statistical results.

**Acceptance checks:** inject failure at each required scalar and matrix save and require the exact nonzero status to propagate, with no successful partial estimate. Include a successful call followed by a failing call to expose stale-state reuse.

**Evidence:** [native reproduction source](G:/GitHub/stata-ctools/temp/audit-four-20261124/cqreg_publication.c). It was compiled with GCC/LTO and run against the included production `store_results()`.

### F04 — cpsmatch can destroy variables that are inputs to the current call

**Locations:** [incomplete input-name guard](G:/GitHub/stata-ctools/build/cpsmatch.ado:48), [output deletion and captured markout](G:/GitHub/stata-ctools/build/cpsmatch.ado:158), and generated `_s_`/`_self_` outcomes later in the wrapper.

**Reproduction:**

```stata
clear
input double(t ps y _y)
0 .1 10 100
0 .8 20 200
1 .2 30 300
1 .7 40 400
end
cpsmatch t, pscore(ps) outcome(y _y)
```

The call returns **success**. The original `_y` outcome, containing `100,200,300,400`, is replaced by the generated match for `y`, containing missing values for the two controls and `10,20` for the treated observations. The requested second outcome is no longer evaluated from its original data.

A second reproduction renames `ps` to `_n1` and calls `cpsmatch t, pscore(_n1) outcome(y)`. The command returns 111 **after deleting the score input**. These are in-memory data changes; no on-disk input file was overwritten by the audit.

**Cause:** the guard covers only a fixed list of reserved names. It misses `_n#`, outcome-derived names, and underlying variable names hidden inside factor-variable expressions. The wrapper drops generated output names before finishing input use. `capture markout` suppresses an error when an outcome has already disappeared.

**Proposed fix:** resolve the complete set of input variables, including factor-variable base names, and calculate every output name for the selected options before any mutation. Reject every input/output collision with an explanatory error and unchanged data, or compute from protected temporary copies if supporting these names is intentional. Validate name lengths and collisions across outcomes, smoothing outputs, and AI outputs together. Do not suppress missing-input errors during sample construction.

**Acceptance checks:** cover treatment, score, covariate, Mahalanobis, and outcome collisions with fixed outputs, `_n1`, `_y`, `_s_y`, and `_self_y`; verify values, variable existence, and sample remain unchanged on rejection. Include `outcome(y _y)` specifically.

**Reference qualification:** `psmatch2` also deletes `_y` in this example, but returns 111; `cpsmatch` continues and reports success. Reference compatibility does not justify destructive handling of an input.

**Evidence:** [initial reproduction](G:/GitHub/stata-ctools/temp/audit-four-20261124/probes_candidate.log), [reference and score-input confirmation](G:/GitHub/stata-ctools/temp/audit-four-20261124/confirm_candidate.log).

### F05 — cpsmatch duplicates donors with ties and multiple neighbors

**Locations:** [psm_scan tie expansion and neighbor loop](G:/GitHub/stata-ctools/src/cpsmatch/cpsmatch_impl.c:232), and the documented [boundary-tie contract](G:/GitHub/stata-ctools/build/cpsmatch.sthlp:137).

**Reproduction:**

```stata
clear
input double(t ps y)
0 .1 10
0 .1 100
0 .1 1000
1 .2 2000
end
cpsmatch t, pscore(ps) outcome(y) ties neighbor(2)
```

There are three equally close, distinct controls. The documented all-boundary-ties behavior should produce three matches with equal weights, a matched outcome of **370**, and ATT **1630**. Instead the command reports `_nn=4`, weights **0.25, 0.5, 0.25**, matched outcome **302.5**, and ATT **1697.5**. With `neighbor(3)` it reports `_nn=5` despite only three controls and yields matched outcome **442**.

A separate case with control scores `0.1,0.2,0.2`, control outcomes `10,100,1000`, and treated score `0.1` shows the other side of the defect: `ties neighbor(2)` includes only the first two controls, omitting the third control tied at the second-neighbor boundary.

**Cause:** the scan first appends ties at the nearest distance, then independently adds `k-1` neighbors starting beside the original nearest donor. With replacement, that second loop can append an already included donor. It also does not extend the final kth-distance boundary.

**Proposed fix:** determine the kth eligible distance, build a list of distinct donor observations through that boundary, and include each tied donor once. Derive counts, weights, matched outcomes, and variance inputs from that same list. Keep tie-breaking and caliper semantics explicit. If strict historical compatibility is required, reject the affected option combination pending a documented policy rather than silently returning duplicate matches.

**Acceptance checks:** assert unique donor IDs and `_nn` no larger than the eligible donor population; use hand-calculated means and weights for ties at the first and kth distances. Exercise ATE and self-matching/AI paths as applicable.

**Reference qualification:** installed `psmatch2` reproduces the `neighbor(2)` error and omitted boundary tie, while it errors at `neighbor(3)` in the first example. These tests must include an analytic oracle, not only reference equality.

**Evidence:** [complete tie probes](G:/GitHub/stata-ctools/temp/audit-four-20261124/ties_candidate.log).

### F06 — cpsmatch Mahalanobis matching is not translation invariant

**Locations:** [rank screening and covariance](G:/GitHub/stata-ctools/build/cpsmatch.ado:931), [nearest-neighbor distance construction](G:/GitHub/stata-ctools/build/cpsmatch.ado:1064), and the corresponding smooth-matching distance calculation in `_cpsmatch_smooth`.

**Reproduction:**

```stata
clear
input double(t x y)
0 0 10
0 2 20
0 5 50
1 .1 100
1 1.9 200
1 4.9 500
end
cpsmatch t, mahal(x) outcome(y)
replace x = x + 1e6
cpsmatch t, mahal(x) outcome(y)
```

Before translation, the matched outcomes are `10,20,50`, giving ATT **240**. After adding `1e6`, every treated observation matches outcome **50**, giving ATT **216.6666666667**. Translation cannot change Mahalanobis distances when the covariance metric is unchanged.

**Cause confirmed for this reproduction:** `_rmdcoll` operates on the uncentered covariate; the audit's direct check returns `o.x` at offsets of `1e5` and above for this dataset, so subsequent `fvrevar` loses the informative dimension. There is also a second numerical weakness in the distance formula `x'Wx - 2x'Wz + z'Wz`. A direct scalar diagnostic at offset `1e8` gives **-2** from that expanded formula where the difference-first squared distance is **0.01**. Fixing rank screening alone would therefore leave cancellation in distance evaluation.

**Proposed fix:** construct centered/scaled temporary matching covariates before rank screening and covariance estimation, retaining a consistent transformation for a supplied `w()` matrix. Compute quadratic distances from the coordinate differences, preferably in a stable transformed metric. Do not treat large negative computed squared distances as valid close matches.

**Acceptance checks:** donor IDs, distances, support, weights, and effects must remain invariant to common translations and appropriate unit changes. Cover nearest-neighbor, kernel/LLR Mahalanobis paths, supplied metrics, and nearly collinear covariates.

**Reference qualification:** installed `psmatch2` shows the same changed matches. Decide whether to correct this behavior or reject ill-conditioned inputs explicitly; matching the reference is insufficient evidence of mathematical correctness.

**Evidence:** [matching comparisons](G:/GitHub/stata-ctools/temp/audit-four-20261124/matching_join_candidate.log), [rank and distance diagnostics](G:/GitHub/stata-ctools/temp/audit-four-20261124/confirm_candidate.log).

### F07 — chdfe's global normalization can erase complete groups

**Locations:** [chdfe_unit_scale](G:/GitHub/stata-ctools/src/chdfe/chdfe_impl.c:29), [weight normalization](G:/GitHub/stata-ctools/src/chdfe/chdfe_impl.c:147), and [slope normalization](G:/GitHub/stata-ctools/src/chdfe/chdfe_impl.c:223), feeding the shared HDFE factors and solver.

**Weight reproduction:**

```stata
clear
set obs 8
generate g = ceil(_n/2)
generate double y = _n
generate double w = cond(_n<=4, 1e-200, 1e200)
chdfe y [aw=w], absorb(g) gen(r_)
```

Weights are constant within each two-observation group, so every group's residuals should be **-0.5, +0.5**. The first two groups instead retain **1,2,3,4**, while the large-weight groups are residualized correctly. The command returns success.

**Slope reproduction:**

```stata
clear
set obs 12
generate g = ceil(_n/3)
generate double x = cond(g<=2, 1e-200, 1e200)*mod(_n,3)
generate double y = mod(_n,3)
chdfe y, absorb(g#c.x) gen(r_) keepsingletons
```

Within each group, `y` is exactly proportional to `x`, so residuals should be zero. The small-slope groups retain **1,2,0**; only the large-slope groups are absorbed.

**Cause:** dividing an entire weight or slope vector by a scale based on its global maximum underflows the small groups to zero. Positive weights are validated before this division, not after it. The solver then operates on effectively absent groups and can satisfy its numerical stopping criterion for the altered problem.

**Proposed fix:** add range-preservation checks around normalization and norm construction. For slopes, group-local scaling preserves each group's projection and can avoid this loss. Weighted absorption needs scaling/accumulation that preserves relative weights in every projection; arbitrary independent group reweighting is not valid for all crossed-FE problems. If the implementation cannot preserve a supplied range, reject it explicitly before output mutation instead of claiming success. Review response-vector normalization for the same failure mechanism.

**Acceptance checks:** analytic group residuals, weighted orthogonality, and idempotence across heterogeneous group scales; test crossed FEs and slope-plus-intercept cases. Include explicit errors when support would be lost by underflow. These are extreme-range cases, which is why this is P2 rather than a general claim that routine absorption fails.

**Evidence:** [weight and slope probes](G:/GitHub/stata-ctools/temp/audit-four-20261124/probes_candidate.log).

### F08 — crangejoin can read a different file after re-expanding a literal filename

**Locations:** [wrapper forwarding of local 0](G:/GitHub/stata-ctools/build/crangejoin.ado:30), [frame/snapshot helper calls](G:/GitHub/stata-ctools/build/crangejoin.ado:164), and [using-file load](G:/GitHub/stata-ctools/build/crangejoin.ado:201).

**Reproduction conditions:** create two files in the same directory: literal `using_$AUDITRJ.dta` containing key 1 and payload **123**, and `using_expanded.dta` containing key 1 and payload **999**. Set global `AUDITRJ` to `expanded`. Construct the first path literally with Mata's `char(36)` and pass it protected by `macval()`:

```stata
global AUDITRJ expanded
mata: st_local("literal", "PATH/using_" + char(36) + "AUDITRJ.dta")
// A direct read proves the literal file is accessible:
use `"`macval(literal)'"', clear
// Reset master to a one-observation dataset containing key=1, then:
crangejoin key 0 0 using `"`macval(literal)'"'
```

The direct `use` reads payload **123**. `crangejoin` returns success with payload **999** in both default and `snapshot` modes. The using argument is re-expanded inside the wrapper/helper chain and the command silently selects the other file. The reproduction does not depend on nonexistent paths or an unprotected caller-side expansion.

**Proposed fix:** preserve literal text across every forwarding stage: the original argument list, parsed using path, helper argument construction, and final `use`. Use `macval()` and suitable compound quoting consistently, or a path transport that does not rescan its contents. Fixing only the last `use` is insufficient if an earlier helper call has already changed the path.

**Acceptance checks:** use a literal filename and an expansion-target decoy file with different sentinel values. Verify correct source selection for default, snapshot, and fast paths; include defined/undefined globals, backticks/apostrophes where the filesystem permits them, spaces, and `keepusing()`. Failure cases must preserve the master according to the selected rollback mode.

**Reference qualification:** installed `rangejoin` also reads the wrong sentinel. This is a file-selection correctness problem, not a native interval-search disagreement, and should be fixed despite reference parity.

**Evidence:** [complete creation/read/join script](G:/GitHub/stata-ctools/temp/audit-four-20261124/matching_join.do), [direct-use control and both join modes](G:/GitHub/stata-ctools/temp/audit-four-20261124/matching_join_candidate.log).

## Audit coverage and verification

The review followed the commands beyond their entry points:

| Area | Reviewed execution paths |
| --- | --- |
| cqreg | Ado syntax, sample/weight/FV preparation and result assembly; native dispatch/state ownership; FN solver, initialization, QR fallback, crossover/basis logic; linear algebra and BLAS fallbacks; weighted/unweighted sparsity and auxiliary fits; IID/robust/cluster/kernel VCE; result publication and delegated `qreg_p` prediction |
| cpsmatch | Ado sample/support construction and outputs; score estimation; native nearest-neighbor scan/index, replacement/ties/ATE/AI paths, kernels/radius, option parsing, statistics and stores; Mata Mahalanobis nearest-neighbor and smoothing/LLR helpers |
| chdfe | Ado option/absorb expansion; compact sample and factor loading; singleton removal, remapping/counting; shared projections/acceleration/convergence, slopes, normalization, nesting/connected components and degrees of freedom; output publication |
| crangejoin | Default frame, snapshot and fast wrappers; using schema and string/strL handling; native using index, interval lookup, output sizing/writes and rollback; stable ordering/grouping and mixed-type key handling |
| Shared code | Relevant filtered-data selection/observation maps, load/store contracts, allocation and bounds checks, grouping/remapping, plugin loading/dispatch, thread setup, and Ado rollback helpers |

This was a time-bounded source audit with targeted execution, not exhaustive branch coverage or a proof that every shared helper is defect-free.

### Existing Stata checks

Run through `validation/run.py --stata cqreg cpsmatch chdfe crangejoin`, first against `build/ctools_windows.plugin`, then against the freshly compiled candidate:

| Suite | Result on each plugin |
| --- | --- |
| cqreg | 300 checks passed |
| cpsmatch | 365 checks passed |
| crangejoin | 102 checks passed |
| cqreg_sample_scale | 371 checks passed |
| cpsmatch_statistics | 511 checks passed |
| chdfe | 131 checks passed |
| cqreg_performance | Passed |
| Selected option/data matrix | 210/210 scenarios passed; log reports 480 assertions passed |

The six counted command/statistics suites total **1,780 checks per plugin**, excluding performance and the separate option matrix. These counts should not be interpreted as independent coverage of all branches.

### Native checks and additional probes

- Existing `cpsmatch_paths`, `cpsmatch_paths_oom`, `hdfe_paths`, `newcommands_paths`, and `cqreg_native` harnesses passed in serial and OpenMP builds with undefined-behavior sanitizer trap mode. The join checks reside in the broader newcommands harness.
- The matching allocation-failure harness exercised 16/18/19/27/25 failure points across its five serial paths and 16/18/19/51/25 in parallel, alongside publication-failure checks. No new failure was observed in those existing checks.
- A deliberate signed-overflow control exited with the expected trap status, confirming the sanitizer mode was active. These runs were **not** AddressSanitizer or ThreadSanitizer runs.
- Dedicated licensed Stata probes reproduced F01/F02/F04–F08 on the supplied and rebuilt implementations; subsequent detailed numeric/reference controls were run against the rebuilt candidate. F03 was demonstrated directly against production source through a mocked host callback.
- The draft combined `targeted.do` and an initial version of `confirm.do` had audit-script errors after some probes. They are not counted as passing runs. Their corrected, completed counterparts are `numeric_candidate.log` and `confirm_candidate.log`.

### Provenance and retained evidence

The workspace HEAD was **`0388e3075907dc080e29849ffcf639d971b9ced8`**, package version **1.2.6**, plugin protocol **3**. Unrelated working-tree changes were present, and additional unrelated edits appeared during the audit; none were modified by this work. The fresh plugin was compiled from the working tree at build time, rather than described as a clean HEAD build. Final verification found no tracked changes in the four audited implementations or the reviewed shared HDFE, data-plan, ordering, and rollback files.

| Artifact | Identity |
| --- | --- |
| Supplied Windows plugin | Embedded revision `9eadc285273add6408b97d27cf9517bfa714b147`; SHA-256 `14b6e1b25fb72598edba2d6f80fb731df05560d8d720abc3f96d3ba661883fa7` |
| Fresh Windows plugin | Embedded revision `audit-four`; SHA-256 `6baad87db2a48fba8bbf983b8ad300c842bd3496e576a7a2fad7750fa75333b8` |
| Execution environment | Windows, StataNow/MP 19.5; GCC 16.2; fresh plugin built with O3, OpenMP, LTO, `-fno-fast-math`, and `-ffp-contract=off` |

Evidence is retained under [temp/audit-four-20261124](G:/GitHub/stata-ctools/temp/audit-four-20261124). It includes the source-build script and manifest, original/candidate suite logs, standalone Stata probes, native harness sources/binaries, and the scalar-publication reproduction. This directory is ignored scratch storage; retain it separately if the report will be reviewed after workspace cleanup. The report itself is the durable audit artifact.

### Limits and release follow-up

Runtime verification was on Windows/Stata 19.5. Linux/macOS builds, older Stata releases (including the pre-frame fallback), address/thread sanitizers, very large observation counts, and every allocation/cancellation point were not exhaustively exercised. Existing native publication/rollback tests passed, but complete rollback after every possible Ado failure is not established by these results.

For the release gate, add the analytic and failure-injection cases specified under each finding, rerun the existing suites against the **exact final release binaries**, and test the other supported platforms. Keep reference parity checks, but supplement them with independent properties: translation/scale invariance, unique donor counts, within-group orthogonality, exact file provenance, and propagation of failed host writes.

Suggested implementation order: prevent input loss and wrong-file selection (F04/F08); normalize and validate cqreg's numerical problem and VCE (F01/F02); repair scalar error propagation (F03); decide and implement the matching semantics (F05/F06); then harden HDFE range handling (F07).

## Implementation follow-up

All eight findings have been addressed in the working tree, with regression coverage and validation against a fresh Windows plugin. Existing unrelated changes were retained. The original findings and audit measurements above describe the pre-fix implementation.

| Finding | Implemented change | Regression coverage |
| --- | --- | --- |
| F01 | Normalize continuous design columns before rank screening, the main and auxiliary fits, and covariance calculations. Transform the entire coefficient vector and covariance back to original units, including intercept and cross-covariances. Preserve the indicator parameterization used for zero-density factor levels. Remove the regularized covariance fallback and reject unrepresentable output with error 498. | Translations through `1e12`, scales `1e-100` and `1e100`, six IID/robust/cluster/density configurations, a full three-parameter covariance transformation, and explicit rejection of unrepresentable units. Existing factor-variable and collinearity checks also pass. |
| F02 | Normalize weights by a common power of two for numerical solving, retaining original weights for the command's statistical semantics. Restore reported objective units and reject unrepresentable weight sums/deviation statistics. | Common iweight and pweight scales from `1e-20` to `1e10`, heterogeneous importance weights, coefficient/objective invariance, and the distinct IID-iweight versus robust/pweight SE behavior. Existing frequency-weight checks pass. |
| F03 | Propagate required scalar and matrix publication errors immediately, including collinearity, basis, and nonuniqueness metadata. Optional timing writes remain best-effort. | A native harness calls the production result-publication function, injects failure at each of its 13 scalar and four covariance writes, verifies the original status and immediate stop, and verifies recovery on a subsequent successful call. |
| F04 | Resolve factor-variable base inputs and enumerate generated names before deleting outputs. Reject aliases involving scores, outcomes, matching indices, smoothed outcomes, and AI outputs. Do not suppress missing-outcome errors during sample marking. | Input preservation for `outcome(y _y)`, `pscore(_n1)`, `c._weight`, `_s_y`, and `_self_y` collisions. |
| F05 | With replacement, walk distinct donors once and include every donor through the kth-distance boundary. Counts, weights, and matched means use that same list. Retain the existing sequential availability policy for no-replacement `neighbor(1)`. | Analytic tied examples for neighbors 1–3 and kth-boundary inclusion, plus an independent exhaustive-distance oracle in the option suite. Existing no-replacement parity and matching-statistics suites pass. |
| F06 | Center and normalize Mahalanobis covariates before rank screening/covariance calculation. Compute distances from coordinate differences and transform supplied metric matrices to the normalized coordinates. | Translations by `1e6` and `1e8` for default and supplied metrics, nearest-neighbor matching, kernel matching, and local-linear matching. |
| F07 | Enforce an explicit numerical-range limit: return error 498 if normalization or weighted squared norms would erase nonzero inputs or entire slope groups. Abort before output publication. | Extreme weight, slope, and response ranges; combined weight/slope norm underflow; input preservation and absence of generated output, including `clear` mode. This is a safe rejection policy, not an implementation of arbitrarily wide-range group-local scaling. |
| F08 | Preserve literal arguments and filenames through the wrappers. Avoid Stata's filename rescan in the varlist-using form for paths containing macro characters; ordinary paths retain selective column loading. | Literal global- and local-macro characters with conflicting decoy files, in default, snapshot, and fast modes, both with and without `keepusing()`. |

The matching changes deliberately diverge from the confirmed `psmatch2` tie and translated-distance defects. This is documented in the command help and compatibility guide. The `cqreg` normalization fixes the audited failures; it does not add a separate general-purpose optimality certificate to the solver.

Two existing `cqreg` tests needed more accurate acceptance criteria. The small `auto` subset now uses the existing alternate-optimum helper, which checks the independently evaluated objective and statistical results instead of requiring identical coordinates for a nonunique optimum. The tiny-response host-posting test now accepts either a correctly transformed nonzero covariance or the explicit posting error; it still rejects a silently discarded covariance.

### Post-fix validation

The final command runs completed with **1,852 counted checks passed, zero failed, zero skipped**, plus the `cqreg_performance` suite and **210/210 option/data cases**. These are the final results per suite; the join suite was rerun after preserving selective loading for ordinary paths and adding the local-macro filename case.

| Suite | Final result |
| --- | --- |
| `cqreg` | 300 passed |
| `cqreg_sample_scale` | 414 passed |
| `cqreg_performance` | Completed successfully |
| `cpsmatch` | 378 passed |
| `cpsmatch_statistics` | 511 passed |
| `chdfe` | 135 passed |
| `crangejoin` | 114 passed |
| Four-command option/data cases | 210/210 passed; the 48 join cases also passed on the final join rerun |
| Native numerical/path checks | cqreg, matching, matching allocation/publication failures, HDFE, and join-containing newcommands harnesses passed in serial and OpenMP builds with UBSan trap mode |
| New cqreg publication harness | All injected scalar/covariance write failures and recovery passed |

Durable regression sources are [the four-command regression cases](G:/GitHub/stata-ctools/validation/support/audit_four_regressions.do), [the publication harness](G:/GitHub/stata-ctools/validation/native/cqreg_publication.c), and [its automatically discovered runner](G:/GitHub/stata-ctools/validation/native_cqreg_publication.py). The public command suites invoke the new Stata cases. The tied-match option oracle is in [coverage_helpers.do](G:/GitHub/stata-ctools/validation/support/coverage_helpers.do).

Execution used Windows, StataNow/MP 19.5, and GCC 16.2. The tested plugin has embedded revision `audit-four-fixes`, SHA-256 **`acc862941a3d7ff2440c50a43371b02c21d2c2152704c1320b090463299a851c`**, and was built with O3, OpenMP, LTO, no fast-math, and floating-point contraction disabled. It is a working-tree build that includes the unrelated source edits present at build time. The supplied release binaries were not overwritten. Linux/macOS builds and older Stata versions still require the normal release validation.

Logs and the isolated candidate are under [temp/audit-four-fixes-20261124](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124): [final command validation](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/validation_final.log), [final join validation](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/join_final.log), [native paths](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/native.log), [native cqreg](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/native_cqreg_final.log), [native HDFE](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/native_hdfe_final.log), and [publication failures](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/publication_final/test_cqreg_publication_native.log). This remains ignored scratch evidence; preserve it separately if needed after cleanup.

Broader repository checks were also attempted and were not all green. Package staging/inventory/Stata metadata/description/licenses, release integrity, documentation, seeds, and benchmark-evidence checks passed. Compiler-probe and runner checks encountered Windows symlink/junction restrictions; make/bash-based checks could not find those tools; the incremental-build check lacked its default compiler on PATH; and the global option-inventory check reported unrelated `cexport` entries (`rename`, `origdbfdate`, `vallabfile`, `version`, `novallabels`). These results are retained in [repository_checks.log](G:/GitHub/stata-ctools/temp/audit-four-fixes-20261124/repository_checks.log). They were not counted as passing or addressed by changing unrelated work.
