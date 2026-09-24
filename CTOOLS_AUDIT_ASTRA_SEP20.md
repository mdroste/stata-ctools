# ctools audit — Astra — September 20, 2026

**Recommendation: do not release this revision.** The audit found reproducible data corruption, incorrect statistical results, incomplete installation metadata, and validation paths that can report success without checking the intended result. Fix the P1 findings before publishing a release advertised as a replacement for the corresponding Stata commands.

Audited checkout: `dev`, commit `2518a9cdf7e4e9b91ae0f7ffaabee43286325da2`. The checkout was clean at the start. This audit changes only this report; it does not change implementation files or validation expectations.

## Scope and evidence

Reviewed the 17 command families, their ado interfaces and postestimation code, shared C infrastructure, package manifest, Makefile, GitHub workflow, validation framework, README, FEATURES, DEVELOPERS, and command documentation. Ran Clang static analysis over all **80 non-vendored C translation units**, built the ARM64 plugin from source, ran every validation script against both the existing plugin and a fresh build in an isolated copy, and constructed small synthetic reproductions of the highest-risk findings. Vendored miniz/libdeflate internals did not receive a separate dependency/security audit.

Evidence labels below distinguish **reproduced** behavior from **source-confirmed** defects and **portability risks** that need another platform. Static analyzer warnings were inspected; they are not treated as bugs merely because the analyzer emitted them.

The native build succeeded with **61 source warnings**, predominantly unused code, plus linker warnings that the installed static libomp was built for macOS 26 while the plugin targets macOS 11. Linux, Windows, Intel macOS, Stata 14/15, and clean-machine installation were not executed. The manifest findings concern the audited checkout; this report does not assume that the remote `main` branch has identical contents.

All Stata runs used `/bin/zsh -lic 'oldstata ...'`. Its required Legacy Time shortcuts account for July 1 timestamps inside Stata logs. All normally completing audit drivers ended with `exit, clear`; the wrapper also restored network time after the one crashed run.

### Validation results

The existing plugin and isolated fresh build produced the same component counts:

| Component | Recorded passes | Recorded failures | Script return code | Interpretation |
|---|---:|---:|---:|---|
| csort | 88 | 0 | 679 | Aborted on web-data access |
| cmerge | 16 | 0 | 679 | Aborted on web-data access |
| cimport | 313 | 0 | 0 | Completed; some web-data cases can be silently skipped |
| cexport | 16 | 0 | 679 | Aborted on web-data access |
| creghdfe | 17 | 0 | 679 | Aborted on web-data access |
| cqreg | 33 | 0 | 679 | Aborted on web-data access |
| civreghdfe | 16 | 0 | 679 | Aborted on web-data access |
| cdecode | 74 | 0 | 679 | Aborted on web-data access |
| cencode | 74 | 0 | 679 | Aborted on web-data access |
| cdestring | 167 | 0 | 0 | Completed |
| csample | 27 | 0 | 0 | Completed |
| cbsample | 29 | 0 | 0 | Completed |
| cbinscatter | 84 | 0 | 679 | Aborted on web-data access |
| cpsmatch | 67 | 0 | 0 | Completed; targeted checks expose missed defects |
| crangestat | 91 | 0 | 679 | Aborted on web-data access |
| cwinsor | 55 | 0 | 0 | Completed |
| cpplmhdfe | 47 | 41 | 0 | Completed with failed assertions |
| **Total** | **1,214** | **41** | — | **7 completed scripts; 10 aborted scripts** |

These are recorded assertion counts, **not a clean bill of health**. In particular, successful script exit is not sufficient: PPML returns zero despite its failed assertions, and the validation helpers have false-positive paths described below. No validation scripts or tolerances were changed.

## Prioritized findings

P1 means a release blocker for the advertised functionality: corruption, wrong estimates, broken installation, crashes on a valid failure path, or an unreliable release gate. P2 means a material defect or release-engineering/documentation gap to address next. P3 denotes refactoring after correctness is secured.

| ID | Priority | Finding | Evidence |
|---|---|---|---|
| A01 | P1 | Installation manifest omits commands and required helpers; three listed binaries are absent | Source/manifest check |
| A02 | P1 | cmerge modifies the wrong frame | Reproduced |
| A03 | P1 | cpsmatch writes results to the wrong variables | Reproduced |
| A04 | P1 | cqreg absorption does not estimate fixed-effect quantile regression | Reproduced |
| A05 | P1 | csample deletes observations outside if/in and overwrites a user variable | Reproduced |
| A06 | P1 | cencode replaces existing value-label mappings | Reproduced |
| A07 | P1 | cbinscatter drops numeric by-groups with ordinary nonconsecutive codes | Reproduced |
| A08 | P1 | cimport silently truncates long strings | Reproduced |
| A09 | P1 | IV residual options are ignored and residual predictions include fixed effects | Reproduced/source |
| A10 | P1 | crangestat variance suffers catastrophic cancellation | Reproduced |
| A11 | P1 | Build flags invalidate compensated arithmetic | Compiled C reproduction |
| A12 | P1 | PPML covariance matrices fail the project's own validation | Reproduced on both builds |
| A13 | P1 | Allocation-failure cleanup frees uninitialized pointers | Source/analyzer |
| A14 | P1 | Validation can report false success and miss missing results | Reproduced/source |
| A15 | P1 | Sampling treats string groups as numeric pointer bits | Reproduced |
| A16 | P1 | Matching without replacement can choose the farther control | Source/targeted trace |
| A17 | P1 | cqreg ignores maxiter and can misreport solver convergence | Reproduced/source |
| A18 | P2 | cencode evaluates literal data as Stata macro syntax | Reproduced |
| A19 | P2 | cencode noextend mishandles zero and negative label codes | Reproduced |
| A20 | P2 | cqreg truncates numeric cluster identifiers to integers | Reproduced |
| A21 | P2 | civreghdfe accepts center but does not apply it | Source-confirmed |
| A22 | P2 | UTF-32 input is silently misread as UTF-16 | Reproduced |
| A23 | P2 | Quoted positional cimport filenames fail | Reproduced |
| A24 | P2 | csort's final native sort ignores if/in | Reproduced |
| A25 | P2 | CI does not gate ado/package changes or run correctness checks | Workflow inspection |
| A26 | P2 | Distribution platform requirements exceed documented requirements | Source/portability risk |
| A27 | P2 | Cross-compiler discovery compares a path-plus-yes string to yes | Source-confirmed |
| A28 | P2 | Public documentation and release identity are inconsistent | Documentation/source |
| A29 | P2 | Shared I/O discards SPI errors and reports success | Source-confirmed |

### A01 — Complete and validate the installable package

**Locations:** `build/ctools.pkg:24–47`; `build/csort.ado:182`; `build/cmerge.ado:757,948`; `build/cexport.ado:304,754`; `build/civreghdfe.ado:1270`.

The manifest omits **11 ado files**: `_ctools_strw.ado`, `cbsample.ado`, `cdecode.ado`, `cdestring.ado`, `cencode.ado`, `civreghdfe_p.ado`, `cpplmhdfe.ado`, `cpsmatch.ado`, `crangestat.ado`, `csample.ado`, and `cwinsor.ado`. It also omits the nine corresponding command help files. `_ctools_strw` is unconditionally called by commands that *are* listed, so installing just the manifest's files breaks basic sorting, merging, and export. IV prediction names an omitted prediction program.

The manifest additionally lists `ctools_linux.plugin`, `ctools_mac_x86.plugin`, and `ctools_windows.plugin`, none of which is present or tracked in this checkout. Both present plugin files are ARM64 Mach-O bundles; the generic `ctools.plugin` is not a cross-platform fallback.

**Recommended fix:** generate or validate the manifest against a single command inventory, include all helpers and prediction programs, and require every manifest file to exist in the release staging directory. Test installation into an empty temporary adopath, with the development directory removed, then smoke-test every advertised command and prediction program. Do not use successful tests against the full `build/` directory as evidence that `net install` works.

### A02 — Restore the actual caller frame in cmerge

**Locations:** `build/cmerge.ado:368–387,407–408,434–435,572–573,771–772,862–864`.

The frame path assumes the master is in `default` and explicitly switches to `default` after reading the using data. It never records the caller's frame. Running a merge from another frame therefore resumes the merge against an unrelated dataset.

**Reproduction:** put `id={1,2}, val=999` in `default`, create a frame `auditmaster` containing only `id`, and run `cmerge 1:1 id using ...` there, where using data contain `val={101,102}`. The command returned zero, left `auditmaster` unmerged, switched the active frame to `default`, and replaced that frame's values with `{101,102}`.

**Recommended fix:** save `c(frame)` before any switch, use a unique temporary frame, and restore the original frame on every success/error path. Add a test that checks both the target frame and an unrelated sentinel frame after normal operation and after an error.

### A03 — Remove the extra index increment in cpsmatch output

**Locations:** `build/cpsmatch.ado:215–278`; `src/cpsmatch/cpsmatch_impl.c:1719–1735`.

The ado passes **one-based** variable positions. The C store loop adds one again to the weight, match, and support destinations. As a result, `_weight` is not populated, weights overwrite the next variable, and subsequent outputs are shifted as well. Ignored SPI return codes conceal this error.

**Reproduction:** with controls `(ps,y)=(.1,10),(.51,20)` and a treated observation `(.5,100)`, `cpsmatch treat, pscore(ps) outcome(y) noreplacement` returns zero, leaves every `_weight` missing, reports no control observations on support, and does not return an ATT. This occurs with both the existing and isolated fresh plugin.

**Recommended fix:** adopt one index convention at the ado/C boundary, remove the inappropriate `+1`, and check every store result. Assert exact `_weight`, support, match IDs, and ATT on a three-row fixture; do not rely solely on aggregate counters or variable existence. Retest all matching modes after the fix.

### A04 — Replace cqreg's least-squares absorption with a valid estimator

**Locations:** `src/cqreg/cqreg_regress.c:835–875,1020–1035`; `src/cqreg/cqreg_hdfe.c`; `build/cqreg.sthlp:63–68,88–89`; `FEATURES.md:111–122`.

The implementation projects y and X off the fixed effects using least-squares partialling, then solves a quantile regression on those residuals. The least-squares partialling identity does not generally extend to quantile loss. This estimates a different object from jointly minimizing quantile loss over slopes and group intercepts.

**Reproduction:** seed 54321; 200 observations; `group=ceil(_n/20)`, `x=rnormal()+group/3`, `y=3*group+2*x+(1+mod(group,3))*rexponential(1)`. At the median, native `qreg y x i.group` gives **2.0647133** for x; `cqreg y x, absorb(group)` gives **2.0258517**. Both return successfully.

**Recommended fix:** implement a quantile fixed-effects algorithm that optimizes the joint objective, or withdraw/restrict `absorb()` until it does. A documented residualized-outcome estimator would need a distinct interpretation and must not be represented as fixed-effect quantile regression. Validate against explicit group indicators on small one-way and two-way examples, including asymmetric and group-dependent errors.

### A05 — Preserve observations outside csample's selection and user variables

**Locations:** `build/csample.ado:132–135,187–194`.

The keep flag starts at zero for all observations, the plugin updates only the selected sample, and the ado drops every remaining zero. Thus an if/in qualifier deletes all observations outside the selection. The hard-coded scratch name `__csample_keep__` also causes an existing user variable with that name to be dropped and lost.

**Reproduction:** on ten observations, `sample 100 if id<=5` leaves all ten; `csample 100 if id<=5` leaves five. A pre-existing `__csample_keep__=123` also disappears.

**Recommended fix:** use `tempvar`, and either initialize excluded observations to keep or condition the final drop on the marked sample. Test if-only, in-only, combined qualifiers, zero/full percentages, and preservation of a user variable with the former scratch name. Existing sampling tests should compare the complement of the selection with native `sample` as well as checking the number sampled.

### A06 — Preserve existing cencode value-label mappings

**Locations:** `build/cencode.ado:225–236,263–269`; `src/cencode/cencode_impl.c:104–119,359–362`.

Except for the `noextend` path, existing labels are not supplied to C. Codes are regenerated from one, and the ado drops the prior value label before loading the new one. This changes the interpretation of every other variable using that label, not just the newly encoded variable.

**Reproduction:** define `existing` as `10 "A" 20 "B" 30 "unused"`, attach it to another numeric variable, and encode strings A/B with `label(existing)`. Native `encode` returns 10/20. `cencode` returns 1/2 and replaces the shared label definition with only those entries.

**Recommended fix:** load existing definitions for both extension and no-extension paths; reuse established codes and append new codes without deleting old entries. Test shared labels, unused existing entries, repeated encodes, and multi-variable encoding with a common label.

### A07 — Remap cbinscatter by-values to dense group IDs

**Locations:** `build/cbinscatter.ado:197–208,314–317`; `src/cbinscatter/cbinscatter_impl.c:321–334,645–684,478–496`.

Only string by-variables are encoded. Numeric values are cast to integers, group count is calculated as `max-min+1`, and observations are assigned using `value-1`; only a minimum of zero is specially handled. Ordinary numeric categories such as 10/20 therefore lose observations, create fictitious groups, and disagree with the ado's result-matrix dimensions. Negative/fractional and large numeric codes need the same review.

**Reproduction:** 40 observations split evenly between group 10 and group 20; `cbinscatter y x, by(group) nquantiles(4) nograph` returns zero and reports **11 groups**. Its 8-row bin matrix contains four bins for group 10 and four missing rows; group 20 is absent. Graph generation on the same data fails.

**Recommended fix:** map distinct values in the estimation sample to 1..G and retain the original values/labels separately. Allocate result matrices from G, not the numeric span; propagate matrix-store failures. Verify every group's observation count and graph legend.

### A08 — Stop silently truncating imported strings

**Locations:** `src/cimport/cimport_impl.c:264–266,296–298,383–384`; `src/ctools_config.h:441`; `build/cimport.ado:530–533,985–988`; `src/cimport/cimport_xlsx.c:2060–2061`.

String widths are capped at 2045, and the CSV extraction path loses an additional byte at the cap. Import succeeds without warning about the discarded content. Excel import has a corresponding fixed-width cap, though the exact Excel truncation boundary was not runtime-tested.

**Reproduction:** a CSV containing a `text` header and one 3,000-character string imports with return code zero into `str2045`; `strlen(text)` is **2044**.

**Recommended fix:** support long strings through an appropriate strL-capable path, or explicitly reject overlong fields before replacing the dataset. If truncation is offered, make it an explicit option with a count/warning. Test 2044/2045/2046-byte boundaries and multibyte UTF-8 boundaries for both CSV and Excel.

### A09 — Implement IV residual storage and correct residual prediction

**Locations:** `build/civreghdfe.ado:44–45,126–130,1270`; `build/civreghdfe_p.ado:50–74`.

`residuals(name)` is parsed but never used to create/store a residual variable or populate `e(resid)`. `residuals2` can delete an existing `_civreghdfe_resid` without replacing it. The prediction program computes `y-xb` for `residuals`, which retains absorbed effects. Its `xbd` branch requires an `e(resid)` that this estimator never supplies.

**Reproduction:** regress `y=10*group+2*x+noise` with IV and `absorb(group) residuals(savedres)`. Estimation returns zero; `savedres` does not exist; residual predictions retain approximately the 10*group component. `predict ..., xbd` reports that residual storage is required even after requesting it.

**Recommended fix:** store the structural residuals for the actual estimation sample, register the residual variable, and use those residuals for predictions involving absorbed effects. If a requested prediction cannot be computed, fail explicitly. Test exclusions, singleton removal, sample marks, and postestimation after saving/restoring estimates.

### A10 — Use stable range-variance calculations

**Locations:** `src/crangestat/crangestat_impl.c:263–344,786–791,1063–1068`.

Variance is computed from uncentered prefix sums as `(sum2-sum*mean)/(n-1)`, and negative results are clamped to zero. With a large location and small variation, cancellation destroys the variance; clamping then produces a plausible but incorrect zero.

**Reproduction:** ten doubles `value=1e12+_n`, `time=_n`; `crangestat (variance) rolling=value, interval(time . .)` returns **zero everywhere**. Native `summarize value` gives variance **9.1666667**.

**Recommended fix:** use centered/compensated aggregates or a stable mergeable-moments representation for range queries; use a stable fallback for ill-conditioned subtraction. Apply the same fix to SD and excludeself. Test invariance to adding a large constant and compare both small and large windows.

### A11 — Remove unsafe floating-point reassociation from accuracy-sensitive code

**Locations:** `Makefile:94–96,166–169,181–194,223–225`; `src/ctools_ols.h:58–79`; `src/ctools_ols.c:174–234,270`; PPML's calls to `dd_add_d`.

Every platform's main build enables `-ffast-math`, while shared numerical routines rely on the exact ordering of floating-point operations for double-double and compensated accumulation. Those assumptions are contradictory.

**Compiled reproduction using the repository's header:** `two_sum(1e16,1)` preserves a low component of 1 with `clang -O3`, but produces a low component of zero with `-O3 -ffast-math`. A non-inlined loop using `dd_add_d` on `[1e16,1,-1e16]` returns **1** in the strict build and **0** in the fast-math build.

**Recommended fix:** compile numerical kernels under strict floating-point semantics, and allow relaxed flags only in individually justified kernels. Account for LTO when isolating strict code. Add cancellation-sensitive arithmetic tests. The audit demonstrates loss of compensation; it does **not** establish that this alone explains the PPML failures or every README precision issue.

### A12 — Resolve PPML covariance discrepancies before release

**Locations:** `src/cpplmhdfe/cpplmhdfe_irls.c:1185–1568`; `validation/validate_cpplmhdfe.do`; `validation/benchmark_helpers.do:2500–2503,2620–2653`.

The unchanged PPML validation script records **41 failures out of 88 checks** with both binaries. All recorded failures concern `e(V)`, spanning offsets/exposures, two-way FE, robust and clustered inference, fweights/pweights, missing values, and larger synthetic samples. Reported agreement ranges roughly from 3.1 to 5.0 significant figures. The PPML comparison helper already uses a **five-significant-figure** default, looser than the shared seven-figure default.

**Recommended fix:** isolate the first synthetic failure and inspect terminal IRLS weights, the final weighted projection, residual construction, scaling, and finite-sample corrections separately. Compare the raw sandwich components with the reference estimator. Record the reference package versions. Do not lower the test threshold or classify these as harmless rounding without establishing the source and size of the discrepancy. The precise root cause remains unresolved by this audit.

### A13 — Make allocation-failure cleanup safe

**Locations:** `src/ctools_sort_radix_lsd.c:87–140,345–385`; `src/ctools_sort_radix_msd.c:114–165`; `src/ctools_types.c:484–508`.

Pointer arrays are allocated with `malloc`, and an earlier allocation can fail before their entries are initialized. The failure path then frees all entries. In `ctools_apply_permutation`, the per-thread allocation loop can break early, but cleanup still visits every thread's uninitialized slots. These are real invalid-free paths under memory pressure, not merely unchecked allocations.

**Recommended fix:** zero-initialize ownership arrays at allocation, or track precisely how many entries were initialized. Audit equivalent cleanup patterns in sample/merge sorting. Add allocator fault injection that fails each allocation in sequence and checks for a clean error rather than crashing Stata. These failure paths were source-confirmed with analyzer traces; system-wide memory exhaustion was not induced.

### A14 — Make validation failures impossible to hide

**Locations:** `validation/validate_all.do:49–59,94–100,124–139`; `validation/validate_setup.do:242–269,314–341`; `validation/validate_cimport.do:1179–1185`; `validation/benchmark_helpers.do:1080–1125`.

Three issues undermine confidence in the current suite:

1. The master runner captures each script's return code but never uses it when declaring PASS. Ten scripts aborted with `r(679)` during this audit with zero recorded assertion failures. Those rows would be printed as PASS by the stock runner.
2. `assert_var_equal` turns a missing computed comparison into 15 significant figures: in Stata, numeric missing compares greater than 15. A one-row test comparing `.` with `42` recorded **one pass and zero failures**.
3. Matrix comparisons do not first assert equal dimensions. Some web-data tests silently skip on failure; missing outputs and skipped cases are not consistently represented in totals. Matching tests explicitly skip the ATT-SE comparison despite compatibility claims.

**Recommended fix:** treat any nonzero script return, missing result, or unexpected missing output as a failure; reset counters before each script; compare dimensions and sample membership explicitly; distinguish PASS/SKIP/FAIL. Make standalone scripts return nonzero on failed assertions. Cache/version necessary datasets or provide explicit offline tests. Add a small self-test for the validation helpers themselves. Keep tolerances unchanged.

### A15 — Compare string groups as strings in csample and cbsample

**Locations:** `src/csample/csample_impl.c:100–114,495`; `src/cbsample/cbsample_impl.c:97–111,552,567`; both ado files accept unrestricted grouping varlists.

Both group-boundary helpers always read `data.dbl`, including when the loaded variable is a string. The union then exposes string-pointer bits as doubles. Group membership becomes unrelated to string equality. Numeric extended missing values are also all collapsed together by the same helper rather than compared by actual code.

**Reproduction:** ten observations, five with string group A and five with B. `csample, count(1) by(group)` keeps **all ten**, rather than one per group. `cbsample 1, strata(group)` likewise yields **ten**, rather than two. Both return zero on the fresh plugin.

**Recommended fix:** use a shared typed key comparator or encode grouping variables to dense IDs before boundary detection. Test string, numeric, mixed multi-key, and distinct extended-missing groups in both commands. Share the corrected implementation rather than maintaining the duplicate helpers.

### A16 — Correct the two-sided nearest-neighbor search

**Locations:** `src/cpsmatch/cpsmatch_impl.c:1245–1265,1340–1415`.

On finding a left candidate, the search sets `best_idx=left` and breaks without decrementing left. The subsequent test checks `best_idx == left+1`, so it takes the wrong branch and does not compare the right candidate. A treated score .50 with available controls .10 and .51 can therefore select .10. The shifted match output in A03's reproduction identifies control observation 1, consistent with this code path. Separately, the treated sorting buffer leaves `obs_idx` uninitialized even though the tie-break comparator reads it.

**Recommended fix:** compare the nearest available candidate on each side explicitly, apply a documented tie rule, and initialize every comparator key. After fixing A03, assert exact matches for nearer-left, nearer-right, ties, exhausted controls, descending order, and calipers. Validate or reject option combinations that the one-neighbor no-replacement branch cannot honor.

### A17 — Honor cqreg iteration limits and actual solver status

**Locations:** `src/cqreg/cqreg_fn.c:32,959,1213–1228,1392`; `src/cqreg/cqreg_regress.c:1000,1038–1041,1066–1068`; `build/cqreg.ado:294,526`.

The default Frisch–Newton solver loops to a hard-coded 500, ignoring `config.maxiter`. It records convergence internally but returns an iteration count. The caller declares convergence whenever that count is positive. A numerical breakdown or exhausted loop after at least one iteration can consequently be posted as successful convergence.

**Reproduction:** a 400-observation noisy regression with `maxiter(1)` returns normal estimates and `e(convcode)=0`. Source inspection confirms that the requested limit is never consulted by this loop.

**Recommended fix:** return an explicit status distinct from iteration count, honor the supplied iteration limit, reject nonfinite termination, and propagate nonconvergence to `e(convcode)` and an appropriate command error policy. Test deliberate exhaustion, singular designs, and numerical failure without relying on a positive iteration count.

### A18 — Preserve literal dollar signs and macro delimiters in cencode labels

**Locations:** `src/ctools_hash.c:761–824`; `build/cencode.ado:267–269`.

Labels are emitted into a do-file and executed. The escaping check handles control characters but does not prevent Stata macro expansion inside the generated compound-quoted strings.

**Reproduction:** set global `AUDIT_WORD` to `changed`, construct a string containing literal `$AUDIT_WORD` with `char(36)+"AUDIT_WORD"`, and cencode it. The resulting label is **changed**, not the source string.

**Recommended fix:** serialize label content so it cannot be interpreted as Stata source syntax. Test dollar signs, local-macro delimiters, nested compound quotes, control characters, Unicode, and decode round trips. This should be treated as data preservation, not just display formatting.

### A19 — Separate hash lookup status from signed label values

**Locations:** `src/ctools_hash.c:702–724`; `src/cencode/cencode_impl.c:200–207,246–252`.

Existing-label parsing treats a negative return from insertion as allocation failure, although a legitimate stored label code can be negative. Lookup uses zero as absence, and the encoding path accepts only `code>0`.

**Reproduction:** `noextend` with label `-1 "A"` produces `r(920)`; with label `0 "A"`, A silently becomes missing. Native encode preserves those label codes.

**Recommended fix:** return a found/error status separately from the signed code. Test negative, zero, positive, and absent labels without conflating valid values with sentinels.

### A20 — Preserve cluster identity in cqreg

**Locations:** `src/cqreg/cqreg_regress.c:192–228`; `src/cqreg/cqreg_vce.c:197–255,812–878`; `build/cqreg.ado:114–129`.

Numeric cluster values are cast directly to `ST_int` before remapping. Distinct fractional categories become the same cluster, and out-of-range double identifiers cannot be represented safely. The ado only creates dense IDs for string variables.

**Reproduction:** 400 observations in 20 categories with `cluster=ceil(_n/20)/100`; clustered cqreg fails with `r(504)` after every category truncates to zero.

**Recommended fix:** remap original numeric values to dense integers using equality-preserving logic, as is already done elsewhere in the project. Verify covariance invariance under one-to-one relabeling of the same groups, including fractional and large double codes.

### A21 — Implement or reject civreghdfe center

**Locations:** `build/civreghdfe.ado:801`; `build/civreghdfe.sthlp:55`; `src/civreghdfe/civreghdfe_impl.c:118,1280`; `src/civreghdfe/civreghdfe_estimate.c:1294–1297`.

The documented HAC centering option is parsed and passed through, then explicitly discarded with `(void)center` and a TODO. Users receive an apparently accepted option without its promised calculation.

**Recommended fix:** center the intended score vectors for every supported affected estimator/VCE combination, or reject the option with a clear unsupported-option error. Add a fixture where centered and uncentered score covariance differ. The nearby Kiefer/FE compatibility TODO at `src/civreghdfe/civreghdfe_estimate.c:1947–1951` also merits a targeted reference comparison; this audit did not establish a separate Kiefer error.

### A22 — Reject or correctly convert UTF-32 input

**Locations:** `src/cimport/cimport_encoding.c:121–130`.

The detector recognizes the UTF-32LE BOM prefix but deliberately treats it as UTF-16LE. Embedded NULs then corrupt field parsing.

**Reproduction:** UTF-32 encoding of `name,value\nAlice,123\n` imports successfully as two observations with two empty string columns named `v` and `v1`.

**Recommended fix:** distinguish UTF-32 BOMs before UTF-16 detection; implement conversion or return a clear unsupported-encoding error. Never assign high-confidence UTF-16 detection to a recognized UTF-32 stream.

### A23 — Strip syntactic quoting from positional import filenames

**Locations:** `build/cimport.ado:24–44`.

The positional `anything` is copied into `using` with its quotes still embedded. The later compound quoting makes those quotes part of the filename. The source comment explicitly promises both positional and using syntax.

**Reproduction:** `cimport delimited "/tmp/.../long.csv", clear` returns `r(601)` for an existing file; `cimport delimited using "/tmp/.../long.csv", clear` works.

**Recommended fix:** parse one filename token with normal Stata filename quoting semantics and reject extra tokens. Cover spaces, quoted paths, both forms, and the Excel dispatcher if it shares this pattern.

### A24 — Define and enforce csort's if/in behavior

**Locations:** `build/csort.ado:13,215–229`; `docs/README_csort.md:15,38–42`.

After the plugin processes a qualified subset, the wrapper runs unqualified `sort varlist, stable` over the entire dataset. Thus the normal command does not preserve the requested scope and incurs another full sort. It can also hide sorting defects in tests that only inspect the final key order.

**Reproduction:** ids 1..4 with keys 4,3,2,1; `csort key if id<=2` reverses all four observations, including those excluded by the condition.

**Recommended fix:** either reject qualifiers consistently with native sort, or define subset ordering and implement it without a final global re-sort. Verify that excluded rows remain unchanged under the chosen contract. Test the C permutation independently of the wrapper's native sort, including non-key columns and `nosortedby`.

### A25 — Make CI a correctness and packaging gate

**Locations:** `.github/workflows/build.yml:3–14,34–52,243–258,273–303`.

The workflow only filters source/Makefile changes, so ado, help, manifest, and validation-only changes do not trigger the ordinary build paths. Push branches are `main` and `develop`, while this checkout's development branch is `dev`. Verification checks binary existence/dependencies, not command behavior. The combined archive copies ado/help/plugins but omits `ctools.pkg` and `stata.toc`, so it is not itself a complete Stata net-install directory. The plugin commit job checks out the current tip of main rather than explicitly binding publication to the source revision that produced its binaries; overlapping pushes can mix revisions.

**Recommended fix:** trigger on every release-relevant file and the actual development branch; gate publication on manifest integrity, deterministic C tests, and available licensed Stata smoke/regression checks. Package the install metadata. Pin artifacts to the build commit and reject publication if main has advanced, or publish immutable versioned releases. Add a concurrency policy. Do not call `file`/`ldd` success a validation pass.

### A26 — State and test the actual platform contract

**Locations:** `Makefile:104–119,136–156,169,184,194,223–269`; `README.md:45–47`; `.github/workflows/build.yml:87–99,131–143,170–204`.

All x86 distribution paths target Haswell, which enables AVX2/FMA/BMI rather than a baseline x86-64 instruction set. Linux dynamically requires libgomp and can dynamically require OpenBLAS. macOS falls back to a Homebrew-path libomp dylib when a static archive is unavailable, while CI merely prints a warning. Windows linkage leaves pthread/runtime deployment to be checked. These conditions are inconsistent with the unconditional “does not require any dependencies” claim.

The local build also links a libomp archive built for macOS 26 into a bundle declared compatible with macOS 11. Successful compilation does not establish that the output runs on macOS 11.

**Recommended fix:** choose a documented minimum CPU/OS baseline, or use runtime dispatch for newer instructions. Build runtimes against that deployment target. Verify dynamic dependencies on clean machines and bundle/install them where appropriate. Treat unsupported dependencies as packaging failures. Windows and old-CPU failures remain risks, not runtime-confirmed findings from this Mac.

### A27 — Fix cross-compiler availability probes

**Locations:** `Makefile:176,188,215–216`.

`$(shell which compiler ... && echo yes || echo no)` emits the compiler path **and** `yes` when installed. Comparing that result with exactly `yes` is false. This prevents selection of the MinGW fallback and the Linux cross-compiler even when present. The MinGW-only Docker route is affected by the same selection logic.

**Recommended fix:** use `command -v ... >/dev/null 2>&1` for a Boolean probe, or retain the discovered executable path and test whether it is nonempty. Exercise compiler-present/compiler-absent branches with a controlled PATH. No installation of additional compilers is needed to test this logic.

### A28 — Reconcile documentation with the executable interface

**Locations and recommended corrections:**

| Location | Stale or misleading statement | Correction |
|---|---|---|
| `README.md:12`; `build/ctools.pkg:18`; `build/ctools.ado:1,105` | Versions disagree: 0.9.1, 1.0.1, and header 1.0.2 | Generate release identity from one source and expose the plugin's identity too |
| `README.md:38` versus `build/cpplmhdfe.ado:3` | Advertised `cppmlhdfe` does not match implemented `cpplmhdfe` | Choose the public spelling, update references, and provide an alias if needed |
| `build/ctools.pkg:20`; ado `version 14.1` statements | Manifest requires 14.0 while commands require 14.1 | Align the minimum version and test that minimum |
| `README.md:40` | Link uses `FEATURES.MD`, but the file is `FEATURES.md` | Correct case for GitHub/case-sensitive filesystems |
| `README.md:88` | Calls the sort memory option `streaming` | Document actual `stream(#)` syntax |
| `docs/README_csort.md:9–23` | Says three algorithms and LSD default | Describe the implemented algorithm set and auto selection |
| `docs/README_cimport.md:40,66,100–107` | Advertises a `fast` option and temporary-DTA path absent from current syntax | Remove the obsolete option/path or implement and test it |
| `build/cimport.sthlp:30,110–112`; `docs/README_cimport.md:21` | Says comma is the default | Explain automatic delimiter detection and explicit overrides |
| `FEATURES.md:93–106`; `validation/validate_csample.do:178,529` | Says native sample does not support by() | Correct the comparison; add native by-group parity tests |
| `README.md:18`; matching/quantile help | Broad syntax/functionality replacement claims | Publish a command-by-command compatibility table, including weights, residual prediction, estimands, and ATT-SE differences |
| `CLAUDE.md` adding-command/build instructions | Refers to per-module Makefile source lists despite automatic discovery | Update repository instructions; the user's current oldstata/Git workflow takes precedence |

Also update the top-level command list and package description for PPML and the omitted newer commands. Add the actual project license file and review third-party notice inclusion when assembling the distribution; README's MIT statement alone is not a complete release artifact inventory. Performance claims should identify benchmark date, hardware, Stata/reference versions, options, and dataset shape.

### A29 — Propagate shared Stata I/O errors

**Locations:** `src/ctools_data_io.c:157–163,214–216,270–273,496–507,526–532`; `src/cimport/cimport_impl.c:1580–1660`; `src/cbinscatter/cbinscatter_impl.c:478–496`.

Core load/store paths ignore the return values of `SF_vdata`, `SF_sdata`, `SF_vstore`, and `SF_sstore`. `store_variable_thread` unconditionally marks success. Import counts attempted rows as stored, and binscatter ignores matrix-store failures. Thus invalid indices, unsupported string representations, and write failures can become apparently successful results; A03 and A07 show why this matters in reachable code.

**Recommended fix:** make shared store functions return status, aggregate the first worker failure safely, and stop before publishing success/results. Validate index and destination dimensions at the boundary even when `SD_FASTMODE` is enabled. Establish an explicit strL policy rather than treating all string reads as fixed-width buffers. Test mocked SPI read/write failures and ensure the ado restores or clearly reports partial state.

## Refactoring opportunities

### R01 — Centralize plugin loading and enforce an ABI/version handshake (P3)

Platform detection and plugin loading are copied across almost every ado, with differing cache mechanisms and fallbacks. Introduce one private loader that returns the resolved path, validates architecture, and checks a plugin-reported API/build version. Use unique names or a defined restart requirement when binaries change; avoid silently falling back to a stale generic plugin. Extend `ctools, environment_check` to check a real call, versions, and required helpers, and return nonzero on failure.

One audit run explicitly loaded the fresh binary while the repository build remained on the adopath; a later command loaded another copy and Stata crashed with SIGSEGV in OpenMP. The crash report shows both plugin images, with frames crossing their OpenMP runtimes. The isolated single-build runs completed. This is evidence for a **mixed-binary/runtime compatibility risk**, not evidence that every clean install crashes. `KMP_DUPLICATE_LIB_OK` in `src/ctools_plugin.c:53–70` suppresses a runtime check; it does not prove that mixed runtimes interoperate.

### R02 — Unify typed grouping, allocation ownership, and cleanup (P3)

Start with the duplicated sampling group helpers exposed by A15. Reuse equality-preserving remapping across cbinscatter, quantile clustering, and HDFE. Give each command a context with explicit ownership and one cleanup path, using zero-initialized pointer arrays. This should also make allocation fault injection and error rollback tractable. Avoid introducing another generic layer before fixing the existing index/type contracts.

### R03 — Consolidate numerical kernels and remove dead implementations (P3)

The source build produces 61 warnings, mostly unused implementations in quantile BLAS/linalg/sparsity, IV VCE/tests, import helpers, and binscatter. There are also nearly empty compatibility translation units such as `creghdfe_ols.c` and `creghdfe_types.c`. Inventory actual call sites, keep one tested implementation per mathematical primitive, and remove obsolete variants/debug scaffolding after parity tests exist. Separate numerical kernels from Stata I/O so they can be tested under sanitizers without a licensed Stata process. Keep strict floating-point policy explicit at the kernel boundary.

### R04 — Generate release metadata and make builds reproducible (P3)

Use one command/release inventory to produce the manifest, command listing, version strings, and documentation checks. Separate release staging from development output; untrack obsolete object/debug products; pin the runtime/toolchain used for distribution. Prefer object-level compilation with correct dependencies over recompiling every module for each phony platform target. Add license/notice and package completeness checks to staging. These changes should follow A01/A25/A26 rather than substitute for their fixes.

## Recommended repair sequence and release criteria

1. **Repair the gate first:** A14 plus CI triggers and deterministic/offline fixtures. A test suite that accepts missing results cannot establish that subsequent fixes work.
2. **Protect user data and prevent crashes:** A02/A03/A05/A06/A07/A08/A13/A15/A18/A19/A29, including allocation failures, error rollback, and non-target state preservation.
3. **Correct the estimators and numerical behavior:** A04/A09/A10/A11/A12/A16/A17/A20/A21. Require reference comparisons of coefficients, sample membership, covariance matrices, and postestimation outputs; validate the intended objective as well as numerical closeness.
4. **Publish a complete, coherent package:** A01/A22–A28. Test a clean installation for each supported platform and the stated Stata/CPU/OS baseline. Publish binaries and ado files from the same revision.
5. **Then refactor:** R01–R04, preserving the new behavioral tests.

Before release, every P1 should have a minimal regression test; all component scripts must finish or explicitly record justified skips; no failed assertion or nonzero script return may be reported as PASS. A full successful run remains outstanding because of the defects above and the web-data limitations. No performance measurements from this audit support updated speedup claims.

## Reproduction material and audit limits

Temporary evidence was kept outside the repository at `/tmp/ctools-audit-astra-sep20/`:

- `build.log` and `analyze.log`; detailed analyzer traces under `analyzer/`.
- `suite.log` and `source-suite.log` for the existing and isolated fresh binaries.
- `repro.log`, `repro2.log`, `repro3.log`, `final-checks.log`, `numeric-checks.log`, and `group-checks.log` for targeted examples; associated `.do` drivers are alongside them.
- `dd.c` and `dd_loop.c` plus strict/fast-math executables for the compensation checks.
- `mirror/` contains unchanged copies of validation/ado files with only the freshly built plugin available in its local build directory.

The early targeted drivers include exploratory cases that used a positional import filename or attempted to read an uncreated IV residual; later drivers isolate those failures explicitly. The conclusions above use the corrected successful reproductions or the identified error itself. Temporary files may be removed by the operating system, so the report records the essential inputs, outputs, locations, and recommended regression checks directly.

The mixed-runtime crash report is `~/Library/Logs/DiagnosticReports/stata-mp-2026-09-20-161449.ips`. Its fault is `EXC_BAD_ACCESS/SIGSEGV` at address `0x8` in OpenMP; the report lists Stata's libomp and two different ctools plugin images. No sanitizer execution, external-platform execution, large-memory stress campaign, full XLSX fuzzing, or remote-release verification was completed. Findings not established by this audit are identified as risks or follow-up checks rather than confirmed defects.
