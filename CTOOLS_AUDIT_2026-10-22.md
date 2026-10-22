# ctools code audit — October 22, 2026

**29 actionable findings: 4 P1, 22 P2, and 3 P3.** Twenty-eight were demonstrated with targeted Stata/native probes or captured generated graph syntax; one Windows-specific issue is source-confirmed and needs execution on Windows. The most consequential findings concern incorrect statistical inference, silent loss of imported data, and unintended data modification.

Reviewed starting revision: `e470fdab096ee3afc14598a7ff52e8ef532c03f4` (version 1.0.2). A fresh Mac ARM plugin was built with identity `e470fda-audit`. Production code and shipped binaries were not modified by this audit. Three parallel reviewers were used, below the requested maximum of 20.

**Priorities:** P1 = high priority; fix before the next release because inference or data integrity can be materially wrong. P2 = important correctness, compatibility, or failure-handling issue. P3 = lower-severity functionality or returned-metadata issue. Priorities consider impact and trigger; they are not claims of exploitability.

The October 21 report was treated as historical context. Its repaired issues were not simply repeated. The XLSX metadata and Windows encode findings identify paths omitted from those earlier fixes.

## Findings at a glance

| ID | Priority | Issue |
|---|---|---|
| A01 | P1 | IV overall F includes the intercept |
| A02 | P1 | Fitted clustered quantile VCE uses the wrong density model |
| A03 | P1 | Corrupt XLSX metadata becomes silent missing data |
| A04 | P1 | SHP record-boundary truncation silently loses observations |
| A05 | P2 | Bootstrap option parsing changes quoted conditions and overwrites data |
| A06 | P2 | Unicode ignore() strips unintended characters |
| A07 | P2 | Destring conversion discards labels and characteristics |
| A08 | P2 | Valid XML text is dropped or truncated |
| A09 | P2 | XLSX control characters and escaped strings do not round-trip |
| A10 | P2 | Range means/sums become wrong at the optimized-path threshold |
| A11 | P2 | IV small option is ignored |
| A12 | P2 | IV no-absorb DOF options do not affect VCE |
| A13 | P2 | No-constant IV uses centered total sum of squares |
| A14 | P2 | Valid small-sample IV models are rejected |
| A15 | P2 | PPML residual prediction rejects accepted time-series outcomes |
| A16 | P2 | Quantile e(rank) counts omitted/base coefficients |
| A17 | P2 | Merge keepusing() mishandles wildcard/range varlists |
| A18 | P2 | Empty-master observation-number merge ignores keepusing() |
| A19 | P2 | Empty-input merge filters away rows before failed assert() |
| A20 | P2 | Binscatter deletes a preexisting user frame |
| A21 | P2 | Binscatter silently drops singleton groups |
| A22 | P2 | Matching memory failures silently alter the estimation sample |
| A23 | P2 | Matching input/output name collisions destroy input on error |
| A24 | P2 | Treatment validation and matching/reporting definitions disagree |
| A25 | P2 | Invalid matching settings silently produce different/no matches |
| A26 | P2 | Encode label files bypass Windows Unicode path handling |
| A27 | P3 | Binscatter connect option does not connect points |
| A28 | P3 | Binscatter reportreg is a no-op |
| A29 | P3 | Binscatter stored metadata names deleted temporary variables |

## Scope, method, and validation

The review covered all public command families and their ado/C boundaries: import/export, merge and sorting, sampling, string conversions, interpolation/splitting, range joins/statistics, winsorization, matching, binscatter, quantile regression, OLS/IV/PPML with fixed effects, shared transport/memory/threading infrastructure, packaging, and CI/release controls. Review depth varied by area; this is a broad source audit with targeted execution, not a line-by-line proof or exhaustive certification of vendored codecs.

Local execution used macOS 27.0 on Apple Silicon, Apple Clang 21.0.0, Python 3.14.5, and StataNow/MP 19.5. The fresh plugin used the supported pthread-only build because the static OpenMP runtime was unavailable.

| Validation | Result |
|---|---|
| Fresh Mac ARM build and dependency/baseline checker | Passed |
| Release metadata, canonical package synchronization, command inventory, package metadata | Passed at initial snapshot |
| Native suite registry | 37 attempted; 26 passed; 11 failed during missing OpenMP setup/build |
| Existing licensed component suites | 2,219 assertions passed; 0 failed; 0 skipped across 12 completed suites |
| Additional rangejoin component | Could not run: reference `rangejoin` was absent |
| New targeted probes | Confirmed the findings below, including native allocation-failure and parser probes with ASan/UBSan |
| Offline official data fixtures | 24 checksum-pinned datasets verified |

The 12 completed components were csort (351), cmerge (230), csample (28), cbsample (29), cencode (239), cdecode (219), cdestring (167), cipolate (31), csplit (132), cimport (360), cexport (213), and cqreg (220). These passing suites do not cover the newly reproduced cases.

The 11 native drivers that could not complete were test_bigdata_native, test_order_parallel, test_p1_native, test_p2_native, test_rowstage_native, test_sort_keys_native, test_transport_adaptive_native, test_transport_native, test_transport_scheduling, test_transport_store_native, and winsor_transport_native. Several ran their serial cases successfully before reaching the missing OpenMP dependency. **The overall native gate failed; it is not reported as green.** These particular failures are setup limitations, not observed numerical failures.

Native Stata built-ins supplied independent references where possible. Third-party rangestat, psmatch2, binscatter and rangejoin were unavailable in this session; installed ivreg2 could not run because ranktest was missing. IV comparisons therefore use built-in ivregress plus explicit Wald/DF identities. Binscatter singleton behavior is assessed against its documented discrete-bin contract; the connection finding checks generated graph syntax rather than a rendered screenshot. Each relevant item states these limits.

No Windows/Linux/Intel execution, minimum supported Stata 14.1 run, minimum-OS execution, or full licensed cross-platform release gate was performed. The existing publication gate was inspected and no new bypass was established. Shared transport, row staging, allocation, ordering and lifecycle paths were reviewed without an additional high-confidence defect beyond A26; that is not a proof of absence.

**Concurrent work:** The checkout was clean when the audit began. Changes to cipolate/csplit/cwinsor source and wrappers, a native test harness, and performance-audit scripts appeared during the audit from other work. They were preserved and were not included in the fresh plugin snapshot or certified by these results. The reported finding locations concern files unchanged by that concurrent work.

## A01. [P1] civreghdfe omnibus F includes the intercept

Locations: [package/ado/civreghdfe.ado:1021](/Users/mike/Documents/GitHub/stata-ctools/package/ado/civreghdfe.ado:1021). Without absorb(), the wrapper appends an explicit _cons to the coefficient vector. Its Wald calculation then tests every coefficient jointly, divides by the complete K, and posts df_m=K. The conventional model F tests slope coefficients, not the intercept. On the 60-row deterministic probe with y=10+.05*x+.4*u+normal error, native `ivregress 2sls y (x=z), small` gives F(1,58)=1.1510264; `civreghdfe y (x=z), small` gives F(2,58)=2312.5382. `test x` immediately after ctools still yields F=1.1510264. This can turn an insignificant model into apparently overwhelming significance solely because its outcome has a nonzero mean.

Fix: construct the joint hypothesis from active nonconstant coefficients, excluding omitted terms; use the hypothesis covariance rank for numerator DF and post/display the same F/df_m. Add explicit-constant no-absorb IV/OLS tests with large intercept and near-zero slopes, plus rank-deficient clustered covariance tests.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.txt).

## A02. [P1] Default fitted cluster VCE incorrectly imposes a common density

Locations: [src/cqreg/cqreg_regress.c:1078](/Users/mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_regress.c:1078); [src/cqreg/cqreg_vce.c:705](/Users/mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_vce.c:705). The per-observation fitted density is computed only when vce_type==ROBUST; CLUSTER goes through the IID scalar sparsity path. Its cluster covariance consequently uses s² (X'X)^-1 M (X'X)^-1 rather than the fitted density Hessian (X'DX)^-1 on both sides of the clustered score meat. The kernel clustered path already uses the latter structure correctly. A common scalar conditional density is not valid for general heteroskedastic clustered data.

Evidence: 2,000 seeded synthetic rows with x~U(0,5), y=2+.5*x+(1+4*x)*normal error and 100 clusters. Native `qreg y x, vce(excluster cl, fitted)` and `cqreg y x, vce(cluster cl)` give identical slope .0661785. Native slope variance .03246172; cqreg .04463079 (ratio1.3748746). With kernel(gaussian) both produce .04926808 exactly, controlling for coefficient solver and cluster mapping. Inspection of installed Stata qreg.ado GetVCE confirms fitted robust/cluster both use per-observation Hessian.

Fix: estimate per-observation fitted densities for both robust and cluster, construct the fitted Hessian once, and combine it with clustered score sums using the intended finite-sample/zero-residual policy. Add heteroskedastic clustered parity tests, including homoskedastic control and singleton-per-observation clusters. This is a systematic covariance error, not a numerical-tolerance difference.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/cqreg.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/cqreg.txt).

## A03. [P1] XLSX integrity protection omits shared strings and other metadata parts

**Locations:** [src/cimport/cimport_xlsx.c:583](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:583) (shared strings), [src/formats/xlsx_xml.c:489](/Users/mike/Documents/GitHub/stata-ctools/src/formats/xlsx_xml.c:489) (generic parser finalization), [src/io/cio_xlsx.c:28](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_xlsx.c:28) (invalid shared-string reference becomes empty). Related metadata paths: [src/cimport/cimport_xlsx.c:259](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:259), `380-408`, `716-727`.

**Trigger:** A structurally sound ZIP and worksheet reference entries 0 and 1, but `xl/sharedStrings.xml` ends at `<sst><si><t>first</t></si><si><t>second` instead of completing the second entry and closing the root. ZIP sizes and CRCs are valid. Truncation can arise while a workbook is produced/repacked.

**Observed impact:** `cimport excel ..., clear` returns 0 and imports `(first, 1)` followed by `("", 2)`. Native `import excel` rejects this workbook with r(603) and retains the prior data. The native sanitizer probe also returns `rc=0 strings=1`. A clean worksheet therefore bypasses the recently added whole-document integrity check, and missing string-table entries are silently converted to empty cells. Truncated styles can analogously lose date classification; relationships parsing additionally ignores parser failure and extraction failure. Those related manifestations are source-confirmed, while the shared-string loss is reproduced end-to-end.

**Proposed fix:** Require structural validation for every consumed XML part, including shared strings, styles, workbook and relationships, before applying callbacks; propagate parse failures instead of treating them as success. Distinguish absent optional parts from present-but-corrupt parts. Validate that each referenced shared-string index exists and reject unresolved references. Preserve the old dataset on validation failure. Add valid-ZIP fixtures with each metadata part truncated at every structural boundary, missing/invalid shared-string indices, and broken relationship XML.

**Relation to previous audit:** October 21 finding 2 was fixed for worksheet XML. This is the remaining metadata path, not a claim that its worksheet fix is absent.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.txt).

## A04. [P1] SHP imports accept files truncated exactly between records

**Locations:** [src/io/cio_shp.c:70](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:70) reads the file header without checking its file-length field; [src/io/cio_shp.c:130](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:130) treats any EOF at a record boundary as successful termination; success is published at `239-248`.

**Trigger:** A valid two-point SHP declares 156 bytes in its main header. Delete the final 28-byte record, leaving 128 bytes and the original header. The remaining record is internally complete.

**Observed impact:** Fresh-plugin `cimport shp ...` returns 0 with one point instead of two. The complete fixture returns 0 with two points. The parser never checks the 32-bit big-endian file-length field at header offset 24, measured in 16-bit words, so losing whole records looks like a successful shorter dataset. Native Stata's behavior for this fixture was not tested; the result is established by direct ctools execution and its file-length contract.

**Proposed fix:** Read and bounds-check the advertised file length using 64-bit arithmetic, require at least the 100-byte header, verify actual stream length matches it, and ensure every record ends within the advertised boundary. Require the final position to match before reporting success. Add truncations at all record boundaries, incomplete final records, overclaimed and underclaimed lengths, and large-file cases.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.txt).

## A05. [P2] cbsample parses weight() inside quoted if text, changing the sample and overwriting unrelated data

Locations: [package/ado/cbsample.ado:27](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cbsample.ado:27) and `177-180`.
The regex scans and removes the first weight(...) substring from the entire original command, including quoted expressions, without requiring it to be an option. It then enables the output-weight branch. An ordinary string comparison can therefore redirect output to an existing variable and silently modify the if condition.

Reproduction: two rows with `s="other"` and `wt=-99`; run `cbsample if s != "weight(wt)"`. Return code is 0, and wt becomes bootstrap frequency counts (observed 1 and 1) even though no weight option was supplied. With `if s == "weight(wt)"`, the quoted text is removed and the command becomes a zero-observation no-op. Separately, valid whitespace syntax `weight( wt )` is rejected with r(198).

Fix: parse the option through Stata syntax (native bsample itself uses Weight(varname)), then copy its local into a private name if necessary. If manual extraction is unavoidable, split the top-level options while respecting quotes and parentheses, and never scan the if/in expression. Add a regression proving wt is untouched without a real weight option and quoted conditions remain unchanged.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.txt).

## A06. [P2] Unicode ignore() silently strips characters the user did not request

Location: [src/cdestring/cdestring_impl.c:130](/Users/mike/Documents/GitHub/stata-ctools/src/cdestring/cdestring_impl.c:130).
`build_ignore_table` records individual UTF-8 bytes rather than complete characters. `convert_value` drops every byte occurring anywhere in that table. Thus the union of bytes for two ignored characters can accidentally remove an entirely different Unicode character. This can turn invalid text into a valid number without `force`, bypassing the intended nonnumeric-value guard.

Reproduction: create `s="12©"`; `destring s, generate(native) ignore("é¢")` correctly declines to create native, while `cdestring s, generate(fast) ignore("é¢")` succeeds and creates fast=12. © is not in the ignore string. `é` supplies A9 and `¢` supplies C2; © is C2 A9. Native Stata 14+ destring defaults to Unicode character semantics (confirmed in installed base/_/_parse_ignore.ado).

Fix: parse ignore() and each source string as complete UTF-8 code points, retain a cheap ASCII table fast path, and remove only requested whole characters. Alternatively delegate non-ASCII ignore strings to native destring. If supporting explicit asbytes, make it opt-in. Regression tests should include overlapping UTF-8 bytes such as this example and both force/non-force modes.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/probes.txt).

## A07. [P2] cdestring discards variable labels and custom characteristics on conversion

Location: [package/ado/cdestring.ado:300](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cdestring.ado:300) (drop/rename at 304-305); destination creation is 132-163.
Successful outputs are new variables with default metadata. replace drops the original before transferring its label/characteristics, so metadata is permanently lost. Generate outputs also lack the source metadata. Native destring copies the variable label and transfers/copies characteristics, including notes, and keeps the output in the source's position.

Reproduction: x is numeric text with variable label "Original variable label" and char x[custom] "custom metadata"; clone y=x; run native `destring x, replace` and `cdestring y, replace`. Native retains both; y has an empty label and empty custom characteristic. Recorded in probes.log.

Fix: copy labels and characteristics to the staged numeric result before dropping the source, and place outputs consistently with native destring. Preserve existing metadata for replace and copy it for generate. Test labels containing quotes, notes/custom characteristics, multiple variables, and source order. This is a data/metadata integrity issue, not merely display-message parity.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/probes.txt).

## A08. [P2] Valid XLSX character data is silently dropped when expressed as CDATA or split by XML comments

**Locations:** [src/cimport/cimport_xlsx.c:1171](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1171) and `1196-1203`; [src/formats/xlsx_xml.c:469](/Users/mike/Documents/GitHub/stata-ctools/src/formats/xlsx_xml.c:469).

**Trigger:** A worksheet inline string `<t><![CDATA[customer & name]]></t>` or `<t>before<!-- comment -->after</t>`. These are ordinary valid XML representations of string content, and native Stata accepts them. Shared strings using CDATA are also affected.

**Observed impact:** Native Stata imports `customer & name`, while ctools stores no value and, when this is the first column's only content, drops that column entirely (2 columns become 1). For the comment case, native imports `beforeafter`; ctools returns `before`. Both commands return success. Native ASan/UBSan probes independently confirm empty inline/shared CDATA and prefix-only commented strings.

**Cause:** The fast worksheet scanner takes only the span before the first `<` inside `<t>`, then ignores subsequent text segments. The generic parser treats CDATA as a declaration and skips it instead of emitting its content.

**Proposed fix:** Retain the optimized path for an uninterrupted text node, but make the fallback concatenate all character-data segments inside the selected `<t>` node; emit CDATA as literal text and skip comments without losing following text. Decode XML entities only in ordinary text, not CDATA. Apply equivalent behavior to shared-string parsing. Add both fixtures and mixed text/CDATA/comment fixtures to native and Stata comparison tests.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.txt).

## A09. [P2] XLSX text does not round-trip OOXML escapes, and exports delete control bytes

**Locations:** [src/cexport/cexport_xlsx.c:300](/Users/mike/Documents/GitHub/stata-ctools/src/cexport/cexport_xlsx.c:300) (filtering controls), [src/cimport/cimport_xlsx.c:533](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:533) (shared-string entity decoding only), [src/cimport/cimport_xlsx.c:1043](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1043) (inline entity decoding only).

**Trigger / reproduced results:**

1. In Stata, `gen str30 s="a"+char(1)+"b"`. Native `export excel` emits `a_x0001_b` in its string table and native re-import restores the control byte. `cimport` of that native workbook instead returns the literal nine-character escape text `a_x0001_b`.
2. `cexport excel` of that same original value emits plain `ab` and returns success. Both native and ctools imports consequently return `ab`, so the control byte is irretrievably lost in the exported file.
3. Native `export excel` of literal `a_x000D_b` emits `a_x005f_x000D_b` to protect the literal escape-looking text. `cimport` returns that protected serialization literally instead of the original text. The same occurs for literal `a_x0041_b`.

**Cause/impact:** XML entity decoding is insufficient for SpreadsheetML's escaped string representation. The writer simply discards C0 controls other than tab/newline/CR, and the reader never handles `_xHHHH_` escapes or escaped underscores. This affects ordinary Stata-produced workbooks and text values, even when every file is structurally valid.

**Proposed fix:** Add an OOXML string codec in addition to XML escaping: preserve supported controls using `_xHHHH_`, protect literal escape-looking substrings with `_x005F_`, and correctly decode the escapes on import without recursively decoding protected literals. Handle surrogate pairs and define/reject values that cannot be represented in Stata text rather than silently deleting bytes. Apply it consistently to shared/inline strings and labels. Add cross-reader round trips for controls, literal escapes, escaped underscores, Unicode supplementary characters, CR/LF, and ordinary XML entities. Note that the observed native inline-string reader did not interpret a raw `_x000D_` test cell, while native shared-string export/import did use the encoding; tests should cover both representations explicitly and avoid assuming native inline behavior is identical.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/io/probe.txt).

## A10. [P2] crangestat prefix sums return wrong means and sums after a large outlier

Locations: [src/crangestat/crangestat_kernels.c:69](/Users/mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_kernels.c:69); threshold selection in [src/crangestat/crangestat_impl.c:798](/Users/mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:798) and [src/crangestat/crangestat_internal.h:50](/Users/mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_internal.h:50).

**Issue and impact.** Every prefix array is centered at the first nonmissing value in its group. Subtracting a very large origin destroys small values, and subtracting two large prefix totals introduces further cancellation. The variance branch has a stable fallback, but mean and sum do not. Results can therefore change when a group reaches the 64-row optimization threshold, even if a requested window contains only one row.

**Reproduction.** Create 64 observations; set `x=_n-1` as double, replace the first value with `1e16`, and set `key=_n`. Run `crangestat (mean) cm=x (sum) cs=x, interval(key 0 0)`. Each interval contains exactly its own observation, so both results must equal x. Instead 62 of 64 results are wrong: the values 1,2,3,5 produce 0,0,4,8 respectively. A four-row public-command control uses the direct path and is correct. A separate native harness reproduces the faulty prefix kernel directly. No external rangestat installation is needed to establish the singleton-window identity.

**Proposed fix.** Make range sums and means cancellation-aware. At minimum return the original valid value for singleton windows, and fall back to a compensated/local-origin direct calculation when prefix subtraction or restoring the origin loses precision. For scalable repeated windows, use a block/tree representation with stable local sums. Merely switching to long double is insufficient on targets where it has double precision. Test singleton and multirow windows, large outliers, mixed signs, excludeself, missing values, and groups just below/at/above 64 rows.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/range64.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/range64.txt).

## A11. [P2] civreghdfe ignores small and always imposes small-sample VCE/inference

Locations: [package/ado/civreghdfe.ado:30](/Users/mike/Documents/GitHub/stata-ctools/package/ado/civreghdfe.ado:30); [src/civreghdfe/civreghdfe_vce.c:220](/Users/mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_vce.c:220); [src/civreghdfe/civreghdfe_impl.c:1446](/Users/mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:1446). SMALL is declared but never referenced by the ado implementation or passed to the plugin. The native covariance always uses N-K absorbed DF denominators. In no-FE IV on 60 observations, native ivregress's default slope variance is .0221972249922; civreghdfe gives .0229626465436, exactly 60/58 times greater. With small, native and civreghdfe match. ctools default and small are also identical for robust and clustered fits. Help explicitly describes small as requesting small-sample corrections ([package/help/civreghdfe.sthlp:46](/Users/mike/Documents/GitHub/stata-ctools/package/help/civreghdfe.sthlp:46)), and without absorb the command claims ordinary IV support.

Fix: pass an explicit small flag through the input/estimator/VCE layer, automatically enabling it when absorb requires upstream behavior. Apply the correct default/asymptotic and small scaling and distribution consistently to VCE, e(df_r), reported tests, RMSE, and CIs. Validate no-absorb default/small for iid, robust, cluster and GMM.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe.txt).

## A12. [P2] no-absorb civreghdfe discards requested dofminus/sdofminus in covariance

Locations: [src/civreghdfe/civreghdfe_impl.c:1450](/Users/mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:1450). The ternaries force sdof_extra=0 and dofminus_vce=0 whenever G==0. Thus `civreghdfe y (x=z), small dofminus(8)` and `..., small sdofminus(8)` produce the exact same variance (.0229626465436) as baseline. dofminus is nevertheless subtracted from reported df_r, so reported DF and covariance scaling disagree. Help promises both options without requiring absorb (help lines47-48).

Fix: retain both user adjustments for no-FE models, combine with automatic absorb/partial adjustments only when present, and centralize effective-DF calculation so covariance, e(df_r), RMSE and diagnostics use consistent values. Test no-FE with nonzero dofminus/sdofminus for all VCE types.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe.txt).

## A13. [P2] civreghdfe uses centered TSS for noconstant models

Location: [src/civreghdfe/civreghdfe_impl.c:1492](/Users/mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:1492). TSS is always sum(w*(y-mean(y))^2), without any constant/centering flag. With noconstant and no absorbed intercept, it should use uncentered sum(w*y^2), or store both definitions and use the correct one in model fit summaries. Probe native uncentered TSS=5661.1831 while ctools e(tss)=80.301186; ctools then reports R²=-75.662598. The same wrong denominator affects e(mss), adjusted R² and display. Reference native ivregress reports R² missing for this particular example because RSS exceeds its uncentered TSS; the primary indisputable mismatch is the TSS denominator.

Fix: carry explicit/absorbed-intercept state through estimation, compute centered and uncentered totals separately as needed, and use the definition required by noconstant consistently. Add no-constant IV and OLS tests with an outcome far from mean zero.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.txt).

## A14. [P2] civreghdfe rejects estimable small samples by adding regressor and instrument counts

Locations: [src/civreghdfe/civreghdfe_impl.c:555](/Users/mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:555). The hard checks require N >= K_total+K_iv+1 before/after singleton removal. These are overlapping spaces, not additive lost degrees of freedom. Four data rows (y,x,z)=(1,1,1),(4,3,2),(2,2,3),(5,4,4) support exactly identified intercept+one-slope IV with 2 residual DF: native ivregress succeeds with slope1.25 and F=32, while civreghdfe returns r(2001). Pre-collinearity checks also count redundant columns rather than effective ranks.

Fix: check identification and effective ranks after collinearity/absorption and require the relevant N/rank/DF conditions for each requested statistic, rather than summing K and L. Allow otherwise estimable coefficients even when a diagnostic cannot be computed, following reference policy. Add four-row valid-IV and many-redundant-regressor tests.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.txt).

## A15. [P2] PPML residual predictions reject time-series dependent variables accepted by estimation

Locations: [package/ado/cpplmhdfe_p.ado:71](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cpplmhdfe_p.ado:71), [package/ado/cpplmhdfe.ado:346](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cpplmhdfe.ado:346). Estimation accepts ts and saves original dependent term in e(depvar), e.g. L.y. Prediction's `confirm variable L.y` returns r(101) because confirm does not accept time-series operators. Probe `cpplmhdfe L.y x, noabsorb` succeeds; predict mu succeeds, but `predict double rr, response` fails r(101). The same shared branch blocks scores, Pearson, deviance, Anscombe and working predictions.

Fix: resolve e(depvar) with tsrevar in the predictor (or validate the underlying TS expression and use it safely), preserving if/in handling and temporary output transaction. Add lag/difference dependent-variable regression tests for every response-based prediction, including estimates restore.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/probe2.txt).

## A16. [P2] cqreg e(rank) counts omitted and factor base columns

Location: [package/ado/cqreg.ado:312](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cqreg.ado:312). K_full=colsof(b) includes explicit zero placeholders for base levels and collinear regressors, and is posted as e(rank), even though df_m and df_r use the actual kept columns. Duplicate x2=x: native qreg rank2, cqreg rank3 while both df_r=1998. Model with x and three-level i.category: native rank4, cqreg rank5. Thus returned rank disagrees with its own covariance and residual DF; downstream postestimation/model-comparison code relying on rank gets the wrong parameter count.

Fix: use Stata's _post_vce_rank or compute covariance rank consistently with the reference. At minimum use K_keep+1 for full-rank covariance, and account for rank deficiency. Add direct e(rank) assertions in collinearity and factor-variable cases; current validation helper compares df_r but omits rank.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/cqreg.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/estimation/cqreg.txt).

## A17. [P2] cmerge keepusing() fails for supported wildcard/range variable lists

Location: [package/ado/cmerge.ado:655](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cmerge.ado:655), especially 658-664.
Selective `use` expands keepusing's varlist, but the wrapper later loops over the raw option string and asks for each token's single-variable type, never normalizing it to concrete names. Wildcards/ranges are therefore treated as literal/ambiguous names. The documented option is keepusing(varlist), as for merge.

Reproduction: using has id,x1,x2,extra; master has id. `merge 1:1 id using ..., keepusing(x*)` succeeds. `cmerge ... keepusing(x*)` fails with r(111), "x ambiguous abbreviation". Logs in merge_probes.log.

Fix: `unab`/normalize keepusing against the loaded using schema, remove duplicates/key variables as appropriate, and derive all names/types/indices/counts from that expanded list. Cover wildcard, range, abbreviation, and overlapping key/payload selections.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.txt).

## A18. [P2] Empty-master observation-number merge ignores keepusing()

Locations: [package/ado/cmerge.ado:553](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cmerge.ado:553) (loads all variables for merge_by_n), `814-815` (saves the entire using dataset), and `896-898` (appends it).
For 1:1 _n the wrapper loads the whole using file even when keepusing is specified. The normal plugin path selects requested payloads later, but the empty-master shortcut appends the unfiltered saved using dataset, importing unwanted variables. The result schema changes depending on whether the master has zero rows.

Reproduction: using has id,x1,x2,extra; zero-row master has id. Native `merge 1:1 _n using ..., keepusing(x1)` yields id,x1,_merge. cmerge yields id,x1,x2,extra,_merge, r(0).

Fix: filter the saved using dataset to the requested concrete variables (plus any required key/master columns under native semantics), or use the same metadata/template path for empty and nonempty merges. Add zero-row-master _n tests for keepusing and shared/excluded columns.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.txt).

## A19. [P2] Empty-input cmerge deletes filtered rows before reporting a failed assert()

Location: [package/ado/cmerge.ado:931](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cmerge.ado:931).
The empty-input shortcut applies keep() before assert(), while the normal cmerge path and native merge assert first. The assertion return code remains correct, but its documented diagnostic result is not: cmerge says "merged result left in memory" after it has already deleted the rows that explain the failure. It also defeats the explicit exception in the cmerge rollback wrapper that commits the full merged result for r(9).

Reproduction: zero-row master with id, using with three rows; run `merge 1:1 id using ..., keep(master) assert(master)`. Native returns r(9) with all three using-only rows still present. cmerge returns r(9) with N=0. An analogous case with nonempty master and empty using can remove the master rows on assertion failure.

Fix: move the empty-branch assert before keep, like the existing normal path; verify r(9) leaves the unfiltered merged result and successful calls still apply keep. This affects the data retained after an assertion error; it does not produce a wrong successful join. The impact is loss of failure-inspection data and inconsistent behavior in scripts that intentionally catch assertion failures.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/merge_probes.txt).

## A20. [P2] cbinscatter unconditionally deletes a preexisting user frame

Location: [package/ado/cbinscatter.ado:42](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cbinscatter.ado:42) (also fixed-name reuse at 565-567 and cleanup at 581).
Every invocation begins with `capture frame drop _ctools_graph_frame`, before option validation and regardless of nograph. If a user already has data in a frame with that legal name, it is silently deleted. This is not limited to cleanup of a frame created by the current call; no ownership is tracked.

Reproduction: create `_ctools_graph_frame`, put one row with protected=99 in it, return to default frame and run `cbinscatter y x, discrete by(g) linetype(none) nograph`. The command returns 0, and the prior frame no longer exists (`USER_FRAME_SURVIVED=0`). Thus even a computation-only command loses unrelated in-memory data. Calling from a current frame with that name would additionally make graph frame creation collide.

Fix: obtain a unique name with `tempname`, create it only when graphing, and clean up only that call-owned frame. Remove the global fixed-name deletion. Test a preexisting similarly named frame, graph failure, nograph, invalid options, and invocation from a nondefault current frame.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.txt).

## A21. [P2] Singleton by-groups silently disappear from bin output

Location: [src/cbinscatter/cbinscatter_impl.c:722](/Users/mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_impl.c:722).
The code skips every group with fewer than two rows before it constructs bins, even in `discrete linetype(none)` mode where a singleton has a perfectly defined mean and no fit needs estimation. The row remains in e(N) but contributes to no returned bin. This violates the documented discrete mode (one bin per unique x value; [package/help/cbinscatter.sthlp:173](/Users/mike/Documents/GitHub/stata-ctools/package/help/cbinscatter.sthlp:173)) and misrepresents small subgroups in exported/charted results.

Reproduction: x=(1,2,3), y=(10,20,30), g=(1,1,2); run `cbinscatter y x, discrete by(g) linetype(none) nograph`. e(N)=3, but e(bindata) contains only group 1's (1,10,n=1) and (2,20,n=1) bins; group 2's (3,30,n=1) is absent and all remaining matrix rows are missing. Return code 0 and no exclusion warning. This was reproduced directly; third-party binscatter was unavailable for comparison.

Fix: construct bins for nonempty singleton groups and independently skip or mark unavailable fit statistics when a fit is underidentified. Quantile binning needs an explicit singleton path too. If some adjusted-model mode cannot support that group, return a clear error or explicitly document/count exclusions rather than silently discarding data. Add discrete and quantile singleton by-group tests with/without fits and verify bin counts sum to e(N).

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.txt).

## A22. [P2] cpsmatch treats allocation failure as absence of a suitable match

Locations: [src/cpsmatch/cpsmatch_impl.c:695](/Users/mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_impl.c:695); workspace-growth failures in [src/cpsmatch/cpsmatch_kernels.c:93](/Users/mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_kernels.c:93).

**Issue and impact.** Nearest-neighbor and kernel helpers return a failure code when their candidate or match workspace cannot grow. The caller combines `match_rc != 0` with the ordinary no-match case, marks the treated observation off support, and continues. The command can return success with fewer matched observations and a changed ATT because of memory pressure, rather than because of caliper/common-support restrictions.

**Reproduction.** A native ASan/UBSan harness supplies one treated observation and 300 controls, all with propensity score .5, and kernel bandwidth 1. With normal allocation, `rc=0 matched=1 treated_support=1`. Injecting failure only in the kernel workspace reallocations produces `rc=0 matched=0 treated_support=0`. The dataset and matching rules are unchanged. This exercises the production matching entry point with mocked SPI transport; it is not a real process-wide out-of-memory experiment.

**Proposed fix.** Collect a parallel-safe fatal-error flag separately from genuine no-match status. If any worker fails, stop/finish outstanding work, clean up, return the memory-error code, and do not publish matching output. Keep the ado output transaction separate so prior data survive the failure. Add fault injection at initial workspace allocation and every growth point for nearest-neighbor, ties, and kernel modes, with serial and OpenMP variants.

**Evidence:** [ASan/UBSan fault-injection harness](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/match_oom.c).

## A23. [P2] cpsmatch deletes a supplied score variable before reading it

Locations: [package/ado/cpsmatch.ado:133](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cpsmatch.ado:133).

**Issue and impact.** The wrapper drops all conventional output names before staging its input variables. If a supplied score is named `_pscore`, it deletes that variable and then attempts `gen double _pscore = _pscore`. Other collisions with treatment/covariate/outcome inputs have the same unsafe ordering. These names are valid Stata variables and the syntax accepts them as inputs.

**Reproduction.** On four rows with a 0/1 treatment, copy a valid score variable to `_pscore` and run `cpsmatch t, pscore(_pscore) outcome(y)`. It returns r(111), and `_pscore` is no longer present. This is an error accompanied by lost input, not just a rejected name. No psmatch2 compatibility comparison is claimed; the input-preservation failure is directly observed.

**Proposed fix.** Validate collisions before any destructive operation, or stage required inputs into private temporary variables first. Construct outputs under temporary names and commit them only on success. If reserved names are intentionally unsupported, reject them explicitly while leaving data unchanged. Test score, treatment, covariate and outcome collisions, plus an estimation error after existing output variables have been created.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/root_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/root_probes.txt).

## A24. [P2] cpsmatch accepts 0/2 treatment coding but omits ATT and miscounts support

Locations: [package/ado/cpsmatch.ado:36](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cpsmatch.ado:36); [src/cpsmatch/cpsmatch_impl.c:172](/Users/mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_impl.c:172).

**Issue and impact.** The advertised binary 0/1 check verifies only that two distinct values exist. C matching treats every nonzero value as treated, but `_treated` retains the original value and the ATT/reporting code selects `_treated == 1`. Thus a coding error passes validation and produces internally inconsistent successful output.

**Reproduction.** Use `(t,ps,y)` rows `(0,.1,2),(0,.2,3),(2,.15,7),(2,.25,8)`. `cpsmatch t, pscore(ps) outcome(y)` returns 0 and reports two treated/two matched observations in r(), but stores `_treated=2` and returns no ATT. The support display also fails to count these treated rows correctly.

**Proposed fix.** Enforce exactly 0 and 1 in the effective matching sample before creating outputs. Alternatively, if arbitrary nonzero treatment is intended, normalize once and use that definition throughout validation, matching, returned variables and ATT. Do not mix the conventions. Test 0/1, 0/2, 1/2, negative codes, missing values and sample restrictions.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/root_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/root_probes.txt).

## A25. [P2] cpsmatch does not validate matching-option domains

Locations: [package/ado/cpsmatch.ado:43](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cpsmatch.ado:43); [src/cpsmatch/cpsmatch_options.c:113](/Users/mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_options.c:113).

**Issue and impact.** Syntactically numeric options are accepted without their required domain checks, and an unknown kernel name silently maps to Epanechnikov. Invalid settings become successful scientific results rather than errors.

**Reproduction.** On four ordinary 0/1 observations, `neighbor(0)` and `neighbor(-1)` both return 0 with every treated row unmatched; `kernel bwidth(0)` does the same. `caliper(-.1)` succeeds and behaves as though no caliper were requested. `kernel kerneltype(nonsense)` succeeds using the default kernel. These are distinct from legitimate no-match results under valid parameters.

**Proposed fix.** Validate neighbor >=1, finite nonnegative caliper, strictly positive finite bandwidth when relevant, and membership in the supported kernel names. Reject contradictory method switches rather than silently selecting one. Perform these checks before dropping or creating any outputs, and repeat essential invariants at the C boundary. Add invalid/boundary-value tests that check both error codes and unchanged data.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/root_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/shared/root_probes.txt).

## A26. [P2] `cencode` label-file writer still bypasses the Windows Unicode path layer

**Locations:** [src/ctools_hash.c:231](/Users/mike/Documents/GitHub/stata-ctools/src/ctools_hash.c:231); [package/ado/cencode.ado:191](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cencode.ado:191); callers [src/cencode/cencode_impl.c:342](/Users/mike/Documents/GitHub/stata-ctools/src/cencode/cencode_impl.c:342) and `738-746`.

**Issue:** The shared label writer calls `fopen(filepath, "wb")` directly. The ado wrapper passes Stata's `tempfile` path as UTF-8 through `__ce_labelfile`; this reaches the narrow Windows CRT API rather than the UTF-8-to-UTF-16 `ctools_fopen` adapter introduced for the other file paths. A Windows account or configured temporary directory containing non-ASCII characters therefore makes ordinary `cencode source, generate(dest)` fail when the process uses a legacy Windows code page. It reports `cencode: failed to write label file`, r(920), even though the cause is a valid filename rather than memory exhaustion. Long temporary paths also miss the adapter's extended-path support.

**Evidence/confidence:** Source-confirmed complete public-command call path; not executed on Windows in this audit. [src/ctools_file.h:3](/Users/mike/Documents/GitHub/stata-ctools/src/ctools_file.h:3) explicitly establishes UTF-8 for all plugin-owned paths, and `42-49` supplies the needed `_wfopen` conversion. Current Windows regression `validation/native/harnesses/oct21_windows.c` tests that adapter/output publication, but never calls `ctools_label_write_file`. This is a missed cencode endpoint of the previous Windows path migration; the previously fixed CSV/binary endpoints were not re-reported.

**Proposed fix:** Include `ctools_file.h` in `ctools_hash.c` and use `ctools_fopen` for the label file. Extend the Windows native test to invoke the actual label writer in a Unicode directory under a legacy code page, and add a licensed cencode test with a Unicode temp directory. Test empty and nonempty label vectors, spaces, Unicode, and a long path. Return a file-write error rather than an allocation error if practical.

**Evidence:** Source-confirmed; Windows execution unavailable.

## A27. [P3] linetype(connect) produces unconnected scatter points

Locations: [package/ado/cbinscatter.ado:76](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cbinscatter.ado:76), `430-431`, `503`.
The option sets linetype_num=0 (correctly disabling the regression fit), but neither single-group nor by-group scatteri specification ever adds a connection option. lcolor() only styles a line; it does not cause scatter points to be connected. Help explicitly promises that connect joins the bin means.

Evidence: a test-only shadow twoway program captured the exact generated graph specification for `cbinscatter y x, nq(2) linetype(connect)`: `(scatteri 15 1.5 30 3, mcolor(navy) lcolor(navy))`, with no connect() or connected recast and no fit plot. Command returned 0. This verifies generated syntax; the probe did not render the resulting graph.

Fix: append `connect(l)` (or recast to connected) when linetype_lower is connect in both branches, ensuring points are ordered by x. Test generated graph syntax plus a graph smoke test for one and multiple groups.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.txt).

## A28. [P3] reportreg is accepted but silently does nothing

Locations: [package/ado/cbinscatter.ado:217](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cbinscatter.ado:217) and [src/cbinscatter/cbinscatter_impl.c:127](/Users/mike/Documents/GitHub/stata-ctools/src/cbinscatter/cbinscatter_impl.c:127).
The documented reportreg flag is parsed and copied into config.reportreg, but no subsequent code ever reads it or reports the underlying regression. Users receive success without the promised table.

Reproduction: `cbinscatter y x, nq(2) nograph reportreg` returns 0 with no regression output; bins_probes.log shows the empty output between the command and RC marker. Repository search finds only config default, parse assignment, and struct declaration, with no consumer.

Fix: implement the regression report using the actual fitted coefficients/statistics and label separate groups correctly, or explicitly reject the option/remove the unsupported promise until implemented. Add an output assertion for reportreg and ensure ordinary runs remain quiet.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_probes.txt).

## A29. [P3] Stored by/control metadata references deleted internal variables

Location: [package/ado/cbinscatter.ado:625](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cbinscatter.ado:625); original-name preservation at 173 and 188, rewrites at 176 and 191.
The wrapper replaces by with a temporary dense group-ID variable and expands factor controls into temporary variables. It correctly keeps original copies in by_orig and controls_original, but ereturn posts the rewritten locals. Stata then drops those temporary variables as the ado returns, leaving e(by) and factor-variable e(controls) invalid. Consumers cannot identify the grouping variable or reconstruct controls from the documented stored metadata.

Reproduction: 20 rows with x=_n, y=x^2+mod(_n,3), g=mod(_n,2), c=mod(_n,3); run `cbinscatter y x, nq(2) by(g) controls(i.c) nograph`. e(by) is __000005; e(controls) is __000002 __000003 __000004. Every recorded name fails confirm variable immediately after return. Fresh plugin probe and output are `bins_metadata.do` and `bins_metadata.log` in the same probe directory.

Fix: post e(by)=by_orig and e(controls)=controls_original (or a deliberately documented permanent expanded-name representation), while keeping internal numeric mappings private. Test a numeric/string by variable and factor interactions/ordinary controls, verifying documented metadata refers to the user's model specification.

**Evidence:** [probe](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_metadata.do) and [results](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/data/bins_metadata.txt).

## Recommended implementation order and coverage improvements

1. Fix the four P1 findings first: model F hypotheses, fitted clustered quantile covariance, validation of every consumed XLSX XML part, and SHP advertised-length validation. Add the smallest reproductions to registered native/licensed suites before changing the implementations.
2. Address unintended modification and loss next: quoted bootstrap parsing, Unicode destring, metadata transfer, frame ownership, matching input collisions, and matching allocation failures. Verify failure paths preserve prior data and do not publish partial successful-looking results.
3. Correct estimator option propagation and result metadata with a shared test matrix covering constant/no-constant, absorb/no-absorb, weights, collinearity, robust/cluster, and factor/time-series terms. Assert model F, e(rank), e(df_r), TSS and prediction behavior as well as coefficient/VCE agreement.
4. Fill the observed boundary gaps: group sizes 63/64/65, empty merge inputs, singleton bins, quoted expressions containing option-like text, Unicode code points with overlapping bytes, and valid XML text represented in different ways. Round trips through the same flawed reader/writer are insufficient; retain independent-reader or exact-content checks.
5. Keep platform certification separate from source confidence. Run A26 on real Windows with a Unicode temporary-directory path and the legacy code-page setup used by the Windows fixture. Install/configure the documented OpenMP dependency for the remaining native runs, and complete the trusted licensed release matrix before publication.

For reproducibility, record versions of installed statistical reference packages and their transitive dependencies alongside the Stata/plugin identity. The current tests can be sensitive to reference behavior and local package availability; this is a test-infrastructure recommendation, not a separately counted runtime bug. Minimum supported Stata/OS execution should remain an explicit release criterion.

## Evidence and reproduction notes

The accompanying evidence directory preserves probe source, complete logs, fixture generators, and native harnesses. These are audit artifacts, not modifications to production code or registered regression tests. Some probes intentionally trigger errors or intercept graph construction; a zero process exit alone is not an assertion of correctness. Read the recorded result markers and the per-finding comparisons.

The scripts retain the original absolute repository and temporary build paths from this session. To replay elsewhere, set their adopath/plugin path to a newly built candidate and adjust fixture-directory constants. The I/O fixture generator recreates the valid/corrupt XLSX and SHP files. No generated binary test executables or private user datasets are included.

[Evidence index and native reproduction commands](/Users/mike/Documents/GitHub/stata-ctools/docs/audit-evidence/2026-10-22/README.md)
