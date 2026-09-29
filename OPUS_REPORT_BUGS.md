# ctools bug report — September 27, 2026 (Claude Opus 5.5)

This report lists defects found by a full audit of the ctools working tree. That tree includes the uncommitted SEP24 repairs and the new modules that no earlier audit covered: `src/io/` (SAS/SPSS/XPORT/XLS/SHP/dBase), `cipolate`, `csplit`, `crangejoin`, and `ctools_order`. Each entry gives the location, the evidence, a trigger or reproduction, the impact and a proposed fix. The audit itself changed no source, ado, help, build or validation file. The eleven P0 findings were fixed afterwards (see [P0 fix status](#p0-fix-status)), and then the 57 P1 findings (see [P1 fix status](#p1-fix-status)).

<a id="p0-fix-status"></a>
## P0 fix status (2026-09-27)

All eleven P0 findings are fixed in the working tree (uncommitted), and `build/ctools_mac_arm.plugin` is rebuilt. Each fix was checked in Stata 18.5 MP on macOS arm64 with the original reproduction and the cases listed. The full `validate_all.do` suite gives the same result as before the fixes: 4,108 passed, with the same three pre-existing failures (the `audit_p1`, `audit_p2` and `sep24` import steps). Windows was checked by cross-compiling (MinGW) and by a macOS test build that mimics `_aligned_malloc`; nothing was run on Windows or Linux. No validation script was changed.

| ID | Fix | Verified by |
|---|---|---|
| SORT-1 (+SORT-6) | Co-rank predicate `>=` in the three timsort parallel merges; the split is sized from the team that performs it | 13 `csort, algorithm(timsort)` cases (N = 200k–1M; 2, 3, 4, 8 and default threads; string and multi-key; `nosortedby`; `stream(1)`) keep every row in stable order. The audit's native suites: timsort 0 failures (was 82 medium, 134 large); random team sizes 0 of 20 bad (was 16 of 20) |
| IMP-1 (+IMP-6) | Strict mode splits only at a newline outside quotes, using the exact quote parity from the data start (counted in parallel). The strict row finder is rewritten (the SIMD version skipped boundaries); `nobind` uses plain line boundaries; boundaries are monotonic, and a search already passed by the previous one is skipped | Strict 60k, 200k-CRLF and long-field files, and the audit's `strict_big*`/`nobind_big*` files (`nobind_big2` gave 272,547 rows, now 150,000): row counts and data match native, and multi-chunk output equals `threads(1)`. `varnames(2)` with strict quotes now matches native |
| EXP-1 | Writes loop until complete (≤ 1 GiB per call, retry on EINTR, fail on 0/−1); each wave's byte count is checked; `ftruncate`/sync failures are fatal; chunk buffers are capped at 16 MB | 600,000 × 220 export (2.37 GB) now succeeds with default threads, byte-identical to `threads(4)`. A short-write/EINTR injection test passes (and fails on the old code). Wide-row exports, including `mmap`, are byte-identical to `export delimited` |
| REG-1 | `unab` the `absorb()` list; the plugin refuses a varlist whose length does not match its layout | `absorb(fe*)` and `absorb(fe1-fe2) resid()` match reghdfe; user variables untouched |
| REG-2 | `vce(cluster a b)` exits r(198) instead of silently clustering on `a`; README and help corrected (multi-way clustering is not implemented) | Error with and without weights and `resid()`; no variables overwritten; one-way clustering unchanged |
| REG-3 | `ctools_aligned_free` for the filtered `obs_map` in creghdfe (3 sites), cpplmhdfe (10) and civreghdfe (19; not in the original finding) | macOS test build with a Windows-style offset allocator: the old code aborts in `creghdfe` ("not an allocated block"); the fixed code passes the full suite |
| IV-1 | Base/omitted re-insertion uses the same rule that selected the plugin columns, with a consistency check | `i.a#c.x`, `i.a#i.b` and `i.a##i.b` (DiD, which gave r(504)): b and V match ivreghdfe to 1e-14 |
| QREG-1 | Sparsity/density failures return qreg's r(498) (920/2000 where qreg does); auxiliary fits use `maxiter()`; no clamping | `auto`, q(.05): r(498) like qreg. Log-normal N = 200k, q(.95): se(x1) = 6792.0567 (qreg 6792.0567; was 0.00106). Robust and residual methods at q = .02/.05/.5/.95: same rc as qreg, same SE at the median |
| BIN-1 | Exact `xtile` cutpoints for every N (radix sort of x); the binsreg adjustments reuse those bin ids; 64-bit cutpoint index and overflow-safe weighted search | Log-normal, integer-tied, outlier, weighted, `by()` and `nquantiles(3/7/9/12)` cases (N a multiple of nq) match `xtile` bin means to 1e-15 |
| SMP-1 | Seed passed as two exact 31-bit halves; strict parsing of `seed`, `count` and `n` | 30 seeds give 30 distinct `csample` draws; 50 `cbsample` draws are all distinct; `by()` draws are distinct; same seed still reproducible |
| RNG-1 | Prefix sums over 64 fixed logical chunks with `omp for` (no per-thread-id chunks) | Unbalanced panels match rangestat; results are identical for `threads(1)`, `threads(3)` and default; audit harness 0 wrong cells (was 11,005 of 12,000) |

Found while fixing (not fixed; outside the P0 scope):
- Native `import delimited` stores an embedded CRLF inside a quoted field as CR only; `cimport` keeps CR+LF (also in single-chunk files).
- `creghdfe` leaks `obs_map`, `cluster_raw_values` and `weighted_counts_orig[]` on its r(2001) all-singletons return and on allocation-failure returns (`src/creghdfe/creghdfe_regress.c`).
- `src/ctools_data_io.c` allocates the zero-observation string pointer array with `calloc` but frees it with `ctools_aligned_free` (unreachable today: loads return early when no rows pass the filter).

<a id="p1-fix-status"></a>
## P1 fix status (2026-09-28)

All 57 P1 findings are fixed in the working tree (uncommitted). IMP-6 was already fixed with IMP-1 in the P0 round. `build/ctools_mac_arm.plugin` is rebuilt with the pinned OpenMP runtime (macOS 11 minimum; dependency contract checked).

**How the fixes were checked:**
- **Stata version:** checks ran in StataNow 19.5 MP, which the project's `stata` alias now starts; the P0 round used Stata 18.5.
- **Comparison do-files:** for each area, a do-file compared ctools with the command it replaces (results below).
- **Full suite:** `validate_all.do` gives 4,161 passed, 0 failed, 48 expected skips.
- **I/O parity runners:** all pass: Excel options 196 + 15, formats 183, parity 32, delimited options 227.
- **Native tests:** every native test tied to these fixes passes. CI now runs all `validation/test_*_native.py`.
- **Windows:** the tree cross-compiles with MinGW.
- **Independent review:** three reviewers checked the changes and their findings are fixed. Two of them (32-bit column offsets in cbinscatter; FE-code overflow in `remap_and_count`) were fixed by a concurrent memory-safety session.

**Caveats on these results:**
- **Concurrent edits:** that concurrent session changed source and validation files in the same working tree while these checks ran, so suite totals reflect its edits too.
- **Validation scripts:** apart from the BLD-1/BLD-2 test ports (which were themselves findings), no validation script was changed for these fixes. The earlier "pre-existing" failures (`audit_p1`, `audit_p2`, `sep24`) no longer occur with the current scripts.

| ID | Fix | Verified by |
|---|---|---|
| SORT-2 | Counting sort checks the key range in double (keys within ±2^53, range below 10^6) before any integer conversion; other keys go to LSD radix | The report's cases (±1e30 at N = 8, 1,000 and 60,000; one-sided 1e19) now sort correctly instead of crashing Stata. `cwinsor`/`crangestat` `by()` on such keys match `egen group()`. 647 ASan/UBSan harness runs |
| SORT-3 | IPS4o uses its parallel partition only with two or more threads | `csort ..., algorithm(ips4o) threads(1)` and `cmerge ..., threads(1)` at 30,000 and 100,000 rows match `sort`/`merge` |
| SORT-4 | MSD radix skips a prefix shared by all keys in one step; past depth 32 it finishes with a stable merge sort instead of insertion sort | 200,000 URL-like keys with a 38-byte common prefix: 0.02 s (the harness took 20 s for 160,000 rows before) |
| SORT-5 | IPS4o and sample-sort string radix loop over shared bytes and recurse only into the smaller buckets (depth ≤ log2 n) | str2045 keys sharing 2,044 bytes, and the "a", "aa", … chains, sort correctly with every algorithm (they crashed Stata); `cmerge` on str2045 keys works |
| SORT-7 | `csort ..., nosortedby` clears Stata's sort marker (one value is changed with `replace` and restored through Mata) | `: sortedby` is empty afterwards, and a later `csort` sorts again |
| MRG-1 | `update`/`update replace` code matched rows 4 or 5 as `merge` does, including extended missing values (`.` ← `.a` and `.a` ← `.` are code 4); `keep()`/`assert()` accept all of merge's words | `update` and `update replace` on the report's data and on a 10-row case with `.a`/`.b`: data and `_merge` identical to `merge` (`cf`). `keep()` words and codes and `assert()` give the same N and return codes |
| MRG-2 | New using variables come from a zero-row template of the using data appended to the master, so formats, variable labels, notes and characteristics follow merge's rules; the using notes of shared variables are appended to the master's | Formats, variable labels, notes, characteristics and `_dta` notes match `merge`, with and without `nonotes`/`nolabel`/`update` |
| MRG-3 | Master value-label definitions win; other using definitions are copied even if no new variable uses them; under `update` a shared variable without a master label takes the using label | Value labels and label definitions match `merge` in the same four option sets |
| MRG-4 | `sorted` is verified, including when either dataset is empty (r(5) "master/using data not sorted"). A sort marker on the keys is only a hint: stale markers are re-sorted, not trusted | Unsorted master with `sorted`: r(5) as with `merge`, with a non-empty and with an empty using |
| MRG-5 | strL variables in the master, or in the using keys or kept variables, hand the command to `merge`, run under the caller's version. A failed `assert()` leaves the merged result in memory, as `merge` does, on both code paths | strL merge identical to `merge`. `keep(match) assert(match master)` gives r(9) with the same data in memory as `merge` |
| DSTR-1 | A variable with any nonnumeric value is left untouched unless `force` is given | Refusal, messages and types match `destring` (116 parity checks) |
| DSTR-2 | `percent` divides every converted value by 100 | Parity checks pass |
| DSTR-3 | Infinite results and values at or above maxdouble are nonnumeric (stored as `.` with `force`). Finite values below −maxdouble stay numbers, as `real()` keeps them. `1d3`-style exponents are read as `real()` reads them | Parity checks pass except by design: with `force`, `destring` stores values such as 9e307 as raw doubles in the missing range (they list as `.e_`/`.z_`); `cdestring` stores `.` |
| DSTR-4 | `dpcomma` treats "." as nonnumeric and the first "," as the decimal point, as `destring` does | Parity checks pass |
| IMP-2 | Delimiter detection counts candidates over the first 50 non-empty lines including the header; ties break comma > tab > pipe > colon > semicolon, as native does | 131 of 135 comparisons with `import delimited` pass; the four others are outside IMP-2..5 (see below) |
| IMP-3 | Blank CRLF lines are skipped like blank LF lines | As above |
| IMP-4 | An explicit variable list renames the columns without forcing `varnames(nonames)` | As above |
| IMP-5 | `rowrange()` counts physical lines as native does for headerless files | As above |
| IMP-6 | Fixed with IMP-1 (P0) | P0 |
| IMP-7 | Windows opens files with `CreateFileW` (UTF-8 converted to UTF-16) and full sharing flags | Cross-compile only |
| EXP-2 | Excel cells are written with the shortest round-trip decimal representation (Schubfach) | `cexport excel` → `import excel` round trip: 0 mismatches (16-digit IDs, fractions, 1e-300). In a harness, 71 million values round-trip exactly and are shortest |
| REG-4 | The robust/cluster F test uses the kept coefficients by position | REG-4/5/6 do-file: 0 failures (F with an omitted first regressor matches `reghdfe`) |
| REG-5 | Singleton removal counts frequency weights | Random fweight panels: N, singletons, df_a and SEs match `reghdfe` |
| REG-6 | Singleton removal iterates to convergence | A 1,000-link chain: N, singletons and SE match `reghdfe`, for creghdfe and civreghdfe |
| REG-7 | 64-bit column offsets in the shared solver, VCE, PPML and cqreg code | Code review; harness above 2^31 elements |
| IV-2 | Tested instruments are validated before `Z_rest` is built; a tested variable that was dropped gives r(198), as ivreghdfe does | civreghdfe P1 do-file (IV-2..IV-7): 45 of 45 checks pass against ivreghdfe |
| IV-3 | `orthog()`/`endogtest()`/`redundant()`/`partial()` are mapped to expanded columns through per-column flag vectors | As above (factor-variable and time-series cases) |
| IV-4 | Regressors and instruments absorbed by the fixed effects are omitted (norm test before equilibration). As in ivreghdfe, an absorbed regressor that is not constant within a single factor still counts in K, the residual df and the J df | As above: K, df_r and SEs match ivreghdfe exactly in the report's case |
| IV-5 | The `partial()` projection is weighted | As above (aw/fw/pw) |
| IV-6 | Collinear `partial()` variables are dropped one by one, not all | As above |
| IV-7 | Hansen J uses (w·e)² for aweights/pweights | As above |
| IV-8 | New module `civreghdfe_omega.c`, a port of ivreg2's `m_omega`, builds every S matrix. The C-statistic uses the S block of the full model, the endogeneity test re-estimates the exogenous model, and the redundancy test uses ranktest's LM with the chosen VCE | A second do-file: 156 of 157 cases match ivreghdfe (b and V ~1e-13, statistics ≤ 1e-11) across iid, robust, cluster, two-way, AC/HAC (all kernels), Driscoll-Kraay, Kiefer, aw/pw/fw, one and two endogenous regressors, 2SLS/LIML/Fuller/k-class/gmm2s/cue, and few-cluster cases. The one mismatch is intended (below) |
| IV-9 | Kleibergen-Paap rk LM/Wald for any number of endogenous regressors via ranktest's algorithm; Kiefer uses iid statistics as ivreg2 does | As above (two-endogenous section) |
| IV-10 | Driscoll-Kraay lags are time differences (t − tmin)/tdelta | As above (unbalanced panel with gaps) |
| IV-11 | `bw()` without `robust` is the AC estimator; `kernel()` without `bw()` is r(102) | As above |
| IV-12 | The tsset panel and time variables go to the plugin; kernel options require tsset (r(111)); lags pair observations by time within panels; ivreg2's other kernel errors (r(5), r(101), r(198)) | As above, including 16 error cases with ivreghdfe's return codes |
| IV-13 | gmm2s and cue use the full S (two-way, HAC, Driscoll-Kraay, Kiefer) with ivreg2's dof; CUE's objective uses (w·e)² | As above (cue b within 1e-8: different optimizers) |
| QREG-2 | `denmethod(residual)` implements qreg's residual method (drops the K basis observations, no clamping); residuals are set to zero only within rounding error | 113 comparisons with `qreg` pass. An uncentered year trend matches `qreg` (SE .0045800529) |
| QREG-3 | `denmethod(kernel)` runs the kernel method (all eight kernels; iid, robust, cluster) | As above; kernel `vce(iid)`/`vce(robust)` V within 1e-8 of `qreg` |
| BIN-2 | Classic fit lines are estimated on residualized data and shifted to the plotted means | 80 of 80 checks against `binscatter`, `binsreg`, `reghdfe` and exact references pass |
| BIN-3 | `method(binsreg)` lines use a joint fit with the controls and absorbed effects, evaluated at the control means | As above |
| BIN-4 | Covariate-adjusted binsreg means go to the correct bins when quantile bins are empty (ties, mass points) | As above |
| BIN-5 | `absorb()` codes are remapped (negative, fractional, 11-digit and string codes) | As above |
| BIN-6 | 64-bit quantile arithmetic; weighted cutpoints | As above (bins with fweights equal those of expanded data) |
| XIO-1 | New XLSX worksheet scanner: no dependence on `<dimension>`; namespace prefixes, either quote style, missing `r=`, comments; storage grows on demand; parallel chunks store only their own rows (serial rescan otherwise) | Comparison do-file: 64 of 73 outcomes as expected; the other 9 are the deliberate differences below. Serial and parallel scans agree on 1,200 random sheets under ASan/UBSan; ThreadSanitizer is clean on a sheet with duplicated row numbers |
| XIO-2 | Shared strings are read in full and entity-decoded once; phonetic runs are ignored; sheet names and relationship targets are entity-decoded ("P&L") | As above |
| XIO-3 | `.xls` formula results keep their type (text, boolean, error, empty) | As above |
| XIO-4 | dBase numbers are written the way Stata writes them, with native field widths (byte 3, int 5, long 10); a field widens by one column only for values native would truncate | As above; files byte-identical to `export dbase` apart from header byte 29 and the final byte |
| XIO-5 | Worksheet parts are resolved through the workbook relationships; without `sheet()` the first worksheet is imported; chart, dialog and macro sheets give r(601); a workbook of chart sheets imports as empty | As above |
| XIO-6 | One unit per column (datetime if any cell has a time, else date; every numeric cell converted); booleans are ignored when choosing formats; built-in number formats map to Stata formats as native does | As above |
| PSM-1 | `cpsmatch` uses psmatch2's estimation sample (marksample, outcome, covariates or pscore, `e(sample)`) | ATT and support counts match `psmatch2` with a missing outcome and with missing covariates |
| WIN-1 | The SIMD trim path keeps extended missing values | `.a` values survive `trim` |
| WIN-2 | Shared quickselect with a three-way partition | 1M rows with 50% zeros: 0.009 s (was 21.7 s). cwinsor do-file: 0 failures |
| RNG-2 | Windows with rows but no nonmissing values give count 0 and sum 0. Windows without rows, including inverted intervals, give missing, on both code paths | Matches `rangestat` for all-missing windows, groups of 63/64/65 rows, and an inverted interval at 20 and 80 rows |
| RNG-9 | Variable positions come from Mata `st_varindex` | A 32-character variable name no longer breaks `crangestat` |
| BLD-1 | CI path filters and the package step include the codec licence and every non-plugin manifest file; the release-gate test copies them | `test_release_gate.py` passes |
| BLD-2 | The stale native tests are ported to the current cimport API, and CI now runs every `validation/test_*_native.py` | All native tests pass locally with clang |

Also fixed while reviewing the P1 changes:
- **cmerge sort marker:** after `update`, cmerge now clears the sort marker when it rewrote master values.
- **Sort markers in other commands:** an in-place `cwinsor` clears the marker when it changes a sort variable.
- **`ctools_invsym`:** it now reproduces Mata's `invsym()` (dynamic largest-pivot order, a per-column 1e-9 relative and 1e-19 absolute tolerance, continuing past dependent columns). The old version was wrong for any matrix with off-diagonal entries.

Remaining differences (deliberate, or outside P1):
- **civreghdfe:**
  - `cue` with a rank-deficient S stops with r(506); ivreghdfe fails inside its optimizer with r(430).
  - Without `absorb()`, civreghdfe always applies small-sample corrections (pre-existing).
  - `xtset` with only a panel variable, plus `bw()`: civreghdfe gives r(111); ivreghdfe gives r(451).
  - When an absorbed regressor is not removed exactly, ivreghdfe estimates a meaningless coefficient for it; civreghdfe omits it.
- **XLSX/XLS — native loses data or rejects valid files; ctools keeps the data:**
  - Rows and cells without `r=` after explicit references: native drops, duplicates or moves values.
  - Characters beyond U+FFFF: native keeps only their low 16 bits.
  - Long cells: native cuts every cell at 32,766 bytes, even mid-character.
  - Worksheet targets starting with `./` or `../`, or outside `xl/`: native lists these sheets without a name and imports them as empty.
  - First tab a chart, dialog or macro sheet: native imports an empty dataset when `sheet()` is not given; ctools imports the first worksheet, as Stata's documentation describes.
- **dBase:** header byte 29 (native 0x57) and the final byte (native 0x0A, ctools 0x1A) differ.
- **cimport delimited:**
  - `emptylines(include)`: native keeps one more trailing blank line.
  - After a variable-list error, native leaves no data in memory; cimport leaves the dataset unchanged.
- **cqreg:**
  - `vce(cluster) denmethod(kernel)`: SEs can differ from `qreg` by about 1% when basis observations have residuals of exactly zero. The sign of a zero residual decides their score; this is rounding-level.
  - Residual method on tied data that `qreg` computes with rounding noise: `qreg` reports tiny SEs, and cqreg stops with r(498).
- **`validation/validate_cqreg_performance.do`, line 42:** it expects `denmethod(residual)` to succeed on `y = mod(_n,3)-1`, but native `qreg` returns r(498) there, and so does cqreg now. The script was left unchanged pending a decision.
- **Not verified:**
  - `cbinscatter, method(binsreg) absorb()` may differ from binsreg in the line's level when singleton groups exist (binsreg's regression drops singletons; cbinscatter does not).
  - IPS4o re-partitions an all-equal bucket up to depth 20 (about 1 s per 100,000 duplicated str2045 keys).
  - Merge sort is slow on long shared prefixes.
  - Counting and sample sorts start two threads under `threads(1)` for 50k–100k rows.

## How the audit was done

1. **Code review by subsystem.** Eighteen parallel reviews covered the core runtime, the data I/O layer, the sort engines, every command family, the import/export formats, and build/packaging/validation. Each review read every assigned file in full. Suspected defects were then confirmed with **native C harnesses**, which compile the real `src/` files against a mocked Stata plugin interface. The harnesses ran randomized differential tests against reference implementations, fault injection, and AddressSanitizer/UndefinedBehaviorSanitizer runs.
2. **Differential testing in Stata.** Stata 18.5 MP on macOS/Apple Silicon ran against the checked-in `build/ctools_mac_arm.plugin`, which was built after the last source change. Randomized ctools-vs-native comparisons covered:
   - cipolate, csplit and crangejoin (which match their references in all tested cases);
   - csort, cmerge, cencode/cdecode/cdestring, cwinsor and crangestat;
   - cimport and cexport (CSV and Excel).

   Most high-impact reviewer claims were also reproduced directly, against native Stata or the installed references: reghdfe 6.13.1, ivreghdfe/ivreg2/ranktest, ppmlhdfe, binscatter/binsreg, rangestat, rangejoin, winsor2 and psmatch2. The exceptions rest on harness or code evidence and are marked as such in the index. They are platform-specific (Windows/Linux) or need multi-gigabyte or very large inputs (for example EXP-1, REG-3, REG-7, BIN-6).
3. **Consolidation.** Duplicate reports from different reviewers are merged, and related minor issues are grouped into single P3 entries, so one entry may cover several individual defects. Items fixed by the September 20/22/24 audits are not repeated. Where a finding shows one of those fixes to be incomplete or regressed, the entry says so, citing the old item number (for example "SEP24 #6").

**Evidence levels** (the index shows the strongest level for each entry):
- *Stata-verified*: reproduced in Stata; the numbers quoted come from those runs.
- *Harness*: reproduced with a native C harness running the real code.
- *Code*: established by reading the source and, where relevant, the native command's own ado or Mata source. "Plausible" marks the few items whose runtime effect could not be demonstrated.

**Priority scale:**
- **P0 — fix immediately.** Silent data corruption or loss, a Stata crash, or grossly wrong statistical output in common, default or documented usage.
- **P1 — fix before the next release.** Impact of the same kind but with narrower triggers, or failures of core functionality in common scenarios.
- **P2 — next correctness/robustness cycle.** Narrow-trigger correctness issues, divergence from the native command in edge cases, error-path problems, and performance cliffs.
- **P3 — minor or latent.** Documentation, cosmetic and maintenance issues.

**Limits.**
- Only macOS/arm64 was executed. Windows, Linux and Intel behaviour (REG-3, IMP-7, EXP-10, XIO-12) is established by code reading.
- The `oldstata` clock-shift wrapper invalidates timers near session start, so timings are reported only from runs where they were reliable.
- Vendored libraries (ReadStat, libxls, miniz, libdeflate) were reviewed at their integration points and for local patches only.
- **Stata crashed twice during testing** (SORT-2 and SORT-5). The first crash has macOS crash report `stata-mp-2026-09-27-155502.ips`.
- **Concurrent edits.** Another session on this machine was editing the repository while this audit ran. Between 16:37 and 17:10 it changed `Makefile`, `validation/test_p2_native.py`, `validation/test_transport_scheduling.py`, `docs/PERFORMANCE_TRANSPORT.md` and `src/ctools_data_io.c`, and it rebuilt `build/ctools_mac_arm.plugin` at 17:09, after the last Stata verification run. Line numbers in those files refer to the tree as it was reviewed and may have shifted.

## Executive summary

The main risks fall into five groups:

- **Silent corruption of user data in data-management commands:**
  - `csort, algorithm(timsort)` duplicates and deletes observations at ≥ 100k rows (SORT-1).
  - `cimport delimited, bindquotes(strict)` loses rows in files ≥ 2 MiB (IMP-1).
  - Large CSV exports can contain NUL-filled holes on Linux (EXP-1).
  - `creghdfe` writes residuals into user variables when `absorb()` uses wildcards or `vce(cluster a b)` is given (REG-1, REG-2).
  - `cdestring` destroys non-numeric strings and mis-scales values (DSTR-1/2/4).
  - `cmerge` drops formats and variable labels, overwrites master value labels, and trusts `sorted` (MRG-2/3/4).
  - XLSX/XLS import silently empties or corrupts common workbooks (XIO-1/2/3).
- **Silently wrong statistics:**
  - `civreghdfe` reports interaction coefficients under the wrong labels (IV-1).
  - `cqreg` tail-quantile standard errors are off by ~10⁶ with rc 0 (QREG-1).
  - `cbinscatter` builds histogram bins instead of quantile bins at ≥ 50k observations per group (BIN-1).
  - `crangestat` returns garbage count/sum/mean/SD for large groups in unbalanced panels (RNG-1).
  - `csample`/`cbsample` can draw only ~9 distinct samples (SMP-1).
  - Many IV diagnostics ignore the chosen VCE or weights (IV-7…13).
  - The creghdfe F-statistic, singleton rules and df are wrong in common cases (REG-4/5/6/9).
  - The binsreg fit lines are wrong (BIN-2/3/4).
- **Crashes and memory safety:**
  - On Windows, aligned buffers are released with `free()` in creghdfe and cpplmhdfe on every call (REG-3).
  - The default numeric sort crashes on keys holding very large values (SORT-2).
  - IPS4o overflows the stack on long shared string prefixes, which crashes Stata inside `cmerge` (SORT-5).
  - `orthog()`/`redundant()` overflow the heap (IV-2).
  - 32-bit offsets overflow on large N×K (REG-7, BIN-6).
  - `threads(#)` has no bound (CORE-1).
- **"Drop-in replacement" divergences:**
  - `cmerge update` never produces `_merge` codes 4/5 (MRG-1).
  - `cimport` splits `;`-delimited files with decimal commas on commas, and mishandles CRLF blank lines, explicit varlists and `rowrange()` (IMP-2…5).
  - `cexport excel` loses numeric precision (EXP-2).
  - `cpsmatch` includes observations with missing outcome or propensity score (PSM-1).
  - `cwinsor` is quadratic on ties and destroys extended missing values (WIN-1/2).
  - `csort` leaves a stale sort marker (SORT-7).
  - `crangestat` fails on any dataset containing a variable name of 24 or more characters, and its `count`/`sum` conventions differ from rangestat (RNG-9, RNG-2).
- **A release gate that cannot catch these:**
  - CI is currently red and cannot publish (BLD-1/2).
  - Several validation helpers pass real divergences (BLD-3/5).
  - Only the runner's own binary, and not the I/O suites, is tested (BLD-4).

**Recurring root causes**, which are worth fixing systemically rather than one item at a time:

1. **Positional ado↔plugin protocols.** Token counts and fixed buffers for varlists, indices and metadata cause REG-1, REG-2, IV-2/3, MRG-9, CORE-2, IMP-11, EXP-3 and DSTR-8. Pass explicit per-variable metadata, and have C cross-check the counts against `SF_nvars()`.
2. **Silent numeric fallbacks** replace errors: QREG-1, IV-14, IV-22, BIN-9, EXP-8, PPML-2. Adopt one rule: a failed computation returns an error, never a substitute value.
3. **Shared parsing and formatting helpers** that are not exact: DSTR-3/5, IMP-8, EXP-2. Keep one correctly rounded parser and one shortest-round-trip formatter.
4. **32-bit size and offset arithmetic:** REG-7, BIN-6, IV-23 and SORT-2. Use `size_t` everywhere and turn on `-Wshorten-64-to-32` in CI.
5. **Parallel partitioning** that assumes a clean boundary or a fixed team size: SORT-1/6, IMP-1, EXP-1.
6. **Validation blind spots:**
   - Comparisons skipped when the reference names a statistic differently.
   - 0.5-significant-figure tolerances.
   - Positional coefficient comparison.
   - No large-N, `threads(1)` or Windows runs.

   Together these let most of the P0/P1 items above pass the current 2,877-check suite (BLD-3/4/5).

**Suggested order:**
1. The P0 data-loss items (SORT-1, IMP-1, EXP-1, REG-1, REG-2, SMP-1).
2. The remaining P0 statistical and crash items (IV-1, QREG-1, BIN-1, RNG-1, REG-3).
3. The P1 data-integrity items (MRG-*, DSTR-*, XIO-*, SORT-2/5/7).
4. The P1 statistical items.
5. The gate (BLD-*). Turn each reproduction in this report into a regression test that fails before the fix.

## Summary counts

**172 findings:** 11 P0, 57 P1, 79 P2, 25 P3. Evidence: 83 fully and 11 partly reproduced in Stata, 53 confirmed with native harnesses, 25 established by code reading (5 of them marked plausible).

| Component | P0 | P1 | P2 | P3 | Total |
|---|---:|---:|---:|---:|---:|
| Core infrastructure and shared data layer |  |  | 3 | 4 | 7 |
| Sorting (csort + shared sort engines) | 1 | 5 | 4 | 3 | 13 |
| cmerge |  | 5 | 9 | 2 | 16 |
| cdestring / cencode / cdecode |  | 4 | 5 | 1 | 10 |
| cimport (delimited text) | 1 | 6 | 7 | 1 | 15 |
| cexport (delimited + Excel) | 1 | 1 | 8 | 1 | 11 |
| creghdfe + shared estimation kernels | 3 | 4 | 3 | 2 | 12 |
| civreghdfe | 1 | 12 | 11 | 1 | 25 |
| cqreg | 1 | 2 | 4 | 1 | 8 |
| cpplmhdfe |  |  | 4 | 1 | 5 |
| cbinscatter | 1 | 5 | 5 | 1 | 12 |
| cimport/cexport Excel and statistical formats |  | 6 | 5 | 1 | 12 |
| csample / cbsample / cwinsor / cpsmatch | 1 | 3 | 3 | 3 | 10 |
| crangestat / crangejoin / cipolate / csplit | 1 | 2 | 6 | 1 | 10 |
| Build, packaging, CI and validation |  | 2 | 2 | 2 | 6 |
| **All** | **11** | **57** | **79** | **25** | **172** |

## Prioritized index

Each ID links to its detailed entry below. Evidence: *Stata-verified* = reproduced in Stata 18.5 MP against the current `build/` plugin (*part* = the headline behavior was reproduced; related sub-items rest on harness or code evidence); *Harness* = native C harness compiling the real sources with a mocked Stata interface (usually under ASan/UBSan); *Code* = source inspection (with the native command's ado source where relevant).

### P0 — fix immediately (11)

| ID | Component | Finding | Evidence |
|---|---|---|---|
| [SORT-1](#sort-1) | Sorting (csort + shared sort engines) | `csort, algorithm(timsort)` silently duplicates and deletes observations (N ≥ 100,000) | Stata-verified |
| [IMP-1](#imp-1) | cimport (delimited text) | `cimport delimited, bindquotes(strict)` (and `nobind`) corrupts, drops, or duplicates rows in files ≥ 2 MiB | Stata-verified |
| [EXP-1](#exp-1) | cexport (delimited + Excel) | Parallel CSV writer issues > 2 GiB `pwritev` batches and ignores short writes: Linux silently publishes CSVs with NUL-filled holes; macOS large exports fail with r(693) | Harness |
| [REG-1](#reg-1) | creghdfe + shared estimation kernels | `creghdfe ..., absorb(fe*)` or `absorb(fe1-fe3)`: only the first FE is absorbed and residuals/sample flags are written into user variables | Stata-verified |
| [REG-2](#reg-2) | creghdfe + shared estimation kernels | `vce(cluster a b)` (advertised as `cluster varlist`) silently clusters on `a` only; with weights `b` is used as the weight; with `resid()`/`groupvar()`/`savefe` the outp… | Stata-verified |
| [REG-3](#reg-3) | creghdfe + shared estimation kernels | Windows: `obs_map` allocated with `_aligned_malloc` is released with plain `free()` in creghdfe and cpplmhdfe (heap corruption on essentially every call) | Code |
| [IV-1](#iv-1) | civreghdfe | Factor-variable base/omitted re-insertion misaligns `e(b)`/`e(V)`: interaction coefficients are silently reported under the wrong labels (and `i.a##i.b` fails) | Stata-verified |
| [QREG-1](#qreg-1) | cqreg | cqreg silently uses sparsity = 1 (or f_i = 1, IQR/0.5) when τ±h falls outside (0,1) or an auxiliary quantile solve fails: standard errors off by up to ~10⁶ with rc 0 (qr… | Stata-verified |
| [BIN-1](#bin-1) | cbinscatter | With ≥ 50,000 observations in a group, cbinscatter bins by a 4,096-bucket equal-width histogram instead of quantiles: skewed x collapses into a few bins | Stata-verified |
| [SMP-1](#smp-1) | csample / cbsample / cwinsor / cpsmatch | `csample`/`cbsample` seeds are truncated to one decimal digit: only ~9 distinct samples/bootstrap draws exist per dataset | Stata-verified |
| [RNG-1](#rng-1) | crangestat / crangejoin / cipolate / csplit | crangestat builds prefix sums inside a nested OpenMP region that assumes a full thread team: count/sum/mean/sd/variance are garbage for every by-group with ≥ 10,000 obse… | Stata-verified |

### P1 — fix before the next release (57)

| ID | Component | Finding | Evidence |
|---|---|---|---|
| [SORT-2](#sort-2) | Sorting (csort + shared sort engines) | Default numeric sort (counting sort) crashes Stata or silently loses rows when a key holds values beyond ±2^63 (e.g. ±1e30 sentinels) | Stata-verified |
| [SORT-3](#sort-3) | Sorting (csort + shared sort engines) | IPS4o fails whenever only one thread is available (N ≥ 30,000): `cmerge ..., threads(1)` and every cmerge on a 1-CPU machine fail | Stata-verified |
| [SORT-4](#sort-4) | Sorting (csort + shared sort engines) | Default string sort (MSD radix) degrades to quadratic time when keys share a ≥32-byte prefix (URLs, paths, padded IDs) | Stata-verified |
| [SORT-5](#sort-5) | Sorting (csort + shared sort engines) | IPS4o and sample-sort string radix recurse once per shared character: stack overflow crashes Stata on long repeated str2045 keys (including every cmerge on such keys) | Stata-verified |
| [SORT-7](#sort-7) | Sorting (csort + shared sort engines) | `csort ..., nosortedby` leaves Stata's previous sort marker in place; `by` groups and csort's own early exit then trust unsorted data | Stata-verified |
| [MRG-1](#mrg-1) | cmerge | `cmerge ..., update` never produces `_merge` codes 4 (missing updated) or 5 (nonmissing conflict), so `keep()`/`assert()` select different observations than `merge` | Stata-verified |
| [MRG-2](#mrg-2) | cmerge | New variables from the using data lose display formats, variable labels, and characteristics (dates appear as raw numbers) | Stata-verified |
| [MRG-3](#mrg-3) | cmerge | Using value-label definitions overwrite same-named master labels, so existing master variables silently display the wrong text | Stata-verified |
| [MRG-4](#mrg-4) | cmerge | `sorted` is trusted without verification: an unsorted master silently mis-joins (native errors r(5)) | Stata-verified |
| [MRG-5](#mrg-5) | cmerge | strL variables: any non-identity merge fails with a misleading r(5); on the "identity" fast path strL keepusing values are written from worker threads (concurrent strL S… | Stata-verified |
| [DSTR-1](#dstr-1) | cdestring / cencode / cdecode | `cdestring` converts (and with `replace` destroys) variables that native `destring` refuses because they contain non-numeric text, even without `force` | Stata-verified |
| [DSTR-2](#dstr-2) | cdestring / cencode / cdecode | `cdestring, percent` divides only the observations that contain "%"; native divides the whole variable by 100 | Stata-verified |
| [DSTR-3](#dstr-3) | cdestring / cencode / cdecode | `cdestring` writes invalid values into the dataset: ±infinity for "1e400"/"-1e400", and raw doubles above Stata's maximum (which display as `.z_`/garbage missing codes) | Stata-verified |
| [DSTR-4](#dstr-4) | cdestring / cencode / cdecode | `cdestring, dpcomma` treats "." as a thousands separator anywhere: "3.5" → 35, "-.5" → −5, "1.25" → 125, with no warning | Stata-verified |
| [IMP-2](#imp-2) | cimport (delimited text) | Delimiter auto-detection ignores the header line and breaks ties in favor of the comma: `;`-delimited files with decimal commas (e.g., `abc;2,5;3,0`) are split on commas | Stata-verified |
| [IMP-3](#imp-3) | cimport (delimited text) | Blank lines in CRLF (Windows) files become all-missing observations | Stata-verified |
| [IMP-4](#imp-4) | cimport (delimited text) | An explicit variable list (`cimport delimited a b using f`) forces `varnames(nonames)`: the header row is imported as data and every column becomes a string | Stata-verified |
| [IMP-5](#imp-5) | cimport (delimited text) | `rowrange()` starting at line 1 is mistranslated for headerless files | Stata-verified |
| [IMP-6](#imp-6) | cimport (delimited text) | `varnames(#≥2)` with `bindquotes(strict)` loses rows and takes the header from the wrong line (SIMD strict row finder skips row boundaries) | Stata-verified |
| [IMP-7](#imp-7) | cimport (delimited text) | Windows: files are opened with `CreateFileA` on a UTF-8 path (non-ASCII paths fail) and with `FILE_SHARE_READ` only (files open in Excel fail) | Code |
| [EXP-2](#exp-2) | cexport (delimited + Excel) | `cexport excel` writes numbers with at most 15 (often 11–14) significant digits: 16-digit IDs collapse and fractions lose precision | Stata-verified |
| [REG-4](#reg-4) | creghdfe + shared estimation kernels | Robust/cluster F statistic uses the first `K_keep` entries of the full `b`/`V`: wrong (down to F = 0) whenever an omitted regressor precedes a kept one (SEP24 #14 fixed… | Stata-verified |
| [REG-5](#reg-5) | creghdfe + shared estimation kernels | Singleton removal ignores fweights: single-row FE levels with fweight ≥ 2 are dropped (reghdfe drops a level only when its total fweight is 1) | Stata-verified |
| [REG-6](#reg-6) | creghdfe + shared estimation kernels | Singleton removal silently stops after 100 passes: long singleton chains survive (creghdfe and civreghdfe) | Stata-verified |
| [REG-7](#reg-7) | creghdfe + shared estimation kernels | 32-bit column offsets (`k*N`, `j*N`) overflow once #regressors × N ≥ 2^31 in the shared solver, VCE and PPML/cqreg code (memory corruption on large data) | Code |
| [IV-2](#iv-2) | civreghdfe | `orthog()`/`redundant()` write past the `Z_rest` heap buffer when a tested instrument is repeated or follows a dropped collinear instrument | Harness |
| [IV-3](#iv-3) | civreghdfe | `orthog()`, `endogtest()`, `redundant()` and `partial()` indices are positions in the raw (unexpanded) token list: the wrong column is tested/partialled with factor or t… | Harness |
| [IV-4](#iv-4) | civreghdfe | Regressors or instruments fully absorbed by the fixed effects are not omitted: they become unit-variance noise columns with absurd coefficients, and inflate K, J df and… | Stata-verified |
| [IV-5](#iv-5) | civreghdfe | `partial()` with aweights/fweights/pweights uses an unweighted FWL projection, changing the coefficients | Stata-verified |
| [IV-6](#iv-6) | civreghdfe | Collinear `partial()` variables make the P'P Cholesky fail; the loop breaks and **all** partial controls are silently dropped (omitted-variable bias, rc 0) | Stata-verified |
| [IV-7](#iv-7) | civreghdfe | Hansen J under `vce(robust)` with aweights/pweights uses `w·e²` instead of `(w·e)²`: every overidentified pweight model reports a wrong J | Stata-verified |
| [IV-8](#iv-8) | civreghdfe | C-statistic (`orthog`), endogeneity test (`endogtest`) and redundancy test are always homoskedastic, whatever the VCE | Stata-verified |
| [IV-9](#iv-9) | civreghdfe | With more than one endogenous regressor, `e(idstat)` is always the homoskedastic Anderson LM (labelled Kleibergen-Paap), and the KP Wald F ignores HAC/DK kernels and use… | Stata-verified |
| [IV-10](#iv-10) | civreghdfe | Driscoll–Kraay time index follows order of first appearance, not calendar order: DK SEs and DK test statistics are wrong in unbalanced panels | Harness |
| [IV-11](#iv-11) | civreghdfe | `bw()`/`kernel()` without `robust` give plain iid SEs or silently switch to HAC; ivreg2 computes the AC estimator | Stata-verified |
| [IV-12](#iv-12) | civreghdfe | HAC panel/time structure is not passed to the plugin: panel HAC uses the first `absorb()` variable as the panel id; lags are row distances (gaps and missing rows mis-pai… | Stata-verified (part) |
| [IV-13](#iv-13) | civreghdfe | `gmm2s`/`cue` ignore the second cluster dimension and HAC/DK/Kiefer structure in both the weighting matrix and the VCE; CUE's robust objective uses `w·e²` for aw/pw | Harness |
| [QREG-2](#qreg-2) | cqreg | `denmethod(residual)` does not implement qreg's residual method: it drops every \|r\| < 1e-8 instead of the K basis observations and clamps s to 1e-10 (SEs ~1e-12 on discr… | Stata-verified |
| [QREG-3](#qreg-3) | cqreg | `denmethod(kernel)` silently runs the residual method (the kernel code is dead) while `e(denmethod)` reports "kernel" | Stata-verified |
| [BIN-2](#bin-2) | cbinscatter | Classic method with `controls()`/`absorb()`: the fit line and `e(coefs)` are estimated on mean-zero residuals but drawn against mean-restored dots (line ~130 units below… | Stata-verified |
| [BIN-3](#bin-3) | cbinscatter | `method(binsreg)`: the fit line ignores controls and absorbed effects (raw y on raw x) | Stata-verified |
| [BIN-4](#bin-4) | cbinscatter | `method(binsreg)` writes covariate-adjusted means to the wrong bins whenever a quantile bin is empty (ties/mass points); with `absorb()` the data are re-binned different… | Harness |
| [BIN-5](#bin-5) | cbinscatter | `absorb()` values are used directly as array indices: negative codes silently dropped, fractional codes merged, codes ≥ 2^31 undefined/dropped (r(2001)), large codes blo… | Stata-verified |
| [BIN-6](#bin-6) | cbinscatter | 32-bit overflow in quantile arithmetic for large N × nquantiles (e.g., nq = 100 and N ≥ 21.5M): bin 1 receives the top 3.4% of the data | Harness |
| [XIO-1](#xio-1) | cimport/cexport Excel and statistical formats | XLSX sheets without a usable `<dimension>` element (none, `ref="A1"`, single-cell `ref="B3"`), namespace-prefixed worksheets (`<x:sheetData>`), cells without `r=`, or si… | Stata-verified |
| [XIO-2](#xio-2) | cimport/cexport Excel and statistical formats | XLSX shared strings are corrupted: truncated at 4,095 raw bytes (before entity decoding, possibly mid-UTF-8), whitespace-only text dropped ("Hello World" → "HelloWorld",… | Stata-verified |
| [XIO-3](#xio-3) | cimport/cexport Excel and statistical formats | Legacy `.xls` import corrupts formula and boolean/error cells: text formulas → 0, `#N/A` → 42, `""` formulas → 0, TRUE/FALSE → the string "bool", error cells → "error" | Stata-verified |
| [XIO-4](#xio-4) | cimport/cexport Excel and statistical formats | `cexport dbase` fails with r(108) on ordinary values (0.05, −0.1, 0.001, 1.5e-5, float −0.3, bytes ≤ −100) | Stata-verified |
| [XIO-5](#xio-5) | cimport/cexport Excel and statistical formats | Worksheet parts are chosen by tab position (`sheet{N}.xml`) instead of workbook relationships: with a chart sheet or reordered parts the wrong sheet is imported silently… | Stata-verified |
| [XIO-6](#xio-6) | cimport/cexport Excel and statistical formats | Excel columns mixing date, datetime and General cells import with mixed units (days and milliseconds under one format; General cells left as raw serials); the XLS reader… | Stata-verified |
| [PSM-1](#psm-1) | csample / cbsample / cwinsor / cpsmatch | `cpsmatch` matches observations psmatch2 excludes (missing outcome, missing covariates/pscore, outside `e(sample)`): wrong ATT | Stata-verified |
| [WIN-1](#win-1) | csample / cbsample / cwinsor / cpsmatch | `cwinsor ..., trim` converts extended missing values (.a–.z) to `.` on Apple Silicon (SIMD path), destroying them in the default in-place mode | Stata-verified |
| [WIN-2](#win-2) | csample / cbsample / cwinsor / cpsmatch | cwinsor's percentile selection is quadratic on tied data (zero-inflated or dummy variables): 235× slower than winsor2 at 1M rows, effectively hanging at 10M | Stata-verified |
| [RNG-2](#rng-2) | crangestat / crangejoin / cipolate / csplit | crangestat `count`/`sum` for empty or all-missing windows differ from rangestat, and `count` depends on which path runs (group size < 64 vs ≥ 64) | Stata-verified |
| [RNG-9](#rng-9) | crangestat / crangejoin / cipolate / csplit | crangestat fails with r(198) whenever *any* variable in the dataset has a name of 24 or more characters, because it creates a `__varpos_<varname>` local for every variab… | Stata-verified |
| [BLD-1](#bld-1) | Build, packaging, CI and validation | CI cannot stage or publish a package: the new codec licence file (`ctools-codecs-LICENSE.txt`) is omitted from the CI copy step and the release-gate test (and is untrack… | Harness |
| [BLD-2](#bld-2) | Build, packaging, CI and validation | Stale native regression tests no longer compile (`test_p2_native.py`, `test_cimport_options_native.py`); one of them breaks the Linux CI gate and the SEP20 A08/A22/A29 r… | Harness |

### P2 — next correctness/robustness cycle (79)

| ID | Component | Finding | Evidence |
|---|---|---|---|
| [CORE-1](#core-1) | Core infrastructure and shared data layer | Unbounded `threads(#)` can terminate Stata: the OpenMP runtime aborts the process when thread creation fails | Harness |
| [CORE-2](#core-2) | Core infrastructure and shared data layer | The `__ctools_strw` width hint is read into a fixed 16 KiB buffer; a width cut mid-number makes valid loads fail on very wide datasets (SEP24 #13 analogue) | Harness |
| [CORE-3](#core-3) | Core infrastructure and shared data layer | `_ctools_load` reads the plugin identity but never enforces ado/plugin compatibility; a stale or pinned older plugin silently runs against the new argument protocol | Code |
| [SORT-6](#sort-6) | Sorting (csort + shared sort engines) | Timsort's parallel merge sizes its work split from a team measured in a different parallel region (SEP24 #9 defect class, missed for timsort) | Harness |
| [SORT-8](#sort-8) | Sorting (csort + shared sort engines) | `csort x x, stream(#)` (repeated key) overflows the heap and leaves the key sorted but all other columns unsorted | Harness |
| [SORT-9](#sort-9) | Sorting (csort + shared sort engines) | Streaming csort writes sorted keys before the payload permutation can fail; a later failure leaves every row scrambled (and returns rc 0 in non-OpenMP builds) | Harness |
| [SORT-10](#sort-10) | Sorting (csort + shared sort engines) | csort rejects any dataset containing a strL variable (even a non-key payload), after loading and sorting everything, with a misleading r(920) | Stata-verified |
| [MRG-6](#mrg-6) | cmerge | Repeating a key variable (`cmerge 1:1 id id using ...`) scrambles the master key and silently mis-joins | Stata-verified |
| [MRG-7](#mrg-7) | cmerge | `replace` without `update` is accepted and overwrites master values (native rejects it) | Stata-verified |
| [MRG-8](#mrg-8) | cmerge | Observation order after the merge differs from `merge`, and the sort marker is not managed | Stata-verified |
| [MRG-9](#mrg-9) | cmerge | Fixed plugin argument buffers (8 KB / 4 KB) truncate silently on wide datasets; shared variables are then treated as new and master values are overwritten | Harness |
| [MRG-10](#mrg-10) | cmerge | The hard-coded scratch variable `_merge_temp_filter` can receive `_merge` codes, filter on user data, and then be dropped | Code |
| [MRG-11](#mrg-11) | cmerge | `generate`/`_merge` existence check matches abbreviations: a master variable such as `_merge_prev` causes a spurious r(110) | Code |
| [MRG-12](#mrg-12) | cmerge | Variable notes are copied incorrectly: new notes are invisible, and a shared variable's first master note is overwritten; using variables with names ≥ 27 characters buil… | Code (plausible) |
| [MRG-13](#mrg-13) | cmerge | Valid native syntax is rejected: `keepusing()` wildcards/abbreviations and keep/assert words or abbreviations (`mat`, `mas`, `3 4`, `match_update`) | Stata-verified (part) |
| [MRG-14](#mrg-14) | cmerge | Undocumented hard limits: at most 1,024 non-key using variables and 32 keys (r(198)) | Code |
| [DSTR-5](#dstr-5) | cdestring / cencode / cdecode | The separator-aware number parser (cdestring `dpcomma`, cimport `decimalseparator()`/`groupseparator()`) mis-scales values with ≥18 significant digits and is not correct… | Harness |
| [DSTR-6](#dstr-6) | cdestring / cencode / cdecode | `cdestring` turns ".a"–".z" into system missing and reports them as non-numeric; "NA"/"NaN" (any case) are silently accepted as missing | Stata-verified |
| [DSTR-7](#dstr-7) | cdestring / cencode / cdecode | `cdestring` and `cencode` drop variable labels and characteristics; `replace` moves the variable to the end of the dataset | Stata-verified |
| [DSTR-8](#dstr-8) | cdestring / cencode / cdecode | cencode's label-file path is split at whitespace in C and `run` unquoted: cencode fails (or writes a stray file) when Stata's temp directory contains a space (SEP24 #6 a… | Harness |
| [DSTR-9](#dstr-9) | cdestring / cencode / cdecode | `cdestring` has no rollback; a plugin failure leaves new `generate()` variables behind (e.g., strL values over 2,045 bytes, which native handles) | Code |
| [IMP-8](#imp-8) | cimport (delimited text) | `NaN`/`Infinity` spellings make the whole column a string (native: numeric missing); overflow values such as `1e999` are stored as raw ±inf | Harness |
| [IMP-9](#imp-9) | cimport (delimited text) | Leading blank line or header-only file: header handling diverges from native | Stata-verified |
| [IMP-10](#imp-10) | cimport (delimited text) | CR-only (classic Mac) line endings are not recognized | Harness |
| [IMP-11](#imp-11) | cimport (delimited text) | Fixed-size metadata buffers: `numericcols()`/`stringcols()` silently truncated at 4 KB; wide files (≥ ~2,000–5,000 columns) fail r(198); header lines longer than one chu… | Harness |
| [IMP-12](#imp-12) | cimport (delimited text) | A failed scan leaks the `CIMPORT_NUMCOLS`/`CIMPORT_STRCOLS` globals into later imports | Code |
| [IMP-13](#imp-13) | cimport (delimited text) | Encoding conversion allocates 4× the file size, the zero-copy ASCII path is too narrow (numeric-only files are detected as ISO-8859-9), and on Windows every conversion o… | Harness |
| [IMP-14](#imp-14) | cimport (delimited text) | The automatic header decision inspects only the first 100 data rows; native decides from full-file column types | Code (plausible) |
| [EXP-3](#exp-3) | cexport (delimited + Excel) | `delimiter(" ")` is lost and the next option is swallowed: output is comma-separated and `replace`/`novarnames`/`quote` are ignored (SEP24 #6 remnant) | Stata-verified |
| [EXP-4](#exp-4) | cexport (delimited + Excel) | Any value-labeled variable whose name has ≥ 24 characters makes `cexport delimited` and `cexport excel` fail with r(198) | Stata-verified |
| [EXP-5](#exp-5) | cexport (delimited + Excel) | Empty XLSX export with two threads writes a corrupt worksheet (bad CRC); empty exports without a header get dimension `A1:A0` | Stata-verified |
| [EXP-6](#exp-6) | cexport (delimited + Excel) | Temp-file publication changes file identity: new outputs are mode 0600, `replace` resets permissions/ACLs, replaces symlinks instead of writing through, and overwrites r… | Harness |
| [EXP-7](#exp-7) | cexport (delimited + Excel) | XLSX strings with invalid UTF-8 produce malformed XML (unreadable workbook); control characters are silently deleted; leading/trailing blanks, literal `_xHHHH_`, and CR… | Harness |
| [EXP-8](#exp-8) | cexport (delimited + Excel) | XLSX chunk formatting ignores allocation failures: rows silently missing from a "successful" export | Harness |
| [EXP-9](#exp-9) | cexport (delimited + Excel) | Memory blow-ups: decoded label columns are recast to `str2045` (N × 2,045 bytes each), and when any column is strL every `str#` cell gets its own 2,046-byte malloc | Harness |
| [EXP-10](#exp-10) | cexport (delimited + Excel) | Windows: file APIs receive Stata's UTF-8 paths in the ANSI code page; for XLSX the reserved temp and the miniz-written workbook can be different files (empty file publis… | Code (plausible) |
| [REG-8](#reg-8) | creghdfe + shared estimation kernels | When every regressor is omitted (`K_keep == 0`) creghdfe fails with undefined-macro errors (r(111) under clustering) and the C path returns the wrong RSS | Stata-verified |
| [REG-9](#reg-9) | creghdfe + shared estimation kernels | Absorbed degrees of freedom (`e(df_a)`) diverge from reghdfe's `estimate_dof()` for G ≥ 3 FEs, cluster-nested FEs and the `dof()` options | Stata-verified |
| [REG-10](#reg-10) | creghdfe + shared estimation kernels | Shared `detect_collinearity()` uses an absolute 1e-14 pivot tolerance: legitimate small-scale regressors are dropped (creghdfe `nostandardize`, cqreg, civreghdfe) | Harness |
| [IV-14](#iv-14) | civreghdfe | gmm2s/cue with a singular moment covariance (few clusters): gmm2s silently returns 2SLS labelled gmm2s; cue posts V = 0 or ~1e-10 SEs | Stata-verified (part) |
| [IV-15](#iv-15) | civreghdfe | Jacobi eigen and `M^{-1/2}` routines stop after 100 single rotations: inaccurate for K ≥ 10 (LIML λ with ≥ 9 endogenous regressors; KP with ≥ 10 excluded instruments) | Harness |
| [IV-16](#iv-16) | civreghdfe | `partial()` bookkeeping: phantom trailing `e(b)`/`e(V)` columns; df_r, V, df_m and rmse ignore the partialled count; `nopartialsmall` works in the wrong direction | Harness |
| [IV-17](#iv-17) | civreghdfe | Display/option handling changes stored results: `nofooter` suppresses posting of every diagnostic `e()`; `rf` leaves creghdfe's reduced-form results in `e()` (and builds… | Stata-verified |
| [IV-18](#iv-18) | civreghdfe | HAC test statistics are O(N²): minutes to hours for time series or badly grouped panels, with no way to interrupt | Harness |
| [IV-19](#iv-19) | civreghdfe | J uses the Bartlett kernel regardless of `kernel()` and the HC S under two-way clustering; `vce(cluster tvar) bw()` and `vce(cluster a b) bw()` silently ignore the kernel | Harness |
| [IV-20](#iv-20) | civreghdfe | Documented DWH test `e(endog_chi2)`/`e(endog_p)` uses `(v'v)⁻¹` instead of `[(X_a'X_a)⁻¹]_vv` and ignores weights/VCE (≈10× inflated) | Harness |
| [IV-21](#iv-21) | civreghdfe | `dofminus()`/`sdofminus()` never reach V or the tests; `small` is a no-op; without `absorb()` small-sample scaling is always applied (diverges from ivreg2's large-sample… | Stata-verified |
| [IV-22](#iv-22) | civreghdfe | IV VCE and test failures are still silent (SEP24 #11 only partly propagated): two-way and Kiefer VCE are `void`; HAC buffer failure returns rc 0 with V = 0; thread-local… | Harness |
| [IV-23](#iv-23) | civreghdfe | Integer overflow and large-data hazards: `int` N·K allocation sizes can wrap (> ~537M obs with K = 8); fweights summing to > 2^31−1 overflow `(ST_int)N_eff` (negative V)… | Code |
| [IV-24](#iv-24) | civreghdfe | Two-way clustering: an FE nested in the *second* cluster variable is not treated as redundant (canonical `absorb(firm year) vce(cluster firm year)` design) | Code |
| [QREG-4](#qreg-4) | cqreg | The Frisch–Newton corrector omits the /x and /s scaling of the second-order Mehrotra term: 2–5× more iterations (r(430) at the default limit on tail quantiles) and a tri… | Harness |
| [QREG-5](#qreg-5) | cqreg | `vce(cluster a b)` silently clusters on the a×b intersection | Stata-verified |
| [QREG-6](#qreg-6) | cqreg | If every regressor is dropped (`K_keep = 0`), `e(V)` is all zeros and SE(_cons) = 0 | Stata-verified |
| [QREG-7](#qreg-7) | cqreg | Absolute tolerances break scale-equivariance: auxiliary gap tolerance 1e-6 (SEs 22× too small for y ~ 1e-8), robust regularized inverse gives ~1e10 SEs, tiny-scale regre… | Harness |
| [PPML-1](#ppml-1) | cpplmhdfe | Default `use_exact_partial(1)` (ppmlhdfe defaults to 0) plus an early-stopping inner CG makes IRLS oscillate: r(430) or ~26× more iterations on poorly connected two-way… | Harness |
| [PPML-2](#ppml-2) | cpplmhdfe | False IRLS convergence for large-scale outcomes (0.1 floor on the standardized deviance): wrong `e(ll)`, `e(deviance)`, `e(r2_p)` and slightly wrong V with `e(converged)… | Harness |
| [PPML-3](#ppml-3) | cpplmhdfe | `savefe`/`absorb(name=var)` uses unaccelerated sweeps capped at `iterate()`: a successful fit turns into r(430) with no message | Harness |
| [PPML-4](#ppml-4) | cpplmhdfe | Factor variables are expanded before the FE/singleton screen and cluster markout, so base/omitted levels (and `_b[#.var]`) can differ from ppmlhdfe | Code (plausible) |
| [BIN-7](#bin-7) | cbinscatter | By-group bookkeeping: the plugin still derives the group count as max − min + 1 (SEP20 A07 fix incomplete), so emptying the lowest group silently drops the highest one;… | Harness |
| [BIN-8](#bin-8) | cbinscatter | binscatter's automatic discrete mode (#unique(x) ≤ nquantiles) is missing (Likert 1–7 → 6 bins); `discrete` with > 500 distinct values fails with a plugin error | Harness |
| [BIN-9](#bin-9) | cbinscatter | HDFE residualization that never converges is silently accepted (10,000 sweeps; message only under `verbose`; batch solver checks convergence on y only) — the cbinscatter… | Harness |
| [BIN-10](#bin-10) | cbinscatter | Polynomial fits use uncentered raw moments with plain Cholesky: `qfit` with binary x fails the whole command (r(199)); fits on dates/datetimes are visibly wrong; constan… | Harness |
| [BIN-11](#bin-11) | cbinscatter | `method(binsreg)` diverges from binsreg with `absorb()` (singleton FE groups kept → every dot shifted) and with factor-variable controls (evaluated at their means instea… | Code (plausible) |
| [XIO-7](#xio-7) | cimport/cexport Excel and statistical formats | Unparseable numeric `<v>` text (including ISO dates in `t="d"` cells) becomes missing, and in string columns becomes the literal text "8.98846567431158e+307" | Stata-verified |
| [XIO-8](#xio-8) | cimport/cexport Excel and statistical formats | SPSS/SAS-catalogue value labels on non-integer keys silently overwrite the integer key's label (1.5 "one and a half" replaces 1 "one") | Stata-verified |
| [XIO-9](#xio-9) | cimport/cexport Excel and statistical formats | Left-justified date formats (`%-td`, `%-tc`) are not recognized as dates on export: SPSS/SAS XPORT/dBase files contain raw day or millisecond counts | Stata-verified |
| [XIO-10](#xio-10) | cimport/cexport Excel and statistical formats | SAS XPORT naming: `cexport sasxport8` fails with r(610) when the file stem is not a SAS name (`sales-2024`, `2024sales`, `my data`); `cexport sasxport5` writes value-lab… | Stata-verified |
| [XIO-11](#xio-11) | cimport/cexport Excel and statistical formats | Smaller XLSX-import defects | Stata-verified (part) |
| [PSM-2](#psm-2) | csample / cbsample / cwinsor / cpsmatch | cpsmatch tie-breaking and `noreplacement` algorithm differ from psmatch2 (different matched controls and ATT) | Stata-verified |
| [SMP-2](#smp-2) | csample / cbsample / cwinsor / cpsmatch | The shared counting sort (by()/strata()/cluster() grouping) corrupts groups for integer keys beyond ±2^63 (see SORT-2); csample/cwinsor allocate threads × N work buffers… | Harness |
| [WIN-3](#win-3) | csample / cbsample / cwinsor / cpsmatch | cwinsor's default in-place mode stores fractional bounds into integer variables without promotion (bounds 1.5/99.5 become 1/99) | Stata-verified |
| [RNG-3](#rng-3) | crangestat / crangejoin / cipolate / csplit | crangestat `first`/`last` return the first/last *non-missing* value (identical to `firstnm`/`lastnm`), unlike rangestat and crangestat's own help | Stata-verified |
| [RNG-4](#rng-4) | crangestat / crangejoin / cipolate / csplit | crangestat default result names are `stat_var` (e.g., `mean_x`); rangestat creates `var_stat` (`x_mean`) | Stata-verified |
| [RNG-5](#rng-5) | crangestat / crangejoin / cipolate / csplit | Reversed bounds (`low > high`) give plausible-looking garbage (count ≈ 1.8e19, sd ≈ 9.5e153) in groups with ≥ 64 observations instead of missing | Stata-verified |
| [RNG-6](#rng-6) | crangestat / crangejoin / cipolate / csplit | crangestat skewness/kurtosis: wrong small-n cut-offs (missing for n = 2/3 where rangestat reports values) and constant windows produce fabricated skewness ±1 / kurtosis… | Stata-verified |
| [RNG-7](#rng-7) | crangestat / crangejoin / cipolate / csplit | crangestat precision: prefix-path variance/SD keep only ~7 significant digits on long trending series (SEP20 A10 fix incomplete), skewness/kurtosis lose digits with larg… | Harness |
| [RNG-8](#rng-8) | crangestat / crangejoin / cipolate / csplit | crangestat percentiles and `iqr` use linear (type-7) interpolation, not Stata's `_pctile`/`summarize` definition, and the help does not say so | Stata-verified |
| [BLD-3](#bld-3) | Build, packaging, CI and validation | Validation helpers let real divergences pass: `benchmark_rangestat` records a PASS whenever `rangestat` errors (9 percentile/IQR gate tests can only fail if all output i… | Code |
| [BLD-4](#bld-4) | Build, packaging, CI and validation | The "complete" release gate omits the I/O, transport, big-data and optimization suites, validates only the runner's own architecture's plugin, and UBSan findings cannot… | Code |

### P3 — minor, latent, or documentation (25)

| ID | Component | Finding | Evidence |
|---|---|---|---|
| [CORE-4](#core-4) | Core infrastructure and shared data layer | Thread-count policy ignores `OMP_NUM_THREADS`, CPU affinity/cgroup limits and Stata's `set processors`, and overwrites the OpenMP ICV on every call | Code |
| [CORE-5](#core-5) | Core infrastructure and shared data layer | `ctools, update` reports "already up to date" when the host cannot be reached (r(631)) | Code |
| [CORE-6](#core-6) | Core infrastructure and shared data layer | Stale-state cleanup compares against the command from two calls earlier, so caches from an interrupted command can survive indefinitely | Harness |
| [CORE-7](#core-7) | Core infrastructure and shared data layer | Latent shared-layer hazards | Harness |
| [SORT-11](#sort-11) | Sorting (csort + shared sort engines) | Pairs-sort fallback re-runs scatter tasks already executed after a partial batch submission failure (heap overflow / duplicated rows) | Harness |
| [SORT-12](#sort-12) | Sorting (csort + shared sort engines) | Minor sort-engine defects: `threads(1)` not honored by sample/counting sort at 50k ≤ N < 100k; pairs radix sorts only 32-bit keys; timsort cleanup reads uninitialized `i… | Harness |
| [SORT-13](#sort-13) | Sorting (csort + shared sort engines) | csort compatibility and documentation gaps versus native `sort` | Code |
| [MRG-15](#mrg-15) | cmerge | Error codes and SPI status: plugin codes 1–5 surface as unrelated Stata errors (r(1) "Break", r(5) "not sorted"); keepusing and `_merge` writers ignore SPI return codes | Stata-verified (part) |
| [MRG-16](#mrg-16) | cmerge | Presentation and edge-case divergences: `_merge` has no value/variable label; new variables follow `keepusing()` token order; a failed `assert()` restores the master (na… | Stata-verified (part) |
| [DSTR-10](#dstr-10) | cdestring / cencode / cdecode | Other cdestring/cencode/cdecode divergences and help errors | Stata-verified (part) |
| [IMP-15](#imp-15) | cimport (delimited text) | Smaller divergences from `import delimited` | Stata-verified (part) |
| [EXP-11](#exp-11) | cexport (delimited + Excel) | Other export divergences and defects | Stata-verified (part) |
| [REG-11](#reg-11) | creghdfe + shared estimation kernels | Wrapper compatibility and documentation divergences from reghdfe | Code |
| [REG-12](#reg-12) | creghdfe + shared estimation kernels | Error-path and robustness defects in the creghdfe C core | Harness |
| [IV-25](#iv-25) | civreghdfe | Other ivreghdfe-compatibility defects | Code |
| [QREG-8](#qreg-8) | cqreg | Other cqreg defects: VCE numerics (uncentered X'X, +1e-10 regularization → wrong/negative/asymmetric V for offset regressors; VCE failure posts zeros with rc 0); non-ver… | Harness |
| [PPML-5](#ppml-5) | cpplmhdfe | Smaller ppmlhdfe divergences | Harness |
| [BIN-12](#bin-12) | cbinscatter | Options that do nothing or fail late, result-metadata hygiene, scalability and documentation | Code |
| [XIO-12](#xio-12) | cimport/cexport Excel and statistical formats | Smaller statistical-format and `.xls` export defects | Harness |
| [PSM-3](#psm-3) | csample / cbsample / cwinsor / cpsmatch | Smaller cpsmatch/psmatch2 semantic differences | Harness |
| [SMP-3](#smp-3) | csample / cbsample / cwinsor / cpsmatch | Sampling-command divergences from `sample`/`bsample` | Harness |
| [WIN-4](#win-4) | csample / cbsample / cwinsor / cpsmatch | cwinsor option/type divergences from winsor2 | Stata-verified (part) |
| [RNG-10](#rng-10) | crangestat / crangejoin / cipolate / csplit | Other crangestat/csplit/crangejoin divergences, help errors and resource issues | Stata-verified (part) |
| [BLD-5](#bld-5) | Build, packaging, CI and validation | Other validation-framework false-pass paths | Code |
| [BLD-6](#bld-6) | Build, packaging, CI and validation | Help files document spellings the wrappers reject, and other packaging/build hygiene | Stata-verified |

## Detailed findings

### P0 findings

<a id="sort-1"></a>
#### SORT-1 [P0] — `csort, algorithm(timsort)` silently duplicates and deletes observations (N ≥ 100,000)
*Component: Sorting (csort + shared sort engines)*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified. On 250,000 rows, all 8 timsort key sets tested (numeric, string, multi-key; default threads, `threads(3)`, `threads(4)`) returned rc 0 with `isid id` failing and payload mismatches; the 200,000-row repro below fails `isid id` with r(459). Native harness (A03): 99,999 rows lost at N=200k/4 threads, 936,512 lost at N=1M/8 threads.
- **Location:** `src/ctools_sort_timsort.c:627` (numeric co-rank), `:689-690` (string), `:753-754` (multi-key); merges ≥ `PARALLEL_MERGE_THRESHOLD` (100000, `:580`) take this path (`:785`, `:802`, `:819`).
- **Problem:** Each thread locates its start in the two runs with a co-rank binary search whose predicate is `keys[temp2[mid2-1]] > keys[temp1[mid1]]`. The sequential merge is stable (left run wins ties), so the predicate must be `>=`. With a tied key straddling a chunk boundary, a thread starts with `i1` too small and `i2` too large: it re-emits left-run elements (duplicates) and never emits a slice of the right run (lost rows). The permutation is then applied to every variable, so whole observations are duplicated and others vanish.
- **Trigger:**
  ```stata
  clear
  set obs 200000
  gen long id = _n
  gen byte g = mod(_n, 3)
  csort g, algorithm(timsort) threads(4)
  isid id            // r(459): rows duplicated/lost, csort returned rc 0
  ```
- **Impact:** Silent, irreversible corruption of the whole dataset for a documented algorithm on any realistic key with ties. In stream mode the lost rows also leave inverse-permutation entries uninitialized, which are then used as buffer offsets (`csort_stream.c:245,288,516-527`) — out-of-bounds writes. The validation suite only runs timsort at N ≤ 5,000.
- **Proposed fix:** Change the three predicates to `>=` (`compare_string(...) >= 0`, `compare_multikey(...) >= 0`), mirroring the correct co-rank in `src/ctools_order.c:101`; also apply SORT-6. Add regressions at N=200k with 2/4/8 threads, `nosortedby`, and `stream()`, asserting `isid id`.

<a id="imp-1"></a>
#### IMP-1 [P0] — `cimport delimited, bindquotes(strict)` (and `nobind`) corrupts, drops, or duplicates rows in files ≥ 2 MiB
*Component: cimport (delimited text)*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: a 60,000-row file with multi-line quoted text in every third row (written by native `export delimited, quote`) imports as **46,776 observations** with `id` turned into a strL of merged rows (native: 60,000, `id` long), rc 0. Harness (A06a): strict 120,000-row file → 129,056 obs (duplicates); nobind 150,000-row file with two inch-mark quotes → 272,547 obs. Files < 2 MiB (single chunk) are correct.
- **Location:** `src/cimport/cimport_impl.c:520-556` (`cimport_find_safe_boundary`), `:1042-1056` (boundaries, no monotonicity), `:1058-1092` (chunk workers); `src/cimport/cimport_parse.c:235-417`.
- **Problem:** Files ≥ 2 MiB are split at `i*chunk_size`. In strict/nobind mode the boundary search guesses quote parity from a ~10 KB window and restarts the row finder with `in_quotes = false`; when the target is inside a quoted field, parity is inverted and the chunk starts mid-field, so real row separators are treated as quoted and rows merge. Boundaries are not forced to be monotonic, so a scan that runs to EOF makes neighbouring chunks re-parse the same bytes (duplicates). `nobind` needlessly takes the quote-parity path.
- **Impact:** Silent corruption of whole chunks (row counts, merged fields, wrong types) in exactly the mode the documentation recommends for multi-line text.
- **Proposed fix:** Treat `nobind` like loose (any previous `\n` is a valid boundary). For strict, compute the true quote state at each target from the file start (parity of `"` count; one SIMD popcount pass or per-chunk parity + prefix sum), scan forward to the first unquoted `\n`, and enforce `b[i] = max(b[i], b[i-1])`. Simplest safe interim fix: parse strict mode in a single chunk.

<a id="exp-1"></a>
#### EXP-1 [P0] — Parallel CSV writer issues > 2 GiB `pwritev` batches and ignores short writes: Linux silently publishes CSVs with NUL-filled holes; macOS large exports fail with r(693)
*Component: cexport (delimited + Excel)*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Harness running the real `cexport_main` streaming pipeline with a mock SPI (A07 `h9_csv_big_wave.c`): 500,000 rows × 250 doubles on 12 threads → "write error during parallel I/O" r(693); the same with 4 threads succeeds (2.17 GB). Linux short-write behavior (writes capped at 0x7ffff000 bytes, partial count returned) simulated with `capped_pwritev.c`: rc 0 and the published CSV contains 9.2 million NUL bytes (3,076 of 30,001 lines intact). Not reproduced in Stata here (requires multi-GB exports).
- **Location:** `src/cexport/cexport_impl.c:145-164` (50,000-row chunks), `:480-483` (a wave = `max_threads` chunks), `:657-670` (only `batch_written < 0` checked, offset advanced by the full wave); `src/cexport/cexport_io.c:716-729`, `:772-839` (short counts returned as-is), `:885-899` (`ftruncate`/`fdatasync` results ignored).
- **Problem:** Each wave writes `threads × 50,000 × row_bytes` in one call. Linux truncates any write to ~2 GiB and returns the partial count, which is treated as success; macOS returns EINVAL for > INT_MAX. ENOSPC or NFS short writes are also published as success.
- **Impact:** Silent corruption of large CSV exports on many-core Linux servers (e.g., 5M rows × 80 doubles on 32 threads); default large exports fail on macOS; chunk buffers can reach tens of GB before the write.
- **Proposed fix:** Write in a loop until all bytes are written (retry on short counts/EINTR, fail on 0/−1); cap each call at ≤ 1 GiB by splitting iovecs; check `fdatasync`/`ftruncate`; bound chunk buffers by bytes (e.g., 8–32 MB) rather than rows.

<a id="reg-1"></a>
#### REG-1 [P0] — `creghdfe ..., absorb(fe*)` or `absorb(fe1-fe3)`: only the first FE is absorbed and residuals/sample flags are written into user variables
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: `creghdfe y x, absorb(fe*)` → β = 4.509 (reghdfe 1.999), header "Absorbing 1 HDFE group", `count if e(sample)` = 0, and **1,891 values of the user variable `fe2` were overwritten**; `absorb(fe1-fe2) resid(r)` overwrote **all 2,000 values of `fe2` with residuals**, rc 0. (civreghdfe, cbinscatter, cpplmhdfe handle the same wildcard correctly — verified.)
- **Location:** `build/creghdfe.ado:186` (`local nfe : word count `absorb'`), `:325`, `:336`, `:348`, `:363`, `:374`; `src/creghdfe/creghdfe_regress.c:254`, `:381` (C assumes exactly G FE columns).
- **Problem:** `absorb()` is kept as raw text; its *token* count is sent as G and used to compute every output index, but `markout` and `plugin call` expand the varlist. C absorbs only the first FE, reads cluster/weight columns from the wrong positions, and writes outputs into shifted positions (user variables). No `unab`, no `SF_nvars()` cross-check.
- **Impact:** Wrong estimates plus silent, permanent corruption of input variables (the wrapper ends with `restore, not`), violating the SEP24 #2 guarantee.
- **Proposed fix:** `unab absorb : `absorb'` before counting and building the plugin varlist; in C compute the expected varlist length from K, G, cluster, weights and outputs and return 198 on mismatch with `SF_nvars()`.

<a id="reg-2"></a>
#### REG-2 [P0] — `vce(cluster a b)` (advertised as `cluster varlist`) silently clusters on `a` only; with weights `b` is used as the weight; with `resid()`/`groupvar()`/`savefe` the output overwrites `b` or the weight variable
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: `vce(cluster firm year)` → creghdfe `N_clust = 40` (firm only), reghdfe two-way `N_clust = 15`; with `[aw=w]` creghdfe β = 0.96037 vs reghdfe 0.96554 (weights were read from `year`); with `resid(r)` **all 1,000 values of `year` were overwritten**, rc 0.
- **Location:** `build/creghdfe.ado:205-234` (`gettoken vce_type clustervar : vce` keeps all remaining tokens), `:327`, `:336`, `:348`, `:363`; `src/creghdfe/creghdfe_regress.c:254`, `:381`, `:421-427`; `docs/README_creghdfe.md:30`.
- **Impact:** A standard reghdfe idiom returns wrong SEs; with weights, wrong coefficients; with outputs, corrupted user data.
- **Proposed fix:** Until multi-way clustering is implemented (the cpplmhdfe inclusion–exclusion code could be shared), `unab` the cluster list and exit 198 when it has more than one variable; add the `SF_nvars()` guard from REG-1; fix the README.

<a id="reg-3"></a>
#### REG-3 [P0] — Windows: `obs_map` allocated with `_aligned_malloc` is released with plain `free()` in creghdfe and cpplmhdfe (heap corruption on essentially every call)
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Source-confirmed: `obs_map` comes from `ctools_safe_cacheline_alloc2` → `ctools_aligned_alloc` → `_aligned_malloc` on `_WIN32` (`src/ctools_config.h:350-371`, which itself states it "MUST be freed with ctools_aligned_free()"; `src/ctools_data_io.c:1053`, `:1070`, `:1115`), but is freed with `free(obs_map)` at `src/creghdfe/creghdfe_regress.c:1032`, `:1162`, `:1775` and `src/cpplmhdfe/cpplmhdfe_irls.c:746`, `:847`, `:880`, `:904`, `:920`, `:961`, `:1058`, `:1109`, `:1640`, `:2351`. The same pattern exists in the last commit. Not executed on Windows (validation runs only on macOS, where `posix_memalign` memory may be passed to `free()`); csample/cbsample handle ownership correctly.
- **Location:** allocation `src/ctools_data_io.c:1053`, `:1070`, `:1115` (`ctools_safe_cacheline_alloc2` → `ctools_aligned_alloc`, `src/ctools_config.h:350-371`); plain `free(obs_map)` at `src/creghdfe/creghdfe_regress.c:1032`, `:1162`, `:1775` and `src/cpplmhdfe/cpplmhdfe_irls.c:746`, `:847`, `:880`, `:904`, `:920`, `:961`, `:1058`, `:1109`, `:1640`, `:2351`.
- **Impact:** On Windows, freeing a non-heap-block pointer is undefined behavior; with 64-bit heap termination-on-corruption this typically terminates Stata, otherwise corrupts the heap. Every creghdfe/cpplmhdfe call reaches these frees.
- **Proposed fix:** Replace these `free(obs_map)` calls with `ctools_aligned_free` (or let `ctools_filtered_data_free` own the map); add a Windows CI smoke test that actually runs each command in Stata, and a lint rule (grep) for `free(` on known aligned buffers.

<a id="iv-1"></a>
#### IV-1 [P0] — Factor-variable base/omitted re-insertion misaligns `e(b)`/`e(V)`: interaction coefficients are silently reported under the wrong labels (and `i.a##i.b` fails)
*Component: civreghdfe*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified. `ivreghdfe price (mpg=length) i.foreign i.foreign#c.weight, absorb(rep78)` → `0b.foreign#c.weight = 8.84`, `1.foreign#c.weight = 19.20`; civreghdfe → `0`, `8.84` (the base group's slope is shown as the foreign slope; 19.20 is lost). `i.a#i.b`: every civreghdfe cell estimate is shifted one position (ivreghdfe `1 1 = .101, 2 0 = .453, …, 3 1 = −.181`; civreghdfe `1 1` omitted, `2 0 = .101, …`, and −.181 is lost). `i.treat##i.post` (DiD): ivreghdfe runs, civreghdfe r(504).
- **Location:** `build/civreghdfe.ado:264-329` (base/omitted classification), `:1111-1151`, `:1153-1193` (re-insertion into b and V).
- **Problem:** The re-insertion logic classifies interaction terms by pattern (base-level markers) and misidentifies which expanded columns are base/omitted, so the compact coefficient vector from C is scattered into the wrong positions.
- **Impact:** Silently wrong inference for common factor syntax (group-specific slopes, cell effects, DiD interactions); the existing `i.race#i.union` validation compares only model-level statistics and misses it.
- **Proposed fix:** Build the column map from `fvexpand` output plus `_ms_parse_parts` (or `_rmcoll`'s `o.`/`b.` flags) on the *estimation* sample, carry an explicit index vector (expanded position → compact position or omitted) through the plugin, and assemble b/V from that map; add coefficient-level (not model-level) comparisons against ivreghdfe for `i.a#c.x`, `i.a#i.b`, `i.a##i.b`.

<a id="qreg-1"></a>
#### QREG-1 [P0] — cqreg silently uses sparsity = 1 (or f_i = 1, IQR/0.5) when τ±h falls outside (0,1) or an auxiliary quantile solve fails: standard errors off by up to ~10⁶ with rc 0 (qreg exits r(498))
*Component: cqreg*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified. Realistic case (log-normal outcome, N = 200,000, 6 regressors, `quantile(.95)`): qreg `se(x1) = 6,792.06`, `e(sparsity) = 6,387,376.5`; cqreg `se(x1) = 0.00106`, `e(sparsity) = 1`, rc 0. `sysuse auto`, `quantile(.05)`: qreg r(498) "VCE computation failed", cqreg rc 0 with sparsity 1.
- **Location:** `src/cqreg/cqreg_regress.c:260-351`, `:286-293`, `:449-529`, `:478-485`; `src/cqreg/cqreg_sparsity.c:697-741`; validation counts these as PASS (`validation/validate_cqreg.do:586-627`, `:1449-1478`).
- **Problem:** Every numeric failure in sparsity estimation is replaced by a fallback value instead of an error. The auxiliary solves use a hard-coded 100-iteration limit (the SEP20 A17 fix covers only the main solve), and tail quantiles of skewed data routinely need more (made worse by QREG-4), so the fallback triggers in realistic large-N tail regressions.
- **Impact:** Silently wrong (by many orders of magnitude) standard errors, t-statistics and confidence intervals.
- **Proposed fix:** Replace every fallback with an error status propagated to `cqreg_full_regression`; return r(498) with qreg's message for τ±h outside (0,1), auxiliary non-convergence, allocation failure, or s < DBL_EPSILON; give auxiliary solves the user's `iterate()` (and fix QREG-4); make the validation require rc parity with qreg for extreme quantiles.

<a id="bin-1"></a>
#### BIN-1 [P0] — With ≥ 50,000 observations in a group, cbinscatter bins by a 4,096-bucket equal-width histogram instead of quantiles: skewed x collapses into a few bins
*Component: cbinscatter*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: 60,000 obs, `x = exp(2*rnormal())`, `nquantiles(20)` → only **9 non-empty bins**, the first holding **27,865** observations (expected 20 bins × 3,000; `xtile` gives 20). Harness (A13-1): a single outlier collapses everything into one bin; integer x shifts every bin; normal x is off by 1–2%.
- **Location:** `src/cbinscatter/cbinscatter_bins.c:26-27` (threshold), `:279-395` (unweighted), `:432-549` (weighted); used from `src/cbinscatter/cbinscatter_impl.c:752` and the binsreg+absorb path (`:771-831`).
- **Problem:** Above the threshold, observations are bucketed by `(x−min)·4095/range` and each bucket gets the bin of its end cumulative count — not xtile/`_pctile` cutpoints (what binscatter's `fastxtile` uses).
- **Impact:** Silently wrong binned scatter, `savedata()` and `e(bindata)` for the command's main use case (large N), catastrophic for skewed variables (income, firm size, prices) or any outlier.
- **Proposed fix:** Compute exact cutpoints (nth_element/quickselect per cutpoint or a sort) and assign with the xtile rule (`<= cutpoint` → lower bin), weighted via cumulative weights, for all N; if a histogram is kept for speed, use it only to locate candidate buckets and refine exactly inside boundary buckets. Fix the help/README claims.

<a id="smp-1"></a>
#### SMP-1 [P0] — `csample`/`cbsample` seeds are truncated to one decimal digit: only ~9 distinct samples/bootstrap draws exist per dataset
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: across 30 different `set seed` values, `csample 50` on 1,000 rows produced **6 distinct samples** (native `sample`: 30); 50 consecutive `cbsample` replicates produced **5 distinct** bootstrap samples. Harness (A15-1): seeds `2.30584300921e+18`, `2.99999999999e+18` and `2.00000000001e+18` give byte-identical samples.
- **Location:** `build/csample.ado:72`, `:118`; `build/cbsample.ado:85`, `:144`; `src/csample/csample_impl.c:235-245`, `:294`, `:581`; `src/cbsample/cbsample_impl.c:224-234`, `:284`, `:634`, `:661`.
- **Problem:** The wrappers build `floor(runiform()*2147483647) + floor(runiform()*2147483647)*2147483648` (≈1e17–4.6e18). Stata stores that in a macro as `2.30584300921e+18`; C parses it with `strtoull` and accepts the partial parse, so the seed becomes the leading digit (1–9). Same-seed reproducibility still "works", which is why the tests pass.
- **Impact:** Bootstrap/subsampling loops built on cbsample/csample have ≤ 9 support points: bootstrap SEs and confidence intervals are silently wrong; ~20% of consecutive calls return exactly the previous sample.
- **Proposed fix:** Pass two exact 31-bit integers (e.g., `runiformint(0, 2147483646)` twice) and combine in C, or format with `%20.0f`; make `parse_size_option` reject partially consumed tokens; add a test that distinct seeds give distinct samples and that B bootstrap draws yield ≈ B distinct replicates.

<a id="rng-1"></a>
#### RNG-1 [P0] — crangestat builds prefix sums inside a nested OpenMP region that assumes a full thread team: count/sum/mean/sd/variance are garbage for every by-group with ≥ 10,000 observations in unbalanced panels
*Component: crangestat / crangejoin / cipolate / csplit*

- **Status:** fixed on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: the repro below gives 11,005 of 12,000 big-group observations wrong (e.g., t = 996: count 1.845e19 vs 11, sum −49,483 vs 561, mean −2.7e−15 vs 51); with `threads(1)` 0 are wrong. Harness with the real kernel and libomp (A14-1): 41 groups (one with 12,000 rows, the rest ~10), `interval(t -5 5)` → 11,005 of 12,000 big-group observations wrong with 12 threads (count `1.8446744073709552e+19` or 0, sum −49,483 vs true 561), 6,005 wrong with `threads(2)`, 0 wrong with `threads(1)`.
- **Location:** `src/crangestat/crangestat_impl.c:249-340` (`build_prefix_arrays` parallel scan; `double chunk_sum[64]` `:262`; `omp_get_thread_num()` `:273`, `:322`; phase-2 loop `:305`); cross-group dispatch `:1743-1744`; called from inside `#pragma omp parallel for` at `:1797-1831`.
- **Problem:** For a group with ≥ 10,000 rows the scan splits the group into `ctools_get_max_threads()` fixed chunks and relies on thread `tid` to fill chunk `tid`. In the cross-group path (≥ 8 groups, average size < 500) the function runs inside the outer parallel loop, so the nested region has one thread (nested parallelism is off by default): only chunk 0 is computed, phase 2 reads uninitialized stack slots, and the rest of the prefix arrays keep stale heap contents reused across groups. The same assumption breaks top-level calls when the runtime supplies a smaller team (`OMP_THREAD_LIMIT`, `OMP_DYNAMIC`). SEP24 #9 fixed this pattern in the sort engines only.
- **Trigger:**
  ```stata
  clear
  set obs 12400
  gen long   g = cond(_n <= 12000, 0, 1 + floor((_n - 12001)/10))
  gen double t = cond(_n <= 12000, _n, mod(_n - 12001, 10) + 1)
  gen double x = mod((_n - 1)*37, 101)
  crangestat (count) cn=x (sum) cs=x (mean) cm=x (sd) csd=x, interval(t -5 5) by(g)
  rangestat  (count) x (sum) x (mean) x (sd) x, interval(t -5 5) by(g)
  count if cn != x_count | reldif(cs, x_sum) > 1e-9 | reldif(cm, x_mean) > 1e-9   // thousands
  ```
- **Impact:** Silently wrong values for the most common statistics on unbalanced panels (one large firm/county among many small ones); results vary with thread count and heap contents; reads uninitialized memory.
- **Proposed fix:** Partition into a fixed number of logical chunks processed with `#pragma omp for` (not `tid`), zero-initialize chunk totals, or take the sequential branch when `omp_in_parallel()`; add a regression with ≥ 8 small groups plus one ≥ 10,000-row group checked against `threads(1)`.

### P1 findings

<a id="sort-2"></a>
#### SORT-2 [P1] — Default numeric sort (counting sort) crashes Stata or silently loses rows when a key holds values beyond ±2^63 (e.g. ±1e30 sentinels)
*Component: Sorting (csort + shared sort engines)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `csort x` on 60,000 rows with `x` containing −1e30 and 1e30 **crashed Stata** (SIGSEGV in `counting_sort_numeric_parallel.omp_outlined`; crash report `~/Library/Logs/DiagnosticReports/stata-mp-2026-09-27-155502.ips`). One-sided values (1e19…3e19) returned unsorted with `nosortedby`. Native harness (A03/A15): ASan heap-buffer-overflow at `ctools_sort_counting.c:139`; 47,811 duplicated indices at N=60k with rc 0.
- **Location:** `src/ctools_sort_counting.c:110-116`, `:136-152`, `:171` (sequential) and `:319-326`, `:385-401`, `:465` (parallel); selected by default for numeric keys at `src/ctools_types.c:517`; called explicitly by `crangestat_impl.c:1497`, `csample_impl.c:454`, `cbsample_impl.c:493`, and via AUTO by `cwinsor_impl.c:402`.
- **Problem:** `min_val = (int64_t)min_val_d; range = (size_t)((int64_t)max_val_d - min_val + 1)` is undefined for |v| ≥ 2^63 (saturates on ARM64, INT64_MIN on x86-64). The subtraction overflows, `range` wraps to 0/1 and passes the `COUNTING_SORT_MAX_RANGE` guard, and bucket indices then read/write far outside `counts`/`offsets`. When all huge values are on one side, every value lands in bucket 0 and nothing is sorted (rc 0).
- **Trigger:** `set obs 60000`, `gen double x = mod(_n,5)`, `replace x = -1e30 in 1`, `replace x = 1e30 in 2`, `csort x` → Stata crashes (smaller N: heap corruption or duplicated/lost rows). One-sided: `gen double g = 1e19*(1+mod(_n,3))`, `csort g, nosortedby` → not sorted; `cwinsor v, by(g)` then treats each row as its own group.
- **Impact:** Process crash (unsaved work lost), silent row duplication/loss, or wrong by-group results in cwinsor/crangestat/csample/cbsample for data containing large sentinel or hashed-ID values.
- **Proposed fix:** Before any integer conversion, require `min_val_d >= -2^53 && max_val_d <= 2^53` and compute `max_val_d - min_val_d` in double against `COUNTING_SORT_MAX_RANGE`; otherwise return `STATA_ERR_UNSUPPORTED_TYPE` so the dispatcher falls back to LSD radix. Apply to both paths; compute bucket offsets in unsigned arithmetic only after the check.

<a id="sort-3"></a>
#### SORT-3 [P1] — IPS4o fails whenever only one thread is available (N ≥ 30,000): `cmerge ..., threads(1)` and every cmerge on a 1-CPU machine fail
*Component: Sorting (csort + shared sort engines)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `cmerge 1:1 id using ..., threads(1)` with 30,000 unsorted using rows → r(459) "Failed to sort using data"; `csort x, algorithm(ips4o) threads(1)` → r(920). Harness: N=29,999 ok, N=30,000 fails.
- **Location:** `src/ctools_sort_ips4o.c:521-538`, `:560-562` (numeric), `:712-743` (string); parallel path chosen at `:843`, `:895` purely on `nobs >= 30000`; `src/ctools_arena.c:76-79` returns NULL for size 0; callers `src/cmerge/cmerge_impl.c:342-350`, `:749-757`.
- **Problem:** With one thread `num_buckets = 1`, so `ctools_arena_alloc((num_buckets-1)*sizeof(uint64_t))` requests 0 bytes, gets NULL, and returns `STATA_ERR_MEMORY`. When the master-side sort fails, cmerge returns raw C code 1, which Stata reports as r(1) (the Break code).
- **Trigger:** see evidence; also any cmerge with an unsorted side ≥ 30k rows inside a 1-vCPU container/VM (`ctools_threads.c` reports 1 CPU).
- **Impact:** cmerge unusable on single-CPU environments and with `threads(1)`; misleading error codes.
- **Proposed fix:** Route `ctools_get_openmp_threads() <= 1` to the sequential IPS4o path (or skip the splitter allocation when `num_buckets == 1`); map `STATA_ERR_*` to Stata codes on cmerge's master path; add `threads(1)` cases to native and Stata tests (the SEP24 #9 native test mocks 4 threads).

<a id="sort-4"></a>
#### SORT-4 [P1] — Default string sort (MSD radix) degrades to quadratic time when keys share a ≥32-byte prefix (URLs, paths, padded IDs)
*Component: Sorting (csort + shared sort engines)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: 100,000 rows of a 40-byte URL key with 10 distinct suffixes → default `csort s` 7.58 s vs `csort s, algorithm(lsd)` 0.015 s (output correct). Native harness with production flags: 0.07 s at 10k rows, 4.85 s at 80k, 21.1 s at 160k (quadratic).
- **Location:** `src/ctools_sort_radix_msd.c:57` (depth cap 32), `:820-827` (insertion-sort fallback), `:839-843` (single-bucket branch still increments depth); AUTO picks MSD for string keys (`src/ctools_types.c:506-512`); also reached by `cwinsor ..., by(strvar)`.
- **Problem:** `depth` counts every character position, including positions where all strings share the byte. After 32 shared bytes the whole bucket goes to `insertion_sort_string`, which is O(n²) for unsorted input.
- **Trigger:** `set obs 200000`, `gen str40 s = "https://www.example.com/products/item/" + string(runiformint(0,9))`, `csort s` (tens of seconds; `algorithm(lsd)` is instant). Extrapolates to ~14 min at 1M rows.
- **Impact:** Default `csort`/`cwinsor by()` effectively hang on common string keys.
- **Proposed fix:** Treat the single-bucket case as a loop that advances `char_pos` without incrementing depth; replace the depth-limit insertion sort with an O(n log n) stable fallback (e.g., merge sort comparing with `memcmp` from `char_pos`) or an explicit heap work stack.

<a id="sort-5"></a>
#### SORT-5 [P1] — IPS4o and sample-sort string radix recurse once per shared character: stack overflow crashes Stata on long repeated str2045 keys (including every cmerge on such keys)
*Component: Sorting (csort + shared sort engines)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified crash: `cmerge 1:1 k using u` with 40 str2045 keys sharing a 2,040-byte prefix (unsorted using data) **killed Stata** (exit status 139 / SIGSEGV) on two consecutive runs. Native harness with production flags (A03): SIGSEGV for shared prefixes ≥ ~1,980 bytes and for the 2,045 strings `"a"`, `"aa"`, …; MSD/LSD/timsort/merge handle them.
- **Location:** `src/ctools_sort_ips4o.c:252-311` (self-call at `:280` not tail-call-eliminated; ~4.2 KB frames), callers `:423-426`, `:803-809`; `src/ctools_sort_sample.c:299-372` (non-tail recursion `:365-371`).
- **Problem:** Recursion depth equals shared-prefix length; with an 8 MB stack ~1,970 levels overflow. Strings up to 2,045 bytes are loaded (str2045 and short strL).
- **Trigger:** `gen str2045 s = "x"*2044 + cond(mod(_n,2),"y","z")`, `csort s, algorithm(ips4o)` → crash; the same key in `cmerge` (which always uses IPS4o for string keys) → crash.
- **Impact:** Stata process crash, unsaved work lost.
- **Proposed fix:** Loop instead of recursing when all strings fall in one bucket; bound depth with an O(n log n) stable fallback (as MSD should, SORT-4) or use a heap-allocated work stack; move the 257-entry count/offset arrays to heap scratch.

<a id="sort-7"></a>
#### SORT-7 [P1] — `csort ..., nosortedby` leaves Stata's previous sort marker in place; `by` groups and csort's own early exit then trust unsorted data
*Component: Sorting (csort + shared sort engines)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: after `sort id` then `csort x, nosortedby`, `: sortedby` still returns `id`; `by id: assert _n==1` runs without the "not sorted" error; a subsequent `csort id` early-exits and leaves `id[1] = 6` (unsorted).
- **Location:** `build/csort.ado:170-174` (final `sort` skipped under `nosortedby`), `:83-88` (early exit trusts `: sortedby`).
- **Problem:** The plugin permutes rows through `SF_vstore/SF_sstore`, which does not clear Stata's sort marker; only the skipped final `sort` would have reset it.
- **Trigger:** `clear`, `set obs 6`, `gen id = 7 - _n`, `gen x = _n`, `sort id`, `csort x, nosortedby`, `di "`: sortedby'"` → `id`.
- **Impact:** Silent wrong `by:` results, merges that skip required sorting (`cmerge.ado:266-280` trusts the marker), and `csort` becoming a no-op on a false marker. The same stale-marker question applies to every plugin that rewrites existing variables in place (cwinsor/crangestat replace paths, cmerge).
- **Proposed fix:** Whenever the plugin permutes rows and the final `sort` is skipped, invalidate the marker from Stata (e.g., touch a sort variable with a no-op `replace`, or run `sort` on a tempvar and drop it) and assert `"`: sortedby'" == ""`; or always run the final `sort ..., stable` if a marker existed beforehand.

<a id="mrg-1"></a>
#### MRG-1 [P1] — `cmerge ..., update` never produces `_merge` codes 4 (missing updated) or 5 (nonmissing conflict), so `keep()`/`assert()` select different observations than `merge`
*Component: cmerge*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified. Controlled 5-row example: native `merge 1:1 k using u, update` → one 3, two 4s, one 5; `cmerge` → all four matched rows coded 3 (updated values themselves agree). `update keep(match)`: native keeps 1 row, cmerge keeps 4; `update assert(master using match)`: native r(9), cmerge rc 0; `keep(match_update)`, `keep(match_conflict)`, `keep(4 5)` → cmerge r(198). Also found independently by harness (A05-1).
- **Location:** `src/cmerge/cmerge_impl.c:1445-1450` (writes only the 1/2/3 join result), `src/cmerge/cmerge_io.c:36-106` (applies update/replace values but never reports 4/5); `build/cmerge.ado:176-199` (keep/assert word parsing), `:1063-1137`.
- **Problem:** The update/replace writer decides per shared variable whether to overwrite but does not aggregate "a missing master value was filled" (4) or "a nonmissing master value differed from a nonmissing using value" (5) into the row's merge result.
- **Impact:** Any workflow that inspects or filters on update conflicts (`tab _merge`, `keep(match)`, `assert(match)`) silently gets different rows than native `merge`. The validation suite misses this because every update test uses `nogenerate`.
- **Proposed fix:** In the keepusing writer, per matched row and shared variable, compare master and using values (numeric: both nonmissing and unequal → conflict; master missing and using nonmissing → updated; strings: "" as missing) and fold into a per-row code with precedence 5 > 4 > 3 (atomic max or per-thread arrays); write `_merge` after all variables are processed. Accept native keep/assert words, abbreviations and numeric codes (reuse native `merge`'s mapping).

<a id="mrg-2"></a>
#### MRG-2 [P1] — New variables from the using data lose display formats, variable labels, and characteristics (dates appear as raw numbers)
*Component: cmerge*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: using variable `bdate` with `%td` and label "Birth date" → after native `merge`: `%td Birth date`; after `cmerge`: `%10.0g` and no label.
- **Location:** `build/cmerge.ado:857-910` (placeholders built from a bare template), `:938-940`, `:651-665` (empty-side branches).
- **Problem:** New variables are created from a type-only template; formats, variable labels and characteristics of the using variables are never applied.
- **Impact:** Every merged-in date/time variable displays as a number, labels vanish; `cf` does not compare metadata, so validation cannot catch it.
- **Proposed fix:** After creating the placeholders, apply `format`, `label variable`, and characteristics from the using file (e.g., captured with `describe, replace`/Mata `st_varformat`/`st_varlabel` while the using data are loaded); share the step with the empty-dataset branches.

<a id="mrg-3"></a>
#### MRG-3 [P1] — Using value-label definitions overwrite same-named master labels, so existing master variables silently display the wrong text
*Component: cmerge*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: master `sex` labeled with `sexlbl` (1 "Male" 2 "Female"); using `grp` labeled with a different `sexlbl` (1 "Treatment" 2 "Control"). After native `merge`, `decode sex` → Male/Female; after `cmerge` → Treatment/Control.
- **Location:** `build/cmerge.ado:547-565` (using labels saved with `label save` and re-run with `, modify`), `:1048-1057` (shared master variables re-pointed to the using's label).
- **Problem:** Native `merge` keeps the master's definition when label names collide (unless `nolabel` changes that); cmerge re-runs the using definitions with `modify`, replacing master mappings, and re-points shared variables.
- **Impact:** Silent change of the meaning of master variables' labels (e.g., standardized names like `yesno`, `lbl`); `decode` produces wrong text. The label do-file round trip may also expand `$` in label text.
- **Proposed fix:** Only define using labels whose names do not exist in the master (native behavior), never re-point shared master variables, and attach using labels only to new variables; transfer label text through Mata (`st_vlload`/`st_vlmodify`) instead of a do-file.

<a id="mrg-4"></a>
#### MRG-4 [P1] — `sorted` is trusted without verification: an unsorted master silently mis-joins (native errors r(5))
*Component: cmerge*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: master ids 3,1,2 (unsorted), using ids 1,2,3: `cmerge 1:1 id using u, sorted` → rc 0 with only 1 match, 2 master-only and 2 using-only rows; native `merge ..., sorted` → "master data not sorted", r(5).
- **Location:** `build/cmerge.ado:150-158`, `:792`, `:986` (only uniqueness is checked; order is not).
- **Impact:** Silent wrong merge for a documented option that users apply when they believe data are sorted.
- **Proposed fix:** Verify sortedness in C (one linear pass) or rely on `: sortedby` as native does, and error r(5) if either side is not sorted.

<a id="mrg-5"></a>
#### MRG-5 [P1] — strL variables: any non-identity merge fails with a misleading r(5); on the "identity" fast path strL keepusing values are written from worker threads (concurrent strL SPI access)
*Component: cmerge*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: master `k L(strL) x` merged with using `k y` → "ctools: writing strL is unsupported ... cmerge: C plugin merge failed (error 5)", r(5) ("not sorted"); master data are restored. Harness (A05-4): with a key-sorted, fully matched master, 8 of 8 strL `SF_sstore` calls ran on pool worker threads, bypassing the shared layer's strL ban (`src/ctools_data_io.c:54-65,99-101` notes concurrent strL access previously caused an allocator abort inside Stata).
- **Location:** `src/cmerge/cmerge_impl.c:1064-1065`, `:1486-1502`; `src/cmerge/cmerge_io.c:66-111` (unchecked `SF_sstore`/`SF_sdata` with a 2,049-byte stack buffer for strL destinations).
- **Impact:** cmerge cannot merge datasets that contain strL variables (common for text fields), with an error that suggests a sort problem; on the fast path, possible crash/corruption and order-dependent behavior.
- **Proposed fix:** Detect strL (master and using) in the ado before any mutation; fall back to native `merge` or convert to str# when lossless; if strL output is supported, write it sequentially on the calling thread using `SF_sdatalen`/`SF_strldata`; map plugin codes to proper Stata errors (MRG-15).

<a id="dstr-1"></a>
#### DSTR-1 [P1] — `cdestring` converts (and with `replace` destroys) variables that native `destring` refuses because they contain non-numeric text, even without `force`
*Component: cdestring / cencode / cdecode*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified. Data `"1" "abc" "2"`: native `destring s, gen(n)` prints "contains nonnumeric characters; no generate" and creates nothing; `cdestring s, gen(n)` creates `n` (byte) with `abc` → missing, rc 0. `cdestring s, replace` turns `s` into a byte variable and the string "abc" is gone; native leaves `s` as str8.
- **Location:** `src/cdestring/cdestring_impl.c:360-370` (unparsable → missing, still stored); `build/cdestring.ado:262-272` (drops the source whenever the plugin returns 0), `:309-312` (only a note). Native reference: `destring.ado:254-266`.
- **Problem:** Without `force`, cdestring behaves as if `force` were specified (plus a message). The help documents this ("the conversion proceeds but a warning is issued", `cdestring.sthlp:94`) while calling the command a drop-in replacement.
- **Impact:** Irreversible loss of every non-numeric string ("n/a", "<5", "10-20", typos) in scripts that rely on destring's refusal as a safety check.
- **Proposed fix:** Match native semantics: if any value fails to parse and `force` is absent, leave that variable untouched (drop its staged output, print native's message) and continue with the others; keep current behavior only under `force`.

<a id="dstr-2"></a>
#### DSTR-2 [P1] — `cdestring, percent` divides only the observations that contain "%"; native divides the whole variable by 100
*Component: cdestring / cencode / cdecode*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: with `force percent`, native gives "1"→.01, "3.5"→.035, "1e3"→10, "0012"→.12; cdestring gives 1, 3.5, 1000, 12 (12 of 26 test values differ). A04 harness: `"0.25"` → .25 vs native .0025.
- **Location:** `src/cdestring/cdestring_impl.c:337-353` (`has_percent` per observation). Native: `destring.ado:186-193,278-284`.
- **Impact:** Silently wrong values (factor of 100) for any column where not every value carries a "%".
- **Proposed fix:** Count "%"-bearing observations per variable in the C pass and divide the whole variable when the count is positive (native semantics).

<a id="dstr-3"></a>
#### DSTR-3 [P1] — `cdestring` writes invalid values into the dataset: ±infinity for "1e400"/"-1e400", and raw doubles above Stata's maximum (which display as `.z_`/garbage missing codes)
*Component: cdestring / cencode / cdecode*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: with `force percent`, `"1.5e308"` was stored as a value that lists as `.z_` (a bit pattern beyond `.z`; native gives `.`). Harnesses (A01, A02, A04): `ctools_parse_double_fast("1e400")` → +inf with ok=1; `"-1e400"` → −inf; `"9e307"` → raw value in Stata's missing range; `"1e+"`/`"1e-"` → 1.
- **Location:** `src/ctools_types.c:566-724` (`ctools_parse_double_fast`: strtod fallback `:690-723` accepts `HUGE_VAL`; exponent `:638-657` accepts a sign with no digits); stored unchecked at `src/cdestring/cdestring_impl.c:348-367`; separator parser `:846` (`pow(10,e)`).
- **Problem:** No finiteness or `< SV_missval` range check after parsing; malformed exponents accepted. Native `real()` returns missing for all of these.
- **Impact:** −inf sorts/sums as a number; +inf and >8.988e307 values are neither valid numbers nor valid missing codes, corrupting summaries and comparisons.
- **Proposed fix:** In `ctools_parse_double_fast` (shared with cimport) return failure when the result is not finite or `>= 8.988465674311579e307`; require ≥1 exponent digit after `e[+-]`. Add the same check to the separator parser.

<a id="dstr-4"></a>
#### DSTR-4 [P1] — `cdestring, dpcomma` treats "." as a thousands separator anywhere: "3.5" → 35, "-.5" → −5, "1.25" → 125, with no warning
*Component: cdestring / cencode / cdecode*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: with `force dpcomma`, native gives missing for "3.5", "-.5", "5." (non-numeric under dpcomma); cdestring gives 35, −5, 5. A04 harness: "1.234" → 1234, "1..2,3" → 12.3, all counted as successful conversions.
- **Location:** `src/cdestring/cdestring_impl.c:239-240` (`grp_sep='.'`); `src/ctools_types.c:807-809` (group separator skipped anywhere, never validated). Native: `destring.ado:228-236`.
- **Impact:** Silent power-of-ten errors when a dot-decimal value appears in a comma-decimal column (common when merging sources).
- **Proposed fix:** Follow native: under `dpcomma` only the first "," becomes the decimal point and any "." makes the value non-numeric unless `ignore(".")` is given; if a grouping mode is desired, make it an explicit option that validates group positions.

<a id="imp-2"></a>
#### IMP-2 [P1] — Delimiter auto-detection ignores the header line and breaks ties in favor of the comma: `;`-delimited files with decimal commas (e.g., `abc;2,5;3,0`) are split on commas
*Component: cimport (delimited text)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: file `a;b / 1,5;x / 2,5;y` → native (default options) splits on `;` (`v1 = "1,5"`, `v2 = "x"`, 3 obs); cimport splits on `,` → variables `ab` (byte) and `v2` (str3) holding `1`/`5;x` and `2`/`5;y`, rc 0. Code: candidates are tried in the order tab, comma, semicolon, colon, pipe over the data rows only, and a later candidate replaces the current best only with a strictly larger consistent field count.
- **Location:** `src/cimport/cimport_impl.c:675-745` (`cimport_auto_detect_delimiter`: candidate order `:677`, header excluded via `data_start = skip_rows + has_header` `:702`, strict `>` tie rule `:731-735`); called at `:892`; `build/cimport.ado:112-132` (default `auto`).
- **Problem:** Whenever every data row contains as many commas as semicolons (typical for European files where k−1 of k columns hold decimal-comma numbers, such as `name;price;qty` rows like `abc;2,5;3,0`), both candidates give the same consistent count and the comma wins; native uses the header line (which has no decimal commas) and chooses `;`.
- **Impact:** European "CSV" exports are silently split in the wrong places, mangling numbers and shifting columns.
- **Proposed fix:** Include the header line in the consistency test (as native does), break ties by preferring the candidate that splits the header into the same number of fields as the data, and treat a comma between digits as a likely decimal separator when `;` is also consistent; add regression files with `;` + decimal commas.

<a id="imp-3"></a>
#### IMP-3 [P1] — Blank lines in CRLF (Windows) files become all-missing observations
*Component: cimport (delimited text)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `id,name,value\r\n1,Alpha,100\r\n\r\n3,Gamma,300\r\n\r\n` → native 2 obs, cimport 4 obs (rows 2 and 4 all missing). LF-only blank lines are handled correctly.
- **Location:** `src/cimport/cimport_impl.c:589`, `:978` (`is_empty_row` requires field length 0); `src/cimport/cimport_parse.c:462-474` (`\r` stays in the field).
- **Impact:** Silent extra observations in common Windows files (including a trailing blank line); a blank CRLF first line also defeats header detection (IMP-9).
- **Proposed fix:** Treat a single-field row whose content is empty after stripping a trailing `\r` as empty; use the same test in header selection.

<a id="imp-4"></a>
#### IMP-4 [P1] — An explicit variable list (`cimport delimited a b using f`) forces `varnames(nonames)`: the header row is imported as data and every column becomes a string
*Component: cimport (delimited text)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `id,value / 1,2.5 / 3,4.5` → native `import delimited a b using f` gives 2 obs (`a` byte, `b` float); cimport gives 3 obs (`a` str2, `b` str5). Harness: explicit names are also case-folded (`ID Val` → `id val`).
- **Location:** `build/cimport.ado:77-81`, `:109`, `:424-433`.
- **Proposed fix:** Keep header auto-detection for an extvarlist and rename variables 1..k after creation without case conversion (native `import_delimited.ado:313-329,851-864`); mirror native's r(103)/r(198) errors.

<a id="imp-5"></a>
#### IMP-5 [P1] — `rowrange()` starting at line 1 is mistranslated for headerless files
*Component: cimport (delimited text)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified on `1,10 / 2,20 / 3,30 / 4,40`: `rowrange(1:2)` native 2 vs cimport 1; `(1:1)` 1 vs 3; `(1)` 4 vs 3; `(:1)` 1 vs 4.
- **Location:** `build/cimport.ado:203-229` (`startrow = max(1, startrow - headerrow)`, `header_only`); `src/cimport/cimport_delimiters.inc:113-131`.
- **Problem:** The ado converts file lines to data rows assuming a header on line 1; C adds one back when no header is detected, and `header_only` with `last=0` removes the upper bound.
- **Impact:** Silently wrong subsets of headerless files.
- **Proposed fix:** Pass original line numbers to C and convert after the header decision.

<a id="imp-6"></a>
#### IMP-6 [P1] — `varnames(#≥2)` with `bindquotes(strict)` loses rows and takes the header from the wrong line (SIMD strict row finder skips row boundaries)
*Component: cimport (delimited text)*

- **Status:** fixed with IMP-1 on 2026-09-27; see [P0 fix status](#p0-fix-status).
- **Evidence:** Stata-verified: `"title line" / "name","value" / "a",1 / "b",2 / "c",3` with `varnames(2) bindquotes(strict)` → native 3 obs; cimport 1 obs (names `b v2`). Harness isolates `cimport_find_next_row_strict`.
- **Location:** `src/cimport/cimport_parse.c:267-286` (AVX2), `:313-326` (SSE2), `:380-401` (NEON); callers `src/cimport/cimport_impl.c:128-134,231-233,254,966-968,1036-1039,1552,1574-1575`.
- **Problem:** In the in-quotes branch the SIMD code jumps to the last quote of a block with odd parity and never examines a newline that follows an earlier closing quote in the same block. The scalar tail is correct, hiding the bug on short inputs.
- **Impact:** Silent row loss and wrong headers; also moves the auto-header probe and variable-label line, and feeds IMP-1.
- **Proposed fix:** Process the first quote in the block (ctz), handle `""` escapes, and continue with out-of-quote logic for the remainder; or use the scalar loop for these rare calls.

<a id="imp-7"></a>
#### IMP-7 [P1] — Windows: files are opened with `CreateFileA` on a UTF-8 path (non-ASCII paths fail) and with `FILE_SHARE_READ` only (files open in Excel fail)
*Component: cimport (delimited text)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Code reading of Win32 semantics (not executed on Windows). The Excel/ReadStat paths already convert with `MultiByteToWideChar(CP_UTF8)`.
- **Location:** `src/cimport/cimport_mmap.c:19-20`.
- **Impact:** `cimport delimited using "C:\Users\José\data.csv"` fails with r(601) after `confirm file` succeeds; common for non-English Windows users.
- **Proposed fix:** `MultiByteToWideChar(CP_UTF8)` + `CreateFileW(..., FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE, ...)`; consider `\\?\` for long paths. (The same issue exists in cexport, see EXP-10.)

<a id="exp-2"></a>
#### EXP-2 [P1] — `cexport excel` writes numbers with at most 15 (often 11–14) significant digits: 16-digit IDs collapse and fractions lose precision
*Component: cexport (delimited + Excel)*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `id = 1234567890123456 + _n` (three distinct values) → after `cexport excel` + `import excel`, all three read back as 1234567890123460; native gives 1234567890123457/…58/…59. `r = _n/7000` → .00014285714286 (native .0001428571428571429); 0.9999999999999999 → 0.99999999999999001 (native 0.99999999999999989). Harness: 537,312 of 580,000 random doubles do not round-trip.
- **Location:** `src/cexport/cexport_xlsx.c:852-1016` (`xlsx_fast_dtoa`): `:891` (`frac_digits = 15 - int_digits`), `:906-908` (rounding carry clamped to `scale-1`), `:872`, `:938` (integers ≥ 1e15 → 15-digit scientific), `:943-948`; call site `:1131-1135`.
- **Impact:** Silent value changes; double-stored identifiers in [1e15, 2^53] become duplicates, breaking keys/merges after an Excel round trip.
- **Proposed fix:** Emit shortest round-trip text (Ryu/Grisu or `%.17g`-then-shorten), use the integer path for all integral |v| ≤ 2^53, and add a `strtod(text) == v` property test across magnitudes.

<a id="reg-4"></a>
#### REG-4 [P1] — Robust/cluster F statistic uses the first `K_keep` entries of the full `b`/`V`: wrong (down to F = 0) whenever an omitted regressor precedes a kept one (SEP24 #14 fixed only for PPML)
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `sysuse auto`, `creghdfe price foreign mpg weight, absorb(foreign) vce(robust)` → e(F) = 0.0366; reghdfe e(F) = 22.10 = `test mpg weight`. (civreghdfe's F is correct in the same setup — verified.)
- **Location:** `build/creghdfe.ado:594-604`.
- **Impact:** Wrong displayed F, `e(F)` and Prob > F under robust/cluster VCE whenever a time-invariant regressor listed first is absorbed — a very common case.
- **Proposed fix:** Build the Wald test from the retained-column index map (or the compact `__creghdfe_beta_1..K_keep` and top-left block of V); set F missing when the compact V is rank-deficient.

<a id="reg-5"></a>
#### REG-5 [P1] — Singleton removal ignores fweights: single-row FE levels with fweight ≥ 2 are dropped (reghdfe drops a level only when its total fweight is 1)
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: reghdfe N = 86, singletons 1, df_a 23, df_r 62, se .15066; creghdfe N = 79, singletons 4, df_a 20, df_r 58, se .15577.
- **Location:** `src/creghdfe/creghdfe_regress.c:483-499`; `src/ctools_hdfe_utils.c:18-137` (unweighted `counts[level] == 1`); civreghdfe calls the same kernel.
- **Impact:** Wrong N, DoF and SEs for fweighted/collapsed data (where single-row cells with large counts are typical).
- **Proposed fix:** Pass compact fweights to `ctools_remove_singletons` and drop a level only when its summed weight equals 1 (unweighted behavior for aw/pw).

<a id="reg-6"></a>
#### REG-6 [P1] — Singleton removal silently stops after 100 passes: long singleton chains survive (creghdfe and civreghdfe)
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: 2×2 core plus a 1,000-edge chain → reghdfe N = 20, se .2958; creghdfe N = 821, se 1.895 (6.4× larger).
- **Location:** `src/creghdfe/creghdfe_regress.c:126`, `:491`; `src/ctools_hdfe_utils.c:107-122`; `src/civreghdfe/civreghdfe_impl.c:373-375`.
- **Impact:** Wrong N, DoF, R² and (drastically) robust/cluster SEs on sparse bipartite data (worker–firm, trade networks), with no warning.
- **Proposed fix:** Iterate to a fixed point (total work is bounded by the number of removals) or use queue-based peeling; if a cap is kept, error when it is hit.

<a id="reg-7"></a>
#### REG-7 [P1] — 32-bit column offsets (`k*N`, `j*N`) overflow once #regressors × N ≥ 2^31 in the shared solver, VCE and PPML/cqreg code (memory corruption on large data)
*Component: creghdfe + shared estimation kernels*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Code reading and UB analysis (A11-2, A12-2): e.g., 30M observations × 72 regressors; cqreg's constant-column fill would write ~17 GB before `X`.
- **Location:** `src/creghdfe/creghdfe_solver.c:421`, `:434`; `src/ctools_matrix.c:233`, `:237`, `:360`; `src/cpplmhdfe/cpplmhdfe_irls.c:647`, `:967-968`, `:1292`–`:1940`; `src/cqreg/cqreg_regress.c:169`, `:876-884`; `src/cqreg/cqreg_blas.c:112-113` (~65 sites in cqreg).
- **Impact:** Heap corruption or crash for large-but-supported problems (the package targets big data).
- **Proposed fix:** Compute all offsets in `size_t` (`(size_t)k * N`), audit with `-Wshorten-64-to-32`/`-Wconversion`, and add a native test with N·K just above 2^31 (can use a mock with sparse allocation).

<a id="iv-2"></a>
#### IV-2 [P1] — `orthog()`/`redundant()` write past the `Z_rest` heap buffer when a tested instrument is repeated or follows a dropped collinear instrument
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Harness with ASan (A10a-5, A10b-1): `heap-buffer-overflow WRITE` in `civreghdfe_compute_cstat` (`civreghdfe_tests.c:2464`) and `civreghdfe_compute_redundant` (`:2928`). In Stata, 20 consecutive `orthog(z3 z3)` calls completed and each displayed a C statistic (3.114) — the corruption is silent rather than an immediate crash.
- **Location:** `src/civreghdfe/civreghdfe_tests.c:2427-2468`, `:2895-2937`; `build/civreghdfe.ado:815-905` (index construction); `src/civreghdfe/civreghdfe_impl.c:1209-1256` (instrument compaction that the indices ignore).
- **Problem:** `Z_rest` is sized `N*(K_iv - n_orthog)` but filled with every column whose mask bit is 0. Duplicated names (`orthog(z3 z3)`) count twice in `n_orthog` but mask one column; indices of instruments after a dropped collinear instrument point beyond the compacted `K_iv`. Either way the fill loop writes ≥ N×8 bytes past the allocation; `redundant()` also reads an uninitialized tail of `Z_test`.
- **Trigger:** `civreghdfe y (x = z1 z2 z3) w, absorb(id) orthog(z3 z3)`; or `gen z1b = 2*z2` and `civreghdfe y (x = z1b z2 z3) w, absorb(id) orthog(z3)`.
- **Impact:** Heap corruption inside Stata (crash or silent corruption); garbage `e(cstat)`/`e(redund)`. The collinear case needs no typo (instrument sets with collinear dummies are common with FE).
- **Proposed fix:** De-duplicate and validate the lists in the ado; pass per-column 0/1 flags aligned with the *expanded* instrument list and compact them together with Z in C; compute `K_rest` from the mask popcount and error if it differs from the requested count.

<a id="iv-3"></a>
#### IV-3 [P1] — `orthog()`, `endogtest()`, `redundant()` and `partial()` indices are positions in the raw (unexpanded) token list: the wrong column is tested/partialled with factor or time-series operators or after collinear drops
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Harness (A10a-6, A10b-2): with `i.grp w, partial(w)` the `3.grp` dummy is partialled instead of `w` and labels shift (`4.grp` shows w's coefficient); `(x = i.grp z3) orthog(z3)` tests `2.grp`; out-of-range `endogtest()` indices silently post 0 with p = 1.
- **Location:** `build/civreghdfe.ado:815-939`; `src/civreghdfe/civreghdfe_tests.c:2449-2454`, `:2690-2695`, `:2907-2912`; `src/civreghdfe/civreghdfe_impl.c:1153-1256`.
- **Impact:** A plausible-looking statistic for the wrong instrument/regressor; wrong partialled model.
- **Proposed fix:** Match names after `fvexpand`/`tsrevar` against the expanded lists, pass per-column flags, remap through the C-side collinearity compaction, and error if a tested column was dropped (shares the fix with IV-2).

<a id="iv-4"></a>
#### IV-4 [P1] — Regressors or instruments fully absorbed by the fixed effects are not omitted: they become unit-variance noise columns with absurd coefficients, and inflate K, J df and CD F
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `age = year − birth_year` with `absorb(id year)` in an unbalanced panel → civreghdfe `b_age = 94,660,130`, `se_age = 1.66e10`, `K = 3`; ivreghdfe omits age. An absorbed instrument (`z2` constant within groups): civreghdfe `sargan_df = 2`, `cd_f = 32.9` vs ivreghdfe `sargandf = 1`, `cdf = .`.
- **Location:** `src/civreghdfe/civreghdfe_impl.c:1021-1121`, `:1209-1255` (equilibration after partialling rescales near-zero columns to unit variance); `src/ctools_ols.c:402-404` (`detect_collinearity` resets stage-1 flags).
- **Impact:** Garbage coefficients/SEs displayed for absorbed regressors, other coefficients perturbed, wrong df and diagnostics.
- **Proposed fix:** Test each partialled column's norm relative to its pre-partialling norm (e.g., `‖x̃‖² < 1e-9·‖x−x̄‖²` → omitted) before equilibration; never rescale such columns; keep stage-1 flags when running `detect_collinearity`.

<a id="iv-5"></a>
#### IV-5 [P1] — `partial()` with aweights/fweights/pweights uses an unweighted FWL projection, changing the coefficients
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `[aw=wt]` model — civreghdfe without partial `b[x] = 2.159455870002` (= ivreghdfe with `partial(w2)`); civreghdfe with `partial(w2)` `b[x] = 2.159219514191`.
- **Location:** `src/civreghdfe/civreghdfe_impl.c:744-797`.
- **Proposed fix:** Weight the partialling regressions (P'WP, P'Wy) consistently with the main estimator.

<a id="iv-6"></a>
#### IV-6 [P1] — Collinear `partial()` variables make the P'P Cholesky fail; the loop breaks and **all** partial controls are silently dropped (omitted-variable bias, rc 0)
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: with `w3 = 2*w2`, `partial(w2 w3)` gives `b[x] = 1.862494` — identical to omitting w2 entirely; ivreghdfe `partial(w2 w3)` gives 1.9214563 (= including w2).
- **Location:** `src/civreghdfe/civreghdfe_impl.c:748-755`, `:846-893`.
- **Proposed fix:** Use a rank-revealing (pivoted) solve that drops only redundant partial columns (as ivreg2 does), or error out; never continue with an unpartialled model.

<a id="iv-7"></a>
#### IV-7 [P1] — Hansen J under `vce(robust)` with aweights/pweights uses `w·e²` instead of `(w·e)²`: every overidentified pweight model reports a wrong J
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: overidentified model with `[pw=pw]` → ivreghdfe `e(j)` = 47.20, civreghdfe `e(sargan)` = 59.09. Harness (A10b-3): civreghdfe J = 9.105 vs reference ivreg2 formula 6.763 (unweighted J matches exactly). Validation never compares J under robust VCE because ivreghdfe posts it as `e(j)` while the helper compares only `e(sargan)`.
- **Location:** `src/civreghdfe/civreghdfe_vce.c:49-50` (`e2 = w * resid[i] * resid[i]`), called from `src/civreghdfe/civreghdfe_tests.c:2074-2077` (also drives the efficient-GMM re-estimation inside the J computation, `:2092-2188`).
- **Proposed fix:** Use `(w·e)²` for aw/pw and `w·e²` for fw in `ivvce_compute_ZOmegaZ_robust`; add a Stata validation comparing `e(sargan)` with ivreghdfe `e(j)` under robust/cluster/HAC.

<a id="iv-8"></a>
#### IV-8 [P1] — C-statistic (`orthog`), endogeneity test (`endogtest`) and redundancy test are always homoskedastic, whatever the VCE
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `orthog(z3)` under robust VCE → ivreghdfe `e(cstat)` = 42.09, civreghdfe 26.61; `endogtest(x2)` robust → ivreghdfe `e(estat)` = 147.89, civreghdfe 162.94. Harness (A10b-4): under robust VCE civreghdfe C = 1.283 vs reference 0.683 (J_full agrees exactly); endogtest/redundant by code reading (`(void)vce_type; (void)cluster_ids;`).
- **Location:** `src/civreghdfe/civreghdfe_tests.c:2436-2438`, `:2580-2610`, `:2635-3123`; `civreghdfe_tests.h:233-282`; callers `src/civreghdfe/civreghdfe_estimate.c:1829-1834`, `:1863-1869`, `:1898-1903`.
- **Problem:** The C-stat subtracts a homoskedastic restricted Sargan from a robust J (an invalid mixture, clamped at 0); endogtest/redundant ignore the VCE entirely, yet are displayed as tests for the chosen VCE.
- **Impact:** Wrong exogeneity/endogeneity/redundancy inference in the most common configuration (`orthog()` with `robust`/`cluster`).
- **Proposed fix:** Pass VCE type, cluster ids (both dimensions), kernel and panel ids; compute C = J_full − J_restricted with the restricted efficient GMM using the relevant block of the full S (ivreg2's `smatrix` approach); endogtest via the same C-stat on the "exogenous" model; redundancy via the robust KP LM machinery.

<a id="iv-9"></a>
#### IV-9 [P1] — With more than one endogenous regressor, `e(idstat)` is always the homoskedastic Anderson LM (labelled Kleibergen-Paap), and the KP Wald F ignores HAC/DK kernels and uses the wrong statistic
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: two endogenous regressors, robust VCE → ivreghdfe `e(idstat)` = 496.40, civreghdfe 724.85 — identical to civreghdfe's iid value (724.85). Harness (A10b-5): `idstat = 5.193145` identical for unadjusted, robust, robust+Bartlett and cluster VCE.
- **Location:** `src/civreghdfe/civreghdfe_tests.c:1493-1496` (never overwritten for K_endog > 1), `:1597-1643`, `:1744-1801`; `build/civreghdfe.ado:1502-1508` (label).
- **Problem:** The KP rk statistic is only implemented for K = 1; the K > 1 "KP Wald" is the minimum eigenvalue of a contracted matrix, not ranktest's statistic (root cause of the known multi-endogenous KP F divergence); Kiefer uses a cluster-robust shat0 (ivreghdfe uses IID).
- **Proposed fix:** Port ranktest's `s_rkstat` for general K (SVD of Θ, λ_q projection, Ω_q from the robust/cluster/2-way/DK/HAC shat0); use IID for Kiefer.

<a id="iv-10"></a>
#### IV-10 [P1] — Driscoll–Kraay time index follows order of first appearance, not calendar order: DK SEs and DK test statistics are wrong in unbalanced panels
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Harness (A10a-7, A10b-6) with the real ID pipeline: DK V = 0.28754 vs 0.36720 with true time order (ratio 0.78); SE 0.14367 vs 0.14117 vs a time-ordered numpy reference. Depends on the first panel's first period.
- **Location:** `src/civreghdfe/civreghdfe_impl.c:201-208`, `:918-932`; `src/ctools_hdfe_utils.c:142-180` (`ctools_remap_cluster_ids` renumbers by first appearance); kernels over ID differences in `civreghdfe_vce.c:687-720`, `civreghdfe_tests.c:634-650,1081-1096,1951-1981`.
- **Proposed fix:** Pass the time variable (or `(t - tmin)/tdelta`) and index the time series by actual period (empty periods contribute zero sums), so lag = time difference.

<a id="iv-11"></a>
#### IV-11 [P1] — `bw()`/`kernel()` without `robust` give plain iid SEs or silently switch to HAC; ivreg2 computes the AC estimator
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: civreghdfe with `bw(3)` alone gives `mreldif(V, V_iid) = 0` (option ignored), while ivreghdfe's AC V differs. Code reading + ivreg2 documentation (A10b-7): `bw(3)` alone leaves `vce_type = 0` and V identical to the unadjusted fit (yet `verbose` prints "HAC: kernel=…").
- **Location:** `build/civreghdfe.ado:485-490`, `:750-760`; `src/civreghdfe/civreghdfe_vce.c:615-622`.
- **Proposed fix:** Treat any `bw`/`kernel` without robust/cluster as AC (generalize `ivvce_compute_kiefer`'s homoskedastic-kernel block to any kernel/bandwidth); set robust only when requested.

<a id="iv-12"></a>
#### IV-12 [P1] — HAC panel/time structure is not passed to the plugin: panel HAC uses the first `absorb()` variable as the panel id; lags are row distances (gaps and missing rows mis-paired); without `absorb()` no panel structure is used
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: after dropping period 10 from every panel, civreghdfe `vce(robust) bw(3)` V differs from ivreghdfe (mreldif 7.5e-6 vs ~1e-18 for the same model without HAC). Harness (A10a-15): SEs change with `absorb(id t)` vs `absorb(t id)`; code reading (A10b-8/9).
- **Location:** `src/civreghdfe/civreghdfe_impl.c:1300-1305` (`hac_panel_ids = fe_levels_c[0]`); `build/civreghdfe.ado:492-504`, `:953-966` (tsset panel/time vars never passed); `src/civreghdfe/civreghdfe_vce.c:455-487`, `:809-906`; `src/civreghdfe/civreghdfe_tests.c:268-310`, `:722-807`, `:1110-1198`, `:1989-2072`.
- **Impact:** Silently wrong HAC/Kiefer SEs and HAC test statistics for common panel data (unbalanced, gaps, or absorb order).
- **Proposed fix:** Pass panel and time variables; build within-panel lag pairs by time value (as `m_omega` does with `tmatrix`); require `tsset` for kernel/bw options.

<a id="iv-13"></a>
#### IV-13 [P1] — `gmm2s`/`cue` ignore the second cluster dimension and HAC/DK/Kiefer structure in both the weighting matrix and the VCE; CUE's robust objective uses `w·e²` for aw/pw
*Component: civreghdfe*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Harness (A10a-8): two-way-clustered gmm2s SEs ≈ 2.4× too small; code reading (A10b-20).
- **Location:** `src/civreghdfe/civreghdfe_estimate.c:69-138`, `:386-392`, `:1336-1444`, `:1486-1538`.
- **Proposed fix:** Build the GMM/CUE S with the same S builder used for the VCE (HAC, DK, two-way, Kiefer) with matching dof; use `(w·e)²` for aw/pw in `cue_objective`.

<a id="qreg-2"></a>
#### QREG-2 [P1] — `denmethod(residual)` does not implement qreg's residual method: it drops every |r| < 1e-8 instead of the K basis observations and clamps s to 1e-10 (SEs ~1e-12 on discrete data); it also accepts `vce(robust)`, which qreg rejects
*Component: cqreg*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: integer outcome with `i.g` → qreg `vce(iid, residual)` sparsity 13.968, se(2.g) 0.155; cqreg `denmethod(residual)` sparsity 1.0e-10, se(2.g) 1.1e-12.
- **Location:** `src/cqreg/cqreg_sparsity.c:641-741`; `src/cqreg/cqreg_vce.c:938-944`.
- **Proposed fix:** Have the crossover return the K basis indices and drop exactly those; error 498 when s < DBL_EPSILON; reject `vce(robust|cluster)` with `denmethod(residual)` (r(184)) as qreg does; add `vce(iid, residual)` comparisons.

<a id="qreg-3"></a>
#### QREG-3 [P1] — `denmethod(kernel)` silently runs the residual method (the kernel code is dead) while `e(denmethod)` reports "kernel"
*Component: cqreg*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `mreldif(V_kernel, V_residual) = 0`; `mreldif(V_kernel, qreg vce(iid, kernel)) = 0.603`. (Known to the maintainer as intentionally kept code, but still user-visible wrong results.)
- **Location:** `src/cqreg/cqreg_regress.c:1036-1064`; `src/cqreg/cqreg_sparsity.c:143-634` (unused kernel statics).
- **Proposed fix:** Implement qreg's kernel VCE (iid and robust; set `e(kernel)`, `e(kbwidth)`) or reject `denmethod(kernel)` with r(198) until it exists.

<a id="bin-2"></a>
#### BIN-2 [P1] — Classic method with `controls()`/`absorb()`: the fit line and `e(coefs)` are estimated on mean-zero residuals but drawn against mean-restored dots (line ~130 units below the points in the example)
*Component: cbinscatter*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `y = 100 + 2x + 3z + e`, `controls(z)` → binscatter `e(y1_coefs)` = (x 1.994, _cons 130.10); cbinscatter `e(coefs)` = (8.7e-15, 1.994).
- **Location:** `src/cbinscatter/cbinscatter_impl.c:536-569`, `:761-768` (means added back to bin means only), `:849-853` (fit on de-meaned arrays); `build/cbinscatter.ado:436-453`, `:503-519`.
- **Impact:** Default usage with controls/absorb produces a graph whose fit line is far from the points; intercepts and qfit/cubic shapes are wrong.
- **Proposed fix:** Add the raw means back to the residual arrays themselves before binning and fitting (binscatter's approach), or re-express the fitted coefficients in shifted coordinates.

<a id="bin-3"></a>
#### BIN-3 [P1] — `method(binsreg)`: the fit line ignores controls and absorbed effects (raw y on raw x)
*Component: cbinscatter*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: same data with `controls(z) method(binsreg)` → `e(coefs)` slope 4.932 (the unadjusted `reg y x` slope) while the dots follow the partial slope ≈ 1.99.
- **Location:** `src/cbinscatter/cbinscatter_impl.c:592-596`, `:849-853`.
- **Proposed fix:** Fit the polynomial jointly with controls/FE and evaluate at w̄ (binsreg's `polyreg()` convention), or draw no line for method(binsreg) and document it.

<a id="bin-4"></a>
#### BIN-4 [P1] — `method(binsreg)` writes covariate-adjusted means to the wrong bins whenever a quantile bin is empty (ties/mass points); with `absorb()` the data are re-binned differently from the x-means; with `discrete`, controls are silently ignored
*Component: cbinscatter*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Harness (A13-4/5): x ∈ {1,2,3} with 50/25/25%, nq(4) → cbinscatter y-means 11.37, 1.36, 21.43 vs binsreg-style 11.82, 21.88, 31.82; `discrete` + `controls(z)` returns exactly the raw per-value means.
- **Location:** `src/cbinscatter/cbinscatter_bins.c:937-949`; `src/cbinscatter/cbinscatter_binsreg.c:47-53`, `:139-153`, `:186-188`, `:209-216`, `:297-306`, `:363-365`, `:421-440`; `src/cbinscatter/cbinscatter_impl.c:742-754`, `:771-831`.
- **Impact:** Silently wrong adjusted means (errors of the order of the y range) whenever x has mass points (zeros, censoring, rounding), or with absorb() at N ≥ 50k; covariate adjustment silently skipped for discrete x.
- **Proposed fix:** Return the exact `bin_ids` used for the x-means from the binning routine, map `bin_id → dense index` once, use it in both adjusters (and for discrete ids); never recompute bins in `impl.c`.

<a id="bin-5"></a>
#### BIN-5 [P1] — `absorb()` values are used directly as array indices: negative codes silently dropped, fractional codes merged, codes ≥ 2^31 undefined/dropped (r(2001)), large codes blow up memory; string FE variables mark out every observation
*Component: cbinscatter*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: `fe = runiformint(-2,2)` → `e(N) = 597` of 1,000 (all negative-code rows dropped as "missing"); 11-digit census-tract codes → r(2001).
- **Location:** `src/cbinscatter/cbinscatter_impl.c:294-303`, `:364-372`; `src/cbinscatter/cbinscatter_resid.c:251-284`, `:304-319`, `:614-615`; `build/cbinscatter.ado:181`.
- **Impact:** Silent data loss / wrong residualization with ordinary FE codings (event time, −1 = unknown); OOM or hangs with large IDs.
- **Proposed fix:** Map each absorb variable to dense 1..L over the estimation sample (in the ado with `egen group()` tempvars, which also handles strings, or in C via sort/hash); drop the `< 1` validity test; `markout` with `strok`.

<a id="bin-6"></a>
#### BIN-6 [P1] — 32-bit overflow in quantile arithmetic for large N × nquantiles (e.g., nq = 100 and N ≥ 21.5M): bin 1 receives the top 3.4% of the data
*Component: cbinscatter*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Harness (A13-7): nq = 1000, N = 3M → 284/999 cutpoints wrong; nq = 100, N = 22M → bin 1 holds 741,393 obs with x-mean 0.703.
- **Location:** `src/cbinscatter/cbinscatter_bins.c:374`, `:712`.
- **Proposed fix:** Compute `(int64_t)q * N` and `(int64_t)(cum-1) * nq`.

<a id="xio-1"></a>
#### XIO-1 [P1] — XLSX sheets without a usable `<dimension>` element (none, `ref="A1"`, single-cell `ref="B3"`), namespace-prefixed worksheets (`<x:sheetData>`), cells without `r=`, or single-quoted attributes import as an **empty dataset with rc 0**; an understated dimension silently drops cells
*Component: cimport/cexport Excel and statistical formats*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified with crafted fixtures (native `import excel` vs `cimport excel`): no dimension → native N=3 k=2, cimport N=0 k=0 rc 0; `ref="A1"` (Apache POI SXSSF style) → same; Excel single-cell `B3` → native N=1, cimport empty; `x:`-prefixed worksheet (.NET Open XML SDK default) → native N=2, cimport N=0; single-quoted attributes and cells without `r` → native N=2, cimport N=0. Harness (A06b-1/2/6/12, A08-1): openpyxl write-only files trigger it; an understated dimension drops all cells outside it.
- **Location:** `src/cimport/cimport_xlsx.c:776-801`, `:834-868`, `:933-935`, `:1051-1065`, `:1114-1135`, `:1420-1546`, `:1521`; `src/io/cio_xlsx.inc:4-48`, `:108-123` (the live route is `cimport.ado` → `_cio_import xlsx` → `scan_xlsx()`).
- **Problem:** Cell storage is sized from `<dimension>` and cells outside it (or all cells, when it is absent/degenerate) are discarded; the fast scanner requires exact, unprefixed, double-quoted, `r`-bearing markup.
- **Impact:** Silent total or partial data loss for valid workbooks produced by common libraries/tools.
- **Proposed fix:** Treat `<dimension>` only as a hint: grow column-major arrays on demand or pre-scan extents; make the scanner prefix-tolerant, accept either quote and any XML whitespace, and infer missing `r` from position; error (never return an empty dataset) when cells were seen but none stored.

<a id="xio-2"></a>
#### XIO-2 [P1] — XLSX shared strings are corrupted: truncated at 4,095 raw bytes (before entity decoding, possibly mid-UTF-8), whitespace-only text dropped ("Hello World" → "HelloWorld", " " → ""), and Japanese phonetic `<rPh>` text appended ("東京" → "東京トウキョウ")
*Component: cimport/cexport Excel and statistical formats*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: 10,000-character cell → 4,095 (native 10,000); 1,000 × `&amp;` → 819 characters (native 1,000); "Hello World" → "HelloWorld"; " " → ""; "東京" → "東京トウキョウ". This leaves SEP20 A08 (string truncation) incomplete for Excel.
- **Location:** `src/cimport/cimport_xlsx.c:62`, `:368-385`, `:410`; `src/cimport/cimport_xlsx_xml.c:169-184`.
- **Proposed fix:** Accumulate `<si>` text in a growable buffer and decode entities once on the full string; deliver all text events (filter by element, not by whitespace); ignore `<t>` inside `<rPh>`; add boundary tests (2045/4095/4096/32767 bytes, multibyte cuts).

<a id="xio-3"></a>
#### XIO-3 [P1] — Legacy `.xls` import corrupts formula and boolean/error cells: text formulas → 0, `#N/A` → 42, `""` formulas → 0, TRUE/FALSE → the string "bool", error cells → "error"
*Component: cimport/cexport Excel and statistical formats*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified on `formulas.xls`: native A = "hello"/"world", B = ./2, C = ./3, D = 1/0, E = ./5; cimport A = 0/0, B = 42/2, C = 0/3, D = "bool"/"bool" (str4), E = "error"/"5" (str5).
- **Location:** `src/io/cio_xls.inc:5-11`, `:161-171`, `:187-234`.
- **Proposed fix:** Classify cells explicitly: FORMULA with string result (resid 0/3) → string; boolean result → 0/1; error result → missing; BOOLERR → 0/1 or missing; only NUMBER/RK/MULRK and numeric formula results are numeric.

<a id="xio-4"></a>
#### XIO-4 [P1] — `cexport dbase` fails with r(108) on ordinary values (0.05, −0.1, 0.001, 1.5e-5, float −0.3, bytes ≤ −100)
*Component: cimport/cexport Excel and statistical formats*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: native `export dbase` rc 0, `cexport dbase` r(108) for both the rate and byte fixtures.
- **Location:** `src/io/cio.c:1025-1028`, `:1135-1141` (`"%.18g"` text overflows the fixed width-20 / width-3 fields).
- **Proposed fix:** Emulate native's `%20.18g` fitting (drop the leading 0 for |x| < 1, switch to exponent form or reduce precision until it fits) or size fields from the data; give bytes width 4.

<a id="xio-5"></a>
#### XIO-5 [P1] — Worksheet parts are chosen by tab position (`sheet{N}.xml`) instead of workbook relationships: with a chart sheet or reordered parts the wrong sheet is imported silently or `sheet()`/`describe` fail with r(610)
*Component: cimport/cexport Excel and statistical formats*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: chartsheet workbook `sheet("Other")` → native N=1, cimport r(610); `cimport excel ..., describe` r(610). Harness (A06b-3): POI-style reordered parts import another sheet.
- **Location:** `src/cimport/cimport_xlsx.c:284`, `:1893-1894`; `src/cimport/cimport_xlsx_zip.h:85-88`; `src/io/cio_excel.inc:68-94`.
- **Proposed fix:** Parse `workbook.xml.rels`, map `r:id` → target, accept only worksheet relationships, and give clear errors for chart/dialog/macro sheets.

<a id="xio-6"></a>
#### XIO-6 [P1] — Excel columns mixing date, datetime and General cells import with mixed units (days and milliseconds under one format; General cells left as raw serials); the XLS reader takes the date/datetime kind from the last date cell
*Component: cimport/cexport Excel and statistical formats*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified on `mixed_dates.xlsx`/`.xls`: native shows `15mar2023 00:00:00`, `16mar2023 12:00:00` (mixed date/datetime column), `15mar2023`/`16mar2023` (date then General) and `15mar2023 06:00:00`/`16mar2023 00:00:00`; cimport XLSX shows `1/1/1960 12:00` (a day count under a datetime format), `17mar2083` (a General serial treated as days) and `1.99e+12`; cimport XLS shows `3/15/2023` for all three columns (time of day dropped). Harness (A06b-14, A08-4).
- **Location:** `src/io/cio_xlsx.inc:152-164`, `:181-219`; `src/io/cio_xls.inc:161-171`, `:224-233`.
- **Proposed fix:** One decision per column (datetime if any datetime cell), convert every numeric cell in the column to that unit with the same epoch logic, shared between the XLS and XLSX bridges.

<a id="psm-1"></a>
#### PSM-1 [P1] — `cpsmatch` matches observations psmatch2 excludes (missing outcome, missing covariates/pscore, outside `e(sample)`): wrong ATT
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: 4-row example (one treated unit with missing outcome) → psmatch2 ATT = 10, cpsmatch ATT = −40. With 40 missing covariate values, 18 treated rows with missing `_pscore` are marked on support.
- **Location:** `build/cpsmatch.ado:139-152`, `:163-176`, `:211`, `:236-261`; `src/cpsmatch/cpsmatch_impl.c:776-829`, `:838-855`, `:1178-1199` (only missing *treatment* is skipped; missing pscore = `SV_missval` treated as a huge number).
- **Impact:** Silent wrong ATT, `_weight`, `_support` whenever outcome, covariates or a supplied pscore have missing values — very common.
- **Proposed fix:** Build `touse` like psmatch2 (`marksample` over treatment + covariates, `markout` outcome, `replace touse = e(sample)` after estimation, `markout` `_pscore`) and pass `if `touse'`; skip missing pscores defensively in C.

<a id="win-1"></a>
#### WIN-1 [P1] — `cwinsor ..., trim` converts extended missing values (.a–.z) to `.` on Apple Silicon (SIMD path), destroying them in the default in-place mode
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: 100 rows with `.a` in rows 1–10 → after `cwinsor x, cuts(5 95) trim`, `count if x == .a` = 0; `winsor2 ..., trim replace` keeps all 10. Harness (A01-3, A15-3): position/platform-dependent (NEON/AVX2 lanes only; x86 release builds without AVX2 are unaffected).
- **Location:** `src/ctools_simd.h:271-328` (`ctools_simd_replace_oob`: AVX2 `:277-297`, NEON `:298-318` vs scalar tail `:322-327`); called from `src/cwinsor/cwinsor_impl.c:195-198`, `:569-570`, `:611-612`; default replace target `build/cwinsor.ado:103`.
- **Proposed fix:** Replace only valid out-of-bounds lanes (`keep = !valid | in_bounds`), mirroring the scalar tail; add a `.a/.b` trim test run on ARM.

<a id="win-2"></a>
#### WIN-2 [P1] — cwinsor's percentile selection is quadratic on tied data (zero-inflated or dummy variables): 235× slower than winsor2 at 1M rows, effectively hanging at 10M
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: 1,000,000 rows with 50% zeros → `cwinsor` 18.8 s vs `winsor2` 0.08 s. Harness: 2M rows with 30% zeros 28.5 s; extrapolates to ~1 hour per variable at 10M. The plugin cannot be interrupted.
- **Location:** `src/ctools_select.h:35-89` (Lomuto partition with `<=` peels ties 1–2 per pass); callers `src/cwinsor/cwinsor_impl.c:57-160`; a private copy with the same problem in `src/crangestat/crangestat_impl.c:128-168`. (`src/cqreg/cqreg_linalg.c:273` already uses a 3-way partition.)
- **Impact:** The canonical winsorizing inputs (financial ratios with many zeros, counts, dummies) make the "fast" command orders of magnitude slower than native.
- **Proposed fix:** Three-way (Dutch-flag) partitioning or Hoare partitioning that stops on equal keys (plus an introselect fallback) in the shared helper and crangestat's copy; add a large all-ties performance test.

<a id="rng-2"></a>
#### RNG-2 [P1] — crangestat `count`/`sum` for empty or all-missing windows differ from rangestat, and `count` depends on which path runs (group size < 64 vs ≥ 64)
*Component: crangestat / crangejoin / cipolate / csplit*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: all-missing look-back windows in a 20-row dataset → crangestat count/sum `.`/`.`, rangestat 0/0; a 64-row group with `interval(t -3 -1)` → obs 1 count 0, rangestat `.`. Harness (A14-2).
- **Location:** `src/crangestat/crangestat_impl.c:773-777` (prefix path), `:1057` (simple path); path selection `:1818`, `:1962`, `:2078`, `:2185`.
- **Problem:** rangestat leaves results missing only when the window has no rows; for rows that are all missing it reports count 0 and sum 0. crangestat returns count 0 for empty windows in groups ≥ 64 rows, missing count/sum for all-missing windows in smaller groups, and missing sum in both.
- **Trigger:** `set obs 20`, `gen t = _n`, `gen x = cond(_n <= 5, ., 1)`, `crangestat (count) cn=x (sum) cs=x, interval(t -2 0)` → obs 1–3 `.`/`.` (rangestat 0/0); with 64 rows and `interval(t -3 -1)`, obs 1 count 0 (rangestat `.`).
- **Impact:** Wrong values in look-back and leave-one-out designs (first period of each panel, singleton groups); small and large groups follow different conventions within one call; `if cn > 0` silently changes meaning because missing > 0.
- **Proposed fix:** Track the window row count separately from the non-missing count (empty window ⇒ all missing; otherwise count = non-missing count, sum = 0 when none) in every path; test group sizes 63/64/65.

<a id="rng-9"></a>
#### RNG-9 [P1] — crangestat fails with r(198) whenever *any* variable in the dataset has a name of 24 or more characters, because it creates a `__varpos_<varname>` local for every variable
*Component: crangestat / crangejoin / cipolate / csplit*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Stata-verified: adding an unrelated variable `a_very_long_variable_name_here_x` to a 20-row dataset makes `crangestat (mean) m=v, interval(k -1 1)` fail with "__varpos_a_very_long_variable_name_here_x invalid name", r(198).
- **Location:** `build/crangestat.ado:263-291` (loop over all dataset variables building `local _varpos_`v'`).
- **Problem:** Local macro names are limited to 31 characters; `_varpos_` (8 characters) plus a variable name of 24 or more characters exceeds the limit.
- **Impact:** crangestat is unusable on any dataset containing a long variable name — common after `import` of files with descriptive headers — even if that variable is not involved in the call.
- **Proposed fix:** Compute plugin positions only for the variables passed (e.g., `: list posof "`v'" in allvars`) instead of creating a local per dataset variable; `unab` key/source names.

<a id="bld-1"></a>
#### BLD-1 [P1] — CI cannot stage or publish a package: the new codec licence file (`ctools-codecs-LICENSE.txt`) is omitted from the CI copy step and the release-gate test (and is untracked in git)
*Component: Build, packaging, CI and validation*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Reproduced (A16-1): `python3 -B validation/test_release_gate.py` → `ValueError: missing or invalid package entry: ctools-codecs-LICENSE.txt`; a replica of the CI package step fails `check_package.py` with "missing package file".
- **Location:** `build/ctools.pkg:40`; `.github/workflows/build.yml:87`, `:259-264`; `validation/test_release_gate.py:56-63`.
- **Impact:** build-linux fails, so `package`, `commit-plugins` and `release` never run; merging dev to main while CI is red ships new .ado files with stale x86-mac/Windows/Linux plugins (and there is no version handshake, CORE-3).
- **Proposed fix:** Stage with `scripts/stage_package.py --platform all`; copy every manifest `f` entry in the test; commit the licence file and include it in CI path filters (or generate it via `sync_release.py`).

<a id="bld-2"></a>
#### BLD-2 [P1] — Stale native regression tests no longer compile (`test_p2_native.py`, `test_cimport_options_native.py`); one of them breaks the Linux CI gate and the SEP20 A08/A22/A29 regressions are dead
*Component: Build, packaging, CI and validation*

- **Status:** fixed on 2026-09-28; see [P1 fix status](#p1-fix-status).
- **Evidence:** Compiled locally (A16-2, A06a-24): `unknown type name 'CImportSPIStoreTask'`, `no member named 'fields' in 'CImportParsedRow'`.
- **Location:** `validation/test_p2_native.py:187-214`; `validation/test_cimport_options_native.py:115`; `.github/workflows/build.yml:79`; current APIs `src/cimport/cimport_impl.c:1614-1616`, `src/cimport/cimport_parse.h:33-54`; `docs/IO_PARITY.md:177` still claims they pass.
- **Proposed fix:** Port to `CImportStoreTile`/`cimport_store_tile_worker`/`cimport_fill_row()`; run every `validation/test_*_native.py` in CI via discovery so API drift fails immediately.

### P2 findings

<a id="core-1"></a>
#### CORE-1 [P2] — Unbounded `threads(#)` can terminate Stata: the OpenMP runtime aborts the process when thread creation fails
*Component: Core infrastructure and shared data layer*

- **Evidence:** Harness against the project's libomp (A01-2): `omp_set_num_threads(13000)` + one parallel region → "OMP: Error #34: System unable to allocate necessary resources for OMP thread", SIGABRT (this Mac's `kern.num_taskthreads` is 12,288). The ctools path is by code reading.
- **Location:** `src/ctools_plugin.c:114-168`, `:292-294` (`parse_threads_arg` accepts 1..INT_MAX); `src/ctools_threads.c:93-104`, `:456-458`; every wrapper forwards the value unchanged.
- **Impact:** A mistyped option (`threads(20000)`, or `threads(5000)` on a Linux node with `ulimit -u 4096`) kills the Stata session and unsaved work; large values also inflate per-thread allocations.
- **Proposed fix:** Clamp once in `parse_threads_arg`/`ctools_set_max_threads` (e.g., to a small multiple of the detected CPUs) or reject values above a documented maximum with r(198); size the persistent pool from the clamped value.

<a id="core-2"></a>
#### CORE-2 [P2] — The `__ctools_strw` width hint is read into a fixed 16 KiB buffer; a width cut mid-number makes valid loads fail on very wide datasets (SEP24 #13 analogue)
*Component: Core infrastructure and shared data layer*

- **Evidence:** Harness under truncating `SF_macro_use` semantics (A02-2, A04-12): 8,190 numeric variables then a `str244` whose "244" straddles the cut → hint 24 → "failed to read Stata data", r(920). Affects csort, cmerge, cexport excel, cencode, cdestring. No corruption (widths are validated).
- **Location:** `src/ctools_data_io.c:164-187` (`char buf[16384]`, no explicit terminator); consumers `:262-291`, `:1470-1497`, `:1799-1809`; producer `build/_ctools_strw.ado`.
- **Proposed fix:** Allocate from the macro's actual length (query first or read in chunks) or pass widths per variable; treat truncated metadata as "no hint".

<a id="core-3"></a>
#### CORE-3 [P2] — `_ctools_load` reads the plugin identity but never enforces ado/plugin compatibility; a stale or pinned older plugin silently runs against the new argument protocol
*Component: Core infrastructure and shared data layer*

- **Evidence:** Code reading (A01-10, A16-8); `DEVELOPERS.md:1059` claims identity checks exist; only cpplmhdfe checks an API number.
- **Location:** `build/_ctools_load.ado:49-60`; `src/ctools_plugin.c:75-86` (image pinned with RTLD_NODELETE, so `ctools, update`/`discard` keep the old plugin resident).
- **Impact:** After an update, an old plugin can misinterpret post-SEP24 tokens (e.g., treat a delimiter token as a filename) instead of failing cleanly; very old plugins fail with the misleading "unknown command 'version'".
- **Proposed fix:** Compare the plugin's version and an explicit protocol number with the wrapper's expected value in `_ctools_load` and exit with a clear "restart Stata / reinstall" message on mismatch.

<a id="sort-6"></a>
#### SORT-6 [P2] — Timsort's parallel merge sizes its work split from a team measured in a different parallel region (SEP24 #9 defect class, missed for timsort)
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Native harness with `OMP_DYNAMIC=true KMP_DYNAMIC_MODE=random`: 16/20 runs at N=200k (all-distinct keys) produced duplicated indices; 0/20 after the fix. Plausible but not observed with default OpenMP settings.
- **Location:** `src/ctools_sort_timsort.c:600-615`, `:665-677`, `:729-741`.
- **Problem:** `num_threads` is read in one `parallel`/`single` region; a second `#pragma omp parallel` assigns chunk `tid` of `total/num_threads`. A smaller second team leaves chunks unwritten (stale `order[]` → duplicates/lost rows).
- **Impact:** Same corruption as SORT-1 under dynamic OpenMP team sizing.
- **Proposed fix:** Read `omp_get_num_threads()` inside the partitioning region, or use `#pragma omp for` over a fixed number of logical chunks.

<a id="sort-8"></a>
#### SORT-8 [P2] — `csort x x, stream(#)` (repeated key) overflows the heap and leaves the key sorted but all other columns unsorted
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Native harness (A04): ASan heap-buffer-overflow WRITE at `csort_stream.c:963`; without ASan rc=920 with 1000/1000 rows misaligned (same data with `nostream`: correct). Not run in Stata (heap corruption).
- **Location:** `src/csort/csort_stream.c:943-977`; `build/csort.ado:13,115-124` (duplicate names passed as duplicate indices).
- **Problem:** `nvars_nonkey = nvars - nkeys` counts duplicate keys, so the non-key arrays are one element short per duplicate; the fill loop writes past both mallocs; the sanity check runs after the overflow and after Phase 4 already wrote the sorted keys (SORT-9).
- **Trigger:** `set obs 1000`, `gen double x = 1001-_n`, `gen long id = _n`, `gen str8 s = "r"+string(_n)`, `csort x x, stream(2)` → r(920) with `x` sorted and `id`/`s` not. Also reachable via auto-streaming (N > 10M, > 30 vars) with programmatic varlists containing a repeat.
- **Impact:** Heap corruption in Stata's process plus a scrambled dataset.
- **Proposed fix:** De-duplicate sort keys in the ado (`list uniq`) and in `parse_sort_vars`; compute non-keys with an `is_key[]` bitmap; validate before writing anything.

<a id="sort-9"></a>
#### SORT-9 [P2] — Streaming csort writes sorted keys before the payload permutation can fail; a later failure leaves every row scrambled (and returns rc 0 in non-OpenMP builds)
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Native harness (A04 `h_csort_strfail.c`): failing one string flat-buffer `calloc` → OpenMP build r(920) with the string column misaligned on 1000/1000 rows; non-OpenMP build **rc=0** with the same misalignment.
- **Location:** `src/csort/csort_stream.c:885-905` (Phase 4 writes keys) precede Phase-5 allocations at `:943-951`, `:378-382`, `:450-469`, `:240-241`; `:743-748` ignores `stream_permute_string_var` failure in the non-OpenMP branch; `src/csort/csort_impl.c:335-347`.
- **Problem:** Streaming mode (explicit `stream()` or automatic for N > 10M with > 30 variables) is non-transactional; memory exhaustion after Phase 4 leaves keys permuted and other columns in original order.
- **Impact:** A user who sees "r(920) insufficient memory" reasonably assumes data are untouched, but every row's keys now belong to other observations.
- **Proposed fix:** Allocate and validate all Phase-5 resources before writing; write keys last; on mid-way failure re-apply the inverse permutation to already-written columns; never ignore per-column failures; report whether data were modified.

<a id="sort-10"></a>
#### SORT-10 [P2] — csort rejects any dataset containing a strL variable (even a non-key payload), after loading and sorting everything, with a misleading r(920)
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Stata-verified: with one `strL` payload variable every `csort` call (all algorithms, numeric or string keys) fails: "ctools: writing strL is unsupported; recast to str2045 when lossless / csort: failed to store data (error 5)" → r(920). Data were left unchanged (verified with `cf`). Native `sort` handles strL.
- **Location:** `src/ctools_data_io.c:746-748` (strL rejected at store validation); `build/csort.ado` (no up-front check); `src/csort/csort_impl.c:335-347` (maps to 920).
- **Impact:** csort is unusable on any dataset with a strL column (common for free-text survey data), after wasting a full load+sort; r(920) wrongly suggests a memory problem.
- **Proposed fix:** Detect strL in the ado before the plugin call; either (a) fall back to native `sort` (or sort with the plugin into a permutation variable and apply it in Stata for strL columns), or (b) reject up front with a type error (r(109)/r(198)) and a clear message. Document the limitation in `csort.sthlp`.

<a id="mrg-6"></a>
#### MRG-6 [P2] — Repeating a key variable (`cmerge 1:1 id id using ...`) scrambles the master key and silently mis-joins
*Component: cmerge*

- **Evidence:** Stata-verified: 5 master and 5 using rows with identical ids → native `merge 1:1 id id` matches all 5; `cmerge 1:1 id id` reports 2 matched, 3 master-only, 3 using-only, rc 0.
- **Location:** `src/cmerge/cmerge_impl.c:772-796` (the master key permutation is applied once per listed key).
- **Proposed fix:** De-duplicate the key list in the ado (`list uniq`) and reject duplicates at the C boundary.

<a id="mrg-7"></a>
#### MRG-7 [P2] — `replace` without `update` is accepted and overwrites master values (native rejects it)
*Component: cmerge*

- **Evidence:** Stata-verified: `merge 1:1 id using u, replace` → r(198); `cmerge 1:1 id using u, replace` → rc 0 (harness: master x=5 becomes 99).
- **Location:** `build/cmerge.ado:992-994`; `src/cmerge/cmerge_io.c:42-52,86-99`.
- **Proposed fix:** `if "`replace'" != "" & "`update'" == ""` exit 198 with native's message.

<a id="mrg-8"></a>
#### MRG-8 [P2] — Observation order after the merge differs from `merge`, and the sort marker is not managed
*Component: cmerge*

- **Evidence:** Stata-verified: in 19 of the 27 merge configurations whose content matched native `merge`, the observation order differed; in 6 configurations native left the data sorted by the key (`: sortedby` = key) while cmerge left no sort marker. Harness (A05-5): cmerge interleaves using-only rows by key; native appends them after the first N (master) rows.
- **Location:** `src/cmerge/cmerge_join.c:152-306`; `src/cmerge/cmerge_impl.c:974-1033` (`preserve_order` defaults to 0); `build/cmerge.ado:959-963,1034-1037`.
- **Impact:** Code relying on native's documented layout ("the original master data can be found in the first N observations") or on `_n` after merge behaves differently; a stale non-key sort marker could survive plugin writes (see SORT-7).
- **Proposed fix:** Emit master rows in sorted-master order followed by using-only rows (native layout); for 1:m/m:m put each master row's first match in the first N rows; keep key-interleaving as an opt-in; set/clear the sort marker exactly as `merge.ado:225-230` does.

<a id="mrg-9"></a>
#### MRG-9 [P2] — Fixed plugin argument buffers (8 KB / 4 KB) truncate silently on wide datasets; shared variables are then treated as new and master values are overwritten
*Component: cmerge*

- **Evidence:** Harness (A05-6): master with 16,026 variables, 924 new + 100 shared using variables → 68 of 100 shared variables had master values replaced; 4 KB load buffer leaves uninitialized payload indices (r(459)).
- **Location:** `src/cmerge/cmerge_impl.c:486-488` (`static char args_copy[8192]`), `:619-632`, `:496,503` (static arrays retain previous values), `:143-150,194-200` (`char args_copy[4096]`).
- **Impact:** Silent corruption of master data in very wide merges (≥ ~10,000 master variables with ~1,000 kept variables).
- **Proposed fix:** Parse the caller's (already private) buffer directly or `strdup`; cross-check parsed counts against declared counts and fail on mismatch; zero-initialize per-call arrays; better, pass index lists through Stata macros/Mata.

<a id="mrg-10"></a>
#### MRG-10 [P2] — The hard-coded scratch variable `_merge_temp_filter` can receive `_merge` codes, filter on user data, and then be dropped
*Component: cmerge*

- **Evidence:** Code reasoning (A05-9).
- **Location:** `build/cmerge.ado:881-887` (fixed name), `:948-952` (`merge_var_idx = c(k)` assumes it is last), `:1079-1098` (`keep if _merge_temp_filter == …`, then `drop`).
- **Trigger:** master already has `_merge_temp_filter`; `cmerge 1:1 id using u, nogenerate keep(match)` keeps the wrong rows and deletes the user's variable.
- **Proposed fix:** Use a `tempvar` and pass its index from `st_varindex` instead of assuming `c(k)`.

<a id="mrg-11"></a>
#### MRG-11 [P2] — `generate`/`_merge` existence check matches abbreviations: a master variable such as `_merge_prev` causes a spurious r(110)
*Component: cmerge*

- **Evidence:** Code + [P] confirm semantics (A05-10).
- **Location:** `build/cmerge.ado:165-173` (`capture confirm variable`), and the key check at `:375`.
- **Proposed fix:** `confirm new variable` / `confirm variable ..., exact`.

<a id="mrg-12"></a>
#### MRG-12 [P2] — Variable notes are copied incorrectly: new notes are invisible, and a shared variable's first master note is overwritten; using variables with names ≥ 27 characters build an illegal local name
*Component: cmerge*

- **Evidence:** Code + `notes.ado` semantics (A05-11); name-length issue PLAUSIBLE (A05-12, not yet run).
- **Location:** `build/cmerge.ado:533-538`, `:1039-1046` (only `char var[note1]`, no `note0`, `_dta` notes ignored); `:536`, `:1042-1043` (`local `var'_note`).
- **Proposed fix:** Enumerate using notes (`note0..noteN`, including `_dta`) and append missing ones via `notes`; store them in indexed locals or Mata rather than name-derived locals.

<a id="mrg-13"></a>
#### MRG-13 [P2] — Valid native syntax is rejected: `keepusing()` wildcards/abbreviations and keep/assert words or abbreviations (`mat`, `mas`, `3 4`, `match_update`)
*Component: cmerge*

- **Evidence:** Stata-verified for `keep(match_update)`, `keep(match_conflict)`, `keep(4 5)` (r(198)); code for wildcards (A05-14).
- **Location:** `build/cmerge.ado:353-354`, `:364-365`, `:435-447`, `:176-199`.
- **Proposed fix:** Expand `keepusing` with `unab` against the loaded using frame (dropping keys, de-duplicating); reuse native's result-word mapping.

<a id="mrg-14"></a>
#### MRG-14 [P2] — Undocumented hard limits: at most 1,024 non-key using variables and 32 keys (r(198))
*Component: cmerge*

- **Evidence:** Code (A05-15).
- **Location:** `src/cmerge/cmerge_impl.c:59-60,166-169,184-187,526-530,567-571`; `build/cmerge.ado:448-470`.
- **Proposed fix:** Allocate index arrays dynamically (after fixing MRG-9) or document the limits and fall back to native `merge`.

<a id="dstr-5"></a>
#### DSTR-5 [P2] — The separator-aware number parser (cdestring `dpcomma`, cimport `decimalseparator()`/`groupseparator()`) mis-scales values with ≥18 significant digits and is not correctly rounded
*Component: cdestring / cencode / cdecode*

- **Evidence:** Harness (A01 `sep_harness.c`, 2M random inputs vs strtod): 209,587 wrong by ≥1 order of magnitude, 184,634 off by ~1 ulp; UBSan signed overflow on long exponents. Examples: `12345678901234567890,5` → 1.23e16 (true 1.23e19); `123456789012345678,99` → 1.23e15 (true 1.23e17); `0,000000000000000001` → 0.
- **Location:** `src/ctools_types.c:733-896` (`ctools_parse_double_with_separators`): 18-digit cap `:798`, leading zeros consume the budget `:798-801`, unchecked exponent accumulation `:825`, exponent path `:832-846`, no-exponent path `:860-883`. Callers: `src/cdestring/cdestring_impl.c:239-240,348,356`; `src/cimport/cimport_parse.c:735`.
- **Impact:** Silently wrong magnitudes for long-digit values with European formatting or thousands separators; 1-ulp divergence from native for 17-18-digit values.
- **Proposed fix:** Normalize into a '.'-decimal buffer (drop group separators, map the decimal separator) and call the correctly rounded `ctools_parse_double_fast`, as `src/cimport/cimport_locale.inc:163` already does; delete the bespoke scaling code.

<a id="dstr-6"></a>
#### DSTR-6 [P2] — `cdestring` turns ".a"–".z" into system missing and reports them as non-numeric; "NA"/"NaN" (any case) are silently accepted as missing
*Component: cdestring / cencode / cdecode*

- **Evidence:** Stata-verified: native `destring` keeps ".a" as `.a` (with or without `force`); cdestring stores `.` and counts it as a failure ("1 observation(s) contained nonnumeric characters; converted to missing"). Harness: "NA", "nan" → `.` without being counted.
- **Location:** `src/ctools_types.c:584-600` (only "."/"NA"/"NaN" recognized); `src/cdestring/cdestring_impl.c:360-370`. Native: `destring.ado:224-226`.
- **Impact:** Loss of survey extended-missing codes (refused / don't know); NA strings converted where native would stop.
- **Proposed fix:** Recognize exactly native's missing tokens ("", ".", ".a"–".z" mapped to the corresponding extended missing values); treat "NA"/"NaN" as non-numeric (or behind an explicit option); `src/cimport/cimport_parse.c:723-731` already has the extended-missing mapping.

<a id="dstr-7"></a>
#### DSTR-7 [P2] — `cdestring` and `cencode` drop variable labels and characteristics; `replace` moves the variable to the end of the dataset
*Component: cdestring / cencode / cdecode*

- **Evidence:** Stata-verified: `cdestring mpg_s, replace` loses the label "Mileage as text"; `cencode make, gen(mk)` has no variable label while `encode` gives "Make and model". Source comparison with native `destring.ado:268,324-337` and [D] encode (A04-9).
- **Location:** `build/cdestring.ado:132-160` (outputs generated at the end), `:262-272` (`drop` + `rename`, no label/char copy); `build/cencode.ado:124-128,161,222-225` (C path: `generate long destvar = .`, no label copy).
- **Trigger:** `sysuse auto`, `tostring mpg, gen(mpg_s)`, `label variable mpg_s "Mileage as text"`, `cdestring mpg_s, replace` → empty label and `mpg_s` moved last; `cencode make, gen(mk)` → no variable label (native: "Make and model").
- **Impact:** Silent loss of documentation metadata; positional varlists (`v1-v20`) change meaning.
- **Proposed fix:** Order outputs after the source, copy the variable label, and `char rename`/copy characteristics as native does; for cencode's C path copy the label and keep position on `replace` (as cdecode does).

<a id="dstr-8"></a>
#### DSTR-8 [P2] — cencode's label-file path is split at whitespace in C and `run` unquoted: cencode fails (or writes a stray file) when Stata's temp directory contains a space (SEP24 #6 analogue)
*Component: cdestring / cencode / cdecode*

- **Evidence:** Harness (A01 `misc_harness.c`, A04 `h_cencode_path.c`): with path `.../tmp dir/labels.do` the plugin returns rc 0, the intended file does not exist, and a stray file `.../tmp` containing the label do-file is created.
- **Location:** `build/cencode.ado:185-186,196` (`labelfile=`path'` inside the tokenized argument string), `:212` (`run `__labelfile'` unquoted); `src/cencode/cencode_impl.c:91-92`; `src/ctools_parse.c:40-68` (value copied up to first space, truncation still returns success).
- **Trigger:** Start Stata with `STATATMP="/tmp/ctools tmp"`; `sysuse auto`; `cencode make, gen(mk)` → r(601)/r(920) and a stray `/tmp/ctools` file.
- **Impact:** cencode unusable on such systems (Windows profiles with spaces when TEMP is not 8.3); may overwrite an unrelated file.
- **Proposed fix:** Pass the path out of band (a local read with `SF_macro_use`, as the SEP24 #6 fix does for export) and quote it in the ado (`run `"`__labelfile'"'`); make `ctools_parse_string_option` fail on truncation or whitespace.

<a id="dstr-9"></a>
#### DSTR-9 [P2] — `cdestring` has no rollback; a plugin failure leaves new `generate()` variables behind (e.g., strL values over 2,045 bytes, which native handles)
*Component: cdestring / cencode / cdecode*

- **Evidence:** Source control-flow reading (A04-11).
- **Location:** `build/cdestring.ado:148-160` (outputs created before the plugin), `:247-250` (`plugin call` not captured, so the `plugin_rc` handling at `:273-278`, `:374-376` is dead); `src/ctools_data_io.c:87-97,123-127` (strL > 2045 rejected).
- **Trigger:** `gen strL s = "1"`, `replace s = "2" + 3000*" " in 2`, `cdestring s, gen(n)` → r(920) and a leftover all-missing `n`; retry → r(110).
- **Impact:** Half-completed state; the help's strL support claim holds only for ≤2,045-byte values.
- **Proposed fix:** Wrap in `preserve`/`capture noisily`/`restore` like cencode/cdecode (or stage as tempvars, rename on success); fall back to native `destring` for long strL; fix the help.

<a id="imp-8"></a>
#### IMP-8 [P2] — `NaN`/`Infinity` spellings make the whole column a string (native: numeric missing); overflow values such as `1e999` are stored as raw ±inf
*Component: cimport (delimited text)*

- **Evidence:** Harness (A06a-7): `0.5, NaN, -Infinity, 2.25` → `str9` (native numeric with `.`); a `byte` column received `-inf`.
- **Location:** `src/cimport/cimport_parse.c:641-703`; `src/cimport/cimport_impl.c:370-376`, `:1641-1649`.
- **Proposed fix:** Accept Java's special spellings as numeric system missing; map any non-finite parsed value to `SV_missval` before storing (see DSTR-3 for the shared parser).

<a id="imp-9"></a>
#### IMP-9 [P2] — Leading blank line or header-only file: header handling diverges from native
*Component: cimport (delimited text)*

- **Evidence:** Stata-verified for a header-only file (`a,b,c\n`): native creates `a b c` (byte) with 0 obs; cimport creates `v1 v2 v3` (str1) with 1 observation containing the header text. Harness (A06a-8): `\nname,value\na,1\nb,2\n` with `varnames(1)` → names lost, 2 obs of `v1 v2`; default → 3 string obs (native: `name value`, 2 obs).
- **Location:** `src/cimport/cimport_impl.c:226-239`, `:595`, `:984`, `:1573-1579`.
- **Proposed fix:** Choose the header as the first non-empty row (including `\r`-only rows) and record its file offset; treat a single row as a header when it is the only row and native would.

<a id="imp-10"></a>
#### IMP-10 [P2] — CR-only (classic Mac) line endings are not recognized
*Component: cimport (delimited text)*

- **Evidence:** Harness (A06a-9): `id,name,value\r1,Alpha,100\r2,Beta,200\r` → 1 obs × 7 vars (native 2 × 3).
- **Location:** `src/cimport/cimport_parse.c:432-531`, `:167-233`.
- **Proposed fix:** Treat a bare `\r` not followed by `\n` as a row terminator outside quotes (parser, row finders, boundary search), or normalize lone CRs in a pre-pass.

<a id="imp-11"></a>
#### IMP-11 [P2] — Fixed-size metadata buffers: `numericcols()`/`stringcols()` silently truncated at 4 KB; wide files (≥ ~2,000–5,000 columns) fail r(198); header lines longer than one chunk are parsed as data
*Component: cimport (delimited text)*

- **Evidence:** Harness (A06a-10/11/12): `stringcols(1/1100)` (4,392 chars) → columns 1041–1100 import as numeric (leading zeros lost); 6,000 columns → only 4,682 names transferred then r(198); a 260 KB header on a 2 MB file → all columns string and load r(459). Possible 1-byte stack overflow (`sizeof(buf)` passed to `SF_macro_use`).
- **Location:** `src/cimport/cimport_impl.c:1415-1463` (`char buf[4096]`), `:1473` (`macro_val[65536]`), `:1519-1562`, `:527-534`, `:595`, `:1109-1111`, `:1760-1771`.
- **Proposed fix:** Length-checked dynamic transfer of lists and per-column metadata (as `cexport_parse.c:171-186` does); locate the header end before chunking; pass `sizeof-1` to `SF_macro_use`; assert row counts before storing.

<a id="imp-12"></a>
#### IMP-12 [P2] — A failed scan leaks the `CIMPORT_NUMCOLS`/`CIMPORT_STRCOLS` globals into later imports
*Component: cimport (delimited text)*

- **Evidence:** Code reading (A06a-13).
- **Location:** `build/cimport.ado:327-332` (set), `:355-359` (early exit on scan error), `:376` (dropped only on success); C re-reads them every scan (`cimport_impl.c:1479-1480`).
- **Trigger:** `cimport delimited using f, stringcols(2) encoding(bogus) clear` (fails) → next plain `cimport delimited using f` imports column 2 as string.
- **Proposed fix:** Always assign both (empty when unused) and drop them immediately after the plugin call, or pass as locals.

<a id="imp-13"></a>
#### IMP-13 [P2] — Encoding conversion allocates 4× the file size, the zero-copy ASCII path is too narrow (numeric-only files are detected as ISO-8859-9), and on Windows every conversion of a > 2 GiB file fails; `encoding(latin1)` is rejected on Windows
*Component: cimport (delimited text)*

- **Evidence:** Stata logs show numeric files auto-detected as ISO-8859-9/-2; harness and code reading (A06a-14/15); Windows behavior not executed.
- **Location:** `src/cimport/cimport_impl.c:833-867`; `src/cimport/cimport_encoding.c:611-626` (`capacity=size*4+32`; only latin2/latin9 canonicalized); `src/io/cio_iconv_win.h:13-58,83-85`.
- **Impact:** Large imports fail with r(920)/r(198) or use far more memory than needed; a common Stata idiom fails on Windows only.
- **Proposed fix:** Skip conversion for pure-ASCII data with ASCII-compatible charsets; size output per charset bound or grow in blocks; feed the Windows shim ≤ INT_MAX blocks; canonicalize Java/ICU aliases (latin1, l1, cp819, …).

<a id="imp-14"></a>
#### IMP-14 [P2] — The automatic header decision inspects only the first 100 data rows; native decides from full-file column types
*Component: cimport (delimited text)*

- **Evidence:** PLAUSIBLE (A06a-16; native behavior from bytecode).
- **Location:** `src/cimport/cimport_impl.c:246-274`.
- **Trigger:** `name,score` + 150 numeric scores + a late `n/a` → native: no header (152 obs); cimport: header (151 obs).
- **Proposed fix:** Decide the header after full inference using per-chunk stats with the first row tracked separately.

<a id="exp-3"></a>
#### EXP-3 [P2] — `delimiter(" ")` is lost and the next option is swallowed: output is comma-separated and `replace`/`novarnames`/`quote` are ignored (SEP24 #6 remnant)
*Component: cexport (delimited + Excel)*

- **Evidence:** Stata-verified: `cexport delimited using f, delimiter(" ") replace` writes `id,s` / `1,a` with commas (native: `id s`). Harness: `"replace noheader"` → delimiter `,`, replace=0.
- **Location:** `build/cexport.ado:101-114`, `:270-271` (delimiter embedded as a bare token); `src/cexport/cexport_parse.c:110-119` (first token is always the delimiter slot).
- **Proposed fix:** Transfer the delimiter as a length-checked local (like the filename) or as `delim=<ascii code>`; parse options by keyword.

<a id="exp-4"></a>
#### EXP-4 [P2] — Any value-labeled variable whose name has ≥ 24 characters makes `cexport delimited` and `cexport excel` fail with r(198)
*Component: cexport (delimited + Excel)*

- **Evidence:** Stata-verified: `respondent_employment_status` with a value label → "_decoded_respondent_employment_status invalid name", r(198) for both CSV and XLSX; native export works.
- **Location:** `build/cexport.ado:127`, `:485` (`tempvar decoded_`var'` exceeds the 31-character local-name limit).
- **Proposed fix:** Index temporaries by position (as the date path does with `fmtvar_`export_idx'`) or keep a list of tempvar names.

<a id="exp-5"></a>
#### EXP-5 [P2] — Empty XLSX export with two threads writes a corrupt worksheet (bad CRC); empty exports without a header get dimension `A1:A0`
*Component: cexport (delimited + Excel)*

- **Evidence:** Stata-verified: `sysuse auto` + `cexport excel using f.xlsx if price < 0, firstrow(variables) replace threads(2)` → success message; `unzip -t` reports "bad CRC 0c0a30fc (should be 7d33ef06)" for `xl/worksheets/sheet1.xml`; native `import excel` → r(603).
- **Location:** `src/cexport/cexport_xlsx.c:691-710` (a task whose last segment is empty never emits `TDEFL_SYNC_FLUSH`), `:748-759`, `:1256`, `:1212-1215`.
- **Impact:** Unreadable workbook reported as saved; 2-CPU machines hit it without `threads(2)`.
- **Proposed fix:** Flush based on the last non-empty segment (or always flush at the end of each task); skip empty segments when partitioning; omit `<dimension>` when there are no rows.

<a id="exp-6"></a>
#### EXP-6 [P2] — Temp-file publication changes file identity: new outputs are mode 0600, `replace` resets permissions/ACLs, replaces symlinks instead of writing through, and overwrites read-only files; exports without `replace` fail on filesystems without hard links (exFAT, SMB, FUSE)
*Component: cexport (delimited + Excel)*

- **Evidence:** Harness (A07-5/6); the repo's own `temp/` outputs from real Stata runs are `-rw-------` for cexport and `-rw-r--r--` for native. exFAT disk image: new destination without `replace` → r(603) "Operation not supported".
- **Location:** `src/cexport/cexport_parse.c:227-266` (`mkstemp` → 0600, then `rename`/`link`); XLSX uses the same path (`cexport_xlsx.c:1686-1698`).
- **Impact:** Outputs unreadable to collaborators on shared servers; stale symlink targets; protection bypassed; exports to USB/network drives fail.
- **Proposed fix:** `fchmod` the temp to `0666 & ~umask` (or the existing destination's mode/owner); resolve symlinks before creating the sibling; fail with 603 when the destination exists but is not writable; fall back to `RENAME_NOREPLACE`/`renamex_np(RENAME_EXCL)` or an `O_EXCL` claim when `link()` is unsupported.

<a id="exp-7"></a>
#### EXP-7 [P2] — XLSX strings with invalid UTF-8 produce malformed XML (unreadable workbook); control characters are silently deleted; leading/trailing blanks, literal `_xHHHH_`, and CR are not preserved
*Component: cexport (delimited + Excel)*

- **Evidence:** Harness (A07-8/13): `"caf\xe9"` makes `sheet1.xml` fail XML parsing; `"a\x01b\x1fc"` is written as `abc`; no `xml:space="preserve"` (native writes it).
- **Location:** `src/cexport/cexport_xlsx.c:160-330` (`xlsx_escape_xml`), `:1083-1090`, `:1225-1228`, `:139-150`; `cexport_xlsx_edit.inc:532-541`.
- **Impact:** Legacy (pre-Unicode) Stata data export to a workbook Excel reports as corrupt; text silently altered.
- **Proposed fix:** Validate UTF-8 while escaping (replace invalid bytes or transcode from Latin-1, or fail clearly); encode C0 controls, CR and literal `_x` sequences with `_xHHHH_`; add `xml:space="preserve"` when needed.

<a id="exp-8"></a>
#### EXP-8 [P2] — XLSX chunk formatting ignores allocation failures: rows silently missing from a "successful" export
*Component: cexport (delimited + Excel)*

- **Evidence:** Fault injection (A07-9): one failed `realloc` > 1 MiB → rc 0, valid archive with `<dimension ref="A1:A3001">` and 0 of 3,000 data rows.
- **Location:** `src/cexport/cexport_xlsx.c:1049` (`(void*)-1` ignored at `:1318-1327`), `:534-583`, `:1354-1361`, `:693`; `cexport_xlsx_edit.inc:596-603`.
- **Proposed fix:** Add a `failed` flag to the chunk args, check all chunks before building segments, and return 920.

<a id="exp-9"></a>
#### EXP-9 [P2] — Memory blow-ups: decoded label columns are recast to `str2045` (N × 2,045 bytes each), and when any column is strL every `str#` cell gets its own 2,046-byte malloc
*Component: cexport (delimited + Excel)*

- **Evidence:** Reasoning (A07-10); harness (A07-11): 100,000 rows of a 1-character `str5` column → 204.8 MB.
- **Location:** `build/cexport.ado:127-132`, `:485-488` (`recast str2045`); `src/cexport/cexport_xlsx.c:1574-1578` (`malloc(length+1)` with `length = 2045` for non-strL), used by CSV (`cexport_impl.c:778`) and XLSX (`:1666`).
- **Impact:** Large labeled or strL-containing exports fail with r(920) or thrash where native succeeds.
- **Proposed fix:** Drop the recast (let `replace` promote widths); load `str#` via a stack buffer and copy the actual length or use the flat-buffer loader; stream strL rows.

<a id="exp-10"></a>
#### EXP-10 [P2] — Windows: file APIs receive Stata's UTF-8 paths in the ANSI code page; for XLSX the reserved temp and the miniz-written workbook can be different files (empty file published, success reported)
*Component: cexport (delimited + Excel)*

- **Evidence:** PLAUSIBLE (code reading; Windows not executed).
- **Location:** `src/cexport/cexport_parse.c:235-240`, `:253-255` (`_mktemp_s`/`_open`/`MoveFileExA`); `src/cexport/cexport_io.c:173` (`CreateFileA`); `src/cexport/cexport_impl.c:311`; `src/cexport/cexport_xlsx.c:1620`; contrast `src/cimport/miniz/miniz_zip.c:49-60` (UTF-16).
- **Proposed fix:** Convert with `MultiByteToWideChar(CP_UTF8)` and use the `W` APIs everywhere (see IMP-7).

<a id="reg-8"></a>
#### REG-8 [P2] — When every regressor is omitted (`K_keep == 0`) creghdfe fails with undefined-macro errors (r(111) under clustering) and the C path returns the wrong RSS
*Component: creghdfe + shared estimation kernels*

- **Evidence:** Stata-verified: `sysuse auto`, `creghdfe price foreign, absorb(foreign)` → r(133) (reghdfe runs with `foreign` omitted). Code reading and harness (A09-5): C stores rss = total SS (89,951) instead of the within SS (146.7); the validation test "all zeros in X" expects the error.
- **Location:** `build/creghdfe.ado:466`, `:545`; `src/creghdfe/creghdfe_regress.c:1135-1165`.
- **Trigger:** `sysuse auto`, `creghdfe price foreign, absorb(foreign)` (reghdfe runs with `foreign` omitted).
- **Proposed fix:** Let `K_keep == 0` continue through the normal path with an empty X; guard the ado's scalar reads; set F missing; update the validation expectation.

<a id="reg-9"></a>
#### REG-9 [P2] — Absorbed degrees of freedom (`e(df_a)`) diverge from reghdfe's `estimate_dof()` for G ≥ 3 FEs, cluster-nested FEs and the `dof()` options
*Component: creghdfe + shared estimation kernels*

- **Evidence:** Stata-verified: `absorb(school year district)` (school nested in district) → reghdfe `e(df_a)` = 64, se .044362; creghdfe 73, se .044740. Harness transcription of reghdfe 6.13.1's DoF logic (A09-6): `absorb(school year district)` 73 vs 64; `absorb(school year) vce(cluster school) dof(pairwise)` 4 vs 64; `dof(firstpair)` with G ≥ 3 and `dof(none)` also differ.
- **Location:** `src/creghdfe/creghdfe_regress.c:549-673`; `build/creghdfe.ado:74-87`, `:421-448`, `:766-815`; `src/ctools_hdfe_utils.c:421-452` (`ctools_compute_hdfe_dof`, shared with civreghdfe).
- **Impact:** Wrong df_a/df_r and therefore SEs, rmse, r2_a (≈1% in the example, larger with few observations per FE).
- **Proposed fix:** Port `estimate_dof()` exactly into shared code (nested detection only when "clusters" is among the adjustments; M = 1 for each remaining intercept; pairwise/firstpair component maxima), mirror `ParseDOF` in the ado, and derive the display table from the same vector.

<a id="reg-10"></a>
#### REG-10 [P2] — Shared `detect_collinearity()` uses an absolute 1e-14 pivot tolerance: legitimate small-scale regressors are dropped (creghdfe `nostandardize`, cqreg, civreghdfe)
*Component: creghdfe + shared estimation kernels*

- **Evidence:** Harness (A09-8): `x = 1e-9*cos(·)` is kept by default but dropped under `nostandardize` (then fails as REG-8); cqreg expected to drop it too (A12-8).
- **Location:** `src/ctools_ols.c:392`, `:422`; `src/creghdfe/creghdfe_regress.c:1013`, `:1089-1098`; `src/cqreg/cqreg_regress.c:814-846`.
- **Proposed fix:** Make the test purely relative (scale by the original diagonal / correlation form); keep an absolute test only for exactly-zero columns.

<a id="iv-14"></a>
#### IV-14 [P2] — gmm2s/cue with a singular moment covariance (few clusters): gmm2s silently returns 2SLS labelled gmm2s; cue posts V = 0 or ~1e-10 SEs
*Component: civreghdfe*

- **Evidence:** Harness (A10a-9); Stata-verified that ivreghdfe errors r(506) ("estimated covariance matrix of moment conditions not of full rank") on the repro.
- **Location:** `src/civreghdfe/civreghdfe_estimate.c:143-154`, `:1368-1374`, `:401-406`, `:537-582`, `:738-847`, `:1515-1538`; `src/ctools_ols.c:34`.
- **Proposed fix:** Detect rank(S) < K_iv; follow ivreg2 (error or generalized inverse with warning); propagate CUE non-convergence; never post V from an uninitialized Hessian.

<a id="iv-15"></a>
#### IV-15 [P2] — Jacobi eigen and `M^{-1/2}` routines stop after 100 single rotations: inaccurate for K ≥ 10 (LIML λ with ≥ 9 endogenous regressors; KP with ≥ 10 excluded instruments)
*Component: civreghdfe*

- **Evidence:** Harness vs numpy (A10a-10).
- **Location:** `src/civreghdfe/civreghdfe_matrix.c:87-145`, `:189-256`.
- **Proposed fix:** Iterate by full sweeps to a relative off-diagonal tolerance (cap ~100 sweeps) or use a LAPACK-style symmetric eigensolver; report failure instead of returning λ = 1.

<a id="iv-16"></a>
#### IV-16 [P2] — `partial()` bookkeeping: phantom trailing `e(b)`/`e(V)` columns; df_r, V, df_m and rmse ignore the partialled count; `nopartialsmall` works in the wrong direction
*Component: civreghdfe*

- **Evidence:** Harness (A10a-11): df_r 548 vs ivreghdfe 547; `matrix colnames` replicates the last name for missing columns.
- **Location:** `build/civreghdfe.ado:941-948`, `:1102-1104`, `:1196-1199`; `src/civreghdfe/civreghdfe_impl.c:1319`, `:1390`, `:1405`, `:1412`; VCE df in `civreghdfe_vce.c`.
- **Proposed fix:** Size b/V to the post-partial K; include n_partial in the dof (unless `nopartialsmall`); set `e(df_m)` accordingly.

<a id="iv-17"></a>
#### IV-17 [P2] — Display/option handling changes stored results: `nofooter` suppresses posting of every diagnostic `e()`; `rf` leaves creghdfe's reduced-form results in `e()` (and builds an invalid weight clause `[aweightw]`); `liml fuller(#)` silently runs plain LIML
*Component: civreghdfe*

- **Evidence:** Stata-verified: with `nofooter`, `e(idstat)`, `e(widstat)`, `e(sargan)`, `e(cstat)` are all missing; `liml fuller(1)` gives b = 1.0583843 (= plain LIML) vs ivreghdfe Fuller 1.058794. Code reading (A10a-12/13/14, A10b-13/17).
- **Location:** `build/civreghdfe.ado:1485-1865` (all `ereturn scalar` diagnostics inside `if "`footer'" == ""`), `:1584-1655` (rf), `:561-597` (liml checked before fuller).
- **Impact:** Scripts using `nofooter` (typical for tables) silently lose J/weak-ID/underID/C statistics; after `rf`, `e(cmd)`, `e(b)`, `e(V)`, `e(sample)` belong to the reduced form; Fuller requests get LIML.
- **Proposed fix:** Post diagnostics unconditionally (only the display is conditional); wrap the RF regression in `_estimates hold/unhold` and build `[`weight'=`exp']`; treat `fuller != 0` as Fuller (ivreg2 `ivparse` rules).

<a id="iv-18"></a>
#### IV-18 [P2] — HAC test statistics are O(N²): minutes to hours for time series or badly grouped panels, with no way to interrupt
*Component: civreghdfe*

- **Evidence:** Harness timings (A10b-11): KP LM+Wald 0.41 s at N=10k, 7.15 s at 40k; extrapolates to ~20 min at 400k.
- **Location:** `src/civreghdfe/civreghdfe_tests.c:790-807`, `:1181-1198`, `:2056-2072` (and panel variants `:757-783`, `:1146-1172`, `:2024-2049`), `civreghdfe_vce.c:840-865`.
- **Proposed fix:** Loop over lags τ = 0..bw and only over partners at lag τ (as `ivvce_compute_full` does); share one HAC-S builder between VCE and tests.

<a id="iv-19"></a>
#### IV-19 [P2] — J uses the Bartlett kernel regardless of `kernel()` and the HC S under two-way clustering; `vce(cluster tvar) bw()` and `vce(cluster a b) bw()` silently ignore the kernel
*Component: civreghdfe*

- **Evidence:** Harness (A10b-10): J = 6.21087732 for all kernels while V changes; two-way J equals the HC-robust J. Code reading (A10b-18).
- **Location:** `src/civreghdfe/civreghdfe_tests.c:1888-1889` (`test_kernel = 1`), `:1943-2077`; `src/civreghdfe/civreghdfe_estimate.c:1539-1557`, `:1775-1783`; `src/civreghdfe/civreghdfe_vce.c:746-778`.
- **Proposed fix:** Use the user's kernel for J (only ranktest id/wid stats default to Bartlett); pass cluster2 ids and use the CGM S for two-way J; route `cluster(tvar) bw()` to the DK path or reject cluster+bw until implemented.

<a id="iv-20"></a>
#### IV-20 [P2] — Documented DWH test `e(endog_chi2)`/`e(endog_p)` uses `(v'v)⁻¹` instead of `[(X_a'X_a)⁻¹]_vv` and ignores weights/VCE (≈10× inflated)
*Component: civreghdfe*

- **Evidence:** Harness (A10b-12): 85.39 vs correct 7.27.
- **Location:** `src/civreghdfe/civreghdfe_tests.c:2307-2376`; posted at `build/civreghdfe.ado:1700-1705`; documented `build/civreghdfe.sthlp:472-475`.
- **Proposed fix:** Use the correct augmented-regression Wald with weights (and sandwich under robust/cluster), or drop it in favor of `endogtest()`.

<a id="iv-21"></a>
#### IV-21 [P2] — `dofminus()`/`sdofminus()` never reach V or the tests; `small` is a no-op; without `absorb()` small-sample scaling is always applied (diverges from ivreg2's large-sample default); no-absorb F/df_m include `_cons`
*Component: civreghdfe*

- **Evidence:** Stata-verified: `dofminus(50)` → civreghdfe V differs from ivreghdfe by mreldif 1.4e-5 (identical, 3.9e-18, without the option), i.e. the option is ignored in V. Code reading vs ivreghdfe/ivreg2 (A10b-14, A10a-19).
- **Location:** `src/civreghdfe/civreghdfe_impl.c:1319`, `:1390-1406`; `src/civreghdfe/civreghdfe_vce.c:612-622`, `:764`, `:929`, `:945`; `build/civreghdfe.ado:31`, `:350-361`, `:1201-1207`, `:1320`, `:1340-1356`.
- **Proposed fix:** Pass `dofminus`, `sdofminus` (+n_partial) and a `small` flag to the estimator/VCE/tests; with G = 0 and no `small`, use dof_adj = 1 and σ² = rss/(N − dofminus); exclude `_cons` from F/df_m.

<a id="iv-22"></a>
#### IV-22 [P2] — IV VCE and test failures are still silent (SEP24 #11 only partly propagated): two-way and Kiefer VCE are `void`; HAC buffer failure returns rc 0 with V = 0; thread-local `calloc` failures silently drop terms (V can go negative); the dense `(G1+1)(G2+1)` pair table can need tens of GB; failed diagnostics are posted as 0 (p = 1)
*Component: civreghdfe*

- **Evidence:** Fault-injection harness (A10b-15): V = 0 or V = −0.00939 posted with rc 0; code reading (A10a-20, A10b-22).
- **Location:** `src/civreghdfe/civreghdfe_vce.c:43-66`, `:93-114`, `:155-176`, `:196-327`, `:348-552`, `:660-681`, `:914-921`; `src/civreghdfe/civreghdfe_tests.c:50-70`, `:363-366`, `:611-687`, `:1891-1892`, `:2440-2441`; `src/civreghdfe/civreghdfe_estimate.c:1549-1587`; `build/civreghdfe.ado:1487-1494`, `:1684-1695` (`widstat = cd_f` when `kp_f <= 0`).
- **Proposed fix:** Return status from all VCE builders and propagate 920; allocate thread-local buffers before parallel regions (or stream by sorted cluster); build intersection ids by sorting pairs; initialize diagnostic outputs to missing and post missing with a warning on failure.

<a id="iv-23"></a>
#### IV-23 [P2] — Integer overflow and large-data hazards: `int` N·K allocation sizes can wrap (> ~537M obs with K = 8); fweights summing to > 2^31−1 overflow `(ST_int)N_eff` (negative V); singleton removal returning −1 on allocation failure makes `N = N_valid + 1` and overruns buffers
*Component: civreghdfe*

- **Evidence:** Code reading (A10b-16, A10b-19, A09-10).
- **Location:** `src/civreghdfe/civreghdfe_vce.c:218`, `:625`, `:790`; `src/civreghdfe/civreghdfe_tests.c:400-402`, `:1296-1297`, `:1514`, `:2256`, `:2275`, `:2976-2977`; `src/civreghdfe/civreghdfe_impl.c:373-378`, `:459-464`, `:1311`.
- **Proposed fix:** Use `ctools_safe_calloc3`/`size_t` arithmetic; carry `N_eff` as double; treat a negative singleton return as 920.

<a id="iv-24"></a>
#### IV-24 [P2] — Two-way clustering: an FE nested in the *second* cluster variable is not treated as redundant (canonical `absorb(firm year) vce(cluster firm year)` design)
*Component: civreghdfe*

- **Evidence:** Source comparison with reghdfe.mata:4311-4341 (A10a-16).
- **Location:** `src/civreghdfe/civreghdfe_impl.c:962-982`; `build/civreghdfe.ado:1034-1039`.
- **Proposed fix:** Run the nesting check against every cluster dimension (and intersections) and export per-FE flags for vce_type 3.

<a id="qreg-4"></a>
#### QREG-4 [P2] — The Frisch–Newton corrector omits the /x and /s scaling of the second-order Mehrotra term: 2–5× more iterations (r(430) at the default limit on tail quantiles) and a trigger for QREG-1
*Component: cqreg*

- **Evidence:** Harness with a patched copy (A12-4): 1,329 vs 249 iterations at q = .99, N = 1M; identical solutions after the fix.
- **Location:** `src/cqreg/cqreg_fn.c:179-180`, `:194-195`, `:300-301` (`corrector_loop1_{neon,scalar}`).
- **Proposed fix:** `dxdz[i] = dx[i]*dz[i]/x_p[i]`, `dsdw[i] = ds[i]*dw[i]/s[i]`.

<a id="qreg-5"></a>
#### QREG-5 [P2] — `vce(cluster a b)` silently clusters on the a×b intersection
*Component: cqreg*

- **Evidence:** Stata-verified: `webuse nlswork`, `cqreg ln_wage age tenure, vce(cluster idcode year)` → `e(N_clust) = 28,101` (idcode-year cells), rc 0.
- **Location:** `build/cqreg.ado:107-121`.
- **Proposed fix:** Accept a single cluster variable (`syntax varname`) and error otherwise, or implement multi-way clustering.

<a id="qreg-6"></a>
#### QREG-6 [P2] — If every regressor is dropped (`K_keep = 0`), `e(V)` is all zeros and SE(_cons) = 0
*Component: cqreg*

- **Evidence:** Stata-verified: constant `x` → qreg `V[_cons,_cons] = .0168`; cqreg `0`.
- **Location:** `build/cqreg.ado:305-329` (the `_cons` variance is copied only inside `if K_keep > 0`).
- **Proposed fix:** Move the `_cons` variance assignment out of the conditional.

<a id="qreg-7"></a>
#### QREG-7 [P2] — Absolute tolerances break scale-equivariance: auxiliary gap tolerance 1e-6 (SEs 22× too small for y ~ 1e-8), robust regularized inverse gives ~1e10 SEs, tiny-scale regressors dropped as collinear, main-solve/pre-check tolerances give wrong β or spurious r(198) at y ~ 1e-12
*Component: cqreg*

- **Evidence:** Harness (A12-7, -8, -14).
- **Location:** `src/cqreg/cqreg_regress.c:289-291`, `:481-483`, `:784`, `:814-846`; `src/cqreg/cqreg_vce.c:695-707`; `src/cqreg/cqreg_fn.c:1010-1021`, `:1413`; `src/ctools_ols.c:392-423` (see REG-10).
- **Proposed fix:** Use relative gap/objective criteria; sweep/generalized inverse instead of +1e-10 regularization; correlation-form collinearity test; error when every f_i is 0.

<a id="ppml-1"></a>
#### PPML-1 [P2] — Default `use_exact_partial(1)` (ppmlhdfe defaults to 0) plus an early-stopping inner CG makes IRLS oscillate: r(430) or ~26× more iterations on poorly connected two-way FE data
*Component: cpplmhdfe*

- **Evidence:** Harness (A11-3); reference side plausible.
- **Location:** `build/cpplmhdfe.ado:35`; `src/cpplmhdfe/cpplmhdfe_irls.c:1378-1396`; `src/creghdfe/creghdfe_solver.c:366-373`.
- **Proposed fix:** Default `use_exact_partial` to 0; make the CG stop more robust (several consecutive small residual changes or a true-residual check); detect IRLS oscillation and fall back to warm-started partialling.

<a id="ppml-2"></a>
#### PPML-2 [P2] — False IRLS convergence for large-scale outcomes (0.1 floor on the standardized deviance): wrong `e(ll)`, `e(deviance)`, `e(r2_p)` and slightly wrong V with `e(converged) = 1`
*Component: cpplmhdfe*

- **Evidence:** Harness vs exact MLE (A11-4): `e(deviance)` 17,126 vs true 393.5, `e(ll)` off by 8,000, SEs +0.4%. ppmlhdfe uses the same formula, so the reference may share the flaw.
- **Location:** `src/cpplmhdfe/cpplmhdfe_irls.c:1544-1546`.
- **Proposed fix:** Use a scale-free denominator (e.g., relative to Σw·y or 0.1/stdev_y) and additionally require a small max |Δη| before accepting convergence.

<a id="ppml-3"></a>
#### PPML-3 [P2] — `savefe`/`absorb(name=var)` uses unaccelerated sweeps capped at `iterate()`: a successful fit turns into r(430) with no message
*Component: cpplmhdfe*

- **Evidence:** Harness (A11-5).
- **Location:** `src/cpplmhdfe/cpplmhdfe_irls.c:184-230`, `:2201`.
- **Proposed fix:** Recover FEs with the same CG/LSMR used for partialling (D'WDα = D'Wd), emit a message on failure, and consider posting results with a warning when only FE recovery fails.

<a id="ppml-4"></a>
#### PPML-4 [P2] — Factor variables are expanded before the FE/singleton screen and cluster markout, so base/omitted levels (and `_b[#.var]`) can differ from ppmlhdfe
*Component: cpplmhdfe*

- **Evidence:** PLAUSIBLE (source reasoning, A11-6).
- **Location:** `build/cpplmhdfe.ado:366-395`, `:431-439`.
- **Proposed fix:** Mark out cluster parts before `fvexpand`; re-expand on the final sample after the first plugin pass and rerun if the level set changed.

<a id="bin-7"></a>
#### BIN-7 [P2] — By-group bookkeeping: the plugin still derives the group count as max − min + 1 (SEP20 A07 fix incomplete), so emptying the lowest group silently drops the highest one; by() bins are computed within groups while binscatter pools bins
*Component: cbinscatter*

- **Evidence:** Harness (A13-8/9): with group 1 emptied by BIN-5, `num_groups = 2` and group 3's 300 valid observations are discarded (still counted in `e(N)`).
- **Location:** `src/cbinscatter/cbinscatter_impl.c:631-681`, `:676-754`, `:723`; `src/cbinscatter/cbinscatter_fit.c:421-461`; `build/cbinscatter.ado:184-192`, `:463`; `build/cbinscatter.sthlp:113-114`.
- **Proposed fix:** Pass G from the ado and index groups 1..G regardless of emptiness; for method(classic) compute cutpoints on the pooled sample (binscatter behavior) and keep per-group bins for method(binsreg).

<a id="bin-8"></a>
#### BIN-8 [P2] — binscatter's automatic discrete mode (#unique(x) ≤ nquantiles) is missing (Likert 1–7 → 6 bins); `discrete` with > 500 distinct values fails with a plugin error
*Component: cbinscatter*

- **Evidence:** Harness (A13-10/13).
- **Location:** `src/cbinscatter/cbinscatter_impl.c:458-471`, `:742-754`; `src/cbinscatter/cbinscatter_bins.c:238-600`; `build/cbinscatter.ado:244-248` (`max_bins = min(nobs, 500)`).
- **Proposed fix:** Count distinct x and switch to discrete when ≤ nquantiles; size `e(bindata)` from a pre-pass count and enforce binscatter's N/2 rule.

<a id="bin-9"></a>
#### BIN-9 [P2] — HDFE residualization that never converges is silently accepted (10,000 sweeps; message only under `verbose`; batch solver checks convergence on y only) — the cbinscatter counterpart of SEP24 #4
*Component: cbinscatter*

- **Evidence:** Harness (A13-11): bin means off by 0.023–0.027 on a y-range of ~3 for a weakly connected two-way FE design.
- **Location:** `src/cbinscatter/cbinscatter_resid.c:401-457`, `:624-662`; `build/cbinscatter.ado:216-217`.
- **Proposed fix:** Return r(430) (or warn and post `e(converged) = 0`) on exhaustion; test all batch columns; consider the shared accelerated solver.

<a id="bin-10"></a>
#### BIN-10 [P2] — Polynomial fits use uncentered raw moments with plain Cholesky: `qfit` with binary x fails the whole command (r(199)); fits on dates/datetimes are visibly wrong; constant x gives a garbage slope
*Component: cbinscatter*

- **Evidence:** Harness (A13-12): qfit on %td dates within a month → quadratic coefficient 0.0253 vs `reg` 0.0100.
- **Location:** `src/cbinscatter/cbinscatter_fit.c:84-95`, `:127-200`, `:240-331`, `:402-414`; `src/cbinscatter/cbinscatter_impl.c:852-853`.
- **Proposed fix:** Center/scale x before forming moments, solve with the collinearity-aware solver, zero omitted terms, transform back; never abort the whole command for a fit problem.

<a id="bin-11"></a>
#### BIN-11 [P2] — `method(binsreg)` diverges from binsreg with `absorb()` (singleton FE groups kept → every dot shifted) and with factor-variable controls (evaluated at their means instead of at zero)
*Component: cbinscatter*

- **Evidence:** PLAUSIBLE/source comparison (A13-14/15); validation compares binsreg output at 0.5 significant figures (see BLD-3), which hides these shifts.
- **Location:** `src/cbinscatter/cbinscatter_binsreg.c:183-216`, `:261-273`, `:297-310`, `:425-440`; `build/cbinscatter.ado:174-177`.
- **Proposed fix:** Drop singletons iteratively before computing ȳ/bin shares/control means; flag fvrevar indicator columns and evaluate them at 0 (binsreg rule), or document.

<a id="xio-7"></a>
#### XIO-7 [P2] — Unparseable numeric `<v>` text (including ISO dates in `t="d"` cells) becomes missing, and in string columns becomes the literal text "8.98846567431158e+307"
*Component: cimport/cexport Excel and statistical formats*

- **Evidence:** Stata-verified: native imports `t="d"` cells as 13jul1905 etc.; cimport gives missing, and with `allstring` the text "8.98846567431158e+307".
- **Location:** `src/cimport/cimport_xlsx.c:952-973`, `:1017-1029`; `src/io/cio_xlsx.inc:185-201`.
- **Proposed fix:** Parse ISO-8601 `t="d"` to serials (respecting date1904); store unparseable numeric text as a string cell; never format `SV_missval` as text.

<a id="xio-8"></a>
#### XIO-8 [P2] — SPSS/SAS-catalogue value labels on non-integer keys silently overwrite the integer key's label (1.5 "one and a half" replaces 1 "one")
*Component: cimport/cexport Excel and statistical formats*

- **Evidence:** Stata-verified: native `import spss` skips 1.5/2.5 with notes and keeps "one"/"two"; `cimport spss` labels 1 as "one and a half" and 2 as "two and a half".
- **Location:** `src/io/cio.c:358-380`, `:853-873`; `build/_cio_import.ado:336-344`.
- **Proposed fix:** Skip (with native's note) any key that is not an integer in Stata's label range or is system missing.

<a id="xio-9"></a>
#### XIO-9 [P2] — Left-justified date formats (`%-td`, `%-tc`) are not recognized as dates on export: SPSS/SAS XPORT/dBase files contain raw day or millisecond counts
*Component: cimport/cexport Excel and statistical formats*

- **Evidence:** Stata-verified: a `%-td` date re-imported from native SPSS output shows 01-Jan-2020; from cexport output shows 21915 (same for sasxport8/sasxport5).
- **Location:** `src/io/cio.c:1016`, `:1298-1301`; the ado datafmt regex.
- **Proposed fix:** Normalize formats (`%-` → `%`) with one helper that parses `%[-]t[dcC]`/`%[-]d` and use it everywhere.

<a id="xio-10"></a>
#### XIO-10 [P2] — SAS XPORT naming: `cexport sasxport8` fails with r(610) when the file stem is not a SAS name (`sales-2024`, `2024sales`, `my data`); `cexport sasxport5` writes value-label names through ReadStat's format grammar (names ending in digits split, reserved names like `month` written, `b_` fails); `cimport sasxport5` turns format names with widths into label names
*Component: cimport/cexport Excel and statistical formats*

- **Evidence:** Stata-verified: native `export sasxport8` accepts all three stems (rc 0), cexport r(610); native `export sasxport5` rejects labels `q1`, `labels12`, `month` (r(100)) while cexport writes them (rc 0) and fails on `b_` (r(610)). Harness for the import side (A08-5).
- **Location:** `build/_cio_export.ado:54`; `src/io/cio.c:288-297`, `:909-968`, `:1220-1226`, `:1329-1333`; `build/_cio_import.ado:274-285`.
- **Proposed fix:** Sanitize member names (map invalid characters to `_`, leading letter, truncate to 8/32, fall back to `DATASET`); write NFORM directly and enforce native label-name rules; build import label names from the NFORM name without width.

<a id="xio-11"></a>
#### XIO-11 [P2] — Smaller XLSX-import defects
*Component: cimport/cexport Excel and statistical formats*

- **Evidence:** Stata-verified for sheet names, cellrange validation, header collision and `.xlsm`; harness for the rest (A06b-7…25).
- **Items:**
  - `sheet("P&L")` → r(601) (names not entity-decoded; `describe` shows "P&amp;L"); names > 63 bytes truncated; sheets after the 256th dropped (`cimport_xlsx.c:263`, `:271-274`; `cimport_xlsx.h:23-24`).
  - Invalid `cellrange()` accepted silently: `foo`, `A1:B`, `A0:B2`, `A1:B2junk` import data; `B5:A1`, `AAAA1`, `Sheet1!A1:B2` give an empty dataset; `A1:A2000000` creates 2,000,000 rows (native r(198) for all) (`cimport_xlsx.c:1973-2017`; `cio_xlsx.inc:82-85`).
  - A header equal to another column's fallback letter name makes the import fail with r(110) (native succeeds); `case()` folds ASCII only (`_cio_import.ado:251-271`).
  - `.xlsm` workbooks rejected with r(198) (native imports) (`cimport.ado:711-716`).
  - Trailing blank formatted columns are dropped (native keeps 4 columns, cimport 2); a leading empty `<row>` makes the header row data under `firstrow`; with `cellrange()` a self-closing `<row/>` skips the next row (`cimport_xlsx.c:1063`, `:1203-1210`, `:1510-1520`; `cio_xlsx.inc:103-106`).
  - Rich-text inline strings keep only one run; datetime text can read "24:00:00"; `_xHHHH_` escapes kept literally; booleans/errors differ from native (`cimport_xlsx.c:410`, `:981-993`, `:1034-1047`; `src/io/cio.c:629-635`).
  - `describe` returns garbage ranges ("B2147483647:A1") and one unresolvable sheet aborts it (`cio_excel.inc:2-15`, `:71-94`).
  - Serial-60 policy disagrees across code: the live readers map serial 60 to 28feb1900 (matching native — Stata-verified), while the dead helper returns missing and `validation/validate_sep24.do:247` asserts missing (the test should now fail).
  - Robustness: 1-byte heap over-read at the end of a stored worksheet part; signed overflow in cell-reference/entity parsing (a single row-less `r="B"` empties the import); parallel scanner restarts the implicit row counter per chunk; CDATA lost; attribute values truncated at 255 bytes / 16 attributes; OOM in callbacks reported as success; no CRC check on the fast path; a huge `<dimension>` pre-allocates rows×cols (`cimport_xlsx.c:744-771`, `:934`, `:1166-1210`, `:1190`, `:1496-1498`, `:1561-1563`; `cimport_xlsx_xml.c:385-399`, `:447-489`; `cimport_xlsx_zip.c:140-153`, `:221-233`).
- **Proposed fix:** Entity-decode attributes and make sheet lists dynamic; validate `cellrange()` as the XLS reader does; unique fallback names; accept `.xlsm` by content signature; derive the used range from stored cells; concatenate rich-text runs; carry 24:00:00 into the next day; one shared serial-date function and a decided serial-60 policy (align test and help); bound all parsers and NUL-terminate buffers; verify CRCs.

<a id="psm-2"></a>
#### PSM-2 [P2] — cpsmatch tie-breaking and `noreplacement` algorithm differ from psmatch2 (different matched controls and ATT)
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Evidence:** Stata-verified: tie example → psmatch2 ATT 50, cpsmatch −50; `noreplacement` example → psmatch2 −30, cpsmatch 20. Harness (A15-6/7): left-side ties → psmatch2 ATT 50, cpsmatch −50; noreplacement example → psmatch2 −30, cpsmatch 20 (cpsmatch uses nearest-available greedy matching; psmatch2 scans forward only).
- **Location:** `src/cpsmatch/cpsmatch_impl.c:222-232`, `:470-499`, `:1009-1082`; help `build/cpsmatch.sthlp`.
- **Proposed fix:** Pick the first element of an equal-pscore run (psmatch2 order); implement psmatch2's forward scan for noreplacement or document the deliberate algorithmic difference prominently.

<a id="smp-2"></a>
#### SMP-2 [P2] — The shared counting sort (by()/strata()/cluster() grouping) corrupts groups for integer keys beyond ±2^63 (see SORT-2); csample/cwinsor allocate threads × N work buffers even when one thread does the work
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Evidence:** Harness (A15-9): `csample, count(1) by(g)` with g ∈ {−1e19, 1e19, 5, 7} keeps 3 rows instead of 4; arithmetic for buffers (A15-8, PLAUSIBLE): 300M obs × 32 threads = 76.8 GB single request.
- **Location:** `src/ctools_sort_counting.c:110-152`, `:319-320`; `src/csample/csample_impl.c:517-534`; `src/cwinsor/cwinsor_impl.c:499-537`.
- **Proposed fix:** Fix SORT-2; size per-thread buffers to `min(threads, ngroups)` (or allocate lazily), and sample without an O(group) index array (Floyd's algorithm).

<a id="win-3"></a>
#### WIN-3 [P2] — cwinsor's default in-place mode stores fractional bounds into integer variables without promotion (bounds 1.5/99.5 become 1/99)
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Evidence:** Stata-verified: `gen int x = _n` (1..100), `cwinsor x, cuts(1 99)` → rows 1 and 100 become 1 and 99 and `x` stays `int`; `winsor2 y, cuts(1 99) replace` → 1.5 and 99.5 with `y` promoted to float.
- **Location:** `build/cwinsor.ado:103-121`, `:152-153`; `src/cwinsor/cwinsor_impl.c:635-644`.
- **Proposed fix:** Recast integer replace targets to double (or float/double by range) before the plugin call when non-integer bounds are possible.

<a id="rng-3"></a>
#### RNG-3 [P2] — crangestat `first`/`last` return the first/last *non-missing* value (identical to `firstnm`/`lastnm`), unlike rangestat and crangestat's own help
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Stata-verified: with `x[2]` missing and `interval(t -1 1)`, obs 1 `last` = 10 (rangestat `.`) and obs 3 `first` = 30 (rangestat `.`); this is also the `(first)` mismatch seen in the randomized differential run. Harness (A14-3) with a rangestat reference.
- **Location:** `src/crangestat/crangestat_impl.c:1023-1039` (SIMD path), `:1042-1054` (scalar path), `:1084-1090`; `build/crangestat.sthlp:128-129`.
- **Problem:** rangestat returns the source value of the first/last observation of the window (by key, then original order) even when it is missing, and with `excludeself` blanks the self row but keeps it; crangestat skips missing values and self.
- **Trigger:** `set obs 10`, `gen t = _n`, `gen x = 10*_n`, `replace x = . in 2`, `crangestat (first) cf=x (last) cl=x, interval(t -1 1)` → obs 1 `last` = 10 (rangestat `.`), obs 3 `first` = 30 (rangestat `.`); with `excludeself` crangestat returns the next row's value where rangestat returns missing.
- **Impact:** Silent wrong values whenever a window edge holds a missing value or self sits at the edge under `excludeself`.
- **Proposed fix:** Return `data[win_start]`/`data[win_end-1]` for FIRST/LAST (missing if that row is self under `excludeself`); keep the non-missing scans for FIRSTNM/LASTNM; add tests with missing edge values.

<a id="rng-4"></a>
#### RNG-4 [P2] — crangestat default result names are `stat_var` (e.g., `mean_x`); rangestat creates `var_stat` (`x_mean`)
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Stata-verified: `crangestat (mean) v` creates `mean_v`, `rangestat (mean) v` creates `v_mean`. Source comparison (`build/crangestat.ado:196-200` vs `rangestat.ado:146`); the maintainer's own notes record rangestat's convention as `varname_stat`; validation always passes explicit names.
- **Location:** `build/crangestat.ado:196-200`.
- **Impact:** Ported rangestat scripts fail later with r(111), or silently pick up an unrelated same-named variable.
- **Proposed fix:** `local result_var "`source_var'_`stat_name'"` + `confirm name`; add an auto-naming test.

<a id="rng-5"></a>
#### RNG-5 [P2] — Reversed bounds (`low > high`) give plausible-looking garbage (count ≈ 1.8e19, sd ≈ 9.5e153) in groups with ≥ 64 observations instead of missing
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Stata-verified: 64 rows with `interval(t 1 -1)` → every row count 1.84e19 and sd 9.5e153. Harness (A14-5): 63-row group → all missing (correct); 64-row group → 256/256 garbage cells.
- **Location:** `src/crangestat/crangestat_impl.c:366-376`, `:760-801`; windows from `:1865-1869`, `:2118-2119`, `:2251-2255`; `build/crangestat.ado:68-96`.
- **Proposed fix:** Treat `win_end <= win_start` as an empty window before any prefix/sparse query; optionally warn when `low > high`.

<a id="rng-6"></a>
#### RNG-6 [P2] — crangestat skewness/kurtosis: wrong small-n cut-offs (missing for n = 2/3 where rangestat reports values) and constant windows produce fabricated skewness ±1 / kurtosis 1 instead of missing
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Stata-verified: two-value windows → rangestat skewness 0/kurtosis 1, crangestat missing; constant 0.1 windows (`interval(t -5 5)`, 40 rows) → crangestat skewness 1/kurtosis 1 at obs 1, 2, 39, 40, while rangestat gives missing at obs 2 and 39 (and its own rounding artifact −1/1 at obs 1 and 40). Harness (A14-6): the naive mean is off by one ulp, so m2 ≈ 1e-36 > 0.
- **Location:** `src/crangestat/crangestat_impl.c:1092-1132`.
- **Proposed fix:** Drop the n ≥ 3/4 thresholds; compute the mean as `origin + mean(x − origin)` so constant windows give exactly zero deviations; return missing only when the central second moment is zero.

<a id="rng-7"></a>
#### RNG-7 [P2] — crangestat precision: prefix-path variance/SD keep only ~7 significant digits on long trending series (SEP20 A10 fix incomplete), skewness/kurtosis lose digits with large offsets, and an outlier in a group's first row degrades rolling means/sums for the whole group
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Harness against compensated references (A14-6): N = 200,000 trending series, `interval(t -50 49)` → variance minimum 6.82 significant digits (176 observations below 7); N = 1M → 6.70; skewness 5.5 digits at offset 1e9. Draft harness (A14 `rs_origin.c`): first value 1e12 in a 100k-row group → rolling mean 6.2373 vs 5.97797.
- **Location:** `src/crangestat/crangestat_impl.c:225-363` (uncompensated prefix sums centered on the group's first value; chunked accumulation makes rounding thread-dependent), `:789-800` (fallback only when `ss <= sqrt(DBL_EPSILON)*scale`), `:1094-1131`.
- **Impact:** Fails the project's own 7-significant-figure standard; last digits vary with `threads()`; sentinel/outlier first rows silently corrupt rolling means.
- **Proposed fix:** Double-double (two-sum) prefix sums or block-local origins combined with Chan/Welford formulas; a rigorous error-bound fallback; deterministic logical chunks; origin-shifted means for higher moments.

<a id="rng-8"></a>
#### RNG-8 [P2] — crangestat percentiles and `iqr` use linear (type-7) interpolation, not Stata's `_pctile`/`summarize` definition, and the help does not say so
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Stata-verified: x = 1,2,3,4 → crangestat p25 1.75, p75 3.25, iqr 1.5; `_pctile` gives 1.5 and 3.5 (iqr 2). Harness (A14-10): p1 1.03 vs 1. (The validation suite cannot catch this; see BLD-3.)
- **Location:** `src/crangestat/crangestat_impl.c:170-193`; `build/crangestat.sthlp:131-139`.
- **Proposed fix:** Implement Stata's definition (average of x(i) and x(i+1) when n·p/100 is an integer i, else x(ceil)), or document the rule prominently.

<a id="bld-3"></a>
#### BLD-3 [P2] — Validation helpers let real divergences pass: `benchmark_rangestat` records a PASS whenever `rangestat` errors (9 percentile/IQR gate tests can only fail if all output is missing); `benchmark_ivreghdfe` never compares Hansen J under robust/cluster/HAC (ivreghdfe posts `e(j)`, not `e(sargan)`) nor any `ffirst` partial R²; cbinscatter `method(binsreg)` is compared at 0.5 significant figures (≈32% relative error)
*Component: Build, packaging, CI and validation*

- **Evidence:** Source (A16-4/5/6); these gaps are why IV-7, IV-8, BIN-4 and BIN-11 pass the current suite.
- **Location:** `validation/benchmark_helpers.do:1306`, `:1321-1324`, `:1355`, `:1370-1373`, `:1426-1442`, `:1580-1597`, `:2237-2269`; `validation/validate_crangestat.do:164,198,205,267,274,1301-1322`; `validation/validate_civreghdfe.do:542-549`, `:1811`; `validation/validate_cbinscatter.do:460-485`, `:696`.
- **Proposed fix:** Compute independent references (or record explicit, counted exclusions) instead of passing; compare ivreghdfe `e(j)/e(jp)/e(jdf)` and `e(first)` rows; compare binsreg output on identical bins at `$DEFAULT_SIGFIGS` (the 0.5-sigfig tolerance is the kind of ad hoc tolerance CLAUDE.md prohibits).

<a id="bld-4"></a>
#### BLD-4 [P2] — The "complete" release gate omits the I/O, transport, big-data and optimization suites, validates only the runner's own architecture's plugin, and UBSan findings cannot fail six native runners
*Component: Build, packaging, CI and validation*

- **Evidence:** Source/workflow (A16-3/7); a scratch UBSan program exits 0 without `halt_on_error`.
- **Location:** `validation/run_stata_audit.py:11`; `validation/validate_all.do:30` (unregistered: `validate_io_formats.do`, `validate_io_parity.do`, `validate_io_excel_options.do`, `validate_io_excel_formulas.do`, `validate_io_delimited_options.do`, `validate_transport_strl.do`, `validate_command_optimizations.do`, `validate_bigdata.do`, `validate_cqreg_performance.do`); `validation/test_p1_native.py:49-50,66`, `test_p2_native.py:16-17,32`, `test_sep22_native.py:88`, `test_newcommands_native.py:130-135`, `test_split_tokens_native.py:47-51`, `test_order_parallel.py:84`; `.github/workflows/build.yml:78-86`.
- **Impact:** SAS/SPSS/XLS/SHP readers and writers, Excel options, strL transport and all Windows/Linux/Intel binaries ship ungated (e.g., REG-3 could not have been caught).
- **Proposed fix:** Fold the I/O/transport/optimization suites into the release driver; run all native runners in CI with `-fno-sanitize-recover=all` / `UBSAN_OPTIONS=halt_on_error=1`; run the Stata gate once per published architecture or scope publication/docs to what is validated.

### P3 findings

<a id="core-4"></a>
#### CORE-4 [P3] — Thread-count policy ignores `OMP_NUM_THREADS`, CPU affinity/cgroup limits and Stata's `set processors`, and overwrites the OpenMP ICV on every call
*Component: Core infrastructure and shared data layer*

- **Evidence:** Code reading (A01-7); the header comment claims OMP_NUM_THREADS is respected.
- **Location:** `src/ctools_threads.c:41-57`, `:75-86`, `:111-122`; `src/ctools_config.h:391-406`.
- **Impact:** Heavy oversubscription on HPC/shared servers (e.g., 64 threads inside a 4-CPU SLURM allocation); on Linux the ICV write can leak into other OpenMP code in Stata's process.
- **Proposed fix:** Default to min(online CPUs, affinity-mask count), honor `OMP_NUM_THREADS`/`OMP_THREAD_LIMIT` (capture `omp_get_max_threads()` before the first `omp_set_num_threads`), optionally pass `c(processors)` from the ado; fix the comment.

<a id="core-5"></a>
#### CORE-5 [P3] — `ctools, update` reports "already up to date" when the host cannot be reached (r(631))
*Component: Core infrastructure and shared data layer*

- **Evidence:** Code reading; Stata maps 631 to "host not found" (A01-8, A16-15).
- **Location:** `build/ctools.ado:13-23`.
- **Proposed fix:** Treat 631 as a network failure; compare remote and local versions if an "up to date" message is wanted.

<a id="core-6"></a>
#### CORE-6 [P3] — Stale-state cleanup compares against the command from two calls earlier, so caches from an interrupted command can survive indefinitely
*Component: Core infrastructure and shared data layer*

- **Evidence:** Harness (A01-9): csort → interrupted cimport → csort leaves the cimport cache resident.
- **Location:** `src/ctools_runtime.c:226-246`; call order `src/ctools_plugin.c:300-302`.
- **Proposed fix:** Compare with the immediately previous command (or set the current command first).

<a id="core-7"></a>
#### CORE-7 [P3] — Latent shared-layer hazards
*Component: Core infrastructure and shared data layer*

- **Evidence:** Code reading/harness (A01-11/12, A02-6…11).
- **Items:**
  - Arena reset reuse is forward-only: mixed-size reset cycles retain ~15× the requested capacity (the SEP24 #19 "bounded capacity" claim does not hold; no production caller) — `src/ctools_arena.c:86-96`, `:142-158`, `:205-217`.
  - Option-string helpers match the first occurrence anywhere (`ctools_parse_bool_option("ignore=verbose verbose","verbose")` → 0; `p=` matches inside `gap=`; locale-dependent `strtod`) — `src/ctools_parse.c:20-131`.
  - Unreachable empty-load branches allocate the string array with `calloc` but free it as aligned memory (heap corruption on Windows if ever reached) — `src/ctools_data_io.c:235`, `:1184`, `:1451` vs `src/ctools_types.c:95`.
  - Empty results are inconsistent (NULL `vars` for an empty range; string columns typed numeric with NULL data when the `if` excludes all rows; storing an empty result errors; `str_maxlen` reported as 2045) — `src/ctools_data_io.c:739`, `:1031-1040`, `:1200`, `:1775-1783`.
  - `str_widths` documented as "one entry per variable position" but indexed by plugin variable index; `obs_start`/`obs_end` of 0 means "full range" so (1, 0) returns every row; `ctools_stream_var_permuted` has no callers and its documented row numbering does not match the code; public store helpers call `SF_error` from OpenMP worker threads (cwinsor `:639`, cdestring `:395`) — `src/ctools_types.h:128-133`, `:496-515`; `src/ctools_data_io.c:497`, `:922-989`, `:1027-1044`, `:1802-1893`.
- **Proposed fix:** Free successor blocks on reset (or first-fit search), or delete the unused API; tokenize options once and match whole tokens; use `ctools_aligned_free` consistently; return well-formed empty results; fix the docs/contracts; collect errors from workers and report them on the main thread.

<a id="sort-11"></a>
#### SORT-11 [P3] — Pairs-sort fallback re-runs scatter tasks already executed after a partial batch submission failure (heap overflow / duplicated rows)
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Fault-injection harness (A01, A03): failing one 24-byte work-item allocation mid-batch → ASan heap-buffer-overflow at `ctools_sort_pairs.c:79`, or 16,597 of 200,000 rows duplicated/lost.
- **Location:** `src/ctools_threads.c:298-320` (`submit_batch` drains already-queued items, returns −1); `src/ctools_sort_pairs.c:206-215` (serial fallback reruns all scatters with advanced offsets). Used by cmerge's preserve-order path (`cmerge_impl.c:990-1003`).
- **Impact:** Memory corruption instead of a clean error under memory pressure.
- **Proposed fix:** Make batch submission all-or-nothing (allocate all work items before enqueuing) or return the number queued; rebuild `all_offsets` from prefix sums before any serial retry, or treat the failure as fatal (920).

<a id="sort-12"></a>
#### SORT-12 [P3] — Minor sort-engine defects: `threads(1)` not honored by sample/counting sort at 50k ≤ N < 100k; pairs radix sorts only 32-bit keys; timsort cleanup reads uninitialized `is_numeric`; merge sort allocates an unused N×4 buffer
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Code reading / harness (A03-8…11).
- **Location:** `src/ctools_sort_sample.c:806-811,848-853`, `src/ctools_sort_counting.c:497-500` (`num_threads = nobs/MIN_OBS_PER_THREAD` → 0 → forced to 2); `src/ctools_sort_pairs.c:131-133` and `ctools_sort_pairs.h:22-37` (4 byte-passes over a `size_t` key; stability relies on input order); `src/ctools_sort_timsort.c:1336,1354-1364,1383-1388` (`malloc`ed `is_numeric` read in cleanup); `src/ctools_sort_merge.c:465,602` (unused `temp`).
- **Impact:** User thread caps violated (results unaffected); latent mis-ordering if a future caller passes ≥2^32 pair keys; UB on an allocation-failure path; 4 bytes/row wasted per merge sort.
- **Proposed fix:** Clamp thread counts to the requested maximum (use the sequential path when it is 1); assert/document pairs-sort preconditions or run enough passes for the key width; `calloc` `is_numeric`; delete the unused buffer.

<a id="sort-13"></a>
#### SORT-13 [P3] — csort compatibility and documentation gaps versus native `sort`
*Component: Sorting (csort + shared sort engines)*

- **Evidence:** Code reading (A04-14); Stata syntax check.
- **Location:** `build/csort.sthlp:56`, `docs/README_csort.md:13` (claim native sort rejects `in` — it accepts `in`); `build/csort.ado:13` (no `stable`, no `in`), `:83-88` (early exit before `return` → empty `r()`), `:330` (`r(algorithm)` reports the request, not the engine used).
- **Impact:** Ported scripts using `sort x in 1/100` or `sort x, stable` fail under the "drop-in replacement"; `r()` is empty when data were already sorted.
- **Proposed fix:** Accept `stable` (csort's final `sort ..., stable` already guarantees stability except SORT-1/6); implement or honestly document `in`; post `r()` on the early-exit path; report the resolved engine.

<a id="mrg-15"></a>
#### MRG-15 [P3] — Error codes and SPI status: plugin codes 1–5 surface as unrelated Stata errors (r(1) "Break", r(5) "not sorted"); keepusing and `_merge` writers ignore SPI return codes
*Component: cmerge*

- **Evidence:** Stata-verified (strL case returns r(5)); code (A05-17).
- **Location:** `src/cmerge/cmerge_impl.c:274,292,699,757,1092,1430,1445-1450`; `src/cmerge/cmerge_io.c:49-110`.
- **Proposed fix:** Map internal codes to 920/198/693 at `cmerge_main`; check every `SF_*` return in `cmerge_io.c` and propagate failures (completes SEP20 A29 for cmerge).

<a id="mrg-16"></a>
#### MRG-16 [P3] — Presentation and edge-case divergences: `_merge` has no value/variable label; new variables follow `keepusing()` token order; a failed `assert()` restores the master (native leaves the merged result); empty-dataset branches ignore `nolabel`/`nonotes` and `_n`+`keepusing()`; help lists `timeit` (rejected) and has a broken Remarks section
*Component: cmerge*

- **Evidence:** Stata-verified for the missing `_merge` label; code for the rest (A05-16, -18, -19, -20).
- **Location:** `build/cmerge.ado` (no `label define _merge`; `:878-910`; `:1133-1136`; `:607-767`); `build/cmerge.sthlp:3-8,45,51,55,101-113,200-204`.
- **Proposed fix:** Replicate native's `_merge` value/variable label block; create new variables in using-file order; `restore, not` before `exit 9` on assert failure (native behavior); pass `nolabel`/`nonotes` in the empty branches; fix the help and document remaining divergences.

<a id="dstr-10"></a>
#### DSTR-10 [P3] — Other cdestring/cencode/cdecode divergences and help errors
*Component: cdestring / cencode / cdecode*

- **Evidence:** Harness/source (A04-13, -15, -16, -18, -19); Stata-verified for "1d3".
- **Items:**
  - "1d3" (Stata's `real()` accepts `d`/`D` exponents → 1000) is non-numeric to cdestring; "1e+"/"2e-" are accepted as 1/2; only space/tab are trimmed (native `ustrtrim`); `ignore()` removes bytes of multibyte characters and native's `ignore("...", asbytes|aschars|illegal)` suboptions are not parsed (the suboption text becomes ignore characters, e.g. "e" turns "1e5" into "15") — `src/ctools_types.c:572-574,639-657`, `src/cdestring/cdestring_impl.c:59-136`, `build/cdestring.ado:15`.
  - `replace` with `if`/`in` silently sets excluded observations to missing (cencode/cdecode/cdestring); cencode computes `touse` but re-evaluates `if` on partially replaced data — `build/cencode.ado:99,124-128,196`, `build/cdecode.ado:55-66`, `build/cdestring.ado:132-147`.
  - cencode help says each variable gets its own label (they share one) and misdescribes extended-missing strings; cencode/cdestring help overstate strL support — `build/cencode.sthlp:56-57,93-95,137-138`, `build/cdestring.sthlp:142-143`.
  - `cdestring ..., verbose` on an empty sample errors r(111) after modifying data; `generate()` counting without a varlist differs from native — `src/cdestring/cdestring_impl.c:261-267`, `build/cdestring.ado:31-49`.
  - `cdecode` exits 108 without a message for numeric-type errors and accepts `maxlength(0)` — `build/cdecode.ado:15-19,35-37`.
  - cencode of > 65,536 distinct values returns r(198) (native r(134)); data are unchanged (Stata-verified).
- **Proposed fix:** Accept `d` exponents and require exponent digits; trim like `ustrtrim`; implement `ignore()` by character with native suboptions; reject or document `replace` + `if/in` and use the precomputed `touse`; correct help text; save timing scalars on every return; print native-style messages; use native error codes.

<a id="imp-15"></a>
#### IMP-15 [P3] — Smaller divergences from `import delimited`
*Component: cimport (delimited text)*

- **Evidence:** Stata-verified for the naming/type items; harness/bytecode for the rest (A06a-17…23).
- **Items:**
  - Variable names: non-ASCII letters are not case-folded (`ÄÖÜ` stays uppercase; native `äöü` and keeps the header as variable label, which cimport drops); fallback names ignore `case(upper)` (`v3`/`v5` vs native `V3`/`V5`); 32-byte (not character) truncation can split UTF-8; leading digits dropped (`1st` → `stplace`) — `src/cimport/cimport_impl.c:168-193`, `build/cimport.ado:426-433`.
  - `numericcols(_all)` storage types differ when a column contains non-numeric text (native int/long, cimport byte; values equal).
  - NUL bytes are kept, truncating strings at the NUL (native drops them with a note).
  - Strict mode keeps embedded CRLF (native keeps `\r`); `maxquotedrows` counted per field, not per row; `rowrange()` counts parsed non-empty rows, native counts physical lines.
  - Internal failures reported as r(601); unmatched-quote row numbers fabricated in multi-chunk mode (`cimport_impl.c:63,598-605,765,1181`).
  - `colrange()` starting beyond the last column fails with r(100) after the full load (`cimport.ado:534-556`).
  - User timers 11/12/13/99 are cleared; `threads()` does not affect chunking; help lists `locale()` for delimited (rejected) and misdescribes `asdouble` (`cimport.ado:341-342,619-622`; `cimport.sthlp`).
  - About 2.2 MB of stack arrays in the scan call chain (risk on 1 MB Windows thread stacks).
- **Proposed fix:** Implement native `makeVarName` rules (Unicode-aware, character counts, keep header as label); strip NULs with native's note; normalize CRLF in strict fields; count lines physically; propagate real error codes; validate `colrange()` before loading; stop clearing user timers; honor `threads()`; move large arrays to the heap; fix help.

<a id="exp-11"></a>
#### EXP-11 [P3] — Other export divergences and defects
*Component: cexport (delimited + Excel)*

- **Evidence:** Stata-verified for the `quote` item; code reading/harness for the rest (A07-14…20).
- **Items:**
  - With `quote`, unlabeled values of a value-labeled numeric variable are quoted (`"2147483620"`, `"-5"`); native writes them unquoted.
  - `sheet(..., modify|replace)` without a sheet name modifies a later sheet named "Sheet1" instead of the first worksheet; `r:id` lookup requires the literal `r` prefix (`cexport_xlsx_edit.inc:477-496`).
  - `sheet(..., modify)` drops attributes of empty formatted rows and drops/misorders cells or rows lacking the optional `r` attribute (`cexport_xlsx_edit.inc:245-323`).
  - `direct` requests unaligned `O_DIRECT`/`FILE_FLAG_NO_BUFFERING` I/O: always fails on Linux/Windows (r(693)), silently ignored on macOS (`cexport_io.c:164-171,591-595`).
  - `keepcellfmt` help describes behavior the code does not implement (`cexport.sthlp:240-245`; `cexport_xlsx.c:73-75,1450-1456`).
  - Several failures return 198/603/610 with no message; XLSX size limits (16,384 columns / 1,048,576 rows) are checked only after the full data load (`cexport_xlsx.c:1629-1693`, `cexport_xlsx_edit.inc:452-458`, `cexport_impl.c:856-857,916`).
  - CSV exponent derived from `floor(log10(x))` rounds up just below powers of ten (999.99999999999989 → `1000`) (`cexport_format.c:177-181`).
  - Unchecked `mz_zip_writer_end`/`msync`; `missing()` numeric detection accepts hex; Windows fallback over 64 batch entries; workbooks without `styles.xml` get an unrelated styles part.
- **Proposed fix:** Match native quoting of unlabeled values; select the first worksheet when `sheet()` is omitted and match relationship ids namespace-agnostically; copy empty-row start tags and track implicit cell/row positions; remove or correctly implement `direct`; fix `keepcellfmt` docs; emit specific messages and validate limits before loading; derive the exponent from the rounded digits; check all results.

<a id="reg-11"></a>
#### REG-11 [P3] — Wrapper compatibility and documentation divergences from reghdfe
*Component: creghdfe + shared estimation kernels*

- **Evidence:** Code reading (A09-9); help-file check (A16-9).
- **Items:** `[pw=w], vce(unadjusted)` resets the pweight-forced robust VCE (silently different SEs); `vce(rob)` abbreviation rejected; string FE variables mark out every observation ("no observations", r(2000)) whereas reghdfe accepts them; `e(N_hdfe_extended)` holds the number of mobility groups; `e(vcetype)` never set (no "Robust" header); bare `residuals` rejected; help advertises `timeit` and `iter:ate` (rejected), says `absorb()` is required, and still says savefe/groupvar "may not contain values".
- **Location:** `build/creghdfe.ado:31-33`, `:91`, `:240-245`, `:715`; `build/creghdfe.sthlp:23`, `:30`, `:53`, `:83-84`, `:96-97`.
- **Proposed fix:** Force robust for pweights (or reject the combination); accept abbreviations and bare `residuals`; encode string FEs with `egen group()`; set `e(N_hdfe_extended)`, `e(vcetype)`; fix the help.

<a id="reg-12"></a>
#### REG-12 [P3] — Error-path and robustness defects in the creghdfe C core
*Component: creghdfe + shared estimation kernels*

- **Evidence:** Code reading; UBSan harness for the extreme-FE case (A09-10).
- **Items:** most early returns leak `obs_map`, `cluster_raw_values`, `weighted_counts_orig[]`, `group_assignments`; `remap_and_count` sets `*num_levels` before its `calloc` and the caller ignores the return (NULL dereference on OOM, `creghdfe_utils.c:135-139,195-196` → `creghdfe_regress.c:402-411,529`); `output_rc` from a failed sample-flag store is overwritten by the VCE status (`:1615`, `:1618`); `compute_vce_unadjusted` is `void` with no finiteness check (saturated models post NaN V → r(504); `creghdfe_vce.c:18-32`); FE values ≥ 2^63 → UB and heap out-of-bounds writes (`creghdfe_utils.c:128,233,243`); `(ST_int)sum_weights` overflows for Σfw > 2^31−1 (`creghdfe_regress.c:1048,1613`); a −1 singleton count on allocation failure is used as a count.
- **Proposed fix:** Single cleanup label; check `remap_and_count`; separate `store_rc`; status + finiteness check for the unadjusted VCE; reject |v| ≥ 2^53 on the counting path and compute ranges unsigned; carry fweight N as double; treat negative singleton returns as 920.

<a id="iv-25"></a>
#### IV-25 [P3] — Other ivreghdfe-compatibility defects
*Component: civreghdfe*

- **Evidence:** Code reading (A10a-17/18/21/22, A10b-21/23).
- **Items:** `e(partial_r2_#)` stores R²_full − R²_reduced instead of the partial R² (0.0219 vs 0.0793) and `F_first#`/`partial_r2_#` labels shift after collinear drops with stale scalars (`civreghdfe_estimate.c:1718-1737`, ado `:1456-1480`); `b0()`, `coviv`, `small`, `first` accepted but ignored; `robust cluster()` rejected; ivreg2 kernel abbreviations (`bar`, `par`, `tru`, `qua`) rejected; fweights accepted with HAC/DK; Kiefer posts `N_clust`/"Cluster"; wrong Stock-Yogo critical values for L > 3; `e(j)`/`e(jdf)`/`e(jp)` never posted (robust J posted as `e(sargan)`); fvexpand/fvrevar run before the final markout (base level chosen on a larger sample); `tsset` re-sorts the user's data without `sortpreserve`; kiefer's `by ivar:` needs sorted data; global temp scalars; string `absorb()` read as doubles; weights not validated.
- **Proposed fix:** `partial_r2 = (r2_full − r2_reduced)/(1 − r2_reduced)` indexed by original position; implement or reject ignored options; accept ivreg2 abbreviations; mirror ivreg2's fweight/HAC rule; correct Stock-Yogo tables; post `e(j)` family; markout before `fvexpand`; add `sortpreserve`; use tempnames; validate absorb types and weights.

<a id="qreg-8"></a>
#### QREG-8 [P3] — Other cqreg defects: VCE numerics (uncentered X'X, +1e-10 regularization → wrong/negative/asymmetric V for offset regressors; VCE failure posts zeros with rc 0); non-vertex solutions on tied data; `nopreprocess(-1)` data race and unverified convergence; O(N·G) cluster id mapping and NaN V for one cluster; stored-result gaps
*Component: cqreg*

- **Evidence:** Harness (A12-10…A12-16).
- **Location:** `src/cqreg/cqreg_vce.c:61-146`, `:197-259`, `:693-745`, `:878`; `src/cqreg/cqreg_regress.c:1073-1077`; `src/cqreg/cqreg_fn.c:1315-1391`; `src/cqreg/cqreg_ipm.c:951`, `:992-1054`; `src/cqreg/cqreg_linalg.c:90-117` (shared loop counter `k` in a parallel loop); `build/cqreg.ado:194`, `:264-271`, `:446`; `build/cqreg.sthlp:100-264`; `docs/README_cqreg.md:63-114`.
- **Items:** `e(rank)` counts base/omitted columns; `e(vcetype)`, `e(kernel)`, `e(title)` missing (help documents `e(title)`); `quantile(25)` percent form rejected; scalars leak and a stale `__cqreg_G` breaks later calls; help/README state the wrong iid method and a 1/n V formula; `is_collinear` leak on one error path; doubled peak memory.
- **Proposed fix:** Center/scale before inverting, symmetrize, pivoted generalized inverse, hard error on VCE failure; select K independent rows by pivoted elimination and verify KKT; make `k` loop-local and verify optimality (or remove the experimental path); map cluster ids directly and error for G < 2; post `e(rank)`, `e(vcetype)`, `e(title)`; accept percent quantiles; clean scalars; fix docs.

<a id="ppml-5"></a>
#### PPML-5 [P3] — Smaller ppmlhdfe divergences
*Component: cpplmhdfe*

- **Evidence:** Harness/code (A11-7…12).
- **Items:** slope-only absorb terms (`g#c.z`) skip the all-zero-group screen (`cpplmhdfe_separation.c:53`); singleton re-screening after a separation refit ignores fweights (`cpplmhdfe_separation.c:52`, ado `:629-636`); an all-base factor varlist (`i.one`) errors r(198) instead of fitting the FE-only model (ado `:367-395`, `cpplmhdfe_irls.c:548`); ReLU hitting its iteration limit silently reports no separation (4 of 60 LP-verified separated designs missed; `cpplmhdfe_irls.c:420-423`); the bare-`verbose` rewrite breaks a variable named `verbose`, `verbose(-1)` enables output, `itolerance()` unclamped, aweights accepted with scale-dependent `e(ll)` (ado `:12`, `:26`; `cpplmhdfe_irls.c:1368`); `obs_map` leak on the remap-failure path, missing messages and an unchecked inversion status (`cpplmhdfe_irls.c:223`, `:705-718`, `:1248`, `:1444`, `:1447`, `:1806`).
- **Proposed fix:** Match the reference screens (or document); Σfw==1 singleton rule on refits; fall through to the zero-regressor branch; warn when ReLU exhausts its limit; restrict the `verbose` rewrite to the options, clamp `itol`, reject or normalize aweights; fix error paths.

<a id="bin-12"></a>
#### BIN-12 [P3] — Options that do nothing or fail late, result-metadata hygiene, scalability and documentation
*Component: cbinscatter*

- **Evidence:** Code reading (A13-16…19).
- **Items:** `linetype(connect)` draws unconnected points; `reportreg` does nothing; `legend()` ignored for a single series; binscatter options `rd()`, `xq()`, `medians`, `noaddmean`, `savegraph()`, `replace` fall into `*` and fail inside `twoway` after all computation (or are ignored under `nograph`); `e(by)`/`e(controls)` hold temporary variable names; `savedata()` always overwrites and stores dense by-codes only; scalars/matrices left behind on error and stale `e()` survives; binsreg+absorb memory/time scale with N × nquantiles (8 GB at 50M obs with 20 bins); help/README claim OpenMP (none), "same HDFE algorithm as creghdfe" (it uses its own projections), equal-sized quantile bins (BIN-1), 7 `e(bindata)` columns (there are 5).
- **Location:** `build/cbinscatter.ado:13-35`, `:76-78`, `:173-188`, `:214`, `:284-299`, `:350-366`, `:412-457`, `:576-642`; `src/cbinscatter/cbinscatter_impl.c:861-871`; `src/cbinscatter/cbinscatter_binsreg.c:277-402`; `build/cbinscatter.sthlp:71-88`, `:109-110`, `:276`, `:283`; `docs/README_cbinscatter.md:7`, `:155`, `:188`, `:196`.
- **Proposed fix:** Add `connect(l)`; implement or reject binscatter options at parse time; honor `legend()`; post user-facing names and original by values; add `replace` to `savedata()`; centralize cleanup; use sparse accumulation for binsreg+absorb; correct docs.

<a id="xio-12"></a>
#### XIO-12 [P3] — Smaller statistical-format and `.xls` export defects
*Component: cimport/cexport Excel and statistical formats*

- **Evidence:** Harness/code (A08-10…19).
- **Items:** `.xls` numeric header cells give labels like "2023.000000"; editing an existing `.xls` with an empty SST produces a file xlrd/pandas cannot open; minimal BIFF8 record set (Excel compatibility unverified); `cexport sasxport5` saturates |x| > ~4.5e74 and drops trailing blanks where native errors; `cimport sasxport8` skips native's `compress`; Windows charset adapter lacks Mac/CJK code pages and replaces invalid bytes silently; dBase date export uses `(int)x` (negative fractional dates shift a day; UB for |x| ≥ 2^31); XPORT5 `._` imported as `.` (native `.u`) and reserved/duplicate names silently renamed; exported files get mode 0600 (EXP-6); import cache not freed on Break; Break reported as r(692); empty `.xls` sheet described as `IW1:A1`.
- **Location:** `src/io/cio_xls.inc:153-160`; `src/io/cexport_xls_edit.inc:119-147`; `src/io/cexport_xls.c:340-417`; `src/io/cio.c:94-104`, `:1122`, `:1439`; `build/_cio_import.ado:221-224`; `src/io/cio_iconv_win.h:13-59`, `:95-126`.
- **Proposed fix:** Format numeric headers with `%.16g`; re-pack SST/EXTSST on edit; emit the standard minimal BIFF8 record set; pre-scan and raise native errors for XPORT limits; add `compress`; extend the code-page table and use `MB_ERR_INVALID_CHARS`; `floor()` with range checks for dBase dates; map `._` to `.u` and display renames.

<a id="psm-3"></a>
#### PSM-3 [P3] — Smaller cpsmatch/psmatch2 semantic differences
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Evidence:** Harness/source (A15-11/15).
- **Items:** caliper inclusive (`<=`, psmatch2 `<`); `radius` without `caliper()` silently uses 0.25; kernel honors `caliper()`; `noreplacement ties` ignores ties; `_weight` 0 (psmatch2 missing) for unused rows; help's "(and vice versa)" for `common` is wrong; ties use a 1e-12 tolerance; noreplacement scans used controls linearly (O(n²) with dense/tied scores).
- **Location:** `src/cpsmatch/cpsmatch_impl.c:492`, `:501`, `:630-632`, `:657`, `:863`, `:970-1082`, `:1042-1043`, `:1365-1374`; `build/cpsmatch.ado:45-50`; `build/cpsmatch.sthlp`.
- **Proposed fix:** Strict `<` caliper; require or mirror psmatch2 for radius; ignore caliper for kernel (or document); implement or reject `noreplacement ties`; store missing for unused `_weight`; fix help; use skip pointers for available controls.

<a id="smp-3"></a>
#### SMP-3 [P3] — Sampling-command divergences from `sample`/`bsample`
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Evidence:** Harness/source (A15-10/12/14).
- **Items:** percent sample size rounded differently at .5 boundaries (`sample 29` on 50 obs keeps 15, csample 14; `csample_impl.c:557-560`); native `sample #, count` syntax fails, `count(0)` rejected, not byable; the keep-flag store return code is ignored (a store failure would silently drop sampled rows; `csample_impl.c:607-611`); cbsample has no singleton-cluster check (bsample r(460)), no `idcluster()`, no expression for the count or `w()` abbreviation, r(198) instead of r(498) for n < 1, and O(strata × clusters) cluster lookup (`cbsample_impl.c:598-649`, `build/cbsample.ado:26-41,61-69`).
- **Proposed fix:** `floor(count*percent/100 + 0.5)`; accept native syntax and `byable(onecall)`; check store rc; add the singleton check, `idcluster()`, expression evaluation, native error codes, and a running cluster cursor.

<a id="win-4"></a>
#### WIN-4 [P3] — cwinsor option/type divergences from winsor2
*Component: csample / cbsample / cwinsor / cpsmatch*

- **Evidence:** Stata-verified for types and `cuts(0 100)`; source for the rest (A15-13).
- **Items:** default is in-place replace (winsor2 creates `_w`/`_tr` variables — porting `winsor2 x` destroys `x`); `cuts(0 100)` accepted (winsor2 refuses); low > high swapped by winsor2 but an error in cwinsor; new variables always `double` without labels (winsor2: source type + label, `label` option unsupported); groups with < 3 non-missing values unchanged; in suffix mode `.a` survives (winsor2 generates `.`).
- **Location:** `build/cwinsor.ado:31`, `:36-53`, `:103-121`; `src/cwinsor/cwinsor_impl.c:37`, `:554-557`, `:596-599`.
- **Proposed fix:** Consider winsor2-compatible defaults (or a prominent note), mirror winsor2's `cuts()` handling, create outputs with the source type and copy labels.

<a id="rng-10"></a>
#### RNG-10 [P3] — Other crangestat/csplit/crangejoin divergences, help errors and resource issues
*Component: crangestat / crangejoin / cipolate / csplit*

- **Evidence:** Code reading/harness (A14-9…14; A02-1 for the quickselect); Stata-verified for the csplit item.
- **Items:**
  - Abbreviated key/source names pass `confirm variable` but have no `_varpos_` entry, so later positional plugin arguments shift and the command fails with r(198) — `build/crangestat.ado:57-66`, `:202-212`, `:263-291`.
  - rangestat syntax not supported: `[if] [in]`, variable-valued interval bounds (r(198)), `(stat) varlist` with several variables, stats `obs`, `missing`, `cov`, `corr`, `reg`, options `casewise`, `describe`, `local()`, abbreviations `i()`/`excl` — `build/crangestat.ado:32-96`, `:120-238`.
  - Help calls kurtosis "excess kurtosis" (it is m4/m2², as in rangestat and `summarize`), and describes first/last as RNG-2's intended behavior — `build/crangestat.sthlp:128-129`, `:141`.
  - median/IQR/percentiles use a private Lomuto quickselect that is O(w²) per window with tied values (see WIN-2) — `src/crangestat/crangestat_impl.c:128-168`.
  - Work arrays of `threads × max_group_size` doubles are always allocated (5 GB for 64 threads × a 10M-row group, even for mean/count); sparse min/max tables need ~320 bytes/row; tables built before a later allocation failure leak — `src/crangestat/crangestat_impl.c:543-603`, `:1713-1733`, `:1836-1845` (plausible).
  - `csplit, parse("")` exits r(198) without a message; native `split` falls back to space parsing (Stata-verified: `split s, parse("")` creates r1–r3) — the ado compares with a one-character `"` instead of `""` — `build/csplit.ado:21`.
  - `crangejoin` rejects strL variables (rangejoin accepts them) and several rejections exit with bare return codes — `build/crangejoin.ado:17`, `:34`, `:50-55`, `:91`, `:94-97`.
- **Proposed fix:** `unab` key/source names and compute positions with `: list posof`; implement or document each syntax gap; fix the help; three-way-partition quickselect; allocate work arrays only for stats that need them and cap sparse tables; compare with `""""` in csplit; document strL limits and print messages before `exit`.

<a id="bld-5"></a>
#### BLD-5 [P3] — Other validation-framework false-pass paths
*Component: Build, packaging, CI and validation*

- **Evidence:** Source (A16-10…13, A16-17).
- **Items:** `sigfigs` treats any two values below 1e-14 as identical (V comparisons for large-scale regressors are vacuous; `validate_setup.do:199-205`); estimator helpers never `ereturn clear`, check `e(cmd)`, or compare coefficient names (label/order defects such as IV-1 pass; `benchmark_qreg` skips V when dimensions differ, `:620`); always-pass or rc-only checks counted as passes (`validate_creghdfe.do:2924-2948`, `:3640-3647`; `validate_cpplmhdfe.do:596-630`, `:801-812`, `:958-970`); `benchmark_psmatch2` ignores partially missing `_pscore`; `benchmark_ppmlhdfe` silently skips missing `ll`/`r2_p`/`N_clust`; `benchmark_export` falls back to 7-sigfig numeric comparison; `validate_cimport.do` renames variables positionally (names never compared); cwinsor is claimed to replace winsor2 but validated only against `gstats winsor`; `validate_setup.do:22-23` lets a sibling `../build` shadow the code under test.
- **Proposed fix:** Remove the absolute floor (compare standardized quantities); `ereturn clear` + `e(cmd)` + `colfullnames` checks; full comparisons or equal-rc requirements in the rc-only tests; compare missing patterns; drop tolerance fallbacks; compare names; add a winsor2 reference; pin the adopath and assert the loaded plugin revision.

<a id="bld-6"></a>
#### BLD-6 [P3] — Help files document spellings the wrappers reject, and other packaging/build hygiene
*Component: Build, packaging, CI and validation*

- **Evidence:** Stata-verified: `creghdfe price mpg, absorb(foreign) iter(50)` and `... timeit` → r(198). Source (A16-9, A16-14, A16-16, A16-18); `timeit`/`iter()`/`locale()` behavior by Stata `syntax` rules.
- **Items:** `timeit` in cmerge/creghdfe help (rejected); creghdfe `iter:ate` (the syntax requires `ITERATE`); cimport delimited `loc:ale`, `decimals:`, `groups:`; cbinscatter `yt:`/`xt:`; civreghdfe `ff:`, `noret:`, `noomit:ted`/`omit:ted` (help markup pasted into `syntax`, `build/civreghdfe.ado:65-66`); cqreg documents `e(title)` but never sets it; `stata.toc` still says Stata 14.0; nothing asserts distribution plugins contain OpenMP (a missing `libomp.a` silently builds serial binaries, `Makefile:107-147`); `make clean` deletes tracked distribution plugins that usually cannot be rebuilt locally and is a no-op under MSYS2 (`Makefile:621-639`); ~340 files under `temp/` are neither tracked nor ignored and required release inputs are untracked; stray `build/test1.csv`/`build/test2.csv`.
- **Proposed fix:** Align help and `syntax` (add an automated abbreviation checker); set `e(title)`; sync `stata.toc` via `sync_release.py`; add `REQUIRE_OPENMP=yes` for CI and assert OpenMP symbols; build into a separate output directory and stage into `build/` only via packaging; ignore or clean `temp/`; track release inputs.

## Areas checked and found correct

These results narrow where fixes are needed. Each was checked by harness or Stata comparison.

- **Command matches in Stata:**
  - `cipolate` matches `ipolate`: 29 configurations including `by()`, `epolate`, ties, float storage and missing x/y.
  - `csplit` matches `split`: 20 configurations including `parse()`, `limit()`, `notrim`, `destring`, `if`/`in`.
  - `crangejoin` matches `rangejoin`: 12 configurations including numeric/variable/missing bounds, `by()`, `prefix()`/`suffix()` and missing keys.
  - The A14 review found these three commands otherwise clean. Native harnesses compared them against emulations of the reference ados: `cipolate` in 600 random cases, `csplit` in 3,000 and `crangejoin` in 400, including the parallel paths, under ASan/UBSan. The only exceptions are RNG-10's `csplit, parse("")` and crangejoin's strL/message gaps.
  - `cwinsor` matches an exact `_pctile` reference apart from WIN-1/WIN-3.
  - `cmerge` join content matches `merge` for 1:1, m:1, 1:m and m:m on numeric, string and composite keys with missing keys. The divergences are the MRG items.
  - `cimport delimited` matches `import delimited` on 19 of 25 torture files: quotes, embedded newlines, ragged rows, BOM, CRLF without blank lines, long strings, leading zeros, big integers, `stringcols()`, `colrange()`, `rowrange()` with headers, `groupseparator()`.
- **Sort engines:** every engine other than timsort (SORT-1/6) orders numeric keys (±0, all 27 missing codes, extreme magnitudes) and string keys (bytes ≥ 0x80, prefixes, 2,045-byte strings) correctly and stably. This held across 1–400,001 rows, thread limits 1/2/3/4/8, and builds with and without OpenMP.
- **Estimation kernels:**
  - The creghdfe 2-FE solver and robust/cluster VCE match a dense dummy-variable regression to ≥ 14 significant figures for all weight types.
  - The civreghdfe 2SLS/LIML/Fuller algebra and the single-endogenous KP statistics match ivreg2/ranktest formulas.
  - cpplmhdfe coefficients, VCE, ll, deviance and separation (default settings) match an exact MLE and LP ground truth in randomized trials.
  - The cqreg main solve returns the exact optimal vertex on continuous data; its bandwidth rules and iid/robust fitted VCE formulas match `qreg.ado`.
- **Core runtime:**
  - The Eisel–Lemire parser and its powers-of-ten table are bit-identical to `strtod` over 160M random inputs.
  - The dispatcher's `threads()` parsing, the thread pool shutdown (SEP24 #18), arena chaining (SEP24 #19), the string hash and the Mata label writer are correct.
  - No Stata display, scalar or macro SPI calls are made from worker threads.
- **Earlier fixes confirmed complete:**
  - SEP20: A02, A03, A04 (absorb rejected in cqreg), A05, A06, A09, A11, A13, A15, A16, A18, A19, A20, A21, A22, A23, A24, A26, A27, and A29 except in cmerge's own writer (MRG-15).
  - SEP22: S01–S14.
  - SEP24: #1, #3, #7, #8, #10, #12, #15, #16, #17, #18, #19. #2 holds apart from the REG-1/REG-2 bypasses. #5 holds apart from the serial-60 inconsistency (XIO-11).
- **Earlier fixes found incomplete by this audit:**
  - SEP20 A07: by-group counting in the plugin (BIN-7).
  - SEP20 A08: XLSX shared strings (XIO-2).
  - SEP20 A17: the cqreg auxiliary solves (QREG-1).
  - SEP24 #4: cbinscatter HDFE and the civreghdfe FWL loop (BIN-9, IV-6).
  - SEP24 #6: cencode's label path and cexport's delimiter (DSTR-8, EXP-3).
  - SEP24 #9: timsort (SORT-6).
  - SEP24 #11: IV two-way/Kiefer VCE and the creghdfe unadjusted VCE (IV-22, REG-12).
  - SEP24 #13: the `__ctools_strw` hint and cimport macros (CORE-2, IMP-11).
  - SEP24 #14: the creghdfe F statistic (REG-4).
  - SEP24 #20: no plugin/ado version handshake (CORE-3).
  - SEP24 #21: documentation still contradicts behaviour in the places listed under BLD-6, SORT-13 and DSTR-10.

## Reproduction material

Every entry includes a reproduction or trigger. The harness sources, fixtures and Stata driver do-files used for this audit are in this session's scratch directory, `/private/tmp/claude-501/-Users-Mike-Documents-GitHub-stata-ctools/e09abaa8-3dcf-4652-8612-b1f6905197af/scratchpad/`:
- `A*_*/`: per-area native harnesses and fixtures.
- `findings/`: the unabridged per-area reviewer write-ups, which carry the IDs cited above as "A09-4", "A10b-3" and so on.
- `stata/`: the Stata do-files and logs behind every "Stata-verified" statement.

That directory is temporary. Copy anything worth keeping, such as the fixtures, into `validation/` when the corresponding fix is made.
