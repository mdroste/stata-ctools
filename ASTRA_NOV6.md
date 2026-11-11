# ctools release audit and remediation plan — November 6, 2026

## November 11 remediation status

Rechecked the report against the current working tree and implemented eleven batches
covering **45 findings**, including all seven P1 implementation defects. Three
findings already had repairs in the tree. **One finding remains open** below.
Existing unrelated changes were preserved. No release binaries were replaced and
no commits or publication were made.

These are implementation statuses, not certification that every acceptance
criterion in the original report has run on every supported platform.

| Batch | Findings | Changes |
|---|---|---|
| Data protection and samples | 01, 03, 06, 10, 16, 17, 34 | Check streaming sort and merge transfers; recover failed streaming sorts with a wrapper snapshot; reject malformed Parquet prefixes; freeze incoming-result samples; use owned bootstrap temporary variables; restore overwritten winsor targets. |
| Result and option contracts | 02, 04, 05, 07, 15, 20, 21, 28, 29, 30, 39, 44, 45, 48 | Complete shortcut/fallback results, retain using labels without frames, guard Parquet row counts before decoding, respect quiet mode, check quantile sample size after rank reduction, validate solver/thread options, handle flexible-only statistics and extended missing bounds, abort failed label backups, correct bootstrap size-expression order, initialize empty sample results, correct Linux diagnostics. |
| Numerical and text handling | 11, 22, 24, 35, 36, 37, 42 | Quote SAS literals; retain double-precision external bin IDs; recompute ill-conditioned local radius means; normalize extreme weights; retain complete ignore lists; publish empty destring timings; allocate delimiter metadata dynamically. |
| Publication and validation | 09, 12, 14, 19, 25, 40, 46, 47 | Propagate required scalar-save failures in matching, encoding, sampling and grouping; repair the remaining Windows runner preparation mock; compare CSV cells; reject unexpected singleton and PPML error codes. |
| HDFE cancellation and PPML follow-up | 26, plus the `hard5` regression | Poll Break only on the calling thread; latch cancellation atomically for workers; keep polling while parallel teams finish; unwind shared/individual/PPML solvers and bootstrap batches with r(1). Use compensated weighted means, working residuals and intercept covariance transformations for PPML. |
| Grouped/parallel command coverage | 31, 33, 43 | Add independent grouped interpolation, split and join oracles; actual OpenMP worker checks; allocation and early/middle/late publication failures; synthetic join output-limit boundaries; licensed staged-output/rollback tests. |
| Destring conversion coverage | 38, plus two parity defects | Exercise the production converter/parser with independent values, UTF-8, long ignore lists, filtering, OpenMP, allocation and publication faults. Preserve valid UTF-8 sequences under illegal-byte removal and hexadecimal missing encodings under force. |
| Binscatter transformation coverage | 23, plus two failure paths | Cover large labels, weighted median ties, sampled quantiles, polynomial predictions, collinear controls and absorption. Propagate median-sort allocation failures and reject empty sampled-quantile inputs. |
| HDFE slopes and DoF metadata | 27, plus constant-slope centering | Assert analytic weighted residuals, orthogonality, idempotence and every native metadata field. Compare per-effect Stata DoF/owner fields. Center slopes relative to an observed group value so constants remain exactly redundant at extreme scales. |
| Label mapping and failure coverage | 41, plus protocol/result checks | Exercise full encoder/decoder mappings, parallel chunks, malformed protocols, transfer/allocation/publication failures. Validate terminators, duplicate codes and decoder integer bounds; propagate decoder result-save errors. |
| Complete Parquet file coverage | 08, plus parser bounds | Exercise complete files from two producers, exact cells, codecs/page versions, malformed metadata/data, allocation faults and byte mutations. Guard row-group count overflow, require leading magic, and reject inconsistent uncompressed sizes before allocation. |
| Already repaired before this pass | 13, 18, 32 | The creghdfe estat handler uses `e(df_m)+1`; obsolete IV saving blocks are absent; interpolation already selects the singular message. These were not reimplemented. |

Additional regressions exposed by validation were corrected: winsor's numeric
`p()` alias now coexists with `prefix()`. CSV fixtures without embedded line breaks
compare exact records while accepting native Windows CRLF terminators. The
nonunique `qreg price mpg weight i.rep78` median fixture now uses the existing
independent coefficient-loss comparison, retaining scalar and covariance checks;
both coefficient vectors attain the same minimum loss. Numerical tolerances were
not relaxed.

### Remaining work

| Finding | Current status and next work |
|---|---|
| 49 | Open: this local serial Windows candidate is not artifact-bound, signed, multiplatform release qualification. |

### Batch eleven verification

- `validation/native/parquet_paths.c` passes UBSan with recovery disabled,
  allocation accounting and tail canaries. It exercises 51 valid full files,
  including independent PyArrow-produced pages and independently compressed
  Zstd data, with exact expected cells and row ranges. Repository-produced
  fixtures cover six codecs, both page versions, nulls and four encoding
  combinations. Every allocation in the 48 small valid files is faulted in turn.
- Nineteen explicit malformed files cover footer/page sizes, dictionary indexes,
  definition levels, delta prefixes, row-group counts and six damaged codecs.
  Another **1,644** byte mutations/truncations run through the complete parser
  and column reader. Test allocation caps bound adversarial runs; they are not
  a claimed universal production memory limit.
- This exposed signed overflow summing row groups, absent leading-magic
  validation, and allocation before validating uncompressed page sizes. All
  three are repaired. Restoring the overflow or missing-magic implementation
  fails the negative controls. The earlier delta-prefix negative control remains.
- The full licensed `s2_parquet.do` passes **172/172**. The new
  `validation/support/parquet_malformed.do` passes **38/38** exact-error and
  data-preservation checks in default/snapshot modes, including the oversized
  row-count guard. Final serial HDFE and label native fixtures also pass.
- Linux's native runner now requests ASan/UBSan for this corpus, and a runner
  check verifies sanitizer selection reaches compilation. Local Zig UBSan is
  verified by a deliberate overflow mutation. Its ASan flag did not detect an
  intentional heap overflow, so the local nominal ASan run is explicitly
  **not counted as ASan evidence**. Real ASan execution remains a release check.
- Logs and build drivers are under `temp/astra-nov11-batch11/`, including
  `native-final.log`, `negative-overflow.log`, `negative-magic.log`,
  `stata-parquet.log`, `rollback.log` and the rejected ASan probe evidence.

### Batches eight through ten verification

- `validation/native/bins_paths.c` passes serial and OpenMP UBSan, with explicit
  float-cast-overflow checking locally. It verifies exact external labels through
  2^52, weighted ties, sampled cutpoints, few-value fits, collinear controls,
  absorbed effects and allocation cleanup. Large fits cross the block-parallel
  threshold. Reintroducing the label cast or ignoring median-sort allocation
  failure makes the tests fail. The full `s2_bin.do` reference suite passes
  **65/65**, including binscatter option/graph comparisons and winsor follow-ups.
- `validation/native/hdfe_paths.c` verifies slope-only/intercept-plus-slope,
  multiple/constant/zero slopes, weighted residuals, orthogonality, idempotence,
  nested/disconnected categories, singleton removal, pool/thread counts and
  1e-200/1e200 scaling. Every native result field is captured and asserted.
  Denominator and metadata mutations fail. The complete Stata `chdfe.do` passes
  **131/131**, now comparing `Kj`, `Mj`, exact/nested flags, owner mappings and
  `N_hdfe_extended`, with an additional analytic extreme constant-slope check.
- `validation/native/labels_paths.c` passes serial and actual OpenMP UBSan.
  It covers bulk/stream encoding, existing-label extension/noextend, smallest-code
  selection for duplicate text, direct-table/binary-search decoding, negative
  and sparse codes, extended missings, Unicode, empty selections, malformed
  mappings, allocation/transfer/writer/scalar/macro failures and cleanup. The
  200,001-row cases assert worker stores and inject early/middle/last failures.
  Negative controls reject unchecked encoder/decoder saves and missing terminators.
- Full Stata encoding and decoding pass **248/248** and **229/229**. The separate
  label fault plugin passes **24/24** checks: encoder failures restore existing
  labels and publish no output; decoder failures discard partial output and
  return complete correct results through its documented Stata fallback.
- Logs and temporary build drivers are under `temp/astra-nov11-batch8/`,
  `temp/astra-nov11-batch9/` and `temp/astra-nov11-batch10/`. These remain local
  candidate builds, not release qualification.

### Batch seven verification

- `validation/native/destring_paths.c` passes serial and actual OpenMP with
  UBSan and assertions enabled. Independent cases cover decimal/d/hex exponents,
  float bounds, percentages, extended missings, escaped/long ignore wires,
  UTF-8 character/byte/illegal modes, filtered/empty samples, one/multiple
  variables, force/nonforce, allocation failures, early/middle/late stores,
  and every result-scalar/macro failure. Cleanup and fresh success are checked.
  Mutations restoring truncated ignore lists, absent empty timings, or broken
  UTF-8 sequence traversal each fail the corresponding test.
- This exposed an additional defect: illegal-byte removal revisited valid UTF-8
  continuation bytes, damaging Unicode whitespace. The scan now preserves whole
  valid sequences. Native hexadecimal missing encodings are also retained under
  force before Stata's storage/compression, while counted as nonnumeric.
- The complete current Stata destring suite passes **229/229**. A separate
  test-only plugin passes **28/28** publication rollback checks in
  `validation/support/destring_publication.do`, covering generate/replace,
  empty samples, every scalar save and the counts macro. Logs and temporary
  build drivers are under `temp/astra-nov11-batch7/`.

### November 11 verification evidence

- New `validation/stata/astra_nov6.do`: **89 focused checks pass** across wrapper,
  contract, numerical, bootstrap and follow-up cases. They include all three
  estimators' pairs/wild bootstrap scratch ownership, post-permutation rollback,
  incoming `r()` samples, old merge-label control flow, and PPML exact errors with
  data and prior-estimate preservation.
- New `validation/native_astra.py`: **eight native harness configurations pass**
  with assertions explicitly enabled (`-UNDEBUG`): Parquet, merge, stream sort,
  encoder/group publication, both samplers, matching, and existing new-command
  allocation/cache tests. Matching also checks the results of its intentional
  workspace-allocation fallback. Pre-fix Parquet/merge/sort implementations each
  fail the corresponding new test. An early assertion-disabled run was discarded.
- Broader licensed Stata testing passed 2,217 checks across sort identity, merge,
  encode, destring, sample, bsample, rangestat, winsor, split and their selected
  contracts. Additional passing runs include creghdfe (385), cpsmatch (314),
  quantile sample/scale (371), and publication/replay cases. After correcting the
  two test comparisons above, standalone cexport **243/243** and cqreg **300/300**
  pass. The later executable option run passes **690/690** cases.
- The broader run is **not wholly green**: cbinscatter's two `binsreg` comparisons
  cannot run because that reference package is absent; 223 other checks pass.
  Its original PPML `native upstream hard5: r(9)` failure was subsequently
  resolved by batch five, with the full current PPML suite passing (see below).
- Runner self-test, all ten runner unit tests, seed checker, documentation,
  release metadata, architecture and `git diff --check` pass. The option inventory
  was regenerated and checks 526 options.
- Evidence is under `temp/astra-nov11/`: `final-validation/`, `final-estimators/`,
  `cexport-final.log`, `cqreg-final.log`, `followup.log`, `native-final.log`, and
  `negative-controls.log`, and `final-smoke.log` (all five focused cases and the
  complete split suite rerun after the last rebuild). Earlier aggregate logs retain failures subsequently
  resolved by the named reruns; they are not represented as successful runs.
- Native testing used a freshly compiled temporary Windows plugin and current ado
  files, StataNow 19, and the locally available Zig C compiler. This is a serial
  build with a local CPU-dispatch compatibility stub, not a distribution binary.
  Repository-wide sanitizer/OpenMP certification, actual Stata 14/15,
  interactive Break and other operating systems remain unverified; batches
  five through eleven add targeted Windows OpenMP/UBSan coverage. The no-frame
  merge branch was forced on Stata 19 for its regression. Batch eleven adds
  a complete oversized-file rejection and rollback fixture.

### Batch five verification

- `validation/native/hdfe_cancel.c` passes serial and actual OpenMP execution,
  including UBSan with recovery disabled. All shared transforms/accelerations,
  weighted/unweighted inputs, individual LSMR, fixed-effect recovery and PPML
  projections return cancellation, and a fresh solve succeeds afterward. Tests
  assert callback thread identity, exercise a host waiting for unfinished workers
  and nested teams, and verify bootstrap cancellation frees all workspaces without
  folding incomplete replicates. The pre-fix shared solver fails the new test.
- `validation/support/hdfe_cancel.do` passes **173 licensed Stata checks** using a
  separate plugin with deterministic host-poll injection. `chdfe` generate/clear
  modes preserve original data and labels across every transform/acceleration,
  with/without weights; OLS, IV and PPML propagate r(1) and can fit again. The
  generator is `validation/support/hdfe_cancel_plugin.py`; production code has no
  scalar-controlled injection hook. Interactive keyboard/UI Break latency was not
  measured.
- The `hard5` sample, coefficients and substantive covariance entries already
  agreed; its near-zero intercept/slope covariance exceeded the existing roundoff
  bound. Compensated means, residuals and covariance transformation resolve this
  without relaxing tolerances. An independent analytic covariance fixture also
  rejects the previous transformation arithmetic.
- Full licensed suites pass: **cppmlhdfe 1,302/1,302**, **chdfe 130/130**, and
  **creghdfe 385/385**. The current PPML suite includes concurrent test maintenance
  in this shared checkout; those edits were preserved, not attributed to this batch.
- Evidence: `temp/astra-nov11-batch5/ppml-final.log`, `chdfe.log`, `creghdfe.log`,
  `rollback.log`, `native-ubsan.log`, `native-omp-ubsan.log`, `negative.log` and
  `negative-vce.log`. Final rebuild smoke checks also pass (`hard5-final.log` and
  the rerun `rollback.log`); seed, documentation, architecture and whitespace
  checks pass. The OpenMP harness uses the installed Intel runtime with a
  local compatibility definition of LLVM's no-op `__kmpc_dispatch_deinit` marker.
  These are temporary test binaries, not packaged release artifacts.

### Batch six verification

- `validation/native_newcommands_paths.py` adds a dynamic fixture with atomic
  allocation/store counters. Serial and actual OpenMP builds pass UBSan with
  recovery disabled. Large split/join tests assert that worker threads performed
  writes. Linear/cubic interpolation checks analytic values, duplicate anchors,
  missing and filtered groups, numeric/string/multiple grouping and exact missing
  counts. Split checks byte-exact UTF-8, overlapping delimiters, trimming, empty
  fields, limits, 2,045-byte inputs and restricted/unrestricted writes.
- Join oracles compare every ordered output field against brute-force expected
  matches, including duplicate keys, interleaved groups, extended-missing group
  keys, absent groups, reversed/unbounded intervals and unmatched rows. Large
  master and using sides exercise separate parallel paths. The existing checked
  offset calculation was extracted without changing its semantics, allowing
  independent `INT_MAX` accepted/rejected-boundary checks without huge output
  allocations.
- Mandatory allocations are swept on small fixtures. Early/middle/last stores
  and final macro failures propagate, release all tracked resources and permit
  subsequent calls. Negative controls fail when interpolation drops its macro
  error, split ignores store errors, or join omits its output-limit check.
- `validation/support/newcommands_publication.do` passes **16 licensed checks**
  with test-only serial SPI injection (`publication_fault.h`): no final split or
  interpolation variables appear on failure, joins restore original data/labels,
  and subsequent commands succeed. The full current crangejoin suite passes
  **102/102** against the temporary production candidate.
- Evidence is under `temp/astra-nov11-batch6/`: `native-serial.log`,
  `native-omp.log`, `negative-interpolate.log`, `negative-split.log`,
  `negative-join.log`, `publication.log` and `crangejoin.log`.

The original November 6 audit and its historical verification record follow.

## Review basis

This is a release-readiness plan, not an implementation change or a release certification. The review covers the working tree at the time of the audit, including its modified and untracked implementation, ado, documentation, and validation files. Base commit: `cccfe84c75f016c467e06c33b6c49b6abcbd038d`. At the time of the audit, the release metadata identified version `1.2.4`, date `20261103`, Stata minimum `14.1`, protocol `2`.

Twenty distinct audit agents reviewed the commands in batches, with at most three reviewers active alongside the coordinating agent because the session supports four concurrent agents. The coordinator checked cross-command contracts, release checks, and the consolidated evidence. Existing source changes were preserved; this audit adds this document only.

Evidence labels distinguish **reproduced** behavior, **source-confirmed** defects with a traced execution path, and **coverage gaps** that require validation rather than asserting a wrong result. Bundled platform plugins are dated November 4 and predate some reviewed source changes. Plugin-dependent execution against those files cannot certify the current source. Ado-only reproductions remain valid where explicitly identified.

Priority: **P1** means silent wrong results, corruption, crashes, or a major release blocker; **P2** means a functional defect or material validation gap; **P3** means a lower-impact compatibility/documentation issue. No P0 claim is made without evidence. Each item below states the work and its completion criteria.

## Release assessment and command coverage

**Recommendation: hold release until the seven P1 items are repaired and validated.** There are **49 action items: 7 P1, 36 P2, and 6 P3**. This count includes test defects and coverage/release-evidence gaps; it is not a count of 49 independently reproduced product bugs. The P1 items are ASTRA-01, 03, 06, 10, 16, 17, and 34. They concern row corruption, a parser crash, incorrect sample selection, or loss of existing data columns.

All 22 public data/estimation commands and the `ctools` management command received an assigned review. Reviewers inspected the wrapper, native implementation, help/documentation, and relevant existing tests. Shared helpers were reviewed with their callers.

| Reviewer | Command(s) / additional scope | Resulting items |
|---|---|---|
| 01 | `csort`; independent second pass of P1 findings | 01-02; verified critical evidence |
| 02 | `cmerge`, native fallback and metadata | 03-05 |
| 03 | `cimport`, including new Parquet reader | 06-08 |
| 04 | `cexport`, binary formats and SAS sidecars | 10-12 |
| 05 | `creghdfe`, prediction, HDFE and inference | 13-17 |
| 06 | `civreghdfe`, auxiliary equations and bootstrap | 16-18 |
| 07 | `cppmlhdfe`, prediction/estat and bootstrap | 16, 19 |
| 08 | `cqreg`, rank checks and solver controls | 20-21 |
| 09 | `cbinscatter`, new binning transformations | 22-23 |
| 10 | `cpsmatch`, matching/statistics/radius paths | 24-25 |
| 11 | `chdfe`, shared solver and DoF | 26-27 |
| 12 | `crangestat`, optimized and flexible statistics | 17, 28-30 |
| 13 | `crangejoin`, rollback, strL, grouped joins | 33; no new functional defect established |
| 14 | `cinterpolate`, linear/cubic/grouped paths | 31-32; no new numerical defect established |
| 15 | `cwinsor`, weighted cutoffs and recovery | 17, 34-35 |
| 16 | `cdestring`, parser, strL and publication | 36-38 |
| 17 | `cencode` and `cdecode`, label mapping/staging | 39-41; no new cdecode output defect established |
| 18 | `csplit`, tokenizer/cache and staged writes | 42-43 |
| 19 | `csample` and `cbsample`, grouping/RNG/sizing | 30, 44-46 |
| 20 | `cdistinct`, `ctools`, loader/protocol, `cgroup`, release gating | 47-49; no new cdistinct functional defect established |

The coordinator added ASTRA-09, confirmed the shared undo-helper part of ASTRA-17, ran repository checks, and consolidated overlapping findings. Documented RNG differences, new-variable requirements for saved residuals/FEs, and intentional strL fallbacks were not automatically classified as bugs. Two investigated suspicions—temporary-frame leakage and early rangejoin schema widening—were excluded after checking Stata behavior/current source.

## Findings and proposed work

### ASTRA-01 — csort streaming transfers discard SPI failures

- **Priority / evidence:** P1; source-confirmed, native fault injection not executed.
- **Location:** `src/csort/csort_stream.c:375-403`, `487-488`, `536-539`, `562-572`, `599-600`, `647-676`, `968-969`.
- **Failure:** `stream_permute_string_var_flat()` and `stream_plan_run()` return `void` and discard numeric/string SPI read and write statuses. `csort key, stream(1)` can return success after failed non-key transfers. Failed reads can introduce empty strings or uninitialized numeric values; failed stores can leave columns partly permuted and rows misaligned. Automatic streaming can reach the same path on large/wide data.
- **Proposed solution:** return and propagate `stata_retcode` from both helpers and entry points; capture the first worker error safely. Preflight reads or retain recovery data and define rollback after writes begin. A nonzero result alone is insufficient if row alignment is already lost.
- **Acceptance:** inject numeric/string read and store failures at first/middle/last rows in serial and OpenMP paths, with keys and mixed payloads. Require a nonzero status, no uninitialized reads, and restoration of the original values, observation order, and metadata. If streaming requires a different failure policy, explicitly decide/document it before release and test it; silent row misalignment is unacceptable. Existing ordered-input tests do not exercise these non-key failure paths.

### ASTRA-02 — csort shortcut omits results and ignores nosortedby

- **Priority / evidence:** P2; reproduced in current ado, independently of the plugin.
- **Location:** `build/csort.ado:95-100`, success result posting at `380-394`; `build/csort.sthlp:94-100`, `217-232`; contrary existing assertion at `validation/stata/csort.do:1621-1632`.
- **Failure:** if `: sortedby` already equals the requested key, the wrapper exits before clearing metadata or posting documented results. With three rows, `sort x` followed by `csort x, nosortedby` leaves the marker `x`, `r(N)=.`, and `r(cmd)` empty.
- **Proposed solution:** route every successful shortcut through a shared result-posting block; clear the sorted marker when requested and define zero/missing timing values consistently.
- **Acceptance:** test repeat calls and pre-marked numeric, fixed-string, and strL data with/without `nosortedby` and verbosity. Assert marker state and all documented `r()` fields. Correct the existing strL assertion that expects the marker to survive.

### ASTRA-03 — cmerge can commit after failed keepusing or merge-marker transfers

- **Priority / evidence:** P1; source-confirmed, fault injection pending.
- **Location:** `src/cmerge/cmerge_io.c:46`, `72`, `98`, `119`, `124`; `src/cmerge/cmerge_impl.c:1492`, `1529`, `1564`.
- **Failure:** numeric/string reads and writes discard SPI status; workers unconditionally set `success=1`, their results are not inspected, and final `_merge` stores are unchecked. A failure can publish uninitialized values or partially written data while returning success and committing the undo copy.
- **Proposed solution:** store first error per worker, stop failed work, inspect every joined result, and pass failure through `cmerge_execute_failed()` so rollback runs. Check sequential fallback and `_merge` stores too.
- **Acceptance:** inject read/write failures in numeric/string keepusing data and final marker stores, at early/middle/last rows with one/multiple workers. Require exact original-master restoration in default and snapshot modes. Existing `fault_store` coverage targets the earlier master-store stage, not these callbacks.

### ASTRA-04 — cmerge strL fallback returns inconsistent category counts

- **Priority / evidence:** P2; reproduced using current ado in StataNow 19.
- **Location:** `build/_cmerge_native.ado:45-51`; main-path contract at `build/cmerge.ado:1653-1658`; incomplete assertions in `validation/stata/s2_strl.do:646-695`.
- **Failure:** the fallback omits `r(N_1)` through `r(N_5)` under `nogenerate`, counts after `keep()` instead of before it, and omits documented missing-valued `r(time)`. For master IDs `{1,2}` and using IDs `{2,3}`, `keep(match)` reports unmatched counts of zero instead of one each.
- **Proposed solution:** collect native category counts before filtering, using a temporary marker even with `nogenerate`; apply output filtering/removal afterward and explicitly return `r(time)=.`.
- **Acceptance:** compare the full `r()` surface between fixed-string/plugin and strL/native paths for normal, `nogenerate`, `keep()`, `update`, and `replace` combinations. Current strL tests compare values/types only.

### ASTRA-05 — cmerge loses using-file label definitions on Stata 14/15

- **Priority / evidence:** P2; source-confirmed branch gap, old-version runtime validation pending.
- **Location:** partial load at `build/cmerge.ado:602-603`; explanation and frames-only reload at `1185-1205`.
- **Failure:** with `keepusing(x)`, only labels attached to retained variables initially load. The complete label reload exists only for the Stata 16+ frames branch. On supported Stata 14/15, definitions attached solely to omitted using variables never reach the master, unlike native merge.
- **Proposed solution:** extract/save all using-file label definitions before narrowing the data, or add a preserve/restore-compatible non-frame reload.
- **Acceptance:** compare `label dir` and label definitions with native merge on Stata 14/15 and current Stata, including selected, omitted, unattached, conflicting, and `nolabel` definitions.

### ASTRA-06 — malformed Parquet DELTA_BYTE_ARRAY data can return success before crashing

- **Priority / evidence:** P1; source-confirmed.
- **Location:** `src/io/parquet_decode.c:283-290`, consumer `emit()` at `713-763`.
- **Failure:** when the prefix delta decoder succeeds but reports fewer values than the page's non-null count, `delta_strings()` jumps to cleanup with `rc == 0`. It returns an empty decoded buffer as success; `emit()` subsequently dereferences the absent string/fixed-value buffer. A malformed imported file can crash the Stata process.
- **Proposed solution:** set a nonzero status on every prefix-count mismatch, matching the suffix-count branch, and validate decoded buffer/count invariants before emission.
- **Acceptance:** zero, one-fewer, and too-many prefix counts for BYTE_ARRAY and fixed-length byte arrays fail cleanly in import and describe, without a crash or changes to the prior dataset; valid delta-encoded pages retain their results.

### ASTRA-07 — Parquet observation guard runs after allocation and decoding

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `build/_cimport_parquet.ado:82-89`; `src/io/cparquet.c:226-296`; allocation at `src/io/parquet_decode.c:1152-1176`.
- **Failure:** `cparquet scan` runs before the wrapper checks the requested range against the SPI limit of 2,147,483,647 observations. Oversized advertised ranges can allocate/decode huge arrays or fail with memory errors before the intended r(920) rejection.
- **Proposed solution:** compute and validate the prepared row range before scanning; enforce the bound independently in the C entry point.
- **Acceptance:** a small synthetic file advertising an oversized range is rejected before column allocation, retaining the original dataset in default and snapshot modes. Test permitted range endpoints as well as the first rejected value.

### ASTRA-08 — Parquet lacks an independent native malformed-input test suite

- **Priority / evidence:** P2; coverage gap, not a claim that every parser path is defective.
- **Location:** `validation/stata/s2_parquet.do:145-162`; `validation/native/import.suite.c`; new `src/io/parquet_*.c` modules.
- **Gap:** current malformed fixtures cover truncation, schema/Arrow metadata, and corrupt Snappy, but not the prefix-count failure in ASTRA-06. The native import suite does not call `pq_open`, `pq_read_column`, or the Parquet decoder. StataNow comparisons alone cannot provide sanitizer coverage or validate the claimed older-Stata support.
- **Proposed solution:** add native parser/decoder tests and sanitizer execution, with fixtures from an independent producer alongside the repository's writer. Cover footer/page sizes, inconsistent delta streams, dictionary indexes, levels, allocation bounds, and row-group totals.
- **Acceptance:** ASTRA-06 fails a test without licensed Stata; malformed inputs produce bounded allocation and explicit errors. Valid coverage spans codecs, v1/v2 pages, nulls, row ranges, and dictionaries.

### ASTRA-09 — validation runner self-test assumes a POSIX launcher on Windows

- **Priority / evidence:** P2; reproduced coordinator check failure and source-confirmed.
- **Location:** `validation/check_runner.py:443-458`; Windows executable requirement at `validation/stata.py:42-50`.
- **Failure:** `validation.check_runner.test_runner()` mocks `subprocess.run` but assumes the command starts with `/bin/zsh -lic`. On this Windows host without `CTOOLS_STATA_EXECUTABLE`, execution fails before the mock with `ValueError: Set CTOOLS_STATA_EXECUTABLE...`. Supplying that variable selects the Windows launcher and still contradicts the mock's POSIX assertion. This is separate from the 10 unittest-discovered tests, which pass.
- **Proposed solution:** mock the launch platform/environment deliberately for the POSIX case and add a Windows case for `/e do`, executable paths with spaces, nonce verification, and cleanup. Do not require a real Stata installation for a mocked runner contract.
- **Acceptance:** the full `test_runner()` contract passes on Windows, Linux, and macOS with and without a caller-provided executable setting. Both launcher forms are checked by negative controls.

### ASTRA-10 — binary exports evaluate if after replacing incoming r() results

- **Priority / evidence:** P1; source-confirmed.
- **Location:** `build/_cio_export.ado:116-117`, `172-174`; correctly frozen sample patterns at `build/cexport.ado:56-62`, `439-445`; current coverage at `validation/stata/results_scope.do:83-98`.
- **Failure:** `label dir` replaces the caller's `r()` before `marksample` evaluates the original condition. SPSS, SAS, XPORT5/8, and dBase exports share this helper. After `summarize w`, an export with `if w < r(mean)` can export every finite row because the missing mean becomes `.`, while `w > r(mean)` can export none.
- **Proposed solution:** freeze `if`/`in` selection immediately after syntax parsing and before metadata commands; use the stable marker for all conversion and export work.
- **Acceptance:** extend incoming-results scope tests to every binary exporter. Compare row IDs and exported counts with explicitly frozen expectations for `<`, `>`, `<=`, and `==` against incoming `r(mean)`/`r(N)`, including labels, dates, and `in` restrictions.

### ASTRA-11 — SAS sidecars do not escape apostrophes

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/io/cio_xpf.c:171-174`, `218`.
- **Failure:** companion SAS code interpolates paths and value-label text directly into single-quoted literals. A label such as `O'Brien`, or an apostrophe in a directory/path, produces invalid code while export reports success. Affected modes include XPORT5 `vallabfile(sascode|both)` and XPORT8 `vallabfile`.
- **Proposed solution:** implement a checked SAS-literal writer that doubles embedded apostrophes; use it for all label and `libname` path literals, retaining grouped file publication and write-error propagation.
- **Acceptance:** check generated code and parsed round-trip text for apostrophes in labels/paths, empty labels, double quotes, semicolons, Unicode, and maximum-length labels, across both XPORT versions and sidecar modes.

### ASTRA-12 — CSV correctness tests accept matching return codes without checking contents

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** `validation/stata/cexport.do:337-352`, `470-480`; custom formatting in `src/cexport/cexport_format.c`.
- **Gap:** extended-missing and extreme-number cases only compare native/ctools return codes. A successful export that erases `.a`/`.z` distinctions or emits inaccurate numbers can pass.
- **Proposed solution:** add CSV cell/token comparisons, with explicit numeric round-trip tolerances where byte equality is intentionally unnecessary.
- **Acceptance:** assert all extended missing codes, system missing, positive/negative extremes, exponents, zero, subnormals, and float/double cases. Exercise buffered, parallel, and mmap formatting paths where applicable.

### ASTRA-13 — creghdfe information criteria undercount fitted parameters

- **Priority / evidence:** P2; source-confirmed against installed reference and Stata ado code.
- **Location:** `build/creghdfe.ado:1060`; execution-only test at `validation/stata/creghdfe.do:3840-3842`.
- **Failure:** the wrapper posts `e(rank)=e(df_m)` without an estat handler or information-criterion df override. Stata therefore uses that rank, while the installed reghdfe estat handler requests `e(df_m)+1`. For the ordinary model with a reported constant, AIC is too small by 2 and BIC by `ln(e(N))`, even if likelihoods agree.
- **Proposed solution:** package a creghdfe estat handler implementing the reference parameter-count rule, or a verified supported-version df override. Do not change covariance rank solely to repair information criteria.
- **Acceptance:** compare the complete `estat ic` `r(S)` matrix, including df/AIC/BIC, with reghdfe for weighted, noabsorb, noconstant, collinear, saturated, and ordinary models. Checking only `_rc==0` is insufficient.

### ASTRA-14 — singleton regression tests turn unexpected estimation failures into passes

- **Priority / evidence:** P2; source-confirmed test defect.
- **Location:** `validation/stata/creghdfe.do:1698`, `1702`, `1727`, `1731`.
- **Gap:** after reghdfe successfully estimates the fixture, a creghdfe error calls `test_pass` as a graceful failure. Reference errors are also treated as passing. These cases can remain green after singleton filtering/sample/fit regressions.
- **Proposed solution:** require success from both commands for intentionally estimable fixtures, and isolate expected failures in separate tests with explicit statuses.
- **Acceptance:** a negative control that forces candidate failure after reference success must fail the case and runner. Retain exact sample, singleton count, observation count, and coefficient comparisons.

### ASTRA-15 — creghdfe quiet mode leaks new covariance warnings

- **Priority / evidence:** P3; source-confirmed.
- **Location:** `build/creghdfe.ado:82`, `1156-1162`.
- **Failure:** Driscoll-Kraay bandwidth and multiway PSD-correction warnings lack the `!__quiet` guard used by other notes, so `verbose(-1)` still emits output.
- **Proposed solution:** apply the same verbosity policy to these warnings.
- **Acceptance:** capture default and `verbose(-1)` logs for fixtures that produce each warning; assert default visibility and quiet suppression.

### ASTRA-16 — bootstrap scratch cleanup deletes existing user variables

- **Priority / evidence:** P1; source-confirmed across all three HDFE estimators.
- **Location:** pairs helper `build/_ctools_hdfe_boot.ado:268-269`, `298`, `309`; wild helpers `build/creghdfe.ado:1450`, `1465`, `1477`; `build/civreghdfe.ado:3099`, `3114`, `3121`; `build/cppmlhdfe.ado:1502`, `1517`, `1524`.
- **Failure:** fixed scratch stubs `__ctools_bsu` and `__ctools_wbv` are cleared with wildcard `capture drop ...*`, without establishing ownership. Default pairs bootstrap uses `draws(stata)` and drops all matching user columns. Wild cleanup runs even on enumeration paths. Successful estimation can silently delete unrelated user data; cppmlhdfe/civreghdfe do not have an enclosing full-data snapshot, and creghdfe's default undo assumes existing columns are untouched.
- **Proposed solution:** generate unique scratch names with tempvars/tempname-derived stubs and track the exact owned variable list. Drop only owned variables on normal completion, error, and interruption. Do not reserve a new user namespace as a substitute for ownership.
- **Acceptance:** for each estimator, create existing numeric/string columns named `__ctools_bsu1`, `__ctools_bsu_user`, `__ctools_wbv1`, and `__ctools_wbv_user`, with labels/characteristics. Assert their exact survival after successful pairs/wild runs, enumerated wild runs, failure, and interruption; verify cleanup leaves no newly owned scratch columns.

### ASTRA-17 — estimation/statistic wrappers evaluate incoming-r() conditions too late

- **Priority / evidence:** P1; source-confirmed; both underlying `r()` losses were reproduced in StataNow 19.5 without the plugin.
- **Location:** `build/civreghdfe.ado:345-350`; `build/creghdfe.ado:21`, `135`; `build/cwinsor.ado:18`, `154`; `build/crangestat.ado:14`, `392`; `_ctools_undo_record` in `build/_ctools_undo.ado:35-36`; missing command coverage in `validation/stata/results_scope.do`.
- **Failure:** IV `tsrevar` replaces incoming `r()` before `marksample`. Default undo setup in creghdfe, cwinsor, and crangestat also calls `label dir` before sample evaluation. After `summarize w`, `... if w<r(mean)` can select the wrong rows. With `w=1,...,5`, the intended two-row condition becomes true for all five after the mean becomes missing; the opposite inequality can exclude everything. The winsorization case can change values on unintended observations.
- **Proposed solution:** preserve incoming r-class state across shared transaction bookkeeping and freeze `if`/`in` immediately after parsing, before r-class operations; apply model-variable missingness afterward. Fixing the IV `tsrevar` order alone does not fix the shared undo-helper path.
- **Acceptance:** compare sample, count, and output values/coefficients with a precomputed marker for all inequality directions, weights, `in`, and time-series dependent variables. Exercise default/snapshot/fast modes for each undo-helper caller that accepts `if`. Retained probes: `temp/astra-nov6/audit06_result.log` and `temp/astra-nov6/root_sample_result.log`.

### ASTRA-18 — obsolete IV saving blocks overwrite valid auxiliary equations

- **Priority / evidence:** P2; source-confirmed.
- **Location:** actual auxiliary storage at `build/civreghdfe.ado:2700-2773`; obsolete structural storage at `2801-2821`; misleading test at `validation/stata/civreghdfe.do:458-470`.
- **Failure:** `savefirst` correctly stores each first stage, then stores the structural estimate as `<prefix>main`. An endogenous variable named `main` has its first-stage equation overwritten. The same happens to a reduced form when the dependent variable is `main` and `saverf` is requested. The late raw name also bypasses prior `strtoname()` handling for long prefixes.
- **Proposed solution:** remove obsolete `<prefix>main` storage/messages and retain the actual auxiliary equation workflow.
- **Acceptance:** restore each saved equation and compare dependent variable, model, coefficients, covariance, and sample with the reference auxiliary fit. Include `main`, long prefixes, repeated calls, replay/drop options. Replace the test that treats restoration of the structural copy as evidence of first-stage saving.

### ASTRA-19 — PPML invalid-input tests accept unrelated failures

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** `validation/stata/cppmlhdfe.do:926-928`, `955-957`, `968-970`; stronger convergence assertion at `899-901`.
- **Gap:** negative-outcome, empty-sample, and negative-exposure cases accept any nonzero status. A plugin-load, memory, or parser failure can pass without reaching the validation being tested. Existing automatic outputs and prior estimates are not checked for preservation.
- **Proposed solution:** assert the intended error code/message and data/output preservation; place nearby valid-input controls alongside the rejection cases.
- **Acceptance:** substituting an unrelated error code must fail each test. Valid neighboring inputs must fit, while invalid inputs preserve existing FE/d outputs and follow the documented estimation-state policy.

### ASTRA-20 — cqreg checks sample size before removing redundant columns

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `build/cqreg.ado:251-254`; later rank screen/compaction at `src/cqreg/cqreg_regress.c:1091-1200`.
- **Failure:** the ado guard compares physical sample size with the raw expanded regressor count. An estimable model with ten rows and ten duplicate regressor columns exits r(2001) before the C collinearity screen can omit redundant columns. Native qreg screens collinearity before fitting.
- **Proposed solution:** retain a minimal empty-sample guard in ado and make the definitive sample/df check after effective-rank determination, either in C or a compatible ado screen.
- **Acceptance:** compare qreg/cqreg with raw column count below/equal to/above N but low effective rank: duplicate, constant, linearly dependent, factor-expanded, and weighted designs. Check omissions, rank, sample, objective, covariance, and expected rejection of genuinely saturated models. Existing small-sample and collinearity tests do not combine these conditions.

### ASTRA-21 — cqreg silently replaces invalid solver settings

- **Priority / evidence:** P3; source-confirmed.
- **Location:** `build/cqreg.ado:43-53`; `src/cqreg/cqreg_fn.c:1243-1247`; `build/cqreg.sthlp:206`; valid-only tolerance tests at `validation/stata/cqreg.do:800-807`.
- **Failure:** `tolerance(0)`, negative tolerances, or values above 1 are silently replaced with `1e-12`. `nopreprocess(-2)` is accepted even though the documented range is nonnegative; only `-1` is rejected.
- **Proposed solution:** define and enforce the public tolerance range at both ado/C boundaries and reject every invalid negative preprocessing value with r(198), retaining the special removed-solver explanation for `-1` if useful.
- **Acceptance:** exercise zero, negative, missing, excessive, and valid boundary values. Invalid options fail before estimation with the documented estimation-state behavior; valid requested settings are honored.

### ASTRA-22 — cbinscatter narrows valid large xq() labels to 32-bit integers

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `build/cbinscatter.ado:179`; `src/cbinscatter/cbinscatter_impl.c:1314-1317`; `src/cbinscatter/cbinscatter_types.h:77`; promised labels at `build/cbinscatter.sthlp:195-201`.
- **Failure:** positive integer doubles are accepted as external bin labels, but the result path casts them to `ST_int`. Repeated labels such as 2,147,483,648 and 2,147,483,649 can group correctly yet yield undefined out-of-range C conversion and wrong IDs in `e(bindata)`/saved data.
- **Proposed solution:** keep dense integer indexes internally and retain external labels as doubles. If an upper restriction is intended, reject it explicitly before C and document it instead.
- **Acceptance:** verify exact IDs, means, and counts for `1`, `INT_MAX`, `INT_MAX+1`, `2^40`, and large exactly represented integers, including `by()`, medians, controls, and saved data; enable float-cast-overflow sanitizer coverage.

### ASTRA-23 — new binscatter transformations are absent from native sanitizer coverage

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** `validation/native/host.suite.c:717-764`; `validation/native/estimators.suite.c:838`, `998-1010`; valid-key tests at `validation/stata/s2_bin.do:419-429`.
- **Gap:** native cases cover publication failures, discrete capacity, and residualization/solver paths, but not xq-label remapping, weighted medians, random-subsample quantiles, or few-value polynomial fitting. Stata comparison tests do not replace native conversion/allocation checks.
- **Proposed solution:** add focused native cases for these transformations using small independently calculated fixtures and allocation failure injection.
- **Acceptance:** cover weighted ties, empty selected subsamples, sparse/extreme keys, few-value fits, allocation failures, and serial/parallel execution. ASTRA-22 must fail before its repair.

### ASTRA-24 — radius matching loses small local outcome means after huge unrelated outcomes

- **Priority / evidence:** P2; source-confirmed; exact arithmetic implementation mirrored in Python, not a compiled reproduction.
- **Location:** `src/cpsmatch/cpsmatch_impl.c:533`, `544-545`; inline `dd_add_d` in `src/ctools_ols.h`.
- **Failure:** local means subtract whole-donor double-double prefix sums. With donor scores `{.1,.2,.8,.9}`, outcomes `{1e100,1e50,1,2}`, and target `.85` with radius `.06`, only the last two donors match. Both relevant prefixes become `(1e100,1e50)`, producing a matched outcome of zero instead of 1.5. Outcomes outside the radius therefore change the local result.
- **Proposed solution:** use local/block/tree range accumulation, or detect cancellation and recompute affected ranges directly with compensated sums.
- **Acceptance:** vary nonmatching outcomes over large positive/negative scales and require unchanged local means; compare direct local sums and equivalent uniform-kernel matching. Include singleton intervals, mixed signs, and ATT/ATE paths.

### ASTRA-25 — cpsmatch ignores required scalar-publication failures

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/cpsmatch/cpsmatch_impl.c:864-878`; consumers at `build/cpsmatch.ado:571-574`, `690-691`; always-successful scalar mock at `validation/native/commands.suite.c:164-284`.
- **Failure:** count/support-bound `SF_scal_save()` statuses are discarded. The wrapper does not clear these scalars between calls; failed publication after a previous successful run can report stale counts/bounds as current results while native rc stays zero. Fresh-session failures can instead surface later as missing-scalar errors.
- **Proposed solution:** check all required saves, stage/clear command-owned results consistently, and propagate failure, preferably before output-column publication. Reuse the checked statistics-publisher pattern.
- **Acceptance:** inject failure at each count/bounds save after a prior success and in a fresh session. Require an explicit failure with no stale successful result contract; cover every required scalar, not only the separately checked 17 statistics.

### ASTRA-26 — long shared HDFE solves cannot cooperatively honor Break

- **Priority / evidence:** P2; source-confirmed; interactive latency not measured.
- **Location:** `build/chdfe.ado:39-40`; `src/hdfe/hdfe_solver.c:360`, `397`, `573`; SPI cancellation mechanisms at `src/stplugin.h:198-199`.
- **Failure:** `maxiterations(.)` maps to 2,147,483,646, but iterative loops do not inspect a cancellation flag or provide calling-thread SPI polling. A slow/stagnating transformation can keep Stata occupied until the C solve returns after the user requests Break. Other consumers of the shared solver need the same cancellation review.
- **Proposed solution:** poll cancellation on the calling thread, propagate it safely to workers, and inspect it at bounded solver intervals. Unwind with cancellation status while retaining staged-output guarantees.
- **Acceptance:** cancel before/during each solver family, including weighted/parallel cases, and assert bounded termination and cleanup. In licensed Stata, interrupt `generate()` and `clear` transformations and verify original-data restoration and no partial output publication.

### ASTRA-27 — chdfe native slopes and detailed DoF metadata are under-asserted

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** `validation/native/count_hdfe.c:40-49`, `214-225`; `validation/stata/chdfe.do:26-29`, `43-46`, `126-128`.
- **Gap:** native slope cases assert the sample marker but skip transformed-value checks. Mock metadata capture omits `_chdfe_*` results, and Stata parity cases omit documented per-effect `e(Kj)`, `e(Mj)`, exact/nested flags, owner mapping, and `e(N_hdfe_extended)`. Stata slope-value checks do not replace independent native math assertions.
- **Proposed solution:** add analytically known weighted/unweighted slope fixtures, orthogonality/idempotence checks, complete mock metadata capture, and full reference DoF comparisons.
- **Acceptance:** cover slope-only/intercept-plus-slope factors, multiple/constant/zero slopes, nested/disconnected categories, singletons, scaling, and pool sizes. Mutating a slope denominator or metadata field must fail the corresponding check.

### ASTRA-28 — crangestat flexible-only verbose calls fail or display stale timings

- **Priority / evidence:** P2; source-confirmed.
- **Location:** native-call condition at `build/crangestat.ado:479`; unconditional verbose scalar reads at `574-578`; flexible-only tests at `validation/stata/regressions.do:2576-2582`.
- **Failure:** flexible-only covariance, correlation, regression, or user-Mata requests do not run the C plugin and create no native timing scalars. `verbose` nevertheless reads them. A fresh session can error after successful computation and roll back the outputs; a prior native call can leave stale times that are printed instead.
- **Proposed solution:** initialize invocation-local timing state, display native fields only when native work ran, and separately time flexible/Mata computation.
- **Acceptance:** run every flexible statistic alone/mixed with native statistics, in fresh sessions and after other calls, with and without verbosity. Outputs match the reference and logs contain only current-invocation timings.

### ASTRA-29 — constant extended-missing interval bounds are rejected inconsistently

- **Priority / evidence:** P2; source-confirmed compatibility defect against installed rangestat source.
- **Location:** `build/crangestat.ado:90`, `108-130`; `src/crangestat/crangestat_impl.c:38-39`; `src/ctools_parse.c:133-144`.
- **Failure:** exact `.` constants are normalized to unbounded, but `.a` through `.z` pass ado numeric validation and are rejected by C. Equivalent variable-valued missing bounds are already normalized, and rangestat treats all numeric missing bounds as unbounded.
- **Proposed solution:** normalize all numeric missing constants to the same unbounded representation before constructing plugin arguments.
- **Acceptance:** compare `.`, `.a`, and `.z` in each low/high position with equivalent bound variables and reference results. Include small/large groups, repeated keys, prefix-query dispatch, and `excludeself`.

### ASTRA-30 — crangestat and csample silently ignore negative threads()

- **Priority / evidence:** P3; source-confirmed.
- **Location:** `build/crangestat.ado:464-466`; `build/crangestat.sthlp:126-127`; `build/csample.ado:153-156`; contrasting explicit check at `build/cbsample.ado:23`.
- **Failure:** both wrappers forward only positive values, so negative values silently choose the default while excessively large positive values are rejected. The crangestat help states 1-256, and the zero/default policy is unclear. cbsample already rejects values outside 0-256 explicitly.
- **Proposed solution:** validate a single documented range in ado/C, explicitly documenting zero if retained as the default sentinel.
- **Acceptance:** test negative values, zero, 1, 256, and 257 with exact statuses and unchanged data on rejection.

### ASTRA-31 — native interpolation tests miss grouped, parallel, and publication failures

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** `validation/native/commands.suite.c:549-551`, `596-612`; parallel thresholds at `src/cinterpolate/cinterpolate_impl.c:101`, `115-116`.
- **Gap:** fixtures exercise only tiny ungrouped linear/cubic cases. The available store-failure switch is never enabled and final macro publication always succeeds. Native tests do not reach grouped/string-key/parallel paths or publication failures, despite substantial Stata numerical coverage.
- **Proposed solution:** add independent grouped fixtures and enough observations to activate OpenMP; inject numeric-store failures and `_cinterpolate_nmissing` macro-save failures.
- **Acceptance:** cover duplicate anchors, empty/all-missing groups, numeric/string grouping, and serial/parallel equality. Early/middle/late store and final-publication errors propagate, release allocations, and never publish the requested wrapper output variable.

### ASTRA-32 — cinterpolate prints a plural message for one missing value

- **Priority / evidence:** P3; source-confirmed, presentation-only.
- **Location:** `build/cinterpolate.ado:89`.
- **Failure:** exactly one missing result prints `(1 missing values generated)` rather than the native singular form.
- **Proposed solution:** choose `value` for one and `values` otherwise.
- **Acceptance:** log assertions for zero/one/multiple missing results, including no message for zero. This item can be deferred independently of correctness work.

### ASTRA-33 — crangejoin native coverage omits grouped parallel joins and output limits

- **Priority / evidence:** P2; confirmed coverage gap; no new functional defect established for this command.
- **Location:** tiny single-thread native fixture at `validation/native/commands.suite.c:553`, `618-628`; grouped/limit paths at `src/crangejoin/crangejoin_impl.c:175-201`; Stata tests at `validation/stata/crangejoin.do:45-60`, `98-110`.
- **Gap:** grouped Stata tests are small, while the large parallel case omits `by()`. A large grouped benchmark exists at `validation/performance/commands.do:167-183`, outside native sanitizer correctness checks. Mixed-group search, cache reuse across thread boundaries, and cumulative output-limit handling lack native coverage.
- **Proposed solution:** add dynamic numeric/string/multiple-key fixtures with serial/OpenMP execution and count-only boundary tests that do not materialize billions of rows.
- **Acceptance:** compare ordered matches across thread counts for repeated/alternating/missing/absent groups, extended missings, duplicate using keys, reversed/unbounded intervals, and unmatched rows. Preparation must reject more than 2,147,483,647 output rows with r(920) before allocation/writing; validate the permitted boundary synthetically. Include allocation/store failures under sanitizers.

### ASTRA-34 — cwinsor default rollback loses overwritten existing output targets

- **Priority / evidence:** P1; reproduced with current ado; the failure hook fires before winsorization.
- **Location:** `build/cwinsor.ado:138`, `260-261`, `281-284`; incomplete rollback fixtures at `validation/stata/cwinsor.do:512-531`.
- **Failure:** existing targets selected with `generate()`/suffix/prefix plus `replace` are recast and overwritten before plugin work. Value undo is enabled only for the option-free in-place mode. With `x=1,...,5`, byte `x_w=99`, and `$CTOOLS_cwinsor_fault="outputs"`, `cwinsor x, suffix(_w) replace` returns r(1) but leaves `x_w=1,...,5` stored as double. Default recovery loses the original column.
- **Proposed solution:** stage replacement targets in temporary columns until success, or retain original values and complete metadata before any mutation. Self-target generation also needs the existing-target recovery policy.
- **Acceptance:** inject outputs/bounds/store/post failures with existing outputs, mixed new/existing targets, and self-target generation. Require exact restoration of values, types, labels, formats, characteristics, order, sorted/changed state. Probe: `temp/astra-nov6/audit15_result.log`.

### ASTRA-35 — weighted winsor cutoffs overflow for finite positive weights

- **Priority / evidence:** P2; source-confirmed, arithmetic mirrored independently.
- **Location:** `src/cwinsor/cwinsor_impl.c:187`, `218-224`; modest-weight fixtures at `validation/stata/s2_bin.do:484-550`.
- **Failure:** the percentile target uses `total * pctl / 100`. For values 1-20 with equal weights `1e306` and `cuts(50 50)`, the finite total `2e307` overflows when multiplied by 50. The search falls through to 20 instead of the median 10.5. Larger totals can overflow during summation too.
- **Proposed solution:** normalize weights by a safe common scale before accumulation and calculate fractional targets without overflowing intermediates. Preserve the specified tie/interpolation rule.
- **Acceptance:** cutoffs remain invariant when all weights are scaled by large/tiny positive constants. Cover near-overflow totals, disparate weights, ties, by-groups, trimming, and each supported weight type.

### ASTRA-36 — cdestring silently truncates ignore() at 255 bytes

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/cdestring/cdestring_impl.c:39`, `50`, `268`.
- **Failure:** a fixed 256-byte buffer silently stops decoding after 255 bytes, although the public option has no such limit. With 255 repeated `a` characters followed by `$` in the ignore list, `$12` stays nonnumeric instead of becoming 12. Unicode can be cut mid-sequence, and the long-strL fallback processes the full list, making behavior path-dependent.
- **Proposed solution:** dynamically store the full decoded ignore string; if a hard bound is unavoidable, reject oversized input explicitly before conversion.
- **Acceptance:** compare native/plugin/strL-fallback results at 254/255/256 bytes and beyond, requiring a character after the boundary. Include repeated characters, Unicode spanning the boundary, escaped whitespace, and backslashes.

### ASTRA-37 — empty verbose destring calls can fail after publishing outputs

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/cdestring/cdestring_impl.c:443-451`; `build/cdestring.ado:507-512`; nonverbose-only empty tests at `validation/stata/cdestring.do:1172-1195`.
- **Failure:** the zero-selected-row branch returns counts without timing scalars. Verbose display reads them unconditionally after output-name publication. A fresh `cdestring s if 0, generate(n) verbose` can error with a published output; subsequent calls can print stale times.
- **Proposed solution:** initialize complete invocation-local counts/timings on every success path, check publication status, and prepare fallible results before output-name commit. Display must not depend on stale command scalars.
- **Acceptance:** fresh/repeated calls on empty datasets and empty selections, generate/replace, verbose/nonverbose, all succeed with the expected missing outputs and current zero/missing timings.

### ASTRA-38 — destring conversion lacks native sanitizer and fault-injection coverage

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** no converter implementation harness in `validation/native/` or the native Python suites; implementation at `src/cdestring/cdestring_impl.c`.
- **Gap:** Stata command/option tests do not exercise C parsing/conversion under sanitizers or injected failures.
- **Proposed solution:** add native ignore-wire, UTF-8, decimal/exponent, float-boundary, percentage, empty-sample, allocation/store, and scalar/macro publication tests.
- **Acceptance:** independent expected values catch ASTRA-36/37 without licensed Stata; failures propagate and release resources. Wrapper-level cases verify no partially published outputs.

### ASTRA-39 — cencode proceeds without a successful label backup

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `build/cencode.ado:14-19`, backup implementation at `444-460`; atomicity promise at `build/cencode.sthlp:105`; later-failure test at `validation/stata/regressions.do:1935-1967`.
- **Failure:** backup errors are captured and recorded as `__saved=0`, but encoding still runs. If labels are extended and a later variable fails, restoration is skipped because no complete backup exists. Memory/Mata failure during backup can thus disable the documented label rollback guarantee.
- **Proposed solution:** abort and propagate the backup error before entering mutation-capable code; clean any incomplete backup state without touching labels.
- **Acceptance:** inject errors before/during backup and verify encoding never begins, all original labels/variables/metadata remain intact, and the intended error survives. Retain later-variable rollback tests as a separate stage.

### ASTRA-40 — cencode can install incorrect labels after failed result-scalar saves

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/cencode/cencode_impl.c:392-396`, `514-516`, `604-606`, `731-733`; ado consumption at `build/cencode.ado:226`, `246`, `266-291`.
- **Failure:** required scalar publication statuses are ignored. When an empty first variable sets `_cencode_n_unique=0` and saving the second variable's new unique count fails, the wrapper can read stale zero, skip the valid label file, and attach an empty label definition to nonempty numeric codes. Other required-save failures have similar stale/missing-result risks.
- **Proposed solution:** check every required result save on every native return path; initialize/stage per-variable results and permit label installation only after complete publication succeeds.
- **Acceptance:** fail each required scalar save in fresh/repeated calls and successive-variable conversions, including empty-then-nonempty input and existing-label extension. Require explicit failure and unchanged original labels/variables, with no published output or stale success.

### ASTRA-41 — native label-file tests do not cover cencode/cdecode implementations

- **Priority / evidence:** P2; confirmed coverage gap; no new cdecode output defect established.
- **Location:** `validation/native_contracts.py:275-277`; label-writer-only case at `validation/native/commands.suite.c:286-336`.
- **Gap:** `test_encode_labels_native()` exercises the label-file writer, not `cencode_main`, `cdecode_main`, or their streaming/mapping implementations. Extensive Stata coverage does not inject failures or sanitize these C paths.
- **Proposed solution:** add native SPI harnesses for streaming encoding, existing-label mapping, and v2 decoding with controllable reads/stores/scalar/macro publication and independent expected mappings.
- **Acceptance:** cover sparse/negative/duplicate-text label codes, empty selections, extended missings, Unicode, malformed mapping protocols, lookup/binary-search dispatch, chunks, threading, and allocation/transfer/publication failures. ASTRA-40 must fail before repair.

### ASTRA-42 — split has an undocumented 2,045-delimiter ceiling

- **Priority / evidence:** P3; current-source guard confirmed; behavior reproduced with current ado and the older bundled plugin.
- **Location:** `src/csplit/csplit_impl.c:87-91`; unrestricted wrapper parsing at `build/csplit.ado:40-53`.
- **Failure:** more than 2,045 delimiter entries is rejected with r(198), although native split accepts the supplied list and current help states no ceiling. A one-row `a:b` input with 2,046 repeated `:` delimiters succeeds natively and fails in csplit.
- **Proposed solution:** allocate delimiter metadata dynamically or deduplicate while preserving precedence. Explicitly documenting/validating a retained ceiling is a compatibility decision, not full parity.
- **Acceptance:** test 2,044/2,045/2,046 entries, repeated/distinct/overlapping alternatives, native matching precedence, atomic rejection, and a successful call after failure. Probe: `temp/astra-nov6/audit18_result.log`.

### ASTRA-43 — native split tests omit publication failures and OpenMP

- **Priority / evidence:** P2; confirmed coverage gap.
- **Location:** `validation/native/commands.suite.c:550-551`, `595`, `613-617`; `validation/native_contracts.py:868-908`; parallel thresholds at `src/csplit/csplit_impl.c:121`, `192`.
- **Gap:** native fixtures use tiny successful writes/allocation failures; string-store failure support is not activated, width-macro publication always succeeds, and token tests compile without OpenMP. Stata parallel comparisons at `validation/stata/contracts.do:192-209` do not supply native sanitizers/fault injection.
- **Proposed solution:** add early/middle/late string-store and width-macro publication failures, plus large OpenMP single-/multiple-delimiter fixtures.
- **Acceptance:** errors propagate, caches release, subsequent calls work, and failed staged writes publish no final output names. Serial/parallel results agree for trimming, overlapping delimiters, empty fields, limits, UTF-8 bytes, filters, and maximum-width input.

### ASTRA-44 — cbsample evaluates subscripted size expressions in the wrong row order

- **Priority / evidence:** P2; reproduced using the current ado helper and native bsample, without the plugin.
- **Location:** `build/cbsample.ado:342-348`; incomplete expression cases at `validation/stata/cbsample.do:583-592`, `786-813`.
- **Failure:** unstratified size expressions are evaluated in original order; native bsample first stably moves eligible rows ahead of excluded rows. For `x=1,...,30`, selecting `x>10` and requesting `x[1]` yields one draw from the current helper versus eleven observations from native bsample.
- **Proposed solution:** evaluate expressions in native eligible-observation order, retain an original-row map, and restore working order after sizing.
- **Acceptance:** compare sizes/statuses for `x[1]`, other subscripts, `_n`, and `_N`, with `if`/`in`, clustering, weight outputs, and stratification. Probe: `temp/astra-nov6/audit19_result.log`.

### ASTRA-45 — empty csample calls leave stale documented count scalars

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `build/csample.ado:68-69`; documented scalars at `build/csample.sthlp:127-137`; return-code-only empty tests at `validation/stata/csample.do:398-407`.
- **Failure:** zero eligible rows return success without updating `_csample_nobs`, `_csample_kept`, `_csample_dropped`, `_csample_ngroups`, or `_csample_threads`. Results are absent in a fresh session or retain unrelated prior counts.
- **Proposed solution:** define and initialize the complete result surface on every success path, including empty selections.
- **Acceptance:** follow successful draws with `if 0` and empty-data calls, with/without grouping/verbosity. Eligible/kept/dropped counts are zero and timing/thread/group fields follow a documented current-invocation policy.

### ASTRA-46 — samplers ignore scalar-publication failures before destructive ado steps

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/csample/csample_impl.c:458-471`; `src/cbsample/cbsample_impl.c:742-743`, `778-790`; consumers at `build/csample.ado:214-215`, `build/cbsample.ado:267-280`.
- **Failure:** required result-scalar saves discard status. csample drops rows before reading possibly stale/missing counts. cbsample clears old scalars but `idcluster()` drops/expands observations before reading ID scalars, so failed publication can cause a late error after resampling.
- **Proposed solution:** check every required native publication and validate all required results before row deletion/expansion. Failure must occur before destructive output publication.
- **Acceptance:** inject every required save failure in fresh/repeated calls, assert explicit failure before resampling publication, and verify cleanup. The sampling mock at `validation/native/commands.suite.c:1207`, `1237` currently always succeeds at scalar saves.

### ASTRA-47 — shared grouping ignores count publication errors

- **Priority / evidence:** P2; source-confirmed.
- **Location:** `src/cgroup/cgroup_impl.c:173`; `build/_ctools_group.ado:99-103`.
- **Failure:** `_ctools_group_n` publication is unchecked. The wrapper reads it outside captured plugin-error handling after generating the output. Missing publication can raise r(111) and bypass cleanup, leaving a new column; an old scalar can instead provide a stale count. This adapter supports multiple public commands even though it is internal.
- **Proposed solution:** check required count publication, preferably before output writes, and include result consumption in cleanup/fallback handling.
- **Acceptance:** fail scalar publication with/without a previous scalar, across count/generate/label and strict/fallback modes. Require explicit failure or correct fallback, no stale count, and no partial output. Add a direct `cgroup_main` harness in addition to existing low-level grouping-utility tests.

### ASTRA-48 — Linux loader troubleshooting incorrectly requires libgomp

- **Priority / evidence:** P3; source-confirmed documentation/diagnostic defect.
- **Location:** `build/_ctools_load.ado:74`; current dependency contract at `docs/INSTALLATION.md:81-85`.
- **Failure:** the loader claims Linux requires `libgomp.so.1` from GCC 6+, while distribution plugins link OpenMP statically and dependency checks prohibit that external runtime. Users receive contradictory advice for unrelated load failures.
- **Proposed solution:** describe the current x86-64/glibc/static-runtime requirements, retain the actual load error, and provide the existing `ldd -r` troubleshooting route. Distinguish custom dynamic builds if mentioned.
- **Acceptance:** a controlled failed-load fixture emits diagnostics consistent with distribution dependency checks and does not claim external OpenMP is mandatory.

### ASTRA-49 — publication is not bound to licensed Stata validation of exact artifacts

- **Priority / evidence:** P2; confirmed release-evidence gap, not a claim that current CI checks failed.
- **Location:** native-only validation at `.github/workflows/build.yml:81-91`, package prerequisites at `368-379`, publication prerequisites at `426`, `470`; absence asserted at `validation/checks.py:72-73`; candidate-testing requirement at `docs/development/PLATFORMS.md:12`.
- **Gap:** builds, signatures, package integrity, and native checks can permit publication without recorded licensed Stata validation of the ado/plugin combination on supported systems. Optional revision checks are useful but do not require this evidence. The reviewed dirty tree also has no fresh matching plugin build established by this audit.
- **Proposed solution:** require licensed-runner results or a deliberate release approval backed by retained test records, bound to candidate source SHA, platform, final staged package hashes, Stata version, reference-package versions, and suite-completion results. A manual licensed validation process can satisfy this without public CI licenses.
- **Acceptance:** missing/failed/skipped/stale-SHA/wrong-platform/wrong-artifact evidence prevents publication. Fresh candidate artifacts pass command, oracle, and executable option suites on supported platforms, including minimum-supported Stata/system checks. Update the harness contract if a new validation job replaces the intentionally native-only workflow.

## Ordered remediation plan

1. **Establish failing regressions for the seven P1 items.** Start with the existing ado fault hooks and tiny standalone cases. Add native SPI-failure mocks for sort/merge and a malformed Parquet prefix fixture. Preserve every reproduction independently of older bundled binaries. Completion: each issue fails a meaningful test for the stated reason before its implementation is changed.
2. **Repair shared data-safety and sample contracts.** Fix scratch-variable ownership (16), preserve incoming r-class state and freeze samples (10/17), restore overwritten winsor outputs (34), propagate sort/merge transfer failures and rollback (01/03), and reject malformed Parquet prefixes (06). Review all callers of changed shared helpers. Completion: no existing user columns/rows/metadata change on protected failure paths, and sample/output parity matches frozen expected data.
3. **Correct the remaining functional issues.** Work in parallel by subsystem: merge metadata/results (04/05), Parquet bounds (07), estimator/postestimation behavior (13/18/20/21), binning/matching/weighted numerics (22/24/35), interval/sample/empty-input behavior (28/29/37/44/45), label/parser behavior (36/39/40), and required publication (25/46/47). Completion: the per-item acceptance criteria pass, including repeated calls and boundary inputs.
4. **Strengthen tests before accepting a green suite.** Repair false-positive assertions (12/14/19) and Windows runner mocks (09). Add missing native adapters and branch coverage (08/23/27/31/33/38/41/43), cancellation tests (26), and required scalar/macro failure injection. Use independent expected answers and negative controls so matching return codes or self-generated fixtures cannot hide incorrect output. Completion: intentional result/error mutations make the appropriate test fail.
5. **Align public validation and diagnostics.** Resolve or explicitly document the low-impact option/message issues (15/21/30/32/42/48). Regenerate inventories only after behavior/help/test evidence agree. Completion: documented options have consistent accepted ranges and diagnostics describe the actual release artifacts.
6. **Build and qualify one frozen release candidate.** Include all intended new files in the candidate source commit, build every platform from that commit, stage the signed packages, and run the artifact-bound validation described in ASTRA-49. Re-run affected focused tests first, then the full native/sanitizer and licensed Stata suites. Completion: archived evidence identifies the exact source and package hashes that will be published.

The work can be divided into wrapper/transaction, import/export, estimator/numerical, and validation/release workstreams. Shared helper changes should land before their callers' final acceptance runs. A failing test must not be made green by weakening its assertion, skipping the branch, or changing the reference result without a documented compatibility decision.

### Release exit criteria

- [ ] All seven P1 items are repaired and their regressions pass against freshly compiled code/current ado files.
- [ ] P2 functional result/data-safety issues are fixed; every remaining P2/P3 item has an explicit fix or documented deferral with scope and rationale. A coverage gap is not silently closed by an inventory check.
- [ ] Full native checks pass with AddressSanitizer/UndefinedBehaviorSanitizer where supported, plus float-cast-overflow coverage for narrowing conversions. Exercise serial/OpenMP and relevant allocation/transfer/publication failures.
- [ ] Licensed Stata command, oracle, bootstrap, incoming-results, and executable option cases pass with no skipped required reference checks. Record Stata and reference package versions.
- [ ] Candidate evidence covers Windows, Linux x86-64, macOS ARM, and macOS Intel; retain minimum-supported Stata/system results for the stated compatibility contract, especially Stata 14/15 non-frame branches.
- [ ] Source SHA, embedded plugin revision, platform, final package hashes, signatures, dependency checks, and full suite completion records agree. Validate the staged installable package in a clean Stata session rather than an unrelated checkout.
- [ ] Release metadata, generated command/option inventories, package manifests, help, licenses, and documentation checks pass on the final tree.

## Verification record

The coordinator ran the following against the reviewed tree:

| Check | Result |
|---|---|
| `scripts/sync_release.py --check` | Pass |
| `scripts/sync_architecture.py --check` | Pass |
| `validation/run.py --checks options` | Pass: 461 inventoried options; this is inventory evidence, not exhaustive interaction coverage |
| `validation/run.py --checks docs` | Pass: 36 files |
| `validation/run.py --checks seeds` | Pass, including negative-control self-tests |
| `validation/run.py --checks package build --metadata-only` | Pass |
| Public command/package/Stata-minimum/license/publication contracts | Pass: 22 public data/estimation commands plus the `ctools` management command reviewed separately |
| `python -m unittest validation.check_runner` | Pass: 10 tests |
| `validation.check_runner.options_main()` and `test_performance()` | Pass |
| `validation.check_runner.test_runner()` | Fails on Windows launcher setup; see ASTRA-09 |

Python used: `C:\Users\Mike\.local\bin\python3.12.exe`, with `PYTHONUTF8=1`. Tests needing scratch space used repository-local `temp/astra-nov6` because the sandbox system-temp path was inaccessible. This environmental issue is not classified as a product defect. No C compiler was found on PATH, so a full fresh native build/sanitizer run remains required.

Targeted Stata probes established the ado-only sort shortcut, native merge fallback results, underlying r-class state loss, winsor rollback, and bootstrap-size-expression findings. The split ceiling probe used the bundled plugin and is explicitly labeled accordingly. Mathematical mirrors support the two extreme-scale numerical findings; they do not constitute a fresh plugin execution. This focused audit did not establish complete numerical parity of every estimator/bootstrap configuration, and a command with no newly confirmed defect is not certified bug-free.
