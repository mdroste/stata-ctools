# C pre-release audit — October 14, 2026

This audit looked for release-blocking defects in the first-party C code under `src/` (vendored code in `src/vendor/` excluded): memory-safety errors, crashes, heap overflows, leaks, races, hangs, and silent data corruption caused by indexing errors. It found **69 distinct defects**. **10 are release-blocking**: 3 P0 and 7 P1. Of the rest, 33 are P2, 25 are P3 (18 of them latent, i.e. no current input reaches them), and 1 report was refuted. Each finding comes with a proposed fix. Full diffs are in the companion file [C_RELEASE_AUDIT_2026-10-14_patches.md](C_RELEASE_AUDIT_2026-10-14_patches.md).

**Status on 15 October:**
- Every P0, P1 and P2 finding is fixed in the working tree (uncommitted), along with two defects found while fixing; see [P0 fix status](#p0-fix-status) and [P1 and P2 fix status](#p1p2-fix-status).
- The 25 P3 findings remain open.

<a id="p0-fix-status"></a>
## P0 fix status (2026-10-14)

The three P0 fixes and the SPI-contract change are applied in the working tree. They are uncommitted, and `build/` has **not** been rebuilt. Every check below ran against a scratch plugin built from the patched `src/` with the repo-local libomp (StataNow MP, macOS arm64).

| ID | Change | Verified by |
|---|---|---|
| F009 (+ F019, F020, csample) | `#define SD_SAFEMODE` before `stplugin.h` in `src/csort/csort_stream.c`, `src/crangestat/crangestat_internal.h`, `src/crangejoin/crangejoin_impl.c` and `src/csample/csample_impl.c` | Post-`predict` `csort id, stream(2)` and `stream(4)`: `cf` passes (before the fix, 71 mismatches in `yhat`). crangestat and csample after `regress` give the same results as after `ereturn clear`. `validate_csort` 351/0, `validate_crangestat` 195/0, `validate_crangejoin` 93/0, `validate_csample` 28/0, native `test_group_windows_native.py` passes |
| F033 | `cio_do_export` keeps `stringneed` apart from `stringcap` (the real size of `s`) and reallocs only when the need exceeds the capacity (`src/io/cio_readstat.c`) | `cexport spss`, `sas` and `sasxport8` of auto plus str20 and two strL variables: rc 0 (before: SIGSEGV). The SPSS file round-trips through native `import spss` with the strL values intact. An empty strL exports. `validate_io_formats` 202/0, `validate_cexport` 213/0 |
| F004 | New `ctools_worker_attr_init()` (at least 8 MiB of stack) in `src/ctools_threads.c/h`, used by the persistent pool and by cimport's two raw-pthread fallbacks. `field_buf` in `cimport_parse_chunk_parallel` is now a heap allocation that sets `error_code` on failure | A 5.8 MB, 120,000-row CSV (multi-chunk) gives `cf` rc 0 against native `import delimited`, both with default threads and with `threads(1)`. `validate_cimport` 360/0. The Intel/Rosetta crash itself was not run: there is no x86 Stata here |

**Full gate (later on 14 October):**
- `build/ctools_mac_arm.plugin` was rebuilt from the patched tree: `make macos-arm` with the repo-local libomp, dependency contract checked, sha256 172a3329…
- `python3 validation/run_stata_audit.py` then passed in full: validate_all 5,151 passed / 0 failed / 0 skipped, io_formats 202/0, delimited options 227, Excel options 196 + 15, and option_coverage 1395/1395. The other extra suites (io_parity, cqreg_performance, bigdata, transport_strl, data_io, refactor_timing, coverage_helpers) also completed without failures.
- That is identical to the pre-rebuild baseline (plugin b1569f1d…).

**Still to do:**
- Rebuild the other three platform plugins. The x86 macOS plugin is from 26 September.
- Confirm F004 on an Intel Mac or under Rosetta.

**Noticed in passing (pre-existing, not fixed):**
- `cexport sasxport5` of plain `sysuse auto` exits r(198) with no message.
- `src/io/cexport_xls.c:372` raises a `-Wdangling-else` warning.

<a id="p1p2-fix-status"></a>
## P1 and P2 fix status (2026-10-15)

All 7 P1 and all 31 P2 findings are fixed in the working tree (uncommitted); F019 and F020 were already fixed with the P0s. Two more defects were found and fixed along the way:
- **ReadStat on Windows:** `src/vendor/readstat/readstat_iconv.h` included `../../cio_iconv_win.h`, a path that broke when the vendored tree moved from `src/io/vendor/` to `src/vendor/`. It now includes `../../io/cio_iconv_win.h`; MinGW could not compile ReadStat before the change.
- **cmerge `_merge` check (MRG-11):** the existence check matched abbreviations, so a master variable named `_merge_temp_filter` or `_merge_prev` made a default `cmerge` fail with r(110). It now passes `exact`.

**How the fixes were checked:**

| Check | Result |
|---|---|
| Plugin | `build/ctools_mac_arm.plugin` rebuilt (sha256 96b3e1c3…, dependency contract checked); `build/*.ado` synced from `package/ado` |
| `python3 validation/run_stata_audit.py` | validate_all 5,151 passed / 0 failed / 0 skipped; option_coverage 1395/1395; every extra suite passes. Identical to the baseline |
| `python3 validation/suite_registry.py --native` | 34 of 35 suites pass, including every ASan/UBSan and fault-injection suite. The one failure is `test_release_gate`'s pre-existing inventory assertion (`encode_labels`, `group_windows` and `winsor_transport` are registered but not public) |
| Per-fix Stata repro do-files | 75 checks, run against a scratch plugin built from the same tree. All pass except three, none of them a defect in the fix: (a) F016's space-in-tmpdir path is a SKIP (it needs Stata started with `STATATMP` pointing at such a directory); (b) the F002 check that uses `cf` runs out of variable room at 8,191 variables, but the same merge checked with `datasignature` matches native `merge`; (c) the F017 `_pctile` reference is miscalibrated: cwinsor's `by()` results are identical to the pre-fix plugin and to `winsor2` |
| Fault-injection and ASan harnesses (scratch) | F005, F006, F011 (every allocation in `csort_stream.c`), F024, F026, F028, F030 (differential fuzzing against the old code: byte-identical output), F031, F045–F050 (allocation-failure sweeps; 180 civreghdfe configurations bit-identical to the old code) |
| Windows-only fixes (F027, F029) | MinGW syntax checks only; not run on Windows |

Fixes that change user-visible behavior (all on paths that crashed, hung or gave wrong results before):

| ID | Change |
|---|---|
| F021 | crangestat: `(stat) new = x y` and `(stat) a b = x` exit r(103), `(stat) = x` exits r(102); the plugin rejects inconsistent argument lists with r(198) |
| F010 | csort: repeated keys are de-duplicated (sorting by `x x` sorts by `x`); an invalid key index in streaming mode is r(198) |
| F026 | cimport: `numericcols()`/`stringcols()` reach C as ado locals with their byte length. An unreadable list is r(198) instead of being silently dropped, and no `CIMPORT_*` global is left behind |
| F016 | cencode: the label-file path reaches C through a local with a length check. **The ado and the plugin must ship together.** |
| F036 | cexport dbase: negative fractional `%td` values floor, as the `%td` display does; dates before year 0 are r(108) |
| F031 | cexport delimited: an observation longer than 2 GiB of text is r(920) with a message |
| F042, F046 | civreghdfe: more than 46,340 regressors or instruments, or a tsset span of 2^31−2 delta units or more, is r(198) |
| F045 | civreghdfe: `e(df_r)` is exact above 2^31 |
| F047, F050, F056, F057 | civreghdfe and cqreg: memory failures stop with r(920) instead of posting fallback, zero or 2SLS-relabeled results. A diagnostic statistic that cannot be computed is posted missing, not 0 |
| F060 | cbinscatter: `discrete method(binsreg) absorb()` with more than 500 bins per group fails at once (rc 503, with a message) instead of building huge matrices first |
| F012 | cmerge: a master variable named `_merge_temp_filter` is no longer overwritten |

Each finding's code change is described by its ID in the P1 and P2 tables below; the original proposals are in the patches file. The fixes as applied differ in a few places. F002 re-reads the hint into a larger heap buffer rather than passing a length macro, because the native transport harness mocks every macro read. F001 is resolved by F002, with no change to `ctools_rowstage.c`.

**Still open:**
- The 25 P3 findings, including F013 and F014 in `csort_stream.c`.
- Notes from the fix work, not fixed:
  - `cimport_do_scan` keeps a 64 KiB `macro_val` buffer that can cut very wide files' variable-name lists.
  - Streaming csort still ignores SPI return codes during payload writes.
  - civreghdfe CUE posts V = 0 if S is singular at the estimate.
  - cpplmhdfe has no K cap for its int K×K arithmetic (F054).
  - `cexport_impl.c`'s mmap and pwrite paths ignore pool submit results, though not silently.
  - `crangestat.sthlp` should state that each statistic takes one source variable.
- Rebuild the other three platform plugins, run the Windows build and tests, and confirm F004 on Intel or Rosetta.
- Commit.

## Release blockers at a glance

| ID | Sev | Command | What happens | Evidence |
|---|---|---|---|---|
| [F009](#f009) | P0 | `csort ..., stream(#)` | After `regress` / `predict`, variables created after the estimation are left unsorted while the rest of the data is sorted. rc 0, silent row misalignment | Reproduced in Stata. **Fixed** ([status](#p0-fix-status)) |
| [F033](#f033) | P0 | `cexport sas/spss/sasxport5/sasxport8` | Any dataset with a strL after a str# variable: write through NULL or a heap overflow, and Stata crashes | Reproduced: Stata exits with SIGSEGV. **Fixed** ([status](#p0-fix-status)) |
| [F004](#f004) | P0 | `cimport delimited` | Any CSV of 2 MiB or more on an Intel Mac, or under Rosetta: a 524 KB stack array on a 512 KiB worker stack, and Stata crashes. On Apple Silicon about 12 KB of headroom is left | Code and disassembly; not run (no x86 Stata here). **Fixed** ([status](#p0-fix-status)) |
| [F021](#f021) | P1 | `crangestat (stat) new = x y` | Overwrites the existing variable `y` with missing values; rc 0 and no undo | Reproduced in Stata. **Fixed** ([status](#p1p2-fix-status)) |
| [F010](#f010) | P1 | `csort x x, stream(#)` | Heap overflow; keys are sorted, all other columns are not | Reproduced (r(198), keys sorted). **Fixed** ([status](#p1p2-fix-status)) |
| [F011](#f011) | P1 | `csort` streaming | An allocation failure after the keys are written leaves rows scrambled | Code trace. **Fixed** ([status](#p1p2-fix-status)) |
| [F005](#f005) | P1 | `cexport delimited/excel` with any strL | 2,046-byte malloc per str# cell: tens of GB, OOM kill or r(920) | Code trace (prior EXP-9). **Fixed** ([status](#p1p2-fix-status)) |
| [F026](#f026) | P1 | `cimport delimited, stringcols()/numericcols()` | Lists longer than 4 KiB are silently cut: ID columns import as numbers | Code trace. **Fixed** ([status](#p1p2-fix-status)) |
| [F028](#f028) | P1 | `cexport excel` | Chunk-format or submit failures are ignored: rows silently dropped or a corrupt DEFLATE stream, reported as saved | Code trace. **Fixed** ([status](#p1p2-fix-status)) |
| [F045](#f045) | P1 | `civreghdfe [fw=]` | Weight totals above 2^31−1 are narrowed to int (UB): wrong SEs and statistics, rc 0 | Code trace. **Fixed** ([status](#p1p2-fix-status)) |

Two more items need a decision before release:
- **Unchecked SPI stores ([F019, F020](#spi-contract)).** crangestat, crangejoin and csample call the unchecked `SF_vstore` directly. DEVELOPERS.md forbids this, and it is the same root cause as F009. A post-estimation crangestat run did *not* reproduce corruption, but the fix is one line per file.
- **Stale Intel plugin.** `build/ctools_mac_x86.plugin` is dated 26 September 2026, so it lacks every fix made since. F004's crash was measured on that binary. Rebuild every platform's plugin before packaging.

## How the audit was done

- **Hunt.** Eighteen read-only reviewers each read one subsystem line by line, together with its ado wrapper. Five more swept the whole tree for one hazard class each: OpenMP data sharing, integer overflow in sizes and offsets, fixed buffers and string handling, error-path cleanup and ownership, and state retained between plugin calls. Each was told to skip anything recorded as fixed in `OPUS_REPORT_BUGS.md`, `docs/audits/MEMORY_SAFETY_AUDIT_2026-09-28.md` and the SEP22/SEP24/P1/P2–P3 fix notes. (Those reports and notes were removed from the repository on 15 October; the September reports remain in git history at the repository root.) Together they filed 106 reports, which de-duplicate to 69 defects.
- **Verify.** 22 findings got two independent reviews: whether the code really does what the report says, and whether a user can reach it and with what impact. The other 47 were checked in seven batches by subsystem. Every reviewer was told to try to refute the claim.
- **Stata runs.** The P0/P1 claims that could run on this Mac (StataNow MP, arm64) were reproduced against `build/ctools_mac_arm.plugin` as rebuilt on 14 October (sha256 b1569f1d…). The F009 fix was checked with a scratch plugin built from a copy of `src/`. `build/` was not touched.

**Limits:**
- Line numbers are from the evening of 14 October. Other sessions were editing `src/civreghdfe/`, `src/cqreg/`, `src/cimport/`, `src/cbinscatter/` and `src/cpsmatch/` during the audit, so locate code by function name.
- The patches have not been compiled, except the one-line F009 fix.
- Windows and Linux were not run. Windows-only findings (F027, F029) come from code reading.
- The cross-cutting sweeps overlapped the subsystem reviews. Their 36 reports mostly merged into existing findings, but there was no second round to look for gaps.

## Release-blocking findings

<a id="f009"></a>
### F009 [P0] — Streaming csort leaves post-estimation variables unsorted (silent row misalignment)

- **Where:** `src/csort/csort_stream.c`, `csort_stream_apply_permutation`: payload reads (`SF_vdata`, ~495/624) and writebacks (`SF_vstore`, ~556–580 and ~690–700).
- **Mechanism:**
  - The build defines `SD_FASTMODE`, and `csort_stream.c` does not define `SD_SAFEMODE`, so `SF_vdata`/`SF_vstore` are the unchecked `(_stata_)->vdata`/`->store`.
  - DEVELOPERS.md ("Numeric variable indices are positions in the plugin varlist…") warns that raw fast-mode calls can use a different index space. Those of cbsample, cpsmatch, creghdfe and cpplmhdfe were already moved to `SD_SAFEMODE` for this reason; this one was not.
  - When an estimation command has left its hidden sample variable in memory, the stream path permutes the wrong column for variables created after it. Separately, every `SF_vstore`/`SF_vdata` return code is discarded.
- **Reproduction** (exact output from this audit):
  ```stata
  sysuse auto, clear
  quietly regress price mpg
  predict yhat
  gen double id = runiform()
  preserve
  sort id
  tempfile ref
  save `ref'
  restore
  csort id, stream(4)
  cf _all using `ref'        // yhat: 71 mismatches, rc 9 (csort itself returned 0)
  ```
  - With `ereturn clear` before `csort`, the same test passes.
  - The non-streaming path also passes. It uses the shared transport.
- **Reach:** the documented `stream(#)` option. Streaming also starts automatically when there are more than 30 variables and more than 10M observations, unless `nostream` is given, which is exactly the big-data case where users won't notice.
- **Fix (verified):** add `#define SD_SAFEMODE` immediately before `#include "stplugin.h"` in `src/csort/csort_stream.c`. With a scratch plugin built from that one-line change, the reproduction gives `cf` rc 0. A fuller patch is in the patches file: it switches the writebacks to `(_stata_)->safestore` and propagates store errors, so a failed write is no longer reported as success.
- **Regression test:** add the reproduction above (stream(2), stream(4), and the automatic path) to `validation/validate_csort.do`.

<a id="f033"></a>
### F033 [P0] — `cexport sas/spss/sasxport*` crashes Stata when a strL follows a str# variable

- **Where:** `src/io/cio_readstat.c`, `cio_do_export`: strL scan ~632–649, capacity bump ~751–752, realloc ~754–763.
- **Mechanism:**
  - `s = NULL` and `stringcap = 0`.
  - In the first per-variable loop, each str# variable raises `stringcap` to `width+1` without reallocating. The real `realloc` happens only after the loop.
  - A strL processed later in the same loop trusts `stringcap`: it calls `SF_strldata(var, obs, s, (int)stringcap)` and writes `s[length] = 0`. If `s` is still NULL, that writes through NULL. If an earlier strL left a 2-byte block, it is a heap overflow.
  - sasxport8 also overflows by one byte with two strLs, because its logical width is `strlen+1`.
- **Reproduction:**
  ```stata
  sysuse auto, clear
  gen strL note = "x"
  cexport spss using t.sav, replace   // Stata process dies with SIGSEGV (exit 139)
  ```
- **Fix:** keep the real allocation size and the needed size apart, and realloc whenever the need exceeds the capacity:
  ```diff
  -    size_t stringcap = 0;
  +    size_t stringcap = 0, stringneed = 0; /* stringcap is always the real size of s */
  ...
  -        if (isstr && (size_t)width + 1 > stringcap)
  -            stringcap = (size_t)width + 1;
  +        if (isstr && (size_t)width + 1 > stringneed)
  +            stringneed = (size_t)width + 1;
       }
  -    if (stringcap)
  +    if (stringneed > stringcap)
       {
  -        char *tmp = realloc(s, stringcap);
  +        char *tmp = realloc(s, stringneed);
           if (!tmp) { rc = 909; goto done; }
           s = tmp;
  +        stringcap = stringneed;
       }
  ```
  The strL branch already grows `s` from `length + 1` and keeps `stringcap` in step, so it is correct once the str# bump stops inflating `stringcap`.
- **Regression test:** add strL cases to the I/O formats suite (`validation/run_io_formats.py`): one strL after a str# variable, two strLs, and an empty strL, for each of spss, sas, sasxport5 and sasxport8.

<a id="f004"></a>
### F004 [P0] — `cimport delimited` overflows the worker thread stack on Intel macOS (and under Rosetta)

- **Where:**
  - `src/cimport/cimport_impl.c:732`, `cimport_parse_chunk_parallel`, declares `CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS]`: 32,767 × 16 bytes = 524,272 bytes.
  - It runs on threads created with `pthread_create(..., NULL, ...)` in `src/ctools_threads.c:233` (`ctools_persistent_pool_init`), and in cimport's own fallbacks at `cimport_impl.c:1255` and `:1563`.
- **Mechanism:**
  - macOS gives secondary threads a default stack of 524,288 bytes on x86_64 and 536,576 bytes on arm64.
  - The function's frame in the shipped `build/ctools_mac_x86.plugin` is 524,416 bytes (`subq $0x80048,%rsp` plus pushes), larger than the whole stack. Its first call faults in the guard page.
  - On arm64 the frame is 524,432 bytes, which leaves about 12 KB for everything above and below it.
  - Every file of 2 MiB or more is split into one chunk per CPU and parsed on the pool, so this is the default path.
- **Reach:**
  - `cimport delimited using big.csv, clear` with default options on any Intel Mac.
  - On Apple Silicon when Stata runs under Rosetta. `_ctools_load` then picks the x86 plugin.
  - Apple Silicon native is one deeper call chain or compiler change away from the same crash.
- **Fix:**
  - Give pool workers an explicit stack size (8 MiB, the main-thread default), and use the same attribute in cimport's two raw `pthread_create` fallbacks.
  - Move `field_buf` to the heap: allocate it once per task and free it at exit, setting `ctx->error_code` on failure.
  ```c
  #define CTOOLS_POOL_STACK_SIZE ((size_t)8 * 1024 * 1024)
  pthread_attr_t attr, *attrp = NULL;
  if (pthread_attr_init(&attr) == 0) {
      size_t cur = 0;
      if (pthread_attr_getstacksize(&attr, &cur) != 0 || cur < CTOOLS_POOL_STACK_SIZE)
          (void)pthread_attr_setstacksize(&attr, CTOOLS_POOL_STACK_SIZE);
      attrp = &attr;
  }
  /* ... */ pthread_create(&pool->workers[i], attrp, persistent_worker_thread, pool) /* ... */
  if (attrp) pthread_attr_destroy(attrp);   /* also on the failure path */
  ```
- **Regression test:** a multi-chunk import (a file larger than 2 MiB) run with the x86_64 plugin under Rosetta in CI.

<a id="f021"></a>
### F021 [P1] — `crangestat (stat) new = x y` silently overwrites `y`

- **Where:**
  - Root cause: `package/ado/crangestat.ado` ~205–244. `confirm variable` accepts a varlist; `source_var_list` gains two tokens while `nsource` grows by one.
  - The C side, `src/crangestat/crangestat_impl.c` `crangestat_main` ~97–172, reads fixed counts with no cross-checks.
- **Mechanism:**
  - The extra source position shifts every later positional argument by one. y's index is read as the result variable, and the real result index is read as the stat code. The stat code is then read as the source map and goes unvalidated.
  - The compute loop skips the out-of-range stat, and the store loop writes its all-missing column into y.
  - The plugin returns 0, so `_ctools_undo` clears its record.
- **Reproduction** (exact output from this audit):
  ```stata
  clear
  set obs 100
  gen t = _n
  gen x = rnormal()
  gen y = rnormal()
  crangestat (mean) m = x y, interval(t -1 1)   // rc 0; afterwards y is missing in all 100 rows
  ```
- **Fix:**
  - In the ado, require exactly one source variable per statistic (`unab`, then `exit 103` on more than one word).
  - In C, as defence in depth, return 198 when a stat code is outside `[0, STAT_COUNT_ALL)`, when a source map is outside `[0, nsource)`, or when a result index equals a key, source or by index.

<a id="f010"></a>
### F010 [P1] — Repeated key in streaming csort: heap overflow, then keys sorted but other columns not (SORT-8, still open)

- **Where:** `src/csort/csort_stream.c`, `csort_stream_sort` ~927–958.
- **Mechanism:**
  - `nvars_nonkey = nvars - nkeys` counts a repeated key twice. `nonkey_var_indices` and `nonkey_str_widths` are sized from that count, but the fill loop writes one entry per real non-key variable, d entries past the end for d repeats.
  - The `nk != nvars_nonkey` check runs only after the overflow, and after Phase 4 has already written the sorted keys back to Stata.
  - With `nkeys > nvars` the subtraction wraps.
- **Reproduction:** `sysuse auto, clear` then `csort price price, stream(2)` gives r(198). `price` is left sorted and every other column is in its original order.
- **Fix:**
  - De-duplicate the sort keys in `csort_main` right after `parse_sort_vars`: sorting by (x, x) equals sorting by x. Alternatively use `list uniq` in `csort.ado`.
  - In `csort_stream_sort`, count non-keys with an `is_key` bitmap and bound-check the fill before anything is written to Stata.

<a id="f011"></a>
### F011 [P1] — Streaming csort writes the keys before allocating the payload buffers (SORT-9, still open)

- **Where:** `src/csort/csort_stream.c`: Phase 4 key stores ~875–883 run before every Phase 5 allocation (`inv_perm` ~360, batch buffers ~432–451, per-string-variable `calloc(nobs, width+1)` ~222).
- **Impact:** if any of those allocations fails, the command returns r(920) with the keys sorted and some or all payload columns in their original order. Each observation's key then belongs to another row. There is no snapshot or undo. Streaming is the mode for memory-constrained wide data, so this is the realistic way it fails.
  - The non-OpenMP branch (~728–732) even returns success after a failed string permutation. That only matters for local non-OpenMP builds.
- **Fix:** allocate everything Phase 5 needs before Phase 4 writes anything:
  - `inv_perm`, the numeric batch buffers, and a fixed pool of string buffers reused across string variables. If the pool does not fit, shrink it down to one thread.
  - Pass them into `csort_stream_apply_permutation`, so that after the first write only SPI errors remain possible, and report those as data-modified errors.

<a id="f005"></a>
### F005 [P1] — `cexport` with any strL allocates 2,046 bytes for every str# cell (EXP-9, still open)

- **Where:** `src/ctools_text.c:58–62`, `ctools_load_text_data`: `length = islong ? SF_sdatalen(...) : 2045; malloc(length + 1)`. It is reached from `cexport_main` and `cexport_xlsx_main` as soon as any exported variable is a strL.
- **Impact:**
  - Resident memory is about 2 KiB × N × (number of str# variables): 2M rows × 5 str10 columns plus one strL is about 20 GB.
  - On Linux the OOM killer kills Stata. On macOS it swaps heavily. On Windows it fails with r(920).
  - Decoded value-label columns are recast to str2045 and count toward this too.
  - The memory is freed correctly afterwards, so this is over-allocation, not a leak.
- **Fix:** for non-strL variables, read into a 2,046-byte local buffer and allocate `strlen + 1`. Better: copy into a per-column `ctools_string_arena`, so cleanup stays O(1). Keep the `SF_sdatalen` sizing only in the strL branch.

<a id="f026"></a>
### F026 [P1] — `cimport delimited` silently truncates `stringcols()`/`numericcols()` lists at 4 KiB

- **Where:** `src/cimport/cimport_impl.c:1587–1592`, `cimport_parse_col_list`, reads the global `CIMPORT_STRCOLS`/`CIMPORT_NUMCOLS` into `char buf[4096]` with no length check. The ado expands the numlist into that global (`cimport.ado` ~68–74 and ~320–326).
- **Trigger:** `cimport delimited using wide.csv, stringcols(1/1100) clear`.
  - The expanded list is 4,392 bytes, so columns past about 1,040 lose their forced type.
  - A token cut mid-number can force the wrong column.
  - ZIP and ID codes then import as numbers, losing leading zeros or rounding values above 2^53, with rc 0.
- **Fix:**
  - Pass the byte length alongside the list (`CIMPORT_STRCOLS_LEN`), allocate `length+1`, and return 198 unless `strlen(buf) == length`. `cexport_read_local` already does this.
  - Also drop the globals on the scan-error exit (`cimport.ado` ~356 exits before the `macro drop`).

<a id="f028"></a>
### F028 [P1] — `cexport excel` reports success after silently dropping rows or corrupting the DEFLATE stream (EXP-8, still open)

- **Where:** `src/cexport/cexport_xlsx.c`:
  - `xlsx_format_chunk` (~1042–1153): an escape-buffer malloc failure, or a strbuf grow failure that frees the chunk.
  - `xlsx_gen_worksheet_segments` (~1321–1365): the submit, wait and serial results are ignored.
  - `xlsx_comp_worker` (~686–705): an empty final segment skips `TDEFL_SYNC_FLUSH`.
  - `xlsx_compress_parallel` (~781–785): the submit result is ignored. The edit path reaches the same code.
- **Impact:** on a big export (for example 1M rows × 150 or more columns of doubles or long strings) under a memory cap:
  - Up to 50,000 rows per failed chunk vanish while `<dimension>` still claims them, and the command says the file was saved.
  - Or the worksheet's DEFLATE stream is left unflushed, and Excel reports the file as damaged.
- **Fix:**
  - Add a `failed` flag to each chunk; check it and every submit/wait return code, and fail the export (r(920)/r(603)) before publication.
  - In `xlsx_comp_worker`, put the closing flush on the last *non-empty* segment.
  - Run a task serially when its submit fails. A full diff is in the patches file.

<a id="f045"></a>
### F045 [P1] — `civreghdfe [fw=...]` with weights summing above 2^31−1: undefined conversion and wrong inference

- **Where:** `src/civreghdfe/civreghdfe_impl.c:1441` passes `(ST_int)N_eff`. `N_eff` is the double sum of the fweights, and nothing bounds it.
  - Every callee takes `ST_int N_eff`: `ivest_compute_2sls` (`civreghdfe_estimate.c:663`), `ivvce_compute_*` (`civreghdfe_vce.c`) and `civreghdfe_compute_underid_test` (`civreghdfe_tests.c`).
  - `impl.c:1527/1535` repeat the narrowing for `df_r`.
- **Impact:** conversion of an out-of-range double is UB.
  - On arm64 it saturates, so the SEs are inflated by √(Σfw / INT_MAX).
  - On x86-64 it gives INT_MIN, and the dof arithmetic then wraps: negative or garbage dof adjustments.
  - `e(N)` is correct, so the results look plausible, with rc 0.
  - Population-weighted cell data easily exceed 2.1e9.
- **Fix:** carry `N_eff` as `ST_double` through the 2SLS, VCE and test signatures, compute the dof in double, and cap only when storing `df_r`. A full diff is in the patches file. The same `(ST_int)N_eff` pattern should be checked in creghdfe and cpplmhdfe; their audits did not flag it.

<a id="spi-contract"></a>
## Unchecked SPI calls outside the shared transport (fix together with F009)

DEVELOPERS.md requires shared checked I/O or `SD_SAFEMODE` for direct SPI calls. These files still call the unchecked `SF_vstore` directly:

| ID | Site | Pattern | Stata check |
|---|---|---|---|
| F019 | `src/crangestat/crangestat_impl.c` ~1161–1184 (store phase) | `omp parallel for` over statistics, one result variable per thread, unchecked `SF_vstore` | Not reproduced: crangestat after `regress` (result and source variables created after the estimation) matched a clean run |
| F020 | `src/crangejoin/crangejoin_impl.c` ~74–86 (`restore_master`), ~228–251 (write) | 4,096-row tiles in parallel over every column once the output exceeds 10,000 rows, unchecked `SF_vstore` | Not run |
| — | `src/csample/csample_impl.c:406` | Serial unchecked `SF_vstore` | Not run |

The project's own notes disagree about the mechanism. The comment above `IO_VSTORE_FN` blames concurrent unchecked stores. F009 shows a hidden-variable index effect. F009 is the only case confirmed in Stata. Recommendation: define `SD_SAFEMODE` at the top of these three files too (before any include that pulls in `stplugin.h`, e.g. `crangestat_internal.h:30`). That matches the four commands already converted. The cost is the per-cell check; it has not been measured for stores, but checked loads cost 10–23% in the transport study.

## P2 findings (33)

**All 33 are fixed** (15 October; see [P1 and P2 fix status](#p1p2-fix-status)). F001 is resolved by the F002 fix; F019 and F020 were fixed with the P0s.

These are real and reachable, but need unusual conditions: allocation failure, sizes near 2^31, malformed input files, Windows-only paths, or rare option combinations. Several are memory-safety defects in parsers of untrusted files, and opening a crafted file should never corrupt the heap: **F066** (repeated `<sst>` in an .xlsx: heap overflow), **F024** (strL cache writes one pointer past `string_data` when chunk 0 has no rows), **F034** (NULL .xls sheet names: crash) and **F039** (uninitialized `xe_record` read). **F006** is a heap overflow after a 24-byte allocation failure. **F029** is a Windows kernel use-after-free on a failed `WriteFile`. **F012** silently writes `_merge` codes into a user variable. These are the P2 items to take first. Entries marked *uncertain* or *not reproduced* had a plausible mechanism that the reviewers could not settle.

| ID | Location | Defect | Proposed fix (summary; full diff in the patches file) |
|---|---|---|---|
| F001 | `src/ctools_rowstage.c:227` `ctools_rowstage_load` | Row-staged direct string reads can silently overwrite already-staged cells when the width hint is narrower than the variable's real width **(uncertain)** | Minimal fix: in the row-major stager, always use the bounced read, so the SPI never writes outside the slot's own width+1 bytes. Replace lines 226-228 with: for (size_t i = 0; i < count && !io_rc; |
| F002 | `src/ctools_transport.c:113` `ctools_io_ctools_read_strw_from_stata` | __ctools_strw width hints are still read into a fixed 16 KiB stack buffer; a cut entry becomes a wrong (too narrow) width Prior: CORE-2. | Treat possibly cut metadata as no hint. After the SF_macro_use call, add: if (rc != 0 \|\| buf[0] == '\0' \|\| strlen(buf) + 2 >= sizeof(buf)) return NULL; so every variable takes the arena path. |
| F006 | `src/ctools_threads.c:313` `ctools_persistent_pool_submit_batch` | submit_batch partially runs a batch before failing; the pairs-sort serial retry then re-runs scatter tasks with advanced offsets (heap overflow) Prior: SORT-11. | Make batch submission all-or-nothing. Allocate and link all `count` work items first; on any allocation failure free the partial list and return -1 without enqueuing anything. |
| F012 | `package/ado/cmerge.ado:1101` `_cmerge_impl` | cmerge's hard-coded `_merge_temp_filter` plus `merge_var_idx = c(k)` makes the plugin write _merge codes into a user variable and skip reordering it (MRG-10 still present) Prior: MRG-10. | In cmerge.ado, use `tempvar __mfilter` for the scratch variable, and after the append pass `merge_var_idx` from `st_varindex("`generate'`__mfilter'")`, not `c(k)`. |
| F016 | `src/cencode/cencode_impl.c:431` `cencode_main` | cencode label-file path is cut at the first space: labels go to a truncated path and cencode fails (DSTR-8 still open) Prior: DSTR-8. | Pass the path outside the argument string, the same way __ce_map is passed. In cencode.ado, replace `labelfile=...` with `local __ce_labelfile `"`__labelfile'"'` before the plugin call, and clear it afterwards. |
| F017 | `src/cwinsor/cwinsor_impl.c:573` `cwinsor_main` | cwinsor allocates one nobs-sized work buffer per CPU even when only nvars threads use them (SMP-2 still open) Prior: SMP-2. | Compute the worker count from the available work: size_t tasks = ngroups == 1 ? nvars : nvars * ngroups; int nwork = (int)(tasks < (size_t)num_threads ? tasks : (size_t)num_threads); |
| F019 | `src/crangestat/crangestat_impl.c:1174` `crangestat_main` | crangestat stores result columns from parallel OpenMP threads through the unchecked SPI store (SF_vstore is (_stata_)->store under SD_FASTMODE) **(not reproduced)** | Use the checked callback in the parallel loop, as the transport does: before the region add `const ST_IIID store_fn = (_stata_)->safestore;` and replace both SF_vstore calls with `store_fn(result_idx, ...)`. |
| F020 | `src/crangejoin/crangejoin_impl.c:251` `crangejoin_main` | crangejoin write and rollback store numeric columns concurrently through the unchecked SF_vstore **(not reproduced)** | Hoist `const ST_IIID store_fn = (_stata_)->safestore;` and use it at lines 82 and 251 (SF_sstore may stay; the transport already writes strings concurrently through sstore). |
| F022 | `src/cpsmatch/cpsmatch_kernels.c:224` `cpsmatch_find_nearest_neighbors_fast` | cpsmatch nearest-neighbour search can index sorted_controls[-1] (or loop forever) when a distance equals the DBL_MAX sentinel | Choose the side explicitly: `int take_left = left >= 0 && (right >= n_controls \|\| dist_left <= dist_right);`. Use it in the first loop. In the ties loop, use `if (take_left && dist_left <= last_dist + 1e-12) ... |
| F024 | `src/cimport/cimport_impl.c:1455` `cimport_build_cache_worker` | strL cache builder writes one pointer past string_data (and rowrange() reads an uninitialized header pointer) when chunk 0 holds no rows | Start chunking at the header row so that chunk 0 always begins with it. Have cimport_parse_header store `header_start - ctx->file_data` in the context. |
| F025 | `src/cimport/cimport_impl.c:321` `cimport_parse_header` | Scan call chain keeps about 2.2 MB of CTOOLS_MAX_COLUMNS stack arrays on the Stata thread **(uncertain)** Prior: IMP-15. | Heap-allocate these buffers, as in the P0 fix. In cimport_parse_header, allocate one `CImportFieldRef *` buffer of 2*CTOOLS_MAX_COLUMNS entries for fields and next_fields, and a `bool *candidates = calloc(ctx->num_columns, 1)`. |
| F027 | `src/io/cio_iconv_win.h:98` `cio_iconv` | Windows iconv shim reports EINVAL for odd-length UTF-16 without converting the complete prefix, so cimport shifts the whole file by one byte and imports garbage | Match iconv's partial-progress semantics. Convert the even prefix and report EINVAL for the leftover byte: size_t usable = *inleft; int partial = 0; if (c->source == 1200 \|\| c->source == 1201) { partial = (int)(usable % 2); |
| F029 | `src/cexport/cexport_io.c:421` `cexport_io_write_batch` | Windows batched async writes return on a synchronous WriteFile error while earlier writes are still in flight; the caller then frees their OVERLAPPED blocks, events and buffers (kernel use-after-free) | Before any early return in the async section, retire the writes already issued. Capture `DWORD err = GetLastError();` first, then call `CancelIoEx(file->handle, NULL);` and, for each issued entry j, `DWORD w; |
| F030 | `src/cexport/cexport_xlsx_formula.inc:110` `xedit_shared_anchors` | sheet(..., modify) shared-formula anchor relocation is quadratic: every anchor rescans the merged sheet from the start and each relocation copies the whole document | Make two linear passes. First, scan `result` once and build a small table keyed by si: whether a non-empty anchor survives, and the first dependent <f> node with its cell coordinates. |
| F031 | `src/cexport/cexport_format.c:429` `cexport_format_row` | CSV row formatters return the row length as int; rows over INT_MAX bytes are misreported, and rows of 4 GiB or more wrap to a small positive length (silent truncation and overlap) | Before returning, reject rows that do not fit the return type: `if (pos > INT_MAX) return -3;` in all three formatters, documented as 'row too large'. format_chunk_thread's `== -1` growth loop then stops immediately. |
| F032 | `src/cexport/cexport_xlsx_edit.inc:424` `xedit_styles` | xedit_styles ignores a failed re-find of cellXfs and passes a node from a freed-later, different buffer to xedit_replace (cross-object pointer subtraction) | `if (!xedit_find(result, result + strlen(result), "cellXfs", &xfs)) { free(result); return NULL; }` The caller already maps NULL to r(920). |
| F034 | `src/io/cio_excel.c:62` `cio_excel_workbook` | NULL libxls sheet names are dereferenced: crash in `cimport excel .xls, describe`, with `sheet()`, and when exporting into an existing .xls | Treat a NULL name as "" everywhere it is read. Add a helper to cio_internal.h, `static inline const char *cio_xls_sheet_name(xlsWorkBook *w, unsigned j) { const char *s = w->sheets.sheet[j].name; return s ? s : ""; |
| F036 | `src/io/cio_dbf.c:417` `cio_export_dbf` | cexport dbase converts out-of-range %td values with (int)x (UB), then overflows cio_days_civil Prior: XIO-12. | `double d = floor(x); if (d < -2000000000.0 \|\| d > 2000000000.0) { rc = 108; goto done; } cio_days_civil((int)d, &year, &month, &day);`. floor() also fixes the off-by-one for negative fractional dates. |
| F039 | `src/io/cexport_xls_edit.inc:152` `xe_globals` | xe_globals/xe_sheet read an uninitialized xe_record when the first record cannot be parsed | Zero-initialize: `xe_record r = {0}, nr = {0}, ...` in xe_globals and `xe_record r = {0};` in xe_sheet, so a failed first read returns 610. |
| F041 | `src/ctools_matrix.c:81` `matmul_xb_blocks` | Row-block count overflows for N within 2047 of INT_MAX, so ctools_matvec_add/ctools_matmul_xb skip every row and the results come out silently wrong | Compute the block count without overflowing: `const ST_int num_blocks = N / CTOOLS_ROW_BLOCK + (N % CTOOLS_ROW_BLOCK != 0);` (or do the addition in size_t/long long). |
| F042 | `src/ctools_matrix.c:61` `ctools_matmul_ata` | K x K sizes and indices are computed in int in the shared OLS/matrix kernels; with 46,341 or more regressors this overflows and corrupts the heap Prior: REG-7. | Index and size the K x K arrays in size_t, for example `C[(size_t)j * K + i]`, `memset(C, 0, (size_t)K1 * K2 * sizeof(ST_double))`, `calloc((size_t)K * K, ...)`, and `ctools_safe_malloc3((size_t)K, (size_t)K, sizeof(ST_double))` i … |
| F043 | `src/cbsample/cbsample_impl.c:442` `cbsample_main` | cbsample with strata() and cluster(): finding each stratum's clusters takes time proportional to strata x clusters (near-hang on large data) Prior: SMP-3. | Keep a running cursor across strata, since both strata and clusters are in sorted order and clusters nest in strata. Declare `size_t cl_cursor = 0;` before the `for (st ...)` loop. |
| F046 | `src/civreghdfe/civreghdfe_impl.c:1395` `do_iv_regression` | tsset time index truncated to int: spans of 2^31 or more delta units wrap tidx, which mis-pairs kernel lags, drops DK periods, and can make pair_next read past om->order | In impl.c, after the tmin/tmax loop and before computing tidx, reject spans that do not fit: if (!((tmax - tmin) / tdelta < (ST_double)INT_MAX - 1.0)) { SF_error("civreghdfe: the tsset time span is too long for kernel estimators ( … |
| F047 | `src/civreghdfe/civreghdfe_estimate.c:82` `ivest_compute_gmm2s` | Unchecked small allocations in GMM2S and the k-class path (NULL dereference); a failed CUE or GMM Hessian allocation silently falls back | In ivest_compute_gmm2s: - Delete WZtX. - After allocating XZW, XZWZX, XZWZy, XZWZX_L and beta_temp, add `if (!XZW \|\| !XZWZX \|\| !XZWZy \|\| !XZWZX_L \|\| !beta_temp) { free(...all...); return 920; }` before first use. |
| F050 | `src/civreghdfe/civreghdfe_tests.c:601` `civreghdfe_compute_underid_test` | Diagnostic test functions are still void: allocation or solve failures post 0 statistics (p = 1) or fall back to the iid Anderson LM labelled as Kleibergen-Paap, with no warning Prior: IV-22. | Make each civreghdfe_compute_* return ST_retcode (0, 920, or -1 for a numerical failure), and set the outputs to SV_missval rather than 0 before any early return. |
| F056 | `src/cqreg/cqreg_linalg.c:185` `cqreg_compute_quantile` | cqreg_compute_quantile swallows an allocation failure and returns 0, so e(q_v), e(sum_rdev) and e(r2_p) are silently wrong with rc 0 | Report the failure instead of returning a substitute value. Change the function to `ST_int cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau, ST_double *out)`, returning -1 when malloc fails. |
| F057 | `src/cqreg/cqreg_vce.c:175` `cqreg_sandwich_product` | cqreg_sandwich_product zeroes V on allocation failure and returns void, so cluster and robust-kernel VCEs post zero standard errors with rc 0 Prior: QREG-8. | Make it `static ST_int cqreg_sandwich_product(...)` returning -1 on overflow or allocation failure and 0 on success. In each caller, write `if (cqreg_sandwich_product(V, H, XtX, K) != 0) goto cleanup;` with rc still -1; |
| F060 | `src/cbinscatter/cbinscatter_binsreg.c:277` `adjust_bins_binsreg_hdfe` | discrete + method(binsreg) + absorb(): dense N x U bin-indicator matrix and (U+K)^2 normal equations are built before the >500-bin capacity check, giving multi-GB allocations and effectively unbounded … Prior: BIN-8 (related; failure noted, not the pre-failure blow-up); BIN-12 (N x nquantiles scalability). | Fail fast in C before any adjustment. In do_compute_bins, read the matrix capacity once on the main thread (ST_int bin_cap = SF_row("__cbinscatter_bins");) and keep ST_int rows_used = 0. |
| F061 | `src/cbinscatter/cbinscatter_resid.c:302` `ols_residualize` | Block-count arithmetic overflows int for N within one block of 2^31: validity flags left uninitialized, or residualization silently skipped | Compute the block count without overflow at both sites, e.g. `const ST_int num_blocks = N / block + (N % block != 0);` (impl.c:351) and `const ST_int num_blocks = N / BLOCK + (N % BLOCK != 0);` (resid.c:302), or use `(ST_int)(((in … |
| F062 | `src/cbinscatter/cbinscatter_resid.c:130` `compute_xtx_weighted` | memset length K*K computed in int: overflows for >= 46,341 control columns and can memset a huge length | Widen before multiplying: `memset(XtX, 0, (size_t)K * (size_t)K * sizeof(ST_double));`. Or delete the memset, since lines 135-155 assign every element of XtX. |
| F066 | `src/cimport/cimport_xlsx.c:557` `shared_strings_callback` | XLSX import: a repeated <sst> element reallocates the shared-string table without resetting its count (heap overflow) | Preallocate only once, and only while empty: `if (unique_count && !ctx->shared_strings && ctx->num_shared_strings == 0)`. Set capacity only after both allocations succeed. |
| F067 | `package/ado/cimport.ado:618` `cimport` | An error after `cimport scan` leaves the whole parsed CSV context (file mmap plus per-row parse arenas) resident until the next ctools call Prior: CORE-6. | In cimport.ado, call `capture plugin call ctools_plugin, "cimport clear"` before `exit 103` and inside the `if `rc'` block before `exit`. |
| F068 | `package/ado/crangejoin.ado:135` `_crangejoin_impl` | crangejoin keeps its full using-data index when _crangejoin_using fails after its plugin call: the clear guard global is set too late | Set `global CTOOLS_crangejoin_plugin "`__ctools_plugin'"` before line 133. A rollback with state.ready==1 frees the state without writing, so this is safe. |

## P3 findings (25: 7 confirmed, 18 latent)

Latent means a real defect that no current input reaches, for example because every caller already guards it. They are worth fixing when the file is next touched.

| ID | Location | Defect | Proposed fix |
|---|---|---|---|
| F003 | `src/ctools_data_store.c:433` `ctools_io_data_store_impl` | Division by zero in the store scheduler when storing a stata_data with zero variables (latent) | In ctools_io_data_store_impl, after resolving nvars, add: if (nvars == 0) return STATA_OK; (nothing to write). Also make ctools_io_chunks_per_column defensive: if (nvars == 0 \|\| min_chunk_rows == 0) return 0; |
| F007 | `src/ctools_types.c:152` `ctools_int64_to_str` | ctools_int64_to_str negates INT64_MIN (signed-overflow UB) (latent) | Negate in unsigned arithmetic: `return 1 + ctools_uint64_to_str(0 - (uint64_t)val, buf + 1);` |
| F008 | `src/ctools_sort_pairs.c:132` `ctools_sort_order_pairs` | Pairs radix sort orders only the low 32 bits of the size_t order_key, and its tie order follows input order rather than orig_idx (latent; SORT-12 unfixed) (latent) Prior: SORT-12. | Compute the maximum key once, then run `passes = max_key ? (bit_width(max_key) + 7) / 8 : 1` passes, capped at 8, instead of a fixed 4. |
| F013 | `src/csort/csort_stream.c:222` `stream_permute_string_var_flat` | Flat string buffer in streaming csort trusts the ado width hint: SF_sdata writes into the slot with no guard tail and no length check (latent) | Read each value into a local `char tmp[2046]`, reject `strlen(tmp) > max_width` with an error, then memcpy into the slot. |
| F014 | `src/csort/csort_stream.c:495` `csort_stream_apply_permutation` | Streaming csort's 16-wide read loop reads observation obs_map[i] + bb instead of obs_map[i + bb] (latent) | Use `(ST_int)obs_map[i + bb]` at lines 495 and 624. |
| F015 | `src/cmerge/cmerge_impl.c:1095` `cmerge_load_using` | cmerge plugin argument parsers do not cross-check counts: negative counts, uninitialized/stale index arrays, and an unchecked keepusing count that indexes the using cache (latent) Prior: MRG-9. | In load_using, reject nkeys < 0 and n_keepusing < 0, and after the token loop return 198 unless arg_idx > nkeys and keepusing_idx == n_keepusing. In execute, memset the three static arrays at entry; |
| F018 | `src/cdestring/cdestring_impl.c:467` `cdestring_main` | cdestring total converted/failed counters are int and overflow (signed UB) when nvars x nobs > INT_MAX | Declare the totals as `double total_converted = 0, total_failed = 0;` (or uint64_t), accumulate with `total_converted += (double)var_converted[v];`, and save them without the int cast. |
| F035 | `src/io/cio_context.c:150` `cio_days_civil` | cio_days_civil signed overflow for Excel date-formatted serials near INT_MAX | Tighten the guard to the domain cio_days_civil supports: `days < (double)INT_MIN + 146097 \|\| days > (double)INT_MAX - 715815` gives ".". Alternatively, compute `z` and `era` in `long long` inside cio_days_civil. |
| F037 | `src/io/cio_shp.c:408` `cio_shp_export` | cexport shp converts an arbitrary _ID double to uint32_t (UB for missing, negative or huge IDs) | Validate before converting: `if (SF_is_missing(id) \|\| id < 0 \|\| id > (double)UINT32_MAX \|\| id != floor(id)) { free(blob); rc = 198; goto done; }`. |
| F038 | `src/io/cexport_xls.c:312` `cexport_xls_main` | cexport excel .xls: cell() column letters accumulate without a bound (signed int overflow); row parsed with atoi | Bail out in the loop as soon as `col0 > 256` (rc 198). Parse the row with `strtol` with errno/range checks (1..65536). Optionally restrict the ado regex to `^[A-Z]{1,3}[0-9]{1,7}$`. |
| F040 | `src/io/cio.c:8` `cio_main` | cio_main leaves `kind`/`index` uninitialized when the argument string has only an operation (latent) | `char op[24] = "", kind[24] = "";`, and return 198 for ops that require a kind when `!*kind`. |
| F044 | `src/creghdfe/creghdfe_regress.c:996` `do_full_regression` | creghdfe inverts the Cholesky factor in place, aliasing L and inv against ctools_invert_from_cholesky's documented contract (latent) | Factor xtx_keep itself, which is dead after this point, and invert into the separate buffer. Drop the memcpy: `if (ctools_cholesky(xtx_keep, K_keep) != 0) { ... |
| F048 | `src/civreghdfe/civreghdfe_impl.c:520` `do_iv_regression` | Signed overflow in the missing-value block count when more than INT_MAX-16383 observations are selected | Replace line 520 with: const ST_int num_blocks = N_total / block + (N_total % block != 0); |
| F049 | `src/civreghdfe/civreghdfe_estimate.c:287` `ivest_compute_cue` | CUE copies beta onto itself with memcpy (overlapping memcpy is UB) | Use `if (beta != initial_beta) memcpy(beta, initial_beta, K_total * sizeof(ST_double));` (or memmove). |
| F051 | `src/civreghdfe/civreghdfe_omega.c:134` `civ_omega_with_kernel` | civ_omega_with_kernel can turn a two-way-with-kernel builder into a kernel-free one whose intersection clusters were never built (cl12 = NULL, ncl12 = 0) (latent) | In civ_omega_with_kernel, keep the original kernel for TWOWAY when the copy would need intersection ids that were not built: `if (c.kind == CIV_OMEGA_TWOWAY && om->kernel > 0 && kernel <= 0) kernel = CIVREGHDFE_KERNEL_BARTLETT;` ( … |
| F052 | `src/civreghdfe/civreghdfe_tests.c:436` `civ_kp_rk` | int products D*D and K1*K3 in the Kleibergen-Paap kernel overflow for K_endog*L > 46340 (latent) | Use size_t loop counters (`for (size_t i = 0; i < (size_t)D * D; i++)`) and `memset(C, 0, (size_t)K1 * K3 * sizeof(ST_double))` in ctools_matmul_ab. Allocate B and T only on the K > 1 path. |
| F053 | `src/cpplmhdfe/cpplmhdfe_irls.c:1574` `do_ppml_regression` | A failed write of the d() variable is reported as success because ppml_save_effects overwrites output_rc (latent) | Keep the first error. Before line 1574 add `if (output_rc) goto ppml_cleanup;`, or write `if (!output_rc) output_rc = ppml_save_effects(...);`. |
| F054 | `src/cpplmhdfe/cpplmhdfe_irls.c:548` `do_ppml_regression` | K_x*K_x is computed in int for matrix buffer sizes and row-major indices (signed overflow for K_x >= 46341) (latent) Prior: REG-7 (related: column offsets k*N were fixed; these K*K products were not). | Widen before multiplying: `malloc((size_t)K_x * (size_t)K_x * sizeof(ST_double))` (likewise for K_keep, K_with_cons, and `(size_t)K * K` in ppml_qrsolve). Index with `(size_t)a * K + b`. |
| F055 | `src/cpplmhdfe/cpplmhdfe_irls.c:889` `do_ppml_regression` | ctools_invert_from_cholesky is called with L aliasing inv, against its documented contract (latent) | Factor into a separate buffer and invert into the destination. The existing workspaces already allow this: run Cholesky on xtx_k (or xtx) itself, which is not reused after the solve, then call `ctools_invert_from_cholesky(xtx_k, K … |
| F058 | `src/cqreg/cqreg_regress.c:1012` `cqreg_full_regression` | Uninitialized u/v arrays are read after a non-converged main solve | Move the `if (!state->converged) { ...; rc = 430; goto cleanup; }` check (lines 1029-1033) directly after line 1008, before the cqreg_ipm_get_objective and cqreg_ipm_get_residuals calls. |
| F059 | `src/cqreg/cqreg_fn.c:1568` `cqreg_fn_solve_core` | C side never checks N > K; with N < K the crossover reads past the N-element idx array (latent) | In cqreg_full_regression, after `ST_int N = (ST_int)filtered.data.nobs;`, add `if (N <= K) { ctools_filtered_data_free(&filtered); ctools_error("cqreg", "not enough observations"); return 2001; }`, matching the ado. |
| F063 | `src/cbinscatter/cbinscatter_impl.c:316` `load_data` | By-group count derived as max-min+1 from a raw (ST_int) cast of the by value: UB/overflow and dropped groups if the dense-id contract from the ado is ever broken (latent) Prior: BIN-7. | Enforce the contract in C. In load_data map any by value that is not an integer in [1, INT_MAX] to -1 (`by_groups[i] = (SF_is_missing(val) \|\| val < 1 \|\| val > INT_MAX \|\| val != floor(val)) ? -1 : (ST_int)val;`). |
| F064 | `src/cbinscatter/cbinscatter_bins.c:179` `compute_bin_statistics` | compute_bin_statistics places the SE double accumulators after the int counts array, so they are misaligned whenever nquantiles is odd (latent) | Put the int array last: sum_y = block; sum_x = sum_y + nq; sum_w = sum_x + nq; if (compute_se) { sum_y2 = sum_w + nq; sum_x2 = sum_y2 + nq; counts = (ST_int *)(sum_x2 + nq); } else counts = (ST_int *)(sum_w + nq);. |
| F065 | `src/ctools_data_store.c:641` `ctools_store_filtered` | Public store helpers call SF_error from OpenMP worker threads (cwinsor/cdestring parallel column stores), reachable with Break (latent) Prior: CORE-7. | Split out a quiet variant (e.g. ctools_store_filtered_noreport) that returns the rc without calling SF_error. Use it inside the parallel loops, and call ctools_io_store_status once on the calling thread after the region. |
| F069 | `src/cimport/cimport_impl.c:2007` `cimport_blob` | cimport blob dereferences a NULL string cache for string columns that were not chosen as strL (latent) | Add `\|\| !ctx->columns[col-1].use_strl` to the argument check. After the cache build, also return 198 if `!ctx->col_cache \|\| !ctx->col_cache[col-1].string_data \|\| !ctx->col_cache[col-1].string_data[row-1]`. |

## Refuted

| ID | Location | Report |
|---|---|---|
| F023 | `src/cpsmatch/cpsmatch_impl.c:535` `cpsmatch_main` | cpsmatch per-thread weight/workspace arrays are sized by ctools_get_max_threads() but used in parallel regions with no num_threads clause |

F023 depends on the OpenMP thread-count ICV exceeding `ctools_get_max_threads()`, and the dispatcher resets the ICV on every call.

## Suggested order of work

1. **F009 plus the SPI contract.** Add `#define SD_SAFEMODE` to `csort_stream.c`, `crangestat`, `crangejoin` and `csample`, and add the post-estimation regression to the csort and crangestat validations.
2. **F033.** Fix the capacity bookkeeping in `cio_do_export`.
3. **F004.** Set an explicit pool thread stack size and move `field_buf` to the heap.
4. **Rebuild all four platform plugins.** The x86 macOS plugin is from 26 September 2026.
5. **The seven P1 items:** F021, F010, F011, F005, F026, F028, F045.
6. **The memory-safety P2 items:** F066, F024, F006, F034, F039, F029, then F012.
