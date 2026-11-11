# cimport / cexport bug audit and fix plan (2026-11-08)

Scope: delimited, Excel and wrapper code for `cimport` and `cexport`, reviewed on the
working tree. That tree includes the uncommitted cimport strL/perf work and the
uncommitted `col_prefix` fix in `cexport_xlsx.c`.

Status column:
- **run** means reproduced in StataNow 19 against native `import`/`export`. Delimited
  and Excel cases ran on the committed 1.2.4 Windows plugin with HEAD's ado files;
  the code involved is unchanged in the working tree.
- **read** means confirmed by tracing the code but not run.
- **plausible** means the code path is clear, but the native behaviour or the
  trigger was not checked.

## Fix status (2026-11-08)

Batches A and B are fixed in the working tree, which is not committed. Every
fix was checked with a zig-built Windows plugin against StataNow 19.

**Validation suite** (`validation/stata.py` `run_groups`)
- **cimport:** 395 passed. The 2 failures are the known `encoding(latin1)`
  baseline cases.
- **s2_strl:** 113 passed, none failed.
- **cexport:** 219 passed. The 9 failures fail the same way on HEAD; they are
  byte comparisons without `crlf`, so Windows line endings break them. With
  `crlf` all of them match native.
- **New regression tests fail on HEAD:**
  - the digits test (B4);
  - the three label tests (B5/B6), which give r(111);
  - the wide shared-string sheet (B1): native `import excel` cannot read
    HEAD's file, r(603).

**A1** was checked separately on a 3.5 MB multi-chunk strL file, imported
through all four strL routes (memory, file, tempdta, strldta) with strict,
loose, `threads(1)` and `threads(3)`. Every result matches native. Builds used:
- tiny limits: `-DCIMPORT_STREAM_THRESHOLD=65536 -DCIMPORT_CACHE_BUDGET=400000
  -DCIMPORT_BLOB_BLOCK_ROWS=50 -DCIMPORT_BLOB_BLOCK_BYTES=4096`;
- the same with `-DCIMPORT_CACHE_BUDGET=1`, where nothing is cached, so
  chunk 0's header skip goes through the split;
- the normal build.

The `CIMPORT_BLOB_BLOCK_*` limits are now `#ifndef` build-time settings for
that purpose.

**How each item was fixed**

| # | Fix |
|---|---|
| A1 | Compact blocks or suffixes over the limits are cut by `cimport_strl_split`, which runs a cursor pre-scan identical to the transfer's re-parse. Normal-sized cached blocks keep the field-ref fast path. The file-route region is capped at 2,147,483,620 bytes, where values fall back to verbatim bytes. |
| A2 | `strofreal(most, "%20.0f")`. |
| A3 | `threads_code` is passed to fulldta, strldta, blobmem, blobs and blobfill. |
| A4 | The home frame is recorded before phase 2, and the restore runs `frame change` back to it first. |
| A5 | `error_code` is reset at the start of the load and of `fulldta`. |
| A6 | `cimport blob` returns 198 for a column without cached strL text. |
| B1 | The uncommitted terminator is kept, with a 40-column regression test. |
| B2 | Only `[h]`, `[hh]`, `[m]`, `[mm]`, `[s]`, `[ss]` count as elapsed. |
| B3 | A trailing field after a final delimiter at EOF is kept (zero-length). |
| B4 | The `floor(log10())` exponent is corrected with exact powers of ten. |
| B5 | An undefined label falls back to unlabeled values, for delimited and excel. |
| B6 | Delimited floats use `%14.0g` of the value rounded to 8 significant digits. Excel uses `%9.0g` text, which is what native writes: the native workbook stores these cells as shared-string text. |
| B7 | The timing table is skipped when no load ran, and the timing scalars are dropped at the start and end of every import. |
| B8 | `DECIMALSeparator`, `GROUPSeparator`, `SHeet`, with the help updated. |

**Batch C** is also fixed and not committed. The suite gives the same results
(cimport 395/2, cexport 222/9, s2_strl 113/0), and its three new batch C tests
pass. Checks on Python-built workbooks compared the new build with HEAD and
with native.

| # | Fix | Evidence |
|---|---|---|
| C1 | The file is pre-sized from the sampled average row; buffers keep the maximum. A failed pre-size is no longer fatal for pwrite (positioned writes extend the file, close truncates), and any mmap pre-size failure falls back to pwrite. | The suite's skewed-row parallel export (120,000 rows, three of them 1,900 bytes) is byte-identical to native. The disk-full path itself was not exercised. |
| C2 | Sheets are matched by relationship type. Without `sheet()` the first worksheet is used; a named chart, dialog or macro sheet gives r(601). | On a chart-first workbook, HEAD's `sheetmodify` failed with a bare r(920) and `sheet("Chart1", replace)` wrote worksheet XML into the chart part. Now `sheetmodify` edits `Data` with the chart part byte-identical, and both named modes give r(601). Native fails on this workbook with r(9901). |
| C3 | A run of the 11 ctools cell formats and the date-time numFmt is reused; formatCode is compared ignoring `\` escapes. | Suite: after 5 sheet additions and a modify, styles.xml still holds 11 cellXfs and 1 numFmt, and the data import intact. |
| C4 | The end pointers are computed once. | Code only; not timed. |
| C5 | In a Strict workbook (main namespace `purl.oclc.org`), the parts written are rewritten to the Strict namespaces. | Modify plus an added sheet: every part parses with only Strict (or OPC) namespaces, and cimport reads both sheets back. |
| C6 | `rel_id[256]`. | Long-id workbook whose tab order differs from its part names: HEAD read sheet A's data for B and the reverse; now it matches native. |
| C7 | Revised after the run against native. Native numbers a row without `r=` by position and ignores its cells' row in `r=`. The audit's proposal (honour the cell's row) was the opposite. Inside such a row a cell now keeps the row's number and takes only its column from `r=`. | Without `cellrange` native gives N=4, HEAD 7, now 4 and identical to native. With `cellrange()` native and cimport both give r(198). |
| C8 | `cimport_parse_row_fast` returns `CIMPORT_FIELD_TOO_LONG` for fields over `CIMPORT_FIELD_MAX_LENGTH` (settable at build time), reported as r(109) with a message. | A test build with a 4000-byte limit gives r(109) and the message for a 5,000-byte field and for an unclosed quote, under loose and strict binding. A 3,000-byte field still loads. |

**Batch D** is fixed in the ado files only, and not committed. Every case
below was run against native StataNow 19. URL cases used a local HTTP server
on 127.0.0.1.

Suite results after batch D:

| Group | Passed | Failed | Notes |
|---|---|---|---|
| cimport | 409 | 2 | the 14 new D tests pass; 2 latin1 baseline failures |
| cexport | 222 | 9 | baseline |
| s2_strl | 113 | 0 | |
| contracts | 18 | 0 | |
| data_io | all | 0 | |
| option_contracts | 99 | 15 | the known CRLF "CSV formatting" baseline |

The package contract and option inventory checks pass. `io_excel_*` and the
other fixture cases cannot run on this machine (symlink privilege).

| # | Fix | Evidence |
|---|---|---|
| D1 | New helper `build/_ctools_home.ado` expands a leading `~` (HOME, or USERPROFILE on Windows), as `_cpq_path` already did for Parquet. It is called for delimited and Excel import, the `_cio_import` formats, and delimited, Excel and `_cio_export` export. It is listed in ctools.pkg and in the developer README. | `~/file.csv` imports as native does; cexport delimited and excel write into the home directory. |
| D2 | A URL is detected before `.csv` is added, and `r(filename)` is the URL. | `…/plain` and `…/plain?format=csv` import as native does; `r(filename)` is the URL. |
| D3 | `macval()` is used through the dispatch, `cimport_excel` and `_cio_import`. One-line `if … local …` forms with user text became blocks: Stata expands the command of a one-line `if` a second time even with `macval()`, as the trace showed. | `sheet("Q1\$Rev")`, `dol$lar.csv`, `dol$lar` (no extension) and `bo$ok` (xlsx, no extension) all match native. |
| D4 | The label and blob globals are dropped at the start, on scan failure, on every early exit and on rollback. Timing scalars were already handled in B7. | After a failed import (r(103)) no `_cimport_*` global is left. |
| D5 | Native's messages for `maxquotedrows()`, `numericcols()`/`stringcols()` (now r(198), not 125/121), `rowrange()`/`colrange()`, `delimiters()`, excel `describe` (now also exclusive of `allstring(fmt)`), `allstring` with `allstring()` (`opts_exclusive`), `allstring(%bad)`, `cellrange(A1:)`, and a varlist with `firstrow` or `case()`. `_cio_import` uses `error 4` instead of a bare `exit 4`. A negative `maxquotedrows()` now means no limit: native accepts it and treats it like 0. | 15 option cases give the same return code as native, and the messages checked match its wording. |
| D6 | Excel URLs are copied to a temporary file that keeps the URL's extension, which is erased on every path. `.xlsm` is accepted as OOXML. | URL `book.xlsx` and `book.xlsm` import as native does. |

The help (`cimport.sthlp`) now covers Excel URLs, `~`, `.xlsm` and negative
`maxquotedrows()`. New suite tests in `validation/stata/cimport.do` cover:
- `$` in file names and in a sheet name;
- `.xlsm`;
- globals after a failed import;
- 7 option return codes;
- `maxquotedrows(-1)`;
- the `~` helper.

URLs and writes to the home directory are left to the separate checks above.

**Batch E** is fixed and not committed. Suite results:

| Group | Passed | Failed | Notes |
|---|---|---|---|
| cimport | 414 | 2 | latin1 baseline |
| cexport | 226 | 9 | baseline |
| s2_strl | 113 | 0 | |
| contracts | 18 | 0 | |
| data_io | all | 0 | |

Native was probed before each change; several audit proposals changed as a
result.

| # | Fix | Evidence |
|---|---|---|
| E1 | Notes record each row's start offset; the first 10 by file position are kept, and they are shown as file line numbers. Native names the line: blank lines and `varnames(#)` rows count. The parallel parse's 50-bytes-per-row guess is gone. | 4 files with blank lines, `varnames(2)` and two notes match native. The suite compares the row numbers. |
| E2 | The audit was wrong here: auto-detection already uses native's rule; the header row did not. Native's `varnames(#)` is the row holding file line # (a row spanning line # is the header), also under strict. `cimport_skip_rows` now skips rows that end before that line; under loose and nobind binding it behaves as before. | Strict `varnames(3)` and `varnames(2)` now take native's header row. **Left open:** on a file whose header and data use different delimiters (`,` vs `;`), auto-detection still picks `;` where native picks `,`. |
| E3 | Auto-detection does not sample past `rowrange()`'s last line after skipping blank lines. | `rowrange(:2)` with a blank line 2 matches native. |
| E4 | `<v>` text has CDATA and entities decoded and surrounding whitespace trimmed (not for `t="str"`). | Native reads `<v> 1 </v>` as TRUE, CDATA `5` as 5, a CDATA shared index and `<v> 7 </v>`; cimport now matches, where it was r(610). |
| E5 | A failed output allocation while inflating, including miniz's `MZ_ZIP_ALLOC_FAILED`, sets the memory refusal, so callers give r(909) "larger than the memory available". | Code only. |
| E6 | One used-range helper for import and describe: formatted BLANK/MULBLANK cells count at every corner. A cell libxls never read is BLANK with format 0. | `.xls` files made by putexcel with a blank at D10 and at A1 match native's import and describe ranges. Before: describe A1:B3, and k=2 instead of 4. |
| E7 | The audit proposal was wrong here: native caps `missing()` at **16 bytes** ("invalid string in option missing() / The string is too long.", r(198)). cexport now does the same, which also covers Excel's 32,767-character limit. | Bisected against native: 16 bytes accepted, 17 refused. |
| E8 | The numeric test for `missing()` uses the plugin's locale-independent decimal parser. | 12 values against native: hex `0x1A` is now text, as native writes it. `inf`, `nan`, `NA` and `1,5` stay text; ` 2`, `1e3`, `+8` and `.5` stay numbers. |
| E9 | `calcPr` is inserted after the last of sheets, functionGroups, externalReferences and definedNames. | Workbook with definedNames and extLst: order is sheets, definedNames, calcPr, extLst. |
| E10 | Proposal changed: native refuses (r(9901)) a sheet name that equals another only ignoring the case of non-ASCII letters, in both add and modify, rather than matching it. cexport refuses too, with a message, using a simple case fold (Latin-1, Latin Extended-A, Greek, Cyrillic). | "DONNÉES" against "Données" now gives r(9901) and no duplicate sheet; before, a second sheet was added. |
| E11 | Every edit failure sets a message: missing or malformed part, memory, or output. A sheet part without sheetData gives r(610) with a message instead of r(920). The relationship id is also found under another prefix (`ns1:id`). Data past XFD or row 1,048,576 gives native's "invalid range in option cell()". | The ns1-prefixed workbook now edits (was r(610)). The suite checks XFD against native. |
| E12 | Native keeps nothing of the old workbook with `replace keepcellfmt`. The help is corrected and the dead `preserved_styles` code removed. | Native output: one cellXfs, no numFmts. |
| E13 | Single-threaded export: a Stata read failure or failed buffer allocation is r(920) with its own message, as in the parallel path. r(693) is kept for write failures. | Code only. |
| E14 | The single-threaded buffer is one worst-case row reserve plus 64 KiB rather than two reserves. Peak memory for a long strL falls by about 2× its length; ordinary rows are unchanged (64 KiB–4 MiB buffer). | Code only; the suite's long-string and strL export tests pass. |
| E15-E17 | `cexport.md` quoting (no quoting for line breaks), `cimport.sthlp` `asdouble` (noninteger columns only) and excel `verbose` (accepted only). | Docs. |

**Found while fixing** (not fixed; see section F)
- the `emptylines()` mismatch;
- Excel CR handling.

## A. Fix before committing the in-flight cimport strL work

| # | Finding | Where | Scenario | Proposed fix | Status |
|---|---|---|---|---|---|
| A1 | Compact chunks with strL columns no longer keep row refs, so strL transfer blocks have no row or byte cap (a regression vs HEAD). | `cimport_impl.c:1684` (HEAD also tested `compact_rows`); `cimport_strl_plan` 2887-2907 | A 3 GB CSV with about 600 KB per strL row, or `threads(1)`: one 4096-row block spans more than 2 GB. `blobmem` refuses with 459 (3189). `blobs` sends a region of up to 4 GB, and Mata's `fread` fails with "file changed while it was imported" (ado 1168-1170). `tempdta` mallocs the whole span for each group. | Split compact blocks and uncached suffixes into sub-blocks of at most `CIMPORT_BLOB_BLOCK_ROWS` and 64 MB, cutting by cached field offsets or a cursor pre-scan. Cap the file-route region at `INT32_MAX`. Simplest alternative: restore HEAD's `compact_rows` retain. | read |
| A2 | `set obs` count formatted with `%9.0g`. | `cimport.ado:1062` | A block of 10^9 or more entries gives `set obs 1.000e+09`. Reachable only through A1. | `strofreal(most, "%20.0f")` | read |
| A3 | `threads()` is not passed to the strL and .dta routes. | `cimport.ado:494, 520, 546, 554, 1069` | `threads(2)` still uses every core in `blobmem`/`blobs`/`fulldta`/`strldta`/`blobfill`. | Prefix those plugin calls with `threads_code`. | read |
| A4 | A Break or out-of-memory inside `_cimport_read_blobmem` leaves the temp frame current, so rollback restores into the wrong frame. | `cimport.ado:1077-1083`, restore at 601 | Break during `st_data` in the temp frame: the backup is copied into the temp frame, and the user's frame keeps the half-imported data. | Save `c(frame)` before phase 2, and `capture frame change` back before restoring. | plausible |
| A5 | A failed malloc in `cimport_context_number` sets the context-wide `error_code=920`, which persists into the fallback load. | `cimport_impl.c:124`; read at 2351/2435/2526 | `fulldta` runs out of memory, the ado falls back to the SPI load, and that load fails with "write failed" on the first text cell in a numeric column. | Reset `error_code` at the start of `cimport_do_load`/`fulldta`, or return the failure instead of storing it. | plausible |
| A6 | The legacy `cimport blob` subcommand dereferences NULL for `str#` columns. | `cimport_impl.c:2767` | `plugin call ctools_plugin, "cimport blob 1 1 0"` on a `str#` column crashes Stata. The ado no longer uses it. | Return 198 unless `use_strl`, or delete the subcommand and its phase entry. | read |

## B. Confirmed wrong results (small, independent fixes)

| # | Finding | Where | Scenario | Proposed fix | Status |
|---|---|---|---|---|---|
| B1 | cexport excel: the cell prefix is not NUL-terminated but is printed with `%s` (in HEAD and the released 1.2.4 plugins). | `cexport_xlsx.c:1332-1347`, `xlsx_shared_cell` 1129-1137 | Any string sent to the shared-string table (one containing `_x`, CR or other control characters, or the `missing()` text) can emit garbage into the sheet XML. It did not show up in a 2-column test because the heap bytes happened to be zero. | The working tree already adds the terminator. Commit it (or use `calloc`), rebuild the plugins, and add a wide-sheet regression with `_x`/CRLF strings. | read (HEAD lacks the line) |
| B2 | cimport excel: `[Magenta]` in a number format is classified as elapsed time. | `src/io/cio_excel_format.h:67-70` | `[Magenta]0.00` on 12.5: native gives double 12.5, cimport gives `%tc` 12jan1900 12:00:00. XLS is affected too. | Treat a bracket as elapsed only when it holds a single repeated `h`/`m`/`s` (`[h]`, `[mm]`, `[ss]`); skip other brackets. | **run** |
| B3 | cimport delimited: a trailing delimiter at EOF with no final newline loses its empty field. | `cimport_parse.c:379-393` | `a,b\n1,2\n,` gives 2 obs native, 1 in cimport (the row is taken as blank). `1,2\n3,4,` with `varnames(nonames)` gives 3 vars native, 2 in cimport. | After the loop, if `field_start == end && end > start && end[-1] == delim && field_count < max_fields`, append a zero-length field. Check that the strL cursor and chunk parsers share this path. | **run** |
| B4 | cexport delimited: `floor(log10())` rounds up just below a power of ten, so a significant digit is dropped. | `cexport_format.c:179-183` | Native writes `99.99999999999999` and `999.9999999999999`; cexport writes `100` and `1000`. | After `floor(log10)`, correct with an exact comparison against `cexport_pow10`: decrement when `magnitude < 10^e`, increment when `magnitude >= 10^(e+1)`. | **run** |
| B5 | cexport delimited: a value label that is attached but not defined aborts the export. | `cexport.ado:131-141` | `label values y nolab` with `nolab` undefined: native writes the numbers, cexport fails with r(111) "value label nolab not found". Check the excel path (`cexport.ado` around 552) too. | Decode only when the label exists (`mata: st_vlexists()`); otherwise treat the variable as unlabeled. | **run** |
| B6 | cexport delimited: unlabeled values of a labeled float variable are formatted with `%18.0g`. | `cexport.ado:143, 153` | A float 0.1 with no label is written `.1` by native and `.1000000014901161` by cexport. | Format the fallback with the writer's own rule (8 significant digits for float, 16 for double). Either use a float-aware format in the ado or leave those cells to the C formatter. | **run** |
| B7 | cimport delimited: `verbose` on a header-only file fails after the import. | `cimport.ado:613-637` (625) | A file containing only `a,b` with `verbose` gives r(111) "_cimport_time_mmap not found". On the `preserve` route that undoes the import. A stale scalar from an earlier run shows the wrong timings. | Skip the timing table when no load ran (`capture confirm scalar`), and drop the timing scalars on every exit. | **run** |
| B8 | Option abbreviations are stricter than native. | `cimport.ado:54` (`DECIMALSEParator`, `GROUPSEParator`); excel `SHEET` at 703 | `decimals(,)` works natively but gives "option decimals() not allowed" r(198) in cimport. Native uses `DECIMALSeparator`, `GROUPSeparator` and `SHeet`. | Use the native capitalisation. | **run** (decimals) / read (sheet) |

## C. Robustness and file-integrity bugs

| # | Finding | Where | Scenario | Proposed fix | Status |
|---|---|---|---|---|---|
| C1 | cexport delimited pre-sizes the file from the largest sampled row (or 2 × the longest strL) and treats a pre-size failure as fatal. | `cexport_impl.c:138-142, 573, 590-600`; `cexport_io.c:218-251` | 20M rows with a mostly-empty str2045 column pre-size to about 52 GB for a 2.4 GB file. Windows `SetEndOfFile` allocates it, and the export fails with r(693) when the disk is short. `mmap` has no fallback. | Estimate from the average row size, keeping the max only for buffers. Make pre-size failure non-fatal for pwrite, and return `UNRESERVED` on Windows so mmap falls back. | read |
| C2 | cexport excel editing treats chart sheets as worksheets. | `cexport_xlsx_edit.inc:572-593, 601-616, 765-773` | A first tab that is a chart, or `sheet("Chart1", replace)`: worksheet XML is written into `xl/chartsheets/…` and Excel drops the part. With `modify` it gives r(920) with no message. | Require a relationship `Type` ending in `/worksheet`. Skip non-worksheets when choosing the default, and return 601 for a named chart, macro or dialog sheet. | read |
| C3 | Every cexport excel edit appends 11 `xf` entries and a new custom `numFmt`. | `cexport_xlsx_edit.inc:428-521`, called at 738 | 300 sheet additions leave 300 duplicate numFmts and 3,300 xfs, past Excel's custom-format limit. | Reuse an existing identical `numFmt`/`xf` before appending. | read |
| C4 | Counting existing shared strings is quadratic. | `cexport_xlsx_edit.inc:689` (`old+strlen(old)` inside the loop) | Editing a workbook with about 200k shared strings takes minutes to hours. | Hoist the end pointer. The rels loop near 603 has the same pattern. | read |
| C5 | cexport excel hard-codes the transitional namespace when editing. | `cexport_xlsx_edit.inc:265-269, 365-368, 414-419, 450-453, 478-482, 512-516, 640-647`; `cexport_xlsx_formula.inc:300-302` | `modify` on a Strict OOXML workbook moves the rows into the wrong namespace, and Excel repairs or drops the sheet data. | Take the namespace from each part's root element, or refuse Strict packages with a clear error. | read |
| C6 | cimport xlsx: relationship ids are truncated to 31 bytes. | `cimport_xlsx.h:57`, `cimport_xlsx.c:453` | Long GUID-style `r:id`s never match, so the reader silently falls back to `sheetN.xml` (wrong sheet) or gives 610. | Make `rel_id` 256 bytes. | read |
| C7 | cimport xlsx: with rows that have no `r=` but cells that do, `cellrange` filters on the wrong row. | `cimport_xlsx.c:1424-1439` with 1237/1243 | `<row><c r="A5">` counted as row 2 under `cellrange(A3:B10)`, so A5 is silently dropped. | Recompute the skip from the cell's resolved row, and move the row counter forward to it. | read |
| C8 | cimport delimited: a field of 2^31 bytes or more overflows into the quote-flag bit. | `cimport_parse.c:328, 360, 389` | An unterminated strict quote early in a 3 GB file is misreported as r(920), or the last field is silently garbled. | Detect spans above `0x7FFFFFFF` and return a clear r(109)-style error. | read |

## D. Wrapper parity and state hygiene

| # | Finding | Where | Proposed fix | Status |
|---|---|---|---|---|
| D1 | A leading `~` is not expanded for delimited and Excel imports or for either cexport format (macOS/Linux). The plugin opens `~/…` literally, so you get r(601)/r(603). | `cimport.ado:107-122`, excel 736-769; `cexport.ado:73-79, 819` | Move `_cpq_path`'s `~`/absolute-path logic (`_cimport_parquet.ado`) into a shared helper and use it everywhere. | plausible (no Mac/Linux run) |
| D2 | `.csv` is appended before the URL test, so `…?format=csv` becomes `…?format=csv.csv`. `r(filename)` returns the temp path. | `cimport.ado:108-120` | Test `pathisurl` first and return the URL. | read |
| D3 | Dispatch to excel, parquet and the other formats uses `` `0' `` instead of `` `macval(0)' ``, so `$name` and backticks in `sheet()` or filenames are expanded twice. | `cimport.ado:25, 31, 38`; `cimport_excel` 749-769 | Use `macval` throughout. | read |
| D4 | Globals and scalars are left behind: `_cimport_label_*` (header text) on several error exits, `_cimport_blob*` after a Mata error, and `_cimport_time_*` after every non-verbose run. | `cimport.ado:393-432, 548-550, 597-605, 663-665` | Drop them on every exit. | read |
| D5 | Invalid options exit 198 with no message: `maxquotedrows`, `_cimport_range`, `_cimport_delimiters`, several excel checks. `numericcols(0)` returns 125 instead of 198. `describe` exclusivity ignores `allstring(fmt)`. | `cimport.ado:74-76, 708-744, 778-784, 1221-1232` | Print native-style messages, and map `numlist` errors to 198. | read |
| D6 | Native `import excel` copies URLs to a temp file; cimport rejects them. `.xlsm` is also rejected. | `cimport.ado:736-762` | Copy URLs to a temp file that keeps the extension, and accept `.xlsm` (or route by signature). | plausible |

## E. Low severity

- cimport delimited:
  - Unmatched-quote warning row numbers are estimated in the parallel parse, as 50 bytes per row (`cimport_impl.c:1032-1039`). Record file offsets and convert them to line numbers after the join.
  - With `varnames(#)` under strict binding and a multi-line quoted field above the header, delimiter auto-detect starts at a physical line while the header parse uses logical rows (`cimport_impl.c:1376`).
  - The auto-detect `last` bound is not enforced after blank lines (1133-1157).
- cimport xlsx/xls:
  - `<v> 1 </v>` booleans read as FALSE, and `<v>` holding CDATA reads as 0 (`cimport_xlsx.c:1377`).
  - Out of memory while inflating a part is reported as r(610), not r(909) (`xlsx_zip.c:268-272`).
  - The XLS used range counts BLANK cells for the bottom-right corner only, and `describe` uses yet another rule (`cio_xls.c:211-231`, `cio_excel.c:51-63`). This needs a native .xls fixture first.
- cexport excel:
  - `missing()` text is not checked against the 32,767-character cell limit (`cexport_xlsx.c:431-437`).
  - The numeric test for `missing()` uses `strtod`, which depends on the locale and accepts hex (479-484).
  - `calcPr` is appended after `extLst`/`pivotCaches` (`cexport_xlsx_formula.inc:297-302`).
  - The sheet-name case match covers ASCII only (`cexport_xlsx_edit.inc:578`).
  - Edit failures give r(610) or r(920) with an empty message (553ff, 739-777), and more than 1,048,576 rows gives a silent r(198).
  - `keepcellfmt` with `replace` does not keep styles as the help says (`preserved_styles` is dead code).
- cexport delimited:
  - Allocation failures and Stata read failures are reported as r(693) "write error" (`cexport_impl.c:283-296, 356-360`).
  - A very long strL needs about 6× its size in buffers (86-113, 258-263).
- Docs:
  - `docs/commands/cexport.md` says default quoting covers embedded newlines. Native and cexport both write `a<LF>b` unquoted (run), so fix the doc.
  - `cimport.sthlp:223` misdescribes `asdouble`: integer columns stay integer.
  - `cimport.sthlp:107` says excel `verbose` shows progress; it does not.

## F. Found while fixing A and B (run against native, not fixed)

| # | Finding | Evidence | Proposed fix |
|---|---|---|---|
| F1 | (Corrected in "F1 and F2 status": the probe ran under `version 16`.) cimport `emptylines()` differs from native in two ways. Native accepts only `skip` and `include`; cimport also accepts `fill`, which native rejects with r(198). Native `include` in StataNow 19 added no observations for blank lines, while cimport adds one per blank line plus one for the final line end. | On `id,b / 1,x / (blank) / 2,"two<LF>lines" / (blank) / 3,z`: strict gives N=3 native vs 6 cimport; loose gives 5 vs 8. LF and CRLF behave alike, and HEAD behaves the same. No validation case covers `emptylines`. | Find out what native `include` actually does (it may only differ for some line shapes), make cimport match it, reject `fill`, and add native-comparison cases. |
| F2 | cexport excel writes CR as `_x000D_`, which keeps it. Native `export excel` writes a raw CR, which XML readers (native `import excel` included) normalise to LF, and it turns CR LF into one character. | `"a" + char(13) + "b"` reads back as `a<LF>b` from the native workbook and `a<CR>b` from cexport's. `a<CR><LF>b` reads back as 3 characters vs 4. | Decide which behaviour to follow. Native parity means writing raw CR, at the cost of losing the CR. Either way, document it. |

### F1 and F2 status (2026-11-09)

Both are fixed and not committed.

**F1.** The original finding came from a probe that ran under `version 16`.
Native `import delimited` picks its Java reader by the caller's version:
`plugins16` below 17, the current reader from 17. Only the old reader reads
`emptylines(include)` as `skip`; from version 17 `include` keeps every empty
line as an observation, which is what cimport did.

The fix follows the caller's version:
- `include` maps to the plugin's fill mode only when `_caller() >= 17`, and to
  skip otherwise.
- `fill`, which cimport accepted but native rejects, now gives native's message
  and r(198).

Checked on 11 file shapes, under both bindings and all three option values:
66 combinations match native under version 19 and 66 under version 16. Suite
tests cover 6 shapes × 2 bindings at the suite's version, plus a version 16
check and the `fill` refusal.

**Probes under `version 16`.** Only `import delimited` switches readers by
version. Every earlier item verified with such a probe (B3, E1, E2,
`maxquotedrows(-1)`) also has a suite test that passes at version 19, and E3
was rerun at 19.

**F2.** Decided for native parity, on both sides:
- **Writer:** cexport writes a CR as it is, as native does (tab and LF already
  were), instead of `_x000D_`. XML readers, Excel and native `import excel`
  read a CR or CR LF there as LF.
- **Reader:** the reader lacked XML's end-of-line handling, so a native-written
  CR stayed CR. `xlsx_xml_decode_text` now reads a literal CR or CR LF as LF
  before entities and `_xHHHH_` are decoded, so `&#13;` and `_x000D_` still
  give CR.

Native and cimport now read identical text from native- and cexport-written
workbooks holding CR, CR LF, LF, a leading CR and trailing CRs. Both helps
describe this.

## Checked and ruled out

- **cexport delimited:** "Windows plugin writes 3-digit exponents (`e+016`)" is false. The committed plugin writes `1.00000000000e-06` and `1.00000000000e+20`, byte-identical to native.
- **cexport excel:** "`missing()` ignored for value-labelled variables" is not a bug. Native also leaves those cells blank.
- **cimport excel:** `t="d"` ISO cells read as their leading number. Native does the same (2024 → 16jul1905).
- **cimport excel:** an empty custom-height row 1 becomes the `firstrow` header. Native does the same.
- **Design note, not a bug:** cexport writes LF by default, while native `export delimited` on Windows writes CRLF. This is documented in `docs/commands/cexport.md`, with `crlf` as the opt-in.

## Suggested order

1. **A, B, C, D, E, F1 and F2:** done (see Fix status). Commit A together with
   the strL work it belongs to.
2. **Residuals:**
   - E2's delimiter auto-detection on files whose header and data use
     different delimiters;
   - cimport's unmatched-quote notes have hard line breaks where native's
     wrap.
