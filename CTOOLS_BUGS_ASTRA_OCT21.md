# ctools pre-release bug report — October 21, 2026

Reviewed revision: `48f72e198bb6e44f39c72edac080e19481af083b` (version 1.0.2).

This report contains **15 findings**: runtime defects, platform compatibility defects, and release validation gaps. The original audit below records the pre-fix behavior. Implementation and regression-test status are tracked in the resolution section. The original release blockers were silent XLSX corruption, partial file publication, Windows path handling, ineffective `snapshot` protection, and a bypassable publication gate.

**Severity:** P1 means a release blocker because of data integrity, loss of expected protection, a major supported-platform failure, or missing release validation. P2 means an important correctness, compatibility, or tooling defect to fix before release if feasible. “Reproduced” means observed locally; “source-confirmed” identifies a definite code path whose target-platform execution was not available.

## Resolution status — implementation follow-up

The P1 fixes were implemented first. Fixes for all 15 findings are committed to `main`. Windows native runtime checks and both Mac CI builds have passed; **licensed cross-platform Stata verification remains pending**. Original reproductions below describe the audited revision, not the corrected behavior. Original location line numbers refer to that revision.

| ID | Status | Implemented fix and evidence |
|---|---|---|
| 1 | Resolved locally | Both libdeflate paths check CRC and exact inflated size. Native sanitizer tests cover both extraction paths, stored/deflated entries, and payload/CRC/length corruption. Fresh-plugin Stata rejects the corrupt workbook and retains the original data. |
| 2 | Resolved locally | A complete structural XML pass precedes the range/parallel scanner. All premature EOF positions, mismatched tags, unclosed attributes/comments, prefixed tags, empty sheets, and a 4,000-row parallel/range scan are tested. Stata rejects truncated XML without replacing data. |
| 3 | Resolved locally | SHP/SHX and XPORT companion outputs use a shared commit group with preflight, backups, and rollback. Native tests inject a second-publication error and verify both original files; licensed regressions verify an existing SHX leaves no SHP and a blocked SAS companion leaves neither XPORT nor XPF output. Destination aliases are rejected before publication. |
| 4 | Windows native checks passed; licensed validation pending | A shared strict UTF-8-to-UTF-16 path layer handles long paths. CSV, binary codecs, XLS, miniz, ReadStat, staging, publication, and cleanup use wide Windows APIs. Native Windows CI passes Unicode and long-path reads, writes, replacement and cleanup in the C locale. Full codec round trips pass on Mac; licensed Windows codec validation is still required. |
| 5 | Resolved locally | Explicit `snapshot` preserves before any destructive clear or scan. Licensed tests cover malformed files in all six binary readers, a post-load fault, retained data/labels/characteristics, and successful snapshot import. |
| 6 | Resolved locally | Trim complete UTF-8 Unicode White_Space sequences at both boundaries. Licensed comparisons with native `destring` pass for ten ASCII/Unicode whitespace cases. |
| 7 | Resolved locally | Predictions use one sample marker and a temporary result, renaming only after success. Licensed regressions cover `if` plus `in`, all twelve prediction options, and an error leaving no requested output variable. |
| 8 | Resolved locally; target mounts pending | Prefer exclusive rename, retain the hard-link fallback, and use an exclusive private copy when both are unsupported. Fault-injected sanitizer tests verify successful publication and no overwrite. Real exFAT/Linux mounts still need execution. |
| 9 | Resolved locally; Linux ACL execution pending | Capture ordinary creation mode without changing process umask; preserve replacement mode, group and POSIX ACLs. Native tests verify new mode 0640, replacement mode 0664, and retention of a Mac extended ACL. |
| 10 | Resolved locally | Short random sibling names no longer extend the destination basename. A 250-character destination component passes the native test. |
| 11 | Windows 64-bit offset check passed; large SHP export pending | SHP offsets and length arithmetic use explicit 64-bit types and seek/tell wrappers (`_ftelli64`/`_fseeki64` on Windows). The Windows CI harness successfully seeks and reports a position beyond 2 GiB without allocating a huge file. An actual multi-gigabyte SHP export still needs licensed validation. |
| 12 | Resolved in workflow; runner setup/execution pending | Candidate packaging is separate from publication. Both publishers require successful licensed validation across a four-platform matrix. Launch supports an explicit executable; CI checks candidate revision/platform. Tests exercise successful, skipped, failed and cancelled gate decisions. |
| 13 | Resolved locally | Both workflow path filters include `package/**`; regression checks cover canonical ado/help/metadata paths. |
| 14 | Resolved locally | Staging and the documented `--integrity` check validate payload hashes, exact inventory and platform scope. Tests reject corrupt plugins/ado files, altered manifests, missing/wrong checksums and scope changes. |
| 15 | Resolved; Windows MSYS2 check passed | A Python cleanup command replaces shell-specific `del`; it removes generic and platform plugins and propagates errors. The real Windows MSYS2 CI job successfully removes stale generic/platform plugins and objects. Local tests also verify cleanup errors propagate. |

Regression entry points: `validation/test_oct21_native.py`, `validation/test_oct21_release.py`, and `validation/test_oct21_stata.py` (all registered). `validation/test_windows_files.py` is registered for execution in the real Windows MSYS2 build job. The licensed tests accept `CTOOLS_STATA_EXECUTABLE` and `CTOOLS_TEST_PLUGIN_DIR` to select a fresh candidate without replacing shipped binaries.

Publication requires repository variable `CTOOLS_STATA_RUNNERS` mapping `mac_arm`, `mac_x86`, `windows`, and `linux` to trusted runner label arrays. Each runner must set `CTOOLS_STATA_EXECUTABLE`, provide Python 3 and the required licensed reference commands. Without that configuration, publication stays blocked.

**Remaining operational limits:** Multi-file rollback handles ordinary errors, not power loss or unsynchronized concurrent writers. If rollback itself fails, the original backup is retained and its path is reported. The last-resort no-hard-link copy preserves no-overwrite semantics but exposes a private destination name during copying; exclusive rename remains the normal path. XML validation bounds nesting at 256 elements. Licensed Windows/Linux, minimum-OS, and Stata 14.1 runs remain release exit criteria.

### Validation of the fixes

- Built a fresh Mac ARM plugin from the working tree in `/tmp/ctools-oct21-fixed-build`, revision label `oct21-fixes`; its macOS dependency/baseline check passed. Shipped plugin binaries were not replaced.
- New licensed regression suite passed against that plugin, including corrupt XLSX rollback, binary snapshots, SHP/XPORT companion failures, Unicode paths through CSV/DBF/XLS/XLSX/SAS/SPSS/XPORT, Unicode whitespace, and all twelve supported PPML prediction options.
- Existing licensed component suites passed **740 checks, 0 failures, 0 skips**: import 360, export 213, destring 167. Checksum-pinned official fixtures were prepared for these runs.
- New native sanitizer and release-contract regressions passed. Existing binary-codec, CSV-parity, XLS/XLSX-edit, I/O CRC/XML/CSV, September 24 atomic-output, release-gate, and build-contract suites passed.
- Package/release synchronization, registry validation, Python compilation, workflow YAML parsing, and whitespace checks passed.
- These are targeted validations, not a completed four-platform release gate. Licensed Windows/Linux execution, a multi-gigabyte Windows SHP export, actual exFAT mounts, the minimum supported OS/Stata versions, and OpenMP execution remain unverified locally. The stricter publication gate intentionally blocks release until configured trusted runners pass.

### GitHub Actions follow-up

- Source fixes were pushed as `e1d637ce368ae7a9f84939677a5b23d5099f0908`; [the first Actions run](https://github.com/mdroste/stata-ctools/actions/runs/39673901895) passed Mac ARM, Mac Intel and Windows plugin builds and dependency contracts. Windows native file API and MSYS2 cleanup regressions passed.
- The initial Linux run failed only in `oct21_native`: its fault-injection wrapper used fixed arguments and conflicted with the included production source's variadic `syscall` declaration. Corrected the wrapper to match that interface. All other registered native suites passed in the first run; the local native sanitizer regression passes after the correction. Linux verification follows in CI.
- Repository inspection confirmed no `CTOOLS_STATA_RUNNERS` variable and zero registered self-hosted runners. The licensed matrix and binary publication remain blocked until trusted licensed runners are configured.

## Scope and validation

The review covered the platform loader and configuration, Windows/POSIX file APIs, allocation and threading infrastructure, import/export engines, selected data-management and prediction paths, package generation, and CI publication.

Local execution used macOS 27.0 on Apple Silicon, Apple Clang 21.0.0, Python 3.14, and StataNow/MP 19.5. A fresh plugin was compiled from the reviewed source into `/tmp/ctools-oct21-build`, with the documented pthread-only fallback because a static OpenMP runtime was unavailable. Targeted Stata probes loaded this fresh plugin and reported its revision as `48f72e198bb6e44f39c72edac080e19481af083b`.

| Check | Result |
|---|---|
| Release identity synchronization | Passed |
| Package source synchronization | Passed |
| Existing package completeness checker | Passed; its limitations are finding 14 |
| Fresh macOS ARM build | Passed, without OpenMP |
| Dependency checks on shipped ARM and Intel Mac plugins | Both passed |
| Registered native suites | 35 attempted; 24 completed successfully; 11 failed during OpenMP setup/build because `brew` or the expected static `libomp.a` was unavailable |
| Targeted Stata reproducers | Confirmed findings 1, 2, 3, 5, 6, 7, and 10; also confirmed new-file mode in finding 9 |
| Focused native/Python/shell probes | Confirmed findings 8, 9, 14, and the shell behavior in 15 |

The native suite's overall result was **failure**, not a green gate. Several of the 11 affected suites passed their serial cases before reaching the missing OpenMP dependency. These setup failures are not evidence of numerical failures. No full licensed Stata correctness gate, Windows execution, Linux execution, Intel execution, minimum-OS run, or Stata 14.1 run was completed. No claim of complete statistical or cross-platform certification is made.

The shipped ARM plugin reported revision `a782a528285564eea06bc2f2de540ff4716aaf8d`, the parent source commit of the reviewed binary-publication commit. That relationship is expected and is **not** reported as a stale-binary bug.

## Findings at a glance

| ID | Priority | Finding | Affected platforms | Evidence |
|---|---|---|---|---|
| 1 | P1 | XLSX worksheet inflation skips CRC checks and accepts changed values | All | Stata reproduction |
| 2 | P1 | Truncated worksheet XML is accepted as a successful partial import | All | Stata reproduction |
| 3 | P1 | SHP/SHX export publishes one file before discovering failure on the other | All | Stata reproduction |
| 4 | P1 | UTF-8 filenames reach ANSI Windows file APIs | Windows with a legacy process/CRT code page | Source-confirmed; Windows run needed |
| 5 | P1 | `snapshot` cannot restore the original data for several binary imports | All | Stata reproduction |
| 6 | P2 | `cdestring` mishandles Unicode whitespace accepted by `destring` | All | Native Stata comparison |
| 7 | P2 | PPML prediction combines two `if` qualifiers and leaves an unfinished output | All | Stata reproduction |
| 8 | P2 | Export without `replace` requires hard-link support | macOS/Linux filesystems without hard links | Native fault injection + platform documentation |
| 9 | P2 | Atomic export replaces file permissions with mode 0600 | macOS/Linux | Native and Stata reproduction |
| 10 | P2 | Temporary-file suffix makes valid destination names too long | All, subject to filesystem/path limits | Native Stata comparison on Mac |
| 11 | P2 | SHP export uses 32-bit file positions on 64-bit Windows | Windows | Source-confirmed; large-file Windows run needed |
| 12 | P1 | Release publication succeeds without licensed Stata validation | All distribution targets | Workflow inspection |
| 13 | P2 | CI does not trigger for changes confined to canonical package sources | All distribution targets | Workflow inspection |
| 14 | P2 | Package validation ignores recorded checksums | All distribution targets | Python reproduction |
| 15 | P2 | Windows `make clean` uses `cmd.exe` commands in the documented MSYS2 shell | Windows/MSYS2 | Windows recipe reproduced under POSIX `sh` |

## 1. [P1] XLSX fast extraction ignores the worksheet checksum

**Locations:** [src/formats/xlsx_zip.c:153](/Users/mike/Documents/GitHub/stata-ctools/src/formats/xlsx_zip.c:153), [src/formats/xlsx_zip.c:238](/Users/mike/Documents/GitHub/stata-ctools/src/formats/xlsx_zip.c:238), [src/cimport/cimport_xlsx.c:1533](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1533).

**Defect and impact.** Both libdeflate extraction paths return successfully decompressed bytes without comparing their CRC-32 with the ZIP directory's `m_crc32`. Successfully inflating DEFLATE data does not establish that its contents match the ZIP checksum. The worksheet reader uses this fast path, so corruption can become ordinary numeric data without an error. Standard miniz extraction verifies the checksum, making behavior dependent on which extraction path is used.

**Reproduction and observed result.** A small XLSX contained numeric cells `123` and `456`. The worksheet was ZIP-DEFLATED with compression level zero, making a single-bit mutation of the first digit possible without breaking the DEFLATE stream. Changing ASCII `1` to `9` by XOR with `8`, while leaving the ZIP checksum unchanged, produced:

```text
Python zipfile: Bad CRC-32 for file 'xl/worksheets/sheet1.xml'
cimport excel: RC=0, N=2, values 923 and 456
```

Changing just the recorded checksum also yielded a successful ctools import. Appendix A generates the fixtures.

**Proposed fix.** In both fast paths, require the decompressed length to equal the advertised uncompressed length and compute/compare CRC-32 before returning the buffer. A checksum or length mismatch must return an import error, release the buffer, and preserve the existing dataset under the normal Excel rollback contract. Reuse the vendored checksum implementation; do not allow a fallback path to bypass a failed integrity check.

**Regression test.** Test valid deflated and stored worksheets, a one-bit payload mutation that still inflates, wrong CRC metadata, and wrong uncompressed length. Exercise direct and fallback paths; every corrupt case must fail before replacing data.

## 2. [P1] The worksheet scanner treats malformed/truncated XML as success

**Locations:** [src/cimport/cimport_xlsx.c:1249](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1249), [src/cimport/cimport_xlsx.c:1425](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1425), [src/cimport/cimport_xlsx.c:1452](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1452), [src/cimport/cimport_xlsx.c:1470](/Users/mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx.c:1470).

**Defect and impact.** Failure to find a closing markup delimiter returns from `xlsx_scan_rows` without setting a parse error. Missing `</sheetData>` falls back to the end of the buffer, and `xlsx_scan_buffer` also returns success for incomplete markup. The final `cells_stored < cells_seen` check cannot catch data never recognized as a complete cell. A damaged or incompletely generated XML part can therefore import missing or incomplete data successfully.

**Reproduction and observed result.** A ZIP with valid CRCs contained a worksheet ending inside its second value:

```xml
...<row r="1"><c r="A1"><v>123</v></c></row>
<row r="2"><c r="A2"><v>45
```

Python's ZIP integrity check passed. `cimport excel ..., clear` returned `RC=0, N=2`, with `A[1]=123` and `A[2]=.`. The XML was not a complete worksheet. This defect is independent of finding 1: checksum validation alone will not fix it.

**Proposed fix.** Add a distinct malformed-XML status, propagate it through cell/row scans and parallel worker summaries, and validate the required enclosing structure and unfinished cells. Distinguish a normal worker chunk boundary or requested cell-range boundary from premature end of the XML part. Reject malformed input before committing the dataset.

**Regression test.** Truncate inside a tag, attribute, numeric value, inline string, row, and worksheet close. Include valid empty sheets, self-closing cells, range-restricted imports, and parallel row boundaries to avoid rejecting valid partial scan regions.

## 3. [P1] SHP/SHX publication is not transactional across the pair

**Location:** [src/io/cio_shp.c:539](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:539).

**Defect and impact.** `cio_export_shp` commits the `.shp` file first and then commits `.shx`. If the second commit fails, the first publication remains. Cleanup removes temporary files only. The command reports failure while leaving a new geometry file paired with an old index, or while having already replaced an existing geometry file.

**Reproduction and observed result.** Import the point fixture from Appendix A. Create `partial.shx` containing a sentinel, ensure `partial.shp` does not exist, and run:

```stata
cexport shp using "/tmp/ctools-oct21-probes/partial.shp", shx
```

The command returned **602**, but `partial.shp` existed afterward and `partial.shx` still contained the old sentinel. This requires no race or injected disk failure.

**Proposed fix.** Prepare and validate all output destinations before publication, then use a coordinated commit with backups and rollback for failures partway through publication. Preflight alone is insufficient because later permission, sharing, or filesystem errors remain possible. Do not remove existing destinations until replacement outputs are complete, and preserve a recoverable copy until the whole operation succeeds. Document any unavoidable crash-consistency limitation of a multi-file commit.

**Regression test.** Cover an existing SHX without `replace`, a blocked SHX destination with `replace`, and injected failure during the second commit. On ordinary error, both destinations must retain their pre-command contents, or both must remain absent. Audit XPORT companion-file publication under the same requirement; it also performs multiple publications, though its commit order differs.

## 4. [P1] Windows file handling does not consistently convert UTF-8 paths

**Locations:** [src/formats/ctools_output.c:28](/Users/mike/Documents/GitHub/stata-ctools/src/formats/ctools_output.c:28), [src/formats/ctools_output.c:52](/Users/mike/Documents/GitHub/stata-ctools/src/formats/ctools_output.c:52), [src/cexport/cexport_io.c:188](/Users/mike/Documents/GitHub/stata-ctools/src/cexport/cexport_io.c:188), [src/cexport/cexport_impl.c:326](/Users/mike/Documents/GitHub/stata-ctools/src/cexport/cexport_impl.c:326), [src/io/cio_dbf.c:27](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_dbf.c:27), [src/io/cio_shp.c:69](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:69), [src/io/cio_xport.c:10](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_xport.c:10), [src/vendor/libxls/src/ole.c:625](/Users/mike/Documents/GitHub/stata-ctools/src/vendor/libxls/src/ole.c:625).

**Defect and impact.** Stata passes UTF-8 paths; the CSV reader explicitly acknowledges this and uses `MultiByteToWideChar` plus `CreateFileW`. The shared export staging layer instead uses narrow `_open`, `MoveFileExA`, and `remove`; other export and import paths use `CreateFileA` or narrow `fopen`. Under a legacy Windows code page, a UTF-8 path can fail to resolve or be interpreted as a different filename. This affects non-ASCII usernames/directories as well as filenames.

The inconsistency is especially hazardous for XLSX export: vendored miniz uses wide-character file opening on the Windows build, while staging and commit use narrow APIs. Those components can operate on different filesystem paths. This consequence is inferred from the source; it was not executed on Windows.

**Trigger to validate on Windows.** With a non-UTF-8 process code page, export and import through paths such as `C:/ctools-test/é.csv` and `C:/ctools-test/測試/data.xlsx`, exercising creation, replacement, workbook modification, and cleanup. Check actual directory entries and file bytes, not merely return codes.

**Proposed fix.** Introduce shared UTF-8 path wrappers that convert once and use `CreateFileW`, `_wopen`/`_wfopen`, `MoveFileExW`, and `_wremove` consistently. Route direct first-party I/O and vendor adapters through those wrappers. Support long paths consistently and return explicit conversion errors. Avoid changing the host process's global code page as a workaround.

**Regression test.** Run the Windows plugin under both legacy and UTF-8 process code pages, with accented, CJK, supplementary-plane, space-containing, and long paths. Verify that failed exports leave no incorrectly named temporary or final files.

**Platform reference:** Windows process UTF-8 behavior requires an explicit configuration; it is not guaranteed for every supported Windows host. [Microsoft: Use UTF-8 code pages in Windows apps](https://learn.microsoft.com/en-us/windows/apps/design/globalizing/use-utf8-code-page).

## 5. [P1] Explicit `snapshot` protection occurs after the original binary-import dataset is cleared

**Locations:** [package/ado/_cio_import.ado:92](/Users/mike/Documents/GitHub/stata-ctools/package/ado/_cio_import.ado:92), [package/ado/_cio_import.ado:217](/Users/mike/Documents/GitHub/stata-ctools/package/ado/_cio_import.ado:217), [package/help/cimport.sthlp:264](/Users/mike/Documents/GitHub/stata-ctools/package/help/cimport.sthlp:264).

**Defect and impact.** The binary import wrapper accepts `snapshot` for SAS, SPSS, XPORT, DBF, and SHP, but clears the original data before loading the plugin and scanning the file. The `preserve` associated with `snapshot` occurs much later, after that destructive clear. It cannot protect the original dataset against a scan failure or a later load failure.

The help explicitly documents destructive default behavior for these binary formats, so that default is not itself the finding. The defect is accepting an explicit protection option whose description promises restoration, while taking its snapshot too late.

**Reproduction and observed result.** With one observation and `sentinel=42` in memory, importing a deliberately invalid DBF yielded:

```stata
capture noisily cimport dbase "/tmp/ctools-oct21-probes/corrupt.dbf", clear snapshot
display _rc " " _N " " c(k)
// Observed: 610 0 0
```

**Proposed fix.** When `snapshot` is explicitly requested, preserve the original dataset before any clear, plugin load, or scan; restore on every failing path and discard the snapshot only after success. Keep the documented native-style destructive default if desired. If protection for these formats is intentionally unsupported, reject `snapshot` before modifying data and explicitly narrow its documentation instead of accepting a misleading option.

**Regression test.** Test malformed input, absent plugin, injected SPI failure, and Break for every binary format with `snapshot`. Verify observations, variables, labels, characteristics, and sorted metadata remain intact after failure.

## 6. [P2] `cdestring` trims ASCII whitespace only

**Location:** [src/cdestring/cdestring_impl.c:75](/Users/mike/Documents/GitHub/stata-ctools/src/cdestring/cdestring_impl.c:75), used at lines 124–125.

**Defect and impact.** `is_trim_space` recognizes only ASCII space and characters 9–13. Native `destring` also accepts Unicode whitespace around numbers. Nonbreaking spaces from spreadsheets and web sources therefore cause valid numbers to be rejected, or converted to missing with `force`.

**Reproduction and observed result:**

```stata
clear
set obs 2
gen str20 s = uchar(160) + "123" + uchar(160)
replace s = uchar(8195) + "123" + uchar(8195) in 2
destring s, generate(native) force
cdestring s, generate(ct) force
list
// native: 123, 123
// ct:       .,   .
```

An ASCII-tab control converted correctly. Without `force`, the discrepancy prevents conversion of the affected variable.

**Proposed fix.** Trim decoded UTF-8 whitespace at both ends using the same character set as Stata's Unicode trim operation, with validation of byte boundaries. Alternatively, normalize through `ustrtrim()` in the wrapper while preserving the fast ASCII path. Do not remove Unicode whitespace inside numbers unless an explicit option requires it.

**Regression test.** Compare native and ctools on nonbreaking space, em space, narrow nonbreaking space, ideographic space, ASCII whitespace, mixed boundaries, and invalid UTF-8, with and without `force`, `replace`, and `dpcomma`.

## 7. [P2] PPML predictions with `if` fail after a no-FE fit with `d()`

**Location:** [package/ado/cpplmhdfe_p.ado:48](/Users/mike/Documents/GitHub/stata-ctools/package/ado/cpplmhdfe_p.ado:48).

**Defect and impact.** When `e(absvars)=="_cons"` and `d()` was saved, prediction executes:

```stata
quietly replace `varlist' = . if missing(`e(d)') `if' `in'
```

A caller's `if` expands into a second `if` qualifier, producing invalid syntax. `_predict` has already created the requested output with the linear predictor, so a failed mean prediction leaves a variable containing unfinished results. Mean, eta/xbd, and residual-like paths reaching this branch are affected; `xb`, `d`, and `stdp` exit earlier.

**Reproduction and observed result:**

```stata
clear
set obs 100
gen double x = _n/100
gen double y = mod(_n,4)
cpplmhdfe y x, d(dsum)
predict double mu_all, mu                 // succeeds
capture noisily predict double mu_if if _n<3, mu
// r(198): invalid syntax
predict double mu_in in 1/2, mu           // succeeds
```

**Proposed fix.** Build one sample marker with `marksample ..., novarlist`, and combine `missing(e(d))` with that marker in a single qualifier. Calculate into a temporary variable and publish the requested name only after all prediction steps succeed.

**Regression test.** Cover no-FE estimation with and without `d()`, absorbed estimation, `if`, `in`, combined qualifiers, default prediction, and every prediction type. Assert that errors do not leave a partial output variable.

## 8. [P2] New exports require hard links on POSIX filesystems

**Location:** [src/formats/ctools_output.c:62](/Users/mike/Documents/GitHub/stata-ctools/src/formats/ctools_output.c:62).

**Defect and impact.** Without `replace`, the common publication function calls `link(temporary, destination)` unconditionally. On a writable filesystem without hard-link support, the entire export can finish writing its temporary file and then fail with 603. With `replace`, the same operation uses `rename`, so behavior changes solely because of that option. This applies to the exporters using the shared helper.

**Evidence.** A native harness compiled the actual source with `link` interposed to return `ENOTSUP`. Observed: `prepare=0, commit=603, destination exists=0`. No physical exFAT volume was mounted for this audit. exFAT and FAT32 lack hard links according to [Microsoft's filesystem comparison](https://learn.microsoft.com/en-us/windows/win32/fileio/filesystem-functionality-comparison).

**Proposed fix.** Use a supported atomic no-replace rename operation when available, such as the appropriate platform-specific exclusive rename. Detect unsupported filesystems and implement an explicitly designed fallback that retains the no-overwrite guarantee. Do not use an existence check followed by ordinary rename, which can overwrite a racing writer. If a fallback cannot preserve full atomicity, document that limit and preserve the complete temporary output on publication failure.

**Regression test.** Exercise `ENOTSUP`/`EPERM`, existing destinations, racing creators, and a real exFAT volume or equivalent mount on both macOS and Linux. Compare default and `replace` behavior.

## 9. [P2] Export publication changes destination permissions to 0600

**Locations:** [src/formats/ctools_output.c:36](/Users/mike/Documents/GitHub/stata-ctools/src/formats/ctools_output.c:36), [src/formats/ctools_output.c:59](/Users/mike/Documents/GitHub/stata-ctools/src/formats/ctools_output.c:59).

**Defect and impact.** `mkstemp` creates a private temporary file. Publishing that inode without applying destination metadata leaves new exports at mode 0600 and replaces an existing file's mode with 0600. Shared-directory workflows can lose group/other access after a successful export. The implementation contains no preservation step for existing permissions or ACLs; only the mode-bit change was tested here.

**Reproduction and observed result.** A native probe created a destination with mode 0640, then called the real prepare/commit functions with `replace=true`:

```text
mode_before=640
replace prepare=0 commit=0 mode_after=600
```

A normal Stata `cexport delimited` also created a mode-0600 file.

**Proposed fix.** Define and implement final-file metadata semantics: preserve an existing destination's required mode/ACL on replacement, and apply ordinary creation permissions subject to the user's umask for a new output. Keep the temporary file private while it is incomplete. Use descriptor-based operations where possible; avoid temporarily changing the host application's process-wide umask in a threaded plugin.

**Regression test.** Check new-file umasks, replacement of group-readable files, directory-inherited ACLs, and metadata-operation failures. Verify that applying metadata does not introduce an overwrite race or expose incomplete content.

## 10. [P2] Appending the temporary suffix rejects valid filenames

**Location:** [src/formats/ctools_output.c:27](/Users/mike/Documents/GitHub/stata-ctools/src/formats/ctools_output.c:27).

**Defect and impact.** Temporary filenames append `.ctools.XXXXXX` to the entire destination basename. This consumes another 14 characters of the component/path limit. A valid final name can therefore be impossible to export through ctools even when the filesystem and native Stata support it.

**Reproduction and observed result.** On the Mac filesystem, a basename consisting of 246 `x` characters plus `.csv` is 250 characters. Native `export delimited ..., replace` succeeded; `cexport delimited` to that same path returned **603**. Its temporary basename would be 264 characters.

**Proposed fix.** Generate a short independent random temporary basename in the destination's parent directory, rather than appending to the destination basename. Keep staging on the same filesystem and use exclusive creation. Handle separators, root directories, Unicode, and Windows extended paths through the common path abstraction from finding 4.

**Regression test.** Test destination components at the filesystem's supported limit, multibyte names, near-limit full paths, replacement, and temporary cleanup. Include CSV and a binary exporter to cover the shared helper.

## 11. [P2] SHP export has an unintended 2 GiB Windows offset limit

**Locations:** [src/io/cio_shp.c:400](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:400), [src/io/cio_shp.c:497](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:497), [src/io/cio_shp.c:504](/Users/mike/Documents/GitHub/stata-ctools/src/io/cio_shp.c:504).

**Defect and impact.** The exporter stores `ftell()` positions in `long` and computes record lengths using `(long)words * 2`. Windows retains a 32-bit `long` on 64-bit builds. The code checks file lengths against a much larger word-based limit, but the position API cannot represent offsets at or beyond 2 GiB. An otherwise permitted large export therefore fails late on Windows, unlike the corresponding 64-bit POSIX path; signed length arithmetic also needs auditing.

This is source-confirmed, not a reproduced multi-gigabyte Windows export. The relevant platform facts are documented in [Microsoft's type ranges](https://learn.microsoft.com/en-us/cpp/cpp/data-type-ranges?view=msvc-170) and [ftell/_ftelli64 reference](https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/ftell-ftelli64?view=vs-2019).

**Proposed fix.** Add checked 64-bit seek/tell wrappers using `_ftelli64`/`_fseeki64` on Windows and `ftello`/`fseeko` on POSIX. Use an explicit 64-bit type for offset and record-length arithmetic. Apply the chosen SHP format/interoperability limit separately, with a clear diagnostic before publication. Audit the XLS writer's `long`/`ftell` usage too, against that format's own limits.

**Regression test.** Use mocked offsets and sparse-file fixtures around 2 GiB, plus a Windows end-to-end export with many individually small records. Ensure oversized outputs fail explicitly and consistently if a deliberate format limit is chosen.

## 12. [P1] Publication treats skipped Stata validation as sufficient

**Locations:** [.github/workflows/build.yml:210](/Users/mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:210), [.github/workflows/build.yml:234](/Users/mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:234), [.github/workflows/build.yml:272](/Users/mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:272), [validation/stata_runner.py:16](/Users/mike/Documents/GitHub/stata-ctools/validation/stata_runner.py:16).

**Defect and impact.** The licensed Stata job is skipped when `CTOOLS_STATA_RUNNER` is unset. The package condition requires only `needs.stata-validation.result != 'failure'`, which accepts `skipped`. The resulting package can be committed to the public installation directory on main or attached to a release. Thus a successful release workflow does not establish that the ado/plugin combination works in Stata. The workflow comments indicate this bypass is intentional; it remains a release-readiness gap, and conflicts with the installation guide's statement that publication requires the licensed correctness gate.

There is also no fresh-binary Windows/Linux licensed validation in this job: it downloads only `ctools_mac_*`. Merely assigning Linux/Windows runner labels would not fix that. The launcher itself hardcodes `/bin/zsh -lic` and a `stata` shell alias.

**Proposed fix.** Separate candidate artifact builds from publication. Require an explicit successful licensed validation result before updating the public install tree or creating release assets. Make missing runner configuration an actionable publication failure. Add configured executable-based launch support and validate fresh binaries for each target platform, with a documented minimum-version/OS policy. Keep untrusted PRs off persistent licensed runners.

**Regression test.** Test the workflow conditions for success, failure, skipped, and cancelled validation. A skipped job may produce a clearly labeled candidate artifact, but must not authorize release publication. Verify the plugin revision loaded by each licensed job matches its candidate build.

## 13. [P2] Canonical package-only edits bypass workflow triggers

**Location:** [.github/workflows/build.yml:6](/Users/mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:6), repeated in the pull-request path filter at line 28.

**Defect and impact.** The editable ado/help/manifest sources are now under `package/`, but both push and pull-request path filters list generated `build/*` files and omit `package/**`. A commit that changes only a canonical source, without regenerating its build copy, starts no workflow. The synchronization check that should reject that exact mistake consequently never runs. The old build copy can remain the installed behavior despite a merged source fix.

**Trigger.** Change only `package/ado/cdestring.ado`, `package/help/cimport.sthlp`, or `package/metadata/ctools.pkg`. None matches the configured workflow path patterns.

**Proposed fix.** Add `package/**` to both filters. Prefer an inexpensive always-running source/package consistency job, and make its check required, so path filtering cannot hide stale generated files.

**Regression test.** Assert that representative canonical ado, help, metadata, and new helper paths trigger the validation workflow. Verify that a deliberate unsynchronized source edit fails `sync_package.py --check` in CI.

## 14. [P2] The documented package checker accepts corrupted payloads

**Locations:** [scripts/package_contract.py:38](/Users/mike/Documents/GitHub/stata-ctools/scripts/package_contract.py:38), [scripts/stage_package.py:56](/Users/mike/Documents/GitHub/stata-ctools/scripts/stage_package.py:56), [docs/INSTALLATION.md:23](/Users/mike/Documents/GitHub/stata-ctools/docs/INSTALLATION.md:23).

**Defect and impact.** Staging writes SHA-256 hashes to `BUILD_INFO.json`, and the installation instructions tell users to validate release directories with `check_package.py`. That checker verifies names and existence only. It never reads the hash manifest, so a damaged plugin or ado file passes the documented validation command.

**Reproduction and observed result.** Stage a Mac ARM package, replace its plugin contents with the bytes `not a plugin`, leave `BUILD_INFO.json` intact, and call `check_package()` on the directory. The returned errors list was **`[]`**. The CLI consequently considers the package complete despite its payload differing from the recorded hash.

**Proposed fix.** Add an explicit release-integrity mode that requires and validates `BUILD_INFO.json`, recalculates every recorded hash, and checks that the hash inventory and platform scope match the package. Keep a separate metadata-only/source-tree mode for pre-build checks. Have the manual installation instructions and release jobs use the integrity mode. Hash checking detects damage relative to the manifest; it should not be described as independent publisher authentication.

**Regression test.** Corrupt a plugin, truncate an ado, modify the manifest, remove a checksum entry, and change the declared scope. Each staged-release validation must fail. Unmodified host and all-platform packages must pass.

## 15. [P2] Windows cleaning silently leaves stale compiled files under MSYS2

**Location:** [Makefile:584](/Users/mike/Documents/GitHub/stata-ctools/Makefile:584).

**Defect and impact.** Windows builds are documented and used in CI through MSYS2, but the Windows `clean` recipe calls `del /Q`, a `cmd.exe` builtin, through make's POSIX shell. The leading `-` suppresses recipe failures. The command prints “Done” and returns success while leaving plugins and loose object files behind. The `2>NUL` redirection can also create a literal `NUL` file under a POSIX shell.

**Reproduction and observed result.** In an isolated temporary directory, invoke the real Makefile's clean target with `DETECTED_OS=Windows`, temporary build/source directories, and POSIX `sh`. With a dummy `ctools_windows.plugin` and `stale.o`, make returned **0**, logged three ignored **127** errors, and retained both files. This tests the recipe's shell mismatch; it is not a claim of a full Windows build test.

**Proposed fix.** Use the Python cleanup already present in this target to remove all intended plugin/object files across platforms, or use MSYS-compatible `rm` consistently. Remove failure suppression for unexpected cleanup errors. Include the generic development `ctools.plugin` in the intended cleanup inventory if `clean` promises to remove all compiled plugins.

**Regression test.** Run `make clean` under the documented MSYS2 shell with stale artifacts present, assert their removal, and assert genuine permission errors produce a nonzero exit. Also verify the generic plugin from `make native` is handled consistently.

## Recommended release exit criteria

1. Fix findings 1–5 and require successful licensed publication validation from finding 12.
2. Add the targeted reproductions as regressions, especially cases where a command currently returns success with incorrect or incomplete data.
3. Run the full native and licensed suites with the actual release toolchains and fresh candidate plugins. Resolve the local OpenMP setup gaps before counting the corresponding tests as passed.
4. Execute a real Windows path/large-file/cleanup campaign and a Linux filesystem campaign. Exercise macOS Intel separately, along with the documented minimum OS and Stata versions; build metadata checks do not substitute for execution.
5. Validate the staged release's hashes, file inventory, dependency baselines, and loaded plugin revisions before publication.

## Appendix A: compact corrupt-XLSX and SHP fixtures

Run this Python code in a disposable directory. It generates only audit fixtures under `/tmp/ctools-oct21-probes`; use a fresh directory or remove earlier audit outputs before rerunning the SHP no-overwrite case.

```python
from pathlib import Path
import struct
import zipfile

p = Path('/tmp/ctools-oct21-probes')
p.mkdir(exist_ok=True)
ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
rel = 'http://schemas.openxmlformats.org/package/2006/relationships'
orel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
files = {
    '[Content_Types].xml': (
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '</Types>'
    ),
    '_rels/.rels': (
        f'<Relationships xmlns="{rel}"><Relationship Id="rId1" '
        f'Type="{orel}/officeDocument" Target="xl/workbook.xml"/></Relationships>'
    ),
    'xl/workbook.xml': (
        f'<workbook xmlns="{ns}" xmlns:r="{orel}"><sheets>'
        '<sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
    ),
    'xl/_rels/workbook.xml.rels': (
        f'<Relationships xmlns="{rel}"><Relationship Id="rId1" '
        f'Type="{orel}/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
    ),
    'xl/worksheets/sheet1.xml': (
        f'<worksheet xmlns="{ns}"><dimension ref="A1:A2"/><sheetData>'
        '<row r="1"><c r="A1"><v>123</v></c></row>'
        '<row r="2"><c r="A2"><v>456</v></c></row>'
        '</sheetData></worksheet>'
    ),
}

def write_xlsx(name, parts, level=6):
    with zipfile.ZipFile(p / name, 'w', compression=zipfile.ZIP_DEFLATED,
                         compresslevel=level) as z:
        for key, value in parts.items():
            z.writestr(key, value)

write_xlsx('valid.xlsx', files)
write_xlsx('bitflip.xlsx', files, level=0)
b = bytearray((p / 'bitflip.xlsx').read_bytes())
b[b.index(b'<v>123</v>') + 3] ^= 8
(p / 'bitflip.xlsx').write_bytes(b)
try:
    with zipfile.ZipFile(p / 'bitflip.xlsx') as z:
        z.read('xl/worksheets/sheet1.xml')
except zipfile.BadZipFile as error:
    print(error)  # Expected CRC failure

truncated = dict(files)
truncated['xl/worksheets/sheet1.xml'] = (
    files['xl/worksheets/sheet1.xml'].split('<row r="2">')[0]
    + '<row r="2"><c r="A2"><v>45'
)
write_xlsx('truncated-xml.xlsx', truncated)
with zipfile.ZipFile(p / 'truncated-xml.xlsx') as z:
    assert z.testzip() is None  # ZIP is intact; contained XML is incomplete

h = bytearray(100)
struct.pack_into('>I', h, 0, 9994)
struct.pack_into('>I', h, 24, 64)
struct.pack_into('<II', h, 28, 1000, 1)
struct.pack_into('<dddd', h, 36, 1, 2, 1, 2)
(p / 'point.shp').write_bytes(
    h + struct.pack('>II', 1, 10) + struct.pack('<Idd', 1, 1, 2)
)
(p / 'partial.shx').write_text('old index sentinel')
(p / 'corrupt.dbf').write_bytes(b'not a dbf')
```

In Stata, use the candidate package's ado files and plugin, then:

```stata
capture noisily cimport excel "/tmp/ctools-oct21-probes/bitflip.xlsx", clear
local rc = _rc
display "bitflip rc=" `rc'
list
// Current bug: success, A[1] == 923

capture noisily cimport excel "/tmp/ctools-oct21-probes/truncated-xml.xlsx", clear
local rc = _rc
display "truncated XML rc=" `rc'
list
// Current bug: success, second row missing

cimport shp "/tmp/ctools-oct21-probes/point.shp", clear
capture noisily cexport shp using "/tmp/ctools-oct21-probes/partial.shp", shx
local rc = _rc
display "SHP export rc=" `rc'
confirm file "/tmp/ctools-oct21-probes/partial.shp"
// Current bug: export failed with 602 but SHP was published
```

## Appendix B: audit evidence and limitations

Local diagnostic artifacts were left outside the repository's tracked sources:

- `/tmp/ctools-oct21-native.log`: all 35 native suite attempts and the 11 OpenMP setup/build failures.
- `/tmp/ctools-oct21-build.log`: successful fresh source build.
- `/tmp/ctools-oct21-probes/stata.log`: fresh-plugin identity, filename, SHP, destring, and snapshot results.
- `/tmp/ctools-oct21-probes/xlsx.log`: one-bit corruption and truncated XML imports.
- `/tmp/ctools-oct21-probes/predict.log`: PPML prediction qualifier failure.
- `/tmp/ctools-oct21-output.c`: focused actual-source hard-link failure and permission probe.

These `/tmp` artifacts are ephemeral; the reproductions and observed outputs above are the durable record. Findings marked source-confirmed still need their listed target-platform regression runs. The review also examined aligned-allocation pairing, worker teardown/fallback, sorting, and numerical code, but absence of a finding in those areas is not a claim of exhaustive verification.
