# ctools release audit — November 16, 2026

## Implementation follow-up — November 16, 2026

RA-01 through RA-04 are fixed in source, with regression coverage. RA-05
(release publication) and RA-06 (metadata-tool encoding) were left unchanged
at the user's request. The original findings and measurements below describe
the audited candidate before these fixes.

| Finding | Implemented change | Verification |
| --- | --- | --- |
| RA-01 | Predict structural values over the union of the estimation and prediction samples, so the partial projection always uses the original estimation sample. | Filtered, `in`, and out-of-sample prediction invariance; weights, factor variables, partialled constants, missing values, and unchanged estimation results. The full model supplies the weighted/FV residual oracle because `ivreg2_p` itself omits estimation weights in its partial-projection branch. |
| RA-02 | Preserve literal text through Excel dispatch, filename parsing, worksheet handling, plugin metadata, and returned filenames. Replace the positional-filename one-line `if`, whose second parse re-expanded protected text. | `.xls` and `.xlsx`, `using` and positional filenames, dollar/backquote characters, worksheet modify/replace, exact returned names, and preservation of a colliding existing workbook. |
| RA-03 | Reject invalid decimal precision/scale, physical-type capacity violations, and conflicting annotations during schema traversal. Widen exponent arithmetic to 64 bits. | 18 added valid decimal fixtures and 29 invalid fixtures; complete native corpus with allocation failures and UBSan; all 48 malformed files preserve data in default and snapshot modes. Validation follows the [Parquet DECIMAL specification](https://github.com/apache/parquet-format/blob/master/LogicalTypes.md#decimal) within the reader's decimal256 capacity. |
| RA-04 | Set x87 extended precision consistently inside each interpolation worker and restore its prior control word afterward. Other architectures retain their existing arithmetic. | Confirmed main/worker control words `0x027f`/`0x037f` in Windows Stata. Exact equality across serial and parallel runs; reused worker pools; restoration of main/worker settings. The new native regression fails against the original source. |

Validation against a locally rebuilt Windows OpenMP plugin passed: `cimport`
536 checks, `cexport` 247, `civreghdfe` 480, `cinterpolate` 65, grouped IV
syntax 17, Parquet parity 172, malformed-file rollback 96, and all 345 selected
option/data cases. Native serial/OpenMP command tests and the complete Parquet
corpus passed with UBSan in trap mode, including 1,644 byte mutations and
truncations. A deliberate overflow negative control verified that UBSan traps
were active. Seed, public-option inventory, documentation, and diff whitespace
checks passed.

The test build and logs are under `temp/audit-fixes/`; the plugin is
`temp/audit-fixes/candidate/ctools_windows.plugin`, built with GCC 16.2.0 and
identified as revision `audit-fixes`. It is a local unsigned test artifact.
Checked-in distribution plugins, signing/publication workflows, and release
metadata were not changed. Linux/macOS and older Stata runtime qualification
remain outside the verification performed here.

## Release assessment

**Address the two P1 defects before release:** filtered residual predictions after `civreghdfe, partial()` are incorrect, and Excel export can write to an unintended filename, including overwriting an existing file. Both were reproduced against the checked-in Windows plugin and compared with the corresponding Stata/reference command.

The audit also confirmed malformed Parquet decimal handling, a failing Windows interpolation reproducibility test, and a Windows Python encoding failure. The release workflow still lacks a requirement for licensed Stata results tied to the final candidate artifacts. The six findings below include proposed fixes and acceptance checks; production code was not changed.

| ID | Priority | Finding | Evidence |
| --- | --- | --- | --- |
| RA-01 | P1 | `civreghdfe` changes partialled residuals when prediction is filtered | Runtime reproduction and source explanation |
| RA-02 | P1 | `cexport excel` reinterprets literal macro characters and can overwrite another destination | Runtime reproduction, workbook inspection, overwrite demonstration |
| RA-03 | P2 | Invalid Parquet decimal metadata produces successful imports with corrupt strings | Runtime reproduction; unchecked arithmetic identified in source |
| RA-04 | P2 | Signed Windows plugin fails the existing cubic interpolation thread-invariance assertion | Existing suite failure and isolated reproduction |
| RA-05 | P2 | Publication does not require licensed Stata validation of candidate artifacts | Workflow and validator inspection; corroborated by RA-04 |
| RA-06 | P3 | Direct release-metadata checks fail under Windows' non-UTF-8 Python defaults | Runtime reproduction |

P1 means fix before release because normal command use can produce wrong results or affect the wrong file. P2 means resolve before release, or explicitly document and accept a bounded limitation. P3 is a development/release-tool compatibility issue with a straightforward workaround. RA-04 concerns small floating-point differences and a failing contract, not demonstrated economically or statistically material error.

## Audited candidate and method

| Property | Value |
| --- | --- |
| Checkout HEAD | `0388e3075907dc080e29849ffcf639d971b9ced8` |
| Release configuration | Version `1.2.6`; date `20261111`; protocol `3`; minimum Stata `14.1` |
| Loaded plugin identity | `1.2.6`, revision `9eadc285273add6408b97d27cf9517bfa714b147`, protocol `3` |
| Windows artifact | `build/ctools_windows.plugin`; Authenticode status `Valid` |
| Artifact SHA-256 | `14b6e1b25fb72598edba2d6f80fb731df05560d8d720abc3f96d3ba661883fa7` |
| Runtime | StataNow/MP 19.5, Windows x86-64, revision October 14, 2026 |
| Relevant reference packages | `ivreg2` 4.1.12; `ivreghdfe` 1.1.4; `reghdfe` 6.13.1; `ppmlhdfe` 2.3.3 |
| Python | 3.12; UTF-8 mode enabled for the broad validation runs |

The checkout commit and the embedded plugin revision are recorded separately; the former is a compiled-plugin update commit. Their difference alone is not reported as a bug.

Review combined source inspection of ado dispatch/prediction, native parsing and numerical paths, release scripts and CI, existing licensed Stata suites, and targeted comparisons against native/reference behavior. All 22 public computational command components were exercised. The prior `ASTRA_NOV6.md` audit and documented compatibility exceptions were consulted to avoid reopening repaired issues or reporting intentional parity behavior as new defects.

### RA-01 — Filtered partialled residual predictions use the wrong projection sample

**Priority: P1. Confidence: high; reproduced.**

**Location:** `build/civreghdfe_p.ado:82–95`, in the `partial_ct > 0` residual-prediction branch. The relevant scope is `partial()` without absorbed fixed effects; the saved-residual branch with `absorb()` follows a different path.

**Trigger and observed behavior:** Estimate `civreghdfe y x z, partial(z)`, then request residuals with an `if` or `in` qualifier. Residuals for the selected observations differ substantially from those returned by unrestricted prediction after the same fit. Requesting only out-of-estimation-sample observations returns `r(2000)` even though unrestricted prediction and `ivreg2` can predict those observations.

Minimal reproduction, with the repository's `build` directory on `adopath`:

```stata
clear
set seed 123
set obs 80
gen long g = ceil(_n/8)
gen double x = rnormal()
gen double z = rnormal()
gen double y = 3*x + 2*z + g + rnormal()

quietly civreghdfe y x z, partial(z)
predict double full, residuals
predict double subset if _n <= 10, residuals
gen double difference = subset - full
summarize difference
```

The ten differences range from **3.580779 to 4.860384**, with mean **4.110960**. For the first observation, unrestricted prediction is approximately `-3.8330298` and filtered prediction is `0.28321917`. Repeating the two predictions after `ivreg2 y x z, partial(z)` gives exactly zero difference on the selected observations. The unrestricted ctools and reference residuals agree to ordinary numerical precision.

A second reproduction estimates on `_n <= 60` and predicts only for `_n > 60`: ctools returns `r(2000)`; `ivreg2` succeeds, and its filtered predictions equal its unrestricted predictions on those rows.

**Cause:** `PredictXB` is called with the prediction qualifiers. Consequently, `ypart = y - xb` is missing outside the prediction subset. Although the subsequent auxiliary regression specifies the saved `e(sample)`, Stata drops those missing observations and estimates the partial projection on the intersection of the estimation and prediction samples. A disjoint prediction sample leaves no observations for that regression.

**Proposed fix:** Compute the structural prediction needed for residualization over at least the union of the original estimation sample and the requested prediction sample. Fit the partial projection on the original estimation sample with the original weights and constant convention. Apply prediction qualifiers only when publishing the requested residuals. Preserve the original estimation results on success and failure.

**Acceptance checks:** After one fit, unrestricted residuals must equal residuals requested with `if`, `in`, and disjoint out-of-sample restrictions wherever both are defined. Cover a partialled constant, factor-variable terms, weights, missing predictors, and observations outside `e(sample)`. Compare with `ivreg2` and assert that prediction leaves the original `e()` results and sample intact. Existing full-sample prediction tests do not establish this invariance.

**Local evidence:** `temp/audit-nov16/probes.do` / `.log`; the out-of-sample extension is in `temp/audit-nov16/final_probes.do` / `.log`.

### RA-02 — Excel export re-expands literal filename and worksheet text

**Priority: P1. Confidence: high; reproduced, including unintended overwrite.**

**Location:** `build/cexport.ado:34` forwards the remaining command as ``cexport_excel `0'``. Additional ordinary macro expansions appear in the Excel handler, including the sheet result at line 472 and filename handling beginning around line 491. The import dispatcher already uses `macval()` to avoid analogous re-expansion.

**Trigger and observed behavior:** Pass an Excel filename or sheet name containing literal Stata macro characters through `macval()`. `export excel` preserves the literal text; `cexport excel` expands it again. With global `AUDITliteral` set to `expanded`, a filename ending in `ctools$AUDITliteral.xlsx` becomes `ctoolsexpanded.xlsx`, and sheet `sales$AUDITliteral` becomes `salesexpanded`.

Run this in a disposable directory:

```stata
clear
set obs 2
gen x = _n
global AUDITliteral expanded
mata: st_local("p", "native" + char(36) + "AUDITliteral.xlsx")
mata: st_local("sheet", "sales" + char(36) + "AUDITliteral")
export excel using `"`macval(p)'"', sheet(`"`macval(sheet)'"') replace

mata: st_local("p", "ctools" + char(36) + "AUDITliteral.xlsx")
cexport excel using `"`macval(p)'"', sheet(`"`macval(sheet)'"') replace
macro drop AUDITliteral
```

Inspection of the actual filenames and `xl/workbook.xml` confirmed both differences. In a separate collision test, an existing `collisionexpanded.xlsx` contained `999`. Exporting `1` to the literal intended filename `collision$AUDITliteral.xlsx` with `replace` overwrote `collisionexpanded.xlsx`; reading it back returned `1`.

**Impact:** The command can create the wrong output file or worksheet. With `replace`, it can overwrite a different existing workbook from the one the caller supplied. This is not merely a cosmetic output-message difference.

**Cause:** Ordinary macro substitution reparses literal user text as Stata macro syntax at command-forwarding and later interpolation boundaries. Protecting only the initial caller expression is insufficient when the wrapper expands it again.

**Proposed fix:** Preserve literal text end to end. Use `macval()` when forwarding command arguments and appropriately protected locals, and retrieve/store text through Mata's `st_local()` interface where feasible instead of inserting it into executable command text. Audit the whole Excel path, including positional filenames, sheet parsing, return values, existence checks, update modes, and plugin argument transport; changing only the first dispatcher line may leave later re-expansion points.

**Acceptance checks:** Compare exact output paths and worksheet names against `export excel` for literal dollar signs, defined/undefined global names, backquote/apostrophe macro sequences, spaces, quotes, and Unicode. Cover `.xls` and `.xlsx`, `replace`, `sheetmodify`, and `sheetreplace`. Place a sentinel workbook at the expanded-name collision destination and assert it is unchanged. Verify workbook contents and XML, not just `_rc` or messages.

**Scope:** The confirmed native-parity difference is Excel export. A CSV probe did not establish the same discrepancy because native CSV export also expanded the tested filename; that case is not included in this finding.

**Local evidence:** `temp/audit-nov16/export_literals.do` / `.log`, generated workbooks in the same directory, and the collision test in `final_probes.do` / `.log`.

### RA-03 — Malformed Parquet decimal schemas are accepted and converted incorrectly

**Priority: P2. Confidence: high for incorrect acceptance/output; integer-overflow risk confirmed by source inspection, not by a successful sanitizer run.**

**Locations:**

- `src/io/parquet_thrift.c:437–442`: `primitive_type()` copies decimal precision/scale and selects decimal128/decimal256 without validating their supported ranges there.
- `src/io/parquet_decode.c:562`: `int adjusted = nd - 1 - scale` can overflow for `scale == INT_MIN`.
- `src/io/cparquet.c:190`: `precision + 2` can overflow and does not establish a safe formatted width for impossible metadata.

**Trigger and observed behavior:** Import an otherwise readable three-row Parquet file with INT32 physical values `1, 2, 3` and invalid DECIMAL metadata. Files were generated using the repository's own `validation.support.parquet_writer`, avoiding a dependency on an external writer accepting these invalid schemas.

| Precision | Scale | Native `import parquet` | `cimport parquet` |
| ---: | ---: | --- | --- |
| 10 | -2147483648 | Rejects schema, `r(3900)` | Success; strings `1E-214748364`, `2E-214748364`, `3E-214748364` |
| 10 | 2147483647 | Rejects schema, `r(3900)` | Success; the same truncated exponent strings |
| 0 | 0 | Rejects schema, `r(3900)` | Success; strings `1`, `2`, `3` |
| 2147483647 | 0 | Rejects schema, `r(3900)` | Later `r(3598)` / “invalid syntax”; original sentinel data restored |

The successful cases used `clear` and replaced the original dataset. Returning success with these strings can conceal a corrupt upstream schema. The maximum-precision case demonstrates a later, less appropriate failure, rather than silent success.

Fixture recipe, run from the repository root:

```python
from pathlib import Path
from validation.support.parquet_writer import Column, INT32, logical, write_table

out = Path("temp/decimal-audit")
out.mkdir(parents=True, exist_ok=True)
for name, precision, scale in [
    ("scale_min", 10, -2147483648),
    ("scale_max", 10, 2147483647),
    ("precision_zero", 0, 0),
    ("precision_max", 2147483647, 0),
]:
    write_table(out / (name + ".parquet"), [
        Column("d", INT32, [1, 2, 3],
               logical_type=logical("DECIMAL", precision=precision, scale=scale))
    ])
```

Compare `import parquet using "temp/decimal-audit/scale_min.parquet", clear` with `cimport parquet` on the same file, and repeat for the other three files.

**Proposed fix:** Validate decimal annotations before selecting the output type or decoding values. Enforce the supported precision/scale domain and its compatibility with the physical representation, including converted and logical annotations. Reject inconsistent metadata with a schema error while retaining the original dataset. Use checked/wider arithmetic defensively for exponent and width calculations, and check formatted lengths before narrowing or publishing strings. Preserve deliberately supported valid decimal cases rather than imposing an unverified blanket rule on all negative scales.

**Acceptance checks:** Add complete-file malformed fixtures for zero/negative/excessive precision, extreme scales, incompatible physical widths/types, and conflicting annotations. Include boundary-valid schemas. Run the decoder under UBSan and ASan where supported. In licensed Stata, require explicit failure and data restoration for invalid schemas; verify full values and storage widths for valid schemas. This audit does not establish an exploitable memory-safety vulnerability.

**Local evidence:** `temp/audit-nov16/parquet_edges.py`, `parquet_invalid.do` / `.log`, and the four generated `.parquet` files.

### RA-04 — Cubic interpolation fails the Windows thread-invariance contract

**Priority: P2. Confidence: high for the failure; root cause unresolved.**

**Locations:** The existing assertion is `validation/stata/cinterpolate.do:253`, `assert got == parallel`. The production cubic calculation uses `long double` in `src/cinterpolate/cinterpolate_impl.c:23–39`; its grouped OpenMP loop begins at line 115.

**Observed behavior:** The unmodified `cinterpolate` component aborts with `r(9)` after 57 counted passes. A trace copy and a separate reproduction identify the serial-versus-parallel cubic equality assertion. The component completion marker is absent, so its `FAILED=0` counter does **not** mean the suite passed; the runner correctly rejects the incomplete result.

Reproduction:

```stata
clear
set seed 87533
set obs 30000
gen byte g = mod(_n, 30)
gen str4 group = "g" + string(g)
gen double x = floor((_n-1)/30) + runiform()/2
gen double y = 2 + x/10 - x^2/1000 + x^3/1000000 + g
replace y = . if mod(floor((_n-1)/30), 3) == 1
cinterpolate y x, gen(serial) by(group) cubic threads(1)
cinterpolate y x, gen(parallel) by(group) cubic threads(4)
count if serial != parallel
gen double difference = abs(serial - parallel)
summarize difference
```

There are **5,869 unequal observations** out of 30,000. The maximum absolute difference is approximately **5.68e-14**, and the maximum relative difference, using the larger absolute result as denominator, is approximately **3.88e-14**. The same count and maxima occurred with 16 threads. Missing/nonmissing status did not differ. These are small rounding differences, but an existing release test requires exact equality.

**Proposed fix:** First determine why the actual Windows artifact differs across thread counts. Inspect main-thread/worker floating-point control state, Windows extended-precision behavior, and compiler/OpenMP code generation around the `long double` calculations. These are hypotheses, not established causes. If exact thread invariance is intended, ensure consistent arithmetic while preserving the host's floating-point environment. If the supported contract permits bounded numerical differences, explicitly document that decision and replace exact equality with a justified error bound while retaining analytic-oracle, missingness, sample, and ordering checks.

**Acceptance checks:** Run the unmodified component through completion against the rebuilt signed artifact. Cover one, four, and sixteen threads, fresh and reused worker pools, numeric/string groups, duplicates, and ill-conditioned anchor spacing. Repeat on each supported platform; a local serial build is insufficient to clear this Windows artifact failure. Do not simply relax every interpolation tolerance.

**Local evidence:** `temp/audit-nov16/baseline2/stata.log`, `interpolate_trace.do` / `.log`, and `cubic_threads.do` / `.log`.

### RA-05 — Release publication is not conditional on candidate Stata results

**Priority: P2. Confidence: high; source-confirmed process gap.**

**Locations:** `.github/workflows/build.yml:74` runs native validation in the Linux build; the package job at line 377 depends on builds/signing; the release job at line 475 depends on packaging. `validation/checks.py:72–73` explicitly asserts that `stata-validation:` and `CTOOLS_STATA_RUNNERS` are absent from the workflow.

**Finding:** Build, dependency, signing, package, and native-test checks do not demonstrate successful ado/plugin behavior inside licensed Stata. The inspected publication path does not require that evidence for the final artifacts. RA-04 is a concrete example: the signed, loadable Windows plugin fails an existing licensed test that the publication dependencies do not run. RA-01 and RA-02 also show gaps that native-only tests cannot close.

This continues the unresolved candidate-qualification concern in `ASTRA_NOV6.md` (ASTRA-49); it is not presented as a newly introduced regression or proof that nobody has performed manual testing.

**Proposed fix:** Require licensed validation evidence before publication, either from suitable runners or a retained manual qualification record checked by the release process. Bind that record to the exact staged ado/plugin package hashes, embedded build revision, platform, Stata version, reference-package versions, completion markers, and failures/skips. Validate the final signed artifacts where applicable. Update the workflow-contract assertions if the gate design changes.

**Acceptance checks:** A missing, failed, interrupted, or stale result must block publication. Changing an artifact after qualification must invalidate its record. Cover Windows, Linux, macOS Intel and Apple Silicon, with the advertised minimum Stata version and a current version; include checks on advertised OS/dependency floors. A smaller documented matrix can be used only with an explicit narrowing of the compatibility claim or a reasoned qualification plan.

**Relationship to this audit:** This audit supplies current Windows/StataNow evidence, not a complete cross-platform release certification. Native sanitizer execution and the other platform/version combinations remain unverified here.

### RA-06 — Release metadata tools depend on the caller enabling UTF-8 mode

**Priority: P3. Confidence: high; reproduced.**

**Location:** `scripts/sync_release.py:41` calls `Path.read_text()` without an encoding; similar defaults occur at lines 11, 70, 74, 92, 108, and 110. `scripts/build.mk:39` exports `PYTHONUTF8=1`, but direct documented Python commands do not inherit that Make-only setting. See `validation/README.md:48`.

**Trigger and observed behavior:** Run the documented check with Windows Python 3.12 using its cp1252 text default:

```powershell
python -X utf8=0 scripts/sync_release.py --check
```

It fails at line 41 with `UnicodeDecodeError: 'charmap' codec can't decode byte 0x90 in position 119`. UTF-8 mode allows the same checkout to pass. The failure occurs while reading repository text, before the script can establish whether release metadata is current. This affects contributor/release tooling, not installed Stata command execution.

**Proposed fix:** Specify `encoding="utf-8"` for repository text reads and writes in the synchronization tool and audit related release/check scripts for the same assumption. Select an appropriate explicit encoding separately for external process logs. Keep the environment workaround for convenience, but do not rely on it as the script's correctness condition.

**Acceptance checks:** Run `--check` on native Windows Python with UTF-8 mode disabled and enabled. In an isolated checkout, run synchronization and verify that UTF-8 content remains unchanged except for intended metadata edits, followed by a successful `--check`. This also protects against accidental re-encoding when synchronization writes files.

**Workaround:** Set `PYTHONUTF8=1` or use `python -X utf8 scripts/sync_release.py --check` until fixed.

**Local evidence:** `temp/audit-nov16/encoding-default.log`; successful metadata checks were run with UTF-8 enabled.

## Validation results and limits

### Existing licensed Stata components

The following are this run's existing suite counters, not counts of newly added tests:

| Component | Result |
| --- | --- |
| `csort` | 380 checks passed |
| `cmerge` | 330 passed |
| `cdecode` | 229 passed |
| `cencode` | 248 passed |
| `cdestring` | 229 passed |
| `csample` | 70 passed |
| `cbsample` | 99 passed |
| `cwinsor` | 134 passed |
| `cinterpolate` | 57 passed before assertion failure; incomplete, `r(9)`; RA-04 |
| `csplit` | 152 passed |
| `cimport` | 536 passed |
| `cexport` | 243 passed |
| `creghdfe` | 385 passed |
| `cqreg` | 300 passed |
| `civreghdfe` | 467 passed |
| `cbinscatter` | 281 passed |
| `cpsmatch` | 365 passed |
| `crangestat` | 306 passed |
| `cppmlhdfe` | 1,302 passed |
| `crangejoin` | 102 passed |
| `cdistinct` | Standalone completion check passed |
| `chdfe` | Standalone completion check passed |

Additional completed cases: sort identity (108), encode labels (25), winsor transport (160), grouped windows (257), IV grouped syntax (17), PPML publication (8), quantile sample/scale (371), and HDFE replay (8). Standalone winsor rounding, quantile performance, and propensity-matching statistics cases also completed successfully. These passes do not cover the new edge cases in RA-01 through RA-03.

All **1,464 executable option/data cases passed**: 561/561 baseline cases, 150/150 import/export cases, and 753/753 estimator cases. Across the 22 public command components, 21 completed successfully and `cinterpolate` aborted at the assertion described in RA-04.

Release metadata, architecture synchronization, package/build metadata checks, documentation checks, public-option inventory, random-seed checks, and the 24 offline dataset integrity checks passed with UTF-8 mode enabled. The option inventory covered 526 public options; inventory coverage is distinct from behavioral validation.

### Reproduction environment and retained evidence

The existing suites were invoked from the repository root in these groups:

```text
python validation/run.py --stata csplit cinterpolate cwinsor cbsample csample cmerge csort cdestring cencode cdecode --output temp/audit-nov16/baseline2
python validation/run.py --stata cimport cexport --output temp/audit-nov16/io
python validation/run.py --stata creghdfe civreghdfe cppmlhdfe cqreg crangestat cbinscatter cpsmatch chdfe cdistinct crangejoin --output temp/audit-nov16/estimators
```

The configured `CTOOLS_STATA_EXECUTABLE` was the local StataNow 19 executable. `TEMP`, `TMP`, and `STATATMP` were directed to a short workspace scratch path, and `PYTHONUTF8=1` was set. An initial run with the inherited sandbox temporary path did not produce a valid Stata completion log; it is excluded from product findings. Options output includes its own results and provenance files.

Raw logs, probes, fixtures, and generated workbooks are retained under `temp/audit-nov16/`, which is ignored by Git. The report embeds the essential reproductions so the findings do not depend on retaining that scratch directory. Some probe logs include exploratory cases that were ruled out; they are not additional findings.

### Limits and release follow-up

- Only the checked-in Windows plugin was exercised inside Stata during this audit. Linux, macOS, and Stata 14.1 runtime behavior were not verified.
- An attempted native/sanitizer build could not complete because the available Zig toolchain failed internally while preparing compiler support files. No native sanitizer pass is claimed, and that toolchain failure is not attributed to ctools.
- This was a source and behavioral audit, not an exhaustive proof of memory safety or numerical correctness. The detailed command suites provide substantial coverage, but the targeted probes demonstrate remaining blind spots.
- Floating-point accumulation and Parquet negative-boundary behaviors that matched the reference or documented compatibility policy were excluded from the findings.
- The pre-existing untracked `a.out` was left untouched. No production source, ado, plugin, release configuration, or existing validation test was modified.

Recommended release sequence: fix RA-01 and RA-02, reject the invalid schemas in RA-03, resolve the numerical contract and Windows failure in RA-04, make the low-risk encoding correction in RA-06, then qualify and retain results for the exact final artifacts as described in RA-05.
