"""Render the concise report from completed, frozen final comparison CSVs."""
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
A='benchmarks/transport_adaptive_20260929'
rows=list(csv.DictReader((ROOT/'final_confirmed/summary.csv').open()))
lookup={(r['case'],r['phase']):r for r in rows}
labels=[('tiny_n100_k20','100 × 20 numeric'),('wide_byte','200K × 64 byte'),
 ('wide_byte128','200K × 128 byte'),('wide_cycle','200K × 256 cycling numeric'),
 ('wide_million','1M × 128 cycling numeric'),('wide_short','50K × 128 cycling numeric'),
 ('filtered_dense','500K × 128, two-thirds selected'),('maxstr_full','200K × 4 full-length str2045'),
 ('maxstr_short','200K × 4 str2045, eight-byte values'),('numeric_control','1M × 20 numeric control'),
 ('mixed_w32','500K × 20 mixed str32'),('large_sorted_control','1M × 20 numeric, sorted writes')]
def row(case,label,lookup,phase='transport'):
 r=lookup[(case,phase)]
 result = f"| {label} | {1000*float(r['baseline_median']):.3f} | {1000*float(r['candidate_median']):.3f} | {float(r['paired_speedup_median']):.2f}× | {float(r['paired_speedup_min']):.2f}–{float(r['paired_speedup_max']):.2f}× |"
 if phase == 'transport':
  result += f" {float(lookup[(case,'load')]['paired_speedup_median']):.2f}× | {float(lookup[(case,'store')]['paired_speedup_median']):.2f}× |"
 return result
table='\n'.join(row(c,l,lookup) for c,l in labels)
command_rows=list(csv.DictReader((ROOT/'commands_final/summary.csv').open()))
clookup={(r['case'],r['phase']):r for r in command_rows}
commands='\n'.join(row(c,l,clookup,'elapsed') for c,l in [
 ('tiny_numeric','100 × 20 numeric'),('narrow_numeric','100K × 20 numeric'),
 ('wide_numeric','200K × 128 numeric'),('str2045_full','100K × 4, ID plus three str2045')])
text=f'''# Adaptive Stata–C transport

The final implementation substantially improves tiny transfers, wide numeric loads, and full-length `str2045`. Final mixed-string, narrow numeric, and large sorted-write controls are near parity. Some store phases still lose performance, and the evidence does **not** establish improvement for every workload.

**Completed:** exact-final-source transport and whole-command benchmarks, same-binary controls, full command acceptance, native UBSan checks, and installation of the tested plugin. The completed measurements below use the frozen production candidate; their provenance will not be replaced by follow-up results.

## Corrected implementation: measured transport

These measurements use the exact final K ≥ 128, selected N ≥ 50,000 policy and final frozen source. Times are milliseconds for load + store + cleanup, including destruction of the production pthread pool. Two matched process blocks each use five timed repetitions after a warmup, plus a separate fully verified pass: **336 records across four processes**. Speedup is the median paired baseline/candidate ratio; ranges are the two block ratios, not confidence intervals. Marginal median times therefore need not divide to the reported paired ratio.

| Workload | Baseline ms | Candidate ms | Total speedup | Block range | Load speedup | Store speedup |
|---|---:|---:|---:|---:|---:|---:|
{table}

The 1M × 128 and 200K × 256 numeric **load** phases improved 4.18× and 3.81×. Full-length str2045 loads improved 1.74×, and cleanup fell from roughly 24 ms to 0.09 ms by avoiding individual fallback allocations. Dense-filtered stores improved 1.19× after vectorizing map validation. The byte64 case now uses column loading; its old 0.75× same-binary tile result motivated the higher gate.

Store phases still cost more in two changed-path cases: 1M × 128 cycling numeric was 0.95× in both blocks, despite 2.21× total transport; 50K × 128 store was 0.78×, with uncertain total transport near parity (0.84–1.15×). The earlier three-block comparison had mixed-str32 and sorted-control total ratios of 0.85× and 0.92×; those losses did not repeat here. Both comparisons remain evidence of variability rather than permission to discard inconvenient results. [Complete final phase results]({A}/final_confirmed/summary.csv) and [process medians]({A}/final_confirmed/process.csv) include every case; the [earlier 64-column-gate comparison]({A}/final_comparison/summary.csv) is preserved as superseded evidence.

A same-binary follow-up (312 records, two process blocks) puts mixed str32 at 0.997×, large sorted writes at 0.998×, and numeric controls at 1.003× using block-level ratios. Mixed str8 remains slightly slower at **0.972× (0.958–0.986×)**. Shared runtime state reduces confounding, but the two compiled code bodies retain different placement and use the final shared helpers. [Block-level control results]({A}/acceptance/dual_impl_control/block_summary.csv) preserve this distinction.

## Whole-command confirmation

`csort id, verbose threads(8) nostream` was measured with one warmup and seven timed calls in each of two reversed process blocks. A monotonic clock includes ado/plugin dispatch and sorting; fixture loading and validation are outside timing. Each result is checked cell by cell and against native Stata's sorted datasignature. The four processes produced **128 records** including warmups.

| Whole csort workload | Baseline ms | Candidate ms | Paired speedup | Block range |
|---|---:|---:|---:|---:|
{commands}

All four whole-command medians improved in both final blocks, including a modest 1.05× narrow numeric gain. The earlier narrow comparison ranged from 0.69× to 1.64× and remains in the archive as evidence of run-to-run variability. These four commands validate that transport gains reach an actual command; they are not estimates for all ctools commands. [Whole-command results]({A}/acceptance/commands_final/summary.csv) preserve both blocks and all phases.

## Selected implementation

- **Small transfers:** estimate work by numeric columns and declared string widths, then use serial load/store when pool setup would dominate. Unknown widths receive a conservative cost. A plain cell-count threshold was rejected because string-heavy cases regressed sharply.
- **Numeric loads:** 512-row tiles for all-numeric varlists with **K ≥ 128 and selected N ≥ 50,000**; no separate cell-count gate. A same-binary scheduler switch showed 64-column byte loads losing, while all five numeric storage types gained at 128/256 columns. Follow-up tests also found losses for compact byte data at two million cells and for very short, wide byte data above six million cells. The final row threshold avoids all measured short-byte losses; it gives up gains for some 64-column int/long/float/double and short-wide cycling datasets.
- **Wide strings:** reserve a larger packed arena for known str2045 values, retain actual-length density, and retry the smaller arena if reservation fails. Flat declared-width slots were rejected because short values in wide declarations lost performance.
- **Unused identity order:** 11 nonordering callers request `CTOOLS_LOAD_NO_SORT_ORDER`. The observed allocation falls from **4N bytes to zero**: 40 MB at 10M observations. A same-binary toggle supports modest numeric gains; strL is approximately unchanged and mixed-string timing is uncertain. This is an allocation saving, not a measured equal fall in peak RSS.
- **Numeric writes:** validate destination bounds and uniqueness before parallel writes; vectorize the common increasing-map check and retain serial semantics for duplicates. Direct sorted gathers apply only when the full source column fits 262,144 doubles (2 MiB); larger sources retain staging.

General store tiling, interleaved scheduling, unconditional OpenMP stores, unrestricted direct gathers, and unweighted small-transfer dispatch were rejected or restricted after regressions. The [detailed evidence]({A}/report_draft.md) records those results, including controls and rejected variants.

## Correctness and profiling

The full baseline, earlier candidate, and accepted final candidate suites each passed **4,161 checks with zero failures and 48 expected method-difference skips**. Four native transport suites passed with UBSan in serial and OpenMP builds: **eight successful runs**, including bounds, allocation failures, strL, sorted/duplicate maps, cancellation, and dispatch boundaries. The ARM plugin passed dependency checks with deployment target macOS 11 and compatible static libomp 21.1.8. These are host build/check results, not new cross-platform runtime measurements.

The initial ASan runtime deadlocked during initialization before `main`; this is not a passing sanitizer run. An independent Homebrew21.1.8 ASan probe also timed out before `main`; neither run tested transport code. Added Linux CI coverage is not a claimed CI pass. [Native results and startup diagnosis]({A}/acceptance/final_native_checks_final/ubsan_results.json), [baseline acceptance]({A}/acceptance/production_baseline/acceptance.log), [candidate acceptance]({A}/acceptance/production_final/acceptance.log), and [dependency metadata]({A}/build_checks_final.json) preserve the evidence.

Separate five-second CPU samples confirm the transition from column-worker loads to tiled OpenMP loads through the same Stata callback. Stores still use the checked callback. Sampling counts pool running and sleeping threads; they are not CPU percentages or hardware cache-miss measurements and do not prove a spin-wait explanation. Profiled timings are excluded from the benchmark summaries.

## Scope, uncertainty, and reproduction

Eight exploratory component campaigns cover N=1–10M and K=1–2,000, all five numeric storage types, mixed fixed strings, bounded textual strL, filtered selections, sorted writes, subsets of a wider host dataset, and 1/8/12 requested threads. Their **9,506 records in 48 processes** include warmups and verification. Subsequent independent integration comparisons and same-binary toggles test the selected changes. This is a targeted matrix, not all possible configurations or a representative workload average.

All Stata launches used `oldstata` through interactive login zsh. Transport uses `CLOCK_MONOTONIC`; it excludes data generation, verification, and ado overhead. The host was Apple M4 Pro, 48 GiB RAM, macOS 26.7, StataNow/MP 18.5. Another user-authorized Stata job ran concurrently. Sequential execution, alternating order, process-level medians, and unchanged controls reduce some confounds but cannot eliminate contention, thermal changes, or compiler layout effects. Exploratory LLVM22/libomp22 timings are kept separate from later comparisons. Earlier manifests recorded `clang` plus libomp21.1.8 without independently resolving the compiler version. The final toolchain record identifies Homebrew Clang21.1.8/libomp21.1.8, and the final whole-command baseline is explicitly rebuilt with the matching compiler. Native final UBSan uses explicit Apple `/usr/bin/clang`; its compiler differs from the performance builds. Early exploratory deployment warnings do not establish macOS 11 compatibility.

The accepted final transport source is `af1442727a737eab81571074f7e6f10e366a9019ecedf70cc329fe661fd30c19`; the full final ARM plugin is `4472bf71a40362d9b0ee18d32ef80a70ceaaff14c0b7199c2e4b94fe72856bb4` (SHA-256). The earlier 12-case comparison retains its separate `c4c184…` source and `e2d01d…` whole-plugin provenance; its results are preserved separately. The final tested plugin was installed only after frozen source/ado equality and dependency checks; [installation record]({A}/install_final.json) preserves the hashes and previous-binary backup. Frozen baseline/candidate snapshots preserve concurrent ownership, cancellation, and empty-string alignment changes. Earlier snapshots and later acceptance are identified separately.

The [archive]({A}/ARCHIVE.md) contains raw CSVs, process summaries, exact manifests and commands, drivers, harnesses, profiles, validation logs, minimal transport source bundles, and the production baseline-to-candidate source patch. Large text logs are losslessly compressed; inventory records original-content and archived-byte SHA-256. Compiled binaries are omitted. Standalone transport builds can be reconstructed from the source bundles; full-plugin builds also need the recorded repository/build dependencies. The [campaign tool](../validation/benchmark_transport_campaign.py) prepares reproducible experiments and refuses incomplete summaries.
'''
(REPO/'docs/PERFORMANCE_TRANSPORT_ADAPTIVE.md').write_text(text)
print(len(text.splitlines()),'lines written')
