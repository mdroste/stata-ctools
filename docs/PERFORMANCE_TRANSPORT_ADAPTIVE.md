# Adaptive Stata–C transport

The final implementation substantially improves tiny transfers, wide numeric loads, and full-length `str2045`. Final mixed-string, narrow numeric, and large sorted-write controls are near parity. Some store phases still lose performance, and the evidence does **not** establish improvement for every workload.

**Completed:** exact-final-source transport and whole-command benchmarks, same-binary controls, full command acceptance, native UBSan checks, and installation of the tested plugin. The completed measurements below use the frozen production candidate; their provenance will not be replaced by follow-up results.

## Corrected implementation: measured transport

These measurements use the exact final K ≥ 128, selected N ≥ 50,000 policy and final frozen source. Times are milliseconds for load + store + cleanup, including destruction of the production pthread pool. Two matched process blocks each use five timed repetitions after a warmup, plus a separate fully verified pass: **336 records across four processes**. Speedup is the median paired baseline/candidate ratio; ranges are the two block ratios, not confidence intervals. Marginal median times therefore need not divide to the reported paired ratio.

| Workload | Baseline ms | Candidate ms | Total speedup | Block range | Load speedup | Store speedup |
|---|---:|---:|---:|---:|---:|---:|
| 100 × 20 numeric | 0.298 | 0.015 | 21.13× | 17.76–24.50× | 23.00× | 11.57× |
| 200K × 64 byte | 12.490 | 12.466 | 1.00× | 0.98–1.02× | 1.03× | 0.97× |
| 200K × 128 byte | 38.777 | 24.565 | 1.58× | 1.57–1.59× | 2.22× | 0.98× |
| 200K × 256 cycling numeric | 112.985 | 45.657 | 2.47× | 2.47–2.48× | 3.81× | 1.02× |
| 1M × 128 cycling numeric | 365.811 | 165.858 | 2.21× | 2.16–2.25× | 4.18× | 0.95× |
| 50K × 128 cycling numeric | 7.170 | 7.347 | 1.00× | 0.84–1.15× | 1.20× | 0.78× |
| 500K × 128, two-thirds selected | 140.889 | 66.564 | 2.12× | 2.10–2.13× | 3.16× | 1.19× |
| 200K × 4 full-length str2045 | 219.230 | 138.887 | 1.58× | 1.56–1.59× | 1.74× | 1.15× |
| 200K × 4 str2045, eight-byte values | 41.627 | 41.697 | 1.00× | 0.99–1.01× | 1.00× | 1.01× |
| 1M × 20 numeric control | 62.081 | 62.241 | 1.00× | 0.99–1.00× | 1.00× | 0.99× |
| 500K × 20 mixed str32 | 24.855 | 24.901 | 1.00× | 0.98–1.02× | 0.98× | 1.02× |
| 1M × 20 numeric, sorted writes | 74.127 | 74.093 | 1.00× | 1.00–1.01× | 1.00× | 1.00× |

The 1M × 128 and 200K × 256 numeric **load** phases improved 4.18× and 3.81×. Full-length str2045 loads improved 1.74×, and cleanup fell from roughly 24 ms to 0.09 ms by avoiding individual fallback allocations. Dense-filtered stores improved 1.19× after vectorizing map validation. The byte64 case now uses column loading; its old 0.75× same-binary tile result motivated the higher gate.

Store phases still cost more in two changed-path cases: 1M × 128 cycling numeric was 0.95× in both blocks, despite 2.21× total transport; 50K × 128 store was 0.78×, with uncertain total transport near parity (0.84–1.15×). The earlier three-block comparison had mixed-str32 and sorted-control total ratios of 0.85× and 0.92×; those losses did not repeat here. Both comparisons remain evidence of variability rather than permission to discard inconvenient results. [Complete final phase results](benchmarks/transport_adaptive_20260929/final_confirmed/summary.csv) and [process medians](benchmarks/transport_adaptive_20260929/final_confirmed/process.csv) include every case; the [earlier 64-column-gate comparison](benchmarks/transport_adaptive_20260929/final_comparison/summary.csv) is preserved as superseded evidence.

A same-binary follow-up (312 records, two process blocks) puts mixed str32 at 0.997×, large sorted writes at 0.998×, and numeric controls at 1.003× using block-level ratios. Mixed str8 remains slightly slower at **0.972× (0.958–0.986×)**. Shared runtime state reduces confounding, but the two compiled code bodies retain different placement and use the final shared helpers. [Block-level control results](benchmarks/transport_adaptive_20260929/acceptance/dual_impl_control/block_summary.csv) preserve this distinction.

## Whole-command confirmation

`csort id, verbose threads(8) nostream` was measured with one warmup and seven timed calls in each of two reversed process blocks. A monotonic clock includes ado/plugin dispatch and sorting; fixture loading and validation are outside timing. Each result is checked cell by cell and against native Stata's sorted datasignature. The four processes produced **128 records** including warmups.

| Whole csort workload | Baseline ms | Candidate ms | Paired speedup | Block range |
|---|---:|---:|---:|---:|
| 100 × 20 numeric | 1.050 | 0.298 | 3.54× | 3.31–3.78× |
| 100K × 20 numeric | 5.345 | 5.078 | 1.05× | 1.01–1.10× |
| 200K × 128 numeric | 71.576 | 44.027 | 1.63× | 1.59–1.66× |
| 100K × 4, ID plus three str2045 | 124.719 | 101.541 | 1.23× | 1.22–1.23× |

All four whole-command medians improved in both final blocks, including a modest 1.05× narrow numeric gain. The earlier narrow comparison ranged from 0.69× to 1.64× and remains in the archive as evidence of run-to-run variability. These four commands validate that transport gains reach an actual command; they are not estimates for all ctools commands. [Whole-command results](benchmarks/transport_adaptive_20260929/acceptance/commands_final/summary.csv) preserve both blocks and all phases.

## Selected implementation

- **Small transfers:** estimate work by numeric columns and declared string widths, then use serial load/store when pool setup would dominate. Unknown widths receive a conservative cost. A plain cell-count threshold was rejected because string-heavy cases regressed sharply.
- **Numeric loads:** 512-row tiles for all-numeric varlists with **K ≥ 128 and selected N ≥ 50,000**; no separate cell-count gate. A same-binary scheduler switch showed 64-column byte loads losing, while all five numeric storage types gained at 128/256 columns. Follow-up tests also found losses for compact byte data at two million cells and for very short, wide byte data above six million cells. The final row threshold avoids all measured short-byte losses; it gives up gains for some 64-column int/long/float/double and short-wide cycling datasets.
- **Wide strings:** reserve a larger packed arena for known str2045 values, retain actual-length density, and retry the smaller arena if reservation fails. Flat declared-width slots were rejected because short values in wide declarations lost performance.
- **Unused identity order:** 11 nonordering callers request `CTOOLS_LOAD_NO_SORT_ORDER`. The observed allocation falls from **4N bytes to zero**: 40 MB at 10M observations. A same-binary toggle supports modest numeric gains; strL is approximately unchanged and mixed-string timing is uncertain. This is an allocation saving, not a measured equal fall in peak RSS.
- **Numeric writes:** validate destination bounds and uniqueness before parallel writes; vectorize the common increasing-map check and retain serial semantics for duplicates. Direct sorted gathers apply only when the full source column fits 262,144 doubles (2 MiB); larger sources retain staging.

General store tiling, interleaved scheduling, unconditional OpenMP stores, unrestricted direct gathers, and unweighted small-transfer dispatch were rejected or restricted after regressions. The [detailed evidence](benchmarks/transport_adaptive_20260929/report_draft.md) records those results, including controls and rejected variants.

## Correctness and profiling

The full baseline, earlier candidate, and accepted final candidate suites each passed **4,161 checks with zero failures and 48 expected method-difference skips**. Four native transport suites passed with UBSan in serial and OpenMP builds: **eight successful runs**, including bounds, allocation failures, strL, sorted/duplicate maps, cancellation, and dispatch boundaries. The ARM plugin passed dependency checks with deployment target macOS 11 and compatible static libomp 21.1.8. These are host build/check results, not new cross-platform runtime measurements.

The initial ASan runtime deadlocked during initialization before `main`; this is not a passing sanitizer run. An independent Homebrew21.1.8 ASan probe also timed out before `main`; neither run tested transport code. Added Linux CI coverage is not a claimed CI pass. [Native results and startup diagnosis](benchmarks/transport_adaptive_20260929/acceptance/final_native_checks_final/ubsan_results.json), [baseline acceptance](benchmarks/transport_adaptive_20260929/acceptance/production_baseline/acceptance.log), [candidate acceptance](benchmarks/transport_adaptive_20260929/acceptance/production_final/acceptance.log), and [dependency metadata](benchmarks/transport_adaptive_20260929/build_checks_final.json) preserve the evidence.

Separate five-second CPU samples confirm the transition from column-worker loads to tiled OpenMP loads through the same Stata callback. Stores still use the checked callback. Sampling counts pool running and sleeping threads; they are not CPU percentages or hardware cache-miss measurements and do not prove a spin-wait explanation. Profiled timings are excluded from the benchmark summaries.

## Scope, uncertainty, and reproduction

Eight exploratory component campaigns cover N=1–10M and K=1–2,000, all five numeric storage types, mixed fixed strings, bounded textual strL, filtered selections, sorted writes, subsets of a wider host dataset, and 1/8/12 requested threads. Their **9,506 records in 48 processes** include warmups and verification. Subsequent independent integration comparisons and same-binary toggles test the selected changes. This is a targeted matrix, not all possible configurations or a representative workload average.

All Stata launches used `oldstata` through interactive login zsh. Transport uses `CLOCK_MONOTONIC`; it excludes data generation, verification, and ado overhead. The host was Apple M4 Pro, 48 GiB RAM, macOS 26.7, StataNow/MP 18.5. Another user-authorized Stata job ran concurrently. Sequential execution, alternating order, process-level medians, and unchanged controls reduce some confounds but cannot eliminate contention, thermal changes, or compiler layout effects. Exploratory LLVM22/libomp22 timings are kept separate from later comparisons. Earlier manifests recorded `clang` plus libomp21.1.8 without independently resolving the compiler version. The final toolchain record identifies Homebrew Clang21.1.8/libomp21.1.8, and the final whole-command baseline is explicitly rebuilt with the matching compiler. Native final UBSan uses explicit Apple `/usr/bin/clang`; its compiler differs from the performance builds. Early exploratory deployment warnings do not establish macOS 11 compatibility.

The accepted final transport source is `af1442727a737eab81571074f7e6f10e366a9019ecedf70cc329fe661fd30c19`; the full final ARM plugin is `4472bf71a40362d9b0ee18d32ef80a70ceaaff14c0b7199c2e4b94fe72856bb4` (SHA-256). The earlier 12-case comparison retains its separate `c4c184…` source and `e2d01d…` whole-plugin provenance; its results are preserved separately. The final tested plugin was installed only after frozen source/ado equality and dependency checks; [installation record](benchmarks/transport_adaptive_20260929/install_final.json) preserves the hashes and previous-binary backup. Frozen baseline/candidate snapshots preserve concurrent ownership, cancellation, and empty-string alignment changes. Earlier snapshots and later acceptance are identified separately.

The [archive](benchmarks/transport_adaptive_20260929/ARCHIVE.md) contains raw CSVs, process summaries, exact manifests and commands, drivers, harnesses, profiles, validation logs, minimal transport source bundles, and the production baseline-to-candidate source patch. Large text logs are losslessly compressed; inventory records original-content and archived-byte SHA-256. Compiled binaries are omitted. Standalone transport builds can be reconstructed from the source bundles; full-plugin builds also need the recorded repository/build dependencies. The [campaign tool](../validation/benchmark_transport_campaign.py) prepares reproducible experiments and refuses incomplete summaries.
