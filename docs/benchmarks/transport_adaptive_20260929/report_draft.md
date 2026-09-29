# Transport optimization: component evidence and final verification

**Evidence complete. Component and integration timings, allocation and scheduler toggles, CPU profiles, whole-command timings, both full command suites, and eight native UBSan configurations are complete. The 128-column gate removes the byte64 regression. Short-byte tests reject both the two-million-cell policy and a six-million-cell alternative. The final policy is K ≥ 128 and selected N ≥ 50,000, without a cell-count check; exact-final-source correctness and transport timings have passed; final whole-command timings and same-binary controls are complete.**

The clearest opportunities are avoiding worker-pool setup for genuinely small transfers and improving row locality when loading many numeric columns. Maximum-width fixed strings also benefit from avoiding repeated fallback allocations. Several plausible write optimizations lost performance on other configurations and were rejected or restricted.

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

Store phases still cost more in two changed-path cases: 1M × 128 cycling numeric was 0.95× in both blocks, despite 2.21× total transport; 50K × 128 store was 0.78×, with uncertain total transport near parity (0.84–1.15×). The earlier three-block comparison had mixed-str32 and sorted-control total ratios of 0.85× and 0.92×; those losses did not repeat here. Both comparisons remain evidence of variability rather than permission to discard inconvenient results. [Complete final phase results](final_confirmed/summary.csv) and [process medians](final_confirmed/process.csv) include every case; the [earlier 64-column-gate comparison](final_comparison/summary.csv) is preserved as superseded evidence.

## Earlier 128-column / two-million-cell integration comparison

Times are milliseconds for load + store + cleanup, including destruction of the production pthread pool. Three matched process blocks each use seven timed repetitions after a warmup, plus a separate fully verified pass: **648 records across six processes**. Speedup is the median paired baseline/candidate ratio; ranges are the three block ratios, not confidence intervals. Marginal median times therefore need not divide to the reported paired ratio.

| Workload | Baseline ms | Candidate ms | Paired speedup | Block range |
|---|---:|---:|---:|---:|
| 100 × 20 numeric | 0.260 | 0.034 | 18.57× | 6.27–72.88× |
| 200K × 64 byte | 14.387 | 12.509 | 1.07× | 1.06–1.20× |
| 200K × 128 byte | 45.988 | 25.121 | 1.74× | 1.62–1.88× |
| 200K × 256 cycling numeric | 129.729 | 52.321 | 2.44× | 2.03–2.49× |
| 1M × 128 cycling numeric | 391.977 | 175.262 | 2.17× | 1.73–2.27× |
| 50K × 128 cycling numeric | 7.333 | 6.580 | 1.08× | 0.91–1.22× |
| 500K × 128, two-thirds selected | 167.430 | 73.084 | 2.03× | 2.00–2.30× |
| 200K × 4 full-length str2045 | 235.058 | 143.358 | 1.52× | 1.47–1.74× |
| 200K × 4 str2045, eight-byte values | 41.867 | 42.684 | 1.00× | 0.98–1.00× |
| 1M × 20 numeric control | 64.185 | 64.893 | 1.00× | 0.99–1.01× |
| 500K × 20 mixed str32 | 25.018 | 29.127 | 0.85× | 0.79–1.03× |
| 1M × 20 numeric, sorted writes | 77.696 | 84.300 | 0.92× | 0.87–1.04× |

The 1M × 128 and 200K × 256 numeric **load** phases improved 4.10× and 3.97×. Full-length str2045 loads improved 1.81×, and cleanup fell from roughly 25 ms to 0.1 ms by avoiding individual fallback allocations. Dense-filtered stores improved 1.30× after vectorizing map validation. The byte64 case now uses column loading; its old 0.75× same-binary tile result motivated the higher gate.

The mixed-str32 total ratio was 0.85×, with two of three blocks below parity; large sorted-write total was 0.92×, also with two blocks below parity. The numeric control was 1.00×. These differences and the wide tiny-case range must remain visible. [Complete final phase results](final_revision/summary.csv) and [process medians](final_revision/process.csv) include every case; the [earlier 64-column-gate comparison](final_comparison/summary.csv) is preserved as superseded evidence.

## Same-binary control comparison

The follow-up plugin links baseline and final data-I/O bodies as separate, namespaced translation units against the same final types, arena, thread pool, and OpenMP runtime. Implementation pointers are selected before timing; each call destroys the pool. Two process blocks alternate AB/BA pairs and reverse case order, with 11 timed pairs, a warmup pair, and a fully verified pair per case: **312 records**. Sorted writes are restored and their datasignature checked before the next implementation sees the fixture.

| Control, total including cleanup | Median block speedup | Two-block range |
|---|---:|---:|
| 500K × 20 mixed str8 | 0.972× | 0.958–0.986× |
| 500K × 20 mixed str32 | 0.997× | 0.988–1.005× |
| 500K × 20 mixed str244 | 1.000× | 0.994–1.005× |
| 1M × 20 numeric | 1.003× | 1.000–1.007× |
| 1M × 20 numeric, sorted writes | 0.998× | 0.992–1.005× |
| 200K × 64 byte | 1.000× | 0.995–1.005× |

The earlier mixed-str32 and large sorted-write losses do not repeat in this control. Mixed str8 remains slightly slower in both blocks, about 0.972× at their median; it is not silently treated as parity. This design reduces process/runtime differences but is not an identical-body flag toggle: the code bodies still have different placement, and both use the final shared helpers. It isolates data-I/O-body behavior rather than every historical source difference.

The reported ratios are computed from the two `block_speedups` entries, not from the summary's 22 individual-pair extrema. [Raw control records](acceptance/dual_impl_control/raw.csv), [original summary](acceptance/dual_impl_control/summary.csv), and [block-level summary](acceptance/dual_impl_control/block_summary.csv) preserve both scales without treating 22 calls as 22 independent processes.

## Whole-command confirmation

`csort id, verbose threads(8) nostream` was measured with one warmup and seven timed calls in each of two reversed process blocks. A monotonic clock includes ado/plugin dispatch and sorting; fixture loading and validation are outside timing. Each result is checked cell by cell and against native Stata's sorted datasignature. The four processes produced **128 records** including warmups.

| Whole csort workload | Baseline ms | Candidate ms | Paired speedup | Block range |
|---|---:|---:|---:|---:|
| 100 × 20 numeric | 1.050 | 0.298 | 3.54× | 3.31–3.78× |
| 100K × 20 numeric | 5.345 | 5.078 | 1.05× | 1.01–1.10× |
| 200K × 128 numeric | 71.576 | 44.027 | 1.63× | 1.59–1.66× |
| 100K × 4, ID plus three str2045 | 124.719 | 101.541 | 1.23× | 1.22–1.23× |

All four whole-command medians improved in both final blocks, including a modest 1.05× narrow numeric gain. The earlier narrow comparison ranged from 0.69× to 1.64× and remains in the archive as evidence of run-to-run variability. These four commands validate that transport gains reach an actual command; they are not estimates for all ctools commands. [Whole-command results](acceptance/commands_final/summary.csv) preserve both blocks and all phases.

## Earlier whole-command comparison

`csort id, verbose threads(8) nostream` was measured with one warmup and seven timed calls in each of two reversed process blocks. A monotonic clock includes ado/plugin dispatch and sorting; fixture loading and validation are outside timing. Each result is checked cell by cell and against native Stata's sorted datasignature. The four processes produced **128 records** including warmups.

| Whole csort workload | Baseline ms | Candidate ms | Paired speedup | Block range |
|---|---:|---:|---:|---:|
| 100 × 20 numeric | 0.669 | 0.307 | 2.25× | 1.65–2.86× |
| 100K × 20 numeric | 8.547 | 8.006 | 1.16× | 0.69–1.64× |
| 200K × 128 numeric | 85.568 | 47.413 | 1.82× | 1.63–2.00× |
| 100K × 4, ID plus three str2045 | 143.154 | 112.566 | 1.27× | 1.19–1.36× |

The wide numeric and full-width-string command improvements occur in both blocks. The narrow numeric result is uncertain, with one block slower. These four commands validate that transport gains reach an actual command; they are not estimates for all ctools commands. [Whole-command results](acceptance/commands/summary.csv) preserve both blocks and all phases.

## Measurement protocol and provenance

The eight completed campaigns contain **9,506 transfer records**, including warmups and verification, from **48 separate Stata processes**. They produce 1,058 process-by-case medians and 340 baseline–variant case comparisons. Every recorded process completed with `TRANSPORT_COMPLETE RC=0`; the summarizer also requires every expected timing, cleanup, and verification record. [Campaign index](campaign_index.csv), [raw records](campaign_raw.csv), [phase summaries](campaign_phase_summary.csv), and [resource records](campaign_resources.csv) contain only these eight campaigns. Prepared but unused experiments and the subsequent integration comparisons are excluded.

Each process loads one benchmark plugin. Variants run in forward order in the first process block and reverse order in the second. Every case has one warmup, 3–11 timed repetitions, and a separate verification pass. The verification compares every loaded and written value, checks the selected-observation map against Stata, and asserts the restored dataset's datasignature. Numeric comparisons are bitwise, including ordinary and extended missing values. The production lifecycle is reproduced by destroying the global pthread pool after each transfer; the older benchmark retained that pool between calls and understated small-transfer overhead.

`CLOCK_MONOTONIC` times loading, writing, and cleanup separately. Load time includes selection construction and allocation; cleanup includes freeing the loaded data and destroying the pool. “Transport” is their sum. It excludes data generation, plugin dispatch overhead, permutation construction, restoration, verification, and Stata ado work. These are in-memory SPI transport measurements, not disk-loading or whole-command timings. All Stata executions use the machine's `oldstata` wrapper; monotonic timing avoids the system-clock changes that wrapper performs.

The tables first take the median of timed calls **within each process**, then summarize the process blocks: two for component campaigns and four for the final transport comparison. The speedup is the median of matched process-block baseline/candidate ratios; 2× means half the measured time. The reported range is the minimum and maximum of the matched block ratios. It is not a confidence interval, and the repeated calls are not treated as independent process replications. The ratio of the displayed marginal medians need not equal the median paired ratio.

Measurements were made on the Apple M4 Pro host with 12 logical cores and 48 GiB RAM, running macOS 26.7 and StataNow/MP 18.5. Another Stata job ran contemporaneously; the user authorized continuing under that contention. Scheduling, thermal state, background work, and retained runtime state therefore remain sources of uncertainty. The experiments themselves ran sequentially. No global OpenMP wait-policy change was adopted.

The component campaigns used Homebrew LLVM Clang 22 and libomp 22. Early campaign link commands specified a macOS 11 deployment target although the installed libomp was built for macOS 26; those builds emitted deployment-target warnings. Later component commands matched the runtime's target. These measurements establish neither macOS 11 compatibility nor behavior on other platforms. The frozen final comparison instead records compiler command **`clang` and the compatible libomp 21.1.8 build** for both variants. Those earlier manifests did not independently record the resolved compiler version, so they must not be described as system-Clang builds. The later final toolchain check resolves `clang` to Homebrew Clang 21.1.8; the final whole-command baseline is rebuilt explicitly with that version to match the final candidate. Results from the two toolchains must not be pooled into a single speedup estimate.

The component baseline is the starting working tree, not Git HEAD. Its `ctools_data_io.c` SHA-256 is `1e42ccdc0242af59a0eb3762b8016de8fa585b7620ddb781e1e671c5c602e188`. Every campaign manifest records the full source-file hashes, harness hash, plugin hash, compiler command, case definitions, and process order. The final comparison uses fresh `integration_baseline` and `integration_candidate` snapshots, preserving concurrent ownership/type and command edits shared by the comparison. At preparation, its candidate transport SHA-256 was `db0adcfde3f0bd484546b1af4957ce4f0b3a556a9106d40c0d1bfae2e6ab0f22`; its shared `ctools_types.c` SHA-256 was `237fe172aa03503315685081c227e461ef2a52dfd36a3199668f71431e894bf7`. The [final comparison manifest](final_comparison/manifest.json) identifies the frozen benchmark build. The subsequently merged validation snapshot includes concurrent cancellation-return-code mapping and empty-string alignment changes, the corrected 128-column gate, and the vectorized observation-map check. Its `ctools_data_io.c` SHA-256 is `c4c18472878f8b89767c4f15bc24448c8798d63f5f504a6197758e4e8381f1f4`. This is distinct from the earlier frozen comparison; its timing and acceptance results are recorded separately, without rewriting that provenance or attributing unmeasured changes to the earlier timings.

## Frozen transport comparison with the provisional 64-column gate

The fresh integration baseline and provisional 64-column-gate candidate completed **29 cases in eight processes**, with seven timed repetitions per case and four matched process blocks: **2,088 records**, 232 process-by-case medians, and 116 paired phase comparisons. All eight processes passed cell verification and datasignature checks. These measurements record compiler command `clang` and compatible libomp 21.1.8; the compiler version was not independently pinned in that earlier manifest. [Raw final records](final_comparison/raw.csv), [complete final summaries](final_comparison/summary.csv), and [process-level medians](final_comparison/process.csv) retain all cases. Allocation opt-out was not requested in this transport matrix; its separate controlled test appears below.

The table reports total measured transport, including cleanup, in milliseconds. It includes gains, regressions, and uncertain controls. It is not an average over a representative distribution of user workloads.

| Workload | Baseline ms | Candidate ms | Paired speedup | Four-block range |
|---|---:|---:|---:|---:|
| 100 × 20 numeric | 0.276 | 0.017 | 16.47× | 14.94–19.79× |
| 1K × 20 numeric | 0.345 | 0.118 | 2.95× | 1.72–3.43× |
| 200K × 64 byte | 14.073 | 16.500 | 0.84× | 0.79–1.18× |
| 200K × 128 double | 71.014 | 31.280 | 2.30× | 1.81–2.93× |
| 200K × 256 cycle | 124.106 | 54.696 | 2.23× | 2.12–2.60× |
| 1M × 128 cycle | 422.477 | 195.966 | 2.18× | 1.74–2.58× |
| 50K × 128 cycle | 7.862 | 6.874 | 1.17× | 0.80–1.55× |
| 10K × 2,000 cycle | 27.296 | 25.923 | 1.05× | 0.86–1.23× |
| 200K × 128 cycle, sorted store | 71.019 | 35.184 | 2.02× | 1.75–2.35× |
| 500K × 20 mixed str32 | 31.892 | 33.880 | 0.91× | 0.76–1.37× |
| 200K × 4 full-length str2045 | 234.351 | 150.801 | 1.55× | 1.50–1.72× |
| 200K × 4 str2045 holding eight-byte values | 42.666 | 43.248 | 0.98× | 0.96–1.06× |
| 500K × 128 cycle, two-thirds selected | 154.921 | 104.136 | 1.45× | 1.35–1.75× |
| 1M × 128 cycle, 1% selected | 24.521 | 23.163 | 1.04× | 1.02–1.08× |
| 1M × 20 numeric control | 67.039 | 68.561 | 0.94× | 0.85–1.07× |
| 100K × 4 str244, no width hints | 34.129 | 33.479 | 1.00× | 0.95–1.04× |

The wide double and mixed-numeric load gains survive the toolchain change: 200K × 128 double loads improved 3.66×, 200K × 256 cycling numeric loads 3.37×, and 1M × 128 cycling loads 3.95×. The 200K × 128 cycling workload also improved at one requested thread (3.95× load; 2.52× transport) and 12 threads (2.93× load; 2.09× transport), supporting a locality benefit beyond additional parallelism. Full-length str2045 loads improved 1.81×, with a tight 1.79–1.86× block range.

**The 64-column byte case regressed.** Its load median increased from 6.594 to 9.106 ms, a paired speedup of 0.714×; total transport was 0.839×. Three of the four load ratios were below one. This reverses the exploratory LLVM22 result. The same-binary test below reproduces the loss, and the corrected implementation raises the minimum to 128 columns. This provisional result is retained to document why the policy changed; it is not a claim about the corrected candidate.

**Dense filtered writes also regressed.** In the 500K × 128 case selecting two-thirds of observations, load improved 2.60× but store fell to 0.744×, with all four store ratios below one (0.668–0.998×). Total transport still improved 1.45×. The vectorized map-validation follow-up below addresses this measured cost; the subsequent 12-case comparison confirms a 1.302× store gain against the original baseline (1.103–1.495× across three blocks). The provisional table is not evidence that every phase improved.

The mixed-str32 control had 0.91× median paired total speedup, with a wide 0.76–1.37× range; the 1M × 20 numeric control was 0.94×, range 0.85–1.07×. Large numeric sorted-store control total was 0.94×, and its store phase was 0.89×. Unknown-width strings, short values in str2045, and checked loads were close to parity. These controls and their variability preclude a universal no-regression claim. Neither low-N width extremes nor the source-array gather cutoff should be declared fully optimized on these measurements alone.

## Completed component results

Selected observations are shown below; milliseconds refer to the median of process medians. “Cycle” rotates through byte, int, long, float, and double storage. Mixed columns alternate numeric and fixed strings. These rows document specific component decisions, not an aggregate package speedup; their baseline and toolchain differ from the frozen final comparison above. [Selected source rows](component_selected.csv) preserve the exact numbers.

| Component and workload | Phase | Baseline ms | Variant ms | Paired speedup | Block-ratio range |
|---|---|---:|---:|---:|---:|
| Adaptive small-transfer dispatch, 1 × 20 numeric | Transport | 0.2225 | 0.0065 | 36.89× | 25.37–48.40× |
| Adaptive small-transfer dispatch, 100 × 20 mixed str32 | Transport | 0.2995 | 0.0370 | 8.09× | 7.16–9.03× |
| Numeric row tiles, 200K × 64 cycle | Load | 25.337 | 9.002 | 2.82× | 2.79–2.84× |
| Numeric row tiles, 200K × 128 double | Load | 57.685 | 16.237 | 3.64× | 3.16–4.12× |
| Numeric row tiles, 200K × 256 double | Load | 121.505 | 27.971 | 4.39× | 4.08–4.70× |
| Numeric row tiles, 1M × 128 cycle | Load | 290.431 | 71.047 | 4.09× | 4.08–4.10× |
| Numeric cell-count gate, 50K × 128 cycle | Load | 10.126 | 3.918 | 2.66× | 2.31–3.01× |
| Packed reservation, 200K × 4 str2045 with full-length values | Load | 105.799 | 71.814 | 1.48× | 1.44–1.52× |
| Direct numeric sorted gather, 200K × 64 cycle | Store | 11.432 | 9.511 | 1.20× | 1.14–1.27× |
| Omit identity-order allocation, 10M × 20 cycle | Load | 108.936 | 117.229 | 0.93× | 0.91–0.94× |

The adaptive small-transfer policy estimates work using both column count and string width. The selected implementation assigns one unit to a numeric column and `2 + ceil(width/8)` units to a string column, treating unknown widths conservatively. It compares total work with `16,384 + 256*K`, using overflow-safe arithmetic. This retains the large absolute savings from avoiding pool setup for tiny datasets while keeping more expensive string transfers parallel. Large ratios at N=1 represent fractions of a millisecond; they are not whole-command speedups. The adaptive screen still contains slower controls and uncertain boundary results, and the subsequent integrated timings and native boundary checks are reported separately.

Numeric tiling changes traversal order to reuse Stata rows across columns. The frozen candidate gate requires an all-numeric varlist with at least 64 variables, 4,096 selected rows, and two million selected cells. The corrected implementation raises the variable-count boundary to 128 after the same-binary follow-up confirms the byte-data regression at 64. The gate expands coverage beyond the former 200K-row cutoff: 50K × 128 improved materially, whereas 10K × 200 gave a small and uncertain gain. The 64-, 128-, and 512-row tile sweep did not establish one winner for every storage layout; 512 rows was retained. The selected gate is based on the passed varlist, not the entire host dataset, so a small subset of a very wide host remains a separate case.

For known `str2045` columns, the packed-arena candidate reserves up to 2 GiB per column but packs actual string bytes densely. This reduces per-string fallback allocations without touching every declared-width slot. The selected revision retries the original 64-byte-per-row arena estimate if the larger reservation fails, then preserves the existing individually allocated fallback. Short values in wide declared columns did not show a general speed gain. Reservation size is not resident-memory consumption; other platforms can account for committed memory differently.

Direct numeric sorted stores are limited to source columns of at most **262,144 rows**, or 2 MiB of doubles. The check uses the full source-column length, not the current worker chunk. The component screen supports a gain at 200K rows but showed a regression at 1M × 20 numeric, so larger random gathers retain their bounded staging buffer. The exact optimal crossover is not established. Native boundary checks passed; the large sorted-write control remains uncertain and receives a further same-binary comparison.

The identity-order opt-out eliminates an otherwise unused `N * sizeof(perm_idx_t)` allocation. Instrumentation observed **4,000,000 bytes versus zero at 1M rows**, and **40,000,000 bytes versus zero at 10M rows**, independently of K. This is an allocation reduction, not a measured equal reduction in peak RSS. Timing was mixed: the separate-plugin screen was about 0.82× for 10M × 1 numeric and 0.93× for 10M × 20 numeric. The subsequent same-binary runtime-flag experiment below removes cross-plugin code-layout differences and resolves the observed numeric timing concern. The 11 callers have now adopted the opt-out; their full command-level correctness suites passed for both baseline and candidate.

The numeric observation-map change is a correctness safeguard: validate all destinations and prove uniqueness before parallel writes. Increasing maps use the existing bounds pass; decreasing or unordered maps require a valid uniqueness proof, and duplicates or unavailable scratch space retain serial semantics. A subsequent vectorized increasing-map check preserves the safeguard while reducing its common-path cost, as measured below.

## Observation-map validation: recovering dense-write performance

The common increasing-map check was rewritten to let the compiler vectorize validation of bounds and strict ordering. The fallback still checks unusual map orders and preserves serial behavior for duplicates. This targets the additional store cost observed in the provisional integration result, while preserving the safety conditions for parallel writes.

The focused scalar-versus-vector comparison completed four cases, three matched process blocks, and nine timed repetitions: **264 records in six processes**. It compares two versions of the candidate, not the original baseline. Dense cases select two-thirds of 500K observations; the sparse control selects 1% of 1M observations. [Raw map records](map_comparison/raw.csv) and [phase summaries](map_comparison/summary.csv) retain every result.

| Selected numeric map | Scalar store ms | Vector store ms | Paired speedup | Three-block range |
|---|---:|---:|---:|---:|
| Dense, K=1 | 1.949 | 1.317 | 1.366× | 1.216–1.698× |
| Dense, K=20 | 10.141 | 7.735 | 1.461× | 1.311–1.475× |
| Dense, K=128 | 63.285 | 51.414 | 1.235× | 1.218–1.317× |
| Sparse, K=128 | 8.785 | 8.837 | 1.006× | 0.977–1.022× |

All three dense configurations improved in every process block; loads were approximately unchanged. The sparse store control remained near parity. These comparisons support the vectorized validation implementation but do not by themselves establish recovery to the original baseline; that is a separate question for the final merged-source comparison.

## Numeric scheduler: controlled policy correction

A second same-binary experiment changes only a runtime switch for numeric row tiling; both branches otherwise execute the same compiled implementation and retain identity-order allocation. It covers 200K × 64/128/256 columns for each homogeneous byte, int, long, float, and double storage type, plus 50K × 128 and 10K × 2,000 cycling types. Two reversed process blocks contain 11 timed pairs per case, with separate warmup and fully verified pairs: **884 records**. [Raw records](numeric_scheduler_toggle/raw.csv) and [summaries](numeric_scheduler_toggle/summary.csv) preserve all cases.

| Load comparison: tiles versus columns | Paired speedup | Matched-block range or scope |
|---|---:|---|
| 200K × 64 byte | 0.754× | 0.729–0.779× |
| 200K × 64 int/long/float/double | 2.17–3.50× | Range of medians across four storage types |
| 200K × 128, all five storage types | 3.21–3.42× | Range of medians across five storage types |
| 200K × 256, all five storage types | 3.02–3.88× | Range of medians across five storage types |
| 50K × 128 cycling types | 1.638× | 1.583–1.693× |
| 10K × 2,000 cycling types | 1.137× | 1.046–1.229× |

The byte-column loss at 64 variables is real in this controlled comparison; it is not just a different-binary or process-order artifact. At this stage the candidate policy used **K ≥ 128, selected N ≥ 4,096, and selected N*K ≥ 2,000,000**, retaining 512-row load tiles; the later short-byte follow-up below supersedes its row/cell conditions. This deliberately gives up gains available for some 64-column int/long/float/double datasets so the generic dispatch avoids the measured byte regression without adding metadata-query overhead or undocumented storage-type assumptions. The later short-byte follow-up gives up some short-wide gains to avoid measured byte-data regressions. The 12-case merged comparison confirms the 128-column boundary, while the earlier 29-case table remains explicitly labeled as the superseded 64-column-gate comparison.

## Compact short-byte boundary follow-up

A six-case same-binary follow-up reused the earlier scheduler-toggle plugin unchanged; it is not a new production build. It contains 312 records over two reversed process blocks, 11 timed pairs per case, and separately verified pairs. The two-million-cell gate is too permissive for compact homogeneous-byte data: at 15,625 × 128 (exactly two million cells), tiles achieved only **0.843× load** and **0.892× load-plus-cleanup**; at 4,096 × 512, the ratios were **0.738×** and **0.780×**, with both blocks below parity. [Raw records](numeric_small_byte_toggle/raw.csv) and [summaries](numeric_small_byte_toggle/summary.csv) preserve these losses.

The larger controls improved: 50K × 128 byte achieved 1.206× load-plus-cleanup, 100K × 128 byte 1.456×, 50K × 256 byte 1.398×, and 50K × 128 double 2.215×. A second follow-up tested 46,875 × 128 byte, 4,096 × 1,536 byte, and 10K × 2,000 byte: 156 records with the same binary and two reversed blocks. Load-plus-cleanup ratios were 1.053×, 0.959×, and 0.928× respectively. Thus a six-million-cell minimum still admits measured short-wide losses. [Six-million-cell records](numeric_sixm_toggle/raw.csv) and [summaries](numeric_sixm_toggle/summary.csv) retain the evidence.

The final policy is simpler and more conservative: **all numeric, K ≥ 128 and selected N ≥ 50,000**, with 512-row tiles and no separate cell-count gate. It retains the dispatch choices for all 12 completed merged-source comparison cases and all four whole-command cases, while avoiding every measured short-byte loss. It deliberately gives up the modest gain for 10K × 2,000 cycling data and any other beneficial short-wide case. The separate `production_final` snapshot has transport source SHA-256 `af1442727a737eab81571074f7e6f10e366a9019ecedf70cc329fe661fd30c19` and full ARM plugin SHA-256 `4472bf71a40362d9b0ee18d32ef80a70ceaaff14c0b7199c2e4b94fe72856bb4`. Its full suite passed 4,161 checks with zero failures and 48 expected method-difference skips; all eight serial/OpenMP native UBSan configurations passed using explicit Apple Clang21. The earlier frozen source and timings remain unchanged. The exact-final 12-case transport comparison is complete (336 records); its table is first in this appendix. Final whole-command timings are complete and reported below.

## Identity-order allocation: controlled follow-up

The separate final allocation matrix completed four cases in four processes, with five timed repetitions: **112 records**. Instrumentation confirmed baseline identity-order allocations of 4,000,000 bytes for the 1M-row cases, 800,000 bytes for the 200K-row fixed-string case, and 400,000 bytes for the 100K-row strL case; the candidate reported zero in every case. [Allocation records](final_allocation/raw.csv) preserve the measurements. Peak RSS remained approximately 2.54 GB in both variants, so these runs establish the removed allocation but not a corresponding peak-RSS reduction.

Those separate binaries also contain the other transport changes: the 1M × 128 total improved 2.53×, short str2045 values were unchanged, the 1M × 20 numeric control improved 1.06×, and the strL control was 0.87× with broad variation. These are combined-implementation comparisons and cannot identify the effect of the allocation flag alone.

The decisive follow-up holds the **binary, source, fixture, and runtime fixed**, toggling only `NO_SORT_ORDER` before the timer. Calls alternate flag order within each repetition; the second process block reverses that order. Eight cases each have 15 timed pairs, a warmup pair, and a separate fully verified pair: **544 records** across two processes. Both flags passed every value check and allocation assertion. The following times include loading and cleanup; they exclude the negligible empty read-only store timer. [Runtime-toggle raw records](no_order_toggle/raw.csv) and [summaries](no_order_toggle/summary.csv) provide the full results.

| Workload | Default ms | No order ms | Paired speedup | Two-block range |
|---|---:|---:|---:|---:|
| 1M × 1 numeric | 2.889 | 2.535 | 1.211× | 1.104–1.318× |
| 1M × 2 numeric | 2.204 | 2.017 | 1.089× | 1.068–1.111× |
| 1M × 20 numeric | 12.169 | 12.017 | 1.012× | 0.998–1.026× |
| 10M × 1 numeric | 8.628 | 7.654 | 1.127× | 1.109–1.146× |
| 10M × 2 numeric | 16.014 | 15.472 | 1.035× | 1.032–1.038× |
| 10M × 20 numeric | 121.413 | 119.075 | 1.019× | 1.013–1.026× |
| 1M × 20 mixed str8 | 67.257 | 68.535 | 0.979× | 0.915–1.042× |
| 1M × 2 strL8 | 170.959 | 170.620 | 1.002× | 1.002–1.002× |

The earlier numeric regressions do not persist in this controlled test. The 10M-row numeric cases improve in both process blocks, while the mixed-string case is slightly slower at the median with a range crossing parity. strL is effectively unchanged. The opt-out has been adopted by 11 callers that do not consume the identity order; it does not establish a runtime gain for every string-heavy call. Numeric saves are modest in absolute time—about 0.5–2.3 ms at 10M rows—alongside the directly observed 40,000,000-byte allocation saving.

## Rejected or restricted approaches

| Approach | Evidence and decision |
|---|---|
| Serial dispatch for every transfer below 65,536 cells | Rejected. It ignored string cost: 100 × 200 str244 fell to 0.33× and 1K × 20 str244 to 0.25×. Replaced by the width- and column-aware work estimate. |
| Flat slots for every `str2045` value | Rejected as a general policy. It helped full-length strings but a single column containing eight-byte values fell to 0.38×. The packed approach retains actual-length density. |
| General numeric store tiling | Rejected. At 200K rows the double-column cases improved about 1.33–1.45×, while byte-column cases were around 0.86–0.96×. |
| Interleaving column tasks by row chunk | Rejected as a general policy. Some wide cases gained, but the 200K × 8 numeric identity and sorted stores both fell to about 0.83×. |
| Moving all stores to OpenMP | Rejected. With the same numeric tiled loader, total transport for 200K × 128 cycle increased from 31.93 to 48.72 ms; the double case increased from 37.99 to 55.77 ms. Avoiding a second scheduler did not reliably offset its cost. |
| Unrestricted direct numeric sorted gathers | Restricted to full source columns no larger than 2 MiB; the 1M × 20 numeric case was 0.95× in the direct-gather screen. |
| Claiming a universal runtime gain from omitting the identity order | The same-binary toggle resolves the numeric regressions and supports modest numeric gains, but mixed strings remain uncertain and strL is essentially unchanged. Retain the direct allocation evidence. |

## CPU profiles

Separate five-second `sample` captures, requested at one-millisecond intervals, observed repeated 200K × 128 cycling-numeric identity transfers with eight requested threads. Both profile workloads completed with return code zero. Their timing records are excluded from all benchmark summaries. [Baseline profile](profile_final_baseline/cpu_profile.txt.gz), [candidate profile](profile_final_candidate/cpu_profile.txt.gz), and [extracted reported leaf counts](profile_leaf_counts.csv) retain the evidence.

The baseline call trees show `persistent_worker_thread → load_variable_thread →` Stata read callbacks; the candidate instead shows `load_tiled_variables.omp_outlined →` the same Stata read entry offset. Stores remain under `persistent_worker_thread → store_variable_thread → store_single_variable →` Stata's checked store path. The top-of-stack listing contains 1,140 samples in baseline `load_variable_thread`, 1,085 in candidate `load_tiled_variables.omp_outlined`, and 448 versus 983 in `store_single_variable`. The Stata read-path instruction at module offset `0x421934` appears 3,712 versus 1,911 times; the store-path instruction at `0x424818` appears 1,148 versus 2,807 times. Those paths are identified from their callers; Stata's internal functions are not symbolicated.

This confirms that the intended load traversal changed and is consistent with stores accounting for more of the observed work after faster loads. It does not establish a per-transfer write slowdown: the faster process can complete more transfers within the same sampling window. The profiles are dominated numerically by sleeping-thread leaves (`__psynch_cvwait`: 70,839 versus 62,359), and OpenMP barrier descendants include suspended condition waits. Counts pool observations from multiple threads, include sleeping threads, and are affected by sampling overhead. They are **not CPU percentages, wall-time shares, or hardware cache-miss measurements**. In particular, these captures do not establish that OpenMP spin waiting is the cause of a measured regression. The runtime scheduling decisions rely on the timed experiments, which rejected the blanket OpenMP-store replacement.

## Coverage and limits

The component matrix spans N=1 through 10M, K=1 through 2,000, all five numeric storage types through cycling, homogeneous byte and double data, fixed declarations of 8/32/244/2045 bytes, actual string lengths from zero to 2045, bounded textual strL, 1/8/12 requested threads, identity and sorted stores, and 1% and two-thirds filtered selections. It includes an eight-variable selection spread through a 256-variable host dataset. These are selected boundary and mechanism tests, not a complete Cartesian product. All component cases supplied width hints. The final comparison adds no-hint fixed strings, empty selections, and all-selected checked loads; additional map-order and transformed-output claims still require separate correctness validation.

Unchanged paths sometimes moved substantially. For example, the numeric screen's eight-of-256 subset load had block ratios of 0.65× and 1.56×, and an unchanged 32-column byte load had an anomalous second-block ratio near 0.23×. Such controls rule out attributing every timing movement to an optimization. They do not erase consistent losses on changed paths or prove those losses are harmless. The final report should identify regressions and uncertainty explicitly, rather than assert universal absence of regressions.

The recorded maximum RSS across the 48 component processes ranged from about 0.62 to 6.67 GB; all resource records reported zero swaps. Each measurement includes Stata, all fixtures in that process, verification, runtime state, and the wrapper. It is neither per-case nor plugin-only peak memory, and zero swaps is not proof of identical memory pressure across runs.

Round-trip tests cannot establish every computed-write behavior: a previous investigation found that identity writes concealed an unsafe callback substitution. This round retains checked stores, but acceptance still requires command-level tests that write transformed values, as well as native error, allocation-failure, duplicate-map, stale-width, ownership, and strL checks. There is no new cross-platform runtime evidence in these campaigns.

## Completed comparisons and acceptance

**COMPLETED — frozen transport, allocation, and allocation-toggle comparisons.** The tables above describe the frozen candidate and explicitly identify its regressions. **COMPLETED — gate disposition.** The final policy uses K ≥ 128 and selected N ≥ 50,000. **COMPLETED — final-source acceptance.** The final snapshot hashes, repeated full command suite, native runs, and final timing reruns are recorded. Preserve the distinction between these frozen measurements and the subsequently merged validation snapshot.

**COMPLETED — allocation timing question.** The [runtime-toggle experiment](no_order_toggle/manifest.json) supports modest numeric gains without the separate-plugin regressions. **COMPLETED — caller acceptance.** The opt-out is merged into 11 callers. Both full baseline and candidate command suites passed 4,161 checks with zero failures and the same 48 expected method-difference skips; exact merged-source hashes are preserved in the acceptance records.

**COMPLETED — CPU profiles.** The section above documents the changed load path and the limits of sampling counts. **COMPLETED — prior native boundary checks.** Numeric tiling, adaptive serial dispatch, and direct-gather checks passed in the tested snapshot. The final 50,000-row constant receives a new native/real-SPI confirmation.

**COMPLETED — numeric scheduler policy test.** The [same-binary experiment](numeric_scheduler_toggle/manifest.json) confirms the 64-column byte loss and supports the corrected 128-column minimum. **COMPLETED — 128-column merged-source comparison.** Its 648 records confirm the byte64 and dense-write fixes. **COMPLETED — final gate selection.** Short-byte follow-ups reject two- and six-million-cell thresholds. Use K ≥ 128 and selected N ≥ 50,000, which keeps all completed 12-case dispatch choices unchanged. **COMPLETED — exact-final-source transport confirmation.** The new snapshot, 336 real-SPI records, binary, full command suite, and eight native UBSan outcomes are recorded separately. Final whole-command timings and same-binary controls are complete.

**COMPLETED — native UBSan acceptance.** The four native suites (`test_transport_native`, `test_transport_scheduling`, `test_transport_adaptive_native`, and `test_transport_store_native`) all passed in serial and OpenMP configurations: **eight successful runs**, against merged candidate transport SHA-256 `c4c18472878f8b89767c4f15bc24448c8798d63f5f504a6197758e4e8381f1f4`. Their checks cover callback bounds, allocation failures, strL, sorted and duplicate maps, cancellation, and dispatch boundaries. [UBSan results](acceptance/final_native_checks/ubsan_results.json) identify the test-source hashes and logs.

**COMPLETED — command acceptance.** The full baseline and candidate suites each completed **4,161 passes, zero failures, and 48 expected method-difference skips**, with return code zero. The skips are expected method differences, not silently omitted failures. The full-plugin macOS ARM dependency check passed with minimum target 11 and compatible statically linked libomp. [Baseline acceptance](acceptance/production_baseline/acceptance.log), [candidate acceptance](acceptance/production_candidate/acceptance.log), and their [baseline](acceptance/production_baseline/source_hashes.json) and [candidate](acceptance/production_candidate/source_hashes.json) source maps preserve the tested snapshots. The earlier whole-command timings are complete. The final source and plugin hashes are recorded above, and the compiler-matched exact-final comparison is complete.

**Sanitizer limitation.** The initial AddressSanitizer runtime on this host stalls before test `main`: its initialization stack enters `InitializeShadowMemory`, dyld, malloc, and recursive sanitizer initialization. This is an unavailable ASan run, not a passing check or evidence that the test body executed. An independent Homebrew Clang21.1.8 probe compiled successfully but timed out after 15 seconds before `main`; no transport test body ran. The planned LLVM22 label was incorrect for the available alternative; the [probe record](acceptance/final_native_checks_final/asan_probe_result.json) captures the actual executable/version and timeout. Both runtime attempts are unavailable results. Linux CI sanitizer coverage has been added; a CI result must be recorded separately before claiming it passed. Preserve the initialization profile with the acceptance artifacts. The separate UBSan runs and both full command suites passed. [ASan diagnosis](acceptance/final_native_checks/asan_runtime_limitation.json) and [startup profile](acceptance/final_native_checks/asan_startup_profile.txt) preserve the failure evidence.

## Exact-final acceptance

The final K ≥ 128 / N ≥ 50,000 source and full plugin passed the full command suite again: **4,161 passes, zero failures, 48 expected method-difference skips**, return code zero. [Final acceptance](acceptance/production_final/acceptance.log) and [source map](acceptance/production_final/source_hashes.json) preserve the snapshot. The four native suites passed serial/OpenMP UBSan again, using explicit Apple `/usr/bin/clang`21.0.0 with the compatible libomp archive: [final UBSan records](acceptance/final_native_checks_final/ubsan_results.json). A comment-only test cleanup after execution is documented with the exact tested source copy. ARM dependency checks and Windows/macOS-x86 syntax checks passed; no cross-platform runtime test is claimed. [Final build checks](build_checks_final.json) and [performance toolchain](final_toolchain.json) distinguish compiler identities. The final performance baseline and candidate use Homebrew21.1.8/libomp21.1.8; the native sanitizer compiler is separate. The exact-final transport table reports the accepted gate and snapshot; earlier gate comparisons remain historical evidence.

## Delivery

The exact tested final ARM plugin was installed after all frozen source and ado files matched the live working tree. The previous binary was preserved as a backup, and the installed plugin passed its dependency check. [Installation metadata](install_final.json) records source and binary hashes and the previous-binary backup path. No later source change is attributed to these frozen results.
