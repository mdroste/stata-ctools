# Encoding and winsorization performance, 2 October 2026

This follow-up implements two measured improvements: `cencode` transfers label
bytes as data and installs the complete label vector in one Mata call;
`cwinsor` loads targets and grouping columns together, sharing selection and
row-order storage. The performance patch leaves calculations and the
existing-label fallback unchanged. The plugin and
`cencode.ado` must be installed together because their private label-file format
changed; the staged build contains the matched pair.

## Integrated checkout status

While the original frozen builds were being measured, concurrent edits changed
`cencode` metadata handling and `cwinsor` storage/missing/trim behavior, along with
other range-statistics code and reference tests. Those changes are preserved.
A second pair of frozen builds includes the same concurrent changes on both
sides, differing only in the performance patch. The following integrated results
therefore isolate the optimization even with the new metadata work present:

| Encoding workload (8 threads) | Before | After | Native encode | Before / after | Native / after |
|---|---:|---:|---:|---:|---:|
| 1,000 categories | 0.11156 | 0.07869 | 0.22776 | 1.42× | 2.89× |
| 10,000 categories | 0.40830 | 0.11122 | 0.42550 | 3.67× | 3.83× |
| 50,000 categories | 1.98293 | 0.53349 | 4.18942 | 3.72× | 7.85× |
| 2,045-byte strings | 3.21503 | 0.32198 | 0.20454 | 9.99× | 0.64× |
| Existing-label control | 0.74448 | 0.74229 | 0.31352 | 1.00× | 0.42× |

Times are seconds. Rows are one million except the 2,045-byte-string case
(100,000 rows). These are the most relevant
encoding timings for the integrated checkout. The earlier broader tables below
isolate the performance patch before the concurrent metadata changes.

The corresponding integrated winsor results (one million rows, eight threads)
retain the same load-phase improvement, although the concurrent ado metadata
work adds overhead to both builds:

| Workload | Before | After | gstats winsor | Before / after |
|---|---:|---:|---:|---:|
| Groups of 200 | 0.06991 | 0.05880 | 0.06803 | 1.19× |
| Groups of 1,000 | 0.06959 | 0.05839 | 0.06111 | 1.19× |
| Eight targets, groups of 200 | 0.43572 | 0.42235 | 0.37125 | 1.03× |
| One group | 0.10972 | 0.09093 | 0.05384 | 1.21× |

The unchanged range-statistics control differs by about 2% in this campaign.

A direct reference confirmation on that integrated snapshot uses two measured
calls after a warmup (the full B–C–C–B campaigns above remain the before/after
estimates):

| Workload | cwinsor | gstats winsor | winsor2 | winsor2 / cwinsor |
|---|---:|---:|---:|---:|
| grouped | 0.06036 | 0.06902 | 0.57900 | 9.59× |
| single | 0.06501 | 0.03210 | 0.77342 | 11.90× |
| wide | 0.41155 | 0.37725 | 6.60103 | 16.04× |
| noby | 0.07582 | 0.03995 | 0.21800 | 2.88× |

**Integrated correctness is not fully green.** The integrated snapshot passes
239 current `cencode` checks plus the 25 new exact-label checks. Its current
`cwinsor` suite passes 62 checks and fails the 100k-row test, where 1,000 values
differ from `winsor2` beyond the unchanged tolerance. The new transport suite
also stops at ungrouped singleton trimming: ctools returns missing while `gstats`
retains 27.473961. Both failures reproduce in the matched integrated baseline
with the performance patch reversed. The optimization is therefore not their
cause. The failing assertions are retained; their logs are included as
`attribution_baseline_*` and `integrated_test_*`. These failures should be
resolved in the concurrent compatibility work before treating the whole checkout
as release-ready.

## Encoding results

Complete-command medians in seconds. All cases below use eight ctools threads
and one million input rows, except the explicitly marked 100k-row long-string
case. The reference is Stata's native `encode`.

| Workload | Before | After | Native encode | Before / after | Native / after |
|---|---:|---:|---:|---:|---:|
| 10 categories (small label-set control) | 0.03669 | 0.03621 | 0.18921 | 1.01× | 5.23× |
| 1,000 categories | 0.07336 | 0.04195 | 0.23590 | 1.75× | 5.62× |
| 10,000 categories | 0.37612 | 0.07521 | 0.44360 | 5.00× | 5.90× |
| 50,000 categories | 1.95674 | 0.49953 | 4.27690 | 3.92× | 8.56× |
| 2,045-byte strings, 1,000 categories (100k rows) | 3.17843 | 0.27019 | 0.19696 | 11.76× | 0.73× |
| Unicode, 1,000 categories | 0.11633 | 0.07665 | 0.35459 | 1.52× | 4.63× |
| Filtered third of rows, 1,000 categories | 0.12759 | 0.09710 | 0.12603 | 1.31× | 1.30× |
| Existing label definition (unchanged fallback) | 0.74887 | 0.75195 | 0.31493 | 1.00× | 0.42× |

The largest gains occur when label creation dominated execution. The original
implementation generated several Stata/Mata statements per label, encoding each
byte as a numeric literal. The new length-prefixed file contains ASCII record
headers and raw label bytes, with no executable source. Mata validates the
complete file, closes it, then calls `st_vlmodify()` once with the vectors.
Quotes, macro-looking text, newlines, control bytes and Unicode remain data.

The C hash, sort and mapping algorithms are unchanged. For example, at 10,000
categories the median collection phase is 28.43 ms before and 28.31 ms after;
label-file writing falls from 8.64 to 1.20 ms. Most of the full-command saving
comes from eliminating the repeated interpreter work after the plugin returns.
For 2,045-byte strings, hashing still consumes about 209 ms, explaining why that
case remains slower than native `encode` despite the large improvement. Existing
label definitions still use the native fallback and its distinct-value count;
this path remains slower than calling `encode` directly.

## Winsorization results

One million rows, eight ctools threads. The reference is `gstats winsor`, a
stronger performance comparator than the original `winsor2` replacement target.

| Workload | Before | After | gstats winsor | Before / after | Reference / after |
|---|---:|---:|---:|---:|---:|
| Groups of 20 | 0.04666 | 0.03534 | 0.06290 | 1.32× | 1.78× |
| Groups of 200 | 0.04476 | 0.03413 | 0.07019 | 1.31× | 2.06× |
| Groups of 1,000 | 0.04430 | 0.03304 | 0.06280 | 1.34× | 1.90× |
| One group | 0.04987 | 0.03910 | 0.03208 | 1.28× | 0.82× |
| Eight targets, groups of 200 | 0.19276 | 0.17624 | 0.39512 | 1.09× | 2.24× |

For groups of 200, median loading time falls from 22.18 to 11.77 ms, while the
quantile phase remains 1.30 ms. Targets and grouping columns now share one
`if`/`in` selection and observation map. Group detection still sees unpermuted
keys, and stores compose the retained permutation with the original row map.
The quantile algorithm, arithmetic and dynamic group scheduler are unchanged.

The change removes an unused target identity permutation, the second observation
map, and the widened copy of the sort permutation: 16 bytes of explicit buffers
per selected row on this 64-bit build (about 15.3 MiB at one million rows).
Without `by()`, the unused identity permutation is also omitted. These are
allocation savings, not measured peak-RSS reductions. Single-group winsorization
still loses to `gstats winsor`; no universal speedup over the reference is claimed.

A separate direct comparison also times the original replacement target,
`winsor2`, on the same generated inputs (four measured calls after a warmup):

| Workload | cwinsor | gstats winsor | winsor2 | winsor2 / cwinsor |
|---|---:|---:|---:|---:|
| grouped | 0.03323 | 0.06970 | 0.61079 | 18.38× |
| single | 0.03960 | 0.03225 | 0.79938 | 20.19× |
| wide | 0.15487 | 0.38455 | 7.54682 | 48.73× |
| noby | 0.07553 | 0.05548 | 0.23623 | 3.13× |

These separate reference runs are not before/after estimates. Background load
increased toward the ungrouped case; all repetitions remain in the data. All
finite outputs and missing locations match. `winsor2` with `suffix()` normalizes
extended missings to ordinary missing, while ctools preserves its existing
extended-missing behavior; the comparison accounts explicitly for that known
representation difference outside the timed section.

At three million rows and eight threads, the encoding improvements persist:

| Workload | Before | After | Before / after | Native / after |
|---|---:|---:|---:|---:|
| 1,000 categories | 0.15635 | 0.12247 | 1.28× | 7.91× |
| 10,000 categories | 0.45631 | 0.15840 | 2.88× | 7.55× |
| Filtered third of rows | 0.16719 | 0.13588 | 1.23× | 2.68× |
| Existing-label control | 1.94124 | 1.96177 | 0.99× | 0.50× |

The existing-label control remains within about 1% at three million rows.

## Method and controls

Frozen baseline and candidate plugin/ado directories use the same compiler,
flags and OpenMP archive. The baseline includes all working-tree changes present
at the start of this follow-up, including the earlier range-statistics work.
In the original comparison, only the four intended C/header files differ
between the source snapshots. The integrated baseline was built from a copied
source tree with only this performance patch reversed; no live files were
reverted.
No unrelated working-tree changes were reverted.

The main campaigns run B–C–C–B in four sequential Stata processes. Each workload
and thread count has one untimed warmup and four recorded repetitions per
process, giving eight measurements per build. Reference commands have one warmup
and two repetitions per process, giving eight measurements overall. A small
`CLOCK_MONOTONIC` plugin times the entire public command; generation, save/use,
comparison and decoding work are outside the timed interval. All observations,
including warmups and extremes, remain in the evidence. No timing outliers are
removed. Plugin phases are supplementary diagnostics rather than command times.

One and eight ctools threads are measured. StataNow/MP 19.5 (12 August 2026
revision) is licensed for four cores. Hardware is an Apple M4 Pro with 12 cores
(8 performance and 4 efficiency), 48 GiB RAM, macOS 26.7.1. Builds use Homebrew
Clang 23.1.2, `-O3`, LTO, strict floating-point flags, `-mcpu=apple-m1`, and the
macOS-11-compatible static OpenMP 21.1.8 archive used elsewhere in this project.
These are workstation measurements with background activity, not isolated lab
measurements or platform-independent guarantees. This investigation starts no concurrent benchmark, compilation, sanitizer or
correctness jobs during its measured campaigns. Other work in the shared checkout
was active; frozen copies prevent code drift within comparisons but cannot
eliminate that work's CPU contention.

At one million rows the unchanged existing-label control differs by less than
1%; the ten-category case differs by about 1%. The unchanged `crangestat` control
in the winsor campaign differs by about 1% at one thread and 7% at eight threads.
The corresponding winsor loading-phase reduction is roughly 10–11 ms in both
baseline blocks, materially larger than the control's absolute drift.
At 1,000 rows, sub-millisecond timing and background variation are material:
the existing-label control varies even though its executed path is unchanged.
The small-data campaign is retained and is not used to claim a regression or
an improvement in that fallback. A follow-up with 200 measured calls per build
finds 2.249 vs 2.272 ms at eight threads and 2.266 vs 2.458 ms at one thread;
these small absolute differences are reported rather than credited to an
optimization of that unchanged path. A pilot had incorrectly dropped the variable
holding the existing label before saving its fixture; the final campaigns keep
that variable. The pilot is exploratory evidence, not the control result.

## Call-stack profiles

Four five-second macOS `sample` captures cover warmed loops of baseline and
optimized encoding and winsorization. The winsor baseline shows two distinct
load/selection branches; the optimized stack has one. The encoding baseline is
dominated by work in Stata's interpreter and string comparison, with little time
inside the plugin; the optimized profile exposes more hashing/copying work as
label interpretation shrinks. Stata's private symbols are mostly unavailable,
so phase timings and the source change provide the more precise attribution.
Waiting-thread samples and differing numbers of completed commands mean the
raw sample counts must not be interpreted as speedup ratios. Complete profiles
and the bounded-loop driver are retained.

## Verification

The original optimized frozen build passes 470 Stata checks without skips: 228 existing
`cencode`, 57 existing `cwinsor`, 25 new label checks, and 160 new winsor transport
checks. This performance change does not alter existing tests or numerical tolerances. New label checks
compare exact codes, entire value-label maps and decoded bytes against native
`encode`, including all nonzero byte values, Unicode, literal macro/quote text,
embedded CR/LF, 2,045-byte strings, filtering, empty samples, and 65,536 distinct
labels. Nine malformed-file cases must fail before installing any labels.
Winsor cases cover single/multiple targets, numeric/string/multiple group keys,
nearby doubles and extended missing keys, sparse samples, trim/winsor, one/eight
threads, and one/two/three-observation groups.

Every timed encoding workload is also compared on every row against both the
frozen baseline and native `encode`, including decoded label bytes. Every winsor
benchmark compares all outputs with the baseline exactly and with the reference
using the existing seven-significant-digit policy. The new native winsor harness
runs the production command, grouping and quantile code with a tracked allocator
and mocked transport, checks sparse row restoration, injects 108 allocation/store
failures, and retries successfully after allocation failures. It passes both
UBSan and combined ASan/UBSan with four OpenMP threads. The later integrated
harness includes concurrent matrix-output checks and passes 145 injections
under both UBSan and combined ASan/UBSan. The label writer passes
byte-framing and failed-open tests under both sanitizers. The existing broader
ASan/UBSan suite also passes its 707 command failure injections.

## Reproduction and evidence

Evidence directory: [encode_winsor_20261002](benchmarks/encode_winsor_20261002/).
It retains individual timings, summaries, phase records, logs, frozen-build hash
manifests, call-stack samples, drivers, and a patch relative to the starting
working tree. Machine-local datasets and frozen executables remain under
`temp/encode_performance_20261002/` and are not committed.

```sh
python3 validation/prepare_encode_performance.py \
  --baseline /absolute/frozen/baseline --candidate /absolute/frozen/candidate \
  --output /absolute/results --rows 1000000 --threads '1 8' --reps 4 --refreps 2
python3 validation/prepare_group_performance.py \
  --baseline /absolute/frozen/baseline --candidate /absolute/frozen/candidate \
  --output /absolute/winsor-results --rows 1000000 --threads '1 8' \
  --cases 'w20 w200 w1000 w_single w_wide rs1000' --reps 4 --refreps 2
```

All Stata launches use the `stata` alias through an interactive login zsh and
end with `exit, clear`. Each benchmark process validates a fresh completion nonce
and checks that neither frozen build changed. One Stata process loads only one
plugin version, avoiding shared-library reuse between builds.
