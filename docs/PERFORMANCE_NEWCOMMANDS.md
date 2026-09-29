# Interpolation, split, and interval-join optimization

Measured September 25, 2026 on an Apple M4 Pro (12 logical processors), macOS
26.7, Stata/MP 18.5. Builds used clang 21.1.8, the repository's `-O3`, LTO and
strict floating-point flags, and the macOS 11-compatible static OpenMP runtime.
These measurements compare the initial C implementations with their optimized
versions, **not with native Stata commands**.

For the README's comparisons with native `ipolate`, native `split`, and SSC
`rangejoin` on 1M–20M observations, see the
[native/SSC benchmark report](BENCHMARKS_README_COMMANDS.md).

## End-to-end results

Seconds below are medians of 12 measured calls per build, workload, and thread
count. The speedup column compares both builds at 12 threads. All observations,
including slow ones, are retained in the linked CSV files.

| Workload | Before, 12 threads | After, 12 threads | Speedup | After, 1 thread | After, 4 threads |
|---|---:|---:|---:|---:|---:|
| Interpolate: ordered numeric groups | 0.0625 | 0.0513 | 1.22× | 0.0693 | 0.0553 |
| Interpolate: shuffled numeric groups | 0.0942 | 0.0605 | 1.56× | 0.1131 | 0.0623 |
| Interpolate: one group | 0.1013 | 0.0682 | 1.48× | 0.0931 | 0.0694 |
| Interpolate: 95% in one group | 0.1078 | 0.0900 | 1.20× | 0.1305 | 0.0859 |
| Interpolate: shuffled string groups | 0.2157 | 0.1550 | 1.39× | 0.4612 | 0.1636 |
| Split: short strings, five fields | 0.1087 | 0.1060 | 1.03× | 0.2267 | 0.1131 |
| Split: long strings, 17 fields | 0.1750 | 0.1321 | 1.33× | 0.3375 | 0.1564 |
| Split: five delimiters, ten fields | 0.3764 | 0.2730 | 1.38× | 0.4815 | 0.2938 |
| Join: ordered numeric groups | 0.1011 | 0.0420 | 2.41× | 0.0741 | 0.0441 |
| Join: shuffled string groups | 0.1792 | 0.0860 | 2.08× | 0.3172 | 0.1068 |
| Join: about 20 matches per master | 0.0758 | 0.0358 | 2.12× | 0.0935 | 0.0416 |
| Join: 12 extra numeric payload columns | 0.1749 | 0.0979 | 1.79× | 0.3212 | 0.1415 |

The short-string result is essentially unchanged. More threads do not always
help: several interpolation workloads saturate around four threads, and a
single large interpolation group still has a serial scan. These are shared
workstation measurements, with noticeable variation between sessions; they are
not guarantees for other machines, data shapes, or memory pressure.

Interpolation workloads contain 1,000,000 rows with 60% missing y values and
`epolate`; grouped cases have 1,000 groups except the deliberately skewed case.
Split workloads contain 1,000,000 `str80` or `str180` rows, or 200,000 `str1800`
rows. Joins have 300,000 rows per input, except the dense case with 100,000 per
input. Output counts are 598,500, 896,865, 2,027,000 and 598,500 respectively.
The deterministic generators and exact options are in
[`benchmark_command_performance.do`](../validation/benchmark_command_performance.do).

## What changed and why

- **Shared ordering:** Already ordered data bypass sorting. Numeric keys use
  stable parallel radix passes, skipping invariant bytes. String and mixed keys
  use a stable merge with independent output tiles, including in the final merge
  stages that previously used very few workers. Numeric keys retain distinct
  extended missing values and treat signed zeros as equal.
- **`cipolate`:** Group discovery and output initialization run in parallel.
  The wrapper clones x or y only when needed because it also appears in `by()`.
  Duplicate-x averaging and the native endpoint-rounding rules are unchanged.
- **`crangejoin`:** A compact group index and contiguous sorted keys replace
  repeated comparisons of all group columns in both interval searches. Adjacent
  master rows in the same group reuse the group lookup. Output is partitioned
  into row tiles, so narrow results and one-to-many expansions use all workers.
  `restore, preserve` removes a redundant full master-file save/reload while
  retaining rollback and variable metadata.
- **`csplit`:** Delimiter lengths and first-byte lookup lists are cached. Parsing
  visits candidate delimiter positions once and preserves the last-listed tie
  rule. The write stage reuses tokens. Single-delimiter input is terminated in
  place; multiple-delimiter input is compacted within the loaded buffer. Both
  scan and write stages parallelize across rows.

The phase profiles support these choices. For example, shuffled numeric
interpolation reduced sorting from roughly 40 ms to roughly 12 ms in the first
profiles. String grouping also benefited from parallel group discovery. Join
searching and writing both improved; eliminating the redundant master-file
round trip materially reduced wrapper time. In split, the initial token-packing
implementation did not improve the long-string end-to-end result, prompting the
single-delimiter in-place path. Full phase observations are retained in
[`command_optimization_profiles.txt`](benchmarks/command_optimization_profiles.txt).

There are memory tradeoffs: numeric radix ordering adds an eight-byte key array
per input row, beyond the existing permutation scratch space. The join index
adds eight bytes per using key plus up to eight bytes per using row for group
boundaries. Split caches four bytes per selected row for token count and initial
offset, without allocating a separate array for every token. Stata allocation
and SPI data transfer remain significant costs.

## Measurement method and reproduction

Each experiment ran separate Stata sessions in baseline/candidate/candidate/
baseline order. Each session warms each workload, checks output, records phase
profiles, and measures six calls at each of 1, 4 and 12 threads after a warmup at
that thread count. Timers surround the entire command, including ado execution,
plugin loading checks, memory allocation, input/output and cleanup; dataset
restoration between calls is outside the timed region.

A benchmark-only plugin uses `CLOCK_MONOTONIC`. Stata's built-in timer was
excluded from the final analysis after the wrapper used for those runs changed
the system date during exploratory measurements. The C phase timers also use a
monotonic clock. Each process loads only one ctools build: attempting to load two
plugins with separate static OpenMP runtimes into one process stalled in an
OpenMP barrier, so that exploratory run was interrupted and discarded before it
produced timing records.

The first experiment measured all 12 workloads (864 observations). A second
experiment measured the final split refinement (216 observations); the table
uses that second experiment for split. Interpolation and join sources, wrappers,
and shared I/O were unchanged between the two optimized builds. All benchmark
output values matched the frozen baseline exactly. Raw observations and build
identities are available here:

- [First optimization experiment](benchmarks/command_optimization_initial.csv)
- [Final split experiment](benchmarks/command_optimization_split_final.csv)
- [SHA-256 identities of the measured builds and relevant sources](benchmarks/command_optimization_builds.json)

The baseline and candidates were built from a frozen source snapshot so
concurrent repository edits could not affect the comparison. Local snapshots
are retained under `/private/tmp/ctools-command-perf/` (`baseline`,
`candidate-v1`, and `candidate`). The final working-tree build also includes
unrelated concurrent changes; its timing is not substituted for the controlled
candidate measurements.

Given frozen `before/build` and `after/build` directories, prepare the experiment:

```sh
python3 validation/prepare_command_performance.py \
    --baseline /absolute/before/build --candidate /absolute/after/build \
    --output /private/tmp/ctools-perf
```

Run the four printed `stata` commands sequentially. Then collect the results:

```sh
python3 validation/summarize_command_performance.py /private/tmp/ctools-perf \
    --csv /private/tmp/ctools-perf/timings.csv
```

Both preparation and summarization accept `--cases sp_short sp_long sp_multi`
for the focused split experiment. Preparation compiles the timer and writes
drivers; it does not launch Stata. Summarization rejects missing completion
markers, failures, missing trials, and duplicate records; it does not trim slow
observations.

## Correctness and limits

The existing 241 command reference checks and 36 additional optimization checks
passed on the final rebuilt working-tree plugin. Additional cases cover radix-size inputs,
`by(x)`/`by(y)` aliases, extended missing groups, delimiter overlaps and limits,
UTF-8 strings, and a single master row spanning many output tiles. The native
UBSan and allocation-failure harness passed. Independent parallel-sort tests
compare against qsort at thread counts 1, 2 and 12, including a runtime-limited
two-thread team, stable incoming ties, signed zero, extended missings, and
partial tiles. The final in-place split path additionally passed native UBSan
checks for whitespace trimming, empty fields, UTF-8, limits, and maximum-width
strings. Linux CI runs the command and token harnesses with ASan and UBSan.

`csplit` retains its native `split` fallback for strL and uses native `destring`
for optional numeric conversion. `crangejoin` supports numeric and fixed-width
strings; strL must be recast losslessly first. No command syntax changed.

The broader working-tree suite recorded 3,931 passing checks, two failures, and
48 expected skips. The failures were in the separate PPML validation component
and the September 22 prediction/weight-expression regression block, while
concurrent PPML changes were in progress. The missing-`absorb()` assertion
expects error 198; an isolated call returns 198 on the untouched baseline and
0 on the current PPML wrapper, without invoking any optimized command. These
unrelated files and assertions were not changed as part of this optimization.

A subsequent isolated run of the complete September 22 prediction/weight block
returned error 111 on the frozen baseline and passed on the then-current working
tree. The full-suite record above remains a failed run; it is not reported as a
green release gate.
