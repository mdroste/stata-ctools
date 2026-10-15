# Range-statistics performance investigation, 2 October 2026

`crangestat` now avoids repeated OpenMP dispatch for medium-sized groups, reuses
results at equal interval keys when `excludeself` is absent, and uses the linear
sliding-window scan when most interval keys are distinct. The final one-million-row
comparison measures 1.19–2.76× faster complete commands on the four targeted
workloads. Relative to installed `rangestat`, those commands run 50–77× faster.
These are workload-specific measurements on this Mac, not platform-independent
speed guarantees.

The investigation reviewed the existing sorting, merging, transport, and I/O
performance evidence, and inspected the range-statistics, winsorization, encoding,
and interpolation implementations. It concentrated implementation work on range
statistics, where phase measurements identified avoidable scheduling and repeated
computation. A guided-scheduling experiment in `cwinsor` was measured and rejected:
its apparent end-to-end gains also appeared on an unchanged code path, while the
winsorization phase showed no convincing benefit. Its original implementation is
retained. Single-group `cwinsor` remains slower than `gstats winsor` in these tests.

## Final measured commands

All rows below have 1,000,000 observations and use eight ctools threads. Seconds
are medians of 16 timed calls per build, collected in two separate processes per
build; warmups are excluded. Reference times pool eight measured calls across the
four processes. The first five rows use `rangestat`; the last uses `gstats winsor`.

| Workload | Before (s) | After (s) | Reference (s) | Before / after | Reference / after |
|---|---:|---:|---:|---:|---:|
| Mean, groups of 200 (control) | 0.02864 | 0.02823 | 2.04185 | 1.01× | 72.34× |
| Mean, groups of 1,000 | 0.07799 | 0.02825 | 2.10379 | 2.76× | 74.48× |
| Mean, groups of 5,000 | 0.04807 | 0.02817 | 2.15655 | 1.71× | 76.55× |
| Mean, one group with unique keys | 0.05494 | 0.04630 | 2.45659 | 1.19× | 53.06× |
| Median, groups of 200, ten equal keys | 0.05720 | 0.03713 | 1.86759 | 1.54× | 50.30× |
| Winsor, one group (unchanged control) | 0.05022 | 0.05024 | 0.03223 | 1.00× | 0.64× |

The unchanged mean and winsor controls differ by approximately 1% or less in this
confirmation. The earlier broad campaign had substantial time-varying background
load: its unchanged single-group winsor control appeared 1.42× faster. That
campaign is retained as exploratory evidence and is not used to credit the
rejected winsor change. No individual timing observations were removed.

The broad campaign includes 13 cases at one and eight threads: groups of 20,
200, 1,000, and 5,000; a single group; repeated keys; rolling medians; a skewed
case with 90% of observations in one group; and winsorization with one or eight
source columns. A preliminary screen used 200,000 rows. See the
[exploratory results](benchmarks/command_performance_20261002/main/summary.csv)
and [final confirmation](benchmarks/command_performance_20261002/confirmation/summary.csv).
The raw CSVs include all repetitions, warmups, and observed min/max times.

At 3,000,000 rows, the same B–C–C–B protocol with four measured calls per
session (eight per build) confirms the gains:

| Workload | Before (s) | After (s) | Before / after | Reference / after |
|---|---:|---:|---:|---:|
| Mean, groups of 1,000 | 0.25641 | 0.10955 | 2.34× | 58.61× |
| Repeated-key median | 0.19450 | 0.13443 | 1.45× | 42.71× |
| Unchanged winsor control | 0.14694 | 0.14657 | 1.00× | 0.63× |

The unchanged winsor control remains within 1%. All output comparisons pass.
[Scaling timings and ranges](benchmarks/command_performance_20261002/scaling/summary.csv).

## Why runtime fell

The original group selector used cross-group parallelism only when the average
group contained fewer than 500 observations. With 1,000-row groups it instead
started a parallel region for each group and dynamically scheduled individual
key runs. The final selector also processes groups together when the largest
group is below the existing 10,000-row within-group parallel threshold and there
are at most 16 source variables. It requires at least eight groups, as before.
The source-count guard preserves the existing min/max acceleration for wider
requests: the cross-group path has a 16-source sparse-table buffer.

For equal keys within a cross-group task, the window and its ordered contents
are identical. Without `excludeself`, the previous row's results can therefore be
reused for every requested statistic. With `excludeself`, each row is still
computed individually. Missing interval keys retain missing outputs, and reuse
never crosses a group boundary.

For large groups with more than one distinct run per two rows, the existing run
inventory is released and execution proceeds through the existing chunked linear
scan. This avoids a pair of binary searches and a dynamically scheduled task for
almost every row. Groups with many duplicates retain the run-based algorithm.
Neither optimization changes the arithmetic within a computed window.

In the final confirmation, median computation-phase times were:

| Workload | Before | After | Computation speedup |
|---|---:|---:|---:|
| Mean, groups of 1,000 | 50.53 ms | 0.79 ms | 63.7× |
| Mean, groups of 5,000 | 20.34 ms | 0.82 ms | 24.8× |
| Mean, one group | 13.02 ms | 4.97 ms | 2.6× |
| Repeated-key median | 29.86 ms | 9.82 ms | 3.0× |

These are internal computation times, not command speedups. Data transfer,
sorting, ado work, and preservation limit the total improvement. Five-second
macOS call-stack samples of a repeated 1M-row/1,000-row-group workload corroborate
the mechanism: the OpenMP dispatch algorithm appeared in 461 leaf samples in the
baseline and six in the final build. Sampling covers different numbers of
completed commands and many waiting threads; those counts are diagnostic, not a
wall-time ratio. Both complete profiles and their driver are retained.

## Correctness and a pre-existing bug

The final build passes:

- 177 existing `crangestat` tests and 57 existing `cwinsor` tests, with no skips.
- 257 new independent comparisons with `rangestat`, with the existing seven-digit
  tolerance unchanged. They cover group sizes around 64, 500, 1,000, and 10,000;
  unique and repeated keys; numeric/string/extended-missing group keys;
  `excludeself`; missing values and keys; reversed, empty, one-sided and unbounded
  windows; both sides of the run-count decision; and 16/17-source min requests.
- Exact all-row comparisons of all 23 statistics with the frozen baseline on
  eight repeated-key fixtures, covering 200/1,000-row groups, one/eight threads,
  and both `excludeself` settings.
- Full-output baseline comparisons and independent reference comparisons in
  every benchmark case/thread warmup, outside the timed region.
- The existing ASan/UBSan memory-safety suite, including 707 allocation/output
  failure injections, plus dedicated sanitized tiny-window moment checks.

The new boundary tests first failed 16 cases on **both** the frozen baseline and
the optimization candidate. `crangestat` unnecessarily required three observations
for skewness and four for kurtosis. The installed `rangestat` computes standardized
population moments: two distinct observations have skewness zero and kurtosis
one; a three-value window can also have defined kurtosis. The minimum-count
checks now allow two values, retaining missing results for zero-variance windows.
The moment formulas and the results for larger windows are unchanged. The baseline
failure log and diagnostic rows are retained; no assertions or tolerances were
relaxed to obtain the final passing result.

## Machine, builds, and timing method

Apple M4 Pro (eight performance and four efficiency cores), 48 GiB RAM, macOS
26.7.1; StataNow/MP 19.5, revision 12 August 2026, licensed for four Stata cores.
The installed references are `rangestat` 1.1.1 and `gstats` ado version 0.5.0.
ctools explicitly uses the thread counts recorded in each CSV; the reference
commands use their normal installed settings. This is not an equal-thread-count
comparison. One-thread ctools results are also available in the broad campaign.

Both plugins were rebuilt from the same uncommitted working-tree snapshot with
Homebrew Clang 23.1.2, `-O3`, LTO, strict floating point (`-fno-fast-math` and
`-ffp-contract=off`), `-mcpu=apple-m1`, and the macOS-11-compatible static OpenMP
archive already present in the checkout. Source, plugin, and ado SHA-256 inventories
identify the actual tested state; a Git revision alone would not describe this
checkout. Only the two range-statistics source files differ in the final build.
The exploratory candidate also contained the subsequently rejected winsor change.

Each baseline/candidate plugin loads in its own Stata process in B–C–C–B order.
The repository launcher uses the `stata` alias through interactive login zsh;
every driver exits with `exit, clear` and the runner verifies a fresh completion
nonce and return code. A separate `CLOCK_MONOTONIC` plugin measures the entire
public command. Input generation, reloads, and comparisons are outside the timer;
ado parsing, preservation, plugin calls/processing/storage, and command cleanup
are inside it. These are warm-cache, synthetic-data benchmarks on an active
workstation. Other desktop jobs were observed consuming CPU. Alternation, retained
raw timings, and unchanged controls help assess variability but cannot establish
isolated-machine performance, cold-storage throughput, peak memory, or results on
Windows/Linux/Intel hardware.

Inputs use seed 928143, shuffled original row order, double outcomes
`100 + 5*rnormal()`, ordinary missing values every 37th row and extended missings
every 107th. Numeric groups have the stated size; interval keys are unique within
groups except the tied workloads. Mean windows use `interval(key -5 5)`.
The median workload has 200 rows per group and ten rows at each interval key.
Winsor controls use cuts 1 and 99 and `suffix(_w)`. The source do-file specifies
all cases, storage types, and options.

## Reproduction and evidence

The [evidence directory](benchmarks/command_performance_20261002/) contains raw
measurements, profiles, successful and initial failing test logs, source hashes,
changed-source archives, the optimization patch, and the compiler/runtime record.
The baseline and final tested plugins remain in
`temp/command_performance_20261002/{baseline,final}`. The final sources use the
ordinary Makefile build path; no experimental scheduler flag is needed.

```sh
python3 validation/prepare_group_performance.py \
  --baseline temp/command_performance_20261002/baseline \
  --candidate temp/command_performance_20261002/final \
  --output temp/group_confirmation \
  --cases 'rs200 rs1000 rs5000 rs_single rs_median w_single' \
  --threads 8 --reps 8 --refreps 1
```

The benchmark helper is currently macOS-specific because it builds the clock as
a bundle. The production optimization itself uses the existing portable C and
OpenMP infrastructure. `validate_group_windows.do` and
`test_group_windows_native.py` are registered with the repository's validation
inventory. No existing validation assertions were changed.

The package-source check, generated command-inventory check, Python lint, and
scoped Git whitespace check pass. The rebuilt `build/ctools_mac_arm.plugin` has
the same SHA-256 as the tested final plugin. A whole-repository validation-registry
check is currently prevented by two unregistered files introduced during
concurrent work: `validation/test_coverage_helpers.do` and
`validation/validate_option_contracts.do`. They were not changed here. The focused
passing results above do not constitute a full release-gate run; the recorded
[registry diagnostic](benchmarks/command_performance_20261002/registry_check.txt)
identifies the outstanding inventory issue.
