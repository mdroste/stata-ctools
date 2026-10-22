# Performance audit, 22 October 2026

The audit benchmarks all 20 public data/estimation commands, plus eight grouping
and grouped-command workloads, with **2 million and 5 million observations and
40 mixed-type variables in memory**. It identifies avoidable full-dataset passes
in the ado wrappers and implements focused changes in `cipolate`, `csplit`, and
`cwinsor`. Separate B–C–C–B campaigns confirm the changes at one and ten threads
(plus an oversubscribed twelve-thread check),
with an additional 10-million-row scaling campaign.

There are **no executable `egen` calls in the production ado files**. That work
was already moved into `_ctools_group` and `src/cgroup/cgroup_impl.c`. References
to `egen` in these files describe compatibility. The helper still has opportunities
to reduce full-row sorting and allocation; they are listed below.

## Confirmed improvements

All rows below use **ten threads on the ten-core host**. Seconds are medians of
eight measured calls per build at 2M/5M and six per build at 10M.

| Workload | Rows | Before (s) | After (s) | Before / after |
|---|---:|---:|---:|---:|
| `cipolate` | 2,000,000 | 0.2934 | 0.0518 | 5.67× |
| `csplit` | 2,000,000 | 0.2716 | 0.0393 | 6.91× |
| `cwinsor`, in place, 3 variables | 2,000,000 | 0.1621 | 0.0870 | 1.86× |
| `cwinsor`, by string group, 7 variables | 2,000,000 | 2.0122 | 0.3087 | 6.52× |
| `cencode` (unchanged control) | 2,000,000 | 0.0250 | 0.0237 | 1.05× |
| `cipolate` | 5,000,000 | 0.7475 | 0.1260 | 5.93× |
| `csplit` | 5,000,000 | 0.7029 | 0.1032 | 6.81× |
| `cwinsor`, in place, 3 variables | 5,000,000 | 0.4022 | 0.2175 | 1.85× |
| `cwinsor`, by string group, 7 variables | 5,000,000 | 5.1495 | 0.9554 | 5.39× |
| `cencode` (unchanged control) | 5,000,000 | 0.0575 | 0.0605 | 0.95× |
| `cipolate` | 10,000,000 | 1.4659 | 0.2601 | 5.64× |
| `csplit` | 10,000,000 | 1.3476 | 0.1928 | 6.99× |
| `cwinsor`, in place, 3 variables | 10,000,000 | 0.8432 | 0.4281 | 1.97× |
| `cwinsor`, by string group, 7 variables | 10,000,000 | 10.9162 | 2.1928 | 4.98× |
| `cencode` (unchanged control) | 10,000,000 | 0.1359 | 0.1472 | 0.92× |

One-thread speedups at 2M/5M were `cipolate` 2.40×/2.45×, `csplit` 2.41×/2.35×,
and grouped seven-variable `cwinsor` 3.46×/3.24×. The unchanged
one-thread encode control varied by 4–12%; at ten threads and 2M/5M it varied
by about ±5%. These controls are why the report does not credit small changes
on unchanged paths. Twelve-thread results are retained as an explicitly
oversubscribed secondary check, not the primary hardware configuration.


These are complete public-command times, including ado parsing, output allocation,
plugin computation, stores and cleanup. They are not isolated kernel speedups.
Each confirmation uses two separate processes per build, discards only the
predeclared warmup, and retains every measured call. `cencode` is an unchanged
control. The detailed CSVs include counts, minima and maxima.

## What changed and why

### `cipolate`

Previously the ado generated a selection variable, excluded missing x values,
initialized the destination, counted selected rows, and later counted missing
results. The plugin then scanned the selection variable again before loading
its inputs. On a wide dataset these are repeated strided passes over Stata's
records, even though the interpolation itself takes only a few milliseconds.

The wrapper now passes the original `if`/`in` directly. An explicit `skipif` flag
allows the loader to avoid predicate calls when there is no `if`. The C algorithm
already leaves missing x values missing; the plugin counts missing results from
its contiguous output and reports that count to the ado. Unrestricted calls
allocate the destination without filling it, since the checked plugin store
writes every row. Restricted destinations remain initialized to missing.

Long-string by-keys retain the earlier selection behavior: missing-x rows must
be excluded before loading a potentially unsupported strL value. Aliases where
x or y also occurs in `by()` retain the existing copies. The interpolation,
duplicate-key averaging and extrapolation formulas are unchanged.

### `csplit`

The wrapper no longer generates/counts a separate sample marker. The plugin reads
the original selection once, excludes empty input strings while tokenizing,
and reports no observations if all selected strings are empty. Whitespace-only
strings remain valid input. Unrestricted output variables are allocated without
an initialization pass. The write phase explicitly fills absent fields and empty
input rows, so **every cell of every uninitialized destination is written**.
Restricted calls retain initialized empty destinations. The native strL fallback
and delimiter/trim/limit/destring behavior remain in place.

### `cwinsor`

The largest problem was after the calculation: for each new output the ado copied
the source to a double, ran the plugin, invoked Stata's `recast` ado, and generated
and dropped an N-row probe merely to retrieve the default display format.
Inspection of native `recast.ado` confirms that it generates another variable,
counts changed values, moves it through the dataset, drops the old variable and
renames the replacement. Seven outputs repeated that process seven times.

For new floating-point winsorized outputs, the plugin now loads the **original
sources** into its existing double buffers and stores into separate destinations
created with their final storage types and default formats. Percentile calculation
therefore retains double precision. Extended missing source values are normalized
to ordinary missing in the new outputs, matching the old conditional generate.
The original variables remain unchanged. In-place calls retain the undo copy;
integer outputs, trimming, and overlapping source/by lists retain the existing
conversion path. Unrestricted loads also skip redundant predicate evaluation.

The verbose post-plugin timer now includes metadata and conversion work. Before
this change it started **after** that work and misleadingly reported almost zero
post-plugin time, even when it accounted for most of the command's runtime.

At 2M rows in the exploratory profiles, interpolation used about 0.005 s of C
computation inside a 0.299 s command; tokenization used about 0.008 s inside a
0.278 s command. Grouped seven-column winsorization used about 0.016 s of quantile
computation inside a 2.064 s command. These observations motivated removing
transport/wrapper work before attempting new numerical algorithms.

## Full command screen

| Command | Baseline, 2M (s) | Baseline, 5M (s) |
|---|---:|---:|
| `cexport` | 0.369 | 1.275 |
| `cimport` | 0.349 | 2.779 |
| `csort` | 0.260 | 0.878 |
| `cmerge` | 0.209 | 0.534 |
| `csample` | 0.167 | 0.427 |
| `cbsample` | 0.276 | 0.705 |
| `cencode` | 0.024 | 0.066 |
| `cdecode` | 0.016 | 0.036 |
| `cdestring` | 0.042 | 0.141 |
| `cipolate` | 0.303 | 0.818 |
| `csplit` | 0.273 | 0.799 |
| `crangejoin` | 1.003 | 3.142 |
| `cwinsor` | 0.170 | 0.413 |
| `crangestat` | 0.074 | 0.183 |
| `cbinscatter` | 0.319 | 0.676 |
| `creghdfe` | 0.382 | 0.928 |
| `civreghdfe` | 0.493 | 1.201 |
| `cpplmhdfe` | 2.739 | 8.438 |
| `cpsmatch` | 2.816 | 7.086 |
| `cqreg` | 2.140 | 5.644 |


These are exploratory baseline medians of three measured calls after one verbose
warmup. They locate costs; they are not controlled speedup estimates. The broad
campaign also ran the candidate and verified 56 matching all-variable data
signatures across 448 successful calls. Unchanged commands such as import/export
and sorting showed substantial timing variation between those sessions. Build
and small diagnostic work overlapped parts of this exploratory screen. The
confirmation campaigns ran sequentially without concurrent compilation or tests.
No observations were removed to improve a result.

## Remaining opportunities, ordered by evidence

| Area | Evidence and next experiment | Constraints |
|---|---|---|
| `cpsmatch` post-processing | At 2M, the plugin took 0.288 s, while ado preparation took 1.599 s and ATT work 0.839 s. The ado repeatedly runs `summarize`, generates squared weights and counts support cells. Fuse sufficient statistics while treatment, outcome and matching weights are already in C memory. | Preserve the weighted control mean, unweighted variances, squared reuse weights, singleton behavior, missingness and displayed results. Native logit/probit estimation remains a separate major cost. |
| `cqreg` inference | At 2M, auxiliary quantile solves took 1.335 s, versus 0.609 s for the main solver. Profile the weighted products and repeated solver loops at larger regressor counts. | The auxiliary fits already reuse the main state, skip unnecessary crossover and choose serial versus concurrent execution based on size/thread count. Simply parallelizing both fits again is not a justified improvement. Keep covariance definitions and convergence criteria. |
| `cpplmhdfe` | Six IRLS iterations dominated the tested workload: 2.352 s in the plugin at 2M versus 0.306 s in setup. Add finer timers around weighted projections, moment formation and mean/deviance updates, then test reuse of per-iteration workspaces. | Separation, singleton removal, numerical scaling and convergence are correctness requirements. No evidence here justifies changing tolerances. |
| `_ctools_group` mixed keys | At 2M, mixed-key numbering/counting cost about 0.28/0.26 s, compared with 0.024 s for numeric IDs. The count-only path still allocates and fills row IDs. Test count-only allocation removal and hashing repeated tuples before sorting unique representatives. | Preserve lexicographic numbering, extended missing exclusion, string byte ordering, label text and first representative semantics. The existing sort dispatcher already uses non-comparison algorithms where appropriate. |
| Shared selection/load transport | Loading and sample preparation remain material for regressions, range statistics and joins, especially when a few active columns sit in wide records. The two new fast paths demonstrate the value of explicit proven-unfiltered loads. | A sample cannot be inferred from a few probes. Never omit checks for a real `if`, missing regressors or invalid weights. Keep checked stores; prior transport work found corruption with unchecked writes. |
| `crangejoin` | At 2M, the verbose profile spent about 0.401 s loading the master, versus 0.056 s matching/counting. Investigate row staging and schema/allocation work before replacing interval search. | Preserve rollback, strL fallback and arbitrary join expansion. The campaign uses one match per row, not explosive joins. |
| `csample`, `cbsample` | The wrappers generate a sample marker and double flags/frequencies; `cbsample` then drops and expands complete wide rows. Optimize the unrestricted wrapper path and measure row movement separately. | Preserve exact seeded draws within each implementation, filtered observations, strata/cluster semantics and weight-only output. |
| `creghdfe`, `civreghdfe` | Dense integer FE remapping is already cheap on this panel; projection and SPI loading dominate. String-cluster preparation adds measurable wrapper work. Noninteger/sparse FE remapping still has full-array `qsort` fallbacks. | Benchmark sparse/noninteger FE codes before replacing these fallbacks; dense integer results do not measure them. |
| `cbinscatter` | Binning, data loading, and graph construction all contribute. `discrete` still sorts a full index array with `qsort`. | The measured default quantile case does not establish a discrete-path speedup. Preserve tie/bin and regression-adjustment semantics. |
| `csort`, `cmerge` | Whole-row transport and Stata sort-marker/schema work remain the main candidates. Existing row staging and adaptive sort selection already address the obvious C bottlenecks. | Stable ordering, missing keys and rollback rule out blindly removing the final Stata sort or snapshot. |
| `cimport`, `cexport` | CSV throughput is already parallel. Wide records, formatting, input preservation and filesystem cache state affect results substantially. | This campaign does not measure binary formats, Excel, compressed input, cold I/O or network filesystems. |
| `cencode`, `cdecode`, `cdestring` | Low-cardinality conversion is already inexpensive. `cencode` hashes strings and sorts only unique labels; `cdecode` has a direct mapped output path. Retain these as controls and test high-cardinality/long-label cases separately. | Do not trade away label, locale, extended-missing or fallback behavior for small baseline savings. |

## Correctness and open findings

Final Stata validation used the rebuilt plugin and synchronized ado files in
`build/`:

- `validate_cipolate.do`: 31 passed, 0 failed, 0 skipped.
- `validate_csplit.do`: 132 passed, 0 failed, 0 skipped.
- `validate_cwinsor.do`: 90 passed, **1 pre-existing failure**, 0 skipped; see below.
- `validate_performance_paths.do`: all 60 compatibility cases passed, covering
  `if`/`in`, empty samples, missing x, aliases in `by()`, long-string by-keys,
  empty/short string fields, output types/formats/labels, extended missing
  values, and rollback after output creation, partial stores and metadata work.
- The standalone winsor macro-rounding probe reproduced the documented difference.

Native `test_newcommands_native.py` and `test_split_tokens_native.py` passed with
AddressSanitizer and UndefinedBehaviorSanitizer. The OpenMP winsor transport
harness passed with UndefinedBehaviorSanitizer and **137 allocation/store failure
injections**, with no leaks reported by its allocation tracker. The new benchmark
evidence parser passed all eight tests for incomplete grids, duplicate/missing
checks, unsuccessful completion, invalid timings and mismatching outputs.
Package synchronization, command inventory and suite registry checks passed.
The numerical parity failure below means the complete winsor suite is **not green**.

The existing `validate_cwinsor.do` suite has one reproducible failure on **both**
the untouched matched baseline and the candidate: its 100K-observation float
fixture differs from installed `winsor2` on 1,000 lower-tail values. The exact
percentile is `-2.3204771280288696`. Stata's numeric local macro renders it as
`-2.32047712802887`, which lies on the opposite side of a float rounding midpoint.
The resulting outputs differ by `2.384185791015625e-7` (one float ULP), just beyond
the suite's seven-significant-digit threshold. The optimization intentionally
preserves existing numerical results; it does not relax this assertion or mask
the failure. See the standalone reproducer in `validation/probes/`.

A smoke probe also found that `creghdfe ..., absorb(string_variable)` returns
r(2000): its ado `markout` treats the string absorb variable as missing. The broad
campaign instead uses supported numeric absorbed effects and explicitly tests a
string **cluster** variable. Supporting string absorbed effects requires more
than removing that `markout`: the FE transport/remapping contract must be reviewed.
This audit does not claim to implement that separate feature/fix.

## Method, provenance and limits

Inputs use seed/sortseed 95183, a 200-period firm panel, continuous outcomes and
regressors, a Poisson count, binary treatment, a labeled category, numeric text,
and delimited text. The 40 columns comprise **32 numeric variables spanning
byte/int/long/float/double and 8 fixed-width strings** (widths 8–40). Unused payload
columns remain resident during every command. The seven-source grouped-winsor
case uses 1,000 string groups; interpolation uses firm groups and 30% missing y.
Input generation, disk save/reload and correctness checks are outside timing.
CSV import/export and join I/O intrinsic to the public commands are timed.

A separate `CLOCK_MONOTONIC` plugin measures elapsed wall time. Each process loads
only one ctools build, avoiding collisions between statically linked OpenMP
runtimes. Thread count is explicit; Stata's licensed processor count is separate.
The initial screen used twelve threads based on older benchmark notes; hardware
inventory established ten cores, so the primary confirmation uses ten threads.
Internal grouping helpers use their default thread setting in the broad screen.
Warm caches are retained. Data signatures cover every output cell and variable
name but do not prove equality of all metadata or estimation coefficients;
targeted independent tests check the changed commands' values and metadata.
Estimator sample counts and coefficient/covariance output are recorded, but no
estimator arithmetic changed.

Both builds use the same Apple Clang compiler, `-O3`, LTO, `-mcpu=apple-m1`, strict
floating-point flags (`-fno-fast-math -ffp-contract=off`), and the repository-pinned
LLVM 21.1.8 static OpenMP runtime targeting macOS 11. The runtime and CMake were
installed only under `temp/performance_audit_20261022/`; the system installation
was not changed. Baseline production sources were frozen before editing. SHA-256
manifests identify each actual plugin, ado file, timer and do-file. The packaged
binary present at task start was not used as a mismatched-compiler baseline.

Machine: Apple M4 Mac mini (Mac16,10), **10 cores: 4 performance and 6 efficiency**,
**16 GB RAM**, macOS 27.0 / Darwin 27.0.0; StataNow/MP 19.5 licensed for
four Stata processors. Apple Clang 21.0.0. The complete campaign contains
**1,112 timed calls: 860 measured repetitions and 252 declared warmups**.

This is an active-workstation, synthetic, warm-cache macOS arm64 campaign.
It does not establish performance on other platforms, arbitrary cardinalities,
very wide/long strings, hundreds of regressors, large join expansion or memory
pressure. Peak RSS was not measured. The correctness runs and native sanitizers
listed above are scoped to the changes; this is not an all-platform release gate.
The macOS ARM plugin in `build/` was rebuilt and verified. Other platform binaries
must be rebuilt from these sources before distributing the updated ado files.

## Reproduction and evidence

- Harness: `validation/benchmark_audit.do` and `validation/run_performance_audit.py`.
- Independent selection/output checks: `validation/validate_performance_paths.do`.
- Evidence: `docs/benchmarks/performance_audit_20261022/`.
- Local frozen builds, complete logs and dependency build logs:
  `temp/performance_audit_20261022/`.

Build before timing, using the same compiler/runtime for both snapshots. For each
thread count, invoke the following with labels `b1`, `c1`, `c2`, `b2`, substituting
the corresponding frozen build; do not run those sessions concurrently:

```sh
python3 validation/run_performance_audit.py \
  --build /absolute/path/to/frozen/build \
  --output /absolute/path/to/campaign --label b1 \
  --sizes 2000000 5000000 --threads 10 --reps 4 \
  --cases ipolate split winsor winsor_by encode
python3 validation/run_performance_audit.py \
  --output /absolute/path/to/campaign --summarize
```

The runner uses the repository's configured Stata launcher, verifies a fresh
completion nonce, freezes its inputs, refuses to overwrite an existing label,
and validates the complete planned measurement grid. Omit `--cases` for the full
28-workload screen. Use `--sizes 10000000 --reps 3` for the scaling confirmation.
