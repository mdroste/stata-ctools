# Large-data csort and cmerge investigation

Baseline: commit `7e5c0c8`. Measurements use an Apple M4 Pro (8 performance and 4 efficiency cores), 48 GiB RAM, Stata/MP 18.5 with twelve processors, and eight plugin threads. These are measurements on this machine, not cross-platform speed guarantees.

The main bottlenecks were moving data between Stata and C, constructing merge output, and duplicate uniqueness checks in the ado wrapper. Sorting and the sorted join accounted for a much smaller share of most numeric workloads.

## Changes

- `csort` in standard mode now gathers small blocks into scratch space and writes them directly to Stata in the computed order. It no longer creates complete permuted copies of numeric and string-pointer columns. Numeric blocks use 8 KiB; string blocks use 64 KiB. Loaded columns and the sort permutation remain unchanged. The existing checked numeric SPI callbacks and validation still apply. The permutation work is now included in the store phase of verbose timings.
- `cmerge` output string arrays borrow their contents from the loaded master and using data, whose lifetimes extend beyond output storage. Only output pointer arrays are allocated; output cleanup runs before either owner is released. Read-ahead prefetching helps the store loop consume these scattered strings.
- Ordinary nonempty merges use the C join's existing complete uniqueness scans, including unmatched groups. The ado wrapper retains `isid` where the join is bypassed for empty inputs, and for explicit `sorted` requests. This removes two expensive duplicate checks from an unsorted 1:1 merge.

An initial direct, unblocked sort writer reduced scratch allocation but did not improve time consistently. It was replaced by the bounded-block implementation. Sort algorithm selection and automatic streaming thresholds were left unchanged.

## Measurement method

`validation/benchmark_bigdata.do` creates deterministic mixed data with byte, int, long, float, double, and fixed-string columns. The inputs contain a shuffled unique ID, 10,000 repeated integer groups, a string group key, a continuous secondary key, and deterministic payloads. A small m:1 lookup matches every master row. The 1:1 using file has the same row count as the master, with a 10% key shift: the result has 10% master-only, 10% using-only, and 90% matched input keys.

Every timed invocation reloads the same saved input first. The outer Stata timer measures the entire command, including ado work, preservation, using-file reads, and the post-sort metadata step. Input generation, the master-file reload, and correctness checks are outside that timer. In particular, `cmerge`'s existing internal verbose timer starts after its old master-side `isid`; relying only on that timer would miss a major baseline cost.

The plugin is warmed up before measurement. Each case runs three times. The first repetition checks every deterministic payload value in every output row; every repetition checks the row count and ordering or merge outcomes. The baseline and final processes run sequentially on identical saved inputs. File caches are warm; these are not cold-disk benchmarks. The desktop remains active, so small differences should be interpreted against the individual timings, not as precise universal gains.

Stata sessions use the `stata` shell alias through an interactive login zsh and exit cleanly. C phase timers use `mach_absolute_time`; the standalone kernels use `CLOCK_MONOTONIC`.

## Results

Medians of three repetitions; speedup is baseline time divided by optimized time. The m:1 workload runs 1.14–1.28× faster (13–22% less time); 1:1 runs 5.2–10.5× faster. `csort` results are mixed: up to 13% less time on the 20M-row integer-group case, near-neutral small cases, and 4–5% **more** time on the 10M-row, 40-column cases. The latter is a real observed tradeoff, not a claimed general sort speedup. Load-phase variability contributes to it, so the individual phase measurements are included.

The blocked writer is retained for its bounded scratch space and gains on some large cases. It should not be presented as faster for every dataset shape. No architecture-specific dispatch threshold was fitted to this small experiment.

Raw data: [end-to-end timings and phases](benchmarks/bigdata_end_to_end.csv), [ablations and configuration checks](benchmarks/bigdata_supplemental.csv), [100M-row kernels](benchmarks/bigdata_native_100m.csv).


### 1,000,000 rows, 20 columns

| Workload | Baseline seconds | Optimized seconds | Speedup |
| --- | ---: | ---: | ---: |
| sort_id | 0.124 | 0.117 | 1.06× |
| sort_group | 0.102 | 0.103 | 0.99× |
| sort_string | 0.138 | 0.136 | 1.01× |
| sort_float | 0.110 | 0.109 | 1.01× |
| merge_group | 0.150 | 0.117 | 1.28× |
| merge_id | 1.385 | 0.173 | 8.01× |

### 10,000,000 rows, 20 columns

| Workload | Baseline seconds | Optimized seconds | Speedup |
| --- | ---: | ---: | ---: |
| sort_id | 1.334 | 1.261 | 1.06× |
| sort_group | 1.213 | 1.197 | 1.01× |
| sort_string | 2.095 | 1.933 | 1.08× |
| sort_float | 1.246 | 1.212 | 1.03× |
| merge_group | 1.740 | 1.521 | 1.14× |
| merge_id | 18.963 | 2.469 | 7.68× |

### 10,000,000 rows, 40 columns

| Workload | Baseline seconds | Optimized seconds | Speedup |
| --- | ---: | ---: | ---: |
| sort_id | 2.535 | 2.632 | 0.96× |
| sort_string | 3.295 | 3.448 | 0.96× |
| merge_group | 3.752 | 2.942 | 1.28× |
| merge_id | 21.043 | 4.064 | 5.18× |

### 20,000,000 rows, 10 columns

| Workload | Baseline seconds | Optimized seconds | Speedup |
| --- | ---: | ---: | ---: |
| sort_group | 1.478 | 1.289 | 1.15× |
| sort_string | 3.329 | 3.152 | 1.06× |
| merge_group | 2.863 | 2.330 | 1.23× |
| merge_id | 37.746 | 3.595 | 10.50× |

### Whole-process memory

| Rows × columns | Baseline peak RSS, GB | Optimized peak RSS, GB |
| --- | ---: | ---: |
| 10M × 20 | 10.32 | 9.67 |
| 20M × 10 | 13.58 | 11.75 |
| 10M × 40 | 16.83 | 15.05 |

These are decimal GB from `/usr/bin/time -l` for the complete Stata process, including setup and all workloads, not a per-command allocation measurement. The 10M×20 inputs existed for both runs; the other baseline processes also generated their master input files. No swapping was observed during the matrix. At 10M rows and eight workers, the old generic permutation pool reserves up to 1.28 GB for numeric and pointer scratch arrays; the new writer uses at most 72 KiB per worker, independently of row count. Input columns still consume memory proportional to dataset size.

### What explains the gains

With the original C plugin and only the ado changes, the 10M×20 1:1 merge takes **2.715 s**, compared with **18.963 s** for the baseline and **2.469 s** for the combined changes. Most of this improvement therefore comes from eliminating redundant `isid` work. The payload changes matter more for m:1 merges, where checking the small 10,000-row lookup was already cheap.

At 20M×40, automatic streaming takes **4.527 s**, versus **6.015 s** with explicit `nostream` using the new writer. Both pass all-row checks. The existing automatic streaming choice remains useful and is unchanged; forcing full loading was about 33% slower in this check.

At 10M×20, the final numeric-group sort medians with 4/8/12 threads were **1.519 / 1.056 / 1.007 s**; m:1 merge medians were **1.884 / 1.339 / 1.347 s**. Eight and twelve threads are close here; four are slower. These configuration tests are separate from the main matrix.

The full matrix was collected before the two alignment fixes described below, on unsorted inputs with no shared payload names. A separate final-source recheck at 10M×20 takes 1.597 s for m:1 and 2.478 s for 1:1, retaining the gains. Those rows are labeled `corrected` in the raw data.


## Isolated 100-million-row kernels

These tests exclude Stata I/O, wide payload columns, file reads, and data generation. Workload codes 0–3 in the raw CSV follow the table order below. Each sort verifies the full permutation for bounds, uniqueness, and monotonic keys. The join verifies every emitted master row, using row, and match code.

| Workload | Median seconds, three repetitions |
| --- | ---: |
| Numeric sort, 10,000 integer groups | 0.0674 |
| Numeric sort, large integer keys | 1.0503 |
| Numeric sort, continuous doubles | 1.6551 |
| Sorted m:1 join, 100M master / 10M using rows | 0.1272 |

The sort uses the existing automatic algorithm choice. The sorted join is serial. These are scaling observations for unchanged kernels, not before/after speedups, and they do not establish full-command performance at 100 million rows.

## Reproduction

Use `validation/benchmark_bigdata.do` with an ado/plugin directory, output directory, label, row count, column count, repetitions, quoted case list, quoted sort options, and thread count. Put the call in a driver do-file so Stata receives quoted multiword arguments correctly:

```stata
capture noisily do "/absolute/repo/validation/benchmark_bigdata.do" ///
    "/absolute/plugin-directory" "/absolute/output-directory" "candidate" ///
    10000000 20 3 "sort_id sort_group sort_string sort_float merge_group merge_id" "" 8
local rc = _rc
if `rc' di "BIGDATA_ERROR,`rc'"
exit, clear
```

Invoke that driver through the machine's wrapper:

```sh
/bin/zsh -lic 'stata -q -b do "/absolute/driver.do"'
```

Cases are `sort_id`, `sort_group`, `sort_string` (string plus continuous key), `sort_float`, `merge_group`, and `merge_id`. `stream(8)` or `nostream` can be passed in the sort-options argument. Plugin variants need separate Stata processes because a loaded plugin is retained for that session.

The standalone kernel benchmark is independent of Stata:

```sh
python3 validation/benchmark_bigdata_native.py --rows 100000000
```

The end-to-end exploratory builds use Homebrew Clang 21.1.8 and the installed static OpenMP 22.1.8 archive, targeting the local macOS 26 machine. Both variants use `-O3 -fPIC -DSYSTEM=APPLEMAC -DSD_FASTMODE -fno-fast-math -ffp-contract=off -funroll-loops -ftree-vectorize -flto -fno-strict-aliasing -arch arm64 -mmacosx-version-min=26.0 -mcpu=apple-m1 -Xpreprocessor -fopenmp`, link the static `libomp.a` and Accelerate as a bundle, and discover all C sources and include directories under `src`. They are temporary benchmark artifacts, not distribution binaries. The final sources also build with the repository's standard macOS 11 Makefile target and pinned compatible OpenMP runtime; the dependency contract is checked separately. Tracked distribution binaries are not replaced by these experiments.

## Validation and scope

The final corrected plugin passes all **278 csort and 110 cmerge existing tests**, with no failures or skips. Additional checks cover unmatched duplicate keys, empty-side uniqueness, borrowed string keys, and 24 combinations of key type, sortedness, shared columns, and update/replace payload handling against native `merge`. The latter use `nogenerate` for update/replace, matching the existing suite's treatment of that option. Native regressions verify serial and parallel sorted stores, multiple block boundaries, unchanged input buffers, NULL/empty strings, invalid mappings, and SPI write failures. The existing P2 native suite also passes under UBSan.

The audit confirmed and fixed two pre-existing row-alignment bugs. A deferred key-only load could incorrectly skip master output when using-only keys inserted rows into sorted data. Such merges now fall back to full master loading and normal output assembly; a true identity output still avoids loading payloads. Shared master columns were also excluded from permutation, leaving values attached to the wrong keys. They now follow the master mapping before using-side update/replace rules are applied. Both minimal baseline reproducers failed comparison with native `merge`; both pass after the fixes.

No changes are made to the sort algorithms, supported merge types, or user-facing options. Automatic streaming remains relevant for more than 30 columns and more than 10 million rows; its permutation pipeline does not use the new standard-mode writer.

Full Stata runs at 100 million rows with 10–40 columns were not attempted on the shared 48 GiB machine. The native 100M-row checks establish kernel behavior only; they do not certify end-to-end memory requirements or speed at that size. Cross-platform benchmarks and release builds for Windows/Linux/Intel remain outside this experiment.
