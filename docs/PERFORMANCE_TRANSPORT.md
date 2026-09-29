# Stata–C transport investigation

This investigation changes the shared transport implementation in `src/ctools_data_io.c`. The strongest improvements come from using row locality for wide fixed-string data, using otherwise idle cores for transfers with few variables, and parallelizing filtered string writes whose destination rows are distinct. Narrow, many-column transfers retain the column scheduler. All numeric accesses still use checked SPI callbacks, and fixed-string reads still validate width hints before copying.

The comparison baseline is the **starting working tree**, including the preceding uncommitted sort/merge optimizations, rather than Git HEAD. The recorded HEAD was `7e5c0c83ca3d19242722806c5add0c56b701bc5a`. Concurrent changes to commands are excluded from the transport comparisons. [Source hashes and machine details](benchmarks/transport_provenance.json) and the [transport-only patch from that baseline](benchmarks/transport_implementation.patch) identify the implementations.

## Measurement protocol

Measurements use real Stata SPI callbacks on an Apple M4 Pro with 12 logical cores and 48 GiB RAM, macOS 26.7, and StataNow/MP 18.5 for Apple Silicon, revision 26 February 2025. Other Stata jobs ran on the machine; the user authorized continuing with timing uncertainty. These results establish performance on this machine, not a guarantee for every architecture, Stata release, or workload.

The standalone plugin times allocation/filter construction/loading separately from writes using `CLOCK_MONOTONIC`. It excludes input generation, permutation construction, restoration, and verification. Each case has a warmup, repeated timed transfers, and a separate verification pass that checks every loaded and written cell. It independently checks the observation map against `SF_ifobs`; the driver checks the restored Stata datasignature after each case. Numeric comparisons are bitwise, including ordinary and extended missing values.

The final matrix alternates baseline/control/candidate and candidate/control/baseline within the same Stata process. The control contains the earlier conservative optimizations and the strL fix, allowing the incremental scheduling changes to be assessed separately. Plugins export only the SPI entry points to avoid cross-plugin symbol interposition. `KMP_BLOCKTIME=0` and `OMP_WAIT_POLICY=PASSIVE` prevent the three statically linked OpenMP runtimes from competing while idle. These are comparison settings, not changes to production defaults. End-to-end command tests use normal runtime defaults and separate processes. Their total elapsed times use `validation/benchmark_clock.c` because Stata timers can jump when another wrapper invocation adjusts the system date.

Most cases have five timed repetitions; the small, empty, and ten-million-row groups have three. Reported speedups are baseline median divided by candidate median: 2× means half the elapsed time. Raw data retain warmups and verification passes, which the summaries exclude. Per-variant ranges and matched-repetition ratio ranges expose timing variation; they are not confidence intervals. Small movements on unchanged paths occurred, so differences near the noise level should not be attributed to an optimization.

Coverage includes 1,000 to 10 million rows; 1, 2, 3, 4, 20, and 40 variables; all five Stata numeric storage types; fixed strings of 1, 8, 32, 64, 244, and 2045 bytes across the final and exploratory matrices; known and absent width hints; ordinary, shuffled-pointer, sorted-gather, filtered, and all-selected checked transfers; dense, sparse, and empty selections; and 1, 8, and 12 threads. These are in-memory transport tests, not disk-loading benchmarks.

## Final transport comparisons

The final matrix completed 229 paired cases (4,281 transfer records including warmups and verification passes). Ten targeted follow-up cases add 510 records, with 15 timed repetitions each. Together with the exploratory runs, the archived evidence contains 8,875 transfer records. [Final raw data](benchmarks/transport_final.csv), [final summaries](benchmarks/transport_final_summary.csv), [follow-up raw data](benchmarks/transport_rechecks.csv), and [follow-up summaries](benchmarks/transport_rechecks_summary.csv) preserve the results.

Eight threads, known width hints; milliseconds are medians. Mixed datasets have alternating numeric and string variables. Filtered cases select two-thirds of rows.

| Input | Transfer | Before (ms) | After (ms) | Speedup | Paired ratio range |
|---|---|---:|---:|---:|---:|
| 1M × 20, mixed str32 | identity load | 82.0 | 26.4 | 3.10× | 2.85–3.60× |
| 1M × 20, mixed str32 | identity store | 61.5 | 27.1 | 2.27× | 2.12–2.58× |
| 1M × 20, mixed str32 | filtered store | 171.7 | 45.1 | 3.81× | 3.42–3.86× |
| 1M × 20, mixed str244 | filtered store | 1105.6 | 157.1 | 7.04× | 6.74–7.08× |
| 2M × 40, mixed str32 | identity load | 405.4 | 105.3 | 3.85× | 3.44–4.13× |
| 2M × 40, mixed str32 | identity store | 220.7 | 112.6 | 1.96× | 1.91–2.24× |
| 1M × 2, mixed str32 | identity load | 14.2 | 4.3 | 3.27× | 2.84–3.78× |
| 1M × 2, mixed str32 | identity store | 19.6 | 3.7 | 5.29× | 4.56–5.91× |
| 10M × 20, mixed str32 | identity load | 971.3 | 257.2 | 3.78× | 3.62–3.89× |
| 10M × 20, mixed str32 | identity store | 514.5 | 230.1 | 2.24× | 2.20–2.43× |
| 10M × 20, mixed str32 | scattered store | 712.7 | 491.6 | 1.45× | 1.40–1.61× |
| 10M × 20, mixed str32 | filtered store | 1767.4 | 372.8 | 4.74× | 4.73–4.77× |

These are selected representative cases, not an aggregate speedup across commands. The 15-repetition follow-ups put single-column loads/writes within about 2% of baseline; 20-column numeric loads improved 3–5% and writes 0–4%. Three-column numeric transfers improved 1.5–2.8×. Earlier slower medians in those short cases did not persist. A few large, unchanged write paths still had roughly 6–7% slower medians in the main matrix. Their repetition variation and concurrent load preclude a universal no-regression claim. The large wide-string gains are substantially larger than that noise.

The [matrix manifest](benchmarks/transport_final_manifest.json) records the exact generator arguments. [Resource measurements](benchmarks/transport_resources.json) record peak process RSS and swaps for each final group; each process includes all three variants, so these are not per-variant memory comparisons.

## Retained changes

| Change | Applicability and safeguards |
|---|---|
| Classify fixed strings and strL outside the cell loop; cache the fixed-string callback | Keeps checked error mapping and bounded reads. Metadata is queried per column rather than per cell. |
| Reserve string bytes directly in worker-private column arenas | Removes atomic reservations only where one worker exclusively owns the arena. General shared arenas remain thread safe; exhaustion still uses the existing individually owned fallback allocations. |
| Write directly from prefetched string pointers | Removes a full-column repack for ordinary string writes. Preserves the former width checks and 256 MiB applicability boundary. Sorted writes retain their bounded gather buffer. |
| Keep a compact selection bitmap after the first excluded row | Avoids evaluating `SF_ifobs` twice for selected rows. Scratch space is about one bit per remaining row and is freed before loading columns. Allocation failure falls back to the previous two-pass behavior. |
| Tile wide fixed-string loads and writes | At least four variables, 200,000 selected rows, known widths below 2045, and at least 256 declared string bytes per row. Loads use 512-row tiles; writes use 4096-row tiles. Workers own disjoint output ranges. Load allocations retain the 2 GiB per-column bound and retry the existing allocator on failure. |
| Use row parallelism when few columns leave cores idle | Loads of two or three variables, and writes of two to four, use spare workers when at least twice as many threads as variables are available and there are at least 200,000 rows. Four-variable loads retain the existing schedule because their results were inconsistent. |
| Parallel filtered string writes | Only for sufficiently large, strictly increasing observation maps. Validate the complete map before writing. Duplicate or reversed destinations retain sequential last-write-wins behavior. |
| Unroll checked numeric column loops four cells at a time | Retains each callback check and the first-error return within the column. This is a modest optimization, not the source of the large string-transfer gains. |

New tiled and few-column scheduling, and parallel filtered string dispatch, are enabled only in OpenMP builds. Builds without OpenMP retain the existing pthread column scheduler.

## Correctness finding: strL

Real-SPI testing exposed incorrect values in the starting implementation's parallel strL reads. An early optimized candidate also triggered an allocator abort inside Stata. Serial reads passed. The final implementation therefore reads strL columns on the calling thread before dispatching fixed-string and numeric work; single-column and few-column paths follow the same rule.

This is a correctness fix, and no speedup is claimed for strL. The existing support boundary remains bounded textual strL up to 2045 bytes; binary or oversized values and strL writes are rejected. The [Stata SPI documentation](https://www.stata.com/plugins/) defines separate length-aware strL accessors; it does not provide a documented bulk array transport interface used by this implementation.

The dedicated real-Stata regression covers empty and UTF-8 text, widths through 2045, reordered mixed varlists, ranges and filters, one to three columns above the production row-parallel threshold, thread settings 1/8/12, oversized read rejection, and write rejection. It checks every value.

## Approaches tested and rejected or restricted

The [experimental raw timings](benchmarks/transport_experiments_raw.csv), [summaries](benchmarks/transport_experiments_summary.csv), and [build commands](benchmarks/transport_experiments_manifest.json) preserve the screening and follow-up results. The screen compared 13 variants in one process; subsequent experiments isolated callback caching, numeric unrolling, filtered writes, conditional tiling, and few-column scheduling.

| Approach | Decision |
|---|---|
| Row-parallel loads for every column in a wide numeric dataset | Rejected as a general policy: sequential parallel-region launches lost to column parallelism. Restricted to two/three columns with spare cores. |
| Row-parallel stores for every dataset | Rejected as a general policy: regressions for many numeric/narrow columns. Restricted by width or low column count. |
| Load tiles of 512, 4096, and 32768 rows | 512 retained behind the wide-data gate. Large tiles and unconditional tiling were inconsistent or slower. |
| Store tiles of 512, 4096, and 32768 rows | 4096 retained behind the dispatch gates. |
| Direct pointer gathers for sorted string stores | Rejected; the existing small gather buffer performed better. |
| Replace zeroed flat buffers with uninitialized allocation | Rejected; no stable improvement justified changing initialization behavior. |
| Parallelize all filtered destination maps | Rejected as unsafe for duplicates. Only strictly increasing maps use parallel writes. |
| Direct SPI reads into slots using untrusted width hints | Not adopted: a stale width can overwrite the next slot. Temporary bounded buffers remain. |

The private-arena change has supporting CPU-profile evidence. In the baseline, worker samples concentrate around the arena reservation loop; disassembly places `casal` at `load_variable_thread+1180`, immediately before a heavily sampled offset. The conservative candidate removes that reservation from private arenas, leaving more samples in Stata string copying and `strlen`/`memcpy`. These are qualitative profiles of an earlier candidate, not percentages of final wall time: [baseline profile](benchmarks/transport_profile_baseline.txt), [candidate profile](benchmarks/transport_profile_candidate.txt), [baseline disassembly](benchmarks/transport_profile_baseline_disassembly.txt).

## End-to-end command checks

Five repetitions per variant on the same saved **two-million-row, 21-variable** dataset (20 mixed payload variables plus an original-position identifier; nine `str32` columns), with eight threads and default OpenMP wait settings. Timings use the independent monotonic clock and include the whole command; loading the saved input and the result assertions are outside the timer. Both variants use identical frozen sort/merge code and ado files, differing only in shared transport. The small smoke cases and all large runs passed ordering, match, and value assertions.

| Command | Before median (s) | After median (s) | Time reduction | Before range (s) | After range (s) |
|---|---:|---:|---:|---:|---:|
| Numeric-key sort, nostream | 0.490 | 0.277 | 43.5% | 0.474–0.728 | 0.253–0.445 |
| String/numeric-key sort, nostream | 0.781 | 0.480 | 38.5% | 0.622–0.852 | 0.468–0.875 |
| String/numeric-key sort, stream(4) | 0.779 | 0.625 | 19.8% | 0.747–0.816 | 0.606–1.251 |
| m:1 merge | 0.644 | 0.469 | 27.1% | 0.549–0.705 | 0.392–0.858 |

[Command timings](benchmarks/transport_commands.csv) and [summaries](benchmarks/transport_commands_summary.csv) retain all repetitions. These command comparisons ran in separate baseline and candidate processes, so they are more exposed to changing background load than the alternating standalone transport matrix. They support practical command-level gains but do not establish identical gains for every ctools command.

## Validation outcome

- The final transport source passed the complete frozen original command suite: **2,877 passed, zero failed, 48 expected skips**.
- The broader repository snapshot produced **3,966 passed, one script failure, 48 expected skips**. The failure occurred at `which rangejoin`: the external SSC reference command was missing, so the `crangejoin` comparison never ran. This is an explicit coverage gap, not a successful test.
- Both new native suites passed with UBSan in serial and OpenMP builds. Existing SPI fault-injection and fused sorted-store regressions also passed against both full source snapshots. Coverage includes bounds, stale widths, allocator failures, selection bit boundaries, reordered variables, UTF-8, ownership, and repeated/reversed filtered destinations.
- The exact final minimal plugin passed the real-Stata strL regression, including large one/two/three-column cases at 1/8/12 threads and rejection checks.
- Both isolated macOS 11 deployment-target plugin builds passed dependency checks. Intel macOS OpenMP syntax checking and Windows x64 compilation without OpenMP passed; these are compile checks, not cross-platform runtime evidence.

[Validation markers](benchmarks/transport_validation.txt) and [full-build source/binary hashes](benchmarks/transport_acceptance_sources.json) preserve the evidence. All transport source and support files in the working tree matched the tested integration snapshot at closeout. Existing validation assertions were unchanged. The source change is applied; tested full plugins were built in isolated directories, and the shared distribution plugin was not overwritten while other command work was ongoing.

## Reproduction

`validation/benchmark_transport.py` builds isolated plugins and writes a do-file; it never invokes Stata itself. Supply frozen source directories for each variant. For the standalone benchmark, the starting baseline can be reconstructed by reversing the recorded transport-only patch in a copy of the final source; the shared transport dependencies were byte-identical in both variants. Verify the recorded source hashes before comparing. For example:

```sh
python3 validation/benchmark_transport.py \
  --variant baseline=/path/to/baseline/src \
  --variant candidate=/path/to/candidate/src \
  --output /private/tmp/transport-check \
  --rows 1000000 --columns 20 --repetitions 5 \
  --threads 8 --shapes numeric mixed --widths 8 32 244
```

On this machine, invoke the generated do-file through the `stata` shell alias:

```sh
KMP_BLOCKTIME=0 OMP_WAIT_POLICY=PASSIVE /bin/zsh -lic \
  'stata -q -b do "/private/tmp/transport-check/run.do"'
```

The driver exits with `exit, clear`. Durations use the monotonic clock. `validation/summarize_transport.py` collects only completed logs. `validation/test_transport_native.py` and `validation/test_transport_scheduling.py` exercise serial and OpenMP paths with UBSan locally and support AddressSanitizer plus UBSan on Linux CI. `validation/benchmark_transport_commands.do` supplies the wide mixed end-to-end sort/merge workload. It expects `clock.plugin` in its output directory, built from `validation/benchmark_clock.c` and `src/stplugin.c` with the platform SPI defines.

AddressSanitizer could not run locally: both available Clang runtimes stalled during sanitizer initialization before reaching the test program. No local ASan pass is claimed. Linux CI coverage was added, but is not a result from this machine. Runtime and performance claims are limited to the tested Apple Silicon/Stata combination.

## September 2026 follow-up: unchecked reads, direct-slot reads, chunked scheduling

A second investigation (September 27, 2026; Apple M4 Pro, 12 logical cores,
48 GiB, StataNow/MP 18.5) revisited three ideas the first campaign had not
tested: unchecked SPI accessors, direct-into-slot fixed-string reads, and
sub-column task scheduling. The final design is in `src/ctools_data_io.c`.
Measurements used `validation/benchmark_transport.py` with the same protocol,
plus single-variant isolated Stata processes for every string-path decision
(see "Measurement caveat" below).

### Retained changes

| Change | Design and safeguards |
|---|---|
| Unchecked SPI reads (`IO_VDATA_FN` → `(_stata_)->vdata`) | Every public entry point validates variable indices and observation ranges before its loops, so the per-cell bounds checks were redundant; per-cell return codes are still checked. Numeric loads improved 1.10–1.23x; strings unaffected. `-DCTOOLS_IO_CHECKED_SPI` restores the checked read callback. The native SPI mock now provides both accessor pairs, exactly like Stata. |
| Direct-slot fixed-string reads (`read_string_slot_direct`), widths >= 32 only | Flat and tiled loads read through SPI directly into the destination slot, skipping the bounce copy. Flat buffers gain whole guard slots (`calloc(nobs + guard_slots, stride)` keeps the fault-injection-visible allocation shape), and each worker's direct reads are bounded by the end of the region it exclusively owns, so a stale width hint can never write outside owned memory; rows too close to a region boundary keep the bounced read, and oversized values still abort the load through the same post-read length check. Isolated str32 tiled identity loads improved up to 1.5x and str244 loads 1.06–1.13x; below width 32 the bounce copy is effectively free and direct reads measured slightly slower, so narrow slots keep the bounce. |
| Chunked column scheduling (`try_chunked_load`, `store_chunk_thread`) | One pool task per column leaves workers idle when the column count divides the workers badly. Loads chunk only all-numeric data with at least 1.5M-row chunks: 13 numeric columns on 8 threads improved 1.2x, while 500K-row chunks were a wash against the allocation pass and task overhead. Stores chunk every type (about `4*threads/nvars` chunks of at least `MIN_OBS_PER_THREAD` rows); isolated sorted stores improved 1.07–1.23x with no string regressions. The tiled and few-column paths keep precedence, strL keeps the serial calling-thread rule, and a failed batch submit aborts rather than re-running chunks. |

### Correctness finding: the unchecked store corrupts columns under concurrent writes

The first candidate also switched writes to the unchecked `(_stata_)->store`.
A single-threaded micro-test showed the two store callbacks bitwise identical,
including float/int truncation and out-of-range-to-missing conversion, and the
standalone transport matrix (which stores loaded values back) verified clean.
Real-Stata validation then failed `cwinsor` and `cipolate`: commands that
store *computed* columns concurrently across variables (cwinsor stores through
an OpenMP parallel-for over `ctools_store_filtered`). Whole columns were
redirected into neighboring variables deterministically — `price` received
`mpg`'s column, `mpg` received `weight`'s — while single-threaded semantics
remained correct. Reverting only the store callback fixed it; unchecked reads
remained clean under every parallel pattern tested. **Writes therefore stay on
`safestore` permanently** (`IO_VSTORE_FN` is not switchable), and the source
carries a warning. The identity/round-trip structure of the standalone
benchmark could not catch this; only command-level validation did.

### End-to-end command effects

Five repetitions per variant, 2M x 21 mixed dataset
(`validation/benchmark_transport_commands.do`), both plugins built from the
same tree differing only in the transport module, separate Stata processes.
Whole-command medians: numeric-key sort 1.13x, string-key sort 1.06x,
stream sort 1.04x, m:1 merge 1.09x. Phase medians: loads 1.18–1.39x,
stores 1.06–1.14x. All ordering/match/value assertions passed, and the full
validation suite matched the baseline exactly (4108 passed; the same 3
pre-existing import/export failures as the unmodified tree).

### Measurement caveat: multi-plugin processes

Alternating several plugins inside one Stata process (the first campaign's
protocol) showed deterministic 1.2–1.5x differences on *identical, unexecuted
code paths* in some string store cells, with the sign depending on the loaded
variant set; single-variant isolated processes showed parity on those same
cells. Same-process comparisons remain useful for tightly-scoped numeric
deltas; every string-path decision above was confirmed or overturned in
isolated processes. One overturned example: same-process runs suggested
chunked and direct-slot str8 loads were fine or better, while isolated runs
showed a consistent 0.85–0.95x regression — which is why chunked loads are
numeric-only and direct reads have the width-32 floor. Wide-string matrices at
2M x 20 x str244 approach memory pressure on a 48 GiB machine and swing
+/-25%; 1M rows is the reliable scale for that width.

### Rejected in this round

- Unchecked stores (see the correctness finding above).
- Direct SPI reads into private arenas and per-thread slabs (width-unknown
  strings): the destination of each read depends on the previous row's length,
  chaining SPI call latency; measured 0.61x on mixed str8 no-hint loads.
- Direct-slot reads below width 32, and chunked loads for any string-bearing
  varlist (isolated 0.85–0.95x on str8).
- Chunking all-numeric loads below 1.5M-row chunks.
- `(_stata_)->data` (direct-return numeric accessor): would drop per-cell
  return codes for at most a marginal gain over the unchecked pointer form.

The native regression suites (`validation/test_transport_native.py`,
`validation/test_transport_scheduling.py`) were extended: the shared mock now
provides the unchecked accessors, a 12-thread case exercises the chunked
store scheduler, and a numeric-only case (with `IO_CHUNK_MIN_NUMERIC_ROWS`
overridden) exercises the chunked load path, identity and filtered, plus
error injection. Both pass with UBSan in serial and OpenMP builds.
