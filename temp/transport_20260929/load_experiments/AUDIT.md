# Shared load transport audit (isolated experiments)

All candidates derive from `temp/transport_20260929/baseline/src`; only isolated copies of `ctools_data_io.c` were changed. The accompanying patches and `manifest.json` identify their changes. Source support files were copied only at the top level, sufficient for standalone transport benchmarks and native tests.

## Candidates

- `tiny_serial`: identity and filtered multi-column loads use the existing serial loop at <=65,536 selected cells. Does not change stores, large loads, or string allocation. Hypothesis: avoid a newly created pthread pool and work queue for submillisecond transfers.
- `numeric_tiles`: permits existing row-tiled loader for all-numeric K>=64; preserves existing N>=200,000 selected-row threshold. No read/write semantics change. Root's preliminary real-Stata result (not independently rerun here) was about3x faster at200Kx128.
- `str2045`: enables full-width flat slots at width2045. **Reject without stronger gating:**100Kx4, actual8-byte contents raised mock process resident memory from9.55MB to824.36MB. Declared-width flat buffers touch nearly every page even when values are empty or short.
- `str2045_packed`: retains2045 in parsed hints, keeps it outside flat slots, and reserves `min(N*2046,2GiB)` in each private arena with `malloc`, packing actual contents consecutively. Single-column, filtered and row-parallel arena paths agree. Other widths unchanged. Native RSS confirms only actual bytes are touched:100Kx4 actual8 bytes remains9.55MB; full2045-byte content is824MB versus828MB baseline. Reserves more virtual address space, and falls back to existing individual allocations if reservation fails or capacity exhausts.

## Correctness evidence

All four candidates passed `test_transport_native.py` and `test_transport_scheduling.py`, serial and OpenMP, UBSan. Dedicated `test_load_edges.py` passed numeric_tiles, str2045 and str2045_packed, covering K={1,2,3,4,20,64}; threads={1,4,12}; all-numeric, mixed, all-string; known2045-byte strings, empty/UTF8; identity/range/filtered loading; every cell/map/order; SPI read errors; invalid variable/range; and failures of the first three numeric tiled allocations with successful exact fallback. Threshold is lowered to512 selected rows in the test to exercise tiled paths with small fixtures.

`probe_packed_memory.py` checks every generated value and records process RSS, virtual size, load time and cleanup time into `memory_probe.json`. This is a native SPI mock allocation experiment; its timings are not Stata transfer measurements.

## Further opportunities and limitations

1. Existing standalone benchmark does not destroy its pool between calls, unlike the real dispatcher. Test cold invocation and teardown separately from warm loops before tuning tiny-work thresholds.
2. Tile size should be swept at K128/K512 and several numeric storage types. Physical host dataset width matters independently of selected K; omitted columns still determine Stata source row spacing. Sparse/random selected observations also differ from the existing dense two-thirds predicate.
3. The N>=200K gate excludes large-cell wide/short datasets; N10KxK2000 still has20M cells. An independently measured total-work gate may cover these without slowing genuinely small inputs.
4. Every load allocates/initializes a4N-byte identity sort order even for regression, export and destring callers that never sort. An opt-out flag can save memory and an OpenMP launch, but requires explicit caller/API changes.
5. `__ctools_strw` reader is limited to16,384 bytes, enough for only about3,276 maximum-width hints. Very large varlists silently lose flat-width optimization. Dynamic metadata or explicit array passing would avoid the cliff.
6. Known maximum widths are not necessarily actual lengths. Unknown/fallback arenas reserve64N bytes, then use per-string `strdup`; long actual strings expose an allocation/cleanup cliff. Packed2045 is a narrow fix; growing private blocks could cover all unknown-width workloads without enormous virtual reservation.
7. `ctools_data_load_ex` header documents an explicit widths array of length nvars, but implementation indexes by plugin-visible `var_indices[j]-1` and requires an SF_nvars-sized array. Current callers/tests use the latter. Correct documentation/API before encouraging external use.
8. Keep strL calling-thread serialization and checked numeric store callbacks. Existing real-SPI evidence documents crashes/corruption when these are parallelized or switched unchecked.

## Follow-up variants and checks

- `packed_retry` supersedes the allocator-only candidate: qualifies platform memory-accounting claims, and retries the original64N-byte arena if the larger reservation fails. `test_packed_retry.py` fails the large allocation and requires exactly one original-capacity retry for row/column, identity/filtered paths, then checks every string. Both original native suites and the dedicated test pass serial/OpenMP UBSan.
- `numeric_tiles_cells` uses K>=64, N_selected>=4096 and at least2M selected cells; existing string tiling remains N_selected>=200K. `test_numeric_cells.py` exercises exact2M-cell and4096-row boundaries, filtered equivalents, all values, SPI read failures and the first three speculative-allocation failures. Serial/OpenMP UBSan pass.
- `no_sort_order` adds opt-in flag0x02 and applies it only to11 demonstrable nonconsumer command families. Core and callers are separately reviewable in `no_sort_order_core.patch` and `no_sort_order_callers.patch`. `test_no_sort_order.py` verifies ordinary/default behavior, opt-outNULL order, numeric/mixed values, filtering/empty cases, repeated cleanup, ordinary stores, sorted-store rejection before writes, and preservation of the opt-out through tiled-allocation fallback. Both original native suites also pass serial/OpenMP UBSan. `width_array_docs.patch` is an optional independent correction to the explicit-width array contract and NULL auto-detection behavior.

Sort-order opt-out saves4N_selected bytes. For a single numeric loaded column, buffers fall from16N bytes (8N values,4N observation map,4N order) to12N bytes, excluding small metadata. Proven consumers excluded from the caller patch include every sort/sample/bsample/rangestat/ipolate load, merge keys, winsor by-groups, and rangejoin using data. Additional proven nonconsumers in mixed-purpose files (winsor target, merge keepusing, rangejoin master) were intentionally left unchanged for easier initial review.

A runtime caveat: unconditional order initialization wakes OpenMP workers before pthread-column loading. With active wait/default blocktime those workers can spin concurrently with pthread workers. No-order may also reduce this contention; test representative command runs under default as well as passive-wait environments.
