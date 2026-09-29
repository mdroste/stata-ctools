# Same-plugin A/B causal control

Preparation freezes baseline and candidate data_io bodies as separate translation
units. All 13 externally defined functions receive base_ or cand_ prefixes.
Both use the final headers, types, arena, thread pool, SPI wrapper, and one OpenMP
runtime. The runtime selector chooses function pointers before load/store timing.

Each block uses one Stata process and the same plugin. Pairs alternate AB/BA;
block 2 reverses both the initial order and case order. Each fixture has one
warmup pair, 11 measured pairs, and one all-cell verification pair. Sorted stores
restore the original fixture outside timing, and Stata checks its datasignature
before the next call. All other cases check signatures after each verification
call and at case completion. Cleanup/free/pool destruction is timed separately.

No load-time environment override is introduced. Source namespacing changes code
placement, so this is not an identical-body flag-toggle control, and shared final
inline helpers mean this isolates the data_io bodies rather than every baseline
source change. Existing process-isolated comparisons remain separate evidence.

Compile only when other timing has stopped:
`python3 temp/transport_20260929/prepare_dual_impl.py build DIRECTORY`

Run the exact oldstata commands in run-plan.tsv sequentially, then:
`python3 temp/transport_20260929/prepare_dual_impl.py summarize DIRECTORY`
