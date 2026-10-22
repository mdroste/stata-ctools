# Performance audit evidence, 22 October 2026

The findings and limitations are in [the audit report](../../PERFORMANCE_AUDIT_20261022.md).
All inputs are synthetic. Every timed workload starts with 40 resident columns:
32 numeric variables across all five storage types and eight fixed-width strings.

| Campaign | Rows | Threads | Order | Calls including warmups | Purpose |
|---|---|---:|---|---:|---|
| `broad` | 2M, 5M | 12 | B, C | 448 | All 20 commands and eight additional workloads |
| `confirmation_t1` | 2M, 5M | 1 | B, C, C, B | 160 | Three changed workloads and unchanged encode control |
| `confirmation_t10` | 2M, 5M | 10 | B, C, C, B | 200 | Primary host-sized confirmation, including in-place winsor |
| `confirmation_t12` | 2M, 5M | 12 | B, C, C, B | 160 | Oversubscribed secondary confirmation |
| `scaling10m` | 10M | 12 | B, C, C, B | 64 | Oversubscribed secondary scaling |
| `scaling10m_t10` | 10M | 10 | B, C, C, B | 80 | Primary host-sized scaling, including in-place winsor |

B is the untouched source baseline at `e470fdab096ee3afc14598a7ff52e8ef532c03f4`;
C is the optimized build. Both use the same compiler, flags and static OpenMP
runtime. Runs in one campaign are sequential. The broad screen is exploratory;
only the separate confirmation campaigns support the reported speedup claims.
Grouping helpers do not accept a thread override and use their hardware default.

Each campaign contains:

- `raw.csv`: every elapsed time, including repetition 0, the declared verbose warmup.
- `summary.csv`: medians, minima and maxima of all measured repetitions; no outliers removed.
- `checks.csv`: complete-data signatures and group counts from every warmup.
  Signatures match across builds for every workload/size/thread combination.
- `profiles.txt`: the verbose warmup sections, including phase timings and, where
  applicable, estimation coefficients/covariances. These matrices are recorded,
  not automatically compared by the timing harness.
- `*.manifest.json`: each process's planned grid and SHA-256 hashes of its frozen
  plugin, ado files, do-file, timer, Python runner and Makefile. Paths identify the
  original local files; they are not dependencies of the CSVs.

`source_manifest.json` identifies the hardware, compiler, source snapshot and
candidate plugin. `correctness/` retains the final independent Stata and native
test outputs, including the pre-existing winsor parity failure on both builds.
The passing rounding probe demonstrates that difference; it does not turn the
failing parity test into a pass.

Complete process logs, frozen builds, datasets and dependency/build logs remain
under the ignored local directory `temp/performance_audit_20261022/`. The benchmark
runner validates fresh launch completion and the full planned grid before
producing summaries. The retained `profiles.txt` files are excerpts, so re-run
the runner's `--summarize` against the original campaign directory or a fresh run,
not this compact evidence directory.

Final targeted native checks:

```sh
CTOOLS_SANITIZERS=address,undefined python3 validation/test_newcommands_native.py
CTOOLS_SANITIZERS=address,undefined python3 validation/test_split_tokens_native.py
CTOOLS_LIBOMP_PREFIX="$PWD/temp/performance_audit_20261022/libomp" \
  python3 validation/test_winsor_transport_native.py
python3 validation/test_performance_audit.py
```

The winsor transport harness uses UndefinedBehaviorSanitizer and allocation
tracking, including 137 injected failures. Other native checks above use both
AddressSanitizer and UndefinedBehaviorSanitizer. This evidence is scoped to the
changed paths and is not an all-platform release certification.
