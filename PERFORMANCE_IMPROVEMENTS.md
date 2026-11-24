# ctools Runtime Performance Audit & Prioritized Optimization Roadmap

**Targets:** `cppmlhdfe`, `creghdfe`, `civreghdfe`, `csort`, `cmerge`, and shared SPI I/O / HDFE / SIMD / `strL` infrastructure\
**Date:** 2026-11-24\
**Measured Platform:** Windows x86-64, 32 logical processors, StataNow 19.5 MP (4-core Stata license), `build/ctools_windows.plugin` (OpenMP enabled, default 16 plugin threads, compiled with `-march=x86-64 -mtune=generic -ffp-contract=off`)

---

## 1. Executive Summary

End-to-end benchmarking at **N = 2,000,000** and **N = 10,000,000** observations combined with full C/ado source auditing and PE binary inspection reveals large, concrete runtime wins across all four target command families:

1. **Direct Answer on SPI I/O (Hypothesis A):** The core row-staged SPI transport ([ctools_rowstage.c](file:///G:/GitHub/stata-ctools/src/ctools_rowstage.c)) and parallel column loader ([ctools_data_io.c](file:///G:/GitHub/stata-ctools/src/ctools_data_io.c)) are **genuinely well-engineered** (unchecked `(_stata_)->vdata` loads, row-tiled parallel `safestore` writebacks). However, **SPI I/O is NOT the primary limiter in several major paths, and multiple call sites fail to use the fast transport**:
   - `cppmlhdfe` executes **40,000,000 serial `SF_vstore()` calls** at N = 10M (30M of which write temporary variables the user never asked for!) plus allocates 180 MB of unused Stata tempvars on every invocation.
   - `cmerge ..., keep(3)` on an already-sorted master loads all 21 master columns into `ctools_rowstage` (8.4× slower load) purely for rollback even when zero master rows are modified, then runs Stata's `keep if inlist(...)` even when zero rows are dropped (+0.55 s at N = 10M, making it **2.8× slower**).
   - `cmerge` writes `keepusing` columns in parallel **across columns on the same row ranges**, causing cache-line false sharing on Stata's row-major table, then makes a separate pass for `_merge` and a serial pass for using-only keys.
   - `csort`'s streaming path ([csort_stream.c](file:///G:/GitHub/stata-ctools/src/csort/csort_stream.c)) enables `SD_SAFEMODE` (checked `SF_vdata`), performs unprefetched indirect scatters, stores keys **single-threaded**, and triggers Stata `preserve` — making it **4.3× slower** than row-staged mode at N = 2M.
   - `csort`'s row-staged path incurs **150–270 ms** of unaccounted OS allocation/deallocation (`VirtualFree`/TLB shootdown) and **up to 530 ms** of cold first-touch page faults on multi-GB buffers at N = 10M.

2. **Direct Answer on `strL` Handling (Hypothesis B):** Because SPI provides no `SF_strlstore`, `csort` unconditionally falls back to native `sort, stable` whenever *any* `strL` exists anywhere in the dataset (3.4× slower at N = 2M), and `cmerge`'s `tempdta` path writes the entire merged dataset to disk and runs Stata `use`:
   - For **master `strL`** (N = 2M), `cmerge, tempdta` beats native `merge` (**1.86 s vs. 2.37 s**), but spends **0.614 s** in single-threaded `cmerge_read_strl` (2–3 SPI calls/obs) and **0.893 s** in Stata `use` + metadata restore.
   - For **using-only `strL`** (N = 2M, 10k-row lookup), `cmerge, tempdta` is **37% slower than native `merge`** (**0.729 s vs. 0.522 s**) because it dumps all 22 master columns to disk just to attach 1 using `strL` column.
   - Concrete fixes exist for both: parallel `strL` staging + direct C `.dta` reading (bypassing Stata `use` for the using file), hybrid SPI + Mata `st_sstore` for narrow using-only `strL`s, and `ctools_dta_write` support for `csort` when `strL` columns are present.

3. **Biggest Cross-Cutting & Estimation Bottlenecks:**
   - **Software `fma()` on all x86-64 builds:** `scripts/build.mk` compiles with `-march=x86-64 -mtune=generic -ffp-contract=off` (SSE2 only, no `-mavx2`/`-mfma`). Every `two_prod()` / `kahan_dot()` / `CTOOLS_MEAT_ADD` calls software-emulated `fma()` (~30–80 cycles/elem instead of 0.5 cycles hardware FMA). At N = 10M (`K = 9`), **`creghdfe, vce(robust)` takes 0.631 s vs. 0.224 s for `vce(cluster firm)` (2.8× slower for simple HC1 than clustered VCE!)**.
   - **Dead AVX2 SIMD code:** `CTOOLS_HAS_AVX2` is `0` in all x86-64 distribution binaries; weighted dot products and fused CG kernels fall back to un-unrolled single-accumulator scalar loops.
   - **`cppmlhdfe` is 17× slower in C than `creghdfe` at N = 10M (18.26 s vs. 1.06 s) and achieves only 1.31× speedup on 16 cores:** Dominated by serial transcendental loops (`ppml_compute_deviance` + double `ppml_compute_loglik` + 3× `exp()`/iter), 40M serial `SF_vstore()` calls, software `fma()` in normal equations/VCE, and 5.31 s of redundant pre-IRLS separation/collinearity CG solves when 0 rows are separated. **Projected speedup: 4×–8×.**
   - **HDFE CG solver (`creghdfe` / `civreghdfe` / `cppmlhdfe`):** Parallelizes only across columns (`min(threads, K)`), leaving 13 of 16 cores idle for `K = 3` (only **1.34× speedup** from 1 to 16 threads at N = 10M!), and uses 6-pass Kaczmarz over N-sized DRAM vectors even for 2-way FE (`G = 2`) where a **dense Schur complement** (small `L_min`) or **bipartite edge-list CG** + **algebraic Kaczmarz fusion** cuts DRAM traffic by 2×–10×.
   - **Single-threaded `qsort()` in FE/cluster remapping:** Wide-range FE IDs (e.g., 9-digit firm IDs) trigger `qsort()`, blowing up `FE remap + counts` by **14.1×** (`0.059 s → 0.834 s` at N = 10M, 46% of C runtime).

---

## 2. Measured Benchmarks (N = 2M and N = 10M)

All benchmarks ran on synthetic panel datasets (21–22 variables: numeric + `str8`/`str40`, 2-way FEs `firm` with 10k/50k levels, `group` with 100 levels, high-dim `firm2` with 10k/50k levels, `count_z` Poisson with ~50% zeros, `count_pos` strictly positive). Warmed-up timings are reported (`StataMP-64.exe` batch mode, default 16 plugin threads unless noted).

### 2.1 `cppmlhdfe` Baselines & Ablations

| Workload | N | Wall (s) | C Plugin (s) | Stata Setup (s) | HDFE Sub-iters | Notes |
|---|---:|---:|---:|---:|---:|---|
| `count_z x1 x2, absorb(firm group)` (default) | 10M | **19.72** | **18.26** | **1.12** | 17 | **17.2× slower in C than `creghdfe` (1.06 s)** |
| Same, `sep(fe)` (skip simplex + ReLU) | 10M | 14.45 | 12.94 | 1.14 | 17 | Simplex + ReLU cost **5.31 s (29.1% of C)** with 0 sep obs |
| Same, `sep(fe) use_exact_partial(0)` | 10M | 13.71 | 12.27 | 1.13 | 11 | Inexact partial cuts CG sub-iters 17 → 11 (−0.68 s) |
| `count_z x1 x2, absorb(firm group)` (default) | 2M | 3.90 | 3.65 | 0.18 | 17 | Warmed-up 16-thread baseline |
| Same, `sep(fe relu)` (skip simplex only) | 2M | 3.50 | 3.25 | 0.18 | 17 | Simplex alone = **0.40 s (11% of C)** |
| Same, `sep(fe)` (skip simplex + ReLU) | 2M | 2.84 | 2.55 | 0.23 | 17 | ReLU alone = **0.70 s (19% of C)** |
| Same, `sep(fe) use_exact_partial(0)` | 2M | 2.74 | 2.46 | 0.22 | 12 | Combined −33% before any code changes |
| `count_pos` (0 zeros, sep checks skipped) | 2M | 3.36 | 3.11 | 0.17 | 19 | Still **3.11 s** in C (11× slower than `creghdfe`!) |
| Thread scaling: `threads(1)` / `(4)` / `(16)` / `(32)` | 2M | 5.07 / 4.01 / 3.90 / 3.90 | 4.79 / 3.75 / 3.65 / 3.65 | ~0.19 | 17 | **Only 1.31× speedup at 16 cores; 0× gain past 4 cores** |

### 2.2 `creghdfe` & `civreghdfe` Baselines (N = 10M, CUE at N = 1M)

| Workload | N | Wall (s) | C Total (s) | HDFE Partial (s) | OLS / IV+VCE (s) | VCE / Remap (s) | Stata Overhead (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `creghdfe y x1 x2, absorb(firm group)` unadj | 10M | 1.49 | 1.06 | 0.578 | 0.226 | VCE 0.026 / Remap 0.059 | 0.42 |
| `creghdfe y x1 x2, absorb(firm group) vce(robust)` | 10M | 1.66 | 1.13 | 0.564 | 0.176 | **VCE 0.153** (6× unadj) | 0.54 |
| `creghdfe y x1 x2, absorb(firm group) vce(cluster firm)` | 10M | 1.77 | 1.18 | 0.576 | 0.228 | **VCE 0.119** (< robust!) | 0.59 |
| `creghdfe y x1 x2, absorb(firm firm2) vce(cluster firm)` | 10M | 1.98 | 1.43 | 0.813 | 0.237 | VCE 0.123 | 0.55 |
| `creghdfe y x1..x8 (K=9), absorb(firm firm2) vce(robust)` | 10M | 4.32 | 3.63 | 2.279 | 0.376 | **VCE 0.631** (2.8× cluster!) | 0.69 |
| `creghdfe y x1..x8 (K=9), absorb(firm firm2) vce(cluster firm)` | 10M | 4.15 | 3.43 | 2.377 | 0.430 | VCE 0.224 | 0.72 |
| `creghdfe y x1 x2, absorb(wide_fe group)` (`firm*1e6`) | 10M | 2.32 | 1.80 | 0.577 | 0.187 | **Remap 0.834** (**14.1×!**) | 0.52 |
| `creghdfe` `firm firm2` K=3: `t=1 / 4 / 16 / 32` | 10M | 2.57 / 1.98 / 1.91 / 1.80 | — | — | — | — | **1.34× at 16 threads** |
| `creghdfe` `firm firm2` K=9: `t=1 / 4 / 16` | 10M | 6.07 / 3.95 / 3.84 | — | — | — | — | **1.58× at 16 threads** |
| `civreghdfe y x2 (x1=z), absorb(firm group)` unadj | 10M | 2.49 | 1.70 | 0.751 | 0.356 | Miss+Sing 0.202 / Stats 0.105 | **0.80 (32%)** |
| `civreghdfe y x2 (x1=z), absorb(firm group) robust` | 10M | 2.69 | 1.94 | 0.758 | **0.619** | Miss+Sing 0.181 / Stats 0.115 | **0.76 (28%)** |
| `civreghdfe y x2 (x1=z), absorb(firm group) cluster(firm)` | 10M | 2.77 | 1.92 | 0.747 | 0.512 | Extr+Miss+Sing 0.259 | **0.85 (31%)** |
| `civreghdfe y x2 (x1=z z2), gmm2s` vs `cue` (1 iter) | 1M | 0.31 vs **0.43** | 0.21 vs **0.33** | 0.095 vs 0.083 | **0.056 vs 0.188** | CUE = 3.3× IV+VCE for 1 iter | 0.10 |

### 2.3 `csort` & `cmerge` Baselines (Including Streaming and `strL`)

| Workload | N | Wall (s) | C Total (s) | Load (s) | Sort / Join (s) | Store / Meta (s) | Stata Overhead (s) | Native Stata (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `csort x1 id` (21 vars, default) | 10M | 1.04 | 0.90 | 0.338 | 0.131 | 0.189 | 0.130 (post-sort 0.128) | 2.42 |
| `csort x1 id, nosortedby` | 10M | 0.67 | 0.67 | 0.236 | 0.083 | 0.199 | 0.002 | 2.42 |
| `csort code` (`str40`, default / `nosortedby`) | 10M | 1.52 / 1.08 | 1.07 / 1.08 | 0.227 | **0.523** | 0.218 | **0.452** / 0.001 | ~9.5 (est.) |
| `csort x1 id` wide `K=35` (cold / warm) | 10M | 1.80 / 1.18 | 1.69 / 1.03 | **0.904** / 0.374 | 0.120 / 0.093 | 0.255 / 0.285 | 0.114 / 0.148 | — |
| `csort x1 id, stream(8)` vs row-staged | 2M | **1.25** vs 0.36 | **1.06** vs 0.27 | keys 0.034 | 0.018 | scatter 0.429 + str 0.294 + keyst 0.064 | 0.089 + ~0.10 `preserve` | 0.74 |
| `csort x1 id` with 1 non-key `strL` (fallback) | 2M | **0.96** | — | — | — | — | Falls back to `sort, stable` | 0.96 |
| `cmerge m:1 firm` unsorted (21 vars) | 10M | 0.92 | 0.68 | 0.246 | 0.049 + 0.044 | 0.209 + 0.071 | 0.140 (inter 0.128) | 3.33 |
| `cmerge m:1 firm` sorted master (default) | 10M | **0.30** | 0.19 | 0.028 | 0.000 + 0.036 | 0.000 + 0.070 | 0.105 (inter 0.096) | — |
| `cmerge m:1 firm` sorted master, `keep(3)` | 10M | **0.85** | **0.43** | **0.236** (**8.4×!**) | 0.000 + 0.036 | 0.000 + 0.101 | **0.333** (post **0.246**) | — |
| `cmerge m:1 firm` sorted master, `keep(3) fast` | 10M | **0.62** | 0.23 | 0.028 | 0.000 + 0.042 | 0.000 + 0.107 | **0.385** (post **0.283**) | — |
| `cmerge 1:1 id` unsorted 2M × 2M using | 2M | 0.34 | 0.23 | u:0.015 m:0.046 | u:0.032 m:0.018 j:0.035 | 0.040 + 0.026 | 0.095 | 0.89 |
| `cmerge m:1` master `strL` (`tempdta`) | 2M | **1.86** | 0.89 | **0.614** (serial) | 0.009 + 0.009 | 0.236 (dta write) | **0.957** (post **0.893**) | 2.37 |
| `cmerge m:1` using `strL` (`tempdta`) | 2M | **0.73** | 0.39 | 0.173 | 0.009 + 0.009 | 0.177 (dta write) | 0.329 (post 0.245) | **0.52** |

---

## 3. Prioritized Master Roadmap

Items are prioritized by **(Wall-Clock Impact × Breadth of Workloads) / Engineering Complexity**:

| Priority | ID | Command(s) | Bottleneck | Measured Cost (N=10M unless noted) | Est. Savings / Speedup | Effort | Primary File(s) |
|---|---|---|---|---|---|---|---|
| **P0** | **P0-1** | `cppmlhdfe`, `creghdfe`, `civreghdfe` | Software `fma()` in `two_prod`, `kahan_dot`, `CTOOLS_MEAT_ADD` (no `-mfma` / runtime dispatch) | VCE robust 0.63 s (K=9), OLS 0.23–0.43 s, PPML normal eqs + quad_dot ~3–5 s | **5×–15×** on OLS/VCE/normal-eq kernels | Low–Med | `ctools_ols.h`, `ctools_ols.c`, `ctools_matrix.c`, `cppmlhdfe_irls.c` |
| **P0** | **P0-2** | `cppmlhdfe` | Unconditional allocation & **40M serial `SF_vstore()`** calls for unused `saved_d`, `tag`, `certificate`, `final_sample` | ~3–5 s (Stata setup 1.12 s + 40M serial SPI writes in C) | **−2.5 to −4.5 s** (15–25% of PPML wall) | **Low** | `build/cppmlhdfe.ado`, `src/cppmlhdfe/cppmlhdfe_irls.c` |
| **P0** | **P0-3** | `cppmlhdfe` | Serial `ppml_compute_deviance` (6×), double serial `ppml_compute_loglik` (20M `lgamma`+`exp`), 3× `exp()`/iter | ~4–6 s (major share of the 10.2 s non-CG IRLS time) | **−4 to −5.5 s** (8×–15× on deviance/ll/exp) | **Low** | `src/cppmlhdfe/cppmlhdfe_moments.c`, `src/cppmlhdfe/cppmlhdfe_irls.c` |
| **P0** | **P0-4** | `cppmlhdfe` | Redundant pre-IRLS CG projections: `simplex` (2× CG) + `relu` (1× CG + serial K=1 CG) + discarded STEP 8b CG | **5.31 s** (simplex+ReLU) + **~0.6 s** (STEP 8b) = **~5.9 s (32% of C)** | **−3.5 to −5.0 s** (reuse basis, fast-path ReLU, merge STEP 8b into iter 0) | Med | `cppmlhdfe_simplex.c`, `cppmlhdfe_relu.c`, `cppmlhdfe_irls.c` |
| **P0** | **P0-5** | `creghdfe`, `civreghdfe`, `cppmlhdfe` | 2-way FE (`G=2`) uses 6-pass N-vector Kaczmarz instead of Schur complement / bipartite edge-list CG + loop fusion | HDFE partial = **0.58–2.38 s** (`creghdfe`), **0.75 s** (`civreghdfe`), **2.1–3.2 s** (`cppmlhdfe`) | **2×–6×** on 2-way FE (zero CG iters for small `L_min`!) | Med–High | `src/hdfe/hdfe_solver.c` |
| **P0** | **P0-6** | `creghdfe`, `civreghdfe`, `cppmlhdfe` | CG solver caps threads at `min(T, K)` (leaves 13/16 cores idle for K=3; serial for K=1 ReLU/refit) | Only **1.34×** speedup at 16 cores for `K=3` | **2×–4×** on low-`K` HDFE solves | Med | `src/hdfe/hdfe_solver.c`, `src/creghdfe/creghdfe_regress.c` |
| **P0** | **P0-7** | `cmerge` | `keep()` triggers `hold_undo` (loads all master vars even when `skip_master_store==1`) + Stata `keep if` even on 0 drops | **+0.55 s** on sorted `m:1 keep(3)` (**2.8× slower**: 0.30 s → 0.85 s) | **−0.50 s** (**2.7× faster** on `keep()`) | **Low** | `build/cmerge.ado`, `src/cmerge/cmerge_impl.c` |
| **P1** | **P1-1** | All HDFE & IV | Dead `#ifdef __AVX2__` in `ctools_simd.h`; weighted dot & fused CG kernels use 1-accum scalar loops | 15–35% overhead across CG vector ops and `matmul_atdb` | **2×–4×** on SIMD vector kernels | Low–Med | `src/ctools_simd.h` |
| **P1** | **P1-2** | `cppmlhdfe` | Default `exact_partial=1` re-partials raw `X` from scratch every IRLS iter instead of warm-starting | +5–6 CG sub-iters (**0.68 s** at 10M; larger on high-dim/3-way FE) | **−15% to −40%** of IRLS CG time | Low | `src/cppmlhdfe/cppmlhdfe_options.c`, `cppmlhdfe_irls.c` |
| **P1** | **P1-3** | `creghdfe`, `civreghdfe`, `cppmlhdfe` | Single-threaded `qsort()` in FE remapping (`wide_fe`) & cluster ID mapping | **0.834 s** vs 0.059 s (**14.1× slower**, 46% of C plugin) | **10×–15×** on wide-range / non-integer FE & cluster IDs | **Low** | `hdfe_remap.c`, `ctools_hdfe_utils.c`, `civreghdfe_omega.c` |
| **P1** | **P1-4** | `civreghdfe` | 7+ redundant `ctools_matmul_atdb` calls (single-threaded when `K1*K2<=4`), double FE remap, unconditional compaction | IV+VCE **0.619 s** + Miss/Sing/Remap **0.212 s** (43% of C) | **−0.45 to −0.55 s** (**1.4×** faster C total) | Low–Med | `src/civreghdfe/civreghdfe_impl.c`, `civreghdfe_estimate.c`, `ctools_matrix.c` |
| **P1** | **P1-5** | `creghdfe`, `cppmlhdfe` | `ctools_vce_cluster` random-gather `long double` via counting sort; `volatile` RSS & serial loops in `creghdfe` | VCE cluster **0.12–0.22 s**; OLS/RSS/post serial loops ~0.08 s | **2×–3×** on cluster VCE & post-OLS | Low–Med | `src/ctools_matrix.c`, `src/creghdfe/creghdfe_regress.c`, `creghdfe_core.c` |
| **P1** | **P1-6** | `csort` | String MSD radix sort (`str17+`) runs pass 0 and large buckets single-threaded with 2-level pointer chasing | `Sort = 0.523 s` on `str40` (**6.3× slower** than numeric sort) | **3×–5×** on `str17+` string sorts | Med | `src/ctools_sort_radix_msd.c` |
| **P1** | **P1-7** | `csort` | Row-staged `csort` has **0.15–0.27 s** unaccounted alloc/free overhead (+0.53 s cold page-faults) + strided key extraction | **0.15–0.27 s** (22–27% of C total at N=10M; +0.53 s cold at K=35) | **−0.15 to −0.35 s** (1.25×–1.5× faster `csort`) | Med | `src/csort/csort_impl.c`, `src/ctools_rowstage.c` |
| **P1** | **P1-8** | `csort` | `csort_stream.c` enables `SD_SAFEMODE`, unprefetched indirect scatter, serial key store, + Stata `preserve` | **4.3× slower** than row-staged (1.25 s vs 0.36 s at N=2M) | **2.5×–3.5×** on streaming mode (wide/huge datasets) | Low–Med | `src/csort/csort_stream.c`, `build/csort.ado` |
| **P1** | **P1-9** | `cmerge` | Multi-pass Stata stores & cross-column false sharing in `keepusing` + separate `_merge` + serial using-only keys | `Write metadata = 0.071–0.107 s` (37% of sorted C total) + extra DRAM pass | **1.5×–2.5×** on store/writeback phase | Med | `src/cmerge/cmerge_impl.c`, `src/cmerge/cmerge_io.c` |
| **P1** | **P1-10** | `cmerge` | `cmerge_load_using` hardcodes `O(N log N)` `ips4o`, double-loads using, eagerly permutes `keepusing` (few threads) | `Phase 1 = 0.050 s` at 2M (~0.28 s at 10M 1:1) | **2×–3×** on large-using `1:1` / `1:m` merges | Low | `src/cmerge/cmerge_impl.c` |
| **P1** | **P1-11** | `cmerge` | `cmerge_plan`, `cmerge_sorted_join`, `orig_rows`, `seen`/`head` partition are 100% single-threaded | `Merge join = 0.036–0.044 s` at 10M (19% of sorted C total) | **3×–5×** on join planning | Med | `src/cmerge/cmerge_join.c`, `src/cmerge/cmerge_impl.c` |
| **P2** | **P2-1** | `cmerge`, `csort` | `strL` acceleration: parallel `cmerge_execute_dta`, direct C `.dta` using reader, hybrid Mata `strL`, `csort` `tempdta` | Master `strL` 1.86 s (`Load 0.61s` serial + `use 0.89s`); Using `strL` 0.73 s vs native 0.52 s; `csort strL` fallback 0.96 s | **1.5×–3×** on `strL` workloads; beats native on using `strL` | Med–High | `cmerge_impl.c`, `cmerge.ado`, `csort_impl.c`, `ctools_dta.c` |
| **P2** | **P2-2** | `creghdfe`, `civreghdfe`, `cppmlhdfe` | Stata `.ado` overhead (`Inter-plugin`, `count if`, `marksample`, `fvrevar`, `ereturn post, esample()`) | `civreghdfe` 0.76–0.85 s (31%), `creghdfe` 0.42–0.72 s (30%), `cmerge` inter-plugin 0.096–0.128 s (32%) | **−0.15 to −0.40 s** on Stata wrapper overhead | Low–Med | `cmerge.ado`, `creghdfe.ado`, `civreghdfe.ado`, `cppmlhdfe.ado` |
| **P2** | **P2-3** | `civreghdfe` | CUE `cue_objective` recomputes `Z'We` in `O(N*K_iv)` per eval; Hessian ignores symmetry (`K^2` vs `K(K+1)/2`) | CUE = 3.3× slower than GMM2S at 1M for 1 iter (`0.188 s` vs `0.056 s`) | **3×–8×** on CUE estimation | Low | `src/civreghdfe/civreghdfe_estimate.c` |
| **P2** | **P2-4** | `creghdfe`, `civreghdfe`, `cppmlhdfe` | Serial standardization, serial FE peeling (`ppml_select_fe_sample_ex`), redundant `hdfe_build_factors` remap on 0 singletons, serial `ctools_remove_singletons` recount | ~0.10–0.40 s across setup phases | **3×–8×** on setup/factor-build phases | Low | `hdfe_sample.c`, `ctools_hdfe_utils.c`, `cppmlhdfe_separation.c`, `cppmlhdfe_irls.c` |
| **P3** | **P3-1** | `csort` | Post-plugin `sort varlist, stable` solely to set Stata's `sortedby` marker | **0.128 s** (`x1 id`) to **0.451 s** (`str40`, 30% of wall!) at N=10M | Document/promote `nosortedby`; avoid in internal callers | Low | `build/csort.ado` |

---

## 4. Detailed Findings & Recommended Fixes

---

### Section A: Cross-Cutting Build, SIMD & Numerical Infrastructure

#### [P0-1] Software-Emulated `fma()` in Compensated Dot Products & Sandwich VCE
- **Locations:**
  - [scripts/build.mk](file:///G:/GitHub/stata-ctools/scripts/build.mk) (`-march=x86-64 -mtune=generic -ffp-contract=off` on Windows/Linux/macOS Intel)
  - [src/ctools_ols.h:76-80](file:///G:/GitHub/stata-ctools/src/ctools_ols.h#L76-L80) (`two_prod` calls `fma(a, b, -p)`)
  - [src/ctools_ols.c:149](file:///G:/GitHub/stata-ctools/src/ctools_ols.c#L149) (`kahan_dot`, called for every `(j,k)` pair in `compute_xtx_xty`)
  - [src/ctools_matrix.c:314-321](file:///G:/GitHub/stata-ctools/src/ctools_matrix.c#L314-L321) (`CTOOLS_MEAT_ADD` in `ctools_vce_robust`)
  - [src/cppmlhdfe/cppmlhdfe_irls.c:160-220](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L160-L220) (`dd_dot_weighted_multi` / `ppml_normal_equations`) and [src/cppmlhdfe/cppmlhdfe_solvers.c](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_solvers.c) (`ppml_quad_dot` in CG during simplex/ReLU)
- **Root Cause:** Without `-mfma`, `fma(a, b, -p)` compiles to a static libm software function call (`libmingwex.a` on Windows, `libm` on Linux/macOS) executing ~30–80 cycles per element (versus 0.5 cycles throughput for hardware `vfmadd231pd` processing 4 doubles simultaneously!).
- **Measured Evidence:**
  - At N = 10M, `K = 9`: `creghdfe ..., vce(robust)` VCE takes **0.631 s** while `vce(cluster firm)` (which does a counting sort + random gather in `long double`, *without* `fma()`) takes only **0.224 s** — simple HC1 is **2.8× slower** than clustered VCE!
  - `creghdfe` `OLS estimation` for just 3 columns (240 MB) takes **0.226 s** (vs. ~15 ms memory bandwidth).
  - In `cppmlhdfe`, `dd_dot_weighted_multi` executes on every IRLS iteration (`K*(K+1)/2` pairs) and `ppml_quad_dot` executes 3–4× per CG iteration in `ppml_reference_partial`.
- **Recommended Fix:**
  1. Add **runtime CPUID multiversioned kernels** using `__attribute__((target("avx2,fma")))` guarded by `__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")` (the exact pattern already used for AVX-512 in [csort_stream.c:101](file:///G:/GitHub/stata-ctools/src/csort/csort_stream.c#L101)).
  2. Vectorize compensated `two_prod` + `two_sum` across 4 lanes (`__m256d`) with `_mm256_fmsub_pd(a, b, p)` (exact hardware FMA error term in 1 instruction!) and 2–4 independent vector accumulator lanes to hide latency.
  3. For pre-Haswell x86-64 CPUs lacking FMA, replace software `fma()` in `two_prod` with pure-FP **Dekker/Veltkamp splitting** (`c = 134217729.0 * a; a_hi = c - (c - a); a_lo = a - a_hi; ...`) guarded against overflow (which is ~10× faster than `libm`'s `fma()` and vectorizes with SSE2).
  4. Parallelize `ctools_vce_robust` and `ppml_normal_equations` across **row blocks** when `K` is small (`K <= num_threads`), summing thread-local `(hi, lo)` double-double pairs deterministically in thread order (`0..T-1`) so results remain reproducible while utilizing all 16 cores.

#### [P1-1] Dead AVX2 Code & Un-Unrolled Scalar Fallbacks in `ctools_simd.h`
- **Location:** [src/ctools_simd.h:497-560, 799-980](file:///G:/GitHub/stata-ctools/src/ctools_simd.h#L497-L560)
- **Root Cause:** Because x86-64 builds compile with `-march=x86-64`, `__AVX2__` is never defined (`CTOOLS_HAS_AVX2 == 0`). Worse, `ctools_simd_dot_weighted`, `ctools_simd_fused_axpy2_dot`, `ctools_simd_fused_axpy2_wdot`, `ctools_simd_axpy`, and `ctools_simd_axpby` have **no SSE2 path** and fall through to single-accumulator scalar loops (e.g., `sum += w[i] * a[i] * b[i]`), serialized by a 4-cycle FP add dependency chain because `-fno-fast-math` forbids compiler reassociation.
- **Recommended Fix:** Dispatch these hot kernels through function pointers initialized once at plugin load via CPUID (`avx2+fma` > `sse2` 4-way unrolled). This immediately speeds up every CG iteration (`hdfe_solver.c`) and every `ctools_matmul_atdb` call (`civreghdfe`) by 2×–4×.

---

### Section B: `cppmlhdfe` (Highest Overall ROI — Projected 4×–8× Faster)

At N = 10M (`count_z x1 x2, absorb(firm group)`), `cppmlhdfe` takes **19.72 s** (**18.26 s** in C, **1.12 s** Stata setup) — **17.2× slower** in C than `creghdfe` on the same dataset — and achieves only **1.31× speedup on 16 cores**. Eliminating the bottlenecks below can bring this benchmark from **~19.7 s down to ~2.5–4.0 s**.

#### [P0-2] Eliminate Unconditional Allocation & 40M Serial `SF_vstore()` Calls
- **Locations:**
  - [build/cppmlhdfe.ado:874-918, 959](file:///G:/GitHub/stata-ctools/build/cppmlhdfe.ado#L874-L918)
  - [src/cppmlhdfe/cppmlhdfe_irls.c:1508-1511, 2093-2099, 2101-2106, 2181-2184](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L2091-L2185)
- **Root Cause:**
  1. `cppmlhdfe.ado` unconditionally creates **4 full-dataset temporary variables** (`final_sample`, `double saved_d`, `byte separation_tag`, `double separation_certificate` = **180 MB** at N = 10M) and re-zeros `final_sample` on line 959 even when the user never requested `d()`, `tagsep()`, or `zvarname()`.
  2. Because their variable indices are non-zero, `cppmlhdfe_irls.c` unconditionally runs:
     - **20M serial `SF_vstore()` calls** writing `tag_idx` and `certificate_idx` (lines 1508–1511).
     - **10M serial `SF_vstore()` calls** computing and writing `saved_d` (lines 2093–2099).
     - **10M serial `SF_vstore()` calls** writing `final_sample` (lines 2181–2184).
     - An unconditional 10M-row serial loop filling `separation_certificate` for `ppml_save_effects` even when no FE is being saved (lines 2101–2106).
- **Measured Evidence:** `cppmlhdfe` Stata setup is **1.12 s** at N = 10M (vs. `creghdfe`'s **0.23 s**), and 40M serial `SF_vstore()` calls add ~2–3 s of single-threaded SPI latency inside C.
- **Recommended Fix:**
  1. In `cppmlhdfe.ado`, only `generate` `saved_d` when `"`d'" != ""` (else set `__cppmlhdfe_d_idx = 0`), `separation_tag` when `"`tagsep'" != ""`, and `separation_certificate` when `"`zvarname'" != ""`. Skip `replace final_sample = 0` on the first pass (it was just generated as `0`).
  2. In `cppmlhdfe_irls.c`, guard lines 2101–2106 on whether any `save_fe_idx[g] > 0`, and replace the serial `sample_var_idx` loop (line 2182) with parallel `ctools_store_filtered_rowpar` (or parallel `(_stata_)->safestore` over disjoint chunks). When `N == N_orig` (no dropped obs), `cppmlhdfe.ado` can reuse `touse` directly as `esample(touse)` without touching a single row!

#### [P0-3] Parallelize & Deduplicate Transcendental Loops (`exp`, `log`, `lgamma`) in Deviance, IRLS, and Log-Likelihood
- **Locations:**
  - [src/cppmlhdfe/cppmlhdfe_irls.c:467-516, 577](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L467-L516)
  - [src/cppmlhdfe/cppmlhdfe_moments.c:87-139](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_moments.c#L87-L139)
  - [src/cppmlhdfe/cppmlhdfe_irls.c:2074-2082, 2164-2177](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L2074-L2178)
- **Root Cause:**
  1. **3× `exp()` per observation per IRLS iteration:** Pass 1 (`irls.c:471`) calls `exp(xb_work[i] + off)` solely to test overflow; Pass 2 (`irls.c:510`) computes `mu[i] = exp(eta[i] + off)`; Pass 3 (`moments.c:89`, called from `irls.c:577`) calls `exp(eta_total)` a **3rd time** inside `ppml_compute_deviance` — **100% single-threaded!**
  2. **Repeated `log(y[i])` in deviance:** `ppml_unit_deviance` recomputes `log(y[i])` on every IRLS iteration even though `y` never changes.
  3. **Double single-threaded `ppml_compute_loglik` at end (`irls.c:2164, 2177`):** Calls `lgamma(y[i] + 1)` + `log(y[i])` + `exp(...)` **single-threaded TWICE** (20M `lgamma` + 20M `exp` at N = 10M!). For the null model `ll_0`, `eta` is a constant `log(y_mean)`, yet it fills `fe_work[i] = log(y_mean)` and evaluates `exp(log(y_mean))` and `ppml_saturated_loglik(y[i])` 10 million times!
  4. **Serial `log(mu[i])` in STEP 12 (`irls.c:2079`):** Calls `log(mu[i])` serially on all N rows when `log(mu[i]) == eta[i] + off` (except on floored rows where `mu[i] == PPML_MU_FLOOR`).
- **Recommended Fix:**
  1. Replace the Pass 1 overflow check `isfinite(exp(x))` (`irls.c:471, 484`) with the exact threshold check `x > 709.782712893384` (zero `exp()` calls!).
  2. Pass already-computed `mu[i]` (before flooring, or compute `mu_raw` + unit deviance **fused inside Pass 2's OpenMP parallel loop** with thread-local `dd_real` reduction) — cutting `exp()` calls by **3×** and parallelizing deviance across all 16 cores!
  3. Precompute `ylog_y` (or use a fast branch/table for small integer counts `y ∈ {0..31}`) once before IRLS, and parallelize `ppml_compute_loglik` with OpenMP.
  4. Compute `ll` and `ll_0` in a **single parallel pass** (or derive `ll = sat_ll - 0.5 * deviance_unstd` directly from the already-converged `deviance`, and `ll_0 = sat_ll + sum_wy * (log(y_mean) - 1.0) - sum_w * (y_mean - ...)` with **zero** `exp()` calls and only **one** parallel evaluation of `ppml_saturated_loglik(y)` using a 10-entry lookup table for integer `y ∈ {0..9}`!).
  5. Parallelize STEP 12's `constant` / `mean_d` reduction and use `eta[i] + off` whenever `mu[i] > PPML_MU_FLOOR`.

#### [P0-4] Eliminate Redundant Pre-IRLS CG Projections (Simplex, ReLU, STEP 8b)
- **Locations:**
  - [src/cppmlhdfe/cppmlhdfe_simplex.c:169-185](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_simplex.c#L169-L185)
  - [src/cppmlhdfe/cppmlhdfe_relu.c:139-220](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_relu.c#L139-L220)
  - [src/cppmlhdfe/cppmlhdfe_irls.c:1590-1637](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L1590-L1637)
- **Root Cause:** On any dataset with `y == 0` rows (default `sep(fe simplex relu)`), `cppmlhdfe` executes **four full CG solves of `X` before IRLS iter 1 even starts**:
  1. `ppml_simplex` runs `ppml_reference_partial` (`tol=1e-9`, using software `fma()` in `ppml_quad_dot`!) AND THEN immediately calls `partial_out_columns` (`tol=1e-13`) on the same `X`, plus 2-pass MGS with `ppml_quad_dot`.
  2. `ppml_relu_ex2` runs `ppml_relu_basis` — **another** full `ppml_reference_partial` on `X` (plus `hypot()` via external `msvcrt.dll` and two `sqrt(S->weights[i])` calls per row even though weights only take 2 distinct values `{1.0, 67108864.0}`!), followed by `relu_iter` calling **single-threaded (`K=1`)** `ppml_reference_partial(S, r, N, 1)`.
  3. **STEP 8b (`irls.c:1590-1637`)** runs **another** full CG solve on `[0, X]` with initial weights `irls_w` at full `tolerance` (`1e-8`) to detect FE-collinear regressors, and then **throws the partialled `X~` away** — IRLS iter 0 immediately copies raw `X` into `aug_copy` and re-partials `X` with the **exact same `irls_w`**!
- **Measured Evidence:** At N = 10M, `simplex + relu` costs **5.31 s (29.1% of C time)** even when **0 observations are separated** (`simplex` = 11%, `relu` = 19%). STEP 8b adds another ~0.6 s.
- **Recommended Fix:**
  1. **Eliminate STEP 8b's duplicate CG solve:** Keep the partialled `X~` from STEP 8b directly for IRLS iter 0 (only partial `z_orig` in iter 0, or better yet, put `z_orig` into column 0 during STEP 8b so STEP 8b *is* IRLS iter 0!). Instant 100% elimination of STEP 8b's CG solve.
  2. **Deduplicate Simplex & ReLU basis construction:** When both `simplex` and `relu` are enabled, reuse the weighted-projection / orthonormalized `X` basis rather than projecting `X` three separate times, and replace `ppml_reference_partial`'s slow `ppml_quad_dot` with hardware-FMA SIMD dots (P0-1) + row-parallel `K=1` CG (P0-6).
  3. **Fast-path `sqrt(w)` in ReLU:** Since ReLU weights are binary (`1.0` for `y == 0` vs `K_relu` for `y > 0`), replace `sqrt(S->weights[i])` per element with a 2-value branch/lookup (`1.0` vs `inv_sqrt_K`).
  4. **Quick early-exit check in ReLU/Simplex:** Verify whether a single coarse-tolerance (`1e-4`) pass of `r` already proves `max_{y=0} r_i <= 0` (no separation certificate exists) before paying for `1e-13` MGS refinement.

#### [P1-2] Warm-Started (`use_exact_partial(0)`) Partialling During IRLS
- **Location:** [src/cppmlhdfe/cppmlhdfe_options.c:19](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_options.c#L19), [src/cppmlhdfe/cppmlhdfe_irls.c:373-375](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L373-L375)
- **Root Cause:** `options->exact_partial` defaults to `1`, overwriting `aug_copy` with raw `X` and re-partialling from scratch on every IRLS iteration. Because `M_{W_new}(X - P_{W_old} X) = M_{W_new}(X)`, warm-starting from the previous iteration's `X~` is mathematically exact up to floating-point drift (which is already guarded by a final exact pass when needed).
- **Measured Evidence:** Switching to `use_exact_partial(0)` reduces HDFE sub-iterations from **17 → 11** (−35%) at N = 10M with identical coefficients/SEs to 7+ digits, saving **0.68 s** on easy 2-way FE and **2×–3× more** on high-dim or 3-way FEs.
- **Recommended Fix:** Make warm-started partialling the default for intermediate IRLS iterations (`irls_iter > 0`), with an automatic exact refresh only on the final convergence check (or if `max_eta_change` / weight ratio jumps sharply).

---

### Section C: `creghdfe`, `civreghdfe`, and Shared HDFE Solver (`src/hdfe/`)

#### [P0-5] 2-Way FE (`G = 2`) Schur Complement Solver, Bipartite Edge-List CG, and Kaczmarz Loop Fusion
- **Location:** [src/hdfe/hdfe_solver.c:160-425](file:///G:/GitHub/stata-ctools/src/hdfe/hdfe_solver.c#L160-L425)
- **Root Cause:**
  1. **6-Pass DRAM-Bound Kaczmarz on `G = 2`:** Every CG iteration of `transform_sym_kaczmarz_threaded` makes **6 full passes over `N`-sized arrays** (`y`, `ans`, `f0.levels`, `f1.levels` = ~480 MB/column/iter at N = 10M), plus 4 more passes for `dot(u,v)`, `fused_axpy2_dot`, and `axpby`.
  2. **Wasted `axpby` on the final CG iteration (`hdfe_solver.c:377, 415`):** Computes `u = r + beta * u` (240 MB/col at N = 10M) **before** checking convergence (`recent_ssr <= tol_sq * improvement_potential`), discarding `u` on the final iteration (25% of all iters in `creghdfe`, 50% in early `cppmlhdfe` IRLS!).
  3. **No parameter-space solver for `G = 2`:** Two-way FE (`firm + year`, `firm + group`, `worker + firm`) is the overwhelmingly dominant empirical specification, yet ctools iterates over `N` observations instead of `L_min = min(L0, L1)` parameters.
- **Recommended Fixes (tiered from quick win to algorithmic breakthrough):**
  1. **Quick Win 1 (3 lines):** Move the convergence check **before** `ctools_simd_axpby(u, r, 1.0, beta, N)` in `cg_solve_column_threaded` (lines 377–387 and 415–425).
  2. **Quick Win 2 (`G = 2` Algebraic Kaczmarz + `uv` Fusion):** For `G = 2` (no slopes), use the exact identity:
     `a = P0(u)`, `b = P1(u - a)`, `c = P0(b)`, `v = (a - c) + b`.
     - Pass 1: scatter `u` → `m0` (read `u, l0`; **no write to `v`**)
     - Pass 2: scatter `u - m0[l0]` → `m1` (read `u, l0, l1`; **no write to `v`**)
     - Pass 3: scatter `m1[l1]` → `m0_b` (read **only `l0, l1`** — 8 bytes/obs, zero reads of `u` or `v`!)
     - Pass 4 (fused with `uv = dot(u, v)`): `v[i] = (m0[l0] - m0_b[l0]) + m1[l1]; uv += w[i] * u[i] * v[i];`
     - Cuts `G = 2` Kaczmarz + dot DRAM traffic by **~50%**.
  3. **Algorithmic Breakthrough A — Dense Schur Complement for `L_min <= 2000` (e.g., `firm + group`, `firm + year`):**
     - Eliminating the larger FE (`L0`, e.g., 50k firms) gives an `L_min × L_min` symmetric positive semidefinite normal system:
       `(W1 - B' W0^{-1} B) α1 = rhs1`, where `B = D0' W D1` is the sparse cross-tab.
     - Build `S = W1 - B' W0^{-1} B` **once** (shared across all `K` columns!), Cholesky-factor `S` (pinning 1 redundant level per connected component) in **<1 ms** for `L_min = 100`, collect `rhs` for all `K` columns in **1 pass**, solve for `(α0, α1)`, and subtract `α0[l0] + α1[l1]` in **1 pass**.
     - **Result: ZERO CG iterations!** For `firm + group` at N = 10M (`K = 3`), HDFE partial-out drops from **0.58 s → ~0.08 s** (~7× faster); for `K = 9`, from **~1.5 s → ~0.15 s** (~10× faster)!
  4. **Algorithmic Breakthrough B — Bipartite Edge-List Parameter-Space CG for large `L_min > 2000` (e.g., `firm + firm2`, AKM panels):**
     - Deduplicate `(l0, l1)` pairs into a compact CSR edge list of size `E <= N` (with summed edge weights `w_e`), and run CG directly on `(W1 - B' W0^{-1} B) α1 = rhs1` in `L_min`-sized vectors (**400 KB** in L2 cache instead of **320 MB** of N-vectors in DRAM!), updating all `K` columns simultaneously so each CSR edge load is amortized across `K` columns.

#### [P0-6] Row-Parallel & Hybrid 2D Parallelism in HDFE CG Solver
- **Locations:** [src/hdfe/hdfe_solver.c:465-495](file:///G:/GitHub/stata-ctools/src/hdfe/hdfe_solver.c#L465-L495), [src/creghdfe/creghdfe_regress.c:308](file:///G:/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c#L308), [src/ctools_hdfe_utils.c:718](file:///G:/GitHub/stata-ctools/src/ctools_hdfe_utils.c#L718)
- **Root Cause:** `num_threads = min(max_threads, K)`. For typical specifications (`K = 3`: `y x1 x2`), **13 of 16 threads sit idle** (and 29 of 32!). For `K = 1` (`ppml_relu`, `ppml_refit_fixed_effects`), the CG solver runs **100% single-threaded**.
- **Measured Evidence:** At N = 10M on `creghdfe y x1 x2, absorb(firm firm2)`: `threads(1) = 2.57 s`, `threads(4) = 1.98 s`, `threads(16) = 1.91 s` (**only 1.34× speedup on 16 cores!**).
- **Recommended Fix:**
  - When `K < max_threads`, assign `T_col = max_threads / K` threads per column (or run columns sequentially / in small groups with all `T` threads over rows `N`).
  - For each factor `g`: if `num_levels[g] * sizeof(double) <= 2 MB` (fits in L2/L3 cache, e.g., up to 250k levels), give each thread a private `means[num_levels]` buffer, scatter over `N / T_col` rows in parallel, reduce the `T_col` buffers in parallel, and gather/subtract in parallel. For huge `num_levels > 250k`, partition by `levels[i] % T_col` or sort/bucket observations by factor-0 level once during setup (which also makes factor-0 scatter sequential!).
  - Parallelize all CG vector ops (`dot`, `fused_axpy2_dot`, `axpby`) over row blocks of `N`.

#### [P1-3] Replace Single-Threaded `qsort()` in FE & Cluster Remapping with Parallel Radix Sort
- **Locations:**
  - [src/hdfe/hdfe_remap.c:195](file:///G:/GitHub/stata-ctools/src/hdfe/hdfe_remap.c#L195) (`remap_sort_impl`)
  - [src/ctools_hdfe_utils.c:235, 270, 349, 908](file:///G:/GitHub/stata-ctools/src/ctools_hdfe_utils.c#L235) (`ctools_strings_to_cluster_ids`, `ctools_pair_ids`, `ctools_numeric_to_cluster_ids`)
  - [src/civreghdfe/civreghdfe_omega.c:124](file:///G:/GitHub/stata-ctools/src/civreghdfe/civreghdfe_omega.c#L124)
- **Root Cause:** Whenever FE or cluster IDs have a range `> 4 * N` (extremely common in practice: 9-digit EINs, CUSIPs, ZIP+4, timestamps, or non-integer IDs), remapping falls back to C stdlib's single-threaded `qsort()` (`msvcrt.dll` on Windows).
- **Measured Evidence:** At N = 10M, changing `absorb(firm group)` to `absorb(wide_fe group)` (`firm * 1e6`) increases `FE remap + counts` from **0.059 s → 0.834 s (14.1× slower — 46% of the entire C plugin time!)**.
- **Recommended Fix:** Replace `qsort()` in `remap_sort_impl` and `ctools_numeric_to_cluster_ids` with either:
  1. A **parallel open-addressing hash table** (or partitioned thread-local hash tables -> merge) when `N` is large (O(N), zero sorting, cache-friendly when `num_levels << N`), OR
  2. ctools' existing **parallel LSD radix sort** on `(ctools_double_to_sortable(val), idx)` pairs (10×–15× faster than `qsort()` even when `num_levels ≈ N`).
  - Also parallelize the min/max scan + counting remap pass inside `remap_and_count` across row chunks!

#### [P1-4] Eliminate `civreghdfe`'s 7× Redundant Matrix Multiplies, Single-Threaded `matmul_atdb`, Double FE Remap & Unconditional Compaction
- **Locations:**
  - [src/ctools_matrix.c:38, 56, 169, 180](file:///G:/GitHub/stata-ctools/src/ctools_matrix.c#L169-L190) (`ctools_matmul_atdb` / `atb` / `ata`)
  - [src/civreghdfe/civreghdfe_impl.c:784-804, 875-920, 1518, 1628, 1691-1693](file:///G:/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c#L1518-L1694)
  - [src/civreghdfe/civreghdfe_estimate.c](file:///G:/GitHub/stata-ctools/src/civreghdfe/civreghdfe_estimate.c) and [src/civreghdfe/civreghdfe_tests.c](file:///G:/GitHub/stata-ctools/src/civreghdfe/civreghdfe_tests.c)
- **Root Cause:**
  1. **`ctools_matmul_atdb` is 100% single-threaded for `K1 * K2 <= 4`!** (`if((size_t)K1 * K2 > 4)`). For `y x2 (x1 = z)` (`K_total = 2, K_iv = 2`), every Gram matrix runs on **1 thread** using the un-unrolled scalar `ctools_simd_dot_weighted` (P1-1), AND symmetric `A'WA` computes both triangles (2× wasted work).
  2. **Redundant Gram computations:** `civreghdfe_impl.c` computes `X'WX` (line 1518, freed), `Z'WZ` (line 1628, freed), `Z'WZ` **again** + `Z'WX_endog` + `X_endog'WX_endog` (lines 1691–1693, freed); then `ivest_compute_2sls` recomputes `Z'WZ`, `Z'WX`, `Z'Wy` a **3rd time**; then `civreghdfe_tests.c` recomputes first-stage cross-products **again**!
  3. **Double FE remap & unconditional copy:** Lines 784–802 serially compact & `remap_and_count` each FE, **free `counts_g`**, and later remap all FEs a second time. Lines 876–920 allocate and gather-copy every array even when `N == N_total` (0 missing values, 0 singletons).
- **Measured Evidence:** At N = 10M (`robust`), `IV estimation + VCE = 0.619 s` (32% of C) and `Missing + Singleton + Remap = 0.212 s` (vs. `creghdfe`'s `0.077 s`).
- **Recommended Fix:**
  1. Compute the single joint weighted Gram matrix `[y, X_exog, X_endog, Z_excl]' W [y, X_exog, X_endog, Z_excl]` **once** in a row-blocked parallel AVX2+FMA kernel (exploiting symmetry!), and extract all sub-blocks (`X'WX`, `Z'WZ`, `Z'WX`, `Z'Wy`, first-stage `Z'WX_endog`) from that tiny `K × K` matrix in microseconds!
  2. Parallelize `ctools_matmul_atdb` / `atb` / `ata` across **row blocks** whenever `(size_t)N * K1 * K2 > 100000` regardless of how small `K1 * K2` is.
  3. Reuse `fe_levels_c` + `counts_g` directly in `HDFE_State` and skip `civ_gather` allocation/copying when `N == N_total` (just alias `filtered.data.vars[...].data.dbl` directly).

#### [P1-5] `ctools_vce_cluster` Cache-Friendly Scatter-Add & `creghdfe` Serial Loops / `volatile` RSS
- **Locations:**
  - [src/ctools_matrix.c:395-507](file:///G:/GitHub/stata-ctools/src/ctools_matrix.c#L395-L507) (`ctools_vce_cluster`)
  - [src/creghdfe/creghdfe_regress.c:1314-1325](file:///G:/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c#L1314-L1325) and [src/creghdfe/creghdfe_core.c:231, 250](file:///G:/GitHub/stata-ctools/src/creghdfe/creghdfe_core.c#L231) (`volatile ST_double sq`)
- **Root Cause:**
  1. `ctools_vce_cluster` performs a serial counting sort (`vce_sort_by_cluster`) to build `sort_perm[N]`, then gathers `e_w[sort_perm[pos]] * xk[sort_perm[pos]]` randomly in unvectorizable x87 `long double` (parallelized only over `K`).
  2. `creghdfe` marks `volatile ST_double sq` inside the N-row RSS loop (forcing memory spills every observation and blocking vectorization!), and runs Step 6 (TSS/stdevs), Step 7 (standardize/`norm2`), and `sample_var_idx`/`resid`/`savefe` writeback serially.
- **Recommended Fix:**
  1. In `ctools_vce_cluster`, when `num_clusters * K * sizeof(double) <= 16 MB` (e.g., 50k clusters × 10 vars = 4 MB — covers almost all workloads!), replace `vce_sort_by_cluster` + random gather with a **sequential pass** over `i = 0..N-1` scattering `e_w[i] * xk[i]` into `cluster_sums[cluster_ids[i] * K + k]` (parallelized over row chunks with thread-local `cluster_sums` when `num_clusters * K <= 256k`, or over columns `k` otherwise), using double-double `(hi, lo)` or compensated accumulation instead of x87 `long double`.
  2. Remove `volatile` from `sq` in `creghdfe_regress.c` and `creghdfe_core.c` (use `-ffp-contract=off` / `two_prod` if reproducibility is desired, not `volatile`!), parallelize Steps 6–7 & RSS, and use `ctools_store_filtered_rowpar` for `sample_var_idx`/`resid`/`savefe`.

#### [P2-3] `civreghdfe` CUE Estimation: `O(K_iv * K)` Moments & Symmetric Hessian
- **Location:** [src/civreghdfe/civreghdfe_estimate.c](file:///G:/GitHub/stata-ctools/src/civreghdfe/civreghdfe_estimate.c) (`ivest_compute_cue`, `cue_objective`)
- **Root Cause:** Every call to `cue_objective(beta)` recomputes `e = y - X*beta` and `g = Z'We` in `O(N * K_iv)` time, and finite-difference Hessian evaluation computes all `K^2` entries instead of the `K*(K+1)/2` upper triangle.
- **Measured Evidence:** At N = 1M, even 1 CUE Newton iteration takes **0.188 s** in `IV estimation + VCE` vs. **0.056 s** for GMM2S (3.3× slower).
- **Recommended Fix:** Precompute `Z'Wy` and `Z'WX` once so `g(beta) = Z'Wy - (Z'WX)*beta` takes `O(K_iv * K)` (~nanoseconds!); exploit symmetry `H[j,k] == H[k,j]` (2× fewer objective evaluations); and for `vce(robust)` CUE, precompute the `O(K_iv^2 * (K+1)^2)` fourth-moment tensor when `K` and `K_iv` are small (e.g., `K <= 10, K_iv <= 10`) so the **entire CUE Newton optimization runs in `O(1)` in `N`** (zero passes over `N` inside the optimizer!).

#### [P2-4] Setup & Factor-Build Micro-Optimizations Across HDFE Commands
- **Locations:** [src/hdfe/hdfe_sample.c:65-132](file:///G:/GitHub/stata-ctools/src/hdfe/hdfe_sample.c#L65-L132), [src/ctools_hdfe_utils.c:51-120](file:///G:/GitHub/stata-ctools/src/ctools_hdfe_utils.c#L51-L120), [src/cppmlhdfe/cppmlhdfe_separation.c:38-66](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_separation.c#L38-L66), [src/cppmlhdfe/cppmlhdfe_moments.c:36-81](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_moments.c#L36-L81), [src/cppmlhdfe/cppmlhdfe_irls.c:1413-1473](file:///G:/GitHub/stata-ctools/src/cppmlhdfe/cppmlhdfe_irls.c#L1413-L1473)
- **Fixes:**
  1. In `hdfe_build_factors`, when `num_singletons == 0` (`N == N_in`), skip re-allocating, re-remapping, and re-counting all `G` factors — move ownership of `levels` and `counts` directly!
  2. Pass already-known `num_levels` and `counts` into `ctools_remove_singletons` to avoid re-scanning `N` twice per factor.
  3. Replace iterative full-scan peeling in `ppml_select_fe_sample_ex` with queue-based peeling (as in `ctools_remove_singletons`).
  4. Parallelize `cppmlhdfe` standardization (`irls.c:1413-1473`) across rows and collapse its 5 serial passes into 2 parallel passes; parallelize `ppml_update_fe_weights` over rows (not just `G = 2` threads) and hoist the `f->slope` check outside the `N`-row loop.

---

### Section D: `csort` & SPI I/O Verification

#### SPI I/O Audit Verdict
At N = 10M (21 variables, ~1.6 GB in Stata), `csort x1 id, nosortedby` achieves:
- `Data load`: **0.236 s** (~6.8 GB/s through unchecked `(_stata_)->vdata`)
- `Sort`: **0.083 s** (parallel LSD/counting radix sort)
- `Data store`: **0.199 s** (~8.0 GB/s through row-tiled `(_stata_)->safestore`)

The core `ctools_rowstage` SPI loops are close to hardware/SPI throughput limits. **However, four major non-SPI and suboptimal-SPI bottlenecks still dominate specific `csort` regimes:**

#### [P1-6] String MSD Radix Sort: Single-Threaded Top Pass & 2-Level Pointer Chasing
- **Location:** [src/ctools_sort_radix_msd.c:803-969](file:///G:/GitHub/stata-ctools/src/ctools_sort_radix_msd.c#L803-L969)
- **Root Cause:** For any string wider than 16 bytes (`str17+`), `auto_select_algorithm` picks `SORT_ALG_MSD`. Unlike numeric MSD, `msd_sort_by_string_var` has **no parallel top-level pass**: it enters `#pragma omp single` and runs `msd_sort_pass_string_seq` **100% single-threaded** on all N rows! When the first character has low cardinality (e.g., `"East"/"North"/"South"/"West"` in `code`), the first 5–8 passes run on only 1–4 threads while chasing 2 pointer indirections (`strings[order[i]][char_pos]`) with no prefetching.
- **Measured Evidence:** At N = 10M, `csort code` spends **0.523 s** in `Sort` (**6.3× slower** than numeric `x1 id` at 0.083 s, and 49% of total C runtime!).
- **Recommended Fix:**
  1. Add a **parallel histogram + scatter pass** (`msd_sort_pass_string_parallel`, using OpenMP — not `pthread_create`!) whenever `end - start >= 65536`.
  2. **8-byte prefix packing (huge cache win):** Pack `(uint64_t prefix_be, perm_idx_t row_idx)` (12–16 bytes, aligned) in parallel during load/key-extraction so the first 8 characters sort via a **zero-indirection parallel 11-bit or 8-bit radix sort** (just like numeric keys!), and only fall back to pointer dereferencing (`char_pos >= 8`) inside equal-prefix buckets! Even better, for fixed-width Stata strings (`str#`) staged in `ctools_rowstage`, read the 8-byte prefix directly from the flat row buffer without ever allocating `char *` pointers or calling `strlen()`!

#### [P1-7] Eliminate `csort_rowstaged`'s 0.15–0.27 s Unaccounted Alloc/Free Gap & Strided Key Extraction
- **Locations:** [src/csort/csort_impl.c:394-433](file:///G:/GitHub/stata-ctools/src/csort/csort_impl.c#L394-L433), [src/ctools_rowstage.c](file:///G:/GitHub/stata-ctools/src/ctools_rowstage.c)
- **Root Cause:**
  1. **Measured 0.15–0.27 s unaccounted C gap:** At N = 10M (`x1 id, nosortedby`), `Load (0.236) + Sort (0.083) + Store (0.199) = 0.518 s`, yet `C plugin total = 0.666 s` (**+0.148 s hidden**, rising to **+0.273 s** at K = 35!). On cold allocation (`K = 35`, 2.8 GB `staged.rows`), `Load` spikes from `0.374 s → 0.904 s` (**+0.530 s** of Windows OS zero-fill page faults), and `ctools_rowstage_free` (outside any sub-timer) pays **150–270 ms** of `VirtualFree` page-table teardown & TLB shootdowns!
  2. **Strided post-load key extraction & serial `sort_order` init:** `csort_rowstaged` initializes `keys.sort_order[i] = i` serially (line 407), then calls `ctools_rowstage_column` (line 410), which re-reads the entire 1.6–2.8 GB `staged.rows` buffer with huge row strides just to copy the 1–2 key columns!
- **Recommended Fix:**
  1. **Fuse key extraction into `ctools_rowstage_load`:** Pass an optional `extract_keys` list into `ctools_rowstage_load` so key doubles (and 8-byte string prefixes!) are written to contiguous key arrays *while each row tile is hot in L1 cache*, eliminating a full 1.6–2.8 GB strided pass over `staged.rows`. Parallelize `keys.sort_order` initialization (or let the first radix pass treat `sort_order == NULL` as implicit identity `0..N-1`, saving 40 MB of writes!).
  2. **Retain/pool the large staging buffer** (optional bounded cache, e.g., freed on plugin unload or after idle/OOM, or `VirtualAlloc` + `MEM_RESET` on Windows / `madvise(MADV_DONTNEED)` on Linux only when shrinking) to eliminate the 150–270 ms `VirtualFree` and 530 ms cold page-fault tax on back-to-back sorts/merges.
  3. Report `ctools_rowstage_free` inside `timer->cleanup_time` so verbose diagnostics never hide deallocation cost.

#### [P1-8] Fix `csort_stream.c` Performance Regressions (4.3× Slower Than Row-Staged)
- **Locations:** [src/csort/csort_stream.c:49, 56, 200-350, 688](file:///G:/GitHub/stata-ctools/src/csort/csort_stream.c#L56), [build/csort.ado:12](file:///G:/GitHub/stata-ctools/build/csort.ado#L12)
- **Root Cause:** When streaming mode is used (`stream()` or auto-triggered at `K > 30 & N > 10M`), it is **4.3× slower** (`1.25 s` vs `0.36 s` at N = 2M) due to five concrete issues:
  1. `#define SD_SAFEMODE` (line 56) forces checked `SF_vdata` instead of unchecked `(_stata_)->vdata` (~15–23% slower loads).
  2. `scatter_8_scalar` is called via a non-inlined function pointer every 8 doubles on non-AVX-512 CPUs, with **zero prefetching** (`CTOOLS_PREFETCH_NTA` is defined at line 49 and never used!) during random scatter across 80 MB+ buffers.
  3. Phase 4 (`csort_stream.c:688`) stores sorted keys via **single-threaded `ctools_store_filtered`** (taking 64 ms for just 2 keys at N = 2M — 63% as long as storing all 22 vars in parallel!).
  4. `csort.ado:12` executes Stata's **`preserve`** whenever streaming triggers, copying the entire multi-GB dataset to RAM/disk *before* C starts (~0.10 s at 2M, ~1.5–3 s at 10M+!).
  5. Streaming loads/stores variables in column batches instead of **gather-reading** via `perm[i]` sequentially (gather + prefetch is dramatically faster than unprefetched random scatter on DRAM!).
- **Recommended Fix:**
  1. Remove `#define SD_SAFEMODE` (use unchecked `(_stata_)->vdata` for loads, matching `ctools_rowstage.c`).
  2. Switch Phase 4 key writeback to `ctools_store_filtered_rowpar`.
  3. Replace random scatter (`buf[inv_perm[i]] = val`) with **prefetched sequential gather** (`buf[i] = raw[perm[i]]` using `CTOOLS_PREFETCH` 16–32 elements ahead, or bucketed/tiled scatter in L2-sized chunks so writes stay cache-resident!).
  4. Inline the non-AVX-512 loop instead of an indirect function-pointer call every 8 elements.
  5. Raise the auto-stream trigger threshold based on actual dataset bytes vs. available RAM (e.g., `dataset_bytes > 0.35 * free_ram`) instead of the fixed `K > 30 & N > 10M` heuristic, and document/avoid automatic `preserve` when `fast` is specified.

#### [P3-1] Post-Plugin `sort varlist, stable` Marker Cost (`nosortedby`)
- **Location:** [build/csort.ado:289](file:///G:/GitHub/stata-ctools/build/csort.ado#L289)
- **Measured Evidence:** Setting Stata's `sortedby` marker via `sort varlist, stable` after C finishes costs **0.128 s** on `x1 id` (12% of wall) and **0.451 s** on `str40 code` (**30% of total wall time!**) at N = 10M.
- **Note:** Because Stata's SPI exposes no API to set `sortedby` directly without invoking `sort`, `csort` already offers `nosortedby`. Recommend: (a) using `nosortedby` in any internal ctools ado helpers that don't require the marker, (b) noting in docs/verbose output how much time `nosortedby` saves on string keys, and (c) noting that `tempdta` mode (Section E) *does* set `sortedby` natively in the `.dta` header at zero extra cost!

---

### Section E: `cmerge` & `strL` Architecture

#### [P0-7] Eliminate `keep()`'s `hold_undo` Full-Master Load & Redundant Stata `keep if`
- **Locations:** [build/cmerge.ado:1362-1369, 1485, 1723-1736](file:///G:/GitHub/stata-ctools/build/cmerge.ado#L1485), [src/cmerge/cmerge_impl.c:1933-1941](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_impl.c#L1933-L1941)
- **Root Cause:**
  1. `cmerge.ado:1485` passes `hold_undo` whenever `keep()` is specified, forcing `cmerge_execute` (`cmerge_impl.c:1940`) to load **every column of master** into `ctools_rowstage` even when `skip_master_store == 1` (master already sorted, zero master cells modified!).
  2. `cmerge.ado:1362-1369` creates a temporary `_merge` variable even with `nogen`, writes all N codes over SPI (+31 ms at 10M), and then runs `qui keep if inlist(_merge, ...)` + `drop` in Stata (+246–283 ms at 10M) — **even when 100% of rows already match `keep()`!**
- **Measured Evidence:** On sorted `m:1` at N = 10M (`keep(3)` where all 10M rows match), total time jumps from **0.304 s → 0.852 s (2.80× slower!)**: `Load master` jumps **8.4×** (`0.028 s → 0.236 s`), `Write metadata` adds `+0.031 s`, and `Post-plugin cleanup` jumps from `0.0002 s → 0.246 s`.
- **Recommended Fix:**
  1. **Zero-drop fast path in `cmerge_plan` / `cmerge.ado` (immediate, trivial fix):** After `cmerge_plan` computes `_cmerge_N_1..5` (and if `update` isn't splitting 3 into 4/5 differently), check if the count of rows to drop (`total_obs - kept_count`) is `0`! If `0`: skip `hold_undo`, skip creating `temp_merge_var`, and skip `qui keep if inlist(...)` completely! That alone turns `0.852 s → 0.304 s` (**2.8× faster**)!
  2. **C-side `keep()` filtering when rows ARE dropped:** Pass `keep_mask` (5-bit mask for codes 1..5) into `cmerge_plan`/`cmerge_execute` (when `assert()` passes — which only depends on the pre-filter counts `_cmerge_N_1..5`!). Filter `output_specs` in C before writing to Stata so: (a) discarded rows are never written over SPI, (b) `nogen` never needs a temp `_merge` column, and (c) Stata only runs `qui drop in <kept+1>/l` (which `cmerge.ado:1701` already notes is 3×+ faster than `keep if`!).

#### [P1-9] Fuse Master Store, Using-Only Keys, `keepusing`, and `_merge` into One Row-Tiled Pass
- **Locations:** [src/cmerge/cmerge_impl.c:1984-2170](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_impl.c#L1984-L2170), [src/cmerge/cmerge_io.c](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_io.c)
- **Root Cause:** `cmerge_execute` currently hits Stata's row-major table in **3–4 separate passes**:
  1. Pass 1 (`:1992`): `ctools_rowstage_store` writes master rows.
  2. Pass 2 (`:1996-2013`): **Single-threaded** loop writes using-only key values.
  3. Pass 3 (`:2051-2127`): Spawns `n_keepusing * ku_chunks` tasks where **different threads write different `keepusing` columns to the exact same row ranges simultaneously** — causing severe cache-line false sharing / ping-pong on Stata's row-major memory layout!
  4. Pass 4 (`:2150-2169`): Separate OpenMP loop writes `_merge`.
- **Measured Evidence:** `Write metadata` (`keepusing` + `_merge`) takes **0.070–0.107 s** at N = 10M (37% of C time on sorted `m:1`, for writing only 2–3 variables!).
- **Recommended Fix:** Partition `[0, output_nobs)` into disjoint contiguous row tiles (`[r0, r1)`, 1 tile per thread). Inside each tile, for each row `i`, one thread writes:
  - staged master row (if `write_rows` and `master_orig_rows[i] >= 0`), OR using-only key columns (if `master_orig_rows[i] < 0`),
  - all `keepusing` columns (tracking per-row `update`/`replace` flags in a stack variable — eliminating the `output_nobs`-byte `row_flags` heap allocation and its serial scan at lines 2140–2149!),
  - `_merge` code (if `merge_var_idx > 0`).
  Each Stata observation cache line is touched **once** by **one** thread, eliminating false sharing and cutting DRAM write traffic substantially.

#### [P1-10] `cmerge_load_using`: Replace `ips4o` with `SORT_ALG_AUTO`, Combine Loads, Skip Eager Permute
- **Location:** [src/cmerge/cmerge_impl.c:535-665](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_impl.c#L535-L665)
- **Root Cause:**
  1. Calls `ctools_data_load` **twice** (keys, then `keepusing`) instead of a single combined load.
  2. Hardcodes `ctools_sort_ips4o_with_perm` (`O(N log N)` comparison sort) instead of `ctools_sort_dispatch(..., SORT_ALG_AUTO)` (parallel counting/radix `O(N)` sort — which master already uses!).
  3. Eagerly permutes every `keepusing` column in parallel only over `v < n_keepusing` (using only 2–3 threads when `n_keepusing = 2–3`).
- **Measured Evidence:** In A2 (`1:1` 2M × 2M unsorted), `Phase 1` spent `0.0184 s` sorting using with `ips4o` + `0.0140 s` permuting 3 keepusing cols = `0.0324 s` (vs. `0.0176 s` for master's radix sort without eager permute). On string keys (`str40`), `ips4o` is **6.6× slower** than MSD radix (`0.404 s` vs `0.061 s` at N = 2M!).
- **Recommended Fix:** Use `ctools_sort_dispatch(..., SORT_ALG_AUTO)` to get `using_perm[sorted_pos] = orig_row`, store `using_perm[using_sorted_row]` directly in `output_specs[i]` during the join, and **delete eager `keepusing` permutation entirely**! Also load keys + `keepusing` in a single `ctools_data_load` call.

#### [P1-11] Parallelize `cmerge_plan` and `cmerge_sorted_join`
- **Locations:** [src/cmerge/cmerge_join.c](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_join.c), [src/cmerge/cmerge_impl.c:332, 1064-1166](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_impl.c#L1064-L1166)
- **Root Cause:** `cmerge_keys_in_order` (sortedness check), uniqueness pre-checks, `cmerge_sorted_join`, `orig_rows` construction, `seen`/`head` stable partition, and result counting are all **100% single-threaded**.
- **Measured Evidence:** At N = 10M, `Merge join = 0.036–0.044 s` (19% of sorted `m:1` C time; at 2M `1:1`, `Merge join = 0.035 s` — **twice as slow as sorting 2M rows!**).
- **Recommended Fix:**
  1. Parallelize `cmerge_keys_in_order` and uniqueness checks across `T` contiguous chunks with boundary checks.
  2. **Parallel merge-join via binary-search partitioning (Merge-Path):** Split the larger side (`master`) into `T` chunks at key group boundaries, binary-search each chunk's boundary key in `using` (O(log N_using)), run independent 2-pointer joins in parallel into thread-local count/output ranges, prefix-sum output offsets, and scatter `output_specs` in parallel!
  3. Parallelize `orig_rows` and `seen`/`head` counting with thread-local prefix sums.

#### [P2-1] `strL` Acceleration for `cmerge` and `csort`
- **Locations:** [src/cmerge/cmerge_impl.c:146-219, 1497-1616](file:///G:/GitHub/stata-ctools/src/cmerge/cmerge_impl.c#L146-L219), [build/cmerge.ado:134, 1675-1687](file:///G:/GitHub/stata-ctools/build/cmerge.ado#L1675-L1687), [build/csort.ado:105-112, 158](file:///G:/GitHub/stata-ctools/build/csort.ado#L105-L112), [src/formats/ctools_dta.c](file:///G:/GitHub/stata-ctools/src/formats/ctools_dta.c)
- **Measured Evidence (N = 2M):**
  - `cmerge m:1` with **master `strL`**: native `merge` = **2.374 s** vs. `cmerge, tempdta` = **1.865 s** (`Load master = 0.614 s` single-threaded `cmerge_read_strl`, `Store/dta write = 0.236 s`, `Post-plugin use + meta = 0.893 s`).
  - `cmerge m:1` with **using `strL`** (10k lookup): native `merge` = **0.522 s** vs. `cmerge, tempdta` = **0.729 s** (37% **slower** than native!).
  - `csort x1 id` with **1 non-key `strL`**: **0.961 s** (unconditional fallback to native `sort, stable`, 3.4× slower than `csort`).
- **Recommended Fixes (4 high-impact upgrades):**
  1. **Hybrid SPI + Mata fast path for using-only `strL`s in `cmerge` (beats native `merge` cleanly!):** When master has **no `strL`s** and only `keepusing` has `strL`s (the `lookup_strl` case!), do **not** use `tempdta` (which dumps all 22 master columns to disk!). Instead:
     - Run the normal fast SPI path for all non-`strL` columns (`0.13 s` at 2M).
     - For `m:1` (or `1:1`) using `strL` columns: either have C write an `int32` index column (`using_row`) to Stata and let **Mata** populate the `strL` column via `st_sstore(., var, u_strl[idx])` (avoiding disk I/O of all 22 master columns!), OR check if `strL` max length ≤ 2045 and text-only (fast `str#` SPI write + `recast strL` if lossless!).
  2. **Parallelize `cmerge_execute_dta` and optimize `cmerge_read_strl`:**
     - Currently `cmerge_execute_dta` (`cmerge_impl.c:1497-1616`) builds `mrow`/`urow`, keys, shared vars, and `_merge` **single-threaded** before handing off to `ctools_dta_write`. Parallelize all array construction.
     - In `cmerge_read_strl` (`:146-219`), test whether `SF_strldata` is thread-safe across disjoint variables/chunks (if SPI requires serial `strL` reads, at minimum eliminate the redundant `SF_var_is_binary` + `SF_sdatalen` calls when `SF_strldata` already returns exact length / null-terminated ASCII!).
     - Fix O(K²) row-appending `C = C \ (...)` in Mata's `_cmerge_dta_meta_save` (`cmerge.ado:134`) by preallocating `J(total_chars, 3, "")`.
  3. **Direct C `.dta` Reader for `cmerge using` (biggest architectural win for `cmerge` with `strL` and large using files!):**
     - Currently `cmerge.ado` opens the `using` `.dta` file in a Stata frame (`use`), then copies keys + `keepusing` out via SPI, then drops the frame (`Pre-plugin1 + Inter-plugin` overhead!).
     - Implementing a lightweight native C `.dta` 117–119 reader (mirroring the existing [ctools_dta.c](file:///G:/GitHub/stata-ctools/src/formats/ctools_dta.c) writer, or leveraging `cio_readstat.c` for non-strL + direct GSO reader for `<strls>`) lets `cmerge` read ONLY the key + `keepusing` columns directly from the `.dta` file in C — **completely bypassing Stata `use`, frame creation/drop, AND serial `SF_strldata` reads for using `strL`s**, while also reading `.dta` value labels and variable metadata directly!
  4. **Enable `tempdta` (and optional `recast`) path for `csort` with `strL`:**
     - `csort` currently refuses any dataset containing a `strL` even if `strL` isn't a sort key. Add a `tempdta` path reusing `ctools_dta_write` (setting `<sortlist>` directly in the `.dta` header so no post-plugin `sort` is needed!) for large datasets where parallel C sorting + parallel `.dta` write + `use` beats single-threaded/4-core Stata `sort, stable` (especially on string keys or wide datasets with 1 `strL`), AND a zero-disk fast path when `strL` values are ≤ 2045 bytes non-binary (or Mata-only permutation of the `strL` column using `sort_order`!).

#### [P2-2] Trim Stata `.ado` Wrapper Overhead (`cmerge`, `creghdfe`, `civreghdfe`, `cppmlhdfe`)
- **Measured Evidence at N = 10M:**
  - `civreghdfe`: **0.75–0.85 s (28–32% of wall time)** (`Setup 0.30–0.36 s`, `Plugin call overhead 0.27–0.29 s`, `Post-processing 0.18–0.20 s`)
  - `creghdfe`: **0.42–0.72 s (28–35% of wall time)** (`Setup 0.21–0.26 s`, `Plugin call overhead 0.12–0.28 s`, `ereturn posting 0.16–0.18 s`)
  - `cppmlhdfe`: **1.12 s Setup + 0.34 s ereturn posting** (dominated by P0-2's 4 unnecessary tempvars + `qui count if depvar < 0 & touse` + `qui count if touse`!)
  - `cmerge`: **`Inter-plugin work = 0.096–0.128 s`** (32% of sorted `m:1` wall time!)
- **Recommended Fixes:**
  1. **`cppmlhdfe.ado` / `civreghdfe.ado`:** Move `count if depvar < 0 & touse` (`cppmlhdfe.ado:729`) and `count if touse` (`cppmlhdfe.ado:784`, `civreghdfe.ado:1186`) into the C plugin (which already scans `y` and `touse` in parallel in ~5 ms instead of ~150 ms in Stata!).
  2. **`creghdfe.ado` / `civreghdfe.ado` / `cppmlhdfe.ado`:** Avoid `fvrevar` when `varlist` contains no factor/time-series operators (already fast-checkable via `strpos` for `.` and `#`!), and skip `ereturn post ..., esample(...)` tempvar copying when `nosample` is set or reuse `touse` directly when 0 singletons/separated obs are dropped!
  3. **`cmerge.ado` Inter-plugin (`0.096–0.128 s`):** Creating a temp `.dta` on disk (`_cmerge_template`) and running `qui append using empty_vars_template` on every merge costs ~40–70 ms of disk + Stata append overhead even when adding 2 numeric variables with no labels/notes. When `new_var_names` has no special value labels/notes/chars (the common case), creating 1–5 variables via `qui gen` / `mata: st_addvar()` directly in memory avoids writing & reading a temporary `.dta` file on disk. (Note: the comment at line 1439 cites an older benchmark vs. `st_addvar` for *many* variables — switch on `n_new_vars <= 8 && !has_notes_or_chars`!)

---

## 5. Recommended Implementation Order

### Wave 1 — High-Impact Quick Wins (1–3 days, Low Risk, Massive Measured Gains)
1. **[P0-2]** `cppmlhdfe`: Gate `saved_d`, `separation_tag`, `separation_certificate` tempvars & their 30M serial `SF_vstore()` calls; parallelize `final_sample` store (or reuse `touse` when 0 obs dropped); move `count if` checks to C.
2. **[P0-3]** `cppmlhdfe`: Replace Pass 1 `exp()` overflow check with `> 709.7827`; fuse & parallelize `ppml_compute_deviance` with Pass 2 `mu`; parallelize & deduplicate `ll` / `ll_0` with integer `y` LUT and closed-form `ll_0`.
3. **[P0-7]** `cmerge`: Skip `hold_undo`, temp `_merge`, and Stata `keep if` when `cmerge_plan` shows 0 rows dropped; pass `keep_mask` into C to filter dropped rows before SPI writeback.
4. **[P0-1 + P1-1]** Cross-cutting SIMD/FMA: Add runtime CPUID AVX2+FMA dispatch (`__attribute__((target("avx2,fma")))`) for `two_prod`, `kahan_dot`, `CTOOLS_MEAT_ADD`, `dd_dot_weighted_multi`, `ctools_simd_dot_weighted`, `fused_axpy2_dot/wdot`, and `axpby` + Dekker fallback.
5. **[P0-5 Quick Win 1 & 2]** `hdfe_solver.c`: Move convergence check before `axpby`; implement 4-pass algebraic `G=2` Kaczmarz with fused `uv` dot.
6. **[P1-3]** Replace `qsort()` in `hdfe_remap.c` and `ctools_hdfe_utils.c` with parallel radix sort / hash remap.
7. **[P1-10]** `cmerge_load_using`: Switch `ips4o` → `SORT_ALG_AUTO`, combine loads, drop eager `keepusing` permutation.

### Wave 2 — Core Engine Parallelism & Kernel Fusion (1–2 weeks)
8. **[P0-4 + P1-2]** `cppmlhdfe`: Merge STEP 8b into IRLS iter 0, deduplicate simplex/ReLU projections, fast-path ReLU binary weights, and evaluate warm-started (`exact_partial=0`) intermediate IRLS iterations with final exact verification.
9. **[P0-6]** `hdfe_solver.c`: Add row-parallel / hybrid 2D column+row parallelism to CG solver so `K = 1..4` utilizes all 16–32 cores.
10. **[P1-4 + P1-5 + P2-3]** `civreghdfe` & `creghdfe`: Single joint parallel Gram matrix in `civreghdfe`, row-blocked `ctools_matmul_atdb`, zero-drop compaction bypass, cache-friendly scatter-add `ctools_vce_cluster`, remove `volatile` RSS, optimize CUE moments/Hessian.
11. **[P1-6 + P1-7 + P1-8]** `csort`: Parallel 8-byte-prefix string MSD sort; fuse key extraction into `ctools_rowstage_load`; scratch buffer reuse for multi-GB `rs->rows`; fix `csort_stream.c` (`SD_SAFEMODE`, prefetched gather, parallel key store).
12. **[P1-9 + P1-11]** `cmerge`: Single-pass row-tiled store fusion (master + using keys + `keepusing` + `_merge`) and parallel Merge-Path join.

### Wave 3 — Architectural Upgrades (2–4 weeks)
13. **[P0-5 Breakthrough A & B]** `hdfe_solver.c`: Dense Schur-complement direct solver for 2-way FE with `L_min <= 2000` (zero CG iterations!) + bipartite edge-list CSR CG in parameter space for large `L_min > 2000`.
14. **[P2-1]** `strL` & `.dta` architecture: Hybrid SPI + Mata `strL` path for using-only `strL` merges, parallel `cmerge_execute_dta`, native C `.dta` 117–119 using reader, and `strL`-capable `csort`.

---

## 6. Verification Checklist for Implementing These Changes

- **Numerical Equivalence:** Verify coefficients, SEs, `t`/`z`, `R²`, `deviance`, `ll`, `ll_0`, `df_a`, and `e(sample)` against existing test suites (`validation/`) — especially for AVX2+FMA `two_prod`, `G=2` Schur complement / algebraic Kaczmarz, `cppmlhdfe` `use_exact_partial(0)`, and `ctools_vce_cluster` scatter-add.
- **Cross-Platform Binary Safety:** Confirm `__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")` builds cleanly on MinGW-w64 (Windows), GCC (Linux), and Clang (macOS x86-64) while keeping base `-march=x86-64` portability for older CPUs and ARM64 NEON FMA intact.
- **Re-Benchmark at N = 2M & N = 10M:** Track both wall clock and C internal phase breakdowns (including `csort` `cleanup_time` and `threads(1/4/16/32)` scaling) after each wave.
