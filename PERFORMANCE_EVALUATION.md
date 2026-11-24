# Evaluation of PERFORMANCE_IMPROVEMENTS.md

**Date:** 2026-11-24 · **Branch:** `perf/audit-eval` (from `dev` at a954de3)
**Machine:** the audit's machine (Ryzen 9 9950X, 16C/32T, 94 GB, dual-channel DDR5), Windows 11, StataNow 19.5 MP (4-core)
**Toolchain:** MSYS2 MINGW64 GCC 16.2.0 (the CI toolchain), plus clang 22.1.8 + LLVM libomp for comparison

Every proposal in the audit was checked against measurements, not taken on trust. Low-effort items were implemented and A/B-tested in Stata. The rest were measured by phase profile or option ablation.

## Bottom line

| Command (N = 10M) | Before | After | Speedup | Result change |
|---|---:|---:|---:|---|
| `cppmlhdfe count_z x1 x2, absorb(firm group)` | 21.2 s | 14.7 s | **1.44x** | ≤ 4e-16 relative |
| `cppmlhdfe count_pos x1 x2, absorb(firm group)` | 18.0 s | 9.3 s | **1.94x** | ≤ 7e-16 relative |
| `cppmlhdfe count_z ..., vce(cluster firm)` | 21.9 s | 14.3 s | **1.53x** | ≤ 4e-16 relative |
| `cmerge m:1` sorted master, `keep(3)` (± `nogen`) | 0.90 s | 0.33 s | **2.9x** | identical data signature |
| `creghdfe` with sparse FE ids (`firm*1000003`) | 2.81 s | 1.84 s | **1.54x** | bitwise identical coefficients |
| `civreghdfe` (unadj / robust / cluster) | 2.75–3.10 s | 2.63–3.02 s | 1.00–1.05x (medians 1.05–1.07x) | ≤ 9e-16 relative |
| `creghdfe` (compact FE ids, K=3 or 9) | 1.79–4.37 s | 1.79–4.23 s | 1.00–1.04x | ≤ 1.3e-15 relative |
| `csort`, unsorted `cmerge` | — | — | unchanged | not touched |

Times are the best of 4 runs (3 per Stata session × 2 interleaved sessions, first run of each session dropped). Run-to-run noise on this machine is ±5–10%, so effects under ~5% are not distinguishable.

**The audit's attributions are mostly right in direction but often wrong in size.** The biggest real wins were not its headline P0-1 (software `fma`) or P0-6 (row-parallel CG). They were the serial transcendental and double-double loops in `cppmlhdfe` (P0-3), the `cmerge keep()` path (P0-7), and the `qsort` remap (P1-3).

## Verdict on every item

| ID | Audit claim | Verdict | Measured evidence | Action |
|---|---|---|---|---|
| P0-1 | Software `fma()` makes OLS/VCE kernels 5–15× slower | **Overstated; but found a real accuracy bug** | Microbenchmark: MinGW `fma()` 1.96 ns/elem vs hardware 0.90 (2.2×, not 5–15×). End-to-end `-mfma` build: 0–7%. **MinGW's `fma()` is not correctly rounded: it disagreed with hardware FMA in 8.5% of 20M random products**, so `two_prod` was not error-free on Windows | **Fixed**: guarded Dekker split product (exact, bitwise equal to hardware FMA; 0 mismatches in 20M), with the FMA instruction for out-of-range operands |
| P0-2 | 4 unused tempvars + 40M serial `SF_vstore` cost 2.5–4.5 s | **Real but ~15× overstated** | ado-side gap 0.52–0.57 → 0.45–0.51 s, in-plugin store phase 0.43 → 0.31 s: ≈0.2 s total | **Fixed**: `d`, `tagsep` and `zvarname` variables are created only when requested |
| P0-3 | Serial deviance/log-likelihood with redundant `exp` | **Confirmed — the largest PPML win** | deviance 2.46 → 0.37 s; ll + ll_0 1.49 → 0.39 s; final solve 1.20 → 0.54 s | **Fixed**: fixed-block parallel double-double sums (deterministic, thread-count independent); `exp()` skipped in the overflow pre-check below 709 |
| P0-4 | Simplex + ReLU cost 5.3 s with 0 separated obs; STEP 8b 0.6 s | **Confirmed** | Ablation: simplex 2.74 s, ReLU 3.10 s, STEP 8b 0.60 s | **Partly**: parallel 4-lane `ppml_quad_dot` cut simplex + ReLU to ~4.3 s. The algorithmic dedup (reusing the basis, merging STEP 8b) is not done: it changes iteration paths |
| P0-5 | Quick win 1: `axpby` before convergence test | Real but negligible (~1%) | — | **Done** (bitwise identical) |
| P0-5 | Schur complement / edge-list CG for G=2 | Plausible, untested | HDFE partial is 0.62–0.81 s for K=3 and 2.2 s for K=9 | Not done: a new solver with different (more exact) results |
| P0-6 | Row-parallel CG gives 2–4× for small K | **Mostly refuted on this hardware** | K=9: 1 thread 7.14 s, 3 → 4.46, 9 → 3.80, 16 → 3.86. Past 3–4 concurrent column solves the CG is memory-bandwidth-bound | Not done. Expect ≤ 15–25% of the CG phase here; may pay on many-channel servers |
| P0-7 | `keep()` loads the full master + runs `keep if` with 0 drops | **Confirmed exactly** | 0.90 → 0.33 s (audit predicted 0.85 → 0.30) | **Fixed**: the plan step drops the undo copy when keep() can drop nothing; the ado skips `keep if` when counts show no drops |
| P1-1 | Dead AVX2 code; single-accumulator scalar fallbacks | **Partly** | `-march=x86-64-v3` adds 3–11% on regressions | **Fixed (baseline path)**: 8-accumulator fallbacks in the AVX2 lane order for dot, weighted dot, fused axpy-dot kernels |
| P1-2 | Exact partialling costs ~0.7 s | Confirmed | `use_exact_partial(0)` saves 0.34 s (before) / 0.70 s (after) | Not changed: it is a default-numerics decision |
| P1-3 | `qsort` remap makes sparse FE ids 14× slower | **Confirmed** | FE remap 1.05 s vs 0.067 s for compact ids | **Fixed**: O(N) hash remap; levels in first-appearance order like the counting path; results bitwise identical |
| P1-4 | `civreghdfe` redundant/serial Gram matrices | Plausible | IV + VCE 0.36–0.60 s | Partly via P1-1 (+5–7%); restructuring not done |
| P1-5 | Cluster VCE gather; `volatile` RSS | `volatile` refuted as meaningful (~7 ms); VCE 0.13–0.40 s plausible | — | Not changed |
| P1-6 | `str17+` MSD sort is serial | Confirmed in direction | `str40` sort 0.345 s vs numeric 0.057 s | Not changed |
| P1-7 | 0.15–0.27 s unaccounted alloc/free in `csort` | **Overstated** | gap ≈ 0.05 s at 22 variables | Not changed |
| P1-8 | `stream()` mode much slower | Confirmed | `stream(8)` 3.40 s vs default 1.42 s (2.4×) | Not changed. `SD_SAFEMODE` there is deliberate per its source comment, so the audit's first fix is unsafe as stated |
| P1-9 | Multi-pass `cmerge` stores | Confirmed | write metadata 0.10 s | Not changed |
| P1-10 | `ips4o` + eager permute for using | Confirmed | 1:1 10M: sort using 0.104 s + permutation 0.074 s | Not changed |
| P1-11 | Serial join/plan | Confirmed | join 0.037 s (m:1) / 0.197 s (1:1) | Not changed |
| P2-1 | `strL` paths | Not evaluated | — | — |
| P2-2 | ado overhead | **Confirmed and large** | `creghdfe`: 0.6–0.86 s of ~2.0 s wall in Stata (setup 0.25, plugin-call overhead 0.15, `ereturn` 0.26); `civreghdfe` 0.95 s; `cppmlhdfe` setup 1.2–1.6 s | Not changed. For `creghdfe` this is now a bigger lever than any C kernel |
| P2-3 | CUE moments | Not evaluated | — | — |
| P2-4 | Setup micro-optimizations | Low value | remap + counts 0.05–0.07 s; singleton passes small | Not changed |
| P3-1 | `sortedby` marker cost | Confirmed | 0.17 s (numeric keys), 0.43 s (`str40`, 35% of wall) | Not changed (documented `nosortedby` exists) |

## Where `cppmlhdfe` time goes now (10M, `count_z`, in-plugin)

| Phase | Before | After |
|---|---:|---:|
| IRLS CG partial (6 iterations, 18 sub-iterations) | 3.18 s | 3.16 s |
| ReLU + simplex separation | 5.56 s | 4.26 s |
| Deviance (6×) | 2.46 s | 0.37 s |
| ll + ll_0 | 1.49 s | 0.39 s |
| WLS normal equations (6×) | 1.13 s | 0.15 s |
| Final projection + solve | 1.20 s | 0.53 s |
| STEP 8b collinearity CG | 0.60 s | 0.64 s |
| VCE | 0.42 s | 0.43 s |
| **Total in plugin** | **18.0 s** | **11.4 s** |

The remaining big items are the separation checks (P0-4, algorithmic) and the CG itself (P0-5 Schur complement).

## Plugin build options (all platforms)

Measured on Windows at 10M (regressions, IV, PPML), relative to the CI flags `-O3 -march=x86-64 -mtune=generic -ffp-contract=off -flto`:

| Variant | Effect | Results |
|---|---|---|
| `-march=x86-64-v2` | none measurable | bitwise identical |
| `-mfma` only | 0–7% | bitwise identical |
| `-mavx2 -mfma` / `-march=x86-64-v3` | 3–11% on regressions | ~1e-16 relative (FMA in reductions) |
| clang 22 + static LLVM libomp (Windows) | 0.87–1.05× (no better) | ~1e-15 relative |
| `-fno-math-errno -fno-trapping-math` | none measurable | bitwise identical |

Recommendations:

- **Windows.** Keep MSYS2 GCC + static libgomp. Clang + libomp is not faster. The real Windows-specific problem was MinGW's inaccurate software `fma()`, now handled in source. Switching to the UCRT runtime would give a correctly rounded `fma` and a newer libm, but it changes `printf`/`strtod` behavior under every I/O command. Not worth it without a full I/O validation pass.
- **Linux and macOS Intel.** Keep the `x86-64` baseline. `v3` buys only 3–11% and would drop pre-Haswell CPUs, some VMs, and Rosetta 2 before macOS 15. If those few percent matter, dispatch at run time with `__attribute__((target("avx2")))` plus `__builtin_cpu_supports` on the reduction kernels. The new 8-accumulator fallbacks use the AVX2 lane order, so an AVX2 clone without FMA would give **bitwise identical results on every CPU**.
  - glibc's and Apple's `fma()` were already correctly rounded; the Dekker path gives them the same exact results without a library call.
  - `-fno-semantic-interposition` is likely a no-op given LTO plus the export map; I could not measure Linux here.
  - Moving to `manylinux_2_28` for a newer GCC would raise the glibc floor from 2.17 to 2.28.
- **macOS arm64.** Already well configured (`-mcpu=apple-m1`, NEON paths, hardware FMA). No change.
- **CI build time.** The critical path is the Linux native validation step (~400 s of a ~565 s job), not compilation.
  - The macOS jobs rebuild LLVM OpenMP from source every run (~47 s on Intel). It could be cached like the Linux runtime.
  - Parallel compilation in `scripts/build_plugin.py` would conflict with `validation/check_build.py`'s "abort on first compile failure" contract, so I left it out.

## Changes on this branch

- `src/ctools_ols.h`, `src/ctools_ols.c`, `src/ctools_matrix.c`: exact `CTOOLS_PROD_ERR` (P0-1). Also `dd_dot_weighted_multi_blocked`, a deterministic block-parallel double-double multi-dot.
- `src/hdfe/hdfe_remap.c`: hash remap for sparse or non-integer FE ids (P1-3).
- `src/hdfe/hdfe_solver.c`: CG convergence test before the direction update (P0-5 QW1).
- `src/ctools_simd.h`: 8-accumulator baseline fallbacks (P1-1).
- `src/cppmlhdfe/cppmlhdfe_moments.c`, `cppmlhdfe_irls.c`, `cppmlhdfe_solvers.c`: block-parallel deviance, log-likelihood, normal equations, RSS and constant sums; parallel final `exp` loop; `exp`-free overflow pre-check; block-parallel 4-lane `ppml_quad_dot` (P0-3, P0-4 partial).
- `build/cppmlhdfe.ado`: optional output variables only when requested (P0-2).
- `src/cmerge/cmerge_impl.c`, `build/cmerge.ado`: `keep()` zero-drop fast path (P0-7).

The ado/plugin protocol is unchanged. The new `keep_mask` plan token is ignored by older plugins, and old ado files still work with the new plugin.

**Numerical note.** Parallel sums use fixed 32,768-row blocks combined in block order, so results do not depend on the thread count or machine. They can differ from the previous release in the last bit or two: ≤ 1.3e-15 relative on every benchmark.

## Validation

These were run with the rebuilt Windows plugin and the changed ado files.

- **Stata suites: all pass.**
  - `cppmlhdfe` 1302 checks, `civreghdfe` 480, `s2_iv` 499, `creghdfe` 385, `cmerge` 362, `cbinscatter` 281, `s2_hdfe` 102, `threads` 24.
  - `cppmlhdfe_publication` 8, `creghdfe_replay` 8, `civreghdfe_grouped_syntax` 17.
  - `chdfe`: 585/585 cases.
- **Bootstrap and wild suites: pass, except two pre-existing nondeterministic checks.**
  - Passing: `creghdfe_bootstrap` 17, `creghdfe_fe_bootstrap` 7, `creghdfe_wild` 20, `civreghdfe_bootstrap` 24, `civreghdfe_fe_bootstrap` 7, `civreghdfe_savefe` 7, `civreghdfe_wild` 31, `cppmlhdfe_wild` 17, `bootstrap_helpers` 13, `bootstrap_options` 35, `bootstrap_reference` 19, `wild_draws` 16.
  - **Not caused by this branch:** `cppmlhdfe_bootstrap` failed 1 of 14 checks and `cppmlhdfe_fe_bootstrap` 1 of 8. Both are exact-equality checks.
    - `draws(fast)`: the thread-count invariance check fails **at HEAD** too (e(V) differs by 1.3e-18).
    - "shared replicates": FE standard errors differ between runs **at HEAD** in 3 of 3 standalone repetitions, by up to 5e-10 relative.
    - With `threads(1)` both checks are exact on HEAD and on this branch. Multi-threaded `cppmlhdfe` bootstrap is therefore not reproducible run to run on Windows: a pre-existing bug to fix separately.
- **Native C tests on MinGW: pass.** These are `test_hdfe_paths_native`, `test_hdfe_cancel_native`, `test_bins_paths_native`, `test_bootstrap_engine_native` and `test_wild_engine_native`. UBSan ran in trap mode, since MinGW has no libubsan.
  - The remaining native suites need ASan, `dlfcn.h` or `posix_memalign`, so they run only on Linux CI.
- **Exactness:** the Dekker product error matched hardware FMA bitwise on 20M random operand pairs across exponent ranges [-1020, 1020].
- **Not run locally:** Linux and macOS builds. CI will cover them.

## Reproducing

The local harness lives in `temp/perf/` (gitignored):

- `gen.do`: synthetic 22-variable panels at 2M/10M rows.
- `bench.do`: workloads.
- `run.py`: interleaved A/B runner with result comparison.
- `build_variant.sh`: plugin variants into `temp/bench/<name>/`.
- `micro/`: `fma` and exactness microbenchmarks.
- `cppmlhdfe_irls.instrumented.c`: the phase profiler used for the PPML tables.
