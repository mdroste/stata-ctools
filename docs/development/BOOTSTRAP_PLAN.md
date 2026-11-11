# Bootstrap inference for creghdfe, cppmlhdfe and civreghdfe

Implementation plan. Status: Phase 0 (reference behavior, including boottest)
done; see Appendix A. Phases 1, 2a and 2b (creghdfe `vce(bootstrap)`) done; see
§10. Written against the
1.2.4 sources; file and line references point at that tree.

## Contents

1. [Summary](#1-summary)
2. [Goals and non-goals](#2-goals-and-non-goals)
3. [Where the code stands today](#3-where-the-code-stands-today)
4. [User-facing design](#4-user-facing-design)
5. [Statistical design](#5-statistical-design)
6. [Architecture](#6-architecture)
7. [creghdfe](#7-creghdfe)
8. [cppmlhdfe](#8-cppmlhdfe)
9. [civreghdfe](#9-civreghdfe)
10. [Phases and milestones](#10-phases-and-milestones)
11. [Validation and benchmarks](#11-validation-and-benchmarks)
12. [Documentation and release plumbing](#12-documentation-and-release-plumbing)
13. [Risks](#13-risks)
14. [Open questions](#14-open-questions)
15. [Later work](#15-later-work)
16. [Appendix A. Phase 0 findings](#appendix-a-phase-0-findings)

---

## 1. Summary

Add a nonparametric (pairs) bootstrap to the three HDFE estimators. Observations
or clusters are resampled, optionally within strata. Inference can target:

- **the regressors:** `vce(bootstrap, ...)` replaces the analytic e(V) with the
  bootstrap covariance of e(b);
- **the fixed effects the user saves:** a new `fevce(bootstrap, ...)` writes a
  bootstrap standard error next to each saved fixed effect;
- **both**, computed from the same replicates in one pass.

A second method, the **wild cluster bootstrap**, is available for the
coefficients only, through a separate `wildbootstrap()` option (§4.7, §5.9).
It keeps the analytic cluster-robust e(V) and adds bootstrap p-values and
confidence intervals for coefficient tests. It never applies to saved fixed
effects; `fevce()` always uses the pairs bootstrap. The method by command:

- creghdfe: wild restricted (WCR) bootstrap;
- civreghdfe: wild restricted efficient (WRE) bootstrap;
- cppmlhdfe: the score-based analog for nonlinear models (Kline–Santos).

With two or more absorbed effects, each fixed effect is identified only up to a
normalization. Whenever the command writes per-level fixed-effect standard
errors in that case, it prints a warning (§4.8). The warning explains the
caveat and points to `dse()`, the standard error of the summed effect, which
does not depend on the normalization.

Every replicate runs inside the plugin on data that was loaded once. No
replicate goes back to Stata. A bootstrap draw becomes an integer multiplicity
per row, and each replicate is a weighted re-fit of the retained full sample.
The shared fixed-effect solver already supports observation weights.
Replicates run in parallel, one per thread, each with its own `HDFE_State`. The
full-sample fit supplies warm starts for the projection and, in PPML, for the
IRLS iterations.

**Random numbers.** Stata's random-number generator drives every draw, so
`set seed`, `seed()`, `set rng` and `set rngstream` work as they do for Stata's
`bootstrap`, and e(rngstate) is posted.

- **Default (`draws(stata)`, decided in Phase 5):** every draw is made with
  Stata's generator, reproducing the replicates of `bootstrap … : cmd` for the
  same seed (§5.8).
- **Opt-in (`draws(fast)`):** Stata's RNG seeds a per-replicate in-plugin
  generator. Results are reproducible, and they do not depend on `threads()` or
  on how replicates are split into chunks.

The work is staged:

1. A shared engine in `src/hdfe/`.
2. creghdfe: pairs bootstrap of the coefficients, then of the fixed effects
   (with the normalization warning), then the wild cluster bootstrap.
3. cppmlhdfe: pairs bootstrap, then the score wild bootstrap.
4. civreghdfe: fixed-effect saving first, which it does not have today; then
   the pairs bootstrap; then the WRE bootstrap.

## 2. Goals and non-goals

### Goals

| # | Goal |
|---|------|
| G1 | `vce(bootstrap[, bootopts])` on all three commands: an observation, cluster or stratified-cluster pairs bootstrap of e(b). It posts the bootstrap e(V) and the result set that `estat bootstrap` reads. |
| G2 | Bootstrap standard errors for every fixed effect the user saves: `absorb(..., savefe)` in all three commands, plus `absorb(newvar=absvar)` and `d()` in cppmlhdfe. Optional percentile bounds and valid-replicate counts. |
| G3 | Each target can be used alone or with the other. In FE-only mode the coefficient VCE stays analytic (for example `vce(cluster firm)`). |
| G4 | Much faster than `bootstrap: <cmd>`. One data load, parallel replicates and warm starts; no per-replicate ado work, SPI transfer, remap or write-back. |
| G5 | Stata's random-number generator drives every draw, and the commands treat it as Stata's `bootstrap` does: `set seed`, `seed()`, `set rng` and `set rngstream` are honored, e(rngstate) is posted, and Stata's RNG is left advanced past the draws taken. Identical seeds give identical draws, and identical results for a fixed `threads()`, whatever the chunking. Draws do not depend on `threads()`. `draws(stata)` reproduces the exact replicates of Stata's `bootstrap` prefix for the same seed (§5.8). |
| G6 | No change when no bootstrap is requested. e(b), e(V) and every other result stay bitwise identical. |
| G7 | Interruptible. Break between chunks leaves the data and e() as the command's existing failure path leaves them, and frees the retained plugin state. |
| G8 | `wildbootstrap([wildopts])` on all three commands, for coefficients only. It gives bootstrap p-values and test-inversion confidence intervals for each tested coefficient, with the null imposed by default. creghdfe uses WCR, civreghdfe uses WRE (2SLS, LIML, Fuller, k-class), and cppmlhdfe uses the score bootstrap. When every absorbed effect is nested within the bootstrap clusters, a replicate costs O(C·K) and does not depend on N. |
| G9 | With two or more absorbed effects, writing per-level fixed-effect SEs always prints a warning about the normalization caveat and points to `dse()`. The warning is recorded in e(fe_norm_warning). |

### Non-goals for this plan

These are listed in [Later work](#15-later-work):

- the wild bootstrap for fixed effects (excluded by design: `fevce()` is pairs
  only);
- joint (multi-coefficient) wild bootstrap tests, multiway wild bootstrap, and
  the wild bootstrap for civreghdfe's GMM2S and CUE;
- jackknife, bootstrap-t, BCa intervals and m-out-of-n resampling;
- multiway "pigeonhole" resampling and time-series block bootstrap (other than
  resampling panels as clusters);
- bootstrapping civreghdfe's diagnostic statistics;
- user-supplied replicate weights (`svy bootstrap`).

In the `draws(fast)` mode, the pairs bootstrap's draws are a
reproducible function of Stata's seed but are not `bsample`'s draws. cbsample
already documents the same difference in `docs/COMPATIBILITY.md`.
`draws(stata)` gives `bsample`'s draws (§5.8).

## 3. Where the code stands today

### Common to all three commands

- **vce(bootstrap) is recognized and then rejected.** `_ctools_vce` classifies
  `boot*` as `bootstrap` (`build/_ctools_vce.ado:18`). Each command then rejects
  it:
  - creghdfe: "unrecognized vce() option", r(9) (`build/creghdfe.ado:207-209`);
  - cppmlhdfe: r(9) (`build/cppmlhdfe.ado:424-426`);
  - civreghdfe: r(198) (`build/civreghdfe.ado:348-350`).

  `validation/stata/creghdfe.do:3328` asserts that creghdfe rejects it as
  reghdfe does. That test changes with this work.
- **The `bootstrap:` prefix already works**, because each command is eclass and
  posts e(b), e(V) and e(sample). Every replicate, though, repeats ado parsing,
  fvexpand, `bsample`, a full SPI load, the FE remap, the singleton pass and the
  DoF counting. This plan removes that overhead. The prefix stays supported and
  is the statistical oracle in §11.
- **Each estimator is one plugin call that loads, estimates, writes and frees:**
  - creghdfe: `do_full_regression`. It starts at
    `src/creghdfe/creghdfe_regress.c:216` and frees `g_state` at line 1587.
  - cppmlhdfe: `do_ppml_regression`. It starts at
    `src/cppmlhdfe/cppmlhdfe_irls.c:247` and uses `g_ppml_state` (line 227).
  - civreghdfe: `do_iv_regression` (`src/civreghdfe/civreghdfe_impl.c:330`).

  None of them has phases in `ctools_command_begin` (`src/ctools_runtime.c:178-240`).
  cmerge, csplit, crangejoin, cwinsor, cimport and cio already use that
  multi-call mechanism, which §6.3 reuses.
- **The fixed-effect projection is ready for concurrency but runs column by
  column.** `partial_out_columns` (`src/hdfe/hdfe_solver.c:437`) parallelizes over
  columns, and `ctools_hdfe_alloc_buffers` caps the thread count at the number of
  columns (`src/ctools_hdfe_utils.c:714-720`). With the usual 2–10 columns, most
  cores sit idle.
  - The solver keeps all scratch in the `HDFE_State`, indexed by thread, and has
    no static data. Separate states can therefore run concurrently.
  - Levels with zero weight get an inverse count of 0, so a replicate can leave
    absent levels in place.
  - The CG stopping rule is relative to the (weighted) norm of the input column
    (`hdfe_solver.c:375-382`). §6.7 relies on this.
- **RNG.** `src/ctools_rng.h` provides xoshiro256** seeded through SplitMix64.
  - cbsample takes its seed from Stata's `runiform()` as two exact 31-bit halves
    (`build/cbsample.ado:128-129`, parsed by `ctools_parse_seed_option`). It seeds
    each stratum separately (`src/cbsample/cbsample_impl.c:669,700`).
  - The resampling code is inline in `cbsample_main` and cannot be reused as it
    stands. There is no binomial sampler.
- **Estimators call the Stata SPI mid-computation.** The SPI may only be called
  from the Stata thread, so these calls must move out of any code that runs on a
  worker:
  - creghdfe's `store_fixed_effects` reads y and X back with `SF_vdata`
    (`creghdfe_regress.c:54-67`).
  - cppmlhdfe reads scalars and options inside IRLS (`SF_scal_use`, `SF_vdata`,
    `ppml_option`).
  - `civreghdfe_estimate.c` makes about 85 `SF_*`/`ctools_scal_save` calls.

### creghdfe

- **Saved fixed effects.** `absorb(..., savefe)` creates `__hdfe#__`
  (`build/creghdfe.ado:105-111, 354-365`). There is no `newvar=` syntax: the
  shared absorb parser rejects it (`build/_ctools_hdfe_absorb.ado:21-23`).
- **Normalization.** `store_fixed_effects` (`creghdfe_regress.c:26-127`) runs
  forward Kaczmarz sweeps from zero with steepest-descent steps, as in reghdfe's
  `store_alphas`. Every effect after the first then has a zero weighted sum within
  each connected component it forms with the previous effect.
- **savefe tightens the projection tolerance** to 1e-12 for the whole regression
  (`src/creghdfe/creghdfe_input.c:108-110`). This matters for the cost of an FE
  bootstrap.
- **One cluster dimension only.** The ado rejects multiway clustering
  (`build/creghdfe.ado:182-185`).
- **Refits.** When singleton removal changes the factor-variable terms, the ado
  reruns the whole command (`build/creghdfe.ado:410-425`). The bootstrap must run
  only on the final fit.

### cppmlhdfe

- **Fixed-effect outputs.**
  - `savefe` and `absorb(newvar=absvar)` go through `ppml_save_effects`
    (`cppmlhdfe_irls.c:73-156`). It runs CGLS from zero, which gives the
    minimum-norm normalization.
  - `d()` (lines 1832-1840) is z − resid − Xb − mean_d. The IRLS-weighted mean
    mean_d is folded into `_cons`.
- **Clustering.** Multiway clustering is supported, through inclusion–exclusion.
- **Sample screening.** Singleton and all-zero-level screening counts rows, not
  weights (`ppml_select_fe_sample_ex`, `cppmlhdfe_separation.c:20`). The
  separation checks are simplex, ReLU and the mu check inside IRLS.
- **Starting values.** `guess(simple|ols|variable)` sets the initial mu (lines
  782-819). There is no start from η or β. Within IRLS, the previous partialled
  columns are reused.
- **Refit loop.** The ado reruns the plugin in a loop when a check drops rows
  (`build/cppmlhdfe.ado:634-667`).

### civreghdfe

- **No fixed-effect saving.** `absorb()` takes no suboptions
  (`build/civreghdfe.ado:299-302`). `predict, d` is y − xb − resid and needs
  `residuals()`.
- **Estimators.** All six run inside `ivest_compute_2sls`
  (`civreghdfe_estimate.c:759`): 2SLS, LIML, Fuller, k-class, GMM2S (line 42) and
  CUE (line 289).
  - GMM2S and CUE depend on the VCE choice through the moment covariance S
    (`civreghdfe_omega.c`).
  - Diagnostics run on every call (`civreghdfe_estimate.c:1658-2007`): first-stage
    F, KP, Hansen J, C, endogeneity and redundancy tests.
- **Restrictions and reruns.**
  - fweights are refused with tsset kernel VCEs (`build/civreghdfe.ado:596-599`).
  - A second full run, `civrerun`, follows when the plugin drops terms
    (`build/civreghdfe.ado:1024-1091`).
  - `first`, `rf` and `savefirst` rerun creghdfe under `_estimates hold`.

## 4. User-facing design

### 4.1 Syntax

```
creghdfe   depvar indepvars [if] [in] [weight], absorb(absvars [, savefe]) ///
           [vce(bootstrap [, bootopts]) | vce(robust|cluster c) wildbootstrap([wildopts])] ///
           [fevce(bootstrap [, bootopts feopts])] ...
cppmlhdfe  depvar indepvars ..., absorb(absvars [, savefe]) [d(newvar)] ///
           [vce(bootstrap [, bootopts]) | vce(robust|cluster c) wildbootstrap([wildopts])] ///
           [fevce(bootstrap [, bootopts feopts])] ...
civreghdfe depvar (endog = iv) exog ..., absorb(absvars [, savefe]) ///
           [vce(bootstrap [, bootopts]) | vce(robust|cluster c) wildbootstrap([wildopts])] ///
           [fevce(bootstrap [, bootopts feopts])] ...
```

`vce()` keeps its usual meaning: the covariance of e(b). `fevce()` is new and
governs inference on the saved fixed effects. Because there are two options,
each target can be requested on its own. `wildbootstrap()` adds wild-bootstrap
tests on the coefficients on top of an analytic `vce()` (§4.7).

| Mode | Specification | e(V) | Coefficient tests | FE output |
|------|---------------|------|-------------------|-----------|
| Coefficients only (pairs) | `vce(bootstrap, reps(500) cluster(firm))` | pairs bootstrap | normal and bootstrap intervals | point estimates only, if saved |
| Coefficients only (wild) | `vce(cluster state) wildbootstrap(reps(9999))` | analytic cluster | analytic, plus a wild-bootstrap table | point estimates only, if saved |
| Fixed effects only | `vce(cluster firm) fevce(bootstrap, reps(500) cluster(firm))` | analytic cluster | analytic | estimates plus pairs-bootstrap SE |
| Both (pairs) | `vce(bootstrap, reps(500) cluster(firm)) fevce(bootstrap)` | pairs bootstrap | normal and bootstrap intervals | estimates plus pairs-bootstrap SE |
| Both (wild coefficients, pairs FE) | `vce(cluster firm) wildbootstrap(reps(999)) fevce(bootstrap, reps(500) cluster(firm))` | analytic cluster | analytic, plus a wild-bootstrap table | estimates plus pairs-bootstrap SE |

- **Both:** `fevce(bootstrap)` with no design options inherits the design from
  `vce(bootstrap, ...)`, and one set of replicates serves both targets. If the two
  specify different designs (reps, cluster, strata or seed), the command stops
  with r(198): "vce(bootstrap) and fevce(bootstrap) share one set of replicates;
  specify the resampling options once".
- **FE only:** if `vce()` clusters on a variable and `fevce(bootstrap)` has no
  `cluster()`, a note says the fixed effects resample observations while the
  coefficients cluster on that variable.
- **Requirement:** `fevce(bootstrap)` needs a saved fixed effect
  (`absorb(..., savefe)`, cppmlhdfe's `newvar=` terms or `d()`). Without one it
  stops with r(198).
- **Wild and FE together:** `fevce()` always uses the pairs bootstrap.
  - With `wildbootstrap()`, `fevce(bootstrap)` inherits nothing and must give
    its own `reps()`, and `cluster()` if wanted. Otherwise it stops with
    r(198): "fevce() uses the pairs bootstrap; specify its reps() and cluster()".
  - `fevce(wildbootstrap ...)` or `fevce(bootstrap, wild)` stops with r(198):
    "the wild bootstrap applies to coefficients only".

Considered and rejected:

- A single `bootstrap(... target(b|fe|both))` option. It hides the coefficient VCE
  from the `vce()` option every Stata user already reads.
- An absorb suboption (`absorb(..., savefe bootse)`). It would also need a home for
  the design options in FE-only mode.
- `vce(bootstrap, wild ...)` for the wild bootstrap. A wild bootstrap is a test
  and interval procedure, and the e(V) it studentizes is the analytic
  cluster-robust one. Its own covariance of b* (unrestricted) is, in
  expectation, the CV0 sandwich, so it adds nothing as an e(V). Putting it
  under `vce()` would suggest that e(V) changed.

### 4.2 Bootstrap options (bootopts)

One parser, `build/_ctools_hdfe_boot.ado`, serves all three commands.

| Option | Default | Meaning |
|--------|---------|---------|
| `reps(#)` | 50 (Stata's default; see §14) | Number of replicates, at least 2 |
| `seed(#)` | current RNG state | Run `set seed #` before the first draw, as `bootstrap`'s `seed()` does |
| `draws(stata\|fast)` | `stata` (§14) | How Stata's random-number stream becomes replicate draws (§5.8). `fast`: Stata's RNG seeds an in-plugin generator per replicate; reproducible and scalable, but not `bsample`'s draws. `stata`: every draw is made with Stata's generator and reproduces the replicates of `bootstrap … : cmd` for the same seed. |
| `cluster(term)` | none: resample observations | Resample clusters. One term: a variable, `i.var` or `var1#var2`, parsed by `_ctools_hdfe_vce` |
| `strata(varlist)` | none | Resample independently within strata. A cluster that appears in several strata is a separate unit in each, as `bsample` numbers clusters within strata |
| `mse` | off | Center the variance on e(b) instead of the replicate mean |
| `level(#)` | `c(level)` | Level of the posted e(ci_*) intervals |
| `nodots`, `dots(#)` | dots every replicate | Progress display, as in `bootstrap` |
| `saving(filename [, replace every(#)])` | none | Coefficient replicates written as a dataset (§6.9) |

An `idcluster()` option is not needed. §5.2 explains why the multiplicity
representation already treats duplicated clusters correctly for the
coefficients, and how the one estimator where it matters (civreghdfe's
cluster-robust S) is handled.

### 4.3 Fixed-effect options (feopts, `fevce()` only)

| Option | Effect |
|--------|--------|
| (default) | For each saved FE variable `v`, a new double variable with the bootstrap SE of `v`'s level on every e(sample) row. Named `__fese#__` for `__hdfe#__` (not `__hdfe#_se__`: a later `reghdfe` run in the same data, as in the tests, runs `cap drop __hdfe*`), and `v_se` for a user-named `v` (cppmlhdfe `newvar=`). All names are checked with `_ctools_newvars` before anything is created. |
| `sestub(stub)` | Name the SE variables `stub1`, `stub2`, … instead |
| `ci` | Also save percentile bounds `..._lb`, `..._ub` at `level()`. Needs the replicate draws kept in memory; refused above a memory cap (§6.5). Delivered in Phase 5. |
| `nreps` | Also save the number of replicates in which the level was estimable: `__fenrep#__`, or `stub#_nrep` |
| `dse(newvar)` | SE of the summed fixed effect d per row. This is the identified quantity when there are two or more effects (§5.7). The warning in §4.8 recommends it. It is implied for cppmlhdfe's `d()`. |

### 4.4 Results

**Coefficient target.** The command posts what Stata's own `vce(bootstrap)`
posts, so that `estat bootstrap` (normal, percentile and bias-corrected intervals)
works. The list below is what `regress, vce(bootstrap)` posts in StataNow 19
(Appendix A). Stata 14.1 itself was not available to compare.

- **Locals:**
  - e(vce) = `bootstrap`, e(vcetype) = `Bootstrap`, e(prefix) = `bootstrap`;
  - e(rngstate): c(rngstate) after `seed()` and before the first draw;
  - e(chi2type) = `Wald`;
  - e(exp1)…e(expK) = `_b[name]`;
  - e(command): the command without its `vce()`; e(cmdname);
  - e(cluster_bs) and e(strata_bs), when used.
- **Scalars:**
  - e(N_reps), e(N_misreps), e(level), e(bs_version) = 3;
  - e(k_exp) = 0, e(k_eexp) = colsof(e(b)), e(k_extra) = 0, e(k_eq) = 1;
  - e(chi2), e(p).
- **Matrices:**
  - e(b_bs): mean of the replicates;
  - e(reps): complete replicates, the same for every coefficient (§5.5);
  - e(bias) = e(b_bs) − e(b);
  - e(se), e(z0);
  - e(ci_normal), e(ci_percentile), e(ci_bc);
  - e(V_modelbased): the analytic covariance the command would have posted
    without `vce(bootstrap)`.
- **Dropped:** e(df_r) is not posted, so `ereturn display` reports z statistics
  and normal intervals, as `regress, vce(bootstrap)` does. e(F) is replaced by
  e(chi2) (the Wald statistic of the kept regressors) with e(df_m).
- **Kept:** the DoF results that describe the absorbed effects (e(df_a),
  e(dof_table) and so on) are unchanged.
- **Command-specific:** e(boot_status), a 1 × 5 matrix counting replicates by
  outcome: ok, collinear, no observations left, not converged, other.

**FE target.**

- **Locals:** e(fevce) = `bootstrap`; e(fe_se), the SE variables; plus e(fe_lb),
  e(fe_ub) and e(fe_nrep) when requested.
- **Scalars:** e(fe_N_reps).
- **Matrices:** e(fe_bs_table), with one row per effect and columns: levels,
  levels never estimable, levels estimable in every replicate, minimum valid
  replicates, median valid replicates.
- **Normalization warning:** e(fe_norm_warning) is 1 when the warning in §4.8
  was issued, 0 otherwise. e(fe_split_reps) counts the replicates in which a
  connected set split (§5.7).

**Wild-bootstrap results** are in §4.7. With `wildbootstrap()`, e(V), e(vce)
and the analytic test results are those of the analytic `vce()`.

### 4.5 Display

- **Header:**
  - "Wald chi2(k)" and "Prob > chi2" replace F;
  - a "Replications" line is added;
  - "(Replications based on C clusters in var)" is shown, as Stata does.
- **Table:** "Bootstrap Std. Err.", z, and "Normal-based [95% Conf. Interval]".
  `_creghdfe_display` (`build/creghdfe.ado:734-779`) and its counterparts branch
  on e(vce), so replay works.
- **Footnote:** when e(N_misreps) > 0, "Note: N replications were incomplete or
  failed; see e(boot_status)". In FE mode, one line gives the per-effect estimable
  summary from e(fe_bs_table).
- **Two or more absorbed effects with per-level FE SEs:** the normalization
  warning in §4.8.
- **wildbootstrap():** the wild-bootstrap table follows the analytic one
  (§4.7).

### 4.6 Restrictions

| Combination | Behavior |
|-------------|----------|
| A bootstrap `cluster()` with more than one term | r(198): one resampling dimension only. A multiway analytic `vce(cluster a b)` (cppmlhdfe, civreghdfe) stays allowed in FE-only mode. |
| civreghdfe `vce(bootstrap)` with `kernel()`, `bw()`, `dkraay()` or `kiefer` | r(198): resample panels with `cluster(panelvar)` instead |
| civreghdfe `b0()` | r(198) with `vce(bootstrap)` |
| `fevce(bootstrap)` without a saved fixed effect | r(198) |
| pweights | Allowed. The note says this is the naive with-replacement bootstrap, not `svy bootstrap`. |
| Fewer than 2 resampling units in the sample or in a stratum | r(2001) or r(459) |
| Fewer than about 30 clusters | With `vce(bootstrap)`, a warning that the pairs cluster bootstrap is unreliable with few clusters, suggesting `wildbootstrap()` |
| `wildbootstrap()` with `vce(bootstrap)` | r(198): choose one bootstrap method for the coefficients |
| `wildbootstrap()` with `vce(unadjusted)` or no `vce()` | r(198): the wild bootstrap studentizes with `vce(robust)` or `vce(cluster)` |
| `wildbootstrap()` with a multiway `vce(cluster a b)`, or a multi-term `bootcluster()` | r(198): one clustering dimension in the first version (§15) |
| `wildbootstrap()` with an observation-level wild bootstrap (`vce(robust)`, no `bootcluster()`) and fweights | r(198): an fweighted row stands for several observations, and each would need its own draw |
| civreghdfe `wildbootstrap()` with `gmm2s` or `cue` | r(198): WRE covers 2SLS, LIML, Fuller and k-class |
| civreghdfe `wildbootstrap()` with `kernel()`, `bw()`, `dkraay()` or `kiefer` | r(198) |

### 4.7 Wild bootstrap (coefficients only)

```
creghdfe   y x1 x2, absorb(firm year) vce(cluster state) wildbootstrap(reps(9999) weighttype(webb))
civreghdfe y (x = z1 z2) w, absorb(firm year) vce(cluster state) wildbootstrap
cppmlhdfe  trade fta, absorb(exp#year imp#year exp#imp) vce(cluster pair) wildbootstrap(coefs(fta) ci)
```

`wildbootstrap()` runs after the analytic fit and leaves e(b), e(V) and the
analytic table unchanged. For each tested coefficient it adds:

- a bootstrap p-value for H0: β_k = 0;
- a confidence interval found by inverting that test (§5.9).

The bootstrap clusters are the `vce(cluster)` variable, or `bootcluster()` when
given. With `vce(robust)` the wild bootstrap draws one weight per observation.

**wildopts**

| Option | Default | Meaning |
|--------|---------|---------|
| `reps(#)` | 999 (as boottest) | Replications. Reduced to 2^C when the sign vectors are enumerated (below). |
| `seed(#)` | current RNG state | `set seed #` first. The weights are always drawn with Stata's generator (§5.8). |
| `weighttype(name)` | `rademacher` | `rademacher`, `mammen`, `webb` (six-point), `normal` or `gamma` |
| `nonull` | off: null imposed (WCR, WRE) | Unrestricted (WCU) bootstrap |
| `coefs(varlist)` | every kept regressor except `_cons` | Coefficients to test. Fewer coefficients cost less when the effects are not nested (§5.9). |
| `ptype(name)` | `symmetric` (as boottest) | `symmetric`, `equaltail`, `lower` or `upper` |
| `noenumerate` | off: enumerate when possible, as boottest | Always draw, even when all 2^C Rademacher sign vectors fit in reps (see below) |
| `fedgp(full\|partialled)` | `full` (as boottest) | How effects not nested in the bootstrap clusters enter the bootstrap DGP (§5.9). `partialled` matches Stata's `wildbootstrap areg` and is faster. |
| `noci` | intervals on, except in cppmlhdfe | Skip test inversion. cppmlhdfe intervals are opt-in (`ci`), because each trial value needs a restricted PPML fit (§8.4). |
| `level(#)` | `c(level)` | Level of the intervals |
| `bootcluster(term)` | the `vce()` cluster | Cluster at which weights are drawn: one term, coarser or different from the `vce()` cluster, as in boottest |
| `nodots`, `dots(#)` | dots | Progress |

**Reference.** The wild bootstrap's defaults and conventions follow boottest
4.5.5, the established tool after reghdfe and ivreghdfe (Appendix A, checks
12–16). A result can therefore be reproduced with boottest wherever boottest
runs: one absorbed effect, and no PPML.

**Enumeration and few clusters.** With Rademacher weights and 2^C ≤ reps, the
command uses all 2^C sign vectors instead of drawing, as boottest does. It sets
e(wb_reps) to 2^C and e(wb_enumerated) to 1, and the p-values are then exact
for the bootstrap distribution. `noenumerate` draws instead, as Stata's
`wildbootstrap` always does. With C ≤ 11 clusters and Rademacher weights, a
note suggests `weighttype(webb)`, because Rademacher then allows only 2^C
distinct draws (MacKinnon and Webb).

**Results.** All results use the prefix `wb_`:

- **Locals:**
  - e(wb_method): `wcr`, `wcu`, `wre` or `score`;
  - e(wb_weighttype), e(wb_ptype);
  - e(wb_bootclust), e(wb_coefs);
  - e(wb_fe_path): `nested` or `projected` (§5.9).
- **Scalars:** e(wb_reps), e(wb_misreps), e(wb_enumerated), e(wb_N_clust),
  e(wb_level).
- **Matrices:**
  - e(wb_t): the observed t statistics;
  - e(wb_p): bootstrap p-values (1 × k);
  - e(wb_ci): test-inversion intervals (k × 2), with missing bounds when the
    confidence set is unbounded on that side.

**Display.** After the analytic table, and on replay:

```
Wild cluster bootstrap: restricted (WCR), Rademacher weights
Replications = 9,999           Bootstrap clusters = 42 (state)
---------------------------------------------------------------
           y |      Coef.          t   P>|t|   [95% Conf. Interval]
-------------+-------------------------------------------------
          x1 |   .1234567       2.31   0.0410   .0061234   .2398765
          x2 |  -.0456789      -0.88   0.3912  -.1543210   .0634567
---------------------------------------------------------------
P-values are symmetric bootstrap-t; intervals invert the bootstrap test.
```

The numbers above are placeholders; only the layout is specified. If a
confidence set is unbounded or not an interval, which can happen with very
few clusters, a footnote names the coefficient and reports the open side.

### 4.8 Fixed-effect normalization warning

**When.** The warning is printed whenever all of these hold:

- `fevce(bootstrap)` writes per-level fixed-effect SE variables (the default
  output of §4.3);
- the model absorbs two or more effects with intercepts. In cppmlhdfe, effects
  that only carry slopes do not count.

**Where.** At estimation it appears after the coefficient table, as a warning
in the style of reghdfe's notes. On replay it appears in the footnote. It is
not suppressed by `nofootnote` at estimation, because it qualifies the new
variables.

```
warning: with 2 absorbed effects (firm year), individual fixed effects are
         identified only up to a normalization within each connected set of
         levels. The bootstrap standard errors in __fese1__ __fese2__
         depend on that normalization and can be large where the sets are thinly
         connected. The summed fixed effect d does not depend on it: specify
         fevce(bootstrap, dse(newvar)) for its standard errors.
```

**Variants.**

- When `dse()` was given, the last sentence instead reads "The summed fixed
  effect d does not depend on it; its standard errors are in `newvar`."
- When connected sets split in some replicates (§5.7), a second line is added:
  "In R of N_reps replicates a connected set split; levels in a split set were
  treated as not estimable in those replicates (see e(fe_bs_table))."

**Storage and documentation.** e(fe_norm_warning) is set to 1 and
e(fe_split_reps) to R. The help files carry the same caveat in a "Remarks on
fixed-effect standard errors" section. With one absorbed effect the warning is
not shown, because the per-level SEs are then free of this problem (§5.7).

## 5. Statistical design

### 5.1 Resampling

- **Units.** A resampling unit is a row of the final estimation sample, or a
  cluster when `cluster()` is given.
- **Draws.** Replicate r draws, independently in each stratum s, `n_s` units with
  replacement from that stratum's `n_s` units. This is a Multinomial(n_s, 1/n_s)
  count `m_u` per unit.
- **Row multiplicity.** A row's multiplicity is its unit's count.

**When resampling happens.** Resampling acts on the full-sample estimation
sample: after `if`/`in`, missing values, singletons, separation and the
command's other exclusions. The `bootstrap:` prefix instead resamples the rows
it is given, and the command then screens each replicate again. Phase 0 found
that the prefix first restricts those rows to e(sample), so the frames agree
(Appendix A, check 2). Here, replicate
screening (§5.4) re-applies the command's rules to the multiplicity-weighted
replicate. Each replicate therefore sees the exclusions the command would apply
to the expanded data, and singletons, separation and collinearity are judged
in each replicate as in the prefix.

**Strata** are evaluated on e(sample). A cluster that crosses strata is a separate unit in each stratum, as in `bsample` (Phase 2b found that parity needs this).

### 5.2 Multiplicities as weights

Each replicate is the full-sample model re-fit with row weights ω_i = m_i · w_i,
where w_i is the user's weight (1 if none). Rows with m_i = 0 are compacted away;
about 63% of rows remain in an observation bootstrap.

**Exact equivalences.** For coefficients and fixed effects, these weighted fits
give exactly what the command gives on the physically expanded data, `expand m`:

- **OLS with fixed effects.** Partialling and the normal equations are weighted
  sums, and a row repeated m times contributes m times.
- **IV.** 2SLS, LIML, Fuller and k-class depend on the data only through weighted
  cross-product matrices.
- **PPML.** The score and IRLS steps are weighted sums over rows.
- **aweights and pweights.** These are relative weights. Duplicating a row with
  weight w equals one row with weight m·w, and normalizing to the replicate's
  row count does not change b.

**Duplicated clusters.** Duplicating a cluster with `idcluster()` and keeping its
fixed-effect level gives the same coefficients as weighting it m times:

- **Effect nested in the cluster** (cluster on firm, absorb firm): each copy is
  demeaned by its own copy of the level mean, which is the same mean.
- **Crossed effects** (years with firm clusters): the copies share levels under
  either representation.

So the coefficient draws match `bootstrap, cluster(c) idcluster(nc): cmd` with
the original absorb() variables.

**The exception: civreghdfe's cluster-robust moment covariance S.** GMM2S and
CUE weight with S, and so do diagnostics, which are not bootstrapped. Under
`idcluster()` each copy of a cluster contributes its own outer product: the
contribution is m_c · s_c s_c′. Weighting rows by m_c instead gives m_c² · s_c
s_c′. §9 adds a "copy-aware" S builder that takes unit multiplicities. The robust
(non-cluster) S is already right with fweight semantics: each row contributes
m_i · (w_i e_i)² z_i z_i′ for aweights and pweights, m_i e_i² z_i z_i′
unweighted. The builder must take the multiplicity and the analytic weight
separately, because today it takes one weight vector of one type.

### 5.3 User weights

| Weight | Resampling | Replicate weight |
|--------|-----------|------------------|
| none | rows or clusters | m_i |
| aweight | rows or clusters | m_i · aw_i (renormalized as the command already does) |
| pweight | rows or clusters | m_i · pw_i |
| fweight | **row resampling:** the units are the Σ fw underlying observations. Draw Multinomial(Σ fw, fw_i / Σ fw) counts per row. **Cluster resampling:** clusters are drawn, and the weight is m_c · fw_i. | the drawn counts, or m_c · fw_i |

Treating an fweighted row as fw_i observations is the correct bootstrap of
frequency-weighted data. Resampling rows would understate the variance. The
row-level fweight draw uses conditional binomials, which is O(N) per replicate
whatever Σ fw is. That needs a `ctools_rng_binomial` (BTPE, with inversion for
small n·p) in `src/ctools_rng.h`. Stata's `vce(bootstrap)` refuses weights
altogether (Appendix A, check 5), so weights here are an extension that
`COMPATIBILITY.md` documents (§14, question 4).

### 5.4 Replicate sample rules

Each replicate re-applies the command's own sample rules to the multiplicity
data, in the units the command would see after `expand`:

- **Singletons.** A level whose total multiplicity is 1 is a singleton. This is
  the rule creghdfe already applies to fweights (`ctools_remove_singletons` with
  fweights). `keepsingletons` turns it off.
  - cppmlhdfe screens by row counts today, so it gets a multiplicity-aware
    variant of `ppml_select_fe_sample_ex`. The full-sample path is unchanged.
- **Separation (cppmlhdfe).**
  - A level whose drawn rows all have y = 0 is separated and dropped.
  - The simplex and ReLU checks, if enabled for the fit, run on the replicate
    with the same settings, warm-started where possible (§8).
- **Collinearity.**
  - The collinear-with-FE rule and the `invsym` rule run on the replicate as on
    the full sample.
  - A regressor that is kept in the full sample but omitted in a replicate makes
    that coefficient missing in that replicate, never zero.
  - A regressor omitted in the full sample stays omitted, and its column of
    e(V) stays zero.
- **Depvar explained by the fixed effects.** The replicate is recorded and
  flagged in e(boot_status). It is not dropped, because it is a legitimate
  value.

### 5.5 Failed and incomplete replicates

| Outcome | Coefficients | Fixed effects |
|---------|--------------|---------------|
| All rows singletons or separated | replicate missing | replicate missing |
| Projection or IRLS not converged (r(430) in a full fit) | replicate missing | replicate missing |
| Omitted-coefficient pattern differs from the full sample (a kept coefficient omitted, or an omitted one estimable) | replicate missing: every coefficient, as Stata's prefix posts it | missing (2c; see "Departures" under Phase 2c) |
| FE not converged (FE recovery hits iterate()) | unaffected | missing |

e(N_misreps) counts replicates with missing coefficients, and e(V) uses
complete replicates only (casewise), so e(reps) is the same for every
coefficient. Both are Stata's conventions, confirmed in Phase 0. Stata then
prints "Note: One or more parameters could not be estimated in N bootstrap
replicates; standard-error estimates include only complete replications."
The commands print the same note.

If fewer than two complete replicates remain, the command stops with r(2000),
"insufficient observations to compute bootstrap standard errors", and posts
nothing.

### 5.6 Summary statistics

R is the number of complete replicates, b̂ is e(b), and b^r is replicate r.

| Statistic | Formula |
|-----------|---------|
| b̄ | mean of b^r |
| e(V), default | V = (R−1)⁻¹ Σ (b^r − b̄)(b^r − b̄)′ |
| e(V), `mse` | V = R⁻¹ Σ (b^r − b̂)(b^r − b̂)′ |
| Normal CI | b̂ ± z_{1−α/2} · se |
| Percentile CI | empirical α/2 and 1−α/2 quantiles, using Stata's percentile definition (match `_pctile`) |
| BC CI | z0 = Φ⁻¹(#{b^r ≤ b̂}/R), missing when the count is 0 or R; percentiles at Φ(2z0 ± z_{1−α/2}) (Stata's `_bs_sum`; ≤, not <) |

Stata posts its replicates as floats unless `saving(…, double)` is given, and
computes e(V) and the intervals from the posted values (Appendix A). The
commands summarize in double precision. They therefore agree with Stata's
default output to about 1e-7 relative, and with `saving(…, double)` output to
rounding.

e(chi2) = b̂′ V⁻ b̂ over the kept regressors, excluding `_cons`. It is computed
with `invsym` and is missing when V is rank deficient, as the analytic F is now.

### 5.7 Fixed-effect targets

**Estimability.** In replicate r, a level of effect g is estimable when both
hold:

1. its replicate weight is positive after replicate screening; and
2. for G ≥ 2, the connected component that the command's normalization is
   defined over has not split. A component splits when the levels that appear in
   the replicate's component did not all come from one full-sample component.

A level that is not estimable is recorded as missing (NaN) for that replicate,
never 0. This matters because both recovery routines start from zero and leave
absent levels at zero. The SE of a level uses its estimable replicates only.
`nreps` and e(fe_bs_table) report how many there were.

**Normalization.** Each replicate is normalized by the same deterministic rule as
the point estimate, applied to the replicate's own data:

- creghdfe and civreghdfe: reghdfe's convention (§3);
- cppmlhdfe: minimum norm.

This bootstraps the estimator as it is defined. With one absorbed effect it is
also the natural target. The level effect is then (level mean of y − Xb) minus the
overall weighted mean, and the bootstrap correctly includes the sampling noise
of that overall mean.

**The caveat with two or more effects.** Individual effects are identified only
up to constants within components. The per-level SE therefore depends on the
normalization, and it can be large in sparsely connected designs (workers and
firms with little mobility). Two outputs help:

- The SE of d, the summed effect per row (`dse()`), is free of this problem. In
  replicate r, d is computed for every full-sample row whose levels are all
  estimable and lie in one replicate component, including rows that were not
  drawn. Its SE is the identified quantity.
- The command reports this caveat as a warning, pointing to `dse()`, whenever
  it writes per-level SEs with G ≥ 2 (§4.8). The help files say the same.

Anchoring each replicate to a full-sample normalization is left to later work
(§15).

**Accumulation.** Replicate draws map back to full-sample level indices. §6.2
shows why no renumbering is needed. They are accumulated with Welford updates, in
replicate order. Percentile bounds (`ci`) keep the draws as float32
(levels × replicates).

### 5.8 Seeds and Stata's random-number generator

Stata's RNG is the only source of randomness. All three commands treat it as
Stata's `bootstrap` prefix does:

- `set seed #` before the command, or the `seed(#)` option (which runs
  `set seed #` first), makes results reproducible.
- `set rng` (mt64, kiss32, …), `set rngstream` and `version` control are
  honored. Every random number is taken from Stata's current generator and
  stream, through `runiform()`, `runiformint()` or `rnormal()` in Mata.
- The command never reseeds Stata or restores an earlier RNG state. Like
  `bootstrap`, it leaves Stata's RNG advanced past the draws it took. A
  following command therefore gets fresh numbers, and a `set seed` before the
  whole do-file reproduces everything.
- e(rngstate) holds c(rngstate) from before the first draw, as `bootstrap`
  posts it. e(seed) holds `seed()` when given. e(draws) records the draw mode.
- No random number is taken when no bootstrap is requested. The analytic
  commands leave c(rngstate) untouched, as they do now.

Two draw modes share these rules. They differ in *how* Stata's stream becomes
replicate draws.

**`draws(fast)`** (the default until Phase 5, now opt-in).

1. The ado takes the master seed from Stata's stream as two exact 31-bit halves
   with `runiform()`, as cbsample does (`build/cbsample.ado:128-129`). That is
   exactly two draws.
2. Replicate r gets its own in-plugin generator: xoshiro256** seeded through
   SplitMix64 from (master seed, r). Within a replicate the strata are drawn in a
   fixed order.

The multiplicities are a pure function of Stata's RNG state at entry and of r.
`set seed` therefore reproduces them exactly, and they do not depend on
`threads()`, the thread schedule or the chunk boundaries. They are not the
draws that Stata's `bootstrap` would make for the same seed. This is the only
mode that scales to observation-level resampling of very large data, because no
per-row draw crosses the SPI.

**`draws(stata)`: Stata's own bootstrap draws.**

- **What.** Each replicate's sample is drawn with Stata's generator,
  reproducing `bsample` exactly: the same units, the same order of uniform
  draws, the same mapping from draws to units within strata, and one call per
  replicate in order, as the `bootstrap` prefix makes them. With the same seed,
  `creghdfe ..., vce(bootstrap, reps(R) seed(s) cluster(c) draws(stata))` then
  uses the replicates that
  `bootstrap, reps(R) seed(s) cluster(c) idcluster(nc): reghdfe ...` uses. It
  leaves c(rngstate) where that command would leave it.
- **What Phase 0 established** (Appendix A, all verified exactly):
  - The prefix resamples **e(sample)** of the full fit, so a command's own drops
    (reghdfe's singletons) are outside the frame. That is the frame of §5.1.
  - Its replicates are R successive `bsample` calls after `set seed`, with no
    other draws in between. It ends with the same c(rngstate) as those calls.
  - Units are rows in data order, or clusters numbered in sorted cluster order
    within strata sorted the same way (`sort strata cluster, stable`).
  - Each replicate takes one uniform per unit for the draw, then one per unit
    for a tie-break, stratum by stratum. The counts are attached to units
    through a sort of (draw, tie-break) and a fixed permutation rule.
    `idcluster()` takes no extra draws.
  - `validation/support/bootstrap_reference.do` reproduces this in about 100
    lines of Mata, with no internal built-ins involved, so no fallback is
    needed.
- **How.**
  - Stata draws, the plugin maps. Per chunk, the ado draws the 2·U uniforms per
    replicate with Mata's `runiform()` on Stata's stream, in `bsample`'s order,
    and passes them in temporary variables.
  - `boot_chunk` applies `bsample`'s draw-to-unit mapping in C (a sort per
    stratum) to get the multiplicities.
  - The units are numbered in `bsample`'s order. The frame is e(sample), which
    the retained model already holds.
- **Cost.** One Mata `runiform()` call and one SPI transfer of 2·U doubles per
  replicate.
  - Cluster bootstraps: negligible (U = C).
  - Observation bootstraps: about two column loads per replicate. That is the
    main cost of choosing `draws(stata)` on very large data, and the reason
    `draws(fast)` exists.
- **Rounding.** The prefix posts replicates as floats by default (§5.6), so
  `draws(stata)` matches its replicates to rounding and its default e(V) to
  about 1e-7.

**Wild bootstrap: always from Stata's generator.** The wild weights are only
C × B values. The ado always draws them in Mata from Stata's stream, in chunks,
and passes them to the plugin the same way (temporary variables, C rows per
replicate). There is no fast mode for the wild bootstrap.

- **Which generator calls.** Each weight type reproduces boottest's
  `MakeWildWeights` (`boottest.mata`), so that the same `seed()` gives the same
  weights as boottest.
  - boottest draws one clusters × (reps + 1) matrix in a single call. Its first
    column becomes the original sample (all ones), and the remaining columns
    are the replicates. It splits the matrix into pieces only when it would
    exceed its memory cap, which the ado must mirror.
  - Rademacher: (u < ½)·(−2) + 1. Webb: √(2·⌈3u₁⌉)·((u₂ ≥ ½) − ½), with the two
    uniform matrices drawn in that order. Mammen: `rdiscrete()` over its two
    points. Normal: `rnormal()`. Gamma: `rgamma(4, ½)` − 2.
  - The exact expressions (including the scaling boottest applies outside WRE)
    are copied from boottest's behavior, not its code. Phase 1w confirms them
    by exact p-value comparison with boottest on the same `seed()` and
    non-enumerated designs. That comparison is the gate for this mapping.
  - Stata's `wildbootstrap` draws differently (Appendix A, check 16), so the
    same seed will not reproduce its results. `fedgp(partialled)`
    `noenumerate` `ptype(equaltail)` `reps(1000)` reproduces its conventions,
    but not its draws.
- **Enumeration.** In enumeration mode (§4.7) no random number is taken: the
  sign vectors are generated in the plugin in boottest's order (its
  `count_binary`), and the RNG state is unchanged. The help file notes this.

**Validation.** `draws(stata)` gives exact comparison with Stata's own commands
for free: the same seed must give the same replicates. §11.1 uses it alongside
the hidden `_bsdraws()` option.

### 5.9 Wild bootstrap (coefficients only)

No rows are resampled. The sample, singletons, collinearity pattern and
fixed-effect nesting are those of the full fit, and a replicate can only fail
numerically: se* = 0 gives a missing t*, counted in e(wb_misreps).

**Notation.**

- W: the full-sample weights.
- M_D: the W-weighted fixed-effect annihilator.
- X̃ = M_D X and ỹ = M_D y: the partialled data, which the fit retains.
- A = (X̃′WX̃)⁻¹.
- β̂ and û = ỹ − X̃β̂: the estimates and the partialled residuals.
- Bootstrap clusters c = 1..C, with weights v_c that have mean 0 and variance 1.

**Weight distributions.**

| Type | Values |
|------|--------|
| Rademacher | ±1, each with probability ½ |
| Mammen | −(√5−1)/2 with probability (√5+1)/(2√5), else (√5+1)/2 |
| Webb | ±√½, ±1, ±√(3/2), each with probability ⅙ |
| normal | N(0, 1) |
| gamma | Gamma(4, ½) − 2 |

#### Linear models (creghdfe): WCR and WCU

For H0: β_k = β0:

- **Restricted residuals by FWL.** Let r_k = M_{−k} x̃_k, the W-weighted residual
  of x̃_k on the other partialled regressors. The restricted residual is
  ũ(β0) = û + (β̂_k − β0) r_k, and the restricted estimate β̃(β0) has β̃_k = β0.
  WCU (`nonull`) uses ũ = û and β̃ = β̂.
- **Bootstrap data.** y*_i = x_i′β̃ + d̃_i + ũ_i v_{c(i)}. Here d̃ is the
  restricted fixed-effect component, so the bootstrap DGP is the full model with
  its fixed effects.
- **Bootstrap statistic.** The command's own estimator is applied to y*: the same
  projection, the same weights and the same analytic cluster-robust V. That gives
  t*_k = (β*_k − β0) / se*_k. The analytic small-sample factor multiplies both t
  and t*, so it cancels in every p-value and interval.
- **The two identities used everywhere:**
  - β* − β̃ = A X̃′W(ũ ⊙ v). This needs no re-partialling, because X̃′W M_D = X̃′W.
  - û* = M_D(ũ ⊙ v) − X̃(β* − β̃). Then
    se*_k² ∝ Σ_g (a_k′ X̃_g′W_g û*_g)² with a_k = A e_k, where g runs over the
    `vce()` clusters.

**Nested path.** Use this path when every absorbed effect is nested within the
bootstrap clusters (e(wb_fe_path) = `nested`). Then M_D(ũ ⊙ v) = ũ ⊙ v: each
level's weighted residual sum is zero, and all its rows share one v_c.

1. Form aggregates on the intersection of the `vce()` clusters and the bootstrap
   clusters (`ctools_intersect_cluster_ids`):
   - S_g = X̃_g′W_g ũ_g, which is linear in β0 through the û and r_k parts;
   - H_g = X̃_g′W_g X̃_g.
2. The statistic then needs only C-sized sums. The setup is O(N·K²); each
   replicate costs O(C·K) per tested coefficient, whatever N is. This is the
   algorithm of Roodman, MacKinnon, Nielsen and Webb (2019). It covers the usual
   `absorb(firm) vce(cluster firm)` and `absorb(county) vce(cluster state)`
   designs.

**Effects not nested in the clusters: two conventions.** Take year effects with
firm clusters. Phase 0 found that the two reference tools disagree here
(Appendix A, checks 9 and 13):

- boottest after reghdfe bootstraps the **full model**: it gives exactly what
  `regress` with dummies gives.
- Stata's `wildbootstrap areg` bootstraps the **partialled** regression: it
  treats ỹ = X̃β + u as the model and never projects the perturbed residuals.

The commands follow boottest:

- **Default, `fedgp(full)`:** the bootstrap of the full model with its fixed
  effects, as the command would estimate it on y*. It matches boottest and
  `regress … i.fe`. It uses the projected path below, and e(wb_fe_path) =
  `projected` when an effect is not nested.
- **Option, `fedgp(partialled)`:** Stata's `wildbootstrap` convention.
  û* = ũ ⊙ v − X̃(β* − β̃), so every design runs on the aggregate path above at
  O(C·K) per replicate. e(wb_fe_path) = `partialled` when an effect is not
  nested. When a design would be slow on the projected path (large N, many
  replicates, many tested coefficients), a note mentions this option.

When every effect is nested the two conventions coincide, and the nested path is
used either way.

**Projected path (`fedgp(full)` with an effect that is not nested).** M_D(ũ ⊙ v)
is needed per replicate. By
linearity in β0:

M_D(ũ(β0) ⊙ v) = M_D(û ⊙ v) + (β̂_k − β0) · M_D(r_k ⊙ v).

- Each replicate needs one projection of û ⊙ v plus one per tested coefficient.
  That is k_tested + 1 columns, through the retained full-sample `HDFE_State`
  (same D, same W; nothing is rebuilt).
- The columns are then aggregated to clusters, and the statistic is evaluated as
  on the nested path.
- With G = 1 the projection is a demean. `coefs()` limits the cost.

**p-values.** B is the number of feasible replicates and t is the observed
statistic. These are boottest's formulas (`boottest.mata`, `getp`). Stata's
`wildbootstrap` uses the same ones. Both have strict inequalities and
denominator B.

| ptype | p |
|-------|---|
| `symmetric` (default) | #{\|t*_b\| > \|t\|} / B |
| `equaltail` | 2 · min(#{t*_b < t}, #{t*_b > t}) / B |
| `lower` | #{t*_b < t} / B |
| `upper` | #{t*_b > t} / B |

**Intervals by test inversion (WCR).**

- The interval is the set of β0 with p(β0) > α.
- For each replicate, the numerator of t*_b(β0) is linear in β0 and the
  denominator is the square root of a quadratic in β0. The engine therefore keeps
  five numbers per (replicate, coefficient), and evaluates p(β0) exactly in
  O(B) for any β0.
- Bounds:
  1. Bracket the bounds from a grid around the analytic interval, widened until
     p < α.
  2. Refine by bisection to 1e-6 · se.
  3. Report a side as open when no bound is found within ±1000 · se.
- With WCU the bootstrap distribution does not depend on β0, and the interval is
  the percentile-t interval.

**User weights.** W appears in every formula, and the weight v is drawn per
bootstrap cluster. Observation-level wild bootstrap with fweights is refused
(§4.6).

#### IV (civreghdfe): WRE

This is the wild restricted efficient bootstrap of Davidson and MacKinnon
(2010), for 2SLS, LIML, Fuller and k-class. For H0 on one coefficient, with
value β0:

1. Restricted structural fit on the partialled data, giving ũ(β0).
2. Efficient restricted first stage: regress Ỹ₂ on Z̃ (the partialled exogenous
   variables and instruments) and ũ(β0), and set ṽ(β0) = Ỹ₂ − Z̃π̃.
3. Generate Y₂* = Z̃π̃ + ṽ ⊙ v and y* from the restricted structural equation
   plus ũ ⊙ v, with the **same** v.
4. Re-estimate with the user's estimator and analytic VCE to get t*.

**Fixed effects and cost.**

- ũ(β0) and ṽ(β0) are combinations of a small basis that does not depend on
  β0. ũ is linear in β0, and ṽ adds the projections of that basis on Z̃ with
  scalar coefficients that are rational in β0.
- **Nested path:** everything reduces to cluster aggregates of the cross
  products of Z̃, X̃, ỹ and the basis.
- **Projected path:** each basis vector ⊙ v is partialled once per replicate,
  which is a few columns per endogenous regressor.
- Either way, any β0 is then evaluated with small dense algebra: a (K₂+1)
  eigenproblem per replicate for LIML and Fuller. That makes test inversion
  affordable.
- Phase 0 pins the hypotheses boottest's WRE supports: endogenous and exogenous
  coefficients, and more than one endogenous regressor. civreghdfe supports the
  same set and refuses the rest with r(198).

**Scope.** The diagnostics (first-stage F, KP, J) are unchanged and analytic.
GMM2S and CUE are excluded (§4.6, §15).

#### PPML (cppmlhdfe): score bootstrap

This is the Kline and Santos (2012) score bootstrap, as boottest implements it
for ML models.

- **Restricted fit.** The fit under H0: β_k = β0 is a PPML fit with x_k·β0 as an
  offset and x_k dropped. It is warm-started from the unrestricted fit (§6.7).
- **Statistic.** At the restricted μ̃:
  - the regressors are partialled once under the IRLS weights μ̃ (and W);
  - the efficient cluster scores are g_c = r̃_k,c′ W_c (y_c − μ̃_c), where r̃_k is
    x_k partialled on the effects and the other regressors;
  - each replicate perturbs the scores, g*_c = v_c g_c, and studentizes them with
    their cluster-robust variance.

  A replicate costs O(C) after the restricted fit, and nothing is re-estimated or
  re-projected per replicate. Phase 0 pins the observed statistic's exact form
  against boottest.
- **Cost.** p-values need one restricted fit per tested coefficient. Intervals
  need restricted fits at the trial β0 values (typically 10–20 per coefficient,
  each warm-started from the previous one), so they are opt-in (`ci`).

## 6. Architecture

### 6.1 Overview

```
 ado: parse ──► plugin "<cmd> full_regression retain" ──► refit check (ado)
                  │ loads once, fits, writes outputs,      │ refit? → state reset,
                  │ keeps compact state (phase 1)          │ rerun as today
                  ▼                                        ▼
            retained model ◄────────── plugin "<cmd> boot_chunk first(r0) count(n) seedhi() seedlo()"
            (Σ ≈ 2·N·K doubles)          │ parallel replicate fits; SPI-free; draws ← (seed, r)
                  │                      │ returns chunk statuses (+ b matrix for saving())
                  │                      ▼
                  │                ado: dots, Break check, saving() posts, next chunk …
                  ▼
            plugin "<cmd> boot_finish"  ── writes FE SE variables (Stata thread),
                  │                         V, b_bs, se, bias, z0, CIs → matrices
                  ▼
            state freed (phase 0) ──► ado: ereturn post / repost, display
```

### 6.2 Shared engine: `src/hdfe/hdfe_boot.{c,h}`

The engine is new code. It owns resampling, scheduling and accumulation. The
estimator owns the replicate fit. A sketch of the header:

```c
/* Resampling design over the N rows of the full-sample estimation sample
 * (compact order). unit[i] in 0..U-1 (NULL: each row is a unit);
 * unit_stratum[u] in 0..S-1 (NULL: one stratum); row_freq: fweights for
 * row resampling (NULL otherwise). All arrays owned by the caller. */
typedef struct {
    ST_int N, U, S;
    const ST_int *unit, *unit_stratum;
    const ST_double *row_freq;
    uint64_t seed;
} hdfe_boot_design;

/* Multiplicities of replicate r: row_mult[N], unit_mult[U]. A pure function
 * of (design, r). Thread safe; no allocation. */
ST_retcode hdfe_boot_draw(const hdfe_boot_design *d, ST_int r,
                          ST_int *row_mult, ST_int *unit_mult);

typedef struct {
    ST_double *b;          /* K_b slots (kept regressors + _cons), filled by fit */
    unsigned char *b_ok;   /* K_b flags */
    ST_double *fe;         /* L_total slots in full-sample level order, NaN = not
                              estimable; NULL when the FE target is off */
    ST_double *d;          /* N slots for dse(), NaN where not identified; or NULL */
    ST_int status;         /* HDFE_BOOT_OK, _COLLINEAR, _EMPTY, _NOCONV, _OTHER */
    ST_int iterations;
} hdfe_boot_rep;

/* The estimator's side. fit() runs on a worker thread: it must not call the
 * SPI, print, or touch globals; all inputs come from model, scratch from ws. */
typedef struct {
    size_t     (*workspace_bytes)(const void *model, ST_int inner_threads);
    void      *(*workspace_new)(const void *model, ST_int inner_threads);
    void       (*workspace_free)(void *ws);
    ST_retcode (*fit)(const void *model, void *ws, const ST_int *row_mult,
                      const ST_int *unit_mult, hdfe_boot_rep *out);
} hdfe_boot_model_ops;

/* Runs replicates [first, first+count), folds results into acc in replicate
 * order. Chooses outer (replicate) vs inner (column) parallelism from threads
 * and mem_budget (§6.5). */
ST_retcode hdfe_boot_run(const hdfe_boot_design *, const hdfe_boot_model_ops *,
                         const void *model, ST_int first, ST_int count,
                         ST_int threads, size_t mem_budget, hdfe_boot_acc *acc);

/* V, b_bs, se, bias, z0, normal/percentile/BC intervals; FE SE and bounds. */
ST_retcode hdfe_boot_summarize(const hdfe_boot_acc *, const ST_double *b_hat,
                               int mse, ST_double level, hdfe_boot_summary *);
```

**Draw sources.** `hdfe_boot_draw` is the `draws(fast)` source. Under
`draws(stata)`, and under the hidden `_bsdraws()` validation option, the chunk's
counts arrive from the ado instead. They come as temporary variables, U values
per replicate, read on the Stata thread in `boot_chunk` before the workers
start. The workers see the same `row_mult` and `unit_mult` arrays either way,
so the replicate fits do not know which source was used.

`hdfe_boot_acc` holds:

- an R × K_b replicate matrix with flags (small; kept whole);
- per-level FE Welford accumulators (count, mean, M2) over L_total = Σ L_g;
- per-row d accumulators when `dse()` is on;
- optional float32 FE draws for `ci`;
- status counts.

**Shared building blocks to add.**

- **`hdfe_build_factors_on_levels`.** A replicate builds factors on the
  full-sample level codes with zero-weight levels left in place:
  `num_levels = L_g`, `weighted_counts` from ω, `inv_weighted_counts = 0` for
  absent levels. No renumbering means a replicate's α arrays are already in
  full-sample indexing. `hdfe_build_factors` renumbers by first appearance, so
  this needs a sibling function in `hdfe_sample.c` that skips the renumbering.
- **`hdfe_fe_recover.{c,h}`.** An SPI-free form of creghdfe's
  `store_fixed_effects`, taking the FE component d as an in-memory array and
  returning α arrays. creghdfe's write-back becomes a thin SPI wrapper around it.
  civreghdfe reuses it. cppmlhdfe gets the same split for `ppml_save_effects`.
- **`hdfe_components`.** Union-find over the levels of the effect pairs a
  normalization uses, with a full-sample component id per level and a split test
  per replicate. It reuses `ctools_UnionFind` from `ctools_hdfe_utils.c`.
- **`ctools_rng_binomial`** in `ctools_rng.h`, for fweights (§5.3).

### 6.3 Retained state and the phase protocol

**Retained model.** The full-sample call keeps a compact model when the ado asks
for it (a `retain` token). The model is a plain struct owned by the command's
module, holding:

- the raw and partialled columns;
- the user weights as given;
- the factor levels and L_g;
- the unit and stratum ids;
- the mask and `obs_map` (for write-back);
- b̂, keep_idx, the collinearity flags and `_cons`;
- α̂ when FE are saved;
- the scaling (stdevs, means);
- a frozen copy of every option the replicate fit reads (tolerance, maxiter,
  standardize, collinear_tol, weight_type, singleton rule, separation settings,
  estimator, k-class parameters).

**Phases.** `ctools_command_begin` (`src/ctools_runtime.c`) gains cases for
`CTOOLS_CMD_CREGHDFE`, `CTOOLS_CMD_CPPMLHDFE` and `CTOOLS_CMD_CIVREGHDFE`:

| Subcommand | required phase | next phase |
|------------|----------------|------------|
| `full_regression retain` (`iv_regression retain`) | — | 1 (state kept) |
| `boot_chunk first(#) count(#) seedhi(#) seedlo(#)` | 1 | 1 |
| `boot_finish` | 1 | 0 (state freed) |
| `boot_release` | — | 0 (frees; called by the ado's error and Break path) |

**Errors.** `ctools_command_finish` already resets on any nonzero return code,
so a failing chunk frees the state. The command's cleanup hook (in
`ctools_commands.def`) frees the retained model as well as the existing state.

**Guard between calls.** The finish call writes to the observations that
`obs_map` points at, so the ado must not sort, add or drop observations between
the calls. Creating variables is fine. The retained model stores `SF_nobs()` and
the sort order (`c(sortedby)` is passed in); `boot_chunk` and `boot_finish` check
both and fail with r(459) otherwise.

**Protocol bump.** The ado-to-plugin contract changes, so `config/release.json`
moves to protocol 3.

### 6.4 SPI-free replicate fits

The SPI is not thread safe, and a replicate runs on a worker. The rule: all
inputs come from the retained model, all outputs go to `hdfe_boot_rep`, and
replicates never call `SF_*`, `ctools_scal_save` or `SF_display`, or write
globals.

- In each command, every option read becomes a field filled once on the Stata
  thread: `creghdfe_input` already does this for creghdfe; cppmlhdfe needs a
  `ppml_config`; civreghdfe needs an `iv_config`.
- Debug builds get an assertion in the SPI wrappers (`ctools_spi.c`) that the
  caller is the thread that entered `stata_call`.
- Verbose output describes the full-sample fit only. The replicates report
  through e(boot_status) and the iteration totals.

### 6.5 Scheduling and memory

The scheduler has two modes:

- **Outer (default).** `T_out` replicates run concurrently, each with a
  single-threaded `HDFE_State` (`num_threads = 1`) and its own workspace. A small
  K no longer limits parallelism.
- **Inner.** When the memory budget allows only one replicate at a time, the
  replicates run one after another, each with column parallelism, as the main fit
  does.

`T_out = min(threads, count, floor(mem_budget / workspace_bytes))`. If `T_out` is
1, the inner mode is used.

**Workspace size.** For creghdfe with K columns and G effects:

- N_r ≈ 0.632 · N rows in an observation bootstrap;
- columns: ≈ N_r · (K partialled + K raw when the FE target is on) doubles;
- ≈ N_r · (weights + 3 CG + residual) doubles;
- G · N_r level ints;
- per-effect means.

That is roughly N_r · (2K + 5) · 8 + G · N_r · 4 bytes. For example, N = 10M,
K = 5, G = 2 gives about 0.8 GB per concurrent replicate.

**Budget.** The default is half of `ctools_physical_memory_bytes()`
(`src/ctools_platform.h:15`), minus the retained model. A hidden option
`_bsmemory(#)` (MB) overrides it for tests and benchmarks. If even one replicate
does not fit, the command stops with r(909) before any replicate runs.

**Nested parallelism.** Combining outer and inner parallelism (threads split
between them) needs nested OpenMP regions, which are disabled today. It is not in
the first version (see §15).

### 6.6 Determinism

- **Draws.** Bitwise identical for any `threads()` and chunking.
  - Under `draws(fast)` they are a pure function of the master seed (taken from
    Stata's RNG) and r.
  - Under `draws(stata)` and for the wild weights, they are Stata's own sequence.
    The ado draws them in replicate order on the Stata thread, and chunking only
    decides when each chunk's draws are taken, never their order.
- **Replicate estimates.** Bitwise identical for a fixed `threads()` and chunking.
  The outer mode is single-threaded per replicate. The inner mode parallelizes
  only over columns, and each column is computed by one thread.
  - Any OpenMP reduction in the dense steps a replicate uses (cross products,
    sums) must be order-fixed. Audit `compute_xtx_xty*` and the `omp parallel
    for` loops in `cppmlhdfe_irls.c` and `civreghdfe_estimate.c` when moving them
    into replicate code.
  - Across `threads()` settings the estimates may differ in the last bits. Tests
    allow 1e-12 relative.
- **Accumulation.** Always in replicate order. Workers write to per-slot buffers,
  and the Stata thread folds a batch in index order. Floating-point sums in the
  accumulators therefore do not depend on the schedule.

### 6.7 Warm starts

**Projection: a Frisch–Waugh–Lovell start.** Let M_r be the weighted FE
annihilator for replicate r. For any column x and any full-sample FE fit Dα,
M_r(x − Dα) = M_r x, because M_r D = 0 for any weights. So each replicate
partials the full-sample *partialled* column x̃ = x − Dα̂_x instead of the raw x.
The result is exact up to rounding, and it is already almost orthogonal to D
under the replicate weights.

The CG stopping rule is relative to the norm of the input
(`hdfe_solver.c:375-382`). The smaller input therefore makes the absolute
tolerance tighter, not looser. The first version keeps the rule as it is. A
later change could pass the raw column's norm as the reference, to stop earlier
at equal absolute accuracy. G = 1 demeans directly, so the warm start costs
nothing there.

**IRLS (cppmlhdfe).** Start replicate r from the full-sample linear predictor η̂
(including offset) on its rows, with μ⁰ = exp(η̂) and the existing floor. The
first projection inside IRLS uses the full-sample final partialled columns as
the FWL start above, and later iterations reuse the previous partialled columns
as today.

**FE recovery: no warm start.** Both recovery routines fix the normalization by
starting from zero. A warm start from α̂ would change which solution is reached,
so recovery starts from zero. Phase 2 measures the cost and revisits this only
with an explicit re-normalization step.

### 6.8 Progress and Break

The ado runs the replicates in chunks. A chunk is `T_out × k` replicates, sized
so that a chunk takes about a second (measured on the first chunk). After each
chunk the ado:

- prints the chunk's dots with `_dots`, using its status vector: `.` ok, `x`
  failed, `e` collinear;
- posts the chunk's b rows for `saving()`;
- lets Break through. `_<cmd>_impl` already runs inside `capture noisily break`.

On Break or error, the ado calls `plugin ... "<cmd> boot_release"` and then its
existing rollback (`_ctools_undo`), which drops the new SE variables and
restores e(). `SF_poll` (`src/stplugin.h:198`) is not used as a Break signal
inside a chunk.

### 6.9 Ado helper: `build/_ctools_hdfe_boot.ado`

The helper is shared by the three commands. It has five entry points:

- **parse.** It turns `vce()` and `fevce()` text into a design: reps, seed,
  cluster term (via `_ctools_hdfe_vce`), strata, mse, level, dots, saving and FE
  options. It checks the cross-option rules in §4.1 and §4.6 and returns
  `r()` locals the command passes on.
- **ids.** It prepares the resampling-unit and stratum variables as plugin
  inputs. Strings and interactions are numbered with `_ctools_group`, and the
  nesting of clusters within strata is checked.
- **run.** It holds:
  - the RNG handling of §5.8: `seed()`, e(rngstate), and the two-draw master
    seed for `draws(fast)`;
  - the Mata draw routines: `bsample`-equivalent multiplicities for
    `draws(stata)`, and wild weights by type, written per chunk to temporary
    variables;
  - the chunk loop, the dots, the Break-safe release and the `saving()`
    postfile.
- **post.** It reposts e(V) and adds the e() results in §4.4. The results are
  computed in C at `boot_finish` and returned as K × 2 interval matrices and so
  on, so Stata never holds the R × K draws. That avoids matsize limits on Stata
  14.1. A chunk's b rows are at most `count × K`, and the ado keeps `count`
  within c(matsize).
- **display.** It prints the header lines and footnotes shared by the three
  commands.

The design must work in Stata 14.1, so it cannot use frames. `saving()` uses
`postfile`.

The same helper parses `wildbootstrap()` and runs it, through a `wild` entry
point with its own option rules (§4.6–4.7).

### 6.10 Wild-bootstrap engine: `src/hdfe/hdfe_wild.{c,h}`

The wild bootstrap is a second replicate type next to the reweighting one. It
perturbs residuals and never rebuilds factors or the sample. It shares the
retained model, the phase protocol, the seeds, the chunking, dots, Break
handling and the determinism rules (§6.3–6.8).

**Draws.** The ado draws the weights v with Stata's generator (§5.8) and passes
each chunk as C values per replicate. The plugin keeps the chunk's B × C matrix.
In enumeration mode the plugin generates the sign vectors itself, from r's
binary digits, and no weights are passed.

**Aggregate kernel.** The nested path and the score bootstrap use it.

- **Inputs:** per-group aggregates on the intersection of the `vce()` clusters
  and the bootstrap clusters, as a group → bootstrap-cluster map, S families
  (K-vectors linear in β0) and H_g.
- **Per replicate:** for each tested coefficient, the five coefficients of
  t*_b(β0)'s numerator and denominator (§5.9).
- **Parallelism:** replicates are independent and their outputs are fixed-size
  slots, so the scheduler runs them over all threads, B/threads each, with no
  memory concern.

**Projection kernel.** The projected path uses it.

- **Per batch:** for a batch of b replicates, build the N × b(k+1) columns
  (û ⊙ v and r_k ⊙ v) and partial them with `partial_out_columns` on the
  retained full-sample `HDFE_State`. This is the solver's native column
  parallelism, with a fixed state, so no per-replicate factor build is needed.
  Then aggregate the projected columns to clusters and hand them to the
  aggregate kernel.
- **Batch size:** b comes from the memory budget, at N · b(k+1) · 8 bytes.
- **The solver's per-thread buffer cap** (`ctools_hdfe_alloc_buffers` caps
  threads at `max_columns`) must be sized for b(k+1) columns, not K. The
  retained state therefore allocates its buffers for the larger count when a
  wild bootstrap is requested.
- **Tolerance:** projections use the fit's tolerance. With savefe the fit's 1e-12
  tolerance is not needed here, so the wild bootstrap uses the user's
  `tolerance()` (§7.5).

**IV basis.** WRE (§5.9) uses the aggregate kernel on a few more aggregate
families, and per (replicate, β0) a small dense solve. LIML and Fuller need a
(K₂+1)-dimensional eigenproblem through the existing
`civreghdfe_compute_liml_lambda` math, moved into an SPI-free helper.

Because the restricted first stage depends on β0 nonlinearly, t*_b(β0) no
longer reduces to five numbers per replicate. Each evaluation of p(β0)
therefore recomputes every replicate from the cluster aggregates and its
weights v. The plugin keeps the weight matrix the ado supplied:

- one byte per weight for Rademacher, Mammen and Webb (an index into the point
  set);
- doubles for normal and gamma.

If the matrix exceeds the memory budget, the command stops with r(909) before
any replicate runs. The cost is O(B · C · L²) per evaluation, where L is the
number of exogenous variables and instruments. On the projected path the per-replicate projected
basis columns are reduced to cluster aggregates once and kept, at
B × C × (basis size) doubles. When that exceeds the memory budget, inversion
falls back to re-projecting per evaluation, with a note on cost.

**Inversion.** One routine, `hdfe_wild_invert`, serves all three methods. It
takes a callback p(β0) and does the grid bracketing, bisection and open-side
detection.

**Phases.** Wild runs as `wild_chunk first(#) count(#) …` and `wild_finish`
calls, in phase 1, alongside `boot_chunk`. When the user also asks for
`fevce(bootstrap)`, both run in one plugin session: the wild bootstrap first,
then the pairs FE replicates, then `boot_finish`.

## 7. creghdfe

### 7.1 Refactor (no behavior change)

Split `do_full_regression` (`src/creghdfe/creghdfe_regress.c:216-1589`, one
1,370-line function) along its existing STEP comments into:

- `creghdfe_prepare`: load, remap, singletons, factors, DoF and compaction
  (STEPS 1–6);
- `creghdfe_fit`: standardize, partial, collinearity, OLS and residuals
  (STEPS 7–8 and the residual part of 9);
- `creghdfe_vce`: STEP 9;
- `creghdfe_publish`: the write-back and scalars of STEP 10, and the call to
  `store_fixed_effects`.

`creghdfe_fit` takes `(columns, weights, factors, config)` and returns
`(beta_keep, keep_idx, collinear flags, resid, rss/tss_within)` without any
SPI call. The full-sample path calls the same functions in the same order with
the same arithmetic (`docs/development/STYLE.md`: numerical refactoring keeps
arithmetic order). The SPI writes of collinearity flags now inside STEP 8 move
into `creghdfe_publish`.

Gate: e(b), e(V) and every e() scalar must be bitwise identical before and after,
on the whole creghdfe validation suite (§11.1).

### 7.2 Retained model and replicate fit

**Retained.** With `retain`, the call keeps:

- the compact raw columns (unstandardized);
- the partialled standardized columns (after STEP 7, the warm start);
- stdevs and means;
- the user weights;
- the factor levels;
- the unit and stratum ids, compacted;
- mask and `obs_map`;
- `creghdfe_input`;
- b̂, keep_idx and `_cons`;
- α̂ when savefe is on.

**Replicate fit** (`creghdfe_boot_fit`):

1. Rows with ω > 0; total multiplicity per level; singleton removal on
   multiplicities (`ctools_remove_singletons` with the multiplicities as fweights).
2. Factors on full-sample codes (`hdfe_build_factors_on_levels`) and
   `ctools_hdfe_alloc_buffers(…, K)` into the workspace's `HDFE_State`.
3. Compact the warm-start columns. Compute the norm2 values for the collinear-
   with-FE rule from the raw columns times the multiplicities. reghdfe's rule
   uses unweighted sums of squares, and after `expand` those are multiplicity
   weighted.
4. `partial_out_columns` (single-threaded in the outer mode), then the full
   sample's collinearity rules, then Cholesky and β. Destandardize and recompute
   `_cons` from the replicate's weighted means.
5. With the FE target, form d = y − `_cons` − Xβ − resid from the raw columns,
   run `hdfe_fe_recover`, set the levels that are not estimable to NaN (§5.7), and
   compute d per row for `dse()`.

**FE tolerance.** savefe already forces a 1e-12 projection tolerance
(`creghdfe_input.c:108-110`). Replicates use the same tolerance, which roughly
doubles to triples the CG iterations relative to a coefficient-only bootstrap.
The benchmarks report it separately (§11.4).

### 7.3 Ado changes (`build/creghdfe.ado`)

- **vce parsing (lines 173-211).** A `bootstrap` kind goes to
  `_ctools_hdfe_boot parse`. Add `FEVCE(string)` to the syntax (line 59).
  - With `vce(bootstrap, cluster(c))`, the plugin also loads c (its position is
    unchanged) so that ids exist for resampling. The analytic `vcetype` passed to
    the plugin becomes 1 (robust), because the full-sample V is replaced anyway.
  - In FE-only mode, the analytic vcetype stays the user's.
  - Strata go in as one extra variable after the cluster variable. The expected
    variable count check (`creghdfe_regress.c:315-327`) gains the new positions.
- **Output names.** The SE, `lb`/`ub`/`nrep` and `dse()` variables are added to
  `output_names` before `_ctools_newvars` (line 289), so the existing undo drops
  them on failure.
- **Order.** Plugin `retain` call, then the refit check (line 410), then the
  bootstrap chunks and `boot_finish`, then `ereturn post`. Posting after the
  bootstrap means a failed bootstrap leaves e() cleared, as any failure does
  today.
- **Refit.** The recursive refit call starts from a fresh `full_regression`.
  `ctools_command_begin` resets the earlier state; the protocol must allow this.
- **Display.** `_creghdfe_display` gains the bootstrap header and footnotes
  (§4.5).
- **noconstant.** `_cons` is bootstrapped and dropped with its row and column, as
  b and V are now (line 572).

### 7.4 What does not change

`groupvar()`, `residuals()`, e(df_a), e(dof_table), e(r2)… and the analytic F
in FE-only mode all describe the full-sample fit as they do today.

### 7.5 Wild bootstrap (WCR and WCU)

- **Retained inputs.** The `retain` call already has what the wild bootstrap
  needs (§7.2): X̃ and ỹ (partialled, standardized), the residuals, the weights,
  the cluster ids and the `HDFE_State`. It also keeps:
  - the `vce()` cluster ids and the `bootcluster()` ids, and their intersection;
  - the per-effect nesting flags. STEP 9 already computes them against the
    `vce()` cluster with `ctools_fe_nested_in_cluster`
    (`creghdfe_regress.c:1182-1262`). They are computed again against the
    bootstrap clusters when `bootcluster()` differs.
- **Setup** (Stata thread, before any chunk):
  1. r_k for each tested coefficient: the FWL residual of x̃_k on the other kept
     columns, from the retained cross products.
  2. The observed t from the analytic V the fit already computed.
  3. The aggregate families S (for û and each r_k) and H.
  4. The path choice: nested if every effect is nested in the bootstrap
     clusters, otherwise projected.
- **Units.** Standardization is undone in the aggregates so that t and t*
  match the posted e(b) and e(V). t is scale invariant, but the inversion bounds
  are in the units of β.
- **Tolerance.** On the projected path with `savefe`, the retained state's
  tolerance is the forced 1e-12 (`creghdfe_input.c:108-110`). The wild
  projections reset it to the user's `tolerance()`, because the wild statistics
  do not need fixed-effect accuracy.
- **Ado.**
  - `wildbootstrap()` is parsed after `vce()`.
  - The `wild_chunk` loop runs after the refit check and before `ereturn post`.
    The results are posted with the analytic ones.
  - `_creghdfe_display` prints the wild table (§4.7) after the analytic table.
  - The `noheader`, `notable` and `nofootnote` rules apply to it as to the
    analytic parts, and `nowildtable` hides it.

## 8. cppmlhdfe

### 8.1 Refactor

- **Options.** Move every `ppml_option(...)` and `SF_scal_use(...)` read inside
  IRLS and separation (`cppmlhdfe_irls.c:733, 921` and the option helper) into a
  `ppml_config` filled by `cppmlhdfe_input.c`.
- **Separation of side effects from IRLS.** Factor
  `do_ppml_regression` (lines 247 to about 1980) into prepare, screen, separation,
  IRLS, final fit, VCE and publish. Put the IRLS core (lines 936-1288) and the
  final fit (1345-1502) behind a function that takes columns, weights, factors,
  the config and a starting η, and returns β, η, the converged flag and the
  iteration count. The verbose `SF_display` calls (lines 513, 771, 884, 1091,
  1186, 1275, 1333) move to the caller, which receives an iteration log.
- **FE recovery.** Split `ppml_save_effects` (lines 73-156) into an SPI-free
  CGLS core and its write-back. Make the `d` computation (lines 1832-1840) a
  function.

Gate: the cppmlhdfe suite is bitwise identical before and after, including the
refit loop and the `guess()` paths.

### 8.2 Replicate fit

1. **Screening.** A multiplicity-aware `ppml_select_fe_sample_ex`: levels with
   total multiplicity 1 are singletons, and levels whose drawn rows all have
   y = 0 are separated. Repeat until stable, as the full sample does.
2. **Separation checks.** Simplex and ReLU run if the full fit ran them, on the
   replicate's rows and weights. This is the most expensive part of a PPML
   replicate. Phase 3 measures it and reports how often a replicate separates
   when the full sample does not. Skipping the checks when they found nothing in
   the full sample is **not** the default. A missed separation makes IRLS
   diverge, and the mu check (lines 1132-1139) and `irls_maxiter` must catch that
   as a failed replicate.
3. **IRLS.** Warm-started from η̂ (§6.7), with weights ω times the IRLS working
   weights, and the same convergence rules. Non-convergence marks the replicate
   failed.
4. **Fixed effects.**
   - The same normalization as the point estimate: minimum norm through CGLS
     from zero. The estimability mask here is over the joint connectivity of all
     effects, because the minimum-norm solution is global.
   - `d` per row for `d()` and `dse()`, with mean_d folded into `_cons` as in
     the full fit.
5. **Standardization.** Fixed to the full sample's `input_scales`
   (`scale_fixed`), so replicates use the same units.

### 8.3 Ado changes (`build/cppmlhdfe.ado`)

- **vce.** `vce(bootstrap)` at lines 367-426. A bootstrap `cluster()` with more
  than one term is r(198). A multiway analytic `vce(cluster a b)` stays allowed
  when only the fixed effects are bootstrapped.
- **Order.** The bootstrap runs once, after the refit `while` loop (lines
  634-667) has finished.
- **Outputs.** SE, `dse()` and `d()`-SE names go into the output validation and
  `_cppmlhdfe_publish`, inside the existing `nobreak`.
- **vce(unadjusted).** It maps to robust today (line 421). With `vce(bootstrap)`
  the analytic V is replaced, so the base vcetype does not matter.

### 8.4 Wild bootstrap (score bootstrap)

- **Restricted fits.** One per tested coefficient: H0: β_k = 0 means x_k·0 as
  the offset and x_k dropped.
  - They run through the SPI-free IRLS core from §8.1, warm-started from the
    unrestricted η̂.
  - They use e(sample), the full fit's screened sample, so the test compares
    the two models on one sample. Dropping a regressor cannot separate a row
    that was not separated before, so the restricted estimate exists on that
    sample.
  - A restricted fit that does not converge makes that coefficient's wild
    result missing, with a note. The analytic results stand.
- **Scores.** At each restricted μ̃: partial the regressors under μ̃·W with the
  retained factors, form r̃_k and the efficient cluster scores g_c, then run the
  aggregate kernel (§6.10). Per replicate this is O(C).
- **Intervals (`ci`).** `hdfe_wild_invert` calls p(β0). Each call is a
  restricted fit at β0 (offset x_k·β0), warm-started from the nearest evaluated
  β0. Expect 10–20 restricted fits per coefficient. The help file states the
  cost.
- **Multiway clustering.** Multiway `vce(cluster)` with `wildbootstrap()` is
  r(198) (§4.6).

## 9. civreghdfe

### 9.1 Phase 4a: fixed-effect saving (no bootstrap yet)

This is needed before any FE bootstrap, and is a feature in its own right.

- **Syntax.** `absorb(absvars, savefe)`. Replace the "no suboptions" error
  (`build/civreghdfe.ado:299-302`) with creghdfe's `_creghdfe_absorb_opts`
  pattern. It creates `__hdfe#__` and checks the names first.
- **Recovery.** d = y − `_cons` − Xβ̂ − ê, where X holds the actual endogenous
  and exogenous regressors (not the first-stage fitted values) and ê is the
  structural residual on the partialled data. This is the same formula as
  creghdfe's, so it calls `hdfe_fe_recover` with the same 1e-12 tolerance rule.
- **`partial()`.** The partialled-out exogenous regressors have no coefficients
  in e(b), and d needs them. The first version refuses `savefe` with
  `partial()` (r(198)). A follow-up recovers them by FWL (regress y − X₂β̂₂ on
  the partialled block).
- **Oracle.** `reghdfe (y − Xβ̂), absorb(..., savefe)` with no regressors
  reproduces the effects exactly. Also ivreghdfe, if it supports savefe in the
  installed version.

### 9.2 Phase 4b: refactor

- **Coefficient-only core.** Add `ivest_fit_coefs()` beside
  `ivest_compute_2sls`. It returns β for 2SLS, LIML, Fuller, k-class, GMM2S and
  CUE, with no diagnostics and no SPI calls. `__civreghdfe_coviv` (`estimate.c`)
  and similar option reads go into an `iv_config` read once.
- **Main path.** The full-sample path keeps calling `ivest_compute_2sls`. It is
  then refactored to call the core and then the diagnostics, with the bitwise
  gate on the civreghdfe suite.
- **Copy-aware S.** An S builder variant for replicates that takes the
  multiplicities and the analytic weights separately (§5.2):
  - cluster: Σ_c m_c s_c s_c′;
  - robust: Σ_i m_i (w_i e_i)² z_i z_i′ for aweights and pweights, and
    Σ_i m_i fw_i e_i² z_i z_i′ for fweights.

### 9.3 Phase 4c: bootstrap

- **Base VCE.** GMM2S and CUE need an S, and the diagnostics need a VCE, even
  when e(V) comes from the bootstrap. Default: cluster on the bootstrap's
  `cluster()` if given, otherwise robust. The diagnostics then use that base VCE.
  The footnote says so: "Diagnostics use the cluster-robust covariance;
  coefficient standard errors are bootstrap."
  - In FE-only mode the base is the user's `vce()`.
  - Whether to expose the base as an option is in §14.
- **Replicate fit.** Partial y, X_endog, X_exog and Z with the FWL warm start,
  then the FWL `partial()`, the collinearity rules (including endogenous in Z),
  `ivest_fit_coefs()`, and FE recovery when requested.
- **Inference is z-based.** `small` has no effect on a bootstrap V. e(df_r) is
  not posted, and the coefficient table and Wald statistic use z and chi2. The
  diagnostics keep their own small-sample conventions unchanged.
- **Auxiliary regressions.** `first`, `rf`, `savefirst` and `saverf` run the
  creghdfe auxiliary regressions with the base VCE, not the bootstrap.
- **Order.** The bootstrap runs once after any `civrerun` (lines 1024-1091).
- **Restrictions.** See §4.6. CUE is allowed with a note on cost, since it is an
  iterative optimization in every replicate.

### 9.4 Phase 4d: wild bootstrap (WRE)

- **Estimators.** 2SLS, LIML, Fuller and k-class, through an SPI-free
  restricted-estimation helper beside `ivest_fit_coefs()` (§9.2).
- **Inputs.** It uses the retained partialled ỹ, Ỹ₂, X̃₁ and Z̃ after
  `partial()`, the cluster ids and the nesting flags.
- **Studentization.** The analytic VCE of the user's `vce()` (robust or one-way
  cluster) with civreghdfe's small-sample conventions. They cancel in p and in
  the interval, as in §5.9.
- **Hypotheses.** Each coefficient in `coefs()`, endogenous or exogenous, within
  the set boottest's WRE supports (Phase 0, §5.9). The rest are refused with
  r(198).
- **Results.** The diagnostics, e(V) and the first-stage results are those of the
  analytic fit. The wild table follows the main table.

## 10. Phases and milestones

Sizes are relative: S is days, M is about a week, L is several weeks. Each phase
ends green on `python3 validation/run.py` for the touched commands, and on
`--checks options|docs|seeds`.

| Phase | Content | Exit criteria | Size |
|-------|---------|---------------|------|
| 0. Reference behavior (**done** on StataNow 19 with boottest 4.5.5, Appendix A; Stata 14.1 runs outstanding) | Record, on Stata 14.1 and a current Stata, what `regress, vce(bootstrap)` and the `bootstrap:` prefix do: default reps; e() list; casewise or per-statistic V; omitted coefficients; fweights and pweights; display lines. Read `bsample.ado` and the prefix's replicate loop for the `draws(stata)` sequence, frame and unit order (§5.8). Read boottest's source, and Stata's `wildbootstrap` output where available, for the weight-to-generator mapping, enumeration rule, p-value ties and denominator, the WRE hypotheses it supports, the score-bootstrap statistic, and how it treats effects not nested in the clusters. Build the oracle loops (§11.1, §11.1b). Decide the §14 questions. | Notes added to this document; oracles run. | M |
| 1. Engine (**done**, see below) | `hdfe_boot.{c,h}`, `ctools_rng_binomial`, the accumulators and summary; `_ctools_hdfe_boot.ado` (parse, seed, Stata-RNG uniforms). | Native tests (§11.3) pass; a toy model through the engine gives the closed-form V; `draws(stata)` counts equal `bsample, weight()` for the same seed. | M |
| 1w. Wild engine (**done**, see below) | `hdfe_wild.{c,h}`: aggregate kernel, projection kernel, enumeration, `hdfe_wild_invert`; the Mata weight draws. | Native tests: hand-computed t* on toy data; enumeration; inversion on a closed-form case. | M |
| 2a. creghdfe refactor (**done**, see below) | §7.1 split; `hdfe_fe_recover` extraction. | Bitwise e() gate on the creghdfe suite and benchmarks unchanged. | M |
| 2b. creghdfe coefficients (**done**, see below) | §7.2–7.3 without the FE target: `vce(bootstrap)` with cluster, strata, weights, mse, saving, dots and both draw modes. Also the pieces moved from Phase 1 because only a command exercises them: `hdfe_build_factors_on_levels`; the ado helper's ids, run (chunk loop, dots, Break release, `saving()`), post and display; phase-protocol support in `ctools_runtime.c`; and the protocol bump in `config/release.json`. `hdfe_components` and the normalization warning come with 2c. | Exact oracle; prefix parity with `draws(stata)`; statistical oracle; determinism and RNG tests; benchmark entry. | M |
| 2c. creghdfe fixed effects (**done**; `ci` in Phase 5, see below) | `fevce(bootstrap)`, estimability masks, SE/ci/nreps/dse outputs, e(fe_*), the normalization warning (§4.8). | FE exact oracle; G = 1, 2, 3 cases; warning tests; memory guard for `ci`. | M |
| 2d. creghdfe wild (**done** except the benchmark, see below) | §7.5: WCR and WCU, nested and projected paths, inversion, display. | Manual-loop and boottest oracles (§11.1b); path agreement; benchmark entry. | M |
| 3. cppmlhdfe (**done** except the benchmark; slopes in Phase 5, see below) | §8.1–8.3: refactor gate, then coefficients and FE (`savefe`, `newvar=`, `d()`), with the warning. | Exact oracle against ppmlhdfe with fweight multiplicities; prefix parity; separation-in-replicate cases; warm and cold starts agree within irls_tol. | L |
| 3w. cppmlhdfe wild (**done**, see below) | §8.4: score bootstrap, restricted fits, opt-in intervals. | boottest oracle after ppmlhdfe; manual-loop oracle. | M |
| 4a. civreghdfe savefe (**done**, see below) | §9.1. | savefe oracle; docs. | S–M |
| 4b. civreghdfe refactor (**done**, see below) | §9.2. | Bitwise gate on the civreghdfe suite. | L |
| 4c. civreghdfe bootstrap (**done**, see below) | §9.3, with the warning. | Exact oracle against ivreghdfe for 2SLS and LIML; GMM2S cluster oracle on idcluster-expanded data; diagnostics unchanged. | M |
| 4d. civreghdfe wild (**done**, see below) | §9.4: WRE for 2SLS, LIML, Fuller and k-class. | Manual-loop and boottest oracles (boottest after ivreg2 with the effects as dummies; see below); diagnostics unchanged. | L |
| 5. Docs and release (**done** except the protocol bump and the §14 confirmations, see below) | §12, and the options the first phases refused: `saving()`, `_bsdraws()`, `ci`, fweighted `draws(stata)`, slopes, civreghdfe savefe with `partial()`. | Docs check, option coverage, help files. | S |

**Phase 1 delivered.**

- **`src/hdfe/hdfe_boot.{c,h}`.**
  - The resampling frame (`hdfe_boot_design`).
  - `draws(fast)`: one keyed generator per replicate, multinomial counts
    within strata, and frequency weights resampled as units by conditional
    binomials.
  - `hdfe_boot_map_bsample`: bsample's draw-to-unit mapping in C, for
    observations and clusters, within strata. It reproduces Stata's
    `int(u·n + 1)` rounding.
  - The replicate scheduler (`hdfe_boot_run`: memory budget, a workspace per
    replicate in flight, OpenMP over replicates, folding in replicate order).
  - The accumulators: replicate matrix, statuses, and Welford sums for fixed
    effects and d.
  - The summary, as `_bs_sum` computes it: casewise V, `mse`, se and normal
    intervals about the replicate mean, `_pctile` percentiles, and
    bias-corrected intervals from #{b^r ≤ b̂}.
  - The fixed-effect SEs.
- **`src/ctools_rng.h`.** `ctools_rng_double` and `ctools_rng_binomial`
  (inversion below a mean of 10, otherwise BTRD).
- **`validation/native/boot.suite.c`** (`test_bootstrap_engine_native`).
  - `bsample` fixtures generated by Stata for four designs.
  - Draw properties.
  - Binomial moments.
  - A toy model through the scheduler: replicates and summary recomputed
    independently; invariance to chunking and threads; supplied draws;
    casewise drops; 198, 909 and 2000 contracts.
  - Percentile and normal-quantile checks.
  - It passes at -O0 and -O2. The engine compiles warning-free for Linux,
    macOS and Windows; this machine built with zig and without OpenMP.
- **`build/_ctools_hdfe_boot.ado`** (`parse`, `seed`, `uniforms`), listed in
  `build/ctools.pkg`.
  - The `bootstrap_helpers` case of `validation/stata/bootstrap_reference.do`
    checks the option parsing, the seed handling, and that the helper's
    uniforms reproduce `bsample`'s draws and its final RNG state (13 checks).
- **Not yet built.**
  - The OpenMP path of the scheduler is untested here. The native test runs
    it with OpenMP on Linux and macOS.
  - The wild-bootstrap weight draws (boottest's mapping) belong to Phase 1w.

**Phase 2a delivered.**

- **What was extracted.** The refactor took out what replicates need, rather
  than the full four-way split of §7.1.
  - `hdfe_fe_recover` (`src/hdfe/hdfe_fe_recover.{c,h}`): creghdfe's
    fixed-effect recovery without the SPI. `store_fixed_effects` is now a thin
    wrapper around it.
  - `creghdfe_collinear` and `creghdfe_ols` (`src/creghdfe/creghdfe_core.{c,h}`):
    STEP 8's collinearity rules and OLS, without the SPI. The scalar writes stay
    in `do_full_regression`.
- **Gate.** A zig-built plugin reproduced 13 creghdfe specifications bit for
  bit, at 17 significant digits: every e() scalar, e(b), e(V), e(dof_table),
  and the residual, savefe and groupvar variables. The specifications covered
  weights of each type, factor variables, `nostandardize quad`, `noabsorb`,
  `keepsingletons`, clusters, and a three-effect panel. The gate passed again
  with the Phase 2b plugin.

**Phase 2b delivered.**

- **Plugin**, in `src/creghdfe/creghdfe_boot.{c,h}`:
  - `full_regression retain` keeps the compact raw and partialled columns, the
    weights, the factor levels, the collinearity pattern and b̂, and builds the
    resampling frame in `bsample`'s unit order.
  - `boot_chunk` runs replicates with `hdfe_boot_run`. Each replicate applies
    the singleton rule to its multiplicities, partials the full-sample
    partialled columns under the replicate weights, applies creghdfe's
    collinearity rules (counting duplicated rows; a changed pattern makes the
    replicate missing), and fits by OLS.
  - `boot_finish` posts the `_bs_sum` summary. `boot_release` frees the model.
  - `ctools_runtime.c` holds creghdfe in phase 1 between these calls.
  - The fit and its cleanup reach the module only through two hooks
    (`creghdfe_hdfe.h`), which the dispatcher installs, so the native harnesses
    that compile the fit alone still link.
- **Ado.**
  - `build/creghdfe.ado` parses `vce(bootstrap, …)`, marks out and numbers the
    cluster and strata variables, runs the replicates in chunks of 50
    (`_creghdfe_bootstrap`: dots, both draw modes), and posts what
    `vce(bootstrap)` posts (`_creghdfe_bootstrap_post`), including
    e(V_modelbased).
  - The header shows replications and a Wald chi2.
  - On failure the ado releases the plugin's model.
  - `estat bootstrap`, replay and `mse` work.
- **Validation.**
  - The `creghdfe_bootstrap` case of `validation/stata/bootstrap_reference.do`
    has 15 checks. With `draws(stata)`, creghdfe reproduces
    `bootstrap _b: reghdfe` to 1e-12 or better, with the same N_reps,
    N_misreps and final RNG state, for five designs:
    - observations;
    - clusters;
    - strata with clusters that cross them and singletons outside e(sample);
    - a dummy that drops out of a replicate;
    - aweights.
  - The case also checks `draws(fast)` reproducibility, the posted e(), `mse`,
    replay, `estat bootstrap` and the option errors.
  - The old creghdfe test that expected `vce(bootstrap)` to be refused now
    checks that it runs.
  - The seeds check counts `vce(boot…)`, `fevce(boot…)` and `wildbootstrap` on
    the three estimators as draws unless they carry `seed()`; its self-test
    covers both directions.
  - The cases run against a test plugin through `CTOOLS_TEST_PLUGIN_DIR`.
- **Documentation:** `build/creghdfe.sthlp` and `docs/commands/creghdfe.md`.

**Departures from the plan.**

- **Clusters across strata.** Clusters that cross strata are units within each
  stratum, not an error (§5.1).
- **Protocol bump deferred.** The protocol in `config/release.json` stays 2
  until release. New ados with an old plugin fail only when `vce(bootstrap)` is
  asked for (an unknown subcommand), whereas a bump would stop every ctools
  command against the shipped plugins in this shared tree.
- **The fit's own VCE.** It is the unadjusted one, posted as e(V_modelbased)
  (robust with pweights).

**Not yet done in 2b.**

- `saving()` and the hidden `_bsdraws()`.
- `draws(stata)` with fweights, which needs `bsample` on the expanded data.
- A Break/fault-hook test of the release path.
- Benchmarks.
- A run of the OpenMP replicate path, which zig cannot build here.
- The Linux and macOS plugins (compile-checked only).

**Phase 2c delivered.**

- **Plugin** (`src/creghdfe/creghdfe_boot.c`):
  - `boot_chunk ... fe [dse]` makes each replicate also recover its fixed
    effects: the replicate's error y − _cons − Xb − e goes through
    `hdfe_fe_recover` on the replicate's factors, the replicate's levels map
    back to full-sample levels, and levels that are not estimable stay NaN.
    d, the summed effect, is computed for every full-sample row whose levels
    are all estimable (drawn or not).
  - Estimability uses the full sample's connected sets. With G ≥ 2 the sets
    are the components of the graph joining all G effects; consecutive pairs
    are also checked, because the recovery normalizes pair by pair. A set
    splits in a replicate when its drawn levels fall into more than one of the
    replicate's components; all its levels are then not estimable, and the
    replicate counts towards e(fe_split_reps).
  - `boot_finish ... fe [nreps] [dse]` writes the SE variables, the counts and
    d's SE on e(sample), e(fe_bs_table)'s matrix and the two scalars. `coef`
    adds the coefficient summary only when `vce(bootstrap)` was asked for, so
    the two targets share replicates or run alone.
  - The module defines `SD_SAFEMODE`: with the fast-mode SPI, `SF_vstore` into
    variables the ado had just created failed silently.
- **Ado.** `fevce(bootstrap[, sestub() nreps dse()])` needs `savefe`; with
  `vce(bootstrap)` the resampling options must agree (they may be given in
  either). Alone it takes the design from `fevce()` and leaves e(V) analytic.
  `_creghdfe_fe_post` posts the e(fe_*) results of §4.4 plus e(fe_rngstate),
  e(fe_draws), e(fe_cluster_bs) and e(fe_strata_bs). The §4.8 warning shows
  with every display, estimation or replay, before the footnote.
- **Validation.** The `creghdfe_fe_bootstrap` case of
  `validation/stata/bootstrap_reference.do` (8 checks) runs an exact oracle
  with `draws(stata)`: per replicate, the data are expanded by the
  multiplicities and `reghdfe ..., absorb(..., savefe)` is fitted; its level
  effects and d are merged back and their standard deviations compared.
  - One effect (auto, rep78) with a cluster bootstrap: per-level SEs, counts
    and d's SE match.
  - Two effects with aweights (a 40 × 8 panel, no splits): per-level SEs, d's
    SE and the counts match; e(fe_norm_warning) is 1.
  - `vce(bootstrap) fevce(bootstrap)` gives the same coefficients as
    `vce(bootstrap)` alone, and FEs as `fevce()` alone.
  - `sestub()`, `nreps`, and the option errors.
  - The oracle expands rows rather than using `[fw=m]`/`[aw=m·w]`, because
    reghdfe counts singletons by rows: an aweighted row drawn twice is not a
    singleton in the replicate's fweight semantics but is one to reghdfe.
  - Reruns with the 2c plugin: the 2a bitwise snapshot (13 specifications)
    is identical, and the reference (19), helper (13) and 2b (15) cases pass.
  - The new and changed C files compile without warnings (`-Wall -Wextra`) for
    x86_64 Linux (glibc 2.17) and arm64 macOS.
- **Documentation:** `fevce()` in `build/creghdfe.sthlp` and
  `docs/commands/creghdfe.md`, with the caveat and the stored results.

**Departures in 2c.**

- **Collinear replicates.** A replicate whose omitted-coefficient pattern
  differs leaves the fixed effects missing too (§5.5 planned to keep them).
  Recovering them needs the replicate's own coefficients, which the collinear
  path does not produce; keeping them would also make the two targets count
  different replicates.
- **Default names** are `__fese#__` and `__fenrep#__` (§4.3).
- **`ci`** is refused. It needs per-level replicate storage and the memory cap
  of §6.5.
- **Only the creghdfe path.** cppmlhdfe and civreghdfe follow in Phases 3
  and 4.

**Not yet done in 2c.** Oracle cases with G = 3 and with sets that split in
some replicates (§11.4 asks for both); a test of the warning's text; `ci`;
benchmarks; the OpenMP path.

**Phase 1w delivered.**

- **Engine** (`src/hdfe/hdfe_wild.{c,h}`, no SPI):
  - `hdfe_wild_setup` takes the partialled X̃, weights, residuals, A, β̂, the
    `vce()` and bootstrap cluster ids (NULL for one per row), the tested
    columns, WCR or WCU, and the path. It builds the intersection groups
    itself (a counting sort by `vce()` cluster), r_k = X̃a_k / a_kk without a
    regression, the aggregate families S_u and S_r per group, and
    d_g = A X̃_g′W_g ρ_g per `vce()` cluster.
  - `hdfe_wild_run` fills the five numbers (n0, n1, D00, D01, D11) of each
    replicate and tested coefficient, in chunks, with OpenMP over replicates.
    The aggregate path (nested effects, or `fedgp(partialled)`) costs
    O(H·K·(1 + k_tested) + G·K·k_tested) per replicate. The projected path
    (`fedgp(full)`) partials û ⊙ v and r_k ⊙ v through a `project()`
    callback in batches sized by a memory budget (909 if one replicate does
    not fit). WCU stores n1 = D01 = D11 = 0.
  - The observed statistic is the replicate with v = 1, computed by the same
    kernel (on the projected path, at the first run). The replicates with
    v = ±1 then tie with it exactly, as in boottest, so enumerated p-values
    do not depend on rounding.
  - `hdfe_wild_pvalue` follows boottest's `getp` (strict inequalities, over
    the feasible replicates; se* = 0 gives NaN and is left out).
  - `hdfe_wild_enumerable` and `hdfe_wild_enum_weights` follow boottest's
    rule and `count_binary` order.
  - `hdfe_wild_invert` is method-neutral (a p(β0) callback). It uses a
    65-point grid over ±2z·se, widens outward by doubling up to 1000·se, and
    bisects to 1e-6·se. It reports open sides, disjoint sets (lo and hi then
    bound the whole set) and empty sets.
- **Weights** (`build/_ctools_hdfe_boot.ado`): `wilddraw` draws boottest's
  matrices from Stata's generator, `wildweights` writes chunks of replicates
  to variables for the plugin, and `wildclear` drops them.
  - Mata fills matrices row by row (checked), so a chunk of columns is not a
    slice of the stream. `wilddraw` therefore draws the whole
    C × (reps + 1) matrix (two for Webb) as boottest does and drops the
    first column, at 8·C·(reps + 1) bytes. boottest's `matsizegb()`
    splitting is not mirrored.
  - The weights are kept at unit variance (boottest's WRE expressions, with
    m = 1), not the scaled values boottest uses outside WRE. The linear t is
    scale-invariant, while WRE (Phase 4d) needs unit variance, so one draw
    serves every method.
  - Enumeration draws nothing and leaves the RNG state unchanged.
- **Validation.**
  - `wild_engine` in `validation/native/boot.suite.c`
    (`test_wild_engine_native`, with an OpenMP build where available). It
    compares t*_b(β0) for every replicate, both tested coefficients and four
    β0 values with a direct dense computation: the explicit restricted
    regression, the perturbed data partialled or not, OLS, and the
    cluster-robust V. Agreement is to 1e-9 for WCR and WCU, on both paths,
    in four designs: an effect nested in the clusters, an effect crossing
    them, coarser bootstrap clusters, and robust. The case also checks:
    - all four p-value types;
    - that the two paths agree on a nested design (1e-10);
    - determinism across chunks, threads and projection batches (bitwise);
    - exact ties at v = ±1;
    - enumeration in `count_binary` order;
    - the error codes;
    - inversion on a closed-form p-function and on open, disjoint and empty
      sets, and the engine's own bounds.
    It passes at -O0 and -O2.
  - The gate: the `wild_draws` case of `validation/stata/bootstrap_reference.do`
    (16 checks). With the same seed, a replicate-by-replicate WCR computed in
    Mata (`_bsref_wild_t`) from the `wilddraw` weights reproduces boottest's
    per-replicate statistics (`svmat`, to 1e-8; differences seen were near 1e-13) and its p-value. This holds
    for every weight type (34 clusters, 19 replicates), and the final RNG
    state is boottest's. It also holds for enumeration (8 clusters, 256
    vectors in boottest's order, no draw) and for the observation-level
    bootstrap after `vce(robust)` (rows in data order). The case also covers
    `noenumerate`, the moments and support of each weight type,
    `wildweights` and `wildclear`, and the errors.
  - `hdfe_wild.c` compiles without warnings (`-Wall -Wextra -Wpedantic`,
    C11) for x86_64 Linux, arm64 macOS and Windows. The seeds, docs and
    options checks, the architecture and release checks, the helper case
    and the pairs engine suite still pass.

**Departures and notes from 1w.**

- **The wild weights stay in Mata, not as bytes in the plugin.** The linear
  path keeps only the five numbers per replicate, so the plugin need not keep
  the weights. The one-byte storage of §6.10 is for WRE (Phase 4d).
- **S_r is stored per group** (k_tested × H × K doubles). That is large when
  H ≈ N (observation-level bootstrap with many tested coefficients); 2d
  should compute those sums from the rows in that case.
- **Inversion grid.** It covers ±2z·se, twice the analytic interval, before
  widening.
- **Not in 1w:** the plugin and ado plumbing (`wild_chunk`, `wild_finish`;
  2d), the IV basis and LIML eigenproblem (4d), the score bootstrap (3w),
  and a run of the OpenMP path, which zig cannot build here.

**Phase 2d delivered.**

- **Plugin** (`src/creghdfe/creghdfe_boot.c`, on the fit kept by
  `full_regression retain`):
  - `wild_setup` rebuilds the partialled kept regressors and y in the units
    of the data (the retained columns times their scales), the residuals and
    A = (X̃′WX̃)⁻¹ (Cholesky). It numbers the `vce()` and bootstrap clusters
    in ascending order on the estimation sample, maps the tested regressor
    positions (matrix `__creghdfe_wb_test`) to kept columns, and checks
    whether every absorbed effect is nested in the bootstrap clusters. The
    path is nested (aggregate), projected (`fedgp(full)` with an effect not
    nested) or partialled (`fedgp(partialled)`).
  - `wild_chunk` builds the engine at the first chunk. On the projected
    path it also builds the full sample's factors with the fit's weights,
    once, and partials batches of replicates' columns with
    `partial_out_columns` (memory budget a quarter of RAM, at most 2 GB).
    The weights come in as variables, 50 replicates per chunk, or are
    enumerated in the plugin.
  - `wild_finish` writes p-values at β0 = 0, the test-inversion sets and
    their status, and e(wb_misreps) (B less the fewest feasible replicates
    over the tested coefficients). With `keep` it holds the fit for
    `fevce(bootstrap)`, whose pairs replicates follow (`ctools_runtime.c`:
    `wild_setup`, `wild_chunk`, and `wild_finish keep` stay in phase 1).
- **Ado** (`build/creghdfe.ado`):
  - `wildbootstrap` or `wildbootstrap(wildopts)` is taken from the free
    options (`_creghdfe_wild_opt`) and parsed by `_ctools_hdfe_boot
    wildparse`, shared with the other estimators. `noci` is stripped before
    `syntax`, which cannot declare both `ci` and `noci`.
  - The refusals of §4.6 apply, and `bootcluster()` is marked out like the
    `vce()` cluster.
  - `_creghdfe_wild` runs `wild_setup`, `wilddraw`, the chunk loop with
    dots, and `wild_finish`. `_creghdfe_wild_post` posts the e(wb_*)
    results of §4.7, plus e(wb_cistatus) (open lower, open upper, disjoint,
    empty, failed) and e(wb_rngstate). e(wb_t) is the analytic t from e(b) and e(V).
  - `_creghdfe_wild_table` prints the table of §4.7 after the coefficient
    table, also on replay, unless `nowildtable`. Its notes cover
    enumeration, Rademacher weights with 11 or fewer clusters,
    `fedgp(partialled)`, open or disjoint sets, and infeasible replicates.
- **Validation** (the `creghdfe_wild` case of
  `validation/stata/bootstrap_reference.do`, 20 checks; boottest 4.5.5 with
  the same seed):
  - p-values equal boottest's exactly after `reghdfe`, for both
    regressors, in six designs. These are a nested effect (Rademacher), a
    crossed effect on the projected path (Webb), aweights (Mammen), WCU
    (`nonull`), the observation-level bootstrap under `vce(robust)`
    (normal), and equal-tailed p-values (gamma).
  - The interval bounds are checked through boottest: its p-value 1e-5·se
    inside each bound exceeds α, and 1e-5·se outside it does not.
  - Two absorbed effects match boottest after `areg` with the second effect
    as dummies, on the projected path.
  - `fedgp(partialled)` matches boottest after `regress` on the demeaned
    data with `nocons`.
  - Enumeration (8 clusters, 256 vectors) matches boottest's p-value and
    leaves the RNG state unchanged.
  - `bootcluster()`: 15 coarser clusters match a replicate-by-replicate
    oracle in Mata (`_bsref_wild_t2`).
  - The hidden `_projected` forces the projected path on a nested design:
    the p-values are identical, and the intervals agree to 1e-8.
  - e(b) and e(V) are unchanged by the wild bootstrap. `coefs()`, `noci`,
    `nowildtable` and replay work.
  - With `fevce(bootstrap)`, each target equals its result when run alone.
  - The seven refusals give r(198).
  - Reruns with the 2d plugin: the 2a bitwise snapshot is identical, and
    the 2b (15), 2c (8), wild draws (16), helpers (13) and reference (19)
    cases pass. The changed C files compile without warnings for Linux,
    macOS and Windows.
- **Documentation:** `wildbootstrap()` and `nowildtable` in
  `build/creghdfe.sthlp` and `docs/commands/creghdfe.md`.
  `validation/option_evidence.json` covers `nowildtable`.

**Departures and findings in 2d.**

- **Interval bounds are exact, not boottest's.** boottest's search
  (`search()` in `boottest.mata`) stops once p is within one replicate step
  (1/B, or 2/B equal-tailed) of α. Its bounds can therefore differ from the
  exact crossing by up to a few percent of a standard error at B = 199. The
  commands bisect to the crossing of p(β0), which boottest's own p-values
  confirm.
- **boottest cannot check `bootcluster()`.** boottest 4.5.5 stops with
  r(3301) in `crosstabCapstarMinus` whenever `bootcluster()` is coarser
  than the clusters, after `reghdfe` and after `regress`. The Mata oracle
  takes its place.
- **`wildbootstrap` is not in the option inventory.** It is read from the
  free options, because `syntax` cannot accept both `wildbootstrap` and
  `wildbootstrap(...)`, so `--checks options` does not list it. Its
  coverage is the `creghdfe_wild` case.
- **`docs/OPTION_COVERAGE.md` was not regenerated.** `--checks options
  --write` refuses while another session's `crangestat.casewise` gap is
  open, so the `nowildtable` row appears at the next regeneration.
- **Not in 2d:** the benchmark entry, a run of the OpenMP paths, and
  computing the S_r sums from the rows for observation-level bootstraps
  with many tested coefficients (Phase 1w note).

**Phase 3 delivered: coefficients (3a).**

- **Refactor** (`src/cppmlhdfe/cppmlhdfe_irls.c`, `cppmlhdfe_fit.h`). The
  IRLS loop (STEP 9), the final projection and collinearity check, and the
  final solve (STEP 10) moved into `ppml_irls_loop`, `ppml_final_collinear`
  and `ppml_final_solve`. They share one `ppml_irls_io` struct.
  - The bodies moved verbatim behind local aliases.
  - Messages print only with `report`.
  - The SPI reads (`ppml_option("log")`, `__cppmlhdfe_inner_tol`) and the
    scalar saves (collinearity flags, e(rss)) stayed with the caller.
  - `ppml_refit_fixed_effects` gained a `report` flag.
  - **Gate:** a zig-built plugin reproduces 15 cppmlhdfe specifications bit
    for bit before and after: every e() scalar, e(b), e(V), and the savefe
    and d() variables. They cover weights of each type, offset and
    exposure, a collinear regressor, a regressor that separates (the refit
    loop), multiway clusters, `guess(ols|simple)`, a slope effect and a
    tight tolerance.
  - This is narrower than §8.1. Preparation, screening, the VCE and
    publishing stay inside `do_ppml_regression`, and the options are still
    read from Stata scalars there. The replicates take a copy of the struct.
- **Replicates** (`src/cppmlhdfe/cppmlhdfe_boot.{c,h}`):
  - `full_regression retain` makes the fit hand its standardized data,
    factor levels, final η, scales, options and b̂ to a hook (so the native
    harnesses still link), after the fixed effects are saved and before
    the un-standardization.
  - The retained separation checks are the user's (`__cppmlhdfe_bs_check_*`),
    not the ones the refit loop turned off.
  - Each replicate:
    1. screens the drawn rows as `ppml_select_fe_sample_ex` does, counting
       each row as often as it was drawn (rows, not fweights, as cppmlhdfe
       counts);
    2. builds its factors and runs the simplex and ReLU checks with the
       multiplicities as frequency weights;
    3. warm-starts the IRLS from η̂ in the full fit's units, with the
       collinearity pre-pass;
    4. runs the mu check.
  - Rows a check finds leave the replicate, which is screened and refit
    (at most 10 passes), as the full fit's refit loop does. Non-convergence
    (r(430) included) and a changed omission pattern make a replicate
    missing. The constant follows the full fit's formula.
  - `boot_chunk`, `boot_finish` and `boot_release` mirror creghdfe's;
    `ctools_runtime.c` holds cppmlhdfe in phase 1 between them.
- **Ado** (`build/cppmlhdfe.ado`):
  - `vce(bootstrap, …)` is parsed before the analytic `vce()`; the fit uses
    robust.
  - It refuses heterogeneous slopes, `tagsep()`/`zvarname()`, a model
    without regressors, `saving()` and `_bsdraws()`.
  - The replicates run once, after the refit loop, and a failure releases
    the retained model.
  - e(V) comes from the replicates through `_ctools_ebv`, e(V_modelbased)
    is the robust V, the Wald chi2 uses the bootstrap V, and the
    `vce(bootstrap)` e() list is posted (`_cppmlhdfe_bootstrap_post`).
  - The header shows Replications. `estat bootstrap` goes through
    `_cppmlhdfe_estat` to `estat_default`.
- **Validation** (the `cppmlhdfe_bootstrap` case of
  `validation/stata/bootstrap_reference.do`, 14 checks):
  - With `draws(stata)`, e(V), N_reps, N_misreps and the final RNG state
    equal `bootstrap _b: cppmlhdfe` (replicates saved in double) to 1e-8.
    That covers five designs: observations; clusters; strata with
    clusters; an offset with a rare dummy that drops out of some
    replicates; and a regressor that separates rows in some replicates.
    The cluster design agrees to about 1e-14 despite the warm start.
  - The prefix on `ppmlhdfe`, an independent estimator, agrees to 1e-6.
  - `draws(fast)` reproducibility across thread counts, the posted e(),
    `estat bootstrap`, replay, `mse`, and five refusals.
  - Reruns with the Phase 3 plugin: both bitwise snapshots (creghdfe and
    cppmlhdfe) are identical, and the creghdfe cases (2b, 2c, wild) and the
    helper case pass. The new and changed C files compile without warnings
    for Linux, macOS and Windows, and the native harnesses still link.
- **Documentation:** `vce(bootstrap)` in `build/cppmlhdfe.sthlp` and
  `docs/commands/cppmlhdfe.md`.

**Departures in Phase 3.**

- **Fweights.** cppmlhdfe resamples rows also with fweights, as the prefix
  does on unexpanded data. creghdfe resamples the observations they stand
  for (`draws(fast)`).
- **Separation frequency.** How often a replicate separates when the full
  sample does not is not measured yet. The separating-regressor design
  shows the path works and matches the prefix.

**Phase 3b delivered: fixed effects.**

- **Refactor.** `ppml_save_effects` is split into `ppml_fe_minnorm` (CGLS from
  zero, SPI-free, messages with `report`) and its write-back. The 15-spec
  cppmlhdfe snapshot is still bit for bit identical.
- **Replicates** (`cppmlhdfe_boot.c`, `replicate_effects`):
  - After a successful replicate, d = z − resid − Xb − mean(d) is formed
    as the full fit forms it, on the replicate's rows with its IRLS
    weights.
  - `ppml_fe_minnorm` gives its minimum-norm coefficients, which map to
    the full-sample levels.
  - Estimability is the plan's: drawn, and the level's full-sample
    connected set over all effects jointly did not split. The set is
    computed once at retain, and per replicate with union-find. Levels in
    a split set are NaN, and the replicate counts towards
    e(fe_split_reps).
  - d is the sum of the effects on every full-sample row whose levels are
    all estimable. A recovery that does not converge leaves the
    replicate's effects NaN (its coefficients stand).
  - `boot_chunk … fe dse` turns the targets on. `boot_finish` writes the SEs
    of the targeted effects (matrix `__cppmlhdfe_bs_fe_targets`), their
    counts, d's SE, e(fe_bs_table)'s matrix and the scalars. It writes the
    coefficient summary only with `coef`.
- **Ado.**
  - `fevce(bootstrap[, sestub() nreps dse()])` needs saved effects or `d()`.
    With `vce(bootstrap)` the resampling options must agree.
  - Names: `newvar_se`; `__fese#__` for `savefe`'s `__hdfe#__`, replaced on
    reruns like it; `d_se` for `d()` unless `dse()` names it; and
    `sestub()`, with `_nrep` counts.
  - The outputs are filled in temporary variables and published with the
    others through `_cppmlhdfe_publish`, in the existing transaction.
  - `_cppmlhdfe_fe_post` posts the e(fe_*) results.
    `_cppmlhdfe_fe_warning` shows with two or more effects, at estimation
    and on replay.
- **A fix in both commands.** e(rngstate), e(fe_rngstate) and
  e(wb_rngstate) were posted wrapped in literal compound quotes, because the
  post programs declared `RNGSTATE(string asis)`. They now hold
  c(rngstate) as set by `seed()`, which the new case checks for creghdfe and
  cppmlhdfe.
- **Validation** (the `cppmlhdfe_fe_bootstrap` case, 9 checks):
  - Two exact oracles with `draws(stata)`: in each replicate, the expanded
    data are refit with cppmlhdfe, `absorb(f=…)`, and the effects and their
    sum are merged back. Per-level SEs, counts and d's SE agree to 1e-5
    for one effect with a cluster bootstrap, and for two effects with
    `savefe` and `d()` (no splits; the warning is set).
  - Shared replicates with `vce(bootstrap)` give the same results as
    either alone.
  - Automatic outputs are replaced on a rerun, `sestub()` works, the
    posted random-number states are correct, and five refusals give
    r(198).
  - Reruns with the 3b plugin: both bitwise snapshots are identical, and
    all eight bootstrap cases pass (2b 15, 2c 8, wild 20, wild draws 16,
    helpers 13, reference 19, cppmlhdfe coefficients 14, fixed effects 9).
    The native harnesses link.
- **Documentation:** `fevce()` in `build/cppmlhdfe.sthlp` and
  `docs/commands/cppmlhdfe.md`. `validation/option_evidence.json` covers
  `cppmlhdfe.fevce`.

**Still to do in Phase 3.** Heterogeneous slopes (with `vce()` or
`fevce()`), measuring how often replicates separate, the benchmark, and an
OpenMP run.

**Phase 3w delivered: the score bootstrap.**

- **boottest's statistic, pinned.** After ML, boottest refits under the null
  with `constraints()`. It then refits without them from that point with
  `iterate(0)`, to take the unrestricted Hessian A = e(V_modelbased) and the
  scores at the restricted estimate. For H0: b_k = b0, a row's
  contribution is s_i'A e_k, which is proportional to (y_i − μ̃_i) r̃_i, with
  r̃ the residual of x_k on the other regressors and the effects under
  weights w·μ̃. With E_c its sums per bootstrap cluster, a replicate's
  statistic is N*/sqrt(D*): N* = Σ v_c E_c, and
  D* = Σ_g (Σ_{c in g} v_c E_c − share_g N*)², where share_g is cluster g's
  share of the observations (of the weights when weighted). The observed
  statistic is the replicate with v = 1. A Mata computation of this formula
  reproduces boottest's `svmat` statistics to 1.6e-12, unweighted.
- **A departure: weighted scores.** boottest sums the ML scores unweighted
  (`_panelsum` without weights). The weights enter only through A and the
  shares. With fweights its statistics then differ from its own on the
  expanded data (t*/t by up to 0.14 in the probe), which is inconsistent.
  cppmlhdfe weights the contributions, w (y − μ̃) r̃, and reproduces boottest
  on the expanded data. Unweighted, the two are identical.
- **Engine** (`src/hdfe/hdfe_wild.{c,h}`): `hdfe_wild_score_*` builds the
  vce() × bootstrap cluster groups and shares, stores the weights (or
  enumerates), and computes the observed and replicate statistics
  (OpenMP over replicates). Ties at v = ±1 are exact.
  `hdfe_wild_pvalue_t` gives boottest's p-value from such a vector.
  `hdfe_wild_invert_ex` takes the grid size and the bisection tolerance, for
  costly p-functions.
- **Plugin** (`src/cppmlhdfe/cppmlhdfe_boot.c`):
  - `restricted_fit` runs the SPI-free IRLS (`ppml_irls_loop`) on the
    retained sample, with x_k dropped and x_k·b0 added to the offset. It is
    warm-started from the full fit's μ̂, or from the previous trial value
    in the interval search.
  - It then partials x_k on the effects and the other kept regressors at
    the restricted means and forms the contributions.
  - `wild_setup` numbers the clusters and maps the tested positions.
    `wild_chunk` stores the weights. `wild_finish` gives the p-value at
    b0 = 0 per coefficient, and with `ci` inverts the test: 9 grid points
    over ±2z·se, then bisection to 1e-4·se, about 30–50 restricted fits
    per coefficient.
  - A restricted fit that fails leaves that coefficient missing
    (e(wb_fail)). `wild_finish keep` holds the fit for `fevce()`.
- **Ado** (`build/cppmlhdfe.ado`; table in the shared
  `build/_ctools_wild_table.ado`, in `ctools.pkg`):
  - `wildbootstrap[(wildopts)]` is taken from the free options.
  - It refuses `nonull`, `vce(bootstrap)`, multiway clusters, and robust
    with fweights without `bootcluster()`, besides the bootstrap
    refusals (slopes, `tagsep()`, no regressors).
  - The analytic standard errors from `__cppmlhdfe_V_full` and
    `__cppmlhdfe_scales` set the scale of the interval search.
  - The run happens before the pairs replicates of `fevce()`.
  - e(wb_method) = `score` and the other e(wb_*) results are posted, and
    the table shows after the coefficient table and on replay, unless
    `nowildtable`.
  - Intervals are opt-in (`ci`).
- **Validation** (the `cppmlhdfe_wild` case, 17 checks; boottest 4.5.5 after
  `poisson` with dummies, same seed):
  - p-values are equal for both regressors in five designs: clusters
    (Rademacher), two effects (Webb), an offset (normal), a cluster per
    observation (Mammen) and equal-tailed (gamma).
  - `vce(robust)` gives exactly what `vce(cluster obs)` gives. boottest
    4.5.5 fails with r(3200) after `poisson, vce(robust)`, so robust is
    checked that way.
  - fweights agree with boottest on the expanded data.
  - Enumeration (8 clusters, 256 vectors) matches, with the RNG state
    unchanged.
  - The `ci` bounds sit where boottest's p-value crosses the level (1e-3·se
    inside and outside). boottest gives no intervals after ML.
  - e(b) and e(V) are unchanged; `coefs()`, no intervals by default,
    e(wb_rngstate), replay and `nowildtable` work.
  - With `fevce()` each target gives the same result as alone, and six
    refusals give r(198).
  - Native: the `wild_engine` suite adds the score kernel against its
    definition (three designs), enumeration, exact ties, p-values and
    errors; it passes at -O0 and -O2.
- **Documentation:** `wildbootstrap()` in `build/cppmlhdfe.sthlp` and
  `docs/commands/cppmlhdfe.md`, and `nowildtable` evidence in
  `validation/option_evidence.json`.

**Not in 3w:** `nonull` (the unrestricted score bootstrap), multiway
clustering, the benchmark, and an OpenMP run. creghdfe keeps its own copy of
the table (`_creghdfe_wild_table`); it could move to `_ctools_wild_table`.

**Phase 4a delivered: civreghdfe saves fixed effects (IV-03).**

- **Syntax** (`build/civreghdfe.ado`): `absorb(…, savefe)` and
  `absorb(newvar=absvar …)`, as reghdfe; `_civreghdfe_absorb_opts` takes
  `savefe` and refuses other suboptions.
  - The names must be new and are checked before estimation.
  - `savefe` names the unnamed effects `__hdfe#__`.
  - It is refused with `partial()` (r(198)), and when a term is not one
    effect.
  - Because the recovery goes through creghdfe's `__hdfe#__`, an existing
    `__hdfe#__` of another name is r(110).
- **Recovery, in the ado rather than C.** Once e() is final, y − Xb is formed
  with `matrix score` on e(b) (the actual regressors). creghdfe then absorbs
  it with `savefe keepsingletons` on e(sample), with the fit's weights,
  under `_estimates hold`.
  - creghdfe's residual is then the IV structural residual on the
    partialled data. Its effects are those of §9.1's formula, with
    reghdfe's normalization and creghdfe's 1e-12 savefe tolerance.
  - ivreghdfe partials the constant out with the effects and posts none, so
    it keeps the constant in the first effect. civreghdfe does the same
    (creghdfe's _cons added to the first effect), so the effects sum to
    y − Xb − e.
  - The collinearity rerun (`civrerun`) passes the options on, so the
    second pass saves.
  - e(savefe) lists the variables.
  - The C-level recovery (`hdfe_fe_recover`) the plan describes is left for
    4c, where the replicates need it.
- **Validation** (the `civreghdfe_savefe` case of
  `validation/stata/bootstrap_reference.do`, 8 checks):
  - 2SLS with clusters matches reghdfe on y − Xb (its constant added to the
    first effect) and ivreghdfe's saved effects, to 1e-6 (2e-8 seen).
  - LIML with aweights matches ivreghdfe.
  - e(b) and e(V) are unchanged, only a named effect is saved (labeled),
    and saving works through the collinearity rerun.
  - The refusals are tested.
  - The `oct27` case of `validation/stata/regressions.do`, which asserted
    that civreghdfe refuses `fe=` and `savefe`, now asserts that it saves
    them (and still refuses other suboptions); it passes.
- **Documentation:** `build/civreghdfe.sthlp` and
  `docs/commands/civreghdfe.md`.
- **Not in 4a:** heterogeneous slopes in civreghdfe's absorb() (IV-02), and
  savefe with partial() (§9.1's FWL follow-up).

**Phase 4b delivered: the SPI-free coefficient core and the copy-aware S.**

- **`ivest_fit_core(civ_fit_io *)`** (`src/civreghdfe/civreghdfe_estimate.c`)
  is §9.2's `ivest_fit_coefs()`. It holds the coefficient phase of
  `ivest_compute_2sls` moved verbatim:
  - k (LIML's λ, Fuller, k-class) and the k-class solve;
  - `wmatrix()`;
  - the moment covariance S of the chosen VCE (`bw(auto)` included);
  - GMM2S and CUE.
  - Messages print only with `io->report`.
  - The scalars it used to save (`__civreghdfe_k`, `__civreghdfe_bw_auto`)
    come back in `io` and are saved by `ivest_compute_2sls`, which fills
    `civ_fit_io`, calls the core and then runs the VCE and diagnostics as
    before.
  - Its arrays are handed back also on failure, for the caller to free.
- **Copy-aware S** (`civ_omega.mult`, `mult_groups`; `civ_fit_io` passes
  them through). A replicate's row copies are separate observations:
  - the heteroskedastic and centering terms use w²/m for aweights and
    pweights (fweights and no weights are already linear in m);
  - with `mult_groups` (clusters are the resampled units) each copy is its
    own cluster, so a cluster's lag-0 outer product is divided by its count.
  - NULL/0 for the full fit, so nothing changes there.
- **Option reads.** The `iv_config` struct the plan mentions was not needed:
  the core takes everything through `civ_fit_io`, and the SPI reads stay in
  `ivest_compute_2sls` and `civreghdfe_impl.c`.
- **Bitwise gate.** A 16-spec snapshot covers:
  - 2SLS, LIML (aw), Fuller (fw), k-class, GMM2S (robust and pw cluster),
    CUE;
  - two-way cluster LIML, `bw()`, `dkraay()`, `bw(auto)` GMM2S;
  - orthog/endog/redundant, savefe with residuals, partial() and b0().
  - Every e() scalar and matrix to %21.17g is identical before and after.
  - The civreghdfe suites pass: 463, grouped syntax 17, cs 8.
- **Not in 4b:** a native test of the copy-aware S; Phase 4c validates it
  end to end against idcluster-expanded data.

**Phase 4c delivered: civreghdfe `vce(bootstrap)` and `fevce(bootstrap)`.**

- **Plugin** (`src/civreghdfe/civreghdfe_boot.{c,h}`, new):
  - `iv_regression retain` keeps the model after the full fit:
    - the raw compact columns and the fixed-effect-partialled ones (before
      `partial()`), weights and levels;
    - the ado's `partial()` flags, the fit's omitted regressors and
      instruments, and its count of endogenous regressors collinear with Z;
    - b in e(b)'s layout, the estimator settings, and the GMM moment
      covariance's kind.
  - `civreghdfe_main` moved to `civreghdfe_boot.c`. `civreghdfe_impl.c`
    calls the retain and release hooks, so the native harness, which
    includes `civreghdfe_impl.c`, links without the bootstrap. `civ_wss`,
    `civ_wss_dev`, `civ_invsym_drops` and `civ_fwl_partial_out` (messages
    only with `report`) are now extern.
  - `ctools_runtime.c` has the civreghdfe phase rules.
  - The bitwise snapshot of 4b is unchanged.
- **Replicate fit:**
  - singletons on the expanded counts;
  - the partialled columns projected again under the replicate weights,
    with the fit's power-of-two column scaling;
  - aweights and pweights normalized to the expanded count, in the
    projection *and* the factors' weighted counts (a mismatch there biased
    single-effect fits, which demean in one pass);
  - `partial()` by FWL;
  - the fit's collinearity rules (absorbed columns, invsym() of X'WX and
    Z'WZ, endogenous regressors in Z). A replicate whose pattern differs is
    COLLINEAR.
  - Then `ivest_fit_core` with the copy-aware S. That S is the fit's
    robust or one-way cluster S for GMM2S and CUE; iid otherwise, where b
    does not depend on S. With the bootstrap's clusters as S's clusters,
    each drawn copy is a cluster of its own.
- **fweights:** rows resampled one by one draw the fweighted observations
  (`unit_freq`), each row weighing the number drawn. With clusters, the
  weight is copies times fweight. `draws(stata)` with fweighted rows is
  refused, as in creghdfe.
- **Fixed effects:** y - Xb on the raw columns less the structural residual,
  less its weighted mean (the constant). Then `hdfe_fe_recover`, with the
  constant added to the first effect as 4a saves it. The connected-set
  split rules and accumulation are creghdfe's; this is the third copy of
  that code (§15).
- **Ado** (`build/civreghdfe.ado`):
  - `vce(bootstrap, ...)` replaces vce() by the base VCE: cluster on the
    bootstrap's `cluster()`, otherwise robust. The diagnostics, the
    first-stage and reduced-form regressions and e(V_modelbased) use it,
    and a note under the table names it.
  - `fevce(bootstrap, ...)` with `sestub()`, `nreps`, `dse()`; the names
    `__fese#__`, `v_se` and `stub#` are checked first.
  - The replicates run once, in the final pass after the collinearity
    rerun.
  - e() as Stata's vce(bootstrap) posts it: z-based, no e(df_r);
    e(chi2), e(p) and e(chi2type) replace e(F), and the header shows the
    Wald chi2 and the replications.
  - The fevce() e() results and the normalization warning, as creghdfe's.
  - Refused (r(198)) with vce(bootstrap): `bw()`, `kernel()`, `dkraay()`,
    `kiefer`. With any bootstrap: `b0()`, `smatrix()`, `wmatrix()`, and
    GMM2S or CUE with a kernel or two-way S.
- **Validation** (`validation/stata/bootstrap_reference.do`):
  - **Oracle.** `ivreghdfe` cannot run under the `bootstrap` prefix: the
    prefix's held estimates are gone when it finishes ("last estimates not
    found", r(301)), even through a wrapper. So the oracle draws bsample's
    counts with `bsref_draws`, expands the data and runs ivreghdfe on each
    replicate (with `idcluster`-style copy ids for clustered GMM2S).
  - `civreghdfe_bootstrap`, 23 checks: draws(stata) equals that oracle for
    - 2SLS on observations (5e-17);
    - LIML with clusters (4e-9, LIML's eigenvalue);
    - GMM2S with clusters (4e-17, the copy-aware S);
    - GMM2S with aweights (5e-18);
    - Fuller with partial(), strata and clusters (4e-15);
    - CUE with pweights on a well-identified model (1e-5 tolerance).
  - It also checks:
    - fweighted clusters against the expanded data's bootstrap (exact);
    - draws(fast) reproducibility;
    - e() under both base VCEs, with e(V_modelbased) and the diagnostics
      equal to the plain fit's;
    - `estat bootstrap` and replay;
    - the collinearity rerun;
    - the refusals.
  - `civreghdfe_fe_bootstrap`, 8 checks:
    - per-level SEs, d's SE and replicate counts equal those of civreghdfe
      `savefe` on each expanded replicate, with one and two effects;
    - shared replicates with vce(bootstrap);
    - names, e() and refusals.
  - `validation/option_evidence.json` has civreghdfe `fevce`.
  - All bootstrap cases, civreghdfe (463), grouped syntax (17), savefe (8)
    and oct27 (9) pass. The native estimators harness compiles and links
    for Linux; it was not run here.
- **creghdfe fix found here:** with fweights and rows resampled
  (`unit_freq`), creghdfe weighted a replicate row by draws times fweight.
  Its replicates were then centered on the fweight-squared estimator:
  - with fweights of 1 and 30 on auto, the replicate mean was -396 with
    e(b) -302 (the fweight-squared estimate is -526);
  - the SE was 31.7, against 42.6 on the expanded data.
  - The row now weighs the observations drawn. A `creghdfe_bootstrap` check
    (16 now) fails on the old plugin and passes on the new.
- **Found, not fixed:**
  - civreghdfe's CUE optimizer can stop at a worse local optimum than
    ivreghdfe's: two of ten replicates of a weakly identified model had
    J = 21.9 against 18.2, and 35.9 against 26.5. The full fit on those
    expanded data does the same, so this is the CUE optimizer, not the
    bootstrap. (Fixed after Phase 4d; see below.)
  - cppmlhdfe resamples fweighted rows as single units, which deviates from
    §5.3. (Fixed after Phase 4d, see below; the singleton counts with
    clusters were right.)
- **Not in 4c:** `smatrix()` and `wmatrix()` with the bootstrap; `saving()`;
  the FE `ci` option; a native test of the copy-aware S (the oracle covers
  it); sharing the connected-set code among the three commands.

**Phase 4d delivered: civreghdfe `wildbootstrap()` (WRE).**

- **What is computed.** boottest's WRE, read from `boottest.mata` 4.5.5
  (`InitVars`, `Estimate`, `MakeResiduals`, `PrepWRE`, `MakeWREStats`):
  - The DGP's structural fit is always LIML (`DGP.LIML = 1`), restricted
    under the null (an OLS fit when no endogenous regressor is left) and
    unrestricted with `nonull`, whatever the user's estimator.
  - The efficient first stage regresses all endogenous regressors on the
    instruments and that residual u; U2 = Y2 - X Pi keeps u's part.
  - Y2* = X Pi + U2 v and y* = X1 g + Y2* b + u v, with the same v.
  - Each replicate applies the user's k-class estimator (LIML's kappa
    recomputed; Fuller with the fit's constant alpha/(N_eff - K_iv)) and
    its sandwich, with P_X Z in the middle as ivreg2 and civreghdfe compute
    it. The statistic is (b*_k - b0)/se*; with `nonull`, (b*_k - b^_k)/se*.
- **Engine** (`src/civreghdfe/civreghdfe_wre.{c,h}`, new, no SPI). Every
  column is a combination of the basis [y, Y2, X] (X = [X1, X2]):
  - The DGP at b0 needs only the basis's cross products (LIML's
    eigenproblem, the first stage by FWL on u), and gives the replicate's
    columns as X Theta + Q C, with Q the perturbed basis.
  - A replicate needs X'WQ, Q'WQ and, for the sandwich, each cluster's
    X_g'W_g Q_g. On the nested path these come from the sums kept at setup
    per bootstrap cluster and per (cluster x bootstrap cluster) group, or
    from the rows when the groups are small (robust); a point then costs
    O(B (H L nb + K^3)). On the projected path the perturbed basis is
    partialled by the fixed effects once per replicate, and its cluster sums
    (or rows) are kept when they fit in a quarter of memory; otherwise they
    are projected again at each evaluation.
  - X'WQ and Q'WQ of every replicate are kept (B (L nb + nb^2) doubles).
    Discrete weights are kept as one byte each (doubles for normal and
    gamma), and 909 stops a run whose weights do not fit.
  - The observed statistic is replicate 0 (v = 1) through the same code, so
    that a replicate with every weight 1 ties with it exactly.
  - LIML's lambda is the smallest root of |A - lambda B| after scaling by
    diag(B)^(-1/2), as civreghdfe computes it: closed form for one or two
    endogenous columns, Jacobi beyond.
  - The interval search is `hdfe_wild_invert_batch` (33-point grid, bisection
    to 1e-6 se), evaluating each batch of points in one pass over the
    replicates.
- **Model columns.** From the retained fit (4c's retain, now also N_eff and
  the kept instrument count): the fixed-effect-partialled columns, the fit's
  omitted regressors and instruments dropped. `partial()` variables return
  as exogenous regressors (redundant ones left out by weighted
  Gram-Schmidt), because each replicate's perturbed residuals must be
  partialled on them as on any regressor; the results then do not depend on
  `partial()` (checked).
- **Plugin and ado.** `wild_setup`, `wild_chunk` and `wild_finish [keep]`
  in phase 1 as creghdfe's, the weights from `_ctools_hdfe_boot wilddraw`,
  e(wb_*) as creghdfe's with e(wb_method) `wre` or `wue`, the table by
  `_ctools_wild_table` (now labelling WRE and WUE), `nowildtable`, and the
  refusals of §4.6 (unadjusted VCE, `vce(bootstrap)`, `gmm2s`, `cue`,
  kernels, two-way clusters, `b0()`, `smatrix()`, `wmatrix()`, robust with
  fweights and no `bootcluster()`, an endogenous regressor collinear with
  the instruments).
- **Phase 0 correction: boottest after ivreghdfe with `absorb()` is not a
  WRE reference.** With firm effects nested in firm clusters, boottest's
  statistic for the exogenous coefficient was t = 29.9 against the
  regression's 2.0, and the endogenous coefficient's p-value differed
  (13/512 against 16/512). boottest after `ivreg2`, `ivregress` or
  `ivreghdfe` without `absorb()`, with the effects as dummies, gives 16/512
  and 52/512, as do civreghdfe and a direct Mata implementation. The
  reference is therefore boottest after `ivreg2` with dummies. boottest also
  refuses k-class and Fuller fits ("Only for use with OLS, 2SLS, ...").
- **boottest's interval bounds** stop at the first point whose p-value is
  within one replicate step of the level (`search()`), so they differ from
  the exact jump, by 5e-3 se in the 9-cluster example. With LIML, boottest's
  p-value function also steps up to about 1e-4 se away from the exact one.
  civreghdfe's bounds are the exact jumps (checked against the Mata WRE to
  1e-6 se).
- **Validation:** `civreghdfe_wild` in `validation/stata/bootstrap_reference.do`,
  31 checks:
  - p-values equal boottest's after `ivreg2` with dummies (1e-12) in 12
    designs, with interval bounds checked through boottest's p-value just
    inside and just outside each bound (1e-5 se; 1e-3 se with LIML):
    - nested Rademacher, crossed (projected) Webb, LIML Mammen, aweights
      with normal weights, robust observation-level (projected) gamma;
    - two endogenous regressors (also with LIML on the projected path),
      equal-tailed p-values, two absorbed effects, pweights, `nonull`;
    - `partial(wx)`, which equals boottest with wx as a regressor.
  - Enumeration (9 clusters, 512 vectors, RNG state unchanged) equals
    boottest.
  - Fuller(1), k-class(0.8) and coarser `bootcluster()` (8 groups) equal the
    replicate-by-replicate Mata WRE (`_bsref_wre_t`, in the support file)
    on the same draws, and the LIML two-endogenous bounds are exact against
    it.
  - Projected and nested paths agree on a nested design; `fedgp(partialled)`;
    e(b), e(V) and e(j) unchanged; `coefs()`, `noci`, `nowildtable`, replay;
    the WUE label; with `fevce(bootstrap)` each as alone; the refusals.
- **Fix found here:** `wildbootstrap` without options failed in all three
  commands ("option wildbootstrap not allowed"): `syntax` leaves a bare
  `wildbootstrap` among the other options when it is declared with an
  argument. The option helpers now look for it there in a second pass; a
  `civreghdfe_wild` check covers creghdfe, cppmlhdfe and civreghdfe.
- **Speed:** 10,000 rows, 200 clusters, 999 replicates, two tested
  coefficients with intervals: 0.3 s nested, 0.5 s projected.
- **Checks run:** every bootstrap case, civreghdfe (463), grouped syntax
  (17), the creghdfe and cppmlhdfe wild cases, and the options, docs, seeds,
  architecture and release checks pass. The new engine compiles without
  warnings (`-Wall -Wextra`, C11) for Windows; the OpenMP build was not run
  here, and there is no native test of the engine yet (the Mata WRE covers
  it in Stata).
- **Not in 4d:** a native `wre` engine test; tests of hypotheses on more than
  one coefficient (boottest's multi-restriction WRE); `bootcluster()`
  coarser than the vce() clusters on the projected path beyond the Mata
  check; `smatrix()`, `wmatrix()`, GMM2S and CUE (§15).

**CUE optimizer fix (after Phase 4d).**

- **Cause.** civreghdfe's own objective at ivreghdfe's estimate was lower
  (J = 18.229 against the 21.934 it stopped at), so the objective agreed and
  the optimizer did not. From the 2SLS start the numerical Hessian was not
  positive definite, and the fallback was a normalized steepest-descent
  step: steps of about 0.008 of a coefficient unit, for all 100 iterations,
  stopping far from the minimum without a message.
- **Fix** (`ivest_compute_cue`, `civreghdfe_estimate.c`). Where Cholesky
  fails, the step is Newton's on the Hessian with its eigenvalues replaced
  by their absolute values (floored at 1e-8 of the largest), in the
  coefficients' natural units (rms(y)/rms(x_k)), with a cyclic Jacobi
  eigensystem; steepest descent remains only for a zero Hessian. Where the
  Hessian is positive definite nothing changes.
- **Results.**
  - The ten bootstrap replicates now all equal ivreghdfe's J and b.
  - On 114 weakly identified simulated fits (robust, cluster and HAC CUE),
    the old optimizer differed from ivreghdfe in 8, 7 of them at a worse J;
    the new one differs in 1, where its J is lower than ivreghdfe's by 4e-11
    on a flat objective.
  - New checks: three nonconvex-start fits against ivreghdfe in
    `validation/stata/civreghdfe.do` (robust, cluster, HAC; six significant
    figures, since ivreg2's optimize() stops at its own tolerance on these
    flat objectives: b within 1e-7, civreghdfe's J lower by 1e-14), and the
    weakly identified model's bootstrap against the expanded-data oracle in
    `civreghdfe_bootstrap`. All four fail on the old plugin and pass on the
    new.
  - civreghdfe (466), civreghdfe_bootstrap (24), the wild, savefe and
    grouped-syntax cases, `s2_iv` (290), `nov04` (13) and the repo checks
    pass.

**cppmlhdfe fweight resampling fix (after Phase 4d).**

- **Bug.** Without `cluster()`, cppmlhdfe drew fweighted rows as single
  units and weighted each by draws times fweight, against §5.3. With
  fweights of 1 and 30 on auto (`price mpg weight, absorb(foreign)`), the
  bootstrap SE was 0.0287 against 0.0078 on the expanded data: 3.7 times
  too large, not too small as the 4c note guessed.
- **Fix** (`cppmlhdfe_boot.c`), as creghdfe and civreghdfe do it: the frame
  sets `unit_freq` to the fweights, so the engine draws the Σ fw
  observations, and a row's replicate weight is the count drawn. The SE is
  now 0.0076 (0.98 times the expanded data's, within Monte Carlo error).
  `draws(stata)` with fweighted rows is refused, in the ado and the
  plugin, as in civreghdfe.
- **Singletons with clusters were right.** ppmlhdfe and cppmlhdfe count an
  fweighted row once for the singleton rule (a lone row with fw = 3 is
  dropped; reghdfe and creghdfe keep it), and §5.4 applies the command's own
  rule to the drawn data, so each drawn copy counts once, as before.
- **Checks** (`cppmlhdfe_bootstrap`, now 17): fweighted clusters equal the
  bootstrap prefix on the expanded data (draws(stata), same seed); fweighted
  rows centered on e(b) with SEs within 20% of the expanded data's (fails on
  the old plugin); the draws(stata) refusal. The cppmlhdfe FE and wild cases
  pass.

**Phase 5 delivered: the refused options, docs and plumbing.**

- **`saving(filename [, replace double every(#)])`** in `vce(bootstrap)` of
  all three commands (in `fevce()` when it runs alone).
  - `boot_chunk ... rows` stores the chunk's coefficients
    (`hdfe_boot_store_rows`) in matrix `__<cmd>_bs_rows`; the helper's `run`
    posts them to a temporary postfile (`_r1 ... _rK`, double); after e(b)
    is posted, `savefinish` writes the file as `bootstrap` does: one
    variable per e(b) term, named by `_prefix_expand`'s rule (`_b_name`, or
    `_bs_#`), its labels and characteristics (`_dta[N]`, `[N_cluster]`,
    `[N_strata]`, `[seed]`, `[bs_version]`, `[observed]`, ...), floats
    unless `double`. A map of each e(b) term to its replicate column (0 for
    omitted and base terms) comes from `_ctools_ebv` on diag(1..K).
  - The chunk is capped at `c(matsize)` rows when saving.
- **The hidden `_bsdraws(varlist)`**: one multiplicity variable per
  replicate, read per chunk (`hdfe_boot_read_draws`, supplied mode; the
  rows of a unit must agree, nonnegative integers), not at retain as §11.1
  first said. `e(draws)` is `supplied`; `draws()` and `_bsdraws()` may not
  be combined.
- **draws(stata) with fweighted rows** (all three commands):
  `hdfe_boot_map_stata` maps bsample's uniforms on the fweight-expanded
  frame, in `expand`'s order (within a stratum, each row's first copy, then
  each row's further copies, row by row); `__<cmd>_bs_U` is then Σ fw. The
  helper caps a chunk at about 1000 variables of uniforms and stops with
  r(900) when one replicate needs more.
- **A draws(stata) bug found and fixed: strata with observations.** The
  prefix calls `bsample` without `weight()`. Its stratified branch
  (`StrSRSWR`) runs `by strata: assert` on the data `expand` left unsorted,
  the `by` fails, and it draws U more uniforms to reorder the sample (which
  does not change it). `draws(stata)` matched the first replicate only; the
  helper now takes and discards those uniforms (`uniforms ... strataobs`),
  and the replicates and the final RNG state equal the prefix's. (Review
  2026-11-08, BOOTSTRAP_REVIEW.md M4: `expand` leaves the data sorted when no
  row was drawn twice, and `bsample` then takes no extra uniforms; the helper
  now takes them only when a stratum drew one of its units twice.) Clusters
  within strata (`StrClustSRSWR`) were already right, and `bsample,
  weight()` (bsref_draws) takes no extra uniforms.
- **`fevce(bootstrap, ci)`**: `hdfe_boot_acc_keep_fe` keeps every
  replicate's effects as floats (L x R, level-major; r(909) above the memory
  budget), and `hdfe_boot_fe_ci` gives `_pctile`'s bounds at `level()` over
  the replicates that estimate a level. Outputs `__felb#__`/`__feub#__` (or
  `stub#_lb`/`_ub`, `v_lb`/`v_ub`), e(fe_lb), e(fe_ub), e(fe_level).
- **Heterogeneous slopes** with `vce(bootstrap)`, `fevce()` and
  `wildbootstrap()`:
  - creghdfe: the retained model keeps each slope factor's variable and the
    fit's slope map; replicates (and the wild projection) set the slope
    factors up with `creghdfe_setup_slopes` under the replicate weights,
    and recover effects with `creghdfe_fe_recover_slopes`.
  - cppmlhdfe: the same slope factors (centered when the fit centers them;
    `ppml_update_fe_weights` computes the centers); slope-only factors are
    not screened for separation, as in the fit; `ppml_fe_minnorm` already
    handled slopes. The wild restricted fits get them too.
  - Connected sets and the normalization warning concern the intercepts; a
    slope is estimable in a level where the replicate's slope variable
    varies (nonzero without an intercept); d sums intercepts and slopes
    times their variable. Output names follow savefe's: `__fese#__Slope#`.
  - Found while testing (not changed): with slopes and a second absorbed
    effect, creghdfe's full-sample savefe intercepts differ from reghdfe's
    by a constant per connected set (another normalization; d agrees), and
    cppmlhdfe's `d()` with slopes differs from ppmlhdfe's by up to 0.07 on
    the test data. Both predate this work.
- **civreghdfe savefe and fevce() with partial()**: the ado recovers the
  effects with the partialled regressors as regressors of creghdfe (FWL);
  each replicate keeps the partialled columns' cross products before its
  FWL step and solves for their coefficients in `civb_effects`. The effects
  and SEs equal those of the model without partial() (1e-13).
- **Shared code**: `src/hdfe/hdfe_boot_stata.{c,h}` (draw reading, row
  storing); the chunk loop moved from the three ados into
  `_ctools_hdfe_boot run`, which loads `ctools_plugin` itself (a plugin
  program is local to the ado that defines it).
- **Validation**: the `bootstrap_options` case (25 checks): saving() against
  the prefix's file (names, labels, characteristics, replicates; bstat),
  double/every()/chunks/failed replicates, noconstant, cppmlhdfe and
  civreghdfe; strata and fweighted draws(stata) against the prefix (with
  RNG states); _bsdraws(); ci bounds against reghdfe, ppmlhdfe and
  civreghdfe on the same draws; slopes (e(V), intercept/slope/d SEs, slopes
  alone, wild against the slopes as regressors on both paths, cppmlhdfe);
  partial(). The old refusal checks were replaced. Every bootstrap case,
  creghdfe (385), cppmlhdfe, civreghdfe (466), s2_hdfe (64), s2_iv (290)
  and the regression cases pass.
- **Docs**: the help files and command pages, COMPATIBILITY.md, option
  evidence (ci), and a developer-guide section in
  `docs/development/README.md`.
- **Plan text superseded by the code**: §6.2's `hdfe_build_factors_on_levels`
  and `hdfe_components` were not built (replicates call
  `hdfe_build_factors` with the full fit's level counts, and each command
  runs its own union-find); §6.3's guard checks `SF_nobs()` only, not the
  sort order; §6.5's nested parallelism exists (`hdfe_boot_run` gives
  leftover threads to each replicate); §6.8's chunks are a multiple of the
  replicates in flight near 50 (at most 400), and the dots show `.` or `x`.
- **Not done**: the protocol bump to 3 in `config/release.json` (deferred to
  the release, see Phase 2b), the §14 decisions (recorded there as
  implemented, for confirmation), benchmarks, native tests of
  `hdfe_boot_map_stata`, `hdfe_boot_fe_ci` and the WRE engine, and an
  OpenMP run.

Phases 2b–2d ship creghdfe before PPML and IV are ready. Each later phase only
adds a `model_ops` implementation, or the restricted-fit and aggregate
callbacks for the wild bootstrap. Phase 1w can run in parallel with 2a–2c, but
2d needs both.

## 11. Validation and benchmarks

### 11.1 Exact oracle: same draws, independent estimator

A hidden option `_bsdraws(varlist)` takes R multiplicity variables m1…mR instead
of the internal RNG. The plugin reads each chunk's variables in `boot_chunk`, in
place of the draws (Phase 5). The reference side exists since Phase 0: `bsref_draws` in
`validation/support/bootstrap_reference.do` reproduces the prefix's draws. On
StataNow, `rwgen bsample` plus `bootstrap, fweights()` is an independent
second oracle. The comparison case extends `validation/stata/bootstrap_reference.do`
and does the following:

1. `set seed`, then generate m1…mR with `bsref_draws` (or `bsample, weight(m_r)`),
   with `cluster()` and `strata()` as needed.
2. For each r, with the same weight type and the multiplicity folded in (fweight:
   `[fw = m_r · fw]`; aweight: `[aw = m_r · aw]`; and so on):
   - creghdfe → `reghdfe ... if m_r > 0, absorb(...)`;
   - cppmlhdfe → `ppmlhdfe`;
   - civreghdfe 2SLS and LIML → `ivreghdfe`;
   - civreghdfe GMM2S with cluster → `expand m_r`, give each copy a new cluster
     id, and run `ivreghdfe` (the idcluster form, because the fweight shortcut is
     not equivalent there).
3. Store the replicate's b and FE estimates (`absorb(fe1=id)` or savefe on the
   expanded or weighted data), and compute V, the SEs and the FE SEs in Mata.
4. Compare with `<cmd> ..., vce(bootstrap, _bsdraws(m1-mR)) fevce(bootstrap)`.
   Tolerances: replicates 1e-8 relative, V 1e-8, FE SE 1e-7 (the FE solver
   tolerance).

**Prefix parity through `draws(stata)`.** With the same seed,
`<cmd> ..., vce(bootstrap, reps(R) seed(s) cluster(c) draws(stata))` must reproduce
`bootstrap _b, reps(R) seed(s) cluster(c) idcluster(nc): <reference> ...`:

- every replicate's b (through `saving()` on both sides), to 1e-8;
- e(V), to 1e-8;
- c(rngstate) after both commands, exactly.

This covers observation and cluster bootstraps, strata, and each weight type
the prefix accepts. It also checks that `draws(fast)` and an
analytic fit leave c(rngstate) as §5.8 states (two draws, and none).

### 11.1b Wild bootstrap oracles

- **Manual loop, exact.** For each tested coefficient, Stata builds the
  restricted fit (reghdfe without x_k; for WRE, the restricted structural and
  efficient first-stage regressions), then y* (and Y₂*) from the same weights v
  the command used. The weights are written out by the hidden `_wbdraws()`
  option, or regenerated from e(rngstate). It then runs reghdfe or ivreghdfe on
  each replicate with the user's analytic `vce()`. The test compares t*_b per
  replicate (1e-8) and the p-values exactly, on both the nested and the
  projected paths.
- **boottest, the reference.**
  - Rademacher enumeration with C ≤ 12 is deterministic on both sides. Compare
    p-values and intervals exactly against boottest after reghdfe (one
    absorbed effect), on nested and non-nested designs, with the default
    `fedgp(full)`. For IV the reference is boottest after `ivreg2` with the
    effects as dummies (Phase 4d: boottest after `ivreghdfe` with `absorb()`
    does not reproduce the WRE statistic).
  - Non-enumerated runs, for every weight type, are compared with the same
    seed. This is the Phase 1w gate for the draw mapping (§5.8).
  - Two or more absorbed effects use boottest after `regress` or `ivreg2` with
    dummies.
  - PPML (cppmlhdfe) uses boottest after `poisson` with dummies.
- **Stata's `wildbootstrap`** (Stata 18 and later): `fedgp(partialled)
  noenumerate ptype(equaltail)` against `wildbootstrap areg` or `regress` on
  conventions only, within Monte Carlo error, since the draws differ.
- **Fast path against projected path.** A hidden switch forces the projected
  path on a nested design. The results must agree to 1e-10.

### 11.2 Statistical oracle

Compare `bootstrap _b, reps(1000) seed(s) cluster(c) idcluster(nc): reghdfe ...`
with `creghdfe ..., vce(bootstrap, reps(1000) seed(s) cluster(c))`. Both seeds
are fixed, so the test is deterministic. The SEs must agree within the Monte
Carlo band: the relative difference between two independent bootstrap SEs has a
standard deviation of about 1/√R. Use 4/√R, about 13% at R = 1000. Run the same
comparison for ppmlhdfe and ivreghdfe.

### 11.3 Native tests

Add these to `validation/native/estimators.suite.c` and a new `boot.suite.c`:

- **Draws (`draws(fast)`):**
  - draws for (seed, r) do not depend on call order, chunking or threads;
  - counts sum to n_s in every stratum;
  - cluster draws give equal multiplicity within a cluster;
  - fweight row draws sum to Σ fw.
- **Binomial sampler:** mean and variance χ² tests; edge cases n = 0, p ∈ {0, 1}.
- **FWL warm start:** partialling x̃ under random weights equals partialling x to
  1e-10.
- **Factors on full-sample levels:** zero-weight levels leave the projection
  unchanged.
- **Accumulators:** Welford against a two-pass computation; percentile and BC
  against hand-computed values.
- **Component split detection:** constructed bipartite graphs.
- **SPI guard (debug build):** calling an SPI wrapper from a worker aborts the
  test.
- **Wild kernels:**
  - the aggregate kernel's t*_b against a direct dense computation on toy data,
    for every weight type;
  - the projection kernel against the aggregate kernel on a nested design;
  - the five-coefficient form of t*_b(β0) against direct recomputation at
    random β0;
  - enumeration covering all 2^C sign vectors;
  - `hdfe_wild_invert` on a case with a closed-form bound, an unbounded case,
    and a disjoint case.

### 11.4 Stata suites

Add these to `validation/stata/{creghdfe,cppmlhdfe,civreghdfe}.do`:

- **Determinism:** identical e(V) and SE variables with `threads(1)` and
  `threads(8)` (within 1e-12), with chunk sizes forced to 1, 7 and all, and
  across reruns with the same seed (bitwise).
- **Modes:** coefficients only, FE only (e(V) equals the analytic run), and both
  (e(V) equals coefficients-only, and the FE SEs equal FE-only).
- **Edge cases:**
  - 2 clusters;
  - strata with one cluster;
  - cluster equal to an absorbed effect;
  - heavy singletons;
  - a dummy that is all zero in some replicates (N_misreps, e(reps));
  - all-singleton replicates;
  - non-convergence;
  - factor variables with bases;
  - the refit path (a level only on singleton rows);
  - noconstant;
  - each weight type.
- **Errors (§4.6):** each restriction; `fevce(bootstrap)` without savefe; output
  name collisions; a sort between calls (hidden fault hook).
- **Break and failure:** fault hooks like the existing `$CTOOLS_creghdfe_fault`
  (`outputs`, `post`) plus a new `boot` hook that fails mid-chunk. The data,
  c(changed) and e() must be as before, and the retained state must be released.
  `ctools` memory reporting, or a second command that would fail its phase check,
  shows the release.
- **Compatibility:**
  - `bootstrap: <cmd>` still works;
  - replay;
  - `estat bootstrap, all`;
  - `predict` after a bootstrap fit.
- **Existing test:** replace the error-match test at
  `validation/stata/creghdfe.do:3328`. reghdfe still rejects `vce(bootstrap)`,
  so the new test checks creghdfe's results instead.
- **Seeds check:** `validation/checks.py` (`seeds_*`) must treat `vce(boot…)`,
  `fevce(boot…)` and `wildbootstrap(` without `seed()` as random draws that need
  a preceding `set seed`. Enumeration needs no seed, but the check stays simple
  and asks for one anyway.
- **RNG behavior:**
  - `set rngstream 1` and `set rngstream 2` give different draws;
  - `set rng kiss32` under `version 13` works;
  - `seed()` equals `set seed` before the command;
  - an analytic fit leaves c(rngstate) unchanged;
  - after an enumerated wild bootstrap, c(rngstate) is unchanged.
- **Wild bootstrap:**
  - enumeration and the Webb note at small C;
  - an unbounded interval with C = 4;
  - `nonull`;
  - `coefs()`;
  - `bootcluster()` coarser than the `vce()` cluster;
  - observation-level wild bootstrap under `vce(robust)`;
  - each weight type;
  - `wildbootstrap()` combined with `fevce(bootstrap, ...)`;
  - each error in §4.6;
  - replay of the wild table.
- **Normalization warning (§4.8):**
  - G = 1: no warning, e(fe_norm_warning) = 0;
  - G = 2 and G = 3: the warning text, naming the effects and the SE variables;
  - the `dse()` variant of the last sentence;
  - the split-replicate line on a thinly connected design;
  - replay shows the warning in the footnote;
  - cppmlhdfe slope-only effects do not count toward G.

### 11.5 Benchmarks

Add profiles to `validation/bench_profiles.py` and `validation/performance/`:

- **Comparison:** `bootstrap, reps(100) cluster(): reghdfe` and
  `bootstrap ...: creghdfe` against `creghdfe, vce(bootstrap, reps(100) cluster())`
  on the HDFE datasets (`docs/performance/HDFE_DATASETS.md`), at threads 1, 4 and
  8.
- **Reported:**
  - per-replicate time relative to one full fit;
  - CG iterations, warm against cold;
  - FE-target overhead;
  - PPML IRLS iterations, warm against cold;
  - peak memory.

- **Wild bootstrap:**
  - creghdfe and civreghdfe `wildbootstrap(reps(9999))` against boottest, on
    nested (`absorb(firm) cluster(firm)`) and projected
    (`absorb(firm year) cluster(firm)`) designs;
  - time per replicate as C, K and N grow;
  - the projection batch size chosen by the budget.
- **Draw modes:** `draws(stata)` against `draws(fast)`, for cluster and
  observation bootstraps, to show the cost of exact parity.

Targets to confirm, not promises:

- one replicate costs at most 0.5× a full creghdfe call at equal threads;
- end to end, at least 5× faster than `bootstrap: creghdfe` at 8 threads;
- a nested-path wild bootstrap costs about the same for N = 10⁵ and N = 10⁷ once
  the aggregates are formed.

## 12. Documentation and release plumbing

- **Command pages and help:** `docs/commands/{creghdfe,cppmlhdfe,civreghdfe}.md`
  and `build/{creghdfe,cppmlhdfe,civreghdfe}.sthlp`:
  - the options, including `wildbootstrap()` and `draws()`;
  - the modes (§4.1);
  - a "Remarks on fixed-effect standard errors" section with the G ≥ 2 caveat
    (§5.7), the warning text and when it appears (§4.8), and why `dse()` is the
    identified alternative;
  - random numbers: what `set seed`, `seed()`, `set rng` and `set rngstream` do,
    how many draws each mode takes, and when `draws(stata)` matches `bootstrap`
    (§5.8);
  - the wild bootstrap: method per command, nested and projected paths, few
    clusters and Webb weights, and intervals by inversion;
  - weights semantics;
  - e() results;
  - examples, including `estat bootstrap` and the wild table.
- **Compatibility:** `docs/COMPATIBILITY.md`.
  - reghdfe, ppmlhdfe and ivreghdfe do not offer `vce(bootstrap)` or a built-in
    wild bootstrap; ctools adds them.
  - `draws(fast)` draws differ from `bsample`'s; `draws(stata)` matches them.
  - Wild results match boottest where Phase 0 confirmed its conventions.
  - Weights are allowed, although Stata's `vce(bootstrap)` refuses them, and
    fweights are resampled as units.
  - The wild bootstrap follows boottest: its defaults, its weight draws for the
    same seed, and the full-model treatment of effects that are not nested. It
    differs from Stata's `wildbootstrap`, whose conventions `fedgp(partialled)
    noenumerate ptype(equaltail) reps(1000)` reproduce, though not its draws.
    It also goes beyond boottest's scope: several absorbed effects, and PPML
    after the command itself.
- **Option coverage:** `docs/OPTION_COVERAGE.md` and
  `validation/option_evidence.json` get the new options and evidence.
- **Developer guide:** `docs/development/README.md`, a section on the bootstrap
  engine and the phase protocol for the estimators. Also
  `docs/development/COMMAND_INVENTORY.md`, regenerated by
  `scripts/sync_architecture.py` if `ctools_commands.def` changes.
- **Release:** `config/release.json` protocol 3, stamped by
  `scripts/sync_release.py`.
- **This document:** update it with the Phase 0 findings and decisions, and
  remove it, or turn it into the developer-guide section, when the work ships.

## 13. Risks

| Risk | Mitigation |
|------|------------|
| Splitting 1,400–2,000-line estimator functions changes numerics | Bitwise e() gates per suite before any bootstrap code (Phases 2a, 3, 4b); no algorithm edits during the moves |
| SPI called from a worker thread (crash or corrupted state) | SPI-free `fit()` contract, config structs, debug-build thread assertion |
| Memory: retained model plus T concurrent workspaces | Explicit per-replicate estimate, budget, fallback to the inner mode, r(909) before work starts, `ci` cap |
| FE SEs with G ≥ 2 are normalization dependent and can mislead | A runtime warning whenever per-level SEs are written with G ≥ 2 (§4.8), pointing to `dse()`; estimability masks and counts; help-file remarks |
| Replicate-level separation in PPML is slow to detect | Measure in Phase 3; separation checks stay on by default; failed replicates are counted, not silently kept |
| Users expect `bsample`-identical draws | `draws(stata)`, the default, reproduces them; `draws(fast)` is documented as for cbsample |
| boottest's weight draws are not reproduced exactly (scaling, matrix layout, memory-cap splitting) | Phase 1w gates on exact p-value equality with boottest for the same seed on non-enumerated designs; enumerated designs need no draws |
| Default `fedgp(full)` wild bootstrap is slow when an effect is not nested and B is large (k+1 projections per replicate) | `coefs()` to test fewer coefficients; batched projections; a note pointing to `fedgp(partialled)`; nested designs (the common `absorb(firm) cluster(firm)`) stay on the O(C·K) path |
| Long runs that cannot be interrupted | Chunked calls with Break between chunks, chunks sized to about a second |
| Ado reorders or changes observations between calls | `SF_nobs` and sort-order guard with r(459) |
| Few-cluster designs give poor pairs-bootstrap inference | A warning that suggests `wildbootstrap()` |

## 14. Open questions

### Decided

- **Wild cluster bootstrap: in scope, for coefficients only.** It is never
  applied to saved fixed effects (§4.7, §5.9, G8).
- **FE target with two or more effects.** Each replicate is normalized by the
  estimator's own rule, with estimability masks (§5.7). Per-level SEs are still
  written, and the command always warns about the normalization caveat and
  points to `dse()` (§4.8, G9). Anchoring to a full-sample normalization is
  later work.
- **Random numbers.** Stata's RNG drives every draw and is treated as Stata's
  `bootstrap` treats it, with `draws(stata)` for exact parity (§5.8, G5).
- **The wild bootstrap's reference is boottest.**
  - Effects not nested in the clusters default to the full-model bootstrap,
    `fedgp(full)`, with `fedgp(partialled)` as the faster option that matches
    Stata's `wildbootstrap areg` (§5.9).
  - Defaults are boottest's: 999 reps, symmetric p-values, enumeration when
    2^C ≤ reps (`noenumerate` to draw instead), and boottest's weight draws for
    the same seed (§4.7, §5.8).
  - This settles former questions 6 and 8.

### Still open

Each question below carries what Phase 0 found and a recommendation. The code
implements one answer to each (noted as **Implemented**), to be confirmed
before the release.

1. **Option names.** `vce()` + `fevce()` + `wildbootstrap()` (§4.1, §4.7), the
   new `fedgp()` (§5.9), and the SE variable names (§4.3).
   - *Implemented:* as named there; slope outputs follow savefe's names
     (`__fese#__Slope#`), bounds `__felb#__`/`__feub#__`.
2. **Default draws mode for the pairs bootstrap.**
   - *Found:* `draws(stata)` is exact and needs no fallback. Its cost is one
     Mata `runiform()` and a transfer of 2·U doubles per replicate (§5.8). That
     is negligible for cluster bootstraps and about two column loads for
     observation bootstraps.
   - *Recommendation:* make `draws(stata)` the default, so results reproduce
     `bootstrap … : cmd`. Keep `draws(fast)` for large observation-level runs,
     and benchmark the gap in Phase 2b before fixing it.
   - *Decided (2026-11-06):* `draws(stata)` is the default, as recommended
     (it also covers fweighted rows and observations within strata);
     `draws(fast)` is opt-in. A chunk's uniforms are capped at about 2^26
     values (512 MB), at least one replicate.
3. **Default reps for the pairs bootstrap.**
   - *Found:* Stata's `vce(bootstrap)` uses 50. The wild bootstrap's 999 is
     decided above.
   - *Recommendation:* 50, with a note recommending more for percentile
     intervals and FE SEs.
   - *Implemented:* 50.
4. **Weights.**
   - *Found:* Stata refuses every weight type with `vce(bootstrap)` and with
     the prefix (r(101)), unless the prefix gets `force`, which resamples rows.
     There is no Stata convention to match.
   - *Recommendation:* allow weights as a documented extension. Resample
     fweights as units (§5.3) in both draw modes; under `draws(stata)` that is
     `bsample` applied to the fweight-expanded data.
   - *Implemented:* as recommended, in both draw modes (Phase 5).
5. **civreghdfe base VCE.** Default it from the bootstrap design (§9.3), or add a
   suboption such as `vce(bootstrap, base(robust | cluster c))`? No new
   evidence.
   - *Implemented:* the default from the design, without a suboption.
6. **cppmlhdfe wild bootstrap.** The score bootstrap is PPML's wild analog
   (§8.4, Phase 3w). Keep it in this project, or ship creghdfe and civreghdfe
   first? boottest cannot run after ppmlhdfe; the oracle is boottest after
   `poisson` with dummies (check 15), so it can only be validated on small
   designs.
   - *Implemented:* kept (Phase 3w).
7. **Precision.** Summarize replicates in double (more accurate, about 1e-7 from
   Stata's default output), or round them to float as the prefix does (bitwise
   parity)? Recommendation: double, with float rounding only behind a hidden
   validation switch.
   - *Implemented:* double; `saving()` writes floats unless `double`, as the
     prefix does.

## 15. Later work

- **Wild bootstrap extensions:**
  - joint (multi-coefficient) hypotheses;
  - multiway clustering;
  - GMM2S and CUE in civreghdfe;
  - a `test()`-style option for arbitrary linear restrictions, as boottest has.
- **Anchored FE normalization** for G ≥ 2 (§5.7).
- **Jackknife.** Delete-one-cluster replicate weights through the same
  `model_ops`. This also enables BCa intervals.
- **Bootstrap-t and percentile-t,** which need replicate VCEs.
- **Combined replicate and column parallelism** with nested OpenMP regions.
- **Raw-norm CG stopping rule** for warm-started columns (§6.7).
- **FE percentile bounds without keeping the draws:** a second deterministic pass
  over the replicates. The draws can be regenerated from (seed, r) under
  `draws(fast)`, or by replaying Stata's RNG from e(rngstate) under
  `draws(stata)`.
- **`savefe` with `partial()` in civreghdfe** (§9.1).
- **reghdfe's `absorb(newvar=absvar)`** syntax in creghdfe and civreghdfe.

## Appendix A. Phase 0 findings

**Setup.** Recorded on StataNow 19 (MP), Windows 11, with reghdfe, ivreghdfe and
ppmlhdfe installed, and boottest 4.5.5 from SSC. Stata 14.1, the ctools
minimum, was not available on this machine.

**Evidence.**

- `validation/stata/bootstrap_reference.do` (case `bootstrap_reference`) checks
  the findings marked ✓ against Stata's own commands and boottest. All 19
  checks pass; the boottest checks skip when it is not installed.
- `validation/support/bootstrap_reference.do` holds `bsref_draws`, the Mata
  emulation of `bsample`.
- Findings without ✓ come from reading Stata's ado and Mata source in
  `ado/base/b`, `ado/base/_` and `ado/base/w`.

### Pairs bootstrap (`bootstrap` prefix, `vce(bootstrap)`, `bsample`)

1. ✓ **Replicates.** The prefix's replicates are R successive `bsample` calls
   after `set seed`, with nothing else drawing in between (`_loop_bs.ado`). It
   leaves c(rngstate) where those calls do.
2. ✓ **Frame.** The prefix resamples e(sample) of the full fit
   (`_prefix_check4esample` → `keep if e(sample)`). reghdfe's singletons are
   therefore never drawn.
3. ✓ **Emulation.** `bsample`'s draws are reproduced exactly in Mata (§5.8) for
   observation, strata, cluster and strata-plus-cluster designs, and so is the
   final RNG state. `regress … [fw=m_r]` on the emulated draws reproduces each
   prefix replicate to 1e-12; `reghdfe` does too, to 1e-9.
   - `bsample` is `version 8`, but its `uniform()` calls follow the caller's
     RNG (mt64 here).
   - Behavior under `version 13` (kiss32) is untested.
4. ✓ **Replay with supplied weights (StataNow).** `rwgen bsample` plus
   `bootstrap, fweights()` replays the prefix's replicates. With a cluster, the
   replay expands clusters first. This is a second exact oracle, available
   only in StataNow.
5. ✓ **Weights.** Every weight type is refused with `vce(bootstrap)` and with
   the prefix (r(101) "weights not allowed with bootstrap estimation"), unless
   the prefix gets `force`.
6. ✓ **Failed replicates.** A replicate whose omitted-coefficient pattern differs
   from the full sample is posted as missing in every statistic, with the note
   quoted in §5.5. e(V) is casewise over complete replicates, and e(reps) is the
   same for every coefficient.
7. ✓ **Float storage.** Replicates are posted as floats unless `saving(…,
   double)` is given, and e(V) and the intervals are computed from the posted
   values. `mse` gives Σ(b^r − b̂)²/R.
8. **Posted results and display.**
   - e() list and display: §4.4 and §4.5. The table headers are "Observed
     coefficient", "Bootstrap std. err.", z, "P>|z|" and "Normal-based [95%
     conf. interval]".
   - BC uses #{b^r ≤ b̂} (`_bs_sum.ado`), with a `ties` variant that counts half
     the ties.
   - reghdfe itself refuses `vce(bootstrap)` with r(9).

### Wild bootstrap (Stata's `wildbootstrap`, Stata 18 and later)

9. ✓ **Absorbed effects.** For an effect nested in the clusters, `wildbootstrap
   areg` equals `regress` with dummies. For an effect crossing the clusters,
   `areg` equals a wild bootstrap of the group-demeaned regression: p-value and
   bounds agree to 1e-7. Dummies give a different answer: p = 0.156 against
   0.162 in the check's data. Its convention is therefore the partialled
   regression (§5.9).
10. ✓ **Defaults:** 1000 reps, Rademacher weights, equal-tailed p-values for one
    coefficient (the F upper tail for several).
11. **Rules from the source** (`wildboot_regress.mata`,
    `wild_cluster_bootstrap.ado`):
    - p-values use strict inequalities with denominator B;
    - no enumeration: with 6 clusters it draws 1000 vectors and warns about ties;
    - reps is raised so that α·reps is an integer when intervals are requested;
    - weights are drawn in blocks of min(1000, reps) from a reps × clusters
      matrix, with an undocumented `btsequence` option for boottest's order;
    - only `regress`, `areg` (one absorbed variable) and `xtreg, fe` are
      supported, so there is no IV or ML reference.
    - Results are e(wboot) (estimate, t, p, lb, ub), e(wb_pci), e(N_wbreps),
      e(wb_weight), e(wb_ptype) and e(wb_rseed).

### boottest 4.5.5 (SSC)

12. **Defaults** (`boottest.ado`, `boottest.mata`):
    - 999 reps, Rademacher weights, symmetric p-values;
    - p = #{t* beyond t} / B with strict inequalities, over the feasible
      replicates;
    - enumeration of all 2^C Rademacher vectors when 2^C ≤ reps, which reports
      r(reps) = 2^C.
13. ✓ **Absorbed effects.** For an effect crossing the clusters, boottest after
    `reghdfe` or `areg` equals `regress` with dummies, so it bootstraps the
    **full model** (p = 0.414 in the check's data, with the demeaned
    regression at 0.391). That is the opposite of Stata's `wildbootstrap`
    (check 9). For a nested effect the two conventions coincide, and so do
    the two programs.
14. **Scope.**
    - ✓ After `reghdfe` and `ivreghdfe` it accepts only one absorbed effect:
      r(198) for two.
    - GMM is no longer supported.
    - ✓ WRE runs after `ivreghdfe` for endogenous and exogenous coefficients,
      for LIML, and with two endogenous regressors. (Phase 4d: with
      `absorb()` its statistic is not the WRE's; after `ivreg2` with the
      effects as dummies it is. boottest refuses k-class and Fuller fits.)
15. **PPML.** boottest cannot run after `ppmlhdfe`. The null-imposed score
    bootstrap needs `constraints()`, and `nonull` needs `predict, score`;
    `ppmlhdfe` offers neither. It does work after `poisson … i.fe`, which is
    the PPML reference for small designs.
16. **Draws.** boottest draws a clusters × (reps + 1) matrix. With the same seed
    its p-values differ from `wildbootstrap`'s, with or without the
    undocumented `btsequence`. The two programs do not share draws.

**Consequences for the plan.**

- The two references disagree on effects not nested in the clusters. The
  decision (§14, Decided) makes boottest the reference: `fedgp(full)` is the
  default and matches boottest; `fedgp(partialled)` matches Stata's
  `wildbootstrap`.
- boottest serves as an exact oracle for one absorbed effect: enumerated runs
  without draws, and non-enumerated runs with the same seed once the draws are
  reproduced (Phase 1w).
  For two or more effects the oracle is the manual loop (§11.1b), or boottest
  after `regress` or `ivreg2` with dummies on small designs.
- The PPML oracle is boottest after `poisson` with dummies.

### Not done in Phase 0

- **Stata 14.1 runs.** Checks 1–8 should be repeated under the minimum version,
  or under `version 13`/`version 14.1` for the RNG question in check 3.
- **The full validation runner.** It stops before running cases until the
  validation data are fetched (`scripts/fetch_validation_data.py`). The
  reference case was run directly instead:
  `do "validation/stata/bootstrap_reference.do" bootstrap_reference`, from the
  repository root.
