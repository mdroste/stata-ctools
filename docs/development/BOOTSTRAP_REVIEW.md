# Bootstrap review: bugs found in the new bootstrap code

Review date: 2026-11-08. It covers the uncommitted working tree on top of
`cccfe84`: the pairs bootstrap (`vce(bootstrap)`, `fevce(bootstrap)`) and the
wild bootstrap (`wildbootstrap()`) in creghdfe, cppmlhdfe and civreghdfe. The
findings were then fixed the same day: the High and Medium ones, and then the
Low ones; see [Fixes applied](#fixes-applied-2026-11-08). The plan items that
are not implemented are still open.

**How it was done.** The shared parts were read line by line:

- the engine, `src/hdfe/hdfe_boot.c`;
- the draw reader, `hdfe_boot_stata.c`;
- the binomial sampler in `ctools_rng.h`;
- the plugin phase table in `ctools_runtime.c`;
- the ado helper, `_ctools_hdfe_boot.ado`.

They were checked against Stata's own `bsample.ado`, `_bs_sum.ado`,
`_loop_bs.ado` and `bootstrap.ado` (StataNow 19). Four reviewers then traced
the code for each command and for the wild-bootstrap engines.

Every High and Medium finding was checked again by reading its code path.
Finding M4 was also confirmed by running Stata. Low findings that were not
re-checked are marked "reviewer-traced".

## Summary

| # | Severity | Command | Finding |
|---|----------|---------|---------|
| H1 | High | all (wild) | The confidence-interval search can loop forever inside the plugin |
| H2 | High | civreghdfe | `fevce()` with a `partial()` variable absorbed by the fixed effects gives fixed effects of about 1e14 |
| H3 | High | creghdfe | Slopes-only models: the replicate `_cons` is not zero, so `fevce()` and `dse()` SEs are wrong |
| H4 | High | civreghdfe | GMM2S/CUE with `center`: the centering in the bootstrap's moment covariance S is not scaled by the multiplicity, so e(V) is wrong |
| M1 | Medium | creghdfe, cppmlhdfe | A slope that does not vary in a replicate still leaves its level's intercept marked estimable |
| M2 | Medium | cppmlhdfe | One replicate's separation-check failure aborts the whole bootstrap |
| M3 | Medium | creghdfe | A singleton refit together with `fevce(bootstrap)` always fails with r(110) |
| M4 | Medium | all (pairs) | `draws(stata)` + `strata()` with observations resampled: the draws get out of step with `bsample` |
| M5 | Medium | cppmlhdfe, civreghdfe (wild) | A failed evaluation in the CI search counts as a rejection, so a false finite bound is reported |
| M6 | Medium | civreghdfe | GMM2S failing in a replicate calls the SPI on a worker thread and silently uses 2SLS |
| M7 | Medium | all (pairs) | `fevce(bootstrap, level() mse dots())` is silently ignored when `vce(bootstrap)` is also given |
| M8 | Medium | creghdfe, civreghdfe (wild) | The projected path stops with r(909) on large data, whatever the RAM |
| L1 | Low | all (pairs) | Under `mse`, e(vcetype) is "Bootstrap", not "Bstrap *", and e(mse) is missing |
| L2 | Low | cppmlhdfe | e(p) is not posted; the normalization warning prints when no SE variable exists |
| L3 | Low | civreghdfe | `partial()` with more than 64 columns plus `fevce()` fails with an unexplained r(198) |
| L4 | Low | civreghdfe | Replicate S depends on the thread layout; nested thread oversubscription |
| L5 | Low | civreghdfe | Error paths after `ereturn post` leave e() and variables inconsistent |
| L6 | Low | all | State left behind after Break or error |
| L7 | Low | cppmlhdfe | Replicate details: simplex weights, memory estimate, 10-pass cap |
| L8 | Low | all (wild) | Small robustness gaps in the wild bootstrap |
| L9 | Low | all (pairs) | The displayed interval level may not follow e(level) (check against Stata) |
| L10 | Low | all (pairs) | fweight totals above INT_MAX |

Also: plan items not implemented, and what was checked and found correct.

## Fixes applied (2026-11-08)

All findings are fixed in the working tree. "Caught" means the new test was
also run against a plugin with that fix reverted and failed there.

| # | Fix | Test |
|---|-----|------|
| H1 | `outward()` and `bound_query()` stop on the stored step (`dist >= limit`), with a cap of 64 doublings as a backstop (`hdfe_wild.c`) | native `wild_engine`: unbounded sets at the two failing (center, scale) pairs, single and batch searches. Caught: the old stop rule hangs |
| H2 | `civb_effects` drops the `partial()` columns that the replicate's FWL step dropped (`isp == 2`) before its Cholesky (`civreghdfe_boot.c`) | `bootstrap_options`: `partial()` of a variable constant within year equals the model without it. Caught: off by about 9e12 |
| H3 | Replicate `_cons` is 0 without an absorbed intercept, as in the full fit (`creghdfe_boot.c`) | `bootstrap_options`: slopes-only `fevce()` slope and `dse()` SEs from reghdfe on the same draws. Caught: 0.13 vs below 1e-6 |
| H4 | `omega_groups` centers each copy: `- mult * md` (`civreghdfe_omega.c`) | `civreghdfe_bootstrap`: GMM2S `center` against the expanded data, for clusters (idcluster) and for observations. Caught: V off by 8e-6 vs below 1e-9 |
| M1 | Where the full fit estimates a level's slope but a replicate cannot, that replicate's intercept for the level is NaN too (`creghdfe_boot.c`, `cppmlhdfe_boot.c`; full-fit estimability kept at retain) | none new; existing slope tests (cluster draws) unchanged |
| M2 | A failed simplex or ReLU check makes the replicate `HDFE_BOOT_NOCONV`; only 920 aborts (`cppmlhdfe_boot.c`) | not reproduced in a test |
| M3 | The singleton refit drops `fevce()`'s SE, nrep, lb and ub variables before rerunning (`creghdfe.ado`) | `bootstrap_options`: `fevce()` with a factor level only on singleton rows |
| M4 | At retain, the plugin saves the strata of two or more draw units (`hdfe_boot_save_dup_strata`). `_ctools_boot_uniforms` takes the U extra uniforms only when one of them drew a unit twice; with more than 64 such strata it always takes them, since the alternative has probability below 2^-64 | `bootstrap_options`: three strata of two rows, 50 replicates, e(V) and c(rngstate) equal the prefix's |
| M5 | A NaN p counts as not rejected (so a failure can only widen the set) and is counted in `hdfe_wild_ci.failed`. `e(wb_cistatus)` gets a fifth column, `failed`, and the tables print a note | native `wild_engine`: NaN p above 3 leaves the upper side open and flagged, single and batch searches |
| M6 | `IVEstContext.report` guards GMM2S's messages; in a replicate (`report == 0`) a GMM2S failure fails the replicate instead of falling back to 2SLS (`civreghdfe_estimate.c/.h`) | not reproduced in a test |
| M7 | `_ctools_hdfe_boot parse` returns `r(given)` (level, mse, dots). With `vce(bootstrap)`, `fevce()` may repeat those options but not change them: r(198) (all three ados) | `bootstrap_options`: conflicting `level()` and `mse` refused; repeated or omitted accepted |
| M8 | Callers pass a quarter of physical memory without the 2/4 GB clamp. The engines batch at `HDFE_WILD_BATCH_BYTES` (2 GB) or `CIV_WRE_SOFT_BYTES` (4 GB), and use more only for one replicate that needs it (`hdfe_wild.c/.h`, `civreghdfe_wre.c/.h`, `creghdfe_boot.c`, `civreghdfe_boot.c`) | not reproduced (needs about 2 GB of columns per replicate) |

**Low findings:**

| # | Fix | Test |
|---|-----|------|
| L1 | With `mse`, the three `_bootstrap_post` programs post `e(vcetype)` `Bstrap *` and `e(mse)` `mse`, as `_bs_sum`; `ereturn display` then drops the "Normal-based" header, as Stata's table does | `bootstrap_options`: all three commands |
| L2 | cppmlhdfe posts `e(p)` = `chi2tail(e(df_m), e(chi2))` under `vce(bootstrap)` (not otherwise, as ppmlhdfe). `e(fe_norm_warning)` needs per-level SE variables | `bootstrap_options`: `e(p)`; `d()` alone gives no warning |
| L3 | `civb_effects` allocates its `partial()` solve per replicate (`pgw`, `pdrop` in the replicate buffer), so the 64-column limit is gone | `bootstrap_options`: `partial()` of 69 dummies with `fevce()` equals the model with them as regressors |
| L4 | A replicate's moment covariance is summed serially (`omega_hetero` and the AC lags when `om->mult`), and CUE runs one thread in a replicate | existing GMM2S and CUE bootstrap tests; not a dedicated thread-layout test |
| L5 | civreghdfe clears e() when FE recovery fails or is interrupted, and drops this run's effect and SE variables when `savefinish` fails | not reproduced in a test |
| L6 | `_ctools_hdfe_boot run` and the three wild chunk loops run inside `capture noisily`: any error or Break drops the uniform or weight variables, the rows matrix, the postfile and the Mata draws. creghdfe's error branch also runs `_ctools_cleanup` and the drops. cppmlhdfe releases the retained fit on its rebuild exit and when its output variables cannot be created | not reproduced in a test (needs a Break) |
| L7 | cppmlhdfe's simplex check gets positive-row weights mult × (fweight or 1), as on the expanded data (`ppml_simplex` takes `w_pos`; the full fit passes NULL). The memory estimate counts the simplex's two copies of X, (4K + 24) row vectors. The pass cap is gone: passes run until no rows are dropped (at most N1 + 1) | existing cppmlhdfe bootstrap tests |
| L8 | The tables note empty sets (`e(wb_cistatus)[., 4]`), and civreghdfe marks an interval it cannot search as empty. The score model counts its weight rows (`got`): in order, and all before the statistic. The WRE checks B = 2^C when enumerating. Enumeration needs C ≤ 30. `WildParse` reads `coefs()` and `bootcluster()` before removing `noci`. cppmlhdfe posts `e(wb_score)`, the score statistic the bootstrap studentizes, and the table shows it (`e(wb_t)` stays the Wald z, as documented) | native `wild_engine`: weight order and completeness. `bootstrap_options`: `coefs(noci)` |
| L9 | Checked against `regress, vce(bootstrap ...)`: the bootstrap's `level()` sets e(level), e(ci_*) and the table at estimation (a replay shows c(level)); the command's own `level()` sets the bootstrap's; two different levels are r(198). New helper `_ctools_hdfe_boot level`, used by all three commands | `bootstrap_options`: the two spellings give the same e(ci_normal), conflicting levels are refused, equal ones accepted |
| L10 | A unit stands for at most INT_MAX / 2 observations, but totals above INT_MAX are valid: draws(fast) draws them. `hdfe_boot_draw_units` returns -1 and `hdfe_boot_draw_total` gives the double total, posted as `__<cmd>_bs_U`, so draws(stata) stops in the ado with its "use draws(fast)" message | native `boot_engine`: a 3e9 total, its draws summing to it, and the unit limit |

**Validation** (Windows plugin built with zig, no OpenMP; StataNow 19 MP):

- Native `boot_engine` and `wild_engine` pass.
- The bootstrap Stata cases pass: `bootstrap_reference`, `bootstrap_helpers`,
  `bootstrap_options` (35 checks with the Low fixes' tests), `wild_draws`, and
  the `*_bootstrap`, `*_fe_bootstrap` and `*_wild` cases of the three
  commands, plus `civreghdfe_savefe`.
- The main suites `creghdfe`, `civreghdfe`, `s2_iv` and `s2_hdfe` pass.
- `cppmlhdfe` passes 1216 of 1217 checks. The failure, `native upstream hard5`,
  is a full-sample fit that none of these fixes reach. A plugin built earlier
  the same day gives identical results: coefficients of about 1e-17 that the
  significant-figure check rejects.
- Not run: the OpenMP builds, and Linux and macOS (only compile-checked for
  macOS).

Help files, `docs/commands/*.md` and the design notes describe the changed
behavior: `fevce()` options, the flat-slope intercept, e(wb_cistatus), the
stratified uniforms, and for the Low fixes `mse`'s labels, `level()`, cppmlhdfe's
`e(p)` and `e(wb_score)`.

---

## High

### H1. The wild-bootstrap confidence-interval search can loop forever (all three commands)

**Where:** `src/hdfe/hdfe_wild.c:745` (`outward`, used by creghdfe) and
`src/hdfe/hdfe_wild.c:864` (`bound_query`, used by civreghdfe and cppmlhdfe).

```c
const double limit = 1000.0 * scale;
for (;;) {
    if (fabs(x_acc - center) >= limit) return NAN;   /* :745 */
    double next = 2.0 * dist;
    if (next > limit) next = limit;
    const double x = center + side * next;
    if (!accepted(p, ctx, x, alpha)) return bisect(...);
    x_acc = x;  dist = next;
}
```

**What happens:**
- Once the step is capped, the query point is `fl(center + side*limit)`.
- If that point is accepted, the next stop test recomputes the distance as
  `fl(x_acc - center)`. Rounding can make this a hair below `limit`.
- The loop then asks for the same point again, accepts it again, and never
  ends.
- This runs inside a plugin call, so Break cannot stop it. In
  `hdfe_wild_invert_batch` the outer `for (;;)` (around line 957) never reaches
  `n == 0`.
- Each pass re-evaluates p: cheap for WCR, but a full WRE evaluation in
  civreghdfe and a full restricted PPML fit in cppmlhdfe.

**Verified numerically (IEEE double):**
- `center = 0.3808`, `scale = 0.064`: `limit = 64`, but
  `(0.3808 + 64) - 0.3808 = 63.99999999999999`.
- `center = 0.598104`, `scale = 0.1278` also hangs on the upper side.
- Over random (center, scale) pairs, about 0.1–0.6% per side are affected.

**When the search reaches the cap:**
- Always with `ptype(lower)` or `ptype(upper)`: one side is then unbounded.
- Whenever the confidence set is unbounded on a side: few clusters, or weak
  instruments with WRE.
- Intervals are on by default in creghdfe and civreghdfe.

**Fix:** test the stored step, which is set to exactly `limit` when capped,
not the recomputed difference:

```c
/* outward() */
for (;;) {
    if (dist >= limit) return NAN;
    ...
}
/* bound_query() */
if (b->dist >= limit) { b->result = NAN; b->active = 0; return 0; }
```

Add an iteration cap (for example 64 doublings) as a backstop. Add a native
test with `center = 0.3808`, `scale = 0.064` and a p-function that accepts
everything.

### H2. civreghdfe: `fevce(bootstrap)` with an absorbed `partial()` variable gives fixed effects of about 1e14

**Where:** `src/civreghdfe/civreghdfe_boot.c:643-668` (`civb_effects`), used
at `:702`.

```c
ST_double dj = WW[j * Kp + j];
const ST_double scale = dj;                 /* relative to its own diagonal */
...
if (dropped[j] || !(dj > 1e-12 * scale) || !(scale > 0.0))
```

**What happens:**
- `WW` holds the cross products of the partialled W columns
  (`b.pw`, built at `:1005-1035` from the FE-partialled `b.A`).
- The tolerance is relative to the column's own partialled diagonal. So it only
  drops W columns that are collinear with earlier W columns.
- A column absorbed by the fixed effects is left as rounding noise
  (W̃ ≈ 1e-16·|W|). WW[j][j] is then about 1e-32·|W|²·N, so the column is not
  dropped.
- γ_j = rhs/WW comes out enormous. Line 702 (`d -= X_exog * gamma`) puts W·γ
  into d, and the recovery folds it into α.
- The replicate's own FWL step flags the column (`b.isp[k] == 2`,
  `civ_fwl_partial_out` at `:1054`), and so does the full fit. But
  `civb_effects` never reads that flag. The comment at `:643-644` ("the fit's
  FWL step drops them the same way") is wrong.

**Scenario:** `civreghdfe y (x = z) w female, absorb(id, savefe) partial(female) fevce(bootstrap)`
with `female` constant within id.
- Correction found while testing the fix: with id the only or first absorbed
  effect, such a column is partialled to exactly 0, and `!(scale > 0.0)` already
  dropped it. The bug needs rounding noise: a column absorbed by a later effect
  or by a combination, for example constant within year in `absorb(idcode
  year)`. With that column, the fixed effects were off by about 9e12.
- Every replicate's `__fese#__` and `dse()` come out around 1e14.
- `partial(w)` with w varying in only a few clusters breaks the cluster
  replicates that drop those clusters.
- The coefficients are unaffected. No test covers `fevce` with `partial()`.

**Fix:** seed `dropped[]` from the replicate's FWL flags before the
factorization:

```c
for (ST_int k = 0, p = 0; k < m->Kx; k++)
    if (m->partial[k]) dropped[p++] = (b->isp[k] == 2);
```

Better still, use the FWL step's own absolute rule: diagonal ≤ `absorb_tol`
times the column's pre-partial sum of squares. Then this step and the
coefficient path cannot disagree.

### H3. creghdfe: slopes-only models give wrong `fevce()` and `dse()` standard errors

**Where:** `src/creghdfe/creghdfe_boot.c:877`, used at `:572`.

```c
out->b[ols.K_keep] = means[0] - xb;            /* :877, always */
...
ST_double e = m->raw[i] - out->b[Kk] - resid;  /* :572, replicate_effects */
```

**What happens:**
- With only slope effects (`absorb(id#c.z)`, so `no_intercept`), the full fit
  posts `_cons = 0`: `creghdfe_regress.c:997` and `:1386`
  (`beta_keep[K_keep] = no_intercept ? 0.0 : ...`).
- Its savefe therefore subtracts 0.
- Each replicate subtracts c_r = ȳ_r − x̄_r·b_r ≠ 0.
- With no intercept, `creghdfe_fe_recover_slopes` has no level mean to absorb
  c_r (`xbar` is filled only when `intercept >= 0`,
  `creghdfe_hdfe.c:289-299`).
- So every slope estimate is shifted by about −c_r·Σx/Σx², differently in each
  replicate.

**Scenario:** `creghdfe y x, absorb(id#c.z, savefe) fevce(bootstrap, dse(d))`.
- The per-level slope SEs are inflated by the spread of c_r.
- The `dse()` values are wrong in the same way.
- e(V) is unaffected: the ado drops `_cons` for slopes-only models
  (`build/creghdfe.ado:1006`).
- `bootstrap_reference.do:2919-2930` tests slopes-only models only with
  `vce(bootstrap)`.

**Fix:** use
`out->b[ols.K_keep] = m->n_intercepts > 0 ? means[0] - xb : 0.0;`.
`n_intercepts` is set in `full_sample_sets` (`:308`) before its early return.

### H4. civreghdfe: GMM2S/CUE with `center`, wrong cluster S in the bootstrap

**Where:** `src/civreghdfe/civreghdfe_omega.c:274` (`omega_groups`).

```c
ST_double wi = om->w ? om->w[i] : 1.0;     /* in a replicate: mult * w */
hd[g] += wi * (e[i] * z[i]) - md;          /* md subtracted once per row */
```

**What happens:**
- In the expanded data (and under `idcluster()`), each of a row's m copies
  subtracts the centering mean. The correct cluster sum is therefore
  h_g = Σ m_i (w_i u_i − md).
- The code subtracts md once per row.
- When clusters are resampled, `mult_groups` then divides the outer product by
  m. This gives (mΣwu − n_g·md)²/m instead of m(Σwu − n_g·md)².
- `omega_center` (`:187-205`) and the robust branch are correct.

**Scenario:** `civreghdfe ..., gmm2s center vce(bootstrap, cluster(c))`, or the
same with CUE.
- Every replicate's weighting matrix differs from the expanded-data S, so the
  replicate coefficients, e(V), e(b_bs) and the CIs are wrong.
- The size of the error grows with the centering term.
- `center` is no longer refused, and the bootstrap restrictions
  (`civreghdfe.ado:1114-1136`) do not exclude it.

**Fix:** use
`hd[g] += wi * (e[i] * z[i]) - (om->mult ? om->mult[i] : 1.0) * md;`.
Alternatively, refuse `center` with `vce(bootstrap)` and GMM2S/CUE until it is
validated against `bootstrap, idcluster(): ivreghdfe ..., center`.

---

## Medium

### M1. A slope that is flat in a replicate leaves its level's intercept marked estimable (creghdfe, cppmlhdfe)

**Where:**
- cppmlhdfe: `src/cppmlhdfe/cppmlhdfe_boot.c:524-531` sets only the slope
  node to NaN.
- creghdfe: `src/creghdfe/creghdfe_boot.c:617` skips only the slope, when the
  slope factor's `counts` (Σ(x − x̄_l)², `creghdfe_hdfe.c:132-145`) is 0.

**What happens:**
- With `absorb(id##c.year, savefe)`, take a replicate whose drawn rows for a
  level all share one x0, for example one row drawn twice, which survives the
  singleton rule.
- The slope cannot be estimated, so the recovery sets it to 0
  (`creghdfe_hdfe.c:346`; `cppmlhdfe_irls.c:134-145`).
- The raw intercept then equals the level mean, a_l + b_l·x0, not a_l.
- Only the sum is identified, yet the intercept is recorded as estimable.
- With a year-valued slope variable and b_l = 0.05, such replicates are off by
  about 100. Small levels hit this often in an observation bootstrap, and the
  intercept SE explodes without a warning.

**Fix:** where the slope is marked not estimable, also set to NaN the
intercept node of the same absorbed term:
- creghdfe: the factor h with `is_intercept[h]` and
  `slope_in.fe_absvar[h] == slope_in.fe_absvar[g]`;
- cppmlhdfe: the factor with matching levels, found as
  `cppmlhdfe_irls.c:137-141` finds it.

d for that level's rows then becomes NaN as well, which is right: it is not
identified either.

### M2. cppmlhdfe: one replicate's separation-check failure aborts the whole bootstrap

**Where:** `src/cppmlhdfe/cppmlhdfe_boot.c:737-742` and `:754-759`.

```c
int nsep = ppml_simplex(&options, &S, y, X, w, 2, N2, K_x, sep_m);
if (nsep < 0) { rc = -nsep; goto done; }        /* same for ppml_relu_ex */
```

**What happens:**
- Return codes:
  - `ppml_simplex` returns −430 at its iteration limit
    (`cppmlhdfe_simplex.c:127`).
  - `ppml_relu_ex` returns −(projection status), or −9010 under
    `relu_strict(1)` (`cppmlhdfe_relu.c:18, 234`).
- These become `fit()` errors, which abort the whole run (`hdfe_boot.h:107-109`
  reserves errors for that).
- The prefix would instead count such a replicate as failed.
- Nothing is printed, because `report` is 0 in replicates. The user gets a bare
  r(430) or r(9010) after all the work so far.

**Fix:** in both places, use
`if (nsep == -920) { rc = 920; goto done; } out->status = HDFE_BOOT_NOCONV; goto done;`.

### M3. creghdfe: a singleton refit together with `fevce(bootstrap)` always fails with r(110)

**Where:** `build/creghdfe.ado:671-674` (SE, nrep, lb and ub variables created
with savefe) and `:751-765` (refit).

**What happens:**
- The refit drops the residual, `groupvar()` and the savefe variables, but not
  `fe_sevars`, `fe_nrepvars`, `fe_lbvars` or `fe_ubvars`.
- The recursive `_creghdfe_impl` (`:763`) then runs `_ctools_newvars`
  (`:599`), which fails because `__fese1__` already exists.

**Scenario:** `creghdfe y i.cat x, absorb(id, savefe) fevce(bootstrap)`, where
a level of `cat` occurs only on singleton rows.

**Fix:** in the refit block, add
``if `__fe' capture drop `fe_sevars' `fe_nrepvars' `fe_lbvars' `fe_ubvars'``.
The `dse()` variable is created later (`:788`). Add a test.

### M4. `draws(stata)` + `strata()`, observations resampled: the draws get out of step with `bsample`

**Where:** `build/_ctools_hdfe_boot.ado:595`, in `_ctools_boot_uniforms`:
`if (extra) skip = runiform(units, 1)`.

**What happens:**
- `bsample`'s `StrSRSWR` (without `weight()`, as the prefix calls it) runs
  `capture by touse strata: assert ...` after `expand`.
- It takes U more uniforms only when that check fails, which happens when the
  expanded data are no longer sorted.
- `expand` adds rows, and so unsorts the data, only when some unit is drawn
  twice within a stratum. If every stratum's draw is a permutation, nothing is
  expanded, the data stay sorted, the check passes, and no extra uniforms are
  taken.
- The helper always discards U, so every later replicate is shifted.

**Verified in StataNow 19** (scripts in the session scratchpad):
- After `expand` with all multiplicities 1, `sortedby` is kept and the
  `by: assert` returns 0.
- With 3 strata of 2 observations, 30 of 200 seeds took 2U uniforms, not 3U.
- Predicting "extras are taken if and only if some stratum has a repeated draw
  value" matched Stata's actual RNG state for 300 of 300 seeds.

**Impact:** the chance that every stratum is a permutation is the product of
n_s!/n_s^n_s, so this only matters when all strata are tiny, for example matched
pairs with few strata. The existing test uses `strata(rep78)` on auto, where it
never triggers. The documented guarantee ("reproduces the replicates of
`bootstrap … : cmd`") is then silently broken.

**Fix:**
- The plugin posts the stratum sizes of the draw frame (expanded sizes with
  fweights), for example a 1 × S matrix `__<cmd>_bs_nstrata`.
- Mata draws the extras only when a stratum has a repeated draw value:

```mata
ud = runiform(units, 1); ut = runiform(units, 1)
/* store ud, ut as now */
if (extra) {
    dup = 0; lo = 1
    for (s = 1; s <= cols(ns) & !dup; s++) {
        n = ns[s]
        if (n > 1) {
            d = sort(floor(ud[|lo \ lo + n - 1|] :* n :+ 1), 1)
            dup = any(d[|2 \ n|] :== d[|1 \ n - 1|])
        }
        lo = lo + n
    }
    if (dup) skip = runiform(units, 1)
}
```

`floor(u*n + 1)` in double is the same value as `bsample`'s
`int(uniform()*_N + 1)`. Add a regression test: three strata of two rows,
`reps(50)`, comparing e(V) and c(rngstate) with the prefix.

### M5. Wild CI: a failed evaluation counts as a rejection, so a false finite bound is reported (cppmlhdfe, civreghdfe)

**Where:** `src/hdfe/hdfe_wild.c:975`, `bound_update(&b[who[k]], pq[k] > alpha)`
(and `accepted()` at `:716`).

**What happens:**
- The p-value is left NaN when `restricted_fit` fails
  (`cppmlhdfe_boot.c:1993-2001`), or when `wre_dgp_build` fails in civreghdfe.
- `NaN > alpha` is false, so the point counts as rejected.
- The search then bisects toward it and reports a finite bound with
  `open = 0` and no note. e(wb_fail) only reflects the b0 = 0 fit.

**Where it bites:** cppmlhdfe `wildbootstrap(ci)`, at the far outward points
(up to 1000·se), where the offset x_k·b0 is extreme and IRLS is most likely to
fail.

**Fix:** keep failures separate from rejections. Count them in the context,
add a status column (for example a fifth `e(wb_cistatus)` column, "failed"),
and print a note. At minimum, treat NaN as accepted in the outward phase, so
the side is reported as open instead of as a false bound.

### M6. civreghdfe: GMM2S failing in a replicate calls the SPI on a worker thread and silently uses 2SLS

**Where:**
- `src/civreghdfe/civreghdfe_estimate.c:69, 75, 129`: unguarded `SF_error`
  calls in `ivest_compute_gmm2s`.
- `:1475-1477`: the caller falls back to 2SLS and returns success.

**What happens:**
- A replicate whose S or X'ZWZ'X fails Cholesky calls `SF_error` from an
  OpenMP worker. The SPI is not thread-safe.
- The replicate is then posted as OK with 2SLS coefficients, mixed into the
  GMM2S distribution.
- This is rare: only for a non-positive pivot.

**Fix:**
- Pass `report` into `IVEstContext` and guard the three calls.
- When `!report` (a replicate), treat any `gmm_rc` as a failure
  (`est_rc = gmm_rc; goto est_cleanup;`), so the replicate is recorded as
  missing.

### M7. `fevce(bootstrap, level() mse dots())` is silently ignored when `vce(bootstrap)` is also given (all three commands)

**Where:**
- `build/creghdfe.ado:296-317`
- `build/cppmlhdfe.ado:427-448`
- `build/civreghdfe.ado:479-495`

**What happens:**
- With `__bs`, only `r(designopts)` (reps, seed, cluster, strata, draws) is
  compared.
- fevce's `level`, `mse` and `dots` are read only in the FE-only branch.
- So `vce(bootstrap) fevce(bootstrap, ci level(90))` writes 95% bounds and posts
  `e(fe_level) = 95`, with no message.

**Fix:** with `__bs`, stop with r(198) when fevce's `r(level)` differs from
`bs_level`, or when fevce sets `mse` or `dots()` ("specify level() once").
Alternatively, honor a separate FE level by passing `felevel=` to
`boot_finish`.

### M8. Wild bootstrap, projected path: r(909) on large data, whatever the RAM (creghdfe, civreghdfe)

**Where:**
- `src/creghdfe/creghdfe_boot.c:1597-1600` caps the budget at 2 GB.
- `src/hdfe/hdfe_wild.c:524-526` returns 909 when
  `N·(1 + n_test)·8 > budget`.
- civreghdfe caps at 4 GB (`civreghdfe_boot.c:1663`), and
  `civreghdfe_wre.c:684-686, 1362-1364` stop the same way.

**Scenario:**
- This is the default `fedgp(full)` with an effect not nested in the clusters
  (year effects with firm clusters).
- With N = 30M and 9 tested coefficients, creghdfe stops with "not enough
  memory for one wild bootstrap replicate", even on a 512 GB machine.

**Fix:** the u⊙v and r_k⊙v columns are independent projections, so partial
them in sub-batches when one replicate's columns exceed the budget. Or derive
the budget from physical memory without the fixed cap.

---

## Low

### L1. Under `mse`, e(vcetype) and e(mse) differ from Stata (all three commands)

**Where:** `build/creghdfe.ado:1664`, `build/cppmlhdfe.ado:1693` and
`build/civreghdfe.ado:2962` always post `ereturn local vcetype Bootstrap`.
e(mse) is never posted.

**What Stata does:** `_bs_sum.ado` returns `vcetype "Bstrap *"` and `mse mse`
under mse.

**Fix:** pass `mse()` to the `_bootstrap_post` programs, then post
`ereturn local vcetype = cond(`mse', "Bstrap *", "Bootstrap")`, and
`ereturn local mse mse` when mse is set.

### L2. cppmlhdfe: e(p) is missing, and the normalization warning can print with no SE variables

- `build/cppmlhdfe.ado:1223` posts e(chi2) but not e(p). creghdfe and
  civreghdfe post both, and plan §4.4 lists e(p). Add
  ``ereturn scalar p = chi2tail(`df_m', e(chi2))``.
- `:1601` sets ``fe_norm_warning = (`nint' >= 2)`` even when only `d()` is
  bootstrapped and no per-level SE variables exist. The warning then talks
  about variables that do not exist (plan §4.8). Use
  ``(`nint' >= 2 & "`sevars'" != "")``.

### L3. civreghdfe: `partial()` with more than 64 columns plus `fevce()` fails with an unexplained r(198)

**Where:** `civreghdfe_boot.c:610` (`if (Kp > 64) return 198;`), with stack
arrays at `:606` and `:617`.

**What happens:** the 198 is returned as a `fit()` error, so the first replicate
aborts the run with no message. This is reachable with `partial(i.year)` and 70
years.

**Fix:** allocate WW, rhs and gamma on the heap, or refuse at retain with a
message.

### L4. civreghdfe: replicate S depends on the thread layout; nested oversubscription (reviewer-traced)

**Where:**
- `civreghdfe_omega.c:334` (`omega_hetero`) uses the global
  `ctools_get_openmp_threads()` for a row-partitioned reduction.
- CUE uses `ctools_get_max_threads()` (`civreghdfe_estimate.c:448`).

**What happens:**
- The summation order, so the last bits of a GMM2S/CUE replicate, depends on
  whether the replicate runs nested. That depends on chunk composition, which
  breaks the engine's "results do not depend on threads" contract (about 1e-15
  relative).
- With nested levels on, each of T workers opens a full-size team (T × nt
  threads).
- CUE's per-thread residual buffers are not counted in `civb_fit_bytes`.

**Fix:** use one thread (or fixed-size blocks) when `om->mult` is set, and
`omp_get_max_threads()` in CUE.

### L5. civreghdfe: error paths after `ereturn post` leave e() and variables inconsistent (reviewer-traced)

- FE-recovery failure or Break (`civreghdfe.ado:2318-2324`) exits with e()
  already posted. e(V) is then the bootstrap V, while e(vce) and e(vcetype) are
  the base ones.
- A `savefinish` failure (`:2399-2403`) clears e() but leaves `__hdfe#__`,
  `__fese#__` and the `dse()` variables in the data.
- **Fix:** `ereturn clear` (or hold and restore e()) on those paths, and create
  the SE variables only after `savefinish` succeeds, or drop them on that path.

### L6. State left behind after Break or error (all three commands)

- The uniform variables `__ctools_bsu*` (draws(stata)) are dropped only on a
  plugin error or at the normal end (`_ctools_hdfe_boot.ado:298, 309`).
  - A Break in `Uniforms`, `post` or `_dots` leaves them in the data under the
    `fast` option, which undoes nothing.
  - The default undo mode rolls them back.
- The Mata matrix `_ctools_wild_v` (C × B doubles; gigabytes when
  `vce(robust)` draws per observation) is dropped only on success or a failed
  chunk call. The error and Break paths do not drop it, and neither does
  `_ctools_cleanup`.
- The error path in `creghdfe.ado:24-30` does not run `_ctools_cleanup`, so the
  `__creghdfe_bs_*` matrices and scalars stay until the next run.
- cppmlhdfe (reviewer-traced): the retained model is not released on the
  multiway "rebuild" exit (`cppmlhdfe.ado:931-940`) or if tempvar generation
  fails (`:960-988`). It stays resident until the next plugin call.
- **Fix:** in each command's error branch, run
  `capture drop __ctools_bsu*`, `capture _ctools_hdfe_boot wildclear`,
  `capture _<cmd>_bs_release` and `_ctools_cleanup <cmd>`.

### L7. cppmlhdfe replicate details (reviewer-traced)

- **Simplex weights:** the replicate passes `w = m·w_user` with `weight_type = 2`
  to `ppml_simplex` (`cppmlhdfe_boot.c:737`). With aw/pw, the full fit gives its
  positive rows weight 1, and the expanded-data equivalent is m. This changes
  only borderline separation decisions near the scale-dependent thresholds.
  Pass `mult · (fw ? fw : 1)`.
- **Memory estimate:** `fit_bytes` (`:1096-1101`) estimates N·(3K + 24)
  doubles. The simplex adds 2·K·N2, so the peak is about (4K + 15)·N, and the
  budget is overrun by about 25% at K = 50.
- **Pass cap:** `PPML_BOOT_PASSES 10` marks a replicate failed when it needs
  more screen-and-refit passes. The full fit has no cap.

### L8. Small robustness gaps in the wild bootstrap (reviewer-traced unless noted)

- **Empty or skipped confidence sets:** shown as two missing bounds with no
  note. The table never reads `e(wb_cistatus)[., 4]`
  (`_ctools_wild_table.ado:69-70, 94`), and civreghdfe's skip path
  (`civreghdfe_boot.c:2163-2164`) sets no status.
- **`cppmlhdfe_wild_finish`:** does not check that all B weight rows arrived,
  and `m->v` is uninitialized (`hdfe_wild.c:1030`). civreghdfe and creghdfe
  check `got == B` and `done == B`.
- **`civ_wre_weights` enumeration** (`civreghdfe_wre.c:502-507`) does not check
  that B == 2^C.
- **C = 31:** the ado enumerates when 31·ln 2 < ln(reps) + 1e-6, but the plugin
  supports C ≤ 30 and `total=2^31` overflows int parsing. Add
  `clusters <= 30` to the rule at `_ctools_hdfe_boot.ado:508`.
- **`noci` removal** (`_ctools_hdfe_boot.ado:422`): removes a regressor
  literally named `noci` from `coefs()`. Names like `nocity` are safe.
- **cppmlhdfe e(wb_t):** posts the Wald z (`_b/_se`, `cppmlhdfe.ado:1550-1553`),
  but the bootstrap p-value is for the score statistic. Label it or post the
  score statistic.

### L9. The displayed interval level may not follow e(level)

`vce(bootstrap, level(90))` posts e(level) = 90 and 90% e(ci_*). The tables
(`_creghdfe_display`, cppmlhdfe's `ereturn display` at `:9` and `:1360`) use
only the command-line `level()`, which defaults to c(level). Conversely,
civreghdfe's main `level()` does not reach e(ci_*). Check what
`regress, vce(bootstrap, level(90))` displays and posts, then make the tables
and e() agree with it.

### L10. fweight totals above INT_MAX (shared engine)

`hdfe_boot_design_check` limits the frequency total per stratum, but
`hdfe_boot_draw_units` (`hdfe_boot.c:186-194`) sums all strata into an
`ST_int`. Above INT_MAX the cast is undefined, and the value is posted as
`__<cmd>_bs_U`. Without strata, creghdfe reports only "Error in regression
(rc=198)", even for `draws(fast)`, which could handle it. Check the overall
total in the design check, and give a message.

---

## Plan items not implemented

None of these is wrong output, but each is a documented promise
(BOOTSTRAP_PLAN.md §4.1, §4.5, §4.6, §5.4):

- The FE-only note when `vce()` clusters but `fevce()` resamples observations.
- The r(198) for `wildbootstrap()` + `fevce(bootstrap)` without its own
  `reps()`. Today fevce silently uses reps(50) with observation resampling.
- The few-clusters (< ~30) warning with `vce(bootstrap)`.
- The pweights "naive bootstrap, not svy" note.
- The "(Replications based on C clusters in var)" header line, and e(N_clust),
  e(cluster), e(N_strata) and e(strata) as `bstat` posts them. Only
  e(cluster_bs) and e(strata_bs) exist.
- "Depvar explained by the fixed effects" flagged in e(boot_status):
  `depvar_explained` is computed at `creghdfe_boot.c:835` and dropped.
- e(command) should omit `vce()` (§4.4). It contains the full option.

## Checked and found correct

- **Engine (`hdfe_boot.c`):**
  - `bsample`'s observation and cluster mappings, including the permutation
    rule for clusters and the stratified variants, checked against
    `bsample.ado` line by line;
  - the fweight-expanded frame order;
  - per-replicate seeding independent of chunks and threads;
  - the scheduling window and the folding of results in replicate order;
  - Welford accumulation;
  - `_bs_sum` parity: casewise complete replicates; V about the mean with
    n − 1, or with mse about b̂ with n (`_bs_sum`'s `r(V)*(n-1)/n`); se and the
    normal CI from the SD about the mean even under mse; the `_pctile`
    definition; z0 from ≤ with missing at 0 or n; the BC percentiles.
- **`ctools_rng_binomial`:** inversion, and BTRD steps 0–3.4 with the Stirling
  table, match Hörmann (1993).
- **Draw reader (`hdfe_boot_stata.c`):** packing positions match
  `_ctools_boot_store`; `_bsdraws()` rows are checked for consistency within
  units.
- **The prefix's draw sequence:** `bsample` is called without `weight()`
  (`_loop_bs.ado:144`), and the `saving()` "observed" characteristic is written
  as `bootstrap.ado:840` writes it.
- **Phase protocol (`ctools_runtime.c`):** retain → chunk → finish, `wild_finish
  keep` → `boot_chunk`, reset on error, and `threads()` stripped before the
  operation token is parsed.
- **Per command (reviewers):**
  - multiplicity × weight for every weight type, with `unit_freq` only for
    observation resampling;
  - the singleton rule judged on copies;
  - collinearity patterns compared with the full fit;
  - destandardization and `_cons` (except H3);
  - FE recovery mapped back to full-sample levels, with NaN presets and
    connected-set split detection;
  - no SPI calls in `fit()` (except M6);
  - frees on every exit path;
  - unit order matching `sort strata cluster, stable`;
  - nothing sorts between retain and finish.
- **Wild bootstrap (reviewer):**
  - weight draws reproduce boottest's `MakeWildWeights` (types, draw order, the
    dropped first column);
  - enumeration order matches `count_binary`;
  - the p-value formulas match `boottest::getp` (strict inequalities,
    denominator = feasible replicates);
  - the WCR/WCU algebra and the intersection of vce() and bootstrap clusters;
  - the WRE restricted-LIML DGP and efficient first stage;
  - the score bootstrap's recentering.

## Suggested order of work

1. **H1:** small fix; prevents a hang that Break cannot stop.
2. **H2, H3, H4:** wrong numbers with no warning.
3. **M1–M3, M6, M7:** wrong numbers or hard failures in reachable
   configurations.
4. **M4:** parity guarantee; add the tiny-strata prefix test.
5. **M5, M8**, then the Low items and the plan items.

Each fix above names the test that would have caught it. None of these
configurations is covered today:
- `fevce` with `partial()`;
- `fevce` with slopes only;
- `fevce` with a singleton refit;
- `center` with the bootstrap;
- tiny strata under `draws(stata)`;
- `ptype(upper)` intervals.
