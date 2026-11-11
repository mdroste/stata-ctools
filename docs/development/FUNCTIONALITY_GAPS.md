# Functionality gaps: ctools vs. the commands it replaces

Audit date: 2026-11-04 (repo at `cccfe84`, release 1.2.4).

This document lists every place found where a ctools command lacks functionality that the command it
replaces has. It is a work plan: each gap has an ID, a severity, and a status so it can be fixed,
documented, or explicitly declined.

**In scope:** missing options and suboptions, syntax forms, weights, prefixes (`by:`, `mi`, `bayes:`...),
data types, stored results, postestimation, and output, plus options that ctools accepts but ignores or
misreads.

**Out of scope:** features ctools adds beyond the original, and numerical or solver-path differences that
do not change what a user can do. These are covered in [COMPATIBILITY.md](../COMPATIBILITY.md).

## How this was produced

Each ctools command's `syntax` statement, help file, and C implementation (`build/*.ado`, `build/*.sthlp`,
`src/`) was compared option by option with the reference implementation installed on the audit machine:

| ctools | Reference (version audited) | Location |
|---|---|---|
| `creghdfe`, `chdfe` | `reghdfe` 6.13.1; `hdfe` 3.2.7 (GitHub tag, not installed locally) | `~/ado/plus/r/` |
| `civreghdfe` | `ivreghdfe` + `ivreg2` | `~/ado/plus/i/` |
| `cppmlhdfe` | `ppmlhdfe` 2.3.3 | `~/ado/plus/p/` |
| `cqreg` | `qreg` (StataNow 19.5) | `C:/Program Files/StataNow19/ado/base/q/` |
| `cimport`, `cexport` | `import *`, `export *` (StataNow 19.5) | `ado/base/i/`, `ado/base/e/` |
| `cmerge`, `csort`, `cencode`, `cdecode`, `cdestring`, `csplit`, `csample`, `cbsample`, `cinterpolate` | `merge`, `sort`, `encode`, `decode`, `destring`, `split`, `sample`, `bsample`, `ipolate` | `ado/base/` |
| `cbinscatter` | `binscatter` 7.02; `binsreg` (not installed; checked from its documentation) | `~/ado/plus/b/` |
| `cpsmatch` | `psmatch2` 4.0.12 | `~/ado/plus/p/` |
| `cdistinct`, `crangejoin`, `crangestat` | `distinct` 1.2.1, `rangejoin` 1.1.3, `rangestat` 1.1.1 | `~/ado/plus/` |
| `cwinsor` | `winsor2`, `gstats winsor` (gtools) | `~/ado/plus/` |

The `cppmlhdfe` and `cqreg` findings were mostly confirmed by running StataNow 19.5 in batch mode. All
other findings come from reading the code. Items marked **Run** in the Status column still need a Stata run
to confirm; see [Open questions](#open-questions-need-a-stata-run).

### Severity

| Level | Meaning |
|---|---|
| **S1** | **Silently wrong.** ctools accepts the option but ignores or misreads it, or produces wrong output, with no error. Fix or reject first. |
| **S2** | **Missing feature.** Rejected with an error. Blocks an analysis the original supports. |
| **S3** | **Compatibility friction.** Abbreviations, legacy or undocumented syntax, programmer options, different error codes. Existing scripts break, but a workaround is easy. |
| **S4** | **Stored results or output.** Missing `e()`/`r()` items, display content, labels, or messages. |

"Documented" means [COMPATIBILITY.md](../COMPATIBILITY.md) or the command's help file already lists the gap.

---

## 1. Priority list: silently wrong behavior (S1)

Fix these first. Each one either changes results without warning or makes an option a no-op.

**Status, 4 Nov 2026:** all fixed or, for the three Excel items, verified not to be gaps. Regression
contracts are in case `nov04` of `validation/stata/regressions.do` (`python3 validation/run.py nov04`).
They compare with the replaced command where possible: `destring`, `ivreg2`/`ivreghdfe`, `reghdfe`
and `ppmlhdfe` output, and Stata's `import excel`.

| ID | Command | Gap | Status |
|---|---|---|---|
| DESTR-01 | `cdestring` | `ignore("chars", illegal\|asbytes\|aschars)`: the suboptions are not parsed. The whole string, including `, illegal`, is taken as characters to strip, so `"1e5"` loses its `e` and becomes 15. | Fixed: parsed with `_parse_ignore`; the plugin handles `asbytes` and `illegal` |
| HDFE-10 | `creghdfe` | `group(var)` (reghdfe's group-FE option) is matched by ctools' own `GROUPvar()` option. With a new name it silently creates a mobility-group variable instead. | Fixed: `group()`, `individual()`, `aggregation()` exit with r(198); `groupvar()` needs `groupv`, as in reghdfe |
| PPML-01 | `cppmlhdfe` | `verbose(-1)` is meant to silence output but turns on the per-iteration trace, because C tests `if (verbose)`. | Fixed: debug traces only at `verbose(1+)`; `verbose(-1)` is quiet |
| PPML-02 | `cppmlhdfe` | `log`/`nolog` are accepted and never read; there is no iteration log to show or hide. | Fixed: ppmlhdfe-style IRLS log by default; `nolog` hides it |
| QREG-09 | `cqreg` | `log`/`nolog` are accepted but do nothing; there is no iteration log even with `log`. | Fixed: qreg's final log line (iterations, minimized sum); `nolog`/`set iterlog off` hide it |
| IV-17 | `civreghdfe` | `savefprefix()`/`saverfprefix()` without `savefirst`/`saverf` are ignored. In ivreg2 the prefix implies saving. | Fixed |
| IV-19 | `civreghdfe` | `first`, `rf`, `savefirst` and `saverf` regressions under `kernel()`/`dkraay`/`kiefer` use robust or unadjusted VCE instead of the main model's HAC VCE. | Fixed: auxiliary fits are `civreghdfe` fits with the main VCE; stored equations match ivreg2 to 1e-7 |
| IV-20 | `civreghdfe` | First-stage and reduced-form regressions ignore `noconstant` and `partial()`: a constant is added and the partialled variables appear. | Fixed (same change) |
| IV-33 | `civreghdfe` | Stock-Yogo critical values: with more than 3 excluded instruments the 1-instrument row is printed (wrong values). LIML and Fuller get the 2SLS table, and there are no relative-bias rows. | Fixed: ivreg2's full table (from its installed Mata library) by estimator; otherwise the built-in cells or `<not available>` |
| IV-32 | `civreghdfe` | The footer always prints `Partialled-out: _cons`, even with `partial(x)` or nothing partialled. Its uncentered lines show centered SS/R2. | Fixed: ivreg2's footer rules; `e(yy)`, `e(yyc)`, `e(r2u)`, `e(r2c)` match |
| IV-36 | `civreghdfe` | The `rf` table hard-codes 95% CIs (ignores `level()`, `eform`, display options) and uses the IV model's df. `rf nofooter` shows nothing. | Fixed: shown ahead of the main table via `ereturn display`, as ivreg2 |
| BIN-05 | `cbinscatter` | `replace` is accepted and ignored. `savedata()` always overwrites, where binscatter refuses without `replace`. | Fixed: r(602) unless `replace` |
| IO-08 | `cimport excel` | `locale()` is parsed and never used. It has no effect on Windows Stata but does on macOS/Linux. | Not a gap: Stata ignores it on Windows; ctools decodes Unicode directly |
| IO-09 | `cimport excel` | `detail` is parsed and never used. | Not a gap: `detail` only adds `r(N)`/`r(k)`, which `cimport` always returns |
| IO-17 | `cexport excel` | `locale()` is parsed and never used. | Not a gap: ctools writes Unicode (UTF-8 xlsx, UTF-16 xls) |
| HDFE-17 | `creghdfe` | `verbose(-1)` is accepted, but the singleton and convergence lines are still printed. | Fixed: notes, singleton and convergence lines follow reghdfe's verbosity rules |
| HDFE-27 | `creghdfe` | Multi-way clustering where an absorbed effect is nested in one cluster dimension and another has 2 clusters (auto: `price weight length, absorb(turn) vce(cluster turn foreign)`): e(V) is all zero where reghdfe's is nonzero. Found 9 Nov 2026. | **Fixed** (9 Nov): the plugin's e(V) equalled reghdfe's to 1e-13. reghdfe_fix_psd() clips the full matrix, then replaces its coefficient block with that block clipped alone, which can leave V indefinite (here a scaled eigenvalue of -0.12). Stata's `ereturn post` sets such a V to zero or keeps it depending on rounding, which is also why reghdfe gives zero with `vce(cluster foreign turn)`. creghdfe still posts reghdfe's matrix when Stata keeps it; when Stata would zero it, creghdfe posts a positive semidefinite version with the same coefficient block (reghdfe's e(V) with `noconstant`), the constant's covariances reduced to their part in that block's column space and its variance raised to the least PSD value (`src/creghdfe/creghdfe_vce.c`, repair_constant; also for dkraay). Tested in s2_hdfe.do "HDFE-27" |
| IV-31 | `civreghdfe` | `nooutput` hides only the coefficient table and first-stage tables. The header, ffirst table, diagnostics and DoF table still print. | Fixed: nothing is shown; results are still stored |

---

## 2. Gaps by command

### creghdfe → reghdfe 6.13.1

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| HDFE-01 | S2 | Multi-way clustering `vce(cluster a b ...)` | r198 (documented) | **Fixed** | `build/creghdfe.ado:182-185` |
| HDFE-02 | S2 | Heterogeneous slopes in `absorb()`: `id#c.x`, `id##c.x`, `id##c.(x z)` | r198 (documented) | **Fixed** | `build/_ctools_hdfe_absorb.ado:21-24` |
| HDFE-03 | S2 | `newvar=absvar` in `absorb()` to save one FE under a chosen name (and `<name>Slope#`) | r198 (documented) | **Fixed** | `_ctools_hdfe_absorb.ado:21` |
| HDFE-04 | S2 | `vce(dkraay [#])` (Driscoll-Kraay) and its `e(dkraay_bw)`, title4 | r9 "unrecognized vce()" | **Fixed** | `_ctools_vce.ado:15-18` |
| HDFE-05 | S2 | Group and individual FEs: `group()`, `individual()`, `aggregation(mean\|sum)` | r198. See HDFE-10 for `group()` | **Fixed** (9 Nov 2026): LSMR projection (`src/hdfe/hdfe_indiv.c`), reghdfe's singleton and DoF rules; matches reghdfe without weights; with aw/pw the weighted fit, which reghdfe 6.13.1 misses. Refused with `individual()`: fw, slopes, saved effects, bootstraps, `groupvar()` | `validation/stata/s2_hdfe.do` |
| HDFE-06 | S3 | `absorb()` suboptions `generate` (synonym for savefe), `replace`, `noisily` | r198; only `savefe` | **Fixed** (9 Nov): generate/g/savefe, replace, noisily; other suboptions r(198) | s2_hdfe.do "HDFE-06/07" |
| HDFE-07 | S3 | Re-running with `savefe`: reghdfe drops `__hdfe*` first, so reruns and `bootstrap:` work | r110 on the second run | **Fixed** (9 Nov): every run drops `__hdfe*` first, as reghdfe | s2_hdfe.do "HDFE-06/07" |
| HDFE-08 | S3 | Bare `resid` overwrites `_reghdfe_resid` (documented reghdfe behavior) | r110 if it exists; breaks reruns and `bootstrap:` | **Fixed** (9 Nov): bare residuals replace `_reghdfe_resid` | s2_hdfe.do "HDFE-08/09" |
| HDFE-09 | S3 | `residuals` abbreviations (`res`, `resi`) and bare `residuals` | `res` is r198; bare `residuals` probably errors | **Fixed** (9 Nov; Run: reghdfe takes res/resi/resid/residuals, bare or named, and refuses resid with residuals()) | s2_hdfe.do "HDFE-08/09" |
| HDFE-10 | **S1** | `group(var)`: see the priority list | Wrong meaning | **Fixed** | `creghdfe.ado:62` |
| HDFE-11 | S3 | vce abbreviations `vce(rob)`, `vce(un)`, other `r*`/`un*` prefixes | r9 | **Fixed** (9 Nov): reghdfe's case-sensitive vce() rules; `vce(conv)` and `vce(Robust)` are now r(9), as reghdfe | s2_hdfe.do |
| HDFE-12 | S3 | `nosample` | r198 | **Fixed** (9 Nov) | s2_hdfe.do "HDFE-12/13/14/21" |
| HDFE-13 | S3 | `timeit` (timer table, `r(timer)`) | r198 | **Fixed** (9 Nov): reghdfe's timer tables and r(timer) | s2_hdfe.do "HDFE-12/13/14/21" |
| HDFE-14 | S3 | `version(3\|5)`, `nowarn` | **Fixed**: nowarn. version(3|5) won't do: historical reghdfe; r(198) (reghdfe 6.13.1's own version(3) fails with r(3499)) | s2_hdfe.do "HDFE-12/13/14/21" | `creghdfe.ado:59-69` |
| HDFE-15 | S3 | Solver options `technique()`, `transform()`, `acceleration()`, `preconditioner()`, `prune`, `poolsize()`, `compact`, `parallel()`, `fastregress` | r198 (documented) | **Fixed** (9 Nov): transform() and acceleration() run reghdfe's MAP variants; technique() and preconditioner() checked as reghdfe; prune, poolsize(), compact, parallel(), fastregress accepted with no effect | s2_hdfe.do "HDFE-15 solver options" |
| HDFE-16 | S3 | Undocumented solver properties (`abort(0)` to continue past maxiter with `e(converged)=0`, `miniter()`, `accel_start()`...) | r198; nonconvergence always exits r430 | **Fixed** (9 Nov): abort(0) as reghdfe's MAP; miniter(), min_ok(), accel_start(), accel_freq() accepted with no effect | s2_hdfe.do |
| HDFE-17 | **S1** | `verbose(-1)` quiet mode | Accepted; lines still print | **Fixed** | `creghdfe.ado:667-670` |
| HDFE-18 | S3 | Programmer options `nopartialout`, `noregress`, `keepmata`, `varlist_is_touse` (other packages rely on them) | r198 | Won't do: reghdfe's Mata HDFE object; r(198) with a message | s2_hdfe.do "HDFE-18 programming options refused" |
| HDFE-19 | S3 | Factor-variable depvar (`creghdfe 1.foreign x, a(id)`) | `tsrevar` probably rejects it | **Fixed** (9 Nov; Run: reghdfe accepts `1.foreign` and refuses `i.foreign` with r(9)) | s2_hdfe.do |
| HDFE-20 | S3 | `ib#.var` in `absorb()` | Probably rejected; only `i.` is stripped | **Fixed** (9 Nov; Run: reghdfe accepts `ib2.`/`ibn.` terms in absorb()) | s2_hdfe.do |
| HDFE-21 | S3 | Display abbreviation `nohe`/`nohea` | r198 (`noHEADer`) | **Fixed** (9 Nov) | s2_hdfe.do "HDFE-12/13/14/21" |
| HDFE-22 | S4 | Header: "Number of clusters (var)" line and "Statistics robust to heteroskedasticity" subtitle | Not shown | **Fixed** (9 Nov): reghdfe_header's layout and reghdfe's order of notes; logged output identical for 8 models | s2_hdfe.do |
| HDFE-23 | S4 | Footnote "?" marker for inexact redundant counts; `e(dof_table)` Inexact column (always missing) and column names | Never computed | **Fixed** (9 Nov): reghdfe's e(dof_table) row and column names and footnote table | s2_hdfe.do |
| HDFE-24 | S4 | `e(estat_cmd)` not set, so `estat ic` uses `e(rank)` instead of `e(df_m)+1` and AIC/BIC differ | Different AIC/BIC | **Fixed** (9 Nov; Run: AIC/BIC differed, df 2 vs 3): `_creghdfe_estat` uses e(df_m)+1 | `build/_creghdfe_estat.ado`; s2_hdfe.do |
| HDFE-25 | S4 | `e(footnote)` (used by bootstrap/jackknife displays), `e(N_full)`, `e(title3)`, `e(title4)` | Not posted | **Fixed** (9 Nov): e(footnote), e(N_full), reghdfe's e(indepvars) | `build/_creghdfe_footnote.ado`; s2_hdfe.do |
| HDFE-26 | S4 | `predict` variable label for `d`/`xbd` under `noabsorb` | Empty label | **Fixed** (9 Nov) | `build/creghdfe_p.ado` |

### chdfe → hdfe 3.x

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| CHDFE-01 | S3 | `g()` abbreviation of `generate()` | r198 | **Fixed** (9 Nov) | chdfe.do |
| CHDFE-02 | S3 | `tra()` abbreviation of `transform()` | r198 | **Fixed** (9 Nov) | chdfe.do |
| CHDFE-03 | S3 | Experimental `transform(random_kaczmarz)` and `acceleration(hybrid)` | r198 | **Fixed** (9 Nov): random_kaczmarz and hybrid, 7 significant figures against hdfe (fixed-seed order) | chdfe.do |
| CHDFE-04 | S3 | `transform(kaczmarz)` with `acceleration(cg)` (hdfe warns and runs) | r198 | **Fixed** (9 Nov): Kaczmarz with CG runs with hdfe's warning | chdfe.do |
| CHDFE-05 | S3 | Slopes with time-series operators (`id#c.L.x`) | `unab L.x` probably errors | Not a gap: hdfe 3.2.9 also gives r(101) | `chdfe.ado:117-125` |
| CHDFE-06 | S4 | Variables kept by `clear` mode (absorbed slope variables etc.) may differ from hdfe | Unverified | **Fixed** (9 Nov; Run: hdfe keeps neither absorbed nor cluster variables) | chdfe.do |

### civreghdfe → ivreghdfe / ivreg2

**Options and syntax**

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| IV-01 | S2 | Replay: bare `civreghdfe`, `estimates replay`, and replay options (`first ffirst rf sfirst drop*`) | Errors (varlist required) | **Fixed** | `build/civreghdfe.ado:3-26` (no replay branch) |
| IV-02 | S2 | Heterogeneous slopes in `absorb()` | r198 (documented) | **Fixed** (9 Nov 2026): matches ivreghdfe, saved slopes included; not with `vce(bootstrap)`, `fevce()`, `wildbootstrap()` | `validation/stata/s2_iv.do` |
| IV-03 | S2 | Saving FEs: `absorb(fe=id)`, `absorb(..., savefe)` | r198 | **Fixed** (bootstrap session, Phase 4a) | `civreghdfe.ado:299-302` |
| IV-04 | S2 | `center` (centered moments in S) | r198 (documented) | **Fixed** | `civreghdfe.ado:181-184` |
| IV-05 | S2 | `sfirst`, `savesfirst`, `savesfprefix()` | r198 (documented for `sfirst`) | **Fixed** (stored system omits base-level columns) | `civreghdfe.ado:188-191` |
| IV-06 | S2 | `bw(auto)` (automatic Newey-West bandwidth, `e(bwchoice)`) | r198 (documented) | **Fixed** (refused with LIML, k-class, b0()) | `civreghdfe.ado:45` |
| IV-07 | S2 | Kernels Tukey-Hamming, Daniell, Tent | r198 (documented) | **Fixed** | `civreghdfe.ado:527-552` |
| IV-08 | S2 | `smatrix()` / `wmatrix()` (user S or W matrix, `gmmw` model) | r198 | **Fixed** | `civreghdfe.ado:26-88` |
| IV-09 | S2 | `first`/`rf`/`savefirst`/`saverf` with two-way clustering | r198 | **Fixed** | `civreghdfe.ado:356-359` |
| IV-10 | S2 | `orthog()` on included exogenous regressors | r198 ("excluded instruments only") | **Fixed** | `civreghdfe.ado:822-830` |
| IV-11 | S3 | `iweight`s (without absorb) | r101 | **Fixed** (9 Nov): as ivreg2 (no absorb(); homoskedastic VCE; e(N) the integer part of the sum; ivreg2 also truncates N in its id statistics with non-integer sums) | s2_iv.do "IV-11" |
| IV-12 | S3 | `psd0`/`psda` | r198 | **Fixed** (9 Nov): S made PSD in the plugin; e(vce) psd0/psda; not with the bootstraps | s2_iv.do "IV-12" |
| IV-13 | S3 | `sw` (Stock-Watson SEs) | r198 | Deferred: ivreg2's sw is a BETA option computed in compiled livreg2.mlib (no source to port or verify against); still r(198) | `civreghdfe.ado:26-88` |
| IV-14 | S3 | `tvar()`/`ivar()` (ivreg2 runs `tsset` itself) | r198; data must already be tsset | **Fixed** (9 Nov; the gap text was wrong: ivreg2 only accepts the tsset variables, r(5) otherwise) | s2_iv.do "IV-14" |
| IV-15 | S3 | `bvclean`, `fvall`, `fvsep`, `useqr` | r198 | **Fixed** (9 Nov): bvclean (implied by partial()); fvall/fvsep/useqr accepted (identical results in ivreg2) | s2_iv.do "IV-15" |
| IV-16 | S3 | `partial(_cons)` / `fwl(_cons)` | Fails varlist parsing | **Fixed** (9 Nov): partial(_cons), fwl(_cons ...), partial(_all); r(198) with absorb()/noconstant as ivreghdfe/ivreg2 | s2_iv.do "IV-16" |
| IV-17 | **S1** | `savefprefix()`/`saverfprefix()` without `savefirst`/`saverf` | Silently ignored | **Fixed** | `civreghdfe.ado:2205-2214` |
| IV-18 | S3 | reghdfe pass-through options: `groupvar()`, group/individual FE, `nosample`, `technique()`, `accel()`, `transform()`, `prune`, `compact`, `poolsize()`, `parallel()`, `timeit`, `verbose(#)`, `nowarn` | r198 | **Fixed** (9 Nov): groupvar(), solver options (validated, unused), prune, poolsize(), nosample, timeit, nowarn, verbose(#); group()/individual(), compact, parallel() not gaps (ivreghdfe fails: r(2001), r(3500), r(3200)) | s2_iv.do "IV-18" |
| IV-19 | **S1** | First-stage/RF VCE under kernel, dkraay or kiefer | Not HAC | **Fixed** | `civreghdfe.ado:1566,1582-1583,1925-1926,2225-2226` |
| IV-20 | **S1** | First-stage/RF ignore `noconstant` and `partial()` | Wrong auxiliary models | **Fixed** | `civreghdfe.ado:1582,1925,2225` |
| IV-21 | S3 | Abbreviations `sm` (small), `nohe` (noheader), `nofo` (nofooter) | r198 | **Fixed** (9 Nov) | s2_iv.do "IV-21" |
| IV-22 | S3 | `residuals2` variable name: ivreghdfe uses `_reghdfe_resid` | Uses `_civreghdfe_resid` | **Fixed** (9 Nov) | s2_iv.do "IV-22 residuals2" |

**Stored results**

| ID | Sev | Gap | Status | Evidence |
|---|---|---|---|---|
| IV-23 | S4 | `e(S)`, `e(W)` | **Fixed** (9 Nov) | s2_iv.do "IV-23/24/28/29/30 stored results" |
| IV-24 | S4 | `e(vce)` for the main model | **Fixed** (9 Nov) | s2_iv.do "IV-24 e(vce)" |
| IV-25 | S2 | Weak-instrument statistics: `e(first)` matrix (Shea partial R2, SW and AP F/chi2), Anderson-Rubin `e(arf ...)`, Stock-Wright `e(sstat ...)`. Not computed anywhere. | **Fixed** | `civreghdfe.ado:1747-1786` |
| IV-26 | S4 | LIML Anderson-Rubin overid test `e(arubin*)` | **Fixed** (9 Nov; also the AR footer line) | s2_iv.do "IV-23/24/28/29/30 stored results liml" |
| IV-27 | S4 | ivreg2 result names: `e(sargandf)`, `e(cstatp/cstatdf/clist)`, `e(estat*)`/`e(elist)`, `e(redstat/redp/reddf/redlist)`, `e(rkf)`, `e(cd)`, `e(ccev)`, `e(cdev)`. ctools posts its own names. | **Fixed** (9 Nov) | s2_iv.do "IV-27 test names" |
| IV-28 | S4 | `e(condzz)`, `e(condxx)`, `e(center)`, `e(partialsmall)`, `e(bwchoice)`; `e(dofminus)` only when nonzero. (`e(r2u)`, `e(yy)`, `e(yyc)` and `e(r2c)` under noconstant were fixed with IV-32.) | **Fixed** (9 Nov) | s2_iv.do "IV-23/24/28/29/30 stored results" |
| IV-29 | S4 | 0/1 varlists `e(depvar0/1)`, `e(instd0/1)`, `e(inexog0/1)`, `e(exexog0/1)`, `e(insts0/1)`, `e(partial1)` | **Fixed** (9 Nov) | s2_iv.do "IV-23/24/28/29/30 stored results" |
| IV-30 | S4 | `e(model)` values (`iv`, `ols`, `liml`, `gmmw`). ctools posts `2sls`, `fuller`, and `2sls` for OLS. | **Fixed** (9 Nov) | s2_iv.do "IV-23/24/28/29/30 stored results" |
| IV-30b | S4 | `e(title)` ignores `title()`; `e(subtitle)` and `e(hacsubtitle*)` missing | **Fixed** (9 Nov) | s2_iv.do "IV-23/24/28/29/30 stored results" |
| IV-30c | S4 | `e(firsteqs)`, `e(rfeq)`, `e(sfirsteq)`; default prefix `_ivreg2_` (ctools uses `_civreghdfe_`) | **Fixed** (9 Nov): _ivreg2_ prefixes; the extra <prefix>main stores removed | s2_iv.do "IV-30c stored names" |
| IV-30d | S4 | Saved first-stage/RF equations lack `e(sample)`, `e(S)`, `e(k_eq)` and kernel info | **Fixed** (9 Nov): e(S), e(k_eq) | s2_iv.do "IV-30d stored equation" |
| IV-30e | S4 | reghdfe footnote results: `e(absvars)`, `e(extended_absvars)`, `e(dof_table)`, `e(df_a_initial/redundant/nested)`, `e(dofmethod)`, `e(N_hdfe_extended)`, `e(drop_singletons)`, `e(title2)` | **Fixed** (9 Nov) | s2_iv.do "IV-30e absorb stored results" |

**Output and postestimation**

| ID | Sev | Gap | Status | Evidence |
|---|---|---|---|---|
| IV-31 | **S1** | `nooutput` is only partial | **Fixed** | `civreghdfe.ado:1573,1636,1725,1747,1791` |
| IV-32 | **S1** | Footer partialled-out list and uncentered SS/R2 are wrong | **Fixed** | `civreghdfe.ado:1715-1720,2150-2156` |
| IV-33 | **S1** | Stock-Yogo critical value table incomplete or wrong | **Fixed** | `civreghdfe.ado:1840-1891` |
| IV-34 | S4 | Header title ignores the estimator: always "IV (2SLS)" instead of LIML, GMM2s, CUE, k-class or OLS | **Fixed** | `civreghdfe.ado:1638-1643` |
| IV-35 | S4 | ffirst/first summary lacks SW/AP F, Shea R2, under-identification and weak-ID tests, AR/SW inference (depends on IV-25) | **Fixed** | `civreghdfe.ado:1747-1786` |
| IV-36 | **S1** | `rf` display problems (fixed 95% CI, df, `nofooter`) | **Fixed** | `civreghdfe.ado:1791-1792,1898-1973` |
| IV-37 | S4 | No restricted-J line with `orthog()` | **Fixed** (9 Nov): footer as ivreghdfe (restricted J, -orthog-/-redundant-/-endog- blocks, LIML AR) | s2_iv.do "IV-37 footer" |
| IV-38 | S3 | `predict, residuals` without absorb is limited to `e(sample)`; ivreg2 predicts out of sample | **Fixed** (9 Nov) | s2_iv.do "IV-38 out-of-sample residuals" |
| IV-39 | S3 | Replay and e(sample)-based commands on restored first-stage/RF equations | **Fixed** (9 Nov) | s2_iv.do "IV-39 replay and e(sample)" |

### cppmlhdfe → ppmlhdfe 2.3.3

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| PPML-01 | **S1** | `verbose(-1)` turns output on | Wrong | **Fixed** | `src/cppmlhdfe/cppmlhdfe_input.c:44-46`; `cppmlhdfe_irls.c:509,768` |
| PPML-02 | **S1** | No IRLS iteration log; `log`/`nolog` are no-ops | Ignored | **Fixed** | `build/cppmlhdfe.ado:49,890-895` |
| PPML-03 | S3 | Replay shows only the coefficient table (no header or footnote); `cppmlhdfe, eform` and `, noheader` fail | r198 | **Fixed** (9 Nov): header, table, footnote; eform/irr, noheader, notable, display options; r(table) | `cppmlhdfe.ado` (_cppmlhdfe_replay); cppmlhdfe.do "PPML-03" |
| PPML-04 | S2 | Alternative slope syntax `absorb(id##(c.x c.z))` (documented in ppmlhdfe help) | r101 | **Fixed** | `cppmlhdfe.ado:205-235` |
| PPML-05 | S3 | `absorb(..., maxiterations() verbose())` suboptions | r198 | **Fixed** (9 Nov): maxiterations() is iterate(); verbose() is the HDFE verbosity (-1 hides the singleton note) | `_cppmlhdfe_abs_options.ado`; cppmlhdfe.do "PPML-05/06" |
| PPML-06 | S3 | reghdfe pass-through: `prune` (named in ppmlhdfe help; `e(pruned)`), `condition` (`e(finite_condition)`), `rre()`, `algorithm()`, `slopemethod()`, `precondition` | r198 | **Fixed** except e(pruned)/e(num_pruned): condition (two effects) posts e(finite_condition)=-1 as ppmlhdfe; rre()/algorithm()/slopemethod()/prune accepted, estimates unchanged; precondition refused. e(pruned) deferred: ppmlhdfe's prune needs moremata (r(3499) here) | `cppmlhdfe.ado`; cppmlhdfe.do "PPML-06" |
| PPML-07 | S3 | `acceleration()`/`transform()` values: case-insensitive prefixes and `accel(off)` | r198 | **Fixed** (9 Nov): any case, leading parts, cg/sd/off, as reghdfe5_parse | `cppmlhdfe.ado` |
| PPML-08 | S3 | `r2` with `tagsep` (regress on the separation certificate) | r198 (documented) | **Fixed** (9 Nov): the certificate is regressed with creghdfe (ppmlhdfe: reghdfe); relu_report_r2(1) synonym | cppmlhdfe.do "PPML-08" |
| PPML-09 | S3 | Undocumented GLM properties (`target_inner_tol()`, `relu_debug()`...) | r198 | **Fixed** (9 Nov): target_inner_tol(), relu_debug(), relu_report_r2(); every property also inside absorb() by exact name | `_cppmlhdfe_abs_options.ado`; cppmlhdfe.do "PPML-09" |
| PPML-10 | S3 | `ppmlhdfe, version` / `reload` | Error | **Fixed** (9 Nov): version shows which and ctools, version; reload prints how to reinstall (ppmlhdfe's reinstalls from its author's local path) | `cppmlhdfe.ado` |
| PPML-11 | S4 | Header lacks Residual df, Deviance and per-dimension cluster counts. Log pseudolikelihood overflows its column (`-2.760e+0` wraps). | Stata | **Fixed** (9 Nov): ppmlhdfe_header layout |
| PPML-12 | S4 | Footnote is a raw `matrix list e(dof_table)`, not reghdfe's formatted table; prints even with noabsorb | Stata | **Fixed** (9 Nov): reghdfe's table, none without absorb() |
| PPML-13 | S4 | `e(title3)` not posted with cluster VCE | Stata | **Fixed** (9 Nov) |
| PPML-14 | S4 | `e(N_hdfe)`/`e(title2)` count unique categorical keys, not absvar terms | Stata | **Fixed** (9 Nov): absorb() terms |
| PPML-15 | S4 | `e(indepvars)` omits `_cons` | Stata | **Fixed** (9 Nov): ppmlhdfe's bn list and "  _cons" |
| PPML-16 | S4 | `predict`: no variable labels and no "(option mu assumed...)" note | Code | **Fixed** (9 Nov) |
| PPML-17 | S4 | Missing warnings: keepsingletons/separation, low standardized depvar, "IRLS step detected N separated observations" | Code | **Fixed** (9 Nov): all of ppmlhdfe's notes and warnings, in its order |

### cqreg → qreg (StataNow 19.5)

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| QREG-01 | S2 | Weights `fweight`, `iweight`, `pweight` (pweight implies robust); `e(wtype)`, `e(wexp)`, `e(sum_w)` | r101; `e(sum_w)` always N | **Fixed**: fw/iw/pw as qreg; negative iweights refused (r(402)) | `build/cqreg.ado:30` |
| QREG-02 | S2 | `absorb()` | Rejected (documented) | Not a gap: qreg has no `absorb()` | `docs/commands/cqreg.md:27` |
| QREG-03 | S2 | `by:` prefix (qreg is `byable(onecall)`) | r190 | **Fixed** | `cqreg.ado:3` |
| QREG-04 | S2 | `mi estimate`, `stepwise`, `nestreg` (qreg declares `prop(sw mi)`) | r198 / r199 | **Fixed**: `prop(sw mi)` | `cqreg.ado:3` |
| QREG-05 | S2 | `bayes:` prefix | r198 | Won't fix: `bayes:` supports a fixed list of official commands | — |
| QREG-06 | S3 | `quantile(25)` meaning .25 (percent form) | r498 | **Fixed** | `cqreg.ado:53-56` |
| QREG-07 | S3 | qreg's `vce(type, density bandwidth)` suboption syntax, e.g. `vce(robust, kernel(gaussian) bofinger)`, `vce(iid, residual)` | r198; only cqreg's own `denmethod()`/`bwmethod()` | **Fixed** (9 Nov): _cqreg_vceparse (also excluster, vce(, residual)); must agree with denmethod()/bwmethod() | `cqreg.ado`; s2_qreg.do |
| QREG-08 | S3 | Abbreviations `vce(rob)`, `bo`/`hs`/`ch`, `fit`, `res`, `residuals` | r198 | **Fixed** (9 Nov): rob/r, bo/hs/ch, fit/res/residuals/ker() in vce(), denmethod(), bwmethod() | s2_qreg.do "QREG-08" |
| QREG-09 | **S1** | No iteration log; `log`/`nolog` are no-ops | Ignored | **Fixed** | `cqreg.ado:32-35` |
| QREG-10 | S3 | `wlsiter(#)` | r198 (documented) | **Fixed** (9 Nov): accepted and checked as qreg (r(198) for # < 1); no effect on the interior-point solver | s2_qreg.do "QREG-10/11" |
| QREG-11 | S3 | `iterate(#)`, `trace` | r198 (has `maxiter()` instead) | **Fixed** (9 Nov): iterate() checked as qreg, no effect (the limit is maxiter()); trace shows the solution's coefficients | s2_qreg.do "QREG-10/11" |
| QREG-12 | S2 | Bootstrap VCE (qreg via `bootstrap:`/`bsqreg`) | Documented as unsupported, but `bootstrap: cqreg` worked in testing. Reconcile the docs. | Not a gap: `bootstrap: cqreg` works; qreg has no `vce(bootstrap)` (docs corrected) | `docs/COMPATIBILITY.md` |
| QREG-13 | S4 | Hidden `e(basis)`, `e(tolerance)`; legacy `$S_E_*` globals | Not set | **Fixed** (9 Nov): e(basis) (qreg's set at a unique optimum), e(tolerance) = tolerance(), $S_E_* ($S_E_cmd cqreg) | `cqreg.ado`; `src/cqreg/cqreg_regress.c`; s2_qreg.do "QREG-13" |
| QREG-14 | S4 | `e(cmd)` is `cqreg`, so code branching on `e(cmd)=="qreg"` misses it. predict, margins, estat, linktest and bootstrap all work. | — | Not a gap (by design): mfp rejects cqreg by command name whatever e(cmd) (r(198)); bayesboot: cqreg works; replay needs e(cmd) cqreg | — |
| QREG-15 | S3 | `version 12:` behavior (aweights, old SE formula) | None | Deferred: qreg's version-12 SEs use _qreg's built-in, undocumented density f_r (it matches no current method); cqreg notes the difference under version 12 or earlier | — |

### cimport / cexport → import / export (StataNow 19.5)

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| IO-01 | S2 | `import parquet` (namelist, `describe`, `rowrange()`, `favormemory`) | r198 "unknown subcommand" (documented) | **Fixed** natively; Brotli too (9 Nov 2026, `src/io/parquet_brotli.c`) | `build/cimport.ado:35-38` |
| IO-02 | S2 | `import fred` (all options) | r198 (documented) | **Fixed**: runs Stata's import fred | `cimport.ado:19-22` |
| IO-03 | S2 | `import haver` / `haverdirect` (Windows) | r198 (documented) | **Fixed**: runs Stata's import haver | `cimport.ado:12-22` |
| IO-04 | S2 | `spshape2dta` (shp + dbf to linked .dta pair) | No equivalent; `cimport shp`/`dbase` cover only the raw pieces | Not a gap: spshape2dta is not a command ctools replaces | — |
| IO-05 | S2 | URL file names (`import delimited using https://...`; check other formats) | r601 from `confirm file` | **Fixed** for delimited; Stata's excel/sas import take no URLs either (r(603)) | `cimport.ado:104,776` |
| IO-06 | S3 | delimited: `decimals()`, `groups()` abbreviations | r198 | **Fixed**: DECIMALSeparator/GROUPSeparator (audit B8) | cimport.do "decimals() and groups() abbreviations" |
| IO-07 | S3 | delimited/sas/spss: `encoding()` names Java accepts but the Windows iconv table lacks (KOI8-U, Big5-HKSCS, EUC-TW, x-MacCyrillic, `latin1`, SAS MAC* encodings...) | r198 or iconv failure | **Fixed** (9 Nov): Java charset names and aliases; single-byte charsets decoded from Java's tables on every platform; malformed UTF-8 bytes dropped. On Windows only, Big5-HKSCS, x-EUC-TW, x-ISCII91, ISO-2022-JP-2, CESU-8 and the IBM EBCDIC CJK pages still give r(198) (no Windows code page), and invalid multi-byte sequences become replacement characters | `src/cimport/cimport_java_charsets.inc`; cimport.do "IO-07" |
| IO-08 | **S1** | excel import `locale()` | Ignored (help says "accepted for compatibility") | Not a gap | `cimport.ado:737,741,780-802` |
| IO-09 | **S1** | excel import `detail` | Ignored | Not a gap | `cimport.ado:737,741` |
| IO-10 | S3 | excel import `sh()`, `all`, `all(fmt)` abbreviations | r198 | **Fixed** (9 Nov): all/allstri; all(fmt) read as neither, as native | cimport.do "IO-10" |
| IO-11 | S3 | `.xlsm`/`.xltx` files: import rejects them; export renames `file.xlsm` to `file.xlsm.xlsx` | Different | **Fixed** (9 Nov): import picks the kind by content (.xltx, .xltm, any name); export keeps any name, adds .xlsx only without an extension | cimport.do and cexport.do "IO-11" |
| IO-12 | S3 | sas/spss namelist and `if`/`in` applied after a full load, so files that only fit when subset fail. `keep` is case-sensitive. | Partial | **Fixed** (9 Nov): the name list and in applied while reading; native name rules (exact, a-b is a and b, unknown names ignored, r(111) without a match) | `src/io/cio_readstat.c`; cimport.do "IO-12" |
| IO-13 | S4 | sasxport5 `describe, member()`: `r(n_members)`/`r(members)` not set; `r(N)`/`r(k)` come from the last member | Partial | Not a gap: r(N)/r(k)/r(size)/r(n_members)/r(members) equal native's with and without member() | cimport.do "IO-13" |
| IO-14 | S3 | `clear` semantics: ctools refuses (r4) whenever data exist; Stata may allow import over unchanged data (`c(changed)`) | Different | Not a gap: native refuses with r(4) with unchanged data for delimited/excel/spss/sasxport8/dbase; only sasxport5 checks c(changed), as cimport does | cimport.do "IO-14" |
| IO-15 | S3 | delimited export `q` (quote), `novar` (novarnames) abbreviations | r198 | **Fixed** (9 Nov) | cexport.do "IO-15" |
| IO-16 | S3 | excel export `sh()`, `date()`, `miss()` abbreviations | r198 | **Fixed** (9 Nov) | cexport.do "IO-16" |
| IO-17 | **S1** | excel export `locale()` | Ignored | Not a gap | `cexport.ado:438` |
| IO-18 | S3 | sasxport5 export `val()` abbreviation of `vallabfile()` | r198 | **Fixed** (9 Nov) | cexport.do "IO-18" |
| IO-19 | S4 | delimited export line endings: ctools always writes LF unless `crlf`; Stata on Windows may write CRLF | Unverified | Not a gap (by design): native writes CRLF on Windows; cexport writes LF unless crlf (documented) | cexport baseline byte checks |

### cmerge → merge

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| MERGE-01 | S3 | Old (pre-v11) syntax `merge [varlist] using f1 f2 ...` including multiple using files, `unique`, `uniqmaster` | r111 "invalid merge type" | **Fixed** (9 Nov): runs Stata's merge (merge_10) | cmerge.do "MERGE-01" |
| MERGE-02 | S3 | Invalid match type error code (merge r198, cmerge r111) | Different rc | **Fixed** (9 Nov): r(198), merge's message | cmerge.do "MERGE-02" |
| MERGE-03 | S3 | `debug` (undocumented) | r198 | **Fixed** (9 Nov): accepted; shows the plugin commands | cmerge.do "MERGE-03 debug" |
| MERGE-04 | S2 | `force` with string/numeric key mismatch | r106 (documented); confirm merge's own behavior | Not a gap: `merge, force` also exits r(106) | `cmerge.ado:653-662` |
| MERGE-05 | S4 | `force` warning for string/numeric mismatch on non-key variables | Silent | **Fixed** (9 Nov) | cmerge.do "MERGE-05/06" |
| MERGE-06 | S4 | "(variable x was int, now float ...)" notes when storage types are widened | Silent | **Fixed** (9 Nov): merge's notes, in the using data's variable order | cmerge.do "MERGE-05/06" |
| MERGE-07 | S3 | `cmerge 1:1 _n using f` with a variable-less master | Probably r111 from `unab _all` | **Fixed** (9 Nov): r(3) "no dataset in use", as merge | cmerge.do "MERGE-07" |
| MERGE-08 | S2 | strL output (updates or using-only strL values) | Unsupported (documented) | **Fixed**: native merge path | `docs/commands/cmerge.md:161` |

### csort → sort

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| SORT-01 | S3 | `sort varlist in range` | r101. The help wrongly says native sort rejects `in`. | **Fixed** (9 Nov): in through sort, stable | csort.do "SORT-01" |
| SORT-02 | S3 | `stable` option (csort is always stable but rejects the word) | r198 | **Fixed** (9 Nov) | csort.do "SORT-02" |
| SORT-03 | S2 | Any strL variable in the dataset, even a non-key one | r109 (documented) | **Fixed**: `sort, stable` with strL data | `csort.ado:81-88` |
| SORT-04 | S4 | `r()` preserved: sort leaves the caller's `r()` alone; csort is rclass and overwrites it | Different | Deferred: csort's r() results are documented and read by contracts.do, s2_strl.do, csort.do and comparisons.do; needs a design decision | `docs/commands/csort.md` |
| SORT-05 | S3 | Frame alias variables (Stata 18+) | Untested | Not a gap: alias variables as key or non-key give sort, stable's order, and the alias stays an alias | csort.do "SORT-05" |

### cencode / cdecode → encode / decode

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| ENC-01 | S4 | `r()` preserved (encode posts nothing; cencode overwrites `r()`) | Different | **Fixed** (9 Nov): the caller's r() is held and returned, plus r(N_unique) and r(N_vars), also on failure | `validation/stata/cencode.do` "ENC-01" |
| ENC-02 | S2 | More than 65,536 unique values: check whether Stata 19 `encode` allows more | r134 | Not a gap: `encode` also exits r(134) | `src/cencode/cencode_impl.c:283-390` |
| DEC-01 | S4 | Non-numeric source: decode prints "not numeric"; cdecode exits r108 silently | No message | **Fixed** (9 Nov): "not possible with string variable" r(108); "x not labeled" r(182) | `validation/stata/cdecode.do` "DEC-01" |
| DEC-02 | S4 | `r()` preserved (cdecode overwrites it) | Different | **Fixed** (9 Nov): as ENC-01, plus r(N_vars) and r(N_plugin) | `validation/stata/cdecode.do` "DEC-02" |

### cdestring → destring

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| DESTR-01 | **S1** | `ignore()` suboptions `illegal`, `asbytes`, `aschars` | Misread as characters to strip | **Fixed** | `build/cdestring.ado:15,164-174`; `src/cdestring/cdestring_impl.c:68-73,226-252` |
| DESTR-02 | S3 | `generate(n1-n3)` variable ranges | r198 (counts words) | **Fixed** (9 Nov): generate() is a newvarlist, with destring's errors | `validation/stata/cdestring.do` "DESTR-04 messages" |
| DESTR-03 | S3 | No varlist means all variables. `generate()` must name every variable, and numeric ones report "already numeric". | Uses string variables only: r198 on a valid count; r111 on an all-numeric dataset | **Fixed** (9 Nov): no varlist was already all variables; numeric ones are now reported "already numeric", and a dataset with no variables gives r(100) | `validation/stata/cdestring.do` "DESTR-03" |
| DESTR-04 | S4 | Per-variable messages: characters removed, "(N missing values generated)", `force` "contains nonnumeric characters" | One combined line | **Fixed** (9 Nov): destring's per-variable messages verbatim | `validation/stata/cdestring.do` "DESTR-04 messages" |
| DESTR-05 | S4 | "x already numeric; no replace" wording | "x: already byte" | **Fixed** (9 Nov) | `validation/stata/cdestring.do` "DESTR-04 messages" |
| DESTR-06 | S2 | strL values over 2,045 bytes or binary | Error (documented) | **Fixed** | `src/ctools_transport.c:37-51` |
| DESTR-07 | S3 | Stata hex number format (`1.8X+004`) | Probably treated as nonnumeric | **Fixed** (9 Nov; confirmed: real() and destring read %21x, e.g. 1.8X+004 = 24) | `src/cdestring/cdestring_impl.c` parse_stata_hex; cdestring.do "DESTR-07" |

### csplit → split

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| SPLIT-01 | S4 | `destring` output (which pieces stayed string, and why) | `quietly destring`; silent | **Fixed** (9 Nov): destring's messages shown | `validation/stata/csplit.do` "SPLIT-01/02" |
| SPLIT-02 | S4 | "variable(s) created as string:" plus `ds` listing | One "variables created:" line | **Fixed** (9 Nov): split's report, the ds listing, r(varlist) with its trailing space | `validation/stata/csplit.do` "SPLIT-01/02" |

### csample / cbsample → sample / bsample

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| SAMP-01 | S3 | `sample # , count` syntax | Needs `count(#)` (documented) | **Fixed** (9 Nov): `#, count` (and `count(#)`), sample's messages, codes and no-op cases | `validation/stata/csample.do` |
| SAMP-02 | S3 | `pduplicates()` (undocumented) | r198 | **Fixed** (9 Nov): accepted, no effect (the ctools draw has no ties) | `validation/stata/csample.do` |
| SAMP-03 | S3 | Data left sorted by `by()` vars with sort marker set | Original order kept; no marker | **Fixed** (9 Nov): sorted by by() (stable) | `validation/stata/csample.do` "SAMP-03" |
| SAMP-04 | S3 | `version 6` algorithm | None | **Fixed** (9 Nov): under version 6 or earlier csample calls sample | `validation/stata/csample.do` |
| SAMP-05 | S4 | "(1 observation deleted)" singular | Always plural | **Fixed** (9 Nov) | `validation/stata/csample.do` |
| SAMP-06 | S2 | strL `by()` values over 2,045 bytes or binary | Error | **Fixed** | `src/ctools_transport.c:37-51` |
| BSAMP-01 | S2 | strL `strata()`/`cluster()` values over 2,045 bytes or binary | Error | **Fixed** | `build/cbsample.ado:232-233` |
| BSAMP-02 | S3 | `version 7` algorithm (`bsample_7`) | None | **Fixed** (9 Nov): under version 7 or earlier cbsample calls bsample_7 | `validation/stata/cbsample.do` |

### cinterpolate → ipolate

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| IPOL-01 | S4 | `r(N)` (missing values generated) left after the command | Not set | **Fixed** (9 Nov): r(N) is the number of missing values generated | `validation/stata/cinterpolate.do` "IPOL-01/02" |
| IPOL-02 | S4 | "(1 missing value generated)" singular | Always plural | **Fixed** (9 Nov) | `validation/stata/cinterpolate.do` "IPOL-01/02" |
| IPOL-03 | S2 | strL `by()` values over 2,045 bytes | Error (documented) | **Fixed** | `docs/USAGE.md:16` |
| IPOL-04 | — | SSC `cipolate` (cubic) edge cases are unverified because the package is not installed | — | Deferred: SSC `cipolate` is not installed (`~/ado/plus/c/cipolate.ado` is a stale ctools 1.2.1 copy) | — |

### cbinscatter → binscatter 7.02 / binsreg

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| BIN-01 | S2 | Several y variables (`binscatter y1 y2 x`) | Rejected (`max=2`) | **Fixed** (one graph; bins match binscatter per y) | `build/cbinscatter.ado:13` |
| BIN-02 | S2 | `medians`, `noaddmean`, `rd()`, `xq()`, `mcolors()`, `lcolors()`, `savegraph()`, `randvar()`, `randcut()`, `randn()` | r198 (documented) | **Fixed** (refused with method(binsreg)) | `cbinscatter.ado:818-837` |
| BIN-03 | S2 | `genxq()` (bin variable) | r198 (documented) | **Fixed** | `cbinscatter.ado:129-132` |
| BIN-04 | S3 | `savedata()` format: binscatter writes `.csv` plus a `.do` to redraw the graph; ctools writes one `.dta` | Different | **Fixed** (9 Nov): without an extension, `savedata()` also writes binscatter's `.csv` and `.do` (default method), checked line by line against binscatter; the `.dta` is kept | s2_bin.do "BIN-04" |
| BIN-05 | **S1** | `replace` ignored; `savedata()` always overwrites | Ignored | **Fixed** | `cbinscatter.ado:815-822,880` |
| BIN-06 | S4 | `e(sample)`, `e(graphcmd)`, `e(y#_coefs)`, `e(byvalues)`, `e(rdintervals)`. ctools clears `e()` and posts `e(coefs)` with a different layout. | Missing | **Fixed** (9 Nov): `e(sample)`, `e(y#_coefs)` (binscatter's layout and names), `e(byvalues)`; `e(graphcmd)`, `e(rdintervals)` and `e(coefs)` kept | s2_bin.do "BIN-06" |
| BIN-07 | S3 | Abbreviations `n()` (nquantiles), `m()`/`ms()` (msymbols) | Fall through to twoway and error | **Fixed** (9 Nov; binscatter accepts `n()`, `m()`, `ms()` with r(0)) | s2_bin.do "BIN-07" |
| BIN-08 | S3 | Legacy `nbins()`, `create_xq`, `x_q()`, `symbols()`, `unique()`, `linetype(noline)` | Error | **Fixed** (9 Nov; binscatter's notes; both forms given: r(198) where binscatter gives r(0)) | s2_bin.do "BIN-08" |
| BIN-09 | S3 | `nquantiles()` above 1000 | r198 | **Fixed** (9 Nov): binscatter itself stops with r(198) above 100 quantile bins (fastxtile) and allows more only for a discrete x; cbinscatter is bounded by `c(max_matdim)` | s2_bin.do "BIN-09" |
| BIN-10 | S2 | binsreg: IMSE-optimal `nbins()` and `binsmethod()`, `nbinsrot()`, `binspos()`, `pselect()`/`sselect()` | Fixed quantile bins only | Won't do now: binsreg is not installed to verify against | `cbinscatter.ado:96-110` |
| BIN-11 | S2 | binsreg: `dots(p s)` with p>0, `line()`, `ci()`, `cb()`, `deriv()`, `at()`, `vce()`, `level()`, `nsims()` and related options (confidence intervals, bands, tests) | Missing (documented) | Won't do now: binsreg is not installed to verify against | — |
| BIN-12 | S3 | binsreg: `polyreg(p)` above 3 | Max cubic | **Documented** (help "Differences from binsreg"); not until binsreg is installed | — |
| BIN-13 | S3 | binsreg: `samebinsby`, `bycolors()`, `bysymbols()`, `bylpatterns()` | Missing; bins per group only | **Documented** (help "Differences from binsreg"); not until binsreg is installed | `cbinscatter.sthlp:108-110` |
| BIN-14 | S3 | binsreg: `masspoints()`, `dfcheck()`, `noplot`, `plotxrange()`/`plotyrange()`, `reghdfeopt()`, binsreg-layout `savedata()` | Missing | **Documented** (help "Differences from binsreg"); not until binsreg is installed | — |
| BIN-15 | S3 | binsreg: positional controls (`binsreg y x w1 w2`) | Needs `controls()` | **Documented** (help "Differences from binsreg"); not until binsreg is installed | — |
| BIN-16 | S4 | binsreg `e()` results | Missing | **Documented** (help "Differences from binsreg"); not until binsreg is installed | — |

### cpsmatch → psmatch2 4.0.12

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| PSM-01 | S2 | Several outcomes `outcome(varlist)` | One variable only | **Fixed** | `build/cpsmatch.ado:14` |
| PSM-02 | S2 | `ate` (ATU/ATE estimates, table rows, `r(atu)`, `r(ate)`, `r(seatu)`, `r(seate)`) | r198 (documented) | **Fixed** | `cpsmatch.ado:13-19` |
| PSM-03 | S2 | `mahalanobis()` and `w()` (covariate-metric matching) | r198 | **Fixed** | — |
| PSM-04 | S2 | `llr` (local linear regression matching) | r198 | **Fixed** | — |
| PSM-05 | S2 | `spline` / `nknots()` | r198 | **Partly fixed**: calls the `spline` package as psmatch2 does; unverified here (package not installed) | — |
| PSM-06 | S2 | `ai()`, `population`, `altvariance` (Abadie-Imbens variance) | r198 (documented) | **Fixed** | — |
| PSM-07 | S2 | `trim()` (density-based common support) | r198 | **Fixed** | — |
| PSM-08 | S3 | `odds`, `index` (match on log odds or linear index) | r198 | **Fixed** | — |
| PSM-09 | S3 | `quietly` (hide the score model output) | r198 | **Fixed** | `cpsmatch.ado:183-186` |
| PSM-10 | S3 | `warnings` (duplicate pscore warning) | r198 | **Fixed** | — |
| PSM-11 | S2 | `noreplacement` with `neighbor(k>1)` | Requires `neighbor(1)` (documented) | **Fixed**: `neighbor(k>1)` with `noreplacement` is refused, as psmatch2 refuses it | — |
| PSM-12 | S2 | Neighbor ID variables `_n1`...`_nK` and psmatch2's `_id`. ctools' `_id` means something else (the obs number of the first match). | Missing or different; breaks downstream code | **Fixed** | `cpsmatch.ado:203-209` |
| PSM-13 | S4 | `_pdif` (distance to nearest match), `_<outcome>` (matched counterfactual) | Not created | **Fixed** | — |
| PSM-14 | S4 | Value labels on `_treated`/`_support`, variable labels, `char _weight[Type]`; `_nn = 0` for on-support controls | Missing / missing value | **Fixed** (labels read `cpsmatch:`) | `cpsmatch.ado:202-206` |
| PSM-15 | S4 | `r(depvar)`, `r(exog)` | `r(outcome)`/`r(treatvar)` instead | **Fixed** | `cpsmatch.ado:453-457` |
| PSM-16 | S3 | `radius` without `caliper()` matches every control in psmatch2 | Uses 0.25 (documented) | **Fixed** (9 Nov): every control, with psmatch2's missing-bandwidth weights and SE | s2_psm.do "PSM-16" |

### cdistinct → distinct 1.2.1

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| DIST-01 | S2 | strL values over 2,045 bytes or binary | Error (documented) | **Fixed** | `src/ctools_transport.c:37-45` |

### crangejoin → rangejoin 1.1.3

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| RJOIN-01 | S2 | strL anywhere in the master or kept using variables, even when not involved in the join | r109 (documented) | **Fixed** (strL `by()` keys refused, as rangejoin's merge refuses them) | `build/crangejoin.ado:113-118,211-214` |
| RJOIN-02 | S3 | Output above 2,147,483,647 observations | Capped (documented) | Deferred: the plugin API's 32-bit observation indices; needs a non-plugin path and a test join of more than 2^31 rows; errors with r(920) | — |

### crangestat → rangestat 1.1.1

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| RSTAT-01 | S2 | Variables as interval bounds `interval(key lowvar highvar)` | r198 (documented) | **Fixed** (C, per-observation windows) | `build/crangestat.ado:75-99` |
| RSTAT-02 | S2 | `cov`, `corr`, `reg` statistics | r3499 (documented) | **Fixed** (Mata, rangestat's names) | `crangestat.ado:150-180` |
| RSTAT-03 | S2 | User-written Mata functions `(myfunc) varlist` | r3499 (documented) | **Fixed** (Mata) | `crangestat.ado:177-180` |
| RSTAT-04 | S3 | `casewise`, `describe`, `local()`, `nocheck` | r198 (documented) | **Fixed** (9 Nov): casewise for every statistic, describe, local(), nocheck | crangestat.do "RSTAT-04" |
| RSTAT-05 | S4 | `r(newvars)`, `r(rangestat_version)` | Not rclass | **Fixed** (9 Nov): r(newvars), r(rangestat_version) | crangestat.do "RSTAT-05" |
| RSTAT-06 | S3 | Variable range with spaces (`x1 - x3`) | r111 | **Fixed** (9 Nov) | crangestat.do "RSTAT-06" |
| RSTAT-07 | S2 | strL `by()` values over 2,045 bytes or binary | Error | **Fixed** | `src/ctools_transport.c:37-45` |

### cwinsor → winsor2 / gstats winsor

| ID | Sev | Gap | ctools behavior | Status | Evidence |
|---|---|---|---|---|---|
| WINS-01 | S3 | `label` (append "-Winsor(p,p)"/"-Trim(...)" to the variable label) | r198 | **Fixed** | `build/cwinsor.ado:74` |
| WINS-02 | S3 | Abbreviations `s()` (suffix), `t` (trim), `p()` (prefix in gstats), `gen()`, `l` | r198, or a different meaning for `p()` (errors, not silent) | **Fixed** (9 Nov; `p(#)` holding a number stays the lower percentile; `prefix()` with `suffix()` allowed as in gstats) | cwinsor.do "WINS-02" |
| WINS-03 | S3 | Reversed `cuts(99 1)`, which winsor2 swaps | r198 | Not a gap: winsor2 1.1 stops with r(111) on reversed cuts; cwinsor's r(198) is gstats winsor's | `cwinsor.ado:94-97` |
| WINS-04 | S2 | gstats: weights `[aw fw iw pw]` | r101 | **Fixed** | `cwinsor.ado:74` |
| WINS-05 | S3 | gstats: `generate(namelist)` | r198 | **Fixed** | — |
| WINS-06 | S3 | gstats: `replace` with `suffix()`/`prefix()` to overwrite existing targets | r198 / r110 | **Fixed** | `cwinsor.ado:104-107,185-189` |
| WINS-07 | S3 | gstats: `nomissing` (drop groups with missing by-keys) | r198 | **Fixed** | — |
| WINS-08 | S4 | gstats: `r(J)`, `r(minJ)`, `r(maxJ)`, `r(cutlow)`, `r(cuthigh)` | Not returned | **Fixed** | — |
| WINS-09 | S4 | gstats: target keeps the source display format | Default format (follows winsor2) | **Fixed** (9 Nov; the source's format with `prefix()`, `generate()` or weights) | cwinsor.do "WINS-09" |
| WINS-10 | S3 | gstats-only `compress`, `forcestrl`, `hashmethod()`, `oncollision()`, `benchmark` | r198 | **Fixed** (9 Nov; accepted and checked, no effect; benchmark shows the verbose timing) | cwinsor.do "WINS-10" |

---

## 3. Differences by design (not counted as gaps)

These are documented intentional differences in default behavior, not missing capability. They are kept
here so they are not reopened as gaps:

- `cwinsor` without `suffix()`/`prefix()` replaces variables in place; winsor2 and gstats create `x_w`.
- `csample`/`cbsample` draws differ from `sample`/`bsample` with the same seed.
- `cmerge` widens long/float combinations to double.
- `civreghdfe` marks absorbed regressors `o.`, where ivreghdfe leaves them unmarked with zero coefficients.
- `civreghdfe` errors on `fuller()`/`kclass()` with `cue`, and on `b0()` with `fuller()`.
- `cppmlhdfe_p` rejects `cooksd`, `hat` and `likelihood`, but so does ppmlhdfe 2.3.3 (confirmed in Stata).

## 4. Documentation fixes found during the audit

These are not functionality gaps but should be corrected:

- `docs/COMPATIBILITY.md` says `cqreg` has "no bootstrap VCE", but `bootstrap: cqreg` worked in testing
  (QREG-12). State what is actually unsupported.
- The compatibility table omits many S1/S3 items above. Once triage is done, either fix them or list them.
- Environment note, not a ctools bug: `C:/Users/Mike/ado/plus/c/cipolate.ado` is a stale ctools 1.2.1 copy
  that calls a plugin subcommand that no longer exists. Delete it.

## 5. Fixed along the way (4 Nov 2026)

Bugs found while fixing the S1 and S2 items, not functionality gaps:

- (9 Nov) `cdestring` with generate() and a repeated source (`u u, gen(a1 a2)`) left the second output
  all missing; outputs are now paired by position and ordered as destring orders them.

- `civreghdfe` (any command using `_ctools_ebv`) aborted with "0bo: operator invalid" when a factor
  interaction cell with a base level was omitted: the omitted name must keep base parts unchanged
  (`3o.race#0b.union`). This stopped the civreghdfe suite at its `i.race#i.union` test.

- `cdestring` printed a stray "invalid syntax" whenever a variable followed the converted one (an inline
  `if` expanded `: word 0 of`), though the results were right.
- `cbinscatter` rejected any `savedata()` path with a directory on Windows (`confirm file "dir/."`
  fails for directories there).
- `cexport excel` could write corrupt cell references (`A01`, read back as column `AHAR`) for strings
  stored as shared strings: the cell prefix was printed with `%s` but never NUL-terminated, so it
  depended on uninitialized heap bytes. The CI plugin happened to read zeros; a zig build failed the
  `oct22` A08/A09 round-trip. Fixed in `src/cexport/cexport_xlsx.c`.

## Open questions (need a Stata run)

Each of these decides whether an item above is a gap, or how severe it is. A single do-file in
`validation/` could settle them all by running the reference and ctools commands side by side.

1. `creghdfe ..., residuals` (bare) and `res` (HDFE-09). Factor-variable depvar (HDFE-19). `ib#.var` in
   `absorb()` (HDFE-20).
   Answered 9 Nov 2026: see HDFE-09, HDFE-19, HDFE-20 (all fixed).
2. `estat ic` after `creghdfe` vs `reghdfe`: do AIC/BIC differ (HDFE-24)? Does the missing `e(footnote)`
   change `bootstrap: creghdfe` output?
   Answered 9 Nov 2026: AIC/BIC differed; fixed with HDFE-24 and HDFE-25.
3. `chdfe` slopes with `L.` operators (CHDFE-05), and the variables kept by `clear` mode (CHDFE-06).
   Answered 9 Nov 2026: CHDFE-05 is not a gap (hdfe 3.2.9 also gives r(101)); CHDFE-06 fixed.
4. `civreghdfe` bare-call error code (IV-01), `partial(_cons)` rc (IV-16). Does ivreghdfe with `absorb()`
   reject iweights (scope of IV-11)? Does ivreghdfe pass group/individual FE through to reghdfe (scope of
   IV-18)?
   Answered 9 Nov 2026: partial(_cons) works in ivreg2 and is r(198) in ivreghdfe with absorb(); ivreghdfe rejects iweights with absorb() (r(101)); ivreghdfe passes group()/individual() to reghdfe but fails (r(2001)). See IV-11, IV-16, IV-18.
5. Does `boottest` accept
   `e(cmd)=="cppmlhdfe"`? Do `mfp` and `bayesboot` accept `cqreg`?
   Answered 9 Nov 2026: boottest fails the same way after ppmlhdfe and cppmlhdfe (r(198), constraints()); bayesboot accepts cqreg; mfp does not (it checks the command name). See QREG-14.
6. Import: URLs for each format (IO-05); `.xlsm`/`.xltx` (IO-11); `clear` with `c(changed)==0` (IO-14);
   sasxport5 `describe, member()` results (IO-13); SAS namelist case sensitivity (IO-12); CRLF on Windows
   export (IO-19).
   Answered 9 Nov 2026: see IO-05, IO-11, IO-12, IO-13, IO-14 and IO-19.
7. `merge ..., force` with a string/numeric key mismatch (MERGE-04); merge's widening message (MERGE-06);
   `_n` merge with an empty master (MERGE-07).
   Answered 9 Nov 2026: see MERGE-04 (not a gap), MERGE-06 and MERGE-07 (fixed).
8. `encode` limit above 65,536 values in Stata 19 (ENC-02); `destring` of hex numbers (DESTR-07); csort
   with frame alias variables (SORT-05).
   Answered 9 Nov 2026: ENC-02 is not a gap; DESTR-07 was a gap and is fixed; SORT-05 is not a gap.
9. `cbinscatter` `n()`, `ms()`, `nbins()` error codes (BIN-07) and binscatter's `nquantiles` range
   (BIN-09). Check the binsreg items against an installed binsreg (BIN-10 to BIN-16).
   Answered 9 Nov 2026: binscatter accepts `n()`, `m()` and `ms()` (r(0)); `nquantiles()` must be at
   least 2, and binscatter fails with r(198) above 100 quantile bins. binsreg items are documented only.
10. `cpsmatch` `_nn` for on-support controls (PSM-14); whether psmatch2's `r()` survives the command.
11. Long-strL errors in `cbinscatter absorb()`, `crangestat by()`, `cwinsor by()`, `csample by()`,
    `cbsample strata()/cluster()`.

## Suggested work order

1. **Done (4 Nov 2026):** fix or reject every S1 item (section 1). Most need a small ado change: either implement the option
   or make it error. DESTR-01, HDFE-10, PPML-01 and the civreghdfe auxiliary-regression issues (IV-19,
   IV-20) change results and come first.
2. **Run the open-questions do-file** and reclassify, close, or confirm the **Run** items.
3. **Cheap S3 fixes in bulk:** option abbreviations (HDFE-09/11/21, CHDFE-01/02, IV-21, IO-06/10/15/16/18,
   QREG-08, BIN-07, WINS-02), `stable`/`in` for csort, `quantile()` as a percent, the vce suboption syntax
   for cqreg, and `replace`/overwrite semantics for `savefe`/`resid`. Add each to `validation/` option
   coverage.
4. **S2 status (9 Nov 2026):** every S2 item is fixed, reclassified, or listed here with the reason it
   remains. HDFE-05 (group and individual fixed effects), IV-02 (civreghdfe slopes) and IO-01's Brotli
   codec were done on 9 Nov. Left: BIN-10/11 (binsreg is not installed to verify against) and PSM-05
   (the `spline` package psmatch2 calls, STB-24 `snp7_1`, is not installed here). Regression
   contracts are the `s2_*` cases (validation/stata/s2_*.do) and `nov04_s2`.
   slopes (HDFE-01/02/03), cqreg weights and `by:` (QREG-01/03), civreghdfe weak-ID statistics and replay
   (IV-25, IV-01), cpsmatch multiple outcomes, `ate`, neighbor IDs (PSM-01/02/12), crangestat variable
   bounds and cov/corr/reg (RSTAT-01/02), the long-strL plugin limit (one shared fix in
   `ctools_transport.c` covers DIST/RSTAT/SAMP/BSAMP/DESTR/IPOL), cbinscatter multiple y. Mark the ones
   declined as "won't do" in COMPATIBILITY.md.
5. **S4 stored results and output:** do these alongside the related S2 work. Prioritize results other
   packages read (`e(estat_cmd)`, `e(footnote)`, ivreg2 `e()` names, psmatch2 variables).
