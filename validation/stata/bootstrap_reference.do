* bootstrap_reference: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "bootstrap_reference"

if "`__case'" == "bootstrap_reference" {
* Reference behavior that the bootstrap options of creghdfe, cppmlhdfe and
* civreghdfe are built against (docs/development/BOOTSTRAP_PLAN.md, Phase 0).
* These tests run Stata's own commands only: they pin the conventions the
* plan relies on, so that a change in Stata shows up here first.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running bootstrap reference checks..."

* SECTION 1: the Mata emulation reproduces bsample's draws
print_section "bsample draws reproduced in Mata"

foreach design in obs strata cluster strata_cluster {
    sysuse auto, clear
    keep price mpg foreign
    gen long id = _n
    gen cl = mod(_n * 7, 11)
    local opts
    if inlist("`design'", "strata", "strata_cluster") local opts `opts' strata(foreign)
    if inlist("`design'", "cluster", "strata_cluster") local opts `opts' cluster(cl)
    local bsopts = cond("`opts'" == "", "", ", `opts'")
    set seed 4242
    bsref_draws, reps(3) generate(em) `opts'
    local state_emulated = c(rngstate)
    set seed 4242
    forvalues r = 1/3 {
        gen double bm`r' = .
        sort id
        local wopt = cond("`opts'" == "", ", weight(bm`r')", ", `opts' weight(bm`r')")
        bsample `wopt'
    }
    local state_bsample = c(rngstate)
    sort id
    local differ = 0
    forvalues r = 1/3 {
        count if em`r' != bm`r'
        local differ = `differ' + r(N)
    }
    if `differ' == 0 & "`state_emulated'" == "`state_bsample'" ///
        test_pass "bsample `design' draws and final RNG state reproduced"
    else test_fail "bsample `design' draws" "`differ' weights differ; RNG states equal: `=("`state_emulated'" == "`state_bsample'")'"
}

* SECTION 2: the bootstrap prefix draws by successive bsample calls on
* e(sample), so frequency-weighted fits on the emulated draws reproduce its
* replicates
print_section "bootstrap prefix replicates"

sysuse auto, clear
keep price mpg weight
gen cl = mod(_n * 7, 11)
tempfile prefix_obs prefix_cl
set seed 777
bootstrap _b, reps(4) saving(`prefix_obs', replace double) nodots: regress price mpg
local state_prefix = c(rngstate)
set seed 777
bsref_draws, reps(4) generate(m)
local state_emulated = c(rngstate)
local maxdiff = 0
forvalues r = 1/4 {
    regress price mpg [fw=m`r']
    local b`r' = _b[mpg]
}
use `prefix_obs', clear
forvalues r = 1/4 {
    local maxdiff = max(`maxdiff', reldif(_b_mpg[`r'], `b`r''))
}
if `maxdiff' < 1e-12 & "`state_prefix'" == "`state_emulated'" ///
    test_pass "observation bootstrap: prefix replicates and RNG state reproduced"
else test_fail "observation bootstrap prefix" "max reldif `maxdiff'; RNG states equal: `=("`state_prefix'" == "`state_emulated'")'"

sysuse auto, clear
keep price mpg weight
gen cl = mod(_n * 7, 11)
set seed 778
bootstrap _b, reps(4) cluster(cl) idcluster(nc) saving(`prefix_cl', replace double) nodots: regress price mpg
set seed 778
bsref_draws, reps(4) generate(m) cluster(cl)
local maxdiff = 0
forvalues r = 1/4 {
    regress price mpg [fw=m`r']
    local b`r' = _b[mpg]
}
use `prefix_cl', clear
forvalues r = 1/4 {
    local maxdiff = max(`maxdiff', reldif(_b_mpg[`r'], `b`r''))
}
if `maxdiff' < 1e-12 test_pass "cluster bootstrap: prefix replicates reproduced"
else test_fail "cluster bootstrap prefix" "max reldif `maxdiff'"

* The prefix resamples e(sample): rows a command drops (here reghdfe's
* singletons) are not in the resampling frame
capture which reghdfe
if _rc test_skip "reghdfe prefix frame" "reghdfe is not installed"
else {
    sysuse auto, clear
    keep price mpg
    gen grp = cond(_n <= 3, 100 + _n, mod(_n, 9))
    tempfile prefix_hdfe
    set seed 779
    bootstrap _b, reps(3) saving(`prefix_hdfe', replace double) nodots: reghdfe price mpg, absorb(grp)
    reghdfe price mpg, absorb(grp)
    keep if e(sample)
    set seed 779
    bsref_draws, reps(3) generate(m)
    local maxdiff = 0
    forvalues r = 1/3 {
        reghdfe price mpg [fw=m`r'], absorb(grp)
        local b`r' = _b[mpg]
    }
    use `prefix_hdfe', clear
    forvalues r = 1/3 {
        local maxdiff = max(`maxdiff', reldif(_b_mpg[`r'], `b`r''))
    }
    if `maxdiff' < 1e-9 test_pass "reghdfe prefix resamples e(sample) (singletons excluded)"
    else test_fail "reghdfe prefix frame" "max reldif `maxdiff'"
}

* Replication weights from rwgen, replayed with bootstrap's fweights()
capture which rwgen
if _rc test_skip "rwgen replay" "rwgen needs StataNow"
else {
    sysuse auto, clear
    keep price mpg
    gen cl = mod(_n * 7, 11)
    tempfile replay
    rwgen bsample rw, reps(4) rseed(778) cluster(cl)
    bootstrap _b, fweights(rw*) saving(`replay', replace double) nodots: regress price mpg
    use `replay', clear
    rename _b_mpg replay_mpg
    merge 1:1 _n using `prefix_cl', nogenerate keepusing(_b_mpg)
    count if reldif(replay_mpg, _b_mpg) > 1e-12
    if r(N) == 0 test_pass "rwgen bsample + bootstrap fweights() replays the prefix"
    else test_fail "rwgen replay" "`r(N)' replicates differ"
}

* SECTION 3: what vce(bootstrap) posts and how it treats failures
print_section "vce(bootstrap) conventions"

foreach w in fw aw pw {
    sysuse auto, clear
    capture regress price mpg [`w'=rep78], vce(bootstrap, reps(10) seed(1))
    if _rc == 101 test_pass "vce(bootstrap) refuses `w' weights (r(101))"
    else test_fail "vce(bootstrap) with `w' weights" "rc=`=_rc', expected 101"
}

* A replicate whose omitted-coefficient pattern differs from the full sample
* is dropped whole; e(V) uses complete replicates only
sysuse auto, clear
gen rare = _n <= 2
tempfile reps_rare
regress price mpg rare, vce(bootstrap, reps(50) seed(1) saving(`reps_rare', replace))
matrix V = e(V)
matrix R = e(reps)
local nreps = e(N_reps)
local nmis = e(N_misreps)
local has_vmb = 0
capture confirm matrix e(V_modelbased)
if !_rc local has_vmb = 1
local allsame = (R[1,1] == `nreps' & R[1,2] == `nreps' & R[1,3] == `nreps')
use `reps_rare', clear
egen nmiss = rowmiss(_b_mpg _b_rare _b_cons)
count if nmiss > 0 & nmiss < 3
local partial = r(N)
drop if nmiss
correlate _b_mpg _b_rare _b_cons, covariance
matrix C = r(C)
local vdiff = mreldif(C, V)
if `nmis' > 0 & `partial' == 0 & `allsame' & `vdiff' < 1e-10 & `has_vmb' ///
    test_pass "vce(bootstrap): whole replicates dropped, casewise e(V), e(V_modelbased) posted"
else test_fail "vce(bootstrap) failures" "N_misreps=`nmis' partial=`partial' allsame=`allsame' mreldif=`vdiff' V_modelbased=`has_vmb'"

* Replicates are posted as floats unless saving(, double) is given, and e(V)
* is computed from what was posted: the default e(V) is the covariance of the
* float-rounded replicates
sysuse auto, clear
tempfile reps_double
regress price mpg weight, vce(bootstrap, reps(20) seed(3))
matrix V_default = e(V)
regress price mpg weight, vce(bootstrap, reps(20) seed(3) saving(`reps_double', replace double))
matrix V_double = e(V)
use `reps_double', clear
foreach v of varlist _b_* {
    generate float f`v' = `v'
}
correlate f_b_*, covariance
matrix C = r(C)
correlate _b_*, covariance
matrix D = r(C)
if mreldif(C, V_default) < 1e-12 & mreldif(D, V_double) < 1e-12 ///
    test_pass "vce(bootstrap): e(V) is computed from float replicates unless saving(, double)"
else test_fail "float replicates" "float mreldif `=mreldif(C, V_default)'; double mreldif `=mreldif(D, V_double)'"

* mse centers the variance on e(b) and divides by the number of replicates
sysuse auto, clear
tempfile reps_mse
regress price mpg, vce(bootstrap, reps(20) seed(2) mse saving(`reps_mse', replace double))
local se_mse = _se[mpg]
local b_hat = _b[mpg]
use `reps_mse', clear
gen double dev2 = (_b_mpg - `b_hat')^2
summarize dev2, meanonly
if reldif(sqrt(r(mean)), `se_mse') < 1e-10 test_pass "mse: sqrt(mean((b_r - b)^2))"
else test_fail "mse" "`=sqrt(r(mean))' vs `se_mse'"

* SECTION 4: Stata's wild cluster bootstrap
print_section "wildbootstrap conventions"

if c(stata_version) < 18 test_skip "wildbootstrap" "needs Stata 18"
else {
    clear
    set seed 2024
    set obs 600
    gen g = ceil(_n / 20)
    gen h = ceil(_n / 60)
    gen a = mod(_n, 17) + 1
    gen x = rnormal() + 0.1 * a
    gen y = 0.04 * x + 0.05 * a + rnormal() + 0.5 * rnormal() * (g / 30)
    foreach v in y x {
        bysort a: egen double m_`v' = mean(`v')
        gen double d_`v' = `v' - m_`v'
    }
    sort g a, stable

    wildbootstrap regress y x, cluster(g) rseed(1) noci
    if e(N_wbreps) == 1000 & "`e(wb_weight)'" == "rademacher" & "`e(wb_ptype)'" == "equal-tailed" ///
        test_pass "wildbootstrap defaults: 1000 reps, Rademacher, equal-tailed"
    else test_fail "wildbootstrap defaults" "reps=`e(N_wbreps)' weight=`e(wb_weight)' ptype=`e(wb_ptype)'"

    * Effect nested in the clusters: absorbing it equals adding dummies
    wildbootstrap areg y x, absorb(g) cluster(h) rseed(1) reps(1000)
    matrix A = e(wboot)
    wildbootstrap regress y x i.g, cluster(h) coefficients(x) rseed(1) reps(1000)
    matrix D = e(wboot)
    if mreldif(A, D) < 1e-7 test_pass "wildbootstrap areg, nested effect: equals dummies"
    else test_fail "wildbootstrap nested" "mreldif `=mreldif(A, D)'"

    * Effect not nested in the clusters: areg bootstraps the partialled
    * regression (the FE projection is not applied to the perturbed
    * residuals); dummies would give the full-model bootstrap instead
    wildbootstrap areg y x, absorb(a) cluster(g) rseed(1) reps(1000)
    matrix A = e(wboot)
    wildbootstrap regress d_y d_x, cluster(g) coefficients(d_x) rseed(1) reps(1000)
    matrix F = e(wboot)
    * columns 3-5: p-value and bounds (t differs by the small-sample factor)
    matrix A = A[1, 3..5]
    matrix F = F[1, 3..5]
    if mreldif(A, F) < 1e-7 ///
        test_pass "wildbootstrap areg, effect crossing clusters: bootstraps the partialled regression"
    else test_fail "wildbootstrap non-nested" "mreldif `=mreldif(A, F)'"
}

* SECTION 5: boottest (SSC). With 8 clusters and Rademacher weights it
* enumerates all 2^8 sign vectors, so these checks need no seed.
print_section "boottest conventions"

capture which boottest
local has_bt = !_rc
capture which reghdfe
local has_hdfe = !_rc
capture which ivreghdfe
local has_iv = !_rc
if !`has_bt' | !`has_hdfe' test_skip "boottest" "boottest or reghdfe is not installed"
else {
    clear
    set seed 2024
    set obs 400
    gen g = ceil(_n / 50)
    gen a = mod(_n, 13) + 1
    gen x = rnormal() + 0.1 * a
    gen y = 0.05 * x + 0.05 * a + rnormal() * (1 + g / 8)
    sort g a, stable
    regress y x i.a, cluster(g)
    boottest x, reps(999) noci
    local p_dummies = r(p)
    local enum = (r(reps) == 256)
    reghdfe y x, absorb(a) cluster(g)
    boottest x, reps(999) noci
    local p_hdfe = r(p)
    * Effect crossing the clusters: boottest bootstraps the full model (as
    * with dummies), unlike Stata's wildbootstrap areg (SECTION 4)
    if `enum' & reldif(`p_dummies', `p_hdfe') < 1e-12 ///
        test_pass "boottest after reghdfe, effect crossing clusters: full-model bootstrap, enumerated"
    else test_fail "boottest non-nested" "enumerated=`enum' p dummies=`p_dummies' reghdfe=`p_hdfe'"

    if `has_iv' {
        gen z1 = rnormal()
        gen z2 = rnormal()
        gen w = rnormal()
        gen xe = 0.5 * z1 + 0.3 * z2 + 0.2 * w + rnormal()
        ivreghdfe y w (xe = z1 z2), absorb(a) cluster(g)
        capture boottest xe, reps(999) noci
        local rc_endog = _rc
        capture boottest w, reps(999) noci
        local rc_exog = _rc
        ivreghdfe y w (xe = z1 z2), absorb(a) cluster(g) liml
        capture boottest xe, reps(999) noci
        local rc_liml = _rc
        ivreghdfe y w (xe = z1 z2), absorb(a g) cluster(g)
        capture boottest xe, reps(999) noci
        local rc_two = _rc
        if !`rc_endog' & !`rc_exog' & !`rc_liml' & `rc_two' == 198 ///
            test_pass "boottest WRE after ivreghdfe: endogenous, exogenous, LIML; refuses two absorbed effects"
        else test_fail "boottest WRE" "rc endog=`rc_endog' exog=`rc_exog' liml=`rc_liml' two effects=`rc_two'"
    }
    else test_skip "boottest WRE" "ivreghdfe is not installed"
}

noi print_summary "bootstrap_reference"
}

global CTOOLS_COMPONENT_COMPLETE "bootstrap_reference"
    global CTOOLS_CASE_COMPLETE "bootstrap_reference"
    exit
}

if "`__case'" == "bootstrap_helpers" {
* The shared bootstrap helper (build/_ctools_hdfe_boot.ado): option parsing,
* seeds, and the Stata-RNG uniforms behind draws(stata). The uniforms are
* mapped to multiplicities with the reference mapping, which the plugin's
* C mapping matches (validation/native/boot.suite.c).

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running bootstrap helper checks..."

print_section "Option parsing"

_ctools_hdfe_boot parse, spec(bootstrap)
if "`r(kind)'" == "bootstrap" & r(reps) == 50 & "`r(draws)'" == "stata" & r(dots) == 1 & ///
    r(level) == c(level) & r(mse) == 0 & r(designopts) == 0 ///
    test_pass "parse: defaults (50 reps, draws(stata))"
else test_fail "parse defaults" "kind=`r(kind)' reps=`r(reps)' draws=`r(draws)' designopts=`r(designopts)'"

sysuse auto, clear
_ctools_hdfe_boot parse, spec(`"bootstrap, reps(200) seed(12) cluster(foreign) strata(rep78) mse level(90) nodots draws(stata)"')
if r(reps) == 200 & "`r(seed)'" == "12" & "`r(cluster)'" == "foreign" & "`r(strata)'" == "rep78" & ///
    r(mse) == 1 & r(level) == 90 & r(dots) == 0 & "`r(draws)'" == "stata" & r(designopts) == 1 ///
    test_pass "parse: every bootstrap option"
else test_fail "parse options" "reps=`r(reps)' seed=`r(seed)' cluster=`r(cluster)' strata=`r(strata)' draws=`r(draws)'"

_ctools_hdfe_boot parse, spec(robust)
if "`r(kind)'" == "robust" test_pass "parse: other vce() kinds pass through"
else test_fail "parse robust" "kind=`r(kind)'"

_ctools_hdfe_boot parse, spec(`"boot, sestub(se_) ci nreps dse(dse)"') fe
if "`r(sestub)'" == "se_" & r(ci) == 1 & r(nreps) == 1 & "`r(dse)'" == "dse" ///
    test_pass "parse: fixed-effect options with fe"
else test_fail "parse fe options" "sestub=`r(sestub)' ci=`r(ci)' dse=`r(dse)'"

foreach bad in "bootstrap, reps(1)" "bootstrap, cluster(foreign rep78)" "bootstrap, draws(slow)" ///
    "bootstrap, sestub(se_)" "bootstrap, nodots dots(5)" "bootstrap reps(5)" {
    capture _ctools_hdfe_boot parse, spec(`"`bad'"')
    if _rc == 198 test_pass "parse refuses: `bad'"
    else test_fail "parse refuses: `bad'" "rc=`=_rc'"
}

print_section "Seeds and Stata-RNG uniforms"

* seed() sets the seed before r(rngstate) is taken; fast takes two draws
set seed 5
local after_seed = c(rngstate)
scalar x1 = runiform()
scalar x2 = runiform()
local after_two = c(rngstate)
set seed 999
_ctools_hdfe_boot seed, seed(5) fast
if `"`r(rngstate)'"' == `"`after_seed'"' & `"`c(rngstate)'"' == `"`after_two'"' & ///
    "`r(seedhi)'" == string(floor(scalar(x1) * 2147483647), "%10.0f") & ///
    "`r(seedlo)'" == string(floor(scalar(x2) * 2147483647), "%10.0f") ///
    test_pass "seed: seed() first, r(rngstate) before draws, master seed from two draws"
else test_fail "seed" "rngstate or master seed differs"

* The uniforms, mapped as bsample maps them, give bsample's draws and leave
* the generator where bsample leaves it
foreach design in obs cluster {
    sysuse auto, clear
    keep price mpg foreign
    gen long id = _n
    gen cl = mod(_n * 7, 11)
    local opts = cond("`design'" == "cluster", "cluster(cl)", "")
    if "`design'" == "cluster" {
        sort cl, stable
        by cl: gen double cidx = (_n == 1)
        quietly count if cidx
        local units = r(N)
        replace cidx = sum(cidx)
    }
    else local units = _N
    gen double one = 1
    set seed 6161
    _ctools_hdfe_boot uniforms, units(`units') reps(3) generate(u)
    local state_helper = c(rngstate)
    * packed: replicate r's draws from position 2 (r - 1) units, then its
    * tie-breaks (in the cluster design, both in u1; for observations, one
    * variable each)
    forvalues r = 1/3 {
        local qd = 2 * (`r' - 1) * `units'
        local qt = `qd' + `units'
        if "`design'" == "cluster" {
            mata: st_store(., st_addvar("double", "hm`r'"), _bsref_clust(st_data(., "one"), ///
                st_data(., "cidx"), _bsref_packed("u", `qd', `units'), _bsref_packed("u", `qt', `units')))
        }
        else {
            mata: st_store(., st_addvar("double", "hm`r'"), _bsref_obs(st_data(., "one"), ///
                _bsref_packed("u", `qd', `units'), _bsref_packed("u", `qt', `units')))
        }
    }
    set seed 6161
    forvalues r = 1/3 {
        gen double bw`r' = .
        sort id
        if "`design'" == "cluster" sort cl, stable
        bsample, `opts' weight(bw`r')
    }
    local state_bsample = c(rngstate)
    local differ = 0
    forvalues r = 1/3 {
        count if hm`r' != bw`r'
        local differ = `differ' + r(N)
    }
    if `differ' == 0 & "`state_helper'" == "`state_bsample'" ///
        test_pass "uniforms: `design' bootstrap reproduces bsample and its RNG state"
    else test_fail "uniforms `design'" "`differ' weights differ"
}

noi print_summary "bootstrap_helpers"
}

global CTOOLS_COMPONENT_COMPLETE "bootstrap_helpers"
    global CTOOLS_CASE_COMPLETE "bootstrap_helpers"
    exit
}

if "`__case'" == "creghdfe_bootstrap" {
* creghdfe, vce(bootstrap) (docs/development/BOOTSTRAP_PLAN.md, Phase 2b).
* With draws(stata) and the same seed, creghdfe must reproduce Stata's
* bootstrap prefix around reghdfe: every summary, the number of complete
* replicates, and the random-number state it leaves.

do "validation/support/validate_setup.do"

capture program drop _bs_compare
program define _bs_compare
    * Compares the creghdfe results (estimates ct) with the prefix's (pf)
    args name state_ct state_pf
    tempname a b
    local worst = 0
    foreach m in V b_bs se bias z0 ci_normal ci_percentile ci_bc {
        estimates restore ct
        matrix `a' = e(`m')
        estimates restore pf
        matrix `b' = e(`m')
        local worst = max(`worst', mreldif(`a', `b'))
    }
    estimates restore ct
    local n_ct = e(N_reps)
    local mis_ct = e(N_misreps)
    estimates restore pf
    if `worst' < 1e-10 & `n_ct' == e(N_reps) & `mis_ct' == e(N_misreps) & ///
        `"`state_ct'"' == `"`state_pf'"' ///
        test_pass "draws(stata) reproduces bootstrap: reghdfe (`name'; max reldif `worst')"
    else test_fail "draws(stata) `name'" "max reldif `worst'; N_reps `n_ct'/`=e(N_reps)'; N_misreps `mis_ct'/`=e(N_misreps)'; RNG states equal: `=(`"`state_ct'"' == `"`state_pf'"')'"
end

quietly {

noi di as text "Running creghdfe vce(bootstrap) checks..."


capture which reghdfe
if _rc test_skip "creghdfe bootstrap parity" "reghdfe is not installed"
else {
    print_section "Parity with the bootstrap prefix"

    tempfile pf

    sysuse auto, clear
    creghdfe price mpg weight, absorb(foreign rep78) vce(bootstrap, reps(30) seed(5) draws(stata) nodots)
    local s_ct = c(rngstate)
    estimates store ct
    bootstrap _b, reps(30) seed(5) saving(`pf', replace double) nodots: reghdfe price mpg weight, absorb(foreign rep78)
    local s_pf = c(rngstate)
    estimates store pf
    _bs_compare "observations" `"`s_ct'"' `"`s_pf'"'

    sysuse auto, clear
    gen cl = mod(_n * 7, 11)
    creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(30) seed(6) cluster(cl) draws(stata) nodots)
    local s_ct = c(rngstate)
    estimates store ct
    bootstrap _b, reps(30) seed(6) cluster(cl) idcluster(nc) saving(`pf', replace double) nodots: reghdfe price mpg weight, absorb(foreign)
    local s_pf = c(rngstate)
    estimates store pf
    _bs_compare "clusters" `"`s_ct'"' `"`s_pf'"'

    * strata, clusters that cross strata, and singletons outside e(sample)
    sysuse auto, clear
    gen cl = mod(_n * 7, 11)
    gen grp = cond(_n <= 3, 100 + _n, mod(_n, 9))
    creghdfe price mpg weight, absorb(grp) vce(bootstrap, reps(25) seed(7) cluster(cl) strata(foreign) draws(stata) nodots)
    local s_ct = c(rngstate)
    estimates store ct
    bootstrap _b, reps(25) seed(7) cluster(cl) strata(foreign) saving(`pf', replace double) nodots: reghdfe price mpg weight, absorb(grp)
    local s_pf = c(rngstate)
    estimates store pf
    _bs_compare "strata, clusters, singletons" `"`s_ct'"' `"`s_pf'"'

    * a dummy that drops out of some replicates: those are missing
    sysuse auto, clear
    gen rare = _n <= 3
    creghdfe price mpg rare, absorb(foreign) vce(bootstrap, reps(40) seed(8) draws(stata) nodots)
    local s_ct = c(rngstate)
    estimates store ct
    bootstrap _b, reps(40) seed(8) saving(`pf', replace double) nodots: reghdfe price mpg rare, absorb(foreign)
    local s_pf = c(rngstate)
    estimates store pf
    _bs_compare "collinear replicates" `"`s_ct'"' `"`s_pf'"'

    * aweights (the prefix needs force)
    sysuse auto, clear
    creghdfe price mpg weight [aw=displacement], absorb(foreign) vce(bootstrap, reps(20) seed(9) draws(stata) nodots)
    local s_ct = c(rngstate)
    estimates store ct
    bootstrap _b, reps(20) seed(9) force saving(`pf', replace double) nodots: reghdfe price mpg weight [aw=displacement], absorb(foreign)
    local s_pf = c(rngstate)
    estimates store pf
    _bs_compare "aweights" `"`s_ct'"' `"`s_pf'"'
}

print_section "Results and reproducibility"

* draws(fast): the same seed gives the same results; another seed does not
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(40) seed(10) draws(fast) nodots)
matrix A = e(V)
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(40) seed(10) draws(fast) nodots)
matrix B = e(V)
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(40) seed(11) draws(fast) nodots)
if mreldif(A, B) == 0 & mreldif(A, e(V)) > 0 test_pass "draws(fast): reproducible from the seed"
else test_fail "draws(fast) seed" "same-seed reldif `=mreldif(A, B)'"

* what vce(bootstrap) posts, and the fit's own covariance as e(V_modelbased)
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign)
matrix Vols = e(V)
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(20) seed(12) nodots)
local ok = ("`e(vce)'" == "bootstrap" & "`e(vcetype)'" == "Bootstrap" & "`e(prefix)'" == "bootstrap" ///
    & e(N_reps) == 20 & e(N_misreps) == 0 & e(df_r) == . & e(chi2) < . & e(bs_version) == 3 ///
    & "`e(rngstate)'" != "" & "`e(draws)'" == "stata")
foreach m in b_bs reps bias se z0 ci_normal ci_percentile ci_bc boot_status {
    capture confirm matrix e(`m')
    if _rc local ok = 0
}
if `ok' & mreldif(e(V_modelbased), Vols) < 1e-12 ///
    test_pass "vce(bootstrap) posts e() as Stata's bootstrap does, with e(V_modelbased)"
else test_fail "vce(bootstrap) e() results" "vce=`e(vce)' vcetype=`e(vcetype)' N_reps=`e(N_reps)'"

capture noisily estat bootstrap, all
if !_rc test_pass "estat bootstrap after creghdfe, vce(bootstrap)"
else test_fail "estat bootstrap" "rc=`=_rc'"
capture creghdfe
if !_rc & "`e(vce)'" == "bootstrap" test_pass "replay after vce(bootstrap)"
else test_fail "replay" "rc=`=_rc'"

* mse: V about e(b), divided by the number of replicates
sysuse auto, clear
creghdfe price mpg, absorb(foreign) vce(bootstrap, reps(30) seed(13) nodots)
local se_mean = _se[mpg]
creghdfe price mpg, absorb(foreign) vce(bootstrap, reps(30) seed(13) mse nodots)
local se_mse = _se[mpg]
local bias = el(e(bias), 1, 1)
* sum((b_r - b)^2)/R = (R-1)/R * sd^2 + bias^2
if reldif((`se_mse')^2, 29/30 * (`se_mean')^2 + (`bias')^2) < 1e-10 test_pass "mse: centered on e(b), divided by R"
else test_fail "mse" "`se_mse' vs `se_mean' and bias `bias'"

print_section "Errors"

sysuse auto, clear
gen cl = mod(_n, 5)
foreach bad in "bootstrap, reps(1)" "bootstrap, cluster(cl foreign)" ///
    "bootstrap, draws(slow)" "bootstrap, reps(3) _bsdraws(cl foreign)" ///
    "bootstrap, reps(2) draws(stata) _bsdraws(cl foreign)" {
    capture creghdfe price mpg, absorb(foreign) vce(`bad')
    if _rc == 198 test_pass "creghdfe refuses vce(`bad')"
    else test_fail "creghdfe vce(`bad')" "rc=`=_rc'"
}
capture creghdfe price mpg, absorb(foreign, savefe) vce(bootstrap, reps(5) seed(1) saving(x)) fevce(bootstrap, saving(y))
if _rc == 198 test_pass "saving() goes in vce(bootstrap) only, when both are given"
else test_fail "saving() twice" "rc=`=_rc'"

* fweighted rows: the observations they stand for are drawn, and a row
* weighs what was drawn of it (a weight of draws times fweight centered the
* replicates on the fweight-squared estimator, -526 here against -302)
sysuse auto, clear
gen fw = cond(mod(_n, 4) == 0, 30, 1)
creghdfe price mpg weight [fw=fw], absorb(foreign) vce(bootstrap, reps(300) seed(1) nodots)
local bias = el(e(bias), 1, 1)
local se = _se[mpg]
expand fw
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(300) seed(1) nodots)
if abs(`bias') < 4 * `se' / sqrt(300) & abs(`se' / _se[mpg] - 1) < 0.2 ///
    test_pass "fweighted rows: replicates centered on e(b), SEs as the expanded data's"
else test_fail "fweighted rows" "bias `bias' se `se' vs expanded `=_se[mpg]'"

noi print_summary "creghdfe_bootstrap"
}

global CTOOLS_COMPONENT_COMPLETE "creghdfe_bootstrap"
    global CTOOLS_CASE_COMPLETE "creghdfe_bootstrap"
    exit
}

if "`__case'" == "creghdfe_fe_bootstrap" {
* creghdfe, fevce(bootstrap) (docs/development/BOOTSTRAP_PLAN.md, Phase 2c).
* Oracle: the same draws (draws(stata), bsref_draws), reghdfe with
* absorb(fe=...) on each replicate's data expanded by its multiplicities, and
* the standard deviation of each level's effect over the replicates that
* estimate it.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

capture program drop _fe_oracle
program define _fe_oracle
    * Per observation: the SD over replicates of each effect's level value
    * (o_se1..), the number of replicates estimating it (o_n1..), and the SD
    * of their sum (o_dse). Data: e(sample) only; m1..mR the draws.
    syntax, reps(integer) absorb(string) model(string) [weight(string)]
    local G : word count `absorb'
    local wexp = cond("`weight'" == "", "", "[aw=`weight']")
    forvalues r = 1/`reps' {
        local fes
        forvalues g = 1/`G' {
            local fes `fes' f`g'_`r'=`: word `g' of `absorb''
        }
        * The replicate as the bootstrap prefix builds it: each row repeated
        * as often as it was drawn (so a level drawn once, twice, is no
        * singleton), then each level's effect back on every row of it
        tempfile full expanded
        forvalues g = 1/`G' {
            tempfile lv`g'
        }
        quietly save `full'
        quietly keep if m`r' > 0
        quietly expand m`r'
        quietly reghdfe `model' `wexp', absorb(`fes') tolerance(1e-12)
        quietly save `expanded'
        forvalues g = 1/`G' {
            local v : word `g' of `absorb'
            use `v' f`g'_`r' using `expanded', clear
            quietly collapse (firstnm) f`g'_`r', by(`v')
            quietly save `lv`g''
        }
        use `full', clear
        local dsum
        forvalues g = 1/`G' {
            local v : word `g' of `absorb'
            quietly merge m:1 `v' using `lv`g'', nogenerate
            local dsum `dsum' + f`g'_`r'
        }
        quietly generate double d_`r' = 0 `dsum'
    }
    forvalues g = 1/`G' {
        quietly egen double o_se`g' = rowsd(f`g'_*)
        quietly egen double o_n`g' = rownonmiss(f`g'_*)
    }
    quietly egen double o_dse = rowsd(d_*)
end

quietly {

noi di as text "Running creghdfe fevce(bootstrap) checks..."

capture which reghdfe
if _rc test_skip "creghdfe fevce oracle" "reghdfe is not installed"
else {
    print_section "Fixed-effect standard errors against the oracle"

    * One effect, cluster bootstrap; some levels drop out of replicates
    sysuse auto, clear
    gen cl = mod(_n * 7, 11)
    creghdfe price mpg weight, absorb(rep78, savefe) ///
        fevce(bootstrap, reps(25) seed(21) cluster(cl) draws(stata) nreps dse(dse) nodots)
    local fe_reps = e(fe_N_reps)
    keep if e(sample)
    keep price mpg weight rep78 cl __hdfe1__ __fese1__ __fenrep1__ dse
    set seed 21
    bsref_draws, reps(25) generate(m) cluster(cl)
    _fe_oracle, reps(25) absorb(rep78) model(price mpg weight)
    count if reldif(__fese1__, o_se1) > 1e-6 | __fenrep1__ != o_n1 | reldif(dse, o_dse) > 1e-6
    if r(N) == 0 & `fe_reps' == 25 ///
        test_pass "one effect: SE, valid replicates and d's SE equal the oracle"
    else test_fail "one effect" "`r(N)' rows differ"

    * Two crossed effects, observation bootstrap, aweights
    clear
    set seed 99
    set obs 320
    gen id = ceil(_n / 8)
    gen yr = mod(_n - 1, 8) + 1
    gen x = rnormal() + id / 40
    gen w = 0.5 + runiform()
    gen y = 0.4 * x + id / 10 + yr / 5 + rnormal()
    creghdfe y x [aw=w], absorb(id yr, savefe) ///
        fevce(bootstrap, reps(20) seed(22) draws(stata) nreps dse(dse) nodots)
    local warned = e(fe_norm_warning)
    local splits = e(fe_split_reps)
    keep if e(sample)
    set seed 22
    bsref_draws, reps(20) generate(m)
    _fe_oracle, reps(20) absorb(id yr) model(y x) weight(w)
    count if reldif(__fese1__, o_se1) > 1e-6 | reldif(__fese2__, o_se2) > 1e-6 | ///
        __fenrep1__ != o_n1 | __fenrep2__ != o_n2 | reldif(dse, o_dse) > 1e-6
    if r(N) == 0 & `warned' == 1 & `splits' == 0 ///
        test_pass "two effects: per-level and d SEs equal the oracle; normalization warning set"
    else test_fail "two effects" "`r(N)' rows differ; warning `warned'; splits `splits'"
}

print_section "Modes, outputs and errors"

* Both targets from one set of replicates: e(V) as with vce(bootstrap) alone,
* the SEs as with fevce(bootstrap) alone
sysuse auto, clear
creghdfe price mpg weight, absorb(rep78, savefe) vce(bootstrap, reps(20) seed(23) nodots) fevce(bootstrap)
matrix Vboth = e(V)
rename __fese1__ se_both
drop __hdfe1__
creghdfe price mpg weight, absorb(rep78) vce(bootstrap, reps(20) seed(23) nodots)
matrix Vcoef = e(V)
creghdfe price mpg weight, absorb(rep78, savefe) vce(robust) fevce(bootstrap, reps(20) seed(23) nodots)
count if reldif(se_both, __fese1__) > 1e-12
if mreldif(Vboth, Vcoef) < 1e-12 & r(N) == 0 & "`e(vce)'" == "robust" ///
    test_pass "vce() and fevce() share replicates; fevce alone keeps the analytic e(V)"
else test_fail "shared replicates" "mreldif `=mreldif(Vboth, Vcoef)'; vce=`e(vce)'"

* Names: sestub(), and no warning with one effect
sysuse auto, clear
creghdfe price mpg, absorb(rep78, savefe) fevce(bootstrap, reps(10) seed(24) sestub(fese) nreps nodots)
capture confirm variable fese1 fese1_nrep
if !_rc & e(fe_norm_warning) == 0 & "`e(fe_se)'" == "fese1" test_pass "sestub() names; no warning with one effect"
else test_fail "sestub" "rc=`=_rc' fe_se=`e(fe_se)'"

foreach bad in "fevce(bootstrap, reps(5) seed(1))" ///
    "absorb(rep78, savefe) fevce(robust)" ///
    "absorb(rep78, savefe) vce(bootstrap, reps(10) seed(1)) fevce(bootstrap, reps(20))" {
    sysuse auto, clear
    local abs = cond(strpos(`"`bad'"', "absorb("), "", "absorb(rep78)")
    capture creghdfe price mpg, `abs' `bad'
    if _rc == 198 test_pass "refuses: `bad'"
    else test_fail "refuses: `bad'" "rc=`=_rc'"
}

noi print_summary "creghdfe_fe_bootstrap"
}

global CTOOLS_COMPONENT_COMPLETE "creghdfe_fe_bootstrap"
    global CTOOLS_CASE_COMPLETE "creghdfe_fe_bootstrap"
    exit
}

if "`__case'" == "wild_draws" {
* The wild bootstrap's weights (build/_ctools_hdfe_boot.ado wilddraw,
* docs/development/BOOTSTRAP_PLAN.md, Phase 1w): for the same seed they are
* boottest's. Gate: a replicate-by-replicate WCR computed in Mata from the
* drawn weights reproduces boottest's bootstrap statistics (svmat, in
* replicate order) and its p-value, for every weight type, for enumeration
* and for the observation-level bootstrap.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running wild bootstrap draw checks..."

capture which boottest
if _rc {
    test_skip "wild draws against boottest" "boottest is not installed"
}
else {
    print_section "Weights against boottest's statistics"

    * Oracle ratios t*_b / t against boottest's (its small-sample factor
    * multiplies both); p by boottest's rule, ties within rounding not counted
    capture program drop _wild_vs_boottest
    program define _wild_vs_boottest, rclass
        syntax, k(integer) cl(varname)
        tempname dist
        matrix `dist' = r(dist)
        local t_bt = r(t)
        local p_bt = r(p)
        local reps_bt = r(reps)
        mata: _wv_X = st_data(., "tenure ttl_exp"), J(st_nobs(), 1, 1)
        mata: _wv_t = _bsref_wild_t(_wv_X, st_data(., "ln_wage"), st_data(., "`cl'"), `k', _wv_V)
        mata: _wv_d = st_matrix("`dist'")
        mata: st_local("maxdiff", strofreal(max(abs(_wv_t[|2 \ .|] :/ _wv_t[1] - _wv_d :/ `t_bt')), "%21.0g"))
        mata: st_local("p", strofreal(sum(abs(_wv_t[|2 \ .|]) :> abs(_wv_t[1]) * (1 + 1e-10)) / (rows(_wv_t) - 1), "%21.0g"))
        mata: mata drop _wv_X _wv_t _wv_d
        return scalar maxdiff = `maxdiff'
        return scalar p = `p'
        return scalar p_bt = `p_bt'
        return scalar reps_bt = `reps_bt'
    end

    webuse nlswork, clear
    keep if inrange(idcode, 1, 40) & !missing(ln_wage, tenure, ttl_exp)
    sort idcode year
    egen long cl = group(idcode)
    summarize cl, meanonly
    local C = r(max)
    regress ln_wage tenure ttl_exp, vce(cluster idcode)
    foreach type in rademacher mammen webb normal gamma {
        set seed 8
        boottest tenure, reps(19) seed(8) weighttype(`type') svmat nograph noci
        local state_bt = c(rngstate)
        set seed 8
        _ctools_hdfe_boot wilddraw, clusters(`C') reps(19) type(`type')
        local state_wd = c(rngstate)
        local enum = r(enumerate)
        mata: _wv_V = _ctools_wild_v
        boottest tenure, reps(19) seed(8) weighttype(`type') svmat nograph noci
        _wild_vs_boottest, k(1) cl(cl)
        if `enum' == 0 & r(maxdiff) < 1e-8 & abs(r(p) - r(p_bt)) < 1e-12 & "`state_wd'" == "`state_bt'" ///
            test_pass "`type': each replicate's statistic and the p-value equal boottest's (`C' clusters)"
        else test_fail "`type' weights" "max diff `r(maxdiff)'; p `r(p)' vs `r(p_bt)'; same state `=("`state_wd'" == "`state_bt'")'"
    }

    * Enumeration: no draws, boottest's order of the 2^C sign vectors
    keep if cl <= 8
    regress ln_wage tenure ttl_exp, vce(cluster idcode)
    set seed 9
    local state0 = c(rngstate)
    _ctools_hdfe_boot wilddraw, clusters(8) reps(999) type(rademacher)
    local enum = r(enumerate)
    local reps = r(reps)
    local same = (c(rngstate) == "`state0'")
    mata: _wv_V = _bsref_wild_enum(8)
    boottest ttl_exp, reps(999) seed(9) svmat nograph noci
    _wild_vs_boottest, k(2) cl(cl)
    if `enum' == 1 & `reps' == 256 & r(reps_bt) == 256 & `same' & r(maxdiff) < 1e-8 & abs(r(p) - r(p_bt)) < 1e-12 ///
        test_pass "enumeration: 256 sign vectors in boottest's order, no draw, same p-value"
    else test_fail "enumeration" "enum `enum' reps `reps'/`r(reps_bt)' same state `same' diff `r(maxdiff)' p `r(p)' vs `r(p_bt)'"

    * noenumerate draws instead
    set seed 9
    _ctools_hdfe_boot wilddraw, clusters(8) reps(999) type(rademacher) noenumerate
    local enum = r(enumerate)
    local reps = r(reps)
    mata: st_local("cols", strofreal(cols(_ctools_wild_v)))
    if `enum' == 0 & `reps' == 999 & `cols' == 999 & c(rngstate) != "`state0'" ///
        test_pass "noenumerate draws 999 replicates"
    else test_fail "noenumerate" "enumerate `enum' reps `reps' cols `cols'"

    * Observation-level bootstrap (vce(robust)): one weight per row, in data order
    webuse nlswork, clear
    keep if inrange(idcode, 1, 12) & !missing(ln_wage, tenure, ttl_exp)
    sort idcode year
    gen long row = _n
    regress ln_wage tenure ttl_exp, vce(robust)
    set seed 10
    _ctools_hdfe_boot wilddraw, clusters(`=_N') reps(29) type(webb)
    mata: _wv_V = _ctools_wild_v
    boottest tenure, reps(29) seed(10) weighttype(webb) svmat nograph noci
    _wild_vs_boottest, k(1) cl(row)
    if r(maxdiff) < 1e-8 & abs(r(p) - r(p_bt)) < 1e-12 ///
        test_pass "robust: one weight per observation, in data order"
    else test_fail "robust" "max diff `r(maxdiff)'; p `r(p)' vs `r(p_bt)'"
    mata: mata drop _wv_V
}

print_section "Helper subcommands"

clear
set obs 10
set seed 11
_ctools_hdfe_boot wilddraw, clusters(7) reps(12) type(normal)
_ctools_hdfe_boot wildweights, first(4) count(3) generate(wv)
local wvars `r(varlist)'
* packed replicate after replicate: 21 values in 3 variables of 10
* observations, replicate 5 spanning wv1 and wv2
mata: st_local("ok", strofreal(_bsref_packed("wv", 0, 21) == vec(_ctools_wild_v[., 4..6])))
capture confirm variable wv4
if `ok' == 1 & "`wvars'" == "wv1 wv2 wv3" & _rc ///
    test_pass "wildweights: replicates 4-6 packed in order into wv1-wv3"
else test_fail "wildweights" "ok `ok' vars `wvars'"

* Five standard errors, derived from sample size and each distribution's
* fourth central moment, replace a fixed 1% allowance. For Rademacher the
* variance is determined by the mean and the exact {-1, 1} support.
foreach type in mammen webb gamma rademacher normal {
    set seed 13
    _ctools_hdfe_boot wilddraw, clusters(5000) reps(400) type(`type') noenumerate
    mata: st_local("m", strofreal(mean(vec(_ctools_wild_v)), "%24.17g"))
    mata: st_local("v", strofreal(variance(vec(_ctools_wild_v)), "%24.17g"))
    mata: st_local("n", strofreal(length(_ctools_wild_v)))
    local mu4 = cond("`type'" == "mammen", 2, cond("`type'" == "webb", 7/6, cond("`type'" == "gamma", 4.5, 3)))
    local var_se = sqrt((`mu4' - (`n'-3)/(`n'-1))/`n')
    local var_ok = abs(`v'-1) <= 5*`var_se'
    local support_ok = 1
    if "`type'" == "mammen" {
        local roundoff = 16*c(epsdouble)
        mata: st_local("support_ok", strofreal(all(abs(_ctools_wild_v :- (1-sqrt(5))/2) :< `roundoff' :| abs(_ctools_wild_v :- (1+sqrt(5))/2) :< `roundoff')))
    }
    if "`type'" == "rademacher" {
        mata: st_local("support_ok", strofreal(all(abs(_ctools_wild_v) :== 1)))
        * |mean| <= 5/sqrt(n) implies |s^2-1| <= (1+5^2)/(n-1).
        local var_ok = abs(`v'-1) <= 26/(`n'-1)
    }
    if `n' == 5000*400 & abs(`m') <= 5/sqrt(`n') & `var_ok' & `support_ok' ///
        test_pass "`type': support and sampling-error bounds for moments"
    else test_fail "`type' moments" "mean `m' variance `v'"
}

_ctools_hdfe_boot wildclear
capture mata: mata describe _ctools_wild_v
local cleared = _rc
capture _ctools_hdfe_boot wildweights, first(1) count(1) generate(zz)
if `cleared' & _rc == 198 test_pass "wildclear drops the draws; wildweights then refuses"
else test_fail "wildclear" "describe rc `cleared'; wildweights rc `=_rc'"

capture _ctools_hdfe_boot wilddraw, clusters(5) reps(10) type(uniform)
local rc1 = _rc
set seed 14
_ctools_hdfe_boot wilddraw, clusters(5) reps(10) type(normal)
capture _ctools_hdfe_boot wildweights, first(9) count(3) generate(zz)
if `rc1' == 198 & _rc == 198 test_pass "wilddraw and wildweights refuse bad input"
else test_fail "wild input errors" "type rc `rc1'; range rc `=_rc'"
_ctools_hdfe_boot wildclear

noi print_summary "wild_draws"
}

global CTOOLS_COMPONENT_COMPLETE "wild_draws"
    global CTOOLS_CASE_COMPLETE "wild_draws"
    exit
}

if "`__case'" == "creghdfe_wild" {
* creghdfe, wildbootstrap() (docs/development/BOOTSTRAP_PLAN.md, Phase 2d).
* Reference: boottest 4.5.5, with the same seed. p-values must be equal;
* interval bounds are checked through boottest's own p-value just inside and
* just outside each bound (boottest's bounds stop within one replicate step
* of the level, so they are not compared directly).

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running creghdfe wildbootstrap() checks..."

* boottest's p-value for H0: coef = value, after the current estimates
capture program drop _wb_bt
program define _wb_bt, rclass
    syntax, coef(string) [value(real 0) opts(string)]
    capture noisily boottest (`coef' = `value'), nograph noci `opts'
    if _rc {
        return scalar p = .
        return scalar reps = .
        exit
    }
    return scalar p = r(p)
    return scalar reps = r(reps)
end

* The data: 60 women of nlswork, no singleton panels
capture program drop _wb_data
program define _wb_data
    webuse nlswork, clear
    keep if inrange(idcode, 1, 60) & !missing(ln_wage, tenure, ttl_exp, year)
    bysort idcode: keep if _N >= 2
    sort idcode year
    gen w = 0.5 + mod(idcode * 7 + year, 11) / 10
    gen grp = ceil(idcode / 4)
end

capture which boottest
local have_bt = !_rc
capture which reghdfe
local have_bt = `have_bt' & !_rc
if !`have_bt' {
    test_skip "wildbootstrap() against boottest" "boottest or reghdfe is not installed"
}
else {
    print_section "p-values and intervals against boottest"

    * designs: label | creghdfe and reghdfe options | boottest options | ptype
    local d1 "nested|absorb(idcode) vce(cluster idcode)|weighttype(rademacher)"
    local d2 "crossed (projected)|absorb(year) vce(cluster idcode)|weighttype(webb)"
    local d3 "aweights|absorb(idcode) vce(cluster idcode) [aw=w]|weighttype(mammen)"
    local d4 "unrestricted|absorb(year) vce(cluster idcode)|weighttype(rademacher) nonull"
    local d5 "robust|absorb(year) vce(robust)|weighttype(normal)"
    local d6 "equal-tailed|absorb(year) vce(cluster idcode)|weighttype(gamma) ptype(equaltail)"
    forvalues d = 1/6 {
        gettoken label rest : d`d', parse("|")
        gettoken bar rest : rest, parse("|")
        gettoken model rest : rest, parse("|")
        gettoken bar btopts : rest, parse("|")
        local wexp
        if strpos("`model'", "[aw=w]") {
            local wexp [aw=w]
            local model = subinstr("`model'", "[aw=w]", "", 1)
        }
        _wb_data
        creghdfe ln_wage tenure ttl_exp `wexp', `model' wildbootstrap(reps(199) seed(31) `btopts' nodots)
        matrix P = e(wb_p)
        matrix CI = e(wb_ci)
        matrix t = e(wb_t)
        local path `e(wb_fe_path)'
        local nonull = strpos("`btopts'", "nonull") > 0
        reghdfe ln_wage tenure ttl_exp `wexp', `model'
        local bad
        local j 0
        foreach c in tenure ttl_exp {
            local ++j
            _wb_bt, coef(`c') opts(reps(199) seed(31) `btopts')
            if abs(r(p) - P[1, `j']) > 1e-12 local bad `bad' p(`c')=`=P[1,`j']'/`r(p)'
            * the interval: boottest's p just inside and just outside each bound
            if !`nonull' {
                local se = _se[`c']
                foreach side in 1 2 {
                    local x = CI[`j', `side']
                    local dir = cond(`side' == 1, 1, -1)
                    _wb_bt, coef(`c') value(`=string(`x' + `dir' * 1e-5 * `se', "%21.17g")') opts(reps(199) seed(31) `btopts')
                    local pin = r(p)
                    _wb_bt, coef(`c') value(`=string(`x' - `dir' * 1e-5 * `se', "%21.17g")') opts(reps(199) seed(31) `btopts')
                    local pout = r(p)
                    if !(`pin' > 0.05 & `pout' <= 0.05) local bad `bad' ci(`c',`side')=`pin'/`pout'
                }
            }
        }
        if "`bad'" == "" test_pass "`label' (`path'): p-values and interval bounds agree with boottest"
        else test_fail "`label'" "`bad'"
    }

    * Coarser bootstrap clusters (bootcluster()): boottest 4.5.5 fails with
    * r(3301) in crosstabCapstarMinus whenever bootcluster() is coarser than
    * the clusters, so the oracle is the replicate-by-replicate WCR in Mata
    * (the effect as dummies) on the weights drawn for the same seed
    _wb_data
    creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode) ///
        wildbootstrap(reps(199) seed(39) bootcluster(grp) nodots)
    matrix P = e(wb_p)
    local C = e(wb_N_clust)
    egen long cl = group(idcode)
    egen long bc = group(grp)
    quietly tabulate idcode, generate(D_)
    set seed 39
    _ctools_hdfe_boot wilddraw, clusters(`C') reps(199) type(rademacher)
    mata: _wb_X = st_data(., "tenure ttl_exp"), st_data(., "D_*")
    local bad
    forvalues k = 1/2 {
        mata: _wb_t = _bsref_wild_t2(_wb_X, st_data(., "ln_wage"), st_data(., "cl"), st_data(., "bc"), `k', _ctools_wild_v)
        mata: st_local("p", strofreal(sum(abs(_wb_t[|2 \ .|]) :> abs(_wb_t[1]) * (1 + 1e-10)) / 199, "%21.17g"))
        if abs(`p' - P[1, `k']) > 1e-12 local bad `bad' p`k'=`=P[1,`k']'/`p'
    }
    mata: mata drop _wb_X _wb_t
    _ctools_hdfe_boot wildclear
    if "`bad'" == "" & `C' == 15 test_pass "bootcluster(): 15 coarser clusters, p-values equal the Mata oracle's"
    else test_fail "bootcluster" "`bad'; clusters `C'"

    * Two absorbed effects: boottest after areg with the second as dummies
    _wb_data
    creghdfe ln_wage tenure ttl_exp, absorb(idcode year) vce(cluster idcode) ///
        wildbootstrap(reps(199) seed(32) weighttype(webb) nodots)
    matrix P = e(wb_p)
    local path `e(wb_fe_path)'
    areg ln_wage tenure ttl_exp i.year, absorb(idcode) vce(cluster idcode)
    _wb_bt, coef(tenure) opts(reps(199) seed(32) weighttype(webb))
    local p1 = r(p)
    _wb_bt, coef(ttl_exp) opts(reps(199) seed(32) weighttype(webb))
    if abs(`p1' - P[1, 1]) < 1e-12 & abs(r(p) - P[1, 2]) < 1e-12 & "`path'" == "projected" ///
        test_pass "two effects (projected): p-values equal boottest after areg with year dummies"
    else test_fail "two effects" "p `=P[1,1]' `=P[1,2]' vs `p1' `r(p)'; path `path'"

    * fedgp(partialled): the bootstrap of the demeaned regression
    _wb_data
    creghdfe ln_wage tenure ttl_exp, absorb(year) vce(cluster idcode) ///
        wildbootstrap(reps(199) seed(33) fedgp(partialled) nodots)
    matrix P = e(wb_p)
    local path `e(wb_fe_path)'
    foreach v in ln_wage tenure ttl_exp {
        egen double m_`v' = mean(`v'), by(year)
        gen double d_`v' = `v' - m_`v'
    }
    regress d_ln_wage d_tenure d_ttl_exp, nocons vce(cluster idcode)
    _wb_bt, coef(d_tenure) opts(reps(199) seed(33))
    local p1 = r(p)
    _wb_bt, coef(d_ttl_exp) opts(reps(199) seed(33))
    if abs(`p1' - P[1, 1]) < 1e-12 & abs(r(p) - P[1, 2]) < 1e-12 & "`path'" == "partialled" ///
        test_pass "fedgp(partialled): p-values equal boottest on the demeaned regression"
    else test_fail "fedgp(partialled)" "p `=P[1,1]' `=P[1,2]' vs `p1' `r(p)'; path `path'"

    * Enumeration with 8 clusters: 256 sign vectors, no draw
    _wb_data
    keep if idcode <= 9
    set seed 34
    local state0 = c(rngstate)
    creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode) wildbootstrap(nodots)
    local same = (c(rngstate) == "`state0'")
    matrix P = e(wb_p)
    local reps = e(wb_reps)
    local enum = e(wb_enumerated)
    local C = e(wb_N_clust)
    reghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode)
    _wb_bt, coef(tenure) opts(reps(999))
    if `enum' == 1 & `reps' == 2^`C' & r(reps) == `reps' & `same' & abs(r(p) - P[1, 1]) < 1e-12 ///
        test_pass "enumeration: all 2^`C' sign vectors, RNG state unchanged, p equals boottest's"
    else test_fail "enumeration" "enum `enum' reps `reps'/`r(reps)' same `same' p `=P[1,1]' vs `r(p)'"
}

print_section "Paths, outputs and errors"

* The projected path on a nested design gives the aggregate path's results
_wb_data
creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode) wildbootstrap(reps(99) seed(35) nodots)
matrix Pa = e(wb_p)
matrix Ca = e(wb_ci)
local pa `e(wb_fe_path)'
creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode) wildbootstrap(reps(99) seed(35) nodots _projected)
matrix Pp = e(wb_p)
matrix Cp = e(wb_ci)
local pp `e(wb_fe_path)'
if mreldif(Pa, Pp) == 0 & mreldif(Ca, Cp) < 1e-8 & "`pa'" == "nested" & "`pp'" == "projected" ///
    test_pass "projected and aggregate paths agree on a nested design"
else test_fail "path agreement" "p `=mreldif(Pa, Pp)' ci `=mreldif(Ca, Cp)' paths `pa' `pp'"

* e(b), e(V) and the analytic table are unchanged; replay and nowildtable
creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode)
matrix b0 = e(b)
matrix V0 = e(V)
creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode) ///
    wildbootstrap(reps(49) seed(36) coefs(ttl_exp) noci nodots) nowildtable
local ok = mreldif(b0, e(b)) == 0 & mreldif(V0, e(V)) == 0 & "`e(wb_coefs)'" == "ttl_exp" & ///
    colsof(e(wb_p)) == 1 & "`e(wb_method)'" == "wcr" & e(wb_reps) == 49
capture confirm matrix e(wb_ci)
local noci = _rc != 0
capture noisily creghdfe
if `ok' & `noci' & _rc == 0 test_pass "e(b), e(V) unchanged; coefs(), noci, nowildtable, replay"
else test_fail "outputs" "ok `ok' noci `noci' replay rc `=_rc'"

* With fevce(bootstrap): both run on one fit, each as alone
creghdfe ln_wage tenure ttl_exp, absorb(idcode, savefe) vce(cluster idcode) ///
    fevce(bootstrap, reps(20) seed(37) nodots)
rename __fese1__ se_alone
drop __hdfe1__
creghdfe ln_wage tenure ttl_exp, absorb(idcode) vce(cluster idcode) wildbootstrap(reps(99) seed(38) nodots)
matrix Pw = e(wb_p)
creghdfe ln_wage tenure ttl_exp, absorb(idcode, savefe) vce(cluster idcode) ///
    wildbootstrap(reps(99) seed(38) nodots) fevce(bootstrap, reps(20) seed(37) nodots)
count if reldif(__fese1__, se_alone) > 1e-12
if r(N) == 0 & mreldif(Pw, e(wb_p)) == 0 & "`e(fevce)'" == "bootstrap" ///
    test_pass "wildbootstrap() with fevce(bootstrap): each as when run alone"
else test_fail "wild with fevce" "SE rows differ `r(N)'; p `=mreldif(Pw, e(wb_p))'"

* Refusals (plan 4.6)
gen fw = 1 + mod(_n, 3)
local e1 "absorb(idcode) wildbootstrap(seed(1))"
local e2 "absorb(idcode) vce(bootstrap, reps(10) seed(1)) wildbootstrap(seed(1))"
local e3 "absorb(year) vce(robust) wildbootstrap(seed(1)) [fw=fw]"
local e4 "absorb(idcode) vce(cluster idcode) wildbootstrap(seed(1) coefs(hours))"
local e5 "absorb(idcode) vce(cluster idcode) wildbootstrap(seed(1) weighttype(uniform))"
local e6 "absorb(idcode) vce(cluster idcode) wildbootstrap(seed(1) bootcluster(grp year))"
local e7 "absorb(idcode) vce(cluster idcode) wildbootstrap(seed(1) fedgp(other))"
forvalues e = 1/7 {
    local opts `e`e''
    local wexp
    if strpos("`opts'", "[fw=fw]") {
        local wexp [fw=fw]
        local opts = subinstr("`opts'", "[fw=fw]", "", 1)
    }
    capture creghdfe ln_wage tenure ttl_exp `wexp', `opts'
    if _rc == 198 test_pass "refuses: `e`e''"
    else test_fail "refuses: `e`e''" "rc=`=_rc'"
}

noi print_summary "creghdfe_wild"
}

global CTOOLS_COMPONENT_COMPLETE "creghdfe_wild"
    global CTOOLS_CASE_COMPLETE "creghdfe_wild"
    exit
}

if "`__case'" == "cppmlhdfe_bootstrap" {
* cppmlhdfe, vce(bootstrap) (docs/development/BOOTSTRAP_PLAN.md, Phase 3).
* Oracle: the bootstrap prefix on cppmlhdfe itself with the same seed
* (draws(stata) reproduces bsample's draws), saving the replicates in
* double. The prefix refits each replicate from scratch; the plugin
* warm-starts from the full fit, so the two agree to the IRLS tolerance.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running cppmlhdfe vce(bootstrap) checks..."

capture program drop _pb_data
program define _pb_data
    clear
    set seed 77
    set obs 720
    gen id = ceil(_n / 6)
    gen t = mod(_n - 1, 6) + 1
    gen cl = ceil(id / 4)
    gen str = mod(cl, 3)
    gen x1 = rnormal()
    gen x2 = rnormal() + id / 100
    gen off = rnormal() / 5
    gen double y = rpoisson(exp(0.3 * x1 - 0.2 * x2 + mod(id, 7) / 5 + t / 10 + off))
    * ids with a single positive outcome: replicates that miss it lose the id
    replace y = 0 if mod(id, 17) == 0 & t > 1
    * a dummy present in few rows: some replicates omit it
    gen rare = inlist(_n, 3, 200, 450)
    * a regressor positive only where y = 0 in a few rows (ReLU separation
    * in replicates that draw those rows without the others)
    gen sepx = (y == 0 & mod(_n, 97) == 0)
end

print_section "Parity with the bootstrap prefix (draws(stata))"

local d1 "observations|y x1 x2, absorb(id t)|"
local d2 "clusters|y x1 x2, absorb(id t)|cluster(cl)"
local d3 "strata and clusters|y x1 x2, absorb(id t)|strata(str) cluster(cl)"
local d4 "offset, rare dummy|y x1 x2 rare, absorb(id t) offset(off)|cluster(cl)"
local d5 "separating regressor|y x1 sepx, absorb(id t)|"
forvalues d = 1/5 {
    gettoken label rest : d`d', parse("|")
    gettoken bar rest : rest, parse("|")
    gettoken model rest : rest, parse("|")
    gettoken bar design : rest, parse("|")
    local design = strtrim(`"`design'"')
    _pb_data
    cppmlhdfe `model' vce(bootstrap, reps(25) seed(41) draws(stata) `design' nodots)
    matrix V1 = e(V)
    matrix b1 = e(b)
    local reps1 = e(N_reps)
    local mis1 = e(N_misreps)
    local state1 = c(rngstate)
    tempfile bsr
    local popts = cond(`"`design'"' == "", "", `"`design'"')
    bootstrap _b, reps(25) seed(41) `popts' nodots saving(`bsr', double): cppmlhdfe `model' nolog
    matrix V2 = e(V)
    local reps2 = e(N_reps)
    local mis2 = e(N_misreps)
    local state2 = c(rngstate)
    if mreldif(V1, V2) < 1e-8 & `reps1' == `reps2' & `mis1' == `mis2' & "`state1'" == "`state2'" ///
        test_pass "`label': e(V), N_reps, N_misreps and RNG state equal the prefix's"
    else test_fail "`label'" "mreldif `=mreldif(V1, V2)'; reps `reps1'/`reps2'; misreps `mis1'/`mis2'; state `=("`state1'"=="`state2'")'"
}

* An independent estimator: the prefix on ppmlhdfe (its own IRLS, cold
* starts and convergence rule), equal to its tolerance
capture which ppmlhdfe
if _rc test_skip "parity with ppmlhdfe" "ppmlhdfe is not installed"
else {
    _pb_data
    cppmlhdfe y x1 x2, absorb(id t) vce(bootstrap, reps(25) seed(44) draws(stata) cluster(cl) nodots)
    matrix V1 = e(V)
    tempfile bsr
    bootstrap _b, reps(25) seed(44) cluster(cl) nodots saving(`bsr', double): ppmlhdfe y x1 x2, absorb(id t)
    matrix V2 = e(V)
    if mreldif(V1, V2) < 1e-6 test_pass "clusters: e(V) equals the prefix on ppmlhdfe (mreldif `=string(mreldif(V1, V2), "%9.2e")')"
    else test_fail "ppmlhdfe parity" "mreldif `=mreldif(V1, V2)'"
}

print_section "Outputs, modes and errors"

* draws(fast): reproducible from the seed, independent of threads
_pb_data
cppmlhdfe y x1 x2, absorb(id t) vce(bootstrap, reps(30) seed(42) draws(fast) nodots)
matrix Va = e(V)
cppmlhdfe y x1 x2, absorb(id t) vce(bootstrap, reps(30) seed(42) draws(fast) nodots) threads(1)
matrix Vb = e(V)
if mreldif(Va, Vb) == 0 test_pass "draws(fast): same seed, same e(V), with any thread count"
else test_fail "draws(fast)" "mreldif `=mreldif(Va, Vb)'"

* e(): what vce(bootstrap) posts; e(b) is the fit's
cppmlhdfe y x1 x2, absorb(id t)
matrix b0 = e(b)
matrix V0 = e(V)
cppmlhdfe y x1 x2, absorb(id t) vce(bootstrap, reps(30) seed(43) nodots)
local ok = mreldif(b0, e(b)) == 0 & mreldif(V0, e(V_modelbased)) == 0 & "`e(vce)'" == "bootstrap" & ///
    "`e(vcetype)'" == "Bootstrap" & e(N_reps) == 30 & "`e(prefix)'" == "bootstrap" & "`e(rngstate)'" != ""
capture noisily estat bootstrap, percentile
local rc1 = _rc
capture noisily cppmlhdfe
if `ok' & `rc1' == 0 & _rc == 0 test_pass "e(b) the fit's, e(V_modelbased) the robust V; estat bootstrap and replay"
else test_fail "posted results" "ok `ok' estat `rc1' replay `=_rc'"

* mse: the variance about e(b)
cppmlhdfe y x1 x2, absorb(id t) vce(bootstrap, reps(30) seed(43) nodots)
matrix Vd = e(V)
cppmlhdfe y x1 x2, absorb(id t) vce(bootstrap, reps(30) seed(43) mse nodots)
matrix Vm = e(V)
matrix bias = e(bias)
local x = Vm[1, 1] - Vd[1, 1] * 29 / 30 - bias[1, 1]^2
if abs(`x') < 1e-12 * Vm[1, 1] test_pass "mse: V = (R-1)/R V + bias^2"
else test_fail "mse" "difference `x'"

* fweights with clusters: the bootstrap of the expanded data, whose clusters
* are drawn with all their copies (the prefix, same seed and draws). Each
* drawn copy of a row counts once for the singleton rule, as cppmlhdfe and
* ppmlhdfe count fweighted rows.
_pb_data
gen fw = 1 + mod(_n, 3)
cppmlhdfe y x1 x2 [fw=fw], absorb(id t) vce(bootstrap, reps(25) seed(45) draws(stata) cluster(cl) nodots)
matrix V1 = e(V)
expand fw
tempfile bsr
bootstrap _b, reps(25) seed(45) cluster(cl) nodots saving(`bsr', double): cppmlhdfe y x1 x2, absorb(id t) nolog
matrix V2 = e(V)
if mreldif(V1, V2) < 1e-8 test_pass "fweighted clusters: e(V) equals the prefix on the expanded data"
else test_fail "fweighted clusters" "mreldif `=mreldif(V1, V2)'"

* fweighted rows: the observations they stand for are drawn, and a row
* weighs what was drawn of it (drawing each row as one unit weighted by its
* fweight gave an SE 3.7 times the expanded data's here)
sysuse auto, clear
gen fw = cond(mod(_n, 4) == 0, 30, 1)
cppmlhdfe price mpg weight [fw=fw], absorb(foreign) vce(bootstrap, reps(300) seed(1) nodots)
local bias = el(e(bias), 1, 1)
local se = _se[mpg]
expand fw
cppmlhdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(300) seed(1) nodots)
if abs(`bias') < 4 * `se' / sqrt(300) & abs(`se' / _se[mpg] - 1) < 0.2 ///
    test_pass "fweighted rows: replicates centered on e(b), SEs as the expanded data's"
else test_fail "fweighted rows" "bias `bias' se `se' vs expanded `=_se[mpg]'"

* Refusals
foreach bad in "absorb(id t) vce(bootstrap, seed(1)) tagsep(ts)" ///
    "absorb(id t) vce(bootstrap, seed(1) cluster(cl t))" {
    _pb_data
    capture cppmlhdfe y x1 x2, `bad'
    if _rc == 198 test_pass "refuses: `bad'"
    else test_fail "refuses: `bad'" "rc=`=_rc'"
}
_pb_data
capture cppmlhdfe y, absorb(id t) vce(bootstrap, seed(1))
if _rc == 198 test_pass "refuses: no regressor"
else test_fail "refuses: no regressor" "rc=`=_rc'"

noi print_summary "cppmlhdfe_bootstrap"
}

global CTOOLS_COMPONENT_COMPLETE "cppmlhdfe_bootstrap"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_bootstrap"
    exit
}

if "`__case'" == "cppmlhdfe_fe_bootstrap" {
* cppmlhdfe, fevce(bootstrap) (docs/development/BOOTSTRAP_PLAN.md, Phase 3b).
* Oracle: the same draws (draws(stata), bsref_draws); in each replicate the
* rows are expanded by their multiplicities and cppmlhdfe refits with
* absorb(f1=... f2=...); each level's effect goes back to every row of it,
* and the standard deviation over the replicates that estimate it is
* compared, as is that of their sum d. The plugin warm-starts its IRLS, so
* the two agree to the solver tolerances.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running cppmlhdfe fevce(bootstrap) checks..."

capture program drop _pf_data
program define _pf_data
    clear
    set seed 78
    set obs 480
    gen id = ceil(_n / 6)
    gen t = mod(_n - 1, 6) + 1
    gen cl = ceil(id / 4)
    gen x1 = rnormal()
    gen double y = rpoisson(exp(0.3 * x1 + mod(id, 7) / 5 + t / 10))
end

capture program drop _pf_oracle
program define _pf_oracle
    * o_se1.., o_n1.., o_dse on e(sample) rows, m1..mR the draws
    syntax, reps(integer) absorb(string) model(string)
    local G : word count `absorb'
    forvalues r = 1/`reps' {
        local fes
        forvalues g = 1/`G' {
            local fes `fes' f`g'_`r'=`: word `g' of `absorb''
        }
        tempfile full expanded
        forvalues g = 1/`G' {
            tempfile lv`g'
        }
        quietly save `full'
        quietly keep if m`r' > 0
        quietly expand m`r'
        quietly cppmlhdfe `model', absorb(`fes') nolog
        quietly keep if e(sample)
        quietly save `expanded'
        forvalues g = 1/`G' {
            local v : word `g' of `absorb'
            use `v' f`g'_`r' using `expanded', clear
            quietly collapse (firstnm) f`g'_`r', by(`v')
            quietly save `lv`g''
        }
        use `full', clear
        local dsum
        forvalues g = 1/`G' {
            local v : word `g' of `absorb'
            quietly merge m:1 `v' using `lv`g'', nogenerate
            local dsum `dsum' + f`g'_`r'
        }
        quietly generate double d_`r' = 0 `dsum'
    }
    forvalues g = 1/`G' {
        quietly egen double o_se`g' = rowsd(f`g'_*)
        quietly egen double o_n`g' = rownonmiss(f`g'_*)
    }
    quietly egen double o_dse = rowsd(d_*)
end

print_section "Fixed-effect standard errors against the oracle"

* One effect, cluster bootstrap
_pf_data
cppmlhdfe y x1, absorb(fid=id) fevce(bootstrap, reps(20) seed(51) cluster(cl) draws(stata) nreps dse(dse) nodots)
local reps = e(fe_N_reps)
local warn = e(fe_norm_warning)
keep if e(sample)
keep y x1 id t cl fid fid_se fid_nrep dse
set seed 51
bsref_draws, reps(20) generate(m) cluster(cl)
_pf_oracle, reps(20) absorb(id) model(y x1)
count if reldif(fid_se, o_se1) > 1e-5 | fid_nrep != o_n1 | reldif(dse, o_dse) > 1e-5
if r(N) == 0 & `reps' == 20 & `warn' == 0 ///
    test_pass "one effect: SE, valid replicates and d's SE equal the oracle; no warning"
else test_fail "one effect" "`r(N)' rows differ; warning `warn'"

* Two effects, observation bootstrap, savefe and d()
_pf_data
cppmlhdfe y x1, absorb(id t, savefe) d(dd) fevce(bootstrap, reps(20) seed(52) draws(stata) nreps nodots)
local warn = e(fe_norm_warning)
local splits = e(fe_split_reps)
keep if e(sample)
keep y x1 id t __fese1__ __fese2__ __fenrep1__ __fenrep2__ dd_se
set seed 52
bsref_draws, reps(20) generate(m)
_pf_oracle, reps(20) absorb(id t) model(y x1)
count if reldif(__fese1__, o_se1) > 1e-5 | reldif(__fese2__, o_se2) > 1e-5 | ///
    __fenrep1__ != o_n1 | __fenrep2__ != o_n2 | reldif(dd_se, o_dse) > 1e-5
if r(N) == 0 & `warn' == 1 & `splits' == 0 ///
    test_pass "two effects: per-level and d SEs equal the oracle; normalization warning set"
else test_fail "two effects" "`r(N)' rows differ; warning `warn'; splits `splits'"

print_section "Modes, outputs and errors"

* Both targets on one set of replicates, each as alone
_pf_data
cppmlhdfe y x1, absorb(fid=id t) vce(bootstrap, reps(20) seed(53) nodots) fevce(bootstrap)
matrix Vboth = e(V)
rename fid_se se_both
drop fid
cppmlhdfe y x1, absorb(id t) vce(bootstrap, reps(20) seed(53) nodots)
matrix Vcoef = e(V)
cppmlhdfe y x1, absorb(fid=id t) fevce(bootstrap, reps(20) seed(53) nodots)
count if reldif(se_both, fid_se) > 1e-12
if mreldif(Vboth, Vcoef) < 1e-12 & r(N) == 0 & "`e(vce)'" == "robust" ///
    test_pass "vce() and fevce() share replicates; fevce alone keeps the robust e(V)"
else test_fail "shared replicates" "mreldif `=mreldif(Vboth, Vcoef)'; rows `r(N)'; vce `e(vce)'"

* Names: automatic outputs are replaced on a rerun, sestub(), the posted
* random-number states
_pf_data
cppmlhdfe y x1, absorb(id, savefe) fevce(bootstrap, reps(10) seed(54) nodots)
capture noisily cppmlhdfe y x1, absorb(id, savefe) fevce(bootstrap, reps(10) seed(54) nodots)
local rerun = _rc
cppmlhdfe y x1, absorb(fx=id) fevce(bootstrap, reps(10) seed(54) sestub(sx) nreps nodots)
capture confirm variable sx1 sx1_nrep
local stub = _rc
set seed 54
local s54 = c(rngstate)
local ok = "`e(fe_rngstate)'" == "`s54'"
cppmlhdfe y x1, absorb(id) vce(bootstrap, reps(10) seed(54) nodots)
local ok = `ok' & "`e(rngstate)'" == "`s54'"
creghdfe y x1, absorb(id) vce(bootstrap, reps(10) seed(54) nodots)
local ok = `ok' & "`e(rngstate)'" == "`s54'"
if `rerun' == 0 & `stub' == 0 & `ok' ///
    test_pass "automatic outputs replaced on rerun; sestub(); e(rngstate) is the seeded state"
else test_fail "names and states" "rerun `rerun' stub `stub' states `ok'"

* Refusals
foreach bad in "absorb(id) fevce(bootstrap, seed(1))" ///
    "absorb(id, savefe) fevce(robust)" ///
    "absorb(id, savefe) vce(bootstrap, reps(10) seed(1)) fevce(bootstrap, reps(20))" ///
    "absorb(id) d(dd) fevce(bootstrap, ci seed(1))" {
    _pf_data
    capture cppmlhdfe y x1, `bad'
    if _rc == 198 test_pass "refuses: `bad'"
    else test_fail "refuses: `bad'" "rc=`=_rc'"
}

noi print_summary "cppmlhdfe_fe_bootstrap"
}

global CTOOLS_COMPONENT_COMPLETE "cppmlhdfe_fe_bootstrap"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_fe_bootstrap"
    exit
}

if "`__case'" == "cppmlhdfe_wild" {
* cppmlhdfe, wildbootstrap(): the score bootstrap (docs/development/
* BOOTSTRAP_PLAN.md, Phase 3w). Reference: boottest 4.5.5 after poisson
* with the effects as dummies, same seed: p-values must be equal. boottest
* takes the ML scores unweighted, which with fweights disagrees with its
* own result on the expanded data; cppmlhdfe weights them, and is compared
* with boottest on the expanded data. boottest computes no intervals after
* ML; the bounds are checked through its p-value just inside and outside.

do "validation/support/validate_setup.do"

quietly {

noi di as text "Running cppmlhdfe wildbootstrap() checks..."

capture program drop _pw_bt
program define _pw_bt, rclass
    syntax, coef(string) [value(real 0) opts(string)]
    capture noisily boottest (`coef' = `value'), nograph noci `opts'
    return scalar p = cond(_rc, ., r(p))
    return scalar reps = cond(_rc, ., r(reps))
end

capture program drop _pw_data
program define _pw_data
    clear
    set seed 82
    set obs 360
    gen id = ceil(_n / 6)
    gen t = mod(_n - 1, 6) + 1
    gen cl = ceil(id / 2)
    gen x1 = rnormal()
    gen x2 = rnormal() + id / 60
    gen off = rnormal() / 5
    gen fw = 1 + mod(_n, 3)
    gen obs = _n
    gen double y = rpoisson(exp(0.08 * x1 - 0.04 * x2 + mod(id, 5) / 4 + t / 10 + off))
end

capture which boottest
if _rc {
    test_skip "wildbootstrap() against boottest" "boottest is not installed"
}
else {
    print_section "p-values against boottest after poisson"

    * label | cppmlhdfe model | poisson model | options for both
    local d1 "clusters|absorb(id) vce(cluster cl)|i.id, vce(cluster cl)|weighttype(rademacher)"
    local d2 "two effects, Webb|absorb(id t) vce(cluster cl)|i.id i.t, vce(cluster cl)|weighttype(webb)"
    local d3 "offset, normal|absorb(id) vce(cluster cl) offset(off)|i.id, vce(cluster cl) offset(off)|weighttype(normal)"
    * (boottest 4.5.5 fails with r(3200) after poisson, vce(robust); a
    * cluster per observation is the same bootstrap)
    local d4 "observation clusters, Mammen|absorb(t) vce(cluster obs)|i.t, vce(cluster obs)|weighttype(mammen)"
    local d5 "equal-tailed, gamma|absorb(id) vce(cluster cl)|i.id, vce(cluster cl)|weighttype(gamma) ptype(equaltail)"
    forvalues d = 1/5 {
        gettoken label rest : d`d', parse("|")
        gettoken bar rest : rest, parse("|")
        gettoken cmodel rest : rest, parse("|")
        gettoken bar rest : rest, parse("|")
        gettoken pmodel rest : rest, parse("|")
        gettoken bar opts : rest, parse("|")
        _pw_data
        cppmlhdfe y x1 x2, `cmodel' wildbootstrap(reps(199) seed(61) `opts' nodots)
        matrix P = e(wb_p)
        poisson y x1 x2 `pmodel'
        _pw_bt, coef(x1) opts(reps(199) seed(61) `opts')
        local p1 = r(p)
        _pw_bt, coef(x2) opts(reps(199) seed(61) `opts')
        if abs(`p1' - P[1, 1]) < 1e-12 & abs(r(p) - P[1, 2]) < 1e-12 ///
            test_pass "`label': p-values equal boottest's"
        else test_fail "`label'" "p `=P[1,1]' `=P[1,2]' vs `p1' `r(p)'"
    }

    * vce(robust) is the bootstrap with a cluster per observation
    _pw_data
    cppmlhdfe y x1 x2, absorb(t) vce(robust) wildbootstrap(reps(199) seed(61) weighttype(mammen) nodots)
    matrix Pr = e(wb_p)
    cppmlhdfe y x1 x2, absorb(t) vce(cluster obs) wildbootstrap(reps(199) seed(61) weighttype(mammen) nodots)
    if mreldif(Pr, e(wb_p)) == 0 & e(wb_N_clust) == _N test_pass "vce(robust): one weight per observation, as vce(cluster obs)"
    else test_fail "robust" "mreldif `=mreldif(Pr, e(wb_p))'"

    * fweights: boottest on the expanded data
    _pw_data
    cppmlhdfe y x1 x2 [fw=fw], absorb(id) vce(cluster cl) wildbootstrap(reps(199) seed(62) nodots)
    matrix P = e(wb_p)
    expand fw
    poisson y x1 x2 i.id, vce(cluster cl)
    _pw_bt, coef(x1) opts(reps(199) seed(62))
    local p1 = r(p)
    _pw_bt, coef(x2) opts(reps(199) seed(62))
    if abs(`p1' - P[1, 1]) < 1e-12 & abs(r(p) - P[1, 2]) < 1e-12 ///
        test_pass "fweights: p-values equal boottest's on the expanded data"
    else test_fail "fweights" "p `=P[1,1]' `=P[1,2]' vs `p1' `r(p)'"

    * Enumeration: 8 clusters, all 256 sign vectors, no draw
    _pw_data
    keep if cl <= 8
    set seed 63
    local state0 = c(rngstate)
    cppmlhdfe y x1 x2, absorb(id) vce(cluster cl) wildbootstrap(nodots)
    local same = (c(rngstate) == "`state0'")
    matrix P = e(wb_p)
    local reps = e(wb_reps)
    local enum = e(wb_enumerated)
    poisson y x1 x2 i.id, vce(cluster cl)
    _pw_bt, coef(x1) opts(reps(999))
    if `enum' == 1 & `reps' == 256 & r(reps) == 256 & `same' & abs(r(p) - P[1, 1]) < 1e-12 ///
        test_pass "enumeration: 256 sign vectors, RNG state unchanged, p equals boottest's"
    else test_fail "enumeration" "enum `enum' reps `reps'/`r(reps)' same `same' p `=P[1,1]' vs `r(p)'"

    * Intervals (ci): boottest's p-value 1e-3 se inside each bound exceeds
    * the level, and 1e-3 se outside does not
    _pw_data
    cppmlhdfe y x1 x2, absorb(id) vce(cluster cl) wildbootstrap(reps(199) seed(64) ci coefs(x1) nodots)
    matrix CI = e(wb_ci)
    local se = _se[x1]
    poisson y x1 x2 i.id, vce(cluster cl)
    local bad
    foreach side in 1 2 {
        local x = CI[1, `side']
        local dir = cond(`side' == 1, 1, -1)
        _pw_bt, coef(x1) value(`=string(`x' + `dir' * 1e-3 * `se', "%21.17g")') opts(reps(199) seed(64))
        local pin = r(p)
        _pw_bt, coef(x1) value(`=string(`x' - `dir' * 1e-3 * `se', "%21.17g")') opts(reps(199) seed(64))
        if !(`pin' > 0.05 & r(p) <= 0.05) local bad `bad' side`side'=`pin'/`r(p)'
    }
    if "`bad'" == "" test_pass "ci: the bounds are where boottest's p-value crosses the level"
    else test_fail "ci bounds" "`bad'"
}

print_section "Outputs, modes and errors"

_pw_data
cppmlhdfe y x1 x2, absorb(id) vce(cluster cl)
matrix b0 = e(b)
matrix V0 = e(V)
cppmlhdfe y x1 x2, absorb(id) vce(cluster cl) wildbootstrap(reps(49) seed(65) coefs(x2) nodots) nowildtable
local ok = mreldif(b0, e(b)) == 0 & mreldif(V0, e(V)) == 0 & "`e(wb_coefs)'" == "x2" & ///
    colsof(e(wb_p)) == 1 & "`e(wb_method)'" == "score" & e(wb_reps) == 49
capture confirm matrix e(wb_ci)
local noci = _rc != 0
set seed 65
local s65 = c(rngstate)
local ok = `ok' & "`e(wb_rngstate)'" == "`s65'"
capture noisily cppmlhdfe
local rc1 = _rc
capture noisily cppmlhdfe, nowildtable
if `ok' & `noci' & `rc1' == 0 & _rc == 0 ///
    test_pass "e(b), e(V) unchanged; coefs(); no intervals by default; e(wb_rngstate); replay"
else test_fail "outputs" "ok `ok' noci `noci' replay `rc1' `=_rc'"

* With fevce(bootstrap): the score bootstrap first, then the pairs
* replicates, each as alone
cppmlhdfe y x1 x2, absorb(fid=id) vce(cluster cl) fevce(bootstrap, reps(20) seed(66) nodots)
rename fid_se se_alone
drop fid
cppmlhdfe y x1 x2, absorb(id) vce(cluster cl) wildbootstrap(reps(99) seed(67) nodots)
matrix Pw = e(wb_p)
cppmlhdfe y x1 x2, absorb(fid=id) vce(cluster cl) wildbootstrap(reps(99) seed(67) nodots) ///
    fevce(bootstrap, reps(20) seed(66) nodots)
count if reldif(fid_se, se_alone) > 1e-12
if r(N) == 0 & mreldif(Pw, e(wb_p)) == 0 test_pass "wildbootstrap() with fevce(bootstrap): each as alone"
else test_fail "wild with fevce" "rows `r(N)'; p `=mreldif(Pw, e(wb_p))'"

* Refusals
local e1 "absorb(id) vce(bootstrap, reps(10) seed(1)) wildbootstrap(seed(1))"
local e2 "absorb(id) vce(cluster cl) wildbootstrap(seed(1) nonull)"
local e3 "absorb(id) vce(cluster cl t) wildbootstrap(seed(1))"
local e4 "absorb(t) vce(robust) wildbootstrap(seed(1)) [fw=fw]"
local e5 "absorb(id) vce(cluster cl) wildbootstrap(seed(1) weighttype(uniform))"
local e6 "absorb(id) vce(cluster cl) wildbootstrap(seed(1) coefs(off))"
forvalues e = 1/6 {
    local opts `e`e''
    local wexp
    if strpos("`opts'", "[fw=fw]") {
        local wexp [fw=fw]
        local opts = subinstr("`opts'", "[fw=fw]", "", 1)
    }
    _pw_data
    capture cppmlhdfe y x1 x2 `wexp', `opts'
    if _rc == 198 test_pass "refuses: `e`e''"
    else test_fail "refuses: `e`e''" "rc=`=_rc'"
}

noi print_summary "cppmlhdfe_wild"
}

global CTOOLS_COMPONENT_COMPLETE "cppmlhdfe_wild"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_wild"
    exit
}

if "`__case'" == "civreghdfe_savefe" {
* civreghdfe, absorb(..., savefe) and absorb(newvar=absvar)
* (docs/development/BOOTSTRAP_PLAN.md, Phase 4a). The effects are those of
* d = y - Xb - e (the actual regressors, the IV residual), normalized as
* reghdfe's, with the constant in the first effect as ivreghdfe saves them.

do "validation/support/validate_setup.do"

quietly {

noi di as text "Running civreghdfe saved fixed-effect checks..."

capture program drop _cs_data
program define _cs_data
    webuse nlswork, clear
    keep if inrange(idcode, 1, 300) & !missing(ln_wage, tenure, ttl_exp, hours, union)
    gen w = 0.5 + mod(idcode, 7) / 5
end

print_section "Effects against the oracles"

capture which ivreghdfe
local have_iv = !_rc
capture which reghdfe
local have_rh = !_rc

* 2SLS, clusters: ivreghdfe's saved effects, and reghdfe on y - Xb
_cs_data
civreghdfe ln_wage (tenure = union hours) ttl_exp, absorb(fid=idcode year, savefe) vce(cluster idcode)
local saved `e(savefe)'
gen double yd = ln_wage - _b[tenure] * tenure - _b[ttl_exp] * ttl_exp if e(sample)
gen byte smp = e(sample)
rename __hdfe2__ fyr
if `have_rh' {
    reghdfe yd if smp, absorb(o1=idcode o2=year)
    replace o1 = o1 + _b[_cons]
    count if smp & (reldif(fid, o1) > 1e-6 | reldif(fyr, o2) > 1e-6)
    if r(N) == 0 & "`saved'" == "fid __hdfe2__" ///
        test_pass "2SLS: effects equal reghdfe's on y - Xb (constant in the first); e(savefe)"
    else test_fail "2SLS reghdfe oracle" "`r(N)' rows differ; e(savefe) `saved'"
}
else test_skip "reghdfe oracle" "reghdfe is not installed"
if `have_iv' {
    ivreghdfe ln_wage (tenure = union hours) ttl_exp, absorb(i1=idcode i2=year) cluster(idcode)
    count if smp & (reldif(fid, i1) > 1e-6 | reldif(fyr, i2) > 1e-6)
    if r(N) == 0 test_pass "2SLS: effects equal ivreghdfe's"
    else test_fail "2SLS ivreghdfe" "`r(N)' rows differ"

    * LIML with aweights
    _cs_data
    civreghdfe ln_wage (tenure = union hours) ttl_exp [aw=w], absorb(idcode year, savefe) liml
    rename (__hdfe1__ __hdfe2__) (c1 c2)
    ivreghdfe ln_wage (tenure = union hours) ttl_exp [aw=w], absorb(i1=idcode i2=year) liml
    count if e(sample) & (reldif(c1, i1) > 1e-6 | reldif(c2, i2) > 1e-6)
    if r(N) == 0 test_pass "LIML with aweights: effects equal ivreghdfe's"
    else test_fail "LIML aweights" "`r(N)' rows differ"
}
else test_skip "ivreghdfe oracle" "ivreghdfe is not installed"

print_section "Names, results and errors"

* e(b) and e(V) are those of the fit without savefe; one named effect
_cs_data
civreghdfe ln_wage (tenure = union hours) ttl_exp, absorb(idcode year) vce(cluster idcode)
matrix b0 = e(b)
matrix V0 = e(V)
civreghdfe ln_wage (tenure = union hours) ttl_exp, absorb(fid=idcode year) vce(cluster idcode)
capture confirm variable __hdfe2__
local none2 = _rc != 0
local lab : variable label fid
if mreldif(b0, e(b)) == 0 & mreldif(V0, e(V)) == 0 & `none2' & "`e(savefe)'" == "fid" & "`lab'" == "FE for idcode" ///
    test_pass "e(b), e(V) unchanged; only the named effect is saved, labeled"
else test_fail "names" "b `=mreldif(b0, e(b))' none2 `none2' savefe `e(savefe)' label `lab'"

* The second pass after ivreg2's collinearity rules (a duplicated regressor)
_cs_data
gen double ttl2 = 2 * ttl_exp
civreghdfe ln_wage (tenure = union hours) ttl_exp ttl2, absorb(idcode year, savefe)
capture confirm variable __hdfe1__ __hdfe2__
local rc1 = _rc
local sv `e(savefe)'
rename (__hdfe1__ __hdfe2__) (a1 a2)
civreghdfe ln_wage (tenure = union hours) ttl_exp, absorb(idcode year, savefe)
count if e(sample) & (reldif(a1, __hdfe1__) > 1e-9 | reldif(a2, __hdfe2__) > 1e-9)
if `rc1' == 0 & r(N) == 0 & "`sv'" == "__hdfe1__ __hdfe2__" ///
    test_pass "savefe through the collinearity rerun"
else test_fail "rerun" "rc `rc1' rows `r(N)' savefe `sv'"

* Refusals
_cs_data
gen fid = 0
foreach bad in "absorb(idcode year, other)" "absorb(fid=idcode year)" {
    capture civreghdfe ln_wage (tenure = union hours) ttl_exp, `bad'
    local want = cond(strpos("`bad'", "fid="), 110, 198)
    if _rc == `want' test_pass "refuses: `bad'"
    else test_fail "refuses: `bad'" "rc=`=_rc'"
}

noi print_summary "civreghdfe_savefe"
}

global CTOOLS_COMPONENT_COMPLETE "civreghdfe_savefe"
    global CTOOLS_CASE_COMPLETE "civreghdfe_savefe"
    exit
}

if "`__case'" == "civreghdfe_bootstrap" {
* civreghdfe, vce(bootstrap) (docs/development/BOOTSTRAP_PLAN.md, Phase 4c).
* The oracle is Stata's bootstrap on the expanded data: with draws(stata)
* and the same seed, each replicate is ivreghdfe on bsample's draw expanded
* (ivreghdfe cannot run under the bootstrap prefix itself, whose held
* estimates it discards). Clustered GMM2S has each drawn copy of a cluster
* as its own cluster, as idcluster() makes it.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

capture program drop _civbs_data
program define _civbs_data
    webuse nlswork, clear
    xtset, clear
    keep if inrange(idcode, 1, 300) & !missing(ln_wage, tenure, ttl_exp, hours, union, wks_ue, age)
    gen w = 0.5 + mod(idcode, 7) / 5
end

* The replicates on the expanded data from the current estimation sample:
* r(V), r(b_bs) and r(N) (complete replicates) as bootstrap summarizes them
capture program drop _civbs_oracle
program define _civbs_oracle, rclass
    syntax, Reps(integer) SEED(string) SPEC(string asis) [CLuster(varname) STRata(varname) ///
        IDCL OPTS(string asis)]
    set seed `seed'
    local dopts
    if "`cluster'" != "" local dopts `dopts' cluster(`cluster')
    if "`strata'" != "" local dopts `dopts' strata(`strata')
    bsref_draws, reps(`reps') generate(__mm) `dopts'
    gen long __oid = _n
    tempname B V M
    forvalues r = 1/`reps' {
        preserve
        quietly drop if __mm`r' == 0
        quietly expand __mm`r'
        if "`idcl'" != "" {
            bysort __oid: gen long __k = _n
            egen long nid = group(`cluster' __k)
        }
        capture quietly ivreghdfe `spec', `opts'
        local rc = _rc
        if !`rc' matrix `B' = nullmat(`B') \ e(b)
        restore
        if `rc' matrix `B' = nullmat(`B') \ J(1, colsof(`B'), .)
    }
    drop __mm* __oid
    mata: _civbs_B = st_matrix("`B'"); _civbs_B = select(_civbs_B, rowmissing(_civbs_B) :== 0)
    mata: st_matrix("`V'", variance(_civbs_B)); st_matrix("`M'", mean(_civbs_B)); st_local("n", strofreal(rows(_civbs_B)))
    return matrix V = `V'
    return matrix b_bs = `M'
    return scalar N = `n'
end

* civreghdfe's results (estimates ct) against the oracle's r()
capture program drop _civbs_check
program define _civbs_check
    args name tol
    tempname V M
    matrix `V' = r(V)
    matrix `M' = r(b_bs)
    local n = r(N)
    estimates restore ct
    local dv = mreldif(e(V), `V')
    local dm = mreldif(e(b_bs), `M')
    if `dv' < `tol' & `dm' < `tol' & e(N_reps) == `n' & e(N_misreps) == 0 ///
        test_pass "draws(stata) reproduces the bootstrap on the expanded data: `name' (V `dv', mean `dm')"
    else test_fail "expanded-data oracle: `name'" "V `dv' mean `dm' reps `=e(N_reps)'/`n' misreps `=e(N_misreps)'"
end

quietly {

noi di as text "Running civreghdfe vce(bootstrap) checks..."

local s1 "ln_wage (tenure = union hours) ttl_exp"
local s2 "ln_wage (tenure = union hours wks_ue) ttl_exp"

capture which ivreghdfe
if _rc test_skip "civreghdfe bootstrap oracles" "ivreghdfe is not installed"
else {
    print_section "Against the bootstrap on the expanded data"

    _civbs_data
    civreghdfe `s1', absorb(idcode year) vce(bootstrap, reps(20) seed(1) draws(stata) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(1) spec(`s1') opts(absorb(idcode year))
    _civbs_check "2SLS, observations" 1e-9

    _civbs_data
    civreghdfe `s1', absorb(idcode year) liml vce(bootstrap, reps(20) seed(2) draws(stata) cluster(idcode) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(2) spec(`s1') cluster(idcode) opts(absorb(idcode year) liml)
    * (LIML's eigenvalue to about 1e-9 on data held differently)
    _civbs_check "LIML, clusters" 1e-7

    * the moment covariance on the drawn copies, each its own cluster
    _civbs_data
    civreghdfe `s2', absorb(idcode) gmm2s vce(bootstrap, reps(20) seed(3) draws(stata) cluster(idcode) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(3) spec(`s2') cluster(idcode) idcl opts(absorb(idcode) gmm2s cluster(nid))
    _civbs_check "GMM2S, clusters (idcluster)" 1e-9

    * center: each drawn copy of a cluster is centered (BOOTSTRAP_REVIEW.md, H4)
    _civbs_data
    civreghdfe `s2', absorb(idcode) gmm2s center vce(bootstrap, reps(20) seed(3) draws(stata) cluster(idcode) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(3) spec(`s2') cluster(idcode) idcl opts(absorb(idcode) gmm2s cluster(nid) center)
    _civbs_check "GMM2S, center, clusters (idcluster)" 1e-9

    _civbs_data
    civreghdfe `s2', absorb(idcode) gmm2s center vce(bootstrap, reps(20) seed(5) draws(stata) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(5) spec(`s2') opts(absorb(idcode) gmm2s robust center)
    _civbs_check "GMM2S, center, observations" 1e-9

    * aweights, one absorbed effect (the weights the projection and the
    * effects' counts use must agree), robust moment covariance
    _civbs_data
    civreghdfe `s2' [aw=w], absorb(idcode) gmm2s vce(bootstrap, reps(20) seed(4) draws(stata) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(4) spec(`s2' [aw=w]) opts(absorb(idcode) gmm2s robust)
    _civbs_check "GMM2S, aweights, observations" 1e-9

    _civbs_data
    gen byte st = mod(idcode, 2)
    civreghdfe ln_wage (tenure = union hours) ttl_exp age, partial(age) absorb(idcode) fuller(1) ///
        vce(bootstrap, reps(20) seed(5) draws(stata) cluster(idcode) strata(st) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(20) seed(5) spec(ln_wage (tenure = union hours) ttl_exp age) cluster(idcode) ///
        strata(st) opts(partial(age) absorb(idcode) fuller(1))
    _civbs_check "Fuller, partial(), strata and clusters" 1e-9

    _civbs_data
    civreghdfe ln_wage (ttl_exp = age hours) tenure [pw=w], absorb(idcode) cue ///
        vce(bootstrap, reps(10) seed(6) draws(stata) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(10) seed(6) spec(ln_wage (ttl_exp = age hours) tenure [pw=w]) opts(absorb(idcode) cue robust)
    * (two optimizers)
    _civbs_check "CUE, pweights" 1e-5

    * A weakly identified CUE model: some replicates start where the
    * objective is not convex (two of these ten ended at worse optima before
    * the absolute-eigenvalue Newton step)
    _civbs_data
    civreghdfe `s2', absorb(idcode) cue vce(bootstrap, reps(10) seed(6) draws(stata) nodots)
    estimates store ct
    keep if e(sample)
    _civbs_oracle, reps(10) seed(6) spec(`s2') opts(absorb(idcode) cue robust)
    _civbs_check "CUE, weakly identified" 1e-5
}

print_section "Results and reproducibility"

* fweights with clusters: the clusters drawn are those of the expanded data
_civbs_data
gen fw = 1 + mod(_n, 3)
civreghdfe `s1' [fw=fw], absorb(idcode year) vce(bootstrap, reps(20) seed(9) cluster(idcode) nodots)
matrix V1 = e(V)
expand fw
civreghdfe `s1', absorb(idcode year) vce(bootstrap, reps(20) seed(9) cluster(idcode) nodots)
if mreldif(V1, e(V)) < 1e-10 test_pass "fweights with clusters: the expanded data's bootstrap"
else test_fail "fweights, clusters" "reldif `=mreldif(V1, e(V))'"

* draws(fast): reproducible from the seed
_civbs_data
civreghdfe `s1', absorb(idcode) vce(bootstrap, reps(20) seed(10) draws(fast) nodots)
matrix A = e(V)
civreghdfe `s1', absorb(idcode) vce(bootstrap, reps(20) seed(10) draws(fast) nodots)
matrix B = e(V)
civreghdfe `s1', absorb(idcode) vce(bootstrap, reps(20) seed(11) draws(fast) nodots)
if mreldif(A, B) == 0 & mreldif(A, e(V)) > 0 test_pass "draws(fast): reproducible from the seed"
else test_fail "draws(fast) seed" "same-seed reldif `=mreldif(A, B)'"

* what vce(bootstrap) posts; the base VCE's covariance in e(V_modelbased)
* and its diagnostics (cluster on the bootstrap's cluster(), else robust)
_civbs_data
foreach base in "robust" "cluster idcode" {
    local bopt = cond("`base'" == "robust", "", "cluster(idcode)")
    civreghdfe `s1', absorb(idcode year) vce(`base')
    matrix Vb = e(V)
    local idst = e(idstat)
    local jst = e(j)
    local kpf = e(widstat)
    civreghdfe `s1', absorb(idcode year) vce(bootstrap, reps(20) seed(12) `bopt' nodots)
    local ok = ("`e(vce)'" == "bootstrap" & "`e(vcetype)'" == "Bootstrap" & "`e(prefix)'" == "bootstrap" ///
        & e(N_reps) == 20 & e(N_misreps) == 0 & e(df_r) == . & e(chi2) < . & e(bs_version) == 3 ///
        & "`e(rngstate)'" != "" & "`e(draws)'" == "stata" & "`e(chi2type)'" == "Wald")
    foreach m in b_bs reps bias se z0 ci_normal ci_percentile ci_bc boot_status {
        capture confirm matrix e(`m')
        if _rc local ok = 0
    }
    if `ok' & mreldif(e(V_modelbased), Vb) < 1e-12 & reldif(e(idstat), `idst') < 1e-12 & ///
        reldif(e(j), `jst') < 1e-12 & reldif(e(widstat), `kpf') < 1e-12 ///
        test_pass "vce(bootstrap) posts e() as Stata's bootstrap does; base vce(`base') in e(V_modelbased) and the diagnostics"
    else test_fail "e() with base vce(`base')" "vce=`e(vce)' vcetype=`e(vcetype)' N_reps=`e(N_reps)' df_r=`e(df_r)'"
}
capture noisily estat bootstrap, all
if !_rc test_pass "estat bootstrap after civreghdfe, vce(bootstrap)"
else test_fail "estat bootstrap" "rc=`=_rc'"
capture noisily civreghdfe
if !_rc & "`e(vce)'" == "bootstrap" test_pass "replay after vce(bootstrap)"
else test_fail "replay" "rc=`=_rc'"

* ivreg2's collinearity rerun bootstraps the final model once
_civbs_data
gen double ttl2 = 2 * ttl_exp
civreghdfe `s1' ttl2, absorb(idcode) vce(bootstrap, reps(20) seed(13) nodots)
local a1 = _se[tenure]
local a2 = _se[ttl_exp]
local cl `e(collin)'
civreghdfe `s1', absorb(idcode) vce(bootstrap, reps(20) seed(13) nodots)
local d = max(reldif(`a1', _se[tenure]), reldif(`a2', _se[ttl_exp]))
if `d' < 1e-12 & "`cl'" == "ttl2" test_pass "the collinearity rerun bootstraps the final model"
else test_fail "rerun" "reldif `d' collin `cl'"

print_section "Errors"

_civbs_data
tsset idcode year
matrix b0 = (0.1, 0.1)
foreach bad in "vce(bootstrap) bw(2)" "vce(bootstrap) dkraay(2)" "vce(bootstrap) kiefer" ///
    "vce(bootstrap) b0(b0)" "vce(bootstrap, reps(1))" "vce(bootstrap, cluster(idcode year))" ///
    "robust vce(bootstrap)" {
    capture civreghdfe `s1', absorb(idcode) `bad'
    if _rc == 198 test_pass "civreghdfe refuses `bad'"
    else test_fail "civreghdfe `bad'" "rc=`=_rc'"
}
* GMM2S needs a moment covariance each replicate can re-estimate
capture civreghdfe `s1', absorb(idcode, savefe) gmm2s vce(cluster idcode year) fevce(bootstrap)
if _rc == 198 test_pass "fevce(bootstrap) with GMM2S and two-way clustering is refused"
else test_fail "gmm2s two-way fevce" "rc=`=_rc'"
capture drop __hdfe1__

noi print_summary "civreghdfe_bootstrap"
}

global CTOOLS_COMPONENT_COMPLETE "civreghdfe_bootstrap"
    global CTOOLS_CASE_COMPLETE "civreghdfe_bootstrap"
    exit
}

if "`__case'" == "civreghdfe_fe_bootstrap" {
* civreghdfe, fevce(bootstrap) (docs/development/BOOTSTRAP_PLAN.md, Phase 4c):
* the bootstrap standard errors of the saved effects and of their sum,
* against the effects civreghdfe saves on each replicate's expanded data.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

capture program drop _civfe_data
program define _civfe_data
    webuse nlswork, clear
    xtset, clear
    keep if inrange(idcode, 1, 300) & !missing(ln_wage, tenure, ttl_exp, hours, union)
end

quietly {

noi di as text "Running civreghdfe fevce(bootstrap) checks..."

local s1 "ln_wage (tenure = union hours) ttl_exp"

print_section "Against the expanded data"

foreach ab in "idcode" "idcode year" {
    _civfe_data
    local nab : word count `ab'
    civreghdfe `s1', absorb(`ab', savefe) fevce(bootstrap, reps(12) seed(8) draws(stata) cluster(idcode) ///
        dse(dse) nreps nodots)
    local fese `e(fe_se)'
    local split = e(fe_split_reps)
    keep if e(sample)
    gen long __oid = _n
    set seed 8
    bsref_draws, reps(12) generate(__mm) cluster(idcode)
    forvalues r = 1/12 {
        preserve
        drop if __mm`r' == 0
        expand __mm`r'
        drop __hdfe*
        civreghdfe `s1', absorb(`ab', savefe)
        gen double __d = __hdfe1__
        if `nab' == 2 replace __d = __d + __hdfe2__
        collapse (first) r1_`r' = __hdfe1__ rd_`r' = __d, by(__oid)
        tempfile f`r'
        save `f`r''
        restore
        merge 1:1 __oid using `f`r'', nogenerate
    }
    egen double o1 = rowsd(r1_*)
    egen double od = rowsd(rd_*)
    egen n1 = rownonmiss(r1_*)
    count if n1 >= 2 & (reldif(o1, __fese1__) > 1e-6 | reldif(od, dse) > 1e-6)
    local bad = r(N)
    count if n1 != __fenrep1__
    if `bad' == 0 & r(N) == 0 & "`fese'" == "__fese1__`=cond(`nab' == 2, " __fese2__", "")'" & `split' == 0 ///
        test_pass "absorb(`ab'): the effects' and d's SEs and replicate counts equal the expanded data's"
    else test_fail "absorb(`ab') FE SEs" "`bad' rows differ, `r(N)' counts differ; e(fe_se) `fese'; split `split'"
}

print_section "Shared replicates, names and errors"

* vce(bootstrap) with fevce(bootstrap): one set of replicates for both
_civfe_data
civreghdfe `s1', absorb(idcode year, savefe) vce(bootstrap, reps(15) seed(20) cluster(idcode) nodots) fevce(bootstrap)
matrix V1 = e(V)
rename (__fese1__ __fese2__) (a1 a2)
drop __hdfe*
civreghdfe `s1', absorb(idcode year) vce(bootstrap, reps(15) seed(20) cluster(idcode) nodots)
matrix V2 = e(V)
civreghdfe `s1', absorb(idcode year, savefe) vce(cluster idcode) fevce(bootstrap, reps(15) seed(20) cluster(idcode) nodots)
count if e(sample) & (reldif(a1, __fese1__) > 1e-12 | reldif(a2, __fese2__) > 1e-12)
if mreldif(V1, V2) == 0 & r(N) == 0 test_pass "vce() and fevce() share the replicates"
else test_fail "shared replicates" "V `=mreldif(V1, V2)', `r(N)' rows differ"

* names: v_se for newvar=, sestub(); e() and the normalization warning
_civfe_data
civreghdfe `s1', absorb(fid=idcode year) fevce(bootstrap, reps(10) seed(21) nodots)
capture confirm variable fid_se
local r1 = _rc
local lab : variable label fid_se
civreghdfe `s1', absorb(idcode, savefe) fevce(bootstrap, reps(10) seed(21) sestub(sefe) nodots)
capture confirm variable sefe1
if `r1' == 0 & _rc == 0 & "`lab'" == "Bootstrap SE of FE for idcode" & "`e(fevce)'" == "bootstrap" & ///
    e(fe_N_reps) == 10 & e(fe_norm_warning) == 0 & rowsof(e(fe_bs_table)) == 1 ///
    test_pass "SE variable names (v_se, sestub), labels and e() results"
else test_fail "names" "fid_se rc `r1', sefe1 rc `=_rc', label `lab'"

_civfe_data
foreach bad in "absorb(idcode) fevce(bootstrap)" ///
    "absorb(idcode, savefe) fevce(robust)" ///
    "absorb(idcode, savefe) vce(bootstrap, reps(5)) fevce(bootstrap, reps(6))" {
    capture civreghdfe `s1', `bad'
    if _rc == 198 test_pass "civreghdfe refuses `bad'"
    else test_fail "civreghdfe `bad'" "rc=`=_rc'"
}

noi print_summary "civreghdfe_fe_bootstrap"
}

global CTOOLS_COMPONENT_COMPLETE "civreghdfe_fe_bootstrap"
    global CTOOLS_CASE_COMPLETE "civreghdfe_fe_bootstrap"
    exit
}

if "`__case'" == "civreghdfe_wild" {
* civreghdfe, wildbootstrap(): the WRE bootstrap (docs/development/
* BOOTSTRAP_PLAN.md, Phase 4d). Reference: boottest 4.5.5 after ivreg2 with
* the absorbed effects as dummies, with the same seed: p-values must be
* equal, and interval bounds are checked through boottest's p-value just
* inside and just outside each bound (boottest stops its own search within
* one replicate step of the level). boottest after ivreghdfe with absorb()
* is not a reference: its WRE statistic for an exogenous coefficient is not
* the regression's. With LIML, boottest's p-value function steps up to
* about 1e-4 standard errors away from the exact one (its LIML replicates
* are less precise), so its bounds are checked at 1e-3 standard errors and
* the Mata WRE checks them at 1e-6. Fuller, k-class and coarser bootstrap
* clusters, which boottest does not take, use the replicate-by-replicate
* WRE in Mata (_bsref_wre_t) on the same draws.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

quietly {

noi di as text "Running civreghdfe wildbootstrap() checks..."

capture program drop _civwb_bt
program define _civwb_bt, rclass
    syntax, coef(string) [value(real 0) opts(string)]
    capture noisily boottest (`coef' = `value'), nograph noci `opts'
    if _rc {
        return scalar p = .
        return scalar reps = .
        exit
    }
    return scalar p = r(p)
    return scalar reps = r(reps)
end

* 30 firms of 20 rows, years crossing them, weak effects
capture program drop _civwb_data
program define _civwb_data
    clear
    set seed 2027
    set obs 600
    gen firm = ceil(_n / 20)
    gen year = mod(_n * 7, 11) + 1
    gen grp = ceil(firm / 4)
    gen uu = rnormal()
    gen ff = rnormal()
    gen yy = rnormal() * .5
    bysort year: replace yy = yy[1]
    gen z1 = rnormal() + .3 * ff
    gen z2 = rnormal()
    gen z3 = rnormal()
    gen wx = rnormal() + .2 * ff
    gen x = .6 * z1 + .4 * z2 + .2 * wx + .5 * uu + rnormal() + ff + yy
    gen x2 = .5 * z2 - .4 * z3 + .3 * z1 + .4 * uu + rnormal() - .5 * ff
    gen y = 1 + .03 * x + .02 * wx + 2 * uu + ff + yy
    gen y2 = 1 + .03 * x - .02 * x2 + .02 * wx + 2 * uu + ff + yy
    gen aw = .5 + mod(_n * 13, 7) / 4
    gen fw = 1 + mod(_n, 3)
    sort firm year, stable
    drop uu ff yy
end

* One design: civreghdfe against boottest after the reference fit. Returns
* r(bad) (empty when every p-value and bound agrees) and r(path).
capture program drop _civwb_design
program define _civwb_design, rclass
    syntax, civ(string asis) ref(string asis) btopts(string) [noci tol(real 1e-5)]
    civreghdfe `civ'
    matrix P = e(wb_p)
    local names : colnames P
    local path `e(wb_fe_path)'
    local hasci = "`ci'" == ""
    if `hasci' matrix CI = e(wb_ci)
    `ref'
    local bad
    local j 0
    foreach c of local names {
        local ++j
        _civwb_bt, coef(`c') opts(`btopts')
        if !(abs(r(p) - P[1, `j']) < 1e-12) local bad `bad' p(`c')=`=P[1,`j']'/`r(p)'
        if !`hasci' continue
        local s = _se[`c']
        foreach side in 1 2 {
            local x = CI[`j', `side']
            local dir = cond(`side' == 1, 1, -1)
            _civwb_bt, coef(`c') value(`=string(`x' + `dir' * `tol' * `s', "%21.17g")') opts(`btopts')
            local pin = r(p)
            _civwb_bt, coef(`c') value(`=string(`x' - `dir' * `tol' * `s', "%21.17g")') opts(`btopts')
            local pout = r(p)
            if !(`pin' > 0.05 & `pout' <= 0.05) local bad `bad' ci(`c',`side')=`pin'/`pout'
        }
    }
    return local bad `bad'
    return local path `path'
end

capture which boottest
local have_bt = !_rc
capture which ivreg2
local have_bt = `have_bt' & !_rc
if !`have_bt' {
    test_skip "wildbootstrap() against boottest" "boottest or ivreg2 is not installed"
}
else {
    print_section "p-values and intervals against boottest (effects as dummies)"

    local m1 "y (x = z1 z2) wx"
    local m2 "y2 (x x2 = z1 z2 z3) wx"
    * label | civreghdfe | ivreg2 | boottest options | noci
    local d1 "nested, Rademacher|`m1', absorb(firm) vce(cluster firm)|`m1' i.firm, cluster(firm)|reps(199) seed(1)|"
    local d2 "crossed (projected), Webb|`m1', absorb(year) vce(cluster firm)|`m1' i.year, cluster(firm)|reps(199) seed(2) weighttype(webb)|"
    local d3 "LIML, Mammen|`m1', absorb(firm) vce(cluster firm) liml|`m1' i.firm, cluster(firm) liml|reps(199) seed(3) weighttype(mammen)|"
    local d4 "aweights, normal|`m1' [aw=aw], absorb(firm) vce(cluster firm)|`m1' i.firm [aw=aw], cluster(firm)|reps(199) seed(4) weighttype(normal)|"
    local d5 "robust, each observation (projected), gamma|`m1', absorb(firm) vce(robust)|`m1' i.firm, robust|reps(99) seed(5) weighttype(gamma)|"
    local d6 "two endogenous|`m2', absorb(firm) vce(cluster firm)|`m2' i.firm, cluster(firm)|reps(199) seed(8)|"
    // Check the LIML p-values against boottest here. Its test inversion can
    // cross at different points with two endogenous regressors. The Mata WRE
    // below checks both sides of all four interval bounds directly.
    local d7 "two endogenous, LIML (projected)|`m2', absorb(year) vce(cluster firm) liml|`m2' i.year, cluster(firm) liml|reps(199) seed(11)|noci"
    local d8 "equal-tailed|`m1', absorb(year) vce(cluster firm)|`m1' i.year, cluster(firm)|reps(199) seed(7) ptype(equaltail)|"
    local d9 "two absorbed effects|`m1', absorb(firm year) vce(cluster firm)|`m1' i.firm i.year, cluster(firm)|reps(199) seed(9)|noci"
    local d10 "pweights|`m1' [pw=aw], absorb(year) vce(cluster firm)|`m1' i.year [pw=aw], cluster(firm)|reps(199) seed(12)|noci"
    local d11 "unrestricted (WUE)|`m1', absorb(firm) vce(cluster firm)|`m1' i.firm, cluster(firm)|reps(199) seed(6) nonull|noci"
    local d12 "partial(wx), as with wx a regressor|`m1', absorb(year) vce(cluster firm) partial(wx)|`m1' i.year, cluster(firm)|reps(199) seed(2) weighttype(webb)|noci"
    forvalues d = 1/12 {
        gettoken lab rest : d`d', parse("|")
        gettoken bar rest : rest, parse("|")
        gettoken civ rest : rest, parse("|")
        gettoken bar rest : rest, parse("|")
        gettoken ref rest : rest, parse("|")
        gettoken bar rest : rest, parse("|")
        gettoken bto rest : rest, parse("|")
        gettoken bar noci : rest, parse("|")
        _civwb_data
        local tol = cond(strpos("`civ'", "liml"), 1e-3, 1e-5)
        _civwb_design, civ(`civ' wildbootstrap(`bto' nodots)) ref(quietly ivreg2 `ref') btopts(`bto') `noci' tol(`tol')
        if "`r(bad)'" == "" test_pass "`lab' (`r(path)'): p-values`=cond("`noci'" == "", " and interval bounds", "")' equal boottest's"
        else test_fail "`lab'" "`r(bad)'"
    }

    * Enumeration with 9 clusters: 512 sign vectors, no draw
    _civwb_data
    keep if firm <= 9
    set seed 34
    local state0 = c(rngstate)
    civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) wildbootstrap(nodots)
    local same = (c(rngstate) == "`state0'")
    matrix P = e(wb_p)
    local reps = e(wb_reps)
    local enum = e(wb_enumerated)
    local C = e(wb_N_clust)
    quietly ivreg2 y (x = z1 z2) wx i.firm, cluster(firm)
    _civwb_bt, coef(x) opts(reps(999))
    local p1 = r(p)
    local r1 = r(reps)
    _civwb_bt, coef(wx) opts(reps(999))
    if `enum' == 1 & `reps' == 2^`C' & `r1' == `reps' & `same' & abs(`p1' - P[1, 1]) < 1e-12 & ///
        abs(r(p) - P[1, 2]) < 1e-12 ///
        test_pass "enumeration: all 2^`C' sign vectors, RNG state unchanged, p equals boottest's"
    else test_fail "enumeration" "enum `enum' reps `reps'/`r1' same `same' p `=P[1,1]' `=P[1,2]' vs `p1' `r(p)'"
}

print_section "Fuller, k-class and coarser bootstrap clusters against the Mata WRE"

* The firm effects as dummies, the draws regenerated from the seed
* label | options | est | kc | bootcluster
local o1 "Fuller(1)|fuller(1)|2|`=1 / (600 - 3)'|"
local o2 "k-class(0.8)|kclass(.8)|3|.8|"
local o3 "LIML, bootcluster(grp): 8 groups|liml|1|0|grp"
forvalues o = 1/3 {
    gettoken lab rest : o`o', parse("|")
    gettoken bar rest : rest, parse("|")
    gettoken eopt rest : rest, parse("|")
    gettoken bar rest : rest, parse("|")
    gettoken est rest : rest, parse("|")
    gettoken bar rest : rest, parse("|")
    gettoken kc rest : rest, parse("|")
    gettoken bar bcv : rest, parse("|")
    local bcopt = cond("`bcv'" != "", "bootcluster(`bcv')", "")
    _civwb_data
    civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) `eopt' ///
        wildbootstrap(reps(199) seed(41) weighttype(webb) `bcopt' nodots)
    matrix P = e(wb_p)
    local C = e(wb_N_clust)
    capture drop D_*
    quietly tabulate firm, generate(D_)
    if "`bcv'" == "" local bcv firm
    set seed 41
    _ctools_hdfe_boot wilddraw, clusters(`C') reps(199) type(webb)
    local bad
    forvalues k = 1/2 {
        mata: _wb_t = _bsref_wre_t(st_data(., "y"), st_data(., "x"), (st_data(., "wx"), st_data(., "D_*")), ///
            st_data(., "z1 z2"), J(st_nobs(), 1, 1), st_data(., "firm"), st_data(., "`bcv'"), `k', 0, `est', `kc', _ctools_wild_v)
        mata: st_local("p", strofreal(sum(abs(_wb_t[|2 \ .|]) :> abs(_wb_t[1])) / 199, "%21.17g"))
        if !(abs(`p' - P[1, `k']) < 1e-12) local bad `bad' p`k'=`=P[1,`k']'/`p'
    }
    mata: mata drop _wb_t
    _ctools_hdfe_boot wildclear
    if "`bad'" == "" test_pass "`lab': p-values equal the Mata WRE's"
    else test_fail "`lab'" "`bad'; clusters `C'"
}

* Interval bounds of the two-endogenous LIML design (projected), exactly:
* the Mata WRE's p-value just inside and just outside each bound
_civwb_data
civreghdfe y2 (x x2 = z1 z2 z3) wx, absorb(year) vce(cluster firm) liml ///
    wildbootstrap(reps(199) seed(11) nodots)
matrix CI = e(wb_ci)
local C = e(wb_N_clust)
capture drop D_*
quietly tabulate year, generate(D_)
set seed 11
_ctools_hdfe_boot wilddraw, clusters(`C') reps(199) type(rademacher)
local bad
forvalues k = 1/2 {
    local s = _se[`: word `k' of x x2']
    foreach side in 1 2 {
        local dir = cond(`side' == 1, 1, -1)
        foreach w in in out {
            local b0 = CI[`k', `side'] + cond("`w'" == "in", 1, -1) * `dir' * 1e-6 * `s'
            mata: _wb_t = _bsref_wre_t(st_data(., "y2"), st_data(., "x x2"), (st_data(., "wx"), st_data(., "D_*")), ///
                st_data(., "z1 z2 z3"), J(st_nobs(), 1, 1), st_data(., "firm"), st_data(., "firm"), `k', `b0', 1, 0, _ctools_wild_v)
            mata: st_local("p`w'", strofreal(sum(abs(_wb_t[|2 \ .|]) :> abs(_wb_t[1])) / 199, "%21.17g"))
        }
        if !(`pin' > 0.05 & `pout' <= 0.05) local bad `bad' ci(`k',`side')=`pin'/`pout'
    }
}
mata: mata drop _wb_t
_ctools_hdfe_boot wildclear
if "`bad'" == "" test_pass "LIML, two endogenous (projected): interval bounds exact against the Mata WRE (1e-6 se)"
else test_fail "LIML interval bounds" "`bad'"

print_section "Paths, outputs and errors"

* The projected path on a nested design gives the aggregate path's results
_civwb_data
civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) wildbootstrap(reps(99) seed(35) nodots)
matrix Pa = e(wb_p)
matrix Ca = e(wb_ci)
local pa `e(wb_fe_path)'
civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) wildbootstrap(reps(99) seed(35) nodots _projected)
matrix Pp = e(wb_p)
matrix Cp = e(wb_ci)
local pp `e(wb_fe_path)'
if mreldif(Pa, Pp) == 0 & mreldif(Ca, Cp) < 1e-8 & "`pa'" == "nested" & "`pp'" == "projected" ///
    test_pass "projected and aggregate paths agree on a nested design"
else test_fail "path agreement" "p `=mreldif(Pa, Pp)' ci `=mreldif(Ca, Cp)' paths `pa' `pp'"

* fedgp(partialled) on a crossed design takes the aggregate path
civreghdfe y (x = z1 z2) wx, absorb(year) vce(cluster firm) wildbootstrap(reps(49) seed(36) fedgp(partialled) noci nodots)
if "`e(wb_fe_path)'" == "partialled" & colsof(e(wb_p)) == 2 test_pass "fedgp(partialled): the partialled path"
else test_fail "fedgp(partialled)" "path `e(wb_fe_path)'"

* e(b), e(V) and the diagnostics are unchanged; coefs(), noci, replay
civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm)
matrix b0 = e(b)
matrix V0 = e(V)
local j0 = e(j)
civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) ///
    wildbootstrap(reps(49) seed(36) coefs(wx) noci nodots) nowildtable
local ok = mreldif(b0, e(b)) == 0 & mreldif(V0, e(V)) == 0 & reldif(`j0', e(j)) == 0 & ///
    "`e(wb_coefs)'" == "wx" & colsof(e(wb_p)) == 1 & "`e(wb_method)'" == "wre" & e(wb_reps) == 49 & ///
    reldif(el(e(wb_t), 1, 1), _b[wx] / _se[wx]) < 1e-12
capture confirm matrix e(wb_ci)
local noci = _rc != 0
capture noisily civreghdfe
if `ok' & `noci' & _rc == 0 test_pass "e(b), e(V), e(j) unchanged; coefs(), noci, nowildtable, replay"
else test_fail "outputs" "ok `ok' noci `noci' replay rc `=_rc'"

* nonull posts the WUE label
civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) wildbootstrap(reps(49) seed(36) nonull noci nodots)
if "`e(wb_method)'" == "wue" test_pass "nonull: e(wb_method) wue"
else test_fail "nonull label" "`e(wb_method)'"

* With fevce(bootstrap): both run on one fit, each as alone
civreghdfe y (x = z1 z2) wx, absorb(firm, savefe) vce(cluster firm) fevce(bootstrap, reps(20) seed(37) nodots)
rename __fese1__ se_alone
drop __hdfe1__
civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) wildbootstrap(reps(99) seed(38) nodots)
matrix Pw = e(wb_p)
civreghdfe y (x = z1 z2) wx, absorb(firm, savefe) vce(cluster firm) ///
    wildbootstrap(reps(99) seed(38) nodots) fevce(bootstrap, reps(20) seed(37) nodots)
count if reldif(__fese1__, se_alone) > 1e-12
if r(N) == 0 & mreldif(Pw, e(wb_p)) == 0 & "`e(fevce)'" == "bootstrap" ///
    test_pass "wildbootstrap() with fevce(bootstrap): each as when run alone"
else test_fail "wild with fevce" "SE rows differ `r(N)'; p `=mreldif(Pw, e(wb_p))'"
drop __hdfe1__ __fese1__ se_alone

* wildbootstrap without options, in all three commands
_civwb_data
capture civreghdfe y (x = z1 z2) wx, absorb(firm) vce(cluster firm) wildbootstrap
local r1 = _rc
local m1 `e(wb_method)'
capture creghdfe y x wx, absorb(firm) vce(cluster firm) wildbootstrap
local r2 = _rc
local m2 `e(wb_method)'
gen ey = exp(y / 4)
capture cppmlhdfe ey x wx, absorb(firm) vce(cluster firm) wildbootstrap
local r3 = _rc
local m3 `e(wb_method)'
if `r1' == 0 & `r2' == 0 & `r3' == 0 & "`m1' `m2' `m3'" == "wre wcr score" ///
    test_pass "wildbootstrap without options: civreghdfe, creghdfe, cppmlhdfe"
else test_fail "bare wildbootstrap" "rc `r1' `r2' `r3'; methods `m1' `m2' `m3'"

* Refusals (plan 4.6)
local e1 "absorb(firm) wildbootstrap(seed(1))"
local e2 "absorb(firm) vce(bootstrap, reps(10) seed(1)) wildbootstrap(seed(1))"
local e3 "absorb(year) vce(robust) wildbootstrap(seed(1)) [fw=fw]"
local e4 "absorb(firm) vce(cluster firm) wildbootstrap(seed(1) coefs(z1))"
local e5 "absorb(firm) vce(cluster firm) gmm2s wildbootstrap(seed(1))"
local e6 "absorb(firm) vce(cluster firm) cue wildbootstrap(seed(1))"
local e7 "absorb(firm) vce(cluster firm year) wildbootstrap(seed(1))"
local e8 "absorb(firm) vce(cluster firm) wildbootstrap(seed(1) bootcluster(grp year))"
forvalues e = 1/8 {
    local opts `e`e''
    local wexp
    if strpos("`opts'", "[fw=fw]") {
        local wexp [fw=fw]
        local opts = subinstr("`opts'", "[fw=fw]", "", 1)
    }
    capture civreghdfe y (x = z1 z2) wx `wexp', `opts'
    if _rc == 198 test_pass "refuses: `e`e''"
    else test_fail "refuses: `e`e''" "rc=`=_rc'"
}

noi print_summary "civreghdfe_wild"
}

global CTOOLS_COMPONENT_COMPLETE "civreghdfe_wild"
    global CTOOLS_CASE_COMPLETE "civreghdfe_wild"
    exit
}

if "`__case'" == "bootstrap_options" {
* The pairs-bootstrap options added in Phase 5 (docs/development/BOOTSTRAP_PLAN.md):
* saving(), the hidden _bsdraws(), draws(stata) with fweighted rows and with
* observations resampled within strata, fevce(bootstrap, ci), heterogeneous
* slopes in absorb() for every bootstrap, and civreghdfe's savefe and fevce()
* with partial(). Oracles: the bootstrap prefix (same seed), reghdfe and
* ppmlhdfe fitted on each replicate's draws, and the same model without the
* option.

do "validation/support/validate_setup.do"
do "validation/support/bootstrap_reference.do"

capture program drop _bo_pf_data
program define _bo_pf_data
    clear
    set seed 78
    set obs 480
    gen id = ceil(_n / 6)
    gen t = mod(_n - 1, 6) + 1
    gen cl = ceil(id / 4)
    gen x1 = rnormal()
    gen double y = rpoisson(exp(0.3 * x1 + mod(id, 7) / 5 + t / 10))
end

capture program drop _bo_pctile
program define _bo_pctile
    * Percentile bounds over the replicates f_1..f_R of each row (as
    * e(ci_percentile): _pctile at (100 - level) / 2 and (100 + level) / 2,
    * over the replicates that are not missing) into o_lb and o_ub
    syntax, stub(name) reps(integer) [level(real 95)]
    quietly generate double o_lb = .
    quietly generate double o_ub = .
    mata: _bo_pctile_rows("`stub'", `reps', `level')
end
capture mata: mata drop _bo_pctile_rows()
mata:
void _bo_pctile_rows(string scalar stub, real scalar R, real scalar level)
{
    real matrix F
    real colvector x
    real scalar i, n, p, k
    real rowvector q
    string rowvector names

    names = J(1, R, "")
    for (k = 1; k <= R; k++) names[k] = stub + strofreal(k)
    F = st_data(., names)
    q = ((100 - level) / 2, (100 + level) / 2)
    for (i = 1; i <= rows(F); i++) {
        x = select(F[i, .], F[i, .] :< .)'
        n = rows(x)
        if (n < 2) continue
        x = sort(x, 1)
        // _pctile's default definition
        for (k = 1; k <= 2; k++) {
            p = n * q[k] / 100
            if (p == floor(p)) st_store(i, (k == 1 ? "o_lb" : "o_ub"), (x[p] + x[p + 1]) / 2)
            else st_store(i, (k == 1 ? "o_lb" : "o_ub"), x[floor(p) + 1])
        }
    }
}
end

quietly {

noi di as text "Running bootstrap option checks..."

local have_reghdfe = 1
capture which reghdfe
if _rc local have_reghdfe = 0
local have_ppmlhdfe = 1
capture which ppmlhdfe
if _rc local have_ppmlhdfe = 0

tempfile f1 f2

print_section "saving()"

if `have_reghdfe' {
    * The replicates, names, labels and characteristics bootstrap writes
    * (an omitted regressor included), and bstat reads the file
    sysuse auto, clear
    drop if missing(rep78)
    gen m2 = 2 * mpg
    creghdfe price mpg m2 weight, absorb(rep78) ///
        vce(bootstrap, reps(6) seed(3) draws(stata) cluster(rep78) strata(foreign) saving(`f1', replace) nodots)
    matrix V1 = e(V)
    bootstrap _b, reps(6) seed(3) cluster(rep78) strata(foreign) nodots saving(`f2', replace): ///
        reghdfe price mpg m2 weight, absorb(rep78)
    matrix V2 = e(V)
    preserve
    use `f1', clear
    local n1 = _N
    local v1 : char _dta[N_cluster]
    local s1 : char _dta[N_strata]
    local l1 : variable label _bs_2
    local c1 : char _b_weight[colname]
    unab names1 : _all
    mkmat _all, matrix(R1)
    use `f2', clear
    unab names2 : _all
    mkmat _all, matrix(R2)
    restore
    local ok = ("`names1'" == "`names2'" & `n1' == 6 & "`v1'" == "8" & "`s1'" == "2" & ///
        "`l1'" == "_b[o.m2]" & "`c1'" == "weight")
    if `ok' & mreldif(R1, R2) < 1e-6 & mreldif(V1, V2) < 1e-6 ///
        test_pass "creghdfe saving(): the prefix's replicates, names, labels and characteristics"
    else test_fail "creghdfe saving()" "names `names1' | `names2', N `n1', clusters `v1', strata `s1', label `l1', reps `=mreldif(R1, R2)'"
    capture noisily bstat using `f1'
    if !_rc & mreldif(e(V), V1) < 1e-6 test_pass "bstat reads the saving() file"
    else test_fail "bstat" "rc=`=_rc'"
}

* double: the file's covariance is e(V); every(); more replicates than a
* chunk; a regressor dropped in some replicates is missing there
sysuse auto, clear
gen rare = _n <= 3
creghdfe price mpg rare, absorb(foreign) vce(bootstrap, reps(120) seed(4) saving(`f1', replace double every(7)) nodots)
matrix V1 = e(V)
local nm = e(N_misreps)
preserve
use `f1', clear
local n1 = _N
count if missing(_b_rare) & missing(_b_mpg) & missing(_b_cons)
local nmiss = r(N)
local type : type _b_mpg
drop if missing(_b_rare)
mata: st_matrix("Vf", variance(st_data(., .)))
restore
if `n1' == 120 & `nmiss' == `nm' & `nm' > 0 & "`type'" == "double" & mreldif(V1, Vf) < 1e-12 ///
    test_pass "saving(, double every()): 120 replicates over chunks, failed ones missing, variance = e(V)"
else test_fail "saving(double)" "N `n1' missing `nmiss' vs `nm' type `type' V `=mreldif(V1, Vf)'"
* without replace, an existing file is an error before any work
capture creghdfe price mpg, absorb(foreign) vce(bootstrap, reps(5) seed(4) saving(`f1'))
if _rc == 602 test_pass "saving() without replace refuses an existing file"
else test_fail "saving() replace" "rc=`=_rc'"
* noconstant: no _cons in e(b) or in the file
creghdfe price mpg weight, absorb(foreign) noconstant vce(bootstrap, reps(5) seed(4) saving(`f1', replace) nodots)
local cn : colnames e(b)
preserve
use `f1', clear
unab vars : _all
restore
if "`cn'" == "mpg weight" & "`vars'" == "_b_mpg _b_weight" test_pass "saving() with noconstant"
else test_fail "saving() noconstant" "e(b) `cn', file `vars'"

if `have_ppmlhdfe' {
    sysuse auto, clear
    drop if missing(rep78)
    cppmlhdfe price mpg weight, absorb(rep78) ///
        vce(bootstrap, reps(5) seed(3) draws(stata) strata(foreign) saving(`f1', replace double) nodots)
    matrix V1 = e(V)
    local s1 = c(rngstate)
    bootstrap _b, reps(5) seed(3) strata(foreign) nodots saving(`f2', replace double): ///
        ppmlhdfe price mpg weight, absorb(rep78)
    matrix V2 = e(V)
    local s2 = c(rngstate)
    preserve
    use `f1', clear
    mkmat _all, matrix(R1)
    use `f2', clear
    mkmat _all, matrix(R2)
    restore
    if mreldif(R1, R2) < 1e-7 & mreldif(V1, V2) < 1e-7 & "`s1'" == "`s2'" ///
        test_pass "cppmlhdfe saving(): the prefix's replicates and RNG state"
    else test_fail "cppmlhdfe saving()" "reps `=mreldif(R1, R2)' V `=mreldif(V1, V2)'"
}

webuse nlswork, clear
keep if idcode <= 300
civreghdfe ln_wage (ttl_exp = age hours) tenure, absorb(idcode) ///
    vce(bootstrap, reps(8) seed(8) draws(stata) cluster(idcode) saving(`f1', replace double) nodots)
matrix V1 = e(V)
local cn : colnames e(b)
preserve
use `f1', clear
unab vars : _all
mata: st_matrix("Vf", variance(st_data(., .)))
restore
if "`vars'" == "_b_ttl_exp _b_tenure" & "`cn'" == "ttl_exp tenure" & mreldif(V1, Vf) < 1e-12 ///
    test_pass "civreghdfe saving(): e(b)'s terms, variance = e(V)"
else test_fail "civreghdfe saving()" "vars `vars' V `=mreldif(V1, Vf)'"

print_section "draws(stata): strata and fweighted rows"

if `have_reghdfe' {
    * Within strata, bsample as the prefix calls it takes U more uniforms
    * after a replicate's draws when a stratum drew a row twice; the
    * replicates and the final RNG state match from the second replicate on
    * only when they are taken too
    sysuse auto, clear
    drop if missing(rep78)
    creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(8) seed(11) draws(stata) strata(rep78) nodots)
    matrix V1 = e(V)
    local s1 = c(rngstate)
    bootstrap _b, reps(8) seed(11) strata(rep78) nodots: reghdfe price mpg weight, absorb(foreign)
    matrix V2 = e(V)
    local s2 = c(rngstate)
    if mreldif(V1, V2) < 1e-6 & "`s1'" == "`s2'" test_pass "observations within strata: the prefix's replicates and RNG state"
    else test_fail "strata draws(stata)" "V `=mreldif(V1, V2)' rng `=("`s1'" == "`s2'")'"

    * The default is draws(stata): without draws() the results are the
    * prefix's for the same seed
    sysuse auto, clear
    creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(10) seed(14) nodots)
    matrix V1 = e(V)
    local dr `e(draws)'
    local s1 = c(rngstate)
    bootstrap _b, reps(10) seed(14) nodots: reghdfe price mpg weight, absorb(foreign)
    matrix V2 = e(V)
    local s2 = c(rngstate)
    if "`dr'" == "stata" & mreldif(V1, V2) < 1e-6 & "`s1'" == "`s2'" ///
        test_pass "default draws(stata): the prefix's replicates and RNG state without draws()"
    else test_fail "default draws" "draws `dr' V `=mreldif(V1, V2)' rng `=("`s1'" == "`s2'")'"

    * fweighted rows: bsample on the fweight-expanded data, in expand's
    * order; large weights split the replicates into several chunks
    foreach design in "" "strata(rep78)" {
        sysuse auto, clear
        drop if missing(rep78)
        gen fw = cond(mod(_n, 5) == 0, 120, 1 + mod(_n, 3))
        creghdfe price mpg weight [fw=fw], absorb(foreign) vce(bootstrap, reps(12) seed(12) draws(stata) `design' nodots)
        matrix V1 = e(V)
        local s1 = c(rngstate)
        expand fw
        bootstrap _b, reps(12) seed(12) `design' nodots: reghdfe price mpg weight, absorb(foreign)
        matrix V2 = e(V)
        local s2 = c(rngstate)
        if mreldif(V1, V2) < 1e-6 & "`s1'" == "`s2'" ///
            test_pass "creghdfe fweighted rows `design': the prefix on the expanded data"
        else test_fail "creghdfe fweighted rows `design'" "V `=mreldif(V1, V2)' rng `=("`s1'" == "`s2'")'"
    }
}
if `have_ppmlhdfe' {
    sysuse auto, clear
    gen fw = 1 + mod(_n, 3)
    cppmlhdfe price mpg weight [fw=fw], absorb(foreign) vce(bootstrap, reps(6) seed(5) draws(stata) nodots)
    matrix V1 = e(V)
    expand fw
    bootstrap _b, reps(6) seed(5) nodots: ppmlhdfe price mpg weight, absorb(foreign)
    matrix V2 = e(V)
    if mreldif(V1, V2) < 1e-6 test_pass "cppmlhdfe fweighted rows: the prefix on the expanded data"
    else test_fail "cppmlhdfe fweighted rows" "V `=mreldif(V1, V2)'"
}
* civreghdfe (which the prefix cannot run): the same draws through
* bsample on the expanded data, passed as _bsdraws()
webuse nlswork, clear
keep if idcode <= 300
gen w = 1 + mod(idcode, 3)
civreghdfe ln_wage (ttl_exp = age hours) tenure [fw=w], absorb(idcode) vce(bootstrap, reps(4) seed(8) draws(stata) nodots)
matrix V1 = e(V)
keep if e(sample)
expand w
gen long obs = _n
set seed 8
forvalues r = 1/4 {
    sort obs
    gen m`r' = .
    bsample, weight(m`r')
}
sort obs
civreghdfe ln_wage (ttl_exp = age hours) tenure, absorb(idcode) vce(bootstrap, reps(4) _bsdraws(m1-m4) nodots)
if mreldif(V1, e(V)) < 1e-10 test_pass "civreghdfe fweighted rows: bsample on the expanded data"
else test_fail "civreghdfe fweighted rows" "V `=mreldif(V1, e(V))'"

print_section "_bsdraws()"

* The prefix's cluster draws (bsref_draws) give draws(stata)'s replicates
sysuse auto, clear
gen cl = mod(_n * 7, 11)
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(10) seed(13) draws(stata) cluster(cl) nodots)
matrix V1 = e(V)
set seed 13
bsref_draws, reps(10) generate(m) cluster(cl)
creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(10) _bsdraws(m1-m10) cluster(cl) nodots)
local dr `e(draws)'
if mreldif(V1, e(V)) < 1e-12 & "`dr'" == "supplied" test_pass "_bsdraws(): the multiplicities of draws(stata)"
else test_fail "_bsdraws()" "V `=mreldif(V1, e(V))' draws `dr'"
replace m1 = m1 + 1 in 1
capture creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(10) _bsdraws(m1-m10) cluster(cl) nodots)
local rc1 = _rc
replace m1 = -1 in 1
capture creghdfe price mpg weight, absorb(foreign) vce(bootstrap, reps(10) _bsdraws(m1-m10) nodots)
if `rc1' == 198 & _rc == 198 test_pass "_bsdraws(): a cluster's rows must agree, multiplicities nonnegative"
else test_fail "_bsdraws() checks" "rc `rc1' `=_rc'"

print_section "fevce(bootstrap, ci)"

if `have_reghdfe' {
    sysuse auto, clear
    drop if missing(rep78)
    gen cl = mod(_n * 7, 11)
    set seed 21
    bsref_draws, reps(30) generate(m) cluster(cl)
    creghdfe price mpg weight, absorb(rep78, savefe) ///
        fevce(bootstrap, reps(30) _bsdraws(m1-m30) cluster(cl) ci level(90) nodots)
    local names `e(fe_lb)' `e(fe_ub)'
    local lev = e(fe_level)
    gen long obs = _n
    forvalues r = 1/30 {
        preserve
        quietly keep if m`r' > 0
        quietly expand m`r'
        quietly reghdfe price mpg weight, absorb(f=rep78) tolerance(1e-12)
        quietly collapse (firstnm) f, by(rep78)
        rename f f_`r'
        tempfile lv
        quietly save `lv'
        restore
        quietly merge m:1 rep78 using `lv', nogenerate
    }
    sort obs
    _bo_pctile, stub(f_) reps(30) level(90)
    gen double d1 = abs(o_lb - __felb1__) / (1 + abs(o_lb))
    gen double d2 = abs(o_ub - __feub1__) / (1 + abs(o_ub))
    sum d1
    local m1 = r(max)
    sum d2
    local m2 = r(max)
    if "`names'" == "__felb1__ __feub1__" & `lev' == 90 & `m1' < 1e-6 & `m2' < 1e-6 ///
        test_pass "creghdfe ci: percentile bounds of reghdfe's effects on the same draws"
    else test_fail "creghdfe ci" "names `names' level `lev' lb `m1' ub `m2'"
}
if `have_ppmlhdfe' {
    _bo_pf_data
    quietly cppmlhdfe y x1, absorb(id)
    keep if e(sample)
    set seed 22
    bsref_draws, reps(20) generate(m) cluster(cl)
    cppmlhdfe y x1, absorb(fid=id) fevce(bootstrap, reps(20) _bsdraws(m1-m20) cluster(cl) ci nodots)
    local names `e(fe_lb)' `e(fe_ub)'
    gen long obs = _n
    forvalues r = 1/20 {
        preserve
        quietly keep if m`r' > 0
        quietly ppmlhdfe y x1 [fw=m`r'], absorb(f=id) tolerance(1e-12)
        quietly collapse (firstnm) f, by(id)
        rename f f_`r'
        tempfile lv
        quietly save `lv'
        restore
        quietly merge m:1 id using `lv', nogenerate
    }
    sort obs
    _bo_pctile, stub(f_) reps(20)
    gen double d1 = abs(o_lb - fid_lb) / (1 + abs(o_lb))
    gen double d2 = abs(o_ub - fid_ub) / (1 + abs(o_ub))
    sum d1
    local m1 = r(max)
    sum d2
    local m2 = r(max)
    if "`names'" == "fid_lb fid_ub" & `m1' < 1e-5 & `m2' < 1e-5 ///
        test_pass "cppmlhdfe ci: percentile bounds of ppmlhdfe's effects on the same draws"
    else test_fail "cppmlhdfe ci" "names `names' lb `m1' ub `m2'"
}
* civreghdfe: its own savefe on each replicate (fweights)
webuse nlswork, clear
keep if idcode <= 200
quietly civreghdfe ln_wage (ttl_exp = age hours) tenure, absorb(idcode)
keep if e(sample)
keep ln_wage ttl_exp age hours tenure idcode
gen grp = ceil(idcode / 10)
set seed 23
bsref_draws, reps(15) generate(m) cluster(grp)
civreghdfe ln_wage (ttl_exp = age hours) tenure, absorb(idcode, savefe) ///
    fevce(bootstrap, reps(15) _bsdraws(m1-m15) cluster(grp) ci sestub(cs) nodots)
local names `e(fe_lb)' `e(fe_ub)'
rename __hdfe1__ fe_main
gen long obs = _n
forvalues r = 1/15 {
    preserve
    quietly keep if m`r' > 0
    quietly civreghdfe ln_wage (ttl_exp = age hours) tenure [fw=m`r'], absorb(f=idcode)
    quietly collapse (firstnm) f, by(idcode)
    rename f f_`r'
    tempfile lv
    quietly save `lv'
    restore
    quietly merge m:1 idcode using `lv', nogenerate
}
sort obs
_bo_pctile, stub(f_) reps(15)
gen double d1 = abs(o_lb - cs1_lb) / (1 + abs(o_lb))
gen double d2 = abs(o_ub - cs1_ub) / (1 + abs(o_ub))
sum d1
local m1 = r(max)
sum d2
local m2 = r(max)
if "`names'" == "cs1_lb cs1_ub" & `m1' < 1e-6 & `m2' < 1e-6 ///
    test_pass "civreghdfe ci: percentile bounds of its own savefe on the same draws"
else test_fail "civreghdfe ci" "names `names' lb `m1' ub `m2'"

print_section "Heterogeneous slopes"

if `have_reghdfe' {
    webuse nlswork, clear
    keep if idcode <= 120
    keep ln_wage age ttl_exp tenure idcode year
    drop if missing(ln_wage, age, ttl_exp, tenure)
    quietly creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year)
    keep if e(sample)
    gen long obs = _n
    set seed 77
    bsref_draws, reps(6) generate(m) cluster(idcode)
    creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) vce(bootstrap, reps(6) _bsdraws(m1-m6) cluster(idcode) nodots)
    matrix Vc = e(V)
    matrix B = J(6, 3, .)
    forvalues r = 1/6 {
        quietly reghdfe ln_wage age tenure [fw=m`r'] if m`r' > 0, absorb(idcode##c.ttl_exp year) tolerance(1e-12)
        matrix B[`r', 1] = e(b)
    }
    mata: st_matrix("Vo", variance(st_matrix("B")))
    if mreldif(Vc, Vo) < 1e-8 test_pass "creghdfe slopes: e(V) from reghdfe on the same draws"
    else test_fail "creghdfe slopes V" "`=mreldif(Vc, Vo)'"

    * One effect with a slope: per-level intercept and slope SEs and d's SE
    quietly creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp)
    keep if e(sample)
    creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp, savefe) ///
        fevce(bootstrap, reps(6) _bsdraws(m1-m6) cluster(idcode) dse(dse) nodots)
    local names `e(fe_se)'
    local warn = e(fe_norm_warning)
    forvalues r = 1/6 {
        quietly reghdfe ln_wage age tenure [fw=m`r'] if m`r' > 0, absorb(a`r'=idcode##c.ttl_exp) tolerance(1e-12)
        gen double d`r' = a`r' + a`r'Slope1 * ttl_exp
    }
    egen double o_a = rowsd(a1 a2 a3 a4 a5 a6)
    egen double o_s = rowsd(a1Slope1 a2Slope1 a3Slope1 a4Slope1 a5Slope1 a6Slope1)
    egen double o_d = rowsd(d1 d2 d3 d4 d5 d6)
    local worst = 0
    foreach p in "o_a __fese1__" "o_s __fese1__Slope1" "o_d dse" {
        gettoken o c : p
        gen double dd = abs(`o' - `c') / (1 + abs(`o'))
        sum dd
        local worst = max(`worst', r(max))
        drop dd
    }
    if "`names'" == "__fese1__ __fese1__Slope1" & `warn' == 0 & `worst' < 1e-6 ///
        test_pass "creghdfe slopes: intercept, slope and d SEs from reghdfe on the same draws"
    else test_fail "creghdfe slopes fevce" "names `names' warn `warn' worst `worst'"

    * Slopes alone (no constant)
    creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(bootstrap, reps(6) _bsdraws(m1-m6) cluster(idcode) nodots)
    matrix Vc = e(V)
    local cn : colnames e(b)
    matrix B = J(6, 2, .)
    forvalues r = 1/6 {
        quietly reghdfe ln_wage age tenure [fw=m`r'] if m`r' > 0, absorb(idcode#c.ttl_exp) tolerance(1e-12)
        matrix B[`r', 1] = e(b)
    }
    mata: st_matrix("Vo", variance(st_matrix("B")))
    if "`cn'" == "age tenure" & mreldif(Vc, Vo) < 1e-8 test_pass "creghdfe slopes alone: no constant, e(V) from reghdfe"
    else test_fail "creghdfe slopes alone" "`cn' `=mreldif(Vc, Vo)'"

    * Slopes alone with fevce(): the replicates subtract no constant before
    * recovering the slopes, as the full fit (BOOTSTRAP_REVIEW.md, H3)
    capture drop __hdfe*
    capture drop __fese*
    creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp, savefe) ///
        fevce(bootstrap, reps(6) _bsdraws(m1-m6) cluster(idcode) dse(dse0) nodots)
    local se1 : word 1 of `e(fe_se)'
    forvalues r = 1/6 {
        quietly reghdfe ln_wage age tenure [fw=m`r'] if m`r' > 0, absorb(b`r'=idcode#c.ttl_exp) tolerance(1e-12)
        gen double e`r' = b`r'Slope1 * ttl_exp
    }
    egen double o_b = rowsd(b1Slope1 b2Slope1 b3Slope1 b4Slope1 b5Slope1 b6Slope1)
    egen double o_e = rowsd(e1 e2 e3 e4 e5 e6)
    local worst = 0
    foreach p in "o_b `se1'" "o_e dse0" {
        gettoken o c : p
        gen double dd = abs(`o' - `c') / (1 + abs(`o'))
        sum dd
        local worst = max(`worst', r(max))
        drop dd
    }
    if "`se1'" != "" & `worst' < 1e-6 test_pass "creghdfe slopes alone: slope and d SEs from reghdfe on the same draws"
    else test_fail "creghdfe slopes alone fevce" "`se1' worst `worst'"
}

* Wild bootstrap: absorbing the slopes or adding them as regressors gives
* the same restricted model, so the same p-values and intervals (nested
* and projected paths)
webuse nlswork, clear
keep if idcode <= 120
drop if missing(ln_wage, age, ttl_exp, tenure)
local a1 idcode##c.ttl_exp
local x1 i.idcode#c.ttl_exp
local e1 idcode
local a2 year##c.ttl_exp
local x2 i.year#c.ttl_exp
local e2 year
forvalues j = 1/2 {
    creghdfe ln_wage age tenure, absorb(`a`j'') vce(cluster idcode) wildbootstrap(reps(199) seed(3) coefs(age) nodots)
    matrix P1 = e(wb_p), e(wb_ci)
    creghdfe ln_wage age tenure `x`j'', absorb(`e`j'') vce(cluster idcode) wildbootstrap(reps(199) seed(3) coefs(age) nodots)
    matrix P2 = e(wb_p), e(wb_ci)
    if mreldif(P1, P2) < 1e-7 test_pass "creghdfe wild, absorb(`a`j''): as the slopes as regressors"
    else test_fail "creghdfe wild `a`j''" "`=mreldif(P1, P2)'"
}

if `have_ppmlhdfe' {
    clear
    set seed 31
    set obs 600
    gen id = ceil(_n / 6)
    gen t = mod(_n - 1, 6) + 1
    gen z = rnormal()
    gen x1 = rnormal()
    gen x2 = rnormal() + id / 100
    gen double y = rpoisson(exp(0.03 * x1 - 0.2 * x2 + mod(id, 7) / 5 + (mod(id, 3) - 1) * z / 4 + t / 10))
    gen cl = ceil(id / 4)
    quietly cppmlhdfe y x1 x2, absorb(id##c.z)
    keep if e(sample)
    set seed 9
    bsref_draws, reps(5) generate(m) cluster(cl)
    cppmlhdfe y x1 x2, absorb(id##c.z, savefe) d(dd) vce(bootstrap, reps(5) _bsdraws(m1-m5) cluster(cl) nodots) fevce(bootstrap)
    matrix Vc = e(V)
    local names `e(fe_se)'
    matrix B = J(5, 3, .)
    forvalues r = 1/5 {
        quietly ppmlhdfe y x1 x2 [fw=m`r'] if m`r' > 0, absorb(a`r'=id##c.z) tolerance(1e-12)
        matrix B[`r', 1] = e(b)
        gen double d`r' = a`r' + a`r'Slope1 * z
    }
    mata: st_matrix("Vo", variance(st_matrix("B")))
    egen double o_a = rowsd(a1 a2 a3 a4 a5)
    egen double o_s = rowsd(a1Slope1 a2Slope1 a3Slope1 a4Slope1 a5Slope1)
    egen double o_d = rowsd(d1 d2 d3 d4 d5)
    local worst = 0
    foreach p in "o_a __fese1__" "o_s __fese1__Slope1" "o_d dd_se" {
        gettoken o c : p
        gen double dd2 = abs(`o' - `c') / (1 + abs(`o'))
        sum dd2
        local worst = max(`worst', r(max))
        drop dd2
    }
    if "`names'" == "__fese1__ __fese1__Slope1" & mreldif(Vc, Vo) < 1e-6 & `worst' < 1e-4 ///
        test_pass "cppmlhdfe slopes: e(V), intercept, slope and d SEs from ppmlhdfe on the same draws"
    else test_fail "cppmlhdfe slopes" "names `names' V `=mreldif(Vc, Vo)' worst `worst'"
    cppmlhdfe y x1 x2, absorb(id##c.z) vce(cluster cl) wildbootstrap(reps(199) seed(3) coefs(x1) nodots)
    matrix P1 = e(wb_p)
    cppmlhdfe y x1 x2 i.id#c.z, absorb(id) vce(cluster cl) wildbootstrap(reps(199) seed(3) coefs(x1) nodots)
    matrix P2 = e(wb_p)
    if mreldif(P1, P2) < 1e-7 & P1[1, 1] > 0.01 test_pass "cppmlhdfe score bootstrap with slopes: as the slopes as regressors"
    else test_fail "cppmlhdfe wild slopes" "p `=P1[1, 1]' vs `=P2[1, 1]'"
}

print_section "civreghdfe: savefe and fevce() with partial()"

webuse nlswork, clear
keep if idcode <= 400
civreghdfe ln_wage (ttl_exp = hours wks_ue) tenure age, absorb(idcode year, savefe) ///
    fevce(bootstrap, reps(10) seed(3) sestub(s0) nodots)
rename (__hdfe1__ __hdfe2__) (a1 a2)
civreghdfe ln_wage (ttl_exp = hours wks_ue) tenure age, absorb(idcode year, savefe) partial(age) ///
    fevce(bootstrap, reps(10) seed(3) sestub(s1) nodots)
local worst = 0
foreach p in "a1 __hdfe1__" "a2 __hdfe2__" "s01 s11" "s02 s12" {
    gettoken o c : p
    gen double dd = abs(`o' - `c') / (1 + abs(`o'))
    sum dd
    local worst = max(`worst', r(max))
    drop dd
}
if `worst' < 1e-9 test_pass "partial(): the effects and their bootstrap SEs of the model without it"
else test_fail "civreghdfe partial() savefe" "worst `worst'"

* A partial() variable the effects absorb gets no coefficient in the
* replicates' effects either (BOOTSTRAP_REVIEW.md, H2). Constant within
* year, the second effect, it is partialled to rounding noise rather than
* exactly zero (a variable constant within the first effect is exactly
* zero after the first sweep)
webuse nlswork, clear
keep if idcode <= 400
gen byte grp = mod(year, 3) == 0
civreghdfe ln_wage (ttl_exp = hours wks_ue) tenure age, absorb(idcode year, savefe) ///
    fevce(bootstrap, reps(10) seed(3) sestub(t0) nodots)
rename (__hdfe1__ __hdfe2__) (b1 b2)
civreghdfe ln_wage (ttl_exp = hours wks_ue) tenure age grp, absorb(idcode year, savefe) partial(grp) ///
    fevce(bootstrap, reps(10) seed(3) sestub(t1) nodots)
local worst = 0
foreach p in "b1 __hdfe1__" "b2 __hdfe2__" "t01 t11" "t02 t12" {
    gettoken o c : p
    gen double dd = abs(`o' - `c') / (1 + abs(`o'))
    sum dd
    local worst = max(`worst', r(max))
    drop dd
}
if `worst' < 1e-9 test_pass "partial() of an absorbed variable: the effects and SEs of the model without it"
else test_fail "civreghdfe absorbed partial() fevce" "worst `worst'"

print_section "creghdfe: fevce() after the singleton refit, and its options"

* A factor level only on singleton rows makes creghdfe refit on the final
* sample; fevce()'s outputs are created again (BOOTSTRAP_REVIEW.md, M3)
clear
set seed 12
set obs 200
gen id = ceil(_n / 5)
gen x = rnormal()
gen cat = mod(_n, 3)
replace id = 1000 in 200
replace cat = 3 in 200
gen y = x + id / 50 + rnormal()
capture noisily creghdfe y x i.cat, absorb(id, savefe) fevce(bootstrap, reps(5) seed(1) nodots)
local rc = _rc
capture confirm variable __fese1__
if !`rc' & !_rc & e(fe_N_reps) == 5 test_pass "fevce(bootstrap) after the singleton refit"
else test_fail "creghdfe refit fevce" "rc `rc'"

* fevce() may repeat vce(bootstrap)'s level(), mse and dots(), not change
* them (BOOTSTRAP_REVIEW.md, M7)
capture drop __hdfe*
capture drop __fese*
capture creghdfe y x, absorb(id, savefe) vce(bootstrap, reps(5) seed(1) level(90) nodots) fevce(bootstrap, ci level(95))
local rc1 = _rc
capture drop __hdfe*
capture drop __fese*
capture drop __felb*
capture drop __feub*
capture creghdfe y x, absorb(id, savefe) vce(bootstrap, reps(5) seed(1) level(90) nodots) fevce(bootstrap, ci level(90))
local rc2 = _rc
local lev = e(fe_level)
capture drop __hdfe*
capture drop __fese*
capture drop __felb*
capture drop __feub*
capture creghdfe y x, absorb(id, savefe) vce(bootstrap, reps(5) seed(1) nodots) fevce(bootstrap, mse)
local rc3 = _rc
capture drop __hdfe*
capture drop __fese*
capture creghdfe y x, absorb(id, savefe) vce(bootstrap, reps(5) seed(1) level(90) nodots) fevce(bootstrap)
local rc4 = _rc
if `rc1' == 198 & `rc2' == 0 & `lev' == 90 & `rc3' == 198 & `rc4' == 0 ///
    test_pass "fevce(bootstrap) with vce(bootstrap): level() and mse given once"
else test_fail "fevce() options with vce(bootstrap)" "rc `rc1' `rc2' `rc3' `rc4' level `lev'"

if `have_reghdfe' {
    print_section "draws(stata): strata of two rows"

    * bsample takes its extra uniforms within strata only when a stratum
    * drew one of its rows twice; with three strata of two rows, about one
    * replicate in eight draws none twice (BOOTSTRAP_REVIEW.md, M4)
    clear
    set obs 6
    gen st = ceil(_n / 2)
    gen x = _n
    gen y = 2 * x + mod(_n * 7, 5)
    gen g = 1
    creghdfe y x, absorb(g) vce(bootstrap, reps(50) seed(11) draws(stata) strata(st) nodots)
    matrix V1 = e(V)
    local s1 `"`c(rngstate)'"'
    bootstrap _b, reps(50) seed(11) strata(st) nodots: reghdfe y x, absorb(g)
    matrix V2 = e(V)
    local s2 `"`c(rngstate)'"'
    local same = (`"`s1'"' == `"`s2'"')
    if mreldif(V1, V2) < 1e-6 & `same' test_pass "draws(stata), strata of two rows: the prefix's replicates and RNG state"
    else test_fail "tiny strata draws(stata)" "V `=mreldif(V1, V2)' rng `same'"
}

print_section "Results as bootstrap posts them (BOOTSTRAP_REVIEW.md, L1, L2, L9)"

* mse: e(vcetype) Bstrap * and e(mse), as _bs_sum; the level() of the
* command sets the bootstrap's, as for bootstrap
webuse nlswork, clear
keep if idcode <= 200
local bad
foreach cmd in creghdfe cppmlhdfe civreghdfe {
    if "`cmd'" == "civreghdfe" local spec ln_wage (tenure = union) age
    else if "`cmd'" == "cppmlhdfe" local spec hours age tenure
    else local spec ln_wage age tenure
    quietly `cmd' `spec', absorb(idcode) vce(bootstrap, reps(5) seed(2) mse nodots)
    if "`e(vcetype)'" != "Bstrap *" | "`e(mse)'" != "mse" local bad `bad' `cmd':mse
    quietly `cmd' `spec', absorb(idcode) vce(bootstrap, reps(5) seed(2) level(90) nodots)
    matrix C1 = e(ci_normal)
    local l1 = e(level)
    quietly `cmd' `spec', absorb(idcode) vce(bootstrap, reps(5) seed(2) nodots) level(90)
    matrix C2 = e(ci_normal)
    if `l1' != 90 | e(level) != 90 | mreldif(C1, C2) > 1e-12 local bad `bad' `cmd':level
    capture `cmd' `spec', absorb(idcode) vce(bootstrap, reps(5) seed(2) level(80) nodots) level(90)
    if _rc != 198 local bad `bad' `cmd':conflict
    capture `cmd' `spec', absorb(idcode) vce(bootstrap, reps(5) seed(2) level(90) nodots) level(90)
    if _rc local bad `bad' `cmd':same
    if "`cmd'" == "cppmlhdfe" {
        quietly cppmlhdfe `spec', absorb(idcode) vce(bootstrap, reps(5) seed(2) nodots)
        if reldif(e(p), chi2tail(e(df_m), e(chi2))) > 1e-12 | missing(e(p)) local bad `bad' cppmlhdfe:p
    }
}
if "`bad'" == "" test_pass "mse labels, e(p), and level() as bootstrap treats them"
else test_fail "bootstrap e() parity" "`bad'"

* cppmlhdfe: d() alone has no normalization caveat
quietly cppmlhdfe hours age tenure, absorb(idcode year) d(dd) fevce(bootstrap, reps(4) seed(1) nodots)
if e(fe_norm_warning) == 0 test_pass "cppmlhdfe fevce() with d() alone: no normalization warning"
else test_fail "cppmlhdfe d() warning" "fe_norm_warning `=e(fe_norm_warning)'"
capture drop dd dd_se

* wildbootstrap(coefs()) may name a regressor noci (BOOTSTRAP_REVIEW.md, L8)
gen double noci = age
quietly creghdfe ln_wage noci tenure, absorb(idcode) vce(cluster idcode) ///
    wildbootstrap(reps(19) seed(1) coefs(noci) nodots)
if "`e(wb_coefs)'" == "noci" & colsof(e(wb_p)) == 1 test_pass "wildbootstrap(coefs(noci)): a regressor named noci"
else test_fail "coefs(noci)" "coefs `e(wb_coefs)'"
drop noci

* civreghdfe: partial() with more than 64 columns and fevce() (L3): the
* effects and SEs of the model with them as regressors
clear
set seed 4
set obs 3000
gen id = ceil(_n / 30)
gen g = 1 + mod(_n * 7, 70)
forvalues j = 1/70 {
    gen byte p`j' = g == `j'
}
gen double z = rnormal()
gen double x = z + rnormal()
gen double y = x + g / 70 + id / 100 + rnormal()
civreghdfe y (x = z) p2-p70, absorb(id, savefe) fevce(bootstrap, reps(4) seed(3) sestub(q0) nodots)
rename __hdfe1__ c1
capture noisily civreghdfe y (x = z) p2-p70, absorb(id, savefe) partial(p2-p70) ///
    fevce(bootstrap, reps(4) seed(3) sestub(q1) nodots)
local rc = _rc
local worst = .
if !`rc' {
    local worst = 0
    foreach p in "c1 __hdfe1__" "q01 q11" {
        gettoken o c : p
        gen double dd = abs(`o' - `c') / (1 + abs(`o'))
        sum dd
        local worst = max(`worst', r(max))
        drop dd
    }
}
if !`rc' & `worst' < 1e-8 test_pass "partial() of 69 columns with fevce(): the effects and SEs of the model without it"
else test_fail "civreghdfe partial() > 64 fevce" "rc `rc' worst `worst'"

noi print_summary "bootstrap_options"
}

global CTOOLS_COMPONENT_COMPLETE "bootstrap_options"
    global CTOOLS_CASE_COMPLETE "bootstrap_options"
    exit
}

di as error "Unknown bootstrap_reference case: `__case'"
exit 198
