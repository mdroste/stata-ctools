* s2_psm: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_psm"

if "`__case'" == "s2_psm" {
do validation/support/validate_setup.do
set seed 20261104

/* _s2psm_cmp runs psmatch2 and then cpsmatch with the same arguments (the
   command line follows the test name and a colon) on the same data and sort
   order. It compares the return code, every r() scalar and macro psmatch2
   returns, and every variable psmatch2 creates (values, variable labels with
   cpsmatch: for psmatch2:, the storage type of _id and _n#), plus
   char _weight[Type]. Values agree within 1e-9 (reldif); both missing is
   agreement. */
capture program drop _s2psm_cmp
program define _s2psm_cmp
    gettoken name 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    local name = strtrim(`"`name'"')
    local tol 1e-9
    capture drop _*
    capture quietly psmatch2 `0'
    local rc_ref = _rc
    local scalars : r(scalars)
    foreach s of local scalars {
        tempname R_`s'
        scalar `R_`s'' = r(`s')
    }
    local macros depvar exog
    foreach m of local macros {
        local M_`m' `"`r(`m')'"'
    }
    local vars
    if !`rc_ref' {
        quietly ds _*
        local vars `r(varlist)'
        foreach v of local vars {
            tempvar REF`v'
            quietly gen double `REF`v'' = `v'
            local lab_`v' : variable label `v'
            local typ_`v' : type `v'
        }
        local char_ref : char _weight[Type]
    }
    * (not drop _*, which would match the temporary variables)
    capture drop `vars'
    capture quietly cpsmatch `0'
    local rc_ct = _rc
    local bad
    if `rc_ref' != `rc_ct' local bad rc:`rc_ref'/`rc_ct'
    else if !`rc_ref' {
        foreach s of local scalars {
            local a = scalar(`R_`s'')
            if !(`a' == r(`s') | reldif(`a', r(`s')) <= `tol') local bad `bad' r(`s'):`a'/`=r(`s')'
        }
        foreach m of local macros {
            if `"`r(`m')'"' != `"`M_`m''"' local bad `bad' r(`m')
        }
        foreach v of local vars {
            capture confirm variable `v'
            if _rc {
                local bad `bad' `v':absent
                continue
            }
            quietly count if !(`v' == `REF`v'' | reldif(`v', `REF`v'') <= `tol')
            if r(N) local bad `bad' `v':`r(N)'
            local lab : variable label `v'
            if subinstr(`"`lab_`v''"', "psmatch2:", "cpsmatch:", 1) != `"`lab'"' local bad `bad' `v':label
            if "`v'" == "_id" | regexm("`v'", "^_n[0-9]+$") {
                local t : type `v'
                if "`t'" != "`typ_`v''" local bad `bad' `v':type
            }
        }
        local char_ct : char _weight[Type]
        if `"`char_ct'"' != `"`char_ref'"' local bad `bad' char
    }
    if "`bad'" == "" & `rc_ref' == 0 test_pass "`name'"
    else if "`bad'" == "" test_pass "`name': same error r(`rc_ref')"
    else test_fail "`name'" "`bad'"
end

* Continuous covariates, a binary covariate and a tied (rounded) score; the
* data in a fixed random order. y2 misses a value (the outcomes' common
* sample), half selects a subsample.
clear
quietly set obs 400
gen double x1 = rnormal()
gen double x2 = rnormal()
gen byte xb = runiform() < .4
gen byte t = (x1 + .5*xb + rnormal() > 0)
gen double y = x1 + t + rnormal()
gen double y2 = x2 + 2*t + rnormal()
quietly replace y2 = . in 7
gen double ps = round(normal(x1/2), .02)
gen byte half = mod(_n, 3) > 0
gen double u = runiform()
sort u

* PSM-12/13/14: nearest neighbor k = 1, 2, 3, with and without ties, on a
* continuous and a tied score: _id, _n#, _nn, _pdif, _<outcome>, labels
foreach o in "" "neighbor(2)" "neighbor(3)" "ties" "neighbor(2) ties" "neighbor(3) ties" ///
    "caliper(.005)" "neighbor(3) caliper(.01)" {
    _s2psm_cmp PSM-12 nearest `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-12 nearest `o', tied score: t, pscore(ps) outcome(y y2) `o'
    _s2psm_cmp PSM-12 nearest `o', if: t x1 x2 if half, outcome(y) `o'
}
_s2psm_cmp PSM-12 nearest, no outcome: t x1 x2
_s2psm_cmp PSM-12 nearest, in range: t x1 x2 in 20/350, outcome(y2) neighbor(2)

* PSM-11: without replacement (1-to-1 only, as psmatch2), descending
foreach o in "noreplacement" "noreplacement ties" "noreplacement descending" ///
    "noreplacement descending ties" "noreplacement caliper(.01)" {
    _s2psm_cmp PSM-11 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-11 `o', tied score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-11 noreplacement neighbor(2) rejected: t x1 x2, outcome(y) noreplacement neighbor(2)
_s2psm_cmp PSM-11 noreplacement kernel rejected: t x1 x2, outcome(y) noreplacement kernel
_s2psm_cmp PSM-11 noreplacement mahalanobis rejected: t, mahalanobis(x1 x2) outcome(y) noreplacement
_s2psm_cmp PSM-11 descending without noreplacement: t x1 x2, outcome(y) descending
* psmatch2 accepts noreplacement with radius, which ignores it
_s2psm_cmp PSM-11 noreplacement radius: t x1 x2, outcome(y) noreplacement radius caliper(.02)
* Treated outnumbering the controls: without replacement the controls run out
gen byte t2 = (x1 + rnormal() > -1)
_s2psm_cmp PSM-11 controls run out: t2, pscore(ps) outcome(y) noreplacement
_s2psm_cmp PSM-11 controls run out, descending: t2 x1 x2, outcome(y) noreplacement descending

* PSM-02: ate (ATU, ATE, their r() and _<outcome>, _weight, _support)
foreach o in "ate" "ate neighbor(3)" "ate ties neighbor(2)" "ate caliper(.005)" "ate noreplacement" ///
    "ate noreplacement descending" "ate common" "ate kernel" "ate kernel kerneltype(normal)" ///
    "ate radius caliper(.02)" "ate llr" "ate llr kerneltype(biweight)" {
    _s2psm_cmp PSM-02 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-02 `o', tied score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-02 ate, many neighbors: t2 x1 x2, outcome(y) ate neighbor(10)
_s2psm_cmp PSM-02 ate, no outcome: t x1 x2, ate

* PSM-16: radius without caliper() (or with caliper(0)) matches every control,
* as psmatch2's missing caliper does
foreach o in "radius" "radius caliper(0)" "radius common" "radius ate" "radius mahalanobis(x1 x2)" {
    _s2psm_cmp PSM-16 `o', estimated: t x1 x2, outcome(y y2) `o'
    _s2psm_cmp PSM-16 `o', tied score: t, pscore(ps) outcome(y) `o'
}

* Radius and kernel matching, each kernel type, and common support
foreach o in "radius caliper(.01)" "radius caliper(.05)" "kernel" "kernel kerneltype(normal)" ///
    "kernel kerneltype(biweight)" "kernel kerneltype(uniform)" "kernel kerneltype(tricube)" ///
    "kernel kerneltype(epan) bwidth(.1)" "kernel kerneltype(normal) bwidth(.02)" ///
    "common" "common neighbor(2)" "common kernel" "common radius caliper(.02)" {
    _s2psm_cmp radius/kernel `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp radius/kernel `o', tied score: t, pscore(ps) outcome(y) `o'
}

* PSM-04: llr, the Epanechnikov kernel (lpoly-smoothed outcomes matched to
* the nearest neighbor, _s_<outcome>) and the other kernels
foreach o in "llr" "llr bwidth(.1)" "llr kerneltype(normal)" "llr kerneltype(biweight)" ///
    "llr kerneltype(uniform) bwidth(.1)" "llr kerneltype(tricube)" "llr neighbor(2)" {
    _s2psm_cmp PSM-04 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-04 `o', tied score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-04 llr without outcome: t x1 x2, llr
_s2psm_cmp PSM-04 llr with kernel rejected: t x1 x2, outcome(y) llr kernel

* PSM-03: mahalanobis() (nearest neighbor, kernel, llr, radius) and w()
foreach o in "" "neighbor(3)" "caliper(.5)" "ate" "ate neighbor(2) caliper(.5)" "kernel" ///
    "kernel kerneltype(normal) bwidth(.5)" "kernel ate bwidth(.5)" "llr bwidth(.5)" ///
    "llr ate kerneltype(normal) bwidth(1)" "radius caliper(.3)" "ties" {
    _s2psm_cmp PSM-03 mahalanobis `o': t, mahalanobis(x1 x2) outcome(y y2) `o'
    _s2psm_cmp PSM-03 mahalanobis `o', with score: t x1 x2, mahalanobis(x1 xb) outcome(y) `o'
}
_s2psm_cmp PSM-03 mahalanobis common: t x1, mahalanobis(x1 x2) outcome(y) common
_s2psm_cmp PSM-03 mahalanobis common without score: t, mahalanobis(x1 x2) outcome(y) common
matrix S2W = (2, .5 \ .5, 1)
matrix colnames S2W = x1 x2
matrix rownames S2W = x1 x2
_s2psm_cmp PSM-03 w(): t, mahalanobis(x1 x2) w(S2W) outcome(y) neighbor(2)
_s2psm_cmp PSM-03 w() kernel: t, mahalanobis(x1 x2) w(S2W) outcome(y) kernel bwidth(.5)
_s2psm_cmp PSM-03 w() ate: t, mahalanobis(x1 x2) w(S2W) outcome(y) ate

* PSM-06: ai(#) with population and altvariance (Abadie-Imbens variance),
* on the score and on the Mahalanobis metric. psmatch2's altvariance fails
* with several outcomes, so it has one.
foreach o in "ai(1)" "ai(2)" "ai(2) ate" "ai(2) population" "ai(3) population ate neighbor(2)" ///
    "ai(2) ties" "ai(2) caliper(.01) ate" "ai(2) noreplacement" {
    _s2psm_cmp PSM-06 `o', score: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-06 `o', mahalanobis: t, mahalanobis(x1 x2) outcome(y y2) `o'
}
foreach o in "ai(2) altvariance" "ai(3) altvariance ate" "ai(2) altvariance population ate" ///
    "ai(2) altvariance ties" {
    _s2psm_cmp PSM-06 `o', score: t, pscore(ps) outcome(y) `o'
    _s2psm_cmp PSM-06 `o', mahalanobis: t, mahalanobis(x1 x2) outcome(y) `o'
}
_s2psm_cmp PSM-06 ai with kernel rejected: t x1 x2, outcome(y) ai(2) kernel
_s2psm_cmp PSM-06 ai with radius: t x1 x2, outcome(y) ai(2) radius caliper(.1)

* PSM-07: trim(#), with ate and with the other methods
foreach o in "trim(5)" "trim(0)" "trim(20) ate" "trim(10) kernel" "trim(10) common neighbor(2)" ///
    "trim(5) ate ai(2)" {
    _s2psm_cmp PSM-07 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-07 `o', given score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-07 trim out of range: t, pscore(ps) outcome(y) trim(-5)
_s2psm_cmp PSM-07 trim above 100 with common: t, pscore(ps) outcome(y) trim(150) common
_s2psm_cmp PSM-07 trim without a score: t, mahalanobis(x1) outcome(y) trim(5)

* S3: odds, index, quietly, warnings, logit
foreach o in "odds" "odds logit" "odds logit trim(5) ate" "index" "index common" "quietly" ///
    "quietly logit kernel" "warnings" "warnings ate" "warnings neighbor(2) ties" {
    _s2psm_cmp S3 `o': t x1 x2 xb, outcome(y y2) `o'
}
_s2psm_cmp S3 warnings, tied score: t, pscore(ps) outcome(y) warnings
_s2psm_cmp S3 warnings, tied score, ate: t, pscore(ps) outcome(y) warnings ate

* PSM-05: spline needs the -spline- package, as in psmatch2
capture which spline
if _rc {
    _s2psm_cmp PSM-05 spline without the spline package: t x1 x2, outcome(y) spline
}
else {
    _s2psm_cmp PSM-05 spline: t x1 x2, outcome(y) spline
    _s2psm_cmp PSM-05 spline nknots(): t x1 x2, outcome(y) spline nknots(3) ate
}

* Other treatment codes and one-arm samples
gen byte t3 = cond(mod(_n, 5) == 0, 2, t)
foreach o in "" "ate" "kernel ate" "ai(2)" "neighbor(2) ties" {
    _s2psm_cmp other codes `o': t3, pscore(ps) outcome(y) `o'
}
_s2psm_cmp other codes, mahalanobis ate: t3, mahalanobis(x1 x2) outcome(y) ate
gen byte t0 = 0
gen byte t1 = 1
_s2psm_cmp no treated, ate: t0, pscore(ps) outcome(y) ate
_s2psm_cmp no controls, ate kernel: t1, pscore(ps) outcome(y) ate kernel

* Option checks in psmatch2's order
_s2psm_cmp checks: no score: t, outcome(y)
_s2psm_cmp checks: score and covariates: t x1, pscore(ps) outcome(y)
_s2psm_cmp checks: kerneltype: t, pscore(ps) outcome(y) kernel kerneltype(foo)
_s2psm_cmp checks: neighbor(0): t, pscore(ps) outcome(y) neighbor(0)

* Contracts: value labels, _weight's characteristic, _nn on support
capture noisily {
    capture drop _*
    quietly cpsmatch t x1 x2, outcome(y) ate
    assert "`: value label _treated'" == "_treated" & "`: value label _support'" == "_support"
    assert "`: label _treated 0'" == "Untreated" & "`: label _treated 1'" == "Treated"
    assert "`: label _support 0'" == "Off support" & "`: label _support 1'" == "On support"
    assert "`: char _weight[Type]'" == "aweight"
    assert _nn < . if _support == 1
    assert !missing(r(atu), r(ate), r(atu_y), r(ate_y))
    assert missing(r(seatu), r(seate))
    quietly cpsmatch t x1 x2, outcome(y) ai(2) ate
    assert !missing(r(seatu), r(seate), r(seatu_y), r(seate_y))
    assert r(method) == "nearest" & r(metric) == "pscore"
}
if _rc test_fail "PSM-14 labels, characteristic and r()" "r(`=_rc')"
else test_pass "PSM-14 labels, characteristic and r()"

* Several outcomes: each outcome's matched values and effects equal a
* run with that outcome alone, over the common sample
capture noisily {
    capture drop _*
    quietly cpsmatch t x1 x2 if !missing(y2), outcome(y) ai(2) ate
    local a = r(ate)
    local s = r(seate)
    quietly gen double one_y = _y
    capture drop _*
    quietly cpsmatch t x1 x2, outcome(y y2) ai(2) ate
    assert reldif(r(ate_y), `a') < 1e-12 & reldif(r(seate_y), `s') < 1e-12
    assert _y == one_y | (missing(_y) & missing(one_y))
    drop one_y
}
if _rc test_fail "PSM-01 several outcomes with ate and ai()" "r(`=_rc')"
else test_pass "PSM-01 several outcomes with ate and ai()"

capture drop _*
print_summary s2_psm
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_psm"
    global CTOOLS_CASE_COMPLETE "s2_psm"
    exit
}

di as error "Unknown s2_psm case: `__case'"
exit 198
