* Test creghdfe weight support
* Compare results with reghdfe

clear all
set more off

* Load ctools from build directory
adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

* Load test dataset
sysuse auto, clear

* Create weight variable (handle missing rep78)
gen w = rep78
replace w = 3 if missing(w)

* Display header
di as text ""
di as text "{hline 70}"
di as text "Testing creghdfe weight support"
di as text "{hline 70}"
di as text ""

* ============================================================================
* Test 1: Unweighted baseline (should match exactly)
* ============================================================================
di as text "{hline 70}"
di as text "Test 1: Unweighted baseline"
di as text "{hline 70}"

quietly reghdfe price mpg weight, absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_b_weight = _b[weight]
local reghdfe_se_mpg = _se[mpg]
di as text "reghdfe: mpg = " %9.4f `reghdfe_b_mpg' " (se = " %9.4f `reghdfe_se_mpg' ")"

quietly creghdfe price mpg weight, absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_b_weight = _b[weight]
local creghdfe_se_mpg = _se[mpg]
di as text "creghdfe: mpg = " %9.4f `creghdfe_b_mpg' " (se = " %9.4f `creghdfe_se_mpg' ")"

local diff = abs(`reghdfe_b_mpg' - `creghdfe_b_mpg')
if `diff' < 0.0001 {
    di as result "PASS: Unweighted coefficients match"
}
else {
    di as error "FAIL: Coefficient difference = " `diff'
}

* ============================================================================
* Test 2: Analytic weights (aweight)
* ============================================================================
di as text ""
di as text "{hline 70}"
di as text "Test 2: Analytic weights [aweight=w]"
di as text "{hline 70}"

quietly reghdfe price mpg weight [aw=w], absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
di as text "reghdfe: mpg = " %9.4f `reghdfe_b_mpg' " (se = " %9.4f `reghdfe_se_mpg' ")"

quietly creghdfe price mpg weight [aw=w], absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
di as text "creghdfe: mpg = " %9.4f `creghdfe_b_mpg' " (se = " %9.4f `creghdfe_se_mpg' ")"

local diff = abs(`reghdfe_b_mpg' - `creghdfe_b_mpg')
if `diff' < 0.0001 {
    di as result "PASS: aweight coefficients match"
}
else {
    di as error "FAIL: Coefficient difference = " `diff'
}

local se_diff = abs(`reghdfe_se_mpg' - `creghdfe_se_mpg')
if `se_diff' < 0.0001 {
    di as result "PASS: aweight standard errors match"
}
else {
    di as error "FAIL: SE difference = " `se_diff'
}

* ============================================================================
* Test 3: Frequency weights (fweight)
* ============================================================================
di as text ""
di as text "{hline 70}"
di as text "Test 3: Frequency weights [fweight=rep78]"
di as text "{hline 70}"

* Note: fweight requires integer values and no missing
preserve
drop if missing(rep78)

quietly reghdfe price mpg weight [fw=rep78], absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
di as text "reghdfe: mpg = " %9.4f `reghdfe_b_mpg' " (se = " %9.4f `reghdfe_se_mpg' ")"

quietly creghdfe price mpg weight [fw=rep78], absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
di as text "creghdfe: mpg = " %9.4f `creghdfe_b_mpg' " (se = " %9.4f `creghdfe_se_mpg' ")"

local diff = abs(`reghdfe_b_mpg' - `creghdfe_b_mpg')
if `diff' < 0.0001 {
    di as result "PASS: fweight coefficients match"
}
else {
    di as error "FAIL: Coefficient difference = " `diff'
}

restore

* ============================================================================
* Test 4: Probability weights (pweight) - forces robust
* ============================================================================
di as text ""
di as text "{hline 70}"
di as text "Test 4: Probability weights [pweight=w]"
di as text "{hline 70}"

quietly reghdfe price mpg weight [pw=w], absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
di as text "reghdfe: mpg = " %9.4f `reghdfe_b_mpg' " (se = " %9.4f `reghdfe_se_mpg' ")"

quietly creghdfe price mpg weight [pw=w], absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
di as text "creghdfe: mpg = " %9.4f `creghdfe_b_mpg' " (se = " %9.4f `creghdfe_se_mpg' ")"

local diff = abs(`reghdfe_b_mpg' - `creghdfe_b_mpg')
if `diff' < 0.0001 {
    di as result "PASS: pweight coefficients match"
}
else {
    di as error "FAIL: Coefficient difference = " `diff'
}

* ============================================================================
* Test 5: Weights with robust VCE
* ============================================================================
di as text ""
di as text "{hline 70}"
di as text "Test 5: aweight with robust VCE"
di as text "{hline 70}"

quietly reghdfe price mpg weight [aw=w], absorb(foreign) vce(robust)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
di as text "reghdfe: mpg = " %9.4f `reghdfe_b_mpg' " (se = " %9.4f `reghdfe_se_mpg' ")"

quietly creghdfe price mpg weight [aw=w], absorb(foreign) vce(robust)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
di as text "creghdfe: mpg = " %9.4f `creghdfe_b_mpg' " (se = " %9.4f `creghdfe_se_mpg' ")"

local diff = abs(`reghdfe_b_mpg' - `creghdfe_b_mpg')
if `diff' < 0.0001 {
    di as result "PASS: aweight+robust coefficients match"
}
else {
    di as error "FAIL: Coefficient difference = " `diff'
}

* ============================================================================
* Test 6: Weights with cluster VCE
* ============================================================================
di as text ""
di as text "{hline 70}"
di as text "Test 6: aweight with cluster VCE"
di as text "{hline 70}"

quietly reghdfe price mpg weight [aw=w], absorb(foreign) vce(cluster foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
di as text "reghdfe: mpg = " %9.4f `reghdfe_b_mpg' " (se = " %9.4f `reghdfe_se_mpg' ")"

quietly creghdfe price mpg weight [aw=w], absorb(foreign) vce(cluster foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
di as text "creghdfe: mpg = " %9.4f `creghdfe_b_mpg' " (se = " %9.4f `creghdfe_se_mpg' ")"

local diff = abs(`reghdfe_b_mpg' - `creghdfe_b_mpg')
if `diff' < 0.0001 {
    di as result "PASS: aweight+cluster coefficients match"
}
else {
    di as error "FAIL: Coefficient difference = " `diff'
}

* ============================================================================
* Test 7: weights=1 should match unweighted
* ============================================================================
di as text ""
di as text "{hline 70}"
di as text "Test 7: All weights=1 should match unweighted"
di as text "{hline 70}"

gen w1 = 1

quietly creghdfe price mpg weight, absorb(foreign)
local unweighted_b = _b[mpg]

quietly creghdfe price mpg weight [aw=w1], absorb(foreign)
local weighted1_b = _b[mpg]

local diff = abs(`unweighted_b' - `weighted1_b')
if `diff' < 0.0001 {
    di as result "PASS: weights=1 matches unweighted (diff = " %9.6f `diff' ")"
}
else {
    di as error "FAIL: weights=1 differs from unweighted by " `diff'
}

di as text ""
di as text "{hline 70}"
di as text "Weight tests complete"
di as text "{hline 70}"
