/*******************************************************************************
 * validate_civreghdfe_gmm_cue.do
 *
 * Validation tests for civreghdfe GMM2S and CUE options
 * All tests use 1e-8 tolerance for PASS criteria
 ******************************************************************************/

* Load setup
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

local tol = 1e-8

di as text ""
di as text "======================================================================"
di as text "         GMM2S AND CUE VALIDATION TESTS (tol=1e-8)"
di as text "======================================================================"

/*******************************************************************************
 * Check prerequisites
 ******************************************************************************/
noi print_section "Prerequisite Check"

capture which ivreghdfe
if _rc != 0 {
    noi di as error "ivreghdfe not installed. Please run: ssc install ivreghdfe"
    exit 1
}
noi test_pass "ivreghdfe is installed"

/*******************************************************************************
 * SECTION 1: GMM2S Coefficient and VCE Tests
 ******************************************************************************/
noi print_section "GMM2S Coefficient/VCE Tests"

clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* GMM2S basic
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s
local iv_b = _b[x_endog]
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s
local c_b = _b[x_endog]
local c_V11 = e(V)[1,1]

local b_diff = abs(`iv_b' - `c_b')
local V_diff = abs(`c_V11' - `iv_V11')

if `b_diff' < `tol' & `V_diff' < `tol' {
    noi test_pass "gmm2s basic (coef diff=`b_diff', VCE diff=`V_diff')"
}
else {
    noi test_fail "gmm2s basic" "coef diff=`b_diff', VCE diff=`V_diff'"
}

* GMM2S + robust
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust)
local iv_b = _b[x_endog]
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust)
local c_b = _b[x_endog]
local c_V11 = e(V)[1,1]

local b_diff = abs(`iv_b' - `c_b')
local V_diff = abs(`c_V11' - `iv_V11')

if `b_diff' < `tol' & `V_diff' < `tol' {
    noi test_pass "gmm2s + robust (coef diff=`b_diff', VCE diff=`V_diff')"
}
else {
    noi test_fail "gmm2s + robust" "coef diff=`b_diff', VCE diff=`V_diff'"
}

* GMM2S + cluster
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local iv_b = _b[x_endog]
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_b = _b[x_endog]
local c_V11 = e(V)[1,1]

local b_diff = abs(`iv_b' - `c_b')
local V_diff = abs(`c_V11' - `iv_V11')

if `b_diff' < `tol' & `V_diff' < `tol' {
    noi test_pass "gmm2s + cluster (coef diff=`b_diff', VCE diff=`V_diff')"
}
else {
    noi test_fail "gmm2s + cluster" "coef diff=`b_diff', VCE diff=`V_diff'"
}

/*******************************************************************************
 * SECTION 2: CUE Coefficient and VCE Tests
 ******************************************************************************/
noi print_section "CUE Coefficient/VCE Tests"

* CUE basic
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue
local iv_b = _b[x_endog]
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue
local c_b = _b[x_endog]
local c_V11 = e(V)[1,1]

local b_diff = abs(`iv_b' - `c_b')
local V_diff = abs(`c_V11' - `iv_V11')

if `b_diff' < `tol' & `V_diff' < `tol' {
    noi test_pass "cue basic (coef diff=`b_diff', VCE diff=`V_diff')"
}
else {
    noi test_fail "cue basic" "coef diff=`b_diff', VCE diff=`V_diff'"
}

* CUE + robust
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(robust)
local iv_b = _b[x_endog]
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(robust)
local c_b = _b[x_endog]
local c_V11 = e(V)[1,1]

local b_diff = abs(`iv_b' - `c_b')
local V_diff = abs(`c_V11' - `iv_V11')

if `b_diff' < `tol' & `V_diff' < `tol' {
    noi test_pass "cue + robust (coef diff=`b_diff', VCE diff=`V_diff')"
}
else {
    noi test_fail "cue + robust" "coef diff=`b_diff', VCE diff=`V_diff'"
}

* CUE + cluster
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(cluster firm)
local iv_b = _b[x_endog]
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(cluster firm)
local c_b = _b[x_endog]
local c_V11 = e(V)[1,1]

local b_diff = abs(`iv_b' - `c_b')
local V_diff = abs(`c_V11' - `iv_V11')

if `b_diff' < `tol' & `V_diff' < `tol' {
    noi test_pass "cue + cluster (coef diff=`b_diff', VCE diff=`V_diff')"
}
else {
    noi test_fail "cue + cluster" "coef diff=`b_diff', VCE diff=`V_diff'"
}

/*******************************************************************************
 * SECTION 3: Just-Identified Equivalence
 ******************************************************************************/
noi print_section "Just-Identified Case"

clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen x_endog = 0.5*z1 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* 2SLS
civreghdfe y (x_endog = z1) x_exog, absorb(firm)
local b_2sls = _b[x_endog]

* GMM2S should equal 2SLS when just-identified
civreghdfe y (x_endog = z1) x_exog, absorb(firm) gmm2s
local b_gmm = _b[x_endog]

* CUE should equal 2SLS when just-identified
civreghdfe y (x_endog = z1) x_exog, absorb(firm) cue
local b_cue = _b[x_endog]

local diff_gmm = abs(`b_2sls' - `b_gmm')
local diff_cue = abs(`b_2sls' - `b_cue')

if `diff_gmm' < `tol' {
    noi test_pass "gmm2s = 2sls when just-identified (diff=`diff_gmm')"
}
else {
    noi test_fail "gmm2s just-identified" "diff=`diff_gmm'"
}

if `diff_cue' < `tol' {
    noi test_pass "cue = 2sls when just-identified (diff=`diff_cue')"
}
else {
    noi test_fail "cue just-identified" "diff=`diff_cue'"
}

/*******************************************************************************
 * SECTION 4: Multiple Endogenous Variables
 ******************************************************************************/
noi print_section "Multiple Endogenous"

clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen z4 = rnormal() + 0.2*z1
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z3 + 0.3*z4 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog1 + 1.5*x_endog2 + 1.0*x_exog + rnormal()

* GMM2S two endog - check coefficients
qui ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) gmm2s
local iv_b1 = _b[x_endog1]
local iv_b2 = _b[x_endog2]

qui civreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) gmm2s
local c_b1 = _b[x_endog1]
local c_b2 = _b[x_endog2]

local diff1 = abs(`iv_b1' - `c_b1')
local diff2 = abs(`iv_b2' - `c_b2')

if `diff1' < `tol' & `diff2' < `tol' {
    noi test_pass "gmm2s two endog (diff1=`diff1', diff2=`diff2')"
}
else {
    noi test_fail "gmm2s two endog" "diff1=`diff1', diff2=`diff2'"
}

* CUE two endog - check coefficients
qui ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) cue
local iv_b1 = _b[x_endog1]
local iv_b2 = _b[x_endog2]

qui civreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) cue
local c_b1 = _b[x_endog1]
local c_b2 = _b[x_endog2]

local diff1 = abs(`iv_b1' - `c_b1')
local diff2 = abs(`iv_b2' - `c_b2')

if `diff1' < `tol' & `diff2' < `tol' {
    noi test_pass "cue two endog (diff1=`diff1', diff2=`diff2')"
}
else {
    noi test_fail "cue two endog" "diff1=`diff1', diff2=`diff2'"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "civreghdfe GMM2S and CUE"

if $TESTS_FAILED > 0 {
    exit 1
}

}
