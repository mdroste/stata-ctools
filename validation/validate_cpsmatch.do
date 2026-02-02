/*******************************************************************************
 * validate_cpsmatch.do
 *
 * Validation tests for cpsmatch propensity score matching
 * Compares cpsmatch results with psmatch2 for all matching methods
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for cpsmatch..."

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
print_section "Plugin Check"

* Check psmatch2 is installed
capture which psmatch2
if _rc != 0 {
    noi di as error "psmatch2 not installed - skipping comparison tests"
    noi di as error "Install with: ssc install psmatch2"
    test_fail "psmatch2 availability" "not installed"
    noi print_summary "cpsmatch"
    exit 1
}
test_pass "psmatch2 available"

* Plugin test
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.5 + 0.3*x1 + 0.2*x2 + rnormal() > 0.5)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

capture noisily cpsmatch treat x1 x2
if _rc != 0 {
    local errc = _rc
    test_fail "cpsmatch plugin load" "plugin returned error `errc'"
    noi print_summary "cpsmatch"
    exit 1
}
test_pass "cpsmatch plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic functionality and return values
 ******************************************************************************/
print_section "Basic Functionality & Return Values"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.5 + 0.3*x1 + 0.2*x2 + rnormal() > 0.5)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

cpsmatch treat x1 x2, outcome(y)

* Store ALL r() values before they get overwritten
local cps_N = r(N)
local cps_n_treated = r(n_treated)
local cps_n_controls = r(n_controls)
local cps_n_matched = r(n_matched)
local cps_n_off_support = r(n_off_support)
local cps_common_min = r(common_min)
local cps_common_max = r(common_max)
local cps_att = r(att)
local cps_att_se = r(att_se)
local cps_att_t = r(att_t)
local cps_method = "`r(method)'"
local cps_treatvar = "`r(treatvar)'"
local cps_outcome = "`r(outcome)'"

* Test variable creation
capture confirm numeric variable _pscore
if _rc == 0 {
    test_pass "creates _pscore variable"
}
else {
    test_fail "creates _pscore" "variable not created"
}

capture confirm numeric variable _weight
if _rc == 0 {
    test_pass "creates _weight variable"
}
else {
    test_fail "creates _weight" "variable not created"
}

capture confirm numeric variable _support
if _rc == 0 {
    test_pass "creates _support variable"
}
else {
    test_fail "creates _support" "variable not created"
}

sum _pscore, meanonly
local pmin = r(min)
local pmax = r(max)
if `pmin' >= 0 & `pmax' <= 1 {
    test_pass "pscores in [0,1] range"
}
else {
    test_fail "pscores range" "min=`pmin', max=`pmax'"
}

* Test ALL scalar return values
local all_scalars_ok = 1

* r(N) - should equal number of observations
if `cps_N' != 500 {
    local all_scalars_ok = 0
}

* r(n_treated) + r(n_controls) should equal N (approximately, excluding missings)
if `cps_n_treated' <= 0 | `cps_n_controls' <= 0 {
    local all_scalars_ok = 0
}

* r(n_matched) should be positive
if `cps_n_matched' <= 0 {
    local all_scalars_ok = 0
}

* r(n_off_support) should be non-negative
if `cps_n_off_support' < 0 | `cps_n_off_support' == . {
    local all_scalars_ok = 0
}

* r(common_min) and r(common_max) should be valid pscores
if `cps_common_min' < 0 | `cps_common_min' > 1 | `cps_common_max' < 0 | `cps_common_max' > 1 {
    local all_scalars_ok = 0
}

* r(att) should be defined (we specified outcome)
if `cps_att' == . {
    local all_scalars_ok = 0
}

* r(att_se) should be positive
if `cps_att_se' <= 0 | `cps_att_se' == . {
    local all_scalars_ok = 0
}

* r(att_t) should be defined
if `cps_att_t' == . {
    local all_scalars_ok = 0
}

if `all_scalars_ok' {
    test_pass "returns all scalar values (N, n_treated, n_controls, n_matched, n_off_support, common_min, common_max, att, att_se, att_t)"
}
else {
    test_fail "scalar returns" "N=`cps_N' n_treated=`cps_n_treated' n_controls=`cps_n_controls' n_matched=`cps_n_matched' n_off=`cps_n_off_support' att=`cps_att'"
}

* Test local return values
local all_locals_ok = 1

* r(method) should be "nearest" (default)
if "`cps_method'" != "nearest" {
    local all_locals_ok = 0
}

* r(treatvar) should be "treat"
if "`cps_treatvar'" != "treat" {
    local all_locals_ok = 0
}

* r(outcome) should be "y"
if "`cps_outcome'" != "y" {
    local all_locals_ok = 0
}

if `all_locals_ok' {
    test_pass "returns all local values (method, treatvar, outcome)"
}
else {
    test_fail "local returns" "method=`cps_method' treatvar=`cps_treatvar' outcome=`cps_outcome'"
}

/*******************************************************************************
 * SECTION 3: psmatch2 comparison - Nearest Neighbor
 ******************************************************************************/
print_section "Nearest Neighbor (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Run psmatch2
psmatch2 treat x1 x2, outcome(y)
local psm2_att = r(att)
tempvar pscore_psm2
gen double `pscore_psm2' = _pscore

* Clean up psmatch2 variables
capture capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

* Run cpsmatch
cpsmatch treat x1 x2, outcome(y)
local cps_att = r(att)

* Compare ATT
sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "NN: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "NN: ATT matches psmatch2" "sigfigs=`sf_fmt'"
}

* Compare propensity scores
gen double _pdiff = abs(_pscore - `pscore_psm2')
sum _pdiff, meanonly
local maxdiff = r(max)
if `maxdiff' < 1e-6 {
    test_pass "NN: pscores match psmatch2"
}
else {
    test_fail "NN: pscores match" "max diff=`maxdiff'"
}
drop _pdiff

/*******************************************************************************
 * SECTION 4: psmatch2 comparison - Multiple Neighbors
 ******************************************************************************/
print_section "Multiple Neighbors (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y) neighbor(3)
local psm2_att = r(att)
capture capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) neighbor(3)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "NN(3): ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "NN(3): ATT matches" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 5: psmatch2 comparison - Kernel Matching
 ******************************************************************************/
print_section "Kernel Matching (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y) kernel
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id

cpsmatch treat x1 x2, outcome(y) kernel
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "Kernel: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "Kernel: ATT matches" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 6: psmatch2 comparison - Radius Matching
 ******************************************************************************/
print_section "Radius Matching (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y) radius caliper(0.1)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id

cpsmatch treat x1 x2, outcome(y) radius caliper(0.1)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "Radius: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "Radius: ATT matches" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 7: psmatch2 comparison - Common Support
 ******************************************************************************/
print_section "Common Support (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y) common
local psm2_att = r(att)
capture capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) common
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "Common support: ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "Common support: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 8: psmatch2 comparison - Without Replacement
 ******************************************************************************/
print_section "Without Replacement (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y) noreplacement
local psm2_att = r(att)
capture capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) noreplacement
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "Noreplacement: ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "Noreplacement: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 9: psmatch2 comparison - Logit
 ******************************************************************************/
print_section "Logit Estimation (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y) logit
local psm2_att = r(att)
capture capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) logit
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "Logit: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "Logit: ATT matches" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 10: Large datasets
 ******************************************************************************/
print_section "Large Datasets"

clear
set obs 10000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

psmatch2 treat x1 x2, outcome(y)
local psm2_att = r(att)
capture capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "10K obs: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "10K obs: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 11: cattaneo2 Example (from help file)
 ******************************************************************************/
print_section "cattaneo2 Example"

* Basic example from help file
webuse cattaneo2, clear
capture drop _pscore _weight _id _support _treated _nn

psmatch2 mbsmoke mage medu, outcome(bweight)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch mbsmoke mage medu, outcome(bweight)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "cattaneo2: basic example ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "cattaneo2: basic example" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 12: if/in Conditions
 ******************************************************************************/
print_section "if/in Conditions"

* Test with if condition
clear
set obs 1000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()
gen group = mod(_n, 2)

* Test if condition
preserve
keep if group == 1
psmatch2 treat x1 x2, outcome(y)
local psm2_att_if = r(att)
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2 if group == 1, outcome(y)
local cps_att_if = r(att)
capture drop _pscore _weight _treated _support _nn _id

sigfigs `psm2_att_if' `cps_att_if'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "if condition: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "if condition: ATT" "sigfigs=`sf_fmt'"
}

* Test in condition
preserve
keep in 1/400
psmatch2 treat x1 x2, outcome(y)
local psm2_att_in = r(att)
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2 in 1/400, outcome(y)
local cps_att_in = r(att)
capture drop _pscore _weight _treated _support _nn _id

sigfigs `psm2_att_in' `cps_att_in'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "in condition: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "in condition: ATT" "sigfigs=`sf_fmt'"
}

* Test if and in combined
preserve
keep if group == 1 in 1/600
psmatch2 treat x1 x2, outcome(y)
local psm2_att_ifin = r(att)
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2 if group == 1 in 1/600, outcome(y)
local cps_att_ifin = r(att)
capture drop _pscore _weight _treated _support _nn _id

sigfigs `psm2_att_ifin' `cps_att_ifin'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "if/in combined: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "if/in combined: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 13: Pre-estimated Propensity Score
 ******************************************************************************/
print_section "Pre-estimated Propensity Score"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Estimate propensity score externally
probit treat x1 x2
predict double pscore_ext, pr

psmatch2 treat, pscore(pscore_ext) outcome(y)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat, pscore(pscore_ext) outcome(y)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "pscore(): ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "pscore(): ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 14: Kernel Types
 ******************************************************************************/
print_section "Kernel Types"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Normal/Gaussian kernel
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel kerneltype(normal)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel kerneltype(normal)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "kernel(normal): ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "kernel(normal): ATT" "sigfigs=`sf_fmt'"
}

* Biweight kernel
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel kerneltype(biweight)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel kerneltype(biweight)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "kernel(biweight): ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "kernel(biweight): ATT" "sigfigs=`sf_fmt'"
}

* Uniform kernel
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel kerneltype(uniform)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel kerneltype(uniform)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "kernel(uniform): ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "kernel(uniform): ATT" "sigfigs=`sf_fmt'"
}

* Tricube kernel
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel kerneltype(tricube)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel kerneltype(tricube)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "kernel(tricube): ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "kernel(tricube): ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 15: Bandwidth
 ******************************************************************************/
print_section "Bandwidth"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Custom bandwidth
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel bwidth(0.08)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel bwidth(0.08)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "bwidth(0.08): ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "bwidth(0.08): ATT" "sigfigs=`sf_fmt'"
}

* Smaller bandwidth
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel bwidth(0.03)
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel bwidth(0.03)
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "bwidth(0.03): ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "bwidth(0.03): ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 16: Ties
 ******************************************************************************/
print_section "Ties"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Test ties option
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) ties
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) ties
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "ties: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "ties: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 17: Descending
 ******************************************************************************/
print_section "Descending"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Test descending option (affects matching order for noreplacement)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) noreplacement descending
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) noreplacement descending
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "descending: ATT matches psmatch2"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "descending: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 18: Combined Options
 ******************************************************************************/
print_section "Combined Options"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* NN with multiple options: neighbor(3), common, logit
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) neighbor(3) common logit
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) neighbor(3) common logit
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "neighbor(3) common logit: ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "neighbor(3) common logit: ATT" "sigfigs=`sf_fmt'"
}

* Kernel with common support
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) kernel common
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) kernel common
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "kernel common: ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "kernel common: ATT" "sigfigs=`sf_fmt'"
}

* Radius with caliper and common support
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif
psmatch2 treat x1 x2, outcome(y) radius caliper(0.05) common
local psm2_att = r(att)
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y) radius caliper(0.05) common
local cps_att = r(att)

sigfigs `psm2_att' `cps_att'
local sf = r(sigfigs)
if `sf' >= 7 {
    test_pass "radius caliper(0.05) common: ATT matches"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "radius caliper(0.05) common: ATT" "sigfigs=`sf_fmt'"
}

/*******************************************************************************
 * SECTION 19: cpsmatch-only Features
 ******************************************************************************/
print_section "cpsmatch-only Features"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen treat = (0.3*x1 + rnormal() > 0)

capture cpsmatch treat x1, verbose
if _rc == 0 {
    test_pass "verbose option"
}
else {
    test_fail "verbose option" "rc=`=_rc'"
}

capture drop _pscore _weight _treated _support _nn _id

capture cpsmatch treat x1, threads(2)
if _rc == 0 {
    test_pass "threads(2) option"
}
else {
    test_fail "threads option" "rc=`=_rc'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cpsmatch"
}
