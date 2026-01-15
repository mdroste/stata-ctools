/*******************************************************************************
 * validate_civreghdfe.do
 *
 * Comprehensive validation tests for civreghdfe vs ivreghdfe
 * Tests all options: absorb, vce, weights, first, small, tolerance, maxiter
 *
 * Requires ivreghdfe to be installed (ssc install ivreghdfe)
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CIVREGHDFE VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * Check if ivreghdfe is installed
 ******************************************************************************/
noi print_section "Prerequisite Check"

capture which ivreghdfe
if _rc != 0 {
    noi di as error "ivreghdfe not installed. Please run: ssc install ivreghdfe"
    noi di as error "Skipping civreghdfe validation tests"
    noi test_fail "ivreghdfe installation" "ivreghdfe not found"
    noi print_summary "civreghdfe"
    exit 0
}
noi test_pass "ivreghdfe is installed"

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)
if _rc != 0 {
    noi test_fail "civreghdfe plugin load" "returned error `=_rc'"
    noi print_summary "civreghdfe"
    exit 1
}
noi test_pass "civreghdfe plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic IV tests (auto)
 ******************************************************************************/
noi print_section "Basic IV Tests (auto)"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("single endog, single instr")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("single endog, two instr")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) testname("with exog var")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight length) turn displacement, absorb(foreign) testname("with multiple exog")

/*******************************************************************************
 * SECTION 3: Multiple endogenous variables
 ******************************************************************************/
noi print_section "Multiple Endogenous Variables"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg weight = length turn displacement), absorb(foreign) testname("two endog")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg weight = length turn displacement headroom), absorb(foreign) testname("two endog, more instr")

/*******************************************************************************
 * SECTION 4: VCE options
 ******************************************************************************/
noi print_section "VCE Options"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("vce(robust)")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("vce(cluster)")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("multi-instr + robust")

/*******************************************************************************
 * SECTION 5: Two-way fixed effects
 ******************************************************************************/
noi print_section "Two-Way Fixed Effects"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) vce(robust) testname("two-way + robust")

/*******************************************************************************
 * SECTION 6: Census dataset
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
noi benchmark_ivreghdfe pop (medage = death), absorb(region) testname("census basic")

sysuse census, clear
noi benchmark_ivreghdfe pop (medage = death marriage), absorb(region) testname("census two instr")

sysuse census, clear
noi benchmark_ivreghdfe pop (medage = death) divorce, absorb(region) testname("census with exog")

/*******************************************************************************
 * SECTION 7: Panel data (nlswork)
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) testname("nlswork basic")

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age ttl_exp), absorb(idcode) testname("nlswork two instr")

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(robust) testname("nlswork robust")

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(cluster idcode) testname("nlswork cluster")

/*******************************************************************************
 * SECTION 8: Large dataset
 ******************************************************************************/
noi print_section "Large Dataset"

clear
set seed 12345
set obs 20000
gen id = runiformint(1, 200)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id)
matrix ivreghdfe_b = e(b)
local ivreghdfe_N = e(N)

civreghdfe y (x_endog = z1 z2) x_exog, absorb(id)
matrix civreghdfe_b = e(b)
local civreghdfe_N = e(N)

if `ivreghdfe_N' == `civreghdfe_N' {
    noi test_pass "20K dataset: N matches"
}
else {
    noi test_fail "20K dataset" "N differs"
}

/*******************************************************************************
 * SECTION 9: first option
 ******************************************************************************/
noi print_section "first Option"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) first
if _rc == 0 {
    noi test_pass "first option accepted"
}
else {
    noi test_fail "first option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 10: small option
 ******************************************************************************/
noi print_section "small Option"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) small
if _rc == 0 {
    noi test_pass "small option accepted"
}
else {
    noi test_fail "small option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: tolerance/maxiter options
 ******************************************************************************/
noi print_section "tolerance/maxiter Options"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10)
if _rc == 0 {
    noi test_pass "tolerance(1e-10) accepted"
}
else {
    noi test_fail "tolerance option" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) maxiter(1000)
if _rc == 0 {
    noi test_pass "maxiter(1000) accepted"
}
else {
    noi test_fail "maxiter option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 12: verbose/timeit options
 ******************************************************************************/
noi print_section "verbose/timeit Options"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) timeit
if _rc == 0 {
    noi test_pass "timeit option accepted"
}
else {
    noi test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 13: if/in conditions
 ******************************************************************************/
noi print_section "if/in Conditions"

sysuse auto, clear
ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign)
local ivreghdfe_N = e(N)

civreghdfe price (mpg = weight) if price > 5000, absorb(foreign)
local civreghdfe_N = e(N)

if `ivreghdfe_N' == `civreghdfe_N' {
    noi test_pass "if condition: N matches"
}
else {
    noi test_fail "if condition" "N differs"
}

sysuse auto, clear
ivreghdfe price (mpg = weight) in 1/50, absorb(foreign)
local ivreghdfe_N = e(N)

civreghdfe price (mpg = weight) in 1/50, absorb(foreign)
local civreghdfe_N = e(N)

if `ivreghdfe_N' == `civreghdfe_N' {
    noi test_pass "in condition: N matches"
}
else {
    noi test_fail "in condition" "N differs"
}

/*******************************************************************************
 * SECTION 14: Coefficient comparison
 ******************************************************************************/
noi print_section "Coefficient Comparison"

sysuse auto, clear

ivreghdfe price (mpg = weight length), absorb(foreign)
matrix ivreghdfe_b = e(b)

civreghdfe price (mpg = weight length), absorb(foreign)
matrix civreghdfe_b = e(b)

tempname diff
matrix `diff' = ivreghdfe_b - civreghdfe_b
local cols = colsof(`diff')
local maxdiff = 0
forvalues j = 1/`cols' {
    local d = abs(`diff'[1, `j'])
    if `d' > `maxdiff' local maxdiff = `d'
}

if `maxdiff' < 1e-6 {
    noi test_pass "coefficients match (maxdiff=`maxdiff')"
}
else {
    noi test_fail "coefficients" "maxdiff=`maxdiff' exceeds tolerance"
}

/*******************************************************************************
 * SECTION 15: Weights
 ******************************************************************************/
noi print_section "Weights"

sysuse auto, clear
capture civreghdfe price (mpg = weight) [aw=weight], absorb(foreign)
if _rc == 0 {
    noi test_pass "aweight accepted"
}
else {
    noi test_fail "aweight" "returned error `=_rc'"
}

sysuse auto, clear
gen int fw = ceil(mpg/5)
capture civreghdfe price (mpg = weight) [fw=fw], absorb(foreign)
if _rc == 0 {
    noi test_pass "fweight accepted"
}
else {
    noi test_fail "fweight" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight) [pw=weight], absorb(foreign)
if _rc == 0 {
    noi test_pass "pweight accepted"
}
else {
    noi test_fail "pweight" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 16: Display options
 ******************************************************************************/
noi print_section "Display Options"

* Test level()
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
if e(level) == 90 {
    noi test_pass "level(90) sets correct confidence level"
}
else {
    noi test_fail "level(90)" "e(level) = `e(level)' instead of 90"
}

* Test noheader - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader
if _rc == 0 {
    noi test_pass "noheader option accepted"
}
else {
    noi test_fail "noheader" "returned error `=_rc'"
}

* Test nofooter - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nofooter
if _rc == 0 {
    noi test_pass "nofooter option accepted"
}
else {
    noi test_fail "nofooter" "returned error `=_rc'"
}

* Test nooutput - verify it runs and stores results
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nooutput
if _rc == 0 & e(N) > 0 {
    noi test_pass "nooutput option accepted and stores e(N)"
}
else {
    noi test_fail "nooutput" "returned error or missing e(N)"
}

* Test title() - verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) title("Custom Title Test")
if _rc == 0 {
    noi test_pass "title() option accepted"
}
else {
    noi test_fail "title()" "returned error `=_rc'"
}

* Test depname() - verify it sets e(depvar)
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) depname("outcome_var")
if "`e(depvar)'" == "outcome_var" {
    noi test_pass "depname() sets e(depvar) correctly"
}
else {
    noi test_fail "depname()" "e(depvar) = `e(depvar)' instead of outcome_var"
}

* Test noid - verify underid test is suppressed but e(idstat) still stored
sysuse auto, clear
civreghdfe price (mpg = weight length), absorb(foreign) noid
if e(idstat) != . {
    noi test_pass "noid suppresses display but e(idstat) still computed"
}
else {
    noi test_fail "noid" "e(idstat) not stored"
}

* Test combined options
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)
if _rc == 0 & e(level) == 99 {
    noi test_pass "combined display options work together"
}
else {
    noi test_fail "combined options" "returned error or wrong level"
}

/*******************************************************************************
 * SECTION 17: Two-way clustering
 ******************************************************************************/
noi print_section "Two-Way Clustering"

* Create test data with firm-time structure
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* Test two-way clustering syntax
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
if _rc == 0 {
    noi test_pass "two-way clustering syntax accepted"
}
else {
    noi test_fail "two-way clustering" "returned error `=_rc'"
}

* Test that both cluster counts are stored
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
if e(N_clust1) > 0 & e(N_clust2) > 0 {
    noi test_pass "two-way clustering stores N_clust1 and N_clust2"
}
else {
    noi test_fail "two-way N_clust" "N_clust1 or N_clust2 not stored"
}

* Test that vcetype is set correctly
if "`e(vcetype)'" == "Two-way Cluster" {
    noi test_pass "two-way clustering sets vcetype correctly"
}
else {
    noi test_fail "two-way vcetype" "vcetype = `e(vcetype)' instead of Two-way Cluster"
}

* Test that two-way SE >= max(one-way SEs)
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
local se_twoway = sqrt(e(V)[1,1])

civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm)
local se_firm = sqrt(e(V)[1,1])

civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster time)
local se_time = sqrt(e(V)[1,1])

* Two-way SE should be >= max of one-way SEs (with small tolerance)
if `se_twoway' >= max(`se_firm', `se_time') - 0.001 {
    noi test_pass "two-way SE >= max(one-way SEs)"
}
else {
    noi test_fail "two-way SE" "se_twoway=`se_twoway' < max(se_firm=`se_firm', se_time=`se_time')"
}

* Test that clustvar1 and clustvar2 are stored
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
if "`e(clustvar1)'" == "firm" & "`e(clustvar2)'" == "time" {
    noi test_pass "two-way clustering stores clustvar1 and clustvar2"
}
else {
    noi test_fail "two-way clustvars" "clustvar1=`e(clustvar1)', clustvar2=`e(clustvar2)'"
}

/*******************************************************************************
 * SECTION 18: Diagnostic test options (orthog, endogtest, redundant)
 ******************************************************************************/
noi print_section "Diagnostic Test Options"

* Create test data with multiple instruments
scalar drop _all
matrix drop _all
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

* Test orthog() option
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3)
if _rc == 0 {
    noi test_pass "orthog() option accepted"
}
else {
    noi test_fail "orthog()" "returned error `=_rc'"
}

* Test orthog() stores results
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3)
if e(cstat) != . & e(cstat_df) != . & e(cstat_p) != . {
    noi test_pass "orthog() stores cstat, cstat_df, cstat_p"
}
else {
    noi test_fail "orthog() results" "cstat=`e(cstat)', df=`e(cstat_df)', p=`e(cstat_p)'"
}

* Test C-stat is non-negative and p-value is in [0,1]
if e(cstat) >= 0 & e(cstat_p) >= 0 & e(cstat_p) <= 1 {
    noi test_pass "orthog() produces valid statistics"
}
else {
    noi test_fail "orthog() values" "invalid cstat or p-value"
}

* Test noid option
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) noid
if _rc == 0 {
    noi test_pass "noid option accepted"
}
else {
    noi test_fail "noid" "returned error `=_rc'"
}

* Test endogtest() with two endogenous variables
drop x_endog y
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
gen y = 2*x_endog1 + 1.5*x_endog2 + 1.0*x_exog + rnormal()

capture civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2)
if _rc == 0 {
    noi test_pass "endogtest() option accepted"
}
else {
    noi test_fail "endogtest()" "returned error `=_rc'"
}

* Test endogtest() stores results
civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2)
if e(endogtest) != . & e(endogtest_df) != . & e(endogtest_p) != . {
    noi test_pass "endogtest() stores endogtest, endogtest_df, endogtest_p"
}
else {
    noi test_fail "endogtest() results" "missing stored results"
}

* Test redundant() option
capture civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3)
if _rc == 0 {
    noi test_pass "redundant() option accepted"
}
else {
    noi test_fail "redundant()" "returned error `=_rc'"
}

* Test redundant() stores results
civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3)
if e(redund) != . & e(redund_df) != . & e(redund_p) != . {
    noi test_pass "redundant() stores redund, redund_df, redund_p"
}
else {
    noi test_fail "redundant() results" "missing stored results"
}

/*******************************************************************************
 * SECTION 19: partial() option (FWL partialling)
 ******************************************************************************/
noi print_section "partial() option"

* Create test data
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog1 = runiform()
gen x_exog2 = rnormal()
gen y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + rnormal()

* Test basic partial() works
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
if _rc == 0 {
    noi test_pass "partial() option accepted"
}
else {
    noi test_fail "partial()" "returned error `=_rc'"
}

* Test that partial() excludes variable from coefficient table
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
capture local test_coef = _b[x_exog2]
if _rc != 0 {
    noi test_pass "partial() excludes variable from coefficient table"
}
else {
    noi test_fail "partial() exclusion" "partialled variable still in output"
}

* Test FWL theorem: coefficients on non-partialled vars should match
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm)
local b_full = _b[x_endog]
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
local b_partial = _b[x_endog]
local diff = abs(`b_full' - `b_partial')
if `diff' < 0.01 {
    noi test_pass "partial() coefficients match full model (FWL theorem)"
}
else {
    noi test_fail "partial() FWL" "coefficient diff = `diff' (should be < 0.01)"
}

* Test partial() with multiple variables
gen x_exog3 = rnormal()
replace y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + 0.3*x_exog3 + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2 x_exog3, absorb(firm) partial(x_exog2 x_exog3)
if _rc == 0 {
    noi test_pass "partial() with multiple variables"
}
else {
    noi test_fail "partial() multiple vars" "returned error `=_rc'"
}

* Test error handling: partial non-existent variable
capture civreghdfe y (x_endog = z1 z2) x_exog1, absorb(firm) partial(nonexistent)
if _rc != 0 {
    noi test_pass "partial() rejects non-existent variable"
}
else {
    noi test_fail "partial() error handling" "should reject non-existent variable"
}

* Test error handling: partial endogenous variable
capture civreghdfe y (x_endog = z1 z2) x_exog1, absorb(firm) partial(x_endog)
if _rc != 0 {
    noi test_pass "partial() rejects endogenous variable"
}
else {
    noi test_fail "partial() error handling" "should reject endogenous variable"
}

/*******************************************************************************
 * SECTION 20: ffirst option (extended first-stage)
 ******************************************************************************/
noi print_section "ffirst option"

* Create test data
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.0*x_exog + rnormal()

* Test basic ffirst works
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst
if _rc == 0 {
    noi test_pass "ffirst option accepted"
}
else {
    noi test_fail "ffirst" "returned error `=_rc'"
}

* Test ffirst stores partial R²
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst
if e(partial_r2_1) != . {
    noi test_pass "ffirst stores partial_r2_1"
}
else {
    noi test_fail "ffirst partial_r2" "partial_r2_1 not stored"
}

* Test ffirst with multiple endogenous vars
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
replace y = 2*x_endog + 1.5*x_endog2 + 1.0*x_exog + rnormal()
civreghdfe y (x_endog x_endog2 = z1 z2 z3) x_exog, absorb(firm) ffirst
if e(partial_r2_1) != . & e(partial_r2_2) != . {
    noi test_pass "ffirst stores partial_r2 for multiple endogenous"
}
else {
    noi test_fail "ffirst multiple endogenous" "partial_r2 not stored correctly"
}

/*******************************************************************************
 * SECTION 21: Save options (savefirst, saverf)
 ******************************************************************************/
noi print_section "Save options"

* Create test data
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.0*x_exog + rnormal()

* Test savefirst option
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst
if _rc == 0 {
    noi test_pass "savefirst option accepted"
}
else {
    noi test_fail "savefirst" "returned error `=_rc'"
}

* Test savefirst stores estimates
capture estimates restore _civreghdfe_main
if _rc == 0 {
    noi test_pass "savefirst stores estimates correctly"
}
else {
    noi test_fail "savefirst storage" "could not restore saved estimates"
}

* Test savefirst with custom prefix
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst savefprefix(test_)
if _rc == 0 {
    capture estimates restore test_main
    if _rc == 0 {
        noi test_pass "savefprefix() with custom prefix"
    }
    else {
        noi test_fail "savefprefix()" "custom prefix not applied"
    }
}
else {
    noi test_fail "savefprefix()" "returned error `=_rc'"
}

* Test saverf option
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) saverf
if _rc == 0 {
    noi test_pass "saverf option accepted"
}
else {
    noi test_fail "saverf" "returned error `=_rc'"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "civreghdfe"

if $TESTS_FAILED > 0 {
    exit 1
}

}
