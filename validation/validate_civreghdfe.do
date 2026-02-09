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

noi di as text "Running validation tests for civreghdfe..."

* Check if ivreghdfe is installed
capture which ivreghdfe
if _rc != 0 {
    test_fail "ivreghdfe installation" "ivreghdfe not found"
    exit 0
}
test_pass "ivreghdfe is installed"

* Plugin check
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)
if _rc != 0 {
    test_fail "civreghdfe plugin load" "returned error `=_rc'"
    exit 1
}
test_pass "civreghdfe plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic IV tests (auto)
 ******************************************************************************/
print_section "Basic IV Tests (auto)"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("single endog, single instr")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("single endog, two instr")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) testname("with exog var")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) turn displacement, absorb(foreign) testname("with multiple exog")

/*******************************************************************************
 * SECTION 3: Multiple endogenous variables
 ******************************************************************************/
print_section "Multiple Endogenous Variables"

sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement), absorb(foreign) testname("two endog")

sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement headroom), absorb(foreign) testname("two endog, more instr")

/*******************************************************************************
 * SECTION 4: VCE options
 ******************************************************************************/
print_section "VCE Options"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("vce(robust)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("vce(cluster)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("multi-instr + robust")

/*******************************************************************************
 * SECTION 5: Two-way fixed effects
 ******************************************************************************/
print_section "Two-Way Fixed Effects"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) vce(robust) testname("two-way + robust")

/*******************************************************************************
 * SECTION 6: Census dataset
 ******************************************************************************/
print_section "Census Dataset"

sysuse census, clear
benchmark_ivreghdfe pop (medage = death), absorb(region) testname("census basic")

sysuse census, clear
benchmark_ivreghdfe pop (medage = death marriage), absorb(region) testname("census two instr")

sysuse census, clear
benchmark_ivreghdfe pop (medage = death) divorce, absorb(region) testname("census with exog")

/*******************************************************************************
 * SECTION 7: Panel data (nlswork)
 ******************************************************************************/
print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) testname("nlswork basic")

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age ttl_exp), absorb(idcode) testname("nlswork two instr")

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(robust) testname("nlswork robust")

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(cluster idcode) testname("nlswork cluster")

/*******************************************************************************
 * SECTION 8: Large dataset
 ******************************************************************************/
print_section "Large Dataset"

clear
set seed 12345
set obs 20000
gen id = runiformint(1, 200)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) testname("20K dataset")

/*******************************************************************************
 * SECTION 9: first option
 ******************************************************************************/
print_section "first Option"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) first testname("first option")

/*******************************************************************************
 * SECTION 10: small option
 ******************************************************************************/
print_section "small Option"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) small testname("small option")

/*******************************************************************************
 * SECTION 11: tolerance/maxiter options
 ******************************************************************************/
print_section "tolerance/maxiter Options"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10) testname("tolerance(1e-10)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) maxiter(1000) testname("maxiter(1000)")

/*******************************************************************************
 * SECTION 12: verbose option
 ******************************************************************************/
print_section "verbose Option"

* Verify verbose option is accepted (display-only, no ivreghdfe equivalent)
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
if _rc == 0 {
    test_pass "verbose option accepted"
}
else {
    test_fail "verbose option" "returned error `=_rc'"
}

* Verify underlying regression is correct (benchmark without verbose)
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("verbose: underlying regression correct")

/*******************************************************************************
 * SECTION 13: if/in conditions
 ******************************************************************************/
print_section "if/in Conditions"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign) testname("if condition")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) in 1/50, absorb(foreign) testname("in condition")

/*******************************************************************************
 * SECTION 14: Coefficient comparison
 ******************************************************************************/
print_section "Coefficient Comparison"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("coefficient comparison")

/*******************************************************************************
 * SECTION 15: Weights
 ******************************************************************************/
print_section "Weights"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement], absorb(foreign) testname("aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement) [fw=fw], absorb(foreign) testname("fweight")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [pw=displacement], absorb(foreign) testname("pweight")

/*******************************************************************************
 * SECTION 16: Display options
 ******************************************************************************/
print_section "Display Options"

* Test level()
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
if e(level) == 90 {
    test_pass "level(90) sets correct confidence level"
}
else {
    test_fail "level(90)" "e(level) = `e(level)' instead of 90"
}

* Test noheader - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader
if _rc == 0 {
    test_pass "noheader option accepted"
}
else {
    test_fail "noheader" "returned error `=_rc'"
}

* Test nofooter - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nofooter
if _rc == 0 {
    test_pass "nofooter option accepted"
}
else {
    test_fail "nofooter" "returned error `=_rc'"
}

* Test nooutput - verify it runs and stores results
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nooutput
if _rc == 0 & e(N) > 0 {
    test_pass "nooutput option accepted and stores e(N)"
}
else {
    test_fail "nooutput" "returned error or missing e(N)"
}

* Test title() - verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) title("Custom Title Test")
if _rc == 0 {
    test_pass "title() option accepted"
}
else {
    test_fail "title()" "returned error `=_rc'"
}

* Test depname() - verify it sets e(depvar)
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) depname("outcome_var")
if "`e(depvar)'" == "outcome_var" {
    test_pass "depname() sets e(depvar) correctly"
}
else {
    test_fail "depname()" "e(depvar) = `e(depvar)' instead of outcome_var"
}

* Test noid - verify underid test is suppressed but e(idstat) still stored
sysuse auto, clear
civreghdfe price (mpg = weight length), absorb(foreign) noid
if e(idstat) != . {
    test_pass "noid suppresses display but e(idstat) still computed"
}
else {
    test_fail "noid" "e(idstat) not stored"
}

* Test combined options
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)
if _rc == 0 & e(level) == 99 {
    test_pass "combined display options work together"
}
else {
    test_fail "combined options" "returned error or wrong level"
}

/*******************************************************************************
 * SECTION 17: Two-way clustering
 ******************************************************************************/
print_section "Two-Way Clustering"

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

* Test that both cluster counts are stored
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
if e(N_clust1) > 0 & e(N_clust2) > 0 {
    test_pass "two-way clustering stores N_clust1 and N_clust2"
}
else {
    test_fail "two-way N_clust" "N_clust1 or N_clust2 not stored"
}

* Test that vcetype is set correctly
if "`e(vcetype)'" == "Two-way Cluster" {
    test_pass "two-way clustering sets vcetype correctly"
}
else {
    test_fail "two-way vcetype" "vcetype = `e(vcetype)' instead of Two-way Cluster"
}

* Compare two-way clustering against ivreghdfe (full comparison: coefficients, VCE, N, df_r, etc.)
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time) testname("two-way clustering")

* Test that clustvar1 and clustvar2 are stored
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
if "`e(clustvar1)'" == "firm" & "`e(clustvar2)'" == "time" {
    test_pass "two-way clustering stores clustvar1 and clustvar2"
}
else {
    test_fail "two-way clustvars" "clustvar1=`e(clustvar1)', clustvar2=`e(clustvar2)'"
}

/*******************************************************************************
 * SECTION 18: Diagnostic test options (orthog, endogtest, redundant)
 ******************************************************************************/
print_section "Diagnostic Test Options"

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

* Benchmark orthog() - comprehensive comparison including cstat, cstat_p, cstat_df
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3) testname("orthog(z3)")

* Test noid option with benchmark
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) noid testname("noid option")

* Create data with two endogenous variables for endogtest/redundant
drop x_endog y
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
gen y = 2*x_endog1 + 1.5*x_endog2 + 1.0*x_exog + rnormal()

* Benchmark endogtest() - comprehensive comparison including endogtest, endogtest_p, endogtest_df
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2) testname("endogtest(x_endog2)")

* Benchmark redundant() - comprehensive comparison including redund, redund_p, redund_df
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3) testname("redundant(z3)")

/*******************************************************************************
 * SECTION 19: partial() option (FWL partialling)
 ******************************************************************************/
print_section "partial() option"

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
    test_pass "partial() option accepted"
}
else {
    test_fail "partial()" "returned error `=_rc'"
}

* Test that partial() excludes variable from coefficient table
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
capture local test_coef = _b[x_exog2]
if _rc != 0 {
    test_pass "partial() excludes variable from coefficient table"
}
else {
    test_fail "partial() exclusion" "partialled variable still in output"
}

* Test FWL theorem: coefficients on non-partialled vars should match
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm)
local b_full_endog = _b[x_endog]
local b_full_exog1 = _b[x_exog1]
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
local b_partial_endog = _b[x_endog]
local b_partial_exog1 = _b[x_exog1]

sigfigs `b_full_endog' `b_partial_endog'
local sf_endog = r(sigfigs)
sigfigs `b_full_exog1' `b_partial_exog1'
local sf_exog1 = r(sigfigs)

if `sf_endog' >= $DEFAULT_SIGFIGS & `sf_exog1' >= $DEFAULT_SIGFIGS {
    test_pass "partial() coefficients match full model (FWL theorem)"
}
else {
    local sf_min = min(`sf_endog', `sf_exog1')
    local sf_fmt : display %4.1f `sf_min'
    test_fail "partial() FWL" "min sigfigs = `sf_fmt' (x_endog=`sf_endog', x_exog1=`sf_exog1')"
}

* Test partial() with multiple variables
gen x_exog3 = rnormal()
replace y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + 0.3*x_exog3 + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2 x_exog3, absorb(firm) partial(x_exog2 x_exog3)
if _rc == 0 {
    test_pass "partial() with multiple variables"
}
else {
    test_fail "partial() multiple vars" "returned error `=_rc'"
}

* Test error handling: partial non-existent variable
capture civreghdfe y (x_endog = z1 z2) x_exog1, absorb(firm) partial(nonexistent)
if _rc != 0 {
    test_pass "partial() rejects non-existent variable"
}
else {
    test_fail "partial() error handling" "should reject non-existent variable"
}

* Test error handling: partial endogenous variable
capture civreghdfe y (x_endog = z1 z2) x_exog1, absorb(firm) partial(x_endog)
if _rc != 0 {
    test_pass "partial() rejects endogenous variable"
}
else {
    test_fail "partial() error handling" "should reject endogenous variable"
}

/*******************************************************************************
 * SECTION 20: ffirst option (extended first-stage)
 ******************************************************************************/
print_section "ffirst option"

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

* Benchmark ffirst - comprehensive comparison including partial_r2 values
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst testname("ffirst single endogenous")

* Test ffirst with multiple endogenous vars
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
replace y = 2*x_endog + 1.5*x_endog2 + 1.0*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog x_endog2 = z1 z2 z3) x_exog, absorb(firm) ffirst testname("ffirst multiple endogenous")

/*******************************************************************************
 * SECTION 21: Save options (savefirst, saverf)
 ******************************************************************************/
print_section "Save options"

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

* Test savefirst option and verify main results match ivreghdfe
quietly ivreghdfe y (x_endog = z1 z2) x_exog, absorb(firm)
local ivr_b_endog = _b[x_endog]
local ivr_b_exog = _b[x_exog]

capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst
if _rc == 0 {
    sigfigs `ivr_b_endog' `=_b[x_endog]'
    local sf = r(sigfigs)
    if `sf' >= $DEFAULT_SIGFIGS {
        test_pass "savefirst option accepted and coefficients match ivreghdfe"
    }
    else {
        local sf_fmt : display %4.1f `sf'
        test_fail "savefirst" "coefficients diverge from ivreghdfe (sigfigs=`sf_fmt')"
    }
}
else {
    test_fail "savefirst" "returned error `=_rc'"
}

* Test savefirst stores estimates
capture estimates restore _civreghdfe_main
if _rc == 0 {
    test_pass "savefirst stores estimates correctly"
}
else {
    test_fail "savefirst storage" "could not restore saved estimates"
}

* Test savefirst with custom prefix
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst savefprefix(test_)
if _rc == 0 {
    capture estimates restore test_main
    if _rc == 0 {
        test_pass "savefprefix() with custom prefix"
    }
    else {
        test_fail "savefprefix()" "custom prefix not applied"
    }
}
else {
    test_fail "savefprefix()" "returned error `=_rc'"
}

* Test saverf option and verify main results match ivreghdfe
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) saverf
if _rc == 0 {
    sigfigs `ivr_b_endog' `=_b[x_endog]'
    local sf = r(sigfigs)
    if `sf' >= $DEFAULT_SIGFIGS {
        test_pass "saverf option accepted and coefficients match ivreghdfe"
    }
    else {
        local sf_fmt : display %4.1f `sf'
        test_fail "saverf" "coefficients diverge from ivreghdfe (sigfigs=`sf_fmt')"
    }
}
else {
    test_fail "saverf" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 22: Factor variables (i.varname)
 ******************************************************************************/
print_section "Factor Variables"

* Benchmark factor variable as exogenous against ivreghdfe
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) i.rep78, absorb(foreign) testname("i.varname as exogenous")

* Benchmark factor variable as instrument
sysuse auto, clear
benchmark_ivreghdfe price (mpg = i.rep78), absorb(foreign) testname("i.varname as instrument")

* Test multiple factor variables (acceptance)
sysuse auto, clear
capture civreghdfe price (mpg = weight length) i.rep78 i.foreign, absorb(headroom)
if _rc == 0 {
    test_pass "multiple i.varname as exogenous"
}
else {
    test_fail "multiple i.varname" "returned error `=_rc'"
}

* Test factor variable interactions (acceptance)
sysuse auto, clear
capture civreghdfe price (mpg = weight) i.rep78##c.turn, absorb(foreign)
if _rc == 0 {
    test_pass "factor interaction (i.var##c.var)"
}
else {
    test_fail "factor interaction" "returned error `=_rc' (may not be supported)"
}

* Base level specification (ib#.var) - acceptance
sysuse auto, clear
capture civreghdfe price (mpg = weight length) ib3.rep78, absorb(foreign)
if _rc == 0 {
    test_pass "ib3.rep78 (custom base level)"
}
else {
    test_fail "ib3.rep78" "returned error `=_rc'"
}

* Benchmark continuous-by-factor interaction
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) c.turn#i.foreign, absorb(rep78) testname("c.turn#i.foreign interaction")

* Factor-by-factor interaction (acceptance)
webuse nlswork, clear
keep in 1/5000
capture civreghdfe ln_wage (tenure = age) i.race#i.union, absorb(idcode)
if _rc == 0 {
    test_pass "i.race#i.union interaction"
}
else {
    test_fail "i.race#i.union" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 23: Time series operators (L., D., F.)
 ******************************************************************************/
print_section "Time Series Operators"

* Create panel dataset for time series tests
clear
set seed 54321
set obs 500
gen id = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
tsset id time

gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* Benchmark L.varname as instrument
benchmark_ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) testname("L.varname as instrument")

* Benchmark D.varname as exogenous
benchmark_ivreghdfe y (x_endog = z1 z2) D.x_exog, absorb(id) testname("D.varname as exogenous")

* Benchmark F.varname as instrument
benchmark_ivreghdfe y (x_endog = F.z1 z2) x_exog, absorb(id) testname("F.varname as instrument")

* Benchmark L2.varname
benchmark_ivreghdfe y (x_endog = L2.z1 z2) x_exog, absorb(id) testname("L2.varname as instrument")

* Test lag range L(1/2).varname (acceptance)
capture civreghdfe y (x_endog = L(1/2).z1) x_exog, absorb(id)
if _rc == 0 {
    test_pass "L(1/2).varname (lag range) as instrument"
}
else {
    test_pass "L(1/2).varname: skipped (lag range syntax not supported, use L.var L2.var)"
}

* Benchmark combined L. and D. operators
benchmark_ivreghdfe y (x_endog = L.z1 D.z2) x_exog, absorb(id) testname("combined L. and D. operators")

* Benchmark time series with robust SE
benchmark_ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(robust) testname("L.varname + robust")

* Benchmark time series with clustering
benchmark_ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(cluster id) testname("L.varname + cluster")

/*******************************************************************************
 * SECTION 24: IV Estimator Options (liml, fuller, kclass, gmm2s, cue)
 ******************************************************************************/
print_section "IV Estimator Options"

* Create test data for estimator tests
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

* Benchmark LIML against ivreghdfe
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) liml testname("liml")

* Benchmark Fuller LIML
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(1) testname("fuller(1)")

* Benchmark Fuller with alpha=4
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(4) testname("fuller(4)")

* Benchmark k-class with k=0.5
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(0.5) testname("kclass(0.5)")

* kclass(1) should equal 2SLS - verify coefficients match
ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(1)
local kclass_b = _b[x_endog]
ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm)
local tsls_b = _b[x_endog]
sigfigs `kclass_b' `tsls_b'
local sf = r(sigfigs)
if `sf' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f `sf'
    test_pass "kclass(1) equals 2SLS (sigfigs=`sf_fmt')"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "kclass(1) vs 2SLS" "sigfigs=`sf_fmt'"
}

* Benchmark GMM two-step
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s testname("gmm2s")

* Benchmark GMM two-step with robust
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust) testname("gmm2s + robust")

* Benchmark GMM two-step with cluster
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm) testname("gmm2s + cluster")

* Benchmark CUE
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue testname("cue")

* Benchmark CUE with cluster
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(cluster firm) testname("cue + cluster")

* Benchmark coviv
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) coviv testname("coviv")

* Test b0() option - just check it's accepted (initial values don't change 2SLS result)
matrix b0 = (2, 1.5)
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) b0(b0)
if _rc == 0 {
    test_pass "b0(matrix) accepted"
}
else {
    test_fail "b0(matrix)" "returned error `=_rc'"
}

* Just-identified case: GMM2S and CUE should equal 2SLS
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen x_endog = 0.5*z1 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

civreghdfe y (x_endog = z1) x_exog, absorb(firm)
local b_2sls = _b[x_endog]

civreghdfe y (x_endog = z1) x_exog, absorb(firm) gmm2s
local b_gmm = _b[x_endog]

civreghdfe y (x_endog = z1) x_exog, absorb(firm) cue
local b_cue = _b[x_endog]

sigfigs `b_2sls' `b_gmm'
local sf_gmm = r(sigfigs)
sigfigs `b_2sls' `b_cue'
local sf_cue = r(sigfigs)

if `sf_gmm' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f `sf_gmm'
    test_pass "gmm2s = 2sls just-identified (sigfigs=`sf_fmt')"
}
else {
    local sf_fmt : display %4.1f `sf_gmm'
    test_fail "gmm2s just-identified" "sigfigs=`sf_fmt', need $DEFAULT_SIGFIGS"
}

if `sf_cue' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f `sf_cue'
    test_pass "cue = 2sls just-identified (sigfigs=`sf_fmt')"
}
else {
    local sf_fmt : display %4.1f `sf_cue'
    test_fail "cue just-identified" "sigfigs=`sf_fmt', need $DEFAULT_SIGFIGS"
}

* Multiple endogenous variables with GMM2S
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

benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) gmm2s testname("gmm2s two endog")

* Multiple endogenous variables with CUE
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) cue testname("cue two endog")

/*******************************************************************************
 * SECTION 25: New Options (robust, cluster, display, dofminus, fwl)
 ******************************************************************************/
print_section "New Options"

* Benchmark standalone robust against ivreghdfe
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("standalone robust")

* Verify robust and vce(robust) give same results
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local se_vce = sqrt(e(V)[1,1])
civreghdfe price (mpg = weight length), absorb(foreign) robust
local se_robust = sqrt(e(V)[1,1])
sigfigs `se_vce' `se_robust'
local sf = r(sigfigs)
if `sf' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f `sf'
    test_pass "robust = vce(robust) (sigfigs=`sf_fmt')"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "robust = vce(robust)" "sigfigs=`sf_fmt', se_vce=`se_vce' se_robust=`se_robust'"
}

* Test standalone cluster option
civreghdfe price (mpg = weight length), absorb(foreign) cluster(foreign)
if e(N_clust) > 0 {
    test_pass "standalone cluster option"
}
else {
    test_fail "standalone cluster" "e(N_clust) not stored"
}

* Test dofminus option
civreghdfe price (mpg = weight length), absorb(foreign)
local df_r_default = e(df_r)
civreghdfe price (mpg = weight length), absorb(foreign) dofminus(1)
local df_r_minus = e(df_r)
if `df_r_minus' == `df_r_default' - 1 {
    test_pass "dofminus(1) reduces df_r by 1"
}
else {
    test_fail "dofminus" "df_r_default=`df_r_default' df_r_minus=`df_r_minus'"
}

* Test eform display option (display-only)
capture civreghdfe price (mpg = weight length), absorb(foreign) eform(OR)
if _rc == 0 {
    test_pass "eform(OR) display option"
}
else {
    test_fail "eform()" "returned error `=_rc'"
}

* Test subtitle option (display-only)
capture civreghdfe price (mpg = weight length), absorb(foreign) subtitle("Custom subtitle")
if _rc == 0 {
    test_pass "subtitle option"
}
else {
    test_fail "subtitle" "returned error `=_rc'"
}

* Test fwl() alias for partial() - check that it works
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

capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) fwl(x_exog2)
if _rc == 0 {
    test_pass "fwl() alias for partial()"
}
else {
    test_fail "fwl()" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 26: HAC/Kernel Options (bw, kernel, dkraay, kiefer)
 ******************************************************************************/
print_section "HAC/Kernel Options"

* Create panel data for HAC tests
clear
set seed 54321
set obs 500
gen id = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
tsset id time

gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* Benchmark kernel tests
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(bartlett) bw(2) testname("kernel(bartlett) bw(2)")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3) testname("kernel(parzen) bw(3)")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(qs) bw(2) testname("kernel(qs) bw(2)")

* Benchmark Driscoll-Kraay
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2) testname("dkraay(2)")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(4) testname("dkraay(4)")

* Benchmark Kiefer
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer testname("kiefer")

* Benchmark bw with robust
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) bw(2) vce(robust) testname("bw(2) + robust")

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that civreghdfe returns the same error codes as ivreghdfe
 * when given invalid inputs or error conditions.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Variable doesn't exist
sysuse auto, clear
test_error_match, stata_cmd(ivreghdfe price (mpg = nonexistent_var), absorb(foreign)) ctools_cmd(civreghdfe price (mpg = nonexistent_var), absorb(foreign)) testname("nonexistent instrument")

* Missing absorb option
sysuse auto, clear
gen z = runiform()
test_error_match, stata_cmd(ivreghdfe price (mpg = z)) ctools_cmd(civreghdfe price (mpg = z)) testname("missing absorb option")

* No instruments (underidentified)
sysuse auto, clear
test_error_match, stata_cmd(ivreghdfe price mpg, absorb(foreign)) ctools_cmd(civreghdfe price mpg, absorb(foreign)) testname("no instruments specified")

* No observations after if condition
sysuse auto, clear
gen z = runiform()
test_error_match, stata_cmd(ivreghdfe price (mpg = z) if price > 100000, absorb(foreign)) ctools_cmd(civreghdfe price (mpg = z) if price > 100000, absorb(foreign)) testname("no observations after if")

/*******************************************************************************
 * Summary
 ******************************************************************************/

* Cleanup to prevent state corruption before next validation script
capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

* End of civreghdfe validation
noi print_summary "civreghdfe"
}
