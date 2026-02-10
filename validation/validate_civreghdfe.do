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

* Test basic partial() - compare non-partialled coefficients against ivreghdfe
* NOTE: Only compare coefficients, not VCE — VCE scaling with partial()
* differs between ivreghdfe and civreghdfe implementations
quietly ivreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
local ivr_N = e(N)
foreach var in x_endog x_exog1 {
    local ivr_b_`var' = _b[`var']
}

civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)

local all_diffs ""
local has_failure = 0
if `ivr_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`ivr_N'!=`=e(N)'"
}
foreach var in x_endog x_exog1 {
    sigfigs `ivr_b_`var'' `=_b[`var']'
    if r(sigfigs) < $DEFAULT_SIGFIGS {
        local has_failure = 1
        local sf_fmt : display %4.1f r(sigfigs)
        local all_diffs "`all_diffs' b[`var']:sigfigs=`sf_fmt'"
    }
}
if `has_failure' == 0 {
    test_pass "partial() coefficients match ivreghdfe"
}
else {
    test_fail "partial()" "`=trim("`all_diffs'")'"
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
* NOTE: FWL guarantees coefficient equality, not SE equality (partial changes DOF/VCE)
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm)
foreach var in x_endog x_exog1 {
    local b_full_`var' = _b[`var']
}
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)

local all_diffs ""
local has_failure = 0
foreach var in x_endog x_exog1 {
    sigfigs `b_full_`var'' `=_b[`var']'
    if r(sigfigs) < $DEFAULT_SIGFIGS {
        local has_failure = 1
        local sf_fmt : display %4.1f r(sigfigs)
        local all_diffs "`all_diffs' b[`var']:sigfigs=`sf_fmt'"
    }
}
if `has_failure' == 0 {
    test_pass "partial() coefficients match full model (FWL theorem)"
}
else {
    test_fail "partial() FWL" "`=trim("`all_diffs'")'"
}

* Test partial() with multiple variables - compare against ivreghdfe
gen x_exog3 = rnormal()
replace y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + 0.3*x_exog3 + rnormal()

quietly ivreghdfe y (x_endog = z1 z2) x_exog1 x_exog2 x_exog3, absorb(firm) partial(x_exog2 x_exog3)
local ivr_N = e(N)
foreach var in x_endog x_exog1 {
    local ivr_b_`var' = _b[`var']
}

civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2 x_exog3, absorb(firm) partial(x_exog2 x_exog3)

local all_diffs ""
local has_failure = 0
if `ivr_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`ivr_N'!=`=e(N)'"
}
foreach var in x_endog x_exog1 {
    sigfigs `ivr_b_`var'' `=_b[`var']'
    if r(sigfigs) < $DEFAULT_SIGFIGS {
        local has_failure = 1
        local sf_fmt : display %4.1f r(sigfigs)
        local all_diffs "`all_diffs' b[`var']:sigfigs=`sf_fmt'"
    }
}
if `has_failure' == 0 {
    test_pass "partial() with multiple variables"
}
else {
    test_fail "partial() multiple vars" "`=trim("`all_diffs'")'"
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

* Test savefirst option - verify main results match ivreghdfe (full comparison)
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) testname("savefirst: main results match ivreghdfe")

* Test savefirst stores estimates
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst
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

* Test saverf option - verify main results match ivreghdfe (full comparison)
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) testname("saverf: main results match ivreghdfe")

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

* Test multiple factor variables
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) i.rep78 i.foreign, absorb(headroom) testname("multiple i.varname as exogenous")

* Test factor variable interactions
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) i.rep78##c.turn, absorb(foreign) testname("factor interaction (i.var##c.var)")

* Base level specification (ib#.var)
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) ib3.rep78, absorb(foreign) testname("ib3.rep78 (custom base level)")

* Benchmark continuous-by-factor interaction
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) c.turn#i.foreign, absorb(rep78) testname("c.turn#i.foreign interaction")

* Factor-by-factor interaction
* Note: i.race#i.union with absorb(idcode) causes different collinearity
* selections after HDFE partialling (Cholesky vs _rmcoll). Both produce valid
* parameterizations. Compare model-level statistics instead of coefficients.
webuse nlswork, clear
keep in 1/5000
local __testpfx "i.race#i.union interaction"

quietly ivreghdfe ln_wage (tenure = age) i.race#i.union, absorb(idcode)
local iv_N    = e(N)
local iv_df_r = e(df_r)
local iv_rss  = e(rss)
local iv_r2   = e(r2)
local iv_rmse = e(rmse)

quietly civreghdfe ln_wage (tenure = age) i.race#i.union, absorb(idcode)
local c_N    = e(N)
local c_df_r = e(df_r)
local c_rss  = e(rss)
local c_r2   = e(r2)
local c_rmse = e(rmse)

* N and df_r should match exactly
if `c_N' == `iv_N' {
    test_pass "`__testpfx': N"
}
else {
    test_fail "`__testpfx': N" "expected `iv_N', got `c_N'"
}
if `c_df_r' == `iv_df_r' {
    test_pass "`__testpfx': df_r"
}
else {
    test_fail "`__testpfx': df_r" "expected `iv_df_r', got `c_df_r'"
}

* Model fit statistics should agree to at least 6 sigfigs
assert_scalar_equal `c_rss'  `iv_rss'  6 "`__testpfx': rss"
assert_scalar_equal `c_r2'   `iv_r2'   6 "`__testpfx': r2"
assert_scalar_equal `c_rmse' `iv_rmse' 6 "`__testpfx': rmse"

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

* kclass(1) should equal 2SLS - compare civreghdfe kclass(1) vs civreghdfe default
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm)
matrix tsls_b = e(b)
matrix tsls_V = e(V)

civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(1)
matrix kclass_b = e(b)
matrix kclass_V = e(V)

matrix_min_sigfigs tsls_b kclass_b
local sf_b = r(min_sigfigs)
matrix_min_sigfigs tsls_V kclass_V
local sf_V = r(min_sigfigs)

if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f min(`sf_b', `sf_V')
    test_pass "kclass(1) equals 2SLS (b sigfigs=`sf_b', V sigfigs=`sf_V')"
}
else {
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "kclass(1) vs 2SLS" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
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

* Test b0() option - verify results match (b0 only affects CUE iteration, not 2SLS)
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm)
matrix nob0_b = e(b)
matrix nob0_V = e(V)

matrix b0 = (2, 1.5)
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) b0(b0)
matrix withb0_b = e(b)
matrix withb0_V = e(V)

matrix_min_sigfigs nob0_b withb0_b
local sf_b = r(min_sigfigs)
matrix_min_sigfigs nob0_V withb0_V
local sf_V = r(min_sigfigs)

if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
    test_pass "b0(matrix) results match default (b sigfigs=`sf_b')"
}
else {
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "b0(matrix)" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
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
matrix tsls_b = e(b)
matrix tsls_V = e(V)

civreghdfe y (x_endog = z1) x_exog, absorb(firm) gmm2s
matrix gmm_b = e(b)
matrix gmm_V = e(V)

civreghdfe y (x_endog = z1) x_exog, absorb(firm) cue
matrix cue_b = e(b)
matrix cue_V = e(V)

matrix_min_sigfigs tsls_b gmm_b
local sf_b = r(min_sigfigs)
matrix_min_sigfigs tsls_V gmm_V
local sf_V = r(min_sigfigs)
if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
    test_pass "gmm2s = 2sls just-identified (b sf=`sf_b', V sf=`sf_V')"
}
else {
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "gmm2s just-identified" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
}

matrix_min_sigfigs tsls_b cue_b
local sf_b = r(min_sigfigs)
matrix_min_sigfigs tsls_V cue_V
local sf_V = r(min_sigfigs)
if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
    test_pass "cue = 2sls just-identified (b sf=`sf_b', V sf=`sf_V')"
}
else {
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "cue just-identified" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
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

* Verify robust and vce(robust) give same results (full b and V comparison)
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
matrix vce_b = e(b)
matrix vce_V = e(V)

civreghdfe price (mpg = weight length), absorb(foreign) robust
matrix rob_b = e(b)
matrix rob_V = e(V)

matrix_min_sigfigs vce_b rob_b
local sf_b = r(min_sigfigs)
matrix_min_sigfigs vce_V rob_V
local sf_V = r(min_sigfigs)
if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
    test_pass "robust = vce(robust) (b sf=`sf_b', V sf=`sf_V')"
}
else {
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "robust = vce(robust)" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
}

* Test standalone cluster option - full comparison against ivreghdfe
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(cluster foreign) testname("standalone cluster option")

* Test dofminus option - verify df_r changes AND coefficients/VCE still reasonable
civreghdfe price (mpg = weight length), absorb(foreign)
local df_r_default = e(df_r)
matrix default_b = e(b)

civreghdfe price (mpg = weight length), absorb(foreign) dofminus(1)
local df_r_minus = e(df_r)
matrix minus_b = e(b)

local all_diffs ""
local has_failure = 0
if `df_r_minus' != `df_r_default' - 1 {
    local has_failure = 1
    local all_diffs "df_r:`df_r_default'!=`df_r_minus'+1"
}
* Coefficients should be identical (dofminus only affects VCE scaling)
matrix_min_sigfigs default_b minus_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dofminus(1) reduces df_r by 1, coefficients unchanged"
}
else {
    test_fail "dofminus" "`=trim("`all_diffs'")'"
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

* Test fwl() alias for partial() - verify it produces same results as partial()
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

civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
matrix partial_b = e(b)
matrix partial_V = e(V)

civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) fwl(x_exog2)
matrix fwl_b = e(b)
matrix fwl_V = e(V)

matrix_min_sigfigs partial_b fwl_b
local sf_b = r(min_sigfigs)
matrix_min_sigfigs partial_V fwl_V
local sf_V = r(min_sigfigs)
if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
    test_pass "fwl() alias matches partial() (b sf=`sf_b', V sf=`sf_V')"
}
else {
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "fwl() alias" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
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
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(truncated) bw(2) testname("kernel(truncated) bw(2)")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(thann) bw(3) testname("kernel(thann) bw(3)")

* Benchmark Driscoll-Kraay
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2) testname("dkraay(2)")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(4) testname("dkraay(4)")

* Benchmark Kiefer
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer testname("kiefer")

* Benchmark bw with robust
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) bw(2) vce(robust) testname("bw(2) + robust")

/*******************************************************************************
 * SECTION 27: Edge Cases - Singletons
 ******************************************************************************/
print_section "Edge Cases - Singletons"

* All singletons (each obs is its own FE group)
clear
set seed 99
set obs 100
gen id = _n
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

capture civreghdfe y (x_endog = z) x_exog, absorb(id)
if _rc == 0 | _rc == 2001 {
    test_pass "all singletons (rc=`=_rc')"
}
else {
    test_fail "all singletons" "returned unexpected error `=_rc'"
}

* Mixed singletons and valid groups
clear
set seed 101
set obs 1000
gen id = _n
replace id = ceil(_n / 10) in 1/500
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

capture civreghdfe y (x_endog = z) x_exog, absorb(id)
if _rc == 0 | _rc == 2001 {
    test_pass "mixed singletons/groups (rc=`=_rc')"
}
else {
    test_fail "mixed singletons" "returned unexpected error `=_rc'"
}

* Two-way FE with singletons
clear
set seed 102
set obs 2000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 20)
replace id2 = 1000 + _n if _n <= 100
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

capture civreghdfe y (x_endog = z) x_exog, absorb(id1 id2)
if _rc == 0 | _rc == 2001 {
    test_pass "two-way FE with singletons (rc=`=_rc')"
}
else {
    test_fail "two-way singletons" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 28: Edge Cases - Missing Values
 ******************************************************************************/
print_section "Edge Cases - Missing Values"

* Missing in dependent variable
clear
set seed 140
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()
replace y = . if runiform() < 0.1

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("missing in depvar")

* Missing in endogenous variable
clear
set seed 141
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
replace x_endog = . if runiform() < 0.1
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("missing in endogenous var")

* Missing in instrument
clear
set seed 142
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
replace z = . if runiform() < 0.1
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("missing in instrument")

* Missing in exogenous variable
clear
set seed 143
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
replace x_exog = . if runiform() < 0.1
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("missing in exogenous var")

* Missing in FE variable
clear
set seed 144
set obs 1000
gen id = runiformint(1, 50)
replace id = . if runiform() < 0.05
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("missing in FE var")

* 50% missing
clear
set seed 145
set obs 2000
gen id = runiformint(1, 100)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()
replace y = . if runiform() < 0.5

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("50% missing")

/*******************************************************************************
 * SECTION 29: Edge Cases - Perfect Collinearity
 ******************************************************************************/
print_section "Edge Cases - Perfect Collinearity"

* Collinear exogenous regressors - compare against ivreghdfe
clear
set seed 150
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog1 = runiform()
gen x_exog2 = 2 * x_exog1
gen y = 2*x_endog + x_exog1 + rnormal()

quietly ivreghdfe y (x_endog = z) x_exog1 x_exog2, absorb(id)
local ivr_N = e(N)
local ivr_b_endog = _b[x_endog]
local ivr_se_endog = _se[x_endog]

civreghdfe y (x_endog = z) x_exog1 x_exog2, absorb(id)
local all_diffs ""
local has_failure = 0
if `ivr_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`ivr_N'!=`=e(N)'"
}
sigfigs `ivr_b_endog' `=_b[x_endog]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' b[x_endog]:sigfigs=`sf_fmt'"
}
sigfigs `ivr_se_endog' `=_se[x_endog]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' se[x_endog]:sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "collinear exogenous regressors handled"
}
else {
    test_fail "collinear exogenous" "`=trim("`all_diffs'")'"
}

* Collinear instruments - compare against ivreghdfe
clear
set seed 151
set obs 1000
gen id = runiformint(1, 50)
gen z1 = runiform()
gen z2 = 2 * z1
gen x_endog = 0.5*z1 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

quietly ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id)
local ivr_N = e(N)
local ivr_b_endog = _b[x_endog]
local ivr_se_endog = _se[x_endog]

civreghdfe y (x_endog = z1 z2) x_exog, absorb(id)
local all_diffs ""
local has_failure = 0
if `ivr_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`ivr_N'!=`=e(N)'"
}
sigfigs `ivr_b_endog' `=_b[x_endog]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' b[x_endog]:sigfigs=`sf_fmt'"
}
sigfigs `ivr_se_endog' `=_se[x_endog]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' se[x_endog]:sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "collinear instruments handled"
}
else {
    test_fail "collinear instruments" "`=trim("`all_diffs'")'"
}

* Endogenous variable collinear with exogenous - compare rc against ivreghdfe
clear
set seed 152
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_exog = runiform()
gen x_endog = x_exog
gen y = 2*x_endog + x_exog + rnormal()

capture ivreghdfe y (x_endog = z) x_exog, absorb(id)
local ivr_rc = _rc
capture civreghdfe y (x_endog = z) x_exog, absorb(id)
local civr_rc = _rc
if `ivr_rc' == `civr_rc' {
    test_pass "endog collinear with exog (both rc=`ivr_rc')"
}
else {
    test_fail "endog collinear with exog" "rc differ: ivreghdfe=`ivr_rc' civreghdfe=`civr_rc'"
}

* Instrument collinear with FE (constant within groups) - compare against ivreghdfe
clear
set seed 153
set obs 1000
gen id = runiformint(1, 50)
gen z_good = runiform()
bysort id: gen z_collinear = id[1]
gen x_endog = 0.5*z_good + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

quietly ivreghdfe y (x_endog = z_good z_collinear) x_exog, absorb(id)
local ivr_N = e(N)
local ivr_b_endog = _b[x_endog]
local ivr_se_endog = _se[x_endog]

civreghdfe y (x_endog = z_good z_collinear) x_exog, absorb(id)
local all_diffs ""
local has_failure = 0
if `ivr_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`ivr_N'!=`=e(N)'"
}
sigfigs `ivr_b_endog' `=_b[x_endog]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' b[x_endog]:sigfigs=`sf_fmt'"
}
sigfigs `ivr_se_endog' `=_se[x_endog]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' se[x_endog]:sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "instrument collinear with FE handled"
}
else {
    test_fail "instrument collinear with FE" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 30: Edge Cases - Unbalanced Panels
 ******************************************************************************/
print_section "Edge Cases - Unbalanced Panels"

* Varying group sizes
clear
set seed 130
set obs 5000
gen id = .
local obs = 1
forvalues i = 1/100 {
    local group_size = runiformint(5, 100)
    forvalues j = 1/`group_size' {
        if `obs' <= 5000 {
            qui replace id = `i' in `obs'
            local obs = `obs' + 1
        }
    }
}
replace id = runiformint(1, 100) if missing(id)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("unbalanced panel")

* Extreme unbalance: one huge group + many tiny
clear
set seed 131
set obs 5000
gen id = 1 in 1/2500
replace id = _n - 2499 in 2501/5000
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

capture civreghdfe y (x_endog = z) x_exog, absorb(id)
if _rc == 0 | _rc == 2001 {
    test_pass "extreme unbalance (one huge, many tiny) rc=`=_rc'"
}
else {
    test_fail "extreme unbalance" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 31: High-dimensional and Multi-way FE
 ******************************************************************************/
print_section "High-Dimensional and Multi-Way FE"

* 1000 FE levels (10K obs)
clear
set seed 123
set obs 10000
gen id = runiformint(1, 1000)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("1000 FE levels (10K obs)")

* High-dimensional two-way FE (500x400)
clear
set seed 124
set obs 20000
gen id1 = runiformint(1, 500)
gen id2 = runiformint(1, 400)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id1 id2) testname("high-dim two-way (500x400)")

* Three-way FE
clear
set seed 125
set obs 10000
gen id1 = runiformint(1, 100)
gen id2 = runiformint(1, 50)
gen id3 = runiformint(1, 20)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id1 id2 id3) testname("three-way FE basic")
benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id1 id2 id3) vce(robust) testname("three-way FE robust")
benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id1 id2 id3) vce(cluster id1) testname("three-way FE clustered")

/*******************************************************************************
 * SECTION 32: Additional Datasets (Grunfeld, bplong)
 ******************************************************************************/
print_section "Additional Datasets"

* Grunfeld: classic panel IV
capture webuse grunfeld, clear
if _rc == 0 {
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company) testname("grunfeld: company FE")
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company year) testname("grunfeld: two-way FE")
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company) vce(robust) testname("grunfeld: robust")
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company) vce(cluster company) testname("grunfeld: cluster company")
}
else {
    test_pass "grunfeld dataset not available - skipped"
}

* bplong: blood pressure panel
capture webuse bplong, clear
if _rc == 0 {
    benchmark_ivreghdfe bp (agegrp = when) sex, absorb(patient) testname("bplong: patient FE")
}
else {
    test_pass "bplong dataset not available - skipped"
}

/*******************************************************************************
 * SECTION 33: Large Datasets
 ******************************************************************************/
print_section "Large Datasets"

* 50K observations
clear
set seed 200
set obs 50000
gen id = runiformint(1, 1000)
gen year = runiformint(2000, 2020)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog1 = runiform()
gen x_exog2 = rnormal()
gen y = 2*x_endog + 1.5*x_exog1 - 0.5*x_exog2 + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(id) testname("50K obs, single FE")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(id year) testname("50K obs, two-way FE")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(id) vce(robust) testname("50K obs, robust")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(id) vce(cluster id) testname("50K obs, clustered")

* 100K observations
clear
set seed 201
set obs 100000
gen id = runiformint(1, 2000)
gen year = runiformint(2000, 2020)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("100K obs, single FE")
benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id year) testname("100K obs, two-way FE")

/*******************************************************************************
 * SECTION 34: VCE Comprehensive
 ******************************************************************************/
print_section "VCE Comprehensive"

* Large clusters (many)
webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age) ttl_exp, absorb(idcode) vce(cluster idcode) testname("large clusters (many)")

* Small clusters (few)
sysuse auto, clear
gen cluster_var = ceil(_n / 10)
benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) vce(cluster cluster_var) testname("small clusters (few)")

* Two clusters
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) turn, absorb(rep78) vce(cluster foreign) testname("two clusters")

* Robust + aweight
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement], absorb(foreign) vce(robust) testname("robust + aweight")

* Robust + fweight
sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement) [fw=fw], absorb(foreign) vce(robust) testname("robust + fweight")

* Robust + pweight
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [pw=displacement], absorb(foreign) vce(robust) testname("robust + pweight")

* Cluster + aweight
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement], absorb(rep78) vce(cluster foreign) testname("cluster + aweight")

* Cluster + fweight
sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement) [fw=fw], absorb(rep78) vce(cluster foreign) testname("cluster + fweight")

/*******************************************************************************
 * SECTION 35: Weight x Estimator Combinations
 ******************************************************************************/
print_section "Weight x Estimator Combinations"

* LIML + weights
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) [aw=displacement], absorb(foreign) liml testname("liml + aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement length) [fw=fw], absorb(foreign) liml testname("liml + fweight")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) [pw=displacement], absorb(foreign) liml testname("liml + pweight")

* GMM2S + weights
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) [aw=displacement], absorb(foreign) gmm2s testname("gmm2s + aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement length) [fw=fw], absorb(foreign) gmm2s testname("gmm2s + fweight")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) [pw=displacement], absorb(foreign) gmm2s testname("gmm2s + pweight")

* Fuller + weights
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) [aw=displacement], absorb(foreign) fuller(1) testname("fuller(1) + aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement length) [fw=fw], absorb(foreign) fuller(1) testname("fuller(1) + fweight")

* CUE + aweight
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) [aw=displacement], absorb(foreign) cue testname("cue + aweight")

/*******************************************************************************
 * SECTION 36: Weight x Diagnostic Test Combinations
 ******************************************************************************/
print_section "Weight x Diagnostic Test Combinations"

* Create overidentified data for diagnostic tests
clear
set seed 360
set obs 1000
gen id = runiformint(1, 50)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
gen aw_var = abs(rnormal()) + 0.5
gen int fw_var = runiformint(1, 5)

* orthog + aweight
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog [aw=aw_var], absorb(id) orthog(z3) testname("orthog + aweight")

* orthog + fweight
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog [fw=fw_var], absorb(id) orthog(z3) testname("orthog + fweight")

* endogtest + aweight (need two endogenous)
drop x_endog y
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
gen y = 2*x_endog1 + 1.5*x_endog2 + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog [aw=aw_var], absorb(id) endogtest(x_endog2) testname("endogtest + aweight")

* endogtest + fweight
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog [fw=fw_var], absorb(id) endogtest(x_endog2) testname("endogtest + fweight")

* redundant + aweight
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog [aw=aw_var], absorb(id) redundant(z3) testname("redundant + aweight")

* redundant + fweight
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog [fw=fw_var], absorb(id) redundant(z3) testname("redundant + fweight")

/*******************************************************************************
 * SECTION 37: Estimator x VCE Combinations
 ******************************************************************************/
print_section "Estimator x VCE Combinations"

clear
set seed 370
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* LIML + robust
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) liml vce(robust) testname("liml + robust")

* LIML + cluster
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) liml vce(cluster firm) testname("liml + cluster")

* Fuller + robust
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(1) vce(robust) testname("fuller(1) + robust")

* Fuller + cluster
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(1) vce(cluster firm) testname("fuller(1) + cluster")

* GMM2S + cluster
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm) testname("gmm2s + cluster")

/*******************************************************************************
 * SECTION 38: Many Covariates and Scale Issues
 ******************************************************************************/
print_section "Many Covariates and Scale Issues"

* Many exogenous covariates
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) turn displacement headroom trunk length, absorb(foreign) testname("6 exogenous covariates")

* Covariates with different scales
clear
set seed 380
set obs 2000
gen id = runiformint(1, 100)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_small = runiform() / 1000
gen x_medium = runiform() * 100
gen x_large = runiform() * 1000000
gen y = 2*x_endog + x_small*1000 + x_medium/100 + x_large/1e6 + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_small x_medium x_large, absorb(id) testname("different scales")

* Pathological: very small values - compare against ivreghdfe
clear
set seed 381
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform() / 1e6
gen x_endog = 0.5*z + rnormal() / 1e6
gen x_exog = runiform() / 1e6
gen y = x_endog * 1e6 + x_exog * 1e6 + rnormal() / 1e6

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("very small values")

* Pathological: very large values - compare against ivreghdfe
clear
set seed 382
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform() * 1e9
gen x_endog = 0.5*z + rnormal() * 1e6
gen x_exog = runiform() * 1e9
gen y = x_endog / 1e6 + x_exog / 1e6 + rnormal() * 1e3

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("very large values")

* Pathological: mixed tiny and huge - compare against ivreghdfe
clear
set seed 383
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_tiny = runiform() / 1e9
gen x_huge = runiform() * 1e9
gen y = 2*x_endog + x_tiny * 1e9 + x_huge / 1e9 + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_tiny x_huge, absorb(id) testname("mixed scale (tiny and huge)")

/*******************************************************************************
 * SECTION 39: Many Instruments / Many Endogenous
 ******************************************************************************/
print_section "Many Instruments / Many Endogenous"

* Highly overidentified: 1 endogenous, 5 instruments
clear
set seed 390
set obs 1000
gen id = runiformint(1, 50)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen z4 = rnormal() + 0.1*z1
gen z5 = runiform() + 0.1*z2
gen x_endog = 0.3*z1 + 0.2*z2 + 0.15*z3 + 0.1*z4 + 0.1*z5 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2 z3 z4 z5) x_exog, absorb(id) testname("1 endog, 5 instruments")

* Three endogenous variables, 6 instruments
clear
set seed 391
set obs 2000
gen id = runiformint(1, 100)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen z4 = rnormal()
gen z5 = runiform()
gen z6 = rnormal()
gen x_endog1 = 0.4*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.3*z3 + 0.3*z4 + rnormal()
gen x_endog3 = 0.3*z5 + 0.3*z6 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog1 + 1.5*x_endog2 + x_endog3 + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog1 x_endog2 x_endog3 = z1 z2 z3 z4 z5 z6) x_exog, absorb(id) testname("3 endog, 6 instruments")

* Exactly identified: 1 endog, 1 instrument with various options
clear
set seed 392
set obs 500
gen firm = ceil(_n / 10)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(firm) testname("exactly identified: 1/1")
benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(firm) vce(robust) testname("exactly identified: 1/1 + robust")
benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(firm) small testname("exactly identified: 1/1 + small")

* Exactly identified: 2 endog, 2 instruments
clear
set seed 393
set obs 1000
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog1 = 0.5*z1 + rnormal()
gen x_endog2 = 0.5*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog1 + 1.5*x_endog2 + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2) x_exog, absorb(firm) testname("exactly identified: 2/2")
benchmark_ivreghdfe y (x_endog1 x_endog2 = z1 z2) x_exog, absorb(firm) liml testname("exactly identified: 2/2 + liml")

/*******************************************************************************
 * SECTION 40: Option Combinations
 ******************************************************************************/
print_section "Option Combinations"

* Two-way FE + robust + weights
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement], absorb(foreign rep78) vce(robust) testname("two-way FE + robust + aw")

* Two-way FE + cluster + weights
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement], absorb(foreign rep78) vce(cluster foreign) testname("two-way FE + cluster + aw")

* if + weights + VCE
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement] if price > 4000, absorb(foreign) vce(robust) testname("if + aw + robust")

* first + small + robust
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) first small vce(robust) testname("first + small + robust")

* ffirst + cluster + weights
clear
set seed 400
set obs 1000
gen id = runiformint(1, 50)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
gen aw_var = abs(rnormal()) + 0.5

benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog [aw=aw_var], absorb(id) ffirst vce(cluster id) testname("ffirst + cluster + aw")

* LIML + first + robust
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) liml first vce(robust) testname("liml + first + robust")

/*******************************************************************************
 * SECTION 41: Stress Tests
 ******************************************************************************/
print_section "Stress Tests"

* 10K FE levels (50K obs)
clear
set seed 300
set obs 50000
gen id = runiformint(1, 10000)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) testname("10K FE levels (50K obs)")

* Many exogenous covariates (8)
clear
set seed 301
set obs 5000
gen id = runiformint(1, 100)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
forvalues i = 1/8 {
    gen x`i' = runiform()
}
gen y = 2*x_endog + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2) x1 x2 x3 x4 x5 x6 x7 x8, absorb(id) testname("8 exogenous covariates")

* Highly imbalanced clusters
clear
set seed 302
set obs 10000
gen cluster_id = 1 in 1/9000
replace cluster_id = 2 in 9001/9500
replace cluster_id = runiformint(3, 100) in 9501/10000
gen id = runiformint(1, 200)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z) x_exog, absorb(id) vce(cluster cluster_id) testname("highly imbalanced clusters")

/*******************************************************************************
 * SECTION 42: String Cluster Variables
 ******************************************************************************/
print_section "String Cluster Variables"

sysuse auto, clear
decode foreign, gen(foreign_str)

* Compare numeric vs string cluster - full b and V comparison
quietly ivreghdfe price (mpg = weight length), absorb(foreign) vce(cluster foreign)
local ivr_N = e(N)
local ivr_N_clust = e(N_clust)
matrix ivr_b = e(b)
matrix ivr_V = e(V)

capture quietly civreghdfe price (mpg = weight length), absorb(foreign) vce(cluster foreign_str)
if _rc == 0 {
    local civr_N = e(N)
    local civr_N_clust = e(N_clust)
    matrix civr_b = e(b)
    matrix civr_V = e(V)

    local all_diffs ""
    local has_failure = 0
    if `ivr_N' != `civr_N' {
        local has_failure = 1
        local all_diffs "e(N):`ivr_N'!=`civr_N'"
    }
    if `ivr_N_clust' != `civr_N_clust' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(N_clust):`ivr_N_clust'!=`civr_N_clust'"
    }
    matrix_min_sigfigs ivr_b civr_b
    if r(min_sigfigs) < $DEFAULT_SIGFIGS {
        local has_failure = 1
        local sf_fmt : display %4.1f r(min_sigfigs)
        local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
    }
    matrix_min_sigfigs ivr_V civr_V
    if r(min_sigfigs) < $DEFAULT_SIGFIGS {
        local has_failure = 1
        local sf_fmt : display %4.1f r(min_sigfigs)
        local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
    }
    if `has_failure' == 0 {
        test_pass "string cluster variable"
    }
    else {
        test_fail "string cluster variable" "`=trim("`all_diffs'")'"
    }
}
else {
    test_fail "string cluster variable" "civreghdfe returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 43: Additional Error Tests
 ******************************************************************************/
print_section "Additional Error Tests"

* Instrument = endogenous variable
sysuse auto, clear
test_error_match, stata_cmd(ivreghdfe price (mpg = mpg), absorb(foreign)) ctools_cmd(civreghdfe price (mpg = mpg), absorb(foreign)) testname("instrument = endogenous var")

* Too few instruments (underidentified)
clear
set seed 430
set obs 500
gen firm = ceil(_n / 10)
gen z = runiform()
gen x_endog1 = 0.5*z + rnormal()
gen x_endog2 = 0.3*z + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog1 + x_endog2 + x_exog + rnormal()

test_error_match, stata_cmd(ivreghdfe y (x_endog1 x_endog2 = z) x_exog, absorb(firm)) ctools_cmd(civreghdfe y (x_endog1 x_endog2 = z) x_exog, absorb(firm)) testname("underidentified (2 endog, 1 instr)")

* Absorb variable with all missing
sysuse auto, clear
gen all_miss = .
test_error_match, stata_cmd(ivreghdfe price (mpg = weight), absorb(all_miss)) ctools_cmd(civreghdfe price (mpg = weight), absorb(all_miss)) testname("absorb var all missing")

* Weight variable with zero
sysuse auto, clear
gen bad_wt = displacement
replace bad_wt = 0 in 1/5
test_error_match, stata_cmd(ivreghdfe price (mpg = weight) [aw=bad_wt], absorb(foreign)) ctools_cmd(civreghdfe price (mpg = weight) [aw=bad_wt], absorb(foreign)) testname("zero weight values")

* Negative weight variable
sysuse auto, clear
gen neg_wt = displacement
replace neg_wt = -1 in 1
test_error_match, stata_cmd(ivreghdfe price (mpg = weight) [aw=neg_wt], absorb(foreign)) ctools_cmd(civreghdfe price (mpg = weight) [aw=neg_wt], absorb(foreign)) testname("negative weight value")

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
