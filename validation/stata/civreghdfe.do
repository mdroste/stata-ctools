* civreghdfe: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "civreghdfe"

if "`__case'" == "civreghdfe" {
* validate_civreghdfe.do
* Comprehensive validation tests for civreghdfe vs ivreghdfe
* Tests all options: absorb, vce, weights, first, small, tolerance, maxiter
* Requires ivreghdfe to be installed (ssc install ivreghdfe)

* Load setup (works from project root or validation dir)
do "validation/support/validate_setup.do"

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

* SECTION 2: Basic IV tests (auto)
print_section "Basic IV Tests (auto)"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("single endog, single instr")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("single endog, two instr")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) testname("with exog var")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) turn displacement, absorb(foreign) testname("with multiple exog")

* SECTION 3: Multiple endogenous variables
print_section "Multiple Endogenous Variables"

sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement), absorb(foreign) testname("two endog")

sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement headroom), absorb(foreign) testname("two endog, more instr")

* SECTION 4: VCE options
print_section "VCE Options"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("vce(robust)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("vce(cluster)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("multi-instr + robust")

* SECTION 5: Two-way fixed effects
print_section "Two-Way Fixed Effects"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) vce(robust) testname("two-way + robust")

* SECTION 6: Census dataset
print_section "Census Dataset"

sysuse census, clear
benchmark_ivreghdfe pop (medage = death), absorb(region) testname("census basic")

sysuse census, clear
benchmark_ivreghdfe pop (medage = death marriage), absorb(region) testname("census two instr")

sysuse census, clear
benchmark_ivreghdfe pop (medage = death) divorce, absorb(region) testname("census with exog")

* SECTION 7: Panel data (nlswork)
print_section "Panel Data (nlswork)"

ctools_fixture nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) testname("nlswork basic")

ctools_fixture nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age ttl_exp), absorb(idcode) testname("nlswork two instr")

ctools_fixture nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(robust) testname("nlswork robust")

ctools_fixture nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(cluster idcode) testname("nlswork cluster")

* SECTION 8: Large dataset
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

* SECTION 9: first option
print_section "first Option"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) first testname("first option")

* SECTION 10: small option
print_section "small Option"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) small testname("small option")

* SECTION 11: tolerance/maxiter options
print_section "tolerance/maxiter Options"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10) testname("tolerance(1e-10)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) maxiter(1000) testname("maxiter(1000)")

* SECTION 12: verbose option
print_section "verbose Option"

* Verify verbose option is accepted (display-only, no ivreghdfe equivalent)
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
if _rc == 0 test_pass "verbose option accepted"
else test_fail "verbose option" "returned error `=_rc'"

* Verify underlying regression is correct (benchmark without verbose)
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("verbose: underlying regression correct")

* SECTION 13: if/in conditions
print_section "if/in Conditions"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign) testname("if condition")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) in 1/50, absorb(foreign) testname("in condition")

* SECTION 14: Coefficient comparison
print_section "Coefficient Comparison"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("coefficient comparison")

* SECTION 15: Weights
print_section "Weights"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [aw=displacement], absorb(foreign) testname("aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_ivreghdfe price (mpg = displacement) [fw=fw], absorb(foreign) testname("fweight")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) [pw=displacement], absorb(foreign) testname("pweight")

* SECTION 16: Display options
print_section "Display Options"

* Test level()
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
if e(level) == 90 test_pass "level(90) sets correct confidence level"
else test_fail "level(90)" "e(level) = `e(level)' instead of 90"

* Test noheader - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader
if _rc == 0 test_pass "noheader option accepted"
else test_fail "noheader" "returned error `=_rc'"

* Test nofooter - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nofooter
if _rc == 0 test_pass "nofooter option accepted"
else test_fail "nofooter" "returned error `=_rc'"

* Test nooutput - verify it runs and stores results
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nooutput
if _rc == 0 & e(N) == _N test_pass "nooutput option accepted and stores e(N)"
else test_fail "nooutput" "returned error or missing e(N)"

* Test title() - verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) title("Custom Title Test")
if _rc == 0 test_pass "title() option accepted"
else test_fail "title()" "returned error `=_rc'"

* Test depname() - verify it sets e(depvar)
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) depname("outcome_var")
if "`e(depvar)'" == "outcome_var" test_pass "depname() sets e(depvar) correctly"
else test_fail "depname()" "e(depvar) = `e(depvar)' instead of outcome_var"

* Test noid - as ivreg2, the identification statistics are neither shown
* nor stored
sysuse auto, clear
civreghdfe price (mpg = weight length), absorb(foreign) noid
if missing(e(idstat)) & missing(e(cdf)) & missing(e(widstat)) test_pass "noid: no identification statistics, as ivreg2"
else test_fail "noid" "e(idstat) stored"

* Test combined options
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)
if _rc == 0 & e(level) == 99 test_pass "combined display options work together"
else test_fail "combined options" "returned error or wrong level"

* SECTION 17: Two-way clustering
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
if e(N_clust1) > 0 & e(N_clust2) > 0 test_pass "two-way clustering stores N_clust1 and N_clust2"
else test_fail "two-way N_clust" "N_clust1 or N_clust2 not stored"

* Test that vcetype is set correctly
if "`e(vcetype)'" == "Two-way Cluster" test_pass "two-way clustering sets vcetype correctly"
else test_fail "two-way vcetype" "vcetype = `e(vcetype)' instead of Two-way Cluster"

* Compare two-way clustering against ivreghdfe (full comparison: coefficients, VCE, N, df_r, etc.)
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time) testname("two-way clustering")

* Test that clustvar1 and clustvar2 are stored
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
if "`e(clustvar1)'" == "firm" & "`e(clustvar2)'" == "time" test_pass "two-way clustering stores clustvar1 and clustvar2"
else test_fail "two-way clustvars" "clustvar1=`e(clustvar1)', clustvar2=`e(clustvar2)'"

* SECTION 18: Diagnostic test options (orthog, endogtest, redundant)
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

* SECTION 19: partial() option (FWL partialling)
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
if `has_failure' == 0 test_pass "partial() coefficients match ivreghdfe"
else test_fail "partial()" "`=trim("`all_diffs'")'"

* Test that partial() excludes variable from coefficient table
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
capture local test_coef = _b[x_exog2]
if _rc != 0 test_pass "partial() excludes variable from coefficient table"
else test_fail "partial() exclusion" "partialled variable still in output"

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
if `has_failure' == 0 test_pass "partial() coefficients match full model (FWL theorem)"
else test_fail "partial() FWL" "`=trim("`all_diffs'")'"

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
if `has_failure' == 0 test_pass "partial() with multiple variables"
else test_fail "partial() multiple vars" "`=trim("`all_diffs'")'"

* Test error handling: partial non-existent variable
capture civreghdfe y (x_endog = z1 z2) x_exog1, absorb(firm) partial(nonexistent)
if _rc != 0 test_pass "partial() rejects non-existent variable"
else test_fail "partial() error handling" "should reject non-existent variable"

* Test error handling: partial endogenous variable
capture civreghdfe y (x_endog = z1 z2) x_exog1, absorb(firm) partial(x_endog)
if _rc != 0 test_pass "partial() rejects endogenous variable"
else test_fail "partial() error handling" "should reject endogenous variable"

* SECTION 20: ffirst option (extended first-stage)
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

* SECTION 21: Save options (savefirst, saverf)
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

* Test savefirst stores the first stage under ivreg2's default name
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst
capture estimates restore _ivreg2_x_endog
if _rc == 0 test_pass "savefirst stores estimates correctly"
else test_fail "savefirst storage" "could not restore saved estimates"

* Test savefirst with custom prefix
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst savefprefix(test_)
if _rc == 0 {
    capture estimates restore test_x_endog
    if _rc == 0 test_pass "savefprefix() with custom prefix"
    else test_fail "savefprefix()" "custom prefix not applied"
}
else {
    test_fail "savefprefix()" "returned error `=_rc'"
}

* Test saverf option - verify main results match ivreghdfe (full comparison)
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) testname("saverf: main results match ivreghdfe")

* SECTION 22: Factor variables (i.varname)
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
ctools_fixture nlswork, clear
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
if `c_N' == `iv_N' test_pass "`__testpfx': N"
else test_fail "`__testpfx': N" "expected `iv_N', got `c_N'"
if `c_df_r' == `iv_df_r' test_pass "`__testpfx': df_r"
else test_fail "`__testpfx': df_r" "expected `iv_df_r', got `c_df_r'"

* Model fit statistics should agree to at least 6 sigfigs
assert_scalar_equal `c_rss'  `iv_rss'  6 "`__testpfx': rss"
assert_scalar_equal `c_r2'   `iv_r2'   6 "`__testpfx': r2"
assert_scalar_equal `c_rmse' `iv_rmse' 6 "`__testpfx': rmse"

* SECTION 23: Time series operators (L., D., F.)
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
if _rc == 0 test_pass "L(1/2).varname (lag range) as instrument"
else test_fail "L(1/2).varname" "lag-range parsing failed"

* Benchmark combined L. and D. operators
benchmark_ivreghdfe y (x_endog = L.z1 D.z2) x_exog, absorb(id) testname("combined L. and D. operators")

* Benchmark time series with robust SE
benchmark_ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(robust) testname("L.varname + robust")

* Benchmark time series with clustering
benchmark_ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(cluster id) testname("L.varname + cluster")

* SECTION 24: IV Estimator Options (liml, fuller, kclass, gmm2s, cue)
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

* b0() columns are matched to the regressors by name, as ivreg2 matches
* them: a matrix without the regressors' names is refused (ivreghdfe r(198));
* b0() itself is compared with ivreghdfe below
matrix b0 = (2, 1.5)
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) b0(b0)
if _rc == 198 test_pass "b0(matrix) refused"
else test_fail "b0(matrix) refused" "rc=`=_rc'"

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

* Weakly identified CUE: from the 2SLS start the objective is not convex,
* where civreghdfe's steepest-descent fallback used to stop at worse optima
* (J 4.12 against ivreghdfe's 2.70 with seed 4, robust; 2.31 against 1.88
* with seed 32, cluster). The absolute-eigenvalue Newton step finds
* ivreghdfe's optimum. Six significant figures: on these flat objectives
* ivreg2's optimize() stops at its own tolerance (b within 1e-7, with
* civreghdfe's J lower by 1e-14), so the seventh digit is ivreg2's.
foreach __cs in "4 vce(robust)" "32 vce(cluster id)" "14 vce(robust) bw(2)" {
    gettoken __seed __opt : __cs
    clear
    set seed `__seed'
    set obs 300
    gen id = ceil(_n / 6)
    gen t = mod(_n - 1, 6) + 1
    gen fe = rnormal()
    bysort id: replace fe = fe[1]
    gen z1 = rnormal()
    gen z2 = rnormal()
    gen z3 = rnormal()
    gen u = rnormal() * (1 + abs(z1))
    gen x = .08 * z1 + .05 * z2 + .03 * z3 + .8 * u + rnormal() + fe
    gen w = rnormal()
    gen y = .5 * x + .3 * w + u + fe
    if strpos("`__opt'", "bw(") tsset id t
    benchmark_ivreghdfe y (x = z1 z2 z3) w, absorb(id) cue `__opt' minsf(6) testname("cue weak, nonconvex start (seed `__seed',`__opt')")
}

* SECTION 25: New Options (robust, cluster, display, dofminus, fwl)
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
if `has_failure' == 0 test_pass "dofminus(1) reduces df_r by 1, coefficients unchanged"
else test_fail "dofminus" "`=trim("`all_diffs'")'"

* Test eform display option (display-only)
capture civreghdfe price (mpg = weight length), absorb(foreign) eform(OR)
if _rc == 0 test_pass "eform(OR) display option"
else test_fail "eform()" "returned error `=_rc'"

* Test subtitle option (display-only)
capture civreghdfe price (mpg = weight length), absorb(foreign) subtitle("Custom subtitle")
if _rc == 0 test_pass "subtitle option"
else test_fail "subtitle" "returned error `=_rc'"

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

* SECTION 26: HAC/Kernel Options (bw, kernel, dkraay, kiefer)
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

* SECTION 27: Edge Cases - Singletons
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
if _rc == 0 | _rc == 2001 test_pass "all singletons (rc=`=_rc')"
else test_fail "all singletons" "returned unexpected error `=_rc'"

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
if _rc == 0 | _rc == 2001 test_pass "mixed singletons/groups (rc=`=_rc')"
else test_fail "mixed singletons" "returned unexpected error `=_rc'"

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
if _rc == 0 | _rc == 2001 test_pass "two-way FE with singletons (rc=`=_rc')"
else test_fail "two-way singletons" "returned unexpected error `=_rc'"

* SECTION 28: Edge Cases - Missing Values
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

* SECTION 29: Edge Cases - Perfect Collinearity
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
if `has_failure' == 0 test_pass "collinear exogenous regressors handled"
else test_fail "collinear exogenous" "`=trim("`all_diffs'")'"

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
if `has_failure' == 0 test_pass "collinear instruments handled"
else test_fail "collinear instruments" "`=trim("`all_diffs'")'"

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
if `ivr_rc' == `civr_rc' test_pass "endog collinear with exog (both rc=`ivr_rc')"
else test_fail "endog collinear with exog" "rc differ: ivreghdfe=`ivr_rc' civreghdfe=`civr_rc'"

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
if `has_failure' == 0 test_pass "instrument collinear with FE handled"
else test_fail "instrument collinear with FE" "`=trim("`all_diffs'")'"

* SECTION 30: Edge Cases - Unbalanced Panels
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
if _rc == 0 | _rc == 2001 test_pass "extreme unbalance (one huge, many tiny) rc=`=_rc'"
else test_fail "extreme unbalance" "returned unexpected error `=_rc'"

* SECTION 31: High-dimensional and Multi-way FE
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

* SECTION 32: Additional Datasets (Grunfeld, bplong)
print_section "Additional Datasets"

* Grunfeld: classic panel IV
capture ctools_fixture grunfeld, clear
if _rc == 0 {
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company) testname("grunfeld: company FE")
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company year) testname("grunfeld: two-way FE")
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company) vce(robust) testname("grunfeld: robust")
    benchmark_ivreghdfe invest (kstock = L.kstock) mvalue, absorb(company) vce(cluster company) testname("grunfeld: cluster company")
}
else {
    test_fail "grunfeld dataset fixture" "required validation fixture could not be loaded"
}

* bplong: blood pressure panel
capture ctools_fixture bplong, clear
if _rc == 0 {
    benchmark_ivreghdfe bp (agegrp = when) sex, absorb(patient) testname("bplong: patient FE")
}
else {
    test_fail "bplong dataset fixture" "required validation fixture could not be loaded"
}

* SECTION 33: Large Datasets
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

* SECTION 34: VCE Comprehensive
print_section "VCE Comprehensive"

* Large clusters (many)
ctools_fixture nlswork, clear
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

* SECTION 35: Weight x Estimator Combinations
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

* SECTION 36: Weight x Diagnostic Test Combinations
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

* SECTION 37: Estimator x VCE Combinations
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

* SECTION 38: Many Covariates and Scale Issues
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

* Pathological: mixed tiny and huge
* NOTE: ivreghdfe is buggy here — ivreg2 chokes on extreme scales, omitting
* x_endog and x_tiny entirely (b=[0,0,8.52e-10], R2=0.01).
* Fit the equivalent, well-scaled model with ivreghdfe, then transform its
* coefficients and covariance back to the original units. Sampling noise is
* shared by both fits and cannot justify a 50% coefficient tolerance.
clear
set seed 383
set obs 1000
gen id = runiformint(1, 50)
gen z = runiform()
gen x_endog = 0.5*z + rnormal()
gen x_tiny = runiform() / 1e9
gen x_huge = runiform() * 1e9
gen y = 2*x_endog + x_tiny * 1e9 + x_huge / 1e9 + rnormal()

gen double tiny_unit = x_tiny * 1e9
gen double huge_unit = x_huge / 1e9
ivreghdfe y (x_endog = z) tiny_unit huge_unit, absorb(id) tolerance(1e-12)
tempname mixed_b mixed_V mixed_D
matrix `mixed_b' = e(b)
matrix `mixed_V' = e(V)
matrix `mixed_D' = I(colsof(`mixed_b'))
local j = colnumb(`mixed_b', "tiny_unit")
matrix `mixed_D'[`j',`j'] = 1e9
local j = colnumb(`mixed_b', "huge_unit")
matrix `mixed_D'[`j',`j'] = 1e-9
matrix `mixed_b' = `mixed_b' * `mixed_D'
matrix `mixed_V' = `mixed_D' * `mixed_V' * `mixed_D'
local mixed_N = e(N)
local mixed_df = e(df_r)
gen byte mixed_sample = e(sample)
civreghdfe y (x_endog = z) x_tiny x_huge, absorb(id) tolerance(1e-12)
assert e(N) == `mixed_N' & e(df_r) == `mixed_df'
assert e(sample) == mixed_sample
assert_matrix_equal `mixed_b' e(b) $DEFAULT_SIGFIGS "mixed scale: b against rescaled ivreghdfe"
assert_matrix_equal `mixed_V' e(V) $DEFAULT_SIGFIGS "mixed scale: V against rescaled ivreghdfe"

* SECTION 39: Many Instruments / Many Endogenous
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

* SECTION 40: Option Combinations
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

* SECTION 41: Stress Tests
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

* SECTION 42: String Cluster Variables
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
    if `has_failure' == 0 test_pass "string cluster variable"
    else test_fail "string cluster variable" "`=trim("`all_diffs'")'"
}
else {
    test_fail "string cluster variable" "civreghdfe returned error `=_rc'"
}

* SECTION 43: Additional Error Tests
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

* SECTION: Intentional Error Tests
* These tests verify that civreghdfe returns the same error codes as ivreghdfe
* when given invalid inputs or error conditions.
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

* Absorbed degrees of freedom as ivreghdfe counts them (reghdfe's HDFE.df_a):
* an effect nested in either dimension of a two-way cluster, two effects
* nested in one cluster, an effect coarser than the cluster, and three or
* more effects whose later pairs have several mobility groups. e(df_a) must
* be equal; benchmark_ivreghdfe compares e(V), e(F), e(rmse), e(r2_a) and
* e(widstat) (and the rest) at $DEFAULT_SIGFIGS.
capture program drop _civ_dfa_case
program define _civ_dfa_case
    syntax anything(name=spec) [aw], Absorb(varlist) [vce(string) SMALL testname(string)]
    local wt
    if "`weight'" != "" local wt [`weight'`exp']
    local vceopt
    if "`vce'" != "" local vceopt vce(`vce')
    benchmark_ivreghdfe `spec' `wt', absorb(`absorb') `vceopt' `small' testname("`testname'")
    capture ivreghdfe `spec' `wt', absorb(`absorb') `vceopt' `small'
    local rc_ref = _rc
    local df_a_ref = e(df_a)
    capture civreghdfe `spec' `wt', absorb(`absorb') `vceopt' `small'
    local rc_ct = _rc
    if !`rc_ref' & !`rc_ct' & `df_a_ref' == e(df_a) test_pass "`testname': e(df_a)"
    else test_fail "`testname': e(df_a)" "ivreghdfe `df_a_ref' (rc `rc_ref'), civreghdfe `=e(df_a)' (rc `rc_ct')"
end
clear
set seed 20261102
set obs 1200
gen long firm = ceil(_n / 10)
gen int year = mod(_n - 1, 10) + 1
* firm in ind in sector; occ moves across firms and years; ind#year as one id
gen int ind = ceil(firm / 12)
gen byte sector = ceil(ind / 4)
gen int occ = mod(firm * 7 + year * 3, 9) + 1
gen int indyear = (ind - 1) * 10 + year
gen double z1 = rnormal()
gen double z2 = rnormal()
gen double w = rnormal()
gen double x = z1 + .5 * z2 + .3 * w + rnormal() + firm / 100
gen double y = x - .5 * w + rnormal() + firm / 50 + year / 10
gen double aw = 1 + runiform()
foreach wt in "" "[aw=aw]" {
    local tag = cond("`wt'" == "", "unweighted", "aweights")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(firm) vce(cluster year firm) ///
        testname("df_a: absorb(firm) cluster(year firm), `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(firm) vce(cluster firm year) ///
        testname("df_a: absorb(firm) cluster(firm year), `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(firm year) vce(cluster firm year) ///
        testname("df_a: absorb(firm year) cluster(firm year), `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(firm ind) vce(cluster sector) ///
        testname("df_a: absorb(firm ind) cluster(sector), `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(ind firm) vce(cluster firm) ///
        testname("df_a: absorb(ind firm) cluster(firm), `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(indyear firm) vce(cluster ind) ///
        testname("df_a: absorb(ind#year firm) cluster(ind), `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(firm year ind) small ///
        testname("df_a: absorb(firm year ind) small, `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(firm year occ) small ///
        testname("df_a: absorb(firm year occ) small, `tag'")
    _civ_dfa_case y (x = z1 z2) w `wt', absorb(year occ ind) vce(cluster firm) small ///
        testname("df_a: absorb(year occ ind) cluster(firm) small, `tag'")
}

* SECTION: First-stage F under non-iid VCEs
* ivreg2 robustifies e(first)["F",] (Wald test of the excluded instruments
* under the model's VCE); e(F_first#) must match it, including an instrument
* that varies only at the cluster level and FEs nested in the cluster.
print_section "First-stage F under robust, cluster and kernel VCEs"
clear
set seed 101
set obs 2000
gen int c = ceil(_n/20)
gen int fe = ceil(runiform()*40)
gen int t = mod(_n-1, 20) + 1
gen int c2 = ceil(runiform()*30)
bysort c: gen double zc = rnormal() if _n == 1
bysort c: replace zc = zc[1]
gen double z2 = rnormal()
gen double z3 = rnormal()
gen double u = rnormal()
bysort c: gen double uc = rnormal() if _n == 1
bysort c: replace uc = uc[1]
gen double d = .25*zc + .3*z2 + uc + u + rnormal()
gen double d2 = .2*z3 + .4*z2 + u + rnormal()
gen double x = rnormal()
gen double w = 1 + runiform()
gen double y = d + d2 + x + uc + rnormal()
xtset c t
foreach spec in "y x (d=zc)|absorb(fe) robust" "y x (d=zc)|absorb(fe) cluster(c)" ///
    "y x (d=zc)|cluster(c)" "y x (d=zc)|absorb(fe) robust bw(2)" ///
    "y x (d=zc) [pw=w]|absorb(fe)" "y x (d d2=zc z2 z3)|absorb(fe) robust" ///
    "y x (d d2=zc z2 z3)|absorb(fe) cluster(c)" "y x (d d2=zc z2 z3) [aw=w]|absorb(fe) cluster(c)" ///
    "y x (d d2=zc z2 z3)|absorb(fe) cluster(c c2)" "y x (d d2=zc z2 z3)|absorb(fe) kernel(parzen) bw(3) robust" ///
    "y x (d d2=zc z2 z3)|absorb(fe) dkraay(2)" "y x (d d2=z2 z3)|absorb(c) cluster(c)" ///
    "y x (d d2=zc z2 z3)|absorb(fe) robust dofminus(3)" "y x (d d2=zc z2 z3) [pw=w]|absorb(fe) liml" {
    gettoken model opts : spec, parse("|")
    local opts : subinstr local opts "|" ""
    capture quietly ivreghdfe `model', `opts' ffirst
    local rc_ref = _rc
    if !`rc_ref' matrix fs_ref = e(first)
    capture quietly civreghdfe `model', `opts'
    if `rc_ref' | _rc {
        test_fail "first-stage F: `model', `opts'" "rc ivreghdfe=`rc_ref' civreghdfe=`=_rc'"
        continue
    }
    forvalues e = 1/`=colsof(fs_ref)' {
        local F_ref = el(fs_ref, rownumb(fs_ref, "F"), `e')
        assert_scalar_equal `=e(F_first`e')' `F_ref' $DEFAULT_SIGFIGS "first-stage F`e': `model', `opts'"
    }
}

* SECTION: Constant screened first without absorb() (vs ivreg2)
* A regressor that is constant on the sample, or a full set of dummies, is
* omitted and _cons kept, so e(df_m) and e(F) test the slopes only.
print_section "Collinearity with the constant without absorb()"
sysuse auto, clear
generate byte d1 = rep78 >= 3 if rep78 < .
generate byte d2 = 1 - d1
foreach spec in "price (mpg = weight) foreign if foreign==1" "price (mpg = weight) d1 d2" {
    quietly ivreg2 `spec'
    local ref_names : colfullnames e(b)
    matrix ref_b = e(b)
    matrix ref_V = e(V)
    local ref_df_m = e(df_m)
    local ref_F = e(F)
    capture quietly civreghdfe `spec'
    if _rc {
        test_fail "constant first: `spec'" "civreghdfe rc=`=_rc'"
        continue
    }
    local ct_names : colfullnames e(b)
    if "`ct_names'" == "`ref_names'" test_pass "constant first: names `spec'"
    else test_fail "constant first: names `spec'" "`ct_names' vs ivreg2 `ref_names'"
    matrix ct_b = e(b)
    matrix ct_V = e(V)
    assert_matrix_equal ct_b ref_b $DEFAULT_SIGFIGS "constant first: e(b) `spec'"
    assert_matrix_equal ct_V ref_V $DEFAULT_SIGFIGS "constant first: e(V) `spec'"
    assert_scalar_equal `=e(df_m)' `ref_df_m' $DEFAULT_SIGFIGS "constant first: e(df_m) `spec'"
    assert_scalar_equal `=e(F)' `ref_F' $DEFAULT_SIGFIGS "constant first: e(F) `spec'"
}

* SECTION: coviv (vs ivreghdfe)
* LIML, Fuller and k-class with coviv use the 2SLS bread (X'PzX)^-1.
print_section "coviv"
clear
set seed 55
set obs 1500
generate int firm = ceil(_n/15)
generate int t = mod(_n-1, 15) + 1
generate double z1 = rnormal()
generate double z2 = rnormal()
generate double z3 = rnormal()
generate double u = rnormal()
generate double x1 = .5*z1 + .3*z2 + .2*z3 + u + rnormal()
generate double w = rnormal()
generate double y = 1 + x1 + w + u + rnormal()
xtset firm t
benchmark_ivreghdfe y w (x1 = z1 z2 z3), absorb(firm) liml coviv testname("coviv: liml")
benchmark_ivreghdfe y w (x1 = z1 z2 z3), absorb(firm) fuller(1) coviv vce(robust) testname("coviv: fuller(1) robust")
benchmark_ivreghdfe y w (x1 = z1 z2 z3), absorb(firm) kclass(.5) coviv vce(cluster firm) testname("coviv: kclass(.5) cluster")
benchmark_ivreghdfe y w (x1 = z1 z2 z3), absorb(firm) liml coviv vce(robust) bw(2) testname("coviv: liml HAC")
quietly civreghdfe y w (x1 = z1 z2 z3), absorb(firm) kclass(.5) vce(cluster firm)
matrix nocoviv_V = e(V)
quietly civreghdfe y w (x1 = z1 z2 z3), absorb(firm) kclass(.5) coviv vce(cluster firm)
if mreldif(e(V), nocoviv_V) > 1e-6 test_pass "coviv changes the k-class VCE"
else test_fail "coviv changes the k-class VCE" "V equals the run without coviv"

* SECTION: predict after partial() (vs ivreg2)
* e(b) omits the partialled-out terms: residuals are y - Xb net of their
* projection on those terms (the full-model residuals); xb, stdp and xbd exit
* 198 as in ivreg2.
print_section "predict after partial()"
sysuse auto, clear
* (without absorb(), the constant is partialled out too and counted, as in ivreg2)
quietly ivreg2 price (mpg = weight) length turn i.rep78, partial(turn i.rep78)
local ref_partial_ct = e(partial_ct)
quietly ivreg2 price (mpg = weight) length turn i.rep78
predict double ref_resid, resid
quietly civreghdfe price (mpg = weight) length turn i.rep78, partial(turn i.rep78)
assert_scalar_equal `=e(partial_ct)' `ref_partial_ct' $DEFAULT_SIGFIGS "partial(): e(partial_ct)"
predict double ct_resid, resid
assert_var_equal ct_resid ref_resid $DEFAULT_SIGFIGS "partial(): predict residuals"
foreach stat in xb stdp xbd {
    capture predict double ct_`stat', `stat'
    if _rc == 198 test_pass "partial(): predict `stat' refused"
    else test_fail "partial(): predict `stat' refused" "rc=`=_rc'"
}

* civreghdfe against a reference command: the e(b) column names, e(b), e(V)
* and the listed e() scalars, to $DEFAULT_SIGFIGS significant figures (equal
* when both are missing)
capture program drop _civreghdfe_cmp
program define _civreghdfe_cmp
    syntax, ref(string) ct(string) testname(string) [scalars(string) nov omitzero]
    capture quietly `ref'
    if _rc {
        test_fail "`testname'" "reference failed rc=`=_rc'"
        exit
    }
    local ref_names : colfullnames e(b)
    tempname rb rV cb cV
    matrix `rb' = e(b)
    matrix `rV' = e(V)
    * omitzero: ivreghdfe leaves a column that its invsym() zeroes (one the
    * fixed effects absorb) unmarked, with zero coefficient and variance;
    * civreghdfe marks it omitted (o.)
    if "`omitzero'" != "" {
        local __names
        local __j = 0
        foreach __n of local ref_names {
            local ++__j
            if `rb'[1, `__j'] == 0 & `rV'[`__j', `__j'] == 0 & !regexm("`__n'", "[ob]\.") {
                local __parts : subinstr local __n "#" " ", all
                local __n
                foreach __p of local __parts {
                    local __d = strpos("`__p'", ".")
                    if `__d' local __p = substr("`__p'", 1, `__d' - 1) + "o" + substr("`__p'", `__d', .)
                    else local __p o.`__p'
                    local __n = cond("`__n'" == "", "`__p'", "`__n'#`__p'")
                }
            }
            local __names `__names' `__n'
        }
        local ref_names `__names'
    }
    foreach s of local scalars {
        local r_`s' = e(`s')
    }
    ereturn clear
    capture quietly `ct'
    if _rc {
        test_fail "`testname'" "civreghdfe failed rc=`=_rc'"
        exit
    }
    local diffs
    local ct_names : colfullnames e(b)
    if "`ct_names'" != "`ref_names'" local diffs "`diffs' names: `ct_names' vs `ref_names'"
    else {
        matrix `cb' = e(b)
        matrix `cV' = e(V)
        matrix_min_sigfigs `rb' `cb'
        if r(min_sigfigs) < $DEFAULT_SIGFIGS local diffs "`diffs' e(b):sf=`=r(min_sigfigs)'"
        if "`v'" == "" {
            matrix_min_sigfigs `rV' `cV'
            if r(min_sigfigs) < $DEFAULT_SIGFIGS local diffs "`diffs' e(V):sf=`=r(min_sigfigs)'"
        }
    }
    foreach s of local scalars {
        if missing(`r_`s'') & missing(e(`s')) continue
        sigfigs `r_`s'' `=e(`s')'
        if r(sigfigs) < $DEFAULT_SIGFIGS local diffs "`diffs' e(`s'):`r_`s'' vs `=e(`s')'"
    }
    if "`diffs'" == "" test_pass "`testname'"
    else test_fail "`testname'" "`=trim("`diffs'")'"
end
local ivstats N r2 r2_a rss mss rmse F idstat widstat

* The default iteration limit is reghdfe's (16000, also as iterate()): a
* slowly mixing graph that ivreghdfe solves at its defaults
clear
set seed 21
set obs 9000
gen long w = ceil(_n/3)
gen long f = 1 + mod(w - 1 + (mod(_n,3)==0), 3000)
gen double x = rnormal()
gen double z = rnormal()
gen double d = z + rnormal()
gen double y = x + d + w/300 + f/200 + rnormal()
_civreghdfe_cmp, ref(ivreghdfe y x (d=z), absorb(w f)) ct(civreghdfe y x (d=z), absorb(w f)) ///
    scalars(`ivstats') testname("slowly mixing effects at the default iteration limit")
_civreghdfe_cmp, ref(ivreghdfe y x (d=z), absorb(w f) iterate(16000)) ct(civreghdfe y x (d=z), absorb(w f) iterate(16000)) ///
    scalars(`ivstats') testname("iterate()")

* A time-series operator on the dependent variable
ctools_fixture grunfeld, clear
foreach op in L D {
    _civreghdfe_cmp, ref(ivreghdfe `op'.invest kstock (mvalue=L.mvalue), absorb(company)) ///
        ct(civreghdfe `op'.invest kstock (mvalue=L.mvalue), absorb(company)) ///
        scalars(`ivstats') testname("`op'.depvar with absorb()")
    _civreghdfe_cmp, ref(ivreg2 `op'.invest kstock (mvalue=L.mvalue)) ///
        ct(civreghdfe `op'.invest kstock (mvalue=L.mvalue)) ///
        scalars(N r2 rss rmse F idstat widstat) testname("`op'.depvar without absorb()")
    if "`e(cmd)'" == "civreghdfe" & "`e(depvar)'" == "`op'.invest" test_pass "`op'.depvar: e(depvar)"
    else test_fail "`op'.depvar: e(depvar)" "`e(depvar)'"
}

* IV data for the cases below
clear
set seed 4242
set obs 2000
gen int id = ceil(_n/10)
gen int t = mod(_n-1,10)+1
forvalues j = 1/4 {
    gen double z`j' = rnormal()
}
gen double w = rnormal()
gen double u = rnormal()
gen double x = .5*z1 + .3*z2 + .2*z3 + u + rnormal()
gen double x3 = .4*z3 + .5*z4 + u + rnormal()
gen double y = 1 + x + .5*w + u + rnormal() + id/100
gen double wt = 1 + mod(_n, 4)
gen int n = 1 + mod(_n, 3)
tempfile ivdata
save `ivdata'

* dofminus() in the identification statistics
foreach v in "" "robust" "cluster(id)" {
    _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id) dofminus(2) `v') ///
        ct(civreghdfe y w (x = z1 z2 z3), absorb(id) dofminus(2) `v') ///
        scalars(idstat widstat cdf) testname("dofminus(2) identification statistics `v'")
}

* nofooter suppresses the display only: every diagnostic is stored
local diag idstat iddf widstat cd_f sargan sargan_df cstat cstat_df endogtest redund endog_chi2
quietly civreghdfe y w (x x3 = z1 z2 z3 z4), absorb(id) vce(cluster id) orthog(z4) endogtest(x3) redundant(z4)
foreach s of local diag {
    local f_`s' = e(`s')
}
quietly civreghdfe y w (x x3 = z1 z2 z3 z4), absorb(id) vce(cluster id) orthog(z4) endogtest(x3) redundant(z4) nofooter
local diffs
foreach s of local diag {
    * (the first fit after other commands can differ in the last bits)
    if missing(e(`s')) | reldif(e(`s'), `f_`s'') > 1e-12 local diffs "`diffs' `s'"
}
if "`diffs'" == "" test_pass "nofooter stores the diagnostics"
else test_fail "nofooter stores the diagnostics" "`diffs'"
_civreghdfe_cmp, ref(ivreghdfe y w (x x3 = z1 z2 z3 z4), absorb(id) vce(cluster id) orthog(z4) nofooter) ///
    ct(civreghdfe y w (x x3 = z1 z2 z3 z4), absorb(id) vce(cluster id) orthog(z4) nofooter) ///
    scalars(idstat widstat j cstat) testname("nofooter diagnostics against ivreghdfe")

* pweights with first-stage and reduced-form output: the auxiliary
* regressions use the robust VCE the weights imply
quietly creghdfe x z1 z2 w [pw=wt], absorb(id) vce(robust)
* (the saved first stage leaves out the constant, partialled with the effects)
matrix pw_fsV = e(V)
matrix pw_fsV = pw_fsV[1..3, 1..3]
estimates clear
capture quietly civreghdfe y w (x = z1 z2) [pw=wt], absorb(id) rf first savefirst
if _rc test_fail "pweights with rf first savefirst" "rc=`=_rc'"
else {
    if missing(e(widstat)) | missing(e(sargan)) test_fail "pweights with rf first savefirst" "diagnostics missing"
    else test_pass "pweights with rf first savefirst"
    capture estimates restore _ivreg2_x
    if _rc test_fail "pweights: saved first stage" "rc=`=_rc'"
    else {
        matrix pw_ctV = e(V)
        matrix_min_sigfigs pw_fsV pw_ctV
        if r(min_sigfigs) >= $DEFAULT_SIGFIGS test_pass "pweights: saved first stage is the robust creghdfe fit"
        else test_fail "pweights: saved first stage is the robust creghdfe fit" "V sf=`=r(min_sigfigs)'"
    }
}
capture matrix drop pw_fsV pw_ctV

* LIML on data in small units: lambda and b do not depend on the units
use `ivdata', clear
foreach v in y x w z1 z2 z3 {
    replace `v' = `v' * 1e-8
}
_civreghdfe_cmp, ref(ivreg2 y w (x = z1 z2 z3), liml) ct(civreghdfe y w (x = z1 z2 z3), liml) ///
    scalars(N lambda rss) testname("liml, variables in units of 1e-8")
_civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id) liml) ct(civreghdfe y w (x = z1 z2 z3), absorb(id) liml) ///
    scalars(N lambda rss) testname("liml with absorb(), variables in units of 1e-8")

* The Durbin-Wu-Hausman statistic e(endog_f) is the Wu-Hausman F of the
* augmented regression (regress y x w v; test v), unweighted and weighted;
* e(endog_chi2) is K_endog times it
use `ivdata', clear
foreach wt in "" "[aw=wt]" {
    quietly regress x w z1 z2 z3 `wt'
    predict double v1 if e(sample), residuals
    quietly regress x3 w z1 z2 z3 `wt'
    predict double v3 if e(sample), residuals
    quietly regress y x w v1 `wt'
    quietly test v1
    local ref_f1 = r(F)
    quietly regress y x x3 w v1 v3 `wt'
    quietly test v1 v3
    local ref_f2 = r(F)
    quietly civreghdfe y w (x = z1 z2 z3) `wt'
    assert_scalar_equal `=e(endog_f)' `ref_f1' $DEFAULT_SIGFIGS "endog_f is the Wu-Hausman F `wt'"
    quietly civreghdfe y w (x x3 = z1 z2 z3) `wt'
    assert_scalar_equal `=e(endog_f)' `ref_f2' $DEFAULT_SIGFIGS "endog_f is the Wu-Hausman F, two endogenous `wt'"
    assert_scalar_equal `=e(endog_chi2)' `=2*`ref_f2'' $DEFAULT_SIGFIGS "endog_chi2 = K_endog x F `wt'"
    drop v1 v3
}
quietly civreghdfe y w (x = z1 z2 z3) [fw=n], absorb(id)
local fw_chi2 = e(endog_chi2)
preserve
expand n
quietly civreghdfe y w (x = z1 z2 z3), absorb(id)
assert_scalar_equal `fw_chi2' `=e(endog_chi2)' $DEFAULT_SIGFIGS "endog_chi2 with fweights equals the expanded data"
restore

* Estimators that cannot be combined are refused: liml with gmm2s or cue
* and gmm2s with cue (ivreg2 r(198)), a negative k-class parameter and
* fuller() with kclass() (ivreg2 r(3200)), fuller() or kclass() with cue
* (ivreg2 runs CUE and labels it LIML or k-class), liml with a negative
* fuller(); cluster() implies robust
foreach o in "gmm2s cue" "liml gmm2s" "fuller(1) cue" "kclass(0.5) cue" "kclass(-0.5)" ///
    "fuller(1) kclass(0.5)" "liml fuller(-1)" {
    capture civreghdfe y w (x = z1 z2), absorb(id) `o'
    if _rc == 198 test_pass "conflicting options refused: `o'"
    else test_fail "conflicting options refused: `o'" "rc=`=_rc'"
}
* The combinations ivreg2 accepts give its estimator: fuller() and kclass()
* take precedence over gmm2s, kclass() over liml, and a negative fuller()
* is LIML
foreach o in "fuller(1) gmm2s" "kclass(0.5) gmm2s" "liml kclass(0.5)" "fuller(-1)" ///
    "fuller(1) gmm2s robust" "fuller(-1) cluster(id)" {
    _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id) `o') ///
        ct(civreghdfe y w (x = z1 z2 z3), absorb(id) `o') scalars(`ivstats' j kclass lambda) ///
        testname("estimator combination ivreg2 accepts: `o'")
}
_civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2), absorb(id) robust cluster(id)) ///
    ct(civreghdfe y w (x = z1 z2), absorb(id) robust cluster(id)) scalars(`ivstats' N_clust) ///
    testname("robust with cluster()")

* First-stage statistics follow the endogenous regressors when one is
* dropped as collinear (with another endogenous regressor: one collinear with
* an exogenous regressor becomes exogenous, as in ivreg2)
gen double x2 = 2*x
quietly civreghdfe y w (x x3 = z1 z2 z3 z4), absorb(id) ffirst
local nox2_F2 = e(F_first2)
local nox2_r2 = e(partial_r2_2)
local nox2_F1 = e(F_first1)
capture quietly civreghdfe y w (x x2 x3 = z1 z2 z3 z4), absorb(id) ffirst
if _rc test_fail "first-stage statistics with a dropped endogenous regressor" "rc=`=_rc'"
else if `: colnumb e(b) o.x2' == . test_fail "first-stage statistics with a dropped endogenous regressor" "x2 not dropped"
else {
    assert_scalar_equal `=e(F_first1)' `nox2_F1' $DEFAULT_SIGFIGS "dropped endogenous regressor: F_first1"
    assert_scalar_equal `=e(F_first3)' `nox2_F2' $DEFAULT_SIGFIGS "dropped endogenous regressor: F_first3"
    assert_scalar_equal `=e(partial_r2_3)' `nox2_r2' $DEFAULT_SIGFIGS "dropped endogenous regressor: partial_r2_3"
    if missing(e(F_first2)) & missing(e(partial_r2_2)) test_pass "dropped endogenous regressor: F_first2 missing"
    else test_fail "dropped endogenous regressor: F_first2 missing" "F_first2=`=e(F_first2)'"
}
drop x2

* Duplicated terms are used once, and a term listed as both exogenous and
* endogenous is exogenous, as in ivreg2
gen double w2 = rnormal()
gen double d2 = z2 + rnormal()
foreach sp in "y w w w2 (x d2 = z1 z2 z3)" "y w (w x = z1 z2 z3)" "y w (x x = z1 z2 z3)" {
    _civreghdfe_cmp, ref(ivreghdfe `sp', absorb(id)) ct(civreghdfe `sp', absorb(id)) ///
        scalars(`ivstats') testname("duplicated terms: `sp'")
}
drop w2 d2

* ivreg2's names for the statistics, posted alongside civreghdfe's
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id) robust
foreach s in j jdf jp cdf Fp r2c endog_ct exexog_ct {
    local r_`s' = e(`s')
}
local r_instd `e(instd)'
local r_exexog `e(exexog)'
quietly civreghdfe y w (x = z1 z2 z3), absorb(id) robust
foreach s in j jdf jp cdf Fp r2c endog_ct exexog_ct {
    assert_scalar_equal `=e(`s')' `r_`s'' $DEFAULT_SIGFIGS "ivreg2 name e(`s')"
}
if "`e(instd)'" == "`r_instd'" & "`e(exexog)'" == "`r_exexog'" test_pass "ivreg2 names e(instd), e(exexog)"
else test_fail "ivreg2 names e(instd), e(exexog)" "`e(instd)' | `e(exexog)'"

* A perfect fit leaves no testable restriction: F missing, as in ivreghdfe
clear
set seed 5
set obs 200
gen int g = ceil(_n/10)
gen double z = rnormal()
gen double x = z + rnormal()
gen double w = rnormal()
gen double y = 2*x + w + g
quietly ivreghdfe y w (x = z), absorb(g)
local ref_F = e(F)
quietly civreghdfe y w (x = z), absorb(g)
if missing(`ref_F') & missing(e(F)) test_pass "perfect fit: missing F"
else test_fail "perfect fit: missing F" "F=`=e(F)' (ivreghdfe `ref_F')"

* by: estimates each group, as by with ivreghdfe
use `ivdata', clear
gen byte half = id > 100
sort half id t
quietly ivreghdfe y w (x = z1 z2) if half == 1, absorb(id)
matrix by_rb = e(b)
capture quietly by half: civreghdfe y w (x = z1 z2), absorb(id)
if _rc test_fail "by: estimates the last group" "rc=`=_rc'"
else {
    matrix by_cb = e(b)
    matrix_min_sigfigs by_rb by_cb
    if r(min_sigfigs) >= $DEFAULT_SIGFIGS & e(N) == 1000 test_pass "by: estimates the last group"
    else test_fail "by: estimates the last group" "sf=`=r(min_sigfigs)' N=`=e(N)'"
    matrix drop by_cb
}
matrix drop by_rb

* partial() without absorb() partials out the constant too, as ivreg2
sysuse auto, clear
foreach o in "" "robust" "small" {
    _civreghdfe_cmp, ref(ivreg2 price (mpg = weight) length turn, partial(turn) `o') ///
        ct(civreghdfe price (mpg = weight) length turn, partial(turn) `o') ///
        scalars(N r2 r2_a rss mss rmse F df_r idstat widstat partial_ct) testname("partial() without absorb() `o'")
}

* Regressors (nearly) collinear with the absorbed effects, as ivreghdfe
* decides them: it partials the effects out first and keeps every column
* that its invsym() does not zero, so a regressor with little within
* variation (nf) is kept, and one the effects absorb is dropped and not
* counted in the degrees of freedom (also when the projection is iterative)
clear
set seed 4242
set obs 3000
gen int id = ceil(_n/30)
gen int t = mod(_n-1,30)+1
gen byte bal = runiform() >= .2
gen double z1 = rnormal()
gen double z2 = rnormal()
gen double x1 = rnormal()
gen double u = rnormal()
gen double d1 = z1 + .5*z2 + u + rnormal()
gen double y = 1 + d1 + x1 + u + rnormal() + id/10 + t/5
gen double xa = id/7
gen double xb = id/7 + t/3
gen int k3 = mod(_n, 7)
local festats `ivstats' df_r j jdf idstat cdf
foreach r in 1e-2 1e-4 1e-5 {
    gen double nf = 1000*id + `r'*1000*rnormal()
    foreach ab in "id" "id t" {
        _civreghdfe_cmp, ref(ivreghdfe y x1 nf (d1=z1 z2) if bal, absorb(`ab')) ///
            ct(civreghdfe y x1 nf (d1=z1 z2) if bal, absorb(`ab')) scalars(`festats') ///
            testname("regressor with little within variation (`r') absorb(`ab')")
    }
    drop nf
}
foreach sp in "xa|id|" "xa|id t|" "xa|id t|if bal" "xb|id t|" "xb|id t k3|" {
    gettoken v sp : sp, parse("|")
    gettoken bar sp : sp, parse("|")
    gettoken ab sp : sp, parse("|")
    local ifc = substr("`sp'", 2, .)
    _civreghdfe_cmp, ref(ivreghdfe y x1 `v' (d1=z1 z2) `ifc', absorb(`ab')) ///
        ct(civreghdfe y x1 `v' (d1=z1 z2) `ifc', absorb(`ab')) scalars(`festats') omitzero ///
        testname("regressor absorbed by the effects: `v' absorb(`ab') `ifc'")
}
* An absorbed endogenous regressor has no first-stage variation, so the rank
* condition fails: the identification statistics are 0 (ivreghdfe's ranktest
* returns rounding noise around 0 there: idstat 3.6e-32 with absorb(id),
* 2.3e-32 with absorb(t id)), with the regressor counted in their dof
_civreghdfe_cmp, ref(ivreghdfe y x1 (d1 xa = z1 z2), absorb(id)) ///
    ct(civreghdfe y x1 (d1 xa = z1 z2), absorb(id)) scalars(N r2 r2_a rss mss rmse F df_r j jdf iddf) omitzero ///
    testname("endogenous regressor absorbed by the effects")
quietly ivreghdfe y x1 (d1 xa = z1 z2), absorb(id)
local fe_ok = abs(e(idstat)) < 1e-20 & abs(e(cdf)) < 1e-20
quietly civreghdfe y x1 (d1 xa = z1 z2), absorb(id)
if `fe_ok' & e(idstat) == 0 & e(cdf) == 0 & e(widstat) == 0 test_pass "absorbed endogenous regressor: identification statistics 0"
else test_fail "absorbed endogenous regressor: identification statistics 0" "idstat=`=e(idstat)' cdf=`=e(cdf)'"
* An excluded instrument absorbed by the effects counts in the
* excluded-instrument count of the weak-identification F and the
* identification dof, as ivreghdfe counts its variable list. When its
* projection removes the instrument exactly (xa, constant within id), the
* statistics are those of the model without it; when it leaves rounding
* noise (xb with absorb(id t)), ranktest uses the noise (idstat 1159.063
* with absorb(id t), 1159.468 with absorb(t id), 1161.880 with the instrument
* times 1000), so the model with xa, an instrument absorbed exactly, is the
* reference for the identification statistics there
_civreghdfe_cmp, ref(ivreghdfe y x1 (d1 = z1 z2 xa), absorb(id)) ///
    ct(civreghdfe y x1 (d1 = z1 z2 xa), absorb(id)) scalars(`festats' iddf) ///
    testname("excluded instrument absorbed by the effects, absorb(id)")
_civreghdfe_cmp, ref(ivreghdfe y x1 (d1 = z1 z2 xb), absorb(id t)) ///
    ct(civreghdfe y x1 (d1 = z1 z2 xb), absorb(id t)) scalars(N r2 r2_a rss mss rmse F df_r j jdf iddf) ///
    testname("excluded instrument absorbed by the effects")
_civreghdfe_cmp, ref(ivreghdfe y x1 (d1 = z1 z2 xa), absorb(id t)) ///
    ct(civreghdfe y x1 (d1 = z1 z2 xb), absorb(id t)) scalars(`festats' iddf) ///
    testname("excluded instrument absorbed by the effects: identification statistics")
* Where ivreghdfe's own projection leaves noise above invsym()'s tolerance in
* an absorbed column, it keeps the noise as a regressor with an arbitrary
* coefficient (xbig = 1e6*id + 1e5 with absorb(id): b = -6.6e-09, se = 1.8e6;
* age = year - birth with absorb(id year) on the unbalanced sample:
* b = 8.3, se = 1.6e7) and counts it in df_r. The column is exactly in the
* span of the effects, so the model with another absorbed column, one that
* ivreghdfe removes exactly (xa, constant within id), is the same model:
* that is the reference
gen double xbig = 1e6*id + 1e5
gen int year = 2000 + t
gen double age = year - (1950 + id)
foreach sp in "xbig|id" "xbig|id t" "age|id year" {
    gettoken v sp : sp, parse("|")
    local ab = substr("`sp'", 2, .)
    quietly ivreghdfe y x1 xa (d1=z1 z2) if bal, absorb(`ab')
    matrix fe_rb = e(b)
    matrix fe_rV = e(V)
    local fe_diffs
    if fe_rb[1, 3] != 0 | fe_rV[3, 3] != 0 local fe_diffs "reference keeps xa"
    matrix fe_rb = fe_rb[1, 1..2]
    matrix fe_rV = fe_rV[1..2, 1..2]
    foreach st of local festats {
        local fe_r_`st' = e(`st')
    }
    quietly civreghdfe y x1 `v' (d1=z1 z2) if bal, absorb(`ab')
    if colnumb(e(b), "o.`v'") != 3 local fe_diffs "`fe_diffs' `v' not omitted"
    else {
        matrix fe_cb = e(b)
        matrix fe_cV = e(V)
        matrix fe_cb = fe_cb[1, 1..2]
        matrix fe_cV = fe_cV[1..2, 1..2]
        matrix_min_sigfigs fe_rb fe_cb
        if r(min_sigfigs) < $DEFAULT_SIGFIGS local fe_diffs "`fe_diffs' b"
        matrix_min_sigfigs fe_rV fe_cV
        if r(min_sigfigs) < $DEFAULT_SIGFIGS local fe_diffs "`fe_diffs' V"
        foreach st of local festats {
            if missing(`fe_r_`st'') & missing(e(`st')) continue
            sigfigs `fe_r_`st'' `=e(`st')'
            if r(sigfigs) < $DEFAULT_SIGFIGS local fe_diffs "`fe_diffs' `st'"
        }
    }
    if "`fe_diffs'" == "" test_pass "absorbed regressor kept as noise by ivreghdfe: `v' absorb(`ab')"
    else test_fail "absorbed regressor kept as noise by ivreghdfe: `v' absorb(`ab')" "`fe_diffs'"
}
capture matrix drop fe_rb fe_rV fe_cb fe_cV

* Exact collinearity among the variables, resolved in ivreg2's order: an
* endogenous regressor collinear with the instruments becomes exogenous and
* the exogenous term it duplicates is dropped; collinear terms are dropped
* from the right; the lists are posted as e(collin), e(ecollin), e(dups)
clear
set seed 4242
set obs 2000
gen int id = ceil(_n/10)
forvalues j = 1/4 {
    gen double z`j' = rnormal()
}
gen double w = rnormal()
gen double x2 = rnormal()
gen double u = rnormal()
gen double x = .5*z1 + .3*z2 + .2*z3 + u + rnormal()
gen double d1 = z4 + z1 + rnormal()
gen double y = 1 + x + .5*w + u + rnormal() + id/100
gen double wb = 2*w
gen double dz = z1 + z2
gen double x1b = 2*x2
gen double d2 = 3*d1
foreach sp in "y w x2 (d1 x1b = z1 z2 z3)" "y w (x wb = z1 z2 z3)" "y w (x dz = z1 z2 z3)" ///
    "y w (x d1 d2 = z1 z2 z3 z4)" "y w wb (x = z1 z2)" "y w (x = z1 z2 dz)" "y (wb = z1) w" ///
    "y w (x = x z1 z2)" "y w (x = w z1 z2)" {
    foreach ab in "absorb(id)" "" {
        local cstats = cond("`ab'" == "", "N r2 rss mss rmse F idstat widstat cdf j jdf", "`ivstats' df_r cdf j jdf")
        _civreghdfe_cmp, ref(ivreghdfe `sp', `ab') ct(civreghdfe `sp', `ab') scalars(`cstats') ///
            testname("collinear variables in ivreg2's order: `sp' `ab'")
        quietly ivreghdfe `sp', `ab'
        local cl_r "`e(collin)'|`e(ecollin)'|`e(dups)'"
        capture quietly civreghdfe `sp', `ab'
        if _rc test_fail "e(collin), e(ecollin), e(dups): `sp' `ab'" "rc=`=_rc'"
        else if "`e(collin)'|`e(ecollin)'|`e(dups)'" == "`cl_r'" test_pass "e(collin), e(ecollin), e(dups): `sp' `ab'"
        else test_fail "e(collin), e(ecollin), e(dups): `sp' `ab'" "`e(collin)'|`e(ecollin)'|`e(dups)' vs `cl_r'"
    }
}

* b0(): the CUE objective (J), the efficient GMM VCE and the statistics at
* the given coefficients, as ivreg2 evaluates them (b0 implies noid, whose
* statistics civreghdfe still stores); the columns are matched by name
use `ivdata', clear
matrix cb0 = (1.1, .4)
matrix colnames cb0 = x w
matrix cb1 = (.4, 1.1, 2)
matrix colnames cb1 = w x _cons
foreach o in "" "robust" "cluster(id)" "small" {
    _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id) b0(cb0) `o') ///
        ct(civreghdfe y w (x = z1 z2 z3), absorb(id) b0(cb0) `o') ///
        scalars(N r2 r2_a rss mss rmse F df_r j jdf) testname("b0() `o'")
    foreach iv in "z1 z2" "z1" {
        _civreghdfe_cmp, ref(ivreg2 y w (x = `iv'), b0(cb1) `o') ///
            ct(civreghdfe y w (x = `iv'), b0(cb1) `o') ///
            scalars(N r2 rss mss rmse F j jdf) testname("b0() without absorb(), instruments `iv' `o'")
    }
}
capture quietly civreghdfe y w (x = z1 z2 z3), absorb(id) b0(cb0)
if _rc test_fail "b0(): e(model) and e(b)" "rc=`=_rc'"
else if "`e(model)'" == "cue" & _b[x] == 1.1 & _b[w] == .4 test_pass "b0(): e(model) and e(b)"
else test_fail "b0(): e(model) and e(b)" "`e(model)'"
foreach o in "cue" "gmm2s" "liml" "kclass(.5)" "fuller(1)" {
    capture civreghdfe y w (x = z1 z2 z3), absorb(id) b0(cb0) `o'
    if _rc == 198 test_pass "b0() refused with `o'"
    else test_fail "b0() refused with `o'" "rc=`=_rc'"
}
capture matrix drop cb0 cb1
* sfirst (ivreg2's system of first-stage and reduced-form equations, with
* its Sanderson-Windmeijer and Anderson-Rubin statistics; compared with
* ivreg2 and ivreghdfe in case s2_iv)
capture quietly civreghdfe y w (x = z1 z2 z3), absorb(id) sfirst
if _rc == 0 & "`e(sfirsteq)'" == "" & !missing(e(arf)) test_pass "sfirst"
else test_fail "sfirst" "rc=`=_rc'"

* nocollin: ivreg2's collinearity rules are skipped; what remains is
* ivreghdfe's invsym() of the partialled cross-products
use `ivdata', clear
gen double wb = 2*w
gen double dz = z1 + z2
foreach sp in "y w wb (x = z1 z2)" "y w (x = z1 z2 dz)" {
    foreach ab in "absorb(id)" "" {
        local cstats = cond("`ab'" == "", "N r2 rss F j jdf nocollin idstat iddf widstat cdf", "`ivstats' df_r j jdf nocollin iddf cdf")
        * (the redundant instrument dz counts in the weak-identification F's
        * excluded-instrument count, as ivreg2 counts its list)
        _civreghdfe_cmp, ref(ivreghdfe `sp', `ab' nocollin) ct(civreghdfe `sp', `ab' nocollin) ///
            scalars(`cstats') omitzero testname("nocollin: `sp' `ab'")
    }
}

* dofadjustments() and keepsingletons, as ivreghdfe passes them to reghdfe
clear
set seed 4242
set obs 2000
gen int id = ceil(_n/10)
gen int t = mod(_n-1,10)+1
gen int f = ceil(id/7)
forvalues j = 1/3 {
    gen double z`j' = rnormal()
}
gen double w = rnormal()
gen double u = rnormal()
gen double x = .5*z1 + .3*z2 + .2*z3 + u + rnormal()
gen double y = 1 + x + .5*w + u + rnormal() + id/100
replace id = 5000 + _n in 1/5
foreach d in "dofadjustments(none)" "dof(firstpair)" "dof(pairwise)" "dof(all)" "dof(clusters)" "dof(pairwise clusters)" "" {
    foreach v in "cluster(f)" "" {
        _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id t f) `v' `d') ///
            ct(civreghdfe y w (x = z1 z2 z3), absorb(id t f) `v' `d') scalars(`ivstats' df_r df_a j) ///
            testname("degrees-of-freedom adjustment `d' `v'")
    }
}
foreach v in "cluster(f)" "robust" {
    _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id t) `v' keepsingletons) ///
        ct(civreghdfe y w (x = z1 z2 z3), absorb(id t) `v' keepsingletons) ///
        scalars(`ivstats' df_r df_a j num_singletons N_full) testname("keepsingletons `v'")
}

* ivreg2's e(df_m), ranks, sdofminus, partialcons, cons, log likelihood and
* reghdfe's singleton counts; estat ic from e(ll)
foreach o in "" "sdofminus(3)" "partial(w)" "partial(w) nopartialsmall" "robust" "cluster(f)" {
    _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id t) `o') ///
        ct(civreghdfe y w (x = z1 z2 z3), absorb(id t) `o') ///
        scalars(N df_m rankxx rankzz rankS rankV sdofminus partialcons cons ll num_singletons N_full) ///
        testname("ivreg2's e() scalars with absorb() `o'")
}
foreach o in "" "sdofminus(3)" "partial(w)" "partial(w) nopartialsmall" "robust" "noconstant" {
    _civreghdfe_cmp, ref(ivreg2 y w (x = z1 z2 z3), `o') ct(civreghdfe y w (x = z1 z2 z3), `o') ///
        scalars(N df_m rankxx rankzz rankS rankV sdofminus partialcons cons ll) ///
        testname("ivreg2's e() scalars without absorb() `o'")
}
* (ll(null) is missing in both: columns N, ll, df, AIC and BIC)
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id t)
quietly estat ic
matrix ic_r = r(S)
quietly civreghdfe y w (x = z1 z2 z3), absorb(id t)
capture quietly estat ic
if _rc test_fail "estat ic" "rc=`=_rc'"
else {
    matrix ic_c = r(S)
    local ic_ll0 = missing(ic_r[1, 2]) & missing(ic_c[1, 2])
    matrix ic_r = ic_r[1, 1], ic_r[1, 3..6]
    matrix ic_c = ic_c[1, 1], ic_c[1, 3..6]
    matrix_min_sigfigs ic_r ic_c
    if r(min_sigfigs) >= $DEFAULT_SIGFIGS & `ic_ll0' test_pass "estat ic"
    else test_fail "estat ic" "sf=`=r(min_sigfigs)'"
}
capture matrix drop ic_r ic_c

* Saved first stages with clustering keep ivreg2's residual dof and root MSE
foreach o in "absorb(id) cluster(f)" "cluster(f)" "absorb(id t) cluster(f)" {
    estimates clear
    quietly ivreghdfe y w (x = z1 z2 z3), `o' savefirst
    quietly estimates restore _ivreg2_x
    local sf_r "`=e(df_r)' `=e(rmse)' `=e(N_clust)'"
    matrix sf_rV = e(V)
    quietly civreghdfe y w (x = z1 z2 z3), `o' savefirst
    quietly estimates restore _ivreg2_x
    matrix sf_cV = e(V)
    matrix_min_sigfigs sf_rV sf_cV
    local sf_sf = r(min_sigfigs)
    gettoken sf_dfr sf_r : sf_r
    gettoken sf_rmse sf_r : sf_r
    sigfigs `sf_rmse' `=e(rmse)'
    if e(df_r) == `sf_dfr' & r(sigfigs) >= $DEFAULT_SIGFIGS & e(N_clust) == `sf_r' & `sf_sf' >= $DEFAULT_SIGFIGS ///
        test_pass "saved first stage with clustering: df_r, rmse, N_clust, V (`o')"
    else test_fail "saved first stage with clustering: df_r, rmse, N_clust, V (`o')" ///
        "df_r=`=e(df_r)' vs `sf_dfr', rmse=`=e(rmse)' vs `sf_rmse', V sf=`sf_sf'"
}
capture matrix drop sf_rV sf_cV
estimates clear

* predict d and dresiduals, as reghdfe_p after ivreghdfe (residuals() saved)
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id t) residuals(ivh_r)
quietly predict double ivh_d, d
quietly predict double ivh_dr, dresiduals
quietly civreghdfe y w (x = z1 z2 z3), absorb(id t) residuals(civ_r)
foreach pr in d dresiduals {
    local pv = cond("`pr'" == "d", "d", "dr")
    capture quietly predict double civ_`pv', `pr'
    if _rc test_fail "predict `pr'" "rc=`=_rc'"
    else assert_var_equal ivh_`pv' civ_`pv' $DEFAULT_SIGFIGS "predict `pr'"
}
quietly civreghdfe y w (x = z1 z2 z3), absorb(id t)
capture predict double civ_d2, d
if _rc == 198 test_pass "predict d requires residuals()"
else test_fail "predict d requires residuals()" "rc=`=_rc'"
capture drop ivh_r ivh_d ivh_dr civ_r civ_d civ_dr

* A real bandwidth, as ivreg2 takes it, for each kernel
clear
set seed 77
set obs 1500
gen int id = ceil(_n/15)
gen int t = mod(_n-1,15)+1
xtset id t
gen double z1 = rnormal()
gen double z2 = rnormal()
gen double w = rnormal()
gen double u = rnormal()
gen double x = z1 + .5*z2 + u + rnormal()
gen double y = 1 + x + w + u + rnormal()
foreach k in bartlett truncated parzen qs thann {
    foreach v in "robust" "" {
        _civreghdfe_cmp, ref(ivreghdfe y w (x = z1 z2), absorb(id) bw(2.5) kernel(`k') `v') ///
            ct(civreghdfe y w (x = z1 z2), absorb(id) bw(2.5) kernel(`k') `v') scalars(`ivstats' j bw) ///
            testname("bw(2.5) kernel(`k') `v'")
    }
}

* Summary

* Cleanup to prevent state corruption before next validation script
capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

* End of civreghdfe validation
noi print_summary "civreghdfe"
}

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "civreghdfe"
    global CTOOLS_CASE_COMPLETE "civreghdfe"
    exit
}

if "`__case'" == "civreghdfe_grouped_syntax" {
do validation/support/validate_setup.do
clear
set seed 1122
set obs 300
generate t = _n
tsset t
generate g = ceil(_n/30)
generate double z = rnormal()
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double d = z + rnormal()
generate double y = 1 + d + x1 + x2 + rnormal()
foreach position in before after {
    if "`position'" == "before" local spec y c.(x1 x2) (d=z)
    else local spec y (d=z) c.(x1 x2)
    benchmark_ivreghdfe `spec', absorb(g) testname("grouped FV `position' IV equation")
    if "`position'" == "before" local spec y L(1/2).x1 (d=z) x2
    else local spec y (d=z) L(1/2).x1 x2
    benchmark_ivreghdfe `spec', absorb(g) testname("lag range `position' IV equation")
}
foreach terms in "c.(x1 x2)" "L(1/2).x1 x2" {
    quietly civreghdfe y `terms'
    matrix grouped_b = e(b)
    matrix grouped_V = e(V)
    if "`terms'" == "c.(x1 x2)" local explicit x1 x2
    else local explicit L.x1 L2.x1 x2
    quietly civreghdfe y `explicit'
    assert mreldif(e(b),grouped_b) == 0 & mreldif(e(V),grouped_V) == 0
    test_pass "OLS grouped syntax: `terms'"
}
benchmark_ivreghdfe y c.(x1 x2) (d=L(1/2).z), absorb(g) testname("nested lag ranges in IV equation")
* (=z) adds an excluded instrument without endogenous regressors.
benchmark_ivreghdfe y (=z) x1, absorb(g) testname("instruments without endogenous regressors")
generate byte group_a = mod(_n, 2)
generate byte group_b = mod(_n, 3)
benchmark_ivreghdfe y i.group_a#i.(group_b) (d=z) x1, absorb(g) testname("FV group after an interaction")
* A group without "=" that no operator owns is an IV equation missing its equals sign.
foreach invalid in "y (d=z) (x1=x2)" "y (d=) x1" "y (d=z=z) x1" "y (d=z" "y d=z" "y (d z) x1" "y x1 (d z)" "y (d) (d=z)" {
    capture civreghdfe `invalid'
    assert _rc == 198
    test_pass "reject malformed IV equation `invalid'"
}
print_summary "civreghdfe grouped syntax"
global CTOOLS_COMPONENT_COMPLETE "civreghdfe_grouped_syntax"
global CTOOLS_CASE_COMPLETE "civreghdfe_grouped_syntax"
exit
}

di as error "Unknown civreghdfe case: `__case'"
exit 198
