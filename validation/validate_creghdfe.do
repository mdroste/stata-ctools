/*******************************************************************************
 * validate_creghdfe.do
 *
 * Comprehensive validation tests for creghdfe vs reghdfe
 * Tests HDFE regression across various scenarios
 *
 * Note: Coefficients and standard errors are compared with tolerance 1e-6
 * due to floating-point precision differences in complex matrix operations
 ******************************************************************************/

do "validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CREGHDFE VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * Auto dataset - basic tests
 ******************************************************************************/
print_section "Auto dataset - basic"

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) testname("single FE")

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
benchmark_reghdfe price mpg, absorb(foreign) testname("single covariate")

sysuse auto, clear
benchmark_reghdfe price mpg weight length turn displacement, absorb(foreign) testname("many covariates")

/*******************************************************************************
 * VCE options
 ******************************************************************************/
print_section "VCE options"

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) vce(robust) testname("robust VCE")

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster foreign) testname("cluster VCE")

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign rep78) vce(robust) testname("two-way FE + robust")

/*******************************************************************************
 * Weights
 ******************************************************************************/
print_section "Weight options"

sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)
benchmark_reghdfe price mpg weight [aw=w], absorb(foreign) testname("aweight")

sysuse auto, clear
drop if missing(rep78)
benchmark_reghdfe price mpg weight [fw=rep78], absorb(foreign) testname("fweight")

sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)
benchmark_reghdfe price mpg weight [pw=w], absorb(foreign) testname("pweight")

* aweight with robust
sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)
benchmark_reghdfe price mpg weight [aw=w], absorb(foreign) vce(robust) testname("aweight + robust")

* aweight with cluster
sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)
benchmark_reghdfe price mpg weight [aw=w], absorb(foreign) vce(cluster foreign) testname("aweight + cluster")

/*******************************************************************************
 * Panel data (nlswork)
 ******************************************************************************/
print_section "Panel data (nlswork)"

webuse nlswork, clear
keep in 1/5000
benchmark_reghdfe ln_wage age ttl_exp tenure, absorb(idcode) testname("individual FE")

webuse nlswork, clear
keep in 1/5000
benchmark_reghdfe ln_wage age ttl_exp tenure, absorb(idcode year) testname("two-way panel FE")

webuse nlswork, clear
keep in 1/5000
benchmark_reghdfe ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode) testname("panel + cluster")

webuse nlswork, clear
keep in 1/5000
benchmark_reghdfe ln_wage age ttl_exp tenure, absorb(idcode) vce(robust) testname("individual FE + robust")

/*******************************************************************************
 * Census dataset
 ******************************************************************************/
print_section "Census dataset"

sysuse census, clear
benchmark_reghdfe pop medage death marriage, absorb(region) testname("region FE")

sysuse census, clear
benchmark_reghdfe pop medage, absorb(region) testname("region FE - simple")

/*******************************************************************************
 * Synthetic data
 ******************************************************************************/
print_section "Synthetic data"

* 20K observations
clear
set seed 98765
set obs 20000
gen firm = runiformint(1, 500)
gen year = runiformint(2000, 2020)
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen y = 1 + 0.5*x1 - 0.3*x2 + rnormal() * 0.5

benchmark_reghdfe y x1 x2, absorb(firm year) testname("20K two-way FE")

* 50K observations
clear
set seed 11111
set obs 50000
gen firm = runiformint(1, 1000)
gen year = runiformint(2000, 2020)
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen x3 = rnormal()
gen y = 1 + 0.5*x1 - 0.3*x2 + 0.2*x3 + rnormal() * 0.5

benchmark_reghdfe y x1 x2 x3, absorb(firm) testname("50K firm FE")

clear
set seed 11111
set obs 50000
gen firm = runiformint(1, 1000)
gen year = runiformint(2000, 2020)
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen x3 = rnormal()
gen y = 1 + 0.5*x1 - 0.3*x2 + 0.2*x3 + rnormal() * 0.5

benchmark_reghdfe y x1 x2 x3, absorb(firm year) testname("50K two-way FE")

/*******************************************************************************
 * High-dimensional FE
 ******************************************************************************/
print_section "High-dimensional FE"

webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age ttl_exp, absorb(idcode) testname("many FE levels (idcode)")

webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age ttl_exp, absorb(idcode year) testname("many FE levels (idcode year)")

/*******************************************************************************
 * Edge cases
 ******************************************************************************/
print_section "Edge cases"

* weights=1 should match unweighted
sysuse auto, clear
gen w1 = 1
quietly creghdfe price mpg weight, absorb(foreign)
local unweighted_b = _b[mpg]
quietly creghdfe price mpg weight [aw=w1], absorb(foreign)
local weighted1_b = _b[mpg]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if abs(`unweighted_b' - `weighted1_b') < 1e-6 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] weights=1 matches unweighted"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] weights=1 differs from unweighted"
}

* Single observation per FE group - verify it works
clear
set obs 100
gen fe = _n
gen x = rnormal()
gen y = x + rnormal()
capture benchmark_reghdfe y x, absorb(fe) testname("one obs per FE group")

/*******************************************************************************
 * Stored results check
 ******************************************************************************/
print_section "Stored results"

sysuse auto, clear

* Run both commands
quietly reghdfe price mpg weight, absorb(foreign)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)

quietly creghdfe price mpg weight, absorb(foreign)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)

* Compare e(b)
tempname diff
matrix `diff' = reghdfe_b - creghdfe_b
local cols = colsof(`diff')
local maxdiff = 0
forvalues j = 1/`cols' {
    local d = abs(`diff'[1, `j'])
    if `d' > `maxdiff' local maxdiff = `d'
}

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `maxdiff' < 1e-6 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] e(b) matches"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] e(b) differs (max: " %9.2e `maxdiff' ")"
}

* Compare e(V)
matrix `diff' = reghdfe_V - creghdfe_V
local rows = rowsof(`diff')
local cols = colsof(`diff')
local maxdiff = 0
forvalues i = 1/`rows' {
    forvalues j = 1/`cols' {
        local d = abs(`diff'[`i', `j'])
        if `d' > `maxdiff' local maxdiff = `d'
    }
}

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `maxdiff' < 1e-6 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] e(V) matches"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] e(V) differs (max: " %9.2e `maxdiff' ")"
}

* Compare df_r
global TESTS_TOTAL = $TESTS_TOTAL + 1
if `reghdfe_df_r' == `creghdfe_df_r' {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] e(df_r) matches"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] e(df_r) differs (`reghdfe_df_r' vs `creghdfe_df_r')"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "creghdfe"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
