/*******************************************************************************
 * validate_creghdfe.do
 *
 * Comprehensive validation tests for creghdfe vs reghdfe
 * Tests all options: absorb, vce, weights, tolerance, maxiter, resid
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CREGHDFE VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign)
if _rc != 0 {
    noi test_fail "creghdfe plugin load" "returned error `=_rc'"
    noi print_summary "creghdfe"
    exit 1
}
noi test_pass "creghdfe plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Auto dataset - basic FE tests
 ******************************************************************************/
noi print_section "Basic Fixed Effects (auto)"

sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign) testname("single FE")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
noi benchmark_reghdfe price mpg, absorb(foreign) testname("single covariate")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight length, absorb(foreign) testname("three covariates")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight length turn displacement, absorb(foreign) testname("many covariates")

/*******************************************************************************
 * SECTION 3: VCE options
 ******************************************************************************/
noi print_section "VCE Options"

sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign) vce(robust) testname("vce(robust)")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster foreign) testname("vce(cluster)")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign rep78) vce(robust) testname("two-way + robust")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign rep78) vce(cluster foreign) testname("two-way + cluster")

/*******************************************************************************
 * SECTION 4: Weight tests
 ******************************************************************************/
noi print_section "Weights"

sysuse auto, clear
noi benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign) testname("aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
noi benchmark_reghdfe price mpg weight [fw=fw], absorb(foreign) testname("fweight")

sysuse auto, clear
noi benchmark_reghdfe price mpg weight [pw=weight], absorb(foreign) testname("pweight")

/*******************************************************************************
 * SECTION 5: Census dataset
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
noi benchmark_reghdfe pop medage, absorb(region) testname("single FE")

sysuse census, clear
noi benchmark_reghdfe pop medage death, absorb(region) testname("two covariates")

sysuse census, clear
noi benchmark_reghdfe pop medage death marriage divorce, absorb(region) testname("many covariates")

sysuse census, clear
noi benchmark_reghdfe pop medage, absorb(region) vce(robust) testname("robust")

/*******************************************************************************
 * SECTION 6: nlswork panel data
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000
noi benchmark_reghdfe ln_wage age tenure, absorb(idcode) testname("individual FE")

webuse nlswork, clear
keep in 1/10000
noi benchmark_reghdfe ln_wage age tenure, absorb(idcode year) testname("two-way FE")

webuse nlswork, clear
keep in 1/10000
noi benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) testname("three covariates")

webuse nlswork, clear
keep in 1/10000
noi benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(robust) testname("robust")

webuse nlswork, clear
keep in 1/10000
noi benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(cluster idcode) testname("cluster idcode")

/*******************************************************************************
 * SECTION 7: Large dataset
 ******************************************************************************/
noi print_section "Large Dataset"

clear
set seed 12345
set obs 50000
gen id = runiformint(1, 500)
gen year = 2000 + runiformint(0, 10)
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 100)
gen y = 2*x1 + 3*x2 - 0.5*x3 + rnormal()

noi benchmark_reghdfe y x1 x2 x3, absorb(id) testname("50K single FE")

noi benchmark_reghdfe y x1 x2 x3, absorb(id year) testname("50K two-way FE")

noi benchmark_reghdfe y x1 x2 x3, absorb(id) vce(robust) testname("50K robust")

/*******************************************************************************
 * SECTION 8: tolerance/maxiter options
 ******************************************************************************/
noi print_section "tolerance/maxiter Options"

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) tolerance(1e-10)
if _rc == 0 {
    noi test_pass "tolerance(1e-10) accepted"
}
else {
    noi test_fail "tolerance option" "returned error `=_rc'"
}

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) maxiter(1000)
if _rc == 0 {
    noi test_pass "maxiter(1000) accepted"
}
else {
    noi test_fail "maxiter option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 9: resid option
 ******************************************************************************/
noi print_section "resid Option"

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) resid
if _rc == 0 {
    capture confirm variable _reghdfe_resid
    if _rc == 0 {
        noi test_pass "resid option creates residual variable"
    }
    else {
        noi test_fail "resid option" "residual variable not created"
    }
}
else {
    noi test_fail "resid option" "returned error `=_rc'"
}

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) resid2(myresid)
if _rc == 0 {
    capture confirm variable myresid
    if _rc == 0 {
        noi test_pass "resid2(name) creates named residual"
    }
    else {
        noi test_fail "resid2(name)" "residual variable not created"
    }
}
else {
    noi test_fail "resid2 option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 10: verbose/timeit options
 ******************************************************************************/
noi print_section "verbose/timeit Options"

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) verbose
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) timeit
if _rc == 0 {
    noi test_pass "timeit option accepted"
}
else {
    noi test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: if/in conditions
 ******************************************************************************/
noi print_section "if/in Conditions"

sysuse auto, clear
reghdfe price mpg weight if price > 5000, absorb(foreign)
matrix reghdfe_b = e(b)
local reghdfe_N = e(N)

creghdfe price mpg weight if price > 5000, absorb(foreign)
matrix creghdfe_b = e(b)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "if condition: N matches"
}
else {
    noi test_fail "if condition" "N differs"
}

sysuse auto, clear
reghdfe price mpg weight in 1/50, absorb(foreign)
local reghdfe_N = e(N)

creghdfe price mpg weight in 1/50, absorb(foreign)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "in condition: N matches"
}
else {
    noi test_fail "in condition" "N differs"
}

/*******************************************************************************
 * SECTION 12: nostandardize option
 ******************************************************************************/
noi print_section "nostandardize Option"

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) nostandardize
if _rc == 0 {
    noi test_pass "nostandardize option accepted"
}
else {
    noi test_fail "nostandardize option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 13: Edge cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Single observation per FE group
clear
set seed 99
set obs 100
gen id = _n
gen x = runiform()
gen y = x + rnormal()

capture creghdfe y x, absorb(id)
* With 100 obs and 100 unique IDs, all are singletons - expect error 2001 or graceful handling
if _rc == 0 | _rc == 2001 {
    noi test_pass "singleton FE handling (rc=`=_rc')"
}
else {
    noi test_fail "singleton FE" "returned unexpected error `=_rc'"
}

* Many FE levels
clear
set seed 123
set obs 10000
gen id = runiformint(1, 1000)
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()

reghdfe y x1 x2, absorb(id)
local reghdfe_N = e(N)

creghdfe y x1 x2, absorb(id)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "many FE levels (1000): N matches"
}
else {
    noi test_fail "many FE levels" "N differs"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "creghdfe"

if $TESTS_FAILED > 0 {
    exit 1
}

}
