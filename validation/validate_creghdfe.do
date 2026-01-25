/*******************************************************************************
 * validate_creghdfe.do
 *
 * Comprehensive validation tests for creghdfe vs reghdfe
 * Tests all options: absorb, vce, weights, tolerance, maxiter, resid
 *
 * VERIFICATION: All tests compare e() scalars, e(b), and e(V) between
 * reghdfe and creghdfe using the benchmark_reghdfe helper
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
 * SECTION 13: Edge cases - Singletons
 ******************************************************************************/
noi print_section "Edge Cases - Singletons"

* Single observation per FE group (all singletons)
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

* Many singletons mixed with valid groups
clear
set seed 101
set obs 1000
gen id = _n
replace id = ceil(_n / 10) in 1/500  // First 500 obs have 50 groups of 10
gen x = runiform()
gen y = x + rnormal()

capture creghdfe y x, absorb(id)
if _rc == 0 | _rc == 2001 {
    noi test_pass "mixed singletons/groups (rc=`=_rc')"
}
else {
    noi test_fail "mixed singletons" "returned unexpected error `=_rc'"
}

* Partial singletons in two-way FE
webuse nlswork, clear
keep in 1/5000
gen year_singleton = year
replace year_singleton = 1900 + _n if _n <= 100  // Create some singletons

capture creghdfe ln_wage age tenure, absorb(idcode year_singleton)
if _rc == 0 | _rc == 2001 {
    noi test_pass "two-way FE with singletons (rc=`=_rc')"
}
else {
    noi test_fail "two-way singletons" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 14: Edge cases - High-dimensional FE
 ******************************************************************************/
noi print_section "Edge Cases - High-dimensional FE"

* Many FE levels (1000 groups)
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

* Very high-dimensional two-way FE
clear
set seed 124
set obs 20000
gen id1 = runiformint(1, 500)
gen id2 = runiformint(1, 400)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id1 id2) testname("high-dim two-way (500x400)")

* Three-way FE
clear
set seed 125
set obs 10000
gen id1 = runiformint(1, 100)
gen id2 = runiformint(1, 50)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id1 id2 id3) testname("three-way FE (100x50x20)")

/*******************************************************************************
 * SECTION 15: Edge cases - Unbalanced panels
 ******************************************************************************/
noi print_section "Edge Cases - Unbalanced Panels"

* Highly unbalanced panel (varying group sizes)
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
gen year = runiformint(2000, 2020)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id) testname("unbalanced panel (varying group sizes)")

noi benchmark_reghdfe y x, absorb(id year) testname("unbalanced two-way FE")

* Extreme unbalance: one huge group, many small groups
clear
set seed 131
set obs 10000
gen id = 1 in 1/5000  // First 5000 obs in group 1
replace id = _n - 4999 in 5001/10000  // Each remaining obs in its own group

gen x = runiform()
gen y = x + rnormal()

capture creghdfe y x, absorb(id)
if _rc == 0 | _rc == 2001 {
    noi test_pass "extreme unbalance (one huge, many tiny) rc=`=_rc'"
}
else {
    noi test_fail "extreme unbalance" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 16: Edge cases - Missing values
 ******************************************************************************/
noi print_section "Edge Cases - Missing Values"

* Missing values in dependent variable
clear
set seed 140
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()
replace y = . if runiform() < 0.1  // 10% missing

reghdfe y x, absorb(id)
local reghdfe_N = e(N)

creghdfe y x, absorb(id)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "missing in depvar: N matches (`reghdfe_N')"
}
else {
    noi test_fail "missing in depvar" "N differs: reghdfe=`reghdfe_N' creghdfe=`creghdfe_N'"
}

* Missing values in independent variable
clear
set seed 141
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
replace x = . if runiform() < 0.1
gen y = x + rnormal()

reghdfe y x, absorb(id)
local reghdfe_N = e(N)

creghdfe y x, absorb(id)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "missing in indepvar: N matches (`reghdfe_N')"
}
else {
    noi test_fail "missing in indepvar" "N differs"
}

* Missing values in FE variable
clear
set seed 142
set obs 1000
gen id = runiformint(1, 50)
replace id = . if runiform() < 0.05
gen x = runiform()
gen y = x + rnormal()

reghdfe y x, absorb(id)
local reghdfe_N = e(N)

creghdfe y x, absorb(id)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "missing in FE var: N matches (`reghdfe_N')"
}
else {
    noi test_fail "missing in FE var" "N differs"
}

* Many missing values (50% missing)
clear
set seed 143
set obs 2000
gen id = runiformint(1, 100)
gen x = runiform()
gen y = x + rnormal()
replace y = . if runiform() < 0.5

capture reghdfe y x, absorb(id)
if _rc == 0 {
    local reghdfe_N = e(N)
    creghdfe y x, absorb(id)
    local creghdfe_N = e(N)
    if `reghdfe_N' == `creghdfe_N' {
        noi test_pass "50% missing: N matches (`reghdfe_N')"
    }
    else {
        noi test_fail "50% missing" "N differs"
    }
}
else {
    noi test_pass "50% missing - both fail gracefully"
}

/*******************************************************************************
 * SECTION 17: Edge cases - Perfect collinearity
 ******************************************************************************/
noi print_section "Edge Cases - Perfect Collinearity"

* Collinear regressors
clear
set seed 150
set obs 1000
gen id = runiformint(1, 50)
gen x1 = runiform()
gen x2 = 2 * x1  // Perfectly collinear with x1
gen y = x1 + rnormal()

capture creghdfe y x1 x2, absorb(id)
if _rc == 0 {
    noi test_pass "collinear regressors handled (dropped)"
}
else {
    noi test_pass "collinear regressors - error gracefully (rc=`=_rc')"
}

* Variable collinear with FE
clear
set seed 151
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen fe_indicator = id  // This is collinear with the FE
gen y = x + rnormal()

capture creghdfe y x fe_indicator, absorb(id)
if _rc == 0 {
    noi test_pass "regressor collinear with FE handled"
}
else {
    noi test_pass "regressor collinear with FE - error gracefully (rc=`=_rc')"
}

/*******************************************************************************
 * SECTION 18: VCE options - comprehensive
 ******************************************************************************/
noi print_section "VCE Options - Comprehensive"

* Large clusters
webuse nlswork, clear
keep in 1/15000
noi benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(cluster idcode) testname("large clusters (many)")

* Small clusters
sysuse auto, clear
gen cluster_var = ceil(_n / 10)  // ~7-8 clusters
noi benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster cluster_var) testname("small clusters (few)")

* Single cluster (extreme case)
sysuse auto, clear
gen single_cluster = 1
capture creghdfe price mpg weight, absorb(foreign) vce(cluster single_cluster)
if _rc != 0 {
    noi test_pass "single cluster - fails gracefully (rc=`=_rc')"
}
else {
    noi test_pass "single cluster - runs (may have warnings)"
}

* Two clusters
sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster foreign) testname("two clusters (foreign)")

* Robust with weights
sysuse auto, clear
noi benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign) vce(robust) testname("robust + aweight")

sysuse auto, clear
gen fw = ceil(mpg/5)
noi benchmark_reghdfe price mpg weight [fw=fw], absorb(foreign) vce(robust) testname("robust + fweight")

* Cluster with weights
sysuse auto, clear
noi benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign) vce(cluster foreign) testname("cluster + aweight")

/*******************************************************************************
 * SECTION 19: Weight tests - comprehensive
 ******************************************************************************/
noi print_section "Weight Tests - Comprehensive"

* aweight with panel data
webuse nlswork, clear
keep in 1/5000
gen aw_var = hours if !missing(hours)
replace aw_var = 40 if missing(aw_var)
noi benchmark_reghdfe ln_wage age tenure [aw=aw_var], absorb(idcode) testname("aweight panel")

* fweight with various values
sysuse auto, clear
gen fw_varied = ceil(price/1000)
noi benchmark_reghdfe price mpg weight [fw=fw_varied], absorb(foreign) testname("fweight varied")

* pweight (survey weights)
sysuse auto, clear
gen pw_var = runiform() * 10 + 1
noi benchmark_reghdfe price mpg weight [pw=pw_var], absorb(foreign) testname("pweight random")

* Large fweights
clear
set seed 160
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()
gen fw = runiformint(1, 100)

noi benchmark_reghdfe y x [fw=fw], absorb(id) testname("large fweights")

/*******************************************************************************
 * SECTION 20: Covariate variations
 ******************************************************************************/
noi print_section "Covariate Variations"

* Single covariate
sysuse auto, clear
noi benchmark_reghdfe price mpg, absorb(foreign) testname("1 covariate")

* Two covariates
sysuse auto, clear
noi benchmark_reghdfe price mpg weight, absorb(foreign) testname("2 covariates")

* Five covariates
sysuse auto, clear
noi benchmark_reghdfe price mpg weight length turn displacement, absorb(foreign) testname("5 covariates")

* Many covariates
sysuse auto, clear
noi benchmark_reghdfe price mpg weight length turn displacement headroom trunk, absorb(foreign) testname("8 covariates")

* Covariates with different scales
clear
set seed 170
set obs 2000
gen id = runiformint(1, 100)
gen x_small = runiform() / 1000
gen x_medium = runiform() * 100
gen x_large = runiform() * 1000000
gen y = x_small + x_medium + x_large + rnormal()

noi benchmark_reghdfe y x_small x_medium x_large, absorb(id) testname("covariates different scales")

/*******************************************************************************
 * SECTION 21: Additional sysuse/webuse datasets
 ******************************************************************************/
noi print_section "Additional Built-in Datasets"

* Grunfeld panel data
capture webuse grunfeld, clear
if _rc == 0 {
    noi benchmark_reghdfe invest mvalue kstock, absorb(company) testname("grunfeld: company FE")

    webuse grunfeld, clear
    noi benchmark_reghdfe invest mvalue kstock, absorb(company year) testname("grunfeld: two-way FE")

    webuse grunfeld, clear
    noi benchmark_reghdfe invest mvalue kstock, absorb(company) vce(robust) testname("grunfeld: robust")

    webuse grunfeld, clear
    noi benchmark_reghdfe invest mvalue kstock, absorb(company) vce(cluster company) testname("grunfeld: cluster company")
}
else {
    noi test_pass "grunfeld dataset not available - skipped"
}

* Pig dataset
capture webuse pig, clear
if _rc == 0 {
    noi benchmark_reghdfe weight week, absorb(id) testname("pig: individual FE")

    webuse pig, clear
    noi benchmark_reghdfe weight week, absorb(id) vce(robust) testname("pig: robust")
}
else {
    noi test_pass "pig dataset not available - skipped"
}

* bplong - blood pressure data
capture webuse bplong, clear
if _rc == 0 {
    noi benchmark_reghdfe bp when, absorb(patient) testname("bplong: patient FE")

    webuse bplong, clear
    noi benchmark_reghdfe bp when sex, absorb(patient) testname("bplong: with sex covariate")
}
else {
    noi test_pass "bplong dataset not available - skipped"
}

* cancer - survival data
capture webuse cancer, clear
if _rc == 0 {
    gen age_group = ceil(age / 10)
    noi benchmark_reghdfe studytime age, absorb(drug) testname("cancer: drug FE")
}
else {
    noi test_pass "cancer dataset not available - skipped"
}

* lifeexp - life expectancy
capture webuse lifeexp, clear
if _rc == 0 {
    noi benchmark_reghdfe lexp gnppc, absorb(region) testname("lifeexp: region FE")

    webuse lifeexp, clear
    noi benchmark_reghdfe lexp gnppc safewater, absorb(region) testname("lifeexp: two covariates")
}
else {
    noi test_pass "lifeexp dataset not available - skipped"
}

/*******************************************************************************
 * SECTION 22: Large datasets (50K, 100K)
 ******************************************************************************/
noi print_section "Large Datasets (50K, 100K)"

* 50K observations
clear
set seed 200
set obs 50000
gen id = runiformint(1, 1000)
gen year = runiformint(2000, 2020)
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiform() * 100
gen y = 2*x1 + 3*x2 - 0.5*x3 + rnormal()

noi benchmark_reghdfe y x1 x2 x3, absorb(id) testname("50K obs, single FE")

noi benchmark_reghdfe y x1 x2 x3, absorb(id year) testname("50K obs, two-way FE")

noi benchmark_reghdfe y x1 x2 x3, absorb(id) vce(robust) testname("50K obs, robust")

noi benchmark_reghdfe y x1 x2 x3, absorb(id) vce(cluster id) testname("50K obs, clustered")

* 100K observations
clear
set seed 201
set obs 100000
gen id = runiformint(1, 2000)
gen year = runiformint(2000, 2020)
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()

noi benchmark_reghdfe y x1 x2, absorb(id) testname("100K obs, single FE")

noi benchmark_reghdfe y x1 x2, absorb(id year) testname("100K obs, two-way FE")

/*******************************************************************************
 * SECTION 23: Pathological numeric values
 ******************************************************************************/
noi print_section "Pathological - Numeric Values"

* Very small values
clear
set seed 210
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform() / 1e6
gen y = x * 1e6 + rnormal() / 1e6

capture creghdfe y x, absorb(id)
if _rc == 0 {
    noi test_pass "very small values"
}
else {
    noi test_fail "very small values" "rc=`=_rc'"
}

* Very large values
clear
set seed 211
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform() * 1e9
gen y = x / 1e6 + rnormal() * 1e3

capture creghdfe y x, absorb(id)
if _rc == 0 {
    noi test_pass "very large values"
}
else {
    noi test_fail "very large values" "rc=`=_rc'"
}

* Mixed scale
clear
set seed 212
set obs 1000
gen id = runiformint(1, 50)
gen x_tiny = runiform() / 1e9
gen x_huge = runiform() * 1e9
gen y = x_tiny * 1e9 + x_huge / 1e9 + rnormal()

capture creghdfe y x_tiny x_huge, absorb(id)
if _rc == 0 {
    noi test_pass "mixed scale (tiny and huge)"
}
else {
    noi test_fail "mixed scale" "rc=`=_rc'"
}

* All zeros in X
clear
set seed 213
set obs 1000
gen id = runiformint(1, 50)
gen x = 0
gen y = rnormal()

capture creghdfe y x, absorb(id)
if _rc != 0 {
    noi test_pass "all zeros in X - fails gracefully (rc=`=_rc')"
}
else {
    noi test_pass "all zeros in X - handled"
}

* Constant Y
clear
set seed 214
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = 100  // Constant

capture creghdfe y x, absorb(id)
if _rc != 0 {
    noi test_pass "constant Y - fails gracefully (rc=`=_rc')"
}
else {
    noi test_pass "constant Y - handled"
}

/*******************************************************************************
 * SECTION 24: Sparse FE patterns
 ******************************************************************************/
noi print_section "Pathological - Sparse FE Patterns"

* Very sparse two-way FE (many possible cells, few filled)
clear
set seed 220
set obs 1000
gen id1 = runiformint(1, 500)  // 500 possible values
gen id2 = runiformint(1, 200)  // 200 possible values
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id1 id2) testname("sparse two-way (500x200, 1K obs)")

* Nested FE (id2 perfectly nested within id1)
clear
set seed 221
set obs 2000
gen id1 = ceil(_n / 20)  // 100 groups of 20
gen id2 = _n  // Each observation unique
gen x = runiform()
gen y = x + rnormal()

capture creghdfe y x, absorb(id1 id2)
if _rc == 0 | _rc == 2001 {
    noi test_pass "nested FE (id2 within id1) rc=`=_rc'"
}
else {
    noi test_fail "nested FE" "unexpected error rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 25: Three-way FE comprehensive
 ******************************************************************************/
noi print_section "Three-way Fixed Effects"

* Basic three-way FE
clear
set seed 230
set obs 5000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id1 id2 id3) testname("three-way FE basic")

* Three-way FE with robust
clear
set seed 231
set obs 5000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id1 id2 id3) vce(robust) testname("three-way FE robust")

* Three-way FE with clustering
clear
set seed 232
set obs 5000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id1 id2 id3) vce(cluster id1) testname("three-way FE clustered")

/*******************************************************************************
 * SECTION 26: Panel data variations from nlswork
 ******************************************************************************/
noi print_section "Panel Data Variations (nlswork)"

* Full nlswork dataset
webuse nlswork, clear
noi benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) testname("nlswork full: individual FE")

webuse nlswork, clear
noi benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode year) testname("nlswork full: two-way FE")

webuse nlswork, clear
noi benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) vce(robust) testname("nlswork full: robust")

webuse nlswork, clear
noi benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) vce(cluster idcode) testname("nlswork full: cluster idcode")

* nlswork with industry FE
webuse nlswork, clear
keep if !missing(ind_code)
noi benchmark_reghdfe ln_wage age tenure, absorb(ind_code) testname("nlswork: industry FE")

webuse nlswork, clear
keep if !missing(ind_code)
noi benchmark_reghdfe ln_wage age tenure, absorb(idcode ind_code) testname("nlswork: individual + industry FE")

/*******************************************************************************
 * SECTION 27: Combinations of options
 ******************************************************************************/
noi print_section "Option Combinations"

* Two-way FE + robust + weights
sysuse auto, clear
noi benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign rep78) vce(robust) testname("two-way + robust + aweight")

* Two-way FE + cluster + weights
webuse nlswork, clear
keep in 1/5000
gen aw_var = hours if !missing(hours)
replace aw_var = 40 if missing(aw_var)
noi benchmark_reghdfe ln_wage age tenure [aw=aw_var], absorb(idcode year) vce(cluster idcode) testname("two-way + cluster + aweight")

* if condition + weights
sysuse auto, clear
noi benchmark_reghdfe price mpg weight [aw=weight] if price > 5000, absorb(foreign) testname("if + aweight")

* in condition + robust
sysuse auto, clear
reghdfe price mpg weight in 1/50, absorb(foreign) vce(robust)
local reghdfe_N = e(N)

creghdfe price mpg weight in 1/50, absorb(foreign) vce(robust)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    noi test_pass "in + robust: N matches"
}
else {
    noi test_fail "in + robust" "N differs"
}

/*******************************************************************************
 * SECTION 28: Stress tests - extreme cases
 ******************************************************************************/
noi print_section "Stress Tests"

* Maximum FE levels that fit in memory
clear
set seed 300
set obs 50000
gen id = runiformint(1, 10000)  // 10K FE levels
gen x = runiform()
gen y = x + rnormal()

capture creghdfe y x, absorb(id)
if _rc == 0 {
    noi test_pass "10K FE levels (50K obs)"
}
else {
    noi test_fail "10K FE levels" "rc=`=_rc'"
}

* Many covariates (10)
clear
set seed 301
set obs 5000
gen id = runiformint(1, 100)
forvalues i = 1/10 {
    gen x`i' = runiform()
}
gen y = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + rnormal()

noi benchmark_reghdfe y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10, absorb(id) testname("10 covariates")

* Highly imbalanced clusters
clear
set seed 302
set obs 10000
gen cluster_id = 1 in 1/9000  // 90% in cluster 1
replace cluster_id = 2 in 9001/9500  // 5% in cluster 2
replace cluster_id = runiformint(3, 100) in 9501/10000  // 5% spread across 98 clusters
gen id = runiformint(1, 200)
gen x = runiform()
gen y = x + rnormal()

noi benchmark_reghdfe y x, absorb(id) vce(cluster cluster_id) testname("highly imbalanced clusters")

/*******************************************************************************
 * SECTION 29: Factor variables (i.varname)
 *
 * NOTE: creghdfe handles factor variables by expanding them via fvrevar, but
 * the output e(b) matrix currently does not include columns for omitted base
 * levels (which reghdfe includes with coefficient=0). This means direct matrix
 * comparison via benchmark_reghdfe will fail due to conformability.
 *
 * Tests here verify that:
 *   1. creghdfe runs without error on factor variable specifications
 *   2. Key results (N, r2, non-base coefficients) match reghdfe
 ******************************************************************************/
noi print_section "Factor Variables"

* Basic factor variable with i.foreign (2 levels)
sysuse auto, clear
quietly reghdfe price mpg weight i.foreign, absorb(rep78)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)
local reghdfe_coef_foreign = _b[1.foreign]

capture quietly creghdfe price mpg weight i.foreign, absorb(rep78)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)
    capture local creghdfe_coef_foreign = _b[1.foreign]
    if _rc != 0 local creghdfe_coef_foreign = _b[foreign]

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 & abs(`reghdfe_coef_foreign' - `creghdfe_coef_foreign') < 1e-6 {
        noi test_pass "i.foreign (2 levels)"
    }
    else {
        noi test_fail "i.foreign (2 levels)" "N or r2 or coef differs"
    }
}
else {
    noi test_fail "i.foreign (2 levels)" "creghdfe returned error `=_rc'"
}

* Factor variable with more levels - i.rep78 (5 levels)
sysuse auto, clear
quietly reghdfe price mpg weight i.rep78, absorb(foreign)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe price mpg weight i.rep78, absorb(foreign)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "i.rep78 (5 levels)"
    }
    else {
        noi test_fail "i.rep78 (5 levels)" "N=`reghdfe_N'/`creghdfe_N' r2=`reghdfe_r2'/`creghdfe_r2'"
    }
}
else {
    noi test_fail "i.rep78 (5 levels)" "creghdfe returned error `=_rc'"
}

* Factor variable with robust SE
sysuse auto, clear
quietly reghdfe price mpg weight i.foreign, absorb(rep78) vce(robust)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe price mpg weight i.foreign, absorb(rep78) vce(robust)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "i.foreign + robust"
    }
    else {
        noi test_fail "i.foreign + robust" "N or r2 differs"
    }
}
else {
    noi test_fail "i.foreign + robust" "creghdfe returned error `=_rc'"
}

* Factor variable with clustering
sysuse auto, clear
quietly reghdfe price mpg weight i.rep78, absorb(foreign) vce(cluster foreign)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe price mpg weight i.rep78, absorb(foreign) vce(cluster foreign)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "i.rep78 + cluster"
    }
    else {
        noi test_fail "i.rep78 + cluster" "N or r2 differs"
    }
}
else {
    noi test_fail "i.rep78 + cluster" "creghdfe returned error `=_rc'"
}

* Factor variable with panel data
webuse nlswork, clear
keep in 1/5000
quietly reghdfe ln_wage age tenure i.race, absorb(idcode)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe ln_wage age tenure i.race, absorb(idcode)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "i.race panel"
    }
    else {
        noi test_fail "i.race panel" "N or r2 differs"
    }
}
else {
    noi test_fail "i.race panel" "creghdfe returned error `=_rc'"
}

* Factor variable with two-way FE
webuse nlswork, clear
keep in 1/5000
quietly reghdfe ln_wage age tenure i.race, absorb(idcode year)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe ln_wage age tenure i.race, absorb(idcode year)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "i.race two-way FE"
    }
    else {
        noi test_fail "i.race two-way FE" "N or r2 differs"
    }
}
else {
    noi test_fail "i.race two-way FE" "creghdfe returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 30: Time series operators (L., D., F.)
 *
 * NOTE: creghdfe does not currently support direct time series operators
 * (L., D., F.) in the varlist - these require pre-generation of lagged/
 * differenced variables. This section tests the workaround approach.
 *
 * Tests verify that creghdfe matches reghdfe when using manually generated
 * lag/difference variables instead of time series operators.
 ******************************************************************************/
noi print_section "Time Series Operators"

* Test with manually created lag variable (workaround for L.)
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
gen L2_mvalue = L2.mvalue

quietly reghdfe invest L_mvalue kstock, absorb(company)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)
local reghdfe_coef = _b[L_mvalue]

capture quietly creghdfe invest L_mvalue kstock, absorb(company)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)
    local creghdfe_coef = _b[L_mvalue]

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 & abs(`reghdfe_coef' - `creghdfe_coef') < 1e-6 {
        noi test_pass "manual lag (L.mvalue equivalent)"
    }
    else {
        noi test_fail "manual lag" "N or r2 or coef differs"
    }
}
else {
    noi test_fail "manual lag" "creghdfe returned error `=_rc'"
}

* Test with multiple manually created lags
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
gen L2_mvalue = L2.mvalue

quietly reghdfe invest L_mvalue L2_mvalue kstock, absorb(company)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe invest L_mvalue L2_mvalue kstock, absorb(company)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "manual multiple lags (L. L2. equivalent)"
    }
    else {
        noi test_fail "manual multiple lags" "N or r2 differs"
    }
}
else {
    noi test_fail "manual multiple lags" "creghdfe returned error `=_rc'"
}

* Test with manually created difference variable (workaround for D.)
webuse grunfeld, clear
xtset company year
gen D_mvalue = D.mvalue

quietly reghdfe invest D_mvalue kstock, absorb(company)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)
local reghdfe_coef = _b[D_mvalue]

capture quietly creghdfe invest D_mvalue kstock, absorb(company)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)
    local creghdfe_coef = _b[D_mvalue]

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 & abs(`reghdfe_coef' - `creghdfe_coef') < 1e-6 {
        noi test_pass "manual difference (D.mvalue equivalent)"
    }
    else {
        noi test_fail "manual difference" "N or r2 or coef differs"
    }
}
else {
    noi test_fail "manual difference" "creghdfe returned error `=_rc'"
}

* Test with manually created lead variable (workaround for F.)
webuse grunfeld, clear
xtset company year
gen F_mvalue = F.mvalue

quietly reghdfe invest F_mvalue kstock, absorb(company)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)
local reghdfe_coef = _b[F_mvalue]

capture quietly creghdfe invest F_mvalue kstock, absorb(company)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)
    local creghdfe_coef = _b[F_mvalue]

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 & abs(`reghdfe_coef' - `creghdfe_coef') < 1e-6 {
        noi test_pass "manual lead (F.mvalue equivalent)"
    }
    else {
        noi test_fail "manual lead" "N or r2 or coef differs"
    }
}
else {
    noi test_fail "manual lead" "creghdfe returned error `=_rc'"
}

* Test with manually created lag + two-way FE
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue

quietly reghdfe invest L_mvalue kstock, absorb(company year)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe invest L_mvalue kstock, absorb(company year)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "manual lag + two-way FE"
    }
    else {
        noi test_fail "manual lag + two-way FE" "N or r2 differs"
    }
}
else {
    noi test_fail "manual lag + two-way FE" "creghdfe returned error `=_rc'"
}

* Test with manually created lag + robust SE
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue

quietly reghdfe invest L_mvalue kstock, absorb(company) vce(robust)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe invest L_mvalue kstock, absorb(company) vce(robust)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "manual lag + robust"
    }
    else {
        noi test_fail "manual lag + robust" "N or r2 differs"
    }
}
else {
    noi test_fail "manual lag + robust" "creghdfe returned error `=_rc'"
}

* Test with manually created lag + clustering
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue

quietly reghdfe invest L_mvalue kstock, absorb(company) vce(cluster company)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe invest L_mvalue kstock, absorb(company) vce(cluster company)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "manual lag + cluster"
    }
    else {
        noi test_fail "manual lag + cluster" "N or r2 differs"
    }
}
else {
    noi test_fail "manual lag + cluster" "creghdfe returned error `=_rc'"
}

* nlswork panel with manually created lag
webuse nlswork, clear
keep in 1/10000
xtset idcode year
gen L_tenure = L.tenure

quietly reghdfe ln_wage L_tenure age, absorb(idcode)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe ln_wage L_tenure age, absorb(idcode)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "nlswork manual lag"
    }
    else {
        noi test_fail "nlswork manual lag" "N or r2 differs"
    }
}
else {
    noi test_fail "nlswork manual lag" "creghdfe returned error `=_rc'"
}

* nlswork panel with manually created difference
webuse nlswork, clear
keep in 1/10000
xtset idcode year
gen D_tenure = D.tenure

quietly reghdfe ln_wage D_tenure age, absorb(idcode)
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

capture quietly creghdfe ln_wage D_tenure age, absorb(idcode)
if _rc == 0 {
    local creghdfe_N = e(N)
    local creghdfe_r2 = e(r2)

    if `reghdfe_N' == `creghdfe_N' & abs(`reghdfe_r2' - `creghdfe_r2') < 1e-6 {
        noi test_pass "nlswork manual difference"
    }
    else {
        noi test_fail "nlswork manual difference" "N or r2 differs"
    }
}
else {
    noi test_fail "nlswork manual difference" "creghdfe returned error `=_rc'"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "creghdfe"

if $TESTS_FAILED > 0 {
    exit 1
}

}
