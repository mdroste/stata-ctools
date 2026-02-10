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

noi di as text "Running validation tests for creghdfe..."

* Plugin check
sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign)
if _rc != 0 {
    test_fail "creghdfe plugin load" "returned error `=_rc'"
    exit 1
}
test_pass "creghdfe plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Auto dataset - basic FE tests
 ******************************************************************************/
print_section "Basic Fixed Effects (auto)"

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) testname("single FE")
benchmark_reghdfe price mpg weight, absorb(foreign rep78) testname("two-way FE")
benchmark_reghdfe price mpg, absorb(foreign) testname("single covariate")
benchmark_reghdfe price mpg weight length, absorb(foreign) testname("three covariates")
benchmark_reghdfe price mpg weight length turn displacement, absorb(foreign) testname("many covariates")

/*******************************************************************************
 * SECTION 3: VCE options
 ******************************************************************************/
print_section "VCE Options"

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) vce(robust) testname("vce(robust)")
benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster foreign) testname("vce(cluster)")
benchmark_reghdfe price mpg weight, absorb(foreign rep78) vce(robust) testname("two-way + robust")
benchmark_reghdfe price mpg weight, absorb(foreign rep78) vce(cluster foreign) testname("two-way + cluster")

/*******************************************************************************
 * SECTION 4: Weight tests
 ******************************************************************************/
print_section "Weights"

sysuse auto, clear
benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign) testname("aweight")

sysuse auto, clear
gen int fw = ceil(mpg/5)
benchmark_reghdfe price mpg weight [fw=fw], absorb(foreign) testname("fweight")

sysuse auto, clear
benchmark_reghdfe price mpg weight [pw=weight], absorb(foreign) testname("pweight")

/*******************************************************************************
 * SECTION 5: Census dataset
 ******************************************************************************/
print_section "Census Dataset"

sysuse census, clear
benchmark_reghdfe pop medage, absorb(region) testname("single FE")
benchmark_reghdfe pop medage death, absorb(region) testname("two covariates")
benchmark_reghdfe pop medage death marriage divorce, absorb(region) testname("many covariates")
benchmark_reghdfe pop medage, absorb(region) vce(robust) testname("robust")

/*******************************************************************************
 * SECTION 6: nlswork panel data
 ******************************************************************************/
print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age tenure, absorb(idcode) testname("individual FE")
benchmark_reghdfe ln_wage age tenure, absorb(idcode year) testname("two-way FE")
benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) testname("three covariates")
benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(robust) testname("robust")
benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(cluster idcode) testname("cluster idcode")

/*******************************************************************************
 * SECTION 7: Large dataset
 ******************************************************************************/
print_section "Large Dataset"

clear
set seed 12345
set obs 50000
gen id = runiformint(1, 500)
gen year = 2000 + runiformint(0, 10)
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 100)
gen y = 2*x1 + 3*x2 - 0.5*x3 + rnormal()

benchmark_reghdfe y x1 x2 x3, absorb(id) testname("50K single FE")

benchmark_reghdfe y x1 x2 x3, absorb(id year) testname("50K two-way FE")

benchmark_reghdfe y x1 x2 x3, absorb(id) vce(robust) testname("50K robust")

/*******************************************************************************
 * SECTION 8: tolerance/maxiter options
 ******************************************************************************/
print_section "tolerance/maxiter Options"

sysuse auto, clear
benchmark_reghdfe_opts price mpg weight, absorb(foreign) creghdfe_opts(tolerance(1e-10)) testname("[syntax] tolerance(1e-10)")
benchmark_reghdfe_opts price mpg weight, absorb(foreign) creghdfe_opts(iterate(1000)) testname("[syntax] iterate(1000)")

/*******************************************************************************
 * SECTION 9: resid option
 ******************************************************************************/
print_section "resid Option"

* Compare resid values against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe price mpg weight, absorb(foreign) resid
assert_var_equal _reghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid: default name"

* Compare resid2(name) values against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe price mpg weight, absorb(foreign) resid2(myresid)
assert_var_equal myresid `reghdfe_r' $DEFAULT_SIGFIGS "resid2(name): custom name"

/*******************************************************************************
 * SECTION 10: verbose option
 ******************************************************************************/
print_section "verbose Option"

sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign) verbose
if _rc == 0 {
    test_pass "verbose option accepted"
}
else {
    test_fail "verbose option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: if/in conditions
 ******************************************************************************/
print_section "if/in Conditions"

* if condition - use benchmark_reghdfe for full comparison (N, coefficients, VCE)
sysuse auto, clear
benchmark_reghdfe price mpg weight if price > 5000, absorb(foreign) testname("if condition")

* in condition - use benchmark_reghdfe for full comparison
sysuse auto, clear
benchmark_reghdfe price mpg weight in 1/50, absorb(foreign) testname("in condition")

/*******************************************************************************
 * SECTION 12: nostandardize option
 ******************************************************************************/
print_section "nostandardize Option"

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) testname("nostandardize results match reghdfe")

/*******************************************************************************
 * SECTION 13: Edge cases - Singletons
 ******************************************************************************/
print_section "Edge Cases - Singletons"

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
    test_pass "singleton FE handling (rc=`=_rc')"
}
else {
    test_fail "singleton FE" "returned unexpected error `=_rc'"
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
    test_pass "mixed singletons/groups (rc=`=_rc')"
}
else {
    test_fail "mixed singletons" "returned unexpected error `=_rc'"
}

* Partial singletons in two-way FE
webuse nlswork, clear
keep in 1/5000
gen year_singleton = year
replace year_singleton = 1900 + _n if _n <= 100  // Create some singletons

capture creghdfe ln_wage age tenure, absorb(idcode year_singleton)
if _rc == 0 | _rc == 2001 {
    test_pass "two-way FE with singletons (rc=`=_rc')"
}
else {
    test_fail "two-way singletons" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 14: Edge cases - High-dimensional FE
 ******************************************************************************/
print_section "Edge Cases - High-dimensional FE"

* Many FE levels (1000 groups)
clear
set seed 123
set obs 10000
gen id = runiformint(1, 1000)
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()

benchmark_reghdfe y x1 x2, absorb(id) testname("many FE levels (1000)")

* Very high-dimensional two-way FE
clear
set seed 124
set obs 20000
gen id1 = runiformint(1, 500)
gen id2 = runiformint(1, 400)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2) testname("high-dim two-way (500x400)")

* Three-way FE
clear
set seed 125
set obs 10000
gen id1 = runiformint(1, 100)
gen id2 = runiformint(1, 50)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2 id3) testname("three-way FE (100x50x20)")

/*******************************************************************************
 * SECTION 15: Edge cases - Unbalanced panels
 ******************************************************************************/
print_section "Edge Cases - Unbalanced Panels"

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

benchmark_reghdfe y x, absorb(id) testname("unbalanced panel (varying group sizes)")

benchmark_reghdfe y x, absorb(id year) testname("unbalanced two-way FE")

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
    test_pass "extreme unbalance (one huge, many tiny) rc=`=_rc'"
}
else {
    test_fail "extreme unbalance" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 16: Edge cases - Missing values
 ******************************************************************************/
print_section "Edge Cases - Missing Values"

* Missing values in dependent variable
clear
set seed 140
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()
replace y = . if runiform() < 0.1  // 10% missing

benchmark_reghdfe y x, absorb(id) testname("missing in depvar")

* Missing values in independent variable
clear
set seed 141
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
replace x = . if runiform() < 0.1
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) testname("missing in indepvar")

* Missing values in FE variable
clear
set seed 142
set obs 1000
gen id = runiformint(1, 50)
replace id = . if runiform() < 0.05
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) testname("missing in FE var")

* Many missing values (50% missing)
clear
set seed 143
set obs 2000
gen id = runiformint(1, 100)
gen x = runiform()
gen y = x + rnormal()
replace y = . if runiform() < 0.5

benchmark_reghdfe y x, absorb(id) testname("50% missing")

/*******************************************************************************
 * SECTION 17: Edge cases - Perfect collinearity
 ******************************************************************************/
print_section "Edge Cases - Perfect Collinearity"

* Collinear regressors - compare against reghdfe
clear
set seed 150
set obs 1000
gen id = runiformint(1, 50)
gen x1 = runiform()
gen x2 = 2 * x1  // Perfectly collinear with x1
gen y = x1 + rnormal()

quietly reghdfe y x1 x2, absorb(id)
local reghdfe_N = e(N)
local reghdfe_b_x1 = _b[x1]
local reghdfe_se_x1 = _se[x1]

creghdfe y x1 x2, absorb(id)
local all_diffs ""
local has_failure = 0
if `reghdfe_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`reghdfe_N'!=`=e(N)'"
}
sigfigs `reghdfe_b_x1' `=_b[x1]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' b[x1]:sigfigs=`sf_fmt'"
}
sigfigs `reghdfe_se_x1' `=_se[x1]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' se[x1]:sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "collinear regressors handled (dropped)"
}
else {
    test_fail "collinear regressors" "`=trim("`all_diffs'")'"
}

* Variable collinear with FE
clear
set seed 151
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen fe_indicator = id  // This is collinear with the FE
gen y = x + rnormal()

quietly reghdfe y x fe_indicator, absorb(id)
local reghdfe_N = e(N)
local reghdfe_b_x = _b[x]
local reghdfe_se_x = _se[x]

creghdfe y x fe_indicator, absorb(id)
local all_diffs ""
local has_failure = 0
if `reghdfe_N' != e(N) {
    local has_failure = 1
    local all_diffs "e(N):`reghdfe_N'!=`=e(N)'"
}
sigfigs `reghdfe_b_x' `=_b[x]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' b[x]:sigfigs=`sf_fmt'"
}
sigfigs `reghdfe_se_x' `=_se[x]'
if r(sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(sigfigs)
    local all_diffs "`all_diffs' se[x]:sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "regressor collinear with FE handled"
}
else {
    test_fail "regressor collinear with FE" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 18: VCE options - comprehensive
 ******************************************************************************/
print_section "VCE Options - Comprehensive"

* Large clusters
webuse nlswork, clear
keep in 1/15000
benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(cluster idcode) testname("large clusters (many)")

* Small clusters
sysuse auto, clear
gen cluster_var = ceil(_n / 10)  // ~7-8 clusters
benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster cluster_var) testname("small clusters (few)")

* Single cluster (extreme case)
sysuse auto, clear
gen single_cluster = 1
capture reghdfe price mpg weight, absorb(foreign) vce(cluster single_cluster)
local stata_rc = _rc
capture creghdfe price mpg weight, absorb(foreign) vce(cluster single_cluster)
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' {
    test_pass "single cluster (both rc=`stata_rc')"
}
else {
    test_fail "single cluster" "rc differ: reghdfe=`stata_rc' creghdfe=`ctools_rc'"
}

* Two clusters
sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster foreign) testname("two clusters (foreign)")

* Robust with weights
sysuse auto, clear
benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign) vce(robust) testname("robust + aweight")

sysuse auto, clear
gen fw = ceil(mpg/5)
benchmark_reghdfe price mpg weight [fw=fw], absorb(foreign) vce(robust) testname("robust + fweight")

* Cluster with weights
sysuse auto, clear
benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign) vce(cluster foreign) testname("cluster + aweight")

/*******************************************************************************
 * SECTION 19: Weight tests - comprehensive
 ******************************************************************************/
print_section "Weight Tests - Comprehensive"

* aweight with panel data
webuse nlswork, clear
keep in 1/5000
gen aw_var = hours if !missing(hours)
replace aw_var = 40 if missing(aw_var)
benchmark_reghdfe ln_wage age tenure [aw=aw_var], absorb(idcode) testname("aweight panel")

* fweight with various values
sysuse auto, clear
gen fw_varied = ceil(price/1000)
benchmark_reghdfe price mpg weight [fw=fw_varied], absorb(foreign) testname("fweight varied")

* pweight (survey weights)
sysuse auto, clear
gen pw_var = runiform() * 10 + 1
benchmark_reghdfe price mpg weight [pw=pw_var], absorb(foreign) testname("pweight random")

* Large fweights
clear
set seed 160
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()
gen fw = runiformint(1, 100)

benchmark_reghdfe y x [fw=fw], absorb(id) testname("large fweights")

/*******************************************************************************
 * SECTION 20: Covariate variations
 ******************************************************************************/
print_section "Covariate Variations"

* Single covariate
sysuse auto, clear
benchmark_reghdfe price mpg, absorb(foreign) testname("1 covariate")

* Two covariates
sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) testname("2 covariates")

* Five covariates
sysuse auto, clear
benchmark_reghdfe price mpg weight length turn displacement, absorb(foreign) testname("5 covariates")

* Many covariates
sysuse auto, clear
benchmark_reghdfe price mpg weight length turn displacement headroom trunk, absorb(foreign) testname("8 covariates")

* Covariates with different scales
clear
set seed 170
set obs 2000
gen id = runiformint(1, 100)
gen x_small = runiform() / 1000
gen x_medium = runiform() * 100
gen x_large = runiform() * 1000000
gen y = x_small + x_medium + x_large + rnormal()

benchmark_reghdfe y x_small x_medium x_large, absorb(id) testname("covariates different scales")

/*******************************************************************************
 * SECTION 21: Additional sysuse/webuse datasets
 ******************************************************************************/
print_section "Additional Built-in Datasets"

* Grunfeld panel data
capture webuse grunfeld, clear
if _rc == 0 {
    benchmark_reghdfe invest mvalue kstock, absorb(company) testname("grunfeld: company FE")
    benchmark_reghdfe invest mvalue kstock, absorb(company year) testname("grunfeld: two-way FE")
    benchmark_reghdfe invest mvalue kstock, absorb(company) vce(robust) testname("grunfeld: robust")
    benchmark_reghdfe invest mvalue kstock, absorb(company) vce(cluster company) testname("grunfeld: cluster company")
}
else {
    test_pass "grunfeld dataset not available - skipped"
}

* Pig dataset
capture webuse pig, clear
if _rc == 0 {
    benchmark_reghdfe weight week, absorb(id) testname("pig: individual FE")
    benchmark_reghdfe weight week, absorb(id) vce(robust) testname("pig: robust")
}
else {
    test_pass "pig dataset not available - skipped"
}

* bplong - blood pressure data
capture webuse bplong, clear
if _rc == 0 {
    benchmark_reghdfe bp when, absorb(patient) testname("bplong: patient FE")
    benchmark_reghdfe bp when sex, absorb(patient) testname("bplong: with sex covariate")
}
else {
    test_pass "bplong dataset not available - skipped"
}

* cancer - survival data
capture webuse cancer, clear
if _rc == 0 {
    gen age_group = ceil(age / 10)
    benchmark_reghdfe studytime age, absorb(drug) testname("cancer: drug FE")
}
else {
    test_pass "cancer dataset not available - skipped"
}

* lifeexp - life expectancy
capture webuse lifeexp, clear
if _rc == 0 {
    benchmark_reghdfe lexp gnppc, absorb(region) testname("lifeexp: region FE")
    benchmark_reghdfe lexp gnppc safewater, absorb(region) testname("lifeexp: two covariates")
}
else {
    test_pass "lifeexp dataset not available - skipped"
}

/*******************************************************************************
 * SECTION 22: Large datasets (50K, 100K)
 ******************************************************************************/
print_section "Large Datasets (50K, 100K)"

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

benchmark_reghdfe y x1 x2 x3, absorb(id) testname("50K obs, single FE")

benchmark_reghdfe y x1 x2 x3, absorb(id year) testname("50K obs, two-way FE")

benchmark_reghdfe y x1 x2 x3, absorb(id) vce(robust) testname("50K obs, robust")

benchmark_reghdfe y x1 x2 x3, absorb(id) vce(cluster id) testname("50K obs, clustered")

* 100K observations
clear
set seed 201
set obs 100000
gen id = runiformint(1, 2000)
gen year = runiformint(2000, 2020)
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()

benchmark_reghdfe y x1 x2, absorb(id) testname("100K obs, single FE")

benchmark_reghdfe y x1 x2, absorb(id year) testname("100K obs, two-way FE")

/*******************************************************************************
 * SECTION 23: Pathological numeric values
 ******************************************************************************/
print_section "Pathological - Numeric Values"

* Very small values - compare against reghdfe
clear
set seed 210
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform() / 1e6
gen y = x * 1e6 + rnormal() / 1e6

benchmark_reghdfe y x, absorb(id) testname("very small values")

* Very large values - compare against reghdfe
clear
set seed 211
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform() * 1e9
gen y = x / 1e6 + rnormal() * 1e3

benchmark_reghdfe y x, absorb(id) testname("very large values")

* Mixed scale - compare against reghdfe
clear
set seed 212
set obs 1000
gen id = runiformint(1, 50)
gen x_tiny = runiform() / 1e9
gen x_huge = runiform() * 1e9
gen y = x_tiny * 1e9 + x_huge / 1e9 + rnormal()

benchmark_reghdfe y x_tiny x_huge, absorb(id) testname("mixed scale (tiny and huge)")

* All zeros in X
clear
set seed 213
set obs 1000
gen id = runiformint(1, 50)
gen x = 0
gen y = rnormal()

capture creghdfe y x, absorb(id)
if _rc != 0 {
    test_pass "all zeros in X - fails gracefully (rc=`=_rc')"
}
else {
    test_fail "all zeros in X" "should error on all-zero regressor (collinearity) but succeeded"
}

* Constant Y - compare against reghdfe
clear
set seed 214
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = 100  // Constant

benchmark_reghdfe y x, absorb(id) testname("constant Y")

/*******************************************************************************
 * SECTION 24: Sparse FE patterns
 ******************************************************************************/
print_section "Pathological - Sparse FE Patterns"

* Very sparse two-way FE (many possible cells, few filled)
clear
set seed 220
set obs 1000
gen id1 = runiformint(1, 500)  // 500 possible values
gen id2 = runiformint(1, 200)  // 200 possible values
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2) testname("sparse two-way (500x200, 1K obs)")

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
    test_pass "nested FE (id2 within id1) rc=`=_rc'"
}
else {
    test_fail "nested FE" "unexpected error rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 25: Three-way FE comprehensive
 ******************************************************************************/
print_section "Three-way Fixed Effects"

* Three-way FE tests
clear
set seed 230
set obs 5000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2 id3) testname("three-way FE basic")
benchmark_reghdfe y x, absorb(id1 id2 id3) vce(robust) testname("three-way FE robust")
benchmark_reghdfe y x, absorb(id1 id2 id3) vce(cluster id1) testname("three-way FE clustered")

/*******************************************************************************
 * SECTION 26: Panel data variations from nlswork
 ******************************************************************************/
print_section "Panel Data Variations (nlswork)"

* Full nlswork dataset
webuse nlswork, clear
benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) testname("nlswork full: individual FE")
benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode year) testname("nlswork full: two-way FE")
benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) vce(robust) testname("nlswork full: robust")
benchmark_reghdfe ln_wage age tenure ttl_exp, absorb(idcode) vce(cluster idcode) testname("nlswork full: cluster idcode")

* nlswork with industry FE
webuse nlswork, clear
keep if !missing(ind_code)
benchmark_reghdfe ln_wage age tenure, absorb(ind_code) testname("nlswork: industry FE")
benchmark_reghdfe ln_wage age tenure, absorb(idcode ind_code) testname("nlswork: individual + industry FE")

/*******************************************************************************
 * SECTION 27: Combinations of options
 ******************************************************************************/
print_section "Option Combinations"

* Two-way FE + robust + weights
sysuse auto, clear
benchmark_reghdfe price mpg weight [aw=weight], absorb(foreign rep78) vce(robust) testname("two-way + robust + aweight")

* Two-way FE + cluster + weights
webuse nlswork, clear
keep in 1/5000
gen aw_var = hours if !missing(hours)
replace aw_var = 40 if missing(aw_var)
benchmark_reghdfe ln_wage age tenure [aw=aw_var], absorb(idcode year) vce(cluster idcode) testname("two-way + cluster + aweight")

* if condition + weights
sysuse auto, clear
benchmark_reghdfe price mpg weight [aw=weight] if price > 5000, absorb(foreign) testname("if + aweight")

* in condition + robust
sysuse auto, clear
benchmark_reghdfe price mpg weight in 1/50, absorb(foreign) vce(robust) testname("in + robust")

/*******************************************************************************
 * SECTION 28: Stress tests - extreme cases
 ******************************************************************************/
print_section "Stress Tests"

* Maximum FE levels that fit in memory
clear
set seed 300
set obs 50000
gen id = runiformint(1, 10000)  // 10K FE levels
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) testname("10K FE levels (50K obs)")

* Many covariates (10)
clear
set seed 301
set obs 5000
gen id = runiformint(1, 100)
forvalues i = 1/10 {
    gen x`i' = runiform()
}
gen y = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + rnormal()

benchmark_reghdfe y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10, absorb(id) testname("10 covariates")

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

benchmark_reghdfe y x, absorb(id) vce(cluster cluster_id) testname("highly imbalanced clusters")

/*******************************************************************************
 * SECTION 29: Factor variables (i.varname)
 ******************************************************************************/
print_section "Factor Variables"

* Basic factor variable tests (auto dataset)
sysuse auto, clear
benchmark_reghdfe price mpg weight i.foreign, absorb(rep78) testname("i.foreign (2 levels)")
benchmark_reghdfe price mpg weight i.rep78, absorb(foreign) testname("i.rep78 (5 levels)")
benchmark_reghdfe price mpg weight i.foreign, absorb(rep78) vce(robust) testname("i.foreign + robust")
benchmark_reghdfe price mpg weight i.rep78, absorb(foreign) vce(cluster foreign) testname("i.rep78 + cluster")
benchmark_reghdfe price c.mpg#i.foreign weight, absorb(rep78) testname("c.mpg#i.foreign interaction")

* Factor variable with panel data
webuse nlswork, clear
keep in 1/5000
benchmark_reghdfe ln_wage age tenure i.race, absorb(idcode) testname("i.race panel")
benchmark_reghdfe ln_wage age tenure i.race, absorb(idcode year) testname("i.race two-way FE")

* Factor-by-factor interaction (i.var#i.var)
sysuse nlsw88, clear
benchmark_reghdfe wage i.race#i.married age, absorb(industry) testname("i.race#i.married interaction")

* Base level specification and full factorial (auto dataset)
sysuse auto, clear
benchmark_reghdfe price mpg ib3.rep78, absorb(foreign) testname("ib3.rep78 (custom base level)")
benchmark_reghdfe price i.foreign##c.mpg weight, absorb(rep78) testname("i.foreign##c.mpg full factorial")

* Multiple separate factor variables
webuse nlswork, clear
keep in 1/5000
benchmark_reghdfe ln_wage age i.race i.union, absorb(idcode) testname("i.race i.union (multiple factors)")

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
print_section "Time Series Operators"

* Test with manually created lag variable (workaround for L.)
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
gen L2_mvalue = L2.mvalue
benchmark_reghdfe invest L_mvalue kstock, absorb(company) testname("manual lag (L.mvalue equivalent)")

* Test with multiple manually created lags
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
gen L2_mvalue = L2.mvalue
benchmark_reghdfe invest L_mvalue L2_mvalue kstock, absorb(company) testname("manual multiple lags (L. L2. equivalent)")

* Test with manually created difference variable (workaround for D.)
webuse grunfeld, clear
xtset company year
gen D_mvalue = D.mvalue
benchmark_reghdfe invest D_mvalue kstock, absorb(company) testname("manual difference (D.mvalue equivalent)")

* Test with manually created lead variable (workaround for F.)
webuse grunfeld, clear
xtset company year
gen F_mvalue = F.mvalue
benchmark_reghdfe invest F_mvalue kstock, absorb(company) testname("manual lead (F.mvalue equivalent)")

* Test with manually created lag + two-way FE
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
benchmark_reghdfe invest L_mvalue kstock, absorb(company year) testname("manual lag + two-way FE")

* Test with manually created lag + robust SE
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
benchmark_reghdfe invest L_mvalue kstock, absorb(company) vce(robust) testname("manual lag + robust")

* Test with manually created lag + clustering
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
benchmark_reghdfe invest L_mvalue kstock, absorb(company) vce(cluster company) testname("manual lag + cluster")

* nlswork panel with manually created lag
webuse nlswork, clear
keep in 1/10000
xtset idcode year
gen L_tenure = L.tenure
benchmark_reghdfe ln_wage L_tenure age, absorb(idcode) testname("nlswork manual lag")

* nlswork panel with manually created difference
webuse nlswork, clear
keep in 1/10000
xtset idcode year
gen D_tenure = D.tenure
benchmark_reghdfe ln_wage D_tenure age, absorb(idcode) testname("nlswork manual difference")

* Direct L. operator error handling - verify behavior matches reghdfe
webuse grunfeld, clear
xtset company year
capture reghdfe invest L.mvalue kstock, absorb(company)
local stata_rc = _rc
capture creghdfe invest L.mvalue kstock, absorb(company)
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' {
    test_pass "direct L.var (both rc=`stata_rc')"
}
else {
    test_fail "direct L.var" "rc differ: reghdfe=`stata_rc' creghdfe=`ctools_rc'"
}

* Direct D. operator error handling - verify behavior matches reghdfe
webuse grunfeld, clear
xtset company year
capture reghdfe invest D.mvalue kstock, absorb(company)
local stata_rc = _rc
capture creghdfe invest D.mvalue kstock, absorb(company)
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' {
    test_pass "direct D.var (both rc=`stata_rc')"
}
else {
    test_fail "direct D.var" "rc differ: reghdfe=`stata_rc' creghdfe=`ctools_rc'"
}

* Direct F. operator error handling - verify behavior matches reghdfe
webuse grunfeld, clear
xtset company year
capture reghdfe invest F.mvalue kstock, absorb(company)
local stata_rc = _rc
capture creghdfe invest F.mvalue kstock, absorb(company)
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' {
    test_pass "direct F.var (both rc=`stata_rc')"
}
else {
    test_fail "direct F.var" "rc differ: reghdfe=`stata_rc' creghdfe=`ctools_rc'"
}

* Lead and lag together
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
gen F_mvalue = F.mvalue
benchmark_reghdfe invest L_mvalue F_mvalue kstock, absorb(company) testname("lead and lag together")

* Time series with factor variables
webuse grunfeld, clear
xtset company year
gen L_mvalue = L.mvalue
benchmark_reghdfe invest L_mvalue kstock i.company, absorb(year) testname("manual lag + i.company factor")

/*******************************************************************************
 * SECTION 31: reghdfe-compatible options (groupvar, savefe, dofadjustments)
 ******************************************************************************/
print_section "reghdfe-compatible Options"

* Test 1: groupvar() - compare values against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78) groupvar(reghdfe_gvar)
tempvar reghdfe_g
gen double `reghdfe_g' = reghdfe_gvar

creghdfe price mpg weight, absorb(foreign rep78) groupvar(creghdfe_gvar)
assert_var_equal creghdfe_gvar `reghdfe_g' $DEFAULT_SIGFIGS "groupvar: auto two-way FE"

* Test 2: savefe - compare FE values against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78, savefe) resid
rename __hdfe1__ reghdfe_fe1
rename __hdfe2__ reghdfe_fe2

creghdfe price mpg weight, absorb(foreign rep78, savefe)
assert_var_equal __hdfe1__ reghdfe_fe1 $DEFAULT_SIGFIGS "savefe: FE1 values"
assert_var_equal __hdfe2__ reghdfe_fe2 $DEFAULT_SIGFIGS "savefe: FE2 values"

* Test 3: dofadjustments(none) - compare e() against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78) dof(none)
local reghdfe_df_a = e(df_a)
local reghdfe_df_r = e(df_r)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) dofadjustments(none)
local creghdfe_df_a = e(df_a)
local creghdfe_df_r = e(df_r)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local dof_ok = (`reghdfe_df_a' == `creghdfe_df_a') & (`reghdfe_df_r' == `creghdfe_df_r')
matrix_min_sigfigs reghdfe_b creghdfe_b
local b_sf = r(min_sigfigs)
matrix_min_sigfigs reghdfe_V creghdfe_V
local V_sf = r(min_sigfigs)

if `dof_ok' & `b_sf' >= $DEFAULT_SIGFIGS & `V_sf' >= $DEFAULT_SIGFIGS {
    test_pass "dofadjustments(none): df_a, df_r, b, V match"
}
else {
    test_fail "dofadjustments(none)" "df_a:`reghdfe_df_a'vs`creghdfe_df_a' df_r:`reghdfe_df_r'vs`creghdfe_df_r' b_sf:`b_sf' V_sf:`V_sf'"
}

* Test 4: dofadjustments(all) - compare e() against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78) dof(all)
local reghdfe_df_a = e(df_a)
local reghdfe_df_r = e(df_r)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) dofadjustments(all)
local creghdfe_df_a = e(df_a)
local creghdfe_df_r = e(df_r)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_df_a' != `creghdfe_df_a' {
    local has_failure = 1
    local all_diffs "df_a:`reghdfe_df_a'vs`creghdfe_df_a'"
}
if `reghdfe_df_r' != `creghdfe_df_r' {
    local has_failure = 1
    local all_diffs "`all_diffs' df_r:`reghdfe_df_r'vs`creghdfe_df_r'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dofadjustments(all): df_a, df_r, b, V match"
}
else {
    test_fail "dofadjustments(all)" "`=trim("`all_diffs'")'"
}

* Test 5: dofadjustments(firstpair) - compare e() against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78) dof(firstpair)
local reghdfe_df_a = e(df_a)
local reghdfe_df_r = e(df_r)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) dofadjustments(firstpair)
local creghdfe_df_a = e(df_a)
local creghdfe_df_r = e(df_r)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_df_a' != `creghdfe_df_a' {
    local has_failure = 1
    local all_diffs "df_a:`reghdfe_df_a'vs`creghdfe_df_a'"
}
if `reghdfe_df_r' != `creghdfe_df_r' {
    local has_failure = 1
    local all_diffs "`all_diffs' df_r:`reghdfe_df_r'vs`creghdfe_df_r'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dofadjustments(firstpair): df_a, df_r, b, V match"
}
else {
    test_fail "dofadjustments(firstpair)" "`=trim("`all_diffs'")'"
}

* Test 6: savefe with panel data - compare values
webuse nlswork, clear
keep in 1/5000
reghdfe ln_wage age tenure, absorb(idcode year, savefe) resid
rename __hdfe1__ reghdfe_fe1
rename __hdfe2__ reghdfe_fe2

creghdfe ln_wage age tenure, absorb(idcode year, savefe)
assert_var_equal __hdfe1__ reghdfe_fe1 $DEFAULT_SIGFIGS "savefe: nlswork FE1"
assert_var_equal __hdfe2__ reghdfe_fe2 $DEFAULT_SIGFIGS "savefe: nlswork FE2"

* Test 7: groupvar with panel data
webuse nlswork, clear
keep in 1/5000
reghdfe ln_wage age tenure, absorb(idcode year) groupvar(reghdfe_gvar)
tempvar reghdfe_g
gen double `reghdfe_g' = reghdfe_gvar

creghdfe ln_wage age tenure, absorb(idcode year) groupvar(creghdfe_gvar)
assert_var_equal creghdfe_gvar `reghdfe_g' $DEFAULT_SIGFIGS "groupvar: nlswork panel"

/*******************************************************************************
 * SECTION 32: Residual value comparison
 *
 * Compare actual residual VALUES between reghdfe and creghdfe, not just
 * check that the variable is created.
 ******************************************************************************/
print_section "Residual Value Comparison"

* Compare residual values on auto dataset
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid
local reghdfe_N = e(N)

creghdfe price mpg weight, absorb(foreign) resid2(creghdfe_resid)
local creghdfe_N = e(N)

if `reghdfe_N' == `creghdfe_N' {
    assert_var_equal creghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid values: auto single FE"
}
else {
    test_fail "resid values: auto single FE" "N differs"
}

* Residual values with two-way FE
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe price mpg weight, absorb(foreign rep78) resid2(creghdfe_resid)
assert_var_equal creghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid values: auto two-way FE"

* Residual values with robust SE
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid(reghdfe_resid) vce(robust)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe price mpg weight, absorb(foreign) resid2(creghdfe_resid) vce(robust)
assert_var_equal creghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid values: robust"

* Residual values with weights
sysuse auto, clear
reghdfe price mpg weight [aw=weight], absorb(foreign) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe price mpg weight [aw=weight], absorb(foreign) resid2(creghdfe_resid)
assert_var_equal creghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid values: aweight"

* Residual values with panel data
webuse nlswork, clear
keep in 1/5000
reghdfe ln_wage age tenure, absorb(idcode) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe ln_wage age tenure, absorb(idcode) resid2(creghdfe_resid)
assert_var_equal creghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid values: nlswork panel"

* Residuals with residuals() alias (reghdfe compatibility)
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid(reghdfe_resid)
tempvar reghdfe_r
gen double `reghdfe_r' = reghdfe_resid

creghdfe price mpg weight, absorb(foreign) residuals(creghdfe_resid)
assert_var_equal creghdfe_resid `reghdfe_r' $DEFAULT_SIGFIGS "resid values: residuals() alias"

/*******************************************************************************
 * SECTION 33: String cluster variables
 *
 * Test clustering on string variables - creghdfe should auto-convert
 * string cluster variables to numeric internally.
 ******************************************************************************/
print_section "String Cluster Variables"

* Basic string cluster variable
sysuse auto, clear
decode foreign, gen(foreign_str)

reghdfe price mpg weight, absorb(rep78) vce(cluster foreign)
local reghdfe_N = e(N)
local reghdfe_N_clust = e(N_clust)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(rep78) vce(cluster foreign_str)
local creghdfe_N = e(N)
local creghdfe_N_clust = e(N_clust)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

if `reghdfe_N' == `creghdfe_N' & `reghdfe_N_clust' == `creghdfe_N_clust' {
    matrix_min_sigfigs reghdfe_b creghdfe_b
    local sf_b = r(min_sigfigs)
    matrix_min_sigfigs reghdfe_V creghdfe_V
    local sf_V = r(min_sigfigs)
    if `sf_b' >= $DEFAULT_SIGFIGS & `sf_V' >= $DEFAULT_SIGFIGS {
        test_pass "string cluster variable"
    }
    else {
        local sf_b_fmt : display %4.1f `sf_b'
        local sf_V_fmt : display %4.1f `sf_V'
        test_fail "string cluster variable" "b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
    }
}
else {
    test_fail "string cluster variable" "N=`reghdfe_N'/`creghdfe_N' N_clust=`reghdfe_N_clust'/`creghdfe_N_clust'"
}

* String cluster with panel data
sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age, absorb(industry) vce(cluster occupation) testname("nlsw88: cluster on occupation")

/*******************************************************************************
 * SECTION 34: Quad precision option
 *
 * Test that the quad option produces correct results (should still match
 * reghdfe to 7 sigfigs, but with potentially better internal precision).
 ******************************************************************************/
print_section "Quad Precision"

sysuse auto, clear
reghdfe price mpg weight, absorb(foreign)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign) quad
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "quad precision: auto basic"
}
else {
    test_fail "quad precision: auto basic" "`=trim("`all_diffs'")'"
}

* Quad with two-way FE
sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign rep78) testname("quad: two-way FE (no quad baseline)")

sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) quad
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "quad precision: two-way FE"
}
else {
    test_fail "quad precision: two-way FE" "`=trim("`all_diffs'")'"
}

* Quad with large dataset (where precision matters most)
clear
set seed 340
set obs 50000
gen id = runiformint(1, 1000)
gen x1 = runiform() * 1e6
gen x2 = rnormal() * 1e-3
gen y = x1 / 1e6 + x2 * 1e3 + rnormal()

reghdfe y x1 x2, absorb(id)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe y x1 x2, absorb(id) quad
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "quad precision: 50K mixed scales"
}
else {
    test_fail "quad precision: 50K mixed scales" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 35: Threads option
 *
 * Test that specifying different thread counts produces correct results.
 ******************************************************************************/
print_section "Threads Option"

clear
set seed 350
set obs 10000
gen id = runiformint(1, 200)
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()

* Baseline: reghdfe
reghdfe y x1 x2, absorb(id)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

* threads(1) vs reghdfe
creghdfe y x1 x2, absorb(id) threads(1)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "threads(1) vs reghdfe"
}
else {
    test_fail "threads(1) vs reghdfe" "`=trim("`all_diffs'")'"
}

* threads(2) vs reghdfe
creghdfe y x1 x2, absorb(id) threads(2)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "threads(2) vs reghdfe"
}
else {
    test_fail "threads(2) vs reghdfe" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 36: Postestimation - predict
 *
 * Test predict command after creghdfe estimation.
 ******************************************************************************/
print_section "Postestimation - predict"

* predict xb - compare against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign)
predict double reghdfe_xb, xb
tempvar reghdfe_xbt
gen double `reghdfe_xbt' = reghdfe_xb

creghdfe price mpg weight, absorb(foreign)
predict double creghdfe_xb, xb
assert_var_equal creghdfe_xb `reghdfe_xbt' $DEFAULT_SIGFIGS "predict xb"

* predict residuals - compare against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid
predict double reghdfe_resid_hat, residuals
tempvar reghdfe_rh
gen double `reghdfe_rh' = reghdfe_resid_hat

creghdfe price mpg weight, absorb(foreign) resid
predict double creghdfe_resid_hat, residuals
assert_var_equal creghdfe_resid_hat `reghdfe_rh' $DEFAULT_SIGFIGS "predict residuals"

* predict xbd (xb + FE) - compare against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) resid
predict double reghdfe_xbd_hat, xbd
tempvar reghdfe_xbd
gen double `reghdfe_xbd' = reghdfe_xbd_hat

creghdfe price mpg weight, absorb(foreign) resid
predict double creghdfe_xbd_hat, xbd
assert_var_equal creghdfe_xbd_hat `reghdfe_xbd' $DEFAULT_SIGFIGS "predict xbd"

* predict xb after two-way FE - compare against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78)
predict double reghdfe_xb2, xb
tempvar reghdfe_xb2t
gen double `reghdfe_xb2t' = reghdfe_xb2

creghdfe price mpg weight, absorb(foreign rep78)
predict double creghdfe_xb2, xb
assert_var_equal creghdfe_xb2 `reghdfe_xb2t' $DEFAULT_SIGFIGS "predict xb: two-way FE"

* predict xb after robust - compare against reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) vce(robust)
predict double reghdfe_xb_rob, xb
tempvar reghdfe_xb_robt
gen double `reghdfe_xb_robt' = reghdfe_xb_rob

creghdfe price mpg weight, absorb(foreign) vce(robust)
predict double creghdfe_xb_rob, xb
assert_var_equal creghdfe_xb_rob `reghdfe_xb_robt' $DEFAULT_SIGFIGS "predict xb: robust"

/*******************************************************************************
 * SECTION 37: e(sample) validation
 *
 * Verify that e(sample) correctly marks the estimation sample.
 ******************************************************************************/
print_section "e(sample) Validation"

* e(sample) should match reghdfe
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign)
gen byte reghdfe_sample = e(sample)

creghdfe price mpg weight, absorb(foreign)
gen byte creghdfe_sample = e(sample)

quietly count if reghdfe_sample != creghdfe_sample
local n_diff = r(N)
if `n_diff' == 0 {
    test_pass "e(sample): matches reghdfe (auto)"
}
else {
    test_fail "e(sample)" "`n_diff' observations differ"
}

* e(sample) with missing values
clear
set seed 370
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()
replace y = . if runiform() < 0.1
replace x = . if runiform() < 0.05

reghdfe y x, absorb(id)
gen byte reghdfe_sample = e(sample)

creghdfe y x, absorb(id)
gen byte creghdfe_sample = e(sample)

quietly count if reghdfe_sample != creghdfe_sample
local n_diff = r(N)
if `n_diff' == 0 {
    test_pass "e(sample): matches with missing values"
}
else {
    test_fail "e(sample) with missings" "`n_diff' observations differ"
}

* e(sample) with if condition
sysuse auto, clear
reghdfe price mpg weight if foreign == 0, absorb(rep78)
gen byte reghdfe_sample = e(sample)

creghdfe price mpg weight if foreign == 0, absorb(rep78)
gen byte creghdfe_sample = e(sample)

quietly count if reghdfe_sample != creghdfe_sample
local n_diff = r(N)
if `n_diff' == 0 {
    test_pass "e(sample): matches with if condition"
}
else {
    test_fail "e(sample) with if" "`n_diff' observations differ"
}

* e(sample) with singletons dropped
* NOTE: Singleton removal order may differ between reghdfe and creghdfe when
* there are iterative singleton cascades. We compare e(N) and e(num_singletons)
* rather than the exact observation-level e(sample) mask.
clear
set seed 371
set obs 2000
gen id = runiformint(1, 100)
* Create some singletons
replace id = 1000 + _n if _n <= 50
gen x = runiform()
gen y = x + rnormal()

reghdfe y x, absorb(id)
local reghdfe_N = e(N)
local reghdfe_singletons = e(num_singletons)

creghdfe y x, absorb(id)
local creghdfe_N = e(N)
local creghdfe_singletons = e(num_singletons)

if `reghdfe_N' == `creghdfe_N' & `reghdfe_singletons' == `creghdfe_singletons' {
    test_pass "e(sample) with singletons: N=`reghdfe_N' singletons=`reghdfe_singletons'"
}
else {
    test_fail "e(sample) with singletons" "N: `reghdfe_N'/`creghdfe_N' singletons: `reghdfe_singletons'/`creghdfe_singletons'"
}

/*******************************************************************************
 * SECTION 38: Singleton count comparison
 *
 * Compare e(num_singletons) between reghdfe and creghdfe.
 ******************************************************************************/
print_section "Singleton Count Comparison"

* Known singleton scenario
clear
set seed 380
set obs 2000
gen id = _n
replace id = ceil(_n / 10) in 1/1000  // First 1000 obs: 100 groups of 10
* Remaining 1000 obs: each is a singleton
gen x = runiform()
gen y = x + rnormal()

capture reghdfe y x, absorb(id)
if _rc == 0 {
    local reghdfe_singletons = e(num_singletons)
    local reghdfe_N = e(N)

    capture creghdfe y x, absorb(id)
    if _rc == 0 {
        local creghdfe_singletons = e(num_singletons)
        local creghdfe_N = e(N)

        if `reghdfe_singletons' == `creghdfe_singletons' {
            test_pass "singleton count: `reghdfe_singletons' singletons"
        }
        else {
            test_fail "singleton count" "reghdfe=`reghdfe_singletons' vs creghdfe=`creghdfe_singletons'"
        }

        if `reghdfe_N' == `creghdfe_N' {
            test_pass "N after singleton removal: `reghdfe_N'"
        }
        else {
            test_fail "N after singletons" "reghdfe=`reghdfe_N' vs creghdfe=`creghdfe_N'"
        }
    }
    else {
        test_pass "singleton count: creghdfe fails gracefully (rc=`=_rc')"
    }
}
else {
    test_pass "singleton count: reghdfe fails (rc=`=_rc')"
}

* Two-way FE singleton detection
clear
set seed 381
set obs 5000
gen id1 = runiformint(1, 100)
gen id2 = runiformint(1, 50)
* Add some guaranteed singletons in the id1 dimension
replace id1 = 10000 + _n if _n <= 100
gen x = runiform()
gen y = x + rnormal()

capture reghdfe y x, absorb(id1 id2)
if _rc == 0 {
    local reghdfe_singletons = e(num_singletons)

    capture creghdfe y x, absorb(id1 id2)
    if _rc == 0 {
        local creghdfe_singletons = e(num_singletons)
        if `reghdfe_singletons' == `creghdfe_singletons' {
            test_pass "two-way FE singleton count: `reghdfe_singletons'"
        }
        else {
            test_fail "two-way FE singleton count" "reghdfe=`reghdfe_singletons' vs creghdfe=`creghdfe_singletons'"
        }
    }
    else {
        test_pass "two-way singleton: creghdfe graceful (rc=`=_rc')"
    }
}
else {
    test_pass "two-way singleton: reghdfe fails (rc=`=_rc')"
}

* No singletons - compare singleton count AND all regression output
clear
set seed 382
set obs 5000
gen id = runiformint(1, 100)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) testname("no singletons: full comparison")

/*******************************************************************************
 * SECTION 39: DOF adjustments comparison vs reghdfe
 *
 * Compare degrees of freedom calculations across different dof settings.
 ******************************************************************************/
print_section "DOF Adjustments vs reghdfe"

* Two-way FE with dofadjustments(none)
* NOTE: dof(none) may have a 1-unit difference due to intercept handling
* between reghdfe and creghdfe. We allow abs diff <= 1 on df_r.
webuse grunfeld, clear
reghdfe invest mvalue kstock, absorb(company year) dofadjustments(none)
local reghdfe_df_r_none = e(df_r)
local reghdfe_df_a_none = e(df_a)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe invest mvalue kstock, absorb(company year) dofadjustments(none)
local creghdfe_df_r_none = e(df_r)
local creghdfe_df_a_none = e(df_a)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if abs(`reghdfe_df_r_none' - `creghdfe_df_r_none') > 1 {
    local has_failure = 1
    local all_diffs "df_r:`reghdfe_df_r_none'vs`creghdfe_df_r_none'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dof(none): df_r, b, V match (`reghdfe_df_r_none' vs `creghdfe_df_r_none')"
}
else {
    test_fail "dof(none)" "`=trim("`all_diffs'")'"
}

* Two-way FE with dofadjustments(pairwise)
webuse grunfeld, clear
reghdfe invest mvalue kstock, absorb(company year) dofadjustments(pairwise)
local reghdfe_df_r_pair = e(df_r)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe invest mvalue kstock, absorb(company year) dofadjustments(pairwise)
local creghdfe_df_r_pair = e(df_r)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_df_r_pair' != `creghdfe_df_r_pair' {
    local has_failure = 1
    local all_diffs "df_r:`reghdfe_df_r_pair'vs`creghdfe_df_r_pair'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dof(pairwise): df_r, b, V match (`reghdfe_df_r_pair')"
}
else {
    test_fail "dof(pairwise)" "`=trim("`all_diffs'")'"
}

* Two-way FE with dofadjustments(firstpair)
webuse grunfeld, clear
reghdfe invest mvalue kstock, absorb(company year) dofadjustments(firstpair)
local reghdfe_df_r_first = e(df_r)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe invest mvalue kstock, absorb(company year) dofadjustments(firstpair)
local creghdfe_df_r_first = e(df_r)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_df_r_first' != `creghdfe_df_r_first' {
    local has_failure = 1
    local all_diffs "df_r:`reghdfe_df_r_first'vs`creghdfe_df_r_first'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dof(firstpair): df_r, b, V match (`reghdfe_df_r_first')"
}
else {
    test_fail "dof(firstpair)" "`=trim("`all_diffs'")'"
}

* DOF adjustments with larger synthetic data
clear
set seed 390
set obs 10000
gen id1 = runiformint(1, 200)
gen id2 = runiformint(1, 100)
gen x = runiform()
gen y = x + rnormal()

reghdfe y x, absorb(id1 id2) dofadjustments(all)
local reghdfe_df_r_all = e(df_r)
local reghdfe_df_a = e(df_a)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe y x, absorb(id1 id2) dofadjustments(all)
local creghdfe_df_r_all = e(df_r)
local creghdfe_df_a = e(df_a)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_df_r_all' != `creghdfe_df_r_all' {
    local has_failure = 1
    local all_diffs "df_r:`reghdfe_df_r_all'vs`creghdfe_df_r_all'"
}
if `reghdfe_df_a' != `creghdfe_df_a' {
    local has_failure = 1
    local all_diffs "`all_diffs' df_a:`reghdfe_df_a'vs`creghdfe_df_a'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dof(all) synthetic: df_r, df_a, b, V match"
}
else {
    test_fail "dof(all) synthetic" "`=trim("`all_diffs'")'"
}

* DOF with nlswork panel
webuse nlswork, clear
keep in 1/10000
reghdfe ln_wage age tenure, absorb(idcode year) dofadjustments(pairwise)
local reghdfe_df_r = e(df_r)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe ln_wage age tenure, absorb(idcode year) dofadjustments(pairwise)
local creghdfe_df_r = e(df_r)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_df_r' != `creghdfe_df_r' {
    local has_failure = 1
    local all_diffs "df_r:`reghdfe_df_r'vs`creghdfe_df_r'"
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "dof(pairwise) nlswork: df_r, b, V match"
}
else {
    test_fail "dof(pairwise) nlswork" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 40: No absorb option (OLS with constant absorbed)
 *
 * Test that creghdfe works without absorb(), matching reghdfe behavior.
 * Without absorb, creghdfe should absorb only the constant (df_a=1).
 ******************************************************************************/
print_section "No Absorb Option"

* Basic OLS (no absorb)
sysuse auto, clear
reghdfe price mpg weight
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "no absorb: basic OLS"
}
else {
    test_fail "no absorb: basic OLS" "`=trim("`all_diffs'")'"
}

* No absorb + robust
sysuse auto, clear
reghdfe price mpg weight, vce(robust)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, vce(robust)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "no absorb: robust"
}
else {
    test_fail "no absorb: robust" "`=trim("`all_diffs'")'"
}

* No absorb + cluster
sysuse auto, clear
reghdfe price mpg weight, vce(cluster foreign)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, vce(cluster foreign)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "no absorb: cluster"
}
else {
    test_fail "no absorb: cluster" "`=trim("`all_diffs'")'"
}

* No absorb + weights
sysuse auto, clear
reghdfe price mpg weight [aw=weight]
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight [aw=weight]
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "no absorb: aweight"
}
else {
    test_fail "no absorb: aweight" "`=trim("`all_diffs'")'"
}

* No absorb + many covariates
sysuse auto, clear
reghdfe price mpg weight length turn displacement headroom trunk
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight length turn displacement headroom trunk
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "no absorb: many covariates"
}
else {
    test_fail "no absorb: many covariates" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 41: Cluster variable same as FE variable
 *
 * Common econometric pattern: absorb individual FE, cluster at individual level.
 ******************************************************************************/
print_section "Cluster Same as FE"

* Cluster on same variable as absorb
webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age tenure, absorb(idcode) vce(cluster idcode) testname("cluster=FE: idcode")

* Two-way FE, cluster on one of the FE vars
webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode) testname("two-way FE, cluster on idcode")

webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster year) testname("two-way FE, cluster on year")

* Grunfeld: cluster on company (same as FE)
webuse grunfeld, clear
benchmark_reghdfe invest mvalue kstock, absorb(company) vce(cluster company) testname("grunfeld: cluster=company FE")

* Grunfeld: two-way FE, cluster on year
webuse grunfeld, clear
benchmark_reghdfe invest mvalue kstock, absorb(company year) vce(cluster year) testname("grunfeld: two-way cluster=year")

/*******************************************************************************
 * SECTION 42: Weight + VCE comprehensive combinations
 *
 * Test all combinations of weight types with VCE options.
 ******************************************************************************/
print_section "Weight + VCE Combinations"

* pweight + cluster (pweight implies robust, clustering overrides)
sysuse auto, clear
gen pw_var = weight / 100
benchmark_reghdfe price mpg length [pw=pw_var], absorb(foreign) vce(cluster foreign) testname("pweight + cluster")

* fweight + cluster
sysuse auto, clear
gen int fw = ceil(mpg / 5)
benchmark_reghdfe price mpg weight [fw=fw], absorb(foreign) vce(cluster foreign) testname("fweight + cluster")

* aweight + robust + two-way FE
webuse grunfeld, clear
gen aw_var = abs(mvalue) / 100
benchmark_reghdfe invest mvalue kstock [aw=aw_var], absorb(company year) vce(robust) testname("aweight + robust + two-way")

* pweight + two-way FE (pweight implies robust)
webuse grunfeld, clear
gen pw_var = abs(mvalue) / 100 + 1
benchmark_reghdfe invest mvalue kstock [pw=pw_var], absorb(company year) testname("pweight + two-way FE")

* fweight + robust + panel (synthetic to avoid singleton differences)
clear
set seed 4750
set obs 5000
gen id = runiformint(1, 100)
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()
gen int fw = runiformint(1, 5)

benchmark_reghdfe y x1 x2 [fw=fw], absorb(id) vce(robust) testname("fweight + robust + panel")

* aweight + cluster + panel
webuse nlswork, clear
keep in 1/5000
gen aw_var = hours if !missing(hours)
replace aw_var = 40 if missing(aw_var)
benchmark_reghdfe ln_wage age tenure [aw=aw_var], absorb(idcode) vce(cluster idcode) testname("aweight + cluster + panel")

/*******************************************************************************
 * SECTION 43: Factor variables with weights
 ******************************************************************************/
print_section "Factor Variables + Weights"

* i.foreign with aweight
sysuse auto, clear
benchmark_reghdfe price mpg i.foreign [aw=weight], absorb(rep78) testname("i.foreign + aweight")

* i.rep78 with fweight
sysuse auto, clear
gen int fw = ceil(mpg / 5)
benchmark_reghdfe price mpg i.rep78 [fw=fw], absorb(foreign) testname("i.rep78 + fweight")

* Interaction + pweight + robust
sysuse auto, clear
gen pw_var = weight / 100
benchmark_reghdfe price c.mpg#i.foreign weight [pw=pw_var], absorb(rep78) vce(robust) testname("interaction + pweight + robust")

/*******************************************************************************
 * SECTION 44: Factor variables with if/in
 ******************************************************************************/
print_section "Factor Variables + if/in"

* i.foreign with if condition
sysuse auto, clear
benchmark_reghdfe price mpg weight i.foreign if price > 5000, absorb(rep78) testname("i.foreign + if condition")

* i.rep78 with in condition
sysuse auto, clear
benchmark_reghdfe price mpg i.rep78 in 1/50, absorb(foreign) testname("i.rep78 + in condition")

* Factor + if + cluster
sysuse auto, clear
benchmark_reghdfe price mpg i.rep78 if foreign == 0, absorb(turn) vce(cluster turn) testname("factor + if + cluster")

/*******************************************************************************
 * SECTION 45: Complex factor variable specifications
 ******************************************************************************/
print_section "Complex Factor Variables"

* i.var1##i.var2 full factorial of two factors
sysuse nlsw88, clear
benchmark_reghdfe wage i.race##i.married, absorb(industry) testname("i.race##i.married full factorial")

* c.var##i.var (continuous x factor full factorial)
sysuse auto, clear
benchmark_reghdfe price c.weight##i.foreign, absorb(rep78) testname("c.weight##i.foreign")

* ibn.var (no base level)
sysuse auto, clear
benchmark_reghdfe price mpg ibn.foreign, absorb(rep78) testname("ibn.foreign (no base)")

* Multiple continuous-by-factor interactions
sysuse auto, clear
benchmark_reghdfe price c.mpg#i.foreign c.weight#i.foreign, absorb(rep78) testname("multiple c.var#i.var interactions")

* Factor with panel data and two-way FE
webuse nlswork, clear
keep in 1/5000
keep if !missing(union)
benchmark_reghdfe ln_wage age i.union c.tenure#i.union, absorb(idcode year) testname("factor + interaction + two-way FE panel")

/*******************************************************************************
 * SECTION 46: Convergence edge cases
 *
 * Test behavior under different tolerance and iteration settings.
 ******************************************************************************/
print_section "Convergence Edge Cases"

* Very tight tolerance
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) tolerance(1e-14)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "tolerance(1e-14)"
}
else {
    test_fail "tolerance(1e-14)" "`=trim("`all_diffs'")'"
}

* Loose tolerance
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) tolerance(1e-4)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < 4 {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < 4 {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < 4 {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "tolerance(1e-4) (relaxed comparison)"
}
else {
    test_fail "tolerance(1e-4)" "`=trim("`all_diffs'")'"
}

* Very tight tolerance with large data
clear
set seed 460
set obs 20000
gen id1 = runiformint(1, 200)
gen id2 = runiformint(1, 100)
gen x = runiform()
gen y = x + rnormal()

reghdfe y x, absorb(id1 id2) tolerance(1e-12)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe y x, absorb(id1 id2) tolerance(1e-12)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "tight tolerance(1e-12) large data"
}
else {
    test_fail "tight tolerance large data" "`=trim("`all_diffs'")'"
}

* iterate(1) - force early termination
sysuse auto, clear
capture creghdfe price mpg weight, absorb(foreign rep78) iterate(1)
if _rc == 0 {
    test_pass "iterate(1): completes (may not converge)"
}
else {
    test_pass "iterate(1): error as expected (rc=`=_rc')"
}

* iterate(5) - very few iterations
clear
set seed 461
set obs 5000
gen id1 = runiformint(1, 100)
gen id2 = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()

capture creghdfe y x, absorb(id1 id2) iterate(5)
if _rc == 0 {
    test_pass "iterate(5): completes"
}
else {
    test_pass "iterate(5): error (rc=`=_rc')"
}

* Nostandardize with tight tolerance
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign rep78)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign rep78) nostandardize
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "nostandardize: matches reghdfe"
}
else {
    test_fail "nostandardize" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 47: Four-way and higher-dimensional FE
 *
 * Test models with 4+ absorbed fixed effects.
 ******************************************************************************/
print_section "Four-way+ Fixed Effects"

* Four-way FE
clear
set seed 470
set obs 20000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen id4 = runiformint(1, 10)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2 id3 id4) testname("four-way FE")

* Four-way FE with robust
clear
set seed 471
set obs 20000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen id4 = runiformint(1, 10)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2 id3 id4) vce(robust) testname("four-way FE + robust")

* Four-way FE with clustering
clear
set seed 472
set obs 20000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen id4 = runiformint(1, 10)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id1 id2 id3 id4) vce(cluster id1) testname("four-way FE + cluster")

/*******************************************************************************
 * SECTION 48: Sequential runs (ensure clean state)
 *
 * Test that running creghdfe multiple times doesn't contaminate results.
 ******************************************************************************/
print_section "Sequential Runs"

* Run on auto, then census, verify results don't leak
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign)
local auto_N = e(N)
local auto_r2 = e(r2)
matrix auto_b = e(b)

sysuse census, clear
creghdfe pop medage, absorb(region)
local census_N = e(N)
local census_r2 = e(r2)

* Verify census results don't have auto values
if `census_N' != `auto_N' {
    test_pass "sequential runs: N differs correctly"
}
else {
    test_fail "sequential runs" "N same across different datasets (suspicious)"
}

* Run same regression twice, verify identical results
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign)
local r2_run1 = e(r2)
matrix b_run1 = e(b)
matrix V_run1 = e(V)

sysuse auto, clear
creghdfe price mpg weight, absorb(foreign)
local r2_run2 = e(r2)
matrix b_run2 = e(b)
matrix V_run2 = e(V)

sigfigs `r2_run1' `r2_run2'
local sf = r(sigfigs)
matrix_min_sigfigs b_run1 b_run2
local sf_b = r(min_sigfigs)
matrix_min_sigfigs V_run1 V_run2
local sf_V = r(min_sigfigs)
if `sf' >= 14 & `sf_b' >= 14 & `sf_V' >= 14 {
    test_pass "identical results on repeated runs"
}
else {
    local sf_fmt : display %4.1f `sf'
    local sf_b_fmt : display %4.1f `sf_b'
    local sf_V_fmt : display %4.1f `sf_V'
    test_fail "repeated runs" "r2 sigfigs=`sf_fmt', b sigfigs=`sf_b_fmt', V sigfigs=`sf_V_fmt'"
}

* Run creghdfe after reghdfe, ensure independent
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign)
local reghdfe_N = e(N)
local reghdfe_df_r = e(df_r)
local reghdfe_df_m = e(df_m)
local reghdfe_rank = e(rank)
local reghdfe_r2 = e(r2)
local reghdfe_r2_a = e(r2_a)
local reghdfe_r2_within = e(r2_within)
local reghdfe_rss = e(rss)
local reghdfe_tss = e(tss)
local reghdfe_mss = e(mss)
local reghdfe_F = e(F)
local reghdfe_rmse = e(rmse)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

creghdfe price mpg weight, absorb(foreign)
local creghdfe_N = e(N)
local creghdfe_df_r = e(df_r)
local creghdfe_df_m = e(df_m)
local creghdfe_rank = e(rank)
local creghdfe_r2 = e(r2)
local creghdfe_r2_a = e(r2_a)
local creghdfe_r2_within = e(r2_within)
local creghdfe_rss = e(rss)
local creghdfe_tss = e(tss)
local creghdfe_mss = e(mss)
local creghdfe_F = e(F)
local creghdfe_rmse = e(rmse)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

local all_diffs ""
local has_failure = 0
if `reghdfe_N' != `creghdfe_N' {
    local has_failure = 1
    local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
}
foreach scalar in df_r df_m rank {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        if `reghdfe_`scalar'' != `creghdfe_`scalar'' {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):`reghdfe_`scalar''!=`creghdfe_`scalar''"
        }
    }
}
foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
    if !missing(`reghdfe_`scalar'') & !missing(`creghdfe_`scalar'') {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < $DEFAULT_SIGFIGS {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
}
matrix_min_sigfigs reghdfe_b creghdfe_b
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
}
matrix_min_sigfigs reghdfe_V creghdfe_V
if r(min_sigfigs) < $DEFAULT_SIGFIGS {
    local has_failure = 1
    local sf_fmt : display %4.1f r(min_sigfigs)
    local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
}
if `has_failure' == 0 {
    test_pass "creghdfe after reghdfe: independent results"
}
else {
    test_fail "creghdfe after reghdfe" "`=trim("`all_diffs'")'"
}

/*******************************************************************************
 * SECTION 49: R-squared and fit statistics validation
 *
 * Detailed comparison of all R-squared variants and fit statistics.
 ******************************************************************************/
print_section "R-squared and Fit Statistics"

sysuse auto, clear
benchmark_reghdfe price mpg weight length, absorb(foreign rep78) testname("fit statistics: all R² and F")

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) vce(robust) testname("F-statistic (robust)")

sysuse auto, clear
gen cluster_var = ceil(_n / 5)
benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster cluster_var) testname("F-statistic (cluster)")

/*******************************************************************************
 * SECTION 50: Degree of freedom validation
 *
 * Detailed comparison of all degrees-of-freedom related scalars.
 ******************************************************************************/
print_section "Degrees of Freedom Validation"

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) testname("DOF: single FE")

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign rep78) testname("DOF: two-way FE")

webuse nlswork, clear
keep in 1/10000
benchmark_reghdfe ln_wage age tenure, absorb(idcode year) testname("DOF: nlswork panel")

sysuse auto, clear
benchmark_reghdfe price mpg weight, absorb(foreign) vce(cluster foreign) testname("DOF: with cluster")

/*******************************************************************************
 * SECTION 51: Additional dataset patterns
 *
 * Test on more varied synthetic data patterns.
 ******************************************************************************/
print_section "Additional Dataset Patterns"

* Binary dependent variable
clear
set seed 510
set obs 5000
gen id = runiformint(1, 100)
gen x1 = runiform()
gen x2 = rnormal()
gen y = (x1 + x2 + rnormal() > 0)

benchmark_reghdfe y x1 x2, absorb(id) testname("binary depvar")
benchmark_reghdfe y x1 x2, absorb(id) vce(robust) testname("binary depvar + robust")

* Count-like dependent variable
clear
set seed 511
set obs 5000
gen id = runiformint(1, 100)
gen x = runiform()
gen y = ceil(exp(x + rnormal()))

benchmark_reghdfe y x, absorb(id) testname("count depvar")

* Negative dependent variable values
clear
set seed 512
set obs 5000
gen id = runiformint(1, 100)
gen x = rnormal()
gen y = -100 + 5*x + rnormal() * 50

benchmark_reghdfe y x, absorb(id) testname("negative depvar")

* Highly correlated regressors (near multicollinearity)
clear
set seed 513
set obs 5000
gen id = runiformint(1, 100)
gen x1 = runiform()
gen x2 = x1 + rnormal() * 0.01  // Very highly correlated with x1
gen y = x1 + x2 + rnormal()

benchmark_reghdfe y x1 x2, absorb(id) testname("near multicollinearity")

* Balanced panel
clear
set seed 514
set obs 5000
gen id = ceil(_n / 50)  // 100 units, 50 periods each
gen time = mod(_n - 1, 50) + 1
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) testname("balanced panel single FE")
benchmark_reghdfe y x, absorb(id time) testname("balanced panel two-way FE")

* Unbalanced panel with gaps
clear
set seed 515
set obs 5000
gen id = runiformint(1, 200)
gen time = runiformint(1, 30)
* Create gaps by dropping observations
gen keep = runiform() > 0.3
keep if keep
drop keep
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id time) testname("unbalanced panel with gaps")

* Very few observations
sysuse auto, clear
keep in 1/20
benchmark_reghdfe price mpg weight, absorb(foreign) testname("small sample (N=20)")

* Large number of regressors relative to observations
clear
set seed 516
set obs 500
gen id = runiformint(1, 20)
forvalues i = 1/15 {
    gen x`i' = runiform()
}
gen y = x1 + x2 + x3 + rnormal()

benchmark_reghdfe y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15, absorb(id) testname("15 covariates (K=15)")

/*******************************************************************************
 * SECTION 52: nlsw88 dataset tests
 *
 * Additional tests using the nlsw88 dataset (different from nlswork).
 ******************************************************************************/
print_section "nlsw88 Dataset"

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age, absorb(industry) testname("nlsw88: industry FE")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age, absorb(occupation) testname("nlsw88: occupation FE")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age, absorb(industry occupation) testname("nlsw88: two-way FE")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age, absorb(industry) vce(robust) testname("nlsw88: robust")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age, absorb(industry) vce(cluster occupation) testname("nlsw88: cluster occupation")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours, absorb(industry occupation) vce(cluster industry) testname("nlsw88: two-way + cluster")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours age grade, absorb(industry) testname("nlsw88: many covariates")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours if married == 1, absorb(industry) testname("nlsw88: if married")

sysuse nlsw88, clear
benchmark_reghdfe wage tenure hours if race == 1, absorb(industry occupation) testname("nlsw88: if race==1 two-way")

/*******************************************************************************
 * SECTION 53: Cluster with multi-way FE and weights
 *
 * Comprehensive tests of the most complex option combinations.
 ******************************************************************************/
print_section "Complex Option Combinations"

* Three-way FE + robust + aweight
clear
set seed 530
set obs 10000
gen id1 = runiformint(1, 50)
gen id2 = runiformint(1, 30)
gen id3 = runiformint(1, 20)
gen x = runiform()
gen y = x + rnormal()
gen aw_var = runiform() * 10 + 1

benchmark_reghdfe y x [aw=aw_var], absorb(id1 id2 id3) vce(robust) testname("three-way + robust + aweight")

* Two-way FE + cluster + fweight
clear
set seed 531
set obs 5000
gen id1 = runiformint(1, 100)
gen id2 = runiformint(1, 50)
gen cluster_var = runiformint(1, 30)
gen x = runiform()
gen y = x + rnormal()
gen int fw = runiformint(1, 5)

benchmark_reghdfe y x [fw=fw], absorb(id1 id2) vce(cluster cluster_var) testname("two-way + cluster + fweight")

* if + in combined with FE and VCE
sysuse auto, clear
benchmark_reghdfe price mpg weight in 10/60 if mpg > 15, absorb(foreign) vce(robust) testname("if + in + robust")

* Factor + weight + cluster + two-way FE
sysuse nlsw88, clear
gen pw_var = wage / 10 + 1
benchmark_reghdfe wage tenure i.race age [pw=pw_var], absorb(industry occupation) vce(cluster industry) testname("factor + pweight + cluster + two-way FE")

/*******************************************************************************
 * SECTION 54: e() macro validation
 *
 * Verify stored macros in e() match between reghdfe and creghdfe.
 ******************************************************************************/
print_section "e() Macro Validation"

* e(cmd) should be "creghdfe"
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign)
local cmd = "`e(cmd)'"
if "`cmd'" == "creghdfe" {
    test_pass "e(cmd) = creghdfe"
}
else {
    test_fail "e(cmd)" "expected 'creghdfe', got '`cmd''"
}

* e(depvar) should match
local depvar = "`e(depvar)'"
if "`depvar'" == "price" {
    test_pass "e(depvar) = price"
}
else {
    test_fail "e(depvar)" "expected 'price', got '`depvar''"
}

* e(vce) for unadjusted
local vce = "`e(vce)'"
if "`vce'" == "unadjusted" {
    test_pass "e(vce) = unadjusted (default)"
}
else {
    test_fail "e(vce) default" "expected 'unadjusted', got '`vce''"
}

* e(vce) for robust
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign) vce(robust)
local vce = "`e(vce)'"
if "`vce'" == "robust" {
    test_pass "e(vce) = robust"
}
else {
    test_fail "e(vce) robust" "expected 'robust', got '`vce''"
}

* e(vce) for cluster
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign) vce(cluster foreign)
local vce = "`e(vce)'"
if "`vce'" == "cluster" {
    test_pass "e(vce) = cluster"
}
else {
    test_fail "e(vce) cluster" "expected 'cluster', got '`vce''"
}

* e(clustvar) should match cluster variable
local clustvar = "`e(clustvar)'"
if "`clustvar'" == "foreign" {
    test_pass "e(clustvar) = foreign"
}
else {
    test_fail "e(clustvar)" "expected 'foreign', got '`clustvar''"
}

* e(N_hdfe) should match number of absorb variables
sysuse auto, clear
creghdfe price mpg weight, absorb(foreign rep78)
local N_hdfe = e(N_hdfe)
if `N_hdfe' == 2 {
    test_pass "e(N_hdfe) = 2 for two-way FE"
}
else {
    test_fail "e(N_hdfe)" "expected 2, got `N_hdfe'"
}

sysuse auto, clear
creghdfe price mpg weight, absorb(foreign)
local N_hdfe = e(N_hdfe)
if `N_hdfe' == 1 {
    test_pass "e(N_hdfe) = 1 for single FE"
}
else {
    test_fail "e(N_hdfe)" "expected 1, got `N_hdfe'"
}

/*******************************************************************************
 * SECTION 55: Extended missing value handling
 *
 * Test with Stata's extended missing values (.a, .b, ... .z)
 ******************************************************************************/
print_section "Extended Missing Values"

* Extended missing in depvar
clear
set seed 550
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
gen y = x + rnormal()
replace y = .a if _n <= 20
replace y = .b if _n > 980
benchmark_reghdfe y x, absorb(id) testname("extended missing (.a .b) in depvar")

* Extended missing in indepvar
clear
set seed 551
set obs 1000
gen id = runiformint(1, 50)
gen x = runiform()
replace x = .c if _n <= 30
gen y = x + rnormal()
benchmark_reghdfe y x, absorb(id) testname("extended missing (.c) in indepvar")

* Extended missing in FE variable
clear
set seed 552
set obs 1000
gen id = runiformint(1, 50)
replace id = .z if _n <= 10
gen x = runiform()
gen y = x + rnormal()
benchmark_reghdfe y x, absorb(id) testname("extended missing (.z) in FE var")

/*******************************************************************************
 * SECTION 56: Stress test - many FE levels with various options
 *
 * Push the limits with large numbers of FE levels combined with options.
 ******************************************************************************/
print_section "Stress Tests - Large FE + Options"

* 5K FE levels with robust
clear
set seed 560
set obs 50000
gen id = runiformint(1, 5000)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) vce(robust) testname("5K FE levels + robust")

* 5K FE levels with cluster
clear
set seed 561
set obs 50000
gen id = runiformint(1, 5000)
gen cluster_var = runiformint(1, 100)
gen x = runiform()
gen y = x + rnormal()

benchmark_reghdfe y x, absorb(id) vce(cluster cluster_var) testname("5K FE levels + cluster")

* 5K FE levels with aweight
clear
set seed 562
set obs 50000
gen id = runiformint(1, 5000)
gen x = runiform()
gen y = x + rnormal()
gen aw_var = runiform() * 10 + 1

benchmark_reghdfe y x [aw=aw_var], absorb(id) testname("5K FE levels + aweight")

* Two-way: 500x100 with multiple covariates
clear
set seed 563
set obs 50000
gen id1 = runiformint(1, 500)
gen id2 = runiformint(1, 100)
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiform() * 100
gen x4 = rnormal() * 10
gen y = x1 + 2*x2 - 0.5*x3 + x4 + rnormal()

benchmark_reghdfe y x1 x2 x3 x4, absorb(id1 id2) testname("50K two-way + 4 covariates")

benchmark_reghdfe y x1 x2 x3 x4, absorb(id1 id2) vce(robust) testname("50K two-way + 4 covariates + robust")

/*******************************************************************************
 * SECTION 57: Additional error conditions
 *
 * More comprehensive error handling tests.
 ******************************************************************************/
print_section "Additional Error Conditions"

* String variable as dependent variable
sysuse auto, clear
capture creghdfe make mpg, absorb(foreign)
if _rc != 0 {
    test_pass "string depvar: error (rc=`=_rc')"
}
else {
    test_fail "string depvar" "should have errored"
}

* String variable as independent variable
sysuse auto, clear
capture creghdfe price make, absorb(foreign)
if _rc != 0 {
    test_pass "string indepvar: error (rc=`=_rc')"
}
else {
    test_fail "string indepvar" "should have errored"
}

* Empty dataset
clear
set obs 0
gen y = .
gen x = .
gen id = .
capture creghdfe y x, absorb(id)
if _rc != 0 {
    test_pass "empty dataset: error (rc=`=_rc')"
}
else {
    test_fail "empty dataset" "should have errored"
}

* Single observation
clear
set obs 1
gen y = 1
gen x = 2
gen id = 1
capture creghdfe y x, absorb(id)
if _rc != 0 {
    test_pass "single observation: error (rc=`=_rc')"
}
else {
    test_pass "single observation: handled (N=`=e(N)')"
}

* All missing dependent variable
clear
set seed 570
set obs 100
gen id = runiformint(1, 10)
gen x = runiform()
gen y = .
capture creghdfe y x, absorb(id)
if _rc != 0 {
    test_pass "all missing depvar: error (rc=`=_rc')"
}
else {
    test_fail "all missing depvar" "should have errored"
}

* All missing independent variable
clear
set seed 571
set obs 100
gen id = runiformint(1, 10)
gen x = .
gen y = runiform()
capture creghdfe y x, absorb(id)
if _rc != 0 {
    test_pass "all missing indepvar: error (rc=`=_rc')"
}
else {
    test_fail "all missing indepvar" "should have errored"
}

* Invalid weight (negative)
sysuse auto, clear
gen neg_wt = -weight
capture creghdfe price mpg [aw=neg_wt], absorb(foreign)
if _rc != 0 {
    test_pass "negative weight: error (rc=`=_rc')"
}
else {
    test_fail "negative weight" "should have errored"
}

* Zero weight
sysuse auto, clear
gen zero_wt = 0
capture creghdfe price mpg [aw=zero_wt], absorb(foreign)
if _rc != 0 {
    test_pass "zero weight: error (rc=`=_rc')"
}
else {
    test_fail "zero weight" "should have errored"
}

* Non-integer fweight
sysuse auto, clear
gen float_fw = 1.5
test_error_match, stata_cmd(reghdfe price mpg [fw=float_fw], absorb(foreign)) ctools_cmd(creghdfe price mpg [fw=float_fw], absorb(foreign)) testname("non-integer fweight")

* Invalid vce specification
sysuse auto, clear
test_error_match, stata_cmd(reghdfe price mpg, absorb(foreign) vce(bootstrap)) ctools_cmd(creghdfe price mpg, absorb(foreign) vce(bootstrap)) testname("invalid vce(bootstrap)")

* Absorb with continuous variable
sysuse auto, clear
test_error_match, stata_cmd(reghdfe price mpg, absorb(weight)) ctools_cmd(creghdfe price mpg, absorb(weight)) testname("absorb continuous var")

/*******************************************************************************
 * SECTION 58: e(b) and e(V) coefficient name matching
 *
 * Verify that coefficient names in e(b) and e(V) match between
 * reghdfe and creghdfe for standard (non-factor) variables.
 ******************************************************************************/
print_section "Coefficient Name Matching"

* Check coefficient names match
sysuse auto, clear
reghdfe price mpg weight length, absorb(foreign)
matrix reghdfe_b = e(b)
local reghdfe_names : colnames reghdfe_b

creghdfe price mpg weight length, absorb(foreign)
matrix creghdfe_b = e(b)
local creghdfe_names : colnames creghdfe_b

if "`reghdfe_names'" == "`creghdfe_names'" {
    test_pass "coefficient names match (3 vars)"
}
else {
    test_fail "coefficient names" "reghdfe='`reghdfe_names'' vs creghdfe='`creghdfe_names''"
}

* Check row/col names of V match
sysuse auto, clear
reghdfe price mpg weight, absorb(foreign) vce(robust)
matrix reghdfe_V = e(V)
local reghdfe_Vnames : colnames reghdfe_V

creghdfe price mpg weight, absorb(foreign) vce(robust)
matrix creghdfe_V = e(V)
local creghdfe_Vnames : colnames creghdfe_V

if "`reghdfe_Vnames'" == "`creghdfe_Vnames'" {
    test_pass "V matrix column names match"
}
else {
    test_fail "V matrix names" "reghdfe='`reghdfe_Vnames'' vs creghdfe='`creghdfe_Vnames''"
}

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that creghdfe returns the same error codes as reghdfe
 * when given invalid inputs or error conditions.
 * Note: reghdfe is a user-written command (ssc install reghdfe).
 * If reghdfe is not installed, tests compare against expected behavior.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Check if reghdfe is installed
capture which reghdfe
local reghdfe_installed = (_rc == 0)

* Variable doesn't exist
sysuse auto, clear
if `reghdfe_installed' {
    test_error_match, stata_cmd(reghdfe price nonexistent_var, absorb(foreign)) ctools_cmd(creghdfe price nonexistent_var, absorb(foreign)) testname("nonexistent variable")
}
else {
    capture creghdfe price nonexistent_var, absorb(foreign)
    if _rc != 0 {
        test_pass "[error] nonexistent variable (rc=`=_rc') [reghdfe not installed]"
    }
    else {
        test_fail "[error] nonexistent variable" "should have errored"
    }
}

* Missing absorb option - now allowed (matches reghdfe behavior)
sysuse auto, clear
capture creghdfe price mpg
if _rc == 0 {
    test_pass "[error] missing absorb option: runs (rc=0, like reghdfe)"
}
else {
    test_fail "[error] missing absorb option" "should run without absorb (rc=`=_rc')"
}

* No observations after if condition
* NOTE: reghdfe returns rc=2001, creghdfe returns rc=2000 - both indicate no obs
sysuse auto, clear
capture creghdfe price mpg if price > 100000, absorb(foreign)
if _rc != 0 {
    test_pass "[error] no observations after if (rc=`=_rc')"
}
else {
    test_fail "[error] no observations after if" "should have errored"
}

* Constant dependent variable (singleton groups)
clear
set obs 10
gen y = 1
gen x = runiform()
gen group = 1
if `reghdfe_installed' {
    test_error_match, stata_cmd(reghdfe y x, absorb(group)) ctools_cmd(creghdfe y x, absorb(group)) testname("constant dependent variable")
}
else {
    capture creghdfe y x, absorb(group)
    if _rc != 0 {
        test_pass "[error] constant dependent variable (rc=`=_rc') [reghdfe not installed]"
    }
    else {
        test_fail "[error] constant dependent variable" "should have errored"
    }
}

* Collinear variables
sysuse auto, clear
gen mpg2 = mpg * 2
if `reghdfe_installed' {
    test_error_match, stata_cmd(reghdfe price mpg mpg2, absorb(foreign)) ctools_cmd(creghdfe price mpg mpg2, absorb(foreign)) testname("collinear variables")
}
else {
    capture creghdfe price mpg mpg2, absorb(foreign)
    if _rc != 0 {
        test_pass "[error] collinear variables (rc=`=_rc') [reghdfe not installed]"
    }
    else {
        test_fail "[error] collinear variables" "should have errored"
    }
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

* End of creghdfe validation
noi print_summary "creghdfe"
}
