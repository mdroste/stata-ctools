* cppmlhdfe: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cppmlhdfe"

if "`__case'" == "cppmlhdfe" {
* validate_cppmlhdfe.do
* Comprehensive validation tests comparing cppmlhdfe vs ppmlhdfe
* Tests: basic Poisson, multi-way FE, weights, cluster VCE, separation,
* offset, exposure, factor variables, if/in, missing data,
* collinearity, large datasets, pathological values, error handling
* VERIFICATION: All tests compare e() scalars, e(b), and e(V) between
* ppmlhdfe and cppmlhdfe using the benchmark_ppmlhdfe helper

do validation/support/validate_setup.do

di as text ""
di as text "===== Validation: cppmlhdfe vs ppmlhdfe ====="
di as text ""

* Check ppmlhdfe is installed
capture which ppmlhdfe
if _rc != 0 {
    di as error "ppmlhdfe not installed. Install with: ssc install ppmlhdfe"
    di as error "Also need: ssc install ftools; ssc install reghdfe"
    exit 601
}

* Use the pinned official ships fixture on every run.
capture program drop load_ships_data
program define load_ships_data
    ctools_fixture ships, clear
end

quietly {

* Plugin check
noi di as text "--- Plugin check ---"
load_ships_data
capture cppmlhdfe accident op_75_79, absorb(ship) vce(robust)
if _rc != 0 {
    test_fail "cppmlhdfe plugin load" "returned error `=_rc'"
    exit 1
}
test_pass "cppmlhdfe plugin loads and runs"

* SECTION 1: Basic FE (ships data)
noi di as text "--- Section 1: Basic Fixed Effects (ships) ---"

load_ships_data

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) vce(robust) testname("ships: single FE robust")

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) testname("ships: single FE unadjusted")

benchmark_ppmlhdfe accident op_75_79, ///
    absorb(ship) vce(robust) testname("ships: single covariate")

benchmark_ppmlhdfe accident co_65_69 co_70_74 co_75_79, ///
    absorb(ship) vce(robust) testname("ships: three covariates")

* SECTION 2: Exposure and offset
noi di as text "--- Section 2: Exposure and Offset ---"

load_ships_data

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) exposure(service) vce(robust) testname("ships: exposure")

* Offset with synthetic data
clear
set obs 1000
set seed 50
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen offset_var = rnormal() * 0.5
gen mu = exp(0.3 * x + offset_var + 0.05 * fe1)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) offset(offset_var) vce(robust) testname("synth: offset")

* Exposure with synthetic data
clear
set obs 1000
set seed 51
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen exposure_var = 1 + runiform() * 10
gen mu = exp(0.3 * x + 0.05 * fe1) * exposure_var
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) exposure(exposure_var) vce(robust) testname("synth: exposure")

* SECTION 3: Two-way FE
noi di as text "--- Section 3: Two-way Fixed Effects ---"

clear
set obs 2000
set seed 99
gen fe1 = mod(_n-1, 20) + 1
gen fe2 = mod(_n-1, 10) + 1
gen x = rnormal()
gen mu = exp(0.3 * x + 0.05 * fe1 + 0.1 * fe2)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) vce(robust) testname("synth: two-way FE robust")

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) testname("synth: two-way FE unadjusted")

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) vce(cluster fe1) testname("synth: two-way FE cluster")

* SECTION 4: Three-way FE
noi di as text "--- Section 4: Three-way Fixed Effects ---"

clear
set obs 5000
set seed 100
gen fe1 = mod(_n-1, 20) + 1
gen fe2 = mod(_n-1, 10) + 1
gen fe3 = mod(_n-1, 5) + 1
gen x = rnormal()
gen mu = exp(0.3 * x + 0.05 * fe1 + 0.1 * fe2 + 0.15 * fe3)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2 fe3) vce(robust) testname("synth: three-way FE robust")

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2 fe3) vce(cluster fe1) testname("synth: three-way FE cluster")

* SECTION 5: VCE options
noi di as text "--- Section 5: VCE Options ---"

load_ships_data

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) vce(robust) testname("ships: vce(robust)")

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) vce(cluster ship) testname("ships: vce(cluster ship)")

* Cluster variable different from FE
clear
set obs 2000
set seed 101
gen cluster_id = mod(_n-1, 50) + 1
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x + 0.05 * fe1)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(cluster cluster_id) testname("synth: cluster != FE")

* Many small clusters
clear
set obs 3000
set seed 102
gen cluster_id = mod(_n-1, 300) + 1
gen fe1 = mod(_n-1, 30) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(cluster cluster_id) testname("synth: many small clusters (300)")

* Few large clusters
clear
set obs 2000
set seed 103
gen cluster_id = mod(_n-1, 5) + 1
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(cluster cluster_id) testname("synth: few large clusters (5)")

* SECTION 6: Frequency weights
noi di as text "--- Section 6: Frequency Weights ---"

clear
set obs 500
set seed 123
gen fe1 = mod(_n-1, 20) + 1
gen x1 = rnormal()
gen x2 = rnormal()
gen mu = exp(0.3 * x1 - 0.2 * x2)
gen y = rpoisson(mu)
gen fw = 1 + int(runiform()*3)

benchmark_ppmlhdfe y x1 x2 [fw=fw], ///
    absorb(fe1) vce(robust) testname("synth: fweight robust")

benchmark_ppmlhdfe y x1 x2 [fw=fw], ///
    absorb(fe1) testname("synth: fweight unadjusted")

benchmark_ppmlhdfe y x1 x2 [fw=fw], ///
    absorb(fe1) vce(cluster fe1) testname("synth: fweight cluster")

* Large frequency weights
clear
set obs 200
set seed 124
gen fe1 = mod(_n-1, 10) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
gen fw = 10 + int(runiform()*90)

benchmark_ppmlhdfe y x [fw=fw], ///
    absorb(fe1) vce(robust) testname("synth: large fweights")

* SECTION 7: Analytic weights
* NOTE: ppmlhdfe does NOT support aweights (error 101). Test that
* cppmlhdfe at least runs without error as a smoke test.
noi di as text "--- Section 7: Analytic Weights (smoke test) ---"

clear
set obs 1000
set seed 130
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
gen aw_var = 1 + runiform() * 5

capture cppmlhdfe y x [aw=aw_var], absorb(fe1) vce(robust)
if _rc == 0 test_pass "synth: aweight robust (smoke test)"
else test_fail "synth: aweight robust" "cppmlhdfe error rc=`=_rc'"

capture cppmlhdfe y x [aw=aw_var], absorb(fe1) vce(cluster fe1)
if _rc == 0 test_pass "synth: aweight cluster (smoke test)"
else test_fail "synth: aweight cluster" "cppmlhdfe error rc=`=_rc'"

* SECTION 8: Probability weights
noi di as text "--- Section 8: Probability Weights ---"

clear
set obs 1000
set seed 135
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
gen pw_var = 1 + runiform() * 5

benchmark_ppmlhdfe y x [pw=pw_var], ///
    absorb(fe1) vce(robust) testname("synth: pweight robust")

benchmark_ppmlhdfe y x [pw=pw_var], ///
    absorb(fe1) vce(cluster fe1) testname("synth: pweight cluster")

* SECTION 9: Weight + VCE combinations
noi di as text "--- Section 9: Weight + VCE Combinations ---"

clear
set obs 1000
set seed 140
gen cluster_id = mod(_n-1, 50) + 1
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
gen fw = 1 + int(runiform()*3)
gen aw_var = 1 + runiform() * 5
gen pw_var = 1 + runiform() * 5

benchmark_ppmlhdfe y x [fw=fw], ///
    absorb(fe1) vce(cluster cluster_id) testname("synth: fweight + cluster")

* aweight + cluster: ppmlhdfe doesn't support aw, smoke test only
capture cppmlhdfe y x [aw=aw_var], absorb(fe1) vce(cluster cluster_id)
if _rc == 0 test_pass "synth: aweight + cluster (smoke test)"
else test_fail "synth: aweight + cluster" "cppmlhdfe error rc=`=_rc'"

benchmark_ppmlhdfe y x [pw=pw_var], ///
    absorb(fe1) vce(cluster cluster_id) testname("synth: pweight + cluster")

* Weights + two-way FE + robust
clear
set obs 1000
set seed 141
gen fe1 = mod(_n-1, 20) + 1
gen fe2 = mod(_n-1, 10) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
gen fw = 1 + int(runiform()*3)

benchmark_ppmlhdfe y x [fw=fw], ///
    absorb(fe1 fe2) vce(robust) testname("synth: fweight + two-way + robust")

* SECTION 10: Factor variables
noi di as text "--- Section 10: Factor Variables ---"

* Auto dataset (mpg as count-like)
sysuse auto, clear
replace mpg = int(mpg)

benchmark_ppmlhdfe mpg i.foreign price weight, ///
    absorb(rep78) vce(robust) testname("auto: i.foreign")

* Synthetic data with multi-level factor
clear
set obs 1000
set seed 150
gen fe1 = mod(_n-1, 10) + 1
gen cat1 = mod(_n-1, 4) + 1
gen x = rnormal()
gen mu = exp(0.3 * x + 0.2 * (cat1==2) + 0.4 * (cat1==3) + 0.6 * (cat1==4))
gen y = rpoisson(mu)

benchmark_ppmlhdfe y i.cat1 x, ///
    absorb(fe1) vce(robust) testname("synth: i.cat1 (4 levels)")

* Factor with continuous interaction
sysuse auto, clear
replace mpg = int(mpg)

benchmark_ppmlhdfe mpg i.foreign##c.weight, ///
    absorb(rep78) vce(robust) testname("auto: i.foreign##c.weight")

* SECTION 11: Multiple covariates
noi di as text "--- Section 11: Multiple Covariates ---"

clear
set obs 2000
set seed 160
gen fe1 = mod(_n-1, 50) + 1
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen mu = exp(0.3*x1 - 0.2*x2 + 0.1*x3)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x1, ///
    absorb(fe1) vce(robust) testname("synth: 1 covariate")

benchmark_ppmlhdfe y x1 x2, ///
    absorb(fe1) vce(robust) testname("synth: 2 covariates")

benchmark_ppmlhdfe y x1 x2 x3, ///
    absorb(fe1) vce(robust) testname("synth: 3 covariates")

benchmark_ppmlhdfe y x1 x2 x3 x4 x5, ///
    absorb(fe1) vce(robust) testname("synth: 5 covariates")

* Covariates with different scales
clear
set obs 2000
set seed 161
gen fe1 = mod(_n-1, 50) + 1
gen x_small = rnormal() / 1000
gen x_medium = rnormal() * 10
gen x_large = rnormal() * 100000
gen mu = exp(100 * x_small + 0.01 * x_medium + 0.000001 * x_large)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x_small x_medium x_large, ///
    absorb(fe1) vce(robust) testname("synth: different covariate scales")

* SECTION 12: if/in conditions
noi di as text "--- Section 12: if/in Conditions ---"

clear
set obs 2000
set seed 170
gen fe1 = mod(_n-1, 20) + 1
gen group = mod(_n-1, 3)
gen x = rnormal()
gen mu = exp(0.3 * x + 0.05 * fe1)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x if group != 0, ///
    absorb(fe1) vce(robust) testname("synth: if condition")

benchmark_ppmlhdfe y x in 1/1000, ///
    absorb(fe1) vce(robust) testname("synth: in condition")

benchmark_ppmlhdfe y x if y > 0, ///
    absorb(fe1) vce(robust) testname("synth: if y > 0")

* if + weights
gen fw = 1 + int(runiform()*3)
benchmark_ppmlhdfe y x [fw=fw] if group == 1, ///
    absorb(fe1) vce(robust) testname("synth: if + fweight")

* in + cluster
benchmark_ppmlhdfe y x in 1/1500, ///
    absorb(fe1) vce(cluster fe1) testname("synth: in + cluster")

* SECTION 13: Separation detection
noi di as text "--- Section 13: Separation Detection ---"

* One FE group all zeros
clear
set obs 1000
set seed 42
gen fe1 = mod(_n-1, 10) + 1
gen x = rnormal()
gen mu = exp(0.5 * x + 0.1 * fe1)
gen y = rpoisson(mu)
replace y = 0 if fe1 == 10

benchmark_ppmlhdfe y x, absorb(fe1) vce(robust) tolerance(1e-12) ///
    testname("separation: single group, full results and sample")
* The helper restores the input data, which invalidates e(sample). Refit to
* test the exact retained observations as well as the reference comparison.
cppmlhdfe y x, absorb(fe1) vce(robust) tolerance(1e-12)
assert e(N) == 900
assert e(sample) == (fe1 != 10)

* Multiple FE groups all zeros
clear
set obs 1000
set seed 43
gen fe1 = mod(_n-1, 10) + 1
gen x = rnormal()
gen mu = exp(0.5 * x + 0.1 * fe1)
gen y = rpoisson(mu)
replace y = 0 if fe1 >= 9

benchmark_ppmlhdfe y x, absorb(fe1) vce(robust) tolerance(1e-12) ///
    testname("separation: two groups, full results and sample")
cppmlhdfe y x, absorb(fe1) vce(robust) tolerance(1e-12)
assert e(N) == 800
assert e(sample) == (fe1 < 9)

* SECTION 14: Missing values
noi di as text "--- Section 14: Missing Values ---"

* Missing in dependent variable
clear
set obs 1000
set seed 180
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
replace y = . in 1/50

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("missing: in depvar (5%)")

* Missing in independent variable
clear
set obs 1000
set seed 181
gen fe1 = mod(_n-1, 20) + 1
gen x1 = rnormal()
gen x2 = rnormal()
replace x2 = . in 51/100
gen mu = exp(0.3 * x1)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x1 x2, ///
    absorb(fe1) vce(robust) testname("missing: in indepvar (5%)")

* Missing in FE variable
clear
set obs 1000
set seed 182
gen fe1 = mod(_n-1, 20) + 1
replace fe1 = . in 901/950
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("missing: in FE var (5%)")

* Heavy missingness (50%)
clear
set obs 2000
set seed 183
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)
replace y = . in 1/1000

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("missing: 50% of depvar")

* SECTION 15: Collinearity
noi di as text "--- Section 15: Collinearity ---"

* Perfectly collinear regressors
clear
set obs 1000
set seed 190
gen fe1 = mod(_n-1, 20) + 1
gen x1 = rnormal()
gen x2 = 2 * x1
gen x3 = rnormal()
gen mu = exp(0.3 * x1 + 0.1 * x3)
gen y = rpoisson(mu)

capture ppmlhdfe y x1 x2 x3, absorb(fe1) vce(robust)
local ppml_rc = _rc
capture cppmlhdfe y x1 x2 x3, absorb(fe1) vce(robust)
local cppml_rc = _rc

if `ppml_rc' == 0 & `cppml_rc' == 0 {
    benchmark_ppmlhdfe y x1 x2 x3, absorb(fe1) vce(robust) tolerance(1e-12) ///
        testname("collinearity: perfect, full results")
}
else if `ppml_rc' == `cppml_rc' {
    test_pass "collinearity: both error (ppml rc=`ppml_rc', cppml rc=`cppml_rc')"
}
else {
    test_fail "collinearity: perfect" "ppml rc=`ppml_rc', cppml rc=`cppml_rc'"
}

* Near-collinear regressors
clear
set obs 2000
set seed 191
gen fe1 = mod(_n-1, 20) + 1
gen x1 = rnormal()
gen x2 = x1 + rnormal() * 1e-6
gen x3 = rnormal()
gen mu = exp(0.3 * x1 + 0.1 * x3)
gen y = rpoisson(mu)

capture ppmlhdfe y x1 x2 x3, absorb(fe1) vce(robust)
local ppml_rc = _rc
capture cppmlhdfe y x1 x2 x3, absorb(fe1) vce(robust)
local cppml_rc = _rc

if `ppml_rc' == 0 & `cppml_rc' == 0 {
    benchmark_ppmlhdfe y x1 x2 x3, absorb(fe1) vce(robust) tolerance(1e-12) ///
        testname("collinearity: near-collinear, full results")
}
else if `ppml_rc' == `cppml_rc' {
    test_pass "collinearity: near-collinear (both error)"
}
else {
    test_fail "collinearity: near-collinear" "ppml rc=`ppml_rc', cppml rc=`cppml_rc'"
}

* SECTION 16: Large synthetic datasets
noi di as text "--- Section 16: Large Datasets ---"

* 10K observations
clear
set seed 200
set obs 10000
gen fe1 = mod(_n-1, 100) + 1
gen x1 = rnormal()
gen x2 = rnormal()
gen mu = exp(0.3*x1 - 0.2*x2 + 0.02*fe1)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x1 x2, ///
    absorb(fe1) vce(robust) testname("10K obs: single FE robust")

benchmark_ppmlhdfe y x1 x2, ///
    absorb(fe1) vce(cluster fe1) testname("10K obs: single FE cluster")

* 10K with two-way FE
clear
set seed 201
set obs 10000
gen fe1 = mod(_n-1, 100) + 1
gen fe2 = mod(_n-1, 50) + 1
gen x = rnormal()
gen mu = exp(0.3*x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) vce(robust) testname("10K obs: two-way FE")

* 50K observations
clear
set seed 202
set obs 50000
gen fe1 = mod(_n-1, 500) + 1
gen x = rnormal()
gen mu = exp(0.3*x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("50K obs: single FE")

* SECTION 17: High-dimensional FE
noi di as text "--- Section 17: High-Dimensional FE ---"

* 1000 FE levels
clear
set seed 210
set obs 10000
gen fe1 = mod(_n-1, 1000) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("1000 FE levels (10K obs)")

* 500 x 200 two-way
clear
set seed 211
set obs 10000
gen fe1 = mod(_n-1, 500) + 1
gen fe2 = mod(_n-1, 200) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) vce(robust) testname("500x200 two-way FE (10K obs)")

* 100 x 50 x 20 three-way
clear
set seed 212
set obs 10000
gen fe1 = mod(_n-1, 100) + 1
gen fe2 = mod(_n-1, 50) + 1
gen fe3 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2 fe3) vce(robust) testname("100x50x20 three-way FE (10K obs)")

* SECTION 18: Pathological numeric values
noi di as text "--- Section 18: Pathological Numeric Values ---"

* Small coefficients (weak effect)
clear
set seed 220
set obs 2000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal() * 100
gen mu = exp(0.001 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("small coefficients")

* Large Y values (high mean)
clear
set seed 221
set obs 1000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(3 + 0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("large Y values (mean ~20)")

* Zero-inflated Y (many zeros)
clear
set seed 222
set obs 2000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(-2 + 0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("zero-inflated Y")

* Mixed covariate scales
clear
set seed 223
set obs 2000
gen fe1 = mod(_n-1, 50) + 1
gen x_tiny = rnormal() / 1e6
gen x_huge = rnormal() * 1e6
gen mu = exp(1e6 * x_tiny + 1e-6 * x_huge)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x_tiny x_huge, ///
    absorb(fe1) vce(robust) testname("mixed scale covariates")

* Constant Y (all same positive value)
clear
set seed 224
set obs 1000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen y = 5

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("constant Y = 5")

* All Y = 0 (degenerate case — all obs separated)
* All observations are separated: the return code must match ppmlhdfe.
clear
set seed 225
set obs 500
gen fe1 = mod(_n-1, 10) + 1
gen x = rnormal()
gen y = 0

capture ppmlhdfe y x, absorb(fe1) vce(robust)
local ppml_rc = _rc
capture cppmlhdfe y x, absorb(fe1) vce(robust)
local cppml_rc = _rc

if `ppml_rc' != 0 & `ppml_rc' == `cppml_rc' {
    test_pass "all Y=0: both error (ppml rc=`ppml_rc', cppml rc=`cppml_rc')"
}
else if `ppml_rc' == 0 & `cppml_rc' == 0 {
    benchmark_ppmlhdfe y x, absorb(fe1) vce(robust) testname("all Y=0: full results")
}
else {
    test_fail "all Y=0" "ppml rc=`ppml_rc', cppml rc=`cppml_rc'"
}

* SECTION 19: Sparse FE / Unbalanced groups
noi di as text "--- Section 19: Sparse FE and Unbalanced Groups ---"

* Very unbalanced group sizes
clear
set obs 2000
set seed 230
gen fe1 = 1 in 1/1000
replace fe1 = 2 in 1001/1500
replace fe1 = runiformint(3, 50) in 1501/2000
gen x = rnormal()
gen mu = exp(0.3 * x + 0.05 * fe1)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("unbalanced groups (50% in group 1)")

* Sparse two-way FE
clear
set seed 231
set obs 1000
gen fe1 = runiformint(1, 200)
gen fe2 = runiformint(1, 100)
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) vce(robust) testname("sparse two-way (200x100, 1K obs)")

* Extreme unbalance (one huge group, many tiny)
clear
set obs 3000
set seed 232
gen fe1 = 1 in 1/2500
replace fe1 = _n - 2499 in 2501/3000
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("extreme unbalance")

* SECTION 20: Singleton patterns
noi di as text "--- Section 20: Singleton Patterns ---"

* Mixed singletons and groups
clear
set obs 500
set seed 240
gen fe1 = _n in 1/100
replace fe1 = runiformint(101, 120) in 101/500
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("mixed singletons (100 of 500)")

* All singletons
clear
set obs 100
set seed 241
gen fe1 = _n
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

capture ppmlhdfe y x, absorb(fe1) vce(robust)
local ppml_rc = _rc
capture cppmlhdfe y x, absorb(fe1) vce(robust)
local cppml_rc = _rc

if `ppml_rc' == 0 & `cppml_rc' == 0 {
    benchmark_ppmlhdfe y x, absorb(fe1) vce(robust) testname("all singletons: full results")
}
else if `ppml_rc' == `cppml_rc' {
    test_pass "all singletons: matching error r(`ppml_rc')"
}
else {
    test_fail "all singletons" "ppml rc=`ppml_rc', cppml rc=`cppml_rc'"
}

* Two-way FE with singletons
clear
set obs 1000
set seed 242
gen fe1 = mod(_n-1, 50) + 1
gen fe2 = _n in 1/200
replace fe2 = runiformint(201, 220) in 201/1000
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1 fe2) vce(robust) testname("two-way FE with singletons")

* SECTION 21: Convergence and solver options
noi di as text "--- Section 21: Convergence Options ---"

clear
set obs 1000
set seed 250
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

* Tight IRLS tolerance
ppmlhdfe y x, absorb(fe1) vce(robust) tolerance(1e-12)
local ppml_ll = e(ll)
tempname ppml_b_tol
matrix `ppml_b_tol' = e(b)

cppmlhdfe y x, absorb(fe1) vce(robust) irlstolerance(1e-12) tolerance(1e-12)
local cppml_ll = e(ll)
tempname cppml_b_tol
matrix `cppml_b_tol' = e(b)

local ppml_bk = `ppml_b_tol'[1,1]
local cppml_bk = `cppml_b_tol'[1,1]
assert_matrix_equal `ppml_b_tol' `cppml_b_tol' $DEFAULT_SIGFIGS "tight tolerance: full b"
assert_scalar_equal `cppml_ll' `ppml_ll' $DEFAULT_SIGFIGS "tight tolerance: log likelihood"

* This fixture cannot converge within two IRLS iterations.
capture cppmlhdfe y x, absorb(fe1) vce(robust) irlsmaxiter(2)
if _rc == 430 test_pass "irlsmaxiter(2): reports non-convergence"
else test_fail "irlsmaxiter(2)" "expected non-convergence r(430), got r(`=_rc')"

* Tight CG tolerance
benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(robust) testname("default tolerance")

* SECTION 22: Verbose option
noi di as text "--- Section 22: Verbose Option ---"

load_ships_data
capture cppmlhdfe accident op_75_79 co_65_69, absorb(ship) vce(robust) verbose
if _rc == 0 test_pass "verbose option accepted"
else test_fail "verbose option" "returned rc=`=_rc'"

* SECTION 23: Error handling
noi di as text "--- Section 23: Error Handling ---"

* Negative dependent variable
clear
set obs 100
set seed 260
gen fe1 = mod(_n-1, 5) + 1
gen x = rnormal()
gen y = rnormal()

capture cppmlhdfe y x, absorb(fe1) vce(robust)
if _rc == 198 test_pass "error: negative Y rejected (rc=`=_rc')"
else test_fail "error: negative Y" "should error on negative depvar"

* No absorbed effects: the C engine fits the constant directly.
clear
set obs 100
set seed 2611
gen x = rnormal()
gen y = rpoisson(exp(0.3*x))
capture noisily {
    quietly ppmlhdfe y x
    matrix noabs_b=e(b)
    matrix noabs_V=e(V)
    quietly cppmlhdfe y x
    assert mreldif(e(b),noabs_b)<1e-6
    assert mreldif(e(V),noabs_V)<1e-6
}
if _rc test_fail "no absorb(): reference parity" "r(`=_rc')"
else test_pass "no absorb(): reference parity"

* No observations after if condition
clear
set obs 100
set seed 261
gen fe1 = mod(_n-1, 5) + 1
gen x = rnormal()
gen y = rpoisson(exp(0.3 * x))

capture cppmlhdfe y x if y > 9999, absorb(fe1) vce(robust)
if _rc == 2000 test_pass "error: no observations (rc=`=_rc')"
else test_fail "error: no observations" "should error when no obs match"

* Negative exposure variable
clear
set obs 100
set seed 262
gen fe1 = mod(_n-1, 5) + 1
gen x = rnormal()
gen y = rpoisson(exp(0.3 * x))
gen bad_exp = rnormal()

capture cppmlhdfe y x, absorb(fe1) exposure(bad_exp)
if _rc == 198 test_pass "error: negative exposure rejected (rc=`=_rc')"
else test_fail "error: negative exposure" "should error on non-positive exposure"

* SECTION 24: Additional real datasets
noi di as text "--- Section 24: Additional Datasets ---"

* Auto dataset variants
sysuse auto, clear
replace mpg = int(mpg)

benchmark_ppmlhdfe mpg weight length, ///
    absorb(foreign) vce(robust) testname("auto: basic")

benchmark_ppmlhdfe mpg weight length, ///
    absorb(foreign rep78) vce(robust) testname("auto: two-way FE")

benchmark_ppmlhdfe mpg weight length, ///
    absorb(foreign) vce(cluster foreign) testname("auto: cluster")

benchmark_ppmlhdfe mpg weight length turn displacement, ///
    absorb(foreign) vce(robust) testname("auto: many covariates")

* Ships dataset variants
load_ships_data

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) vce(cluster ship) testname("ships: cluster ship")

benchmark_ppmlhdfe accident op_75_79 co_65_69 co_70_74 co_75_79, ///
    absorb(ship) exposure(service) vce(cluster ship) testname("ships: exposure + cluster")

* SECTION 25: Stress tests
noi di as text "--- Section 25: Stress Tests ---"

* Many covariates (8)
clear
set seed 300
set obs 5000
gen fe1 = mod(_n-1, 50) + 1
forvalues i = 1/8 {
    gen x`i' = rnormal()
}
gen mu = exp(0.3*x1 - 0.2*x2 + 0.1*x3 + 0.05*x4)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x1 x2 x3 x4 x5 x6 x7 x8, ///
    absorb(fe1) vce(robust) testname("stress: 8 covariates")

* Imbalanced clusters
clear
set seed 301
set obs 5000
gen cluster_id = 1 in 1/4500
replace cluster_id = 2 in 4501/4750
replace cluster_id = runiformint(3, 50) in 4751/5000
gen fe1 = mod(_n-1, 50) + 1
gen x = rnormal()
gen mu = exp(0.3 * x)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x, ///
    absorb(fe1) vce(cluster cluster_id) testname("stress: imbalanced clusters (90% in one)")

* High-dim FE + many covariates
clear
set seed 302
set obs 20000
gen fe1 = mod(_n-1, 2000) + 1
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen mu = exp(0.3*x1 - 0.2*x2 + 0.1*x3)
gen y = rpoisson(mu)

benchmark_ppmlhdfe y x1 x2 x3, ///
    absorb(fe1) vce(robust) testname("stress: 2000 FE levels (20K obs)")

* Weights + two-way FE + cluster + exposure
clear
set seed 303
set obs 3000
gen fe1 = mod(_n-1, 30) + 1
gen fe2 = mod(_n-1, 15) + 1
gen cluster_id = mod(_n-1, 100) + 1
gen x = rnormal()
gen exp_var = 1 + runiform() * 10
gen mu = exp(0.3 * x) * exp_var
gen y = rpoisson(mu)
gen fw = 1 + int(runiform()*3)

benchmark_ppmlhdfe y x [fw=fw], ///
    absorb(fe1 fe2) exposure(exp_var) vce(cluster cluster_id) ///
    testname("stress: fw + two-way + exposure + cluster")

do validation/stata/cppmlhdfe.do cppmlhdfe_precision_cases
do validation/stata/cppmlhdfe.do cppmlhdfe_compatibility_cases
do validation/stata/cppmlhdfe.do cppmlhdfe_upstream_cases
do validation/stata/cppmlhdfe.do cppmlhdfe_advanced_cases

* Wide offsets: rows whose fitted means fall far below ppmlhdfe's 2^-46 mean
* floor (in units of sd(y)). ppmlhdfe's own e(ll)/e(deviance) lose most digits
* there (rounding of the huge IRLS working variable and the floored means),
* and its b at +/-400, so the references are well-conditioned equivalents.
* (a) Offsets +/-s alternating within each absorbed group: raising s by D
* moves the negative-offset means by exp(-2D) and leaves the other rows and b
* unchanged up to exp(-2s), so well-conditioned ppmlhdfe fits give the
* references: ll(s) = ll(10) - 2(s-10)S and dev(s) = dev(10) + 4(s-10)S, with
* S the sum of y over the negative-offset rows, and e(b) from s = 20 (where
* ppmlhdfe's b is still exact; _cons, the weighted mean of the absorbed
* effects, moves by -(s-20), and at s = 10 it is not yet settled to 7 sf).
clear
set seed 504
set obs 1000
gen i = ceil(_n/20)
gen x = rnormal()
gen z = runiform()
gen y = rpoisson(exp(.3 + .5*x - .3*z))
gen double ov = cond(mod(_n, 2), 10, -10)
quietly ppmlhdfe y x z, absorb(i) offset(ov)
local wide_ll10 = e(ll)
local wide_dev10 = e(deviance)
quietly summarize y if ov < 0 & e(sample), meanonly
local wide_S = r(sum)
replace ov = cond(mod(_n, 2), 20, -20)
quietly ppmlhdfe y x z, absorb(i) offset(ov)
matrix wide_b20 = e(b)
foreach s in 20 50 100 400 {
    replace ov = cond(mod(_n, 2), `s', -`s')
    capture cppmlhdfe y x z, absorb(i) offset(ov)
    if _rc {
        test_fail "wide offset +/-`s'" "rc=`=_rc'"
        continue
    }
    matrix wide_b = e(b)
    matrix wide_bs = wide_b20
    matrix wide_bs[1, 3] = wide_b20[1, 3] - (`s' - 20)
    assert_matrix_equal wide_b wide_bs $DEFAULT_SIGFIGS "wide offset +/-`s': e(b)"
    assert_scalar_equal `=e(ll)' `=`wide_ll10' - 2*(`s'-10)*`wide_S'' $DEFAULT_SIGFIGS "wide offset +/-`s': e(ll)"
    assert_scalar_equal `=e(deviance)' `=`wide_dev10' + 4*(`s'-10)*`wide_S'' $DEFAULT_SIGFIGS "wide offset +/-`s': e(deviance)"
}
* (b) Normal offsets without absorb(), where ppmlhdfe's b also drifts (4 sf at
* SD 30): the references are Stata's poisson at tight tolerances for y unrelated
* to the offset, and glm's IRLS for y generated with it (where poisson does not
* converge at tight tolerances); the deviance is 2*(saturated ll - ll), with the
* saturated part by Stirling's series for y >= 10.
foreach design in "unrelated 7" "unrelated 30" "fitted 8" {
    gettoken kind sd : design
    clear
    set seed 504
    set obs 1000
    gen x = rnormal()
    gen z = runiform()
    gen y = rpoisson(exp(.3 + .5*x - .3*z))
    set seed 505
    gen double ov = rnormal()*`sd'
    if "`kind'" == "fitted" replace y = rpoisson(exp(.3 + .5*x - .3*z + ov))
    quietly generate double wide_sat = cond(y == 0, 0, cond(y < 10, ///
        y*ln(y) - y - lngamma(y + 1), ///
        -0.5*ln(2*_pi*y) - 1/(12*y) + 1/(360*y^3) - 1/(1260*y^5) + 1/(1680*y^7)))
    quietly summarize wide_sat if !missing(y), meanonly
    local wide_satsum = r(sum)
    if "`kind'" == "unrelated" {
        quietly poisson y x z, offset(ov) tolerance(1e-13) ltolerance(1e-15) nrtolerance(1e-13)
        local wide_llref = e(ll)
        local wide_devref = 2*(`wide_satsum' - e(ll))
    }
    else {
        quietly glm y x z, family(poisson) offset(ov) irls
        local wide_devref = e(deviance)
        local wide_llref = `wide_satsum' - e(deviance)/2
    }
    matrix wide_bref = e(b)
    capture cppmlhdfe y x z, offset(ov)
    if _rc {
        test_fail "normal offset `design'" "rc=`=_rc'"
        continue
    }
    matrix wide_b = e(b)
    assert_matrix_equal wide_b wide_bref $DEFAULT_SIGFIGS "normal offset `design': e(b)"
    assert_scalar_equal `=e(ll)' `wide_llref' $DEFAULT_SIGFIGS "normal offset `design': e(ll)"
    assert_scalar_equal `=e(deviance)' `wide_devref' $DEFAULT_SIGFIGS "normal offset `design': e(deviance)"
}

* Low-mobility worker-firm panels: six periods, a few percent of moves.
* Some zero outcomes head to a large negative eta without being separated;
* the IRLS must still converge (it stalled in a deviance cycle when each
* iteration re-projected the whole working response), and saved individual
* effects must be recovered.
capture program drop ppml_mobility_panel
program define ppml_mobility_panel
    args seed workers firms moves
    clear
    set seed `seed'
    quietly set obs `workers'
    gen w = _n
    quietly expand 6
    bysort w: gen t = _n
    gen long f = 1 + mod(w, `firms')
    bysort w (t): gen byte mv = runiform() < `moves'
    bysort w (t): replace f = cond(mv, 1 + floor(runiform()*`firms'), f[_n-1]) if _n > 1
    bysort w (t): gen wa = rnormal() if _n == 1
    bysort w (t): replace wa = wa[_n-1] if _n > 1
    gen x = rnormal()
    gen y = rpoisson(exp(1 + .3*x + .2*wa + mod(f,7)/7))
end
* The zero outcomes heading to a large negative eta are separated in the
* limit: there the fit equals the fit without them. ppmlhdfe stops short in
* that direction (its _cons differs from the limit by 6e-8 to 1.2e-7 and its
* e(ll) by 4e-8 to 7e-8, at tolerance 1e-8 and 1e-12 alike), so the
* reference is ppmlhdfe at a tight tolerance on the sample without those
* rows: a y = 0 row whose fitted mean is below 1e-2, where every other fitted
* mean is at least 100 times larger. cppmlhdfe keeps the rows and must reach
* that limit (and e(N) of the full sample): its coefficients at the default
* tolerance and at tolerance(1e-12), with e(ll) and e(r2_p).
foreach design in "2 1000 100 .01" "3 2000 200 .01" "77 400 40 .02" {
    quietly ppml_mobility_panel `design'
    capture cppmlhdfe y x, absorb(w f)
    if _rc {
        test_fail "low-mobility panel `design'" "rc=`=_rc'"
        continue
    }
    matrix mobility_default_b = e(b)
    capture cppmlhdfe y x, absorb(w f) tolerance(1e-12)
    if _rc {
        test_fail "low-mobility panel `design', tolerance(1e-12)" "rc=`=_rc'"
        continue
    }
    matrix mobility_b = e(b)
    matrix mobility_V = e(V)
    local mobility_N = e(N)
    local mobility_ll = e(ll)
    local mobility_r2_p = e(r2_p)
    quietly ppmlhdfe y x, absorb(w f) d(mobility_d) tolerance(1e-12)
    local mobility_ref_N = e(N)
    local mobility_ll_0 = e(ll_0)
    quietly predict double mobility_mu, mu
    quietly summarize mobility_mu if mobility_mu >= 1e-2, meanonly
    local mobility_gap = r(min)
    quietly count if mobility_mu < 1e-2 & y != 0
    local mobility_positive = r(N)
    quietly count if mobility_mu < 1e-2
    local mobility_limit = r(N)
    quietly summarize mobility_mu if mobility_mu < 1e-2, meanonly
    if `mobility_limit' == 0 | `mobility_positive' | r(max) * 100 > `mobility_gap' {
        test_fail "low-mobility panel `design'" "reference limit rows not identified"
        drop mobility_d mobility_mu
        continue
    }
    quietly ppmlhdfe y x if mobility_mu >= 1e-2, absorb(w f) tolerance(1e-12)
    matrix mobility_ref_b = e(b)
    matrix mobility_ref_V = e(V)
    assert_matrix_equal mobility_default_b mobility_ref_b $DEFAULT_SIGFIGS "low-mobility panel `design': e(b) at the default tolerance"
    assert_matrix_equal mobility_b mobility_ref_b $DEFAULT_SIGFIGS "low-mobility panel `design': e(b)"
    * The covariance depends on the fitted means to first order: the IRLS
    * stops with the separated rows' means near 1e-6 (ppmlhdfe near 1e-4),
    * which moves e(V) by up to 2e-7 relative from the limit. With
    * separation(mu), ppmlhdfe's method for such rows, they are dropped once
    * their means fall below mu_tol(), giving the limit's sample.
    capture cppmlhdfe y x, absorb(w f) tolerance(1e-12) separation(fe simplex relu mu)
    if _rc test_fail "low-mobility panel `design', separation(mu)" "rc=`=_rc'"
    else {
        matrix mobility_mu_V = e(V)
        matrix mobility_mu_b = e(b)
        assert_matrix_equal mobility_mu_b mobility_ref_b $DEFAULT_SIGFIGS "low-mobility panel `design', separation(mu): e(b)"
        assert_matrix_equal mobility_mu_V mobility_ref_V $DEFAULT_SIGFIGS "low-mobility panel `design', separation(mu): e(V)"
        if e(N) == `mobility_ref_N' - `mobility_limit' test_pass "low-mobility panel `design', separation(mu): e(N)"
        else test_fail "low-mobility panel `design', separation(mu): e(N)" "`=e(N)'"
    }
    assert_scalar_equal `mobility_ll' `=e(ll)' $DEFAULT_SIGFIGS "low-mobility panel `design': e(ll)"
    assert_scalar_equal `mobility_r2_p' `=1 - e(ll)/`mobility_ll_0'' $DEFAULT_SIGFIGS "low-mobility panel `design': e(r2_p)"
    if `mobility_N' == `mobility_ref_N' test_pass "low-mobility panel `design': e(N)"
    else test_fail "low-mobility panel `design': e(N)" "`mobility_N' vs `mobility_ref_N'"
    drop mobility_d mobility_mu
}
* Saved individual effects on such a panel: recovered (they stopped with
* r(430)), and they reconstruct the saved sum d(). The effects are
* identified through the IRLS weights (the fitted means), so the check is the
* weighted norm of fe1 + fe2 - d() relative to that of d(), at a tight
* tolerance. (Per row, the sum in the rows heading to the separated limit is
* only as accurate as the IRLS projection there.) ppmlhdfe's own saved
* effects differ from its own d() by up to 3e-5 relative per row even at
* tolerance(1e-12), so they are not the reference.
foreach design in "1 1000 100 .02" "2 1000 100 .01" {
    quietly ppml_mobility_panel `design'
    capture cppmlhdfe y x, absorb(fe1=w fe2=f)
    if _rc test_fail "saved effects, low-mobility panel `design'" "rc=`=_rc'"
    else test_pass "saved effects, low-mobility panel `design'"
    capture drop fe1 fe2
    capture cppmlhdfe y x, absorb(fe1=w fe2=f) d(saved_sum) tolerance(1e-12)
    if _rc {
        test_fail "saved effects reconstruct d(), panel `design'" "rc=`=_rc'"
        continue
    }
    quietly predict double saved_mu, mu
    quietly generate double saved_error = saved_mu * (fe1 + fe2 - saved_sum)^2 if e(sample)
    quietly generate double saved_norm = saved_mu * saved_sum^2 if e(sample)
    quietly summarize saved_error, meanonly
    local saved_error = r(sum)
    quietly summarize saved_norm, meanonly
    local saved_sf = -log10(sqrt(`saved_error' / r(sum)))
    if `saved_sf' >= $DEFAULT_SIGFIGS test_pass "saved effects reconstruct d(), panel `design'"
    else test_fail "saved effects reconstruct d(), panel `design'" "weighted sigfigs=`saved_sf'"
}

* Automatic outputs must not replace a variable the model uses
clear
set seed 2411
set obs 200
generate double x = rnormal()
generate g = ceil(_n/10)
generate double y = rpoisson(exp(.4 + .3*x + .02*g))
generate double _ppmlhdfe_d = rnormal()
generate double __hdfe1__ = g
generate double automatic_copy = _ppmlhdfe_d
foreach model in "y x _ppmlhdfe_d, absorb(g) d2" "y x [pw=exp(_ppmlhdfe_d)], absorb(g) d2" ///
    "y x, absorb(g) vce(cluster _ppmlhdfe_d) d2" "y x, absorb(__hdfe1__, savefe)" {
    capture cppmlhdfe `model'
    local automatic_rc = _rc
    capture assert _ppmlhdfe_d == automatic_copy & __hdfe1__ == g
    if `automatic_rc' == 198 & !_rc test_pass "automatic output used in the model: `model'"
    else test_fail "automatic output used in the model: `model'" "rc=`automatic_rc'"
}
drop _ppmlhdfe_d __hdfe1__ automatic_copy

* e(title) and e(marginsnotok) as ppmlhdfe posts them
foreach model in "absorb(g)" "noabsorb" {
    quietly ppmlhdfe y x, `model'
    local ref_title "`e(title)'"
    local ref_marginsnotok "`e(marginsnotok)'"
    quietly cppmlhdfe y x, `model'
    if "`model'" == "absorb(g)" local ref_title "HDFE PPML regression"
    if "`e(title)'" == "`ref_title'" & "`e(marginsnotok)'" == "`ref_marginsnotok'" ///
        test_pass "e(title) and e(marginsnotok): `model'"
    else test_fail "e(title) and e(marginsnotok): `model'" "[`e(title)'] [`e(marginsnotok)']"
}

* Normal offsets with absorb(), SD 6 to 100, stored as float and as double.
* The IRLS must converge (it stopped with r(430): the deviance cycled from
* the rounding of the working values of positive outcomes with means below
* the floor, X'WX was numerically singular far from the optimum, or a full
* step overflowed a mean). ppmlhdfe's b drifts to 3-5 sf at SD 30-100 and its
* e(ll) keeps 4-6 sf, so the reference is Stata's poisson with explicit group
* dummies at tight tolerances, started at ppmlhdfe's b with the group effects
* concentrated out, a_g = ln(Y_g) - ln(sum over g of exp(ov + xb)); it agrees
* with an independent Newton solution of the concentrated likelihood to at
* least 13 sf. The deviance is 2*(saturated ll - ll), as above.
foreach type in float double {
    clear
    set seed 504
    set obs 1000
    gen i = ceil(_n/20)
    gen x = rnormal()
    gen z = runiform()
    gen y = rpoisson(exp(.3 + .5*x - .3*z))
    gen `type' ov = .
    quietly generate double wide_sat = cond(y == 0, 0, cond(y < 10, ///
        y*ln(y) - y - lngamma(y + 1), ///
        -0.5*ln(2*_pi*y) - 1/(12*y) + 1/(360*y^3) - 1/(1260*y^5) + 1/(1680*y^7)))
    quietly summarize wide_sat, meanonly
    local wide_satsum = r(sum)
    sort i, stable
    quietly by i: egen double wide_Y = total(y)
    foreach sd in 6 7 8 10 30 60 100 {
        quietly replace ov = rnormal()*`sd'
        capture ppmlhdfe y x z, absorb(i) offset(ov)
        if _rc {
            test_fail "normal offset SD `sd' (`type'), absorb(): reference" "ppmlhdfe rc=`=_rc'"
            continue
        }
        matrix wide_from = e(b)
        matrix wide_from = wide_from[1, 1..2]
        capture drop wide_v wide_vmax wide_sv wide_a
        quietly generate double wide_v = ov + wide_from[1, 1]*x + wide_from[1, 2]*z
        quietly by i: egen double wide_vmax = max(wide_v)
        quietly by i: egen double wide_sv = total(exp(wide_v - wide_vmax))
        quietly generate double wide_a = ln(wide_Y) - wide_vmax - ln(wide_sv)
        forvalues g = 1/50 {
            quietly summarize wide_a if i == `g', meanonly
            matrix wide_from = wide_from, r(mean)
        }
        capture poisson y x z ibn.i, noconstant offset(ov) from(wide_from, copy) ///
            tolerance(1e-13) ltolerance(1e-15) nrtolerance(1e-13) iterate(100)
        if _rc | e(converged) != 1 {
            test_fail "normal offset SD `sd' (`type'), absorb(): reference" "poisson rc=`=_rc'"
            continue
        }
        matrix wide_bref = (_b[x], _b[z])
        local wide_llref = e(ll)
        capture cppmlhdfe y x z, absorb(i) offset(ov)
        if _rc {
            test_fail "normal offset SD `sd' (`type'), absorb()" "rc=`=_rc'"
            continue
        }
        matrix wide_b = (_b[x], _b[z])
        assert_matrix_equal wide_b wide_bref $DEFAULT_SIGFIGS "normal offset SD `sd' (`type'), absorb(): e(b)"
        assert_scalar_equal `=e(ll)' `wide_llref' $DEFAULT_SIGFIGS "normal offset SD `sd' (`type'), absorb(): e(ll)"
        assert_scalar_equal `=e(deviance)' `=2*(`wide_satsum' - `wide_llref')' $DEFAULT_SIGFIGS "normal offset SD `sd' (`type'), absorb(): e(deviance)"
    }
}

* e(rss), e(rmse), e(ic2), e(title2) and estat ic as after ppmlhdfe; e(footnote)
* and e(estat_cmd) name cppmlhdfe's own programs (ppmlhdfe's are reghdfe's,
* which refuse other commands). e(ic2) counts the partialling subiterations,
* which only one-way demeaning (one per IRLS iteration) shares with ppmlhdfe.
sysuse auto, clear
generate w = round(turn/5)
foreach model in "|absorb(rep78)" "|noabsorb" "|absorb(rep78 foreign)" "|absorb(rep78##c.mpg)" ///
    "|absorb(rep78) vce(cluster foreign)" "[fw=w]|absorb(foreign)" "[pw=turn]|absorb(rep78) vce(robust)" {
    local wgt = substr("`model'", 1, strpos("`model'", "|") - 1)
    local opts = substr("`model'", strpos("`model'", "|") + 1, .)
    quietly ppmlhdfe price weight `wgt', `opts'
    local ref_rss = e(rss)
    local ref_rmse = e(rmse)
    local ref_ic2 = e(ic2)
    local ref_title2 "`e(title2)'"
    quietly estat ic
    matrix ref_ic = r(S)
    quietly cppmlhdfe price weight `wgt', `opts'
    assert_scalar_equal `=e(rss)' `ref_rss' $DEFAULT_SIGFIGS "e(rss): `model'"
    assert_scalar_equal `=e(rmse)' `ref_rmse' $DEFAULT_SIGFIGS "e(rmse): `model'"
    if "`e(title2)'" == "`ref_title2'" test_pass "e(title2): `model'"
    else test_fail "e(title2): `model'" "[`e(title2)'] vs [`ref_title2']"
    if strpos("`model'", "rep78 foreign") | strpos("`model'", "c.mpg") {
        if e(ic2) >= e(ic) & e(ic2) < . test_pass "e(ic2): `model'"
        else test_fail "e(ic2): `model'" "`=e(ic2)'"
    }
    else assert_scalar_equal `=e(ic2)' `ref_ic2' $DEFAULT_SIGFIGS "e(ic2): `model'"
    capture estat ic
    if _rc test_fail "estat ic: `model'" "rc=`=_rc'"
    else {
        matrix cand_ic = r(S)
        assert_matrix_equal cand_ic ref_ic $DEFAULT_SIGFIGS "estat ic: `model'"
    }
    capture `e(footnote)'
    local foot_rc = _rc
    capture estat summarize
    if "`e(footnote)'" != "" & "`e(estat_cmd)'" != "" & !`foot_rc' & !_rc test_pass "e(footnote), e(estat_cmd): `model'"
    else test_fail "e(footnote), e(estat_cmd): `model'" "[`e(footnote)'] [`e(estat_cmd)'] rc=`foot_rc'/`=_rc'"
}

} // end quietly

* (outside quietly: these compare the displayed lines, captured in logs)
* ---------------------------------------------------------------------
* PPML-03..17 (docs/development/FUNCTIONALITY_GAPS.md): replay and its
* options, the header and footnote, messages, e() results, absorb()
* suboptions, partialling and GLM options, tagsep's r2, version/reload and
* predict, each against ppmlhdfe 2.3.3
* ---------------------------------------------------------------------

* The lines a command displays (between markers, in a text log)
capture program drop _pq_lines
program define _pq_lines, rclass
    args file cmd
    capture log close _pqlog
    quietly log using `"`file'"', text replace name(_pqlog)
    noisily di "@@BEGIN"
    capture noisily `cmd'
    local rc = _rc
    noisily di "@@END"
    quietly log close _pqlog
    return scalar rc = `rc'
end

* The next line between the markers of an open log that begins with one of
* the words in keep (any line when keep is empty): r(line), r(eof)
capture program drop _pq_next
program define _pq_next, rclass
    args h keep
    while 1 {
        file read `h' line
        if r(eof) | `"`macval(line)'"' == "@@END" {
            return scalar eof = 1
            exit
        }
        local ok = (`"`keep'"' == "")
        local trimmed = strtrim(`"`macval(line)'"')
        foreach p of local keep {
            if strpos(`"`macval(trimmed)'"', `"`p'"') == 1 local ok 1
        }
        if `ok' {
            return local line = strrtrim(`"`macval(line)'"')
            return scalar eof = 0
            exit
        }
    }
end

* Skip an open log to its first marker
capture program drop _pq_start
program define _pq_start
    args h
    while 1 {
        file read `h' line
        if r(eof) | `"`macval(line)'"' == "@@BEGIN" exit
    }
end

* Whether the lines between the markers of two logs are the same (only
* those beginning with a word of keep, when given): r(same), r(n) lines
capture program drop _pq_same
program define _pq_same, rclass
    args f1 f2 keep
    tempname h1 h2
    file open `h1' using `"`f1'"', read text
    file open `h2' using `"`f2'"', read text
    _pq_start `h1'
    _pq_start `h2'
    local same 1
    local n 0
    while 1 {
        _pq_next `h1' `"`keep'"'
        local e1 = r(eof)
        local l1 `"`r(line)'"'
        _pq_next `h2' `"`keep'"'
        local e2 = r(eof)
        local l2 `"`r(line)'"'
        if `e1' & `e2' continue, break
        if `e1' != `e2' | `"`macval(l1)'"' != `"`macval(l2)'"' {
            noisily di as text `"REF: `macval(l1)'"'
            noisily di as text `"CT:  `macval(l2)'"'
            local same 0
            continue, break
        }
        local ++n
    }
    file close `h1'
    file close `h2'
    return scalar same = `same'
    return scalar n = `n'
end

* The number of lines between the markers of a log that contain text
capture program drop _pq_count
program define _pq_count, rclass
    args file text
    tempname h
    file open `h' using `"`file'"', read text
    _pq_start `h'
    local n 0
    while 1 {
        _pq_next `h' ""
        if r(eof) continue, break
        if strpos(`"`r(line)'"', `"`text'"') local ++n
    }
    file close `h'
    return scalar n = `n'
end

* The replays of ppmlhdfe and cppmlhdfe after the same model, with the
* same replay options, show the same lines
capture program drop _pq_replay
program define _pq_replay
    syntax, MODEL(string) NAME(string) [ROPTions(string)]
    tempfile f1 f2
    capture noisily {
        quietly ppmlhdfe `model'
        _pq_lines `"`f1'"' "ppmlhdfe, `roptions'"
        local rc1 = r(rc)
        quietly cppmlhdfe `model'
        _pq_lines `"`f2'"' "cppmlhdfe, `roptions'"
        local rc2 = r(rc)
        _pq_same `"`f1'"' `"`f2'"' ""
        local same = r(same) & r(n) > 0
    }
    if _rc test_fail "`name'" "r(`=_rc')"
    else if `rc1' | `rc2' | !`same' test_fail "`name'" "rc `rc1'/`rc2', same `same'"
    else test_pass "`name'"
end

* The notes and warnings ppmlhdfe and cppmlhdfe print while fitting a model
capture program drop _pq_notes
program define _pq_notes
    syntax, MODEL(string) NAME(string) [NONE]
    tempfile f1 f2
    local keep "(dropped (ReLU (simplex (IRLS warning: (identifying"
    capture noisily {
        _pq_lines `"`f1'"' "ppmlhdfe `model'"
        local rc1 = r(rc)
        _pq_lines `"`f2'"' "cppmlhdfe `model'"
        local rc2 = r(rc)
        _pq_same `"`f1'"' `"`f2'"' `"`keep'"'
        local same = r(same) & (r(n) > 0) == ("`none'" == "")
    }
    if _rc test_fail "`name'" "r(`=_rc')"
    else if `rc1' | `rc2' | !`same' test_fail "`name'" "rc `rc1'/`rc2', same `same'"
    else test_pass "`name'"
end

clear
set seed 12345
set obs 2000
generate id1 = ceil(_n/40)
generate id2 = mod(_n, 20) + 1
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double y = rpoisson(exp(0.3*x1 - 0.2*x2 + 0.1*mod(id1,5)))
replace y = 0 if id1 == 3
generate double yl = y * 1e-8
replace yl = 50 in 1
generate double rrev = x1
generate byte sepx = (id1 == 7 & y == 0)

* PPML-03, PPML-11, PPML-12: replay with header and footnote, and its
* eform, noheader, notable and display options
foreach model in "y x1 x2, absorb(id1 id2)" "y x1 x2, absorb(id1) vce(cluster id2)" ///
    "y x1 x2, absorb(id1) vce(cluster id1 id2)" "y x1 x2" "y x1 x2, absorb(id1##c.x2 id2)" ///
    "y x1 x2 [fw=id2], absorb(id1)" {
    _pq_replay, model(`model') name("PPML-03/11/12 replay: `model'")
}
foreach ropt in eform noheader notable "level(90)" "eform level(80) noheader" {
    _pq_replay, model(y x1 x2, absorb(id1 id2)) roptions(`ropt') name("PPML-03 replay, `ropt'")
}
* irr at estimation shows exp(b), as ppmlhdfe's
_pq_replay, model(y x1 x2, absorb(id1 id2) irr) name("PPML-03 irr shows exp(b)")
capture noisily {
    quietly cppmlhdfe y x1 x2, absorb(id1 id2)
    cppmlhdfe, eform
    confirm matrix r(table)
    assert colsof(r(table)) == 3
}
if _rc test_fail "PPML-03 replay returns r(table)" "r(`=_rc')"
else test_pass "PPML-03 replay returns r(table)"

* PPML-13, PPML-14, PPML-15: e(title2), e(title3), e(N_hdfe) (absorb()
* terms), e(N_hdfe_extended) and e(indepvars) (with _cons), as ppmlhdfe
foreach model in "y x1 x2, absorb(id1 id2)" "y x1 x2, absorb(id1) vce(cluster id2)" ///
    "y x1 x2, absorb(i.id1#i.id2 id1)" "y x1 x2, absorb(id1 id1#c.x2)" "y x1 x2, absorb(id1##c.x2 id2)" ///
    "y x1, absorb(id2##c.(x2) id1)" "y x1 x2" "y i.id2 x1, absorb(id1)" "y x1, absorb(id1 id2)" {
    quietly ppmlhdfe `model'
    local r_title2 `"`e(title2)'"'
    local r_title3 `"`e(title3)'"'
    local r_indep `"`e(indepvars)'"'
    local r_N_hdfe = e(N_hdfe)
    local r_ext = e(N_hdfe_extended)
    quietly cppmlhdfe `model'
    if `"`r_title2'"' == `"`e(title2)'"' & `"`r_title3'"' == `"`e(title3)'"' & ///
        `"`r_indep'"' == `"`e(indepvars)'"' & `r_N_hdfe' == e(N_hdfe) & `r_ext' == e(N_hdfe_extended) ///
        test_pass "PPML-13/14/15 e() names: `model'"
    else test_fail "PPML-13/14/15 e() names: `model'" ///
        `"[`e(title2)'|`e(title3)'|`e(indepvars)'|`=e(N_hdfe)'] vs [`r_title2'|`r_title3'|`r_indep'|`r_N_hdfe']"'
}

* PPML-17 and the notes: dropped singletons, keepsingletons, a low
* standardized outcome, ReLU, simplex and IRLS separation, as ppmlhdfe
foreach model in "y x1 x2, absorb(id1 id2)" "y x1 x2, absorb(id1 id2, keepsingletons)" ///
    "yl x1 x2, absorb(id1 id2)" "y x1 sepx, absorb(id2) separation(mu)" ///
    "y x1 x2, absorb(id1 id2) separation(relu)" "y x1 sepx, absorb(id2) separation(simplex)" ///
    "y x1 sepx, absorb(id2) separation(relu)" ///
    "y x1 x2, absorb(id1 id2) separation(relu) keepsingletons" {
    _pq_notes, model(`model') name("PPML-17 notes: `model'")
}

* verbose(-1), and the HDFE object's verbose(-1), hide them
foreach model in "y x1 x2, absorb(id1 id2, verbose(-1))" "y x1 x2, absorb(id1 id2) verbose(-1)" {
    _pq_notes, model(`model') name("PPML-17 notes: `model'") none
}

* PPML-05: absorb(, maxiterations() verbose()) as ppmlhdfe accepts them;
* PPML-06: prune, condition, rre(), algorithm(), slopemethod() leave the
* estimates as they are; condition posts ppmlhdfe's e(finite_condition)
quietly cppmlhdfe y x1 x2, absorb(id1 id2)
matrix pq_b0 = e(b)
foreach opts in "absorb(id1 id2, maxiterations(500) verbose(0))" "absorb(id1 id2, maxit(500))" ///
    "absorb(id1 id2) condition" "absorb(id1 id2, condition)" "absorb(id1 id2) rre(rrev)" ///
    "absorb(id1 id2) algorithm(cg)" "absorb(id1 id2) slopemethod(foo)" "absorb(id1 id2, algorithm(map))" {
    capture quietly ppmlhdfe y x1 x2, `opts'
    local rc1 = _rc
    local r_fc = e(finite_condition)
    capture quietly cppmlhdfe y x1 x2, `opts'
    local rc2 = _rc
    if !`rc1' & !`rc2' & mreldif(e(b), pq_b0) < 1e-7 & ///
        (e(finite_condition) == `r_fc' | (missing(e(finite_condition)) & missing(`r_fc'))) ///
        test_pass "PPML-05/06 `opts'"
    else test_fail "PPML-05/06 `opts'" "rc `rc1'/`rc2'"
}
* prune is accepted (ppmlhdfe's needs moremata, not installed here)
capture quietly cppmlhdfe y x1 x2, absorb(id1 id2) prune
if !_rc & mreldif(e(b), pq_b0) < 1e-7 test_pass "PPML-06 prune"
else test_fail "PPML-06 prune" "r(`=_rc')"
* condition requires two absorbed effects; precondition and residuals()
* are refused, as by ppmlhdfe (whose codes come from Mata: 3598, 3498)
foreach opts in "absorb(id1) condition" "noabsorb condition" "absorb(id1 id2) precondition" ///
    "absorb(id1 id2) residuals(rr)" "absorb(id1 id2, verbose(6))" "absorb(id1 id2) rre(nosuchvar)" {
    capture quietly ppmlhdfe y x1 x2, `opts'
    local rc1 = _rc
    capture quietly cppmlhdfe y x1 x2, `opts'
    local rc2 = _rc
    if `rc1' & `rc2' & (`rc2' == `rc1' | inlist(`rc1', 3010, 3498, 3598)) test_pass "PPML-06 refused: `opts'"
    else test_fail "PPML-06 refused: `opts'" "rc `rc1'/`rc2'"
}

* PPML-07: acceleration() and transform() values as reghdfe reads them
foreach a in off CG Conj sd stee ai none lsmr steepest_descent {
    capture quietly ppmlhdfe y x1 x2, absorb(id1 id2) acceleration(`a')
    local rc1 = _rc
    capture quietly cppmlhdfe y x1 x2, absorb(id1 id2) acceleration(`a')
    local rc2 = _rc
    if !`rc1' & !`rc2' & mreldif(e(b), pq_b0) < 1e-6 test_pass "PPML-07 acceleration(`a')"
    else test_fail "PPML-07 acceleration(`a')" "rc `rc1'/`rc2'"
}
foreach t in CIM kac sym Symmetric_Kaczmarz cimmino {
    capture quietly ppmlhdfe y x1 x2, absorb(id1 id2) transform(`t')
    local rc1 = _rc
    capture quietly cppmlhdfe y x1 x2, absorb(id1 id2) transform(`t')
    local rc2 = _rc
    if !`rc1' & !`rc2' & mreldif(e(b), pq_b0) < 1e-6 test_pass "PPML-07 transform(`t')"
    else test_fail "PPML-07 transform(`t')" "rc `rc1'/`rc2'"
}
foreach o in "acceleration(bogus)" "transform(bogus)" "transform(rand)" {
    capture quietly ppmlhdfe y x1 x2, absorb(id1 id2) `o'
    local rc1 = _rc
    capture quietly cppmlhdfe y x1 x2, absorb(id1 id2) `o'
    local rc2 = _rc
    if `rc1' & `rc2' & (`rc2' == `rc1' | inlist(`rc1', 3010, 3498, 3598)) test_pass "PPML-07 refused: `o'"
    else test_fail "PPML-07 refused: `o'" "rc `rc1'/`rc2'"
}

* PPML-09: ppmlhdfe's GLM properties, at the top level and in absorb()
foreach prop in "target_inner_tol(1e-10)" "relu_debug(0)" "relu_report_r2(0)" "min_ok(2)" ///
    "standardize_data(1)" "use_exact_partial(1)" "simplex_tol(1e-11)" "mu_tol(1e-7)" {
    capture quietly ppmlhdfe y x1 x2, absorb(id1 id2) `prop'
    local rc1 = _rc
    matrix pq_r = e(b)
    capture quietly cppmlhdfe y x1 x2, absorb(id1 id2) `prop'
    local rc2 = _rc
    matrix pq_c = e(b)
    capture quietly cppmlhdfe y x1 x2, absorb(id1 id2, `prop')
    local rc3 = _rc
    if !`rc1' & !`rc2' & !`rc3' & mreldif(pq_c, pq_r) < 1e-5 & mreldif(e(b), pq_c) == 0 ///
        test_pass "PPML-09 `prop'"
    else test_fail "PPML-09 `prop'" "rc `rc1'/`rc2'/`rc3'"
}
* target_inner_tol() is itolerance()
quietly cppmlhdfe y x1 x2, absorb(id1 id2) itolerance(1e-11)
matrix pq_c = e(b)
quietly cppmlhdfe y x1 x2, absorb(id1 id2) target_inner_tol(1e-11)
if mreldif(e(b), pq_c) == 0 & e(ic2) < . test_pass "PPML-09 target_inner_tol() is itolerance()"
else test_fail "PPML-09 target_inner_tol() is itolerance()" "b differs"
foreach o in "absorb(id1 id2, bogus_prop(1))" "absorb(id1 id2, min_o(2))" "absorb(id1 id2) bogus_prop(1)" {
    capture quietly ppmlhdfe y x1 x2, `o'
    local rc1 = _rc
    capture quietly cppmlhdfe y x1 x2, `o'
    local rc2 = _rc
    if `rc1' & `rc2' == `rc1' test_pass "PPML-09 refused: `o'"
    else test_fail "PPML-09 refused: `o'" "rc `rc1'/`rc2'"
}

* PPML-08: r2 with tagsep() and zvarname() regresses the certificate on the
* model (creghdfe here, reghdfe in ppmlhdfe): R-squared 1
tempfile pq_f1 pq_f2
capture drop pq_sep pq_z
_pq_lines `"`pq_f1'"' "ppmlhdfe y x1 x2, absorb(id1 id2) tagsep(pq_sep) zvarname(pq_z) r2"
local rc1 = r(rc)
local r_r2 = e(r2)
quietly count if pq_sep == 1
local r_nsep = r(N)
capture drop pq_sep pq_z
_pq_lines `"`pq_f2'"' "cppmlhdfe y x1 x2, absorb(id1 id2) tagsep(pq_sep) zvarname(pq_z) r2"
local rc2 = r(rc)
local c_r2 = e(r2)
quietly count if pq_sep == 1
local c_nsep = r(N)
_pq_same `"`pq_f1'"' `"`pq_f2'"' "(identifying (ReLU Verifying"
if !`rc1' & !`rc2' & r(same) & r(n) == 3 & `r_nsep' == `c_nsep' & `c_nsep' > 0 & abs(`c_r2' - 1) < 1e-8 & ///
    abs(`r_r2' - 1) < 1e-8 & "`e(cmd)'" == "creghdfe" test_pass "PPML-08 r2 with tagsep()"
else test_fail "PPML-08 r2 with tagsep()" "rc `rc1'/`rc2' sep `r_nsep'/`c_nsep' r2 `r_r2'/`c_r2'"
capture drop pq_sep pq_z
* relu_report_r2(1) is r2; without zvarname() nothing is regressed
capture drop pq_sep pq_z
capture noisily quietly cppmlhdfe y x1 x2, absorb(id1 id2) tagsep(pq_sep) zvarname(pq_z) relu_report_r2(1)
if !_rc & "`e(cmd)'" == "creghdfe" test_pass "PPML-08 relu_report_r2(1)"
else test_fail "PPML-08 relu_report_r2(1)" "r(`=_rc') `e(cmd)'"
capture drop pq_sep pq_z
ereturn clear
capture quietly cppmlhdfe y x1 x2, absorb(id1 id2) tagsep(pq_sep) r2
if !_rc & "`e(cmd)'" == "" test_pass "PPML-08 r2 without zvarname()"
else test_fail "PPML-08 r2 without zvarname()" "r(`=_rc') `e(cmd)'"
capture drop pq_sep pq_z

* PPML-10: cppmlhdfe, version and cppmlhdfe, reload (standalone)
ereturn clear
_pq_lines `"`pq_f2'"' "cppmlhdfe, version"
local rc_v = r(rc)
_pq_count `"`pq_f2'"' "Required packages installed?"
local pq_which = r(n)
capture quietly cppmlhdfe, reload
local rc_r = _rc
capture quietly cppmlhdfe, version level(90)
local rc_bad = _rc
if !`rc_v' & `pq_which' == 1 & !`rc_r' & `rc_bad' == 301 test_pass "PPML-10 version and reload"
else test_fail "PPML-10 version and reload" "rc `rc_v'/`rc_r'/`rc_bad' which `pq_which'"

* PPML-16: predict's statistics, labels and the default's note, as
* ppmlhdfe_p
local pq_stats "default mu eta xbd xb d response pearson deviance anscombe working scores stdp"
foreach model in "y x1 x2, absorb(id1 id2) d" "y x1 x2" "y x1 x2 if id2 > 3, absorb(id1 id2) d" {
    quietly ppmlhdfe `model'
    local i 0
    foreach stat of local pq_stats {
        local ++i
        local o = cond("`stat'" == "default", "", ", `stat'")
        capture drop pq_r`i'
        * (ppmlhdfe_p's d fails without absorbed effects; cppmlhdfe_p gives 0)
        capture quietly predict pq_r`i' `o'
        local r_rc`i' = _rc
        if !_rc local r_lab`i' : variable label pq_r`i'
    }
    quietly cppmlhdfe `model'
    local fails
    local i 0
    foreach stat of local pq_stats {
        local ++i
        local o = cond("`stat'" == "default", "", ", `stat'")
        capture drop pq_c
        quietly predict pq_c `o'
        if `r_rc`i'' {
            drop pq_c
            continue
        }
        local c_lab : variable label pq_c
        capture assert reldif(pq_r`i', pq_c) < 1e-6 | (missing(pq_r`i') & missing(pq_c))
        if _rc | `"`r_lab`i''"' != `"`c_lab'"' local fails `fails' `stat'
        drop pq_r`i' pq_c
    }
    if "`fails'" == "" test_pass "PPML-16 predict labels and values: `model'"
    else test_fail "PPML-16 predict labels and values: `model'" "`fails'"
    capture drop _ppmlhdfe_d
}
quietly ppmlhdfe y x1 x2, absorb(id1 id2) d(pq_d1)
_pq_lines `"`pq_f1'"' "predict pq_m1"
quietly cppmlhdfe y x1 x2, absorb(id1 id2) d(pq_d2)
_pq_lines `"`pq_f2'"' "predict pq_m2"
_pq_same `"`pq_f1'"' `"`pq_f2'"' ""
if r(same) & r(n) == 1 test_pass "PPML-16 predict's default note"
else test_fail "PPML-16 predict's default note" "lines differ"
capture drop pq_d1 pq_d2 pq_m1 pq_m2

* Summary
print_summary "cppmlhdfe"

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "cppmlhdfe"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe"
    exit
}

if "`__case'" == "cppmlhdfe_advanced_cases" {
* Explicit reference controls and initialization methods.
clear
set seed 84123
set obs 700
gen g=runiformint(1,20)
gen h=runiformint(1,15)
gen double x=rnormal()
gen double z=rnormal()
gen double y=rpoisson(exp(.4*x-.2*z+.03*g))
gen double initial=exp(.3*x+.02*g)
foreach option in "guess(simple)" "guess(ols)" "guess(variable initial)" ///
    "standardize_data(0)" "min_ok(3)" "use_exact_solver(1)" ///
    "use_exact_partial(0)" "use_heuristic_tol(0)" "start_inner_tol(.001)" ///
    "use_step_halving(1) step_halving_memory(.8) max_step_halving(4)" ///
    "remove_collinear_variables(0)" "separation(mu)" "separation(fe mu)" ///
    "relu_tol(1e-6) relu_zero_tol(1e-10) relu_maxiter(200) relu_strict(1)" {
    ppml_native_compare y x z, absorb(g h) extra("`option'") name("advanced `option'") predict
}
* A separated dummy exercises the explicit mu-only removal and native refit.
gen byte separated=(_n<=70)
replace y=0 if separated
ppml_native_compare y x z separated, absorb(g h) extra("separation(mu)") name("mu-only separation")

* Named effects and absorb suboptions must reconstruct the saved sum in C.
capture noisily {
    cppmlhdfe y x z, absorb(firm_fe=g year_fe=h) d(total_fe) tol(1e-12)
    assert missing(firm_fe)==!e(sample)
    assert missing(year_fe)==!e(sample)
    assert abs(firm_fe+year_fe-total_fe)<1e-7 if e(sample)
    bysort g: assert firm_fe==firm_fe[1] if e(sample) & !missing(firm_fe[1])
    drop firm_fe year_fe total_fe
    cppmlhdfe y x, absorb(trend_fe=g##c.z year_fe=h, tolerance(1e-12)) d(total_fe) tol(1e-12)
    assert abs(trend_fe+trend_feSlope1*z+year_fe-total_fe)<1e-7 if e(sample)
    drop trend_fe trend_feSlope1 year_fe total_fe
    cppmlhdfe y x, absorb(g h, savefe keepsingletons) d(total_fe) tol(1e-12)
    assert abs(__hdfe1__+__hdfe2__-total_fe)<1e-7 if e(sample)
    assert e(drop_singletons)==0
    cppmlhdfe y x, absorb(g h, savefe) tol(1e-12)
    assert !missing(__hdfe1__,__hdfe2__) if e(sample)
    drop total_fe
    gen group_one=g
    gen group_two=h
    cppmlhdfe y x, absorb(group_*, savefe) d(total_fe) tol(1e-12)
    assert abs(__hdfe1__+__hdfe2__-total_fe)<1e-7 if e(sample)
    drop total_fe
    cppmlhdfe y x, absorb(trend_fe=g##c.(z initial) year_fe=h) d(total_fe) tol(1e-12)
    assert abs(trend_fe+trend_feSlope1*z+trend_feSlope2*initial+year_fe-total_fe)<1e-7 if e(sample)
}
if _rc test_fail "native named and saved fixed effects" "r(`=_rc')"
else test_pass "native named and saved fixed effects"

* Individual separation methods are independent options.
foreach method in simplex relu "fe simplex" "fe relu" "fe simplex relu mu" {
    ppml_native_compare y x z separated, absorb(g h) extra("separation(`method')") name("separation `method'")
}

foreach option in "accel(cg)" "accel(none)" "accel(sd)" "accel(aitken)" ///
    "accel(hybrid)" "accel(lsmr) btol(1e-12)" "accel(sd) transform(kaczmarz)" ///
    "transform(cimmino)" "poolsize(1)" {
    ppml_native_compare y x z separated, absorb(g h) nativeopts("`option'") name("projection `option'")
}
capture noisily {
    cppmlhdfe y x z separated, absorb(g h) separation(relu) ///
        tagsep(is_separated) zvarname(certificate) relu_accelerate(1) tol(1e-12)
    assert is_separated==separated
    assert abs(certificate)<1e-10 if y>0
    assert certificate>0 if separated
    assert certificate>=0 if y==0
}
if _rc test_fail "native ReLU tags and certificate" "r(`=_rc')"
else test_pass "native ReLU tags and certificate"

* tagsep selects diagnostic-only ReLU and overrides FE/simplex screening.
preserve
capture noisily {
    replace y=0 if g==1
    ppmlhdfe y x z separated, absorb(g h) separation(fe relu) tagsep(ref_tag) zvarname(ref_certificate)
    cppmlhdfe y x z separated, absorb(g h) separation(fe relu) tagsep(native_tag) zvarname(native_certificate)
    assert native_tag==ref_tag
    assert native_tag==1 if g==1
    assert missing(native_certificate)==missing(ref_certificate)
    assert native_certificate==0 if y>0 & !missing(ref_certificate)
    assert native_certificate>0 if native_tag==1
}
if _rc test_fail "native ReLU diagnostic sample" "r(`=_rc')"
else test_pass "native ReLU diagnostic sample"
restore

* Solver controls also interact with weights, offset, slopes, and warm partials.
clear
set seed 88231
set obs 1200
gen g=runiformint(1,20)
gen h=runiformint(1,15)
gen double x=rnormal()
gen double z=rnormal()
gen double off=.2*rnormal()
gen double pw=exp(rnormal())
gen double y=rpoisson(exp(.4*x+.01*g*z+.02*h+off))
gen double poor_initial=1000+_n
foreach method in cg none sd aitken hybrid lsmr {
    ppml_native_compare y x [pw=pw], absorb(g##c.z h) offset(off) vce(cluster g h) ///
        nativeopts("accel(`method') use_exact_partial(0) use_exact_solver(1) btol(1e-12)") ///
        name("weighted slope projection `method'") predict
}
ppml_native_compare y x z [pw=pw], absorb(g h) offset(off) ///
    extra("guess(variable poor_initial) use_step_halving(1) max_step_halving(8)") ///
    name("poor initial means with step halving") predict

capture noisily {
    foreach invalid in "tolerance(0)" "iterate(0)" "poolsize(-1)" {
        capture cppmlhdfe y x, absorb(g h, `invalid')
        assert _rc==198
    }
    foreach invalid in "btol(.)" "conlim(.)" "relu_tol(0)" "guess(nonsense)" "relu_maxiter(.)" "threads(-1)" {
        capture cppmlhdfe y x, absorb(g h) `invalid'
        assert _rc==198
    }
}
if _rc test_fail "advanced invalid controls" "r(`=_rc')"
else test_pass "advanced invalid controls"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_advanced_cases"
    exit
}

if "`__case'" == "cppmlhdfe_compatibility_cases" {
* Coefficient and covariance agreement in significant figures. An entry
* that is zero in exact arithmetic (a symmetric design, an exactly fitted
* model) is reported at roundoff by both programs, and relative digits of two
* roundoff values carry no information. Such an entry agrees when both values
* lie below a roundoff bound in the entry's own units, as in
* sigfigs_cancellation (16*epsdouble):
*   b_j:  16*eps*SE_j, roundoff at the coefficient's sampling scale.
*   V_ij: max(16*eps*SE_i*SE_j, r_i*r_j) with r_j = 16*eps*|b_j|. The first
*         term is the benchmark_helpers covariance bound; the second treats a
*         standard error below 16 ulps of its coefficient as unresolvable.
* Scales come from the reference (or analytic oracle) b and V. Every other
* entry must agree to the requested significant figures.
capture program drop ppml_min_sigfigs
program define ppml_min_sigfigs, rclass
    args A B kind refb refV
    tempname a c bound
    local min_sf 15
    local tol=16*c(epsdouble)
    forvalues i=1/`=rowsof(`A')' {
        forvalues j=1/`=colsof(`A')' {
            scalar `a'=`A'[`i',`j']
            scalar `c'=`B'[`i',`j']
            if "`kind'"=="b" scalar `bound'=`tol'*sqrt(`refV'[`j',`j'])
            else scalar `bound'=max(`tol'*sqrt(`refV'[`i',`i'])*sqrt(`refV'[`j',`j']), ///
                `tol'^2*abs(`refb'[1,`i']*`refb'[1,`j']))
            if missing(`bound') scalar `bound'=0
            if !missing(`a',`c') & abs(`a')<=`bound' & abs(`c')<=`bound' continue
            sigfigs `a' `c'
            local min_sf=min(`min_sf',r(sigfigs))
        }
    }
    return scalar min_sigfigs=`min_sf'
end

* Native PPML compatibility: compare full results and retained rows, then
* exercise prediction after estimates store/restore. No reference fallback.
capture program drop ppml_native_compare
program define ppml_native_compare
    syntax varlist(fv ts) [fw pw] [if] [in], NAME(string) [ABSorb(string) VCE(string) OFFset(varname) EXTRA(string) REFerenceopts(string) NATIVEopts(string) PREDICT VCEREF(name)]
    local wt
    if "`weight'"!="" local wt [`weight'`exp']
    local opts `extra'
    if "`absorb'"!="" local opts `opts' absorb(`absorb')
    if "`vce'"!="" local opts `opts' vce(`vce')
    if "`offset'"!="" local opts `opts' offset(`offset')
    local precise_ref
    foreach control in use_exact_partial use_exact_solver min_ok {
        if !strpos("`opts' `referenceopts'","`control'(") {
            local value=cond("`control'"=="min_ok",3,1)
            local precise_ref `precise_ref' `control'(`value')
        }
    }
    preserve
    capture noisily {
        tempvar refsample rd cd
        ereturn clear
        quietly ppmlhdfe `varlist' `wt' `if' `in', `opts' d(`rd') tol(1e-14) ///
            `precise_ref' `referenceopts'
        matrix rb=e(b)
        matrix rv=e(V)
        * An explicit independent oracle is needed when Stata rejects the
        * upstream covariance at posting. All other result checks still apply.
        if "`vceref'" != "" matrix rv=`vceref'
        gen byte `refsample'=e(sample)
        foreach stat in N df df_a df_a_initial df_a_redundant df_a_nested df_m rank num_singletons num_separated {
            scalar ref_`stat'=e(`stat')
        }
        local prediction_types xb mu eta xbd d stdp response scores pearson deviance anscombe working
        if "`predict'"!="" {
            foreach opt of local prediction_types {
                tempvar rp_`opt'
                quietly predict double `rp_`opt'', `opt'
            }
        }
        ereturn clear
        quietly cppmlhdfe `varlist' `wt' `if' `in', `opts' `nativeopts' d(`cd') tol(1e-14)
        matrix cb=e(b)
        matrix cv=e(V)
        assert `refsample'==e(sample)
        foreach stat in N df df_a df_a_initial df_a_redundant df_a_nested df_m rank num_singletons num_separated {
            if reldif(e(`stat'),ref_`stat')>=1e-12 noi di "MISMATCH `name': `stat' reference=" ref_`stat' " native=" e(`stat')
            assert reldif(e(`stat'),ref_`stat')<1e-12
        }
        assert colsof(rb)==colsof(cb)
        if mreldif(rb,cb)>=1e-7 noi di "MISMATCH `name': coefficients " mreldif(rb,cb)
        ppml_min_sigfigs rb cb b rb rv
        assert r(min_sigfigs)>=$DEFAULT_SIGFIGS
        ppml_min_sigfigs rv cv V rb rv
        assert r(min_sigfigs)>=$DEFAULT_SIGFIGS
        forvalues j=1/`=colsof(cb)' {
            assert abs(sqrt(cv[`j',`j'])-sqrt(rv[`j',`j'])) < 1e-6*sqrt(rv[`j',`j'])+1e-12
            if "`vceref'" != "" & rv[`j',`j'] > 0 {
                assert abs(cv[`j',`j']/rv[`j',`j']-1) < 1e-7
            }
        }
        if "`predict'"!="" {
            estimates store native_restore
            quietly regress `: word 1 of `varlist''
            estimates restore native_restore
            * Predictions must agree to 7 significant figures (1e-7) relative
            * to max(|entry|, S). The linear index xb + d = eta is built from
            * global quantities (coefficients, fixed effects) whose errors do
            * not shrink where a row's sum crosses zero, so for xb, d, eta and
            * xbd S is their common sup norm over the estimation sample.
            * Without fixed effects d is identically zero and pure roundoff in
            * both programs. Residual-type predictions are differences of y
            * and mu terms, and stdp is the root of a signed quadratic form
            * x'Vx; each uses its own sup norm (normwise relative error). mu
            * has no cancellation and keeps elementwise relative error (S=0).
            local linear_scale 0
            foreach comp in xb d eta {
                quietly summarize `rp_`comp'' if `refsample', meanonly
                if r(N) local linear_scale=max(`linear_scale',abs(r(min)),abs(r(max)))
            }
            foreach opt of local prediction_types {
                local scale 0
                if inlist("`opt'","xb","d","eta","xbd") local scale `linear_scale'
                else if "`opt'"!="mu" {
                    quietly summarize `rp_`opt'' if `refsample', meanonly
                    if r(N) local scale=max(abs(r(min)),abs(r(max)))
                }
                tempvar cp
                quietly predict double `cp', `opt'
                quietly count if missing(`rp_`opt'')!=missing(`cp') | ///
                    (abs(`cp'-`rp_`opt'')>=1e-6*(1+abs(`rp_`opt'')) & !missing(`cp'))
                if r(N) noi di "MISMATCH `name': prediction `opt' in " r(N) " rows"
                assert missing(`rp_`opt'')==missing(`cp')
                assert `cp'==`rp_`opt'' | abs(`cp'-`rp_`opt'') <= 1e-7*max(abs(`cp'),abs(`rp_`opt''),`scale') if !missing(`cp')
            }
            estimates drop native_restore
        }
    }
    local rc=_rc
    restore
    if `rc' test_fail "native `name'" "r(`rc')"
    else test_pass "native `name'"
end

noi di "--- Native parity: separation, absorption, and prediction ---"
forvalues seed=1/10 {
    clear
    set seed `seed'
    set obs 1000
    gen g=runiformint(1,20)
    gen h=runiformint(1,10)
    gen double x=rnormal()
    gen double z=rnormal()
    gen double off=rnormal()*.3
    gen double y=rpoisson(exp(.5*x-.3*z+.03*g+off))
    gen fw=runiformint(1,4)
    gen double pw=exp(rnormal())
    gen byte sep=(_n<=100)
    replace y=0 if sep
    gen double combo=x+sep
    foreach wt in "" "[fw=fw]" "[pw=pw]" {
        ppml_native_compare y x z combo `wt', absorb(g h) vce(cluster g h) offset(off) ///
            name("separating combination seed `seed' `wt'") predict
    }
    drop combo sep
    ppml_native_compare y x z, absorb(g#h) offset(off) name("categorical interaction `seed'") predict
    ppml_native_compare y x, absorb(g##c.z h) vce(cluster g) offset(off) name("heterogeneous slopes `seed'") predict
    ppml_native_compare y x, absorb(g#c.z h) offset(off) name("slope-only component `seed'") predict
    ppml_native_compare y, absorb(g h) offset(off) name("pure fixed effects `seed'") predict
    ppml_native_compare y x z, offset(off) name("no fixed effects `seed'") predict
    ppml_native_compare y, offset(off) name("constant only `seed'") predict
}

* The FE design itself separates one off-diagonal cell. No individual level
* is all-zero, so an all-zero-FE screen cannot find this case.
clear
set seed 2201
set obs 240
gen g=1+(_n>80)
gen h=1+(_n>160)
gen double x=rnormal()
gen double y=cond(inrange(_n,81,160),0,1+runiformint(0,3))
ppml_native_compare y x, absorb(g h) name("separation from FE combination") predict
ppml_native_compare y, absorb(g h) name("pure-FE separation") predict

* Base levels, interaction omissions, disconnected FE graphs, and missing rows.
clear
set seed 9121
set obs 800
gen g=runiformint(1,12)
gen h=runiformint(1,8)+20*(g>6)
gen t=runiformint(1,5)+10*(g>6)
gen double x=rnormal()
gen double z=rnormal()
gen byte cat=runiformint(1,3)
gen double y=rpoisson(exp(.4*x+.2*z+.03*g))
replace y=0 if g==1
replace x=. in 101/105
ppml_native_compare y x z, absorb(g h t) name("disconnected three-way DOF") predict
ppml_native_compare y x z, absorb(g h t) vce(cluster g) name("nested and connected DOF") predict
ppml_native_compare y i.cat##c.x z, absorb(g h) name("factor bases and interactions") predict
ppml_native_compare y x z, absorb(g h) vce(cluster g#h) name("interaction cluster") predict
ppml_native_compare y x z, absorb(g h) extra("keepsingletons") name("keepsingletons") predict

* Exposure must be logged exactly once, including after estimates restore.
gen double exposure=exp(z)
ppml_native_compare y x, absorb(g h) extra("exposure(exposure)") name("exposure prediction") predict
ppml_native_compare y x, extra("exposure(exposure)") name("no-FE exposure prediction") predict
* The installed reference's preconditioner errors for multiple slopes in a
* single factor. Its equivalent explicit slope design remains a valid oracle.
capture noisily {
    tempvar rd cd rs rm cm
    quietly ppmlhdfe y x ibn.g#c.z ibn.g#c.exposure, absorb(g h) d(`rd') ///
        tol(1e-14) use_exact_partial(1) use_exact_solver(1) min_ok(3)
    scalar slope_b=_b[x]
    scalar slope_se=_se[x]
    gen byte `rs'=e(sample)
    quietly predict double `rm', mu
    quietly cppmlhdfe y x, absorb(g##c.(z exposure) h) d(`cd') tol(1e-14)
    assert e(sample)==`rs'
    assert reldif(_b[x],slope_b)<1e-7
    assert abs(_se[x]-slope_se)<1e-6*slope_se+1e-12
    quietly predict double `cm', mu
    assert abs(`cm'-`rm')<1e-6*(1+abs(`rm')) if e(sample)
}
if _rc test_fail "native multiple heterogeneous slopes" "r(`=_rc')"
else test_pass "native multiple heterogeneous slopes"

foreach dof in none firstpair pairwise clusters continuous "pairwise clusters" "firstpair continuous" all {
    ppml_native_compare y x z, absorb(g h t) vce(cluster g) extra("dof(`dof')") ///
        name("DOF adjustment `dof'")
}
capture noisily {
    foreach option in "dof(bogus)" "guess(bogus)" "relu_maxiter(0)" {
        capture cppmlhdfe y x, absorb(g h) `option'
        assert _rc==198
    }
    local scratch : all scalars
    foreach name of local scratch {
        assert substr("`name'",1,12)!="__cppmlhdfe_"
    }
}
if _rc test_fail "native invalid options and scratch cleanup" "r(`=_rc')"
else test_pass "native invalid options and scratch cleanup"

* vce() as ppmlhdfe's ms_parse_vce reads it: ols, un and r prefixes and
* conventional; iid is rejected by both
clear
set seed 9122
set obs 600
gen g=runiformint(1,10)
gen h=runiformint(1,6)
gen double x=rnormal()
gen double z=rnormal()
gen double y=rpoisson(exp(.3*x-.2*z+.05*g))
foreach v in ols conventional un unadj rob r {
    ppml_native_compare y x z, absorb(g h) vce(`v') name("vce(`v')")
}
test_error_match, stata_cmd(ppmlhdfe y x z, absorb(g h) vce(iid)) ///
    ctools_cmd(cppmlhdfe y x z, absorb(g h) vce(iid)) testname("vce(iid) rejected")
* A variable named verbose is data, never the legacy bare verbose option
gen double verbose=rnormal()
ppml_native_compare y verbose x, absorb(g h) name("regressor named verbose")
ppml_native_compare y x if verbose > 0, absorb(g h) name("if condition on verbose")
replace verbose=runiformint(1,4)
ppml_native_compare y x z, absorb(g verbose h) name("absorb() term named verbose")
drop verbose
capture quietly cppmlhdfe y x z, absorb(g h) verbose
local rc=_rc
if `rc'==0 & `"`e(cmdline)'"'=="cppmlhdfe y x z, absorb(g h) verbose" test_pass "bare verbose option; e(cmdline) as typed"
else test_fail "bare verbose option; e(cmdline) as typed" `"rc=`rc' cmdline=`e(cmdline)'"'
* guess(variable): a missing value in the estimation sample is reported
gen double m=1+y
replace m=. in 5
capture cppmlhdfe y x z, absorb(g h) guess(variable m)
local rc=_rc
local scratch : all scalars "__cppmlhdfe_*"
if `rc'==198 & "`scratch'"=="" test_pass "guess() variable with a missing value rejected"
else test_fail "guess() variable with a missing value rejected" "rc=`rc' scratch=`scratch'"
quietly cppmlhdfe y x z if _n!=5, absorb(g h) guess(simple) tol(1e-12)
matrix guess_ref=e(b)
matrix guess_V=e(V)
capture quietly cppmlhdfe y x z if _n!=5, absorb(g h) guess(variable m) tol(1e-12)
local rc=_rc
if `rc'==0 {
    matrix guess_got=e(b)
    ppml_min_sigfigs guess_got guess_ref b guess_ref guess_V
    if r(min_sigfigs)>=$DEFAULT_SIGFIGS test_pass "guess() variable missing only outside the sample"
    else test_fail "guess() variable missing only outside the sample" "sigfigs=`r(min_sigfigs)'"
}
else test_fail "guess() variable missing only outside the sample" "rc=`rc'"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_compatibility_cases"
    exit
}

if "`__case'" == "cppmlhdfe_precision_cases" {
* At default tolerance compare SEs and scale covariance errors by marginal SEs.
* Tight fits below also check relative precision of every covariance entry.
capture program drop ppml_default_check
program define ppml_default_check
    syntax varlist [fw pw], absorb(varlist) offset(varname) vce(string) name(string)
    local wt
    if "`weight'"!="" local wt [`weight'`exp']
    preserve
    capture noisily {
        quietly ppmlhdfe `varlist' `wt', absorb(`absorb') offset(`offset') vce(`vce')
        matrix default_b=e(b)
        matrix default_V=e(V)
        tempvar reference_sample
        gen byte `reference_sample'=e(sample)
        quietly cppmlhdfe `varlist' `wt', absorb(`absorb') offset(`offset') vce(`vce')
        assert `reference_sample'==e(sample)
        matrix actual_b=e(b)
        matrix actual_V=e(V)
        forvalues i=1/2 {
            scalar relative_se=abs(sqrt(actual_V[`i',`i']/default_V[`i',`i'])-1)
            scalar PPML_MAX_DEFAULT_SE=max(PPML_MAX_DEFAULT_SE,relative_se)
            assert relative_se < 1e-5
            assert abs(actual_b[1,`i']-default_b[1,`i']) < 1e-5*max(abs(default_b[1,`i']),1e-3)
            forvalues j=1/2 {
                assert abs(actual_V[`i',`j']-default_V[`i',`j']) < ///
                    1e-5*sqrt(default_V[`i',`i']*default_V[`j',`j'])
            }
        }
    }
    local rc=_rc
    restore
    if `rc' test_fail "default `name'" "r(`rc')"
    else test_pass "default `name'"
end
scalar PPML_MAX_DEFAULT_SE=0

* Invoked by validate_cppmlhdfe.do after loading the common test helpers.
* Tight reference fits distinguish estimator errors from stopping-rule error.
noi di as text "--- Precision: randomized weights, VCE, and FE dimensions ---"
forvalues seed = 1/40 {
    clear
    set seed `seed'
    set obs 600
    gen g = runiformint(1,20)
    gen h = runiformint(1,15)
    gen t = runiformint(1,8)
    gen c = runiformint(1,25)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double off = rnormal()*1.5
    gen double y = rpoisson(exp(-1+.7*x-.3*z+.03*g-.02*h+off))
    gen fw = runiformint(1,10)
    gen double pw = exp(rnormal())
    local fes g
    if mod(`seed',3)>0 local fes g h
    if mod(`seed',3)>1 local fes g h t
    foreach w in none fw pw {
        local wt
        if "`w'"!="none" local wt [`w'=`w']
        foreach vc in robust "cluster g" "cluster c" {
            ppml_default_check y x z `wt', absorb(`fes') offset(off) ///
                vce(`vc') name("seed `seed' `w' `vc'")
            benchmark_ppmlhdfe y x z `wt', absorb(`fes') offset(off) ///
                vce(`vc') tolerance(1e-14) ///
                referenceopts("use_exact_partial(1) use_exact_solver(1) min_ok(2)") ///
                testname("precision seed `seed' `w' `vc'")
        }
    }
}

noi di as text "--- Precision: multiway clustering and sample alignment ---"
forvalues seed = 1/10 {
    clear
    set seed `seed'
    set obs 700
    gen g = runiformint(1,25)
    gen h = runiformint(1,17)
    gen t = runiformint(1,9)
    gen c = runiformint(1,40)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double off = rnormal()*.4
    gen double y = rpoisson(exp(.5*x-.3*z+.04*g+off))
    gen fw = runiformint(1,5)
    gen double pw = exp(rnormal())
    * Drop blocks and missing cluster values, not just a prefix of the data.
    replace y = 0 if g==3 | g==9
    replace c = . in 351/360
    gen str8 sh = string(h)
    gen str8 sg = string(g)
    * Multiway covariances combine cluster meats with opposite signs, so a
    * small entry inherits the absolute error of the fixed-effect solve. The
    * reference's conjugate-gradient partial stops at its 1e-12 inner floor
    * with ~1e-11 FE error (seen at the 6th digit of V[x,z] for seed 6 pw
    * "cluster c sh"); its unaccelerated projections converge to the same
    * point as the native solve to ~12 digits, so use them for the oracle.
    foreach w in none fw pw {
        local wt
        if "`w'"!="none" local wt [`w'=`w']
        foreach vc in "cluster g h" "cluster c sh" "cluster g h t" {
            benchmark_ppmlhdfe y x z `wt', absorb(sg t) offset(off) ///
                vce(`vc') tolerance(1e-14) ///
                referenceopts("use_exact_partial(1) use_exact_solver(1) min_ok(2) accel(none)") ///
                testname("multiway seed `seed' `w' `vc'")
        }
    }
}

noi di as text "--- Precision: weight scaling and inference ---"
capture noisily {
    clear
    set seed 99
    set obs 1000
    gen g=mod(_n,20)
    gen h=mod(_n,7)
    gen double x=rnormal()
    gen double z=rnormal()
    gen double y=rpoisson(exp(.5*x-.3*z+.1*g))
    gen double pw=exp(rnormal())
    gen double small=pw*1e-12
    gen double large=pw*1e12
    foreach vc in robust "cluster g" "cluster g h" {
        ppmlhdfe y x z [pw=pw], absorb(g h) vce(`vc') tolerance(1e-14)
        matrix expected_b=e(b)
        matrix expected_V=e(V)
        scalar expected_chi2=e(chi2)
        foreach w in pw small large {
            cppmlhdfe y x z [pw=`w'], absorb(g h) vce(`vc') tolerance(1e-14)
            matrix actual_b=e(b)
            matrix actual_V=e(V)
            assert mreldif(actual_b,expected_b) < 1e-8
            assert mreldif(actual_V,expected_V) < 1e-9
            assert reldif(e(chi2),expected_chi2) < 1e-7
            assert missing(e(df_r))
            test x z
            assert reldif(r(chi2),e(chi2)) < 1e-10
        }
    }
}
if _rc test_fail "PPML weight-scale invariance and Wald inference" "r(`=_rc')"
else test_pass "PPML weight-scale invariance and Wald inference"

capture noisily {
    * The reference may itself omit heavily shifted regressors. Check the
    * estimator's shift invariance against the well-conditioned reference.
    ppmlhdfe y x z, absorb(g) tolerance(1e-14)
    matrix expected_b=e(b)
    matrix expected_V=e(V)
    gen double shifted_x=x+1e6
    gen double shifted_z=z+1e6
    cppmlhdfe y shifted_x shifted_z, absorb(g) tolerance(1e-14)
    matrix shifted_b=e(b)
    matrix shifted_V=e(V)
    matrix shifted_b=shifted_b[1,1..2]
    matrix shifted_V=shifted_V[1..2,1..2]
    matrix expected_b=expected_b[1,1..2]
    matrix expected_V=expected_V[1..2,1..2]
    assert mreldif(shifted_b,expected_b) < 1e-7
    assert mreldif(shifted_V,expected_V) < 1e-9
}
if _rc test_fail "PPML large regressor means" "r(`=_rc')"
else test_pass "PPML large regressor means"

capture noisily {
    gen one=1
    capture cppmlhdfe y x, absorb(g) vce(cluster g one)
    assert _rc==459
    capture cppmlhdfe y x, absorb(g, bogus_option)
    assert _rc==198
    capture cppmlhdfe y x, absorb(g) offset(pw) exposure(pw)
    assert _rc==198
}
if _rc test_fail "PPML invalid options fail explicitly" "r(`=_rc')"
else test_pass "PPML invalid options fail explicitly"

noi di "PPML_MAX_DEFAULT_RELATIVE_SE_ERROR=" %12.5g PPML_MAX_DEFAULT_SE
capture noisily {
    * Scaling the response changes the absorbed intercept only.
    foreach vc in robust "cluster g" "cluster g h" {
        ppmlhdfe y x z, absorb(g h) vce(`vc') tolerance(1e-14)
        matrix expected_b=e(b)
        matrix expected_V=e(V)
        foreach scale in 1e-10 1e10 {
            tempvar scaled_y
            gen double `scaled_y'=y*`scale'
            cppmlhdfe `scaled_y' x z, absorb(g h) vce(`vc') tolerance(1e-14)
            matrix scaled_b=e(b)
            matrix scaled_b[1,colsof(scaled_b)]=scaled_b[1,colsof(scaled_b)]-ln(`scale')
            assert mreldif(scaled_b,expected_b)<1e-8
            assert mreldif(e(V),expected_V)<1e-9
        }
    }
}
if _rc test_fail "PPML response-scale invariance" "r(`=_rc')"
else test_pass "PPML response-scale invariance"

* General separation must select the same sample and refit in C.
capture noisily {
    gen byte separates=(y==0 & x<0)
    ppmlhdfe y x z separates, absorb(g h) tol(1e-12)
    gen byte reference_sample=e(sample)
    matrix reference_b=e(b)
    matrix reference_V=e(V)
    cppmlhdfe y x z separates, absorb(g h) tol(1e-12)
    assert e(sample)==reference_sample
    assert mreldif(e(b),reference_b)<1e-7
    assert mreldif(e(V),reference_V)<1e-7
}
if _rc test_fail "PPML native general separation" "r(`=_rc')"
else test_pass "PPML native general separation"
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_precision_cases"
    exit
}

if "`__case'" == "cppmlhdfe_upstream_cases" {
* Numerical examples adapted from the upstream ppmlhdfe test suite.
* Source: https://github.com/sergiocorreia/ppmlhdfe/tree/master/test

clear
input byte(y id1 id2)
0 1 1
1 1 1
0 2 1
0 2 2
1 2 2
end
ppml_native_compare y, absorb(id1 id2) name("upstream hard1")

clear
input byte FTA int(ikt_fe jkt_fe ijk_fe)
1 1 3 1
1 2 4 1
1 1 5 2
1 2 6 2
1 3 1 3
1 4 2 3
1 3 7 4
1 4 8 4
1 5 1 5
1 6 2 5
0 5 3 6
1 6 4 6
0 5 5 7
1 6 6 7
1 5 7 8
1 6 8 8
end
ppml_native_compare FTA, absorb(ijk_fe jkt_fe ikt_fe) extra("keepsingletons") name("upstream hard2")

clear
input byte(y x1 x2 x3 x4)
0  5 11  2  2
0  5  2 11  2
0  5  2  2 11
0  0 -1 -1 -1
1  5  5  5  5
2  4  4  4  4
3  3  3  3  3
4  2  2  2  2
5  1  1  1  1
end
ppml_native_compare y x1 x2 x3 x4, name("upstream hard4")

clear
input byte(y x1 x2 x3 x4 x5)
0  0  0  0  0 1
0  0  0  0  0 2
0  5 11  2  2 0
0  5  2 11  2 0
0  5  2  2 11 0
0  0 -1 -1 -1 0
1  5  5  5  5 0
2  4  4  4  4 0
3  3  3  3  3 0
4  2  2  2  2 0
5  1  1  1  1 0
end
ppml_native_compare y x1 x2 x3 x4 x5, name("upstream hard5")

clear
input double trade byte(fta exporter trend c)
35412.88183544995 0  8 15 1
619124.7429669127 0  8 16 1
 645947.355868591 0 35  0 1
2113023.937040806 1 35 12 1
end
ppml_native_compare trade fta, absorb(exporter##c.trend c) referenceopts("guess(variable trade)") name("upstream slopes1")

clear
input double(id1 y) byte(x z c)
1   67793.6092504498 0 -1 1
1              68774 0  0 1
2 29460066.818006903 1 -1 1
2  40370377.65078497 1 10 1
3  4215.552981376648 0  7 1
3  4288.904005765915 0 15 1
4  483882.9337809086 0  7 1
4 306164.57000368834 0  8 1
5 29557.571687340736 0  6 1
5  170794.3896008134 0 12 1
6 1963894.9483253956 0  6 1
6 3364687.2930059433 0 12 1
end
ppml_native_compare y x, absorb(id1##c.z c) name("upstream slopes2")

clear
input byte(index id1 id2 x z) double y
	1 1 1 0 -1  3369071263
	2 1 1 0  2  2535746605
	3 1 1 1 15 15241548895
	4 1 2 0 11        1040
	5 1 2 0 13         166
	6 1 2 0 15    21682500
	7 1 3 0 -1           0
	8 1 3 0  7      764800
end
* At tight tolerance, Stata rejects the upstream covariance as nonsymmetric
* or highly singular and posts zeros. Use an analytic sandwich oracle. The rows
* in id2=1 are exactly fitted (so the slope variance is zero). In id2=2 the
* fitted means form (a, a*r, a*r^2); solve the count and trend score equations.
* The normalized constant's variance is 6/5 * sum((y-mu)^2) / sum(y)^2.
scalar slopes3_total = y[4]+y[5]+y[6]
scalar slopes3_m = (y[5]+2*y[6])/slopes3_total
scalar slopes3_r = (slopes3_m-1 + sqrt((1-slopes3_m)^2 + ///
    4*(2-slopes3_m)*slopes3_m))/(2*(2-slopes3_m))
scalar slopes3_a = slopes3_total/(1+slopes3_r+slopes3_r^2)
scalar slopes3_v = (6/5)*((y[4]-slopes3_a)^2 + ///
    (y[5]-slopes3_a*slopes3_r)^2 + (y[6]-slopes3_a*slopes3_r^2)^2) / ///
    (y[1]+y[2]+y[3]+slopes3_total)^2
matrix slopes3_V = (0,0 \ 0,slopes3_v)
ppml_native_compare y x, absorb(id1 id2##c.z) name("upstream slopes3") vceref(slopes3_V)
matrix drop slopes3_V
scalar drop slopes3_total slopes3_m slopes3_r slopes3_a slopes3_v

clear
input long y byte(x z i j)
	      0 0 -1  1  1
	   5211 0 -1  1  4
	  10498 0  2  2  4
	   1310 0  2  2  5
	    189 0  3  3  2
	      0 0  3  3  6
	   1468 0  6  4  5
	     34 0  6  4  7
	  12000 0  7  5  3
	   1036 0  7  5  5
	      0 0  8  6  3
	   1247 0  8  6  7
	    460 0  9  7  1
	     70 0  9  7  7
	     28 0 14  8  2
	      4 0 14  8  3
	     17 0 14  8  6
	      0 0 15  9  1
	      8 0 15  9  6
	 203886 0  3 10  8
	3840369 1  3 10  9
	3961747 1  8 11  9
	 710771 0  8 11 10
	 810565 0 13 12  8
	 249014 0 13 12 10
	 819845 0 16 13  8
	5589636 1 16 13  9
end
* Match the reference's default stopping rule and retained rows, including
* the ill-conditioned separation direction that it treats as numerical zero.
ppml_native_compare y x, absorb(i j##c.z) name("upstream slopes4 default sample")
ppml_native_compare y x, absorb(i j##c.z) extra("relu_accelerate(1)") name("upstream slopes4 accelerated ReLU")
capture noisily cppmlhdfe y x, absorb(i j##c.z) separation(relu) relu_maxiter(1) relu_strict(1)
if _rc==9010 test_pass "native strict ReLU iteration limit"
else test_fail "native strict ReLU iteration limit" "expected r(9010), got r(`=_rc')"

clear
set obs 6
gen y=max(_n-3,0)
gen double x1=100*(_n<3)+(_n==6)
ppml_native_compare y x1, name("upstream ill_conditioned")
clear
drawnorm x1, n(1000) seed(101010) double
gen double u=rpoisson(1)
gen double y=exp(1+10*x1)*u
gen double x2=(y==0)
ppml_native_compare y x1 x2, name("upstream collinear2")
clear
drawnorm u x1 x2, n(1000) seed(101010) double
gen double y=exp(40+x1+x2+u)
ppml_native_compare y x1 x2, name("upstream large_y")
    global CTOOLS_CASE_COMPLETE "cppmlhdfe_upstream_cases"
    exit
}

if "`__case'" == "cppmlhdfe_publication" {
do validation/support/validate_setup.do
clear
set seed 2411
set obs 200
generate double x = rnormal()
generate g = ceil(_n/10)
generate double y = rpoisson(exp(.4 + .3*x + .02*g))
generate double _ppmlhdfe_d = 777 + _n
generate double __hdfe1__ = -_n
label variable _ppmlhdfe_d "previous fixed effects"
format _ppmlhdfe_d %20.8f
label variable __hdfe1__ "previous individual effects"
sort x
tempfile original before_failure
save `original'
foreach failure in negative absorb vce empty nonconvergence publish {
    use `original', clear
    local opts absorb(g, savefe) d2
    local qualifier
    if "`failure'" == "negative" replace y = -1 in 1
    if "`failure'" == "absorb" local opts absorb(no_such_effect) d2
    if "`failure'" == "vce" local opts `opts' vce(no_such_vce)
    if "`failure'" == "empty" local qualifier if 0
    if "`failure'" == "nonconvergence" local opts `opts' irlsmaxiter(1)
    if "`failure'" == "publish" global CTOOLS_cppmlhdfe_fault publish
    unab names_before : _all
    local sort_before : sortedby
    save `before_failure', replace
    capture noisily cppmlhdfe y x `qualifier', `opts'
    local failed_rc = _rc
    macro drop CTOOLS_cppmlhdfe_fault
    assert `failed_rc' != 0
    if "`failure'" == "publish" assert `failed_rc' == 459
    cf _all using `before_failure', all
    unab names_after : _all
    local sort_after : sortedby
    assert "`names_before'" == "`names_after'"
    assert "`sort_before'" == "`sort_after'"
    assert "`: variable label _ppmlhdfe_d'" == "previous fixed effects"
    assert "`: format _ppmlhdfe_d'" == "%20.8f"
    assert "`: type _ppmlhdfe_d'" == "double"
    assert "`: variable label __hdfe1__'" == "previous individual effects"
    test_pass "failed `failure' preserves automatic PPML outputs"
}
use `original', clear
quietly cppmlhdfe y x, absorb(g, savefe) d2
assert "`e(d)'" == "_ppmlhdfe_d"
assert _ppmlhdfe_d != 777 + _n if e(sample)
quietly predict double fitted, mu
assert fitted > 0 & !missing(fitted) if e(sample)
quietly cppmlhdfe y x, absorb(g, savefe) d2
test_pass "successful reruns replace automatic PPML outputs"
* An interruption after posting (a Break during the display, simulated by
* the display fault) leaves the new e() with the new automatic outputs
quietly cppmlhdfe y x, absorb(g, savefe) d2
generate double first_d = _ppmlhdfe_d
global CTOOLS_cppmlhdfe_fault display
capture cppmlhdfe y x if _n > 20, absorb(g, savefe) d2
local display_rc = _rc
macro drop CTOOLS_cppmlhdfe_fault
quietly count if _ppmlhdfe_d != first_d & !missing(_ppmlhdfe_d)
local display_changed = r(N)
quietly count if !missing(_ppmlhdfe_d)
if `display_rc' == 1 & e(N) == 180 & `display_changed' == 180 & r(N) == 180 ///
    test_pass "interrupted display: new e() with the new automatic outputs"
else test_fail "interrupted display: new e() with the new automatic outputs" ///
    "rc=`display_rc' N=`=e(N)' changed=`display_changed'"
print_summary "cppmlhdfe output publication"
global CTOOLS_COMPONENT_COMPLETE "cppmlhdfe_publication"
global CTOOLS_CASE_COMPLETE "cppmlhdfe_publication"
exit
}

di as error "Unknown cppmlhdfe case: `__case'"
exit 198
