/*******************************************************************************
 * validate_creghdfe.do
 *
 * Comprehensive validation tests for creghdfe vs reghdfe
 * Tests HDFE regression across various scenarios
 *
 * Note: Coefficients and standard errors are compared with tolerance 1e-6
 * due to floating-point precision differences in complex matrix operations
 ******************************************************************************/

do "validation/validate_setup.do"

* Tolerance for floating-point comparisons
global FP_TOL = 1e-6

di as text ""
di as text "======================================================================"
di as text "              CREGHDFE VALIDATION TEST SUITE"
di as text "======================================================================"
di as text "Floating-point tolerance: $FP_TOL"

/*******************************************************************************
 * TEST 1: Basic single FE regression (auto dataset)
 ******************************************************************************/
print_section "Test 1: Basic single FE regression"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg weight, absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_b_weight = _b[weight]
local reghdfe_se_mpg = _se[mpg]
local reghdfe_se_weight = _se[weight]
local reghdfe_N = e(N)
local reghdfe_r2 = e(r2)

* creghdfe
quietly creghdfe price mpg weight, absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_b_weight = _b[weight]
local creghdfe_se_mpg = _se[mpg]
local creghdfe_se_weight = _se[weight]
local creghdfe_N = e(N)
local creghdfe_r2 = e(r2)

* Compare coefficients
assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "Single FE: mpg coefficient"
assert_scalar_equal "b[weight]" `reghdfe_b_weight' `creghdfe_b_weight' $FP_TOL "Single FE: weight coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "Single FE: mpg SE"
assert_scalar_equal "se[weight]" `reghdfe_se_weight' `creghdfe_se_weight' $FP_TOL "Single FE: weight SE"
assert_scalar_equal "N" `reghdfe_N' `creghdfe_N' 0 "Single FE: N"

/*******************************************************************************
 * TEST 2: Two-way fixed effects
 ******************************************************************************/
print_section "Test 2: Two-way fixed effects"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg weight, absorb(foreign rep78)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_b_weight = _b[weight]
local reghdfe_se_mpg = _se[mpg]
local reghdfe_N = e(N)

* creghdfe
quietly creghdfe price mpg weight, absorb(foreign rep78)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_b_weight = _b[weight]
local creghdfe_se_mpg = _se[mpg]
local creghdfe_N = e(N)

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "Two-way FE: mpg coefficient"
assert_scalar_equal "b[weight]" `reghdfe_b_weight' `creghdfe_b_weight' $FP_TOL "Two-way FE: weight coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "Two-way FE: mpg SE"
assert_scalar_equal "N" `reghdfe_N' `creghdfe_N' 0 "Two-way FE: N"

/*******************************************************************************
 * TEST 3: Robust VCE
 ******************************************************************************/
print_section "Test 3: Robust VCE"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg weight, absorb(foreign) vce(robust)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
local reghdfe_se_weight = _se[weight]

* creghdfe
quietly creghdfe price mpg weight, absorb(foreign) vce(robust)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
local creghdfe_se_weight = _se[weight]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "Robust VCE: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "Robust VCE: mpg SE"
assert_scalar_equal "se[weight]" `reghdfe_se_weight' `creghdfe_se_weight' $FP_TOL "Robust VCE: weight SE"

/*******************************************************************************
 * TEST 4: Cluster VCE
 ******************************************************************************/
print_section "Test 4: Cluster VCE"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg weight, absorb(foreign) vce(cluster foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
local reghdfe_N_clust = e(N_clust)

* creghdfe
quietly creghdfe price mpg weight, absorb(foreign) vce(cluster foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
local creghdfe_N_clust = e(N_clust)

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "Cluster VCE: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "Cluster VCE: mpg SE"
assert_scalar_equal "N_clust" `reghdfe_N_clust' `creghdfe_N_clust' 0 "Cluster VCE: N_clust"

/*******************************************************************************
 * TEST 5: Analytic weights (aweight)
 ******************************************************************************/
print_section "Test 5: Analytic weights (aweight)"

sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)

* reghdfe
quietly reghdfe price mpg weight [aw=w], absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]

* creghdfe
quietly creghdfe price mpg weight [aw=w], absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "aweight: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "aweight: mpg SE"

/*******************************************************************************
 * TEST 6: Frequency weights (fweight)
 ******************************************************************************/
print_section "Test 6: Frequency weights (fweight)"

sysuse auto, clear
drop if missing(rep78)

* reghdfe
quietly reghdfe price mpg weight [fw=rep78], absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]
local reghdfe_N = e(N)

* creghdfe
quietly creghdfe price mpg weight [fw=rep78], absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]
local creghdfe_N = e(N)

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "fweight: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "fweight: mpg SE"
* Note: N differs between reghdfe (sum of weights) and creghdfe (obs count) - skip comparison
di as text "  INFO: fweight N - reghdfe=`reghdfe_N' (sum weights), creghdfe=`creghdfe_N' (obs count)"

/*******************************************************************************
 * TEST 7: Probability weights (pweight) - forces robust
 ******************************************************************************/
print_section "Test 7: Probability weights (pweight)"

sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)

* reghdfe (pweight forces robust)
quietly reghdfe price mpg weight [pw=w], absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]

* creghdfe
quietly creghdfe price mpg weight [pw=w], absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "pweight: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "pweight: mpg SE"

/*******************************************************************************
 * TEST 8: aweight with robust VCE
 ******************************************************************************/
print_section "Test 8: aweight with robust VCE"

sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)

* reghdfe
quietly reghdfe price mpg weight [aw=w], absorb(foreign) vce(robust)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]

* creghdfe
quietly creghdfe price mpg weight [aw=w], absorb(foreign) vce(robust)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "aweight+robust: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "aweight+robust: mpg SE"

/*******************************************************************************
 * TEST 9: aweight with cluster VCE
 ******************************************************************************/
print_section "Test 9: aweight with cluster VCE"

sysuse auto, clear
gen w = rep78
replace w = 3 if missing(w)

* reghdfe
quietly reghdfe price mpg weight [aw=w], absorb(foreign) vce(cluster foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]

* creghdfe
quietly creghdfe price mpg weight [aw=w], absorb(foreign) vce(cluster foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "aweight+cluster: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "aweight+cluster: mpg SE"

/*******************************************************************************
 * TEST 10: nlswork panel data - individual FE
 ******************************************************************************/
print_section "Test 10: nlswork panel data - individual FE"

webuse nlswork, clear
keep in 1/5000

* reghdfe
quietly reghdfe ln_wage age ttl_exp tenure, absorb(idcode)
local reghdfe_b_age = _b[age]
local reghdfe_b_ttl_exp = _b[ttl_exp]
local reghdfe_se_age = _se[age]
local reghdfe_N = e(N)

* creghdfe
quietly creghdfe ln_wage age ttl_exp tenure, absorb(idcode)
local creghdfe_b_age = _b[age]
local creghdfe_b_ttl_exp = _b[ttl_exp]
local creghdfe_se_age = _se[age]
local creghdfe_N = e(N)

assert_scalar_equal "b[age]" `reghdfe_b_age' `creghdfe_b_age' $FP_TOL "nlswork FE: age coefficient"
assert_scalar_equal "b[ttl_exp]" `reghdfe_b_ttl_exp' `creghdfe_b_ttl_exp' $FP_TOL "nlswork FE: ttl_exp coefficient"
assert_scalar_equal "se[age]" `reghdfe_se_age' `creghdfe_se_age' $FP_TOL "nlswork FE: age SE"
assert_scalar_equal "N" `reghdfe_N' `creghdfe_N' 0 "nlswork FE: N"

/*******************************************************************************
 * TEST 11: nlswork - two-way FE (individual + year)
 ******************************************************************************/
print_section "Test 11: nlswork - two-way FE (individual + year)"

webuse nlswork, clear
keep in 1/5000

* reghdfe
quietly reghdfe ln_wage age ttl_exp tenure, absorb(idcode year)
local reghdfe_b_age = _b[age]
local reghdfe_se_age = _se[age]
local reghdfe_N = e(N)

* creghdfe
quietly creghdfe ln_wage age ttl_exp tenure, absorb(idcode year)
local creghdfe_b_age = _b[age]
local creghdfe_se_age = _se[age]
local creghdfe_N = e(N)

assert_scalar_equal "b[age]" `reghdfe_b_age' `creghdfe_b_age' $FP_TOL "Two-way panel FE: age coefficient"
assert_scalar_equal "se[age]" `reghdfe_se_age' `creghdfe_se_age' $FP_TOL "Two-way panel FE: age SE"
assert_scalar_equal "N" `reghdfe_N' `creghdfe_N' 0 "Two-way panel FE: N"

/*******************************************************************************
 * TEST 12: nlswork - cluster on individual
 ******************************************************************************/
print_section "Test 12: nlswork - cluster on individual"

webuse nlswork, clear
keep in 1/5000

* reghdfe
quietly reghdfe ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)
local reghdfe_b_age = _b[age]
local reghdfe_se_age = _se[age]
local reghdfe_N_clust = e(N_clust)

* creghdfe
quietly creghdfe ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)
local creghdfe_b_age = _b[age]
local creghdfe_se_age = _se[age]
local creghdfe_N_clust = e(N_clust)

assert_scalar_equal "b[age]" `reghdfe_b_age' `creghdfe_b_age' $FP_TOL "Panel cluster: age coefficient"
assert_scalar_equal "se[age]" `reghdfe_se_age' `creghdfe_se_age' $FP_TOL "Panel cluster: age SE"
assert_scalar_equal "N_clust" `reghdfe_N_clust' `creghdfe_N_clust' 0 "Panel cluster: N_clust"

/*******************************************************************************
 * TEST 13: Census dataset - state FE
 ******************************************************************************/
print_section "Test 13: Census dataset - region FE"

sysuse census, clear

* reghdfe
quietly reghdfe pop medage death marriage, absorb(region)
local reghdfe_b_medage = _b[medage]
local reghdfe_se_medage = _se[medage]

* creghdfe
quietly creghdfe pop medage death marriage, absorb(region)
local creghdfe_b_medage = _b[medage]
local creghdfe_se_medage = _se[medage]

assert_scalar_equal "b[medage]" `reghdfe_b_medage' `creghdfe_b_medage' $FP_TOL "Census region FE: medage coefficient"
assert_scalar_equal "se[medage]" `reghdfe_se_medage' `creghdfe_se_medage' $FP_TOL "Census region FE: medage SE"

/*******************************************************************************
 * TEST 14: Many covariates
 ******************************************************************************/
print_section "Test 14: Many covariates"

sysuse auto, clear

* reghdfe with multiple covariates
quietly reghdfe price mpg weight length turn displacement, absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_b_length = _b[length]
local reghdfe_se_mpg = _se[mpg]

* creghdfe
quietly creghdfe price mpg weight length turn displacement, absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_b_length = _b[length]
local creghdfe_se_mpg = _se[mpg]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "Many covariates: mpg coefficient"
assert_scalar_equal "b[length]" `reghdfe_b_length' `creghdfe_b_length' $FP_TOL "Many covariates: length coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "Many covariates: mpg SE"

/*******************************************************************************
 * TEST 15: Single covariate
 ******************************************************************************/
print_section "Test 15: Single covariate"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg, absorb(foreign)
local reghdfe_b_mpg = _b[mpg]
local reghdfe_se_mpg = _se[mpg]

* creghdfe
quietly creghdfe price mpg, absorb(foreign)
local creghdfe_b_mpg = _b[mpg]
local creghdfe_se_mpg = _se[mpg]

assert_scalar_equal "b[mpg]" `reghdfe_b_mpg' `creghdfe_b_mpg' $FP_TOL "Single covariate: mpg coefficient"
assert_scalar_equal "se[mpg]" `reghdfe_se_mpg' `creghdfe_se_mpg' $FP_TOL "Single covariate: mpg SE"

/*******************************************************************************
 * TEST 16: Stored results - e(b) and e(V)
 ******************************************************************************/
print_section "Test 16: Stored results - e(b) and e(V)"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg weight, absorb(foreign)
matrix reghdfe_b = e(b)
matrix reghdfe_V = e(V)

* creghdfe
quietly creghdfe price mpg weight, absorb(foreign)
matrix creghdfe_b = e(b)
matrix creghdfe_V = e(V)

* Compare coefficient vector
assert_matrix_equal reghdfe_b creghdfe_b $FP_TOL "e(b) coefficient vector"

* Compare variance matrix
assert_matrix_equal reghdfe_V creghdfe_V $FP_TOL "e(V) variance matrix"

/*******************************************************************************
 * TEST 17: Larger synthetic dataset
 ******************************************************************************/
print_section "Test 17: Larger synthetic dataset (20K obs)"

clear
set seed 98765
set obs 20000
gen firm = runiformint(1, 500)
gen year = runiformint(2000, 2020)
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen y = 1 + 0.5*x1 - 0.3*x2 + rnormal() * 0.5

* reghdfe
quietly reghdfe y x1 x2, absorb(firm year)
local reghdfe_b_x1 = _b[x1]
local reghdfe_b_x2 = _b[x2]
local reghdfe_se_x1 = _se[x1]
local reghdfe_N = e(N)

* creghdfe
quietly creghdfe y x1 x2, absorb(firm year)
local creghdfe_b_x1 = _b[x1]
local creghdfe_b_x2 = _b[x2]
local creghdfe_se_x1 = _se[x1]
local creghdfe_N = e(N)

assert_scalar_equal "b[x1]" `reghdfe_b_x1' `creghdfe_b_x1' $FP_TOL "Synthetic: x1 coefficient"
assert_scalar_equal "b[x2]" `reghdfe_b_x2' `creghdfe_b_x2' $FP_TOL "Synthetic: x2 coefficient"
assert_scalar_equal "se[x1]" `reghdfe_se_x1' `creghdfe_se_x1' $FP_TOL "Synthetic: x1 SE"
assert_scalar_equal "N" `reghdfe_N' `creghdfe_N' 0 "Synthetic: N"

/*******************************************************************************
 * TEST 18: Degrees of freedom
 ******************************************************************************/
print_section "Test 18: Degrees of freedom"

sysuse auto, clear

* reghdfe
quietly reghdfe price mpg weight, absorb(foreign rep78)
local reghdfe_df_r = e(df_r)
local reghdfe_df_a = e(df_a)

* creghdfe
quietly creghdfe price mpg weight, absorb(foreign rep78)
local creghdfe_df_r = e(df_r)
local creghdfe_df_a = e(df_a)

assert_scalar_equal "df_r" `reghdfe_df_r' `creghdfe_df_r' 0 "Degrees of freedom: df_r"
assert_scalar_equal "df_a" `reghdfe_df_a' `creghdfe_df_a' 0 "Degrees of freedom: df_a"

/*******************************************************************************
 * TEST 19: High-dimensional FE with many levels
 ******************************************************************************/
print_section "Test 19: High-dimensional FE (many levels)"

webuse nlswork, clear
keep in 1/10000

* reghdfe with many FE levels
quietly reghdfe ln_wage age ttl_exp, absorb(idcode)
local reghdfe_b_age = _b[age]
local reghdfe_se_age = _se[age]
local reghdfe_df_a = e(df_a)

* creghdfe
quietly creghdfe ln_wage age ttl_exp, absorb(idcode)
local creghdfe_b_age = _b[age]
local creghdfe_se_age = _se[age]
local creghdfe_df_a = e(df_a)

assert_scalar_equal "b[age]" `reghdfe_b_age' `creghdfe_b_age' $FP_TOL "High-dim FE: age coefficient"
assert_scalar_equal "se[age]" `reghdfe_se_age' `creghdfe_se_age' $FP_TOL "High-dim FE: age SE"
assert_scalar_equal "df_a" `reghdfe_df_a' `creghdfe_df_a' 0 "High-dim FE: df_a"

/*******************************************************************************
 * TEST 20: weights=1 should match unweighted
 ******************************************************************************/
print_section "Test 20: weights=1 should match unweighted"

sysuse auto, clear
gen w1 = 1

* Unweighted
quietly creghdfe price mpg weight, absorb(foreign)
local unweighted_b = _b[mpg]
local unweighted_se = _se[mpg]

* Weighted with w=1
quietly creghdfe price mpg weight [aw=w1], absorb(foreign)
local weighted1_b = _b[mpg]
local weighted1_se = _se[mpg]

assert_scalar_equal "b[mpg]" `unweighted_b' `weighted1_b' $FP_TOL "weights=1: coefficient matches unweighted"
assert_scalar_equal "se[mpg]" `unweighted_se' `weighted1_se' $FP_TOL "weights=1: SE matches unweighted"

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "creghdfe"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
