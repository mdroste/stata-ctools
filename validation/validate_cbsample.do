/*******************************************************************************
 * validate_cbsample.do
 *
 * Validation tests for cbsample command
 * Tests bootstrap sampling with replacement functionality
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text ""
noi di as text "======================================================================"
noi di as text "              CBSAMPLE VALIDATION TEST SUITE"
noi di as text "======================================================================"

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
noi print_section "Plugin Check"

clear
set obs 100
gen id = _n
set seed 12345
capture noisily cbsample
if _rc != 0 {
    noi test_fail "cbsample plugin load" "plugin returned error `=_rc'"
    noi print_summary "cbsample"
    exit 1
}
noi test_pass "cbsample plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic bootstrap sample
 ******************************************************************************/
noi print_section "Basic Bootstrap Sample"

* Default n (same as original)
clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N
cbsample
local N_after = _N

* With bootstrap, some obs get dropped (sampled 0 times)
if `N_after' > 0 & `N_after' <= `N_before' {
    noi test_pass "basic bootstrap sample (N=`N_after' from `N_before')"
}
else {
    noi test_fail "basic bootstrap sample" "unexpected N=`N_after'"
}

* Explicit n (positional argument, like bsample)
clear
set obs 100
gen id = _n
set seed 12345
cbsample 50
if _N > 0 & _N <= 100 {
    noi test_pass "bootstrap with n=50"
}
else {
    noi test_fail "bootstrap with n=50" "N=`=_N'"
}

/*******************************************************************************
 * SECTION 3: Weight variable (matches bsample syntax: requires existing var)
 ******************************************************************************/
noi print_section "Weight Variable"

clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
cbsample, weight(bsweight)

* Should keep all obs and have weight variable populated
if _N == 100 {
    noi test_pass "weight option keeps all obs"
}
else {
    noi test_fail "weight option keeps all obs" "N=`=_N'"
}

* Check weight variable was populated (not all missing)
count if bsweight != .
if r(N) > 0 {
    noi test_pass "weight variable populated"
}
else {
    noi test_fail "weight variable populated" "all values missing"
}

* Check total weight sums to ~100
summarize bsweight, meanonly
local total_weight = r(sum)
if abs(`total_weight' - 100) < 10 {
    noi test_pass "total weight ~100 (got `total_weight')"
}
else {
    noi test_fail "total weight ~100" "got `total_weight'"
}

/*******************************************************************************
 * SECTION 4: Stratified bootstrap
 ******************************************************************************/
noi print_section "Stratified Bootstrap"

clear
set obs 1000
gen stratum = mod(_n, 10)
gen id = _n
gen bsweight = .
set seed 12345
cbsample, strata(stratum) weight(bsweight)

* Check each stratum has some weight
local all_ok = 1
forvalues s = 0/9 {
    summarize bsweight if stratum == `s', meanonly
    if r(sum) == 0 {
        local all_ok = 0
    }
}
if `all_ok' {
    noi test_pass "stratified bootstrap has weight in all strata"
}
else {
    noi test_fail "stratified bootstrap" "some strata have zero weight"
}

/*******************************************************************************
 * SECTION 5: Cluster bootstrap
 ******************************************************************************/
noi print_section "Cluster Bootstrap"

clear
set obs 1000
gen cluster_id = ceil(_n / 10)  /* 100 clusters of 10 each */
gen id = _n
gen bsweight = .
set seed 12345
cbsample, cluster(cluster_id) weight(bsweight)

* Within each cluster, all obs should have same weight
local all_ok = 1
levelsof cluster_id, local(clusters)
foreach cl of local clusters {
    summarize bsweight if cluster_id == `cl'
    if r(sd) > 0.001 & r(sd) != . {
        local all_ok = 0
    }
}
if `all_ok' {
    noi test_pass "cluster bootstrap: uniform weights within clusters"
}
else {
    noi test_fail "cluster bootstrap" "weights vary within some clusters"
}

/*******************************************************************************
 * SECTION 6: Reproducibility with set seed
 ******************************************************************************/
noi print_section "Reproducibility"

clear
set obs 100
gen id = _n
gen w1 = .
set seed 99999
cbsample, weight(w1)

tempfile sample1
save `sample1'

clear
set obs 100
gen id = _n
gen w2 = .
set seed 99999
cbsample, weight(w2)

* Merge and compare
merge 1:1 id using `sample1'
count if abs(w1 - w2) > 0.001
local n_diff = r(N)

if `n_diff' == 0 {
    noi test_pass "same seed yields identical weights"
}
else {
    noi test_fail "same seed yields identical weights" "`n_diff' weights differ"
}

/*******************************************************************************
 * SECTION 7: Edge cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Small dataset
clear
set obs 10
gen id = _n
gen bsweight = .
set seed 12345
cbsample, weight(bsweight)
summarize bsweight, meanonly
local total = r(sum)
if abs(`total' - 10) < 3 {
    noi test_pass "small dataset bootstrap (total weight ~10)"
}
else {
    noi test_fail "small dataset bootstrap" "total weight=`total'"
}

* Large n (positional argument, like bsample)
clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
cbsample 200, weight(bsweight)
summarize bsweight, meanonly
local total = r(sum)
if abs(`total' - 200) < 20 {
    noi test_pass "n > nobs bootstrap (total weight ~200)"
}
else {
    noi test_fail "n > nobs bootstrap" "total weight=`total'"
}

} /* end quietly */

/*******************************************************************************
 * Summary
 ******************************************************************************/
print_summary "cbsample"

if $TESTS_FAILED > 0 {
    exit 1
}
