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

noi di as text "Running validation tests for cbsample..."

* Plugin functionality check
clear
set obs 100
gen id = _n
set seed 12345
capture cbsample
if _rc != 0 {
    test_fail "cbsample plugin load" "plugin returned error `=_rc'"
    exit 1
}
test_pass "cbsample plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic bootstrap sample
 ******************************************************************************/
print_section "Basic Bootstrap Sample"

* Default n (same as original)
clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N
cbsample
local N_after = _N

* With bootstrap expand/drop, N should be approximately N_before
if `N_after' > 0 & `N_after' <= 2 * `N_before' {
    test_pass "basic bootstrap sample (N=`N_after' from `N_before')"
}
else {
    test_fail "basic bootstrap sample" "N=`N_after' from `N_before' (unexpected size)"
}

* Explicit n (positional argument, like bsample)
clear
set obs 100
gen id = _n
set seed 12345
cbsample 50
if _N == 50 {
    test_pass "bootstrap with n=50 (N=`=_N')"
}
else {
    test_fail "bootstrap with n=50" "N=`=_N', expected 50"
}

/*******************************************************************************
 * SECTION 2b: Deterministic comparison with bsample
 ******************************************************************************/
print_section "Deterministic Comparison with bsample"

* Compare cbsample to Stata's bsample using the same seed
clear
set obs 100
gen id = _n

set seed 12345
preserve
bsample
local stata_N = _N
restore

clear
set obs 100
gen id = _n

set seed 12345
preserve
cbsample
local ctools_N = _N
restore

* Both should produce same N with same seed
if `stata_N' == `ctools_N' {
    test_pass "cbsample matches bsample with same seed (N=`ctools_N')"
}
else {
    test_fail "cbsample vs bsample" "N differ: bsample=`stata_N' cbsample=`ctools_N'"
}

* Compare with explicit n argument
clear
set obs 100
gen id = _n

set seed 54321
preserve
bsample 50
local stata_N = _N
restore

clear
set obs 100
gen id = _n

set seed 54321
preserve
cbsample 50
local ctools_N = _N
restore

if `stata_N' == `ctools_N' {
    test_pass "cbsample 50 matches bsample 50 with same seed (N=`ctools_N')"
}
else {
    test_fail "cbsample 50 vs bsample 50" "N differ: bsample=`stata_N' cbsample=`ctools_N'"
}

/*******************************************************************************
 * SECTION 3: Weight variable (matches bsample syntax: requires existing var)
 ******************************************************************************/
print_section "Weight Variable"

clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
cbsample, weight(bsweight)

* Should keep all obs and have weight variable populated
if _N == 100 {
    test_pass "weight option keeps all obs"
}
else {
    test_fail "weight option keeps all obs" "N=`=_N'"
}

* Check weight variable was populated (not all missing)
count if bsweight != .
if r(N) > 0 {
    test_pass "weight variable populated"
}
else {
    test_fail "weight variable populated" "all values missing"
}

* Check total weight sums to ~100
summarize bsweight, meanonly
local total_weight = r(sum)
if abs(`total_weight' - 100) < 5 {
    test_pass "total weight ~100 (got `total_weight')"
}
else {
    test_fail "total weight ~100" "got `total_weight'"
}

/*******************************************************************************
 * SECTION 4: Stratified bootstrap
 ******************************************************************************/
print_section "Stratified Bootstrap"

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
    test_pass "stratified bootstrap has weight in all strata"
}
else {
    test_fail "stratified bootstrap" "some strata have zero weight"
}

/*******************************************************************************
 * SECTION 5: Cluster bootstrap
 ******************************************************************************/
print_section "Cluster Bootstrap"

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
    test_pass "cluster bootstrap: uniform weights within clusters"
}
else {
    test_fail "cluster bootstrap" "weights vary within some clusters"
}

/*******************************************************************************
 * SECTION 6: Reproducibility with set seed
 ******************************************************************************/
print_section "Reproducibility"

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

* Merge and compare (deterministic: same seed must produce exact same weights)
merge 1:1 id using `sample1'
count if w1 != w2
local n_diff = r(N)

if `n_diff' == 0 {
    test_pass "same seed yields identical weights"
}
else {
    test_fail "same seed yields identical weights" "`n_diff' weights differ"
}

/*******************************************************************************
 * SECTION 7: Edge cases
 ******************************************************************************/
print_section "Edge Cases"

* Small dataset
clear
set obs 10
gen id = _n
gen bsweight = .
set seed 12345
cbsample, weight(bsweight)
summarize bsweight, meanonly
local total = r(sum)
* Small N means more variance; use 2-unit tolerance (~20% for N=10)
if abs(`total' - 10) < 2 {
    test_pass "small dataset bootstrap (total weight ~10)"
}
else {
    test_fail "small dataset bootstrap" "total weight=`total'"
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
if abs(`total' - 200) < 10 {
    test_pass "n > nobs bootstrap (total weight ~200)"
}
else {
    test_fail "n > nobs bootstrap" "total weight=`total'"
}

/*******************************************************************************
 * SECTION 8: If/In Conditions
 ******************************************************************************/
print_section "If/In Conditions"

* if condition test
clear
set obs 200
gen id = _n
gen bsweight = .
gen keep_flag = id <= 100

set seed 12345
cbsample if keep_flag == 1, weight(bsweight)

* All 200 obs should remain, but only first 100 should have non-zero weight
count if bsweight > 0 & id > 100
local n_wrong = r(N)
count if bsweight > 0 & id <= 100
local n_sampled = r(N)

if `n_wrong' == 0 & `n_sampled' > 0 {
    test_pass "if condition respected (weights only in valid range)"
}
else {
    test_fail "if condition respected" "n_wrong=`n_wrong', n_sampled=`n_sampled'"
}

* in condition test
clear
set obs 200
gen id = _n
gen bsweight = .

set seed 12345
cbsample in 1/100, weight(bsweight)

* Only first 100 should have non-zero weight
count if bsweight > 0 & id > 100
local n_wrong = r(N)
count if bsweight > 0 & id <= 100
local n_sampled = r(N)

if `n_wrong' == 0 & `n_sampled' > 0 {
    test_pass "in condition respected (weights only in valid range)"
}
else {
    test_fail "in condition respected" "n_wrong=`n_wrong', n_sampled=`n_sampled'"
}

* Combined if/in test
clear
set obs 300
gen id = _n
gen bsweight = .
gen group = mod(_n, 3)

set seed 12345
cbsample if group == 1 in 1/200, weight(bsweight)

* Only obs with id<=200 and group==1 should have non-zero weight
count if bsweight > 0 & id > 200
local n_wrong_in = r(N)
count if bsweight > 0 & group != 1
local n_wrong_if = r(N)
count if bsweight > 0 & id <= 200 & group == 1
local n_sampled = r(N)

if `n_wrong_in' == 0 & `n_wrong_if' == 0 & `n_sampled' > 0 {
    test_pass "combined if/in conditions"
}
else {
    test_fail "combined if/in conditions" "n_sampled=`n_sampled', wrong_in=`n_wrong_in', wrong_if=`n_wrong_if'"
}

/*******************************************************************************
 * SECTION 9: Multiple strata/cluster variables
 ******************************************************************************/
print_section "Multiple Strata/Cluster Variables"

* Multiple strata variables
clear
set obs 1000
gen strata1 = mod(_n, 5)
gen strata2 = mod(_n, 4)
gen id = _n
gen bsweight = .
set seed 12345

cbsample, strata(strata1 strata2) weight(bsweight)

* Check each stratum combination has weight
local all_ok = 1
forvalues s1 = 0/4 {
    forvalues s2 = 0/3 {
        summarize bsweight if strata1 == `s1' & strata2 == `s2', meanonly
        if r(sum) == 0 | r(sum) == . {
            local all_ok = 0
        }
    }
}
if `all_ok' {
    test_pass "multiple strata variables: all combinations sampled"
}
else {
    test_fail "multiple strata variables" "some combinations have zero weight"
}

* Multiple cluster variables (treated as compound cluster ID)
clear
set obs 1000
gen cluster1 = ceil(_n / 20)   /* 50 first-level clusters */
gen cluster2 = mod(_n, 2)       /* 2 second-level clusters within each */
gen id = _n
gen bsweight = .
set seed 12345

cbsample, cluster(cluster1 cluster2) weight(bsweight)

* Within each compound cluster, weights should be uniform
local all_ok = 1
forvalues c1 = 1/50 {
    forvalues c2 = 0/1 {
        summarize bsweight if cluster1 == `c1' & cluster2 == `c2'
        if r(sd) > 0.001 & r(sd) != . {
            local all_ok = 0
        }
    }
}
if `all_ok' {
    test_pass "multiple cluster variables: uniform weights within compound clusters"
}
else {
    test_fail "multiple cluster variables" "weights vary within some compound clusters"
}

/*******************************************************************************
 * SECTION 10: Combined strata and cluster
 ******************************************************************************/
print_section "Combined Strata and Cluster"

* Create data where clusters are NESTED within strata
* 4 strata, each with 25 clusters of 20 obs each = 2000 obs
clear
set obs 2000
gen stratum = floor((_n - 1) / 500)           /* 4 strata of 500 obs each */
gen cluster_id = floor((_n - 1) / 20)          /* 100 clusters of 20 obs each */
gen id = _n
gen bsweight = .
set seed 12345

cbsample, strata(stratum) cluster(cluster_id) weight(bsweight)

* Within each cluster, weights should be uniform
local cluster_ok = 1
levelsof cluster_id, local(clusters)
foreach cl of local clusters {
    summarize bsweight if cluster_id == `cl'
    if r(sd) > 0.001 & r(sd) != . {
        local cluster_ok = 0
    }
}

* Each stratum should have some weight
local strata_ok = 1
forvalues s = 0/3 {
    summarize bsweight if stratum == `s', meanonly
    if r(sum) == 0 | r(sum) == . {
        local strata_ok = 0
    }
}

if `cluster_ok' & `strata_ok' {
    test_pass "combined strata+cluster: both respected"
}
else {
    test_fail "combined strata+cluster" "cluster_ok=`cluster_ok', strata_ok=`strata_ok'"
}

/*******************************************************************************
 * SECTION 11: Without weight option (expand behavior)
 ******************************************************************************/
print_section "Without Weight Option (Expand)"

clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N

cbsample

* Without weight(), data should be expanded/contracted via expand
* Total obs should be approximately equal to N_before (may differ due to sampling)
local N_after = _N
if `N_after' > 0 & `N_after' <= 2 * `N_before' {
    test_pass "without weight(): data modified (N went from `N_before' to `N_after')"
}
else {
    test_fail "without weight(): data modified" "N=`N_after' from `N_before' (unexpected size)"
}

* With explicit n, expanded data should sum to approximately n
clear
set obs 100
gen id = _n
set seed 12345
cbsample 150

local N_after = _N
* Should have approximately 150 obs (some variance expected)
if `N_after' > 110 & `N_after' < 190 {
    test_pass "without weight() n=150: expanded to ~150 (got `N_after')"
}
else {
    test_fail "without weight() n=150: expanded to ~150" "got `N_after'"
}

/*******************************************************************************
 * SECTION 12: Strata independence verification
 ******************************************************************************/
print_section "Strata Independence Verification"

* Verify that strata are truly sampled independently
* Each stratum should have total weight approximately equal to its size
clear
set obs 1000
gen stratum = cond(_n <= 200, 0, cond(_n <= 600, 1, 2))
* Stratum 0: 200 obs, Stratum 1: 400 obs, Stratum 2: 400 obs
gen id = _n
gen bsweight = .
set seed 12345

cbsample, strata(stratum) weight(bsweight)

* Total weight per stratum should be close to stratum size
local all_ok = 1
summarize bsweight if stratum == 0, meanonly
if abs(r(sum) - 200) > 20 {
    local all_ok = 0
}
summarize bsweight if stratum == 1, meanonly
if abs(r(sum) - 400) > 30 {
    local all_ok = 0
}
summarize bsweight if stratum == 2, meanonly
if abs(r(sum) - 400) > 30 {
    local all_ok = 0
}

if `all_ok' {
    test_pass "strata independence: each stratum weighted ~proportionally"
}
else {
    test_fail "strata independence" "stratum weights not proportional to size"
}

/*******************************************************************************
 * SECTION 13: Cluster sampling verification
 ******************************************************************************/
print_section "Cluster Sampling Verification"

* Verify clusters are sampled as units (all or nothing within cluster)
clear
set obs 500
gen cluster_id = ceil(_n / 5)  /* 100 clusters of 5 each */
gen id = _n
gen bsweight = .
set seed 12345

cbsample, cluster(cluster_id) weight(bsweight)

* For each cluster, either all have same weight or all have zero
* (With bootstrap, a cluster may be sampled multiple times, so weight can vary
*  but within a cluster, all obs should have the same weight)
local all_ok = 1
forvalues cl = 1/100 {
    summarize bsweight if cluster_id == `cl'
    if r(sd) > 0.001 & r(sd) != . {
        local all_ok = 0
    }
}
if `all_ok' {
    test_pass "cluster sampling: all obs in cluster have same weight"
}
else {
    test_fail "cluster sampling" "weights vary within some clusters"
}

* Verify some clusters have zero weight (weren't selected)
count if bsweight == 0
local n_zero = r(N)
if `n_zero' > 0 {
    test_pass "cluster sampling: some obs have zero weight (clusters not selected)"
}
else {
    test_fail "cluster sampling" "all obs have non-zero weight (unlikely for n=100 clusters)"
}

/*******************************************************************************
 * SECTION 14: Verbose and threads options
 ******************************************************************************/
print_section "Verbose and Threads Options"

* Verbose option should not crash
clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
capture noisily cbsample, weight(bsweight) verbose
if _rc == 0 {
    test_pass "verbose option runs without error"
}
else {
    test_fail "verbose option runs without error" "rc=`=_rc'"
}

* Threads option should not crash
clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
capture cbsample, weight(bsweight) threads(2)
if _rc == 0 {
    test_pass "threads option runs without error"
}
else {
    test_fail "threads option runs without error" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cbsample returns the same error codes as Stata's bsample
 * when given invalid inputs or error conditions.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Weight variable doesn't exist
clear
set obs 100
gen id = _n
test_error_match, stata_cmd(bsample, weight(nonexistent_var)) ctools_cmd(cbsample, weight(nonexistent_var)) testname("nonexistent weight variable")

* Empty dataset
clear
set obs 0
gen id = .
gen bsweight = .
capture bsample, weight(bsweight)
local stata_rc = _rc
capture cbsample, weight(bsweight)
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' {
    test_pass "[error] empty dataset (rc=`stata_rc')"
}
else {
    test_fail "[error] empty dataset" "stata rc=`stata_rc', cbsample rc=`ctools_rc'"
}

* Strata variable doesn't exist
clear
set obs 100
gen id = _n
gen bsweight = .
test_error_match, stata_cmd(bsample, weight(bsweight) strata(nonexistent_var)) ctools_cmd(cbsample, weight(bsweight) strata(nonexistent_var)) testname("nonexistent strata variable")

* End of cbsample validation
noi print_summary "cbsample"
}
