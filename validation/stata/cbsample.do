* cbsample: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cbsample"

if "`__case'" == "cbsample" {
* validate_cbsample.do
* Validation tests for cbsample command
* Tests bootstrap sampling with replacement functionality

* Load setup (works from project root or validation dir)
do "validation/support/validate_setup.do"

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

* SECTION 2: Basic bootstrap sample
print_section "Basic Bootstrap Sample"

* Default n (same as original)
clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N
cbsample
local N_after = _N

* An unstratified bootstrap draws exactly the original sample size.
if `N_after' == `N_before' test_pass "basic bootstrap sample (N=`N_after' from `N_before')"
else test_fail "basic bootstrap sample" "N=`N_after' from `N_before' (unexpected size)"

* Explicit n (positional argument, like bsample)
clear
set obs 100
gen id = _n
set seed 12345
cbsample 50
if _N == 50 test_pass "bootstrap with n=50 (N=`=_N')"
else test_fail "bootstrap with n=50" "N=`=_N', expected 50"

* SECTION 2b: Self-consistency with set seed
print_section "Self-Consistency with set seed"

* Check every bootstrap weight, not the unchanged number of source rows.
clear
set obs 100
gen id = _n
gen w = .
set seed 12345
preserve
cbsample, weight(w)
tempfile weights_first
save `weights_first'
restore
set seed 12345
cbsample, weight(w)
capture cf _all using `weights_first'
if _rc == 0 test_pass "cbsample reproducible weights with same seed"
else test_fail "cbsample reproducible weights with same seed" "bootstrap weights differ"

* Compare the sampled observations, not the fixed requested sample size.
clear
set obs 100
gen id = _n
set seed 54321
preserve
cbsample 50
tempfile sample_first
save `sample_first'
restore
set seed 54321
cbsample 50
capture cf _all using `sample_first'
if _rc == 0 test_pass "cbsample 50 reproducible observations with same seed"
else test_fail "cbsample 50 reproducible observations with same seed" "sampled observations differ"

* SECTION 3: Weight variable (matches bsample syntax: requires existing var)
print_section "Weight Variable"

clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
cbsample, weight(bsweight)

* Should keep all obs and have weight variable populated
if _N == 100 test_pass "weight option keeps all obs"
else test_fail "weight option keeps all obs" "N=`=_N'"

* Bootstrap multiplicities are nonnegative integers for every input row.
count if missing(bsweight) | bsweight < 0 | bsweight != floor(bsweight)
if r(N) == 0 test_pass "weight variable contains nonnegative integer counts"
else test_fail "weight variable populated" "`r(N)' invalid weights"

* Randomness changes which rows are selected, not the number of draws.
summarize bsweight, meanonly
local total_weight = r(sum)
if `total_weight' == 100 test_pass "total weight exactly 100"
else test_fail "total weight exactly 100" "got `total_weight'"

* SECTION 4: Stratified bootstrap
print_section "Stratified Bootstrap"

clear
set obs 1000
gen stratum = mod(_n, 10)
gen id = _n
gen bsweight = .
set seed 12345
cbsample, strata(stratum) weight(bsweight)

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* Every stratum draws exactly its original number of observations.
local all_ok = 1
forvalues s = 0/9 {
    summarize bsweight if stratum == `s', meanonly
    if r(N) != 100 | r(sum) != 100 {
        local all_ok = 0
    }
}
if `all_ok' test_pass "stratified bootstrap has exactly 100 draws in every stratum"
else test_fail "stratified bootstrap" "incorrect per-stratum draw count"

* SECTION 5: Cluster bootstrap
print_section "Cluster Bootstrap"

clear
set obs 1000
gen cluster_id = ceil(_n / 10)  /* 100 clusters of 10 each */
gen id = _n
gen bsweight = .
set seed 12345
cbsample, cluster(cluster_id) weight(bsweight)

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* Within each cluster, all obs should have same weight
local all_ok = 1
levelsof cluster_id, local(clusters)
foreach cl of local clusters {
    summarize bsweight if cluster_id == `cl'
    if r(N) != 10 | r(min) != r(max) {
        local all_ok = 0
    }
}
if `all_ok' test_pass "cluster bootstrap: uniform weights within clusters"
else test_fail "cluster bootstrap" "weights vary within some clusters"

* SECTION 6: Reproducibility with set seed
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

if `n_diff' == 0 test_pass "same seed yields identical weights"
else test_fail "same seed yields identical weights" "`n_diff' weights differ"

* SECTION 7: Edge cases
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
assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
if `total' == 10 test_pass "small dataset bootstrap (total weight exactly 10)"
else test_fail "small dataset bootstrap" "total weight=`total'"

* Large n (positional argument, like bsample)
clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
test_error_match, stata_cmd(bsample 200, weight(bsweight)) ctools_cmd(cbsample 200, weight(bsweight)) testname("oversized bootstrap draw rejected")

* SECTION 8: If/In Conditions
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

if `n_wrong' == 0 & `n_sampled' > 0 test_pass "if condition respected (weights only in valid range)"
else test_fail "if condition respected" "n_wrong=`n_wrong', n_sampled=`n_sampled'"

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

if `n_wrong' == 0 & `n_sampled' > 0 test_pass "in condition respected (weights only in valid range)"
else test_fail "in condition respected" "n_wrong=`n_wrong', n_sampled=`n_sampled'"

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

if `n_wrong_in' == 0 & `n_wrong_if' == 0 & `n_sampled' > 0 test_pass "combined if/in conditions"
else test_fail "combined if/in conditions" "n_sampled=`n_sampled', wrong_in=`n_wrong_in', wrong_if=`n_wrong_if'"

* SECTION 9: Multiple strata/cluster variables
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

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* Each of the 20 combinations has exactly 50 observations and draws.
local all_ok = 1
forvalues s1 = 0/4 {
    forvalues s2 = 0/3 {
        summarize bsweight if strata1 == `s1' & strata2 == `s2', meanonly
        if r(N) != 50 | r(sum) != 50 {
            local all_ok = 0
        }
    }
}
if `all_ok' test_pass "multiple strata variables: exactly 50 draws per combination"
else test_fail "multiple strata variables" "incorrect per-combination draw count"

* Multiple cluster variables (treated as compound cluster ID)
clear
set obs 1000
gen cluster1 = ceil(_n / 20)   /* 50 first-level clusters */
gen cluster2 = mod(_n, 2)       /* 2 second-level clusters within each */
gen id = _n
gen bsweight = .
set seed 12345

cbsample, cluster(cluster1 cluster2) weight(bsweight)

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* Within each compound cluster, weights should be uniform
local all_ok = 1
forvalues c1 = 1/50 {
    forvalues c2 = 0/1 {
        summarize bsweight if cluster1 == `c1' & cluster2 == `c2'
        if r(N) != 10 | r(min) != r(max) {
            local all_ok = 0
        }
    }
}
if `all_ok' test_pass "multiple cluster variables: uniform weights within compound clusters"
else test_fail "multiple cluster variables" "weights vary within some compound clusters"

* SECTION 10: Combined strata and cluster
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

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* Within each cluster, weights should be uniform
local cluster_ok = 1
levelsof cluster_id, local(clusters)
foreach cl of local clusters {
    summarize bsweight if cluster_id == `cl'
    if r(N) != 20 | r(min) != r(max) {
        local cluster_ok = 0
    }
}

* Equal-size clusters imply exactly 25*20 weighted rows in each stratum.
local strata_ok = 1
forvalues s = 0/3 {
    summarize bsweight if stratum == `s', meanonly
    if r(N) != 500 | r(sum) != 500 {
        local strata_ok = 0
    }
}

if `cluster_ok' & `strata_ok' test_pass "combined strata+cluster: both respected"
else test_fail "combined strata+cluster" "cluster_ok=`cluster_ok', strata_ok=`strata_ok'"

* SECTION 11: Without weight option (expand behavior)
print_section "Without Weight Option (Expand)"

clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N

cbsample

* Without clusters, the draw count equals the original number of rows.
local N_after = _N
if `N_after' == `N_before' test_pass "without weight(): exact bootstrap draw count"
else test_fail "without weight(): data modified" "N=`N_after' from `N_before' (unexpected size)"

* With explicit n, expanded data must contain exactly n observations.
clear
set obs 100
gen id = _n
set seed 12345
cbsample 50
assert _N == 50
test_pass "explicit count expands to exactly 50 observations"

* SECTION 12: Strata independence verification
print_section "Strata Independence Verification"

* Verify that strata are truly sampled independently
* Each stratum must have total weight exactly equal to its size.
clear
set obs 1000
gen stratum = cond(_n <= 200, 0, cond(_n <= 600, 1, 2))
* Stratum 0: 200 obs, Stratum 1: 400 obs, Stratum 2: 400 obs
gen id = _n
gen bsweight = .
set seed 12345

cbsample, strata(stratum) weight(bsweight)

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* These are fixed draw counts, not estimated population proportions.
local all_ok = 1
summarize bsweight if stratum == 0, meanonly
if r(N) != 200 | r(sum) != 200 {
    local all_ok = 0
}
summarize bsweight if stratum == 1, meanonly
if r(N) != 400 | r(sum) != 400 {
    local all_ok = 0
}
summarize bsweight if stratum == 2, meanonly
if r(N) != 400 | r(sum) != 400 {
    local all_ok = 0
}

if `all_ok' test_pass "strata independence: exact draw count in each stratum"
else test_fail "strata independence" "incorrect per-stratum draw count"

* SECTION 13: Cluster sampling verification
print_section "Cluster Sampling Verification"

* Verify clusters are sampled as units (all or nothing within cluster)
clear
set obs 500
gen cluster_id = ceil(_n / 5)  /* 100 clusters of 5 each */
gen id = _n
gen bsweight = .
set seed 12345

cbsample, cluster(cluster_id) weight(bsweight)

assert !missing(bsweight) & bsweight >= 0 & bsweight == floor(bsweight)
* For each cluster, either all have same weight or all have zero
* (With bootstrap, a cluster may be sampled multiple times, so weight can vary
*  but within a cluster, all obs should have the same weight)
local all_ok = 1
forvalues cl = 1/100 {
    summarize bsweight if cluster_id == `cl'
    if r(N) != 5 | r(min) != r(max) {
        local all_ok = 0
    }
}
if `all_ok' test_pass "cluster sampling: all obs in cluster have same weight"
else test_fail "cluster sampling" "weights vary within some clusters"

* Verify some clusters have zero weight (weren't selected)
count if bsweight == 0
local n_zero = r(N)
if `n_zero' > 0 test_pass "cluster sampling: some obs have zero weight (clusters not selected)"
else test_fail "cluster sampling" "all obs have non-zero weight (unlikely for n=100 clusters)"

* SECTION 14: Verbose and threads options
print_section "Verbose and Threads Options"

* Verbose option should not crash
clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
capture noisily cbsample, weight(bsweight) verbose
if _rc == 0 test_pass "verbose option runs without error"
else test_fail "verbose option runs without error" "rc=`=_rc'"

* Threads option should not crash
clear
set obs 100
gen id = _n
gen bsweight = .
set seed 12345
capture cbsample, weight(bsweight) threads(2)
if _rc == 0 test_pass "threads option runs without error"
else test_fail "threads option runs without error" "rc=`=_rc'"

* SECTION: Intentional Error Tests
* These tests verify that cbsample returns the same error codes as Stata's bsample
* when given invalid inputs or error conditions.
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
if `stata_rc' == `ctools_rc' test_pass "[error] empty dataset (rc=`stata_rc')"
else test_fail "[error] empty dataset" "stata rc=`stata_rc', cbsample rc=`ctools_rc'"

* Strata variable doesn't exist
clear
set obs 100
gen id = _n
gen bsweight = .
test_error_match, stata_cmd(bsample, weight(bsweight) strata(nonexistent_var)) ctools_cmd(cbsample, weight(bsweight) strata(nonexistent_var)) testname("nonexistent strata variable")

* Empty selection semantics differ from empty data and empty strata.
foreach mode in empty ifzero nstrata sstrata {
    foreach weighted in 0 1 {
        clear
        if "`mode'" != "empty" set obs 10
        gen id = _n
        gen w = 99
        gen sn = .
        gen str1 ss = ""
        local select
        local opts
        if "`mode'" == "ifzero" local select if 0
        if "`mode'" == "nstrata" local opts strata(sn)
        if "`mode'" == "sstrata" local opts strata(ss)
        if `weighted' local opts `opts' weight(w)
        preserve
        capture bsample `select', `opts'
        local reference_rc = _rc
        * Native stratified bsample sorts ties; compare by original row ID.
        sort id
        tempfile expected
        save `expected'
        restore
        * Prior successful-call output must not masquerade as this draw.
        scalar _cbsample_total_weight = 123
        local rng_before `c(rngstate)'
        capture cbsample `select', `opts'
        local rc = _rc
        assert `rc' == `reference_rc'
        cf _all using `expected'
        assert "`c(rngstate)'" == "`rng_before'"
        if `rc' {
            capture confirm scalar _cbsample_total_weight
            assert _rc == 111
        }
        else assert scalar(_cbsample_total_weight) == 0
        test_pass "empty sample parity: `mode', weight=`weighted'"
    }
}
clear
set obs 3
gen id = _n
gen str1 w = "a"
foreach count in 0 invalid {
    capture cbsample `count' if 0
    assert _rc == 198
    assert _N == 3 & id == _n & w == "a"
}
capture cbsample if 0, weight(w)
assert _rc == 7
assert _N == 3 & id == _n & w == "a"
test_pass "empty selection still validates count and weight type"

* Even on zero rows, native's explicit selection/count paths can error.
foreach command in "if 0" "if 1" "1" ", strata(sn)" ", cluster(id)" {
    clear
    gen id = .
    gen sn = .
    preserve
    capture bsample `command'
    local reference_rc = _rc
    restore
    capture cbsample `command'
    local rc = _rc
    assert `rc' == `reference_rc'
    assert _N == 0 & c(k) == 2
    test_pass "zero-row explicit syntax parity: `command'"
}

* DATA-06: the resample size is an expression; DAT-06: singleton clusters
clear
set obs 60
gen id = _n
gen s = mod(_n, 3)
gen cl = floor(_n / 4)
replace cl = 999 if s == 2
gen one = 1
foreach call in "_N/2" "round(0.1*_N)" "5*4" "0" "2.5" "70" "_N/2, strata(s)" ///
    ", cluster(cl) strata(s)" ", cluster(one)" "3, cluster(one)" {
    preserve
    capture bsample `call'
    local rc_native = _rc
    local n_native = _N
    restore, preserve
    capture cbsample `call'
    local rc_ctools = _rc
    if `rc_native' == `rc_ctools' & _N == `n_native' test_pass "DATA-06/DAT-06 cbsample `call' (r(`rc_native'), N=`n_native')"
    else test_fail "DATA-06/DAT-06 cbsample `call'" "bsample r(`rc_native') N=`n_native', cbsample r(`rc_ctools') N=`=_N'"
    restore
}

* P2-2: idcluster(), bsample's new cluster id. The draws differ from
* bsample's (another RNG), so every output is checked against bsample's rules:
* each id is one whole copy of one cluster in observation order; ids run
* 1, 2, ... by stratum (on across strata under version 14+, within strata
* under version 13); rows sorted by (strata, id, observation); copies numbered
* in cluster order for full-size draws; the same draw as weight() and the same
* RNG use; and, against bsample on the same data, the variable's type, format,
* label, position and the sort order.
capture program drop _cbs_idcheck
program define _cbs_idcheck, rclass
    args strata cluster id full uniq clsize
    local bad
    tempvar o n
    gen long `o' = _n
    local newstratum = cond("`strata'" != "", "(`strata' != `strata'[_n-1]) |", "")
    capture assert `newstratum' (`id' > `id'[_n-1]) | (`id' == `id'[_n-1] & obs > obs[_n-1]) if _n > 1
    if _rc local bad `bad' order
    bysort `strata' `id' (`o'): gen `n' = _N
    capture by `strata' `id': assert `cluster' == `cluster'[1]
    if _rc local bad `bad' onecluster
    capture assert `n' == `clsize'
    if _rc local bad `bad' wholecluster
    sort `o'
    if `uniq' | "`strata'" == "" capture assert `id' == 1 if _n == 1
    else capture assert `id' == 1 if _n == 1 | `strata' != `strata'[_n-1]
    if _rc local bad `bad' firstid
    local same = cond("`strata'" != "", "& `strata' == `strata'[_n-1]", "")
    if `uniq' | "`strata'" == "" capture assert inlist(`id' - `id'[_n-1], 0, 1) if _n > 1
    else capture assert inlist(`id' - `id'[_n-1], 0, 1) if _n > 1 `same'
    if _rc local bad `bad' ids
    if `full' {
        capture assert `cluster' >= `cluster'[_n-1] if _n > 1 `same'
        if _rc local bad `bad' clusterorder
    }
    return local bad `bad'
end
foreach cfg in "" ", strata(s)" ", strata(ss)" "5, strata(s)" "4" "if keep, strata(s)" ///
    "sz, strata(s)" "v13, strata(s)" "v13 sz, strata(s)" {
    local fails
    forvalues seed = 1/20 {
        clear
        set seed `seed'
        set obs 120
        gen long obs = _n
        gen s = mod(_n, 3) + 1
        gen str3 ss = cond(s == 1, "b", cond(s == 2, "a", "c c"))
        gen cl = ceil(runiform() * 20)
        replace cl = . if mod(_n, 17) == 0
        replace cl = .a if mod(_n, 19) == 0
        gen keep = runiform() < .8
        gen sz = cond(s == 1, 3, cond(s == 2, 7, 5))
        local call `cfg'
        local ver
        if strpos("`call'", "v13") {
            local ver "version 13:"
            local call = subinstr("`call'", "v13", "", 1)
        }
        local st = cond(strpos("`call'", "ss)"), "ss", cond(strpos("`call'", "s)"), "s", ""))
        if strpos("`call'", "if keep") bysort `st' cl: egen clsize = total(keep)
        else bysort `st' cl: gen clsize = _N
        sort obs
        local full = !regexm("`call'", "^ *[0-9]") & !strpos("`call'", "sz")
        local uniq = "`ver'" == ""
        local call = cond(strpos("`call'", ","), "`call' cluster(cl)", "`call', cluster(cl)")
        preserve
        `ver' bsample `call' idcluster(nid)
        local meta_b "`: type nid' `: format nid' `: variable label nid' [`: sortedby']"
        unab vars_b : *
        _cbs_idcheck "`st'" cl nid `full' `uniq' clsize
        if "`r(bad)'" != "" local fails `fails' bsample`seed':`r(bad)'
        restore, preserve
        set seed `seed'
        gen double w = .
        `ver' cbsample `call' weight(w)
        local rng_w `c(rngstate)'
        tempfile draws
        keep obs w
        save `draws'
        restore, preserve
        set seed `seed'
        `ver' cbsample `call' idcluster(nid)
        local rng_id `c(rngstate)'
        local meta_c "`: type nid' `: format nid' `: variable label nid' [`: sortedby']"
        unab vars_c : *
        _cbs_idcheck "`st'" cl nid `full' `uniq' clsize
        local bad `r(bad)'
        if `"`meta_b' `vars_b'"' != `"`meta_c' `vars_c'"' local bad `bad' meta
        if "`rng_w'" != "`rng_id'" local bad `bad' rng
        bysort obs: gen copies = _N
        by obs: keep if _n == 1
        merge 1:1 obs using `draws', nogenerate
        capture assert (w == copies) | (w == 0 & missing(copies))
        if _rc local bad `bad' draw
        if "`bad'" != "" local fails `fails' cbsample`seed':`bad'
        restore
    }
    if "`fails'" == "" test_pass "P2-2 idcluster rules vs bsample: `cfg' (20 seeds)"
    else test_fail "P2-2 idcluster rules vs bsample: `cfg'" "`fails'"
}
* Full-size draws number the copies deterministically: the exact dataset
* bsample's NewClustID builds from the same draw (taken from weight())
local fails
forvalues seed = 1/20 {
    foreach ver in 13 19 {
        clear
        set seed `seed'
        set obs 90
        gen long obs = _n
        gen s = mod(_n, 3) + 1
        gen cl = ceil(runiform() * 12)
        preserve
        gen double w = .
        set seed `seed'
        cbsample, cluster(cl) strata(s) weight(w)
        keep if w > 0
        sort s cl obs
        expand w
        bysort obs: gen copy = _n
        sort s cl copy obs
        by s cl copy: gen byte first = _n == 1
        if `ver' > 13 gen nid = sum(first)
        else by s: gen nid = sum(first)
        sort s nid obs
        drop w copy first
        tempfile expected
        save `expected'
        restore, preserve
        set seed `seed'
        version `ver': cbsample, cluster(cl) strata(s) idcluster(nid)
        capture cf _all using `expected'
        if _rc local fails `fails' `seed'/v`ver'
        restore
    }
}
if "`fails'" == "" test_pass "P2-2 idcluster ids equal NewClustID on the same draw (20 seeds, v13/v19)"
else test_fail "P2-2 idcluster exact ids" "seeds `fails'"
* Storage type (bsample's byte id promoted by replace), errors, empty sample
local bad
foreach n in 100 101 32740 32741 {
    clear
    set seed 4
    set obs `n'
    gen id = _n
    cbsample, cluster(id) idcluster(nid)
    local want = cond(`n' <= 100, "byte", cond(`n' <= 32740, "int", "long"))
    if "`: type nid'" != "`want'" | "`: format nid'" != "%8.0g" local bad `bad' `n':`: type nid'
}
if "`bad'" == "" test_pass "P2-2 idcluster storage type and format"
else test_fail "P2-2 idcluster storage type" "`bad'"
clear
set obs 10
gen id = _n
gen w = .
gen s = 1
foreach o in "idcluster(x)" "strata(id) idcluster(x)" "cluster(id) idcluster(x) weight(w)" ///
    "cluster(id) idcluster(1x)" "cluster(id) idcluster(w)" "if 0, cluster(id) strata(s) idcluster(nid)" ///
    "if 0, cluster(id) idcluster(nid)" "in 1, cluster(id) idcluster(nid)" {
    local opt = cond(strpos("`o'", ","), "`o'", ", `o'")
    preserve
    capture bsample `opt'
    local rc_b = _rc
    unab vars_b : *
    local meta_b "`c(N)' [`: sortedby']"
    capture local meta_b "`meta_b' `: type nid' `: variable label nid'"
    restore, preserve
    capture cbsample `opt'
    local rc_c = _rc
    unab vars_c : *
    local meta_c "`c(N)' [`: sortedby']"
    capture local meta_c "`meta_c' `: type nid' `: variable label nid'"
    restore
    if `rc_b' == `rc_c' & "`vars_b' `meta_b'" == "`vars_c' `meta_c'" test_pass "P2-2 idcluster parity: `o' (r(`rc_b'))"
    else test_fail "P2-2 idcluster parity: `o'" "bsample r(`rc_b') `vars_b' `meta_b'; cbsample r(`rc_c') `vars_c' `meta_c'"
}

* P3-13: resample sizes that differ across strata, as bsample draws them:
* the same rc, and the same number of draws in every stratum
clear
set seed 2026
set obs 90
gen long obs = _n
gen cl = ceil(_n / 3)
gen s = mod(cl, 3) + 1
gen double w = .
gen sz = cond(s == 1, 3, cond(s == 2, 7, 5))
gen szfrac = sz + (s == 2) / 2
gen szzero = sz * (s != 3)
gen szbig = sz + 40 * (s == 1)
gen s2 = _n <= 30
foreach call in "sz, strata(s)" "sz, strata(s) weight(w)" "sz, strata(s) cluster(cl)" ///
    "sz, strata(s) cluster(cl) weight(w)" "_N/2, strata(s2)" "_N/3, strata(s2) cluster(cl)" ///
    "szfrac, strata(s)" "szzero, strata(s)" "szbig, strata(s)" "szbig, strata(s) cluster(cl)" ///
    "sz if obs > 9, strata(s)" {
    local draws_bsample
    local draws_cbsample
    foreach cmd in bsample cbsample {
        preserve
        set seed 11
        capture `cmd' `call'
        local rc_`cmd' = _rc
        if !_rc {
            local st = cond(strpos("`call'", "s2)"), "s2", "s")
            * clusters drawn per stratum: weight() counts each copy once per
            * cluster, expanded data once per row
            if strpos("`call'", "weight(w)") {
                if strpos("`call'", "cluster(cl)") bysort `st' cl: keep if _n == 1
                gen double drawn = w
            }
            else gen double drawn = cond(strpos("`call'", "cluster(cl)"), 1 / 3, 1)
            collapse (sum) drawn, by(`st')
            mata: st_local("draws_`cmd'", invtokens(strofreal(round(st_data(., "drawn"))')))
        }
        restore
    }
    if `rc_bsample' == `rc_cbsample' & "`draws_bsample'" == "`draws_cbsample'" ///
        test_pass "P3-13 per-stratum sizes: `call' (r(`rc_bsample')) draws `draws_bsample'"
    else test_fail "P3-13 per-stratum sizes: `call'" "bsample r(`rc_bsample') [`draws_bsample'], cbsample r(`rc_cbsample') [`draws_cbsample']"
}

* _ct_output runs a command and returns its output (r(text), lines joined
* by |) and return code (r(rc)). _ct_same_output runs a native and a ctools
* command on copies of the data in memory and passes when the output, return
* code, variable names and types, and (unless "nodata") values agree.
capture program drop _ct_output
program define _ct_output, rclass
    tempfile f
    tempname fh
    capture log close _ctout
    quietly log using `"`f'"', text replace name(_ctout)
    capture noisily `0'
    local rc = _rc
    quietly log close _ctout
    local text
    file open `fh' using `"`f'"', read text
    file read `fh' line
    while !r(eof) {
        if `"`macval(line)'"' != "" & substr(`"`macval(line)'"', 1, 5) != "r; t=" {
            local text `"`macval(text)'`macval(line)'|"'
        }
        file read `fh' line
    }
    file close `fh'
    return local text `"`macval(text)'"'
    return scalar rc = `rc'
end
capture program drop _ct_same_output
program define _ct_same_output
    args name native ctools nodata
    tempfile want
    preserve
    _ct_output `native'
    local want_text `"`r(text)'"'
    local want_rc = r(rc)
    local want_vars
    local want_types
    forvalues j = 1/`c(k)' {
        mata: st_local("v", st_varname(`j'))
        local want_vars `want_vars' `v'
        local want_types `want_types' `: type `v''
    }
    capture quietly save `"`want'"', emptyok
    restore, preserve
    _ct_output `ctools'
    local got_text `"`r(text)'"'
    local got_rc = r(rc)
    local got_vars
    local got_types
    forvalues j = 1/`c(k)' {
        mata: st_local("v", st_varname(`j'))
        local got_vars `got_vars' `v'
        local got_types `got_types' `: type `v''
    }
    local same_data = 1
    if "`nodata'" == "" & "`got_vars'" != "" {
        capture cf _all using `"`want'"'
        local same_data = !_rc
    }
    restore
    if "`nodata'" != "" {
        local want_vars
        local got_vars
        local want_types
        local got_types
    }
    if `"`macval(want_text)'"' == `"`macval(got_text)'"' & `want_rc' == `got_rc' & ///
        "`want_vars'" == "`got_vars'" & "`want_types'" == "`got_types'" & `same_data' {
        test_pass "`name'"
    }
    else {
        noisily di as error `"`name': ctools [`macval(got_text)'] rc `got_rc' vars `got_vars' (`got_types')"'
        noisily di as error `"`name': native [`macval(want_text)'] rc `want_rc' vars `want_vars' (`want_types') same data `same_data'"'
        test_fail "`name'" "output, return code, variables or data differ from the native command"
    }
end
* BSAMP-02: under version 7 or earlier, bsample's old algorithm (bsample_7)
* draws the same resample, with its messages and return codes
clear
set obs 30
gen long id = _n
gen byte c = mod(_n, 6)
gen byte s = mod(_n, 2)
foreach pair in `""' `"20"' `", cluster(c)"' `"4, cluster(c)"' `"4, cluster(c) idcluster(nc)"' ///
    `", cluster(c) idcluster(id)"' `"40"' `"0"' `"7, cluster(c)"' `", idcluster(nc)"' `"2.5"' ///
    `", strata(s)"' `"10 if s"' `", weight(id)"' `"20, threads(2) verbose"' {
    tempfile want
    local npair : subinstr local pair "threads(2) verbose" ""
    preserve
    set seed 23
    _ct_output version 7: bsample `npair'
    local want_text `"`r(text)'"'
    local want_rc = r(rc)
    quietly save `want'
    restore, preserve
    set seed 23
    _ct_output version 7: cbsample `pair'
    local got_text `"`r(text)'"'
    local got_rc = r(rc)
    capture cf _all using `want'
    local same = !_rc
    restore
    if `same' & `got_rc' == `want_rc' & `"`got_text'"' == `"`want_text'"' {
        test_pass "BSAMP-02 version 7 resample: `pair'"
    }
    else test_fail "BSAMP-02 version 7 resample: `pair'" "rc `got_rc' vs `want_rc', same data `same'"
}
* Version 8 and later use the ctools draw; version 7 does not leak into it
set seed 23
version 8: cbsample 10
if _N == 10 test_pass "BSAMP-02 version 8 uses cbsample"
else test_fail "BSAMP-02 version 8 uses cbsample" "N=`=_N'"
* End of cbsample validation
noi print_summary "cbsample"
}

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "cbsample"
    global CTOOLS_CASE_COMPLETE "cbsample"
    exit
}

di as error "Unknown cbsample case: `__case'"
exit 198
