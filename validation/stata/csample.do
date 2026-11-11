* csample: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "csample"

if "`__case'" == "csample" {
* validate_csample.do
* Validation tests for csample command
* Tests random sampling without replacement functionality

* Load setup (works from project root or validation dir)
do "validation/support/validate_setup.do"

quietly {

noi di as text "Running validation tests for csample..."

* Plugin functionality check
clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N
capture csample 50
if _rc != 0 {
    test_fail "csample plugin load" "plugin returned error `=_rc'"
    exit 1
}
if _N == 50 {
    test_pass "csample plugin loads and runs (50% of 100 = `=_N' obs)"
}
else {
    test_fail "csample plugin load" "expected 50 observations, got `=_N'"
    exit 1
}

* SECTION 2: Basic percentage sampling
print_section "Basic Percentage Sampling"

* 50% sample
clear
set obs 1000
gen id = _n
set seed 12345
csample 50
local n_after = _N
if `n_after' == 500 test_pass "50% sample yields exactly 500 obs"
else test_fail "50% sample yields exactly 500 obs" "got `n_after'"

* 25% sample
clear
set obs 1000
gen id = _n
set seed 12345
csample 25
local n_after = _N
if `n_after' == 250 test_pass "25% sample yields exactly 250 obs"
else test_fail "25% sample yields exactly 250 obs" "got `n_after'"

* SECTION 2b: Self-consistency with set seed
print_section "Self-Consistency with set seed"

* csample with same seed should produce identical results
clear
set obs 1000
gen id = _n

set seed 54321
preserve
csample 50
local N1 = _N
tempfile csample_result
save `csample_result'
restore

set seed 54321
preserve
csample 50
local N2 = _N

if `N1' == `N2' {
    merge 1:1 id using `csample_result'
    quietly count if _merge != 3
    if r(N) == 0 test_pass "csample reproducible with same seed (50%)"
    else test_fail "csample reproducible with same seed" "`=r(N)' observations differ"
}
else {
    test_fail "csample reproducible with same seed" "N differ: run1=`N1' run2=`N2'"
}
restore

* SECTION 3: Fixed count sampling
print_section "Fixed Count Sampling"

clear
set obs 1000
gen id = _n
set seed 12345
csample, count(100)
if _N == 100 test_pass "count(100) yields exactly 100 obs"
else test_fail "count(100) yields exactly 100 obs" "got `=_N'"

clear
set obs 500
gen id = _n
set seed 12345
csample, count(50)
if _N == 50 test_pass "count(50) yields exactly 50 obs"
else test_fail "count(50) yields exactly 50 obs" "got `=_N'"

* Deterministic comparison: count(100) csample vs itself
clear
set obs 1000
gen id = _n

set seed 99999
preserve
csample, count(100)
sort id
tempfile csample_result
save `csample_result'
restore

set seed 99999
preserve
csample, count(100)
sort id

merge 1:1 id using `csample_result'
quietly count if _merge != 3
if r(N) == 0 test_pass "count(100) reproducible with same seed"
else test_fail "count(100) reproducible with same seed" "`=r(N)' observations differ"
restore

* SECTION 4: By-group sampling
print_section "By-Group Sampling"

* By-group percentage
clear
set obs 1000
gen group = mod(_n, 10)
gen id = _n
set seed 12345
csample 50, by(group)

* The selected identities are random; the sample size in each group is fixed.
local all_ok = 1
forvalues g = 0/9 {
    count if group == `g'
    local n_g = r(N)
    if `n_g' != 50 {
        local all_ok = 0
    }
}
if `all_ok' test_pass "50% by-group sample (exactly 50 per group)"
else test_fail "50% by-group sample" "some groups do not contain exactly 50 observations"

* By-group fixed count
clear
set obs 1000
gen group = mod(_n, 10)
gen id = _n
set seed 12345
csample, count(10) by(group)

local all_ok = 1
forvalues g = 0/9 {
    count if group == `g'
    local n_g = r(N)
    if `n_g' != 10 {
        local all_ok = 0
    }
}
if `all_ok' & _N == 100 test_pass "count(10) by-group yields 10 per group, 100 total"
else test_fail "count(10) by-group" "expected 100 total, got `=_N'"

* SECTION 5: Reproducibility with set seed
print_section "Reproducibility"

clear
set obs 1000
gen id = _n
set seed 99999
csample 30
tempfile sample1
save `sample1'
local N1 = _N

clear
set obs 1000
gen id = _n
set seed 99999
csample 30
tempfile sample2
save `sample2'
local N2 = _N

* Check same count
if `N1' == `N2' test_pass "same seed yields same count"
else test_fail "same seed yields same count" "N1=`N1' N2=`N2'"

* Check same observations
use `sample1', clear
merge 1:1 id using `sample2'
count if _merge != 3
local n_diff = r(N)
if `n_diff' == 0 test_pass "same seed yields identical samples"
else test_fail "same seed yields identical samples" "`n_diff' obs differ"

* SECTION 6: Edge cases
print_section "Edge Cases"

* 100% sample
clear
set obs 100
gen id = _n
set seed 12345
csample 100
if _N == 100 test_pass "100% sample keeps all obs"
else test_fail "100% sample keeps all obs" "got `=_N'"

* 0% sample
clear
set obs 100
gen id = _n
set seed 12345
csample 0
if _N == 0 test_pass "0% sample drops all obs"
else test_fail "0% sample drops all obs" "got `=_N'"

* Small dataset
clear
set obs 10
gen id = _n
set seed 12345
csample 50
if _N == 5 test_pass "small dataset sample (N=10, 50%)"
else test_fail "small dataset sample" "got `=_N', expected 5"

* SECTION 7: if/in conditions
print_section "If/In Conditions"

* if condition test
clear
set obs 200
gen id = _n
gen x = runiform()
gen keep_flag = id <= 100

preserve
sample 50 if keep_flag == 1
local expected_N = _N
count if id > 100
local expected_outside = r(N)
restore
set seed 12345
csample 50 if keep_flag == 1

local n_kept = _N
count if id > 100
local n_wrong = r(N)

if `n_wrong' == `expected_outside' & `n_kept' == `expected_N' test_pass "if condition respected"
else test_fail "if condition respected" "n_wrong=`n_wrong', n_kept=`n_kept'"

* in condition test
clear
set obs 200
gen id = _n
preserve
sample 50 in 1/100
local expected_N = _N
count if id > 100
local expected_outside = r(N)
restore
set seed 12345
csample 50 in 1/100

local n_kept = _N
count if id > 100
local n_wrong = r(N)

if `n_wrong' == `expected_outside' & `n_kept' == `expected_N' test_pass "in condition respected"
else test_fail "in condition respected" "n_wrong=`n_wrong', n_kept=`n_kept'"

* Combined if/in test
clear
set obs 300
gen id = _n
gen group = mod(_n, 3)

preserve
sample 50 if group == 1 in 1/200
local expected_N = _N
count if id > 200
local expected_outside_in = r(N)
count if group != 1
local expected_outside_if = r(N)
restore
set seed 12345
csample 50 if group == 1 in 1/200

local n_kept = _N
count if id > 200
local n_wrong_in = r(N)
count if group != 1
local n_wrong_if = r(N)

if `n_wrong_in' == `expected_outside_in' & `n_wrong_if' == `expected_outside_if' & `n_kept' == `expected_N' test_pass "combined if/in conditions"
else test_fail "combined if/in conditions" "n_kept=`n_kept', wrong_in=`n_wrong_in', wrong_if=`n_wrong_if'"

* SECTION 8: Multiple by-variables
print_section "Multiple By-Variables"

clear
set obs 1000
gen group1 = mod(_n, 5)
gen group2 = mod(_n, 4)
gen id = _n
set seed 12345

* There are 20 unique combinations (5*4), each with ~50 obs
csample, count(5) by(group1 group2)

* Check total count is 20 groups * 5 = 100
if _N == 100 test_pass "multiple by-vars: correct total count"
else test_fail "multiple by-vars: correct total count" "expected 100, got `=_N'"

* Verify each combination has exactly 5 obs
local all_ok = 1
forvalues g1 = 0/4 {
    forvalues g2 = 0/3 {
        count if group1 == `g1' & group2 == `g2'
        if r(N) != 5 {
            local all_ok = 0
        }
    }
}
if `all_ok' test_pass "multiple by-vars: 5 per cell"
else test_fail "multiple by-vars: 5 per cell" "some cells don't have 5"

* SECTION 9: Verbose and threads options
print_section "Verbose and Threads Options"

* Verbose option should not crash
clear
set obs 100
gen id = _n
set seed 12345
capture noisily csample 50, verbose
if _rc == 0 test_pass "verbose option runs without error"
else test_fail "verbose option runs without error" "rc=`=_rc'"

* Threads option should not crash
clear
set obs 100
gen id = _n
set seed 12345
capture csample 50, threads(2)
if _rc == 0 test_pass "threads option runs without error"
else test_fail "threads option runs without error" "rc=`=_rc'"

* SECTION 10: Error handling
print_section "Error Handling"

* count > N should still work (samples all available)
clear
set obs 50
gen id = _n
set seed 12345
capture csample, count(100) by(id)
* This should succeed but only get 50 (1 per group, 50 groups)
if _rc == 0 & _N == 50 test_pass "count > group size samples all available"
else test_fail "count > group size samples all available" "rc=`=_rc', N=`=_N'"

* Invalid percent should error
clear
set obs 100
gen id = _n
capture csample 150
if _rc != 0 test_pass "invalid percent (>100) errors correctly"
else test_fail "invalid percent (>100) errors correctly" "should have errored"

* Negative percent should error
clear
set obs 100
gen id = _n
capture csample -10
if _rc != 0 test_pass "negative percent errors correctly"
else test_fail "negative percent errors correctly" "should have errored"

* SECTION: Additional Intentional Error Tests
* These tests verify that csample returns the same error codes as Stata's sample
* when given invalid inputs or error conditions.
print_section "Additional Intentional Error Tests"

* Empty dataset
clear
set obs 0
gen id = .
capture sample 50
local stata_rc = _rc
capture csample 50
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' test_pass "[error] empty dataset (rc=`stata_rc')"
else test_fail "[error] empty dataset" "stata rc=`stata_rc', csample rc=`ctools_rc'"

* Zero percent
clear
set obs 100
gen id = _n
capture sample 0
local stata_rc = _rc
local stata_n = _N
clear
set obs 100
gen id = _n
capture csample 0
local ctools_rc = _rc
local ctools_n = _N
if `stata_rc' == `ctools_rc' test_pass "[error] zero percent (rc=`stata_rc')"
else test_fail "[error] zero percent" "stata rc=`stata_rc', csample rc=`ctools_rc'"

* By variable doesn't exist
* Native sample also supports by(); these cases test csample group sizes.
clear
set obs 100
gen id = _n
capture csample 50, by(nonexistent_var)
if _rc != 0 test_pass "[error] nonexistent by variable (rc=`=_rc')"
else test_fail "[error] nonexistent by variable" "should have errored"

* By-groups . and .a stay distinct: count(1) keeps one observation per group
clear
set obs 30
gen double g = cond(mod(_n, 3) == 0, ., cond(mod(_n, 3) == 1, .a, 1))
set seed 20260927
capture csample, count(1) by(g)
local rc = _rc
if `rc' == 0 & _N == 3 test_pass "count(1) by(g) with . and .a groups"
else test_fail "count(1) by(g) with . and .a groups" "rc=`rc', kept `=_N' obs (expected 3)"

* End of csample validation
* DATA-07: the by prefix and count(0), as sample
clear
set obs 100
gen g = mod(_n, 4)
sort g
set seed 9
preserve
by g: sample 10
local n_native = _N
restore, preserve
by g: csample 10
if _N == `n_native' test_pass "DATA-07 by: prefix"
else test_fail "DATA-07 by: prefix" "N=`=_N' vs `n_native'"
restore, preserve
sample 0, count
local n_native = _N
restore, preserve
capture csample, count(0)
if !_rc & _N == `n_native' test_pass "DATA-07 count(0)"
else test_fail "DATA-07 count(0)" "rc=`=_rc' N=`=_N' vs `n_native'"
restore
capture by g: csample 10, by(g)
if _rc == 190 test_pass "DATA-07 by: with by() refused (r(190))"
else test_fail "DATA-07 by: with by() refused" "r(`=_rc')"

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
* SAMP-01 to SAMP-05: sample's syntax (# with count, pduplicates()), its
* messages and return codes, sort order after by(), and version 6 draws
clear
set obs 50
gen long id = _n
gen byte g = mod(_n, 3)
gen double u = runiform()
sort u
drop u
local i = 0
foreach pair in `"3, count"' `"3, c"' `"2, count by(g)"' `"29"' `"29, by(g)"' `"35 if g == 1"' ///
    `"50, pdup(1e-6)"' `"50, pduplicates(.5) by(g)"' `"1, pdup(abc)"' `"100"' `"101"' `"-1"' ///
    `"50, count"' `"60, count"' `"5.5, count"' `"-1, count"' `"0, count by(g)"' `"abc"' ///
    `"10 20"' `"98"' `"0"' {
    local ++i
    _ct_same_output "SAMP-01 output `i'" `"sample `pair'"' `"csample `pair'"' nodata
}
* Observation counts (draws differ by design; the numbers kept do not)
foreach pair in `"3, count"' `"2, count by(g)"' `"29"' `"29, by(g)"' `"35 if g == 1"' `"-1, count"' {
    preserve
    sample `pair'
    local want = _N
    restore, preserve
    csample `pair'
    local got = _N
    restore
    if `got' == `want' test_pass "SAMP-01 kept count: `pair'"
    else test_fail "SAMP-01 kept count: `pair'" "kept `got', sample keeps `want'"
}
* SAMP-03: sample leaves the data sorted by the by() variables
preserve
sample 50, by(g)
local want : sortedby
restore, preserve
csample 50, by(g)
local got : sortedby
local ok = "`got'" == "`want'" & "`got'" == "g"
sort g, stable
restore
if `ok' test_pass "SAMP-03 sorted by by() variables"
else test_fail "SAMP-03 sorted by by() variables" "sortedby [`got'] vs [`want']"
* the kept observations keep their order within each group
preserve
gen long orig = _n
csample 50, by(g)
gen long n = _n
sort g orig
assert n == _n
restore
test_pass "SAMP-03 order within groups"
* SAMP-05: one observation is singular
clear
set obs 2
_ct_same_output "SAMP-05 one observation deleted" "sample 50" "csample 50" nodata
* SAMP-04: under version 6 sample's old algorithm draws the same sample
clear
set obs 40
gen long id = _n
gen byte g = mod(_n, 4)
sort g
foreach pair in `"50"' `"25 if g < 2"' `"30 in 5/30"' `"50, by(g)"' `"50, count"' `"101"' ///
    `"50, by(g) threads(2) verbose"' {
    tempfile want
    local npair : subinstr local pair " threads(2) verbose" ""
    preserve
    set seed 11
    capture noisily version 6: sample `npair'
    local want_rc = _rc
    quietly save `want'
    restore, preserve
    set seed 11
    capture noisily version 6: csample `pair'
    local got_rc = _rc
    capture cf _all using `want'
    local same = !_rc & _N > 0
    restore
    if `same' & `got_rc' == `want_rc' test_pass "SAMP-04 version 6 draw: `pair'"
    else test_fail "SAMP-04 version 6 draw: `pair'" "rc `got_rc' vs `want_rc', same data `same'"
}
_ct_same_output "SAMP-04 version 6 by prefix" "version 6: by g: sample 50" "version 6: by g: csample 50" nodata
_ct_same_output "SAMP-04 version 6 count" "version 6: sample 50, count" "version 6: csample 50, count" nodata
noi print_summary "csample"
}

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "csample"
    global CTOOLS_CASE_COMPLETE "csample"
    exit
}

di as error "Unknown csample case: `__case'"
exit 198
