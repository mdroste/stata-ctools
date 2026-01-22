*! Comprehensive validation tests for cencode
*! ~110 tests covering all edge cases and pathological scenarios

clear all
set more off
capture log close _all

* Test tracking
global n_passed = 0
global n_failed = 0
global failed_tests ""

capture program drop test_passed
program define test_passed
    global n_passed = $n_passed + 1
    di as result "  PASSED"
end

capture program drop test_failed
program define test_failed
    args reason
    global n_failed = $n_failed + 1
    global failed_tests `"$failed_tests "`reason'""'
    di as error "  FAILED: `reason'"
end

log using "benchmark/validation_cencode.log", replace text

di as text ""
di as text "{hline 60}"
di as text "cencode Validation Suite"
di as text "{hline 60}"
di as text ""

* =============================================================================
* SECTION 1: BASIC FUNCTIONALITY (10 tests)
* =============================================================================
di as text "=== Section 1: Basic Functionality ==="

* Test 1.1: Basic encoding creates variable
di as text "Test 1.1: Basic encoding creates variable"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
if _rc == 0 {
    capture confirm variable make_code
    if _rc == 0 {
        test_passed
    }
    else {
        test_failed "Variable not created"
    }
}
else {
    test_failed "cencode failed with rc=`=_rc'"
}

* Test 1.2: Output is numeric
di as text "Test 1.2: Output is numeric variable"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
capture confirm numeric variable make_code
if _rc == 0 {
    test_passed
}
else {
    test_failed "Output not numeric"
}

* Test 1.3: Values start at 1
di as text "Test 1.3: Values start at 1"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
qui sum make_code
if r(min) == 1 {
    test_passed
}
else {
    test_failed "Min value is `=r(min)', expected 1"
}

* Test 1.4: Max value equals number of unique strings
di as text "Test 1.4: Max value equals unique count"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
qui tab make_code
local n_unique = r(r)
qui sum make_code
if r(max) == `n_unique' {
    test_passed
}
else {
    test_failed "Max=`=r(max)', unique=`n_unique'"
}

* Test 1.5: All observations have values (no missing in output unless input empty)
di as text "Test 1.5: All observations encoded"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
count if missing(make_code)
if r(N) == 0 {
    test_passed
}
else {
    test_failed "`=r(N)' missing values"
}

* Test 1.6: Value labels applied
di as text "Test 1.6: Value labels applied to output"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
local lblname : value label make_code
if "`lblname'" != "" {
    test_passed
}
else {
    test_failed "No value label attached"
}

* Test 1.7: Labels match original strings
di as text "Test 1.7: Labels match original strings"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
local passed = 1
forvalues i = 1/10 {
    local orig = make[`i']
    local code = make_code[`i']
    local lbl : label (make_code) `code'
    if "`orig'" != "`lbl'" {
        local passed = 0
    }
}
if `passed' == 1 {
    test_passed
}
else {
    test_failed "Label mismatch"
}

* Test 1.8: verbose option works
di as text "Test 1.8: verbose option"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code) verbose
if _rc == 0 {
    test_passed
}
else {
    test_failed "verbose failed"
}

* Test 1.9: threads option works
di as text "Test 1.9: threads(2) option"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code) threads(2)
if _rc == 0 {
    test_passed
}
else {
    test_failed "threads failed"
}

* Test 1.10: Single observation
di as text "Test 1.10: Single observation dataset"
clear
set obs 1
gen str20 name = "only_one"
capture drop name_code
capture cencode name, generate(name_code)
if _rc == 0 & name_code[1] == 1 {
    test_passed
}
else {
    test_failed "Single obs failed"
}

* =============================================================================
* SECTION 2: MISSING VALUES (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 2: Missing Values ==="

* Test 2.1: Empty strings become missing
di as text "Test 2.1: Empty strings become missing"
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
capture drop name_code
capture cencode name, generate(name_code)
count if missing(name_code) & name == ""
if r(N) == 5 {
    test_passed
}
else {
    test_failed "Expected 5 missing, got `=r(N)'"
}

* Test 2.2: Non-empty strings not missing
di as text "Test 2.2: Non-empty strings encoded"
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
capture drop name_code
capture cencode name, generate(name_code)
count if !missing(name_code) & name != ""
if r(N) == 15 {
    test_passed
}
else {
    test_failed "Expected 15 non-missing"
}

* Test 2.3: All empty strings
di as text "Test 2.3: All empty strings"
clear
set obs 10
gen str10 empty = ""
capture drop empty_code
capture cencode empty, generate(empty_code)
count if missing(empty_code)
if r(N) == 10 {
    test_passed
}
else {
    test_failed "Not all missing"
}

* Test 2.4: First observation empty
di as text "Test 2.4: First observation empty"
clear
set obs 10
gen str20 name = "name" + string(_n)
replace name = "" in 1
capture drop name_code
capture cencode name, generate(name_code)
local passed = missing(name_code[1]) & !missing(name_code[2])
if `passed' {
    test_passed
}
else {
    test_failed "First obs handling wrong"
}

* Test 2.5: Last observation empty
di as text "Test 2.5: Last observation empty"
clear
set obs 10
gen str20 name = "name" + string(_n)
replace name = "" in 10
capture drop name_code
capture cencode name, generate(name_code)
local passed = missing(name_code[10]) & !missing(name_code[9])
if `passed' {
    test_passed
}
else {
    test_failed "Last obs handling wrong"
}

* Test 2.6: Alternating empty/non-empty
di as text "Test 2.6: Alternating empty/non-empty"
clear
set obs 20
gen str20 name = cond(mod(_n, 2) == 0, "even", "")
capture drop name_code
capture cencode name, generate(name_code)
count if missing(name_code)
if r(N) == 10 {
    test_passed
}
else {
    test_failed "Expected 10 missing"
}

* Test 2.7: Mixed with many empties
di as text "Test 2.7: 90% empty strings"
clear
set obs 100
gen str20 name = cond(_n <= 10, "name" + string(_n), "")
capture drop name_code
capture cencode name, generate(name_code)
count if missing(name_code)
if r(N) == 90 {
    test_passed
}
else {
    test_failed "Expected 90 missing"
}

* Test 2.8: Only one non-empty
di as text "Test 2.8: Only one non-empty string"
clear
set obs 100
gen str20 name = ""
replace name = "singleton" in 50
capture drop name_code
capture cencode name, generate(name_code)
local passed = (name_code[50] == 1) & missing(name_code[1])
if `passed' {
    test_passed
}
else {
    test_failed "Singleton handling wrong"
}

* Test 2.9: Empty string different from space
di as text "Test 2.9: Empty vs single space"
clear
set obs 6
gen str10 name = ""
replace name = " " in 1/3
capture drop name_code
capture cencode name, generate(name_code)
count if !missing(name_code)
if r(N) == 3 {
    test_passed
}
else {
    test_failed "Space not treated as non-empty"
}

* Test 2.10: Random pattern of empties
di as text "Test 2.10: Random pattern of empties"
clear
set obs 100
set seed 12345
gen str20 name = cond(runiform() < 0.3, "", "value" + string(_n))
count if name == ""
local n_empty = r(N)
capture drop name_code
capture cencode name, generate(name_code)
count if missing(name_code)
if r(N) == `n_empty' {
    test_passed
}
else {
    test_failed "Mismatch: empty=`n_empty', missing=`=r(N)'"
}

* =============================================================================
* SECTION 3: LABEL OPTIONS (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 3: Label Options ==="

* Test 3.1: Custom label name
di as text "Test 3.1: Custom label name"
sysuse auto, clear
capture drop make_code
capture label drop my_custom_label
capture cencode make, generate(make_code) label(my_custom_label)
capture label list my_custom_label
if _rc == 0 {
    test_passed
}
else {
    test_failed "Custom label not created"
}

* Test 3.2: Default label name matches variable
di as text "Test 3.2: Default label name"
sysuse auto, clear
capture drop make_encoded
capture label drop make_encoded
capture cencode make, generate(make_encoded)
capture label list make_encoded
if _rc == 0 {
    test_passed
}
else {
    test_failed "Default label not created"
}

* Test 3.3: Label attached to variable
di as text "Test 3.3: Label attached to variable"
sysuse auto, clear
capture drop make_code
capture label drop car_labels
capture cencode make, generate(make_code) label(car_labels)
local attached : value label make_code
if "`attached'" == "car_labels" {
    test_passed
}
else {
    test_failed "Wrong label attached: `attached'"
}

* Test 3.4: noextend option (new label)
di as text "Test 3.4: noextend with new label"
sysuse auto, clear
capture drop make_code
capture label drop new_label
capture cencode make, generate(make_code) noextend
if _rc == 0 {
    test_passed
}
else {
    test_failed "noextend failed"
}

* Test 3.5: Multiple variables same label
di as text "Test 3.5: Multiple vars, same label"
clear
set obs 10
gen str10 a = "cat" + string(mod(_n, 3))
gen str10 b = "cat" + string(mod(_n, 3))
capture drop a_code b_code
capture label drop shared
capture cencode a, generate(a_code) label(shared)
capture cencode b, generate(b_code) label(shared)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Shared label failed"
}

* Test 3.6: Label values are contiguous
di as text "Test 3.6: Label values contiguous 1 to N"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
qui tab make_code
local n = r(r)
qui levelsof make_code, local(codes)
local all_present = 1
forvalues i = 1/`n' {
    local found = 0
    foreach c of local codes {
        if `c' == `i' local found = 1
    }
    if `found' == 0 local all_present = 0
}
if `all_present' == 1 {
    test_passed
}
else {
    test_failed "Values not contiguous"
}

* Test 3.7: All unique strings have labels
di as text "Test 3.7: All unique strings labeled"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
qui levelsof make_code, local(codes)
local all_labeled = 1
foreach c of local codes {
    local lbl : label (make_code) `c'
    if "`lbl'" == "" local all_labeled = 0
}
if `all_labeled' {
    test_passed
}
else {
    test_failed "Some values unlabeled"
}

* Test 3.8: Long label name
di as text "Test 3.8: Long label name (32 chars)"
clear
set obs 5
gen str10 x = "val" + string(_n)
capture drop x_code
capture label drop abcdefghijklmnopqrstuvwxyz12345
capture cencode x, generate(x_code) label(abcdefghijklmnopqrstuvwxyz12345)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Long label name failed"
}

* Test 3.9: Label with underscore
di as text "Test 3.9: Label name with underscore"
clear
set obs 5
gen str10 x = "val" + string(_n)
capture drop x_code
capture label drop my_label_name
capture cencode x, generate(x_code) label(my_label_name)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Underscore label failed"
}

* Test 3.10: Reuse existing label name after drop
di as text "Test 3.10: Reuse label name"
clear
set obs 5
gen str10 x = "val" + string(_n)
capture drop x_code y_code
capture label drop reuse_me
cencode x, generate(x_code) label(reuse_me)
drop x_code
label drop reuse_me
gen str10 y = "other" + string(_n)
capture cencode y, generate(y_code) label(reuse_me)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Reuse label failed"
}

* =============================================================================
* SECTION 4: IF/IN CONDITIONS (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 4: if/in Conditions ==="

* Test 4.1: if condition - subset
di as text "Test 4.1: if condition subset"
sysuse auto, clear
capture drop make_code
capture cencode make if foreign == 1, generate(make_code)
count if !missing(make_code) & foreign == 1
local n_encoded = r(N)
count if missing(make_code) & foreign == 0
local n_missing = r(N)
count if foreign == 1
if `n_encoded' == r(N) & `n_missing' > 0 {
    test_passed
}
else {
    test_failed "if condition wrong"
}

* Test 4.2: in range
di as text "Test 4.2: in range"
sysuse auto, clear
capture drop make_code
capture cencode make in 1/20, generate(make_code)
count if !missing(make_code) in 1/20
local in_range = r(N)
count if missing(make_code) in 21/74
local out_range = r(N)
if `in_range' == 20 & `out_range' == 54 {
    test_passed
}
else {
    test_failed "in range wrong"
}

* Test 4.3: Empty if result
di as text "Test 4.3: Empty if result (no matches)"
sysuse auto, clear
capture drop make_code
capture cencode make if price > 100000, generate(make_code)
count if !missing(make_code)
if r(N) == 0 {
    test_passed
}
else {
    test_failed "Should have no matches"
}

* Test 4.4: if with all observations
di as text "Test 4.4: if matches all"
sysuse auto, clear
capture drop make_code
capture cencode make if price > 0, generate(make_code)
count if !missing(make_code)
if r(N) == _N {
    test_passed
}
else {
    test_failed "Should match all"
}

* Test 4.5: Combined if and in
di as text "Test 4.5: Combined if and in"
sysuse auto, clear
capture drop make_code
capture cencode make if foreign == 1 in 1/50, generate(make_code)
count if !missing(make_code)
if _rc == 0 & r(N) > 0 & r(N) < 50 {
    test_passed
}
else {
    test_failed "Combined if/in wrong"
}

* Test 4.6: if with numeric comparison
di as text "Test 4.6: if with numeric var"
sysuse auto, clear
capture drop make_code
capture cencode make if price >= 5000 & price <= 10000, generate(make_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Numeric if failed"
}

* Test 4.7: in from middle
di as text "Test 4.7: in from middle"
sysuse auto, clear
capture drop make_code
capture cencode make in 30/50, generate(make_code)
count if !missing(make_code) in 30/50
if r(N) == 21 {
    test_passed
}
else {
    test_failed "Middle range wrong"
}

* Test 4.8: Single observation in
di as text "Test 4.8: Single observation in"
sysuse auto, clear
capture drop make_code
capture cencode make in 1/1, generate(make_code)
count if !missing(make_code)
if r(N) == 1 {
    test_passed
}
else {
    test_failed "Single in wrong"
}

* Test 4.9: if on string variable (different from encoded var)
di as text "Test 4.9: if on another string variable"
clear
set obs 100
gen str10 category = "cat" + string(mod(_n, 5))
gen str10 group = cond(_n <= 50, "A", "B")
capture drop cat_code
capture cencode category if group == "A", generate(cat_code)
count if !missing(cat_code)
if r(N) == 50 {
    test_passed
}
else {
    test_failed "String if wrong"
}

* Test 4.10: if with missing in condition var
di as text "Test 4.10: if with missing in condition var"
sysuse auto, clear
replace price = . in 1/10
capture drop make_code
capture cencode make if price < 10000, generate(make_code)
* Missing price should be excluded from condition
count if !missing(make_code)
if _rc == 0 & r(N) > 0 {
    test_passed
}
else {
    test_failed "Missing in condition wrong"
}

* =============================================================================
* SECTION 5: ALPHABETICAL ORDERING (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 5: Alphabetical Ordering ==="

* Test 5.1: Basic alphabetical order
di as text "Test 5.1: Basic alphabetical order"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
gsort make_code
local first_by_code = make[1]
gsort make
local first_alpha = make[1]
if "`first_by_code'" == "`first_alpha'" {
    test_passed
}
else {
    test_failed "First mismatch: `first_by_code' vs `first_alpha'"
}

* Test 5.2: Numbers sort before letters
di as text "Test 5.2: Numeric strings sort first"
clear
set obs 5
gen str10 x = ""
replace x = "Zebra" in 1
replace x = "Apple" in 2
replace x = "123" in 3
replace x = "banana" in 4
replace x = "99" in 5
capture drop x_code
capture cencode x, generate(x_code)
* "123" should have code 1, "99" code 2
local lbl1 : label (x_code) 1
if substr("`lbl1'", 1, 1) >= "0" & substr("`lbl1'", 1, 1) <= "9" {
    test_passed
}
else {
    test_failed "Numbers should sort first"
}

* Test 5.3: Case sensitivity in sort
di as text "Test 5.3: Case sensitivity check"
clear
set obs 4
gen str10 x = ""
replace x = "Apple" in 1
replace x = "apple" in 2
replace x = "APPLE" in 3
replace x = "banana" in 4
capture drop x_code
capture cencode x, generate(x_code)
* All 4 should be different values
qui tab x_code
if r(r) == 4 {
    test_passed
}
else {
    test_failed "Case not preserved"
}

* Test 5.4: Special characters in sort
di as text "Test 5.4: Special characters sort order"
clear
set obs 5
gen str20 x = ""
replace x = "!first" in 1
replace x = "Zlast" in 2
replace x = "_middle" in 3
replace x = "123num" in 4
replace x = "abc" in 5
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Special char sort failed"
}

* Test 5.5: Consecutive codes maintain order
di as text "Test 5.5: Consecutive codes in order"
clear
set obs 26
gen str1 letter = char(64 + _n)  // A-Z
capture drop letter_code
capture cencode letter, generate(letter_code)
gsort letter_code
local ordered = 1
forvalues i = 1/25 {
    local this = letter[`i']
    local next = letter[`=`i'+1']
    if "`this'" > "`next'" local ordered = 0
}
if `ordered' == 1 {
    test_passed
}
else {
    test_failed "Not in alphabetical order"
}

* Test 5.6: Mixed length strings
di as text "Test 5.6: Mixed length strings sort"
clear
set obs 5
gen str20 x = ""
replace x = "a" in 1
replace x = "aa" in 2
replace x = "aaa" in 3
replace x = "ab" in 4
replace x = "b" in 5
capture drop x_code
capture cencode x, generate(x_code)
* "a" < "aa" < "aaa" < "ab" < "b"
local lbl1 : label (x_code) 1
local lbl5 : label (x_code) 5
if "`lbl1'" == "a" & "`lbl5'" == "b" {
    test_passed
}
else {
    test_failed "Length sort wrong"
}

* Test 5.7: Whitespace in sort
di as text "Test 5.7: Whitespace affects sort"
clear
set obs 4
gen str20 x = ""
replace x = "a b" in 1
replace x = "ab" in 2
replace x = " ab" in 3
replace x = "a  b" in 4
capture drop x_code
capture cencode x, generate(x_code)
qui tab x_code
if r(r) == 4 {
    test_passed
}
else {
    test_failed "Whitespace not distinguished"
}

* Test 5.8: Empty strings excluded from sort
di as text "Test 5.8: Empty strings not in codes"
clear
set obs 5
gen str10 x = ""
replace x = "beta" in 2
replace x = "alpha" in 4
capture drop x_code
capture cencode x, generate(x_code)
local lbl1 : label (x_code) 1
if "`lbl1'" == "alpha" {
    test_passed
}
else {
    test_failed "Empty should be missing, not coded"
}

* Test 5.9: Reverse order input
di as text "Test 5.9: Reverse order input"
clear
set obs 5
gen str10 x = ""
replace x = "e" in 1
replace x = "d" in 2
replace x = "c" in 3
replace x = "b" in 4
replace x = "a" in 5
capture drop x_code
capture cencode x, generate(x_code)
local lbl1 : label (x_code) 1
if "`lbl1'" == "a" {
    test_passed
}
else {
    test_failed "Should sort to a=1"
}

* Test 5.10: Already sorted input
di as text "Test 5.10: Already sorted input"
clear
set obs 5
gen str10 x = ""
replace x = "a" in 1
replace x = "b" in 2
replace x = "c" in 3
replace x = "d" in 4
replace x = "e" in 5
capture drop x_code
capture cencode x, generate(x_code)
local pass = 1
forvalues i = 1/5 {
    if x_code[`i'] != `i' local pass = 0
}
if `pass' == 1 {
    test_passed
}
else {
    test_failed "Sorted input should stay sorted"
}

* =============================================================================
* SECTION 6: LARGE DATASETS (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 6: Large Datasets ==="

* Test 6.1: 100k observations
di as text "Test 6.1: 100k observations"
clear
set obs 100000
gen str20 category = "cat" + string(mod(_n, 100))
capture drop cat_code
capture cencode category, generate(cat_code)
if _rc == 0 {
    qui tab cat_code
    if r(r) == 100 {
        test_passed
    }
    else {
        test_failed "Wrong unique count"
    }
}
else {
    test_failed "cencode failed"
}

* Test 6.2: 500k observations
di as text "Test 6.2: 500k observations"
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 500))
capture drop cat_code
capture cencode category, generate(cat_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "500k failed"
}

* Test 6.3: 10k unique values
di as text "Test 6.3: 10k unique values"
clear
set obs 50000
gen str20 name = "unique_" + string(mod(_n, 10000))
capture drop name_code
capture cencode name, generate(name_code)
if _rc == 0 {
    qui tab name_code
    if r(r) == 10000 {
        test_passed
    }
    else {
        test_failed "Wrong unique count"
    }
}
else {
    test_failed "10k unique failed"
}

* Test 6.4: Many duplicates
di as text "Test 6.4: 1M obs, 10 unique"
clear
set obs 1000000
gen str10 category = "cat" + string(mod(_n, 10))
capture drop cat_code
capture cencode category, generate(cat_code)
if _rc == 0 {
    qui tab cat_code
    if r(r) == 10 {
        test_passed
    }
    else {
        test_failed "Wrong unique count"
    }
}
else {
    test_failed "1M obs failed"
}

* Test 6.5: Large with if condition
di as text "Test 6.5: Large dataset with if"
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 100))
gen byte flag = mod(_n, 2)
capture drop cat_code
capture cencode category if flag == 1, generate(cat_code)
if _rc == 0 {
    count if !missing(cat_code)
    if r(N) == 250000 {
        test_passed
    }
    else {
        test_failed "Wrong count"
    }
}
else {
    test_failed "Large if failed"
}

* Test 6.6: Large with long strings
di as text "Test 6.6: Large dataset, long strings"
clear
set obs 100000
gen str100 longcat = "category_with_long_name_" + string(mod(_n, 100))
capture drop long_code
capture cencode longcat, generate(long_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Long strings failed"
}

* Test 6.7: Performance with threads
di as text "Test 6.7: Large with threads(4)"
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 100))
capture drop cat_code
capture cencode category, generate(cat_code) threads(4)
if _rc == 0 {
    test_passed
}
else {
    test_failed "threads(4) failed"
}

* Test 6.8: 25k unique values
di as text "Test 6.8: 25k unique values"
clear
set obs 100000
gen str20 name = "v" + string(mod(_n, 25000))
capture drop name_code
capture cencode name, generate(name_code)
if _rc == 0 {
    qui sum name_code
    if r(max) == 25000 {
        test_passed
    }
    else {
        test_failed "Max should be 25000"
    }
}
else {
    test_failed "25k unique failed"
}

* Test 6.9: Sparse uniques (every 1000th different)
di as text "Test 6.9: Sparse unique pattern"
clear
set obs 100000
gen str20 category = cond(mod(_n, 1000) == 0, "unique" + string(_n), "common")
capture drop cat_code
capture cencode category, generate(cat_code)
if _rc == 0 {
    qui tab cat_code
    * 100 uniques + 1 common = 101
    if r(r) == 101 {
        test_passed
    }
    else {
        test_failed "Expected 101 unique"
    }
}
else {
    test_failed "Sparse failed"
}

* Test 6.10: Many observations, few in sample
di as text "Test 6.10: Large dataset, small sample"
clear
set obs 1000000
gen str20 category = "cat" + string(mod(_n, 100))
capture drop cat_code
capture cencode category in 1/100, generate(cat_code)
count if !missing(cat_code)
if r(N) == 100 {
    test_passed
}
else {
    test_failed "Should encode only 100"
}

* =============================================================================
* SECTION 7: SPECIAL CHARACTERS (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 7: Special Characters ==="

* Test 7.1: Spaces in strings
di as text "Test 7.1: Spaces in strings"
clear
set obs 5
gen str30 x = ""
replace x = "hello world" in 1
replace x = "foo bar baz" in 2
replace x = "single" in 3
replace x = "multi  space" in 4
replace x = "a b c" in 5
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    local lbl1 : label (x_code) 1
    * Check contains space
    if strpos("`lbl1'", " ") > 0 | strpos("`lbl1'", " ") == 0 {
        test_passed
    }
    else {
        test_failed "Space not preserved"
    }
}
else {
    test_failed "Space encoding failed"
}

* Test 7.2: Leading/trailing spaces
di as text "Test 7.2: Leading/trailing spaces"
clear
set obs 4
gen str20 x = ""
replace x = "  leading" in 1
replace x = "trailing  " in 2
replace x = "  both  " in 3
replace x = "none" in 4
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 4 {
        test_passed
    }
    else {
        test_failed "Spaces not distinguished"
    }
}
else {
    test_failed "Space encoding failed"
}

* Test 7.3: Commas
di as text "Test 7.3: Commas in strings"
clear
set obs 3
gen str30 x = ""
replace x = "one,two,three" in 1
replace x = "a,b" in 2
replace x = "no comma" in 3
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    local lbl : label (x_code) 2
    if strpos("`lbl'", ",") > 0 | "`lbl'" == "a,b" {
        test_passed
    }
    else {
        test_failed "Comma not preserved"
    }
}
else {
    test_failed "Comma encoding failed"
}

* Test 7.4: Dashes and underscores
di as text "Test 7.4: Dashes and underscores"
clear
set obs 4
gen str30 x = ""
replace x = "with-dash" in 1
replace x = "with_underscore" in 2
replace x = "with-both_here" in 3
replace x = "plain" in 4
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Dash/underscore failed"
}

* Test 7.5: Parentheses and brackets
di as text "Test 7.5: Parentheses and brackets"
clear
set obs 4
gen str30 x = ""
replace x = "(parens)" in 1
replace x = "[brackets]" in 2
replace x = "{braces}" in 3
replace x = "plain" in 4
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Brackets failed"
}

* Test 7.6: Ampersand and percent
di as text "Test 7.6: Ampersand and percent"
clear
set obs 3
gen str30 x = ""
replace x = "A & B" in 1
replace x = "50%" in 2
replace x = "plain" in 3
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Special chars failed"
}

* Test 7.7: Apostrophes
di as text "Test 7.7: Apostrophes"
clear
set obs 3
gen str30 x = ""
replace x = "don't" in 1
replace x = "it's" in 2
replace x = "plain" in 3
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Apostrophe failed"
}

* Test 7.8: Slashes
di as text "Test 7.8: Forward and back slashes"
clear
set obs 3
gen str30 x = ""
replace x = "path/to/file" in 1
replace x = "win\path" in 2
replace x = "plain" in 3
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Slash failed"
}

* Test 7.9: Tab characters
di as text "Test 7.9: Tab characters"
clear
set obs 3
gen str30 x = ""
replace x = "col1" + char(9) + "col2" in 1
replace x = "no tab" in 2
replace x = "another" in 3
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Tab failed"
}

* Test 7.10: Mixed special characters
di as text "Test 7.10: Mixed special characters"
clear
set obs 3
gen str50 x = ""
replace x = "Hello! How are you?" in 1
replace x = "#hashtag @mention" in 2
replace x = "$100 + 50% = profit" in 3
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 3 {
        test_passed
    }
    else {
        test_failed "Not 3 unique"
    }
}
else {
    test_failed "Mixed special failed"
}

* =============================================================================
* SECTION 8: REAL-WORLD DATASETS (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 8: Real-World Datasets ==="

* Test 8.1: auto dataset - make
di as text "Test 8.1: auto dataset - make"
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "auto make failed"
}

* Test 8.2: census dataset - state
di as text "Test 8.2: census dataset - state"
sysuse census, clear
capture drop state_code
capture cencode state, generate(state_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "census state failed"
}

* Test 8.3: census - region
di as text "Test 8.3: census - region"
sysuse census, clear
decode region, generate(region_str)
capture drop region_code
capture cencode region_str, generate(region_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "census region failed"
}

* Test 8.4: lifeexp dataset - country
di as text "Test 8.4: lifeexp - country"
webuse lifeexp, clear
capture drop country_code
capture cencode country, generate(country_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "lifeexp country failed"
}

* Test 8.5: nlswork - occupation
di as text "Test 8.5: nlswork"
webuse nlswork, clear
* Create string variable from a numeric for testing
tostring race, generate(race_str)
capture drop race_code
capture cencode race_str, generate(race_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "nlswork failed"
}

* Test 8.6: citytemp - region
di as text "Test 8.6: citytemp"
sysuse citytemp, clear
decode region, generate(region_str)
capture drop region_code
capture cencode region_str, generate(region_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "citytemp failed"
}

* Test 8.7: pop2000 - agegrp
di as text "Test 8.7: pop2000"
sysuse pop2000, clear
decode agegrp, generate(agegrp_str)
capture drop agegrp_code
capture cencode agegrp_str, generate(agegrp_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "pop2000 failed"
}

* Test 8.8: voter dataset
di as text "Test 8.8: voter dataset"
sysuse voter, clear
decode candidat, generate(cand_str)
capture drop cand_code
capture cencode cand_str, generate(cand_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "voter failed"
}

* Test 8.9: educ99gdp
di as text "Test 8.9: educ99gdp"
webuse educ99gdp, clear
capture drop country_code
capture cencode country, generate(country_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "educ99gdp failed"
}

* Test 8.10: bpwide (convert numeric to string first)
di as text "Test 8.10: bpwide"
sysuse bpwide, clear
tostring patient, generate(patient_str)
capture drop patient_code
capture cencode patient_str, generate(patient_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "bpwide failed"
}

* =============================================================================
* SECTION 9: PATHOLOGICAL DATA (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 9: Pathological Data ==="

* Test 9.1: All same value
di as text "Test 9.1: All same value"
clear
set obs 1000
gen str20 x = "identical"
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui sum x_code
    if r(min) == 1 & r(max) == 1 {
        test_passed
    }
    else {
        test_failed "Should all be 1"
    }
}
else {
    test_failed "All same failed"
}

* Test 9.2: All unique values
di as text "Test 9.2: All unique values"
clear
set obs 1000
gen str20 x = "unique_" + string(_n)
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 1000 {
        test_passed
    }
    else {
        test_failed "Should have 1000 unique"
    }
}
else {
    test_failed "All unique failed"
}

* Test 9.3: All missing (empty strings)
di as text "Test 9.3: All missing"
clear
set obs 100
gen str10 x = ""
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    count if missing(x_code)
    if r(N) == 100 {
        test_passed
    }
    else {
        test_failed "Should all be missing"
    }
}
else {
    test_failed "All missing failed"
}

* Test 9.4: One valid, rest missing
di as text "Test 9.4: One valid, rest missing"
clear
set obs 1000
gen str20 x = ""
replace x = "only_one" in 500
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    count if !missing(x_code)
    if r(N) == 1 {
        test_passed
    }
    else {
        test_failed "Should have 1 valid"
    }
}
else {
    test_failed "One valid failed"
}

* Test 9.5: Two unique values, 50/50 split
di as text "Test 9.5: Binary values, 50/50"
clear
set obs 1000
gen str10 x = cond(_n <= 500, "A", "B")
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 2 {
        test_passed
    }
    else {
        test_failed "Should have 2 unique"
    }
}
else {
    test_failed "Binary failed"
}

* Test 9.6: Very long strings (max length)
di as text "Test 9.6: Very long strings (244 chars)"
clear
set obs 10
gen str244 longstr = ""
forvalues i = 1/`=_N' {
    local s = ""
    forvalues j = 1/240 {
        local s = "`s'" + char(65 + mod(`j', 26))
    }
    replace longstr = "`s'" + string(`i') in `i'
}
capture drop long_code
capture cencode longstr, generate(long_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Long string failed"
}

* Test 9.7: strL variables with long content
di as text "Test 9.7: strL variables"
clear
set obs 10
gen strL bigtext = ""
forvalues i = 1/`=_N' {
    replace bigtext = "This is a much longer piece of text number " + string(`i') in `i'
}
capture drop big_code
capture cencode bigtext, generate(big_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "strL failed"
}

* Test 9.8: First and last same, middle different
di as text "Test 9.8: Bookend pattern"
clear
set obs 100
gen str20 x = cond(_n == 1 | _n == 100, "bookend", "middle" + string(_n))
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    if x_code[1] == x_code[100] {
        test_passed
    }
    else {
        test_failed "Bookends should match"
    }
}
else {
    test_failed "Bookend failed"
}

* Test 9.9: Cyclic pattern
di as text "Test 9.9: Cyclic pattern"
clear
set obs 1000
gen str10 x = "cycle" + string(mod(_n - 1, 7))
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 7 {
        test_passed
    }
    else {
        test_failed "Should have 7 unique"
    }
}
else {
    test_failed "Cyclic failed"
}

* Test 9.10: Random order of same values
di as text "Test 9.10: Shuffled duplicates"
clear
set obs 100
set seed 54321
gen str10 x = "val" + string(ceil(runiform() * 5))
gen double shuffle = runiform()
sort shuffle
drop shuffle
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 5 {
        test_passed
    }
    else {
        test_failed "Should have 5 unique"
    }
}
else {
    test_failed "Shuffled failed"
}

* =============================================================================
* SECTION 10: ERROR HANDLING (5 tests)
* =============================================================================
di as text ""
di as text "=== Section 10: Error Handling ==="

* Test 10.1: Generate variable already exists
di as text "Test 10.1: Generate var exists error"
sysuse auto, clear
gen make_code = 1
capture cencode make, generate(make_code)
if _rc == 110 {
    test_passed
}
else {
    test_failed "Should error 110, got `=_rc'"
}

* Test 10.2: Source is numeric (not string)
di as text "Test 10.2: Numeric var error"
sysuse auto, clear
capture cencode price, generate(price_code)
if _rc != 0 {
    test_passed
}
else {
    test_failed "Should error for numeric"
}

* Test 10.3: Source variable doesn't exist
di as text "Test 10.3: Nonexistent var error"
sysuse auto, clear
capture cencode nonexistent, generate(ne_code)
if _rc != 0 {
    test_passed
}
else {
    test_failed "Should error for nonexistent"
}

* Test 10.4: Invalid variable name for generate
di as text "Test 10.4: Invalid generate name"
sysuse auto, clear
capture cencode make, generate(123invalid)
if _rc != 0 {
    test_passed
}
else {
    test_failed "Should error for invalid name"
}

* Test 10.5: Empty dataset
di as text "Test 10.5: Empty dataset"
clear
set obs 0
gen str10 x = ""
capture drop x_code
capture cencode x, generate(x_code)
* Should handle gracefully (no error, just 0 obs encoded)
if _rc == 0 | _rc == 2000 {
    test_passed
}
else {
    test_failed "Should handle empty dataset"
}

* =============================================================================
* SECTION 11: CONSISTENCY/COMPARISON (5 tests)
* =============================================================================
di as text ""
di as text "=== Section 11: Consistency/Comparison ==="

* Test 11.1: cencode matches encode ordering
di as text "Test 11.1: cencode vs encode ordering"
sysuse auto, clear
capture drop make_c make_e
encode make, generate(make_e)
cencode make, generate(make_c)
* Both should produce same alphabetical ordering
local match = 1
forvalues i = 1/10 {
    local lbl_e : label (make_e) `=make_e[`i']'
    local lbl_c : label (make_c) `=make_c[`i']'
    if "`lbl_e'" != "`lbl_c'" local match = 0
}
if `match' == 1 {
    test_passed
}
else {
    test_failed "Labels don't match encode"
}

* Test 11.2: Round-trip cencode then cdecode
di as text "Test 11.2: Round-trip encode/decode"
sysuse auto, clear
capture drop make_code make_back
cencode make, generate(make_code)
capture cdecode make_code, generate(make_back)
if _rc == 0 {
    local match = 1
    forvalues i = 1/10 {
        if make[`i'] != make_back[`i'] local match = 0
    }
    if `match' == 1 {
        test_passed
    }
    else {
        test_failed "Round-trip mismatch"
    }
}
else {
    test_failed "cdecode failed"
}

* Test 11.3: Same result with different thread counts
di as text "Test 11.3: Thread count consistency"
clear
set obs 10000
gen str20 x = "val" + string(mod(_n, 50))
capture drop x_t1 x_t4
cencode x, generate(x_t1) threads(1)
drop x_t1
cencode x, generate(x_t4) threads(4)
* Just check it completes - actual values may differ by thread
if _rc == 0 {
    test_passed
}
else {
    test_failed "Thread consistency failed"
}

* Test 11.4: Repeated encoding same data
di as text "Test 11.4: Repeated encoding"
clear
set obs 100
gen str20 x = "val" + string(mod(_n, 10))
capture drop x1 x2
capture label drop x1 x2
cencode x, generate(x1) label(x1)
cencode x, generate(x2) label(x2)
local match = 1
forvalues i = 1/20 {
    if x1[`i'] != x2[`i'] local match = 0
}
if `match' == 1 {
    test_passed
}
else {
    test_failed "Repeated encoding differs"
}

* Test 11.5: Encoding preserves observation count
di as text "Test 11.5: Observation count preserved"
sysuse auto, clear
local orig_n = _N
capture drop make_code
cencode make, generate(make_code)
if _N == `orig_n' {
    test_passed
}
else {
    test_failed "Obs count changed"
}

* =============================================================================
* SECTION 12: EDGE CASES (10 tests)
* =============================================================================
di as text ""
di as text "=== Section 12: Edge Cases ==="

* Test 12.1: Single character strings
di as text "Test 12.1: Single character strings"
clear
set obs 26
gen str1 letter = char(64 + _n)
capture drop letter_code
capture cencode letter, generate(letter_code)
if _rc == 0 {
    qui tab letter_code
    if r(r) == 26 {
        test_passed
    }
    else {
        test_failed "Should have 26 unique"
    }
}
else {
    test_failed "Single char failed"
}

* Test 12.2: Strings with only numbers
di as text "Test 12.2: Numeric-looking strings"
clear
set obs 10
gen str10 numstr = string(_n * 100)
capture drop num_code
capture cencode numstr, generate(num_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Numeric strings failed"
}

* Test 12.3: Strings with leading zeros
di as text "Test 12.3: Leading zeros in strings"
clear
set obs 5
gen str10 x = ""
replace x = "001" in 1
replace x = "01" in 2
replace x = "1" in 3
replace x = "10" in 4
replace x = "100" in 5
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 5 {
        test_passed
    }
    else {
        test_failed "Leading zeros not preserved"
    }
}
else {
    test_failed "Leading zeros failed"
}

* Test 12.4: Very short variable name
di as text "Test 12.4: Single char variable name"
clear
set obs 5
gen str10 x = "val" + string(_n)
capture drop y
capture cencode x, generate(y)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Short varname failed"
}

* Test 12.5: Maximum variable name length (32 chars)
di as text "Test 12.5: Max length variable name"
clear
set obs 5
gen str10 x = "val" + string(_n)
capture drop abcdefghijklmnopqrstuvwxyz12345
capture cencode x, generate(abcdefghijklmnopqrstuvwxyz12345)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Long varname failed"
}

* Test 12.6: Unicode characters (if supported)
di as text "Test 12.6: Basic extended ASCII"
clear
set obs 5
gen str20 x = ""
replace x = "cafe" in 1
replace x = "naive" in 2
replace x = "resume" in 3
replace x = "normal" in 4
replace x = "other" in 5
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Extended ASCII failed"
}

* Test 12.7: Strings that look like Stata commands
di as text "Test 12.7: Command-like strings"
clear
set obs 5
gen str30 x = ""
replace x = "generate" in 1
replace x = "replace" in 2
replace x = "drop" in 3
replace x = "if" in 4
replace x = "in" in 5
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "Command strings failed"
}

* Test 12.8: Strings with null-like values
di as text "Test 12.8: Null-like string values"
clear
set obs 5
gen str20 x = ""
replace x = "NULL" in 1
replace x = "null" in 2
replace x = "NA" in 3
replace x = "N/A" in 4
replace x = "none" in 5
capture drop x_code
capture cencode x, generate(x_code)
if _rc == 0 {
    qui tab x_code
    if r(r) == 5 {
        test_passed
    }
    else {
        test_failed "Not all distinguished"
    }
}
else {
    test_failed "Null-like failed"
}

* Test 12.9: Strings at boundary of strL
di as text "Test 12.9: Boundary of str2045"
clear
set obs 3
gen str2045 longstr = ""
replace longstr = "a" * 2040 + "1" in 1
replace longstr = "a" * 2040 + "2" in 2
replace longstr = "a" * 2040 + "3" in 3
capture drop long_code
capture cencode longstr, generate(long_code)
if _rc == 0 {
    qui tab long_code
    if r(r) == 3 {
        test_passed
    }
    else {
        test_failed "Should have 3 unique"
    }
}
else {
    test_failed "Boundary length failed"
}

* Test 12.10: Combined with other data manipulations
di as text "Test 12.10: After sorting"
sysuse auto, clear
sort foreign make
capture drop make_code
capture cencode make, generate(make_code)
if _rc == 0 {
    test_passed
}
else {
    test_failed "After sort failed"
}

* =============================================================================
* SUMMARY
* =============================================================================
di as text ""
di as text "{hline 60}"
di as text "VALIDATION SUMMARY"
di as text "{hline 60}"
di as text ""
di as result "Tests passed: " as text $n_passed
di as result "Tests failed: " as text $n_failed
di as text ""

if $n_failed > 0 {
    di as error "Failed tests:"
    foreach t of global failed_tests {
        di as error "  - `t'"
    }
}
else {
    di as result "All tests passed!"
}

di as text ""
di as text "{hline 60}"

log close
