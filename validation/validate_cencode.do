/*******************************************************************************
 * validate_cencode.do
 *
 * Comprehensive validation tests for cencode vs native Stata encode
 * Tests all options, data types, and edge cases
 *
 * cencode: C-accelerated string to numeric encoding
 * Syntax: cencode varname [if] [in], generate(newvar) [label(name) noextend verbose threads(#)]
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text ""
noi di as text "======================================================================"
noi di as text "              CENCODE VALIDATION TEST SUITE"
noi di as text "======================================================================"

/*******************************************************************************
 * Helper: Compare cencode vs encode
 ******************************************************************************/
capture program drop benchmark_encode
program define benchmark_encode
    syntax varname, testname(string) [GENerate(name) LABel(name) noextend if2(string) in2(string)]

    * Set default generate name
    if "`generate'" == "" local generate = "encoded"

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build options
    local opts "generate(`generate'_c)"
    if "`label'" != "" local opts "`opts' label(`label'_c)"
    if "`noextend'" != "" local opts "`opts' noextend"

    local opts_stata "generate(`generate'_s)"
    if "`label'" != "" local opts_stata "`opts_stata' label(`label'_s)"
    if "`noextend'" != "" local opts_stata "`opts_stata' noextend"

    * Run cencode
    capture drop `generate'_c
    capture label drop `generate'_c
    if "`label'" != "" capture label drop `label'_c
    capture cencode `varlist' `ifin', `opts'
    local rc_c = _rc

    * Run native encode
    capture drop `generate'_s
    capture label drop `generate'_s
    if "`label'" != "" capture label drop `label'_s
    capture encode `varlist' `ifin', `opts_stata'
    local rc_s = _rc

    * Check both succeeded or both failed
    if `rc_c' != `rc_s' {
        noi test_fail "`testname'" "cencode rc=`rc_c', encode rc=`rc_s'"
        exit
    }

    if `rc_c' != 0 {
        noi test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare numeric values
    quietly count if `generate'_c != `generate'_s & !missing(`generate'_c) & !missing(`generate'_s)
    local ndiff = r(N)

    * Compare missing patterns
    quietly count if missing(`generate'_c) != missing(`generate'_s)
    local nmiss_diff = r(N)

    if `ndiff' > 0 | `nmiss_diff' > 0 {
        noi test_fail "`testname'" "`ndiff' value diffs, `nmiss_diff' missing diffs"
    }
    else {
        * Also verify labels match
        local all_match = 1
        quietly levelsof `generate'_c, local(codes)
        foreach c of local codes {
            local lbl_c : label (`generate'_c) `c'
            local lbl_s : label (`generate'_s) `c'
            if `"`lbl_c'"' != `"`lbl_s'"' {
                local all_match = 0
            }
        }
        if `all_match' {
            noi test_pass "`testname'"
        }
        else {
            noi test_fail "`testname'" "label text mismatch"
        }
    }

    * Cleanup
    capture drop `generate'_c `generate'_s
    capture label drop `generate'_c `generate'_s
    if "`label'" != "" {
        capture label drop `label'_c `label'_s
    }
end

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture noisily cencode make, generate(make_code)
if _rc != 0 {
    noi test_fail "cencode plugin load" "plugin returned error `=_rc'"
    noi print_summary "cencode"
    exit 1
}
noi test_pass "cencode plugin loads and runs"
drop make_code

/*******************************************************************************
 * SECTION 2: Basic functionality (10 tests)
 ******************************************************************************/
noi print_section "Basic Functionality"

* Test 2.1: Basic encoding creates variable
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
capture confirm numeric variable make_code
if _rc == 0 {
    noi test_pass "basic encoding creates numeric variable"
}
else {
    noi test_fail "basic encoding" "variable not numeric"
}

* Test 2.2: Values start at 1
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
quietly sum make_code
if r(min) == 1 {
    noi test_pass "values start at 1"
}
else {
    noi test_fail "values start at 1" "min=`=r(min)'"
}

* Test 2.3: Max equals unique count
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
quietly tab make_code
local n_unique = r(r)
quietly sum make_code
if r(max) == `n_unique' {
    noi test_pass "max value equals unique count"
}
else {
    noi test_fail "max value" "max=`=r(max)', unique=`n_unique'"
}

* Test 2.4: Compare with encode
sysuse auto, clear
noi benchmark_encode make, testname("vs encode: auto make")

* Test 2.5: Value labels applied
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
local lblname : value label make_code
if "`lblname'" != "" {
    noi test_pass "value labels applied"
}
else {
    noi test_fail "value labels" "no label attached"
}

* Test 2.6: Labels match original strings
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
local pass = 1
forvalues i = 1/10 {
    local orig = make[`i']
    local code = make_code[`i']
    local lbl : label (make_code) `code'
    if "`orig'" != "`lbl'" {
        local pass = 0
    }
}
if `pass' {
    noi test_pass "labels match original strings"
}
else {
    noi test_fail "labels match" "mismatch found"
}
drop make_code

* Test 2.7: Single observation
clear
set obs 1
gen str20 name = "only_one"
capture drop name_code
cencode name, generate(name_code)
if name_code[1] == 1 {
    noi test_pass "single observation"
}
else {
    noi test_fail "single observation" "value=`=name_code[1]'"
}

* Test 2.8: Two observations
clear
set obs 2
gen str10 x = cond(_n == 1, "apple", "banana")
cencode x, generate(x_code)
if x_code[1] == 1 & x_code[2] == 2 {
    noi test_pass "two observations alphabetical"
}
else {
    noi test_fail "two observations" "order wrong"
}

* Test 2.9: verbose option
sysuse auto, clear
capture cencode make, generate(make_code) verbose
if _rc == 0 {
    noi test_pass "verbose option"
}
else {
    noi test_fail "verbose option" "rc=`=_rc'"
}
capture drop make_code

* Test 2.10: threads option
sysuse auto, clear
capture cencode make, generate(make_code) threads(2)
if _rc == 0 {
    noi test_pass "threads(2) option"
}
else {
    noi test_fail "threads option" "rc=`=_rc'"
}
capture drop make_code

/*******************************************************************************
 * SECTION 3: Missing values (10 tests)
 ******************************************************************************/
noi print_section "Missing Values"

* Test 3.1: Empty strings become missing
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
cencode name, generate(name_code)
count if missing(name_code) & name == ""
if r(N) == 5 {
    noi test_pass "empty strings become missing"
}
else {
    noi test_fail "empty strings" "expected 5, got `=r(N)'"
}

* Test 3.2: Non-empty strings encoded
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
cencode name, generate(name_code)
count if !missing(name_code) & name != ""
if r(N) == 15 {
    noi test_pass "non-empty strings encoded"
}
else {
    noi test_fail "non-empty strings" "expected 15, got `=r(N)'"
}

* Test 3.3: All empty strings
clear
set obs 10
gen str10 empty = ""
cencode empty, generate(empty_code)
count if missing(empty_code)
if r(N) == 10 {
    noi test_pass "all empty strings -> all missing"
}
else {
    noi test_fail "all empty" "not all missing"
}

* Test 3.4: First observation empty
clear
set obs 10
gen str20 name = "name" + string(_n)
replace name = "" in 1
cencode name, generate(name_code)
local pass = missing(name_code[1]) & !missing(name_code[2])
if `pass' {
    noi test_pass "first observation empty"
}
else {
    noi test_fail "first empty" "handling wrong"
}

* Test 3.5: Last observation empty
clear
set obs 10
gen str20 name = "name" + string(_n)
replace name = "" in 10
cencode name, generate(name_code)
local pass = missing(name_code[10]) & !missing(name_code[9])
if `pass' {
    noi test_pass "last observation empty"
}
else {
    noi test_fail "last empty" "handling wrong"
}

* Test 3.6: Alternating empty/non-empty
clear
set obs 20
gen str20 name = cond(mod(_n, 2) == 0, "even", "")
cencode name, generate(name_code)
count if missing(name_code)
if r(N) == 10 {
    noi test_pass "alternating empty/non-empty"
}
else {
    noi test_fail "alternating" "expected 10 missing"
}

* Test 3.7: 90% empty strings
clear
set obs 100
gen str20 name = cond(_n <= 10, "name" + string(_n), "")
cencode name, generate(name_code)
count if missing(name_code)
if r(N) == 90 {
    noi test_pass "90% empty strings"
}
else {
    noi test_fail "90% empty" "expected 90 missing"
}

* Test 3.8: Only one non-empty
clear
set obs 100
gen str20 name = ""
replace name = "singleton" in 50
cencode name, generate(name_code)
local pass = (name_code[50] == 1) & missing(name_code[1])
if `pass' {
    noi test_pass "only one non-empty"
}
else {
    noi test_fail "singleton" "handling wrong"
}

* Test 3.9: Empty vs single space
clear
set obs 6
gen str10 name = ""
replace name = " " in 1/3
cencode name, generate(name_code)
count if !missing(name_code)
if r(N) == 3 {
    noi test_pass "space treated as non-empty"
}
else {
    noi test_fail "space vs empty" "space not distinguished"
}

* Test 3.10: Compare missing handling with encode
clear
set obs 50
gen str20 name = cond(mod(_n, 3) == 0, "", "val" + string(_n))
noi benchmark_encode name, testname("vs encode: missing pattern")

/*******************************************************************************
 * SECTION 4: Label options (10 tests)
 ******************************************************************************/
noi print_section "Label Options"

* Test 4.1: Custom label name
sysuse auto, clear
capture label drop my_labels
cencode make, generate(make_code) label(my_labels)
capture label list my_labels
if _rc == 0 {
    noi test_pass "custom label name"
}
else {
    noi test_fail "custom label" "label not created"
}
drop make_code
capture label drop my_labels

* Test 4.2: Default label name
sysuse auto, clear
capture label drop make_encoded
cencode make, generate(make_encoded)
capture label list make_encoded
if _rc == 0 {
    noi test_pass "default label name matches varname"
}
else {
    noi test_fail "default label" "label not created"
}
drop make_encoded
capture label drop make_encoded

* Test 4.3: Label attached to variable
sysuse auto, clear
capture label drop car_labels
cencode make, generate(make_code) label(car_labels)
local attached : value label make_code
if "`attached'" == "car_labels" {
    noi test_pass "label attached to variable"
}
else {
    noi test_fail "label attached" "wrong label: `attached'"
}
drop make_code
capture label drop car_labels

* Test 4.4: noextend option
sysuse auto, clear
capture label drop new_label
cencode make, generate(make_code) noextend
if _rc == 0 {
    noi test_pass "noextend option"
}
else {
    noi test_fail "noextend" "rc=`=_rc'"
}
capture drop make_code

* Test 4.5: Multiple vars same label
clear
set obs 10
gen str10 var_a = "cat" + string(mod(_n, 3))
gen str10 var_b = "cat" + string(mod(_n, 3))
capture label drop shared
capture cencode var_a, generate(a_code) label(shared)
local rc1 = _rc
capture cencode var_b, generate(b_code) label(shared)
local rc2 = _rc
if `rc1' == 0 & `rc2' == 0 {
    noi test_pass "multiple vars same label"
}
else {
    noi test_fail "shared label" "rc1=`rc1', rc2=`rc2'"
}

* Test 4.6: Label values contiguous
sysuse auto, clear
cencode make, generate(make_code)
quietly tab make_code
local n = r(r)
quietly levelsof make_code, local(codes)
local all_present = 1
forvalues i = 1/`n' {
    local found = 0
    foreach c of local codes {
        if `c' == `i' local found = 1
    }
    if `found' == 0 local all_present = 0
}
if `all_present' {
    noi test_pass "label values contiguous 1 to N"
}
else {
    noi test_fail "contiguous" "gaps in values"
}
drop make_code

* Test 4.7: All unique strings labeled
sysuse auto, clear
cencode make, generate(make_code)
quietly levelsof make_code, local(codes)
local all_labeled = 1
foreach c of local codes {
    local lbl : label (make_code) `c'
    if "`lbl'" == "" local all_labeled = 0
}
if `all_labeled' {
    noi test_pass "all unique strings labeled"
}
else {
    noi test_fail "all labeled" "some unlabeled"
}
drop make_code

* Test 4.8: Long label name (32 chars)
clear
set obs 5
gen str10 x = "val" + string(_n)
capture label drop abcdefghijklmnopqrstuvwxyz12345
capture cencode x, generate(x_code) label(abcdefghijklmnopqrstuvwxyz12345)
if _rc == 0 {
    noi test_pass "long label name (32 chars)"
}
else {
    noi test_fail "long label name" "rc=`=_rc'"
}

* Test 4.9: Label with underscore
clear
set obs 5
gen str10 x = "val" + string(_n)
capture label drop my_label_name
capture cencode x, generate(x_code) label(my_label_name)
if _rc == 0 {
    noi test_pass "label name with underscore"
}
else {
    noi test_fail "underscore label" "rc=`=_rc'"
}

* Test 4.10: Compare labels with encode
sysuse auto, clear
noi benchmark_encode make, testname("vs encode: labels") label(test_lbl)

/*******************************************************************************
 * SECTION 5: if/in conditions (10 tests)
 ******************************************************************************/
noi print_section "if/in Conditions"

* Test 5.1: if condition subset
sysuse auto, clear
cencode make if foreign == 1, generate(make_code)
count if !missing(make_code) & foreign == 1
local n_enc = r(N)
count if missing(make_code) & foreign == 0
local n_miss = r(N)
count if foreign == 1
if `n_enc' == r(N) & `n_miss' > 0 {
    noi test_pass "if condition subset"
}
else {
    noi test_fail "if condition" "wrong encoding pattern"
}
drop make_code

* Test 5.2: in range
sysuse auto, clear
cencode make in 1/20, generate(make_code)
count if !missing(make_code) in 1/20
local in_range = r(N)
count if missing(make_code) in 21/74
local out_range = r(N)
if `in_range' == 20 & `out_range' == 54 {
    noi test_pass "in range"
}
else {
    noi test_fail "in range" "wrong counts"
}
drop make_code

* Test 5.3: Empty if result
sysuse auto, clear
cencode make if price > 100000, generate(make_code)
count if !missing(make_code)
if r(N) == 0 {
    noi test_pass "empty if result"
}
else {
    noi test_fail "empty if" "should have no matches"
}
drop make_code

* Test 5.4: if matches all
sysuse auto, clear
cencode make if price > 0, generate(make_code)
count if !missing(make_code)
if r(N) == _N {
    noi test_pass "if matches all"
}
else {
    noi test_fail "if all" "should match all"
}
drop make_code

* Test 5.5: Combined if and in
sysuse auto, clear
cencode make if foreign == 1 in 1/50, generate(make_code)
count if !missing(make_code)
if r(N) > 0 & r(N) < 50 {
    noi test_pass "combined if and in"
}
else {
    noi test_fail "combined if/in" "wrong count"
}
drop make_code

* Test 5.6: Compare if with encode
sysuse auto, clear
noi benchmark_encode make, testname("vs encode: if foreign") if2("foreign == 1")

* Test 5.7: Compare in with encode
sysuse auto, clear
noi benchmark_encode make, testname("vs encode: in 1/30") in2("1/30")

* Test 5.8: Single observation in
sysuse auto, clear
cencode make in 1/1, generate(make_code)
count if !missing(make_code)
if r(N) == 1 {
    noi test_pass "single observation in"
}
else {
    noi test_fail "single in" "wrong count"
}
drop make_code

* Test 5.9: if on another string var
clear
set obs 100
gen str10 category = "cat" + string(mod(_n, 5))
gen str10 group = cond(_n <= 50, "A", "B")
cencode category if group == "A", generate(cat_code)
count if !missing(cat_code)
if r(N) == 50 {
    noi test_pass "if on string variable"
}
else {
    noi test_fail "string if" "wrong count"
}

* Test 5.10: in from middle
sysuse auto, clear
cencode make in 30/50, generate(make_code)
count if !missing(make_code) in 30/50
if r(N) == 21 {
    noi test_pass "in from middle"
}
else {
    noi test_fail "middle in" "wrong count"
}
drop make_code

/*******************************************************************************
 * SECTION 6: Alphabetical ordering (10 tests)
 ******************************************************************************/
noi print_section "Alphabetical Ordering"

* Test 6.1: Basic alphabetical order
sysuse auto, clear
cencode make, generate(make_code)
gsort make_code
local first_by_code = make[1]
gsort make
local first_alpha = make[1]
if "`first_by_code'" == "`first_alpha'" {
    noi test_pass "basic alphabetical order"
}
else {
    noi test_fail "alphabetical" "`first_by_code' vs `first_alpha'"
}
drop make_code

* Test 6.2: Numbers sort before letters
clear
set obs 5
gen str10 x = ""
replace x = "Zebra" in 1
replace x = "Apple" in 2
replace x = "123" in 3
replace x = "banana" in 4
replace x = "99" in 5
cencode x, generate(x_code)
local lbl1 : label (x_code) 1
if substr("`lbl1'", 1, 1) >= "0" & substr("`lbl1'", 1, 1) <= "9" {
    noi test_pass "numbers sort before letters"
}
else {
    noi test_fail "number sort" "first label: `lbl1'"
}

* Test 6.3: Case sensitivity
clear
set obs 4
gen str10 x = ""
replace x = "Apple" in 1
replace x = "apple" in 2
replace x = "APPLE" in 3
replace x = "banana" in 4
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 4 {
    noi test_pass "case sensitivity preserved"
}
else {
    noi test_fail "case sensitivity" "not 4 unique"
}

* Test 6.4: Consecutive codes A-Z
clear
set obs 26
gen str1 letter = char(64 + _n)
cencode letter, generate(letter_code)
gsort letter_code
local ordered = 1
forvalues i = 1/25 {
    local this = letter[`i']
    local next = letter[`=`i'+1']
    if "`this'" > "`next'" local ordered = 0
}
if `ordered' {
    noi test_pass "consecutive codes in order"
}
else {
    noi test_fail "consecutive" "not alphabetical"
}

* Test 6.5: Mixed length strings
clear
set obs 5
gen str20 x = ""
replace x = "a" in 1
replace x = "aa" in 2
replace x = "aaa" in 3
replace x = "ab" in 4
replace x = "b" in 5
cencode x, generate(x_code)
local lbl1 : label (x_code) 1
local lbl5 : label (x_code) 5
if "`lbl1'" == "a" & "`lbl5'" == "b" {
    noi test_pass "mixed length strings sort"
}
else {
    noi test_fail "mixed length" "wrong order"
}

* Test 6.6: Whitespace affects sort
clear
set obs 4
gen str20 x = ""
replace x = "a b" in 1
replace x = "ab" in 2
replace x = " ab" in 3
replace x = "a  b" in 4
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 4 {
    noi test_pass "whitespace distinguished"
}
else {
    noi test_fail "whitespace" "not 4 unique"
}

* Test 6.7: Empty strings excluded from codes
clear
set obs 5
gen str10 x = ""
replace x = "beta" in 2
replace x = "alpha" in 4
cencode x, generate(x_code)
local lbl1 : label (x_code) 1
if "`lbl1'" == "alpha" {
    noi test_pass "empty strings not coded"
}
else {
    noi test_fail "empty strings" "wrong first label"
}

* Test 6.8: Reverse order input
clear
set obs 5
gen str10 x = ""
replace x = "e" in 1
replace x = "d" in 2
replace x = "c" in 3
replace x = "b" in 4
replace x = "a" in 5
cencode x, generate(x_code)
local lbl1 : label (x_code) 1
if "`lbl1'" == "a" {
    noi test_pass "reverse order input sorted"
}
else {
    noi test_fail "reverse order" "wrong first: `lbl1'"
}

* Test 6.9: Already sorted input
clear
set obs 5
gen str10 x = ""
replace x = "a" in 1
replace x = "b" in 2
replace x = "c" in 3
replace x = "d" in 4
replace x = "e" in 5
cencode x, generate(x_code)
local pass = 1
forvalues i = 1/5 {
    if x_code[`i'] != `i' local pass = 0
}
if `pass' {
    noi test_pass "already sorted stays sorted"
}
else {
    noi test_fail "already sorted" "order changed"
}

* Test 6.10: Compare order with encode
sysuse auto, clear
noi benchmark_encode make, testname("vs encode: alphabetical order")

/*******************************************************************************
 * SECTION 7: Large datasets (10 tests)
 ******************************************************************************/
noi print_section "Large Datasets"

* Test 7.1: 100k observations
clear
set obs 100000
gen str20 category = "cat" + string(mod(_n, 100))
cencode category, generate(cat_code)
quietly tab cat_code
if r(r) == 100 {
    noi test_pass "100k observations, 100 unique"
}
else {
    noi test_fail "100k obs" "wrong unique count"
}

* Test 7.2: 500k observations
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 500))
capture cencode category, generate(cat_code)
if _rc == 0 {
    noi test_pass "500k observations"
}
else {
    noi test_fail "500k obs" "rc=`=_rc'"
}

* Test 7.3: 10k unique values
clear
set obs 50000
gen str20 name = "unique_" + string(mod(_n, 10000))
cencode name, generate(name_code)
quietly tab name_code
if r(r) == 10000 {
    noi test_pass "10k unique values"
}
else {
    noi test_fail "10k unique" "wrong count"
}

* Test 7.4: 1M obs, 10 unique
clear
set obs 1000000
gen str10 category = "cat" + string(mod(_n, 10))
capture cencode category, generate(cat_code)
if _rc == 0 {
    quietly tab cat_code
    if r(r) == 10 {
        noi test_pass "1M obs, 10 unique"
    }
    else {
        noi test_fail "1M obs" "wrong unique count"
    }
}
else {
    noi test_fail "1M obs" "rc=`=_rc'"
}

* Test 7.5: Large with if condition
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 100))
gen byte flag = mod(_n, 2)
cencode category if flag == 1, generate(cat_code)
count if !missing(cat_code)
if r(N) == 250000 {
    noi test_pass "large with if condition"
}
else {
    noi test_fail "large if" "wrong count"
}

* Test 7.6: Large with long strings
clear
set obs 100000
gen str100 longcat = "category_with_long_name_" + string(mod(_n, 100))
capture cencode longcat, generate(long_code)
if _rc == 0 {
    noi test_pass "large with long strings"
}
else {
    noi test_fail "long strings" "rc=`=_rc'"
}

* Test 7.7: threads(4)
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 100))
capture cencode category, generate(cat_code) threads(4)
if _rc == 0 {
    noi test_pass "large with threads(4)"
}
else {
    noi test_fail "threads(4)" "rc=`=_rc'"
}

* Test 7.8: 25k unique values
clear
set obs 100000
gen str20 name = "v" + string(mod(_n, 25000))
cencode name, generate(name_code)
quietly sum name_code
if r(max) == 25000 {
    noi test_pass "25k unique values"
}
else {
    noi test_fail "25k unique" "max=`=r(max)'"
}

* Test 7.9: Sparse unique pattern
clear
set obs 100000
gen str20 category = cond(mod(_n, 1000) == 0, "unique" + string(_n), "common")
cencode category, generate(cat_code)
quietly tab cat_code
if r(r) == 101 {
    noi test_pass "sparse unique pattern"
}
else {
    noi test_fail "sparse" "expected 101 unique"
}

* Test 7.10: Large dataset, small sample
clear
set obs 1000000
gen str20 category = "cat" + string(mod(_n, 100))
cencode category in 1/100, generate(cat_code)
count if !missing(cat_code)
if r(N) == 100 {
    noi test_pass "large dataset, small sample"
}
else {
    noi test_fail "small sample" "wrong count"
}

/*******************************************************************************
 * SECTION 8: Special characters (10 tests)
 ******************************************************************************/
noi print_section "Special Characters"

* Test 8.1: Spaces in strings
clear
set obs 5
gen str30 x = ""
replace x = "hello world" in 1
replace x = "foo bar baz" in 2
replace x = "single" in 3
replace x = "multi  space" in 4
replace x = "a b c" in 5
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "spaces in strings"
}
else {
    noi test_fail "spaces" "rc=`=_rc'"
}

* Test 8.2: Leading/trailing spaces
clear
set obs 4
gen str20 x = ""
replace x = "  leading" in 1
replace x = "trailing  " in 2
replace x = "  both  " in 3
replace x = "none" in 4
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 4 {
    noi test_pass "leading/trailing spaces distinguished"
}
else {
    noi test_fail "leading/trailing" "not 4 unique"
}

* Test 8.3: Commas
clear
set obs 3
gen str30 x = ""
replace x = "one,two,three" in 1
replace x = "a,b" in 2
replace x = "no comma" in 3
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "commas in strings"
}
else {
    noi test_fail "commas" "rc=`=_rc'"
}

* Test 8.4: Dashes and underscores
clear
set obs 4
gen str30 x = ""
replace x = "with-dash" in 1
replace x = "with_underscore" in 2
replace x = "with-both_here" in 3
replace x = "plain" in 4
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "dashes and underscores"
}
else {
    noi test_fail "dash/underscore" "rc=`=_rc'"
}

* Test 8.5: Parentheses and brackets
clear
set obs 4
gen str30 x = ""
replace x = "(parens)" in 1
replace x = "[brackets]" in 2
replace x = "{braces}" in 3
replace x = "plain" in 4
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "parentheses and brackets"
}
else {
    noi test_fail "brackets" "rc=`=_rc'"
}

* Test 8.6: Ampersand and percent
clear
set obs 3
gen str30 x = ""
replace x = "A & B" in 1
replace x = "50%" in 2
replace x = "plain" in 3
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "ampersand and percent"
}
else {
    noi test_fail "amp/percent" "rc=`=_rc'"
}

* Test 8.7: Apostrophes
clear
set obs 3
gen str30 x = ""
replace x = "don't" in 1
replace x = "it's" in 2
replace x = "plain" in 3
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "apostrophes"
}
else {
    noi test_fail "apostrophes" "rc=`=_rc'"
}

* Test 8.8: Slashes
clear
set obs 3
gen str30 x = ""
replace x = "path/to/file" in 1
replace x = "win\path" in 2
replace x = "plain" in 3
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "forward and back slashes"
}
else {
    noi test_fail "slashes" "rc=`=_rc'"
}

* Test 8.9: Tab characters
clear
set obs 3
gen str30 x = ""
replace x = "col1" + char(9) + "col2" in 1
replace x = "no tab" in 2
replace x = "another" in 3
capture cencode x, generate(x_code)
if _rc == 0 {
    noi test_pass "tab characters"
}
else {
    noi test_fail "tabs" "rc=`=_rc'"
}

* Test 8.10: Compare special chars with encode
clear
set obs 5
gen str50 x = ""
replace x = "Hello! How are you?" in 1
replace x = "#hashtag @mention" in 2
replace x = "$100 + 50%" in 3
replace x = "plain" in 4
replace x = "another" in 5
noi benchmark_encode x, testname("vs encode: special chars")

/*******************************************************************************
 * SECTION 9: Real-world datasets (10 tests)
 ******************************************************************************/
noi print_section "Real-World Datasets"

* Test 9.1: auto - make
sysuse auto, clear
noi benchmark_encode make, testname("auto: make")

* Test 9.2: census - state
sysuse census, clear
noi benchmark_encode state, testname("census: state")

* Test 9.3: census - region (decode first)
sysuse census, clear
decode region, generate(region_str)
noi benchmark_encode region_str, testname("census: region")

* Test 9.4: lifeexp - country
webuse lifeexp, clear
noi benchmark_encode country, testname("lifeexp: country")

* Test 9.5: nlswork - race (tostring first)
webuse nlswork, clear
tostring race, generate(race_str)
noi benchmark_encode race_str, testname("nlswork: race")

* Test 9.6: citytemp - region
sysuse citytemp, clear
decode region, generate(region_str)
noi benchmark_encode region_str, testname("citytemp: region")

* Test 9.7: pop2000 - agegrp
sysuse pop2000, clear
decode agegrp, generate(agegrp_str)
noi benchmark_encode agegrp_str, testname("pop2000: agegrp")

* Test 9.8: voter - candidat
sysuse voter, clear
decode candidat, generate(cand_str)
noi benchmark_encode cand_str, testname("voter: candidat")

* Test 9.9: educ99gdp - country
webuse educ99gdp, clear
noi benchmark_encode country, testname("educ99gdp: country")

* Test 9.10: bpwide - patient
sysuse bpwide, clear
tostring patient, generate(patient_str)
noi benchmark_encode patient_str, testname("bpwide: patient")

/*******************************************************************************
 * SECTION 10: Pathological data (10 tests)
 ******************************************************************************/
noi print_section "Pathological Data"

* Test 10.1: All same value
clear
set obs 1000
gen str20 x = "identical"
cencode x, generate(x_code)
quietly sum x_code
if r(min) == 1 & r(max) == 1 {
    noi test_pass "all same value"
}
else {
    noi test_fail "all same" "not all 1"
}

* Test 10.2: All unique values
clear
set obs 1000
gen str20 x = "unique_" + string(_n)
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 1000 {
    noi test_pass "all unique values"
}
else {
    noi test_fail "all unique" "not 1000 unique"
}

* Test 10.3: Binary values 50/50
clear
set obs 1000
gen str10 x = cond(_n <= 500, "A", "B")
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 2 {
    noi test_pass "binary values 50/50"
}
else {
    noi test_fail "binary" "not 2 unique"
}

* Test 10.4: Very long strings (244 chars)
clear
set obs 10
gen str244 longstr = "a" * 240 + string(_n)
capture cencode longstr, generate(long_code)
if _rc == 0 {
    noi test_pass "very long strings (244 chars)"
}
else {
    noi test_fail "long strings" "rc=`=_rc'"
}

* Test 10.5: strL variables
clear
set obs 10
gen strL bigtext = "This is a longer piece of text number " + string(_n)
capture cencode bigtext, generate(big_code)
if _rc == 0 {
    noi test_pass "strL variables"
}
else {
    noi test_fail "strL" "rc=`=_rc'"
}

* Test 10.6: Bookend pattern (first and last same)
clear
set obs 100
gen str20 x = cond(_n == 1 | _n == 100, "bookend", "middle" + string(_n))
cencode x, generate(x_code)
if x_code[1] == x_code[100] {
    noi test_pass "bookend pattern"
}
else {
    noi test_fail "bookend" "first and last differ"
}

* Test 10.7: Cyclic pattern
clear
set obs 1000
gen str10 x = "cycle" + string(mod(_n - 1, 7))
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 7 {
    noi test_pass "cyclic pattern"
}
else {
    noi test_fail "cyclic" "not 7 unique"
}

* Test 10.8: Shuffled duplicates
clear
set obs 100
set seed 54321
gen str10 x = "val" + string(ceil(runiform() * 5))
gen double shuffle = runiform()
sort shuffle
drop shuffle
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 5 {
    noi test_pass "shuffled duplicates"
}
else {
    noi test_fail "shuffled" "not 5 unique"
}

* Test 10.9: Strings at str2045 boundary
clear
set obs 3
gen str2045 longstr = "a" * 2040 + string(_n)
capture cencode longstr, generate(long_code)
if _rc == 0 {
    quietly tab long_code
    if r(r) == 3 {
        noi test_pass "str2045 boundary"
    }
    else {
        noi test_fail "str2045" "not 3 unique"
    }
}
else {
    noi test_fail "str2045" "rc=`=_rc'"
}

* Test 10.10: Compare pathological with encode
clear
set obs 50
gen str20 x = cond(mod(_n, 10) == 0, "rare", "common")
noi benchmark_encode x, testname("vs encode: pathological")

/*******************************************************************************
 * SECTION 11: Error handling (5 tests)
 ******************************************************************************/
noi print_section "Error Handling"

* Test 11.1: Generate var exists
sysuse auto, clear
gen make_code = 1
capture cencode make, generate(make_code)
if _rc == 110 {
    noi test_pass "error: generate var exists"
}
else {
    noi test_fail "var exists" "expected rc=110, got `=_rc'"
}

* Test 11.2: Source is numeric
sysuse auto, clear
capture cencode price, generate(price_code)
if _rc != 0 {
    noi test_pass "error: numeric source"
}
else {
    noi test_fail "numeric source" "should error"
}

* Test 11.3: Source doesn't exist
sysuse auto, clear
capture cencode nonexistent, generate(ne_code)
if _rc != 0 {
    noi test_pass "error: nonexistent source"
}
else {
    noi test_fail "nonexistent" "should error"
}

* Test 11.4: Invalid generate name
sysuse auto, clear
capture cencode make, generate(123invalid)
if _rc != 0 {
    noi test_pass "error: invalid generate name"
}
else {
    noi test_fail "invalid name" "should error"
}

* Test 11.5: Empty dataset
clear
set obs 0
gen str10 x = ""
capture cencode x, generate(x_code)
if _rc == 0 | _rc == 2000 {
    noi test_pass "empty dataset handled"
}
else {
    noi test_fail "empty dataset" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 12: Consistency/comparison (5 tests)
 ******************************************************************************/
noi print_section "Consistency/Comparison"

* Test 12.1: cencode vs encode ordering
sysuse auto, clear
encode make, generate(make_e)
cencode make, generate(make_c)
local match = 1
forvalues i = 1/10 {
    local lbl_e : label (make_e) `=make_e[`i']'
    local lbl_c : label (make_c) `=make_c[`i']'
    if "`lbl_e'" != "`lbl_c'" local match = 0
}
if `match' {
    noi test_pass "cencode vs encode ordering"
}
else {
    noi test_fail "ordering" "labels mismatch"
}

* Test 12.2: Round-trip cencode then cdecode
sysuse auto, clear
cencode make, generate(make_code)
capture cdecode make_code, generate(make_back)
if _rc == 0 {
    local match = 1
    forvalues i = 1/10 {
        if make[`i'] != make_back[`i'] local match = 0
    }
    if `match' {
        noi test_pass "round-trip cencode/cdecode"
    }
    else {
        noi test_fail "round-trip" "mismatch"
    }
}
else {
    noi test_fail "round-trip" "cdecode failed"
}

* Test 12.3: Repeated encoding same data
clear
set obs 100
gen str20 x = "val" + string(mod(_n, 10))
capture label drop x1 x2
cencode x, generate(x1) label(x1)
cencode x, generate(x2) label(x2)
local match = 1
forvalues i = 1/20 {
    if x1[`i'] != x2[`i'] local match = 0
}
if `match' {
    noi test_pass "repeated encoding same data"
}
else {
    noi test_fail "repeated" "values differ"
}

* Test 12.4: Observation count preserved
sysuse auto, clear
local orig_n = _N
cencode make, generate(make_code)
if _N == `orig_n' {
    noi test_pass "observation count preserved"
}
else {
    noi test_fail "obs count" "count changed"
}

* Test 12.5: Compare full dataset
sysuse census, clear
noi benchmark_encode state, testname("full comparison: census state")

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cencode"

if $TESTS_FAILED > 0 {
    exit 1
}

}
