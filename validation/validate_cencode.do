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

noi di as text "Running validation tests for cencode..."

/*******************************************************************************
 * Helper: Compare cencode vs encode
 ******************************************************************************/
capture program drop benchmark_encode
program define benchmark_encode
    syntax varname, testname(string) [GENerate(name) LABel(name) noextend if2(string) in2(string) cencodeopts(string)]

    * Set default generate name
    if "`generate'" == "" local generate = "encoded"

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build options for cencode (includes cencodeopts like verbose, threads)
    local opts "generate(`generate'_c)"
    if "`label'" != "" local opts "`opts' label(`label'_c)"
    if "`noextend'" != "" local opts "`opts' noextend"
    if "`cencodeopts'" != "" local opts "`opts' `cencodeopts'"

    * Build options for Stata encode (no cencodeopts - Stata doesn't support them)
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
        test_fail "`testname'" "cencode rc=`rc_c', encode rc=`rc_s'"
        exit
    }

    if `rc_c' != 0 {
        test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare numeric values
    quietly count if `generate'_c != `generate'_s & !missing(`generate'_c) & !missing(`generate'_s)
    local ndiff = r(N)

    * Compare missing patterns
    quietly count if missing(`generate'_c) != missing(`generate'_s)
    local nmiss_diff = r(N)

    if `ndiff' > 0 | `nmiss_diff' > 0 {
        test_fail "`testname'" "`ndiff' value diffs, `nmiss_diff' missing diffs"
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
            test_pass "`testname'"
        }
        else {
            test_fail "`testname'" "label text mismatch"
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
print_section "Plugin Check"

sysuse auto, clear
capture noisily cencode make, generate(make_code)
if _rc != 0 {
    test_fail "cencode plugin load" "plugin returned error `=_rc'"
    print_summary "cencode"
    exit 1
}
test_pass "cencode plugin loads and runs"
drop make_code

/*******************************************************************************
 * SECTION 2: Basic functionality (10 tests)
 ******************************************************************************/
print_section "Basic Functionality"

* Test 2.1: Basic encoding creates variable
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
capture confirm numeric variable make_code
if _rc == 0 {
    test_pass "basic encoding creates numeric variable"
}
else {
    test_fail "basic encoding" "variable not numeric"
}

* Test 2.2: Values start at 1
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
quietly sum make_code
if r(min) == 1 {
    test_pass "values start at 1"
}
else {
    test_fail "values start at 1" "min=`=r(min)'"
}

* Test 2.3: Max equals unique count
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
quietly tab make_code
local n_unique = r(r)
quietly sum make_code
if r(max) == `n_unique' {
    test_pass "max value equals unique count"
}
else {
    test_fail "max value" "max=`=r(max)', unique=`n_unique'"
}

* Test 2.4: Compare with encode
sysuse auto, clear
benchmark_encode make, testname("vs encode: auto make")

* Test 2.5: Value labels applied
sysuse auto, clear
capture drop make_code
cencode make, generate(make_code)
local lblname : value label make_code
if "`lblname'" != "" {
    test_pass "value labels applied"
}
else {
    test_fail "value labels" "no label attached"
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
    test_pass "labels match original strings"
}
else {
    test_fail "labels match" "mismatch found"
}
drop make_code

* Test 2.7: Single observation
clear
set obs 1
gen str20 name = "only_one"
capture drop name_code
cencode name, generate(name_code)
if name_code[1] == 1 {
    test_pass "single observation"
}
else {
    test_fail "single observation" "value=`=name_code[1]'"
}

* Test 2.8: Two observations
clear
set obs 2
gen str10 x = cond(_n == 1, "apple", "banana")
cencode x, generate(x_code)
if x_code[1] == 1 & x_code[2] == 2 {
    test_pass "two observations alphabetical"
}
else {
    test_fail "two observations" "order wrong"
}

* Test 2.9: verbose option - verify correctness with verbose enabled
sysuse auto, clear
benchmark_encode make, testname("[syntax] verbose option") cencodeopts(verbose)

* Test 2.10: threads option - verify correctness with threads(2)
sysuse auto, clear
benchmark_encode make, testname("[syntax] threads(2) option") cencodeopts(threads(2))

/*******************************************************************************
 * SECTION 3: Missing values (10 tests)
 ******************************************************************************/
print_section "Missing Values"

* Test 3.1: Empty strings become missing
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
cencode name, generate(name_code)
count if missing(name_code) & name == ""
if r(N) == 5 {
    test_pass "empty strings become missing"
}
else {
    test_fail "empty strings" "expected 5, got `=r(N)'"
}

* Test 3.2: Non-empty strings encoded
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
cencode name, generate(name_code)
count if !missing(name_code) & name != ""
if r(N) == 15 {
    test_pass "non-empty strings encoded"
}
else {
    test_fail "non-empty strings" "expected 15, got `=r(N)'"
}

* Test 3.3: All empty strings
clear
set obs 10
gen str10 empty = ""
cencode empty, generate(empty_code)
count if missing(empty_code)
if r(N) == 10 {
    test_pass "all empty strings -> all missing"
}
else {
    test_fail "all empty" "not all missing"
}

* Test 3.4: First observation empty
clear
set obs 10
gen str20 name = "name" + string(_n)
replace name = "" in 1
cencode name, generate(name_code)
local pass = missing(name_code[1]) & !missing(name_code[2])
if `pass' {
    test_pass "first observation empty"
}
else {
    test_fail "first empty" "handling wrong"
}

* Test 3.5: Last observation empty
clear
set obs 10
gen str20 name = "name" + string(_n)
replace name = "" in 10
cencode name, generate(name_code)
local pass = missing(name_code[10]) & !missing(name_code[9])
if `pass' {
    test_pass "last observation empty"
}
else {
    test_fail "last empty" "handling wrong"
}

* Test 3.6: Alternating empty/non-empty
clear
set obs 20
gen str20 name = cond(mod(_n, 2) == 0, "even", "")
cencode name, generate(name_code)
count if missing(name_code)
if r(N) == 10 {
    test_pass "alternating empty/non-empty"
}
else {
    test_fail "alternating" "expected 10 missing"
}

* Test 3.7: 90% empty strings
clear
set obs 100
gen str20 name = cond(_n <= 10, "name" + string(_n), "")
cencode name, generate(name_code)
count if missing(name_code)
if r(N) == 90 {
    test_pass "90% empty strings"
}
else {
    test_fail "90% empty" "expected 90 missing"
}

* Test 3.8: Only one non-empty
clear
set obs 100
gen str20 name = ""
replace name = "singleton" in 50
cencode name, generate(name_code)
local pass = (name_code[50] == 1) & missing(name_code[1])
if `pass' {
    test_pass "only one non-empty"
}
else {
    test_fail "singleton" "handling wrong"
}

* Test 3.9: Empty vs single space
clear
set obs 6
gen str10 name = ""
replace name = " " in 1/3
cencode name, generate(name_code)
count if !missing(name_code)
if r(N) == 3 {
    test_pass "space treated as non-empty"
}
else {
    test_fail "space vs empty" "space not distinguished"
}

* Test 3.10: Compare missing handling with encode
clear
set obs 50
gen str20 name = cond(mod(_n, 3) == 0, "", "val" + string(_n))
benchmark_encode name, testname("vs encode: missing pattern")

/*******************************************************************************
 * SECTION 4: Label options (10 tests)
 ******************************************************************************/
print_section "Label Options"

* Test 4.1: Custom label name
sysuse auto, clear
capture label drop my_labels
cencode make, generate(make_code) label(my_labels)
capture label list my_labels
if _rc == 0 {
    test_pass "custom label name"
}
else {
    test_fail "custom label" "label not created"
}
drop make_code
capture label drop my_labels

* Test 4.2: Default label name
sysuse auto, clear
capture label drop make_encoded
cencode make, generate(make_encoded)
capture label list make_encoded
if _rc == 0 {
    test_pass "default label name matches varname"
}
else {
    test_fail "default label" "label not created"
}
drop make_encoded
capture label drop make_encoded

* Test 4.3: Label attached to variable
sysuse auto, clear
capture label drop car_labels
cencode make, generate(make_code) label(car_labels)
local attached : value label make_code
if "`attached'" == "car_labels" {
    test_pass "label attached to variable"
}
else {
    test_fail "label attached" "wrong label: `attached'"
}
drop make_code
capture label drop car_labels

* Test 4.4: noextend option
sysuse auto, clear
benchmark_encode make, testname("noextend option") noextend

* Test 4.5: Multiple vars same label
clear
set obs 10
gen str10 var_a = "cat" + string(mod(_n, 3))
gen str10 var_b = "cat" + string(mod(_n, 3))
benchmark_encode var_a, testname("multiple vars same label - var_a") label(shared)
benchmark_encode var_b, testname("multiple vars same label - var_b") label(shared2)

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
    test_pass "label values contiguous 1 to N"
}
else {
    test_fail "contiguous" "gaps in values"
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
    test_pass "all unique strings labeled"
}
else {
    test_fail "all labeled" "some unlabeled"
}
drop make_code

* Test 4.8: Long label name (32 chars)
clear
set obs 5
gen str10 x = "val" + string(_n)
benchmark_encode x, testname("long label name (32 chars)") label(abcdefghijklmnopqrstuv12345)

* Test 4.9: Label with underscore
clear
set obs 5
gen str10 x = "val" + string(_n)
benchmark_encode x, testname("label name with underscore") label(my_label_name)

* Test 4.10: Compare labels with encode
sysuse auto, clear
benchmark_encode make, testname("vs encode: labels") label(test_lbl)

/*******************************************************************************
 * SECTION 5: if/in conditions (10 tests)
 ******************************************************************************/
print_section "if/in Conditions"

* Test 5.1: if condition subset
sysuse auto, clear
cencode make if foreign == 1, generate(make_code)
count if !missing(make_code) & foreign == 1
local n_enc = r(N)
count if missing(make_code) & foreign == 0
local n_miss = r(N)
count if foreign == 1
if `n_enc' == r(N) & `n_miss' > 0 {
    test_pass "if condition subset"
}
else {
    test_fail "if condition" "wrong encoding pattern"
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
    test_pass "in range"
}
else {
    test_fail "in range" "wrong counts"
}
drop make_code

* Test 5.3: Empty if result
sysuse auto, clear
cencode make if price > 100000, generate(make_code)
count if !missing(make_code)
if r(N) == 0 {
    test_pass "empty if result"
}
else {
    test_fail "empty if" "should have no matches"
}
drop make_code

* Test 5.4: if matches all
sysuse auto, clear
cencode make if price > 0, generate(make_code)
count if !missing(make_code)
if r(N) == _N {
    test_pass "if matches all"
}
else {
    test_fail "if all" "should match all"
}
drop make_code

* Test 5.5: Compare if with encode
sysuse auto, clear
benchmark_encode make, testname("vs encode: if foreign") if2("foreign == 1")

* Test 5.7: Compare in with encode
sysuse auto, clear
benchmark_encode make, testname("vs encode: in 1/30") in2("1/30")

* Test 5.8: Single observation in
sysuse auto, clear
cencode make in 1/1, generate(make_code)
count if !missing(make_code)
if r(N) == 1 {
    test_pass "single observation in"
}
else {
    test_fail "single in" "wrong count"
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
    test_pass "if on string variable"
}
else {
    test_fail "string if" "wrong count"
}

* Test 5.10: in from middle
sysuse auto, clear
cencode make in 30/50, generate(make_code)
count if !missing(make_code) in 30/50
if r(N) == 21 {
    test_pass "in from middle"
}
else {
    test_fail "middle in" "wrong count"
}
drop make_code

/*******************************************************************************
 * SECTION 6: Alphabetical ordering (10 tests)
 ******************************************************************************/
print_section "Alphabetical Ordering"

* Test 6.1: Basic alphabetical order
sysuse auto, clear
cencode make, generate(make_code)
gsort make_code
local first_by_code = make[1]
gsort make
local first_alpha = make[1]
if "`first_by_code'" == "`first_alpha'" {
    test_pass "basic alphabetical order"
}
else {
    test_fail "alphabetical" "`first_by_code' vs `first_alpha'"
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
    test_pass "numbers sort before letters"
}
else {
    test_fail "number sort" "first label: `lbl1'"
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
    test_pass "case sensitivity preserved"
}
else {
    test_fail "case sensitivity" "not 4 unique"
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
    test_pass "consecutive codes in order"
}
else {
    test_fail "consecutive" "not alphabetical"
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
    test_pass "mixed length strings sort"
}
else {
    test_fail "mixed length" "wrong order"
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
    test_pass "whitespace distinguished"
}
else {
    test_fail "whitespace" "not 4 unique"
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
    test_pass "empty strings not coded"
}
else {
    test_fail "empty strings" "wrong first label"
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
    test_pass "reverse order input sorted"
}
else {
    test_fail "reverse order" "wrong first: `lbl1'"
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
    test_pass "already sorted stays sorted"
}
else {
    test_fail "already sorted" "order changed"
}

* Test 6.10: Compare order with encode
sysuse auto, clear
benchmark_encode make, testname("vs encode: alphabetical order")

/*******************************************************************************
 * SECTION 7: Large datasets (10 tests)
 ******************************************************************************/
print_section "Large Datasets"

* Test 7.1: 100k observations
clear
set obs 100000
gen str20 category = "cat" + string(mod(_n, 100))
cencode category, generate(cat_code)
quietly tab cat_code
if r(r) == 100 {
    test_pass "100k observations, 100 unique"
}
else {
    test_fail "100k obs" "wrong unique count"
}

* Test 7.2: 500k observations - verify correctness
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 500))
benchmark_encode category, testname("500k observations")

* Test 7.3: 10k unique values
clear
set obs 50000
gen str20 name = "unique_" + string(mod(_n, 10000))
cencode name, generate(name_code)
quietly tab name_code
if r(r) == 10000 {
    test_pass "10k unique values"
}
else {
    test_fail "10k unique" "wrong count"
}

* Test 7.4: 1M obs, 10 unique - compare on subset
clear
set obs 1000000
gen str10 category = "cat" + string(mod(_n, 10))
benchmark_encode category, testname("1M obs, 10 unique") generate(cat_code)

* Test 7.5: Large with if condition
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 100))
gen byte flag = mod(_n, 2)
benchmark_encode category, testname("large with if condition") if2("flag == 1")

* Test 7.6: Large with long strings
clear
set obs 100000
gen str100 longcat = "category_with_long_name_" + string(mod(_n, 100))
benchmark_encode longcat, testname("large with long strings") generate(long_code)

* Test 7.7: threads(4)
clear
set obs 500000
gen str20 category = "cat" + string(mod(_n, 100))
benchmark_encode category, testname("large with threads(4)") cencodeopts(threads(4))

* Test 7.8: 25k unique values
clear
set obs 100000
gen str20 name = "v" + string(mod(_n, 25000))
benchmark_encode name, testname("25k unique values")

* Test 7.9: Sparse unique pattern
clear
set obs 100000
gen str20 category = cond(mod(_n, 1000) == 0, "unique" + string(_n), "common")
benchmark_encode category, testname("sparse unique pattern")

* Test 7.10: Large dataset, small sample
clear
set obs 1000000
gen str20 category = "cat" + string(mod(_n, 100))
benchmark_encode category, testname("large dataset, small sample") in2("1/100")

/*******************************************************************************
 * SECTION 8: Special characters (10 tests)
 ******************************************************************************/
print_section "Special Characters"

* Test 8.1: Spaces in strings
clear
set obs 5
gen str30 x = ""
replace x = "hello world" in 1
replace x = "foo bar baz" in 2
replace x = "single" in 3
replace x = "multi  space" in 4
replace x = "a b c" in 5
benchmark_encode x, testname("spaces in strings")

* Test 8.2: Leading/trailing spaces
clear
set obs 4
gen str20 x = ""
replace x = "  leading" in 1
replace x = "trailing  " in 2
replace x = "  both  " in 3
replace x = "none" in 4
benchmark_encode x, testname("leading/trailing spaces distinguished")

* Test 8.3: Commas
clear
set obs 3
gen str30 x = ""
replace x = "one,two,three" in 1
replace x = "a,b" in 2
replace x = "no comma" in 3
benchmark_encode x, testname("commas in strings")

* Test 8.4: Dashes and underscores
clear
set obs 4
gen str30 x = ""
replace x = "with-dash" in 1
replace x = "with_underscore" in 2
replace x = "with-both_here" in 3
replace x = "plain" in 4
benchmark_encode x, testname("dashes and underscores")

* Test 8.5: Parentheses and brackets
clear
set obs 4
gen str30 x = ""
replace x = "(parens)" in 1
replace x = "[brackets]" in 2
replace x = "{braces}" in 3
replace x = "plain" in 4
benchmark_encode x, testname("parentheses and brackets")

* Test 8.6: Ampersand and percent
clear
set obs 3
gen str30 x = ""
replace x = "A & B" in 1
replace x = "50%" in 2
replace x = "plain" in 3
benchmark_encode x, testname("ampersand and percent")

* Test 8.7: Apostrophes
clear
set obs 3
gen str30 x = ""
replace x = "don't" in 1
replace x = "it's" in 2
replace x = "plain" in 3
benchmark_encode x, testname("apostrophes")

* Test 8.8: Slashes
clear
set obs 3
gen str30 x = ""
replace x = "path/to/file" in 1
replace x = "win\path" in 2
replace x = "plain" in 3
benchmark_encode x, testname("forward and back slashes")

* Test 8.9: Tab characters
clear
set obs 3
gen str30 x = ""
replace x = "col1" + char(9) + "col2" in 1
replace x = "no tab" in 2
replace x = "another" in 3
benchmark_encode x, testname("tab characters")

* Test 8.10: Compare special chars with encode
clear
set obs 5
gen str50 x = ""
replace x = "Hello! How are you?" in 1
replace x = "#hashtag @mention" in 2
replace x = "$100 + 50%" in 3
replace x = "plain" in 4
replace x = "another" in 5
benchmark_encode x, testname("vs encode: special chars")

/*******************************************************************************
 * SECTION 9: Real-world datasets (10 tests)
 ******************************************************************************/
print_section "Real-World Datasets"

* Test 9.1: auto - make
sysuse auto, clear
benchmark_encode make, testname("auto: make")

* Test 9.2: census - state
sysuse census, clear
benchmark_encode state, testname("census: state")

* Test 9.3: census - region (decode first)
sysuse census, clear
decode region, generate(region_str)
benchmark_encode region_str, testname("census: region")

* Test 9.4: lifeexp - country
webuse lifeexp, clear
benchmark_encode country, testname("lifeexp: country")

* Test 9.5: nlswork - race (tostring first)
webuse nlswork, clear
tostring race, generate(race_str)
benchmark_encode race_str, testname("nlswork: race")

* Test 9.6: citytemp - region
sysuse citytemp, clear
decode region, generate(region_str)
benchmark_encode region_str, testname("citytemp: region")

* Test 9.7: pop2000 - agegrp
sysuse pop2000, clear
decode agegrp, generate(agegrp_str)
benchmark_encode agegrp_str, testname("pop2000: agegrp")

* Test 9.8: voter - candidat
sysuse voter, clear
decode candidat, generate(cand_str)
benchmark_encode cand_str, testname("voter: candidat")

* Test 9.9: educ99gdp - country
webuse educ99gdp, clear
benchmark_encode country, testname("educ99gdp: country")

* Test 9.10: bpwide - patient
sysuse bpwide, clear
tostring patient, generate(patient_str)
benchmark_encode patient_str, testname("bpwide: patient")

/*******************************************************************************
 * SECTION 10: Pathological data (10 tests)
 ******************************************************************************/
print_section "Pathological Data"

* Test 10.1: All same value
clear
set obs 1000
gen str20 x = "identical"
cencode x, generate(x_code)
quietly sum x_code
if r(min) == 1 & r(max) == 1 {
    test_pass "all same value"
}
else {
    test_fail "all same" "not all 1"
}

* Test 10.2: All unique values
clear
set obs 1000
gen str20 x = "unique_" + string(_n)
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 1000 {
    test_pass "all unique values"
}
else {
    test_fail "all unique" "not 1000 unique"
}

* Test 10.3: Binary values 50/50
clear
set obs 1000
gen str10 x = cond(_n <= 500, "A", "B")
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 2 {
    test_pass "binary values 50/50"
}
else {
    test_fail "binary" "not 2 unique"
}

* Test 10.4: Very long strings (244 chars)
clear
set obs 10
gen str244 longstr = "a" * 240 + string(_n)
capture cencode longstr, generate(long_code_c)
local rc_c = _rc
capture encode longstr, generate(long_code_s)
local rc_s = _rc
if `rc_c' != 0 | `rc_s' != 0 {
    test_fail "very long strings (244 chars)" "cencode rc=`rc_c', encode rc=`rc_s'"
}
else {
    assert_var_equal long_code_c long_code_s $DEFAULT_SIGFIGS "very long strings (244 chars)"
}

* Test 10.5: strL variables
clear
set obs 10
gen strL bigtext = "This is a longer piece of text number " + string(_n)
capture cencode bigtext, generate(big_code_c)
local rc_c = _rc
capture encode bigtext, generate(big_code_s)
local rc_s = _rc
if `rc_c' != 0 | `rc_s' != 0 {
    test_fail "strL variables" "cencode rc=`rc_c', encode rc=`rc_s'"
}
else {
    assert_var_equal big_code_c big_code_s $DEFAULT_SIGFIGS "strL variables"
}

* Test 10.6: Bookend pattern (first and last same)
clear
set obs 100
gen str20 x = cond(_n == 1 | _n == 100, "bookend", "middle" + string(_n))
cencode x, generate(x_code)
if x_code[1] == x_code[100] {
    test_pass "bookend pattern"
}
else {
    test_fail "bookend" "first and last differ"
}

* Test 10.7: Cyclic pattern
clear
set obs 1000
gen str10 x = "cycle" + string(mod(_n - 1, 7))
cencode x, generate(x_code)
quietly tab x_code
if r(r) == 7 {
    test_pass "cyclic pattern"
}
else {
    test_fail "cyclic" "not 7 unique"
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
    test_pass "shuffled duplicates"
}
else {
    test_fail "shuffled" "not 5 unique"
}

* Test 10.9: Strings at str2045 boundary
clear
set obs 3
gen str2045 longstr = "a" * 2040 + string(_n)
capture cencode longstr, generate(long_code)
if _rc == 0 {
    quietly tab long_code
    if r(r) == 3 {
        test_pass "str2045 boundary"
    }
    else {
        test_fail "str2045" "not 3 unique"
    }
}
else {
    test_fail "str2045" "rc=`=_rc'"
}

* Test 10.10: Compare pathological with encode
clear
set obs 50
gen str20 x = cond(mod(_n, 10) == 0, "rare", "common")
benchmark_encode x, testname("vs encode: pathological")

/*******************************************************************************
 * SECTION 11: Error handling (5 tests)
 ******************************************************************************/
print_section "Error Handling"

* Test 11.1: Generate var exists
sysuse auto, clear
gen make_code = 1
capture cencode make, generate(make_code)
if _rc == 110 {
    test_pass "error: generate var exists"
}
else {
    test_fail "var exists" "expected rc=110, got `=_rc'"
}

* Test 11.2: Source is numeric
sysuse auto, clear
capture cencode price, generate(price_code)
if _rc != 0 {
    test_pass "error: numeric source"
}
else {
    test_fail "numeric source" "should error"
}

* Test 11.3: Source doesn't exist
sysuse auto, clear
capture cencode nonexistent, generate(ne_code)
if _rc != 0 {
    test_pass "error: nonexistent source"
}
else {
    test_fail "nonexistent" "should error"
}

* Test 11.4: Invalid generate name
sysuse auto, clear
capture cencode make, generate(123invalid)
if _rc != 0 {
    test_pass "error: invalid generate name"
}
else {
    test_fail "invalid name" "should error"
}

* Test 11.5: Empty dataset - compare with Stata's encode
clear
set obs 0
gen str10 x = ""
capture encode x, generate(x_enc)
local stata_rc = _rc
capture cencode x, generate(x_code)
local cencode_rc = _rc
if `stata_rc' == `cencode_rc' {
    test_pass "empty dataset - matches Stata behavior"
}
else {
    test_fail "empty dataset" "cencode rc=`cencode_rc' but encode rc=`stata_rc'"
}

/*******************************************************************************
 * SECTION 12: Consistency/comparison (5 tests)
 ******************************************************************************/
print_section "Consistency/Comparison"

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
    test_pass "cencode vs encode ordering"
}
else {
    test_fail "ordering" "labels mismatch"
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
        test_pass "round-trip cencode/cdecode"
    }
    else {
        test_fail "round-trip" "mismatch"
    }
}
else {
    test_fail "round-trip" "cdecode failed"
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
    test_pass "repeated encoding same data"
}
else {
    test_fail "repeated" "values differ"
}

* Test 12.4: Observation count preserved
sysuse auto, clear
local orig_n = _N
cencode make, generate(make_code)
if _N == `orig_n' {
    test_pass "observation count preserved"
}
else {
    test_fail "obs count" "count changed"
}

* Test 12.5: Compare full dataset
sysuse census, clear
benchmark_encode state, testname("full comparison: census state")

/*******************************************************************************
 * SECTION 13: Extended sysuse/webuse datasets
 ******************************************************************************/
print_section "Extended sysuse/webuse Datasets"

* Test 13.1: auto - make (string variable)
sysuse auto, clear
benchmark_encode make, testname("sysuse auto: make")

* Test 13.2: census - state (string variable)
sysuse census, clear
benchmark_encode state, testname("sysuse census: state")

* Test 13.3: census - state2 (string variable)
sysuse census, clear
benchmark_encode state2, testname("sysuse census: state2")

* Test 13.4: lifeexp - country (string variable)
webuse lifeexp, clear
benchmark_encode country, testname("webuse lifeexp: country")

* Test 13.5: nlsw88 - occupation (decode first since labeled)
sysuse nlsw88, clear
decode occupation, generate(occ_str)
benchmark_encode occ_str, testname("sysuse nlsw88: occupation")

* Test 13.6: nlsw88 - industry (decode first since labeled)
sysuse nlsw88, clear
decode industry, generate(ind_str)
benchmark_encode ind_str, testname("sysuse nlsw88: industry")

* Test 13.7: nlswork - occ (tostring since numeric)
webuse nlswork, clear
tostring occ, generate(occ_str)
benchmark_encode occ_str, testname("webuse nlswork: occ")

* Test 13.8: nlswork - ind (tostring since numeric)
webuse nlswork, clear
tostring ind, generate(ind_str)
benchmark_encode ind_str, testname("webuse nlswork: ind")

* Test 13.9: grunfeld - company (tostring since numeric)
webuse grunfeld, clear
tostring company, generate(company_str)
benchmark_encode company_str, testname("webuse grunfeld: company")

* Test 13.10: citytemp - division (decode since labeled)
sysuse citytemp, clear
decode division, generate(div_str)
benchmark_encode div_str, testname("sysuse citytemp: division")

* Test 13.11: pop2000 - agegrp (decode since labeled)
sysuse pop2000, clear
decode agegrp, generate(age_str)
benchmark_encode age_str, testname("sysuse pop2000: agegrp")

* Test 13.12: bpwide - patient (tostring since numeric)
sysuse bpwide, clear
tostring patient, generate(pat_str)
benchmark_encode pat_str, testname("sysuse bpwide: patient")

* Test 13.13: educ99gdp - country
webuse educ99gdp, clear
benchmark_encode country, testname("webuse educ99gdp: country")

/*******************************************************************************
 * SECTION 14: Empty string edge cases
 ******************************************************************************/
print_section "Empty String Edge Cases"

* Test 14.1: All empty strings
clear
set obs 50
gen str20 allEmpty = ""
cencode allEmpty, generate(code)
count if missing(code)
if r(N) == 50 {
    test_pass "all empty: all become missing"
}
else {
    test_fail "all empty" "not all missing"
}

* Test 14.2: Mixed empty/non-empty (50/50)
clear
set obs 100
gen str20 mixed = cond(mod(_n, 2) == 0, "value" + string(_n), "")
cencode mixed, generate(code)
count if missing(code) & mixed == ""
local n_miss = r(N)
count if !missing(code) & mixed != ""
local n_enc = r(N)
if `n_miss' == 50 & `n_enc' == 50 {
    test_pass "mixed 50/50: correct split"
}
else {
    test_fail "mixed 50/50" "miss=`n_miss', enc=`n_enc'"
}

* Test 14.3: First 10 empty, rest non-empty
clear
set obs 100
gen str20 firstEmpty = cond(_n <= 10, "", "val" + string(_n))
cencode firstEmpty, generate(code)
count if missing(code) in 1/10
local first_miss = r(N)
count if !missing(code) in 11/100
local rest_enc = r(N)
if `first_miss' == 10 & `rest_enc' == 90 {
    test_pass "first 10 empty"
}
else {
    test_fail "first 10 empty" "wrong pattern"
}

* Test 14.4: Last 10 empty, rest non-empty
clear
set obs 100
gen str20 lastEmpty = cond(_n > 90, "", "val" + string(_n))
cencode lastEmpty, generate(code)
count if missing(code) in 91/100
local last_miss = r(N)
count if !missing(code) in 1/90
local first_enc = r(N)
if `last_miss' == 10 & `first_enc' == 90 {
    test_pass "last 10 empty"
}
else {
    test_fail "last 10 empty" "wrong pattern"
}

* Test 14.5: Only first observation empty
clear
set obs 100
gen str20 x = "val" + string(_n)
replace x = "" in 1
cencode x, generate(code)
if missing(code[1]) & !missing(code[2]) {
    test_pass "only first empty"
}
else {
    test_fail "only first empty" "wrong pattern"
}

* Test 14.6: Only last observation empty
clear
set obs 100
gen str20 x = "val" + string(_n)
replace x = "" in 100
cencode x, generate(code)
if missing(code[100]) & !missing(code[99]) {
    test_pass "only last empty"
}
else {
    test_fail "only last empty" "wrong pattern"
}

* Test 14.7: Random sparse empty (10% empty)
clear
set obs 1000
set seed 12345
gen str20 x = cond(runiform() < 0.1, "", "val" + string(_n))
cencode x, generate(code)
count if missing(code) & x == ""
local correct_miss = r(N)
count if !missing(code) & x != ""
local correct_enc = r(N)
count if x == ""
local total_empty = r(N)
if `correct_miss' == `total_empty' {
    test_pass "10% sparse empty"
}
else {
    test_fail "sparse empty" "miss pattern wrong"
}

* Test 14.8: Compare with encode - mixed empty
clear
set obs 50
gen str20 x = cond(mod(_n, 5) == 0, "", "val" + string(_n))
benchmark_encode x, testname("vs encode: mixed empty")

/*******************************************************************************
 * SECTION 15: Single character strings
 ******************************************************************************/
print_section "Single Character Strings"

* Test 15.1: Single character A-Z
clear
set obs 26
gen str1 letter = char(64 + _n)
cencode letter, generate(code)
quietly tab code
if r(r) == 26 {
    test_pass "26 single letters A-Z"
}
else {
    test_fail "A-Z" "not 26 unique"
}

* Test 15.2: Single character a-z
clear
set obs 26
gen str1 letter = char(96 + _n)
cencode letter, generate(code)
quietly tab code
if r(r) == 26 {
    test_pass "26 single letters a-z"
}
else {
    test_fail "a-z" "not 26 unique"
}

* Test 15.3: Single digits 0-9
clear
set obs 10
gen str1 digit = string(_n - 1)
cencode digit, generate(code)
quietly tab code
if r(r) == 10 {
    test_pass "10 single digits 0-9"
}
else {
    test_fail "0-9" "not 10 unique"
}

* Test 15.4: Single special characters
clear
set obs 10
gen str1 special = ""
replace special = "!" in 1
replace special = "@" in 2
replace special = "#" in 3
replace special = "$" in 4
replace special = "%" in 5
replace special = "^" in 6
replace special = "&" in 7
replace special = "*" in 8
replace special = "(" in 9
replace special = ")" in 10
* Run Stata's encode as reference
capture encode special, generate(stata_result)
local stata_rc = _rc
* Run cencode
capture cencode special, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "10 single special chars" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "10 single special chars (both error rc=`stata_rc')"
}
else {
    * Compare results
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "10 single special chars"
    }
    else {
        test_fail "10 single special chars" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 15.5: Single space character
clear
set obs 10
gen str1 sp = " "
cencode sp, generate(code)
quietly sum code
if r(min) == 1 & r(max) == 1 {
    test_pass "single space - all same code"
}
else {
    test_fail "single space" "not all same"
}

* Test 15.6: Mix of single chars with duplicates
clear
set obs 100
gen str1 x = char(65 + mod(_n - 1, 5))
cencode x, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "5 repeated single chars"
}
else {
    test_fail "repeated chars" "not 5 unique"
}

* Test 15.7: Compare single chars with encode
clear
set obs 26
gen str1 letter = char(64 + _n)
benchmark_encode letter, testname("vs encode: single chars A-Z")

/*******************************************************************************
 * SECTION 16: Very long strings
 ******************************************************************************/
print_section "Very Long Strings"

* Test 16.1: 100 character strings
clear
set obs 10
gen str100 long100 = "a" * 95 + string(_n)
capture encode long100, generate(stata_result)
local stata_rc = _rc
capture cencode long100, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "100 char strings" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "100 char strings (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "100 char strings"
    }
    else {
        test_fail "100 char strings" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.2: 200 character strings
clear
set obs 10
gen str200 long200 = "b" * 195 + string(_n)
capture encode long200, generate(stata_result)
local stata_rc = _rc
capture cencode long200, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "200 char strings" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "200 char strings (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "200 char strings"
    }
    else {
        test_fail "200 char strings" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.3: Maximum str244 strings
clear
set obs 5
gen str244 max244 = "c" * 240 + string(_n)
capture encode max244, generate(stata_result)
local stata_rc = _rc
capture cencode max244, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "244 char strings (max str244)" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "244 char strings (max str244) (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "244 char strings (max str244)"
    }
    else {
        test_fail "244 char strings (max str244)" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.4: strL with 500 characters
clear
set obs 5
gen strL longL = "d" * 495 + string(_n)
capture encode longL, generate(stata_result)
local stata_rc = _rc
capture cencode longL, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strL 500 char" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strL 500 char (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strL 500 char"
    }
    else {
        test_fail "strL 500 char" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.5: strL with 1000 characters
clear
set obs 5
gen strL longL2 = "e" * 995 + string(_n)
capture encode longL2, generate(stata_result)
local stata_rc = _rc
capture cencode longL2, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strL 1000 char" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strL 1000 char (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strL 1000 char"
    }
    else {
        test_fail "strL 1000 char" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.6: strL with 2000 characters
clear
set obs 3
gen strL longL3 = "f" * 1995 + string(_n)
capture encode longL3, generate(stata_result)
local stata_rc = _rc
capture cencode longL3, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strL 2000 char" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strL 2000 char (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strL 2000 char"
    }
    else {
        test_fail "strL 2000 char" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.7: Very long strings differing only at end
clear
set obs 5
gen str244 diffEnd = "x" * 240
replace diffEnd = diffEnd + "1" in 1
replace diffEnd = diffEnd + "2" in 2
replace diffEnd = diffEnd + "3" in 3
replace diffEnd = diffEnd + "4" in 4
replace diffEnd = diffEnd + "5" in 5
capture encode diffEnd, generate(stata_result)
local stata_rc = _rc
capture cencode diffEnd, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "long strings differing at end" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "long strings differing at end (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "long strings differing at end"
    }
    else {
        test_fail "long strings differing at end" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.8: Very long strings differing only at start
clear
set obs 5
gen str244 diffStart = ""
replace diffStart = "1" + "y" * 240 in 1
replace diffStart = "2" + "y" * 240 in 2
replace diffStart = "3" + "y" * 240 in 3
replace diffStart = "4" + "y" * 240 in 4
replace diffStart = "5" + "y" * 240 in 5
capture encode diffStart, generate(stata_result)
local stata_rc = _rc
capture cencode diffStart, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "long strings differing at start" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "long strings differing at start (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "long strings differing at start"
    }
    else {
        test_fail "long strings differing at start" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.9: Mixed long and short strings
clear
set obs 10
gen str244 mixed = cond(mod(_n, 2) == 0, "short" + string(_n), "z" * 200 + string(_n))
capture encode mixed, generate(stata_result)
local stata_rc = _rc
capture cencode mixed, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "mixed long and short strings" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "mixed long and short strings (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "mixed long and short strings"
    }
    else {
        test_fail "mixed long and short strings" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 16.10: Compare long strings with encode
clear
set obs 10
gen str100 x = "a" * 95 + string(_n)
benchmark_encode x, testname("vs encode: 100 char strings")

/*******************************************************************************
 * SECTION 17: Strings with special characters
 ******************************************************************************/
print_section "Strings with Special Characters"

* Test 17.1: Strings with double quotes
clear
set obs 5
gen str50 withQuotes = ""
replace withQuotes = `"He said "hello""' in 1
replace withQuotes = `"She replied "goodbye""' in 2
replace withQuotes = `"Quote: "test""' in 3
replace withQuotes = `"Multiple "quotes" here"' in 4
replace withQuotes = "no quotes" in 5
capture encode withQuotes, generate(stata_result)
local stata_rc = _rc
capture cencode withQuotes, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with double quotes" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with double quotes (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with double quotes"
    }
    else {
        test_fail "strings with double quotes" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.2: Strings with commas
clear
set obs 5
gen str50 withCommas = ""
replace withCommas = "Smith, John" in 1
replace withCommas = "one, two, three" in 2
replace withCommas = "a,b,c,d,e" in 3
replace withCommas = "comma,at,end," in 4
replace withCommas = "no comma" in 5
capture encode withCommas, generate(stata_result)
local stata_rc = _rc
capture cencode withCommas, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with commas" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with commas (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with commas"
    }
    else {
        test_fail "strings with commas" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.3: Strings with tabs
clear
set obs 5
gen str50 withTabs = ""
replace withTabs = "col1" + char(9) + "col2" in 1
replace withTabs = "a" + char(9) + "b" + char(9) + "c" in 2
replace withTabs = char(9) + "leading tab" in 3
replace withTabs = "trailing tab" + char(9) in 4
replace withTabs = "no tab" in 5
capture encode withTabs, generate(stata_result)
local stata_rc = _rc
capture cencode withTabs, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with tabs" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with tabs (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with tabs"
    }
    else {
        test_fail "strings with tabs" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.4: Strings with newlines
clear
set obs 5
gen str50 withNewlines = ""
replace withNewlines = "line1" + char(10) + "line2" in 1
replace withNewlines = "a" + char(13) + "b" in 2
replace withNewlines = "cr" + char(13) + char(10) + "lf" in 3
replace withNewlines = char(10) + "leading" in 4
replace withNewlines = "no newline" in 5
capture encode withNewlines, generate(stata_result)
local stata_rc = _rc
capture cencode withNewlines, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with newlines" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with newlines (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with newlines"
    }
    else {
        test_fail "strings with newlines" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.5: Strings with apostrophes
clear
set obs 5
gen str50 withApos = ""
replace withApos = "don't" in 1
replace withApos = "it's" in 2
replace withApos = "O'Brien" in 3
replace withApos = "couldn't've" in 4
replace withApos = "no apostrophe" in 5
capture encode withApos, generate(stata_result)
local stata_rc = _rc
capture cencode withApos, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with apostrophes" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with apostrophes (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with apostrophes"
    }
    else {
        test_fail "strings with apostrophes" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.6: Strings with backslashes
clear
set obs 5
gen str50 withSlash = ""
replace withSlash = "path\to\file" in 1
replace withSlash = "C:\Windows\System32" in 2
replace withSlash = "\\server\share" in 3
replace withSlash = "escape\\n" in 4
replace withSlash = "no backslash" in 5
capture encode withSlash, generate(stata_result)
local stata_rc = _rc
capture cencode withSlash, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with backslashes" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with backslashes (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with backslashes"
    }
    else {
        test_fail "strings with backslashes" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.7: Strings with various punctuation
clear
set obs 10
gen str50 punct = ""
replace punct = "hello!" in 1
replace punct = "what?" in 2
replace punct = "stop." in 3
replace punct = "item: value" in 4
replace punct = "semi;colon" in 5
replace punct = "50%" in 6
replace punct = "$100" in 7
replace punct = "@mention" in 8
replace punct = "#hashtag" in 9
replace punct = "a & b" in 10
capture encode punct, generate(stata_result)
local stata_rc = _rc
capture cencode punct, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "various punctuation" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "various punctuation (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "various punctuation"
    }
    else {
        test_fail "various punctuation" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.8: Strings with brackets
clear
set obs 6
gen str50 brackets = ""
replace brackets = "(parens)" in 1
replace brackets = "[square]" in 2
replace brackets = "{curly}" in 3
replace brackets = "<angle>" in 4
replace brackets = "((nested))" in 5
replace brackets = "no brackets" in 6
capture encode brackets, generate(stata_result)
local stata_rc = _rc
capture cencode brackets, generate(ctools_result)
local ctools_rc = _rc
if `stata_rc' != `ctools_rc' {
    test_fail "strings with brackets" "rc differ: stata=`stata_rc' cencode=`ctools_rc'"
}
else if `stata_rc' != 0 {
    test_pass "strings with brackets (both error rc=`stata_rc')"
}
else {
    quietly count if stata_result != ctools_result & !missing(stata_result) & !missing(ctools_result)
    local ndiff = r(N)
    quietly count if missing(stata_result) != missing(ctools_result)
    local nmiss = r(N)
    if `ndiff' == 0 & `nmiss' == 0 {
        test_pass "strings with brackets"
    }
    else {
        test_fail "strings with brackets" "`ndiff' values differ, `nmiss' missing pattern mismatches"
    }
    capture drop stata_result ctools_result
}

* Test 17.9: Compare special chars with encode
clear
set obs 5
gen str50 x = ""
replace x = "A,B,C" in 1
replace x = "X;Y;Z" in 2
replace x = "a	b	c" in 3
replace x = "1|2|3" in 4
replace x = "plain" in 5
benchmark_encode x, testname("vs encode: special chars mix")

/*******************************************************************************
 * SECTION 18: Leading/trailing spaces
 ******************************************************************************/
print_section "Leading/Trailing Spaces"

* Test 18.1: Leading spaces
clear
set obs 5
gen str20 lead = ""
replace lead = " one" in 1
replace lead = "  two" in 2
replace lead = "   three" in 3
replace lead = "    four" in 4
replace lead = "five" in 5
cencode lead, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "leading spaces distinguished"
}
else {
    test_fail "leading spaces" "not 5 unique"
}

* Test 18.2: Trailing spaces
clear
set obs 5
gen str20 trail = ""
replace trail = "one " in 1
replace trail = "two  " in 2
replace trail = "three   " in 3
replace trail = "four    " in 4
replace trail = "five" in 5
cencode trail, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "trailing spaces distinguished"
}
else {
    test_fail "trailing spaces" "not 5 unique"
}

* Test 18.3: Both leading and trailing
clear
set obs 4
gen str20 both = ""
replace both = " x " in 1
replace both = "  x  " in 2
replace both = "x" in 3
replace both = "   x   " in 4
cencode both, generate(code)
quietly tab code
if r(r) == 4 {
    test_pass "leading and trailing distinguished"
}
else {
    test_fail "both spaces" "not 4 unique"
}

* Test 18.4: Strings differing only by whitespace
clear
set obs 5
gen str20 ws = ""
replace ws = "abc" in 1
replace ws = " abc" in 2
replace ws = "abc " in 3
replace ws = " abc " in 4
replace ws = "a bc" in 5
cencode ws, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "differ only by whitespace"
}
else {
    test_fail "whitespace diff" "not 5 unique"
}

* Test 18.5: Only whitespace strings
clear
set obs 5
gen str10 onlyWs = ""
replace onlyWs = " " in 1
replace onlyWs = "  " in 2
replace onlyWs = "   " in 3
replace onlyWs = char(9) in 4
replace onlyWs = " " + char(9) + " " in 5
cencode onlyWs, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "only whitespace strings distinguished"
}
else {
    test_fail "only whitespace" "not 5 unique"
}

* Test 18.6: Compare whitespace with encode
clear
set obs 4
gen str20 x = ""
replace x = "test" in 1
replace x = " test" in 2
replace x = "test " in 3
replace x = " test " in 4
benchmark_encode x, testname("vs encode: whitespace variations")

/*******************************************************************************
 * SECTION 19: Numeric-looking strings
 ******************************************************************************/
print_section "Numeric-Looking Strings"

* Test 19.1: Leading zeros
clear
set obs 5
gen str10 zeros = ""
replace zeros = "001" in 1
replace zeros = "01" in 2
replace zeros = "1" in 3
replace zeros = "0001" in 4
replace zeros = "00001" in 5
cencode zeros, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "leading zeros distinguished"
}
else {
    test_fail "leading zeros" "not 5 unique"
}

* Test 19.2: Zip codes with leading zeros
clear
set obs 5
gen str5 zip = ""
replace zip = "01234" in 1
replace zip = "00501" in 2
replace zip = "00101" in 3
replace zip = "00001" in 4
replace zip = "90210" in 5
cencode zip, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "zip codes preserved"
}
else {
    test_fail "zip codes" "not 5 unique"
}

* Test 19.3: Decimal number strings
clear
set obs 5
gen str10 decs = ""
replace decs = "1.5" in 1
replace decs = "1.50" in 2
replace decs = "1.500" in 3
replace decs = "01.5" in 4
replace decs = "1.05" in 5
cencode decs, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "decimal strings distinguished"
}
else {
    test_fail "decimal strings" "not 5 unique"
}

* Test 19.4: Negative number strings
clear
set obs 5
gen str10 negs = ""
replace negs = "-1" in 1
replace negs = "-01" in 2
replace negs = "-001" in 3
replace negs = "-1.0" in 4
replace negs = "-1.00" in 5
cencode negs, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "negative strings distinguished"
}
else {
    test_fail "negative strings" "not 5 unique"
}

* Test 19.5: Scientific notation strings
clear
set obs 5
gen str15 sci = ""
replace sci = "1e5" in 1
replace sci = "1E5" in 2
replace sci = "1.0e5" in 3
replace sci = "1e+5" in 4
replace sci = "1e-5" in 5
cencode sci, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "scientific notation distinguished"
}
else {
    test_fail "scientific notation" "not 5 unique"
}

* Test 19.6: Phone numbers
clear
set obs 5
gen str15 phone = ""
replace phone = "555-1234" in 1
replace phone = "555-1234" in 2
replace phone = "(555) 123-4567" in 3
replace phone = "5551234567" in 4
replace phone = "+1-555-123-4567" in 5
cencode phone, generate(code)
quietly tab code
if r(r) == 4 {
    test_pass "phone numbers (4 unique)"
}
else {
    test_fail "phone numbers" "not 4 unique"
}

* Test 19.7: Compare numeric strings with encode
clear
set obs 5
gen str10 x = ""
replace x = "001" in 1
replace x = "01" in 2
replace x = "1" in 3
replace x = "1.0" in 4
replace x = "1.00" in 5
benchmark_encode x, testname("vs encode: numeric strings")

/*******************************************************************************
 * SECTION 20: Unicode/international characters
 ******************************************************************************/
print_section "Unicode/International Characters"

* Test 20.1: Basic Latin Extended
clear
set obs 5
gen str50 latin = ""
replace latin = "Cafe" in 1
replace latin = "Jose" in 2
replace latin = "Francois" in 3
replace latin = "Munchen" in 4
replace latin = "Zurich" in 5
capture cencode latin, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 5 {
        test_pass "basic Latin extended"
    }
    else {
        test_fail "Latin extended" "not 5 unique"
    }
}
else {
    test_fail "Latin extended" "rc=`=_rc'"
}

* Test 20.2: Mixed ASCII and extended
clear
set obs 5
gen str50 mixed = ""
replace mixed = "normal" in 1
replace mixed = "Cafe" in 2
replace mixed = "test123" in 3
replace mixed = "Jose Garcia" in 4
replace mixed = "plain text" in 5
capture cencode mixed, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 5 {
        test_pass "mixed ASCII and extended"
    }
    else {
        test_fail "mixed ASCII" "not 5 unique"
    }
}
else {
    test_fail "mixed ASCII" "rc=`=_rc'"
}

* Test 20.3: Compare international with encode
clear
set obs 5
gen str50 x = ""
replace x = "Garcia" in 1
replace x = "Muller" in 2
replace x = "Jensen" in 3
replace x = "normal" in 4
replace x = "plain" in 5
benchmark_encode x, testname("vs encode: international chars")

/*******************************************************************************
 * SECTION 21: Case sensitivity tests
 ******************************************************************************/
print_section "Case Sensitivity"

* Test 21.1: Same word different cases
clear
set obs 5
gen str10 cases = ""
replace cases = "Apple" in 1
replace cases = "apple" in 2
replace cases = "APPLE" in 3
replace cases = "aPpLe" in 4
replace cases = "ApPlE" in 5
cencode cases, generate(code)
quietly tab code
if r(r) == 5 {
    test_pass "case variations distinguished"
}
else {
    test_fail "case variations" "not 5 unique"
}

* Test 21.2: Uppercase vs lowercase letters
clear
set obs 52
gen str1 letter = ""
replace letter = char(65 + _n - 1) in 1/26
replace letter = char(97 + _n - 27) in 27/52
cencode letter, generate(code)
quietly tab code
if r(r) == 52 {
    test_pass "52 letters (A-Z, a-z)"
}
else {
    test_fail "52 letters" "not 52 unique"
}

* Test 21.3: Case difference with duplicates
clear
set obs 100
gen str10 x = cond(mod(_n, 4) == 0, "UPPER", cond(mod(_n, 4) == 1, "lower", cond(mod(_n, 4) == 2, "Mixed", "mIxEd")))
cencode x, generate(code)
quietly tab code
if r(r) == 4 {
    test_pass "4 case variations with duplicates"
}
else {
    test_fail "case duplicates" "not 4 unique"
}

* Test 21.4: Case sort order - compare with Stata's encode
clear
set obs 4
gen str10 x = ""
replace x = "Zebra" in 1
replace x = "apple" in 2
replace x = "APPLE" in 3
replace x = "zebra" in 4
encode x, generate(code_stata)
cencode x, generate(code_ctools)
gsort code_stata
local stata_first = x[1]
gsort code_ctools
local ctools_first = x[1]
if "`stata_first'" == "`ctools_first'" {
    test_pass "case sort order matches Stata"
}
else {
    test_fail "case sort order" "cencode first=`ctools_first' but encode first=`stata_first'"
}

* Test 21.5: Compare case sensitivity with encode
clear
set obs 4
gen str10 x = ""
replace x = "Test" in 1
replace x = "test" in 2
replace x = "TEST" in 3
replace x = "tEsT" in 4
benchmark_encode x, testname("vs encode: case sensitivity")

/*******************************************************************************
 * SECTION 22: Label edge cases
 ******************************************************************************/
print_section "Label Edge Cases"

* Test 22.1: Using existing label with noextend
clear
set obs 10
gen str10 x = cond(mod(_n, 2) == 0, "even", "odd")
capture label drop mylab
label define mylab 1 "odd" 2 "even"
cencode x, generate(code) label(mylab) noextend
local lbl_1 : label mylab 1
local lbl_2 : label mylab 2
if "`lbl_1'" == "odd" & "`lbl_2'" == "even" {
    test_pass "noextend uses existing label"
}
else {
    test_fail "noextend label" "expected 1=odd,2=even but got 1=`lbl_1',2=`lbl_2'"
}
capture label drop mylab

* Test 22.2: Extending existing labels
clear
set obs 15
gen str10 x = "cat" + string(mod(_n, 5) + 1)
capture label drop extlab
label define extlab 1 "cat1" 2 "cat2"
capture cencode x, generate(code) label(extlab)
if _rc == 0 {
    * Verify all 5 categories are encoded
    quietly tab code
    if r(r) == 5 {
        test_pass "extending existing labels"
    }
    else {
        test_fail "extending labels" "expected 5 categories, got `=r(r)'"
    }
}
else {
    test_fail "extending labels" "cencode failed with rc=`=_rc'"
}
capture label drop extlab

* Test 22.3: Very long label text (244 chars in value) - compare with Stata
clear
set obs 3
gen str244 x = ""
replace x = "a" * 240 in 1
replace x = "b" * 240 in 2
replace x = "c" * 240 in 3
capture encode x, generate(code_stata)
local stata_rc = _rc
local stata_len = 0
if `stata_rc' == 0 {
    local stata_lbl : label (code_stata) 1
    local stata_len = strlen("`stata_lbl'")
}
capture cencode x, generate(code_ctools)
local cencode_rc = _rc
if `stata_rc' == `cencode_rc' {
    if `cencode_rc' == 0 {
        local ctools_lbl : label (code_ctools) 1
        local ctools_len = strlen("`ctools_lbl'")
        if `ctools_len' == `stata_len' {
            test_pass "very long label text matches Stata"
        }
        else {
            test_fail "long label text" "cencode len=`ctools_len' but encode len=`stata_len'"
        }
    }
    else {
        test_pass "very long label text - both error (as expected)"
    }
}
else {
    test_fail "long label text" "cencode rc=`cencode_rc' but encode rc=`stata_rc'"
}

* Test 22.4: 100 unique values
clear
set obs 100
gen str15 x = "category_" + string(_n)
capture cencode x, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 100 {
        test_pass "100 unique values encoded"
    }
    else {
        test_fail "100 unique" "not 100"
    }
}
else {
    test_fail "100 unique" "rc=`=_rc'"
}

* Test 22.5: 500 unique values
clear
set obs 500
gen str15 x = "cat_" + string(_n)
capture cencode x, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 500 {
        test_pass "500 unique values encoded"
    }
    else {
        test_fail "500 unique" "not 500"
    }
}
else {
    test_fail "500 unique" "rc=`=_rc'"
}

* Test 22.6: 1000 unique values
clear
set obs 1000
gen str15 x = "v" + string(_n)
capture cencode x, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 1000 {
        test_pass "1000 unique values encoded"
    }
    else {
        test_fail "1000 unique" "not 1000"
    }
}
else {
    test_fail "1000 unique" "rc=`=_rc'"
}

* Test 22.7: 5000 unique values
clear
set obs 5000
gen str15 x = "val_" + string(_n)
capture cencode x, generate(code)
if _rc == 0 {
    quietly sum code
    if r(max) == 5000 {
        test_pass "5000 unique values encoded"
    }
    else {
        test_fail "5000 unique" "max not 5000"
    }
}
else {
    test_fail "5000 unique" "rc=`=_rc'"
}

* Test 22.8: Label name same as variable name
clear
set obs 10
gen str10 myvar = "cat" + string(mod(_n, 3))
capture label drop myvar
capture cencode myvar, generate(myvar_code) label(myvar)
if _rc == 0 {
    test_pass "label name same as var name"
}
else {
    test_fail "label=var name" "rc=`=_rc'"
}
capture label drop myvar

/*******************************************************************************
 * SECTION 23: Large dataset tests
 ******************************************************************************/
print_section "Large Dataset Tests"

* Test 23.1: 10K observations, 10 unique
clear
set obs 10000
gen str10 x = "cat" + string(mod(_n, 10))
cencode x, generate(code)
quietly tab code
if r(r) == 10 {
    test_pass "10K obs, 10 unique"
}
else {
    test_fail "10K/10" "not 10 unique"
}

* Test 23.2: 50K observations, 50 unique
clear
set obs 50000
gen str10 x = "cat" + string(mod(_n, 50))
cencode x, generate(code)
quietly tab code
if r(r) == 50 {
    test_pass "50K obs, 50 unique"
}
else {
    test_fail "50K/50" "not 50 unique"
}

* Test 23.3: 100K observations, 100 unique
clear
set obs 100000
gen str15 x = "category_" + string(mod(_n, 100))
capture cencode x, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 100 {
        test_pass "100K obs, 100 unique"
    }
    else {
        test_fail "100K/100" "not 100 unique"
    }
}
else {
    test_fail "100K/100" "rc=`=_rc'"
}

* Test 23.4: 100K observations, 1K unique (high cardinality)
clear
set obs 100000
gen str15 x = "val_" + string(mod(_n, 1000))
capture cencode x, generate(code)
if _rc == 0 {
    quietly sum code
    if r(max) == 1000 {
        test_pass "100K obs, 1K unique (high cardinality)"
    }
    else {
        test_fail "100K/1K" "max not 1000"
    }
}
else {
    test_fail "100K/1K" "rc=`=_rc'"
}

* Test 23.5: 100K observations, 2 unique (low cardinality)
clear
set obs 100000
gen str10 x = cond(mod(_n, 2) == 0, "A", "B")
cencode x, generate(code)
quietly tab code
if r(r) == 2 {
    test_pass "100K obs, 2 unique (low cardinality)"
}
else {
    test_fail "100K/2" "not 2 unique"
}

* Test 23.6: 500K observations
clear
set obs 500000
gen str10 x = "cat" + string(mod(_n, 100))
capture cencode x, generate(code)
if _rc == 0 {
    test_pass "500K observations"
}
else {
    test_fail "500K obs" "rc=`=_rc'"
}

* Test 23.8: 1M observations
clear
set obs 1000000
gen str10 x = "cat" + string(mod(_n, 10))
capture cencode x, generate(code)
if _rc == 0 {
    quietly tab code
    if r(r) == 10 {
        test_pass "1M observations"
    }
    else {
        test_fail "1M obs" "wrong unique count"
    }
}
else {
    test_fail "1M obs" "rc=`=_rc'"
}

* Test 23.9: Large with if condition
clear
set obs 500000
gen str10 x = "cat" + string(mod(_n, 50))
gen byte flag = mod(_n, 3) == 0
cencode x if flag == 1, generate(code)
count if !missing(code)
if abs(r(N) - 166667) <= 1 {
    test_pass "large with if condition"
}
else {
    test_fail "large if" "wrong count"
}

* Test 23.10: Large with in range
clear
set obs 500000
gen str10 x = "cat" + string(mod(_n, 50))
cencode x in 1/100000, generate(code)
count if !missing(code)
if r(N) == 100000 {
    test_pass "large with in range"
}
else {
    test_fail "large in" "wrong count"
}

/*******************************************************************************
 * SECTION 24: Extended if/in conditions
 ******************************************************************************/
print_section "Extended if/in Conditions"

* Test 24.1: if on numeric condition
sysuse auto, clear
cencode make if price > 10000, generate(code)
count if !missing(code)
count if price > 10000
if r(N) > 0 {
    test_pass "if numeric condition"
}
else {
    test_fail "if numeric" "wrong count"
}
drop code

* Test 24.2: if on string condition
clear
set obs 100
gen str10 category = "cat" + string(mod(_n, 5))
gen str10 group = cond(_n <= 50, "GroupA", "GroupB")
cencode category if group == "GroupA", generate(code)
count if !missing(code) & group == "GroupA"
if r(N) == 50 {
    test_pass "if string condition"
}
else {
    test_fail "if string" "wrong count"
}

* Test 24.3: Complex if with AND
sysuse auto, clear
cencode make if price > 5000 & foreign == 1, generate(code)
count if !missing(code)
local n_enc = r(N)
count if price > 5000 & foreign == 1
if `n_enc' == r(N) {
    test_pass "if with AND condition"
}
else {
    test_fail "if AND" "wrong count"
}
drop code

* Test 24.4: Complex if with OR
sysuse auto, clear
cencode make if price > 10000 | mpg > 30, generate(code)
count if !missing(code)
local n_enc = r(N)
count if price > 10000 | mpg > 30
if `n_enc' == r(N) {
    test_pass "if with OR condition"
}
else {
    test_fail "if OR" "wrong count"
}
drop code

* Test 24.5: in first few rows
sysuse auto, clear
cencode make in 1/5, generate(code)
count if !missing(code) in 1/5
if r(N) == 5 {
    test_pass "in first 5 rows"
}
else {
    test_fail "in first 5" "wrong count"
}
drop code

* Test 24.6: in last few rows
sysuse auto, clear
local last5 = _N - 4
cencode make in `last5'/`=_N', generate(code)
count if !missing(code)
if r(N) == 5 {
    test_pass "in last 5 rows"
}
else {
    test_fail "in last 5" "wrong count"
}
drop code

* Test 24.7: in middle rows
sysuse auto, clear
cencode make in 30/40, generate(code)
count if !missing(code)
if r(N) == 11 {
    test_pass "in middle rows 30/40"
}
else {
    test_fail "in middle" "wrong count"
}
drop code

* Test 24.8: if that selects no rows
sysuse auto, clear
cencode make if price > 100000, generate(code)
count if !missing(code)
if r(N) == 0 {
    test_pass "if selects no rows"
}
else {
    test_fail "if no rows" "should be 0"
}
drop code

* Test 24.9: if that selects all rows
sysuse auto, clear
cencode make if price > 0, generate(code)
count if !missing(code)
if r(N) == _N {
    test_pass "if selects all rows"
}
else {
    test_fail "if all rows" "should be all"
}
drop code

* Test 24.10: Compare if/in with encode
sysuse auto, clear
benchmark_encode make, testname("vs encode: if price>8000") if2("price > 8000")

sysuse auto, clear
benchmark_encode make, testname("vs encode: in 10/60") in2("10/60")

/*******************************************************************************
 * SECTION 25: Replace option tests
 ******************************************************************************/
print_section "Replace Option"

* Test 25.1: Replace basic functionality
clear
set obs 10
gen str20 x = "category_" + string(mod(_n - 1, 5) + 1)
cencode x, replace
capture confirm numeric variable x
if _rc == 0 {
    test_pass "replace converts to numeric"
}
else {
    test_fail "replace basic" "variable not numeric"
}

* Test 25.2: Replace preserves values
clear
set obs 10
gen str10 x = "cat" + string(mod(_n - 1, 3) + 1)
clonevar x_backup = x
cencode x, replace
* Check that each unique string maps to same code
quietly tab x
if r(r) == 3 {
    test_pass "replace preserves unique values"
}
else {
    test_fail "replace values" "expected 3 unique"
}

* Test 25.3: Replace comparison with encode
clear
set obs 50
gen str20 x = "item_" + string(mod(_n - 1, 10) + 1)
clonevar x_c = x
clonevar x_s = x
cencode x_c, replace
encode x_s, generate(x_s_enc)
drop x_s
rename x_s_enc x_s
* Compare values
count if x_c != x_s
if r(N) == 0 {
    test_pass "replace matches encode"
}
else {
    test_fail "replace vs encode" "`=r(N)' values differ"
}

* Test 25.4: Replace with empty strings
clear
set obs 20
gen str20 x = cond(mod(_n, 3) == 0, "", "val" + string(_n))
clonevar x_backup = x
cencode x, replace
count if missing(x) & x_backup == ""
local n_miss = r(N)
count if !missing(x) & x_backup != ""
local n_enc = r(N)
* Empty strings should become missing (just like encode)
if `n_miss' > 0 {
    test_pass "replace handles empty strings"
}
else {
    test_fail "replace empty" "expected some missing"
}

* Test 25.5: Replace with if condition
clear
set obs 100
gen str10 x = "cat" + string(mod(_n - 1, 5) + 1)
gen flag = mod(_n, 2)
cencode x if flag == 1, replace
capture confirm numeric variable x
if _rc == 0 {
    count if !missing(x) & flag == 1
    local n_conv = r(N)
    count if flag == 1
    if `n_conv' == r(N) {
        test_pass "replace with if condition"
    }
    else {
        test_fail "replace if" "wrong count"
    }
}
else {
    test_fail "replace if" "not converted"
}

* Test 25.6: Replace with in range
clear
set obs 100
gen str10 x = "cat" + string(mod(_n - 1, 5) + 1)
cencode x in 1/50, replace
capture confirm numeric variable x
if _rc == 0 {
    count if !missing(x) in 1/50
    if r(N) == 50 {
        test_pass "replace with in range"
    }
    else {
        test_fail "replace in" "wrong count"
    }
}
else {
    test_fail "replace in" "not converted"
}

* Test 25.7: Replace with label option
clear
set obs 20
gen str10 x = "val" + string(mod(_n - 1, 4) + 1)
capture label drop my_custom_label
cencode x, replace label(my_custom_label)
local attached : value label x
if "`attached'" == "my_custom_label" {
    test_pass "replace with custom label"
}
else {
    test_fail "replace label" "label not my_custom_label"
}
capture label drop my_custom_label

/*******************************************************************************
 * SECTION 26: Multiple variables (varlist) tests
 ******************************************************************************/
print_section "Multiple Variables (varlist)"

* Test 26.1: Two variables with generate
clear
set obs 50
gen str10 a = "cat" + string(mod(_n - 1, 5) + 1)
gen str10 b = "grp" + string(mod(_n - 1, 3) + 1)
cencode a b, generate(a_num b_num)
capture confirm numeric variable a_num
local rc1 = _rc
capture confirm numeric variable b_num
local rc2 = _rc
if `rc1' == 0 & `rc2' == 0 {
    test_pass "two variables with generate"
}
else {
    test_fail "two vars gen" "not created"
}
drop a_num b_num

* Test 26.2: Two variables comparison with encode
clear
set obs 50
gen str10 a = "cat" + string(mod(_n - 1, 5) + 1)
gen str10 b = "grp" + string(mod(_n - 1, 3) + 1)
capture label drop a_c b_c a_s b_s
cencode a b, generate(a_c b_c)
encode a, generate(a_s)
encode b, generate(b_s)
count if a_c != a_s
local n1 = r(N)
count if b_c != b_s
local n2 = r(N)
if `n1' == 0 & `n2' == 0 {
    test_pass "two variables match encode"
}
else {
    test_fail "two vars encode" "a diff=`n1' b diff=`n2'"
}
drop a_c b_c a_s b_s

* Test 26.3: Three variables
clear
set obs 30
gen str10 x = "x_" + string(mod(_n - 1, 3) + 1)
gen str10 y = "y_" + string(mod(_n - 1, 4) + 1)
gen str10 z = "z_" + string(mod(_n - 1, 5) + 1)
capture label drop x_num y_num z_num x_e y_e z_e
cencode x y z, generate(x_num y_num z_num)
encode x, generate(x_e)
encode y, generate(y_e)
encode z, generate(z_e)
count if x_num != x_e
local n1 = r(N)
count if y_num != y_e
local n2 = r(N)
count if z_num != z_e
local n3 = r(N)
if `n1' == 0 & `n2' == 0 & `n3' == 0 {
    test_pass "three variables match encode"
}
else {
    test_fail "three vars" "x diff=`n1' y diff=`n2' z diff=`n3'"
}
drop x_num y_num z_num x_e y_e z_e

* Test 26.4: Multiple variables with replace
clear
set obs 40
gen str10 a = "cat" + string(mod(_n - 1, 4) + 1)
gen str10 b = "grp" + string(mod(_n - 1, 2) + 1)
clonevar a_backup = a
clonevar b_backup = b
encode a_backup, generate(a_expected)
encode b_backup, generate(b_expected)
cencode a b, replace
count if a != a_expected
local n1 = r(N)
count if b != b_expected
local n2 = r(N)
if `n1' == 0 & `n2' == 0 {
    test_pass "multiple variables replace"
}
else {
    test_fail "multi replace" "a diff=`n1' b diff=`n2'"
}

* Test 26.5: Multiple variables with if condition
clear
set obs 50
gen str10 a = "cat" + string(mod(_n - 1, 5) + 1)
gen str10 b = "grp" + string(mod(_n - 1, 3) + 1)
gen flag = mod(_n, 2)
capture label drop a_num b_num
cencode a b if flag == 1, generate(a_num b_num)
count if !missing(a_num) & flag == 1
local n_conv = r(N)
count if flag == 1
if `n_conv' == r(N) {
    test_pass "multiple vars with if"
}
else {
    test_fail "multi if" "wrong count"
}
drop a_num b_num

* Test 26.6: Multiple variables count mismatch error
clear
set obs 10
gen str10 a = "cat" + string(_n)
gen str10 b = "grp" + string(_n)
capture cencode a b, generate(only_one)
if _rc != 0 {
    test_pass "error: generate count mismatch"
}
else {
    test_fail "count mismatch" "should error"
}

/*******************************************************************************
 * SECTION 27: Empty string edge cases comprehensive
 ******************************************************************************/
print_section "Empty String Edge Cases Comprehensive"

* Test 27.1: Single empty string observation
clear
set obs 5
gen str10 x = "val" + string(_n)
replace x = "" in 3
cencode x, generate(x_num)
encode x, generate(x_enc)
count if missing(x_num) != missing(x_enc)
if r(N) == 0 {
    test_pass "single empty matches encode"
}
else {
    test_fail "single empty" "missing pattern differs"
}
drop x_num x_enc

* Test 27.2: First observation empty
clear
set obs 10
gen str10 x = "val" + string(_n)
replace x = "" in 1
cencode x, generate(x_num)
encode x, generate(x_enc)
count if missing(x_num) != missing(x_enc)
if r(N) == 0 {
    test_pass "first empty matches encode"
}
else {
    test_fail "first empty" "missing pattern differs"
}
drop x_num x_enc

* Test 27.3: Last observation empty
clear
set obs 10
gen str10 x = "val" + string(_n)
replace x = "" in 10
cencode x, generate(x_num)
encode x, generate(x_enc)
count if missing(x_num) != missing(x_enc)
if r(N) == 0 {
    test_pass "last empty matches encode"
}
else {
    test_fail "last empty" "missing pattern differs"
}
drop x_num x_enc

* Test 27.4: Consecutive empties
clear
set obs 20
gen str10 x = "val" + string(_n)
replace x = "" in 5/10
cencode x, generate(x_num)
encode x, generate(x_enc)
count if missing(x_num) != missing(x_enc)
if r(N) == 0 {
    test_pass "consecutive empties match encode"
}
else {
    test_fail "consecutive empty" "missing pattern differs"
}
drop x_num x_enc

* Test 27.5: Bookend pattern (first and last empty)
clear
set obs 20
gen str10 x = "val" + string(_n)
replace x = "" in 1
replace x = "" in 20
cencode x, generate(x_num)
encode x, generate(x_enc)
count if missing(x_num) != missing(x_enc)
if r(N) == 0 {
    test_pass "bookend empties match encode"
}
else {
    test_fail "bookend empty" "missing pattern differs"
}
drop x_num x_enc

/*******************************************************************************
 * SECTION 28: Combined options tests
 ******************************************************************************/
print_section "Combined Options"

* Test 28.1: verbose with threads
clear
set obs 1000
gen str20 x = "category_" + string(mod(_n - 1, 10) + 1)
capture cencode x, generate(x_num) verbose threads(2)
if _rc == 0 {
    quietly tab x_num
    if r(r) == 10 {
        test_pass "verbose with threads"
    }
    else {
        test_fail "verbose+threads" "wrong unique count"
    }
}
else {
    test_fail "verbose+threads" "rc=`=_rc'"
}
drop x_num

* Test 28.2: label with noextend (existing label)
clear
set obs 20
gen str10 x = "cat" + string(mod(_n - 1, 3) + 1)
capture label drop test_label
label define test_label 1 "cat1" 2 "cat2" 3 "cat3"
capture cencode x, generate(x_num) label(test_label) noextend
if _rc == 0 {
    local lbl_1 : label test_label 1
    local lbl_2 : label test_label 2
    if "`lbl_1'" == "cat1" & "`lbl_2'" == "cat2" {
        test_pass "label with noextend"
    }
    else {
        test_fail "label+noextend" "labels changed"
    }
}
else {
    test_fail "label+noextend" "rc=`=_rc'"
}
capture label drop test_label

* Test 28.3: if with threads
clear
set obs 1000
gen str10 x = "cat" + string(mod(_n - 1, 5) + 1)
gen flag = mod(_n, 2)
capture cencode x if flag == 1, generate(x_num) threads(2)
if _rc == 0 {
    count if !missing(x_num) & flag == 1
    local n_conv = r(N)
    count if flag == 1
    if `n_conv' == r(N) {
        test_pass "if with threads"
    }
    else {
        test_fail "if+threads" "wrong count"
    }
}
else {
    test_fail "if+threads" "rc=`=_rc'"
}
drop x_num

* Test 28.4: in with verbose
clear
set obs 100
gen str10 x = "cat" + string(mod(_n - 1, 5) + 1)
capture cencode x in 1/50, generate(x_num) verbose
if _rc == 0 {
    count if !missing(x_num) in 1/50
    if r(N) == 50 {
        test_pass "in with verbose"
    }
    else {
        test_fail "in+verbose" "wrong count"
    }
}
else {
    test_fail "in+verbose" "rc=`=_rc'"
}
drop x_num

* Test 28.5: replace with threads
clear
set obs 5000
gen str10 x = "cat" + string(mod(_n - 1, 20) + 1)
capture cencode x, replace threads(4)
if _rc == 0 {
    capture confirm numeric variable x
    if _rc == 0 {
        quietly tab x
        if r(r) == 20 {
            test_pass "replace with threads"
        }
        else {
            test_fail "replace+threads" "wrong unique count"
        }
    }
    else {
        test_fail "replace+threads" "not numeric"
    }
}
else {
    test_fail "replace+threads" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 29: strL variable type tests
 ******************************************************************************/
print_section "strL Variable Tests"

* Test 29.1: Basic strL encoding
clear
set obs 10
gen strL x = "value_" + string(_n)
capture cencode x, generate(x_num)
if _rc == 0 {
    quietly tab x_num
    if r(r) == 10 {
        test_pass "basic strL encoding"
    }
    else {
        test_fail "strL basic" "wrong unique count"
    }
}
else {
    test_fail "strL basic" "rc=`=_rc'"
}

* Test 29.2: strL with long content
clear
set obs 5
gen strL x = ""
replace x = "short" in 1
replace x = "a" * 500 + "_1" in 2
replace x = "b" * 500 + "_2" in 3
replace x = "c" * 500 + "_3" in 4
replace x = "medium length text" in 5
capture cencode x, generate(x_num)
if _rc == 0 {
    quietly tab x_num
    if r(r) == 5 {
        test_pass "strL with long content"
    }
    else {
        test_fail "strL long" "wrong unique count"
    }
}
else {
    test_fail "strL long" "rc=`=_rc'"
}

* Test 29.3: strL comparison with encode
clear
set obs 20
gen strL x = "category_" + string(mod(_n - 1, 5) + 1)
capture label drop x_c x_s
capture cencode x, generate(x_c)
local cencode_rc = _rc
capture encode x, generate(x_s)
local encode_rc = _rc
if `cencode_rc' == `encode_rc' {
    if `cencode_rc' == 0 {
        count if x_c != x_s
        if r(N) == 0 {
            test_pass "strL matches encode"
        }
        else {
            test_fail "strL encode" "`=r(N)' differ"
        }
    }
    else {
        test_pass "strL - both error as expected"
    }
}
else {
    test_fail "strL encode" "cencode rc=`cencode_rc' encode rc=`encode_rc'"
}

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cencode returns the same error codes as Stata's encode
 * when given invalid inputs or error conditions.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Variable doesn't exist
sysuse auto, clear
test_error_match, stata_cmd(encode nonexistent_var, generate(test)) ctools_cmd(cencode nonexistent_var, generate(test)) testname("nonexistent variable")

* Generate variable already exists
sysuse auto, clear
test_error_match, stata_cmd(encode make, generate(price)) ctools_cmd(cencode make, generate(price)) testname("generate var already exists")

* No generate option (missing required option)
sysuse auto, clear
test_error_match, stata_cmd(encode make) ctools_cmd(cencode make) testname("missing generate option")

* Numeric variable (encode expects string)
sysuse auto, clear
test_error_match, stata_cmd(encode price, generate(test)) ctools_cmd(cencode price, generate(test)) testname("numeric variable input")

* Empty dataset
clear
set obs 0
gen str10 x = ""
test_error_match, stata_cmd(encode x, generate(test)) ctools_cmd(cencode x, generate(test)) testname("empty dataset")

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
* End of cencode validation
noi print_summary "cencode"
}
