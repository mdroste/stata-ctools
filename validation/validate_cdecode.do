/*******************************************************************************
 * validate_cdecode.do
 *
 * Comprehensive validation tests for cdecode vs native Stata decode
 * Tests all options, data types, and edge cases
 *
 * cdecode: C-accelerated numeric to string decoding
 * Syntax: cdecode varname [if] [in], generate(newvar) [maxlength(#) verbose threads(#)]
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for cdecode..."

/*******************************************************************************
 * Helper: Compare cdecode vs decode
 ******************************************************************************/
capture program drop benchmark_decode
program define benchmark_decode
    syntax varname, testname(string) [GENerate(name) MAXLength(integer 0) if2(string) in2(string)]

    * Set default generate name
    if "`generate'" == "" local generate = "decoded"

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build options
    local opts "generate(`generate'_c)"
    if `maxlength' > 0 local opts "`opts' maxlength(`maxlength')"

    local opts_stata "generate(`generate'_s)"
    if `maxlength' > 0 local opts_stata "`opts_stata' maxlength(`maxlength')"

    * Run cdecode
    capture drop `generate'_c
    capture cdecode `varlist' `ifin', `opts'
    local rc_c = _rc

    * Run native decode
    capture drop `generate'_s
    capture decode `varlist' `ifin', `opts_stata'
    local rc_s = _rc

    * Check both succeeded or both failed
    if `rc_c' != `rc_s' {
        test_fail "`testname'" "cdecode rc=`rc_c', decode rc=`rc_s'"
        exit
    }

    if `rc_c' != 0 {
        test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare string values
    quietly count if `generate'_c != `generate'_s
    local ndiff = r(N)

    if `ndiff' > 0 {
        test_fail "`testname'" "`ndiff' values differ"
    }
    else {
        test_pass "`testname'"
    }

    * Cleanup
    capture drop `generate'_c `generate'_s
end

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
print_section "Plugin Check"

sysuse auto, clear
capture noisily cdecode foreign, generate(foreign_str)
if _rc != 0 {
    test_fail "cdecode plugin load" "plugin returned error `=_rc'"
    print_summary "cdecode"
    exit 1
}
test_pass "cdecode plugin loads and runs"
drop foreign_str

/*******************************************************************************
 * SECTION 2: Basic functionality (10 tests)
 ******************************************************************************/
print_section "Basic Functionality"

* Test 2.1: Basic decoding creates variable
sysuse auto, clear
capture drop foreign_str
cdecode foreign, generate(foreign_str)
capture confirm string variable foreign_str
if _rc == 0 {
    test_pass "basic decoding creates string variable"
}
else {
    test_fail "basic decoding" "variable not string"
}

* Test 2.2: Compare with decode
sysuse auto, clear
benchmark_decode foreign, testname("vs decode: foreign")

* Test 2.3: Labels correctly decoded
sysuse auto, clear
cdecode foreign, generate(foreign_str)
local pass = 1
forvalues i = 1/10 {
    local code = foreign[`i']
    local expected : label (foreign) `code'
    local actual = foreign_str[`i']
    if "`expected'" != "`actual'" {
        local pass = 0
    }
}
if `pass' {
    test_pass "labels correctly decoded"
}
else {
    test_fail "labels decoded" "mismatch found"
}
drop foreign_str

* Test 2.4: Variable with missing values
clear
set obs 20
gen byte x = _n
replace x = . in 1/5
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten" ///
    11 "eleven" 12 "twelve" 13 "thirteen" 14 "fourteen" 15 "fifteen"
label values x x_lbl
benchmark_decode x, testname("vs decode: var with missings")

* Test 2.5: verbose option - verify correctness with verbose enabled
sysuse auto, clear
benchmark_decode foreign, testname("[syntax] verbose option") cdecodeopts(verbose)

* Test 2.6: threads option - verify correctness with threads(2)
sysuse auto, clear
benchmark_decode foreign, testname("[syntax] threads(2) option") cdecodeopts(threads(2))

* Test 2.7: Single observation
clear
set obs 1
gen byte x = 1
label define x_lbl 1 "one"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "one" {
    test_pass "single observation"
}
else {
    test_fail "single observation" "got `=x_str[1]'"
}

* Test 2.8: Two observations
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "first" 2 "second"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "first" & x_str[2] == "second" {
    test_pass "two observations"
}
else {
    test_fail "two observations" "wrong values"
}

* Test 2.9: Many label values
clear
set obs 100
gen byte x = mod(_n - 1, 10) + 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "one" & x == 1
local n1 = r(N)
count if x_str == "ten" & x == 10
local n10 = r(N)
if `n1' == 10 & `n10' == 10 {
    test_pass "many label values"
}
else {
    test_fail "many labels" "wrong counts"
}

* Test 2.10: Compare many labels with decode
clear
set obs 50
gen byte category = mod(_n - 1, 5) + 1
label define cat_lbl 1 "Alpha" 2 "Beta" 3 "Gamma" 4 "Delta" 5 "Epsilon"
label values category cat_lbl
benchmark_decode category, testname("vs decode: 5 categories")

/*******************************************************************************
 * SECTION 3: Missing values (10 tests)
 ******************************************************************************/
print_section "Missing Values"

* Test 3.1: Missing numeric becomes empty string
clear
set obs 20
gen byte x = _n
replace x = . in 1/5
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten" ///
    11 "eleven" 12 "twelve" 13 "thirteen" 14 "fourteen" 15 "fifteen" ///
    16 "sixteen" 17 "seventeen" 18 "eighteen" 19 "nineteen" 20 "twenty"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "" & missing(x)
local n_empty = r(N)
count if missing(x)
if `n_empty' == r(N) {
    test_pass "missing numeric -> empty string"
}
else {
    test_fail "missing handling" "expected `=r(N)', got `n_empty'"
}

* Test 3.2: Non-missing values decoded
clear
set obs 20
gen byte x = _n
replace x = . in 1/5
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten" ///
    11 "eleven" 12 "twelve" 13 "thirteen" 14 "fourteen" 15 "fifteen" ///
    16 "sixteen" 17 "seventeen" 18 "eighteen" 19 "nineteen" 20 "twenty"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str != "" & !missing(x)
local n_decoded = r(N)
count if !missing(x)
if `n_decoded' == r(N) {
    test_pass "non-missing values decoded"
}
else {
    test_fail "non-missing" "expected `=r(N)', got `n_decoded'"
}

* Test 3.3: All missing values
clear
set obs 10
gen x = .
label define x_lbl 1 "one"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 10 {
    test_pass "all missing -> all empty"
}
else {
    test_fail "all missing" "not all empty"
}

* Test 3.4: First observation missing
clear
set obs 10
gen byte x = _n
replace x = . in 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten"
label values x x_lbl
cdecode x, generate(x_str)
local pass = (x_str[1] == "") & (x_str[2] == "two")
if `pass' {
    test_pass "first observation missing"
}
else {
    test_fail "first missing" "wrong handling"
}

* Test 3.5: Last observation missing
clear
set obs 10
gen byte x = _n
replace x = . in 10
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten"
label values x x_lbl
cdecode x, generate(x_str)
local pass = (x_str[10] == "") & (x_str[9] == "nine")
if `pass' {
    test_pass "last observation missing"
}
else {
    test_fail "last missing" "wrong handling"
}

* Test 3.6: Alternating missing/non-missing
clear
set obs 20
gen byte x = cond(mod(_n, 2) == 0, 1, .)
label define x_lbl 1 "value"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 10 {
    test_pass "alternating missing"
}
else {
    test_fail "alternating" "expected 10 empty"
}

* Test 3.7: Extended missing values (.a, .b, etc.)
clear
set obs 5
gen x = .
replace x = .a in 1
replace x = .b in 2
replace x = .c in 3
replace x = 1 in 4
replace x = 2 in 5
label define x_lbl 1 "one" 2 "two"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "" in 1/3
if r(N) == 3 {
    test_pass "extended missing values"
}
else {
    test_fail "extended missing" "not all empty"
}

* Test 3.8: Compare missing handling with decode
clear
set obs 30
gen byte x = mod(_n, 5) + 1
replace x = . if mod(_n, 6) == 0
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values x x_lbl
benchmark_decode x, testname("vs decode: missing values")

* Test 3.9: 90% missing
clear
set obs 100
gen byte x = cond(_n <= 10, mod(_n - 1, 5) + 1, .)
label define x_lbl 1 "a" 2 "b" 3 "c" 4 "d" 5 "e"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 90 {
    test_pass "90% missing"
}
else {
    test_fail "90% missing" "expected 90 empty"
}

* Test 3.10: Only one non-missing
clear
set obs 100
gen byte x = .
replace x = 1 in 50
label define x_lbl 1 "singleton"
label values x x_lbl
cdecode x, generate(x_str)
local pass = (x_str[50] == "singleton") & (x_str[1] == "")
if `pass' {
    test_pass "only one non-missing"
}
else {
    test_fail "singleton" "wrong handling"
}

/*******************************************************************************
 * SECTION 4: Unlabeled values (10 tests)
 ******************************************************************************/
print_section "Unlabeled Values"

* Test 4.1: Value without label -> empty string
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "one" 3 "three" 5 "five"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "" & (x == 2 | x == 4)
if r(N) == 2 {
    test_pass "unlabeled values -> empty"
}
else {
    test_fail "unlabeled" "not empty"
}

* Test 4.2: Labeled values decoded
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "one" 3 "three" 5 "five"
label values x x_lbl
cdecode x, generate(x_str)
local pass = (x_str[1] == "one") & (x_str[3] == "three") & (x_str[5] == "five")
if `pass' {
    test_pass "labeled values decoded"
}
else {
    test_fail "labeled values" "wrong decoding"
}

* Test 4.3: Compare unlabeled handling with decode
clear
set obs 10
gen byte x = _n
label define x_lbl 1 "one" 5 "five" 10 "ten"
label values x x_lbl
benchmark_decode x, testname("vs decode: sparse labels")

* Test 4.4: All values unlabeled
clear
set obs 5
gen byte x = _n
label define x_lbl 99 "not used"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 5 {
    test_pass "all unlabeled -> all empty"
}
else {
    test_fail "all unlabeled" "not all empty"
}

* Test 4.5: Negative values labeled
clear
set obs 5
gen int x = _n - 3  // -2, -1, 0, 1, 2
label define x_lbl -2 "minus two" -1 "minus one" 0 "zero" 1 "one" 2 "two"
label values x x_lbl
cdecode x, generate(x_str)
local pass = (x_str[1] == "minus two") & (x_str[3] == "zero") & (x_str[5] == "two")
if `pass' {
    test_pass "negative values labeled"
}
else {
    test_fail "negative labels" "wrong decoding"
}

* Test 4.6: Zero labeled
clear
set obs 3
gen byte x = _n - 1  // 0, 1, 2
label define x_lbl 0 "zero" 1 "one" 2 "two"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "zero" {
    test_pass "zero labeled"
}
else {
    test_fail "zero label" "got `=x_str[1]'"
}

* Test 4.7: Large label values
clear
set obs 3
gen int x = 1000 * _n
label define x_lbl 1000 "thousand" 2000 "two thousand" 3000 "three thousand"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "thousand" & x_str[3] == "three thousand" {
    test_pass "large label values"
}
else {
    test_fail "large values" "wrong decoding"
}

* Test 4.8: Mix of labeled, unlabeled, missing
clear
set obs 10
gen byte x = .
replace x = 1 in 1/3
replace x = 2 in 4/6
replace x = 3 in 7/8
label define x_lbl 1 "labeled_one" 3 "labeled_three"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "labeled_one"
local n1 = r(N)
count if x_str == "" & x == 2
local n2 = r(N)
count if x_str == "" & missing(x)
local nm = r(N)
if `n1' == 3 & `n2' == 3 & `nm' == 2 {
    test_pass "mix of labeled/unlabeled/missing"
}
else {
    test_fail "mixed" "wrong counts"
}

* Test 4.9: Compare mixed with decode
clear
set obs 20
gen byte x = mod(_n, 5)
replace x = . if mod(_n, 7) == 0
label define x_lbl 0 "zero" 2 "two" 4 "four"
label values x x_lbl
benchmark_decode x, testname("vs decode: mixed labels")

* Test 4.10: Float values (non-integer)
clear
set obs 3
gen float x = _n + 0.5
label define x_lbl 1 "one" 2 "two" 3 "three"
label values x x_lbl
cdecode x, generate(x_str)
* Non-integer values won't match integer labels
count if x_str == ""
if r(N) == 3 {
    test_pass "float values unlabeled"
}
else {
    test_fail "float values" "expected all empty"
}

/*******************************************************************************
 * SECTION 5: if/in conditions (10 tests)
 ******************************************************************************/
print_section "if/in Conditions"

* Test 5.1: if condition subset
sysuse auto, clear
cdecode foreign if price > 10000, generate(foreign_str)
count if foreign_str != "" & price > 10000
local n_dec = r(N)
count if foreign_str == "" & price <= 10000
local n_empty = r(N)
count if price > 10000
if `n_dec' == r(N) & `n_empty' > 0 {
    test_pass "if condition subset"
}
else {
    test_fail "if condition" "wrong pattern"
}
drop foreign_str

* Test 5.2: in range
sysuse auto, clear
cdecode foreign in 1/20, generate(foreign_str)
count if foreign_str != "" in 1/20
local in_range = r(N)
count if foreign_str == "" in 21/74
local out_range = r(N)
if `in_range' == 20 & `out_range' == 54 {
    test_pass "in range"
}
else {
    test_fail "in range" "wrong counts"
}
drop foreign_str

* Test 5.3: Empty if result
sysuse auto, clear
cdecode foreign if price > 100000, generate(foreign_str)
count if foreign_str != ""
if r(N) == 0 {
    test_pass "empty if result"
}
else {
    test_fail "empty if" "should have no decoded"
}
drop foreign_str

* Test 5.4: if matches all
sysuse auto, clear
cdecode foreign if price > 0, generate(foreign_str)
count if foreign_str != ""
if r(N) == _N {
    test_pass "if matches all"
}
else {
    test_fail "if all" "should decode all"
}
drop foreign_str

* Test 5.5: Combined if and in
sysuse auto, clear
* Use a condition that actually matches some observations
cdecode foreign if mpg > 20 in 1/50, generate(foreign_str)
* Compare with native decode
decode foreign if mpg > 20 in 1/50, generate(foreign_decode)
count if foreign_str != foreign_decode
if r(N) == 0 {
    test_pass "combined if and in"
}
else {
    test_fail "combined if/in" "`=r(N)' values differ"
}
drop foreign_str foreign_decode

* Test 5.6: Compare if with decode
sysuse auto, clear
benchmark_decode foreign, testname("vs decode: if price>8000") if2("price > 8000")

* Test 5.7: Compare in with decode
sysuse auto, clear
benchmark_decode foreign, testname("vs decode: in 1/30") in2("1/30")

* Test 5.8: Single observation in
sysuse auto, clear
cdecode foreign in 1/1, generate(foreign_str)
count if foreign_str != ""
if r(N) == 1 {
    test_pass "single observation in"
}
else {
    test_fail "single in" "wrong count"
}
drop foreign_str

* Test 5.9: if with condition on another variable
sysuse auto, clear
cdecode foreign if price < 10000, generate(foreign_str)
* Should decode foreign where price < 10000
count if foreign_str != "" & price < 10000
if r(N) > 0 {
    test_pass "if with condition on price"
}
else {
    test_fail "price condition if" "no decoded"
}
drop foreign_str

* Test 5.10: in from middle
sysuse auto, clear
cdecode foreign in 30/50, generate(foreign_str)
count if foreign_str != "" in 30/50
if r(N) == 21 {
    test_pass "in from middle"
}
else {
    test_fail "middle in" "wrong count"
}
drop foreign_str

/*******************************************************************************
 * SECTION 6: maxlength option (10 tests)
 ******************************************************************************/
print_section "maxlength Option"

* Test 6.1: maxlength truncates long labels
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "this is a very long label" 2 "short" 3 "another long one here"
label values x x_lbl
cdecode x, generate(x_str) maxlength(10)
if strlen(x_str[1]) <= 10 {
    test_pass "maxlength truncates"
}
else {
    test_fail "maxlength" "string too long: `=strlen(x_str[1])'"
}

* Test 6.2: Short labels unaffected by maxlength
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "abc" 2 "de" 3 "f"
label values x x_lbl
cdecode x, generate(x_str) maxlength(100)
if x_str[1] == "abc" & x_str[2] == "de" & x_str[3] == "f" {
    test_pass "short labels unaffected"
}
else {
    test_fail "short labels" "wrong values"
}

* Test 6.3: maxlength(1)
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "alpha" 2 "beta" 3 "gamma"
label values x x_lbl
cdecode x, generate(x_str) maxlength(1)
if strlen(x_str[1]) == 1 & strlen(x_str[2]) == 1 {
    test_pass "maxlength(1)"
}
else {
    test_fail "maxlength(1)" "not truncated to 1"
}

* Test 6.4: Compare maxlength with decode
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "very long label text here" 2 "another long one" ///
    3 "short" 4 "medium length" 5 "x"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(5)") maxlength(5)

* Test 6.5: maxlength with exact length match
clear
set obs 1
gen byte x = 1
label define x_lbl 1 "exact"
label values x x_lbl
cdecode x, generate(x_str) maxlength(5)
if x_str[1] == "exact" {
    test_pass "maxlength exact match"
}
else {
    test_fail "exact match" "got `=x_str[1]'"
}

* Test 6.6: maxlength preserves prefix
clear
set obs 1
gen byte x = 1
label define x_lbl 1 "abcdefghij"
label values x x_lbl
cdecode x, generate(x_str) maxlength(5)
if x_str[1] == "abcde" {
    test_pass "maxlength preserves prefix"
}
else {
    test_fail "prefix" "got `=x_str[1]'"
}

* Test 6.7: Default determines width from labels
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "short" 2 "this is longer" 3 "x"
label values x x_lbl
cdecode x, generate(x_str)
* Should accommodate longest label
if strlen(x_str[2]) == strlen("this is longer") {
    test_pass "default width from labels"
}
else {
    test_fail "default width" "wrong length"
}

* Test 6.8: maxlength larger than any label
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "a" 2 "bb" 3 "ccc"
label values x x_lbl
cdecode x, generate(x_str) maxlength(1000)
if x_str[1] == "a" & x_str[2] == "bb" & x_str[3] == "ccc" {
    test_pass "maxlength larger than labels"
}
else {
    test_fail "large maxlength" "wrong values"
}

* Test 6.9: Compare various maxlength with decode
clear
set obs 10
gen byte x = mod(_n, 3) + 1
label define x_lbl 1 "category_one" 2 "category_two" 3 "category_three"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(8)") maxlength(8)

* Test 6.10: maxlength with empty strings (missing/unlabeled)
clear
set obs 5
gen byte x = _n
replace x = . in 5
label define x_lbl 1 "label" 3 "another"
label values x x_lbl
cdecode x, generate(x_str) maxlength(3)
* Check empty strings for missing and unlabeled
if x_str[2] == "" & x_str[4] == "" & x_str[5] == "" {
    test_pass "maxlength with empties"
}
else {
    test_fail "maxlength empties" "wrong handling"
}

/*******************************************************************************
 * SECTION 7: Large datasets (10 tests)
 ******************************************************************************/
print_section "Large Datasets"

* Test 7.1: 100k observations
clear
set obs 100000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "Alpha" 2 "Beta" 3 "Gamma" 4 "Delta" 5 "Epsilon"
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == "Alpha"
    if r(N) == 20000 {
        test_pass "100k observations"
    }
    else {
        test_fail "100k obs" "wrong count"
    }
}
else {
    test_fail "100k obs" "rc=`=_rc'"
}

* Test 7.2: 500k observations - verify correctness on sample
clear
set obs 500000
gen byte category = mod(_n, 10) + 1
forvalues i = 1/10 {
    local lbl = "category_`i'"
    if `i' == 1 label define cat_lbl `i' "`lbl'"
    else label define cat_lbl `i' "`lbl'", add
}
label values category cat_lbl
benchmark_decode category, testname("500k observations")

* Test 7.3: 1M observations - verify correctness on sample
clear
set obs 1000000
gen byte category = mod(_n, 3) + 1
label define cat_lbl 1 "One" 2 "Two" 3 "Three"
label values category cat_lbl
benchmark_decode category, testname("1M observations")

* Test 7.4: Large with if condition
clear
set obs 500000
gen byte category = mod(_n, 5) + 1
gen byte flag = mod(_n, 2)
label define cat_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values category cat_lbl
cdecode category if flag == 1, generate(cat_str)
count if cat_str != ""
if r(N) == 250000 {
    test_pass "large with if condition"
}
else {
    test_fail "large if" "wrong count"
}

* Test 7.5: Large with long labels
clear
set obs 100000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "This is a fairly long label for category one" ///
    2 "Another somewhat lengthy label text here" ///
    3 "Short" ///
    4 "Medium length label" ///
    5 "The fifth and final category with long text"
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    test_pass "large with long labels"
}
else {
    test_fail "long labels" "rc=`=_rc'"
}

* Test 7.6: Large with threads(4)
clear
set obs 500000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values category cat_lbl
capture cdecode category, generate(cat_str) threads(4)
if _rc == 0 {
    test_pass "large with threads(4)"
}
else {
    test_fail "threads(4)" "rc=`=_rc'"
}

* Test 7.7: Large with many missing
clear
set obs 500000
gen byte category = cond(mod(_n, 10) == 0, mod(_n, 5) + 1, .)
label define cat_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == ""
    if r(N) == 450000 {
        test_pass "large with many missing"
    }
    else {
        test_fail "many missing" "wrong empty count"
    }
}
else {
    test_fail "many missing" "rc=`=_rc'"
}

* Test 7.8: Large dataset, small sample
clear
set obs 1000000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values category cat_lbl
cdecode category in 1/100, generate(cat_str)
count if cat_str != ""
if r(N) == 100 {
    test_pass "large dataset, small sample"
}
else {
    test_fail "small sample" "wrong count"
}

* Test 7.9: 50 unique labels
clear
set obs 100000
gen byte category = mod(_n, 50) + 1
forvalues i = 1/50 {
    local lbl = "label_`i'"
    if `i' == 1 label define cat_lbl `i' "`lbl'"
    else label define cat_lbl `i' "`lbl'", add
}
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == "label_1"
    if r(N) == 2000 {
        test_pass "50 unique labels"
    }
    else {
        test_fail "50 labels" "wrong count"
    }
}
else {
    test_fail "50 labels" "rc=`=_rc'"
}

* Test 7.10: Compare large dataset with decode
clear
set obs 50000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values category cat_lbl
benchmark_decode category, testname("vs decode: 50k obs")

/*******************************************************************************
 * SECTION 8: Special characters in labels (10 tests)
 ******************************************************************************/
print_section "Special Characters in Labels"

* Test 8.1: Spaces in labels
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "hello world" 2 "foo bar baz" 3 "a b c"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "hello world" {
    test_pass "spaces in labels"
}
else {
    test_fail "spaces" "got `=x_str[1]'"
}

* Test 8.2: Leading/trailing spaces
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "  leading" 2 "trailing  " 3 "  both  "
label values x x_lbl
cdecode x, generate(x_str)
* Note: Stata's decode may trim spaces
capture benchmark_decode x, testname("vs decode: leading/trailing spaces")

* Test 8.3: Commas in labels
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "one, two, three" 2 "a,b"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], ",") > 0 {
    test_pass "commas in labels"
}
else {
    test_fail "commas" "comma not preserved"
}

* Test 8.4: Parentheses and brackets
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "(parens)" 2 "[brackets]" 3 "{braces}"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "(parens)" & x_str[2] == "[brackets]" {
    test_pass "parentheses and brackets"
}
else {
    test_fail "brackets" "not preserved"
}

* Test 8.5: Ampersand and percent
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "A & B" 2 "50%"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], "&") > 0 & strpos(x_str[2], "%") > 0 {
    test_pass "ampersand and percent"
}
else {
    test_fail "amp/percent" "not preserved"
}

* Test 8.6: Apostrophes
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "don't" 2 "it's"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], "'") > 0 {
    test_pass "apostrophes"
}
else {
    test_fail "apostrophes" "not preserved"
}

* Test 8.7: Slashes
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "path/to/file" 2 "a\b"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], "/") > 0 {
    test_pass "forward slashes"
}
else {
    test_fail "slashes" "not preserved"
}

* Test 8.8: Compare special chars with decode
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "Hello!" 2 "#hashtag" 3 "$100" 4 "a+b=c" 5 "x*y"
label values x x_lbl
benchmark_decode x, testname("vs decode: special chars")

* Test 8.9: Dashes and underscores
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "with-dash" 2 "with_underscore" 3 "both-and_here"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "with-dash" & x_str[2] == "with_underscore" {
    test_pass "dashes and underscores"
}
else {
    test_fail "dash/underscore" "not preserved"
}

* Test 8.10: Mixed special characters
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "Item #1 (new)" 2 "50% off - sale!" 3 "A/B test: v2.0"
label values x x_lbl
cdecode x, generate(x_str)
benchmark_decode x, testname("vs decode: mixed special")

/*******************************************************************************
 * SECTION 9: Real-world datasets (10 tests)
 ******************************************************************************/
print_section "Real-World Datasets"

* Test 9.1: auto - foreign
sysuse auto, clear
benchmark_decode foreign, testname("auto: foreign")

* Test 9.2: citytemp - division
sysuse citytemp, clear
benchmark_decode division, testname("citytemp: division")

* Test 9.3: census - region
sysuse census, clear
benchmark_decode region, testname("census: region")

* Test 9.4: nlswork - race
webuse nlswork, clear
benchmark_decode race, testname("nlswork: race")

* Test 9.5: nlswork - msp
webuse nlswork, clear
benchmark_decode msp, testname("nlswork: msp")

* Test 9.6: citytemp - region
sysuse citytemp, clear
benchmark_decode region, testname("citytemp: region")

* Test 9.7: pop2000 - agegrp
sysuse pop2000, clear
benchmark_decode agegrp, testname("pop2000: agegrp")

* Test 9.8: voter - candidat
sysuse voter, clear
benchmark_decode candidat, testname("voter: candidat")

* Test 9.9: voter - inc
sysuse voter, clear
benchmark_decode inc, testname("voter: inc")

* Test 9.10: bpwide - agegrp
sysuse bpwide, clear
benchmark_decode agegrp, testname("bpwide: agegrp")

/*******************************************************************************
 * SECTION 10: Pathological data (10 tests)
 ******************************************************************************/
print_section "Pathological Data"

* Test 10.1: All same value
clear
set obs 1000
gen byte x = 1
label define x_lbl 1 "constant"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "constant"
if r(N) == 1000 {
    test_pass "all same value"
}
else {
    test_fail "all same" "wrong count"
}

* Test 10.2: All unique values (up to label limit)
clear
set obs 100
gen byte x = _n
forvalues i = 1/100 {
    local lbl = "val_`i'"
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
cdecode x, generate(x_str)
count if x_str != ""
if r(N) == 100 {
    test_pass "100 unique values"
}
else {
    test_fail "100 unique" "wrong count"
}

* Test 10.3: Binary values
clear
set obs 1000
gen byte x = mod(_n, 2)
label define x_lbl 0 "No" 1 "Yes"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "No"
local n_no = r(N)
count if x_str == "Yes"
local n_yes = r(N)
if `n_no' == 500 & `n_yes' == 500 {
    test_pass "binary values"
}
else {
    test_fail "binary" "wrong counts"
}

* Test 10.4: Very long label text
clear
set obs 3
gen byte x = _n
local longtext = "a" * 200
label define x_lbl 1 "`longtext'" 2 "short" 3 "medium length"
label values x x_lbl
capture cdecode x, generate(x_str)
if _rc == 0 {
    if strlen(x_str[1]) >= 100 {
        test_pass "very long label text"
    }
    else {
        test_fail "long text" "truncated too much"
    }
}
else {
    test_fail "long text" "rc=`=_rc'"
}

* Test 10.5: Bookend pattern
clear
set obs 100
gen byte x = cond(_n == 1 | _n == 100, 1, 2)
label define x_lbl 1 "bookend" 2 "middle"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "bookend" & x_str[100] == "bookend" & x_str[50] == "middle" {
    test_pass "bookend pattern"
}
else {
    test_fail "bookend" "wrong values"
}

* Test 10.6: Cyclic pattern
clear
set obs 1000
gen byte x = mod(_n - 1, 7) + 1
forvalues i = 1/7 {
    local lbl = "cycle_`i'"
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "cycle_1"
* 1000/7 = 142 or 143
if r(N) >= 142 & r(N) <= 143 {
    test_pass "cyclic pattern"
}
else {
    test_fail "cyclic" "wrong count"
}

* Test 10.7: Sparse data (mostly missing)
clear
set obs 10000
gen byte x = cond(mod(_n, 100) == 0, 1, .)
label define x_lbl 1 "rare"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "rare"
if r(N) == 100 {
    test_pass "sparse data"
}
else {
    test_fail "sparse" "wrong count"
}

* Test 10.8: Shuffled data
clear
set obs 100
set seed 12345
gen byte x = mod(_n, 5) + 1
gen double shuffle = runiform()
sort shuffle
drop shuffle
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "A"
if r(N) == 20 {
    test_pass "shuffled data"
}
else {
    test_fail "shuffled" "wrong count"
}

* Test 10.9: Compare pathological with decode
clear
set obs 50
gen byte x = cond(mod(_n, 10) == 0, 1, cond(mod(_n, 10) == 5, ., 2))
label define x_lbl 1 "rare" 2 "common"
label values x x_lbl
benchmark_decode x, testname("vs decode: pathological")

* Test 10.10: Extreme value range
clear
set obs 5
gen int x = 1
replace x = -1000 in 1
replace x = -1 in 2
replace x = 0 in 3
replace x = 1 in 4
replace x = 1000 in 5
label define x_lbl -1000 "neg thousand" -1 "neg one" 0 "zero" 1 "one" 1000 "thousand"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "neg thousand" & x_str[3] == "zero" & x_str[5] == "thousand" {
    test_pass "extreme value range"
}
else {
    test_fail "extreme range" "wrong values"
}

/*******************************************************************************
 * SECTION 11: Error handling (5 tests)
 ******************************************************************************/
print_section "Error Handling"

* Test 11.1: Generate var exists
sysuse auto, clear
gen foreign_str = "test"
capture cdecode foreign, generate(foreign_str)
if _rc == 110 {
    test_pass "error: generate var exists"
}
else {
    test_fail "var exists" "expected rc=110, got `=_rc'"
}

* Test 11.2: Source is string (not numeric)
sysuse auto, clear
capture cdecode make, generate(make_decoded)
if _rc != 0 {
    test_pass "error: string source"
}
else {
    test_fail "string source" "should error"
}

* Test 11.3: Source has no value label
clear
set obs 10
gen x = _n
capture cdecode x, generate(x_str)
if _rc != 0 {
    test_pass "error: no value label"
}
else {
    test_fail "no label" "should error"
}

* Test 11.4: Source doesn't exist
sysuse auto, clear
capture cdecode nonexistent, generate(ne_str)
if _rc != 0 {
    test_pass "error: nonexistent source"
}
else {
    test_fail "nonexistent" "should error"
}

* Test 11.5: Empty dataset - compare with Stata's decode
clear
set obs 0
gen byte x = .
label define x_lbl 1 "one"
label values x x_lbl
capture decode x, generate(x_dec)
local stata_rc = _rc
capture cdecode x, generate(x_str)
local cdecode_rc = _rc
if `stata_rc' == `cdecode_rc' {
    test_pass "empty dataset - matches Stata behavior"
}
else {
    test_fail "empty dataset" "cdecode rc=`cdecode_rc' but decode rc=`stata_rc'"
}

/*******************************************************************************
 * SECTION 12: Consistency/round-trip (5 tests)
 ******************************************************************************/
print_section "Consistency/Round-trip"

* Test 12.1: Round-trip encode then decode
sysuse auto, clear
encode make, generate(make_num)
decode make_num, generate(make_back)
local match = 1
forvalues i = 1/10 {
    if make[`i'] != make_back[`i'] local match = 0
}
if `match' {
    test_pass "round-trip encode/decode"
}
else {
    test_fail "round-trip" "mismatch"
}

* Test 12.2: Round-trip cencode then cdecode
sysuse auto, clear
cencode make, generate(make_num)
cdecode make_num, generate(make_back)
local match = 1
forvalues i = 1/10 {
    if make[`i'] != make_back[`i'] local match = 0
}
if `match' {
    test_pass "round-trip cencode/cdecode"
}
else {
    test_fail "cencode/cdecode" "mismatch"
}

* Test 12.3: cdecode matches decode
sysuse auto, clear
decode foreign, generate(foreign_d)
cdecode foreign, generate(foreign_c)
count if foreign_d != foreign_c
if r(N) == 0 {
    test_pass "cdecode matches decode"
}
else {
    test_fail "match decode" "`=r(N)' differ"
}

* Test 12.4: Observation count preserved
sysuse auto, clear
local orig_n = _N
cdecode foreign, generate(foreign_str)
if _N == `orig_n' {
    test_pass "observation count preserved"
}
else {
    test_fail "obs count" "changed"
}

* Test 12.5: Full comparison on census
sysuse census, clear
benchmark_decode region, testname("full comparison: census region")

/*******************************************************************************
 * SECTION 13: Comprehensive sysuse/webuse dataset tests
 ******************************************************************************/
print_section "Comprehensive Dataset Tests (sysuse/webuse)"

* auto dataset - foreign and rep78
sysuse auto, clear
benchmark_decode foreign, testname("auto: foreign (2 labels)")
benchmark_decode rep78, testname("auto: rep78 (with missing)")

* census dataset - region
sysuse census, clear
benchmark_decode region, testname("census: region (4 labels)")

* citytemp dataset - division and region
sysuse citytemp, clear
benchmark_decode division, testname("citytemp: division (9 labels)")
benchmark_decode region, testname("citytemp: region (4 labels)")

* pop2000 dataset - agegrp
sysuse pop2000, clear
benchmark_decode agegrp, testname("pop2000: agegrp")

* voter dataset - candidat, inc
sysuse voter, clear
benchmark_decode candidat, testname("voter: candidat")
benchmark_decode inc, testname("voter: inc")

* bpwide dataset - agegrp
sysuse bpwide, clear
benchmark_decode agegrp, testname("bpwide: agegrp")

* educ99gdp dataset
sysuse educ99gdp, clear
capture confirm variable region
if _rc == 0 {
    benchmark_decode region, testname("educ99gdp: region")
}

* nlsw88 dataset - comprehensive labeled variables
webuse nlsw88, clear
benchmark_decode race, testname("nlsw88: race")
benchmark_decode occupation, testname("nlsw88: occupation")
benchmark_decode industry, testname("nlsw88: industry")
benchmark_decode married, testname("nlsw88: married")
benchmark_decode collgrad, testname("nlsw88: collgrad")
benchmark_decode union, testname("nlsw88: union")
benchmark_decode south, testname("nlsw88: south")
benchmark_decode smsa, testname("nlsw88: smsa")
benchmark_decode c_city, testname("nlsw88: c_city")

* nlswork dataset - panel data with labels
webuse nlswork, clear
benchmark_decode race, testname("nlswork: race")
benchmark_decode msp, testname("nlswork: msp")
benchmark_decode nev_mar, testname("nlswork: nev_mar")
benchmark_decode collgrad, testname("nlswork: collgrad")
benchmark_decode union, testname("nlswork: union")
benchmark_decode south, testname("nlswork: south")

* grunfeld dataset
webuse grunfeld, clear
capture confirm numeric variable company
if _rc == 0 {
    capture local lbl : value label company
    if "`lbl'" != "" {
        benchmark_decode company, testname("grunfeld: company")
    }
}

* lifeexp dataset
webuse lifeexp, clear
capture confirm numeric variable region
if _rc == 0 {
    capture local lbl : value label region
    if "`lbl'" != "" {
        benchmark_decode region, testname("lifeexp: region")
    }
}

/*******************************************************************************
 * SECTION 14: Pathological label edge cases - Special characters
 ******************************************************************************/
print_section "Pathological Labels - Special Characters"

* Labels with embedded quotes
clear
set obs 3
gen byte x = _n
label define x_lbl 1 `"He said "hello""' 2 `"She's here"' 3 `"It's "fine""'
label values x x_lbl
capture benchmark_decode x, testname("labels with quotes")

* Labels with commas
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "Smith, John" 2 "Doe, Jane" 3 "one, two, three"
label values x x_lbl
benchmark_decode x, testname("labels with commas")

* Labels with tabs
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "col1	col2" 2 "tab	here" 3 "no tab"
label values x x_lbl
capture benchmark_decode x, testname("labels with tabs")

* Labels with semicolons and colons
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "time: 12:30" 2 "list; item; another" 3 "a:b;c"
label values x x_lbl
benchmark_decode x, testname("labels with semicolons/colons")

* Labels with special symbols
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "100% done" 2 "A & B" 3 "$500" 4 "#hashtag" 5 "@mention"
label values x x_lbl
benchmark_decode x, testname("labels with special symbols")

* Labels with parentheses and brackets
clear
set obs 4
gen byte x = _n
label define x_lbl 1 "(parentheses)" 2 "[brackets]" 3 "{braces}" 4 "<angles>"
label values x x_lbl
benchmark_decode x, testname("labels with brackets")

* Labels with slashes
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "path/to/file" 2 "a\b\c" 3 "yes/no"
label values x x_lbl
benchmark_decode x, testname("labels with slashes")

* Labels with math operators
clear
set obs 4
gen byte x = _n
label define x_lbl 1 "a+b=c" 2 "x-y" 3 "a*b" 4 "x/y"
label values x x_lbl
benchmark_decode x, testname("labels with math operators")

/*******************************************************************************
 * SECTION 15: Pathological labels - Leading/trailing spaces
 ******************************************************************************/
print_section "Pathological Labels - Whitespace"

* Labels with leading spaces
clear
set obs 3
gen byte x = _n
label define x_lbl 1 " leading" 2 "  double leading" 3 "   triple leading"
label values x x_lbl
benchmark_decode x, testname("labels with leading spaces")

* Labels with trailing spaces
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "trailing " 2 "trailing  " 3 "trailing   "
label values x x_lbl
benchmark_decode x, testname("labels with trailing spaces")

* Labels with both leading and trailing spaces
clear
set obs 3
gen byte x = _n
label define x_lbl 1 " both " 2 "  both  " 3 "   both   "
label values x x_lbl
benchmark_decode x, testname("labels with both lead/trail spaces")

* Labels with multiple internal spaces
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "hello  world" 2 "a   b   c" 3 "lots    of    space"
label values x x_lbl
benchmark_decode x, testname("labels with internal spaces")

* Labels that are only whitespace
clear
set obs 3
gen byte x = _n
label define x_lbl 1 " " 2 "  " 3 "   "
label values x x_lbl
capture benchmark_decode x, testname("labels with only whitespace")

/*******************************************************************************
 * SECTION 16: Pathological labels - Long labels
 ******************************************************************************/
print_section "Pathological Labels - Long Labels"

* Labels near 80 character limit
clear
set obs 3
gen byte x = _n
local lbl75 = "a" * 75
local lbl79 = "b" * 79
local lbl80 = "c" * 80
label define x_lbl 1 "`lbl75'" 2 "`lbl79'" 3 "`lbl80'"
label values x x_lbl
cdecode x, generate(x_str)
if strlen(x_str[3]) == 80 {
    test_pass "labels at 80 char limit"
}
else {
    test_fail "80 char labels" "wrong length: `=strlen(x_str[3])'"
}

* Mix of short and long labels
clear
set obs 4
gen byte x = _n
local long70 = "x" * 70
label define x_lbl 1 "a" 2 "`long70'" 3 "bc" 4 "longer label here"
label values x x_lbl
benchmark_decode x, testname("mixed short and long labels")

* Very long labels with special characters
clear
set obs 2
gen byte x = _n
local mixed_long = "a,b c:d;e" * 7
label define x_lbl 1 "`mixed_long'" 2 "short"
label values x x_lbl
capture benchmark_decode x, testname("long labels with special chars")

/*******************************************************************************
 * SECTION 17: Pathological labels - Empty and unusual labels
 ******************************************************************************/
print_section "Pathological Labels - Empty and Unusual"

* Empty label text
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "" 2 "not empty" 3 ""
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "" & x_str[2] == "not empty" & x_str[3] == "" {
    test_pass "empty label text"
}
else {
    test_fail "empty labels" "unexpected values"
}

* Single character labels
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "a" 2 "b" 3 "1" 4 " " 5 "."
label values x x_lbl
benchmark_decode x, testname("single character labels")

* Numeric-looking labels
clear
set obs 4
gen byte x = _n
label define x_lbl 1 "123" 2 "45.67" 3 "-89" 4 "1e10"
label values x x_lbl
benchmark_decode x, testname("numeric-looking labels")

/*******************************************************************************
 * SECTION 18: Pathological labels - Unicode/international characters
 ******************************************************************************/
print_section "Pathological Labels - Unicode/International"

* Basic Latin extended characters
clear
set obs 4
gen byte x = _n
label define x_lbl 1 "Jose Garcia" 2 "Francois Muller" 3 "Soren Jensen" 4 "Cafe Creme"
label values x x_lbl
benchmark_decode x, testname("basic Latin extended")

* German characters
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "Munchen" 2 "Dusseldorf" 3 "Strasse"
label values x x_lbl
capture benchmark_decode x, testname("German characters")

* Mixed ASCII and extended
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "Mix: A-Z and cafe" 2 "123 strasse" 3 "Jose #1"
label values x x_lbl
capture benchmark_decode x, testname("mixed ASCII and extended")

/*******************************************************************************
 * SECTION 19: Pathological labels - Gap values (non-consecutive)
 ******************************************************************************/
print_section "Pathological Labels - Non-Consecutive Values"

* Simple gaps
clear
set obs 5
gen byte x = 1 in 1
replace x = 2 in 2
replace x = 5 in 3
replace x = 10 in 4
replace x = 20 in 5
label define x_lbl 1 "one" 2 "two" 5 "five" 10 "ten" 20 "twenty"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "one" & x_str[3] == "five" & x_str[5] == "twenty" {
    test_pass "labels with gaps (non-consecutive values)"
}
else {
    test_fail "gap values" "wrong decoding"
}

* Large gaps
clear
set obs 4
gen int x = 1 in 1
replace x = 100 in 2
replace x = 1000 in 3
replace x = 10000 in 4
label define x_lbl 1 "tiny" 100 "hundred" 1000 "thousand" 10000 "ten thousand"
label values x x_lbl
benchmark_decode x, testname("large gaps in values")

* Random-looking value assignments
clear
set obs 5
gen int x = 7 in 1
replace x = 23 in 2
replace x = 42 in 3
replace x = 99 in 4
replace x = 137 in 5
label define x_lbl 7 "lucky" 23 "jordan" 42 "answer" 99 "gretzky" 137 "fine structure"
label values x x_lbl
benchmark_decode x, testname("random value assignments")

/*******************************************************************************
 * SECTION 20: Pathological labels - Negative values
 ******************************************************************************/
print_section "Pathological Labels - Negative Values"

* Simple negative values
clear
set obs 5
gen int x = -2 in 1
replace x = -1 in 2
replace x = 0 in 3
replace x = 1 in 4
replace x = 2 in 5
label define x_lbl -2 "minus two" -1 "minus one" 0 "zero" 1 "one" 2 "two"
label values x x_lbl
benchmark_decode x, testname("negative values labeled")

* Large negative values
clear
set obs 4
gen int x = -1000 in 1
replace x = -100 in 2
replace x = -10 in 3
replace x = -1 in 4
label define x_lbl -1000 "neg thousand" -100 "neg hundred" -10 "neg ten" -1 "neg one"
label values x x_lbl
benchmark_decode x, testname("large negative values")

* Mix of negative and positive
clear
set obs 6
gen int x = -999 in 1
replace x = -1 in 2
replace x = 0 in 3
replace x = 1 in 4
replace x = 999 in 5
replace x = . in 6
label define x_lbl -999 "min" -1 "neg" 0 "zero" 1 "pos" 999 "max"
label values x x_lbl
benchmark_decode x, testname("mixed negative/positive with missing")

* Only negative values
clear
set obs 5
gen int x = -1 * _n
label define x_lbl -1 "one" -2 "two" -3 "three" -4 "four" -5 "five"
label values x x_lbl
benchmark_decode x, testname("only negative values")

/*******************************************************************************
 * SECTION 21: Pathological labels - Very large number of labels
 ******************************************************************************/
print_section "Pathological Labels - Many Label Values"

* 100 unique label values
clear
set obs 100
gen byte x = _n
forvalues i = 1/100 {
    local lbl = "label_`i'"
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "label_1"
local n1 = r(N)
count if x_str == "label_100"
local n100 = r(N)
if `n1' == 1 & `n100' == 1 {
    test_pass "100 unique labels"
}
else {
    test_fail "100 labels" "wrong counts"
}

* 500 unique label values (using int to hold values)
clear
set obs 500
gen int x = _n
forvalues i = 1/500 {
    local lbl = "cat_`i'"
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
capture cdecode x, generate(x_str)
if _rc == 0 {
    count if x_str == "cat_1"
    local n1 = r(N)
    count if x_str == "cat_500"
    local n500 = r(N)
    if `n1' == 1 & `n500' == 1 {
        test_pass "500 unique labels"
    }
    else {
        test_fail "500 labels" "wrong counts"
    }
}
else {
    test_fail "500 labels" "rc=`=_rc'"
}

* Compare 200 labels with decode
clear
set obs 200
gen int x = _n
forvalues i = 1/200 {
    local lbl = "item_`i'"
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
benchmark_decode x, testname("vs decode: 200 unique labels")

/*******************************************************************************
 * SECTION 22: Missing value handling - Extended missing
 ******************************************************************************/
print_section "Missing Values - Extended"

* All extended missing values (.a through .z)
clear
set obs 28
gen x = .
replace x = .a in 1
replace x = .b in 2
replace x = .c in 3
replace x = .d in 4
replace x = .e in 5
replace x = .f in 6
replace x = .g in 7
replace x = .h in 8
replace x = .i in 9
replace x = .j in 10
replace x = .k in 11
replace x = .l in 12
replace x = .m in 13
replace x = .n in 14
replace x = .o in 15
replace x = .p in 16
replace x = .q in 17
replace x = .r in 18
replace x = .s in 19
replace x = .t in 20
replace x = .u in 21
replace x = .v in 22
replace x = .w in 23
replace x = .x in 24
replace x = .y in 25
replace x = .z in 26
replace x = 1 in 27
replace x = 2 in 28
label define x_lbl 1 "one" 2 "two"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "" in 1/26
local n_empty = r(N)
count if x_str != "" in 27/28
local n_filled = r(N)
if `n_empty' == 26 & `n_filled' == 2 {
    test_pass "all extended missing values (.a-.z)"
}
else {
    test_fail "extended missing" "expected 26 empty, 2 filled"
}

* Mixed extended missing
clear
set obs 10
gen x = _n
replace x = .a in 2
replace x = .z in 5
replace x = . in 8
label define x_lbl 1 "one" 3 "three" 4 "four" 6 "six" 7 "seven" 9 "nine" 10 "ten"
label values x x_lbl
benchmark_decode x, testname("vs decode: mixed extended missing")

* Alternating missing pattern
clear
set obs 20
gen x = .
forvalues i = 1(2)20 {
    replace x = `i' in `i'
}
forvalues i = 1/10 {
    local lbl = "val_`=`i'*2-1'"
    local val = `i' * 2 - 1
    if `i' == 1 label define x_lbl `val' "`lbl'"
    else label define x_lbl `val' "`lbl'", add
}
label values x x_lbl
benchmark_decode x, testname("alternating missing pattern")

* All missing values (stress test)
clear
set obs 1000
gen x = .
label define x_lbl 1 "one"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 1000 {
    test_pass "all missing (1000 obs)"
}
else {
    test_fail "all missing 1000" "not all empty"
}

/*******************************************************************************
 * SECTION 23: Numeric type coverage
 ******************************************************************************/
print_section "Numeric Type Coverage"

* byte variable with labels
clear
set obs 10
gen byte x = mod(_n - 1, 5) + 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
benchmark_decode x, testname("byte variable")

* int variable with labels
clear
set obs 10
gen int x = mod(_n - 1, 5) + 1000
label define x_lbl 1000 "thousand" 1001 "thousand one" 1002 "thousand two" 1003 "thousand three" 1004 "thousand four"
label values x x_lbl
benchmark_decode x, testname("int variable")

* long variable with labels
clear
set obs 10
gen long x = mod(_n - 1, 5) + 100000
label define x_lbl 100000 "100k" 100001 "100k+1" 100002 "100k+2" 100003 "100k+3" 100004 "100k+4"
label values x x_lbl
benchmark_decode x, testname("long variable")

* float variable (unusual but allowed)
clear
set obs 5
gen float x = _n
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
benchmark_decode x, testname("float variable (integer values)")

* double variable (unusual but allowed)
clear
set obs 5
gen double x = _n
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
benchmark_decode x, testname("double variable (integer values)")

* byte at boundaries
clear
set obs 3
gen byte x = -127 in 1
replace x = 0 in 2
replace x = 100 in 3
label define x_lbl -127 "min byte" 0 "zero" 100 "hundred"
label values x x_lbl
benchmark_decode x, testname("byte at boundaries")

* int at boundaries
clear
set obs 3
gen int x = -32767 in 1
replace x = 0 in 2
replace x = 32740 in 3
label define x_lbl -32767 "min int" 0 "zero" 32740 "near max int"
label values x x_lbl
benchmark_decode x, testname("int at boundaries")

/*******************************************************************************
 * SECTION 24: maxlength option comprehensive tests
 ******************************************************************************/
print_section "maxlength Option Comprehensive"

* maxlength(1) - minimal truncation
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "alpha" 2 "beta" 3 "gamma" 4 "delta" 5 "epsilon"
label values x x_lbl
cdecode x, generate(x_str) maxlength(1)
if x_str[1] == "a" & x_str[2] == "b" & x_str[3] == "g" {
    test_pass "maxlength(1)"
}
else {
    test_fail "maxlength(1)" "wrong truncation"
}

* maxlength(2)
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "alpha" 2 "beta" 3 "gamma" 4 "delta" 5 "epsilon"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(2)") maxlength(2)

* maxlength(3)
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "alpha" 2 "beta" 3 "gamma" 4 "delta" 5 "epsilon"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(3)") maxlength(3)

* maxlength(10)
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "this is very long" 2 "short" 3 "medium text" 4 "another long one" 5 "tiny"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(10)") maxlength(10)

* maxlength(20)
clear
set obs 5
gen byte x = _n
local long30 = "a" * 30
label define x_lbl 1 "`long30'" 2 "short" 3 "medium" 4 "`long30'" 5 "x"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(20)") maxlength(20)

* maxlength(50)
clear
set obs 3
gen byte x = _n
local long70 = "b" * 70
label define x_lbl 1 "`long70'" 2 "short" 3 "`long70'"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(50)") maxlength(50)

* maxlength exactly matching label length
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "exact" 2 "match" 3 "five!"
label values x x_lbl
cdecode x, generate(x_str) maxlength(5)
if x_str[1] == "exact" & x_str[2] == "match" & x_str[3] == "five!" {
    test_pass "maxlength exact match"
}
else {
    test_fail "exact match" "unexpected truncation"
}

* maxlength with empty labels
clear
set obs 4
gen byte x = _n
label define x_lbl 1 "" 2 "a" 3 "" 4 "test"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength with empty") maxlength(2)

* Large maxlength (no truncation)
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "short" 2 "medium length" 3 "this is a bit longer"
label values x x_lbl
benchmark_decode x, testname("vs decode: maxlength(100) no truncation") maxlength(100)

/*******************************************************************************
 * SECTION 25: Large datasets
 ******************************************************************************/
print_section "Large Datasets"

* 10K observations with few unique values
clear
set obs 10000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "Category A" 2 "Category B" 3 "Category C" 4 "Category D" 5 "Category E"
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == "Category A"
    if r(N) == 2000 {
        test_pass "10K obs, 5 unique values"
    }
    else {
        test_fail "10K/5" "wrong count"
    }
}
else {
    test_fail "10K/5" "rc=`=_rc'"
}

* 50K observations
clear
set obs 50000
gen byte category = mod(_n, 10) + 1
forvalues i = 1/10 {
    local lbl = "item_`i'"
    if `i' == 1 label define cat_lbl `i' "`lbl'"
    else label define cat_lbl `i' "`lbl'", add
}
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    test_pass "50K observations"
}
else {
    test_fail "50K obs" "rc=`=_rc'"
}

* 100K observations
clear
set obs 100000
gen byte category = mod(_n, 20) + 1
forvalues i = 1/20 {
    local lbl = "label_`i'"
    if `i' == 1 label define cat_lbl `i' "`lbl'"
    else label define cat_lbl `i' "`lbl'", add
}
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == "label_1"
    if r(N) == 5000 {
        test_pass "100K observations"
    }
    else {
        test_fail "100K" "wrong count"
    }
}
else {
    test_fail "100K obs" "rc=`=_rc'"
}

* Large with many unique values (100 labels)
clear
set obs 50000
gen int category = mod(_n, 100) + 1
forvalues i = 1/100 {
    local lbl = "category_`i'"
    if `i' == 1 label define cat_lbl `i' "`lbl'"
    else label define cat_lbl `i' "`lbl'", add
}
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == "category_1"
    if r(N) == 500 {
        test_pass "50K obs, 100 unique values"
    }
    else {
        test_fail "50K/100" "wrong count"
    }
}
else {
    test_fail "50K/100" "rc=`=_rc'"
}

* Large with missing values
clear
set obs 100000
gen byte category = mod(_n, 10) + 1
replace category = . if mod(_n, 5) == 0
label define cat_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "H" 9 "I" 10 "J"
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    count if cat_str == ""
    if r(N) == 20000 {
        test_pass "100K with 20% missing"
    }
    else {
        test_fail "100K missing" "wrong empty count"
    }
}
else {
    test_fail "100K missing" "rc=`=_rc'"
}

* Large dataset comparison with decode
clear
set obs 20000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "Alpha" 2 "Beta" 3 "Gamma" 4 "Delta" 5 "Epsilon"
label values category cat_lbl
benchmark_decode category, testname("vs decode: 20K observations")

/*******************************************************************************
 * SECTION 26: if/in conditions comprehensive
 ******************************************************************************/
print_section "if/in Conditions Comprehensive"

* if with numeric comparison
clear
set obs 100
gen byte x = mod(_n - 1, 5) + 1
gen value = runiform()
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x if value > 0.5, generate(x_str)
count if x_str != "" & value > 0.5
local n_dec = r(N)
count if x_str == "" & value <= 0.5
local n_empty = r(N)
if `n_dec' > 0 & `n_empty' > 0 {
    test_pass "if with numeric comparison"
}
else {
    test_fail "if numeric" "wrong pattern"
}

* if with string comparison (on another variable)
clear
set obs 50
gen byte x = mod(_n - 1, 5) + 1
gen str10 group = cond(mod(_n, 2) == 0, "even", "odd")
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x if group == "even", generate(x_str)
count if x_str != "" & group == "even"
if r(N) == 25 {
    test_pass "if with string comparison"
}
else {
    test_fail "if string" "wrong count"
}

* in first 10
clear
set obs 100
gen byte x = mod(_n - 1, 5) + 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x in 1/10, generate(x_str)
count if x_str != ""
if r(N) == 10 {
    test_pass "in first 10"
}
else {
    test_fail "in first 10" "wrong count"
}

* in last 10
clear
set obs 100
gen byte x = mod(_n - 1, 5) + 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x in 91/100, generate(x_str)
count if x_str != "" in 91/100
if r(N) == 10 {
    test_pass "in last 10"
}
else {
    test_fail "in last 10" "wrong count"
}

* in middle range
clear
set obs 100
gen byte x = mod(_n - 1, 5) + 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x in 40/60, generate(x_str)
count if x_str != ""
if r(N) == 21 {
    test_pass "in middle range"
}
else {
    test_fail "in middle" "wrong count"
}

* Combined if and in - compare with decode
sysuse auto, clear
benchmark_decode foreign, testname("vs decode: if mpg>20 in 1/50") if2("mpg > 20") in2("1/50")

* Complex if condition
sysuse auto, clear
benchmark_decode foreign, testname("vs decode: complex if") if2("price > 5000 & mpg < 25")

* in with single observation
clear
set obs 100
gen byte x = mod(_n - 1, 5) + 1
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x in 50/50, generate(x_str)
count if x_str != ""
if r(N) == 1 {
    test_pass "in single observation"
}
else {
    test_fail "single in" "wrong count"
}

* if that matches no observations
clear
set obs 50
gen byte x = mod(_n - 1, 5) + 1
gen value = _n
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
label values x x_lbl
cdecode x if value > 1000, generate(x_str)
count if x_str != ""
if r(N) == 0 {
    test_pass "if matches no observations"
}
else {
    test_fail "if none" "should have no decoded"
}

* if and in on large dataset
clear
set obs 50000
gen byte x = mod(_n - 1, 10) + 1
gen value = runiform()
forvalues i = 1/10 {
    if `i' == 1 label define x_lbl `i' "cat_`i'"
    else label define x_lbl `i' "cat_`i'", add
}
label values x x_lbl
cdecode x if value > 0.9 in 1/25000, generate(x_str)
count if x_str != ""
* Should be approximately 2500 (10% of 25000)
if r(N) > 2000 & r(N) < 3000 {
    test_pass "if and in on large dataset"
}
else {
    test_fail "large if/in" "unexpected count: `=r(N)'"
}

/*******************************************************************************
 * SECTION 27: Replace option tests
 ******************************************************************************/
print_section "Replace Option"

* Test 27.1: Replace basic functionality
sysuse auto, clear
capture drop foreign_backup
clonevar foreign_backup = foreign
cdecode foreign, replace
capture confirm string variable foreign
if _rc == 0 {
    test_pass "replace converts to string"
}
else {
    test_fail "replace basic" "variable not string"
}

* Test 27.2: Replace preserves value label text
sysuse auto, clear
local lbl1 : label (foreign) 0
cdecode foreign, replace
count if foreign == "`lbl1'" & foreign_backup == 0 in 1
* All foreign==0 should become "Domestic"
sysuse auto, clear
count if foreign == 0
local n_domestic = r(N)
cdecode foreign, replace
count if foreign == "Domestic"
if r(N) == `n_domestic' {
    test_pass "replace preserves label text"
}
else {
    test_fail "replace labels" "wrong count"
}

* Test 27.3: Replace vs decode comparison
sysuse auto, clear
clonevar foreign_c = foreign
clonevar foreign_s = foreign
cdecode foreign_c, replace
decode foreign_s, generate(foreign_s_str)
drop foreign_s
rename foreign_s_str foreign_s
count if foreign_c != foreign_s
if r(N) == 0 {
    test_pass "replace matches decode"
}
else {
    test_fail "replace vs decode" "`=r(N)' values differ"
}
drop foreign_c foreign_s

* Test 27.4: Replace with missing values
clear
set obs 20
gen byte x = _n
replace x = . in 1/5
label define x_lbl 1 "one" 2 "two" 3 "three" 4 "four" 5 "five" ///
    6 "six" 7 "seven" 8 "eight" 9 "nine" 10 "ten" ///
    11 "eleven" 12 "twelve" 13 "thirteen" 14 "fourteen" 15 "fifteen" ///
    16 "sixteen" 17 "seventeen" 18 "eighteen" 19 "nineteen" 20 "twenty"
label values x x_lbl
clonevar x_backup = x
cdecode x, replace
count if x == "" & missing(x_backup)
if r(N) == 5 {
    test_pass "replace handles missings"
}
else {
    test_fail "replace missings" "expected 5 empty"
}

* Test 27.5: Replace with if condition
sysuse auto, clear
cdecode foreign if price > 10000, replace
capture confirm string variable foreign
if _rc == 0 {
    count if foreign != "" & price > 10000
    local n_conv = r(N)
    count if price > 10000
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

* Test 27.6: Replace with in range
sysuse auto, clear
cdecode foreign in 1/30, replace
capture confirm string variable foreign
if _rc == 0 {
    count if foreign != "" in 1/30
    if r(N) == 30 {
        test_pass "replace with in range"
    }
    else {
        test_fail "replace in" "wrong count"
    }
}
else {
    test_fail "replace in" "not converted"
}

* Test 27.7: Replace with maxlength
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "very long label text here" 2 "another long one" ///
    3 "short" 4 "medium length" 5 "x"
label values x x_lbl
cdecode x, replace maxlength(5)
local pass = 1
forvalues i = 1/5 {
    if strlen(x[`i']) > 5 {
        local pass = 0
    }
}
if `pass' {
    test_pass "replace with maxlength"
}
else {
    test_fail "replace maxlen" "strings too long"
}

/*******************************************************************************
 * SECTION 28: Multiple variables (varlist) tests
 ******************************************************************************/
print_section "Multiple Variables (varlist)"

* Test 28.1: Two variables with generate
sysuse auto, clear
capture drop foreign_str rep78_str
cdecode foreign rep78, generate(foreign_str rep78_str)
capture confirm string variable foreign_str
local rc1 = _rc
capture confirm string variable rep78_str
local rc2 = _rc
if `rc1' == 0 & `rc2' == 0 {
    test_pass "two variables with generate"
}
else {
    test_fail "two vars gen" "not created"
}
capture drop foreign_str rep78_str

* Test 28.2: Two variables comparison with decode
sysuse auto, clear
cdecode foreign rep78, generate(foreign_c rep78_c)
decode foreign, generate(foreign_s)
decode rep78, generate(rep78_s)
count if foreign_c != foreign_s
local n1 = r(N)
count if rep78_c != rep78_s
local n2 = r(N)
if `n1' == 0 & `n2' == 0 {
    test_pass "two variables match decode"
}
else {
    test_fail "two vars decode" "foreign diff=`n1' rep78 diff=`n2'"
}
drop foreign_c rep78_c foreign_s rep78_s

* Test 28.3: Three variables
sysuse citytemp, clear
cdecode region division, generate(region_str division_str)
decode region, generate(region_s)
decode division, generate(division_s)
count if region_str != region_s
local n1 = r(N)
count if division_str != division_s
local n2 = r(N)
if `n1' == 0 & `n2' == 0 {
    test_pass "three variables match decode"
}
else {
    test_fail "three vars" "region diff=`n1' division diff=`n2'"
}
drop region_str division_str region_s division_s

* Test 28.4: Multiple variables with replace
sysuse auto, clear
clonevar foreign_backup = foreign
clonevar rep78_backup = rep78
decode foreign, generate(foreign_expected)
decode rep78, generate(rep78_expected)
cdecode foreign rep78, replace
count if foreign != foreign_expected
local n1 = r(N)
count if rep78 != rep78_expected
local n2 = r(N)
if `n1' == 0 & `n2' == 0 {
    test_pass "multiple variables replace"
}
else {
    test_fail "multi replace" "foreign diff=`n1' rep78 diff=`n2'"
}

* Test 28.5: Multiple variables with if condition
sysuse auto, clear
cdecode foreign rep78 if price > 8000, generate(foreign_str rep78_str)
count if foreign_str != "" & price > 8000
local n_conv = r(N)
count if price > 8000
if `n_conv' == r(N) {
    test_pass "multiple vars with if"
}
else {
    test_fail "multi if" "wrong count"
}
drop foreign_str rep78_str

* Test 28.6: Multiple variables count mismatch error
sysuse auto, clear
capture cdecode foreign rep78, generate(only_one)
if _rc != 0 {
    test_pass "error: generate count mismatch"
}
else {
    test_fail "count mismatch" "should error"
}

/*******************************************************************************
 * SECTION 29: Extended missing values comprehensive
 ******************************************************************************/
print_section "Extended Missing Values Comprehensive"

* Test 29.1: All extended missings with labels defined for them
clear
set obs 30
gen x = .
forvalues i = 1/26 {
    local letter = char(96 + `i')
    replace x = .`letter' in `i'
}
replace x = 1 in 27
replace x = 2 in 28
replace x = . in 29
replace x = . in 30
label define x_lbl 1 "one" 2 "two"
label values x x_lbl
* Compare with decode
clonevar x_c = x
clonevar x_s = x
cdecode x_c, replace
decode x_s, generate(x_s_str)
drop x_s
rename x_s_str x_s
count if x_c != x_s
if r(N) == 0 {
    test_pass "extended missings match decode"
}
else {
    test_fail "extended miss" "`=r(N)' differ"
}

* Test 29.2: Extended missings scattered throughout
clear
set obs 100
gen byte x = mod(_n - 1, 10) + 1
replace x = .a if mod(_n, 11) == 0
replace x = .b if mod(_n, 13) == 0
replace x = .z if mod(_n, 17) == 0
forvalues i = 1/10 {
    local lbl = "label_`i'"
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
benchmark_decode x, testname("scattered extended missings")

* Test 29.3: Only extended missing values
clear
set obs 26
gen x = .
forvalues i = 1/26 {
    local letter = char(96 + `i')
    replace x = .`letter' in `i'
}
label define x_lbl 1 "one"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 26 {
    test_pass "only extended missings -> all empty"
}
else {
    test_fail "only extended" "not all empty"
}

/*******************************************************************************
 * SECTION 30: Combined options tests
 ******************************************************************************/
print_section "Combined Options"

* Test 30.1: maxlength with verbose
clear
set obs 10
gen byte x = _n
forvalues i = 1/10 {
    local lbl = "label_" + string(`i')
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
capture cdecode x, generate(x_str) maxlength(5) verbose
if _rc == 0 {
    local pass = 1
    forvalues i = 1/10 {
        if strlen(x_str[`i']) > 5 {
            local pass = 0
        }
    }
    if `pass' {
        test_pass "maxlength with verbose"
    }
    else {
        test_fail "maxlen+verbose" "strings too long"
    }
}
else {
    test_fail "maxlen+verbose" "rc=`=_rc'"
}

* Test 30.2: maxlength with threads
clear
set obs 1000
gen byte x = mod(_n - 1, 10) + 1
forvalues i = 1/10 {
    local lbl = "category_" + string(`i')
    if `i' == 1 label define x_lbl `i' "`lbl'"
    else label define x_lbl `i' "`lbl'", add
}
label values x x_lbl
capture cdecode x, generate(x_str) maxlength(5) threads(2)
if _rc == 0 {
    local pass = 1
    forvalues i = 1/10 {
        if strlen(x_str[`i']) > 5 {
            local pass = 0
        }
    }
    if `pass' {
        test_pass "maxlength with threads"
    }
    else {
        test_fail "maxlen+threads" "strings too long"
    }
}
else {
    test_fail "maxlen+threads" "rc=`=_rc'"
}

* Test 30.3: if with maxlength
sysuse auto, clear
cdecode foreign if mpg > 20, generate(foreign_str) maxlength(3)
count if foreign_str != "" & mpg > 20
local n_conv = r(N)
count if mpg > 20
if `n_conv' == r(N) {
    local pass = 1
    forvalues i = 1/`=_N' {
        if strlen(foreign_str[`i']) > 3 & foreign_str[`i'] != "" {
            local pass = 0
        }
    }
    if `pass' {
        test_pass "if with maxlength"
    }
    else {
        test_fail "if+maxlen" "strings too long"
    }
}
else {
    test_fail "if+maxlen" "wrong count"
}
drop foreign_str

* Test 30.4: in with maxlength
sysuse auto, clear
cdecode foreign in 1/30, generate(foreign_str) maxlength(4)
count if foreign_str != "" in 1/30
if r(N) == 30 {
    local pass = 1
    forvalues i = 1/30 {
        if strlen(foreign_str[`i']) > 4 {
            local pass = 0
        }
    }
    if `pass' {
        test_pass "in with maxlength"
    }
    else {
        test_fail "in+maxlen" "strings too long"
    }
}
else {
    test_fail "in+maxlen" "wrong count"
}
drop foreign_str

* Test 30.5: replace with maxlength
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "very long label one" 2 "long label two" 3 "short" 4 "medium" 5 "x"
label values x x_lbl
cdecode x, replace maxlength(6)
capture confirm string variable x
if _rc == 0 {
    local pass = 1
    forvalues i = 1/5 {
        if strlen(x[`i']) > 6 {
            local pass = 0
        }
    }
    if `pass' {
        test_pass "replace with maxlength"
    }
    else {
        test_fail "replace+maxlen" "strings too long"
    }
}
else {
    test_fail "replace+maxlen" "not string"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
* End of cdecode validation
noi print_summary "cdecode"
}
