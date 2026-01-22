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

noi di as text ""
noi di as text "======================================================================"
noi di as text "              CDECODE VALIDATION TEST SUITE"
noi di as text "======================================================================"

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
        noi test_fail "`testname'" "cdecode rc=`rc_c', decode rc=`rc_s'"
        exit
    }

    if `rc_c' != 0 {
        noi test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare string values
    quietly count if `generate'_c != `generate'_s
    local ndiff = r(N)

    if `ndiff' > 0 {
        noi test_fail "`testname'" "`ndiff' values differ"
    }
    else {
        noi test_pass "`testname'"
    }

    * Cleanup
    capture drop `generate'_c `generate'_s
end

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture noisily cdecode foreign, generate(foreign_str)
if _rc != 0 {
    noi test_fail "cdecode plugin load" "plugin returned error `=_rc'"
    noi print_summary "cdecode"
    exit 1
}
noi test_pass "cdecode plugin loads and runs"
drop foreign_str

/*******************************************************************************
 * SECTION 2: Basic functionality (10 tests)
 ******************************************************************************/
noi print_section "Basic Functionality"

* Test 2.1: Basic decoding creates variable
sysuse auto, clear
capture drop foreign_str
cdecode foreign, generate(foreign_str)
capture confirm string variable foreign_str
if _rc == 0 {
    noi test_pass "basic decoding creates string variable"
}
else {
    noi test_fail "basic decoding" "variable not string"
}

* Test 2.2: Compare with decode
sysuse auto, clear
noi benchmark_decode foreign, testname("vs decode: foreign")

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
    noi test_pass "labels correctly decoded"
}
else {
    noi test_fail "labels decoded" "mismatch found"
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
noi benchmark_decode x, testname("vs decode: var with missings")

* Test 2.5: verbose option
sysuse auto, clear
capture cdecode foreign, generate(foreign_str) verbose
if _rc == 0 {
    noi test_pass "verbose option"
}
else {
    noi test_fail "verbose option" "rc=`=_rc'"
}
capture drop foreign_str

* Test 2.6: threads option
sysuse auto, clear
capture cdecode foreign, generate(foreign_str) threads(2)
if _rc == 0 {
    noi test_pass "threads(2) option"
}
else {
    noi test_fail "threads option" "rc=`=_rc'"
}
capture drop foreign_str

* Test 2.7: Single observation
clear
set obs 1
gen byte x = 1
label define x_lbl 1 "one"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "one" {
    noi test_pass "single observation"
}
else {
    noi test_fail "single observation" "got `=x_str[1]'"
}

* Test 2.8: Two observations
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "first" 2 "second"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "first" & x_str[2] == "second" {
    noi test_pass "two observations"
}
else {
    noi test_fail "two observations" "wrong values"
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
    noi test_pass "many label values"
}
else {
    noi test_fail "many labels" "wrong counts"
}

* Test 2.10: Compare many labels with decode
clear
set obs 50
gen byte category = mod(_n - 1, 5) + 1
label define cat_lbl 1 "Alpha" 2 "Beta" 3 "Gamma" 4 "Delta" 5 "Epsilon"
label values category cat_lbl
noi benchmark_decode category, testname("vs decode: 5 categories")

/*******************************************************************************
 * SECTION 3: Missing values (10 tests)
 ******************************************************************************/
noi print_section "Missing Values"

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
    noi test_pass "missing numeric -> empty string"
}
else {
    noi test_fail "missing handling" "expected `=r(N)', got `n_empty'"
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
    noi test_pass "non-missing values decoded"
}
else {
    noi test_fail "non-missing" "expected `=r(N)', got `n_decoded'"
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
    noi test_pass "all missing -> all empty"
}
else {
    noi test_fail "all missing" "not all empty"
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
    noi test_pass "first observation missing"
}
else {
    noi test_fail "first missing" "wrong handling"
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
    noi test_pass "last observation missing"
}
else {
    noi test_fail "last missing" "wrong handling"
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
    noi test_pass "alternating missing"
}
else {
    noi test_fail "alternating" "expected 10 empty"
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
    noi test_pass "extended missing values"
}
else {
    noi test_fail "extended missing" "not all empty"
}

* Test 3.8: Compare missing handling with decode
clear
set obs 30
gen byte x = mod(_n, 5) + 1
replace x = . if mod(_n, 6) == 0
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: missing values")

* Test 3.9: 90% missing
clear
set obs 100
gen byte x = cond(_n <= 10, mod(_n - 1, 5) + 1, .)
label define x_lbl 1 "a" 2 "b" 3 "c" 4 "d" 5 "e"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 90 {
    noi test_pass "90% missing"
}
else {
    noi test_fail "90% missing" "expected 90 empty"
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
    noi test_pass "only one non-missing"
}
else {
    noi test_fail "singleton" "wrong handling"
}

/*******************************************************************************
 * SECTION 4: Unlabeled values (10 tests)
 ******************************************************************************/
noi print_section "Unlabeled Values"

* Test 4.1: Value without label -> empty string
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "one" 3 "three" 5 "five"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "" & (x == 2 | x == 4)
if r(N) == 2 {
    noi test_pass "unlabeled values -> empty"
}
else {
    noi test_fail "unlabeled" "not empty"
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
    noi test_pass "labeled values decoded"
}
else {
    noi test_fail "labeled values" "wrong decoding"
}

* Test 4.3: Compare unlabeled handling with decode
clear
set obs 10
gen byte x = _n
label define x_lbl 1 "one" 5 "five" 10 "ten"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: sparse labels")

* Test 4.4: All values unlabeled
clear
set obs 5
gen byte x = _n
label define x_lbl 99 "not used"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == ""
if r(N) == 5 {
    noi test_pass "all unlabeled -> all empty"
}
else {
    noi test_fail "all unlabeled" "not all empty"
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
    noi test_pass "negative values labeled"
}
else {
    noi test_fail "negative labels" "wrong decoding"
}

* Test 4.6: Zero labeled
clear
set obs 3
gen byte x = _n - 1  // 0, 1, 2
label define x_lbl 0 "zero" 1 "one" 2 "two"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "zero" {
    noi test_pass "zero labeled"
}
else {
    noi test_fail "zero label" "got `=x_str[1]'"
}

* Test 4.7: Large label values
clear
set obs 3
gen int x = 1000 * _n
label define x_lbl 1000 "thousand" 2000 "two thousand" 3000 "three thousand"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "thousand" & x_str[3] == "three thousand" {
    noi test_pass "large label values"
}
else {
    noi test_fail "large values" "wrong decoding"
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
    noi test_pass "mix of labeled/unlabeled/missing"
}
else {
    noi test_fail "mixed" "wrong counts"
}

* Test 4.9: Compare mixed with decode
clear
set obs 20
gen byte x = mod(_n, 5)
replace x = . if mod(_n, 7) == 0
label define x_lbl 0 "zero" 2 "two" 4 "four"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: mixed labels")

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
    noi test_pass "float values unlabeled"
}
else {
    noi test_fail "float values" "expected all empty"
}

/*******************************************************************************
 * SECTION 5: if/in conditions (10 tests)
 ******************************************************************************/
noi print_section "if/in Conditions"

* Test 5.1: if condition subset
sysuse auto, clear
cdecode foreign if price > 10000, generate(foreign_str)
count if foreign_str != "" & price > 10000
local n_dec = r(N)
count if foreign_str == "" & price <= 10000
local n_empty = r(N)
count if price > 10000
if `n_dec' == r(N) & `n_empty' > 0 {
    noi test_pass "if condition subset"
}
else {
    noi test_fail "if condition" "wrong pattern"
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
    noi test_pass "in range"
}
else {
    noi test_fail "in range" "wrong counts"
}
drop foreign_str

* Test 5.3: Empty if result
sysuse auto, clear
cdecode foreign if price > 100000, generate(foreign_str)
count if foreign_str != ""
if r(N) == 0 {
    noi test_pass "empty if result"
}
else {
    noi test_fail "empty if" "should have no decoded"
}
drop foreign_str

* Test 5.4: if matches all
sysuse auto, clear
cdecode foreign if price > 0, generate(foreign_str)
count if foreign_str != ""
if r(N) == _N {
    noi test_pass "if matches all"
}
else {
    noi test_fail "if all" "should decode all"
}
drop foreign_str

* Test 5.5: Combined if and in
sysuse auto, clear
cdecode foreign if foreign == 1 in 1/50, generate(foreign_str)
count if foreign_str != ""
if r(N) > 0 & r(N) < 50 {
    noi test_pass "combined if and in"
}
else {
    noi test_fail "combined if/in" "wrong count"
}
drop foreign_str

* Test 5.6: Compare if with decode
sysuse auto, clear
noi benchmark_decode foreign, testname("vs decode: if price>8000") if2("price > 8000")

* Test 5.7: Compare in with decode
sysuse auto, clear
noi benchmark_decode foreign, testname("vs decode: in 1/30") in2("1/30")

* Test 5.8: Single observation in
sysuse auto, clear
cdecode foreign in 1/1, generate(foreign_str)
count if foreign_str != ""
if r(N) == 1 {
    noi test_pass "single observation in"
}
else {
    noi test_fail "single in" "wrong count"
}
drop foreign_str

* Test 5.9: if with condition on another variable
sysuse auto, clear
cdecode foreign if price < 10000, generate(foreign_str)
* Should decode foreign where price < 10000
count if foreign_str != "" & price < 10000
if r(N) > 0 {
    noi test_pass "if with condition on price"
}
else {
    noi test_fail "price condition if" "no decoded"
}
drop foreign_str

* Test 5.10: in from middle
sysuse auto, clear
cdecode foreign in 30/50, generate(foreign_str)
count if foreign_str != "" in 30/50
if r(N) == 21 {
    noi test_pass "in from middle"
}
else {
    noi test_fail "middle in" "wrong count"
}
drop foreign_str

/*******************************************************************************
 * SECTION 6: maxlength option (10 tests)
 ******************************************************************************/
noi print_section "maxlength Option"

* Test 6.1: maxlength truncates long labels
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "this is a very long label" 2 "short" 3 "another long one here"
label values x x_lbl
cdecode x, generate(x_str) maxlength(10)
if strlen(x_str[1]) <= 10 {
    noi test_pass "maxlength truncates"
}
else {
    noi test_fail "maxlength" "string too long: `=strlen(x_str[1])'"
}

* Test 6.2: Short labels unaffected by maxlength
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "abc" 2 "de" 3 "f"
label values x x_lbl
cdecode x, generate(x_str) maxlength(100)
if x_str[1] == "abc" & x_str[2] == "de" & x_str[3] == "f" {
    noi test_pass "short labels unaffected"
}
else {
    noi test_fail "short labels" "wrong values"
}

* Test 6.3: maxlength(1)
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "alpha" 2 "beta" 3 "gamma"
label values x x_lbl
cdecode x, generate(x_str) maxlength(1)
if strlen(x_str[1]) == 1 & strlen(x_str[2]) == 1 {
    noi test_pass "maxlength(1)"
}
else {
    noi test_fail "maxlength(1)" "not truncated to 1"
}

* Test 6.4: Compare maxlength with decode
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "very long label text here" 2 "another long one" ///
    3 "short" 4 "medium length" 5 "x"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: maxlength(5)") maxlength(5)

* Test 6.5: maxlength with exact length match
clear
set obs 1
gen byte x = 1
label define x_lbl 1 "exact"
label values x x_lbl
cdecode x, generate(x_str) maxlength(5)
if x_str[1] == "exact" {
    noi test_pass "maxlength exact match"
}
else {
    noi test_fail "exact match" "got `=x_str[1]'"
}

* Test 6.6: maxlength preserves prefix
clear
set obs 1
gen byte x = 1
label define x_lbl 1 "abcdefghij"
label values x x_lbl
cdecode x, generate(x_str) maxlength(5)
if x_str[1] == "abcde" {
    noi test_pass "maxlength preserves prefix"
}
else {
    noi test_fail "prefix" "got `=x_str[1]'"
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
    noi test_pass "default width from labels"
}
else {
    noi test_fail "default width" "wrong length"
}

* Test 6.8: maxlength larger than any label
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "a" 2 "bb" 3 "ccc"
label values x x_lbl
cdecode x, generate(x_str) maxlength(1000)
if x_str[1] == "a" & x_str[2] == "bb" & x_str[3] == "ccc" {
    noi test_pass "maxlength larger than labels"
}
else {
    noi test_fail "large maxlength" "wrong values"
}

* Test 6.9: Compare various maxlength with decode
clear
set obs 10
gen byte x = mod(_n, 3) + 1
label define x_lbl 1 "category_one" 2 "category_two" 3 "category_three"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: maxlength(8)") maxlength(8)

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
    noi test_pass "maxlength with empties"
}
else {
    noi test_fail "maxlength empties" "wrong handling"
}

/*******************************************************************************
 * SECTION 7: Large datasets (10 tests)
 ******************************************************************************/
noi print_section "Large Datasets"

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
        noi test_pass "100k observations"
    }
    else {
        noi test_fail "100k obs" "wrong count"
    }
}
else {
    noi test_fail "100k obs" "rc=`=_rc'"
}

* Test 7.2: 500k observations
clear
set obs 500000
gen byte category = mod(_n, 10) + 1
forvalues i = 1/10 {
    local lbl = "category_`i'"
    if `i' == 1 label define cat_lbl `i' "`lbl'"
    else label define cat_lbl `i' "`lbl'", add
}
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    noi test_pass "500k observations"
}
else {
    noi test_fail "500k obs" "rc=`=_rc'"
}

* Test 7.3: 1M observations
clear
set obs 1000000
gen byte category = mod(_n, 3) + 1
label define cat_lbl 1 "One" 2 "Two" 3 "Three"
label values category cat_lbl
capture cdecode category, generate(cat_str)
if _rc == 0 {
    noi test_pass "1M observations"
}
else {
    noi test_fail "1M obs" "rc=`=_rc'"
}

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
    noi test_pass "large with if condition"
}
else {
    noi test_fail "large if" "wrong count"
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
    noi test_pass "large with long labels"
}
else {
    noi test_fail "long labels" "rc=`=_rc'"
}

* Test 7.6: Large with threads(4)
clear
set obs 500000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values category cat_lbl
capture cdecode category, generate(cat_str) threads(4)
if _rc == 0 {
    noi test_pass "large with threads(4)"
}
else {
    noi test_fail "threads(4)" "rc=`=_rc'"
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
        noi test_pass "large with many missing"
    }
    else {
        noi test_fail "many missing" "wrong empty count"
    }
}
else {
    noi test_fail "many missing" "rc=`=_rc'"
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
    noi test_pass "large dataset, small sample"
}
else {
    noi test_fail "small sample" "wrong count"
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
        noi test_pass "50 unique labels"
    }
    else {
        noi test_fail "50 labels" "wrong count"
    }
}
else {
    noi test_fail "50 labels" "rc=`=_rc'"
}

* Test 7.10: Compare large dataset with decode
clear
set obs 50000
gen byte category = mod(_n, 5) + 1
label define cat_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values category cat_lbl
noi benchmark_decode category, testname("vs decode: 50k obs")

/*******************************************************************************
 * SECTION 8: Special characters in labels (10 tests)
 ******************************************************************************/
noi print_section "Special Characters in Labels"

* Test 8.1: Spaces in labels
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "hello world" 2 "foo bar baz" 3 "a b c"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "hello world" {
    noi test_pass "spaces in labels"
}
else {
    noi test_fail "spaces" "got `=x_str[1]'"
}

* Test 8.2: Leading/trailing spaces
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "  leading" 2 "trailing  " 3 "  both  "
label values x x_lbl
cdecode x, generate(x_str)
* Note: Stata's decode may trim spaces
capture noi benchmark_decode x, testname("vs decode: leading/trailing spaces")

* Test 8.3: Commas in labels
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "one, two, three" 2 "a,b"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], ",") > 0 {
    noi test_pass "commas in labels"
}
else {
    noi test_fail "commas" "comma not preserved"
}

* Test 8.4: Parentheses and brackets
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "(parens)" 2 "[brackets]" 3 "{braces}"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "(parens)" & x_str[2] == "[brackets]" {
    noi test_pass "parentheses and brackets"
}
else {
    noi test_fail "brackets" "not preserved"
}

* Test 8.5: Ampersand and percent
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "A & B" 2 "50%"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], "&") > 0 & strpos(x_str[2], "%") > 0 {
    noi test_pass "ampersand and percent"
}
else {
    noi test_fail "amp/percent" "not preserved"
}

* Test 8.6: Apostrophes
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "don't" 2 "it's"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], "'") > 0 {
    noi test_pass "apostrophes"
}
else {
    noi test_fail "apostrophes" "not preserved"
}

* Test 8.7: Slashes
clear
set obs 2
gen byte x = _n
label define x_lbl 1 "path/to/file" 2 "a\b"
label values x x_lbl
cdecode x, generate(x_str)
if strpos(x_str[1], "/") > 0 {
    noi test_pass "forward slashes"
}
else {
    noi test_fail "slashes" "not preserved"
}

* Test 8.8: Compare special chars with decode
clear
set obs 5
gen byte x = _n
label define x_lbl 1 "Hello!" 2 "#hashtag" 3 "$100" 4 "a+b=c" 5 "x*y"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: special chars")

* Test 8.9: Dashes and underscores
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "with-dash" 2 "with_underscore" 3 "both-and_here"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "with-dash" & x_str[2] == "with_underscore" {
    noi test_pass "dashes and underscores"
}
else {
    noi test_fail "dash/underscore" "not preserved"
}

* Test 8.10: Mixed special characters
clear
set obs 3
gen byte x = _n
label define x_lbl 1 "Item #1 (new)" 2 "50% off - sale!" 3 "A/B test: v2.0"
label values x x_lbl
cdecode x, generate(x_str)
noi benchmark_decode x, testname("vs decode: mixed special")

/*******************************************************************************
 * SECTION 9: Real-world datasets (10 tests)
 ******************************************************************************/
noi print_section "Real-World Datasets"

* Test 9.1: auto - foreign
sysuse auto, clear
noi benchmark_decode foreign, testname("auto: foreign")

* Test 9.2: citytemp - division
sysuse citytemp, clear
noi benchmark_decode division, testname("citytemp: division")

* Test 9.3: census - region
sysuse census, clear
noi benchmark_decode region, testname("census: region")

* Test 9.4: nlswork - race
webuse nlswork, clear
noi benchmark_decode race, testname("nlswork: race")

* Test 9.5: nlswork - msp
webuse nlswork, clear
noi benchmark_decode msp, testname("nlswork: msp")

* Test 9.6: citytemp - region
sysuse citytemp, clear
noi benchmark_decode region, testname("citytemp: region")

* Test 9.7: pop2000 - agegrp
sysuse pop2000, clear
noi benchmark_decode agegrp, testname("pop2000: agegrp")

* Test 9.8: voter - candidat
sysuse voter, clear
noi benchmark_decode candidat, testname("voter: candidat")

* Test 9.9: voter - inc
sysuse voter, clear
noi benchmark_decode inc, testname("voter: inc")

* Test 9.10: bpwide - agegrp
sysuse bpwide, clear
noi benchmark_decode agegrp, testname("bpwide: agegrp")

/*******************************************************************************
 * SECTION 10: Pathological data (10 tests)
 ******************************************************************************/
noi print_section "Pathological Data"

* Test 10.1: All same value
clear
set obs 1000
gen byte x = 1
label define x_lbl 1 "constant"
label values x x_lbl
cdecode x, generate(x_str)
count if x_str == "constant"
if r(N) == 1000 {
    noi test_pass "all same value"
}
else {
    noi test_fail "all same" "wrong count"
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
    noi test_pass "100 unique values"
}
else {
    noi test_fail "100 unique" "wrong count"
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
    noi test_pass "binary values"
}
else {
    noi test_fail "binary" "wrong counts"
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
        noi test_pass "very long label text"
    }
    else {
        noi test_fail "long text" "truncated too much"
    }
}
else {
    noi test_fail "long text" "rc=`=_rc'"
}

* Test 10.5: Bookend pattern
clear
set obs 100
gen byte x = cond(_n == 1 | _n == 100, 1, 2)
label define x_lbl 1 "bookend" 2 "middle"
label values x x_lbl
cdecode x, generate(x_str)
if x_str[1] == "bookend" & x_str[100] == "bookend" & x_str[50] == "middle" {
    noi test_pass "bookend pattern"
}
else {
    noi test_fail "bookend" "wrong values"
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
    noi test_pass "cyclic pattern"
}
else {
    noi test_fail "cyclic" "wrong count"
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
    noi test_pass "sparse data"
}
else {
    noi test_fail "sparse" "wrong count"
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
    noi test_pass "shuffled data"
}
else {
    noi test_fail "shuffled" "wrong count"
}

* Test 10.9: Compare pathological with decode
clear
set obs 50
gen byte x = cond(mod(_n, 10) == 0, 1, cond(mod(_n, 10) == 5, ., 2))
label define x_lbl 1 "rare" 2 "common"
label values x x_lbl
noi benchmark_decode x, testname("vs decode: pathological")

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
    noi test_pass "extreme value range"
}
else {
    noi test_fail "extreme range" "wrong values"
}

/*******************************************************************************
 * SECTION 11: Error handling (5 tests)
 ******************************************************************************/
noi print_section "Error Handling"

* Test 11.1: Generate var exists
sysuse auto, clear
gen foreign_str = "test"
capture cdecode foreign, generate(foreign_str)
if _rc == 110 {
    noi test_pass "error: generate var exists"
}
else {
    noi test_fail "var exists" "expected rc=110, got `=_rc'"
}

* Test 11.2: Source is string (not numeric)
sysuse auto, clear
capture cdecode make, generate(make_decoded)
if _rc != 0 {
    noi test_pass "error: string source"
}
else {
    noi test_fail "string source" "should error"
}

* Test 11.3: Source has no value label
clear
set obs 10
gen x = _n
capture cdecode x, generate(x_str)
if _rc != 0 {
    noi test_pass "error: no value label"
}
else {
    noi test_fail "no label" "should error"
}

* Test 11.4: Source doesn't exist
sysuse auto, clear
capture cdecode nonexistent, generate(ne_str)
if _rc != 0 {
    noi test_pass "error: nonexistent source"
}
else {
    noi test_fail "nonexistent" "should error"
}

* Test 11.5: Empty dataset
clear
set obs 0
gen byte x = .
label define x_lbl 1 "one"
label values x x_lbl
capture cdecode x, generate(x_str)
if _rc == 0 | _rc == 2000 {
    noi test_pass "empty dataset handled"
}
else {
    noi test_fail "empty dataset" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 12: Consistency/round-trip (5 tests)
 ******************************************************************************/
noi print_section "Consistency/Round-trip"

* Test 12.1: Round-trip encode then decode
sysuse auto, clear
encode make, generate(make_num)
decode make_num, generate(make_back)
local match = 1
forvalues i = 1/10 {
    if make[`i'] != make_back[`i'] local match = 0
}
if `match' {
    noi test_pass "round-trip encode/decode"
}
else {
    noi test_fail "round-trip" "mismatch"
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
    noi test_pass "round-trip cencode/cdecode"
}
else {
    noi test_fail "cencode/cdecode" "mismatch"
}

* Test 12.3: cdecode matches decode
sysuse auto, clear
decode foreign, generate(foreign_d)
cdecode foreign, generate(foreign_c)
count if foreign_d != foreign_c
if r(N) == 0 {
    noi test_pass "cdecode matches decode"
}
else {
    noi test_fail "match decode" "`=r(N)' differ"
}

* Test 12.4: Observation count preserved
sysuse auto, clear
local orig_n = _N
cdecode foreign, generate(foreign_str)
if _N == `orig_n' {
    noi test_pass "observation count preserved"
}
else {
    noi test_fail "obs count" "changed"
}

* Test 12.5: Full comparison on census
sysuse census, clear
noi benchmark_decode region, testname("full comparison: census region")

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cdecode"

if $TESTS_FAILED > 0 {
    exit 1
}

}
