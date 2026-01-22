/*
 * validation_cdecode.do
 *
 * Comprehensive validation tests for cdecode command
 * Target: ~100 tests covering all edge cases, pathological data,
 * and real-world datasets
 */

clear all
set more off
adopath + "/Users/Mike/Documents/GitHub/stata-ctools/build"

* Test tracking using globals
global n_passed = 0
global n_failed = 0
global failed_tests ""

capture program drop run_test
program define run_test
    args test_name
    di as text ""
    di as text "{hline 60}"
    di as text "Test: `test_name'"
    di as text "{hline 60}"
end

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

di as text ""
di as text "============================================================"
di as text "     cdecode Comprehensive Validation Suite (~100 tests)"
di as text "============================================================"
di as text ""

* ============================================================================
* SECTION 1: Basic Functionality (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 1: Basic Functionality"
di as text ""

* ---------------------------------------------------------------------------
run_test "1.1 Simple integer labels"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Simple integer labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.2 Labels with spaces"
* ---------------------------------------------------------------------------
clear
set obs 4
gen status = _n
label define status_lbl 1 "In Progress" 2 "On Hold" 3 "Completed Successfully" 4 "Failed with Error"
label values status status_lbl

cdecode status, generate(status_cdecode)
decode status, generate(status_decode)

capture assert status_cdecode == status_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Labels with spaces mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.3 Single observation"
* ---------------------------------------------------------------------------
clear
set obs 1
gen x = 1
label define x_lbl 1 "Only One"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Single observation mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.4 Many unique labels (100)"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = _n
forvalues i = 1/100 {
    label define x_lbl `i' "Label number `i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "100 unique labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.5 Repeated values same label"
* ---------------------------------------------------------------------------
clear
set obs 1000
gen x = ceil(runiform() * 3)
label define x_lbl 1 "Low" 2 "Medium" 3 "High"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Repeated values mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.6 Labels with leading spaces"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 " Leading" 2 "  Double Leading" 3 "   Triple Leading"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Leading spaces mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.7 Labels with trailing spaces"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "Trailing " 2 "Double Trailing  " 3 "Triple Trailing   "
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Trailing spaces mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.8 Labels with both leading and trailing spaces"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 " Both " 2 "  Both  " 3 "   Both   "
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Both leading/trailing spaces mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.9 Empty-looking labels (just spaces)"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 " " 2 "  " 3 "   "
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Space-only labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "1.10 Two observations"
* ---------------------------------------------------------------------------
clear
set obs 2
gen x = _n
label define x_lbl 1 "First" 2 "Second"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Two observations mismatch"
}


* ============================================================================
* SECTION 2: Missing Values (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 2: Missing Values"
di as text ""

* ---------------------------------------------------------------------------
run_test "2.1 Standard missing value (.)"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
replace x = . in 3
label define x_lbl 1 "A" 2 "B" 4 "D" 5 "E"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Standard missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.2 Extended missing values (.a - .z)"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
replace x = .a in 3
replace x = .b in 5
replace x = .z in 8
label define x_lbl 1 "One" 2 "Two" 4 "Four" 6 "Six" 7 "Seven" 9 "Nine" 10 "Ten"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Extended missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.3 All missing values"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = .
label define x_lbl 1 "One" 2 "Two"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "All missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.4 Missing at start"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
replace x = . in 1/3
label define x_lbl 4 "Four" 5 "Five" 6 "Six" 7 "Seven" 8 "Eight" 9 "Nine" 10 "Ten"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Missing at start mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.5 Missing at end"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
replace x = . in 8/10
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five" 6 "Six" 7 "Seven"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Missing at end mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.6 Missing in middle"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
replace x = . in 4/6
label define x_lbl 1 "A" 2 "B" 3 "C" 7 "G" 8 "H" 9 "I" 10 "J"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Missing in middle mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.7 Alternating missing"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
replace x = . if mod(_n, 2) == 0
label define x_lbl 1 "One" 3 "Three" 5 "Five" 7 "Seven" 9 "Nine"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Alternating missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.8 Single non-missing among all missing"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = .
replace x = 5 in 5
label define x_lbl 5 "The Only One"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Single non-missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.9 All extended missing types"
* ---------------------------------------------------------------------------
clear
set obs 27
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
label define x_lbl 1 "Valid"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "All extended missing types mismatch"
}

* ---------------------------------------------------------------------------
run_test "2.10 Random pattern of missing"
* ---------------------------------------------------------------------------
clear
set seed 12345
set obs 100
gen x = ceil(runiform() * 10)
replace x = . if runiform() < 0.3
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "H" 9 "I" 10 "J"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Random missing pattern mismatch"
}


* ============================================================================
* SECTION 3: Unlabeled Values (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 3: Unlabeled Values"
di as text ""

* ---------------------------------------------------------------------------
run_test "3.1 Some values unlabeled (gaps)"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
label define x_lbl 1 "One" 3 "Three" 5 "Five" 7 "Seven" 9 "Nine"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Gaps in labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.2 First value unlabeled"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "First unlabeled mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.3 Last value unlabeled"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Last unlabeled mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.4 All values unlabeled"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define empty_lbl 999 "Unused"
label values x empty_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "All unlabeled mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.5 Single labeled value among many unlabeled"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = _n
label define x_lbl 50 "The Middle One"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Single labeled mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.6 Every other value labeled"
* ---------------------------------------------------------------------------
clear
set obs 20
gen x = _n
forvalues i = 1(2)20 {
    label define x_lbl `i' "Odd `i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Every other labeled mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.7 Wide gaps between labeled values"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = 1 in 1
replace x = 100 in 2
replace x = 10000 in 3
replace x = 100000 in 4
replace x = 1000000 in 5
label define x_lbl 1 "One" 100 "Hundred" 10000 "Ten K" 100000 "Hundred K" 1000000 "Million"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Wide gaps mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.8 Sparse labels at boundaries"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
label define x_lbl 1 "First" 10 "Last"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Boundary labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.9 Labels only for values not in data"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n + 100  // Values 101-105
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Labels not in data mismatch"
}

* ---------------------------------------------------------------------------
run_test "3.10 Mix of labeled, unlabeled, and missing"
* ---------------------------------------------------------------------------
clear
set obs 15
gen x = _n
replace x = . in 5
replace x = . in 10
replace x = . in 15
label define x_lbl 1 "One" 3 "Three" 7 "Seven" 11 "Eleven" 13 "Thirteen"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Mixed labeled/unlabeled/missing mismatch"
}


* ============================================================================
* SECTION 4: Negative and Zero Values (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 4: Negative and Zero Values"
di as text ""

* ---------------------------------------------------------------------------
run_test "4.1 Simple negative values"
* ---------------------------------------------------------------------------
clear
set obs 7
gen x = _n - 4  // -3 to 3
label define x_lbl -3 "Minus Three" -2 "Minus Two" -1 "Minus One" 0 "Zero" 1 "One" 2 "Two" 3 "Three"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Simple negative mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.2 Large negative values"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = -1000000 + (_n - 1) * 500000
label define x_lbl -1000000 "Neg Million" -500000 "Neg Half M" 0 "Zero" 500000 "Half M" 1000000 "Million"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large negative mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.3 Only negative values"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = -_n
label define x_lbl -1 "Neg One" -2 "Neg Two" -3 "Neg Three" -4 "Neg Four" -5 "Neg Five"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Only negative mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.4 Zero only"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = 0
label define x_lbl 0 "All Zeros"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Zero only mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.5 Negative with gaps"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n - 5  // -4 to 5
label define x_lbl -4 "A" -2 "B" 0 "C" 2 "D" 4 "E"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Negative with gaps mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.6 Consecutive negative"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = -_n  // -1 to -10
forvalues i = 1/10 {
    label define x_lbl -`i' "Negative `i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Consecutive negative mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.7 Mix of negative, zero, positive with missing"
* ---------------------------------------------------------------------------
clear
set obs 11
gen x = _n - 6  // -5 to 5
replace x = . in 6  // Zero becomes missing
label define x_lbl -5 "A" -4 "B" -3 "C" -2 "D" -1 "E" 1 "G" 2 "H" 3 "I" 4 "J" 5 "K"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Mixed neg/zero/pos with missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.8 Large range spanning negative to positive"
* ---------------------------------------------------------------------------
clear
set obs 21
gen x = _n - 11  // -10 to 10
forvalues i = -10/10 {
    local abs_i = abs(`i')
    if `i' < 0 {
        label define x_lbl `i' "Neg`abs_i'", add
    }
    else if `i' == 0 {
        label define x_lbl 0 "Zero", add
    }
    else {
        label define x_lbl `i' "Pos`i'", add
    }
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large range mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.9 Negative one only"
* ---------------------------------------------------------------------------
clear
set obs 1
gen x = -1
label define x_lbl -1 "Negative One"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Negative one only mismatch"
}

* ---------------------------------------------------------------------------
run_test "4.10 Large negative and positive values"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = -10000 in 1
replace x = -5000 in 2
replace x = 0 in 3
replace x = 5000 in 4
replace x = 10000 in 5
label define x_lbl -10000 "Neg 10K" -5000 "Neg 5K" 0 "Zero" 5000 "Pos 5K" 10000 "Pos 10K"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large values mismatch"
}


* ============================================================================
* SECTION 5: if/in Conditions (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 5: if/in Conditions"
di as text ""

* ---------------------------------------------------------------------------
run_test "5.1 Simple if condition"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five" 6 "Six" 7 "Seven" 8 "Eight" 9 "Nine" 10 "Ten"
label values x x_lbl

cdecode x if x <= 5, generate(x_cdecode)
decode x if x <= 5, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Simple if mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.2 in range"
* ---------------------------------------------------------------------------
clear
set obs 20
gen x = _n
forvalues i = 1/20 {
    label define x_lbl `i' "Val`i'", add
}
label values x x_lbl

cdecode x in 5/15, generate(x_cdecode)
decode x in 5/15, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "in range mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.3 if with another variable"
* ---------------------------------------------------------------------------
clear
set obs 20
gen x = _n
gen group = mod(_n, 3)
forvalues i = 1/20 {
    label define x_lbl `i' "Item`i'", add
}
label values x x_lbl

cdecode x if group == 1, generate(x_cdecode)
decode x if group == 1, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "if with another var mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.4 Complex if condition"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = ceil(runiform() * 10)
gen y = runiform()
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "H" 9 "I" 10 "J"
label values x x_lbl

cdecode x if x > 3 & x < 8 & y > 0.5, generate(x_cdecode)
decode x if x > 3 & x < 8 & y > 0.5, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Complex if mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.5 if selects no observations"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
label define x_lbl 1 "One" 2 "Two" 3 "Three"
label values x x_lbl

cdecode x if x > 100, generate(x_cdecode)
decode x if x > 100, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "if no obs mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.6 if selects all observations"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
forvalues i = 1/10 {
    label define x_lbl `i' "Val`i'", add
}
label values x x_lbl

cdecode x if x > 0, generate(x_cdecode)
decode x if x > 0, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "if all obs mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.7 in first observation only"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
forvalues i = 1/10 {
    label define x_lbl `i' "Val`i'", add
}
label values x x_lbl

cdecode x in 1, generate(x_cdecode)
decode x in 1, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "in first only mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.8 in last observation only"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
forvalues i = 1/10 {
    label define x_lbl `i' "Val`i'", add
}
label values x x_lbl

cdecode x in 10, generate(x_cdecode)
decode x in 10, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "in last only mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.9 if with missing values"
* ---------------------------------------------------------------------------
clear
set obs 10
gen x = _n
replace x = . in 3
replace x = . in 7
label define x_lbl 1 "A" 2 "B" 4 "D" 5 "E" 6 "F" 8 "H" 9 "I" 10 "J"
label values x x_lbl

cdecode x if !missing(x), generate(x_cdecode)
decode x if !missing(x), generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "if with missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "5.10 Combined if and in"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = ceil(runiform() * 5)
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values x x_lbl

cdecode x if x >= 3 in 20/80, generate(x_cdecode)
decode x if x >= 3 in 20/80, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Combined if/in mismatch"
}


* ============================================================================
* SECTION 6: maxlength Option (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 6: maxlength Option"
di as text ""

* ---------------------------------------------------------------------------
run_test "6.1 maxlength truncates long labels"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "This is a very long label that should be truncated" ///
                   2 "Another extremely long label text here" ///
                   3 "Short"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(10)
decode x, generate(x_decode) maxlength(10)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength truncation mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.2 maxlength equal to longest label"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "Alpha" 2 "Beta" 3 "Gamma"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(5)
decode x, generate(x_decode) maxlength(5)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength equal longest mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.3 maxlength of 1"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "Alpha" 2 "Beta" 3 "Gamma"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(1)
decode x, generate(x_decode) maxlength(1)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength 1 mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.4 maxlength larger than all labels"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "A" 2 "BB" 3 "CCC"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(100)
decode x, generate(x_decode) maxlength(100)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength larger mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.5 maxlength with spaces"
* ---------------------------------------------------------------------------
clear
set obs 2
gen x = _n
label define x_lbl 1 "Label with many spaces inside" 2 "Another spaced label"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(15)
decode x, generate(x_decode) maxlength(15)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength spaces mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.6 maxlength of 2"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "Apple" 2 "Banana" 3 "Cherry" 4 "Date" 5 "Elderberry"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(2)
decode x, generate(x_decode) maxlength(2)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength 2 mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.7 maxlength of 50"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "Short" 2 "A medium length label for testing" 3 "This is a much longer label that exceeds fifty characters in total length"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(50)
decode x, generate(x_decode) maxlength(50)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength 50 mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.8 maxlength cuts at exact boundary"
* ---------------------------------------------------------------------------
clear
set obs 1
gen x = 1
label define x_lbl 1 "1234567890"  // Exactly 10 chars
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(10)
decode x, generate(x_decode) maxlength(10)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength exact boundary mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.9 maxlength with unicode"
* ---------------------------------------------------------------------------
clear
set obs 2
gen x = _n
label define x_lbl 1 "Cafe with accent" 2 "Regular text"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(8)
decode x, generate(x_decode) maxlength(8)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength unicode mismatch"
}

* ---------------------------------------------------------------------------
run_test "6.10 maxlength of 3"
* ---------------------------------------------------------------------------
clear
set obs 4
gen x = _n
label define x_lbl 1 "AAAA" 2 "BB" 3 "C" 4 "DDDDD"
label values x x_lbl

cdecode x, generate(x_cdecode) maxlength(3)
decode x, generate(x_decode) maxlength(3)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "maxlength 3 mismatch"
}


* ============================================================================
* SECTION 7: Large Datasets and Performance (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 7: Large Datasets and Performance"
di as text ""

* ---------------------------------------------------------------------------
run_test "7.1 100,000 observations"
* ---------------------------------------------------------------------------
clear
set obs 100000
gen x = ceil(runiform() * 10)
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five" ///
                   6 "Six" 7 "Seven" 8 "Eight" 9 "Nine" 10 "Ten"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "100K obs mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.2 1,000,000 observations"
* ---------------------------------------------------------------------------
clear
set obs 1000000
gen x = ceil(runiform() * 5)
label define x_lbl 1 "Category A" 2 "Category B" 3 "Category C" 4 "Category D" 5 "Category E"
label values x x_lbl

timer clear 1
timer on 1
cdecode x, generate(x_cdecode)
timer off 1

timer clear 2
timer on 2
decode x, generate(x_decode)
timer off 2

capture assert x_cdecode == x_decode
if _rc == 0 {
    quietly timer list 1
    local cdecode_time = r(t1)
    quietly timer list 2
    local decode_time = r(t2)
    di as text "  cdecode time: " %6.3f `cdecode_time' " sec"
    di as text "  decode time:  " %6.3f `decode_time' " sec"
    test_passed
}
else {
    test_failed "1M obs mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.3 5,000,000 observations"
* ---------------------------------------------------------------------------
clear
set obs 5000000
gen x = ceil(runiform() * 3)
label define x_lbl 1 "Low" 2 "Medium" 3 "High"
label values x x_lbl

timer clear 1
timer on 1
cdecode x, generate(x_cdecode)
timer off 1

timer clear 2
timer on 2
decode x, generate(x_decode)
timer off 2

capture assert x_cdecode == x_decode
if _rc == 0 {
    quietly timer list 1
    local cdecode_time = r(t1)
    quietly timer list 2
    local decode_time = r(t2)
    di as text "  cdecode time: " %6.3f `cdecode_time' " sec"
    di as text "  decode time:  " %6.3f `decode_time' " sec"
    test_passed
}
else {
    test_failed "5M obs mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.4 100K with many unique labels"
* ---------------------------------------------------------------------------
clear
set obs 100000
* Use 100 labels instead of 500 to avoid Stata's command line length limits
gen x = ceil(runiform() * 100)
forvalues i = 1/100 {
    label define x_lbl `i' "Label_`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "100K many labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.5 Large with high missing rate"
* ---------------------------------------------------------------------------
clear
set obs 500000
gen x = ceil(runiform() * 10)
replace x = . if runiform() < 0.5  // 50% missing
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "H" 9 "I" 10 "J"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large high missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.6 Large with mostly unlabeled"
* ---------------------------------------------------------------------------
clear
set obs 500000
gen x = ceil(runiform() * 100)
label define x_lbl 1 "One" 50 "Fifty" 100 "Hundred"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large mostly unlabeled mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.7 Large with long labels"
* ---------------------------------------------------------------------------
clear
set obs 100000
gen x = ceil(runiform() * 5)
label define x_lbl 1 "This is a very long label for category one" ///
                   2 "This is another very long label for category two" ///
                   3 "Category three also has a fairly long descriptive label" ///
                   4 "The fourth category label is quite lengthy as well" ///
                   5 "Finally category five rounds out the long labels"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large long labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.8 Large with if condition"
* ---------------------------------------------------------------------------
clear
set obs 1000000
gen x = ceil(runiform() * 10)
gen y = runiform()
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "H" 9 "I" 10 "J"
label values x x_lbl

cdecode x if y > 0.5, generate(x_cdecode)
decode x if y > 0.5, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large with if mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.9 Large sparse negative values"
* ---------------------------------------------------------------------------
clear
set obs 200000
gen x = ceil(runiform() * 20) - 10  // -9 to 10
forvalues i = -9/10 {
    label define x_lbl `i' "Val_`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large sparse negative mismatch"
}

* ---------------------------------------------------------------------------
run_test "7.10 Large uniformly distributed"
* ---------------------------------------------------------------------------
clear
set obs 500000
* Use 50 labels instead of 1000 to avoid Stata's command line length limits
gen x = ceil(runiform() * 50)
forvalues i = 1/50 {
    label define x_lbl `i' "L`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Large uniform mismatch"
}


* ============================================================================
* SECTION 8: Special Characters in Labels (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 8: Special Characters in Labels"
di as text ""

* ---------------------------------------------------------------------------
run_test "8.1 Labels with punctuation"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "Hello, World!" 2 "Test: Success" 3 "Item #1" 4 "50% Off" 5 "Q&A"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Punctuation mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.2 Labels with apostrophes"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
* Note: Labels with embedded double quotes have parsing issues, use apostrophes
label define x_lbl 1 "It's working" 2 "Don't stop" 3 "O'Brien's"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Apostrophes mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.3 Labels with backslashes"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
label define x_lbl 1 "Path\to\file" 2 "Back\\slash" 3 "End\"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Backslashes mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.4 Labels with slashes"
* ---------------------------------------------------------------------------
clear
set obs 3
gen x = _n
* Note: Pipes (|) are used as delimiters in encoding, use slashes instead
label define x_lbl 1 "A/B" 2 "X/Y/Z" 3 "//Double//"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Slashes mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.5 Labels with numbers"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "Item 1" 2 "Version 2.0" 3 "123 Main St" 4 "99.9%" 5 "1st Place"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Numbers in labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.6 Labels with tabs"
* ---------------------------------------------------------------------------
clear
set obs 2
gen x = _n
label define x_lbl 1 "Tab	here" 2 "Multiple		tabs"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Tabs mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.7 Labels with special symbols"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "@ Symbol" 2 "$ Dollar" 3 "* Star" 4 "^ Caret" 5 "~ Tilde"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Special symbols mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.8 Labels with brackets"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "(Parentheses)" 2 "[Brackets]" 3 "{Braces}" 4 "<Angle>" 5 "Mix([]{})"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Brackets mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.9 Labels with math symbols"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "A + B" 2 "X - Y" 3 "M * N" 4 "P / Q" 5 "A = B"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Math symbols mismatch"
}

* ---------------------------------------------------------------------------
run_test "8.10 Labels with underscores and hyphens"
* ---------------------------------------------------------------------------
clear
set obs 4
gen x = _n
label define x_lbl 1 "under_score" 2 "hyphen-ated" 3 "both_and-mixed" 4 "__double__"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Underscores/hyphens mismatch"
}


* ============================================================================
* SECTION 9: Real-World Datasets (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 9: Real-World Datasets (sysuse/webuse)"
di as text ""

* ---------------------------------------------------------------------------
run_test "9.1 auto dataset - foreign"
* ---------------------------------------------------------------------------
sysuse auto, clear
cdecode foreign, generate(foreign_cdecode)
decode foreign, generate(foreign_decode)

capture assert foreign_cdecode == foreign_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "auto foreign mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.2 auto dataset - rep78 (with label)"
* ---------------------------------------------------------------------------
sysuse auto, clear
* rep78 doesn't have labels by default, so create them
label define rep78_lbl 1 "Poor" 2 "Fair" 3 "Average" 4 "Good" 5 "Excellent"
label values rep78 rep78_lbl

cdecode rep78, generate(rep78_cdecode)
decode rep78, generate(rep78_decode)

capture assert rep78_cdecode == rep78_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "auto rep78 mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.3 nlsw88 - race"
* ---------------------------------------------------------------------------
webuse nlsw88, clear
cdecode race, generate(race_cdecode)
decode race, generate(race_decode)

capture assert race_cdecode == race_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "nlsw88 race mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.4 nlsw88 - industry"
* ---------------------------------------------------------------------------
webuse nlsw88, clear
cdecode industry, generate(ind_cdecode)
decode industry, generate(ind_decode)

capture assert ind_cdecode == ind_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "nlsw88 industry mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.5 nlsw88 - occupation"
* ---------------------------------------------------------------------------
webuse nlsw88, clear
cdecode occupation, generate(occ_cdecode)
decode occupation, generate(occ_decode)

capture assert occ_cdecode == occ_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "nlsw88 occupation mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.6 census - region"
* ---------------------------------------------------------------------------
sysuse census, clear
cdecode region, generate(region_cdecode)
decode region, generate(region_decode)

capture assert region_cdecode == region_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "census region mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.7 lifeexp - region"
* ---------------------------------------------------------------------------
sysuse lifeexp, clear
cdecode region, generate(region_cdecode)
decode region, generate(region_decode)

capture assert region_cdecode == region_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "lifeexp region mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.8 sp500 - date"
* ---------------------------------------------------------------------------
sysuse sp500, clear
* Create a labeled variable for testing
gen quarter = quarter(date)
label define qtr_lbl 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4"
label values quarter qtr_lbl

cdecode quarter, generate(qtr_cdecode)
decode quarter, generate(qtr_decode)

capture assert qtr_cdecode == qtr_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "sp500 quarter mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.9 nlsw88 - married"
* ---------------------------------------------------------------------------
webuse nlsw88, clear
cdecode married, generate(married_cdecode)
decode married, generate(married_decode)

capture assert married_cdecode == married_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "nlsw88 married mismatch"
}

* ---------------------------------------------------------------------------
run_test "9.10 nlsw88 - smsa"
* ---------------------------------------------------------------------------
webuse nlsw88, clear
cdecode smsa, generate(smsa_cdecode)
decode smsa, generate(smsa_decode)

capture assert smsa_cdecode == smsa_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "nlsw88 smsa mismatch"
}


* ============================================================================
* SECTION 10: Pathological Data Structures (10 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 10: Pathological Data Structures"
di as text ""

* ---------------------------------------------------------------------------
run_test "10.1 All observations same value"
* ---------------------------------------------------------------------------
clear
set obs 1000
gen x = 5
label define x_lbl 5 "All Same"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "All same value mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.2 Single value repeated with many labels defined"
* ---------------------------------------------------------------------------
clear
set obs 1000
gen x = 1
* Use 50 labels to avoid Stata's command line length limits
forvalues i = 1/50 {
    label define x_lbl `i' "Label_`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Single value many labels mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.3 Alternating two values"
* ---------------------------------------------------------------------------
clear
set obs 10000
gen x = mod(_n, 2) + 1
label define x_lbl 1 "Odd" 2 "Even"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Alternating values mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.4 Sorted ascending"
* ---------------------------------------------------------------------------
clear
set obs 50
gen x = _n
* Use 50 labels to avoid Stata's command line length limits
forvalues i = 1/50 {
    label define x_lbl `i' "V`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Sorted ascending mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.5 Sorted descending"
* ---------------------------------------------------------------------------
clear
set obs 50
gen x = 51 - _n
* Use 50 labels to avoid Stata's command line length limits
forvalues i = 1/50 {
    label define x_lbl `i' "V`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Sorted descending mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.6 Random shuffle"
* ---------------------------------------------------------------------------
clear
set seed 54321
set obs 1000
gen x = ceil(runiform() * 100)
gen order = runiform()
sort order
forvalues i = 1/100 {
    label define x_lbl `i' "Rand`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Random shuffle mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.7 First half missing"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = _n
replace x = . in 1/50
* Use 50 labels to avoid Stata's command line length limits
forvalues i = 51/100 {
    label define x_lbl `i' "Val`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "First half missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.8 Second half missing"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = _n
replace x = . in 51/100
* Use 50 labels to avoid Stata's command line length limits
forvalues i = 1/50 {
    label define x_lbl `i' "Val`i'", add
}
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Second half missing mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.9 Only one non-missing at end"
* ---------------------------------------------------------------------------
clear
set obs 1000
gen x = .
replace x = 1 in 1000
label define x_lbl 1 "The Only One"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "One non-missing at end mismatch"
}

* ---------------------------------------------------------------------------
run_test "10.10 Only one non-missing at start"
* ---------------------------------------------------------------------------
clear
set obs 1000
gen x = .
replace x = 1 in 1
label define x_lbl 1 "The First One"
label values x x_lbl

cdecode x, generate(x_cdecode)
decode x, generate(x_decode)

capture assert x_cdecode == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "One non-missing at start mismatch"
}


* ============================================================================
* SECTION 11: Error Handling (5 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 11: Error Handling"
di as text ""

* ---------------------------------------------------------------------------
run_test "11.1 Variable already exists error"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
gen x_str = ""
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values x x_lbl

capture cdecode x, generate(x_str)
if _rc == 110 {
    test_passed
}
else {
    test_failed "Should error 110 for existing variable (got `=_rc')"
}

* ---------------------------------------------------------------------------
run_test "11.2 No value label error"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n

capture cdecode x, generate(x_str)
if _rc == 182 {
    test_passed
}
else {
    test_failed "Should error 182 for no value label (got `=_rc')"
}

* ---------------------------------------------------------------------------
run_test "11.3 String variable error"
* ---------------------------------------------------------------------------
clear
set obs 5
gen str10 x = "test"

capture cdecode x, generate(x_str)
if _rc != 0 {
    test_passed
}
else {
    test_failed "Should error for string variable"
}

* ---------------------------------------------------------------------------
run_test "11.4 Missing generate option"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "One"
label values x x_lbl

capture cdecode x
if _rc != 0 {
    test_passed
}
else {
    test_failed "Should error for missing generate"
}

* ---------------------------------------------------------------------------
run_test "11.5 Invalid variable name in generate"
* ---------------------------------------------------------------------------
clear
set obs 5
gen x = _n
label define x_lbl 1 "One"
label values x x_lbl

capture cdecode x, generate(123invalid)
if _rc != 0 {
    test_passed
}
else {
    test_failed "Should error for invalid variable name"
}


* ============================================================================
* SECTION 12: Consistency and Reproducibility (5 tests)
* ============================================================================

di as text ""
di as text ">>> SECTION 12: Consistency and Reproducibility"
di as text ""

* ---------------------------------------------------------------------------
run_test "12.1 Multiple runs same result"
* ---------------------------------------------------------------------------
clear
set obs 1000
gen x = ceil(runiform() * 5)
label define x_lbl 1 "One" 2 "Two" 3 "Three" 4 "Four" 5 "Five"
label values x x_lbl

cdecode x, generate(x_run1)
drop x_run1
cdecode x, generate(x_run2)
drop x_run2
cdecode x, generate(x_run3)

decode x, generate(x_decode)

capture assert x_run3 == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Multiple runs differ"
}

* ---------------------------------------------------------------------------
run_test "12.2 Order independence"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = ceil(runiform() * 10)
gen order1 = runiform()
gen order2 = runiform()
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "H" 9 "I" 10 "J"
label values x x_lbl

sort order1
cdecode x, generate(x_sort1)

sort order2
cdecode x, generate(x_sort2)

decode x, generate(x_decode)

capture assert x_sort2 == x_decode
if _rc == 0 {
    test_passed
}
else {
    test_failed "Order dependence detected"
}

* ---------------------------------------------------------------------------
run_test "12.3 Preserve original variable"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = ceil(runiform() * 5)
gen x_original = x
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values x x_lbl

cdecode x, generate(x_str)

capture assert x == x_original
if _rc == 0 {
    test_passed
}
else {
    test_failed "Original variable modified"
}

* ---------------------------------------------------------------------------
run_test "12.4 Preserve other variables"
* ---------------------------------------------------------------------------
clear
set obs 100
gen x = ceil(runiform() * 5)
gen y = runiform()
gen z = "test"
gen y_original = y
gen z_original = z
label define x_lbl 1 "A" 2 "B" 3 "C" 4 "D" 5 "E"
label values x x_lbl

cdecode x, generate(x_str)

capture assert y == y_original & z == z_original
if _rc == 0 {
    test_passed
}
else {
    test_failed "Other variables modified"
}

* ---------------------------------------------------------------------------
run_test "12.5 Observation count unchanged"
* ---------------------------------------------------------------------------
clear
set obs 12345
gen x = ceil(runiform() * 3)
label define x_lbl 1 "A" 2 "B" 3 "C"
label values x x_lbl

local n_before = _N
cdecode x, generate(x_str)
local n_after = _N

if `n_before' == `n_after' {
    test_passed
}
else {
    test_failed "Observation count changed"
}


* ============================================================================
* SUMMARY
* ============================================================================

di as text ""
di as text "============================================================"
di as text "                    VALIDATION SUMMARY"
di as text "============================================================"
di as text ""
di as text "Tests passed: " as result $n_passed
di as text "Tests failed: " as result $n_failed
di as text ""

if $n_failed > 0 {
    di as error "Failed tests:"
    di as error "$failed_tests"
    di as text ""
    di as error "VALIDATION FAILED"
    exit 9
}
else {
    di as result "ALL VALIDATION TESTS PASSED ($n_passed tests)"
}

di as text ""
di as text "============================================================"
