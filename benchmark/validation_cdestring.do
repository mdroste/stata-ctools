/*
 * validation_cdestring.do
 *
 * Validation tests for cdestring command
 * Tests all options and compares with Stata's destring
 */

clear all
set more off

* Set working directory and load ctools
cd "${CTOOLS_ROOT}"
if "${CTOOLS_ROOT}" == "" {
    local cwd "`c(pwd)'"
    cd "`cwd'"
}
adopath + "build"

di as text ""
di as text "=========================================="
di as text "   cdestring Validation Tests"
di as text "=========================================="
di as text ""

* ============================================================================
* Test 1: Basic numeric strings
* ============================================================================
di as text "{hline 50}"
di as text "Test 1: Basic numeric strings"
di as text "{hline 50}"

clear
set obs 10
gen str10 num_str = string(_n * 100)

* Test with generate
cdestring num_str, generate(num_val)
assert num_val[1] == 100
assert num_val[5] == 500
assert num_val[10] == 1000
di as result "  PASSED: Basic numeric conversion with generate()"

* Test with replace
gen str10 num_str2 = string(_n * 10)
cdestring num_str2, replace
assert num_str2[1] == 10
assert num_str2[10] == 100
di as result "  PASSED: Basic numeric conversion with replace"

* ============================================================================
* Test 2: Decimal numbers
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 2: Decimal numbers"
di as text "{hline 50}"

clear
set obs 5
gen str15 dec_str = ""
replace dec_str = "123.456" in 1
replace dec_str = "-0.5" in 2
replace dec_str = ".25" in 3
replace dec_str = "0.0" in 4
replace dec_str = "999.999" in 5

cdestring dec_str, generate(dec_val)
assert abs(dec_val[1] - 123.456) < 0.0001
assert abs(dec_val[2] - (-0.5)) < 0.0001
assert abs(dec_val[3] - 0.25) < 0.0001
assert dec_val[4] == 0
assert abs(dec_val[5] - 999.999) < 0.0001
di as result "  PASSED: Decimal number conversion"

* ============================================================================
* Test 3: Scientific notation
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 3: Scientific notation"
di as text "{hline 50}"

clear
set obs 4
gen str15 sci_str = ""
replace sci_str = "1.23e4" in 1
replace sci_str = "1.5E-3" in 2
replace sci_str = "-2e10" in 3
replace sci_str = "5e0" in 4

cdestring sci_str, generate(sci_val)
assert abs(sci_val[1] - 12300) < 1
assert abs(sci_val[2] - 0.0015) < 0.00001
assert abs(sci_val[3] - (-2e10)) < 1e6
assert sci_val[4] == 5
di as result "  PASSED: Scientific notation conversion"

* ============================================================================
* Test 4: Missing values
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 4: Missing values"
di as text "{hline 50}"

clear
set obs 6
gen str10 miss_str = ""
replace miss_str = "" in 1
replace miss_str = "." in 2
replace miss_str = "NA" in 3
replace miss_str = "NaN" in 4
replace miss_str = "   " in 5  // whitespace
replace miss_str = "123" in 6

cdestring miss_str, generate(miss_val)
assert miss_val[1] == .
assert miss_val[2] == .
assert miss_val[3] == .
assert miss_val[4] == .
assert miss_val[5] == .
assert miss_val[6] == 123
di as result "  PASSED: Missing value handling"

* ============================================================================
* Test 5: ignore() option
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 5: ignore() option"
di as text "{hline 50}"

clear
set obs 4
gen str15 curr_str = ""
replace curr_str = "$1,234.56" in 1
replace curr_str = "$500" in 2
replace curr_str = "1,000,000" in 3
replace curr_str = "€99.99" in 4

cdestring curr_str, generate(curr_val) ignore("$,€")
assert abs(curr_val[1] - 1234.56) < 0.01
assert curr_val[2] == 500
assert curr_val[3] == 1000000
assert abs(curr_val[4] - 99.99) < 0.01
di as result "  PASSED: ignore() option"

* ============================================================================
* Test 6: percent option
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 6: percent option"
di as text "{hline 50}"

clear
set obs 4
gen str10 pct_str = ""
replace pct_str = "50%" in 1
replace pct_str = "100%" in 2
replace pct_str = "0.5%" in 3
replace pct_str = "25" in 4  // no percent sign

cdestring pct_str, generate(pct_val) percent
assert abs(pct_val[1] - 0.50) < 0.001
assert abs(pct_val[2] - 1.00) < 0.001
assert abs(pct_val[3] - 0.005) < 0.0001
assert pct_val[4] == 25  // no % sign, no division
di as result "  PASSED: percent option"

* ============================================================================
* Test 7: dpcomma option (European format)
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 7: dpcomma option (European format)"
di as text "{hline 50}"

clear
set obs 3
gen str15 euro_str = ""
replace euro_str = "1.234,56" in 1
replace euro_str = "999,99" in 2
replace euro_str = "1.000.000,00" in 3

cdestring euro_str, generate(euro_val) dpcomma
assert abs(euro_val[1] - 1234.56) < 0.01
assert abs(euro_val[2] - 999.99) < 0.01
assert abs(euro_val[3] - 1000000.00) < 0.01
di as result "  PASSED: dpcomma option"

* ============================================================================
* Test 8: force option with non-numeric data
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 8: force option"
di as text "{hline 50}"

clear
set obs 5
gen str15 mixed_str = ""
replace mixed_str = "123" in 1
replace mixed_str = "abc" in 2
replace mixed_str = "45.6" in 3
replace mixed_str = "N/A" in 4
replace mixed_str = "789" in 5

cdestring mixed_str, generate(mixed_val) force
assert mixed_val[1] == 123
assert mixed_val[2] == .
assert abs(mixed_val[3] - 45.6) < 0.01
assert mixed_val[4] == .
assert mixed_val[5] == 789
di as result "  PASSED: force option"

* ============================================================================
* Test 9: float option
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 9: float option"
di as text "{hline 50}"

clear
set obs 3
gen str10 flt_str = "123.456"

cdestring flt_str, generate(flt_val) float
assert "`: type flt_val'" == "float"
di as result "  PASSED: float option"

* ============================================================================
* Test 10: Multiple variables
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 10: Multiple variables"
di as text "{hline 50}"

clear
set obs 5
gen str10 a_str = string(_n)
gen str10 b_str = string(_n * 10)
gen str10 c_str = string(_n * 100)

cdestring a_str b_str c_str, generate(a b c)
assert a[3] == 3
assert b[3] == 30
assert c[3] == 300
di as result "  PASSED: Multiple variable conversion"

* ============================================================================
* Test 11: if/in restrictions
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 11: if/in restrictions"
di as text "{hline 50}"

clear
set obs 10
gen str10 cond_str = string(_n)
gen flag = mod(_n, 2) == 0

cdestring cond_str if flag, generate(cond_val)
assert cond_val[1] == .  // odd, not converted
assert cond_val[2] == 2  // even, converted
assert cond_val[3] == .  // odd
assert cond_val[4] == 4  // even
di as result "  PASSED: if condition"

clear
set obs 10
gen str10 cond_str = string(_n)

cdestring cond_str in 3/7, generate(cond_val)
assert cond_val[1] == .
assert cond_val[2] == .
assert cond_val[3] == 3
assert cond_val[5] == 5
assert cond_val[7] == 7
assert cond_val[8] == .
di as result "  PASSED: in range"

* ============================================================================
* Test 12: Combined options
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 12: Combined options"
di as text "{hline 50}"

clear
set obs 3
gen str20 combo_str = ""
replace combo_str = "$1,234.56%" in 1
replace combo_str = "€500,00%" in 2
replace combo_str = "invalid" in 3

* Note: dpcomma and ignore with force, percent
* This tests ignore("$€,") with dpcomma and percent
* But dpcomma changes comma to decimal, which conflicts with ignore(",")
* So we use just ignore("$€%") with dpcomma

clear
set obs 2
gen str20 combo_str = ""
replace combo_str = "€1.234,56" in 1
replace combo_str = "€999,00" in 2

cdestring combo_str, generate(combo_val) ignore("€") dpcomma
assert abs(combo_val[1] - 1234.56) < 0.01
assert abs(combo_val[2] - 999.00) < 0.01
di as result "  PASSED: Combined options"

* ============================================================================
* Test 13: Performance comparison (small)
* ============================================================================
di as text ""
di as text "{hline 50}"
di as text "Test 13: Performance test with verbose"
di as text "{hline 50}"

clear
set obs 100000
gen str12 perf_str = string(runiform() * 10000, "%12.4f")

timer clear 1
timer on 1
cdestring perf_str, generate(perf_val) verbose
timer off 1

timer list 1
local cdestring_time = r(t1)
di as result "  cdestring time: `cdestring_time' seconds"

* Compare with native destring
clear
set obs 100000
gen str12 perf_str = string(runiform() * 10000, "%12.4f")

timer clear 2
timer on 2
destring perf_str, generate(perf_val2)
timer off 2

timer list 2
local destring_time = r(t2)
di as result "  destring time:  `destring_time' seconds"

local speedup = `destring_time' / `cdestring_time'
di as result "  Speedup: " %5.2f `speedup' "x"

* ============================================================================
* Summary
* ============================================================================
di as text ""
di as text "=========================================="
di as text "   All cdestring tests PASSED"
di as text "=========================================="
di as text ""
