/*******************************************************************************
 * validate_cdestring.do
 *
 * Comprehensive validation tests for cdestring vs native Stata destring
 * Tests all options, data types, and edge cases
 *
 * cdestring: C-accelerated string to numeric conversion
 * Syntax: cdestring varlist [if] [in], {generate(newvarlist) | replace}
 *         [ignore("chars") force float percent dpcomma verbose threads(#)]
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text ""
noi di as text "======================================================================"
noi di as text "              CDESTRING VALIDATION TEST SUITE"
noi di as text "======================================================================"

/*******************************************************************************
 * Helper: Compare cdestring vs destring
 ******************************************************************************/
capture program drop benchmark_destring
program define benchmark_destring
    syntax varlist(string), testname(string) [GENerate(string) replace ///
        IGnore(string) force float percent dpcomma if2(string) in2(string)]

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build common options
    local common_opts ""
    if `"`ignore'"' != "" local common_opts `"`common_opts' ignore("`ignore'")"'
    if "`force'" != "" local common_opts "`common_opts' force"
    if "`float'" != "" local common_opts "`common_opts' float"
    if "`percent'" != "" local common_opts "`common_opts' percent"
    if "`dpcomma'" != "" local common_opts "`common_opts' dpcomma"

    * Build generate names with _c and _s suffixes
    if "`generate'" != "" {
        local gen_c ""
        local gen_s ""
        foreach v of local generate {
            local gen_c "`gen_c' `v'_c"
            local gen_s "`gen_s' `v'_s"
        }
        local opts_c "generate(`gen_c') `common_opts'"
        local opts_s "generate(`gen_s') `common_opts'"
    }
    else if "`replace'" != "" {
        * For replace, we need to work with copies
        local opts_c "replace `common_opts'"
        local opts_s "replace `common_opts'"
    }
    else {
        noi test_fail "`testname'" "must specify generate or replace"
        exit
    }

    * Handle replace option specially - need to create copies
    if "`replace'" != "" {
        * Create copies of original variables
        local copy_c ""
        local copy_s ""
        foreach v of local varlist {
            capture drop `v'_c `v'_s
            clonevar `v'_c = `v'
            clonevar `v'_s = `v'
        }

        * Run cdestring on copies
        local varlist_c ""
        foreach v of local varlist {
            local varlist_c "`varlist_c' `v'_c"
        }
        capture cdestring `varlist_c' `ifin', replace `common_opts'
        local rc_c = _rc

        * Run destring on copies
        local varlist_s ""
        foreach v of local varlist {
            local varlist_s "`varlist_s' `v'_s"
        }
        capture destring `varlist_s' `ifin', replace `common_opts'
        local rc_s = _rc

        * Check both succeeded or both failed
        if `rc_c' != `rc_s' {
            noi test_fail "`testname'" "cdestring rc=`rc_c', destring rc=`rc_s'"
            exit
        }

        if `rc_c' != 0 {
            noi test_pass "`testname' (both error as expected)"
            exit
        }

        * Compare results (use tolerance for floating-point comparison)
        local ndiff_total = 0
        foreach v of local varlist {
            quietly count if abs(`v'_c - `v'_s) > 1e-10 & !missing(`v'_c) & !missing(`v'_s)
            local ndiff_total = `ndiff_total' + r(N)
            * Also check missing patterns match
            quietly count if missing(`v'_c) != missing(`v'_s)
            local ndiff_total = `ndiff_total' + r(N)
        }

        if `ndiff_total' > 0 {
            noi test_fail "`testname'" "`ndiff_total' values differ"
        }
        else {
            noi test_pass "`testname'"
        }

        * Cleanup
        foreach v of local varlist {
            capture drop `v'_c `v'_s
        }
    }
    else {
        * Generate option - simpler case

        * Drop any existing test variables
        foreach v of local gen_c {
            capture drop `v'
        }
        foreach v of local gen_s {
            capture drop `v'
        }

        * Run cdestring
        capture cdestring `varlist' `ifin', `opts_c'
        local rc_c = _rc

        * Run destring
        capture destring `varlist' `ifin', `opts_s'
        local rc_s = _rc

        * Check both succeeded or both failed
        if `rc_c' != `rc_s' {
            noi test_fail "`testname'" "cdestring rc=`rc_c', destring rc=`rc_s'"
            * Cleanup
            foreach v of local gen_c {
                capture drop `v'
            }
            foreach v of local gen_s {
                capture drop `v'
            }
            exit
        }

        if `rc_c' != 0 {
            noi test_pass "`testname' (both error as expected)"
            exit
        }

        * Compare numeric values
        local ndiff_total = 0
        local i = 1
        foreach v of local generate {
            local vc = "`v'_c"
            local vs = "`v'_s"
            * Compare non-missing values with tolerance
            quietly count if abs(`vc' - `vs') > 1e-10 & !missing(`vc') & !missing(`vs')
            local ndiff_total = `ndiff_total' + r(N)
            * Check missing patterns match
            quietly count if missing(`vc') != missing(`vs')
            local ndiff_total = `ndiff_total' + r(N)
            local ++i
        }

        if `ndiff_total' > 0 {
            noi test_fail "`testname'" "`ndiff_total' values differ"
        }
        else {
            noi test_pass "`testname'"
        }

        * Cleanup
        foreach v of local gen_c {
            capture drop `v'
        }
        foreach v of local gen_s {
            capture drop `v'
        }
    }
end

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
noi print_section "Plugin Check"

clear
set obs 10
gen str10 numstr = string(_n)
capture noisily cdestring numstr, generate(numval)
if _rc != 0 {
    noi test_fail "cdestring plugin load" "plugin returned error `=_rc'"
    noi print_summary "cdestring"
    exit 1
}
noi test_pass "cdestring plugin loads and runs"
drop numval

/*******************************************************************************
 * SECTION 2: Basic functionality (10 tests)
 ******************************************************************************/
noi print_section "Basic Functionality"

* Test 2.1: Basic conversion creates numeric variable
clear
set obs 10
gen str10 x = string(_n)
cdestring x, generate(x_num)
capture confirm numeric variable x_num
if _rc == 0 {
    noi test_pass "basic conversion creates numeric variable"
}
else {
    noi test_fail "basic conversion" "variable not numeric"
}

* Test 2.2: Compare with destring - simple integers
clear
set obs 100
gen str10 x = string(_n)
noi benchmark_destring x, testname("vs destring: simple integers") generate(x)

* Test 2.3: Values correctly converted
clear
set obs 5
gen str10 x = string(_n * 10)
cdestring x, generate(x_num)
local pass = 1
forvalues i = 1/5 {
    local expected = `i' * 10
    local actual = x_num[`i']
    if `expected' != `actual' {
        local pass = 0
    }
}
if `pass' {
    noi test_pass "values correctly converted"
}
else {
    noi test_fail "values converted" "mismatch found"
}

* Test 2.4: Single observation
clear
set obs 1
gen str10 x = "42"
cdestring x, generate(x_num)
if x_num[1] == 42 {
    noi test_pass "single observation"
}
else {
    noi test_fail "single observation" "got `=x_num[1]'"
}

* Test 2.5: verbose option
clear
set obs 10
gen str10 x = string(_n)
capture cdestring x, generate(x_num) verbose
if _rc == 0 {
    noi test_pass "verbose option"
}
else {
    noi test_fail "verbose option" "rc=`=_rc'"
}

* Test 2.6: threads option
clear
set obs 100
gen str10 x = string(_n)
capture cdestring x, generate(x_num) threads(2)
if _rc == 0 {
    noi test_pass "threads(2) option"
}
else {
    noi test_fail "threads option" "rc=`=_rc'"
}

* Test 2.7: Negative numbers
clear
set obs 10
gen str10 x = string(_n - 5)
cdestring x, generate(x_num)
if x_num[1] == -4 & x_num[5] == 0 & x_num[10] == 5 {
    noi test_pass "negative numbers"
}
else {
    noi test_fail "negative numbers" "wrong values"
}

* Test 2.8: Compare negative numbers with destring
clear
set obs 50
gen str10 x = string(_n - 25)
noi benchmark_destring x, testname("vs destring: negative integers") generate(x)

* Test 2.9: Decimals
clear
set obs 5
gen str20 x = ""
replace x = "1.5" in 1
replace x = "2.25" in 2
replace x = "3.125" in 3
replace x = "0.001" in 4
replace x = "-99.99" in 5
cdestring x, generate(x_num)
local pass = (abs(x_num[1] - 1.5) < 1e-10) & (abs(x_num[3] - 3.125) < 1e-10)
if `pass' {
    noi test_pass "decimal numbers"
}
else {
    noi test_fail "decimal numbers" "wrong values"
}

* Test 2.10: Compare decimals with destring
clear
set obs 100
gen str20 x = string(_n / 7, "%20.10f")
noi benchmark_destring x, testname("vs destring: decimal numbers") generate(x)

/*******************************************************************************
 * SECTION 3: Replace option (10 tests)
 ******************************************************************************/
noi print_section "Replace Option"

* Test 3.1: Replace converts in place
clear
set obs 10
gen str10 x = string(_n)
cdestring x, replace
capture confirm numeric variable x
if _rc == 0 {
    noi test_pass "replace converts in place"
}
else {
    noi test_fail "replace option" "variable still string"
}

* Test 3.2: Replace preserves values
clear
set obs 5
gen str10 x = string(_n * 100)
cdestring x, replace
local pass = (x[1] == 100) & (x[3] == 300) & (x[5] == 500)
if `pass' {
    noi test_pass "replace preserves values"
}
else {
    noi test_fail "replace values" "wrong conversion"
}

* Test 3.3: Compare replace with destring
clear
set obs 50
gen str10 x = string(_n)
noi benchmark_destring x, testname("vs destring: replace option") replace

* Test 3.4: Replace with negative numbers
clear
set obs 20
gen str10 x = string(_n - 10)
noi benchmark_destring x, testname("vs destring: replace negative") replace

* Test 3.5: Replace with decimals
clear
set obs 30
gen str20 x = string(_n / 3, "%10.6f")
noi benchmark_destring x, testname("vs destring: replace decimals") replace

* Test 3.6: Replace multiple variables
clear
set obs 20
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
gen str10 c = string(_n * 3)
noi benchmark_destring a b c, testname("vs destring: replace multi") replace

* Test 3.7: Replace with if condition
clear
set obs 50
gen str10 x = string(_n)
gen flag = mod(_n, 2)
cdestring x if flag == 1, replace
* Should only convert odd observations
count if missing(x) & flag == 0
local n_missing = r(N)
* Even obs should still be string (but replace makes them missing if not in sample)
if `n_missing' == 25 {
    noi test_pass "replace with if condition"
}
else {
    noi test_fail "replace if" "unexpected missing count"
}

* Test 3.8: Replace with in range
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 1/50, replace
* Only first 50 should be converted
capture confirm numeric variable x
if _rc == 0 {
    count if !missing(x) in 1/50
    if r(N) == 50 {
        noi test_pass "replace with in range"
    }
    else {
        noi test_fail "replace in" "wrong count"
    }
}
else {
    noi test_fail "replace in" "not converted"
}

* Test 3.9: Replace float option
clear
set obs 10
gen str20 x = string(_n + 0.5)
cdestring x, replace float
local vtype : type x
if "`vtype'" == "float" {
    noi test_pass "replace with float option"
}
else {
    noi test_fail "replace float" "type is `vtype'"
}

* Test 3.10: Compare replace float with destring
clear
set obs 50
gen str20 x = string(_n / 2)
noi benchmark_destring x, testname("vs destring: replace float") replace float

/*******************************************************************************
 * SECTION 4: Ignore option (10 tests)
 ******************************************************************************/
noi print_section "Ignore Option"

* Test 4.1: Ignore single character
clear
set obs 5
gen str20 x = "$" + string(_n * 100)
cdestring x, generate(x_num) ignore("$")
if x_num[1] == 100 & x_num[5] == 500 {
    noi test_pass "ignore single character ($)"
}
else {
    noi test_fail "ignore single" "wrong values"
}

* Test 4.2: Compare ignore with destring
clear
set obs 50
gen str20 x = "$" + string(_n)
noi benchmark_destring x, testname("vs destring: ignore $") generate(x) ignore("$")

* Test 4.3: Ignore multiple characters
clear
set obs 5
gen str30 x = "$" + string(_n * 1000) + ",00"
replace x = subinstr(x, "000", ",000", 1) if _n >= 1
cdestring x, generate(x_num) ignore("$,")
* Should handle currency-like formatting
capture confirm numeric variable x_num
if _rc == 0 {
    noi test_pass "ignore multiple characters"
}
else {
    noi test_fail "ignore multiple" "conversion failed"
}

* Test 4.4: Ignore comma in numbers
clear
set obs 5
gen str20 x = ""
replace x = "1,000" in 1
replace x = "2,500" in 2
replace x = "10,000" in 3
replace x = "100,000" in 4
replace x = "1,234,567" in 5
noi benchmark_destring x, testname("vs destring: ignore comma") generate(x) ignore(",")

* Test 4.5: Ignore space
clear
set obs 5
gen str20 x = ""
replace x = "1 000" in 1
replace x = "2 500" in 2
replace x = "10 000" in 3
replace x = "100 000" in 4
replace x = "1 234 567" in 5
noi benchmark_destring x, testname("vs destring: ignore space") generate(x) ignore(" ")

* Test 4.6: Ignore letters (force may be needed)
clear
set obs 5
gen str20 x = ""
replace x = "100kg" in 1
replace x = "200kg" in 2
replace x = "300kg" in 3
replace x = "400kg" in 4
replace x = "500kg" in 5
noi benchmark_destring x, testname("vs destring: ignore kg") generate(x) ignore("kg")

* Test 4.7: Ignore percent sign
clear
set obs 5
gen str10 x = string(_n * 10) + "%"
noi benchmark_destring x, testname("vs destring: ignore %") generate(x) ignore("%")

* Test 4.8: Ignore with replace
clear
set obs 50
gen str20 x = "$" + string(_n * 100)
noi benchmark_destring x, testname("vs destring: ignore $ replace") replace ignore("$")

* Test 4.9: Ignore parentheses (negative numbers)
clear
set obs 5
gen str20 x = ""
replace x = "(100)" in 1
replace x = "(200)" in 2
replace x = "300" in 3
replace x = "(400)" in 4
replace x = "500" in 5
* Note: destring doesn't interpret () as negative, just strips them
noi benchmark_destring x, testname("vs destring: ignore parens") generate(x) ignore("()")

* Test 4.10: Empty ignore string (no effect)
clear
set obs 20
gen str10 x = string(_n)
noi benchmark_destring x, testname("vs destring: empty ignore") generate(x) ignore("")

/*******************************************************************************
 * SECTION 5: Force option (10 tests)
 ******************************************************************************/
noi print_section "Force Option"

* Test 5.1: Force converts non-numeric to missing
clear
set obs 5
gen str20 x = ""
replace x = "100" in 1
replace x = "abc" in 2
replace x = "300" in 3
replace x = "def" in 4
replace x = "500" in 5
cdestring x, generate(x_num) force
local pass = (x_num[1] == 100) & missing(x_num[2]) & (x_num[3] == 300) & missing(x_num[4])
if `pass' {
    noi test_pass "force converts non-numeric to missing"
}
else {
    noi test_fail "force option" "wrong values"
}

* Test 5.2: Compare force with destring
clear
set obs 50
gen str20 x = cond(mod(_n, 3) == 0, "abc", string(_n))
noi benchmark_destring x, testname("vs destring: force option") generate(x) force

* Test 5.3: Force with all non-numeric
clear
set obs 10
gen str10 x = char(64 + _n)  // A, B, C, ...
cdestring x, generate(x_num) force
count if missing(x_num)
if r(N) == 10 {
    noi test_pass "force with all non-numeric"
}
else {
    noi test_fail "force all non-numeric" "expected 10 missing"
}

* Test 5.4: Force with mixed content
clear
set obs 20
gen str20 x = ""
replace x = string(_n) if mod(_n, 2) == 1
replace x = "text" + string(_n) if mod(_n, 2) == 0
noi benchmark_destring x, testname("vs destring: force mixed") generate(x) force

* Test 5.5: Force with leading/trailing spaces
clear
set obs 5
gen str20 x = ""
replace x = " 100" in 1
replace x = "200 " in 2
replace x = " 300 " in 3
replace x = "  400  " in 4
replace x = "500" in 5
noi benchmark_destring x, testname("vs destring: force spaces") generate(x) force

* Test 5.6: Force with replace
clear
set obs 30
gen str20 x = cond(mod(_n, 4) == 0, "NA", string(_n))
noi benchmark_destring x, testname("vs destring: force replace") replace force

* Test 5.7: Force with ignore
clear
set obs 20
gen str20 x = "$" + cond(mod(_n, 5) == 0, "NA", string(_n))
noi benchmark_destring x, testname("vs destring: force + ignore") generate(x) force ignore("$")

* Test 5.8: Force preserves valid numbers
clear
set obs 100
gen str10 x = string(_n)
cdestring x, generate(x_num) force
count if x_num == real(x)
if r(N) == 100 {
    noi test_pass "force preserves valid numbers"
}
else {
    noi test_fail "force valid" "some values changed"
}

* Test 5.9: Force with empty strings
clear
set obs 10
gen str10 x = cond(mod(_n, 2) == 0, "", string(_n))
noi benchmark_destring x, testname("vs destring: force empty") generate(x) force

* Test 5.10: Force with special characters
clear
set obs 10
gen str20 x = ""
replace x = "100" in 1
replace x = "200!" in 2
replace x = "300" in 3
replace x = "@400" in 4
replace x = "500" in 5
replace x = "600#" in 6
replace x = "700" in 7
replace x = "800*" in 8
replace x = "900" in 9
replace x = "1000" in 10
noi benchmark_destring x, testname("vs destring: force special chars") generate(x) force

/*******************************************************************************
 * SECTION 6: Percent option (10 tests)
 ******************************************************************************/
noi print_section "Percent Option"

* Test 6.1: Percent strips % and divides by 100
clear
set obs 5
gen str10 x = string(_n * 10) + "%"
cdestring x, generate(x_num) percent
local pass = (abs(x_num[1] - 0.1) < 1e-10) & (abs(x_num[5] - 0.5) < 1e-10)
if `pass' {
    noi test_pass "percent option basic"
}
else {
    noi test_fail "percent option" "wrong values"
}

* Test 6.2: Compare percent with destring
clear
set obs 50
gen str10 x = string(_n) + "%"
noi benchmark_destring x, testname("vs destring: percent option") generate(x) percent

* Test 6.3: Percent with decimals
clear
set obs 5
gen str10 x = ""
replace x = "10.5%" in 1
replace x = "25.25%" in 2
replace x = "50.125%" in 3
replace x = "75.0%" in 4
replace x = "99.99%" in 5
noi benchmark_destring x, testname("vs destring: percent decimal") generate(x) percent

* Test 6.4: Percent with replace
clear
set obs 30
gen str10 x = string(_n * 3) + "%"
noi benchmark_destring x, testname("vs destring: percent replace") replace percent

* Test 6.5: Percent without % sign (no transformation)
clear
set obs 10
gen str10 x = string(_n * 10)
cdestring x, generate(x_num) percent
* Values without % should not be divided by 100
if x_num[1] == 10 & x_num[10] == 100 {
    noi test_pass "percent without % sign"
}
else {
    noi test_fail "percent no %" "values changed unexpectedly"
}

* Test 6.6: Compare percent without % with destring
clear
set obs 50
gen str10 x = string(_n)
noi benchmark_destring x, testname("vs destring: percent no %") generate(x) percent

* Test 6.7: Percent with negative values
clear
set obs 5
gen str10 x = ""
replace x = "-10%" in 1
replace x = "-25%" in 2
replace x = "0%" in 3
replace x = "50%" in 4
replace x = "-100%" in 5
noi benchmark_destring x, testname("vs destring: percent negative") generate(x) percent

* Test 6.8: Percent with force
clear
set obs 10
gen str10 x = cond(mod(_n, 3) == 0, "NA%", string(_n * 10) + "%")
noi benchmark_destring x, testname("vs destring: percent force") generate(x) percent force

* Test 6.9: 100% value
clear
set obs 5
gen str10 x = "100%"
cdestring x, generate(x_num) percent
if abs(x_num[1] - 1.0) < 1e-10 {
    noi test_pass "100% converts to 1.0"
}
else {
    noi test_fail "100% value" "got `=x_num[1]'"
}

* Test 6.10: Percent with small values
clear
set obs 5
gen str10 x = ""
replace x = "0.1%" in 1
replace x = "0.01%" in 2
replace x = "0.001%" in 3
replace x = "1%" in 4
replace x = "0.5%" in 5
noi benchmark_destring x, testname("vs destring: percent small") generate(x) percent

/*******************************************************************************
 * SECTION 7: dpcomma option (10 tests)
 ******************************************************************************/
noi print_section "dpcomma Option"

* Test 7.1: dpcomma uses comma as decimal point
clear
set obs 5
gen str10 x = ""
replace x = "1,5" in 1
replace x = "2,25" in 2
replace x = "3,125" in 3
replace x = "10,0" in 4
replace x = "0,5" in 5
cdestring x, generate(x_num) dpcomma
local pass = (abs(x_num[1] - 1.5) < 1e-10) & (abs(x_num[3] - 3.125) < 1e-10)
if `pass' {
    noi test_pass "dpcomma basic"
}
else {
    noi test_fail "dpcomma option" "wrong values"
}

* Test 7.2: Compare dpcomma with destring
clear
set obs 50
gen str20 x = subinstr(string(_n / 3, "%10.3f"), ".", ",", 1)
noi benchmark_destring x, testname("vs destring: dpcomma option") generate(x) dpcomma

* Test 7.3: dpcomma with negative numbers
clear
set obs 5
gen str10 x = ""
replace x = "-1,5" in 1
replace x = "-2,25" in 2
replace x = "0,0" in 3
replace x = "1,5" in 4
replace x = "-0,5" in 5
noi benchmark_destring x, testname("vs destring: dpcomma negative") generate(x) dpcomma

* Test 7.4: dpcomma with replace
clear
set obs 30
gen str20 x = subinstr(string(_n / 7, "%10.4f"), ".", ",", 1)
noi benchmark_destring x, testname("vs destring: dpcomma replace") replace dpcomma

* Test 7.5: dpcomma with integers (no comma)
clear
set obs 20
gen str10 x = string(_n)
cdestring x, generate(x_num) dpcomma
if x_num[1] == 1 & x_num[20] == 20 {
    noi test_pass "dpcomma with integers"
}
else {
    noi test_fail "dpcomma int" "wrong values"
}

* Test 7.6: dpcomma with force
clear
set obs 10
gen str20 x = cond(mod(_n, 3) == 0, "abc", subinstr(string(_n / 2, "%10.2f"), ".", ",", 1))
noi benchmark_destring x, testname("vs destring: dpcomma force") generate(x) dpcomma force

* Test 7.7: dpcomma with ignore (period as thousand separator)
clear
set obs 5
gen str20 x = ""
replace x = "1.000,50" in 1
replace x = "2.500,25" in 2
replace x = "10.000,00" in 3
replace x = "100.000,99" in 4
replace x = "1.234.567,89" in 5
noi benchmark_destring x, testname("vs destring: dpcomma + ignore .") generate(x) dpcomma ignore(".")

* Test 7.8: dpcomma zero
clear
set obs 3
gen str10 x = ""
replace x = "0,0" in 1
replace x = "0,00" in 2
replace x = "0,000" in 3
cdestring x, generate(x_num) dpcomma
if x_num[1] == 0 & x_num[2] == 0 & x_num[3] == 0 {
    noi test_pass "dpcomma zero values"
}
else {
    noi test_fail "dpcomma zero" "not zero"
}

* Test 7.9: dpcomma with percent
clear
set obs 5
gen str10 x = ""
replace x = "10,5%" in 1
replace x = "25,25%" in 2
replace x = "50,0%" in 3
replace x = "75,75%" in 4
replace x = "99,99%" in 5
noi benchmark_destring x, testname("vs destring: dpcomma percent") generate(x) dpcomma percent

* Test 7.10: Compare standard decimal with dpcomma
clear
set obs 30
* Standard decimal
gen str20 x_dot = string(_n / 7, "%10.4f")
* Comma decimal
gen str20 x_comma = subinstr(x_dot, ".", ",", 1)
destring x_dot, generate(x_dot_num)
cdestring x_comma, generate(x_comma_num) dpcomma
count if abs(x_dot_num - x_comma_num) > 1e-10
if r(N) == 0 {
    noi test_pass "dpcomma equals standard decimal"
}
else {
    noi test_fail "dpcomma vs dot" "values differ"
}

/*******************************************************************************
 * SECTION 8: if/in conditions (10 tests)
 ******************************************************************************/
noi print_section "if/in Conditions"

* Test 8.1: if condition subset
clear
set obs 100
gen str10 x = string(_n)
gen flag = (_n > 50)
cdestring x if flag == 1, generate(x_num)
count if !missing(x_num) & flag == 1
local n_conv = r(N)
count if missing(x_num) & flag == 0
local n_miss = r(N)
if `n_conv' == 50 & `n_miss' == 50 {
    noi test_pass "if condition subset"
}
else {
    noi test_fail "if condition" "wrong pattern"
}

* Test 8.2: in range
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 1/30, generate(x_num)
count if !missing(x_num) in 1/30
local in_range = r(N)
count if missing(x_num) in 31/100
local out_range = r(N)
if `in_range' == 30 & `out_range' == 70 {
    noi test_pass "in range"
}
else {
    noi test_fail "in range" "wrong counts"
}

* Test 8.3: if condition works correctly (destring doesn't support if/in)
clear
set obs 100
gen str10 x = string(_n)
gen value = runiform()
set seed 12345
replace value = runiform()
cdestring x if value > 0.5, generate(x_num)
* Count how many observations match the condition
count if value > 0.5
local n_expected = r(N)
count if !missing(x_num)
local n_converted = r(N)
if `n_converted' == `n_expected' {
    noi test_pass "if condition filters observations"
}
else {
    noi test_fail "if condition" "expected `n_expected' converted, got `n_converted'"
}

* Test 8.4: in range works correctly
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 20/80, generate(x_num)
count if !missing(x_num) in 20/80
local in_range = r(N)
count if missing(x_num) in 1/19
local miss_before = r(N)
count if missing(x_num) in 81/100
local miss_after = r(N)
if `in_range' == 61 & `miss_before' == 19 & `miss_after' == 20 {
    noi test_pass "in range works correctly"
}
else {
    noi test_fail "in range" "in=`in_range' miss_before=`miss_before' miss_after=`miss_after'"
}

* Test 8.5: Combined if and in
clear
set obs 100
gen str10 x = string(_n)
gen flag = mod(_n, 3) == 0
cdestring x if flag == 1 in 1/50, generate(x_num)
* Should convert only obs 3,6,9,...,48 (multiples of 3 up to 50) = 16 values
count if !missing(x_num)
local n_converted = r(N)
if `n_converted' == 16 {
    noi test_pass "combined if and in"
}
else {
    noi test_fail "if+in" "expected 16, got `n_converted'"
}

* Test 8.6: Single observation in
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 50/50, generate(x_num)
count if !missing(x_num)
if r(N) == 1 & x_num[50] == 50 {
    noi test_pass "single observation in"
}
else {
    noi test_fail "single in" "wrong result"
}

* Test 8.7: if that matches no observations
clear
set obs 50
gen str10 x = string(_n)
gen value = _n
cdestring x if value > 1000, generate(x_num)
count if !missing(x_num)
if r(N) == 0 {
    noi test_pass "if matches no observations"
}
else {
    noi test_fail "if none" "some values converted"
}

* Test 8.8: if with string comparison on another variable
clear
set obs 50
gen str10 x = string(_n)
gen str10 group = cond(mod(_n, 2) == 0, "even", "odd")
cdestring x if group == "even", generate(x_num)
count if !missing(x_num)
if r(N) == 25 {
    noi test_pass "if with string comparison"
}
else {
    noi test_fail "if string" "wrong count"
}

* Test 8.9: in first observations
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 1/10, generate(x_num)
count if !missing(x_num)
if r(N) == 10 {
    noi test_pass "in first 10"
}
else {
    noi test_fail "in first" "wrong count"
}

* Test 8.10: in last observations
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 91/100, generate(x_num)
count if !missing(x_num) in 91/100
if r(N) == 10 {
    noi test_pass "in last 10"
}
else {
    noi test_fail "in last" "wrong count"
}

/*******************************************************************************
 * SECTION 9: Multiple variables (10 tests)
 ******************************************************************************/
noi print_section "Multiple Variables"

* Test 9.1: Two variables with generate
clear
set obs 50
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
cdestring a b, generate(a_num b_num)
local pass = (a_num[1] == 1) & (b_num[1] == 2) & (a_num[50] == 50) & (b_num[50] == 100)
if `pass' {
    noi test_pass "two variables with generate"
}
else {
    noi test_fail "two vars gen" "wrong values"
}

* Test 9.2: Compare two variables with destring
clear
set obs 50
gen str10 a = string(_n)
gen str10 b = string(_n * 3)
noi benchmark_destring a b, testname("vs destring: two vars") generate(a b)

* Test 9.3: Three variables
clear
set obs 30
gen str10 x = string(_n)
gen str10 y = string(_n * 2)
gen str10 z = string(_n * 3)
noi benchmark_destring x y z, testname("vs destring: three vars") generate(x y z)

* Test 9.4: Multiple variables with replace
clear
set obs 40
gen str10 a = string(_n)
gen str10 b = string(_n + 100)
gen str10 c = string(_n + 200)
noi benchmark_destring a b c, testname("vs destring: multi replace") replace

* Test 9.5: Multiple variables with force
clear
set obs 30
gen str20 a = cond(mod(_n, 3) == 0, "NA", string(_n))
gen str20 b = cond(mod(_n, 4) == 0, "NA", string(_n * 2))
noi benchmark_destring a b, testname("vs destring: multi force") generate(a b) force

* Test 9.6: Multiple variables with ignore
clear
set obs 30
gen str20 a = "$" + string(_n)
gen str20 b = "$" + string(_n * 2)
noi benchmark_destring a b, testname("vs destring: multi ignore") generate(a b) ignore("$")

* Test 9.7: Multiple variables with percent
clear
set obs 30
gen str10 a = string(_n) + "%"
gen str10 b = string(_n * 2) + "%"
noi benchmark_destring a b, testname("vs destring: multi percent") generate(a b) percent

* Test 9.8: Multiple variables with dpcomma
clear
set obs 30
gen str20 a = subinstr(string(_n / 3, "%10.2f"), ".", ",", 1)
gen str20 b = subinstr(string(_n / 7, "%10.2f"), ".", ",", 1)
noi benchmark_destring a b, testname("vs destring: multi dpcomma") generate(a b) dpcomma

* Test 9.9: Five variables
clear
set obs 20
forvalues i = 1/5 {
    gen str10 v`i' = string(_n * `i')
}
noi benchmark_destring v1 v2 v3 v4 v5, testname("vs destring: five vars") generate(v1 v2 v3 v4 v5)

* Test 9.10: Multiple variables with if condition (destring doesn't support if/in)
clear
set obs 50
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
gen flag = mod(_n, 2)
cdestring a b if flag == 1, generate(a_num b_num)
count if !missing(a_num)
local n_a = r(N)
count if !missing(b_num)
local n_b = r(N)
* Should convert 25 observations (where flag==1)
if `n_a' == 25 & `n_b' == 25 {
    noi test_pass "multiple vars with if condition"
}
else {
    noi test_fail "multi if" "expected 25 each, got a=`n_a' b=`n_b'"
}

/*******************************************************************************
 * SECTION 10: Large datasets (10 tests)
 ******************************************************************************/
noi print_section "Large Datasets"

* Test 10.1: 10K observations
clear
set obs 10000
gen str10 x = string(_n)
capture cdestring x, generate(x_num)
if _rc == 0 {
    count if x_num == _n
    if r(N) == 10000 {
        noi test_pass "10K observations"
    }
    else {
        noi test_fail "10K obs" "wrong count"
    }
}
else {
    noi test_fail "10K obs" "rc=`=_rc'"
}

* Test 10.2: 50K observations
clear
set obs 50000
gen str10 x = string(_n)
capture cdestring x, generate(x_num)
if _rc == 0 {
    noi test_pass "50K observations"
}
else {
    noi test_fail "50K obs" "rc=`=_rc'"
}

* Test 10.3: 100K observations
clear
set obs 100000
gen str10 x = string(_n)
capture cdestring x, generate(x_num)
if _rc == 0 {
    noi test_pass "100K observations"
}
else {
    noi test_fail "100K obs" "rc=`=_rc'"
}

* Test 10.4: Large with decimals
clear
set obs 50000
gen str20 x = string(_n / 7, "%15.8f")
capture cdestring x, generate(x_num)
if _rc == 0 {
    noi test_pass "50K with decimals"
}
else {
    noi test_fail "50K decimals" "rc=`=_rc'"
}

* Test 10.5: Large with force
clear
set obs 50000
gen str20 x = cond(mod(_n, 100) == 0, "NA", string(_n))
capture cdestring x, generate(x_num) force
if _rc == 0 {
    count if missing(x_num)
    if r(N) == 500 {
        noi test_pass "50K with force"
    }
    else {
        noi test_fail "50K force" "wrong missing count"
    }
}
else {
    noi test_fail "50K force" "rc=`=_rc'"
}

* Test 10.6: Large with ignore
clear
set obs 50000
gen str20 x = "$" + string(_n)
capture cdestring x, generate(x_num) ignore("$")
if _rc == 0 {
    noi test_pass "50K with ignore"
}
else {
    noi test_fail "50K ignore" "rc=`=_rc'"
}

* Test 10.7: Large with percent
clear
set obs 50000
gen str15 x = string(_n) + "%"
capture cdestring x, generate(x_num) percent
if _rc == 0 {
    noi test_pass "50K with percent"
}
else {
    noi test_fail "50K percent" "rc=`=_rc'"
}

* Test 10.8: Large with threads
clear
set obs 100000
gen str10 x = string(_n)
capture cdestring x, generate(x_num) threads(4)
if _rc == 0 {
    noi test_pass "100K with threads(4)"
}
else {
    noi test_fail "100K threads" "rc=`=_rc'"
}

* Test 10.9: Large dataset comparison with destring
clear
set obs 10000
gen str10 x = string(_n)
noi benchmark_destring x, testname("vs destring: 10K obs") generate(x)

* Test 10.10: Large with multiple variables
clear
set obs 20000
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
gen str10 c = string(_n * 3)
capture cdestring a b c, generate(a_num b_num c_num)
if _rc == 0 {
    noi test_pass "20K with 3 variables"
}
else {
    noi test_fail "20K 3 vars" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 11: Edge cases and special values (10 tests)
 ******************************************************************************/
noi print_section "Edge Cases"

* Test 11.1: Scientific notation
clear
set obs 5
gen str20 x = ""
replace x = "1e5" in 1
replace x = "1E5" in 2
replace x = "1.5e3" in 3
replace x = "2.5E-2" in 4
replace x = "-1e10" in 5
noi benchmark_destring x, testname("vs destring: scientific notation") generate(x)

* Test 11.2: Leading zeros
clear
set obs 5
gen str10 x = ""
replace x = "001" in 1
replace x = "010" in 2
replace x = "0123" in 3
replace x = "00001" in 4
replace x = "000" in 5
noi benchmark_destring x, testname("vs destring: leading zeros") generate(x)

* Test 11.3: Trailing zeros after decimal
clear
set obs 5
gen str10 x = ""
replace x = "1.0" in 1
replace x = "2.00" in 2
replace x = "3.000" in 3
replace x = "4.0000" in 4
replace x = "0.10" in 5
noi benchmark_destring x, testname("vs destring: trailing zeros") generate(x)

* Test 11.4: Just decimal point
clear
set obs 5
gen str10 x = ""
replace x = ".5" in 1
replace x = ".25" in 2
replace x = ".125" in 3
replace x = ".0" in 4
replace x = "-.5" in 5
noi benchmark_destring x, testname("vs destring: leading decimal") generate(x)

* Test 11.5: Very large numbers
clear
set obs 5
gen str30 x = ""
replace x = "999999999999" in 1
replace x = "1234567890123" in 2
replace x = "9999999999999999" in 3
replace x = "-999999999999" in 4
replace x = "123456789012345678" in 5
noi benchmark_destring x, testname("vs destring: large numbers") generate(x)

* Test 11.6: Very small decimals
clear
set obs 5
gen str30 x = ""
replace x = "0.000001" in 1
replace x = "0.0000001" in 2
replace x = "0.00000001" in 3
replace x = "1.23e-10" in 4
replace x = "-0.000001" in 5
noi benchmark_destring x, testname("vs destring: small decimals") generate(x)

* Test 11.7: Plus sign prefix
clear
set obs 5
gen str10 x = ""
replace x = "+1" in 1
replace x = "+100" in 2
replace x = "+0.5" in 3
replace x = "+1e5" in 4
replace x = "5" in 5
noi benchmark_destring x, testname("vs destring: plus prefix") generate(x)

* Test 11.8: Zero variations
clear
set obs 6
gen str10 x = ""
replace x = "0" in 1
replace x = "00" in 2
replace x = "0.0" in 3
replace x = "-0" in 4
replace x = "+0" in 5
replace x = "0.00" in 6
cdestring x, generate(x_num)
count if x_num == 0
if r(N) == 6 {
    noi test_pass "zero variations"
}
else {
    noi test_fail "zeros" "not all zero"
}

* Test 11.9: Empty strings
clear
set obs 10
gen str10 x = cond(mod(_n, 2) == 0, "", string(_n))
cdestring x, generate(x_num) force
count if missing(x_num) & mod(_n, 2) == 0
if r(N) == 5 {
    noi test_pass "empty strings -> missing"
}
else {
    noi test_fail "empty strings" "wrong missing count"
}

* Test 11.10: Whitespace only strings
clear
set obs 5
gen str10 x = ""
replace x = " " in 1
replace x = "  " in 2
replace x = "   " in 3
replace x = " 1 " in 4
replace x = "2" in 5
noi benchmark_destring x, testname("vs destring: whitespace") generate(x) force

/*******************************************************************************
 * SECTION 12: Error handling (5 tests)
 ******************************************************************************/
noi print_section "Error Handling"

* Test 12.1: Generate var exists
clear
set obs 10
gen str10 x = string(_n)
gen x_num = 1
capture cdestring x, generate(x_num)
if _rc == 110 {
    noi test_pass "error: generate var exists"
}
else {
    noi test_fail "var exists" "expected rc=110, got `=_rc'"
}

* Test 12.2: Source is numeric (not string)
clear
set obs 10
gen x = _n
capture cdestring x, generate(x_num)
if _rc != 0 {
    noi test_pass "error: numeric source"
}
else {
    noi test_fail "numeric source" "should error"
}

* Test 12.3: Missing generate and replace
clear
set obs 10
gen str10 x = string(_n)
capture cdestring x
if _rc != 0 {
    noi test_pass "error: missing generate/replace"
}
else {
    noi test_fail "missing option" "should error"
}

* Test 12.4: Both generate and replace
clear
set obs 10
gen str10 x = string(_n)
capture cdestring x, generate(x_num) replace
if _rc != 0 {
    noi test_pass "error: both generate and replace"
}
else {
    noi test_fail "both options" "should error"
}

* Test 12.5: Generate count mismatch
clear
set obs 10
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
capture cdestring a b, generate(x_num)
if _rc != 0 {
    noi test_pass "error: generate count mismatch"
}
else {
    noi test_fail "count mismatch" "should error"
}

/*******************************************************************************
 * SECTION 13: Float option (5 tests)
 ******************************************************************************/
noi print_section "Float Option"

* Test 13.1: Float creates float variable
clear
set obs 10
gen str20 x = string(_n + 0.5)
cdestring x, generate(x_num) float
local vtype : type x_num
if "`vtype'" == "float" {
    noi test_pass "float creates float variable"
}
else {
    noi test_fail "float type" "type is `vtype'"
}

* Test 13.2: Compare float with destring
clear
set obs 50
gen str20 x = string(_n / 3)
noi benchmark_destring x, testname("vs destring: float option") generate(x) float

* Test 13.3: Float with replace
clear
set obs 30
gen str20 x = string(_n / 7)
noi benchmark_destring x, testname("vs destring: float replace") replace float

* Test 13.4: Float with integers (note: compress converts to byte)
clear
set obs 30
gen str10 x = string(_n)
* Both cdestring and destring compress integers to byte even with float option
cdestring x, generate(x_c) float
destring x, generate(x_s) float
local vtype_c : type x_c
local vtype_s : type x_s
* Test passes if cdestring matches destring behavior
if "`vtype_c'" == "`vtype_s'" {
    noi test_pass "float with integers (matches destring)"
}
else {
    noi test_fail "float int" "cdestring type `vtype_c' != destring type `vtype_s'"
}

* Test 13.5: Float precision
clear
set obs 5
gen str30 x = ""
replace x = "1.23456789" in 1
replace x = "9.87654321" in 2
replace x = "0.123456789" in 3
replace x = "123.456789" in 4
replace x = "0.000123456789" in 5
noi benchmark_destring x, testname("vs destring: float precision") generate(x) float

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cdestring"

if $TESTS_FAILED > 0 {
    exit 1
}

}
