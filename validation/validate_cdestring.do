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

noi di as text "Running validation tests for cdestring..."

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
        test_fail "`testname'" "must specify generate or replace"
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
            test_fail "`testname'" "cdestring rc=`rc_c', destring rc=`rc_s'"
            exit
        }

        if `rc_c' != 0 {
            test_pass "`testname' (both error as expected)"
            exit
        }

        * Compare results (use sigfigs-based comparison)
        local ndiff_total = 0
        foreach v of local varlist {
            tempvar _sf_`v'
            quietly {
                gen double `_sf_`v'' = 15 if `v'_c == `v'_s
                replace `_sf_`v'' = 0 if (`_sf_`v'' == .) & ((`v'_c == 0) | (`v'_s == 0))
                replace `_sf_`v'' = -log10(abs(`v'_c - `v'_s) / max(abs(`v'_c), abs(`v'_s))) if `_sf_`v'' == .
                replace `_sf_`v'' = 0 if `_sf_`v'' < 0
                replace `_sf_`v'' = 15 if `_sf_`v'' > 15
            }
            quietly count if `_sf_`v'' < $DEFAULT_SIGFIGS & !missing(`v'_c) & !missing(`v'_s)
            local ndiff_total = `ndiff_total' + r(N)
            * Also check missing patterns match
            quietly count if missing(`v'_c) != missing(`v'_s)
            local ndiff_total = `ndiff_total' + r(N)
        }

        if `ndiff_total' > 0 {
            test_fail "`testname'" "`ndiff_total' values differ"
        }
        else {
            test_pass "`testname'"
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
            test_fail "`testname'" "cdestring rc=`rc_c', destring rc=`rc_s'"
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
            test_pass "`testname' (both error as expected)"
            exit
        }

        * Compare numeric values
        local ndiff_total = 0
        local i = 1
        foreach v of local generate {
            local vc = "`v'_c"
            local vs = "`v'_s"
            * Compare non-missing values with sigfigs
            tempvar _sf_`v'
            quietly {
                gen double `_sf_`v'' = 15 if `vc' == `vs'
                replace `_sf_`v'' = 0 if (`_sf_`v'' == .) & ((`vc' == 0) | (`vs' == 0))
                replace `_sf_`v'' = -log10(abs(`vc' - `vs') / max(abs(`vc'), abs(`vs'))) if `_sf_`v'' == .
                replace `_sf_`v'' = 0 if `_sf_`v'' < 0
                replace `_sf_`v'' = 15 if `_sf_`v'' > 15
            }
            quietly count if `_sf_`v'' < $DEFAULT_SIGFIGS & !missing(`vc') & !missing(`vs')
            local ndiff_total = `ndiff_total' + r(N)
            * Check missing patterns match
            quietly count if missing(`vc') != missing(`vs')
            local ndiff_total = `ndiff_total' + r(N)
            local ++i
        }

        if `ndiff_total' > 0 {
            test_fail "`testname'" "`ndiff_total' values differ"
        }
        else {
            test_pass "`testname'"
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
print_section "Plugin Check"

clear
set obs 10
gen str10 numstr = string(_n)
capture noisily cdestring numstr, generate(numval)
if _rc != 0 {
    test_fail "cdestring plugin load" "plugin returned error `=_rc'"
    print_summary "cdestring"
    exit 1
}
test_pass "cdestring plugin loads and runs"
drop numval

/*******************************************************************************
 * SECTION 2: Basic functionality (10 tests)
 ******************************************************************************/
print_section "Basic Functionality"

* Test 2.1: Basic conversion creates numeric variable
clear
set obs 10
gen str10 x = string(_n)
cdestring x, generate(x_num)
capture confirm numeric variable x_num
if _rc == 0 {
    test_pass "basic conversion creates numeric variable"
}
else {
    test_fail "basic conversion" "variable not numeric"
}

* Test 2.2: Compare with destring - simple integers
clear
set obs 100
gen str10 x = string(_n)
benchmark_destring x, testname("vs destring: simple integers") generate(x)

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
    test_pass "values correctly converted"
}
else {
    test_fail "values converted" "mismatch found"
}

* Test 2.4: Single observation
clear
set obs 1
gen str10 x = "42"
cdestring x, generate(x_num)
if x_num[1] == 42 {
    test_pass "single observation"
}
else {
    test_fail "single observation" "got `=x_num[1]'"
}

* Test 2.5: verbose option
clear
set obs 10
gen str10 x = string(_n)
capture cdestring x, generate(x_num_c) verbose
local rc_c = _rc
capture destring x, generate(x_num_s)
local rc_s = _rc
if `rc_c' != 0 | `rc_s' != 0 {
    test_fail "verbose option" "cdestring rc=`rc_c', destring rc=`rc_s'"
}
else {
    assert_var_equal x_num_c x_num_s $DEFAULT_SIGFIGS "verbose option"
}

* Test 2.6: threads option
clear
set obs 100
gen str10 x = string(_n)
capture cdestring x, generate(x_num_c) threads(2)
local rc_c = _rc
capture destring x, generate(x_num_s)
local rc_s = _rc
if `rc_c' != 0 | `rc_s' != 0 {
    test_fail "threads(2) option" "cdestring rc=`rc_c', destring rc=`rc_s'"
}
else {
    assert_var_equal x_num_c x_num_s $DEFAULT_SIGFIGS "threads(2) option"
}

* Test 2.7: Negative numbers
clear
set obs 10
gen str10 x = string(_n - 5)
cdestring x, generate(x_num)
if x_num[1] == -4 & x_num[5] == 0 & x_num[10] == 5 {
    test_pass "negative numbers"
}
else {
    test_fail "negative numbers" "wrong values"
}

* Test 2.8: Compare negative numbers with destring
clear
set obs 50
gen str10 x = string(_n - 25)
benchmark_destring x, testname("vs destring: negative integers") generate(x)

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
local pass = (abs(x_num[1] - 1.5) < $DEFAULT_TOL) & (abs(x_num[3] - 3.125) < $DEFAULT_TOL)
if `pass' {
    test_pass "decimal numbers"
}
else {
    test_fail "decimal numbers" "wrong values"
}

* Test 2.10: Compare decimals with destring
clear
set obs 100
gen str20 x = string(_n / 7, "%20.10f")
benchmark_destring x, testname("vs destring: decimal numbers") generate(x)

/*******************************************************************************
 * SECTION 3: Replace option (10 tests)
 ******************************************************************************/
print_section "Replace Option"

* Test 3.1: Replace converts in place
clear
set obs 10
gen str10 x = string(_n)
cdestring x, replace
capture confirm numeric variable x
if _rc == 0 {
    test_pass "replace converts in place"
}
else {
    test_fail "replace option" "variable still string"
}

* Test 3.2: Replace preserves values
clear
set obs 5
gen str10 x = string(_n * 100)
cdestring x, replace
local pass = (x[1] == 100) & (x[3] == 300) & (x[5] == 500)
if `pass' {
    test_pass "replace preserves values"
}
else {
    test_fail "replace values" "wrong conversion"
}

* Test 3.3: Compare replace with destring
clear
set obs 50
gen str10 x = string(_n)
benchmark_destring x, testname("vs destring: replace option") replace

* Test 3.4: Replace with negative numbers
clear
set obs 20
gen str10 x = string(_n - 10)
benchmark_destring x, testname("vs destring: replace negative") replace

* Test 3.5: Replace with decimals
clear
set obs 30
gen str20 x = string(_n / 3, "%10.6f")
benchmark_destring x, testname("vs destring: replace decimals") replace

* Test 3.6: Replace multiple variables
clear
set obs 20
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
gen str10 c = string(_n * 3)
benchmark_destring a b c, testname("vs destring: replace multi") replace

* Test 3.7: Replace with if condition
clear
set obs 50
gen str10 x = string(_n)
gen flag = mod(_n, 2)
cdestring x if flag == 1, replace
* Should only convert odd observations (flag==1)
* Even observations (flag==0) become missing after replace
count if missing(x) & flag == 0
local n_missing = r(N)
count if !missing(x) & flag == 1
local n_converted = r(N)
if `n_missing' == 25 & `n_converted' == 25 {
    test_pass "replace with if condition"
}
else {
    test_fail "replace if" "missing=`n_missing' converted=`n_converted'"
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
        test_pass "replace with in range"
    }
    else {
        test_fail "replace in" "wrong count"
    }
}
else {
    test_fail "replace in" "not converted"
}

* Test 3.9: Replace float option
clear
set obs 10
gen str20 x = string(_n + 0.5)
cdestring x, replace float
local vtype : type x
if "`vtype'" == "float" {
    test_pass "replace with float option"
}
else {
    test_fail "replace float" "type is `vtype'"
}

* Test 3.10: Compare replace float with destring
clear
set obs 50
gen str20 x = string(_n / 2)
benchmark_destring x, testname("vs destring: replace float") replace float

/*******************************************************************************
 * SECTION 4: Ignore option (10 tests)
 ******************************************************************************/
print_section "Ignore Option"

* Test 4.1: Ignore single character
clear
set obs 5
gen str20 x = "$" + string(_n * 100)
cdestring x, generate(x_num) ignore("$")
if x_num[1] == 100 & x_num[5] == 500 {
    test_pass "ignore single character ($)"
}
else {
    test_fail "ignore single" "wrong values"
}

* Test 4.2: Compare ignore with destring
clear
set obs 50
gen str20 x = "$" + string(_n)
benchmark_destring x, testname("vs destring: ignore $") generate(x) ignore("$")

* Test 4.3: Ignore multiple characters
clear
set obs 5
gen str30 x = ""
replace x = "$1,000" in 1
replace x = "$2,500" in 2
replace x = "$10,000" in 3
replace x = "$25,000" in 4
replace x = "$100,000" in 5
cdestring x, generate(x_num) ignore("$,")
* Should strip $ and , to get numeric values
local pass = (x_num[1] == 1000) & (x_num[2] == 2500) & (x_num[3] == 10000)
if `pass' {
    test_pass "ignore multiple characters"
}
else {
    test_fail "ignore multiple" "wrong values"
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
benchmark_destring x, testname("vs destring: ignore comma") generate(x) ignore(",")

* Test 4.5: Ignore space
clear
set obs 5
gen str20 x = ""
replace x = "1 000" in 1
replace x = "2 500" in 2
replace x = "10 000" in 3
replace x = "100 000" in 4
replace x = "1 234 567" in 5
benchmark_destring x, testname("vs destring: ignore space") generate(x) ignore(" ")

* Test 4.6: Ignore letters (force may be needed)
clear
set obs 5
gen str20 x = ""
replace x = "100kg" in 1
replace x = "200kg" in 2
replace x = "300kg" in 3
replace x = "400kg" in 4
replace x = "500kg" in 5
benchmark_destring x, testname("vs destring: ignore kg") generate(x) ignore("kg")

* Test 4.7: Ignore percent sign
clear
set obs 5
gen str10 x = string(_n * 10) + "%"
benchmark_destring x, testname("vs destring: ignore %") generate(x) ignore("%")

* Test 4.8: Ignore with replace
clear
set obs 50
gen str20 x = "$" + string(_n * 100)
benchmark_destring x, testname("vs destring: ignore $ replace") replace ignore("$")

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
benchmark_destring x, testname("vs destring: ignore parens") generate(x) ignore("()")

* Test 4.10: Empty ignore string (no effect)
clear
set obs 20
gen str10 x = string(_n)
benchmark_destring x, testname("vs destring: empty ignore") generate(x) ignore("")

/*******************************************************************************
 * SECTION 5: Force option (10 tests)
 ******************************************************************************/
print_section "Force Option"

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
    test_pass "force converts non-numeric to missing"
}
else {
    test_fail "force option" "wrong values"
}

* Test 5.2: Compare force with destring
clear
set obs 50
gen str20 x = cond(mod(_n, 3) == 0, "abc", string(_n))
benchmark_destring x, testname("vs destring: force option") generate(x) force

* Test 5.3: Force with all non-numeric
clear
set obs 10
gen str10 x = char(64 + _n)  // A, B, C, ...
cdestring x, generate(x_num) force
count if missing(x_num)
if r(N) == 10 {
    test_pass "force with all non-numeric"
}
else {
    test_fail "force all non-numeric" "expected 10 missing"
}

* Test 5.4: Force with mixed content
clear
set obs 20
gen str20 x = ""
replace x = string(_n) if mod(_n, 2) == 1
replace x = "text" + string(_n) if mod(_n, 2) == 0
benchmark_destring x, testname("vs destring: force mixed") generate(x) force

* Test 5.5: Force with leading/trailing spaces
clear
set obs 5
gen str20 x = ""
replace x = " 100" in 1
replace x = "200 " in 2
replace x = " 300 " in 3
replace x = "  400  " in 4
replace x = "500" in 5
benchmark_destring x, testname("vs destring: force spaces") generate(x) force

* Test 5.6: Force with replace
clear
set obs 30
gen str20 x = cond(mod(_n, 4) == 0, "NA", string(_n))
benchmark_destring x, testname("vs destring: force replace") replace force

* Test 5.7: Force with ignore
clear
set obs 20
gen str20 x = "$" + cond(mod(_n, 5) == 0, "NA", string(_n))
benchmark_destring x, testname("vs destring: force + ignore") generate(x) force ignore("$")

* Test 5.8: Force preserves valid numbers
clear
set obs 100
gen str10 x = string(_n)
cdestring x, generate(x_num) force
count if x_num == real(x)
if r(N) == 100 {
    test_pass "force preserves valid numbers"
}
else {
    test_fail "force valid" "some values changed"
}

* Test 5.9: Force with empty strings
clear
set obs 10
gen str10 x = cond(mod(_n, 2) == 0, "", string(_n))
benchmark_destring x, testname("vs destring: force empty") generate(x) force

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
benchmark_destring x, testname("vs destring: force special chars") generate(x) force

/*******************************************************************************
 * SECTION 6: Percent option (10 tests)
 ******************************************************************************/
print_section "Percent Option"

* Test 6.1: Percent strips % and divides by 100
clear
set obs 5
gen str10 x = string(_n * 10) + "%"
cdestring x, generate(x_num) percent
local pass = (abs(x_num[1] - 0.1) < $DEFAULT_TOL) & (abs(x_num[5] - 0.5) < $DEFAULT_TOL)
if `pass' {
    test_pass "percent option basic"
}
else {
    test_fail "percent option" "wrong values"
}

* Test 6.2: Compare percent with destring
clear
set obs 50
gen str10 x = string(_n) + "%"
benchmark_destring x, testname("vs destring: percent option") generate(x) percent

* Test 6.3: Percent with decimals
clear
set obs 5
gen str10 x = ""
replace x = "10.5%" in 1
replace x = "25.25%" in 2
replace x = "50.125%" in 3
replace x = "75.0%" in 4
replace x = "99.99%" in 5
benchmark_destring x, testname("vs destring: percent decimal") generate(x) percent

* Test 6.4: Percent with replace
clear
set obs 30
gen str10 x = string(_n * 3) + "%"
benchmark_destring x, testname("vs destring: percent replace") replace percent

* Test 6.5: Percent without % sign (no transformation)
clear
set obs 10
gen str10 x = string(_n * 10)
cdestring x, generate(x_num) percent
* Values without % should not be divided by 100
if x_num[1] == 10 & x_num[10] == 100 {
    test_pass "percent without % sign"
}
else {
    test_fail "percent no %" "values changed unexpectedly"
}

* Test 6.6: Compare percent without % with destring
clear
set obs 50
gen str10 x = string(_n)
benchmark_destring x, testname("vs destring: percent no %") generate(x) percent

* Test 6.7: Percent with negative values
clear
set obs 5
gen str10 x = ""
replace x = "-10%" in 1
replace x = "-25%" in 2
replace x = "0%" in 3
replace x = "50%" in 4
replace x = "-100%" in 5
benchmark_destring x, testname("vs destring: percent negative") generate(x) percent

* Test 6.8: Percent with force
clear
set obs 10
gen str10 x = cond(mod(_n, 3) == 0, "NA%", string(_n * 10) + "%")
benchmark_destring x, testname("vs destring: percent force") generate(x) percent force

* Test 6.9: 100% value
clear
set obs 5
gen str10 x = "100%"
cdestring x, generate(x_num) percent
if abs(x_num[1] - 1.0) < $DEFAULT_TOL {
    test_pass "100% converts to 1.0"
}
else {
    test_fail "100% value" "got `=x_num[1]'"
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
benchmark_destring x, testname("vs destring: percent small") generate(x) percent

/*******************************************************************************
 * SECTION 7: dpcomma option (10 tests)
 ******************************************************************************/
print_section "dpcomma Option"

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
local pass = (abs(x_num[1] - 1.5) < $DEFAULT_TOL) & (abs(x_num[3] - 3.125) < $DEFAULT_TOL)
if `pass' {
    test_pass "dpcomma basic"
}
else {
    test_fail "dpcomma option" "wrong values"
}

* Test 7.2: Compare dpcomma with destring
clear
set obs 50
gen str20 x = subinstr(string(_n / 3, "%10.3f"), ".", ",", 1)
benchmark_destring x, testname("vs destring: dpcomma option") generate(x) dpcomma

* Test 7.3: dpcomma with negative numbers
clear
set obs 5
gen str10 x = ""
replace x = "-1,5" in 1
replace x = "-2,25" in 2
replace x = "0,0" in 3
replace x = "1,5" in 4
replace x = "-0,5" in 5
benchmark_destring x, testname("vs destring: dpcomma negative") generate(x) dpcomma

* Test 7.4: dpcomma with replace
clear
set obs 30
gen str20 x = subinstr(string(_n / 7, "%10.4f"), ".", ",", 1)
benchmark_destring x, testname("vs destring: dpcomma replace") replace dpcomma

* Test 7.5: dpcomma with integers (no comma)
clear
set obs 20
gen str10 x = string(_n)
cdestring x, generate(x_num) dpcomma
if x_num[1] == 1 & x_num[20] == 20 {
    test_pass "dpcomma with integers"
}
else {
    test_fail "dpcomma int" "wrong values"
}

* Test 7.6: dpcomma with force
clear
set obs 10
gen str20 x = cond(mod(_n, 3) == 0, "abc", subinstr(string(_n / 2, "%10.2f"), ".", ",", 1))
benchmark_destring x, testname("vs destring: dpcomma force") generate(x) dpcomma force

* Test 7.7: dpcomma with ignore (period as thousand separator)
clear
set obs 5
gen str20 x = ""
replace x = "1.000,50" in 1
replace x = "2.500,25" in 2
replace x = "10.000,00" in 3
replace x = "100.000,99" in 4
replace x = "1.234.567,89" in 5
benchmark_destring x, testname("vs destring: dpcomma + ignore .") generate(x) dpcomma ignore(".")

* Test 7.8: dpcomma zero
clear
set obs 3
gen str10 x = ""
replace x = "0,0" in 1
replace x = "0,00" in 2
replace x = "0,000" in 3
cdestring x, generate(x_num) dpcomma
if x_num[1] == 0 & x_num[2] == 0 & x_num[3] == 0 {
    test_pass "dpcomma zero values"
}
else {
    test_fail "dpcomma zero" "not zero"
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
benchmark_destring x, testname("vs destring: dpcomma percent") generate(x) dpcomma percent

* Test 7.10: Compare standard decimal with dpcomma
clear
set obs 30
* Standard decimal
gen str20 x_dot = string(_n / 7, "%10.4f")
* Comma decimal
gen str20 x_comma = subinstr(x_dot, ".", ",", 1)
destring x_dot, generate(x_dot_num)
cdestring x_comma, generate(x_comma_num) dpcomma
count if abs(x_dot_num - x_comma_num) > $DEFAULT_TOL
if r(N) == 0 {
    test_pass "dpcomma equals standard decimal"
}
else {
    test_fail "dpcomma vs dot" "values differ"
}

/*******************************************************************************
 * SECTION 8: if/in conditions (10 tests)
 ******************************************************************************/
print_section "if/in Conditions"

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
    test_pass "if condition subset"
}
else {
    test_fail "if condition" "wrong pattern"
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
    test_pass "in range"
}
else {
    test_fail "in range" "wrong counts"
}

* Test 8.3: if condition works correctly (destring doesn't support if/in)
clear
set seed 12345
set obs 100
gen str10 x = string(_n)
gen value = runiform()
cdestring x if value > 0.5, generate(x_num)
* Count how many observations match the condition
count if value > 0.5
local n_expected = r(N)
count if !missing(x_num)
local n_converted = r(N)
if `n_converted' == `n_expected' {
    test_pass "if condition filters observations"
}
else {
    test_fail "if condition" "expected `n_expected' converted, got `n_converted'"
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
    test_pass "in range works correctly"
}
else {
    test_fail "in range" "in=`in_range' miss_before=`miss_before' miss_after=`miss_after'"
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
    test_pass "combined if and in"
}
else {
    test_fail "if+in" "expected 16, got `n_converted'"
}

* Test 8.6: Single observation in
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 50/50, generate(x_num)
count if !missing(x_num)
if r(N) == 1 & x_num[50] == 50 {
    test_pass "single observation in"
}
else {
    test_fail "single in" "wrong result"
}

* Test 8.7: if that matches no observations
clear
set obs 50
gen str10 x = string(_n)
gen value = _n
cdestring x if value > 1000, generate(x_num)
count if !missing(x_num)
if r(N) == 0 {
    test_pass "if matches no observations"
}
else {
    test_fail "if none" "some values converted"
}

* Test 8.8: if with string comparison on another variable
clear
set obs 50
gen str10 x = string(_n)
gen str10 group = cond(mod(_n, 2) == 0, "even", "odd")
cdestring x if group == "even", generate(x_num)
count if !missing(x_num)
if r(N) == 25 {
    test_pass "if with string comparison"
}
else {
    test_fail "if string" "wrong count"
}

* Test 8.9: in first observations
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 1/10, generate(x_num)
count if !missing(x_num)
if r(N) == 10 {
    test_pass "in first 10"
}
else {
    test_fail "in first" "wrong count"
}

* Test 8.10: in last observations
clear
set obs 100
gen str10 x = string(_n)
cdestring x in 91/100, generate(x_num)
count if !missing(x_num) in 91/100
if r(N) == 10 {
    test_pass "in last 10"
}
else {
    test_fail "in last" "wrong count"
}

/*******************************************************************************
 * SECTION 9: Multiple variables (10 tests)
 ******************************************************************************/
print_section "Multiple Variables"

* Test 9.1: Two variables with generate
clear
set obs 50
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
cdestring a b, generate(a_num b_num)
local pass = (a_num[1] == 1) & (b_num[1] == 2) & (a_num[50] == 50) & (b_num[50] == 100)
if `pass' {
    test_pass "two variables with generate"
}
else {
    test_fail "two vars gen" "wrong values"
}

* Test 9.2: Compare two variables with destring
clear
set obs 50
gen str10 a = string(_n)
gen str10 b = string(_n * 3)
benchmark_destring a b, testname("vs destring: two vars") generate(a b)

* Test 9.3: Three variables
clear
set obs 30
gen str10 x = string(_n)
gen str10 y = string(_n * 2)
gen str10 z = string(_n * 3)
benchmark_destring x y z, testname("vs destring: three vars") generate(x y z)

* Test 9.4: Multiple variables with replace
clear
set obs 40
gen str10 a = string(_n)
gen str10 b = string(_n + 100)
gen str10 c = string(_n + 200)
benchmark_destring a b c, testname("vs destring: multi replace") replace

* Test 9.5: Multiple variables with force
clear
set obs 30
gen str20 a = cond(mod(_n, 3) == 0, "NA", string(_n))
gen str20 b = cond(mod(_n, 4) == 0, "NA", string(_n * 2))
benchmark_destring a b, testname("vs destring: multi force") generate(a b) force

* Test 9.6: Multiple variables with ignore
clear
set obs 30
gen str20 a = "$" + string(_n)
gen str20 b = "$" + string(_n * 2)
benchmark_destring a b, testname("vs destring: multi ignore") generate(a b) ignore("$")

* Test 9.7: Multiple variables with percent
clear
set obs 30
gen str10 a = string(_n) + "%"
gen str10 b = string(_n * 2) + "%"
benchmark_destring a b, testname("vs destring: multi percent") generate(a b) percent

* Test 9.8: Multiple variables with dpcomma
clear
set obs 30
gen str20 a = subinstr(string(_n / 3, "%10.2f"), ".", ",", 1)
gen str20 b = subinstr(string(_n / 7, "%10.2f"), ".", ",", 1)
benchmark_destring a b, testname("vs destring: multi dpcomma") generate(a b) dpcomma

* Test 9.9: Five variables
clear
set obs 20
forvalues i = 1/5 {
    gen str10 v`i' = string(_n * `i')
}
benchmark_destring v1 v2 v3 v4 v5, testname("vs destring: five vars") generate(v1 v2 v3 v4 v5)

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
    test_pass "multiple vars with if condition"
}
else {
    test_fail "multi if" "expected 25 each, got a=`n_a' b=`n_b'"
}

/*******************************************************************************
 * SECTION 10: Large datasets (10 tests)
 ******************************************************************************/
print_section "Large Datasets"

* Test 10.1: 10K observations
clear
set obs 10000
gen str10 x = string(_n)
capture cdestring x, generate(x_num)
if _rc == 0 {
    count if x_num == _n
    if r(N) == 10000 {
        test_pass "10K observations"
    }
    else {
        test_fail "10K obs" "wrong count"
    }
}
else {
    test_fail "10K obs" "rc=`=_rc'"
}

* Test 10.2: 50K observations
clear
set obs 50000
gen str10 x = string(_n)
capture cdestring x, generate(x_num)
if _rc == 0 {
    count if x_num == _n
    if r(N) == 50000 {
        test_pass "50K observations"
    }
    else {
        test_fail "50K obs" "only `r(N)' correct"
    }
}
else {
    test_fail "50K obs" "rc=`=_rc'"
}

* Test 10.3: 100K observations
clear
set obs 100000
gen str10 x = string(_n)
capture cdestring x, generate(x_num)
if _rc == 0 {
    count if x_num == _n
    if r(N) == 100000 {
        test_pass "100K observations"
    }
    else {
        test_fail "100K obs" "only `r(N)' correct"
    }
}
else {
    test_fail "100K obs" "rc=`=_rc'"
}

* Test 10.4: Large with decimals
clear
set obs 50000
gen str20 x = string(_n / 7, "%15.8f")
capture cdestring x, generate(x_num)
if _rc == 0 {
    test_pass "50K with decimals"
}
else {
    test_fail "50K decimals" "rc=`=_rc'"
}

* Test 10.5: Large with force
clear
set obs 50000
gen str20 x = cond(mod(_n, 100) == 0, "NA", string(_n))
capture cdestring x, generate(x_num) force
if _rc == 0 {
    count if missing(x_num)
    if r(N) == 500 {
        test_pass "50K with force"
    }
    else {
        test_fail "50K force" "wrong missing count"
    }
}
else {
    test_fail "50K force" "rc=`=_rc'"
}

* Test 10.6: Large with ignore
clear
set obs 50000
gen str20 x = "$" + string(_n)
capture cdestring x, generate(x_num) ignore("$")
if _rc == 0 {
    count if x_num == _n
    if r(N) == 50000 {
        test_pass "50K with ignore"
    }
    else {
        test_fail "50K ignore" "only `r(N)' correct"
    }
}
else {
    test_fail "50K ignore" "rc=`=_rc'"
}

* Test 10.7: Large with percent
clear
set obs 50000
gen str15 x = string(_n) + "%"
capture cdestring x, generate(x_num) percent
if _rc == 0 {
    * Values should be _n / 100
    count if abs(x_num - _n/100) < $DEFAULT_TOL
    if r(N) == 50000 {
        test_pass "50K with percent"
    }
    else {
        test_fail "50K percent" "only `r(N)' correct"
    }
}
else {
    test_fail "50K percent" "rc=`=_rc'"
}

* Test 10.8: Large with threads
clear
set obs 100000
gen str10 x = string(_n)
capture cdestring x, generate(x_num) threads(4)
if _rc == 0 {
    test_pass "100K with threads(4)"
}
else {
    test_fail "100K threads" "rc=`=_rc'"
}

* Test 10.9: Large dataset comparison with destring
clear
set obs 10000
gen str10 x = string(_n)
benchmark_destring x, testname("vs destring: 10K obs") generate(x)

* Test 10.10: Large with multiple variables
clear
set obs 20000
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
gen str10 c = string(_n * 3)
capture cdestring a b c, generate(a_num b_num c_num)
if _rc == 0 {
    test_pass "20K with 3 variables"
}
else {
    test_fail "20K 3 vars" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 11: Edge cases and special values (10 tests)
 ******************************************************************************/
print_section "Edge Cases"

* Test 11.1: Scientific notation
clear
set obs 5
gen str20 x = ""
replace x = "1e5" in 1
replace x = "1E5" in 2
replace x = "1.5e3" in 3
replace x = "2.5E-2" in 4
replace x = "-1e10" in 5
benchmark_destring x, testname("vs destring: scientific notation") generate(x)

* Test 11.2: Leading zeros
clear
set obs 5
gen str10 x = ""
replace x = "001" in 1
replace x = "010" in 2
replace x = "0123" in 3
replace x = "00001" in 4
replace x = "000" in 5
benchmark_destring x, testname("vs destring: leading zeros") generate(x)

* Test 11.3: Trailing zeros after decimal
clear
set obs 5
gen str10 x = ""
replace x = "1.0" in 1
replace x = "2.00" in 2
replace x = "3.000" in 3
replace x = "4.0000" in 4
replace x = "0.10" in 5
benchmark_destring x, testname("vs destring: trailing zeros") generate(x)

* Test 11.4: Just decimal point
clear
set obs 5
gen str10 x = ""
replace x = ".5" in 1
replace x = ".25" in 2
replace x = ".125" in 3
replace x = ".0" in 4
replace x = "-.5" in 5
benchmark_destring x, testname("vs destring: leading decimal") generate(x)

* Test 11.5: Very large numbers
clear
set obs 5
gen str30 x = ""
replace x = "999999999999" in 1
replace x = "1234567890123" in 2
replace x = "9999999999999999" in 3
replace x = "-999999999999" in 4
replace x = "123456789012345678" in 5
benchmark_destring x, testname("vs destring: large numbers") generate(x)

* Test 11.6: Very small decimals
clear
set obs 5
gen str30 x = ""
replace x = "0.000001" in 1
replace x = "0.0000001" in 2
replace x = "0.00000001" in 3
replace x = "1.23e-10" in 4
replace x = "-0.000001" in 5
benchmark_destring x, testname("vs destring: small decimals") generate(x)

* Test 11.7: Plus sign prefix
clear
set obs 5
gen str10 x = ""
replace x = "+1" in 1
replace x = "+100" in 2
replace x = "+0.5" in 3
replace x = "+1e5" in 4
replace x = "5" in 5
benchmark_destring x, testname("vs destring: plus prefix") generate(x)

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
    test_pass "zero variations"
}
else {
    test_fail "zeros" "not all zero"
}

* Test 11.9: Empty strings
clear
set obs 10
gen str10 x = cond(mod(_n, 2) == 0, "", string(_n))
cdestring x, generate(x_num) force
count if missing(x_num) & mod(_n, 2) == 0
if r(N) == 5 {
    test_pass "empty strings -> missing"
}
else {
    test_fail "empty strings" "wrong missing count"
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
benchmark_destring x, testname("vs destring: whitespace") generate(x) force

/*******************************************************************************
 * SECTION 12: Error handling (5 tests)
 ******************************************************************************/
print_section "Error Handling"

* Test 12.1: Generate var exists
clear
set obs 10
gen str10 x = string(_n)
gen x_num = 1
capture cdestring x, generate(x_num)
if _rc == 110 {
    test_pass "error: generate var exists"
}
else {
    test_fail "var exists" "expected rc=110, got `=_rc'"
}

* Test 12.2: Source is numeric (match destring behavior - succeed gracefully)
clear
set obs 10
gen x = _n
capture destring x, generate(x_stata)
local stata_rc = _rc
capture cdestring x, generate(x_ctools)
local ctools_rc = _rc
if `stata_rc' == `ctools_rc' {
    if `ctools_rc' == 0 {
        * destring on numeric input may succeed without creating the variable;
        * check whether both variables exist before comparing values
        capture confirm variable x_stata
        local has_stata = (_rc == 0)
        capture confirm variable x_ctools
        local has_ctools = (_rc == 0)
        if `has_stata' & `has_ctools' {
            quietly count if x_stata != x_ctools & !missing(x_stata) & !missing(x_ctools)
            if r(N) == 0 {
                test_pass "numeric source handled (matches destring)"
            }
            else {
                test_fail "numeric source" "`=r(N)' values differ"
            }
        }
        else if `has_stata' == `has_ctools' {
            test_pass "numeric source handled (matches destring, no var created)"
        }
        else {
            test_fail "numeric source" "var creation mismatch: destring=`has_stata' cdestring=`has_ctools'"
        }
        capture drop x_stata x_ctools
    }
    else {
        test_pass "numeric source (both error rc=`ctools_rc')"
    }
}
else {
    test_fail "numeric source" "rc differ: destring=`stata_rc' cdestring=`ctools_rc'"
}

* Test 12.3: Missing generate and replace
clear
set obs 10
gen str10 x = string(_n)
capture cdestring x
if _rc != 0 {
    test_pass "error: missing generate/replace"
}
else {
    test_fail "missing option" "should error"
}

* Test 12.4: Both generate and replace
clear
set obs 10
gen str10 x = string(_n)
capture cdestring x, generate(x_num) replace
if _rc != 0 {
    test_pass "error: both generate and replace"
}
else {
    test_fail "both options" "should error"
}

* Test 12.5: Generate count mismatch
clear
set obs 10
gen str10 a = string(_n)
gen str10 b = string(_n * 2)
capture cdestring a b, generate(x_num)
if _rc != 0 {
    test_pass "error: generate count mismatch"
}
else {
    test_fail "count mismatch" "should error"
}

/*******************************************************************************
 * SECTION 13: Float option (5 tests)
 ******************************************************************************/
print_section "Float Option"

* Test 13.1: Float creates float variable
clear
set obs 10
gen str20 x = string(_n + 0.5)
cdestring x, generate(x_num) float
local vtype : type x_num
if "`vtype'" == "float" {
    test_pass "float creates float variable"
}
else {
    test_fail "float type" "type is `vtype'"
}

* Test 13.2: Compare float with destring
clear
set obs 50
gen str20 x = string(_n / 3)
benchmark_destring x, testname("vs destring: float option") generate(x) float

* Test 13.3: Float with replace
clear
set obs 30
gen str20 x = string(_n / 7)
benchmark_destring x, testname("vs destring: float replace") replace float

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
    test_pass "float with integers (matches destring)"
}
else {
    test_fail "float int" "cdestring type `vtype_c' != destring type `vtype_s'"
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
benchmark_destring x, testname("vs destring: float precision") generate(x) float

/*******************************************************************************
 * SECTION 14: Empty dataset and boundary tests
 ******************************************************************************/
print_section "Empty Dataset and Boundary Tests"

* Test 14.1: Empty dataset with generate - compare with Stata
clear
set obs 0
gen str10 x = ""
capture destring x, generate(x_num_s)
local stata_rc = _rc
capture cdestring x, generate(x_num_c)
local cdestring_rc = _rc
if `stata_rc' == `cdestring_rc' {
    test_pass "empty dataset with generate - matches Stata"
}
else {
    test_fail "empty gen" "cdestring rc=`cdestring_rc' but destring rc=`stata_rc'"
}

* Test 14.2: Empty dataset with replace - compare with Stata
clear
set obs 0
gen str10 x = ""
capture clonevar x_s = x
capture destring x_s, replace
local stata_rc = _rc
capture clonevar x_c = x
capture cdestring x_c, replace
local cdestring_rc = _rc
if `stata_rc' == `cdestring_rc' {
    test_pass "empty dataset with replace - matches Stata"
}
else {
    test_fail "empty replace" "cdestring rc=`cdestring_rc' but destring rc=`stata_rc'"
}

* Test 14.3: Single observation
clear
set obs 1
gen str10 x = "42"
benchmark_destring x, testname("single observation") generate(x)

* Test 14.4: Single observation with non-numeric
clear
set obs 1
gen str10 x = "abc"
benchmark_destring x, testname("single non-numeric with force") generate(x) force

* Test 14.5: Two observations
clear
set obs 2
gen str10 x = string(_n)
benchmark_destring x, testname("two observations") generate(x)

/*******************************************************************************
 * SECTION 15: Comprehensive option combinations
 ******************************************************************************/
print_section "Option Combinations"

* Test 15.1: ignore + force
clear
set obs 20
gen str20 x = "$" + cond(mod(_n, 4) == 0, "NA", string(_n * 10))
benchmark_destring x, testname("ignore + force") generate(x) ignore("$") force

* Test 15.2: ignore + percent
clear
set obs 20
gen str20 x = "$" + string(_n) + "%"
benchmark_destring x, testname("ignore $ + percent") generate(x) ignore("$") percent

* Test 15.3: dpcomma + force
clear
set obs 20
gen str20 x = cond(mod(_n, 5) == 0, "NA", subinstr(string(_n / 3, "%10.2f"), ".", ",", 1))
benchmark_destring x, testname("dpcomma + force") generate(x) dpcomma force

* Test 15.4: dpcomma + ignore + force
clear
set obs 20
gen str20 x = cond(mod(_n, 5) == 0, "NA", "€" + subinstr(string(_n * 100, "%10.2f"), ".", ",", 1))
benchmark_destring x, testname("dpcomma + ignore + force") generate(x) dpcomma ignore("€") force

* Test 15.5: percent + force
clear
set obs 20
gen str15 x = cond(mod(_n, 4) == 0, "N/A%", string(_n * 5) + "%")
benchmark_destring x, testname("percent + force") generate(x) percent force

* Test 15.6: float + force
clear
set obs 30
gen str20 x = cond(mod(_n, 6) == 0, "missing", string(_n / 7, "%10.4f"))
benchmark_destring x, testname("float + force") generate(x) float force

* Test 15.7: float + ignore
clear
set obs 30
gen str20 x = "$" + string(_n / 3, "%10.3f")
benchmark_destring x, testname("float + ignore") generate(x) float ignore("$")

* Test 15.8: float + percent
clear
set obs 30
gen str15 x = string(_n * 2) + "%"
benchmark_destring x, testname("float + percent") generate(x) float percent

* Test 15.9: Multiple ignore characters + force
clear
set obs 30
gen str30 x = cond(mod(_n, 5) == 0, "N/A", "$" + string(_n * 1000, "%15.0gc"))
benchmark_destring x, testname("ignore $, + force") generate(x) ignore("$,") force

* Test 15.10: All options combined (ignore + force + float)
clear
set obs 30
gen str30 x = cond(mod(_n, 7) == 0, "NA", "£" + string(_n / 5, "%10.2f"))
benchmark_destring x, testname("ignore + force + float") generate(x) ignore("£") force float

/*******************************************************************************
 * SECTION 16: Missing string patterns
 ******************************************************************************/
print_section "Missing String Patterns"

* Test 16.1: All empty strings
clear
set obs 20
gen str10 x = ""
cdestring x, generate(x_num) force
count if missing(x_num)
if r(N) == 20 {
    test_pass "all empty -> all missing"
}
else {
    test_fail "all empty" "not all missing"
}

* Test 16.2: Mixed empty and numeric
clear
set obs 20
gen str10 x = cond(mod(_n, 2) == 0, "", string(_n))
benchmark_destring x, testname("mixed empty and numeric") generate(x) force

* Test 16.3: First observation empty
clear
set obs 10
gen str10 x = string(_n)
replace x = "" in 1
benchmark_destring x, testname("first observation empty") generate(x) force

* Test 16.4: Last observation empty
clear
set obs 10
gen str10 x = string(_n)
replace x = "" in 10
benchmark_destring x, testname("last observation empty") generate(x) force

* Test 16.5: Consecutive empties
clear
set obs 20
gen str10 x = string(_n)
replace x = "" in 5/10
benchmark_destring x, testname("consecutive empties") generate(x) force

* Test 16.6: Only first non-empty
clear
set obs 20
gen str10 x = ""
replace x = "42" in 1
cdestring x, generate(x_num) force
if x_num[1] == 42 & missing(x_num[2]) {
    test_pass "only first non-empty"
}
else {
    test_fail "only first" "wrong pattern"
}

* Test 16.7: Only last non-empty
clear
set obs 20
gen str10 x = ""
replace x = "99" in 20
cdestring x, generate(x_num) force
if x_num[20] == 99 & missing(x_num[1]) {
    test_pass "only last non-empty"
}
else {
    test_fail "only last" "wrong pattern"
}

/*******************************************************************************
 * SECTION 17: Special number formats
 ******************************************************************************/
print_section "Special Number Formats"

* Test 17.1: Hexadecimal-looking strings (treated as non-numeric)
clear
set obs 5
gen str10 x = ""
replace x = "0x10" in 1
replace x = "0xFF" in 2
replace x = "1A2B" in 3
replace x = "CAFE" in 4
replace x = "100" in 5
benchmark_destring x, testname("hex-looking strings") generate(x) force

* Test 17.2: Infinity-looking strings
clear
set obs 5
gen str10 x = ""
replace x = "Inf" in 1
replace x = "-Inf" in 2
replace x = "inf" in 3
replace x = "+Inf" in 4
replace x = "100" in 5
benchmark_destring x, testname("infinity-looking") generate(x) force

* Test 17.3: NaN-looking strings
clear
set obs 5
gen str10 x = ""
replace x = "NaN" in 1
replace x = "nan" in 2
replace x = "NA" in 3
replace x = "N/A" in 4
replace x = "100" in 5
benchmark_destring x, testname("NaN-looking") generate(x) force

* Test 17.4: Currency formats
clear
set obs 10
gen str20 x = ""
replace x = "$100" in 1
replace x = "£200" in 2
replace x = "€300" in 3
replace x = "¥400" in 4
replace x = "500 USD" in 5
replace x = "600" in 6
replace x = "-$100" in 7
replace x = "$-100" in 8
replace x = "($100)" in 9
replace x = "100$" in 10
benchmark_destring x, testname("currency formats") generate(x) force

* Test 17.5: Accounting format (parentheses for negative)
clear
set obs 5
gen str20 x = ""
replace x = "(100)" in 1
replace x = "(200.50)" in 2
replace x = "300" in 3
replace x = "(0)" in 4
replace x = "()" in 5
benchmark_destring x, testname("accounting format") generate(x) ignore("()") force

* Test 17.6: Number with units
clear
set obs 10
gen str20 x = ""
replace x = "100kg" in 1
replace x = "200lb" in 2
replace x = "300m" in 3
replace x = "400km" in 4
replace x = "500ml" in 5
replace x = "600" in 6
replace x = "1.5kg" in 7
replace x = "-2.5lb" in 8
replace x = "0.001m" in 9
replace x = "1e3km" in 10
benchmark_destring x, testname("numbers with units") generate(x) ignore("kglbmm") force

/*******************************************************************************
 * SECTION 18: Verbose and threads options
 ******************************************************************************/
print_section "Verbose and Threads Options"

* Test 18.1: verbose verifies correctness
clear
set obs 100
gen str10 x = string(_n)
cdestring x, generate(x_c) verbose
destring x, generate(x_s)
count if x_c != x_s
if r(N) == 0 {
    test_pass "[syntax] verbose option verifies correctness"
}
else {
    test_fail "verbose" "`=r(N)' differ"
}
drop x_c x_s

* Test 18.2: threads(1) matches destring
clear
set obs 1000
gen str10 x = string(_n)
cdestring x, generate(x_c) threads(1)
destring x, generate(x_s)
count if x_c != x_s
if r(N) == 0 {
    test_pass "threads(1) matches destring"
}
else {
    test_fail "threads(1)" "`=r(N)' differ"
}
drop x_c x_s

* Test 18.3: threads(2) matches destring
clear
set obs 1000
gen str10 x = string(_n)
cdestring x, generate(x_c) threads(2)
destring x, generate(x_s)
count if x_c != x_s
if r(N) == 0 {
    test_pass "threads(2) matches destring"
}
else {
    test_fail "threads(2)" "`=r(N)' differ"
}
drop x_c x_s

* Test 18.4: threads(4) matches destring
clear
set obs 5000
gen str15 x = string(_n / 3, "%10.5f")
cdestring x, generate(x_c) threads(4)
destring x, generate(x_s)
count if abs(x_c - x_s) > $DEFAULT_TOL & !missing(x_c) & !missing(x_s)
if r(N) == 0 {
    test_pass "threads(4) matches destring"
}
else {
    test_fail "threads(4)" "`=r(N)' differ"
}
drop x_c x_s

* Test 18.5: verbose + threads
clear
set obs 500
gen str10 x = string(_n)
capture cdestring x, generate(x_num) verbose threads(2)
if _rc == 0 {
    count if x_num == _n
    if r(N) == 500 {
        test_pass "[syntax] verbose + threads"
    }
    else {
        test_fail "verbose+threads" "wrong values"
    }
}
else {
    test_fail "verbose+threads" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 19: if/in with various options
 ******************************************************************************/
print_section "if/in with Options"

* Test 19.1: if with force (destring doesn't support if/in, just verify cdestring works)
clear
set obs 50
gen str20 x = cond(mod(_n, 3) == 0, "NA", string(_n))
gen flag = mod(_n, 2)
cdestring x if flag == 1, generate(x_num) force
count if !missing(x_num) & flag == 1 & mod(_n, 3) != 0
local n_conv = r(N)
* Should convert numeric values where flag==1 and not "NA"
if `n_conv' > 0 {
    test_pass "if with force"
}
else {
    test_fail "if+force" "no conversions"
}

* Test 19.2: if with ignore
clear
set obs 50
gen str20 x = "$" + string(_n)
gen flag = mod(_n, 2)
cdestring x if flag == 1, generate(x_num) ignore("$")
count if !missing(x_num) & flag == 1
local n_conv = r(N)
count if flag == 1
if `n_conv' == r(N) {
    test_pass "if with ignore"
}
else {
    test_fail "if+ignore" "wrong count"
}

* Test 19.3: if with percent
clear
set obs 50
gen str15 x = string(_n * 2) + "%"
gen flag = mod(_n, 2)
cdestring x if flag == 1, generate(x_num) percent
count if !missing(x_num) & flag == 1
local n_conv = r(N)
count if flag == 1
if `n_conv' == r(N) {
    test_pass "if with percent"
}
else {
    test_fail "if+percent" "wrong count"
}

* Test 19.4: in with force
clear
set obs 100
gen str20 x = cond(mod(_n, 5) == 0, "NA", string(_n))
cdestring x in 1/50, generate(x_num) force
count if !missing(x_num) in 1/50
local n_conv = r(N)
* 50 observations, 10 are "NA" -> 40 should convert
if `n_conv' == 40 {
    test_pass "in with force"
}
else {
    test_fail "in+force" "expected 40, got `n_conv'"
}

* Test 19.5: in with dpcomma
clear
set obs 100
gen str20 x = subinstr(string(_n / 3, "%10.2f"), ".", ",", 1)
cdestring x in 20/70, generate(x_num) dpcomma
count if !missing(x_num) in 20/70
if r(N) == 51 {
    test_pass "in with dpcomma"
}
else {
    test_fail "in+dpcomma" "expected 51, got `=r(N)'"
}

* Test 19.6: Combined if and in with options
clear
set obs 100
gen str20 x = "$" + string(_n)
gen flag = mod(_n, 3) == 0
cdestring x if flag == 1 in 1/60, generate(x_num) ignore("$")
* Should convert multiples of 3 up to 60: 3,6,9,...,60 = 20 values
count if !missing(x_num)
if r(N) == 20 {
    test_pass "if + in + ignore"
}
else {
    test_fail "if+in+ignore" "expected 20, got `=r(N)'"
}

/*******************************************************************************
 * SECTION 20: String variable types
 ******************************************************************************/
print_section "String Variable Types"

* Test 20.1: str1 variable
clear
set obs 10
gen str1 x = string(mod(_n - 1, 10))
benchmark_destring x, testname("str1 variable") generate(x)

* Test 20.2: str10 variable
clear
set obs 100
gen str10 x = string(_n)
benchmark_destring x, testname("str10 variable") generate(x)

* Test 20.3: str50 variable
clear
set obs 50
gen str50 x = string(_n, "%50.10f")
benchmark_destring x, testname("str50 variable") generate(x)

* Test 20.4: str100 variable
clear
set obs 30
gen str100 x = string(_n / 7, "%100.20f")
benchmark_destring x, testname("str100 variable") generate(x)

* Test 20.5: str244 variable (max fixed string)
clear
set obs 20
gen str244 x = string(_n / 3, "%20.10f")
benchmark_destring x, testname("str244 variable") generate(x)

* Test 20.6: strL variable
clear
set obs 20
gen strL x = string(_n * 100)
capture cdestring x, generate(x_c)
local cdestring_rc = _rc
capture destring x, generate(x_s)
local destring_rc = _rc
if `cdestring_rc' == `destring_rc' {
    if `cdestring_rc' == 0 {
        count if x_c != x_s
        if r(N) == 0 {
            test_pass "strL variable"
        }
        else {
            test_fail "strL" "`=r(N)' differ"
        }
    }
    else {
        test_pass "strL - both error as expected"
    }
}
else {
    test_fail "strL" "cdestring rc=`cdestring_rc' destring rc=`destring_rc'"
}

/*******************************************************************************
 * SECTION 21: Replace option comprehensive
 ******************************************************************************/
print_section "Replace Option Comprehensive"

* Test 21.1: Replace with force
clear
set obs 30
gen str20 x = cond(mod(_n, 4) == 0, "NA", string(_n))
benchmark_destring x, testname("replace with force") replace force

* Test 21.2: Replace with ignore + force
clear
set obs 30
gen str20 x = "$" + cond(mod(_n, 5) == 0, "NA", string(_n))
benchmark_destring x, testname("replace ignore + force") replace ignore("$") force

* Test 21.3: Replace with percent
clear
set obs 30
gen str15 x = string(_n * 3) + "%"
benchmark_destring x, testname("replace with percent") replace percent

* Test 21.4: Replace with dpcomma
clear
set obs 30
gen str20 x = subinstr(string(_n / 5, "%10.3f"), ".", ",", 1)
benchmark_destring x, testname("replace with dpcomma") replace dpcomma

* Test 21.5: Replace with float
clear
set obs 30
gen str20 x = string(_n / 7, "%10.5f")
benchmark_destring x, testname("replace with float") replace float

* Test 21.6: Replace multiple variables with force
clear
set obs 30
gen str20 a = cond(mod(_n, 3) == 0, "NA", string(_n))
gen str20 b = cond(mod(_n, 4) == 0, "NA", string(_n * 2))
gen str20 c = cond(mod(_n, 5) == 0, "NA", string(_n * 3))
benchmark_destring a b c, testname("replace multi with force") replace force

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cdestring returns the same error codes as Stata's destring
 * when given invalid inputs or error conditions.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Variable doesn't exist
sysuse auto, clear
test_error_match, stata_cmd(destring nonexistent_var, generate(test)) ctools_cmd(cdestring nonexistent_var, generate(test)) testname("nonexistent variable")

* Already numeric variable
sysuse auto, clear
test_error_match, stata_cmd(destring price, generate(test)) ctools_cmd(cdestring price, generate(test)) testname("numeric variable input")

* Missing generate or replace option
clear
set obs 10
gen str10 x = string(_n)
test_error_match, stata_cmd(destring x) ctools_cmd(cdestring x) testname("missing generate/replace option")

* Generate variable already exists
clear
set obs 10
gen str10 x = string(_n)
gen y = _n
test_error_match, stata_cmd(destring x, generate(y)) ctools_cmd(cdestring x, generate(y)) testname("generate var already exists")

* Non-numeric string without force
clear
set obs 5
gen str10 x = "abc"
test_error_match, stata_cmd(destring x, generate(test)) ctools_cmd(cdestring x, generate(test)) testname("non-numeric string without force")

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
* End of cdestring validation
noi print_summary "cdestring"
}
