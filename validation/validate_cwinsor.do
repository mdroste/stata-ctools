/*******************************************************************************
 * validate_cwinsor.do
 *
 * Comprehensive validation tests for cwinsor vs gstats winsor
 * Tests all options: cuts, p/q, trim, by, suffix/prefix, replace
 *
 * VERIFICATION: All tests compare winsorized/trimmed values between
 * gstats winsor and cwinsor using significant figures comparison
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for cwinsor..."

/*******************************************************************************
 * Helper: Compare cwinsor vs gstats winsor
 ******************************************************************************/
capture program drop benchmark_winsor
program define benchmark_winsor
    syntax varlist(numeric), testname(string) [Cuts(numlist min=2 max=2) P(real -1) Q(real -1) TRIM BY(varlist) SUFfix(string) PREfix(string) cwinsoropts(string)]

    * Set default percentiles
    if `p' < 0 & `q' < 0 & "`cuts'" == "" {
        local p = 1
        local q = 99
    }

    * Build cuts string for gstats
    if "`cuts'" != "" {
        local p : word 1 of `cuts'
        local q : word 2 of `cuts'
    }

    * Build suffix for new variables - always use suffix for comparison
    local suf "_cw"
    if "`suffix'" != "" local suf "`suffix'"
    if "`prefix'" != "" local suf ""

    * Build options for cwinsor - always use suffix or prefix for comparison
    local copts ""
    if "`cuts'" != "" {
        local copts "cuts(`cuts')"
    }
    else {
        local copts "p(`p') q(`q')"
    }
    if "`trim'" != "" local copts "`copts' trim"
    if "`by'" != "" local copts "`copts' by(`by')"
    if "`suffix'" != "" {
        local copts "`copts' suffix(`suffix')"
    }
    else if "`prefix'" != "" {
        local copts "`copts' prefix(`prefix')"
    }
    else {
        * Default: use suffix(_cw) for comparison
        local copts "`copts' suffix(_cw)"
    }
    if "`cwinsoropts'" != "" local copts "`copts' `cwinsoropts'"

    * Build options for gstats winsor
    local gopts "cuts(`p' `q')"
    if "`trim'" != "" local gopts "`gopts' trim"
    if "`by'" != "" local gopts "`gopts' by(`by')"

    * For gstats, always use suffix to create new variable
    local gopts "`gopts' suffix(_gw)"

    * Run cwinsor
    foreach v of local varlist {
        capture drop `v'`suf'
        if "`prefix'" != "" {
            capture drop `prefix'`v'
        }
    }
    capture cwinsor `varlist', `copts'
    local rc_c = _rc

    * Run gstats winsor
    foreach v of local varlist {
        capture drop `v'_gw
    }
    capture gstats winsor `varlist', `gopts'
    local rc_g = _rc

    * Check both succeeded or both failed
    if `rc_c' != `rc_g' {
        test_fail "`testname'" "cwinsor rc=`rc_c', gstats rc=`rc_g'"
        exit
    }

    if `rc_c' != 0 {
        test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare results for each variable
    local all_pass = 1
    local fail_reason ""
    foreach v of local varlist {
        * Determine cwinsor output variable name
        local cvar "`v'`suf'"
        if "`prefix'" != "" local cvar "`prefix'`v'"

        * gstats output variable name
        local gvar "`v'_gw"

        * Compare using sigfigs
        quietly count if `cvar' != `gvar' & !missing(`cvar') & !missing(`gvar')
        local ndiff = r(N)

        if `ndiff' > 0 {
            * Use significant figures comparison for floating point tolerance
            tempvar sf
            quietly {
                gen double `sf' = 15 if `cvar' == `gvar'
                replace `sf' = 0 if (`sf' == .) & ((`cvar' == 0) | (`gvar' == 0))
                replace `sf' = -log10(abs(`cvar' - `gvar') / max(abs(`cvar'), abs(`gvar'))) if `sf' == .
                replace `sf' = 0 if `sf' < 0
                replace `sf' = 15 if `sf' > 15
            }

            quietly count if `sf' < $DEFAULT_SIGFIGS & !missing(`cvar') & !missing(`gvar')
            local nfail = r(N)
            drop `sf'

            if `nfail' > 0 {
                local all_pass = 0
                local fail_reason "`v': `nfail' values differ beyond tolerance"
            }
        }

        * Also check missing value patterns match
        quietly count if missing(`cvar') != missing(`gvar')
        local nmiss_diff = r(N)
        if `nmiss_diff' > 0 {
            local all_pass = 0
            local fail_reason "`v': `nmiss_diff' missing value patterns differ"
        }
    }

    if `all_pass' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "`fail_reason'"
    }

    * Cleanup (use capture to ignore if variable doesn't exist)
    foreach v of local varlist {
        capture drop `v'`suf'
        if "`prefix'" != "" {
            capture drop `prefix'`v'
        }
        capture drop `v'_gw
    }
end

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
print_section "Plugin Check"

sysuse auto, clear
capture noisily cwinsor price, suffix(_w)
if _rc != 0 {
    test_fail "cwinsor plugin load" "plugin returned error `=_rc'"
    print_summary "cwinsor"
    exit 1
}
test_pass "cwinsor plugin loads and runs"
drop price_w

/*******************************************************************************
 * SECTION 2: Basic functionality
 ******************************************************************************/
print_section "Basic Functionality"

* Test 2.1: Basic winsorization creates variable
sysuse auto, clear
capture drop price_w
cwinsor price, suffix(_w)
capture confirm numeric variable price_w
if _rc == 0 {
    test_pass "basic winsorization creates variable"
}
else {
    test_fail "basic winsorization" "variable not created"
}

* Test 2.2: Compare with gstats winsor (default 1/99)
sysuse auto, clear
benchmark_winsor price, testname("vs gstats: price default cuts")

* Test 2.3: Multiple variables
sysuse auto, clear
benchmark_winsor price mpg weight, testname("vs gstats: multiple variables")

* Test 2.4: Variable with missing values
clear
set obs 100
gen x = rnormal()
replace x = . in 1/10
benchmark_winsor x, testname("vs gstats: variable with missings")

* Test 2.5: All missing values
clear
set obs 50
gen x = .
capture cwinsor x, suffix(_w)
if _rc != 0 {
    test_pass "all missing values handled (error expected)"
}
else {
    * Check result is all missing
    quietly count if !missing(x_w)
    if r(N) == 0 {
        test_pass "all missing values handled"
    }
    else {
        test_fail "all missing values" "expected all missing output"
    }
    capture drop x_w
}

/*******************************************************************************
 * SECTION 3: Percentile options
 ******************************************************************************/
print_section "Percentile Options"

* Test 3.1: cuts(5 95)
sysuse auto, clear
benchmark_winsor price, testname("cuts(5 95)") cuts(5 95)

* Test 3.2: cuts(10 90)
sysuse auto, clear
benchmark_winsor price, testname("cuts(10 90)") cuts(10 90)

* Test 3.3: cuts(0.5 99.5)
sysuse auto, clear
benchmark_winsor price, testname("cuts(0.5 99.5)") cuts(0.5 99.5)

* Test 3.4: p() q() syntax
sysuse auto, clear
cwinsor price, p(5) q(95) suffix(_c)
gstats winsor price, cuts(5 95) suffix(_g)
quietly count if price_c != price_g & !missing(price_c) & !missing(price_g)
if r(N) == 0 {
    test_pass "p(5) q(95) syntax"
}
else {
    test_fail "p(5) q(95) syntax" "`=r(N)' values differ"
}
drop price_c price_g

* Test 3.5: Asymmetric percentiles
sysuse auto, clear
benchmark_winsor price, testname("cuts(2 98)") cuts(2 98)

* Test 3.6: Wide percentiles cuts(25 75)
sysuse auto, clear
benchmark_winsor price, testname("cuts(25 75)") cuts(25 75)

* Test 3.7: Narrow percentiles cuts(0.1 99.9)
sysuse auto, clear
benchmark_winsor price, testname("cuts(0.1 99.9)") cuts(0.1 99.9)

/*******************************************************************************
 * SECTION 4: Trim mode
 ******************************************************************************/
print_section "Trim Mode"

* Test 4.1: Basic trim
sysuse auto, clear
benchmark_winsor price, testname("trim mode basic") trim

* Test 4.2: Trim with different percentiles
sysuse auto, clear
benchmark_winsor price, testname("trim cuts(5 95)") trim cuts(5 95)

* Test 4.3: Trim with multiple variables
sysuse auto, clear
benchmark_winsor price mpg, testname("trim multiple vars") trim

* Test 4.4: Trim with missing values
clear
set obs 100
gen x = rnormal()
replace x = . in 1/10
benchmark_winsor x, testname("trim with missings") trim

* Test 4.5: Verify trim actually replaces with missing
sysuse auto, clear
cwinsor price, suffix(_w) trim cuts(5 95)
gstats winsor price, suffix(_g) trim cuts(5 95)
quietly count if missing(price_w) != missing(price_g)
if r(N) == 0 {
    test_pass "trim missing pattern matches gstats"
}
else {
    test_fail "trim missing pattern" "`=r(N)' obs differ in missing status"
}
drop price_w price_g

/*******************************************************************************
 * SECTION 5: By-group processing
 ******************************************************************************/
print_section "By-Group Processing"

* Test 5.1: Single by-variable
sysuse auto, clear
benchmark_winsor price, testname("by(foreign)") by(foreign)

* Test 5.2: Multiple by-variables
sysuse auto, clear
benchmark_winsor price, testname("by(foreign rep78)") by(foreign rep78)

* Test 5.3: By with multiple target variables
sysuse auto, clear
benchmark_winsor price mpg, testname("by(foreign) multiple vars") by(foreign)

* Test 5.4: By with trim
sysuse auto, clear
benchmark_winsor price, testname("by(foreign) trim") by(foreign) trim

* Test 5.5: By with different percentiles
sysuse auto, clear
benchmark_winsor price, testname("by(foreign) cuts(10 90)") by(foreign) cuts(10 90)

* Test 5.6: By-variable with missing values
sysuse auto, clear
replace rep78 = . if _n <= 5
benchmark_winsor price, testname("by(rep78) with missing by-var") by(rep78)

* Test 5.7: Many groups
clear
set obs 10000
gen group = mod(_n, 100)
gen x = rnormal()
benchmark_winsor x, testname("by() with 100 groups") by(group)

* Test 5.8: Single observation per group
clear
set obs 50
gen group = _n
gen x = rnormal()
benchmark_winsor x, testname("by() singleton groups") by(group)

/*******************************************************************************
 * SECTION 6: Suffix and prefix options
 ******************************************************************************/
print_section "Suffix/Prefix Options"

* Test 6.1: Custom suffix
sysuse auto, clear
cwinsor price, suffix(_winsor)
capture confirm variable price_winsor
if _rc == 0 {
    test_pass "custom suffix(_winsor)"
}
else {
    test_fail "custom suffix" "variable not created"
}
drop price_winsor

* Test 6.2: Prefix option
sysuse auto, clear
cwinsor price, prefix(w_)
capture confirm variable w_price
if _rc == 0 {
    test_pass "prefix(w_)"
}
else {
    test_fail "prefix option" "variable not created"
}
drop w_price

* Test 6.3: Suffix with multiple variables
sysuse auto, clear
cwinsor price mpg, suffix(_w)
capture confirm variable price_w
local rc1 = _rc
capture confirm variable mpg_w
local rc2 = _rc
if `rc1' == 0 & `rc2' == 0 {
    test_pass "suffix with multiple variables"
}
else {
    test_fail "suffix multiple vars" "not all variables created"
}
drop price_w mpg_w

/*******************************************************************************
 * SECTION 7: Replace option
 ******************************************************************************/
print_section "Replace Option"

* Test 7.1: Basic replace
sysuse auto, clear
gen price_orig = price
cwinsor price, replace
quietly count if price != price_orig
if r(N) > 0 {
    test_pass "replace modifies original variable"
}
else {
    test_fail "replace option" "variable not modified"
}

* Test 7.2: Replace matches suffix result
sysuse auto, clear
gen price_orig = price
cwinsor price, suffix(_w)
cwinsor price_orig, replace
quietly count if price_w != price_orig
if r(N) == 0 {
    test_pass "replace matches suffix result"
}
else {
    test_fail "replace matches suffix" "`=r(N)' values differ"
}

/*******************************************************************************
 * SECTION 8: Edge cases
 ******************************************************************************/
print_section "Edge Cases"

* Test 8.1: Small dataset
clear
set obs 5
gen x = _n
benchmark_winsor x, testname("small dataset (5 obs)")

* Test 8.2: Large dataset
clear
set obs 100000
gen x = rnormal()
benchmark_winsor x, testname("large dataset (100k obs)")

* Test 8.3: Constant values
clear
set obs 100
gen x = 5
benchmark_winsor x, testname("constant values")

* Test 8.4: Extreme values
clear
set obs 100
gen x = rnormal()
replace x = 1e10 in 1
replace x = -1e10 in 2
benchmark_winsor x, testname("extreme values")

* Test 8.5: Integer variable
sysuse auto, clear
benchmark_winsor rep78, testname("integer variable")

* Test 8.6: Byte variable
clear
set obs 100
gen byte x = mod(_n, 10)
benchmark_winsor x, testname("byte variable")

* Test 8.7: Float variable
clear
set obs 100
gen float x = rnormal()
benchmark_winsor x, testname("float variable")

* Test 8.8: Double variable
clear
set obs 100
gen double x = rnormal() * 1e10
benchmark_winsor x, testname("double variable")

/*******************************************************************************
 * SECTION 9: If/in conditions
 ******************************************************************************/
print_section "If/In Conditions"

* Test 9.1: if condition
sysuse auto, clear
cwinsor price if foreign == 1, suffix(_c)
gstats winsor price if foreign == 1, suffix(_g)
quietly count if price_c != price_g & !missing(price_c) & !missing(price_g) & foreign == 1
if r(N) == 0 {
    test_pass "if condition"
}
else {
    test_fail "if condition" "`=r(N)' values differ"
}
drop price_c price_g

* Test 9.2: in range
sysuse auto, clear
cwinsor price in 1/50, suffix(_c)
gstats winsor price in 1/50, suffix(_g)
quietly count if price_c != price_g & !missing(price_c) & !missing(price_g) in 1/50
if r(N) == 0 {
    test_pass "in range"
}
else {
    test_fail "in range" "`=r(N)' values differ"
}
drop price_c price_g

* Test 9.3: if and in combined
sysuse auto, clear
cwinsor price if mpg > 20 in 1/50, suffix(_c)
gstats winsor price if mpg > 20 in 1/50, suffix(_g)
quietly count if price_c != price_g & !missing(price_c) & !missing(price_g)
if r(N) == 0 {
    test_pass "if and in combined"
}
else {
    test_fail "if and in combined" "`=r(N)' values differ"
}
drop price_c price_g

/*******************************************************************************
 * SECTION 10: Error handling
 ******************************************************************************/
print_section "Error Handling"

* Test 10.1: Invalid percentile (p > q)
sysuse auto, clear
capture cwinsor price, p(99) q(1) suffix(_w)
if _rc != 0 {
    test_pass "error: p > q"
}
else {
    test_fail "error: p > q" "should have failed"
    capture drop price_w
}

* Test 10.2: Percentile out of range (negative)
sysuse auto, clear
capture cwinsor price, p(-1) q(99) suffix(_w)
if _rc != 0 {
    test_pass "error: negative percentile"
}
else {
    test_fail "error: negative percentile" "should have failed"
    capture drop price_w
}

* Test 10.3: Percentile out of range (> 100)
sysuse auto, clear
capture cwinsor price, p(1) q(101) suffix(_w)
if _rc != 0 {
    test_pass "error: percentile > 100"
}
else {
    test_fail "error: percentile > 100" "should have failed"
    capture drop price_w
}

* Test 10.4: String variable
sysuse auto, clear
capture cwinsor make, suffix(_w)
if _rc != 0 {
    test_pass "error: string variable"
}
else {
    test_fail "error: string variable" "should have failed"
    capture drop make_w
}

* Test 10.5: Nonexistent variable
sysuse auto, clear
capture cwinsor nonexistent, suffix(_w)
if _rc != 0 {
    test_pass "error: nonexistent variable"
}
else {
    test_fail "error: nonexistent variable" "should have failed"
}

* Test 10.6: Both suffix and prefix
sysuse auto, clear
capture cwinsor price, suffix(_s) prefix(p_)
if _rc != 0 {
    test_pass "error: both suffix and prefix"
}
else {
    test_fail "error: both suffix and prefix" "should have failed"
    capture drop price_s p_price
}

* Test 10.7: Replace with suffix
sysuse auto, clear
capture cwinsor price, replace suffix(_w)
if _rc != 0 {
    test_pass "error: replace with suffix"
}
else {
    test_fail "error: replace with suffix" "should have failed"
    capture drop price_w
}

* Test 10.8: Variable already exists
sysuse auto, clear
gen price_w = 0
capture cwinsor price, suffix(_w)
if _rc != 0 {
    test_pass "error: variable already exists"
}
else {
    test_fail "error: variable already exists" "should have failed"
}
drop price_w

/*******************************************************************************
 * SECTION 11: Verbose and threads options
 ******************************************************************************/
print_section "Options"

* Test 11.1: Verbose option produces output
sysuse auto, clear
cwinsor price, suffix(_w) verbose
test_pass "[syntax] verbose option"
drop price_w

* Test 11.2: Threads option (correctness check)
sysuse auto, clear
cwinsor price, suffix(_c1) threads(1)
cwinsor price, suffix(_c4) threads(4)
quietly count if price_c1 != price_c4
if r(N) == 0 {
    test_pass "threads(1) vs threads(4) identical results"
}
else {
    test_fail "threads option" "`=r(N)' values differ between thread counts"
}
drop price_c1 price_c4

/*******************************************************************************
 * SECTION 12: Numerical precision
 ******************************************************************************/
print_section "Numerical Precision"

* Test 12.1: Very small values
clear
set obs 100
gen x = rnormal() * 1e-10
benchmark_winsor x, testname("very small values (1e-10)")

* Test 12.2: Very large values
clear
set obs 100
gen x = rnormal() * 1e10
benchmark_winsor x, testname("very large values (1e10)")

* Test 12.3: Mixed magnitude
clear
set obs 100
gen x = rnormal()
replace x = x * 1e10 if _n <= 25
replace x = x * 1e-10 if _n > 75
benchmark_winsor x, testname("mixed magnitude values")

} /* end quietly */

print_summary "cwinsor"
