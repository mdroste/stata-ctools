* cwinsor: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cwinsor"

if "`__case'" == "cwinsor" {
* validate_cwinsor.do
* Comprehensive validation tests for cwinsor vs winsor2 (documented replacement)
* Tests all options: cuts, p/q, trim, by, suffix/prefix, replace
* VERIFICATION: All tests compare winsorized/trimmed values between
* winsor2 and cwinsor using significant figures comparison. Direct gstats
* checks below are supplementary checks on fixtures where the references agree.

* Load setup (works from project root or validation dir)
do "validation/support/validate_setup.do"

quietly {

noi di as text "Running validation tests for cwinsor..."

* SECTION 1: Plugin functionality check
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

* SECTION 2: Basic functionality
print_section "Basic Functionality"

* Test 2.1: Basic winsorization creates variable
sysuse auto, clear
capture drop price_w
cwinsor price, suffix(_w)
capture confirm numeric variable price_w
if _rc == 0 test_pass "basic winsorization creates variable"
else test_fail "basic winsorization" "variable not created"

* Test 2.2: Compare with winsor2 (default 1/99)
sysuse auto, clear
benchmark_winsor price, testname("vs winsor2: price default cuts")

* Test 2.3: Multiple variables
sysuse auto, clear
benchmark_winsor price mpg weight, testname("vs winsor2: multiple variables")

* Test 2.4: Variable with missing values
clear
set obs 100
gen x = rnormal()
replace x = . in 1/10
benchmark_winsor x, testname("vs winsor2: variable with missings")

* Test 2.5: All missing values
clear
set obs 50
gen x = .
capture cwinsor x, suffix(_w)
local ctools_rc = _rc
capture gstats winsor x, suffix(_g) cuts(1 99)
local gstats_rc = _rc
if `ctools_rc' == `gstats_rc' test_pass "all missing values handled (both rc=`ctools_rc')"
else test_fail "all missing values" "rc differ: gstats=`gstats_rc' cwinsor=`ctools_rc'"
capture drop x_w x_g

* SECTION 3: Percentile options
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
assert_var_equal price_c price_g $DEFAULT_SIGFIGS "p(5) q(95) syntax"
capture drop price_c price_g

* Test 3.5: Asymmetric percentiles
sysuse auto, clear
benchmark_winsor price, testname("cuts(2 98)") cuts(2 98)

* Test 3.6: Wide percentiles cuts(25 75)
sysuse auto, clear
benchmark_winsor price, testname("cuts(25 75)") cuts(25 75)

* Test 3.7: Narrow percentiles cuts(0.1 99.9)
sysuse auto, clear
benchmark_winsor price, testname("cuts(0.1 99.9)") cuts(0.1 99.9)

* SECTION 4: Trim mode
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
if r(N) == 0 test_pass "trim missing pattern matches gstats"
else test_fail "trim missing pattern" "`=r(N)' obs differ in missing status"
drop price_w price_g

* SECTION 5: By-group processing
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

* SECTION 6: Suffix and prefix options
print_section "Suffix/Prefix Options"

* Test 6.1: Custom suffix
sysuse auto, clear
cwinsor price, suffix(_winsor)
capture confirm variable price_winsor
if _rc == 0 test_pass "custom suffix(_winsor)"
else test_fail "custom suffix" "variable not created"
drop price_winsor

* Test 6.2: Prefix option
sysuse auto, clear
cwinsor price, prefix(w_)
capture confirm variable w_price
if _rc == 0 test_pass "prefix(w_)"
else test_fail "prefix option" "variable not created"
drop w_price

* Test 6.3: Suffix with multiple variables
sysuse auto, clear
cwinsor price mpg, suffix(_w)
capture confirm variable price_w
local rc1 = _rc
capture confirm variable mpg_w
local rc2 = _rc
if `rc1' == 0 & `rc2' == 0 test_pass "suffix with multiple variables"
else test_fail "suffix multiple vars" "not all variables created"
drop price_w mpg_w

* SECTION 7: Replace option
print_section "Replace Option"

* Test 7.1: Basic replace (use cuts(5 95) so values actually get winsorized)
* Compare against gstats winsor
sysuse auto, clear
gen price_orig = price
gstats winsor price, replace cuts(5 95)
rename price price_gstats
gen price = price_orig
cwinsor price, replace cuts(5 95)
assert_var_equal price price_gstats $DEFAULT_SIGFIGS "replace matches gstats winsor"

* Test 7.2: Replace matches suffix result
sysuse auto, clear
gen price_orig = price
cwinsor price, suffix(_w) cuts(5 95)
cwinsor price_orig, replace cuts(5 95)
quietly count if price_w != price_orig
if r(N) == 0 test_pass "replace matches suffix result"
else test_fail "replace matches suffix" "`=r(N)' values differ"

* SECTION 8: Edge cases
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

* SECTION 9: If/in conditions
print_section "If/In Conditions"

* Test 9.1: if condition
sysuse auto, clear
cwinsor price if foreign == 1, suffix(_c)
gstats winsor price if foreign == 1, suffix(_g)
assert_var_equal price_c price_g $DEFAULT_SIGFIGS "if condition"
capture drop price_c price_g

* Test 9.2: in range
sysuse auto, clear
cwinsor price in 1/50, suffix(_c)
gstats winsor price in 1/50, suffix(_g)
assert_var_equal price_c price_g $DEFAULT_SIGFIGS "in range"
capture drop price_c price_g

* Test 9.3: if and in combined
sysuse auto, clear
cwinsor price if mpg > 20 in 1/50, suffix(_c)
gstats winsor price if mpg > 20 in 1/50, suffix(_g)
assert_var_equal price_c price_g $DEFAULT_SIGFIGS "if and in combined"
capture drop price_c price_g

* SECTION 10: Error handling
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
if _rc != 0 test_pass "error: nonexistent variable"
else test_fail "error: nonexistent variable" "should have failed"

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
if _rc != 0 test_pass "error: variable already exists"
else test_fail "error: variable already exists" "should have failed"
drop price_w

* SECTION 11: Verbose and threads options
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
if r(N) == 0 test_pass "threads(1) vs threads(4) identical results"
else test_fail "threads option" "`=r(N)' values differ between thread counts"
drop price_c1 price_c4

* SECTION 12: Numerical precision
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

* SECTION 13: By-group identity
print_section "By-group Identity"

* Test 13.1: String by-variable
sysuse auto, clear
gen str8 grp = cond(foreign, "foreign", "domestic")
benchmark_winsor price, testname("string by-variable") by(grp) cuts(10 90)

* Test 13.2: . and .a are separate by-groups, interleaved in data order
sysuse auto, clear
gen double g = cond(mod(_n, 3) == 0, ., cond(mod(_n, 3) == 1, .a, 1))
benchmark_winsor price, testname("by-groups . and .a kept distinct") by(g) cuts(10 90)

* winsor2's macro-rounded bounds can trim an observation at a percentile.
* Compare exact missingness and metadata, with no boundary tolerance.
quietly do validation/support/coverage_helpers.do
foreach scale in 1 1e-6 {
    clear
    set obs 8
    gen double x=.
    local i=0
    foreach value in -4 -1.8189731519126051 -1 0 1 4.387736506821946 9 .a {
        local ++i
        replace x=`value' in `i'
    }
    replace x=x*`scale' if !missing(x)
    label variable x "Boundary fixture"
    foreach output in "suffix(_w)" "replace" {
        coverage_data, reference(winsor2 x, cuts(20 80) trim `output') ///
            candidate(cwinsor x, cuts(20 80) trim `output') ///
            name("winsor2 rounded trim boundary scale `scale' `output'")
    }
    gen byte g=mod(_n,2)
    coverage_data, reference(winsor2 x, cuts(20 80) by(g) suffix(_w)) ///
        candidate(cwinsor x, cuts(20 80) by(g) suffix(_w)) ///
        name("winsor2 grouped output format scale `scale'")
}

* Rollback of failed calls
* A failed or interrupted cwinsor leaves the data as they were without a
* snapshot of the dataset; option snapshot preserves it instead and option
* fast undoes nothing (see validation/support/rollback_helpers.do). In place, the
* plugin keeps the original values of the observations it rewrites. The
* hook $CTOOLS_cwinsor_fault fails the call:
* outputs  with r(1) after suffix()/prefix() variables were created
* store    inside the plugin after half the observations were written,
* with the memory error it simulates, r(920)
* post     with r(1) after the targets were written, relabeled and recast
print_section "Rollback of Failed Calls"
capture confirm file "validation/support/rollback_helpers.do"
if !_rc do "validation/support/rollback_helpers.do"
else do "rollback_helpers.do"

clear
set seed 1732
set obs 3000
gen long id = _n
gen int f2 = mod(_n * 7, 23) + 1
gen double y = rnormal() * 10
gen float x1 = rnormal()
gen int x2 = round(rnormal() * 100)
replace y = .a in 5
replace x1 = . in 7
gen byte keep = mod(_n, 9) != 0
gen byte cat = mod(_n, 3)
label define rb_cat 0 "zero" 1 "one" 2 "two"
label values cat rb_cat
label variable y "outcome"
label data "cwinsor rollback data"
note: dataset note
note y: note on y
char _dta[rb_char] dataset characteristic
char x1[rb_char] x1 characteristic
sort x2 id
tempfile rb_copy rb_out
save `rb_copy'

_ctools_rb_case, cmd(cwinsor y x1, cuts(5 95)) ///
    global(CTOOLS_cwinsor_fault) fault(post) testname("[rollback] in place after plugin")
_ctools_rb_case, cmd(cwinsor y x1, cuts(5 95)) ///
    global(CTOOLS_cwinsor_fault) fault(bounds) modes(default snapshot) testname("[rollback] in place between cutoff phases")
_ctools_rb_case, cmd(cwinsor y x1, cuts(5 95) suffix(_w)) ///
    global(CTOOLS_cwinsor_fault) fault(bounds) testname("[rollback] new outputs between cutoff phases")
_ctools_rb_case, cmd(cwinsor y x1, cuts(5 95)) ///
    global(CTOOLS_cwinsor_fault) fault(store) rc(920) testname("[rollback] multiple variables, deferred store fails")
_ctools_rb_case, cmd(cwinsor y x2 if keep, cuts(10 90) by(f2)) ///
    global(CTOOLS_cwinsor_fault) fault(post) testname("[rollback] in place by() if, sort variable")
_ctools_rb_case, cmd(cwinsor y x1 if keep, cuts(10 90) by(f2) replace) ///
    global(CTOOLS_cwinsor_fault) fault(store) rc(920) testname("[rollback] in place, store fails")
_ctools_rb_case, cmd(cwinsor x2, cuts(10 90)) ///
    global(CTOOLS_cwinsor_fault) fault(store) rc(920) testname("[rollback] one variable, store fails")
_ctools_rb_case, cmd(cwinsor y x2, cuts(5 95) trim) ///
    global(CTOOLS_cwinsor_fault) fault(post) testname("[rollback] in place trim after replace")
_ctools_rb_case, cmd(cwinsor y x2, cuts(5 95) suffix(_w)) ///
    global(CTOOLS_cwinsor_fault) fault(outputs) testname("[rollback] suffix outputs created")
_ctools_rb_case, cmd(cwinsor y x1 x2, cuts(5 95) by(f2) trim prefix(w_)) ///
    global(CTOOLS_cwinsor_fault) fault(post) testname("[rollback] prefix after recast")

* Options snapshot and fast give the default's results; together they are refused
foreach spec in "cuts(5 95) by(f2)" "cuts(5 95) trim" "cuts(5 95) suffix(_w)" {
    use `rb_copy', clear
    cwinsor y x1 x2, `spec'
    save `rb_out', replace
    foreach mode in snapshot fast {
        use `rb_copy', clear
        cwinsor y x1 x2, `spec' `mode'
        capture cf _all using `rb_out', all
        if !_rc test_pass "[rollback] option `mode' gives the default results (`spec')"
        else test_fail "[rollback] option `mode' gives the default results (`spec')" "results differ"
    }
}
use `rb_copy', clear
capture cwinsor y, cuts(5 95) snapshot fast
local rc = _rc
capture cf _all using `rb_copy', all
if `rc' == 198 & !_rc test_pass "[rollback] snapshot and fast together refused"
else test_fail "[rollback] snapshot and fast together refused" "rc=`rc', data changed=`=_rc!=0'"

} /* end quietly */

print_summary "cwinsor"

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "cwinsor"
    global CTOOLS_CASE_COMPLETE "cwinsor"
    exit
}

if "`__case'" == "winsor_rounding" {
* Exact winsor2 cutoff expansion, output storage, selection, and phase cleanup.
version 16
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
local cases=0
foreach storage in float double {
    foreach sign in -1 1 {
        foreach scale in 1 1e-20 1e20 {
            clear
            set obs 204
            gen `storage' x=0
            replace x=-3 in 1
            replace x=-2.3206636905670166 in 2
            replace x=-2.3202905654907227 in 3
            replace x=`sign'*`scale'*x
            replace x=. in 201
            replace x=.a in 202
            replace x=.z in 203
            replace x=.b in 204
            gen double y=1.2345678901234567*x
            gen long id=_n
            foreach sample in "" "if id<=200 in 1/203" {
                foreach trim in "" "trim" {
                    foreach cuts in "1 99" "0 99" "1 100" {
                        quietly winsor2 x y `sample', cuts(`cuts') `trim' suffix(_expected)
                        foreach threads in 1 8 {
                            quietly cwinsor x y `sample', cuts(`cuts') `trim' suffix(_actual) threads(`threads')
                            assert x_actual==x_expected & y_actual==y_expected
                            assert "`: type x_actual'"=="`: type x_expected'"
                            drop x_actual y_actual
                            local ++cases
                        }
                        drop x_expected y_expected
                        preserve
                        quietly winsor2 x y `sample', cuts(`cuts') `trim' replace
                        rename x x_expected
                        rename y y_expected
                        tempfile expected
                        quietly save `expected', replace
                        restore
                        preserve
                        quietly cwinsor x y `sample', cuts(`cuts') `trim' replace
                        quietly merge 1:1 id using `expected', nogen
                        assert x==x_expected & y==y_expected
                        restore
                        local ++cases
                    }
                }
            }
        }
    }
}
* A single by-group must use grouped semantics, without numeric macro rounding.
clear
set obs 200
gen float x=0
replace x=-3 in 1
replace x=-2.3206636905670166 in 2
replace x=-2.3202905654907227 in 3
gen g=1
winsor2 x, by(g) suffix(_expected)
cwinsor x, by(g) suffix(_actual)
assert x_actual==x_expected
local ++cases
noi di "WINSOR_ROUNDING_CASES=`cases'"
noi di "WINSOR_ROUNDING_COMPLETE"
    global CTOOLS_CASE_COMPLETE "winsor_rounding"
    exit
}

if "`__case'" == "winsor_transport" {
* Combined target/by loading: index spaces, selection and exact grouping.
capture do "validation/support/validate_setup.do"
if _rc do "validate_setup.do"
quietly {
    foreach n in 1 2 3 1000 20000 {
        clear
        set seed 981723
        set obs `n'
        gen long id = _n
        gen double x = rnormal()*50
        gen double y = rnormal()*2
        gen double z = mod(_n,13)
        replace x = .a if mod(_n,17)==0
        replace y = .z if mod(_n,31)==0
        gen double g = mod(_n,7) + 1e-9*mod(_n,3)
        replace g = .b if mod(_n,23)==0
        gen str20 s = cond(mod(_n,3)==0,"",cond(mod(_n,2),"prefix_a","prefix_b"))
        gen double shuffle=runiform()
        sort shuffle
        drop shuffle
        foreach keys in none g s gs {
            local byopt
            if "`keys'" == "g" local byopt by(g)
            if "`keys'" == "s" local byopt by(s)
            if "`keys'" == "gs" local byopt by(g s)
            foreach threads in 1 8 {
                foreach mode in winsor trim {
                    local trimopt
                    if "`mode'" == "trim" local trimopt trim
                    foreach count in 1 3 {
                        local vars x
                        if `count' == 3 local vars x y z
                        local sample
                        if `n' > 3 local sample if mod(id,5)<2 in 2/`n'
                        cwinsor `vars' `sample', `byopt' `trimopt' cuts(5 95) suffix(_got) threads(`threads')
                        * Ungrouped trimming follows winsor2, which rounds the
                        * bounds through macros (see validate_cwinsor.do), so
                        * gstats can keep an observed extreme that both trim.
                        if "`mode'" == "trim" & "`keys'" == "none" {
                            winsor2 `vars' `sample', trim cuts(5 95) suffix(_ref)
                        }
                        else gstats winsor `vars' `sample', `byopt' `trimopt' cuts(5 95) suffix(_ref)
                        foreach v of local vars {
                            assert missing(`v'_got) == missing(`v'_ref)
                            assert `v'_got == `v'_ref | abs(`v'_got-`v'_ref)/max(abs(`v'_got),abs(`v'_ref)) < $DEFAULT_TOL if !missing(`v'_got)
                        }
                        drop *_got *_ref
                        test_pass "winsor combined transport matches gstats"
                    }
                }
            }
        }
    }
}
print_summary "winsor_transport"
if $TESTS_FAILED exit 9
di "WINSOR_TRANSPORT_COMPLETE"

global CTOOLS_COMPONENT_COMPLETE "winsor_transport"
    global CTOOLS_CASE_COMPLETE "winsor_transport"
    exit
}

di as error "Unknown cwinsor case: `__case'"
exit 198
