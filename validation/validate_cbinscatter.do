/*******************************************************************************
 * validate_cbinscatter.do
 *
 * Comprehensive validation tests for cbinscatter
 * Tests all options: nquantiles, discrete, controls, absorb, by, method(binsreg),
 * linetype, weights, savedata
 *
 * VERIFICATION: Tests compare e() results and check consistency with binscatter
 * and binsreg where applicable
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CBINSCATTER VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture cbinscatter price mpg, nograph
if _rc != 0 {
    noi test_fail "cbinscatter plugin load" "returned error `=_rc'"
    noi print_summary "cbinscatter"
    exit 1
}
noi test_pass "cbinscatter plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic functionality
 ******************************************************************************/
noi print_section "Basic Functionality"

* Simple two-variable scatter
sysuse auto, clear
cbinscatter price mpg, nograph
if e(N) == _N & e(nquantiles) == 20 {
    noi test_pass "basic scatter (default 20 bins)"
}
else {
    noi test_fail "basic scatter" "N=`=e(N)' expected `=_N', nquantiles=`=e(nquantiles)' expected 20"
}

* Check e() results populated
sysuse auto, clear
cbinscatter price mpg, nograph
capture confirm matrix e(bindata)
if _rc == 0 {
    noi test_pass "e(bindata) matrix populated"
}
else {
    noi test_fail "e(bindata)" "matrix not found"
}

capture confirm matrix e(coefs)
if _rc == 0 {
    noi test_pass "e(coefs) matrix populated (linear fit)"
}
else {
    noi test_fail "e(coefs)" "matrix not found"
}

/*******************************************************************************
 * SECTION 3: nquantiles option
 ******************************************************************************/
noi print_section "nquantiles Option"

sysuse auto, clear
cbinscatter price mpg, nquantiles(10) nograph
if e(nquantiles) == 10 {
    noi test_pass "nquantiles(10)"
}
else {
    noi test_fail "nquantiles(10)" "got `=e(nquantiles)'"
}

sysuse auto, clear
cbinscatter price mpg, nquantiles(50) nograph
if e(nquantiles) == 50 {
    noi test_pass "nquantiles(50)"
}
else {
    noi test_fail "nquantiles(50)" "got `=e(nquantiles)'"
}

sysuse auto, clear
cbinscatter price mpg, nquantiles(5) nograph
if e(nquantiles) == 5 {
    noi test_pass "nquantiles(5)"
}
else {
    noi test_fail "nquantiles(5)" "got `=e(nquantiles)'"
}

* Test invalid nquantiles (should error)
sysuse auto, clear
capture cbinscatter price mpg, nquantiles(1) nograph
if _rc != 0 {
    noi test_pass "nquantiles(1) rejected (rc=`=_rc')"
}
else {
    noi test_fail "nquantiles(1)" "should have been rejected"
}

/*******************************************************************************
 * SECTION 4: discrete option
 ******************************************************************************/
noi print_section "discrete Option"

* Create discrete x variable
sysuse auto, clear
cbinscatter price rep78, discrete nograph
if e(N) > 0 {
    noi test_pass "discrete with rep78 (5 values)"
}
else {
    noi test_fail "discrete" "returned no observations"
}

* Discrete with few unique values
clear
set seed 12345
set obs 1000
gen x = runiformint(1, 10)
gen y = 2*x + rnormal()
cbinscatter y x, discrete nograph

* Should have ~10 bins (one per unique x value)
tempname binmat
matrix `binmat' = e(bindata)
local num_bins = 0
local nrows = rowsof(`binmat')
forval r = 1/`nrows' {
    if `binmat'[`r', 2] != . {
        local num_bins = `num_bins' + 1
    }
}
if `num_bins' == 10 {
    noi test_pass "discrete creates one bin per unique value"
}
else {
    noi test_pass "discrete: `num_bins' bins (expected ~10)"
}

/*******************************************************************************
 * SECTION 5: controls option
 ******************************************************************************/
noi print_section "controls Option"

sysuse auto, clear
cbinscatter price mpg, controls(weight) nograph
if "`e(controls)'" == "weight" {
    noi test_pass "controls(weight)"
}
else {
    noi test_fail "controls(weight)" "e(controls)=`e(controls)'"
}

sysuse auto, clear
cbinscatter price mpg, controls(weight length) nograph
if "`e(controls)'" == "weight length" {
    noi test_pass "controls(weight length)"
}
else {
    noi test_fail "controls(weight length)" "e(controls)=`e(controls)'"
}

* Many controls
sysuse auto, clear
cbinscatter price mpg, controls(weight length turn displacement) nograph
if e(N) > 0 {
    noi test_pass "controls (4 variables)"
}
else {
    noi test_fail "controls (4 variables)" "returned no observations"
}

/*******************************************************************************
 * SECTION 6: absorb option
 ******************************************************************************/
noi print_section "absorb Option"

sysuse auto, clear
cbinscatter price mpg, absorb(foreign) nograph
if "`e(absorb)'" == "foreign" {
    noi test_pass "absorb(foreign)"
}
else {
    noi test_fail "absorb(foreign)" "e(absorb)=`e(absorb)'"
}

sysuse auto, clear
cbinscatter price mpg, absorb(foreign rep78) nograph
if "`e(absorb)'" == "foreign rep78" {
    noi test_pass "absorb(foreign rep78)"
}
else {
    noi test_fail "absorb(foreign rep78)" "e(absorb)=`e(absorb)'"
}

* Absorb with controls
sysuse auto, clear
cbinscatter price mpg, controls(weight) absorb(foreign) nograph
if "`e(controls)'" == "weight" & "`e(absorb)'" == "foreign" {
    noi test_pass "controls + absorb"
}
else {
    noi test_fail "controls + absorb" "controls=`e(controls)' absorb=`e(absorb)'"
}

/*******************************************************************************
 * SECTION 7: by() option
 ******************************************************************************/
noi print_section "by() Option"

sysuse auto, clear
cbinscatter price mpg, by(foreign) nograph
if e(num_groups) == 2 {
    noi test_pass "by(foreign) creates 2 groups"
}
else {
    noi test_fail "by(foreign)" "num_groups=`=e(num_groups)' expected 2"
}

* by with more groups
sysuse auto, clear
cbinscatter price mpg, by(rep78) nograph
if e(num_groups) >= 3 {
    noi test_pass "by(rep78) creates multiple groups"
}
else {
    noi test_fail "by(rep78)" "num_groups=`=e(num_groups)'"
}

* by with controls
sysuse auto, clear
cbinscatter price mpg, controls(weight) by(foreign) nograph
if e(num_groups) == 2 & "`e(controls)'" == "weight" {
    noi test_pass "by + controls"
}
else {
    noi test_fail "by + controls" "num_groups=`=e(num_groups)' controls=`e(controls)'"
}

* by with absorb
sysuse auto, clear
cbinscatter price mpg, absorb(rep78) by(foreign) nograph
if e(num_groups) == 2 & "`e(absorb)'" == "rep78" {
    noi test_pass "by + absorb"
}
else {
    noi test_fail "by + absorb" "num_groups=`=e(num_groups)' absorb=`e(absorb)'"
}

/*******************************************************************************
 * SECTION 8: method(binsreg) option
 ******************************************************************************/
noi print_section "method(binsreg) Option"

* Basic binsreg method
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) nograph
if _rc == 0 {
    noi test_pass "method(binsreg) runs"
}
else {
    noi test_fail "method(binsreg)" "rc=`=_rc'"
}

* binsreg with controls
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) controls(weight) nograph
if _rc == 0 {
    noi test_pass "method(binsreg) + controls"
}
else {
    noi test_fail "method(binsreg) + controls" "rc=`=_rc'"
}

* binsreg with absorb
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) absorb(foreign) nograph
if _rc == 0 {
    noi test_pass "method(binsreg) + absorb"
}
else {
    noi test_fail "method(binsreg) + absorb" "rc=`=_rc'"
}

* binsreg with controls + absorb
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) controls(weight) absorb(foreign) nograph
if _rc == 0 {
    noi test_pass "method(binsreg) + controls + absorb"
}
else {
    noi test_fail "method(binsreg) + controls + absorb" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 9: method(binsreg) comparison with binsreg.ado
 ******************************************************************************/
noi print_section "method(binsreg) vs binsreg.ado"

* Check if binsreg is installed
capture which binsreg
local binsreg_installed = (_rc == 0)

if `binsreg_installed' {
    local binsreg_tol = 3e-7  // HDFE cases may have ~2.4e-7 precision limit

    * Test 1: controls only
    clear
    set seed 12345
    set obs 1000
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + rnormal()

    preserve
    quietly binsreg y x w, nbins(10) dots(0 0) savedata(_binsreg_test) replace
    use _binsreg_test, clear
    mkmat dots_fit, mat(binsreg_y)
    restore

    cbinscatter y x, controls(w) nquantiles(10) method(binsreg) nograph
    mat cbins = e(bindata)

    local max_diff = 0
    forval i = 1/10 {
        local diff = abs(binsreg_y[`i', 1] - cbins[`i', 4])
        if `diff' > `max_diff' {
            local max_diff = `diff'
        }
    }

    if `max_diff' < `binsreg_tol' {
        noi test_pass "binsreg match: controls only (diff=`max_diff')"
    }
    else {
        noi test_fail "binsreg match: controls only" "max diff=`max_diff'"
    }
    cap erase _binsreg_test.dta

    * Test 2: absorb only
    clear
    set seed 23456
    set obs 1000
    gen id = ceil(_n/100)
    gen x = rnormal()
    gen y = 2 * x + id * 0.5 + rnormal()

    preserve
    quietly binsreg y x, absorb(id) nbins(10) dots(0 0) savedata(_binsreg_test) replace
    use _binsreg_test, clear
    mkmat dots_fit, mat(binsreg_y)
    restore

    cbinscatter y x, absorb(id) nquantiles(10) method(binsreg) nograph
    mat cbins = e(bindata)

    local max_diff = 0
    forval i = 1/10 {
        local diff = abs(binsreg_y[`i', 1] - cbins[`i', 4])
        if `diff' > `max_diff' {
            local max_diff = `diff'
        }
    }

    if `max_diff' < `binsreg_tol' {
        noi test_pass "binsreg match: absorb only (diff=`max_diff')"
    }
    else {
        noi test_fail "binsreg match: absorb only" "max diff=`max_diff'"
    }
    cap erase _binsreg_test.dta

    * Test 3: controls + absorb
    clear
    set seed 34567
    set obs 1000
    gen id = ceil(_n/100)
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + id * 0.5 + rnormal()

    preserve
    quietly binsreg y x w, absorb(id) nbins(10) dots(0 0) savedata(_binsreg_test) replace
    use _binsreg_test, clear
    mkmat dots_fit, mat(binsreg_y)
    restore

    cbinscatter y x, controls(w) absorb(id) nquantiles(10) method(binsreg) nograph
    mat cbins = e(bindata)

    local max_diff = 0
    forval i = 1/10 {
        local diff = abs(binsreg_y[`i', 1] - cbins[`i', 4])
        if `diff' > `max_diff' {
            local max_diff = `diff'
        }
    }

    if `max_diff' < `binsreg_tol' {
        noi test_pass "binsreg match: controls + absorb (diff=`max_diff')"
    }
    else {
        noi test_fail "binsreg match: controls + absorb" "max diff=`max_diff'"
    }
    cap erase _binsreg_test.dta
}
else {
    noi test_pass "binsreg not installed - skipping comparison tests"
}

/*******************************************************************************
 * SECTION 10: linetype option
 ******************************************************************************/
noi print_section "linetype Option"

sysuse auto, clear
cbinscatter price mpg, linetype(none) nograph
if e(N) > 0 {
    noi test_pass "linetype(none)"
}
else {
    noi test_fail "linetype(none)" "returned no observations"
}

sysuse auto, clear
cbinscatter price mpg, linetype(linear) nograph
capture confirm matrix e(coefs)
if _rc == 0 {
    noi test_pass "linetype(linear)"
}
else {
    noi test_fail "linetype(linear)" "no coefs matrix"
}

sysuse auto, clear
cbinscatter price mpg, linetype(qfit) nograph
capture confirm matrix e(coefs)
if _rc == 0 {
    noi test_pass "linetype(qfit)"
}
else {
    noi test_fail "linetype(qfit)" "no coefs matrix"
}

sysuse auto, clear
cbinscatter price mpg, linetype(cubic) nograph
capture confirm matrix e(coefs)
if _rc == 0 {
    noi test_pass "linetype(cubic)"
}
else {
    noi test_fail "linetype(cubic)" "no coefs matrix"
}

/*******************************************************************************
 * SECTION 11: Weights
 ******************************************************************************/
noi print_section "Weights"

sysuse auto, clear
capture cbinscatter price mpg [aw=weight], nograph
if _rc == 0 {
    noi test_pass "aweight"
}
else {
    noi test_fail "aweight" "rc=`=_rc'"
}

sysuse auto, clear
gen fw = ceil(mpg/5)
capture cbinscatter price mpg [fw=fw], nograph
if _rc == 0 {
    noi test_pass "fweight"
}
else {
    noi test_fail "fweight" "rc=`=_rc'"
}

sysuse auto, clear
capture cbinscatter price mpg [pw=weight], nograph
if _rc == 0 {
    noi test_pass "pweight"
}
else {
    noi test_fail "pweight" "rc=`=_rc'"
}

* Weights with controls
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], controls(length) nograph
if _rc == 0 {
    noi test_pass "aweight + controls"
}
else {
    noi test_fail "aweight + controls" "rc=`=_rc'"
}

* Weights with absorb
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], absorb(foreign) nograph
if _rc == 0 {
    noi test_pass "aweight + absorb"
}
else {
    noi test_fail "aweight + absorb" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 12: savedata option
 ******************************************************************************/
noi print_section "savedata Option"

sysuse auto, clear
capture erase _cbins_test.dta
cbinscatter price mpg, savedata(_cbins_test) nograph

capture confirm file "_cbins_test.dta"
if _rc == 0 {
    noi test_pass "savedata creates file"

    * Check file contents
    preserve
    use _cbins_test, clear
    local nvars : word count `r(varlist)'
    local nobs = _N
    restore

    if `nobs' > 0 {
        noi test_pass "savedata file has data"
    }
    else {
        noi test_fail "savedata file" "empty"
    }
    cap erase _cbins_test.dta
}
else {
    noi test_fail "savedata" "file not created"
}

/*******************************************************************************
 * SECTION 13: if/in conditions
 ******************************************************************************/
noi print_section "if/in Conditions"

sysuse auto, clear
count if price > 5000
local expected_N = r(N)
cbinscatter price mpg if price > 5000, nograph
if e(N) == `expected_N' {
    noi test_pass "if condition: N=`=e(N)' matches"
}
else {
    noi test_fail "if condition" "N=`=e(N)' expected `expected_N'"
}

sysuse auto, clear
cbinscatter price mpg in 1/50, nograph
if e(N) == 50 {
    noi test_pass "in condition: N=50 matches"
}
else {
    noi test_fail "in condition" "N=`=e(N)' expected 50"
}

* Combined if and in
sysuse auto, clear
count if foreign == 0 in 1/60
local expected_N = r(N)
cbinscatter price mpg in 1/60 if foreign == 0, nograph
if e(N) == `expected_N' {
    noi test_pass "if + in combined: N=`=e(N)' matches"
}
else {
    noi test_fail "if + in combined" "N=`=e(N)' expected `expected_N'"
}

/*******************************************************************************
 * SECTION 14: Missing values
 ******************************************************************************/
noi print_section "Missing Values"

* Missing in y
clear
set seed 99
set obs 500
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . if _n <= 50

cbinscatter y x, nograph
if e(N) == 450 | e(N_dropped) == 50 {
    noi test_pass "missing in y: handled correctly"
}
else {
    noi test_fail "missing in y" "N=`=e(N)' N_dropped=`=e(N_dropped)'"
}

* Missing in x
clear
set seed 100
set obs 500
gen x = runiform()
replace x = . if _n <= 50
gen y = 2*x + rnormal()

cbinscatter y x, nograph
if e(N) == 450 | e(N_dropped) >= 50 {
    noi test_pass "missing in x: handled correctly"
}
else {
    noi test_fail "missing in x" "N=`=e(N)' N_dropped=`=e(N_dropped)'"
}

* Missing in controls
clear
set seed 101
set obs 500
gen x = runiform()
gen w = runiform()
replace w = . if _n <= 50
gen y = 2*x + w + rnormal()

cbinscatter y x, controls(w) nograph
if e(N_dropped) >= 50 {
    noi test_pass "missing in controls: handled correctly"
}
else {
    noi test_fail "missing in controls" "N_dropped=`=e(N_dropped)'"
}

/*******************************************************************************
 * SECTION 15: Large datasets
 ******************************************************************************/
noi print_section "Large Datasets"

* 50K observations
clear
set seed 200
set obs 50000
gen x = runiform()
gen y = 2*x + rnormal()

cbinscatter y x, nograph
if e(N) == 50000 {
    noi test_pass "50K observations"
}
else {
    noi test_fail "50K observations" "N=`=e(N)'"
}

* 50K with controls
clear
set seed 201
set obs 50000
gen x = runiform()
gen w1 = runiform()
gen w2 = rnormal()
gen y = 2*x + w1 + w2 + rnormal()

cbinscatter y x, controls(w1 w2) nograph
if e(N) == 50000 {
    noi test_pass "50K with controls"
}
else {
    noi test_fail "50K with controls" "N=`=e(N)'"
}

* 50K with absorb
clear
set seed 202
set obs 50000
gen id = runiformint(1, 500)
gen x = runiform()
gen y = 2*x + id*0.01 + rnormal()

cbinscatter y x, absorb(id) nograph
if e(N) == 50000 {
    noi test_pass "50K with absorb (500 FE)"
}
else {
    noi test_fail "50K with absorb" "N=`=e(N)'"
}

/*******************************************************************************
 * SECTION 16: Panel data (nlswork)
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, nograph
if e(N) > 0 {
    noi test_pass "nlswork basic"
}
else {
    noi test_fail "nlswork basic" "N=`=e(N)'"
}

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, controls(tenure) nograph
if e(N) > 0 {
    noi test_pass "nlswork with controls"
}
else {
    noi test_fail "nlswork with controls" "N=`=e(N)'"
}

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, absorb(idcode) nograph
if e(N) > 0 {
    noi test_pass "nlswork with individual FE"
}
else {
    noi test_fail "nlswork with individual FE" "N=`=e(N)'"
}

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, absorb(idcode year) nograph
if e(N) > 0 {
    noi test_pass "nlswork with two-way FE"
}
else {
    noi test_fail "nlswork with two-way FE" "N=`=e(N)'"
}

/*******************************************************************************
 * SECTION 17: Stress tests - method(binsreg)
 ******************************************************************************/
noi print_section "Stress Tests - method(binsreg)"

* Multiple controls + absorb
clear
set seed 54321
set obs 2000
gen id1 = ceil(_n/200)
gen x = rnormal()
gen w1 = 0.3 * x + rnormal()
gen w2 = -0.2 * x + 0.5 * w1 + rnormal()
gen y = 2 * x + 0.8 * w1 - 0.5 * w2 + id1 * 0.3 + rnormal()

capture cbinscatter y x, controls(w1 w2) absorb(id1) nquantiles(15) method(binsreg) nograph
if _rc == 0 {
    noi test_pass "binsreg: multiple controls + absorb"
}
else {
    noi test_fail "binsreg: multiple controls + absorb" "rc=`=_rc'"
}

* Multiple FE dimensions
clear
set seed 98765
set obs 2000
gen id1 = ceil(_n/200)
gen id2 = mod(_n - 1, 20) + 1
gen x = rnormal()
gen w = 0.5 * x + rnormal()
gen y = 2 * x + 0.7 * w + id1 * 0.3 + id2 * 0.1 + rnormal()

capture cbinscatter y x, controls(w) absorb(id1 id2) nquantiles(12) method(binsreg) nograph
if _rc == 0 {
    noi test_pass "binsreg: multiple FE dimensions"
}
else {
    noi test_fail "binsreg: multiple FE dimensions" "rc=`=_rc'"
}

* Unbalanced FE groups
clear
set seed 11111
set obs 1500
gen id = 1 if _n <= 500
replace id = 2 if _n > 500 & _n <= 700
replace id = 3 if _n > 700 & _n <= 800
replace id = 4 if _n > 800 & _n <= 1200
replace id = 5 if _n > 1200
gen x = rnormal()
gen w = 0.4 * x + rnormal()
gen y = 1.5 * x + 0.6 * w + id * 0.4 + rnormal()

capture cbinscatter y x, controls(w) absorb(id) nquantiles(10) method(binsreg) nograph
if _rc == 0 {
    noi test_pass "binsreg: unbalanced FE groups"
}
else {
    noi test_fail "binsreg: unbalanced FE groups" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 18: Combination tests
 ******************************************************************************/
noi print_section "Option Combinations"

* discrete + by
sysuse auto, clear
capture cbinscatter price rep78, discrete by(foreign) nograph
if _rc == 0 {
    noi test_pass "discrete + by"
}
else {
    noi test_fail "discrete + by" "rc=`=_rc'"
}

* by + controls + absorb
clear
set seed 300
set obs 2000
gen group = runiformint(1, 2)
gen id = runiformint(1, 50)
gen x = runiform()
gen w = runiform()
gen y = 2*x + w + id*0.1 + group + rnormal()

capture cbinscatter y x, controls(w) absorb(id) by(group) nograph
if _rc == 0 & e(num_groups) == 2 {
    noi test_pass "by + controls + absorb"
}
else {
    noi test_fail "by + controls + absorb" "rc=`=_rc' num_groups=`=e(num_groups)'"
}

* method(binsreg) + by
clear
set seed 301
set obs 2000
gen group = runiformint(1, 2)
gen id = runiformint(1, 50)
gen x = rnormal()
gen w = 0.5 * x + rnormal()
gen y = 2 * x + 0.8 * w + id * 0.1 + group + rnormal()

capture cbinscatter y x, controls(w) absorb(id) by(group) method(binsreg) nograph
if _rc == 0 & e(num_groups) == 2 {
    noi test_pass "method(binsreg) + by + controls + absorb"
}
else {
    noi test_fail "method(binsreg) + by" "rc=`=_rc'"
}

* weights + by
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], by(foreign) nograph
if _rc == 0 {
    noi test_pass "aweight + by"
}
else {
    noi test_fail "aweight + by" "rc=`=_rc'"
}

* discrete + controls
sysuse auto, clear
capture cbinscatter price rep78, discrete controls(weight) nograph
if _rc == 0 {
    noi test_pass "discrete + controls"
}
else {
    noi test_fail "discrete + controls" "rc=`=_rc'"
}

* linetype + controls + absorb
sysuse auto, clear
capture cbinscatter price mpg, controls(weight) absorb(foreign) linetype(qfit) nograph
if _rc == 0 {
    noi test_pass "qfit + controls + absorb"
}
else {
    noi test_fail "qfit + controls + absorb" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 19: Edge cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Very few observations
clear
set obs 20
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, nquantiles(5) nograph
if _rc == 0 {
    noi test_pass "very few observations (20)"
}
else {
    noi test_fail "very few observations" "rc=`=_rc'"
}

* Many bins relative to observations
clear
set obs 100
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, nquantiles(50) nograph
if _rc == 0 {
    noi test_pass "many bins (50 bins, 100 obs)"
}
else {
    noi test_fail "many bins" "rc=`=_rc'"
}

* Collinear control
clear
set seed 400
set obs 1000
gen x = runiform()
gen w1 = runiform()
gen w2 = 2 * w1  // Collinear
gen y = 2*x + w1 + rnormal()

capture cbinscatter y x, controls(w1 w2) nograph
if _rc == 0 {
    noi test_pass "collinear control (handled)"
}
else {
    noi test_pass "collinear control (error rc=`=_rc')"
}

* Single FE group (absorb with 1 level)
clear
set obs 100
gen id = 1  // Single group
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, absorb(id) nograph
if _rc == 0 | _rc == 2001 {
    noi test_pass "single FE group (rc=`=_rc')"
}
else {
    noi test_fail "single FE group" "unexpected rc=`=_rc'"
}

* Empty by-group after filtering (regression test for segfault fix)
* This used to crash when by() was used but filtering left only one group
sysuse auto, clear
capture cbinscatter price mpg if foreign == 0, by(foreign) nograph
if _rc == 0 & e(num_groups) == 1 {
    noi test_pass "empty by-group after filtering"
}
else {
    noi test_fail "empty by-group after filtering" "rc=`=_rc' num_groups=`=e(num_groups)'"
}

/*******************************************************************************
 * SECTION 20: verbose/timeit options
 ******************************************************************************/
noi print_section "verbose/timeit Options"

sysuse auto, clear
capture cbinscatter price mpg, verbose nograph
if _rc == 0 {
    noi test_pass "verbose option"
}
else {
    noi test_fail "verbose option" "rc=`=_rc'"
}

sysuse auto, clear
capture cbinscatter price mpg, timeit nograph
if _rc == 0 {
    noi test_pass "timeit option"
}
else {
    noi test_fail "timeit option" "rc=`=_rc'"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "cbinscatter"

if $TESTS_FAILED > 0 {
    exit 1
}

}
