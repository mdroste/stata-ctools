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

noi di as text "Running validation tests for cbinscatter..."

* Plugin check
sysuse auto, clear
capture cbinscatter price mpg, nograph
if _rc != 0 {
    test_fail "cbinscatter plugin load" "returned error `=_rc'"
    exit 1
}
test_pass "cbinscatter plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic functionality
 ******************************************************************************/
print_section "Basic Functionality"

* Simple two-variable scatter
sysuse auto, clear
cbinscatter price mpg, nograph
if e(N) == _N & e(nquantiles) == 20 {
    test_pass "basic scatter (default 20 bins)"
}
else {
    test_fail "basic scatter" "N=`=e(N)' expected `=_N', nquantiles=`=e(nquantiles)' expected 20"
}

* Check e() results populated
sysuse auto, clear
cbinscatter price mpg, nograph
capture confirm matrix e(bindata)
if _rc == 0 {
    test_pass "e(bindata) matrix populated"
}
else {
    test_fail "e(bindata)" "matrix not found"
}

capture confirm matrix e(coefs)
if _rc == 0 {
    test_pass "e(coefs) matrix populated (linear fit)"
}
else {
    test_fail "e(coefs)" "matrix not found"
}

/*******************************************************************************
 * SECTION 3: nquantiles option
 ******************************************************************************/
print_section "nquantiles Option"

sysuse auto, clear
cbinscatter price mpg, nquantiles(10) nograph
if e(nquantiles) == 10 {
    test_pass "nquantiles(10)"
}
else {
    test_fail "nquantiles(10)" "got `=e(nquantiles)'"
}

sysuse auto, clear
cbinscatter price mpg, nquantiles(50) nograph
if e(nquantiles) == 50 {
    test_pass "nquantiles(50)"
}
else {
    test_fail "nquantiles(50)" "got `=e(nquantiles)'"
}

sysuse auto, clear
cbinscatter price mpg, nquantiles(5) nograph
if e(nquantiles) == 5 {
    test_pass "nquantiles(5)"
}
else {
    test_fail "nquantiles(5)" "got `=e(nquantiles)'"
}

* Test invalid nquantiles (should error)
sysuse auto, clear
capture cbinscatter price mpg, nquantiles(1) nograph
if _rc != 0 {
    test_pass "nquantiles(1) rejected (rc=`=_rc')"
}
else {
    test_fail "nquantiles(1)" "should have been rejected"
}

/*******************************************************************************
 * SECTION 4: discrete option
 ******************************************************************************/
print_section "discrete Option"

* Create discrete x variable
sysuse auto, clear
cbinscatter price rep78, discrete nograph
if e(N) > 0 {
    test_pass "discrete with rep78 (5 values)"
}
else {
    test_fail "discrete" "returned no observations"
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
    test_pass "discrete creates one bin per unique value"
}
else {
    test_fail "discrete creates one bin per unique value" "got `num_bins' bins, expected 10"
}

/*******************************************************************************
 * SECTION 5: controls option
 ******************************************************************************/
print_section "controls Option"

sysuse auto, clear
cbinscatter price mpg, controls(weight) nograph
if "`e(controls)'" == "weight" {
    test_pass "controls(weight)"
}
else {
    test_fail "controls(weight)" "e(controls)=`e(controls)'"
}

sysuse auto, clear
cbinscatter price mpg, controls(weight length) nograph
if "`e(controls)'" == "weight length" {
    test_pass "controls(weight length)"
}
else {
    test_fail "controls(weight length)" "e(controls)=`e(controls)'"
}

* Many controls
sysuse auto, clear
cbinscatter price mpg, controls(weight length turn displacement) nograph
if e(N) > 0 {
    test_pass "controls (4 variables)"
}
else {
    test_fail "controls (4 variables)" "returned no observations"
}

/*******************************************************************************
 * SECTION 6: absorb option
 ******************************************************************************/
print_section "absorb Option"

sysuse auto, clear
cbinscatter price mpg, absorb(foreign) nograph
if "`e(absorb)'" == "foreign" {
    test_pass "absorb(foreign)"
}
else {
    test_fail "absorb(foreign)" "e(absorb)=`e(absorb)'"
}

sysuse auto, clear
cbinscatter price mpg, absorb(foreign rep78) nograph
if "`e(absorb)'" == "foreign rep78" {
    test_pass "absorb(foreign rep78)"
}
else {
    test_fail "absorb(foreign rep78)" "e(absorb)=`e(absorb)'"
}

* Absorb with controls
sysuse auto, clear
cbinscatter price mpg, controls(weight) absorb(foreign) nograph
if "`e(controls)'" == "weight" & "`e(absorb)'" == "foreign" {
    test_pass "controls + absorb"
}
else {
    test_fail "controls + absorb" "controls=`e(controls)' absorb=`e(absorb)'"
}

/*******************************************************************************
 * SECTION 7: by() option
 ******************************************************************************/
print_section "by() Option"

sysuse auto, clear
cbinscatter price mpg, by(foreign) nograph
if e(num_groups) == 2 {
    test_pass "by(foreign) creates 2 groups"
}
else {
    test_fail "by(foreign)" "num_groups=`=e(num_groups)' expected 2"
}

* by with more groups
sysuse auto, clear
cbinscatter price mpg, by(rep78) nograph
if e(num_groups) >= 3 {
    test_pass "by(rep78) creates multiple groups"
}
else {
    test_fail "by(rep78)" "num_groups=`=e(num_groups)'"
}

* by with controls
sysuse auto, clear
cbinscatter price mpg, controls(weight) by(foreign) nograph
if e(num_groups) == 2 & "`e(controls)'" == "weight" {
    test_pass "by + controls"
}
else {
    test_fail "by + controls" "num_groups=`=e(num_groups)' controls=`e(controls)'"
}

* by with absorb
sysuse auto, clear
cbinscatter price mpg, absorb(rep78) by(foreign) nograph
if e(num_groups) == 2 & "`e(absorb)'" == "rep78" {
    test_pass "by + absorb"
}
else {
    test_fail "by + absorb" "num_groups=`=e(num_groups)' absorb=`e(absorb)'"
}

/*******************************************************************************
 * SECTION 8: method(binsreg) option
 ******************************************************************************/
print_section "method(binsreg) Option"

* Basic binsreg method
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) == _N {
        test_pass "method(binsreg) runs (N=`=e(N)')"
    }
    else {
        test_fail "method(binsreg)" "bindata missing or N=`=e(N)' expected `=_N'"
    }
}
else {
    test_fail "method(binsreg)" "rc=`=_rc'"
}

* binsreg with controls
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) controls(weight) nograph
if _rc == 0 {
    if "`e(controls)'" == "weight" & e(N) == _N {
        test_pass "method(binsreg) + controls"
    }
    else {
        test_fail "method(binsreg) + controls" "controls=`e(controls)' N=`=e(N)'"
    }
}
else {
    test_fail "method(binsreg) + controls" "rc=`=_rc'"
}

* binsreg with absorb
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) absorb(foreign) nograph
if _rc == 0 {
    if "`e(absorb)'" == "foreign" & e(N) == _N {
        test_pass "method(binsreg) + absorb"
    }
    else {
        test_fail "method(binsreg) + absorb" "absorb=`e(absorb)' N=`=e(N)'"
    }
}
else {
    test_fail "method(binsreg) + absorb" "rc=`=_rc'"
}

* binsreg with controls + absorb
sysuse auto, clear
capture cbinscatter price mpg, method(binsreg) controls(weight) absorb(foreign) nograph
if _rc == 0 {
    if "`e(controls)'" == "weight" & "`e(absorb)'" == "foreign" & e(N) == _N {
        test_pass "method(binsreg) + controls + absorb"
    }
    else {
        test_fail "method(binsreg) + controls + absorb" "controls=`e(controls)' absorb=`e(absorb)' N=`=e(N)'"
    }
}
else {
    test_fail "method(binsreg) + controls + absorb" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 9: Comparison with binscatter.ado
 ******************************************************************************/
print_section "Comparison vs binscatter.ado"

* Check if binscatter is installed
capture which binscatter
local binscatter_installed = (_rc == 0)

if `binscatter_installed' {
    * For binscatter comparison, we compare the actual bin x and y means
    * binscatter saves as CSV, cbinscatter saves as DTA

    * Test 1: Basic scatter - compare bin means
    sysuse auto, clear
    quietly binscatter price mpg, nquantiles(10) savedata(_bs_test) replace

    import delimited _bs_test.csv, clear
    rename mpg bs_x
    rename price bs_y
    sort bs_x
    mkmat bs_x bs_y, mat(bs_data)

    sysuse auto, clear
    cbinscatter price mpg, nquantiles(10) nograph
    mat cb_data = e(bindata)

    * Compare x_mean (column 3 of cbinscatter) with binscatter x
    local min_sf_x = 15
    local min_sf_y = 15
    forval i = 1/10 {
        local bs_x = bs_data[`i', 1]
        local bs_y = bs_data[`i', 2]
        local cb_x = cb_data[`i', 3]
        local cb_y = cb_data[`i', 4]
        sigfigs `bs_x' `cb_x'
        if r(sigfigs) < `min_sf_x' {
            local min_sf_x = r(sigfigs)
        }
        sigfigs `bs_y' `cb_y'
        if r(sigfigs) < `min_sf_y' {
            local min_sf_y = r(sigfigs)
        }
    }

    * Require at least $DEFAULT_SIGFIGS significant figures for bin means
    if `min_sf_x' >= $DEFAULT_SIGFIGS & `min_sf_y' >= $DEFAULT_SIGFIGS {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_pass "binscatter match: basic (x:`x_fmt', y:`y_fmt' sigfigs)"
    }
    else {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_fail "binscatter match: basic" "x:`x_fmt', y:`y_fmt' sigfigs, need $DEFAULT_SIGFIGS"
    }
    cap erase _bs_test.csv
    cap erase _bs_test.do

    * Test 2: With controls - compare bin means
    * Note: binscatter adds sample means back to residuals; cbinscatter doesn't
    * So we compare after adjusting for mean differences
    sysuse auto, clear
    sum mpg, meanonly
    local mean_x = r(mean)
    sum price, meanonly
    local mean_y = r(mean)

    quietly binscatter price mpg, controls(weight) nquantiles(10) savedata(_bs_test) replace

    import delimited _bs_test.csv, clear
    rename mpg bs_x
    rename price bs_y
    sort bs_x
    mkmat bs_x bs_y, mat(bs_data)

    sysuse auto, clear
    cbinscatter price mpg, controls(weight) nquantiles(10) nograph
    mat cb_data = e(bindata)

    * Get mean of cbinscatter output to compute offset
    * Use high-precision format to avoid truncation in offset calculation
    mata: st_local("cb_mean_x", strofreal(mean(st_matrix("cb_data")[., 3]), "%21.15g"))
    mata: st_local("cb_mean_y", strofreal(mean(st_matrix("cb_data")[., 4]), "%21.15g"))
    mata: st_local("bs_mean_x", strofreal(mean(st_matrix("bs_data")[., 1]), "%21.15g"))
    mata: st_local("bs_mean_y", strofreal(mean(st_matrix("bs_data")[., 2]), "%21.15g"))

    local offset_x = `bs_mean_x' - `cb_mean_x'
    local offset_y = `bs_mean_y' - `cb_mean_y'

    local min_sf_x = 15
    local min_sf_y = 15
    forval i = 1/10 {
        local bs_x = bs_data[`i', 1]
        local bs_y = bs_data[`i', 2]
        local cb_x = cb_data[`i', 3] + `offset_x'
        local cb_y = cb_data[`i', 4] + `offset_y'
        sigfigs `bs_x' `cb_x'
        if r(sigfigs) < `min_sf_x' {
            local min_sf_x = r(sigfigs)
        }
        sigfigs `bs_y' `cb_y'
        if r(sigfigs) < `min_sf_y' {
            local min_sf_y = r(sigfigs)
        }
    }

    if `min_sf_x' >= $DEFAULT_SIGFIGS & `min_sf_y' >= $DEFAULT_SIGFIGS {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_pass "binscatter match: controls (x:`x_fmt', y:`y_fmt' sigfigs)"
    }
    else {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_fail "binscatter match: controls" "x:`x_fmt', y:`y_fmt' sigfigs, need $DEFAULT_SIGFIGS"
    }
    cap erase _bs_test.csv
    cap erase _bs_test.do

    * Test 3: With absorb - compare bin means (with mean adjustment)
    * For absorb, use overall raw data means for offset since HDFE produces
    * zero-mean residuals, and binscatter adds means back
    sysuse auto, clear
    sum mpg, meanonly
    local raw_mean_x = r(mean)
    sum price, meanonly
    local raw_mean_y = r(mean)

    quietly binscatter price mpg, absorb(foreign) nquantiles(10) savedata(_bs_test) replace

    import delimited _bs_test.csv, clear
    rename mpg bs_x
    rename price bs_y
    sort bs_x
    mkmat bs_x bs_y, mat(bs_data)

    sysuse auto, clear
    cbinscatter price mpg, absorb(foreign) nquantiles(10) nograph
    mat cb_data = e(bindata)

    * Use raw data means for offset (HDFE residuals are near-zero mean)
    local offset_x = `raw_mean_x'
    local offset_y = `raw_mean_y'

    local min_sf_x = 15
    local min_sf_y = 15
    forval i = 1/10 {
        local bs_x = bs_data[`i', 1]
        local bs_y = bs_data[`i', 2]
        local cb_x = cb_data[`i', 3] + `offset_x'
        local cb_y = cb_data[`i', 4] + `offset_y'
        sigfigs `bs_x' `cb_x'
        if r(sigfigs) < `min_sf_x' {
            local min_sf_x = r(sigfigs)
        }
        sigfigs `bs_y' `cb_y'
        if r(sigfigs) < `min_sf_y' {
            local min_sf_y = r(sigfigs)
        }
    }

    if `min_sf_x' >= $DEFAULT_SIGFIGS & `min_sf_y' >= $DEFAULT_SIGFIGS {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_pass "binscatter match: absorb (x:`x_fmt', y:`y_fmt' sigfigs)"
    }
    else {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_fail "binscatter match: absorb" "x:`x_fmt', y:`y_fmt' sigfigs, need $DEFAULT_SIGFIGS"
    }
    cap erase _bs_test.csv
    cap erase _bs_test.do

    * Test 4: Controls + absorb - compare bin means (with mean adjustment)
    * Use raw data means for offset
    sysuse auto, clear
    sum mpg, meanonly
    local raw_mean_x = r(mean)
    sum price, meanonly
    local raw_mean_y = r(mean)

    quietly binscatter price mpg, controls(weight) absorb(foreign) nquantiles(10) savedata(_bs_test) replace

    import delimited _bs_test.csv, clear
    rename mpg bs_x
    rename price bs_y
    sort bs_x
    mkmat bs_x bs_y, mat(bs_data)

    sysuse auto, clear
    cbinscatter price mpg, controls(weight) absorb(foreign) nquantiles(10) nograph
    mat cb_data = e(bindata)

    * Use raw data means for offset
    local offset_x = `raw_mean_x'
    local offset_y = `raw_mean_y'

    local min_sf_x = 15
    local min_sf_y = 15
    forval i = 1/10 {
        local bs_x = bs_data[`i', 1]
        local bs_y = bs_data[`i', 2]
        local cb_x = cb_data[`i', 3] + `offset_x'
        local cb_y = cb_data[`i', 4] + `offset_y'
        sigfigs `bs_x' `cb_x'
        if r(sigfigs) < `min_sf_x' {
            local min_sf_x = r(sigfigs)
        }
        sigfigs `bs_y' `cb_y'
        if r(sigfigs) < `min_sf_y' {
            local min_sf_y = r(sigfigs)
        }
    }

    if `min_sf_x' >= $DEFAULT_SIGFIGS & `min_sf_y' >= $DEFAULT_SIGFIGS {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_pass "binscatter match: controls+absorb (x:`x_fmt', y:`y_fmt' sigfigs)"
    }
    else {
        local x_fmt : display %4.1f `min_sf_x'
        local y_fmt : display %4.1f `min_sf_y'
        test_fail "binscatter match: controls+absorb" "x:`x_fmt', y:`y_fmt' sigfigs, need $DEFAULT_SIGFIGS"
    }
    cap erase _bs_test.csv
    cap erase _bs_test.do
}
else {
    test_fail "binscatter comparison" "binscatter not installed - cannot validate"
}

/*******************************************************************************
 * SECTION 10: method(binsreg) comparison with binsreg.ado
 ******************************************************************************/
print_section "method(binsreg) vs binsreg.ado"

* Check if binsreg is installed
capture which binsreg
local binsreg_installed = (_rc == 0)

if `binsreg_installed' {
    * Use significant figures comparison for binsreg tests
    * This is more appropriate when comparing two independent implementations
    * as numerical precision can differ slightly due to algorithm differences
    local min_sigfigs = $DEFAULT_SIGFIGS

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
    mkmat dots_x, mat(binsreg_x)
    restore

    cbinscatter y x, controls(w) nquantiles(10) method(binsreg) nograph
    mat cbins = e(bindata)

    * Compare Y values
    local min_sf_y = 15
    forval i = 1/10 {
        local val1 = binsreg_y[`i', 1]
        local val2 = cbins[`i', 4]
        sigfigs `val1' `val2'
        local sf = r(sigfigs)
        if `sf' < `min_sf_y' {
            local min_sf_y = `sf'
        }
    }

    if `min_sf_y' >= `min_sigfigs' {
        local sf_fmt : display %4.1f `min_sf_y'
        test_pass "binsreg match: controls only Y (sigfigs=`sf_fmt')"
    }
    else {
        local sf_fmt : display %4.1f `min_sf_y'
        test_fail "binsreg match: controls only Y" "sigfigs=`sf_fmt', need `min_sigfigs'"
    }

    * Also compare X values (bin centers)
    local min_sf_x = 15
    forval i = 1/10 {
        local val1 = binsreg_x[`i', 1]
        local val2 = cbins[`i', 3]
        sigfigs `val1' `val2'
        local sf = r(sigfigs)
        if `sf' < `min_sf_x' {
            local min_sf_x = `sf'
        }
    }

    if `min_sf_x' >= `min_sigfigs' {
        local sf_fmt : display %4.1f `min_sf_x'
        test_pass "binsreg match: controls only X (sigfigs=`sf_fmt')"
    }
    else {
        local sf_fmt : display %4.1f `min_sf_x'
        test_fail "binsreg match: controls only X" "sigfigs=`sf_fmt', need `min_sigfigs'"
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
    mkmat dots_x, mat(binsreg_x)
    restore

    cbinscatter y x, absorb(id) nquantiles(10) method(binsreg) nograph
    mat cbins = e(bindata)

    * Compare Y values
    local min_sf_y = 15
    forval i = 1/10 {
        local val1 = binsreg_y[`i', 1]
        local val2 = cbins[`i', 4]
        sigfigs `val1' `val2'
        local sf = r(sigfigs)
        if `sf' < `min_sf_y' {
            local min_sf_y = `sf'
        }
    }

    if `min_sf_y' >= `min_sigfigs' {
        local sf_fmt : display %4.1f `min_sf_y'
        test_pass "binsreg match: absorb only Y (sigfigs=`sf_fmt')"
    }
    else {
        local sf_fmt : display %4.1f `min_sf_y'
        test_fail "binsreg match: absorb only Y" "sigfigs=`sf_fmt', need `min_sigfigs'"
    }

    * Also compare X values (bin centers)
    local min_sf_x = 15
    forval i = 1/10 {
        local val1 = binsreg_x[`i', 1]
        local val2 = cbins[`i', 3]
        sigfigs `val1' `val2'
        local sf = r(sigfigs)
        if `sf' < `min_sf_x' {
            local min_sf_x = `sf'
        }
    }

    if `min_sf_x' >= `min_sigfigs' {
        local sf_fmt : display %4.1f `min_sf_x'
        test_pass "binsreg match: absorb only X (sigfigs=`sf_fmt')"
    }
    else {
        local sf_fmt : display %4.1f `min_sf_x'
        test_fail "binsreg match: absorb only X" "sigfigs=`sf_fmt', need `min_sigfigs'"
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
    mkmat dots_x, mat(binsreg_x)
    restore

    cbinscatter y x, controls(w) absorb(id) nquantiles(10) method(binsreg) nograph
    mat cbins = e(bindata)

    * Compare Y values
    local min_sf_y = 15
    forval i = 1/10 {
        local val1 = binsreg_y[`i', 1]
        local val2 = cbins[`i', 4]
        sigfigs `val1' `val2'
        local sf = r(sigfigs)
        if `sf' < `min_sf_y' {
            local min_sf_y = `sf'
        }
    }

    if `min_sf_y' >= `min_sigfigs' {
        local sf_fmt : display %4.1f `min_sf_y'
        test_pass "binsreg match: controls + absorb Y (sigfigs=`sf_fmt')"
    }
    else {
        local sf_fmt : display %4.1f `min_sf_y'
        test_fail "binsreg match: controls + absorb Y" "sigfigs=`sf_fmt', need `min_sigfigs'"
    }

    * Also compare X values (bin centers)
    local min_sf_x = 15
    forval i = 1/10 {
        local val1 = binsreg_x[`i', 1]
        local val2 = cbins[`i', 3]
        sigfigs `val1' `val2'
        local sf = r(sigfigs)
        if `sf' < `min_sf_x' {
            local min_sf_x = `sf'
        }
    }

    if `min_sf_x' >= `min_sigfigs' {
        local sf_fmt : display %4.1f `min_sf_x'
        test_pass "binsreg match: controls + absorb X (sigfigs=`sf_fmt')"
    }
    else {
        local sf_fmt : display %4.1f `min_sf_x'
        test_fail "binsreg match: controls + absorb X" "sigfigs=`sf_fmt', need `min_sigfigs'"
    }
    cap erase _binsreg_test.dta
}
else {
    test_fail "binsreg comparison" "binsreg not installed - cannot validate"
}

/*******************************************************************************
 * SECTION 10: linetype option
 ******************************************************************************/
print_section "linetype Option"

sysuse auto, clear
cbinscatter price mpg, linetype(none) nograph
if e(N) > 0 {
    test_pass "linetype(none)"
}
else {
    test_fail "linetype(none)" "returned no observations"
}

sysuse auto, clear
cbinscatter price mpg, linetype(linear) nograph
capture confirm matrix e(coefs)
if _rc == 0 {
    test_pass "linetype(linear)"
}
else {
    test_fail "linetype(linear)" "no coefs matrix"
}

sysuse auto, clear
cbinscatter price mpg, linetype(qfit) nograph
capture confirm matrix e(coefs)
if _rc == 0 {
    test_pass "linetype(qfit)"
}
else {
    test_fail "linetype(qfit)" "no coefs matrix"
}

sysuse auto, clear
cbinscatter price mpg, linetype(cubic) nograph
capture confirm matrix e(coefs)
if _rc == 0 {
    test_pass "linetype(cubic)"
}
else {
    test_fail "linetype(cubic)" "no coefs matrix"
}

/*******************************************************************************
 * SECTION 11: Weights
 ******************************************************************************/
print_section "Weights"

sysuse auto, clear
capture cbinscatter price mpg [aw=weight], nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) == _N {
        test_pass "aweight (N=`=e(N)')"
    }
    else {
        test_fail "aweight" "bindata missing or N=`=e(N)' expected `=_N'"
    }
}
else {
    test_fail "aweight" "rc=`=_rc'"
}

sysuse auto, clear
gen fw = ceil(mpg/5)
capture cbinscatter price mpg [fw=fw], nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) > 0 {
        test_pass "fweight (N=`=e(N)')"
    }
    else {
        test_fail "fweight" "bindata missing or N=`=e(N)'"
    }
}
else {
    test_fail "fweight" "rc=`=_rc'"
}

sysuse auto, clear
capture cbinscatter price mpg [pw=weight], nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) == _N {
        test_pass "pweight (N=`=e(N)')"
    }
    else {
        test_fail "pweight" "bindata missing or N=`=e(N)' expected `=_N'"
    }
}
else {
    test_fail "pweight" "rc=`=_rc'"
}

* Weights with controls
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], controls(length) nograph
if _rc == 0 {
    if "`e(controls)'" == "length" & e(N) == _N {
        test_pass "aweight + controls"
    }
    else {
        test_fail "aweight + controls" "controls=`e(controls)' N=`=e(N)'"
    }
}
else {
    test_fail "aweight + controls" "rc=`=_rc'"
}

* Weights with absorb
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], absorb(foreign) nograph
if _rc == 0 {
    if "`e(absorb)'" == "foreign" & e(N) == _N {
        test_pass "aweight + absorb"
    }
    else {
        test_fail "aweight + absorb" "absorb=`e(absorb)' N=`=e(N)'"
    }
}
else {
    test_fail "aweight + absorb" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 12: savedata option
 ******************************************************************************/
print_section "savedata Option"

sysuse auto, clear
capture erase _cbins_test.dta
cbinscatter price mpg, savedata(_cbins_test) nograph

capture confirm file "_cbins_test.dta"
if _rc == 0 {
    test_pass "savedata creates file"

    * Check file contents
    preserve
    use _cbins_test, clear
    local nvars : word count `r(varlist)'
    local nobs = _N
    restore

    if `nobs' > 0 {
        test_pass "savedata file has data"
    }
    else {
        test_fail "savedata file" "empty"
    }
    cap erase _cbins_test.dta
}
else {
    test_fail "savedata" "file not created"
}

/*******************************************************************************
 * SECTION 13: if/in conditions
 ******************************************************************************/
print_section "if/in Conditions"

sysuse auto, clear
count if price > 5000
local expected_N = r(N)
cbinscatter price mpg if price > 5000, nograph
if e(N) == `expected_N' {
    test_pass "if condition: N=`=e(N)' matches"
}
else {
    test_fail "if condition" "N=`=e(N)' expected `expected_N'"
}

sysuse auto, clear
cbinscatter price mpg in 1/50, nograph
if e(N) == 50 {
    test_pass "in condition: N=50 matches"
}
else {
    test_fail "in condition" "N=`=e(N)' expected 50"
}

* Combined if and in
sysuse auto, clear
count if foreign == 0 in 1/60
local expected_N = r(N)
cbinscatter price mpg in 1/60 if foreign == 0, nograph
if e(N) == `expected_N' {
    test_pass "if + in combined: N=`=e(N)' matches"
}
else {
    test_fail "if + in combined" "N=`=e(N)' expected `expected_N'"
}

/*******************************************************************************
 * SECTION 14: Missing values
 ******************************************************************************/
print_section "Missing Values"

* Missing in y
clear
set seed 99
set obs 500
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . if _n <= 50

cbinscatter y x, nograph
if e(N) == 450 & e(N_dropped) == 50 {
    test_pass "missing in y: handled correctly"
}
else {
    test_fail "missing in y" "N=`=e(N)' N_dropped=`=e(N_dropped)'"
}

* Missing in x
clear
set seed 100
set obs 500
gen x = runiform()
replace x = . if _n <= 50
gen y = 2*x + rnormal()

cbinscatter y x, nograph
if e(N) == 450 & e(N_dropped) >= 50 {
    test_pass "missing in x: handled correctly"
}
else {
    test_fail "missing in x" "N=`=e(N)' N_dropped=`=e(N_dropped)'"
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
    test_pass "missing in controls: handled correctly"
}
else {
    test_fail "missing in controls" "N_dropped=`=e(N_dropped)'"
}

/*******************************************************************************
 * SECTION 15: Large datasets
 ******************************************************************************/
print_section "Large Datasets"

* 50K observations
clear
set seed 200
set obs 50000
gen x = runiform()
gen y = 2*x + rnormal()

cbinscatter y x, nograph
if e(N) == 50000 {
    test_pass "50K observations"
}
else {
    test_fail "50K observations" "N=`=e(N)'"
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
    test_pass "50K with controls"
}
else {
    test_fail "50K with controls" "N=`=e(N)'"
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
    test_pass "50K with absorb (500 FE)"
}
else {
    test_fail "50K with absorb" "N=`=e(N)'"
}

/*******************************************************************************
 * SECTION 16: Panel data (nlswork)
 ******************************************************************************/
print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, nograph
if e(N) > 0 {
    test_pass "nlswork basic"
}
else {
    test_fail "nlswork basic" "N=`=e(N)'"
}

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, controls(tenure) nograph
if e(N) > 0 {
    test_pass "nlswork with controls"
}
else {
    test_fail "nlswork with controls" "N=`=e(N)'"
}

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, absorb(idcode) nograph
if e(N) > 0 {
    test_pass "nlswork with individual FE"
}
else {
    test_fail "nlswork with individual FE" "N=`=e(N)'"
}

webuse nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, absorb(idcode year) nograph
if e(N) > 0 {
    test_pass "nlswork with two-way FE"
}
else {
    test_fail "nlswork with two-way FE" "N=`=e(N)'"
}

/*******************************************************************************
 * SECTION 17: Stress tests - method(binsreg)
 ******************************************************************************/
print_section "Stress Tests - method(binsreg)"

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
    mat _bd = e(bindata)
    local nrows = rowsof(_bd)
    if `nrows' == 15 & e(N) == 2000 & e(nquantiles) == 15 {
        test_pass "binsreg: multiple controls + absorb (bins=`nrows', N=`=e(N)')"
    }
    else {
        test_fail "binsreg: multiple controls + absorb" "bins=`nrows' expected 15, N=`=e(N)' expected 2000"
    }
}
else {
    test_fail "binsreg: multiple controls + absorb" "rc=`=_rc'"
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
    mat _bd = e(bindata)
    local nrows = rowsof(_bd)
    if `nrows' == 12 & e(N) == 2000 & e(nquantiles) == 12 {
        test_pass "binsreg: multiple FE dimensions (bins=`nrows', N=`=e(N)')"
    }
    else {
        test_fail "binsreg: multiple FE dimensions" "bins=`nrows' expected 12, N=`=e(N)' expected 2000"
    }
}
else {
    test_fail "binsreg: multiple FE dimensions" "rc=`=_rc'"
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
    mat _bd = e(bindata)
    local nrows = rowsof(_bd)
    if `nrows' == 10 & e(N) == 1500 & e(nquantiles) == 10 {
        test_pass "binsreg: unbalanced FE groups (bins=`nrows', N=`=e(N)')"
    }
    else {
        test_fail "binsreg: unbalanced FE groups" "bins=`nrows' expected 10, N=`=e(N)' expected 1500"
    }
}
else {
    test_fail "binsreg: unbalanced FE groups" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 18: Combination tests
 ******************************************************************************/
print_section "Option Combinations"

* discrete + by
sysuse auto, clear
capture cbinscatter price rep78, discrete by(foreign) nograph
if _rc == 0 {
    * Verify output has correct number of groups
    if e(num_groups) == 2 {
        * Check bindata has rows: discrete rep78 has 5 unique values, by(foreign) has 2 groups
        * bindata should have non-missing rows for each group
        mat _bd = e(bindata)
        local nrows = rowsof(_bd)
        if `nrows' >= 5 {
            test_pass "discrete + by (groups=`=e(num_groups)', rows=`nrows')"
        }
        else {
            test_fail "discrete + by" "bindata has only `nrows' rows, expected >= 5"
        }
    }
    else {
        test_fail "discrete + by" "num_groups=`=e(num_groups)' expected 2"
    }
}
else {
    test_fail "discrete + by" "rc=`=_rc'"
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
    test_pass "by + controls + absorb"
}
else {
    test_fail "by + controls + absorb" "rc=`=_rc' num_groups=`=e(num_groups)'"
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
    test_pass "method(binsreg) + by + controls + absorb"
}
else {
    test_fail "method(binsreg) + by" "rc=`=_rc'"
}

* weights + by
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], by(foreign) nograph
if _rc == 0 {
    if e(num_groups) == 2 & e(N) == _N {
        test_pass "aweight + by (groups=`=e(num_groups)', N=`=e(N)')"
    }
    else {
        test_fail "aweight + by" "num_groups=`=e(num_groups)' expected 2, N=`=e(N)' expected `=_N'"
    }
}
else {
    test_fail "aweight + by" "rc=`=_rc'"
}

* discrete + controls
sysuse auto, clear
capture cbinscatter price rep78, discrete controls(weight) nograph
if _rc == 0 {
    * Verify e(controls) is set and bindata has correct discrete bins
    if "`e(controls)'" == "weight" {
        mat _bd = e(bindata)
        local nrows = rowsof(_bd)
        if `nrows' >= 5 {
            test_pass "discrete + controls (controls=`e(controls)', rows=`nrows')"
        }
        else {
            test_fail "discrete + controls" "bindata has only `nrows' rows, expected >= 5"
        }
    }
    else {
        test_fail "discrete + controls" "e(controls)=`e(controls)' expected weight"
    }
}
else {
    test_fail "discrete + controls" "rc=`=_rc'"
}

* linetype + controls + absorb
sysuse auto, clear
capture cbinscatter price mpg, controls(weight) absorb(foreign) linetype(qfit) nograph
if _rc == 0 {
    * Verify coefs matrix exists (qfit should produce polynomial coefficients)
    capture confirm matrix e(coefs)
    if _rc == 0 & "`e(controls)'" == "weight" & "`e(absorb)'" == "foreign" {
        test_pass "qfit + controls + absorb"
    }
    else {
        test_fail "qfit + controls + absorb" "coefs rc=`=_rc' controls=`e(controls)' absorb=`e(absorb)'"
    }
}
else {
    test_fail "qfit + controls + absorb" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 19: Edge cases
 ******************************************************************************/
print_section "Edge Cases"

* Very few observations
clear
set obs 20
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, nquantiles(5) nograph
if _rc == 0 {
    if e(N) == 20 & e(nquantiles) == 5 {
        mat _bd = e(bindata)
        local nrows = rowsof(_bd)
        if `nrows' == 5 {
            test_pass "very few observations (20 obs, 5 bins)"
        }
        else {
            test_fail "very few observations" "bindata rows=`nrows' expected 5"
        }
    }
    else {
        test_fail "very few observations" "N=`=e(N)' expected 20, nquantiles=`=e(nquantiles)' expected 5"
    }
}
else {
    test_fail "very few observations" "rc=`=_rc'"
}

* Many bins relative to observations
clear
set obs 100
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, nquantiles(50) nograph
if _rc == 0 {
    if e(N) == 100 & e(nquantiles) == 50 {
        mat _bd = e(bindata)
        local nrows = rowsof(_bd)
        if `nrows' == 50 {
            test_pass "many bins (50 bins, 100 obs)"
        }
        else {
            test_fail "many bins" "bindata rows=`nrows' expected 50"
        }
    }
    else {
        test_fail "many bins" "N=`=e(N)' expected 100, nquantiles=`=e(nquantiles)' expected 50"
    }
}
else {
    test_fail "many bins" "rc=`=_rc'"
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
    test_pass "collinear control (handled)"
}
else {
    test_fail "collinear control" "unexpected error rc=`=_rc'"
}

* Single FE group (absorb with 1 level)
clear
set obs 100
gen id = 1  // Single group
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, absorb(id) nograph
if _rc == 0 | _rc == 2001 {
    test_pass "single FE group (rc=`=_rc')"
}
else {
    test_fail "single FE group" "unexpected rc=`=_rc'"
}

* Empty by-group after filtering (regression test for segfault fix)
* This used to crash when by() was used but filtering left only one group
sysuse auto, clear
capture cbinscatter price mpg if foreign == 0, by(foreign) nograph
if _rc == 0 & e(num_groups) == 1 {
    test_pass "empty by-group after filtering"
}
else {
    test_fail "empty by-group after filtering" "rc=`=_rc' num_groups=`=e(num_groups)'"
}

/*******************************************************************************
 * SECTION 20: verbose/timeit options
 ******************************************************************************/
print_section "verbose/timeit Options"

sysuse auto, clear
capture cbinscatter price mpg, verbose nograph
if _rc == 0 {
    test_pass "verbose option"
}
else {
    test_fail "verbose option" "rc=`=_rc'"
}

sysuse auto, clear
capture cbinscatter price mpg, timeit nograph
if _rc == 0 {
    test_pass "timeit option"
}
else {
    test_fail "timeit option" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cbinscatter returns the same error codes as binscatter
 * when given invalid inputs or error conditions.
 * Note: binscatter is a user-written command (ssc install binscatter).
 * If binscatter is not installed, tests compare against expected behavior.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Check if binscatter is installed
capture which binscatter
local binscatter_installed = (_rc == 0)

* Variable doesn't exist
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price nonexistent_var, nograph) ctools_cmd(cbinscatter price nonexistent_var, nograph) testname("nonexistent variable")
}
else {
    capture cbinscatter price nonexistent_var, nograph
    if _rc != 0 {
        test_pass "[error] nonexistent variable (rc=`=_rc') [binscatter not installed]"
    }
    else {
        test_fail "[error] nonexistent variable" "should have errored"
    }
}

* String variable as y
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter make mpg, nograph) ctools_cmd(cbinscatter make mpg, nograph) testname("string y variable")
}
else {
    capture cbinscatter make mpg, nograph
    if _rc != 0 {
        test_pass "[error] string y variable (rc=`=_rc') [binscatter not installed]"
    }
    else {
        test_fail "[error] string y variable" "should have errored"
    }
}

* String variable as x
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price make, nograph) ctools_cmd(cbinscatter price make, nograph) testname("string x variable")
}
else {
    capture cbinscatter price make, nograph
    if _rc != 0 {
        test_pass "[error] string x variable (rc=`=_rc') [binscatter not installed]"
    }
    else {
        test_fail "[error] string x variable" "should have errored"
    }
}

* Invalid nquantiles (0)
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price mpg, nquantiles(0) nograph) ctools_cmd(cbinscatter price mpg, nquantiles(0) nograph) testname("nquantiles(0)")
}
else {
    capture cbinscatter price mpg, nquantiles(0) nograph
    if _rc != 0 {
        test_pass "[error] nquantiles(0) (rc=`=_rc') [binscatter not installed]"
    }
    else {
        test_fail "[error] nquantiles(0)" "should have errored"
    }
}

* Invalid nquantiles (negative)
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price mpg, nquantiles(-5) nograph) ctools_cmd(cbinscatter price mpg, nquantiles(-5) nograph) testname("negative nquantiles")
}
else {
    capture cbinscatter price mpg, nquantiles(-5) nograph
    if _rc != 0 {
        test_pass "[error] negative nquantiles (rc=`=_rc') [binscatter not installed]"
    }
    else {
        test_fail "[error] negative nquantiles" "should have errored"
    }
}

* No observations after if
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price mpg if price > 100000, nograph) ctools_cmd(cbinscatter price mpg if price > 100000, nograph) testname("no observations after if")
}
else {
    capture cbinscatter price mpg if price > 100000, nograph
    if _rc != 0 {
        test_pass "[error] no observations after if (rc=`=_rc') [binscatter not installed]"
    }
    else {
        test_fail "[error] no observations after if" "should have errored"
    }
}

* End of cbinscatter validation
noi print_summary "cbinscatter"
}
