* cbinscatter: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cbinscatter"

if "`__case'" == "cbinscatter" {
* validate_cbinscatter.do
* Comprehensive validation tests for cbinscatter
* Tests all options: nquantiles, discrete, controls, absorb, by, method(binsreg),
* linetype, weights, savedata
* VERIFICATION: Tests compare e() results and check consistency with binscatter
* and binsreg where applicable

* Load setup (works from project root or validation dir)
do "validation/support/validate_setup.do"

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

* Check if binscatter is installed (used for value comparison throughout)
capture which binscatter
local binscatter_installed = (_rc == 0)

if `binscatter_installed' {
    * -------------------------------------------------------------------------
    * Helper: compare_binscatter_bins
    *   Runs binscatter and cbinscatter, compares bin X and Y values.
    *   Arguments:
    *     bs_cmd   - full binscatter command (must include savedata(_bs_test) replace)
    *     cb_cmd   - full cbinscatter command (must include nograph)
    *     nbins    - number of bins requested
    *     yvar     - y variable name (for renaming CSV columns)
    *     xvar     - x variable name (for renaming CSV columns)
    *     testname - test label
    *     minsf    - minimum significant figures (default $DEFAULT_SIGFIGS)
    * -------------------------------------------------------------------------
    capture program drop compare_binscatter_bins
    program define compare_binscatter_bins
        syntax, bs_cmd(string asis) cb_cmd(string asis) nbins(integer) yvar(string) xvar(string) testname(string) [minsf(real $DEFAULT_SIGFIGS)]

        local min_sigfigs = `minsf'

        * Run binscatter and extract bin data from CSV
        preserve
        quietly `bs_cmd'
        import delimited _bs_test.csv, clear asdouble
        rename `xvar' bs_x
        rename `yvar' bs_y
        sort bs_x
        mkmat bs_x bs_y, mat(_bs_data)
        local bs_nbins = _N
        restore

        * Run cbinscatter
        `cb_cmd'
        mat _cb_data = e(bindata)

        * Count non-missing cbinscatter bins
        local cb_nbins = 0
        local nrows = rowsof(_cb_data)
        forval r = 1/`nrows' {
            if _cb_data[`r', 2] != . {
                local cb_nbins = `cb_nbins' + 1
            }
        }

        * Check bin counts match before comparing values
        if `bs_nbins' != `cb_nbins' {
            test_fail "binscatter X: `testname'" "bin count mismatch: binscatter=`bs_nbins' cbinscatter=`cb_nbins'"
            test_fail "binscatter Y: `testname'" "bin count mismatch"
            cap erase _bs_test.csv
            cap erase _bs_test.do
            mat drop _bs_data _cb_data
            exit
        }

        * Direct comparison
        local min_sf_x = 15
        forval i = 1/`bs_nbins' {
            local val1 = _bs_data[`i', 1]
            local val2 = _cb_data[`i', 3]
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_x' {
                local min_sf_x = `sf'
            }
        }

        if `min_sf_x' >= `min_sigfigs' {
            local sf_fmt : display %4.1f `min_sf_x'
            test_pass "binscatter X: `testname' (`sf_fmt' sf)"
        }
        else {
            local sf_fmt : display %4.1f `min_sf_x'
            test_fail "binscatter X: `testname'" "sigfigs=`sf_fmt', need `min_sigfigs'"
        }

        * Compare Y values
        local min_sf_y = 15
        forval i = 1/`bs_nbins' {
            local val1 = _bs_data[`i', 2]
            local val2 = _cb_data[`i', 4]
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_y' {
                local min_sf_y = `sf'
            }
        }

        if `min_sf_y' >= `min_sigfigs' {
            local sf_fmt : display %4.1f `min_sf_y'
            test_pass "binscatter Y: `testname' (`sf_fmt' sf)"
        }
        else {
            local sf_fmt : display %4.1f `min_sf_y'
            test_fail "binscatter Y: `testname'" "sigfigs=`sf_fmt', need `min_sigfigs'"
        }

        cap erase _bs_test.csv
        cap erase _bs_test.do
        mat drop _bs_data _cb_data
    end
}

* SECTION 2: Basic functionality
print_section "Basic Functionality"

* Simple two-variable scatter
sysuse auto, clear
cbinscatter price mpg, nograph
if e(N) == _N & e(nquantiles) == 20 test_pass "basic scatter (default 20 bins)"
else test_fail "basic scatter" "N=`=e(N)' expected `=_N', nquantiles=`=e(nquantiles)' expected 20"

* Check e() results populated
sysuse auto, clear
cbinscatter price mpg, nograph
capture confirm matrix e(bindata)
if _rc == 0 test_pass "e(bindata) matrix populated"
else test_fail "e(bindata)" "matrix not found"

capture confirm matrix e(coefs)
if _rc == 0 test_pass "e(coefs) matrix populated (linear fit)"
else test_fail "e(coefs)" "matrix not found"

* Compare actual bin X/Y values against binscatter
if `binscatter_installed' {
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, nquantiles(20) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, nquantiles(20) nograph) ///
        nbins(20) yvar(price) xvar(mpg) testname("basic 20 bins")
}

* SECTION 3: nquantiles option
print_section "nquantiles Option"

sysuse auto, clear
cbinscatter price mpg, nquantiles(10) nograph
if e(nquantiles) == 10 test_pass "nquantiles(10)"
else test_fail "nquantiles(10)" "got `=e(nquantiles)'"

sysuse auto, clear
cbinscatter price mpg, nquantiles(50) nograph
if e(nquantiles) == 50 test_pass "nquantiles(50)"
else test_fail "nquantiles(50)" "got `=e(nquantiles)'"

sysuse auto, clear
cbinscatter price mpg, nquantiles(5) nograph
if e(nquantiles) == 5 test_pass "nquantiles(5)"
else test_fail "nquantiles(5)" "got `=e(nquantiles)'"

* Compare actual bin X/Y values for different nquantiles
if `binscatter_installed' {
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, nquantiles(5) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, nquantiles(5) nograph) ///
        nbins(5) yvar(price) xvar(mpg) testname("nquantiles(5)")

    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("nquantiles(10)")

    * nquantiles(50) on larger dataset (auto has only 21 unique mpg values)
    clear
    set seed 12345
    set obs 2000
    gen x = rnormal()
    gen y = 2*x + rnormal()
    * CSV serialization in binscatter uses a limited-width numeric format.
    * Near-zero means lose significant digits there. Use its actual bin IDs
    * and double-precision Stata means to test the numerical result instead.
    tempvar reference_bin xtile_bin
    tempname reference_x reference_y candidate_x candidate_y
    binscatter y x, nquantiles(50) genxq(`reference_bin')
    xtile `xtile_bin' = x, nq(50)
    assert `reference_bin' == `xtile_bin'
    test_pass "binscatter nquantiles(50): assignments match xtile"
    preserve
    recast double x y
    collapse (mean) x y, by(`reference_bin')
    mkmat x, matrix(`reference_x')
    mkmat y, matrix(`reference_y')
    restore
    cbinscatter y x, nquantiles(50) nograph
    matrix `candidate_x' = e(bindata)[1...,3]
    matrix `candidate_y' = e(bindata)[1...,4]
    assert_matrix_equal `reference_x' `candidate_x' $DEFAULT_SIGFIGS "binscatter X: nquantiles(50)"
    assert_matrix_equal `reference_y' `candidate_y' $DEFAULT_SIGFIGS "binscatter Y: nquantiles(50)"
}

* Test invalid nquantiles (should error)
sysuse auto, clear
capture cbinscatter price mpg, nquantiles(1) nograph
if _rc != 0 test_pass "nquantiles(1) rejected (rc=`=_rc')"
else test_fail "nquantiles(1)" "should have been rejected"

* SECTION 4: discrete option
print_section "discrete Option"

* Create discrete x variable
sysuse auto, clear
count if !missing(price, rep78)
local expected_N = r(N)
cbinscatter price rep78, discrete nograph
if e(N) == `expected_N' test_pass "discrete with rep78 (5 values)"
else test_fail "discrete" "returned no observations"

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
if `num_bins' == 10 test_pass "discrete creates one bin per unique value"
else test_fail "discrete creates one bin per unique value" "got `num_bins' bins, expected 10"

* Compare discrete bin X/Y values against binscatter
if `binscatter_installed' {
    clear
    set seed 12345
    set obs 1000
    gen x = runiformint(1, 10)
    gen y = 2*x + rnormal()
    compare_binscatter_bins, ///
        bs_cmd(binscatter y x, discrete savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter y x, discrete nograph) ///
        nbins(10) yvar(y) xvar(x) testname("discrete 10 values")

    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price rep78, discrete savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78, discrete nograph) ///
        nbins(5) yvar(price) xvar(rep78) testname("discrete rep78")
}

* Without option discrete, binscatter bins an x with at most nquantiles
* distinct values one bin per value (pooled count; each by-group's own
* values), weighted or not. rep78 has 5 values with uneven frequencies.
if `binscatter_installed' {
    * Bins of every by-group, from binscatter's saved CSV (columns
    * <x>_by<g>, <y>_by<g>) and from e(bindata), compared row by row
    capture program drop _cbin_by_compare
    program define _cbin_by_compare
        syntax, bs_cmd(string asis) cb_cmd(string asis) yvar(name) xvar(name) ///
            ngroups(integer) testname(string)
        preserve
        quietly `bs_cmd'
        quietly import delimited _bs_test.csv, clear asdouble
        tempfile wide long
        quietly save `wide'
        forvalues g = 1/`ngroups' {
            quietly use `wide', clear
            keep `xvar'_by`g' `yvar'_by`g'
            rename (`xvar'_by`g' `yvar'_by`g') (bs_x bs_y)
            quietly drop if missing(bs_x)
            gen int bs_g = `g'
            if `g' > 1 quietly append using `long'
            quietly save `long', replace
        }
        sort bs_g bs_x
        mkmat bs_g bs_x bs_y, mat(_bs_data)
        restore
        `cb_cmd'
        mata: _cb = st_matrix("e(bindata)"); _cb = _cb[selectindex(_cb[., 2] :!= .), (1, 3, 4)]
        mata: _cb = sort(_cb, (1, 2)); st_matrix("_cb_data", _cb)
        mata: mata drop _cb
        local ok = rowsof(_bs_data) == rowsof(_cb_data)
        local min_sf = 15
        if `ok' {
            forvalues r = 1/`=rowsof(_bs_data)' {
                if _bs_data[`r', 1] != _cb_data[`r', 1] local ok 0
                forvalues c = 2/3 {
                    sigfigs `=_bs_data[`r', `c']' `=_cb_data[`r', `c']'
                    local min_sf = min(`min_sf', r(sigfigs))
                }
            }
        }
        if `ok' & `min_sf' >= $DEFAULT_SIGFIGS test_pass "`testname'"
        else test_fail "`testname'" "rows binscatter=`=rowsof(_bs_data)' cbinscatter=`=rowsof(_cb_data)', min sigfigs=`min_sf'"
        capture erase _bs_test.csv
        capture erase _bs_test.do
        matrix drop _bs_data _cb_data
    end

    sysuse auto, clear
    gen int fw = 1 + mod(_n, 4)
    compare_binscatter_bins, ///
        bs_cmd(binscatter price rep78, savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78, nograph) ///
        nbins(5) yvar(price) xvar(rep78) testname("automatic discrete rep78")
    compare_binscatter_bins, ///
        bs_cmd(binscatter price rep78 [aw=weight], savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78 [aw=weight], nograph) ///
        nbins(5) yvar(price) xvar(rep78) testname("automatic discrete rep78, aweights")
    compare_binscatter_bins, ///
        bs_cmd(binscatter price rep78 [fw=fw], savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78 [fw=fw], nograph) ///
        nbins(5) yvar(price) xvar(rep78) testname("automatic discrete rep78, fweights")
    compare_binscatter_bins, ///
        bs_cmd(binscatter price rep78, nquantiles(5) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78, nquantiles(5) nograph) ///
        nbins(5) yvar(price) xvar(rep78) testname("automatic discrete rep78, nquantiles(5)")
    _cbin_by_compare, bs_cmd(binscatter price rep78, by(foreign) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78, by(foreign) nograph) ///
        yvar(price) xvar(rep78) ngroups(2) testname("automatic discrete rep78, by(foreign)")
    _cbin_by_compare, bs_cmd(binscatter price rep78 [aw=weight], by(foreign) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price rep78 [aw=weight], by(foreign) nograph) ///
        yvar(price) xvar(rep78) ngroups(2) testname("automatic discrete rep78, by(foreign) aweights")
    * Ten values: discrete with nquantiles(10), quantile bins with nquantiles(9)
    clear
    set seed 4021
    set obs 600
    gen double x = cond(runiform() < .4, 1, runiformint(2, 10))
    gen double y = x + rnormal()
    gen byte g = runiform() < .5
    compare_binscatter_bins, ///
        bs_cmd(binscatter y x, nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter y x, nquantiles(10) nograph) ///
        nbins(10) yvar(y) xvar(x) testname("automatic discrete at nquantiles = distinct values")
    compare_binscatter_bins, ///
        bs_cmd(binscatter y x, nquantiles(9) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter y x, nquantiles(9) nograph) ///
        nbins(9) yvar(y) xvar(x) testname("quantile bins at nquantiles < distinct values")
    _cbin_by_compare, bs_cmd(binscatter y x, by(g) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter y x, by(g) nograph) ///
        yvar(y) xvar(x) ngroups(2) testname("automatic discrete 10 values, by()")
}

* SECTION 5: controls option
print_section "controls Option"

sysuse auto, clear
cbinscatter price mpg, controls(weight) nograph
if "`e(controls)'" == "weight" test_pass "controls(weight)"
else test_fail "controls(weight)" "e(controls)=`e(controls)'"

sysuse auto, clear
cbinscatter price mpg, controls(weight length) nograph
if "`e(controls)'" == "weight length" test_pass "controls(weight length)"
else test_fail "controls(weight length)" "e(controls)=`e(controls)'"

* Many controls
sysuse auto, clear
cbinscatter price mpg, controls(weight length turn displacement) nograph
if e(N) == _N test_pass "controls (4 variables)"
else test_fail "controls (4 variables)" "returned no observations"

* Compare bin X/Y values with controls against binscatter
if `binscatter_installed' {
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls(weight) values")

    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight length) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight length) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls(weight length) values")

    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight length turn displacement) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight length turn displacement) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls(4 vars) values")
}

* SECTION 6: absorb option
print_section "absorb Option"

sysuse auto, clear
cbinscatter price mpg, absorb(foreign) nograph
if "`e(absorb)'" == "foreign" test_pass "absorb(foreign)"
else test_fail "absorb(foreign)" "e(absorb)=`e(absorb)'"

sysuse auto, clear
cbinscatter price mpg, absorb(foreign rep78) nograph
if "`e(absorb)'" == "foreign rep78" test_pass "absorb(foreign rep78)"
else test_fail "absorb(foreign rep78)" "e(absorb)=`e(absorb)'"

* Absorb with controls
sysuse auto, clear
cbinscatter price mpg, controls(weight) absorb(foreign) nograph
if "`e(controls)'" == "weight" & "`e(absorb)'" == "foreign" test_pass "controls + absorb"
else test_fail "controls + absorb" "controls=`e(controls)' absorb=`e(absorb)'"

* Compare bin X/Y values with absorb against binscatter
* Note: binscatter only supports single absorb variable
if `binscatter_installed' {
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("absorb(foreign) values")

    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight) absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight) absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls + absorb values")
}

* SECTION 7: by() option
print_section "by() Option"

sysuse auto, clear
cbinscatter price mpg, by(foreign) nograph
if e(num_groups) == 2 test_pass "by(foreign) creates 2 groups"
else test_fail "by(foreign)" "num_groups=`=e(num_groups)' expected 2"

* by with more groups
sysuse auto, clear
cbinscatter price mpg, by(rep78) nograph
if e(num_groups) >= 3 test_pass "by(rep78) creates multiple groups"
else test_fail "by(rep78)" "num_groups=`=e(num_groups)'"

* by with controls
sysuse auto, clear
cbinscatter price mpg, controls(weight) by(foreign) nograph
if e(num_groups) == 2 & "`e(controls)'" == "weight" test_pass "by + controls"
else test_fail "by + controls" "num_groups=`=e(num_groups)' controls=`e(controls)'"

* by with absorb
sysuse auto, clear
cbinscatter price mpg, absorb(rep78) by(foreign) nograph
if e(num_groups) == 2 & "`e(absorb)'" == "rep78" test_pass "by + absorb"
else test_fail "by + absorb" "num_groups=`=e(num_groups)' absorb=`e(absorb)'"

* SECTION 8: method(binsreg) option
* Compares actual bin X and Y coordinates between cbinscatter method(binsreg)
* and binsreg.ado. Uses price as X (continuous, 73 unique values in auto).
* Compare at the standard significant-figure threshold on identical bins.
* The reference adapter below uses explicit knots between xtile groups to
* avoid binsreg's cutpoint-to-macro rounding, double storage for saved dots,
* and keepsingletons to match cbinscatter's estimation and centering sample.
print_section "method(binsreg) Option"

* Check if binsreg is installed
capture which binsreg
local binsreg_installed = (_rc == 0)

if `binsreg_installed' {

    * Compare the same partition at full precision using binsreg's public options.
    * Its default helper rounds cutpoints through local macros; mid-gap knots
    * preserve xtile's partition. Assert that equality before fitting anything.
    * For degree-zero dots, moving knots inside empty gaps also preserves the
    * fitted model. This adapter is deliberately restricted to dots(0 0).
    * cbinscatter keeps FE singletons and centers over the full eligible sample;
    * reghdfe's default drops them and can shift every adjusted bin vertically.
    capture program drop ctools_binsreg_reference
    program define ctools_binsreg_reference
        syntax varlist(min=2 numeric) [if] [in] [aw fw pw], ///
            NBINS(integer) DOTS(numlist integer min=2 max=2) [ABSorb(varlist) *]
        assert "`dots'" == "0 0"
        marksample touse
        if "`absorb'" != "" markout `touse' `absorb'
        local x : word 2 of `varlist'
        local wt ""
        if "`weight'" != "" local wt "[`weight'`exp']"
        tempvar q refq
        tempname bounds knots
        quietly xtile `q' = `x' if `touse' `wt', nq(`nbins')
        preserve
        quietly keep if `touse'
        quietly collapse (min) lo=`x' (max) hi=`x', by(`q')
        assert _N == `nbins'
        mkmat lo hi, matrix(`bounds')
        restore
        matrix `knots' = J(`nbins'+1, 1, .)
        matrix `knots'[1,1] = `bounds'[1,1]
        matrix `knots'[`nbins'+1,1] = `bounds'[`nbins',2]
        local positions ""
        forvalues b = 2/`nbins' {
            matrix `knots'[`b',1] = `bounds'[`b'-1,2]/2 + `bounds'[`b',1]/2
            local cut : display %24.17g `knots'[`b',1]
            local positions "`positions' `cut'"
        }
        quietly gen long `refq' = .
        quietly binsreg_irecode `x' if `touse', knotmat(`knots') binid(`refq')
        assert `refq' == `q' if `touse'
        local absopt ""
        if "`absorb'" != "" local absopt "absorb(`absorb') reghdfeopt(keepsingletons tolerance(1e-12))"
        * binsreg creates its saved dot variables using Stata's default type.
        * Restore the caller's setting even when the reference command fails.
        local oldtype "`c(type)'"
        set type double
        capture noisily binsreg `varlist' if `touse' `wt', ///
            dots(`dots') binspos(`positions') `absopt' `options'
        local rc = _rc
        set type `oldtype'
        exit `rc'
    end

    * -------------------------------------------------------------------------
    * Helper: compare_binsreg_bins
    *   Runs binsreg and cbinscatter, compares bin X and Y values.
    *   Arguments:
    *     binsreg_cmd  - full binsreg command (must include savedata(_binsreg_test) replace)
    *     cbins_cmd    - full cbinscatter command (must include method(binsreg) nograph)
    *     nbins        - number of bins requested
    *     testname     - test label
    * -------------------------------------------------------------------------
    capture program drop compare_binsreg_bins
    program define compare_binsreg_bins
        syntax, binsreg_cmd(string asis) cbins_cmd(string asis) nbins(integer) testname(string) [minsf(real $DEFAULT_SIGFIGS)]

        local min_sigfigs = `minsf'

        * Run binsreg and extract bin coordinates
        preserve
        gettoken refcommand refargs : binsreg_cmd
        assert "`refcommand'" == "binsreg"
        capture erase _binsreg_test.dta
        quietly ctools_binsreg_reference `refargs'
        use _binsreg_test, clear
        mkmat dots_fit, mat(_br_y)
        mkmat dots_x, mat(_br_x)
        local br_nbins = _N
        restore

        * Run cbinscatter
        `cbins_cmd'
        mat _cb = e(bindata)

        * Count non-missing cbinscatter bins
        local cb_nbins = 0
        local nrows = rowsof(_cb)
        forval r = 1/`nrows' {
            if _cb[`r', 2] != . {
                local cb_nbins = `cb_nbins' + 1
            }
        }

        * Check bin counts match before comparing values
        if `br_nbins' != `cb_nbins' {
            test_fail "binsreg X: `testname'" "bin count mismatch: binsreg=`br_nbins' cbinscatter=`cb_nbins'"
            test_fail "binsreg Y: `testname'" "bin count mismatch: binsreg=`br_nbins' cbinscatter=`cb_nbins'"
            cap erase _binsreg_test.dta
            mat drop _br_y _br_x _cb
            exit
        }

        * Compare X values (bin centers)
        local min_sf_x = 15
        forval i = 1/`br_nbins' {
            local val1 = _br_x[`i', 1]
            local val2 = _cb[`i', 3]
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_x' {
                local min_sf_x = `sf'
            }
        }

        if `min_sf_x' >= `min_sigfigs' {
            local sf_fmt : display %4.1f `min_sf_x'
            test_pass "binsreg X: `testname' (`sf_fmt' sf)"
        }
        else {
            local sf_fmt : display %4.1f `min_sf_x'
            test_fail "binsreg X: `testname'" "sigfigs=`sf_fmt', need `min_sigfigs'"
        }

        * Compare Y values (conditional means)
        local min_sf_y = 15
        forval i = 1/`br_nbins' {
            local val1 = _br_y[`i', 1]
            local val2 = _cb[`i', 4]
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_y' {
                local min_sf_y = `sf'
            }
        }

        if `min_sf_y' >= `min_sigfigs' {
            local sf_fmt : display %4.1f `min_sf_y'
            test_pass "binsreg Y: `testname' (`sf_fmt' sf)"
        }
        else {
            local sf_fmt : display %4.1f `min_sf_y'
            test_fail "binsreg Y: `testname'" "sigfigs=`sf_fmt', need `min_sigfigs'"
        }

        cap erase _binsreg_test.dta
        mat drop _br_y _br_x _cb
    end

    * Basic binsreg method (no controls, no absorb)
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("basic")

    * binsreg with controls
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price length, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, controls(length) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("controls")

    * binsreg with absorb
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price, absorb(foreign) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, absorb(foreign) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("absorb")

    * binsreg with controls + absorb
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price length, absorb(foreign) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, controls(length) absorb(foreign) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("controls + absorb")
}
else {
    test_fail "method(binsreg) comparison" "binsreg not installed - cannot validate"
}

* SECTION 9: Comparison with binscatter.ado
* Compares actual bin X and Y coordinates between cbinscatter and binscatter.
* binscatter saves bin means as CSV via savedata().
* cbinscatter stores bin means in e(bindata)[,3] (x) and e(bindata)[,4] (y).
* MEAN OFFSET: binscatter adds sample means back to residuals after
* partialling out controls/absorb; cbinscatter stores raw residuals. The
* helper computes the mean difference between the two outputs and shifts
* cbinscatter values to match. For the basic case (no controls/absorb),
* the offset is ~0.
print_section "Comparison vs binscatter.ado"

if `binscatter_installed' {

    * Basic scatter
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("basic")

    * With controls
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls(weight)")

    * With absorb
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("absorb(foreign)")

    * Controls + absorb
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight) absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight) absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls + absorb")

    * Multiple controls
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, controls(weight length) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, controls(weight length) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("controls(weight length)")

    * Different nquantiles
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg, nquantiles(20) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg, nquantiles(20) nograph) ///
        nbins(20) yvar(price) xvar(mpg) testname("nquantiles(20)")

    * Different y/x variables
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter weight displacement, nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter weight displacement, nquantiles(10) nograph) ///
        nbins(10) yvar(weight) xvar(displacement) testname("weight~displacement")

    * Controls + absorb + different variables
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter weight displacement, controls(length) absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter weight displacement, controls(length) absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(weight) xvar(displacement) testname("weight~displ ctrl+absorb")
}
else {
    test_fail "binscatter comparison" "binscatter not installed - cannot validate"
}

* SECTION 10: method(binsreg) comparison with binsreg.ado
* Compares actual bin X and Y coordinates between cbinscatter method(binsreg)
* and the binsreg package (Cattaneo, Crump, Farrell, Feng).
* binsreg with dots(0 0) produces degree-0 (constant) fits = conditional means,
* which matches what cbinscatter method(binsreg) computes.
* binsreg savedata() creates a .dta with dots_x (bin center) and dots_fit (y).
* cbinscatter stores bin centers in e(bindata)[,3] and y means in e(bindata)[,4].
* IMPORTANT: Use continuous X variables. binsreg with integer/discrete X
* (e.g. mpg with 21 unique values) creates one bin per unique value regardless
* of nbins(), making comparison impossible. Use price, weight, ln_wage, etc.
* Compare at the standard significant-figure threshold on identical bins.
* The reference adapter below uses explicit knots between xtile groups to
* avoid binsreg's cutpoint-to-macro rounding, double storage for saved dots,
* and keepsingletons to match cbinscatter's estimation and centering sample.
print_section "method(binsreg) vs binsreg.ado"

* Check if binsreg is installed
capture which binsreg
local binsreg_installed = (_rc == 0)

if `binsreg_installed' {

    * =====================================================================
    * AUTO DATASET - use price as X (73 unique values among 74 obs)
    * =====================================================================

    * auto: one control
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price length, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, controls(length) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: ctrl(length)")

    * auto: multiple controls
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price length mpg, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, controls(length mpg) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: ctrl(length mpg)")

    * auto: absorb only
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price, absorb(foreign) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, absorb(foreign) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: absorb(foreign)")

    * auto: controls + absorb
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price length, absorb(foreign) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price, controls(length) absorb(foreign) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: ctrl + absorb")

    * auto: different y/x with controls
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg displacement price length, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter displacement price, controls(length) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: displ~price ctrl(length)")

    * auto: if condition with controls
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price length if price > 4000, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price if price > 4000, controls(length) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: if price>4000 + ctrl")

    * auto: aweights + controls
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price mpg [aw=length], nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price [aw=length], controls(mpg) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: aw + ctrl")

    * auto: aweights + absorb
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price [aw=length], absorb(foreign) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price [aw=length], absorb(foreign) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: aw + absorb")

    * auto: aweights + controls + absorb
    sysuse auto, clear
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg weight price mpg [aw=length], absorb(foreign) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter weight price [aw=length], controls(mpg) absorb(foreign) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("auto: aw + ctrl + absorb")

    * =====================================================================
    * NLSWORK DATASET - pre-drop missing values for consistent samples
    * Use ln_wage as X (continuous with many unique values)
    * =====================================================================

    * nlswork: with control (tenure)
    ctools_fixture nlswork, clear
    keep in 1/5000
    drop if missing(hours) | missing(ln_wage) | missing(tenure)
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg hours ln_wage tenure, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter hours ln_wage, controls(tenure) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("nlswork: ctrl(tenure)")

    * nlswork: absorb(idcode) individual FE
    ctools_fixture nlswork, clear
    keep in 1/5000
    drop if missing(hours) | missing(ln_wage)
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg hours ln_wage, absorb(idcode) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter hours ln_wage, absorb(idcode) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("nlswork: absorb(idcode)")

    * nlswork: controls + absorb
    ctools_fixture nlswork, clear
    keep in 1/5000
    drop if missing(hours) | missing(ln_wage) | missing(tenure)
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg hours ln_wage tenure, absorb(idcode) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter hours ln_wage, controls(tenure) absorb(idcode) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("nlswork: ctrl + absorb(idcode)")

    * nlswork: two-way FE
    ctools_fixture nlswork, clear
    keep in 1/5000
    drop if missing(hours) | missing(ln_wage)
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg hours ln_wage, absorb(idcode year) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter hours ln_wage, absorb(idcode year) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("nlswork: two-way FE")

    * nlswork: controls + two-way FE
    ctools_fixture nlswork, clear
    keep in 1/5000
    drop if missing(hours) | missing(ln_wage) | missing(tenure)
    compare_binsreg_bins, ///
        binsreg_cmd(binsreg hours ln_wage tenure, absorb(idcode year) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter hours ln_wage, controls(tenure) absorb(idcode year) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("nlswork: ctrl + two-way FE")

    * =====================================================================
    * SYNTHETIC DATA (controlled DGP, continuous X, no missing values)
    * =====================================================================

    * synthetic: controls only
    clear
    set seed 12345
    set obs 2000
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x, controls(w) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: controls only")

    * synthetic: absorb only
    clear
    set seed 23456
    set obs 2000
    gen id = ceil(_n/200)
    gen x = rnormal()
    gen y = 2 * x + id * 0.5 + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x, absorb(id) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x, absorb(id) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: absorb only")

    * synthetic: controls + absorb
    clear
    set seed 34567
    set obs 2000
    gen id = ceil(_n/200)
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + id * 0.5 + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w, absorb(id) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x, controls(w) absorb(id) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: controls + absorb")

    * synthetic: multiple controls
    clear
    set seed 45678
    set obs 2000
    gen x = rnormal()
    gen w1 = 0.3 * x + rnormal()
    gen w2 = -0.2 * x + 0.5 * w1 + rnormal()
    gen w3 = rnormal()
    gen y = 2 * x + 0.8 * w1 - 0.5 * w2 + 0.3 * w3 + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w1 w2 w3, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x, controls(w1 w2 w3) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: 3 controls")

    * synthetic: many FE levels (200)
    clear
    set seed 56789
    set obs 5000
    gen id = ceil(_n/25)
    gen x = rnormal()
    gen y = 2 * x + id * 0.1 + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x, absorb(id) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x, absorb(id) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: 200 FE levels")

    * synthetic: two-way FE
    clear
    set seed 98765
    set obs 2000
    gen id1 = ceil(_n/200)
    gen id2 = mod(_n - 1, 20) + 1
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.7 * w + id1 * 0.3 + id2 * 0.1 + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w, absorb(id1 id2) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x, controls(w) absorb(id1 id2) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: ctrl + two-way FE")

    * synthetic: aweights
    clear
    set seed 67890
    set obs 2000
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + rnormal()
    gen wt = abs(rnormal()) + 0.1

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w [aw=wt], nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x [aw=wt], controls(w) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: aw + ctrl")

    * synthetic: fweights + absorb
    clear
    set seed 89012
    set obs 1000
    gen id = ceil(_n/100)
    gen x = rnormal()
    gen y = 2 * x + id * 0.5 + rnormal()
    gen fw = ceil(abs(rnormal()) * 3) + 1

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x [fw=fw], absorb(id) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x [fw=fw], absorb(id) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: fw + absorb")

    * synthetic: if condition
    clear
    set seed 90123
    set obs 2000
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + rnormal()

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w if x > 0, nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x if x > 0, controls(w) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: if condition")

    * synthetic: aweights + absorb
    clear
    set seed 13579
    set obs 2000
    gen id = ceil(_n/200)
    gen x = rnormal()
    gen y = 2 * x + id * 0.5 + rnormal()
    gen wt = abs(rnormal()) + 0.1

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x [aw=wt], absorb(id) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x [aw=wt], absorb(id) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: aw + absorb")

    * synthetic: aweights + controls + absorb
    clear
    set seed 24680
    set obs 2000
    gen id = ceil(_n/200)
    gen x = rnormal()
    gen w = 0.5 * x + rnormal()
    gen y = 2 * x + 0.8 * w + id * 0.5 + rnormal()
    gen wt = abs(rnormal()) + 0.1

    compare_binsreg_bins, ///
        binsreg_cmd(binsreg y x w [aw=wt], absorb(id) nbins(20) dots(0 0) savedata(_binsreg_test) replace) ///
        cbins_cmd(cbinscatter y x [aw=wt], controls(w) absorb(id) nquantiles(20) method(binsreg) nograph) ///
        nbins(20) testname("synthetic: aw + ctrl + absorb")

    * Clean up helper
    capture program drop compare_binsreg_bins
}
else {
    test_fail "binsreg comparison" "binsreg not installed - cannot validate"
}

* SECTION 10: linetype option
print_section "linetype Option"

sysuse auto, clear
cbinscatter price mpg, linetype(none) nograph
if e(N) == _N test_pass "linetype(none)"
else test_fail "linetype(none)" "returned no observations"

sysuse auto, clear
cbinscatter price mpg, linetype(linear) nograph
capture confirm matrix e(coefs)
if _rc == 0 test_pass "linetype(linear)"
else test_fail "linetype(linear)" "no coefs matrix"

sysuse auto, clear
cbinscatter price mpg, linetype(qfit) nograph
capture confirm matrix e(coefs)
if _rc == 0 test_pass "linetype(qfit)"
else test_fail "linetype(qfit)" "no coefs matrix"

sysuse auto, clear
cbinscatter price mpg, linetype(cubic) nograph
capture confirm matrix e(coefs)
if _rc == 0 test_pass "linetype(cubic)"
else test_fail "linetype(cubic)" "no coefs matrix"

* SECTION 11: Weights
print_section "Weights"

sysuse auto, clear
capture cbinscatter price mpg [aw=weight], nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) == _N test_pass "aweight (N=`=e(N)')"
    else test_fail "aweight" "bindata missing or N=`=e(N)' expected `=_N'"
}
else {
    test_fail "aweight" "rc=`=_rc'"
}

sysuse auto, clear
gen fw = ceil(mpg/5)
capture cbinscatter price mpg [fw=fw], nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) > 0 & e(N) < . test_pass "fweight (N=`=e(N)')"
    else test_fail "fweight" "bindata missing or N=`=e(N)'"
}
else {
    test_fail "fweight" "rc=`=_rc'"
}

sysuse auto, clear
capture cbinscatter price mpg [pw=weight], nograph
if _rc == 0 {
    capture confirm matrix e(bindata)
    if _rc == 0 & e(N) == _N test_pass "pweight (N=`=e(N)')"
    else test_fail "pweight" "bindata missing or N=`=e(N)' expected `=_N'"
}
else {
    test_fail "pweight" "rc=`=_rc'"
}

* Weights with controls
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], controls(length) nograph
if _rc == 0 {
    if "`e(controls)'" == "length" & e(N) == _N test_pass "aweight + controls"
    else test_fail "aweight + controls" "controls=`e(controls)' N=`=e(N)'"
}
else {
    test_fail "aweight + controls" "rc=`=_rc'"
}

* Weights with absorb
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], absorb(foreign) nograph
if _rc == 0 {
    if "`e(absorb)'" == "foreign" & e(N) == _N test_pass "aweight + absorb"
    else test_fail "aweight + absorb" "absorb=`e(absorb)' N=`=e(N)'"
}
else {
    test_fail "aweight + absorb" "rc=`=_rc'"
}

* SECTION 11b: Weighted bin comparison vs binscatter.ado
* Compares weighted cbinscatter bins against binscatter [aw=weight].
* binscatter uses xtile [aw=weight] for weighted quantile binning;
* cbinscatter now uses the same cutpoint-based algorithm.
print_section "Weighted Comparison vs binscatter.ado"

if `binscatter_installed' {

    * Analytic weights - basic
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg [aw=weight], nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg [aw=weight], nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("aweight basic")

    * Analytic weights + controls
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg [aw=weight], controls(length) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg [aw=weight], controls(length) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("aweight + controls")

    * Analytic weights + absorb
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg [aw=weight], absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg [aw=weight], absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("aweight + absorb")

    * Frequency weights - basic
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg [fw=rep78], nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg [fw=rep78], nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("fweight basic")

    * Analytic weights + controls + absorb
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg [aw=weight], controls(length) absorb(foreign) nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg [aw=weight], controls(length) absorb(foreign) nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("aweight ctrl + absorb")
}
else {
    test_fail "weighted binscatter comparison" "binscatter not installed - cannot validate"
}

* SECTION 12: savedata option
print_section "savedata Option"

sysuse auto, clear
capture erase _cbins_test.dta
capture erase _cbins_test.csv
capture erase _cbins_test.do
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

    if `nobs' > 0 test_pass "savedata file has data"
    else test_fail "savedata file" "empty"
    cap erase _cbins_test.dta
    * binscatter's .csv and .do, written besides the .dta
    capture confirm file "_cbins_test.csv"
    local csv_rc = _rc
    capture confirm file "_cbins_test.do"
    if !`csv_rc' & !_rc test_pass "savedata also writes binscatter's .csv and .do"
    else test_fail "savedata .csv and .do" "not written"
    cap erase _cbins_test.csv
    cap erase _cbins_test.do
}
else {
    test_fail "savedata" "file not created"
}

* SECTION 13: if/in conditions
print_section "if/in Conditions"

sysuse auto, clear
count if price > 5000
local expected_N = r(N)
cbinscatter price mpg if price > 5000, nograph
if e(N) == `expected_N' test_pass "if condition: N=`=e(N)' matches"
else test_fail "if condition" "N=`=e(N)' expected `expected_N'"

sysuse auto, clear
cbinscatter price mpg in 1/50, nograph
if e(N) == 50 test_pass "in condition: N=50 matches"
else test_fail "in condition" "N=`=e(N)' expected 50"

* Combined if and in
sysuse auto, clear
count if foreign == 0 in 1/60
local expected_N = r(N)
cbinscatter price mpg in 1/60 if foreign == 0, nograph
if e(N) == `expected_N' test_pass "if + in combined: N=`=e(N)' matches"
else test_fail "if + in combined" "N=`=e(N)' expected `expected_N'"

* Compare bin X/Y values with if condition against binscatter
if `binscatter_installed' {
    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg if price > 5000, nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg if price > 5000, nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("if condition values")

    sysuse auto, clear
    compare_binscatter_bins, ///
        bs_cmd(binscatter price mpg if foreign == 0, nquantiles(10) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter price mpg if foreign == 0, nquantiles(10) nograph) ///
        nbins(10) yvar(price) xvar(mpg) testname("if foreign==0 values")
}

* SECTION 14: Missing values
print_section "Missing Values"

* Missing in y
clear
set seed 99
set obs 500
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . if _n <= 50

cbinscatter y x, nograph
if e(N) == 450 & e(N_dropped) == 50 test_pass "missing in y: handled correctly"
else test_fail "missing in y" "N=`=e(N)' N_dropped=`=e(N_dropped)'"

* Missing in x
clear
set seed 100
set obs 500
gen x = runiform()
replace x = . if _n <= 50
gen y = 2*x + rnormal()

cbinscatter y x, nograph
if e(N) == 450 & e(N_dropped) == 50 test_pass "missing in x: handled correctly"
else test_fail "missing in x" "N=`=e(N)' N_dropped=`=e(N_dropped)'"

* Missing in controls
clear
set seed 101
set obs 500
gen x = runiform()
gen w = runiform()
replace w = . if _n <= 50
gen y = 2*x + w + rnormal()

cbinscatter y x, controls(w) nograph
if e(N_dropped) >= 50 test_pass "missing in controls: handled correctly"
else test_fail "missing in controls" "N_dropped=`=e(N_dropped)'"

* SECTION 15: Large datasets
print_section "Large Datasets"

* 50K observations
clear
set seed 200
set obs 50000
gen x = runiform()
gen y = 2*x + rnormal()

cbinscatter y x, nograph
if e(N) == 50000 test_pass "50K observations"
else test_fail "50K observations" "N=`=e(N)'"

* 50K with controls
clear
set seed 201
set obs 50000
gen x = runiform()
gen w1 = runiform()
gen w2 = rnormal()
gen y = 2*x + w1 + w2 + rnormal()

cbinscatter y x, controls(w1 w2) nograph
if e(N) == 50000 test_pass "50K with controls"
else test_fail "50K with controls" "N=`=e(N)'"

* 50K with absorb
clear
set seed 202
set obs 50000
gen id = runiformint(1, 500)
gen x = runiform()
gen y = 2*x + id*0.01 + rnormal()

cbinscatter y x, absorb(id) nograph
if e(N) == 50000 test_pass "50K with absorb (500 FE)"
else test_fail "50K with absorb" "N=`=e(N)'"

* SECTION 16: Panel data (nlswork)
print_section "Panel Data (nlswork)"

ctools_fixture nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, nograph
if e(N) > 0 & e(N) < . test_pass "nlswork basic"
else test_fail "nlswork basic" "N=`=e(N)'"

ctools_fixture nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, controls(tenure) nograph
if e(N) > 0 & e(N) < . test_pass "nlswork with controls"
else test_fail "nlswork with controls" "N=`=e(N)'"

ctools_fixture nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, absorb(idcode) nograph
if e(N) > 0 & e(N) < . test_pass "nlswork with individual FE"
else test_fail "nlswork with individual FE" "N=`=e(N)'"

ctools_fixture nlswork, clear
keep in 1/10000
cbinscatter ln_wage age, absorb(idcode year) nograph
if e(N) > 0 & e(N) < . test_pass "nlswork with two-way FE"
else test_fail "nlswork with two-way FE" "N=`=e(N)'"

* SECTION 17: Stress tests - method(binsreg)
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
    if `nrows' == 15 & e(N) == 2000 & e(nquantiles) == 15 test_pass "binsreg: multiple controls + absorb (bins=`nrows', N=`=e(N)')"
    else test_fail "binsreg: multiple controls + absorb" "bins=`nrows' expected 15, N=`=e(N)' expected 2000"
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
    if `nrows' == 12 & e(N) == 2000 & e(nquantiles) == 12 test_pass "binsreg: multiple FE dimensions (bins=`nrows', N=`=e(N)')"
    else test_fail "binsreg: multiple FE dimensions" "bins=`nrows' expected 12, N=`=e(N)' expected 2000"
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
    if `nrows' == 10 & e(N) == 1500 & e(nquantiles) == 10 test_pass "binsreg: unbalanced FE groups (bins=`nrows', N=`=e(N)')"
    else test_fail "binsreg: unbalanced FE groups" "bins=`nrows' expected 10, N=`=e(N)' expected 1500"
}
else {
    test_fail "binsreg: unbalanced FE groups" "rc=`=_rc'"
}

* SECTION 18: Combination tests
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
        if `nrows' >= 5 test_pass "discrete + by (groups=`=e(num_groups)', rows=`nrows')"
        else test_fail "discrete + by" "bindata has only `nrows' rows, expected >= 5"
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
if _rc == 0 & e(num_groups) == 2 test_pass "by + controls + absorb"
else test_fail "by + controls + absorb" "rc=`=_rc' num_groups=`=e(num_groups)'"

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
if _rc == 0 & e(num_groups) == 2 test_pass "method(binsreg) + by + controls + absorb"
else test_fail "method(binsreg) + by" "rc=`=_rc'"

* weights + by
sysuse auto, clear
capture cbinscatter price mpg [aw=weight], by(foreign) nograph
if _rc == 0 {
    if e(num_groups) == 2 & e(N) == _N test_pass "aweight + by (groups=`=e(num_groups)', N=`=e(N)')"
    else test_fail "aweight + by" "num_groups=`=e(num_groups)' expected 2, N=`=e(N)' expected `=_N'"
}
else {
    test_fail "aweight + by" "rc=`=_rc'"
}

* discrete + controls: the residualized x has more distinct values than
* half the sample, which binscatter refuses (r(134))
sysuse auto, clear
test_error_match, stata_cmd(binscatter price rep78, discrete controls(weight) nodraw) ///
    ctools_cmd(cbinscatter price rep78, discrete controls(weight) nograph) ///
    testname("discrete + controls")

* linetype + controls + absorb
sysuse auto, clear
capture cbinscatter price mpg, controls(weight) absorb(foreign) linetype(qfit) nograph
if _rc == 0 {
    * Verify coefs matrix exists (qfit should produce polynomial coefficients)
    capture confirm matrix e(coefs)
    if _rc == 0 & "`e(controls)'" == "weight" & "`e(absorb)'" == "foreign" test_pass "qfit + controls + absorb"
    else test_fail "qfit + controls + absorb" "coefs rc=`=_rc' controls=`e(controls)' absorb=`e(absorb)'"
}
else {
    test_fail "qfit + controls + absorb" "rc=`=_rc'"
}

* SECTION 19: Edge cases
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
        if `nrows' == 5 test_pass "very few observations (20 obs, 5 bins)"
        else test_fail "very few observations" "bindata rows=`nrows' expected 5"
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
        if `nrows' == 50 test_pass "many bins (50 bins, 100 obs)"
        else test_fail "many bins" "bindata rows=`nrows' expected 50"
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
if _rc == 0 test_pass "collinear control (handled)"
else test_fail "collinear control" "unexpected error rc=`=_rc'"

* Single FE group (absorb with 1 level)
clear
set obs 100
gen id = 1  // Single group
gen x = runiform()
gen y = 2*x + rnormal()

capture cbinscatter y x, absorb(id) nograph
if _rc == 0 | _rc == 2001 test_pass "single FE group (rc=`=_rc')"
else test_fail "single FE group" "unexpected rc=`=_rc'"

* Empty by-group after filtering (regression test for segfault fix)
* This used to crash when by() was used but filtering left only one group
sysuse auto, clear
capture cbinscatter price mpg if foreign == 0, by(foreign) nograph
if _rc == 0 & e(num_groups) == 1 test_pass "empty by-group after filtering"
else test_fail "empty by-group after filtering" "rc=`=_rc' num_groups=`=e(num_groups)'"

* SECTION 20: verbose option
print_section "verbose Option"

sysuse auto, clear
capture cbinscatter price mpg, verbose nograph
if _rc == 0 test_pass "verbose option"
else test_fail "verbose option" "rc=`=_rc'"

* SECTION: Intentional Error Tests
* These tests verify that cbinscatter returns the same error codes as binscatter
* when given invalid inputs or error conditions.
* Note: binscatter is a user-written command (ssc install binscatter).
* If binscatter is not installed, tests compare against expected behavior.
print_section "Intentional Error Tests"

* Variable doesn't exist
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price nonexistent_var, nograph) ctools_cmd(cbinscatter price nonexistent_var, nograph) testname("nonexistent variable")
}
else {
    capture cbinscatter price nonexistent_var, nograph
    if _rc != 0 test_pass "[error] nonexistent variable (rc=`=_rc') [binscatter not installed]"
    else test_fail "[error] nonexistent variable" "should have errored"
}

* String variable as y
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter make mpg, nograph) ctools_cmd(cbinscatter make mpg, nograph) testname("string y variable")
}
else {
    capture cbinscatter make mpg, nograph
    if _rc != 0 test_pass "[error] string y variable (rc=`=_rc') [binscatter not installed]"
    else test_fail "[error] string y variable" "should have errored"
}

* String variable as x
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price make, nograph) ctools_cmd(cbinscatter price make, nograph) testname("string x variable")
}
else {
    capture cbinscatter price make, nograph
    if _rc != 0 test_pass "[error] string x variable (rc=`=_rc') [binscatter not installed]"
    else test_fail "[error] string x variable" "should have errored"
}

* Invalid nquantiles (0)
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price mpg, nquantiles(0) nograph) ctools_cmd(cbinscatter price mpg, nquantiles(0) nograph) testname("nquantiles(0)")
}
else {
    capture cbinscatter price mpg, nquantiles(0) nograph
    if _rc != 0 test_pass "[error] nquantiles(0) (rc=`=_rc') [binscatter not installed]"
    else test_fail "[error] nquantiles(0)" "should have errored"
}

* Invalid nquantiles (negative)
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price mpg, nquantiles(-5) nograph) ctools_cmd(cbinscatter price mpg, nquantiles(-5) nograph) testname("negative nquantiles")
}
else {
    capture cbinscatter price mpg, nquantiles(-5) nograph
    if _rc != 0 test_pass "[error] negative nquantiles (rc=`=_rc') [binscatter not installed]"
    else test_fail "[error] negative nquantiles" "should have errored"
}

* No observations after if
sysuse auto, clear
if `binscatter_installed' {
    test_error_match, stata_cmd(binscatter price mpg if price > 100000, nograph) ctools_cmd(cbinscatter price mpg if price > 100000, nograph) testname("no observations after if")
}
else {
    capture cbinscatter price mpg if price > 100000, nograph
    if _rc != 0 test_pass "[error] no observations after if (rc=`=_rc') [binscatter not installed]"
    else test_fail "[error] no observations after if" "should have errored"
}

* An already-residualized outcome must not stop demeaning the bin indicators.
* Check bin contrasts against an independent regression on explicit FE dummies.
foreach weighted in 0 1 {
    clear
    set obs 48
    gen int obsid = _n - 1
    gen int fe1 = 1 + floor(obsid / 16)
    gen int fe2 = 1 + floor(obsid / 12)
    gen byte x = 1 + mod(mod(obsid, 7), 3)
    gen double raw_y = sin(.7 * obsid) + .05 * obsid
    gen double w = 1 + mod(obsid, 5)
    local wt ""
    if `weighted' local wt "[aw=w]"
    regress raw_y i.fe1 i.fe2 `wt'
    predict double y, residuals
    regress y i.x i.fe1 i.fe2 `wt'
    scalar expected2 = _b[2.x]
    scalar expected3 = _b[3.x]
    cbinscatter y x `wt', discrete absorb(fe1 fe2) method(binsreg) nograph
    matrix batch_bins = e(bindata)
    if abs(batch_bins[2,4] - batch_bins[1,4] - expected2) < 1e-9 & ///
       abs(batch_bins[3,4] - batch_bins[1,4] - expected3) < 1e-9 {
        test_pass "batch FE convergence with residualized y (weighted=`weighted')"
    }
    else {
        test_fail "batch FE convergence with residualized y (weighted=`weighted')" "bin contrasts differ from dummy OLS"
    }

}

* savedata() leaves the data in memory alone
* The bins are built and saved in a frame of their own (Stata 16 or later):
* the data are neither copied nor cleared, and a failed save changes nothing.
print_section "savedata() and the Data in Memory"
capture confirm file "validation/support/rollback_helpers.do"
if !_rc do "validation/support/rollback_helpers.do"
else do "rollback_helpers.do"

sysuse auto, clear
label data "cbinscatter savedata data"
note: dataset note
char _dta[rb_char] dataset characteristic
sort foreign make
tempfile rb_copy rb_bins
tempname saved_bins
save `rb_copy'
_ctools_rb_state
local before `"`r(state)'"'
local frames ""
if c(stata_version) >= 16 {
    quietly frames dir
    local frames `"`r(frames)'"'
}
foreach target in "`rb_bins'" "`c(tmpdir)'/rb_no_such_dir/bins.dta" {
    capture cbinscatter price mpg, nquantiles(5) nograph savedata("`target'")
    local rc = _rc
    if "`target'" == "`rb_bins'" & !`rc' matrix `saved_bins' = e(bindata)
    _ctools_rb_state
    local after `"`r(state)'"'
    local frames_now ""
    if c(stata_version) >= 16 {
        quietly frames dir
        local frames_now `"`r(frames)'"'
    }
    local name = cond("`target'" == "`rb_bins'", "[savedata] data unchanged", "[savedata] failed save changes nothing")
    local ok = cond("`target'" == "`rb_bins'", `rc' == 0, `rc' != 0)
    if `ok' & `"`macval(after)'"' == `"`macval(before)'"' & `"`frames_now'"' == `"`frames'"' {
        test_pass "`name'"
    }
    else test_fail "`name'" "rc=`rc', state or frames changed"
}
use `rb_bins', clear
capture {
    if _N != 5 | c(k) != 5 | rowsof(`saved_bins') != 5 error 9
    local j = 0
    foreach v in by_group bin_id x_mean y_mean n_obs {
        local ++j
        * svmat stores float columns by default: require exactly the values
        * representable in the file's declared storage type, not just rows.
        local type : type `v'
        if "`type'" == "float" assert `v' == float(`saved_bins'[_n, `j'])
        else assert `v' == `saved_bins'[_n, `j']
    }
}
if !_rc test_pass "[savedata] exact bin coordinates and counts saved"
else test_fail "[savedata] bins saved" "dimensions or values differ from e(bindata)"

print_section "binscatter Options, Titles, Discrete Capacity, reportreg, qfit"
if `binscatter_installed' {
    * binscatter's own options are accepted with and without a graph and
    * leave the data unchanged (s2_bin compares their results); with
    * method(binsreg), to which they do not apply, they are refused
    clear
    set seed 4101
    quietly set obs 400
    gen double x = rnormal()
    gen double z = rnormal() + .5*x
    gen double y = 1 + x + z + rnormal()
    gen q = ceil(_n/20)
    gen t = _n
    gen double u = runiform()
    tempfile before
    quietly save `before'
    foreach o in "medians" "noaddmean controls(z)" "rd(0)" "xq(q)" "mcolors(red)" ///
        "lcolors(red)" "randvar(u)" "randcut(.5)" "randn(100)" {
        foreach g in nograph nodraw {
            capture cbinscatter y x, `o' `g'
            local rc = _rc
            capture cf _all using `before'
            if `rc' == 0 & _rc == 0 test_pass "binscatter option `o' accepted (`g')"
            else test_fail "binscatter option `o' accepted (`g')" "rc=`rc'"
        }
        if !inlist("`o'", "mcolors(red)", "lcolors(red)") {
            capture cbinscatter y x, `o' method(binsreg) nograph
            local rc = _rc
            capture cf _all using `before'
            if `rc' == 198 & _rc == 0 test_pass "binscatter option `o' refused with method(binsreg)"
            else test_fail "binscatter option `o' refused with method(binsreg)" "rc=`rc'"
        }
    }
    capture erase _cb_unused.gph
    capture cbinscatter y x, savegraph(_cb_unused.gph) nograph
    local rc = _rc
    capture confirm file _cb_unused.gph
    if `rc' == 198 & _rc test_pass "savegraph() with nograph refused, writes nothing"
    else test_fail "savegraph() with nograph refused, writes nothing" "rc=`rc'"
    * Typos are errors without a graph too, as twoway makes them in binscatter
    foreach o in "bogusoption" "bogus(1)" "xlabel(1, bogus)" {
        test_error_match, stata_cmd(binscatter y x, `o' nodraw) ///
            ctools_cmd(cbinscatter y x, `o' nograph) testname("nograph option `o'")
    }
    * replace and nofastxtile change nothing in cbinscatter; valid graph
    * options are accepted without a graph
    quietly cbinscatter y x, nograph
    matrix _cb_base = e(bindata)
    foreach o in "replace" "nofastxtile" "xlabel(-2(1)2) title(T) ytitle(Y) legend(off)" "name(_cb_unused) saving(_cb_unused)" {
        capture quietly cbinscatter y x, `o' nograph
        local rc = _rc
        if `rc' == 0 & mreldif(e(bindata), _cb_base) == 0 test_pass "option `o' accepted, bins unchanged"
        else test_fail "option `o' accepted, bins unchanged" "rc=`rc'"
    }
    capture confirm file _cb_unused.gph
    if _rc test_pass "saving() without a graph writes nothing"
    else test_fail "saving() without a graph writes nothing" "file written"
    capture graph describe _cb_unused
    if _rc test_pass "name() without a graph stores nothing"
    else test_fail "name() without a graph stores nothing" "graph stored"
    matrix drop _cb_base

    * Title and legend options reach twoway as binscatter passes them:
    * suboptions, several lines, empty titles, legend() with one series.
    * Compared on the graph objects (built with nodraw, never displayed).
    capture program drop _cb_graph_props
    program define _cb_graph_props, rclass
        local p
        foreach item in title.text.arrnels title.style.size.setting ///
            yaxis1.title.text.arrnels xaxis1.title.text.arrnels legend.draw_view.setting {
            local v ""
            capture local v `"`.Graph.`item''"'
            local p `"`p'|`v'"'
        }
        forvalues i = 1/2 {
            foreach item in title yaxis1.title xaxis1.title {
                local v ""
                capture local v `"`.Graph.`item'.text[`i']'"'
                local p `"`p'|`v'"'
            }
        }
        return local props `"`p'"'
    end
    gen bv = mod(_n, 2)
    foreach o in `"title("Line one" "Line two")"' `"title("My title", size(vsmall))"' ///
        `"ytitle("Why", size(large)) xtitle("")"' `"ytitle("")"' `"ytitle("A" "B") xtitle("C" "D")"' ///
        `"title(My plain title)"' `"legend(on)"' `"by(bv) legend(off)"' `"by(bv) legend(rows(2))"' {
        quietly binscatter y x, `o' nodraw
        _cb_graph_props
        local ref `"`r(props)'"'
        capture quietly cbinscatter y x, `o' nodraw
        local rc = _rc
        _cb_graph_props
        local same = `rc' == 0 & `"`r(props)'"' == `"`ref'"'
        * Test names and messages carry no quotation marks
        local label : subinstr local o `"""' "", all
        local got `"`r(props)'"'
        local got : subinstr local got `"""' "", all
        local ref : subinstr local ref `"""' "", all
        if `same' test_pass "graph titles and legend: `label'"
        else test_fail "graph titles and legend: `label'" "rc=`rc' binscatter `ref' cbinscatter `got'"
    }
    capture graph drop Graph
    program drop _cb_graph_props

    * Time-series operators in controls(). binscatter's savedata() CSV
    * holds about 9 significant digits of a value near 1 but fewer of one
    * near 0, so x and y are shifted away from 0 (the same bins, shifted).
    tsset t
    gen double xs = x + 10
    gen double ys = y + 10
    capture cbinscatter ys xs, controls(L.z) nograph
    if _rc test_fail "controls(L.z)" "rc=`=_rc'"
    else compare_binscatter_bins, ///
        bs_cmd(binscatter ys xs, controls(L.z) savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter ys xs, controls(L.z) nograph) ///
        nbins(20) yvar(ys) xvar(xs) testname("controls(L.z)")
    tsset, clear
    drop xs ys

    * discrete: binscatter's (sample size)/2 rule on the residualized x, an
    * exact result size for many values and many by-groups, and the
    * capacity check before the per-bin buffers
    test_error_match, stata_cmd(binscatter y x, discrete controls(z) nodraw) ///
        ctools_cmd(cbinscatter y x, discrete controls(z) nograph) testname("discrete controls(z): too many values")
    test_error_match, stata_cmd(binscatter y x, discrete nodraw) ///
        ctools_cmd(cbinscatter y x, discrete nograph) testname("discrete continuous x: too many values")
    clear
    set seed 4102
    quietly set obs 10000
    gen x = mod(_n, 600)
    gen double y = 10 + x/100 + rnormal()
    capture cbinscatter y x, discrete nograph
    if _rc test_fail "discrete with 600 values" "rc=`=_rc'"
    else compare_binscatter_bins, ///
        bs_cmd(binscatter y x, discrete savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter y x, discrete nograph) ///
        nbins(600) yvar(y) xvar(x) testname("discrete with 600 values")
    clear
    quietly set obs 200000
    gen double x = mod(_n, 70000)
    gen double y = mod(_n, 7)
    tempfile before
    quietly save `before'
    capture cbinscatter y x, discrete nograph
    local rc = _rc
    capture cf _all using `before'
    if `rc' == 503 & _rc == 0 test_pass "discrete beyond c(max_matdim) bins refused, data unchanged"
    else test_fail "discrete beyond c(max_matdim) bins refused, data unchanged" "rc=`rc'"
    clear
    set seed 4103
    quietly set obs 70000
    gen x = mod(_n, 5)
    gen state = mod(floor(_n/5), 140)
    gen double y = x + state/10 + rnormal()
    capture quietly cbinscatter y x, discrete by(state) nograph
    local rc = _rc
    if `rc' == 0 {
        mata: _cb = st_matrix("e(bindata)"); _cb = _cb[selectindex(_cb[., 2] :!= .), (1, 3, 4, 5)]
        mata: st_matrix("_cb_data", sort(_cb, (1, 2)))
        mata: mata drop _cb
        preserve
        collapse (mean) y (count) n=y, by(state x)
        sort state x
        mkmat x y n, mat(_ref)
        restore
        local ok = rowsof(_cb_data) == 700
        local min_sf = 15
        if `ok' {
            forvalues r = 1/700 {
                if _cb_data[`r', 2] != _ref[`r', 1] | _cb_data[`r', 4] != _ref[`r', 3] local ok 0
                sigfigs `=_cb_data[`r', 3]' `=_ref[`r', 2]'
                local min_sf = min(`min_sf', r(sigfigs))
            }
        }
        if `ok' & `min_sf' >= $DEFAULT_SIGFIGS test_pass "discrete by() with 140 groups"
        else test_fail "discrete by() with 140 groups" "rows=`=rowsof(_cb_data)' min sigfigs=`min_sf'"
        matrix drop _cb_data _ref
    }
    else test_fail "discrete by() with 140 groups" "rc=`rc'"

    * reportreg: the fit regressions use cbinscatter's sample, and pweights
    * and iweights work; a failure leaves no scratch results
    sysuse auto, clear
    gen w = 1 + mod(_n, 3)
    capture quietly cbinscatter price rep78, controls(weight) reportreg nograph
    local rc = _rc
    if `rc' == 0 capture confirm matrix e(reportreg_b)
    if `rc' == 0 & _rc local rc = "no e(reportreg_b)"
    if "`rc'" == "0" {
        sigfigs `=e(reportreg_b)[1, colnumb(e(reportreg_b), "_cons")]' `=e(coefs)[1, 1]'
        local sf_cons = r(sigfigs)
        sigfigs `=e(reportreg_b)[1, colnumb(e(reportreg_b), "rep78")]' `=e(coefs)[1, 2]'
        if min(`sf_cons', r(sigfigs)) >= $DEFAULT_SIGFIGS test_pass "reportreg fit on the plugin's sample"
        else test_fail "reportreg fit on the plugin's sample" "sigfigs _cons=`sf_cons' slope=`r(sigfigs)'"
    }
    else test_fail "reportreg fit on the plugin's sample" "rc=`rc'"
    foreach spec in "price mpg [pw=w], controls(weight)" "price mpg [iw=w], absorb(foreign)" ///
        "price mpg [pw=w], absorb(foreign) method(binsreg)" "price mpg [iw=w], absorb(foreign) controls(weight) method(binsreg)" {
        capture quietly cbinscatter `spec' reportreg nograph
        local rc = _rc
        local scratch : all scalars "__cbinscatter*"
        local scratch_m : all matrices "__cbinscatter*"
        if `rc' == 0 & "`scratch'`scratch_m'" == "" test_pass "reportreg with `spec'"
        else test_fail "reportreg with `spec'" "rc=`rc' scratch=`scratch' `scratch_m'"
    }

    * linetype(qfit) with two distinct x values: binscatter's regress keeps
    * x^2 and omits x (e(y1_coefs): x^2, x, _cons; e(coefs): _cons, x, x^2)
    clear
    set seed 4104
    quietly set obs 400
    gen x = runiform() < .4
    gen double y = 2*x + rnormal()
    gen g = mod(_n, 2)
    foreach shift in 0 1 -.5 {
        quietly replace x = x + `shift'
        foreach byopt in "" "by(g)" {
            quietly binscatter y x, linetype(qfit) `byopt' nodraw
            matrix _ref = e(y1_coefs)
            quietly cbinscatter y x, linetype(qfit) `byopt' nograph
            local min_sf = 15
            local ok = rowsof(_ref) == rowsof(e(coefs))
            forvalues gr = 1/`=rowsof(_ref)' {
                forvalues c = 1/3 {
                    local a = _ref[`gr', 4 - `c']
                    local b = e(coefs)[`gr', `c']
                    if `a' == 0 & `b' == 0 continue
                    sigfigs `a' `b'
                    local min_sf = min(`min_sf', r(sigfigs))
                }
            }
            if `ok' & `min_sf' >= $DEFAULT_SIGFIGS test_pass "qfit with two x values, shift `shift' `byopt'"
            else test_fail "qfit with two x values, shift `shift' `byopt'" "min sigfigs=`min_sf'"
            matrix drop _ref
        }
    }
}

* P3-11: by() legend keys are binscatter's byvarname=byvalname: a numeric
* value as a macro expression evaluates it, or its value label; a string's
* text. Compared on the graph objects (built with nodraw, never displayed).
clear
set seed 305
set obs 600
gen double x = rnormal()
gen double y = 1 + x + rnormal()
gen bv = mod(_n, 2)
gen float bf = cond(mod(_n, 3) == 0, .1, cond(mod(_n, 3) == 1, -1.5, 2.25))
gen double bd = cond(mod(_n, 3) == 0, .1, cond(mod(_n, 3) == 1, 1e-7, 123456789))
gen bl = mod(_n, 3)
label define _cb_bll 0 "zero" 1 "one and more"
label values bl _cb_bll
gen str12 bs = cond(mod(_n, 3) == 0, "alpha", cond(mod(_n, 3) == 1, "b c", "z"))
gen long bbig = cond(mod(_n, 2), 2147483000, -5)
foreach b in bv bf bd bl bs bbig {
    local keys_ref
    local keys_ct
    foreach cmd in binscatter cbinscatter {
        local n = 0
        capture quietly `cmd' y x, by(`b') nodraw
        local rc_`cmd' = _rc
        capture local n : di `.Graph.legend.labels.arrnels'
        local keys
        forvalues i = 1/`n' {
            local keys `"`keys'[`.Graph.legend.labels[`i']']"'
        }
        if "`cmd'" == "binscatter" local keys_ref `"`keys'"'
        else local keys_ct `"`keys'"'
    }
    if !`rc_binscatter' & !`rc_cbinscatter' & `"`keys_ref'"' == `"`keys_ct'"' & `"`keys_ref'"' != "" ///
        test_pass "P3-11 by(`b') legend keys `keys_ref'"
    else test_fail "P3-11 by(`b') legend keys" "binscatter `keys_ref' cbinscatter `keys_ct'"
}
capture graph drop Graph
label drop _cb_bll

* P3-12: binscatter switches to discrete bins when the x it stores (predict's
* float residual plus the mean, stored as float) has at most nquantiles
* values; a residual that is exactly constant within absorb() groups, or
* collinear with controls(), carries ~1e-16 noise in double. binscatter's
* float storage of the residualized x and y holds about 7 digits of each
* value, so y and the quantile-binned x4 are shifted away from 0 (the same
* bins, shifted), as for controls(L.z) above.
clear
set seed 412
set obs 600
gen g = ceil(_n / 3)
gen k = mod(_n, 3)
gen double z = rnormal()
gen double x1 = g * .1 + k
gen double x2 = g * .37
gen double x3 = 3 * z + .1
gen double x4 = 1.7 * z + k + 10
gen double y = rnormal() + k + 10
gen w = 1 + mod(_n, 4)
foreach spec in "x1, absorb(g)" "x2, absorb(g)" "x3, controls(z)" "x4, controls(z)" ///
    "x1 [aw=w], absorb(g)" "x1, absorb(g) discrete" "x1, absorb(g) controls(z)" {
    gettoken xv rest : spec, parse(" ,[")
    preserve
    quietly binscatter y `spec' nodraw savedata(_bs_count) replace
    quietly import delimited _bs_count.csv, clear
    local nb = _N
    restore
    capture erase _bs_count.csv
    capture erase _bs_count.do
    compare_binscatter_bins, bs_cmd(binscatter y `spec' savedata(_bs_test) replace) ///
        cb_cmd(cbinscatter y `spec' nograph) nbins(`nb') yvar(y) xvar(`xv') ///
        testname("P3-12 discrete switch: `spec' (`nb' bins)")
}
* by(): the switch counts the pooled x; each group's values are its bins
preserve
quietly binscatter y x1, absorb(g) by(k) nodraw savedata(_bs_count) replace
quietly import delimited _bs_count.csv, clear
local nb = _N
restore
capture erase _bs_count.csv
capture erase _bs_count.do
capture quietly cbinscatter y x1, absorb(g) by(k) nograph
local nc = 0
if !_rc {
    matrix _cb_b = e(bindata)
    forvalues r = 1/`=rowsof(_cb_b)' {
        if _cb_b[`r', 2] < . local ++nc
    }
    matrix drop _cb_b
}
if `nc' == `nb' test_pass "P3-12 discrete switch: by(k) (`nb' bins)"
else test_fail "P3-12 discrete switch: by(k)" "binscatter `nb' bins, cbinscatter `nc'"

* Clean up helpers
capture program drop ctools_binsreg_reference
capture program drop compare_binscatter_bins
capture program drop compare_binsreg_bins

* End of cbinscatter validation
noi print_summary "cbinscatter"
}

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "cbinscatter"
    global CTOOLS_CASE_COMPLETE "cbinscatter"
    exit
}

di as error "Unknown cbinscatter case: `__case'"
exit 198
