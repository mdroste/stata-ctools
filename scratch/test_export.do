********************************************************************************
* test_export.do
* Quick test to verify cexport works correctly on sysuse auto data
********************************************************************************

clear all
set more off

* Add ctools build directory to adopath
adopath ++ "./build"

di as txt _n
di as txt "================================================================================"
di as txt "                    CEXPORT vs EXPORT DELIMITED TEST"
di as txt "================================================================================"

* Load test data
sysuse auto, clear
di as txt "Loaded auto dataset: " _N " observations, " c(k) " variables"

* Define test files
local testfile_stata "scratch/test_export_stata.csv"
local testfile_cexport "scratch/test_export_ctools.csv"

********************************************************************************
* Test 1: Stata's export delimited
********************************************************************************

di as txt _n "{hline 60}"
di as txt "Test 1: Stata's export delimited"
di as txt "{hline 60}"

timer clear 1
timer on 1
export delimited using "`testfile_stata'", replace
timer off 1
quietly timer list 1
local t_stata = r(t1)

di as txt "  Time: " %8.4f `t_stata' " sec"
di as txt "  File: `testfile_stata'"

********************************************************************************
* Test 2: cexport delimited
********************************************************************************

di as txt _n "{hline 60}"
di as txt "Test 2: cexport delimited"
di as txt "{hline 60}"

timer clear 2
timer on 2
capture noisily cexport delimited using "`testfile_cexport'", replace
local cexport_rc = _rc
timer off 2
quietly timer list 2
local t_cexport = r(t2)

if `cexport_rc' != 0 {
    di as error "cexport failed with error code `cexport_rc'"
}
else {
    di as txt "  Time: " %8.4f `t_cexport' " sec"
    di as txt "  File: `testfile_cexport'"
}

********************************************************************************
* Verification: Compare exported files
********************************************************************************

di as txt _n "{hline 60}"
di as txt "Verification: Comparing exported CSV files"
di as txt "{hline 60}"

* Check if cexport file exists
capture confirm file "`testfile_cexport'"
if _rc != 0 {
    di as error "cexport output file does not exist"
    exit
}

* Import both files and compare
di as txt _n "Importing Stata export..."
import delimited using "`testfile_stata'", clear
tempfile stata_data
quietly save `stata_data'
local stata_nobs = _N
local stata_nvars = c(k)

di as txt "Importing cexport output..."
import delimited using "`testfile_cexport'", clear
tempfile cexport_data
quietly save `cexport_data'
local cexport_nobs = _N
local cexport_nvars = c(k)

di as txt _n "Comparison:"
di as txt "  Observations: stata=" `stata_nobs' "  cexport=" `cexport_nobs' ///
    cond(`stata_nobs'==`cexport_nobs', "  [MATCH]", "  [DIFFER]")
di as txt "  Variables:    stata=" `stata_nvars' "  cexport=" `cexport_nvars' ///
    cond(`stata_nvars'==`cexport_nvars', "  [MATCH]", "  [DIFFER]")

* Full data comparison
use `stata_data', clear
capture cf _all using `cexport_data'
if _rc == 0 {
    di as txt _n "  Data verification: {result:PASS} - Files are identical"
    local verification "PASS"
}
else {
    di as txt _n "  Data verification: {error:DIFFER} - Files have differences"
    local verification "DIFFER"

    * Show which variables differ
    di as txt _n "  Checking individual variables:"
    foreach var of varlist * {
        capture cf `var' using `cexport_data'
        if _rc == 0 {
            di as txt "    `var': [MATCH]"
        }
        else {
            di as txt "    `var': [DIFFER]"
        }
    }
}

********************************************************************************
* Summary
********************************************************************************

di as txt _n
di as txt "================================================================================"
di as txt "                              SUMMARY"
di as txt "================================================================================"
di as txt "  export delimited time: " %8.4f `t_stata' " sec"
if `cexport_rc' == 0 {
    di as txt "  cexport time:          " %8.4f `t_cexport' " sec"
    di as txt "  Verification:          `verification'"
}
else {
    di as txt "  cexport:               FAILED (rc=`cexport_rc')"
}
di as txt "================================================================================"

* Clean up test files
capture erase "`testfile_stata'"
capture erase "`testfile_cexport'"
di as txt _n "Test files cleaned up."
