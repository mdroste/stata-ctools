********************************************************************************
* benchmark_cimport.do
* Speed benchmark comparing import delimited and cimport performance
* 50M observations, 10 variables of varying types
********************************************************************************

clear all
set more off

* Add ctools build directory to adopath
adopath ++ "./build"

********************************************************************************
* Configuration
********************************************************************************

local N = 50000000   // 50 million observations
local testfile "benchmark_test_data.csv"

di as txt _n _n
di as txt "================================================================================"
di as txt "              IMPORT DELIMITED vs CIMPORT SPEED BENCHMARK"
di as txt "================================================================================"
di as txt "Date: `c(current_date)' `c(current_time)'"
di as txt "Stata: `c(stata_version)' `c(processors)'-core"
di as txt "Observations: " %12.0fc `N'
di as txt "Variables: 10 (mixed types: double, int, string)"
di as txt "================================================================================"

********************************************************************************
* Generate test CSV file
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Phase 1: Generating test CSV file"
di as txt "{hline 80}"

* First create the data in Stata
set seed 54321
set obs `N'

* Create variables of varying types
* 4 doubles (continuous)
gen double var_dbl1 = rnormal()
gen double var_dbl2 = rnormal() * 1000
gen double var_dbl3 = runiform() * 10000
gen double var_dbl4 = rnormal() + runiform() * 100

* 3 integers
gen int var_int1 = ceil(runiform() * 1000)
gen long var_int2 = ceil(runiform() * 1000000)
gen byte var_int3 = ceil(runiform() * 100)

* 3 strings of varying lengths
gen str10 var_str1 = "str" + string(ceil(runiform() * 10000))
gen str20 var_str2 = "category" + string(ceil(runiform() * 100))
gen str5 var_str3 = substr("ABCDEFGHIJ", ceil(runiform() * 10), 1) + string(ceil(runiform() * 1000))

compress

di as txt "Data generated in memory. Exporting to CSV..."

* Export to CSV using Stata's export delimited
timer clear 99
timer on 99
quietly export delimited using "`testfile'", replace
timer off 99
quietly timer list 99
local export_time = r(t99)

* Get file size
local fsize = 0
tempname fh
file open `fh' using "`testfile'", read binary
file seek `fh' eof
local fsize = r(loc)
file close `fh'
local fsize_mb = `fsize' / 1048576

di as txt "CSV file created: `testfile'"
di as txt "  File size: " %12.1f `fsize_mb' " MB"
di as txt "  Export time: " %8.2f `export_time' " sec"

* Save copy of data for verification
tempfile original_data
quietly save `original_data'

********************************************************************************
* Benchmark 1: Stata's import delimited
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Benchmark 1: Stata's import delimited"
di as txt "{hline 80}"

clear
timer clear 1

di as txt "Running import delimited..."
timer on 1
quietly import delimited using "`testfile'", clear
timer off 1

quietly timer list 1
local t_import = r(t1)
local mbps_import = `fsize_mb' / `t_import'

di as txt "  Time:  " %10.2f `t_import' " sec"
di as txt "  Speed: " %10.1f `mbps_import' " MB/s"

* Store imported data for comparison
tempfile stata_imported
quietly save `stata_imported'

* Store variable info
describe, short
local stata_nobs = r(N)
local stata_nvars = r(k)

* Store checksums for verification
local stata_checksums ""
foreach var of varlist * {
    capture confirm numeric variable `var'
    if !_rc {
        quietly sum `var'
        local stata_checksums "`stata_checksums' `var'=`r(sum)'"
    }
}

********************************************************************************
* Benchmark 2: cimport delimited
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Benchmark 2: cimport delimited"
di as txt "{hline 80}"

clear
timer clear 2

di as txt "Running cimport delimited..."
timer on 2
quietly cimport delimited using "`testfile'", clear
timer off 2

quietly timer list 2
local t_cimport = r(t2)
local mbps_cimport = `fsize_mb' / `t_cimport'

di as txt "  Time:  " %10.2f `t_cimport' " sec"
di as txt "  Speed: " %10.1f `mbps_cimport' " MB/s"

* Store cimport data for comparison
tempfile cimport_data
quietly save `cimport_data'

* Store variable info
describe, short
local cimport_nobs = r(N)
local cimport_nvars = r(k)

* Store checksums for verification
local cimport_checksums ""
foreach var of varlist * {
    capture confirm numeric variable `var'
    if !_rc {
        quietly sum `var'
        local cimport_checksums "`cimport_checksums' `var'=`r(sum)'"
    }
}

********************************************************************************
* Verification: Compare imported datasets
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Verification: Comparing imported datasets"
di as txt "{hline 80}"

di as txt _n "Basic comparison:"
di as txt "  Observations: import=" %12.0fc `stata_nobs' "  cimport=" %12.0fc `cimport_nobs' cond(`stata_nobs'==`cimport_nobs', "  [MATCH]", "  [DIFFER]")
di as txt "  Variables:    import=" %12.0fc `stata_nvars' "  cimport=" %12.0fc `cimport_nvars' cond(`stata_nvars'==`cimport_nvars', "  [MATCH]", "  [DIFFER]")

* Full data comparison using cf
di as txt _n "Full data comparison:"
use `stata_imported', clear
capture cf _all using `cimport_data'
if _rc == 0 {
    di as txt "  Result: {result:PASS} - Datasets are identical"
    local verification = "PASS"
}
else {
    di as txt "  Result: {error:DIFFER} - Datasets have differences"
    local verification = "DIFFER"

    * Show which variables differ
    di as txt _n "  Checking individual variables:"
    foreach var of varlist * {
        capture cf `var' using `cimport_data'
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

di as txt _n _n
di as txt "================================================================================"
di as txt "                           SPEED BENCHMARK SUMMARY"
di as txt "================================================================================"
di as txt _n

local speedup = `t_import' / `t_cimport'

di as txt "{hline 60}"
di as txt %25s "Metric" " | " %15s "import" " | " %15s "cimport"
di as txt "{hline 60}"
di as txt %25s "Time (seconds)" " | " as result %13.2f `t_import' as txt " | " as result %13.2f `t_cimport'
di as txt %25s "Speed (MB/s)" " | " as result %13.1f `mbps_import' as txt " | " as result %13.1f `mbps_cimport'
di as txt %25s "Rows imported" " | " as result %13.0fc `stata_nobs' as txt " | " as result %13.0fc `cimport_nobs'
di as txt "{hline 60}"
di as txt _n

if `t_cimport' < `t_import' {
    di as txt "cimport is " as result %5.2f `speedup' as txt "x faster than import delimited"
}
else {
    local slowdown = `t_cimport' / `t_import'
    di as txt "cimport is " as result %5.2f `slowdown' as txt "x slower than import delimited"
}

di as txt _n "Data verification: `verification'"

di as txt _n "================================================================================"
di as txt "Test file: `testfile' (" %0.1f `fsize_mb' " MB)"
di as txt "================================================================================"

* Clean up test file
capture erase "`testfile'"
di as txt _n "Test file cleaned up."
