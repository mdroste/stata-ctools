********************************************************************************
* benchmark_cexport.do
* Speed benchmark comparing export delimited and cexport performance
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
local testfile_stata "benchmark_export_stata.csv"
local testfile_ctools "benchmark_export_ctools.csv"

di as txt _n _n
di as txt "================================================================================"
di as txt "              EXPORT DELIMITED vs CEXPORT SPEED BENCHMARK"
di as txt "================================================================================"
di as txt "Date: `c(current_date)' `c(current_time)'"
di as txt "Stata: `c(stata_version)' `c(processors)'-core"
di as txt "Observations: " %12.0fc `N'
di as txt "Variables: 10 (mixed types: double, int, string)"
di as txt "================================================================================"

********************************************************************************
* Generate test data
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Phase 1: Generating test data in memory"
di as txt "{hline 80}"

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

describe, short
local data_nobs = r(N)
local data_nvars = r(k)

di as txt "Data generated:"
di as txt "  Observations: " %12.0fc `data_nobs'
di as txt "  Variables:    " %12.0fc `data_nvars'

* Save copy for reuse between benchmarks
tempfile testdata
quietly save `testdata'

********************************************************************************
* Benchmark 1: Stata's export delimited
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Benchmark 1: Stata's export delimited"
di as txt "{hline 80}"

timer clear 1

di as txt "Running export delimited..."
timer on 1
quietly export delimited using "`testfile_stata'", replace
timer off 1

quietly timer list 1
local t_export = r(t1)

* Get file size
local fsize_stata = 0
tempname fh
file open `fh' using "`testfile_stata'", read binary
file seek `fh' eof
local fsize_stata = r(loc)
file close `fh'
local fsize_stata_mb = `fsize_stata' / 1048576

local mbps_export = `fsize_stata_mb' / `t_export'
local rows_per_sec_export = `data_nobs' / `t_export'

di as txt "  Time:       " %10.2f `t_export' " sec"
di as txt "  File size:  " %10.1f `fsize_stata_mb' " MB"
di as txt "  Speed:      " %10.1f `mbps_export' " MB/s"
di as txt "  Throughput: " %10.0fc `rows_per_sec_export' " rows/sec"

********************************************************************************
* Benchmark 2: cexport delimited
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Benchmark 2: cexport delimited"
di as txt "{hline 80}"

* Reload data (in case export modified anything)
use `testdata', clear

timer clear 2

di as txt "Running cexport delimited..."
timer on 2
quietly cexport delimited using "`testfile_ctools'", replace
timer off 2

quietly timer list 2
local t_cexport = r(t2)

* Get file size
local fsize_ctools = 0
file open `fh' using "`testfile_ctools'", read binary
file seek `fh' eof
local fsize_ctools = r(loc)
file close `fh'
local fsize_ctools_mb = `fsize_ctools' / 1048576

local mbps_cexport = `fsize_ctools_mb' / `t_cexport'
local rows_per_sec_cexport = `data_nobs' / `t_cexport'

di as txt "  Time:       " %10.2f `t_cexport' " sec"
di as txt "  File size:  " %10.1f `fsize_ctools_mb' " MB"
di as txt "  Speed:      " %10.1f `mbps_cexport' " MB/s"
di as txt "  Throughput: " %10.0fc `rows_per_sec_cexport' " rows/sec"

********************************************************************************
* Verification: Compare exported files
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Verification: Comparing exported CSV files"
di as txt "{hline 80}"

di as txt _n "File size comparison:"
di as txt "  export delimited: " %12.0fc `fsize_stata' " bytes (" %0.1f `fsize_stata_mb' " MB)"
di as txt "  cexport delimited: " %12.0fc `fsize_ctools' " bytes (" %0.1f `fsize_ctools_mb' " MB)"
if `fsize_stata' == `fsize_ctools' {
    di as txt "  File sizes: [MATCH]"
}
else {
    local pct_diff = abs(`fsize_stata' - `fsize_ctools') / `fsize_stata' * 100
    di as txt "  File sizes: [DIFFER by " %0.2f `pct_diff' "%]"
}

* Import both files and compare the data
di as txt _n "Content verification (re-importing and comparing):"

* Import Stata's export
clear
quietly import delimited using "`testfile_stata'", clear
tempfile reimport_stata
quietly save `reimport_stata'
local stata_reimport_nobs = _N

* Import ctools' export
clear
quietly import delimited using "`testfile_ctools'", clear
tempfile reimport_ctools
quietly save `reimport_ctools'
local ctools_reimport_nobs = _N

di as txt "  Rows in Stata CSV:  " %12.0fc `stata_reimport_nobs'
di as txt "  Rows in ctools CSV: " %12.0fc `ctools_reimport_nobs'

* Compare the re-imported datasets
use `reimport_stata', clear
capture cf _all using `reimport_ctools'
if _rc == 0 {
    di as txt "  Data comparison: {result:PASS} - CSV contents are identical"
    local verification = "PASS"
}
else {
    di as txt "  Data comparison: {error:DIFFER} - CSV contents have differences"
    local verification = "DIFFER"

    * Show which variables differ
    di as txt _n "  Checking individual variables:"
    foreach var of varlist * {
        capture cf `var' using `reimport_ctools'
        if _rc == 0 {
            di as txt "    `var': [MATCH]"
        }
        else {
            di as txt "    `var': [DIFFER]"
        }
    }
}

* Also compare against original data
di as txt _n "Verification against original data:"
use `testdata', clear
capture cf _all using `reimport_stata'
if _rc == 0 {
    di as txt "  Stata export -> import: [MATCH]"
}
else {
    di as txt "  Stata export -> import: [DIFFER]"
}

use `testdata', clear
capture cf _all using `reimport_ctools'
if _rc == 0 {
    di as txt "  ctools export -> import: [MATCH]"
}
else {
    di as txt "  ctools export -> import: [DIFFER]"
}

********************************************************************************
* Summary
********************************************************************************

di as txt _n _n
di as txt "================================================================================"
di as txt "                           SPEED BENCHMARK SUMMARY"
di as txt "================================================================================"
di as txt _n

local speedup = `t_export' / `t_cexport'

di as txt "{hline 70}"
di as txt %25s "Metric" " | " %20s "export delimited" " | " %20s "cexport"
di as txt "{hline 70}"
di as txt %25s "Time (seconds)" " | " as result %18.2f `t_export' as txt " | " as result %18.2f `t_cexport'
di as txt %25s "File size (MB)" " | " as result %18.1f `fsize_stata_mb' as txt " | " as result %18.1f `fsize_ctools_mb'
di as txt %25s "Speed (MB/s)" " | " as result %18.1f `mbps_export' as txt " | " as result %18.1f `mbps_cexport'
di as txt %25s "Throughput (rows/sec)" " | " as result %18.0fc `rows_per_sec_export' as txt " | " as result %18.0fc `rows_per_sec_cexport'
di as txt "{hline 70}"
di as txt _n

if `t_cexport' < `t_export' {
    di as txt "cexport is " as result %5.2f `speedup' as txt "x faster than export delimited"
}
else {
    local slowdown = `t_cexport' / `t_export'
    di as txt "cexport is " as result %5.2f `slowdown' as txt "x slower than export delimited"
}

di as txt _n "Data verification: `verification'"

di as txt _n "================================================================================"
di as txt "Test files:"
di as txt "  Stata:  `testfile_stata' (" %0.1f `fsize_stata_mb' " MB)"
di as txt "  ctools: `testfile_ctools' (" %0.1f `fsize_ctools_mb' " MB)"
di as txt "================================================================================"

* Clean up test files
capture erase "`testfile_stata'"
capture erase "`testfile_ctools'"
di as txt _n "Test files cleaned up."
