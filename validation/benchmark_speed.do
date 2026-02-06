/*******************************************************************************
 * benchmark_speed.do
 *
 * Speed benchmarks comparing ctools commands vs Stata native equivalents
 * Uses 25 million observations with typical use cases
 * Saves timestamped log to validation/ folder
 *
 * NOTE: Uses locals to store times because cimport clears Stata's timers
 ******************************************************************************/

clear all
set more off

* Add build directory with ctools to adopath
capture adopath ++ "build"
capture adopath ++ "../build"

* Create timestamped log
local datetime : di %tcCCYY-NN-DD-HH-MM-SS clock("`c(current_date)' `c(current_time)'", "DMY hms")
local datetime = subinstr("`datetime'", " ", "", .)
local logfile "validation/benchmark_speed_`datetime'.log"
capture log close _all
log using "`logfile'", replace text

di as text "Benchmark started: `c(current_date)' `c(current_time)'"
di as text "Log file: `logfile'"

local N = 25000000
local G1 = 10000       // small groups (for FEs)
local G2 = 10000000    // large groups (for merge keys)

di as text _n "Creating test data with N = `N' observations..."

*-------------------------------------------------------------------------------
* Create main test dataset (10 variables + outcome)
*-------------------------------------------------------------------------------
set obs `N'
forval i = 1/4 {
    gen double x`i' = rnormal()
    gen long g`i' = ceil(runiform()*`G1')
}
gen str6 s1 = "A1B2C3"
replace s1 = char(runiformint(65,90)) + char(runiformint(48,57)) + char(runiformint(65,90)) + char(runiformint(65,90)) + char(runiformint(48,57)) + char(runiformint(65,90))
gen str6 s2 = "A1B2C3"
replace s2 = char(runiformint(65,90)) + char(runiformint(48,57)) + char(runiformint(65,90)) + char(runiformint(65,90)) + char(runiformint(48,57)) + char(runiformint(65,90))
gen double y = 1 + x1 + 2*x2 + 3*x3 + sin(g1) + cos(g2) + rnormal()
compress
tempfile main
save `main', replace

* Merge dataset (m:1 on large int)
clear
set obs `G2'
gen long g1 = _n
gen double z1 = rnormal()
tempfile using1
save `using1', replace

di as text _n "Starting benchmarks..." _n

*===============================================================================
* SORT BENCHMARKS
*===============================================================================
di as text "{hline 60}"
di as text "SORT BENCHMARKS (csort vs sort)"
di as text "{hline 60}"

* Sort on double
use `main', clear
timer clear
timer on 1
sort x1
timer off 1
quietly timer list 1
local t_sort_dbl = r(t1)
use `main', clear
timer clear
timer on 1
csort x1, verbose
timer off 1
quietly timer list 1
local t_csort_dbl = r(t1)

* Sort on long int
use `main', clear
timer clear
timer on 1
sort g1
timer off 1
quietly timer list 1
local t_sort_long = r(t1)
use `main', clear
timer clear
timer on 1
csort g1, verbose
timer off 1
quietly timer list 1
local t_csort_long = r(t1)

* Sort on string
use `main', clear
timer clear
timer on 1
sort s1
timer off 1
quietly timer list 1
local t_sort_str = r(t1)
use `main', clear
timer clear
timer on 1
csort s1, verbose
timer off 1
quietly timer list 1
local t_csort_str = r(t1)

* Sort on multiple keys
use `main', clear
timer clear
timer on 1
sort g1 g2 x1
timer off 1
quietly timer list 1
local t_sort_multi = r(t1)
use `main', clear
timer clear
timer on 1
csort g1 g2 x1, verbose
timer off 1
quietly timer list 1
local t_csort_multi = r(t1)

*===============================================================================
* MERGE BENCHMARKS
*===============================================================================
di as text _n "{hline 60}"
di as text "MERGE BENCHMARKS (cmerge vs merge)"
di as text "{hline 60}"

* m:1 merge on large int
use `main', clear
timer clear
timer on 1
merge m:1 g1 using `using1', nogen
timer off 1
quietly timer list 1
local t_merge = r(t1)
use `main', clear
timer clear
timer on 1
cmerge m:1 g1 using `using1', nogen verbose
timer off 1
quietly timer list 1
local t_cmerge = r(t1)

*===============================================================================
* REGRESSION BENCHMARKS
*===============================================================================
di as text _n "{hline 60}"
di as text "REGRESSION BENCHMARKS"
di as text "{hline 60}"

* reghdfe vs creghdfe
use `main', clear
timer clear
timer on 1
reghdfe y x1 x2 x3, absorb(g1 g2)
timer off 1
quietly timer list 1
local t_reghdfe = r(t1)
timer clear
timer on 1
creghdfe y x1 x2 x3, absorb(g1 g2) verbose
timer off 1
quietly timer list 1
local t_creghdfe = r(t1)

* qreg vs cqreg (1M subset for speed)
use `main', clear
keep in 1/1000000
timer clear
timer on 1
qreg y x1 x2 x3
timer off 1
quietly timer list 1
local t_qreg = r(t1)
timer clear
timer on 1
cqreg y x1 x2 x3, verbose
timer off 1
quietly timer list 1
local t_cqreg = r(t1)

* ivreghdfe vs civreghdfe
use `main', clear
timer clear
timer on 1
ivreghdfe y (x1=x4) x2 x3, absorb(g1 g2)
timer off 1
quietly timer list 1
local t_ivreghdfe = r(t1)
timer clear
timer on 1
civreghdfe y (x1=x4) x2 x3, absorb(g1 g2) verbose
timer off 1
quietly timer list 1
local t_civreghdfe = r(t1)

*===============================================================================
* BINSCATTER BENCHMARKS
*===============================================================================
di as text _n "{hline 60}"
di as text "BINSCATTER BENCHMARKS"
di as text "{hline 60}"

use `main', clear
timer clear
timer on 1
binscatter y x1
timer off 1
quietly timer list 1
local t_binscatter = r(t1)
timer clear
timer on 1
cbinscatter y x1, nograph verbose
timer off 1
quietly timer list 1
local t_cbinscatter = r(t1)

* With controls and absorb
use `main', clear
timer clear
timer on 1
binscatter y x1, controls(x2 x3) absorb(g1)
timer off 1
quietly timer list 1
local t_binscatter2 = r(t1)
timer clear
timer on 1
cbinscatter y x1, controls(x2 x3) absorb(g1) nograph verbose
timer off 1
quietly timer list 1
local t_cbinscatter2 = r(t1)

*===============================================================================
* IMPORT/EXPORT BENCHMARKS
*===============================================================================
di as text _n "{hline 60}"
di as text "IMPORT/EXPORT BENCHMARKS"
di as text "{hline 60}"

* Export (1M rows subset)
use `main', clear
keep in 1/1000000
tempfile exp_stata exp_ctools
timer clear
timer on 1
export delimited using "`exp_stata'.csv", replace
timer off 1
quietly timer list 1
local t_export = r(t1)
timer clear
timer on 1
cexport delimited using "`exp_ctools'.csv", replace verbose
timer off 1
quietly timer list 1
local t_cexport = r(t1)

* Import
timer clear
timer on 1
import delimited using "`exp_stata'.csv", clear
timer off 1
quietly timer list 1
local t_import = r(t1)

* cimport - use wall clock since cimport clears Stata's timers
local t0 = clock("`c(current_date)' `c(current_time)'", "DMY hms")
cimport delimited using "`exp_stata'.csv", clear verbose
local t1 = clock("`c(current_date)' `c(current_time)'", "DMY hms")
local t_cimport = (`t1' - `t0') / 1000

*===============================================================================
* RESULTS SUMMARY
*===============================================================================
di as text _n "{hline 70}"
di as text "BENCHMARK RESULTS SUMMARY (seconds)"
di as text "{hline 70}"
di as text ""
di as text "Command                         Stata      ctools     Speedup"
di as text "{hline 70}"
di as text "SORT (N=25M):"
di as text "  sort (double)             " %8.2f `t_sort_dbl' "    " %8.2f `t_csort_dbl' "    " %5.1f (`t_sort_dbl'/`t_csort_dbl') "x"
di as text "  sort (long)               " %8.2f `t_sort_long' "    " %8.2f `t_csort_long' "    " %5.1f (`t_sort_long'/`t_csort_long') "x"
di as text "  sort (string)             " %8.2f `t_sort_str' "    " %8.2f `t_csort_str' "    " %5.1f (`t_sort_str'/`t_csort_str') "x"
di as text "  sort (multi-key)          " %8.2f `t_sort_multi' "    " %8.2f `t_csort_multi' "    " %5.1f (`t_sort_multi'/`t_csort_multi') "x"
di as text ""
di as text "MERGE (N=25M, m:1):"
di as text "  merge m:1                 " %8.2f `t_merge' "    " %8.2f `t_cmerge' "    " %5.1f (`t_merge'/`t_cmerge') "x"
di as text ""
di as text "REGRESSION (N=25M):"
di as text "  reghdfe                   " %8.2f `t_reghdfe' "    " %8.2f `t_creghdfe' "    " %5.1f (`t_reghdfe'/`t_creghdfe') "x"
di as text "  ivreghdfe                 " %8.2f `t_ivreghdfe' "    " %8.2f `t_civreghdfe' "    " %5.1f (`t_ivreghdfe'/`t_civreghdfe') "x"
di as text ""
di as text "QREG (N=1M):"
di as text "  qreg                      " %8.2f `t_qreg' "    " %8.2f `t_cqreg' "    " %5.1f (`t_qreg'/`t_cqreg') "x"
di as text ""
di as text "BINSCATTER (N=25M):"
di as text "  binscatter                " %8.2f `t_binscatter' "    " %8.2f `t_cbinscatter' "    " %5.1f (`t_binscatter'/`t_cbinscatter') "x"
di as text "  binscatter (controls+abs) " %8.2f `t_binscatter2' "    " %8.2f `t_cbinscatter2' "    " %5.1f (`t_binscatter2'/`t_cbinscatter2') "x"
di as text ""
di as text "IMPORT/EXPORT (N=1M):"
di as text "  export delimited          " %8.2f `t_export' "    " %8.2f `t_cexport' "    " %5.1f (`t_export'/`t_cexport') "x"
di as text "  import delimited          " %8.2f `t_import' "    " %8.2f `t_cimport' "    " %5.1f (`t_import'/`t_cimport') "x"
di as text "{hline 70}"

di as text _n "Benchmark completed: `c(current_date)' `c(current_time)'"
log close
di as text "Log saved to: `logfile'"
