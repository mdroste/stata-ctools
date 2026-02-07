/*******************************************************************************
 * benchmark_speed.do
 *
 * Speed benchmarks comparing ctools commands vs Stata native equivalents
 * Saves timestamped log to validation/ folder
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


quietly {

noi di as text "Benchmark started: `c(current_date)' `c(current_time)'"
noi di as text "Log file: `logfile'"

local N = 10000000
local G1 = 10000       // small groups (for FEs)
local G2 = 10000000    // large groups (for merge keys)

noi di as text _n "Creating test data with N = `N' observations..."

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

* Variables for encode/decode/destring benchmarks
gen str_cat = "category_" + string(ceil(runiform()*1000))
gen str_num = string(x1, "%9.3f")

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

noi di as text "Starting benchmarks..."

*===============================================================================
* SORT BENCHMARKS
*===============================================================================

noi di as text "  Running sort benchmarks"

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

noi di as text "  Running merge benchmarks"

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

noi di as text "  Running reghdfe benchmarks"

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

noi di as text "  Running qreg benchmarks"

* qreg vs cqreg
use `main', clear
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

noi di as text "  Running ivreghdfe benchmarks"

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

noi di as text "  Running binscatter benchmarks"

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

noi di as text "  Running import/export benchmarks"

* Export (1M rows subset)
use `main', clear
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

timer clear
timer on 1
cimport delimited using "`exp_stata'.csv", clear verbose
timer off 1
quietly timer list 1
local t_cimport = r(t1)

*===============================================================================
* EXCEL IMPORT/EXPORT BENCHMARKS
*===============================================================================

noi di as text "  Running Excel import/export benchmarks"

* xlsx has a 1,048,576 row limit per sheet
local N_excel = min(`N', 1000000)

* Export Excel - Stata
use `main' in 1/`N_excel', clear
tempfile exp_xlsx_stata exp_xlsx_ctools
timer clear
timer on 1
export excel using "`exp_xlsx_stata'.xlsx", firstrow(variables) replace
timer off 1
quietly timer list 1
local t_export_xl = r(t1)

* Export Excel - ctools
use `main' in 1/`N_excel', clear
timer clear
timer on 1
cexport excel using "`exp_xlsx_ctools'.xlsx", replace verbose
timer off 1
quietly timer list 1
local t_cexport_xl = r(t1)

* Import Excel - Stata
timer clear
timer on 1
import excel using "`exp_xlsx_stata'.xlsx", firstrow clear
timer off 1
quietly timer list 1
local t_import_xl = r(t1)

* Import Excel - ctools
timer clear
timer on 1
cimport excel using "`exp_xlsx_stata'.xlsx", firstrow clear verbose
timer off 1
quietly timer list 1
local t_cimport_xl = r(t1)

*===============================================================================
* DECODE BENCHMARKS
*===============================================================================

noi di as text "  Running decode benchmarks"

* First encode str_cat to get a labeled numeric variable
use `main', clear
encode str_cat, gen(cat_encoded)
tempfile encoded
save `encoded', replace

* Stata decode
timer clear
timer on 1
decode cat_encoded, gen(decoded_stata)
timer off 1
quietly timer list 1
local t_decode = r(t1)

* ctools cdecode
use `encoded', clear
timer clear
timer on 1
cdecode cat_encoded, gen(decoded_ctools) verbose
timer off 1
quietly timer list 1
local t_cdecode = r(t1)

*===============================================================================
* ENCODE BENCHMARKS
*===============================================================================

noi di as text "  Running encode benchmarks"

use `main', clear
timer clear
timer on 1
encode str_cat, gen(encoded_stata)
timer off 1
quietly timer list 1
local t_encode = r(t1)

use `main', clear
timer clear
timer on 1
cencode str_cat, gen(encoded_ctools) verbose
timer off 1
quietly timer list 1
local t_cencode = r(t1)

*===============================================================================
* DESTRING BENCHMARKS
*===============================================================================

noi di as text "  Running destring benchmarks"

use `main', clear
timer clear
timer on 1
destring str_num, gen(num_stata)
timer off 1
quietly timer list 1
local t_destring = r(t1)

use `main', clear
timer clear
timer on 1
cdestring str_num, gen(num_ctools) verbose
timer off 1
quietly timer list 1
local t_cdestring = r(t1)

*===============================================================================
* WINSOR BENCHMARKS
*===============================================================================

noi di as text "  Running winsor benchmarks"

use `main', clear
timer clear
timer on 1
winsor2 x1 x2 x3 x4, cuts(1 99) replace
timer off 1
quietly timer list 1
local t_winsor = r(t1)

use `main', clear
timer clear
timer on 1
cwinsor x1 x2 x3 x4, cuts(1 99) replace verbose
timer off 1
quietly timer list 1
local t_cwinsor = r(t1)

*===============================================================================
* SAMPLE BENCHMARKS
*===============================================================================

noi di as text "  Running sample benchmarks"

use `main', clear
set seed 12345
timer clear
timer on 1
sample 50
timer off 1
quietly timer list 1
local t_sample = r(t1)

use `main', clear
set seed 12345
timer clear
timer on 1
csample 50, verbose
timer off 1
quietly timer list 1
local t_csample = r(t1)

*===============================================================================
* BSAMPLE BENCHMARKS
*===============================================================================

noi di as text "  Running bsample benchmarks"

use `main', clear
set seed 12345
timer clear
timer on 1
bsample
timer off 1
quietly timer list 1
local t_bsample = r(t1)

use `main', clear
set seed 12345
timer clear
timer on 1
cbsample, verbose
timer off 1
quietly timer list 1
local t_cbsample = r(t1)

*===============================================================================
* RANGESTAT BENCHMARKS
*===============================================================================

noi di as text "  Running rangestat benchmarks"

capture which rangestat
local rangestat_installed = (_rc == 0)

if `rangestat_installed' {
    use `main', clear
    gen long obs = _n
    timer clear
    timer on 1
    rangestat (mean) rs_mean=x1, interval(obs -50 50)
    timer off 1
    quietly timer list 1
    local t_rangestat = r(t1)
}
else {
    noi di as text "    (rangestat not installed - skipping Stata comparison)"
    local t_rangestat = .
}

use `main', clear
gen long obs = _n
timer clear
timer on 1
crangestat (mean) cr_mean=x1, interval(obs -50 50) verbose
timer off 1
quietly timer list 1
local t_crangestat = r(t1)

*===============================================================================
* PSMATCH BENCHMARKS
*===============================================================================

noi di as text "  Running psmatch benchmarks"

capture which psmatch2
local psmatch2_installed = (_rc == 0)

local N_psm = 10000
clear
set obs `N_psm'
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()
tempfile psm_data
save `psm_data', replace

if `psmatch2_installed' {
    timer clear
    timer on 1
    psmatch2 treat x1 x2, outcome(y)
    timer off 1
    quietly timer list 1
    local t_psmatch = r(t1)
}
else {
    noi di as text "    (psmatch2 not installed - skipping Stata comparison)"
    local t_psmatch = .
}

use `psm_data', clear
timer clear
timer on 1
cpsmatch treat x1 x2, outcome(y) verbose
timer off 1
quietly timer list 1
local t_cpsmatch = r(t1)

*===============================================================================
* RESULTS SUMMARY
*===============================================================================

	
noi di as text _n "{hline 70}"
noi di as text "SPEED BENCHMARK RESULTS (seconds)"
noi di as text "{hline 70}"
noi di as text ""
noi di as text "Command                         Stata      ctools     Speedup"
noi di as text "{hline 70}"
noi di as text "SORT (N=`N'):"
noi di as text "  sort (double)             " %8.2f `t_sort_dbl' "    " %8.2f `t_csort_dbl' "    " %5.1f (`t_sort_dbl'/`t_csort_dbl') "x"
noi di as text "  sort (long)               " %8.2f `t_sort_long' "    " %8.2f `t_csort_long' "    " %5.1f (`t_sort_long'/`t_csort_long') "x"
noi di as text "  sort (string)             " %8.2f `t_sort_str' "    " %8.2f `t_csort_str' "    " %5.1f (`t_sort_str'/`t_csort_str') "x"
noi di as text "  sort (multi-key)          " %8.2f `t_sort_multi' "    " %8.2f `t_csort_multi' "    " %5.1f (`t_sort_multi'/`t_csort_multi') "x"
noi di as text ""
noi di as text "MERGE (N=`N', m:1):"
noi di as text "  merge m:1                 " %8.2f `t_merge' "    " %8.2f `t_cmerge' "    " %5.1f (`t_merge'/`t_cmerge') "x"
noi di as text ""
noi di as text "REGRESSION (N=`N'):"
noi di as text "  reghdfe                   " %8.2f `t_reghdfe' "    " %8.2f `t_creghdfe' "    " %5.1f (`t_reghdfe'/`t_creghdfe') "x"
noi di as text "  ivreghdfe                 " %8.2f `t_ivreghdfe' "    " %8.2f `t_civreghdfe' "    " %5.1f (`t_ivreghdfe'/`t_civreghdfe') "x"
noi di as text ""
noi di as text "QREG (N=`N'):"
noi di as text "  qreg                      " %8.2f `t_qreg' "    " %8.2f `t_cqreg' "    " %5.1f (`t_qreg'/`t_cqreg') "x"
noi di as text ""
noi di as text "BINSCATTER (N=`N'):"
noi di as text "  binscatter                " %8.2f `t_binscatter' "    " %8.2f `t_cbinscatter' "    " %5.1f (`t_binscatter'/`t_cbinscatter') "x"
noi di as text "  binscatter (controls+abs) " %8.2f `t_binscatter2' "    " %8.2f `t_cbinscatter2' "    " %5.1f (`t_binscatter2'/`t_cbinscatter2') "x"
noi di as text ""
noi di as text "DECODE (N=`N'):"
noi di as text "  decode                    " %8.2f `t_decode' "    " %8.2f `t_cdecode' "    " %5.1f (`t_decode'/`t_cdecode') "x"
noi di as text ""
noi di as text "ENCODE (N=`N'):"
noi di as text "  encode                    " %8.2f `t_encode' "    " %8.2f `t_cencode' "    " %5.1f (`t_encode'/`t_cencode') "x"
noi di as text ""
noi di as text "DESTRING (N=`N'):"
noi di as text "  destring                  " %8.2f `t_destring' "    " %8.2f `t_cdestring' "    " %5.1f (`t_destring'/`t_cdestring') "x"
noi di as text ""
noi di as text "WINSOR (N=`N'):"
noi di as text "  winsor2                   " %8.2f `t_winsor' "    " %8.2f `t_cwinsor' "    " %5.1f (`t_winsor'/`t_cwinsor') "x"
noi di as text ""
noi di as text "IMPORT/EXPORT DELIMITED (N=`N'):"
noi di as text "  export delimited          " %8.2f `t_export' "    " %8.2f `t_cexport' "    " %5.1f (`t_export'/`t_cexport') "x"
noi di as text "  import delimited          " %8.2f `t_import' "    " %8.2f `t_cimport' "    " %5.1f (`t_import'/`t_cimport') "x"
noi di as text ""
noi di as text "IMPORT/EXPORT EXCEL (N=`N_excel'):"
noi di as text "  export excel              " %8.2f `t_export_xl' "    " %8.2f `t_cexport_xl' "    " %5.1f (`t_export_xl'/`t_cexport_xl') "x"
noi di as text "  import excel              " %8.2f `t_import_xl' "    " %8.2f `t_cimport_xl' "    " %5.1f (`t_import_xl'/`t_cimport_xl') "x"
noi di as text ""
noi di as text "SAMPLE (N=`N'):"
noi di as text "  sample (50%)              " %8.2f `t_sample' "    " %8.2f `t_csample' "    " %5.1f (`t_sample'/`t_csample') "x"
noi di as text ""
noi di as text "BSAMPLE (N=`N'):"
noi di as text "  bsample                   " %8.2f `t_bsample' "    " %8.2f `t_cbsample' "    " %5.1f (`t_bsample'/`t_cbsample') "x"
noi di as text ""
noi di as text "RANGESTAT (N=`N'):"
if `t_rangestat' != . {
    noi di as text "  rangestat (mean)          " %8.2f `t_rangestat' "    " %8.2f `t_crangestat' "    " %5.1f (`t_rangestat'/`t_crangestat') "x"
}
else {
    noi di as text "  crangestat (mean)              N/A    " %8.2f `t_crangestat' "       N/A"
}
noi di as text ""
noi di as text "PSMATCH (N=`N_psm'):"
if `t_psmatch' != . {
    noi di as text "  psmatch2 (NN)             " %8.2f `t_psmatch' "    " %8.2f `t_cpsmatch' "    " %5.1f (`t_psmatch'/`t_cpsmatch') "x"
}
else {
    noi di as text "  cpsmatch (NN)                  N/A    " %8.2f `t_cpsmatch' "       N/A"
}
noi di as text "{hline 70}"

noi di as text _n "Benchmark completed: `c(current_date)' `c(current_time)'"
log close
noi di as text "Log saved to: `logfile'"

}
