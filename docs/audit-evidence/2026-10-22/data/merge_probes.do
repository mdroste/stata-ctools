clear all
set more off
set linesize 255
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
capture log close _all
log using "/tmp/ctools-data-audit/merge_probes.log", text replace
set obs 3
gen id = _n
gen x1 = id*10
gen x2 = id*100
gen extra = id*1000
save "/tmp/ctools-data-audit/using.dta", replace
keep id
save "/tmp/ctools-data-audit/master.dta", replace
capture noisily cmerge 1:1 id using "/tmp/ctools-data-audit/using.dta", keepusing(x*) noreport
di "CMERGE_WILDCARD_RC=" _rc
use "/tmp/ctools-data-audit/master.dta", clear
capture noisily merge 1:1 id using "/tmp/ctools-data-audit/using.dta", keepusing(x*) noreport
di "MERGE_WILDCARD_RC=" _rc
use "/tmp/ctools-data-audit/master.dta", clear
keep if 0
capture noisily cmerge 1:1 _n using "/tmp/ctools-data-audit/using.dta", keepusing(x1) noreport
di "CMERGE_EMPTY_N_RC=" _rc
unab vars: _all
di "CMERGE_EMPTY_N_VARS=`vars'"
use "/tmp/ctools-data-audit/master.dta", clear
keep if 0
capture noisily merge 1:1 _n using "/tmp/ctools-data-audit/using.dta", keepusing(x1) noreport
di "MERGE_EMPTY_N_RC=" _rc
unab vars: _all
di "MERGE_EMPTY_N_VARS=`vars'"
use "/tmp/ctools-data-audit/master.dta", clear
keep if 0
capture noisily cmerge 1:1 id using "/tmp/ctools-data-audit/using.dta", keep(master) assert(master) noreport
di "CMERGE_EMPTY_ASSERT_RC=" _rc " N=" _N
use "/tmp/ctools-data-audit/master.dta", clear
keep if 0
capture noisily merge 1:1 id using "/tmp/ctools-data-audit/using.dta", keep(master) assert(master) noreport
di "MERGE_EMPTY_ASSERT_RC=" _rc " N=" _N
clear
set obs 6
gen x = mod(_n,2)
capture noisily cwinsor x, by(x) cuts(10 90)
di "WINSOR_REPLACE_OVERLAP_RC=" _rc
clear
set obs 2
gen s="other"
gen wt=-99
capture noisily cbsample if s != "weight(wt)"
di "BSAMPLE_QUOTED_IF_RC=" _rc
list
log close
exit, clear
