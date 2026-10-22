clear all
set more off
set linesize 255
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
capture log close _all
log using "/tmp/ctools-data-audit/bins_probes.log", text replace
capture which binscatter
di "BINS_REFERENCE_AVAILABLE=" (_rc==0)
set obs 3
gen x=_n
gen y=x*10
gen g=cond(_n==3,2,1)
frame create _ctools_graph_frame
frame _ctools_graph_frame: set obs 1
frame _ctools_graph_frame: gen protected=99
capture noisily cbinscatter y x, discrete by(g) linetype(none) nograph
di "CBINS_SINGLETON_RC=" _rc
capture confirm frame _ctools_graph_frame
di "USER_FRAME_SURVIVED=" (_rc==0)
matrix list e(bindata)
di "CBINS_REPORTED_N=" e(N)
capture noisily binscatter y x, discrete by(g) linetype(none) savedata("/tmp/ctools-data-audit/native_bins") replace
di "BINS_REFERENCE_RC=" _rc
capture noisily cbinscatter y x, nquantiles(2) nograph reportreg
di "CBINS_REPORTREG_RC=" _rc
program define twoway
    di `"CAPTURED_TWOWAY: `macval(0)'"'
end
capture noisily cbinscatter y x, nquantiles(2) linetype(connect)
di "CBINS_CONNECT_RC=" _rc
log close
exit, clear
