clear all
set more off
set linesize 255
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/candidate"
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/c_external.log", text replace
program perf_clock, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/results/clock.plugin")
foreach case in numeric mixed {
    forvalues rep=0/5 {
        clear
        plugin call perf_clock, "start"
        quietly cimport excel using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/native_`case'.xlsx", firstrow clear
        plugin call perf_clock, "stop"
        di "EXTERNAL,c,`case',`rep'," %12.6f scalar(__ctools_perf_seconds)
        cf _all using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/expected_`case'.dta"
    }
}
di "EXTERNAL_COMPLETE"
log close
exit, clear
