clear all
set more off
set linesize 255
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/baseline"
program perf_clock, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/results/clock.plugin")
foreach case in numeric mixed {
    forvalues rep=0/5 {
        clear
        plugin call perf_clock, "start"
        quietly cimport excel using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/native_`case'.xlsx", firstrow clear
        plugin call perf_clock, "stop"
        di "EXTERNAL,b,`case',`rep'," %12.6f scalar(__ctools_perf_seconds)
        if `rep'==0 quietly save "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/expected_`case'.dta", replace
    }
}
di "EXTERNAL_COMPLETE"
