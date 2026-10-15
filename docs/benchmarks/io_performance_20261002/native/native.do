clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/native.log", text replace
program perf_clock, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/results/clock.plugin")
foreach case in numeric mixed {
    use "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/results/expected_`case'_import_xlsx.dta", clear
    tempfile base
    quietly save `base'
    foreach op in export_csv import_csv export_xlsx import_xlsx {
        forvalues rep=0/2 {
            quietly use `base', clear
            plugin call perf_clock, "start"
            if "`op'"=="export_csv" quietly export delimited using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/native_`case'.csv", replace
            if "`op'"=="import_csv" quietly import delimited using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/results/`case'.csv", clear
            if "`op'"=="export_xlsx" quietly export excel using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/native/native_`case'.xlsx", replace firstrow(variables)
            if "`op'"=="import_xlsx" quietly import excel using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_performance/results/`case'.xlsx", clear firstrow
            plugin call perf_clock, "stop"
            di "NATIVE,`case',`op',`rep'," %12.6f scalar(__ctools_perf_seconds)
        }
    }
}
di "NATIVE_COMPLETE"
log close
exit, clear
