clear all
set more off
set linesize 255
capture log close _all
log using "temp/io_options/formula-probe.log", text replace
do "validation/io_test_helpers.do"
foreach location in C2 C3 C4 D2 D3 D4 A2 {
    foreach prefix in native result {
        clear
        set obs 1
        gen double value=99
        copy "temp/io_options/formulas/base.xlsx" "temp/io_options/formulas/`prefix'_`location'.xlsx", replace
        local command export
        if "`prefix'"=="result" local command cexport
        capture noisily `command' excel value using "temp/io_options/formulas/`prefix'_`location'.xlsx", sheet("Data",modify) cell(`location') keepcellfmt
        local rc=_rc
        di "FORMULA: `prefix' `location' rc=`rc'"
    }
    capture noisily import excel using "temp/io_options/formulas/native_`location'.xlsx", clear
    if !_rc {
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_ref=cio_snapshot_data()
        capture noisily import excel using "temp/io_options/formulas/result_`location'.xlsx", clear
        local rc=_rc
        if !`rc' {
            local datalabel : data label
            local sortedby : sortedby
            quietly label dir
            local labelsets `r(names)'
            mata: cio_compare(cio_ref)
            local rc=cond(`same',0,9)
        }
        cio_record "Formula edit `location' reimport" `rc'
    }
}
log close
exit, clear
