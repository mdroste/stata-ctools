capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/xls-edits.log", text replace
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/helpers.do"
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
set obs 3
gen double x=_n
gen str8 s="new"
gen double day=mdy(6,17,2020)+_n
format day %td
foreach opt in "sheet(New)" "sheet(Sheet1)" "sheet(Sheet1, modify)" "sheet(Sheet1, replace)" "sheet(Other, modify)" "sheet(Other, replace)" "sheet(Sheet1, modify) keepcellfmt" "sheet(Sheet1) sheetmodify" "sheet(Sheet1) sheetreplace" "sheet(Sheet1, modify) replace" "keepcellfmt" "sheet(Sheet1, replace) keepcellfmt" "sheet(,modify)" {
    copy "`root'/styled.xls" "`root'/nedit.xls", replace
    copy "`root'/styled.xls" "`root'/cedit.xls", replace
    capture noisily export excel using "`root'/nedit.xls", `opt' cell(B2)
    local nr=_rc
    capture noisily cexport excel using "`root'/cedit.xls", `opt' cell(B2)
    local cr=_rc
    di "OPTION=`opt' NATIVE=`nr' C=`cr'"
    if `nr'!=`cr' cio_record "edit error `opt'" 9
    else if `nr' cio_record "edit error `opt'" 0
    else {
        preserve
        import excel using "`root'/nedit.xls", describe
        local k=r(N_worksheet)
        forvalues i=1/`k' {
            local sh`i' `"`r(worksheet_`i')'"'
        }
        forvalues i=1/`k' {
            import excel using "`root'/nedit.xls", sheet(`"`sh`i''"') clear
            local datalabel : data label
            local sortedby : sortedby
            quietly label dir
            local labelsets `r(names)'
            mata: cio_ref=cio_snapshot_data()
            import excel using "`root'/cedit.xls", sheet(`"`sh`i''"') clear
            local datalabel : data label
            local sortedby : sortedby
            quietly label dir
            local labelsets `r(names)'
            mata: cio_compare(cio_ref)
            local rc=cond(`same',0,9)
            cio_record "edit `opt' `sh`i''" `rc'
        }
        restore
    }
}
di "EDITS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
log close
exit, clear
