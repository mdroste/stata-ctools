/* Native/C workbook formula edits: the Python runner also checks XML structure. */
args root
if `"`root'"' == "" local root "`c(pwd)'/temp/io_excel_options"
forvalues test=1/15 {
    local kind simple
    local rows 1
    local cols 1
    local action modify
    local location : word `test' of C2 C3 C4 D2 D3 D4 A2 C2 D3 C2 B2 C2 C2 C2 C2
    if `test'>7 local kind grid
    if `test'==10 | `test'==11 local rows 2
    if `test'==10 local cols 2
    if `test'==11 local cols 3
    if `test'==15 local cols 3
    if `test'==12 {
        local rows 3
        local cols 3
    }
    if `test'==13 local action replace
    if `test'==14 local action new
    clear
    set obs `rows'
    forvalues j=1/`cols' {
        gen double value`j'=99+_n+`j'
    }
    foreach prefix in native result {
        copy "`root'/formulas_`kind'.xlsx" "`root'/formula_`prefix'_`test'.xlsx", replace
        local command export
        if "`prefix'"=="result" local command cexport
        local sheetopt sheet("Data",`action') keepcellfmt
        if "`action'"=="new" local sheetopt sheet("New")
        capture noisily `command' excel using "`root'/formula_`prefix'_`test'.xlsx", `sheetopt' cell(`location')
        local `prefix'rc=_rc
    }
    if `nativerc' | `resultrc' {
        cio_record "XLSX formula edit `test' export" 9
        continue
    }
    local sheets Data
    if "`action'"=="new" local sheets Data New
    foreach sheet of local sheets {
        import excel using "`root'/formula_native_`test'.xlsx", sheet("`sheet'") clear
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_ref=cio_snapshot_data()
        import excel using "`root'/formula_result_`test'.xlsx", sheet("`sheet'") clear
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_compare(cio_ref)
        cio_record "XLSX formula edit `test' `sheet' values/metadata" `=cond(`same',0,9)'
    }
}
