/* Differential native options, workbook mutation, and error atomicity. */
do "validation/io_test_helpers.do"
args root
if `"`root'"' == "" local root "`c(pwd)'/temp/io_excel_options"
capture mkdir "`root'"
set obs 4
gen double v=_n/7
gen double day=mdy(6,17,2020)+_n
format day %td
gen str6 text="str"
replace text="" in 2
foreach ext in xlsx xls {
    export excel using "`root'/numbers.`ext'", replace firstrow(variables)
}
foreach ext in xls xlsx {
    foreach opt in "allstring" "allstring(%9.2f)" "allstring(%21x)" "allstring(%td)" "allstring(%9.0gc)" "allstring(%9.2e)" "allstring(%9.0g) cellrange(:C4)" "cellrange(B2)" "cellrange(:B4)" "firstrow case(l)" "firstrow case(pre)" {
        cio_import_test using "`root'/numbers.`ext'", kind(excel) name("`ext' `opt'") opts(`opt')
    }
}
cio_import_test using "`root'/numbers", kind(excel) name("Excel extensionless prefers xls")
foreach ext in xlsx xls {
    clear
    set obs 3
    gen double x=_n
    replace x=. in 2
    gen str6 text="str"
    replace text="" in 2
    foreach missing in 99 NA {
        export excel using "`root'/nm.`ext'", missing(`missing') replace
        cexport excel using "`root'/cm.`ext'", missing(`missing') replace
        preserve
        import excel using "`root'/nm.`ext'", clear
        mata: cio_ref=cio_snapshot_data()
        import excel using "`root'/cm.`ext'", clear
        mata: cio_compare(cio_ref)
        local rc=cond(`same',0,9)
        cio_record "`ext' missing(`missing')" `rc'
        restore
    }
}

clear
set obs 4
gen double x=_n
gen str8 s="old"
gen double day=mdy(1,1,2020)+_n
format day %td
foreach ext in xlsx xls {
    export excel using "`root'/base.`ext'", firstrow(variables) sheet(Original) replace
    export excel using "`root'/base.`ext'", firstrow(variables) sheet(Untouched)
}
replace x=10+_n
replace s="new"
replace x=. in 2
replace s="" in 2
foreach ext in xlsx xls {
    foreach opt in "sheet(New)" "sheet(Original)" "sheet(Original, modify)" "sheet(Original, replace)" "sheet(Other, modify)" "sheet(Other, replace)" "sheet(Original, modify) keepcellfmt" "sheet(Original) sheetmodify" "sheet(Original) sheetreplace" "sheet(Original, modify) replace" "keepcellfmt" "sheet(Original, replace) keepcellfmt" "sheet(,modify)" "sheet(Original, modify) missing(99)" "sheet(Original, modify) missing(NA)" "sheet(Original, modify) firstrow(var)" {
        copy "`root'/base.`ext'" "`root'/native.`ext'", replace
        copy "`root'/base.`ext'" "`root'/result.`ext'", replace
        capture noisily export excel using "`root'/native.`ext'", `opt' cell(B2)
        local nr=_rc
        capture noisily cexport excel using "`root'/result.`ext'", `opt' cell(B2)
        local cr=_rc
        if `nr'!=`cr' cio_record "`ext' edit error `opt'" 9
        else if `nr' cio_record "`ext' edit error `opt'" 0
        else {
            preserve
            import excel using "`root'/native.`ext'", describe
            local k=r(N_worksheet)
            forvalues i=1/`k' {
                local sh`i' `"`r(worksheet_`i')'"'
            }
            forvalues i=1/`k' {
                import excel using "`root'/native.`ext'", sheet(`"`sh`i''"') clear
                local datalabel : data label
                local sortedby : sortedby
                quietly label dir
                local labelsets `r(names)'
                mata: cio_ref=cio_snapshot_data()
                import excel using "`root'/result.`ext'", sheet(`"`sh`i''"') clear
                local datalabel : data label
                local sortedby : sortedby
                quietly label dir
                local labelsets `r(names)'
                mata: cio_compare(cio_ref)
                local rc=cond(`same',0,9)
                cio_record "`ext' edit `opt' `sh`i''" `rc'
            }
            restore
        }
    }
}
clear
set obs 2
forvalues j=1/8 {
    gen v`j'=_n
}
label var v1 "A"
label var v2 "A"
label var v3 "123"
label var v4 "if"
label var v5 "A B"
label var v6 "C"
label var v7 "ΩABC"
label var v8 "ΩABC"
foreach ext in xls xlsx {
    export excel using "`root'/headers.`ext'", firstrow(varlabels) replace
}
foreach ext in xls xlsx {
    cio_import_test using "`root'/headers.`ext'", kind(excel) name("`ext' duplicate/numeric/reserved/Unicode names") opts(firstrow)
}
foreach opt in "allstring allstring(%9.2f)" "allstring(%9s)" "allstring(bad)" "case(p)" "cellrange(B2:)" {
    capture noisily import excel using "`root'/numbers.xlsx", `opt' clear
    local nr=_rc
    capture noisily cimport excel using "`root'/numbers.xlsx", `opt' clear
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "invalid Excel import option `opt'" `rc'
}
clear
set obs 3
gen byte x=_n
cio_export_test, kind(dbase) ext(dbf) name("DBF version(3)") opts(version(3))
foreach ver in 4 III IV iii iv {
    capture noisily export dbase using "`root'/native.dbf", version(`ver') replace
    local nr=_rc
    capture noisily cexport dbase using "`root'/result.dbf", version(`ver') replace
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "DBF rejects version(`ver')" `rc'
}
foreach ext in xls xlsx {
    foreach sel in "id=A note=C" "id=B note=D" "empty=D" "x=A y=A" "x y" "x y z w" {
        foreach opt in "" "allstring" "allstring(%9.2f)" {
            cio_import_test using "`root'/numbers.`ext'", kind(excel) name("`ext' selector `sel' `opt'") select(`sel') opts(`opt')
        }
    }
    foreach opt in "locale(UTF-8)" "locale(ISO-8859-1)" "locale(en_US)" "locale(tr_TR)" "locale(bogus)" {
        cio_import_test using "`root'/headers.`ext'", kind(excel) name("`ext' Unicode header `opt'") opts(firstrow `opt')
    }
    capture noisily import excel id=A using "`root'/numbers.`ext'", firstrow clear
    local nr=_rc
    capture noisily cimport excel id=A using "`root'/numbers.`ext'", firstrow clear
    local cr=_rc
    cio_record "`ext' extvarlist conflicts with firstrow" `=cond(`nr'==`cr',0,9)'
    foreach sel in "empty=Z" "id=3" "x=A x=B" {
        capture noisily import excel `sel' using "`root'/numbers.`ext'", clear
        local nr=_rc
        capture noisily cimport excel `sel' using "`root'/numbers.`ext'", clear
        local cr=_rc
        cio_record "`ext' invalid selector `sel'" `=cond(`nr'==`cr',0,9)'
    }
}
foreach ext in xls xlsx {
    foreach opt in "firstrow" "firstrow allstring" "firstrow allstring(%9.2f)" {
        cio_import_test using "`root'/display_formats.`ext'", kind(excel) name("`ext' built-in/custom date and time formats `opt'") opts(`opt')
    }
}
do "validation/validate_io_excel_formulas.do" "`root'"
foreach ext in xls xlsx {
    foreach opt in "clear" "firstrow" "allstring" "sheet(Formats)" "cellrange(A1)" "detail" {
        capture noisily import excel using "`root'/display_formats.`ext'", describe `opt'
        local nr=_rc
        capture noisily cimport excel using "`root'/display_formats.`ext'", describe `opt'
        cio_record "`ext' describe rejects `opt'" `=cond(`nr'==_rc,0,9)'
    }
    foreach corruption in invalid truncated {
        cio_import_state_test using "`root'/`ext'_`corruption'.`ext'", kind(excel) state(changed) opts(describe)
    }
}
di "IO_OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
if $IO_FORMATS_FAILED exit 9
