capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/import_options.log", text replace
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/helpers.do"
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
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
di "OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
log close
exit, clear
