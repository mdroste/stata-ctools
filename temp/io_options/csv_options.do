capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/csv_options.log", text replace
do "validation/io_test_helpers.do"
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
foreach opt in `"delimiters("|;",collapse)"' `"delimiters("||",asstring)"' `"delimiters("|;")"' `"delimiters("|;") collapsedelimiters"' {
    cio_import_test using "`root'/multidelim.txt", kind(delimited) name("delimiter case") opts(varnames(nonames) `macval(opt)')
}
cio_import_test using "`root'/long.csv", kind(delimited) name("CSV long Unicode strL") opts(encoding(utf8))
cio_import_test using "`root'/long.csv", kind(delimited) name("CSV long Unicode favorstrfixed") opts(encoding(utf8) favorstrfixed)
foreach names in "one" "one two" {
    cio_import_test using "`root'/multiline.csv", kind(delimited) name("CSV explicit names `names'") select(`names')
}
foreach limit in 1 2 3 20 unlimited {
    capture noisily import delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "maxquotedrows(`limit') error" `rc'
}
di "CSV_OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
log close
exit, clear
