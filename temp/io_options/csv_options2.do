capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/csv_options2.log", text replace
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

foreach opt in "varnames(1) rowrange(3:4)" "varnames(1) rowrange(2:3)" "rowrange(3:4)" "varnames(1) rowrange(4)" "varnames(1) rowrange(1:1)" {
    cio_import_test using "`root'/range.csv", kind(delimited) name("CSV range inference `opt'") opts(`opt')
}
import delimited using "`root'/long.csv", encoding(utf8) clear
cio_export_test, kind(delimited) ext(csv) name("CSV long Unicode writer") opts(quote) iopts(encoding(utf8))
import delimited using "`root'/long.csv", encoding(utf8) clear
cio_export_test, kind(delimited) ext(csv) name("CSV long Unicode filtered writer") select(text if id==1) iopts(encoding(utf8))
di "CSV_OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
log close
exit, clear
