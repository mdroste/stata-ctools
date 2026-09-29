clear all
set more off
capture log close
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/csv_extra_probe.log", text replace
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
foreach limit in 1 2 3 20 unlimited {
    capture noisily import delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local rc=_rc
    di "LIMIT=`limit' RC=`rc'"
    list
}
foreach names in "one" "one two" "one two three four" {
    capture noisily import delimited `names' using "`root'/multiline.csv", clear
    local rc=_rc
    di "NAMES=`names' RC=`rc'"
    describe
}
log close
exit, clear
