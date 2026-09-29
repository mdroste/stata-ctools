clear all
set more off
set linesize 255
capture log close _all
log using "temp/io_options/csv-results2-probe.log", text replace
foreach enc in latin1 Latin1 LATIN1 latin2 Latin2 latin5 latin9 unicode UTF16 UTF-16 UTF-32 utf-32le {
    capture noisily import delimited using "temp/io_allformats/memory.csv", encoding(`enc') clear
    local rc=_rc
    di "ENC: `enc' rc=`rc' return=" `"`r(encoding)'"'
}
local j=0
foreach opt in `"delimiters("|;",collapse)"' `"delimiters("||",asstring)"' `"delimiters("|;")"' `"delimiters("tab")"' `"delimiters(" ")"' {
    local ++j
    capture noisily import delimited using "temp/io_delimited_options/multidelim.txt", `opt' clear
    local rc=_rc
    local ret=r(delimiters)
    di `"DELIM: `j' rc=`rc' return=`macval(ret)'"'
}
foreach opt in "" "encoding(utf8)" "encoding(bogus)" {
    capture noisily import delimited using "temp/io_options/empty.csv", `opt' clear
    local rc=_rc
    di "EMPTY: rc=`rc' return=" `"`r(encoding)'"' " delimiters=" `"`r(delimiters)'"'
}
log close
exit, clear
