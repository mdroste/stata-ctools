clear all
set more off
set linesize 255
capture log close _all
log using "temp/io_options/csv-results-probe.log", text replace
adopath ++ "build"
foreach enc in utf8 UTF-8 windows-1252 cp1252 latin1 ISO-8859-1 ascii US-ASCII UTF-16BE bogus {
    capture noisily import delimited using "temp/io_allformats/memory.csv", encoding(`enc') clear
    local rc=_rc
    di "ENC: `enc' rc=`rc' return=" `"`r(encoding)'"' " delimiters=" `"`r(delimiters)'"'
}
foreach pair in "long.csv" "utf32bom.csv" "utf32le.csv" "cp1251.csv" "sjis.csv" "latin2.csv" "range.csv" {
    capture noisily import delimited using "temp/io_delimited_options/`pair'", clear
    local rc=_rc
    di "AUTO: `pair' rc=`rc' return=" `"`r(encoding)'"' " delimiters=" `"`r(delimiters)'"'
}
foreach opt in `"delimiters("|;",collapse)"' `"delimiters("||",asstring)"' `"delimiters("|;")"' `"delimiters("tab")"' `"delimiters(" ")"' {
    capture noisily import delimited using "temp/io_delimited_options/multidelim.txt", `opt' clear
    local rc=_rc
    di "DELIM: `opt' rc=`rc' return=" `"`r(encoding)'"' " delimiters=" `"`r(delimiters)'"'
}
log close
exit, clear
