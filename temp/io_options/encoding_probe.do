capture log close _all
log using "temp/io_options/encoding_probe.log", text replace
foreach enc in utf32le utf32be utf16be cp1251 sjis latin2 {
 capture noisily import delimited using "temp/io_options/`enc'.csv", encoding(`enc') clear
 di "ENC `enc' RC=" _rc
 capture describe
 capture list
 capture return list
}
foreach opt in "encoding(bogus)" "parselocale(bogus)" "locale(bogus)" "locale(tr_TR)" "charset(cp1252) encoding(utf8)" {
 capture noisily import delimited using "temp/io_options/range.csv", `opt' clear
 di "OPTION `opt' RC=" _rc
}
foreach opt in "" "favorstrfixed" {
 import delimited using "temp/io_options/storage.csv", `opt' clear
 describe
}
log close
exit, clear
