capture log close _all
log using "temp/io_options/encoding_probe2.log", text replace
foreach enc in utf-32le UTF-32BE utf-16be utf16 windows-1251 cp1251 SHIFT-JIS shift_jis sjis ISO-8859-2 latin2 {
 local f cp1251
 if strpos(lower("`enc'"),"32le") local f utf32le
 if strpos(lower("`enc'"),"32be") local f utf32be
 if strpos(lower("`enc'"),"16") local f utf16be
 if strpos(lower("`enc'"),"jis") | "`enc'"=="sjis" local f sjis
 if strpos(lower("`enc'"),"8859") | "`enc'"=="latin2" local f latin2
 capture noisily import delimited using "temp/io_options/`f'.csv", encoding(`enc') clear
 local rc=_rc
 di "ENC `enc' RC=`rc'"
 if !`rc' {
  describe
  list
  return list
 }
}
foreach locale in de_DE de fr_FR fr en en_US tr_TR bogus {
 capture noisily import delimited using "temp/io_options/range.csv", parselocale(`locale') clear
 di "LOCALE `locale' RC=" _rc
}
log close
exit, clear
