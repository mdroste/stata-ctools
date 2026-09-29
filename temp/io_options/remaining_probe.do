capture log close _all
log using "temp/io_options/remaining_probe.log", text replace
adopath ++ build
foreach pct in 1 5 10 25 50 80 99 {
 import delimited using "temp/io_options/storage`pct'.csv", clear
 local typ : type text
 di "STORAGE `pct' `typ'"
}
foreach opt in "" "member(SECOND)" {
 import sasxport5 "temp/io_allformats/native/multiple.xpt", describe `opt'
 return list
}
foreach ext in xls xlsx {
 foreach locale in bogus UTF-8 ISO-8859-1 en_US tr_TR {
  capture noisily import excel using "temp/io_excel_options/numbers.`ext'", locale(`locale') clear
  di "EXCEL IMPORT `ext' `locale' RC=" _rc
  capture noisily export excel using "temp/io_options/locale.`ext'", locale(`locale') replace
  di "EXCEL EXPORT `ext' `locale' RC=" _rc
 }
}
log close
exit, clear
