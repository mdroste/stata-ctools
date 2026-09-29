capture log close _all
log using "temp/io_options/storage_probe.log", text replace
foreach width in 20 40 50 100 200 1000 {
 foreach pct in 1 50 70 75 76 77 78 79 80 85 90 95 96 {
  quietly import delimited using "temp/io_options/storagew`width'p`pct'.csv", clear
  local typ : type text
  di "STORAGE `width' `pct' `typ'"
 }
}
log close
exit, clear
