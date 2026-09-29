clear all
log using "temp/io_formats/blob.log", text replace
set obs 1
gen strL s=""
mata: st_sstore(1,"s",char((0,0,0,10,1,0,0,0)))
di "BLOB_LEN=" strlen(s[1])
di "BLOB_TYPE=" strltype(s[1])
log close
exit, clear
