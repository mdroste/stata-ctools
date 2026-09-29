clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/registry_transport/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/build"
program io0, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/registry_transport/current.plugin")
program run_benchmark
version 16
clear
quietly set obs 1
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_numeric_w8_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_numeric_w8_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_numeric_w8_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_numeric_w8_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_numeric_w8_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_numeric_w8_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_numeric_w8_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_numeric_w8_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_numeric_w8_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_numeric_w8_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_numeric_w8_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_numeric_w8_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_numeric_w8_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_numeric_w8_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_numeric_w8_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_numeric_w8_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_numeric_w8_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_numeric_w8_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_numeric_w8_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_numeric_w8_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_numeric_w8_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_numeric_w8_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_numeric_w8_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_numeric_w8_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
clear
quietly set obs 1
local transport_pad "xxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_mixed_w8_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w8_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w8_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w8_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w8_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w8_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w8_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w8_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w8_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w8_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w8_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w8_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_mixed_w8_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w8_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w8_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w8_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w8_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w8_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w8_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w8_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w8_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w8_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w8_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w8_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
clear
quietly set obs 1
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_mixed_w32_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w32_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w32_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w32_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w32_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w32_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w32_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w32_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w32_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w32_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w32_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w32_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_mixed_w32_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w32_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w32_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w32_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w32_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w32_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w32_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w32_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w32_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w32_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w32_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w32_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
clear
quietly set obs 1
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_mixed_w244_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w244_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w244_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w244_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w244_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w244_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w244_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w244_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w244_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w244_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w244_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w244_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_mixed_w244_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w244_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_mixed_w244_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w244_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w244_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_mixed_w244_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_mixed_w244_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w244_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_mixed_w244_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w244_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w244_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_mixed_w244_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
clear
quietly set obs 1
local transport_pad "xxxxxxxx"
quietly gen str8 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_strings_w8_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_strings_w8_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_strings_w8_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w8_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w8_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w8_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w8_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w8_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w8_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w8_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w8_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w8_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_strings_w8_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_strings_w8_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_strings_w8_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w8_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w8_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w8_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w8_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w8_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w8_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w8_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w8_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w8_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
clear
quietly set obs 1
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen str32 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_strings_w32_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_strings_w32_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_strings_w32_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w32_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w32_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w32_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w32_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w32_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w32_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w32_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w32_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w32_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_strings_w32_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_strings_w32_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_strings_w32_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w32_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w32_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w32_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w32_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w32_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w32_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w32_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w32_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w32_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
clear
quietly set obs 1
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen str244 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v*
plugin call io0 v*, "current_strings_w244_h1_r0" "8" "identity" "0"
plugin call io0 v*, "current_strings_w244_h1_r1" "8" "identity" "0"
plugin call io0 v*, "current_strings_w244_h1_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w244_h1_r0" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w244_h1_r1" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w244_h1_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w244_h1_r0" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w244_h1_r1" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w244_h1_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w244_h1_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w244_h1_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w244_h1_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
local __ctools_strw ""
plugin call io0 v*, "current_strings_w244_h0_r0" "8" "identity" "0"
plugin call io0 v*, "current_strings_w244_h0_r1" "8" "identity" "0"
plugin call io0 v*, "current_strings_w244_h0_r2" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w244_h0_r0" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w244_h0_r1" "8" "scattered" "0"
plugin call io0 v*, "current_strings_w244_h0_r2" "8" "scattered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v*, "current_strings_w244_h0_r0" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w244_h0_r1" "8" "sorted" "0"
plugin call io0 v*, "current_strings_w244_h0_r2" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w244_h0_r0" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w244_h0_r1" "8" "filtered" "0"
plugin call io0 v* if mod(_n,3)!=0, "current_strings_w244_h0_r2" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
end
capture noisily run_benchmark
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
