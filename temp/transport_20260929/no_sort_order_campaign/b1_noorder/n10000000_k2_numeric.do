clear
quietly set obs 10000000
local transport_pad "xxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen int v2 = mod(_n+2,101)-50
quietly replace v2 = . if mod(_n,101)==0
quietly replace v2 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "noorder_n10000000_k2_numeric_r0" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n10000000_k2_numeric_r1" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n10000000_k2_numeric_r2" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n10000000_k2_numeric_r3" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n10000000_k2_numeric_r4" "8" "read" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
