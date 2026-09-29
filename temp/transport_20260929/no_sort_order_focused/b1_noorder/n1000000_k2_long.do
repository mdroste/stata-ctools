clear
quietly set obs 1000000
local transport_pad "xxxxxxxx"
quietly gen strL v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "noorder_n1000000_k2_long_r0" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n1000000_k2_long_r1" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n1000000_k2_long_r2" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n1000000_k2_long_r3" "8" "read" "0"
plugin call transport_io v1 v2, "noorder_n1000000_k2_long_r4" "8" "read" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
