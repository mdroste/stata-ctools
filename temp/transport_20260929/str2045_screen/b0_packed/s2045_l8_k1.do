clear
quietly set obs 200000
local transport_pad "xxxxxxxx"
quietly gen str2045 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1
plugin call transport_io v1, "packed_s2045_l8_k1_r0" "8" "identity" "0"
plugin call transport_io v1, "packed_s2045_l8_k1_r1" "8" "identity" "0"
plugin call transport_io v1, "packed_s2045_l8_k1_r2" "8" "identity" "0"
plugin call transport_io v1, "packed_s2045_l8_k1_r3" "8" "identity" "0"
plugin call transport_io v1, "packed_s2045_l8_k1_r4" "8" "identity" "0"
plugin call transport_io v1, "packed_s2045_l8_k1_r5" "8" "identity" "0"
plugin call transport_io v1, "packed_s2045_l8_k1_r6" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
