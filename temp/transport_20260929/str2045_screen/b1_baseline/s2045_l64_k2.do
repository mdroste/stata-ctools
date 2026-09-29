clear
quietly set obs 200000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen str2045 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,64))
quietly gen str2045 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,64))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r0" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r1" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r2" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r3" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r4" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r5" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_s2045_l64_k2_r6" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
