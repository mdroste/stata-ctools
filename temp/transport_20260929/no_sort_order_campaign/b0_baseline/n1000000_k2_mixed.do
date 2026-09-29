clear
quietly set obs 1000000
local transport_pad "xxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen str8 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "baseline_n1000000_k2_mixed_r0" "8" "read" "0"
plugin call transport_io v1 v2, "baseline_n1000000_k2_mixed_r1" "8" "read" "0"
plugin call transport_io v1 v2, "baseline_n1000000_k2_mixed_r2" "8" "read" "0"
plugin call transport_io v1 v2, "baseline_n1000000_k2_mixed_r3" "8" "read" "0"
plugin call transport_io v1 v2, "baseline_n1000000_k2_mixed_r4" "8" "read" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
