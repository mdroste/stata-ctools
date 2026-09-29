clear
quietly set obs 1000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen str244 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r0" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r1" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r2" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r3" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r4" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r5" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r6" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r7" "8" "identity" "0"
plugin call transport_io v1 v2, "baseline_n1000_k2_strings_r8" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
