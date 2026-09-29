clear
quietly set obs 100000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen strL v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen strL v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "candidate_strl_control_r0" "8" "read" "0"
plugin call transport_io v1 v2, "candidate_strl_control_r1" "8" "read" "0"
plugin call transport_io v1 v2, "candidate_strl_control_r2" "8" "read" "0"
plugin call transport_io v1 v2, "candidate_strl_control_r3" "8" "read" "0"
plugin call transport_io v1 v2, "candidate_strl_control_r4" "8" "read" "0"
plugin call transport_io v1 v2, "candidate_strl_control_r5" "8" "read" "0"
plugin call transport_io v1 v2, "candidate_strl_control_r6" "8" "read" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
