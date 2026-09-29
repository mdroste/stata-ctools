clear
quietly set obs 200000
local transport_pad "xxxxxxxx"
quietly gen str2045 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen str2045 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen str2045 v3 = cond(mod(_n,17)==0,"",substr(string(mod(_n+3,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen str2045 v4 = cond(mod(_n,17)==0,"",substr(string(mod(_n+4,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4, "candidate_maxstr_short_r6" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
