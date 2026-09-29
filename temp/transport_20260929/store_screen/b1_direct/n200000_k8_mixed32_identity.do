clear
quietly set obs 200000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen str32 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen long v3 = _n+3
quietly replace v3 = . if mod(_n,101)==0
quietly replace v3 = .z if mod(_n,103)==0
quietly gen str32 v4 = cond(mod(_n,17)==0,"",substr(string(mod(_n+4,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen double v5 = (_n+5)/7
quietly replace v5 = . if mod(_n,101)==0
quietly replace v5 = .z if mod(_n,103)==0
quietly gen str32 v6 = cond(mod(_n,17)==0,"",substr(string(mod(_n+6,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen int v7 = mod(_n+7,101)-50
quietly replace v7 = . if mod(_n,101)==0
quietly replace v7 = .z if mod(_n,103)==0
quietly gen str32 v8 = cond(mod(_n,17)==0,"",substr(string(mod(_n+8,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8, "direct_n200000_k8_mixed32_identity_r6" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
