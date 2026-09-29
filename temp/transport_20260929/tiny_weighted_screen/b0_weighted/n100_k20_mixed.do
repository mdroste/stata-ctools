clear
quietly set obs 100
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
quietly gen float v9 = (_n+9)/7
quietly replace v9 = . if mod(_n,101)==0
quietly replace v9 = .z if mod(_n,103)==0
quietly gen str32 v10 = cond(mod(_n,17)==0,"",substr(string(mod(_n+10,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen byte v11 = mod(_n+11,101)-50
quietly replace v11 = . if mod(_n,101)==0
quietly replace v11 = .z if mod(_n,103)==0
quietly gen str32 v12 = cond(mod(_n,17)==0,"",substr(string(mod(_n+12,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen long v13 = _n+13
quietly replace v13 = . if mod(_n,101)==0
quietly replace v13 = .z if mod(_n,103)==0
quietly gen str32 v14 = cond(mod(_n,17)==0,"",substr(string(mod(_n+14,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen double v15 = (_n+15)/7
quietly replace v15 = . if mod(_n,101)==0
quietly replace v15 = .z if mod(_n,103)==0
quietly gen str32 v16 = cond(mod(_n,17)==0,"",substr(string(mod(_n+16,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen int v17 = mod(_n+17,101)-50
quietly replace v17 = . if mod(_n,101)==0
quietly replace v17 = .z if mod(_n,103)==0
quietly gen str32 v18 = cond(mod(_n,17)==0,"",substr(string(mod(_n+18,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly gen float v19 = (_n+19)/7
quietly replace v19 = . if mod(_n,101)==0
quietly replace v19 = .z if mod(_n,103)==0
quietly gen str32 v20 = cond(mod(_n,17)==0,"",substr(string(mod(_n+20,1000000),"%06.0f")+"`transport_pad'",1,32))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r6" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r7" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r8" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r9" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r10" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r11" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_n100_k20_mixed_r12" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
