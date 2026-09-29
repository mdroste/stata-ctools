clear
quietly set obs 1000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen str244 v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v3 = cond(mod(_n,17)==0,"",substr(string(mod(_n+3,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v4 = cond(mod(_n,17)==0,"",substr(string(mod(_n+4,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v5 = cond(mod(_n,17)==0,"",substr(string(mod(_n+5,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v6 = cond(mod(_n,17)==0,"",substr(string(mod(_n+6,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v7 = cond(mod(_n,17)==0,"",substr(string(mod(_n+7,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v8 = cond(mod(_n,17)==0,"",substr(string(mod(_n+8,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v9 = cond(mod(_n,17)==0,"",substr(string(mod(_n+9,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v10 = cond(mod(_n,17)==0,"",substr(string(mod(_n+10,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v11 = cond(mod(_n,17)==0,"",substr(string(mod(_n+11,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v12 = cond(mod(_n,17)==0,"",substr(string(mod(_n+12,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v13 = cond(mod(_n,17)==0,"",substr(string(mod(_n+13,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v14 = cond(mod(_n,17)==0,"",substr(string(mod(_n+14,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v15 = cond(mod(_n,17)==0,"",substr(string(mod(_n+15,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v16 = cond(mod(_n,17)==0,"",substr(string(mod(_n+16,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v17 = cond(mod(_n,17)==0,"",substr(string(mod(_n+17,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v18 = cond(mod(_n,17)==0,"",substr(string(mod(_n+18,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v19 = cond(mod(_n,17)==0,"",substr(string(mod(_n+19,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly gen str244 v20 = cond(mod(_n,17)==0,"",substr(string(mod(_n+20,1000000),"%06.0f")+"`transport_pad'",1,244))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r6" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r7" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "tiny_n1000_k20_strings_r8" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
