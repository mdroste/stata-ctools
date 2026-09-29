clear
quietly set obs 1000000
local transport_pad "xxxxxxxx"
quietly gen strL v1 = cond(mod(_n,17)==0,"",substr(string(mod(_n+1,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v3 = cond(mod(_n,17)==0,"",substr(string(mod(_n+3,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v4 = cond(mod(_n,17)==0,"",substr(string(mod(_n+4,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v5 = cond(mod(_n,17)==0,"",substr(string(mod(_n+5,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v6 = cond(mod(_n,17)==0,"",substr(string(mod(_n+6,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v7 = cond(mod(_n,17)==0,"",substr(string(mod(_n+7,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v8 = cond(mod(_n,17)==0,"",substr(string(mod(_n+8,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v9 = cond(mod(_n,17)==0,"",substr(string(mod(_n+9,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v10 = cond(mod(_n,17)==0,"",substr(string(mod(_n+10,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v11 = cond(mod(_n,17)==0,"",substr(string(mod(_n+11,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v12 = cond(mod(_n,17)==0,"",substr(string(mod(_n+12,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v13 = cond(mod(_n,17)==0,"",substr(string(mod(_n+13,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v14 = cond(mod(_n,17)==0,"",substr(string(mod(_n+14,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v15 = cond(mod(_n,17)==0,"",substr(string(mod(_n+15,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v16 = cond(mod(_n,17)==0,"",substr(string(mod(_n+16,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v17 = cond(mod(_n,17)==0,"",substr(string(mod(_n+17,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v18 = cond(mod(_n,17)==0,"",substr(string(mod(_n+18,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v19 = cond(mod(_n,17)==0,"",substr(string(mod(_n+19,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen strL v20 = cond(mod(_n,17)==0,"",substr(string(mod(_n+20,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "noorder_n1000000_k20_long_r0" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "noorder_n1000000_k20_long_r1" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "noorder_n1000000_k20_long_r2" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "noorder_n1000000_k20_long_r3" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "noorder_n1000000_k20_long_r4" "8" "read" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
