clear
quietly set obs 3275
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen int v2 = mod(_n+2,101)-50
quietly replace v2 = . if mod(_n,101)==0
quietly replace v2 = .z if mod(_n,103)==0
quietly gen long v3 = _n+3
quietly replace v3 = . if mod(_n,101)==0
quietly replace v3 = .z if mod(_n,103)==0
quietly gen float v4 = (_n+4)/7
quietly replace v4 = . if mod(_n,101)==0
quietly replace v4 = .z if mod(_n,103)==0
quietly gen double v5 = (_n+5)/7
quietly replace v5 = . if mod(_n,101)==0
quietly replace v5 = .z if mod(_n,103)==0
quietly gen byte v6 = mod(_n+6,101)-50
quietly replace v6 = . if mod(_n,101)==0
quietly replace v6 = .z if mod(_n,103)==0
quietly gen int v7 = mod(_n+7,101)-50
quietly replace v7 = . if mod(_n,101)==0
quietly replace v7 = .z if mod(_n,103)==0
quietly gen long v8 = _n+8
quietly replace v8 = . if mod(_n,101)==0
quietly replace v8 = .z if mod(_n,103)==0
quietly gen float v9 = (_n+9)/7
quietly replace v9 = . if mod(_n,101)==0
quietly replace v9 = .z if mod(_n,103)==0
quietly gen double v10 = (_n+10)/7
quietly replace v10 = . if mod(_n,101)==0
quietly replace v10 = .z if mod(_n,103)==0
quietly gen byte v11 = mod(_n+11,101)-50
quietly replace v11 = . if mod(_n,101)==0
quietly replace v11 = .z if mod(_n,103)==0
quietly gen int v12 = mod(_n+12,101)-50
quietly replace v12 = . if mod(_n,101)==0
quietly replace v12 = .z if mod(_n,103)==0
quietly gen long v13 = _n+13
quietly replace v13 = . if mod(_n,101)==0
quietly replace v13 = .z if mod(_n,103)==0
quietly gen float v14 = (_n+14)/7
quietly replace v14 = . if mod(_n,101)==0
quietly replace v14 = .z if mod(_n,103)==0
quietly gen double v15 = (_n+15)/7
quietly replace v15 = . if mod(_n,101)==0
quietly replace v15 = .z if mod(_n,103)==0
quietly gen byte v16 = mod(_n+16,101)-50
quietly replace v16 = . if mod(_n,101)==0
quietly replace v16 = .z if mod(_n,103)==0
quietly gen int v17 = mod(_n+17,101)-50
quietly replace v17 = . if mod(_n,101)==0
quietly replace v17 = .z if mod(_n,103)==0
quietly gen long v18 = _n+18
quietly replace v18 = . if mod(_n,101)==0
quietly replace v18 = .z if mod(_n,103)==0
quietly gen float v19 = (_n+19)/7
quietly replace v19 = . if mod(_n,101)==0
quietly replace v19 = .z if mod(_n,103)==0
quietly gen double v20 = (_n+20)/7
quietly replace v20 = . if mod(_n,101)==0
quietly replace v20 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r6" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r7" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r8" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r9" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r10" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r11" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20, "weighted_boundary_n3275_k20_numeric_r12" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
