clear
quietly set obs 500000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r0" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r1" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r2" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r3" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r4" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r5" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r6" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r7" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r8" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r9" "8" "filtered" "0"
plugin call transport_io v1 if mod(_n,3)!=0, "scalar_dense_k1_r10" "8" "filtered" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
