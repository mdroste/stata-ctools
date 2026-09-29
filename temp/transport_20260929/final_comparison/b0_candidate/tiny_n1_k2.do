clear
quietly set obs 1
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen int v2 = mod(_n+2,101)-50
quietly replace v2 = . if mod(_n,101)==0
quietly replace v2 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r0" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r1" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r2" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r3" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r4" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r5" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r6" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r7" "8" "identity" "0"
plugin call transport_io v1 v2, "candidate_tiny_n1_k2_r8" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
