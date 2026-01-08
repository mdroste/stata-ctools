* Final comparison of qreg vs cqreg
clear all
set more off
adopath + "../build"

set seed 12345
set obs 5000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di "{hline 78}"
di "QREG"
di "{hline 78}"
qreg y x1 x2

di _newline(2) "{hline 78}"
di "CQREG"
di "{hline 78}"
cqreg y x1 x2

di _newline(2) "{hline 78}"
di "COMPARISON OF KEY e() VALUES"
di "{hline 78}"
quietly qreg y x1 x2
local qreg_r2 = e(r2_p)
local qreg_adev = e(sum_adev)
local qreg_rdev = e(sum_rdev)
local qreg_N = e(N)
local qreg_sparsity = e(sparsity)
matrix qreg_b = e(b)

quietly cqreg y x1 x2
local cqreg_r2 = e(r2_p)
local cqreg_adev = e(sum_adev)
local cqreg_rdev = e(sum_rdev)
local cqreg_N = e(N)
local cqreg_sparsity = e(sparsity)
matrix cqreg_b = e(b)

di as text "                        qreg          cqreg         diff"
di as text "{hline 60}"
di as text "N                " as result %12.0f `qreg_N' "  " %12.0f `cqreg_N' "  " %12.0f `qreg_N'-`cqreg_N'
di as text "Pseudo R2        " as result %12.6f `qreg_r2' "  " %12.6f `cqreg_r2' "  " %12.6f `qreg_r2'-`cqreg_r2'
di as text "Sum adev         " as result %12.4f `qreg_adev' "  " %12.4f `cqreg_adev' "  " %12.4f `qreg_adev'-`cqreg_adev'
di as text "Sum rdev         " as result %12.4f `qreg_rdev' "  " %12.4f `cqreg_rdev' "  " %12.4f `qreg_rdev'-`cqreg_rdev'
di as text "Sparsity         " as result %12.6f `qreg_sparsity' "  " %12.6f `cqreg_sparsity' "  " %12.6f `qreg_sparsity'-`cqreg_sparsity'
di as text "{hline 60}"
di as text "Coefficient diffs:"
di as text "  x1:   " as result %12.6f qreg_b[1,1]-cqreg_b[1,1]
di as text "  x2:   " as result %12.6f qreg_b[1,2]-cqreg_b[1,2]
di as text "  cons: " as result %12.6f qreg_b[1,3]-cqreg_b[1,3]
