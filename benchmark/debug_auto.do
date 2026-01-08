* Debug cqreg vs qreg on sysuse auto
clear all
set more off
adopath + "../build"

sysuse auto, clear

di _newline "{hline 78}"
di "QREG price mpg"
di "{hline 78}"
qreg price mpg

matrix qreg_b = e(b)
matrix qreg_V = e(V)
local qreg_sum_adev = e(sum_adev)
local qreg_sum_rdev = e(sum_rdev)
local qreg_r2_p = e(r2_p)
local qreg_q_v = e(q_v)

di _newline(2) "{hline 78}"
di "CQREG price mpg"
di "{hline 78}"
cqreg price mpg

matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
local cqreg_sum_adev = e(sum_adev)
local cqreg_sum_rdev = e(sum_rdev)
local cqreg_r2_p = e(r2_p)
local cqreg_q_v = e(q_v)

di _newline(2) "{hline 78}"
di "COEFFICIENT COMPARISON"
di "{hline 78}"
di "              qreg            cqreg           diff"
di "{hline 60}"
di "mpg:   " %14.8f qreg_b[1,1] "  " %14.8f cqreg_b[1,1] "  " %14.8f qreg_b[1,1]-cqreg_b[1,1]
di "_cons: " %14.8f qreg_b[1,2] "  " %14.8f cqreg_b[1,2] "  " %14.8f qreg_b[1,2]-cqreg_b[1,2]

di _newline "{hline 78}"
di "STANDARD ERROR COMPARISON"
di "{hline 78}"
di "              qreg            cqreg           diff"
di "{hline 60}"
di "mpg:   " %14.8f sqrt(qreg_V[1,1]) "  " %14.8f sqrt(cqreg_V[1,1]) "  " %14.8f sqrt(qreg_V[1,1])-sqrt(cqreg_V[1,1])
di "_cons: " %14.8f sqrt(qreg_V[2,2]) "  " %14.8f sqrt(cqreg_V[2,2]) "  " %14.8f sqrt(qreg_V[2,2])-sqrt(cqreg_V[2,2])

di _newline "{hline 78}"
di "SCALAR COMPARISON"
di "{hline 78}"
di "              qreg            cqreg           diff"
di "{hline 60}"
di "sum_adev: " %12.6f `qreg_sum_adev' "  " %12.6f `cqreg_sum_adev' "  " %12.6f `qreg_sum_adev'-`cqreg_sum_adev'
di "sum_rdev: " %12.6f `qreg_sum_rdev' "  " %12.6f `cqreg_sum_rdev' "  " %12.6f `qreg_sum_rdev'-`cqreg_sum_rdev'
di "r2_p:     " %12.8f `qreg_r2_p' "  " %12.8f `cqreg_r2_p' "  " %12.8f `qreg_r2_p'-`cqreg_r2_p'
di "q_v:      " %12.6f `qreg_q_v' "  " %12.6f `cqreg_q_v' "  " %12.6f `qreg_q_v'-`cqreg_q_v'
