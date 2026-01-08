clear all
set more off
sysuse auto, clear

display "=== IID Residual method ==="
qreg price mpg, vce(iid, residual)
display "sparsity = " e(sparsity)
display "f_r = " e(f_r)
display "bwidth = " e(bwidth)
display "denmethod = " e(denmethod)

display ""
display "=== IID Fitted method (default) ==="
qreg price mpg
display "sparsity = " e(sparsity)
display "f_r = " e(f_r)
display "bwidth = " e(bwidth)
display "denmethod = " e(denmethod)

display ""
display "Note: sparsity = 1/f_r"
display "Residual: 1/f_r = " 1/0.0002241622
display "Fitted: 1/f_r = " 1/0.0001503742
