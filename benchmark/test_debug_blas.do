clear all
set more off
adopath + "../build"
sysuse auto, clear

display "============================================================"
display "Debug: Testing cqreg with fitted method (Siddiqui)"
display "============================================================"

* Test 1: qreg default
display ""
display "=== qreg default ==="
quietly qreg price mpg
display "Coefficient: " _b[mpg]
display "SE: " _se[mpg]
display "Sparsity: " e(sparsity)

* Test 2: cqreg with fitted method (should use Siddiqui)
display ""
display "=== cqreg with fitted method ==="
cqreg price mpg, verbose
display "Coefficient: " _b[mpg]
display "SE: " _se[mpg]
display "Sparsity: " e(sparsity)

* Test 3: cqreg with residual method
display ""
display "=== cqreg with residual method ==="
cqreg price mpg, denmethod(residual) verbose
display "Coefficient: " _b[mpg]
display "SE: " _se[mpg]
display "Sparsity: " e(sparsity)

display ""
display "============================================================"
display "If fitted sparsity != 6650.08, the Siddiqui method is failing"
display "============================================================"
