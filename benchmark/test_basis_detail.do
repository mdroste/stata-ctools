clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid

* Show residuals near 0 with various thresholds
display "=== RESIDUALS NEAR ZERO ==="
sort resid
gen abs_resid = abs(resid)
format resid %15.10f
format abs_resid %15.10f

list resid abs_resid if abs_resid < 1e-6

* Count at different thresholds
foreach thresh in 1e-10 1e-9 1e-8 1e-7 1e-6 {
    count if abs_resid < `thresh'
    display "Threshold `thresh': " r(N) " basis obs"
}

* The actual basis observations in qreg are those with residual = 0 exactly
* In linear programming, the basis is where the constraint is exactly binding
* For qreg: u_i - v_i = residual, and either u_i = 0 or v_i = 0
* Basis observation: residual = 0 (on the hyperplane)

* Show the specific residual values around the median
_pctile resid if abs_resid >= 1e-8, percentile(25.958888 74.041112)
display _n "Percentile at 25.96: " %15.6f r(r1)
display "Percentile at 74.04: " %15.6f r(r2)

* Compare with my indexing
gen seq = _n if abs_resid >= 1e-8
sort resid
* With N_adj = 66, index = (66-1) * q = 65 * q
* q_lo = 0.25959: index = 16.87
* q_hi = 0.74041: index = 48.13

display _n "Manual check of order statistics:"
list resid if seq == 17
list resid if seq == 18
list resid if seq == 48
list resid if seq == 49
