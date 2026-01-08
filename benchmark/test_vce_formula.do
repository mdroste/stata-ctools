clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
matrix V_qreg = e(V)
local N = e(N)
local df_r = e(df_r)
local q = 0.5

* Get residuals
predict double resid, resid

* Compute X'X^{-1}
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)

display "=== MODEL INFO ==="
display "N = " `N'
display "df_r = " `df_r'  // This is N - K
display "K = " (`N' - `df_r')

local K = `N' - `df_r'

* Check different VCE formulas
* Formula 1: V = s² * q(1-q) * (X'X)^{-1}  (no adjustment)
* Formula 2: V = s² * q(1-q) * (X'X)^{-1} * N/(N-K)  (df adjustment)
* Formula 3: V = s² * q(1-q) * (X'X/N)^{-1}  (per-obs X'X)

* From qreg's V[1,1], back-calculate sparsity for each formula
local v11 = V_qreg[1,1]
local xtxinv11 = XtX_inv[1,1]
local tau_factor = `q' * (1 - `q')

* Formula 1: s² = V / (τ(1-τ) * (X'X)^{-1})
local s2_f1 = `v11' / (`tau_factor' * `xtxinv11')
local s_f1 = sqrt(`s2_f1')

* Formula 2: s² = V / (τ(1-τ) * (X'X)^{-1} * N/(N-K))
local s2_f2 = `v11' / (`tau_factor' * `xtxinv11' * `N'/(`N'-`K'))
local s_f2 = sqrt(`s2_f2')

* Formula 3: s² = V * N / (τ(1-τ) * (X'X)^{-1})  [if V uses (X'X/N)^{-1}]
local s2_f3 = `v11' * `N' / (`tau_factor' * `xtxinv11' * `N')
* This is same as f1

display _n "=== BACK-CALCULATED SPARSITY ==="
display "Formula 1 (no adjustment): s = " %12.4f `s_f1'
display "Formula 2 (df adjustment): s = " %12.4f `s_f2'

* Compute DQ sparsity with Hall-Sheather
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

_pctile resid, percentile(`=(`q'-`h')*100' `=(`q'+`h')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local sparsity_dq = (`r_hi' - `r_lo') / (2 * `h')

display _n "=== DIFFERENCE QUOTIENT SPARSITY ==="
display "DQ sparsity (h_hs): " %12.4f `sparsity_dq'

display _n "=== RATIOS ==="
display "qreg_implied / DQ: " %8.6f (`s_f1' / `sparsity_dq')

* Check if Stata divides by df_r instead of N somewhere
local scale_factor = `s_f1' / `sparsity_dq'
display "Scale factor needed: " %8.6f `scale_factor'

* Let's compute what bandwidth would give the qreg sparsity
* sparsity = (r_hi - r_lo) / (2h)
* h = (r_hi - r_lo) / (2 * sparsity)
local h_implied = (`r_hi' - `r_lo') / (2 * `s_f1')
display _n "Implied h that would match qreg: " %12.8f `h_implied'
display "h_hs: " %12.8f `h'
display "h ratio: " %8.6f (`h' / `h_implied')

* Could Stata be using a narrower bandwidth?
* Let's check: if h is 1.47x smaller, what is it?
local h_smaller = `h' / `scale_factor'
display "If qreg uses h = " %12.8f `h_smaller'

* That would correspond to:
local q_lo_new = `q' - `h_smaller'
local q_hi_new = `q' + `h_smaller'
_pctile resid, percentile(`=`q_lo_new'*100' `=`q_hi_new'*100')
display "With this h: r_lo=" %8.2f r(r1) " r_hi=" %8.2f r(r2)
local sparsity_new = (r(r2) - r(r1)) / (2 * `h_smaller')
display "Resulting sparsity: " %12.4f `sparsity_new'
