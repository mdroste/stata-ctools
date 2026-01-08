clear all
set more off
sysuse auto, clear

quietly qreg price mpg
predict resid, residuals

* Count LP basis observations (residual = 0)
count if abs(resid) < 1e-8
local n_basis = r(N)
display "LP basis observations (|r| < 1e-8): " `n_basis'

* Hall-Sheather bandwidth
local N = 74
local q = 0.5
local alpha = 0.05
local z_alpha = invnormal(1 - `alpha'/2)  
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h_prob = `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3) * `N'^(-1/3)
display "h_prob (Hall-Sheather) = " `h_prob'

* Residual quantiles for bandwidth
_pctile resid, p(`=(`q'-`h_prob')*100' `=(`q'+`h_prob')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local h_resid = (`r_hi' - `r_lo') / 2
display "h_resid = " `h_resid'

* Kernel density at 0 (all observations)
gen kernel_all = 0
replace kernel_all = 0.75 * (1 - (resid/`h_resid')^2) if abs(resid) < `h_resid'
qui sum kernel_all
local f0_all = r(sum) / (`N' * `h_resid')
local s_all = 1 / `f0_all'
display "Kernel density at 0 (all obs): f(0)=" `f0_all' ", sparsity=" `s_all'

* Kernel density at 0 (excluding LP basis)
gen kernel_no_basis = 0
replace kernel_no_basis = 0.75 * (1 - (resid/`h_resid')^2) if abs(resid) < `h_resid' & abs(resid) >= 1e-8
qui sum kernel_no_basis
local N_adj = `N' - `n_basis'
local f0_nb = r(sum) / (`N_adj' * `h_resid')
local s_nb = 1 / `f0_nb'
display "Kernel density at 0 (no basis): f(0)=" `f0_nb' ", sparsity=" `s_nb'

* What SE would these give?
* First get (X'X)^{-1}[1,1] from the IID residual result
quietly qreg price mpg, vce(iid, residual)
local se_iid_resid = _se[mpg]
local s_iid_resid = e(sparsity)
local xtx_inv11 = (`se_iid_resid'^2) / (0.25 * `s_iid_resid'^2)
display "(X'X)^{-1}[1,1] = " `xtx_inv11'

local se_all = sqrt(0.25 * `s_all'^2 * `xtx_inv11')
local se_nb = sqrt(0.25 * `s_nb'^2 * `xtx_inv11')
display "Predicted SE (all obs):      " `se_all'
display "Predicted SE (no basis):     " `se_nb'

display ""
display "Actual qreg robust SE:       " 43.452
display "Actual qreg IID residual SE: " 45.117
