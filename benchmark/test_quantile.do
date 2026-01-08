clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid

* Count basis observations with same threshold
gen byte basis = abs(resid) < 1e-8
count if basis
local n_basis = r(N)
local N_full = _N
local N_adj = _N - `n_basis'

display "N = " `N_full'
display "Basis obs = " `n_basis'
display "N_adj = " `N_adj'

* Compute bandwidth
local tau = 0.5
local alpha = 0.05
local z_alpha = invnormal(1 - `alpha'/2)
local z_tau = invnormal(`tau')
local phi_tau = normalden(`z_tau')
local h = `N_adj'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_tau'^2) / (2*`z_tau'^2 + 1))^(1/3)

display "h = " %12.8f `h'

* Compute residual quantiles at tau ± h, excluding basis
local q_lo = `tau' - `h'
local q_hi = `tau' + `h'
display "q_lo = " %12.8f `q_lo'
display "q_hi = " %12.8f `q_hi'

_pctile resid if !basis, percentile(`=`q_lo'*100' `=`q_hi'*100')
local r_lo = r(r1)
local r_hi = r(r2)
display "r_lo = " %12.6f `r_lo'
display "r_hi = " %12.6f `r_hi'

local sparsity = (`r_hi' - `r_lo') / (2 * `h')
display "Sparsity = " %12.4f `sparsity'

* List the non-basis residuals sorted
sort resid
list resid if !basis in 1/5
list resid if !basis in -5/l

* Also show the specific quantiles
local idx_lo = round(`N_adj' * `q_lo')
local idx_hi = round(`N_adj' * `q_hi')
display "Index for q_lo (0-based): " `idx_lo'
display "Index for q_hi (0-based): " `idx_hi'
