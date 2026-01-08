clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

* Run qreg to get basis observations
qreg price mpg weight, quantile(0.5)
predict double resid, resid

* Count basis observations (residual very close to 0)
gen byte basis = abs(resid) < 1e-8
count if basis
local n_basis = r(N)
local N_full = _N
local N_excl = _N - `n_basis'
display "Full N: " `N_full'
display "Basis obs: " `n_basis'
display "N after exclusion: " `N_excl'

* Get bandwidth for N vs N_excl
local tau = 0.5
local alpha = 0.05
local z_alpha = invnormal(1 - `alpha'/2)
local z_tau = invnormal(`tau')
local phi_tau = normalden(`z_tau')

* Hall-Sheather bandwidth for full N
local h_full = `N_full'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_tau'^2) / (2*`z_tau'^2 + 1))^(1/3)

* Hall-Sheather bandwidth for N_excl
local h_excl = `N_excl'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_tau'^2) / (2*`z_tau'^2 + 1))^(1/3)

display _n "=== BANDWIDTH COMPARISON ==="
display "h (full N=" `N_full' "): " %12.8f `h_full'
display "h (excl N=" `N_excl' "): " %12.8f `h_excl'
display "Stata reported h: 0.24041112"

* Compute sparsity with exclusion
_pctile resid if !basis, percentile(`=(`tau'-`h_excl')*100' `=(`tau'+`h_excl')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local sparsity_excl = (`r_hi' - `r_lo') / (2 * `h_excl')

display _n "=== SPARSITY COMPARISON ==="
display "Stata reported sparsity: 4883.3016"
display "Our sparsity (with exclusion, h_excl): " %12.4f `sparsity_excl'

* Also compute without exclusion for comparison
_pctile resid, percentile(`=(`tau'-`h_full')*100' `=(`tau'+`h_full')*100')
local r_lo2 = r(r1)
local r_hi2 = r(r2)
local sparsity_full = (`r_hi2' - `r_lo2') / (2 * `h_full')
display "Our sparsity (no exclusion, h_full): " %12.4f `sparsity_full'

* Show the VCE formula verification
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)
local xtxinv11 = XtX_inv[1,1]
local tau_factor = `tau' * (1 - `tau')

* Expected SE with each sparsity
local se_stata = sqrt(`tau_factor' * 4883.3016^2 * `xtxinv11')
local se_excl = sqrt(`tau_factor' * `sparsity_excl'^2 * `xtxinv11')
local se_full = sqrt(`tau_factor' * `sparsity_full'^2 * `xtxinv11')

display _n "=== SE COMPARISON ==="
display "Stata SE (residual method): 85.176144"
display "Expected SE (Stata sparsity): " %12.6f `se_stata'
display "Expected SE (our sparsity, excl): " %12.6f `se_excl'
display "Expected SE (our sparsity, full): " %12.6f `se_full'
