clear all
set more off
sysuse auto, clear

* Get qreg results
quietly qreg price mpg
local b_cons = _b[_cons]
local b_mpg = _b[mpg]
predict resid, residuals
gen fitted = price - resid

local N = 74
local q = 0.5
local alpha = 0.05

display "============================================================"
display "Understanding Stata's fitted method for sparsity estimation"
display "============================================================"
display ""
display "Coefficients: b[_cons] = " `b_cons' ", b[mpg] = " `b_mpg'
display "N = " `N' ", q = " `q'
display ""

* Stata's reported values
quietly qreg price mpg
display "qreg IID fitted sparsity = " e(sparsity) ", bandwidth = " e(bwidth)

quietly qreg price mpg, vce(iid, residual)
display "qreg IID residual sparsity = " e(sparsity) ", bandwidth = " e(bwidth)

display ""
display "============================================================"
display "Testing different sparsity estimation methods"
display "============================================================"

* Method 1: Difference quotient on residuals (residual method)
local z_alpha = invnormal(1 - `alpha'/2)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h_prob_N = `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3) * `N'^(-1/3)

* Exclude LP basis observations for residual method
count if abs(resid) < 1e-8
local n_basis = r(N)
local N_adj = `N' - `n_basis'
local h_prob_adj = `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3) * `N_adj'^(-1/3)

display ""
display "Hall-Sheather bandwidth (N=" `N' "): " `h_prob_N'
display "Hall-Sheather bandwidth (N_adj=" `N_adj' "): " `h_prob_adj'
display "LP basis observations: " `n_basis'

* Difference quotient on residuals
tempvar resid_nobase
gen `resid_nobase' = resid if abs(resid) >= 1e-8
_pctile `resid_nobase', p(`=(`q'-`h_prob_adj')*100' `=(`q'+`h_prob_adj')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local s_dq_resid = (`r_hi' - `r_lo') / (2 * `h_prob_adj')
display ""
display "Method 1: Difference quotient on residuals (excl. basis)"
display "  r_lo = " `r_lo' ", r_hi = " `r_hi'
display "  sparsity = " `s_dq_resid'

* Method 2: Difference quotient on Y
_pctile price, p(`=(`q'-`h_prob_N')*100' `=(`q'+`h_prob_N')*100')
local y_lo = r(r1)
local y_hi = r(r2)
local s_dq_y = (`y_hi' - `y_lo') / (2 * `h_prob_N')
display ""
display "Method 2: Difference quotient on Y"
display "  y_lo = " `y_lo' ", y_hi = " `y_hi'
display "  sparsity = " `s_dq_y'

* Method 3: Difference quotient on fitted values
_pctile fitted, p(`=(`q'-`h_prob_N')*100' `=(`q'+`h_prob_N')*100')
local f_lo = r(r1)
local f_hi = r(r2)
local s_dq_fit = (`f_hi' - `f_lo') / (2 * `h_prob_N')
display ""
display "Method 3: Difference quotient on fitted values"
display "  f_lo = " `f_lo' ", f_hi = " `f_hi'
display "  sparsity = " `s_dq_fit'

* Method 4: Sort Y by fitted, then take quantiles of Y
gsort fitted
tempvar y_sorted order_idx
gen `order_idx' = _n
gen `y_sorted' = price
_pctile `y_sorted', p(`=(`q'-`h_prob_N')*100' `=(`q'+`h_prob_N')*100')
local ys_lo = r(r1)
local ys_hi = r(r2)
local s_sorted = (`ys_hi' - `ys_lo') / (2 * `h_prob_N')
display ""
display "Method 4: Y quantiles after sorting by fitted"
display "  ys_lo = " `ys_lo' ", ys_hi = " `ys_hi'
display "  sparsity = " `s_sorted'

* Restore sort order
sort `order_idx'

* Method 5: Kernel density on residuals at 0 (Silverman bandwidth)
local r_25 = 0
local r_75 = 0
_pctile resid, p(25 75)
local r_25 = r(r1)
local r_75 = r(r2)
local iqr_r = `r_75' - `r_25'
local h_silverman = 0.9 * `iqr_r' / 1.34 * `N'^(-0.2)
display ""
display "Method 5: Kernel density on residuals (Silverman bandwidth)"
display "  IQR of residuals = " `iqr_r'
display "  Silverman h = " `h_silverman'

gen kernel_val = 0
replace kernel_val = 0.75 * (1 - (resid/`h_silverman')^2) if abs(resid) < `h_silverman'
qui sum kernel_val
local sum_K = r(sum)
local f0_silver = `sum_K' / (`N' * `h_silverman')
local s_silver = 1 / `f0_silver'
display "  f(0) = " `f0_silver'
display "  sparsity = " `s_silver'
drop kernel_val

* Method 6: Kernel density on residuals at 0 (Hall-Sheather bandwidth)
* Convert h_prob to residual scale using residual quantiles
_pctile resid, p(`=(`q'-`h_prob_N')*100' `=(`q'+`h_prob_N')*100')
local h_resid = (r(r2) - r(r1)) / 2
display ""
display "Method 6: Kernel density on residuals (H-S bandwidth)"
display "  h in residual scale = " `h_resid'

gen kernel_val = 0
replace kernel_val = 0.75 * (1 - (resid/`h_resid')^2) if abs(resid) < `h_resid'
qui sum kernel_val
local sum_K = r(sum)
local f0_hs = `sum_K' / (`N' * `h_resid')
local s_hs = 1 / `f0_hs'
display "  f(0) = " `f0_hs'
display "  sparsity = " `s_hs'
drop kernel_val

* Method 7: Kernel density on residuals (H-S bandwidth in FITTED scale)
_pctile fitted, p(`=(`q'-`h_prob_N')*100' `=(`q'+`h_prob_N')*100')
local h_fitted = (r(r2) - r(r1)) / 2
display ""
display "Method 7: Kernel density on residuals (H-S bandwidth in fitted scale)"
display "  h in fitted scale = " `h_fitted'

gen kernel_val = 0
replace kernel_val = 0.75 * (1 - (resid/`h_fitted')^2) if abs(resid) < `h_fitted'
qui sum kernel_val
local sum_K = r(sum)
local f0_fit = `sum_K' / (`N' * `h_fitted')
local s_fit = 1 / `f0_fit'
display "  f(0) = " `f0_fit'
display "  sparsity = " `s_fit'
drop kernel_val

* Method 8: Kernel density on residuals (H-S bandwidth in Y scale)
_pctile price, p(`=(`q'-`h_prob_N')*100' `=(`q'+`h_prob_N')*100')
local h_y = (r(r2) - r(r1)) / 2
display ""
display "Method 8: Kernel density on residuals (H-S bandwidth in Y scale)"
display "  h in Y scale = " `h_y'

gen kernel_val = 0
replace kernel_val = 0.75 * (1 - (resid/`h_y')^2) if abs(resid) < `h_y'
qui sum kernel_val
local sum_K = r(sum)
local f0_y = `sum_K' / (`N' * `h_y')
local s_y = 1 / `f0_y'
display "  f(0) = " `f0_y'
display "  sparsity = " `s_y'
drop kernel_val

* Method 9: Per-observation kernel density on Y at fitted values (Powell)
* For each observation i, estimate f_i = density of Y at yhat_i
display ""
display "Method 9: Per-observation kernel density (Powell method)"
display "  h in Y scale = " `h_y'

gen obs_density = 0
qui forval i = 1/`N' {
    local yhat_i = fitted[`i']
    tempvar kernel_contrib
    gen `kernel_contrib' = 0
    replace `kernel_contrib' = 0.75 * (1 - ((price - `yhat_i')/`h_y')^2) if abs(price - `yhat_i') < `h_y'
    qui sum `kernel_contrib'
    replace obs_density = r(sum) / (`N' * `h_y') in `i'
    drop `kernel_contrib'
}
qui sum obs_density
local avg_f = r(mean)
local s_powell = 1 / `avg_f'
display "  Average f = " `avg_f'
display "  sparsity = " `s_powell'

display ""
display "============================================================"
display "Summary - Target: qreg IID fitted sparsity = 6650.08"
display "============================================================"
display "Method 1 (DQ residuals, excl basis): " %8.2f `s_dq_resid'
display "Method 2 (DQ Y):                     " %8.2f `s_dq_y'
display "Method 3 (DQ fitted):                " %8.2f `s_dq_fit'
display "Method 4 (DQ Y sorted by fitted):    " %8.2f `s_sorted'
display "Method 5 (KDE Silverman):            " %8.2f `s_silver'
display "Method 6 (KDE H-S residual):         " %8.2f `s_hs'
display "Method 7 (KDE H-S fitted):           " %8.2f `s_fit'
display "Method 8 (KDE H-S Y):                " %8.2f `s_y'
display "Method 9 (Powell per-obs):           " %8.2f `s_powell'
