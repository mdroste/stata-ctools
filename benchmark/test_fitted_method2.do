clear all
set more off
sysuse auto, clear

quietly qreg price mpg
predict resid, residuals
gen fitted = price - resid

local N = 74
local q = 0.5
local h_prob = 0.23141524

display "============================================================"
display "More methods to find qreg fitted sparsity = 6650.08"
display "============================================================"

* Method 10: Using full range (y_{q+h} - y_{q-h}) as bandwidth
_pctile price, p(`=(`q'-`h_prob')*100' `=(`q'+`h_prob')*100')
local y_range = r(r2) - r(r1)
display ""
display "Method 10: Kernel on residuals, h = full Y range"
display "  Y range (y_{q+h} - y_{q-h}) = " `y_range'

gen kernel_val = 0
replace kernel_val = 0.75 * (1 - (resid/`y_range')^2) if abs(resid) < `y_range'
qui sum kernel_val
local sum_K = r(sum)
local f0 = `sum_K' / (`N' * `y_range')
local s = 1 / `f0'
display "  f(0) = " `f0'
display "  sparsity = " `s'
drop kernel_val

* Method 11: Try different kernel normalizations
display ""
display "Method 11: Raw difference quotient on Y (no /2h)"
_pctile price, p(`=(`q'-`h_prob')*100' `=(`q'+`h_prob')*100')
local y_lo = r(r1)
local y_hi = r(r2)
local s_raw = `y_hi' - `y_lo'
display "  y_hi - y_lo = " `s_raw'
local s_div_h = `s_raw' / `h_prob'
display "  (y_hi - y_lo) / h = " `s_div_h'

* Method 12: Using 1/(2*h*f(xq)) formula directly
* Let me try computing f(xq) at the median Y value
_pctile price, p(50)
local y_med = r(r1)
display ""
display "Method 12: Kernel density of Y at median Y"
display "  Median Y = " `y_med'

* Bandwidth for Y
local h_y = (6303 - 4296) / 2
display "  h in Y scale = " `h_y'

gen kernel_val = 0
replace kernel_val = 0.75 * (1 - ((price - `y_med')/`h_y')^2) if abs(price - `y_med') < `h_y'
qui sum kernel_val
local sum_K = r(sum)
local f_y_med = `sum_K' / (`N' * `h_y')
local s_y_med = 1 / `f_y_med'
display "  f(y_med) = " `f_y_med'
display "  sparsity = " `s_y_med'
drop kernel_val

* Method 13: Use twice the half-range (full span 2h)
* Sparsity formula: s = (x_{q+h} - x_{q-h}) / (2h) where x is residuals
* But for fitted, maybe s = (y_{q+h} - y_{q-h}) / (2h) where y is price?
display ""
display "Method 13: Y range divided by 2h (probability scale)"
local s_y_2h = `y_range' / (2 * `h_prob')
display "  sparsity = (y_hi - y_lo) / (2h) = " `s_y_2h'

* Method 14: Try using PREDICTD quantile formula from Koenker
* Linear interpolation between order statistics
display ""
display "Method 14: Linear interpolation quantiles"

* Sort and get order statistics
sort price
local n_lo = ceil(`N' * (`q' - `h_prob'))
local n_hi = ceil(`N' * (`q' + `h_prob'))
local y_lo_ord = price[`n_lo']
local y_hi_ord = price[`n_hi']
local s_ord = (`y_hi_ord' - `y_lo_ord') / (2 * `h_prob')
display "  n_lo = " `n_lo' ", n_hi = " `n_hi'
display "  y_lo_ord = " `y_lo_ord' ", y_hi_ord = " `y_hi_ord'
display "  sparsity = " `s_ord'

* Method 15: Try using fitted values range but on different scale
display ""
display "Method 15: Fitted range with Y-to-fitted scaling"
_pctile fitted, p(`=(`q'-`h_prob')*100' `=(`q'+`h_prob')*100')
local f_range = r(r2) - r(r1)
display "  Fitted range = " `f_range'

* Compute ratio of Y spread to fitted spread
_pctile price, p(25 75)
local y_iqr = r(r2) - r(r1)
_pctile fitted, p(25 75)
local f_iqr = r(r2) - r(r1)
local scale = `y_iqr' / `f_iqr'
display "  Y IQR = " `y_iqr' ", Fitted IQR = " `f_iqr' ", scale = " `scale'
local s_scaled = (`f_range' * `scale') / (2 * `h_prob')
display "  Scaled fitted sparsity = " `s_scaled'

* Method 16: Try Gaussian kernel instead of Epanechnikov
display ""
display "Method 16: Gaussian kernel on residuals"
local h_resid = 1041.6667
gen kernel_gauss = normalden(resid / `h_resid')
qui sum kernel_gauss
local sum_K = r(sum)
local f0_gauss = `sum_K' / (`N' * `h_resid')
local s_gauss = 1 / `f0_gauss'
display "  f(0) = " `f0_gauss'
display "  sparsity = " `s_gauss'
drop kernel_gauss

* Method 17: Check if Stata uses a different effective N
* Maybe Stata divides by a different value?
display ""
display "Method 17: Testing different N divisors"
local f0_base = 0.00028277  /* from residual kernel */
forval divisor = 50/100 {
    local s_test = 1 / (`f0_base' * `divisor' / 74)
    if abs(`s_test' - 6650) < 50 {
        display "  N/divisor = " `divisor' ", sparsity = " `s_test'
    }
}

* Method 18: What density gives exactly 6650?
local target_f = 1/6650.08
display ""
display "Method 18: Required density for sparsity = 6650.08"
display "  Required f(0) = " `target_f'
display "  With h=1041.67, required sum_K = " (`target_f' * 74 * 1041.67)

display ""
display "============================================================"
display "Summary"
display "============================================================"
display "Target: qreg IID fitted sparsity = 6650.08"
display "Method 10 (full Y range): " %8.2f `s'
display "Method 11a (raw Y diff): " %8.2f `s_raw'
display "Method 11b (Y diff / h): " %8.2f `s_div_h'
display "Method 12 (KDE Y at median): " %8.2f `s_y_med'
display "Method 13 (Y range / 2h): " %8.2f `s_y_2h'
display "Method 14 (order stat): " %8.2f `s_ord'
display "Method 15 (scaled fitted): " %8.2f `s_scaled'
display "Method 16 (Gaussian): " %8.2f `s_gauss'
