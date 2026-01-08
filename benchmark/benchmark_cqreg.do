* Benchmark cqreg vs qreg - Standard Error Comparison
* Goal: Match standard errors to machine precision (diff < 1E-5)

clear all
set more off
adopath + "../build"

sysuse auto, clear

di _newline "{hline 78}"
di "TEST 1: QREG price mpg"
di "{hline 78}"
qreg price mpg

matrix qreg_b = e(b)
matrix qreg_V = e(V)
local qreg_N = e(N)
local qreg_df_r = e(df_r)
local qreg_q = e(q)
local qreg_q_v = e(q_v)
local qreg_f = e(f)
local qreg_convcode = e(convcode)
local qreg_sum_adev = e(sum_adev)
local qreg_sum_rdev = e(sum_rdev)
local qreg_r2_p = e(r2_p)

di _newline(2) "{hline 78}"
di "CQREG price mpg"
di "{hline 78}"
cqreg price mpg

matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
local cqreg_N = e(N)
local cqreg_df_r = e(df_r)
local cqreg_q = e(q)
local cqreg_q_v = e(q_v)
local cqreg_sum_adev = e(sum_adev)
local cqreg_sum_rdev = e(sum_rdev)
local cqreg_r2_p = e(r2_p)

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
di "VARIANCE-COVARIANCE MATRIX COMPARISON"
di "{hline 78}"
di "qreg V[1,1]:   " %14.6f qreg_V[1,1] "  cqreg V[1,1]:   " %14.6f cqreg_V[1,1] "  diff: " %14.6f qreg_V[1,1]-cqreg_V[1,1]
di "qreg V[1,2]:   " %14.6f qreg_V[1,2] "  cqreg V[1,2]:   " %14.6f cqreg_V[1,2] "  diff: " %14.6f qreg_V[1,2]-cqreg_V[1,2]
di "qreg V[2,2]:   " %14.6f qreg_V[2,2] "  cqreg V[2,2]:   " %14.6f cqreg_V[2,2] "  diff: " %14.6f qreg_V[2,2]-cqreg_V[2,2]

di _newline "{hline 78}"
di "SCALAR COMPARISON"
di "{hline 78}"
di "              qreg            cqreg           diff"
di "{hline 60}"
di "N:        " %12.0f `qreg_N' "  " %12.0f `cqreg_N' "  " %12.0f `qreg_N'-`cqreg_N'
di "df_r:     " %12.0f `qreg_df_r' "  " %12.0f `cqreg_df_r' "  " %12.0f `qreg_df_r'-`cqreg_df_r'
di "q:        " %12.6f `qreg_q' "  " %12.6f `cqreg_q' "  " %12.6f `qreg_q'-`cqreg_q'
di "q_v:      " %12.6f `qreg_q_v' "  " %12.6f `cqreg_q_v' "  " %12.6f `qreg_q_v'-`cqreg_q_v'
di "sum_adev: " %12.6f `qreg_sum_adev' "  " %12.6f `cqreg_sum_adev' "  " %12.6f `qreg_sum_adev'-`cqreg_sum_adev'
di "sum_rdev: " %12.6f `qreg_sum_rdev' "  " %12.6f `cqreg_sum_rdev' "  " %12.6f `qreg_sum_rdev'-`cqreg_sum_rdev'
di "r2_p:     " %12.8f `qreg_r2_p' "  " %12.8f `cqreg_r2_p' "  " %12.8f `qreg_r2_p'-`cqreg_r2_p'

di _newline "{hline 78}"
di "VCE DETAILS FROM QREG"
di "{hline 78}"
di "f (density at median):  " %12.8f `qreg_f'
di "Sparsity (1/f):         " %12.8f 1/`qreg_f'

* Get residuals and estimate sparsity manually
predict double resid_qreg, resid
summ resid_qreg, detail

* Compute sparsity using Stata's method
* qreg uses kdens option, let's check the residual distribution
_pctile resid_qreg, percentile(25 50 75)
local q25 = r(r1)
local q50 = r(r2)
local q75 = r(r3)
di "Residual quartiles: Q25=" %12.4f `q25' " Q50=" %12.4f `q50' " Q75=" %12.4f `q75'

* Show what VCE cqreg should match
di _newline "{hline 78}"
di "MANUAL VCE CALCULATION (should match qreg)"
di "{hline 78}"

* Number of obs
local N = 74
local K = 2  // number of regressors including constant
local q = 0.5

* Hall-Sheather bandwidth (with alpha=0.05)
* h = n^{-1/3} * z_{1-α/2}^{2/3} * (1.5 * φ(z_q)^2 / (2*z_{1-α/2}^2 + 1))^{1/3}
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
di "z_alpha = " %12.8f `z_alpha'
di "z_q = " %12.8f `z_q'
di "phi(z_q) = " %12.8f `phi_zq'

local h_raw = (`z_alpha'^(2/3)) * ((1.5 * `phi_zq'^2) / (2*`z_alpha'^2 + 1))^(1/3) * `N'^(-1/3)
di "h_raw (no scaling) = " %12.8f `h_raw'

* Get residual quantiles for sparsity estimation
local q_lo = `q' - `h_raw'
local q_hi = `q' + `h_raw'
if `q_lo' < 0 local q_lo = 0
if `q_hi' > 1 local q_hi = 1
di "q_lo = " %12.8f `q_lo' " q_hi = " %12.8f `q_hi'

_pctile resid_qreg, percentile(`=`q_lo'*100' `=`q_hi'*100')
local x_lo = r(r1)
local x_hi = r(r2)
di "x_lo = " %12.4f `x_lo' " x_hi = " %12.4f `x_hi'

local sparsity_manual = (`x_hi' - `x_lo') / (`q_hi' - `q_lo')
di "Manual sparsity (no 2x scale): " %12.4f `sparsity_manual'
di "Manual sparsity (with 2x scale): " %12.4f `sparsity_manual' * 2

* Implied sparsity from qreg VCE
* V = sparsity^2 * q*(1-q) * (X'X)^{-1}
* So sparsity = sqrt(V[1,1] / (q*(1-q) * (X'X)^{-1}[1,1]))
* Need to compute (X'X)^{-1}
matrix accum XtX = mpg, noconstant
matrix XtX_inv = invsym(XtX)
di "X'X = " %12.4f XtX[1,1]
di "(X'X)^{-1} = " %12.8f XtX_inv[1,1]

* V[1,1] = sparsity^2 * q*(1-q) * (X'X)^{-1}[1,1]
* Need to compute (X'X)^{-1} with constant term
gen byte _one = 1
matrix accum XtX_full = mpg _one, noconstant
matrix XtX_inv_full = invsym(XtX_full)
matrix list XtX_full
matrix list XtX_inv_full

local scale = `q' * (1-`q')
di "q*(1-q) = " %12.8f `scale'

* Implied sparsity from V[1,1]
local implied_sparsity_mpg = sqrt(qreg_V[1,1] / (`scale' * XtX_inv_full[1,1]))
local implied_sparsity_cons = sqrt(qreg_V[2,2] / (`scale' * XtX_inv_full[2,2]))
di "Implied sparsity from V[1,1] (mpg): " %12.4f `implied_sparsity_mpg'
di "Implied sparsity from V[2,2] (_cons): " %12.4f `implied_sparsity_cons'

drop _one

drop resid_qreg
