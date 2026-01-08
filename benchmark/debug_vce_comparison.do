/*******************************************************************************
* debug_vce_comparison.do
*
* Diagnostic comparison of cqreg vs qreg variance-covariance estimation.
* Run this to understand discrepancies between the two implementations.
*******************************************************************************/

clear all
set more off

* Load test data
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

local N = _N
display "Sample size: `N'"

********************************************************************************
* Part 1: Basic IID VCE Comparison
********************************************************************************

display _n "{hline 70}"
display "PART 1: IID VCE COMPARISON (default)"
display "{hline 70}"

* Run qreg with default VCE
quietly qreg price mpg weight, quantile(0.5)
local b_mpg_qreg = _b[mpg]
local b_wt_qreg = _b[weight]
local b_cons_qreg = _b[_cons]
local se_mpg_qreg = _se[mpg]
local se_wt_qreg = _se[weight]
local se_cons_qreg = _se[_cons]
matrix V_qreg = e(V)
local sum_adev_qreg = e(sum_adev)

display "qreg results:"
display "  b[mpg] = " %12.6f `b_mpg_qreg' "  SE = " %10.6f `se_mpg_qreg'
display "  b[weight] = " %12.6f `b_wt_qreg' "  SE = " %10.6f `se_wt_qreg'
display "  b[_cons] = " %12.6f `b_cons_qreg' "  SE = " %10.6f `se_cons_qreg'
display "  sum_adev = " %12.4f `sum_adev_qreg'

* Run cqreg with verbose output
display _n "Running cqreg..."
cqreg price mpg weight, quantile(0.5) verbose

local b_mpg_cqreg = _b[mpg]
local b_wt_cqreg = _b[weight]
local b_cons_cqreg = _b[_cons]
local se_mpg_cqreg = _se[mpg]
local se_wt_cqreg = _se[weight]
local se_cons_cqreg = _se[_cons]

display _n "cqreg results:"
display "  b[mpg] = " %12.6f `b_mpg_cqreg' "  SE = " %10.6f `se_mpg_cqreg'
display "  b[weight] = " %12.6f `b_wt_cqreg' "  SE = " %10.6f `se_wt_cqreg'
display "  b[_cons] = " %12.6f `b_cons_cqreg' "  SE = " %10.6f `se_cons_cqreg'

* Comparison
display _n "{hline 70}"
display "COMPARISON:"
display "{hline 70}"
display "                 qreg          cqreg       coef_diff     SE_ratio"
display "mpg:       " %12.6f `b_mpg_qreg' %12.6f `b_mpg_cqreg' %12.6f (`b_mpg_cqreg'-`b_mpg_qreg') %10.4f (`se_mpg_cqreg'/`se_mpg_qreg')
display "weight:    " %12.6f `b_wt_qreg' %12.6f `b_wt_cqreg' %12.6f (`b_wt_cqreg'-`b_wt_qreg') %10.4f (`se_wt_cqreg'/`se_wt_qreg')
display "_cons:     " %12.6f `b_cons_qreg' %12.6f `b_cons_cqreg' %12.6f (`b_cons_cqreg'-`b_cons_qreg') %10.4f (`se_cons_cqreg'/`se_cons_qreg')

********************************************************************************
* Part 2: Extract qreg's Implied Sparsity
********************************************************************************

display _n "{hline 70}"
display "PART 2: SPARSITY ANALYSIS"
display "{hline 70}"

* The IID VCE formula is:
*   V = s^2 * tau*(1-tau) * (X'X)^{-1}
*
* Where s = sparsity = 1/f(0)
*
* We can back-calculate s from qreg's V and our computed (X'X)^{-1}

* Compute (X'X)^{-1}
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)

* For tau = 0.5, tau*(1-tau) = 0.25
local tau = 0.5
local tau_factor = `tau' * (1 - `tau')

* Back-calculate sparsity from qreg's variance of first coefficient
* V[1,1] = s^2 * tau*(1-tau) * (X'X)^{-1}[1,1]
* s^2 = V[1,1] / (tau*(1-tau) * (X'X)^{-1}[1,1])
local v11_qreg = V_qreg[1,1]
local xtxinv11 = XtX_inv[1,1]
local implied_s2 = `v11_qreg' / (`tau_factor' * `xtxinv11')
local implied_s = sqrt(`implied_s2')

display "qreg implied sparsity (from mpg variance):"
display "  V[mpg,mpg] = " %12.6f `v11_qreg'
display "  (X'X)^{-1}[1,1] = " %12.6e `xtxinv11'
display "  tau*(1-tau) = " %8.6f `tau_factor'
display "  Implied s^2 = " %12.4f `implied_s2'
display "  Implied s = " %12.4f `implied_s'

* Get cqreg's reported sparsity
capture scalar cqreg_sparsity = __cqreg_sparsity
capture scalar cqreg_bandwidth = __cqreg_bandwidth

display _n "cqreg reported values:"
display "  Sparsity = " %12.4f cqreg_sparsity
display "  Bandwidth = " %12.6f cqreg_bandwidth

local sparsity_ratio = cqreg_sparsity / `implied_s'
display _n "Sparsity ratio (cqreg/qreg implied): " %8.4f `sparsity_ratio'

********************************************************************************
* Part 3: Test Different Quantiles
********************************************************************************

display _n "{hline 70}"
display "PART 3: MULTIPLE QUANTILES"
display "{hline 70}"

foreach q in 0.25 0.5 0.75 {
    display _n "Quantile: `q'"

    quietly qreg price mpg weight, quantile(`q')
    local se_mpg_q = _se[mpg]

    quietly cqreg price mpg weight, quantile(`q')
    local se_mpg_c = _se[mpg]

    local ratio = `se_mpg_c' / `se_mpg_q'
    display "  SE[mpg] qreg: " %10.6f `se_mpg_q' "  cqreg: " %10.6f `se_mpg_c' "  ratio: " %8.4f `ratio'
}

********************************************************************************
* Part 4: Bandwidth Verification
********************************************************************************

display _n "{hline 70}"
display "PART 4: BANDWIDTH CALCULATION VERIFICATION"
display "{hline 70}"

* Hall-Sheather bandwidth formula (from R's quantreg):
* h = n^(-1/3) * qnorm(1 - alpha/2)^(2/3) * ((1.5 * f0^2) / (2*x0^2 + 1))^(1/3)
*
* Where: x0 = qnorm(p), f0 = dnorm(x0), alpha = 0.05

local p = 0.5
local alpha = 0.05
local x0 = invnormal(`p')
local f0 = normalden(`x0')
local z_alpha = invnormal(1 - `alpha'/2)

local h_correct = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `f0'^2) / (2*`x0'^2 + 1))^(1/3)

display "Hall-Sheather bandwidth calculation:"
display "  n = " %d `N'
display "  p = " %6.4f `p'
display "  x0 = qnorm(p) = " %10.6f `x0'
display "  f0 = dnorm(x0) = " %10.6f `f0'
display "  z_{1-alpha/2} = " %10.6f `z_alpha'
display "  CORRECT h = " %10.6f `h_correct'
display "  cqreg h = " %10.6f cqreg_bandwidth
display "  Ratio (cqreg/correct) = " %8.4f (cqreg_bandwidth / `h_correct')

********************************************************************************
* Part 5: Test with Robust VCE
********************************************************************************

display _n "{hline 70}"
display "PART 5: ROBUST VCE COMPARISON"
display "{hline 70}"

capture {
    quietly qreg price mpg weight, quantile(0.5) vce(robust)
    local se_mpg_qreg_robust = _se[mpg]
    local se_wt_qreg_robust = _se[weight]
    local se_cons_qreg_robust = _se[_cons]

    display "qreg vce(robust):"
    display "  SE[mpg] = " %10.6f `se_mpg_qreg_robust'
    display "  SE[weight] = " %10.6f `se_wt_qreg_robust'
    display "  SE[_cons] = " %10.6f `se_cons_qreg_robust'

    quietly cqreg price mpg weight, quantile(0.5) vce(robust)
    local se_mpg_cqreg_robust = _se[mpg]
    local se_wt_cqreg_robust = _se[weight]
    local se_cons_cqreg_robust = _se[_cons]

    display _n "cqreg vce(robust):"
    display "  SE[mpg] = " %10.6f `se_mpg_cqreg_robust'
    display "  SE[weight] = " %10.6f `se_wt_cqreg_robust'
    display "  SE[_cons] = " %10.6f `se_cons_cqreg_robust'

    display _n "Ratios (cqreg/qreg):"
    display "  mpg: " %8.4f (`se_mpg_cqreg_robust'/`se_mpg_qreg_robust')
    display "  weight: " %8.4f (`se_wt_cqreg_robust'/`se_wt_qreg_robust')
    display "  _cons: " %8.4f (`se_cons_cqreg_robust'/`se_cons_qreg_robust')
}

********************************************************************************
* Part 6: Test with Cluster VCE
********************************************************************************

display _n "{hline 70}"
display "PART 6: CLUSTER VCE COMPARISON"
display "{hline 70}"

capture {
    quietly qreg price mpg weight, quantile(0.5) vce(cluster rep78)
    local se_mpg_qreg_cluster = _se[mpg]
    local se_wt_qreg_cluster = _se[weight]
    local se_cons_qreg_cluster = _se[_cons]

    display "qreg vce(cluster rep78):"
    display "  SE[mpg] = " %10.6f `se_mpg_qreg_cluster'
    display "  SE[weight] = " %10.6f `se_wt_qreg_cluster'
    display "  SE[_cons] = " %10.6f `se_cons_qreg_cluster'

    quietly cqreg price mpg weight, quantile(0.5) vce(cluster rep78)
    local se_mpg_cqreg_cluster = _se[mpg]
    local se_wt_cqreg_cluster = _se[weight]
    local se_cons_cqreg_cluster = _se[_cons]

    display _n "cqreg vce(cluster rep78):"
    display "  SE[mpg] = " %10.6f `se_mpg_cqreg_cluster'
    display "  SE[weight] = " %10.6f `se_wt_cqreg_cluster'
    display "  SE[_cons] = " %10.6f `se_cons_cqreg_cluster'

    display _n "Ratios (cqreg/qreg):"
    display "  mpg: " %8.4f (`se_mpg_cqreg_cluster'/`se_mpg_qreg_cluster')
    display "  weight: " %8.4f (`se_wt_cqreg_cluster'/`se_wt_qreg_cluster')
    display "  _cons: " %8.4f (`se_cons_cqreg_cluster'/`se_cons_qreg_cluster')
}

********************************************************************************
* Summary
********************************************************************************

display _n "{hline 70}"
display "SUMMARY"
display "{hline 70}"
display "Key finding: If SE ratios are approximately 2.0, check bandwidth calculation."
display "The h *= 2.0 scaling factor in cqreg_bandwidth_hsheather may be the issue."
display "{hline 70}"
