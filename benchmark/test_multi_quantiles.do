clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display "============================================="
display "TESTING CQREG VS QREG VCE(IID, RESIDUAL)"
display "ACROSS MULTIPLE QUANTILES"
display "============================================="
display ""
display "      q      SE[mpg]_qreg  SE[mpg]_cqreg     ratio"
display "  -------  ------------  -------------  ---------"

local all_match = 1

foreach q in 0.10 0.25 0.50 0.75 0.90 {
    quietly qreg price mpg weight, quantile(`q') vce(iid, residual)
    local se_qreg = _se[mpg]
    
    quietly cqreg price mpg weight, quantile(`q')
    local se_cqreg = _se[mpg]
    
    local ratio = `se_cqreg' / `se_qreg'
    
    display "   " %4.2f `q' "     " %12.4f `se_qreg' "   " %12.4f `se_cqreg' "    " %8.6f `ratio'
    
    if abs(`ratio' - 1.0) > 0.001 {
        local all_match = 0
    }
}

display ""
if `all_match' {
    display "{result}SUCCESS: All quantiles match within 0.1%{txt}"
}
else {
    display "{error}WARNING: Some quantiles differ by more than 0.1%{txt}"
}

display ""
display "============================================="
display "TESTING WITH LARGER DATASET (nlsw88)"
display "============================================="

sysuse nlsw88, clear
drop if missing(wage, ttl_exp, tenure, grade)

display ""
display "      q      SE[exp]_qreg  SE[exp]_cqreg     ratio"
display "  -------  ------------  -------------  ---------"

foreach q in 0.10 0.25 0.50 0.75 0.90 {
    quietly qreg wage ttl_exp tenure grade, quantile(`q') vce(iid, residual)
    local se_qreg = _se[ttl_exp]
    
    quietly cqreg wage ttl_exp tenure grade, quantile(`q')
    local se_cqreg = _se[ttl_exp]
    
    local ratio = `se_cqreg' / `se_qreg'
    
    display "   " %4.2f `q' "     " %12.6f `se_qreg' "   " %12.6f `se_cqreg' "    " %8.6f `ratio'
}

display ""
display "============================================="
