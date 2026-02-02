*! version 0.9.0 26Jan2026
*! cqreg_p: Prediction after cqreg
*! Part of the ctools suite

program define cqreg_p
    version 14.1

    syntax newvarname [if] [in], [XB Residuals STDP]

    * Default is xb (linear prediction)
    local nopt = ("`xb'" != "") + ("`residuals'" != "") + ("`stdp'" != "")
    if `nopt' > 1 {
        di as error "only one statistic may be specified"
        exit 198
    }
    if `nopt' == 0 {
        local xb "xb"
    }

    * Check that cqreg was run
    if "`e(cmd)'" != "cqreg" {
        di as error "last estimates not found or not from cqreg"
        exit 301
    }

    marksample touse, novarlist

    * Get dependent variable name
    local depvar "`e(depvar)'"

    * Generate prediction
    if "`xb'" != "" {
        * Linear prediction: xb = X*beta
        tempvar linpred
        qui matrix score double `linpred' = e(b) if `touse'
        gen `typlist' `varlist' = `linpred' if `touse'
        label variable `varlist' "Linear prediction"
    }
    else if "`residuals'" != "" {
        * Residuals: r = y - xb
        tempvar linpred
        qui matrix score double `linpred' = e(b) if `touse'
        gen `typlist' `varlist' = `depvar' - `linpred' if `touse'
        label variable `varlist' "Residuals"
    }
    else if "`stdp'" != "" {
        * Standard error of prediction: sqrt(X * V * X')
        * For each observation, compute x_i' * V * x_i
        tempname V b
        matrix `V' = e(V)
        matrix `b' = e(b)
        local K = colsof(`b')

        * Get variable names from coefficient matrix
        local varnames : colnames `b'

        * Create temporary variables for X values
        tempvar stdp_val
        gen double `stdp_val' = 0 if `touse'

        * Build X matrix observation by observation (slow but correct)
        * For large datasets, this could be optimized
        local i = 0
        foreach v of local varnames {
            local i = `i' + 1
            if "`v'" == "_cons" {
                * Constant term
                tempvar x`i'
                gen double `x`i'' = 1 if `touse'
            }
            else {
                * Try to use the variable directly
                capture confirm variable `v'
                if _rc == 0 {
                    tempvar x`i'
                    gen double `x`i'' = `v' if `touse'
                }
                else {
                    * Handle factor variables - try original name
                    local vclean = subinstr("`v'", "_", ".", .)
                    capture confirm variable `vclean'
                    if _rc == 0 {
                        tempvar x`i'
                        gen double `x`i'' = `vclean' if `touse'
                    }
                    else {
                        * Use 0 for unavailable variables
                        tempvar x`i'
                        gen double `x`i'' = 0 if `touse'
                    }
                }
            }
        }

        * Compute quadratic form: sum_j sum_k x_j * V[j,k] * x_k
        forval j = 1/`K' {
            forval k = 1/`K' {
                local Vjk = `V'[`j', `k']
                qui replace `stdp_val' = `stdp_val' + `x`j'' * `Vjk' * `x`k'' if `touse'
            }
        }

        * Take square root
        gen `typlist' `varlist' = sqrt(`stdp_val') if `touse'
        label variable `varlist' "Std. error of prediction"
    }

end
