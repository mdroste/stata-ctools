*! version 0.9.0 01Feb2026
*! civreghdfe_p: Predict program for civreghdfe
*! Part of the ctools suite

program define civreghdfe_p, rclass
    version 14.1

    * Handle score option (used by margins)
    capture syntax anything [if] [in], SCore
    local was_score = !c(rc)
    if `was_score' {
        _score_spec `anything', score
        local 0 `s(varlist)' `if' `in', residuals
    }

    syntax newvarname [if] [in] [, XB STDP Residuals XBD]

    * Ensure there is only one option
    opts_exclusive "`xb' `stdp' `residuals' `xbd'"

    * Default option is xb
    capture opts_exclusive "`xb' `stdp' `residuals' `xbd' placeholder"
    if !c(rc) {
        di as text "(option xb assumed; fitted values)"
        local xb "xb"
    }

    * Check that civreghdfe was run
    if "`e(cmd)'" != "civreghdfe" {
        di as error "last estimates not found or not from civreghdfe"
        exit 301
    }

    marksample touse, novarlist

    local depvar "`e(depvar)'"
    local fixed_effects "`e(absorb)'"

    if "`xb'" != "" {
        * xb: linear prediction X*beta (structural equation)
        PredictXB `varlist' `if' `in', xb
        label var `varlist' "Fitted values"
    }
    else if "`stdp'" != "" {
        * stdp: standard error of linear prediction
        PredictXB `varlist' `if' `in', stdp
        label var `varlist' "S.E. of prediction"
    }
    else if "`residuals'" != "" {
        * residuals: y - xb
        tempvar xbvar
        PredictXB `xbvar' `if' `in', xb
        gen double `varlist' = `depvar' - `xbvar' if `touse'
        label var `varlist' "Residuals"
        if `was_score' {
            return local scorevars `varlist'
        }
    }
    else if "`xbd'" != "" {
        * xbd: fitted values including fixed effects (y - residuals)
        * This requires stored residuals from the estimation
        if "`e(resid)'" == "" {
            di as error "xbd requires residuals() option during civreghdfe estimation"
            di as error "Re-run civreghdfe with the residuals() option"
            exit 198
        }
        capture confirm numeric variable `e(resid)', exact
        if _rc {
            di as error "residual variable `e(resid)' not found"
            exit 111
        }
        gen double `varlist' = `depvar' - `e(resid)' if `touse'
        label var `varlist' "Xb + d[`fixed_effects']"
    }
    else {
        error 100
    }
end

program PredictXB
    syntax newvarname [if] [in], [*]
    * Check if there are any regressors
    capture matrix list e(b)
    if c(rc) {
        gen double `varlist' = 0 `if' `in'
    }
    else {
        _predict double `varlist' `if' `in', `options'
    }
end
