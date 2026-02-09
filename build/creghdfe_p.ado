*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

program define creghdfe_p, rclass
    version 14.1

    * Handle score option (used by margins)
    capture syntax anything [if] [in], SCore
    local was_score = !c(rc)
    if `was_score' {
        _score_spec `anything', score
        local 0 `s(varlist)' `if' `in', residuals
    }

    syntax newvarname [if] [in] [, XB STDP Residuals D XBD DResiduals]

    * Ensure there is only one option
    opts_exclusive "`xb' `stdp' `residuals' `d' `xbd' `dresiduals'"

    * Default option is xb
    capture opts_exclusive "`xb' `stdp' `residuals' `d' `xbd' `dresiduals' placeholder"
    if !c(rc) {
        di as text "(option xb assumed; fitted values)"
        local xb "xb"
    }

    local fixed_effects "`e(absorb)'"

    * Except for xb and stdp, we need the previously computed residuals
    if "`xb'" == "" & "`stdp'" == "" {
        if "`e(resid)'" == "" {
            di as error "you must add the {bf:resid} option to creghdfe before running this prediction"
            exit 9
        }
        capture confirm numeric variable `e(resid)', exact
        if _rc {
            di as error "residual variable `e(resid)' not found"
            exit 111
        }
    }

    if "`xb'" != "" | "`stdp'" != "" {
        * xb/stdp: use standard predict
        PredictXB `varlist' `if' `in', `xb' `stdp'
    }
    else if "`residuals'" != "" {
        * residuals: return the preexisting residual variable
        gen double `varlist' = `e(resid)' `if' `in'
        label var `varlist' "Residuals"
        if `was_score' {
            return local scorevars `varlist'
        }
    }
    else if "`d'" != "" {
        * d: y - xb - resid = fixed effects
        tempvar xbvar
        PredictXB `xbvar' `if' `in', xb
        gen double `varlist' = `e(depvar)' - `xbvar' - `e(resid)' `if' `in'
        label var `varlist' "d[`fixed_effects']"
    }
    else if "`xbd'" != "" {
        * xbd: y - resid = xb + d = fitted values including FE
        gen double `varlist' = `e(depvar)' - `e(resid)' `if' `in'
        label var `varlist' "Xb + d[`fixed_effects']"
    }
    else if "`dresiduals'" != "" {
        * dresiduals: y - xb = d + resid
        tempvar xbvar
        PredictXB `xbvar' `if' `in', xb
        gen double `varlist' = `e(depvar)' - `xbvar' `if' `in'
        label var `varlist' "d + e"
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
