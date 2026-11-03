*! version 1.2.4 20261103 github.com/mdroste/stata-ctools

program define civreghdfe_p, rclass
    version 14.1

    * Handle score option (used by margins)
    capture syntax anything [if] [in], SCore
    local was_score = !c(rc)
    if `was_score' {
        _score_spec `anything', score
        local 0 `s(varlist)' `if' `in', residuals
    }

    syntax newvarname [if] [in] [, XB STDP Residuals XBD D DResiduals]

    * Ensure there is only one option
    opts_exclusive "`xb' `stdp' `residuals' `xbd' `d' `dresiduals'"

    * Default option is xb
    capture opts_exclusive "`xb' `stdp' `residuals' `xbd' `d' `dresiduals' placeholder"
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

    * After partial(), e(b) omits the partialled-out terms. Without absorb()
    * only residuals can be predicted, as ivreg2; with absorb() xb and stdp
    * refer to the posted coefficients, as ivreghdfe's reghdfe_p
    local partial_ct = cond(missing(e(partial_ct)), 0, e(partial_ct))
    if `partial_ct' > 0 & "`fixed_effects'" == "" & "`residuals'" == "" {
        di as error "option `xb'`stdp'`xbd'`d'`dresiduals' not supported with the partial() option"
        exit 198
    }

    * d (the sum of the fixed effects, y - xb - residuals) and dresiduals
    * (y - xb) need the residuals saved by residuals(), as in reghdfe_p
    if "`d'`dresiduals'" != "" & "`e(resid)'" == "" {
        di as error "you must add the residuals() option to civreghdfe before running this prediction"
        exit 198
    }

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
        if "`e(resid)'" != "" {
            confirm numeric variable `e(resid)', exact
            gen `typlist' `varlist' = `e(resid)' if `touse' & e(sample)
        }
        else if "`fixed_effects'" != "" {
            di as error "residual predictions with fixed effects require residuals() during estimation"
            exit 198
        }
        else if `partial_ct' > 0 {
            * Frisch-Waugh-Lovell: the residual is y - Xb net of its
            * projection on the partialled-out terms, on the estimation
            * sample with its weights (with the constant: in e(b), or
            * partialled out with the other terms)
            tempvar xbvar ypart esample rpart
            tempname ivest
            PredictXB `xbvar' `if' `in', xb
            quietly gen double `ypart' = `depvar' - `xbvar'
            quietly gen byte `esample' = e(sample)
            local pwt
            if "`e(wtype)'" != "" local pwt [`e(wtype)'`e(wexp)']
            local pcons = cond(colnumb(e(b), "_cons") < ., "noconstant", "")
            * a partialled-out constant (no absorb()) is the regression's own
            local pterms `e(partial)'
            local __cons_term _cons
            local pterms : list pterms - __cons_term
            _estimates hold `ivest', restore
            quietly regress `ypart' `pterms' `pwt' if `esample', `pcons'
            quietly predict double `rpart' if `esample', residuals
            _estimates unhold `ivest'
            gen `typlist' `varlist' = `rpart' if `touse' & e(sample)
        }
        else {
            tempvar xbvar
            PredictXB `xbvar' `if' `in', xb
            gen `typlist' `varlist' = `depvar' - `xbvar' if `touse' & e(sample)
        }
        label var `varlist' "Residuals"
        if `was_score' {
            return local scorevars `varlist'
        }
    }
    else if "`d'" != "" {
        confirm numeric variable `e(resid)', exact
        tempvar xbvar
        PredictXB `xbvar' `if' `in', xb
        gen `typlist' `varlist' = `depvar' - `xbvar' - `e(resid)' if `touse'
        label var `varlist' "d[`fixed_effects']"
    }
    else if "`dresiduals'" != "" {
        confirm numeric variable `e(resid)', exact
        tempvar xbvar
        PredictXB `xbvar' `if' `in', xb
        gen `typlist' `varlist' = `depvar' - `xbvar' if `touse'
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
        gen `typlist' `varlist' = `depvar' - `e(resid)' if `touse' & e(sample)
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
