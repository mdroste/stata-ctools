*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
* predict after cppmlhdfe, with ppmlhdfe_p's statistics, default (mu, with
* its note), variable labels and notes on missing values
program define cppmlhdfe_p
    version 14.1
    if "`e(cmd)'" != "cppmlhdfe" error 301
    capture syntax anything [if] [in], SCORES
    if !_rc {
        _score_spec `anything', score
        local 0 `s(varlist)' `if' `in', scores
    }
    syntax newvarname [if] [in] [, XB XBD D Mu Eta STDP ///
        Anscombe Cooksd DEViance Hat Likelihood Pearson Response Scores Working]
    local opt `xb' `xbd' `d' `mu' `eta' `stdp' `anscombe' `cooksd' ///
        `deviance' `hat' `likelihood' `pearson' `response' `scores' `working'
    opts_exclusive "`opt'"
    if "`opt'" == "" {
        di as text "(option mu assumed; predicted mean of depvar)"
        local opt mu
    }
    if inlist("`opt'", "cooksd", "hat", "likelihood") {
        di as error "option not implemented: `opt'"
        exit 198
    }
    if "`opt'" != "xb" & "`e(absvars)'" != "_cons" {
        if "`e(d)'" == "" {
            di as error "predict `opt' requires the d() option of cppmlhdfe"
            exit 198
        }
        confirm double variable `e(d)', exact
    }
    marksample touse, novarlist
    local if "if `touse'"
    local in ""
    local target `varlist'
    tempvar prediction
    local varlist `prediction'
    if "`opt'" == "stdp" {
        _predict double `varlist' `if' `in', stdp
        rename `varlist' `target'
        exit
    }
    if "`opt'" == "d" {
        if "`e(d)'" == "" generate double `varlist' = 0 `if' `in'
        else generate double `varlist' = `e(d)' `if' `in'
        label variable `varlist' "d[`e(absvars)']"
        rename `varlist' `target'
        exit
    }
    * _predict adds e(offset). Exposure needs its logarithm instead.
    _predict double `varlist' `if' `in', xb nooffset
    if "`e(exposure)'" != "" {
        quietly replace `varlist' = `varlist' + ln(`e(exposure)') `if' `in'
    }
    else if "`e(offset)'" != "" {
        quietly replace `varlist' = `varlist' + `e(offset)' `if' `in'
    }
    if "`opt'" == "xb" {
        label variable `varlist' "Linear prediction"
        rename `varlist' `target'
        exit
    }
    if "`e(absvars)'" != "_cons" {
        quietly replace `varlist' = `varlist' + `e(d)' `if' `in'
        label variable `varlist' "Linear prediction: xb + d[`e(absvars)']"
    }
    else {
        if "`e(d)'" != "" quietly replace `varlist' = . if missing(`e(d)') & `touse'
        label variable `varlist' "Linear prediction: xb"
    }
    if inlist("`opt'", "xbd", "eta") {
        rename `varlist' `target'
        exit
    }
    if "`opt'" == "mu" {
        quietly replace `varlist' = exp(`varlist') `if' `in'
        label variable `varlist' "Predicted mean of `e(depvar)'"
        rename `varlist' `target'
        exit
    }
    tempvar fitted
    generate double `fitted' = exp(`varlist') `if' `in'
    local y `e(depvar)'
    quietly tsrevar `y'
    local y `r(varlist)'
    if "`opt'" == "deviance" {
        quietly replace `varlist' = 2*cond(`y'>0, `fitted'-`y'+`y'*ln(`y'/`fitted'), `fitted') `if' `in'
        label variable `varlist' "deviance residual"
    }
    else if "`opt'" == "response" {
        quietly replace `varlist' = `y' - `fitted' `if' `in'
        label variable `varlist' "response residual"
    }
    else if "`opt'" == "pearson" {
        quietly replace `varlist' = (`y' - `fitted') / sqrt(`fitted') `if' `in'
        label variable `varlist' "Pearson residual"
    }
    else if "`opt'" == "anscombe" {
        quietly replace `varlist' = 1.5*(`y'^(2/3)-`fitted'^(2/3))/`fitted'^(1/6) `if' `in'
        label variable `varlist' "Anscombe residual"
    }
    else if "`opt'" == "scores" {
        quietly replace `varlist' = `y' - `fitted' `if' `in'
        label variable `varlist' "score residual"
    }
    else if "`opt'" == "working" {
        quietly replace `varlist' = (`y' - `fitted') / `fitted' `if' `in'
        label variable `varlist' "working residual"
    }
    rename `varlist' `target'
end
