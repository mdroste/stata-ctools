*! version 1.1.0 25sep2026 github.com/mdroste/stata-ctools
program define cpplmhdfe_p
    version 14.1
    if "`e(cmd)'" != "cpplmhdfe" error 301
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
    if "`opt'" == "" local opt mu
    if inlist("`opt'", "cooksd", "hat", "likelihood") {
        di as error "option not implemented: `opt'"
        exit 198
    }
    if "`opt'" != "xb" & "`e(absvars)'" != "_cons" {
        if "`e(d)'" == "" {
            di as error "predict `opt' requires the d() option of cpplmhdfe"
            exit 198
        }
        confirm double variable `e(d)', exact
    }
    if "`opt'" == "stdp" {
        _predict double `varlist' `if' `in', stdp
        exit
    }
    if "`opt'" == "d" {
        if "`e(d)'" == "" quietly generate double `varlist' = 0 `if' `in'
        else quietly generate double `varlist' = `e(d)' `if' `in'
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
    if "`opt'" == "xb" exit
    if "`e(absvars)'" != "_cons" {
        quietly replace `varlist' = `varlist' + `e(d)' `if' `in'
    }
    else if "`e(d)'" != "" {
        quietly replace `varlist' = . if missing(`e(d)') `if' `in'
    }
    if inlist("`opt'", "xbd", "eta") exit
    quietly replace `varlist' = exp(`varlist') `if' `in'
    if "`opt'" == "mu" exit
    local y `e(depvar)'
    confirm variable `y'
    if inlist("`opt'", "response", "scores") ///
        quietly replace `varlist' = `y' - `varlist' `if' `in'
    else if "`opt'" == "pearson" ///
        quietly replace `varlist' = (`y' - `varlist') / sqrt(`varlist') `if' `in'
    else if "`opt'" == "deviance" ///
        quietly replace `varlist' = 2*cond(`y'>0, `varlist'-`y'+`y'*ln(`y'/`varlist'), `varlist') `if' `in'
    else if "`opt'" == "anscombe" ///
        quietly replace `varlist' = 1.5*(`y'^(2/3)-`varlist'^(2/3))/`varlist'^(1/6) `if' `in'
    else if "`opt'" == "working" ///
        quietly replace `varlist' = (`y'-`varlist')/`varlist' `if' `in'
end
