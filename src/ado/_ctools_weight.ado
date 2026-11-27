*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! ctools weight-expression staging; caller owns generate() tempvar
* Sets the caller's local __ctools_wvar to the variable that holds the
* weights: generate(), filled here and checked. With option marked the
* caller's marksample has seen this weight: the sample then leaves out zero
* and missing weights, and Stata has refused negative weights (r(402)) and
* noninteger frequency weights (r(401)). A weight that is one numeric
* variable is then used as it is, and generate() is not created.
program define _ctools_weight
    version 14.1
    syntax, Generate(name) TOUse(varname) TYPE(string) EXPression(string) [MARKED]
    c_local __ctools_wvar `generate'
    if "`marked'" != "" {
        local wexp = strtrim(substr(strtrim(`"`expression'"'), 2, .))
        if `: word count `wexp'' == 1 {
            capture confirm numeric variable `wexp', exact
            if !_rc {
                c_local __ctools_wvar `wexp'
                exit
            }
        }
    }
    confirm new variable `generate'
    quietly generate double `generate' `expression' if `touse'
    markout `touse' `generate'
    quietly count if `generate' <= 0 & `touse'
    if r(N) {
        di as error "ctools: weights must be positive"
        exit 198
    }
    if "`type'" == "fweight" {
        quietly count if `generate' != floor(`generate') & `touse'
        if r(N) {
            di as error "frequency weights must be integers"
            exit 401
        }
    }
end
