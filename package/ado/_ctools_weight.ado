*! ctools weight-expression staging; caller owns generate() tempvar
program define _ctools_weight
    version 14.1
    syntax, Generate(name) TOUse(varname) TYPE(string) EXPression(string)
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
