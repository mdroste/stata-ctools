*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! Private helper for civreghdfe.
program define _civreghdfe_optflags, rclass
    args optvars all touse
    fvexpand `optvars' if `touse'
    local optnames `r(varlist)'

    * Plugin columns: the terms that are not base/omitted placeholders
    _ctools_fvsplit, names(`all')
    local is_base `r(base)'
    local sall ""
    local plugpos ""
    local np = 0
    foreach v of local all {
        gettoken base is_base : is_base
        _civreghdfe_fvstrip `"`v'"'
        local sall `sall' `r(name)'
        if !`base' {
            local ++np
            local plugpos `plugpos' `np'
        }
        else local plugpos `plugpos' 0
    }

    tempname F
    if `np' > 0 matrix `F' = J(1, `np', 0)
    local n = 0
    local terms ""
    local bad ""
    foreach v of local optnames {
        _civreghdfe_fvstrip `"`v'"'
        local sv `r(name)'
        local j : list posof "`sv'" in sall
        if `j' == 0 {
            local bad `v'
            continue, break
        }
        local terms `terms' `: word `j' of `all''
        local p : word `j' of `plugpos'
        if `p' > 0 {
            if `F'[1, `p'] == 0 local ++n
            matrix `F'[1, `p'] = 1
        }
    }
    local terms : list uniq terms
    return local bad `bad'
    return local terms `terms'
    return scalar n = `n'
    if `n' > 0 return matrix flags = `F'
end
