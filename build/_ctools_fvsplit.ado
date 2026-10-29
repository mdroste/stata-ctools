*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! ctools: estimated terms and base/omitted placeholders of an fvexpand-ed list
* names() is an fvexpand result. A term is a placeholder, not estimated, only
* when every #-component carries a base (b.) or omitted (o.) marker:
* 1b.a#2.b and 0b.foreign#c.weight are estimable, 1b.a#1b.b and 1b.a#co.x
* are not. vars(), when given, is the fvrevar list for the same terms.
* Returns r(base) (one 0/1 flag per term, 1 = placeholder), r(names) and
* r(vars) (the estimated terms and their variables).
program define _ctools_fvsplit, rclass
    version 14.1
    syntax [, NAMES(string) VARS(string)]
    local has_vars = (`"`vars'"' != "")
    local flags
    local kept_names
    local kept_vars
    foreach name of local names {
        if `has_vars' gettoken var vars : vars
        local base = 1
        local parts : subinstr local name "#" " ", all
        foreach part of local parts {
            if !regexm("`part'", "b\.") & !regexm("`part'", "o\.") local base = 0
        }
        local flags `flags' `base'
        if !`base' {
            local kept_names `kept_names' `name'
            if `has_vars' local kept_vars `kept_vars' `var'
        }
    }
    return local base `flags'
    return local names `kept_names'
    return local vars `kept_vars'
end
