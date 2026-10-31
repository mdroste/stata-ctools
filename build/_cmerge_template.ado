*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_template
    args file newvars mergevars
    local keepvars : list newvars - mergevars
    if "`keepvars'" != "" {
        keep `keepvars'
        if _N > 0 {
            quietly keep in 1
            quietly drop in 1
        }
    }
    else drop _all
    foreach v of local mergevars {
        quietly generate byte `v' = .
    }
    quietly save `"`file'"', emptyok replace
end
