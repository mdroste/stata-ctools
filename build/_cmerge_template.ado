*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_template
    args file newvars mergevars
    local keepvars : list newvars - mergevars
    if "`keepvars'" != "" keep `keepvars'
    if _N > 0 {
        quietly keep in 1
        quietly drop in 1
    }
    foreach v of local mergevars {
        quietly generate byte `v' = .
    }
    * keep, unlike drop _all, retains the _dta notes and characteristics
    * that append passes on to the master
    if "`keepvars'" == "" {
        if "`mergevars'" != "" keep `mergevars'
        else drop _all
    }
    quietly save `"`file'"', emptyok replace
end
