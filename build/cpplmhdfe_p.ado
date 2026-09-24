*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools
program define cpplmhdfe_p
    version 14.1
    if "`e(cmd)'" != "cpplmhdfe" error 301
    syntax newvarname [if] [in] [, XB]
    if "`xb'" == "" {
        di as error "cpplmhdfe: specify xb for the slope index; fitted means are not available"
        di as error "xb excludes absorbed effects and offset/exposure contributions"
        exit 198
    }
    _predict `typlist' `varlist' `if' `in', xb
    label variable `varlist' "Slope index (excludes absorbed effects and offsets)"
end
