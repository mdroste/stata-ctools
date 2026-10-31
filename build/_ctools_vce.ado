*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! ctools: classify a vce() specification
* _ctools_vce <kind> <rest> : <specification> sets two locals of the caller:
* <kind> is cluster (any word starting with cl), robust (robust, r),
* unadjusted, ols, iid, bootstrap (any word starting with boot), or empty for
* anything else, case-insensitively; <rest> is the text after the first word.
* Which kinds a command accepts, and its messages, stay with the command.
program define _ctools_vce
    version 14.1
    gettoken kind_macro 0 : 0
    gettoken rest_macro 0 : 0
    gettoken colon 0 : 0, parse(":")
    local lower = lower(trim(`"`0'"'))
    local kind
    if substr(`"`lower'"', 1, 2) == "cl" local kind cluster
    else if inlist(`"`lower'"', "robust", "r") local kind robust
    else if inlist(`"`lower'"', "unadjusted", "ols", "iid") local kind `lower'
    else if substr(`"`lower'"', 1, 4) == "boot" local kind bootstrap
    gettoken first rest : 0
    local rest = trim(`"`rest'"')
    c_local `kind_macro' `kind'
    c_local `rest_macro' `rest'
end
