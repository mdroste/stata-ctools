*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_keys_sorted, rclass
    syntax varlist
    local n : word count `varlist'
    local v : word `n' of `varlist'
    local cond "(`v'[_n-1] <= `v')"
    forvalues j = `=`n'-1'(-1)1 {
        local v : word `j' of `varlist'
        local cond "(`v'[_n-1] < `v' | (`v'[_n-1] == `v' & `cond'))"
    }
    capture assert `cond' in 2/l
    return scalar sorted = (_rc == 0)
end
