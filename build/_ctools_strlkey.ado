*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! _ctools_strlkey — numeric stand-ins for strL keys the plugin cannot read
*!
*! Usage:  _ctools_strlkey varlist, generate(namelist) [missing]
*!         _ctools_strlkey varlist, check
*!
*! The plugin reads strL text of at most 2,045 bytes and no binary values.
*! Each strL variable of varlist that holds another value anywhere in the
*! data gets a stand-in: the corresponding name of generate() becomes a long
*! holding the rank of the value among the variable's distinct values, in
*! Stata's sort order. Grouping by the ranks, ordering by them and counting
*! their distinct values give what the strL values give. With missing, ""
*! (a missing string) gets a missing rank, so that commands that leave out
*! missing keys still leave it out; otherwise "" ranks first, as it sorts.
*!
*! Sets local __ctools_keys in the caller: varlist with each such variable
*! replaced by its stand-in. Variables the plugin can read stay as they are.
*! The order of the observations, the sort marker and r() are left alone.
*! check only sets local __ctools_strl_bad in the caller: the variables of
*! varlist that would get a stand-in.

program define _ctools_strlkey
    version 14.1
    syntax varlist, [GENerate(namelist) MISSing CHECK]
    if "`check'" == "" & `: word count `generate'' != `: word count `varlist'' exit 198

    local strl
    foreach v of local varlist {
        if "`: type `v''" == "strL" local strl `strl' `v'
    }
    local strl : list uniq strl
    if "`strl'" == "" | _N == 0 {
        c_local __ctools_keys `varlist'
        c_local __ctools_strl_bad
        exit
    }

    * strlcheck reads lengths and binary flags only. If the plugin cannot
    * answer, every strL key gets a stand-in, which is exact anyway.
    local bad `strl'
    capture _ctools_load
    if !_rc {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if !_rc | _rc == 110 {
            local ctools_strl_bad ""
            capture plugin call ctools_plugin `strl', "strlcheck"
            if !_rc & `: word count `ctools_strl_bad'' == `: word count `strl'' {
                local bad
                local i = 0
                foreach v of local strl {
                    local ++i
                    if `: word `i' of `ctools_strl_bad'' local bad `bad' `v'
                }
            }
        }
    }

    if "`check'" != "" {
        c_local __ctools_strl_bad `bad'
        exit
    }

    local keys
    local made
    local i = 0
    foreach v of local varlist {
        local ++i
        if !`: list v in bad' {
            local keys `keys' `v'
            continue
        }
        * A variable named twice shares one stand-in
        local at : list posof "`v'" in made
        if `at' {
            local keys `keys' `: word `at' of `madeid''
            continue
        }
        local id : word `i' of `generate'
        _ctools_strlkey_rank `v' `id' `missing'
        local made `made' `v'
        local madeid `madeid' `id'
        local keys `keys' `id'
    }
    c_local __ctools_keys `keys'
end

* Rank the values of v in Stata's sort order (sortpreserve puts the
* observations, and the sort marker, back; on Break as well).
program define _ctools_strlkey_rank, sortpreserve
    args v id missing
    sort `v'
    quietly by `v': generate long `id' = _n == 1
    quietly replace `id' = sum(`id')
    if "`missing'" != "" quietly replace `id' = . if `v' == ""
end
