*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! ctools: i.f indicators built by the plugin instead of fvrevar and fvexpand
* _ctools_fvnative <terms>, touse(varname) plugin(string) [threads(#)]
*
* terms is an independent-variable list as syntax varlist(fv ts) leaves it.
* When every term is a numeric variable or one of i.v, ib#.v, ibn.v,
* ib(first).v and ib(last).v on a numeric v without fvset base (and set
* fvbase is on), and at least one term is such a factor term, the levels of
* the factor variables on the touse rows come from one plugin pass and the
* caller's locals get
*     __fvn_ok      1 (otherwise 0, and nothing else is set)
*     __fvn_names   the names fvexpand <terms> if touse gives
*     __fvn_vars    one variable per name: the term, or its factor variable
*     __fvn_levels  one value per name: its level, or -1 for a plain term
*     __fvn_counts  one value per term: how many names it gave
* Anything else (interactions, c./time-series operators, i(numlist).v,
* ib(freq).v, noninteger or negative values, levels above 65535) leaves the
* caller on Stata's own fvrevar/fvexpand path, which also gives Stata's
* messages. At 10M rows fvrevar and fvexpand of one i.f cost about 0.2 s
* each; the plugin pass reads each factor variable once.
program define _ctools_fvnative
    version 14.1
    syntax anything(name=terms), TOUSE(varname numeric) PLUGIN(string) [THREADS(integer 0)]
    c_local __fvn_ok 0
    if "`c(fvbase)'" == "off" exit
    local kinds
    local tvars
    local fvars
    foreach t of local terms {
        if ustrregexm("`t'", "^i(b[0-9]+|bn|b\(first\)|b\(last\))?\.([A-Za-z_][A-Za-z0-9_]*)$") {
            local spec = ustrregexs(1)
            local v = ustrregexs(2)
            capture confirm numeric variable `v', exact
            if _rc exit
            if `"`: char `v'[_fv_base]'"' != "" exit
            if "`spec'" == "" | "`spec'" == "b(first)" local kind first
            else if "`spec'" == "b(last)" local kind last
            else if "`spec'" == "bn" local kind bn
            else local kind = substr("`spec'", 2, .)
            local kinds `kinds' `kind'
            local tvars `tvars' `v'
            local fvars `fvars' `v'
        }
        else {
            capture confirm numeric variable `t', exact
            if _rc exit
            local kinds `kinds' plain
            local tvars `tvars' `t'
        }
    }
    local fvars : list uniq fvars
    if "`fvars'" == "" exit
    capture program ctools_plugin, plugin using("`plugin'")
    if _rc != 0 & _rc != 110 exit
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    local nf : word count `fvars'
    forvalues j = 1/`nf' {
        local ctools_fvl_`j'
        local ctools_fvl_ok_`j' 0
    }
    capture plugin call ctools_plugin `fvars' `touse', "fvlevels `threadopt'"
    if _rc exit
    forvalues j = 1/`nf' {
        if "`ctools_fvl_ok_`j''" != "1" exit
    }
    local names
    local vars
    local levels
    local i = 0
    foreach t of local terms {
        local ++i
        local kind : word `i' of `kinds'
        if "`kind'" == "plain" {
            local names `names' `t'
            local vars `vars' `t'
            local levels `levels' -1
            local counts `counts' 1
            continue
        }
        local v : word `i' of `tvars'
        local j : list posof "`v'" in fvars
        local all `ctools_fvl_`j''
        local nl : word count `all'
        local base
        if "`kind'" == "first" local base : word 1 of `all'
        else if "`kind'" == "last" local base : word `nl' of `all'
        else if "`kind'" != "bn" {
            * ib#.v: a base that is not a level is still listed, in order
            local base = real("`kind'")
            if missing(`base') | `base' < 0 | `base' > 65535 | `base' != floor(`base') exit
            local base = strofreal(`base', "%12.0f")
            if !`: list base in all' {
                local merged
                local placed = 0
                foreach l of local all {
                    if !`placed' & `l' > `base' {
                        local merged `merged' `base'
                        local placed = 1
                    }
                    local merged `merged' `l'
                }
                if !`placed' local merged `merged' `base'
                local all `merged'
            }
        }
        local first = 1
        foreach l of local all {
            if "`kind'" == "bn" local name = cond(`first', "`l'bn.`v'", "`l'.`v'")
            else local name = cond("`l'" == "`base'", "`l'b.`v'", "`l'.`v'")
            local first = 0
            local names `names' `name'
            local vars `vars' `v'
            local levels `levels' `l'
        }
        local counts `counts' `: word count `all''
    }
    c_local __fvn_names `names'
    c_local __fvn_vars `vars'
    c_local __fvn_levels `levels'
    c_local __fvn_counts `counts'
    c_local __fvn_ok 1
end
