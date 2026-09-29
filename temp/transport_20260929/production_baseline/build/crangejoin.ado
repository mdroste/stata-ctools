*! C-accelerated interval joins (ctools)
program define crangejoin
    version 14.1
    syntax anything(name=interval) using/, [BY(varlist) Keepusing(string) ///
        Prefix(string) Suffix(string) All Verbose THReads(integer 0)]
    if `threads' < 0 exit 198
    if _N > 2147483647 exit 920
    if _N == 0 error 2000
    local interval : subinstr local interval "," " ", all
    if `: word count `interval'' != 3 exit 198
    tokenize `interval'
    local key `1'
    local low `2'
    local high `3'
    confirm name `key'
    capture confirm numeric variable `key', exact
    if _rc == 7 exit 7
    local haskey = _rc == 0
    if `: list key in by' {
        di as error "key variable cannot appear in by()"
        exit 198
    }
    if "`prefix'`suffix'" == "" local suffix _U
    unab mastervars : *
    local nmaster : word count `mastervars'
    tempvar lo hi
    foreach side in low high {
        local dest `lo'
        if "`side'" == "high" local dest `hi'
        capture confirm numeric variable ``side''
        if !_rc quietly gen double `dest' = cond(missing(``side''), ., ``side'')
        else {
            capture confirm number ``side''
            if _rc & "``side''" != "." exit 198
            if `haskey' quietly gen double `dest' = `key' + ``side''
            else quietly gen double `dest' = ``side''
        }
    }
    quietly replace `lo' = c(mindouble) if missing(`lo')
    if `haskey' {
        quietly replace `lo' = 1 if missing(`key')
        quietly replace `hi' = 0 if missing(`key')
    }
    quietly count if `lo' <= `hi'
    if !r(N) {
        di as error "no observation with valid interval bounds to use"
        exit 2000
    }
    * Fixed strings are the writable string type supported by the plugin API.
    foreach v of local mastervars {
        if "`: type `v''" == "strL" {
            di as error "crangejoin: recast strL variables to fixed strings when lossless"
            exit 109
        }
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')
    local timing = "`verbose'" != ""
    local nby : word count `by'
    local byindices
    foreach v of local by {
        local pos : list posof "`v'" in mastervars
        local byindices `byindices' `pos'
        local bytype_`pos' = substr("`: type `v''",1,3) == "str"
    }
    tempfile schema
    preserve
    capture noisily {
        if `"`keepusing'"' != "" quietly use `keepusing' `by' `key' using `"`using'"', clear
        else quietly use `"`using'"', clear
        if _N > 2147483647 exit 920
        confirm numeric variable `key', exact
        foreach v of local by {
            confirm variable `v', exact
            local pos : list posof "`v'" in mastervars
            if (substr("`: type `v''",1,3) == "str") != `bytype_`pos'' exit 106
        }
        quietly drop if missing(`key')
        if _N == 0 error 2000
        unab uvars : *
        local uvars : list uvars - by
        local renamed
        foreach v of local uvars {
            local dest `v'
            if "`all'" != "" | `: list v in mastervars' {
                local dest `prefix'`v'`suffix'
                confirm name `dest'
                if `: list dest in mastervars' exit 110
                rename `v' `dest'
            }
            if "`: type `dest''" == "strL" {
                di as error "crangejoin: recast using strL variables to fixed strings when lossless"
                exit 109
            }
            if "`v'" == "`key'" local keyname `dest'
            local renamed `renamed' `dest'
        }
        local keypos : list posof "`keyname'" in renamed
        local keypos = `keypos' + `nby'
        _ctools_strw `by' `renamed'
        plugin call ctools_plugin `by' `renamed', ///
            "crangejoin `threadopt' using `nby' `keypos' `timing'"
        keep `renamed'
        quietly keep if 0
        quietly save `"`schema'"'
        restore, preserve
        quietly append using `"`schema'"'
        order `mastervars' `lo' `hi' `renamed'
        _ctools_strw `mastervars' `lo' `hi' `renamed'
        plugin call ctools_plugin `mastervars' `lo' `hi' `renamed', ///
            "crangejoin `threadopt' prepare `nmaster' `byindices'"
        quietly set obs `crangejoin_nout'
        plugin call ctools_plugin `mastervars' `lo' `hi' `renamed', ///
            "crangejoin `threadopt' write"
        drop `lo' `hi'
    }
    local rc = _rc
    capture plugin call ctools_plugin, "crangejoin clear"
    if `rc' {
        restore
        exit `rc'
    }
    restore, not
end
