*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! C-accelerated interval joins (ctools)
program define crangejoin
    version 14.1
    syntax anything(name=interval) using/, [BY(varlist) Keepusing(string) ///
        Prefix(string) Suffix(string) All Verbose THReads(integer 0) SNAPshot FAST]
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
    local mode `snapshot'`fast'
    * A failed or interrupted join leaves the master as it was without a
    * snapshot of the dataset (Stata 16 or later): the using data are read in
    * a frame of their own; the master only gains variables, value labels,
    * _dta notes and characteristics (append) and observations at the end;
    * and the plugin keeps the master observations it loaded, writing them
    * back after a failure. _ctools_undo removes the rest and restores the
    * sort marker and c(changed). Option snapshot preserves the whole dataset
    * instead (on disk above c(max_preservemem)), as does every join before
    * Stata 16, where the using data can only be read in place of the master.
    * Option fast keeps no copy and undoes nothing.
    * Break is recognized only inside _crangejoin_impl.
    nobreak {
        macro drop CTOOLS_crangejoin_plugin
        if "`mode'" == "" _ctools_undo record CTOOLS_u_crangejoin_
        capture noisily break _crangejoin_impl `0'
        local rc = _rc
        if `rc' {
            if "$CTOOLS_crangejoin_plugin" != "" {
                capture program ctools_plugin, plugin using("$CTOOLS_crangejoin_plugin")
                if "`mode'" == "" {
                    * Observations first, while the master variables lead
                    local k0 = $CTOOLS_u_crangejoin_k
                    mata: st_local("mastervars", invtokens(st_varname(1..`k0')))
                    capture noisily plugin call ctools_plugin `mastervars', "crangejoin rollback"
                    if _rc di as error "crangejoin: the master observations could not be restored"
                }
                capture plugin call ctools_plugin, "crangejoin clear"
            }
            if "`mode'" == "" _ctools_undo rollback CTOOLS_u_crangejoin_
            macro drop CTOOLS_crangejoin_plugin
            exit `rc'
        }
        * Free the master copy the plugin kept
        if "$CTOOLS_crangejoin_plugin" != "" {
            capture program ctools_plugin, plugin using("$CTOOLS_crangejoin_plugin")
            capture plugin call ctools_plugin, "crangejoin clear"
        }
        if "`mode'" == "" _ctools_undo clear CTOOLS_u_crangejoin_
        macro drop CTOOLS_crangejoin_plugin
    }
end

program define _crangejoin_impl
    version 14.1
    syntax anything(name=interval) using/, [BY(varlist) Keepusing(string) ///
        Prefix(string) Suffix(string) All Verbose THReads(integer 0) SNAPshot FAST]
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
    local usingopts `"key(`key') by(`by') keepusing(`keepusing') prefix(`prefix') suffix(`suffix') `all'"'
    local bytypes
    foreach v of local by {
        local pos : list posof "`v'" in mastervars
        local bytypes `bytypes' `bytype_`pos''
    }
    local usingopts `"`usingopts' mastervars(`mastervars') threadopt(`threadopt') timing(`timing') bytypes(`bytypes')"'
    tempfile schema
    if c(stata_version) >= 16 & "`snapshot'" == "" {
        * The using data go to a frame of their own: the plugin indexes them,
        * and only their empty schema is appended to the master.
        tempname uframe
        frame create `uframe'
        * Set before the plugin indexes the using data, so that the wrapper
        * releases them when a later step of _crangejoin_using fails
        global CTOOLS_crangejoin_plugin "`__ctools_plugin'"
        frame `uframe': _crangejoin_using `"`using'"' `"`schema'"', `usingopts'
        frame drop `uframe'
        _crangejoin_master `"`schema'"', mastervars(`mastervars') lo(`lo') hi(`hi') ///
            renamed(`renamed') nmaster(`nmaster') byindices(`byindices') ///
            threadopt(`threadopt') `fast'
        exit
    }
    * Option snapshot, and every join before Stata 16: read the using data in
    * place of the master, from a snapshot of the whole dataset.
    preserve
    capture noisily {
        _crangejoin_using `"`using'"' `"`schema'"', `usingopts'
        restore, preserve
        _crangejoin_master `"`schema'"', mastervars(`mastervars') lo(`lo') hi(`hi') ///
            renamed(`renamed') nmaster(`nmaster') byindices(`byindices') ///
            threadopt(`threadopt') norestore
    }
    local rc = _rc
    capture plugin call ctools_plugin, "crangejoin clear"
    if `rc' {
        restore
        exit `rc'
    }
    restore, not
end

* Read the using data into the current dataset, check the by() types against
* the master's, rename clashing variables, index them in the plugin and save
* their empty schema. Sets local renamed in the caller.
program define _crangejoin_using
    version 14.1
    syntax anything(name=files), key(name) [by(string) keepusing(string) ///
        prefix(string) suffix(string) all mastervars(string) threadopt(string) ///
        timing(integer 0) bytypes(string)]
    gettoken using files : files
    gettoken schema : files
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    if `"`keepusing'"' != "" quietly use `keepusing' `by' `key' using `"`using'"', clear
    else quietly use `"`using'"', clear
    if _N > 2147483647 exit 920
    confirm numeric variable `key', exact
    local i = 0
    foreach v of local by {
        local ++i
        confirm variable `v', exact
        if (substr("`: type `v''",1,3) == "str") != `: word `i' of `bytypes'' exit 106
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
    local nby : word count `by'
    local keypos : list posof "`keyname'" in renamed
    local keypos = `keypos' + `nby'
    _ctools_strw `by' `renamed'
    plugin call ctools_plugin `by' `renamed', ///
        "crangejoin `threadopt' using `nby' `keypos' `timing'"
    keep `renamed'
    quietly keep if 0
    quietly save `"`schema'"'
    c_local renamed `renamed'
end

* Append the using schema to the master and have the plugin write the join.
program define _crangejoin_master
    version 14.1
    syntax anything(name=schema), mastervars(string) lo(name) hi(name) ///
        renamed(string) nmaster(integer) [byindices(string) threadopt(string) ///
        fast norestore]
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    quietly append using `schema'
    order `mastervars' `lo' `hi' `renamed'
    _ctools_strw `mastervars' `lo' `hi' `renamed'
    plugin call ctools_plugin `mastervars' `lo' `hi' `renamed', ///
        "crangejoin `threadopt' prepare `nmaster' `byindices'"
    quietly set obs `crangejoin_nout'
    * Validation hooks: fail as Break would before the plugin writes, or have
    * it fail after writing half the output rows, or fail once it wrote them
    if "$CTOOLS_crangejoin_fault" == "schema" error 1
    local wopts
    if "$CTOOLS_crangejoin_fault" == "store" local wopts fault_store
    if "`fast'`norestore'" != "" local wopts `wopts' norestore
    plugin call ctools_plugin `mastervars' `lo' `hi' `renamed', ///
        "crangejoin `threadopt' write `wopts'"
    if "$CTOOLS_crangejoin_fault" == "post" error 1
    drop `lo' `hi'
end
