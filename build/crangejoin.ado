*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
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
    * The plugin cannot write strL values: strL variables travel through the
    * join as observation numbers (mobs, uobs) and are copied in afterwards.
    * Break is recognized only inside _crangejoin_impl.
    tempvar mobs uobs
    nobreak {
        macro drop CTOOLS_crangejoin_plugin CTOOLS_crangejoin_mvars
        global CTOOLS_crangejoin_obs `mobs' `uobs'
        if "`mode'" == "" _ctools_undo record CTOOLS_u_crangejoin_
        capture noisily break _crangejoin_impl `0'
        local rc = _rc
        if `rc' {
            if "$CTOOLS_crangejoin_plugin" != "" {
                capture program ctools_plugin, plugin using("$CTOOLS_crangejoin_plugin")
                if "`mode'" == "" {
                    * Observations first, while the master variables lead;
                    * the variables the plugin loaded, when strL ones stayed out
                    local mastervars $CTOOLS_crangejoin_mvars
                    if "`mastervars'" == "" {
                        local k0 = $CTOOLS_u_crangejoin_k
                        mata: st_local("mastervars", invtokens(st_varname(1..`k0')))
                    }
                    capture noisily plugin call ctools_plugin `mastervars', "crangejoin rollback"
                    if _rc di as error "crangejoin: the master observations could not be restored"
                }
                capture plugin call ctools_plugin, "crangejoin clear"
            }
            if "`mode'" == "" _ctools_undo rollback CTOOLS_u_crangejoin_
            macro drop CTOOLS_crangejoin_plugin CTOOLS_crangejoin_mvars CTOOLS_crangejoin_obs
            capture mata: mata drop __crangejoin_U
            exit `rc'
        }
        * Free the master copy the plugin kept
        if "$CTOOLS_crangejoin_plugin" != "" {
            capture program ctools_plugin, plugin using("$CTOOLS_crangejoin_plugin")
            capture plugin call ctools_plugin, "crangejoin clear"
        }
        if "`mode'" == "" _ctools_undo clear CTOOLS_u_crangejoin_
        macro drop CTOOLS_crangejoin_plugin CTOOLS_crangejoin_mvars CTOOLS_crangejoin_obs
        capture mata: mata drop __crangejoin_U
        * rangejoin replaces the dataset, which empties e(sample) and keeps
        * the other estimation results; the plugin rewrote the rows in place,
        * so the hidden e(sample) column would mark unrelated observations.
        local esample : e(functions)
        if `: list posof "sample" in esample' {
            tempname eshold
            tempvar esvar
            _estimates hold `eshold', copy varname(`esvar')
            drop `esvar'
            _estimates unhold `eshold'
        }
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
    * The plugin writes fixed strings only. Master strL variables (other than
    * by() keys, which rangejoin refuses) stay out of it: the plugin carries
    * each observation number (mobs), and the values are copied after the join.
    local obsvars $CTOOLS_crangejoin_obs
    if "`obsvars'" == "" {
        tempvar mobs uobs
        local obsvars `mobs' `uobs'
    }
    gettoken mobs uobs : obsvars
    local mstrl
    foreach v of local mastervars {
        if "`: type `v''" == "strL" & !`: list v in by' local mstrl `mstrl' `v'
    }
    local mplug : list mastervars - mstrl
    if "`mstrl'" != "" local mplug `mplug' `mobs'
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')
    local timing = "`verbose'" != ""
    local nby : word count `by'
    local byindices
    local mbytypes
    foreach v of local by {
        local pos : list posof "`v'" in mplug
        local byindices `byindices' `pos'
        local mbytypes `mbytypes' `: type `v''
    }
    local usingopts `"key(`key') by(`by') keepusing(`keepusing') prefix(`prefix') suffix(`suffix') `all'"'
    local usingopts `"`usingopts' mastervars(`mastervars') threadopt(`threadopt') timing(`timing') mbytypes(`mbytypes') uobs(`uobs')"'
    local masteropts `"mastervars(`mplug') allmaster(`mastervars') mstrl(`mstrl') mobs(`mobs') uobs(`uobs')"'
    local nmaster : word count `mplug'
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
        _crangejoin_master `"`schema'"', `masteropts' lo(`lo') hi(`hi') ///
            renamed(`renamed') ustrl(`ustrl') nmaster(`nmaster') byindices(`byindices') ///
            threadopt(`threadopt') by(`by') usingbytypes(`usingbytypes') `fast'
        exit
    }
    * Option snapshot, and every join before Stata 16: read the using data in
    * place of the master, from a snapshot of the whole dataset.
    preserve
    capture noisily {
        _crangejoin_using `"`using'"' `"`schema'"', `usingopts'
        restore, preserve
        _crangejoin_master `"`schema'"', `masteropts' lo(`lo') hi(`hi') ///
            renamed(`renamed') ustrl(`ustrl') nmaster(`nmaster') byindices(`byindices') ///
            threadopt(`threadopt') by(`by') usingbytypes(`usingbytypes') norestore
    }
    local rc = _rc
    capture plugin call ctools_plugin, "crangejoin clear"
    if `rc' {
        restore
        exit `rc'
    }
    restore, not
end

* Read the using data into the current dataset, rename clashing variables,
* check the by() types against the master's, index the data in the plugin and
* save their empty schema. Using strL variables stay out of the plugin: their
* values are kept in Mata (__crangejoin_U), and uobs, the observation number
* the plugin carries, finds them. Sets locals renamed and ustrl (the strL
* variables among them) in the caller.
program define _crangejoin_using
    version 14.1
    syntax anything(name=files), key(name) uobs(name) [by(string) keepusing(string) ///
        prefix(string) suffix(string) all mastervars(string) threadopt(string) ///
        timing(integer 0) mbytypes(string)]
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
    local keyerror = 0
    foreach v of local by {
        local ++i
        confirm variable `v', exact
        local ut : type `v'
        local mt : word `i' of `mbytypes'
        * Keys merge refuses: string with numeric, or strL on either side
        if (substr("`ut'",1,3) == "str") != (substr("`mt'",1,3) == "str") | ///
            "`ut'" == "strL" | "`mt'" == "strL" local keyerror = 1
        local usingbytypes `usingbytypes' `ut'
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
        if "`v'" == "`key'" local keyname `dest'
        local renamed `renamed' `dest'
        if "`: type `dest''" == "strL" local ustrl `ustrl' `dest'
    }
    * rangejoin pairs the data with merge, which then reports these keys
    if `keyerror' _crangejoin_keyerror, by(`by') mbytypes(`mbytypes')
    local uplug : list renamed - ustrl
    if "`ustrl'" != "" {
        quietly generate double `uobs' = _n
        local uplug `uplug' `uobs'
        mata: __crangejoin_U = st_sdata(., tokens(st_local("ustrl")))
    }
    local nby : word count `by'
    local keypos : list posof "`keyname'" in uplug
    local keypos = `keypos' + `nby'
    _ctools_strw `by' `uplug'
    plugin call ctools_plugin `by' `uplug', ///
        "crangejoin `threadopt' using `nby' `keypos' `timing'"
    keep `renamed' `: list uplug - renamed'
    quietly keep if 0
    quietly save `"`schema'"'
    c_local renamed `renamed'
    c_local ustrl `ustrl'
    c_local usingbytypes `usingbytypes'
end

* merge's messages and error code for by() keys it refuses (strL, or a string
* key paired with a numeric one), from empty data with the master's and the
* using data's key types. Replaces the data in memory.
program define _crangejoin_keyerror
    version 14.1
    syntax, by(string) mbytypes(string)
    tempvar obs
    tempfile keys
    keep `by'
    quietly keep if 0
    quietly generate long `obs' = .
    quietly save `"`keys'"'
    clear
    local i = 0
    foreach v of local by {
        local ++i
        local t : word `i' of `mbytypes'
        if substr("`t'", 1, 3) == "str" quietly generate `t' `v' = ""
        else quietly generate `t' `v' = .
    }
    quietly generate long `obs' = .
    capture noisily quietly merge 1:1 `by' `obs' using `"`keys'"', nogenerate
    local rc = _rc
    if !`rc' local rc = 106
    exit `rc'
end

* Append the using schema to the master and have the plugin write the join.
* mastervars are the master variables the plugin takes (allmaster without the
* strL ones mstrl, plus mobs when there are any); renamed are the using
* variables, ustrl the strL ones among them, which uobs finds.
program define _crangejoin_master
    version 14.1
    syntax anything(name=schema), mastervars(string) allmaster(string) lo(name) hi(name) ///
        renamed(string) nmaster(integer) mobs(name) uobs(name) [mstrl(string) ///
        ustrl(string) byindices(string) threadopt(string) by(string) ///
        usingbytypes(string) fast norestore]
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    if "`mstrl'" != "" {
        quietly generate double `mobs' = _n
        * crangejoin's rollback hands the plugin the variables it loaded
        global CTOOLS_crangejoin_mvars `mastervars'
    }
    quietly append using `schema'
    local uplug : list renamed - ustrl
    if "`ustrl'" != "" local uplug `uplug' `uobs'
    * The master variables keep the first places, where a rollback expects them
    order `allmaster' `lo' `hi' `uplug'
    _ctools_strw `mastervars' `lo' `hi' `uplug'
    plugin call ctools_plugin `mastervars' `lo' `hi' `uplug', ///
        "crangejoin `threadopt' prepare `nmaster' `byindices'"
    quietly set obs `crangejoin_nout'
    * Validation hooks: fail as Break would before the plugin writes, or have
    * it fail after writing half the output rows, or fail once it wrote them
    if "$CTOOLS_crangejoin_fault" == "schema" error 1
    local wopts
    if "$CTOOLS_crangejoin_fault" == "store" local wopts fault_store
    if "`fast'`norestore'" != "" local wopts `wopts' norestore
    plugin call ctools_plugin `mastervars' `lo' `hi' `uplug', ///
        "crangejoin `threadopt' write `wopts'"
    if "$CTOOLS_crangejoin_fault" == "post" error 1
    drop `lo' `hi'
    * Using strL values by the using observation each row joined (none: "")
    if "`ustrl'" != "" mata: _crangejoin_fill(st_local("uobs"), tokens(st_local("ustrl")))
    * Master strL values by each row's master observation, into copies: the
    * master variables stay as they were until the end, for a rollback
    local copies
    foreach v of local mstrl {
        tempvar copy
        quietly generate strL `copy' = `v'[`mobs']
        local copies `copies' `copy'
    }
    * Validation hook: fail as Break would once the strL values are in place
    if "$CTOOLS_crangejoin_fault" == "strl" error 1
    * rangejoin's by() variables take the wider of the master and using
    * storage types (the values are the master's)
    local numeric_types "byte int long float double"
    local i = 0
    foreach v of local by {
        local ++i
        local mt : type `v'
        local ut : word `i' of `usingbytypes'
        local target `mt'
        if substr("`mt'", 1, 3) == "str" {
            local width = max(real(substr("`mt'", 4, .)), real(substr("`ut'", 4, .)))
            if !missing(`width') local target str`width'
        }
        else {
            local mr : list posof "`mt'" in numeric_types
            local ur : list posof "`ut'" in numeric_types
            local target : word `=max(`mr', `ur')' of `numeric_types'
            if (`mr' == 3 & `ur' == 4) | (`mr' == 4 & `ur' == 3) local target double
        }
        if "`target'" != "`mt'" quietly recast `target' `v'
    }
    if "`mstrl'`ustrl'" != "" {
        nobreak {
            * Each copy takes its variable's place: format, label and
            * characteristics (notes included)
            foreach v of local mstrl {
                gettoken copy copies : copies
                format `copy' `: format `v''
                mata: st_varlabel(st_local("copy"), st_varlabel(st_local("v")))
                local chars : char `v'[]
                foreach c of local chars {
                    mata: st_global(st_local("copy") + "[" + st_local("c") + "]", ///
                        st_global(st_local("v") + "[" + st_local("c") + "]"))
                }
                drop `v'
                rename `copy' `v'
            }
            capture drop `mobs'
            capture drop `uobs'
            order `allmaster' `renamed'
        }
    }
end

capture mata: mata drop _crangejoin_fill()
mata:
/* Store the using strL values (__crangejoin_U, one column per variable of
   names) in each joined row, by the using observation uobs holds; rows
   without a match (uobs missing) keep "". In blocks of rows, to bound the
   strings held at once. */
void _crangejoin_fill(string scalar uobs, string rowvector names)
{
    external string matrix __crangejoin_U
    real scalar n, a, b
    real colvector idx, sel
    string matrix values

    n = st_nobs()
    for (a = 1; a <= n; a = b + 1) {
        b = min((a + 99999, n))
        idx = st_data((a, b), uobs)
        values = J(rows(idx), cols(names), "")
        sel = selectindex(idx :< .)
        if (rows(sel)) values[sel, .] = __crangejoin_U[idx[sel], .]
        st_sstore((a, b), names, values)
    }
}
end
