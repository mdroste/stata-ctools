*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! C-accelerated interval joins (ctools)
*! original: rangejoin 1.1.3 (April 13, 2021) by Robert Picard
*! original: see help crangejoin for details
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
        capture noisily break _crangejoin_impl `macval(0)'
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
    * Each bound is a variable, or a number added to the key. A missing low
    * bound is the lowest double, a missing high bound stays missing (the
    * highest), and an observation with a missing key gets the empty
    * interval [1, 0]. One pass builds each bound.
    tempvar lo hi
    foreach side in low high {
        capture confirm numeric variable ``side''
        if !_rc local `side'val ``side''
        else {
            capture confirm number ``side''
            if _rc & "``side''" != "." exit 198
            if `haskey' local `side'val (`key' + ``side'')
            else local `side'val (``side'')
        }
    }
    * (missing() may not be nested: each value is tested once)
    local lowexp cond(missing(`lowval'), c(mindouble), `lowval')
    local highexp cond(missing(`highval'), ., `highval')
    if `haskey' {
        local lowexp cond(missing(`key'), 1, `lowexp')
        local highexp cond(missing(`key'), 0, `highexp')
    }
    quietly gen double `lo' = `lowexp'
    quietly gen double `hi' = `highexp'
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
        frame `uframe': _crangejoin_using `"`macval(using)'"' `"`schema'"', `usingopts'
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
        _crangejoin_using `"`macval(using)'"' `"`schema'"', `usingopts'
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
    * The varlist-using form of use rescans macro characters in its filename.
    * Preserve selective loading for ordinary paths; literal macro characters
    * require loading by filename first, then selecting the columns.
    local literal = strpos(`"`macval(using)'"', char(36)) | strpos(`"`macval(using)'"', char(96))
    if `"`keepusing'"' != "" & !`literal' {
        quietly use `keepusing' `by' `key' using `"`macval(using)'"', clear
    }
    else {
        quietly use `"`macval(using)'"', clear
        if `"`keepusing'"' != "" keep `keepusing' `by' `key'
    }
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
    * Only the schema of the using data is saved: drop the observations
    * before the variables, by range and not by expression
    if _N quietly drop in 1/l
    keep `renamed' `: list uplug - renamed'
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

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* Shared plugin loader and identity contract
program define _ctools_load
    version 14.1
    * Not rclass, and findfile's r() is undone, so the caller's r() survives;
    * identity comes back through c_local.
    local __os = c(os)
    local __machine = c(machine_type)
    local __is_mac = 0
    if "`__os'" == "MacOSX" {
        local __is_mac = 1
    }
    else if strpos(lower("`__machine'"), "mac") > 0 {
        local __is_mac = 1
    }

    local __plugin = ""
    if "`__os'" == "Windows" {
        local __plugin "ctools_windows.plugin"
    }
    else if `__is_mac' {
        local __is_arm = 0
        if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
            local __is_arm = 1
        }
        if `__is_arm' == 0 & !strpos(lower("`__machine'"), "intel") & !strpos(lower("`__machine'"), "x86") {
            tempfile __archfile
            quietly shell uname -m > "`__archfile'" 2>&1
            tempname __fh
            file open `__fh' using "`__archfile'", read text
            file read `__fh' __archline
            file close `__fh'
            capture erase "`__archfile'"
            if strpos("`__archline'", "arm64") > 0 {
                local __is_arm = 1
            }
        }
        if `__is_arm' {
            local __plugin "ctools_mac_arm.plugin"
        }
        else {
            local __plugin "ctools_mac_x86.plugin"
        }
    }
    else if "`__os'" == "Unix" {
        local __plugin "ctools_linux.plugin"
    }
    else {
        local __plugin "ctools.plugin"
    }

    * The ado/plugin argument protocol, stamped from config/release.json.
    local __ado_protocol "4"
    local __reinstall `"net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace"'

    tempname __rhold
    _return hold `__rhold'
    capture findfile `__plugin'
    local __found = (_rc == 0)
    local __path `"`r(fn)'"'
    _return restore `__rhold'
    local __wanted "`__plugin'"
    if !`__found' & "`__plugin'" != "ctools.plugin" local __plugin "ctools.plugin"
    capture program ctools_plugin, plugin using("`__plugin'")
    local __load_rc = _rc
    if `__load_rc' != 0 & `__load_rc' != 110 {
        di as error "ctools: could not load the plugin `__plugin'"
        di as error `"  Stata reports c(os) "`__os'" and c(machine_type) "`__machine'"; ctools selected `__wanted'."'
        if `__found' di as error `"  The plugin file is `__path'."'
        else di as error "  `__wanted' was not found on the ado-path; reinstall ctools."
        di as error "  Stata's plugin loader reported:"
        capture noisily program ctools_plugin, plugin using("`__plugin'")
        if "`__os'" == "Unix" & !`__is_mac' {
            di as error "  The distribution Linux plugin needs x86-64 and glibc 2.17 or later; OpenMP is linked statically."
            di as error "  Run ldd -r on the plugin to inspect missing libraries or symbols."
        }
        else if `__is_mac' {
            di as error "  On macOS, check that the plugin is the notarized ctools release and is not quarantined."
        }
        di as error "  See https://github.com/mdroste/stata-ctools/blob/main/docs/INSTALLATION.md#troubleshooting"
        exit `__load_rc'
    }
    local ctools_plugin_hdfe_screen
    local ctools_plugin_hdfe_markout
    local ctools_plugin_hdfe_fvlevels
    capture plugin call ctools_plugin, "version"
    local __id_rc = _rc
    if `__id_rc' | "`ctools_plugin_protocol'" != "`__ado_protocol'" {
        local __has = cond("`ctools_plugin_protocol'" == "", "none", "`ctools_plugin_protocol'")
        di as error "ctools: the plugin `__plugin' does not match the installed ado files."
        local __ident = cond("`ctools_plugin_version'" == "", "", " (plugin version `ctools_plugin_version', revision `ctools_plugin_revision')")
        di as error "  The ado files need plugin protocol `__ado_protocol'; the loaded plugin reports `__has'`__ident'."
        di as error "  Reinstall ctools so its ado files and plugins come from the same release:"
        di as error `"    `__reinstall'"'
        di as error "  Then restart Stata, because a plugin stays loaded until Stata exits."
        exit 9
    }
    c_local __ctools_plugin "`__plugin'"
    c_local __ctools_plugin_version "`ctools_plugin_version'"
    c_local __ctools_plugin_revision "`ctools_plugin_revision'"
    c_local __ctools_plugin_protocol "`ctools_plugin_protocol'"
    c_local __ctools_plugin_hdfe_screen = ("`ctools_plugin_hdfe_screen'" == "1")
    c_local __ctools_plugin_hdfe_markout = ("`ctools_plugin_hdfe_markout'" == "1")
    c_local __ctools_plugin_hdfe_fvlevels = ("`ctools_plugin_hdfe_fvlevels'" == "1")
    capture confirm number 0
end

* _ctools_strw — set string width metadata for ctools flat buffer optimization
* 
* Usage:  _ctools_strw varlist
* 
* Sets local __ctools_strw in CALLER's scope containing comma-separated
* string variable widths (0 = numeric, >0 = str width e.g. 17 for str17).
* The ctools C plugin reads this local automatically via SF_macro_use
* to enable flat buffer string I/O (eliminates per-string strlen, memcpy,
* and arena atomic CAS overhead).
* 
* Example:
* _ctools_strw `allvars'
* plugin call ctools_plugin `allvars', "mycommand args"
* * C plugin automatically detects __ctools_strw local

program define _ctools_strw
    syntax varlist

    local __strw ""
    local __sep ""
    foreach v of local varlist {
        local __type : type `v'
        if substr("`__type'", 1, 3) == "str" {
            local __w = real(substr("`__type'", 4, .))
            if missing(`__w') local __w = 0
        }
        else {
            local __w = 0
        }
        local __strw "`__strw'`__sep'`__w'"
        local __sep ","
    }
    c_local __ctools_strw `__strw'
end

* Exact undo of a ctools wrapper's changes without a dataset snapshot.
* The wrappers that use it only add observations and variables at the end of
* the dataset, add value labels and _dta characteristics or notes, change
* variable labels they name, and lose the sort marker; the plugin restores any
* values it overwrote first. The record lives in globals named <prefix>*.
*
*   _ctools_undo record <prefix>
*       before the first change: _N, c(k), the sort marker, c(changed), the
*       value-label names and the _dta characteristics (notes included)
*   _ctools_undo varlabels <prefix> <varlist>
*       add the variable labels of varlist to the record
*   _ctools_undo rollback <prefix>
*       drop what was added and put the recorded state back (call under
*       nobreak, after the plugin has restored overwritten values)
*   _ctools_undo clear <prefix>
*       drop the record after success
program define _ctools_undo
    version 14.1
    gettoken sub 0 : 0
    gettoken p 0 : 0
    if "`sub'" == "record" _ctools_undo_record `p'
    else if "`sub'" == "varlabels" _ctools_undo_varlabels `p' `0'
    else if "`sub'" == "rollback" _ctools_undo_rollback `p'
    else if "`sub'" == "clear" macro drop `p'*
    else exit 198
end

program define _ctools_undo_record
    args p
    macro drop `p'*
    global `p'N = _N
    global `p'k = c(k)
    global `p'sortedby : sortedby
    global `p'changed = c(changed)
    tempname caller_r
    _return hold `caller_r'
    quietly label dir
    global `p'labels `r(names)'
    _return restore `caller_r'
    local names : char _dta[]
    global `p'dn `names'
    local i = 0
    foreach c of local names {
        local ++i
        mata: st_global("`p'd`i'", st_global("_dta[`c']"))
    }
    global `p'nl = 0
    global `p'active 1
end

program define _ctools_undo_varlabels
    gettoken p 0 : 0
    if "${`p'active}" == "" exit
    local i = ${`p'nl}
    foreach v of local 0 {
        local ++i
        global `p'lv`i' `v'
        mata: st_global("`p'll`i'", st_varlabel("`v'"))
    }
    global `p'nl = `i'
end

program define _ctools_undo_rollback
    args p
    if "${`p'active}" == "" {
        macro drop `p'*
        exit
    }
    local N0 = ${`p'N}
    local k0 = ${`p'k}
    if _N > `N0' quietly drop in `=`N0'+1'/l
    if c(k) > `k0' {
        mata: st_local("added", invtokens(st_varname((`k0'+1)..st_nvar())))
        drop `added'
    }
    quietly label dir
    local now `r(names)'
    local before ${`p'labels}
    local added : list now - before
    if "`added'" != "" label drop `added'
    local now : char _dta[]
    local old ${`p'dn}
    foreach c of local now {
        if !`: list c in old' char _dta[`c']
    }
    local i = 0
    foreach c of local old {
        local ++i
        mata: st_global("_dta[`c']", st_global("`p'd`i'"))
    }
    forvalues i = 1/${`p'nl} {
        mata: st_varlabel(st_global("`p'lv`i'"), st_global("`p'll`i'"))
    }
    local sortedby ${`p'sortedby}
    if "`: sortedby'" != "`sortedby'" {
        * The observations are back in their order: a stable sort only sets
        * the marker again
        if "`sortedby'" != "" sort `sortedby', stable
        else _ctools_undo_unsort
    }
    mata: st_updata(${`p'changed})
    macro drop `p'*
end

* Drop the sort marker without moving data: change one value of the first sort
* variable with replace, which Stata notices, then put it back through Mata.
program define _ctools_undo_unsort
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_undo_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_undo_hold)
    }
    else {
        mata: __ctools_undo_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_undo_hold)
    }
    mata: mata drop __ctools_undo_hold
end
* ctools-helpers:end
