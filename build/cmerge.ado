*! version 1.2.0 20261029 github.com/mdroste/stata-ctools

* Mata helpers for optimized operations
capture mata: mata drop _cmerge_shared_flags()
capture mata: mata drop _cmerge_chars_record()
capture mata: mata drop _cmerge_chars_restore()
mata:
void _cmerge_shared_flags(string rowvector keepusing, string rowvector master_vars)
{
    real scalar i, j, n_ku, n_mv, is_shared
    string scalar result, vname
    transmorphic A

    n_ku = cols(keepusing)
    n_mv = cols(master_vars)

    /* Build associative array for O(1) lookup */
    A = asarray_create()
    for (j = 1; j <= n_mv; j++) {
        asarray(A, master_vars[j], 1)
    }

    /* Check each keepusing var against the hash */
    result = ""
    for (i = 1; i <= n_ku; i++) {
        is_shared = asarray_contains(A, keepusing[i]) ? 1 : 0
        result = result + (i > 1 ? " " : "") + strofreal(is_shared)
    }

    st_local("shared_var_flags", result)
}

/* Characteristics of evarname ev (a variable or _dta), kept in the globals
   <p>n (the names) and <p>1, <p>2, ... (the values) */
void _cmerge_chars_record(string scalar ev, string scalar p)
{
    string colvector names
    real scalar i

    names = st_dir("char", ev, "*")
    st_global(p + "n", invtokens(names'))
    for (i = 1; i <= rows(names); i++) {
        st_global(p + strofreal(i), st_global(ev + "[" + names[i] + "]"))
    }
}

/* Put back the characteristics _cmerge_chars_record kept: drop the ones added
   since and rewrite the ones whose values changed */
void _cmerge_chars_restore(string scalar ev, string scalar p)
{
    string rowvector old
    string colvector now
    string scalar value
    real scalar i

    old = tokens(st_global(p + "n"))
    now = st_dir("char", ev, "*")
    for (i = 1; i <= rows(now); i++) {
        if (!anyof(old, now[i])) st_global(ev + "[" + now[i] + "]", "")
    }
    for (i = 1; i <= cols(old); i++) {
        value = st_global(p + strofreal(i))
        if (st_global(ev + "[" + old[i] + "]") != value) {
            st_global(ev + "[" + old[i] + "]", value)
        }
    }
}
end

program define cmerge, rclass
    * merge's behavior depends on the caller's version (strL fallback below)
    global CTOOLS_cmerge_caller = _caller()
    global CTOOLS_cmerge_assert_failed
    version 14.1
    local caller_frame ""
    if c(stata_version) >= 16 local caller_frame "`c(frame)'"
    * By default a failed or interrupted merge leaves the master as it was
    * without a snapshot of the whole dataset: _cmerge_impl records what it is
    * about to change (_cmerge_record), the plugin keeps a copy of every master
    * observation before it overwrites one, and _cmerge_rollback undoes both.
    * The few paths that need a snapshot keep it in memory, in this frame.
    * Option snapshot preserves the whole dataset instead (on disk once it
    * exceeds c(max_preservemem)); option fast keeps no copy at all.
    _cmerge_mode `0'
    tempname snapframe
    * Break is recognized only inside _cmerge_impl.
    nobreak {
        macro drop CTOOLS_cmu_* CTOOLS_cmerge_plugin
        global CTOOLS_cmerge_snapframe `snapframe'
        if "`mode'" == "snapshot" preserve
        capture noisily break _cmerge_impl `0'
        local rc = _rc
        global CTOOLS_cmerge_caller
        if "`caller_frame'" != "" frame change `caller_frame'
        * Like merge, a failed assert() leaves the merged result in memory
        if `rc' == 9 & "$CTOOLS_cmerge_assert_failed" == "1" {
            global CTOOLS_cmerge_assert_failed
            _cmerge_commit
            if "`mode'" == "snapshot" restore, not
            exit 9
        }
        global CTOOLS_cmerge_assert_failed
        if `rc' {
            _cmerge_rollback
            if "`mode'" == "snapshot" restore
            exit `rc'
        }
        return add
        _cmerge_commit
        if "`mode'" == "snapshot" restore, not
    }
end

* The undo mode requested in cmerge's options: local mode in the caller is
* snapshot, fast or empty (the default).
program define _cmerge_mode
    gettoken lhs 0 : 0, parse(",")
    syntax [, SNAPshot FAST *]
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
    c_local mode `snapshot'`fast'
end

* Record what _cmerge_impl is about to change in the master, for
* _cmerge_rollback: the observation and variable counts (new variables and
* observations are only ever added at the end), the sort marker, c(changed),
* the value-label names, the storage types and formats of the variables to be
* promoted, and the value labels and characteristics (notes included) of the
* shared variables and of _dta (append copies the using file's).
program define _cmerge_record
    syntax [, Promote(string) Shared(string)]
    macro drop CTOOLS_cmu_*
    global CTOOLS_cmu_N = _N
    global CTOOLS_cmu_k = c(k)
    global CTOOLS_cmu_sortedby : sortedby
    global CTOOLS_cmu_changed = c(changed)
    quietly label dir
    global CTOOLS_cmu_labels `r(names)'
    local i = 0
    foreach v of local promote {
        local ++i
        global CTOOLS_cmu_pv`i' `v'
        global CTOOLS_cmu_pt`i' : type `v'
        global CTOOLS_cmu_pf`i' : format `v'
    }
    global CTOOLS_cmu_np = `i'
    local i = 0
    foreach v of local shared {
        local ++i
        global CTOOLS_cmu_sv`i' `v'
        global CTOOLS_cmu_sl`i' : value label `v'
        mata: _cmerge_chars_record("`v'", "CTOOLS_cmu_c`i'_")
    }
    global CTOOLS_cmu_ns = `i'
    mata: _cmerge_chars_record("_dta", "CTOOLS_cmu_d_")
    global CTOOLS_cmu_active 1
end

* Undo a failed or interrupted merge (called by cmerge under nobreak). The
* plugin first rewrites the original observations if it holds them; then the
* added observations, variables and value labels go, and the recorded types,
* formats, value labels, characteristics, sort marker and c(changed) return.
* A path that kept an in-memory snapshot frame is restored from it instead.
* Without either (failure before any change, or option fast) only the
* plugin's state is freed.
program define _cmerge_rollback
    local active = "$CTOOLS_cmu_active" != ""
    if `active' {
        local N0 = $CTOOLS_cmu_N
        local k0 = $CTOOLS_cmu_k
    }
    if "$CTOOLS_cmerge_plugin" != "" {
        * Observations first, while the original variables still have their
        * positions and the storage types the plugin copied them with
        local vars ""
        if `active' {
            if _N < `N0' quietly set obs `N0'
            if `k0' > 0 mata: st_local("vars", invtokens(st_varname(1..`k0')))
        }
        capture program ctools_plugin, plugin using("$CTOOLS_cmerge_plugin")
        capture noisily plugin call ctools_plugin `vars', "cmerge rollback"
        if _rc di as error "cmerge: the master observations could not be restored"
    }
    if c(stata_version) >= 16 {
        capture confirm frame $CTOOLS_cmerge_snapframe
        if !_rc {
            frame copy $CTOOLS_cmerge_snapframe `c(frame)', replace
            frame drop $CTOOLS_cmerge_snapframe
        }
    }
    if `active' {
        if _N > `N0' quietly drop in `=`N0'+1'/l
        if c(k) > `k0' {
            mata: st_local("added", invtokens(st_varname((`k0'+1)..st_nvar())))
            drop `added'
        }
        forvalues i = 1/$CTOOLS_cmu_np {
            local v ${CTOOLS_cmu_pv`i'}
            capture noisily recast ${CTOOLS_cmu_pt`i'} `v'
            format `v' ${CTOOLS_cmu_pf`i'}
        }
        quietly label dir
        local now `r(names)'
        local before $CTOOLS_cmu_labels
        local added : list now - before
        if "`added'" != "" label drop `added'
        forvalues i = 1/$CTOOLS_cmu_ns {
            local v ${CTOOLS_cmu_sv`i'}
            if "`: value label `v''" != "${CTOOLS_cmu_sl`i'}" {
                label values `v' ${CTOOLS_cmu_sl`i'}
            }
            mata: _cmerge_chars_restore("`v'", "CTOOLS_cmu_c`i'_")
        }
        mata: _cmerge_chars_restore("_dta", "CTOOLS_cmu_d_")
        local sortedby $CTOOLS_cmu_sortedby
        if "`: sortedby'" != "`sortedby'" {
            if "`sortedby'" != "" sort `sortedby', stable
            else _cmerge_clear_sortedby
        }
        mata: st_updata($CTOOLS_cmu_changed)
    }
    macro drop CTOOLS_cmu_* CTOOLS_cmerge_plugin CTOOLS_cmerge_snapframe
end

* Release what the plugin and any snapshot frame hold for a merge that stands
* (success, or a failed assert(), whose merged result stays in memory).
program define _cmerge_commit
    if "$CTOOLS_cmerge_plugin" != "" {
        capture program ctools_plugin, plugin using("$CTOOLS_cmerge_plugin")
        capture plugin call ctools_plugin, "cmerge commit"
    }
    if c(stata_version) >= 16 & "$CTOOLS_cmerge_snapframe" != "" {
        capture frame drop $CTOOLS_cmerge_snapframe
    }
    macro drop CTOOLS_cmu_* CTOOLS_cmerge_plugin CTOOLS_cmerge_snapframe
end

program define _cmerge_impl, rclass
    version 14.1
    local __cmdline `"`0'"'

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Parse merge type and key variables
    gettoken merge_type 0 : 0, parse(" ")

    * Validate merge type
    local merge_code = -1
    if "`merge_type'" == "1:1" {
        local merge_code = 0
    }
    else if "`merge_type'" == "m:1" {
        local merge_code = 1
    }
    else if "`merge_type'" == "1:m" {
        local merge_code = 2
    }
    else if "`merge_type'" == "m:m" {
        local merge_code = 3
    }
    else {
        di as error "cmerge: invalid merge type `merge_type'"
        di as error "Must be one of: 1:1, m:1, 1:m, m:m"
        exit 111
    }

    * Check for _n merge (merge by observation number)
    gettoken first_token rest : 0, parse(" ")
    local merge_by_n = 0
    if "`first_token'" == "_n" {
        * Validate: _n only valid with 1:1 merge
        if "`merge_type'" != "1:1" {
            di as error "cmerge: _n may only be used with 1:1 merge"
            exit 198
        }
        local merge_by_n = 1
        * Remove _n from 0 and continue parsing
        local 0 "`rest'"
    }

    * Parse key variables and using filename
    if `merge_by_n' {
        * For _n merge, no varlist needed
        syntax using/ [, ///
            Keep(string) ///
            ASSert(string) ///
            GENerate(name) ///
            NOGENerate ///
            KEEPUSing(string) ///
            SORTED ///
            FORCE ///
            NOREPort ///
            Verbose ///
            NOLabel ///
            NONotes ///
            UPDATE ///
            REPLACE ///
            PRESERVE_order(integer 0) ///
            THReads(integer 0) ///
            SNAPshot ///
            FAST ///
            ]

        * No key variables for _n merge
        local keyvars ""
        local nkeys = 0
    }
    else {
        * Standard merge with key variables
        syntax varlist using/ [, ///
            Keep(string) ///
            ASSert(string) ///
            GENerate(name) ///
            NOGENerate ///
            KEEPUSing(string) ///
            SORTED ///
            FORCE ///
            NOREPort ///
            Verbose ///
            NOLabel ///
            NONotes ///
            UPDATE ///
            REPLACE ///
            PRESERVE_order(integer 0) ///
            THReads(integer 0) ///
            SNAPshot ///
            FAST ///
            ]

        * Validate key variables exist in master
        local keyvars `varlist'
        local nkeys : word count `keyvars'

        if `nkeys' == 0 {
            di as error "cmerge: no key variables specified"
            exit 198
        }
    }

    * The sorted C join checks every group, including unmatched keys.
    * Retain isid for an explicit sorted promise, which can bypass C sorting.
    if !`merge_by_n' & "`sorted'" != "" & inlist(`merge_code', 0, 2) & _N > 0 {
        capture isid `keyvars', missok
        if _rc {
            di as error "cmerge: key variables do not uniquely identify observations in master"
            exit 459
        }
    }

    * Default generate variable name
    if "`generate'" == "" & "`nogenerate'" == "" {
        local generate "_merge"
    }

    * Check _merge variable doesn't already exist (exact name: _merge_prev
    * must not count as an abbreviation of _merge)
    if "`nogenerate'" == "" {
        capture confirm variable `generate', exact
        if !_rc {
            di as error "variable `generate' already exists"
            di as error "use generate() option or drop `generate' first"
            exit 110
        }
    }

    * replace performs an update, so merge requires update with it
    if "`replace'" != "" & "`update'" == "" {
        di as err "option {bf:replace} not allowed"
        di as err "{p 4 4 2}"
        di as err "option {bf:replace} requires you also specify"
        di as err "option {bf:update}, thus demonstrating your"
        di as err "understanding that an update will be performed"
        di as err "{p_end}"
        exit 198
    }

    * keep() and assert() take merge's result words or the codes 1-5
    * (validated before any data manipulation)
    local keep_codes ""
    if "`keep'" != "" _cmerge_results keep_codes : `keep'
    local assert_codes ""
    if "`assert'" != "" _cmerge_results assert_codes : `assert'

    * Validate using file exists (skip check for web URLs)
    local is_url = 0
    if substr(`"`using'"', 1, 7) == "http://" | substr(`"`using'"', 1, 8) == "https://" {
        local is_url = 1
    }
    if !`is_url' {
        capture confirm file `"`using'"'
        if _rc {
            capture confirm file `"`using'.dta"'
            if _rc {
                di as error `"file `using' not found"'
                exit 601
            }
            local using `"`using'.dta"'
        }
    }

    * Get variable types for key variables in master (skip for _n merge)
    local master_keytypes ""
    if !`merge_by_n' {
        foreach var of local keyvars {
            capture confirm string variable `var'
            if !_rc {
                local master_keytypes "`master_keytypes' str"
            }
            else {
                local master_keytypes "`master_keytypes' num"
            }
        }
    }

    * Store master dataset info BEFORE any modifications
    local master_nobs = _N
    local master_nvars = c(k)
    unab master_varlist : _all
    quietly label dir
    local __master_labels `r(names)'
    local master_storage_types ""
    foreach var of local master_varlist {
        local storage : type `var'
        local master_storage_types `master_storage_types' `storage'
    }

    * The plugin cannot write strL values: merges involving strL variables
    * run through merge itself (checked before any data are changed).
    if `: list posof "strL" in master_storage_types' {
        _cmerge_native `"`__cmdline'"'
        return add
        exit
    }

    * Allow empty master - will just add using-only observations

    * Get master key variable indices (1-based) using Mata - skip for _n merge
    local master_key_indices ""
    local master_sorted = 0
    if !`merge_by_n' {
        * Use Mata st_varindex() for O(1) lookup instead of O(n*m) nested loops
        mata: st_local("master_key_indices", invtokens(strofreal(st_varindex(tokens(st_local("keyvars"))))))
        * Verify all keys found (check for missing values from st_varindex)
        local n_found : word count `master_key_indices'
        if `n_found' != `nkeys' {
            di as error "cmerge: not all key variables found in master"
            exit 111
        }
        * Check for any zeros (variable not found)
        foreach idx of local master_key_indices {
            if `idx' == . {
                di as error "cmerge: not all key variables found in master"
                exit 111
            }
        }

        * Check if master is already sorted on key variables
        * This allows skipping the sort step in the C plugin
        local sortedby : sortedby
        if "`sortedby'" != "" {
            * Check if keyvars are a prefix of sortedby
            local master_sorted = 1
            local i = 1
            foreach k of local keyvars {
                local s : word `i' of `sortedby'
                if "`k'" != "`s'" {
                    local master_sorted = 0
                    continue, break
                }
                local ++i
            }
        }
    }
    else {
        * For _n merge, data is implicitly "sorted" by row number
        local master_sorted = 1
    }

    * Display header if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cmerge: Optimized C-Accelerated Merge"
        di as text "{hline 60}"
        di as text "Merge type:  " as result "`merge_type'"
        if `merge_by_n' {
            di as text "Key vars:    " as result "_n (observation number)"
        }
        else {
            di as text "Key vars:    " as result "`keyvars'"
        }
        di as text "Using file:  " as result `"`using'"'
        di as text "Master obs:  " as result `master_nobs'
        di as text "Master vars: " as result `master_nvars'
        if `master_sorted' {
            di as text "Master sort: " as result "already sorted on keys - skipping sort"
        }
        di as text "{hline 60}"
        di ""
    }

    * Initialize timing (like csort)
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        local __ct_elapsed93 = 0
        local __ct_elapsed94 = 0
        local __ct_elapsed95 = 0
        if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
        quietly _ctools_clock __ct_start90
        if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
        quietly _ctools_clock __ct_start91
    }

    * Load the platform-appropriate ctools plugin (cached after first load)
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0
    * Lets cmerge free or roll back what the plugin holds for this merge
    global CTOOLS_cmerge_plugin "`__ctools_plugin'"

    * =========================================================================
    * Phase 1: Load using dataset - ONLY key + keepusing vars
    * =========================================================================

    if "`verbose'" != "" {
        di as text "Phase 1: Loading using dataset (keys + keepusing only)..."
    }

    * Use frames (Stata 16+) for faster context switching vs preserve/restore
    * Frames avoid copying the entire master dataset twice
    local __use_frames = 0
    tempname __using_frame
    local __master_frame ""
    if c(stata_version) >= 16 {
        local __master_frame "`c(frame)'"
        capture frame create `__using_frame'
        if _rc == 0 {
            local __use_frames = 1
        }
    }

    if `__use_frames' {
        * Load using data into separate frame (master stays in its caller frame)
        * If keepusing specified, load only key + keepusing vars (much faster for wide datasets)
        if "`keepusing'" != "" & !`merge_by_n' {
            frame `__using_frame': qui use `keyvars' `keepusing' using `"`using'"', clear
        }
        else {
            frame `__using_frame': qui use `"`using'"', clear
        }
        frame change `__using_frame'
    }
    else {
        * Fallback: preserve/restore for Stata < 16
        preserve
        if "`keepusing'" != "" & !`merge_by_n' {
            qui use `keyvars' `keepusing' using `"`using'"', clear
        }
        else {
            qui use `"`using'"', clear
        }
    }

    * Check key variables exist in using (skip for _n merge)
    if !`merge_by_n' {
        foreach var of local keyvars {
            capture confirm variable `var'
            if _rc {
                di as error "cmerge: key variable `var' not found in using dataset"
                if `__use_frames' {
                    frame change `__master_frame'
                    frame drop `__using_frame'
                }
                else {
                    restore
                }
                exit 111
            }
        }

        * Empty-master merges bypass the C join and still need uniqueness checks.
        if inlist(`merge_code', 0, 1) & _N > 0 & (`master_nobs' == 0 | "`sorted'" != "") {
            capture isid `keyvars', missok
            if _rc {
                di as error "cmerge: key variables do not uniquely identify observations in using"
                if `__use_frames' {
                    frame change `__master_frame'
                    frame drop `__using_frame'
                }
                else restore
                exit 459
            }
        }

        * Validate key variable types match
        local i = 1
        foreach var of local keyvars {
            local master_type : word `i' of `master_keytypes'
            capture confirm string variable `var'
            if !_rc {
                local using_type "str"
            }
            else {
                local using_type "num"
            }

            if "`master_type'" != "`using_type'" {
                di as error "cmerge: key variable `var' is `master_type' in master but `using_type' in using"
                di as error "key types must match even with force"
                if `__use_frames' {
                    frame change `__master_frame'
                    frame drop `__using_frame'
                }
                else restore
                exit 106
            }
            local ++i
        }
    }

    * Using keys or kept variables stored as strL: run through merge itself
    local __using_strl = 0
    foreach var of varlist _all {
        if "`: type `var''" == "strL" local __using_strl = 1
    }
    if `__using_strl' {
        if `__use_frames' {
            frame change `__master_frame'
            frame drop `__using_frame'
        }
        else restore
        _cmerge_native `"`__cmdline'"'
        return add
        exit
    }

    * Build list of variables to keep: keys + keepusing (no keys for _n merge)
    local using_keep_vars "`keyvars'"
    local keepusing_count = 0
    local keepusing_names ""
    local keepusing_types ""

    if "`keepusing'" != "" {
        * Expand wildcard, abbreviation and range syntax against using data.
        unab keepusing : `keepusing'
        if `merge_by_n' keep `keepusing'
        * Build concrete name/type lists
        foreach var of local keepusing {
            local using_keep_vars "`using_keep_vars' `var'"
            local keepusing_names "`keepusing_names' `var'"
            local ++keepusing_count

            * Capture type for placeholder creation (preserve original type)
            local vtype : type `var'
            local keepusing_types "`keepusing_types' `vtype'"
        }
    }
    else {
        * No keepusing specified - get all non-key variables from using
        * Include variables that exist in master (shared vars) - they'll be
        * written only for using-only rows to avoid overwriting master values
        unab all_using_vars : _all
        foreach var of local all_using_vars {
            local is_key = 0
            foreach k of local keyvars {
                if "`var'" == "`k'" {
                    local is_key = 1
                }
            }
            if !`is_key' {
                local keepusing_names "`keepusing_names' `var'"
                local ++keepusing_count

                * Capture type for placeholder creation (preserve original type)
                local vtype : type `var'
                local keepusing_types "`keepusing_types' `vtype'"
            }
        }
        local using_keep_vars "`keyvars' `keepusing_names'"
    }

    * Resolve all shared destinations, including keys for using-only rows,
    * while both input schemas are available. Recast master only after restore.
    local promote_names ""
    local promote_types ""
    local numeric_types "byte int long float double"
    local schema_vars `keyvars' `keepusing_names'
    local schema_vars : list uniq schema_vars
    foreach var of local schema_vars {
        local master_pos : list posof "`var'" in master_varlist
        if `master_pos' {
            local mt : word `master_pos' of `master_storage_types'
            local ut : type `var'
            local mstr = substr("`mt'",1,3) == "str"
            local ustr = substr("`ut'",1,3) == "str"
            if `mstr' != `ustr' {
                * force applies only to non-key using values: discard values
                * whose representation cannot be stored in the master type.
                local is_key : list var in keyvars
                if "`force'" == "" | `is_key' {
                    di as error "cmerge: incompatible types for shared variable `var'"
                    if `__use_frames' {
                        frame change `__master_frame'
                        frame drop `__using_frame'
                    }
                    else restore
                    exit 106
                }
                drop `var'
                if `mstr' quietly generate `mt' `var' = ""
                else quietly generate `mt' `var' = .
                local ut "`mt'"
            }
            local target "`mt'"
            if `mstr' {
                if "`mt'" == "strL" | "`ut'" == "strL" local target "strL"
                else {
                    local width = max(real(substr("`mt'",4,.)), real(substr("`ut'",4,.)))
                    local target "str`width'"
                }
            }
            else {
                local mr : list posof "`mt'" in numeric_types
                local ur : list posof "`ut'" in numeric_types
                local rank = max(`mr', `ur')
                * float cannot exactly represent every long integer.
                if (`mr' == 3 & `ur' == 4) | (`mr' == 4 & `ur' == 3) local rank = 5
                local target : word `rank' of `numeric_types'
            }
            if "`target'" != "`mt'" {
                local promote_names `promote_names' `var'
                local promote_types `promote_types' `target'
            }
        }
    }
    * A force conversion may have changed the using-side storage metadata.
    local keepusing_types ""
    foreach var of local keepusing_names {
        local storage : type `var'
        local keepusing_types `keepusing_types' `storage'
    }

    * Identify shared variables using Mata (O(n) vs O(n*m) nested loops)
    * A variable is shared if it exists in both master and using
    local shared_var_flags ""
    if `keepusing_count' > 0 {
        mata: _cmerge_shared_flags(tokens(st_local("keepusing_names")), tokens(st_local("master_varlist")))
    }

    local using_nobs = _N

    * Check if using is already sorted on key variables
    local using_sorted = 0
    if !`merge_by_n' {
        local sortedby : sortedby
        if "`sortedby'" != "" {
            local using_sorted = 1
            local i = 1
            foreach k of local keyvars {
                local s : word `i' of `sortedby'
                if "`k'" != "`s'" {
                    local using_sorted = 0
                    continue, break
                }
                local ++i
            }
        }
    }
    else {
        * For _n merge, data is implicitly "sorted" by row number
        local using_sorted = 1
    }

    if "`verbose'" != "" {
        if `merge_by_n' {
            di as text "  Using: " as result `using_nobs' as text " obs, keeping " as result `keepusing_count' as text " keepusing vars (merge by _n)"
        }
        else {
            di as text "  Using: " as result `using_nobs' as text " obs, keeping " as result `nkeys' as text " keys + " as result `keepusing_count' as text " keepusing vars"
            if `using_sorted' {
                di as text "  Using data already sorted on keys - skipping sort"
            }
        }
    }

    * =========================================================================
    * Special handling for empty datasets (handle entirely in Stata)
    * =========================================================================

    if `master_nobs' == 0 | `using_nobs' == 0 {
        * Like merge, the sorted option is checked here too (r(5))
        if "`sorted'" != "" & !`merge_by_n' & `using_nobs' > 1 {
            _cmerge_keys_sorted `keyvars'
            if !r(sorted) {
                if `__use_frames' {
                    frame change `__master_frame'
                    frame drop `__using_frame'
                }
                else restore
                di as error "using data not sorted"
                exit 5
            }
        }
        * Save using data if needed; with no using rows, save the new using
        * variables (and their metadata) as a template instead
        tempfile using_saved
        if `using_nobs' > 0 {
            qui save `using_saved', replace
        }
        else {
            local new_vars ""
            local i = 0
            foreach v of local keepusing_names {
                local ++i
                if !`: word `i' of `shared_var_flags'' local new_vars `new_vars' `v'
            }
            _cmerge_template `"`using_saved'"' "`new_vars'" ""
        }

        if `__use_frames' {
            frame change `__master_frame'
            frame drop `__using_frame'
        }
        else {
            restore
        }

        if "`sorted'" != "" & !`merge_by_n' & `master_nobs' > 1 {
            _cmerge_keys_sorted `keyvars'
            if !r(sorted) {
                di as error "master data not sorted"
                exit 5
            }
        }

        * Empty-using merges also bypass the C join.
        if !`merge_by_n' & `using_nobs' == 0 & `master_nobs' > 0 & ///
            inlist(`merge_code', 0, 2) & "`sorted'" == "" {
            capture isid `keyvars', missok
            if _rc {
                di as error "cmerge: key variables do not uniquely identify observations in master"
                exit 459
            }
        }

        * By default this path is undone from a record, unless it drops master
        * observations (keep()) or may reorder them (1:1 sort): then it keeps
        * an in-memory snapshot frame (preserve before Stata 16)
        local __preserved = 0
        if "`snapshot'`fast'" == "" {
            local __rows_lost = 0
            if `master_nobs' > 0 & "`keep'" != "" {
                local __rows_lost = !`: list posof "1" in keep_codes'
            }
            if `master_nobs' > 1 & `merge_code' == 0 & !`merge_by_n' & ///
                !`preserve_order' & !`master_sorted' {
                local __rows_lost = 1
            }
            if `__rows_lost' & c(stata_version) >= 16 {
                frame copy `c(frame)' $CTOOLS_cmerge_snapframe
            }
            else if `__rows_lost' {
                preserve
                local __preserved = 1
            }
            else {
                local __common `keyvars' `promote_names'
                local __i = 0
                foreach v of local keepusing_names {
                    local ++__i
                    if `: word `__i' of `shared_var_flags'' local __common `__common' `v'
                }
                local __common : list uniq __common
                _cmerge_record, promote(`__common') shared(`__common')
            }
        }
        if "$CTOOLS_cmerge_fault" == "record" exit 1

        local p = 0
        foreach var of local promote_names {
            local ++p
            local storage : word `p' of `promote_types'
            quietly recast `storage' `var'
        }

        * Handle empty master: append using with _merge=2
        if `master_nobs' == 0 & `using_nobs' > 0 {
            qui append using `using_saved', `nolabel' `nonotes'
            if "`nogenerate'" == "" {
                qui gen byte `generate' = 2
                _cmerge_label_mergevar `generate'
            }
            local merge_master = 0
            local merge_using = `using_nobs'
            local merge_matched = 0
        }
        * Handle empty using: keep master with _merge=1
        else if `master_nobs' > 0 & `using_nobs' == 0 {
            * Add the new using variables (no rows, with their metadata)
            qui append using `using_saved', `nolabel' `nonotes'
            if "`nogenerate'" == "" {
                qui gen byte `generate' = 1
                _cmerge_label_mergevar `generate'
            }
            local merge_master = `master_nobs'
            local merge_using = 0
            local merge_matched = 0
        }
        * Handle both empty
        else {
            qui append using `using_saved', `nolabel' `nonotes'
            if "`nogenerate'" == "" {
                qui gen byte `generate' = .
                _cmerge_label_mergevar `generate'
            }
            local merge_master = 0
            local merge_using = 0
            local merge_matched = 0
        }

        * Handle assert() option for empty datasets (a failure keeps the result)
        if "`assert'" != "" {
            capture noisily _cmerge_assert "`assert_codes'" `merge_master' `merge_using' 0 0 0
            if _rc {
                local rc = _rc
                if `__preserved' restore, not
                exit `rc'
            }
        }

        * Handle keep() option for empty datasets
        if "`keep'" != "" {
            local empty_result = cond(`master_nobs' > 0, 1, 2)
            if !`: list empty_result in keep_codes' quietly keep if 0
        }

        * Like merge, a 1:1 result without using-only rows is sorted by keys
        if `merge_code' == 0 & !`merge_by_n' & !`preserve_order' & _N > 0 {
            local kept_using = `merge_using'
            if "`keep'" != "" & !`: list posof "2" in keep_codes' local kept_using = 0
            if `kept_using' == 0 sort `keyvars'
        }

        * Display merge table for empty datasets
        if "`noreport'" == "" {
            _cmerge_table "`update'" `merge_master' `merge_using' `merge_matched' 0 0 "`generate'"
        }
        * Validation hook: fail as Break would once the result is in place
        if "$CTOOLS_cmerge_fault" == "post" exit 1

        * Return results for empty datasets
        return scalar N = _N
        return scalar N_1 = `merge_master'
        return scalar N_2 = `merge_using'
        return scalar N_3 = `merge_matched'
        return scalar time = .
        return local using `"`using'"'
        if `merge_by_n' {
            return local keyvars "_n"
        }
        else {
            return local keyvars "`keyvars'"
        }
        if `__preserved' restore, not
        exit 0
    }

    * Get variable indices in the reduced using dataset using Mata
    local using_key_indices ""
    local using_keepusing_indices ""
    unab using_varlist : _all
    if !`merge_by_n' & `nkeys' > 0 {
        * Use Mata for O(1) key index lookup
        mata: st_local("using_key_indices", invtokens(strofreal(st_varindex(tokens(st_local("keyvars"))))))
    }
    if `keepusing_count' > 0 {
        * Use Mata for O(1) keepusing index lookup
        mata: st_local("using_keepusing_indices", invtokens(strofreal(st_varindex(tokens(st_local("keepusing_names"))))))
    }

    * Build plugin command for load_using
    local plugin_args "load_using `nkeys' `using_key_indices'"
    local plugin_args "`plugin_args' n_keepusing `keepusing_count'"
    if `keepusing_count' > 0 {
        local plugin_args "`plugin_args' keepusing_indices `using_keepusing_indices'"
    }
    if "`verbose'" != "" {
        local plugin_args "`plugin_args' verbose"
    }
    * The sorted option is verified (r(5) if false); a sort marker found on
    * the keys is a hint: the plugin sorts if the data are not in key order
    if "`sorted'" != "" | (`using_sorted' & `merge_by_n') {
        local plugin_args "`plugin_args' sorted"
    }
    else if `using_sorted' {
        local plugin_args "`plugin_args' sorted_hint"
    }
    if `merge_by_n' {
        local plugin_args "`plugin_args' merge_by_n"
    }

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * End pre-plugin1 timer, start plugin1 timer
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
        if "`__ct_elapsed92'" == "" local __ct_elapsed92 = 0
        quietly _ctools_clock __ct_start92
    }

    * Set string width metadata for flat buffer optimization
    _ctools_strw `using_varlist'

    * Call plugin Phase 1 with reduced varlist (using_varlist already computed)
    capture noisily plugin call ctools_plugin `using_varlist', "cmerge `threads_code' `plugin_args'"
    local plugin_rc = _rc

    * End plugin1 timer
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')
    }

    if `plugin_rc' {
        di as error "cmerge: failed to load using data (error `plugin_rc')"
        if `__use_frames' {
            frame change `__master_frame'
            frame drop `__using_frame'
        }
        else {
            restore
        }
        exit `plugin_rc'
    }

    * Start inter-plugin timer
    if `__do_timing' {
        if "`__ct_elapsed93'" == "" local __ct_elapsed93 = 0
        quietly _ctools_clock __ct_start93
    }

    * =========================================================================
    * Phase 2: Execute merge with streaming
    * =========================================================================

    if "`verbose'" != "" {
        di as text ""
        di as text "Phase 2: Executing merge with streaming..."
    }

    * Sub-timers for inter-plugin breakdown
    if `__do_timing' {
        local __ct_elapsed96 = 0
        local __ct_elapsed97 = 0
        local __ct_elapsed98 = 0
        if "`__ct_elapsed97'" == "" local __ct_elapsed97 = 0
        quietly _ctools_clock __ct_start97
    }

    * Build lists of new variables to create (excluding existing shared vars)
    * Use shared_var_flags (already computed) to determine which vars need creation
    local new_var_names ""
    local new_var_types ""
    local placeholder_num = 1
    foreach vtype of local keepusing_types {
        local vname : word `placeholder_num' of `keepusing_names'
        local is_shared : word `placeholder_num' of `shared_var_flags'
        if `is_shared' == 0 {
            local new_var_names "`new_var_names' `vname'"
            local new_var_types "`new_var_types' `vtype'"
        }
        local ++placeholder_num
    }
    * Create merge variable if:
    * 1. nogenerate is not specified, OR
    * 2. keep() is specified (need _merge for filtering, will drop later)
    local need_merge_var = 0
    local temp_merge_var = ""
    if "`nogenerate'" == "" {
        local need_merge_var = 1
        local new_var_names "`new_var_names' `generate'"
        local new_var_types "`new_var_types' byte"
    }
    else if "`keep'" != "" {
        * Need temporary merge variable for keep() filtering; a tempvar
        * cannot be the name of a master variable
        local need_merge_var = 1
        tempvar temp_merge_var
        local new_var_names "`new_var_names' `temp_merge_var'"
        local new_var_types "`new_var_types' byte"
    }

    * Count new variables
    local n_new_vars : word count `new_var_names'

    * Like merge, copy every using value-label definition whose name the
    * master does not define (not only those of new variables), give shared
    * variables without a master label the using label under update, and
    * append the using notes of shared variables to the master's
    tempfile __using_labels
    local __new_labels ""
    if "`nolabel'" == "" {
        quietly label dir
        local __using_label_names `r(names)'
        local __new_labels : list __using_label_names - __master_labels
        if "`__new_labels'" != "" {
            quietly label save `__new_labels' using `"`__using_labels'"', replace
        }
    }
    local __shared_i = 0
    local __shared_vars ""
    local __ki = 0
    foreach v of local keepusing_names {
        local ++__ki
        if !`: word `__ki' of `shared_var_flags'' continue
        local ++__shared_i
        local __shared_vars `__shared_vars' `v'
        local __sv_`__shared_i' `v'
        local __sv_lbl_`__shared_i' : value label `v'
        local __sv_nn_`__shared_i' = 0
        if "`nonotes'" == "" {
            local __n : char `v'[note0]
            capture confirm integer number `__n'
            if !_rc {
                forvalues j = 1/`__n' {
                    local __sv_note_`__shared_i'_`j' : char `v'[note`j']
                }
                local __sv_nn_`__shared_i' = `__n'
            }
        }
    }
    local __n_shared = `__shared_i'

    * Create empty template dataset with new variable definitions
    * Using data is already cached in plugin, so we can clear and reuse memory
    * Benchmark shows template+append is 3.5x faster than Mata st_addvar
    tempfile empty_vars_template
    if `n_new_vars' > 0 {
        _cmerge_template `"`empty_vars_template'"' "`new_var_names'" ///
            "`generate' `temp_merge_var'"
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed97 = `__ct_elapsed97' + (`__ct_now' - `__ct_start97')
        if "`__ct_elapsed96'" == "" local __ct_elapsed96 = 0
        quietly _ctools_clock __ct_start96
    }

    if `__use_frames' {
        frame change `__master_frame'
        frame drop `__using_frame'
    }
    else {
        restore
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed96 = `__ct_elapsed96' + (`__ct_now' - `__ct_start96')
        quietly _ctools_clock __ct_plan_start
    }

    * Plan the join in the untouched master: key, sortedness and uniqueness
    * errors stop here with nothing to undo.
    local plan_args "plan `merge_code' `nkeys' `master_key_indices'"
    local plan_args "`plan_args' master_nobs `master_nobs' preserve_order `preserve_order'"
    if "`sorted'" != "" | (`master_sorted' & `merge_by_n') {
        local plan_args "`plan_args' sorted"
    }
    else if `master_sorted' {
        local plan_args "`plan_args' sorted_hint"
    }
    * keep() drops observations, which the plugin's copy must then cover;
    * snapshot and fast need no copy at all
    if "`keep'" != "" local plan_args "`plan_args' hold_undo"
    if "`snapshot'`fast'" != "" local plan_args "`plan_args' nocopy"
    if `merge_by_n' local plan_args "`plan_args' merge_by_n"
    local __ctools_strw ""
    if `master_nvars' > 0 _ctools_strw `master_varlist'
    capture noisily plugin call ctools_plugin `master_varlist', "cmerge `threads_code' `plan_args'"
    local plugin_rc = _rc
    if `plugin_rc' {
        di as error "cmerge: C plugin merge failed (error `plugin_rc')"
        exit `plugin_rc'
    }
    local __plan_N = scalar(_cmerge_plan_N)
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __time_plan = `__ct_now' - `__ct_plan_start'
        if "`__ct_elapsed97'" == "" local __ct_elapsed97 = 0
        quietly _ctools_clock __ct_start97
    }

    if "`snapshot'`fast'" == "" {
        _cmerge_record, promote(`promote_names') shared(`__shared_vars')
    }
    * Validation hooks: fail as Break would before or after the changes below
    if "$CTOOLS_cmerge_fault" == "record" exit 1

    local p = 0
    foreach var of local promote_names {
        local ++p
        local storage : word `p' of `promote_types'
        quietly recast `storage' `var'
    }

    * Append empty template to add variable definitions (faster than st_addvar).
    * append carries the new variables' formats, labels, notes and
    * characteristics; like merge, a value-label name already defined in the
    * master keeps the master's definition.
    if "`__new_labels'" != "" {
        quietly do `"`__using_labels'"'
    }
    if `n_new_vars' > 0 {
        qui append using `empty_vars_template', `nolabel' `nonotes'
    }
    forvalues i = 1/`__n_shared' {
        local v `__sv_`i''
        if "`update'" != "" & "`nolabel'" == "" & "`__sv_lbl_`i''" != "" & ///
            "`: value label `v''" == "" {
            capture label values `v' `__sv_lbl_`i''
        }
        if `__sv_nn_`i'' > 0 {
            local __n : char `v'[note0]
            capture confirm integer number `__n'
            if _rc local __n = 0
            forvalues j = 1/`__sv_nn_`i'' {
                local ++__n
                char `v'[note`__n'] `"`macval(__sv_note_`i'_`j')'"'
            }
            char `v'[note0] `__n'
        }
    }

    * Compute keepusing variable indices using Mata (O(1) lookup vs O(n*m) loops)
    local keepusing_placeholder_indices ""
    if `keepusing_count' > 0 {
        mata: st_local("keepusing_placeholder_indices", invtokens(strofreal(st_varindex(tokens(st_local("keepusing_names"))))))
    }

    * Get _merge variable index (the plugin writes the merge codes there)
    local merge_var_idx = 0
    if `need_merge_var' {
        local __merge_var `generate'
        if "`temp_merge_var'" != "" local __merge_var `temp_merge_var'
        mata: st_local("merge_var_idx", strofreal(st_varindex("`__merge_var'")))
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed97 = `__ct_elapsed97' + (`__ct_now' - `__ct_start97')
        if "`__ct_elapsed98'" == "" local __ct_elapsed98 = 0
        quietly _ctools_clock __ct_start98
    }

    * Expand the dataset to the planned output size
    local max_output = `__plan_N'
    if `max_output' > _N {
        qui set obs `max_output'
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed98 = `__ct_elapsed98' + (`__ct_now' - `__ct_start98')
    }

    * Build plugin command for execute
    local current_nvars = c(k)
    local plugin_args "execute `merge_code' `nkeys' `master_key_indices'"
    local plugin_args "`plugin_args' master_nobs `master_nobs'"
    local plugin_args "`plugin_args' master_nvars `current_nvars'"
    local plugin_args "`plugin_args' n_keepusing `keepusing_count'"
    if `keepusing_count' > 0 {
        local plugin_args "`plugin_args' keepusing_placeholders `keepusing_placeholder_indices'"
    }
    if `need_merge_var' {
        local plugin_args "`plugin_args' merge_var_idx `merge_var_idx'"
    }
    local plugin_args "`plugin_args' preserve_order `preserve_order'"
    if "`verbose'" != "" {
        local plugin_args "`plugin_args' verbose"
    }
    * The sorted option is verified; a master sort marker is only a hint
    if "`sorted'" != "" | (`master_sorted' & `merge_by_n') {
        local plugin_args "`plugin_args' sorted"
    }
    else if `master_sorted' {
        local plugin_args "`plugin_args' sorted_hint"
    }
    if "`update'" != "" {
        local plugin_args "`plugin_args' update"
    }
    if "`replace'" != "" {
        local plugin_args "`plugin_args' replace"
    }
    if `keepusing_count' > 0 {
        local plugin_args "`plugin_args' shared_flags `shared_var_flags'"
    }
    if `merge_by_n' {
        local plugin_args "`plugin_args' merge_by_n"
    }
    if "$CTOOLS_cmerge_fault" == "schema" exit 1
    if "$CTOOLS_cmerge_fault" == "store" local plugin_args "`plugin_args' fault_store"

    * End inter-plugin timer, start plugin2 timer
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed93 = `__ct_elapsed93' + (`__ct_now' - `__ct_start93')
        if "`__ct_elapsed94'" == "" local __ct_elapsed94 = 0
        quietly _ctools_clock __ct_start94
    }

    * Set string width metadata for flat buffer optimization
    unab __all_phase2_vars : _all
    _ctools_strw `__all_phase2_vars'

    * Pass every current variable: SPI indices are relative to this varlist.
    capture noisily plugin call ctools_plugin `__all_phase2_vars', "cmerge `threads_code' `plugin_args'"
    local plugin_rc = _rc

    * End plugin2 timer, start post-plugin timer
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed94 = `__ct_elapsed94' + (`__ct_now' - `__ct_start94')
        if "`__ct_elapsed95'" == "" local __ct_elapsed95 = 0
        quietly _ctools_clock __ct_start95
    }

    if `plugin_rc' {
        di as error "cmerge: C plugin merge failed (error `plugin_rc')"
        exit `plugin_rc'
    }
    * Validation hook: a Break pressed during the plugin is acted on here
    if "$CTOOLS_cmerge_fault" == "post" exit 1

    * Get merge results from scalars (with update, matched rows split into
    * 3 not updated, 4 missing updated, 5 nonmissing conflict)
    local merge_master = scalar(_cmerge_N_1)
    local merge_using = scalar(_cmerge_N_2)
    local merge_matched = scalar(_cmerge_N_3)
    local merge_updated = scalar(_cmerge_N_4)
    local merge_conflict = scalar(_cmerge_N_5)
    local same_order = scalar(_cmerge_same_order)
    local total_obs = scalar(_cmerge_N)

    * Trim excess observations (drop in is 3x faster than keep in)
    if `total_obs' < _N {
        local __drop_start = `total_obs' + 1
        qui drop in `__drop_start'/l
    }

    if "`nogenerate'" == "" _cmerge_label_mergevar `generate'

    * =========================================================================
    * Phase 3: Apply keep/assert options and display results
    * =========================================================================

    * merge checks assert() before applying keep()
    if "`assert'" != "" {
        _cmerge_assert "`assert_codes'" `merge_master' `merge_using' ///
            `merge_matched' `merge_updated' `merge_conflict'
    }

    * Handle keep option (already validated upfront)
    if "`keep'" != "" {
        * Use the appropriate merge variable for filtering
        local filter_var "`generate'"
        if "`temp_merge_var'" != "" {
            local filter_var "`temp_merge_var'"
        }
        local keep_list : subinstr local keep_codes " " ",", all
        qui keep if inlist(`filter_var', `keep_list')

        * Drop temp filter variable if we created one
        if "`temp_merge_var'" != "" {
            qui drop `temp_merge_var'
        }
    }
    if "$CTOOLS_cmerge_fault" == "keep" exit 1

    * Rows were reordered or appended, so Stata's previous sort marker may no
    * longer hold. Like merge, a 1:1 result without using-only rows is sorted
    * by the keys; otherwise clear the marker. When the output kept the master
    * order the plugin may hold no copy of the master, so the sort (a check of
    * data already in key order) and the marker reset must not be interrupted.
    local kept_using = `merge_using'
    if "`keep'" != "" & !`: list posof "2" in keep_codes' local kept_using = 0
    if `merge_code' == 0 & !`merge_by_n' & `kept_using' == 0 & !`preserve_order' {
        if `same_order' nobreak sort `keyvars'
        else sort `keyvars'
    }
    else if !`same_order' {
        nobreak _cmerge_clear_sortedby
    }

    * End all timers
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed95 = `__ct_elapsed95' + (`__ct_now' - `__ct_start95')
        quietly _ctools_clock __ct_now
        local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

        * Extract ado-file timer values
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
        * Elapsed seconds are held in caller-local timing state.
        local __time_preplugin1 = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_plugin1 = `__ct_elapsed92'
        * Elapsed seconds are held in caller-local timing state.
        * The plan call runs between the two timed plugin calls; count it
        * as plugin time (the plugin's Phase 2 total includes it)
        local __time_interplugin = `__ct_elapsed93' - `__time_plan'
        * Elapsed seconds are held in caller-local timing state.
        local __time_plugin2 = `__ct_elapsed94' + `__time_plan'
        * Elapsed seconds are held in caller-local timing state.
        local __time_postplugin = `__ct_elapsed95'

        * Extract inter-plugin sub-timers
        * Elapsed seconds are held in caller-local timing state.
        local __time_restore = `__ct_elapsed96'
        * Elapsed seconds are held in caller-local timing state.
        local __time_addvar = `__ct_elapsed97'
        * Elapsed seconds are held in caller-local timing state.
        local __time_setobs = `__ct_elapsed98'

        * Get C plugin timing from scalars (Phase 1)
        local __p1_load_keys = scalar(_cmerge_p1_load_keys)
        local __p1_load_keepusing = scalar(_cmerge_p1_load_keepusing)
        local __p1_sort = scalar(_cmerge_p1_sort)
        local __p1_apply_perm = scalar(_cmerge_p1_apply_perm)
        local __p1_total = scalar(_cmerge_p1_total)

        * Get C plugin timing from scalars (Phase 2)
        local __p2_load_master = scalar(_cmerge_p2_load_master)
        local __p2_sort_master = scalar(_cmerge_p2_sort_master)
        local __p2_merge_join = scalar(_cmerge_p2_merge_join)
        local __p2_reorder = scalar(_cmerge_p2_reorder)
        local __p2_permute = scalar(_cmerge_p2_permute)
        local __p2_store = scalar(_cmerge_p2_store)
        local __p2_write_meta = scalar(_cmerge_p2_write_meta)
        local __p2_cleanup = scalar(_cmerge_p2_cleanup)
        local __p2_total = scalar(_cmerge_p2_total)
        local __n_output_vars = scalar(_cmerge_n_output_vars)

        * Calculate overhead
        local __p1_overhead = `__time_plugin1' - `__p1_total'
        local __p2_overhead = `__time_plugin2' - `__p2_total'
        local __c_total = `__p1_total' + `__p2_total'
        local __ado_total = `__time_preplugin1' + `__time_interplugin' + `__time_postplugin'
        local __overhead_total = `__p1_overhead' + `__p2_overhead'
    }

    * Display merge table
    if "`noreport'" == "" {
        * Counts of the categories kept by keep()
        forvalues c = 1/5 {
            local kept`c' = ("`keep'" == "" | `: list c in keep_codes')
        }
        _cmerge_table "`update'" `=`merge_master'*`kept1'' `=`merge_using'*`kept2'' ///
            `=`merge_matched'*`kept3'' `=`merge_updated'*`kept4'' ///
            `=`merge_conflict'*`kept5'' "`generate'"
    }

    * Display comprehensive timing breakdown
    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "cmerge timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals (Phase 1: Load using data):"
        di as text "    Load keys:              " as result %8.4f `__p1_load_keys' " sec"
        di as text "    Load keepusing:         " as result %8.4f `__p1_load_keepusing' " sec"
        di as text "    Sort using:             " as result %8.4f `__p1_sort' " sec"
        di as text "    Apply permutation:      " as result %8.4f `__p1_apply_perm' " sec"
        di as text "  {hline 53}"
        di as text "    Phase 1 C total:        " as result %8.4f `__p1_total' " sec"
        di as text "  {hline 53}"
        di as text "  C plugin internals (Phase 2: Plan and execute merge):"
        di as text "    Load master:            " as result %8.4f `__p2_load_master' " sec"
        di as text "    Sort master:            " as result %8.4f `__p2_sort_master' " sec"
        di as text "    Merge join:             " as result %8.4f `__p2_merge_join' " sec"
        di as text "    Reorder output:         " as result %8.4f `__p2_reorder' " sec"
        di as text "    Permute data:           " as result %8.4f `__p2_permute' " sec"
        di as text "    Store to Stata:         " as result %8.4f `__p2_store' " sec"
        di as text "    Write metadata:         " as result %8.4f `__p2_write_meta' " sec"
        di as text "    Cleanup:                " as result %8.4f `__p2_cleanup' " sec"
        di as text "  {hline 53}"
        di as text "    Phase 2 C total:        " as result %8.4f `__p2_total' " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f `__c_total' " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:"
        di as text "    Pre-plugin1 parsing:    " as result %8.4f `__time_preplugin1' " sec"
        di as text "    Plugin1 call overhead:  " as result %8.4f `__p1_overhead' " sec"
        di as text "    Inter-plugin work:      " as result %8.4f `__time_interplugin' " sec"
        di as text "    Plugin2 call overhead:  " as result %8.4f `__p2_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "  {hline 53}"
        di as text "    Stata overhead total:   " as result %8.4f (`__ado_total' + `__overhead_total') " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = _cmerge_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cmerge_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clear timers
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        local __ct_elapsed93 = 0
        local __ct_elapsed94 = 0
        local __ct_elapsed95 = 0
        local __ct_elapsed96 = 0
        local __ct_elapsed97 = 0
        local __ct_elapsed98 = 0
    }

    * Return results
    if `__do_timing' {
        local elapsed = `__time_total'
    }
    else {
        local elapsed = .
    }
    return scalar N = _N
    return scalar N_1 = `merge_master'
    return scalar N_2 = `merge_using'
    return scalar N_3 = `merge_matched'
    if "`update'" != "" {
        return scalar N_4 = `merge_updated'
        return scalar N_5 = `merge_conflict'
    }
    return scalar time = `elapsed'
    return local using `"`using'"'
    if `merge_by_n' {
        return local keyvars "_n"
    }
    else {
        return local keyvars "`keyvars'"
    }

end

* Result codes for keep() and assert(): the numbers 1-5 or merge's words
* (masters, usings, matches/matched, match_updates, match_conflicts and the
* abbreviations merge accepts).

* merge's message when assert() fails (r(9)); counts are _merge 1..5.

* merge's result table; with update, matched rows are split by _merge 3/4/5.

* Clear Stata's sort marker after the plugin reordered or added rows: a real
* change to a marked variable drops the marker, so change observation 1 of the
* first marked variable and restore it.

* Save the using data's new variables with no rows (plus the byte merge
* variables) as the template the master appends. Run in the using data.

* merge's value and variable labels for the generated _merge variable

* Run the merge with merge itself (strL variables are involved), passing the
* options merge accepts; cmerge-only options (threads, verbose,
* preserve_order) do not apply.

* r(sorted) = 1 when the observations are in ascending order of varlist
* (Stata's sort order: strings bytewise, missing values above numbers)
