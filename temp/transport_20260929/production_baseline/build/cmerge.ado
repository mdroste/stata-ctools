*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

* Mata helpers for optimized operations
capture mata: mata drop _cmerge_shared_flags()
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
end

program define cmerge, rclass
    * merge's behavior depends on the caller's version (strL fallback below)
    global CTOOLS_cmerge_caller = _caller()
    global CTOOLS_cmerge_assert_failed
    version 14.1
    local caller_frame ""
    if c(stata_version) >= 16 local caller_frame "`c(frame)'"
    preserve
    capture noisily _cmerge_impl `0'
    local rc = _rc
    global CTOOLS_cmerge_caller
    if "`caller_frame'" != "" frame change `caller_frame'
    * Like merge, a failed assert() leaves the merged result in memory
    if `rc' == 9 & "$CTOOLS_cmerge_assert_failed" == "1" {
        global CTOOLS_cmerge_assert_failed
        restore, not
        exit 9
    }
    global CTOOLS_cmerge_assert_failed
    if `rc' {
        restore
        exit `rc'
    }
    return add
    restore, not
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

    * Check _merge variable doesn't already exist
    if "`nogenerate'" == "" {
        capture confirm variable `generate'
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
        timer clear 90
        timer clear 91
        timer clear 92
        timer clear 93
        timer clear 94
        timer clear 95
        timer on 90   /* Total wall clock */
        timer on 91   /* Pre-plugin1 */
    }

    * Load the platform-appropriate ctools plugin (cached after first load)
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

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
        * When selective use was done, vars are already validated by Stata's use
        * Just build the name/type lists
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

        * Handle keep() option for empty datasets
        if "`keep'" != "" {
            local empty_result = cond(`master_nobs' > 0, 1, 2)
            if !`: list empty_result in keep_codes' quietly keep if 0
        }

        * Handle assert() option for empty datasets
        if "`assert'" != "" {
            _cmerge_assert "`assert_codes'" `merge_master' `merge_using' 0 0 0
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
        timer off 91
        timer on 92   /* Plugin Phase 1 */
    }

    * Set string width metadata for flat buffer optimization
    _ctools_strw `using_varlist'

    * Call plugin Phase 1 with reduced varlist (using_varlist already computed)
    capture noisily plugin call ctools_plugin `using_varlist', "cmerge `threads_code' `plugin_args'"
    local plugin_rc = _rc

    * End plugin1 timer
    if `__do_timing' {
        timer off 92
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
        timer on 93   /* Inter-plugin (restore, create vars) */
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
        timer clear 96
        timer clear 97
        timer clear 98
        timer on 97   /* create vars template */
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
        * Need temporary merge variable for keep() filtering
        local need_merge_var = 1
        local temp_merge_var "_merge_temp_filter"
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
    local __ki = 0
    foreach v of local keepusing_names {
        local ++__ki
        if !`: word `__ki' of `shared_var_flags'' continue
        local ++__shared_i
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
        timer off 97
        timer on 96   /* frame switch / restore */
    }

    if `__use_frames' {
        frame change `__master_frame'
        frame drop `__using_frame'
    }
    else {
        restore
    }

    if `__do_timing' {
        timer off 96
        timer on 97   /* append vars */
    }

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

    * Get _merge variable index
    local merge_var_idx = 0
    if `need_merge_var' {
        local merge_var_idx = c(k)
    }

    if `__do_timing' {
        timer off 97
        timer on 98   /* set obs */
    }

    * Expand dataset to maximum possible output size
    local max_output = `master_nobs' + `using_nobs'
    if `max_output' > _N {
        qui set obs `max_output'
    }

    if `__do_timing' {
        timer off 98
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

    * End inter-plugin timer, start plugin2 timer
    if `__do_timing' {
        timer off 93
        timer on 94   /* Plugin Phase 2 */
    }

    * Set string width metadata for flat buffer optimization
    unab __all_phase2_vars : _all
    _ctools_strw `__all_phase2_vars'

    * Pass every current variable: SPI indices are relative to this varlist.
    capture noisily plugin call ctools_plugin `__all_phase2_vars', "cmerge `threads_code' `plugin_args'"
    local plugin_rc = _rc

    * End plugin2 timer, start post-plugin timer
    if `__do_timing' {
        timer off 94
        timer on 95   /* Post-plugin */
    }

    if `plugin_rc' {
        di as error "cmerge: C plugin merge failed (error `plugin_rc')"
        exit `plugin_rc'
    }

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

    * Rows were reordered or appended, so Stata's previous sort marker may no
    * longer hold. Like merge, a 1:1 result without using-only rows is sorted
    * by the keys; otherwise clear the marker.
    local kept_using = `merge_using'
    if "`keep'" != "" & !`: list posof "2" in keep_codes' local kept_using = 0
    if `merge_code' == 0 & !`merge_by_n' & `kept_using' == 0 & !`preserve_order' {
        sort `keyvars'
    }
    else if !`same_order' {
        _cmerge_clear_sortedby
    }

    * End all timers
    if `__do_timing' {
        timer off 95   /* Post-plugin */
        timer off 90   /* Total wall clock */

        * Extract ado-file timer values
        quietly timer list 90
        local __time_total = r(t90)
        quietly timer list 91
        local __time_preplugin1 = r(t91)
        quietly timer list 92
        local __time_plugin1 = r(t92)
        quietly timer list 93
        local __time_interplugin = r(t93)
        quietly timer list 94
        local __time_plugin2 = r(t94)
        quietly timer list 95
        local __time_postplugin = r(t95)

        * Extract inter-plugin sub-timers
        quietly timer list 96
        local __time_restore = r(t96)
        quietly timer list 97
        local __time_addvar = r(t97)
        quietly timer list 98
        local __time_setobs = r(t98)

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
        di as text "  C plugin internals (Phase 2: Execute merge):"
        di as text "    Load master keys:       " as result %8.4f `__p2_load_master' " sec"
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
        timer clear 90
        timer clear 91
        timer clear 92
        timer clear 93
        timer clear 94
        timer clear 95
        timer clear 96
        timer clear 97
        timer clear 98
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
program define _cmerge_results
    gettoken macname 0 : 0
    gettoken colon 0 : 0
    local codes ""
    foreach w of local 0 {
        local l = strlen(`"`w'"')
        local c ""
        if inlist(`"`w'"', "1", "2", "3", "4", "5") local c `w'
        else if substr("masters", 1, max(3, `l')) == `"`w'"' local c 1
        else if substr("usings", 1, max(2, `l')) == `"`w'"' local c 2
        else if substr("matches", 1, max(3, `l')) == `"`w'"' local c 3
        else if substr("matched", 1, max(3, `l')) == `"`w'"' local c 3
        else if substr("match_updates", 1, max(8, `l')) == `"`w'"' local c 4
        else if substr("match_conflicts", 1, max(8, `l')) == `"`w'"' local c 5
        else {
            di as err `"`w':  invalid {it:resulttype}"'
            di as err "{p 4 4 2}"
            di as err "results, specified in options {bf:assert()} and {bf:keep()},"
            di as err "are the integers 1 through 5, or"
            di as err "{bf:master} (equivalent to 1),"
            di as err "{bf:using} (2),"
            di as err "{bf:match} (3),"
            di as err "{bf:match_update} (4), or"
            di as err "{bf:match_conflict} (5).  The last two arise only when"
            di as err "option {cmd:update} is specified."
            di as err "{p_end}"
            exit 198
        }
        local codes `codes' `c'
    }
    local codes : list uniq codes
    local codes : list sort codes
    c_local `macname' "`codes'"
end

* merge's message when assert() fails (r(9)); counts are _merge 1..5.
program define _cmerge_assert
    args codes n1 n2 n3 n4 n5
    local failed = 0
    forvalues c = 1/5 {
        if !`: list c in codes' & `n`c'' > 0 local failed = 1
    }
    if !`failed' exit
    local n : word count `codes'
    local comma = cond(`n' > 2, ",", "")
    di as err "{p 0 4 2}"
    di as err "after {bf:merge}, not all observations"
    local i = 0
    foreach c of local codes {
        local ++i
        if `c' == 1 local text "from master"
        else if `c' == 2 local text "from using"
        else if `c' == 3 local text "matched"
        else if `c' == 4 local text "matched and updated"
        else local text "matched and conflicting"
        if `i' == 1 di as err "`text'`comma'"
        else if `i' < `n' di as err "`text'`comma'"
        else di as err "or `text'"
    }
    di as err "{p_end}"
    di as err "(merged result left in memory)"
    global CTOOLS_cmerge_assert_failed 1
    exit 9
end

* merge's result table; with update, matched rows are split by _merge 3/4/5.
program define _cmerge_table
    args isupdate m1 m2 m3 m4 m5 mergevar
    if "`mergevar'" != "" {
        forvalues c = 1/5 {
            local v`c' "(`mergevar'==`c')"
        }
    }
    di
    di as txt _col(5) "Result" _col(33) "Number of obs"
    di as txt _col(5) "{hline 41}"
    di as txt _col(5) "Not matched" _col(30) as res %16.0fc (`m1'+`m2')
    if (`m1' | `m2') {
        di as txt _col(9) "from master" _col(30) as res %16.0fc `m1' as txt "  `v1'"
        di as txt _col(9) "from using" _col(30) as res %16.0fc `m2' as txt "  `v2'"
        di
    }
    if "`isupdate'" == "" {
        di as txt _col(5) "Matched" _col(30) as res %16.0fc `m3' as txt "  `v3'"
    }
    else {
        if (`m1' == 0 & `m2' == 0) di
        di as txt _col(5) "Matched" _col(30) as res %16.0fc (`m3'+`m4'+`m5')
        di as txt _col(9) "not updated" _col(30) as res %16.0fc `m3' as txt "  `v3'"
        di as txt _col(9) "missing updated" _col(30) as res %16.0fc `m4' as txt "  `v4'"
        di as txt _col(9) "nonmissing conflict" _col(30) as res %16.0fc `m5' as txt "  `v5'"
    }
    di as txt _col(5) "{hline 41}"
end

* Clear Stata's sort marker after the plugin reordered or added rows: a real
* change to a marked variable drops the marker, so change observation 1 of the
* first marked variable and restore it.
program define _cmerge_clear_sortedby
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    * Stata drops the sort marker when a sort variable changes: change one
    * value with replace, then put the original back through Mata (no copy
    * of the variable is needed)
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_sb_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_sb_hold)
    }
    else {
        mata: __ctools_sb_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_sb_hold)
    }
    mata: mata drop __ctools_sb_hold
end

* Save the using data's new variables with no rows (plus the byte merge
* variables) as the template the master appends. Run in the using data.
program define _cmerge_template
    args file newvars mergevars
    local keepvars : list newvars - mergevars
    if "`keepvars'" != "" {
        keep `keepvars'
        if _N > 0 {
            quietly keep in 1
            quietly drop in 1
        }
    }
    else drop _all
    foreach v of local mergevars {
        quietly generate byte `v' = .
    }
    quietly save `"`file'"', emptyok replace
end

* merge's value and variable labels for the generated _merge variable
program define _cmerge_label_mergevar
    args mergevar
    capture label list _merge
    if _rc {
        label define _merge 1 "Master only (1)" 2 "Using only (2)" ///
            3 "Matched (3)" 4 "Missing updated (4)" 5 "Nonmissing conflict (5)"
    }
    label values `mergevar' _merge
    label variable `mergevar' "Matching result from merge"
end

* Run the merge with merge itself (strL variables are involved), passing the
* options merge accepts; cmerge-only options (threads, verbose,
* preserve_order) do not apply.
program define _cmerge_native, rclass
    args cmdline
    local 0 `"`cmdline'"'
    gettoken mtype 0 : 0, parse(" ")
    syntax [anything] using/ [, Keep(string) ASSert(string) GENerate(name) ///
        NOGENerate KEEPUSing(string) SORTED FORCE NOREPort Verbose NOLabel ///
        NONotes UPDATE REPLACE PRESERVE_order(integer 0) THReads(integer 0)]
    local opts `nogenerate' `sorted' `force' `noreport' `nolabel' `nonotes' `update' `replace'
    if "`keep'" != "" local opts `opts' keep(`keep')
    if "`assert'" != "" local opts `opts' assert(`assert')
    if "`generate'" != "" local opts `opts' generate(`generate')
    if "`keepusing'" != "" local opts `opts' keepusing(`keepusing')
    di as txt "(note: strL variables are involved; cmerge runs {bf:merge})"
    local caller = cond("$CTOOLS_cmerge_caller" != "", "$CTOOLS_cmerge_caller", "`c(stata_version)'")
    capture noisily version `caller': merge `mtype' `anything' using `"`using'"', `opts'
    if _rc {
        * merge exits 9 only when assert() fails; its result stays in memory
        if _rc == 9 global CTOOLS_cmerge_assert_failed 1
        exit _rc
    }
    local mvar = cond("`generate'" != "", "`generate'", "_merge")
    return scalar N = _N
    if "`nogenerate'" == "" {
        forvalues c = 1/5 {
            quietly count if `mvar' == `c'
            return scalar N_`c' = r(N)
        }
    }
    return local using `"`using'"'
    return local keyvars "`anything'"
end

* r(sorted) = 1 when the observations are in ascending order of varlist
* (Stata's sort order: strings bytewise, missing values above numbers)
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
