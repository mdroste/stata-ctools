*! version 0.9.0 26Jan2026
*! cmerge: Optimized C-accelerated merge for Stata datasets
*! Part of the ctools suite
*!
*! Description:
*!   High-performance merge with MINIMAL data transfer:
*!   - Only loads key + keepusing variables into C memory
*!   - Master non-key variables are STREAMED (never fully loaded)
*!   - Memory: O(nobs * (nkeys + nkeepusing)) instead of O(nobs * all_vars)
*!
*! Syntax:
*!   cmerge 1:1 varlist using filename [, options]
*!   cmerge 1:1 _n using filename [, options]
*!   cmerge m:1 varlist using filename [, options]
*!   cmerge 1:m varlist using filename [, options]
*!   cmerge m:m varlist using filename [, options]

* Mata helpers for optimized operations
capture mata: mata drop _cmerge_addvars()
capture mata: mata drop _cmerge_shared_flags()
mata:
void _cmerge_addvars(string rowvector types, string rowvector names)
{
    real scalar i, n
    n = cols(names)
    for (i = 1; i <= n; i++) {
        (void) st_addvar(types[i], names[i])
    }
}

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
    version 14.1

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
            TIMEit ///
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
            TIMEit ///
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

    * Validate keep() option values upfront (before any data manipulation)
    if "`keep'" != "" {
        foreach k of local keep {
            if "`k'" != "match" & "`k'" != "matched" & "`k'" != "3" ///
             & "`k'" != "master" & "`k'" != "1" ///
             & "`k'" != "using" & "`k'" != "2" {
                di as error "cmerge: invalid keep() value: `k'"
                di as error "  valid values: match (3), master (1), using (2)"
                exit 198
            }
        }
    }

    * Validate assert() option values upfront (before any data manipulation)
    if "`assert'" != "" {
        foreach a of local assert {
            if "`a'" != "match" & "`a'" != "matched" & "`a'" != "3" ///
             & "`a'" != "master" & "`a'" != "1" ///
             & "`a'" != "using" & "`a'" != "2" {
                di as error "cmerge: invalid assert() value: `a'"
                di as error "  valid values: match (3), master (1), using (2)"
                exit 198
            }
        }
    }

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
    local __do_timing = ("`verbose'" != "" | "`timeit'" != "")
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
    capture program list ctools_plugin
    if _rc != 0 {
        local __os = c(os)
        local __machine = c(machine_type)
        local __plugin = ""
        if "`__os'" == "Windows" {
            local __plugin "ctools_windows.plugin"
        }
        else if "`__os'" == "MacOSX" | strpos(lower("`__machine'"), "mac") > 0 {
            * Check cached architecture first (avoids subprocess on repeat calls)
            local __is_arm = 0
            if "$CTOOLS_ARCH_CACHE" == "arm64" {
                local __is_arm = 1
            }
            else if "$CTOOLS_ARCH_CACHE" == "x86_64" {
                local __is_arm = 0
            }
            else {
                * First call - detect and cache
                if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
                    local __is_arm = 1
                    global CTOOLS_ARCH_CACHE "arm64"
                }
                else {
                    tempfile __archfile
                    quietly shell uname -m > "`__archfile'" 2>&1
                    tempname __fh
                    file open `__fh' using "`__archfile'", read text
                    file read `__fh' __archline
                    file close `__fh'
                    capture erase "`__archfile'"
                    if strpos("`__archline'", "arm64") > 0 {
                        local __is_arm = 1
                        global CTOOLS_ARCH_CACHE "arm64"
                    }
                    else {
                        global CTOOLS_ARCH_CACHE "x86_64"
                    }
                }
            }
            local __plugin = cond(`__is_arm', "ctools_mac_arm.plugin", "ctools_mac_x86.plugin")
        }
        else if "`__os'" == "Unix" {
            local __plugin "ctools_linux.plugin"
        }
        else {
            local __plugin "ctools.plugin"
        }
        capture program ctools_plugin, plugin using("`__plugin'")
        if _rc != 0 & _rc != 110 & "`__plugin'" != "ctools.plugin" {
            capture program ctools_plugin, plugin using("ctools.plugin")
        }
        if _rc != 0 & _rc != 110 {
            di as error "cmerge: Could not load ctools plugin"
            exit 601
        }
    }

    * =========================================================================
    * Phase 1: Load using dataset - ONLY key + keepusing vars
    * =========================================================================

    if "`verbose'" != "" {
        di as text "Phase 1: Loading using dataset (keys + keepusing only)..."
    }

    * Use frames (Stata 16+) for faster context switching vs preserve/restore
    * Frames avoid copying the entire master dataset twice
    local __use_frames = 0
    if c(stata_version) >= 16 {
        capture frame create _cmerge_using
        if _rc == 0 {
            local __use_frames = 1
        }
    }

    if `__use_frames' {
        * Load using data into separate frame (master stays in default)
        frame _cmerge_using: qui use `"`using'"', clear
        frame change _cmerge_using
    }
    else {
        * Fallback: preserve/restore for Stata < 16
        preserve
        qui use `"`using'"', clear
    }

    * Check key variables exist in using (skip for _n merge)
    if !`merge_by_n' {
        foreach var of local keyvars {
            capture confirm variable `var'
            if _rc {
                di as error "cmerge: key variable `var' not found in using dataset"
                if `__use_frames' {
                    frame change default
                    frame drop _cmerge_using
                }
                else {
                    restore
                }
                exit 111
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
                if "`force'" == "" {
                    di as error "cmerge: key variable `var' is `master_type' in master but `using_type' in using"
                    di as error "use force option to override"
                    if `__use_frames' {
                        frame change default
                        frame drop _cmerge_using
                    }
                    else {
                        restore
                    }
                    exit 106
                }
            }
            local ++i
        }
    }

    * Build list of variables to keep: keys + keepusing (no keys for _n merge)
    local using_keep_vars "`keyvars'"
    local keepusing_count = 0
    local keepusing_names ""
    local keepusing_types ""

    if "`keepusing'" != "" {
        foreach var of local keepusing {
            capture confirm variable `var'
            if _rc {
                di as error "cmerge: keepusing variable `var' not found in using dataset"
                if `__use_frames' {
                    frame change default
                    frame drop _cmerge_using
                }
                else {
                    restore
                }
                exit 111
            }
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

    * Capture variable notes from using dataset (unless nonotes specified)
    if "`nonotes'" == "" {
        foreach var of local keepusing_names {
            local `var'_note : char `var'[note1]
        }
    }

    * Identify shared variables using Mata (O(n) vs O(n*m) nested loops)
    * A variable is shared if it exists in both master and using
    local shared_var_flags ""
    if `keepusing_count' > 0 {
        mata: _cmerge_shared_flags(tokens(st_local("keepusing_names")), tokens(st_local("master_varlist")))
    }

    * Capture value labels from using dataset (unless nolabel specified)
    local using_label_list ""
    local using_label_varmap ""
    if "`nolabel'" == "" {
        qui label dir
        local using_label_list "`r(names)'"

        foreach var of local keepusing_names {
            local vallbl : value label `var'
            if "`vallbl'" != "" {
                local using_label_varmap "`using_label_varmap' `var':`vallbl'"
            }
        }

        if "`using_label_list'" != "" {
            tempfile using_labels_file
            qui label save `using_label_list' using "`using_labels_file'", replace
        }
    }

    * Keep ONLY the necessary variables
    qui keep `using_keep_vars'

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
        * Save using data if needed
        tempfile using_saved
        if `using_nobs' > 0 {
            qui save `using_saved', replace
        }

        if `__use_frames' {
            frame change default
            frame drop _cmerge_using
        }
        else {
            restore
        }

        * Handle empty master: append using with _merge=2
        if `master_nobs' == 0 & `using_nobs' > 0 {
            qui append using `using_saved'
            if "`nogenerate'" == "" {
                qui gen byte `generate' = 2
            }
            local merge_master = 0
            local merge_using = `using_nobs'
            local merge_matched = 0
        }
        * Handle empty using: keep master with _merge=1
        else if `master_nobs' > 0 & `using_nobs' == 0 {
            * Add placeholder keepusing variables
            foreach vname of local keepusing_names {
                local vtype : word 1 of `keepusing_types'
                local keepusing_types : list keepusing_types - vtype
                capture confirm variable `vname'
                if _rc {
                    if substr("`vtype'", 1, 3) == "str" {
                        qui gen `vtype' `vname' = ""
                    }
                    else {
                        qui gen `vtype' `vname' = .
                    }
                }
            }
            if "`nogenerate'" == "" {
                qui gen byte `generate' = 1
            }
            local merge_master = `master_nobs'
            local merge_using = 0
            local merge_matched = 0
        }
        * Handle both empty
        else {
            if "`nogenerate'" == "" {
                qui gen byte `generate' = .
            }
            local merge_master = 0
            local merge_using = 0
            local merge_matched = 0
        }

        * Handle keep() option for empty datasets
        if "`keep'" != "" {
            local keep_codes ""
            foreach k of local keep {
                if "`k'" == "match" | "`k'" == "matched" | "`k'" == "3" {
                    local keep_codes "`keep_codes' 3"
                }
                else if "`k'" == "master" | "`k'" == "1" {
                    local keep_codes "`keep_codes' 1"
                }
                else if "`k'" == "using" | "`k'" == "2" {
                    local keep_codes "`keep_codes' 2"
                }
            }
            if "`nogenerate'" == "" {
                local keep_expr ""
                foreach code of local keep_codes {
                    if "`keep_expr'" == "" {
                        local keep_expr "`generate' == `code'"
                    }
                    else {
                        local keep_expr "`keep_expr' | `generate' == `code'"
                    }
                }
                qui keep if `keep_expr'
            }
        }

        * Handle assert() option for empty datasets
        if "`assert'" != "" {
            local allow_master = 0
            local allow_using = 0
            local allow_match = 0
            foreach a of local assert {
                if "`a'" == "match" | "`a'" == "matched" | "`a'" == "3" {
                    local allow_match = 1
                }
                else if "`a'" == "master" | "`a'" == "1" {
                    local allow_master = 1
                }
                else if "`a'" == "using" | "`a'" == "2" {
                    local allow_using = 1
                }
            }
            local assert_failed = 0
            if !`allow_master' & `merge_master' > 0 {
                local assert_failed = 1
            }
            if !`allow_using' & `merge_using' > 0 {
                local assert_failed = 1
            }
            if !`allow_match' & `merge_matched' > 0 {
                local assert_failed = 1
            }
            if `assert_failed' {
                di as error "cmerge: assertion failed"
                exit 9
            }
        }

        * Display merge table for empty datasets
        if "`noreport'" == "" {
            local disp_not_matched = `merge_master' + `merge_using'
            di as text ""
            di as text "    Result" _col(33) "Number of obs"
            di as text "    {hline 41}"
            di as text "    Not matched" _col(33) as result %13.0fc `disp_not_matched'
            if `disp_not_matched' > 0 {
                di as text "        from master" _col(33) as result %13.0fc `merge_master' as text "  (_merge==1)"
                di as text "        from using" _col(33) as result %13.0fc `merge_using' as text "  (_merge==2)"
            }
            di as text ""
            di as text "    Matched" _col(33) as result %13.0fc `merge_matched' as text "  (_merge==3)"
            di as text "    {hline 41}"
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
    * Pass sorted flag if user specified OR auto-detected
    if "`sorted'" != "" | `using_sorted' {
        local plugin_args "`plugin_args' sorted"
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
            frame change default
            frame drop _cmerge_using
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
        timer clear 99
        timer on 97   /* create vars template */
    }

    * Build lists of new variables to create (excluding existing shared vars)
    * Use shared_var_flags (already computed) to determine which vars need creation
    local new_var_names "_cmerge_orig_row"
    local new_var_types "long"
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

    * Create empty template dataset with new variable definitions
    * Using data is already cached in plugin, so we can clear and reuse memory
    * Benchmark shows template+append is 3.5x faster than Mata st_addvar
    tempfile empty_vars_template
    if `n_new_vars' > 0 {
        clear
        qui set obs 1
        local var_idx = 1
        foreach vtype of local new_var_types {
            local vname : word `var_idx' of `new_var_names'
            if substr("`vtype'", 1, 3) == "str" {
                qui gen `vtype' `vname' = ""
            }
            else {
                qui gen `vtype' `vname' = .
            }
            local ++var_idx
        }
        qui drop in 1
        qui save `empty_vars_template', emptyok replace
    }

    if `__do_timing' {
        timer off 97
        timer on 96   /* frame switch / restore */
    }

    if `__use_frames' {
        frame change default
        frame drop _cmerge_using
    }
    else {
        restore
    }

    if `__do_timing' {
        timer off 96
        timer on 97   /* append vars */
    }

    * Append empty template to add variable definitions (faster than st_addvar)
    if `n_new_vars' > 0 {
        qui append using `empty_vars_template'
    }

    if `__do_timing' {
        timer off 97
        timer on 99   /* store row index */
    }

    * Fill _cmerge_orig_row with row numbers
    local orig_row_idx = c(k) - `n_new_vars' + 1
    mata: st_store(., `orig_row_idx', (1::st_nobs()))

    if `__do_timing' {
        timer off 99
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
    * IMPORTANT: Use current c(k) not original master_nvars, because we added
    * _cmerge_orig_row and other new variables that need to be loaded
    local current_nvars = c(k)
    local plugin_args "execute `merge_code' `nkeys' `master_key_indices'"
    local plugin_args "`plugin_args' orig_row_idx `orig_row_idx'"
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
    * Pass sorted flag if user specified OR auto-detected for master
    if "`sorted'" != "" | `master_sorted' {
        local plugin_args "`plugin_args' sorted"
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

    * Call plugin Phase 2 (current_varlist already computed)
    capture noisily plugin call ctools_plugin `current_varlist', "cmerge `threads_code' `plugin_args'"
    local plugin_rc = _rc

    * End plugin2 timer, start post-plugin timer
    if `__do_timing' {
        timer off 94
        timer on 95   /* Post-plugin */
    }

    if `plugin_rc' {
        di as error "cmerge: C plugin merge failed (error `plugin_rc')"
        capture drop _cmerge_orig_row
        exit `plugin_rc'
    }

    * Get merge results from scalars
    local merge_master = scalar(_cmerge_N_1)
    local merge_using = scalar(_cmerge_N_2)
    local merge_matched = scalar(_cmerge_N_3)
    local total_obs = scalar(_cmerge_N)

    * Drop temporary row index variable
    qui capture drop _cmerge_orig_row

    * Trim excess observations (drop in is 3x faster than keep in)
    if `total_obs' < _N {
        local __drop_start = `total_obs' + 1
        qui drop in `__drop_start'/l
    }

    * Apply variable notes from using dataset (unless nonotes specified)
    if "`nonotes'" == "" {
        foreach var of local keepusing_names {
            if `"``var'_note'"' != "" {
                char `var'[note1] `"``var'_note'"'
            }
        }
    }

    * Apply value labels from using dataset (unless nolabel specified)
    if "`nolabel'" == "" & "`using_label_list'" != "" {
        qui do "`using_labels_file'"

        foreach mapping of local using_label_varmap {
            gettoken vname mapping : mapping, parse(":")
            local labelname = substr("`mapping'", 2, .)
            capture label values `vname' `labelname'
        }
    }

    * =========================================================================
    * Phase 3: Apply keep/assert options and display results
    * =========================================================================

    * Handle keep option (already validated upfront)
    if "`keep'" != "" {
        local keep_codes ""
        foreach k of local keep {
            if "`k'" == "match" | "`k'" == "matched" | "`k'" == "3" {
                local keep_codes "`keep_codes' 3"
            }
            else if "`k'" == "master" | "`k'" == "1" {
                local keep_codes "`keep_codes' 1"
            }
            else if "`k'" == "using" | "`k'" == "2" {
                local keep_codes "`keep_codes' 2"
            }
        }

        * Use the appropriate merge variable for filtering
        local filter_var "`generate'"
        if "`temp_merge_var'" != "" {
            local filter_var "`temp_merge_var'"
        }

        local keep_expr ""
        foreach code of local keep_codes {
            if "`keep_expr'" == "" {
                local keep_expr "`filter_var' == `code'"
            }
            else {
                local keep_expr "`keep_expr' | `filter_var' == `code'"
            }
        }
        qui keep if `keep_expr'

        * Drop temp filter variable if we created one
        if "`temp_merge_var'" != "" {
            qui drop `temp_merge_var'
        }
    }

    * Handle assert option
    * assert(X Y Z) means all observations must have _merge in {X, Y, Z}
    * I.e., observations NOT in the allowed set cause failure
    if "`assert'" != "" {
        local allow_master = 0
        local allow_using = 0
        local allow_match = 0

        foreach a of local assert {
            if "`a'" == "match" | "`a'" == "matched" | "`a'" == "3" {
                local allow_match = 1
            }
            else if "`a'" == "master" | "`a'" == "1" {
                local allow_master = 1
            }
            else if "`a'" == "using" | "`a'" == "2" {
                local allow_using = 1
            }
        }

        * Check that no disallowed categories exist
        local assert_failed = 0
        if !`allow_master' & `merge_master' > 0 {
            local assert_failed = 1
        }
        if !`allow_using' & `merge_using' > 0 {
            local assert_failed = 1
        }
        if !`allow_match' & `merge_matched' > 0 {
            local assert_failed = 1
        }

        if `assert_failed' {
            di as error "cmerge: assertion failed"
            exit 9
        }
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
        quietly timer list 99
        local __time_storerow = r(t99)
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
        * Calculate display counts based on keep() option
        local disp_master = `merge_master'
        local disp_using = `merge_using'
        local disp_matched = `merge_matched'

        if "`keep'" != "" {
            * Reset counts to 0 for categories not kept
            local disp_master = 0
            local disp_using = 0
            local disp_matched = 0
            foreach code of local keep_codes {
                if `code' == 1 local disp_master = `merge_master'
                if `code' == 2 local disp_using = `merge_using'
                if `code' == 3 local disp_matched = `merge_matched'
            }
        }

        local disp_not_matched = `disp_master' + `disp_using'

        di as text ""
        di as text "    Result" _col(33) "Number of obs"
        di as text "    {hline 41}"
        di as text "    Not matched" _col(33) as result %13.0fc `disp_not_matched'
        * Only show detail rows if there are unmatched observations
        if `disp_not_matched' > 0 {
            di as text "        from master" _col(33) as result %13.0fc `disp_master' as text "  (_merge==1)"
            di as text "        from using" _col(33) as result %13.0fc `disp_using' as text "  (_merge==2)"
        }
        di as text ""
        di as text "    Matched" _col(33) as result %13.0fc `disp_matched' as text "  (_merge==3)"
        di as text "    {hline 41}"
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
        timer clear 99
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
    return scalar time = `elapsed'
    return local using `"`using'"'
    if `merge_by_n' {
        return local keyvars "_n"
    }
    else {
        return local keyvars "`keyvars'"
    }

end
