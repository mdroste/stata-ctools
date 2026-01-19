*! version 1.0.0 17Jan2026
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

program define cmerge, rclass
    version 14.0

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
        exit 198
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

    if `master_nobs' == 0 {
        di as error "cmerge: no observations in master dataset"
        exit 2000
    }

    * Get master key variable indices (1-based) - skip for _n merge
    local master_key_indices ""
    if !`merge_by_n' {
        local pos = 1
        foreach v of local master_varlist {
            foreach var of local keyvars {
                if "`v'" == "`var'" {
                    local master_key_indices "`master_key_indices' `pos'"
                }
            }
            local ++pos
        }
        * Verify all keys found
        local n_found : word count `master_key_indices'
        if `n_found' != `nkeys' {
            di as error "cmerge: not all key variables found in master"
            exit 111
        }
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

    * Load the platform-appropriate ctools plugin
    capture program list ctools_plugin
    if _rc != 0 {
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
            if `__is_arm' == 0 {
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

    preserve

    * Load using dataset
    qui use `"`using'"', clear

    * Check key variables exist in using (skip for _n merge)
    if !`merge_by_n' {
        foreach var of local keyvars {
            capture confirm variable `var'
            if _rc {
                di as error "cmerge: key variable `var' not found in using dataset"
                restore
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
                    restore
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
                restore
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

    * Identify shared variables (exist in both master and using)
    local shared_var_flags ""
    foreach vname of local keepusing_names {
        local is_shared = 0
        foreach mvar of local master_varlist {
            if "`mvar'" == "`vname'" {
                local is_shared = 1
                continue, break
            }
        }
        local shared_var_flags "`shared_var_flags' `is_shared'"
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

    if "`verbose'" != "" {
        if `merge_by_n' {
            di as text "  Using: " as result `using_nobs' as text " obs, keeping " as result `keepusing_count' as text " keepusing vars (merge by _n)"
        }
        else {
            di as text "  Using: " as result `using_nobs' as text " obs, keeping " as result `nkeys' as text " keys + " as result `keepusing_count' as text " keepusing vars"
        }
    }

    * Get variable indices in the reduced using dataset - optimized single pass
    local using_key_indices ""
    local using_keepusing_indices ""
    local pos = 1
    unab using_varlist : _all
    foreach v of local using_varlist {
        local is_key = 0
        if !`merge_by_n' {
            foreach k of local keyvars {
                if "`v'" == "`k'" {
                    local using_key_indices "`using_key_indices' `pos'"
                    local is_key = 1
                }
            }
        }
        if !`is_key' {
            local using_keepusing_indices "`using_keepusing_indices' `pos'"
        }
        local ++pos
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
    if "`sorted'" != "" {
        local plugin_args "`plugin_args' sorted"
    }
    if `merge_by_n' {
        local plugin_args "`plugin_args' merge_by_n"
    }

    * End pre-plugin1 timer, start plugin1 timer
    if `__do_timing' {
        timer off 91
        timer on 92   /* Plugin Phase 1 */
    }

    * Call plugin Phase 1 with reduced varlist (using_varlist already computed)
    capture noisily plugin call ctools_plugin `using_varlist', "cmerge `plugin_args'"
    local plugin_rc = _rc

    * End plugin1 timer
    if `__do_timing' {
        timer off 92
    }

    if `plugin_rc' {
        di as error "cmerge: failed to load using data (error `plugin_rc')"
        restore
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

    * Create empty template dataset with new variable definitions (FAST append trick)
    * This is 4.5x faster than Java/Mata because it only adds metadata, not data
    local n_new_vars : word count `new_var_names'
    tempfile empty_vars_template
    if `n_new_vars' > 0 {
        clear
        qui set obs 1
        local var_idx = 1
        foreach vtype of local new_var_types {
            local vname : word `var_idx' of `new_var_names'
            * String variables need empty string, numeric need missing
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
        timer on 96   /* restore */
    }

    restore

    if `__do_timing' {
        timer off 96
        timer on 97   /* append vars */
    }

    * Append empty template to add variable definitions (FAST - no data initialization)
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

    * Compute variable indices in a single pass
    local keepusing_placeholder_indices ""
    unab current_varlist : _all
    foreach vname of local keepusing_names {
        local pos = 1
        foreach v of local current_varlist {
            if "`v'" == "`vname'" {
                local keepusing_placeholder_indices "`keepusing_placeholder_indices' `pos'"
                continue, break
            }
            local ++pos
        }
    }

    * Get _merge variable index (already in current_varlist from unab)
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
    if "`sorted'" != "" {
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
    capture noisily plugin call ctools_plugin `current_varlist', "cmerge `plugin_args'"
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

    * Trim excess observations
    if `total_obs' < _N {
        qui keep in 1/`total_obs'
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
    if "`assert'" != "" {
        local assert_failed = 0
        foreach a of local assert {
            if "`a'" == "match" | "`a'" == "matched" | "`a'" == "3" {
                if `merge_master' > 0 | `merge_using' > 0 {
                    local assert_failed = 1
                }
            }
            else if "`a'" == "master" | "`a'" == "1" {
                if `merge_using' > 0 | `merge_matched' > 0 {
                    local assert_failed = 1
                }
            }
            else if "`a'" == "using" | "`a'" == "2" {
                if `merge_master' > 0 | `merge_matched' > 0 {
                    local assert_failed = 1
                }
            }
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
