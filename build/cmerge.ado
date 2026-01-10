*! version 3.0.0 10Jan2026
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

    * Parse key variables and using filename
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

    * Validate using file exists
    capture confirm file `"`using'"'
    if _rc {
        capture confirm file `"`using'.dta"'
        if _rc {
            di as error `"file `using' not found"'
            exit 601
        }
        local using `"`using'.dta"'
    }

    * Get variable types for key variables in master
    local master_keytypes ""
    foreach var of local keyvars {
        capture confirm string variable `var'
        if !_rc {
            local master_keytypes "`master_keytypes' str"
        }
        else {
            local master_keytypes "`master_keytypes' num"
        }
    }

    * Store master dataset info BEFORE any modifications
    local master_nobs = _N
    local master_nvars = c(k)
    qui ds
    local master_varlist "`r(varlist)'"

    if `master_nobs' == 0 {
        di as error "cmerge: no observations in master dataset"
        exit 2000
    }

    * Get master key variable indices (1-based)
    local master_key_indices ""
    foreach var of local keyvars {
        local idx = 0
        local pos = 1
        qui ds
        foreach v in `r(varlist)' {
            if "`v'" == "`var'" {
                local idx = `pos'
                continue, break
            }
            local ++pos
        }
        if `idx' == 0 {
            di as error "cmerge: key variable `var' not found"
            exit 111
        }
        local master_key_indices "`master_key_indices' `idx'"
    }

    * Display header if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cmerge: Optimized C-Accelerated Merge"
        di as text "{hline 60}"
        di as text "Merge type:  " as result "`merge_type'"
        di as text "Key vars:    " as result "`keyvars'"
        di as text "Using file:  " as result `"`using'"'
        di as text "Master obs:  " as result `master_nobs'
        di as text "Master vars: " as result `master_nvars'
        di as text "{hline 60}"
        di ""
    }

    timer clear 99
    timer on 99

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

    * Check key variables exist in using
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

    * Build list of variables to keep: keys + keepusing
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

            * Capture type for placeholder creation
            capture confirm string variable `var'
            if _rc == 0 {
                local vtype : type `var'
                local keepusing_types "`keepusing_types' `vtype'"
            }
            else {
                local keepusing_types "`keepusing_types' double"
            }
        }
    }
    else {
        * No keepusing specified - get all non-key variables from using
        * Include variables that exist in master (shared vars) - they'll be
        * written only for using-only rows to avoid overwriting master values
        qui ds
        local all_using_vars `r(varlist)'
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

                capture confirm string variable `var'
                if _rc == 0 {
                    local vtype : type `var'
                    local keepusing_types "`keepusing_types' `vtype'"
                }
                else {
                    local keepusing_types "`keepusing_types' double"
                }
            }
        }
        local using_keep_vars "`keyvars' `keepusing_names'"
    }

    * Keep ONLY the necessary variables
    qui keep `using_keep_vars'

    local using_nobs = _N

    if "`verbose'" != "" {
        di as text "  Using: " as result `using_nobs' as text " obs, keeping " as result `nkeys' as text " keys + " as result `keepusing_count' as text " keepusing vars"
    }

    * Get variable indices in the reduced using dataset
    local using_key_indices ""
    local using_keepusing_indices ""
    local pos = 1
    qui ds
    foreach v in `r(varlist)' {
        local is_key = 0
        foreach k of local keyvars {
            if "`v'" == "`k'" {
                local using_key_indices "`using_key_indices' `pos'"
                local is_key = 1
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

    * Call plugin Phase 1 with reduced varlist
    qui ds
    capture noisily plugin call ctools_plugin `r(varlist)', "cmerge `plugin_args'"
    local plugin_rc = _rc

    if `plugin_rc' {
        di as error "cmerge: failed to load using data (error `plugin_rc')"
        restore
        exit `plugin_rc'
    }

    * =========================================================================
    * Phase 2: Execute merge with streaming
    * =========================================================================

    if "`verbose'" != "" {
        di as text ""
        di as text "Phase 2: Executing merge with streaming..."
    }

    restore

    * Generate original row index BEFORE any modifications
    gen long _cmerge_orig_row = _n
    local orig_row_idx = c(k)

    * Create placeholder variables for keepusing (with correct types)
    local keepusing_placeholder_indices ""
    local placeholder_num = 1
    foreach vtype of local keepusing_types {
        local vname : word `placeholder_num' of `keepusing_names'

        * Check if variable already exists in master (overlapping var)
        capture confirm variable `vname'
        if _rc != 0 {
            * Variable doesn't exist - create placeholder
            if substr("`vtype'", 1, 3) == "str" {
                qui gen `vtype' `vname' = ""
            }
            else {
                qui gen double `vname' = .
            }
        }
        * Get current index (whether new or existing)
        local new_idx = 0
        local pos = 1
        qui ds
        foreach v in `r(varlist)' {
            if "`v'" == "`vname'" {
                local new_idx = `pos'
            }
            local ++pos
        }
        local keepusing_placeholder_indices "`keepusing_placeholder_indices' `new_idx'"
        local ++placeholder_num
    }

    * Create _merge variable
    local merge_var_idx = 0
    if "`nogenerate'" == "" {
        qui gen byte `generate' = .
        local merge_var_idx = c(k)
    }

    * Expand dataset to maximum possible output size
    local max_output = `master_nobs' + `using_nobs'
    if `max_output' > _N {
        qui set obs `max_output'
    }

    * Build plugin command for execute
    local plugin_args "execute `merge_code' `nkeys' `master_key_indices'"
    local plugin_args "`plugin_args' orig_row_idx `orig_row_idx'"
    local plugin_args "`plugin_args' master_nobs `master_nobs'"
    local plugin_args "`plugin_args' master_nvars `master_nvars'"
    local plugin_args "`plugin_args' n_keepusing `keepusing_count'"
    if `keepusing_count' > 0 {
        local plugin_args "`plugin_args' keepusing_placeholders `keepusing_placeholder_indices'"
    }
    if "`nogenerate'" == "" {
        local plugin_args "`plugin_args' merge_var_idx `merge_var_idx'"
    }
    local plugin_args "`plugin_args' preserve_order `preserve_order'"
    if "`verbose'" != "" {
        local plugin_args "`plugin_args' verbose"
    }

    * Call plugin Phase 2
    qui ds
    capture noisily plugin call ctools_plugin `r(varlist)', "cmerge `plugin_args'"
    local plugin_rc = _rc

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
    capture drop _cmerge_orig_row

    * Trim excess observations
    if `total_obs' < _N {
        qui keep in 1/`total_obs'
    }

    * =========================================================================
    * Phase 3: Apply keep/assert options and display results
    * =========================================================================

    * Handle keep option
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
            else {
                di as error "cmerge: invalid keep() value: `k'"
                exit 198
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

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Display merge table
    if "`noreport'" == "" {
        di as text ""
        di as text "    Result" _col(33) "Number of obs"
        di as text "    {hline 41}"
        di as text "    Not matched" _col(33) as result %13.0fc `=`merge_master' + `merge_using''
        if `merge_master' > 0 {
            di as text "        from master" _col(33) as result %13.0fc `merge_master' as text "  (_merge==1)"
        }
        if `merge_using' > 0 {
            di as text "        from using" _col(33) as result %13.0fc `merge_using' as text "  (_merge==2)"
        }
        di as text ""
        di as text "    Matched" _col(33) as result %13.0fc `merge_matched' as text "  (_merge==3)"
        di as text "    {hline 41}"
    }

    * Display timing
    if "`timeit'" != "" | "`verbose'" != "" {
        di as text ""
        di as text "Total time: " as result %9.3f `elapsed' as text " seconds"
    }

    timer clear 99

    * Return results
    return scalar N = _N
    return scalar N_1 = `merge_master'
    return scalar N_2 = `merge_using'
    return scalar N_3 = `merge_matched'
    return scalar time = `elapsed'
    return local using `"`using'"'
    return local keyvars "`keyvars'"

end
