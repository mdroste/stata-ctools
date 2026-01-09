*! version 2.1.0 05Jan2026
*! cmerge: C-accelerated merge for Stata datasets
*! Part of the ctools suite
*!
*! Description:
*!   High-performance drop-in replacement for Stata's merge command using a
*!   C plugin with parallel data loading, radix sort, and parallel merge.
*!   All merge operations are performed entirely in C for maximum speed.
*!
*! Syntax:
*!   cmerge 1:1 varlist using filename [, options]
*!   cmerge m:1 varlist using filename [, options]
*!   cmerge 1:m varlist using filename [, options]
*!   cmerge m:m varlist using filename [, options]
*!
*! Options:
*!   keep(match|master|using|1 2 3)  - Which observations to keep
*!   assert(match|master|using|1 2 3) - Assert merge results
*!   generate(varname)               - Name for _merge variable (default: _merge)
*!   nogenerate                      - Do not create _merge variable
*!   keepusing(varlist)              - Variables to keep from using dataset
*!   sorted                          - Assert data is already sorted
*!   force                           - Force string/numeric mismatches
*!   noreport                        - Do not display merge result table
*!   verbose                         - Display timing and progress info
*!   nolabel                         - Do not copy value labels from using

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
        * Try adding .dta extension
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

    * Check if master is already sorted on key variables
    local master_sorted = 0
    if "`sorted'" != "" {
        local master_sorted = 1
    }
    else {
        * Check sortedby characteristic
        local sortedby : sortedby
        if "`sortedby'" != "" {
            local is_sorted = 1
            local i = 1
            foreach var of local keyvars {
                local sorted_var : word `i' of `sortedby'
                if "`sorted_var'" != "`var'" {
                    local is_sorted = 0
                    continue, break
                }
                local ++i
            }
            if `is_sorted' {
                local master_sorted = 1
            }
        }
    }

    * Store master dataset info
    local master_nobs = _N
    local master_nvars = c(k)

    if `master_nobs' == 0 {
        di as error "cmerge: no observations in master dataset"
        exit 2000
    }

    * Get key variable indices (1-based for plugin)
    local keyvar_indices ""
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
        local keyvar_indices "`keyvar_indices' `idx'"
    }

    * Display header if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cmerge: C-Accelerated Merge (Full Native)"
        di as text "{hline 60}"
        di as text "Merge type:  " as result "`merge_type'"
        di as text "Key vars:    " as result "`keyvars'"
        di as text "Using file:  " as result `"`using'"'
        di as text "Master obs:  " as result `master_nobs'
        di as text "Master vars: " as result `master_nvars'
        di as text "Sorted:      " as result cond(`master_sorted', "yes", "no")
        di as text "{hline 60}"
        di ""
    }

    * Record start time
    timer clear 99
    timer on 99

    * Load the platform-appropriate ctools plugin if not already loaded
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
    * Phase 1: Load using dataset into C plugin memory
    * =========================================================================

    if "`verbose'" != "" {
        di as text "Phase 1: Loading using dataset into C..."
    }

    * Preserve master dataset
    preserve

    * Load using dataset into Stata
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

    * Get using dataset info
    local using_nobs = _N
    local using_nvars = c(k)

    if "`verbose'" != "" {
        di as text "  Using dataset: " as result `using_nobs' as text " obs, " as result `using_nvars' as text " vars"
    }

    * Check if using is sorted on key variables
    local using_sorted = 0
    if "`sorted'" != "" {
        local using_sorted = 1
    }
    else {
        local sortedby : sortedby
        if "`sortedby'" != "" {
            local is_sorted = 1
            local i = 1
            foreach var of local keyvars {
                local sorted_var : word `i' of `sortedby'
                if "`sorted_var'" != "`var'" {
                    local is_sorted = 0
                    continue, break
                }
                local ++i
            }
            if `is_sorted' {
                local using_sorted = 1
            }
        }
    }

    * Get using key variable indices (1-based)
    local using_keyvar_indices ""
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
        local using_keyvar_indices "`using_keyvar_indices' `idx'"
    }

    * Handle keepusing option - keep only key vars and requested vars
    if "`keepusing'" != "" {
        local keepvars "`keyvars' `keepusing'"
        qui keep `keepvars'

        * Recalculate using key indices after keep
        local using_keyvar_indices ""
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
            local using_keyvar_indices "`using_keyvar_indices' `idx'"
        }
    }

    * Check if keep option excludes using-only rows (allows hash merge optimization)
    * Hash merge can only be used when we don't need to include using-only rows
    * This is the case when keep(1 3), keep(match master), or similar is specified
    local keep_excludes_using = 0
    if "`keep'" != "" {
        local keep_excludes_using = 1
        foreach k of local keep {
            if "`k'" == "using" | "`k'" == "2" {
                local keep_excludes_using = 0
            }
        }
    }

    * Build plugin command for load_using
    local plugin_args "load_using `nkeys' `using_keyvar_indices'"
    if `using_sorted' {
        local plugin_args "`plugin_args' sorted"
    }
    if "`verbose'" != "" {
        local plugin_args "`plugin_args' verbose"
    }
    * Pass merge type so plugin can skip sorting for m:1/1:1 (hash-based lookup)
    local plugin_args "`plugin_args' mergetype `merge_code'"
    * Pass keep_excludes_using flag so plugin knows if hash merge will be used
    if `keep_excludes_using' {
        local plugin_args "`plugin_args' keep_excludes_using"
    }

    * Get all variables in using dataset for plugin call
    qui ds
    local using_allvars `r(varlist)'

    * Build list of non-key using variable names and types (in order)
    * These will be used to create correctly-typed placeholder variables
    local using_nonkey_varnames ""
    local using_nonkey_vartypes ""
    foreach v of local using_allvars {
        local is_key = 0
        foreach k of local keyvars {
            if "`v'" == "`k'" {
                local is_key = 1
            }
        }
        if !`is_key' {
            local using_nonkey_varnames "`using_nonkey_varnames' `v'"
            * Determine variable type
            capture confirm string variable `v'
            if _rc == 0 {
                * It's a string - get the storage type
                local vtype : type `v'
                local using_nonkey_vartypes "`using_nonkey_vartypes' `vtype'"
            }
            else {
                local using_nonkey_vartypes "`using_nonkey_vartypes' double"
            }
        }
    }

    * Call plugin to load using data into C memory
    * IMPORTANT: Must pass varlist to plugin so it can access the data
    capture noisily plugin call ctools_plugin `using_allvars', "cmerge `plugin_args'"
    local plugin_rc = _rc

    if `plugin_rc' {
        di as error "cmerge: failed to load using data into C plugin (error `plugin_rc')"
        restore
        exit `plugin_rc'
    }

    * =========================================================================
    * Phase 2: Restore master and execute merge in C
    * =========================================================================

    if "`verbose'" != "" {
        di as text ""
        di as text "Phase 2: Executing merge in C..."
    }

    restore

    * We need to prepare Stata's dataset to receive the output
    * The plugin will resize and populate it

    * Get the number of non-key variables from using
    local using_nonkey_nvars = `using_nvars' - `nkeys'
    if "`keepusing'" != "" {
        local keepusing_count : word count `keepusing'
        local using_nonkey_nvars = `keepusing_count'
    }

    * Calculate expected output variables
    local output_nvars = `master_nvars' + `using_nonkey_nvars'
    if "`nogenerate'" == "" {
        local ++output_nvars
    }

    * Add placeholder variables for using non-key vars and _merge
    * These will be overwritten by the plugin
    * Create placeholders with correct types matching using variables
    local v = 1
    foreach vtype of local using_nonkey_vartypes {
        if substr("`vtype'", 1, 3) == "str" {
            * String variable - create with correct length
            capture gen `vtype' __using_temp`v' = ""
        }
        else {
            * Numeric variable
            capture gen double __using_temp`v' = .
        }
        local ++v
    }
    if "`nogenerate'" == "" {
        capture gen byte `generate' = .
    }

    * For merges that may include using-only rows, expand dataset to accommodate
    * Maximum output size = master_nobs + using_nobs (if all using rows are unmatched)
    * Always expand - the C plugin produces all rows, keep() filtering happens after
    local expanded = 0
    local max_output_nobs = `master_nobs' + `using_nobs'
    if `max_output_nobs' > `master_nobs' {
        qui set obs `max_output_nobs'
        local expanded = 1
    }

    * Build plugin command for execute
    * Pass original master_nobs and master_nvars so plugin knows real data sizes
    local plugin_args "execute `merge_code' `nkeys' `keyvar_indices'"
    if `master_sorted' {
        local plugin_args "`plugin_args' sorted"
    }
    if "`verbose'" != "" {
        local plugin_args "`plugin_args' verbose"
    }
    if "`nogenerate'" != "" {
        local plugin_args "`plugin_args' nogenerate"
    }
    if `keep_excludes_using' {
        local plugin_args "`plugin_args' keep_excludes_using"
    }
    * Pass original master_nobs if we expanded
    if `expanded' {
        local plugin_args "`plugin_args' master_nobs `master_nobs'"
    }
    * Always pass original master_nvars (before adding placeholders)
    local plugin_args "`plugin_args' master_nvars `master_nvars'"

    * Get all variables in master dataset for plugin call
    qui ds
    local master_allvars `r(varlist)'

    * Call plugin to execute merge
    * IMPORTANT: Must pass varlist to plugin so it can access the data
    capture noisily plugin call ctools_plugin `master_allvars', "cmerge `plugin_args'"
    local plugin_rc = _rc
    if `plugin_rc' {
        di as error "cmerge: C plugin merge failed (error `plugin_rc')"
        * Cleanup temp vars
        capture drop __using_temp*
        exit `plugin_rc'
    }

    * Get merge results from scalars saved by plugin
    local merge_master = scalar(_cmerge_N_1)
    local merge_using = scalar(_cmerge_N_2)
    local merge_matched = scalar(_cmerge_N_3)
    local total_obs = scalar(_cmerge_N)

    * Rename placeholder variables to actual using variable names
    local v = 1
    foreach varname of local using_nonkey_varnames {
        * Check if variable already exists in master (skip if so)
        capture confirm variable `varname'
        if _rc != 0 {
            * Variable doesn't exist, rename placeholder to this name
            capture rename __using_temp`v' `varname'
        }
        else {
            * Variable exists in master, just drop the placeholder
            * (master values are kept, like native merge behavior)
            capture drop __using_temp`v'
        }
        local ++v
    }
    * Drop any remaining placeholders (shouldn't happen but just in case)
    capture drop __using_temp*

    * Trim excess observations if we over-expanded
    * (happens when some using rows matched and weren't added as using-only)
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

    * Display merge table if not suppressed
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

    * Display timing if requested
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
