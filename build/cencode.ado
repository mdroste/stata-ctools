*! version 0.9.0 26Jan2026
*! cencode: C-accelerated string encoding for Stata
*! Part of the ctools suite
*!
*! Replaces Stata's built-in encode command with a high-performance
*! C implementation featuring:
*!   - Parallel unique string collection
*!   - Lock-free parallel encoding
*!   - Efficient string interning with arena allocator
*!
*! Syntax: cencode varlist [if] [in], generate(newvarlist) | replace [label(name) noextend Verbose THReads(integer)]
*!
*! The cencode command inherits the same syntax and functionality as Stata's
*! built-in encode command, producing identical results but with better
*! performance on large datasets.
*!
*! ctools extensions (not in Stata's encode):
*!   - varlist support: encode multiple string variables at once
*!   - replace option: replace original variables with encoded versions

program define cencode, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    syntax varlist(string) [if] [in], [Generate(string) replace Label(name) noextend Verbose THReads(integer 0)]

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Check that either generate or replace is specified, but not both
    if "`generate'" == "" & "`replace'" == "" {
        di as error "cencode: must specify either generate() or replace option"
        exit 198
    }
    if "`generate'" != "" & "`replace'" != "" {
        di as error "cencode: cannot specify both generate() and replace options"
        exit 198
    }

    * Count input variables
    local n_vars : word count `varlist'

    * Handle replace vs generate
    local __do_replace = 0
    if "`replace'" != "" {
        local __do_replace = 1
    }
    else {
        * Check generate list has same number of variables
        local n_gen : word count `generate'
        if `n_gen' != `n_vars' {
            di as error "cencode: generate() must specify `n_vars' variable(s) to match varlist"
            exit 198
        }
        * Check that none of the generate variables already exist
        forvalues i = 1/`n_gen' {
            local gvar : word `i' of `generate'
            capture confirm variable `gvar'
            if _rc == 0 {
                di as error "cencode: variable `gvar' already exists"
                exit 110
            }
        }
    }

    * Check that all source variables are string
    foreach v of local varlist {
        capture confirm string variable `v'
        if _rc != 0 {
            di as error "cencode: `v' is not a string variable"
            exit 198
        }
    }

    * =========================================================================
    * END VALIDATION
    * =========================================================================

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        timer clear 90
        timer on 90
    }

    * Mark sample
    marksample touse, strok

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
            di as error "cencode: Could not load ctools plugin"
            exit 601
        }
    }

    * Build threads option
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * =========================================================================
    * Process each variable
    * =========================================================================

    forvalues i = 1/`n_vars' {
        local srcvar : word `i' of `varlist'

        * Determine destination variable name
        if `__do_replace' {
            tempvar __tempgen
            local destvar "`__tempgen'"
            local __final_name "`srcvar'"
        }
        else {
            local destvar : word `i' of `generate'
            local __final_name "`destvar'"
        }

        * Set label name for this variable
        if "`label'" != "" {
            local this_label "`label'"
        }
        else {
            local this_label "`__final_name'"
        }

        * Get variable index for source
        unab allvars : *
        local var_idx = 0
        local idx = 1
        foreach v of local allvars {
            if ("`v'" == "`srcvar'") {
                local var_idx = `idx'
                continue, break
            }
            local ++idx
        }

        if `var_idx' == 0 {
            di as error "cencode: could not find variable `srcvar'"
            exit 111
        }

        * Create the destination variable (numeric, long type for sufficient range)
        quietly generate long `destvar' = .

        * Get index of new variable
        unab allvars : *
        local gen_idx = 0
        local idx = 1
        foreach v of local allvars {
            if ("`v'" == "`destvar'") {
                local gen_idx = `idx'
                continue, break
            }
            local ++idx
        }

        * Build label option
        local label_code "label=`this_label'"

        * Build noextend option and extract existing label mappings if needed
        local noextend_code ""
        local existing_code ""
        if "`extend'" == "noextend" {
            local noextend_code "noextend"
            * Check if the label already exists and extract its mappings
            capture label list `this_label'
            if _rc == 0 {
                * Label exists - extract value->string mappings
                local existing_mappings ""

                * Iterate through reasonable range of label values
                forvalues v = 1/1000 {
                    local lbl : label `this_label' `v', strict
                    if `"`lbl'"' != "" {
                        * This value has a label - escape special characters
                        local lbl_escaped = subinstr(`"`lbl'"', "\", "\\", .)
                        local lbl_escaped = subinstr(`"`lbl_escaped'"', "|", "\|", .)
                        if "`existing_mappings'" == "" {
                            local existing_mappings "`v'|`lbl_escaped'"
                        }
                        else {
                            local existing_mappings "`existing_mappings'||`v'|`lbl_escaped'"
                        }
                    }
                }

                if "`existing_mappings'" != "" {
                    local existing_code "existing=`existing_mappings'"
                }
            }
        }

        * Call the C plugin with ALL variables
        unab allvars : *
        plugin call ctools_plugin `allvars' `if' `in', "cencode `threads_code' `var_idx' `gen_idx' `label_code' `noextend_code' `existing_code'"

        * =====================================================================
        * Create value labels from plugin output
        * =====================================================================

        local n_unique = _cencode_n_unique
        local n_chunks = _cencode_n_chunks

        if `n_unique' > 0 {
            * Check if label already exists
            local label_exists = 0
            if "`extend'" == "noextend" {
                capture label list `this_label'
                if _rc == 0 {
                    local label_exists = 1
                }
            }

            if `label_exists' == 0 {
                * Create new label definition from plugin output
                forvalues chunk = 0/`=`n_chunks'-1' {
                    local labels_chunk "${_cencode_labels_`chunk'}"

                    * Parse the label chunk: format is "value|label||value|label||..."
                    while "`labels_chunk'" != "" {
                        * Find next pair delimiter ||
                        local delim_pos = strpos("`labels_chunk'", "||")
                        if `delim_pos' > 0 {
                            local pair = substr("`labels_chunk'", 1, `delim_pos' - 1)
                            local labels_chunk = substr("`labels_chunk'", `delim_pos' + 2, .)
                        }
                        else {
                            local pair = "`labels_chunk'"
                            local labels_chunk = ""
                        }

                        * Parse pair: value|label
                        local pipe_pos = strpos("`pair'", "|")
                        if `pipe_pos' > 0 {
                            local value = substr("`pair'", 1, `pipe_pos' - 1)
                            local labeltext = substr("`pair'", `pipe_pos' + 1, .)

                            * Unescape special characters
                            local labeltext = subinstr("`labeltext'", "\|", "|", .)
                            local labeltext = subinstr("`labeltext'", "\\", "\", .)
                            local labeltext = subinstr(`"`labeltext'"', `"\""', `"""', .)

                            * Add to label definition
                            capture label define `this_label' `value' `"`labeltext'"', add
                            if _rc != 0 & _rc != 110 {
                                * Try modify if label value already exists
                                capture label define `this_label' `value' `"`labeltext'"', modify
                            }
                        }
                    }
                }
            }

            * Apply label to the new variable
            label values `destvar' `this_label'
        }

        * Clean up macros
        forvalues chunk = 0/`=`n_chunks'-1' {
            capture macro drop _cencode_labels_`chunk'
        }
        capture macro drop _cencode_label_name

        * Handle replace option: drop original var and rename temp var
        if `__do_replace' {
            capture drop `srcvar'
            rename `destvar' `srcvar'
        }

        * Store rclass results for this variable
        return scalar N_unique = `n_unique'
    }

    * Store rclass results
    return scalar N_vars = `n_vars'

    if `__do_timing' {
        timer off 90

        * Extract timer values
        quietly timer list 90
        local __time_total = r(t90)

        di as text ""
        di as text "{hline 55}"
        di as text "cencode timing breakdown:"
        di as text "{hline 55}"
        di as text "  Variables encoded:        " as result `n_vars'
        di as text "  Wall clock total:         " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture confirm scalar __cencode_threads_max
        if _rc == 0 {
            local __threads_max = scalar(__cencode_threads_max)
            capture local __openmp_enabled = scalar(__cencode_openmp_enabled)
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"

            * Clean up scalars
            capture scalar drop __cencode_threads_max
            capture scalar drop __cencode_openmp_enabled
        }
    }

    * Clean up any remaining scalars
    capture scalar drop _cencode_n_unique
    capture scalar drop _cencode_n_chunks
end
