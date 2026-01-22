*! version 1.0.0 21Jan2026
*! cencode: C-accelerated string encoding for Stata
*! Part of the ctools suite
*!
*! Replaces Stata's built-in encode command with a high-performance
*! C implementation featuring:
*!   - Parallel unique string collection
*!   - Lock-free parallel encoding
*!   - Efficient string interning with arena allocator
*!
*! Syntax: cencode varname [if] [in], generate(newvar) [label(name) noextend Verbose THReads(integer)]
*!
*! The cencode command inherits the same syntax and functionality as Stata's
*! built-in encode command, producing identical results but with better
*! performance on large datasets.

program define cencode
    version 14.0

    syntax varname(string) [if] [in], Generate(name) [Label(name) noextend Verbose THReads(integer 0)]

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Check that source variable is string
    capture confirm string variable `varlist'
    if _rc != 0 {
        di as error "cencode: `varlist' is not a string variable"
        exit 198
    }

    * Check that generate variable doesn't already exist
    capture confirm variable `generate'
    if _rc == 0 {
        di as error "cencode: variable `generate' already exists"
        exit 110
    }

    * Set default label name if not specified
    if "`label'" == "" {
        local label "`generate'"
    }

    * =========================================================================
    * END VALIDATION
    * =========================================================================

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        timer clear 90
        timer clear 91
        timer clear 92
        timer clear 93
        timer on 90
        timer on 91
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

    * Get variable index
    unab allvars : *
    local var_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`varlist'") {
            local var_idx = `idx'
            continue, break
        }
        local ++idx
    }

    if `var_idx' == 0 {
        di as error "cencode: could not find variable `varlist'"
        exit 111
    }

    * Create the destination variable (numeric, long type for sufficient range)
    quietly generate long `generate' = .

    * Get index of new variable
    unab allvars : *
    local gen_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`generate'") {
            local gen_idx = `idx'
            continue, break
        }
        local ++idx
    }

    * Build threads option
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Build label option
    local label_code "label=`label'"

    * Build noextend option
    local noextend_code ""
    if "`noextend'" != "" {
        local noextend_code "noextend"
    }

    if `__do_timing' {
        timer off 91
        timer on 92
    }

    * Call the C plugin
    plugin call ctools_plugin `varlist' `generate' `if' `in', "cencode `threads_code' `var_idx' `gen_idx' `label_code' `noextend_code'"

    if `__do_timing' {
        timer off 92
        timer on 93
    }

    * =========================================================================
    * Create value labels from plugin output
    * =========================================================================

    local n_unique = _cencode_n_unique
    local n_chunks = _cencode_n_chunks

    if `n_unique' > 0 {
        * Check if label already exists
        local label_exists = 0
        if "`noextend'" != "" {
            capture label list `label'
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
                        capture label define `label' `value' `"`labeltext'"', add
                        if _rc != 0 & _rc != 110 {
                            * Try modify if label value already exists
                            capture label define `label' `value' `"`labeltext'"', modify
                        }
                    }
                }
            }
        }

        * Apply label to the new variable
        label values `generate' `label'
    }

    * Clean up macros
    forvalues chunk = 0/`=`n_chunks'-1' {
        capture macro drop _cencode_labels_`chunk'
    }
    capture macro drop _cencode_label_name

    if `__do_timing' {
        timer off 93
        timer off 90

        * Extract timer values
        quietly timer list 90
        local __time_total = r(t90)
        quietly timer list 91
        local __time_preplugin = r(t91)
        quietly timer list 92
        local __time_plugin = r(t92)
        quietly timer list 93
        local __time_postplugin = r(t93)

        * Calculate plugin call overhead
        local __plugin_call_overhead = `__time_plugin' - _cencode_time_total

        di as text ""
        di as text "{hline 55}"
        di as text "cencode timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f _cencode_time_parse " sec"
        di as text "    Load strings:           " as result %8.4f _cencode_time_load " sec"
        di as text "    Collect unique:         " as result %8.4f _cencode_time_collect " sec"
        di as text "    Sort unique:            " as result %8.4f _cencode_time_sort " sec"
        di as text "    Encode values:          " as result %8.4f _cencode_time_encode " sec"
        di as text "    Store labels:           " as result %8.4f _cencode_time_labels " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cencode_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin labels:     " as result %8.4f `__time_postplugin' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Unique values encoded:    " as result `n_unique'

        * Display thread diagnostics
        capture local __threads_max = _cencode_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cencode_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }
end
