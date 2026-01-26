*! version 0.9.0 26Jan2026
*! cdestring: C-accelerated string to numeric conversion for Stata
*! Part of the ctools suite
*!
*! Replaces Stata's built-in destring command with a high-performance
*! C implementation featuring:
*!   - Parallel conversion across observations using OpenMP
*!   - Fast numeric parsing optimized for common formats
*!   - Full compatibility with destring syntax and options
*!
*! Syntax: cdestring [varlist] [if] [in], {generate(newvarlist) | replace}
*!         [ignore("chars") force float percent dpcomma Verbose THReads(integer)]
*!
*! The cdestring command inherits the same syntax and functionality as Stata's
*! built-in destring command, producing identical results but with better
*! performance on large datasets.

program define cdestring
    version 14.0

    * Parse syntax - note varlist is optional (defaults to all string vars)
    syntax [varlist(string)] [if] [in], [Generate(string) replace ///
        IGnore(string) force float percent dpcomma Verbose THReads(integer 0)]

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Either generate or replace must be specified, but not both
    if "`generate'" == "" & "`replace'" == "" {
        di as error "cdestring: must specify either generate() or replace option"
        exit 198
    }
    if "`generate'" != "" & "`replace'" != "" {
        di as error "cdestring: cannot specify both generate() and replace options"
        exit 198
    }

    * If no varlist specified, use all string variables
    if "`varlist'" == "" {
        quietly ds, has(type string)
        local varlist "`r(varlist)'"
        if "`varlist'" == "" {
            di as error "cdestring: no string variables found"
            exit 111
        }
    }

    * Count variables
    local nvars : word count `varlist'

    * Validate generate() has same number of variables as varlist
    if "`generate'" != "" {
        local ngen : word count `generate'
        if `ngen' != `nvars' {
            di as error "cdestring: generate() must specify `nvars' variable names"
            exit 198
        }

        * Check that generate variables are valid new names
        foreach v of local generate {
            capture confirm new variable `v'
            if _rc != 0 {
                * Either invalid name or already exists
                capture confirm variable `v'
                if _rc == 0 {
                    di as error "cdestring: variable `v' already exists"
                    exit 110
                }
                else {
                    di as error "cdestring: invalid variable name `v'"
                    exit 198
                }
            }
        }
    }

    * Verify all source variables are string
    foreach v of local varlist {
        capture confirm string variable `v'
        if _rc != 0 {
            di as error "cdestring: `v' is not a string variable"
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
        timer clear 91
        timer clear 92
        timer clear 93
        timer on 90
        timer on 91
    }

    * Mark sample
    marksample touse, strok

    * Load the platform-appropriate ctools plugin
    * Note: program list doesn't work for plugins, so we just try to load
    * and treat _rc=110 (already defined) as success
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
        di as error "cdestring: Could not load ctools plugin"
        exit 601
    }

    * Reset _rc to 0 (110 means plugin already loaded, which is fine)
    capture confirm number 1

    * Create destination variables
    if "`replace'" != "" {
        * For replace, we create temporary variables first, then replace
        local dest_varlist ""
        local tempvars ""
        foreach v of local varlist {
            tempvar tmp_`v'
            if "`float'" != "" {
                quietly generate float `tmp_`v'' = .
            }
            else {
                quietly generate double `tmp_`v'' = .
            }
            local dest_varlist "`dest_varlist' `tmp_`v''"
            local tempvars "`tempvars' `tmp_`v''"
        }
    }
    else {
        * For generate, create the new variables
        local dest_varlist ""
        local i = 1
        foreach gv of local generate {
            if "`float'" != "" {
                quietly generate float `gv' = .
            }
            else {
                quietly generate double `gv' = .
            }
            local dest_varlist "`dest_varlist' `gv'"
            local ++i
        }
    }

    * Get all variable names and build indices for plugin call
    * Pass ALL variables to plugin so indices match dataset positions
    unab allvars : *

    * Build list of source variable indices (position in allvars)
    local src_indices ""
    foreach v of local varlist {
        local idx = 1
        foreach av of local allvars {
            if "`av'" == "`v'" {
                local src_indices "`src_indices' `idx'"
                continue, break
            }
            local ++idx
        }
    }

    * Build list of destination variable indices (position in allvars)
    local dst_indices ""
    foreach v of local dest_varlist {
        local idx = 1
        foreach av of local allvars {
            if "`av'" == "`v'" {
                local dst_indices "`dst_indices' `idx'"
                continue, break
            }
            local ++idx
        }
    }


    * Build option strings for C plugin
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    local ignore_code ""
    if `"`ignore'"' != "" {
        * Escape special characters for passing to C
        * Order matters: first escape backslashes, then encode spaces
        local ignore_escaped = subinstr(`"`ignore'"', "\", "\\", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', " ", "\s", .)
        local ignore_code "ignore=`ignore_escaped'"
    }

    local force_code ""
    if "`force'" != "" {
        local force_code "force"
    }

    local percent_code ""
    if "`percent'" != "" {
        local percent_code "percent"
    }

    local dpcomma_code ""
    if "`dpcomma'" != "" {
        local dpcomma_code "dpcomma"
    }

    local verbose_code ""
    if "`verbose'" != "" {
        local verbose_code "verbose"
    }

    if `__do_timing' {
        timer off 91
        timer on 92
    }

    * Build the full command string
    * Format: src1 dst1 src2 dst2 ... nvars=N [options]
    local cmd_indices ""
    local i = 1
    foreach src_idx of local src_indices {
        local dst_idx : word `i' of `dst_indices'
        local cmd_indices "`cmd_indices' `src_idx' `dst_idx'"
        local ++i
    }

    * Call the C plugin with ALL variables (so indices match dataset positions)
    plugin call ctools_plugin `allvars' `if' `in', ///
        "cdestring `threads_code' `cmd_indices' nvars=`nvars' `ignore_code' `force_code' `percent_code' `dpcomma_code' `verbose_code'"

    local plugin_rc = _rc

    if `__do_timing' {
        timer off 92
        timer on 93
    }

    * Get results from plugin
    local n_converted = _cdestring_n_converted
    local n_failed = _cdestring_n_failed

    * Handle replace option - copy temp vars to originals and drop temps
    if "`replace'" != "" & `plugin_rc' == 0 {
        local i = 1
        foreach v of local varlist {
            local tmpv : word `i' of `tempvars'
            * Drop the original string variable
            drop `v'
            * Rename temp to original name
            rename `tmpv' `v'
            local ++i
        }
    }
    else if "`replace'" != "" {
        * Clean up temp vars on error
        foreach tmpv of local tempvars {
            capture drop `tmpv'
        }
    }

    * Compress the converted variables to optimal storage type (mimics destring)
    if `plugin_rc' == 0 {
        if "`replace'" != "" {
            quietly compress `varlist'
        }
        else {
            quietly compress `generate'
        }
    }

    * Display output similar to destring
    if `plugin_rc' == 0 {
        if "`replace'" != "" {
            foreach v of local varlist {
                local vtype : type `v'
                di as text "`v': " as text "string to " as result "`vtype'"
            }
        }
        else {
            local i = 1
            foreach gv of local generate {
                local sv : word `i' of `varlist'
                local vtype : type `gv'
                di as text "`sv' -> `gv': " as text "string to " as result "`vtype'"
                local ++i
            }
        }
    }

    * Report any failures
    if `n_failed' > 0 & "`force'" == "" {
        di as text "(`n_failed' observation(s) contained nonnumeric characters; converted to missing)"
    }

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
        capture local __plugin_time_total = _cdestring_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "cdestring timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f _cdestring_time_parse " sec"
        di as text "    String conversion:      " as result %8.4f _cdestring_time_convert " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cdestring_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Observations processed:   " as result `n_converted'
        if `n_failed' > 0 {
            di as text "  Non-numeric (-> missing): " as result `n_failed'
        }

        * Display thread diagnostics
        capture local __threads_max = _cdestring_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cdestring_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    if `plugin_rc' != 0 {
        exit `plugin_rc'
    }

    * Ensure clean exit
    exit 0
end
