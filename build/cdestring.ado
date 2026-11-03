*! version 1.2.3 20261030 github.com/mdroste/stata-ctools

program define cdestring
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Parse syntax - note varlist is optional (defaults to all string vars)
    syntax [varlist] [if] [in], [Generate(string) replace ///
        IGnore(string) force float percent dpcomma Verbose THReads(integer 0)]
    * The if is evaluated now, before any helper replaces the r() or e()
    * results it may name.
    local iftouse ""
    if `"`if'"' != "" {
        tempvar iftouse
        mark `iftouse' `if' `in'
        local if "if `iftouse'"
    }
    if !inrange(`threads', 0, 256) exit 198

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

        * Native destring rejects duplicate destinations before creating any.
        local unique : list uniq generate
        if `: word count `unique'' != `ngen' {
            di as error "generate() contains duplicate variable names"
            exit 110
        }
        _ctools_newvars `generate'
    }

    * destring converts a repeated variable once and then finds it numeric.
    * Staging it twice would make the group rename fail.
    if "`replace'" != "" local varlist : list uniq varlist

    * Handle already-numeric variables (match destring behavior: skip gracefully)
    local string_vars ""
    local string_gen ""
    local _var_idx = 0
    foreach v of local varlist {
        local ++_var_idx
        capture confirm string variable `v'
        if _rc != 0 {
            * Variable is already numeric - match destring behavior
            if "`replace'" != "" {
                local vtype : type `v'
                di as text "`v': already " as result "`vtype'"
            }
            else if "`generate'" != "" {
                local vtype : type `v'
                di as text "`v': already " as result "`vtype'"
            }
        }
        else {
            local string_vars "`string_vars' `v'"
            if "`generate'" != "" {
                local gv : word `_var_idx' of `generate'
                local string_gen "`string_gen' `gv'"
            }
        }
    }

    * Update varlist and generate to only include string variables
    local varlist "`string_vars'"
    if "`generate'" != "" {
        local generate "`string_gen'"
    }
    local nvars : word count `varlist'

    * If no string variables remain, exit successfully (like destring)
    if `nvars' == 0 {
        exit 0
    }

    * =========================================================================
    * END VALIDATION
    * =========================================================================

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        local __ct_elapsed93 = 0
        if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
        quietly _ctools_clock __ct_start90
        if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
        quietly _ctools_clock __ct_start91
    }

    * Load the platform-appropriate ctools plugin (cached after first call)
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Both modes stage values and metadata. A plugin error or Break then
    * removes only temporary variables, without publishing partial outputs.
    unab original_order : _all
    local dest_varlist ""
    local i = 0
    foreach v of local varlist {
        local ++i
        tempvar tmp`i'
        if "`float'" != "" {
            quietly generate float `tmp`i'' = .
        }
        else {
            quietly generate double `tmp`i'' = .
        }
        local dest_varlist "`dest_varlist' `tmp`i''"
    }

    * Get all variable names and build indices for plugin call
    * Pass ALL variables to plugin so indices match dataset positions
    unab allvars : *
    local allvars : list allvars - iftouse

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
        * Escape backslashes first so literal \t differs from an actual tab.
        local ignore_escaped = subinstr(`"`ignore'"', "\", "\\", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', " ", "\s", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', char(9), "\t", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', char(10), "\n", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', char(13), "\r", .)
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

    * The plugin needs float for the float range and percent rounding
    local float_code ""
    if "`float'" != "" {
        local float_code "float"
    }

    local verbose_code ""
    if "`verbose'" != "" {
        local verbose_code "verbose"
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
        if "`__ct_elapsed92'" == "" local __ct_elapsed92 = 0
        quietly _ctools_clock __ct_start92
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

    * Set string width metadata for flat buffer optimization
    _ctools_strw `allvars'

    * Stata answers an if one observation at a time. Without one the plugin
    * skips the check; the if marked after syntax is passed last, as a
    * sample marker the plugin reads in parallel.
    local select_code "skipif"
    local marker ""
    local ifcond ""
    if `"`if'"' != "" {
        local marker "`iftouse'"
        local select_code "touse"
    }

    * Call the C plugin with ALL variables (so indices match dataset positions).
    * ignore= goes last so its characters cannot be read as option keywords.
    * The plugin sets local cdestring_failed (nonnumeric count per variable).
    local cdestring_failed ""
    capture noisily plugin call ctools_plugin `allvars' `marker' `ifcond' `in', ///
        "cdestring `threads_code' `cmd_indices' nvars=`nvars' `force_code' `percent_code' `dpcomma_code' `float_code' `verbose_code' `select_code' `ignore_code'"

    local plugin_rc = _rc
    if `plugin_rc' exit `plugin_rc'

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')
        if "`__ct_elapsed93'" == "" local __ct_elapsed93 = 0
        quietly _ctools_clock __ct_start93
    }

    * Get results from plugin
    local n_converted = _cdestring_n_converted
    local n_failed = _cdestring_n_failed

    * Like destring, a variable with any nonnumeric value (after removing
    * ignore() characters) is left untouched unless force is specified: no
    * generate, no replace (the plugin did not store it). Others proceed.
    local genrep = cond("`replace'" != "", "replace", "generate")
    local publish
    local final_names
    local backups
    local i = 1
    foreach v of local varlist {
        local dv : word `i' of `dest_varlist'
        local nbad : word `i' of `cdestring_failed'
        if "`nbad'" == "" local nbad 0
        if `nbad' > 0 & "`force'" == "" {
            drop `dv'
            if `"`ignore'"' != "" | "`percent'" != "" {
                di as text "`v': contains characters not specified in " ///
                    as result "ignore()" as text "; no " as result "`genrep'"
            }
            else {
                di as text "`v': contains nonnumeric characters; no " as result "`genrep'"
            }
        }
        else {
            local vlabel : variable label `v'
            label variable `dv' `"`vlabel'"'
            local characteristics : char `v'[]
            foreach ch of local characteristics {
                local value : char `v'[`ch']
                char `dv'[`ch'] `"`value'"'
            }
            * Compress to the optimal storage type (mimics destring)
            quietly compress `dv'
            local vtype : type `dv'
            local publish `publish' `dv'
            if "`replace'" != "" {
                tempvar backup`i'
                local backups `backups' `backup`i''
                local final_names `final_names' `v'
                di as text "`v': " as text "string to " as result "`vtype'"
            }
            else {
                local final : word `i' of `generate'
                local final_names `final_names' `final'
                di as text "`v' -> `final': " as text "string to " as result "`vtype'"
            }
        }
        local ++i
    }

    * Keep originals until the complete name swap succeeds. Group rename
    * validates its whole mapping; nobreak prevents interruption mid-commit.
    if "`publish'" != "" {
        nobreak {
            if "`replace'" != "" rename (`final_names') (`backups')
            capture rename (`publish') (`final_names')
            local publish_rc = _rc
            if `publish_rc' {
                if "`replace'" != "" rename (`backups') (`final_names')
                exit `publish_rc'
            }
            if "`replace'" != "" {
                drop `backups'
                order `original_order'
            }
        }
    }

    * Report values that force converted to missing
    if `n_failed' > 0 & "`force'" != "" {
        di as text "(`n_failed' observation(s) contained nonnumeric characters; converted to missing)"
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed93 = `__ct_elapsed93' + (`__ct_now' - `__ct_start93')
        quietly _ctools_clock __ct_now
        local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

        * Extract timer values
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
        * Elapsed seconds are held in caller-local timing state.
        local __time_preplugin = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_plugin = `__ct_elapsed92'
        * Elapsed seconds are held in caller-local timing state.
        local __time_postplugin = `__ct_elapsed93'

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
        di as text "    Data load:              " as result %8.4f _cdestring_time_load " sec"
        di as text "    String conversion:      " as result %8.4f _cdestring_time_convert " sec"
        di as text "    Store results:          " as result %8.4f _cdestring_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cdestring_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__plugin_call_overhead' + `__time_postplugin'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Numeric or missing:       " as result `n_converted'
        if `n_failed' > 0 {
            di as text "  Nonnumeric:               " as result `n_failed'
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

    * Ensure clean exit
    exit 0
end
