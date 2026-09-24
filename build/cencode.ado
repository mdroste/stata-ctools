*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

program define cencode, rclass
    version 14.1
    preserve
    capture noisily _cencode_impl `0'
    local rc = _rc
    if `rc' {
        restore
        exit `rc'
    }
    return add
    restore, not
end

program define _cencode_impl, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    syntax varlist [if] [in], [Generate(string) replace Label(name) NOExtend Verbose THReads(integer 0)]

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
        _ctools_newvars `generate'
    }

    * Check that all source variables are string (match encode rc=107)
    foreach v of local varlist {
        capture confirm string variable `v'
        if _rc != 0 {
            di as error "`v' is not a string variable"
            exit 107
        }
    }

    if _N == 0 {
        forvalues i=1/`n_vars' {
            if `__do_replace' {
                local target : word `i' of `varlist'
                drop `target'
            }
            else local target : word `i' of `generate'
            quietly generate long `target' = .
        }
        return scalar N_unique = 0
        return scalar N_vars = `n_vars'
        exit
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
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0
    * Reset _rc (plugin load may leave _rc=110 for already-loaded plugin)
    capture confirm number 0

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

        capture label list `this_label'
        if !_rc {
            quietly encode `srcvar' `if' `in', generate(`destvar') label(`this_label') `noextend'
            tempvar label_tag
            quietly egen byte `label_tag' = tag(`destvar')
            quietly count if `label_tag'
            return scalar N_unique = r(N)
            if `__do_replace' {
                drop `srcvar'
                rename `destvar' `srcvar'
            }
            * No C phases ran for this variable.
            foreach phase in parse load collect sort encode labels total {
                scalar _cencode_time_`phase' = 0
            }
            continue
        }

        * Create the destination variable (numeric, long type for sufficient range)
        quietly generate long `destvar' = .

        * Get variable indices with a single unab + loop
        * (gen_idx is always the last variable since we just created it)
        unab allvars : *
        local var_idx = 0
        local idx = 1
        foreach v of local allvars {
            if ("`v'" == "`srcvar'") {
                local var_idx = `idx'
            }
            local ++idx
        }
        local gen_idx : word count `allvars'

        if `var_idx' == 0 {
            di as error "cencode: could not find variable `srcvar'"
            exit 111
        }

        * Build label option
        local label_code "label=`this_label'"

        * Build noextend option: save existing label to file for C plugin
        local noextend_code ""
        local existfile_code ""
        if "`noextend'" != "" {
            local noextend_code "noextend"
            * Check if the label already exists
            capture label list `this_label'
            if _rc == 0 {
                * Label exists - save to temp file for C plugin to parse
                tempfile __existlblfile
                quietly label save `this_label' using `__existlblfile', replace
                local existfile_code "existfile=`__existlblfile'"
            }
        }

        * Temp file for C plugin to write label definitions
        tempfile __labelfile
        local labelfile_code "labelfile=`__labelfile'"

        * Set string width metadata for flat buffer optimization
        _ctools_strw `allvars'

        * Call the C plugin with ALL variables
        if `__do_timing' {
            timer off 91
            timer on 92
        }
        plugin call ctools_plugin `allvars' `if' `in', "cencode `threads_code' `var_idx' `gen_idx' `label_code' `noextend_code' `existfile_code' `labelfile_code'"
        if `__do_timing' {
            timer off 92
            timer on 93
        }

        * =====================================================================
        * Create value labels from .do file written by C plugin
        * =====================================================================

        local n_unique = _cencode_n_unique

        if `n_unique' > 0 {
            * Drop any pre-existing label with this name for a clean slate
            capture label drop `this_label'
            * Run the .do file written by C plugin (label define commands)
            run `__labelfile'
        }

        * Apply label to variable if it exists
        capture label list `this_label'
        if _rc == 0 {
            label values `destvar' `this_label'
        }

        * Handle replace option: drop original var and rename temp var
        if `__do_replace' {
            capture drop `srcvar'
            rename `destvar' `srcvar'
        }

        * Store rclass results for this variable
        return scalar N_unique = `n_unique'

        if `__do_timing' {
            timer off 93
            timer on 91
        }
    }

    if `__do_timing' {
        timer off 91
    }

    * Store rclass results
    return scalar N_vars = `n_vars'

    if `__do_timing' {
        timer off 90

        * Extract timer values
        quietly timer list 90
        local __time_total = r(t90)

        * Extract sub-timer values
        quietly timer list 91
        local __time_preplugin = r(t91)
        quietly timer list 92
        local __time_plugin = r(t92)
        quietly timer list 93
        local __time_postplugin = r(t93)

        * Calculate plugin call overhead
        capture local __plugin_time_total = _cencode_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "cencode timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f _cencode_time_parse " sec"
        di as text "    Data load:              " as result %8.4f _cencode_time_load " sec"
        di as text "    Collect unique values:  " as result %8.4f _cencode_time_collect " sec"
        di as text "    Sort:                   " as result %8.4f _cencode_time_sort " sec"
        di as text "    Encode:                 " as result %8.4f _cencode_time_encode " sec"
        di as text "    Apply labels:           " as result %8.4f _cencode_time_labels " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cencode_time_total " sec"
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
        di as text "  Variables encoded:        " as result `n_vars'

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

        * Clean up timing scalars
        capture scalar drop _cencode_time_parse _cencode_time_load _cencode_time_collect
        capture scalar drop _cencode_time_sort _cencode_time_encode _cencode_time_labels _cencode_time_total
        capture scalar drop _cencode_threads_max _cencode_openmp_enabled
    }

    * Clean up any remaining scalars
    capture scalar drop _cencode_n_unique
end
