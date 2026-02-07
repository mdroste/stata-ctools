*! version 0.9.1 06Feb2026
*! crangestat: C-accelerated range statistics for Stata
*! Part of the ctools suite

program define crangestat
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Parse the syntax - similar to rangestat
    * crangestat (stat) [newvar=]varname [(stat) [newvar=]varname ...], interval(keyvar low high) [by(varlist) excludeself verbose threads(#)]

    * First, separate the stat specifications from the options
    local 0_orig "`0'"

    * Find the comma to split stats from options
    local comma_pos = strpos("`0'", ",")
    if `comma_pos' == 0 {
        di as error "crangestat: interval() option is required"
        exit 198
    }

    local stat_part = substr("`0'", 1, `comma_pos' - 1)
    local opt_part = substr("`0'", `comma_pos' + 1, .)

    * Parse options
    local 0 ", `opt_part'"
    syntax, INTerval(string) [BY(varlist) EXCLUDEself Verbose THReads(integer 0)]

    * Parse interval specification: keyvar low high
    tokenize "`interval'"
    local keyvar "`1'"
    local int_low "`2'"
    local int_high "`3'"

    if "`keyvar'" == "" | "`int_low'" == "" | "`int_high'" == "" {
        di as error "crangestat: interval() requires keyvar low high"
        exit 198
    }

    * Validate key variable exists and is numeric
    capture confirm variable `keyvar'
    if _rc != 0 {
        di as error "crangestat: variable `keyvar' not found"
        exit 111
    }
    capture confirm numeric variable `keyvar'
    if _rc != 0 {
        di as error "crangestat: key variable `keyvar' must be numeric"
        exit 109
    }

    * Convert interval bounds (handle . for +/- infinity)
    local low_is_missing = 0
    local high_is_missing = 0

    if "`int_low'" == "." {
        local low_is_missing = 1
        local int_low_val = 0
    }
    else {
        capture confirm number `int_low'
        if _rc != 0 {
            di as error "crangestat: invalid interval low bound '`int_low''"
            exit 198
        }
        local int_low_val = `int_low'
    }

    if "`int_high'" == "." {
        local high_is_missing = 1
        local int_high_val = 0
    }
    else {
        capture confirm number `int_high'
        if _rc != 0 {
            di as error "crangestat: invalid interval high bound '`int_high''"
            exit 198
        }
        local int_high_val = `int_high'
    }

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        timer clear 90
        timer clear 91
        timer clear 92
        timer on 90
        timer on 91
    }

    * Parse stat specifications
    * Format: (stat) [newvar=]varname [(stat) [newvar=]varname ...]
    local stat_codes ""
    local source_vars ""
    local result_vars ""
    local source_var_map ""
    local nstats = 0
    local nsource = 0
    local source_var_list ""

    local remaining "`stat_part'"

    while "`remaining'" != "" {
        * Skip leading whitespace
        local remaining = strtrim("`remaining'")
        if "`remaining'" == "" continue, break

        * Expect opening parenthesis
        if substr("`remaining'", 1, 1) != "(" {
            di as error "crangestat: expected '(' before statistic name"
            exit 198
        }

        * Find closing parenthesis
        local paren_end = strpos("`remaining'", ")")
        if `paren_end' == 0 {
            di as error "crangestat: missing ')' after statistic name"
            exit 198
        }

        * Extract stat name
        local stat_name = strtrim(substr("`remaining'", 2, `paren_end' - 2))
        local remaining = strtrim(substr("`remaining'", `paren_end' + 1, .))

        * Validate stat name and get code
        local stat_code = -1
        if "`stat_name'" == "count" local stat_code = 0
        else if "`stat_name'" == "mean" local stat_code = 1
        else if "`stat_name'" == "sum" local stat_code = 2
        else if "`stat_name'" == "min" local stat_code = 3
        else if "`stat_name'" == "max" local stat_code = 4
        else if "`stat_name'" == "sd" local stat_code = 5
        else if "`stat_name'" == "variance" local stat_code = 6
        else if "`stat_name'" == "median" local stat_code = 7
        else if "`stat_name'" == "iqr" local stat_code = 8
        else if "`stat_name'" == "first" local stat_code = 9
        else if "`stat_name'" == "last" local stat_code = 10
        else if "`stat_name'" == "firstnm" local stat_code = 11
        else if "`stat_name'" == "lastnm" local stat_code = 12
        else if "`stat_name'" == "p1" local stat_code = 13
        else if "`stat_name'" == "p5" local stat_code = 14
        else if "`stat_name'" == "p10" local stat_code = 15
        else if "`stat_name'" == "p25" local stat_code = 16
        else if "`stat_name'" == "p75" local stat_code = 17
        else if "`stat_name'" == "p90" local stat_code = 18
        else if "`stat_name'" == "p95" local stat_code = 19
        else if "`stat_name'" == "p99" local stat_code = 20
        else if "`stat_name'" == "skewness" local stat_code = 21
        else if "`stat_name'" == "kurtosis" local stat_code = 22

        if `stat_code' == -1 {
            di as error "crangestat: unknown statistic '`stat_name''"
            exit 3499
        }

        * Parse variable specification: [newvar=]varname
        * Find next ( or end
        local next_paren = strpos("`remaining'", "(")
        if `next_paren' == 0 {
            local var_spec = strtrim("`remaining'")
            local remaining ""
        }
        else {
            local var_spec = strtrim(substr("`remaining'", 1, `next_paren' - 1))
            local remaining = substr("`remaining'", `next_paren', .)
        }

        if "`var_spec'" == "" {
            di as error "crangestat: missing variable after (`stat_name')"
            exit 198
        }

        * Check for newvar=varname format
        local eq_pos = strpos("`var_spec'", "=")
        if `eq_pos' > 0 {
            local result_var = strtrim(substr("`var_spec'", 1, `eq_pos' - 1))
            local source_var = strtrim(substr("`var_spec'", `eq_pos' + 1, .))
        }
        else {
            * Default: result var is stat_name_varname
            local source_var "`var_spec'"
            local result_var "`stat_name'_`source_var'"
        }

        * Validate source variable
        capture confirm variable `source_var'
        if _rc != 0 {
            di as error "crangestat: variable `source_var' not found"
            exit 111
        }
        capture confirm numeric variable `source_var'
        if _rc != 0 {
            di as error "crangestat: source variable `source_var' must be numeric"
            exit 109
        }

        * Create result variable if it doesn't exist
        capture confirm variable `result_var'
        if _rc != 0 {
            quietly gen double `result_var' = .
        }
        else {
            * Variable exists - check it's numeric and replace with missing
            capture confirm numeric variable `result_var'
            if _rc != 0 {
                di as error "crangestat: result variable `result_var' must be numeric"
                exit 109
            }
            quietly replace `result_var' = .
        }

        * Track source variable index
        local src_idx = 0
        local found = 0
        foreach sv of local source_var_list {
            local src_idx = `src_idx' + 1
            if "`sv'" == "`source_var'" {
                local found = 1
                continue, break
            }
        }
        if `found' == 0 {
            local source_var_list "`source_var_list' `source_var'"
            local nsource = `nsource' + 1
            local src_idx = `nsource'
        }

        * Store mapping (0-based for C)
        local source_var_map "`source_var_map' `=`src_idx' - 1'"

        local stat_codes "`stat_codes' `stat_code'"
        local result_vars "`result_vars' `result_var'"
        local nstats = `nstats' + 1
    }

    if `nstats' == 0 {
        di as error "crangestat: no statistics specified"
        exit 198
    }

    * Load plugin
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
            di as error "crangestat: Could not load ctools plugin"
            exit 601
        }
    }

    * Get variable indices
    * Build list of ALL dataset variables - this determines plugin varlist order
    * The plugin interface uses 1-based indices relative to the varlist passed to it
    unab all_vars : *

    * Get the position of each variable in the FULL dataset varlist
    local pos = 1
    foreach v of local all_vars {
        local _varpos_`v' = `pos'
        local ++pos
    }

    * Key variable index in dataset order
    local key_idx = `_varpos_`keyvar''

    * Source variable indices
    local source_indices ""
    foreach sv of local source_var_list {
        local source_indices "`source_indices' `_varpos_`sv''"
    }

    * Result variable indices
    local result_indices ""
    foreach rv of local result_vars {
        local result_indices "`result_indices' `_varpos_`rv''"
    }

    * By-variable indices
    local nby = 0
    local by_indices ""
    if "`by'" != "" {
        local nby : word count `by'
        foreach bv of local by {
            local by_indices "`by_indices' `_varpos_`bv''"
        }
    }

    * Build options string
    local opts ""
    if "`excludeself'" != "" local opts "`opts' excludeself"
    if "`verbose'" != "" local opts "`opts' verbose"
    local opts "`opts' low=`int_low_val' high=`int_high_val'"
    if `low_is_missing' local opts "`opts' low_miss"
    if `high_is_missing' local opts "`opts' high_miss"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Build argument string for plugin
    * Format: nstats nsource nby key_idx source_indices... result_indices... by_indices... stat_types... source_var_map... [options]
    local plugin_args "`nstats' `nsource' `nby' `key_idx' `source_indices' `result_indices' `by_indices' `stat_codes' `source_var_map' `opts'"

    * Call plugin
    if `__do_timing' {
        timer off 91
        timer on 92
    }
    plugin call ctools_plugin `all_vars', "crangestat `threads_opt' `plugin_args'"
    if `__do_timing' {
        timer off 92
    }

    * Display summary or timing
    if `__do_timing' {
        timer off 90
        quietly timer list 90
        local __time_total = r(t90)

        * Extract sub-timer values
        quietly timer list 91
        local __time_preplugin = r(t91)
        quietly timer list 92
        local __time_plugin = r(t92)

        * Calculate plugin call overhead
        capture local __plugin_time_total = _crangestat_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "crangestat timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f _crangestat_time_load " sec"
        di as text "    Sort:                   " as result %8.4f _crangestat_time_sort " sec"
        di as text "    Group detection:        " as result %8.4f _crangestat_time_groups " sec"
        di as text "    Compute statistics:     " as result %8.4f _crangestat_time_compute " sec"
        di as text "    Store results:          " as result %8.4f _crangestat_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _crangestat_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__plugin_call_overhead'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = _crangestat_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _crangestat_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clean up timing scalars
        capture scalar drop _crangestat_time_load _crangestat_time_sort _crangestat_time_groups
        capture scalar drop _crangestat_time_compute _crangestat_time_store _crangestat_time_total
        capture scalar drop _crangestat_threads_max _crangestat_openmp_enabled
    }
    else {
        di as text "(`nstats' range statistics computed)"
    }
end
