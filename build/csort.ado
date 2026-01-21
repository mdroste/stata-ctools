*! version 1.1.0 19Jan2026
*! csort: C-accelerated sorting for Stata datasets
*! Part of the ctools suite
*!
*! Supports multiple sort algorithms:
*!   algorithm(ips4o)    - IPS4o (default) - in-place parallel super scalar samplesort
*!   algorithm(lsd)      - LSD radix sort - best for fixed-width keys
*!   algorithm(msd)      - MSD radix sort - best for variable-length strings
*!   algorithm(timsort)  - Timsort - best for partially sorted data
*!   algorithm(sample)   - Sample sort - best for large datasets with many cores
*!   algorithm(counting) - Counting sort - best for integer data with small range
*!   algorithm(merge)    - Parallel merge sort - stable, predictable O(n log n)
*!
*! Streaming mode (stream option):
*!   Only loads sort key variables into C memory, then streams the permutation
*!   to non-key variables. Reduces memory usage for wide datasets.
*!   Best when: many non-key columns, limited memory, or dataset too large for RAM.
*!
*! nosortedby option:
*!   Skips setting Stata's internal sort order metadata after sorting.
*!   Use when you don't need Stata to recognize the data as sorted (faster).

program define csort
    version 14.0

    syntax varlist [if] [in], [Verbose ALGorithm(string) STReam THReads(integer 0) NOSORTedby]

    * =========================================================================
    * UPFRONT VALIDATION - check all options before any data manipulation
    * =========================================================================

    * Count variables for the sort
    local nvars : word count `varlist'
    if (`nvars' == 0) {
        di as error "csort: no sort variables specified"
        exit 198
    }

    * Validate algorithm option early (before any data operations)
    local alg_code ""
    if "`algorithm'" != "" {
        local algorithm = lower("`algorithm'")
        if "`algorithm'" == "lsd" | "`algorithm'" == "0" {
            local alg_code "alg=lsd"
        }
        else if "`algorithm'" == "msd" | "`algorithm'" == "1" {
            local alg_code "alg=msd"
        }
        else if "`algorithm'" == "timsort" | "`algorithm'" == "tim" | "`algorithm'" == "2" {
            local alg_code "alg=timsort"
        }
        else if "`algorithm'" == "sample" | "`algorithm'" == "3" {
            local alg_code "alg=sample"
        }
        else if "`algorithm'" == "counting" | "`algorithm'" == "count" | "`algorithm'" == "4" {
            local alg_code "alg=counting"
        }
        else if "`algorithm'" == "merge" | "`algorithm'" == "5" {
            local alg_code "alg=merge"
        }
        else if "`algorithm'" == "ips4o" | "`algorithm'" == "6" {
            local alg_code "alg=ips4o"
        }
        else {
            di as error "csort: invalid algorithm '`algorithm''"
            di as error "Valid options: ips4o (default), lsd, msd, timsort, sample, counting, merge"
            exit 198
        }
    }

    * =========================================================================
    * END UPFRONT VALIDATION
    * =========================================================================

    * Start overall wall-clock timer
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
    marksample touse, novarlist

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
            di as error "csort: Could not load ctools plugin"
            exit 601
        }
    }

    * Get ALL permanent variables in the dataset (we need to sort all of them)
    * Use unab to expand varlist and avoid temporary variables
    unab allvars : *

    * Find indices of sort variables within all variables (1-based)
    local var_indices ""
    foreach sortvar of local varlist {
        local idx = 1
        foreach v of local allvars {
            if ("`v'" == "`sortvar'") {
                local var_indices "`var_indices' `idx'"
                continue, break
            }
            local ++idx
        }
    }

    * Build stream option string
    local stream_code ""
    if "`stream'" != "" {
        local stream_code "stream"
    }

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * End pre-plugin timer, start plugin timer
    if `__do_timing' {
        timer off 91
        timer on 92
    }

    * Call the C plugin with ALL variables (so it can sort the entire dataset)
    plugin call ctools_plugin `allvars' `if' `in', "csort `threads_code' `var_indices' `alg_code' `stream_code'"

    * End plugin timer, start post-plugin timer
    if `__do_timing' {
        timer off 92
        timer on 93
    }

    * Set sort order in Stata (use stable to match csort's stable radix sort)
    * Skip if nosortedby option is specified
    if "`nosortedby'" == "" {
        sort `varlist', stable
    }

    * End all timers
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

        * Calculate plugin call overhead (time in plugin call but outside C code)
        local __plugin_call_overhead = `__time_plugin' - _csort_time_total

        local __ado_overhead = `__time_preplugin' + `__time_postplugin'
        local __ado_total = `__ado_overhead' + `__plugin_call_overhead'

        * Check if stream mode was used
        capture local __is_stream = _csort_stream
        if _rc != 0 local __is_stream = 0

        * Check if permutation time is available
        capture local __time_permute = _csort_time_permute
        if _rc != 0 local __time_permute = 0

        * Check if streaming time is available
        capture local __time_stream = _csort_time_stream
        if _rc != 0 local __time_stream = 0

        di as text ""
        di as text "{hline 55}"
        if `__is_stream' {
            di as text "csort timing breakdown (streaming mode):"
        }
        else {
            di as text "csort timing breakdown:"
        }
        di as text "{hline 55}"
        di as text "  C plugin internals:"

        if `__is_stream' {
            * Memory-efficient mode timing display
            di as text "    Load key vars:          " as result %8.4f _csort_time_load " sec"
            di as text "    Sort (compute order):   " as result %8.4f _csort_time_sort " sec"
            di as text "    Permute keys in C:      " as result %8.4f `__time_permute' " sec"
            di as text "    Store sorted keys:      " as result %8.4f _csort_time_store " sec"
            di as text "    Stream non-key vars:    " as result %8.4f `__time_stream' " sec"
        }
        else {
            * Standard mode timing display
            di as text "    Data load:              " as result %8.4f _csort_time_load " sec"
            if `__time_permute' > 0 {
                di as text "    Sort (compute order):   " as result %8.4f _csort_time_sort " sec"
                di as text "    Sort (apply permute):   " as result %8.4f `__time_permute' " sec"
            }
            else {
                di as text "    Sort:                   " as result %8.4f _csort_time_sort " sec"
            }
            di as text "    Data store:             " as result %8.4f _csort_time_store " sec"
            di as text "    Memory cleanup:         " as result %8.4f _csort_time_cleanup " sec"
        }

        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _csort_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:"
        di as text "    Pre-plugin parsing:     " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin sort:       " as result %8.4f `__time_postplugin' " sec"
        di as text "  {hline 53}"
        di as text "    Stata overhead total:   " as result %8.4f `__ado_total' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = _csort_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _csort_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }
end
