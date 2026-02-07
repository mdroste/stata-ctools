*! version 0.9.1 06Feb2026
*! csort: C-accelerated sorting for Stata datasets
*! Part of the ctools suite
*!
*! By default, csort auto-selects the optimal algorithm based on data type:
*!   - String variables: MSD radix sort
*!   - Numeric variables: Counting sort (falls back to LSD radix for large ranges)
*!
*! Explicit algorithm options:
*!   algorithm(auto)     - Auto-select (default)
*!   algorithm(counting) - Counting sort - optimal for small-range integers
*!   algorithm(lsd)      - LSD radix sort - for fixed-width numeric keys
*!   algorithm(msd)      - MSD radix sort - for string variables
*!   algorithm(ips4o)    - IPS4o - parallel samplesort
*!   algorithm(timsort)  - Timsort - for partially sorted data
*!   algorithm(sample)   - Sample sort - for very large datasets
*!   algorithm(merge)    - Parallel merge sort - stable O(n log n)
*!
*! Streaming mode (stream option):
*!   Only loads sort key variables into C memory, then streams the permutation
*!   to non-key variables. Reduces memory usage for wide datasets.
*!   Best when: many non-key columns, limited memory, or dataset too large for RAM.
*!
*! nosortedby option:
*!   Skips setting Stata's internal sort order metadata after sorting.
*!   Use when you don't need Stata to recognize the data as sorted (faster).

program define csort, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    syntax varlist [if] [in], [Verbose ALGorithm(string) STReam(integer 0) THReads(integer 0) NOSORTedby NOSTReam]

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
    * Default is "auto" which intelligently selects based on data characteristics
    local alg_code ""
    local alg_name "auto"
    if "`algorithm'" != "" {
        local algorithm = lower("`algorithm'")
        if "`algorithm'" == "lsd" | "`algorithm'" == "0" {
            local alg_code "alg=lsd"
            local alg_name "lsd"
        }
        else if "`algorithm'" == "msd" | "`algorithm'" == "1" {
            local alg_code "alg=msd"
            local alg_name "msd"
        }
        else if "`algorithm'" == "timsort" | "`algorithm'" == "tim" | "`algorithm'" == "2" {
            local alg_code "alg=timsort"
            local alg_name "timsort"
        }
        else if "`algorithm'" == "sample" | "`algorithm'" == "3" {
            local alg_code "alg=sample"
            local alg_name "sample"
        }
        else if "`algorithm'" == "counting" | "`algorithm'" == "count" | "`algorithm'" == "4" {
            local alg_code "alg=counting"
            local alg_name "counting"
        }
        else if "`algorithm'" == "merge" | "`algorithm'" == "5" {
            local alg_code "alg=merge"
            local alg_name "merge"
        }
        else if "`algorithm'" == "ips4o" | "`algorithm'" == "6" {
            local alg_code "alg=ips4o"
            local alg_name "ips4o"
        }
        else if "`algorithm'" == "auto" | "`algorithm'" == "7" {
            local alg_code "alg=auto"
            local alg_name "auto"
        }
        else {
            di as error "csort: invalid algorithm '`algorithm''"
            di as error "Valid options: auto (default), lsd, msd, timsort, sample, counting, merge, ips4o"
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
    * stream(#) enables streaming mode with # variables loaded at a time (1-16)
    local stream_code ""
    if `stream' > 0 {
        * Validate and constrain batch size
        local batch_size = `stream'
        if `batch_size' > 16 {
            local batch_size = 16
        }
        local stream_code "stream(`batch_size')"
    }

    * Build nostream option string (disables auto-streaming)
    local nostream_code ""
    if "`nostream'" != "" {
        local nostream_code "nostream"
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
    plugin call ctools_plugin `allvars' `if' `in', "csort `threads_code' `var_indices' `alg_code' `stream_code' `nostream_code'"

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

            * Detailed breakdown of stream non-key vars phase
            capture local __stream_build_inv = _csort_stream_time_build_inv
            if _rc == 0 {
                capture local __stream_scatter = _csort_stream_time_scatter
                capture local __stream_writeback = _csort_stream_time_writeback
                capture local __stream_strings = _csort_stream_time_strings
                capture local __n_numeric = _csort_stream_n_numeric
                capture local __n_string = _csort_stream_n_string
                if _rc == 0 {
                    di as text "      ├─ Build inv perm:    " as result %8.4f `__stream_build_inv' " sec"
                    di as text "      ├─ Scatter (read+buf):" as result %8.4f `__stream_scatter' " sec" as text " (" as result `__n_numeric' as text " numeric vars)"
                    di as text "      ├─ Writeback to Stata:" as result %8.4f `__stream_writeback' " sec"
                    if `__n_string' > 0 {
                        di as text "      └─ String vars:       " as result %8.4f `__stream_strings' " sec" as text " (" as result `__n_string' as text " string vars)"
                    }
                }
            }
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

    * =========================================================================
    * STORE RESULTS IN r()
    * =========================================================================

    * Get values from scalars set by C plugin
    capture local __time_load = _csort_time_load
    if _rc != 0 local __time_load = .
    capture local __time_sort = _csort_time_sort
    if _rc != 0 local __time_sort = .
    capture local __time_permute = _csort_time_permute
    if _rc != 0 local __time_permute = .
    capture local __time_store = _csort_time_store
    if _rc != 0 local __time_store = .
    capture local __time_stream = _csort_time_stream
    if _rc != 0 local __time_stream = .
    capture local __time_cleanup = _csort_time_cleanup
    if _rc != 0 local __time_cleanup = .
    capture local __time_total_c = _csort_time_total
    if _rc != 0 local __time_total_c = .
    capture local __is_stream = _csort_stream
    if _rc != 0 local __is_stream = 0
    capture local __threads_max = _csort_threads_max
    if _rc != 0 local __threads_max = .
    capture local __openmp_enabled = _csort_openmp_enabled
    if _rc != 0 local __openmp_enabled = .

    * Count observations in sample
    quietly count `if' `in'
    local __nobs = r(N)

    * Post r() results
    return scalar N = `__nobs'
    return scalar time_load = `__time_load'
    return scalar time_sort = `__time_sort'
    return scalar time_permute = `__time_permute'
    return scalar time_store = `__time_store'
    return scalar time_stream = `__time_stream'
    return scalar time_cleanup = `__time_cleanup'
    return scalar time_total = `__time_total_c'
    return scalar stream = `__is_stream'
    return scalar threads_max = `__threads_max'
    return scalar openmp_enabled = `__openmp_enabled'

    return local algorithm "`alg_name'"
    return local sortvars "`varlist'"
    return local cmd "csort"

end
