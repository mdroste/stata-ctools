*! version 1.4.0 10Jan2026
*! csort: C-accelerated sorting for Stata datasets
*! Part of the ctools suite
*!
*! Supports multiple sort algorithms:
*!   algorithm(lsd)      - LSD radix sort (default) - best for fixed-width keys
*!   algorithm(msd)      - MSD radix sort - best for variable-length strings
*!   algorithm(timsort)  - Timsort - best for partially sorted data
*!   algorithm(sample)   - Sample sort - best for large datasets with many cores
*!   algorithm(counting) - Counting sort - best for integer data with small range
*!   algorithm(merge)    - Parallel merge sort - stable, predictable O(n log n)
*!   algorithm(ips4o)    - IPS4o - in-place parallel super scalar samplesort

program define csort
    version 14.0

    syntax varlist [if] [in], [Verbose TIMEit ALGorithm(string)]

    * Mark sample
    marksample touse, novarlist

    * Count variables for the sort
    local nvars : word count `varlist'

    if (`nvars' == 0) {
        di as error "csort: no sort variables specified"
        exit 198
    }

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

    * Parse algorithm option
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
            di as error "Valid options: lsd (default), msd, timsort, sample, counting, merge, ips4o"
            exit 198
        }
    }

    if ("`verbose'" != "") {
        di as text "csort: Sorting on variables: `varlist'"
        di as text "       All variables: `allvars'"
        di as text "       Sort key indices: `var_indices'"
        if "`alg_code'" != "" {
            di as text "       Algorithm: `algorithm'"
        }
        else {
            di as text "       Algorithm: lsd (default)"
        }
    }

    * Call the C plugin with ALL variables (so it can sort the entire dataset)
    plugin call ctools_plugin `allvars' `if' `in', "csort `var_indices' `alg_code'"

    * Display timing if requested (verbose or timeit)
    if ("`verbose'" != "" | "`timeit'" != "") {
        di as text _n "csort timing breakdown:"
        di as text "  Data load:  " as result %8.4f _csort_time_load " sec"
        di as text "  Sort:       " as result %8.4f _csort_time_sort " sec"
        di as text "  Data store: " as result %8.4f _csort_time_store " sec"
        di as text "  Total:      " as result %8.4f _csort_time_total " sec"
    }

    * Set sort order in Stata (use stable to match csort's stable radix sort)
    sort `varlist', stable
end
