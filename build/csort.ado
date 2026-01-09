*! version 1.1.0 05Jan2026
*! csort: C-accelerated sorting for Stata datasets
*! Part of the ctools suite

program define csort
    version 14.0

    syntax varlist [if] [in], [Verbose TIMEit]

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

    if ("`verbose'" != "") {
        di as text "csort: Sorting on variables: `varlist'"
        di as text "       All variables: `allvars'"
        di as text "       Sort key indices: `var_indices'"
    }

    * Call the C plugin with ALL variables (so it can sort the entire dataset)
    plugin call ctools_plugin `allvars' `if' `in', "csort `var_indices'"

    * Display timing if requested
    if ("`timeit'" != "") {
        di as text _n "Timing breakdown:"
        di as text "  Load:  " as result %8.4f _csort_time_load " sec"
        di as text "  Sort:  " as result %8.4f _csort_time_sort " sec"
        di as text "  Store: " as result %8.4f _csort_time_store " sec"
        di as text "  Total: " as result %8.4f _csort_time_total " sec"
    }

    * Set sort order in Stata (use stable to match csort's stable radix sort)
    sort `varlist', stable
end
