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

    * Get variable indices (1-based for Stata)
    local var_indices ""
    foreach var of local varlist {
        * Find variable index
        local idx = 0
        local allvars : char _dta[__varlist]
        if ("`allvars'" == "") {
            qui ds
            local allvars `r(varlist)'
        }
        local pos = 1
        foreach v of local allvars {
            if ("`v'" == "`var'") {
                local idx = `pos'
                continue, break
            }
            local ++pos
        }
        if (`idx' == 0) {
            * Variable not found in simple list, try describe
            qui describe `var', varlist
            * Just use position
            local idx = `pos'
        }
        local var_indices "`var_indices' `idx'"
    }

    if ("`verbose'" != "") {
        di as text "csort: Sorting on variables: `varlist'"
        di as text "       Variable indices: `var_indices'"
    }

    * Call the C plugin
    plugin call ctools_plugin `varlist' `if' `in', "csort `var_indices'"

    * Display timing if requested
    if ("`timeit'" != "") {
        di as text _n "Timing breakdown:"
        di as text "  Load:  " as result %8.4f _csort_time_load " sec"
        di as text "  Sort:  " as result %8.4f _csort_time_sort " sec"
        di as text "  Store: " as result %8.4f _csort_time_store " sec"
        di as text "  Total: " as result %8.4f _csort_time_total " sec"
    }

    * Set sort order in Stata
    sort `varlist'
end
