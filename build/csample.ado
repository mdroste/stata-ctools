*! version 1.0.1 07Feb2026
*! csample: C-accelerated random sampling for Stata
*! Part of the ctools suite

program define csample
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    timer clear 90
    timer on 90

    syntax [anything(name=percent)] [if] [in], [Count(integer 0) BY(varlist) Verbose THReads(integer 0)]

    local __do_timing = ("`verbose'" != "")

    * Validate percent or count
    local use_percent = 0
    local use_count = 0
    local pct_value = 0

    if "`percent'" != "" {
        capture confirm number `percent'
        if _rc != 0 {
            di as error "csample: percent must be a number"
            exit 198
        }
        if `percent' < 0 | `percent' > 100 {
            di as error "csample: percent must be between 0 and 100"
            exit 198
        }
        local use_percent = 1
        local pct_value = `percent'
    }

    if `count' > 0 {
        local use_count = 1
    }

    if `use_percent' == 0 & `use_count' == 0 {
        di as error "csample: must specify either percent or count()"
        exit 198
    }

    if `use_percent' == 1 & `use_count' == 1 {
        di as error "csample: cannot specify both percent and count()"
        exit 198
    }

    * Mark sample
    marksample touse, novarlist

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
    if `nobs_sample' == 0 {
        exit 0
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
            di as error "csample: Could not load ctools plugin"
            exit 601
        }
    }

    * Generate seed from Stata's RNG (respects -set seed-)
    * Use runiform to generate a large integer seed
    local stata_seed = floor(runiform() * 2147483647) + floor(runiform() * 2147483647) * 2147483648

    * Count by-variables
    local nby = 0
    if "`by'" != "" {
        local nby : word count `by'
    }

    * Create variable for keep flags (will be dropped at end)
    capture drop __csample_keep__
    quietly gen double __csample_keep__ = 0
    local __keep "__csample_keep__"

    * Get ALL variables and find indices (cencode pattern)
    unab allvars : *
    local keep_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "__csample_keep__") {
            local keep_idx = `idx'
        }
        local ++idx
    }

    * Build by-variable indices
    local by_indices ""
    if "`by'" != "" {
        foreach bv of local by {
            local idx = 1
            foreach v of local allvars {
                if ("`v'" == "`bv'") {
                    local by_indices "`by_indices' `idx'"
                    continue, break
                }
                local ++idx
            }
        }
    }

    * Build options
    local opts ""
    if `use_percent' == 1 {
        local opts "percent=`pct_value'"
    }
    else {
        local opts "count=`count'"
    }
    local opts "`opts' seed=`stata_seed'"
    if "`verbose'" != "" local opts "`opts' verbose"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Record pre-plugin time
    timer off 90
    quietly timer list 90
    local __time_preplugin = r(t90)
    timer on 90

    * Call plugin with ALL variables (cencode pattern)
    * Pass keep_idx so C code knows which variable to write to
    plugin call ctools_plugin `allvars' if `touse', ///
        "csample `threads_opt' `keep_idx' `nby' `by_indices' `opts'"

    * Drop observations where __keep == 0
    quietly drop if `__keep' == 0

    * Clean up the keep variable
    capture drop __csample_keep__

    * Finish timing
    timer off 90

    if `__do_timing' {
        quietly timer list 90
        local __time_total = r(t90)

        di as text ""
        di as text "{hline 60}"
        di as text "csample timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:         " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f _csample_time_load " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f _csample_time_sort " sec"
        }
        di as text "    Sampling:               " as result %8.4f _csample_time_sample " sec"
        di as text "    Data store:             " as result %8.4f _csample_time_store " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f _csample_time_total " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Original obs: " as result `nobs_sample'
        di as text "  Kept:         " as result %8.0f _csample_kept
        di as text "  Dropped:      " as result %8.0f _csample_dropped
        di as text "  Groups:       " as result %8.0f _csample_ngroups
        di as text "  Threads:      " as result %8.0f _csample_threads
        di as text "{hline 60}"
    }

    * Summary message
    if "`verbose'" == "" {
        local kept = _csample_kept
        local dropped = _csample_dropped
        di as text "(`dropped' observations deleted)"
    }
end
