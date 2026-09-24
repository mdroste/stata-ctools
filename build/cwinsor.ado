*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

program define cwinsor
    version 14.1
    preserve
    capture noisily _cwinsor_impl `0'
    local rc = _rc
    if `rc' {
        restore
        exit `rc'
    }
    restore, not
end

program define _cwinsor_impl
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start total timing FIRST
    timer clear 90
    timer clear 91
    timer clear 92
    timer on 90

    syntax varlist(numeric) [if] [in], [Cuts(numlist min=2 max=2 >=0 <=100) P(real 1) Q(real 99) TRIM BY(varlist) SUFfix(string) PREfix(string) REPLACE Verbose THReads(integer 0)]

    local __do_timing = ("`verbose'" != "")

    * Handle percentile options
    if "`cuts'" != "" {
        local p : word 1 of `cuts'
        local q : word 2 of `cuts'
    }

    * Validate percentiles
    if `p' < 0 | `p' > 100 {
        di as error "cwinsor: lower percentile must be between 0 and 100"
        exit 198
    }
    if `q' < 0 | `q' > 100 {
        di as error "cwinsor: upper percentile must be between 0 and 100"
        exit 198
    }
    if `p' >= `q' {
        di as error "cwinsor: lower percentile (`p') must be less than upper (`q')"
        exit 198
    }

    * Validate suffix/prefix
    if "`suffix'" != "" & "`prefix'" != "" {
        di as error "cwinsor: cannot specify both suffix() and prefix()"
        exit 198
    }
    if ("`suffix'" != "" | "`prefix'" != "") & "`replace'" != "" {
        di as error "cwinsor: cannot specify replace with suffix() or prefix()"
        exit 198
    }

    * Check numeric variables
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc != 0 {
            di as error "cwinsor: `v' is not numeric"
            exit 109
        }
    }

    local nvars : word count `varlist'

    * Mark sample
    marksample touse, novarlist

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
    if `nobs_sample' == 0 {
        di as error "cwinsor: no observations"
        exit 2000
    }

    * Count by-variables (sorting now done in C for speed)
    timer on 91
    local nby = 0
    if "`by'" != "" {
        local nby : word count `by'
    }
    timer off 91

    * Handle generate new variables
    local target_varlist "`varlist'"
    if "`suffix'" != "" | "`prefix'" != "" {
        local target_varlist ""
        foreach v of local varlist {
            if "`suffix'" != "" {
                local newvar "`v'`suffix'"
            }
            else {
                local newvar "`prefix'`v'"
            }
            capture confirm variable `newvar'
            if _rc == 0 {
                di as error "cwinsor: variable `newvar' already exists"
                exit 110
            }
            quietly gen double `newvar' = `v' if `touse'
            local target_varlist "`target_varlist' `newvar'"
        }
    }

    * Checked SPI stores use positions in the plugin varlist, just like reads.
    * Only target variables need store indices; by-variables are read-only.
    local store_indices ""
    local __pos = 0
    foreach v of local target_varlist {
        local ++__pos
        local store_indices "`store_indices' `__pos'"
    }

    * Build options
    local opts "p=`p' q=`q'"
    if "`trim'" != "" local opts "`opts' trim"
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

    * Call plugin with only needed variables (target + by)
    * C uses plugin-local indices for both loading and writing
    local nvars_target : word count `target_varlist'

    plugin call ctools_plugin `target_varlist' `by' if `touse', ///
        "cwinsor `threads_opt' `nvars_target' `nby' `store_indices' `opts'"

    * Record post-plugin time
    timer on 92

    * Finish timing
    timer off 90
    timer off 92

    if `__do_timing' {
        quietly timer list 90
        local __time_total = r(t90)
        quietly timer list 91
        local __time_sort = r(t91)
        quietly timer list 92
        local __time_postplugin = r(t92)

        * Calculate pre-plugin time (excluding sort which is timed separately)
        local __time_preplugin_other = `__time_preplugin' - `__time_sort'

        di as text ""
        di as text "{hline 60}"
        di as text "cwinsor timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:"
        di as text "    Syntax/validation:      " as result %8.4f `__time_preplugin_other' " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f `__time_sort' " sec"
        }
        di as text "  {hline 58}"
        di as text "    Pre-plugin total:       " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f _cwinsor_time_load " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f _cwinsor_time_sort " sec"
        }
        di as text "    Group detection:        " as result %8.4f _cwinsor_time_groups " sec"
        di as text "    Winsorization:          " as result %8.4f _cwinsor_time_winsor " sec"
        di as text "    Data store:             " as result %8.4f _cwinsor_time_store " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f _cwinsor_time_total " sec"
        di as text "  {hline 58}"
        di as text "  Stata post-plugin:        " as result %8.4f `__time_postplugin' " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Dataset: " as result `nobs_sample' as text " obs x " as result `nvars_target' as text " vars"
        di as text "  Groups:  " as result %8.0f _cwinsor_ngroups
        di as text "  Mode:    " as result cond("`trim'"!="", "trim", "winsorize") as text " at [" as result `p' as text ", " as result `q' as text "] percentiles"
        di as text "  Threads: " as result %8.0f _cwinsor_threads
        di as text "{hline 60}"
    }

    * Summary
    if "`verbose'" == "" {
        local mode_str = cond("`trim'" != "", "trimmed", "winsorized")
        di as text "(`nvars_target' variables `mode_str' at [`p', `q'] percentiles)"
    }
end
