*! version 0.9.0 26Jan2026
*! cbsample: C-accelerated bootstrap sampling for Stata
*! Drop-in replacement for bsample

program define cbsample
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    timer clear 90
    timer on 90

    * Parse syntax to match bsample exactly:
    * bsample [n] [if] [in], [strata(varlist) cluster(varlist) weight(varname)]
    * Note: Must manually extract weight() because Stata treats it specially

    * First, extract weight() option if present (before syntax parsing)
    * Use __wgtvar to avoid conflict with Stata's reserved 'weight' local
    local __wgtvar ""
    local cmdline : copy local 0

    * Use regex-like extraction for weight(varname)
    if regexm(`"`cmdline'"', "weight\(([a-zA-Z_][a-zA-Z0-9_]*)\)") {
        local __wgtvar = regexs(1)
        local cmdline = regexr(`"`cmdline'"', "weight\([a-zA-Z_][a-zA-Z0-9_]*\)", "")
    }
    else if regexm(`"`cmdline'"', "Weight\(([a-zA-Z_][a-zA-Z0-9_]*)\)") {
        local __wgtvar = regexs(1)
        local cmdline = regexr(`"`cmdline'"', "Weight\([a-zA-Z_][a-zA-Z0-9_]*\)", "")
    }
    else if regexm(`"`cmdline'"', "WEIGHT\(([a-zA-Z_][a-zA-Z0-9_]*)\)") {
        local __wgtvar = regexs(1)
        local cmdline = regexr(`"`cmdline'"', "WEIGHT\([a-zA-Z_][a-zA-Z0-9_]*\)", "")
    }

    * Now parse the rest without the weight option
    local 0 `"`cmdline'"'
    syntax [anything(name=n)] [if] [in], [STRata(varlist) CLuster(varlist) Verbose THReads(integer 0)]

    local __do_timing = ("`verbose'" != "")

    * Mark sample
    marksample touse, novarlist

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
    if `nobs_sample' == 0 {
        exit 0
    }

    * Parse n (positional argument, defaults to _N)
    if "`n'" == "" {
        local n = `nobs_sample'
    }
    else {
        capture confirm integer number `n'
        if _rc != 0 {
            di as error "cbsample: n must be a positive integer"
            exit 198
        }
        if `n' < 1 {
            di as error "cbsample: n must be a positive integer"
            exit 198
        }
    }

    * Validate weight variable exists if specified (matches bsample behavior)
    if "`__wgtvar'" != "" {
        capture confirm variable `__wgtvar'
        if _rc != 0 {
            di as error "variable `__wgtvar' not found"
            exit 111
        }
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
            di as error "cbsample: Could not load ctools plugin"
            exit 601
        }
    }

    * Generate seed from Stata's RNG (respects -set seed-)
    local stata_seed = floor(runiform() * 2147483647) + floor(runiform() * 2147483647) * 2147483648

    * Count cluster and strata variables
    local ncluster = 0
    if "`cluster'" != "" {
        local ncluster : word count `cluster'
    }

    local nstrata = 0
    if "`strata'" != "" {
        local nstrata : word count `strata'
    }

    * Create temporary frequency weight variable
    tempvar __cbsample_weight__
    quietly gen double `__cbsample_weight__' = 0

    * Get ALL variables and find indices
    unab allvars : *
    local freq_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`__cbsample_weight__'") {
            local freq_idx = `idx'
        }
        local ++idx
    }

    * Build cluster variable indices
    local cluster_indices ""
    if "`cluster'" != "" {
        foreach cv of local cluster {
            local idx = 1
            foreach v of local allvars {
                if ("`v'" == "`cv'") {
                    local cluster_indices "`cluster_indices' `idx'"
                    continue, break
                }
                local ++idx
            }
        }
    }

    * Build strata variable indices
    local strata_indices ""
    if "`strata'" != "" {
        foreach sv of local strata {
            local idx = 1
            foreach v of local allvars {
                if ("`v'" == "`sv'") {
                    local strata_indices "`strata_indices' `idx'"
                    continue, break
                }
                local ++idx
            }
        }
    }

    * Build options
    local opts "n=`n' seed=`stata_seed'"
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

    * Call plugin with ALL variables
    plugin call ctools_plugin `allvars' if `touse', ///
        "cbsample `threads_opt' `freq_idx' `ncluster' `nstrata' `cluster_indices' `strata_indices' `opts'"

    * Record post-plugin time
    timer off 90
    quietly timer list 90
    local __time_postplugin_start = r(t90)
    timer on 90

    * Handle output based on weight option (match bsample behavior exactly)
    if "`__wgtvar'" != "" {
        * weight() specified: copy weights to user's variable, keep all obs
        quietly replace `__wgtvar' = `__cbsample_weight__'
    }
    else {
        * No weight(): expand data to match bsample behavior
        * Drop observations with zero weight first
        quietly drop if `__cbsample_weight__' == 0
        * Expand based on weight (creates duplicates for weight > 1)
        quietly expand `__cbsample_weight__'
    }

    * Finish timing
    timer off 90

    if `__do_timing' {
        quietly timer list 90
        local __time_total = r(t90)
        local __time_postplugin = `__time_total' - `__time_postplugin_start'

        di as text ""
        di as text "{hline 60}"
        di as text "cbsample timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:         " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f _cbsample_time_load " sec"
        if `ncluster' > 0 | `nstrata' > 0 {
            di as text "    Sort groups:            " as result %8.4f _cbsample_time_sort " sec"
        }
        di as text "    Sampling:               " as result %8.4f _cbsample_time_sample " sec"
        di as text "    Data store:             " as result %8.4f _cbsample_time_store " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f _cbsample_time_total " sec"
        di as text "  {hline 58}"
        di as text "  Stata post-plugin:        " as result %8.4f `__time_postplugin' " sec"
        if "`__wgtvar'" == "" {
            di as text "    (drop zeros + expand)"
        }
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Original obs: " as result `nobs_sample'
        di as text "  Selected:     " as result %8.0f _cbsample_n_selected
        di as text "  Total weight: " as result %8.0f _cbsample_total_weight
        di as text "  Clusters:     " as result %8.0f _cbsample_nclusters
        di as text "  Strata:       " as result %8.0f _cbsample_nstrata
        di as text "  Threads:      " as result %8.0f _cbsample_threads
        di as text "{hline 60}"
    }
end
