*! version 1.1.0 20261016 github.com/mdroste/stata-ctools

program define cbsample
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    local __ct_elapsed90 = 0
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Parse syntax to match bsample exactly:
    * bsample [n] [if] [in], [strata(varlist) cluster(varlist) weight(varname)]
    syntax [anything(name=n)] [if] [in], [STRata(varlist) CLuster(varlist) Weight(varname) Verbose THReads(integer 0)]
    local __wgtvar `weight'
    local weight

    local __do_timing = ("`verbose'" != "")

    * Mark sample
    marksample touse, novarlist
    markout `touse' `strata', strok

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
    if `nobs_sample' == 0 {
        exit 0
    }

    * Parse n (positional argument, defaults to _N)
    if "`n'" == "" {
        local n = 0
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
        confirm numeric variable `__wgtvar'
    }

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Generate seed from Stata's RNG (respects -set seed-). Pass two exact
    * 31-bit halves: a single ~62-bit value loses digits in a macro.
    local seed_hi = string(floor(runiform() * 2147483647), "%10.0f")
    local seed_lo = string(floor(runiform() * 2147483647), "%10.0f")

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

    * Get ALL variables and find indices. The sample marker goes last, in
    * place of an if: Stata answers an if one observation at a time, while
    * the plugin reads the marker in parallel.
    unab allvars : *
    local allvars : list allvars - touse
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
    local n_s = strtrim(string(`n', "%20.0f"))
    local opts "n=`n_s' seedhi=`seed_hi' seedlo=`seed_lo'"
    if "`verbose'" != "" local opts "`opts' verbose"
    if `nobs_sample' == _N local opts "`opts' skipif"
    else local opts "`opts' touse"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Record pre-plugin time
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    * Elapsed seconds are held in caller-local timing state.
    local __time_preplugin = `__ct_elapsed90'
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Call plugin with ALL variables
    plugin call ctools_plugin `allvars' `touse', ///
        "cbsample `threads_opt' `freq_idx' `ncluster' `nstrata' `cluster_indices' `strata_indices' `opts'"

    * Record post-plugin time
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    * Elapsed seconds are held in caller-local timing state.
    local __time_postplugin_start = `__ct_elapsed90'
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

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
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

    if `__do_timing' {
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
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
