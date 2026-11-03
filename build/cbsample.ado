*! version 1.2.3 20261030 github.com/mdroste/stata-ctools

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
    * bsample [exp] [if] [in], [strata(varlist) cluster(varlist) weight(varname)]
    syntax [anything(name=n)] [if] [in], [STRata(varlist) CLuster(varlist) Weight(varname) Verbose THReads(integer 0)]
    local __wgtvar `weight'
    local weight
    _ctools_cleanup cbsample
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }

    local __do_timing = ("`verbose'" != "")

    * With no rows, native bsample distinguishes its fast no-option path
    * from explicit count/if/strata/cluster syntax. Let it validate those
    * empty-data cases; no random draws or C transfers are required.
    if _N == 0 {
        local empty_opts
        if "`strata'" != "" local empty_opts `empty_opts' strata(`strata')
        if "`cluster'" != "" local empty_opts `empty_opts' cluster(`cluster')
        if "`__wgtvar'" != "" local empty_opts `empty_opts' weight(`__wgtvar')
        quietly bsample `n' `if' `in', `empty_opts'
        foreach result in nobs n_selected total_weight nclusters nstrata {
            scalar _cbsample_`result' = 0
        }
        scalar _cbsample_threads = 1
        exit
    }

    local n_exp `"`n'"'

    * Validate weight variable exists if specified (matches bsample behavior)
    if "`__wgtvar'" != "" {
        confirm numeric variable `__wgtvar'
    }

    * A stratified draw with no eligible strata has zero weight. An empty
    * unstratified selection in nonempty data is an error instead.
    marksample touse, novarlist
    markout `touse' `strata', strok
    * bsample leaves the caller's r() in place.
    tempname __rhold
    _return hold `__rhold'
    quietly count if `touse'
    local nobs_sample = r(N)
    _return restore `__rhold'
    * The resample size is an expression, evaluated as bsample does: on the
    * sample, within each stratum under strata() (so _N counts the stratum)
    if `"`n_exp'"' == "" | `nobs_sample' == 0 local n = 0
    else _cbsample_size n : `"`n_exp'"' `touse' "`strata'"
    if `nobs_sample' == 0 {
        if "`strata'" == "" {
            di as error "cbsample: no observations in sampling sample"
            exit 198
        }
        if "`__wgtvar'" != "" quietly replace `__wgtvar' = 0
        else quietly drop if 1
        foreach result in nobs n_selected total_weight nclusters nstrata {
            scalar _cbsample_`result' = 0
        }
        scalar _cbsample_threads = 1
        exit
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

* The resample size bsample takes from [exp]: constant (within each stratum
* under strata()), an integer, and greater than zero. Stata evaluates _N in
* the sample (strata: in the stratum). Sizes that differ across strata are
* not supported.
program define _cbsample_size
    args macname colon exp touse strata
    capture confirm number `exp'
    if !_rc {
        local values `exp'
        local constant 1
    }
    else if "`strata'" == "" {
        tempvar nsamp
        quietly generate double `nsamp' = (`exp') if `touse'
        mata: __cb_v = st_data(., st_local("nsamp"), st_local("touse")); ///
            st_local("constant", strofreal(all(__cb_v :== __cb_v[1]))); ///
            st_local("values", strofreal(__cb_v[1], "%21.0g"))
        mata: mata drop __cb_v
    }
    else {
        preserve
        quietly keep if `touse'
        tempvar nsamp first
        quietly bysort `strata': generate double `nsamp' = (`exp')
        quietly by `strata': generate byte `first' = _n == 1
        quietly by `strata': replace `first' = . if `nsamp' != `nsamp'[1]
        quietly count if missing(`first')
        local constant = r(N) == 0
        mata: __cb_v = uniqrows(st_data(., st_local("nsamp"), st_local("first"))); ///
            st_local("values", invtokens(strofreal(__cb_v', "%21.0g")))
        mata: mata drop __cb_v
        restore
        if !`constant' {
            di as error "expression: `exp': is not constant within strata"
            exit 198
        }
    }
    if !`constant' {
        di as error "expression: `exp': is not constant"
        exit 198
    }
    if `: word count `values'' > 1 {
        di as error "cbsample: resample sizes that differ across strata are not supported; use bsample"
        exit 198
    }
    local value = real("`values'")
    if !(abs(`value' - int(`value')) < 1e-7) {
        di as error "expression: `exp': is not an integer"
        exit 198
    }
    if !(`value' > 0) {
        di as error "resample size must be greater than zero"
        exit 498
    }
    local value = round(`value')
    c_local `macname' `value'
end
