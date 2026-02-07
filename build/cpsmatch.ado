*! version 0.9.1 06Feb2026
*! cpsmatch: C-accelerated propensity score matching for Stata
*! Part of the ctools suite
*!
*! Implements psmatch2-compatible propensity score matching with:
*!   - Nearest neighbor matching (with/without replacement)
*!   - Radius/caliper matching
*!   - Kernel matching (Epanechnikov, Gaussian, biweight, uniform, tricube)
*!   - Common support restriction
*!
*! Syntax: cpsmatch treatvar [varlist] [if] [in] [, options]
*!
*! Options:
*!   OUTcome(varname)     - outcome variable for ATT calculation
*!   Pscore(varname)      - use existing propensity score variable
*!   Logit | Probit       - estimation method (default: probit)
*!   Neighbor(#)          - number of nearest neighbors (default: 1)
*!   Caliper(#)           - caliper width
*!   Radius               - radius matching
*!   Kernel               - kernel matching
*!   Kerneltype(string)   - kernel type: epan, normal, biweight, uniform, tricube
*!   Bwidth(#)            - bandwidth for kernel matching
*!   Common               - impose common support
*!   Noreplacement        - match without replacement
*!   Ties                 - include tied matches
*!   Descending           - match in descending pscore order
*!   Verbose              - display timing breakdown

program define cpsmatch, rclass sortpreserve
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    capture noisily syntax varlist(min=1 numeric) [if] [in], ///
        [OUTcome(varname numeric) Pscore(varname numeric) ///
         Logit Probit ///
         Neighbor(integer 1) Caliper(real 0) Radius Kernel ///
         Kerneltype(string) Bwidth(real 0.06) ///
         Common NOREPlacement Ties Descending ///
         Verbose THReads(integer 0)]
    if _rc == 109 {
        exit 2000
    }
    else if _rc != 0 {
        exit _rc
    }

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Split varlist into treatment variable and covariates
    gettoken depvar treatvar : varlist
    * depvar = treatment variable (first var)
    * treatvar = covariates (remaining vars)

    * Check treatment is binary
    quietly tab `depvar' `if' `in'
    if r(r) != 2 {
        di as error "cpsmatch: treatment variable must be binary (0/1)"
        exit 2000
    }

    * Determine matching method
    local method = 0  /* 0=nearest, 1=radius, 2=kernel, 3=mahal */
    if "`radius'" != "" {
        local method = 1
        if `caliper' == 0 {
            local caliper = 0.25  /* Default radius */
        }
    }
    if "`kernel'" != "" {
        local method = 2
    }

    * Kernel type
    local kernelcode = 0  /* 0=epan, 1=normal, 2=biweight, 3=uniform, 4=tricube */
    if "`kerneltype'" != "" {
        if "`kerneltype'" == "normal" | "`kerneltype'" == "gaussian" {
            local kernelcode = 1
        }
        else if "`kerneltype'" == "biweight" {
            local kernelcode = 2
        }
        else if "`kerneltype'" == "uniform" {
            local kernelcode = 3
        }
        else if "`kerneltype'" == "tricube" {
            local kernelcode = 4
        }
    }

    * Default estimation method
    if "`logit'" == "" & "`probit'" == "" {
        local probit "probit"
    }

    * With/without replacement
    local with_replace = 1
    if "`noreplacement'" != "" {
        local with_replace = 0
    }

    * Ties handling
    local handle_ties = 0
    if "`ties'" != "" {
        local handle_ties = 1
    }

    * Common support
    local common_support = 0
    if "`common'" != "" {
        local common_support = 1
    }

    * Descending order
    local desc_order = 0
    if "`descending'" != "" {
        local desc_order = 1
    }

    * Verbose timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        timer clear 90
        timer clear 91
        timer clear 92
        timer clear 94
        timer on 90
        timer on 91
    }

    * =========================================================================
    * LOAD PLUGIN
    * =========================================================================

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
            di as error "cpsmatch: Could not load ctools plugin"
            exit 601
        }
    }

    * =========================================================================
    * ESTIMATE PROPENSITY SCORE (if not provided)
    * =========================================================================

    * Drop existing output variables if they exist
    capture drop _weight
    capture drop _id
    capture drop _support
    capture drop _pscore
    capture drop _treated
    capture drop _nn

    if "`pscore'" != "" {
        * Use provided propensity score
        quietly gen double _pscore = `pscore' `if' `in'
    }
    else if "`treatvar'" != "" {
        * Estimate propensity score using covariates (show output like psmatch2)
        if "`logit'" != "" {
            logit `depvar' `treatvar' `if' `in', nolog
        }
        else {
            probit `depvar' `treatvar' `if' `in', nolog
        }
        quietly predict double _pscore `if' `in', pr
    }
    else {
        di as error "cpsmatch: must specify either pscore() or covariates for estimation"
        exit 198
    }

    * =========================================================================
    * CREATE OUTPUT VARIABLES
    * =========================================================================

    * Create output variables
    quietly gen double _weight = .
    quietly gen double _id = .
    quietly gen double _support = .
    quietly gen byte _treated = `depvar' `if' `in'
    quietly gen int _nn = .

    * Get variable indices - include ALL variables (like cencode)
    unab allvars : *

    local treat_idx = 0
    local pscore_idx = 0
    local outcome_idx = 0
    local weight_idx = 0
    local match_idx = 0
    local support_idx = 0

    local idx = 1
    foreach v of local allvars {
        if "`v'" == "`depvar'" {
            local treat_idx = `idx'
        }
        if "`v'" == "_pscore" {
            local pscore_idx = `idx'
        }
        if "`outcome'" != "" & "`v'" == "`outcome'" {
            local outcome_idx = `idx'
        }
        if "`v'" == "_weight" {
            local weight_idx = `idx'
        }
        if "`v'" == "_id" {
            local match_idx = `idx'
        }
        if "`v'" == "_support" {
            local support_idx = `idx'
        }
        local ++idx
    }

    * Count observations
    quietly count `if' `in'
    local nobs = r(N)

    if `nobs' == 0 {
        di as error "cpsmatch: no observations"
        exit 2000
    }

    * Build threads option
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    if `__do_timing' {
        timer off 91
        timer on 92
    }

    * =========================================================================
    * CALL C PLUGIN
    * =========================================================================

    * Build argument string
    * Format: treat_idx pscore_idx outcome_idx weight_idx match_idx support_idx
    *         nobs method neighbor caliper common kernel bwidth with_replace ties desc verbose

    local args "`treat_idx' `pscore_idx' `outcome_idx' `weight_idx' `match_idx' `support_idx'"
    local args "`args' `nobs' `method' `neighbor' `caliper' `common_support'"
    local args "`args' `kernelcode' `bwidth' `with_replace' `handle_ties' `desc_order'"

    * Call plugin with ALL variables (like cencode does)
    plugin call ctools_plugin `allvars' `if' `in', "cpsmatch `threads_code' `args'"

    if `__do_timing' {
        timer off 92
    }

    * =========================================================================
    * POST-PROCESSING
    * =========================================================================

    * Update _nn with neighbor count from matching
    quietly replace _nn = `neighbor' if _treated == 1 & _support == 1

    * =========================================================================
    * CALCULATE ATT (if outcome specified)
    * =========================================================================

    if `__do_timing' {
        timer on 94
    }

    local att = .
    local att_se = .
    local att_t = .

    if "`outcome'" != "" {
        * Calculate weighted mean difference
        quietly summarize `outcome' if _treated == 1 & _support == 1 [aw=_weight]
        local y1_bar = r(mean)
        local n1 = r(N)

        quietly summarize `outcome' if _treated == 0 & _support == 1 [aw=_weight]
        local y0_bar = r(mean)
        local n0 = r(N)

        if `n1' > 0 & `n0' > 0 {
            local att = `y1_bar' - `y0_bar'

            * Approximate SE using matched sample variance
            quietly summarize `outcome' if _treated == 1 & _support == 1 [aw=_weight]
            local var1 = r(Var)

            quietly summarize `outcome' if _treated == 0 & _support == 1 [aw=_weight]
            local var0 = r(Var)

            local att_se = sqrt(`var1'/`n1' + `var0'/`n0')
            if `att_se' > 0 {
                local att_t = `att' / `att_se'
            }
        }
    }

    * =========================================================================
    * DISPLAY RESULTS (psmatch2-compatible format)
    * =========================================================================

    * ATT results table (psmatch2 format)
    if "`outcome'" != "" & `att' != . {
        * Calculate unmatched difference
        quietly summarize `outcome' if _treated == 1, meanonly
        local y1_unmatched = r(mean)
        quietly summarize `outcome' if _treated == 0, meanonly
        local y0_unmatched = r(mean)
        local unmatched_diff = `y1_unmatched' - `y0_unmatched'

        * Calculate unmatched SE and t-stat
        quietly summarize `outcome' if _treated == 1
        local var1_u = r(Var)
        local n1_u = r(N)
        quietly summarize `outcome' if _treated == 0
        local var0_u = r(Var)
        local n0_u = r(N)
        local unmatched_se = sqrt(`var1_u'/`n1_u' + `var0_u'/`n0_u')
        local unmatched_t = `unmatched_diff' / `unmatched_se'

        di as text "{hline 78}"
        di as text "        Variable     Sample |    Treated     Controls   Difference         S.E.   T-stat"
        di as text "{hline 28}+{hline 59}"
        di as text %16s "`outcome'" "  Unmatched |" as result %11.6g `y1_unmatched' %12.6g `y0_unmatched' %12.6g `unmatched_diff' %12.6g `unmatched_se' %9.2f `unmatched_t'
        di as text "                        ATT |" as result %11.6g `y1_bar' %12.6g `y0_bar' %12.6g `att' %12.6g `att_se' %9.2f `att_t'
        di as text "{hline 28}+{hline 59}"
        di as text "Note: S.E. does not take into account that the propensity score is estimated."
    }

    if `__do_timing' {
        timer off 94
    }

    * Treatment assignment table (psmatch2 format)
    di as text ""
    quietly count if _treated == 0 & _support == 1
    local n_untreated_on = r(N)
    quietly count if _treated == 0 & _support == 0
    local n_untreated_off = r(N)
    quietly count if _treated == 1 & _support == 1
    local n_treated_on = r(N)
    quietly count if _treated == 1 & _support == 0
    local n_treated_off = r(N)
    local n_untreated_total = `n_untreated_on' + `n_untreated_off'
    local n_treated_total = `n_treated_on' + `n_treated_off'
    local n_on_total = `n_untreated_on' + `n_treated_on'
    local n_off_total = `n_untreated_off' + `n_treated_off'
    local n_total = `n_on_total' + `n_off_total'

    di as text "           | cpsmatch:"
    di as text " cpsmatch: |   Common"
    di as text " Treatment |  support"
    if `n_off_total' > 0 {
        di as text "assignment | On suppor  Off suppo |     Total"
        di as text "{hline 11}+{hline 22}+{hline 10}"
        di as text " Untreated |" as result %10.0f `n_untreated_on' %10.0f `n_untreated_off' as text " |" as result %10.0f `n_untreated_total'
        di as text "   Treated |" as result %10.0f `n_treated_on' %10.0f `n_treated_off' as text " |" as result %10.0f `n_treated_total'
        di as text "{hline 11}+{hline 22}+{hline 10}"
        di as text "     Total |" as result %10.0f `n_on_total' %10.0f `n_off_total' as text " |" as result %10.0f `n_total'
    }
    else {
        di as text "assignment | On suppor |     Total"
        di as text "{hline 11}+{hline 11}+{hline 10}"
        di as text " Untreated |" as result %10.0f `n_untreated_on' as text " |" as result %10.0f `n_untreated_total'
        di as text "   Treated |" as result %10.0f `n_treated_on' as text " |" as result %10.0f `n_treated_total'
        di as text "{hline 11}+{hline 11}+{hline 10}"
        di as text "     Total |" as result %10.0f `n_on_total' as text " |" as result %10.0f `n_total'
    }

    di as text ""

    * =========================================================================
    * TIMING OUTPUT (verbose)
    * =========================================================================

    if `__do_timing' {
        timer off 90

        quietly timer list 90
        local __time_total = r(t90)
        quietly timer list 91
        local __time_preplugin = r(t91)
        quietly timer list 92
        local __time_plugin = r(t92)
        quietly timer list 94
        local __time_att = r(t94)

        di as text "{hline 55}"
        di as text "cpsmatch timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Load data:              " as result %8.4f _cpsmatch_time_load " sec"
        di as text "    Group separation:       " as result %8.4f _cpsmatch_time_setup " sec"
        di as text "    Sort/index build:       " as result %8.4f _cpsmatch_time_sort " sec"
        di as text "    Matching:               " as result %8.4f _cpsmatch_time_match " sec"
        di as text "    Store results:          " as result %8.4f _cpsmatch_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cpsmatch_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Post-plugin (ATT):      " as result %8.4f `__time_att' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__time_att'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Thread diagnostics
        capture local __threads_max = _cpsmatch_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cpsmatch_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * =========================================================================
    * RETURN RESULTS
    * =========================================================================

    return scalar N = `nobs'
    return scalar n_treated = _cpsmatch_n_treated
    return scalar n_controls = _cpsmatch_n_controls
    return scalar n_matched = _cpsmatch_n_matched
    return scalar n_off_support = _cpsmatch_n_off_support
    return scalar common_min = _cpsmatch_common_min
    return scalar common_max = _cpsmatch_common_max

    if "`outcome'" != "" & `att' != . {
        return scalar att = `att'
        return scalar att_se = `att_se'
        return scalar att_t = `att_t'
    }

    return local method = cond(`method'==0, "nearest", cond(`method'==1, "radius", "kernel"))
    return local treatvar "`depvar'"
    if "`outcome'" != "" {
        return local outcome "`outcome'"
    }
end
