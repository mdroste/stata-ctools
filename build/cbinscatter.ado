*! version 1.0.0 09Jan2026
*! cbinscatter: C-accelerated binned scatter plots
*! Part of the ctools suite
*!
*! Description:
*!   High-performance replacement for binscatter using a C plugin
*!   with optimized data processing and optional HDFE residualization.
*!
*! Syntax:
*!   cbinscatter yvar xvar [if] [in] [weight], [options]
*!
*! Options:
*!   nquantiles(#)       - Number of bins (default 20)
*!   controls(varlist)   - Control variables to partial out
*!   absorb(varlist)     - Fixed effects to absorb
*!   by(varname)         - Separate series by group
*!   linetype(string)    - Line fit type: none, linear, quadratic, cubic
*!   discrete            - Treat x as discrete (one bin per unique value)
*!   genxq(varname)      - Generate bin assignment variable
*!   savedata(filename)  - Save bin data to file
*!   reportreg           - Report underlying regression
*!   nograph             - Suppress graph
*!   verbose             - Display progress information
*!   timeit              - Display timing breakdown

program define cbinscatter, eclass sortpreserve
    version 14.0

    syntax varlist(min=2 max=2 fv) [aw fw pw iw] [if] [in], ///
        [NQuantiles(integer 20)] ///
        [Controls(varlist fv)] ///
        [Absorb(varlist)] ///
        [BY(varname)] ///
        [LINEtype(string)] ///
        [DISCrete] ///
        [GENxq(name)] ///
        [SAVEdata(string)] ///
        [REPORTreg] ///
        [NOgraph] ///
        [Verbose] ///
        [TIMEit] ///
        /* Graph options */ ///
        [TITLE(string)] ///
        [YTItle(string)] ///
        [XTItle(string)] ///
        [LEGEND(string asis)] ///
        [COLors(string)] ///
        [MSYMbols(string)] ///
        [MLAbels(string)] ///
        [*]

    * Mark sample - include all relevant variables
    marksample touse
    if "`controls'" != "" {
        markout `touse' `controls'
    }
    if "`absorb'" != "" {
        markout `touse' `absorb'
    }
    if "`by'" != "" {
        markout `touse' `by'
    }

    * Parse weight
    local weight_var ""
    local weight_type = 0   // 0=none, 1=aweight, 2=fweight, 3=pweight, 4=iweight
    if "`weight'" != "" {
        local weight_var "`exp'"
        local weight_var = subinstr("`weight_var'", "=", "", .)
        local weight_var = trim("`weight_var'")

        if "`weight'" == "aweight" {
            local weight_type = 1
        }
        else if "`weight'" == "fweight" {
            local weight_type = 2
        }
        else if "`weight'" == "pweight" {
            local weight_type = 3
        }
        else if "`weight'" == "iweight" {
            local weight_type = 4
        }

        markout `touse' `weight_var'

        * Check for non-positive weights
        qui count if `weight_var' <= 0 & `touse'
        if r(N) > 0 {
            di as error "cbinscatter: weights must be positive"
            exit 198
        }
    }

    * Parse variable list
    gettoken depvar xvar : varlist
    local depvar = trim("`depvar'")
    local xvar = trim("`xvar'")

    * Parse linetype option
    local linetype_num = 1   // Default to linear
    if "`linetype'" == "" {
        local linetype = "linear"
        local linetype_num = 1
    }
    else {
        local linetype_lower = lower("`linetype'")
        if "`linetype_lower'" == "none" | "`linetype_lower'" == "n" {
            local linetype_num = 0
        }
        else if "`linetype_lower'" == "linear" | "`linetype_lower'" == "line" | "`linetype_lower'" == "lfit" {
            local linetype_num = 1
        }
        else if "`linetype_lower'" == "quadratic" | "`linetype_lower'" == "qfit" {
            local linetype_num = 2
        }
        else if "`linetype_lower'" == "cubic" {
            local linetype_num = 3
        }
        else if "`linetype_lower'" == "connect" {
            local linetype_num = 0  // No fit, just connected points
        }
        else {
            di as error "cbinscatter: invalid linetype option"
            exit 198
        }
    }

    * Validate nquantiles
    if `nquantiles' < 2 | `nquantiles' > 1000 {
        di as error "cbinscatter: nquantiles must be between 2 and 1000"
        exit 198
    }

    * Count observations
    qui count if `touse'
    local nobs = r(N)
    if `nobs' == 0 {
        di as error "cbinscatter: no observations"
        exit 2000
    }
    if `nobs' < `nquantiles' & "`discrete'" == "" {
        di as text "(warning: fewer observations than bins, adjusting nquantiles)"
        local nquantiles = `nobs'
    }

    * Count control and absorb variables
    local num_controls = 0
    if "`controls'" != "" {
        local num_controls : word count `controls'
    }
    local num_absorb = 0
    if "`absorb'" != "" {
        local num_absorb : word count `absorb'
    }

    * Count by-groups if specified
    local num_by_groups = 1
    if "`by'" != "" {
        qui tab `by' if `touse'
        local num_by_groups = r(r)
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
            di as error "cbinscatter: Could not load ctools plugin"
            exit 601
        }
    }

    * Set up parameters via Stata scalars
    scalar __cbinscatter_nquantiles = `nquantiles'
    scalar __cbinscatter_linetype = `linetype_num'
    scalar __cbinscatter_compute_se = 0
    scalar __cbinscatter_discrete = ("`discrete'" != "")
    scalar __cbinscatter_has_controls = (`num_controls' > 0)
    scalar __cbinscatter_num_controls = `num_controls'
    scalar __cbinscatter_has_absorb = (`num_absorb' > 0)
    scalar __cbinscatter_num_absorb = `num_absorb'
    scalar __cbinscatter_has_by = ("`by'" != "")
    scalar __cbinscatter_has_weights = (`weight_type' > 0)
    scalar __cbinscatter_weight_type = `weight_type'
    scalar __cbinscatter_verbose = ("`verbose'" != "")
    scalar __cbinscatter_reportreg = ("`reportreg'" != "")
    scalar __cbinscatter_maxiter = 10000
    scalar __cbinscatter_tolerance = 1e-8

    * Display info if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cbinscatter: C-Accelerated Binned Scatter"
        di as text "{hline 60}"
        di as text "Y variable: " as result "`depvar'"
        di as text "X variable: " as result "`xvar'"
        if `num_controls' > 0 {
            di as text "Controls:   " as result "`controls'"
        }
        if `num_absorb' > 0 {
            di as text "Absorb:     " as result "`absorb'"
        }
        if "`by'" != "" {
            di as text "By:         " as result "`by'"
        }
        di as text "Bins:       " as result `nquantiles'
        di as text "Obs:        " as result `nobs'
        di as text "{hline 60}"
        di ""
    }

    * Create matrices for C plugin to fill
    local max_bins = `nquantiles'
    if "`discrete'" != "" {
        * For discrete, could have many more bins - use a reasonable max
        local max_bins = min(`nobs', 500)
    }
    matrix __cbinscatter_bins = J(`max_bins' * `num_by_groups', 7, .)
    if `linetype_num' > 0 {
        matrix __cbinscatter_coefs = J(`num_by_groups', 4, .)
        matrix __cbinscatter_fit_stats = J(`num_by_groups', 2, .)
    }

    * Build varlist for plugin: y x [controls] [absorb] [by] [weight]
    local plugin_varlist `depvar' `xvar'
    if "`controls'" != "" {
        local plugin_varlist `plugin_varlist' `controls'
    }
    if "`absorb'" != "" {
        local plugin_varlist `plugin_varlist' `absorb'
    }
    if "`by'" != "" {
        * Convert by variable to numeric if string
        capture confirm numeric variable `by'
        if _rc != 0 {
            tempvar by_numeric
            qui egen `by_numeric' = group(`by')
            local by_orig "`by'"
            local by "`by_numeric'"
        }
        local plugin_varlist `plugin_varlist' `by'
    }
    if `weight_type' > 0 {
        local plugin_varlist `plugin_varlist' `weight_var'
    }

    * Record start time
    timer clear 99
    timer on 99

    * Call the C plugin
    capture noisily plugin call ctools_plugin `plugin_varlist' if `touse', ///
        "cbinscatter compute_bins"

    local plugin_rc = _rc
    if `plugin_rc' {
        * Clean up scalars on error
        capture scalar drop __cbinscatter_nquantiles __cbinscatter_linetype
        capture scalar drop __cbinscatter_compute_se __cbinscatter_discrete
        capture scalar drop __cbinscatter_has_controls __cbinscatter_num_controls
        capture scalar drop __cbinscatter_has_absorb __cbinscatter_num_absorb
        capture scalar drop __cbinscatter_has_by __cbinscatter_has_weights
        capture scalar drop __cbinscatter_weight_type __cbinscatter_verbose
        capture scalar drop __cbinscatter_reportreg __cbinscatter_maxiter
        capture scalar drop __cbinscatter_tolerance
        capture matrix drop __cbinscatter_bins
        capture matrix drop __cbinscatter_coefs
        capture matrix drop __cbinscatter_fit_stats
        di as error "cbinscatter: plugin error (rc=`plugin_rc')"
        exit `plugin_rc'
    }

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Retrieve results from scalars
    local N_used = __cbinscatter_N
    local N_dropped = __cbinscatter_N_dropped
    local actual_num_groups = __cbinscatter_num_groups

    * Display timing if requested
    if "`timeit'" != "" | "`verbose'" != "" {
        di as text ""
        di as text "Timing breakdown:"
        di as text "  Load:       " as result %8.4f __cbinscatter_time_load " sec"
        di as text "  Residualize:" as result %8.4f __cbinscatter_time_resid " sec"
        di as text "  Bin comp:   " as result %8.4f __cbinscatter_time_bins " sec"
        di as text "  Line fit:   " as result %8.4f __cbinscatter_time_fit " sec"
        di as text "  Total:      " as result %8.4f __cbinscatter_time_total " sec"
    }

    * Display summary
    if `N_dropped' > 0 {
        di as text "(dropped " as result `N_dropped' as text " observations due to missing values)"
    }

    * Save bin data if requested
    if "`savedata'" != "" {
        preserve
        quietly {
            clear
            svmat __cbinscatter_bins, names(col)
            rename c1 by_group
            rename c2 bin_id
            rename c3 x_mean
            rename c4 y_mean
            rename c5 x_se
            rename c6 y_se
            rename c7 n_obs
            * Drop empty rows
            drop if missing(bin_id)
            save "`savedata'", replace
        }
        di as text "(bin data saved to " as result "`savedata'" as text ")"
        restore
    }

    * Generate bin assignment variable if requested
    if "`genxq'" != "" {
        * This would require additional plugin work to return per-obs bin IDs
        * For now, just display a message
        di as text "(genxq option not yet implemented)"
    }

    * Create graph if not suppressed
    if "`nograph'" == "" {
        * Prepare graph data
        preserve
        quietly {
            clear
            svmat __cbinscatter_bins, names(col)
            rename c1 by_group
            rename c2 bin_id
            rename c3 x_mean
            rename c4 y_mean
            rename c5 x_se
            rename c6 y_se
            rename c7 n_obs
            * Drop empty rows
            drop if missing(bin_id)
        }

        * Build scatter command
        local scatter_cmd ""
        local line_cmd ""
        local legend_order ""

        if `actual_num_groups' == 1 {
            * Single group - simple scatter
            local scatter_cmd "(scatter y_mean x_mean)"

            * Add line fit if requested
            if `linetype_num' > 0 {
                * Get fit coefficients
                local coef_const = __cbinscatter_coefs[1, 1]
                local coef_linear = __cbinscatter_coefs[1, 2]

                if `linetype_num' == 1 {
                    * Linear fit
                    qui gen fitted = `coef_const' + `coef_linear' * x_mean
                    local line_cmd "(line fitted x_mean, sort)"
                }
                else if `linetype_num' >= 2 {
                    * Polynomial fit
                    local coef_quad = __cbinscatter_coefs[1, 3]
                    qui gen fitted = `coef_const' + `coef_linear' * x_mean + `coef_quad' * x_mean^2
                    if `linetype_num' >= 3 {
                        local coef_cubic = __cbinscatter_coefs[1, 4]
                        qui replace fitted = fitted + `coef_cubic' * x_mean^3
                    }
                    local line_cmd "(line fitted x_mean, sort)"
                }
            }
        }
        else {
            * Multiple by-groups
            local scatter_cmd ""
            local line_cmd ""
            local idx = 1

            forval g = 1/`actual_num_groups' {
                if `g' == 1 {
                    local scatter_cmd "(scatter y_mean x_mean if by_group == `g')"
                }
                else {
                    local scatter_cmd "`scatter_cmd' (scatter y_mean x_mean if by_group == `g')"
                }

                * Add line fit for this group if requested
                if `linetype_num' > 0 {
                    local coef_const = __cbinscatter_coefs[`g', 1]
                    local coef_linear = __cbinscatter_coefs[`g', 2]

                    if `linetype_num' == 1 {
                        qui gen fitted`g' = `coef_const' + `coef_linear' * x_mean if by_group == `g'
                        local line_cmd "`line_cmd' (line fitted`g' x_mean if by_group == `g', sort)"
                    }
                    else if `linetype_num' >= 2 {
                        local coef_quad = __cbinscatter_coefs[`g', 3]
                        qui gen fitted`g' = `coef_const' + `coef_linear' * x_mean + `coef_quad' * x_mean^2 if by_group == `g'
                        if `linetype_num' >= 3 {
                            local coef_cubic = __cbinscatter_coefs[`g', 4]
                            qui replace fitted`g' = fitted`g' + `coef_cubic' * x_mean^3 if by_group == `g'
                        }
                        local line_cmd "`line_cmd' (line fitted`g' x_mean if by_group == `g', sort)"
                    }
                }
            }
        }

        * Set default titles if not specified
        if `"`title'"' == "" {
            local title "Binned Scatter Plot"
        }
        if `"`ytitle'"' == "" {
            local ytitle "`depvar'"
        }
        if `"`xtitle'"' == "" {
            local xtitle "`xvar'"
        }

        * Execute graph
        twoway `scatter_cmd' `line_cmd', ///
            title(`"`title'"') ///
            ytitle(`"`ytitle'"') ///
            xtitle(`"`xtitle'"') ///
            `legend' ///
            `options'

        restore
    }

    * Post e() results
    ereturn clear
    ereturn scalar N = `N_used'
    ereturn scalar N_dropped = `N_dropped'
    ereturn scalar nquantiles = `nquantiles'
    ereturn scalar num_groups = `actual_num_groups'
    ereturn matrix bindata = __cbinscatter_bins
    if `linetype_num' > 0 {
        ereturn matrix coefs = __cbinscatter_coefs
        ereturn matrix fit_stats = __cbinscatter_fit_stats
    }
    ereturn local depvar "`depvar'"
    ereturn local xvar "`xvar'"
    if "`controls'" != "" {
        ereturn local controls "`controls'"
    }
    if "`absorb'" != "" {
        ereturn local absorb "`absorb'"
    }
    if "`by'" != "" {
        ereturn local by "`by'"
    }
    ereturn local linetype "`linetype'"
    ereturn local cmd "cbinscatter"
    ereturn local cmdline "cbinscatter `0'"

    * Clean up scalars
    capture scalar drop __cbinscatter_nquantiles __cbinscatter_linetype
    capture scalar drop __cbinscatter_compute_se __cbinscatter_discrete
    capture scalar drop __cbinscatter_has_controls __cbinscatter_num_controls
    capture scalar drop __cbinscatter_has_absorb __cbinscatter_num_absorb
    capture scalar drop __cbinscatter_has_by __cbinscatter_has_weights
    capture scalar drop __cbinscatter_weight_type __cbinscatter_verbose
    capture scalar drop __cbinscatter_reportreg __cbinscatter_maxiter
    capture scalar drop __cbinscatter_tolerance
    capture scalar drop __cbinscatter_N __cbinscatter_N_dropped __cbinscatter_num_groups
    capture scalar drop __cbinscatter_time_load __cbinscatter_time_resid
    capture scalar drop __cbinscatter_time_bins __cbinscatter_time_fit __cbinscatter_time_total

    timer clear 99
end
