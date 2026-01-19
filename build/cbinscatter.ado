*! version 1.0.0 17Jan2026
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
        [METHod(string)] ///
        /* Graph options */ ///
        [TITLE(string)] ///
        [YTItle(string)] ///
        [XTItle(string)] ///
        [LEGEND(string asis)] ///
        [COLors(string)] ///
        [MSYMbols(string)] ///
        [MLAbels(string)] ///
        [*]

    * =========================================================================
    * UPFRONT VALIDATION - check all options before any data manipulation
    * This ensures errors don't leave the dataset in a bad state
    * =========================================================================

    * Clean up any leftover frame from a previous failed run
    capture frame drop _ctools_graph_frame

    * Validate nquantiles early
    if `nquantiles' < 2 | `nquantiles' > 1000 {
        di as error "cbinscatter: nquantiles must be between 2 and 1000"
        exit 198
    }

    * Validate linetype option early
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
            di as error "cbinscatter: invalid linetype() option: `linetype'"
            di as error "  valid options: none, linear, qfit, cubic, connect"
            exit 198
        }
    }

    * Validate method option early (Cattaneo et al. "On Binscatter")
    * 0 = classic (default): residualize both Y and X, bin on residualized X
    * 1 = binsreg: bin on raw X, residualize Y only (conditional means)
    local method_num = 0   // Default to classic
    if "`method'" != "" {
        local method_lower = lower("`method'")
        if "`method_lower'" == "classic" | "`method_lower'" == "c" | "`method_lower'" == "residualize" {
            local method_num = 0
        }
        else if "`method_lower'" == "binsreg" | "`method_lower'" == "b" | "`method_lower'" == "conditional" {
            local method_num = 1
        }
        else {
            di as error "cbinscatter: invalid method() option: `method'"
            di as error "  valid options: classic, binsreg"
            exit 198
        }
    }

    * Validate savedata filename if specified (check directory exists)
    if "`savedata'" != "" {
        * Extract directory path
        local savedir = ""
        local lastslash = max(strpos(reverse("`savedata'"), "/"), strpos(reverse("`savedata'"), "\"))
        if `lastslash' > 0 {
            local savedir = substr("`savedata'", 1, strlen("`savedata'") - `lastslash')
        }
        if "`savedir'" != "" {
            capture confirm file "`savedir'/."
            if _rc {
                di as error "cbinscatter: directory does not exist for savedata(): `savedir'"
                exit 601
            }
        }
    }

    * Validate genxq variable name if specified
    if "`genxq'" != "" {
        capture confirm new variable `genxq'
        if _rc {
            di as error "cbinscatter: genxq() variable already exists: `genxq'"
            exit 110
        }
    }

    * =========================================================================
    * END UPFRONT VALIDATION - now proceed with data processing
    * =========================================================================

    * Parse weight (before marksample to include weight var)
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
    }

    * Parse variable list
    gettoken depvar xvar : varlist
    local depvar = trim("`depvar'")
    local xvar = trim("`xvar'")

    * Count control and absorb variables
    local num_controls = 0
    if "`controls'" != "" {
        local num_controls : word count `controls'
    }
    local num_absorb = 0
    if "`absorb'" != "" {
        local num_absorb : word count `absorb'
    }

    * Minimal sample marking - just handle if/in, let plugin handle missing values
    * This avoids expensive markout operations on large datasets
    marksample touse, novarlist

    * Get observation count from Stata (fast, uses existing if/in)
    qui count if `touse'
    local nobs = r(N)
    if `nobs' == 0 {
        di as error "cbinscatter: no observations"
        exit 2000
    }

    * Count by-groups if specified (use levelsof instead of tab - faster)
    local num_by_groups = 1
    if "`by'" != "" {
        * Convert by variable to numeric if string
        capture confirm numeric variable `by'
        if _rc != 0 {
            tempvar by_numeric
            qui egen `by_numeric' = group(`by') if `touse'
            local by_orig "`by'"
            local by "`by_numeric'"
        }
        qui levelsof `by' if `touse', local(bylevels)
        local num_by_groups : word count `bylevels'
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
    scalar __cbinscatter_method = `method_num'
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
        capture scalar drop __cbinscatter_reportreg __cbinscatter_method
        capture scalar drop __cbinscatter_maxiter __cbinscatter_tolerance
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
        local plugin_overhead = r(t99) - __cbinscatter_time_total
        if `plugin_overhead' < 0 local plugin_overhead = 0
        di as text ""
        di as text "{hline 55}"
        di as text "cbinscatter timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f __cbinscatter_time_load " sec"
        di as text "    Residualize:            " as result %8.4f __cbinscatter_time_resid " sec"
        di as text "    Bin computation:        " as result %8.4f __cbinscatter_time_bins " sec"
        di as text "    Line fitting:           " as result %8.4f __cbinscatter_time_fit " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f __cbinscatter_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Plugin call overhead:     " as result %8.4f `plugin_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f r(t99) " sec"
        di as text "{hline 55}"
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
        * Default colors (same as binscatter)
        local colorlist "navy maroon forest_green dkorange teal cranberry lavender khaki sienna emidblue emerald brown erose gold bluishgray"

        * Parse user color options or use defaults
        if "`colors'" == "" {
            local colors "`colorlist'"
        }

        * Get bin data from matrix
        tempname binmat
        matrix `binmat' = __cbinscatter_bins

        * Find x range for fit lines
        local x_min = .
        local x_max = .
        local nrows = rowsof(`binmat')
        forval r = 1/`nrows' {
            local xval = `binmat'[`r', 3]
            if `xval' != . {
                if `x_min' == . | `xval' < `x_min' {
                    local x_min = `xval'
                }
                if `x_max' == . | `xval' > `x_max' {
                    local x_max = `xval'
                }
            }
        }

        * Build scatter commands using scatteri (like binscatter)
        local scatters ""
        local fits ""
        local legend_labels ""
        local legend_order ""

        if `actual_num_groups' == 1 {
            * Single group - build scatteri coordinate list
            local coords ""
            forval r = 1/`nrows' {
                local yval = `binmat'[`r', 4]
                local xval = `binmat'[`r', 3]
                if `yval' != . & `xval' != . {
                    local coords "`coords' `yval' `xval'"
                }
            }

            * Get color for this series
            local thiscolor : word 1 of `colors'

            local scatters "(scatteri `coords', mcolor(`thiscolor') msymbol(O))"

            * Add fit line if requested
            if `linetype_num' > 0 {
                local coef_const = __cbinscatter_coefs[1, 1]
                local coef_linear = __cbinscatter_coefs[1, 2]

                if `linetype_num' == 1 {
                    * Linear fit using function
                    local fits "(function `coef_linear'*x + `coef_const', range(`x_min' `x_max') lcolor(`thiscolor'))"
                }
                else if `linetype_num' >= 2 {
                    local coef_quad = __cbinscatter_coefs[1, 3]
                    if `linetype_num' == 2 {
                        local fits "(function `coef_quad'*x^2 + `coef_linear'*x + `coef_const', range(`x_min' `x_max') lcolor(`thiscolor'))"
                    }
                    else {
                        local coef_cubic = __cbinscatter_coefs[1, 4]
                        local fits "(function `coef_cubic'*x^3 + `coef_quad'*x^2 + `coef_linear'*x + `coef_const', range(`x_min' `x_max') lcolor(`thiscolor'))"
                    }
                }
            }

            * Legend off for single group
            local legend_opt "legend(off)"
        }
        else {
            * Multiple by-groups
            local legend_idx = 1

            forval g = 1/`actual_num_groups' {
                * Build coordinate list for this group
                local coords ""
                local grp_x_min = .
                local grp_x_max = .

                forval r = 1/`nrows' {
                    local grpval = `binmat'[`r', 1]
                    if `grpval' == `g' {
                        local yval = `binmat'[`r', 4]
                        local xval = `binmat'[`r', 3]
                        if `yval' != . & `xval' != . {
                            local coords "`coords' `yval' `xval'"
                            if `grp_x_min' == . | `xval' < `grp_x_min' {
                                local grp_x_min = `xval'
                            }
                            if `grp_x_max' == . | `xval' > `grp_x_max' {
                                local grp_x_max = `xval'
                            }
                        }
                    }
                }

                * Get color for this series
                local thiscolor : word `g' of `colors'
                if "`thiscolor'" == "" {
                    local thiscolor "navy"
                }

                local scatters "`scatters' (scatteri `coords', mcolor(`thiscolor') msymbol(O))"
                local legend_labels `"`legend_labels' label(`legend_idx' "`by'==`g'")"'
                local legend_order "`legend_order' `legend_idx'"
                local legend_idx = `legend_idx' + 1

                * Add fit line for this group if requested
                if `linetype_num' > 0 & "`coords'" != "" {
                    local coef_const = __cbinscatter_coefs[`g', 1]
                    local coef_linear = __cbinscatter_coefs[`g', 2]

                    if `linetype_num' == 1 {
                        local fits "`fits' (function `coef_linear'*x + `coef_const', range(`grp_x_min' `grp_x_max') lcolor(`thiscolor'))"
                    }
                    else if `linetype_num' >= 2 {
                        local coef_quad = __cbinscatter_coefs[`g', 3]
                        if `linetype_num' == 2 {
                            local fits "`fits' (function `coef_quad'*x^2 + `coef_linear'*x + `coef_const', range(`grp_x_min' `grp_x_max') lcolor(`thiscolor'))"
                        }
                        else {
                            local coef_cubic = __cbinscatter_coefs[`g', 4]
                            local fits "`fits' (function `coef_cubic'*x^3 + `coef_quad'*x^2 + `coef_linear'*x + `coef_const', range(`grp_x_min' `grp_x_max') lcolor(`thiscolor'))"
                        }
                    }
                }
            }

            * Build legend option
            if `"`legend'"' == "" {
                local legend_opt `"legend(`legend_labels' order(`legend_order'))"'
            }
            else {
                local legend_opt `"legend(`legend')"'
            }
        }

        * Set default titles if not specified (like binscatter)
        if `"`ytitle'"' == "" {
            local ytitle "`depvar'"
        }
        if `"`xtitle'"' == "" {
            local xtitle "`xvar'"
        }

        * Build title option only if specified
        local title_opt ""
        if `"`title'"' != "" {
            local title_opt `"title(`"`title'"')"'
        }

        * Execute graph (matching binscatter format)
        * Use temporary frame to avoid twoway overhead on large datasets
        * (twoway speed depends on dataset size in memory, even for scatteri)
        local _use_frame = 0
        if c(stata_version) >= 16 {
            local _use_frame = 1
        }

        if `_use_frame' {
            * Stata 16+: use frames for speed (no preserve/restore overhead)
            * Use capture to ensure proper cleanup if twoway fails
            local _orig_frame = c(frame)
            capture frame drop _ctools_graph_frame
            frame create _ctools_graph_frame
            frame change _ctools_graph_frame

            capture noisily twoway `scatters' `fits', ///
                graphregion(fcolor(white)) ///
                ytitle(`"`ytitle'"') ///
                xtitle(`"`xtitle'"') ///
                `legend_opt' ///
                `title_opt' ///
                `options'

            local _graph_rc = _rc

            * Always clean up frame, even on error
            frame change `_orig_frame'
            capture frame drop _ctools_graph_frame

            if `_graph_rc' {
                di as error "cbinscatter: graph generation failed (rc=`_graph_rc')"
                exit `_graph_rc'
            }
        }
        else {
            * Stata 14-15: fall back to preserve/keep/restore
            * Use capture to ensure restore happens even on error
            preserve
            qui keep in 1/1

            capture noisily twoway `scatters' `fits', ///
                graphregion(fcolor(white)) ///
                ytitle(`"`ytitle'"') ///
                xtitle(`"`xtitle'"') ///
                `legend_opt' ///
                `title_opt' ///
                `options'

            local _graph_rc = _rc
            restore

            if `_graph_rc' {
                di as error "cbinscatter: graph generation failed (rc=`_graph_rc')"
                exit `_graph_rc'
            }
        }
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
    capture scalar drop __cbinscatter_reportreg __cbinscatter_method
    capture scalar drop __cbinscatter_maxiter __cbinscatter_tolerance
    capture scalar drop __cbinscatter_N __cbinscatter_N_dropped __cbinscatter_num_groups
    capture scalar drop __cbinscatter_time_load __cbinscatter_time_resid
    capture scalar drop __cbinscatter_time_bins __cbinscatter_time_fit __cbinscatter_time_total

    timer clear 99
end
