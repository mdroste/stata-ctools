* benchmark_creghdfe.do
* Benchmark comparison between reghdfe and creghdfe
* Created: 2026-01-06

clear all
set more off

* Add build directory to adopath
adopath + "../build"

* ===========================================================================
* Program: compare_hdfe
*
* Runs the same regression with reghdfe and creghdfe, compares all stored
* results, and reports PASSED/FAILED based on tolerance threshold.
*
* Syntax:
*   compare_hdfe depvar indepvars, absorb(varlist) [vce(string) tol(real)]
*
* Options:
*   absorb(varlist)  - Fixed effects to absorb (required)
*   vce(string)      - VCE type (optional)
*   tol(real)        - Tolerance for comparison (default: 1e-7)
* ===========================================================================

capture program drop compare_hdfe
program define compare_hdfe, rclass
    version 14.0

    syntax varlist(min=2 fv), Absorb(varlist) [VCE(string) TOLerance(real 1e-7)]

    local tol = `tolerance'

    * Parse depvar and indepvars
    gettoken depvar indepvars : varlist

    di as text ""
    di as text "{hline 70}"
    di as text "BENCHMARK: reghdfe vs creghdfe"
    di as text "{hline 70}"
    di as text "Specification: `depvar' `indepvars', absorb(`absorb')" _c
    if "`vce'" != "" {
        di as text " vce(`vce')"
    }
    else {
        di as text ""
    }
    di as text "Tolerance: `tol'"
    di as text "{hline 70}"

    * =========================================================================
    * Run reghdfe
    * =========================================================================
    di as text _n "Running reghdfe..."
    timer clear 1
    timer on 1

    if "`vce'" != "" {
        quietly reghdfe `varlist', absorb(`absorb') vce(`vce')
    }
    else {
        quietly reghdfe `varlist', absorb(`absorb')
    }

    timer off 1
    quietly timer list 1
    local time_reghdfe = r(t1)

    * Store reghdfe results
    local reghdfe_N = e(N)
    local reghdfe_df_m = e(df_m)
    local reghdfe_df_r = e(df_r)
    local reghdfe_F = e(F)
    local reghdfe_r2 = e(r2)
    local reghdfe_r2_a = e(r2_a)
    local reghdfe_r2_within = e(r2_within)
    local reghdfe_rmse = e(rmse)
    local reghdfe_rss = e(rss)
    local reghdfe_tss = e(tss)
    local reghdfe_df_a = e(df_a)
    local reghdfe_num_singletons = e(num_singletons)

    tempname reghdfe_b reghdfe_V
    matrix `reghdfe_b' = e(b)
    matrix `reghdfe_V' = e(V)

    local nvars = colsof(`reghdfe_b')

    di as text "  Time: " as result %9.4f `time_reghdfe' as text " seconds"

    * =========================================================================
    * Run creghdfe
    * =========================================================================
    di as text _n "Running creghdfe..."
    timer clear 2
    timer on 2

    if "`vce'" != "" {
        quietly creghdfe `varlist', absorb(`absorb') vce(`vce')
    }
    else {
        quietly creghdfe `varlist', absorb(`absorb')
    }

    timer off 2
    quietly timer list 2
    local time_creghdfe = r(t2)

    * Store creghdfe results
    local creghdfe_N = e(N)
    local creghdfe_df_m = e(df_m)
    local creghdfe_df_r = e(df_r)
    local creghdfe_F = e(F)
    local creghdfe_r2 = e(r2)
    local creghdfe_r2_a = e(r2_a)
    local creghdfe_r2_within = e(r2_within)
    local creghdfe_rmse = e(rmse)
    local creghdfe_rss = e(rss)
    local creghdfe_tss = e(tss)
    local creghdfe_df_a = e(df_a)
    local creghdfe_num_singletons = e(num_singletons)

    tempname creghdfe_b creghdfe_V
    matrix `creghdfe_b' = e(b)
    matrix `creghdfe_V' = e(V)

    di as text "  Time: " as result %9.4f `time_creghdfe' as text " seconds"

    * =========================================================================
    * Compare results
    * =========================================================================
    di as text _n "{hline 70}"
    di as text "COMPARISON RESULTS"
    di as text "{hline 70}"
    di as text ""
    di as text %25s "Statistic" " | " %15s "reghdfe" " | " %15s "creghdfe" " | " %12s "Diff" " | " "Status"
    di as text "{hline 25}-+-{hline 15}-+-{hline 15}-+-{hline 12}-+-{hline 8}"

    local all_passed = 1
    local num_comparisons = 0
    local num_passed = 0

    * Compare scalar statistics
    foreach stat in N df_m df_r F r2 r2_a r2_within rmse rss tss df_a num_singletons {
        local val_reghdfe = `reghdfe_`stat''
        local val_creghdfe = `creghdfe_`stat''
        local diff = abs(`val_reghdfe' - `val_creghdfe')

        * Handle relative difference for large values
        if abs(`val_reghdfe') > 1 {
            local rel_diff = `diff' / abs(`val_reghdfe')
            local passed = (`rel_diff' < `tol')
        }
        else {
            local passed = (`diff' < `tol')
        }

        local ++num_comparisons
        if `passed' {
            local ++num_passed
            local status_text "PASS"
            local status_color "as result"
        }
        else {
            local all_passed = 0
            local status_text "FAIL"
            local status_color "as error"
        }

        di as text %25s "`stat'" " | " as result %15.6g `val_reghdfe' as text " | " as result %15.6g `val_creghdfe' as text " | " as result %12.2e `diff' " | " `status_color' "`status_text'"
    }

    * Compare coefficients
    di as text "{hline 25}-+-{hline 15}-+-{hline 15}-+-{hline 12}-+-{hline 8}"
    di as text %25s "COEFFICIENTS" " | " " " " | " " " " | " " " " | "

    local coef_names : colnames `reghdfe_b'
    local k = 1
    foreach coef of local coef_names {
        local val_reghdfe = `reghdfe_b'[1, `k']
        local val_creghdfe = `creghdfe_b'[1, `k']
        local diff = abs(`val_reghdfe' - `val_creghdfe')

        if abs(`val_reghdfe') > 1 {
            local rel_diff = `diff' / abs(`val_reghdfe')
            local passed = (`rel_diff' < `tol')
        }
        else {
            local passed = (`diff' < `tol')
        }

        local ++num_comparisons
        if `passed' {
            local ++num_passed
            local status_text "PASS"
            local status_color "as result"
        }
        else {
            local all_passed = 0
            local status_text "FAIL"
            local status_color "as error"
        }

        di as text %25s "b[`coef']" " | " as result %15.6g `val_reghdfe' as text " | " as result %15.6g `val_creghdfe' as text " | " as result %12.2e `diff' " | " `status_color' "`status_text'"
        local ++k
    }

    * Compare standard errors (diagonal of V)
    di as text "{hline 25}-+-{hline 15}-+-{hline 15}-+-{hline 12}-+-{hline 8}"
    di as text %25s "STANDARD ERRORS" " | " " " " | " " " " | " " " " | "

    local k = 1
    foreach coef of local coef_names {
        local val_reghdfe = sqrt(`reghdfe_V'[`k', `k'])
        local val_creghdfe = sqrt(`creghdfe_V'[`k', `k'])
        local diff = abs(`val_reghdfe' - `val_creghdfe')

        if abs(`val_reghdfe') > 1 {
            local rel_diff = `diff' / abs(`val_reghdfe')
            local passed = (`rel_diff' < `tol')
        }
        else {
            local passed = (`diff' < `tol')
        }

        local ++num_comparisons
        if `passed' {
            local ++num_passed
            local status_text "PASS"
            local status_color "as result"
        }
        else {
            local all_passed = 0
            local status_text "FAIL"
            local status_color "as error"
        }

        di as text %25s "se[`coef']" " | " as result %15.6g `val_reghdfe' as text " | " as result %15.6g `val_creghdfe' as text " | " as result %12.2e `diff' " | " `status_color' "`status_text'"
        local ++k
    }

    * Compare t-statistics
    di as text "{hline 25}-+-{hline 15}-+-{hline 15}-+-{hline 12}-+-{hline 8}"
    di as text %25s "T-STATISTICS" " | " " " " | " " " " | " " " " | "

    local k = 1
    foreach coef of local coef_names {
        local b_reghdfe = `reghdfe_b'[1, `k']
        local se_reghdfe = sqrt(`reghdfe_V'[`k', `k'])
        local val_reghdfe = `b_reghdfe' / `se_reghdfe'

        local b_creghdfe = `creghdfe_b'[1, `k']
        local se_creghdfe = sqrt(`creghdfe_V'[`k', `k'])
        local val_creghdfe = `b_creghdfe' / `se_creghdfe'

        local diff = abs(`val_reghdfe' - `val_creghdfe')

        if abs(`val_reghdfe') > 1 {
            local rel_diff = `diff' / abs(`val_reghdfe')
            local passed = (`rel_diff' < `tol')
        }
        else {
            local passed = (`diff' < `tol')
        }

        local ++num_comparisons
        if `passed' {
            local ++num_passed
            local status_text "PASS"
            local status_color "as result"
        }
        else {
            local all_passed = 0
            local status_text "FAIL"
            local status_color "as error"
        }

        di as text %25s "t[`coef']" " | " as result %15.6g `val_reghdfe' as text " | " as result %15.6g `val_creghdfe' as text " | " as result %12.2e `diff' " | " `status_color' "`status_text'"
        local ++k
    }

    * =========================================================================
    * Summary
    * =========================================================================
    di as text ""
    di as text "{hline 70}"
    di as text "TIMING COMPARISON"
    di as text "{hline 70}"
    di as text "  reghdfe:  " as result %9.4f `time_reghdfe' as text " seconds"
    di as text "  creghdfe: " as result %9.4f `time_creghdfe' as text " seconds"
    local time_diff = `time_reghdfe' - `time_creghdfe'
    if `time_diff' > 0 {
        local speedup = `time_reghdfe' / `time_creghdfe'
        di as text "  Difference: " as result %9.4f `time_diff' as text " seconds (creghdfe is " as result %5.2f `speedup' as text "x faster)"
    }
    else {
        local speedup = `time_creghdfe' / `time_reghdfe'
        di as text "  Difference: " as result %9.4f `time_diff' as text " seconds (reghdfe is " as result %5.2f `speedup' as text "x faster)"
    }

    di as text ""
    di as text "{hline 70}"
    di as text "OVERALL RESULT"
    di as text "{hline 70}"
    di as text "  Comparisons: " as result `num_passed' as text "/" as result `num_comparisons' as text " passed"

    if `all_passed' {
        di as text ""
        di as result "  *** PASSED ***"
        di as text ""
    }
    else {
        di as text ""
        di as error "  *** FAILED ***"
        di as text ""
    }

    * Return results
    return scalar passed = `all_passed'
    return scalar num_comparisons = `num_comparisons'
    return scalar num_passed = `num_passed'
    return scalar time_reghdfe = `time_reghdfe'
    return scalar time_creghdfe = `time_creghdfe'
    return scalar speedup = `time_reghdfe' / `time_creghdfe'

end

* ===========================================================================
* Run benchmarks
* ===========================================================================

sysuse auto, clear

di as text _n _n
di as text "╔══════════════════════════════════════════════════════════════════════════════╗"
di as text "║                         CREGHDFE BENCHMARK SUITE                             ║"
di as text "╚══════════════════════════════════════════════════════════════════════════════╝"

* Test 1: Basic regression with single FE
compare_hdfe price mpg, absorb(foreign)

* Test 2: Multiple FEs with singletons
compare_hdfe price mpg, absorb(foreign trunk)

* Test 3: Multiple FEs with singletons, robust SEs
compare_hdfe price mpg, absorb(foreign trunk) vce(robust)

* Test 4: Multiple FEs with singletons, robust SEs
compare_hdfe price mpg, absorb(foreign trunk) vce(cluster foreign)

* Test 5: Cluster by string variable (make)
compare_hdfe price mpg, absorb(foreign trunk) vce(cluster make)

di as text _n "Benchmark suite complete."
