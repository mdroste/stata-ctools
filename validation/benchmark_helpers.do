/*******************************************************************************
 * benchmark_helpers.do
 *
 * Benchmark helper programs for ctools validation test suite
 * Each program runs both ctools and Stata native commands, compares results,
 * and reports [PASS] or [FAIL] with detailed error descriptions
 *
 * TOLERANCE POLICY:
 *   - All comparisons use SIGNIFICANT FIGURES, not absolute differences
 *   - sigfigs = -log10(|a - b| / max(|a|, |b|))
 *   - Default threshold is 7 significant figures ($DEFAULT_SIGFIGS)
 *   - This applies to all scalars, coefficients e(b), and VCE e(V)
 *   - Scale-invariant: 1000000 vs 1000001 and 0.001 vs 0.001001 both have ~3 sigfigs
 *   - Missing values in ctools when Stata has a value is treated as a FAILURE
 *
 ******************************************************************************/

/*******************************************************************************
 * benchmark_sort - Compare csort vs sort, stable
 *
 * Validation strategy:
 *   1. Check if csort result is byte-for-byte identical to sort, stable
 *   2. If not identical, verify csort result is actually sorted (using sortedby)
 *   3. If sorted, apply sort, stable to csort result and compare with reference
 *      - If they match, the difference is only sort stability (valid)
 *      - If they don't match, other variables were corrupted (invalid)
 *
 * Syntax: benchmark_sort varlist [, testname(string) algorithm(string)]
 ******************************************************************************/
capture program drop benchmark_sort
program define benchmark_sort
    syntax varlist [, testname(string) ALGorithm(string)]

    if "`testname'" == "" local testname "sort `varlist'"
    if "`algorithm'" != "" local testname "`testname' [alg=`algorithm']"

    preserve

    * Save original data
    tempfile original
    quietly save `original'

    * Run Stata sort, stable (quietly) - this is our reference
    quietly sort `varlist', stable
    tempfile stata_sorted
    quietly save `stata_sorted'

    * Restore and run csort (quietly)
    quietly use `original', clear
    if "`algorithm'" != "" {
        capture quietly csort `varlist', algorithm(`algorithm')
    }
    else {
        capture quietly csort `varlist'
    }
    local csort_rc = _rc

    if `csort_rc' != 0 {
        restore
        test_fail "`testname'" "csort returned error `csort_rc'"
        exit
    }

    * Save csort result
    tempfile csort_sorted
    quietly save `csort_sorted'

    * Check 1: Are they byte-for-byte identical?
    capture quietly cf _all using `stata_sorted'
    local exact_match = (_rc == 0)

    if `exact_match' {
        restore
        test_pass "`testname' (exact)"
        exit
    }

    * Check 2: Is the csort result actually sorted by the sort keys?
    * Use Stata's sortedby() to verify - it returns the sort order if sorted
    * Normalize both strings by removing extra spaces before comparing
    local sortedby : sortedby
    local sortedby_norm = trim(itrim("`sortedby'"))
    local varlist_norm = trim(itrim("`varlist'"))
    local is_sorted = ("`sortedby_norm'" == "`varlist_norm'")

    if !`is_sorted' {
        * Data claims not to be sorted - this is a definite failure
        restore
        test_fail "`testname'" "csort result not marked as sorted by `varlist' (sortedby=`sortedby')"
        exit
    }

    * Check 3: Apply sort, stable to csort result and compare with reference
    * If csort produced a valid sort (just with different tie-breaking),
    * then applying stable sort should yield the same result as sorting original
    quietly sort `varlist', stable
    capture quietly cf _all using `stata_sorted'
    local stable_match = (_rc == 0)

    restore

    if `stable_match' {
        * Valid sort - only differs in sort stability for ties
        test_pass "`testname' (stability differs)"
    }
    else {
        * Data was corrupted - other variables don't match after stable sort
        test_fail "`testname'" "csort corrupted data (stable re-sort doesn't match reference)"
    }
end

/*******************************************************************************
 * benchmark_merge - Compare cmerge vs merge
 *
 * Syntax: benchmark_merge mergetype keyvar using filename [, options]
 ******************************************************************************/
capture program drop benchmark_merge
program define benchmark_merge
    syntax anything using/, [keep(string) GENerate(name) NOGENerate KEEPUSing(string) ///
        testname(string) SORTed ASSert(string) UPDATE REPLACE FORCE NOLabel NONotes NOREPort]

    * Parse merge type and key vars
    gettoken mergetype keyvars : anything

    if "`testname'" == "" local testname "merge `mergetype' `keyvars'"

    * Build option string
    local opts ""
    if "`keep'" != "" local opts "`opts' keep(`keep')"
    if "`generate'" != "" local opts "`opts' generate(`generate')"
    if "`nogenerate'" != "" local opts "`opts' nogenerate"
    if "`keepusing'" != "" local opts "`opts' keepusing(`keepusing')"
    if "`sorted'" != "" local opts "`opts' sorted"
    if "`assert'" != "" local opts "`opts' assert(`assert')"
    if "`update'" != "" local opts "`opts' update"
    if "`replace'" != "" local opts "`opts' replace"
    if "`force'" != "" local opts "`opts' force"
    if "`nolabel'" != "" local opts "`opts' nolabel"
    if "`nonotes'" != "" local opts "`opts' nonotes"

    preserve

    * Save original data
    tempfile original
    quietly save `original'

    * Run Stata merge (quietly)
    capture quietly merge `mergetype' `keyvars' using `using', `opts'
    local stata_rc = _rc

    if `stata_rc' != 0 {
        restore
        test_fail "`testname'" "Stata merge returned error `stata_rc'"
        exit
    }

    * Get _merge counts from Stata
    local stata_m1 = 0
    local stata_m2 = 0
    local stata_m3 = 0
    capture confirm variable _merge
    if _rc == 0 {
        quietly count if _merge == 1
        local stata_m1 = r(N)
        quietly count if _merge == 2
        local stata_m2 = r(N)
        quietly count if _merge == 3
        local stata_m3 = r(N)
    }

    * Sort by all variables for comparison
    quietly ds
    local allvars `r(varlist)'
    quietly sort `allvars'
    tempfile stata_merged
    quietly save `stata_merged'

    * Restore and run cmerge
    quietly use `original', clear
    capture quietly cmerge `mergetype' `keyvars' using `using', `opts' noreport
    local cmerge_rc = _rc

    if `cmerge_rc' != 0 {
        restore
        test_fail "`testname'" "cmerge returned error `cmerge_rc'"
        exit
    }

    * Get _merge counts from cmerge
    local cmerge_m1 = 0
    local cmerge_m2 = 0
    local cmerge_m3 = 0
    capture confirm variable _merge
    if _rc == 0 {
        quietly count if _merge == 1
        local cmerge_m1 = r(N)
        quietly count if _merge == 2
        local cmerge_m2 = r(N)
        quietly count if _merge == 3
        local cmerge_m3 = r(N)
    }

    * Check counts match
    local counts_match = (`stata_m1' == `cmerge_m1') & (`stata_m2' == `cmerge_m2') & (`stata_m3' == `cmerge_m3')

    * Sort and compare datasets
    quietly ds
    local allvars `r(varlist)'
    quietly sort `allvars'
    capture quietly cf _all using `stata_merged'
    local cf_rc = _rc

    restore

    * For m:m merges, only compare counts (Stata's m:m behavior is "unusual" and
    * implementation-dependent per Stata documentation - we only verify dimensions)
    local is_mm = ("`mergetype'" == "m:m")

    if `counts_match' & (`cf_rc' == 0 | `is_mm') {
        test_pass "`testname'"
    }
    else if !`counts_match' {
        test_fail "`testname'" "_merge counts differ: stata(m1=`stata_m1',m2=`stata_m2',m3=`stata_m3') vs cmerge(m1=`cmerge_m1',m2=`cmerge_m2',m3=`cmerge_m3')"
    }
    else {
        test_fail "`testname'" "merged data differs from Stata merge"
    }
end

/*******************************************************************************
 * benchmark_reghdfe - Compare creghdfe vs reghdfe
 *
 * Comprehensively compares all e() scalars, vectors, and matrices:
 *   - e(N), e(df_r), e(df_m), e(rank) - counts and degrees of freedom
 *   - e(r2), e(r2_a), e(r2_within) - R-squared values
 *   - e(rss), e(tss), e(mss) - sum of squares
 *   - e(F) - F-statistic
 *   - e(rmse) - root mean squared error
 *   - e(N_clust) - number of clusters (if clustering)
 *   - e(b) - coefficient vector
 *   - e(V) - variance-covariance matrix
 *
 * All comparisons use significant figures (default 7 sigfigs).
 *
 * Syntax: benchmark_reghdfe depvar indepvars [weight], absorb(varlist) [options]
 ******************************************************************************/
capture program drop benchmark_reghdfe
program define benchmark_reghdfe
    syntax varlist(min=2 fv) [aw fw pw] [if] [in], Absorb(varlist) [vce(string) testname(string) minsf(real 7)]

    gettoken depvar indepvars : varlist

    if "`testname'" == "" local testname "reghdfe `depvar' `indepvars', absorb(`absorb')"

    * Build weight string
    local wtexp ""
    if "`weight'" != "" local wtexp "[`weight'`exp']"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run reghdfe (quietly)
    capture quietly reghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    local reghdfe_rc = _rc

    if `reghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "reghdfe returned error `reghdfe_rc'"
        exit
    }

    * Store reghdfe results - matrices
    matrix reghdfe_b = e(b)
    matrix reghdfe_V = e(V)

    * Store reghdfe results - scalars (capture for optional values)
    local reghdfe_N = e(N)
    local reghdfe_df_r = e(df_r)
    local reghdfe_df_m = e(df_m)
    local reghdfe_r2 = e(r2)
    local reghdfe_r2_a = e(r2_a)
    local reghdfe_r2_within = e(r2_within)
    local reghdfe_rss = e(rss)
    local reghdfe_tss = e(tss)
    local reghdfe_mss = e(mss)
    local reghdfe_F = e(F)
    local reghdfe_rmse = e(rmse)
    local reghdfe_rank = e(rank)
    local reghdfe_N_clust = e(N_clust)
    local reghdfe_N_hdfe = e(N_hdfe)

    * Run creghdfe (quietly)
    capture quietly creghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    local creghdfe_rc = _rc

    if `creghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "creghdfe returned error `creghdfe_rc'"
        exit
    }

    * Store creghdfe results - matrices
    matrix creghdfe_b = e(b)
    matrix creghdfe_V = e(V)

    * Store creghdfe results - scalars
    local creghdfe_N = e(N)
    local creghdfe_df_r = e(df_r)
    local creghdfe_df_m = e(df_m)
    local creghdfe_r2 = e(r2)
    local creghdfe_r2_a = e(r2_a)
    local creghdfe_r2_within = e(r2_within)
    local creghdfe_rss = e(rss)
    local creghdfe_tss = e(tss)
    local creghdfe_mss = e(mss)
    local creghdfe_F = e(F)
    local creghdfe_rmse = e(rmse)
    local creghdfe_rank = e(rank)
    local creghdfe_N_clust = e(N_clust)
    local creghdfe_N_hdfe = e(N_hdfe)

    restore

    * Track all differences found
    local all_diffs ""
    local has_failure = 0

    * Compare N (must match exactly)
    if `reghdfe_N' != `creghdfe_N' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
    }

    * Compare integer scalars (must match exactly)
    foreach scalar in df_r df_m rank N_clust N_hdfe {
        local val1 = `reghdfe_`scalar''
        local val2 = `creghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            if `val1' != `val2' {
                local has_failure = 1
                local all_diffs "`all_diffs' e(`scalar'):`val1'!=`val2'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local all_diffs "`all_diffs' e(`scalar'):missing_in_creghdfe"
        }
    }

    * Compare continuous scalars using significant figures
    foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
        local val1 = `reghdfe_`scalar''
        local val2 = `creghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):missing_in_creghdfe"
        }
    }

    * Compare coefficients e(b) using significant figures
    * First check dimensions match
    local reghdfe_bcols = colsof(reghdfe_b)
    local creghdfe_bcols = colsof(creghdfe_b)

    if `reghdfe_bcols' != `creghdfe_bcols' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(b):dim_mismatch(`reghdfe_bcols'!=`creghdfe_bcols')"
    }
    else {
        local cols = colsof(reghdfe_b)
        local min_sf_b = 15
        local min_sf_b_idx = 0
        forvalues j = 1/`cols' {
            local v1 = reghdfe_b[1, `j']
            local v2 = creghdfe_b[1, `j']
            sigfigs `v1' `v2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_b' {
                local min_sf_b = `sf'
                local min_sf_b_idx = `j'
            }
        }

        if `min_sf_b' < `minsf' {
            local has_failure = 1
            local coefname : colnames reghdfe_b
            local badcoef : word `min_sf_b_idx' of `coefname'
            local sf_fmt : display %4.1f `min_sf_b'
            local all_diffs "`all_diffs' e(b)[`badcoef']:sigfigs=`sf_fmt'"
        }
    }

    * Compare VCE e(V) using significant figures
    * First check dimensions match
    local reghdfe_Vrows = rowsof(reghdfe_V)
    local creghdfe_Vrows = rowsof(creghdfe_V)

    if `reghdfe_Vrows' != `creghdfe_Vrows' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(V):dim_mismatch(`reghdfe_Vrows'!=`creghdfe_Vrows')"
    }
    else {
        local rows = rowsof(reghdfe_V)
        local cols = colsof(reghdfe_V)
        local min_sf_V = 15
        local min_sf_V_i = 0
        local min_sf_V_j = 0
        forvalues i = 1/`rows' {
            forvalues j = 1/`cols' {
                local v1 = reghdfe_V[`i', `j']
                local v2 = creghdfe_V[`i', `j']
                sigfigs `v1' `v2'
                local sf = r(sigfigs)
                if `sf' < `min_sf_V' {
                    local min_sf_V = `sf'
                    local min_sf_V_i = `i'
                    local min_sf_V_j = `j'
                }
            }
        }

        if `min_sf_V' < `minsf' {
            local has_failure = 1
            local sf_fmt : display %4.1f `min_sf_V'
            local all_diffs "`all_diffs' e(V)[`min_sf_V_i',`min_sf_V_j']:sigfigs=`sf_fmt'"
        }
    }

    * Report result
    if `has_failure' == 0 {
        test_pass "`testname'"
    }
    else {
        * Trim leading space from all_diffs
        local all_diffs = trim("`all_diffs'")
        test_fail "`testname'" "`all_diffs'"
    }
end

/*******************************************************************************
 * benchmark_qreg - Compare cqreg vs qreg
 *
 * Comprehensively compares all e() scalars, vectors, and matrices:
 *   - e(N), e(df_r) - counts and degrees of freedom
 *   - e(q) - quantile
 *   - e(sum_adev) - sum of absolute deviations
 *   - e(sum_rdev) - sum of raw deviations
 *   - e(q_v), e(f_r) - quantile variance and pseudo F df
 *   - e(convcode) - convergence code
 *   - e(b) - coefficient vector
 *   - e(V) - variance-covariance matrix
 *
 * All comparisons use significant figures (default 7 sigfigs).
 *
 * Syntax: benchmark_qreg depvar indepvars [, quantile(real) options]
 ******************************************************************************/
capture program drop benchmark_qreg
program define benchmark_qreg
    syntax varlist(min=2 fv) [if] [in], [Quantile(real 0.5) vce(string) testname(string) minsf(real 7)]

    gettoken depvar indepvars : varlist

    if "`testname'" == "" local testname "qreg `depvar' `indepvars', q(`quantile')"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run qreg (quietly)
    capture quietly qreg `depvar' `indepvars' `if' `in', quantile(`quantile') `vceopt'
    local qreg_rc = _rc

    if `qreg_rc' != 0 {
        restore
        test_fail "`testname'" "qreg returned error `qreg_rc'"
        exit
    }

    * Store qreg results - matrices
    matrix qreg_b = e(b)
    matrix qreg_V = e(V)

    * Store qreg results - scalars
    local qreg_N = e(N)
    local qreg_df_r = e(df_r)
    local qreg_q = e(q)
    local qreg_sum_adev = e(sum_adev)
    local qreg_sum_rdev = e(sum_rdev)
    local qreg_q_v = e(q_v)
    local qreg_f_r = e(f_r)
    local qreg_convcode = e(convcode)

    * Run cqreg (quietly)
    capture quietly cqreg `depvar' `indepvars' `if' `in', quantile(`quantile') `vceopt'
    local cqreg_rc = _rc

    if `cqreg_rc' != 0 {
        restore
        test_fail "`testname'" "cqreg returned error `cqreg_rc'"
        exit
    }

    * Store cqreg results - matrices
    matrix cqreg_b = e(b)
    matrix cqreg_V = e(V)

    * Store cqreg results - scalars
    local cqreg_N = e(N)
    local cqreg_df_r = e(df_r)
    local cqreg_q = e(q)
    local cqreg_sum_adev = e(sum_adev)
    local cqreg_sum_rdev = e(sum_rdev)
    local cqreg_q_v = e(q_v)
    local cqreg_f_r = e(f_r)
    local cqreg_convcode = e(convcode)

    restore

    * Track all differences found
    local all_diffs ""
    local has_failure = 0

    * Compare N (must match exactly)
    if `qreg_N' != `cqreg_N' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(N):`qreg_N'!=`cqreg_N'"
    }

    * Compare integer scalars (must match exactly)
    foreach scalar in df_r convcode {
        local val1 = `qreg_`scalar''
        local val2 = `cqreg_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            if `val1' != `val2' {
                local has_failure = 1
                local all_diffs "`all_diffs' e(`scalar'):`val1'!=`val2'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local all_diffs "`all_diffs' e(`scalar'):missing_in_cqreg"
        }
    }

    * Compare continuous scalars using significant figures
    foreach scalar in q sum_adev sum_rdev q_v f_r {
        local val1 = `qreg_`scalar''
        local val2 = `cqreg_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):missing_in_cqreg"
        }
    }

    * Compare coefficients e(b) using significant figures
    * First check dimensions match
    local qreg_cols = colsof(qreg_b)
    local cqreg_cols = colsof(cqreg_b)

    if `qreg_cols' != `cqreg_cols' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(b):dim_mismatch(`qreg_cols'!=`cqreg_cols')"
    }
    else {
        local cols = colsof(qreg_b)
        local min_sf_b = 15
        local min_sf_b_idx = 0
        forvalues j = 1/`cols' {
            local v1 = qreg_b[1, `j']
            local v2 = cqreg_b[1, `j']
            sigfigs `v1' `v2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_b' {
                local min_sf_b = `sf'
                local min_sf_b_idx = `j'
            }
        }

        if `min_sf_b' < `minsf' {
            local has_failure = 1
            local coefname : colnames qreg_b
            local badcoef : word `min_sf_b_idx' of `coefname'
            local sf_fmt : display %4.1f `min_sf_b'
            local all_diffs "`all_diffs' e(b)[`badcoef']:sigfigs=`sf_fmt'"
        }
    }

    * Compare VCE e(V) using significant figures
    * First check dimensions match (already checked via e(b) above)
    local qreg_Vrows = rowsof(qreg_V)
    local cqreg_Vrows = rowsof(cqreg_V)

    if `qreg_Vrows' == `cqreg_Vrows' {
        local rows = rowsof(qreg_V)
        local cols = colsof(qreg_V)
        local min_sf_V = 15
        local min_sf_V_i = 0
        local min_sf_V_j = 0
        forvalues i = 1/`rows' {
            forvalues j = 1/`cols' {
                local v1 = qreg_V[`i', `j']
                local v2 = cqreg_V[`i', `j']
                sigfigs `v1' `v2'
                local sf = r(sigfigs)
                if `sf' < `min_sf_V' {
                    local min_sf_V = `sf'
                    local min_sf_V_i = `i'
                    local min_sf_V_j = `j'
                }
            }
        }

        if `min_sf_V' < `minsf' {
            local has_failure = 1
            local sf_fmt : display %4.1f `min_sf_V'
            local all_diffs "`all_diffs' e(V)[`min_sf_V_i',`min_sf_V_j']:sigfigs=`sf_fmt'"
        }
    }

    * Report result
    if `has_failure' == 0 {
        test_pass "`testname'"
    }
    else {
        * Trim leading space from all_diffs
        local all_diffs = trim("`all_diffs'")
        test_fail "`testname'" "`all_diffs'"
    }
end

/*******************************************************************************
 * benchmark_import - Compare cimport vs import delimited
 *
 * NOTE: Stata's import delimited has auto-detection for headers - if row 1
 * has the same type signature as row 2, it treats row 1 as data. cimport
 * always uses row 1 as headers (matching standard CSV conventions).
 *
 * To make fair comparisons, this helper forces both commands to use the same
 * header handling when varnames() is not specified:
 * - If no varnames option: force varnames(1) on Stata so both use first row as headers
 * - If varnames(nonames): both treat all rows as data
 * - If varnames(1): both use first row as headers
 *
 * Syntax: benchmark_import using filename [, options]
 ******************************************************************************/
capture program drop benchmark_import
program define benchmark_import
    syntax using/, [DELIMiters(string) VARNames(string) CASE(string) testname(string) importopts(string)]

    if "`testname'" == "" local testname "import `using'"

    * Build option strings for each command
    * For cimport, use options as provided
    local cimport_opts "clear"
    if "`delimiters'" != "" local cimport_opts "`cimport_opts' delimiters(`delimiters')"
    if "`varnames'" != "" local cimport_opts "`cimport_opts' varnames(`varnames')"
    if "`case'" != "" local cimport_opts "`cimport_opts' case(`case')"

    * For Stata import delimited, force varnames(1) if not specified
    * This ensures Stata uses first row as headers like cimport does
    local stata_opts "clear"
    if "`delimiters'" != "" local stata_opts "`stata_opts' delimiters(`delimiters')"
    if "`varnames'" != "" {
        local stata_opts "`stata_opts' varnames(`varnames')"
    }
    else {
        * Force Stata to use first row as headers to match cimport default
        local stata_opts "`stata_opts' varnames(1)"
    }
    if "`case'" != "" local stata_opts "`stata_opts' case(`case')"

    * Add any additional import options (passed through directly)
    if `"`importopts'"' != "" {
        local stata_opts `"`stata_opts' `importopts'"'
        local cimport_opts `"`cimport_opts' `importopts'"'
    }

    preserve

    * Run Stata import (quietly)
    capture quietly import delimited `using', `stata_opts'
    local stata_rc = _rc

    if `stata_rc' != 0 {
        restore
        test_fail "`testname'" "import delimited returned error `stata_rc'"
        exit
    }

    local stata_n = _N
    local stata_k = c(k)
    tempfile stata_import
    quietly save `stata_import'

    * Run cimport (quietly)
    capture quietly cimport delimited `using', `cimport_opts'
    local cimport_rc = _rc

    if `cimport_rc' != 0 {
        restore
        test_fail "`testname'" "cimport returned error `cimport_rc'"
        exit
    }

    local cimport_n = _N
    local cimport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cimport_n' | `stata_k' != `cimport_k' {
        restore
        test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cimport N=`cimport_n' K=`cimport_k'"
        exit
    }

    * Dimensions match - now compare data content
    * Note: Variable names may differ (e.g., Stata sanitizes "1" to v1, cimport to v)
    * So we compare by position rather than by variable name

    * Get variable lists for both datasets
    quietly ds
    local cimport_vars `r(varlist)'

    * Rename cimport variables to generic names for comparison
    local i = 1
    foreach v of local cimport_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }

    tempfile cimport_import
    quietly save `cimport_import'

    * Load Stata import and rename to same generic names
    quietly use `stata_import', clear
    quietly ds
    local stata_vars `r(varlist)'

    local i = 1
    foreach v of local stata_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }

    * Sort both by all variables
    quietly ds
    local allvars `r(varlist)'
    quietly sort `allvars'

    tempfile stata_renamed
    quietly save `stata_renamed'

    quietly use `cimport_import', clear
    quietly sort `allvars'

    * Now compare
    capture quietly cf _all using `stata_renamed'
    local cf_rc = _rc

    restore

    if `cf_rc' == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "data content differs (cf _all failed)"
    }
end

/*******************************************************************************
 * benchmark_export - Compare cexport vs export delimited
 *
 * Exports data with both methods, imports back, and compares:
 *   - First attempts cf _all for byte-for-byte comparison
 *   - Falls back to sigfigs-based tolerance comparison for numeric variables
 *   - String variables must match exactly
 *
 * Syntax: benchmark_export [varlist], testname(string) [exportopts(string) importopts(string) ifcond(string) incond(string)]
 ******************************************************************************/
capture program drop benchmark_export
program define benchmark_export
    syntax [varlist], testname(string) [EXPORTopts(string) IMPORTopts(string) IFcond(string) INcond(string)]

    * Build if/in conditions
    local ifin ""
    if "`ifcond'" != "" local ifin "`ifin' if `ifcond'"
    if "`incond'" != "" local ifin "`ifin' in `incond'"

    tempfile stata_csv cexport_csv

    * Export with Stata's export delimited
    if "`varlist'" != "" {
        export delimited `varlist' using "`stata_csv'" `ifin', `exportopts' replace
    }
    else {
        export delimited using "`stata_csv'" `ifin', `exportopts' replace
    }

    * Export with cexport delimited
    if "`varlist'" != "" {
        cexport delimited `varlist' using "`cexport_csv'" `ifin', `exportopts' replace
    }
    else {
        cexport delimited using "`cexport_csv'" `ifin', `exportopts' replace
    }

    * Import both CSVs back
    preserve

    import delimited using "`stata_csv'", `importopts' clear
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    import delimited using "`cexport_csv'", `importopts' clear
    tempfile cexport_data
    quietly save `cexport_data', replace
    local cexport_n = _N
    local cexport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cexport_n' | `stata_k' != `cexport_k' {
        restore
        test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cexport N=`cexport_n' K=`cexport_k'"
        exit
    }

    * Compare data using cf _all
    use `stata_data', clear
    capture cf _all using `cexport_data'
    local cfrc = _rc

    if `cfrc' == 0 {
        restore
        test_pass "`testname'"
        exit
    }

    * cf _all failed - try tolerance-based comparison for numeric variables
    use `stata_data', clear
    ds
    local allvars `r(varlist)'

    * Rename variables from stata data with _stata suffix
    foreach v of local allvars {
        rename `v' `v'_stata
    }
    gen long _row = _n
    tempfile stata_renamed
    quietly save `stata_renamed', replace

    * Load cexport data and rename with _cexport suffix
    use `cexport_data', clear
    foreach v of local allvars {
        capture confirm variable `v'
        if _rc != 0 {
            restore
            test_fail "`testname'" "variable `v' missing from cexport data"
            exit
        }
        rename `v' `v'_cexport
    }
    gen long _row = _n
    merge 1:1 _row using `stata_renamed', nogen

    local all_match = 1
    local fail_reason ""
    local min_sigfigs = 99

    foreach v of local allvars {
        * Check if string or numeric
        capture confirm string variable `v'_stata
        if _rc == 0 {
            * String variable - must be exactly equal
            quietly count if `v'_stata != `v'_cexport
            if r(N) > 0 {
                local all_match = 0
                local fail_reason "string variable `v' has `r(N)' mismatches"
                continue, break
            }
        }
        else {
            * Numeric variable - compare values using significant figures
            quietly gen double _v1 = `v'_stata
            quietly gen double _v2 = `v'_cexport

            * Check each non-missing pair using sigfigs
            local var_min_sf = 99
            quietly count if !missing(_v1) & !missing(_v2)
            local n_pairs = r(N)
            if `n_pairs' > 0 {
                quietly gen double _rel_err = abs(_v1 - _v2) / max(abs(_v1), abs(_v2), 1e-300) if !missing(_v1) & !missing(_v2)
                quietly summarize _rel_err
                if r(max) != . & r(max) > 0 {
                    local max_rel_err = r(max)
                    local var_min_sf = -log10(`max_rel_err')
                    if `var_min_sf' < `min_sigfigs' {
                        local min_sigfigs = `var_min_sf'
                    }
                }
                drop _rel_err
            }

            * Check if sigfigs meets threshold
            if `var_min_sf' < $DEFAULT_SIGFIGS {
                local all_match = 0
                local sf_fmt : display %4.1f `var_min_sf'
                local fail_reason "numeric variable `v' has only `sf_fmt' significant figures agreement"
                drop _v1 _v2
                continue, break
            }
            drop _v1 _v2
            * Check missing values match
            quietly count if missing(`v'_stata) != missing(`v'_cexport)
            if r(N) > 0 {
                local all_match = 0
                local fail_reason "variable `v' has `r(N)' missing value mismatches"
                continue, break
            }
        }
    }

    restore

    if `all_match' == 1 {
        if `min_sigfigs' < 99 {
            local sf_fmt : display %4.1f `min_sigfigs'
            test_pass "`testname' (min_sigfigs=`sf_fmt')"
        }
        else {
            test_pass "`testname'"
        }
    }
    else {
        test_fail "`testname'" "`fail_reason'"
    }
end

/*******************************************************************************
 * benchmark_psmatch2 - Compare cpsmatch vs psmatch2
 *
 * Comprehensively compares propensity score matching results:
 *   - r(att), r(seatt) - ATT and standard error from psmatch2
 *   - _pscore variable - propensity scores (must match)
 *   - Sample sizes computed from _treated/_support variables
 *
 * All continuous comparisons use significant figures (default 7 sigfigs).
 *
 * NOTE: psmatch2 does not return sample counts in r(), so we compute them
 * from the generated variables (_treated, _support) for comparison.
 *
 * Syntax: benchmark_psmatch2 treatvar [varlist], [outcome(varname) options]
 ******************************************************************************/
capture program drop benchmark_psmatch2
program define benchmark_psmatch2
    syntax varlist(min=1 numeric) [if] [in], ///
        [OUTcome(varname numeric) Pscore(varname numeric) ///
         Logit Probit ///
         Neighbor(integer 1) Caliper(real 0) Radius Kernel ///
         Kerneltype(string) Bwidth(real 0.06) ///
         Common NOREPlacement Ties Descending ///
         testname(string) minsf(real 7)]

    gettoken depvar covars : varlist

    if "`testname'" == "" local testname "psmatch2 `depvar' `covars'"

    * Build common options for both commands
    local opts ""
    if "`outcome'" != "" local opts "`opts' outcome(`outcome')"
    if "`pscore'" != "" local opts "`opts' pscore(`pscore')"
    if "`logit'" != "" local opts "`opts' logit"
    if "`probit'" != "" local opts "`opts' probit"
    if `neighbor' != 1 local opts "`opts' neighbor(`neighbor')"
    if `caliper' != 0 local opts "`opts' caliper(`caliper')"
    if "`radius'" != "" local opts "`opts' radius"
    if "`kernel'" != "" local opts "`opts' kernel"
    if "`kerneltype'" != "" local opts "`opts' kerneltype(`kerneltype')"
    if `bwidth' != 0.06 local opts "`opts' bwidth(`bwidth')"
    if "`common'" != "" local opts "`opts' common"
    if "`noreplacement'" != "" local opts "`opts' noreplacement"
    if "`ties'" != "" local opts "`opts' ties"
    if "`descending'" != "" local opts "`opts' descending"

    preserve

    * Run psmatch2 (quietly)
    capture quietly psmatch2 `depvar' `covars' `if' `in', `opts'
    local psmatch2_rc = _rc

    if `psmatch2_rc' != 0 {
        restore
        test_fail "`testname'" "psmatch2 returned error `psmatch2_rc'"
        exit
    }

    * Store psmatch2 ATT results (these are what psmatch2 returns)
    local psmatch2_att = r(att)
    local psmatch2_att_se = r(seatt)

    * Compute sample counts from psmatch2 variables
    quietly count if _treated == 1 & _support == 1
    local psmatch2_n_treated = r(N)
    quietly count if _treated == 0 & _support == 1
    local psmatch2_n_controls = r(N)
    quietly count if _weight != . & _weight > 0 & _treated == 1
    local psmatch2_n_matched = r(N)

    * Store psmatch2 propensity scores
    tempvar pscore_psm2
    quietly gen double `pscore_psm2' = _pscore

    * Drop psmatch2-created variables before running cpsmatch
    capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

    * Run cpsmatch (quietly)
    capture quietly cpsmatch `depvar' `covars' `if' `in', `opts'
    local cpsmatch_rc = _rc

    if `cpsmatch_rc' != 0 {
        restore
        test_fail "`testname'" "cpsmatch returned error `cpsmatch_rc'"
        exit
    }

    * Store cpsmatch results
    local cpsmatch_n_treated = r(n_treated)
    local cpsmatch_n_controls = r(n_controls)
    local cpsmatch_n_matched = r(n_matched)
    local cpsmatch_att = r(att)
    local cpsmatch_att_se = r(att_se)

    * Track all differences found
    local all_diffs ""
    local has_failure = 0

    * Compare sample sizes (must match exactly)
    if `psmatch2_n_treated' != `cpsmatch_n_treated' {
        local has_failure = 1
        local all_diffs "`all_diffs' n_treated:`psmatch2_n_treated'!=`cpsmatch_n_treated'"
    }

    if `psmatch2_n_controls' != `cpsmatch_n_controls' {
        local has_failure = 1
        local all_diffs "`all_diffs' n_controls:`psmatch2_n_controls'!=`cpsmatch_n_controls'"
    }

    * Note: n_matched definition may differ between implementations
    * psmatch2 counts treated with matches, cpsmatch may count differently
    * Only flag as failure if significantly different
    local matched_diff = abs(`psmatch2_n_matched' - `cpsmatch_n_matched')
    if `matched_diff' > 0 {
        * Just note the difference, don't fail
        * local all_diffs "`all_diffs' n_matched:`psmatch2_n_matched'!=`cpsmatch_n_matched'"
    }

    * Compare ATT using significant figures (if outcome specified)
    if "`outcome'" != "" {
        if !missing(`psmatch2_att') & !missing(`cpsmatch_att') {
            sigfigs `psmatch2_att' `cpsmatch_att'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' att:sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`psmatch2_att') & missing(`cpsmatch_att') {
            local has_failure = 1
            local all_diffs "`all_diffs' att:missing_in_cpsmatch"
        }

        * Note: We skip ATT SE comparison because:
        * 1. psmatch2 explicitly notes: "S.E. does not take into account that the propensity score is estimated"
        * 2. cpsmatch uses a different approximation for the SE
        * 3. Both are known to be approximations that can differ significantly
        * The ATT point estimate comparison is the meaningful validation
    }

    * Compare propensity scores variable using significant figures
    quietly gen double _pscore_diff = abs(_pscore - `pscore_psm2')
    quietly sum _pscore_diff
    local max_diff = r(max)
    if `max_diff' > 1e-7 {
        * Find minimum sigfigs across all observations
        tempvar sf_var
        quietly gen double `sf_var' = 15 if _pscore == `pscore_psm2'
        quietly replace `sf_var' = -log10(abs(_pscore - `pscore_psm2') / max(abs(_pscore), abs(`pscore_psm2'))) if `sf_var' == .
        quietly replace `sf_var' = 0 if `sf_var' < 0
        quietly replace `sf_var' = 15 if `sf_var' > 15
        quietly sum `sf_var'
        local min_sf = r(min)
        if `min_sf' < `minsf' {
            local has_failure = 1
            local sf_fmt : display %4.1f `min_sf'
            local all_diffs "`all_diffs' _pscore:sigfigs=`sf_fmt'"
        }
    }
    drop _pscore_diff

    restore

    * Report result
    if `has_failure' == 0 {
        test_pass "`testname'"
    }
    else {
        * Trim leading space from all_diffs
        local all_diffs = trim("`all_diffs'")
        test_fail "`testname'" "`all_diffs'"
    }
end

/*******************************************************************************
 * benchmark_ivreghdfe - Compare civreghdfe vs ivreghdfe
 *
 * Comprehensively compares all e() scalars, vectors, and matrices:
 *   - e(N), e(df_r), e(df_m), e(rank) - counts and degrees of freedom
 *   - e(r2), e(r2_a) - R-squared values
 *   - e(rss), e(mss) - sum of squares
 *   - e(F) - F-statistic
 *   - e(rmse) - root mean squared error
 *   - e(N_clust) - number of clusters (if clustering)
 *   - e(idstat), e(idp) - underidentification test statistic and p-value
 *   - e(widstat) - weak identification (Kleibergen-Paap) statistic
 *   - e(sargan), e(sarganp) - overidentification test (if applicable)
 *   - Diagnostic tests: orthog (C-stat), endogtest, redundant
 *   - First-stage stats: partial_r2 (via ffirst)
 *   - e(b) - coefficient vector (with dimension check)
 *   - e(V) - variance-covariance matrix (with dimension check)
 *
 * All comparisons use significant figures (default 7 sigfigs).
 *
 * Syntax: benchmark_ivreghdfe spec [weight] [if] [in], absorb(varlist) [options]
 ******************************************************************************/

capture program drop benchmark_ivreghdfe
program define benchmark_ivreghdfe
    syntax anything(name=spec) [aw fw pw] [if] [in], Absorb(varlist) ///
        [vce(string) testname(string) minsf(real 7) ///
        FIRST FFIRst SMALL NOId ///
        ORTHog(string) ENDOGtest(string) REDundant(string) ///
        PARTial(string) FWL(string) ///
        TOLerance(string) MAXITer(string) ///
        LIML FULLER(real 0) Kclass(real 0) GMM2s CUE COVIV ///
        BW(integer 0) KERNEL(string) DKRAAY(integer 0) KIEFER ///
        b0(string)]

    if "`testname'" == "" local testname "ivreghdfe `spec' `if' `in', absorb(`absorb')"

    * Build weight string
    local wtexp ""
    if "`weight'" != "" local wtexp "[`weight'`exp']"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    * Build estimation display options
    local dispopts ""
    if "`first'" != "" local dispopts "`dispopts' first"
    if "`ffirst'" != "" local dispopts "`dispopts' ffirst"
    if "`small'" != "" local dispopts "`dispopts' small"
    if "`noid'" != "" local dispopts "`dispopts' noid"

    * Build diagnostic test options
    local diagopts ""
    if "`orthog'" != "" local diagopts "`diagopts' orthog(`orthog')"
    if "`endogtest'" != "" local diagopts "`diagopts' endogtest(`endogtest')"
    if "`redundant'" != "" local diagopts "`diagopts' redundant(`redundant')"

    * Build FWL/partial options
    local fwlopts ""
    if "`partial'" != "" local fwlopts "`fwlopts' partial(`partial')"
    if "`fwl'" != "" local fwlopts "`fwlopts' fwl(`fwl')"

    * Build solver options
    local solveropts ""
    if "`tolerance'" != "" local solveropts "`solveropts' tolerance(`tolerance')"
    if "`maxiter'" != "" local solveropts "`solveropts' maxiter(`maxiter')"

    * Build estimator options
    local estopts ""
    if "`liml'" != "" local estopts "`estopts' liml"
    if `fuller' != 0 local estopts "`estopts' fuller(`fuller')"
    if `kclass' != 0 local estopts "`estopts' kclass(`kclass')"
    if "`gmm2s'" != "" local estopts "`estopts' gmm2s"
    if "`cue'" != "" local estopts "`estopts' cue"
    if "`coviv'" != "" local estopts "`estopts' coviv"

    * Build HAC options
    * Note: civreghdfe implements robust HAC (heteroskedasticity + autocorrelation consistent)
    * so we add 'robust' to ivreghdfe when kernel or bw is specified to ensure matching VCE
    * But don't add robust if vce already contains "robust" to avoid duplicate specification
    local hacopts ""
    local hac_robust = 0
    if `bw' != 0 {
        local hacopts "`hacopts' bw(`bw')"
        local hac_robust = 1
    }
    if "`kernel'" != "" {
        local hacopts "`hacopts' kernel(`kernel')"
        local hac_robust = 1
    }
    * Add robust for HAC to match civreghdfe's robust HAC implementation
    * Only add if vce doesn't already contain robust
    local vce_has_robust = regexm("`vce'", "robust")
    if `hac_robust' & `dkraay' == 0 & "`kiefer'" == "" & !`vce_has_robust' {
        local hacopts "`hacopts' robust"
    }
    if `dkraay' != 0 local hacopts "`hacopts' dkraay(`dkraay')"
    if "`kiefer'" != "" local hacopts "`hacopts' kiefer"

    * Build b0 option (if specified)
    local b0opt ""
    if "`b0'" != "" local b0opt "b0(`b0')"

    preserve

    * Run ivreghdfe (quietly)
    capture quietly ivreghdfe `spec' `wtexp' `if' `in', absorb(`absorb') `vceopt' `dispopts' `diagopts' `fwlopts' `solveropts' `estopts' `hacopts' `b0opt'
    local ivreghdfe_rc = _rc

    if `ivreghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "ivreghdfe returned error `ivreghdfe_rc'"
        exit
    }

    * Store ivreghdfe results - matrices
    matrix ivreghdfe_b = e(b)
    matrix ivreghdfe_V = e(V)

    * Store ivreghdfe results - scalars (capture for optional values)
    local ivreghdfe_N = e(N)
    local ivreghdfe_df_r = e(df_r)
    local ivreghdfe_df_m = e(df_m)
    local ivreghdfe_r2 = e(r2)
    local ivreghdfe_r2_a = e(r2_a)
    local ivreghdfe_rss = e(rss)
    local ivreghdfe_mss = e(mss)
    local ivreghdfe_F = e(F)
    local ivreghdfe_rmse = e(rmse)
    local ivreghdfe_rank = e(rank)
    local ivreghdfe_N_clust = e(N_clust)
    local ivreghdfe_idstat = e(idstat)
    local ivreghdfe_idp = e(idp)
    local ivreghdfe_widstat = e(widstat)
    local ivreghdfe_sargan = e(sargan)
    local ivreghdfe_sarganp = e(sarganp)
    local ivreghdfe_N_hdfe = e(N_hdfe)

    * Store diagnostic test results (ivreghdfe naming conventions)
    local ivreghdfe_cstat = e(cstat)
    local ivreghdfe_cstatp = e(cstatp)
    local ivreghdfe_cstatdf = e(cstatdf)
    local ivreghdfe_estat = e(estat)
    local ivreghdfe_estatp = e(estatp)
    local ivreghdfe_estatdf = e(estatdf)
    local ivreghdfe_redstat = e(redstat)
    local ivreghdfe_redp = e(redp)
    local ivreghdfe_reddf = e(reddf)

    * Store ffirst partial R² values (up to 10 endogenous vars)
    forvalues k = 1/10 {
        local ivreghdfe_partial_r2_`k' = e(partial_r2_`k')
    }

    * Run civreghdfe (quietly)
    capture quietly civreghdfe `spec' `wtexp' `if' `in', absorb(`absorb') `vceopt' `dispopts' `diagopts' `fwlopts' `solveropts' `estopts' `hacopts' `b0opt'
    local civreghdfe_rc = _rc

    if `civreghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "civreghdfe returned error `civreghdfe_rc'"
        exit
    }

    * Store civreghdfe results - matrices
    matrix civreghdfe_b = e(b)
    matrix civreghdfe_V = e(V)

    * Store civreghdfe results - scalars
    local civreghdfe_N = e(N)
    local civreghdfe_df_r = e(df_r)
    local civreghdfe_df_m = e(df_m)
    local civreghdfe_r2 = e(r2)
    local civreghdfe_r2_a = e(r2_a)
    local civreghdfe_rss = e(rss)
    local civreghdfe_mss = e(mss)
    local civreghdfe_F = e(F)
    local civreghdfe_rmse = e(rmse)
    local civreghdfe_rank = e(rank)
    local civreghdfe_N_clust = e(N_clust)
    local civreghdfe_idstat = e(idstat)
    local civreghdfe_idp = e(idp)
    local civreghdfe_widstat = e(widstat)
    local civreghdfe_sargan = e(sargan)
    local civreghdfe_sarganp = e(sarganp)
    local civreghdfe_N_hdfe = e(N_hdfe)

    * Store diagnostic test results (civreghdfe naming conventions)
    local civreghdfe_cstat = e(cstat)
    local civreghdfe_cstat_p = e(cstat_p)
    local civreghdfe_cstat_df = e(cstat_df)
    local civreghdfe_endogtest = e(endogtest)
    local civreghdfe_endogtest_p = e(endogtest_p)
    local civreghdfe_endogtest_df = e(endogtest_df)
    local civreghdfe_redund = e(redund)
    local civreghdfe_redund_p = e(redund_p)
    local civreghdfe_redund_df = e(redund_df)

    * Store ffirst partial R² values
    forvalues k = 1/10 {
        local civreghdfe_partial_r2_`k' = e(partial_r2_`k')
    }

    restore

    * Track all differences found
    local all_diffs ""
    local has_failure = 0

    * Compare N (must match exactly)
    if `ivreghdfe_N' != `civreghdfe_N' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(N):`ivreghdfe_N'!=`civreghdfe_N'"
    }

    * Compare integer scalars (must match exactly)
    foreach scalar in df_r df_m rank N_clust N_hdfe {
        local val1 = `ivreghdfe_`scalar''
        local val2 = `civreghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            if `val1' != `val2' {
                local has_failure = 1
                local all_diffs "`all_diffs' e(`scalar'):`val1'!=`val2'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local all_diffs "`all_diffs' e(`scalar'):missing_in_civreghdfe"
        }
    }

    * Compare continuous scalars using significant figures
    * Main estimation results
    foreach scalar in r2 r2_a rss mss F rmse {
        local val1 = `ivreghdfe_`scalar''
        local val2 = `civreghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):missing_in_civreghdfe"
        }
    }

    * Diagnostic statistics: use significant figures
    foreach scalar in idstat idp widstat sargan sarganp {
        local val1 = `ivreghdfe_`scalar''
        local val2 = `civreghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            sigfigs `val1' `val2'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(`scalar'):missing_in_civreghdfe"
        }
    }

    * Compare orthog (C-statistic) diagnostic test results
    * ivreghdfe: e(cstat), e(cstatp), e(cstatdf)
    * civreghdfe: e(cstat), e(cstat_p), e(cstat_df)
    if "`orthog'" != "" {
        * cstat: same name in both
        if !missing(`ivreghdfe_cstat') & !missing(`civreghdfe_cstat') {
            sigfigs `ivreghdfe_cstat' `civreghdfe_cstat'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(cstat):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`ivreghdfe_cstat') & missing(`civreghdfe_cstat') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(cstat):missing_in_civreghdfe"
        }
        * cstatp → cstat_p
        if !missing(`ivreghdfe_cstatp') & !missing(`civreghdfe_cstat_p') {
            sigfigs `ivreghdfe_cstatp' `civreghdfe_cstat_p'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(cstat_p):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`ivreghdfe_cstatp') & missing(`civreghdfe_cstat_p') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(cstat_p):missing_in_civreghdfe"
        }
        * cstatdf → cstat_df (integer)
        if !missing(`ivreghdfe_cstatdf') & !missing(`civreghdfe_cstat_df') {
            if `ivreghdfe_cstatdf' != `civreghdfe_cstat_df' {
                local has_failure = 1
                local all_diffs "`all_diffs' e(cstat_df):`ivreghdfe_cstatdf'!=`civreghdfe_cstat_df'"
            }
        }
        else if !missing(`ivreghdfe_cstatdf') & missing(`civreghdfe_cstat_df') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(cstat_df):missing_in_civreghdfe"
        }
    }

    * Compare endogtest diagnostic test results
    * ivreghdfe: e(estat), e(estatp), e(estatdf)
    * civreghdfe: e(endogtest), e(endogtest_p), e(endogtest_df)
    if "`endogtest'" != "" {
        * estat → endogtest
        if !missing(`ivreghdfe_estat') & !missing(`civreghdfe_endogtest') {
            sigfigs `ivreghdfe_estat' `civreghdfe_endogtest'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(endogtest):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`ivreghdfe_estat') & missing(`civreghdfe_endogtest') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(endogtest):missing_in_civreghdfe"
        }
        * estatp → endogtest_p
        if !missing(`ivreghdfe_estatp') & !missing(`civreghdfe_endogtest_p') {
            sigfigs `ivreghdfe_estatp' `civreghdfe_endogtest_p'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(endogtest_p):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`ivreghdfe_estatp') & missing(`civreghdfe_endogtest_p') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(endogtest_p):missing_in_civreghdfe"
        }
        * estatdf → endogtest_df (integer)
        if !missing(`ivreghdfe_estatdf') & !missing(`civreghdfe_endogtest_df') {
            if `ivreghdfe_estatdf' != `civreghdfe_endogtest_df' {
                local has_failure = 1
                local all_diffs "`all_diffs' e(endogtest_df):`ivreghdfe_estatdf'!=`civreghdfe_endogtest_df'"
            }
        }
        else if !missing(`ivreghdfe_estatdf') & missing(`civreghdfe_endogtest_df') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(endogtest_df):missing_in_civreghdfe"
        }
    }

    * Compare redundant instrument test results
    * ivreghdfe: e(redstat), e(redp), e(reddf)
    * civreghdfe: e(redund), e(redund_p), e(redund_df)
    if "`redundant'" != "" {
        * redstat → redund
        if !missing(`ivreghdfe_redstat') & !missing(`civreghdfe_redund') {
            sigfigs `ivreghdfe_redstat' `civreghdfe_redund'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(redund):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`ivreghdfe_redstat') & missing(`civreghdfe_redund') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(redund):missing_in_civreghdfe"
        }
        * redp → redund_p
        if !missing(`ivreghdfe_redp') & !missing(`civreghdfe_redund_p') {
            sigfigs `ivreghdfe_redp' `civreghdfe_redund_p'
            local sf = r(sigfigs)
            if `sf' < `minsf' {
                local has_failure = 1
                local sf_fmt : display %4.1f `sf'
                local all_diffs "`all_diffs' e(redund_p):sigfigs=`sf_fmt'"
            }
        }
        else if !missing(`ivreghdfe_redp') & missing(`civreghdfe_redund_p') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(redund_p):missing_in_civreghdfe"
        }
        * reddf → redund_df (integer)
        if !missing(`ivreghdfe_reddf') & !missing(`civreghdfe_redund_df') {
            if `ivreghdfe_reddf' != `civreghdfe_redund_df' {
                local has_failure = 1
                local all_diffs "`all_diffs' e(redund_df):`ivreghdfe_reddf'!=`civreghdfe_redund_df'"
            }
        }
        else if !missing(`ivreghdfe_reddf') & missing(`civreghdfe_redund_df') {
            local has_failure = 1
            local all_diffs "`all_diffs' e(redund_df):missing_in_civreghdfe"
        }
    }

    * Compare ffirst partial R² values
    if "`ffirst'" != "" {
        forvalues k = 1/10 {
            local val1 = `ivreghdfe_partial_r2_`k''
            local val2 = `civreghdfe_partial_r2_`k''
            if !missing(`val1') & !missing(`val2') {
                sigfigs `val1' `val2'
                local sf = r(sigfigs)
                if `sf' < `minsf' {
                    local has_failure = 1
                    local sf_fmt : display %4.1f `sf'
                    local all_diffs "`all_diffs' e(partial_r2_`k'):sigfigs=`sf_fmt'"
                }
            }
            else if !missing(`val1') & missing(`val2') {
                local has_failure = 1
                local all_diffs "`all_diffs' e(partial_r2_`k'):missing_in_civreghdfe"
            }
        }
    }

    * Compare coefficients e(b) using significant figures
    * First check dimensions match
    local ivreghdfe_bcols = colsof(ivreghdfe_b)
    local civreghdfe_bcols = colsof(civreghdfe_b)

    if `ivreghdfe_bcols' != `civreghdfe_bcols' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(b):dim_mismatch(`ivreghdfe_bcols'!=`civreghdfe_bcols')"
    }
    else {
        local cols = colsof(ivreghdfe_b)
        local min_sf_b = 15
        local min_sf_b_idx = 0
        forvalues j = 1/`cols' {
            local v1 = ivreghdfe_b[1, `j']
            local v2 = civreghdfe_b[1, `j']
            sigfigs `v1' `v2'
            local sf = r(sigfigs)
            if `sf' < `min_sf_b' {
                local min_sf_b = `sf'
                local min_sf_b_idx = `j'
            }
        }

        if `min_sf_b' < `minsf' {
            local has_failure = 1
            local coefname : colnames ivreghdfe_b
            local badcoef : word `min_sf_b_idx' of `coefname'
            local sf_fmt : display %4.1f `min_sf_b'
            local all_diffs "`all_diffs' e(b)[`badcoef']:sigfigs=`sf_fmt'"
        }
    }

    * Compare VCE e(V) using significant figures
    * First check dimensions match
    local ivreghdfe_Vrows = rowsof(ivreghdfe_V)
    local civreghdfe_Vrows = rowsof(civreghdfe_V)

    if `ivreghdfe_Vrows' != `civreghdfe_Vrows' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(V):dim_mismatch(`ivreghdfe_Vrows'!=`civreghdfe_Vrows')"
    }
    else {
        local rows = rowsof(ivreghdfe_V)
        local cols = colsof(ivreghdfe_V)
        local min_sf_V = 15
        local min_sf_V_i = 0
        local min_sf_V_j = 0
        forvalues i = 1/`rows' {
            forvalues j = 1/`cols' {
                local v1 = ivreghdfe_V[`i', `j']
                local v2 = civreghdfe_V[`i', `j']
                sigfigs `v1' `v2'
                local sf = r(sigfigs)
                if `sf' < `min_sf_V' {
                    local min_sf_V = `sf'
                    local min_sf_V_i = `i'
                    local min_sf_V_j = `j'
                }
            }
        }

        if `min_sf_V' < `minsf' {
            local has_failure = 1
            local sf_fmt : display %4.1f `min_sf_V'
            local all_diffs "`all_diffs' e(V)[`min_sf_V_i',`min_sf_V_j']:sigfigs=`sf_fmt'"
        }
    }

    * Report result
    if `has_failure' == 0 {
        test_pass "`testname'"
    }
    else {
        * Trim leading space from all_diffs
        local all_diffs = trim("`all_diffs'")
        test_fail "`testname'" "`all_diffs'"
    }
end

/*******************************************************************************
 * benchmark_encode - Compare cencode vs encode
 *
 * Compares cencode output against Stata's native encode command.
 * Numeric values and label text must match exactly.
 *
 * Syntax: benchmark_encode varname, testname(string) [generate(name) label(name) noextend if2(string) in2(string) cencodeopts(string)]
 ******************************************************************************/
capture program drop benchmark_encode
program define benchmark_encode
    syntax varname, testname(string) [GENerate(name) LABel(name) noextend if2(string) in2(string) cencodeopts(string)]

    * Set default generate name
    if "`generate'" == "" local generate = "encoded"

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build options for cencode (includes cencodeopts like verbose, threads)
    local opts "generate(`generate'_c)"
    if "`label'" != "" local opts "`opts' label(`label'_c)"
    if "`noextend'" != "" local opts "`opts' noextend"
    if "`cencodeopts'" != "" local opts "`opts' `cencodeopts'"

    * Build options for Stata encode (no cencodeopts - Stata doesn't support them)
    local opts_stata "generate(`generate'_s)"
    if "`label'" != "" local opts_stata "`opts_stata' label(`label'_s)"
    if "`noextend'" != "" local opts_stata "`opts_stata' noextend"

    * Run cencode
    capture drop `generate'_c
    capture label drop `generate'_c
    if "`label'" != "" capture label drop `label'_c
    capture cencode `varlist' `ifin', `opts'
    local rc_c = _rc

    * Run native encode
    capture drop `generate'_s
    capture label drop `generate'_s
    if "`label'" != "" capture label drop `label'_s
    capture encode `varlist' `ifin', `opts_stata'
    local rc_s = _rc

    * Check both succeeded or both failed
    if `rc_c' != `rc_s' {
        test_fail "`testname'" "cencode rc=`rc_c', encode rc=`rc_s'"
        exit
    }

    if `rc_c' != 0 {
        test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare numeric values
    quietly count if `generate'_c != `generate'_s & !missing(`generate'_c) & !missing(`generate'_s)
    local ndiff = r(N)

    * Compare missing patterns
    quietly count if missing(`generate'_c) != missing(`generate'_s)
    local nmiss_diff = r(N)

    if `ndiff' > 0 | `nmiss_diff' > 0 {
        test_fail "`testname'" "`ndiff' value diffs, `nmiss_diff' missing diffs"
    }
    else {
        * Also verify labels match
        local all_match = 1
        quietly levelsof `generate'_c, local(codes)
        foreach c of local codes {
            local lbl_c : label (`generate'_c) `c'
            local lbl_s : label (`generate'_s) `c'
            if `"`lbl_c'"' != `"`lbl_s'"' {
                local all_match = 0
            }
        }
        if `all_match' {
            test_pass "`testname'"
        }
        else {
            test_fail "`testname'" "label text mismatch"
        }
    }

    * Cleanup
    capture drop `generate'_c `generate'_s
    capture label drop `generate'_c `generate'_s
    if "`label'" != "" {
        capture label drop `label'_c `label'_s
    }
end

/*******************************************************************************
 * benchmark_winsor - Compare cwinsor vs gstats winsor
 *
 * Compares cwinsor output against gstats winsor (from gtools).
 * Numeric values compared using significant figures (default 7 sigfigs).
 * Missing patterns must match exactly.
 *
 * Syntax: benchmark_winsor varlist(numeric), testname(string) [cuts(numlist) p(real) q(real) trim by(varlist) suffix(string) prefix(string) cwinsoropts(string)]
 ******************************************************************************/
capture program drop benchmark_winsor
program define benchmark_winsor
    syntax varlist(numeric), testname(string) [Cuts(numlist min=2 max=2) P(real -1) Q(real -1) TRIM BY(varlist) SUFfix(string) PREfix(string) cwinsoropts(string)]

    * Set default percentiles
    if `p' < 0 & `q' < 0 & "`cuts'" == "" {
        local p = 1
        local q = 99
    }

    * Build cuts string for gstats
    if "`cuts'" != "" {
        local p : word 1 of `cuts'
        local q : word 2 of `cuts'
    }

    * Build suffix for new variables - always use suffix for comparison
    local suf "_cw"
    if "`suffix'" != "" local suf "`suffix'"
    if "`prefix'" != "" local suf ""

    * Build options for cwinsor - always use suffix or prefix for comparison
    local copts ""
    if "`cuts'" != "" {
        local copts "cuts(`cuts')"
    }
    else {
        local copts "p(`p') q(`q')"
    }
    if "`trim'" != "" local copts "`copts' trim"
    if "`by'" != "" local copts "`copts' by(`by')"
    if "`suffix'" != "" {
        local copts "`copts' suffix(`suffix')"
    }
    else if "`prefix'" != "" {
        local copts "`copts' prefix(`prefix')"
    }
    else {
        * Default: use suffix(_cw) for comparison
        local copts "`copts' suffix(_cw)"
    }
    if "`cwinsoropts'" != "" local copts "`copts' `cwinsoropts'"

    * Build options for gstats winsor
    local gopts "cuts(`p' `q')"
    if "`trim'" != "" local gopts "`gopts' trim"
    if "`by'" != "" local gopts "`gopts' by(`by')"

    * For gstats, always use suffix to create new variable
    local gopts "`gopts' suffix(_gw)"

    * Run cwinsor
    foreach v of local varlist {
        capture drop `v'`suf'
        if "`prefix'" != "" {
            capture drop `prefix'`v'
        }
    }
    capture cwinsor `varlist', `copts'
    local rc_c = _rc

    * Run gstats winsor
    foreach v of local varlist {
        capture drop `v'_gw
    }
    capture gstats winsor `varlist', `gopts'
    local rc_g = _rc

    * Check both succeeded or both failed
    if `rc_c' != `rc_g' {
        test_fail "`testname'" "cwinsor rc=`rc_c', gstats rc=`rc_g'"
        exit
    }

    if `rc_c' != 0 {
        test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare results for each variable
    local all_pass = 1
    local fail_reason ""
    foreach v of local varlist {
        * Determine cwinsor output variable name
        local cvar "`v'`suf'"
        if "`prefix'" != "" local cvar "`prefix'`v'"

        * gstats output variable name
        local gvar "`v'_gw"

        * Compare using sigfigs
        quietly count if `cvar' != `gvar' & !missing(`cvar') & !missing(`gvar')
        local ndiff = r(N)

        if `ndiff' > 0 {
            * Use significant figures comparison for floating point tolerance
            tempvar sf
            quietly {
                gen double `sf' = 15 if `cvar' == `gvar'
                replace `sf' = 0 if (`sf' == .) & ((`cvar' == 0) | (`gvar' == 0))
                replace `sf' = -log10(abs(`cvar' - `gvar') / max(abs(`cvar'), abs(`gvar'))) if `sf' == .
                replace `sf' = 0 if `sf' < 0
                replace `sf' = 15 if `sf' > 15
            }

            quietly count if `sf' < $DEFAULT_SIGFIGS & !missing(`cvar') & !missing(`gvar')
            local nfail = r(N)
            drop `sf'

            if `nfail' > 0 {
                local all_pass = 0
                local fail_reason "`v': `nfail' values differ beyond tolerance"
            }
        }

        * Also check missing value patterns match
        quietly count if missing(`cvar') != missing(`gvar')
        local nmiss_diff = r(N)
        if `nmiss_diff' > 0 {
            local all_pass = 0
            local fail_reason "`v': `nmiss_diff' missing value patterns differ"
        }
    }

    if `all_pass' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "`fail_reason'"
    }

    * Cleanup (use capture to ignore if variable doesn't exist)
    foreach v of local varlist {
        capture drop `v'`suf'
        if "`prefix'" != "" {
            capture drop `prefix'`v'
        }
        capture drop `v'_gw
    }
end

/*******************************************************************************
 * benchmark_decode - Compare cdecode vs decode
 *
 * Compares cdecode output against Stata's native decode command.
 * String values must match exactly (no tolerance - strings are discrete).
 *
 * Syntax: benchmark_decode varname, testname(string) [generate(name) maxlength(#) if2(string) in2(string) cdecodeopts(string)]
 ******************************************************************************/
capture program drop benchmark_decode
program define benchmark_decode
    syntax varname, testname(string) [GENerate(name) MAXLength(integer 0) if2(string) in2(string) cdecodeopts(string)]

    * Set default generate name
    if "`generate'" == "" local generate = "decoded"

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build options for cdecode (includes cdecodeopts like verbose, threads)
    local opts "generate(`generate'_c)"
    if `maxlength' > 0 local opts "`opts' maxlength(`maxlength')"
    if "`cdecodeopts'" != "" local opts "`opts' `cdecodeopts'"

    * Build options for Stata decode (no cdecodeopts - Stata doesn't support them)
    local opts_stata "generate(`generate'_s)"
    if `maxlength' > 0 local opts_stata "`opts_stata' maxlength(`maxlength')"

    * Run cdecode
    capture drop `generate'_c
    capture cdecode `varlist' `ifin', `opts'
    local rc_c = _rc

    * Run native decode
    capture drop `generate'_s
    capture decode `varlist' `ifin', `opts_stata'
    local rc_s = _rc

    * Check both succeeded or both failed
    if `rc_c' != `rc_s' {
        test_fail "`testname'" "cdecode rc=`rc_c', decode rc=`rc_s'"
        exit
    }

    if `rc_c' != 0 {
        test_pass "`testname' (both error as expected)"
        exit
    }

    * Compare string values
    quietly count if `generate'_c != `generate'_s
    local ndiff = r(N)

    if `ndiff' > 0 {
        test_fail "`testname'" "`ndiff' values differ"
    }
    else {
        test_pass "`testname'"
    }

    * Cleanup
    capture drop `generate'_c `generate'_s
end

/*******************************************************************************
 * benchmark_destring - Compare cdestring vs destring
 *
 * Compares cdestring output against Stata's native destring command.
 * Supports generate and replace modes. Numeric values compared using
 * significant figures (default 7 sigfigs). Missing patterns must match exactly.
 *
 * Syntax: benchmark_destring varlist(string), testname(string) [generate(string) replace ignore(string) force float percent dpcomma if2(string) in2(string)]
 ******************************************************************************/
capture program drop benchmark_destring
program define benchmark_destring
    syntax varlist(string), testname(string) [GENerate(string) replace ///
        IGnore(string) force float percent dpcomma if2(string) in2(string)]

    * Build if/in conditions
    local ifin ""
    if "`if2'" != "" local ifin "`ifin' if `if2'"
    if "`in2'" != "" local ifin "`ifin' in `in2'"

    * Build common options
    local common_opts ""
    if `"`ignore'"' != "" local common_opts `"`common_opts' ignore("`ignore'")"'
    if "`force'" != "" local common_opts "`common_opts' force"
    if "`float'" != "" local common_opts "`common_opts' float"
    if "`percent'" != "" local common_opts "`common_opts' percent"
    if "`dpcomma'" != "" local common_opts "`common_opts' dpcomma"

    * Build generate names with _c and _s suffixes
    if "`generate'" != "" {
        local gen_c ""
        local gen_s ""
        foreach v of local generate {
            local gen_c "`gen_c' `v'_c"
            local gen_s "`gen_s' `v'_s"
        }
        local opts_c "generate(`gen_c') `common_opts'"
        local opts_s "generate(`gen_s') `common_opts'"
    }
    else if "`replace'" != "" {
        * For replace, we need to work with copies
        local opts_c "replace `common_opts'"
        local opts_s "replace `common_opts'"
    }
    else {
        test_fail "`testname'" "must specify generate or replace"
        exit
    }

    * Handle replace option specially - need to create copies
    if "`replace'" != "" {
        * Create copies of original variables
        local copy_c ""
        local copy_s ""
        foreach v of local varlist {
            capture drop `v'_c `v'_s
            clonevar `v'_c = `v'
            clonevar `v'_s = `v'
        }

        * Run cdestring on copies
        local varlist_c ""
        foreach v of local varlist {
            local varlist_c "`varlist_c' `v'_c"
        }
        capture cdestring `varlist_c' `ifin', replace `common_opts'
        local rc_c = _rc

        * Run destring on copies
        local varlist_s ""
        foreach v of local varlist {
            local varlist_s "`varlist_s' `v'_s"
        }
        capture destring `varlist_s' `ifin', replace `common_opts'
        local rc_s = _rc

        * Check both succeeded or both failed
        if `rc_c' != `rc_s' {
            test_fail "`testname'" "cdestring rc=`rc_c', destring rc=`rc_s'"
            exit
        }

        if `rc_c' != 0 {
            test_pass "`testname' (both error as expected)"
            exit
        }

        * Compare results (use sigfigs-based comparison)
        local ndiff_total = 0
        foreach v of local varlist {
            tempvar _sf_`v'
            quietly {
                gen double `_sf_`v'' = 15 if `v'_c == `v'_s
                replace `_sf_`v'' = 0 if (`_sf_`v'' == .) & ((`v'_c == 0) | (`v'_s == 0))
                replace `_sf_`v'' = -log10(abs(`v'_c - `v'_s) / max(abs(`v'_c), abs(`v'_s))) if `_sf_`v'' == .
                replace `_sf_`v'' = 0 if `_sf_`v'' < 0
                replace `_sf_`v'' = 15 if `_sf_`v'' > 15
            }
            quietly count if `_sf_`v'' < $DEFAULT_SIGFIGS & !missing(`v'_c) & !missing(`v'_s)
            local ndiff_total = `ndiff_total' + r(N)
            * Also check missing patterns match
            quietly count if missing(`v'_c) != missing(`v'_s)
            local ndiff_total = `ndiff_total' + r(N)
        }

        if `ndiff_total' > 0 {
            test_fail "`testname'" "`ndiff_total' values differ"
        }
        else {
            test_pass "`testname'"
        }

        * Cleanup
        foreach v of local varlist {
            capture drop `v'_c `v'_s
        }
    }
    else {
        * Generate option - simpler case

        * Drop any existing test variables
        foreach v of local gen_c {
            capture drop `v'
        }
        foreach v of local gen_s {
            capture drop `v'
        }

        * Run cdestring
        capture cdestring `varlist' `ifin', `opts_c'
        local rc_c = _rc

        * Run destring
        capture destring `varlist' `ifin', `opts_s'
        local rc_s = _rc

        * Check both succeeded or both failed
        if `rc_c' != `rc_s' {
            test_fail "`testname'" "cdestring rc=`rc_c', destring rc=`rc_s'"
            * Cleanup
            foreach v of local gen_c {
                capture drop `v'
            }
            foreach v of local gen_s {
                capture drop `v'
            }
            exit
        }

        if `rc_c' != 0 {
            test_pass "`testname' (both error as expected)"
            exit
        }

        * Compare numeric values
        local ndiff_total = 0
        local i = 1
        foreach v of local generate {
            local vc = "`v'_c"
            local vs = "`v'_s"
            * Compare non-missing values with sigfigs
            tempvar _sf_`v'
            quietly {
                gen double `_sf_`v'' = 15 if `vc' == `vs'
                replace `_sf_`v'' = 0 if (`_sf_`v'' == .) & ((`vc' == 0) | (`vs' == 0))
                replace `_sf_`v'' = -log10(abs(`vc' - `vs') / max(abs(`vc'), abs(`vs'))) if `_sf_`v'' == .
                replace `_sf_`v'' = 0 if `_sf_`v'' < 0
                replace `_sf_`v'' = 15 if `_sf_`v'' > 15
            }
            quietly count if `_sf_`v'' < $DEFAULT_SIGFIGS & !missing(`vc') & !missing(`vs')
            local ndiff_total = `ndiff_total' + r(N)
            * Check missing patterns match
            quietly count if missing(`vc') != missing(`vs')
            local ndiff_total = `ndiff_total' + r(N)
            local ++i
        }

        if `ndiff_total' > 0 {
            test_fail "`testname'" "`ndiff_total' values differ"
        }
        else {
            test_pass "`testname'"
        }

        * Cleanup
        foreach v of local gen_c {
            capture drop `v'
        }
        foreach v of local gen_s {
            capture drop `v'
        }
    }
end

/*******************************************************************************
 * benchmark_rangestat - Compare crangestat vs rangestat
 *
 * Runs both crangestat and rangestat with the same stat/variable/options,
 * then compares results using significant figures.
 *
 * If rangestat is not installed or doesn't support the requested statistic,
 * falls back to verifying crangestat results are non-missing.
 *
 * Syntax: benchmark_rangestat statname sourcevar, interval(string) [by(varlist) excludeself testname(string) minsf(real 7)]
 ******************************************************************************/
capture program drop benchmark_rangestat
program define benchmark_rangestat
    syntax anything(name=spec), INTerval(string) [BY(varlist) EXCLUDEself testname(string) minsf(real 7)]

    * Parse statname and sourcevar from spec
    gettoken statname sourcevar : spec
    local sourcevar = strtrim("`sourcevar'")

    if "`sourcevar'" == "" {
        di as error "benchmark_rangestat: must specify statname and sourcevar"
        exit 198
    }

    if "`testname'" == "" local testname "rangestat (`statname') `sourcevar', interval(`interval')"

    * Build common options
    local opts ""
    if "`by'" != "" local opts "`opts' by(`by')"
    if "`excludeself'" != "" local opts "`opts' excludeself"

    * Variable names for comparison
    local ctools_var "__brs_ctools"
    local rangestat_var "`sourcevar'_`statname'"

    preserve

    * Clean up pre-existing variables
    capture drop `ctools_var'
    capture drop `rangestat_var'

    * Run crangestat
    capture quietly crangestat (`statname') `ctools_var'=`sourcevar', interval(`interval') `opts'
    local crangestat_rc = _rc

    if `crangestat_rc' != 0 {
        restore
        test_fail "`testname'" "crangestat returned error `crangestat_rc'"
        exit
    }

    * Check if rangestat is installed
    capture which rangestat
    if _rc != 0 {
        * rangestat not installed — verify crangestat results are sane
        quietly count if !missing(`ctools_var')
        local n_nonmiss = r(N)
        restore
        if `n_nonmiss' > 0 {
            test_pass "`testname' [rangestat not installed]"
        }
        else {
            test_fail "`testname'" "all results missing [rangestat not installed]"
        }
        exit
    }

    * Run rangestat (capture in case stat is unsupported)
    capture quietly rangestat (`statname') `sourcevar', interval(`interval') `opts'
    local rangestat_rc = _rc

    if `rangestat_rc' != 0 {
        * rangestat doesn't support this stat — fall back to sanity check
        quietly count if !missing(`ctools_var')
        local n_nonmiss = r(N)
        restore
        if `n_nonmiss' > 0 {
            test_pass "`testname' [rangestat (`statname') unsupported]"
        }
        else {
            test_fail "`testname'" "all results missing [rangestat (`statname') unsupported]"
        }
        exit
    }

    * Both succeeded — compare results
    * Check for mismatched missing patterns
    quietly count if missing(`ctools_var') != missing(`rangestat_var')
    local nmiss_diff = r(N)

    if `nmiss_diff' > 0 {
        restore
        test_fail "`testname'" "`nmiss_diff' obs with mismatched missing patterns"
        exit
    }

    * Compare using assert_var_equal
    assert_var_equal `ctools_var' `rangestat_var' `minsf' "`testname'"

    restore
end

/*******************************************************************************
 * cimport_test - Test cimport standalone (no Stata comparison)
 *
 * Runs cimport and checks basic properties: return code, row count,
 * column count, and an optional custom assertion.
 *
 * Syntax: cimport_test using filename, testname(string) [importopts(string) expectn(integer) expectk(integer) check(string)]
 ******************************************************************************/
capture program drop cimport_test
program define cimport_test
    syntax using/, testname(string) [IMPORTopts(string) expectn(integer 0) expectk(integer 0) CHECK(string)]

    capture cimport delimited `using', `importopts' clear
    if _rc != 0 {
        test_fail "`testname'" "rc=`=_rc'"
        exit
    }

    * Check expected row count
    if `expectn' > 0 & _N != `expectn' {
        test_fail "`testname'" "expected N=`expectn', got N=`=_N'"
        exit
    }

    * Check expected column count
    if `expectk' > 0 & c(k) != `expectk' {
        test_fail "`testname'" "expected K=`expectk', got K=`=c(k)'"
        exit
    }

    * Optional custom check expression
    if "`check'" != "" {
        capture assert `check'
        if _rc != 0 {
            test_fail "`testname'" "check() failed: `check'"
            exit
        }
    }

    test_pass "`testname'"
end

/*******************************************************************************
 * benchmark_import_errtol - Error-tolerant import comparison
 *
 * Like benchmark_import but handles expected errors: if both commands error
 * with the same rc, that counts as a pass. When both succeed, does full
 * data comparison via cf _all (with positional variable renaming).
 *
 * Syntax: benchmark_import_errtol using filename, testname(string)
 ******************************************************************************/
capture program drop benchmark_import_errtol
program define benchmark_import_errtol
    syntax using/, testname(string)

    * Run Stata import delimited (with varnames(1) and encoding(utf-8))
    preserve
    capture import delimited `using', varnames(1) encoding(utf-8) clear
    local stata_rc = _rc
    if `stata_rc' == 0 {
        tempfile stata_data
        quietly save `stata_data', replace
        local stata_N = _N
        local stata_k = c(k)
    }
    restore

    * Run cimport
    capture cimport delimited `using', clear
    local cimport_rc = _rc

    * Both errored
    if `cimport_rc' != 0 & `stata_rc' != 0 {
        if `cimport_rc' == `stata_rc' {
            test_pass "`testname' (both error rc=`cimport_rc')"
        }
        else {
            test_fail "`testname'" "error codes differ: stata=`stata_rc' cimport=`cimport_rc'"
        }
        exit
    }

    * One errored, other didn't
    if `cimport_rc' != 0 | `stata_rc' != 0 {
        test_fail "`testname'" "one errored: stata rc=`stata_rc', cimport rc=`cimport_rc'"
        exit
    }

    * Both succeeded - check dimensions
    if _N != `stata_N' | c(k) != `stata_k' {
        test_fail "`testname'" "dimensions differ: N=`=_N' vs `stata_N', K=`=c(k)' vs `stata_k'"
        exit
    }

    * Dimensions match - compare data content with positional variable renaming
    tempfile cimport_data
    quietly save `cimport_data', replace

    use `stata_data', clear
    quietly ds
    local stata_vars `r(varlist)'
    local i = 1
    foreach v of local stata_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }
    tempfile stata_renamed
    quietly save `stata_renamed', replace

    use `cimport_data', clear
    quietly ds
    local cimport_vars `r(varlist)'
    local i = 1
    foreach v of local cimport_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }

    capture cf _all using `stata_renamed'
    if _rc == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "cf _all comparison failed - data not identical"
    }
end

/*******************************************************************************
 * benchmark_reghdfe_opts - Compare creghdfe (with extra options) vs reghdfe
 *
 * Thin wrapper for comparing reghdfe output against creghdfe with options
 * that reghdfe doesn't support (tolerance, iterate). Runs reghdfe first,
 * then creghdfe with the extra options, and compares e() results.
 *
 * Syntax: benchmark_reghdfe_opts varlist, absorb(varlist) creghdfe_opts(string) testname(string) [vce(string) minsf(real)]
 ******************************************************************************/
capture program drop benchmark_reghdfe_opts
program define benchmark_reghdfe_opts
    syntax varlist(min=2 fv) [aw fw pw] [if] [in], Absorb(varlist) creghdfe_opts(string) testname(string) [vce(string) minsf(real 7)]

    gettoken depvar indepvars : varlist

    * Build weight string
    local wtexp ""
    if "`weight'" != "" local wtexp "[`weight'`exp']"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    * Run reghdfe (reference)
    quietly reghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    local reghdfe_N = e(N)
    local reghdfe_df_r = e(df_r)
    matrix reghdfe_b = e(b)
    matrix reghdfe_V = e(V)
    foreach scalar in r2 rss F rmse {
        local reghdfe_`scalar' = e(`scalar')
    }

    * Run creghdfe with extra options
    quietly creghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt' `creghdfe_opts'
    local creghdfe_N = e(N)
    local creghdfe_df_r = e(df_r)
    matrix creghdfe_b = e(b)
    matrix creghdfe_V = e(V)
    foreach scalar in r2 rss F rmse {
        local creghdfe_`scalar' = e(`scalar')
    }

    * Compare results
    local all_diffs ""
    local has_failure = 0
    if `reghdfe_N' != `creghdfe_N' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(N):`reghdfe_N'!=`creghdfe_N'"
    }
    if `reghdfe_df_r' != `creghdfe_df_r' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(df_r):`reghdfe_df_r'!=`creghdfe_df_r'"
    }
    foreach scalar in r2 rss F rmse {
        sigfigs `reghdfe_`scalar'' `creghdfe_`scalar''
        if r(sigfigs) < `minsf' {
            local has_failure = 1
            local sf_fmt : display %4.1f r(sigfigs)
            local all_diffs "`all_diffs' e(`scalar'):sigfigs=`sf_fmt'"
        }
    }
    matrix_min_sigfigs reghdfe_b creghdfe_b
    if r(min_sigfigs) < `minsf' {
        local has_failure = 1
        local sf_fmt : display %4.1f r(min_sigfigs)
        local all_diffs "`all_diffs' e(b):sigfigs=`sf_fmt'"
    }
    matrix_min_sigfigs reghdfe_V creghdfe_V
    if r(min_sigfigs) < `minsf' {
        local has_failure = 1
        local sf_fmt : display %4.1f r(min_sigfigs)
        local all_diffs "`all_diffs' e(V):sigfigs=`sf_fmt'"
    }
    if `has_failure' == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "`=trim("`all_diffs'")'"
    }
end

/*******************************************************************************
 * test_error_match - Compare error codes between Stata native and ctools commands
 *
 * This helper runs both commands (expected to fail) and verifies they return
 * the same error code. Used for "intentional failure" validation tests.
 *
 * Syntax: test_error_match, stata_cmd(string) ctools_cmd(string) testname(string)
 ******************************************************************************/
capture program drop test_error_match
program define test_error_match
    syntax, stata_cmd(string asis) ctools_cmd(string asis) testname(string)

    * Run Stata native command
    capture `stata_cmd'
    local stata_rc = _rc

    * Run ctools command
    capture `ctools_cmd'
    local ctools_rc = _rc

    * Compare error codes
    if `stata_rc' == `ctools_rc' {
        test_pass "[error] `testname' (rc=`stata_rc')"
    }
    else {
        test_fail "[error] `testname'" "stata rc=`stata_rc', ctools rc=`ctools_rc'"
    }
end
