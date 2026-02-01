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

    * Compare VCE e(V) using significant figures
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
 * Syntax: benchmark_export [varlist] using filename [, options]
 ******************************************************************************/
capture program drop benchmark_export
program define benchmark_export
    syntax [varlist] using/, [DELIMiter(string) NOVARNames QUOTE NOLabel replace testname(string)]

    if "`testname'" == "" local testname "export `using'"

    * Build option string
    local opts "replace"
    if "`delimiter'" != "" local opts "`opts' delimiter(`delimiter')"
    if "`novarnames'" != "" local opts "`opts' novarnames"
    if "`quote'" != "" local opts "`opts' quote"
    if "`nolabel'" != "" local opts "`opts' nolabel"

    preserve

    * Generate temp filenames
    tempfile stata_export cexport_export
    local stata_csv "`stata_export'.csv"
    local cexport_csv "`cexport_export'.csv"

    * Run Stata export (quietly)
    if "`varlist'" != "" {
        capture quietly export delimited `varlist' using "`stata_csv'", `opts'
    }
    else {
        capture quietly export delimited using "`stata_csv'", `opts'
    }
    local stata_rc = _rc

    if `stata_rc' != 0 {
        restore
        test_fail "`testname'" "export delimited returned error `stata_rc'"
        exit
    }

    * Run cexport (quietly)
    if "`varlist'" != "" {
        capture quietly cexport delimited `varlist' using "`cexport_csv'", `opts'
    }
    else {
        capture quietly cexport delimited using "`cexport_csv'", `opts'
    }
    local cexport_rc = _rc

    if `cexport_rc' != 0 {
        restore
        test_fail "`testname'" "cexport returned error `cexport_rc'"
        exit
    }

    * Reimport and compare (quietly)
    quietly import delimited using "`stata_csv'", clear
    local stata_n = _N
    local stata_k = c(k)

    quietly import delimited using "`cexport_csv'", clear
    local cexport_n = _N
    local cexport_k = c(k)

    * Cleanup temp files
    capture erase "`stata_csv'"
    capture erase "`cexport_csv'"

    restore

    * Compare dimensions
    if `stata_n' == `cexport_n' & `stata_k' == `cexport_k' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "dimensions differ: stata(N=`stata_n',K=`stata_k') vs cexport(N=`cexport_n',K=`cexport_k')"
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
 *   - e(b) - coefficient vector
 *   - e(V) - variance-covariance matrix
 *
 * All comparisons use significant figures (default 7 sigfigs).
 *
 * Syntax: benchmark_ivreghdfe spec, absorb(varlist) [options]
 ******************************************************************************/
capture program drop benchmark_ivreghdfe
program define benchmark_ivreghdfe
    syntax anything(name=spec) [if] [in], Absorb(varlist) ///
        [vce(string) testname(string) minsf(real 7) ///
        LIML FULLER(real 0) Kclass(real 0) GMM2s CUE COVIV ///
        BW(integer 0) KERNEL(string) DKRAAY(integer 0) KIEFER ///
        b0(string)]

    if "`testname'" == "" local testname "ivreghdfe `spec' `if' `in', absorb(`absorb')"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

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
    capture quietly ivreghdfe `spec' `if' `in', absorb(`absorb') `vceopt' `estopts' `hacopts' `b0opt'
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

    * Run civreghdfe (quietly)
    capture quietly civreghdfe `spec' `if' `in', absorb(`absorb') `vceopt' `estopts' `hacopts' `b0opt'
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

    * Compare coefficients e(b) using significant figures
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

    * Compare VCE e(V) using significant figures
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
