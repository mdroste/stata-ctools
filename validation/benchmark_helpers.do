/*******************************************************************************
 * benchmark_helpers.do
 *
 * Benchmark helper programs for ctools validation test suite
 * Each program runs both ctools and Stata native commands, compares results,
 * and reports [PASS] or [FAIL] with detailed error descriptions
 ******************************************************************************/

/*******************************************************************************
 * benchmark_sort - Compare csort vs sort, stable
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

    * Run Stata sort, stable (quietly)
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

    * Compare datasets
    capture quietly cf _all using `stata_sorted'
    local cf_rc = _rc

    restore

    if `cf_rc' == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "sorted data differs from Stata sort"
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

    if `counts_match' & `cf_rc' == 0 {
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
 * Syntax: benchmark_reghdfe depvar indepvars [weight], absorb(varlist) [options]
 ******************************************************************************/
capture program drop benchmark_reghdfe
program define benchmark_reghdfe
    syntax varlist(min=2 fv) [aw fw pw] [if] [in], Absorb(varlist) [vce(string) testname(string) tol(real 1e-6)]

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

    * Compare continuous scalars (use tolerance)
    foreach scalar in r2 r2_a r2_within rss tss mss F rmse {
        local val1 = `reghdfe_`scalar''
        local val2 = `creghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            local diff = abs(`val1' - `val2')
            * Use relative tolerance for larger values
            if abs(`val1') > 1 {
                local reldiff = `diff' / abs(`val1')
                if `reldiff' > `tol' {
                    local has_failure = 1
                    local all_diffs "`all_diffs' e(`scalar'):reldiff=`reldiff'"
                }
            }
            else {
                if `diff' > `tol' {
                    local has_failure = 1
                    local all_diffs "`all_diffs' e(`scalar'):diff=`diff'"
                }
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local all_diffs "`all_diffs' e(`scalar'):missing_in_creghdfe"
        }
    }

    * Compare coefficients e(b)
    tempname diff_b
    matrix `diff_b' = reghdfe_b - creghdfe_b
    local cols = colsof(`diff_b')
    local maxdiff_b = 0
    local maxdiff_b_idx = 0
    forvalues j = 1/`cols' {
        local d = abs(`diff_b'[1, `j'])
        if `d' > `maxdiff_b' {
            local maxdiff_b = `d'
            local maxdiff_b_idx = `j'
        }
    }

    if `maxdiff_b' >= `tol' {
        local has_failure = 1
        local coefname : colnames reghdfe_b
        local badcoef : word `maxdiff_b_idx' of `coefname'
        local all_diffs "`all_diffs' e(b)[`badcoef']:maxdiff=`maxdiff_b'"
    }

    * Compare VCE e(V) using relative tolerance
    tempname diff_V
    matrix `diff_V' = reghdfe_V - creghdfe_V
    local rows = rowsof(`diff_V')
    local cols = colsof(`diff_V')
    local maxreldiff_V = 0
    local maxreldiff_V_i = 0
    local maxreldiff_V_j = 0
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local diff = abs(`diff_V'[`i', `j'])
            local base = abs(reghdfe_V[`i', `j'])
            if `base' > 1e-10 {
                local reldiff = `diff' / `base'
            }
            else {
                local reldiff = `diff'
            }
            if `reldiff' > `maxreldiff_V' {
                local maxreldiff_V = `reldiff'
                local maxreldiff_V_i = `i'
                local maxreldiff_V_j = `j'
            }
        }
    }

    if `maxreldiff_V' >= `tol' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(V)[`maxreldiff_V_i',`maxreldiff_V_j']:maxreldiff=`maxreldiff_V'"
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
 * Syntax: benchmark_qreg depvar indepvars [, quantile(real) options]
 ******************************************************************************/
capture program drop benchmark_qreg
program define benchmark_qreg
    syntax varlist(min=2 fv) [if] [in], [Quantile(real 0.5) vce(string) testname(string) tol(real 1e-6)]

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

    * Compare continuous scalars (use tolerance)
    foreach scalar in q sum_adev sum_rdev q_v f_r {
        local val1 = `qreg_`scalar''
        local val2 = `cqreg_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            local diff = abs(`val1' - `val2')
            * Use relative tolerance for larger values
            if abs(`val1') > 1 {
                local reldiff = `diff' / abs(`val1')
                if `reldiff' > `tol' {
                    local has_failure = 1
                    local all_diffs "`all_diffs' e(`scalar'):reldiff=`reldiff'"
                }
            }
            else {
                if `diff' > `tol' {
                    local has_failure = 1
                    local all_diffs "`all_diffs' e(`scalar'):diff=`diff'"
                }
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local all_diffs "`all_diffs' e(`scalar'):missing_in_cqreg"
        }
    }

    * Compare coefficients e(b)
    tempname diff_b
    matrix `diff_b' = qreg_b - cqreg_b
    local cols = colsof(`diff_b')
    local maxdiff_b = 0
    local maxdiff_b_idx = 0
    forvalues j = 1/`cols' {
        local d = abs(`diff_b'[1, `j'])
        if `d' > `maxdiff_b' {
            local maxdiff_b = `d'
            local maxdiff_b_idx = `j'
        }
    }

    if `maxdiff_b' >= `tol' {
        local has_failure = 1
        local coefname : colnames qreg_b
        local badcoef : word `maxdiff_b_idx' of `coefname'
        local all_diffs "`all_diffs' e(b)[`badcoef']:maxdiff=`maxdiff_b'"
    }

    * Compare VCE e(V) using relative tolerance
    tempname diff_V
    matrix `diff_V' = qreg_V - cqreg_V
    local rows = rowsof(`diff_V')
    local cols = colsof(`diff_V')
    local maxreldiff_V = 0
    local maxreldiff_V_i = 0
    local maxreldiff_V_j = 0
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local diff = abs(`diff_V'[`i', `j'])
            local base = abs(qreg_V[`i', `j'])
            if `base' > 1e-10 {
                local reldiff = `diff' / `base'
            }
            else {
                local reldiff = `diff'
            }
            if `reldiff' > `maxreldiff_V' {
                local maxreldiff_V = `reldiff'
                local maxreldiff_V_i = `i'
                local maxreldiff_V_j = `j'
            }
        }
    }

    if `maxreldiff_V' >= `tol' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(V)[`maxreldiff_V_i',`maxreldiff_V_j']:maxreldiff=`maxreldiff_V'"
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
 * Syntax: benchmark_import using filename [, options]
 ******************************************************************************/
capture program drop benchmark_import
program define benchmark_import
    syntax using/, [DELIMiters(string) VARNames(string) CASE(string) testname(string)]

    if "`testname'" == "" local testname "import `using'"

    * Build option string
    local opts "clear"
    if "`delimiters'" != "" local opts "`opts' delimiters(`delimiters')"
    if "`varnames'" != "" local opts "`opts' varnames(`varnames')"
    if "`case'" != "" local opts "`opts' case(`case')"

    preserve

    * Run Stata import (quietly)
    capture quietly import delimited `using', `opts'
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
    capture quietly cimport delimited `using', `opts'
    local cimport_rc = _rc

    if `cimport_rc' != 0 {
        restore
        test_fail "`testname'" "cimport returned error `cimport_rc'"
        exit
    }

    local cimport_n = _N
    local cimport_k = c(k)

    restore

    * Compare dimensions
    if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "dimensions differ: stata(N=`stata_n',K=`stata_k') vs cimport(N=`cimport_n',K=`cimport_k')"
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
 * Syntax: benchmark_ivreghdfe spec, absorb(varlist) [options]
 ******************************************************************************/
capture program drop benchmark_ivreghdfe
program define benchmark_ivreghdfe
    syntax anything(name=spec), Absorb(varlist) [vce(string) testname(string) tol(real 1e-6)]

    if "`testname'" == "" local testname "ivreghdfe `spec', absorb(`absorb')"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run ivreghdfe (quietly)
    capture quietly ivreghdfe `spec', absorb(`absorb') `vceopt'
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
    capture quietly civreghdfe `spec', absorb(`absorb') `vceopt'
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

    * Compare continuous scalars (use tolerance)
    foreach scalar in r2 r2_a rss mss F rmse idstat idp widstat sargan sarganp {
        local val1 = `ivreghdfe_`scalar''
        local val2 = `civreghdfe_`scalar''
        * Only compare if both are non-missing
        if !missing(`val1') & !missing(`val2') {
            local diff = abs(`val1' - `val2')
            * Use relative tolerance for larger values
            if abs(`val1') > 1 {
                local reldiff = `diff' / abs(`val1')
                if `reldiff' > `tol' {
                    local has_failure = 1
                    local all_diffs "`all_diffs' e(`scalar'):reldiff=`reldiff'"
                }
            }
            else {
                if `diff' > `tol' {
                    local has_failure = 1
                    local all_diffs "`all_diffs' e(`scalar'):diff=`diff'"
                }
            }
        }
        else if !missing(`val1') & missing(`val2') {
            local all_diffs "`all_diffs' e(`scalar'):missing_in_civreghdfe"
        }
    }

    * Compare coefficients e(b)
    tempname diff_b
    matrix `diff_b' = ivreghdfe_b - civreghdfe_b
    local cols = colsof(`diff_b')
    local maxdiff_b = 0
    local maxdiff_b_idx = 0
    forvalues j = 1/`cols' {
        local d = abs(`diff_b'[1, `j'])
        if `d' > `maxdiff_b' {
            local maxdiff_b = `d'
            local maxdiff_b_idx = `j'
        }
    }

    if `maxdiff_b' >= `tol' {
        local has_failure = 1
        local coefname : colnames ivreghdfe_b
        local badcoef : word `maxdiff_b_idx' of `coefname'
        local all_diffs "`all_diffs' e(b)[`badcoef']:maxdiff=`maxdiff_b'"
    }

    * Compare VCE e(V) using relative tolerance
    tempname diff_V
    matrix `diff_V' = ivreghdfe_V - civreghdfe_V
    local rows = rowsof(`diff_V')
    local cols = colsof(`diff_V')
    local maxreldiff_V = 0
    local maxreldiff_V_i = 0
    local maxreldiff_V_j = 0
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local diff = abs(`diff_V'[`i', `j'])
            local base = abs(ivreghdfe_V[`i', `j'])
            if `base' > 1e-10 {
                local reldiff = `diff' / `base'
            }
            else {
                local reldiff = `diff'
            }
            if `reldiff' > `maxreldiff_V' {
                local maxreldiff_V = `reldiff'
                local maxreldiff_V_i = `i'
                local maxreldiff_V_j = `j'
            }
        }
    }

    if `maxreldiff_V' >= `tol' {
        local has_failure = 1
        local all_diffs "`all_diffs' e(V)[`maxreldiff_V_i',`maxreldiff_V_j']:maxreldiff=`maxreldiff_V'"
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
