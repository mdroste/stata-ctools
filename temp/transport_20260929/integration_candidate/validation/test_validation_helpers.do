/* Mutation checks for the validation harness. Called first by audit_p1, after
   validate_setup.do. Each bad implementation must record exactly one failure.
   Drop the stubs before returning so subsequent tests load the real ado files. */

capture program drop helper_expect
program define helper_expect
    args failures passes
    assert $TESTS_FAILED == `failures'
    assert $TESTS_PASSED == `passes'
    assert $TESTS_TOTAL == `failures' + `passes'
    reset_counters
end

assert_scalar_equal 1e-16 9e-16 7 "tiny values must retain relative precision"
helper_expect 1 0
assert_scalar_equal 1e-16 1e-16 7 "identical tiny values"
helper_expect 0 1
sigfigs 1e-16 1.000001e-16
assert abs(r(sigfigs) - 6) < 1e-5
matrix helper_a = (1e-16, 1)
matrix helper_b = (9e-16, 1)
assert_matrix_equal helper_a helper_b 7 "tiny matrix element differs"
helper_expect 1 0

/* Specialized cancellation checks must tolerate only machine roundoff.
   A small covariance in small units is still tested on its own scale. */
sigfigs_cancellation 0 1e-16 1
assert r(sigfigs) == 15
sigfigs_cancellation 0 1e-12 1
assert r(sigfigs) < 7
sigfigs_cancellation 1e-24 9e-24 1e-20
assert r(sigfigs) < 7
sigfigs_cancellation 0 . 1
assert r(sigfigs) < 7
sigfigs_cancellation 0 1e-16 0
assert r(sigfigs) < 7

capture program drop csort
program define csort, rclass
    syntax varlist [, ALGorithm(string) STReam(integer 0) VERBOSE]
    assert `stream' == 2
    if "$HELPER_MUTATION" == "identity" {
        return scalar stream = 1
        exit
    }
    sort `varlist', stable
    if "$HELPER_MUTATION" == "corrupt" replace witness = -witness in 1
    if "$HELPER_MUTATION" == "dropvar" drop witness
    return scalar stream = ("$HELPER_MUTATION" != "ignorestream")
end
clear
set obs 3
gen key = 4 - _n
gen witness = _n
foreach mutation in identity corrupt dropvar ignorestream correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_sort key, stream(2) testname("sort mutation: `mutation'")
    if "`mutation'" == "correct" helper_expect 0 1
    else helper_expect 1 0
}
program drop csort

capture program drop cimport
program define cimport
    gettoken format rest : 0
    import `format' `rest'
    if "$HELPER_MUTATION" == "reorder" gsort -id
    if "$HELPER_MUTATION" == "corrupt" replace value = -value in 1
end
tempfile csv
tempname fh
file open `fh' using "`csv'", write text replace
file write `fh' "id,value" _n "1,100" _n "2,200" _n "3,300" _n
file close `fh'
foreach mutation in reorder corrupt correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_import using "`csv'", testname("import mutation: `mutation'")
    if "`mutation'" == "correct" helper_expect 0 1
    else helper_expect 1 0
}
program drop cimport

/* Use identical coefficients, covariance, and samples, varying only the one
   result being tested. A missing scalar or a tiny coefficient must not escape
   comparison just because all the other output agrees. */
capture program drop helper_ppml_result
program define helper_ppml_result, eclass
    tempname b V
    tempvar sample
    matrix `b' = (1e-16, 2)
    matrix colnames `b' = x _cons
    matrix `V' = (1, 1e-16 \ 1e-16, 1)
    matrix rownames `V' = x _cons
    matrix colnames `V' = x _cons
    if "$HELPER_MUTATION" == "coefficient" matrix `b'[1,1] = 9e-16
    if "$HELPER_MUTATION" == "covariance" {
        matrix `V'[1,2] = 9e-16
        matrix `V'[2,1] = 9e-16
    }
    if "$HELPER_MUTATION" == "covariance_above_roundoff" {
        matrix `V'[1,2] = 1e-12
        matrix `V'[2,1] = 1e-12
    }
    gen byte `sample' = 1
    ereturn post `b' `V', esample(`sample') obs(3)
    ereturn scalar ll = -10
    ereturn scalar r2_p = 0
    ereturn scalar N_clust = 3
    if "$HELPER_MUTATION" == "ll" ereturn scalar ll = .
    if "$HELPER_MUTATION" == "r2_p" ereturn scalar r2_p = .
    if "$HELPER_MUTATION" == "N_clust" ereturn scalar N_clust = .
    if "$HELPER_MUTATION" == "zero" ereturn scalar r2_p = .5
    if "$HELPER_MUTATION" == "r2_small" ereturn scalar r2_p = 1e-12
    if "$HELPER_MUTATION" == "r2_roundoff" ereturn scalar r2_p = -2e-15
end
capture program drop ppmlhdfe
program define ppmlhdfe, eclass
    local mutation "$HELPER_MUTATION"
    global HELPER_MUTATION "correct"
    helper_ppml_result
    global HELPER_MUTATION "`mutation'"
end
capture program drop cpplmhdfe
program define cpplmhdfe, eclass
    helper_ppml_result
end
clear
set obs 3
gen x = _n
gen y = _n
gen g = _n
foreach mutation in ll r2_p N_clust zero r2_small r2_roundoff coefficient covariance correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_ppmlhdfe y x, absorb(g) testname("PPML mutation: `mutation'")
    if inlist("`mutation'", "correct", "r2_roundoff") helper_expect 0 1
    else helper_expect 1 0
}
/* Exercise the quantile helper itself: covariance roundoff may pass, but
   a 1e-12 covariance error or a tiny coefficient error must still fail. */
capture program drop qreg
program define qreg, eclass
    ppmlhdfe `0'
end
capture program drop cqreg
program define cqreg, eclass
    cpplmhdfe `0'
end
foreach mutation in covariance covariance_above_roundoff coefficient correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_qreg y x, testname("quantile mutation: `mutation'")
    if inlist("`mutation'", "covariance", "correct") helper_expect 0 1
    else helper_expect 1 0
}
program drop qreg cqreg ppmlhdfe cpplmhdfe helper_ppml_result

/* A failing reference must not be turned into a nonmissing-value smoke test. */
capture program drop crangestat
program define crangestat
    if strpos("`0'", "__brp_ctools") gen __brp_ctools = 1
    else gen __brs_ctools = 1
end
capture program drop rangestat
program define rangestat
    exit 198
end
benchmark_rangestat mean x, interval(g -1 1) testname("failed range reference")
helper_expect 1 0
benchmark_range_pctile p25 x, interval(g -1 1) testname("incorrect window percentiles")
helper_expect 1 0
program drop crangestat rangestat helper_expect
macro drop HELPER_MUTATION
matrix drop helper_a helper_b
