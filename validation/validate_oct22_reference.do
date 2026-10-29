* Additional parity against installed reference ado commands.
* Install ivreg2/ranktest, binscatter and psmatch2 before running.
version 16
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
local reference : environment CTOOLS_TEST_REFERENCE_DIR
if `"`reference'"' != "" adopath ++ `"`reference'"'
set obs 100
set seed 2345
gen double z=rnormal()
gen double z2=rnormal()
gen double x=z+.5*z2+rnormal()
gen double y=10+.2*x+rnormal()
gen cl=ceil(_n/5)
foreach options in "" "small" "robust" "robust small" "cluster(cl)" "cluster(cl) small" ///
    "dofminus(2)" "small dofminus(2) sdofminus(3)" "robust dofminus(2) sdofminus(3)" ///
    "small robust dofminus(2) sdofminus(3)" "gmm2s robust" "gmm2s robust small" ///
    "gmm2s robust dofminus(2) sdofminus(3)" "gmm2s robust small dofminus(2) sdofminus(3)" ///
    "noconstant small" "cluster(cl) small dofminus(2) sdofminus(3)" {
    quietly ivreg2 y (x=z z2), `options'
    matrix expected=e(V)
    scalar expected_F=e(F)
    scalar expected_rmse=e(rmse)
    scalar expected_r2a=e(r2_a)
    quietly civreghdfe y (x=z z2), `options'
    matrix actual=e(V)
    assert mreldif(actual,expected)<1e-8
    assert reldif(e(F),expected_F)<1e-8
    assert reldif(e(rmse),expected_rmse)<1e-8
    assert reldif(e(r2_a),expected_r2a)<1e-8
}
* Effective sample size and partialled controls must retain the same corrections.
gen fw=1+mod(_n,3)
gen double aw=.2+runiform()
foreach weight in "" "[fw=fw]" "[aw=aw]" {
    foreach options in "partial(z2)" "partial(z2) small" ///
        "partial(z2) robust dofminus(2) sdofminus(3)" ///
        "partial(z2) cluster(cl) small dofminus(2) sdofminus(3)" {
        quietly ivreg2 y z2 (x=z) `weight', `options'
        matrix expected=e(V)
        scalar expected_F=e(F)
        scalar expected_rmse=e(rmse)
        quietly civreghdfe y z2 (x=z) `weight', `options'
        matrix actual=e(V)
        * Compare the common slope: ivreg2 partial() omits the constant from e(V).
        matrix actual=actual[1,1]
        assert mreldif(actual,expected)<1e-8
        assert reldif(e(F),expected_F)<1e-8
        assert reldif(e(rmse),expected_rmse)<1e-8
    }
}
* Covariance rank deficiency does not change ivreg2's slope-count divisor.
gen two_clusters=mod(_n,2)
gen double large_scale=1e6*rnormal()
gen double small_scale=1e-6*rnormal()
foreach rhs in "z2" "large_scale small_scale" {
    quietly ivreg2 y `rhs' (x=z), cluster(two_clusters) small
    scalar expected_F=e(F)
    scalar expected_df=e(df_m)
    quietly civreghdfe y `rhs' (x=z), cluster(two_clusters) small
    assert e(df_m)==expected_df
    assert reldif(e(F),expected_F)<1e-7
}
* Singleton fit requirements differ from retaining singleton scatter points.
clear
set obs 8
gen x=mod(_n,2)
gen y=_n
gen g=(_n==8)
foreach line in none lfit qfit {
    capture quietly binscatter y x, by(g) discrete linetype(`line') nodraw
    local expected_rc=_rc
    capture quietly cbinscatter y x, by(g) discrete linetype(`line') nograph
    assert _rc==`expected_rc'
}
* Reference treatment and option semantics, not assumptions about valid inputs.
foreach options in "neighbor(0)" "neighbor(-1)" "kernel bwidth(0)" "caliper(-1)" "kerneltype(unknown)" {
    clear
    set obs 10
    gen t=mod(_n,2)
    gen p=_n/11
    gen y=_n
    capture quietly psmatch2 t, pscore(p) outcome(y) `options'
    local expected_rc=_rc
    scalar expected_att=r(att)
    capture quietly cpsmatch t, pscore(p) outcome(y) `options'
    assert _rc==`expected_rc'
    if !`expected_rc' assert r(att)==expected_att
}
noi di "OCT22_REFERENCE_COMPLETE"
