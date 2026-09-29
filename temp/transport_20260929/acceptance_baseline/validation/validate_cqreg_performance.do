* Small end-to-end checks for cqreg performance changes (at most 603 obs).
* Optional argument selects an isolated plugin/ado build directory.
args builddir
if `"`builddir'"' == "" local builddir "build"
adopath ++ "`builddir'"
clear all
set more off
set seed 81273
set obs 603
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double y = 2*x1 - 3*x2 + 4 + rnormal()
generate long cluster = ceil(_n/9)
generate byte keep = mod(_n,7) != 0

foreach q in .1 .5 .9 {
    quietly qreg y x1 x2 if keep, quantile(`q')
    matrix expected_b = e(b)
    scalar expected_obj = e(sum_adev)
    scalar expected_q = e(q_v)
    foreach density in fitted residual kernel {
        quietly cqreg y x1 x2 if keep, quantile(`q') denmethod(`density')
        assert e(convcode) == 0
        assert e(N) == 517
        assert mreldif(e(b),expected_b) < 1e-6
        assert abs(e(sum_adev)-expected_obj) < 1e-6
        assert abs(e(q_v)-expected_q) < 1e-12*(1+abs(expected_q))
    }
}
foreach vce in robust "cluster cluster" {
    quietly cqreg y x1 x2, vce(`vce')
    assert e(convcode) == 0
    assert e(N) == 603
    matrix V = e(V)
    forvalues j=1/3 {
        assert V[`j',`j'] > 0 & V[`j',`j'] < .
    }
}
* Repeated responses and residuals exercise tie-aware order statistics.
replace y = mod(_n,3)-1
replace x1 = ceil(_n/3)/100
quietly cqreg y x1, denmethod(residual)
assert abs(_b[_cons]) < 1e-8
assert e(convcode) == 0
capture cqreg y x1 x2, maxiter(1)
assert _rc == 430
display "cqreg performance smoke checks passed"
exit, clear
