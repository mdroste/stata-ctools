* Analytical coefficient oracles: default reghdfe also fails the first case.
* Called noisily because the diagnostic itself is part of the contract.
_ct_loud creghdfe_numerics
capture program drop _ct_creg_weak
program define _ct_creg_weak
    clear
    set obs 6
    gen byte a = cond(inrange(_n,3,4),2,1)
    gen byte b = cond(_n<=4,2,1)
    gen byte u = 2*mod(_n-1,2)-1
    gen double w = cond(a==b, 2, 1e-14)
    gen double x = u + 1e6*(a-b)
    gen double y = 2*x + 1e8*(4*a-5*b) + (a==b)*(2*a-3)*u
    * The error is orthogonal to x and every FE; beta_x is exactly 2.
    tempfile diagnostic
    tempname audit fh
    log using "`diagnostic'", name(`audit') text replace
    capture noisily creghdfe y x [aw=w], absorb(a b)
    local rc = _rc
    log close `audit'
    assert `rc' == 498
    assert "`e(cmd)'" == ""
    file open `fh' using "`diagnostic'", read text
    local found = 0
    file read `fh' line
    while r(eof)==0 {
        if strpos(`"`line'"', "numerical breakdown in fixed-effect projection") local found = 1
        file read `fh' line
    }
    file close `fh'
    assert `found'
    * Changing units or requested tolerance cannot turn breakdown into success.
    foreach power in -20 0 20 {
        gen double xs = x*2^`power'
        gen double ys = y*2^`power'
        foreach opt in "" "tolerance(1e-14)" "quad" "nostandardize" "abort(0)" {
            capture creghdfe ys xs [aw=w], absorb(a b) `opt'
            local rc = _rc
            if `rc' == 0 {
                assert abs(_b[xs]-2)<1e-7 & e(converged)==1
            }
            else {
                assert `rc' == 498
                assert "`e(cmd)'" == ""
            }
        }
        drop xs ys
    }
    * Removing only absorbed components preserves the slope and must succeed.
    gen double xc = x-1e6*(a-b)
    gen double yc = y-2e6*(a-b)-1e8*(4*a-5*b)
    quietly creghdfe yc xc [aw=w], absorb(a b)
    assert abs(_b[xc]-2)<1e-10 & e(converged)==1
end
capture noisily _ct_creg_weak
if _rc test_fail "creghdfe weak connection: accurate or clear numerical rejection" "r(`=_rc')"
else test_pass "creghdfe weak connection: accurate or clear numerical rejection"
program drop _ct_creg_weak

capture program drop _ct_creg_scaled
program define _ct_creg_scaled
    clear
    set obs 8
    gen byte x = 2*mod(_n-1,2)-1
    gen byte v = 2*mod(floor((_n-1)/2),2)-1
    gen byte z = 2*(_n>4)-1
    gen byte g = 1+(x*z>0)
    gen double q = x+v/1024
    gen double y = 2*x+3*q+2^34*z+x*v+g
    * x*v is orthogonal to all regressors and the FE. The inputs are exact.
    foreach opt in "" "quad" "nostandardize" "tolerance(1e-14)" "threads(1)" "threads(4)" {
        quietly creghdfe y x q z, absorb(g) `opt'
        assert e(converged)==1
        assert abs(_b[x]-2)<1e-7 & abs(_b[q]-3)<1e-7
        assert reldif(_b[z],2^34)<1e-14
        assert abs(e(rss)-8)<1e-5
    }
    gen double yc = y-2^34*z
    assert yc == 2*x+3*q+x*v+g
    quietly creghdfe yc x q z, absorb(g)
    assert abs(_b[x]-2)<1e-7 & abs(_b[q]-3)<1e-7
    * A nearly absorbed column can also break the unweighted recurrence.
    clear
    set seed 251125
    set obs 2000
    gen long a = 1+mod(_n-1,200)
    gen long b = runiformint(1,80)
    gen long c = runiformint(1,30)
    gen double z = rnormal()
    gen double x = .5*z+rnormal()
    gen double y = 2*x-.7*z+a/200-b/80+rnormal()
    gen double near = a+.01*x
    quietly creghdfe y x z, absorb(a b) tolerance(1e-12)
    local expected = _b[x]/.01
    capture creghdfe y near z, absorb(a b) tolerance(1e-12)
    local rc = _rc
    if `rc'==0 assert reldif(_b[near],`expected')<1e-7 & e(converged)==1
    else assert `rc'==498 & "`e(cmd)'"==""
    * A nearby well-behaved column must still succeed.
    replace near = a+.001*x
    quietly creghdfe y near z, absorb(a b) tolerance(1e-12)
    assert reldif(_b[near],10*`expected')<1e-7 & e(converged)==1
end
capture noisily _ct_creg_scaled
if _rc test_fail "creghdfe near-collinear coefficients with extreme response scaling" "r(`=_rc')"
else test_pass "creghdfe near-collinear coefficients with extreme response scaling"
program drop _ct_creg_scaled

capture program drop _ct_creg_technique
program define _ct_creg_technique
    sysuse auto, clear
    foreach technique in lsmr lsqr {
        tempfile diagnostic
        tempname audit fh
        log using "`diagnostic'", name(`audit') text replace
        capture noisily creghdfe price weight, absorb(turn foreign) technique(`technique')
        local rc = _rc
        log close `audit'
        assert `rc' == 198 & "`e(cmd)'" == ""
        file open `fh' using "`diagnostic'", read text
        local found = 0
        file read `fh' line
        while r(eof)==0 {
            if strpos(`"`line'"', "technique(`technique') is not implemented for this model") local found = 1
            file read `fh' line
        }
        file close `fh'
        assert `found'
    }
    quietly creghdfe price weight, absorb(turn foreign) technique(map)
    assert e(converged)==1
end
capture noisily _ct_creg_technique
if _rc test_fail "creghdfe refuses unsupported Krylov solvers explicitly" "r(`=_rc')"
else test_pass "creghdfe refuses unsupported Krylov solvers explicitly"
program drop _ct_creg_technique
