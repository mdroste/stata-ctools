/* Memory-audit integration contracts; native fault injection lives in
   test_memory_safety_native.py. Run through the machine's oldstata wrapper. */
do validation/validate_setup.do

capture noisily {
    clear
    set obs 120
    gen double x = _n + sin(_n)
    gen double y = 3 + 2*x + cos(_n)
    gen byte g = mod(_n, 3)
    gen double large_g = cond(g==0, -2^63, cond(g==1, 2^62, 1e30))
    quietly creghdfe y x, absorb(g)
    scalar slope = _b[x]
    quietly creghdfe y x, absorb(large_g)
    assert reldif(_b[x], scalar(slope)) < 1e-10
    replace large_g = g + .25
    quietly creghdfe y x, absorb(large_g)
    assert reldif(_b[x], scalar(slope)) < 1e-10

    gen double w = 300000000
    foreach vce in "" "vce(robust)" "vce(cluster g)" {
        quietly creghdfe y x [fw=w], absorb(g) `vce'
        assert e(N)==36000000000
        assert reldif(_b[x], scalar(slope)) < 1e-10
        assert _se[x]>0 & _se[x]<.
    }
}
if _rc test_fail "large FE labels and frequency weights" "r(`=_rc')"
else test_pass "large FE labels and frequency weights"

capture noisily {
    clear
    set obs 512
    gen long key = _n
    gen byte g = floor((_n-1)/64)
    gen double x = mod(_n, 10)
    crangestat (min) lo=x (max) hi=x (median) med=x (mean) avg=x, ///
        interval(key 0 0) by(g) threads(8)
    assert lo==x & hi==x & med==x & avg==x
    drop lo hi med avg
    crangestat (median) med=x (mean) avg=x, interval(key 0 0) threads(8)
    assert med==x & avg==x
    capture crangestat (mean) invalid=x, interval(key 0 0) threads(2147483647)
    assert _rc==198
}
if _rc test_fail "range scratch paths and thread bound" "r(`=_rc')"
else test_pass "range scratch paths and thread bound"

capture noisily {
    tempfile original selected
    save `original'
    set seed 142
    csample, count(20) by(g) threads(1)
    assert _N==160
    sort key
    save `selected'
    use `original', clear
    set seed 142
    csample, count(20) by(g) threads(8)
    sort key
    cf _all using `selected'
    use `original', clear
    csample if key<=100, count(20) threads(8)
    count if key<=100
    assert r(N)==20
    use `original', clear
    gen double freq=0
    cbsample 20 if key<=100, weight(freq) threads(8)
    quietly summarize freq, meanonly
    assert r(sum)==20
}
if _rc test_fail "sampling grouping and filtered maps" "r(`=_rc')"
else test_pass "sampling grouping and filtered maps"

capture noisily {
    clear
    set obs 1000
    set seed 324
    gen double x = rnormal()
    gen double z = x + rnormal()
    gen double z2 = 2*z
    gen double y = 2*x + 3*z + rnormal()
    gen g = mod(_n,10)
    foreach adjustment in "" "absorb(g)" {
        cbinscatter y x, controls(z) method(binsreg) nograph `adjustment'
        matrix reference=e(bindata)
        cbinscatter y x, controls(z z2) method(binsreg) nograph `adjustment'
        mata: assert(max(abs(st_matrix("reference")-st_matrix("e(bindata)"))) < 1e-7)
    }
}
if _rc test_fail "cbinscatter controls and FE offsets" "r(`=_rc')"
else test_pass "cbinscatter controls and FE offsets"

print_summary "memory safety"
global CTOOLS_COMPONENT_COMPLETE "memory_safety"
