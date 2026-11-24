* Analytic regression checks for the November 2026 four-command audit.
* Called by each public suite; deliberately independent of reference bugs.
args command
if "`command'" == "cqreg" {
    clear
    set seed 71234
    set obs 240
    generate double x0 = rnormal()
    generate double y = 3 + 2*x0 + rnormal()
    generate double x = x0
    generate int cluster = ceil(_n/8)
    generate double w = 1
    foreach vce in "iid" "iid, residual" "iid, kernel(epan)" "robust" "robust, kernel(epan)" "cluster cluster" {
        replace x = x0
        quietly cqreg y x, vce(`vce')
        scalar audit_b = _b[x]
        scalar audit_se = _se[x]
        scalar audit_obj = e(sum_adev)
        foreach offset in 1e6 1e8 1e12 {
            replace x = x0 + `offset'
            quietly cqreg y x, vce(`vce')
            * At 1e12 the double input itself loses about 1e-4 precision.
            local tol = cond(`offset'==1e12, 5e-4, 1e-6)
            assert reldif(_b[x], audit_b) < `tol'
            assert reldif(_se[x], audit_se) < `tol'
            assert reldif(e(sum_adev), audit_obj) < `tol'
            assert e(convcode)==0
            test_pass "cqreg translated design `offset' `vce'"
        }
        foreach scale in 1e-100 1e100 {
            replace x = x0 * `scale'
            quietly cqreg y x, vce(`vce')
            assert reldif(_b[x]*`scale', audit_b) < 1e-6
            assert reldif(_se[x]*`scale', audit_se) < 1e-6
            assert reldif(e(sum_adev), audit_obj) < 1e-6
            test_pass "cqreg design units `scale' `vce'"
        }
    }
    replace x = x0
    foreach weight in iw pw {
        quietly cqreg y x [`weight'=w]
        scalar audit_b = _b[x]
        scalar audit_se = _se[x]
        scalar audit_obj = e(sum_adev)
        foreach scale in 1e-20 1e-10 1 1e10 {
            replace w = `scale'
            quietly cqreg y x [`weight'=w]
            assert reldif(_b[x], audit_b) < 1e-7
            * qreg's IID iweight VCE uses X'diag(w^2)X; retain its units.
            if "`weight'"=="iw" assert reldif(_se[x]*`scale', audit_se) < 1e-7
            else assert reldif(_se[x], audit_se) < 1e-7
            assert reldif(e(sum_adev)/`scale', audit_obj) < 1e-7
            test_pass "cqreg common weight scale `scale' `weight'"
        }
        replace w=1
    }
    replace w=1+mod(_n,7)
    quietly cqreg y x [iw=w], vce(robust)
    scalar audit_b=_b[x]
    scalar audit_se=_se[x]
    scalar audit_obj=e(sum_adev)
    replace w=w*1e-20
    quietly cqreg y x [iw=w], vce(robust)
    assert reldif(_b[x],audit_b)<1e-7
    assert reldif(_se[x],audit_se)<1e-7
    assert reldif(e(sum_adev)/1e-20,audit_obj)<1e-7
    test_pass "cqreg heterogeneous weight scale preserves robust inference"
    generate double z=rnormal()
    quietly cqreg y x z, vce(robust)
    matrix audit_b0=e(b)
    matrix audit_v0=e(V)
    matrix audit_T=I(3)
    matrix audit_T[1,1]=1/100
    matrix audit_T[3,1]=-1e6
    replace x=100*(x0+1e6)
    quietly cqreg y x z, vce(robust)
    matrix audit_expected_b=audit_b0*audit_T'
    matrix audit_expected_v=audit_T*audit_v0*audit_T'
    assert mreldif(e(b),audit_expected_b)<1e-6
    assert mreldif(e(V),audit_expected_v)<1e-6
    test_pass "cqreg full coefficient and covariance unit transformation"
    replace w=1e307
    capture cqreg y x [iw=w]
    assert _rc==498
    test_pass "cqreg rejects unrepresentable weight sums"
    * Unrepresentable covariance must fail without replacing prior estimates.
    quietly cqreg y x0
    matrix audit_prior=e(b)
    foreach scale in 1e-200 1e200 {
        replace x=x0*`scale'
        capture cqreg y x
        assert _rc==498
        assert mreldif(e(b), audit_prior)==0
        test_pass "cqreg rejects unrepresentable regressor units `scale'"
    }
}

if "`command'" == "cpsmatch" {
    clear
    input double(t ps y _y)
    0 .1 10 100
    0 .8 20 200
    1 .2 30 300
    1 .7 40 400
    end
    capture cpsmatch t, pscore(ps) outcome(y _y)
    assert _rc==111
    assert _y==100*_n
    rename ps _n1
    capture cpsmatch t, pscore(_n1) outcome(y)
    assert _rc==111
    assert _n1==cond(_n==1,.1,cond(_n==2,.8,cond(_n==3,.2,.7)))
    rename _n1 ps
    generate double _weight=ps
    capture cpsmatch t c._weight, outcome(y)
    assert _rc==111
    assert _weight==ps
    generate double _s_y=y+7
    capture cpsmatch t, pscore(ps) outcome(y _s_y) llr kerneltype(epan)
    assert _rc==111
    assert _s_y==y+7
    generate double _self_y=y+9
    capture cpsmatch t, mahal(ps) outcome(y _self_y) ai(1)
    assert _rc==111
    assert _self_y==y+9
    test_pass "cpsmatch preserves aliased outcome, score, FV, smooth and AI inputs"

    clear
    input double(t ps y)
    0 .1 10
    0 .1 100
    0 .1 1000
    1 .2 2000
    end
    foreach k in 1 2 3 {
        quietly cpsmatch t, pscore(ps) outcome(y) ties neighbor(`k')
        assert _nn==3 if t==1
        assert abs(_weight-1/3)<1e-12 if t==0
        assert _y==370 if t==1
        assert r(att)==1630
        test_pass "cpsmatch unique tied donors neighbor(`k')"
    }
    replace ps=.2 in 2/3
    replace ps=.1 in 4
    quietly cpsmatch t, pscore(ps) outcome(y) ties neighbor(2)
    assert _nn==3 if t==1
    assert _y==370 if t==1
    test_pass "cpsmatch includes kth boundary ties"

    clear
    input double(t x y)
    0 0 10
    0 2 20
    0 5 50
    1 .1 100
    1 1.9 200
    1 4.9 500
    end
    generate double original=x
    matrix metric=(2)
    foreach opt in "" "w(metric)" "kernel kerneltype(normal) bwidth(2)" "llr kerneltype(normal) bwidth(2)" {
        replace x=original
        quietly cpsmatch t, mahal(x) outcome(y) `opt'
        scalar audit_att=r(att)
        generate double audit_y=_y
        foreach shift in 1e6 1e8 {
            replace x=original+`shift'
            quietly cpsmatch t, mahal(x) outcome(y) `opt'
            assert reldif(r(att), audit_att)<1e-7
            assert reldif(_y,audit_y)<1e-7 if t==1
            test_pass "cpsmatch translated Mahalanobis `shift' `opt'"
        }
        drop audit_y
    }
}

if "`command'" == "chdfe" {
    clear
    set obs 8
    generate g=ceil(_n/2)
    generate double y=_n
    generate double w=cond(_n<=4,1e-200,1e200)
    capture chdfe y [aw=w], absorb(g) gen(r_)
    assert _rc==498
    capture confirm variable r_y
    assert _rc==111
    assert y==_n
    test_pass "chdfe rejects weights that erase small groups"
    clear
    set obs 12
    generate g=ceil(_n/3)
    generate double x=cond(g<=2,1e-200,1e200)*mod(_n,3)
    generate double y=mod(_n,3)
    capture chdfe y, absorb(g#c.x) gen(r_) keepsingletons
    assert _rc==498
    capture confirm variable r_y
    assert _rc==111
    assert y==mod(_n,3)
    test_pass "chdfe rejects slopes that erase small groups"
    generate double w=cond(g<=2,1e-150,1)
    replace x=cond(g<=2,1e-150,1)*mod(_n,3)
    capture chdfe y [aw=w], absorb(g#c.x) gen(r_) keepsingletons
    assert _rc==498
    test_pass "chdfe rejects weighted slope norm underflow"
    replace x=mod(_n,3)
    replace y=cond(g<=2,1e-200,1e200)*mod(_n,3)
    capture chdfe y, absorb(g) clear
    assert _rc==498
    assert _N==12
    assert y==cond(g<=2,1e-200,1e200)*mod(_n,3)
    test_pass "chdfe rejects response underflow before clear publication"
}

if "`command'" == "crangejoin" {
    tempfile stem
    global CTOOLS_AUDIT_PATH expanded
    foreach kind in global local {
    if "`kind'"=="global" {
    mata: st_local("literal", st_local("stem") + "_" + char(36) + "CTOOLS_AUDIT_PATH.dta")
    local expanded `"`stem'_expanded.dta"'
    }
    else {
        mata: st_local("literal", st_local("stem") + "_" + char(96) + "CTOOLS_AUDIT_LOCAL_PATH" + char(39) + ".dta")
        local expanded `"`stem'_.dta"'
    }
    clear
    set obs 1
    generate double key=1
    generate double payload=123
    quietly save `"`macval(literal)'"'
    replace payload=999
    quietly save `"`expanded'"'
    quietly use `"`macval(literal)'"', clear
    assert payload==123
    foreach mode in "" "snapshot" "fast" {
        foreach keep in "" "keepusing(payload)" {
            clear
            set obs 1
            generate double key=1
            crangejoin key 0 0 using `"`macval(literal)'"', `mode' `keep'
            assert payload==123 & _N==1
            test_pass "crangejoin literal `kind' filename `mode' `keep'"
        }
    }
    erase `"`macval(literal)'"'
    erase `"`expanded'"'
    }
    macro drop CTOOLS_AUDIT_PATH
}
