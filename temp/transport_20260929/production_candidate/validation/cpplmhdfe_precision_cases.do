* At default tolerance compare SEs and scale covariance errors by marginal SEs.
* Tight fits below also check relative precision of every covariance entry.
capture program drop ppml_default_check
program define ppml_default_check
    syntax varlist [fw pw], absorb(varlist) offset(varname) vce(string) name(string)
    local wt
    if "`weight'"!="" local wt [`weight'`exp']
    preserve
    capture noisily {
        quietly ppmlhdfe `varlist' `wt', absorb(`absorb') offset(`offset') vce(`vce')
        matrix default_b=e(b)
        matrix default_V=e(V)
        tempvar reference_sample
        gen byte `reference_sample'=e(sample)
        quietly cpplmhdfe `varlist' `wt', absorb(`absorb') offset(`offset') vce(`vce')
        assert `reference_sample'==e(sample)
        matrix actual_b=e(b)
        matrix actual_V=e(V)
        forvalues i=1/2 {
            scalar relative_se=abs(sqrt(actual_V[`i',`i']/default_V[`i',`i'])-1)
            scalar PPML_MAX_DEFAULT_SE=max(PPML_MAX_DEFAULT_SE,relative_se)
            assert relative_se < 1e-5
            assert abs(actual_b[1,`i']-default_b[1,`i']) < 1e-5*max(abs(default_b[1,`i']),1e-3)
            forvalues j=1/2 {
                assert abs(actual_V[`i',`j']-default_V[`i',`j']) < ///
                    1e-5*sqrt(default_V[`i',`i']*default_V[`j',`j'])
            }
        }
    }
    local rc=_rc
    restore
    if `rc' test_fail "default `name'" "r(`rc')"
    else test_pass "default `name'"
end
scalar PPML_MAX_DEFAULT_SE=0

* Invoked by validate_cpplmhdfe.do after loading the common test helpers.
* Tight reference fits distinguish estimator errors from stopping-rule error.
noi di as text "--- Precision: randomized weights, VCE, and FE dimensions ---"
forvalues seed = 1/40 {
    clear
    set seed `seed'
    set obs 600
    gen g = runiformint(1,20)
    gen h = runiformint(1,15)
    gen t = runiformint(1,8)
    gen c = runiformint(1,25)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double off = rnormal()*1.5
    gen double y = rpoisson(exp(-1+.7*x-.3*z+.03*g-.02*h+off))
    gen fw = runiformint(1,10)
    gen double pw = exp(rnormal())
    local fes g
    if mod(`seed',3)>0 local fes g h
    if mod(`seed',3)>1 local fes g h t
    foreach w in none fw pw {
        local wt
        if "`w'"!="none" local wt [`w'=`w']
        foreach vc in robust "cluster g" "cluster c" {
            ppml_default_check y x z `wt', absorb(`fes') offset(off) ///
                vce(`vc') name("seed `seed' `w' `vc'")
            benchmark_ppmlhdfe y x z `wt', absorb(`fes') offset(off) ///
                vce(`vc') tolerance(1e-14) minsf(6) ///
                referenceopts("use_exact_partial(1) use_exact_solver(1) min_ok(2)") ///
                testname("precision seed `seed' `w' `vc'")
        }
    }
}

noi di as text "--- Precision: multiway clustering and sample alignment ---"
forvalues seed = 1/10 {
    clear
    set seed `seed'
    set obs 700
    gen g = runiformint(1,25)
    gen h = runiformint(1,17)
    gen t = runiformint(1,9)
    gen c = runiformint(1,40)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double off = rnormal()*.4
    gen double y = rpoisson(exp(.5*x-.3*z+.04*g+off))
    gen fw = runiformint(1,5)
    gen double pw = exp(rnormal())
    * Drop blocks and missing cluster values, not just a prefix of the data.
    replace y = 0 if g==3 | g==9
    replace c = . in 351/360
    gen str8 sh = string(h)
    gen str8 sg = string(g)
    foreach w in none fw pw {
        local wt
        if "`w'"!="none" local wt [`w'=`w']
        foreach vc in "cluster g h" "cluster c sh" "cluster g h t" {
            benchmark_ppmlhdfe y x z `wt', absorb(sg t) offset(off) ///
                vce(`vc') tolerance(1e-14) minsf(6) ///
                referenceopts("use_exact_partial(1) use_exact_solver(1) min_ok(2)") ///
                testname("multiway seed `seed' `w' `vc'")
        }
    }
}

noi di as text "--- Precision: weight scaling and inference ---"
capture noisily {
    clear
    set seed 99
    set obs 1000
    gen g=mod(_n,20)
    gen h=mod(_n,7)
    gen double x=rnormal()
    gen double z=rnormal()
    gen double y=rpoisson(exp(.5*x-.3*z+.1*g))
    gen double pw=exp(rnormal())
    gen double small=pw*1e-12
    gen double large=pw*1e12
    foreach vc in robust "cluster g" "cluster g h" {
        ppmlhdfe y x z [pw=pw], absorb(g h) vce(`vc') tolerance(1e-14)
        matrix expected_b=e(b)
        matrix expected_V=e(V)
        scalar expected_chi2=e(chi2)
        foreach w in pw small large {
            cpplmhdfe y x z [pw=`w'], absorb(g h) vce(`vc') tolerance(1e-14)
            matrix actual_b=e(b)
            matrix actual_V=e(V)
            assert mreldif(actual_b,expected_b) < 1e-8
            assert mreldif(actual_V,expected_V) < 1e-9
            assert reldif(e(chi2),expected_chi2) < 1e-7
            assert missing(e(df_r))
            test x z
            assert reldif(r(chi2),e(chi2)) < 1e-10
        }
    }
}
if _rc test_fail "PPML weight-scale invariance and Wald inference" "r(`=_rc')"
else test_pass "PPML weight-scale invariance and Wald inference"

capture noisily {
    * The reference may itself omit heavily shifted regressors. Check the
    * estimator's shift invariance against the well-conditioned reference.
    ppmlhdfe y x z, absorb(g) tolerance(1e-14)
    matrix expected_b=e(b)
    matrix expected_V=e(V)
    gen double shifted_x=x+1e6
    gen double shifted_z=z+1e6
    cpplmhdfe y shifted_x shifted_z, absorb(g) tolerance(1e-14)
    matrix shifted_b=e(b)
    matrix shifted_V=e(V)
    matrix shifted_b=shifted_b[1,1..2]
    matrix shifted_V=shifted_V[1..2,1..2]
    matrix expected_b=expected_b[1,1..2]
    matrix expected_V=expected_V[1..2,1..2]
    assert mreldif(shifted_b,expected_b) < 1e-7
    assert mreldif(shifted_V,expected_V) < 1e-9
}
if _rc test_fail "PPML large regressor means" "r(`=_rc')"
else test_pass "PPML large regressor means"

capture noisily {
    gen one=1
    capture cpplmhdfe y x, absorb(g) vce(cluster g one)
    assert _rc==459
    capture cpplmhdfe y x, absorb(g, bogus_option)
    assert _rc==198
    capture cpplmhdfe y x, absorb(g) offset(pw) exposure(pw)
    assert _rc==198
}
if _rc test_fail "PPML invalid options fail explicitly" "r(`=_rc')"
else test_pass "PPML invalid options fail explicitly"

noi di "PPML_MAX_DEFAULT_RELATIVE_SE_ERROR=" %12.5g PPML_MAX_DEFAULT_SE
capture noisily {
    * Scaling the response changes the absorbed intercept only.
    foreach vc in robust "cluster g" "cluster g h" {
        ppmlhdfe y x z, absorb(g h) vce(`vc') tolerance(1e-14)
        matrix expected_b=e(b)
        matrix expected_V=e(V)
        foreach scale in 1e-10 1e10 {
            tempvar scaled_y
            gen double `scaled_y'=y*`scale'
            cpplmhdfe `scaled_y' x z, absorb(g h) vce(`vc') tolerance(1e-14)
            matrix scaled_b=e(b)
            matrix scaled_b[1,colsof(scaled_b)]=scaled_b[1,colsof(scaled_b)]-ln(`scale')
            assert mreldif(scaled_b,expected_b)<1e-8
            assert mreldif(e(V),expected_V)<1e-9
        }
    }
}
if _rc test_fail "PPML response-scale invariance" "r(`=_rc')"
else test_pass "PPML response-scale invariance"

* General separation must select the same sample and refit in C.
capture noisily {
    gen byte separates=(y==0 & x<0)
    ppmlhdfe y x z separates, absorb(g h) tol(1e-12)
    gen byte reference_sample=e(sample)
    matrix reference_b=e(b)
    matrix reference_V=e(V)
    cpplmhdfe y x z separates, absorb(g h) tol(1e-12)
    assert e(sample)==reference_sample
    assert mreldif(e(b),reference_b)<1e-7
    assert mreldif(e(V),reference_V)<1e-7
}
if _rc test_fail "PPML native general separation" "r(`=_rc')"
else test_pass "PPML native general separation"
