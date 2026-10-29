* Coefficient and covariance agreement in significant figures. An entry
* that is zero in exact arithmetic (a symmetric design, an exactly fitted
* model) is reported at roundoff by both programs, and relative digits of two
* roundoff values carry no information. Such an entry agrees when both values
* lie below a roundoff bound in the entry's own units, as in
* sigfigs_cancellation (16*epsdouble):
*   b_j:  16*eps*SE_j, roundoff at the coefficient's sampling scale.
*   V_ij: max(16*eps*SE_i*SE_j, r_i*r_j) with r_j = 16*eps*|b_j|. The first
*         term is the benchmark_helpers covariance bound; the second treats a
*         standard error below 16 ulps of its coefficient as unresolvable.
* Scales come from the reference (or analytic oracle) b and V. Every other
* entry must agree to the requested significant figures.
capture program drop ppml_min_sigfigs
program define ppml_min_sigfigs, rclass
    args A B kind refb refV
    tempname a c bound
    local min_sf 15
    local tol=16*c(epsdouble)
    forvalues i=1/`=rowsof(`A')' {
        forvalues j=1/`=colsof(`A')' {
            scalar `a'=`A'[`i',`j']
            scalar `c'=`B'[`i',`j']
            if "`kind'"=="b" scalar `bound'=`tol'*sqrt(`refV'[`j',`j'])
            else scalar `bound'=max(`tol'*sqrt(`refV'[`i',`i'])*sqrt(`refV'[`j',`j']), ///
                `tol'^2*abs(`refb'[1,`i']*`refb'[1,`j']))
            if missing(`bound') scalar `bound'=0
            if !missing(`a',`c') & abs(`a')<=`bound' & abs(`c')<=`bound' continue
            sigfigs `a' `c'
            local min_sf=min(`min_sf',r(sigfigs))
        }
    }
    return scalar min_sigfigs=`min_sf'
end

* Native PPML compatibility: compare full results and retained rows, then
* exercise prediction after estimates store/restore. No reference fallback.
capture program drop ppml_native_compare
program define ppml_native_compare
    syntax varlist(fv ts) [fw pw] [if] [in], NAME(string) [ABSorb(string) VCE(string) OFFset(varname) EXTRA(string) REFerenceopts(string) NATIVEopts(string) PREDICT VCEREF(name)]
    local wt
    if "`weight'"!="" local wt [`weight'`exp']
    local opts `extra'
    if "`absorb'"!="" local opts `opts' absorb(`absorb')
    if "`vce'"!="" local opts `opts' vce(`vce')
    if "`offset'"!="" local opts `opts' offset(`offset')
    local precise_ref
    foreach control in use_exact_partial use_exact_solver min_ok {
        if !strpos("`opts' `referenceopts'","`control'(") {
            local value=cond("`control'"=="min_ok",3,1)
            local precise_ref `precise_ref' `control'(`value')
        }
    }
    preserve
    capture noisily {
        tempvar refsample rd cd
        quietly ppmlhdfe `varlist' `wt' `if' `in', `opts' d(`rd') tol(1e-14) ///
            `precise_ref' `referenceopts'
        matrix rb=e(b)
        matrix rv=e(V)
        * An explicit independent oracle is needed when Stata rejects the
        * upstream covariance at posting. All other result checks still apply.
        if "`vceref'" != "" matrix rv=`vceref'
        gen byte `refsample'=e(sample)
        foreach stat in N df df_a df_a_initial df_a_redundant df_a_nested df_m rank num_singletons num_separated {
            scalar ref_`stat'=e(`stat')
        }
        local prediction_types xb mu eta xbd d stdp response scores pearson deviance anscombe working
        if "`predict'"!="" {
            foreach opt of local prediction_types {
                tempvar rp_`opt'
                quietly predict double `rp_`opt'', `opt'
            }
        }
        quietly cppmlhdfe `varlist' `wt' `if' `in', `opts' `nativeopts' d(`cd') tol(1e-14)
        matrix cb=e(b)
        matrix cv=e(V)
        assert `refsample'==e(sample)
        foreach stat in N df df_a df_a_initial df_a_redundant df_a_nested df_m rank num_singletons num_separated {
            if reldif(e(`stat'),ref_`stat')>=1e-12 noi di "MISMATCH `name': `stat' reference=" ref_`stat' " native=" e(`stat')
            assert reldif(e(`stat'),ref_`stat')<1e-12
        }
        assert colsof(rb)==colsof(cb)
        if mreldif(rb,cb)>=1e-7 noi di "MISMATCH `name': coefficients " mreldif(rb,cb)
        ppml_min_sigfigs rb cb b rb rv
        assert r(min_sigfigs)>=$DEFAULT_SIGFIGS
        ppml_min_sigfigs rv cv V rb rv
        assert r(min_sigfigs)>=$DEFAULT_SIGFIGS
        forvalues j=1/`=colsof(cb)' {
            assert abs(sqrt(cv[`j',`j'])-sqrt(rv[`j',`j'])) < 1e-6*sqrt(rv[`j',`j'])+1e-12
            if "`vceref'" != "" & rv[`j',`j'] > 0 {
                assert abs(cv[`j',`j']/rv[`j',`j']-1) < 1e-7
            }
        }
        if "`predict'"!="" {
            estimates store native_restore
            quietly regress `: word 1 of `varlist''
            estimates restore native_restore
            * Predictions must agree to 7 significant figures (1e-7) relative
            * to max(|entry|, S). The linear index xb + d = eta is built from
            * global quantities (coefficients, fixed effects) whose errors do
            * not shrink where a row's sum crosses zero, so for xb, d, eta and
            * xbd S is their common sup norm over the estimation sample.
            * Without fixed effects d is identically zero and pure roundoff in
            * both programs. Residual-type predictions are differences of y
            * and mu terms, and stdp is the root of a signed quadratic form
            * x'Vx; each uses its own sup norm (normwise relative error). mu
            * has no cancellation and keeps elementwise relative error (S=0).
            local linear_scale 0
            foreach comp in xb d eta {
                quietly summarize `rp_`comp'' if `refsample', meanonly
                if r(N) local linear_scale=max(`linear_scale',abs(r(min)),abs(r(max)))
            }
            foreach opt of local prediction_types {
                local scale 0
                if inlist("`opt'","xb","d","eta","xbd") local scale `linear_scale'
                else if "`opt'"!="mu" {
                    quietly summarize `rp_`opt'' if `refsample', meanonly
                    if r(N) local scale=max(abs(r(min)),abs(r(max)))
                }
                tempvar cp
                quietly predict double `cp', `opt'
                quietly count if missing(`rp_`opt'')!=missing(`cp') | ///
                    (abs(`cp'-`rp_`opt'')>=1e-6*(1+abs(`rp_`opt'')) & !missing(`cp'))
                if r(N) noi di "MISMATCH `name': prediction `opt' in " r(N) " rows"
                assert missing(`rp_`opt'')==missing(`cp')
                assert `cp'==`rp_`opt'' | abs(`cp'-`rp_`opt'') <= 1e-7*max(abs(`cp'),abs(`rp_`opt''),`scale') if !missing(`cp')
            }
            estimates drop native_restore
        }
    }
    local rc=_rc
    restore
    if `rc' test_fail "native `name'" "r(`rc')"
    else test_pass "native `name'"
end

noi di "--- Native parity: separation, absorption, and prediction ---"
forvalues seed=1/10 {
    clear
    set seed `seed'
    set obs 1000
    gen g=runiformint(1,20)
    gen h=runiformint(1,10)
    gen double x=rnormal()
    gen double z=rnormal()
    gen double off=rnormal()*.3
    gen double y=rpoisson(exp(.5*x-.3*z+.03*g+off))
    gen fw=runiformint(1,4)
    gen double pw=exp(rnormal())
    gen byte sep=(_n<=100)
    replace y=0 if sep
    gen double combo=x+sep
    foreach wt in "" "[fw=fw]" "[pw=pw]" {
        ppml_native_compare y x z combo `wt', absorb(g h) vce(cluster g h) offset(off) ///
            name("separating combination seed `seed' `wt'") predict
    }
    drop combo sep
    ppml_native_compare y x z, absorb(g#h) offset(off) name("categorical interaction `seed'") predict
    ppml_native_compare y x, absorb(g##c.z h) vce(cluster g) offset(off) name("heterogeneous slopes `seed'") predict
    ppml_native_compare y x, absorb(g#c.z h) offset(off) name("slope-only component `seed'") predict
    ppml_native_compare y, absorb(g h) offset(off) name("pure fixed effects `seed'") predict
    ppml_native_compare y x z, offset(off) name("no fixed effects `seed'") predict
    ppml_native_compare y, offset(off) name("constant only `seed'") predict
}

* The FE design itself separates one off-diagonal cell. No individual level
* is all-zero, so an all-zero-FE screen cannot find this case.
clear
set seed 2201
set obs 240
gen g=1+(_n>80)
gen h=1+(_n>160)
gen double x=rnormal()
gen double y=cond(inrange(_n,81,160),0,1+runiformint(0,3))
ppml_native_compare y x, absorb(g h) name("separation from FE combination") predict
ppml_native_compare y, absorb(g h) name("pure-FE separation") predict

* Base levels, interaction omissions, disconnected FE graphs, and missing rows.
clear
set seed 9121
set obs 800
gen g=runiformint(1,12)
gen h=runiformint(1,8)+20*(g>6)
gen t=runiformint(1,5)+10*(g>6)
gen double x=rnormal()
gen double z=rnormal()
gen byte cat=runiformint(1,3)
gen double y=rpoisson(exp(.4*x+.2*z+.03*g))
replace y=0 if g==1
replace x=. in 101/105
ppml_native_compare y x z, absorb(g h t) name("disconnected three-way DOF") predict
ppml_native_compare y x z, absorb(g h t) vce(cluster g) name("nested and connected DOF") predict
ppml_native_compare y i.cat##c.x z, absorb(g h) name("factor bases and interactions") predict
ppml_native_compare y x z, absorb(g h) vce(cluster g#h) name("interaction cluster") predict
ppml_native_compare y x z, absorb(g h) extra("keepsingletons") name("keepsingletons") predict

* Exposure must be logged exactly once, including after estimates restore.
gen double exposure=exp(z)
ppml_native_compare y x, absorb(g h) extra("exposure(exposure)") name("exposure prediction") predict
ppml_native_compare y x, extra("exposure(exposure)") name("no-FE exposure prediction") predict
* The installed reference's preconditioner errors for multiple slopes in a
* single factor. Its equivalent explicit slope design remains a valid oracle.
capture noisily {
    tempvar rd cd rs rm cm
    quietly ppmlhdfe y x ibn.g#c.z ibn.g#c.exposure, absorb(g h) d(`rd') ///
        tol(1e-14) use_exact_partial(1) use_exact_solver(1) min_ok(3)
    scalar slope_b=_b[x]
    scalar slope_se=_se[x]
    gen byte `rs'=e(sample)
    quietly predict double `rm', mu
    quietly cppmlhdfe y x, absorb(g##c.(z exposure) h) d(`cd') tol(1e-14)
    assert e(sample)==`rs'
    assert reldif(_b[x],slope_b)<1e-7
    assert abs(_se[x]-slope_se)<1e-6*slope_se+1e-12
    quietly predict double `cm', mu
    assert abs(`cm'-`rm')<1e-6*(1+abs(`rm')) if e(sample)
}
if _rc test_fail "native multiple heterogeneous slopes" "r(`=_rc')"
else test_pass "native multiple heterogeneous slopes"

foreach dof in none firstpair pairwise clusters continuous "pairwise clusters" "firstpair continuous" all {
    ppml_native_compare y x z, absorb(g h t) vce(cluster g) extra("dof(`dof')") ///
        name("DOF adjustment `dof'")
}
capture noisily {
    foreach option in "dof(bogus)" "guess(bogus)" "relu_maxiter(0)" {
        capture cppmlhdfe y x, absorb(g h) `option'
        assert _rc==198
    }
    local scratch : all scalars
    foreach name of local scratch {
        assert substr("`name'",1,12)!="__cppmlhdfe_"
    }
}
if _rc test_fail "native invalid options and scratch cleanup" "r(`=_rc')"
else test_pass "native invalid options and scratch cleanup"
