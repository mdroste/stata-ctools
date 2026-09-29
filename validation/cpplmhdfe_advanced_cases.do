* Explicit reference controls and initialization methods.
clear
set seed 84123
set obs 700
gen g=runiformint(1,20)
gen h=runiformint(1,15)
gen double x=rnormal()
gen double z=rnormal()
gen double y=rpoisson(exp(.4*x-.2*z+.03*g))
gen double initial=exp(.3*x+.02*g)
foreach option in "guess(simple)" "guess(ols)" "guess(variable initial)" ///
    "standardize_data(0)" "min_ok(3)" "use_exact_solver(1)" ///
    "use_exact_partial(0)" "use_heuristic_tol(0)" "start_inner_tol(.001)" ///
    "use_step_halving(1) step_halving_memory(.8) max_step_halving(4)" ///
    "remove_collinear_variables(0)" "separation(mu)" "separation(fe mu)" ///
    "relu_tol(1e-6) relu_zero_tol(1e-10) relu_maxiter(200) relu_strict(1)" {
    ppml_native_compare y x z, absorb(g h) extra("`option'") name("advanced `option'") predict
}
* A separated dummy exercises the explicit mu-only removal and native refit.
gen byte separated=(_n<=70)
replace y=0 if separated
ppml_native_compare y x z separated, absorb(g h) extra("separation(mu)") name("mu-only separation")

* Named effects and absorb suboptions must reconstruct the saved sum in C.
capture noisily {
    cpplmhdfe y x z, absorb(firm_fe=g year_fe=h) d(total_fe) tol(1e-12)
    assert missing(firm_fe)==!e(sample)
    assert missing(year_fe)==!e(sample)
    assert abs(firm_fe+year_fe-total_fe)<1e-7 if e(sample)
    bysort g: assert firm_fe==firm_fe[1] if e(sample) & !missing(firm_fe[1])
    drop firm_fe year_fe total_fe
    cpplmhdfe y x, absorb(trend_fe=g##c.z year_fe=h, tolerance(1e-12)) d(total_fe) tol(1e-12)
    assert abs(trend_fe+trend_feSlope1*z+year_fe-total_fe)<1e-7 if e(sample)
    drop trend_fe trend_feSlope1 year_fe total_fe
    cpplmhdfe y x, absorb(g h, savefe keepsingletons) d(total_fe) tol(1e-12)
    assert abs(__hdfe1__+__hdfe2__-total_fe)<1e-7 if e(sample)
    assert e(drop_singletons)==0
    cpplmhdfe y x, absorb(g h, savefe) tol(1e-12)
    assert !missing(__hdfe1__,__hdfe2__) if e(sample)
    drop total_fe
    gen group_one=g
    gen group_two=h
    cpplmhdfe y x, absorb(group_*, savefe) d(total_fe) tol(1e-12)
    assert abs(__hdfe1__+__hdfe2__-total_fe)<1e-7 if e(sample)
    drop total_fe
    cpplmhdfe y x, absorb(trend_fe=g##c.(z initial) year_fe=h) d(total_fe) tol(1e-12)
    assert abs(trend_fe+trend_feSlope1*z+trend_feSlope2*initial+year_fe-total_fe)<1e-7 if e(sample)
}
if _rc test_fail "native named and saved fixed effects" "r(`=_rc')"
else test_pass "native named and saved fixed effects"

* Individual separation methods are independent options.
foreach method in simplex relu "fe simplex" "fe relu" "fe simplex relu mu" {
    ppml_native_compare y x z separated, absorb(g h) extra("separation(`method')") name("separation `method'")
}

foreach option in "accel(cg)" "accel(none)" "accel(sd)" "accel(aitken)" ///
    "accel(hybrid)" "accel(lsmr) btol(1e-12)" "accel(sd) transform(kaczmarz)" ///
    "transform(cimmino)" "poolsize(1)" {
    ppml_native_compare y x z separated, absorb(g h) nativeopts("`option'") name("projection `option'")
}
capture noisily {
    cpplmhdfe y x z separated, absorb(g h) separation(relu) ///
        tagsep(is_separated) zvarname(certificate) relu_accelerate(1) tol(1e-12)
    assert is_separated==separated
    assert abs(certificate)<1e-10 if y>0
    assert certificate>0 if separated
    assert certificate>=0 if y==0
}
if _rc test_fail "native ReLU tags and certificate" "r(`=_rc')"
else test_pass "native ReLU tags and certificate"

* tagsep selects diagnostic-only ReLU and overrides FE/simplex screening.
preserve
capture noisily {
    replace y=0 if g==1
    ppmlhdfe y x z separated, absorb(g h) separation(fe relu) tagsep(ref_tag) zvarname(ref_certificate)
    cpplmhdfe y x z separated, absorb(g h) separation(fe relu) tagsep(native_tag) zvarname(native_certificate)
    assert native_tag==ref_tag
    assert native_tag==1 if g==1
    assert missing(native_certificate)==missing(ref_certificate)
    assert native_certificate==0 if y>0 & !missing(ref_certificate)
    assert native_certificate>0 if native_tag==1
}
if _rc test_fail "native ReLU diagnostic sample" "r(`=_rc')"
else test_pass "native ReLU diagnostic sample"
restore

* Solver controls also interact with weights, offset, slopes, and warm partials.
clear
set seed 88231
set obs 1200
gen g=runiformint(1,20)
gen h=runiformint(1,15)
gen double x=rnormal()
gen double z=rnormal()
gen double off=.2*rnormal()
gen double pw=exp(rnormal())
gen double y=rpoisson(exp(.4*x+.01*g*z+.02*h+off))
gen double poor_initial=1000+_n
foreach method in cg none sd aitken hybrid lsmr {
    ppml_native_compare y x [pw=pw], absorb(g##c.z h) offset(off) vce(cluster g h) ///
        nativeopts("accel(`method') use_exact_partial(0) use_exact_solver(1) btol(1e-12)") ///
        name("weighted slope projection `method'") predict
}
ppml_native_compare y x z [pw=pw], absorb(g h) offset(off) ///
    extra("guess(variable poor_initial) use_step_halving(1) max_step_halving(8)") ///
    name("poor initial means with step halving") predict

capture noisily {
    foreach invalid in "tolerance(0)" "iterate(0)" "poolsize(-1)" {
        capture cpplmhdfe y x, absorb(g h, `invalid')
        assert _rc==198
    }
    foreach invalid in "btol(.)" "conlim(.)" "relu_tol(0)" "guess(nonsense)" "relu_maxiter(.)" "threads(-1)" {
        capture cpplmhdfe y x, absorb(g h) `invalid'
        assert _rc==198
    }
}
if _rc test_fail "advanced invalid controls" "r(`=_rc')"
else test_pass "advanced invalid controls"
