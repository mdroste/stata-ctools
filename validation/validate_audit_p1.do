/* Targeted regression tests for audit findings A01--A17. Run from the root. */
do validation/validate_setup.do

program p1_merge
    clear
    set obs 2
    gen id = _n
    gen val = 100 + _n
    tempfile usingdata
    save `usingdata'
    replace val = 999
    tempname master sentinel
    frame copy default `master'
    frame create `sentinel'
    frame change `master'
    drop val
    cmerge 1:1 id using `usingdata'
    assert "`c(frame)'" == "`master'"
    assert val == 100 + id
    frame default: assert val == 999
    drop _merge
    capture cmerge 1:1 id using `usingdata', keepusing(not_a_variable)
    assert _rc != 0
    assert "`c(frame)'" == "`master'"
    frame default: assert val == 999
    frame change default
end

program p1_sample
    foreach qualifier in "if id<=5" "in 1/5" "if id<=7 in 1/5" {
        foreach pct in 0 100 {
            clear
            set obs 10
            gen id = _n
            gen __csample_keep__ = 123
            csample `pct' `qualifier'
            assert __csample_keep__ == 123
            count if id > 5
            assert r(N) == 5
            assert _N == 5 + (`pct' == 100)*5
        }
    }
end

program p1_groups
    foreach cmd in csample cbsample {
        foreach kind in string numeric mixed missing {
            clear
            set obs 10
            gen id = _n
            gen str1 s = cond(_n<=5,"A","B")
            gen double g = cond(_n<=5,.1,1e12)
            if "`kind'" == "missing" replace g = cond(_n<=5,.a,.b)
            local groups s
            if inlist("`kind'","numeric","missing") local groups g
            if "`kind'" == "mixed" local groups s g
            if "`cmd'" == "cbsample" & "`kind'" == "missing" {
                replace g = 1 in 1/5
                cbsample 1, strata(g)
                assert _N == 1 & g == 1
                continue
            }
            if "`cmd'" == "csample" csample, count(1) by(`groups')
            else cbsample 1, strata(`groups')
            assert _N == 2
            isid `groups', missok
        }
    }
end

program p1_matching
    foreach option in "" "noreplacement" "noreplacement descending" {
        clear
        set obs 3
        gen byte treat = _n==3
        gen double ps = cond(_n==1,.1,cond(_n==2,.51,.5))
        gen double y = cond(_n==1,10,cond(_n==2,20,100))
        cpsmatch treat, pscore(ps) outcome(y) `option'
        assert r(att) == 80
        assert _weight[2] == 1 & _weight[3] == 1
        assert _support == 1
        assert _id[3] == 2
    }
    clear
    set obs 4
    gen byte treat = _n>=3
    gen double ps = cond(_n==1,.1,cond(_n==2,.51,cond(_n==3,.11,.5)))
    gen double y = cond(_n==1,10,cond(_n==2,20,cond(_n==3,100,200)))
    cpsmatch treat, pscore(ps) outcome(y) noreplacement
    assert _id[3] == 1 & _id[4] == 2
    assert _weight == 1
    assert r(att) == 135
    cpsmatch treat, pscore(ps) outcome(y) noreplacement caliper(.001)
    assert _support == 0 if treat
    capture cpsmatch treat, pscore(ps) noreplacement neighbor(2)
    assert _rc == 198
    capture cpsmatch treat, pscore(ps) noreplacement kernel
    assert _rc == 198
    capture cpsmatch treat, pscore(ps) noreplacement radius
    assert _rc == 198
    clear
    set obs 5
    gen byte treat = _n>2
    gen double ps = cond(_n==1,.25,cond(_n==2,.75,.5))
    cpsmatch treat, pscore(ps) noreplacement
    assert _id[3]==1 & _id[4]==2
    assert _support[5]==0
    assert r(n_matched)==2
    clear
    set seed 12345
    set obs 500
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen byte treat = runiform() < normal(.3*x1+.2*x2)
    gen double y = 2*treat+x1+x2+rnormal()
    gen double saved_y = y
    cpsmatch treat x1 x2, outcome(y)
    assert y == saved_y
    assert inrange(_pscore,0,1)
    assert _support == 1
    assert !missing(_weight, r(att))
end

program p1_labels
    clear
    set obs 3
    gen str1 s = char(64+_n)
    label define existing 10 "A" 20 "B" 30 "unused"
    gen other = 30
    label values other existing
    encode s, gen(expected) label(existing)
    cencode s, gen(actual) label(existing)
    assert actual == expected
    local unused : label existing 30
    assert "`unused'" == "unused"
    cencode s, gen(again) label(existing)
    assert again == actual
    gen str1 s2 = s
    cencode s s2, gen(multi1 multi2) label(existing)
    assert multi1 == expected & multi2 == expected
end

program p1_bins
    clear
    set obs 40
    gen double x = _n
    gen double y = 2*x + mod(_n,3)
    gen double g = cond(_n<=20,10,20)
    cbinscatter y x, by(g) nquantiles(4) nograph
    assert e(N) == 40 & e(num_groups) == 2
    matrix bins = e(bindata)
    assert rowsof(bins) == 8
    scalar total = 0
    forvalues i=1/8 {
        assert !missing(bins[`i',3],bins[`i',4],bins[`i',5])
        scalar total = total + bins[`i',5]
    }
    assert total == 40
    replace g = cond(_n<=20,-.1,1e12)
    cbinscatter y x, by(g) nquantiles(4) nodraw
    assert e(num_groups) == 2
    gen str1 gs = cond(_n<=20,"A","B")
    cbinscatter y x, by(gs) nquantiles(4) nograph
    assert e(num_groups) == 2
end

program p1_import
    tempfile csv xlsx
    foreach encoding in ascii utf8 {
        foreach width in 2044 2045 2046 3000 {
            foreach format in csv xlsx {
                clear
                set obs 1
                gen strL text = ""
                if "`encoding'" == "ascii" {
                    mata: st_sstore(1,"text",invtokens(J(1,`width',"a"),""))
                }
                else {
                    local symbol = uchar(233)
                    local chars = floor(`width'/2)
                    mata: st_sstore(1,"text",invtokens(J(1,`chars',"`symbol'"),""))
                    if mod(`width',2) replace text = text + "a"
                }
                local expected = text[1]
                if "`format'" == "csv" export delimited using `csv', replace
                else export excel using "`xlsx'.xlsx", firstrow(variables) replace
                clear
                set obs 1
                gen sentinel = 42
                if "`format'" == "csv" {
                    capture noisily cimport delimited using `csv', clear
                    local rc = _rc
                }
                else {
                    capture noisily cimport excel using "`xlsx'.xlsx", firstrow clear
                    local rc = _rc
                }
                if `width' <= 2045 {
                    assert `rc' == 0
                    assert strlen(text) == `width'
                    assert text == "`expected'"
                }
                else {
                    assert `rc' == 198
                    assert _N == 1 & sentinel == 42
                }
            }
        }
    }
    erase "`xlsx'.xlsx"
end

program p1_iv
    clear
    set seed 2050
    set obs 400
    gen group = ceil(_n/20)
    replace group = 999 in 400
    gen double z = rnormal()
    gen double x = z + rnormal()
    gen double y = 10*group + 2*x + rnormal()
    replace z = . in 6
    civreghdfe y (x=z) if _n>3, absorb(group) residuals(savedres)
    count if e(sample)
    assert r(N) == e(N) & r(N) == 395
    assert !missing(savedres) if _n>3 & _n!=6 & _n!=400
    assert missing(savedres) if _n<=3 | _n==6 | _n==400
    assert !e(sample) if _n<=3 | _n==6 | _n==400
    assert "`e(resid)'" == "savedres"
    predict double residual, residuals
    predict double fit, xbd
    assert residual == savedres
    assert abs(y-fit-residual) < 1e-10 if e(sample)
    tempvar meanres
    egen double `meanres' = mean(residual), by(group)
    assert abs(`meanres') < 1e-7 if e(sample)
    estimates store iv_saved
    quietly regress y x
    estimates restore iv_saved
    predict double restored, residuals
    assert restored == savedres
    gen _civreghdfe_resid = 999
    civreghdfe y (x=z), absorb(group) residuals2
    count if e(sample)
    assert r(N) == e(N) & r(N) == 398
    assert !missing(_civreghdfe_resid) if e(sample)
    assert _civreghdfe_resid != 999 if e(sample)
end

program p1_variance
    foreach n in 10 200 {
        foreach exclude in "" "excludeself" {
            clear
            set obs `n'
            gen time = _n
            gen double x = 1e12 + _n
            gen double small = _n
            crangestat (variance) vx=x (variance) vs=small (sd) sx=x, interval(time . .) `exclude'
            assert abs(vx-vs) < 1e-8
            assert abs(sx^2-vx) < 1e-8
            summarize small
            if "`exclude'" == "" assert abs(vx-r(Var)) < 1e-8
            crangestat (variance) wx=x (variance) ws=small, interval(time -2 2) `exclude'
            assert abs(wx-ws) < 1e-8
        }
    }
    clear
    set obs 400000
    gen time = _n
    gen double x = 1e12 + _n
    crangestat (variance) vx=x, interval(time -99 0)
    assert abs(vx-100*101/12) < 1e-8 in 399900/400000

end

program p1_quantile
    clear
    set seed 54321
    set obs 400
    gen group = ceil(_n/20)
    gen double x = rnormal()
    gen double y = 2*x + rnormal()
    capture cqreg y x, absorb(group)
    assert _rc == 198
    capture cqreg y x, maxiter(1)
    assert _rc == 430
    cqreg y x
    assert e(convcode) == 0
    assert e(iterations) >= 0 & e(iterations) <= 200
    qreg y x
    scalar expected = _b[x]
    cqreg y x
    assert reldif(_b[x],expected) < 1e-7
end

/* The helper self-tests intentionally record failures, then reset the counters. */
clear
set obs 1
gen a = .
gen b = 42
assert_var_equal a b 7 "missing mismatch"
assert $TESTS_FAILED == 1
matrix a = (1,2)
matrix b = (1,2,3)
assert_matrix_equal a b 7 "dimension mismatch"
assert $TESTS_FAILED == 2
matrix_min_sigfigs a b
assert r(min_sigfigs) == 0
capture print_summary "expected helper failures"
assert _rc == 9
global TESTS_PASSED = 0
global TESTS_FAILED = 0
global TESTS_TOTAL = 0
global FAILURE_COUNT = 0

foreach case in merge sample groups matching labels bins import iv variance quantile {
    capture noisily p1_`case'
    if _rc test_fail "P1 `case'" "r(`=_rc')"
    else test_pass "P1 `case'"
}
print_summary "audit P1 regressions"

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "audit_p1"
