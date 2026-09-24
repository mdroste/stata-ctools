/* Regression contracts for SEP22_ASTRA_AUDIT.md. No network data required. */
do validation/validate_setup.do

capture noisily {
    tempfile using_data
    foreach relation in 1:1 m:1 1:m {
        foreach side in master using {
            if ("`relation'" == "m:1" & "`side'" == "master") | ("`relation'" == "1:m" & "`side'" == "using") continue
            clear
            set obs 3
            gen double id = _n
            if "`side'" == "using" replace id = 99 in 2/3
            gen str8 us = "using"
            save `using_data', replace
            clear
            set obs 3
            gen double id = _n
            if "`side'" == "master" replace id = 99 in 2/3
            gen long row = _n
            quietly datasignature
            local before "`r(datasignature)'"
            capture cmerge `relation' id using `using_data'
            assert _rc == 459
            quietly datasignature
            assert "`r(datasignature)'" == "`before'"
        }
    }
    * Missing string keys obey the same cardinality contract.
    clear
    set obs 2
    gen str5 key = ""
    gen row = _n
    save `using_data', replace
    capture cmerge 1:1 key using `using_data'
    assert _rc == 459
    assert _N == 2 & row == _n
}
if _rc test_fail "S01 merge cardinality and rollback" "r(`=_rc')"
else test_pass "S01 merge cardinality and rollback"

capture noisily {
    sysuse auto, clear
    gen str8 group = cond(foreign, "foreign", "domestic")
    clonevar price_before = price
    cwinsor price, suffix(_w) by(foreign group)
    assert price == price_before
    assert !missing(price_w)
    drop price_w
    order price, last
    cwinsor price if mpg > 20, suffix(_w)
    assert missing(price_w) if mpg <= 20
    cwinsor price, replace
    assert !missing(price)
    * A later invalid target cannot leave an earlier generated variable behind.
    drop price_w
    gen mpg_w = 1
    capture cwinsor price mpg, suffix(_w)
    assert _rc == 110
    capture confirm variable price_w
    assert _rc == 111
    assert mpg_w == 1
}
if _rc test_fail "S02 winsor destinations and rollback" "r(`=_rc')"
else test_pass "S02 winsor destinations and rollback"

capture noisily {
    clear
    set obs 400
    set seed 45678
    gen g = ceil(_n/20)
    gen double cl = g + cond(mod(_n,2), .1, .2)
    egen cid = group(cl)
    gen double x = rnormal()
    gen double y = .4*x + g/10 + rnormal()
    gen double yp = rpoisson(exp(.3*x+g/100))
    foreach cmd in creghdfe cpplmhdfe {
        local dep = cond("`cmd'" == "creghdfe", "y", "yp")
        `cmd' `dep' x, absorb(g) vce(cluster cl)
        assert e(N_clust) == 40
        matrix V = e(V)
        `cmd' `dep' x, absorb(g) vce(cluster cid)
        assert e(N_clust) == 40
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
}
if _rc test_fail "S03 cluster partition invariance" "r(`=_rc')"
else test_pass "S03 cluster partition invariance"

capture noisily {
    clear
    set obs 400
    set seed 45679
    gen g = ceil(_n/20)
    gen double frac = g/100
    gen double large = 1e12+g
    gen double neg = -g/100
    gen str8 strcl = string(g)
    gen h = mod(_n,17)
    gen double hfrac = h/100
    gen double z = rnormal()
    gen double x = z + rnormal()
    gen double y = .4*x + g/10 + rnormal()
    civreghdfe y (x=z), absorb(g) vce(cluster g)
    matrix V = e(V)
    foreach cl in frac large neg strcl {
        civreghdfe y (x=z), absorb(g) vce(cluster `cl')
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
    civreghdfe y (x=z), absorb(g) vce(cluster g h)
    matrix V = e(V)
    civreghdfe y (x=z), absorb(g) vce(cluster frac hfrac)
    mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    gen byte one_cluster = 1
    capture civreghdfe y (x=z), absorb(g) vce(cluster one_cluster)
    assert _rc == 459
    capture civreghdfe y (x=z), absorb(g) vce(cluster g one_cluster)
    assert _rc == 459
}
if _rc test_fail "S04 IV cluster labels" "r(`=_rc')"
else test_pass "S04 IV cluster labels"

capture noisily {
    clear
    set obs 402
    set seed 45680
    gen g = ceil(_n/20)
    replace g = 99 in 401
    replace g = 98 in 402
    gen double x = rnormal()
    gen double y = .4*x + g/10 + rnormal()
    gen double yp = rpoisson(exp(.3*x))
    creghdfe y x if _n <= 401, absorb(g) residuals(res)
    assert e(N) == 400
    count if e(sample)
    assert r(N) == 400
    assert e(sample) == (_n <= 400)
    assert !missing(res) == e(sample)
    estimates store ols_sample
    cpplmhdfe yp x if _n <= 401, absorb(g)
    assert e(N) == 400
    assert e(sample) == (_n <= 400)
    estimates restore ols_sample
    assert e(sample) == (_n <= 400)
}
if _rc test_fail "S05 final sample and residuals" "r(`=_rc')"
else test_pass "S05 final sample and residuals"

capture noisily {
    foreach zero in 1 10 20 {
        clear
        set obs 402
        set seed 45681
        gen g = ceil(_n/20)
        replace g = 99 in 401
        gen double cl = mod(_n,37)/100
        gen double x = rnormal()
        gen double off = rnormal()/10
        gen double w = 1+runiform()
        gen double y = rpoisson(exp(.3*x+off))
        replace y = 0 if g == `zero'
        cpplmhdfe y x if _n <= 401 [aw=w], absorb(g) vce(cluster cl) offset(off)
        assert e(N) == 380
        assert e(sample) == (_n <= 400 & g != `zero')
        matrix B = e(b)
        matrix V = e(V)
        local nc = e(N_clust)
        drop if _n > 400 | g == `zero'
        cpplmhdfe y x [aw=w], absorb(g) vce(cluster cl) offset(off)
        assert e(N_clust) == `nc'
        mata: assert(mreldif(st_matrix("B"),st_matrix("e(b)")) < 1e-7)
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
}
if _rc test_fail "S06 PPML separation and row alignment" "r(`=_rc')"
else test_pass "S06 PPML separation and row alignment"

capture noisily {
    clear
    set obs 200
    set seed 922
    gen g = ceil(_n/20)
    gen double x = rnormal()
    gen y = (_n < 110)
    gen byte one = 1
    capture cpplmhdfe y x, absorb(g) vce(cluster one)
    assert _rc == 459
    * A deliberately large guard threshold must fail, not count retained rows
    * as repeatedly dropped. A subsequent ordinary fit must remain usable.
    capture cpplmhdfe y x, absorb(g) septolerance(10)
    assert _rc == 430
    cpplmhdfe y x, absorb(g)
    assert e(num_separated) == 80
    assert e(sample) == (g <= 6)
}
if _rc test_fail "S05/S06 PPML guards and recovery" "r(`=_rc')"
else test_pass "S05/S06 PPML guards and recovery"

capture noisily {
    clear
    set obs 8
    gen double x = _n
    replace x = .a in 7
    replace x = .z in 8
    mata: st_vlmodify("lab", (1\2\3\4\5\6\.a\.z), ("plain"\(char(96)+"literal"+char(39))\("quote "+char(34)+" text")\("tabs"+char(9)+"line"+char(10)+char(13))\"café 文"\(3000*"x")\"refused"\"unknown"))
    label values x lab
    decode x, gen(native)
    cdecode x, gen(actual)
    assert r(N_vars) == 1
    assert actual == native
    assert strlen(actual[6]) == 3000
    cdecode x, gen(short) maxlength(4)
    decode x, gen(native_short) maxlength(4)
    assert short == native_short
    cdecode x, replace
    assert x == native
    * All variables are staged before replacing any input.
    gen y = 1
    label values y lab
    gen unlabeled = 2
    capture cdecode y unlabeled, replace
    assert _rc == 182
    confirm numeric variable y
}
if _rc test_fail "S07 literal labels and long strings" "r(`=_rc')"
else test_pass "S07 literal labels and long strings"

capture noisily {
    clear
    set obs 40
    gen cl = ceil(_n/10)
    gen str3 st = cond(cl<=2,"a","b")
    gen freq = 0
    cbsample 1, cluster(cl) weight(freq)
    summarize freq, meanonly
    assert r(sum) == 10
    cbsample 2, strata(cl) weight(freq)
    summarize freq, meanonly
    assert r(sum) == 8
    cbsample, strata(cl) weight(freq)
    summarize freq, meanonly
    assert r(sum) == 40
    cbsample 1, strata(st) cluster(cl) weight(freq)
    bysort st: egen drawn = total(freq)
    assert drawn == 10
    capture cbsample 3, strata(st) cluster(cl) weight(freq)
    assert _rc == 498
    assert _N == 40
    drop drawn
    * Repeated cluster labels in different strata are separate sampling units.
    replace cl = mod(cl,2)
    cbsample 1, strata(st) cluster(cl) weight(freq)
    bysort st: egen drawn = total(freq)
    assert drawn == 10
    cbsample 2, strata(st)
    assert _N == 4
}
if _rc test_fail "S08 bootstrap draw counts" "r(`=_rc')"
else test_pass "S08 bootstrap draw counts"

capture noisily {
    clear
    set obs 400
    set seed 45682
    gen g = ceil(_n/20)
    gen double x = rnormal()
    gen double z = x+rnormal()
    gen double y = .4*x + rnormal()
    gen double yp = rpoisson(exp(.3*x))
    gen double w = 1+runiform()
    gen double twice = 2*w
    foreach cmd in creghdfe cpplmhdfe civreghdfe {
        local dep = cond("`cmd'" == "cpplmhdfe", "yp", "y")
        local rhs = cond("`cmd'" == "civreghdfe", "(x=z)", "x")
        `cmd' `dep' `rhs' [aw=2*w], absorb(g)
        matrix B = e(b)
        matrix V = e(V)
        `cmd' `dep' `rhs' [aw=twice], absorb(g)
        mata: assert(mreldif(st_matrix("B"),st_matrix("e(b)")) < 1e-7)
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
    cpplmhdfe yp x, absorb(g)
    estimates store ppml_prediction
    capture predict double mu
    assert _rc == 198
    capture confirm variable mu
    assert _rc == 111
    predict double xb, xb
    assert reldif(xb, x*_b[x]) < 1e-7
    estimates restore ppml_prediction
    capture predict double mu
    assert _rc == 198
}
if _rc test_fail "S09/S10 predictions and weight expressions" "r(`=_rc')"
else test_pass "S09/S10 predictions and weight expressions"

capture noisily {
    clear
    input byte treat double(ps y)
    0 .5 1
    0 .5 3
    1 .5 7
    end
    cpsmatch treat, pscore(ps) outcome(y) ties
    assert _nn == 2 if treat
    assert missing(_nn) if !treat
    cpsmatch treat, pscore(ps) outcome(y) neighbor(5)
    assert _nn == 2 if treat
    cpsmatch treat, pscore(ps) outcome(y) radius caliper(.1)
    assert _nn == 2 if treat
    cpsmatch treat, pscore(ps) outcome(y) kernel bwidth(.1)
    assert _nn == 2 if treat
}
if _rc test_fail "S11 actual matching neighbor counts" "r(`=_rc')"
else test_pass "S11 actual matching neighbor counts"

capture noisily {
    tempfile empty full
    clear
    set obs 2
    gen id = _n
    gen str10 text = "example"
    save `full', replace
    drop in 1/2
    save `empty', replace
    foreach nogenerate in "" "nogenerate" {
        use `full', clear
        cmerge 1:1 id using `empty', keep(match) `nogenerate'
        assert _N == 0
        use `empty', clear
        cmerge 1:1 id using `full', keep(master) `nogenerate'
        assert _N == 0
        use `empty', clear
        cmerge 1:1 id using `full', keep(using) `nogenerate'
        assert _N == 2
        local typ : type text
        assert "`typ'" == "str10"
    }
    sysuse auto, clear
    capture cbinscatter price mpg, genxq(q) nograph
    assert _rc == 198
    capture confirm variable q
    assert _rc == 111
}
if _rc test_fail "S12/S15 empty-side filtering and unsupported option" "r(`=_rc')"
else test_pass "S12/S15 empty-side filtering and unsupported option"

print_summary sep22
if $TESTS_FAILED > 0 exit 9

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "sep22"
