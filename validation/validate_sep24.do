/* Regression contracts for the 21 SEP24_ASTRA_AUDIT findings. */
do validation/validate_setup.do

capture noisily {
    tempfile using_data
    clear
    set obs 2
    gen long id = 1000 + _n
    gen double value = 1000.125
    gen str20 text = "abcdefghijklmnopqrst"
    save `using_data', replace
    foreach option in "" "update" "update replace" {
        clear
        set obs 1
        gen byte id = 1
        gen byte value = .
        gen str1 text = ""
        cmerge 1:1 id using `using_data', `option'
        assert id == 1000 + _n - 1 if _n > 1
        assert value == 1000.125 if _merge == 2
        assert text == "abcdefghijklmnopqrst" if _merge == 2
        local vt : type value
        local st : type text
        assert "`vt'" == "double" & "`st'" == "str20"
    }
    clear
    set obs 1
    gen str4 id = "1001"
    gen row = _n
    quietly datasignature
    local before "`r(datasignature)'"
    foreach option in "" "force" {
        capture cmerge 1:1 id using `using_data', `option'
        assert _rc == 106
        quietly datasignature
        assert "`before'" == "`r(datasignature)'"
    }
}
if _rc test_fail "A01/A08 merge promotion and type safety" "r(`=_rc')"
else test_pass "A01/A08 merge promotion and type safety"

capture noisily {
    sysuse auto, clear
    quietly datasignature
    local before "`r(datasignature)'"
    capture creghdfe price mpg, residuals(mpg)
    assert _rc == 110
    capture crangestat (mean) mpg=mpg, interval(weight -100 100)
    assert _rc == 110
    capture crangestat (mean) out=mpg (sum) out=price, interval(weight -100 100)
    assert _rc == 198
    quietly datasignature
    assert "`before'" == "`r(datasignature)'"
    capture confirm variable out
    assert _rc == 111
}
if _rc test_fail "A02 output aliases and duplicates preserve data" "r(`=_rc')"
else test_pass "A02 output aliases and duplicates preserve data"

capture noisily {
    clear
    set obs 1000
    set seed 296519
    gen double x = rnormal()
    gen double z = x+rnormal()
    gen double z2 = 2*z
    gen double tiny_z = 1e-9*z
    gen double y = .3*x+4*z+rnormal()
    gen g = mod(_n,10)
    foreach absorb in "" "absorb(g)" {
        cbinscatter y x, controls(z) method(binsreg) nograph `absorb'
        matrix one = e(bindata)
        cbinscatter y x, controls(tiny_z) method(binsreg) nograph `absorb'
        matrix scaled = e(bindata)
        mata: assert(max(abs(st_matrix("one")-st_matrix("scaled"))) < 1e-7)
        cbinscatter y x, controls(z z2) method(binsreg) nograph `absorb'
        matrix dup = e(bindata)
        mata: assert(max(abs(st_matrix("one")-st_matrix("dup"))) < 1e-7)
    }
    gen double w = 1+abs(z)
    cbinscatter y x [aw=w], nograph
    matrix one = e(bindata)
    cbinscatter y x [aw=1+abs(z)], nograph
    matrix dup = e(bindata)
    mata: assert(max(abs(st_matrix("one")-st_matrix("dup"))) < 1e-10)
    gen byte category = mod(_n,4)
    quietly tab category, gen(dummy)
    cbinscatter y x, controls(dummy2 dummy3 dummy4) nograph
    matrix one = e(bindata)
    cbinscatter y x, controls(i.category) nograph
    matrix dup = e(bindata)
    mata: assert(max(abs(st_matrix("one")-st_matrix("dup"))) < 1e-7)
}
if _rc test_fail "A03/A16 redundant and factor controls, weight expressions" "r(`=_rc')"
else test_pass "A03/A16 redundant and factor controls, weight expressions"

capture noisily {
    clear
    set obs 1000
    set seed 742913
    gen g = ceil(_n/10)
    gen h = ceil(runiform()*50)
    gen double z = rnormal()
    gen double x = z+rnormal()+g/10+h/5
    gen double y = .7*x+g/5-h/7+rnormal()
    capture creghdfe y x, absorb(g h) iterate(1) tolerance(1e-14) residuals(res)
    assert _rc == 430
    capture confirm variable res
    assert _rc == 111
    capture civreghdfe y (x=z), absorb(g h) maxiter(1) tolerance(1e-14)
    assert _rc == 430
    foreach option in "iterate(0)" "tolerance(0)" {
        capture creghdfe y x, absorb(g h) `option'
        assert _rc == 198
    }
    creghdfe y x, absorb(g h) tolerance(1e-10)
    assert !missing(_b[x])
    replace y = rpoisson(exp(.3+.05*x+g/200+h/100))
    capture cpplmhdfe y x, absorb(g h) irlsmaxiter(1)
    assert _rc == 430
}
if _rc test_fail "A04 exhausted solvers return errors" "r(`=_rc')"
else test_pass "A04 exhausted solvers return errors"

capture noisily {
    foreach n in 999 1000 1001 2001 {
        clear
        set obs `n'
        gen double key = _n
        replace key = .a in `n'
        gen double v = 1
        replace v = 1001 in `n'
        foreach excl in "" "excludeself" {
            crangestat (mean) actual=v (sum) total=v (count) count=v, interval(key . .) `excl'
            assert actual == 1 if !missing(key)
            assert total == `n'-1-("`excl'"!="") if !missing(key)
            assert count == total if !missing(key)
            drop actual total count
        }
    }
}
if _rc test_fail "A07 missing range keys across thresholds" "r(`=_rc')"
else test_pass "A07 missing range keys across thresholds"

capture noisily {
    foreach n in 31 32 100 30000 {
        clear
        set obs `n'
        gen str8 key = cond(mod(_n,3)==0,"",cond(mod(_n,2),"az","ba"))
        gen long row = _n
        csort key, algorithm(merge) nosortedby threads(2)
        assert key >= key[_n-1] if _n > 1
        assert row > row[_n-1] if _n > 1 & key == key[_n-1]
    }
}
if _rc test_fail "A12 merge lexicographic order and stable ties" "r(`=_rc')"
else test_pass "A12 merge lexicographic order and stable ties"

capture noisily {
    clear
    set obs 600
    set seed 92815
    gen g = ceil(_n/10)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double zero = 0
    gen double y = rpoisson(exp(.3+.4*x+.2*z+g/200))
    gen byte small_cluster = 1+mod(g,2)
    foreach vc in "vce(robust)" "vce(cluster g)" "vce(cluster small_cluster)" {
    foreach vars in "zero x z" "x zero z" "x z zero" {
        cpplmhdfe y `vars', absorb(g) `vc'
        scalar reported = e(F)
        test x z
        assert reldif(reported,r(F)) < 1e-10
        assert e(df_m) == r(df)
    }
    }
}
if _rc test_fail "A14 PPML joint tests select retained regressors" "r(`=_rc')"
else test_pass "A14 PPML joint tests select retained regressors"

capture noisily {
    clear
    set obs 1
    gen str1 a = "x"
    gen strL b = "z"
    forvalues i=1/12 {
        quietly replace b=b+b
    }
    quietly datasignature
    local before "`r(datasignature)'"
    capture cencode a b, replace
    assert _rc != 0
    quietly datasignature
    assert "`before'" == "`r(datasignature)'"
    label define first 1 "old"
    capture cencode a b, gen(first second)
    assert _rc != 0
    local old : label first 1
    local added : label first 2, strict
    assert "`old'" == "old" & "`added'" == ""
    gen byte target=9
    drop in 1
    capture cencode a, gen(target)
    assert _rc == 110
    local vt : type target
    assert "`vt'" == "byte"
    capture cencode a b, gen(new new)
    assert _rc == 198
    capture cencode target, gen(new)
    assert _rc == 107
    cencode a, gen(new)
    confirm numeric variable new
}
if _rc test_fail "A15 encoding failure and empty data contracts" "r(`=_rc')"
else test_pass "A15 encoding failure and empty data contracts"

capture noisily {
    clear
    set obs 3
    gen double day = _n-1
    format day %td
    gen double timestamp = clock("24sep2026 12:34:56.789","DMYhms")+_n
    format timestamp %tc
    gen str10 note = "hello " + string(_n)
    foreach ext in csv xlsx {
        local path "temp/space threads(2) name.`ext'"
        if "`ext'" == "csv" cexport delimited using "`path'", replace
        else cexport excel using "`path'", firstrow(variables) sheet("A sheet") replace
        confirm file "`path'"
        if "`ext'" == "csv" {
            capture cexport delimited using "`path'"
        }
        else {
            capture cexport excel using "`path'"
        }
        assert _rc == 602
    }
    cimport excel using "temp/space threads(2) name.xlsx", firstrow clear
    assert day == _n-1
    assert note == "hello " + string(_n)
    assert abs(timestamp - (clock("24sep2026 12:34:56.789","DMYhms")+_n)) < .1
    cimport excel using "validation/fixtures/sep24_dates1900.xlsx", firstrow clear
    assert day == td(28feb1900) in 1
    assert missing(day) in 2
    assert day == td(01mar1900) in 3
    assert day == 0 in 4
    cimport excel using "validation/fixtures/sep24_dates1904.xlsx", firstrow clear
    assert day == td(01jan1904) in 1
    assert day == 0 in 2
}
if _rc test_fail "A05/A06 dates, datetimes, literal export paths" "r(`=_rc')"
else test_pass "A05/A06 dates, datetimes, literal export paths"

capture noisily {
    clear
    set obs 1
    forvalues i=1/1100 {
        gen double long_variable_name_for_test_`i' = `i'+.125
    }
    cexport delimited using "temp/sep24_wide.csv", replace
    cexport excel using "temp/sep24_wide.xlsx", firstrow(variables) replace
    import delimited using "temp/sep24_wide.csv", clear varnames(1)
    assert c(k) == 1100
    forvalues i=1/1100 {
        assert long_variable_name_for_test_`i' == `i'+.125
    }
    import excel using "temp/sep24_wide.xlsx", clear firstrow
    assert c(k) == 1100
    forvalues i=1/1100 {
        assert long_variable_name_for_test_`i' == `i'+.125
    }
}
if _rc test_fail "A13 complete wide export metadata" "r(`=_rc')"
else test_pass "A13 complete wide export metadata"

capture noisily {
    _ctools_load
    assert "`r(plugin_version)'" != ""
    local revision "`r(plugin_revision)'"
    _ctools_load
    assert "`r(plugin_revision)'" == "`revision'"
}
if _rc test_fail "A20 repeated shared plugin loading" "r(`=_rc')"
else test_pass "A20 repeated shared plugin loading"

print_summary sep24
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "sep24"
