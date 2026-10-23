* Sample selection and complete destination writes for the October 2026 fast paths.
version 16
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
set seed 72391
set obs 20000
gen long id = _n
gen str8 g = "g" + string(mod(id,13))
gen double x = floor(id/13)
gen double y = 3*x + mod(id,13)
replace y = . if mod(id,5)==0
replace x = .a if mod(id,11)==0
replace g = "" if mod(id,31)==0
gen str80 text = cond(mod(id,7)==0,"", cond(mod(id,3)==0,"alpha", ///
    cond(mod(id,3)==1,"  alpha|beta  ","alpha||東京|tail")))
forvalues j=1/8 {
    gen double v`j' = rnormal()
    gen str24 s`j' = "unused payload" + string(mod(id,100))
}
local comparisons=0
foreach selection in "" "if mod(id,3)==0" "in 13/17000" ///
    "if mod(id,3)==1 in 13/17000" "if 0" "if missing(x)" {
    foreach grouping in "" "by(g)" "by(x)" "by(y)" {
        ipolate y x `selection', gen(expected) epolate `grouping'
        cipolate y x `selection', gen(actual) epolate `grouping'
        assert reldif(actual,expected)<1e-12 | (missing(actual)&missing(expected))
        assert id==_n
        drop expected actual
        local ++comparisons
    }
    capture split text `selection', gen(expected) parse("|")
    local expected_rc=_rc
    capture csplit text `selection', gen(actual) parse("|")
    assert _rc==`expected_rc'
    if !`expected_rc' {
        unab expected_vars : expected*
        unab actual_vars : actual*
        local fields : word count `expected_vars'
        assert `fields'==`: word count `actual_vars''
        forvalues j=1/`fields' {
            assert actual`j'==expected`j'
            assert "`: type actual`j''"=="`: type expected`j''"
        }
        drop actual* expected*
    }
    assert id==_n
    local ++comparisons
}
* Empty strings are excluded; whitespace-only strings are valid split inputs.
replace text=""
capture csplit text, gen(actual)
assert _rc==2000
capture confirm variable actual1
assert _rc==111
replace text="   "
split text, gen(expected)
csplit text, gen(actual)
assert actual1==expected1
drop expected* actual*
local comparisons=`comparisons'+2
* Every missing x value must be written even with an uninitialized destination.
replace x=.z
cipolate y x, gen(actual) epolate
assert actual==.
drop actual
local ++comparisons

* Grouped winsor output precision, metadata and missing values match winsor2.
* No-group comparisons use the explicit order-statistic formula, since native
* winsor2 rounds its percentile through a macro (a pre-existing difference).
gen byte vb=mod(id,103)
gen int vi=mod(id*37,20000)
replace vb=.z in 7
replace vi=.a in 9
foreach storage in float double {
    set type `storage'
    foreach selection in "" "if mod(id,3)==0" "in 13/17000" "if id>50 in 13/17000" {
        replace v1=.a in 3
        replace v2=.z in 4
        label variable v1 "first source"
        format v1 %15.4f
        cwinsor v1 v2 vb vi, by(g) suffix(_all)
        winsor2 v1 v2 vb vi, by(g) suffix(_ref)
        foreach v in v1 v2 vb vi {
            assert `v'_all==`v'_ref
            drop `v'_all `v'_ref
        }
        cwinsor v1 v2 vb vi `selection', by(g) suffix(_actual)
        winsor2 v1 v2 vb vi `selection', by(g) suffix(_expected)
        foreach v in v1 v2 vb vi {
            assert `v'_actual==`v'_expected
            assert "`: type `v'_actual'"=="`: type `v'_expected'"
            assert "`: format `v'_actual'"=="`: format `v'_expected'"
            assert `"`: variable label `v'_actual'"'==`"`: variable label `v'_expected'"'
        }
        drop v1_actual v2_actual vb_actual vi_actual v1_expected v2_expected vb_expected vi_expected
        local comparisons=`comparisons'+2
    }
}
set type float
* All outputs belong to the wrapper transaction: direct-store errors must undo them.
foreach fault in outputs store post {
    quietly datasignature
    local before `r(datasignature)'
    global CTOOLS_cwinsor_fault `fault'
    capture cwinsor v1 v2, by(g) suffix(_failure)
    local rc=_rc
    global CTOOLS_cwinsor_fault
    assert `rc'==cond("`fault'"=="store",920,1)
    quietly datasignature
    assert "`before'"=="`r(datasignature)'"
    local ++comparisons
}
gen float vf=v2
foreach selection in "" "if mod(id,3)==0" "in 13/17000" {
    gen byte selected=0
    replace selected=1 `selection'
    foreach source in v1 vf {
        _pctile `source' if selected, p(5 95)
        * winsor2 expands numeric locals, rounding the ungrouped cutoffs.
        local lower=r(r1)
        local upper=r(r2)
        local storage : type `source'
        gen `storage' expected=min(max(`source',`lower'),`upper') ///
            if selected & !missing(`source')
        cwinsor `source' `selection', cuts(5 95) suffix(_direct)
        assert `source'_direct==expected
        drop `source'_direct expected
        local ++comparisons
    }
    drop selected
}
* Missing-x rows must be excluded before attempting SPI reads of long by-keys.
clear
set obs 3
gen double x=_n
gen double y=2*x
replace x=.z in 2
gen strL longkey="g"
forvalues j=1/12 {
    quietly replace longkey=longkey+longkey in 2
}
cipolate y x, by(longkey) gen(actual) epolate
assert actual==y if !missing(x)
assert actual==. if missing(x)
drop actual
replace x=.
cipolate y x, by(longkey) gen(actual)
assert actual==.
local comparisons=`comparisons'+2
di "PERFORMANCE_PATHS_COMPARISONS=`comparisons'"
di "PERFORMANCE_PATHS_COMPLETE"
