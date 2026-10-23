* Exact winsor2 cutoff expansion, output storage, selection, and phase cleanup.
version 16
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
local cases=0
foreach storage in float double {
    foreach sign in -1 1 {
        foreach scale in 1 1e-20 1e20 {
            clear
            set obs 204
            gen `storage' x=0
            replace x=-3 in 1
            replace x=-2.3206636905670166 in 2
            replace x=-2.3202905654907227 in 3
            replace x=`sign'*`scale'*x
            replace x=. in 201
            replace x=.a in 202
            replace x=.z in 203
            replace x=.b in 204
            gen double y=1.2345678901234567*x
            gen long id=_n
            foreach sample in "" "if id<=200 in 1/203" {
                foreach trim in "" "trim" {
                    foreach cuts in "1 99" "0 99" "1 100" {
                        quietly winsor2 x y `sample', cuts(`cuts') `trim' suffix(_expected)
                        foreach threads in 1 8 {
                            quietly cwinsor x y `sample', cuts(`cuts') `trim' suffix(_actual) threads(`threads')
                            assert x_actual==x_expected & y_actual==y_expected
                            assert "`: type x_actual'"=="`: type x_expected'"
                            drop x_actual y_actual
                            local ++cases
                        }
                        drop x_expected y_expected
                        preserve
                        quietly winsor2 x y `sample', cuts(`cuts') `trim' replace
                        rename x x_expected
                        rename y y_expected
                        tempfile expected
                        quietly save `expected', replace
                        restore
                        preserve
                        quietly cwinsor x y `sample', cuts(`cuts') `trim' replace
                        quietly merge 1:1 id using `expected', nogen
                        assert x==x_expected & y==y_expected
                        restore
                        local ++cases
                    }
                }
            }
        }
    }
}
* A single by-group must use grouped semantics, without numeric macro rounding.
clear
set obs 200
gen float x=0
replace x=-3 in 1
replace x=-2.3206636905670166 in 2
replace x=-2.3202905654907227 in 3
gen g=1
winsor2 x, by(g) suffix(_expected)
cwinsor x, by(g) suffix(_actual)
assert x_actual==x_expected
local ++cases
noi di "WINSOR_ROUNDING_CASES=`cases'"
noi di "WINSOR_ROUNDING_COMPLETE"
