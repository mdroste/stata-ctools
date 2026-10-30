* Manual diagnostics; sizes and settings are deliberately unchanged.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}

if "`__workload'" == "cppml_pathology_debug" {
do validation/support/validate_setup.do

capture program drop cmp_case
program define cmp_case
    syntax, name(string)

    matrix pV = e(V)
    matrix pB = e(b)

    capture noisily cppmlhdfe `name'
    matrix cV = e(V)
    matrix cB = e(b)

    local pcols = colsof(pB)
    local ccols = colsof(cB)
    local K = min(`pcols', `ccols')

    local min_sf_b = 15
    forvalues j = 1/`K' {
        sigfigs pB[1,`j'] cB[1,`j']
        if r(sigfigs) < `min_sf_b' local min_sf_b = r(sigfigs)
    }

    local min_sf_v = 15
    forvalues i = 1/`K' {
        forvalues j = 1/`K' {
            sigfigs pV[`i',`j'] cV[`i',`j']
            if r(sigfigs) < `min_sf_v' local min_sf_v = r(sigfigs)
        }
    }

    di as txt "case=`name' min_sf_b=" as res %6.3f `min_sf_b' as txt " min_sf_v=" as res %6.3f `min_sf_v'
end

* Case 1: zero-inflated Y
clear
set seed 222
set obs 2000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(-2 + 0.3 * x)
gen y = rpoisson(mu)
quietly ppmlhdfe y x, absorb(fe1) vce(robust)

cmp_case, name("y x, absorb(fe1) vce(robust)")
cmp_case, name("y x, absorb(fe1) vce(robust) irlstolerance(1e-14) tolerance(1e-12)")
cmp_case, name("y x, absorb(fe1) vce(robust) irlstolerance(1e-16) tolerance(1e-14)")

* Case 2: mixed covariate scales
clear
set seed 223
set obs 2000
gen fe1 = mod(_n-1, 50) + 1
gen x_tiny = rnormal() / 1e6
gen x_huge = rnormal() * 1e6
gen mu = exp(1e6 * x_tiny + 1e-6 * x_huge)
gen y = rpoisson(mu)
quietly ppmlhdfe y x_tiny x_huge, absorb(fe1) vce(robust)

cmp_case, name("y x_tiny x_huge, absorb(fe1) vce(robust)")
cmp_case, name("y x_tiny x_huge, absorb(fe1) vce(robust) irlstolerance(1e-14) tolerance(1e-12)")
cmp_case, name("y x_tiny x_huge, absorb(fe1) vce(robust) irlstolerance(1e-16) tolerance(1e-14)")
exit
}

if "`__workload'" == "cppml_pathology_minchange" {
do validation/support/validate_setup.do

capture program drop cmp_case2
program define cmp_case2
    syntax, cmd(string)
    matrix pV = e(V)
    matrix pB = e(b)
    capture noisily cppmlhdfe `cmd'
    matrix cV = e(V)
    matrix cB = e(b)
    local K = min(colsof(pB), colsof(cB))
    local min_sf_v = 15
    forvalues i=1/`K' {
        forvalues j=1/`K' {
            sigfigs pV[`i',`j'] cV[`i',`j']
            if r(sigfigs) < `min_sf_v' local min_sf_v = r(sigfigs)
        }
    }
    di as txt "cmd=" as txt "`cmd'" as txt " min_sf_v=" as res %6.3f `min_sf_v'
end

* zero-inflated
clear
set seed 222
set obs 2000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(-2 + 0.3 * x)
gen y = rpoisson(mu)
qui ppmlhdfe y x, absorb(fe1) vce(robust)
cmp_case2, cmd("y x, absorb(fe1) vce(robust)")
cmp_case2, cmd("y x, absorb(fe1) vce(robust) tolerance(1e-12)")
cmp_case2, cmd("y x, absorb(fe1) vce(robust) irlstolerance(1e-14)")
cmp_case2, cmd("y x, absorb(fe1) vce(robust) tolerance(1e-12) irlstolerance(1e-14)")

* mixed scale
clear
set seed 223
set obs 2000
gen fe1 = mod(_n-1, 50) + 1
gen x_tiny = rnormal() / 1e6
gen x_huge = rnormal() * 1e6
gen mu = exp(1e6 * x_tiny + 1e-6 * x_huge)
gen y = rpoisson(mu)
qui ppmlhdfe y x_tiny x_huge, absorb(fe1) vce(robust)
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust)")
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust) tolerance(1e-12)")
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust) irlstolerance(1e-14)")
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust) tolerance(1e-12) irlstolerance(1e-14)")
exit
}

if "`__workload'" == "cwinsor_macro_rounding" {
* Regression for the former one-float-ULP cutoff rounding difference.
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
set obs 200
gen float x=0
replace x=-3 in 1
replace x=-2.3206636905670166 in 2
replace x=-2.3202905654907227 in 3
_pctile x, p(1 99)
scalar exact = r(r1)
local through_macro = scalar(exact)
di "exact cutoff: " %24.17g scalar(exact)
di "macro cutoff: " %24.17g `through_macro'
cwinsor x, suffix(_c)
winsor2 x, suffix(_ref)
format x x_c x_ref %24.17g
list in 1/3
assert x_c==x_ref
di "CWINSOR_MACRO_ROUNDING_FIXED"
exit
}

if "`__workload'" == "merge" {
* single variable merge key

* the NOBATCH key is important!
* global ctools_nobatch 1
clear all
set obs 25000000
gen long i = _n
forval i=1/10 {
	gen x`i' = rnormal()
}
sort x1
tempfile t1
save `t1'

clear all
set obs 25000000
gen long i = _n
forval i=1/10 {
	gen y`i' = rnormal()
}
sort y1
tempfile t2
save `t2'

use `t1', clear
merge 1:1 i using `t2'


use `t1', clear
cmerge 1:1 i using `t2', verbose
exit
}

if "`__workload'" == "nobatch" {
clear all
set obs 25000000
forval i=1/10 {
	gen x`i' = rnormal()
}

sort x1
csort x2, verbose
csort x3, verbose
global ctools_nobatch 1
csort x4, verbose
exit
}

if "`__workload'" == "which_cppml" {
do validation/support/validate_setup.do
which cppmlhdfe
type "build/cppmlhdfe.ado", lines(35)
exit
}

if "`__workload'" == "settings" {
clear all
adopath + build
set obs 10000000
forval i=1/10 {
	gen x`i' = rnormal()
}

sort x1
csort x2, verbose
csort x3, verbose
exit
}

exit 198
