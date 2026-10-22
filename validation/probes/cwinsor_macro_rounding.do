* Reproduces a pre-existing one-float-ULP difference; not a passing parity test.
* Both the pre-audit and optimized ctools builds give the same result.
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
assert x_c[1]!=x_ref[1]
assert abs(x_c[1]-x_ref[1])==2^-22
di "CWINSOR_MACRO_ROUNDING_REPRODUCED"
