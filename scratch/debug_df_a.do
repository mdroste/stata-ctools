* Debug df_a calculation
clear all
set more off
adopath + "../build"

sysuse auto

* Test 4: cluster by foreign (which is also an absorbed FE)
di _n "=== reghdfe with cluster(foreign) ==="
reghdfe price mpg, absorb(foreign trunk) vce(cluster foreign)
di "e(df_a) = " e(df_a)
di "e(df_a_nested) = " e(df_a_nested)
di "e(df_a_initial) = " e(df_a_initial)
di "e(df_a_redundant) = " e(df_a_redundant)
ereturn list

di _n "=== creghdfe with cluster(foreign) ==="
creghdfe price mpg, absorb(foreign trunk) vce(cluster foreign)
ereturn list
