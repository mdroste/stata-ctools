* Debug: check exact interaction behavior for failing tests
capture adopath ++ "build"

* Test 1: i.race#i.married from nlsw88
sysuse nlsw88, clear
di "=== reghdfe: i.race#i.married ==="
quietly reghdfe wage i.race#i.married age, absorb(industry)
matrix list e(b)
di "df_m = " e(df_m) " rank = " e(rank)

di ""
di "=== fvexpand i.race#i.married ==="
fvexpand i.race#i.married
di "`r(varlist)'"

di ""
di "=== creghdfe: i.race#i.married ==="
quietly creghdfe wage i.race#i.married age, absorb(industry)
matrix list e(b)
di "df_m = " e(df_m) " rank = " e(rank)

* Test 2: c.mpg#i.foreign from auto
sysuse auto, clear
di ""
di "=== reghdfe: c.mpg#i.foreign ==="
quietly reghdfe price c.mpg#i.foreign weight, absorb(rep78)
matrix list e(b)
di "df_m = " e(df_m) " rank = " e(rank)

di ""
di "=== creghdfe: c.mpg#i.foreign ==="
quietly creghdfe price c.mpg#i.foreign weight, absorb(rep78)
matrix list e(b)
di "df_m = " e(df_m) " rank = " e(rank)
