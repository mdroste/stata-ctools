* Debug: verify reghdfe behavior for different interaction types
capture adopath ++ "build"

* Pure indicator interaction
webuse nlswork, clear
quietly reghdfe ln_wage i.race#i.collgrad age, absorb(ind_code)
di "=== i.race#i.collgrad ==="
matrix list e(b)
di "df_m = " e(df_m)

* Continuous x indicator interaction
sysuse auto, clear
quietly reghdfe price c.mpg#i.foreign weight, absorb(rep78)
di "=== c.mpg#i.foreign ==="
matrix list e(b)
di "df_m = " e(df_m)
