* Debug factor variable base level handling
capture adopath ++ "build"

sysuse auto, clear

* Show what fvexpand gives for a simple interaction
di "=== fvexpand c.mpg#i.foreign ==="
fvexpand c.mpg#i.foreign
di "`r(varlist)'"

* Show what fvrevar gives
di "=== fvrevar c.mpg#i.foreign ==="
fvrevar c.mpg#i.foreign
di "`r(varlist)'"

* Run reghdfe and show e(b) column names
di "=== reghdfe ==="
quietly reghdfe price c.mpg#i.foreign weight, absorb(rep78)
matrix list e(b)
di "df_m = " e(df_m) " rank = " e(rank)

* Run creghdfe and show e(b) column names
di "=== creghdfe ==="
quietly creghdfe price c.mpg#i.foreign weight, absorb(rep78)
matrix list e(b)
di "df_m = " e(df_m) " rank = " e(rank)
