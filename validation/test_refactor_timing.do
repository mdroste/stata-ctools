adopath ++ build
clear
set obs 200
gen double x = mod(_n, 17)
gen double y = _n
gen int g = mod(_n, 4)
timer clear 99
timer on 99
quietly summarize y
local mean = r(mean)
quietly _ctools_clock now
assert r(mean) == `mean'
assert `now' > 0
csort x, verbose
creghdfe y x, absorb(g) verbose
timer off 99
quietly timer list 99
assert r(nt99) == 1
assert r(t99) > 0
local prior = r(t99)
csort y, verbose
quietly timer list 99
assert r(nt99) == 1
assert r(t99) == `prior'
di "CTOOLS_REFACTOR_TIMING_COMPLETE"
