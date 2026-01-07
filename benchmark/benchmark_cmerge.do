* Add ctools to adopath
adopath + "/Users/Mike/Documents/GitHub/stata-ctools/build"

*-------------------------------------------------------------
* Simulated data speed benchmark
*-------------------------------------------------------------

clear all
set obs 25000000
forval i=1/10 {
	gen x`i' = rnormal()
}
gen key = ceil(runiform()*10000)
tempfile t1
save `t1'
clear all
set obs 10000
gen key = _n
gen y = rnormal()
gen x = rnormal()
sort x
tempfile t2
save `t2'

* Benchmark cmerge and merge for m:1 merge
di as text "{hline 60}"
di as text "Benchmark: Stata merge vs cmerge"
di as text "{hline 60}"

* Benchmark merge
use `t1', clear
timer clear 1
timer on 1
merge m:1 key using `t2'
timer off 1
tempfile m1
save `m1'
drop _merge

* Display merge benchmark result
qui timer list 1
local stata_time = r(t1)
di as text "Stata merge time: " as result %9.3f `stata_time' as text " seconds"

* Benchmark cmerge
use `t1', clear
timer clear 2
timer on 2
cmerge m:1 key using `t2'
timer off 2
tempfile m2
save `m2'

* Display cmerge benchmark result
qui timer list 2
local cmerge_time = r(t2)
di as text "cmerge time:      " as result %9.3f `cmerge_time' as text " seconds"

* Summary
di as text "{hline 60}"
local speedup = `stata_time' / `cmerge_time'
di as text "Speedup: " as result %5.2f `speedup' as text "x"
di as text "{hline 60}"

* Checking if data is the same
use `m1'
cf _all using `m2'
