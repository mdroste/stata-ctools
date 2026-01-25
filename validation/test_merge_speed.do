*-------------------------------------------------------------------------------
* Speed benchmark, 1:1 and m:1 merge vs. cmerge, large dataset, one key
*-------------------------------------------------------------------------------


local N 25000000
local K1 10
local K2 10
local G1 10000

clear
set obs `N'
forval i=1/`K1' {
	gen x`i' = rnormal()
}
gen long i = _n
gen g = ceil(runiform()*`G1')
sort x1
tempfile t1
save `t1'
sort i
tempfile t1s
save `t1s'

clear
set obs `N'
forval i=1/`K2' {
	gen y`i'=rnormal()
}
gen long i = _n
sort y1
tempfile t2
save `t2'
sort i
tempfile t2s
save `t2s'

clear
set obs `G1'
gen g = _n
forval i=1/10 {
	gen z`i' = rnormal()
}
tempfile t3
save `t3'

* Test 1:1 merge, unsorted datasets
use `t1', clear
merge 1:1 i using `t2'
use `t1', clear
cmerge 1:1 i using `t2', verbose

* Test 1:1 merge, sorted datasets
use `t1s', clear
merge 1:1 i using `t2s'
use `t1s', clear
cmerge 1:1 i using `t2s', verbose

* Test m:1 merge
use `t1', clear
merge m:1 g using `t3'
use `t1', clear
cmerge m:1 g using `t3', verbose
