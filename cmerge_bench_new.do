clear all
adopath + build
set obs 50000000
gen i = _n
forval i=1/10 {
gen x`i' = rnormal()
}
gen g = ceil(runiform()*10000000)
tempfile t1
save `t1'
clear
set obs 10000000
gen g = _n
forval i=1/10 {
gen y`i' = rnormal()
}
sort y1
tempfile t2
save `t2'
use `t1', clear
merge m:1 g using `t2'
use `t1', clear
cmerge m:1 g using `t2', verbose
