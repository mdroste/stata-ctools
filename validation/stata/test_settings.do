clear all
adopath + build
set obs 10000000
forval i=1/10 {
	gen x`i' = rnormal()
}

sort x1
csort x2, verbose
csort x3, verbose
