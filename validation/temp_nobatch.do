clear all
set obs 25000000
forval i=1/10 {
	gen x`i' = rnormal()
}

sort x1
csort x2, verbose
csort x3, verbose
global ctools_nobatch 1
csort x4, verbose
