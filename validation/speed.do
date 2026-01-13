*==============================================================================
* Benchmarking times for Stata
*==============================================================================

*------------------------------------------------------------------------------
* Benchmark csort
*------------------------------------------------------------------------------

timer clear
clear
set obs 50000000
forval i=1/8 {
	gen f`i' = rnormal()
}
forval i=1/2 {
	gen g`i' = ceil(runiform()*1000)
}
gen s1 = char(runiformint(65,90)) + string(runiformint(0,9)) ///
		 + char(runiformint(65,90)) + char(runiformint(65,90)) ///
         + string(runiformint(0,9)) + char(runiformint(65,90))
gen y = 0+f1+2*f2+3*f3+4*f4+5*f5+rnormal()

* Float sort comparison
timer on 1
sort f1
timer off 1
timer on 2
csort f2
timer off 2

* Integer sort comparison
timer on 3
sort g1
timer off 3
timer on 4
csort g2
timer off 4

* String sort comparison
timer on 5
sort s1
timer off 5
csort f1
timer on 6
csort s1
timer off 6


*------------------------------------------------------------------------------
* Benchmark cmerge
*------------------------------------------------------------------------------

* Create data for 1:1 merge
gen long id = _n
csort f1
preserve
keep id f*
tempfile t1
save `t1'
restore, preserve
keep id g1 g2 s1 y
tempfile t2
save `t2'

* Merge 1:1 test
use `t1', clear
timer on 7
merge 1:1 id using `t2'
timer off 7
restore, preserve
timer on 8
cmerge 1:1 id using `t2'
timer off 8
restore, preserve

* Merge m:1 test data setup
tempfile t3
save `t3'
clear
set obs 1000
gen g1 = _n
gen z = rnormal()
tempfile t4
save `t4'
use `t3', clear
timer on 9
merge m:1 g1 using `t4'
timer off 9
use `t3', clear
timer on 10
cmerge m:1 g1 using `t4'
timer off 10

* Restore data now
restore

*------------------------------------------------------------------------------
* Benchmark cqreg
*------------------------------------------------------------------------------

/*
timer on 11
qreg y f1 f2 f3 f4 f5
timer off 11
timer on 12
cqreg y f1 f2 f3 f4 f5
timer off 12
*/

*------------------------------------------------------------------------------
* Benchmark creghdfe
*------------------------------------------------------------------------------


timer on 13
reghdfe y f1 f2 f3 f4 f5, absorb(g1 g2)
timer off 13
timer on 14
creghdfe y f1 f2 f3 f4 f5, absorb(g1 g2)
timer off 14

*------------------------------------------------------------------------------
* Benchmark cbinscatter
*------------------------------------------------------------------------------


* Binscatter without controls
timer on 15
binscatter y f1
timer off 15
timer on 16
cbinscatter y f1
timer off 16

* Binscatter with controls
timer on 17
binscatter y f1, controls(f2 f3 f4) absorb(g1)
timer off 17
timer on 18
cbinscatter y f1, controls(f2 f3 f4) absorb(g1)
timer off 18

timer list
timer clear
clear
set obs 18
gen n1 = _n
gen ctools = mod(n1,2)==0
gen time=.
forval i=1/18 {
replace time = r(t`i') in `i'
}
generate str operation = "sort float" in 1
replace operation = "sort float" in 2
replace operation = "sort int" in 3
replace operation = "sort int" in 4
replace operation = "sort str" in 5
replace operation = "sort str" in 6
replace operation = "merge 1:1" in 7
replace operation = "merge 1:1" in 8
replace operation = "merge m:1" in 9
replace operation = "merge m:1" in 10
replace operation = "qreg" in 11
replace operation = "qreg" in 12
replace operation = "reghdfe" in 13
replace operation = "reghdfe" in 14
replace operation = "binscatter" in 15
replace operation = "binscatter" in 16
replace operation = "binscatter w/ controls" in 17
replace operation = "binscatter w/ controls" in 18

encode operation, gen(op)

twoway (bar time op if ctools==0, color(maroon) barwidth(0.4)) (bar time op if ctools==1, color(ebblue) barwidth(0.4))
gen time1 = time if ctools==1
gen time2 = time if ctools==0
graph hbar time1 time2, over(operation)  legend( ring(0) pos(2) label(1 "ctools") label(2 "Stata default") )
