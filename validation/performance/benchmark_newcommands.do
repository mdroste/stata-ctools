* Run via a Stata driver that captures errors and exits cleanly.
* Reference rangejoin must be on adopath. Timings are end-to-end, one run each.
clear all
set more off
adopath ++ "build"
ctools, version
creturn list
which ipolate
which split
which rangejoin
local N = 1000000
set seed 95183
clear
set obs `N'
gen long id=_n
gen int g=mod(_n,1000)
gen double x=floor((_n-1)/1000)
gen double y=rnormal()
replace y=. if runiform()<.6
tempfile base
save `base'
timer clear
timer on 1
ipolate y x, gen(ref) by(g) epolate
timer off 1
sort id
tempfile interpolated
save `interpolated'
use `base', clear
timer on 2
cipolate y x, gen(ref) by(g) epolate verbose
timer off 2
cf _all using `interpolated'
clear
set obs `N'
gen str80 text="NY|123|business|active|2024-01-31"
timer on 3
split text, gen(ref) parse("|")
timer off 3
timer on 4
csplit text, gen(got) parse("|") verbose
timer off 4
forvalues j=1/5 {
    assert ref`j'==got`j'
}
clear
set obs 200000
gen long id=_n
gen long g=floor((_n-1)/200)
gen double key=mod(_n-1,200)
gen str20 text="row"+string(_n)
tempfile using master joined
save `using'
replace key=key+.5
save `master'
timer on 5
rangejoin key -1 1 using `using', by(g)
timer off 5
save `joined'
use `master', clear
timer on 6
crangejoin key -1 1 using `using', by(g) verbose
timer off 6
cf _all using `joined'
quietly timer list 1
local ipolate=r(t1)
quietly timer list 2
local cipolate=r(t2)
quietly timer list 3
local split=r(t3)
quietly timer list 4
local csplit=r(t4)
quietly timer list 5
local rangejoin=r(t5)
quietly timer list 6
local crangejoin=r(t6)
di "BENCH ipolate rows=`N' native=`ipolate' ctools=`cipolate'"
di "BENCH split rows=`N' native=`split' ctools=`csplit'"
di "BENCH rangejoin rows=200000+200000 output=`=_N' native=`rangejoin' ctools=`crangejoin'"
