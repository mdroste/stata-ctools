* Test cmerge with string key - partial match case
clear all
set more off
adopath + build

* Create master dataset (fewer cars)
sysuse auto, clear
keep if foreign == 0   // Domestic only = 52 cars
keep make price mpg
sort make
tempfile master
save `master'

* Create using dataset (all cars)
sysuse auto, clear
keep make weight length
sort make
tempfile using
save `using'

* Load master and merge
use `master', clear
di "Master observations (domestic only): " _N

* Merge with keep(match using)
cmerge 1:1 make using `using', verbose keep(match using)

* Check results
di ""
di "After merge observations: " _N
di ""
tab _merge

di ""
di "First 10 observations:"
list make price mpg weight length _merge in 1/10

di ""
di "Observations from using only (should be foreign cars):"
list make weight length _merge if _merge == 2

* Verify no duplicates in make
duplicates report make
