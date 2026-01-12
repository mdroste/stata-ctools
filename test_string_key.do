* Test cmerge with string key (auto dataset)
clear all
set more off
adopath + build

* Load auto dataset and save master
sysuse auto, clear
keep make price mpg
sort make
tempfile master
save `master'

* Create using dataset with different variables
sysuse auto, clear
keep make weight length
sort make
tempfile using
save `using'

* Load master and merge
use `master', clear
di "Master observations: " _N
list make in 1/10

* Merge
cmerge 1:1 make using `using', verbose

* Check results
di ""
di "After merge observations: " _N
di ""
di "Checking make variable for duplicates/corruption:"
tab _merge
di ""
di "First 20 observations:"
list make price mpg weight length _merge in 1/20

* Verify no duplicates in make
duplicates report make
