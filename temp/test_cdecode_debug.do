capture do "validation/validate_setup.do"
clear
set obs 3
gen int x = _n
label define xlbl 1 "A" 2 "B" 3 "C"
label values x xlbl
list

cdecode x, generate(y)
list
