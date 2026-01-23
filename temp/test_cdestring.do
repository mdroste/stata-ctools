capture do "validation/validate_setup.do"
clear
set obs 3
gen str5 x = ""
replace x = "1" in 1
replace x = "2" in 2
replace x = "3" in 3
list

cdestring x, generate(y)
list
