capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

* Create a test scenario where we know the order
clear
set obs 3
gen byte first = 1
gen str10 x = "A"
gen byte third = 3

di "Dataset variables before cencode:"
describe

* cencode will create y. At that point:
* first=1, x=2, third=3, y=4 in the dataset
* But plugin call is: plugin call x y, ...
* So in plugin varlist: x=1, y=2

cencode x, generate(y)
list
