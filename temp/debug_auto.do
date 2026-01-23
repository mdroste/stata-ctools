adopath ++ "./build"
sysuse auto, clear
di "Number of variables: " _N
describe
capture drop make_code
cencode make, generate(make_code) verbose
di "After cencode:"
describe make_code
sum make_code
list make make_code in 1/10
