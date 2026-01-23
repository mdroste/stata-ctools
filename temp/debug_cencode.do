* Debug cencode
clear all

* Add build directory to adopath
adopath + "./build"

sysuse auto, clear
list make in 1/5
describe make

di "About to run cencode..."
capture drop make_code
cencode make, generate(make_code) verbose

di "After cencode:"
list make make_code in 1/5

sum make_code
