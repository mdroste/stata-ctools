* Minimal test for data I/O refactoring
clear all
adopath ++ "build"
set obs 10

* Test 1: cencode with 1 variable (string)
gen str10 city = ""
replace city = "NYC" in 1
replace city = "LA" in 2
replace city = "NYC" in 3
replace city = "Chicago" in 4
replace city = "LA" in 5
replace city = "Boston" in 6
replace city = "NYC" in 7
replace city = "Chicago" in 8
replace city = "LA" in 9
replace city = "Boston" in 10

di "About to call cencode..."
cencode city, gen(city_num) verbose
di "cencode done: " _rc
list city city_num in 1/5

* Test 2: cdecode with 1 variable (numeric)
di "About to call cdecode..."
cdecode city_num, gen(city_str) verbose
di "cdecode done: " _rc
list city_num city_str in 1/5

* Test 3: cdestring with 1 variable
clear
set obs 5
gen str10 x = ""
replace x = "1.5" in 1
replace x = "2.3" in 2
replace x = "3.7" in 3
replace x = "4.1" in 4
replace x = "5.9" in 5

di "About to call cdestring with generate..."
cdestring x, gen(x_num) verbose
di "cdestring generate done: " _rc
list x x_num

di "ALL TESTS PASSED"
