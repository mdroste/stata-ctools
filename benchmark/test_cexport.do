* Test cexport functionality
* Add build directory to adopath to find ctools commands
adopath + "`c(pwd)'/build"
clear all
set obs 10000
gen id = _n
gen x = rnormal()
gen y = runiform()
gen str20 name = "test" + string(id)

* Test 1: Basic export
timer clear 1
timer on 1
cexport delimited id x y using "/tmp/test_cexport1.csv", replace verbose
timer off 1
timer list 1

* Verify file was created and has content
shell wc -l /tmp/test_cexport1.csv
shell head -3 /tmp/test_cexport1.csv

* Test 2: Export with strings
timer clear 2
timer on 2
cexport delimited id name x y using "/tmp/test_cexport2.csv", replace verbose
timer off 2
timer list 2

* Verify
shell wc -l /tmp/test_cexport2.csv
shell head -3 /tmp/test_cexport2.csv

* Test 3: Large dataset
clear
set obs 100000
gen id = _n
gen x = rnormal()
gen y = runiform()
gen z = rnormal()

timer clear 3
timer on 3
cexport delimited id x y z using "/tmp/test_cexport3.csv", replace verbose
timer off 3
timer list 3

shell wc -l /tmp/test_cexport3.csv

di "All tests completed successfully!"
