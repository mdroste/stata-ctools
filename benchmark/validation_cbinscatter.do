*! Validation tests for cbinscatter - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cbinscatter.log", replace text
validation_header "cbinscatter"

sysuse auto, clear

* === Basic binscatter ===
capture cbinscatter price mpg, nograph
local passed = (_rc == 0 & e(N) > 0)
report_test "Basic cbinscatter" `passed'

* === nquantiles options ===
foreach nq in 5 10 15 20 {
    capture cbinscatter price mpg, nquantiles(`nq') nograph
    local passed = (_rc == 0 & e(nquantiles) == `nq')
    report_test "nquantiles(`nq')" `passed'
}

* === Controls and absorb ===
run_test "controls()" cbinscatter price mpg, controls(weight length) nograph
run_test "absorb()" cbinscatter price mpg, absorb(foreign) nograph
run_test "controls+absorb" cbinscatter price mpg, controls(weight) absorb(foreign) nograph

* === By option ===
capture cbinscatter price mpg, by(foreign) nograph
local passed = (_rc == 0 & e(num_groups) == 2)
report_test "by(foreign)" `passed'

capture gen rep_group = rep78
capture replace rep_group = 3 if mi(rep_group)
run_test "by() multi-groups" cbinscatter price mpg, by(rep_group) nograph

* === Linetype options ===
foreach lt in none linear qfit cubic {
    run_test "linetype(`lt')" cbinscatter price mpg, linetype(`lt') nograph
}

* === Discrete option ===
run_test "discrete" cbinscatter price rep78, discrete nograph

* === Savedata option ===
tempfile bindata
capture cbinscatter price mpg, savedata(`bindata') nograph
local p = (_rc == 0)
capture confirm file `bindata'.dta
local passed = (`p' & _rc == 0)
report_test "savedata()" `passed'

* === Method options ===
run_test "method(classic)" cbinscatter price mpg, method(classic) controls(weight) nograph
run_test "method(binsreg)" cbinscatter price mpg, method(binsreg) controls(weight) nograph

* === Output options ===
run_test "nograph" cbinscatter price mpg, nograph
run_test "verbose" cbinscatter price mpg, verbose nograph
run_test "timeit" cbinscatter price mpg, timeit nograph
run_test "threads(2)" cbinscatter price mpg, threads(2) nograph

* === Weights ===
run_test "aweight" cbinscatter price mpg [aw=weight], nograph

capture gen int_weight = round(weight/100)
run_test "fweight" cbinscatter price mpg [fw=int_weight], nograph

* === Conditionals ===
run_test "if condition" cbinscatter price mpg if foreign == 1, nograph
run_test "in range" cbinscatter price mpg in 1/50, nograph

* === Title options ===
run_test "title()" cbinscatter price mpg, title("Test Title") nograph
run_test "ytitle()" cbinscatter price mpg, ytitle("Y Axis") nograph
run_test "xtitle()" cbinscatter price mpg, xtitle("X Axis") nograph

* === Colors option ===
run_test "colors()" cbinscatter price mpg, by(foreign) colors(red blue) nograph

* === Stored results ===
capture cbinscatter price mpg, nograph
local passed = (e(N) > 0 & e(nquantiles) > 0)
report_test "e() results" `passed'

* === Combined options ===
run_test "Combined options" cbinscatter price mpg, nquantiles(15) controls(weight) absorb(foreign) linetype(qfit) method(classic) verbose timeit threads(2) nograph

* === Error case ===
capture cbinscatter price mpg, nquantiles(1) nograph
local passed = (_rc == 198)
report_test "nquantiles(1) error" `passed'

* === Large dataset ===
clear
set obs 10000
gen y = rnormal()
gen x = rnormal() + 0.5 * y
run_test "Large dataset (10k)" cbinscatter y x, nograph

validation_summary
log close
