*! Validation tests for csort - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_csort.log", replace text
validation_header "csort"

* === Basic functionality ===
sysuse auto, clear
capture csort price
local p1 = (_rc == 0)
local sorted = 1
forvalues i = 2/`=_N' {
    if price[`i'] < price[`i'-1] local sorted = 0
}
local passed = (`p1' & `sorted')
report_test "Basic numeric sort" `passed'

sysuse auto, clear
run_test "String sort" csort make

sysuse auto, clear
run_test "Multiple keys" csort foreign mpg

* === Algorithm options ===
foreach alg in ips4o lsd msd timsort sample counting merge {
    sysuse auto, clear
    run_test "algorithm(`alg')" csort price, algorithm(`alg')
}

* Algorithm abbreviations
sysuse auto, clear
run_test "algorithm(tim)" csort price, algorithm(tim)
run_test "algorithm(count)" csort price, algorithm(count)

* === Mode options ===
sysuse auto, clear
run_test "stream mode" csort price, stream

* === Other options ===
sysuse auto, clear
run_test "verbose" csort price, verbose

sysuse auto, clear
run_test "threads(2)" csort price, threads(2)

sysuse auto, clear
run_test "nosortedby" csort price, nosortedby

* === Conditionals ===
sysuse auto, clear
run_test "if condition" csort price if foreign == 1

sysuse auto, clear
run_test "in range" csort price in 1/50

* === Combined options ===
sysuse auto, clear
run_test "Combined options" csort price mpg, algorithm(lsd) threads(2) verbose

* === Stability test ===
clear
set obs 100
gen group = mod(_n - 1, 10) + 1
gen orig = _n
capture csort group
local p = (_rc == 0)
if `p' {
    bysort group: gen check = orig >= orig[_n-1] if _n > 1
    egen minc = min(check), by(group)
    count if minc == 0
    local p = (r(N) == 0)
}
report_test "Sort stability" `p'

* === Large dataset ===
clear
set obs 100000
gen x = runiform()
run_test "Large dataset (100k)" csort x

* === Missing values ===
clear
set obs 100
gen x = _n
replace x = . in 1/10
capture csort x
local passed = (_rc == 0 & mi(x[_N]))
report_test "Missing values" `passed'

* === Error cases ===
sysuse auto, clear
run_test_error "Invalid algorithm error" 198 csort price, algorithm(invalid)

* === Data integrity ===
sysuse auto, clear
local orig_n = _N
sum price, meanonly
local orig_sum = r(sum)
capture csort make
local p = (_rc == 0)
if `p' {
    sum price, meanonly
    local p = (`orig_n' == _N & abs(`orig_sum' - r(sum)) < 0.001)
}
report_test "Data integrity" `p'

validation_summary
log close
