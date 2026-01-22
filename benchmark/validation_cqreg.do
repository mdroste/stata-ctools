*! Validation tests for cqreg - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cqreg.log", replace text
validation_header "cqreg"

sysuse auto, clear

* === Basic quantile regression ===
capture cqreg price mpg weight
local passed = (_rc == 0 & e(N) > 0)
report_test "Basic cqreg (q=0.5)" `passed'

* === Quantile options ===
foreach q in 0.10 0.25 0.50 0.75 0.90 {
    run_test "quantile(`q')" cqreg price mpg weight, quantile(`q')
}

* === Extreme quantiles ===
run_test "quantile(0.05)" cqreg price mpg weight, quantile(0.05)
run_test "quantile(0.95)" cqreg price mpg weight, quantile(0.95)

* === VCE options ===
run_test "vce(iid)" cqreg price mpg weight, vce(iid)
run_test "vce(robust)" cqreg price mpg weight, vce(robust)
run_test "vce(cluster)" cqreg price mpg weight, vce(cluster foreign)

* === Absorb ===
run_test "absorb()" cqreg price mpg weight, absorb(foreign)
run_test "absorb+quantile" cqreg price mpg weight, absorb(foreign) quantile(0.25)

* === Density estimation methods ===
foreach dm in fitted residual kernel {
    run_test "denmethod(`dm')" cqreg price mpg weight, denmethod(`dm')
}

* === Bandwidth methods ===
foreach bw in hsheather bofinger chamberlain {
    run_test "bwmethod(`bw')" cqreg price mpg weight, bwmethod(`bw')
}

* === Algorithm options ===
run_test "tolerance(1e-6)" cqreg price mpg weight, tolerance(1e-6)
run_test "maxiter(1000)" cqreg price mpg weight, maxiter(1000)
run_test "nopreprocess(1)" cqreg price mpg weight, nopreprocess(1)

* === Output options ===
run_test "verbose" cqreg price mpg weight, verbose
run_test "threads(2)" cqreg price mpg weight, threads(2)

* === Conditionals ===
run_test "if condition" cqreg price mpg weight if foreign == 1
run_test "in range" cqreg price mpg weight in 1/50

* === Multiple regressors ===
run_test "Multiple regressors" cqreg price mpg weight length headroom trunk

* === Stored results ===
capture cqreg price mpg weight
local passed = (e(N) > 0 & e(q) == 0.5)
report_test "e() results" `passed'

* === Combined options ===
run_test "Combined options" cqreg price mpg weight, quantile(0.75) vce(robust) denmethod(kernel) bwmethod(hsheather) verbose tolerance(1e-7) threads(2)

* === Coefficient checks ===
capture cqreg price mpg weight
local passed = (_b[weight] > 0)
report_test "Coef sign (weight>0)" `passed'

* === Different quantiles give different coefficients ===
capture cqreg price mpg weight, quantile(0.25)
local b_q25 = _b[weight]
capture cqreg price mpg weight, quantile(0.75)
local b_q75 = _b[weight]
local passed = (abs(`b_q25' - `b_q75') > 0.001)
report_test "Quantiles differ" `passed'

* === Large dataset ===
clear
set obs 5000
gen price = 1000 + 100 * rnormal()
gen mpg = 20 + 5 * rnormal()
gen weight = 3000 + 500 * rnormal()
run_test "Large dataset (5k)" cqreg price mpg weight

validation_summary
log close
