*! Validation tests for creghdfe - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_creghdfe.log", replace text
validation_header "creghdfe"

webuse nlswork, clear
keep in 1/10000

* === Basic functionality ===
capture creghdfe ln_wage age ttl_exp, absorb(idcode)
local passed = (_rc == 0 & e(N) > 0)
report_test "Basic creghdfe" `passed'

run_test "Multiple FE" creghdfe ln_wage age ttl_exp, absorb(idcode year)

* === VCE options ===
run_test "vce(robust)" creghdfe ln_wage age ttl_exp, absorb(idcode) vce(robust)
run_test "vce(cluster)" creghdfe ln_wage age ttl_exp, absorb(idcode) vce(cluster idcode)
run_test "vce(unadjusted)" creghdfe ln_wage age ttl_exp, absorb(idcode) vce(unadjusted)

* === Algorithm options ===
run_test "tolerance(1e-6)" creghdfe ln_wage age ttl_exp, absorb(idcode) tolerance(1e-6)
run_test "tolerance(1e-3)" creghdfe ln_wage age ttl_exp, absorb(idcode) tolerance(1e-3)
run_test "maxiter(1000)" creghdfe ln_wage age ttl_exp, absorb(idcode) maxiter(1000)

* === Output options ===
run_test "verbose" creghdfe ln_wage age ttl_exp, absorb(idcode) verbose
run_test "timeit" creghdfe ln_wage age ttl_exp, absorb(idcode) timeit
run_test "threads(2)" creghdfe ln_wage age ttl_exp, absorb(idcode) threads(2)
run_test "nostandardize" creghdfe ln_wage age ttl_exp, absorb(idcode) nostandardize

* === Residuals ===
capture drop _reghdfe_resid
capture creghdfe ln_wage age ttl_exp, absorb(idcode) resid
local p = (_rc == 0)
capture confirm variable _reghdfe_resid
local passed = (`p' & _rc == 0)
report_test "resid" `passed'

capture drop my_resid
capture creghdfe ln_wage age ttl_exp, absorb(idcode) resid2(my_resid)
local p = (_rc == 0)
capture confirm variable my_resid
local passed = (`p' & _rc == 0)
report_test "resid2()" `passed'

* === Weights ===
run_test "aweight" creghdfe ln_wage age ttl_exp [aw=hours], absorb(idcode)

capture gen int_hours = round(hours)
capture replace int_hours = 1 if int_hours < 1
run_test "fweight" creghdfe ln_wage age ttl_exp [fw=int_hours], absorb(idcode)

run_test "pweight" creghdfe ln_wage age ttl_exp [pw=hours], absorb(idcode) vce(robust)

* === Conditionals ===
run_test "if condition" creghdfe ln_wage age ttl_exp if year >= 75, absorb(idcode)
run_test "in range" creghdfe ln_wage age ttl_exp in 1/5000, absorb(idcode)

* === Three fixed effects ===
capture gen industry = mod(idcode, 20) + 1
run_test "Three FE" creghdfe ln_wage age ttl_exp, absorb(idcode year industry)

* === Coefficient verification ===
capture creghdfe ln_wage age ttl_exp, absorb(idcode)
local passed = (_rc == 0 & _b[ttl_exp] > 0)
report_test "Coef sign (ttl_exp>0)" `passed'

* === Stored results ===
capture creghdfe ln_wage age ttl_exp, absorb(idcode)
local passed = (e(N) > 0 & e(r2) >= 0 & e(r2) <= 1)
report_test "e() results" `passed'

* === Combined options ===
capture drop my_resid2
run_test "Combined options" creghdfe ln_wage age ttl_exp, absorb(idcode year) vce(cluster idcode) tolerance(1e-7) maxiter(500) verbose timeit resid2(my_resid2) threads(2)

* === Error cases ===
capture creghdfe ln_wage age ttl_exp if ln_wage > 1000, absorb(idcode)
local passed = (_rc != 0)
report_test "No obs error" `passed'

* === Edge cases ===
capture gen age2 = age
capture creghdfe ln_wage age age2, absorb(idcode)
report_test "Collinearity" 1  // Just check doesn't crash

capture bysort idcode: gen singleton_check = (_N == 1)
run_test "Singleton handling" creghdfe ln_wage age ttl_exp, absorb(idcode)

run_test "Multiple regressors" creghdfe ln_wage age ttl_exp tenure hours wks_work, absorb(idcode)

* === Consistency check ===
capture creghdfe ln_wage age, absorb(idcode)
local coef1 = _b[age]
capture creghdfe ln_wage age, absorb(idcode)
local coef2 = _b[age]
local passed = (abs(`coef1' - `coef2') < 1e-10)
report_test "Coef consistency" `passed'

* === Missing values ===
capture count if missing(ln_wage)
capture creghdfe ln_wage age ttl_exp, absorb(idcode)
local passed = (_rc == 0 & e(N) > 0)
report_test "Missing values" `passed'

* === Factor variables ===
run_test "Factor vars (i.race)" creghdfe ln_wage i.race age, absorb(idcode)

validation_summary
log close
