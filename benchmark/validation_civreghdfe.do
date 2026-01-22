*! Validation tests for civreghdfe - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_civreghdfe.log", replace text
validation_header "civreghdfe"

webuse nlswork, clear
keep in 1/5000

* === Basic 2SLS ===
capture civreghdfe ln_wage (tenure = union) age, absorb(idcode)
local passed = (_rc == 0 & e(N) > 0)
report_test "Basic 2SLS" `passed'

* === Endogenous and instruments ===
run_test "Multiple endogenous" civreghdfe ln_wage (tenure hours = union south) age, absorb(idcode)
run_test "Multiple instruments" civreghdfe ln_wage (tenure = union south wks_ue) age, absorb(idcode)

* === Fixed effects ===
run_test "Multiple FE" civreghdfe ln_wage (tenure = union) age, absorb(idcode year)

* === VCE options ===
run_test "vce(robust)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) vce(robust)
run_test "vce(cluster)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) vce(cluster idcode)
run_test "vce(unadjusted)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) vce(unadjusted)
run_test "Two-way cluster" civreghdfe ln_wage (tenure = union) age, absorb(idcode) vce(cluster idcode year)

* === Estimation methods ===
run_test "liml" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) liml
run_test "fuller(1)" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) fuller(1)
run_test "kclass(0.9)" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) kclass(0.9)
run_test "gmm2s" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) gmm2s vce(robust)

* === First stage options ===
run_test "first" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) first
run_test "ffirst" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) ffirst

* === Algorithm options ===
run_test "tolerance(1e-6)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) tolerance(1e-6)
run_test "maxiter(1000)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) maxiter(1000)

* === Output options ===
run_test "verbose" civreghdfe ln_wage (tenure = union) age, absorb(idcode) verbose
run_test "timeit" civreghdfe ln_wage (tenure = union) age, absorb(idcode) timeit
run_test "threads(2)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) threads(2)
run_test "noheader" civreghdfe ln_wage (tenure = union) age, absorb(idcode) noheader
run_test "nofooter" civreghdfe ln_wage (tenure = union) age, absorb(idcode) nofooter
run_test "nooutput" civreghdfe ln_wage (tenure = union) age, absorb(idcode) nooutput
run_test "level(90)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) level(90)
run_test "noid" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode) noid
run_test "small" civreghdfe ln_wage (tenure = union) age, absorb(idcode) small
run_test "title()" civreghdfe ln_wage (tenure = union) age, absorb(idcode) title("Custom Title")
run_test "depname()" civreghdfe ln_wage (tenure = union) age, absorb(idcode) depname("Log Wage")

* === Diagnostic options ===
run_test "orthog()" civreghdfe ln_wage (tenure = union wks_ue south) age, absorb(idcode) orthog(south)
run_test "endogtest()" civreghdfe ln_wage (tenure hours = union wks_ue south) age, absorb(idcode) endogtest(hours)
run_test "redundant()" civreghdfe ln_wage (tenure = union wks_ue south) age, absorb(idcode) redundant(south)
run_test "partial()" civreghdfe ln_wage (tenure = union wks_ue) age ttl_exp, absorb(idcode) partial(ttl_exp)

* === HAC standard errors ===
run_test "HAC (bartlett)" civreghdfe ln_wage (tenure = union) age, absorb(idcode) kernel(bartlett) bw(3)

* === Weights ===
run_test "aweight" civreghdfe ln_wage (tenure = union) age [aw=hours], absorb(idcode)

* === Stored results ===
capture civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode)
local passed = (e(N) > 0 & e(K_endog) == 1 & e(sargan) >= 0)
report_test "e() results" `passed'

* === Combined options ===
run_test "Combined options" civreghdfe ln_wage (tenure = union wks_ue) age, absorb(idcode year) vce(cluster idcode) liml verbose timeit threads(2) ffirst

validation_summary
log close
