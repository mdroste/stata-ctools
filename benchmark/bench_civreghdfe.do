* Benchmark script for civreghdfe optimization
* Tests with 10 million observations

clear all
set more off
set seed 12345

* Configuration
local N = 10000000  // 10 million observations
local n_fe = 10000  // Number of fixed effect groups

di as text "=============================================="
di as text "civreghdfe Benchmark"
di as text "Observations: `N'"
di as text "FE groups: `n_fe'"
di as text "=============================================="

* Generate random data
di as text _n "Generating data..."
timer clear 1
timer on 1

set obs `N'

* Exogenous variables
gen x1 = rnormal()
gen x2 = rnormal()

* Instruments (for endogenous variable)
gen z1 = rnormal()
gen z2 = rnormal()

* Endogenous variable (correlated with error through omitted var)
gen omitted = rnormal()
gen endog = 0.5*z1 + 0.5*z2 + 0.3*omitted + rnormal()*0.5

* Fixed effects
gen fe1 = ceil(runiform() * `n_fe')

* Dependent variable
gen y = 1 + 2*x1 + 3*x2 + 4*endog + fe1/`n_fe' + omitted + rnormal()

timer off 1
timer list 1
di as text "Data generation time: " r(t1) " seconds"

* Compress to minimize memory
compress

di as text _n "Data summary:"
summarize y x1 x2 endog z1 z2 fe1

* ============================================
* Benchmark civreghdfe
* ============================================

di as text _n "=============================================="
di as text "Running civreghdfe benchmark"
di as text "=============================================="

* Warm-up run (to load plugin)
di as text _n "Warm-up run (1M obs subset)..."
preserve
keep if _n <= 1000000
timer clear 2
timer on 2
quietly civreghdfe y x1 x2 (endog = z1 z2), absorb(fe1) verbose
timer off 2
timer list 2
di as text "Warm-up time: " r(t2) " seconds"
restore

* Full benchmark run
di as text _n "Full benchmark (10M observations)..."
timer clear 3
timer on 3
civreghdfe y x1 x2 (endog = z1 z2), absorb(fe1) verbose
timer off 3
timer list 3

di as text _n "=============================================="
di as text "BENCHMARK RESULT"
di as text "civreghdfe total time: " r(t3) " seconds"
di as text "=============================================="

* Store timing for comparison
scalar civreghdfe_time = r(t3)

* Save benchmark results
file open results using "benchmark/civreghdfe_bench_results.txt", write replace
file write results "civreghdfe_baseline," (civreghdfe_time) _n
file close results

di as text _n "Results saved to benchmark/civreghdfe_bench_results.txt"
di as text "Benchmark complete."
