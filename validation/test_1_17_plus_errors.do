/*******************************************************************************
 * test_1_17_plus_errors.do - Sections 1-17 then error tests then clear all
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running sections 1-17 then error tests then clear all..."

* Plugin check
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)

* Sections 2-8: Basic tests
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("basic 1")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("basic 2")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) testname("basic 3")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) turn displacement, absorb(foreign) testname("basic 4")
sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement), absorb(foreign) testname("multi endog 1")
sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement headroom), absorb(foreign) testname("multi endog 2")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("vce robust")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("vce cluster")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("multi-instr robust")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) testname("two-way FE")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) vce(robust) testname("two-way robust")
sysuse census, clear
benchmark_ivreghdfe pop (medage = death), absorb(region) testname("census 1")
sysuse census, clear
benchmark_ivreghdfe pop (medage = death marriage), absorb(region) testname("census 2")
sysuse census, clear
benchmark_ivreghdfe pop (medage = death) divorce, absorb(region) testname("census 3")
webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) testname("nlswork 1")
webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age ttl_exp), absorb(idcode) testname("nlswork 2")
webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(robust) testname("nlswork robust")
webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(cluster idcode) testname("nlswork cluster")
clear
set seed 12345
set obs 20000
gen id = runiformint(1, 200)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) testname("20K dataset")

* Options
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) first
capture civreghdfe price (mpg = weight), absorb(foreign) small
capture civreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10)
capture civreghdfe price (mpg = weight), absorb(foreign) maxiter(1000)
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
capture civreghdfe price (mpg = weight), absorb(foreign) timeit

* if/in
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign) testname("if condition")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) in 1/50, absorb(foreign) testname("in condition")

* Coefficient comparison
sysuse auto, clear
ivreghdfe price (mpg = weight length), absorb(foreign)
matrix ivreghdfe_b = e(b)
civreghdfe price (mpg = weight length), absorb(foreign)
matrix civreghdfe_b = e(b)

* Weights
sysuse auto, clear
capture civreghdfe price (mpg = weight) [aw=weight], absorb(foreign)
sysuse auto, clear
gen int fw = ceil(mpg/5)
capture civreghdfe price (mpg = weight) [fw=fw], absorb(foreign)
sysuse auto, clear
capture civreghdfe price (mpg = weight) [pw=weight], absorb(foreign)

* Display options
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)

* Two-way clustering
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm)
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster time)
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)

noi di as text _n "Sections 1-17 complete. Now running error tests..."

/*******************************************************************************
 * ERROR TESTS
 ******************************************************************************/
* Variable doesn't exist
sysuse auto, clear
capture civreghdfe price (mpg = nonexistent_var), absorb(foreign)

* Missing absorb option
sysuse auto, clear
gen z = runiform()
capture civreghdfe price (mpg = z)

* No instruments (underidentified)
sysuse auto, clear
capture civreghdfe price mpg, absorb(foreign)

* No observations after if condition
sysuse auto, clear
gen z2 = runiform()
capture civreghdfe price (mpg = z2) if price > 100000, absorb(foreign)

* String variable as dependent
sysuse auto, clear
gen z3 = runiform()
capture civreghdfe make (mpg = z3), absorb(foreign)

noi di as text _n "Error tests complete. Now trying clear all..."

capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

}

di as text _n "Now running: clear all"
clear all

di as text "SUCCESS: No crash after sections 1-17 + error tests!"
