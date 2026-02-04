/*******************************************************************************
 * test_bisect_1_24.do - Run sections 1-24 (skip error tests), then clear all
 *
 * This is the full validation minus the "Intentional Error Tests" section.
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running sections 1-24 (no error tests) then clear all..."

* Check if ivreghdfe is installed
capture which ivreghdfe
if _rc != 0 {
    noi di as error "ivreghdfe not found"
    exit 0
}

* Plugin check
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)
if _rc != 0 {
    noi di as error "civreghdfe plugin load failed"
    exit 1
}

/*******************************************************************************
 * SECTION 2-8: Basic tests
 ******************************************************************************/
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

/*******************************************************************************
 * SECTION 9-12: Options
 ******************************************************************************/
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) first
capture civreghdfe price (mpg = weight), absorb(foreign) small
capture civreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10)
capture civreghdfe price (mpg = weight), absorb(foreign) maxiter(1000)
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
capture civreghdfe price (mpg = weight), absorb(foreign) timeit

/*******************************************************************************
 * SECTION 13: if/in
 ******************************************************************************/
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign) testname("if condition")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) in 1/50, absorb(foreign) testname("in condition")

/*******************************************************************************
 * SECTION 14: Coefficient comparison
 ******************************************************************************/
sysuse auto, clear
ivreghdfe price (mpg = weight length), absorb(foreign)
matrix ivreghdfe_b = e(b)
civreghdfe price (mpg = weight length), absorb(foreign)
matrix civreghdfe_b = e(b)
matrix_min_sigfigs ivreghdfe_b civreghdfe_b

/*******************************************************************************
 * SECTION 15: Weights
 ******************************************************************************/
sysuse auto, clear
capture civreghdfe price (mpg = weight) [aw=weight], absorb(foreign)
sysuse auto, clear
gen int fw = ceil(mpg/5)
capture civreghdfe price (mpg = weight) [fw=fw], absorb(foreign)
sysuse auto, clear
capture civreghdfe price (mpg = weight) [pw=weight], absorb(foreign)

/*******************************************************************************
 * SECTION 16: Display options
 ******************************************************************************/
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
capture civreghdfe price (mpg = weight), absorb(foreign) noheader
capture civreghdfe price (mpg = weight), absorb(foreign) nofooter
capture civreghdfe price (mpg = weight), absorb(foreign) nooutput
capture civreghdfe price (mpg = weight), absorb(foreign) title("Test")
civreghdfe price (mpg = weight), absorb(foreign) depname("outcome")
civreghdfe price (mpg = weight length), absorb(foreign) noid
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)

/*******************************************************************************
 * SECTION 17: Two-way clustering
 ******************************************************************************/
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

/*******************************************************************************
 * SECTION 18: Diagnostic tests
 ******************************************************************************/
scalar drop _all
matrix drop _all
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3)
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3)
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) noid
drop x_endog y
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
gen y = 2*x_endog1 + 1.5*x_endog2 + 1.0*x_exog + rnormal()
capture civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2)
civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2)
capture civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3)
civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3)

/*******************************************************************************
 * SECTION 19: partial()
 ******************************************************************************/
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog1 = runiform()
gen x_exog2 = rnormal()
gen y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm)
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
gen x_exog3 = rnormal()
replace y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + 0.3*x_exog3 + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2 x_exog3, absorb(firm) partial(x_exog2 x_exog3)

/*******************************************************************************
 * SECTION 20: ffirst
 ******************************************************************************/
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.0*x_exog + rnormal()
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
replace y = 2*x_endog + 1.5*x_endog2 + 1.0*x_exog + rnormal()
civreghdfe y (x_endog x_endog2 = z1 z2 z3) x_exog, absorb(firm) ffirst

/*******************************************************************************
 * SECTION 21: Save options
 ******************************************************************************/
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.0*x_exog + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst
capture estimates restore _civreghdfe_main
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst savefprefix(test_)
capture estimates restore test_main
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) saverf

/*******************************************************************************
 * SECTION 22: Factor variables
 ******************************************************************************/
sysuse auto, clear
capture civreghdfe price (mpg = weight) i.rep78, absorb(foreign)
capture ivreghdfe price (mpg = weight) i.rep78, absorb(foreign)
if _rc == 0 {
    civreghdfe price (mpg = weight) i.rep78, absorb(foreign)
}
capture civreghdfe price (mpg = i.rep78), absorb(foreign)
capture civreghdfe price (mpg = weight length) i.rep78 i.foreign, absorb(headroom)
capture civreghdfe price (mpg = weight) i.rep78##c.turn, absorb(foreign)
capture civreghdfe price (mpg = weight length) ib3.rep78, absorb(foreign)
capture civreghdfe price (mpg = weight) c.turn#i.foreign, absorb(rep78)
webuse nlswork, clear
keep in 1/5000
capture civreghdfe ln_wage (tenure = age) i.race#i.union, absorb(idcode)
sysuse auto, clear
capture ivreghdfe price (mpg = weight length) c.turn#i.foreign, absorb(rep78)
if _rc == 0 {
    capture civreghdfe price (mpg = weight length) c.turn#i.foreign, absorb(rep78)
}

/*******************************************************************************
 * SECTION 23: Time series
 ******************************************************************************/
clear
set seed 54321
set obs 500
gen id = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
tsset id time
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
capture civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = z1 z2) D.x_exog, absorb(id)
capture civreghdfe y (x_endog = F.z1 z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = L2.z1 z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = L(1/2).z1) x_exog, absorb(id)
capture ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id)
if _rc == 0 {
    civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id)
}
capture civreghdfe y (x_endog = L.z1 D.z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(robust)
capture civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(cluster id)

/*******************************************************************************
 * SECTION 24: IV Estimator Options
 ******************************************************************************/
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) liml testname("liml")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(1) testname("fuller(1)")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(4) testname("fuller(4)")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(0.5) testname("kclass(0.5)")
ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(1)
ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm)
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s testname("gmm2s")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust) testname("gmm2s + robust")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue testname("cue")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) coviv testname("coviv")
matrix b0 = (2, 1.5)
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) b0(b0)

/*******************************************************************************
 * SECTION 24b: GMM2S/CUE tests
 ******************************************************************************/
local strict_sigfigs = $DEFAULT_SIGFIGS

qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust)
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust)
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(cluster firm)
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue vce(cluster firm)

* Just-identified case
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen x_endog = 0.5*z1 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()
civreghdfe y (x_endog = z1) x_exog, absorb(firm)
civreghdfe y (x_endog = z1) x_exog, absorb(firm) gmm2s
civreghdfe y (x_endog = z1) x_exog, absorb(firm) cue

* Multiple endogenous with GMM2S/CUE
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen z4 = rnormal() + 0.2*z1
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z3 + 0.3*z4 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog1 + 1.5*x_endog2 + 1.0*x_exog + rnormal()
qui ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) gmm2s
qui civreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) gmm2s
qui ivreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) cue
qui civreghdfe y (x_endog1 x_endog2 = z1 z2 z3 z4) x_exog, absorb(firm) cue

/*******************************************************************************
 * SECTION 24c: New Options
 ******************************************************************************/
sysuse auto, clear
civreghdfe price (mpg = weight length), absorb(foreign) robust
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
civreghdfe price (mpg = weight length), absorb(foreign) cluster(foreign)
civreghdfe price (mpg = weight length), absorb(foreign)
civreghdfe price (mpg = weight length), absorb(foreign) dofminus(1)
capture civreghdfe price (mpg = weight length), absorb(foreign) eform(OR)
capture civreghdfe price (mpg = weight length), absorb(foreign) subtitle("Custom")

clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog1 = runiform()
gen x_exog2 = rnormal()
gen y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) fwl(x_exog2)

/*******************************************************************************
 * SECTION 25: HAC/Kernel
 ******************************************************************************/
clear
set seed 54321
set obs 500
gen id = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
tsset id time
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(bartlett) bw(2) testname("kernel bartlett")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3) testname("kernel parzen")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(qs) bw(2) testname("kernel qs")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2) testname("dkraay 2")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(4) testname("dkraay 4")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer testname("kiefer")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) bw(2) vce(robust) testname("bw robust")

/*******************************************************************************
 * SKIPPING: Intentional Error Tests
 ******************************************************************************/

noi di as text _n "Sections 1-25 (no error tests) complete. Now trying clear all..."

capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

}

di as text _n "Now running: clear all"
clear all

di as text "SUCCESS: No crash after sections 1-25 (no error tests)!"
