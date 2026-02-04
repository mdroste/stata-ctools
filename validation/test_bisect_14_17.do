/*******************************************************************************
 * test_bisect_14_17.do - Test sections 14-17
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running sections 14-17 then clear all..."

/*******************************************************************************
 * SECTION 14: Coefficient comparison
 ******************************************************************************/
print_section "Coefficient Comparison"
sysuse auto, clear
ivreghdfe price (mpg = weight length), absorb(foreign)
matrix ivreghdfe_b = e(b)
civreghdfe price (mpg = weight length), absorb(foreign)
matrix civreghdfe_b = e(b)
matrix_min_sigfigs ivreghdfe_b civreghdfe_b
local min_sf = r(min_sigfigs)
if `min_sf' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f `min_sf'
    test_pass "coefficients match (sigfigs=`sf_fmt')"
}
else {
    local sf_fmt : display %4.1f `min_sf'
    test_fail "coefficients" "sigfigs=`sf_fmt' below threshold"
}

/*******************************************************************************
 * SECTION 15: Weights
 ******************************************************************************/
print_section "Weights"
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
print_section "Display Options"
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nofooter
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nooutput
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) title("Custom Title Test")
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) depname("outcome_var")
sysuse auto, clear
civreghdfe price (mpg = weight length), absorb(foreign) noid
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)

/*******************************************************************************
 * SECTION 17: Two-way clustering
 ******************************************************************************/
print_section "Two-Way Clustering"
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
local se_twoway = sqrt(e(V)[1,1])
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm)
local se_firm = sqrt(e(V)[1,1])
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster time)
local se_time = sqrt(e(V)[1,1])
civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) vce(cluster firm time)

noi di as text _n "Sections 14-17 complete. Now trying clear all..."

capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

}

di as text _n "Now running: clear all"
clear all

di as text "SUCCESS: No crash after sections 14-17!"
