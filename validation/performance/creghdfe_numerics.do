* Before/after benchmark for CG guards and response scaling.
* Select a plugin/ado directory with CTOOLS_TEST_PLUGIN_DIR in the caller.
* Run each build in a fresh Stata process; compare alternating A/B/B/A runs.
args output label smoke
if "`output'" == "" local output "temp/creghdfe-numerics"
if "`label'" == "" local label candidate
set seed 251125
set more off
tempname timings
postfile `timings' str32 path long N byte threads int trial double seconds using "`output'/timing-`label'.dta", replace
local paths no_fe one_fe two_fe three_fe weighted unscaled quad cimmino kaczmarz saved_fe near_fe slopes iv bins partial poisson poisson_slopes bootstrap
local trials = cond("`smoke'"=="smoke",1,3)
foreach n in 2000 200000 {
    clear
    set obs `n'
    gen long a = 1+mod(_n-1,200)
    gen long b = runiformint(1,80)
    gen long c = runiformint(1,30)
    gen double z = rnormal()
    gen double x = .5*z+rnormal()
    gen double y = 2*x-.7*z+a/200-b/80+rnormal()
    gen double w = 1+runiform()
    gen double near = a+.001*x
    gen long count = rpoisson(exp(.1*x+.1*z+a/200-b/80))
    local repeats = cond(`n'==2000,20,3)
    if "`smoke'"=="smoke" local repeats = 1
    foreach nt in 1 4 {
        local cmd_no_fe "creghdfe y x z, noabsorb threads(`nt')"
        local cmd_one_fe "creghdfe y x z, absorb(a) threads(`nt')"
        local cmd_two_fe "creghdfe y x z, absorb(a b) threads(`nt')"
        local cmd_three_fe "creghdfe y x z, absorb(a b c) threads(`nt')"
        local cmd_weighted "creghdfe y x z [aw=w], absorb(a b) threads(`nt')"
        local cmd_unscaled "creghdfe y x z, absorb(a b) nostandardize threads(`nt')"
        local cmd_quad "creghdfe y x z [aw=w], absorb(a b) quad threads(`nt')"
        local cmd_cimmino "creghdfe y x z, absorb(a b) transform(cimmino) threads(`nt')"
        local cmd_kaczmarz "creghdfe y x z, absorb(a b) transform(kaczmarz) threads(`nt')"
        local cmd_saved_fe "creghdfe y x z, absorb(a b, savefe) threads(`nt')"
        local cmd_near_fe "creghdfe y near z, absorb(a b) tolerance(1e-12) threads(`nt')"
        local cmd_slopes "creghdfe y x, absorb(a##c.z b) threads(`nt')"
        local cmd_iv "civreghdfe y (x=z), absorb(a b) threads(`nt')"
        local cmd_bins "cbinscatter y x, controls(z) absorb(a b) nograph threads(`nt')"
        local cmd_partial "chdfe y x z, absorb(a b) generate(res_) threads(`nt')"
        local cmd_poisson "cppmlhdfe count x z, absorb(a b) threads(`nt')"
        local cmd_poisson_slopes "cppmlhdfe count x, absorb(a##c.z) threads(`nt')"
        local cmd_bootstrap "creghdfe y x z, absorb(a b) vce(bootstrap, reps(10) seed(251125)) threads(`nt')"
        foreach path of local paths {
            if "`path'"=="bootstrap" & `n'>2000 continue
            capture drop res_*
            quietly `cmd_`path''
            forvalues trial=1/`trials' {
                timer clear 1
                timer on 1
                forvalues repeat=1/`repeats' {
                    if "`path'"=="partial" quietly drop res_*
                    quietly `cmd_`path''
                }
                timer off 1
                quietly timer list 1
                post `timings' ("`path'") (`n') (`nt') (`trial') (r(t1)/`repeats')
            }
            noisily display "TIMED `label' `path' N=`n' threads=`nt'"
        }
    }
}
* The individual-effects solver shares response scaling, but uses real LSMR.
clear
set seed 4242
set obs 600
gen long patent = _n
gen byte year = 1+mod(_n,12)
gen double x = rnormal()
gen double y = 2*x+year/10+rnormal()
expand 3
gen int inventor = runiformint(1,300)
bysort patent inventor: keep if _n==1
foreach nt in 1 4 {
    quietly creghdfe y x, absorb(year inventor) group(patent) individual(inventor) threads(`nt')
    local repeats = cond("`smoke'"=="smoke",1,20)
    forvalues trial=1/`trials' {
        timer clear 1
        timer on 1
        forvalues repeat=1/`repeats' {
            quietly creghdfe y x, absorb(year inventor) group(patent) individual(inventor) threads(`nt')
        }
        timer off 1
        quietly timer list 1
        post `timings' ("individual") (600) (`nt') (`trial') (r(t1)/`repeats')
    }
}
postclose `timings'
use "`output'/timing-`label'.dta", clear
export delimited using "`output'/timing-`label'.csv", replace
display "NUMERICS_BENCH_COMPLETE `label'"
