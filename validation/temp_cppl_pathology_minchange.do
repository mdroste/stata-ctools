do validation/validate_setup.do

capture program drop cmp_case2
program define cmp_case2
    syntax, cmd(string)
    matrix pV = e(V)
    matrix pB = e(b)
    capture noisily cpplmhdfe `cmd'
    matrix cV = e(V)
    matrix cB = e(b)
    local K = min(colsof(pB), colsof(cB))
    local min_sf_v = 15
    forvalues i=1/`K' {
        forvalues j=1/`K' {
            sigfigs pV[`i',`j'] cV[`i',`j']
            if r(sigfigs) < `min_sf_v' local min_sf_v = r(sigfigs)
        }
    }
    di as txt "cmd=" as txt "`cmd'" as txt " min_sf_v=" as res %6.3f `min_sf_v'
end

* zero-inflated
clear
set seed 222
set obs 2000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(-2 + 0.3 * x)
gen y = rpoisson(mu)
qui ppmlhdfe y x, absorb(fe1) vce(robust)
cmp_case2, cmd("y x, absorb(fe1) vce(robust)")
cmp_case2, cmd("y x, absorb(fe1) vce(robust) tolerance(1e-12)")
cmp_case2, cmd("y x, absorb(fe1) vce(robust) irlstolerance(1e-14)")
cmp_case2, cmd("y x, absorb(fe1) vce(robust) tolerance(1e-12) irlstolerance(1e-14)")

* mixed scale
clear
set seed 223
set obs 2000
gen fe1 = mod(_n-1, 50) + 1
gen x_tiny = rnormal() / 1e6
gen x_huge = rnormal() * 1e6
gen mu = exp(1e6 * x_tiny + 1e-6 * x_huge)
gen y = rpoisson(mu)
qui ppmlhdfe y x_tiny x_huge, absorb(fe1) vce(robust)
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust)")
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust) tolerance(1e-12)")
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust) irlstolerance(1e-14)")
cmp_case2, cmd("y x_tiny x_huge, absorb(fe1) vce(robust) tolerance(1e-12) irlstolerance(1e-14)")
