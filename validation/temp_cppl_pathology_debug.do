do validation/validate_setup.do

capture program drop cmp_case
program define cmp_case
    syntax, name(string)

    matrix pV = e(V)
    matrix pB = e(b)

    capture noisily cpplmhdfe `name'
    matrix cV = e(V)
    matrix cB = e(b)

    local pcols = colsof(pB)
    local ccols = colsof(cB)
    local K = min(`pcols', `ccols')

    local min_sf_b = 15
    forvalues j = 1/`K' {
        sigfigs pB[1,`j'] cB[1,`j']
        if r(sigfigs) < `min_sf_b' local min_sf_b = r(sigfigs)
    }

    local min_sf_v = 15
    forvalues i = 1/`K' {
        forvalues j = 1/`K' {
            sigfigs pV[`i',`j'] cV[`i',`j']
            if r(sigfigs) < `min_sf_v' local min_sf_v = r(sigfigs)
        }
    }

    di as txt "case=`name' min_sf_b=" as res %6.3f `min_sf_b' as txt " min_sf_v=" as res %6.3f `min_sf_v'
end

* Case 1: zero-inflated Y
clear
set seed 222
set obs 2000
gen fe1 = mod(_n-1, 20) + 1
gen x = rnormal()
gen mu = exp(-2 + 0.3 * x)
gen y = rpoisson(mu)
quietly ppmlhdfe y x, absorb(fe1) vce(robust)

cmp_case, name("y x, absorb(fe1) vce(robust)")
cmp_case, name("y x, absorb(fe1) vce(robust) irlstolerance(1e-14) tolerance(1e-12)")
cmp_case, name("y x, absorb(fe1) vce(robust) irlstolerance(1e-16) tolerance(1e-14)")

* Case 2: mixed covariate scales
clear
set seed 223
set obs 2000
gen fe1 = mod(_n-1, 50) + 1
gen x_tiny = rnormal() / 1e6
gen x_huge = rnormal() * 1e6
gen mu = exp(1e6 * x_tiny + 1e-6 * x_huge)
gen y = rpoisson(mu)
quietly ppmlhdfe y x_tiny x_huge, absorb(fe1) vce(robust)

cmp_case, name("y x_tiny x_huge, absorb(fe1) vce(robust)")
cmp_case, name("y x_tiny x_huge, absorb(fe1) vce(robust) irlstolerance(1e-14) tolerance(1e-12)")
cmp_case, name("y x_tiny x_huge, absorb(fe1) vce(robust) irlstolerance(1e-16) tolerance(1e-14)")
