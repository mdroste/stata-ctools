* Manually compute canonical correlations for Test 4
clear all
set more off
adopath ++ "build"

sysuse auto, clear

* Absorb foreign FE manually (demean by foreign)
foreach v in price mpg weight length turn displacement {
    bys foreign: egen mean_`v' = mean(`v')
    gen `v'_dm = `v' - mean_`v'
}

* Y = endogenous (mpg, weight) after demeaning
* Z = excluded instruments (length, turn, displacement) after demeaning

* Compute Y'Y
mat accum YtY = mpg_dm weight_dm, noconstant
mat list YtY

* Compute Z'Z
mat accum ZtZ = length_dm turn_dm displacement_dm, noconstant
mat list ZtZ

* Compute Y'Z
mat accum YtZ_full = mpg_dm weight_dm length_dm turn_dm displacement_dm, noconstant
mat YtZ = YtZ_full[1..2, 3..5]
mat list YtZ

* Compute eigenvalues: inv(Y'Y) * Y'Z * inv(Z'Z) * Z'Y
mat ZtZ_inv = inv(ZtZ)
mat YtZ_ZtZ_inv = YtZ * ZtZ_inv
mat temp = YtZ_ZtZ_inv * YtZ'
mat C = inv(YtY) * temp
mat list C

* Eigenvalues of C
mat symeigen eigvec eigval = C
mat list eigval
mat list eigvec

di "Sum of eigenvalues: " eigval[1,1] + eigval[1,2]
di "Expected underid_stat ≈ N * sum(eigenvalues) = " _N * (eigval[1,1] + eigval[1,2])

* Also run ivreghdfe to compare
ivreghdfe price (mpg weight = length turn displacement), absorb(foreign)
di "ivreghdfe idstat = " e(idstat)
di "ivreghdfe widstat = " e(widstat)
