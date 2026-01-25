adopath ++ build
sysuse auto, clear

* Clustered VCE
di "=== ivreghdfe with cluster ==="
ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign)
local iv_widstat = e(widstat)
local iv_idstat = e(idstat)
local iv_kp_f = e(kp_f)

sysuse auto, clear
di ""
di "=== civreghdfe with cluster ==="
civreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign)
local civ_widstat = e(widstat)
local civ_idstat = e(idstat)
capture local civ_kp_f = e(kp_f)
if missing(`civ_kp_f') local civ_kp_f = .

di ""
di "=== Comparison ==="
di "ivreghdfe widstat = `iv_widstat'"
di "civreghdfe widstat = `civ_widstat'"
di "Ratio = " `civ_widstat' / `iv_widstat'
di ""
di "ivreghdfe idstat = `iv_idstat'"
di "civreghdfe idstat = `civ_idstat'"
di ""
di "ivreghdfe kp_f = `iv_kp_f'"
di "civreghdfe kp_f = `civ_kp_f'"
