clear all
set more off
set linesize 255
capture log close _all
log using "/tmp/ctools-regression-audit/probe2.log", text replace
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
clear
set obs 60
set seed 8754
gen double z = rnormal()
gen double u = rnormal()
gen double x = z + u
gen double y = 10+ .05*x + .4*u + rnormal()
quietly ivregress 2sls y (x=z), small
scalar rf=e(F)
scalar rdf=e(df_m)
quietly civreghdfe y (x=z), small
scalar cf=e(F)
scalar cdf=e(df_m)
quietly test x
di "IV OMNIBUS: native F=" rf " df_m=" rdf " c F=" cf " df_m=" cdf " test x F=" r(F)
foreach opts in "noconstant small" "noconstant small robust" {
 quietly ivregress 2sls y (x=z), `opts'
 scalar rtss=e(tss)
 scalar rr2=e(r2)
 scalar rr2a=e(r2_a)
 quietly civreghdfe y (x=z), `opts'
 di "IV NOCONS `opts': native tss=" rtss " c tss=" e(tss) " native R2=" rr2 " c R2=" e(r2) " native R2adj=" rr2a " c=" e(r2_a)
}
* Missing time-series dependency after estimation
clear
set obs 100
gen t=_n
tsset t
set seed 884
gen double x=rnormal()
gen double y=exp(.2*x+rnormal())
quietly cpplmhdfe L.y x, noabsorb
local depvar `e(depvar)'
capture predict double rr, response
di "PPML TS RESPONSE: e(depvar)=`depvar' rc=" _rc
capture predict double mu, mu
di "PPML TS MU: rc=" _rc
* An allowed positive-DF model rejected by overly strong N cutoff
clear
input double(y x z)
1 1 1
4 3 2
2 2 3
5 4 4
end
capture noisily ivregress 2sls y (x=z), small
local ref_rc=_rc
capture noisily civreghdfe y (x=z), small
di "SMALL VALID MODEL: native_rc=`ref_rc' civ_rc=" _rc
log close
exit, clear
