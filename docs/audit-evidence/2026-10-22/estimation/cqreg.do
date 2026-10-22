clear all
set more off
set linesize 255
capture log close _all
log using "/tmp/ctools-regression-audit/cqreg.log", text replace
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
clear
set obs 2000
set seed 9854
gen double x=runiform()*5
gen double y=2+ .5*x +(1+4*x)*rnormal()
gen int cl=ceil(_n/20)
quietly qreg y x, vce(excluster cl, fitted)
matrix Vr=e(V)
matrix Br=e(b)
local rdf=e(df_r)
local rr=e(rank)
quietly cqreg y x, vce(cluster cl)
matrix Vc=e(V)
matrix Bc=e(b)
di "CLUSTER FITTED: qreg V=" Vr[1,1] " cqreg V=" Vc[1,1] " ratio=" Vc[1,1]/Vr[1,1] " b ref=" Br[1,1] " b c=" Bc[1,1]
di "CLUSTER DF: ref=`rdf' cqreg=" e(df_r) " rank ref=`rr' c=" e(rank)
quietly qreg y x, vce(excluster cl, kernel(gaussian))
matrix Vr=e(V)
quietly cqreg y x, vce(cluster cl) denmethod(kernel(gaussian))
matrix Vc=e(V)
di "CLUSTER KERNEL: qreg V=" Vr[1,1] " cqreg V=" Vc[1,1] " ratio=" Vc[1,1]/Vr[1,1]
gen double x2=x
quietly qreg y x x2
local rr=e(rank)
local rdf=e(df_r)
quietly cqreg y x x2
di "COLLINEAR RANK: qreg rank=`rr' cqreg=" e(rank) " refdf=`rdf' cdf=" e(df_r)
gen int category=mod(_n,3)
quietly qreg y x i.category
local rr=e(rank)
quietly cqreg y x i.category
di "FACTOR RANK: qreg rank=`rr' cqreg=" e(rank)
replace cl=1
capture noisily cqreg y x, vce(cluster cl)
di "ONE CLUSTER: rc=" _rc
replace y=y*1e-14
capture noisily cqreg y x
di "TINY OUTCOME: rc=" _rc
di "FINISHED"
log close
exit, clear
