clear all
set more off
set linesize 255
capture log close _all
log using "/tmp/ctools-regression-audit/probe.log", text replace
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
which civreghdfe
which cpplmhdfe
capture which ppmlhdfe
capture which ivreghdfe
_ctools_load
return list
clear
set obs 60
set seed 8754
gen double z = rnormal()
gen double u = rnormal()
gen double x = z + u
gen double y = 2*x + .4*u + rnormal()
gen int fe = ceil(_n/10)
foreach model in "" "small" {
 quietly ivregress 2sls y (x=z), `model'
 matrix vr=e(V)
 scalar refrmse=e(rmse)
 local refdf=e(df_r)
 quietly civreghdfe y (x=z), `model'
 matrix vc=e(V)
 di "IV MODEL `model': V ref=" %16.12g vr[1,1] " c=" %16.12g vc[1,1] " df ref=`refdf' c=" e(df_r) " rmse ref=" refrmse " c=" e(rmse)
}
foreach model in "dofminus(8)" "small dofminus(8)" "small sdofminus(8)" "vce(robust)" "vce(robust) small" "vce(cluster fe)" "vce(cluster fe) small" {
 quietly civreghdfe y (x=z), `model'
 matrix vc=e(V)
 di "CIV EXTRA `model': V=" %16.12g vc[1,1]
}
clear
set obs 100
set seed 8754
gen double x=rnormal()
gen byte fe=mod(_n,5)
gen double y=rpoisson(exp(.1 + .3*x))
quietly cpplmhdfe y x, noabsorb
predict double mu_c, mu
predict double dev_c, deviance
gen double dev_exp=sign(y-mu_c)*sqrt(2*cond(y>0,mu_c-y+y*ln(y/mu_c),mu_c))
gen double difference=abs(dev_exp-dev_c)
quietly summarize difference
local devdiff=r(max)
quietly count if dev_c<0
local neg_c=r(N)
quietly count if dev_exp<0
di "PPML DEVIANCE: maxdiff=`devdiff' negative_c=`neg_c' negative_expected=" r(N)
list y mu_c dev_c dev_exp in 1/5
capture predict float xb_f, xb
local floatrc=_rc
local xbt : type xb_f
di "PPML FLOAT: rc=`floatrc' type=`xbt'"
quietly cpplmhdfe y x, absorb(fe) d(d)
predict double mud, mu
capture predict float xb_f2, xb
local floatrc=_rc
local xbt : type xb_f2
di "PPML ABS FLOAT: rc=`floatrc' type=`xbt'"
quietly creghdfe y x, absorb(fe) resid
capture predict float xbr, xb
local floatrc=_rc
capture local xbt : type xbr
di "REG FLOAT: rc=`floatrc' type=`xbt'"
di "FINISHED"
log close
exit, clear
