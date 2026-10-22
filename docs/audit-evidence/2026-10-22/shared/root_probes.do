clear all
set more off
set linesize 220
adopath ++ /Users/mike/Documents/GitHub/stata-ctools/package/ado
adopath ++ /tmp/ctools-audit-oct22/build
capture log close _all
log using /tmp/ctools-audit-oct22/root_probes.log, text replace
ctools, version
clear
input double x
1e16
1
2
3
end
gen key=_n
crangestat (mean) cm=x (sum) cs=x, interval(key 0 0)
capture noisily rangestat (mean) rm=x (sum) rs=x, interval(key 0 0)
list, noobs
clear
input double (t ps y)
0 .1 2
0 .2 3
2 .15 7
2 .25 8
end
capture noisily cpsmatch t, pscore(ps) outcome(y)
di "NONBINARY_RC=" _rc
return list
list, noobs
capture noisily psmatch2 t, pscore(ps) outcome(y)
di "NONBINARY_NATIVE_RC=" _rc
return list
clear
input double (t ps y)
0 .1 2
0 .2 3
1 .15 7
1 .25 8
end
foreach opts in "neighbor(0)" "neighbor(-1)" "caliper(-.1)" "kernel bwidth(0)" "kernel kerneltype(nonsense)" {
 capture noisily cpsmatch t, pscore(ps) outcome(y) `opts'
 di "OPTION=`opts' RC=" _rc
 return list
}
clear
input double (t ps y)
0 .1 2
0 .2 3
1 .15 7
1 .25 8
end
clonevar _pscore=ps
capture noisily cpsmatch t, pscore(_pscore) outcome(y)
di "RESERVED_INPUT_RC=" _rc
capture confirm variable _pscore
di "RESERVED_INPUT_EXISTS_RC=" _rc
capture drop _pscore
clonevar _pscore=ps
capture noisily psmatch2 t, pscore(_pscore) outcome(y)
di "RESERVED_NATIVE_RC=" _rc
log close
exit, clear
