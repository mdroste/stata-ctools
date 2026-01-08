* Compare e() contents between qreg and cqreg
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di "{hline 60}"
di "QREG OUTPUT"
di "{hline 60}"
qreg y x1 x2

di _newline "e() scalars:"
ereturn list

di _newline "{hline 60}"
di "CQREG OUTPUT"
di "{hline 60}"
cqreg y x1 x2

di _newline "e() scalars:"
ereturn list
