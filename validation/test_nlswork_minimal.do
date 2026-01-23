* Minimal test for nlswork crash
clear all
webuse nlswork, clear
keep in 1/10000

* Check the data
describe ln_wage age tenure idcode
summarize idcode, detail
count if missing(ln_wage) | missing(age) | missing(tenure) | missing(idcode)

* Try creghdfe with verbose
creghdfe ln_wage age tenure, absorb(idcode) verbose
