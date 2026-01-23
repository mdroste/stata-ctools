* Check state2 in census
sysuse census, clear
describe state2
list state2 in 1/5

* Check lifeexp country
webuse lifeexp, clear
describe country
list country in 1/5

* Try cencode on state2
sysuse census, clear
capture cencode state2, generate(test)
di "state2 cencode rc = " _rc

* Check what's different
webuse lifeexp, clear
capture cencode country, generate(test)
di "country cencode rc = " _rc

* Simple if test
sysuse auto, clear
cencode make if foreign, generate(make_code)
tab make_code
