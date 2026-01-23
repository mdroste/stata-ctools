* Load setup first
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

* Try cencode on state2
sysuse census, clear
capture noisily cencode state2, generate(test)
di "state2 cencode rc = " _rc

* Check what's different about lifeexp
webuse lifeexp, clear
capture noisily cencode country, generate(test)
di "country cencode rc = " _rc

* Simple if test
sysuse auto, clear
capture noisily cencode make if foreign, generate(make_code)
di "make if foreign rc = " _rc
tab make_code

* Compare with encode
sysuse auto, clear
capture noisily encode make if foreign, generate(make_code2)
di "encode make if foreign rc = " _rc
tab make_code2
