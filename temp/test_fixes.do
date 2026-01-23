capture do "validation/validate_setup.do"

di as text "TEST 1: state2 variable (previously failing with rc=198)"
sysuse census, clear
capture noisily cencode state2, generate(test)
di "cencode state2 rc = " _rc
if _rc == 0 {
    di as result "[PASS] state2 encoded successfully"
    tab test in 1/10
}
else {
    di as error "[FAIL] state2 failed"
}
drop test

di ""
di as text "TEST 2: if condition (previously encoded all 74 cars)"
sysuse auto, clear
capture drop make_code
cencode make if foreign, generate(make_code)
count if !missing(make_code)
local n_encoded = r(N)
di "Observations encoded: " `n_encoded'
if `n_encoded' == 22 {
    di as result "[PASS] if condition works - only 22 foreign cars encoded"
}
else {
    di as error "[FAIL] if condition wrong - expected 22, got " `n_encoded'
}

di ""
di as text "TEST 3: Compare cencode vs encode with if condition"
sysuse auto, clear
encode make if foreign, generate(make_encode)
cencode make if foreign, generate(make_cencode)

* Both should have same count of non-missing
count if !missing(make_encode)
local n_encode = r(N)
count if !missing(make_cencode)
local n_cencode = r(N)
di "encode count: `n_encode', cencode count: `n_cencode'"

if `n_encode' == `n_cencode' {
    di as result "[PASS] Same number of observations encoded"
}
else {
    di as error "[FAIL] Different counts"
}

di ""
di as text "TEST 4: webuse lifeexp country"
webuse lifeexp, clear
capture noisily cencode country, generate(test)
if _rc == 0 {
    di as result "[PASS] lifeexp country encoded"
    count if !missing(test)
    di "Encoded " r(N) " observations"
}
else {
    di as error "[FAIL] lifeexp country failed with rc=" _rc
}

di ""
di as text "TEST 5: webuse grunfeld company"
webuse grunfeld, clear
capture noisily cencode company, generate(test)
if _rc == 0 {
    di as result "[PASS] grunfeld company encoded"
}
else {
    di as error "[FAIL] grunfeld company failed with rc=" _rc
}
