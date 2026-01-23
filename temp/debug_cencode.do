capture do "validation/validate_setup.do"

di as text "DEBUG: Basic encoding without if"
sysuse auto, clear
cencode make, generate(make_code) verbose
count if !missing(make_code)
di "Non-missing count: " r(N)
tab make_code in 1/5
