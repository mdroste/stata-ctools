capture do "validation/validate_setup.do"

sysuse auto, clear
cencode make, generate(make_code)
count if !missing(make_code)
di "Non-missing: " r(N)
tab make_code in 1/5
