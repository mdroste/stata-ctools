capture do "validation/validate_setup.do"

sysuse auto, clear
cdecode rep78, generate(rep_str)
count if !missing(rep_str)
di "Non-missing: " r(N)
list rep78 rep_str in 1/10
