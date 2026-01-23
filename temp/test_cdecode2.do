capture do "validation/validate_setup.do"

sysuse auto, clear
cdecode foreign, generate(foreign_str)
count if !missing(foreign_str)
di "Non-missing: " r(N)
list foreign foreign_str in 1/10
