adopath ++ "./build"
sysuse auto, clear
di "Total obs: " _N
count if foreign == 1
di "Foreign cars: " r(N)
count if foreign == 1 in 1/50
di "Foreign cars in 1/50: " r(N)

* Test native encode
encode make if foreign == 1 in 1/50, generate(make_native)
count if !missing(make_native)
di "encode non-missing: " r(N)

drop make_native

* Test cencode
cencode make if foreign == 1 in 1/50, generate(make_code) verbose
count if !missing(make_code)
di "cencode non-missing: " r(N)
list make foreign make_code in 1/20
