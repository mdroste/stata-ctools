* Debug delimiter edge cases section

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* Load helpers
quietly do "validation/validate_setup.do"

di "=== Testing pipe delimiter ==="

* Create test file
file open fh using "temp/pipe_delim.csv", write replace
file write fh "id|name|value" _n
file write fh "1|Alpha|100" _n
file write fh "2|Beta|200" _n
file write fh "3|Gamma|300" _n
file close fh

benchmark_import using "temp/pipe_delim.csv", testname("pipe delimiter") importopts(`"delimiters("|")"')

di "=== Testing space delimiter ==="

* Create test file
file open fh using "temp/space_delim.csv", write replace
file write fh "id name value" _n
file write fh "1 Alpha 100" _n
file write fh "2 Beta 200" _n
file write fh "3 Gamma 300" _n
file close fh

benchmark_import using "temp/space_delim.csv", testname("space delimiter") importopts(`"delimiters(" ")"')

di "=== Done ==="
