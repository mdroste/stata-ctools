* Debug pipe character in syntax parsing

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* Define a simple test program
capture program drop test_opts
program define test_opts
    syntax using/, testname(string) [IMPORTopts(string)]
    di "testname = `testname'"
    di "importopts = `importopts'"
end

* Test 1: Simple options
di "=== Test 1: Simple options ==="
test_opts using "test.csv", testname("test1") importopts("clear")

* Test 2: Pipe character in quotes - compound quotes
di _n "=== Test 2a: Pipe with compound quotes ==="
capture test_opts using "test.csv", testname("test2a") importopts(`"delimiters("|")"')
di "rc = " _rc

* Test 2b: Double quotes only
di _n "=== Test 2b: Double quotes only ==="
capture test_opts using "test.csv", testname("test2b") importopts("delimiters(|)")
di "rc = " _rc

* Test 2c: Try with escaped pipe
di _n "=== Test 2c: Escaped pipe ==="
capture test_opts using "test.csv", testname("test2c") importopts(`"delimiters("\|")"')
di "rc = " _rc

* Test 2d: Simple variable approach
di _n "=== Test 2d: Via local variable ==="
local myopts `"delimiters("|")"'
di "myopts = `myopts'"
capture test_opts using "test.csv", testname("test2d") importopts(`"`myopts'"')
di "rc = " _rc
