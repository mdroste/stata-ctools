/**
 * test_cleanup.do
 *
 * Test the cleanup system for persistent state between plugin calls.
 * Verifies that switching between commands doesn't cause memory issues.
 */

clear all
set more off

* Test: Run multiple commands in sequence to ensure cleanup works
* This should not crash or leak memory

* 1. Create some test data
set obs 1000
gen id = _n
gen x = runiform()
gen y = rnormal()
gen z = rnormal()

* 2. Run csort
di "Testing csort..."
csort id
assert id[1] == 1
assert id[_N] == 1000
di "  csort passed"

* 3. Run csort again (same command, state preserved)
di "Testing csort again (same command)..."
gsort -id
csort id
assert id[1] == 1
di "  csort (repeated) passed"

* 4. Run cmerge load_using (this creates cache)
di "Testing cmerge load_using..."
preserve
keep id x
tempfile using_data
save `using_data'
restore
cmerge load_using `using_data', verbose
di "  cmerge load_using passed"

* 5. Run csort (different command - should clean cmerge cache)
di "Testing csort (switching from cmerge)..."
gsort -id
csort id
assert id[1] == 1
di "  csort after cmerge passed"

* 6. Run csort multiple times in a loop to check for memory leaks
di "Testing repeated csort (100 iterations)..."
forvalues i = 1/100 {
    gsort -id
    csort id
}
di "  repeated csort passed"

* 7. Final verification
assert id[1] == 1
assert id[_N] == 1000
count
assert r(N) == 1000

di ""
di "========================================="
di "  CLEANUP SYSTEM TESTS PASSED"
di "========================================="
