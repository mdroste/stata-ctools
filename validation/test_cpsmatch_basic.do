* test_cpsmatch_basic.do - Basic cpsmatch validation
clear all
set seed 12345
discard

* Add build directory to adopath
local cwd = c(pwd)
adopath ++ "`cwd'/build"

* Load the plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

* Create test data
set obs 200
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

di ""
di "=========================================="
di "TEST 1: Basic nearest neighbor matching"
di "=========================================="
cpsmatch treat x1 x2, outcome(y)
assert _pscore > 0 & _pscore < 1
assert _support == 1

di ""
di "=========================================="
di "TEST 2: Multiple neighbors"
di "=========================================="
cpsmatch treat x1 x2, neighbor(3) outcome(y)
assert _pscore > 0 & _pscore < 1

di ""
di "=========================================="
di "TEST 3: Kernel matching"
di "=========================================="
cpsmatch treat x1 x2, kernel outcome(y)
assert _pscore > 0 & _pscore < 1

di ""
di "=========================================="
di "TEST 4: Without replacement"
di "=========================================="
cpsmatch treat x1 x2, noreplacement outcome(y)
assert _pscore > 0 & _pscore < 1

di ""
di "=========================================="
di "TEST 5: With caliper"
di "=========================================="
cpsmatch treat x1 x2, caliper(0.1) outcome(y)
assert _pscore > 0 & _pscore < 1

di ""
di "=========================================="
di "All tests passed!"
di "=========================================="
