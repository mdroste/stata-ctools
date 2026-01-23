* Test basic cqreg without HDFE
capture adopath ++ "build"

sysuse auto, clear
di "Testing basic cqreg..."
cqreg price mpg weight
di "Basic test passed!"

di "Testing with timeit option..."
cqreg price mpg weight, timeit
di "Timeit test passed!"

di "Testing vce(robust)..."
cqreg price mpg weight, vce(robust)
di "Robust test passed!"

di "All tests completed!"
