* Test cqreg with absorb option
capture adopath ++ "build"

sysuse auto, clear
di "Testing absorb(foreign)..."
cqreg price mpg weight, absorb(foreign)

di ""
di "Testing absorb(rep78)..."
sysuse auto, clear
cqreg price mpg weight, absorb(rep78)

di ""
di "Testing with if condition..."
sysuse auto, clear
cqreg price mpg weight if foreign==0, absorb(rep78)

di ""
di "All tests completed successfully!"
