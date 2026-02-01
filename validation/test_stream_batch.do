* Test streaming mode with batch size parameter
clear all
adopath ++ "build"

sysuse auto, clear
expand 100

di "Test 1: stream(1) - process 1 variable at a time"
csort price, stream(1) verbose

sysuse auto, clear
expand 100

di _n "Test 2: stream(4) - process 4 variables at a time"
csort make, stream(4) verbose

sysuse auto, clear
expand 100

di _n "Test 3: stream(16) - maximum batch size"
csort foreign mpg, stream(16) verbose

sysuse auto, clear
expand 100

di _n "Test 4: stream(8) with multiple sort keys"
csort price mpg, stream(8) verbose

di _n "All streaming batch tests passed!"
