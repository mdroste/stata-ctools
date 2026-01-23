* Simple test to debug creghdfe crash
clear
sysuse auto, clear
di "Running creghdfe with verbose mode..."
creghdfe price mpg weight, absorb(foreign) verbose
di "Test completed successfully"
