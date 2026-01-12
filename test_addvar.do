clear all
set obs 100

di "Before javacall:"
desc

javacall AddVar3 create, classpath(/Users/Mike/Documents/GitHub/stata-ctools) args(double testvar)

di "After javacall:"
desc
