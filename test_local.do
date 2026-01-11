adopath + "./build"

* Create local test files instead of web
sysuse auto, clear
keep make price mpg
tempfile master
save `master'

sysuse auto, clear
keep make weight length
tempfile using_file
save `using_file'

use `master', clear
di "Master loaded: " _N " obs"
desc, short

di "Calling cmerge with string key..."
cmerge 1:1 make using `using_file', verbose
di "cmerge completed!"
list in 1/5
