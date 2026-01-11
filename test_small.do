adopath + "./build"

* Use small dataset like webuse autosize
webuse autosize, clear
desc, short
tempfile master
save `master'

webuse autoexpense, clear
desc, short
tempfile using_file
save `using_file'

use `master', clear
di "Master loaded: " _N " obs"

di "Calling cmerge..."
cmerge 1:1 make using `using_file', verbose
di "cmerge completed!"
list
