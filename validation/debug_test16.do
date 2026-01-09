* Debug Test 16 - Merge with no matches

version 14.0
clear all
set more off

adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

clear
input int id float value
1 10
2 20
3 30
end
tempfile master
quietly save `master'

clear
input int id float other
100 1.5
200 2.5
300 3.5
end
tempfile using_data
quietly save `using_data'

di "=== STATA MERGE ==="
use `master', clear
merge 1:1 id using `using_data'
list, sep(0)
describe

tempfile stata_merged
quietly save `stata_merged'

di ""
di "=== CMERGE ==="
use `master', clear
cmerge 1:1 id using `using_data'
list, sep(0)
describe

tempfile cmerge_merged
quietly save `cmerge_merged'

di ""
di "=== COMPARISON ==="
use `stata_merged', clear
cf _all using `cmerge_merged', all verbose
