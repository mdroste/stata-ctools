* Debug Test 19b - More detailed comparison

version 14.0
clear all
set more off

adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

sysuse census, clear
keep state region pop
tempfile master
quietly save `master'

clear
input byte region str20 region_name
1 "NE"
2 "N Cntrl"
3 "South"
4 "West"
end
tempfile using_data
quietly save `using_data'

di "=== STATA MERGE (first 15 rows) ==="
use `master', clear
merge m:1 region using `using_data'
list state region pop region_name _merge in 1/15, sep(0)

tempfile stata_merged
quietly save `stata_merged'

di ""
di "=== CMERGE (first 15 rows) ==="
use `master', clear
cmerge m:1 region using `using_data'
list state region pop region_name _merge in 1/15, sep(0)

tempfile cmerge_merged
quietly save `cmerge_merged'
