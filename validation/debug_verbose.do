* Debug with verbose

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

use `master', clear
cmerge m:1 region using `using_data', verbose
list state region pop region_name _merge in 1/10, sep(0)
