* Debug with simple data

version 14.0
clear all
set more off

adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

* Create simple master data (10 rows, 3 groups)
clear
input int id int grp float value
1 2 10
2 1 20
3 3 30
4 1 40
5 2 50
6 3 60
7 1 70
8 2 80
9 3 90
10 1 100
end

di "=== ORIGINAL MASTER DATA ==="
list, sep(0)

tempfile master
quietly save `master'

* Create using data
clear
input int grp str10 label
1 "GroupA"
2 "GroupB"
3 "GroupC"
end
tempfile using_data
quietly save `using_data'

di ""
di "=== STATA m:1 MERGE ==="
use `master', clear
merge m:1 grp using `using_data'
list id grp value label _merge, sep(0)

di ""
di "=== CMERGE m:1 (verbose) ==="
use `master', clear
cmerge m:1 grp using `using_data', verbose
list id grp value label _merge, sep(0)
