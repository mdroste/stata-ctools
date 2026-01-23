adopath ++ build

// Create master dataset with extended missing values
clear
input float key
1
.a
.b
.z
end
save temp/master.dta, replace

// Create using dataset with some matching extended missing values
clear
input float key
1
.a
.c
end
save temp/using.dta, replace

// Test Stata's native merge
use temp/master.dta, clear
di "=== Stata native merge ==="
merge 1:1 key using temp/using.dta
tab _merge
list key _merge
save temp/stata_merge.dta, replace

// Test cmerge
use temp/master.dta, clear
di "=== cmerge ==="
cmerge 1:1 key using temp/using.dta
tab _merge
list key _merge

