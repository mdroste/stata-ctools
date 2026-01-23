capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

clear
set obs 3
gen str10 x = "A"
replace x = "B" in 2
replace x = "C" in 3

* After this, x=1, y=2 in dataset
gen long y = .

* Get dataset indices
unab allvars : *
di "allvars = `allvars'"
local x_idx = 1
local y_idx = 2
foreach v of local allvars {
    if "`v'" == "x" local x_idx = `i'
    if "`v'" == "y" local y_idx = `i'
    local ++i
}
di "x_idx=`x_idx', y_idx=`y_idx'"

* Direct plugin call with dataset indices
plugin call ctools_plugin x y, "cencode `x_idx' `y_idx' label=y"
list
