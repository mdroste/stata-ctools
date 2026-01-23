* Load setup
capture do "validation/validate_setup.do"

* Check what variable index we're passing
sysuse census, clear
unab allvars : *
di "All vars: `allvars'"

local var_idx = 0
local idx = 1
foreach v of local allvars {
    if ("`v'" == "state2") {
        local var_idx = `idx'
        di "Found state2 at index: `var_idx'"
        continue, break
    }
    local ++idx
}
di "Final var_idx for state2 = `var_idx'"

* What's at that index?
di "Variable 1 is: " _byindex(1)
di "Variable `var_idx' is: " _byindex(`var_idx')
