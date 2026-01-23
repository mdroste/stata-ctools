capture do "validation/validate_setup.do"

webuse lifeexp, clear
describe
di "Checking country..."

* Get all variables
unab allvars : *
di "All vars: `allvars'"

* Find country index
local var_idx = 0
local idx = 1
foreach v of local allvars {
    if ("`v'" == "country") {
        local var_idx = `idx'
        di "Found country at index: `var_idx'"
        continue, break
    }
    local ++idx
}

* Try cencode
capture noisily cencode country, generate(test)
di "rc = " _rc
