capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

webuse lifeexp, clear
di "Initial variables:"
describe

* Simulate marksample
tempvar touse
gen byte `touse' = 1

di ""
di "After marksample (tempvar created):"
describe

* Find country index
unab allvars : *
di ""
di "allvars = `allvars'"
local var_idx = 0
local idx = 1
foreach v of local allvars {
    if "`v'" == "country" {
        local var_idx = `idx'
        continue, break
    }
    local ++idx
}
di "var_idx for country = `var_idx'"

* Create destination
gen long test = .

di ""
di "After creating test:"
describe

* Find test index
unab allvars : *
di ""
di "allvars = `allvars'"
local gen_idx = 0
local idx = 1
foreach v of local allvars {
    if "`v'" == "test" {
        local gen_idx = `idx'
        continue, break
    }
    local ++idx
}
di "gen_idx for test = `gen_idx'"

di ""
di "Would pass: var_idx=`var_idx', gen_idx=`gen_idx'"
di "country is at dataset position 2, test is at position 8"
di "With marksample tempvar at position 7"
