adopath ++ "./build"
webuse lifeexp, clear
di "Variables in dataset:"
describe

* Check what cencode passes to plugin
gen long country_code = .
di "After generate:"
unab allvars : *
local idx = 1
foreach v of local allvars {
    di "  `idx': `v'"
    local ++idx
}
drop country_code
