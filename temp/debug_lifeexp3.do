adopath ++ "./build"
webuse lifeexp, clear

* Replicate what cencode.ado does
local varlist "country"
local generate "country_code"

* Get variable index
unab allvars : *
local var_idx = 0
local idx = 1
foreach v of local allvars {
    if ("`v'" == "`varlist'") {
        local var_idx = `idx'
        continue, break
    }
    local ++idx
}
di "var_idx for `varlist' = `var_idx'"

* Create generate var
quietly generate long `generate' = .

* Get index of new variable
unab allvars : *
local gen_idx = 0
local idx = 1
foreach v of local allvars {
    if ("`v'" == "`generate'") {
        local gen_idx = `idx'
        continue, break
    }
    local ++idx
}
di "gen_idx for `generate' = `gen_idx'"

* Check if source is string
confirm string variable `varlist'
di "Source variable is string: yes (no error)"

di "Arguments to plugin would be:"
di `"cencode  `var_idx' `gen_idx' label=country_code"'
