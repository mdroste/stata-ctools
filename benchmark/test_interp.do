clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid

* Exclude basis
gen byte basis = abs(resid) < 1e-8

* Sort non-basis residuals
preserve
keep if !basis
sort resid
gen obs_num = _n

* Get the relevant observations
local N = _N
local q_lo = 0.25958888
local q_hi = 0.74041112

* Different index calculations
local idx_type7_lo = (`N' - 1) * `q_lo'  // Type 7: (n-1)*p, 0-based
local idx_type7_hi = (`N' - 1) * `q_hi'

local idx_ceil_lo = ceil(`N' * `q_lo')   // Ceiling
local idx_ceil_hi = ceil(`N' * `q_hi')

local idx_floor_lo = floor(`N' * `q_lo') // Floor  
local idx_floor_hi = floor(`N' * `q_hi')

local idx_round_lo = round(`N' * `q_lo') // Round
local idx_round_hi = round(`N' * `q_hi')

display "N = " `N'
display _n "=== INDEX CALCULATIONS ==="
display "Type 7 (R default): lo=" %6.2f `idx_type7_lo' " hi=" %6.2f `idx_type7_hi'
display "Ceiling: lo=" `idx_ceil_lo' " hi=" `idx_ceil_hi'
display "Floor: lo=" `idx_floor_lo' " hi=" `idx_floor_hi'
display "Round: lo=" `idx_round_lo' " hi=" `idx_round_hi'

display _n "=== VALUES AT INDICES ==="
display "Stata returns: lo=-731.00 hi=1617.00"
foreach i in 17 18 48 49 50 {
    sum resid if obs_num == `i', meanonly
    display "Obs `i': " %12.4f r(mean)
}

* Type 7 interpolation
local lo_floor = floor(`idx_type7_lo')
local lo_ceil = `lo_floor' + 1
local frac_lo = `idx_type7_lo' - `lo_floor'
sum resid if obs_num == `lo_floor' + 1, meanonly  // 1-based
local val_lo_floor = r(mean)
sum resid if obs_num == `lo_ceil' + 1, meanonly
local val_lo_ceil = r(mean)
local interp_lo = `val_lo_floor' * (1 - `frac_lo') + `val_lo_ceil' * `frac_lo'

local hi_floor = floor(`idx_type7_hi')
local hi_ceil = `hi_floor' + 1
local frac_hi = `idx_type7_hi' - `hi_floor'
sum resid if obs_num == `hi_floor' + 1, meanonly
local val_hi_floor = r(mean)
sum resid if obs_num == `hi_ceil' + 1, meanonly
local val_hi_ceil = r(mean)
local interp_hi = `val_hi_floor' * (1 - `frac_hi') + `val_hi_ceil' * `frac_hi'

display _n "=== INTERPOLATED VALUES ==="
display "Type 7 lo: " %12.4f `interp_lo' " (floor=" %12.4f `val_lo_floor' " ceil=" %12.4f `val_lo_ceil' " frac=" %6.4f `frac_lo' ")"
display "Type 7 hi: " %12.4f `interp_hi' " (floor=" %12.4f `val_hi_floor' " ceil=" %12.4f `val_hi_ceil' " frac=" %6.4f `frac_hi' ")"

* Check what Stata's _pctile actually returns
_pctile resid, percentile(`=`q_lo'*100' `=`q_hi'*100')
display _n "Stata _pctile: lo=" %12.4f r(r1) " hi=" %12.4f r(r2)

restore
