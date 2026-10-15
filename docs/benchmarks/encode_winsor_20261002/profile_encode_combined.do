clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/profile_encode_combined.log", text replace
program define run_profile
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/combined"
program perf_clock, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_1m/clock.plugin")
set obs 200000
gen str12 x="value_"+string(mod(_n-1,10000), "%05.0f")
quietly cencode x, gen(y) threads(8)
drop y
label drop y
di "PROFILE_READY"
plugin call perf_clock, "start"
scalar profile_elapsed = 0
while scalar(profile_elapsed) < 10 {
    quietly cencode x, gen(y) threads(8)
    drop y
label drop y
    plugin call perf_clock, "stop"
    scalar profile_elapsed = scalar(__ctools_perf_seconds)
}
end
capture noisily run_profile
local rc=_rc
di "PROFILE_COMPLETE_RC=`rc'"
log close
exit, clear
