clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/profile_winsor_baseline.log", text replace
program define run_profile
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/baseline"
program perf_clock, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_1m/clock.plugin")
set obs 1000000
set seed 928143
gen double x=rnormal()
gen long g=floor((_n-1)/1000)
gen double shuffle=runiform()
sort shuffle
drop shuffle
quietly cwinsor x, by(g) suffix(_w) threads(8)
drop x_w
di "PROFILE_READY"
plugin call perf_clock, "start"
scalar profile_elapsed = 0
while scalar(profile_elapsed) < 10 {
    quietly cwinsor x, by(g) suffix(_w) threads(8)
    drop x_w
    plugin call perf_clock, "stop"
    scalar profile_elapsed = scalar(__ctools_perf_seconds)
}
end
capture noisily run_profile
local rc=_rc
di "PROFILE_COMPLETE_RC=`rc'"
log close
exit, clear
