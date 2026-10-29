*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! ctools: seconds since the previous lap of a caller-local wall clock
* _ctools_lap <target> <clock> <clock value> stores the elapsed seconds in the
* caller's local <target> and restarts the clock, a local the caller starts
* with _ctools_clock <clock> and passes back expanded:
*     _ctools_clock __ct_lap
*     ...
*     _ctools_lap __time_setup __ct_lap `__ct_lap'
program define _ctools_lap
    version 14.1
    args target clock start
    quietly _ctools_clock now
    c_local `target' = `now' - `start'
    c_local `clock' `now'
end
