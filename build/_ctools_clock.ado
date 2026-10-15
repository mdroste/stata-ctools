*! Private monotonic timestamp in seconds; no public Stata timer IDs.
program define _ctools_clock
    version 14.1
    args target
    tempname held
    _return hold `held'
    capture noisily {
        _ctools_load
        capture program ctools_clock_plugin, plugin using("`__ctools_plugin'")
        if _rc != 0 & _rc != 110 exit _rc
        plugin call ctools_clock_plugin, "clock"
    }
    local rc = _rc
    _return restore `held'
    if `rc' exit `rc'
    c_local `target' `ctools_clock_value'
end
