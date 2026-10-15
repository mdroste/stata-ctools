*! _ctools_strw — set string width metadata for ctools flat buffer optimization
*!
*! Usage:  _ctools_strw varlist
*!
*! Sets local __ctools_strw in CALLER's scope containing comma-separated
*! string variable widths (0 = numeric, >0 = str width e.g. 17 for str17).
*! The ctools C plugin reads this local automatically via SF_macro_use
*! to enable flat buffer string I/O (eliminates per-string strlen, memcpy,
*! and arena atomic CAS overhead).
*!
*! Example:
*!   _ctools_strw `allvars'
*!   plugin call ctools_plugin `allvars', "mycommand args"
*!   * C plugin automatically detects __ctools_strw local

program define _ctools_strw
    syntax varlist

    local __strw ""
    local __sep ""
    foreach v of local varlist {
        local __type : type `v'
        if substr("`__type'", 1, 3) == "str" {
            local __w = real(substr("`__type'", 4, .))
            if missing(`__w') local __w = 0
        }
        else {
            local __w = 0
        }
        local __strw "`__strw'`__sep'`__w'"
        local __sep ","
    }
    c_local __ctools_strw `__strw'
end
