*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! ctools: drop the scalars and matrices a command and its plugin leave behind
* _ctools_cleanup <command> [<command> ...] drops every scalar and matrix
* named __<command>_* or _<command>_*.
program define _ctools_cleanup
    version 14.1
    foreach command of local 0 {
        foreach pattern in __`command'_* _`command'_* {
            local found : all scalars "`pattern'"
            if "`found'" != "" capture scalar drop `found'
            local found : all matrices "`pattern'"
            if "`found'" != "" capture matrix drop `found'
        }
    }
end
