*! ctools: drop the scalars and matrices a command and its plugin leave behind
* _ctools_cleanup <command> [<command> ...] drops every scalar and matrix
* named __<command>_* or _<command>_*.
program define _ctools_cleanup
    version 14.1
    foreach command of local 0 {
        foreach pattern in __`command'_* _`command'_* {
            local found : all scalars "`pattern'"
            foreach name of local found {
                capture scalar drop `name'
            }
            local found : all matrices "`pattern'"
            foreach name of local found {
                capture matrix drop `name'
            }
        }
    }
end
