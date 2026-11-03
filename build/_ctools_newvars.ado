*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
*! ctools output-name validation
program define _ctools_newvars
    version 14.1
    syntax anything(name=outputs)
    local seen ""
    foreach name of local outputs {
        if `: list name in seen' {
            di as error "ctools: duplicate output variable `name'"
            exit 198
        }
        confirm new variable `name'
        local seen `seen' `name'
    }
end
