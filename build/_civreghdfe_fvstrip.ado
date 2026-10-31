*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! Private helper for civreghdfe.
program define _civreghdfe_fvstrip, rclass
    args s
    while regexm(`"`s'"', "([0-9]+)(bn|b|o)\.") {
        local mk = regexs(0)
        local lev = regexs(1)
        local s : subinstr local s "`mk'" "`lev'.", all
    }
    local s : subinstr local s "#co." "#c.", all
    local s : subinstr local s "#o." "#", all
    if substr(`"`s'"', 1, 3) == "co." local s = "c." + substr(`"`s'"', 4, .)
    else if substr(`"`s'"', 1, 2) == "o." local s = substr(`"`s'"', 3, .)
    return local name `"`s'"'
end
