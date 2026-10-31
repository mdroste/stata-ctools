*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! Absorbed effects of creghdfe and civreghdfe.
* _ctools_hdfe_absorb terms, touse(varname) cmd(name)
*
* Each term of absorb() is a variable list (abbreviations, wildcards and
* ranges expand to one effect per variable), an i.-prefixed variable, or an
* interaction var1#var2 of such variables, as in reghdfe; string variables
* are allowed. Heterogeneous slopes (c.), ##, newvar= targets and other
* operators are rejected. Marks out touse for missing values (empty strings
* included) and returns r(n) effects with, for g = 1..r(n):
*   r(name_g)   the effect as reghdfe names it: rep78, i.rep78, rep78#foreign
*   r(vars_g)   its variables
*   r(group_g)  1 when they must be numbered with _ctools_group first (a
*               string variable or an interaction), 0 when the numeric
*               variable goes to the plugin as it is
program define _ctools_hdfe_absorb, rclass
    version 14.1
    syntax anything(name=terms id="absorb()" equalok), TOUSE(varname) CMD(name)
    local n = 0
    foreach term of local terms {
        if strpos(`"`term'"', "c.") | strpos(`"`term'"', "##") | strpos(`"`term'"', "=") | ///
            strpos(`"`term'"', "(") | strpos(`"`term'"', ",") {
            Unsupported `cmd' `"`term'"'
        }
        if !strpos("`term'", "#") & !strpos("`term'", ".") {
            unab expanded : `term', name(absorb())
            foreach v of local expanded {
                local ++n
                capture confirm numeric variable `v'
                return scalar group_`n' = (_rc != 0)
                return local name_`n' `v'
                return local vars_`n' `v'
                markout `touse' `v', strok
            }
            continue
        }
        local parts : subinstr local term "#" " ", all
        local group = (`: word count `parts'' > 1)
        local name
        local vars
        foreach part of local parts {
            local prefix
            if substr("`part'", 1, 2) == "i." {
                local prefix "i."
                local part = substr("`part'", 3, .)
            }
            if "`part'" == "" | strpos("`part'", ".") Unsupported `cmd' `"`term'"'
            unab v : `part', max(1) name(absorb())
            capture confirm numeric variable `v'
            if _rc local group = 1
            local vars `vars' `v'
            if "`name'" != "" local name "`name'#"
            local name "`name'`prefix'`v'"
        }
        local ++n
        return scalar group_`n' = `group'
        return local name_`n' `name'
        return local vars_`n' `vars'
        markout `touse' `vars', strok
    }
    return scalar n = `n'
end

program define Unsupported
    args cmd term
    di as error `"`cmd': absorb() term `term' is not supported"'
    di as error "specify variables, i.var, or interactions such as var1#var2"
    exit 198
end
