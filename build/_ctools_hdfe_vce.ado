*! version 1.2.1 20261029 github.com/mdroste/stata-ctools
*! vce() of the HDFE estimators.
* _ctools_hdfe_vce, spec(string) touse(varname) cmd(name)
*
* Classifies spec with _ctools_vce: r(kind) is cluster, robust, unadjusted,
* ols, iid, bootstrap, or empty for anything else. Which kinds and how many
* cluster dimensions a command accepts, and its messages, stay with the
* command. For cluster (vce(cluster a b) or vce(cluster(a b))), each term
* is a dimension: a variable (abbreviations and wildcards expand to one
* dimension per variable), an i.-prefixed variable, or an interaction
* var1#var2 of such variables, as in reghdfe; string variables are allowed.
* Marks out touse for missing cluster values (empty strings included) and
* returns r(ndims) dimensions with, for d = 1..r(ndims):
*   r(name_d)   the dimension as reghdfe names it: turn, turn#foreign
*   r(vars_d)   its variables
*   r(group_d)  1 when they must be numbered with _ctools_group first (a
*               string variable or an interaction), 0 when the numeric
*               variable goes to the plugin as it is
program define _ctools_hdfe_vce, rclass
    version 14.1
    syntax, TOUSE(varname) CMD(name) [SPEC(string)]
    _ctools_vce kind unused : `spec'
    return local kind `kind'
    local n = 0
    if "`kind'" == "cluster" {
        gettoken first terms : spec, parse(" (")
        local terms : subinstr local terms "(" " ", all
        local terms : subinstr local terms ")" " ", all
        foreach term of local terms {
            if !strpos("`term'", "#") & !strpos("`term'", ".") {
                unab expanded : `term', name(vce(cluster))
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
                if substr("`part'", 1, 2) == "i." local part = substr("`part'", 3, .)
                if "`part'" == "" | strpos("`part'", ".") {
                    di as error `"`cmd': vce(cluster) term `term' is not supported"'
                    di as error "specify variables, i.var, or interactions such as var1#var2"
                    exit 198
                }
                unab v : `part', max(1) name(vce(cluster))
                capture confirm numeric variable `v'
                if _rc local group = 1
                local vars `vars' `v'
                if "`name'" != "" local name "`name'#"
                local name "`name'`v'"
            }
            local ++n
            return scalar group_`n' = `group'
            return local name_`n' `name'
            return local vars_`n' `vars'
            markout `touse' `vars', strok
        }
    }
    return scalar ndims = `n'
end
