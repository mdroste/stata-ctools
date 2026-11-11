*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! Absorbed effects of creghdfe and civreghdfe.
* _ctools_hdfe_absorb terms, touse(varname) cmd(name) [slopes]
*
* Each term of absorb() is a variable list (abbreviations, wildcards and
* ranges expand to one effect per variable), an i.-prefixed variable (ib#.
* and ibn. as well: the base level does not matter to an absorbed effect),
* or an interaction var1#var2 of such variables, as in reghdfe; string variables
* are allowed; a##b of categorical variables is the interaction a#b, as in
* reghdfe. Without option slopes, heterogeneous slopes (c.), newvar=
* targets and other operators are rejected. Marks out touse for missing
* values (empty strings included) and returns r(n) effects with, for
* g = 1..r(n):
*   r(name_g)   the effect as reghdfe names it: rep78, i.rep78, rep78#foreign
*   r(vars_g)   its variables
*   r(group_g)  1 when they must be numbered with _ctools_group first (a
*               string variable or an interaction), 0 when the numeric
*               variable goes to the plugin as it is
* Option slopes also accepts reghdfe's heterogeneous slopes, id#c.x (slope
* only), id##c.x (intercept and slope), id##c.(x z) and id##(c.x c.z), and
* a newvar= target before any term (with option replace it may name an
* existing variable, which the caller replaces), and returns in addition
*   r(intercept_g)  1 when the effect has an intercept, 0 for slopes only
*   r(cvars_g)      its slope variables, in order (empty without slopes)
*   r(target_g)     the newvar= target, or empty
program define _ctools_hdfe_absorb, rclass
    version 14.1
    syntax anything(name=terms id="absorb()" equalok), TOUSE(varname) CMD(name) [SLOPES REPLACE]
    local n = 0
    if "`slopes'" != "" {
        * Terms are bound by parentheses, c.(x z); a target may be written
        * with spaces around the equal sign, as reghdfe allows
        local terms = ustrregexra(`"`terms'"', "\s*=\s*", "=")
        while `"`terms'"' != "" {
            gettoken term terms : terms, bind
            local target
            if strpos(`"`term'"', "=") {
                gettoken target term : term, parse("=")
                gettoken equal term : term, parse("=")
                if "`target'" == "" | `"`term'"' == "" | strpos(`"`term'"', "=") Unsupported `cmd' `"`target'=`term'"'
                if "`replace'" != "" confirm name `target'
                else confirm new variable `target'
            }
            * a##b of categorical variables is the interaction a#b, as in reghdfe
            if strpos(`"`term'"', "##") & !strpos(`"`term'"', "c.") {
                local term : subinstr local term "##" "#", all
            }
            if strpos(`"`term'"', "c.") {
                SlopeTerm `cmd' `"`term'"'
                local ++n
                return scalar group_`n' = (`: word count `s_ivars'' > 1)
                return local name_`n' `"`s_name'"'
                return local vars_`n' `s_ivars'
                return scalar intercept_`n' = `s_intercept'
                return local cvars_`n' `s_cvars'
                return local target_`n' `target'
                markout `touse' `s_ivars' `s_cvars'
                continue
            }
            local before = `n'
            PlainTerm `cmd' `"`term'"' `touse' `n'
            local n = `r_n'
            forvalues g = `=`before'+1'/`n' {
                return scalar group_`g' = `r_group_`g''
                return local name_`g' `r_name_`g''
                return local vars_`g' `r_vars_`g''
                return scalar intercept_`g' = 1
                return local cvars_`g'
                return local target_`g' `target'
            }
            if "`target'" != "" & `n' != `before' + 1 {
                di as error `"`cmd': absorb() target `target'= must name one absorbed effect"'
                exit 198
            }
        }
        return scalar n = `n'
        exit
    }
    foreach term of local terms {
        * a##b of categorical variables is the interaction a#b, as in reghdfe
        if strpos(`"`term'"', "##") & !strpos(`"`term'"', "c.") {
            local term : subinstr local term "##" "#", all
        }
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
            else if regexm("`part'", "^(ib(n|[0-9]+)[.])(.+)$") {
                local prefix = regexs(1)
                local part = regexs(3)
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

* One term without slopes (the default parser above), for the slopes mode:
* caller locals r_n and r_name_g, r_vars_g, r_group_g for the new effects
program define PlainTerm
    args cmd term touse n
    if strpos(`"`term'"', "##") | strpos(`"`term'"', "(") | strpos(`"`term'"', ",") {
        Unsupported `cmd' `"`term'"'
    }
    _ctools_hdfe_absorb `term', touse(`touse') cmd(`cmd')
    forvalues g = 1/`r(n)' {
        local ++n
        c_local r_name_`n' `r(name_`g')'
        c_local r_vars_`n' `r(vars_`g')'
        c_local r_group_`n' = r(group_`g')
    }
    c_local r_n `n'
end

* A term with heterogeneous slopes. Factor-variable syntax expands it, as
* it does for reghdfe: id##c.(x z) is i.id x z i.id#c.x i.id#c.z. Its
* indicator variables (i.) form the effect; its continuous ones (c.) are the
* slopes; it has an intercept when a term of the expansion has no c. part.
* Caller locals s_name (the term with its variables unabbreviated, as
* e(absvars) shows it), s_ivars, s_cvars and s_intercept.
program define SlopeTerm
    args cmd term
    ExpandTerm `"`term'"'
    local expanded `r_expanded'
    local ivars
    local cvars
    local intercept = 0
    foreach factor of local expanded {
        if !strpos("`factor'", ".") continue
        local parts : subinstr local factor "#" " ", all
        local has_c = 0
        foreach part of local parts {
            local v
            if inlist(substr("`part'", 1, 2), "c.", "i.") local v = substr("`part'", 3, .)
            else if regexm("`part'", "^ib(n|[0-9]+)[.](.+)$") local v = regexs(2)
            if "`v'" == "" | strpos("`v'", ".") Unsupported `cmd' `"`term'"'
            if substr("`part'", 1, 2) == "c." {
                local cvars `cvars' `v'
                local has_c = 1
            }
            else local ivars `ivars' `v'
        }
        if !`has_c' local intercept = 1
    }
    local ivars : list uniq ivars
    local unique : list uniq cvars
    if "`ivars'" == "" | "`cvars'" == "" Unsupported `cmd' `"`term'"'
    if !`: list unique == cvars' {
        di as error `"`cmd': absorb() term `term' repeats a slope variable"'
        exit 198
    }
    c_local s_ivars `ivars'
    c_local s_cvars `cvars'
    c_local s_intercept `intercept'
    c_local s_name `"`r_name'"'
end

* The expansion of a slope term by factor-variable syntax (r_expanded) and
* the term itself with every variable name unabbreviated (r_name)
program define ExpandTerm
    args term
    local rest `"`term'"'
    local name
    while `"`rest'"' != "" {
        gettoken token rest : rest, parse(" #.()")
        local next = substr(`"`rest'"', 1, 1)
        if inlist(`"`token'"', "#", ".", "(", ")") | "`next'" == "." local name `"`name'`token'"'
        else {
            unab vars : `token', name(absorb())
            local name `"`name'`vars'"'
        }
        if "`next'" == " " local name `"`name' "'
    }
    local 0 `"`name'"'
    syntax varlist(numeric fv)
    c_local r_expanded `varlist'
    c_local r_name `"`name'"'
end

program define Unsupported
    args cmd term
    di as error `"`cmd': absorb() term `term' is not supported"'
    di as error "specify variables, i.var, or interactions such as var1#var2"
    exit 198
end
