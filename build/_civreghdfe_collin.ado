*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
*! civreghdfe: ivreg2's collinearity rules for the variable lists
* _civreghdfe_collin, endo() endovars() inexog() inexogvars() exexog()
*     exexogvars() touse() [wvar() wtype() noconstant]
*
* ivreg2's CheckDupsCollin, applied to the variables before the fixed effects
* are partialled out: endo(), inexog() and exexog() hold the estimated terms
* of each list (no base levels, no constant) and *vars() their variables.
* Each pass sweeps invsym() in order over the (weighted, centred unless
* noconstant) cross-products, so a term collinear with earlier ones is
* dropped: first among the endogenous regressors, then over
* [inexog exexog endo]. An endogenous regressor dropped there is collinear
* with the instruments: it becomes the first included exogenous regressor
* and the pass is repeated, which drops the exogenous term it duplicates.
* Returns r(endo), r(inexog), r(exexog): the lists in ivreg2's order with
* dropped terms omitted (o.), r(collin) the dropped terms, r(ecollin) the
* reclassified endogenous regressors and r(changed) (0/1).
program define _civreghdfe_collin, rclass
    version 14.1
    syntax, touse(varname) [endo(string) endovars(string) inexog(string) ///
        inexogvars(string) exexog(string) exexogvars(string) wvar(string) ///
        wtype(string) noCONStant]
    local cons = ("`constant'" == "")
    mata: _civreghdfe_collin_rules()
    return local endo `"`r_endo'"'
    return local inexog `"`r_inexog'"'
    return local exexog `"`r_exexog'"'
    return local collin `"`r_collin'"'
    return local ecollin `"`r_ecollin'"'
    return scalar changed = `r_changed'
end

capture mata: mata drop _civreghdfe_collin_rules()
capture mata: mata drop _civreghdfe_collin_omit()
capture mata: mata drop _civreghdfe_collin_pass()
mata:
/* The omitted form of a term (o.x, 2o.g, co.x#co.z), as _ctools_ebv names it */
string scalar _civreghdfe_collin_omit(string scalar name)
{
    string rowvector parts
    real scalar i, dot

    parts = tokens(subinstr(name, "#", " "))
    for (i = 1; i <= cols(parts); i++) {
        dot = strpos(parts[i], ".")
        if (dot) parts[i] = substr(parts[i], 1, dot - 1) + "o" + substr(parts[i], dot, .)
        else parts[i] = "o." + parts[i]
    }
    return(invtokens(parts, "#"))
}

/* ivreg2_rmcollright2: invsym() swept in order over the columns idx; 1 marks
   a dropped column */
real rowvector _civreghdfe_collin_pass(real matrix XX, real rowvector idx)
{
    if (cols(idx) == 0) return(J(1, 0, 0))
    return((diagonal(invsym(XX[idx, idx], 1..cols(idx))) :== 0)')
}

void _civreghdfe_collin_rules()
{
    string rowvector en, in, xn, vars, names, ecollin, collin
    real rowvector ei, ii, xi, keep_e, om, idx
    real matrix X, XX
    real colvector w
    real scalar ne, ni, nx, i, changed

    en = tokens(st_local("endo"));   ne = cols(en)
    in = tokens(st_local("inexog")); ni = cols(in)
    xn = tokens(st_local("exexog")); nx = cols(xn)
    vars = (tokens(st_local("inexogvars")), tokens(st_local("exexogvars")), tokens(st_local("endovars")))
    names = (in, xn, en)
    /* column positions of the three lists in X */
    ii = (ni ? 1..ni : J(1, 0, 0))
    xi = (nx ? (ni + 1)..(ni + nx) : J(1, 0, 0))
    ei = (ne ? (ni + nx + 1)..(ni + nx + ne) : J(1, 0, 0))

    X = st_data(., vars, st_local("touse"))
    w = 1
    if (st_local("wvar") != "") {
        w = st_data(., st_local("wvar"), st_local("touse"))
        if (st_local("wtype") == "aweight" | st_local("wtype") == "pweight") w = w * rows(w) / quadsum(w)
    }
    if (strtoreal(st_local("cons"))) XX = quadcrossdev(X, mean(X, w), w, X, mean(X, w))
    else XX = quadcross(X, w, X)

    collin = J(1, 0, "")
    ecollin = J(1, 0, "")
    om = J(1, cols(names), 0)
    if (ne == 0) {
        idx = (ii, xi)
        om[idx] = _civreghdfe_collin_pass(XX, idx)
    }
    else {
        /* 1: collinear endogenous regressors */
        om[ei] = _civreghdfe_collin_pass(XX, ei)
        keep_e = select(ei, !om[ei])
        if (cols(keep_e) == 0) keep_e = J(1, 0, 0)
        /* 2: [inexog exexog endo]; an endogenous regressor dropped here is
           reclassified as the first exogenous regressor and the pass redone */
        idx = (ii, xi, keep_e)
        om[idx] = _civreghdfe_collin_pass(XX, idx)
        if (any(om[keep_e])) {
            idx = select(keep_e, om[keep_e])
            ecollin = names[idx]
            ei = (cols(ei) ? select(ei, (rowsum(ei' :== idx) :== 0)') : ei)
            if (cols(ei) == 0) ei = J(1, 0, 0)
            ii = (idx, ii)
            keep_e = select(keep_e, !om[keep_e])
            if (cols(keep_e) == 0) keep_e = J(1, 0, 0)
            om[ii] = J(1, cols(ii), 0)
            idx = (ii, xi, keep_e)
            om[idx] = _civreghdfe_collin_pass(XX, idx)
        }
    }
    changed = any(om) | cols(ecollin) > 0
    for (i = 1; i <= cols(names); i++) {
        if (om[i]) {
            collin = (collin, names[i])
            names[i] = _civreghdfe_collin_omit(names[i])
        }
    }
    st_local("r_endo", cols(ei) ? invtokens(names[ei]) : "")
    st_local("r_inexog", cols(ii) ? invtokens(names[ii]) : "")
    st_local("r_exexog", cols(xi) ? invtokens(names[xi]) : "")
    st_local("r_collin", invtokens(collin))
    st_local("r_ecollin", invtokens(ecollin))
    st_local("r_changed", strofreal(changed))
}
end
