*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! dofadjustments() of the HDFE estimators.
* _ctools_hdfe_dof [, all none firstpair pairwise clusters continuous]
*
* all (the default), none, or any of firstpair|pairwise, clusters and
* continuous, with reghdfe's abbreviations and exclusions; all stands for
* pairwise clusters continuous. Returns
*   s(flags)     pairwise = 1, firstpair = 2, clusters = 4, continuous = 8
*   s(methods)   the selected adjustments in that order, as e(dofmethod)
*   s(type)      0 = all, 1 = none, 2 = firstpair, 3 = pairwise
*   s(clusters)  1 when the cluster-nesting adjustment applies
program define _ctools_hdfe_dof, sclass
    version 14.1
    syntax [, ALL NONE FIRSTpair PAIRwise CLusters CONTinuous]
    if "`none'`firstpair'`pairwise'`clusters'`continuous'" == "" local all all
    opts_exclusive "`all' `none' `firstpair' `pairwise'" dofadjustments
    opts_exclusive "`all' `none' `clusters'" dofadjustments
    opts_exclusive "`all' `none' `continuous'" dofadjustments
    if "`all'" != "" {
        local pairwise pairwise
        local clusters clusters
        local continuous continuous
    }
    sreturn local flags = ("`pairwise'" != "") + 2 * ("`firstpair'" != "") + ///
        4 * ("`clusters'" != "") + 8 * ("`continuous'" != "")
    sreturn local methods `pairwise' `firstpair' `clusters' `continuous'
    sreturn local type = cond("`all'" != "", 0, cond("`firstpair'" != "", 2, cond("`pairwise'" != "", 3, 1)))
    sreturn local clusters = "`all'`clusters'" != ""
end
