*! e() results shared by creghdfe and cpplmhdfe; maintained in package/ado.
* _ctools_hdfe_post, stub(name) absvars(string) nhdfe(#) vce(name)
*     [extended(string) clustvars(string) nclust(numlist) dofmethod(string)]
*
* Call after ereturn post, while the plugin's per-effect scalars
* <stub>_num_levels_g, <stub>_fe_redundant_g and <stub>_fe_nested_g (and
* <stub>_fe_exact_g where the plugin knows it) exist, g = 1..the number of
* words in extended(). Posts, as reghdfe does:
*   e(absvars), e(extended_absvars) (default: absvars without i. prefixes),
*   e(N_hdfe) = nhdfe, e(N_hdfe_extended)
*   e(dof_table): Categories, Redundant, Num_Coefs, Inexact (missing where
*   unknown) and Nested per effect; e(df_a_initial) and e(df_a_redundant)
*   e(dofmethod); e(vce), and e(vcetype) Robust unless vce is unadjusted
*   with clustvars: e(clustvar), e(N_clustervars), e(clustvar#), e(N_clust#)
*   from nclust (one count per variable), e(N_clust) their minimum
program define _ctools_hdfe_post, eclass
    version 14.1
    syntax, STUB(name) ABSvars(string) NHDFE(integer) VCE(name) ///
        [EXTended(string) CLUSTvars(string) NCLUST(numlist integer) DOFmethod(string)]
    if `"`extended'"' == "" local extended : subinstr local absvars "i." "", all
    local nfe : word count `extended'
    tempname table
    matrix `table' = J(`nfe', 5, .)
    local initial = 0
    local redundant = 0
    forvalues g = 1/`nfe' {
        matrix `table'[`g', 1] = scalar(`stub'_num_levels_`g')
        matrix `table'[`g', 2] = scalar(`stub'_fe_redundant_`g')
        matrix `table'[`g', 3] = `table'[`g', 1] - `table'[`g', 2]
        capture confirm scalar `stub'_fe_exact_`g'
        if !_rc matrix `table'[`g', 4] = 1 - scalar(`stub'_fe_exact_`g')
        matrix `table'[`g', 5] = scalar(`stub'_fe_nested_`g')
        local initial = `initial' + `table'[`g', 1]
        local redundant = `redundant' + `table'[`g', 2]
    }
    matrix colnames `table' = Categories Redundant Num_Coefs Inexact Nested
    matrix rownames `table' = `extended'

    ereturn local absvars `"`absvars'"'
    ereturn local extended_absvars `"`extended'"'
    ereturn scalar N_hdfe = `nhdfe'
    ereturn scalar N_hdfe_extended = `nfe'
    ereturn scalar df_a_initial = `initial'
    ereturn scalar df_a_redundant = `redundant'
    ereturn matrix dof_table = `table'
    ereturn local dofmethod `"`dofmethod'"'
    ereturn local vce `vce'
    if "`vce'" != "unadjusted" ereturn local vcetype "Robust"

    local ndims : word count `clustvars'
    if `ndims' {
        if `: word count `nclust'' != `ndims' {
            di as error "_ctools_hdfe_post: one cluster count per cluster variable"
            exit 198
        }
        local fewest .
        forvalues d = 1/`ndims' {
            local count : word `d' of `nclust'
            ereturn local clustvar`d' `: word `d' of `clustvars''
            ereturn scalar N_clust`d' = `count'
            local fewest = min(`fewest', `count')
        }
        ereturn local clustvar `"`clustvars'"'
        ereturn scalar N_clustervars = `ndims'
        ereturn scalar N_clust = `fewest'
    }
end
