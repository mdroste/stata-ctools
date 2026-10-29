/* Regression contracts for the HDFE wrapper fixes of 27 Oct 2026. */
do validation/validate_setup.do

capture noisily {
    sysuse auto, clear
    creghdfe price weight, absorb(rep78 foreign)
    assert e(N_hdfe) == 2 & e(N_hdfe_extended) == 2
    assert "`e(vcetype)'" == ""
    creghdfe price weight, absorb(rep78) vce(robust)
    assert "`e(vcetype)'" == "Robust" & "`e(vce)'" == "robust"
    creghdfe price weight, absorb(rep78) vce(cluster turn)
    assert "`e(vcetype)'" == "Robust" & "`e(vce)'" == "cluster"
    cppmlhdfe price weight, absorb(rep78)
    assert "`e(vcetype)'" == "Robust"
    cppmlhdfe price weight, absorb(rep78) vce(cluster turn)
    assert "`e(vcetype)'" == "Robust"
    * The first-stage regressions run by creghdfe carry the label as well
    civreghdfe price (weight = length), absorb(rep78) vce(robust) savefirst
    estimates restore _civreghdfe_weight
    assert "`e(vcetype)'" == "Robust"
    capture which reghdfe
    if !_rc {
        reghdfe price weight, absorb(rep78 foreign)
        assert e(N_hdfe_extended) == 2
    }
}
if _rc test_fail "N_hdfe_extended and Robust vcetype" "r(`=_rc')"
else test_pass "N_hdfe_extended and Robust vcetype"

capture noisily {
    sysuse auto, clear
    gen str5 srep = string(rep78) if !missing(rep78)
    gen str1 sfor = cond(foreign, "F", "D")
    gen int cell = 10 * rep78 + foreign
    creghdfe price weight, absorb(rep78 foreign) vce(cluster turn)
    matrix b0 = e(b)
    matrix V0 = e(V)
    local N0 = e(N)
    creghdfe price weight, absorb(srep sfor) vce(cluster turn)
    assert mreldif(b0, e(b)) < 1e-12 & mreldif(V0, e(V)) < 1e-12 & e(N) == `N0'
    assert "`e(absorb)'" == "srep sfor"
    creghdfe price weight, absorb(i.rep78 i.foreign) vce(cluster turn)
    assert mreldif(b0, e(b)) < 1e-12 & "`e(absorb)'" == "i.rep78 i.foreign"
    creghdfe price weight length, absorb(cell)
    matrix b1 = e(b)
    matrix V1 = e(V)
    creghdfe price weight length, absorb(srep#sfor)
    assert mreldif(b1, e(b)) < 1e-12 & mreldif(V1, e(V)) < 1e-12
    assert "`e(absorb)'" == "srep#sfor" & e(N_hdfe) == 1
    creghdfe price weight length, absorb(i.rep78#i.foreign)
    assert mreldif(b1, e(b)) < 1e-12 & "`e(absorb)'" == "i.rep78#i.foreign"
    capture which reghdfe
    if !_rc {
        reghdfe price weight length, absorb(rep78#foreign)
        assert mreldif(b1, e(b)) < 1e-10 & mreldif(V1, e(V)) < 1e-10
    }
    civreghdfe price (weight = length) turn, absorb(rep78 foreign) vce(cluster turn)
    matrix b2 = e(b)
    matrix V2 = e(V)
    civreghdfe price (weight = length) turn, absorb(srep sfor) vce(cluster turn)
    assert mreldif(b2, e(b)) < 1e-12 & mreldif(V2, e(V)) < 1e-12
    assert "`e(absorb)'" == "srep sfor"
    civreghdfe price (weight = length) turn, absorb(cell)
    matrix b3 = e(b)
    matrix V3 = e(V)
    civreghdfe price (weight = length) turn, absorb(srep#sfor)
    assert mreldif(b3, e(b)) < 1e-12 & mreldif(V3, e(V)) < 1e-12
    assert "`e(absorb)'" == "srep#sfor"
}
if _rc test_fail "string, i. and # absorbed effects" "r(`=_rc')"
else test_pass "string, i. and # absorbed effects"

capture noisily {
    sysuse auto, clear
    foreach bad in "c.length#i.rep78" "rep78##foreign" "L.rep78" "fe=rep78" {
        capture creghdfe price weight, absorb(`bad')
        assert _rc == 198
        capture civreghdfe price (weight = length), absorb(`bad')
        assert _rc == 198
    }
    * savefe is reghdfe's only absorb() suboption; civreghdfe takes none
    capture creghdfe price weight, absorb(rep78, keepsingletons)
    assert _rc == 198
    capture civreghdfe price (weight = length), absorb(rep78, savefe)
    assert _rc == 198
    creghdfe price weight, absorb(rep78, savefe)
    confirm numeric variable __hdfe1__
}
if _rc test_fail "unsupported absorb() terms and suboptions" "r(`=_rc')"
else test_pass "unsupported absorb() terms and suboptions"

capture noisily {
    sysuse auto, clear
    civreghdfe price (weight = length) if foreign == 0 [aw=mpg], absorb(rep78) vce(robust) small
    local cmdline `"`e(cmdline)'"'
    assert strpos(`"`cmdline'"', "civreghdfe price (weight = length)") == 1
    assert strpos(`"`cmdline'"', "if foreign == 0") > 0
    assert strpos(`"`cmdline'"', "[aw=mpg]") > 0
    assert strpos(`"`cmdline'"', "vce(robust) small") > 0
}
if _rc test_fail "civreghdfe e(cmdline) keeps if, weights and options" "r(`=_rc')"
else test_pass "civreghdfe e(cmdline) keeps if, weights and options"

capture noisily {
    sysuse auto, clear
    gen byte cl2 = mod(_n, 7)
    creghdfe price weight, absorb(rep78)
    foreach option in first rf savefirst saverf {
        capture civreghdfe price (weight = length), absorb(rep78) vce(cluster turn cl2) `option'
        assert _rc == 198
        assert "`e(cmd)'" == "creghdfe"
    }
    civreghdfe price (weight = length), absorb(rep78) vce(cluster turn cl2)
    assert "`e(vcetype)'" == "Two-way Cluster"
}
if _rc test_fail "two-way clustering refuses auxiliary regressions up front" "r(`=_rc')"
else test_pass "two-way clustering refuses auxiliary regressions up front"

capture noisily {
    sysuse auto, clear
    cppmlhdfe price weight, absorb(rep78 foreign) dof(pairwise clusters)
    matrix b0 = e(b)
    matrix V0 = e(V)
    local df_a = e(df_a)
    assert "`e(dofmethod)'" == "pairwise clusters"
    cppmlhdfe price weight, absorb(rep78 foreign) dof(pair clus)
    assert mreldif(b0, e(b)) == 0 & mreldif(V0, e(V)) == 0 & e(df_a) == `df_a'
    assert "`e(dofmethod)'" == "pairwise clusters"
    cppmlhdfe price weight, absorb(rep78 foreign) dof(firstp)
    assert "`e(dofmethod)'" == "firstpair"
    cppmlhdfe price weight, absorb(rep78 foreign)
    assert "`e(dofmethod)'" == "pairwise clusters continuous"
    foreach bad in "pair first" "all clusters" "junk" {
        capture cppmlhdfe price weight, absorb(rep78 foreign) dof(`bad')
        assert _rc == 198
    }
    creghdfe price weight, absorb(rep78 foreign) dof(pair clus)
}
if _rc test_fail "dofadjustments() abbreviations" "r(`=_rc')"
else test_pass "dofadjustments() abbreviations"

capture noisily {
    sysuse auto, clear
    civreghdfe price (weight = length), absorb(rep78) noreturn
    mata: st_local("left", invtokens((st_dir("global", "numscalar", "__civreghdfe_*") \ ///
        st_dir("global", "matrix", "__civreghdfe_*"))'))
    assert "`left'" == ""
}
if _rc test_fail "civreghdfe noreturn drops its plugin scalars" "r(`=_rc')"
else test_pass "civreghdfe noreturn drops its plugin scalars"

capture noisily {
    sysuse auto, clear
    gen byte cl2 = mod(_n, 7)
    gen str3 smk = substr(make, 1, 3)
    gen int cell = 100 * turn + foreign
    * String, interaction and i. cluster dimensions, as reghdfe names them
    creghdfe price weight, absorb(rep78) vce(cluster cell)
    matrix b0 = e(b)
    matrix V0 = e(V)
    local G0 = e(N_clust)
    creghdfe price weight, absorb(rep78) vce(cluster turn#foreign)
    assert mreldif(b0, e(b)) < 1e-12 & mreldif(V0, e(V)) < 1e-12 & e(N_clust) == `G0'
    assert "`e(clustvar)'" == "turn#foreign" & "`e(clustvar1)'" == "turn#foreign"
    creghdfe price weight, absorb(rep78) vce(cluster i.turn)
    assert "`e(clustvar)'" == "turn" & e(N_clustervars) == 1 & e(N_clust1) == e(N_clust)
    creghdfe price weight, absorb(rep78) vce(cluster smk)
    assert "`e(clustvar)'" == "smk" & e(N_clust) == 22
    matrix bs = e(b)
    matrix Vs = e(V)
    civreghdfe price (weight = length), absorb(rep78) vce(cluster cell)
    matrix b1 = e(b)
    matrix V1 = e(V)
    civreghdfe price (weight = length), absorb(rep78) vce(cluster turn#foreign)
    assert mreldif(b1, e(b)) < 1e-12 & mreldif(V1, e(V)) < 1e-12
    civreghdfe price (weight = length), absorb(rep78) vce(cluster turn cl2)
    assert "`e(clustvar)'" == "turn cl2" & "`e(clustvar1)'" == "turn" & "`e(clustvar2)'" == "cl2"
    cppmlhdfe price weight, absorb(rep78) vce(cluster i.turn foreign#rep78)
    assert "`e(clustvar)'" == "turn foreign#rep78" & e(N_clustervars) == 2
    assert e(N_clust) == min(e(N_clust1), e(N_clust2))
    foreach bad in "L.turn" "c.turn" {
        capture creghdfe price weight, absorb(rep78) vce(cluster `bad')
        assert _rc == 198
    }
    * pweights need robust or clustered standard errors, as in reghdfe
    capture creghdfe price weight [pw=mpg], absorb(rep78) vce(ols)
    assert _rc == 9
    creghdfe price weight [pw=mpg], absorb(rep78)
    assert "`e(vce)'" == "robust"
    capture which reghdfe
    if !_rc {
        reghdfe price weight, absorb(rep78) vce(cluster smk)
        assert mreldif(e(b), bs) < 1e-10 & mreldif(e(V), Vs) < 1e-10 & e(N_clust) == 22
    }
}
if _rc test_fail "shared vce() parsing" "r(`=_rc')"
else test_pass "shared vce() parsing"

capture noisily {
    sysuse auto, clear
    creghdfe price weight length, absorb(i.rep78 foreign#turn) dof(firstpair clusters) vce(cluster turn)
    assert "`e(absvars)'" == "i.rep78 foreign#turn"
    assert "`e(extended_absvars)'" == "rep78 foreign#turn"
    assert "`e(title)'" == "HDFE Linear regression" & "`e(title2)'" == "Absorbing 2 HDFE groups"
    assert "`e(dofmethod)'" == "firstpair clusters"
    assert e(N_hdfe_extended) == 2 & e(drop_singletons) == 1
    matrix t = e(dof_table)
    assert rowsof(t) == 2 & colsof(t) == 5
    assert e(df_a_initial) == t[1,1] + t[2,1] & e(df_a_redundant) == t[1,2] + t[2,2]
    assert t[2,5] == 1 & missing(t[1,4])
    local c_initial = e(df_a_initial)
    local c_redundant = e(df_a_redundant)
    local c_df_a = e(df_a)
    creghdfe price weight
    assert "`e(absvars)'" == "_cons" & e(df_a_initial) == 1 & e(df_a_redundant) == 0
    cppmlhdfe price weight, absorb(rep78 foreign)
    matrix t = e(dof_table)
    local columns : colnames t
    assert "`columns'" == "Categories Redundant Num_Coefs Inexact Nested"
    assert e(df_a_initial) == t[1,1] + t[2,1] & e(df_a_redundant) == t[1,2] + t[2,2]
    capture which reghdfe
    if !_rc {
        reghdfe price weight length, absorb(i.rep78 foreign#turn) dof(firstpair clusters) vce(cluster turn)
        assert e(df_a_initial) == `c_initial' & e(df_a_redundant) == `c_redundant' & e(df_a) == `c_df_a'
        * Cluster variables are marked out before the factor variables are
        * expanded, so both commands name the same coefficients
        gen cl3 = turn if foreign == 0
        reghdfe price weight i.foreign, absorb(rep78) vce(cluster cl3)
        local names_r : colfullnames e(b)
        creghdfe price weight i.foreign, absorb(rep78) vce(cluster cl3)
        local names_c : colfullnames e(b)
        assert "`names_r'" == "`names_c'"
    }
}
if _rc test_fail "shared reghdfe e() results" "r(`=_rc')"
else test_pass "shared reghdfe e() results"

print_summary oct27
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "oct27"
