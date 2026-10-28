/*******************************************************************************
 * validate_hdfe_datasets.do
 *
 * creghdfe vs reghdfe on Sergio Correia's HDFE benchmark collection
 * (https://www.scorreia.com/data/hdfe/): 15 real and synthetic two-way graphs,
 * from complete bipartite to path-like, several with extra identifiers.
 *
 * Opt-in: the collection is about 227 MB and is not part of the offline gate.
 *   python3 scripts/fetch_validation_data.py --collection hdfe
 *   python3 validation/run_hdfe_datasets.py
 *
 * VERIFICATION: every case runs both commands and compares the estimation
 * sample, singleton count, absorbed and residual degrees of freedom exactly,
 * and fit statistics, e(b) and e(V) at $DEFAULT_SIGFIGS. Robust and cluster
 * cases are tested on converged solutions (see hdfe_case); two small graphs
 * and the path graph are also checked against exact direct solutions. Every
 * comparison is written to <outdir>/results.csv and oracle.csv for
 * docs/HDFE_BENCHMARK_DATASETS.md. Timings are single wall-clock runs recorded
 * for context, never asserted.
 ******************************************************************************/

* args: output directory; directory holding the plugin and ado files to test
args outdir plugin_dir
if `"`outdir'"' == "" local outdir "temp/hdfe_datasets"

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}
if `"`plugin_dir'"' != "" adopath ++ `"`plugin_dir'"'

capture mkdir "temp"
capture mkdir "`outdir'"

/*******************************************************************************
 * hdfe_use - Load one dataset from the pinned cache (never downloads)
 ******************************************************************************/
capture program drop hdfe_use
program define hdfe_use
    args dataset
    local path "validation/fixtures/data/hdfe/`dataset'.dta"
    capture confirm file "`path'"
    if _rc local path "fixtures/data/hdfe/`dataset'.dta"
    confirm file "`path'"
    use "`path'", clear
    global HDFE_DATASET "`dataset'"
    di as text "  `dataset' (" as result _N as text " obs)"
end

/*******************************************************************************
 * hdfe_same_partition - Do two group identifiers define the same partition?
 *
 * Returns r(different), rows whose group is split or merged by the other
 * identifier or missing in only one, and r(groups), distinct values of the first.
 ******************************************************************************/
capture program drop hdfe_same_partition
program define hdfe_same_partition, rclass
    args g1 g2
    preserve
    keep `g1' `g2'
    quietly count if missing(`g1') != missing(`g2')
    local different = r(N)
    quietly keep if !missing(`g1', `g2')
    tempvar flag
    foreach pair in "`g1' `g2'" "`g2' `g1'" {
        gettoken a b : pair
        quietly bysort `a' (`b'): gen byte `flag' = `b'[1] != `b'[_N]
        quietly count if `flag'
        local different = `different' + r(N)
        drop `flag'
    }
    quietly bysort `g1': gen byte `flag' = _n == 1
    quietly count if `flag'
    return scalar groups = r(N)
    return scalar different = `different'
    restore
end

/*******************************************************************************
 * hdfe_run - Run reghdfe and creghdfe once each and compare; records no test
 *
 * Exact: e(N), e(num_singletons), e(df_a), e(df_r), e(df_m), e(rank),
 *        e(N_hdfe), e(N_clust). With two effects and no clustering, also the
 *        number of mobility groups: the redundant coefficients of the second
 *        effect in each command's e(dof_table). That count is not the
 *        mobility groups in other models.
 * $DEFAULT_SIGFIGS: e(r2), e(r2_a), e(r2_within), e(rss), e(tss), e(mss),
 *        e(F), e(rmse), and every element of e(b) and e(V).
 *
 * Returns r(differences), empty when everything agrees; r(rc_reghdfe) and
 * r(rc_creghdfe); minimum significant figures r(sf_b), r(sf_v), r(sf_stats);
 * reghdfe's r(n), r(singletons), r(groups), r(iterations); both commands'
 * r(df_a_*) and r(df_r_*); and wall-clock seconds r(t_reghdfe), r(t_creghdfe).
 *
 * Syntax: hdfe_run depvar indepvars, absorb(varlist) [vce(string)
 *             reghdfe(options) creghdfe(options)]
 ******************************************************************************/
capture program drop hdfe_run
program define hdfe_run, rclass
    syntax varlist(min=2 numeric), Absorb(varlist) [vce(string) reghdfe(string) creghdfe(string)]

    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"
    local nfe : word count `absorb'
    local stats_exact "N num_singletons df_a df_r df_m rank N_hdfe N_clust"
    local stats_float "r2 r2_a r2_within rss tss mss F rmse"
    local options_r "`reghdfe'"
    local options_c "`creghdfe'"

    tempname b_r V_r b_c V_c dof
    local r_groups = .
    local c_groups = .
    foreach cmd in r c {
        local command = cond("`cmd'" == "r", "reghdfe", "creghdfe")
        timer clear 91
        timer on 91
        capture quietly `command' `varlist', absorb(`absorb') `vceopt' `options_`cmd''
        local rc_`cmd' = _rc
        timer off 91
        quietly timer list 91
        local t_`cmd' = r(t91)
        if `rc_`cmd'' continue
        matrix `b_`cmd'' = e(b)
        matrix `V_`cmd'' = e(V)
        foreach stat in `stats_exact' `stats_float' ic {
            local `cmd'_`stat' = e(`stat')
        }
        if `nfe' == 2 & !strmatch("`vce'", "cl*") {
            matrix `dof' = e(dof_table)
            local `cmd'_groups = `dof'[2, 2]
        }
    }

    return scalar rc_reghdfe = `rc_r'
    return scalar rc_creghdfe = `rc_c'
    return scalar t_reghdfe = `t_r'
    return scalar t_creghdfe = `t_c'
    if `rc_r' | `rc_c' {
        return local differences "reghdfe rc=`rc_r', creghdfe rc=`rc_c'"
        exit
    }

    local diffs ""
    foreach stat in `stats_exact' groups {
        if !missing(`r_`stat'') {
            if missing(`c_`stat'') local diffs "`diffs' `stat':missing_in_creghdfe"
            else if `r_`stat'' != `c_`stat'' local diffs "`diffs' `stat':`r_`stat''!=`c_`stat''"
        }
    }

    local sf_stats = 15
    foreach stat in `stats_float' {
        if !missing(`r_`stat'') {
            if missing(`c_`stat'') local diffs "`diffs' `stat':missing_in_creghdfe"
            else {
                sigfigs `r_`stat'' `c_`stat''
                local sf = r(sigfigs)
                local sf_stats = min(`sf_stats', `sf')
                if `sf' < $DEFAULT_SIGFIGS {
                    local sf_fmt : display %4.1f `sf'
                    local diffs "`diffs' `stat':sigfigs=`sf_fmt'"
                }
            }
        }
    }

    local sf_b = .
    local sf_V = .
    local names_r : colfullnames `b_r'
    local names_c : colfullnames `b_c'
    if "`names_r'" != "`names_c'" {
        local diffs "`diffs' b:names(`names_r'|`names_c')"
    }
    else {
        foreach m in b V {
            matrix_min_sigfigs ``m'_r' ``m'_c'
            local sf_`m' = r(min_sigfigs)
            if `sf_`m'' < $DEFAULT_SIGFIGS {
                local sf_fmt : display %4.1f `sf_`m''
                local diffs "`diffs' `m':sigfigs=`sf_fmt'"
            }
        }
    }

    return scalar n = `r_N'
    return scalar singletons = `r_num_singletons'
    return scalar groups = `r_groups'
    return scalar iterations = `r_ic'
    return scalar df_a_reghdfe = `r_df_a'
    return scalar df_a_creghdfe = `c_df_a'
    return scalar df_r_reghdfe = `r_df_r'
    return scalar df_r_creghdfe = `c_df_r'
    return scalar sf_b = `sf_b'
    return scalar sf_v = `sf_V'
    return scalar sf_stats = `sf_stats'
    return local differences "`=trim("`diffs'")'"
end

/*******************************************************************************
 * hdfe_case - One test: reghdfe vs creghdfe with their default solver settings
 *
 * groups: also request groupvar() from both and require the same partition
 *     (a second test).
 * converged: for sandwich covariances. The robust and cluster "meat" is
 *     first-order sensitive to projection error, so at tolerance(1e-8) both
 *     commands carry it to only five to seven digits on weakly connected
 *     graphs, and reghdfe's default accelerated projections stop improving
 *     at six to eight. The default-settings comparison is then only
 *     recorded, and the test is the comparison of converged solutions under
 *     the same $DEFAULT_SIGFIGS rule: $HDFE_CONVERGED_REGHDFE against
 *     $HDFE_CONVERGED_CREGHDFE.
 *
 * Syntax: hdfe_case depvar indepvars, absorb(varlist) [vce(string) groups converged]
 ******************************************************************************/
global HDFE_CONVERGED_REGHDFE "technique(lsmr) tolerance(1e-13)"
global HDFE_CONVERGED_CREGHDFE "tolerance(1e-12)"

capture program drop hdfe_case
program define hdfe_case
    syntax varlist(min=2 numeric), Absorb(varlist) [vce(string) GROUPs CONVerged]

    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"
    local testname = trim("$HDFE_DATASET: `varlist', absorb(`absorb') `vceopt'")
    local nfull = _N
    local modes "default"
    if "`converged'" != "" local modes "default converged"
    local asserted : word `: word count `modes'' of `modes'
    local numbers "n singletons groups df_a_reghdfe df_a_creghdfe df_r_reghdfe df_r_creghdfe"
    local numbers "`numbers' sf_b sf_v sf_stats iterations t_reghdfe t_creghdfe"

    foreach mode of local modes {
        local options_r ""
        local options_c ""
        if "`mode'" == "converged" {
            local options_r "$HDFE_CONVERGED_REGHDFE"
            local options_c "$HDFE_CONVERGED_CREGHDFE"
        }
        else if "`groups'" != "" {
            local options_r "groupvar(__hdfe_group_r)"
            local options_c "groupvar(__hdfe_group_c)"
        }
        foreach cmd in r c {
            local pass_`cmd' ""
            if "`options_`cmd''" != "" {
                local pass_`cmd' = cond("`cmd'" == "r", "reghdfe", "creghdfe") + "(`options_`cmd'')"
            }
        }
        hdfe_run `varlist', absorb(`absorb') `vceopt' `pass_r' `pass_c'
        local differences "`r(differences)'"
        local failed = r(rc_reghdfe) | r(rc_creghdfe)
        foreach number of local numbers {
            local v_`number' = r(`number')
        }
        if "`mode'" == "default" local group_count = `v_groups'

        local status "recorded"
        if `failed' {
            * Either command failing is a difference under any settings.
            local status "error"
            test_fail "`testname' [`mode' settings]" "`differences'"
        }
        else if "`mode'" == "`asserted'" {
            local status = cond("`differences'" == "", "pass", "fail")
            if "`status'" == "pass" test_pass "`testname'"
            else test_fail "`testname'" "`differences'"
        }

        post hdfe_results ("$HDFE_SECTION") ("$HDFE_DATASET") ("`varlist'") ("`absorb'") ("`vce'") ("`mode'") ///
            (`nfull') (`v_n') (`v_singletons') (`v_groups') (`v_df_a_reghdfe') (`v_df_a_creghdfe') ///
            (`v_df_r_reghdfe') (`v_df_r_creghdfe') (`v_sf_b') (`v_sf_v') (`v_sf_stats') (`v_iterations') ///
            (`v_t_reghdfe') (`v_t_creghdfe') ("`status'") ("`differences'")
        di as text "    " %-9s "`status'" %-10s "`mode'" %-44s "absorb(`absorb') `vceopt'" ///
            " N=" %7.0f `v_n' " df_a=" %7.0f `v_df_a_reghdfe' " sf(b)=" %4.1f `v_sf_b' " sf(V)=" %4.1f `v_sf_v' ///
            "  reghdfe " %6.2f `v_t_reghdfe' "s  creghdfe " %6.2f `v_t_creghdfe' "s"
    }

    capture confirm variable __hdfe_group_r __hdfe_group_c
    if "`groups'" != "" & _rc == 0 {
        hdfe_same_partition __hdfe_group_r __hdfe_group_c
        if r(different) == 0 & r(groups) == `group_count' {
            test_pass "`testname' groupvar()"
        }
        else {
            test_fail "`testname' groupvar()" "`=r(different)' rows disagree; `=r(groups)' groups, expected `group_count'"
        }
    }
    capture drop __hdfe_group_r
    capture drop __hdfe_group_c
end

/*******************************************************************************
 * hdfe_common - The specifications every dataset supports
 ******************************************************************************/
capture program drop hdfe_common
program define hdfe_common
    global HDFE_SECTION "oneway"
    hdfe_case y x1 x2, absorb(id1)
    hdfe_case y x1 x2, absorb(id2)
    global HDFE_SECTION "twoway"
    hdfe_case y x1 x2, absorb(id1 id2) groups
    global HDFE_SECTION "vce"
    hdfe_case y x1 x2, absorb(id1 id2) vce(robust) converged
    hdfe_case y x1 x2, absorb(id1 id2) vce(cluster id1) converged
    hdfe_case y x1 x2, absorb(id1 id2) vce(cluster id2) converged
end

/*******************************************************************************
 * hdfe_oracle - Check both commands against an exact direct solution
 *
 * oracle() is a dummy-variable regression listing the same regressors first.
 * Its slopes, their covariance block, e(rss) and e(df_r) equal the absorbed
 * model's without iteration. reghdfe and creghdfe at default tolerance are
 * recorded against it; the asserted run is creghdfe, with tolerance() when
 * given, at $DEFAULT_SIGFIGS.
 *
 * Syntax: hdfe_oracle depvar indepvars, absorb(varlist) oracle(command) [tolerance(#)]
 ******************************************************************************/
capture program drop hdfe_oracle
program define hdfe_oracle
    syntax varlist(min=2 numeric), Absorb(varlist) oracle(string asis) [TOLerance(string)]

    gettoken depvar indepvars : varlist
    local k : word count `indepvars'
    local testname "$HDFE_DATASET: `varlist', absorb(`absorb') vs exact solution"

    tempname b_o V_o b V
    capture quietly `oracle'
    if _rc {
        test_fail "`testname'" "exact solution returned error `=_rc'"
        exit
    }
    matrix `b_o' = e(b)
    matrix `b_o' = `b_o'[1, 1..`k']
    matrix `V_o' = e(V)
    matrix `V_o' = `V_o'[1..`k', 1..`k']
    local o_rss = e(rss)
    local o_df_r = e(df_r)

    local runs "reghdfe creghdfe"
    if "`tolerance'" != "" local runs "`runs' creghdfe_tol"
    local asserted : word `: word count `runs'' of `runs'
    foreach run of local runs {
        local command = cond("`run'" == "reghdfe", "reghdfe", "creghdfe")
        local opt ""
        if "`run'" == "creghdfe_tol" local opt "tolerance(`tolerance')"
        local sf_b = .
        local sf_V = .
        local sf_rss = .
        local df_r = .
        capture quietly `command' `varlist', absorb(`absorb') `opt'
        local rc = _rc
        if `rc' == 0 {
            matrix `b' = e(b)
            matrix `b' = `b'[1, 1..`k']
            matrix `V' = e(V)
            matrix `V' = `V'[1..`k', 1..`k']
            local rss = e(rss)
            local df_r = e(df_r)
            matrix_min_sigfigs `b' `b_o'
            local sf_b = r(min_sigfigs)
            matrix_min_sigfigs `V' `V_o'
            local sf_V = r(min_sigfigs)
            sigfigs `rss' `o_rss'
            local sf_rss = r(sigfigs)
        }
        local status "recorded"
        if "`run'" == "`asserted'" {
            if `rc' == 0 & `df_r' == `o_df_r' & min(`sf_b', `sf_V', `sf_rss') >= $DEFAULT_SIGFIGS {
                local status "pass"
                test_pass "`testname'"
            }
            else {
                local status "fail"
                local sf_fmt : display %4.1f min(`sf_b', `sf_V', `sf_rss')
                test_fail "`testname'" "`command' `opt' rc=`rc' df_r:`df_r' vs `o_df_r' min sigfigs=`sf_fmt'"
            }
        }
        post hdfe_oracle ("$HDFE_DATASET") ("`varlist'") ("`absorb'") ("`command' `opt'") ///
            (`rc') (`df_r') (`o_df_r') (`sf_b') (`sf_V') (`sf_rss') ("`status'")
        di as text "    " %-9s "`status'" %-10s "exact" %-44s "`command' `opt'" ///
            " sf(b)=" %4.1f `sf_b' " sf(V)=" %4.1f `sf_V' " sf(rss)=" %4.1f `sf_rss'
    }
end

di as text "Running creghdfe vs reghdfe on the HDFE benchmark collection..."

* Provenance for the published tables
di "HDFE_ENV date=`c(current_date)'"
di "HDFE_ENV stata=`c(stata_version)' `c(edition_real)' (born `c(born_date)')"
di "HDFE_ENV processors=`c(processors)'"
foreach package in reghdfe ftools creghdfe {
    local version_line ""
    capture findfile `package'.ado
    if !_rc {
        tempname handle
        file open `handle' using `"`r(fn)'"', read text
        file read `handle' version_line
        file close `handle'
    }
    di `"HDFE_ENV `package'=`version_line'"'
}
_ctools_load
capture findfile `__ctools_plugin'
di `"HDFE_ENV plugin=`r(fn)'"'

capture postclose hdfe_results
postfile hdfe_results str12 section str28 dataset str40 varlist str32 absorb str16 vce str10 settings ///
    double(n_full n singletons groups df_a df_a_creghdfe df_r df_r_creghdfe sf_b sf_v sf_stats ///
    iterations t_reghdfe t_creghdfe) str8 status str200 differences using "`outdir'/results.dta", replace
capture postclose hdfe_oracle
postfile hdfe_oracle str28 dataset str40 varlist str32 absorb str40 run ///
    double(rc df_r df_r_exact sf_b sf_v sf_rss) str8 status using "`outdir'/oracle.dta", replace

/*******************************************************************************
 * Credit: borrower (id1) by month (id2), one connected component.
 * ruc is the raw 11-digit borrower number behind id1; t is the month string.
 ******************************************************************************/
foreach d in credit2 credit {
    hdfe_use `d'
    hdfe_common
    quietly gen int year = real(substr(t, 1, 4))
    quietly egen long id1_year = group(id1 year)
    global HDFE_SECTION "regroup"
    hdfe_case y x1 x2, absorb(ruc id2)
    hdfe_case y x1 x2, absorb(id1 year)
    hdfe_case y x1 x2, absorb(id1_year id2)
    global HDFE_SECTION "multiway"
    hdfe_case y x1 x2, absorb(id1 id2 year)
    if "`d'" == "credit2" {
        hdfe_oracle y x1 x2, absorb(id1 id2) oracle(areg y x1 x2 i.id2, absorb(id1))
    }
}

/*******************************************************************************
 * Soccer and the synthetic complete/uniform graphs: id1 and id2 only.
 ******************************************************************************/
hdfe_use soccer
hdfe_common
hdfe_oracle y x1 x2, absorb(id1 id2) oracle(areg y x1 x2 i.id2, absorb(id1))

foreach d in synthetic-complete synthetic-uniform-easy synthetic-uniform-hard synthetic-uniform-harder {
    hdfe_use `d'
    hdfe_common
}

/*******************************************************************************
 * Synthetic assortative matching: CEO (id1), firm (id2), year and two types.
 ******************************************************************************/
hdfe_use synthetic-assortative
hdfe_common
global HDFE_SECTION "regroup"
hdfe_case y x1 x2, absorb(ceo_type firm_type)
hdfe_case y x1 x2, absorb(id1 year)
hdfe_case y x1 x2, absorb(id2 year)
global HDFE_SECTION "multiway"
hdfe_case y x1 x2, absorb(id1 id2 year)
hdfe_case y x1 x2, absorb(id1 id2 year) vce(cluster id2) converged

/*******************************************************************************
 * Synthetic zigzag: a path graph, so the two-way projection equals the
 * projection on (id1, id2) cells and a one-way regression is exact. Two cells
 * hold two rows each: with x1 and x2 the two-way model is saturated.
 ******************************************************************************/
hdfe_use synthetic-zigzag
global HDFE_SECTION "oneway"
hdfe_case y x1 x2, absorb(id1)
hdfe_case y x1 x2, absorb(id2)
quietly egen long cell = group(id1 id2)
hdfe_oracle y x1, absorb(id1 id2) oracle(areg y x1, absorb(cell)) tolerance(1e-10)

* Documented difference: reghdfe posts estimates with zero residual degrees of
* freedom and a missing covariance; creghdfe returns r(498).
capture quietly reghdfe y x1 x2, absorb(id1 id2)
local saturated_rc = _rc
local saturated_df_r = e(df_r)
capture quietly creghdfe y x1 x2, absorb(id1 id2)
if `saturated_rc' == 0 & `saturated_df_r' == 0 & _rc == 498 {
    test_pass "synthetic-zigzag: saturated model returns r(498)"
}
else {
    test_fail "synthetic-zigzag: saturated model returns r(498)" ///
        "reghdfe rc=`saturated_rc' df_r=`saturated_df_r', creghdfe rc=`=_rc'"
}

/*******************************************************************************
 * Enron, patents, directors: id1 and id2 only, many components and singletons.
 ******************************************************************************/
foreach d in enron patents directors {
    hdfe_use `d'
    hdfe_common
}

/*******************************************************************************
 * GitHub: id1, id2, a third identifier and a commit timestamp string.
 ******************************************************************************/
hdfe_use github
hdfe_common
quietly gen int month = mofd(dofc(clock(t, "DMYhms")))
global HDFE_SECTION "regroup"
hdfe_case y x1 x2 len_commit_msg num_files, absorb(id1 id2)
hdfe_case y x1 x2, absorb(id3 month)
global HDFE_SECTION "multiway"
hdfe_case y x1 x2, absorb(id1 id2 id3)
hdfe_case y x1 x2, absorb(id1 id2 month)
hdfe_case y x1 x2, absorb(id1 id2 id3 month)

/*******************************************************************************
 * Workers: firm (id1), worker (id2), job title (id3), year (id4).
 ******************************************************************************/
hdfe_use workers
hdfe_common
quietly egen long firm_year = group(id1 id4)
global HDFE_SECTION "regroup"
hdfe_case y x1 x2, absorb(id1 id3)
hdfe_case y x1 x2, absorb(id2 firm_year)
global HDFE_SECTION "multiway"
hdfe_case y x1 x2, absorb(id1 id2 id4)
hdfe_case y x1 x2, absorb(id1 id2 id3)
hdfe_case y x1 x2, absorb(id1 id2 id3 id4)
hdfe_case y x1 x2, absorb(id1 id2 id3 id4) vce(cluster id1) converged

/*******************************************************************************
 * Schools: student (id1), teacher (id2), year (id3).
 ******************************************************************************/
hdfe_use schools
hdfe_common
quietly egen long teacher_year = group(id2 id3)
global HDFE_SECTION "regroup"
hdfe_case y x1 x2, absorb(id1 teacher_year)
global HDFE_SECTION "multiway"
hdfe_case y x1 x2, absorb(id1 id2 id3)
hdfe_case y x1 x2, absorb(id1 id2 id3) vce(cluster id2) converged

postclose hdfe_results
postclose hdfe_oracle
foreach table in results oracle {
    use "`outdir'/`table'.dta", clear
    quietly export delimited using "`outdir'/`table'.csv", replace
}

print_summary "hdfe_datasets"

* print_summary exits with an error when a check failed.
di "HDFE_DATASETS_COMPLETE"
