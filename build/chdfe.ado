*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! Remove absorbed fixed effects from a list of variables
*! original: hdfe 3.2.9 (February 21, 2016; SSC distribution of August 6, 2016) by Sergio Correia
*! original: see help chdfe for details
program define chdfe, eclass
    version 14.1
    capture syntax, Version
    if !_rc {
        ctools, version
        ereturn clear
        ereturn local version "`r(version)'"
        exit
    }
    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    local __chdfe_args `"`0'"'
    * hdfe accepts verbose(#); a bare verbose also requests ctools timing.
    * verbose(#) arrives among the other options and is parsed on its own,
    * so no text of the command line is rewritten.
    syntax varlist(numeric fv ts) [if] [in] [aw fw pw], Absorb(string) ///
        [Generate(name) SAMPLE(name) CLEAR KEEPIDs KEEPVars(varlist) ///
         CLUSTERVars(string) Verbose VERBOSITY(integer -1) TIMEit ///
         TOLerance(real 1e-8) MAXITerations(real 10000) POOLsize(integer 10) ///
         ACCELeration(string) TRAnsform(string) DOFadjustments(string) ///
         GROUPVar(name) KEEPSINgletons THReads(integer 0) *]
    * verbosity(#) is the same setting
    local __verbosity = `verbosity'
    _chdfe_verbosity, `options'
    if `__verbosity' >= 0 local verbosity = `__verbosity'
    if ("`generate'" != "") + ("`clear'" != "") != 1 {
        di as error "specify exactly one of generate() and clear"
        exit 198
    }
    if "`clear'" != "" & "`sample'" != "" exit 198
    if "`generate'" != "" & "`keepids'`keepvars'" != "" exit 198
    if !inrange(`threads',0,256) | !inrange(`verbosity',0,5) | `poolsize'<1 | ///
        `tolerance'<1e-16 | `tolerance'>=1 | missing(`tolerance') exit 198
    if missing(`maxiterations') local maxiterations 2147483646
    if `maxiterations'<1 | `maxiterations'!=floor(`maxiterations') | `maxiterations'>2147483646 exit 198
    if `verbosity' | "`timeit'" != "" local verbose verbose
    * hdfe's transforms and accelerations, with its abbreviations; the
    * experimental random Kaczmarz and hybrid acceleration included. Either
    * option given runs that iteration on the rows, at hdfe's defaults too;
    * without them the plugin solves for the coefficients of the effects
    local rowspace = cond(`"`transform'`acceleration'"' != "", "rowspace", "")
    local transform = lower("`transform'")
    local trans = 0
    if "`transform'" != "" {
        if strpos("symmetric_kaczmarz","`transform'")==1 local trans 0
        else if strpos("kaczmarz","`transform'")==1 local trans 1
        else if strpos("cimmino","`transform'")==1 local trans 2
        else if strpos("rand_kaczmarz","`transform'")==1 | strpos("random_kaczmarz","`transform'")==1 local trans 3
        else exit 198
    }
    local acceleration = lower("`acceleration'")
    local accel = 1
    if "`acceleration'" != "" {
        if inlist("`acceleration'","no","none","off") local accel 0
        else if "`acceleration'"=="cg" | strpos("conjugate_gradient","`acceleration'")==1 local accel 1
        else if "`acceleration'"=="sd" | strpos("steepest_descent","`acceleration'")==1 | strpos("steep_descent","`acceleration'")==1 local accel 2
        else if strpos("aitken","`acceleration'")==1 local accel 3
        else if strpos("hybrid","`acceleration'")==1 local accel 4
        else exit 198
    }
    _ctools_hdfe_dof, `dofadjustments'
    local flags `s(flags)'
    local methods `s(methods)'
    if "`groupvar'" != "" & !mod(`flags',4) exit 198
    foreach v in `sample' `groupvar' {
        confirm new variable `v'
    }
    if "`sample'" != "" & "`sample'" == "`groupvar'" exit 198
    * chdfe ends by clearing e(). Released before marksample allocates its
    * marker, the old e(sample) storage is reused (as in creghdfe). Kept
    * when the qualifiers, weights or options can refer to e() or to
    * coefficients (_b[], _se[], _coef[] as a token); a failed run after
    * this point leaves e() empty.
    if !ustrregexm(`"`__chdfe_args'"', "(^|[^A-Za-z0-9_])(e[(]|_b\[|_se\[|_coef\[)") {
        ereturn clear
    }
    marksample touse, novarlist
    tempvar w used mobility
    local weightvar
    local weightkeep
    local wt = 0
    if "`weight'" != "" {
        * marksample has seen the weight: the sample leaves out zero and
        * missing weights, and Stata has refused negative weights (r(402))
        * and noninteger frequency weights (r(401)). A weight that is one
        * numeric variable goes to the plugin as it is.
        local expression = strtrim(substr(`"`exp'"',2,.))
        capture confirm numeric variable `expression', exact
        if !_rc & `: word count `expression'' == 1 {
            local weightkeep `expression'
            local weightvar `expression'
        }
        else {
            quietly generate double `w' `exp' if `touse'
            quietly count if `touse' & `w'<0
            if r(N) exit 402
            if "`weight'"=="fweight" {
                quietly count if `touse' & `w'<. & `w'!=floor(`w')
                if r(N) exit 401
            }
            quietly replace `touse'=0 if missing(`w') | `w'==0
            local weightvar `w'
        }
        local wt = cond("`weight'"=="fweight",2,cond("`weight'"=="pweight",3,1))
    }
    * Parse absorb() in Stata; numbering categories and all projections use C.
    local pending `absorb'
    local g = 0
    local original = 0
    local original_names
    local extended
    local fevars
    local slopevars
    local slopecodes
    local idvars
    while `"`pending'"' != "" {
        gettoken term pending : pending, bind
        if strpos(`"`term'"',"=") | strpos(`"`term'"',",") exit 198
        local items `term'
        if !strpos("`term'","#") & !strpos("`term'",".") unab items : `term'
        while `"`items'"' != "" {
            gettoken item items : items, bind
            local ++original
            local intercept = strpos("`item'","##")>0 | !strpos("`item'","c.")
            local parts : subinstr local item "#" " ", all
            local parts : subinstr local parts "(" " ", all
            local parts : subinstr local parts ")" " ", all
            local categories
            local continuous
            local incontinuous = 0
            foreach part of local parts {
                if substr("`part'",1,2)=="c." {
                    local incontinuous = 1
                    local part = substr("`part'",3,.)
                    if "`part'"=="" continue
                }
                if `incontinuous' {
                    unab cv : `part', max(1)
                    confirm numeric variable `cv'
                    local continuous `continuous' `cv'
                }
                else {
                    if substr("`part'",1,2)=="i." local part = substr("`part'",3,.)
                    unab cv : `part', max(1)
                    local categories `categories' `cv'
                }
            }
            if "`categories'"=="" exit 198
            if `: word count `continuous''>1 & !strpos("`item'","(") exit 198
            if strpos("`item'","##") & "`continuous'"=="" exit 198
            * Marked out lazily: before any statement that reads the marker
            * (_ctools_group, fvrevar of an operator term), else by the
            * plugin's screen of every variable it loads
            local __pend_fe `__pend_fe' `categories' `continuous'
            local key : subinstr local categories " " "#", all
            local name `key'
            if "`continuous'"!="" {
                local op = cond(`intercept',"##","#")
                if `: word count `continuous''==1 local name `key'`op'c.`continuous'
                else local name `key'`op'c.(`continuous')
            }
            local original_names `original_names' `name'
            * A numeric variable goes to the plugin as it is, which numbers
            * the levels of every effect on its sample. A string variable or
            * an interaction is numbered first; so is every effect with
            * keepids, where the plugin stores its identifiers in place.
            local id
            if "`keepids'" == "" & `: word count `categories'' == 1 {
                capture confirm numeric variable `categories', exact
                if !_rc local id `categories'
            }
            if "`id'" == "" {
                if "`__pend_fe'" != "" markout `touse' `__pend_fe', strok
                local __pend_fe
                tempvar id
                quietly _ctools_group `categories' if `touse', gen(`id') strict
            }
            if `intercept' {
                local ++g
                local fevars `fevars' `id'
                local slopecodes `slopecodes' 0
                local extended `extended' `key'
                local owner_`g' `original'
            }
            local keptid `id'
            local firstslope = 1
            foreach cv of local continuous {
                * With keepids the plugin stores an identifier in each
                * effect's variable, so each slope has its own copy.
                if "`keepids'" != "" {
                    tempvar sid sx
                    quietly clonevar `sid' = `id'
                    quietly clonevar `sx' = `cv'
                }
                else {
                    local sid `id'
                    local sx `cv'
                }
                if !`intercept' & `firstslope' local keptid `sid'
                local firstslope = 0
                local ++g
                local fevars `fevars' `sid'
                local slopevars `slopevars' `sx'
                local code = 1 + `intercept'
                local slopecodes `slopecodes' `code'
                local extended `extended' `key'#c.`cv'
                local owner_`g' `original'
            }
            local idvars `idvars' `keptid'
        }
    }
    if !`g' exit 198
    local clusters
    local clusterkeep
    local nc = 0
    foreach term of local clustervars {
        local keys : subinstr local term "#" " ", all
        local keys : subinstr local keys "i." "", all
        unab keys : `keys'
        local __pend_fe `__pend_fe' `keys'
        local clusterkeep `clusterkeep' `keys'
        local cid
        if "`keepids'" == "" & `: word count `keys'' == 1 {
            capture confirm numeric variable `keys'
            if !_rc local cid `keys'
        }
        if "`cid'" == "" {
            markout `touse' `__pend_fe', strok
            local __pend_fe
            tempvar cid
            quietly _ctools_group `keys' if `touse', gen(`cid') strict
        }
        local clusters `clusters' `cid'
        local ++nc
    }
    * Use fvrevar only for syntax expansion. Base and omitted indicators are
    * omitted from the output, as in hdfe.
    local inputs
    local outputs
    local destinations
    local cachevars
    local k = 0
    foreach term of local varlist {
        * an operator term's expansion depends on the marker
        capture confirm numeric variable `term', exact
        if _rc & "`__pend_fe'" != "" {
            markout `touse' `__pend_fe', strok
            local __pend_fe
        }
        fvrevar `term' if `touse'
        local expanded `r(varlist)'
        foreach v of local expanded {
            local fv : char `v'[fvrevar]
            local ts : char `v'[tsrevar]
            local name `v'
            local label `v'
            if "`fv'`ts'"!="" & substr("`v'",1,2)=="__" {
                if regexm("`fv'","[0-9]+o[.]") continue
                if regexm("`fv'","[0-9]+b[n]?[.]") & !strpos("`fv'","c.") continue
                local label `fv'`ts'
                local terms : subinstr local label "#" " ", all
                foreach expandedterm of local terms {
                    local base = regexr("`expandedterm'", "^.*[.]", "")
                    capture confirm variable `base', exact
                    if !_rc local cachevars `cachevars' `base'
                }
                if "`ts'"!="" {
                    quietly tsset
                    local cachevars `cachevars' `r(timevar)' `r(panelvar)'
                }
                local name : subinstr local label "." "__", all
                local name : subinstr local name "#" "_X_", all
                local name : permname __`name', length(30)
            }
            local dest `generate'`name'
            confirm name `dest'
            if "`clear'"=="" confirm new variable `dest'
            if "`dest'"=="`sample'" | "`dest'"=="`groupvar'" exit 198
            if `: list dest in destinations' exit 198
            * The plugin only reads its inputs, so the variable itself is
            * passed; the outputs are created together below
            tempvar result
            local inputs `inputs' `v'
            local outputs `outputs' `result'
            local destinations `destinations' `dest'
            local ++k
            local outlabel_`k' `label'
        }
    }
    if !`k' exit 102
    * The rows missing in an input, or in an effect, slope or cluster
    * variable not yet marked out, are dropped by the plugin (screen), which
    * loads all of them; that is the sample these markouts would leave
    local __screen
    if "`__pend_fe'`inputs'" != "" local __screen screen
    * The outputs are added without filling them: the plugin writes every
    * observation (fillall), missing outside the sample
    mata: (void) st_addvar("double", tokens(st_local("outputs")), 1)
    forvalues j = 1/`k' {
        local result : word `j' of `outputs'
        label variable `result' "Residuals: `outlabel_`j''"
        char `result'[name] `outlabel_`j''
    }
    * Without sample() the sample marker itself returns the rows used: the
    * plugin sets it to 0 on the rows it drops, one store per dropped row.
    * sample(newvar) keeps its own variable, created here so that it follows
    * the outputs in the dataset.
    local prefilled
    if "`sample'" == "" {
        local used `touse'
        local prefilled prefilled
    }
    else mata: (void) st_addvar("byte", st_local("used"), 1)
    local groupoutput
    if "`groupvar'"!="" {
        mata: (void) st_addvar("long", st_local("mobility"), 1)
        local groupoutput `mobility'
    }
    local hasgroup = "`groupvar'"!=""
    local keeps = "`keepsingletons'"!=""
    local threadopt
    if `threads' local threadopt threads(`threads')
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local storeids
    if "`keepids'"!="" local storeids storeids
    capture noisily plugin call ctools_plugin `inputs' `fevars' `slopevars' `clusters' `weightvar' ///
        `outputs' `used' `groupoutput' `touse', ///
        "chdfe `k' `g' `nc' `wt' `keeps' `maxiterations' `flags' `poolsize' `trans' `accel' `hasgroup' `slopecodes' tol=`tolerance' `threadopt' `verbose' `storeids' `prefilled' `rowspace' fillall `__screen'"
    local rc = _rc
    * a plugin that did not update the marker would report dropped rows as
    * used, and one that did not fill every row would leave the unfilled
    * outputs undefined outside the sample
    if !`rc' & (("`prefilled'" != "" & "`chdfe_marked'" != "1") | "`chdfe_filled'" != "1") {
        di as error "chdfe: rebuild the ctools plugin alongside the ado files"
        local rc = 601
    }
    * hdfe's messages: the singletons dropped, its warning on Kaczmarz with
    * conjugate gradient (which it runs anyway), and the iterations
    if `rc' {
        if `trans' == 1 & `accel' == 1 di as error "(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)"
        exit `rc'
    }
    if `chdfe_singletons' > 0 di as text "(dropped `chdfe_singletons' singleton observations)"
    if `trans' == 1 & `accel' == 1 di as error "(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)"
    di as text "(converged in `chdfe_iterations' iterations)"
    * No user variables change until all C work succeeds. Clear mode also
    * restores the original dataset if publication fails.
    if "`clear'"!="" preserve
    capture noisily {
        if "`clear'"!="" {
            * As hdfe: the cluster variables are not kept (only as
            * __CL#__ with keepids)
            local retain `outputs' `used' `groupoutput' `keepvars' `weightkeep' `cachevars'
            if "`keepids'"!="" local retain `retain' `idvars' `clusters'
            * every observation in the sample: nothing to drop (a pass)
            if `chdfe_N' < _N quietly keep if `used'
            keep `retain'
            foreach dest of local destinations {
                capture drop `dest'
            }
        }
        forvalues j=1/`k' {
            local source : word `j' of `outputs'
            local dest : word `j' of `destinations'
            rename `source' `dest'
        }
        if "`sample'"!="" {
            rename `used' `sample'
            label variable `sample' "[HDFE Sample]"
        }
        if "`groupvar'"!="" & `chdfe_group_saved' rename `mobility' `groupvar'
        if "`clear'`keepids'"=="clearkeepids" {
            forvalues j=1/`original' {
                local source : word `j' of `idvars'
                capture drop __ID`j'__
                rename `source' __ID`j'__
            }
            forvalues j=1/`nc' {
                local source : word `j' of `clusters'
                capture drop __CL`j'__
                rename `source' __CL`j'__
            }
        }
    }
    local rc = _rc
    if "`clear'"!="" {
        if `rc' restore
        else restore, not
    }
    if `rc' {
        if "`clear'"=="" {
            foreach v in `destinations' `sample' `groupvar' {
                capture drop `v'
            }
        }
        exit `rc'
    }
    ereturn clear
    ereturn scalar N = `chdfe_N'
    ereturn scalar df_a = `chdfe_df'
    ereturn scalar N_hdfe = `original'
    ereturn scalar N_hdfe_extended = `g'
    ereturn scalar mobility = `chdfe_mobility'
    ereturn scalar M_due_to_nested = `chdfe_nested'
    ereturn scalar num_singletons = `chdfe_singletons'
    ereturn scalar iteration = `chdfe_iterations'
    ereturn scalar converged = 1
    forvalues j=1/`g' {
        ereturn scalar K`j' = `chdfe_K`j''
        ereturn scalar M`j' = `chdfe_M`j''
        ereturn scalar M`j'_exact = `chdfe_exact`j''
        ereturn scalar M`j'_nested = `chdfe_nested`j''
        ereturn scalar G`j' = `owner_`j''
    }
    ereturn local absvars `original_names'
    ereturn local extended_absvars `extended'
    ereturn local dofmethod `methods'
    ereturn local cmd chdfe
end

* verbose(#), the only option left after chdfe's own: the caller's local
* verbosity (0 without it); anything else is an error
program define _chdfe_verbosity
    version 14.1
    syntax [, Verbose(integer 0)]
    c_local verbosity `verbose'
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* _ctools_group — group identifiers from the ctools plugin, in place of egen group()
* 
* Usage:  _ctools_group varlist [if] [in], generate(newvar) [label]
* _ctools_group varlist [if] [in], count [strict]
* strict returns plugin errors instead of sorting in Stata.
* 
* Numbers the distinct values of varlist 1, 2, ... in sort order, as
* egen group() does: observations excluded by if/in, or with a missing value
* in any variable of varlist (string "" counts as missing), get missing. The
* new variable is a long labeled "group(varlist)". label also defines a value
* label named after the new variable, with egen group(, label)'s text for
* each group. count only counts the groups. Returns r(N_groups), and with
* label r(reps), one observation number of each group in group order.
* 
* The plugin reads, sorts and numbers the keys. Stata numbers them by sorting
* only when the plugin is unavailable or cannot read a key (strL text over
* 2,045 bytes).

program define _ctools_group, rclass
    version 14.1
    syntax varlist [if] [in], [GENerate(name) LABel COUNT STRICT]
    if ("`generate'" == "") == ("`count'" == "") | ("`label'" != "" & "`generate'" == "") {
        di as error "_ctools_group: specify generate() [label] or count"
        exit 198
    }
    if "`generate'" != "" {
        confirm new variable `generate'
        if "`label'" != "" {
            capture label list `generate'
            if !_rc {
                di as error "default label `generate' already defined"
                exit 110
            }
        }
    }
    local nkeys : word count `varlist'

    local rc = 199
    capture _ctools_load
    if !_rc {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if _rc == 0 | _rc == 110 local rc = 0
    }
    capture confirm number 0
    if !`rc' {
        _ctools_strw `varlist'
        local out ""
        if "`generate'" != "" {
            * The plugin writes every observation of the range, so without in
            * the variable is created without initial values.
            if `"`in'"' == "" mata: (void) st_addvar("long", st_local("generate"), 1)
            else quietly generate long `generate' = .
            local out "out"
        }
        local reps = cond("`label'" != "", "reps", "")
        * Stata answers an if one observation at a time. Without one the
        * plugin skips the check; an if that names a numeric variable (the
        * callers' touse) passes that variable last, as a sample marker the
        * plugin reads in parallel. Any other if is marked first, unless in
        * limits the range (marking covers every observation); that marker is
        * the last variable until it is dropped, which is then free.
        local select "skipif"
        local marker ""
        local marked ""
        local ifcond ""
        if `"`if'"' != "" {
            local select "touse"
            local marker = strtrim(substr(`"`if'"', 3, .))
            capture confirm numeric variable `marker', exact
            if _rc | wordcount(`"`marker'"') != 1 {
                local marker ""
                if `"`in'"' == "" {
                    tempvar marker
                    capture mark `marker' `if'
                    local rc = _rc
                    local marked `marker'
                }
                else {
                    local select ""
                    local ifcond `"`if'"'
                }
            }
        }
        local ctools_group_reps ""
        if !`rc' {
            capture plugin call ctools_plugin `varlist' `generate' `marker' `ifcond' `in', ///
                "group nkeys=`nkeys' `out' `reps' `select'"
            local rc = _rc
        }
        if "`marked'" != "" capture drop `marked'
        * Every failure, Break included, drops the unfilled output variable.
        if `rc' & "`generate'" != "" capture drop `generate'
        if `rc' == 1 exit 1
        if `rc' == 134 {
            di as error "too many values for {bf:label define}"
            exit 134
        }
        if !`rc' {
            local n_groups = scalar(_ctools_group_n)
            local group_reps `"`ctools_group_reps'"'
        }
        capture scalar drop _ctools_group_n
    }
    if `rc' & "`strict'" != "" exit `rc'
    if `rc' {
        _ctools_group_sort `varlist' `if' `in', generate(`generate') `count' `label'
        local n_groups = r(N_groups)
        local group_reps `"`r(reps)'"'
    }

    if "`generate'" != "" {
        if udstrlen("group(`varlist')") > 80 {
            note `generate' : group(`varlist')
            label variable `generate' "see notes"
        }
        else {
            label variable `generate' "group(`varlist')"
        }
        if "`label'" != "" {
            if `n_groups' > 0 {
                mata: _ctools_group_label(st_local("varlist"), st_local("group_reps"), ///
                    st_local("generate"))
            }
            capture label list `generate'
            if !_rc label values `generate' `generate'
        }
    }
    return scalar N_groups = `n_groups'
    * label: the observation number of one observation of each group
    if "`label'" != "" return local reps `"`group_reps'"'
end

* Number groups by sorting in Stata (the data and their sort order are restored).
program define _ctools_group_sort, rclass sortpreserve
    syntax varlist [if] [in], [GENerate(name) COUNT LABel]
    marksample touse, strok
    tempvar order g
    quietly generate long `order' = _n
    sort `touse' `varlist' `order'
    quietly by `touse' `varlist': generate long `g' = (_n == 1) if `touse'
    quietly replace `g' = sum(`g') if `touse'
    quietly summarize `g', meanonly
    return scalar N_groups = cond(r(N), r(max), 0)
    if "`label'" != "" {
        mata: _ctools_group_firsts(st_local("g"), st_local("order"))
        return local reps `"`__reps'"'
    }
    if "`generate'" != "" rename `g' `generate'
end

mata:
/* Original observation number of the first observation of each group. */
void _ctools_group_firsts(string scalar g, string scalar order)
{
    real colvector id, obs, previous, first

    id = st_data(., g)
    obs = st_data(., order)
    st_local("__reps", "")
    if (rows(id) == 0) return
    previous = (rows(id) > 1 ? (. \ id[|1 \ rows(id) - 1|]) : .)
    first = selectindex(id :< . :& (id :!= previous))
    if (rows(first)) st_local("__reps", invtokens(strofreal(obs[first], "%12.0g")'))
}

/*
    Value label of egen group(varlist), label for groups 1..G, from one
    observation of each group (reps, in group order). Each variable
    contributes, separated by spaces: a string's value; a numeric value's
    label text when its value label defines one, otherwise the value
    ("%18.0g" for integers, "%9.0g" or "%9,0g" otherwise). Text is cut to 32
    characters and trimmed, and a string's NUL bytes become \0, which is what
    egen group() does when called under version 16 or earlier, as the ctools
    ado files are.
*/
void _ctools_group_label(string scalar varlist, string scalar reps, string scalar lname)
{
    real rowvector v
    real colvector obs, s
    real scalar i, j
    string scalar fmt, vl
    string colvector x, y

    v = st_varindex(tokens(varlist))
    obs = strtoreal(tokens(reps))'
    if (rows(obs) > 65536) {
        errprintf("too many values for {bf:label define}\n")
        exit(134)
    }
    fmt = (st_global("c(dp)") == "comma" ? "%9,0g" : "%9.0g")
    for (j = 1; j <= cols(v); j++) {
        if (st_isnumvar(v[j])) {
            s = st_data(obs, v[j])
            x = J(rows(s), 1, "")
            vl = st_varvaluelabel(v[j])
            if (vl != "") {
                if (st_vlexists(vl)) x = ustrtrim(usubstr(st_vlmap(vl, s), 1, 32))
            }
            for (i = 1; i <= rows(s); i++) {
                if (x[i] == "") x[i] = strofreal(s[i], (s[i] != trunc(s[i]) ? fmt : "%18.0g"))
            }
        }
        else {
            x = st_sdata(obs, v[j])
            x = usubinstr(x, char(0), "\0", .)
            x = ustrtrim(usubstr(x, 1, 32))
        }
        y = (j == 1 ? x : y :+ " " :+ x)
    }
    st_vlmodify(lname, 1::rows(y), y)
}
end

* dofadjustments() of the HDFE estimators.
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

* Shared plugin loader and identity contract
program define _ctools_load
    version 14.1
    * Not rclass, and findfile's r() is undone, so the caller's r() survives;
    * identity comes back through c_local.
    local __os = c(os)
    local __machine = c(machine_type)
    local __is_mac = 0
    if "`__os'" == "MacOSX" {
        local __is_mac = 1
    }
    else if strpos(lower("`__machine'"), "mac") > 0 {
        local __is_mac = 1
    }

    local __plugin = ""
    if "`__os'" == "Windows" {
        local __plugin "ctools_windows.plugin"
    }
    else if `__is_mac' {
        local __is_arm = 0
        if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
            local __is_arm = 1
        }
        if `__is_arm' == 0 & !strpos(lower("`__machine'"), "intel") & !strpos(lower("`__machine'"), "x86") {
            tempfile __archfile
            quietly shell uname -m > "`__archfile'" 2>&1
            tempname __fh
            file open `__fh' using "`__archfile'", read text
            file read `__fh' __archline
            file close `__fh'
            capture erase "`__archfile'"
            if strpos("`__archline'", "arm64") > 0 {
                local __is_arm = 1
            }
        }
        if `__is_arm' {
            local __plugin "ctools_mac_arm.plugin"
        }
        else {
            local __plugin "ctools_mac_x86.plugin"
        }
    }
    else if "`__os'" == "Unix" {
        local __plugin "ctools_linux.plugin"
    }
    else {
        local __plugin "ctools.plugin"
    }

    * The ado/plugin argument protocol, stamped from config/release.json.
    local __ado_protocol "4"
    local __reinstall `"net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace"'

    tempname __rhold
    _return hold `__rhold'
    capture findfile `__plugin'
    local __found = (_rc == 0)
    local __path `"`r(fn)'"'
    _return restore `__rhold'
    local __wanted "`__plugin'"
    if !`__found' & "`__plugin'" != "ctools.plugin" local __plugin "ctools.plugin"
    capture program ctools_plugin, plugin using("`__plugin'")
    local __load_rc = _rc
    if `__load_rc' != 0 & `__load_rc' != 110 {
        di as error "ctools: could not load the plugin `__plugin'"
        di as error `"  Stata reports c(os) "`__os'" and c(machine_type) "`__machine'"; ctools selected `__wanted'."'
        if `__found' di as error `"  The plugin file is `__path'."'
        else di as error "  `__wanted' was not found on the ado-path; reinstall ctools."
        di as error "  Stata's plugin loader reported:"
        capture noisily program ctools_plugin, plugin using("`__plugin'")
        if "`__os'" == "Unix" & !`__is_mac' {
            di as error "  The distribution Linux plugin needs x86-64 and glibc 2.17 or later; OpenMP is linked statically."
            di as error "  Run ldd -r on the plugin to inspect missing libraries or symbols."
        }
        else if `__is_mac' {
            di as error "  On macOS, check that the plugin is the notarized ctools release and is not quarantined."
        }
        di as error "  See https://github.com/mdroste/stata-ctools/blob/main/docs/INSTALLATION.md#troubleshooting"
        exit `__load_rc'
    }
    local ctools_plugin_hdfe_screen
    local ctools_plugin_hdfe_markout
    local ctools_plugin_hdfe_fvlevels
    capture plugin call ctools_plugin, "version"
    local __id_rc = _rc
    if `__id_rc' | "`ctools_plugin_protocol'" != "`__ado_protocol'" {
        local __has = cond("`ctools_plugin_protocol'" == "", "none", "`ctools_plugin_protocol'")
        di as error "ctools: the plugin `__plugin' does not match the installed ado files."
        local __ident = cond("`ctools_plugin_version'" == "", "", " (plugin version `ctools_plugin_version', revision `ctools_plugin_revision')")
        di as error "  The ado files need plugin protocol `__ado_protocol'; the loaded plugin reports `__has'`__ident'."
        di as error "  Reinstall ctools so its ado files and plugins come from the same release:"
        di as error `"    `__reinstall'"'
        di as error "  Then restart Stata, because a plugin stays loaded until Stata exits."
        exit 9
    }
    c_local __ctools_plugin "`__plugin'"
    c_local __ctools_plugin_version "`ctools_plugin_version'"
    c_local __ctools_plugin_revision "`ctools_plugin_revision'"
    c_local __ctools_plugin_protocol "`ctools_plugin_protocol'"
    c_local __ctools_plugin_hdfe_screen = ("`ctools_plugin_hdfe_screen'" == "1")
    c_local __ctools_plugin_hdfe_markout = ("`ctools_plugin_hdfe_markout'" == "1")
    c_local __ctools_plugin_hdfe_fvlevels = ("`ctools_plugin_hdfe_fvlevels'" == "1")
    capture confirm number 0
end

* _ctools_strw — set string width metadata for ctools flat buffer optimization
* 
* Usage:  _ctools_strw varlist
* 
* Sets local __ctools_strw in CALLER's scope containing comma-separated
* string variable widths (0 = numeric, >0 = str width e.g. 17 for str17).
* The ctools C plugin reads this local automatically via SF_macro_use
* to enable flat buffer string I/O (eliminates per-string strlen, memcpy,
* and arena atomic CAS overhead).
* 
* Example:
* _ctools_strw `allvars'
* plugin call ctools_plugin `allvars', "mycommand args"
* * C plugin automatically detects __ctools_strw local

program define _ctools_strw
    syntax varlist

    local __strw ""
    local __sep ""
    foreach v of local varlist {
        local __type : type `v'
        if substr("`__type'", 1, 3) == "str" {
            local __w = real(substr("`__type'", 4, .))
            if missing(`__w') local __w = 0
        }
        else {
            local __w = 0
        }
        local __strw "`__strw'`__sep'`__w'"
        local __sep ","
    }
    c_local __ctools_strw `__strw'
end
* ctools-helpers:end
