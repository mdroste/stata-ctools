*! version 1.3.1 20261124 github.com/mdroste/stata-ctools
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
