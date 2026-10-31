*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! Remove absorbed fixed effects from a list of variables
program define chdfe, eclass
    version 14.1
    capture syntax, Version
    if !_rc {
        ctools, version
        ereturn clear
        ereturn local version "`r(version)'"
        exit
    }
    * hdfe accepts verbose(#); a bare verbose also requests ctools timing.
    local 0 : subinstr local 0 "verbose(" "verbosity(", all
    local 0 : subinstr local 0 "v(" "verbosity(", all
    syntax varlist(numeric fv ts) [if] [in] [aw fw pw], Absorb(string) ///
        [GENerate(name) SAMPLE(name) CLEAR KEEPIDs KEEPVars(varlist) ///
         CLUSTERVars(string) Verbose VERBOSITY(integer 0) TIMEit ///
         TOLerance(real 1e-8) MAXITerations(real 10000) POOLsize(integer 10) ///
         ACCELeration(string) TRANSform(string) DOFadjustments(string) ///
         GROUPVar(name) KEEPSINgletons THReads(integer 0)]
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
    local transform = lower("`transform'")
    local trans = 0
    if "`transform'" != "" {
        if strpos("symmetric_kaczmarz","`transform'")==1 local trans 0
        else if strpos("kaczmarz","`transform'")==1 local trans 1
        else if strpos("cimmino","`transform'")==1 local trans 2
        else exit 198
    }
    local acceleration = lower("`acceleration'")
    local accel = 1
    if "`acceleration'" != "" {
        if inlist("`acceleration'","no","none","off") local accel 0
        else if "`acceleration'"=="cg" | strpos("conjugate_gradient","`acceleration'")==1 local accel 1
        else if "`acceleration'"=="sd" | strpos("steepest_descent","`acceleration'")==1 | strpos("steep_descent","`acceleration'")==1 local accel 2
        else if strpos("aitken","`acceleration'")==1 local accel 3
        else exit 198
    }
    if `trans'==1 & `accel'==1 {
        di as error "conjugate gradient requires a symmetric transform; use transform(sym) or acceleration(none)"
        exit 198
    }
    _ctools_hdfe_dof, `dofadjustments'
    local flags `s(flags)'
    local methods `s(methods)'
    if "`groupvar'" != "" & !mod(`flags',4) exit 198
    foreach v in `sample' `groupvar' {
        confirm new variable `v'
    }
    if "`sample'" != "" & "`sample'" == "`groupvar'" exit 198
    marksample touse, novarlist
    tempvar w used mobility
    local weightvar
    local weightkeep
    local wt = 0
    if "`weight'" != "" {
        quietly generate double `w' `exp' if `touse'
        quietly count if `touse' & `w'<0
        if r(N) exit 402
        if "`weight'"=="fweight" {
            quietly count if `touse' & `w'<. & `w'!=floor(`w')
            if r(N) exit 401
        }
        quietly replace `touse'=0 if missing(`w') | `w'==0
        local weightvar `w'
        local expression = strtrim(substr(`"`exp'"',2,.))
        capture confirm numeric variable `expression', exact
        if !_rc local weightkeep `expression'
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
            markout `touse' `categories' `continuous', strok
            local key : subinstr local categories " " "#", all
            local name `key'
            if "`continuous'"!="" {
                local op = cond(`intercept',"##","#")
                if `: word count `continuous''==1 local name `key'`op'c.`continuous'
                else local name `key'`op'c.(`continuous')
            }
            local original_names `original_names' `name'
            tempvar id
            quietly _ctools_group `categories' if `touse', gen(`id') strict
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
                * Plugin variables cannot repeat. Clone category and slope inputs.
                tempvar sid sx
                quietly clonevar `sid' = `id'
                if !`intercept' & `firstslope' local keptid `sid'
                local firstslope = 0
                quietly clonevar `sx' = `cv'
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
        markout `touse' `keys', strok
        local clusterkeep `clusterkeep' `keys'
        tempvar cid
        quietly _ctools_group `keys' if `touse', gen(`cid') strict
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
            tempvar x result
            quietly clonevar `x' = `v'
            quietly generate double `result' = .
            label variable `result' "Residuals: `label'"
            char `result'[name] `label'
            local inputs `inputs' `x'
            local outputs `outputs' `result'
            local destinations `destinations' `dest'
            local ++k
        }
    }
    if !`k' exit 102
    markout `touse' `inputs'
    quietly generate byte `used' = 0
    local groupoutput
    if "`groupvar'"!="" {
        quietly generate long `mobility' = .
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
    plugin call ctools_plugin `inputs' `fevars' `slopevars' `clusters' `weightvar' ///
        `outputs' `used' `groupoutput' `touse', ///
        "chdfe `k' `g' `nc' `wt' `keeps' `maxiterations' `flags' `poolsize' `trans' `accel' `hasgroup' `slopecodes' tol=`tolerance' `threadopt' `verbose' `storeids'"
    * No user variables change until all C work succeeds. Clear mode also
    * restores the original dataset if publication fails.
    if "`clear'"!="" preserve
    capture noisily {
        if "`clear'"!="" {
            local retain `outputs' `used' `groupoutput' `keepvars' `clusterkeep' `weightkeep' `cachevars'
            if "`keepids'"!="" local retain `retain' `idvars' `clusters'
            quietly keep if `used'
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
