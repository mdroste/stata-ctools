*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! Bootstrap options and Stata-RNG draws of the HDFE estimators.
* Design: docs/development/BOOTSTRAP_PLAN.md (sections 4.2-4.3, 4.7 and 5.8).
*
* _ctools_hdfe_boot parse, spec(string) [fe]
*     spec is the text of vce() or fevce(): bootstrap[, bootopts]. Returns
*     r(kind) = bootstrap, or empty when spec is not a bootstrap, and
*     r(reps) r(seed) r(cluster) r(strata) r(mse) r(level) r(dots)
*     r(saving) r(draws) r(bsdraws) r(designopts) (1 when any of reps(),
*     seed(), cluster(), strata() or draws() was given), r(given) (which of
*     level, mse and dots, for dots() or nodots, were given). With fe, also the
*     fixed-effect options r(sestub) r(ci) r(nreps) r(dse); without it they
*     are errors. saving(filename [, replace double every(#)]) is checked as
*     bootstrap checks it and returned split: r(save_file) r(save_double)
*     r(save_every) r(save_replace). The hidden _bsdraws(varlist) takes one
*     multiplicity variable per replicate in place of the draws; r(draws)
*     is then supplied.
*
* _ctools_hdfe_boot run, cmd(name) reps(#) draws(fast|stata|supplied)
*     dots(#) kb(#) [seed(string) threads(string) targets(string)
*     bsdraws(varlist) savefile(string) every(#) strataobs]
*     The replicates of the fit cmd's plugin retained, in chunks of about
*     50 (dots and Break between chunks; a multiple of the replicates the
*     plugin runs at once, scalar __cmd_bs_par, so that no thread idles at
*     a chunk's end): plugin call "cmd threads boot_chunk ... targets".
*     draws(stata) draws bsample's uniforms with Stata's generator for the
*     scalar __cmd_bs_U units and the plugin maps them as bsample does;
*     draws(fast) takes a master seed from Stata's generator
*     (docs/development/BOOTSTRAP_PLAN.md, 5.8); draws(supplied) passes the
*     bsdraws() variables. With savefile(), the kb coefficients of each
*     replicate are posted to that file (variables _r1 ... _rkb, double),
*     for savefinish. Returns r(rngstate), the random-number state before
*     the first draw.
*
* _ctools_hdfe_boot savefinish, tmp(file) file(filename) map(matrix)
*     [double replace command() seed() n(#) cluster() nclust(#) strata()
*     nstrata(#) cmdname(name)]
*     After e(b) is posted: the replicates in tmp written to file as
*     bootstrap's saving() writes them (one variable per e(b) term, named
*     and labeled as bootstrap names and labels them, and its
*     characteristics, so that bstat reads the file). map is 1 x colsof(e(b)):
*     the replicate column of each term, or 0 for a term fixed at zero
*     (omitted and base levels). Floats unless double.
*
* _ctools_hdfe_boot level, cmd(name) bslevel(#) bsgiven(0|1) [diopts(string)
*     cmdlevel(string)]
*     The level of vce(bootstrap)'s results when the command's own level()
*     can also be given (in diopts(), the display options, or as
*     cmdlevel()), as Stata's bootstrap treats them: the command's level()
*     sets the bootstrap's unless vce(bootstrap, level()) gave another, an
*     error (r(198)). Returns r(level) (bslevel when the command gave none)
*     and r(diopts), diopts with level(r(level)), so that the table at
*     estimation shows the bootstrap's intervals (a replay shows c(level)).
*     cmdlevel() is checked as level() is checked.
*
* _ctools_hdfe_boot seed [, seed(string) fast]
*     Runs set seed when seed() is given and returns r(rngstate), the
*     random-number state before any draw (posted as e(rngstate)). With
*     fast, also takes the draws(fast) master seed from Stata's generator
*     as two exact 31-bit halves, r(seedhi) and r(seedlo), as cbsample does.
*
* _ctools_hdfe_boot uniforms, units(#) reps(#) generate(stub) [strataobs(cmd)]
*     For draws(stata): the uniforms bsample takes, in its order, for reps
*     successive replicates of a frame of units units: per replicate, units
*     uniforms for the draws, then units for the tie-breaks. With
*     strataobs (observations resampled within strata), bsample as the
*     bootstrap prefix calls it (without weight()) then takes units more
*     that do not change the sample when expand added rows, that is when a
*     stratum drew one of its units twice: its by strata: check then fails
*     on the expanded data, which expand left unsorted, and it redraws an
*     order (StrSRSWR). They are taken and discarded. The strata that can
*     draw a unit twice are those the plugin of cmd saved at retain
*     (scalars __cmd_bs_dupK, __cmd_bs_dupo#, __cmd_bs_dupn#; dupK missing:
*     so many that a replicate without a repeated draw has probability
*     below 2^-64, and the extras are always taken). Packed in that
*     order into double variables stub1, stub2, ... (value q, 0-based, in
*     variable floor(q / _N) + 1, observation mod(q, _N) + 1), as few as hold
*     the 2 * reps * units values; variables left by an earlier call are
*     reused. Returns r(varlist), the variables of this call's values. The
*     plugin applies bsample's mapping to them.
*
* _ctools_hdfe_boot wildparse, spec(string)
*     spec is the text of wildbootstrap() (its options, or empty). Returns
*     r(reps) (default 999) r(seed) r(weighttype) r(null) r(coefs) r(ptype)
*     r(ptype_code) (0 symmetric, 1 equaltail, 2 lower, 3 upper)
*     r(enumerate) r(fedgp) (full or partialled) r(ci) r(noci) r(level)
*     r(bootcluster) r(dots) r(projected) (hidden _projected: the projected
*     path even for nested effects, for validation).
*
* _ctools_hdfe_boot wilddraw, clusters(#) reps(#) type(name) [noenumerate]
*     The wild bootstrap's weights for clusters bootstrap clusters (in
*     sorted order) and reps replicates, drawn from Stata's generator
*     exactly as boottest draws them: one clusters x (reps + 1) matrix (two
*     for webb) whose first column is discarded, so that the same seed gives
*     boottest's replicates. type is rademacher, mammen, webb, normal or
*     gamma; the weights have mean 0 and variance 1. With rademacher and
*     2^clusters <= reps (boottest's rule), all 2^clusters sign vectors are
*     used instead and nothing is drawn (unless noenumerate); the plugin
*     generates them. Returns r(enumerate) and r(reps) (2^clusters when
*     enumerating). The draws are kept in Mata (_ctools_wild_v) until
*     wildclear.
*
* _ctools_hdfe_boot wildweights, first(#) count(#) generate(stub)
*     Replicates first .. first + count - 1 (1-based) of the last wilddraw,
*     replicate after replicate (clusters values each), packed into double
*     variables stub1, stub2, ... as uniforms packs its values, and returned
*     in r(varlist).
*
* _ctools_hdfe_boot wildclear
*     Drops the draws.

program define _ctools_hdfe_boot, rclass
    version 14.1
    gettoken sub 0 : 0, parse(" ,")
    if "`sub'" == "parse" Parse `0'
    else if "`sub'" == "seed" Seed `0'
    else if "`sub'" == "uniforms" Uniforms `0'
    else if "`sub'" == "wildparse" WildParse `0'
    else if "`sub'" == "wilddraw" WildDraw `0'
    else if "`sub'" == "wildweights" WildWeights `0'
    else if "`sub'" == "wildclear" WildClear `0'
    else if "`sub'" == "run" Run `0'
    else if "`sub'" == "savefinish" SaveFinish `0'
    else if "`sub'" == "level" BootLevel `0'
    else {
        di as error "_ctools_hdfe_boot: unknown subcommand `sub'"
        exit 198
    }
    return add
end

program define Parse, rclass
    syntax, [SPEC(string) FE]
    gettoken word rest : spec, parse(" ,")
    _ctools_vce kind unused : `word'
    if "`kind'" != "bootstrap" {
        return local kind `kind'
        exit
    }
    local rest = strtrim(`"`rest'"')
    if `"`rest'"' != "" {
        gettoken comma rest : rest, parse(",")
        if "`comma'" != "," {
            di as error `"bootstrap options must follow a comma: `spec'"'
            exit 198
        }
    }
    local 0 `", `rest'"'
    local feopts
    if "`fe'" != "" local feopts SEstub(name) CI NREPS DSE(name)
    syntax [, Reps(integer -1) SEED(string) CLuster(string) STRata(varlist) MSE ///
        Level(string) NODOTS DOTS(integer -1) SAVing(string asis) DRAWs(string) ///
        _BSDRAWS(varlist numeric) `feopts']
    * level() checked as cilevel checks it, and remembered as given
    local given
    if `"`level'"' != "" {
        local given level
        local 0 `", level(`level')"'
        syntax [, Level(cilevel)]
    }
    else local level = c(level)
    if "`mse'" != "" local given `given' mse
    if "`nodots'" != "" | `dots' != -1 local given `given' dots

    local designopts = (`reps' != -1 | `"`seed'"' != "" | `"`cluster'"' != "" | ///
        "`strata'" != "" | "`draws'" != "")
    if `reps' == -1 local reps 50
    if `reps' < 2 {
        di as error "bootstrap option reps() requires an integer greater than 1"
        exit 198
    }
    if `"`seed'"' != "" {
        * Validate as set seed would, without changing the state yet
        capture confirm integer number `seed'
        if _rc & substr(`"`seed'"', 1, 1) != "X" {
            di as error `"bootstrap option seed(`seed') invalid"'
            exit 198
        }
    }
    if `: word count `cluster'' > 1 {
        di as error "bootstrap option cluster() takes one term: one resampling dimension"
        exit 198
    }
    local draws = lower(strtrim("`draws'"))
    if "`_bsdraws'" != "" {
        if "`draws'" != "" {
            di as error "bootstrap options draws() and _bsdraws() may not be combined"
            exit 198
        }
        if `: word count `_bsdraws'' != `reps' {
            di as error "bootstrap option _bsdraws() takes one variable per replicate (reps(`reps'))"
            exit 198
        }
        local draws supplied
    }
    * draws(stata), bsample's draws, by default (plan 14, question 2): the
    * results reproduce bootstrap ...: cmd for the same seed
    if "`draws'" == "" local draws stata
    if !inlist("`draws'", "fast", "stata", "supplied") {
        di as error "bootstrap option draws() must be fast or stata"
        exit 198
    }
    if `"`saving'"' != "" {
        _prefix_saving `saving'
        return local save_file `"`s(filename)'"'
        return local save_double `s(double)'
        return local save_every = cond(`"`s(every)'"' != "", substr(`"`s(every)'"', 7, strlen(`"`s(every)'"') - 7), "")
        return local save_replace `s(replace)'
    }
    if "`nodots'" != "" & `dots' != -1 {
        di as error "bootstrap options nodots and dots() may not be combined"
        exit 198
    }
    if `dots' == -1 local dots = cond("`nodots'" != "", 0, 1)
    if `dots' < 0 {
        di as error "bootstrap option dots() requires a nonnegative integer"
        exit 198
    }

    return local kind bootstrap
    return scalar reps = `reps'
    return scalar designopts = `designopts'
    return local seed `"`seed'"'
    return local cluster `"`cluster'"'
    return local strata `"`strata'"'
    return scalar mse = ("`mse'" != "")
    return scalar level = `level'
    return scalar dots = `dots'
    return local given `given'
    return local saving `"`saving'"'
    return local draws `draws'
    return local bsdraws `_bsdraws'
    if "`fe'" != "" {
        return local sestub `sestub'
        return scalar ci = ("`ci'" != "")
        return scalar nreps = ("`nreps'" != "")
        return local dse `dse'
    }
end

program define Run, rclass
    syntax, CMD(name) REPS(integer) DRAWS(name) DOTS(integer) KB(integer) ///
        [SEED(string) THREADS(string) TARGETS(string) BSDRAWS(varlist numeric) ///
        SAVEFILE(string) EVERY(string) STRATAOBS]
    local fast = cond("`draws'" == "fast", "fast", "")
    Seed, seed(`"`seed'"') `fast'
    local rngstate `"`r(rngstate)'"'
    local how
    if "`draws'" == "fast" local how " seedhi=`r(seedhi)' seedlo=`r(seedlo)'"
    else if "`draws'" == "stata" local how " stata"
    else if "`draws'" == "supplied" {
        local how " supplied"
        if `: word count `bsdraws'' != `reps' {
            di as error "`cmd': _bsdraws() takes one variable per replicate"
            exit 198
        }
    }
    local U = scalar(__`cmd'_bs_U)
    local par = scalar(__`cmd'_bs_par)
    if missing(`par') | `par' < 1 local par 1
    local chunk = min(`par' * ceil(50 / `par'), 400)
    if "`draws'" == "stata" {
        * At most about 1000 variables of uniforms per chunk (2 U values per
        * replicate; U exceeds _N when fweighted rows stand for more
        * observations), and about 2^26 values (512 MB) unless one replicate
        * needs more
        local vmax = min(1000, c(maxvar) - c(k) - 10)
        local nmax = floor(`vmax' * _N / (2 * `U'))
        if `nmax' >= 1 local nmax = min(`nmax', max(1, floor(2^26 / (2 * `U'))))
        if `nmax' < 1 {
            di as error "`cmd': draws(stata) needs " ceil(2 * `U' / _N) " variables of uniforms " ///
                "for one replicate (the frequency weights stand for `U' observations); use draws(fast)"
            exit 900
        }
        local chunk = min(`chunk', `nmax')
    }
    local rows
    if `"`savefile'"' != "" {
        local chunk = min(`chunk', c(matsize))
        local rows " rows"
        tempname post
        local vars
        forvalues k = 1/`kb' {
            local vars `vars' _r`k'
        }
        local everyopt
        if "`every'" != "" local everyopt every(`every')
        postfile `post' `vars' using `"`savefile'"', double `everyopt'
    }
    * The plugin, as this file sees it (a plugin program is local to the
    * ado-file that defines it; the library and its retained fit are shared)
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if `dots' _dots 0, title(Bootstrap replications) reps(`reps')
    * draws(stata): each chunk's uniforms, packed into variables reused
    * across chunks
    * Own every scratch variable individually; tempvars also clean up on Break.
    local uniform_vars
    local first 0
    * Any failure or Break between chunks (the uniforms, post, _dots) leaves
    * neither the uniform variables nor the rows matrix behind
    capture noisily {
        if "`draws'" == "stata" {
            forvalues j = 1/`=ceil(2 * `chunk' * `U' / _N)' {
                tempvar uniform
                quietly generate double `uniform' = .
                local uniform_vars `uniform_vars' `uniform'
            }
        }
        while `first' < `reps' {
            local n = min(`chunk', `reps' - `first')
            matrix __`cmd'_bs_status = J(1, `n', .)
            if "`rows'" != "" matrix __`cmd'_bs_rows = J(`n', `kb', .)
            local vars
            if "`draws'" == "stata" {
                Uniforms, units(`U') reps(`n') variables(`uniform_vars') ///
                    `=cond("`strataobs'" != "", "strataobs(`cmd')", "")'
                local vars `r(varlist)'
            }
            else if "`draws'" == "supplied" {
                forvalues r = `=`first' + 1'/`=`first' + `n'' {
                    local vars `vars' `: word `r' of `bsdraws''
                }
            }
            capture noisily plugin call ctools_plugin `vars', ///
                "`cmd' `threads' boot_chunk first=`first' count=`n' total=`reps'`how'`targets'`rows'"
            local rc = _rc
            if !`rc' & "`rows'" != "" {
                forvalues r = 1/`n' {
                    local exp
                    forvalues k = 1/`kb' {
                        local exp `exp' (__`cmd'_bs_rows[`r', `k'])
                    }
                    post `post' `exp'
                }
            }
            if `rc' {
                if "`uniform_vars'" != "" capture drop `uniform_vars'
                capture matrix drop __`cmd'_bs_rows
                exit `rc'
            }
            if `dots' {
                forvalues r = 1/`n' {
                    _dots `=`first' + `r'' `=(__`cmd'_bs_status[1, `r'] != 0)'
                }
            }
            local first = `first' + `n'
        }
    }
    local rc = _rc
    if `rc' {
        capture drop `uniform_vars'
        capture matrix drop __`cmd'_bs_rows
        if "`rows'" != "" capture postclose `post'
        exit `rc'
    }
    if "`uniform_vars'" != "" capture drop `uniform_vars'
    if "`rows'" != "" {
        postclose `post'
        capture matrix drop __`cmd'_bs_rows
    }
    if `dots' di as text ""
    return local rngstate `"`rngstate'"'
end

program define SaveFinish
    syntax, TMP(string) FILE(string) MAP(name) [DOUBLE REPLACE COMMAND(string) ///
        SEED(string) N(string) CLUSTER(string) NCLUST(string) STRATA(string) ///
        NSTRATA(string) CMDNAME(name)]
    tempname b
    matrix `b' = e(b)
    local K = colsof(`b')
    if colsof(`map') != `K' {
        di as error "_ctools_hdfe_boot savefinish: map() does not match e(b)"
        exit 198
    }
    local colnames : colnames `b'
    * Variable names as bootstrap gives them (_prefix_expand): _b_name, or
    * _bs_# where that is not a valid name
    local names
    forvalues j = 1/`K' {
        local na : word `j' of `colnames'
        if bsubstr("`na'", 1, 1) != "_" local na _`na'
        local id _b`na'
        if `: list id in names' local id `id'`j'
        capture confirm name `id'
        if _rc local id _bs_`j'
        local names `names' `id'
    }
    local type = cond("`double'" != "", "double", "float")
    preserve
    quietly use `"`tmp'"', clear
    forvalues j = 1/`K' {
        local id : word `j' of `names'
        local src = `map'[1, `j']
        * A replicate that failed is missing in every term; the others hold
        * zero for the terms fixed at zero
        if `src' > 0 quietly generate `type' `id' = _r`src'
        else quietly generate `type' `id' = cond(missing(_r1), ., 0)
    }
    keep `names'
    order `names'
    label data `"bootstrap: `cmdname'"'
    char _dta[command] `"`command'"'
    char _dta[seed] `"`seed'"'
    if "`nclust'" != "" char _dta[N_cluster] `nclust'
    if "`cluster'" != "" char _dta[cluster] `cluster'
    if "`strata'" != "" char _dta[strata] `strata'
    if "`nstrata'" != "" char _dta[N_strata] `nstrata'
    char _dta[N] `n'
    char _dta[bs_version] 3
    char _dta[k_eq] 1
    char _dta[k_extra] 0
    forvalues j = 1/`K' {
        local id : word `j' of `names'
        local na : word `j' of `colnames'
        char `id'[observed] `=`b'[1, `j']'
        label variable `id' `"`=usubstr("_b[`na']", 1, 80)'"'
        char `id'[expression] `"_b[`na']"'
        char `id'[coleq] _
        char `id'[colname] `na'
        char `id'[is_eexp] 1
    }
    if "`replace'" != "" quietly save `"`file'"', replace
    else quietly save `"`file'"'
    restore
end

program define BootLevel, rclass
    syntax, CMD(name) BSLEVEL(real) BSGIVEN(integer) [DIOPTS(string) CMDLEVEL(string)]
    local 0 `", `diopts'"'
    syntax [, Level(string) *]
    if `"`cmdlevel'"' == "" local cmdlevel `"`level'"'
    local rest `"`options'"'
    local lev `bslevel'
    if `"`cmdlevel'"' != "" {
        local 0 `", level(`cmdlevel')"'
        syntax [, Level(cilevel)]
        if `bsgiven' & `level' != `bslevel' {
            di as error "option level() of vce(bootstrap) conflicts with option level() of `cmd'"
            exit 198
        }
        local lev `level'
    }
    return scalar level = `lev'
    return local diopts `"`rest' level(`lev')"'
end

program define Seed, rclass
    syntax [, SEED(string) FAST]
    if `"`seed'"' != "" set seed `seed'
    return local rngstate `"`c(rngstate)'"'
    if "`fast'" != "" {
        return local seedhi = string(floor(runiform() * 2147483647), "%10.0f")
        return local seedlo = string(floor(runiform() * 2147483647), "%10.0f")
    }
end

program define Uniforms, rclass
    syntax, UNITS(integer) REPS(integer) [GENerate(name) VARIABLES(varlist numeric) STRATAOBS(name)]
    if `units' < 1 | `reps' < 1 {
        di as error "_ctools_hdfe_boot uniforms: units() and reps() must be positive"
        exit 198
    }
    Packed `generate' `=2 * `reps' * `units'', variables(`variables')
    local vars `r(varlist)'
    * Stata's sequence: per replicate the draws, then the tie-breaks, then
    * within strata the uniforms bsample takes when a unit was drawn twice
    mata: _ctools_boot_uniforms("`vars'", `units', `reps', "`strataobs'")
    return local varlist `vars'
end

* Double variables stub1 ... stubV holding n packed values (V = ceil(n /
* _N)): those that exist already are reused, the others created. Returns
* r(varlist).
program define Packed, rclass
    syntax [anything], [VARIABLES(varlist numeric)]
    if "`variables'" != "" {
        local n `anything'
        local V = ceil(`n' / _N)
        if `: word count `variables'' < `V' exit 198
        local vars
        forvalues v = 1/`V' {
            local vars `vars' `: word `v' of `variables''
        }
    }
    else {
        gettoken stub n : anything
        if "`stub'" == "" | "`n'" == "" exit 198
        local V = ceil(`n' / _N)
        local vars
        forvalues v = 1/`V' {
            capture confirm variable `stub'`v', exact
            if _rc quietly generate double `stub'`v' = .
            local vars `vars' `stub'`v'
        }
    }
    return local varlist `vars'
end

program define WildParse, rclass
    syntax, [SPEC(string)]
    * noci and ci are both options, which syntax cannot declare as a pair:
    * noci is taken out of the options as a word, after coefs() and
    * bootcluster(), whose variable names could include one called noci
    local 0 `", `spec'"'
    syntax [, COEFs(string) BOOTCLuster(string) *]
    local noci : subinstr local options "noci" "", word count(local n_noci)
    local 0 `", `noci'"'
    local noci = cond(`n_noci', "noci", "")
    syntax [, Reps(integer 999) SEED(string) WEIGHTtype(string) NONULL ///
        Ptype(string) NOENUMerate FEDGP(string) CI Level(cilevel) ///
        NODOTS DOTS(integer -1) _PROJected]
    if `reps' < 1 {
        di as error "wildbootstrap option reps() requires a positive integer"
        exit 198
    }
    if `"`seed'"' != "" {
        capture confirm integer number `seed'
        if _rc & substr(`"`seed'"', 1, 1) != "X" {
            di as error `"wildbootstrap option seed(`seed') invalid"'
            exit 198
        }
    }
    local type = lower(strtrim("`weighttype'"))
    if "`type'" == "" local type rademacher
    local found
    foreach t in rademacher mammen webb normal gamma {
        if strlen("`type'") >= 3 & "`type'" == substr("`t'", 1, strlen("`type'")) local found `t'
    }
    if "`found'" == "" {
        di as error "wildbootstrap option weighttype() must be rademacher, mammen, webb, normal or gamma"
        exit 198
    }
    local ptype = lower(strtrim("`ptype'"))
    if "`ptype'" == "" local ptype symmetric
    local code = cond("`ptype'" == "symmetric", 0, cond("`ptype'" == "equaltail", 1, ///
        cond("`ptype'" == "lower", 2, cond("`ptype'" == "upper", 3, .))))
    if `code' == . {
        di as error "wildbootstrap option ptype() must be symmetric, equaltail, lower or upper"
        exit 198
    }
    local fedgp = lower(strtrim("`fedgp'"))
    if "`fedgp'" == "" local fedgp full
    if "`fedgp'" == substr("partialled", 1, max(4, strlen("`fedgp'"))) local fedgp partialled
    if !inlist("`fedgp'", "full", "partialled") {
        di as error "wildbootstrap option fedgp() must be full or partialled"
        exit 198
    }
    if "`ci'" != "" & "`noci'" != "" {
        di as error "wildbootstrap options ci and noci may not be combined"
        exit 198
    }
    if `: word count `bootcluster'' > 1 {
        di as error "wildbootstrap option bootcluster() takes one term: one clustering dimension"
        exit 198
    }
    if "`nodots'" != "" & `dots' != -1 {
        di as error "wildbootstrap options nodots and dots() may not be combined"
        exit 198
    }
    if `dots' == -1 local dots = cond("`nodots'" != "", 0, 1)
    return scalar reps = `reps'
    return local seed `"`seed'"'
    return local weighttype `found'
    return scalar null = ("`nonull'" == "")
    return local coefs `"`coefs'"'
    return local ptype `ptype'
    return scalar ptype_code = `code'
    return scalar enumerate = ("`noenumerate'" == "")
    return local fedgp `fedgp'
    return scalar ci = ("`ci'" != "")
    return scalar noci = ("`noci'" != "")
    return scalar level = `level'
    return local bootcluster `"`bootcluster'"'
    return scalar dots = `dots'
    return scalar projected = ("`_projected'" != "")
end

program define WildDraw, rclass
    syntax, CLusters(integer) REPS(integer) TYPE(name) [NOENUMerate]
    local type = lower("`type'")
    if !inlist("`type'", "rademacher", "mammen", "webb", "normal", "gamma") {
        di as error "wild bootstrap weight type `type' not allowed: rademacher, mammen, webb, normal or gamma"
        exit 198
    }
    if `clusters' < 1 | `reps' < 1 {
        di as error "_ctools_hdfe_boot wilddraw: clusters() and reps() must be positive"
        exit 198
    }
    capture mata: mata drop _ctools_wild_v
    * boottest: enumerate when clusters * ln(2) < ln(reps) + 1e-6
    * (at most 30 clusters: the plugin enumerates up to 2^30 vectors)
    if "`type'" == "rademacher" & "`noenumerate'" == "" & `clusters' <= 30 & ///
        `clusters' * ln(2) < ln(`reps') + 1e-6 {
        return scalar enumerate = 1
        return scalar reps = 2^`clusters'
        exit
    }
    local C `clusters'
    local B1 = `reps' + 1
    * boottest's MakeWildWeights for its WRE (unit-variance) weights, m = 1
    if "`type'" == "rademacher" {
        mata: _ctools_wild_v = (runiform(`C', `B1') :< .5) * -2 :+ 1
    }
    else if "`type'" == "webb" {
        mata: _ctools_wild_v = sqrt(2 * ceil(runiform(`C', `B1') * 3)) :* ((runiform(`C', `B1') :>= .5) :- .5)
    }
    else if "`type'" == "mammen" {
        mata: _ctools_wild_v = (rdiscrete(`C', `B1', 1.727c9716ffb76X-001 \ 1.1b06d1d200914X-002) :- 1.5) * 1.1e3779b97f4a8X+001 :+ .5
    }
    else if "`type'" == "normal" {
        mata: _ctools_wild_v = rnormal(`C', `B1', 0, 1)
    }
    else {
        mata: _ctools_wild_v = rgamma(`C', `B1', 4, .5) :- 2
    }
    mata: _ctools_wild_v = _ctools_wild_v[|1, 2 \ ., .|]
    return scalar enumerate = 0
    return scalar reps = `reps'
end

program define WildWeights, rclass
    syntax, FIRST(integer) COUNT(integer) [GENerate(name) VARIABLES(varlist numeric)]
    capture mata: st_local("C", strofreal(rows(_ctools_wild_v)))
    if _rc {
        di as error "_ctools_hdfe_boot wildweights: no draws (wilddraw first)"
        exit 198
    }
    mata: st_local("B", strofreal(cols(_ctools_wild_v)))
    if `first' < 1 | `count' < 1 | `first' + `count' - 1 > `B' {
        di as error "_ctools_hdfe_boot wildweights: replicates `first' to `=`first'+`count'-1' out of range 1 to `B'"
        exit 198
    }
    if `C' > _N {
        di as error "_ctools_hdfe_boot wildweights: more clusters than observations"
        exit 198
    }
    Packed `generate' `=`count' * `C'', variables(`variables')
    local vars `r(varlist)'
    mata: _ctools_boot_wildweights("`vars'", `first', `count')
    return local varlist `vars'
end

program define WildClear
    syntax
    capture mata: mata drop _ctools_wild_v
end

mata:
/* Store x at the 0-based packed position q0 of the supplied variable list
 * (value q in variable floor(q / N) + 1, observation mod(q, N) + 1) */
void _ctools_boot_store(string scalar stub, real scalar q0, real colvector x)
{
    real scalar N, q, done, v, o, take
    string rowvector names

    names = tokens(stub)

    N = st_nobs()
    q = q0
    done = 0
    while (done < rows(x)) {
        v = floor(q / N) + 1
        o = mod(q, N) + 1
        take = min((rows(x) - done, N - o + 1))
        st_store((o, o + take - 1), names[v], x[|done + 1 \ done + take|])
        done = done + take
        q = q + take
    }
}

/* bsample's uniforms for reps replicates of units units, in Stata's order.
 * With cmd (observations resampled within strata), the units more the
 * prefix's bsample takes when a stratum drew a unit twice are drawn and
 * discarded: the strata are those the plugin of cmd saved (dupK of them;
 * missing: always) */
void _ctools_boot_uniforms(string scalar stub, real scalar units, real scalar reps,
    string scalar cmd)
{
    real scalar r, K, j, o, n, dup
    real colvector u, d, skip
    real rowvector off, size

    K = 0
    if (cmd != "") {
        K = st_numscalar("__" + cmd + "_bs_dupK")
        if (rows(K) == 0) K = .
        if (K < .) {
            off = J(1, K, .)
            size = J(1, K, .)
            for (j = 1; j <= K; j++) {
                off[j] = st_numscalar("__" + cmd + "_bs_dupo" + strofreal(j))
                size[j] = st_numscalar("__" + cmd + "_bs_dupn" + strofreal(j))
            }
        }
    }
    for (r = 0; r < reps; r++) {
        u = runiform(units, 1)
        _ctools_boot_store(stub, 2 * r * units, u)
        _ctools_boot_store(stub, (2 * r + 1) * units, runiform(units, 1))
        if (cmd == "") continue
        /* bsample's draw within a stratum of n: int(uniform() * _N + 1) */
        dup = (K == .)
        for (j = 1; j <= K & !dup; j++) {
            o = off[j]
            n = size[j]
            d = sort(floor(u[|o + 1 \ o + n|] :* n :+ 1), 1)
            dup = any(d[|2 \ n|] :== d[|1 \ n - 1|])
        }
        if (dup) skip = runiform(units, 1)
    }
}

/* Columns first .. first + count - 1 of _ctools_wild_v, one after another */
void _ctools_boot_wildweights(string scalar stub, real scalar first, real scalar count)
{
    external real matrix _ctools_wild_v
    real scalar C, j

    C = rows(_ctools_wild_v)
    for (j = 0; j < count; j++) {
        _ctools_boot_store(stub, j * C, _ctools_wild_v[., first + j])
    }
}
end
