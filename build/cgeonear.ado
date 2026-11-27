*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! C-accelerated nearest neighbors and geodetic distances (ctools)
*! original: geonear 2.0.3 (September 4, 2019) by Robert Picard
*! original: geodist 1.1.0 (June 18, 2019) by Robert Picard
*! original: see help cgeonear for details

* Mata helpers: what the plugin cannot reach through Stata's plugin
* interface (variable creation, labels, characteristics and notes). The
* distances and the neighbors are found in the plugin.
capture mata: mata drop _cgeonear_stash()
capture mata: mata drop _cgeonear_fetch()
capture mata: mata drop _cgeonear_chars()
capture mata: mata drop _cgeonear_set_chars()
capture mata: mata drop _cgeonear_merge_dta()
capture mata: mata drop _cgeonear_hold_labels()
capture mata: mata drop _cgeonear_define_labels()
capture mata: mata drop _cgeonear_hold_var()
capture mata: mata drop _cgeonear_apply_var()
mata:
/* Metadata held between the data of the two files, in one external */
void _cgeonear_stash(string scalar key, transmorphic value)
{
    pointer() scalar p

    if ((p = findexternal("__cgeonear_meta")) == NULL) {
        p = crexternal("__cgeonear_meta")
        *p = asarray_create()
    }
    asarray(*p, key, value)
}

transmorphic _cgeonear_fetch(string scalar key)
{
    pointer() scalar p

    p = findexternal("__cgeonear_meta")
    if (p == NULL) return(J(0, 2, ""))
    if (!asarray_contains(*p, key)) return(J(0, 2, ""))
    return(asarray(*p, key))
}

/* The characteristics of name (_dta or a variable), in Stata's order */
string matrix _cgeonear_chars(string scalar name)
{
    string rowvector c
    string matrix out
    real scalar j

    c = tokens(st_macroexpand("`" + ": char " + name + "[]" + "'"))
    out = J(cols(c), 2, "")
    for (j = 1; j <= cols(c); j++) {
        out[j, 1] = c[j]
        out[j, 2] = st_global(name + "[" + c[j] + "]")
    }
    return(out)
}

/* Define characteristics; Stata lists them newest first, so defining them
   from the end of the list leaves the list as it was */
void _cgeonear_set_chars(string scalar name, string matrix c)
{
    real scalar j

    for (j = rows(c); j >= 1; j--) st_global(name + "[" + c[j, 1] + "]", c[j, 2])
}

/* What merging the neighbor file's dataset characteristics does to those
   in memory (geonear's last step): its notes are appended unless the same
   text is already a note, its other characteristics are added where the
   data have none of that name */
void _cgeonear_merge_dta()
{
    string matrix c
    string scalar name, text
    real scalar j, k, n, isnote, dup

    c = _cgeonear_fetch("n_dta")
    n = strtoreal(st_global("_dta[note0]"))
    if (missing(n)) n = 0
    for (j = 1; j <= rows(c); j++) {
        name = c[j, 1]
        isnote = regexm(name, "^note[0-9]+$")
        if (name == "note0") continue
        if (isnote) {
            text = c[j, 2]
            dup = 0
            for (k = 1; k <= n; k++) {
                if (st_global("_dta[note" + strofreal(k) + "]") == text) {
                    dup = 1
                    break
                }
            }
            if (dup) continue
            n++
            st_global("_dta[note" + strofreal(n) + "]", text)
            st_global("_dta[note0]", strofreal(n))
            continue
        }
        if (st_global("_dta[" + name + "]") == "") st_global("_dta[" + name + "]", c[j, 2])
    }
}

/* The value labels attached to vars, kept under key */
void _cgeonear_hold_labels(string scalar key, string rowvector vars)
{
    string colvector names, text
    real colvector values
    string scalar lab
    real scalar j, n

    names = J(cols(vars), 1, "")
    n = 0
    for (j = 1; j <= cols(vars); j++) {
        lab = st_varvaluelabel(vars[j])
        if (lab == "" | anyof(names, lab)) continue
        if (!st_vlexists(lab)) continue
        names[++n] = lab
        st_vlload(lab, values, text)
        _cgeonear_stash(key + ":v:" + lab, values)
        _cgeonear_stash(key + ":t:" + lab, text)
    }
    _cgeonear_stash(key, n ? names[|1 \ n|] : J(0, 1, ""))
}

/* Define the labels held under key that the data do not have */
void _cgeonear_define_labels(string scalar key)
{
    string colvector names
    real scalar j

    names = _cgeonear_fetch(key)
    for (j = 1; j <= rows(names); j++) {
        if (st_vlexists(names[j])) continue
        st_vlmodify(names[j], _cgeonear_fetch(key + ":v:" + names[j]),
            _cgeonear_fetch(key + ":t:" + names[j]))
    }
}

/* What describes variable v: type, format, label, value label and
   characteristics, kept under key */
void _cgeonear_hold_var(string scalar key, string scalar v)
{
    _cgeonear_stash(key + ":desc", (st_vartype(v), st_varformat(v), st_varlabel(v),
        st_varvaluelabel(v)))
    _cgeonear_stash(key + ":chars", _cgeonear_chars(v))
}

void _cgeonear_apply_var(string scalar key, string scalar v)
{
    string rowvector desc

    desc = _cgeonear_fetch(key + ":desc")
    st_varformat(v, desc[2])
    st_varlabel(v, desc[3])
    if (desc[4] != "") st_varvaluelabel(v, desc[4])
    _cgeonear_set_chars(v, _cgeonear_fetch(key + ":chars"))
}
end

program define cgeonear, rclass
    version 14.1
    * Two forms. With a using file (or neighbors()), geonear's: the nearest
    * neighbors of each base location. With four coordinates, geodist's: the
    * distance between two points, in each observation or once.
    capture syntax anything [if] [in] [using/] [, Neighbors(string) *]
    if _rc local near = strpos(`" `macval(0)' "', " using ") > 0
    * neighbors() without a using file is geonear's form unless four
    * coordinates are given: geodist's syntax then rejects the option
    else local near = `"`using'"' != "" | (`"`neighbors'"' != "" & `: word count `anything'' != 4)
    if `near' {
        * Break is recognized only inside _cgeonear_impl. A failure leaves
        * the data as they were: the plugin frees what it holds, and the
        * wide-form variables made for the results are dropped.
        nobreak {
            macro drop CTOOLS_cgeonear_plugin
            local cgeonear_newvars
            capture noisily break _cgeonear_impl `macval(0)'
            local rc = _rc
            if `rc' {
                if "$CTOOLS_cgeonear_plugin" != "" {
                    capture program ctools_plugin, plugin using("$CTOOLS_cgeonear_plugin")
                    capture plugin call ctools_plugin, "cgeonear clear"
                }
                if "`cgeonear_newvars'" != "" capture drop `cgeonear_newvars'
            }
            macro drop CTOOLS_cgeonear_plugin
            capture mata: mata drop __cgeonear_meta
            _ctools_cleanup cgeonear
            if `rc' exit `rc'
        }
        exit
    }
    _ctools_cleanup cgeonear
    capture noisily _cgeonear_dist `macval(0)'
    local rc = _rc
    if !`rc' & "`cgeonear_dscalar'" != "" {
        return scalar distance = scalar(__cgeonear_dd)
        if "`cgeonear_diter'" != "" return scalar iterations = `cgeonear_diter'
    }
    _ctools_cleanup cgeonear
    if `rc' exit `rc'
end

program define _cgeonear_impl
    version 14.1
    syntax varlist(min=3 max=3) using/ , Neighbors(namelist min=3 max=3) ///
        [ ///
        MIles Ellipsoid a(string) f(string) Ignoreself RAdius(string) ///
        NEArcount(integer 1) Genstub(name) REport(numlist max=1 >=0) ///
        MAxbases(string) ///
        WIDe LOng Ops(numlist max=1 integer >0) WIThin(numlist max=1 >0) ///
        LImit(numlist max=1 integer >0) ///
        Verbose THReads(integer 0) ///
        ]
    quietly _ctools_clock __ct_start
    local __ct_lap `__ct_start'
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    _ctools_cleanup cgeonear

    * geonear's checks, in its order and with its messages. report(), ops()
    * and maxbases() steer geonear's regions; they are accepted and change
    * nothing here.
    if "`wide'" != "" & "`long'" != "" {
        dis as err "wide and long cannot be combined"
        exit 198
    }
    local wide = cond("`long'" == "", 1, 0)
    if `wide' {
        if `nearcount' < 1 {
            dis as err "value must be >= 1 : option nearcount(`nearcount')"
            exit 198
        }
        if "`genstub'" == "" local genstub nid
        if "`within'`limit'" != "" {
            dis as err "long form option(s) not allowed"
            exit 198
        }
    }
    else {
        if `nearcount' < 0 {
            dis as err "value must be >= 0 : option nearcount(`nearcount')"
            exit 198
        }
        if "`genstub'" != "" {
            dis as err "genstub(`genstub') not allowed in long form mode"
            exit 198
        }
    }
    local limitopt
    if "`limit'" != "" {
        if `limit' < `nearcount' {
            dis as err "limit(`limit') < nearcount(`nearcount')"
            exit 198
        }
        local limitopt limit=`limit'
    }
    if "`miles'" != "" local units mi
    else local units km
    local mi2km 1.609344
    local ignore = "`ignoreself'" != ""
    local minnbor = `nearcount' + `ignore'

    if "`ellipsoid'`a'`f'" != "" {
        if "`radius'" != "" {
            dis as err "Sphere and ellipsoid model cannot be combined"
            exit 198
        }
        if "`a'`f'" != "" {
            capture confirm number `a'
            local rc = _rc
            capture confirm number `f'
            if _rc | `rc' {
                dis as err "reference ellipsoid parameters incorrectly specified"
                dis as err "use option a() for the semi-major axis, in meters"
                dis as err "use option f() for the reciprocal of flattening ratio"
                dis as err "e.g. a(6378249.145) f(293.465) for Clarke 1880"
                exit 7
            }
            if `a' <= 0 {
                dis as err "semi-major axis (in meters) <= 0 : option a(`a')"
                exit 198
            }
            if `f' <= 0 {
                dis as err "reciprocal of flattening ratio <= 0 : option f(`f')"
                exit 198
            }
            if `f' <= 1 {
                dis as err "reciprocal of flattening ratio <= 1 : option f(`f')"
                exit 198
            }
            scalar __cgeonear_a = `a'
            scalar __cgeonear_f = 1 / `f'
        }
        else {
            * WGS 1984: semi-major axis and flattening
            scalar __cgeonear_a = 6378137
            scalar __cgeonear_f = 1 / 298.257223563
        }
        local sphere 0
        local modelopt
    }
    else {
        if "`radius'" == "" scalar __cgeonear_R = 6371
        else {
            capture confirm number `radius'
            if _rc {
                dis as err "was expecting a number: option radius(`radius')"
                exit 7
            }
            scalar __cgeonear_R = `radius'
            if scalar(__cgeonear_R) <= 0 {
                dis as err "radius must be > 0"
                exit 198
            }
        }
        local sphere 1
        local modelopt sphere
    }
    if "`within'" == "" scalar __cgeonear_within = 0
    else scalar __cgeonear_within = `within'
    if "`units'" == "mi" scalar __cgeonear_within = scalar(__cgeonear_within) * `mi2km'

    tokenize `varlist'
    local baseid `1'
    local baselat `2'
    local baselon `3'
    if "`: type `baseid''" == "strL" {
        dis as err "cgeonear: strL identifiers are not supported; recast `baseid' to str#"
        exit 109
    }

    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    global CTOOLS_cgeonear_plugin "`__ctools_plugin'"
    local threadopt
    if `threads' local threadopt threads(`threads')
    local timing = "`verbose'" != ""
    quietly _ctools_lap __time_setup __ct_lap `__ct_lap'

    * 1. The base locations: one per identifier with coordinates
    local cgeonear_fail
    _ctools_strw `varlist'
    capture plugin call ctools_plugin `varlist', "cgeonear `threadopt' base"
    if _rc _cgeonear_fail `=_rc' "`cgeonear_fail'" base `varlist'
    local baseN = scalar(__cgeonear_nbase)
    local basestr = substr("`: type `baseid''", 1, 3) == "str"
    quietly _ctools_lap __time_base __ct_lap `__ct_lap'

    * 2. The neighbor locations, read from the file in a frame of their own
    * (Stata 16 or later) or in place of the master data
    local longopt
    if !`wide' local longopt long baseid(`baseid')
    local usingopts names(`neighbors') threadopt(`threadopt') `longopt'
    if c(stata_version) >= 16 {
        tempname uframe
        frame create `uframe'
        frame `uframe': _cgeonear_using `"`macval(using)'"', `usingopts'
        frame drop `uframe'
    }
    else {
        preserve
        _cgeonear_using `"`macval(using)'"', `usingopts'
        restore
    }
    local nborN = scalar(__cgeonear_nnbor)
    if `nborN' < `minnbor' {
        dis as err "not enough observations to satisfy nearcount"
        dis as err "number of nbor locations = " `nborN'
        error 2001
    }
    * geonear compares the identifiers under ignoreself
    if `ignore' & `basestr' != `nborstr' error 109
    quietly _ctools_lap __time_using __ct_lap `__ct_lap'

    local copts `modelopt' nearcount=`nearcount' `limitopt'
    if `ignore' local copts `copts' ignoreself
    if "`units'" == "mi" local copts `copts' miles
    if `timing' local copts `copts' verbose

    if `wide' {
        * 3. Wide form: two new variables per neighbor, written for every
        * observation of the master data, which is then ordered by baseid
        local newvars
        local newtypes
        forvalues i = 1/`nearcount' {
            local j
            if `nearcount' > 1 local j `i'
            foreach v in `genstub'`j' `units'_to_`genstub'`j' {
                capture confirm name `v'
                if _rc {
                    dis as err "`v' invalid name"
                    exit 198
                }
                confirm new variable `v'
                if `: list v in newvars' {
                    dis as err "variable `v' already defined"
                    exit 110
                }
                local newvars `newvars' `v'
            }
            local newtypes `newtypes' `nbortype' double
        }
        mata: (void) st_addvar(tokens(st_local("newtypes")), tokens(st_local("newvars")), 1)
        c_local cgeonear_newvars `newvars'
        _ctools_strw `newvars'
        plugin call ctools_plugin `newvars', "cgeonear `threadopt' compute wide `copts'"
        quietly _ctools_lap __time_compute __ct_lap `__ct_lap'
        * geonear orders the data by baseid; csort (stable) only when the
        * observations are not in that order already
        if !scalar(__cgeonear_inorder) quietly csort `baseid'
        mata: _cgeonear_merge_dta()
        c_local cgeonear_newvars
        quietly _ctools_lap __time_store __ct_lap `__ct_lap'
    }
    else {
        * 3. Long form: one observation per base and neighbor pair, in place
        * of the master data
        local dname `units'_to_`nborid'
        capture confirm name `dname'
        if _rc {
            dis as err "`dname' invalid name"
            exit 198
        }
        if "`dname'" == "`baseid'" {
            dis as err "variable `dname' already defined"
            exit 110
        }
        plugin call ctools_plugin, "cgeonear `threadopt' compute long `copts'"
        local nout = scalar(__cgeonear_nout)
        quietly _ctools_lap __time_compute __ct_lap `__ct_lap'
        * What geonear's result keeps of the base data: the identifier with
        * everything that describes it, the data label, the dataset's
        * characteristics and notes. Of the value labels only baseid's:
        * geonear's intermediate files drop labels no variable uses
        mata: _cgeonear_hold_var("b_id", st_local("baseid"))
        mata: _cgeonear_stash("b_dta", _cgeonear_chars("_dta"))
        local keeplabs : value label `baseid'
        local dlabel : data label
        local btype : type `baseid'
        nobreak {
            quietly drop _all
            quietly label dir
            foreach l in `r(names)' {
                if !`: list l in keeplabs' label drop `l'
            }
            mata: _cgeonear_define_labels("n_labels")
            local rtypes `btype' `nbortype' double
            local rnames `baseid' `nborid' `dname'
            mata: (void) st_addvar(tokens(st_local("rtypes")), tokens(st_local("rnames")), 1)
            * the plugin writes the pairs in this order
            sort `baseid' `dname' `nborid'
            quietly set obs `nout'
            _ctools_strw `rnames'
            capture noisily plugin call ctools_plugin `rnames', "cgeonear `threadopt' store"
            local rc = _rc
            mata: _cgeonear_apply_var("b_id", st_local("baseid"))
            mata: _cgeonear_apply_var("n_id", st_local("nborid"))
            label data `"`dlabel'"'
            mata: _cgeonear_set_chars("_dta", _cgeonear_fetch("b_dta"))
            mata: _cgeonear_merge_dta()
            global S_FN
            global S_FNDATE
            if `rc' exit `rc'
        }
        quietly _ctools_lap __time_store __ct_lap `__ct_lap'
    }

    local ndist = scalar(__cgeonear_ndist)
    local nthreads = scalar(__cgeonear_threads)
    quietly _ctools_lap __time_end __ct_start `__ct_start'
    dis as txt _n "{hline 79}"
    dis as txt "Unique base locations   = " ///
        as res %-10.0fc `baseN' ///
        _col(42) as txt "Unique neighbor locations = " ///
        as res %-10.0fc `nborN'
    dis as txt "Bases * Neighbors  " _cont
    if `baseN' * `nborN' > 1e6 ///
        dis as txt "(M)  = " ///
        as res %-12.1fc `baseN' * `nborN' / 1e6 _cont
    else dis as txt "     = " ///
        as res %-12.0fc `baseN' * `nborN' _cont
    dis _col(42) as txt "Threads                   = " ///
        as res %-10.0fc `nthreads'
    dis as txt "Computed distances " _cont
    if `ndist' > 1e6 ///
        dis as txt "(M)  = " as res %-12.2fc `ndist' / 1e6 _cont
    else dis as txt "     = " as res %-12.0fc `ndist' _cont
    dis _col(42) as txt "Total run time (seconds)  = " ///
        as res %-10.3f `__time_end'
    dis "{hline 79}"

    if `timing' {
        di as text ""
        di as text "{hline 60}"
        di as text "cgeonear timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
        di as text "  Base locations:           " as result %8.4f `__time_base' " sec"
        di as text "  Neighbor file:            " as result %8.4f `__time_using' " sec"
        di as text "  Nearest neighbors:        " as result %8.4f `__time_compute' " sec"
        di as text "    Index:                  " as result %8.4f scalar(__cgeonear_time_index) " sec"
        di as text "    Search:                 " as result %8.4f scalar(__cgeonear_time_search) " sec"
        if `wide' {
            di as text "    Write:                  " as result %8.4f scalar(__cgeonear_time_store) " sec"
            di as text "  Order by baseid:          " as result %8.4f `__time_store' " sec"
        }
        else {
            di as text "  Result dataset:           " as result %8.4f `__time_store' " sec"
        }
        di as text "{hline 60}"
        di as text "  Total:                    " as result %8.4f `__time_end' " sec"
        di as text "{hline 60}"
    }
end

* Read the neighbor locations from the file, in the current (empty or
* preserved) data: geonear's checks, the plugin's location table, and what
* the result needs from the file. Sets nborid, nbortype and nborstr in the
* caller.
program define _cgeonear_using
    version 14.1
    syntax anything(name=files), names(string) [threadopt(string) long baseid(string)]
    gettoken using : files
    * The varlist-using form of use rescans macro characters in its file
    * name: such names are loaded whole, then the variables are selected
    local literal = strpos(`"`macval(using)'"', char(36)) | strpos(`"`macval(using)'"', char(96))
    if !`literal' quietly use `names' using `"`macval(using)'"', clear
    else {
        quietly use `"`macval(using)'"', clear
        keep `names'
    }
    unab full : `names'
    tokenize `full'
    local nid `1'
    local nlat `2'
    local nlon `3'
    if "`long'" != "" & "`nid'" == "`baseid'" {
        dis as err "identifier variables names must not be the same"
        exit 198
    }
    if "`: type `nid''" == "strL" {
        dis as err "cgeonear: strL identifiers are not supported; recast `nid' to str#"
        exit 109
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local cgeonear_fail
    _ctools_strw `nid' `nlat' `nlon'
    capture plugin call ctools_plugin `nid' `nlat' `nlon', "cgeonear `threadopt' using"
    if _rc {
        * geonear names the variables as the neighbors() option gave them
        * in wide form, and in full in long form
        local shown `names'
        if "`long'" != "" local shown `full'
        _cgeonear_fail `=_rc' "`cgeonear_fail'" nbor `shown'
    }
    mata: _cgeonear_stash("n_dta", _cgeonear_chars("_dta"))
    if "`long'" != "" {
        mata: _cgeonear_hold_var("n_id", st_local("nid"))
        mata: _cgeonear_hold_labels("n_labels", st_local("nid"))
    }
    c_local nborid `nid'
    c_local nbortype `: type `nid''
    c_local nborstr = substr("`: type `nid''", 1, 3) == "str"
end

* geonear_prep_locations' message for a check the plugin reports
program define _cgeonear_fail
    version 14.1
    args rc why loctype id lat lon
    if "`why'" == "notconst" dis as err "`lat' or `lon' not constant within `id' group"
    else if "`why'" == "noobs" dis as err "no observation for `loctype' locations"
    else if "`why'" == "string" dis as err "lat/lon variables not numeric for `loctype' locations"
    else if "`why'" == "lat" dis as err "`loctype' latitude var `lat' must be between -90 and 90"
    else if "`why'" == "lon" dis as err "`loctype' longitude var `lon' must be between -180 and 180"
    else error `rc'
    exit `rc'
end

* geodist's form: lat1 lon1 lat2 lon2, each a variable, a scalar or a
* number. With a variable, one distance per observation into generate();
* otherwise one distance, displayed and returned in r(). The plugin call
* is "cgeonear dist ..." (src/cgeonear/cgeonear_dist.c).
program define _cgeonear_dist
    version 14.1
    syntax anything [in] [if], [Generate(name) MIles Radius(string) Sphere ///
        Ellipsoid(string) Maxiter(integer 25) Verbose THReads(integer 0)]
    * The if is evaluated now, before anything replaces the r() or e()
    * results it may name
    local iftouse
    if `"`if'"' != "" {
        tempvar iftouse
        mark `iftouse' `if' `in'
    }
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    local timing = "`verbose'" != ""
    if `timing' quietly _ctools_clock __ct_start

    if ("`sphere'" != "" | "`radius'" != "") & "`ellipsoid'" != "" {
        dis as err "You must choose either a sphere or ellipsoid model"
        exit 198
    }
    * each lat1 lon1 lat2 lon2 is either a variable, a scalar, or a number;
    * coordinates are in signed decimal degrees, east and north positive
    tokenize `anything'
    if "`5'" != "" {
        dis as err "unexpected extra text: `5'"
        exit 198
    }
    if "`4'" == "" {
        dis as err "You must provide lat1 lon1 lat2 lon2"
        exit 198
    }

    * geodist checks the coordinates in order and stops at the first that
    * fails. A variable's range is checked on the observations still in the
    * sample once the variables before it are marked out, which the plugin
    * does in one pass; the first other failure is recorded here (errpos)
    * and reported after the checks of the variables before it.
    local hasvar 0
    local allmissing 0
    local kinds
    local vars
    local errpos 0
    local j 0
    forvalues i = 1/2 {
        foreach l in lat lon {
            local ++j
            if `errpos' continue
            local tok ``j''
            capture confirm string var `tok'
            if !_rc {
                local errpos `j'
                local errrc 109
                local errmsg "`tok' is a string variable"
                continue
            }
            capture confirm numeric var `tok'
            if !_rc {
                local hasvar 1
                unab full : `tok', max(1)
                local vars `vars' `full'
                local kinds `kinds'V
                continue
            }
            capture confirm number `tok'
            local isnum = _rc == 0
            capture confirm scalar `tok'
            local isscalar = _rc == 0
            if `isscalar' {
                tempname what
                capture scalar `what' = `tok' / 1
                if _rc {
                    local errpos `j'
                    local errrc 109
                    local errmsg "`tok' is a string scalar"
                    continue
                }
            }
            capture local is_missing = mi(`tok')
            if _rc local is_missing 0
            if !(`isnum' | `isscalar' | `is_missing') {
                local errpos `j'
                local errrc 198
                local errmsg "Was expecting a number or a numeric scalar instead of -`tok'-"
                continue
            }
            if `is_missing' local allmissing 1
            else {
                if "`l'" == "lat" & abs(`tok') > 90 {
                    local errpos `j'
                    local errrc 198
                    local errmsg "latitude `i' must be between -90 and 90"
                    continue
                }
                if "`l'" == "lon" & abs(`tok') > 180 {
                    local errpos `j'
                    local errrc 198
                    local errmsg "longitude `i' must be between -180 and 180"
                    continue
                }
            }
            scalar __cgeonear_dc`j' = `tok'
            local kinds `kinds'C
        }
    }

    if "`miles'" != "" {
        local km_to_miles / 1.609344
        local units miles
    }
    else local units km
    if "`sphere'" != "" | "`radius'" != "" {
        if "`radius'" == "" local radius 6371
        local model sphere
    }
    else {
        tempname a f
        if "`ellipsoid'" != "" {
            tokenize "`ellipsoid'", parse(" ,")
            capture confirm number `1'
            local rc = _rc
            if "`2'" == "," local 2 `3'
            capture confirm number `2'
            if _rc | `rc' {
                dis as err "the reference ellipsoid parameters (a,f^-1) are not numbers"
                dis as err "a = semi-major axis, in meters"
                dis as err "f^-1 = reciprocal of flattening ratio"
                dis as err "e.g. ellipsoid(6378249.145, 293.465) for Clarke 1880"
                exit 198
            }
            scalar `a' = `1'
            scalar `f' = 1 / `2'
        }
        else {
            * WGS 1984: semi-major axis and flattening
            scalar `a' = 6378137
            scalar `f' = 1 / 298.257223563
            local ename "WGS 1984"
        }
    }

    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')
    local modelopt `model' `=cond("`units'" == "miles", "miles", "")' bodies=`=`maxiter' - 1'
    if "`model'" == "sphere" scalar __cgeonear_dfac = `radius' `km_to_miles'
    else {
        scalar __cgeonear_da = `a'
        scalar __cgeonear_df = `f'
    }

    * The sample geodist starts from: if and in
    local markopt
    local inopt `in'
    if "`iftouse'" != "" {
        local markopt marker
        local markvar `iftouse'
    }

    if `errpos' {
        * the variables before the failing coordinate are checked first
        if `hasvar' {
            * the coordinates from the failing one on are not read
            local pad = 4 - length("`kinds'")
            forvalues k = `=5 - `pad''/4 {
                scalar __cgeonear_dc`k' = 0
            }
            local kinds `kinds'`=substr("CCCC", 1, `pad')'
            _ctools_strw `vars'
            plugin call ctools_plugin `vars' `markvar' `inopt', ///
                "cgeonear `threadopt' dist vars kinds=`kinds' check=`=`errpos' - 1' `modelopt' `markopt'"
            _cgeonear_dist_range `cgeonear_dviol'
        }
        dis as err "`errmsg'"
        exit `errrc'
    }

    if !`hasvar' {
        if "`if'`in'" != "" {
            dis as err "lat/lon are not variables, if or in option not allowed"
            exit 198
        }
        if "`generate'" != "" {
            dis as err "lat/lon are not variables, nothing to generate"
            exit 198
        }
        c_local cgeonear_dscalar 1
        if `allmissing' {
            scalar __cgeonear_dd = .
            exit 0
        }
        if "`model'" != "sphere" & `maxiter' <= 1 {
            * geodist reads sin_sigma, which no iteration has set
            dis as err "sin_sigma not found"
            exit 111
        }
        plugin call ctools_plugin, "cgeonear `threadopt' dist scalar `modelopt'"
        if "`model'" == "sphere" {
            dis as txt "Great-circle distance " ///
                "(haversine formula, radius of `radius'km) = " ///
                as res scalar(__cgeonear_dd) " `units'"
        }
        else {
            if "`cgeonear_dfailed'" == "1" {
                dis as err "Warning: failed to converge due to near-antipodal points"
            }
            else {
                dis as txt "`ename' ellipsoid(`=`a'',`=1/`f'') distance = " ///
                    as res scalar(__cgeonear_dd) " `units'"
            }
            c_local cgeonear_diter = scalar(__cgeonear_diter)
        }
        exit 0
    }

    * coordinates in variables: one distance per observation of the sample
    local out out
    local fail
    if "`generate'" == "" local out
    else {
        capture confirm new variable `generate'
        if _rc local out
    }
    local allopt
    if `allmissing' local allopt allmissing
    if "`out'" != "" {
        mata: (void) st_addvar("double", st_local("generate"), 1)
        local outvar `generate'
    }
    _ctools_strw `vars'
    capture noisily plugin call ctools_plugin `vars' `outvar' `markvar' `inopt', ///
        "cgeonear `threadopt' dist vars kinds=`kinds' check=4 `modelopt' `out' `markopt' `allopt'"
    local rc = _rc
    if `rc' | "`cgeonear_dviol'" != "0" {
        if "`out'" != "" capture drop `generate'
        if `rc' exit `rc'
        _cgeonear_dist_range `cgeonear_dviol'
    }
    if "`generate'" == "" {
        dis as err "you must specify gen(newvar) option"
        exit 198
    }
    if "`out'" == "" {
        * geodist's generate fails here: the name is in use
        confirm new variable `generate'
    }
    if `timing' {
        quietly _ctools_lap __time_total __ct_start `__ct_start'
        di as text "cgeonear: " as result %8.4f `__time_total' as text " sec"
    }
    if "`model'" != "sphere" & "`cgeonear_dnfail'" != "" & "`cgeonear_dnfail'" != "0" {
        dis as err "Warning: failed to converge due to near-antipodal points"
        dis as err "Replaced distance(s) with missing value .a"
        dis as err "Number of distance(s) affected = " `cgeonear_dnfail'
    }
end

* geodist's message for a variable outside its range: coordinate j (1-4)
program define _cgeonear_dist_range
    version 14.1
    args j
    if "`j'" == "" | "`j'" == "0" exit 0
    local i = ceil(`j' / 2)
    if mod(`j', 2) dis as err "latitude `i' must be between -90 and 90"
    else dis as err "longitude `i' must be between -180 and 180"
    exit 198
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* ctools: drop the scalars and matrices a command and its plugin leave behind
* _ctools_cleanup <command> [<command> ...] drops every scalar and matrix
* named __<command>_* or _<command>_*.
program define _ctools_cleanup
    version 14.1
    foreach command of local 0 {
        foreach pattern in __`command'_* _`command'_* {
            local found : all scalars "`pattern'"
            if "`found'" != "" capture scalar drop `found'
            local found : all matrices "`pattern'"
            if "`found'" != "" capture matrix drop `found'
        }
    }
end

* Private monotonic timestamp in seconds; no public Stata timer IDs.
program define _ctools_clock
    version 14.1
    args target
    tempname held
    _return hold `held'
    capture noisily {
        * Once registered, this plugin stays loaded for the Stata session.
        * Resolve its path and check its protocol only on first use (or
        * after discard), rather than searching the ado-path at every lap.
        capture plugin call ctools_clock_plugin, "clock"
        if _rc {
            _ctools_load
            capture program ctools_clock_plugin, plugin using("`__ctools_plugin'")
            if _rc != 0 & _rc != 110 exit _rc
            plugin call ctools_clock_plugin, "clock"
        }
    }
    local rc = _rc
    _return restore `held'
    if `rc' exit `rc'
    c_local `target' `ctools_clock_value'
end

* ctools: seconds since the previous lap of a caller-local wall clock
* _ctools_lap <target> <clock> <clock value> stores the elapsed seconds in the
* caller's local <target> and restarts the clock, a local the caller starts
* with _ctools_clock <clock> and passes back expanded:
*     _ctools_clock __ct_lap
*     ...
*     _ctools_lap __time_setup __ct_lap `__ct_lap'
program define _ctools_lap
    version 14.1
    args target clock start
    quietly _ctools_clock now
    c_local `target' = `now' - `start'
    c_local `clock' `now'
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
