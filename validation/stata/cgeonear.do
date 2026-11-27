* cgeonear (nearest neighbors) and cgeonear_dist (distances): named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cgeonear"


if "`__case'" == "cgeonear" {
* Reference tests require SSC geonear 2.0.3.
do "validation/support/validate_setup.do"
do "validation/support/esample_helpers.do"
do "validation/support/cgeonear_helpers.do"
which geonear
which cgeonear
tempfile bg cell glob_b glob_n dup panel sbase snbor ref tmp
quietly {

* ---------------------------------------------------------------------------
* Helpers
* ---------------------------------------------------------------------------
* Wide form against geonear. Rows of one baseid may come in any order from
* geonear; cgeonear keeps their order. Distances are compared exactly, or to
* tol (relative) when given (ellipsoid: geonear measures some pairs from the
* neighbor to the base).
capture program drop _cgn_wide
program define _cgn_wide
    syntax, master(string) using(string) vars(string) n(string) name(string) ///
        [opts(string) tol(real 0)]
    tempfile ref out
    gettoken baseid : vars
    use `"`master'"', clear
    gen long __obs = _n
    quietly geonear `vars' using `"`using'"', n(`n') `opts'
    mata: st_local("sig_ref", _cgn_sig(0))
    sort `baseid' __obs
    save `ref'
    use `"`master'"', clear
    gen long __obs = _n
    capture cgeonear `vars' using `"`using'"', n(`n') `opts'
    if _rc {
        test_fail "`name'" "cgeonear r(`=_rc')"
        exit
    }
    mata: st_local("sig", _cgn_sig(0))
    * ordered by baseid, observations of one baseid in their original order
    local sorted : sortedby
    tempvar pos
    gen long `pos' = _n
    sort `baseid' __obs
    capture assert `pos' == _n
    * marked as sorted by baseid when cgeonear reordered the data
    local order_ok = !_rc & inlist("`: word 1 of `sorted''", "", "`baseid'")
    if `"`sig'"' != `"`sig_ref'"' {
        test_fail "`name'" "metadata differ: `sig' vs `sig_ref'"
        exit
    }
    if !`order_ok' {
        test_fail "`name'" "not ordered by `baseid' with ties in their order (sortedby `sorted')"
        exit
    }
    drop `pos'
    if `tol' == 0 {
        capture cf _all using `ref', all
        if _rc test_fail "`name'" "results differ from geonear"
        else test_pass "`name'"
        exit
    }
    * identifiers exactly, distances to tol
    local dvars
    foreach stub in km_to_* mi_to_* {
        capture unab more : `stub'
        if !_rc local dvars `dvars' `more'
    }
    local ok 1
    foreach v of local dvars {
        rename `v' __c_`v'
    }
    merge 1:1 `baseid' __obs using `ref', nogenerate
    foreach v of local dvars {
        capture assert reldif(__c_`v', `v') <= `tol' | (missing(__c_`v') & missing(`v'))
        if _rc local ok 0
        drop `v' __c_`v'
    }
    * every other variable exactly (cf ignores the distances left in ref)
    capture cf _all using `ref'
    if _rc local ok 0
    if `ok' test_pass "`name'"
    else test_fail "`name'" "results differ from geonear beyond `tol'"
end

* Long form against geonear: same observations in the same order, the same
* metadata and sort marker
capture program drop _cgn_long
program define _cgn_long
    syntax, master(string) using(string) vars(string) n(string) name(string) [opts(string)]
    tempfile ref
    use `"`master'"', clear
    quietly geonear `vars' using `"`using'"', n(`n') long `opts'
    mata: st_local("sig_ref", _cgn_sig(1))
    local N_ref = _N
    save `ref'
    use `"`master'"', clear
    capture cgeonear `vars' using `"`using'"', n(`n') long `opts'
    if _rc {
        test_fail "`name'" "cgeonear r(`=_rc')"
        exit
    }
    mata: st_local("sig", _cgn_sig(1))
    if `"`sig'"' != `"`sig_ref'"' {
        test_fail "`name'" "metadata differ: `sig' vs `sig_ref'"
        exit
    }
    if _N != `N_ref' {
        test_fail "`name'" "`=_N' observations, geonear `N_ref'"
        exit
    }
    capture cf _all using `ref', all
    if _rc test_fail "`name'" "results differ from geonear"
    else test_pass "`name'"
end

* The same error code as geonear, and the data left as they were
capture program drop _cgn_error
program define _cgn_error
    syntax, master(string) cmd(string asis) name(string)
    use `"`master'"', clear
    capture geonear `cmd'
    local rc_ref = _rc
    use `"`master'"', clear
    tempfile before
    save `before'
    capture cgeonear `cmd'
    local rc = _rc
    if `rc' != `rc_ref' test_fail "`name'" "r(`rc'), geonear r(`rc_ref')"
    else if `rc' == 0 test_fail "`name'" "no error"
    else {
        capture cf _all using `before', all
        if _rc test_fail "`name'" "data changed after the error"
        else test_pass "`name'"
    }
end

* ---------------------------------------------------------------------------
* Data: geonear's help example (block groups, cell towers) with co-located
* towers, and a global set including the poles and the date line
* ---------------------------------------------------------------------------
clear
set seed 123456
set obs 3278
gen bgid = _n
gen double bglat = 37 + (41 - 37) * runiform()
gen double bglon = -109 + (109 - 102) * runiform()
save `bg'
clear
set obs 1000
gen ctid = _n
gen double ctlat = 37 + (41 - 37) * runiform()
gen double ctlon = -109 + (109 - 102) * runiform()
* three towers at one place, two more at another
replace ctlat = ctlat[10] in 11/12
replace ctlon = ctlon[10] in 11/12
replace ctlat = ctlat[500] in 900
replace ctlon = ctlon[500] in 900
label variable ctid "cell tower"
save `cell'

clear
set seed 2468
set obs 4000
gen long bid = _n
gen double blat = asin(2 * runiform() - 1) * 180 / _pi
gen double blon = 360 * runiform() - 180
replace blat = 90 in 1
replace blat = -90 in 2
replace blon = 180 in 3
replace blon = -180 in 4
replace blat = 0 in 5/6
save `glob_b'
clear
set obs 25000
gen long nid0 = 100000 + _n
gen double nlat = asin(2 * runiform() - 1) * 180 / _pi
gen double nlon = 360 * runiform() - 180
replace nlat = 89.9 in 1
replace nlat = -89.95 in 2
replace nlon = 179.99 in 3
replace nlon = -179.99 in 4
replace nlat = 0 in 5/9
save `glob_n'

* ---------------------------------------------------------------------------
* Wide form
* ---------------------------------------------------------------------------
_cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    name("wide: nearest neighbor")
foreach k in 2 3 5 {
    _cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
        opts(near(`k')) name("wide: near(`k')")
}
_cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    opts(near(3) miles genstub(tower)) name("wide: miles genstub()")
_cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    opts(radius(6378.137) near(2)) name("wide: radius()")
_cgn_wide, master(`cell') using(`cell') vars(ctid ctlat ctlon) n(ctid ctlat ctlon) ///
    opts(ignoreself) name("wide: ignoreself, co-located towers")
_cgn_wide, master(`cell') using(`cell') vars(ctid ctlat ctlon) n(ctid ctlat ctlon) ///
    opts(ignoreself near(4)) name("wide: ignoreself near(4)")
_cgn_wide, master(`cell') using(`cell') vars(ctid ctlat ctlon) n(ctid ctlat ctlon) ///
    opts(near(3)) name("wide: self included, ties at distance 0")
_cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    opts(ellipsoid) tol(1e-12) name("wide: ellipsoid")
_cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    opts(near(3) a(6377397.155) f(299.1528128) miles) tol(1e-12) name("wide: custom ellipsoid")
_cgn_wide, master(`glob_b') using(`glob_n') vars(bid blat blon) n(nid0 nlat nlon) ///
    opts(near(5)) name("wide: global, poles and date line")
_cgn_wide, master(`glob_b') using(`glob_n') vars(bid blat blon) n(nid0 nlat nlon) ///
    opts(near(2) ellipsoid) tol(1e-12) name("wide: global ellipsoid")
* compatibility options are accepted
_cgn_wide, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    opts(ops(500) report(1) maxbases(10)) name("wide: ops() report() maxbases()")

* ---------------------------------------------------------------------------
* Long form
* ---------------------------------------------------------------------------
use `cell', clear
rename ctid ctid0
label variable ctid0 "base tower"
save `tmp'
foreach opts in "" "ignoreself" "near(3)" "near(0) within(50) ignoreself" "within(20) near(2)" ///
    "within(20) near(2) limit(3)" "within(30) near(0) limit(5)" "within(25) miles" ///
    "near(4) ellipsoid" "within(30) ellipsoid miles" "radius(6000) within(40)" {
    _cgn_long, master(`tmp') using(`cell') vars(ctid0 ctlat ctlon) n(ctid ctlat ctlon) ///
        opts(`opts') name("long: `opts'")
}
_cgn_long, master(`bg') using(`cell') vars(bgid bglat bglon) n(ctid ctlat ctlon) ///
    opts(within(5)) name("long: nearest and all within 5 km")
_cgn_long, master(`glob_b') using(`glob_n') vars(bid blat blon) n(nid0 nlat nlon) ///
    opts(within(250) near(2)) name("long: global within(250)")
* every pair: a within() larger than any distance
use `cell' in 1/40, clear
rename ctid ctid0
save `ref', replace
use `cell' in 1/60, clear
save `dup', replace
_cgn_long, master(`ref') using(`dup') vars(ctid0 ctlat ctlon) n(ctid ctlat ctlon) ///
    opts(within(25000) near(0)) name("long: all pairs")

* ---------------------------------------------------------------------------
* Identifiers: repeated, missing, string, labelled; coordinates as float
* ---------------------------------------------------------------------------
clear
set seed 77
set obs 900
gen int pid = ceil(_n / 3)
gen float plat = 37 + 4 * runiform()
gen float plon = -109 + 7 * runiform()
bysort pid: replace plat = plat[1]
bysort pid: replace plon = plon[1]
gen year = 2000 + mod(_n, 3)
replace plat = . in 4
replace plon = . in 7/9
replace pid = . in 10/11
replace pid = .a in 12
gsort -year pid
label define pidlab 1 "first" 2 "second"
label values pid pidlab
label variable pid "panel id"
notes pid: a note
notes: base note
label data "panel"
char _dta[origin] base
save `panel'
_cgn_wide, master(`panel') using(`cell') vars(pid plat plon) n(ctid ctlat ctlon) ///
    opts(near(2)) name("wide: repeated, missing and labelled ids")
_cgn_long, master(`panel') using(`cell') vars(pid plat plon) n(ctid ctlat ctlon) ///
    opts(near(2) within(3)) name("long: repeated, missing and labelled ids")

use `cell', clear
gen str8 sid = "t" + string(ctid, "%04.0f")
replace sid = "" in 5
label variable sid "string tower id"
notes: nbor note
notes: base note
char _dta[origin] nbor
char _dta[extra] nbor-only
save `snbor'
use `bg', clear
gen str6 sbg = "b" + string(bgid)
save `sbase'
_cgn_wide, master(`sbase') using(`snbor') vars(sbg bglat bglon) n(sid ctlat ctlon) ///
    opts(near(2)) name("wide: string ids")
_cgn_wide, master(`bg') using(`snbor') vars(bgid bglat bglon) n(sid ctlat ctlon) ///
    name("wide: numeric base, string neighbor ids")
_cgn_long, master(`sbase') using(`snbor') vars(sbg bglat bglon) n(sid ctlat ctlon) ///
    opts(within(4)) name("long: string ids")
_cgn_long, master(`panel') using(`snbor') vars(pid plat plon) n(sid ctlat ctlon) ///
    opts(within(4) near(1)) name("long: notes and characteristics of both files")
use `snbor', clear
rename sid sid0
save `tmp', replace
_cgn_long, master(`tmp') using(`snbor') vars(sid0 ctlat ctlon) n(sid ctlat ctlon) ///
    opts(ignoreself within(2)) name("long: string ids, ignoreself")
* abbreviated neighbor variables, file name with spaces
use `cell', clear
local spaced "`c(tmpdir)'/cgeonear test cell.dta"
save "`spaced'", replace
_cgn_wide, master(`bg') using(`spaced') vars(bgid bglat bglon) n(ctid ctla ctlo) ///
    name("wide: abbreviations and a file name with spaces")
use `bg', clear
rename bgid ctidx
save `tmp', replace
_cgn_long, master(`tmp') using(`spaced') vars(ctidx bglat bglon) n(cti ctla ctlo) ///
    opts(within(3)) name("long: abbreviations")
erase "`spaced'"

* ---------------------------------------------------------------------------
* Errors: geonear's return codes, data unchanged
* ---------------------------------------------------------------------------
local errs `""bgid bglat bglon using `cell', n(ctid ctlat ctlon) near(2000)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) wide long""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) within(5)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) limit(5)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) long genstub(x)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) long near(3) limit(2)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) near(0)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) long near(-1)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) radius(6000) ellipsoid""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) a(6378137)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) a(-1) f(298)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) a(6378137) f(-3)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) radius(abc)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) radius(-1)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(bgid ctlat ctlon) long""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon ctid)""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctidx ctlat ctlon)""'
local errs `"`errs' "bgid bglat bglon using `cell'_nofile, n(ctid ctlat ctlon)""'
local errs `"`errs' "bgid bglat bglon using `snbor', n(sid ctlat ctlon) ignoreself""'
local errs `"`errs' "bgid bglat bglon using `cell', n(ctid ctlat ctlon) ignoreself near(1000)""'
local i 0
foreach e of local errs {
    local ++i
    _cgn_error, master(`bg') cmd(`e') name("error `i': `e'")
}
foreach change in "replace bglat = 95 in 5" "replace bglon = -181 in 9" ///
    "replace bglat = ." "replace bgid = 2 in 3" "drop if _n" "tostring bglon, replace force" {
    use `bg', clear
    `change'
    save `tmp', replace
    _cgn_error, master(`tmp') cmd(bgid bglat bglon using `cell', n(ctid ctlat ctlon)) ///
        name("base data: `change'")
}
foreach change in "replace ctlat = -91 in 5" "replace ctlon = 181 in 9" ///
    "replace ctlat = ." "replace ctid = 2 in 3" "replace ctid = ." {
    use `cell', clear
    `change'
    save `tmp', replace
    _cgn_error, master(`bg') cmd(bgid bglat bglon using `tmp', n(ctid ctlat ctlon)) ///
        name("neighbor data: `change'")
}
use `cell', clear
tostring ctlat, replace force
save `tmp', replace
_cgn_error, master(`bg') cmd(bgid bglat bglon using `tmp', n(ctid ctlat ctlon)) ///
    name("neighbor data: string latitude")

* ---------------------------------------------------------------------------
* Messages
* ---------------------------------------------------------------------------
use `bg', clear
replace bglat = 95 in 5
capture noisily cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon)
if _rc == 198 test_pass "message: latitude range"
else test_fail "message: latitude range" "r(`=_rc')"

* ---------------------------------------------------------------------------
* Where cgeonear differs (documented)
* ---------------------------------------------------------------------------
use `bg', clear
gen nid = -1
tempfile before
save `before'
capture cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon)
local rc = _rc
capture cf _all using `before', all
if `rc' == 110 & !_rc & c(k) == 4 test_pass "difference: an existing result variable is r(110)"
else test_fail "difference: an existing result variable is r(110)" "r(`rc')"
use `bg', clear
gen _merge = 1
set seed 5
gen double __u = runiform()
sort __u
drop __u
capture cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon)
if !_rc & "`: sortedby'" == "bgid" test_pass "difference: _merge allowed, reordered data sorted by baseid"
else test_fail "difference: _merge allowed, reordered data sorted by baseid" "r(`=_rc') sortedby `: sortedby'"
* neighbors() without a using file is the nearest-neighbor form: geonear's error
use `bg', clear
capture geonear bgid bglat bglon, n(ctid ctlat ctlon)
local rc_ref = _rc
capture cgeonear bgid bglat bglon, n(ctid ctlat ctlon)
if _rc == `rc_ref' & _rc test_pass "form: neighbors() without using"
else test_fail "form: neighbors() without using" "r(`=_rc'), geonear r(`rc_ref')"
* data already in baseid order are not reordered, and keep their sort marker
use `bg', clear
sort bgid bglat
capture cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon)
local ok = !_rc & "`: sortedby'" == "bgid bglat"
capture assert bgid == _n
if `ok' & !_rc test_pass "data in baseid order stay as they are"
else test_fail "data in baseid order stay as they are" "sortedby `: sortedby'"
use `cell', clear
generate strL sl = string(ctid)
save `tmp', replace
use `bg', clear
capture cgeonear bgid bglat bglon using `tmp', n(sl ctlat ctlon)
if _rc == 109 test_pass "difference: strL identifiers are r(109)"
else test_fail "difference: strL identifiers are r(109)" "r(`=_rc')"

* ---------------------------------------------------------------------------
* Threads and verbose: the same results at any thread count
* ---------------------------------------------------------------------------
foreach form in "" "long within(30)" {
    use `glob_b', clear
    cgeonear bid blat blon using `glob_n', n(nid0 nlat nlon) near(3) `form' threads(1)
    save `ref', replace
    use `glob_b', clear
    noisily cgeonear bid blat blon using `glob_n', n(nid0 nlat nlon) near(3) `form' threads(4) verbose
    capture cf _all using `ref', all
    if _rc test_fail "threads(1) and threads(4): `form'" "results differ"
    else test_pass "threads(1) and threads(4): `form'"
}
use `bg', clear
capture cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon) threads(300)
if _rc == 198 test_pass "threads() range"
else test_fail "threads() range" "r(`=_rc')"

* ---------------------------------------------------------------------------
* e(sample) as geonear leaves it
* ---------------------------------------------------------------------------
global CTOOLS_cgn_master `bg'
global CTOOLS_cgn_using `cell'
capture program drop _cgn_es_setup
program define _cgn_es_setup
    use "$CTOOLS_cgn_master", clear
    gen long __obs = _n
    quietly regress bglat bglon if mod(bgid, 3)
    estimates store m2
    quietly regress bglat bglon if mod(bgid, 2)
    estimates store m1
end
_ctools_esample_case, setup(_cgn_es_setup) ///
    native(geonear bgid bglat bglon using "$CTOOLS_cgn_using", n(ctid ctlat ctlon)) ///
    ctools(cgeonear bgid bglat bglon using "$CTOOLS_cgn_using", n(ctid ctlat ctlon)) ///
    sortby(bgid __obs) testname("e(sample) after cgeonear (wide)")
* the same with data that must be reordered by baseid
capture program drop _cgn_es_shuffled
program define _cgn_es_shuffled
    _cgn_es_setup
    set seed 11
    gen double __u = runiform()
    sort __u
    drop __u
end
_ctools_esample_case, setup(_cgn_es_shuffled) ///
    native(geonear bgid bglat bglon using "$CTOOLS_cgn_using", n(ctid ctlat ctlon) near(2)) ///
    ctools(cgeonear bgid bglat bglon using "$CTOOLS_cgn_using", n(ctid ctlat ctlon) near(2)) ///
    sortby(bgid __obs) testname("e(sample) after cgeonear (wide, reordered)")
global CTOOLS_cgn_master
global CTOOLS_cgn_using

* ---------------------------------------------------------------------------
* Break and failures leave the data as they were; the plugin keeps nothing
* ---------------------------------------------------------------------------
use `bg', clear
capture cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon) near(3) genstub(1x)
if _rc test_pass "invalid genstub() refused"
else test_fail "invalid genstub() refused" "no error"
capture cgeonear bgid bglat bglon using `cell', n(ctid ctlat ctlon)
if !_rc & c(k) == 5 test_pass "a call after a failure"
else test_fail "a call after a failure" "r(`=_rc') k=`c(k)'"
}
print_summary "cgeonear"
global CTOOLS_COMPONENT_COMPLETE "cgeonear"
    global CTOOLS_CASE_COMPLETE "cgeonear"
    exit
}

if "`__case'" == "cgeonear_dist" {
* cgeonear's distance form (four coordinates) against SSC geodist 1.1.0.
do "validation/support/validate_setup.do"
which geodist
which cgeonear
quietly {

* Variables: geodist and cgeonear into g and c. Distances agree to the last
* digit except where Stata's atan2 and the runtime's differ in the last bit
* (a few distances in 10,000, by at most a few units in the last place).
capture program drop _cgd_vars
program define _cgd_vars
    syntax anything, name(string) [opts(string) copts(string) ifin(string asis)]
    capture drop __g
    capture drop __c
    capture geodist `anything' `ifin', gen(__g) `opts'
    local rc_ref = _rc
    capture cgeonear `anything' `ifin', gen(__c) `opts' `copts'
    local rc = _rc
    if `rc' != `rc_ref' {
        test_fail "`name'" "r(`rc'), geodist r(`rc_ref')"
        exit
    }
    if `rc' {
        test_pass "`name'"
        exit
    }
    quietly count if __g != __c
    local ndiff = r(N)
    capture assert reldif(__g, __c) < 1e-14 | __g == __c
    local close = !_rc
    quietly count if __g == .a
    local na_ref = r(N)
    quietly count if __c == .a
    if !`close' | r(N) != `na_ref' | `ndiff' > ceil(_N / 1000) {
        test_fail "`name'" "`ndiff' distances differ (.a: `r(N)' vs `na_ref')"
    }
    else test_pass "`name'"
end

* The same return code as geodist and, when there is none, the same r()
capture program drop _cgd_scalar
program define _cgd_scalar
    syntax anything, name(string) [opts(string)]
    tempname d_ref
    capture geodist `anything', `opts'
    local rc_ref = _rc
    scalar `d_ref' = r(distance)
    local i_ref "`r(iterations)'"
    local r_ref : r(scalars)
    capture cgeonear `anything', `opts'
    local rc = _rc
    local same = (r(distance) == scalar(`d_ref')) | (missing(r(distance)) & ///
        missing(scalar(`d_ref')) & string(r(distance)) == string(scalar(`d_ref')))
    if `rc' != `rc_ref' test_fail "`name'" "r(`rc'), geodist r(`rc_ref')"
    else if !`rc' & (!`same' | "`r(iterations)'" != "`i_ref'" | "`: r(scalars)'" != "`r_ref'") {
        test_fail "`name'" "r(distance) `=r(distance)' vs `=scalar(`d_ref')', r(iterations) `r(iterations)' vs `i_ref'"
    }
    else test_pass "`name'"
end

clear
set seed 4321
set obs 30000
gen double lat1 = asin(2 * runiform() - 1) * 180 / _pi
gen double lon1 = 360 * runiform() - 180
gen double lat2 = lat1 + (runiform() - .5) * 30
replace lat2 = 90 * sign(lat2) if abs(lat2) > 90
gen double lon2 = lon1 + (runiform() - .5) * 60
replace lon2 = lon2 - 360 if lon2 > 180
replace lon2 = lon2 + 360 if lon2 < -180
replace lat1 = . in 1/10
replace lon2 = .b in 11/20
* antipodal points (Vincenty fails), coincident points, the equator, poles
replace lat2 = -lat1 in 21/60
replace lon2 = cond(lon1 > 0, lon1 - 180, lon1 + 180) in 21/60
replace lat2 = lat1 in 61/70
replace lon2 = lon1 in 61/70
replace lat1 = 0 in 71/80
replace lat2 = 0 in 71/80
replace lat1 = 90 in 81/85
replace lat2 = -90 in 86/90
gen float flat = lat1
gen float flon = lon1
foreach opts in "" "miles" "sphere" "sphere miles" "radius(6000)" "radius(6000) miles" ///
    "ellipsoid(6378249.145, 293.465)" "ellipsoid(6377563.396,299.3249646) miles" ///
    "maxiter(4)" "maxiter(2)" "maxiter(1)" {
    _cgd_vars lat1 lon1 lat2 lon2, opts(`opts') name("variables: `opts'")
}
_cgd_vars flat flon lat2 lon2, name("float variables")
_cgd_vars 42.265837 -83.748696 lat2 lon2, name("constants and variables")
_cgd_vars lat1 lon1 -33.9 151.2, opts(sphere) name("variables and constants, sphere")
scalar __s1 = 51.47803
scalar __s2 = -0.00145
_cgd_vars __s1 __s2 lat2 lon2, name("scalars and variables")
_cgd_vars lat1 lon1 lat2 lon2, ifin(if lon1 > 0 in 50/20000) name("if and in")
_cgd_vars lat1 lon1 lat2 lon2, ifin(in 100/200) name("in")
_cgd_vars lat1 lon1 lat2 lon2, ifin(if 0) name("empty sample")
_cgd_vars lat1 lon1 . lon2, name("a missing constant")
_cgd_vars lat1 lon1 lat2 lon2, opts(sphere) ifin(if lat1 < 0) name("if, sphere")
_cgd_vars lat1 lon1 lat2 lon2, copts(threads(1)) name("threads(1)")
_cgd_vars lat1 lon1 lat2 lon2, opts(sphere) copts(threads(3) verbose) name("threads(3) verbose")

* geodist's range checks: each variable on the observations whose
* variables before it are nonmissing
preserve
replace lat2 = 95 in 1
_cgd_vars lat1 lon1 lat2 lon2, name("out of range where lat1 is missing: no error")
replace lat2 = 95 in 30000
_cgd_vars lat1 lon1 lat2 lon2, name("lat2 out of range: error")
restore
preserve
replace lat1 = 100 in 11
_cgd_vars lat1 lon1 lat2 lon2, name("lat1 out of range where lon2 is missing: error")
restore
preserve
replace lon1 = 200 in 5
_cgd_vars lat1 lon1 lat2 lon2, name("lon1 out of range where lat1 is missing: no error")
replace lon2 = -181 in 22
replace lat1 = -91 in 23
_cgd_vars lat1 lon1 lat2 lon2, name("first failing coordinate reported")
restore
gen str3 s = "1"
scalar __str = "a"
foreach call in "lat1 lon1 lat2 s" "lat1 lon1 lat2 __str" "lat1 lon1 lat2 abc" "lat1 lon1 lat2 200" ///
    "lat1 lon1 95 lon2" "lat1 lon1 lat2" "lat1 lon1 lat2 lon2 lon1" "s lon1 lat2 lon2" {
    _cgd_vars `call', name("error: `call'")
}
_cgd_vars lat1 lon1 lat2 lon2, opts(sphere ellipsoid(1,2)) name("error: sphere and ellipsoid")
_cgd_vars lat1 lon1 lat2 lon2, opts(ellipsoid(a,b)) name("error: ellipsoid() not numbers")
capture cgeonear lat1 lon1 lat2 lon2
if _rc == 198 test_pass "error: no generate()"
else test_fail "error: no generate()" "r(`=_rc')"
capture cgeonear lat1 lon1 lat2 lon2, gen(lat1)
if _rc == 110 test_pass "error: generate() exists"
else test_fail "error: generate() exists" "r(`=_rc')"
capture cgeonear lat1 lon1 lat2 lon2, gen(new) threads(-1)
if _rc == 198 test_pass "error: threads()"
else test_fail "error: threads()" "r(`=_rc')"

* Scalars and numbers: displayed and returned in r()
_cgd_scalar 42.265837 -83.748696 42.207667 -83.356022, name("scalar: ellipsoid")
_cgd_scalar 42.265837 -83.748696 42.207667 -83.356022, opts(miles) name("scalar: miles")
_cgd_scalar 42.265837 -83.748696 42.207667 -83.356022, opts(sphere) name("scalar: sphere")
_cgd_scalar 42.265837 -83.748696 42.207667 -83.356022, opts(radius(6000) miles) name("scalar: radius()")
_cgd_scalar 51.47803 -0.00145 51.47803 0, name("scalar: Greenwich")
_cgd_scalar __s1 __s2 51.47803 0, opts(ellipsoid(6378249.145,293.465)) name("scalar: scalars, Clarke 1880")
_cgd_scalar 10 20 10 20, name("scalar: coincident points")
_cgd_scalar 0 10 0 50, name("scalar: equator")
_cgd_scalar 10 20 -10 -160, name("scalar: antipodal (no convergence)")
_cgd_scalar 10 20 -10 -159.5, name("scalar: nearly antipodal")
forvalues m = 2/6 {
    _cgd_scalar 42.265837 -83.748696 -12.5 100, opts(maxiter(`m')) name("scalar: maxiter(`m')")
}
_cgd_scalar 10 20 . 30, name("scalar: a missing coordinate")
_cgd_scalar 10 20 95 30, name("scalar: latitude out of range")
_cgd_scalar 10 200 5 30, name("scalar: longitude out of range")
_cgd_scalar 10 20 5 abc, name("scalar: not a number")
capture cgeonear 10 20 30 40 in 1/5
if _rc == 198 test_pass "scalar: in not allowed"
else test_fail "scalar: in not allowed" "r(`=_rc')"
capture cgeonear 10 20 30 40, gen(x)
if _rc == 198 test_pass "scalar: nothing to generate"
else test_fail "scalar: nothing to generate" "r(`=_rc')"
}
print_summary "cgeonear_dist"
global CTOOLS_COMPONENT_COMPLETE "cgeonear_dist"
    global CTOOLS_CASE_COMPLETE "cgeonear_dist"
    exit
}

di as error "Unknown cgeonear case: `__case'"
exit 198
