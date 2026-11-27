* cegen: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cegen"

if "`__case'" == "cegen" {
* Reference comparisons of cegen with egen: every function and option, with
* and without by-groups, if and in, on sorted and unsorted data. Each case
* compares the values, the storage type, the display format, the variable
* label and the value label of the new variable, or the return code when
* egen refuses the call.
do "validation/support/validate_setup.do"

* _cegen_cmp [by varlist] [tol(#)] : <egen arguments, NEW for the new
* variable and NEWL for a value label it names>
* tol(0), the default, requires identical values; tol(#) allows that
* relative difference where the order of a sum is not egen's.
capture program drop _cegen_cmp
program define _cegen_cmp
    gettoken head 0 : 0, parse(":")
    if `"`head'"' == ":" local head ""
    else gettoken colon 0 : 0, parse(":")
    local tol 0
    local pre ""
    foreach w of local head {
        if regexm("`w'", "^tol\((.*)\)$") local tol = regexs(1)
        else local pre `pre' `w'
    }
    if "`pre'" != "" local pre "`pre':"
    local a : subinstr local 0 "NEWL" "e_lab", all
    local b : subinstr local 0 "NEWL" "c_lab", all
    local a : subinstr local a "NEW" "e_out", all
    local b : subinstr local b "NEW" "c_out", all
    foreach v in e_out c_out {
        capture drop `v'
        capture label drop `v'
    }
    capture label drop e_lab
    capture label drop c_lab
    capture `pre' egen `a'
    local rc_e = _rc
    capture `pre' cegen `b'
    local rc_c = _rc
    local why ""
    if `rc_e' != `rc_c' local why "rc egen=`rc_e' cegen=`rc_c'"
    else if `rc_e' == 0 {
        if "`: type e_out'" != "`: type c_out'" local why "`why' type `: type e_out'/`: type c_out'"
        if "`: format e_out'" != "`: format c_out'" local why "`why' format"
        if `"`: variable label e_out'"' != `"`: variable label c_out'"' local why "`why' label"
        capture confirm string variable e_out
        if _rc == 0 {
            capture confirm string variable c_out
            if _rc local why "`why' kind"
            else {
                quietly count if e_out != c_out
                if r(N) local why "`why' `r(N)' values"
            }
        }
        else {
            capture confirm numeric variable c_out
            if _rc local why "`why' kind"
            else {
                if `tol' == 0 quietly count if e_out != c_out
                else quietly count if (reldif(e_out, c_out) > `tol' & !missing(e_out) & !missing(c_out)) | (missing(e_out) != missing(c_out))
                if r(N) local why "`why' `r(N)' values"
            }
        }
        local le : value label e_out
        local lc : value label c_out
        if ("`le'" == "") != ("`lc'" == "") local why "`why' value label"
        else if "`le'" != "" {
            capture decode e_out, gen(__de)
            capture decode c_out, gen(__dc)
            capture count if __de != __dc
            if _rc == 0 & r(N) local why "`why' `r(N)' label texts"
            capture drop __de
            capture drop __dc
        }
    }
    local shown : subinstr local 0 `"""' "'", all
    if `"`why'"' == "" test_pass "cegen `pre'`shown'"
    else test_fail "cegen `pre'`shown'" "`why'"
end

* _cegen_data #: the fixture, # observations in id order
capture program drop _cegen_data
program define _cegen_data
    args n
    clear
    set seed 20261126
    set obs `n'
    gen long id = _n
    gen int g = ceil(runiform() * 40)
    gen byte h = mod(_n, 3)
    gen double y = rnormal() * 100
    replace y = . if runiform() < .05
    replace y = .a if runiform() < .01
    gen float yf = y
    gen int k = ceil(runiform() * 6)
    replace k = . if runiform() < .03
    replace k = .b if runiform() < .01
    gen double z = round(runiform() * 10)
    gen str8 s = "s" + string(ceil(runiform() * 12))
    replace s = "" if runiform() < .04
    gen str12 name = word("alpha beta gamma delta", ceil(runiform() * 4)) + " " + word("one two three", h + 1)
    gen double a1 = rnormal()
    gen double a2 = rnormal()
    gen double a3 = rnormal()
    gen float a4 = rnormal()
    replace a1 = . if runiform() < .2
    replace a2 = .b if runiform() < .2
    replace a3 = . if runiform() < .5
    replace a4 = . if runiform() < .7
    gen byte i1 = ceil(runiform() * 5)
    gen byte i2 = ceil(runiform() * 5)
    gen byte i3 = ceil(runiform() * 5)
    replace i2 = . if runiform() < .1
    * one group whose values are all missing
    replace y = . if g == 7
    replace yf = . if g == 7
    label define klab 1 "one" 2 "two" 3 "three"
    label values k klab
end

quietly {
_cegen_data 5000

* ---------------------------------------------------------------- statistics
* Without by, egen takes sums from summarize: equal to rounding error.
foreach f in count mean total sd min max median iqr mad mdev kurt skew pc std {
    _cegen_cmp tol(1e-12) : NEW = `f'(y)
    _cegen_cmp tol(1e-12) : double NEW = `f'(y)
    _cegen_cmp tol(1e-12) : double NEW = `f'(y) if h == 1
    _cegen_cmp tol(1e-12) : double NEW = `f'(y) in 100/3000
    _cegen_cmp tol(1e-12) : double NEW = `f'(y * 2 + 1) if z > 3 in 50/4000
}
_cegen_cmp : NEW = count(s)
_cegen_cmp : NEW = count(s), by(g)
_cegen_cmp : NEW = count(s + "x") if h
_cegen_cmp : double NEW = min(y), missing
_cegen_cmp : double NEW = max(y), missing
_cegen_cmp : double NEW = min(y), missing by(g)
_cegen_cmp : double NEW = max(y), missing by(g)
_cegen_cmp : double NEW = max(y) if h == 2, missing by(g)
_cegen_cmp : double NEW = max(y), by(k)
_cegen_cmp tol(1e-12) : double NEW = total(y), missing
_cegen_cmp tol(1e-12) : double NEW = total(y) if g == 7, missing
_cegen_cmp tol(1e-9) : double NEW = total(y), missing by(g)
_cegen_cmp tol(1e-12) : double NEW = pc(y), prop
_cegen_cmp tol(1e-12) : double NEW = std(y), mean(10) sd(3)
_cegen_cmp : double NEW = std(y), sd(1) std(2)
_cegen_cmp : byte NEW = mean(y)
_cegen_cmp : int NEW = total(z)
_cegen_cmp : str5 NEW = mean(y)
foreach p in 1 5 10 25 33.3 50 66.6 75 90 99 99.9 {
    _cegen_cmp : double NEW = pctile(y), p(`p')
    _cegen_cmp : double NEW = pctile(y), p(`p') by(g)
    _cegen_cmp : double NEW = pctile(z), p(`p') by(g h)
}
_cegen_cmp : double NEW = pctile(y), p(0)
_cegen_cmp : double NEW = pctile(y), p(100)
_cegen_cmp : NEW = median(z), autotype
_cegen_cmp : NEW = median(y), autotype by(g)
_cegen_cmp : NEW = pctile(z), p(30) autotype

* With by on data sorted by the by-variables the values are egen's bit for bit.
sort g id
foreach f in count mean total sd min max median iqr mad mdev kurt skew pc std {
    _cegen_cmp : double NEW = `f'(y), by(g)
    _cegen_cmp : NEW = `f'(yf), by(g)
    _cegen_cmp : double NEW = `f'(y) if h == 1, by(g)
    _cegen_cmp by g : double NEW = `f'(y)
}
_cegen_cmp : double NEW = total(y[_n-1]), by(g)
_cegen_cmp : double NEW = mean(y) if _n <= 5, by(g)
_cegen_cmp by g : double NEW = total(_n)
_cegen_cmp by g (id) : double NEW = mean(y*_N)
_cegen_cmp : NEW = seq() if h == 1, from(1) to(3) by(g)
* Unsorted data: egen sorts first, so its sums run in another order.
sort id
foreach f in count mean total sd min max median iqr mad mdev kurt skew pc std {
    _cegen_cmp tol(1e-9) : double NEW = `f'(y), by(g)
    _cegen_cmp tol(1e-9) : double NEW = `f'(y), by(g h)
    _cegen_cmp tol(1e-9) : double NEW = `f'(y), by(k)
    _cegen_cmp tol(1e-9) : double NEW = `f'(y) if h == 1, by(s)
    _cegen_cmp tol(1e-9) : double NEW = `f'(y) in 200/4500, by(g)
    _cegen_cmp tol(1e-9) : double NEW = `f'(y + z) if h, by(s g)
}

* ---------------------------------------------------------------- rank, mode, seq
foreach o in "" field track {
    _cegen_cmp : NEW = rank(y), `o'
    _cegen_cmp : double NEW = rank(z), `o'
    _cegen_cmp : double NEW = rank(z) if h == 1, `o'
    _cegen_cmp : double NEW = rank(z), `o' by(g)
    _cegen_cmp : double NEW = rank(-z) if h, `o' by(g h)
    _cegen_cmp : double NEW = rank(y) in 10/4000, `o' by(s)
    _cegen_cmp : double NEW = rank(z), `o' by(k)
}
_cegen_cmp : double NEW = rank(id), unique
_cegen_cmp : double NEW = rank(a1), unique by(g)
_cegen_cmp : NEW = rank(y), field track
* unique ranks of tied values: egen breaks ties at random, cegen in order
cegen double u1 = rank(z), unique
cegen double u2 = rank(z), track
sort z id
by z: assert u1 == u2 + _n - 1 if !missing(z)
assert missing(u1) == missing(z)
sort id
test_pass "cegen rank unique: ties in observation order"
drop u1 u2
foreach o in "" minmode maxmode "nummode(1)" "nummode(2)" "nummode(40)" missing "missing minmode" "missing maxmode" {
    _cegen_cmp : NEW = mode(k), `o'
    _cegen_cmp : NEW = mode(z), `o'
    _cegen_cmp : NEW = mode(k), `o' by(g)
    _cegen_cmp : NEW = mode(z) if h == 1, `o' by(g)
    _cegen_cmp : NEW = mode(s), `o'
    _cegen_cmp : NEW = mode(s), `o' by(g)
    _cegen_cmp : NEW = mode(y), `o' by(g)
}
_cegen_cmp : NEW = mode(k), minmode maxmode
_cegen_cmp : NEW = seq()
_cegen_cmp : NEW = seq(), from(5) to(9)
_cegen_cmp : NEW = seq(), from(9) to(5) block(3)
_cegen_cmp : NEW = seq(), by(g)
_cegen_cmp : NEW = seq() in 10/200, from(250) to(700)
_cegen_cmp : NEW = seq(), from(-5) to(40000) block(2)
_cegen_cmp : NEW = seq(), block(0)
_cegen_cmp by g : NEW = seq(), to(4)

* ---------------------------------------------------------------- across variables
foreach f in rowtotal rowmean rowsd rowmin rowmax rowfirst rowlast rowmiss rownonmiss rowmedian {
    _cegen_cmp : NEW = `f'(a1 a2 a3 a4)
    _cegen_cmp : double NEW = `f'(a1 a2 a3 a4)
    _cegen_cmp : double NEW = `f'(a1-a4) if h == 1
    _cegen_cmp : double NEW = `f'(a*) in 5/2000
    _cegen_cmp : double NEW = `f'(a1 i1 a4 i2 y)
}
* Integer types: generate truncates and replace widens, step by step
foreach f in rowtotal rowmax rowmin rowfirst rowlast {
    foreach t in byte int long {
        _cegen_cmp : `t' NEW = `f'(a1 a2 a3 a4)
        _cegen_cmp : `t' NEW = `f'(i1 i2 i3 id)
        _cegen_cmp : `t' NEW = `f'(i1 a1 id y)
        _cegen_cmp : `t' NEW = `f'(i1 i2 i3)
    }
}
_cegen_cmp : NEW = rowtotal(a1 a2 a3 a4), missing
_cegen_cmp : double NEW = rowtotal(a1 a2 a3 a4), missing
_cegen_cmp : byte NEW = rowtotal(a1 a2 a3 a4), missing
_cegen_cmp : int NEW = rowtotal(i1 i2 id y), missing
_cegen_cmp : NEW = rowmiss(a1 s name)
_cegen_cmp : NEW = rownonmiss(a1 s name), strok
_cegen_cmp : NEW = rownonmiss(a1 s name)
_cegen_cmp : NEW = rowsd(a1)
foreach p in 1 10 25 50 60 75 99 {
    _cegen_cmp : double NEW = rowpctile(a1 a2 a3 a4 y z), p(`p')
}
_cegen_cmp : NEW = rowpctile(a1 a2), p(0)
_cegen_cmp : NEW = rowmedian(a1 s)
_cegen_cmp : NEW = anycount(i1 i2 i3), values(1 2)
_cegen_cmp : NEW = anycount(i1 i2 i3), values(1 2 2 5)
_cegen_cmp : NEW = anymatch(i1 i2 i3), values(3)
_cegen_cmp : NEW = anymatch(i1 i2 i3) if h, values(3/5)
_cegen_cmp : NEW = anyvalue(i1), values(2 4)
_cegen_cmp : double NEW = anyvalue(k) if h, values(1/3 6)
_cegen_cmp : NEW = anycount(i1 i2 i3)
_cegen_cmp : NEW = diff(i1 i2 i3)
_cegen_cmp : NEW = diff(i1 i1)
_cegen_cmp : NEW = diff(s name) in 1/1000
_cegen_cmp : NEW = diff(i1)
_cegen_cmp : NEW = rowmean(a1 a2), by(g)
* the names egen kept from before Stata 9
_cegen_cmp : NEW = rsum(a1 a2)
_cegen_cmp : NEW = rmean(a1 a2)
_cegen_cmp : NEW = robs(a1 a2)
_cegen_cmp : NEW = rmiss(a1 a2)
_cegen_cmp : NEW = rmax(a1 a2)
_cegen_cmp : NEW = rmin(a1 a2)
_cegen_cmp : NEW = rsd(a1 a2 a3)
_cegen_cmp : NEW = rfirst(a1 a2)
_cegen_cmp : NEW = rlast(a1 a2)
_cegen_cmp : NEW = neqany(i1 i2), values(1)
_cegen_cmp : NEW = eqany(i1 i2), values(1)
_cegen_cmp : NEW = any(i1), values(1 3)
_cegen_cmp tol(1e-12) : NEW = sum(y)

* ---------------------------------------------------------------- group, tag
_cegen_cmp : NEW = group(g)
_cegen_cmp : NEW = group(g h)
_cegen_cmp : long NEW = group(s g)
_cegen_cmp : NEW = group(k s), missing
_cegen_cmp : NEW = group(k y), missing
_cegen_cmp : NEW = group(k)
_cegen_cmp : NEW = group(g) if h == 1
_cegen_cmp : NEW = group(s) in 100/900
_cegen_cmp : NEW = group(k), label
_cegen_cmp : NEW = group(k s), label
_cegen_cmp : NEW = group(name h), label(NEWL)
_cegen_cmp : NEW = group(name h), label(NEWL, replace truncate(4))
_cegen_cmp : NEW = group(k), lname(NEWL)
_cegen_cmp : NEW = group(k), truncate(3)
_cegen_cmp : NEW = group(g), autotype
_cegen_cmp : NEW = group(id), autotype
_cegen_cmp : byte NEW = group(id)
_cegen_cmp : NEW = group()
_cegen_cmp : NEW = group(y), label
_cegen_cmp : NEW = tag(g)
_cegen_cmp : NEW = tag(g h s)
_cegen_cmp : NEW = tag(k s), missing
_cegen_cmp : NEW = tag(k)
_cegen_cmp : NEW = tag(g) if h == 2
_cegen_cmp : NEW = tag(s) in 1000/1100
_cegen_cmp by g : NEW = tag(g)
sort g h id
_cegen_cmp : NEW = group(g h)
_cegen_cmp : NEW = tag(g h)
sort id

* ---------------------------------------------------------------- cut, fill, ma
_cegen_cmp : NEW = cut(y), at(-300 -100 0 100 300)
_cegen_cmp : NEW = cut(y), at(-300 -100 0 100 300) icodes
_cegen_cmp : NEW = cut(y), at(-300 -100 0 100 300) label
_cegen_cmp : NEW = cut(z), at(0 .7 2.1 5 9.3)
_cegen_cmp : double NEW = cut(z), at(0 .7 2.1 5 9.3)
_cegen_cmp : byte NEW = cut(z), at(0 .7 2.1 5 9.3)
_cegen_cmp : NEW = cut(y), group(4)
_cegen_cmp : NEW = cut(y), group(7) label
_cegen_cmp : NEW = cut(z), group(3)
_cegen_cmp : NEW = cut(z), group(5) icodes
_cegen_cmp : NEW = cut(y) if h == 1, group(5) label
_cegen_cmp : NEW = cut(y) in 1/60, group(10)
_cegen_cmp : NEW = cut(y), group(1)
_cegen_cmp : NEW = cut(y)
_cegen_cmp : NEW = cut(y), at(1 2) group(2)
_cegen_cmp : NEW = cut(y*2), group(2)
forvalues n = 3/40 {
    _cegen_cmp : NEW = cut(id) in 1/`n', group(3)
    _cegen_cmp : NEW = cut(id) in 1/`n', group(4)
}
_cegen_cmp : NEW = fill(1 2)
_cegen_cmp : NEW = fill(10 8 6)
_cegen_cmp : NEW = fill(1 1 2 2)
_cegen_cmp : NEW = fill(0.1 0.2)
_cegen_cmp : double NEW = fill(0.1 0.2)
_cegen_cmp : byte NEW = fill(1 2)
_cegen_cmp : NEW = fill(1 3 5 9 2)
_cegen_cmp : NEW = fill(1)
_cegen_cmp : NEW = ma(z)
_cegen_cmp : NEW = ma(z), t(5)
_cegen_cmp : double NEW = ma(z), t(7) nomiss
_cegen_cmp : NEW = ma(z) in 10/400, t(5)
_cegen_cmp : NEW = ma(z*1.5), nomiss
_cegen_cmp : NEW = ma(y)
_cegen_cmp : NEW = ma(z), t(4)

* ---------------------------------------------------------------- strings
_cegen_cmp : NEW = concat(s name)
_cegen_cmp : NEW = concat(s name), punct(" ")
_cegen_cmp : NEW = concat(s name g), punct(", ")
_cegen_cmp : NEW = concat(g h k)
_cegen_cmp : NEW = concat(g h k), punct(-)
_cegen_cmp : NEW = concat(g k), decode punct(":")
_cegen_cmp : NEW = concat(g k), decode maxlength(2) punct(":")
_cegen_cmp : NEW = concat(y z), punct("|")
_cegen_cmp : NEW = concat(y z), format(%9.2f) punct("|")
_cegen_cmp : NEW = concat(s g) if h == 1, punct("_")
_cegen_cmp : NEW = concat(s id) in 20/300
_cegen_cmp : NEW = concat(id y)
_cegen_cmp : NEW = concat(s), format(%zz)
_cegen_cmp : NEW = ends(name)
_cegen_cmp : NEW = ends(name), head
_cegen_cmp : NEW = ends(name), tail
_cegen_cmp : NEW = ends(name), last
_cegen_cmp : NEW = ends(name), punct(a)
_cegen_cmp : NEW = ends(name), punct(a) tail
_cegen_cmp : NEW = ends(name), punct(a) last
_cegen_cmp : NEW = ends(name), punct(e) trim
_cegen_cmp : NEW = ends(name), punct(e) trim tail
_cegen_cmp : NEW = ends(name) if h == 1, punct("a ") last
_cegen_cmp : NEW = ends(name), head tail
_cegen_cmp : NEW = ends(g)
_cegen_cmp : NEW = ends(s) in 5/50

* ---------------------------------------------------------------- errors, egen's own functions
_cegen_cmp : NEW = nosuchfunction(y)
_cegen_cmp : NEW = mean(y), by(g) nosuchoption
_cegen_cmp : NEW = mtr(1990 y)
_cegen_cmp : NEW = mean(nosuchvar)
_cegen_cmp : g = mean(y)
capture cegen bad = mean(y), threads(-1)
assert _rc == 198
capture confirm variable bad
assert _rc == 111
test_pass "cegen: a refused call creates no variable"

* cegen leaves the order of the data and its sort marker alone
sort g id
cegen double m = mean(y), by(h)
assert "`: sortedby'" == "g id"
sort id
cegen double m2 = median(y), by(g)
assert id == _n
assert "`: sortedby'" == "id"
test_pass "cegen: observation order and sort marker unchanged"
drop m m2

* strL sources the plugin cannot read are egen's
gen strL big = s
replace big = 3000 * "x" in 1
_cegen_cmp : NEW = count(big), by(g)
_cegen_cmp : NEW = group(big)
_cegen_cmp : NEW = concat(big s)
drop big

* An empty dataset
clear
set obs 0
gen double y = .
gen str3 s = ""
cegen m = mean(y)
cegen c = concat(s)
assert _N == 0
confirm float variable m
confirm string variable c
test_pass "cegen: empty dataset"

* ---------------------------------------------------------------- parallel paths
* 250,000 observations run the threaded loops. Results must not depend on
* the number of threads, and must still be egen's.
_cegen_data 250000
foreach spec in "mean(y)" "sd(y)" "median(y)" "std(y)" "rowtotal(a1 a2 a3 a4)" ///
    "rowmedian(a1 a2 a3 a4)" "concat(s name)" "ends(name)" "cut(y), group(7)" "rank(z)" ///
    "tag(g h)" "group(s g)" "group(g h)" "mode(k)" {
    local comma = cond(strpos("`spec'", ","), "", ",")
    _cegen_cmp tol(1e-12) : double NEW = `spec'
    capture drop t1
    capture drop t4
    cegen double t1 = `spec'`comma' threads(1)
    cegen double t4 = `spec'`comma' threads(4)
    capture confirm string variable t1
    if _rc == 0 assert t1 == t4
    else assert t1 == t4
    test_pass "cegen threads(1) and threads(4) identical: `spec'"
}
foreach spec in "mean(y)" "total(y)" "sd(y)" "skew(y)" "median(y)" "pctile(y), p(10)" ///
    "iqr(y)" "mad(y)" "std(y)" "pc(y)" "min(y)" "count(y)" "rank(y)" "rank(z), field" ///
    "mode(k), minmode" "seq(), to(7)" {
    foreach groups in "g" "g h" "k" "s" {
        capture drop t1
        capture drop t4
        local comma = cond(strpos("`spec'", ","), "", ",")
        cegen double t1 = `spec'`comma' by(`groups') threads(1)
        cegen double t4 = `spec'`comma' by(`groups') threads(4)
        assert t1 == t4
        test_pass "cegen threads(1) and threads(4) identical: `spec' by `groups'"
    }
    _cegen_cmp tol(1e-9) : double NEW = `spec'`comma' by(g)
    _cegen_cmp tol(1e-9) : double NEW = `spec'`comma' by(s h)
}
sort g id
foreach spec in "mean(y)" "total(y)" "sd(y)" "skew(y)" "kurt(y)" "mdev(y)" "median(y)" ///
    "pctile(y), p(10)" "iqr(y)" "mad(y)" "std(y)" "pc(y)" "rank(y)" "mode(k), minmode" {
    local comma = cond(strpos("`spec'", ","), "", ",")
    _cegen_cmp : double NEW = `spec'`comma' by(g)
}
}
print_summary "cegen"
global CTOOLS_COMPONENT_COMPLETE "cegen"
    global CTOOLS_CASE_COMPLETE "cegen"
    exit
}

if "`__case'" == "cegen_extra" {
* More comparisons of cegen with egen: expressions as arguments, if with _n
* inside by-groups, if and in with by, string keys (unicode, quotes, empty
* values, strL, str2045), extended missing codes, extreme magnitudes,
* explicit storage types and their overflow, notes, error parity with the
* data left as it was, many groups, tiny and empty samples, abbreviations,
* _all, egen's version control, and results that do not depend on the
* order of the observations.
do "validation/support/validate_setup.do"

* _cegen_xcmp [by varlist] [tol(#)] [keeplab] : <egen arguments, NEW for
* the new variable and NEWL for a value label it names>
* As _cegen_cmp in the cegen case, and also: the sort marker, _N and the
* number of variables after cegen; on an error, no new variable and the
* data unchanged; notes; the number of value-label entries. keeplab leaves
* existing value labels e_lab, c_lab, e_out and c_out in place.
capture program drop _cegen_xcmp
program define _cegen_xcmp
    gettoken head 0 : 0, parse(":")
    if `"`head'"' == ":" local head ""
    else gettoken colon 0 : 0, parse(":")
    local tol 0
    local keeplab 0
    local pre ""
    foreach w of local head {
        if regexm("`w'", "^tol\((.*)\)$") local tol = regexs(1)
        else if "`w'" == "keeplab" local keeplab 1
        else local pre `pre' `w'
    }
    if "`pre'" != "" local pre "`pre':"
    local a : subinstr local 0 "NEWL" "e_lab", all
    local b : subinstr local 0 "NEWL" "c_lab", all
    local a : subinstr local a "NEW" "e_out", all
    local b : subinstr local b "NEW" "c_out", all
    foreach v in e_out c_out {
        capture drop `v'
        if !`keeplab' capture label drop `v'
    }
    if !`keeplab' {
        capture label drop e_lab
        capture label drop c_lab
    }
    capture `pre' egen `a'
    local rc_e = _rc
    quietly datasignature
    local sig1 "`r(datasignature)'"
    local sort1 : sortedby
    local n1 = _N
    local k1 = c(k)
    capture `pre' cegen `b'
    local rc_c = _rc
    local why ""
    if "`: sortedby'" != "`sort1'" local why "sortedby `sort1' to `: sortedby'"
    if _N != `n1' local why "`why' _N"
    if `rc_e' != `rc_c' local why "`why' rc egen=`rc_e' cegen=`rc_c'"
    else if `rc_e' {
        if c(k) != `k1' local why "`why' variables left behind"
        quietly datasignature
        if "`r(datasignature)'" != "`sig1'" local why "`why' data changed"
    }
    else {
        if c(k) != `k1' + 1 local why "`why' `=c(k) - `k1'' variables added"
        if "`: type e_out'" != "`: type c_out'" local why "`why' type `: type e_out'/`: type c_out'"
        if "`: format e_out'" != "`: format c_out'" local why "`why' format `: format e_out'/`: format c_out'"
        if `"`: variable label e_out'"' != `"`: variable label c_out'"' local why "`why' label"
        forvalues j = 0/2 {
            if `"`: char e_out[note`j']'"' != `"`: char c_out[note`j']'"' local why "`why' note`j'"
        }
        capture confirm string variable e_out
        if _rc == 0 {
            capture confirm string variable c_out
            if _rc local why "`why' kind"
            else {
                quietly count if e_out != c_out
                if r(N) local why "`why' `r(N)' values"
            }
        }
        else {
            capture confirm numeric variable c_out
            if _rc local why "`why' kind"
            else {
                if `tol' == 0 quietly count if e_out != c_out
                else quietly count if (reldif(e_out, c_out) > `tol' & !missing(e_out) & !missing(c_out)) | (missing(e_out) != missing(c_out)) | (missing(e_out) & e_out != c_out)
                if r(N) local why "`why' `r(N)' values"
            }
        }
        local le : value label e_out
        local lc : value label c_out
        if ("`le'" == "") != ("`lc'" == "") local why "`why' value label"
        else if "`le'" != "" {
            capture label list `le'
            local ke = r(k)
            capture label list `lc'
            if r(k) != `ke' local why "`why' label entries `ke'/`r(k)'"
            capture drop __de
            capture drop __dc
            capture decode e_out, gen(__de)
            capture decode c_out, gen(__dc)
            capture count if __de != __dc
            if _rc == 0 & r(N) local why "`why' `r(N)' label texts"
            capture drop __de
            capture drop __dc
        }
    }
    local shown : subinstr local 0 `"""' "'", all
    if `"`why'"' == "" test_pass "cegen_extra `pre'`shown'"
    else test_fail "cegen_extra `pre'`shown'" "`why'"
end

* _cegen_xall "fcn(args)[, opts]": egen and cegen with _all or * in the
* arguments, each run while the other's result does not exist
capture program drop _cegen_xall
program define _cegen_xall
    args spec
    capture drop e_out
    capture drop c_out
    tempfile ef
    local k0 = c(k)
    capture egen e_out = `spec'
    local rc_e = _rc
    if `rc_e' == 0 {
        preserve
        keep e_out
        quietly save `ef'
        restore
        drop e_out
    }
    capture cegen c_out = `spec'
    local rc_c = _rc
    local why ""
    if `rc_e' != `rc_c' local why "rc egen=`rc_e' cegen=`rc_c'"
    else if `rc_e' {
        if c(k) != `k0' local why "variables left behind"
    }
    else {
        quietly merge 1:1 _n using `ef', nogenerate
        if "`: type e_out'" != "`: type c_out'" local why "type `: type e_out'/`: type c_out'"
        quietly count if e_out != c_out
        if r(N) local why "`why' `r(N)' values"
    }
    capture drop e_out
    capture drop c_out
    if "`why'" == "" test_pass "cegen_extra `spec'"
    else test_fail "cegen_extra `spec'" "`why'"
end

* _cegen_xdata # seed: the fixture of this case
capture program drop _cegen_xdata
program define _cegen_xdata
    args n seed
    clear
    set seed `seed'
    set obs `n'
    gen long id = _n
    gen int g = ceil(runiform() * 40)
    gen byte h = mod(_n, 3)
    gen double y = rnormal() * 100
    replace y = . if runiform() < .05
    replace y = .a if runiform() < .01
    gen float yf = y
    gen double z = round(runiform() * 10)
    gen byte kx = ceil(runiform() * 6)
    foreach c in a b c d e f g h i j k l m n o p q r s t u v w x y z {
        replace kx = .`c' if runiform() < .003
    }
    replace kx = . if runiform() < .02
    label define kxl 1 "uno" 2 "dos" 6 "seis" .a "falta", replace
    label values kx kxl
    gen str10 s = "s" + string(ceil(runiform() * 15))
    replace s = "" if runiform() < .1
    gen str24 su = word("é ü ñ 東京 Ωmega naïve", ceil(runiform() * 6)) + " " + word("à ß ç", ceil(runiform() * 3))
    replace su = "" if runiform() < .05
    gen str12 sq = "q" + char(34) + string(ceil(runiform() * 3))
    replace sq = char(96) + "t" + char(39) + string(ceil(runiform() * 2)) if runiform() < .3
    gen strL sl = s
    replace sl = sl + " L" if h == 1
    gen str12 name = word("alpha beta gamma", ceil(runiform() * 3)) + " " + word("one two", ceil(runiform() * 2))
    gen double big = rnormal() * 1e306
    gen double tiny = rnormal() * 1e-300
    gen double cst = 5
    gen double shift = 1e9 + rnormal() * 1e-3
    gen byte bl = -127 + floor(runiform() * 228)
    gen int il = -32767 + floor(runiform() * 65508)
    gen long ll = -2147483647 + floor(runiform() * 4294967268)
    gen double bigint = cond(runiform() < .5, 1, -1) * (12345678 + floor(runiform() * 100))
    gen double pm = cond(mod(_n, 2), 1.5, -1.5)
    gen long pair = ceil(_n / 2)
    gen long key = ceil(runiform() * max(1, `n' / 2))
    gen double a1 = rnormal()
    gen double a2 = rnormal()
    gen float a3 = rnormal()
    replace a1 = . if runiform() < .2
    replace a2 = .c if runiform() < .2
    gen byte i1 = ceil(runiform() * 5)
    gen byte i2 = ceil(runiform() * 5)
    gen byte i3 = ceil(runiform() * 5)
    replace i2 = . if runiform() < .1
    gen double yvalue = y
    gen int grp = g
    gen byte a_very_long_variable_name_numb_1 = mod(_n, 2)
    gen byte a_very_long_variable_name_numb_2 = mod(_n, 3)
    gen byte a_very_long_variable_name_numb_3 = mod(_n, 5)
end

quietly {
_cegen_xdata 5000 20261127

* ---------------------------------------------------------------- expressions
* On data sorted by the by-variable the sums run in egen's order
sort g id
foreach f in count mean total sd min max median iqr mad mdev kurt skew pc std {
    _cegen_xcmp : double NEW = `f'(y*2 + z), by(g)
}
* Unsorted, by a unicode string key with empty values, if and in
sort id
foreach f in count mean total sd min max median iqr mad mdev kurt skew pc std {
    _cegen_xcmp tol(1e-9) : double NEW = `f'(abs(y) - z/3) if h in 50/4900, by(su)
}
_cegen_xcmp : double NEW = pctile(y^2 - z), p(17) by(sq)
_cegen_xcmp : NEW = count(cond(h, y, .)), by(s)
_cegen_xcmp : NEW = rank( y + z ), field
_cegen_xcmp : NEW = rank(y + z) if h, track by(g)
_cegen_xcmp : NEW = std( y / 2 ), mean(1) sd(2)
_cegen_xcmp tol(1e-12) : NEW = std((y)), by(g)
_cegen_xcmp : double NEW = rank(-yf), unique

* if with _n is evaluated within each group (sorted by the group)
sort g id
foreach f in count mean total sd min max median pctile iqr mad mdev kurt skew pc std rank {
    _cegen_xcmp : double NEW = `f'(y) if _n > 3 & _n < _N, by(g)
}
_cegen_xcmp : NEW = mode(kx) if _n <= 20, maxmode by(g)
_cegen_xcmp : NEW = seq() if _n > 2, by(g)
_cegen_xcmp by g : double NEW = total(y * _n) if mod(_n, 2)
_cegen_xcmp : double NEW = mean(y[_n+1] - y) in 10/4900, by(g)
_cegen_xcmp : double NEW = max(_n) if _n < 30 in 1/4000, by(g)

* ---------------------------------------------------------------- if and in with by
sort id
_cegen_xcmp tol(1e-9) : double NEW = mean(y) if h in 100/4000, by(s)
_cegen_xcmp : NEW = rank(y) if !h in 1/3000, track by(su s)
_cegen_xcmp : NEW = mode(z) if h in 1/4000, nummode(2) by(s)
_cegen_xcmp : NEW = count(y), by(sl)
_cegen_xcmp tol(1e-9) : double NEW = mean(y), by(sl kx)
_cegen_xcmp : double NEW = median(y), by(kx)
sort s id
_cegen_xcmp : double NEW = sd(y) if h in 100/4000, by(s)
_cegen_xcmp by s : double NEW = median(y) if h in 100/4000
_cegen_xcmp by s : NEW = seq() if h in 30/4000, block(2)
sort kx id
_cegen_xcmp : double NEW = mean(y), by(kx)
_cegen_xcmp by kx : double NEW = sd(y) if h
sort id
_cegen_xcmp by g : NEW = mean(y)

* ---------------------------------------------------------------- strings
_cegen_xcmp : NEW = concat(su s sq), punct(" | ")
_cegen_xcmp : NEW = concat(sq su)
_cegen_xcmp : NEW = concat(bigint il ll bl kx)
_cegen_xcmp : NEW = concat(yf z), punct(;)
_cegen_xcmp : NEW = concat(kx g), decode punct(/)
_cegen_xcmp : NEW = concat(sl s)
_cegen_xcmp : str3 NEW = concat(s name)
_cegen_xcmp : strL NEW = concat(s name)
_cegen_xcmp : double NEW = concat(s name)
_cegen_xcmp : NEW = ends(su)
_cegen_xcmp : NEW = ends(su), tail
_cegen_xcmp : NEW = ends(su), punct("ü") last
_cegen_xcmp : NEW = ends(sq), p(t) tr tail
_cegen_xcmp : str2 NEW = ends(name), last
_cegen_xcmp : NEW = ends(sl)
_cegen_xcmp : NEW = group(s)
_cegen_xcmp : NEW = group(s), missing
_cegen_xcmp : NEW = group(s su), missing label
_cegen_xcmp : NEW = group(su), label truncate(2)
_cegen_xcmp : NEW = group(sq g), label
_cegen_xcmp : NEW = group(sl h), label
_cegen_xcmp : NEW = group(kx), missing label
_cegen_xcmp : NEW = tag(s)
_cegen_xcmp : NEW = tag(s su), missing
_cegen_xcmp : NEW = tag(sl kx), missing
_cegen_xcmp : NEW = tag(kx)
_cegen_xcmp : NEW = mode(su), by(g) minmode
_cegen_xcmp : NEW = mode(sl), missing maxmode
_cegen_xcmp : NEW = mode(kx), missing nummode(3)
_cegen_xcmp : int NEW = mode(s)
_cegen_xcmp : str3 NEW = mode(kx)
_cegen_xcmp : NEW = rank(s)
_cegen_xcmp : NEW = mean(s)
_cegen_xcmp : NEW = mean(s + "1")
_cegen_xcmp : NEW = rowfirst(s su)
_cegen_xcmp : NEW = rowlast(s su sq)
_cegen_xcmp : str2 NEW = rowfirst(s name)
_cegen_xcmp : strL NEW = rowlast(s su)
_cegen_xcmp : NEW = rowfirst(s y)
_cegen_xcmp : NEW = rowmiss(s su sl)
_cegen_xcmp : NEW = rownonmiss(s su sl), strok
_cegen_xcmp : NEW = rownonmiss(s su)
_cegen_xcmp : NEW = rowtotal(s y)
_cegen_xcmp : NEW = diff(s su)
_cegen_xcmp : NEW = diff(kx i1)

* ---------------------------------------------------------------- missing option, missing codes
_cegen_xcmp : NEW = rowtotal(y kx), missing
_cegen_xcmp : NEW = rowtotal(y kx) if h, m
_cegen_xcmp tol(1e-9) : double NEW = total(y) if g == 7, missing by(h)
_cegen_xcmp : NEW = rowmax(kx a2)
_cegen_xcmp : NEW = rowmin(kx a2) in 1/2000
_cegen_xcmp : NEW = anyvalue(kx), values(1/3)
_cegen_xcmp : NEW = anymatch(kx i1), values(-1 0 2)
_cegen_xcmp : NEW = cut(kx), at(0 2 4 7)
_cegen_xcmp : NEW = cut(kx), group(3) label

* ---------------------------------------------------------------- extreme magnitudes
foreach f in mean total sd min max median skew kurt mdev {
    _cegen_xcmp tol(1e-12) : double NEW = `f'(big)
    _cegen_xcmp tol(1e-12) : double NEW = `f'(tiny)
}
_cegen_xcmp tol(1e-12) : double NEW = sd(shift)
sort g id
foreach f in mean total sd {
    if "`f'" != "sd" _cegen_xcmp : double NEW = `f'(big), by(g)
    _cegen_xcmp : double NEW = `f'(shift), by(g)
}
* sd() of values near 1e306 by groups (a documented difference): the
* squares overflow; egen's sum() skips them as missing and gives 0, cegen
* gives missing
capture drop c_out
cegen double c_out = sd(big), by(g)
capture assert missing(c_out)
if _rc == 0 test_pass "cegen_extra double NEW = sd(big), by(g): missing, as documented"
else test_fail "cegen_extra double NEW = sd(big), by(g): missing, as documented" "nonmissing values"
drop c_out
sort id
_cegen_xcmp : NEW = rowtotal(big big big)
_cegen_xcmp : double NEW = rowmean(big big)
_cegen_xcmp : double NEW = rowsd(tiny tiny y)
_cegen_xcmp : NEW = mean(ll)
_cegen_xcmp : long NEW = max(ll)
_cegen_xcmp : NEW = rowtotal(bl il ll)
_cegen_xcmp : long NEW = rowtotal(ll ll)
_cegen_xcmp : double NEW = cut(big), at(-1e307 0 1e307)
_cegen_xcmp : NEW = pc(pm), by(pair)
_cegen_xcmp : NEW = kurt(cst), by(g)
_cegen_xcmp : NEW = std(cst)
_cegen_xcmp : NEW = skew(z) if z == 3, by(g)

* ---------------------------------------------------------------- explicit types and overflow
_cegen_xcmp : byte NEW = count(y)
_cegen_xcmp : int NEW = count(y), by(h)
_cegen_xcmp : byte NEW = max(bl)
_cegen_xcmp : byte NEW = max(il)
_cegen_xcmp : int NEW = mean(y), by(g)
_cegen_xcmp : long NEW = total(big)
_cegen_xcmp : byte NEW = rank(y)
_cegen_xcmp : int NEW = rank(y), field
_cegen_xcmp : byte NEW = seq()
_cegen_xcmp : int NEW = seq(), from(32700)
_cegen_xcmp : byte NEW = tag(g)
_cegen_xcmp : double NEW = tag(g)
_cegen_xcmp : int NEW = group(id)
_cegen_xcmp : float NEW = group(id)
_cegen_xcmp : byte NEW = anycount(il bl), values(100)
_cegen_xcmp : int NEW = fill(0.5 1)
_cegen_xcmp : long NEW = fill(2147483000 2147483300)
_cegen_xcmp : byte NEW = cut(il), group(5)
_cegen_xcmp : byte NEW = pctile(il), p(90)
_cegen_xcmp : int NEW = mode(il), minmode
_cegen_xcmp : byte NEW = ma(il), t(5)
_cegen_xcmp : bogus NEW = mean(y)

* ---------------------------------------------------------------- errors: data left as it was
_cegen_xcmp : y = mean(z)
_cegen_xcmp : double y = mean(z), by(g)
sort g id
_cegen_xcmp by g : id = seq()
sort id
_cegen_xcmp by g : NEW = mean(y)
_cegen_xcmp : NEW = tag(g), by(h)
_cegen_xcmp : NEW = group(g), by(h)
_cegen_xcmp : NEW = rowtotal(y z), by(g)
_cegen_xcmp : NEW = cut(y), group(3) by(g)
_cegen_xcmp : NEW = concat(s su), by(g)
_cegen_xcmp : NEW = fill(1 2), by(g)
_cegen_xcmp : NEW = fill(1 2) in 1/10
_cegen_xcmp : NEW = mean(y), by(nosuch)
_cegen_xcmp : NEW = pctile(y), p(101)
_cegen_xcmp : NEW = cut(y), at(3 1)
_cegen_xcmp : NEW = cut(y), group(2.5)
* egen (pctile) limits group() by _N + 1, not by the size of the sample
_cegen_xcmp : NEW = cut(y) in 1/1, group(3)
_cegen_xcmp : NEW = cut(y) in 1/2, group(5)
_cegen_xcmp : NEW = cut(y) if id <= 3, group(9)
_cegen_xcmp : NEW = ma(z) if h
_cegen_xcmp : NEW = ma(z), t(2)
_cegen_xcmp : NEW = std(y), sd(-1)
_cegen_xcmp : NEW = rank(y), field track
_cegen_xcmp : NEW = mean(y) in 10/5
_cegen_xcmp : NEW = mean(y) in 1/6000
_cegen_xcmp : NEW = group(kx), lname(NEWL) label(NEWL)
_cegen_xcmp : NEW = mean(y y)
_cegen_xcmp : NEW = mean()
_cegen_xcmp : NEW = rowtotal()
_cegen_xcmp : NEW = mean y
sort g id
_cegen_xcmp by g : NEW = mean(y), by(h)
sort id
* value labels that already exist
label define e_lab 1 "old", replace
label define c_lab 1 "old", replace
_cegen_xcmp keeplab : NEW = group(s), label(NEWL)
_cegen_xcmp keeplab : NEW = group(s), label(NEWL, replace)
capture label drop e_out
capture label drop c_out
label define e_out 1 "old"
label define c_out 1 "old"
_cegen_xcmp keeplab : NEW = group(s), label
_cegen_xcmp keeplab : NEW = cut(y), group(3) label
capture label drop e_lab
capture label drop c_lab
capture label drop e_out
capture label drop c_out

* ---------------------------------------------------------------- notes, long labels
_cegen_xcmp : NEW = tag(a_very_long_variable_name_numb_1 a_very_long_variable_name_numb_2 a_very_long_variable_name_numb_3)
_cegen_xcmp : NEW = group(a_very_long_variable_name_numb_1 a_very_long_variable_name_numb_2 a_very_long_variable_name_numb_3)
_cegen_xcmp : NEW = anycount(a_very_long_variable_name_numb_1 a_very_long_variable_name_numb_2 a_very_long_variable_name_numb_3), values(1)
_cegen_xcmp : NEW = anyvalue(a_very_long_variable_name_numb_3), values(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)

* ---------------------------------------------------------------- abbreviations
_cegen_xcmp tol(1e-9) : NEW = mean(yval), by(gr)
_cegen_xcmp : NEW = mode(kx), minm
_cegen_xcmp : NEW = mode(z), n(2)
_cegen_xcmp : NEW = rank(y), f
_cegen_xcmp : NEW = rank(y), t by(g)
_cegen_xcmp : NEW = seq(), b(3) f(2) t(9)
_cegen_xcmp : NEW = cut(y), g(4) ic
_cegen_xcmp : NEW = ma(z), t(5) nom
_cegen_xcmp : NEW = rownonmiss(s y), s
_cegen_xcmp : NEW = anycount(i*), v(1 2)
_cegen_xcmp : NEW = concat(s g), p(-) d maxl(3) f(%3.0f)
_cegen_xcmp : NEW = group(s g), l lname(NEWL) t(3)
_cegen_xcmp : NEW = rowtotal(i1-i3), m
_cegen_xcmp : NEW = total(y), m by(h)
_cegen_xcmp : NEW = max(y), miss by(h)

* ---------------------------------------------------------------- _all and *
preserve
keep id g h s y kx
_cegen_xall "rowmiss(_all)"
_cegen_xall "rownonmiss(*), strok"
_cegen_xall "concat(_all), punct(-)"
_cegen_xall "rowmax(*)"
restore

* ---------------------------------------------------------------- version control
_cegen_xcmp : NEW = group(name h), label truncate(3)
version 16: _cegen_xcmp : NEW = group(name h), label
version 16: _cegen_xcmp : NEW = group(su s), label
version 16.1: _cegen_xcmp : byte NEW = group(id)
version 16: _cegen_xcmp : NEW = concat(g kx), decode

* ---------------------------------------------------------------- verbose and threads change nothing
foreach spec in "mean(y), by(su)" "median(y), by(sl)" "group(su s), missing" "concat(su sq)" "rank(y), by(s) field" {
    capture drop v1
    capture drop v2
    capture drop v3
    local comma = cond(strpos("`spec'", ","), "", ",")
    cegen v1 = `spec'
    cegen v2 = `spec'`comma' verbose threads(3)
    cegen v3 = `spec'`comma' threads(1)
    capture assert v1 == v2 & v1 == v3
    if _rc == 0 test_pass "cegen_extra verbose threads(): `spec'"
    else test_fail "cegen_extra verbose threads(): `spec'" "results differ"
    drop v1 v2 v3
}

* ---------------------------------------------------------------- long strings
_cegen_xdata 300 20261131
gen str2045 lw = s + 1990 * "w" + string(h)
gen strL lb = lw + lw
_cegen_xcmp : NEW = concat(lw s)
_cegen_xcmp : NEW = concat(lw), punct(x)
_cegen_xcmp : NEW = group(lw h), label
version 16: _cegen_xcmp : NEW = group(lw h), label
_cegen_xcmp : NEW = ends(lw), punct(w) tail
_cegen_xcmp : NEW = tag(lw)
_cegen_xcmp : NEW = mode(lw), maxmode
_cegen_xcmp : NEW = rownonmiss(lw s), strok
_cegen_xcmp tol(1e-12) : double NEW = mean(y), by(lw)
_cegen_xcmp : NEW = group(lb)
_cegen_xcmp : NEW = count(y), by(lb)
drop lw lb

* ---------------------------------------------------------------- samples with no observations
foreach spec in "mean(y)" "total(y), missing" "count(y)" "median(y)" "rank(y)" "mode(z)" ///
    "group(g)" "tag(g)" "seq()" "concat(s)" "rowtotal(a1 a2)" "cut(y), group(3)" "std(y)" ///
    "max(y), by(g)" "pctile(y), by(g)" "anymatch(i1), values(1)" "ends(name)" "sd(y), by(s)" {
    gettoken fn opts : spec, parse(",")
    _cegen_xcmp : NEW = `fn' if h == 99 `opts'
}

* ---------------------------------------------------------------- one and two observations
foreach n in 1 2 {
    _cegen_xdata `n' 2026113`n'
    foreach spec in "mean(y)" "sd(y)" "median(y)" "kurt(y)" "rank(y)" "std(y)" "pc(y)" ///
        "mode(z)" "group(s)" "tag(s)" "seq()" "concat(s su)" "ends(name)" "rowsd(a1 a2)" ///
        "cut(id), group(2)" "ma(z)" "fill(1 2)" "mean(y), by(g)" "iqr(y), by(s)" "group(s), label" {
        _cegen_xcmp : NEW = `spec'
    }
}
clear
set obs 0
gen double y = .
gen int g = .
gen str3 s = ""
_cegen_xcmp : NEW = group(s g), label
* egen itself stops with r(111) on an empty dataset with by(); cegen
* creates the empty variable, as it does without by
foreach spec in "mean(y), by(g)" "rank(y), by(s)" "median(y), by(g s)" {
    capture drop c_out
    capture cegen c_out = `spec'
    local rc = _rc
    capture confirm float variable c_out
    if `rc' == 0 & _rc == 0 & _N == 0 test_pass "cegen_extra empty dataset: `spec'"
    else test_fail "cegen_extra empty dataset: `spec'" "rc `rc'"
}
capture drop c_out

* ---------------------------------------------------------------- many groups
_cegen_xdata 120000 20261128
_cegen_xcmp : NEW = group(key)
_cegen_xcmp : NEW = group(key), label
_cegen_xcmp : int NEW = group(key)
_cegen_xcmp : NEW = tag(key s)
_cegen_xcmp : double NEW = median(y), by(key)
_cegen_xcmp : double NEW = rank(y), by(key)
_cegen_xcmp tol(1e-12) : double NEW = mean(y), by(key)
_cegen_xcmp : double NEW = sd(y), by(id)
_cegen_xcmp : NEW = mode(z), by(key) maxmode
_cegen_xcmp : NEW = count(y), by(key su)
_cegen_xcmp : NEW = group(su key s), missing
_cegen_xcmp : NEW = concat(key s), punct(_)
_cegen_xcmp : int NEW = rank(id)
_cegen_xcmp : NEW = seq(), by(key)
_cegen_xcmp : NEW = cut(y), group(50) label

* ---------------------------------------------------------------- order of the observations
* Order statistics, groups and counts do not depend on it
set seed 20261201
foreach spec in "median(y), by(key)" "rank(y), field by(key)" "pctile(y), p(10) by(g)" ///
    "group(s key)" "count(y), by(su)" "min(y), by(s g)" "mode(z), minmode by(key)" ///
    "iqr(y), by(g h)" "rowmedian(a1 a2 a3)" "cut(y), group(9)" "mad(y), by(sl)" {
    capture drop o1
    capture drop o2
    capture drop u
    sort id
    cegen double o1 = `spec'
    gen double u = runiform()
    sort u
    cegen double o2 = `spec'
    sort id
    capture assert o1 == o2
    if _rc == 0 test_pass "cegen_extra order of observations: `spec'"
    else test_fail "cegen_extra order of observations: `spec'" "results differ"
    drop o1 o2 u
}
}
print_summary "cegen_extra"
global CTOOLS_COMPONENT_COMPLETE "cegen_extra"
    global CTOOLS_CASE_COMPLETE "cegen_extra"
    exit
}

di as error "Unknown cegen case: `__case'"
exit 198
