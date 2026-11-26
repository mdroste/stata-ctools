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

di as error "Unknown cegen case: `__case'"
exit 198
