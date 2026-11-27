* ccollapse against collapse: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
if "`__case'" == "" local __case ccollapse
if "`__case'" == "ccollapse" {
do "validation/support/validate_setup.do"

* every dataset below is drawn from this seed
set seed 20261126

* The collapsed dataset as describe shows it, with every characteristic
capture program drop _cc_meta
program define _cc_meta, rclass
    local sig `"N=`=_N' k=`=c(k)' sort=<`: sortedby'> label=<`: data label'> dta=<`: char _dta[]'>"'
    foreach v of varlist _all {
        local sig `"`sig' | `v' `: type `v'' `: format `v'' <`: variable label `v''> <`: value label `v''> <`: char `v'[]'>"'
    }
    return local sig `"`sig'"'
end

* _cc_pair "name" tol : clist [if] [in] [weight] [, options]
* runs collapse and ccollapse on the data in memory and compares the return
* codes and, after success, the two datasets in full: observations, order,
* names, storage types, formats, labels, characteristics, sort marker, values.
* tol 0 asks for identical values; tol > 0 allows that relative difference in
* numeric values (at least float rounding for float variables). Standard
* deviations and standard errors of the mean are compared to 1e-13 unless
* global CC_EXACT_SD is set: Stata/MP sums the squared deviations of a group
* in as many pieces as it has processors, so collapse's own last bits change
* with set processors. Extra words before the colon are options of ccollapse
* only.
capture program drop _cc_pair
program define _cc_pair
    gettoken name 0 : 0
    gettoken tol 0 : 0
    gettoken extra 0 : 0, parse(":")
    if `"`extra'"' == ":" local extra
    else gettoken colon 0 : 0, parse(":")
    local cmd `"`0'"'
    local why
    preserve
    capture collapse `cmd'
    local rc_ref = _rc
    if !`rc_ref' {
        _cc_meta
        local meta_ref `"`r(sig)'"'
        local i 0
        foreach v of varlist _all {
            local ++i
            rename `v' __r`i'
        }
        tempfile ref
        quietly save `"`ref'"'
    }
    restore, preserve
    if `"`extra'"' != "" & strpos(`"`cmd'"', ",") == 0 local extra `", `extra'"'
    capture ccollapse `cmd' `extra'
    local rc_new = _rc
    if `rc_ref' != `rc_new' local why "return codes `rc_ref' (collapse) and `rc_new' (ccollapse)"
    else if !`rc_ref' {
        _cc_meta
        if `"`r(sig)'"' != `"`meta_ref'"' {
            local why `"description: `r(sig)' instead of `meta_ref'"'
        }
        else {
            unab mine : _all
            quietly merge 1:1 _n using `"`ref'"', nogenerate
            local i 0
            foreach v of local mine {
                local ++i
                local t `tol'
                local label : variable label `v'
                if "$CC_EXACT_SD" == "" & regexm(`"`label'"', "^(\(sd\)|\(semean\)|sd of|semean of)") {
                    local t = max(`t', 1e-13)
                }
                if `t' > 0 & "`: type `v''" == "float" local t = max(`t', 2e-7)
                capture confirm string variable `v'
                if !_rc | `t' == 0 quietly count if `v' != __r`i'
                else quietly count if `v' != __r`i' & !(reldif(`v', __r`i') <= `t')
                if r(N) & `"`why'"' == "" local why "`r(N)' values of `v' differ"
            }
        }
    }
    restore
    if `"`why'"' == "" test_pass `"`name'"'
    else test_fail `"`name'"' `"`why'"'
end

* What a command prints, and its return code. Must run with output on.
capture program drop _cc_output
program define _cc_output, rclass
    _ct_loud _cc_output
    tempfile log
    tempname handle
    quietly log using `"`log'"', text replace name(__cc_output) nomsg
    capture noisily `0'
    local rc = _rc
    quietly log close __cc_output
    file open `handle' using `"`log'"', read text
    file read `handle' line
    local text
    while r(eof) == 0 {
        local text `"`text';`macval(line)'"'
        file read `handle' line
    }
    file close `handle'
    return local text `"`text'"'
    return scalar rc = `rc'
end

* _cc_error "name" : arguments   both commands fail alike: return code, the
* message, and the data left as they were
capture program drop _cc_error
program define _cc_error
    _ct_loud _cc_error
    gettoken name 0 : 0
    gettoken colon 0 : 0, parse(":")
    local cmd `"`0'"'
    preserve
    noisily _cc_output collapse `cmd'
    local rc_ref = r(rc)
    local text_ref `"`r(text)'"'
    restore, preserve
    quietly datasignature
    local before `r(datasignature)'
    noisily _cc_output ccollapse `cmd'
    local rc_new = r(rc)
    local text_new `"`r(text)'"'
    quietly datasignature
    local after `r(datasignature)'
    restore
    if `rc_ref' == 0 | `rc_ref' != `rc_new' {
        test_fail `"`name'"' "return codes `rc_ref' (collapse) and `rc_new' (ccollapse)"
    }
    else if `"`text_ref'"' != `"`text_new'"' {
        test_fail `"`name'"' `"message `text_new' instead of `text_ref'"'
    }
    else if "`before'" != "`after'" test_fail `"`name'"' "the data changed"
    else test_pass `"`name'"'
end

* Mixed data: skewed integer groups, a second key with missing values, a
* string key with empty values, numbers with ordinary and extended missing
* values, values with exact sums, labels, formats and characteristics.
capture program drop _cc_data
program define _cc_data
    args n
    clear
    set obs `n'
    gen long id = _n
    gen int g = ceil(40 * runiform()^2)
    gen byte g2 = floor(4 * runiform())
    replace g2 = . if runiform() < .04
    replace g2 = .b if runiform() < .02
    gen str3 s = char(97 + floor(5 * runiform()))
    replace s = "" if runiform() < .05
    replace s = s + "x" if runiform() < .3
    gen double x = cond(runiform() < .3, rnormal() * 1e6, cond(runiform() < .5, runiform(), rnormal() * 1e-3))
    replace x = . if runiform() < .05
    replace x = .a if runiform() < .02
    gen double r = floor(64 * runiform()) / 8 - 2
    replace r = . if runiform() < .03
    gen byte b = runiform() < .3
    replace b = . if runiform() < .03
    gen int k = rpoisson(3)
    gen float f = runiform() * 100
    gen long big = floor(runiform() * 1e9)
    gen byte wi = floor(6 * runiform())
    replace wi = . if runiform() < .02
    gen double wr = 3 * runiform()
    replace wr = 0 if runiform() < .02
    replace wr = . if runiform() < .01
    gen str8 name = "n" + string(floor(1000 * runiform()))
    replace name = "" if runiform() < .2
    gen int day = 20000 + floor(365 * runiform())
    format day %td
    format x %12.3f
    label define yesno 0 "no" 1 "yes"
    label values b yesno
    label define keys 0 "zero" 1 "one" 2 "two" 3 "three"
    label values g2 keys
    label variable x "a measurement"
    label variable g "group"
    char x[unit] "kg"
    char g[source] "survey"
    note g: a note on g
    note x: a note on x
    label data "ccollapse validation data"
    note: a note on the data
    char _dta[origin] "generated"
    * an all-missing group, and a group with one observation
    replace x = cond(mod(_n, 2), ., .c) if g == 7
    replace r = . if g == 7
    replace b = . if g == 7
    quietly count if g == 40
    if r(N) == 0 replace g = 40 in 1
    quietly replace g = 39 if g == 40 & _n > 1
end


* -----------------------------------------------------------------------
print_section "Statistics, data sorted by the by() variables"
* -----------------------------------------------------------------------
* On sorted data collapse sums each group in the order of the observations,
* as ccollapse does on any data: every statistic must be identical.
_cc_data 3000
sort g g2 id
local allstats mean median p1 p10 p25 p33 p75 p90 p99 sd semean sepoisson sum rawsum count percent max min iqr first last firstnm lastnm
foreach st of local allstats {
    * collapse's sepoisson is not a number once a group total reaches 2^31
    local vars x r k f big
    if "`st'" == "sepoisson" local vars x r k f
    _cc_pair "(`st') by(g)" 0 : (`st') `vars', by(g)
}
_cc_pair "(sebinomial) by(g)" 0 : (sebinomial) b (sebinomial) sk = k (sebinomial) sr = r, by(g)
_cc_pair "sepoisson with negative values" 0 : (sepoisson) x r, by(g)
_cc_pair "all statistics of one variable" 0 : (mean) m=x (median) md=x (p5) p5=x (sd) sd=x (semean) se=x (sum) s=x (rawsum) rs=x (count) n=x (percent) pc=x (max) mx=x (min) mn=x (iqr) iq=x (first) fi=x (last) la=x (firstnm) fn=x (lastnm) ln=x, by(g)
_cc_pair "two keys" 0 : (mean) x r (sd) sx=x (p50) k (max) f (count) n=x, by(g g2)
_cc_pair "no by()" 0 : (mean) x r (sd) sx=x (p50) k (max) f (count) n=x (sum) sr=r (first) id (last) lid=id
_cc_pair "default statistic" 0 : x r k
_cc_pair "default statistic, by" 0 : x r k, by(g)
_cc_pair "statistic carried over a target" 0 : (sum) a=x b2=r k (max) x
_cc_pair "target = source without spaces" 0 : (sum) a=x (mean) c =r d= k (max) e = f, by(g)
_cc_pair "one variable under several statistics" 0 : (mean) x (sd) sx=x (min) mnx=x (mean) mr=r (sum) r, by(g)
_cc_pair "varlist ranges and wildcards" 0 : (mean) x-k (max) w*, by(g)
_cc_pair "abbreviated variable names" 0 : (mean) bi (sum) da, by(g)
_cc_pair "sem, seb and sep abbreviations" 0 : (sem) x (seb) b (sep) k (semea) s2=x (sebinom) b2=b (sepois) k2=k, by(g)
_cc_pair "leading zero in p#" 0 : (p05) x (p5) q=x, by(g)
_cc_pair "integer and long sources" 0 : (mean) k big id (sum) sk=k sb=big (median) mk=k mb=big, by(g)
_cc_pair "float source" 0 : (mean) f (sd) sf=f (sum) tf=f (p25) pf=f (min) mnf=f (max) mxf=f, by(g)

sort s id
_cc_pair "string key" 0 : (mean) x r (sd) sx=x (p50) k (count) n=x (first) id, by(s)
sort s g2 id
_cc_pair "string and numeric keys" 0 : (mean) x r (sd) sx=x (p50) k (count) n=x, by(s g2)
_cc_pair "keys in another order than the data" 0 : (sum) r k (count) n=x, by(g2 s)
sort day id
_cc_pair "formatted key" 0 : (mean) x (max) mid=id, by(day)

* with one processor collapse sums the squared deviations in order too
local processors = c(processors)
set processors 1
global CC_EXACT_SD 1
sort g g2 id
_cc_pair "one processor: sd and semean" 0 : (sd) x r k f (semean) sx=x sr=r sk=k sf=f, by(g)
_cc_pair "one processor: semean without sd" 0 : (semean) x r k f (mean) mx=x, by(g)
_cc_pair "one processor: sd, two keys" 0 : (sd) x r (semean) k f, by(g g2)
_cc_pair "one processor: sd, no by()" 0 : (sd) x r (semean) k f (semean) sx=x
_cc_pair "one processor: sd [fw]" 0 : (sd) x r k f (semean) sx=x sr=r [fw=wi], by(g)
_cc_pair "one processor: semean without sd [fw]" 0 : (semean) x r k f [fw=wi], by(g)
_cc_pair "one processor: sd [iw]" 0 : (sd) x r k f [iw=wr], by(g)
global CC_EXACT_SD
set processors `processors'

* -----------------------------------------------------------------------
print_section "Weights"
* -----------------------------------------------------------------------
* Weighted standard deviations may differ from collapse's in the last bit
_cc_data 3000
sort g g2 id
local plain mean median p10 p25 p75 p90 sum rawsum count percent max min iqr first last firstnm lastnm
foreach wt in fw aw iw pw {
    foreach st of local plain {
        _cc_pair "(`st') [`wt'=wi]" 0 : (`st') x r k [`wt'=wi], by(g)
    }
    if "`wt'" != "fw" {
        foreach st of local plain {
            _cc_pair "(`st') [`wt'=wr]" 0 : (`st') x r k [`wt'=wr], by(g)
        }
    }
}
foreach wt in fw aw iw {
    _cc_pair "(sd) [`wt'=wi]" 1e-14 : (sd) x r k f [`wt'=wi], by(g)
}
_cc_pair "(sd) [aw=wr]" 1e-14 : (sd) x r k f [aw=wr], by(g)
_cc_pair "(sd) [iw=wr]" 1e-14 : (sd) x r k f [iw=wr], by(g)
_cc_pair "(semean) [fw=wi]" 1e-14 : (semean) x r k f [fw=wi], by(g)
_cc_pair "(semean) [aw=wr]" 1e-14 : (semean) x r k f [aw=wr], by(g)
_cc_pair "(sebinomial) (sepoisson) [fw=wi]" 0 : (sebinomial) b (sepoisson) k r, by(g)
_cc_pair "weight expression" 0 : (mean) x (sum) r (count) n=k [aw=wr*2+1], by(g)
_cc_pair "weight expression with missing values" 0 : (mean) x (sum) r (count) n=k [pw=wr/wi], by(g)
_cc_pair "weight that is a source" 0 : (mean) wr (sum) x [aw=wr], by(g)
_cc_pair "negative iweights" 1e-14 : (mean) x (sum) sx=x (count) n=x (sd) sd=x (p50) md=x [iw=r], by(g)
_cc_pair "weights, no by()" 0 : (mean) x r (sum) sx=x (count) n=x (p50) md=x (percent) pc=x [aw=wr]
_cc_pair "weights, if and in" 0 : (mean) x r (sum) sx=x (count) n=x (p50) md=x [fw=wi] if f > 20 in 100/2800, by(g)
_cc_pair "weights, two keys" 0 : (mean) x r (sum) sx=x (count) n=x (p50) md=x [pw=wr], by(g g2)

* -----------------------------------------------------------------------
print_section "Sample: if, in, cw"
* -----------------------------------------------------------------------
_cc_data 3000
sort g g2 id
_cc_pair "if" 0 : (mean) x r (sd) sx=x (count) n=x if f > 50, by(g)
_cc_pair "in" 0 : (mean) x r (sd) sx=x (count) n=x in 500/2500, by(g)
_cc_pair "if and in" 0 : (mean) x r (sd) sx=x (count) n=x if f > 50 in 500/2500, by(g)
_cc_pair "in, no by()" 0 : (mean) x r (sd) sx=x (count) n=x in 2/40
_cc_pair "if selecting one observation" 0 : (mean) x r (sd) sx=x (count) n=x (p50) md=r if id == 1234, by(g)
_cc_pair "if on a key" 0 : (mean) x r (count) n=x if g2 >= 2, by(g g2)
_cc_pair "if with missing values in the keys" 0 : (sum) r k (count) n=x if g2 >= ., by(g2)
_cc_pair "cw" 0 : (mean) x r b (count) nx=x nr=r nb=b, by(g) cw
_cc_pair "cw, if, weights" 0 : (mean) x r b (count) nx=x nr=r nb=b [aw=wr] if f > 10, by(g) cw
_cc_pair "cw, no by()" 0 : (mean) x r b (sum) sx=x cw
_cc_pair "fast" 0 : (mean) x r, by(g) fast
_cc_pair "clabel" 0 : (mean) x (p50) md=r (sum) k, by(g) clabel
_cc_pair "frame()" 0 : (mean) x, by(g) frame(somewhere)

* -----------------------------------------------------------------------
print_section "Strings, missing values, storage types"
* -----------------------------------------------------------------------
_cc_data 3000
sort g g2 id
_cc_pair "string first and last" 0 : (first) name s (last) lname=name (firstnm) fname=name (lastnm) lnname=name, by(g)
_cc_pair "string statistics, no by()" 0 : (first) name (last) lname=name (firstnm) fname=name (lastnm) lnname=name
_cc_pair "strings next to numbers, weights" 0 : (first) name (mean) x (lastnm) ln=name (max) r [aw=wr], by(g g2)
_cc_pair "storage types" 0 : (max) b k big f x day (min) nb=b nk=k nbig=big nf=f nx=x (first) fb=b fk=k (last) lbig=big (firstnm) fnf=f (lastnm) lnx=x, by(g)
_cc_pair "labelled sources" 0 : (mean) b (max) mb=b (first) fb=b (count) nb=b, by(g2)
sort g2 id
_cc_pair "missing values in the key" 0 : (mean) x (count) n=x, by(g2)
clear
input g x str3 s
1 .a ""
1 .b ""
2 .b "q"
2 3  ""
2 .  "z"
3 .  ""
4 .z "a"
4 .  "b"
4 .c ""
5 -0 "m"
5 0  "m"
end
_cc_pair "groups of missing values" 0 : (first) f=x (firstnm) fn=x (last) l=x (lastnm) ln=x (max) mx=x (min) mn=x (first) fs=s (firstnm) fns=s (last) ls=s (lastnm) lns=s (count) c=x (sum) sm=x (mean) m=x (sd) sd=x (median) md=x (iqr) i=x (semean) se=x (percent) pc=x (rawsum) rs=x, by(g)
_cc_pair "all values missing" 0 : (mean) x (sum) s2=x (count) n=x (percent) pc=x (max) mx=x (min) mn=x (p50) md=x if g == 1 | g == 3, by(g)
clear
input g str3 s x
.a "b" 1
.  "a" 2
1  ""  3
.b "B" 4
1  "a" 5
-1 "a" 6
.  "a" 7
0  "Z" 8
-0 "Z" 9
end
_cc_pair "order of missing and string keys" 0 : (sum) x, by(g s)
_cc_pair "order of keys, string first" 0 : (sum) x, by(s g)
set type double
_cc_data 500
sort g id
_cc_pair "set type double" 0 : (mean) k b f (sd) sk=k (p50) mf=f, by(g)
set type float

* -----------------------------------------------------------------------
print_section "Data in any order"
* -----------------------------------------------------------------------
* collapse's sums depend on its own order of the observations within a
* group when the data are not sorted. Variables whose partial sums are exact
* (whole numbers, eighths) give identical results in any order; other sums
* and every standard deviation agree to rounding.
_cc_data 3000
gen double u = runiform()
sort u
foreach st in mean median p10 p90 sum rawsum count percent max min iqr {
    _cc_pair "unsorted (`st')" 0 : (`st') r k big, by(g)
}
_cc_pair "unsorted (sepoisson)" 0 : (sepoisson) r k, by(g)
_cc_pair "unsorted sums of inexact values" 1e-9 : (mean) f (sum) sf=f (sd) sdf=f sdr=r (semean) sef=f, by(g)
_cc_pair "unsorted first and last" 0 : (first) id name (last) lid=id lname=name (firstnm) fx=x (lastnm) lx=x, by(g)
_cc_pair "unsorted, two keys" 0 : (mean) r k (max) x (p50) md=x (count) n=x, by(g g2)
_cc_pair "unsorted, string key" 0 : (mean) r k (max) x (p50) md=x (count) n=x, by(s)
_cc_pair "unsorted, string and numeric keys" 0 : (mean) r k (max) x (p50) md=x (count) n=x (first) id, by(s g2)
_cc_pair "unsorted, weights" 0 : (mean) r k (sum) sr=r (count) n=x (max) x (p50) md=x [fw=wi], by(g)
_cc_pair "unsorted, real weights" 1e-12 : (mean) r k (sum) sr=r (count) n=x (percent) pc=x [aw=wr], by(g g2)
_cc_pair "unsorted, if in cw" 0 : (mean) r k (count) nr=r nx=x if f > 30 in 20/2900, by(g2 g) cw
_cc_pair "unsorted, non-integer key" 0 : (sum) k (count) n=x (max) r, by(f)
_cc_pair "unsorted, wide integer key" 0 : (sum) k (count) n=x (max) r, by(big)
_cc_pair "unsorted, by the observation number" 0 : (sum) k x (first) name, by(id)

* -----------------------------------------------------------------------
print_section "Large data: every grouping method, in parallel"
* -----------------------------------------------------------------------
clear
set obs 400000
gen long id = _n
gen long firm = ceil(_n / 37)
gen int year = 1990 + mod(_n, 30)
gen double x = rnormal()
replace x = . if runiform() < .02
gen double r = floor(64 * runiform()) / 8
gen int k = rpoisson(4)
gen byte wi = 1 + floor(4 * runiform())
gen double wr = runiform() + .1
gen str6 code = "c" + string(floor(300 * runiform()))
gen double odd = floor(runiform() * 1e6) / 7
gen long sparse = floor(runiform() * 2000) * 1000003
local big (mean) x r k (sd) sx=x (sum) tr=r tx=x (count) n=x (p50) md=x (p90) p90=r (iqr) iq=x (min) mnx=x (max) mxx=x (first) fid=id (last) lid=id (firstnm) fx=x (lastnm) lx=x (percent) pc=x

_cc_pair "ordered rows: one key" 0 : `big', by(firm)
if scalar(_ccollapse_method) == 2 test_pass "ordered rows are recognized"
else test_fail "ordered rows are recognized" "method `=scalar(_ccollapse_method)'"
_cc_pair "ordered rows: fweights" 1e-14 : `big' [fw=wi], by(firm)
_cc_pair "ordered rows: aweights" 1e-14 : `big' [aw=wr], by(firm)
_cc_pair "ordered rows: if" 0 : `big' if k > 2, by(firm)
_cc_pair "one group" 0 : `big'
if scalar(_ccollapse_method) == 1 test_pass "no by() is one group"
else test_fail "no by() is one group" "method `=scalar(_ccollapse_method)'"
_cc_pair "one group, weights" 1e-14 : `big' [aw=wr]

* the same on shuffled data, statistics that do not depend on the order
gen double u = runiform()
sort u
local anyorder (mean) r k (sum) tr=r tk=k (count) n=x (p50) md=x (p90) p90=r (iqr) iq=x (min) mnx=x (max) mxx=x (percent) pc=x (first) fid=id (last) lid=id (firstnm) fx=x (lastnm) lx=x
_cc_pair "dense keys: one integer key" 0 : `anyorder', by(firm)
if scalar(_ccollapse_method) == 3 test_pass "dense integer keys are recognized"
else test_fail "dense integer keys are recognized" "method `=scalar(_ccollapse_method)'"
_cc_pair "dense keys: two integer keys" 0 : `anyorder', by(year wi)
_cc_pair "dense keys: weights" 0 : `anyorder' [fw=wi], by(year)
_cc_pair "hashed keys: strings" 0 : `anyorder', by(code)
if scalar(_ccollapse_method) == 4 test_pass "string keys are hashed"
else test_fail "string keys are hashed" "method `=scalar(_ccollapse_method)'"
_cc_pair "hashed keys: string and number" 0 : `anyorder', by(code year)
_cc_pair "hashed keys: sparse integers" 0 : `anyorder', by(sparse)
_cc_pair "hashed keys: weights and if" 0 : `anyorder' [fw=wi] if k > 1, by(code)
_cc_pair "sorted keys: many distinct values" 0 : `anyorder', by(odd)
if scalar(_ccollapse_method) == 5 test_pass "many distinct keys are sorted"
else test_fail "many distinct keys are sorted" "method `=scalar(_ccollapse_method)'"
_cc_pair "sorted keys: two keys" 0 : `anyorder', by(odd code)
_cc_pair "inexact sums in any order" 1e-9 : (mean) x (sum) tx=x (sd) sx=x sr=r, by(year)

* results do not depend on the number of threads
foreach keys in firm code odd "year wi" "" {
    local opt
    if "`keys'" != "" local opt by(`keys')
    preserve
    ccollapse `big' [aw=wr], `opt' threads(1)
    tempfile one
    quietly save `"`one'"'
    restore, preserve
    ccollapse `big' [aw=wr], `opt' threads(7)
    capture cf _all using `"`one'"'
    if _rc == 0 test_pass "threads(1) and threads(7) agree, by(`keys')"
    else test_fail "threads(1) and threads(7) agree, by(`keys')" "cf returned `=_rc'"
    restore
}

* groups of a million rows and more: percentiles selected by all threads
clear
set obs 2600000
gen byte half = _n > 1400000
gen byte lopsided = _n > 2590000
gen double v = rnormal()
replace v = . if runiform() < .01
gen int tied = floor(20 * runiform())
gen double ramp = _n / 8
gen float fl = runiform()
local pcts (p1) a1=v (p25) a25=v (median) a50=v (p75) a75=v (p99) a99=v (iqr) aiq=v (p50) t50=tied (p10) t10=tied (p90) t90=tied (iqr) tiq=tied (p33) r33=ramp (p50) r50=ramp (p50) f50=fl (count) n=v
_cc_pair "large group: no by()" 0 : `pcts'
_cc_pair "large groups: two ordered groups" 0 : `pcts', by(half)
_cc_pair "large and small group" 0 : `pcts', by(lopsided)
_cc_pair "large group: if" 0 : `pcts' if fl > .3
_cc_pair "large group: every percentile" 0 : (p1) v (p2) p2=v (p3) p3=v (p5) p5=v (p10) p10=v (p20) p20=v (p30) p30=v (p40) p40=v (p50) p50=v (p60) p60=v (p70) p70=v (p80) p80=v (p90) p90=v (p95) p95=v (p97) p97=v (p98) p98=v (p99) p99=v
gen double u = runiform()
sort u
_cc_pair "large groups: shuffled rows" 0 : `pcts', by(half)
_cc_pair "large group: one thread" 0 threads(1) : `pcts'

* -----------------------------------------------------------------------
print_section "Time-series operators and strL variables"
* -----------------------------------------------------------------------
clear
set obs 600
gen int panel = ceil(_n / 20)
bysort panel: gen int t = _n
gen double y = rnormal()
gen int z = rpoisson(5)
format y %10.2f
quietly xtset panel t
_cc_pair "time-series operators" 0 : (mean) L.y D.z (sum) F2.y (max) y, by(panel)
_cc_pair "time-series operator as a source" 0 : (mean) ly = L.y (sd) dz = D.z (count) n = L2.y, by(panel)
_cc_pair "time-series operators, no by()" 0 : (mean) L.y l2.z

_cc_data 400
sort g id
gen strL long1 = name + "-" + string(id)
gen strL key = s
_cc_pair "strL source" 0 : (first) long1 (last) l=long1 (mean) x, by(g)
_cc_pair "strL key" 0 : (mean) x r (count) n=x, by(key)
_cc_pair "strL key and source" 0 : (firstnm) long1 (sum) r, by(key g2)

* -----------------------------------------------------------------------
print_section "Errors"
* -----------------------------------------------------------------------
_cc_data 300
sort g id
foreach bad in "(bogus) x" "(p0) x" "(p100) x" "(p50.5) x" "(p) x" "(MEAN) x" "(se) x" "(mean median) x" "() x" ///
    "" "(mean)" "(mean) x (sd)" "(mean) nosuchvar" "(mean) x, by(nosuchvar)" "(mean) 1x = x" "(mean) a = nosuchvar" ///
    "(mean) x (sum) x, by(g)" "(mean) a=x (p50) a=r" "(mean) x x" "(mean) g, by(g)" "(mean) q=x (sum) g2=r, by(g g2)" ///
    "(mean) name, by(g)" "(p50) name" "(count) name, by(g)" "(sd) x [pw=wr], by(g)" "(semean) x [pw=wr]" ///
    "(sebinomial) b [pw=wr]" "(sepoisson) k [pw=wr]" "(semean) x [iw=wr]" "(sebinomial) b [iw=wr]" ///
    "(sepoisson) k [iw=wr]" "(sebinomial) b [aw=wr]" "(sepoisson) k [aw=wr]" "(mean) x [fw=wr], by(g)" ///
    "(mean) x [aw=r], by(g)" "(mean) x [fw=r], by(g)" "(mean) x [pw=r]" "(mean) x if g > 1000, by(g)" ///
    "(mean) x in 1/5 if g > 1000" "(first) name (mean) x, by(g) cw" "(mean) x [aw=name]" "(mean) x, bogus" ///
    "(mean) x if" "(mean) x [aw=wr] if g > 1000, by(g)" "(mean) x (sum) x if g > 1000" "(sd) x [pw=wr] if g > 1000" ///
    "(mean) name if g > 1000" "(mean) x (sum) x [aw=r]" "(first) name, cw" "(mean) x [iw=wr*0]" {
    noisily _cc_error `"error: collapse `bad'"' : `bad'
}
clear
noisily _cc_error "error: no data" : (mean) x
set obs 5
gen x = _n
noisily _cc_error "error: no variable" : (mean) y
capture ccollapse (mean) x, threads(-1)
if _rc == 198 test_pass "threads(-1) is rejected"
else test_fail "threads(-1) is rejected" "rc `=_rc'"
capture ccollapse (mean) x, threads(300)
if _rc == 198 test_pass "threads(300) is rejected"
else test_fail "threads(300) is rejected" "rc `=_rc'"

* -----------------------------------------------------------------------
print_section "ccollapse alone"
* -----------------------------------------------------------------------
* The note on aweights and counts, as collapse prints it
_cc_data 300
preserve
noisily _cc_output collapse (count) x (mean) r [aw=wr], by(g)
local text_ref `"`r(text)'"'
restore, preserve
noisily _cc_output ccollapse (count) x (mean) r [aw=wr], by(g)
if `"`r(text)'"' == `"`text_ref'"' & strpos(`"`text_ref'"', "not used to compute") test_pass "note on aweights and counts"
else test_fail "note on aweights and counts" `"`r(text)' instead of `text_ref'"'
restore
* verbose reports and changes nothing else
preserve
ccollapse (mean) x r (p50) k, by(g)
tempfile quiet
quietly save `"`quiet'"'
restore, preserve
noisily _cc_output ccollapse (mean) x r (p50) k, by(g) verbose threads(2)
local text `"`r(text)'"'
capture cf _all using `"`quiet'"'
if _rc == 0 & strpos(`"`text'"', "ccollapse timing breakdown") & scalar(_ccollapse_threads) == 2 test_pass "verbose and threads()"
else test_fail "verbose and threads()" `"cf `=_rc', threads `=scalar(_ccollapse_threads)'"'
restore
* stored results
preserve
quietly count if !missing(wr) & wr != 0 & f > 10
local expected = r(N)
ccollapse (mean) x [aw=wr] if f > 10, by(g)
if scalar(_ccollapse_nobs) == `expected' & scalar(_ccollapse_ngroups) == _N test_pass "stored results"
else test_fail "stored results" "nobs `=scalar(_ccollapse_nobs)' ngroups `=scalar(_ccollapse_ngroups)'"
restore
* sepoisson where the group total reaches 2^31: collapse stores a value
* that is not a number; ccollapse stores sqrt(total)/n
preserve
collapse (sum) total = big (count) n = big, by(g)
tempfile totals
quietly save `"`totals'"'
restore, preserve
ccollapse (sepoisson) se = big, by(g)
quietly merge 1:1 g using `"`totals'"', nogenerate
quietly count if total >= 2^31 & n > 0
local large = r(N)
quietly count if reldif(se, sqrt(total) / n) > 1e-6 & n > 0
if r(N) == 0 & `large' > 0 test_pass "sepoisson of large totals"
else test_fail "sepoisson of large totals" "`r(N)' of `large' groups differ"
restore
* a by() variable may also be a source (collapse stops with r(198))
preserve
quietly generate long gcopy = g
collapse (mean) m = gcopy (count) n = gcopy (max) x, by(g)
tempfile ref
quietly save `"`ref'"'
restore, preserve
capture noisily ccollapse (mean) m = g (count) n = g (max) x, by(g)
if _rc == 0 {
    label variable m "(mean) gcopy"
    label variable n "(count) gcopy"
    format m n %8.0g
    capture cf _all using `"`ref'"'
}
if _rc == 0 test_pass "a by() variable as a source"
else test_fail "a by() variable as a source" "rc `=_rc'"
restore
* S_FN and S_FNDATE are cleared, as by collapse
preserve
global S_FN "somefile.dta"
global S_FNDATE "1 Jan 2026"
ccollapse (mean) x, by(g)
if "$S_FN" == "" & "$S_FNDATE" == "" test_pass "S_FN and S_FNDATE cleared"
else test_fail "S_FN and S_FNDATE cleared" "$S_FN"
restore
* nothing is left behind besides the documented scalars
preserve
ccollapse (mean) x (p50) r [aw=wr] if f > 10, by(g s)
capture confirm variable __000000
local leftovers : all globals "__ct*"
if "`leftovers'" == "" test_pass "no temporary globals"
else test_fail "no temporary globals" "`leftovers'"
restore
* Break and errors before the dataset is replaced leave the data alone
preserve
quietly datasignature
local before `r(datasignature)'
capture ccollapse (mean) x [fw=wr], by(g)
local rc = _rc
quietly datasignature
if `rc' == 401 & "`r(datasignature)'" == "`before'" test_pass "an error leaves the data as they were"
else test_fail "an error leaves the data as they were" "rc `rc'"
restore
* version control: under version 4 collapse's old command runs
preserve
capture noisily version 4: ccollapse x, by(g)
local rc_new = _rc
restore, preserve
capture noisily version 4: collapse x, by(g)
if `rc_new' == _rc test_pass "version 4 goes to collapse"
else test_fail "version 4 goes to collapse" "rc `rc_new' and `=_rc'"
restore

noisily print_summary "ccollapse"
global CTOOLS_CASE_COMPLETE "ccollapse"
}
* -----------------------------------------------------------------------
* ccollapse_extra: weights with every statistic, percentile ties, extreme
* values and type limits, keys that stress each grouping method, metadata,
* clist parsing, error order, large data, threads
* -----------------------------------------------------------------------
if "`__case'" == "ccollapse_extra" {
do "validation/support/validate_setup.do"

* every dataset below is drawn from this seed
set seed 20261127

* The collapsed dataset as describe shows it, with every characteristic
capture program drop _cx_meta
program define _cx_meta, rclass
    local sig `"N=`=_N' k=`=c(k)' sort=<`: sortedby'> label=<`: data label'> dta=<`: char _dta[]'>"'
    foreach v of varlist _all {
        local sig `"`sig' | `v' `: type `v'' `: format `v'' <`: variable label `v''> <`: value label `v''> <`: char `v'[]'>"'
    }
    return local sig `"`sig'"'
end

* _cx_pair "name" tol [ccollapse option] : clist [if] [in] [weight] [, options]
* collapse and ccollapse on the data in memory: the same return code and,
* after success, the same dataset (description and values). tol as in the
* ccollapse case: 0 asks for identical values, sd and semean get 1e-13.
capture program drop _cx_pair
program define _cx_pair
    gettoken name 0 : 0
    gettoken tol 0 : 0
    gettoken extra 0 : 0, parse(":")
    if `"`extra'"' == ":" local extra
    else gettoken colon 0 : 0, parse(":")
    local cmd `"`0'"'
    local why
    preserve
    capture collapse `cmd'
    local rc_ref = _rc
    if !`rc_ref' {
        _cx_meta
        local meta_ref `"`r(sig)'"'
        local i 0
        foreach v of varlist _all {
            local ++i
            rename `v' __r`i'
        }
        tempfile ref
        quietly save `"`ref'"'
    }
    restore, preserve
    if `"`extra'"' != "" & strpos(`"`cmd'"', ",") == 0 local extra `", `extra'"'
    capture ccollapse `cmd' `extra'
    local rc_new = _rc
    if `rc_ref' != `rc_new' local why "return codes `rc_ref' (collapse) and `rc_new' (ccollapse)"
    else if !`rc_ref' {
        _cx_meta
        if `"`r(sig)'"' != `"`meta_ref'"' {
            local why `"description: `r(sig)' instead of `meta_ref'"'
        }
        else {
            unab mine : _all
            quietly merge 1:1 _n using `"`ref'"', nogenerate
            local i 0
            foreach v of local mine {
                local ++i
                local t `tol'
                local label : variable label `v'
                if regexm(`"`label'"', "^(\(sd\)|\(semean\)|sd of|semean of)") local t = max(`t', 1e-13)
                if `t' > 0 & "`: type `v''" == "float" local t = max(`t', 2e-7)
                capture confirm string variable `v'
                if !_rc | `t' == 0 quietly count if `v' != __r`i'
                else quietly count if `v' != __r`i' & !(reldif(`v', __r`i') <= `t')
                if r(N) & `"`why'"' == "" local why "`r(N)' values of `v' differ"
            }
        }
    }
    restore
    if `"`why'"' == "" test_pass `"`name'"'
    else test_fail `"`name'"' `"`why'"'
end

* What a command prints, and its return code. Must run with output on.
capture program drop _cx_output
program define _cx_output, rclass
    _ct_loud _cx_output
    tempfile log
    tempname handle
    quietly log using `"`log'"', text replace name(__cx_output) nomsg
    capture noisily `0'
    local rc = _rc
    quietly log close __cx_output
    file open `handle' using `"`log'"', read text
    file read `handle' line
    local text
    while r(eof) == 0 {
        local text `"`text';`macval(line)'"'
        file read `handle' line
    }
    file close `handle'
    return local text `"`text'"'
    return scalar rc = `rc'
end

* _cx_error "name" : arguments   both commands fail alike: return code, the
* message, and the data left as they were
capture program drop _cx_error
program define _cx_error
    _ct_loud _cx_error
    gettoken name 0 : 0
    gettoken colon 0 : 0, parse(":")
    local cmd `"`0'"'
    preserve
    noisily _cx_output collapse `cmd'
    local rc_ref = r(rc)
    local text_ref `"`r(text)'"'
    restore, preserve
    quietly datasignature
    local before `r(datasignature)'
    noisily _cx_output ccollapse `cmd'
    local rc_new = r(rc)
    local text_new `"`r(text)'"'
    quietly datasignature
    local after `r(datasignature)'
    restore
    if `rc_ref' == 0 | `rc_ref' != `rc_new' {
        test_fail `"`name'"' "return codes `rc_ref' (collapse) and `rc_new' (ccollapse)"
    }
    else if `"`text_ref'"' != `"`text_new'"' {
        test_fail `"`name'"' `"message `text_new' instead of `text_ref'"'
    }
    else if "`before'" != "`after'" test_fail `"`name'"' "the data changed"
    else test_pass `"`name'"'
end

* _cx_same "name" "options 1" "options 2" : clist ...   two ccollapse runs
* with different ccollapse options give the same dataset (cf _all)
capture program drop _cx_same
program define _cx_same
    gettoken name 0 : 0
    gettoken o1 0 : 0
    gettoken o2 0 : 0
    gettoken colon 0 : 0, parse(":")
    local cmd `"`0'"'
    if strpos(`"`cmd'"', ",") == 0 {
        local o1 `", `o1'"'
        local o2 `", `o2'"'
    }
    local why
    preserve
    capture ccollapse `cmd' `o1'
    local rc1 = _rc
    tempfile one
    if !`rc1' quietly save `"`one'"'
    restore, preserve
    capture ccollapse `cmd' `o2'
    local rc2 = _rc
    if `rc1' | `rc2' local why "return codes `rc1' and `rc2'"
    else {
        capture cf _all using `"`one'"'
        if _rc local why "cf returned `=_rc'"
    }
    restore
    if `"`why'"' == "" test_pass `"`name'"'
    else test_fail `"`name'"' `"`why'"'
end

* -----------------------------------------------------------------------
print_section "Every weight type with the statistics the Weights section leaves out"
* -----------------------------------------------------------------------
* positive whole weights (fw) and exact binary weights (aw, iw, pw); the
* standard errors collapse refuses for some weights must fail alike
clear
set obs 2000
gen int g = ceil(25 * runiform())
gen double x = rnormal() * 10
replace x = . if runiform() < .05
replace x = .b if runiform() < .02
gen double r = floor(32 * runiform()) / 4 - 3
gen int k = rpoisson(4)
replace k = . if runiform() < .03
gen byte b = runiform() < .4
replace b = . if runiform() < .03
gen float f = runiform() * 10
gen byte wp = 1 + floor(5 * runiform())
gen double wq = .25 + floor(16 * runiform()) / 4
gen str5 name = "n" + string(floor(100 * runiform()))
gen long id = _n
sort g id
foreach wt in fw aw iw pw {
    local w wq
    if "`wt'" == "fw" local w wp
    _cx_pair "(p1) (p33) (p99) [`wt']" 0 : (p1) x r k (p33) x33=x r33=r k33=k (p99) x99=x r99=r f99=f [`wt'=`w'], by(g)
    foreach st in sd semean sebinomial sepoisson {
        _cx_pair "(`st') [`wt']" 0 : (`st') x r k b f [`wt'=`w'], by(g)
    }
}
* a whole group whose weights are all zero drops out
gen double wz = cond(g == 3, 0, wq)
gen byte wzp = cond(g == 3, 0, wp)
foreach wt in fw aw iw pw {
    local w wz
    if "`wt'" == "fw" local w wzp
    _cx_pair "group with zero weights [`wt']" 0 : (mean) x (sum) sr=r (count) n=x (p50) md=r (percent) pc=x (first) id [`wt'=`w'], by(g)
}
_cx_pair "percent with negative iweights" 0 : (percent) x (count) n=x (sum) s=r [iw=r], by(g)
_cx_pair "percent, if, no by()" 0 : (percent) x (count) n=x [fw=wp] if g < 5
_cx_pair "weight that is the by() variable" 0 : (mean) x (count) n=x (sum) s=r (p50) md=r [fw=g], by(g)
_cx_pair "float weight variable" 0 : (mean) x (p50) md=r (sum) s=r (percent) pc=x [aw=f], by(g)
* hundreds of targets
preserve
clear
set obs 500
gen int g = ceil(_n / 50)
forvalues j = 1/600 {
    gen float v`j' = runiform()
}
_cx_pair "600 targets" 0 : (mean) v1-v600 (p50) m1=v1 m600=v600 (sd) s600=v600, by(g)
_cx_pair "600 percentile targets" 0 : (p75) v1-v600, by(g)
restore

* ccollapse's own results do not depend on the incoming order for statistics
* that collapse also computes order-free
preserve
ccollapse (count) x (sum) r (min) mn=x (max) mx=x (p50) md=x (p90) p90=r (iqr) iq=k (percent) pc=x [fw=wp], by(g)
tempfile ordered
quietly save `"`ordered'"'
restore, preserve
gen double u = runiform()
sort u
ccollapse (count) x (sum) r (min) mn=x (max) mx=x (p50) md=x (p90) p90=r (iqr) iq=k (percent) pc=x [fw=wp], by(g)
capture cf _all using `"`ordered'"'
if _rc == 0 test_pass "same result from shuffled data"
else test_fail "same result from shuffled data" "cf returned `=_rc'"
restore

* -----------------------------------------------------------------------
print_section "Error order and error parity"
* -----------------------------------------------------------------------
foreach bad in "(sd) x (mean) name [pw=wq], by(g)" "(mean) name (sd) x [pw=wq], by(g)" ///
    "(mean) x [aw=nosuchvar], by(g)" `"(mean) x [aw="a"], by(g)"' "(mean) x [aw=wq+], by(g)" ///
    "(mean) x in 1/99999, by(g)" "(mean) x [fw=wp] if 0, by(g)" "(mean) x [fw=0*wp]" ///
    "(mean) x, by(g) cw nosuchopt" "(mean) x [aw=wq], by(g x)" "(sebinomial) b [aw=wq] if g > 100" ///
    "(mean) x [aw=-wq], by(g)" "(sd) x (mean) x [pw=wq]" "(sum) x [iw=x*.] if g < 3" ///
    "(mean) x (p50) n2=x (sum) n2=r" "(mean) x [fw=wq], by(g)" {
    noisily _cx_error `"error: collapse `bad'"' : `bad'
}
* empty data with variables
preserve
quietly drop in 1/l
noisily _cx_error "error: zero observations" : (mean) x, by(g)
noisily _cx_error "error: zero observations, bad weight" : (mean) x [aw=nosuchvar], by(g)
noisily _cx_error "error: zero observations, conflict" : (mean) x (sum) x
restore
* cw with a variable missing everywhere
gen double allmiss = .
noisily _cx_error "error: cw with an all-missing source" : (mean) x a=allmiss, by(g) cw
drop allmiss

* -----------------------------------------------------------------------
print_section "Parsing of the clist"
* -----------------------------------------------------------------------
_cx_pair "no space after (stat)" 0 : (mean)x (sum)sr=r, by(g)
_cx_pair "no space before (stat)" 0 : (mean) x(sum) sr=r, by(g)
_cx_pair "spaces inside the parentheses" 0 : ( mean ) x (median ) md=r, by(g)
_cx_pair "p# written oddly" 0 : (p050) x (p5.0) y=r (p1e1) z=k, by(g)
_cx_pair "spaces in by()" 0 : (mean) x (count) n=x, by( g )
_cx_pair "weight spelled out" 0 : (mean) x (sum) sr=r [aweight = wq], by(g)
_cx_pair "targets swap names with sources" 0 : (mean) x=r (mean) r=x (sum) k=x (max) f=k, by(g)
_cx_pair "default statistic after named ones" 0 : x (count) n=x r, by(g)
_cx_pair "in with f and l, if" 0 : (mean) x (count) n=x in f/l if g == 3, by(g)
_cx_pair "in -#/l" 0 : (mean) x (count) n=x in -50/l, by(g)
_cx_pair "in #/# reversed" 0 : (mean) x (count) n=x in 50/10, by(g)

* -----------------------------------------------------------------------
print_section "Percentiles with ties and weights at exact boundaries"
* -----------------------------------------------------------------------
clear
input g x w
1 1 1
1 2 1
1 3 1
1 4 1
2 2 1
2 1 1
2 2 1
2 1 1
3 5 3
3 5 1
3 5 2
4 3 1
4 1 2
4 2 1
5 1 2
5 2 1
5 3 1
5 4 4
6 7 5
7 1 0
7 2 4
7 3 0
7 . 3
8 .a 2
8 . 1
9 2 3
9 1 1
9 2 1
9 1 3
end
gen double w10 = w / 10
local pq (p1) p1=x (p10) p10=x (p25) p25=x (p33) p33=x (p50) p50=x (p66) p66=x (p75) p75=x (p90) p90=x (p99) p99=x (iqr) iq=x (median) md=x (count) n=x (sum) s=x (percent) pc=x
_cx_pair "percentile ties, no weights" 0 : `pq', by(g)
foreach wt in fw aw iw pw {
    _cx_pair "percentile ties [`wt']" 0 : `pq' [`wt'=w], by(g)
}
foreach wt in aw iw pw {
    _cx_pair "percentile ties, inexact weights [`wt']" 0 : `pq' [`wt'=w10], by(g)
}

* -----------------------------------------------------------------------
print_section "Extreme values"
* -----------------------------------------------------------------------
clear
input g double x
1 1e300
1 1e300
1 -1e300
3 -8.9e307
3 8.9e307
4 1e-300
4 1e-310
4 -1e-320
5 1e300
5 1e-300
5 -1e300
7 0
7 -0
8 1e200
8 -1e200
8 1e200
end
gen double wbig = cond(mod(_n, 2), 1e300, 1e-300)
gen double wnum = 1 + mod(_n, 3)
_cx_pair "huge and tiny values: sums and means" 0 : (sum) s=x (rawsum) rs=x (mean) m=x (count) n=x (percent) pc=x, by(g)
_cx_pair "huge and tiny values: sd and semean" 0 : (sd) sd=x (semean) se=x, by(g)
_cx_pair "huge and tiny values: order statistics" 0 : (min) mn=x (max) mx=x (median) md=x (p1) p1=x (p99) p99=x (first) fi=x (last) la=x (firstnm) fn=x (lastnm) ln=x, by(g)
_cx_pair "huge and tiny values: sebinomial and sepoisson" 0 : (sebinomial) sb=x (sepoisson) sp=x, by(g)
_cx_pair "huge and tiny weights [aw]" 0 : (mean) m=x (sum) s=x (sd) sd=x (median) md=x (percent) pc=x [aw=wbig], by(g)
_cx_pair "huge and tiny weights [iw]" 0 : (mean) m=x (sum) s=x (count) n=x (median) md=x (percent) pc=x [iw=wbig], by(g)
* the running total of the weights stays at P*W/100 over observations whose
* weights are too small to change it: collapse takes the last of them
preserve
clear
input g x double w
1 1 1
1 2 1e-20
1 3 1
2 1 1
2 2 1e-20
2 3 1e-20
2 4 1
end
foreach wt in aw iw pw {
    _cx_pair "median where weights do not move the total [`wt']" 0 : (p50) md=x [`wt'=w], by(g)
}
restore
_cx_pair "huge and tiny values, weights, no by()" 0 : (mean) m=x (sd) sd=x (p50) md=x (count) n=x [fw=wnum]
* sums beyond Stata's largest double: collapse stores the raw value (it
* lists as .z_), ccollapse stores a missing value
clear
input g double x
1 4.4e307
1 4.4e307
1 4.4e307
2 8.9e307
2 -8.9e307
end
preserve
ccollapse (sum) s=x (mean) m=x (iqr) q=x, by(g)
if s[1] == . & reldif(m[1], 4.4e307) < 1e-15 & q[1] == 0 & s[2] == 0 & q[2] == . test_pass "results beyond the largest double are missing"
else test_fail "results beyond the largest double are missing" "s `=s[1]' m `=m[1]' q `=q[2]'"
restore

* storage types at their limits
clear
set obs 15
gen byte g = ceil(_n / 3)
gen byte b = cond(mod(_n, 3) == 0, 100, cond(mod(_n, 3) == 1, -127, 0))
gen int i = cond(mod(_n, 3) == 0, 32740, cond(mod(_n, 3) == 1, -32767, 1))
gen long l = cond(mod(_n, 3) == 0, 2147483620, cond(mod(_n, 3) == 1, -2147483647, 7))
gen float f = cond(mod(_n, 3) == 0, 1.70141e38, cond(mod(_n, 3) == 1, -1.70141e38, 1e-38))
gen double d = cond(mod(_n, 3) == 0, 8.98846567431e307, cond(mod(_n, 3) == 1, -8.98846567431e307, 5e-324))
replace b = 100 if g == 2
replace i = 32740 if g == 2
replace l = 2147483620 if g == 2
replace f = 1.70141e38 if g == 2
foreach v in b i l f d {
    replace `v' = .z in 13
    replace `v' = .a in 14
}
format i %9.0gc
format l %12.0fc
_cx_pair "type limits: mean, sum, count" 0 : (mean) b i l f (sum) sb=b si=i sl=l sf=f (count) nb=b nl=l nd=d (rawsum) rl=l, by(g)
_cx_pair "type limits: order statistics" 0 : (max) b i l f d (min) nb=b ni=i nl=l nf=f nd=d (first) fb=b fl=l (lastnm) lf=f ld=d (median) mb=b mi=i ml=l mf=f, by(g)
_cx_pair "type limits: sd, semean, percent" 0 : (sd) b i l f (semean) sb=b sl=l (percent) pb=b pf=f, by(g)
set type double
_cx_pair "type limits, set type double" 0 : (mean) b i l f (p50) mb=b mf=f (iqr) qi=i, by(g)
set type float
* weights near the top of a long: counts and sums of weights pass 2^31
clear
set obs 9
gen byte g = ceil(_n / 3)
gen long w = 2000000000
replace w = 1 in 9
gen double x = _n
gen int y = _n * 1000
_cx_pair "fweights summing past 2^31" 0 : (count) n=x (sum) s=x (mean) m=x (percent) pc=x (median) md=x (rawsum) rs=x (p90) p90=y (sd) sd=y [fw=w], by(g)
_cx_pair "fweights summing past 2^31, no by()" 0 : (count) n=x (sum) s=x (mean) m=x (median) md=x (semean) se=y [fw=w]
_cx_pair "iweights of 1e15" 0 : (count) n=x (sum) s=x (mean) m=x (median) md=x [iw=w*500000], by(g)

* -----------------------------------------------------------------------
print_section "Missing values with weights"
* -----------------------------------------------------------------------
clear
input g x w
1 .a 1
1 .z 2
1 .  3
2 .b 0
2 3  1
2 .c 1
3 1  .
3 .d 2
4 .  .
4 2  4
4 .e 1
end
foreach wt in fw aw iw pw {
    _cx_pair "missing codes [`wt']" 0 : (max) mx=x (min) mn=x (first) fi=x (last) la=x (firstnm) fn=x (lastnm) ln=x (count) n=x (sum) s=x (mean) m=x (percent) pc=x (median) md=x [`wt'=w], by(g)
}

* -----------------------------------------------------------------------
print_section "Keys"
* -----------------------------------------------------------------------
* strings that differ only in spaces, case, accents, quotes, tabs and bytes
clear
set obs 400
gen str20 s = ""
replace s = "a" if mod(_n, 12) == 0
replace s = "a " if mod(_n, 12) == 1
replace s = " a" if mod(_n, 12) == 2
replace s = "A" if mod(_n, 12) == 3
replace s = uchar(233) if mod(_n, 12) == 4
replace s = "e" + uchar(769) if mod(_n, 12) == 5
replace s = uchar(26085) + uchar(26412) if mod(_n, 12) == 6
replace s = `"q"uote"' if mod(_n, 12) == 7
replace s = "tab" + char(9) if mod(_n, 12) == 8
replace s = char(255) + "z" if mod(_n, 12) == 9
replace s = "`" + "x'" if mod(_n, 12) == 10
gen str20 s2 = s
gen byte g2 = mod(_n, 3)
gen double x = _n
gen long id = _n
gen double u = runiform()
sort u
_cx_pair "string keys with spaces, case and bytes" 0 : (count) n=x (sum) x (min) mn=id (first) fi=id (last) la=id (firstnm) fs=s2, by(s)
_cx_pair "string keys with spaces, two keys" 0 : (count) n=x (sum) x (first) fi=id, by(g2 s)
* the longest str#: values that differ only in the last byte
clear
set obs 600
gen str2045 s = "a" * 2044 + string(mod(_n, 7))
replace s = "b" * 2045 if mod(_n, 11) == 0
replace s = "b" * 2044 if mod(_n, 13) == 0
gen int g = mod(_n, 5)
gen double x = _n
gen double u = runiform()
sort u
_cx_pair "str2045 key" 0 : (count) n=x (sum) x (first) fi=x (last) la=x, by(s)
sort g x
_cx_pair "str2045 source" 0 : (first) s (last) ls=s (firstnm) fs=s (lastnm) lns=s (count) n=x, by(g)
* doubles one ulp apart, the largest and smallest magnitudes, missing codes
clear
set obs 3000
gen double key = cond(mod(_n, 6) == 0, 1, cond(mod(_n, 6) == 1, 1 + epsdouble(), cond(mod(_n, 6) == 2, 1 - epsdouble() / 2, cond(mod(_n, 6) == 3, 8.98846567431e307, cond(mod(_n, 6) == 4, -8.98846567431e307, .)))))
replace key = .a if mod(_n, 7) == 0
replace key = 5e-324 if mod(_n, 11) == 0
replace key = 0 if mod(_n, 13) == 0
replace key = -1e-300 if mod(_n, 17) == 0
gen float fkey = 1 + floor(8 * runiform()) * epsfloat()
gen double x = _n
gen double u = runiform()
sort u
_cx_pair "double keys one ulp apart" 0 : (count) n=x (sum) x (first) fi=x (max) mx=x, by(key)
_cx_pair "float keys one ulp apart" 0 : (count) n=x (sum) x (first) fi=x, by(fkey)
* whole-number keys with negative values and missing codes, unsorted
clear
set obs 200000
gen int a = floor(41 * runiform()) - 20
replace a = . if runiform() < .01
replace a = .q if runiform() < .01
gen byte c = floor(3 * runiform()) - 1
replace c = .a if runiform() < .005
gen double r = floor(16 * runiform()) / 4
gen long id = _n
_cx_pair "whole-number keys with negatives and missing" 0 : (count) n=r (sum) r (min) mn=id (max) mx=id (first) fi=id (last) la=id (p50) md=r, by(a)
_cx_pair "two whole-number keys with negatives and missing" 0 : (count) n=r (sum) r (min) mn=id (last) la=id (p50) md=r, by(a c)
* whole-number keys whose ranges multiply past 2^64
clear
set obs 30000
gen long k1 = cond(runiform() < .5, 0, 4194303)
gen long k2 = cond(runiform() < .5, 0, 4194303)
gen long k3 = cond(runiform() < .5, 0, 4194303)
gen long m1 = cond(runiform() < .5, -2147483647, 2147483620)
gen long m2 = cond(runiform() < .5, -2147483647, 2147483620)
gen long q1 = cond(runiform() < .5, 0, 65535)
gen long q2 = cond(runiform() < .5, 0, 65535)
gen long q3 = cond(runiform() < .5, 0, 65535)
gen long q4 = cond(runiform() < .5, 0, 65535)
gen double r = floor(16 * runiform()) / 4
gen long id = _n
_cx_pair "three keys, ranges of 2^22" 0 : (count) n=r (sum) r (first) fi=id (last) la=id, by(k1 k2 k3)
_cx_pair "two keys spanning a long" 0 : (count) n=r (sum) r (first) fi=id (last) la=id, by(m1 m2)
_cx_pair "four keys, ranges of 2^16" 0 : (count) n=r (sum) r (first) fi=id (last) la=id, by(q1 q2 q3 q4)
replace k1 = floor(4194304 * runiform())
replace k2 = floor(4194304 * runiform())
replace k3 = floor(4194304 * runiform())
_cx_pair "three keys, ranges of 2^22, many groups" 0 : (count) n=r (sum) r (first) fi=id, by(k1 k2 k3)
_cx_pair "four mixed keys" 0 : (count) n=r (sum) r (first) fi=id, by(q1 m2 k1 q4)
* sorted runs that start again: the order check must look across chunks
clear
set obs 400000
gen long id = _n
gen double r = floor(8 * runiform()) / 8
foreach runs in 2 4 8 {
    capture drop key
    gen long key = mod(_n - 1, 400000 / `runs')
    _cx_pair "`runs' sorted runs, threads(`runs')" 0 threads(`runs') : (count) n=r (sum) r (first) fi=id (last) la=id, by(key)
}
replace key = _n
replace key = 0 in 400000
_cx_pair "sorted but for the last row" 0 : (count) n=r (sum) r (first) fi=id, by(key)
gsort -id
replace key = id
_cx_pair "sorted descending" 0 : (count) n=r (sum) r (first) fi=id (last) la=id, by(key)

* -----------------------------------------------------------------------
print_section "Large data and threads"
* -----------------------------------------------------------------------
clear
set obs 1000000
gen long a = floor(500000 * runiform())
gen double x = rnormal()
replace x = . if runiform() < .01
gen double r = floor(64 * runiform()) / 8
gen byte wi = 1 + floor(3 * runiform())
gen double wr = runiform() + .1
gen str5 code = string(floor(20000 * runiform()), "%05.0f")
gen long id = _n
local free (count) r (sum) sr=r (min) mn=x (max) mx=x (p50) md=x (p10) p10=r (iqr) iq=x (first) fi=id (last) la=id (percent) pc=x
_cx_pair "1e6 rows, 430,000 groups" 0 : `free', by(a)
_cx_pair "1e6 rows, string key, fweights" 0 : `free' [fw=wi], by(code)
_cx_pair "1e6 groups of one observation" 0 : (mean) x (sd) sd=x (semean) se=x (count) n=x (first) fr=r (p50) md=r, by(id)
_cx_same "threads(1) and threads(16): sort, aweights" "threads(1)" "threads(16)" : (mean) x r (sd) sx=x (p25) p25=x (p75) p75=r (median) md=x (count) n=x (first) fi=id [aw=wr], by(a)
_cx_same "threads(1) and threads(16): string key" "threads(1)" "threads(16)" : (mean) x (p99) p99=x (lastnm) ln=x (firstnm) fn=code (count) n=x [fw=wi], by(code)
_cx_same "threads(1) and threads(3): no by()" "threads(1)" "threads(3)" : (mean) x (sd) sx=x (p1) p1=x (p99) p99=x (sum) sr=r [aw=wr]

* -----------------------------------------------------------------------
print_section "Metadata, names and variables of other kinds"
* -----------------------------------------------------------------------
clear
set obs 300
gen int g = ceil(10 * runiform())
gen double abcdefghijklmnopqrstuvwxyz012345 = rnormal()
gen int kk = rpoisson(2)
gen double y = runiform()
label define lab1 1 "one"
label define lab2 2 "two"
label values kk lab1
label define unused 3 "three"
sort g
_cx_pair "32-character names" 0 : (mean) abcdefghijklmnopqrstuvwxyz012345 (sebinomial) zyxwvutsrqponmlkjihgfedcba543210 = abcdefghijklmnopqrstuvwxyz012345 (sum) kk, by(g)
_cx_pair "32-character names, clabel" 0 : (sepoisson) abcdefghijklmnopqrstuvwxyz012345 (p50) zyxwvutsrqponmlkjihgfedcba543210 = kk, by(g) clabel
* value label definitions stay in memory
preserve
collapse (mean) y, by(g)
quietly label dir
local labels_ref `r(names)'
restore, preserve
ccollapse (mean) y, by(g)
quietly label dir
local labels_new `r(names)'
restore
if "`labels_ref'" == "`labels_new'" test_pass "value label definitions kept"
else test_fail "value label definitions kept" "`labels_new' instead of `labels_ref'"
* variables named like Stata's temporary variables, with if and a weight expression
gen double __000000 = y * 2
gen double __000001 = kk + 1
gen double __000002 = 1
_cx_pair "sources named like temporary variables" 0 : (mean) __000000 (sum) __000001 (count) n=__000000 [aw=__000002 * 2] if __000000 > .2, by(g)
drop __00000*
* alias variables from a linked frame
capture frame drop __cx_link
frame create __cx_link
frame __cx_link {
    quietly set obs 10
    gen int g = _n
    gen double asrc = _n * 1.5
    gen int akey = mod(_n, 3)
}
quietly frlink m:1 g, frame(__cx_link)
quietly fralias add asrc akey, from(__cx_link)
_cx_pair "alias variable as a source" 0 : (mean) asrc (sum) s=asrc (max) y (p50) md=asrc, by(g)
* in the order of the key, where collapse too adds each group's values in the
* order of the observations (on other orders its sums vary in the last digits)
gen long __cx_o = _n
sort akey __cx_o
_cx_pair "alias variable as a key" 0 : (mean) y (count) n=y (sum) s=asrc, by(akey)
sort g __cx_o
drop __cx_o __cx_link asrc akey
frame drop __cx_link
* strL: collapse does the work, options passed on
gen strL sl = "v" + string(kk)
_cx_pair "strL source with clabel and fast" 0 : (first) sl (mean) y (count) n=y, by(g) clabel fast
_cx_pair "strL key with weights and cw" 0 : (mean) y (sum) s=kk [aw=y] if g > 2, by(sl) cw
drop sl
* time-series operators with weights and if
sort g
by g: gen int t = _n
quietly xtset g t
_cx_pair "time-series operators, weights, if" 0 : (mean) L.y D.kk (sum) s=L2.y [aw=y] if t > 3, by(g)
_cx_pair "operator groups" 0 : (mean) L(1/2).y (max) LD.(y kk), by(g)

* -----------------------------------------------------------------------
print_section "Long clists and the order of tied weighted values"
* -----------------------------------------------------------------------
* a target list far longer than a plugin's argument string
clear
set obs 200
gen int g = ceil(_n / 20)
forvalues j = 1/700 {
    gen float v`j' = runiform()
}
local cl (mean) v1-v700 (p50)
forvalues j = 1/700 {
    local cl `cl' m`j'=v`j'
}
local cl `cl' (sebinomial)
forvalues j = 1/700 {
    local cl `cl' b`j'=v`j'
}
_cx_pair "2100 targets" 0 : `cl', by(g)
_cx_pair "2100 targets, no by()" 0 : `cl'
* equal values are taken in the order of their weights; a running total that
* stays at P*W/100 over weights too small to change it (positive weights
* only: with a negative weight collapse takes the first of them)
clear
set obs 20000
gen long g = ceil(_n / 5)
gen double x = floor(5 * runiform())
gen double u = runiform()
gen double wp = cond(u < .3, 1, cond(u < .5, 2, cond(u < .7, 3, cond(u < .85, 1e-20, .5))))
gen double wn = cond(u < .15, -1, cond(u < .25, -2, cond(u < .5, 1, cond(u < .65, 2, cond(u < .8, 3, cond(u < .9, 1e-20, .5))))))
* groups whose iweights add to zero or less are kept: collapse rescales
* their weights to w*n/W (which turns the running total around) and gives
* missing when W is 0. A total that is nothing but the tiny weights (1 - 1
* + 1e-20) is 0 or 1e-20 depending on the order collapse adds in, which
* its result then follows, so there the tiny weights become .5
bysort g: egen double tw = total(cond(wn == 1e-20, 0, wn))
replace wn = .5 if wn == 1e-20 & tw == 0
drop tw u
gen long id = _n
gen double r = runiform()
sort r
local pq (p10) p10=x (p20) p20=x (p25) p25=x (p40) p40=x (p50) p50=x (p60) p60=x (p75) p75=x (p80) p80=x (p90) p90=x (iqr) iq=x
foreach wt in aw iw pw {
    _cx_pair "tied values, tiny weights [`wt']" 0 : `pq' [`wt'=wp], by(g)
}
_cx_pair "tied values, negative and tiny iweights" 0 : `pq' [iw=wn], by(g)
* iweights whose group totals are negative, zero or small, in exact
* arithmetic (integers and halves, so no rounding decides a boundary)
clear
set seed 20261127
set obs 50000
gen long g = ceil(_n / 10)
drop if runiform() < .4
gen double x = floor(8 * runiform())
gen double wi = (floor(7 * runiform()) - 3) / 2
replace wi = 1 if wi == 0
bysort g: egen double tw = total(wi)
quietly count if tw < 0
local nneg = r(N)
quietly count if tw == 0
local nzero = r(N)
drop tw
gen double r = runiform()
sort r
drop r
if `nneg' > 1000 & `nzero' > 100 test_pass "exact iweight data has negative (`nneg') and zero (`nzero') group totals"
else test_fail "exact iweight data has negative and zero group totals" "negative `nneg' zero `nzero'"
local pq2 (p1) p1=x (p5) p5=x (p10) p10=x (p25) p25=x (p33) p33=x (p50) p50=x (p67) p67=x (p75) p75=x (p90) p90=x (p99) p99=x (median) md=x (iqr) iq=x
_cx_pair "percentiles with negative, zero and small iweight totals" 0 : `pq2' [iw=wi], by(g)

noisily print_summary "ccollapse_extra"
global CTOOLS_COMPONENT_COMPLETE "ccollapse_extra"
global CTOOLS_CASE_COMPLETE "ccollapse_extra"
exit
}
