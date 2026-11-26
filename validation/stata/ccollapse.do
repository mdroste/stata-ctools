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
