* cfuzzyjoin against reclink, matchit and strgroup. Run from the repository root:
*     do validation/performance/cfuzzyjoin.do cfuzzyjoin_performance [sizes] [large]
* sizes: observation counts of the comparisons (default "1000 2000 5000");
* large: observation counts timed for cfuzzyjoin alone (default "50000").
* Each line FUZZY_PERF,<task>,<N>,<method>,<seconds>,<check> is one timed run;
* the reference commands run once, cfuzzyjoin best of ten (best of three at
* the large sizes). check 0 means
* the two results agree (reclink: pairs and scores; matchit: pairs and scores;
* strgroup: group numbers).
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}
if "`__workload'" == "cfuzzyjoin_performance" {
args sizes large
if `"`sizes'"' == "" local sizes "1000 2000 5000"
if `"`large'"' == "" local large "50000"
clear all
set more off
set linesize 255
adopath ++ "build"
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
about
ctools, version
foreach cmd in cfuzzyjoin reclink matchit freqindex strgroup {
    which `cmd'
}

* Company-like names: two or three words from a vocabulary of 25 syllables,
* with a typo (deleted, doubled or swapped byte) in 40% of them, upper case
* in 5%, and a city and a ZIP code for the record linkage.
capture program drop _fz_names
program define _fz_names
    args n seed
    clear
    set seed `seed'
    quietly set obs `n'
    local syl "an ber cal dor el fin gar hol is jo kar lim mon nor ol pet quin ros sal tor ul vin wex yor zan"
    quietly gen str80 name = ""
    forvalues w = 1/3 {
        quietly gen int k`w' = 1 + floor(runiform() * 25)
        quietly gen int k`w'b = 1 + floor(runiform() * 25)
        quietly gen int l`w' = 1 + floor(runiform() * 3)
        forvalues s = 1/3 {
            quietly replace name = name + word("`syl'", mod(k`w' * `s' + k`w'b, 25) + 1) if l`w' >= `s'
        }
        if `w' < 3 quietly replace name = name + " " if runiform() < .85
        drop k`w' k`w'b l`w'
    }
    quietly replace name = strtrim(stritrim(name))
    quietly gen double u = runiform()
    quietly gen int p = 1 + floor(runiform() * max(1, length(name) - 1))
    quietly replace name = substr(name, 1, p - 1) + substr(name, p + 1, .) if u < .15
    quietly replace name = substr(name, 1, p) + substr(name, p, .) if u >= .15 & u < .3
    quietly replace name = substr(name, 1, p - 1) + substr(name, p + 1, 1) + substr(name, p, 1) ///
        + substr(name, p + 2, .) if u >= .3 & u < .4
    quietly replace name = upper(name) if runiform() < .05
    quietly gen long id = _n
    quietly gen str12 city = word("boston newton cambridge somerville quincy medford salem lynn", 1 + floor(runiform() * 8))
    quietly gen long zip = 2100 + floor(runiform() * 60)
    drop u p
    quietly compress
end

capture program drop _fz_time
program define _fz_time, rclass
    * _fz_time reps : command  -> r(seconds) the best of reps runs
    gettoken reps 0 : 0
    gettoken colon 0 : 0
    local best .
    forvalues r = 1/`reps' {
        preserve
        timer clear 99
        timer on 99
        quietly `0'
        timer off 99
        quietly timer list 99
        local best = min(`best', r(t99))
        if `r' < `reps' restore
        else restore, not
    }
    return scalar seconds = `best'
end

foreach N of local sizes {
    * ------------------------------------------------------------ strgroup
    _fz_names `N' 4101
    tempfile g ref
    quietly save `g'
    timer clear 99
    timer on 99
    quietly strgroup name, gen(grp) threshold(0.2) force
    timer off 99
    quietly timer list 99
    local tref = r(t99)
    keep id grp
    quietly save `ref'
    quietly use `g', clear
    _fz_time 10 : cfuzzyjoin group name, gen(grp) threshold(0.2)
    local tc = r(seconds)
    quietly use `g', clear
    quietly cfuzzyjoin group name, gen(grp2) threshold(0.2)
    quietly merge 1:1 id using `ref', nogen
    capture assert grp2 == grp
    local chk = _rc
    di "FUZZY_PERF,group,`N',ref,`tref',`chk'"
    di "FUZZY_PERF,group,`N',c,`tc',`chk'"

    * ------------------------------------------------------------ matchit
    _fz_names `N' 4201
    keep id name
    tempfile mm uu
    quietly save `mm'
    _fz_names `N' 4202
    keep id name
    rename (id name) (uid uname)
    quietly save `uu'
    quietly use `mm', clear
    timer clear 99
    timer on 99
    quietly matchit id name using `uu', idu(uid) txtu(uname) threshold(0.6)
    timer off 99
    quietly timer list 99
    local tref = r(t99)
    sort id uid
    quietly save `ref', replace
    quietly use `mm', clear
    _fz_time 10 : cfuzzyjoin match id name using `uu', idu(uid) txtu(uname) threshold(0.6)
    local tc = r(seconds)
    quietly use `mm', clear
    quietly cfuzzyjoin match id name using `uu', idu(uid) txtu(uname) threshold(0.6)
    sort id uid
    capture cf _all using `ref'
    local chk = _rc
    di "FUZZY_PERF,match,`N',ref,`tref',`chk'"
    di "FUZZY_PERF,match,`N',c,`tc',`chk'"

    * ------------------------------------------------------------ reclink
    _fz_names `N' 4301
    tempfile lm lu
    quietly save `lm'
    _fz_names `N' 4302
    rename id uid
    * one record per name and city (see help cfuzzyjoin on duplicated records)
    quietly bysort name city zip (uid): keep if _n == 1
    quietly save `lu'
    quietly use `lm', clear
    timer clear 99
    timer on 99
    quietly reclink name city zip using `lu', idm(id) idu(uid) gen(score) wmatch(10 2 4) wnomatch(8 1 2)
    timer off 99
    quietly timer list 99
    local tref = r(t99)
    sort id uid
    quietly save `ref', replace
    quietly use `lm', clear
    _fz_time 10 : cfuzzyjoin link name city zip using `lu', idm(id) idu(uid) gen(score) wmatch(10 2 4) wnomatch(8 1 2)
    local tc = r(seconds)
    quietly use `lm', clear
    quietly cfuzzyjoin link name city zip using `lu', idm(id) idu(uid) gen(score) wmatch(10 2 4) wnomatch(8 1 2)
    sort id uid
    capture cf _all using `ref'
    local chk = _rc
    di "FUZZY_PERF,link,`N',ref,`tref',`chk'"
    di "FUZZY_PERF,link,`N',c,`tc',`chk'"
}

* cfuzzyjoin alone at larger sizes
foreach N of local large {
    _fz_names `N' 4401
    _fz_time 3 : cfuzzyjoin group name, gen(grp) threshold(0.2)
    di "FUZZY_PERF,group,`N',c,`r(seconds)',."
    _fz_names `N' 4402
    keep id name
    tempfile mm uu
    quietly save `mm'
    _fz_names `N' 4403
    keep id name
    rename (id name) (uid uname)
    quietly save `uu'
    quietly use `mm', clear
    _fz_time 3 : cfuzzyjoin match id name using `uu', idu(uid) txtu(uname) threshold(0.6)
    di "FUZZY_PERF,match,`N',c,`r(seconds)',."
    _fz_names `N' 4404
    tempfile lm lu
    quietly save `lm'
    _fz_names `N' 4405
    rename id uid
    quietly bysort name city zip (uid): keep if _n == 1
    quietly save `lu'
    quietly use `lm', clear
    _fz_time 3 : cfuzzyjoin link name city zip using `lu', idm(id) idu(uid) gen(score) required(city)
    di "FUZZY_PERF,link_required,`N',c,`r(seconds)',."
}
di "FUZZY_PERF_COMPLETE"
exit
}
