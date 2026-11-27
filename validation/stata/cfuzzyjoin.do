* cfuzzyjoin: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cfuzzyjoin"

if "`__case'" == "cfuzzyjoin" {
* Reference tests compare cfuzzyjoin link with reclink 1.7, cfuzzyjoin match
* with matchit 1.5.2 (and freqindex 1.3.1), and cfuzzyjoin group with
* strgroup 1.0.4, all from SSC.
do "validation/support/validate_setup.do"
foreach ref in reclink matchit freqindex strgroup {
    capture which `ref'
    if _rc {
        di as error "cfuzzyjoin suite: `ref' is not installed. Install the references from SSC:"
        di as error "    ssc install reclink"
        di as error "    ssc install freqindex"
        di as error "    ssc install matchit"
        di as error "    ssc install strgroup"
        exit 601
    }
}
which reclink
which matchit
which strgroup

* -----------------------------------------------------------------------------
* Data: names built from syllables, with typos, case changes, blanks
* -----------------------------------------------------------------------------
capture program drop _cfj_names
program define _cfj_names
    * _cfj_names n seed [long]: variable name (str), id 1..n
    args n seed long
    clear
    set seed `seed'
    quietly set obs `n'
    local syl "an ber cal dor el fin gar hol is jo kar lim mon nor ol pet quin ros sal tor ul vin wex yor zan"
    quietly gen str200 name = ""
    forvalues w = 1/3 {
        quietly gen int k`w' = 1 + floor(runiform() * 25)
        quietly gen int l`w' = 1 + floor(runiform() * 3)
        forvalues s = 1/3 {
            quietly replace name = name + word("`syl'", mod(k`w' * `s' + `w', 25) + 1) if l`w' >= `s'
        }
        if `w' < 3 quietly replace name = name + " " if runiform() < .8
        drop k`w' l`w'
    }
    if "`long'" != "" {
        * strings beyond 64 bytes take the banded distance
        quietly replace name = name + " " + name + " " + name + " " + name if runiform() < .5
    }
    quietly gen double u = runiform()
    quietly gen int p = 1 + floor(runiform() * max(1, length(name) - 1))
    quietly replace name = substr(name, 1, p - 1) + substr(name, p + 1, .) if u < .15
    quietly replace name = substr(name, 1, p) + substr(name, p, .) if u >= .15 & u < .3
    quietly replace name = substr(name, 1, p - 1) + substr(name, p + 1, 1) + substr(name, p, 1) ///
        + substr(name, p + 2, .) if u >= .3 & u < .4
    quietly replace name = upper(name) if runiform() < .05
    quietly replace name = "" if runiform() < .02
    quietly gen long id = _n
    drop u p
    quietly compress name
end

* _cfj_rows v1 [v2 ...] : "a" "b" ... makes a dataset of string variables
* v1 v2 ... from the quoted values, row by row (input cannot be used inside a
* program, and its end line would end one)
capture program drop _cfj_rows
program define _cfj_rows
    gettoken vars 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    local k : word count `vars'
    clear
    foreach v of local vars {
        quietly gen str244 `v' = ""
    }
    local i 0
    while `"`0'"' != "" {
        local ++i
        quietly set obs `i'
        foreach v of local vars {
            gettoken s 0 : 0
            quietly replace `v' = `"`s'"' in `i'
        }
    }
    quietly compress
end

* Odd strings: quotes, tabs, punctuation, digits, bytes above 127, short ones
capture program drop _cfj_odd
program define _cfj_odd
    _cfj_rows s1 s2 : "john smith" "smith, john" `"say "hi there" now"' `"say hi there now"' ///
        `"a"b c"d"' `"a"b"' "x" "x" "" "" "" "ab" "a" "ab" "ab" "ab" "abc" "abd" ///
        "MacDonald" "McDonald" "Knight" "Night" "Schmidt" "Smith" "Philips" "Phillips" ///
        "O'Hara" "OHara" "12 main st." "12 main street" "élan vital" "elan vital" ///
        "tab_inside" "tab inside" "  padded  " "padded" "eve" "eave" "ashcraft" "ashcroft" ///
        "wbwb bhb" "wb bh" "kn kt nd rd" "ee ie dt" "sch ph pf mac" "s f p m" ///
        "Johnatan Smithsonian" "Jonathan Smithson" ///
        `"`"x"' y"' "x y"
    quietly replace s1 = subinstr(s1, "_", char(9), .)
end

* =============================================================================
* group: strgroup
* =============================================================================
capture program drop _cfj_group_parity
program define _cfj_group_parity
    * strgroup and cfuzzyjoin group on the data in memory: same return code
    * and the same group numbers
    syntax anything(name=var) [if] [in], testname(string) [opts(string asis) by(varlist)]
    tempfile data
    quietly save `data'
    local byp
    if "`by'" != "" local byp by `by':
    capture `byp' strgroup `var' `if' `in', gen(__g1) `opts' force
    local rc1 = _rc
    quietly use `data', clear
    capture `byp' cfuzzyjoin group `var' `if' `in', gen(__g2) `opts'
    local rc2 = _rc
    if `rc1' != `rc2' {
        test_fail "`testname'" "strgroup r(`rc1'), cfuzzyjoin r(`rc2')"
        exit
    }
    if `rc1' {
        test_pass "`testname'"
        exit
    }
    quietly merge 1:1 _n using `data', nogen
    capture confirm variable __g1
    if _rc {
        quietly use `data', clear
        quietly `byp' strgroup `var' `if' `in', gen(__g1) `opts' force
        tempfile ref
        quietly save `ref'
        quietly use `data', clear
        quietly `byp' cfuzzyjoin group `var' `if' `in', gen(__g2) `opts'
        quietly merge 1:1 _n using `ref', nogen
    }
    capture assert __g1 == __g2
    if _rc {
        quietly count if __g1 != __g2
        test_fail "`testname'" "`r(N)' group numbers differ from strgroup"
    }
    else test_pass "`testname'"
end

_cfj_names 2500 9101
tempfile names
quietly save `names'
foreach o in "threshold(0.2)" "threshold(0.35) normalize(longer)" "threshold(0.15) normalize(shorter)" ///
    "threshold(2) normalize(none)" "threshold(1) normalize(none)" "threshold(0.3) first" ///
    "threshold(0.5) noclean" "threshold(1)" "threshold(0.05)" "threshold(3)" "threshold(.4) NORMalize(LONGER)" {
    quietly use `names', clear
    _cfj_group_parity name, testname("group equals strgroup: `o'") opts(`o')
}
quietly use `names', clear
_cfj_group_parity name if id > 200 in 1/2000, testname("group equals strgroup: if in") opts(threshold(0.25))
quietly use `names', clear
quietly gen byte grp = mod(id, 4)
sort grp id
_cfj_group_parity name, by(grp) testname("group equals strgroup: by") opts(threshold(0.25))
quietly use `names', clear
quietly gen byte grp = mod(id, 3)
quietly gen str1 grs = char(97 + mod(id, 2))
sort grp grs
_cfj_group_parity name, by(grp grs) testname("group equals strgroup: by two variables, first") ///
    opts(threshold(0.3) first)
* long strings (banded distance)
_cfj_names 1500 9102 long
foreach o in "threshold(0.1)" "threshold(0.25) normalize(longer)" "threshold(12) normalize(none)" {
    _cfj_group_parity name, testname("group equals strgroup, long strings: `o'") opts(`o')
    _cfj_names 1500 9102 long
}
* duplicates and thresholds at the boundary (1/7 = edit 1 of 7 bytes)
_cfj_rows s : "widgets" "widgetts" "widgets" "wodgets" "gadget" "gadgets" "gadgets" "x" "y" "xy" "abcdefg" "abcdefh"
foreach t in 0.142857142857 0.1428571428571429 0.14285714285714285 0.143 0.5 1 {
    preserve
    _cfj_group_parity s, testname("group boundary threshold `t'") opts(threshold(`t'))
    restore
}
* odd strings
_cfj_odd
foreach t in 0.3 1 {
    preserve
    _cfj_group_parity s1, testname("group odd strings `t'") opts(threshold(`t') noclean)
    restore
}
* errors: same return codes as strgroup
quietly use `names', clear
foreach o in "threshold(0)" "threshold(-1)" "threshold(1.5) normalize(none)" "threshold(0.2) normalize(middle)" ///
    "" {
    _cfj_group_parity name, testname("group error `o'") opts(`o')
}
quietly gen long n2 = id
capture cfuzzyjoin group id, gen(g) threshold(.2)
local rc2 = _rc
capture strgroup id, gen(g) threshold(.2)
if _rc == `rc2' & `rc2' test_pass "group of a numeric variable refused"
else test_fail "group of a numeric variable refused" "r(`rc2') vs strgroup r(`=_rc')"
capture cfuzzyjoin group name, gen(n2) threshold(.2)
if _rc == 110 test_pass "group into an existing variable refused"
else test_fail "group into an existing variable refused" "r(`=_rc')"
* blanks: trimmed unless noclean (strgroup 1.0.4 never trims)
_cfj_rows s : "abc" " abc" "abc  " "xyz" "   "
quietly cfuzzyjoin group s, gen(g1) threshold(0.01)
quietly cfuzzyjoin group s, gen(g2) threshold(0.01) noclean
quietly strgroup s, gen(g3) threshold(0.01)
capture assert g1 == 1 in 1/3
local a = _rc
capture assert g1[4] == 2 & missing(g1[5])
local b = _rc
capture assert g2 == g3
if !`a' & !`b' & !_rc test_pass "group trims blanks unless noclean"
else test_fail "group trims blanks unless noclean" "clean r(`a') r(`b'), noclean vs strgroup r(`=_rc')"
* force: accepted, nothing changes
quietly use `names', clear
quietly cfuzzyjoin group name, gen(g1) threshold(0.2) force
quietly cfuzzyjoin group name, gen(g2) threshold(0.2)
capture assert g1 == g2
if _rc test_fail "group force is accepted" "groups differ"
else test_pass "group force is accepted"
* stored results
quietly use `names', clear
quietly cfuzzyjoin group name, gen(g) threshold(0.2)
quietly count if !missing(g)
local n = r(N)
quietly summarize g
if r(N) == `n' & `=r(max)' > 0 test_pass "group stored results"
else test_fail "group stored results" "N"
quietly cfuzzyjoin group name, gen(g4) threshold(0.2)
local groups = r(groups)
quietly summarize g4
if `groups' == r(max) test_pass "group r(groups)"
else test_fail "group r(groups)" "`groups' vs `=r(max)'"

* =============================================================================
* match: matchit
* =============================================================================
capture program drop _cfj_cols_parity
program define _cfj_cols_parity
    * matchit and cfuzzyjoin match on two columns: the same scores
    syntax varlist(min=2 max=2), testname(string) [opts(string asis) tol(real 0)]
    tempfile data
    quietly save `data'
    capture matchit `varlist', `opts' gen(__s1)
    local rc1 = _rc
    if `rc1' {
        quietly use `data', clear
        capture cfuzzyjoin match `varlist', `opts' gen(__s2)
        if _rc == `rc1' | (_rc & `rc1') test_pass "`testname'"
        else test_fail "`testname'" "matchit r(`rc1'), cfuzzyjoin r(`=_rc')"
        quietly use `data', clear
        exit
    }
    keep __s1
    tempfile ref
    quietly save `ref'
    quietly use `data', clear
    capture cfuzzyjoin match `varlist', `opts' gen(__s2)
    if _rc {
        test_fail "`testname'" "cfuzzyjoin r(`=_rc')"
        exit
    }
    quietly merge 1:1 _n using `ref', nogen
    if `tol' == 0 capture assert __s1 == __s2
    else capture assert reldif(__s1, __s2) <= `tol'
    if _rc {
        quietly count if reldif(__s1, __s2) > `tol'
        test_fail "`testname'" "`r(N)' scores differ from matchit"
    }
    else test_pass "`testname'"
    quietly use `data', clear
end

_cfj_names 1500 9201
rename name name1
tempfile cols
quietly save `cols'
_cfj_names 1500 9202
rename name name2
quietly merge 1:1 id using `cols', nogen
quietly save `cols', replace
local methods `""bigram" "ngram, 3" "ngram, 1" "ngram_circ, 3" "ngram_circ, 1" "token" "cotoken" "scotoken" "firstgram, 2" "tkngram, 2" "soundex" "soundex_nara" "soundex_fk" "soundex_ext" "nysiis_fk" "token_soundex" `"tokenwrap, "soundex""' `"tokenwrap, "soundex_nara""' `"tokenwrap, "soundex_fk""' `"tokenwrap, "soundex_ext""' `"tokenwrap, "nysiis_fk""' "bigram, 1" "ngram, 2, 1" "token, 1" "tkngram, 3, 1" "token_soundex, 1""'
foreach sim of local methods {
    foreach sc in jaccard simple minsimple {
        quietly use `cols', clear
        _cfj_cols_parity name1 name2, testname(`"match columns equal matchit: `sim' `sc'"') opts(sim(`sim') score(`sc'))
    }
}
foreach sim in "bigram" "token" "soundex_ext" "ngram, 3" {
    foreach w in simple root log {
        quietly use `cols', clear
        _cfj_cols_parity name1 name2, testname("match columns weights equal matchit: `sim' `w'") ///
            opts(sim(`sim') weights(`w')) tol(1e-13)
    }
}
* odd strings, every method
_cfj_odd
tempfile odd
quietly save `odd'
foreach sim of local methods {
    quietly use `odd', clear
    _cfj_cols_parity s1 s2, testname(`"match columns odd strings: `sim'"') opts(sim(`sim'))
}
* weights from a file (grams missing from it count as 1)
_cfj_rows grams freq : "an" "3" "er" "10" "on" "1.5" "ol" "40"
quietly destring freq, replace
tempfile wfile
quietly save `wfile'
quietly use `cols', clear
_cfj_cols_parity name1 name2, testname("match columns wgtfile") opts(weights(log) wgtfile(`wfile')) tol(1e-13)
quietly use `cols', clear
_cfj_cols_parity name1 name2, testname("match columns wgtfile simple") opts(weights(simple) wgtfile(`wfile')) tol(1e-13)
* errors
quietly use `cols', clear
foreach o in "sim(nosuch)" "score(nosuch)" "weights(nosuch)" {
    _cfj_cols_parity name1 name2, testname("match columns error `o'") opts(`o')
}
capture cfuzzyjoin match name1 id
if _rc == 7 test_pass "match columns of a number refused"
else test_fail "match columns of a number refused" "r(`=_rc')"

capture program drop _cfj_index_parity
program define _cfj_index_parity
    * matchit and cfuzzyjoin match with using: the same pairs, scores,
    * variable names and storage types
    syntax anything(name=vars), using(string) testname(string) [opts(string asis) tol(real 0)]
    tempfile data ref
    quietly save `data'
    capture matchit `vars' using `"`using'"', `opts' override
    local rc1 = _rc
    if !`rc1' {
        quietly describe, varlist
        local v1 `r(varlist)'
        local t1
        foreach v of local v1 {
            local t1 `t1' `: type `v'':`: format `v''
        }
        local k1 = c(k)
        quietly ds
        local sortv : word 1 of `v1'
        local sortu : word 3 of `v1'
        sort `sortv' `sortu'
        quietly save `ref'
    }
    quietly use `data', clear
    capture cfuzzyjoin match `vars' using `"`using'"', `opts' override
    local rc2 = _rc
    if `rc1' != `rc2' {
        test_fail "`testname'" "matchit r(`rc1'), cfuzzyjoin r(`rc2')"
        quietly use `data', clear
        exit
    }
    if `rc1' {
        test_pass "`testname'"
        quietly use `data', clear
        exit
    }
    quietly describe, varlist
    local v2 `r(varlist)'
    local t2
    foreach v of local v2 {
        local t2 `t2' `: type `v'':`: format `v''
    }
    if "`v1'" != "`v2'" | "`t1'" != "`t2'" {
        test_fail "`testname'" "variables [`v1'] [`t1'] vs [`v2'] [`t2']"
        quietly use `data', clear
        exit
    }
    sort `sortv' `sortu'
    if `tol' == 0 {
        capture cf _all using `ref'
        local rc = _rc
    }
    else {
        local sc : word 5 of `v2'
        rename `sc' __mine
        quietly merge 1:1 `sortv' `sortu' using `ref'
        capture assert _merge == 3 & reldif(__mine, `sc') <= `tol'
        local rc = _rc
    }
    if `rc' test_fail "`testname'" "pairs or scores differ from matchit"
    else test_pass "`testname'"
    quietly use `data', clear
end

_cfj_names 900 9301
keep id name
tempfile mm uu
quietly save `mm'
_cfj_names 1100 9302
keep id name
rename (id name) (uid uname)
quietly save `uu'
foreach o in "" "threshold(0.7)" "threshold(0)" "threshold(1)" "sim(token) threshold(0.3)" "score(minsimple)" ///
    "score(simple) sim(ngram, 3)" "sim(soundex_fk)" `"sim(tokenwrap, "nysiis_fk") threshold(.4)"' ///
    "sim(scotoken) threshold(.2)" "sim(bigram, 1) threshold(.6)" "generate(sim)" "diagnose" "time flag(10) keepmata" {
    quietly use `mm', clear
    _cfj_index_parity id name, using(`uu') testname("match equals matchit: `o'") opts(idu(uid) txtu(uname) `o')
}
* matchit's weighted path calls a Mata function for every posting it reads,
* minutes per run at the size above: these cases use fewer records
quietly use `uu', clear
quietly keep in 1/300
tempfile uus
quietly save `uus'
* a weights file with every gram of both datasets (freqindex, as matchit's
* help suggests); grams missing from a file make matchit fail (below)
rename uname name
keep name
tempfile both wfull
quietly save `both'
quietly use `mm', clear
quietly keep in 1/250
keep name
quietly append using `both'
quietly freqindex name, sim(bigram)
quietly save `wfull'
foreach o in "weights(log)" "weights(simple) threshold(0.3)" "weights(root) sim(token) threshold(.3)" ///
    "stopwordsauto" "stopwordsauto swthreshold(0.05)" "stopwordsauto swthreshold(0.02) weights(root)" ///
    "weights(log) wgtfile(`wfull')" "weights(log) wgtfile(`wfile') stopwordsauto" {
    quietly use `mm', clear
    quietly keep in 1/250
    _cfj_index_parity id name, using(`uus') testname("match equals matchit: `o'") opts(idu(uid) txtu(uname) `o') tol(1e-13)
}
quietly use `mm', clear
quietly keep in 1/250
capture matchit id name using `uus', idu(uid) txtu(uname) weights(log) wgtfile(`wfile') override
local rc1 = _rc
quietly use `mm', clear
quietly keep in 1/250
capture cfuzzyjoin match id name using `uus', idu(uid) txtu(uname) weights(log) wgtfile(`wfile') override
if `rc1' & !_rc & _N > 0 test_pass "match wgtfile with grams the file lacks (matchit fails)"
else test_fail "match wgtfile with grams the file lacks (matchit fails)" "matchit r(`rc1'), cfuzzyjoin r(`=_rc')"
* clashing names: matchit renames the using variables and the score
quietly use `uu', clear
rename (uid uname) (id name)
tempfile uu2
quietly save `uu2'
quietly use `mm', clear
_cfj_index_parity id name, using(`uu2') testname("match renames clashing variables") opts(idu(id) txtu(name))
quietly use `uu', clear
rename uname similscore
tempfile uu3
quietly save `uu3'
* (matchit never returns here: its loop that renames the score compares the
* requested name, not the candidate, with txtusing())
quietly use `mm', clear
capture cfuzzyjoin match id name using `uu3', idu(uid) txtu(similscore) override
local rc = _rc
capture confirm numeric variable similscore1
local a = _rc
capture confirm string variable similscore
if !`rc' & !`a' & !_rc test_pass "match renames the score"
else test_fail "match renames the score" "r(`rc')"
* no pairs: the result keeps the storage types and has no observations
quietly use `mm', clear
_cfj_index_parity id name, using(`uu') testname("match without pairs") opts(idu(uid) txtu(uname) sim(soundex) threshold(1.1))
* odd strings
_cfj_odd
quietly gen long id = _n
keep id s1
tempfile om ou
quietly save `om'
_cfj_odd
quietly gen long uid = _n
keep uid s2
quietly save `ou'
foreach sim in "bigram" "token" "cotoken" "soundex" "nysiis_fk" "ngram_circ, 2" {
    quietly use `om', clear
    _cfj_index_parity id s1, using(`ou') testname("match odd strings: `sim'") opts(idu(uid) txtu(s2) sim(`sim') threshold(0))
}
* unsaved changes need override
quietly use `mm', clear
quietly replace name = name + "x" in 1
capture cfuzzyjoin match id name using `uu', idu(uid) txtu(uname)
local rc2 = _rc
capture matchit id name using `uu', idu(uid) txtu(uname)
if `rc2' == 4 & _rc == 4 test_pass "match refuses to replace changed data"
else test_fail "match refuses to replace changed data" "r(`rc2') vs matchit r(`=_rc')"
* the data stay when the using file is wrong
quietly use `mm', clear
capture cfuzzyjoin match id name using `uu', idu(nosuch) txtu(uname)
local rc2 = _rc
capture confirm variable name
if `rc2' == 111 & !_rc & _N == 900 test_pass "match error leaves the data"
else test_fail "match error leaves the data" "r(`rc2')"
quietly use `mm', clear
capture cfuzzyjoin match id name using `uu', idu(uid) txtu(uname) sim(ngram, 0)
if _rc == 198 & _N == 900 test_pass "match refuses an invalid gram size"
else test_fail "match refuses an invalid gram size" "r(`=_rc')"

* =============================================================================
* link: reclink
* =============================================================================
capture program drop _cfj_link_parity
program define _cfj_link_parity
    * reclink and cfuzzyjoin link: the same return code, variables (order,
    * types, formats, labels) and observations
    syntax anything(name=vars), using(string) testname(string) opts(string asis) ///
        idm(name) idu(name) [nodescribe]
    tempfile data ref
    quietly save `data'
    capture reclink `vars' using `"`using'"', idm(`idm') idu(`idu') `opts'
    local rc1 = _rc
    if !`rc1' {
        _cfj_describe
        local d1 `"`r(desc)'"'
        sort `idm' `idu'
        quietly save `ref'
    }
    quietly use `data', clear
    capture cfuzzyjoin link `vars' using `"`using'"', idm(`idm') idu(`idu') `opts'
    local rc2 = _rc
    if `rc1' != `rc2' {
        test_fail "`testname'" "reclink r(`rc1'), cfuzzyjoin r(`rc2')"
        quietly use `data', clear
        exit
    }
    if `rc1' {
        test_pass "`testname'"
        quietly use `data', clear
        exit
    }
    _cfj_describe
    local d2 `"`r(desc)'"'
    if "`describe'" == "" & `"`d1'"' != `"`d2'"' {
        test_fail "`testname'" `"variables differ: [`d1'] [`d2']"'
        quietly use `data', clear
        exit
    }
    sort `idm' `idu'
    capture cf _all using `ref'
    if _rc test_fail "`testname'" "observations differ from reclink"
    else test_pass "`testname'"
    quietly use `data', clear
end

capture program drop _cfj_describe
program define _cfj_describe, rclass
    * names, storage types, formats, value labels, variable labels, then the
    * characteristics (notes included) of each variable, in sorted order
    local d
    foreach v of varlist _all {
        local d `d' `v':`: type `v'':`: format `v'':`: value label `v'':`: var label `v''
        local chars : char `v'[]
        local chars : list sort chars
        foreach c of local chars {
            local d `d' `v'[`c']=`: char `v'[`c']'
        }
    }
    quietly label dir
    local labels `r(names)'
    local labels : list sort labels
    local d `d' labels:`labels'
    return local desc `"`d'"'
end

_cfj_names 600 9401
keep id name
quietly gen str12 city = word("boston newton cambridge somerville quincy", 1 + mod(id, 5))
quietly replace city = "" if mod(id, 23) == 0
quietly gen zip = 2100 + mod(id, 37)
quietly replace zip = . if mod(id, 31) == 0
quietly gen int year = 1990 + mod(id, 9)
label variable city "city of the master"
tempfile lm lu
quietly save `lm'
_cfj_names 800 9402
keep id name
rename id uid
* one record per name: where several master records share the best match among
* duplicated using records, reclink pairs them off in an arbitrary order
* (cfuzzyjoin keeps every pair; tested below)
quietly bysort name (uid): keep if _n == 1
sort uid
quietly gen str12 city = word("boston newton cambridge somerville quincy", 1 + mod(uid, 5))
quietly replace city = "brookline" if mod(uid, 17) == 0
quietly gen zip = 2100 + mod(uid, 37)
quietly replace zip = . if mod(uid, 29) == 0
quietly gen int year = 1990 + mod(uid, 7)
quietly gen str6 extra = "x" + string(mod(uid, 7))
label variable extra "a using variable"
label define yl 1990 "first"
label values year yl
notes extra: a note on extra
char extra[mychar] hello
quietly save `lu'
foreach spec in "name" "name city" "name city zip" "name city zip|wmatch(10 2 4) wnomatch(8 1 2)" ///
    "name city zip|required(city)" "name city zip year|orblock(city zip)" "name city zip|minscore(0.8) minbigram(0.7)" ///
    "name city|exactstr(city)" "name city zip year|orblock(none)" "name zip|uprefix(V_)" ///
    "name city zip year" "name city zip|required(city zip)" "name city|_merge(how)" "name city zip|minscore(0.3)" ///
    "name city|wmatch(5 1) wnomatch(0.5 2)" {
    gettoken vars opts : spec, parse("|")
    local opts = substr(`"`opts'"', 2, .)
    quietly use `lm', clear
    _cfj_link_parity `vars', using(`lu') idm(id) idu(uid) testname("link equals reclink: `spec'") opts(gen(score) `opts')
}
* or-blocking on name: reclink fails (syntax error) on a master record whose
* name is empty, so these records are left out here
quietly use `lm', clear
quietly drop if name == ""
_cfj_link_parity name year, using(`lu') idm(id) idu(uid) testname("link equals reclink: required() and orblock()") ///
    opts(gen(score) required(year) orblock(name))
quietly use `lm', clear
capture reclink name year using `lu', idm(id) idu(uid) gen(score) required(year) orblock(name)
local rc1 = _rc
quietly use `lm', clear
capture cfuzzyjoin link name year using `lu', idm(id) idu(uid) gen(score) required(year) orblock(name)
if `rc1' & !_rc test_pass "link or-blocking with empty master strings (reclink fails)"
else test_fail "link or-blocking with empty master strings (reclink fails)" "reclink r(`rc1'), cfuzzyjoin r(`=_rc')"
* uvarlist: other names in the using data
quietly use `lu', clear
rename (name city) (uname ucity)
tempfile lu2
quietly save `lu2'
quietly use `lm', clear
_cfj_link_parity name city, using(`lu2') idm(id) idu(uid) testname("link equals reclink: uvarlist") ///
    opts(gen(score) uvarlist(uname ucity))
* exclude(): earlier matches
quietly use `lm', clear
quietly reclink name city using `lu', idm(id) idu(uid) gen(score) minscore(.9)
quietly keep if _merge == 3
keep id uid
tempfile ex
quietly save `ex'
quietly use `lm', clear
_cfj_link_parity name city, using(`lu') idm(id) idu(uid) testname("link equals reclink: exclude") ///
    opts(gen(score) exclude(`ex'))
* debug(): the pairs file
tempfile dbg1 dbg2
quietly use `lm', clear
quietly reclink name city using `lu', idm(id) idu(uid) gen(score) debug(`dbg1')
quietly use `lm', clear
quietly cfuzzyjoin link name city using `lu', idm(id) idu(uid) gen(score) debug(`dbg2')
quietly use `dbg2', clear
sort id uid
quietly save `dbg2', replace
quietly use `dbg1', clear
sort id uid
capture cf _all using `dbg2'
if _rc test_fail "link debug file equals reclink's" "cf r(`=_rc')"
else test_pass "link debug file equals reclink's"
* ties and duplicated using records
_cfj_rows name : "abcd" "wxyz" "qqqq"
quietly gen long id = _n
tempfile tm tu
quietly save `tm'
_cfj_rows name : "abce" "abcf" "wxyz" "wxyz" "qqqr"
quietly gen long uid = _n + 4
quietly save `tu'
quietly use `tm', clear
_cfj_link_parity name, using(`tu') idm(id) idu(uid) testname("link ties and duplicates equal reclink") ///
    opts(gen(score) minscore(.1))
* two master records and two duplicated using records: every pair
_cfj_rows name : "wxyy" "wxyy"
quietly gen long id = _n
tempfile tm2
quietly save `tm2'
quietly cfuzzyjoin link name using `tu', idm(id) idu(uid) gen(score)
quietly count if _merge == 3
local n3 = r(N)
quietly count if (id == 1 | id == 2) & (uid == 7 | uid == 8)
if `n3' == 4 & r(N) == 4 test_pass "link gives every duplicate using record to every master record"
else test_fail "link gives every duplicate using record to every master record" "`n3' matched rows"
* blanks in master strings reach the comparator as reclink's syntax passes
* them: spaces trimmed and runs of them made one, tabs kept (reclink fails on
* double quotes)
_cfj_rows name : "ab  cd" " ab cd" "ab_cd" "abcd  "
quietly replace name = subinstr(name, "_", char(9), .)
quietly gen long id = _n
quietly save `tm2', replace
_cfj_rows name : "ab cd" "abxcd"
quietly gen long uid = _n + 4
tempfile tu2
quietly save `tu2'
quietly use `tm2', clear
_cfj_link_parity name, using(`tu2') idm(id) idu(uid) testname("link master blanks equal reclink") ///
    opts(gen(score) minscore(.1))
* string identifiers
quietly use `lm', clear
quietly tostring id, replace format(%05.0f)
quietly replace id = "m" + id
tempfile sm su
quietly save `sm'
quietly use `lu', clear
quietly tostring uid, replace
quietly replace uid = "u" + uid
quietly save `su'
quietly use `sm', clear
capture reclink name city using `su', idm(id) idu(uid) gen(score)
local rc1 = _rc
quietly use `sm', clear
capture cfuzzyjoin link name city using `su', idm(id) idu(uid) gen(score)
local rc2 = _rc
keep id uid score
quietly replace id = substr(id, 2, .)
quietly destring id, replace
quietly replace uid = substr(uid, 2, .)
quietly destring uid, replace
sort id uid
tempfile sres
quietly save `sres'
quietly use `lm', clear
quietly cfuzzyjoin link name city using `lu', idm(id) idu(uid) gen(score)
keep id uid score
sort id uid
capture cf _all using `sres'
if `rc1' & !`rc2' & !_rc test_pass "link with string identifiers (reclink fails)"
else test_fail "link with string identifiers (reclink fails)" "reclink r(`rc1'), cfuzzyjoin r(`rc2'), cf r(`=_rc')"
* errors
quietly use `lm', clear
foreach spec in "name|wmatch(0.5)" "name city|wmatch(1)" "name|uvarlist(name city)" "name|required(city)" ///
    "name|orblock(city)" "name|gen(name)" {
    gettoken vars opts : spec, parse("|")
    local opts = substr(`"`opts'"', 2, .)
    if !strpos(`"`opts'"', "gen(") local opts `opts' gen(score)
    _cfj_link_parity `vars', using(`lu') idm(id) idu(uid) testname("link error equals reclink: `spec'") opts(`opts')
}
quietly use `lm', clear
quietly replace id = 1 in 2
_cfj_link_parity name, using(`lu') idm(id) idu(uid) testname("link error: idmaster not unique") opts(gen(score))
quietly use `lm', clear
quietly gen byte _merge = 1
_cfj_link_parity name, using(`lu') idm(id) idu(uid) testname("link error: _merge exists") opts(gen(score))
quietly use `lu', clear
quietly replace uid = 3 in 5
tempfile ludup
quietly save `ludup'
quietly use `lm', clear
_cfj_link_parity name, using(`ludup') idm(id) idu(uid) testname("link error: idusing not unique") opts(gen(score))
quietly use `lm', clear
capture cfuzzyjoin link name using `lu', idm(id) idu(uid) gen(score) threads(1)
if !_rc test_pass "link with one thread"
else test_fail "link with one thread" "r(`=_rc')"

* =============================================================================
* Results do not depend on the thread count
* =============================================================================
quietly use `names', clear
quietly cfuzzyjoin group name, gen(g1) threshold(0.3) threads(1)
quietly cfuzzyjoin group name, gen(g2) threshold(0.3) threads(4)
capture assert g1 == g2
if _rc test_fail "group threads(1) equals threads(4)" ""
else test_pass "group threads(1) equals threads(4)"
foreach o in "" "weights(log)" {
    quietly use `mm', clear
    quietly cfuzzyjoin match id name using `uu', idu(uid) txtu(uname) `o' threads(1)
    tempfile t1
    quietly save `t1'
    quietly use `mm', clear
    quietly cfuzzyjoin match id name using `uu', idu(uid) txtu(uname) `o' threads(4) verbose
    capture cf _all using `t1'
    if _rc test_fail "match threads(1) equals threads(4) `o'" "cf r(`=_rc')"
    else test_pass "match threads(1) equals threads(4) `o'"
}
quietly use `lm', clear
quietly cfuzzyjoin link name city zip using `lu', idm(id) idu(uid) gen(score) threads(1)
tempfile t1
quietly save `t1'
quietly use `lm', clear
quietly cfuzzyjoin link name city zip using `lu', idm(id) idu(uid) gen(score) threads(4) verbose
capture cf _all using `t1'
if _rc test_fail "link threads(1) equals threads(4)" "cf r(`=_rc')"
else test_pass "link threads(1) equals threads(4)"
quietly use `cols', clear
quietly cfuzzyjoin match name1 name2, gen(a) weights(root) threads(1)
quietly cfuzzyjoin match name1 name2, gen(b) weights(root) threads(4) verbose
capture assert a == b
if _rc test_fail "match columns threads(1) equals threads(4)" ""
else test_pass "match columns threads(1) equals threads(4)"
quietly use `names', clear
quietly cfuzzyjoin group name, gen(g1) threshold(0.3) verbose
if !_rc test_pass "group verbose"
capture cfuzzyjoin group name, gen(g9) threshold(0.3) threads(300)
if _rc == 198 test_pass "threads() out of range refused"
else test_fail "threads() out of range refused" "r(`=_rc')"
capture cfuzzyjoin nosuch name
if _rc == 198 test_pass "unknown subcommand refused"
else test_fail "unknown subcommand refused" "r(`=_rc')"
quietly use `lm', clear
quietly gen byte grp = mod(id, 2)
sort grp
capture by grp: cfuzzyjoin link name using `lu', idm(id) idu(uid) gen(score)
if _rc == 190 test_pass "link with by refused"
else test_fail "link with by refused" "r(`=_rc')"

print_summary "cfuzzyjoin"
global CTOOLS_COMPONENT_COMPLETE "cfuzzyjoin"
    global CTOOLS_CASE_COMPLETE "cfuzzyjoin"
    exit
}

di as error "Unknown cfuzzyjoin case: `__case'"
exit 198
