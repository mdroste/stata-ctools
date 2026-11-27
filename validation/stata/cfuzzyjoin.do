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

if "`__case'" == "cfuzzyjoin_extra" {
* Extra checks (2026-11-26): edge cases of the three subcommands against
* strgroup 1.0.4, matchit 1.5.2 and reclink 1.7: tiny and empty data,
* unicode, case and punctuation, very long strings and strL, thresholds at
* the cutoff, ties and many-to-many pairs, gram sizes, weights and scores,
* abbreviations, data left unchanged by errors, sort order and threads.

* Checks that could bring Stata down (data with no observations) run in a
* separate Stata process: the suite calls itself as
*     do <this file> cfuzzyjoin_extra child <check> <file>
* and the child writes "1" (or "0 <why>") to <file>, so a crash fails that
* check instead of the whole suite.
if "`1'" == "child" {
    local check `2'
    local resfile `"`3'"'
    local ok 0
    local why "unknown check"
    tempfile m u
    clear
    if "`check'" == "group" {
        quietly set obs 0
        quietly gen str5 s = ""
        * (strgroup stops with r(198) here, from its duplicates warning)
        capture cfuzzyjoin group s, gen(g2) threshold(0.2)
        local rc2 = _rc
        capture confirm numeric variable g2
        local ok = !`rc2' & !_rc
        local why "cfuzzyjoin r(`rc2')"
    }
    else if "`check'" == "columns" {
        quietly set obs 0
        quietly gen str5 s1 = ""
        quietly gen str5 s2 = ""
        capture cfuzzyjoin match s1 s2, gen(sc)
        local rc2 = _rc
        capture confirm numeric variable sc
        local ok = !`rc2' & !_rc & _N == 0
        local why "cfuzzyjoin r(`rc2'), N=`=_N' (matchit r(0))"
    }
    else if inlist("`check'", "match_master", "match_using", "link_master", "link_using") {
        quietly set obs 1
        quietly gen str20 name = "hello world"
        quietly gen long id = 1
        if "`check'" == "match_master" | "`check'" == "link_master" quietly drop in 1
        quietly save `m', emptyok
        clear
        quietly set obs 1
        quietly gen str20 uname = "hello word"
        quietly gen long uid = 1
        if "`check'" == "match_using" | "`check'" == "link_using" quietly drop in 1
        quietly save `u', emptyok
        quietly use `m', clear
        if substr("`check'", 1, 5) == "match" {
            capture matchit id name using `u', idu(uid) txtu(uname) override
            local rc1 = _rc
            local n1 = _N
            local k1 = c(k)
            quietly use `m', clear
            capture cfuzzyjoin match id name using `u', idu(uid) txtu(uname) override
        }
        else {
            capture reclink name using `u', idm(id) idu(uid) gen(score) uvarlist(uname)
            local rc1 = _rc
            local n1 = _N
            local k1 = c(k)
            quietly use `m', clear
            capture cfuzzyjoin link name using `u', idm(id) idu(uid) gen(score) uvarlist(uname)
        }
        local rc2 = _rc
        local ok = `rc1' == `rc2' & `n1' == _N & `k1' == c(k)
        local why "reference r(`rc1') N=`n1' k=`k1', cfuzzyjoin r(`rc2') N=`=_N' k=`=c(k)'"
    }
    tempname fh
    file open `fh' using `"`resfile'"', write text replace
    file write `fh' `"`ok' `why'"' _n
    file close `fh'
    exit
}

do "validation/support/validate_setup.do"
foreach ref in reclink matchit strgroup {
    capture which `ref'
    if _rc {
        di as error "cfuzzyjoin_extra: `ref' is not installed (ssc install `ref')"
        exit 601
    }
}

* _cfx_rows v1 [v2 ...] : "a" "b" ... : string variables from quoted values,
* row by row
capture program drop _cfx_rows
program define _cfx_rows
    gettoken vars 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
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

* _cfx_names n seed: name (str) and id 1..n from a small vocabulary with
* accented syllables (many near names), with typos made by characters
capture program drop _cfx_names
program define _cfx_names
    args n seed
    clear
    set seed `seed'
    quietly set obs `n'
    local syl "an bér cal dör el fin gär hol ïs jo"
    quietly gen str200 name = ""
    forvalues w = 1/2 {
        quietly gen int k = 1 + floor(runiform() * 10)
        quietly gen int l = 1 + floor(runiform() * 2)
        forvalues s = 1/2 {
            quietly replace name = name + word("`syl'", mod(k * `s' + `w', 10) + 1) if l >= `s'
        }
        if `w' < 2 quietly replace name = name + " " if runiform() < .8
        drop k l
    }
    quietly gen double u = runiform()
    quietly gen int p = 1 + floor(runiform() * max(1, ustrlen(name) - 1))
    quietly replace name = usubstr(name, 1, p - 1) + usubstr(name, p + 1, .) if u < .2
    quietly replace name = usubstr(name, 1, p) + usubstr(name, p, .) if u >= .2 & u < .35
    quietly replace name = ustrupper(name) if runiform() < .05
    quietly replace name = upper(name) if runiform() < .05
    quietly replace name = "" if runiform() < .03
    quietly gen long id = _n
    drop u p
    quietly compress name
end

capture program drop _cfx_gpar
program define _cfx_gpar
    * strgroup and cfuzzyjoin group: same return code and group numbers
    syntax anything(name=var) [if] [in], testname(string) [opts(string asis) by(varlist)]
    tempfile data ref
    quietly save `data'
    local byp
    if "`by'" != "" local byp by `by':
    capture `byp' strgroup `var' `if' `in', gen(__g1) `opts' force
    local rc1 = _rc
    if !`rc1' {
        keep __g1
        quietly save `ref'
    }
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
    quietly merge 1:1 _n using `ref', nogen
    capture assert __g1 == __g2
    if _rc {
        quietly count if __g1 != __g2
        test_fail "`testname'" "`r(N)' group numbers differ from strgroup"
    }
    else test_pass "`testname'"
end

capture program drop _cfx_cpar
program define _cfx_cpar
    * matchit and cfuzzyjoin match on two columns: the same scores
    syntax varlist(min=2 max=2), testname(string) [opts(string asis) tol(real 0)]
    tempfile data ref
    quietly save `data'
    capture matchit `varlist', `opts' gen(__s1)
    local rc1 = _rc
    if `rc1' {
        quietly use `data', clear
        capture cfuzzyjoin match `varlist', `opts' gen(__s2)
        if _rc == `rc1' test_pass "`testname'"
        else test_fail "`testname'" "matchit r(`rc1'), cfuzzyjoin r(`=_rc')"
        quietly use `data', clear
        exit
    }
    keep __s1
    quietly save `ref'
    quietly use `data', clear
    capture cfuzzyjoin match `varlist', `opts' gen(__s2)
    if _rc {
        test_fail "`testname'" "cfuzzyjoin r(`=_rc')"
        quietly use `data', clear
        exit
    }
    quietly merge 1:1 _n using `ref', nogen
    if `tol' == 0 capture assert __s1 == __s2
    else capture assert reldif(__s1, __s2) <= `tol'
    if _rc {
        quietly count if !(reldif(__s1, __s2) <= `tol')
        test_fail "`testname'" "`r(N)' scores differ from matchit"
    }
    else test_pass "`testname'"
    quietly use `data', clear
end

capture program drop _cfx_ipar
program define _cfx_ipar
    * matchit and cfuzzyjoin match with using: the same return code, pairs,
    * scores, variable names, storage types and formats (rows sorted on all
    * variables, so duplicated identifiers are compared too)
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
        local N1 = _N
        local sortv : word 1 of `v1'
        local sortu : word 3 of `v1'
        if `tol' == 0 sort _all
        else sort `sortv' `sortu'
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
    if "`v1'" != "`v2'" | "`t1'" != "`t2'" | `N1' != _N {
        test_fail "`testname'" "variables [`v1'] [`t1'] N=`N1' vs [`v2'] [`t2'] N=`=_N'"
        quietly use `data', clear
        exit
    }
    if `tol' == 0 {
        sort _all
        capture cf _all using `ref'
        local rc = _rc
    }
    else {
        sort `sortv' `sortu'
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

capture program drop _cfx_lpar
program define _cfx_lpar
    * reclink and cfuzzyjoin link: the same return code, variables (order,
    * types, formats, labels) and observations
    syntax anything(name=vars), using(string) testname(string) opts(string asis) ///
        idm(name) idu(name)
    tempfile data ref
    quietly save `data'
    capture reclink `vars' using `"`using'"', idm(`idm') idu(`idu') `opts'
    local rc1 = _rc
    if !`rc1' {
        _cfx_describe
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
    _cfx_describe
    local d2 `"`r(desc)'"'
    if `"`d1'"' != `"`d2'"' {
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

capture program drop _cfx_child
program define _cfx_child
    * _cfx_child check testname: runs a check of the child section above in a
    * separate Stata process (same executable, working directory and ado-path)
    args check testname
    local exe : env CTOOLS_STATA_EXECUTABLE
    if `"`exe'"' == "" {
        local fl = cond(c(MP), "MP", cond(c(SE), "SE", "BE"))
        if c(os) == "Windows" local exe `"`c(sysdir_stata)'Stata`fl'-64.exe"'
        else local exe `"`c(sysdir_stata)'stata-`=lower("`fl'")'"'
    }
    capture confirm file `"`exe'"'
    if _rc {
        test_fail "`testname'" `"no Stata executable to run the check (`exe')"'
        exit
    }
    tempfile res
    tempname fh
    local tmp = c(tmpdir)
    if inlist(substr("`tmp'", -1, 1), "/", "\") local tmp = substr("`tmp'", 1, length("`tmp'") - 1)
    local wrap `"`tmp'/cfx_child_`check'.do"'
    file open `fh' using `"`wrap'"', write text replace
    file write `fh' `"global S_ADO `"`c(adopath)'"'"' _n
    file write `fh' `"cd "`c(pwd)'""' _n
    * the path sits in a macro so the runner does not take this call for a
    * nested case (validation/stata.py cases()), which would hide the case
    local self "validation/stata/cfuzzyjoin.do"
    file write `fh' `"do "`self'" cfuzzyjoin_extra child `check' "`res'""' _n
    file close `fh'
    if c(os) == "Windows" quietly shell cd /d "`tmp'" && "`exe'" /e do "`wrap'"
    else quietly shell cd "`tmp'" && "`exe'" -b do "`wrap'"
    capture erase `"`wrap'"'
    capture erase `"`tmp'/cfx_child_`check'.log"'
    capture confirm file `"`res'"'
    if _rc {
        test_fail "`testname'" "Stata stopped (crashed) before the check finished"
        exit
    }
    file open `fh' using `"`res'"', read text
    file read `fh' line
    file close `fh'
    gettoken ok why : line
    if "`ok'" == "1" test_pass "`testname'"
    else test_fail "`testname'" `"`why'"'
end

capture program drop _cfx_describe
program define _cfx_describe, rclass
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

* =============================================================================
* group: strgroup
* =============================================================================
* tiny data: nothing selected, one string, all empty, all equal
* (strgroup stops with r(198) when if/in selects nothing, from its warning
* about duplicates, which reads results that duplicates report did not leave:
* the contract is that of all-empty strings, a missing group for each)
_cfx_rows s : "abc" "abd" "xyz"
capture cfuzzyjoin group s if s == "none", gen(g) threshold(0.4)
local rc = _rc
capture assert missing(g)
if !`rc' & !_rc test_pass "extra group: no observation selected gives missing groups"
else test_fail "extra group: no observation selected gives missing groups" "r(`rc')"
_cfx_child group "extra group without observations"
_cfx_rows s : "abc"
_cfx_gpar s, testname("extra group: one observation") opts(threshold(0.4))
_cfx_rows s : "" "" ""
_cfx_gpar s, testname("extra group: all strings empty") opts(threshold(0.4))
_cfx_rows s : "same" "same" "same" "same"
_cfx_gpar s, testname("extra group: identical strings") opts(threshold(0.01))
quietly cfuzzyjoin group s, gen(g) threshold(0.01)
capture assert g == 1
if _rc == 0 & r(groups) == 1 & r(distinct) == 1 & r(N) == 4 test_pass "extra group: identical strings form one group"
else test_fail "extra group: identical strings form one group" "r(groups)=`r(groups)'"

* case, punctuation and unicode (distances count bytes)
_cfx_rows s : "Smith" "smith" "SMITH" "smiTh" "Smyth" "O'Brien" "OBrien" "o brien" "A.B.C." "ABC" ///
    "a-b-c" "abc" "café" "cafe" "cafè" "CAFÉ" "naïve" "naive" "Zürich" "Zurich" "東京" "東京都" ///
    "Straße" "Strasse" "ñandú" "nandu" "😀 smile" "😀 smiles" "x,y;z" "x y z"
tempfile gu
quietly save `gu'
foreach o in "threshold(0.2)" "threshold(0.34)" "threshold(0.5) normalize(longer)" ///
    "threshold(2) normalize(none)" "threshold(0.5) first" {
    quietly use `gu', clear
    _cfx_gpar s, testname("extra group: case, punctuation, unicode: `o'") opts(`o')
}
* first: strings that start with digits, punctuation, blanks, upper case and
* bytes above 127
_cfx_rows s : "1abc" "1abd" "_x" "-x" "Abc" "abc" "Zed" "zed" "zee" "éa" "éb" "[a]" "a" " a" " b"
_cfx_gpar s, testname("extra group: first with non-letters") opts(threshold(0.5) first noclean)

* very long strings (the banded distance), up to 2,045 bytes
clear
quietly set obs 12
quietly gen str2045 s = substr("abcdefghij" * 205, 1, 2000)
quietly replace s = substr(s, 1, 10 * _n - 1) + "X" + substr(s, 10 * _n + 1, .) if inrange(_n, 2, 6)
quietly replace s = substr(s, 1, 1999) if _n == 7
quietly replace s = substr(s, 1, 500) + "YY" + substr(s, 503, .) if _n == 8
quietly replace s = s + "Z" if _n == 9
quietly replace s = "Q" + substr(s, 2, .) if _n == 10
quietly replace s = substr("abcdefghij" * 205, 1, 2045) if _n >= 11
quietly replace s = substr(s, 1, 2044) + "W" if _n == 12
* strgroup's plugin reads only the first 244 bytes of each string, so these
* are checked against the groups of the exact distances: rows 2-6, 9 and 10
* are 1 edit from row 1 (2,000 bytes), row 7 is row 1 less its last byte
* (1/1999 > 0.0005), row 8 is 2 edits from row 1, rows 11 and 12 (2,045
* bytes) are 1 edit apart and 45 from row 1
foreach spec in "threshold(0.0005)|1 1 1 1 1 1 2 3 1 1 4 4" ///
    "threshold(0.0004999)|1 2 3 4 5 6 7 8 9 10 11 11" ///
    "threshold(0.001) normalize(longer)|1 1 1 1 1 1 1 1 1 1 2 2" ///
    "threshold(1) normalize(none)|1 1 1 1 1 1 1 2 1 1 3 3" ///
    "threshold(45) normalize(none)|1 1 1 1 1 1 1 1 1 1 1 1" ///
    "threshold(44) normalize(none)|1 1 1 1 1 1 1 1 1 1 2 2" {
    gettoken o want : spec, parse("|")
    local want = substr("`want'", 2, .)
    capture drop g
    quietly cfuzzyjoin group s, gen(g) `o'
    local got
    forvalues i = 1/12 {
        local got `got' `=g[`i']'
    }
    if "`got'" == "`want'" test_pass "extra group: long strings, exact distances: `o'"
    else test_fail "extra group: long strings, exact distances: `o'" "groups `got', expected `want'"
}
* strL: the same groups as the str# copy
quietly use `gu', clear
quietly cfuzzyjoin group s, gen(g1) threshold(0.34)
quietly recast strL s
capture noisily cfuzzyjoin group s, gen(g2) threshold(0.34)
if _rc test_fail "extra group: strL equals str#" "r(`=_rc')"
else {
    capture assert g1 == g2
    if _rc test_fail "extra group: strL equals str#" "groups differ"
    else test_pass "extra group: strL equals str#"
}
* strL beyond 2,045 bytes: compared on the first 2,045 bytes (documented)
clear
quietly set obs 3
quietly gen strL s = "abcdefghij" * 300
quietly replace s = s + "tail" in 2
quietly replace s = "Q" + s in 3
capture noisily cfuzzyjoin group s, gen(g) threshold(0.0001)
if _rc test_fail "extra group: strL compared on 2,045 bytes" "r(`=_rc')"
else {
    capture assert g[1] == g[2] & g[3] != g[1]
    if _rc test_fail "extra group: strL compared on 2,045 bytes" "groups `=g[1]' `=g[2]' `=g[3]'"
    else test_pass "extra group: strL compared on 2,045 bytes"
}

* thresholds at the cutoff (normalize(longer): 1/4 and 1/8; none: 1 edit;
* shorter: 1/3)
_cfx_rows s : "abcd" "abce" "wxyz" "wxy" "pqrstuvw" "pqrstuvx" "klm" "klmn"
tempfile gb
quietly save `gb'
foreach t in 0.25 0.2499999 0.24999995 0.24999994 0.2500001 0.125 0.1249999 0.12499995 0.12499994 {
    quietly use `gb', clear
    _cfx_gpar s, testname("extra group: longer at the cutoff `t'") opts(threshold(`t') normalize(longer))
}
foreach t in 1 0.9999999 0.99999995 0.99999994 1.0000001 {
    quietly use `gb', clear
    _cfx_gpar s, testname("extra group: none at the cutoff `t'") opts(threshold(`t') normalize(none))
}
foreach t in 0.3333333 0.33333333 0.3333334 {
    quietly use `gb', clear
    _cfx_gpar s, testname("extra group: shorter at the cutoff `t'") opts(threshold(`t'))
}

* random data over several sizes and seeds
local gopts `""threshold(0.2)" "threshold(0.4) normalize(longer)" "threshold(1) normalize(none) first""'
foreach spec in "1 7201 1" "7 7202 2" "60 7203 3" "400 7204 1" "3000 7205 1 2 3" {
    gettoken n spec : spec
    gettoken seed which : spec
    _cfx_names `n' `seed'
    tempfile gr
    quietly save `gr'
    foreach j of local which {
        local o : word `j' of `gopts'
        quietly use `gr', clear
        _cfx_gpar name, testname("extra group: random n=`n': `o'") opts(`o')
    }
}
* shuffled order: groups are numbered by first observation
_cfx_names 2000 7206
set seed 7207
quietly gen double u = runiform()
sort u
drop u
_cfx_gpar name, testname("extra group: shuffled order") opts(threshold(0.25))
* by with if, an all-empty by-group and a one-string by-group
_cfx_names 1500 7208
quietly gen byte grp = mod(id, 3)
quietly replace grp = 9 if id <= 20
quietly replace name = "" if grp == 9
quietly replace grp = 8 if id == 21
sort grp id
_cfx_gpar name if id > 10, by(grp) testname("extra group: by with if, empty and single by-groups") ///
    opts(threshold(0.3))
* abbreviations and a threshold larger than any distance
_cfx_names 500 7209
tempfile ga
quietly save `ga'
_cfx_gpar name, testname("extra group: abbreviated options") opts(thresh(.3) norm(Longer))
quietly use `ga', clear
_cfx_gpar name, testname("extra group: threshold above every distance") opts(threshold(1000))
* errors leave the data as they were
quietly use `ga', clear
quietly datasignature
local sig `r(datasignature)'
foreach o in "threshold(-1)" "threshold(0.2) normalize(bad)" "threshold(0.2) threads(999)" "" {
    capture cfuzzyjoin group name, gen(g) `o'
    local rc = _rc
    quietly datasignature
    local sig2 `r(datasignature)'
    capture confirm variable g
    if `rc' & _rc & "`sig2'" == "`sig'" test_pass "extra group: error leaves the data: `o'"
    else test_fail "extra group: error leaves the data: `o'" "r(`rc')"
}
* stored results: strings, distinct strings, groups numbered 1..r(groups)
quietly use `ga', clear
quietly cfuzzyjoin group name, gen(g) threshold(0.3)
local N = r(N)
local d = r(distinct)
local G = r(groups)
quietly count if name != ""
local eN = r(N)
quietly egen byte __t = tag(name) if name != ""
quietly count if __t == 1
local ed = r(N)
quietly egen byte __tg = tag(g) if !missing(g)
quietly count if __tg == 1
local eg = r(N)
quietly summarize g
if `N' == `eN' & `d' == `ed' & `G' == `eg' & r(max) == `G' & r(min) == 1 test_pass "extra group: r(N), r(distinct), r(groups)"
else test_fail "extra group: r(N), r(distinct), r(groups)" "`N'/`eN' `d'/`ed' `G'/`eg'"
* names of temporary variables already in the data (Stata's tempvar does not
* skip them: both commands stop with r(110))
_cfx_rows s : "abc" "abd" "xyz"
quietly gen byte __000000 = 7
_cfx_gpar s, testname("extra group: temporary names in the data, as strgroup") opts(threshold(0.4))
* threads on long strings and a large dataset
_cfx_names 3000 7210
quietly replace name = name + " " + name + " " + name + " " + name + " " + name if name != ""
quietly cfuzzyjoin group name, gen(g1) threshold(0.15) threads(1)
quietly cfuzzyjoin group name, gen(g2) threshold(0.15) threads(16)
quietly cfuzzyjoin group name, gen(g3) threshold(3) normalize(none) threads(1)
quietly cfuzzyjoin group name, gen(g4) threshold(3) normalize(none) threads(16)
capture assert g1 == g2 & g3 == g4
if _rc test_fail "extra group: threads(1) equals threads(16), long strings" ""
else test_pass "extra group: threads(1) equals threads(16), long strings"
_cfx_names 100000 7211
quietly cfuzzyjoin group name, gen(g1) threshold(0.2) threads(1)
quietly cfuzzyjoin group name, gen(g2) threshold(0.2)
local G = r(groups)
quietly summarize g2
capture assert g1 == g2
if !_rc & r(max) == `G' test_pass "extra group: 100,000 strings, threads give the same groups"
else test_fail "extra group: 100,000 strings, threads give the same groups" ""

* =============================================================================
* match: matchit, two columns
* =============================================================================
_cfx_rows s1 s2 : "Zürich" "Zurich" "東京" "東京都" "naïve café" "naive cafe" "" "abc" "abc" "" "" "" ///
    "ÉCOLE" "école" "😀 smile" "😀 smiles" "straße" "STRASSE" "same text" "same text" "Ab" "ab" "a" "a" ///
    "O'Brien & Co." "OBrien and Co" "x,y;z" "x y z" "  lead" "lead  " "MÜLLER" "Mueller" ///
    "the the the" "the" "aaaa" "aa"
tempfile cu
quietly save `cu'
foreach sim in "bigram" "ngram, 3" "ngram_circ, 3" "token" "cotoken" "scotoken" "firstgram, 2" ///
    "tkngram, 2" "soundex" "soundex_nara" "soundex_fk" "soundex_ext" "nysiis_fk" "token_soundex" ///
    `"tokenwrap, "nysiis_fk""' {
    quietly use `cu', clear
    _cfx_cpar s1 s2, testname(`"extra match columns unicode, case, punctuation: `sim' minsimple"') ///
        opts(sim(`sim') score(minsimple))
}
quietly use `cu', clear
quietly cfuzzyjoin match s1 s2, gen(b)
quietly cfuzzyjoin match s1 s2, gen(t) sim(token)
capture assert b == 1 & t == 1 if s1 == s2 & length(s1) > 1
if _rc test_fail "extra match columns: identical strings score 1" ""
else test_pass "extra match columns: identical strings score 1"
* strL columns
quietly use `cu', clear
quietly recast strL s1 s2
foreach sim in "bigram" "token" {
    _cfx_cpar s1 s2, testname("extra match columns strL: `sim'") opts(sim(`sim'))
}
* strL values beyond 2,045 bytes: compared on their first 2,045 (documented)
clear
quietly set obs 2
quietly gen strL s1 = "abcdefghij" * 300
quietly gen strL s2 = s1 + "tail"
quietly replace s2 = "zz" in 2
capture noisily cfuzzyjoin match s1 s2, gen(sc)
if _rc test_fail "extra match columns: strL compared on 2,045 bytes" "r(`=_rc')"
else {
    capture assert sc[1] == 1 & sc[2] == 0
    if _rc test_fail "extra match columns: strL compared on 2,045 bytes" "scores `=sc[1]' `=sc[2]'"
    else test_pass "extra match columns: strL compared on 2,045 bytes"
}
* all strings empty
clear
quietly set obs 3
quietly gen str1 s1 = ""
quietly gen str1 s2 = ""
foreach sim in "bigram" "token" "soundex" {
    _cfx_cpar s1 s2, testname("extra match columns all empty: `sim'") opts(sim(`sim'))
}
* very long strings
clear
quietly set obs 6
quietly gen str2045 s1 = substr("lorem ipsum dolor sit amet " * 80, 1, 2000)
quietly gen str2045 s2 = s1
quietly replace s2 = substr(s2, 1, 999) + "Q" + substr(s2, 1001, .) in 2
quietly replace s2 = upper(s2) in 3
quietly replace s2 = substr(s2, 1, 1500) in 4
quietly replace s2 = "lorem" in 5
quietly replace s1 = substr("abcdefghijklmnopqrstuvwxyz" * 80, 1, 2045) in 6
quietly replace s2 = substr("bcdefghijklmnopqrstuvwxyza" * 80, 1, 2045) in 6
tempfile clong
quietly save `clong'
foreach sim in "bigram" "ngram, 4" "token" "soundex_ext" "tkngram, 3" {
    quietly use `clong', clear
    _cfx_cpar s1 s2, testname("extra match columns long strings: `sim'") opts(sim(`sim'))
}
* gram sizes, weights with every score, abbreviations
_cfx_names 800 7301
rename name name1
tempfile c1
quietly save `c1'
_cfx_names 800 7302
rename name name2
quietly merge 1:1 id using `c1', nogen
tempfile cols
quietly save `cols'
foreach sim in "ngram, 4" "ngram, 6" "ngram_circ, 2" "ngram_circ, 5" "firstgram, 1" "firstgram, 4" ///
    "tkngram, 1" "tkngram, 4" {
    quietly use `cols', clear
    _cfx_cpar name1 name2, testname("extra match columns gram size: `sim'") opts(sim(`sim'))
}
foreach o in "weights(log) score(simple)" "sim(token) weights(root) score(minsimple)" ///
    "sim(ngram, 4) weights(simple) score(minsimple)" "sim(soundex_fk) weights(log) score(simple)" ///
    "sim(tkngram, 2) weights(root)" "sim(token) s(minsimple) w(log)" {
    quietly use `cols', clear
    _cfx_cpar name1 name2, testname("extra match columns: `o'") opts(`o') tol(1e-13)
}
quietly use `cols', clear
capture matchit name1 name2, gen(name1)
local rc1 = _rc
capture cfuzzyjoin match name1 name2, gen(name1)
if _rc == `rc1' & `rc1' test_pass "extra match columns: existing score variable refused"
else test_fail "extra match columns: existing score variable refused" "r(`=_rc') vs matchit r(`rc1')"

* =============================================================================
* match: matchit, two datasets
* =============================================================================
* scores at the threshold: "abc"/"abd" 0.5 (jaccard), "abc"/"abcd" 0.8
* (simple), "abc"/"abcde" 1 (minsimple)
_cfx_rows name : "abc" "abcd" "wxyz"
quietly gen long id = _n
tempfile mb ub
quietly save `mb'
_cfx_rows uname : "abd" "abc" "abcde" "wxyy" "zzzz"
quietly gen long uid = _n + 10
quietly save `ub'
foreach o in "threshold(0.5)" "threshold(0.5000001)" "threshold(0.4999999)" "score(simple) threshold(0.8)" ///
    "score(simple) threshold(0.8000001)" "score(minsimple) threshold(1)" {
    quietly use `mb', clear
    _cfx_ipar id name, using(`ub') testname("extra match at the threshold: `o'") opts(idu(uid) txtu(uname) `o')
}
* many-to-many: duplicated strings and identifiers on both sides, empty strings
_cfx_rows name : "acme corp" "acme corp" "acme co" "beta" "beta" "" "gamma"
quietly gen long id = _n
quietly replace id = 1 in 2
tempfile mm uu
quietly save `mm'
_cfx_rows uname : "acme corp" "acme corp" "acme corporation" "beta" "betas" "" "beta"
quietly gen long uid = _n
quietly replace uid = 2 in 2
quietly replace uid = 4 in 7
quietly save `uu'
foreach o in "threshold(0)" "threshold(0.5)" "sim(token) threshold(0)" {
    quietly use `mm', clear
    _cfx_ipar id name, using(`uu') testname("extra match many-to-many: `o'") opts(idu(uid) txtu(uname) `o')
}
* one record each; no master records and no using records (in a separate
* Stata process, see _cfx_child)
_cfx_rows name : "hello world"
quietly gen long id = 1
tempfile m1 u1
quietly save `m1'
_cfx_rows uname : "hello word"
quietly gen long uid = 1
quietly save `u1'
quietly use `m1', clear
_cfx_ipar id name, using(`u1') testname("extra match one record each") opts(idu(uid) txtu(uname))
_cfx_rows name : "" "" ""
quietly gen long id = _n
tempfile me
quietly save `me'
_cfx_ipar id name, using(`u1') testname("extra match all master strings empty") opts(idu(uid) txtu(uname))
quietly use `me', clear
rename (id name) (uid uname)
tempfile ue
quietly save `ue'
quietly use `m1', clear
_cfx_ipar id name, using(`ue') testname("extra match all using strings empty") opts(idu(uid) txtu(uname))
_cfx_child match_master "extra match without master records equals matchit"
_cfx_child match_using "extra match without using records equals matchit"
_cfx_child columns "extra match columns without observations"
* unicode, case and punctuation
quietly use `cu', clear
quietly gen long id = _n
keep id s1
tempfile mu uu2
quietly save `mu'
quietly use `cu', clear
quietly gen long uid = _n
keep uid s2
quietly save `uu2'
foreach sim in "bigram" "token" "soundex" "nysiis_fk" "ngram_circ, 3" "tkngram, 2" {
    quietly use `mu', clear
    _cfx_ipar id s1, using(`uu2') testname("extra match unicode, case, punctuation: `sim'") ///
        opts(idu(uid) txtu(s2) sim(`sim') threshold(0))
}
* strL text variables (the result holds str#)
quietly use `cu', clear
quietly gen long uid = _n
keep uid s2
quietly recast strL s2
tempfile uul
quietly save `uul'
quietly use `mu', clear
quietly recast strL s1
_cfx_ipar id s1, using(`uul') testname("extra match strL text variables") opts(idu(uid) txtu(s2) threshold(0))
* long strings
quietly use `clong', clear
quietly gen long id = _n
keep id s1
tempfile ml ul
quietly save `ml'
quietly use `clong', clear
quietly gen long uid = _n
keep uid s2
quietly save `ul'
foreach sim in "bigram" "token" {
    quietly use `ml', clear
    _cfx_ipar id s1, using(`ul') testname("extra match long strings: `sim'") opts(idu(uid) txtu(s2) sim(`sim') threshold(0))
}
* other similarity functions and scores on random names
_cfx_names 300 7303
keep id name
tempfile mr ur
quietly save `mr'
_cfx_names 300 7304
keep id name
rename (id name) (uid uname)
quietly save `ur'
foreach o in "sim(tkngram, 2) threshold(0.3)" "sim(firstgram, 3) score(simple) threshold(0.3)" ///
    "sim(token_soundex) score(minsimple)" "sim(soundex_ext) threshold(0.6)" "sim(cotoken) threshold(0.2)" ///
    "sim(ngram_circ, 4) score(simple) threshold(0.4)" "sim(ngram, 1) score(minsimple) threshold(0.9)" {
    quietly use `mr', clear
    _cfx_ipar id name, using(`ur') testname("extra match random names: `o'") opts(idu(uid) txtu(uname) `o')
}
* weights and stopwords with other scores (matchit's weighted path is slow:
* few records)
quietly use `ur', clear
quietly keep in 1/150
tempfile urs
quietly save `urs'
foreach o in "weights(log) score(minsimple) threshold(0.4)" "weights(root) score(simple) sim(token) threshold(0.3)" ///
    "sim(token) stopwordsauto swthreshold(0.3) threshold(0.3)" {
    quietly use `mr', clear
    quietly keep in 1/150
    _cfx_ipar id name, using(`urs') testname("extra match: `o'") opts(idu(uid) txtu(uname) `o') tol(1e-13)
}
* abbreviated options and a float identifier that is not an integer
quietly use `mr', clear
_cfx_ipar id name, using(`ur') testname("extra match abbreviated options") ///
    opts(idu(uid) txtu(uname) sim(token) s(minsimple) t(0.4) gen(sc))
quietly use `mr', clear
quietly gen float fid = id + 0.5
drop id
order fid
_cfx_ipar fid name, using(`ur') testname("extra match non-integer float identifier") opts(idu(uid) txtu(uname))
* a string identifier is refused as by matchit
quietly use `mr', clear
quietly tostring id, replace
_cfx_ipar id name, using(`ur') testname("extra match string identifier refused") opts(idu(uid) txtu(uname))
* the result does not depend on the order of the master data
quietly use `mr', clear
quietly cfuzzyjoin match id name using `ur', idu(uid) txtu(uname) threshold(0.4)
sort _all
tempfile o1
quietly save `o1'
quietly use `mr', clear
set seed 7305
quietly gen double u = runiform()
sort u
drop u
quietly cfuzzyjoin match id name using `ur', idu(uid) txtu(uname) threshold(0.4) override
sort _all
capture cf _all using `o1'
if _rc test_fail "extra match: master order does not change the pairs" "cf r(`=_rc')"
else test_pass "extra match: master order does not change the pairs"
* errors leave the data as they were
quietly use `mr', clear
quietly datasignature
local sig `r(datasignature)'
foreach o in "idu(uid) txtu(uname) sim(nosuch)" "idu(uid) txtu(uname) sim(ngram, 0)" "idu(uid) txtu(nosuch)" ///
    "idu(uid) txtu(uname) score(bad)" "idu(uid) txtu(uname) weights(log) wgtfile(nosuchfile)" {
    capture cfuzzyjoin match id name using `ur', `o' override
    local rc = _rc
    quietly datasignature
    if `rc' & "`r(datasignature)'" == "`sig'" test_pass "extra match: error leaves the data: `o'"
    else test_fail "extra match: error leaves the data: `o'" "r(`rc')"
}
* the index scores equal the column scores of the same pairs, and threads
* give the same pairs (larger data)
_cfx_names 5000 7306
keep id name
tempfile mbig ubig
quietly save `mbig'
_cfx_names 5000 7307
keep id name
rename (id name) (uid uname)
quietly save `ubig'
quietly use `mbig', clear
quietly cfuzzyjoin match id name using `ubig', idu(uid) txtu(uname) threshold(0.7) threads(1)
tempfile b1
quietly save `b1'
quietly cfuzzyjoin match name uname, gen(chk)
capture assert chk == similscore
local a = _rc
quietly use `mbig', clear
quietly cfuzzyjoin match id name using `ubig', idu(uid) txtu(uname) threshold(0.7) threads(16)
capture cf _all using `b1'
if !`a' & !_rc test_pass "extra match: index scores equal column scores; threads(1) equals threads(16)"
else test_fail "extra match: index scores equal column scores; threads(1) equals threads(16)" "columns r(`a'), cf r(`=_rc')"

* =============================================================================
* link: reclink
* =============================================================================
_cfx_names 300 7401
keep id name
quietly bysort name (id): keep if _n == 1
sort id
quietly gen str12 city = word("boston newton cambridge", 1 + mod(id, 3))
quietly replace city = upper(city) if mod(id, 11) == 0
quietly gen double zip = 2100 + mod(id, 13)
quietly replace zip = .a if mod(id, 17) == 0
quietly replace zip = .b if mod(id, 19) == 0
tempfile xm xu
quietly save `xm'
_cfx_names 400 7402
keep id name
rename id uid
quietly bysort name (uid): keep if _n == 1
sort uid
quietly gen str12 city = word("boston newton cambridge", 1 + mod(uid, 3))
quietly gen double zip = 2100 + mod(uid, 13)
quietly replace zip = .a if mod(uid, 7) == 0
quietly replace zip = .b if mod(uid, 23) == 0
quietly save `xu'
foreach spec in "name|minscore(0.3)" "name city|minbigram(0.3)" "name city|minbigram(0.95)" ///
    "name city zip|wmatch(1.5 3.25 2) wnomatch(0.25 1 0.5)" "name city|exactstr(city)" ///
    "name city zip|required(zip)" "name city zip|orblock(zip city)" "name city|minscore(0.99999)" {
    gettoken vars opts : spec, parse("|")
    local opts = substr(`"`opts'"', 2, .)
    quietly use `xm', clear
    _cfx_lpar `vars', using(`xu') idm(id) idu(uid) testname("extra link: `spec'") opts(gen(score) `opts')
}
* scores at minscore(): the lowest fuzzy score of a first run as the cutoff
quietly use `xm', clear
quietly reclink name city using `xu', idm(id) idu(uid) gen(score) minscore(0.3)
quietly summarize score if _merge == 3 & score < 1
local s : di %21.17g r(min)
local s2 : di %21.17g r(min) - 1e-7
foreach m in `s' `s2' {
    quietly use `xm', clear
    _cfx_lpar name city, using(`xu') idm(id) idu(uid) testname("extra link: minscore at a score (`m')") ///
        opts(gen(score) minscore(`m'))
}
* identical datasets: every record a perfect match
quietly use `xm', clear
rename id uid
tempfile xsame
quietly save `xsame'
quietly use `xm', clear
_cfx_lpar name city zip, using(`xsame') idm(id) idu(uid) testname("extra link: identical datasets") opts(gen(score))
quietly use `xm', clear
quietly cfuzzyjoin link name city zip using `xsame', idm(id) idu(uid) gen(score)
local nm = r(matched)
local np = r(perfect)
capture assert score == 1 & _merge == 3 & uid == id
if !_rc & `nm' == _N & `np' == _N test_pass "extra link: identical datasets are perfect matches"
else test_fail "extra link: identical datasets are perfect matches" "matched `nm', perfect `np' of `=_N'"
* one master record; no master records, no using records (separate process)
quietly use `xm', clear
quietly keep in 5
_cfx_lpar name city, using(`xu') idm(id) idu(uid) testname("extra link: one master record") opts(gen(score) minscore(0.3))
quietly use `xm', clear
quietly replace name = ""
_cfx_lpar name city, using(`xu') idm(id) idu(uid) testname("extra link: all master names empty") opts(gen(score) minscore(0.3))
_cfx_child link_master "extra link without master records equals reclink"
_cfx_child link_using "extra link without using records equals reclink"
* unicode, case and punctuation
_cfx_rows name city : "Zürich Versicherung" "Zürich" "MÜLLER GmbH" "Köln" "O'Brien & Co." "Dublin" ///
    "café de flore" "Paris" "東京電力" "東京" "Straße 5" "Berlin"
quietly gen long id = _n
tempfile um uu3
quietly save `um'
_cfx_rows name city : "Zurich Versicherung" "Zürich" "Mueller GmbH" "Köln" "OBrien and Co" "Dublin" ///
    "Cafe de Flore" "Paris" "東京電力株式会社" "東京" "Strasse 5" "Berlin" "zürich versicherung" "Zürich"
quietly gen long uid = _n + 100
quietly save `uu3'
foreach o in "minscore(0.2)" "minscore(0.2) required(city)" {
    quietly use `um', clear
    _cfx_lpar name city, using(`uu3') idm(id) idu(uid) testname("extra link: unicode, case, punctuation: `o'") ///
        opts(gen(score) `o')
}
* long strings
clear
quietly set obs 4
quietly gen str700 name = substr("acme widget holdings international " * 20, 1, 650)
quietly replace name = name + string(_n)
quietly gen long id = _n
tempfile lgm lgu
quietly save `lgm'
quietly replace name = substr(name, 1, 300) + "X" + substr(name, 302, .) if _n == 2
quietly replace name = upper(name) if _n == 3
quietly rename id uid
quietly replace uid = uid + 10
quietly save `lgu'
quietly use `lgm', clear
_cfx_lpar name, using(`lgu') idm(id) idu(uid) testname("extra link: long strings") opts(gen(score) minscore(0.1))
* a using variable the master has and a strL using variable; a prefixed name
* the using data already has (reclink stops with r(110) on its rename)
quietly use `xu', clear
quietly gen strL note = "note " + string(uid)
tempfile xu1 xu2 xu3
quietly save `xu1'
quietly gen str6 extra = "u" + string(uid)
quietly save `xu2'
quietly gen str5 Uname = "taken"
quietly save `xu3'
quietly use `xm', clear
_cfx_lpar name city, using(`xu1') idm(id) idu(uid) testname("extra link: strL using variable") ///
    opts(gen(score) minscore(0.3))
quietly use `xm', clear
quietly gen str2 extra = "m"
_cfx_lpar name city, using(`xu2') idm(id) idu(uid) testname("extra link: using variable the master has (wider)") ///
    opts(gen(score) minscore(0.3))
quietly use `xm', clear
_cfx_lpar name city, using(`xu3') idm(id) idu(uid) testname("extra link: prefixed name taken in the using data") ///
    opts(gen(score) minscore(0.3))
* the order of the master and using data does not change the result
quietly use `xm', clear
quietly cfuzzyjoin link name city zip using `xu', idm(id) idu(uid) gen(score) minscore(0.3)
sort id uid
tempfile r0
quietly save `r0'
quietly use `xu', clear
set seed 7403
quietly gen double u = runiform()
sort u
drop u
tempfile xush
quietly save `xush'
quietly use `xm', clear
quietly gen double u = runiform()
sort u
drop u
quietly cfuzzyjoin link name city zip using `xush', idm(id) idu(uid) gen(score) minscore(0.3)
sort id uid
capture cf _all using `r0'
if _rc test_fail "extra link: data order does not change the result" "cf r(`=_rc')"
else test_pass "extra link: data order does not change the result"
* stored results
quietly use `xm', clear
quietly cfuzzyjoin link name city using `xu', idm(id) idu(uid) gen(score) minscore(0.3)
local rN = r(N)
local rm = r(matched)
local ru = r(unmatched)
quietly egen byte __t = tag(id) if _merge == 3
quietly count if __t == 1
local em = r(N)
quietly count if _merge == 1
if `rN' == _N & `rm' == `em' & `ru' == r(N) test_pass "extra link: r(N), r(matched), r(unmatched)"
else test_fail "extra link: r(N), r(matched), r(unmatched)" "`rN' `rm'/`em' `ru'/`=r(N)'"
* errors leave the data as they were
quietly use `xm', clear
quietly datasignature
local sig `r(datasignature)'
foreach o in "gen(score) wmatch(0.5 1)" "gen(score) exclude(nosuchfile)" "gen(score) uvarlist(name nosuch)" ///
    "gen(score) required(zip)" "gen(score) threads(999)" {
    capture cfuzzyjoin link name city using `xu', idm(id) idu(uid) `o'
    local rc = _rc
    quietly datasignature
    if `rc' & "`r(datasignature)'" == "`sig'" test_pass "extra link: error leaves the data: `o'"
    else test_fail "extra link: error leaves the data: `o'" "r(`rc')"
}
* threads on larger data with or-blocking
_cfx_names 4000 7404
keep id name
quietly gen str12 city = word("boston newton cambridge somerville", 1 + mod(id, 4))
quietly gen zip = 2100 + mod(id, 37)
tempfile bm bu
quietly save `bm'
_cfx_names 4000 7405
keep id name
rename id uid
quietly gen str12 city = word("boston newton cambridge somerville", 1 + mod(uid, 4))
quietly gen zip = 2100 + mod(uid, 29)
quietly save `bu'
quietly use `bm', clear
quietly cfuzzyjoin link name city zip using `bu', idm(id) idu(uid) gen(score) orblock(city zip) threads(1)
tempfile t1
quietly save `t1'
quietly use `bm', clear
quietly cfuzzyjoin link name city zip using `bu', idm(id) idu(uid) gen(score) orblock(city zip) threads(16)
capture cf _all using `t1'
if _rc test_fail "extra link: threads(1) equals threads(16), or-blocking" "cf r(`=_rc')"
else test_pass "extra link: threads(1) equals threads(16), or-blocking"

print_summary "cfuzzyjoin_extra"
global CTOOLS_COMPONENT_COMPLETE "cfuzzyjoin_extra"
    global CTOOLS_CASE_COMPLETE "cfuzzyjoin_extra"
    exit
}

di as error "Unknown cfuzzyjoin case: `__case'"
exit 198
