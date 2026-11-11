* s2_iv: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_iv"

if "`__case'" == "s2_iv" {
do validation/support/validate_setup.do
set seed 20261104

/* civreghdfe against ivreg2 (no absorb()) and ivreghdfe (absorb()) for the
   S2 gaps IV-01 (replay), IV-04 (center), IV-05 (sfirst), IV-06 (bw(auto)),
   IV-07 (kernels), IV-08 (smatrix/wmatrix), IV-10 (orthog() on included
   instruments), IV-25/IV-35 (first-stage and weak-instrument statistics).
   Tolerance: mreldif 1e-7 on e(b), e(V), every e() statistic compared and
   e(first); the displayed blocks must be identical lines. */

* s2iv_ts: a time series with autocorrelated errors (tsset t)
capture program drop s2iv_ts
program define s2iv_ts
    clear
    set obs 150
    gen t = _n
    tsset t
    forvalues j = 1/4 {
        gen z`j' = rnormal()
    }
    gen w = rnormal()
    gen w2 = rnormal()
    gen u = rnormal()
    replace u = u + .6 * L.u if _n > 1
    gen x = z1 + .5*z2 + .3*z3 + .2*z4 + u + rnormal()
    gen x2 = z2 - z1 + .4*z4 + rnormal()
    gen y = 1 + x + .5*w - .3*w2 + 2*u + rnormal()
    gen cl = ceil(_n / 6)
    gen wt = runiform() + .5
    gen fw = 1 + mod(_n, 3)
end

* s2iv_panel: 40 panels of 15 periods (xtset id t)
capture program drop s2iv_panel
program define s2iv_panel
    clear
    set obs 600
    gen id = ceil(_n / 15)
    bysort id: gen t = _n
    xtset id t
    forvalues j = 1/4 {
        gen z`j' = rnormal() + (mod(id, 3) == 0) * .3
    }
    gen w = rnormal()
    gen u = rnormal()
    by id: replace u = u + .5 * u[_n-1] if _n > 1
    gen x = z1 + .5*z2 + .3*z3 + .2*z4 + u + rnormal() + id / 40
    gen x2 = z2 - z1 + .4*z4 + rnormal()
    gen y = 1 + x + .5*w + 2*u + rnormal() + id / 20
    gen cl = mod(id, 9)
    gen cl2 = mod(t, 5)
    gen wt = runiform() + .5
end

* s2iv_cmp: estimate with ref() and ct(); compare e(b), e(V), the scalars
* (ref=ct name pairs allowed, as cstatdf=cstat_df), the matrices and the
* macros at mreldif 1e-7
capture program drop s2iv_cmp
program define s2iv_cmp
    syntax, ref(string) ct(string) testname(string) [scalars(string) mats(string) ///
        macros(string) nob tol(real 1e-7)]
    ereturn clear
    capture quietly `ref'
    if _rc {
        test_fail "`testname'" "reference failed rc=`=_rc'"
        exit
    }
    tempname rb rV
    matrix `rb' = e(b)
    matrix `rV' = e(V)
    local ref_names : colfullnames e(b)
    foreach s of local scalars {
        gettoken rs cs : s, parse("=")
        local cs : subinstr local cs "=" ""
        if "`cs'" == "" local cs `rs'
        local r_`rs' = e(`rs')
    }
    local i = 0
    foreach m of local mats {
        local ++i
        tempname rm`i'
        capture matrix `rm`i'' = e(`m')
        if _rc {
            test_fail "`testname'" "reference has no e(`m')"
            exit
        }
    }
    foreach m of local macros {
        local rmac_`m' `"`e(`m')'"'
    }
    ereturn clear
    capture quietly `ct'
    if _rc {
        test_fail "`testname'" "civreghdfe failed rc=`=_rc'"
        exit
    }
    local diffs
    if "`b'" == "" {
        local ct_names : colfullnames e(b)
        if "`ct_names'" != "`ref_names'" local diffs "`diffs' names: `ct_names' vs `ref_names'"
        else {
            if mreldif(e(b), `rb') > `tol' local diffs "`diffs' e(b):`=mreldif(e(b), `rb')'"
            if mreldif(e(V), `rV') > `tol' local diffs "`diffs' e(V):`=mreldif(e(V), `rV')'"
        }
    }
    foreach s of local scalars {
        gettoken rs cs : s, parse("=")
        local cs : subinstr local cs "=" ""
        if "`cs'" == "" local cs `rs'
        if missing(`r_`rs'') & missing(e(`cs')) continue
        if missing(`r_`rs'') | missing(e(`cs')) | reldif(`r_`rs'', e(`cs')) > `tol' {
            local diffs "`diffs' e(`cs'):`r_`rs'' vs `=e(`cs')'"
        }
    }
    local i = 0
    foreach m of local mats {
        local ++i
        tempname cm
        capture matrix `cm' = e(`m')
        if _rc local diffs "`diffs' no e(`m')"
        else if rowsof(`cm') != rowsof(`rm`i'') | colsof(`cm') != colsof(`rm`i'') {
            local diffs "`diffs' e(`m') dims"
        }
        else {
            if mreldif(`cm', `rm`i'') > `tol' local diffs "`diffs' e(`m'):`=mreldif(`cm', `rm`i'')'"
            if "`: rownames `cm''" != "`: rownames `rm`i'''" | "`: colfullnames `cm''" != "`: colfullnames `rm`i'''" {
                local diffs "`diffs' e(`m') names"
            }
        }
    }
    foreach m of local macros {
        if `"`e(`m')'"' != `"`rmac_`m''"' local diffs `"`diffs' e(`m'): "`e(`m')'" vs "`rmac_`m''""'
    }
    if "`diffs'" == "" test_pass "`testname'"
    else test_fail "`testname'" "`=trim("`diffs'")'"
end

* s2iv_log <file> : <command> runs the command with its output written to
* a text log; c_local __logrc its return code
capture program drop s2iv_log
program define s2iv_log
    gettoken file 0 : 0
    gettoken colon 0 : 0
    quietly log using `"`file'"', text replace name(s2iv)
    capture noisily `0'
    local rc = _rc
    quietly log close s2iv
    c_local __logrc = `rc'
end

* s2iv_block <file> <from> <to>: r(text) the lines of the log from the
* first one containing <from> up to, not including, the next one
* containing <to> ("" = to the end); r(n) their count. A marker that starts
* with ^ must start the line. The command echo lines (". ") are skipped.
capture program drop s2iv_block
program define s2iv_block, rclass
    args file from to
    mata: s2iv_block(st_local("file"), st_local("from"), st_local("to"))
    return local text `"`macval(__text)'"'
    return scalar n = `__n'
end
mata:
void s2iv_block(string scalar file, string scalar from, string scalar to)
{
    real scalar fh, on, n, f1, t1
    string scalar line, text
    fh = fopen(file, "r")
    on = 0; n = 0; text = ""
    f1 = (substr(from, 1, 1) == "^")
    t1 = (substr(to, 1, 1) == "^")
    if (f1) from = substr(from, 2, .)
    if (t1) to = substr(to, 2, .)
    while ((line = fget(fh)) != J(0, 0, "")) {
        if (substr(line, 1, 2) == ". " | substr(line, 1, 2) == "> ") continue
        if (!on & (f1 ? strpos(line, from) == 1 : strpos(line, from) > 0)) on = 1
        else if (on & to != "" & (t1 ? strpos(line, to) == 1 : strpos(line, to) > 0)) break
        if (on) {
            text = text + strrtrim(line) + "|"
            n++
        }
    }
    fclose(fh)
    st_local("__text", text)
    st_local("__n", strofreal(n))
}
void s2iv_firstdiff(string scalar a, string scalar b)
{
    string rowvector x, y
    real scalar i
    x = ustrsplit(a, "[|]"); y = ustrsplit(b, "[|]")
    for (i = 1; i <= min((cols(x), cols(y))); i++) {
        if (x[i] != y[i]) {
            printf("{txt}REF: %s\n{txt}CT:  %s\n", x[i], y[i])
            return
        }
    }
    printf("{txt}lengths %g vs %g\n", cols(x), cols(y))
}
end

* s2iv_disp: the displayed blocks (from|to pairs) of ref() and ct() are
* the same lines
capture program drop s2iv_disp
program define s2iv_disp
    syntax, ref(string) ct(string) testname(string) blocks(string asis)
    tempfile f1 f2
    s2iv_log `"`f1'"' : `ref'
    local rc1 = `__logrc'
    s2iv_log `"`f2'"' : `ct'
    local rc2 = `__logrc'
    if `rc1' | `rc2' {
        test_fail "`testname'" "rc ref=`rc1' ct=`rc2'"
        exit
    }
    local diffs
    foreach blk of local blocks {
        gettoken from to : blk, parse("|")
        local to : subinstr local to "|" ""
        s2iv_block `"`f1'"' `"`from'"' `"`to'"'
        local t1 `"`r(text)'"'
        local n1 = r(n)
        s2iv_block `"`f2'"' `"`from'"' `"`to'"'
        if `n1' == 0 | `"`r(text)'"' != `"`t1'"' {
            local diffs `"`diffs' [`from']"'
            mata: s2iv_firstdiff(st_local("t1"), st_global("r(text)"))
        }
    }
    if `"`diffs'"' == "" test_pass "`testname'"
    else test_fail "`testname'" `"blocks differ:`diffs'"'
end

* s3_pre <vars> : <command> drops the variables (saved effects of an
* earlier run) and runs the command
capture program drop s3_pre
program define s3_pre
    gettoken vars 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    foreach v of local vars {
        capture drop `v'
    }
    `0'
end
* s3_rest <name> : <command> runs the command and restores the stored
* estimates <name>
capture program drop s3_rest
program define s3_rest
    gettoken eq 0 : 0
    gettoken colon 0 : 0
    quietly `0'
    quietly estimates restore `eq'
end

local ivstats "N rss F j jdf idstat iddf cdf widstat"
local ffstats "arf arfp archi2 archi2p ardf ardf_r sstat sstatp sstatdf"
local ffblock `""Summary results for first-stage regressions|Number of excluded instruments""'

capture noisily {

* ---------------------------------------------------------------------
* IV-25/IV-35: e(first), Anderson-Rubin, Stock-Wright and the ffirst
* display, against ivreg2 and ivreghdfe
* ---------------------------------------------------------------------
s2iv_ts
foreach v in "" "robust" "cluster(cl)" "bw(3)" "robust bw(4) kernel(parzen)" "liml" "gmm2s robust" ///
    "cue robust" "partial(w)" "noconstant robust" "dofminus(2) robust" "sdofminus(3)" "robust small" {
    foreach m in "(x = z1 z2 z3)" "(x x2 = z1 z2 z3 z4)" {
        s2iv_cmp, ref(ivreg2 y w w2 `m', ffirst `v') ct(civreghdfe y w w2 `m', ffirst `v') ///
            scalars(`ivstats' `ffstats') mats(first) testname("IV-25 ffirst `m' `v'")
    }
    s2iv_disp, ref(ivreg2 y w w2 (x x2 = z1 z2 z3 z4), ffirst `v') ///
        ct(civreghdfe y w w2 (x x2 = z1 z2 z3 z4), ffirst `v') blocks(`ffblock') ///
        testname("IV-35 ffirst display `v'")
}
foreach wt in "[aw=wt]" "[fw=fw]" "[pw=wt]" {
    s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4) `wt', ffirst) ct(civreghdfe y w (x x2 = z1 z2 z3 z4) `wt', ffirst) ///
        scalars(`ivstats' `ffstats') mats(first) testname("IV-25 ffirst weights `wt'")
}
s2iv_panel
* (ivreghdfe leaves the first-stage statistics missing with two-way
* clustering: they are compared with ivreg2 below)
foreach v in "" "robust" "cluster(id)" "dkraay(2)" "kiefer" "bw(2) robust" ///
    "liml cluster(id)" "partial(w) robust" {
    s2iv_cmp, ref(ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) ffirst `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) ffirst `v') ///
        scalars(`ivstats' `ffstats') mats(first) testname("IV-25 ffirst absorb `v'")
    s2iv_disp, ref(ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) ffirst `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) ffirst `v') blocks(`ffblock') ///
        testname("IV-35 ffirst display absorb `v'")
}
foreach v in "cluster(cl cl2)" "cluster(id t)" "dkraay(2)" "kiefer" {
    s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), ffirst `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), ffirst `v') ///
        scalars(`ivstats' `ffstats') mats(first) testname("IV-25 ffirst panel `v'")
}
foreach wt in "[aw=wt]" "[pw=wt]" {
    s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3) `wt', absorb(id) ffirst cluster(id)) ///
        ct(civreghdfe y w (x = z1 z2 z3) `wt', absorb(id) ffirst cluster(id)) ///
        scalars(`ivstats' `ffstats') mats(first) testname("IV-25 ffirst absorb weights `wt'")
}

* ---------------------------------------------------------------------
* IV-35: first displays each first stage with ivreg2's F and SW lines
* ---------------------------------------------------------------------
s2iv_ts
foreach v in "" "robust" "cluster(cl)" "bw(2)" {
    s2iv_disp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), first `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), first `v') ///
        blocks(`"First-stage regressions|Summary results"' `ffblock') ///
        testname("IV-35 first display `v'")
}
s2iv_panel
s2iv_disp, ref(ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) first cluster(id)) ///
    ct(civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) first cluster(id)) ///
    blocks(`"First-stage regressions|Summary results"' `ffblock') ///
    testname("IV-35 first display absorb cluster")

* ---------------------------------------------------------------------
* IV-05: sfirst, savesfirst, savesfprefix()
* ---------------------------------------------------------------------
s2iv_ts
foreach v in "" "robust" "cluster(cl)" "bw(3) robust" "partial(w)" "dofminus(1)" {
    s2iv_disp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), sfirst `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), sfirst `v') ///
        blocks(`"System of first-stage|Summary results"' `ffblock') ///
        testname("IV-05 sfirst display `v'")
    * the stored system: e(b), e(V), e(S), dof
    * (ivreg2 and civreghdfe store it under the same name)
    capture estimates drop _ivreg2_sfirst_y
    s2iv_cmp, ref(s3_rest _ivreg2_sfirst_y : ivreg2 y w (x x2 = z1 z2 z3 z4), savesfirst `v') ///
        ct(s3_rest _ivreg2_sfirst_y : civreghdfe y w (x x2 = z1 z2 z3 z4), savesfirst `v') ///
        scalars(N df_r k_eq) mats(S) macros(model depvar) testname("IV-05 savesfirst `v'")
}
s2iv_panel
foreach v in "" "cluster(id)" "cluster(cl cl2)" "dkraay(2)" {
    capture estimates drop _ivreg2_sfirst_y
    s2iv_cmp, ref(s3_rest _ivreg2_sfirst_y : ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) savesfirst `v') ///
        ct(s3_rest _ivreg2_sfirst_y : civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) savesfirst `v') ///
        scalars(N df_r k_eq) mats(S) macros(model depvar) testname("IV-05 savesfirst absorb `v'")
}
quietly civreghdfe y w (x = z1 z2 z3), absorb(id) savesfprefix(my_)
if "`e(sfirsteq)'" == "my_sfirst_y" & "`e(model)'" == "iv" {
    quietly estimates restore my_sfirst_y
    if "`e(model)'" == "sfirst" & colsof(e(b)) == 8 test_pass "IV-05 savesfprefix() implies savesfirst"
    else test_fail "IV-05 savesfprefix() implies savesfirst" "`e(model)' `=colsof(e(b))'"
}
else test_fail "IV-05 savesfprefix() implies savesfirst" "`e(sfirsteq)'"
capture estimates drop my_sfirst_y
quietly civreghdfe y w (x = z1 z2 z3), absorb(id) sfirst
capture estimates describe _ivreg2_sfirst_y
if _rc & "`e(sfirsteq)'" == "" test_pass "IV-05 sfirst alone leaves no stored system"
else test_fail "IV-05 sfirst alone leaves no stored system" "rc=`=_rc'"

* ---------------------------------------------------------------------
* IV-07: Tukey-Hamming, Daniell and Tent kernels (with every VCE)
* ---------------------------------------------------------------------
s2iv_ts
foreach k in thamm dan ten tukey-hamming danielle tent {
    foreach v in "" "robust" "cluster(t)" "gmm2s" "liml robust" {
        s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), bw(4) kernel(`k') `v') ///
            ct(civreghdfe y w (x = z1 z2 z3), bw(4) kernel(`k') `v') ///
            scalars(`ivstats' bw) macros(kernel) testname("IV-07 kernel(`k') `v'")
    }
}
s2iv_panel
foreach k in thamm dan ten {
    foreach v in "robust" "cluster(id t)" "dkraay(3)" {
        s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id) bw(3) kernel(`k') `v') ///
            ct(civreghdfe y w (x = z1 z2 z3), absorb(id) bw(3) kernel(`k') `v') ///
            scalars(`ivstats' bw) macros(kernel) testname("IV-07 absorb kernel(`k') `v'")
    }
    s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3) [aw=wt], absorb(id) bw(2.5) kernel(`k') robust) ///
        ct(civreghdfe y w (x = z1 z2 z3) [aw=wt], absorb(id) bw(2.5) kernel(`k') robust) ///
        scalars(`ivstats' bw) macros(kernel) testname("IV-07 absorb aw kernel(`k')")
}
foreach k in daniell tukey hamming {
    capture ivreg2 y w (x = z1 z2 z3), bw(3) kernel(`k')
    local rc_ref = _rc
    capture civreghdfe y w (x = z1 z2 z3), bw(3) kernel(`k')
    if _rc == `rc_ref' | ("`k'" == "tukey" & _rc == 0) test_pass "IV-07 kernel(`k') accepted as ivreg2"
    else test_fail "IV-07 kernel(`k') accepted as ivreg2" "rc=`=_rc' ref=`rc_ref'"
}

* ---------------------------------------------------------------------
* IV-06: bw(auto), Newey-West (1994)
* ---------------------------------------------------------------------
s2iv_ts
foreach v in "robust" "" "kernel(par) robust" "kernel(qs) robust" "gmm2s" "kernel(qs) gmm2s robust" ///
    "cue robust" "partial(w) robust" "noconstant robust" "cluster(t)" "kernel(thamm) robust" {
    s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), bw(auto) `v') ct(civreghdfe y w (x = z1 z2 z3), bw(auto) `v') ///
        scalars(`ivstats' bw) macros(kernel bwchoice) testname("IV-06 bw(auto) `v'")
}
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3) [aw=wt] if t > 4, bw(auto) robust) ///
    ct(civreghdfe y w (x = z1 z2 z3) [aw=wt] if t > 4, bw(auto) robust) ///
    scalars(`ivstats' bw) macros(bwchoice) testname("IV-06 bw(auto) aw if")
s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(cl) bw(auto) robust) ///
    ct(civreghdfe y w (x = z1 z2 z3), absorb(cl) bw(auto) robust) ///
    scalars(`ivstats' bw) macros(bwchoice) testname("IV-06 bw(auto) absorb")
s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), bw(auto) robust ffirst) ///
    ct(civreghdfe y w (x x2 = z1 z2 z3 z4), bw(auto) robust ffirst) ///
    scalars(`ivstats' `ffstats' bw) mats(first) testname("IV-06 bw(auto) ffirst")
s2iv_disp, ref(ivreg2 y w (x = z1 z2 z3), bw(auto) robust first) ///
    ct(civreghdfe y w (x = z1 z2 z3), bw(auto) robust first) ///
    blocks(`"  kernel=|time variable"' `"Statistics robust|Number"' `"First-stage regressions|Summary results"') ///
    testname("IV-06 bw(auto) display")
s2iv_disp, ref(ivreg2 y w (x = z1 z2 z3), bw(auto) kernel(dan) robust) ///
    ct(civreghdfe y w (x = z1 z2 z3), bw(auto) kernel(dan) robust) ///
    blocks(`"Error: kernel|Number"') testname("IV-06 bw(auto) incompatible kernel display")
foreach b in "0" "-1" "abc" "AUTO" {
    capture ivreg2 y w (x = z1 z2 z3), bw(`b') robust
    local rc_ref = _rc
    capture civreghdfe y w (x = z1 z2 z3), bw(`b') robust
    if _rc == `rc_ref' & _rc != 0 test_pass "IV-06 bw(`b') error"
    else test_fail "IV-06 bw(`b') error" "rc=`=_rc' ref=`rc_ref'"
}
s2iv_panel
capture ivreg2 y w (x = z1 z2 z3), bw(auto) robust
local rc_ref = _rc
capture civreghdfe y w (x = z1 z2 z3), bw(auto) robust
if _rc == `rc_ref' & _rc == 198 test_pass "IV-06 bw(auto) refused for panel data"
else test_fail "IV-06 bw(auto) refused for panel data" "rc=`=_rc' ref=`rc_ref'"

* ---------------------------------------------------------------------
* IV-04: center
* ---------------------------------------------------------------------
s2iv_ts
foreach v in "robust" "gmm2s robust" "cue robust" "liml robust" "fuller(1) robust" "kclass(.5) robust" ///
    "cluster(cl)" "gmm2s cluster(cl)" "robust bw(3)" "gmm2s robust bw(3) kernel(qs)" "gmm2s bw(3)" ///
    "gmm2s" "cluster(t) bw(2)" "gmm2s robust partial(w)" {
    s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), center `v') ct(civreghdfe y w (x x2 = z1 z2 z3 z4), center `v') ///
        scalars(`ivstats' center) testname("IV-04 center `v'")
}
foreach wt in "[aw=wt]" "[pw=wt]" "[fw=fw]" {
    foreach v in "gmm2s robust" "gmm2s cluster(cl)" {
        s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3) `wt', center `v' ffirst) ///
            ct(civreghdfe y w (x = z1 z2 z3) `wt', center `v' ffirst) ///
            scalars(`ivstats' `ffstats') mats(first) testname("IV-04 center `wt' `v'")
    }
}
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), center gmm2s robust orthog(z3)) ///
    ct(civreghdfe y w (x = z1 z2 z3), center gmm2s robust orthog(z3)) ///
    scalars(`ivstats' cstat cstatdf=cstat_df) testname("IV-04 center orthog()")
s2iv_panel
foreach v in "gmm2s robust" "gmm2s cluster(id)" "cluster(cl cl2)" "cue cluster(id)" {
    s2iv_cmp, ref(ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) center `v' ffirst) ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) center `v' ffirst) ///
        scalars(`ivstats' `ffstats' center) mats(first) testname("IV-04 center absorb `v'")
}

* ---------------------------------------------------------------------
* IV-08: smatrix(), wmatrix() (the gmmw model)
* ---------------------------------------------------------------------
s2iv_ts
quietly ivreg2 y w (x = z1 z2 z3), robust
matrix S1 = e(S)
matrix W1 = invsym(S1)
matrix Wi = I(5)
matrix rownames Wi = z1 z2 z3 w _cons
matrix colnames Wi = z1 z2 z3 w _cons
* a superset in another order: rows and columns are picked by name
matrix Sbig = (S1, J(5, 1, 0) \ J(1, 6, 0))
matrix rownames Sbig = z1 z2 z3 w _cons extra
matrix colnames Sbig = z1 z2 z3 w _cons extra
foreach v in "" "robust" "gmm2s" "gmm2s robust" "cluster(cl)" "gmm2s cluster(cl)" "robust bw(3)" "liml robust" "cue robust" {
    foreach mo in "wmatrix(Wi)" "smatrix(S1)" "wmatrix(W1)" "smatrix(Sbig)" "smatrix(S1) wmatrix(Wi)" {
        s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), `mo' `v') ct(civreghdfe y w (x = z1 z2 z3), `mo' `v') ///
            scalars(`ivstats') testname("IV-08 `mo' `v'")
    }
}
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3) [aw=wt], wmatrix(Wi) robust) ct(civreghdfe y w (x = z1 z2 z3) [aw=wt], wmatrix(Wi) robust) ///
    scalars(`ivstats') testname("IV-08 wmatrix() aw")
quietly civreghdfe y w (x = z1 z2 z3), wmatrix(Wi) robust
if "`e(model)'" == "gmmw" test_pass "IV-08 e(model) gmmw"
else test_fail "IV-08 e(model) gmmw" "`e(model)'"
s2iv_disp, ref(ivreg2 y w (x = z1 z2 z3), wmatrix(Wi) robust) ct(civreghdfe y w (x = z1 z2 z3), wmatrix(Wi) robust) ///
    blocks(`"GMM estimation|Number of obs"') testname("IV-08 gmmw title display")
s2iv_disp, ref(ivreg2 y w (x = z1 z2 z3), wmatrix(Wi) gmm2s) ct(civreghdfe y w (x = z1 z2 z3), wmatrix(Wi) gmm2s) ///
    blocks(`"2-Step GMM|Number of obs"') testname("IV-08 gmm2s wmatrix title display")
matrix Sasym = S1
matrix Sasym[1, 2] = 5
foreach mo in "smatrix(S1[1..4,1..4])" "smatrix(nosuch)" "smatrix(Sasym)" "wmatrix(Sasym)" {
    capture ivreg2 y w (x = z1 z2 z3), `mo'
    local rc_ref = _rc
    capture civreghdfe y w (x = z1 z2 z3), `mo'
    if _rc == `rc_ref' & _rc != 0 test_pass "IV-08 `mo' error"
    else test_fail "IV-08 `mo' error" "rc=`=_rc' ref=`rc_ref'"
}
s2iv_panel
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id) cluster(id)
matrix SA = e(S)
matrix WA = I(4)
matrix rownames WA = z1 z2 z3 w
matrix colnames WA = z1 z2 z3 w
foreach v in "cluster(id)" "gmm2s cluster(id)" "" {
    foreach mo in "smatrix(SA)" "wmatrix(WA)" {
        s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id) `mo' `v') ct(civreghdfe y w (x = z1 z2 z3), absorb(id) `mo' `v') ///
            scalars(`ivstats') testname("IV-08 absorb `mo' `v'")
    }
}

* ---------------------------------------------------------------------
* IV-10: orthog() on included (and excluded) instruments
* ---------------------------------------------------------------------
s2iv_ts
foreach o in "w" "w2" "w z3" "w w2" "z4" {
    foreach v in "" "robust" "gmm2s robust" "cluster(cl)" "liml" "bw(2) robust" {
        s2iv_cmp, ref(ivreg2 y w w2 (x = z1 z2 z3 z4), orthog(`o') `v') ///
            ct(civreghdfe y w w2 (x = z1 z2 z3 z4), orthog(`o') `v') ///
            scalars(`ivstats' cstat cstatdf=cstat_df) testname("IV-10 orthog(`o') `v'")
    }
}
s2iv_cmp, ref(ivreg2 y w w2 (x = z1 z2 z3 z4) [aw=wt], orthog(w2) robust) ///
    ct(civreghdfe y w w2 (x = z1 z2 z3 z4) [aw=wt], orthog(w2) robust) ///
    scalars(`ivstats' cstat cstatdf=cstat_df) testname("IV-10 orthog() aw")
s2iv_panel
foreach v in "" "cluster(id)" "gmm2s cluster(id)" {
    s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3 z4), absorb(id) orthog(w) `v') ///
        ct(civreghdfe y w (x = z1 z2 z3 z4), absorb(id) orthog(w) `v') ///
        scalars(`ivstats' cstat cstatdf=cstat_df) testname("IV-10 absorb orthog(w) `v'")
}
capture civreghdfe y w (x = z1 z2 z3 z4), absorb(id) orthog(nosuchvar)
if _rc test_pass "IV-10 orthog() of a non-instrument refused"
else test_fail "IV-10 orthog() of a non-instrument refused" "accepted"

* ---------------------------------------------------------------------
* IV-01: replay
* ---------------------------------------------------------------------
s2iv_ts
foreach v in "robust" "cluster(cl) liml" "bw(2)" {
    s2iv_disp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), `v') ///
        blocks(`"^Estimates|^Instrumented:"') testname("IV-01 main display matches ivreg2 `v'")
    quietly ivreg2 y w (x x2 = z1 z2 z3 z4), `v' savefirst saverf savesfirst
    tempfile r1 r2 r3
    s2iv_log `"`r1'"' : ivreg2, first ffirst rf sfirst level(90) eform(exp)
    quietly civreghdfe y w (x x2 = z1 z2 z3 z4), `v' savefirst saverf savesfirst
    s2iv_log `"`r2'"' : civreghdfe, first ffirst rf sfirst level(90) eform(exp)
    local fails
    foreach blk in `"^Reduced-form regression|^First-stage regressions"' `"^First-stage regressions|^System of"' ///
        `"^System of first-stage|^Summary results"' `ffblock' `"^Estimates|^Instrumented:"' {
        gettoken from to : blk, parse("|")
        local to : subinstr local to "|" ""
        s2iv_block `"`r1'"' `"`from'"' `"`to'"'
        local t1 `"`r(text)'"'
        local n1 = r(n)
        s2iv_block `"`r2'"' `"`from'"' `"`to'"'
        if `n1' == 0 | `"`r(text)'"' != `"`t1'"' {
            local fails `"`fails' [`from']"'
            mata: s2iv_firstdiff(st_local("t1"), st_global("r(text)"))
        }
    }
    if `"`fails'"' == "" test_pass "IV-01 replay first ffirst rf sfirst level() eform() `v'"
    else test_fail "IV-01 replay first ffirst rf sfirst level() eform() `v'" `"`fails'"'
    * replay of the results equals the estimation display
    quietly civreghdfe y w (x x2 = z1 z2 z3 z4), `v'
    s2iv_log `"`r1'"' : civreghdfe y w (x x2 = z1 z2 z3 z4), `v'
    s2iv_log `"`r2'"' : civreghdfe
    s2iv_block `"`r1'"' "estimation" ""
    local t1 `"`r(text)'"'
    s2iv_block `"`r2'"' "estimation" ""
    if `"`r(text)'"' == `"`t1'"' test_pass "IV-01 bare replay equals the estimation display `v'"
    else test_fail "IV-01 bare replay equals the estimation display `v'" "differs"
    s2iv_log `"`r3'"' : estimates replay
    s2iv_block `"`r3'"' "estimation" ""
    if `"`r(text)'"' == `"`t1'"' test_pass "IV-01 estimates replay `v'"
    else test_fail "IV-01 estimates replay `v'" "differs"
    foreach eq in _ivreg2_x _ivreg2_x2 _ivreg2_y _ivreg2_sfirst_y {
        capture estimates drop `eq'
    }
}
* replay options noheader, nofooter, plus
quietly civreghdfe y w (x = z1 z2 z3), robust
tempfile r1 r2
s2iv_log `"`r1'"' : civreghdfe, noheader nofooter plus
s2iv_block `"`r1'"' "IV (2SLS) estimation" ""
local n_head = r(n)
s2iv_block `"`r1'"' "Hansen J" ""
local n_foot = r(n)
s2iv_block `"`r1'"' "Coefficient" ""
if `n_head' == 0 & `n_foot' == 0 & r(n) > 0 test_pass "IV-01 replay noheader nofooter"
else test_fail "IV-01 replay noheader nofooter" "head=`n_head' foot=`n_foot'"
* drop* options drop the stored equations
quietly civreghdfe y w (x = z1 z2 z3), robust savefirst saverf savesfirst
quietly civreghdfe, dropfirst droprf dropsfirst
local ok = ("`e(firsteqs)'`e(rfeq)'`e(sfirsteq)'" == "")
foreach eq in _ivreg2_x _ivreg2_y _ivreg2_sfirst_y {
    capture estimates describe `eq'
    if !_rc local ok = 0
}
if `ok' test_pass "IV-01 replay drop* options"
else test_fail "IV-01 replay drop* options" "equations kept"
* by: with replay, and replay after another command
sort cl
capture noisily by cl: civreghdfe
local rc = _rc
sort t
if `rc' == 601 test_pass "IV-01 by: replay refused"
else test_fail "IV-01 by: replay refused" "rc=`rc'"
quietly regress y x
capture civreghdfe
if _rc == 301 test_pass "IV-01 replay after another command"
else test_fail "IV-01 replay after another command" "rc=`=_rc'"
* replay of first without stored first stages
quietly ivreg2 y w (x = z1 z2 z3), ffirst
tempfile r1 r2
s2iv_log `"`r1'"' : ivreg2, first
quietly civreghdfe y w (x = z1 z2 z3), ffirst
s2iv_log `"`r2'"' : civreghdfe, first
s2iv_block `"`r1'"' "First-stage regressions" "estimation"
local t1 `"`r(text)'"'
s2iv_block `"`r2'"' "First-stage regressions" "estimation"
if `"`r(text)'"' == `"`t1'"' & `"`t1'"' != "" test_pass "IV-01 replay first without stored equations"
else test_fail "IV-01 replay first without stored equations" "differs"
s2iv_panel
quietly ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) cluster(id) savefirst saverf savesfirst
tempfile r1 r2
s2iv_log `"`r1'"' : ivreghdfe, first ffirst sfirst
quietly civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) cluster(id) savefirst saverf savesfirst
s2iv_log `"`r2'"' : civreghdfe, first ffirst sfirst
local fails
foreach blk in `"^First-stage regressions|^System of"' `"^System of first-stage|^Summary results"' `ffblock' {
    gettoken from to : blk, parse("|")
    local to : subinstr local to "|" ""
    s2iv_block `"`r1'"' `"`from'"' `"`to'"'
    local t1 `"`r(text)'"'
    local n1 = r(n)
    s2iv_block `"`r2'"' `"`from'"' `"`to'"'
    if `n1' == 0 | `"`r(text)'"' != `"`t1'"' {
        local fails `"`fails' [`from']"'
        mata: s2iv_firstdiff(st_local("t1"), st_global("r(text)"))
    }
}
if `"`fails'"' == "" test_pass "IV-01 replay absorb first ffirst sfirst"
else test_fail "IV-01 replay absorb first ffirst sfirst" `"`fails'"'

/* IV-02: heterogeneous slopes in absorb(), against ivreghdfe (reghdfe's
   slopes): e(b), e(V), the absorbed degrees of freedom, the statistics,
   e(extended_absvars), the saved effects and the absorbed DoF table */
clear
set obs 3000
gen id = 1 + mod(_n, 40)
gen t = 1 + mod(floor(_n / 7), 25)
gen grp = 1 + mod(_n, 3)
gen cl = 1 + mod(_n, 60)
gen x1 = rnormal() + id / 20
gen x2 = rnormal()
gen x3 = rnormal()
forvalues j = 1/3 {
    gen z`j' = rnormal()
}
gen d = .5*z1 + .3*z2 - .2*z3 + .4*x1 + rnormal()
gen y = 1 + 2*d - x1 + .5*x2 - .4*x3 + id/10 + (id/40)*x1 + t/30*x2 + rnormal() * (1 + abs(x1))
gen w = .5 + runiform()
gen fw = 1 + floor(runiform() * 3)
gen xc = cond(id <= 5, 1, x1)
local iv02_sc N df_a df_r rmse rss widstat idstat j cdf df_m F N_hdfe N_hdfe_extended
local iv02_specs `""absorb(id##c.x1)" "absorb(id#c.x1 t) cluster(id)" "absorb(id##c.(x1 x2) t)""'
local iv02_specs `"`iv02_specs' "absorb(id##c.x1 t) robust" "absorb(id##c.x1 t) cluster(cl t)""'
local iv02_specs `"`iv02_specs' "absorb(id##c.x1 t) liml" "absorb(id##c.x1 t) gmm2s robust""'
local iv02_specs `"`iv02_specs' "absorb(id grp#t##c.x2)" "absorb(t#c.x1 id)" "absorb(id##c.x1 t grp)""'
local iv02_specs `"`iv02_specs' "absorb(id##c.x1 t) dofadjustments(pairwise)" "absorb(id##c.xc t)""'
local iv02_specs `"`iv02_specs' "absorb(id#c.xc t)" "absorb(id#c.x2)""'
foreach s of local iv02_specs {
    s2iv_cmp, ref(ivreghdfe y (d = z1 z2 z3) x3, `s') ct(civreghdfe y (d = z1 z2 z3) x3, `s') ///
        scalars(`iv02_sc') macros(extended_absvars) testname("IV-02 slopes `s'")
}
foreach wt in "aw=w" "fw=fw" "pw=w" {
    s2iv_cmp, ref(ivreghdfe y (d = z1 z2 z3) x3 [`wt'], absorb(id##c.x1 t)) ///
        ct(civreghdfe y (d = z1 z2 z3) x3 [`wt'], absorb(id##c.x1 t)) ///
        scalars(`iv02_sc') macros(extended_absvars) testname("IV-02 slopes [`wt']")
}
s2iv_disp, ref(ivreghdfe y (d = z1 z2 z3) x3, absorb(id##c.(x1 x2) t grp)) ///
    ct(civreghdfe y (d = z1 z2 z3) x3, absorb(id##c.(x1 x2) t grp)) ///
    blocks(`"^Absorbed degrees of freedom|^? = number"') testname("IV-02 slopes absorbed DoF table")
* Saved effects: newvar= and savefe, intercepts and slopes
quietly ivreghdfe y (d = z1 z2 z3) x3, absorb(fe=id##c.x2 t, savefe)
rename (fe feSlope1 __hdfe2__) (rfe rfeSlope1 rfe2)
quietly civreghdfe y (d = z1 z2 z3) x3, absorb(fe=id##c.x2 t, savefe)
local fails
if "`e(savefe)'" != "fe feSlope1 __hdfe2__" local fails "`fails' e(savefe)=`e(savefe)'"
foreach p in "fe rfe" "feSlope1 rfeSlope1" "__hdfe2__ rfe2" {
    gettoken a b : p
    tempvar dd
    quietly gen double `dd' = abs(`a' - `b')
    summarize `dd', meanonly
    if !(r(max) < 1e-6) local fails "`fails' `a':`r(max)'"
}
if "`fails'" == "" test_pass "IV-02 slopes saved effects"
else test_fail "IV-02 slopes saved effects" "`fails'"
drop fe feSlope1 __hdfe2__ rfe rfeSlope1 rfe2
capture civreghdfe y (d = z1 z2 z3) x2, absorb(id##c.x2 t) vce(bootstrap, reps(5) seed(1))
if _rc == 198 test_pass "IV-02 slopes refuse vce(bootstrap)"
else test_fail "IV-02 slopes refuse vce(bootstrap)" "rc=`=_rc'"

/* S3/S4 gaps IV-11 to IV-39: ivreg2's and ivreghdfe's stored results
   (e(S), e(W), e(vce), e(model), titles, 0/1 varlists, test names, the
   footnote results), options (iweights, psd0/psda, tvar()/ivar(),
   bvclean, partial(_cons), reghdfe's options, abbreviations), the
   stored equations, their replay, and out-of-sample residuals */
local s3sc N j jdf sargandf cd rkf condzz condxx dofminus center inexog_ct lambda kclass ///
    arubin arubin_lin arubinp arubin_linp arubindf idstat cdf widstat r2_a F
local s3mac model vce title hacsubtitleB hacsubtitleV hacsubtitleV2 partialsmall depvar0 depvar1 ///
    instd instd0 instd1 inexog inexog0 inexog1 exexog exexog0 exexog1 insts insts0 insts1 partial partial1
local s3foot `""^Underidentification test|^Instrumented:""'
s2iv_ts
foreach v in "" "robust" "liml" "fuller(1)" "kclass(.5)" "liml kclass(.8)" "gmm2s" "gmm2s robust" ///
    "cue" "cue robust" "cluster(cl)" "small dofminus(2)" "robust dofminus(2)" "partial(w)" ///
    "noconstant robust" "bw(3)" "bw(3) robust" "bw(2) kernel(qs) robust gmm2s" "noid" "center robust" {
    local mats S W
    if strpos("`v'", "liml") | strpos("`v'", "fuller") | strpos("`v'", "kclass") local mats S
    if !strpos("`v'", "noid") local mats `mats' ccev cdev
    s2iv_cmp, ref(ivreg2 y w w2 (x x2 = z1 z2 z3 z4), `v') ct(civreghdfe y w w2 (x x2 = z1 z2 z3 z4), `v') ///
        scalars(`s3sc') mats(`mats') macros(`s3mac') testname("IV-23/24/28/29/30 stored results `v'")
}
* without endogenous regressors (OLS) and with one
s2iv_cmp, ref(ivreg2 y w w2 x) ct(civreghdfe y w w2 x) scalars(`s3sc') mats(S W) macros(`s3mac') ///
    testname("IV-23/30 stored results OLS")
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), robust) ct(civreghdfe y w (x = z1 z2 z3), robust) ///
    scalars(`s3sc') mats(S W ccev cdev) macros(`s3mac') testname("IV-23 stored results one endogenous")
foreach wt in "[aw=wt]" "[fw=fw]" "[pw=wt]" {
    s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4) `wt') ct(civreghdfe y w (x x2 = z1 z2 z3 z4) `wt') ///
        scalars(`s3sc') mats(S W ccev cdev) macros(`s3mac') testname("IV-23 stored results `wt'")
}
* wmatrix() and smatrix(): e(W) and e(S) as given
quietly ivreg2 y w (x = z1 z2 z3)
matrix s3W = 2 * e(W)
foreach v in "wmatrix(s3W)" "gmm2s wmatrix(s3W)" "smatrix(s3W)" {
    s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), `v') ct(civreghdfe y w (x = z1 z2 z3), `v') ///
        scalars(`s3sc') mats(S W) macros(`s3mac') testname("IV-23 e(S) e(W) `v'")
}
* kernels: e(vce) with the kernel and bandwidth, kiefer, dkraay
foreach v in "bw(3) kernel(par) robust" "bw(2)" "cluster(t) bw(2)" {
    s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), `v') ct(civreghdfe y w (x = z1 z2 z3), `v') ///
        scalars(`s3sc') mats(S W) macros(`s3mac') testname("IV-24 e(vce) `v'")
}
s2iv_panel
foreach v in "kiefer" "cluster(id t) bw(2)" "cluster(cl cl2)" "dkraay(2)" "dkraay(2) robust" {
    s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), `v') ct(civreghdfe y w (x = z1 z2 z3), `v') ///
        scalars(`s3sc') mats(S W) macros(`s3mac') testname("IV-24 e(vce) `v'")
}
* absorb(): against ivreghdfe, with reghdfe's footnote results (IV-30e)
local s3fsc `s3sc' df_a df_a_initial df_a_redundant df_a_nested drop_singletons N_hdfe N_hdfe_extended
local s3fmac `s3mac' absvars extended_absvars title2 dofmethod
foreach v in "absorb(id)" "absorb(id t)" "absorb(i.id i.t) robust" "absorb(id t) cluster(cl)" ///
    "absorb(id t) cluster(id)" "absorb(id) gmm2s robust" "absorb(id) liml" "absorb(id t) keepsingletons" ///
    "absorb(id t) dofadjustments(none)" "absorb(id t) dofadjustments(pairwise)" "absorb(id#cl2 t)" ///
    "absorb(fe1=id t)" "absorb(id##c.wt)" {
    local mats S W ccev cdev dof_table
    if strpos("`v'", "liml") local mats S ccev cdev dof_table
    s2iv_cmp, ref(s3_pre fe1 : ivreghdfe y w (x x2 = z1 z2 z3 z4), `v') ///
        ct(s3_pre fe1 : civreghdfe y w (x x2 = z1 z2 z3 z4), `v') ///
        scalars(`s3fsc') mats(`mats') macros(`s3fmac') testname("IV-30e absorb stored results `v'")
    capture drop fe1
}

* IV-27: ivreg2's names for the C, endogeneity and redundancy tests, and
* IV-37: the first footer (restricted J with orthog(), redundancy and
* endogeneity blocks, LIML's Anderson-Rubin statistic) as ivreg2 shows it
local s3tsc cstat cstatp cstatdf estat estatp estatdf redstat redp reddf j jdf jp
s2iv_ts
foreach v in "orthog(z3)" "orthog(z3 z4 w) robust" "endog(x)" "endog(x x2) robust" "redundant(z4)" ///
    "orthog(z4) endog(x2) redundant(z3)" "liml" "fuller(1) orthog(z3)" "gmm2s robust orthog(w2)" "noid redundant(z4)" {
    s2iv_cmp, ref(ivreg2 y w w2 (x x2 = z1 z2 z3 z4), `v') ct(civreghdfe y w w2 (x x2 = z1 z2 z3 z4), `v') ///
        scalars(`s3tsc') macros(clist elist redlist) testname("IV-27 test names `v'")
    local blk `"`s3foot'"'
    if strpos("`v'", "noid") local blk `""^Sargan statistic|^Instrumented:""'
    s2iv_disp, ref(ivreg2 y w w2 (x x2 = z1 z2 z3 z4), `v') ct(civreghdfe y w w2 (x x2 = z1 z2 z3 z4), `v') ///
        blocks(`blk') testname("IV-37 footer `v'")
}
s2iv_disp, ref(ivreg2 y w w2 x) ct(civreghdfe y w w2 x) blocks(`"^OLS estimation|^Included instruments"') ///
    testname("IV-37 OLS footer")
s2iv_panel
foreach v in "orthog(z3)" "endog(x) cluster(cl)" "redundant(z4) robust" "liml orthog(z4)" {
    s2iv_cmp, ref(ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) `v') ///
        scalars(`s3tsc') macros(clist elist redlist) testname("IV-27 test names absorb `v'")
    s2iv_disp, ref(ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) `v') ///
        ct(civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) `v') blocks(`s3foot') testname("IV-37 footer absorb `v'")
}

* IV-12: psd0 and psda (kernels whose S need not be positive semidefinite)
s2iv_ts
foreach k in tent thann tru {
    foreach p in psd0 psda {
        local vs `""robust" """'
        if "`p'" == "psda" local vs `"`vs' "gmm2s robust""'
        foreach v of local vs {
            s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4), bw(4) kernel(`k') `v' `p') ///
                ct(civreghdfe y w (x x2 = z1 z2 z3 z4), bw(4) kernel(`k') `v' `p') ///
                scalars(`s3sc' rankS) mats(S) macros(vce) testname("IV-12 `p' kernel(`k') `v'")
        }
    }
}
capture civreghdfe y w (x = z1 z2 z3), bw(3) robust psd0 psda
local rc1 = _rc
capture ivreg2 y w (x = z1 z2 z3), bw(3) robust psd0 psda
if `rc1' == _rc & `rc1' == 198 test_pass "IV-12 psd0 with psda refused"
else test_fail "IV-12 psd0 with psda refused" "rc `rc1' vs `=_rc'"

* IV-11: iweights (homoskedastic VCE, no absorb()); with a non-integer sum
* ivreg2 truncates N in later statistics, so only e(b), e(V), N and the
* residual dof are compared there
gen int s3iw = 1 + mod(_n, 4)
gen double s3iwr = .5 + runiform()
foreach v in "" "small" "liml" "partial(w)" {
    s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4) [iw=s3iw], `v') ct(civreghdfe y w (x x2 = z1 z2 z3 z4) [iw=s3iw], `v') ///
        scalars(`s3sc') mats(S) macros(`s3mac') testname("IV-11 iweights `v'")
}
s2iv_cmp, ref(ivreg2 y w (x x2 = z1 z2 z3 z4) [iw=s3iwr], small) ct(civreghdfe y w (x x2 = z1 z2 z3 z4) [iw=s3iwr], small) ///
    scalars(N df_r rss j) testname("IV-11 iweights non-integer")
foreach v in "robust" "gmm2s" "bw(2)" {
    capture civreghdfe y w (x = z1 z2 z3) [iw=s3iw], `v'
    local rc1 = _rc
    capture ivreg2 y w (x = z1 z2 z3) [iw=s3iw], `v'
    if `rc1' == _rc & `rc1' == 101 test_pass "IV-11 iweights refused with `v'"
    else test_fail "IV-11 iweights refused with `v'" "rc `rc1' vs `=_rc'"
}
gen id = ceil(_n / 10)
capture civreghdfe y w (x = z1 z2 z3) [iw=s3iw], absorb(id)
local rc1 = _rc
capture ivreghdfe y w (x = z1 z2 z3) [iw=s3iw], absorb(id)
if `rc1' == _rc & `rc1' == 101 test_pass "IV-11 iweights refused with absorb()"
else test_fail "IV-11 iweights refused with absorb()" "rc `rc1' vs `=_rc'"
drop id s3iw s3iwr

* IV-14: tvar() and ivar() name the tsset variables, as ivreg2 requires
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), bw(2) tvar(t)) ct(civreghdfe y w (x = z1 z2 z3), bw(2) tvar(t)) ///
    scalars(`s3sc') macros(vce) testname("IV-14 tvar()")
gen t2 = t
foreach o in "tvar(t2)" "ivar(t2)" {
    capture civreghdfe y w (x = z1 z2 z3), bw(2) `o'
    local rc1 = _rc
    capture ivreg2 y w (x = z1 z2 z3), bw(2) `o'
    if `rc1' == _rc & `rc1' == 5 test_pass "IV-14 `o' refused"
    else test_fail "IV-14 `o' refused" "rc `rc1' vs `=_rc'"
}
drop t2
s2iv_panel
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), bw(2) robust tvar(t) ivar(id)) ///
    ct(civreghdfe y w (x = z1 z2 z3), bw(2) robust tvar(t) ivar(id)) scalars(`s3sc') macros(vce) ///
    testname("IV-14 tvar() ivar() panel")

* IV-15: bvclean (e(b) without base levels), fvall, fvsep and useqr
gen a = mod(_n, 4)
gen b = mod(_n, 3)
foreach v in "bvclean" "bvclean robust" "fvall" "fvsep" "useqr" "fvall robust" {
    s2iv_cmp, ref(ivreg2 y i.a w (x = z1 z2 i.b), `v') ct(civreghdfe y i.a w (x = z1 z2 i.b), `v') ///
        scalars(`s3sc') macros(`s3mac') testname("IV-15 `v'")
}
s2iv_cmp, ref(ivreg2 y i.a##c.w (x = z1 z2 i.a#c.z3), fvsep) ct(civreghdfe y i.a##c.w (x = z1 z2 i.a#c.z3), fvsep) ///
    scalars(`s3sc') testname("IV-15 fvsep interactions")
s2iv_cmp, ref(ivreghdfe y i.a w (x = z1 z2 i.b), absorb(id) bvclean) ///
    ct(civreghdfe y i.a w (x = z1 z2 i.b), absorb(id) bvclean) scalars(`s3sc') macros(`s3mac') ///
    testname("IV-15 bvclean absorb")
s2iv_cmp, ref(ivreg2 y i.a w (x = z1 z2 i.b), partial(w)) ct(civreghdfe y i.a w (x = z1 z2 i.b), partial(w)) ///
    scalars(`s3sc') macros(`s3mac') testname("IV-15 partial() implies bvclean")

* IV-16: partial(_cons) and fwl(_cons ...)
s2iv_ts
foreach v in "partial(_cons)" "partial(_cons) robust" "fwl(_cons w)" "partial(w _cons) liml" "partial(_all)" {
    s2iv_cmp, ref(ivreg2 y w w2 (x = z1 z2 z3), `v') ct(civreghdfe y w w2 (x = z1 z2 z3), `v') ///
        scalars(`s3sc' partial_ct partialcons cons) mats(S) macros(`s3mac') testname("IV-16 `v'")
}
gen id = ceil(_n / 10)
foreach v in "absorb(id) partial(_cons)" "noconstant partial(_cons)" "partial(_cons _cons)" {
    capture civreghdfe y w (x = z1 z2 z3), `v'
    local rc1 = _rc
    capture ivreghdfe y w (x = z1 z2 z3), `v'
    if `rc1' == _rc & `rc1' == 198 test_pass "IV-16 `v' refused"
    else test_fail "IV-16 `v' refused" "rc `rc1' vs `=_rc'"
}
drop id

* IV-18: reghdfe's options, which ivreghdfe passes on: groupvar(), the
* solver options (same estimates), nosample, timeit, nowarn, verbose(#)
s2iv_panel
gen fe2 = mod(_n, 13)
foreach o in "technique(lsmr)" "acceleration(sd) transform(cim)" "transform(kac) acceleration(none)" ///
    "preconditioner(diag) technique(lsqr)" "prune poolsize(5)" "nosample" "nowarn" "verbose(-1)" "timeit" {
    s2iv_cmp, ref(ivreghdfe y w (x = z1 z2 z3), absorb(id fe2) `o') ct(civreghdfe y w (x = z1 z2 z3), absorb(id fe2) `o') ///
        scalars(N df_a j idstat) testname("IV-18 `o'")
}
foreach o in "technique(xyz)" "transform(foo)" "acceleration(bar)" {
    capture civreghdfe y w (x = z1 z2 z3), absorb(id fe2) `o'
    local rc1 = _rc
    capture ivreghdfe y w (x = z1 z2 z3), absorb(id fe2) `o'
    if `rc1' == 198 & _rc != 0 test_pass "IV-18 invalid `o' refused"
    else test_fail "IV-18 invalid `o' refused" "rc `rc1' vs `=_rc'"
}
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id fe2) groupvar(s3g_ref)
quietly civreghdfe y w (x = z1 z2 z3), absorb(id fe2) groupvar(s3g_ct)
quietly count if s3g_ref != s3g_ct
local nd = r(N)
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id) groupvar(s3g1_ref)
quietly civreghdfe y w (x = z1 z2 z3), absorb(id) groupvar(s3g1_ct)
capture confirm variable s3g1_ref
local r1 = _rc
capture confirm variable s3g1_ct
if `nd' == 0 & `r1' == _rc test_pass "IV-18 groupvar()"
else test_fail "IV-18 groupvar()" "`nd' values differ; one effect rc `r1' vs `=_rc'"
capture drop s3g_ref s3g_ct s3g1_ref s3g1_ct
drop fe2

* IV-21: ivreg2's abbreviations sm, nohe, nofo
s2iv_cmp, ref(ivreg2 y w (x = z1 z2 z3), sm) ct(civreghdfe y w (x = z1 z2 z3), sm) ///
    scalars(N df_r F) testname("IV-21 sm")
s2iv_disp, ref(ivreg2 y w (x = z1 z2 z3), nohe nofo) ct(civreghdfe y w (x = z1 z2 z3), nohe nofo) ///
    blocks(`"Coefficient|s3_no_such_line"') testname("IV-21 nohe nofo")

* IV-22: residuals2 saves _reghdfe_resid, replaced on a rerun
capture drop _reghdfe_resid
quietly ivreghdfe y w (x = z1 z2 z3), absorb(id) residuals2
rename _reghdfe_resid s3r_ref
quietly civreghdfe y w (x = z1 z2 z3), absorb(id) residuals2
capture quietly civreghdfe y w (x = z1 z2 z3), absorb(id) residuals2
local rc2 = _rc
quietly count if reldif(s3r_ref, _reghdfe_resid) > 1e-7 & !missing(s3r_ref)
if r(N) == 0 & `rc2' == 0 test_pass "IV-22 residuals2"
else test_fail "IV-22 residuals2" "`r(N)' differ; rerun rc `rc2'"
drop s3r_ref _reghdfe_resid

* IV-30c/30d: stored equations under ivreg2's names, with e(S); IV-39:
* their replay and e(sample)
s2iv_ts
foreach v in "" "robust" "cluster(cl)" "bw(2) robust" "partial(w)" {
    capture estimates drop _ivreg2_*
    quietly ivreg2 y w (x x2 = z1 z2 z3 z4), savefirst saverf `v'
    local rnames `e(firsteqs)' `e(rfeq)'
    quietly civreghdfe y w (x x2 = z1 z2 z3 z4), savefirst saverf `v'
    if "`e(firsteqs)' `e(rfeq)'" == "`rnames'" test_pass "IV-30c stored names `v'"
    else test_fail "IV-30c stored names `v'" "`e(firsteqs)' `e(rfeq)' vs `rnames'"
    foreach eq of local rnames {
        s2iv_cmp, ref(s3_rest `eq' : ivreg2 y w (x x2 = z1 z2 z3 z4), savefirst saverf `v') ///
            ct(s3_rest `eq' : civreghdfe y w (x x2 = z1 z2 z3 z4), savefirst saverf `v') nob ///
            scalars(N df_r df_m k_eq rmse) mats(S) macros(model depvar vce hacsubtitleV) ///
            testname("IV-30d stored equation `eq' `v'")
        tempfile r1 r2
        s3_rest `eq' : ivreg2 y w (x x2 = z1 z2 z3 z4), savefirst saverf `v'
        quietly count if e(sample)
        local nref = r(N)
        s2iv_log `"`r1'"' : ivreg2
        s3_rest `eq' : civreghdfe y w (x x2 = z1 z2 z3 z4), savefirst saverf `v'
        quietly count if e(sample)
        local nct = r(N)
        s2iv_log `"`r2'"' : civreghdfe
        s2iv_block `"`r1'"' "regression" ""
        local t1 `"`r(text)'"'
        s2iv_block `"`r2'"' "regression" ""
        if `"`r(text)'"' == `"`t1'"' & `"`t1'"' != "" & `nref' == `nct' test_pass "IV-39 replay and e(sample) `eq' `v'"
        else {
            test_fail "IV-39 replay and e(sample) `eq' `v'" "e(sample) `nct' vs `nref'"
            mata: s2iv_firstdiff(st_local("t1"), st_global("r(text)"))
        }
    }
}
capture estimates drop _ivreg2_*
quietly ivreg2 y w (x x2 = z1 z2 z3 z4), savesfirst robust
quietly estimates restore _ivreg2_sfirst_y
tempfile r1 r2
s2iv_log `"`r1'"' : ivreg2
quietly civreghdfe y w (x x2 = z1 z2 z3 z4), savesfirst robust
quietly estimates restore _ivreg2_sfirst_y
s2iv_log `"`r2'"' : civreghdfe
s2iv_block `"`r1'"' "System of" ""
local t1 `"`r(text)'"'
s2iv_block `"`r2'"' "System of" ""
if `"`r(text)'"' == `"`t1'"' & `"`t1'"' != "" test_pass "IV-39 replay of the stored system"
else {
    test_fail "IV-39 replay of the stored system" "differs"
    mata: s2iv_firstdiff(st_local("t1"), st_global("r(text)"))
}
s2iv_panel
foreach v in "" "cluster(id)" {
    capture estimates drop _ivreg2_*
    quietly ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) savefirst saverf `v'
    local rnames `e(firsteqs)' `e(rfeq)'
    quietly civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) savefirst saverf `v'
    if "`e(firsteqs)' `e(rfeq)'" == "`rnames'" test_pass "IV-30c stored names absorb `v'"
    else test_fail "IV-30c stored names absorb `v'" "`e(firsteqs)' `e(rfeq)' vs `rnames'"
    foreach eq of local rnames {
        s2iv_cmp, ref(s3_rest `eq' : ivreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) savefirst saverf `v') ///
            ct(s3_rest `eq' : civreghdfe y w (x x2 = z1 z2 z3 z4), absorb(id) savefirst saverf `v') nob ///
            scalars(N df_r k_eq rmse) mats(S) macros(model depvar) testname("IV-30d stored equation absorb `eq' `v'")
    }
}
capture estimates drop _ivreg2_*

* IV-38: residuals predicted out of the estimation sample, as ivreg2_p
s2iv_ts
foreach v in "" "robust" "partial(w)" "partial(_cons)" "noconstant partial(w)" {
    quietly ivreg2 y w w2 (x = z1 z2 z3) if t <= 100, `v'
    quietly predict double s3r1, residuals
    quietly civreghdfe y w w2 (x = z1 z2 z3) if t <= 100, `v'
    quietly predict double s3r2, residuals
    quietly count if reldif(s3r1, s3r2) > 1e-7 | missing(s3r1) != missing(s3r2)
    local nd = r(N)
    quietly count if !missing(s3r2) & t > 100
    if `nd' == 0 & r(N) == 50 test_pass "IV-38 out-of-sample residuals `v'"
    else test_fail "IV-38 out-of-sample residuals `v'" "`nd' differ; `r(N)' out of sample"
    drop s3r1 s3r2
}

}
if _rc test_fail "s2_iv aborted" "r(`=_rc')"

print_summary s2_iv
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_iv"
    global CTOOLS_CASE_COMPLETE "s2_iv"
    exit
}

di as error "Unknown s2_iv case: `__case'"
exit 198
