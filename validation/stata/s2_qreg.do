* s2_qreg: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_qreg"

if "`__case'" == "s2_qreg" {
do validation/support/validate_setup.do
set seed 20261104

* QREG-01: cqreg accepts fweights, iweights and pweights as qreg does.
* Each comparison runs qreg and then cqreg and checks e(b), e(V), every
* scalar qreg posts and the weight and VCE macros at 1e-7 relative
* difference. The data have continuous responses and weights, so the
* weighted optimum is unique and the coefficients themselves are compared.

* Store the reference results of the last qreg
capture program drop _s2q_save
program define _s2q_save
    matrix __s2q_b = e(b)
    matrix __s2q_V = e(V)
    foreach s in sum_adev sum_rdev q_v N sum_w sparsity bwidth kbwidth f_r ///
        df_r df_m rank r2_p q N_clust {
        scalar __s2q_`s' = e(`s')
    }
    foreach m in wtype wexp vce vcetype denmethod bwmethod kernel clustvar {
        global S2Q_m_`m' `"`e(`m')'"'
    }
end

* Compare the last cqreg with the stored qreg results; r(msg) lists the
* results that differ. nob skips e(b) (a non-unique optimum).
capture program drop _s2q_check
program define _s2q_check, rclass
    syntax [, tol(real 1e-7) nob]
    local msg ""
    if "`b'" == "" {
        if mreldif(__s2q_b, e(b)) > `tol' local msg "`msg' e(b)"
        if "`: colfullnames __s2q_b'" != "`: colfullnames e(b)'" local msg "`msg' names"
        if mreldif(__s2q_V, e(V)) > `tol' local msg "`msg' e(V)"
    }
    foreach s in sum_adev sum_rdev q_v N sum_w sparsity bwidth kbwidth f_r ///
        df_r df_m rank r2_p q N_clust {
        if "`b'" != "" & inlist("`s'", "sparsity", "f_r", "kbwidth") continue
        local a = __s2q_`s'
        local c = e(`s')
        if missing(`a') != missing(`c') local msg "`msg' e(`s')"
        else if !missing(`a') & reldif(`a', `c') > `tol' local msg "`msg' e(`s')"
    }
    foreach m in wtype wexp vce vcetype denmethod bwmethod kernel clustvar {
        if `"${S2Q_m_`m'}"' != `"`e(`m')'"' local msg "`msg' e(`m')"
    }
    return local msg "`msg'"
end

* Run "qreg ..." and "cqreg ...", compare, record the test
capture program drop _s2q_run
program define _s2q_run
    args name qcmd ccmd extra
    capture noisily {
        ereturn clear
        quietly `qcmd'
        _s2q_save
        ereturn clear
        quietly `ccmd'
        _s2q_check, `extra'
    }
    if _rc test_fail "`name'" "r(`=_rc')"
    else if "`r(msg)'" != "" test_fail "`name'" "differs:`r(msg)'"
    else test_pass "`name'"
end

* Error code of a command (0 when it runs)
capture program drop _s2q_rc
program define _s2q_rc, rclass
    args cmd
    capture quietly `cmd'
    return scalar rc = _rc
end

sysuse auto, clear
gen fw = 1 + floor(runiform()*4)
gen double iw = 0.2 + runiform()*3
gen grp = mod(_n, 7)

* -------------------------------------------------------------------------
* Weight type x VCE x density method x bandwidth
* -------------------------------------------------------------------------
foreach w in fw iw pw {
    local wv = cond("`w'" == "fw", "fw", "iw")
    foreach v in iid robust {
        if "`w'" == "pw" & "`v'" == "iid" continue
        foreach d in fitted residual kernel {
            if "`d'" == "residual" & "`v'" != "iid" continue
            foreach bw in hsheather bofinger chamberlain {
                _s2q_run "QREG-01 [`w'] vce(`v') denmethod(`d') bwmethod(`bw')" ///
                    "qreg price weight length [`w'=`wv'], vce(`v', `d' `bw')" ///
                    "cqreg price weight length [`w'=`wv'], vce(`v') denmethod(`d') bwmethod(`bw')"
            }
        }
        foreach k in epanechnikov epan2 biweight cosine gaussian parzen rectangle triangle {
            _s2q_run "QREG-01 [`w'] vce(`v') denmethod(kernel(`k'))" ///
                "qreg price weight length [`w'=`wv'], vce(`v', kernel(`k'))" ///
                "cqreg price weight length [`w'=`wv'], vce(`v') denmethod(kernel(`k'))"
        }
    }
    foreach q in .25 .75 {
        _s2q_run "QREG-01 [`w'] quantile(`q')" ///
            "qreg price weight length [`w'=`wv'], q(`q')" ///
            "cqreg price weight length [`w'=`wv'], q(`q')"
    }
    * the constant-only model: e(q_v), e(sum_rdev) and the weighted quantile
    _s2q_run "QREG-01 [`w'] constant only" ///
        "qreg price [`w'=`wv'], q(.4)" "cqreg price [`w'=`wv'], q(.4)"
    * qreg's vce(excluster) is cqreg's vce(cluster)
    foreach d in fitted kernel {
        _s2q_run "QREG-01 [`w'] vce(cluster) denmethod(`d')" ///
            "qreg price weight length [`w'=`wv'], vce(excluster grp, `d')" ///
            "cqreg price weight length [`w'=`wv'], vce(cluster grp) denmethod(`d')"
    }
    _s2q_run "QREG-01 [`w'] vce(cluster) bwmethod(chamberlain)" ///
        "qreg price weight length [`w'=`wv'], vce(excluster grp, chamberlain)" ///
        "cqreg price weight length [`w'=`wv'], vce(cluster grp) bwmethod(chamberlain)"
}

* pweights imply vce(robust)
_s2q_run "QREG-01 [pw] default vce(robust)" ///
    "qreg price weight length [pw=iw]" "cqreg price weight length [pw=iw]"
* weight expressions and abbreviations
_s2q_run "QREG-01 [fweight=] expression" ///
    "qreg price weight length [fweight=fw*2]" "cqreg price weight length [fweight=fw*2]"
_s2q_run "QREG-01 [iweight=] expression" ///
    "qreg price weight length [iweight=iw/10]" "cqreg price weight length [iweight=iw/10]"
* zero weights are excluded from the sample (e(N) counts the rest)
gen zw = fw
replace zw = 0 in 1/5
_s2q_run "QREG-01 [fw] zero weights" ///
    "qreg price weight length [fw=zw]" "cqreg price weight length [fw=zw]"
_s2q_run "QREG-01 [iw] zero weights" ///
    "qreg price weight length [iw=zw], vce(robust)" "cqreg price weight length [iw=zw], vce(robust)"
* weights on a half-integer grid: cumulative weights reach the percentile
* targets of the residual and kernel methods exactly
gen double hw = 0.5 * (1 + floor(runiform()*4))
foreach d in residual kernel {
    _s2q_run "QREG-01 [iw] half-integer weights denmethod(`d')" ///
        "qreg price weight length [iw=hw], vce(iid, `d')" ///
        "cqreg price weight length [iw=hw], denmethod(`d')"
}
_s2q_run "QREG-01 [iw] half-integer weights constant only residual" ///
    "qreg price [iw=hw], vce(iid, residual)" "cqreg price [iw=hw], denmethod(residual)"
* iweights over eight orders of magnitude
gen double wide = 10^(4*(2*runiform() - 1))
_s2q_run "QREG-01 [iw] weights from 1e-4 to 1e4" ///
    "qreg price weight length [iw=wide], vce(robust)" "cqreg price weight length [iw=wide], vce(robust)"

* A non-unique weighted optimum: indicator regressors with integer
* fweights. qreg and cqreg can report different optimal vertices, so the
* minimized loss is compared, and cqreg's e(sum_adev) must be the weighted
* loss of its own e(b).
capture noisily {
    quietly qreg price i.rep78 [fw=fw]
    _s2q_save
    ereturn clear
    quietly cqreg price i.rep78 [fw=fw]
    _s2q_check, nob
    if "`r(msg)'" != "" {
        di as error "differs:`r(msg)'"
        error 9
    }
    tempvar xb loss
    quietly predict double `xb' if e(sample), xb
    quietly gen double `loss' = fw * cond(price - `xb' >= 0, .5, -.5) * (price - `xb') if e(sample)
    quietly summarize `loss', meanonly
    if reldif(r(sum), e(sum_adev)) > 1e-7 error 9
}
if _rc test_fail "QREG-01 [fw] non-unique optimum: e(sum_adev)" "r(`=_rc')"
else test_pass "QREG-01 [fw] non-unique optimum: e(sum_adev)"

* -------------------------------------------------------------------------
* Rejected combinations, with qreg's error codes
* -------------------------------------------------------------------------
foreach spec in ///
    "[pw=iw], vce(iid)|[pw=iw], vce(iid)|404" ///
    "[pw=iw], vce(iid, residual)|[pw=iw], vce(iid) denmethod(residual)|404" ///
    "[pw=iw], vce(, residual)|[pw=iw], denmethod(residual)|404" ///
    "[pw=iw], vce(robust, residual)|[pw=iw], vce(robust) denmethod(residual)|184" ///
    "[iw=iw], vce(robust, residual)|[iw=iw], vce(robust) denmethod(residual)|184" ///
    "[aw=iw]|[aw=iw]|101" ///
    "[fw=iw]|[fw=iw]|401" {
    gettoken qs rest : spec, parse("|")
    gettoken bar rest : rest, parse("|")
    gettoken cs rest : rest, parse("|")
    local code = substr("`rest'", 2, .)
    _s2q_rc "qreg price weight `qs'"
    local rq = r(rc)
    _s2q_rc "cqreg price weight `cs'"
    local rc = r(rc)
    if `rq' == `code' & `rc' == `code' test_pass "QREG-01 rejects `cs' (r(`code'))"
    else test_fail "QREG-01 rejects `cs' (r(`code'))" "qreg r(`rq'), cqreg r(`rc')"
}
* pweights with a negative weight: r(402) in both
gen double npw = iw - 1
_s2q_rc "qreg price weight [pw=npw]"
local rq = r(rc)
_s2q_rc "cqreg price weight [pw=npw]"
if `rq' == 402 & r(rc) == 402 test_pass "QREG-01 negative pweights r(402)"
else test_fail "QREG-01 negative pweights r(402)" "qreg r(`rq'), cqreg r(`=r(rc)')"
* Contract: a negative iweight makes the weighted check loss non-convex.
* qreg's simplex still reports a vertex; cqreg refuses (r(402)).
_s2q_rc "cqreg price weight [iw=npw]"
if r(rc) == 402 test_pass "QREG-01 contract: negative iweights r(402)"
else test_fail "QREG-01 contract: negative iweights r(402)" "r(`=r(rc)')"
* a failed weighted call leaves no plugin scalars behind
local left : all scalars "__cqreg_*"
if "`left'" == "" test_pass "QREG-01 no scalars left after a rejected weight"
else test_fail "QREG-01 no scalars left after a rejected weight" "`left'"

* -------------------------------------------------------------------------
* predict after a weighted fit
* -------------------------------------------------------------------------
capture noisily {
    quietly qreg price weight length [pw=iw]
    quietly predict double q_xb, xb
    quietly predict double q_se, stdp
    quietly predict double q_r, residuals
    quietly cqreg price weight length [pw=iw]
    quietly predict double c_xb, xb
    quietly predict double c_se, stdp
    quietly predict double c_r, residuals
    foreach v in xb se r {
        tempvar d
        quietly gen double `d' = reldif(q_`v', c_`v')
        quietly summarize `d'
        if r(max) > 1e-7 {
            di as error "predict `v' differs: " r(max)
            error 9
        }
    }
}
if _rc test_fail "QREG-01 predict xb, stdp, residuals after [pw]" "r(`=_rc')"
else test_pass "QREG-01 predict xb, stdp, residuals after [pw]"

* -------------------------------------------------------------------------
* by: and mi estimate with weights
* -------------------------------------------------------------------------
foreach w in fw iw pw {
    local wv = cond("`w'" == "fw", "fw", "iw")
    capture noisily {
        sort foreign, stable
        quietly by foreign: qreg price weight [`w'=`wv']
        _s2q_save
        quietly by foreign: cqreg price weight [`w'=`wv']
        _s2q_check
        if "`r(msg)'" != "" {
            di as error "differs:`r(msg)'"
            error 9
        }
    }
    if _rc test_fail "QREG-01 by: with [`w']" "r(`=_rc')"
    else test_pass "QREG-01 by: with [`w']"
}

preserve
gen double x = mpg
replace x = . if runiform() < 0.2
quietly mi set mlong
quietly mi register imputed x
quietly mi impute regress x price weight, add(3) rseed(5)
foreach w in fw iw pw {
    local wv = cond("`w'" == "fw", "fw", "iw")
    capture noisily {
        quietly mi estimate: qreg price x weight [`w'=`wv']
        tempname qb qV
        matrix `qb' = e(b_mi)
        matrix `qV' = e(V_mi)
        quietly mi estimate: cqreg price x weight [`w'=`wv']
        if mreldif(`qb', e(b_mi)) > 1e-7 | mreldif(`qV', e(V_mi)) > 1e-7 error 9
    }
    if _rc test_fail "QREG-01 mi estimate with [`w']" "r(`=_rc')"
    else test_pass "QREG-01 mi estimate with [`w']"
}
restore

* -------------------------------------------------------------------------
* A larger panel: nlswork
* -------------------------------------------------------------------------
ctools_fixture nlswork, clear
keep if !missing(ln_wage, age, tenure, ttl_exp, grade)
gen fw = 1 + floor(runiform()*5)
gen double iw = exp(rnormal())
foreach w in fw iw pw {
    local wv = cond("`w'" == "fw", "fw", "iw")
    foreach q in .5 .9 {
        if "`w'" != "pw" {
            _s2q_run "QREG-01 nlswork [`w'] q(`q') vce(iid)" ///
                "qreg ln_wage age tenure ttl_exp grade [`w'=`wv'], q(`q')" ///
                "cqreg ln_wage age tenure ttl_exp grade [`w'=`wv'], q(`q')"
        }
        _s2q_run "QREG-01 nlswork [`w'] q(`q') vce(robust, kernel)" ///
            "qreg ln_wage age tenure ttl_exp grade [`w'=`wv'], q(`q') vce(robust, kernel)" ///
            "cqreg ln_wage age tenure ttl_exp grade [`w'=`wv'], q(`q') vce(robust) denmethod(kernel)"
    }
    _s2q_run "QREG-01 nlswork [`w'] vce(cluster idcode)" ///
        "qreg ln_wage age tenure ttl_exp grade [`w'=`wv'], vce(excluster idcode)" ///
        "cqreg ln_wage age tenure ttl_exp grade [`w'=`wv'], vce(cluster idcode)"
}

* -------------------------------------------------------------------------
* QREG-07/08/10/11/13 (docs/development/FUNCTIONALITY_GAPS.md): qreg's
* vce(type, density bandwidth) syntax and abbreviations, wlsiter(),
* iterate() and trace, and the hidden e(basis), e(tolerance) and the $S_E_*
* globals, against qreg
* -------------------------------------------------------------------------
sysuse auto, clear
set seed 20261109
generate grp = mod(_n, 7)
generate fw = 1 + floor(runiform()*4)
generate double iw = 0.2 + runiform()*3
foreach v in "rob" "r" "robust, bo" "robust, hs" "robust, ch" "iid, res" "iid, residuals" ///
    "iid, fit" "robust, kernel(gaussian) bofinger" "robust, ker(gau)" ", residual" ", kernel" ///
    "iid, kernel" "robust, kernel(epan2) chamberlain" "excluster grp, fit" "excl grp, ker(tri) bo" ///
    "iid, residual bofinger" ", ch" "iid, kernel(parzen) chamberlain" {
    _s2q_run "QREG-07/08 vce(`v')" "qreg price weight length, vce(`v')" ///
        "cqreg price weight length, vce(`v')"
}
* cqreg's own denmethod() and bwmethod() take qreg's names and abbreviations
foreach pair in "iid, residual bofinger|denmethod(res) bwmethod(bo)" ///
    "iid, residual|denmethod(residuals)" "robust, kernel(gaussian)|vce(robust) denmethod(ker(gau))" ///
    "robust, chamberlain|vce(r) bwmethod(ch)" "iid, hsheather|bwmethod(HS)" ///
    "robust, kernel bofinger|vce(robust, kernel) bwmethod(bofinger)" {
    gettoken q c : pair, parse("|")
    local c : subinstr local c "|" ""
    _s2q_run "QREG-08 `c'" "qreg price weight length, vce(`q')" "cqreg price weight length, `c'"
}
* The same error codes as qreg for invalid specifications
foreach v in "robust iid" "robust, fit res" "robust, bo hs" "robust, bogus" "robust, kernel(bogus)" ///
    "iid foo" "robust, residual" "cluster" "bogus" {
    _s2q_rc "qreg price weight length, vce(`v')"
    local qrc = r(rc)
    _s2q_rc "cqreg price weight length, vce(`v')"
    if `qrc' & r(rc) == `qrc' test_pass "QREG-07 vce(`v') refused, r(`qrc')"
    else test_fail "QREG-07 vce(`v') refused, r(`qrc')" "r(`=r(rc)')"
}
_s2q_rc "qreg price weight length [pw=iw], vce(, residual)"
local qrc = r(rc)
_s2q_rc "cqreg price weight length [pw=iw], vce(, residual)"
if `qrc' == 404 & r(rc) == 404 test_pass "QREG-07 vce(, residual) with pweights, r(404)"
else test_fail "QREG-07 vce(, residual) with pweights, r(404)" "r(`qrc') r(`=r(rc)')"
_s2q_rc "cqreg price weight length, vce(robust, bo) bwmethod(hs)"
if r(rc) == 198 test_pass "QREG-07 conflicting vce() and bwmethod() bandwidths"
else test_fail "QREG-07 conflicting vce() and bwmethod() bandwidths" "r(`=r(rc)')"

* QREG-10/11: wlsiter(), iterate() and trace are accepted and checked as
* qreg checks them (the interior-point solver's limit is maxiter())
foreach o in "wlsiter(3)" "iterate(50)" "trace" "wlsiter(2) iterate(1000) trace" {
    _s2q_run "QREG-10/11 `o'" "qreg price weight length, `o'" "cqreg price weight length, `o'"
}
foreach o in "wlsiter(0)" "iterate(0)" "iterate(-2)" "iterate(1.5)" "wlsiter(-1)" {
    _s2q_rc "qreg price weight length, `o'"
    local qrc = r(rc)
    _s2q_rc "cqreg price weight length, `o'"
    if `qrc' & r(rc) == `qrc' test_pass "QREG-10/11 `o' refused, r(`qrc')"
    else test_fail "QREG-10/11 `o' refused, r(`qrc')" "r(`=r(rc)')"
}

* QREG-13: e(basis) is the same set of observations as qreg's (unique
* solutions), e(tolerance) is posted, and $S_E_* hold qreg's values
foreach m in "price weight length" "price weight length [fw=fw]" "price weight length, q(.25)" ///
    "price weight length, vce(robust)" {
    capture noisily {
        quietly qreg `m'
        local q_basis `e(basis)'
        local q_basis : list sort q_basis
        local q_globals
        foreach g in S_E_mdf S_E_tdf S_E_rd S_E_nobs S_E_q S_E_rq S_E_rsd S_E_msd S_E_vl S_E_frc {
            local q_`g' `"${`g'}"'
        }
        quietly cqreg `m'
        local c_basis `e(basis)'
        local c_basis : list sort c_basis
        local bad
        if "`q_basis'" != "`c_basis'" | `: word count `c_basis'' != colsof(e(b)) local bad `bad' basis
        if !(e(tolerance) > 0 & e(tolerance) < 1) local bad `bad' tolerance
        foreach g in S_E_mdf S_E_tdf S_E_nobs S_E_q S_E_rq S_E_frc {
            if `"`q_`g''"' != `"${`g'}"' local bad `bad' `g'
        }
        foreach g in S_E_rd S_E_rsd S_E_msd {
            if reldif(`q_`g'', ${`g'}) > 1e-7 local bad `bad' `g'
        }
        if `"`q_S_E_vl'"' != `"$S_E_vl"' local bad `bad' S_E_vl
        if "$S_E_cmd" != "cqreg" local bad `bad' S_E_cmd
    }
    if _rc test_fail "QREG-13 hidden results: `m'" "r(`=_rc')"
    else if "`bad'" != "" test_fail "QREG-13 hidden results: `m'" "`bad' [`q_basis'] [`c_basis']"
    else test_pass "QREG-13 hidden results: `m'"
}

* The header as qreg's Replay prints it (estimation and replay)
capture program drop _s2q_header
program define _s2q_header, rclass
    args file cmd
    capture log close _s2qh
    quietly log using `"`file'"', text replace name(_s2qh)
    noisily `cmd'
    quietly log close _s2qh
    tempname h
    local text
    file open `h' using `"`file'"', read text
    while 1 {
        file read `h' line
        if r(eof) continue, break
        if strpos(`"`macval(line)'"', "regression") | strpos(`"`macval(line)'"', "sum of deviations") ///
            local text `"`macval(text)'|`macval(line)'"'
    }
    file close `h'
    return local text `"`macval(text)'"'
end
tempfile s2q_h1 s2q_h2
foreach m in "price weight length" "price weight length, q(.25)" "price weight length, q(.1234)" ///
    "price weight length [fw=fw], q(75)" "mpg weight" {
    quietly qreg `m'
    _s2q_header `"`s2q_h1'"' "qreg"
    local t1 `"`r(text)'"'
    quietly cqreg `m'
    _s2q_header `"`s2q_h2'"' "cqreg"
    if `"`t1'"' == `"`r(text)'"' & `"`t1'"' != "" test_pass "QREG header: `m'"
    else test_fail "QREG header: `m'" `"[`r(text)'] vs [`t1']"'
}

capture scalar drop __s2q_*
capture matrix drop __s2q_b __s2q_V
macro drop S2Q_m_*

print_summary s2_qreg
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_qreg"
    global CTOOLS_CASE_COMPLETE "s2_qreg"
    exit
}

di as error "Unknown s2_qreg case: `__case'"
exit 198
