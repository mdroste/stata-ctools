* Regression coverage for audit A01-A29. Run from the repository root.
version 16
clear all
set more off
set seed 9854
adopath ++ package/ado
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
args fixtures
if `"`fixtures'"' == "" local fixtures "validation/fixtures/oct22"

* A01, A11-A14: independent IV reference and positive-DF tiny model.
set obs 100
gen double z=rnormal()
gen double x=z+rnormal()
gen double y=10+.05*x+rnormal()
foreach options in "small" "small robust" "noconstant small" {
    quietly ivregress 2sls y (x=z), `options'
    matrix expected=e(V)
    quietly test x
    scalar expected_F=r(F)
    quietly summarize y
    scalar expected_tss=(r(N)-1)*r(Var)
    if strpos("`options'","noconstant") scalar expected_tss=expected_tss+r(N)*r(mean)^2
    quietly civreghdfe y (x=z), `options'
    matrix actual=e(V)
    noi di "IV `options' Vdiff=" mreldif(actual,expected) " F=" e(F) " expected=" expected_F " TSSdiff=" reldif(e(tss),expected_tss)
    assert mreldif(actual,expected)<1e-8
    assert reldif(e(F),expected_F)<1e-8
    assert reldif(e(tss),expected_tss)<1e-12
    assert e(df_m)==1
}
quietly civreghdfe y (x=z)
matrix large=e(V)
scalar rms=e(rmse)
quietly civreghdfe y (x=z), small
matrix small=e(V)
assert reldif(small[1,1]/large[1,1],100/98)<1e-10
quietly civreghdfe y (x=z), small dofminus(2) sdofminus(3)
matrix adjusted=e(V)
assert e(df_r)==93
assert reldif(adjusted[1,1]/small[1,1],98/93)<1e-10
quietly civreghdfe y (x=z), robust
matrix large=e(V)
quietly civreghdfe y (x=z), robust dofminus(2)
matrix adjusted=e(V)
assert reldif(adjusted[1,1]/large[1,1],100/98)<1e-10
clear
input double(y x z)
1 1 1
4 3 2
2 2 3
5 4 4
end
quietly civreghdfe y (x=z), small
assert e(N)==4 & e(df_r)==2

* A02/A16: fitted and kernel cluster covariance, omitted/base rank.
clear
set obs 2000
set seed 9854
gen double x=5*runiform()
gen double y=2+.5*x+(1+4*x)*rnormal()
gen cl=ceil(_n/20)
foreach density in fitted kernel(gaussian) {
    quietly qreg y x, vce(excluster cl, `density')
    matrix expected=e(V)
    quietly cqreg y x, vce(cluster cl) denmethod(`density')
    matrix actual=e(V)
    assert mreldif(actual,expected)<1e-7
}
gen double duplicate=x
gen category=mod(_n,3)
foreach rhs in "x duplicate" "x i.category" {
    quietly qreg y `rhs'
    scalar expected_rank=e(rank)
    quietly cqreg y `rhs'
    assert e(rank)==expected_rank
}
replace cl=1
quietly cqreg y x, vce(cluster cl)
matrix actual=e(V)
scalar undefined_v=actual[1,1]
assert missing(undefined_v)

* A03/A04: invalid metadata/indexes and record-boundary SHP truncation.
foreach bad in truncated_sst truncated_styles truncated_book truncated_rels invalid_index {
    clear
    set obs 1
    gen sentinel=99
    capture cimport excel "`fixtures'/`bad'.xlsx", clear
    assert _rc!=0
    assert sentinel==99 & _N==1
}
capture cimport shp "`fixtures'/truncated.shp", clear
assert _rc!=0
assert _N==0
cimport shp "`fixtures'/good.shp", clear
assert _N==2

* A08/A09: XML text nodes and OOXML strings, native interchange both ways.
foreach textcase in cdata comment shared_cdata shared_comment escaped {
    cimport excel "`fixtures'/`textcase'.xlsx", clear
    if inlist("`textcase'","cdata","shared_cdata") assert A=="customer & name"
    if inlist("`textcase'","comment","shared_comment") assert A=="beforeafter"
    if "`textcase'"=="escaped" assert A=="a_x000D_b"
}
clear
set obs 4
gen str80 text=cond(_n==1,"a"+char(1)+"b", cond(_n==2,"a"+char(13)+"b",cond(_n==3,"_x0041_","_x005F_x0041_")))
tempfile exported native original
save `original'
cexport excel using "`exported'.xlsx", replace
import excel using "`exported'.xlsx", clear
rename A text
cf _all using `original'
use `original', clear
export excel using "`native'.xlsx", replace
cimport excel using "`native'.xlsx", clear
rename A text
cf _all using `original'

* The existing workbook's shared-string indices survive appends/updates.
use `original', clear
cexport excel using "`exported'.xlsx", sheet(Sheet1, modify) cell(B1)
import excel using "`exported'.xlsx", clear
assert A==B
rename A text
drop B
cf _all using `original'
erase "`exported'.xlsx"
erase "`native'.xlsx"

* A05: weight() inside a quoted condition must remain data.
clear
set obs 6
gen str30 s=cond(_n<=3,"weight(w)","other")
gen w=42
cbsample if s=="weight(w)"
assert _N==3 & w==42 & s=="weight(w)"
cbsample, weight(w)
quietly summarize w
assert r(sum)==3

* A06/A07: Unicode codepoint removal and metadata for generate/replace.
clear
set obs 2
gen str20 s=cond(_n==1,"12©","12é")
clonevar expected=s
clonevar actual=s
quietly destring expected, ignore("é¢") replace
quietly cdestring actual, ignore("é¢") replace
assert actual==expected
replace s="12é"
label var s "Original label"
char s[note1] "Original characteristic"
quietly destring s, gen(expected_num) ignore("é¢")
quietly cdestring s, gen(actual_num) ignore("é¢")
assert actual_num==expected_num
assert `"`: variable label actual_num'"'==`"`: variable label expected_num'"'
assert `"`: char actual_num[note1]'"'==`"`: char expected_num[note1]'"'
quietly cdestring s, replace ignore("é¢")
assert s==12
assert `"`: variable label s'"'=="Original label"
assert `"`: char s[note1]'"'=="Original characteristic"

* A10: prefix threshold, multi-row windows and excludeself cancellation.
clear
set obs 64
gen double x=_n-1
replace x=1e16 in 1
gen key=_n
crangestat (mean) m=x (sum) s=x, interval(key 0 0)
assert m==x & s==x
crangestat (mean) m2=x (sum) s2=x, interval(key -1 1) excludeself
assert m2==x & s2==2*x if inrange(_n,3,63)

* A15: every outcome-dependent PPML prediction accepts TS expressions.
clear
set obs 100
gen time=_n
tsset time
gen double x=rnormal()
gen double y=exp(.2*x+rnormal())
quietly cpplmhdfe L.y x, noabsorb
predict double mu, mu
foreach option in response scores pearson deviance anscombe working {
    predict double result, `option'
    assert !missing(result) if e(sample)
    if "`option'"=="response" assert reldif(result,L.y-mu)<1e-12 if e(sample)
    drop result
}

* A17-A19: expanded using varlists, empty _n merge and assert before keep.
clear
set obs 3
gen id=_n
gen a1=10*_n
gen a2=20*_n
gen extra=99
tempfile usingfile
save `usingfile'
foreach keepvars in "a*" "a1-a2" {
    clear
    set obs 3
    gen id=_n
    cmerge 1:1 id using `usingfile', keepusing(`keepvars')
    assert a1==10*id & a2==20*id
    capture confirm variable extra
    assert _rc==111
}
clear
gen master=.
cmerge 1:1 _n using `usingfile', keepusing(a1)
assert _N==3
capture confirm variable extra
assert _rc==111
capture confirm variable id
assert _rc==111
clear
set obs 0
gen id=.
capture cmerge 1:1 id using `usingfile', assert(match) keep(match)
assert _rc==9
assert _N==3 & _merge==2

* A20/A21/A27-A29: user frame ownership, singleton points, graph and metadata.
clear
set obs 8
gen x=mod(_n,2)
gen y=_n
gen str1 g=cond(_n==8,"b","a")
gen category=mod(_n,3)
frame create _ctools_graph_frame
frame _ctools_graph_frame: set obs 1
frame _ctools_graph_frame: gen protected=99
cbinscatter y x, discrete by(g) linetype(none) nograph
frame _ctools_graph_frame: assert protected==99
matrix bins=e(bindata)
matrix list bins
assert rowsof(bins)>=3
assert bins[3,4]==8 & bins[3,5]==1
assert "`e(by)'"=="g"
foreach line in lfit qfit {
    capture cbinscatter y x, discrete by(g) linetype(`line') nograph
    assert _rc==2001
}
cbinscatter y x, controls(i.category) nograph
assert "`e(controls)'"=="i.category"
program define twoway
    assert strpos(`"`0'"',"connect(l)")>0
end
cbinscatter y x, linetype(connect)
cbinscatter y x, by(g) linetype(connect)
program drop twoway
foreach line in none connect {
    capture cbinscatter y x, linetype(`line') reportreg nograph
    assert _rc==198
}
cbinscatter y x, reportreg nograph
assert "`e(cmd)'"=="cbinscatter"
frame drop _ctools_graph_frame

* Reporting must also work with residualization, groups, weights and polynomials.
clear
set obs 100
gen double x=rnormal()
gen double c=rnormal()
gen double y=2*x+c+rnormal()
gen group=mod(_n,2)
gen fe=mod(_n,10)
gen w=1+mod(_n,3)
foreach options in "controls(c)" "controls(c) absorb(fe)" ///
    "controls(c) by(group) linetype(qfit)" "method(binsreg) controls(c) absorb(fe)" {
    quietly cbinscatter y x [aw=w], `options' nograph
    matrix expected_bins=e(bindata)
    cbinscatter y x [aw=w], `options' reportreg nograph
    matrix actual_bins=e(bindata)
    assert mreldif(actual_bins,expected_bins)==0
    assert "`e(cmd)'"=="cbinscatter"
}

* A23-A25: preserve collided inputs; reference domains and treatment coding.
clear
set obs 10
gen t=mod(_n,2)
gen double p=_n/11
gen y=_n
gen double _pscore=p
capture cpsmatch t, pscore(_pscore) outcome(y)
assert _rc==111
assert _pscore==p
capture cpsmatch t, pscore(p) neighbor(0)
assert _rc==111
capture cpsmatch t, pscore(p) kerneltype(unknown)
assert _rc==198
cpsmatch t, pscore(p) outcome(y) caliper(-1)
assert missing(r(att))
assert _nn==0 if t==1
cpsmatch t, pscore(p) outcome(y) kernel bwidth(0)
assert missing(r(att))
replace t=2*t
cpsmatch t, pscore(p) outcome(y)
assert missing(r(att))
assert missing(_nn)
assert _weight==0
noi di "OCT22_REGRESSIONS_COMPLETE"
