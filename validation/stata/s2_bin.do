* s2_bin: named cases share this file, each do-call has its own local scope.
* BIN-02/BIN-03 (cbinscatter vs binscatter 7.02) and WINS-01/04..08 (cwinsor vs
* gstats winsor and winsor2).
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_bin"

if "`__case'" == "s2_bin" {
do validation/support/validate_setup.do
set seed 20261104

* ---------------------------------------------------------------------------
* Helpers. A graph is compared on binscatter's e(graphcmd) and cbinscatter's:
* the scatteri coordinates, each series' marker and line colors and symbols,
* each fit line's range and its values at the range's ends and middle, and
* the rd() xline(), with relative tolerance tol (reldif). With svg, the two
* graphs exported as SVG must also be identical.
* ---------------------------------------------------------------------------
mata:
void _s2b_split(string scalar g, string colvector plots, string scalar opts)
{
    real scalar i, n, depth, start
    string scalar ch
    plots = J(0, 1, "")
    opts = ""
    n = strlen(g)
    depth = 0
    start = 0
    for (i = 1; i <= n; i++) {
        ch = substr(g, i, 1)
        if (ch == "(") {
            if (depth == 0) start = i
            depth++
        }
        else if (ch == ")") {
            depth--
            if (depth == 0) plots = plots \ substr(g, start + 1, i - start - 1)
        }
        else if (ch == "," & depth == 0) {
            opts = substr(g, i + 1, .)
            break
        }
    }
}

real scalar _s2b_comma(string scalar s)
{
    real scalar i, depth
    string scalar ch
    depth = 0
    for (i = 1; i <= strlen(s); i++) {
        ch = substr(s, i, 1)
        if (ch == "(") depth++
        else if (ch == ")") depth--
        else if (ch == "," & depth == 0) return(i)
    }
    return(0)
}

string scalar _s2b_opt(string scalar s0, string scalar name)
{
    real scalar p, i, depth
    string scalar s, ch
    s = " " + s0
    p = strpos(s, " " + name + "(")
    if (p == 0) return("<none>")
    p = p + strlen(name) + 2
    depth = 1
    for (i = p; i <= strlen(s); i++) {
        ch = substr(s, i, 1)
        if (ch == "(") depth++
        else if (ch == ")") {
            depth--
            if (depth == 0) return(strtrim(substr(s, p, i - p)))
        }
    }
    return("<open>")
}

real scalar _s2b_eval(string scalar expr, real scalar v)
{
    stata("scalar __s2b_v = " + subinstr(expr, "x", "(" + strofreal(v, "%21.0g") + ")"))
    return(st_numscalar("__s2b_v"))
}

real scalar _s2b_rd(real scalar a, real scalar b)
{
    if (a == b) return(0)
    return(abs(a - b) / (abs(b) + 1))
}

string scalar _s2b_gcmp(string scalar ga, string scalar gb, real scalar tol)
{
    string colvector pa, pb
    string scalar oa, ob, sa, sb, ea, eb, ra, rb
    string rowvector ta, tb
    real rowvector na, nb, rna, rnb
    real scalar k, c, j, v
    string rowvector names
    _s2b_split(ga, pa, oa)
    _s2b_split(gb, pb, ob)
    if (rows(pa) != rows(pb)) return("series " + strofreal(rows(pa)) + " vs " + strofreal(rows(pb)))
    if (_s2b_opt(oa, "xline") != _s2b_opt(ob, "xline")) return("xline " + _s2b_opt(oa, "xline") + " vs " + _s2b_opt(ob, "xline"))
    for (k = 1; k <= rows(pa); k++) {
        sa = strtrim(pa[k])
        sb = strtrim(pb[k])
        if (tokens(sa)[1] != tokens(sb)[1]) return("series " + strofreal(k) + " type")
        c = _s2b_comma(sa)
        ea = substr(sa, 1, c - 1)
        ra = substr(sa, c + 1, .)
        c = _s2b_comma(sb)
        eb = substr(sb, 1, c - 1)
        rb = substr(sb, c + 1, .)
        names = ("mcolor", "lcolor", "msymbol")
        for (j = 1; j <= 3; j++) {
            if (_s2b_opt(ra, names[j]) != _s2b_opt(rb, names[j])) {
                if (!(tokens(sa)[1] == "function" & names[j] != "lcolor"))
                    return("series " + strofreal(k) + " " + names[j] + " " + _s2b_opt(ra, names[j]) + " vs " + _s2b_opt(rb, names[j]))
            }
        }
        if (tokens(sa)[1] == "scatteri") {
            na = strtoreal(tokens(ea)[|2 \ .|])
            nb = strtoreal(tokens(eb)[|2 \ .|])
            if (cols(na) != cols(nb)) return("series " + strofreal(k) + " points")
            for (j = 1; j <= cols(na); j++) {
                if (_s2b_rd(na[j], nb[j]) > tol) return("series " + strofreal(k) + " coordinate " + strofreal(j) + ": " + strofreal(na[j], "%21.0g") + " vs " + strofreal(nb[j], "%21.0g"))
            }
        }
        else {
            rna = strtoreal(tokens(_s2b_opt(ra, "range")))
            rnb = strtoreal(tokens(_s2b_opt(rb, "range")))
            if (cols(rna) != 2 | cols(rnb) != 2) return("series " + strofreal(k) + " range")
            if (_s2b_rd(rna[1], rnb[1]) > tol | _s2b_rd(rna[2], rnb[2]) > tol) return("series " + strofreal(k) + " range " + _s2b_opt(ra, "range") + " vs " + _s2b_opt(rb, "range"))
            ea = substr(ea, strpos(ea, "function") + 8, .)
            eb = substr(eb, strpos(eb, "function") + 8, .)
            for (j = 0; j <= 2; j++) {
                v = rna[1] + j * (rna[2] - rna[1]) / 2
                if (_s2b_rd(_s2b_eval(ea, v), _s2b_eval(eb, v)) > tol) return("series " + strofreal(k) + " fit at " + strofreal(v))
            }
        }
    }
    return("")
}

real scalar _s2b_samefile(string scalar a, string scalar b)
{
    real scalar fa, fb, same
    string scalar sa, sb
    fa = fopen(a, "r")
    fb = fopen(b, "r")
    same = 1
    while (1) {
        sa = fread(fa, 1048576)
        sb = fread(fb, 1048576)
        if (rows(sa) != rows(sb) | cols(sa) != cols(sb)) {
            same = 0
            break
        }
        if (!rows(sa)) break
        if (sa != sb) {
            same = 0
            break
        }
    }
    fclose(fa)
    fclose(fb)
    return(same)
}
end

* Run the binscatter and the cbinscatter command, each after the same seed,
* and compare their graphs (see above). r(msg) is empty when they agree.
capture program drop _s2b_graphs
program define _s2b_graphs, rclass
    syntax, REFerence(string asis) CANDidate(string asis) [TOL(real 1e-9) svg SEED(integer 4242)]
    tempfile f
    local svga "`f'_a.svg"
    local svgb "`f'_b.svg"
    set seed `seed'
    capture graph drop _all
    ereturn clear
    quietly `reference'
    local ga `"`e(graphcmd)'"'
    if "`svg'" != "" quietly graph export "`svga'", as(svg) replace
    set seed `seed'
    capture graph drop _all
    ereturn clear
    quietly `candidate'
    local gb `"`e(graphcmd)'"'
    if "`svg'" != "" quietly graph export "`svgb'", as(svg) replace
    mata: st_local("msg", _s2b_gcmp(st_local("ga"), st_local("gb"), `tol'))
    if `"`ga'"' == "" | `"`gb'"' == "" local msg "no e(graphcmd)"
    if "`msg'" == "" & "`svg'" != "" {
        mata: st_local("same", strofreal(_s2b_samefile(st_local("svga"), st_local("svgb"))))
        if !`same' local msg "SVG differs"
    }
    capture erase "`svga'"
    capture erase "`svgb'"
    return local msg `"`msg'"'
end

* One y and one series: binscatter's savedata() CSV (x, y) against
* cbinscatter's e(bindata) at the CSV's precision
capture program drop _s2b_csv
program define _s2b_csv, rclass
    syntax, REFerence(string asis) CANDidate(string asis) [SEED(integer 4242)]
    tempfile f
    set seed `seed'
    ereturn clear
    quietly `reference' savedata(`"`f'"') replace
    set seed `seed'
    ereturn clear
    quietly `candidate'
    tempname cb bs
    matrix `cb' = e(bindata)
    preserve
    quietly import delimited using `"`f'.csv"', clear asdouble
    quietly ds
    local xv : word 1 of `r(varlist)'
    local yv : word 2 of `r(varlist)'
    mkmat `xv' `yv', matrix(`bs')
    restore
    capture erase `"`f'.csv"'
    capture erase `"`f'.do"'
    local msg
    if rowsof(`bs') != rowsof(`cb') local msg "bins `=rowsof(`bs')' vs `=rowsof(`cb')'"
    else {
        forvalues r = 1/`=rowsof(`bs')' {
            if reldif(`bs'[`r', 1], `cb'[`r', 3]) > 1e-6 | reldif(`bs'[`r', 2], `cb'[`r', 4]) > 1e-6 {
                local msg "bin `r'"
                continue, break
            }
        }
    }
    return local msg `"`msg'"'
end

capture program drop _s2b_data
program define _s2b_data
    clear
    set seed 917
    quietly set obs 1200
    gen double x = rnormal()
    gen double z = rnormal() + .3*x
    gen double y = 1 + x + .5*x^2 + z + rnormal()
    gen double y2 = 2 - x + rnormal()
    gen double w = runiform()*3
    gen f = 1 + mod(_n, 4)
    gen q = ceil(_n/100)
    gen double u = runiform()
    gen g = mod(_n, 3)
    gen a = mod(_n, 7)
    gen xd = mod(_n, 6)
    gen double yr = round(y, .5)
end

* ===========================================================================
* BIN-02/BIN-03: graphs (bins, fits, colors, rd lines) against binscatter
* ===========================================================================
_s2b_data
local fails
* Linear fits without residualization: the same graph, byte for byte
foreach o in "medians" "medians by(g)" "[aw=w]|medians" "[fw=f]|medians" "rd(0)" ///
    "rd(-1 .5) by(g)" "rd(-5 0)" "xq(q)" "xq(q) by(g)" "xq(q) medians" ///
    "randvar(u) randcut(.3)" "randcut(.4)" "randn(150)" "[aw=w]|randcut(.5) rd(0)" ///
    "randcut(.5) by(g)" "mcolors(red blue) lcolors(green)" ///
    "by(g) mcolors(red blue green) lcolors(black gs8)" "linetype(connect) mcolors(red) lcolors(blue)" ///
    "medians linetype(connect) by(g)" "noaddmean" {
    gettoken wt opts : o, parse("|")
    if `"`opts'"' == "" {
        local opts `"`o'"'
        local wt
    }
    else local opts = substr(`"`opts'"', 2, .)
    capture noisily _s2b_graphs, reference(binscatter y x `wt', `opts') ///
        candidate(cbinscatter y x `wt', `opts') svg
    if _rc local fails `"`fails' [`o': r(`=_rc')]"'
    else if `"`r(msg)'"' != "" local fails `"`fails' [`o': `r(msg)']"'
}
if `"`fails'"' != "" test_fail "BIN-02 graph equals binscatter's (SVG)" `"`fails'"'
else test_pass "BIN-02 graph equals binscatter's (SVG)"

* Several y, with the options (common bins, random subsample and colors)
local fails
foreach o in "medians" "rd(0 1)" "xq(q)" "randcut(.5)" "randn(200) by(g)" ///
    "by(g) rd(0) mcolors(red blue green) lcolors(black gs8)" "medians xq(q) by(g)" {
    capture noisily _s2b_graphs, reference(binscatter y y2 x, `o') candidate(cbinscatter y y2 x, `o') svg
    if _rc local fails `"`fails' [`o': r(`=_rc')]"'
    else if `"`r(msg)'"' != "" local fails `"`fails' [`o': `r(msg)']"'
}
if `"`fails'"' != "" test_fail "BIN-02 several y equal binscatter's (SVG)" `"`fails'"'
else test_pass "BIN-02 several y equal binscatter's (SVG)"

* With controls()/absorb(), binscatter stores the residuals as float, and its
* qfit regresses on a float x^2: compared at 1e-6
local fails
foreach o in "medians controls(z)" "noaddmean controls(z)" "noaddmean absorb(a)" ///
    "noaddmean controls(z) by(g) rd(0)" "medians absorb(a) controls(z) rd(.2)" ///
    "xq(q) absorb(a)" "randcut(.5) controls(z)" "rd(-1 .5) linetype(qfit)" ///
    "rd(0) linetype(qfit) controls(z) by(g)" "medians linetype(qfit)" {
    capture noisily _s2b_graphs, reference(binscatter y x, `o') candidate(cbinscatter y x, `o') tol(1e-6)
    if _rc local fails `"`fails' [`o': r(`=_rc')]"'
    else if `"`r(msg)'"' != "" local fails `"`fails' [`o': `r(msg)']"'
    capture noisily _s2b_graphs, reference(binscatter y y2 x, `o') candidate(cbinscatter y y2 x, `o') tol(1e-6)
    if _rc local fails `"`fails' [y y2 `o': r(`=_rc')]"'
    else if `"`r(msg)'"' != "" local fails `"`fails' [y y2 `o': `r(msg)']"'
}
if `"`fails'"' != "" test_fail "BIN-02 controls/absorb/qfit equal binscatter's" `"`fails'"'
else test_pass "BIN-02 controls/absorb/qfit equal binscatter's"

* Bin coordinates from binscatter's savedata() CSV
local fails
foreach o in "medians" "xq(q)" "randcut(.3)" "randvar(u) randcut(.6)" "medians discrete" {
    local x x
    if strpos("`o'", "discrete") local x xd
    capture noisily _s2b_csv, reference(binscatter yr `x', `o' nodraw) candidate(cbinscatter yr `x', `o' nograph)
    if _rc local fails `"`fails' [`o': r(`=_rc')]"'
    else if `"`r(msg)'"' != "" local fails `"`fails' [`o': `r(msg)']"'
}
if `"`fails'"' != "" test_fail "BIN-02 savedata CSV bins equal binscatter's" `"`fails'"'
else test_pass "BIN-02 savedata CSV bins equal binscatter's"

* medians with ties and even counts: _pctile's average of the middle values
* (yr is rounded to .5), weighted ties included
capture noisily {
    _s2b_data
    gen fw2 = 1 + (_n > 600)
    foreach wt in "" "[fw=fw2]" "[aw=fw2]" {
        tempvar b
        quietly binscatter yr x `wt', medians genxq(`b') nodraw
        quietly cbinscatter yr x `wt', medians nograph nquantiles(20)
        tempname cb
        matrix `cb' = e(bindata)
        forvalues k = 1/20 {
            quietly _pctile yr `wt' if `b' == `k', p(50)
            assert reldif(r(r1), `cb'[`k', 4]) < 1e-12
            quietly _pctile x `wt' if `b' == `k', p(50)
            assert reldif(r(r1), `cb'[`k', 3]) < 1e-12
        }
        drop `b'
    }
}
if _rc test_fail "BIN-02 medians are _pctile medians of binscatter's bins" "r(`=_rc')"
else test_pass "BIN-02 medians are _pctile medians of binscatter's bins"

* rd(): fit coefficients by interval, e(rdintervals), reportreg, the message
* for a side without observations and r(2001) for a side with one
capture noisily {
    _s2b_data
    foreach o in "rd(0)" "rd(-1 .5) by(g)" "rd(0) controls(z)" "rd(-5 0) linetype(qfit)" {
        quietly binscatter y x, `o' nodraw
        tempname bsb bsr
        matrix `bsb' = e(y1_coefs)
        matrix `bsr' = e(rdintervals)
        quietly cbinscatter y x, `o' nograph
        assert mreldif(e(rdintervals), `bsr') == 0
        assert rowsof(e(coefs)) == rowsof(`bsb')
        local ord = colsof(`bsb') - 1
        forvalues r = 1/`=rowsof(`bsb')' {
            if missing(`bsb'[`r', 1]) {
                assert missing(e(coefs)[`r', 1])
                continue
            }
            * binscatter: x^2, x, _cons; cbinscatter: _cons, x, x^2
            forvalues j = 0/`ord' {
                assert reldif(`bsb'[`r', `ord' + 1 - `j'], e(coefs)[`r', `j' + 1]) < 1e-6
            }
        }
    }
    quietly cbinscatter y x, rd(0) reportreg controls(z) nograph
    tempname rb cf
    matrix `rb' = e(reportreg_b)
    matrix `cf' = e(coefs)
    assert rowsof(`rb') == 2
    assert reldif(`rb'[1, 1], `cf'[1, 2]) < 1e-6 & reldif(`rb'[2, 1], `cf'[2, 2]) < 1e-6
    sort x
    local c1 = (x[1] + x[2]) / 2
    capture binscatter y x, rd(`c1') nodraw
    local brc = _rc
    capture cbinscatter y x, rd(`c1') nograph
    assert _rc == `brc' & `brc' == 2001
}
if _rc test_fail "BIN-02 rd() fits, e(rdintervals) and errors match binscatter" "r(`=_rc')"
else test_pass "BIN-02 rd() fits, e(rdintervals) and errors match binscatter"

* savegraph(): .gph saved, other extensions exported, replace needed
capture noisily {
    tempfile sg
    * graph export takes the format from the extension: drop the .tmp
    local sg = substr(`"`sg'"', 1, strlen(`"`sg'"') - 4) + "_s2b"
    local gph "`sg'.gph"
    local png "`sg'.png"
    _s2b_data
    cbinscatter y x, savegraph("`gph'") nodraw
    confirm file "`gph'"
    capture cbinscatter y x, savegraph("`gph'") nodraw
    local crc = _rc
    capture binscatter y x, savegraph("`gph'") nodraw
    assert _rc == 602 & `crc' == 602
    cbinscatter y x, savegraph("`gph'") replace nodraw
    cbinscatter y y2 x, savegraph("`png'")
    confirm file "`png'"
    capture cbinscatter y x, savegraph("`sg'_new.gph") nograph
    assert _rc == 198
    capture confirm file "`sg'_new.gph"
    assert _rc
    erase "`gph'"
    erase "`png'"
}
if _rc test_fail "BIN-02 savegraph() as binscatter's" "r(`=_rc')"
else test_pass "BIN-02 savegraph() as binscatter's"

* The random subsample: the same random-number state afterwards, none drawn
* when x is discrete, and fastxtile's errors
capture noisily {
    _s2b_data
    set seed 33
    quietly binscatter y x, randcut(.5) nodraw
    local after_bs = runiform()
    set seed 33
    quietly cbinscatter y x, randcut(.5) nograph
    assert runiform() == `after_bs'
    set seed 33
    quietly cbinscatter y y2 x, randcut(.5) nograph
    assert runiform() == `after_bs'
    local before `"`c(rngstate)'"'
    quietly cbinscatter y xd, randcut(.5) nograph
    assert `"`c(rngstate)'"' == `"`before'"'
    foreach o in "randn(0)" "randcut(1.5)" "randcut(.5) randn(10)" "randvar(y) randn(10)" ///
        "nquantiles(500) randcut(.2)" {
        capture binscatter y x, `o' nodraw
        local brc = _rc
        capture cbinscatter y x, `o' nograph
        assert _rc == `brc' & `brc' == 198
    }
    * binscatter stops with r(0) after these messages; cbinscatter r(198)
    foreach o in "nofastxtile randcut(.5)" "xq(q) nquantiles(10)" "xq(q) genxq(b)" "xq(q) discrete" {
        capture cbinscatter y x, `o' nograph
        assert _rc == 198
    }
    gen double qbad = q + .5
    capture cbinscatter y x, xq(qbad) nograph
    assert _rc == 198
    gen qmany = _n
    capture binscatter y x, xq(qmany) nodraw
    local brc = _rc
    capture cbinscatter y x, xq(qmany) nograph
    assert _rc == `brc' & `brc' == 134
    foreach o in medians noaddmean "rd(0)" "xq(q)" "randcut(.5)" "randvar(u)" "randn(10)" {
        capture cbinscatter y x, `o' method(binsreg) nograph
        assert _rc == 198
    }
}
if _rc test_fail "BIN-02 random subsample state and option errors" "r(`=_rc')"
else test_pass "BIN-02 random subsample state and option errors"

* BIN-03: genxq() equals binscatter's variable
local fails
_s2b_data
foreach o in "" "by(g)" "controls(z)" "absorb(a) by(g)" "[aw=w]|" "[fw=f]|by(g)" ///
    "randcut(.3)" "randvar(u) randcut(.5) by(g)" "nquantiles(80)" "medians rd(0)" {
    gettoken wt opts : o, parse("|")
    if `"`opts'"' == "" {
        local opts `"`o'"'
        local wt
    }
    else local opts = substr(`"`opts'"', 2, .)
    capture drop bq cq
    set seed 51
    quietly binscatter y x `wt', `opts' genxq(bq) nodraw
    set seed 51
    capture quietly cbinscatter y x `wt', `opts' genxq(cq) nograph
    if _rc {
        local fails `"`fails' [`o': r(`=_rc')]"'
        continue
    }
    quietly count if bq != cq
    if r(N) local fails `"`fails' [`o': `r(N)' differ]"'
    if "`: type bq'" != "`: type cq'" local fails `"`fails' [`o': type `: type cq'']"'
}
* Several y: one variable, from the common bins
capture drop bq cq
quietly binscatter y y2 x, genxq(bq) nodraw
quietly cbinscatter y y2 x, genxq(cq) nograph
quietly count if bq != cq
if r(N) local fails `"`fails' [y y2: `r(N)' differ]"'
* x with few values: discrete, no variable (with binscatter's note)
quietly binscatter y xd, genxq(bd) nodraw
quietly cbinscatter y xd, genxq(cd) nograph
capture confirm variable bd
local brc = _rc
capture confirm variable cd
if _rc != 111 | `brc' != 111 local fails `"`fails' [discrete x: variable made]"'
capture cbinscatter y x, genxq(y) nograph
if _rc != 110 local fails `"`fails' [existing name: r(`=_rc')]"'
if `"`fails'"' != "" test_fail "BIN-03 genxq() equals binscatter's" `"`fails'"'
else test_pass "BIN-03 genxq() equals binscatter's"

* ===========================================================================
* WINS-04: weights, as gstats winsor computes the weighted cutoffs
* ===========================================================================
clear
set seed 2604
quietly set obs 3000
gen double x = rnormal()
replace x = round(x, .1) in 1/800
replace x = . in 1/25
gen float xf = rnormal()*10
gen g = mod(_n, 7)
replace g = . in 30/45
gen double aw = runiform()*3
replace aw = 0 in 50/60
gen fw = 1 + mod(_n, 4)
gen double iw = runiform()
gen double pw = runiform() + .1
* tied cumulative weights at the cut: 0.1 aweights
gen double tie = .1
local fails
local n = 0
foreach w in "aw=aw" "fw=fw" "iw=iw" "pw=pw" "aw=tie" {
    foreach c in "1 99" "5 95" "10 90" "25 75" "50 50" "0.5 99.5" "20 80" {
        foreach t in "" "trim" {
            foreach b in "" "by(g)" "by(g) nomissing" {
                foreach v in x xf {
                    capture drop `v'_g `v'_c
                    quietly gstats winsor `v' [`w'], cuts(`c') `t' `b' gen(`v'_g)
                    local J = r(J)
                    local minJ = r(minJ)
                    local maxJ = r(maxJ)
                    capture quietly cwinsor `v' [`w'], cuts(`c') `t' `b' generate(`v'_c)
                    if _rc {
                        local fails `"`fails' [`w' `c' `t' `b' `v': r(`=_rc')]"'
                        continue
                    }
                    local ok = r(J) == `J' & r(minJ) == `minJ' & r(maxJ) == `maxJ'
                    * winsorized by() targets are float in cwinsor (winsor2's type)
                    if "`b'" == "" | "`t'" == "trim" quietly count if `v'_g != `v'_c & !(missing(`v'_g) & missing(`v'_c))
                    else quietly count if float(`v'_g) != `v'_c & !(missing(`v'_g) & missing(`v'_c))
                    if r(N) | !`ok' local fails `"`fails' [`w' `c' `t' `b' `v': `r(N)' differ, J `ok']"'
                    local ++n
                }
            }
        }
    }
}
if `"`fails'"' != "" test_fail "WINS-04 weighted cutoffs equal gstats winsor's" `"`fails'"'
else test_pass "WINS-04 weighted cutoffs equal gstats winsor's (`n' cases)"

* In place (cwinsor's default) with weights gives the generated values;
* zero-weight observations are left unchanged; errors as gstats
capture noisily {
    gen double x1 = x
    quietly cwinsor x1 [aw=aw], cuts(5 95) by(g)
    capture drop x_g
    quietly gstats winsor x [aw=aw], cuts(5 95) by(g) gen(x_g)
    assert x1 == x_g | (missing(x_g) & x1 == x) | (missing(x1) & missing(x_g))
    assert x1 == x in 50/60
    gen double wn = -aw
    capture gstats winsor x [aw=wn], gen(neg)
    local grc = _rc
    capture cwinsor x [aw=wn], generate(neg)
    assert _rc == `grc' & `grc' == 402
    capture cwinsor x [iw=wn], generate(neg)
    assert _rc == 402
    gen double fwx = 1.5
    capture gstats winsor x [fw=fwx], gen(neg)
    local grc = _rc
    capture cwinsor x [fw=fwx], generate(neg)
    assert _rc == `grc' & `grc' == 401
}
if _rc test_fail "WINS-04 in place and weight errors" "r(`=_rc')"
else test_pass "WINS-04 in place and weight errors"

* WINS-05 generate(), WINS-06 replace with targets, WINS-07 nomissing,
* WINS-08 r(): against gstats winsor
capture noisily {
    sysuse auto, clear
    replace rep78 = . in 1/5
    gen double pr = price
    gen double mp = mpg
    quietly gstats winsor pr mp, gen(gp gm) by(rep78) nomissing cuts(5 95)
    foreach r in J minJ maxJ cutlow cuthigh {
        local g_`r' = r(`r')
    }
    quietly cwinsor pr mp, generate(cp cm) by(rep78) nomissing cuts(5 95)
    foreach r in J minJ maxJ cutlow cuthigh {
        assert r(`r') == `g_`r''
    }
    assert (float(gp) == cp | missing(gp) & missing(cp)) & (float(gm) == cm | missing(gm) & missing(cm))
    assert missing(cp) if missing(rep78)
    quietly gstats winsor price, gen(gp0)
    local gJ = r(J)
    local gmin = r(minJ)
    quietly cwinsor price, generate(cp0)
    assert r(J) == `gJ' & r(minJ) == `gmin' & r(maxJ) == _N
    * replace overwrites existing targets (keeping their place); without
    * replace an existing one is r(110)
    gen byte price_w = 1
    gen byte p_mpg = 1
    order price_w, first
    quietly gstats winsor price if foreign, suffix(_g) cuts(10 90)
    quietly cwinsor price if foreign, suffix(_w) replace cuts(10 90)
    assert reldif(price_w, price_g) < 1e-12 | (missing(price_w) & missing(price_g))
    confirm int variable price_w
    unab all : _all
    assert "`: word 1 of `all''" == "price_w"
    capture cwinsor price, suffix(_w)
    assert _rc == 110
    quietly cwinsor mpg, prefix(p_) replace
    capture cwinsor price mpg, generate(one)
    assert _rc == 198
    capture cwinsor price, generate(one) suffix(_x)
    assert _rc == 198
    capture cwinsor price mpg, generate(mpg two) replace
    assert _rc == 198
}
if _rc test_fail "WINS-05..08 generate(), replace, nomissing, r() as gstats winsor" "r(`=_rc')"
else test_pass "WINS-05..08 generate(), replace, nomissing, r() as gstats winsor"

* WINS-01: label, as winsor2 writes it
capture noisily {
    sysuse auto, clear
    label variable mpg ""
    * winsor2 adds another 0 before a lower percentile below 1 for each
    * further variable (a bug cwinsor does not copy): one variable at a time
    foreach o in "cuts(0.5 99)" "cuts(5 95) trim" "cuts(0 90)" "cuts(1 99) by(foreign)" {
        foreach v in price mpg {
            quietly winsor2 `v', `o' label suffix(_v)
            quietly cwinsor `v', `o' label suffix(_c)
            assert `"`: variable label `v'_v'"' == `"`: variable label `v'_c'"'
        }
        drop *_v *_c
    }
    quietly winsor2 price, cuts(5 95) label replace
    local wl : variable label price
    sysuse auto, clear
    quietly cwinsor price, cuts(5 95) label
    assert `"`: variable label price'"' == `"`wl'"'
}
if _rc test_fail "WINS-01 label as winsor2's" "r(`=_rc')"
else test_pass "WINS-01 label as winsor2's"

* ---------------------------------------------------------------------------
* BIN-04/06/07/08/09: cbinscatter against binscatter 7.02. _s2b_s3 runs
* binscatter and then cbinscatter with the same arguments on the same data.
* It compares the return code, the variables either one leaves behind,
* e(N), e(sample), e(y#_coefs) (values within 1e-6 and row and column
* names), e(byvalues) and e(rdintervals), and, with savedata(), the .csv
* (header, rows, values within 1e-6) and the .do (every line, with the
* numbers masked: binscatter takes the fits from regress).
* ---------------------------------------------------------------------------
capture program drop _s2b_s3
program define _s2b_s3
    gettoken name 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    local name = strtrim(`"`name'"')
    local save = strpos(`"`0'"', "savedata(_s2b_sd)") > 0
    capture erase _s2b_sd.csv
    capture erase _s2b_sd.do
    capture erase _s2b_sd.dta
    preserve
    quietly ds
    local vars0 `r(varlist)'
    capture noisily binscatter `0'
    local rc_ref = _rc
    quietly ds
    local vars_ref `r(varlist)'
    tempvar es_ref
    capture generate byte `es_ref' = e(sample)
    local N_ref = e(N)
    local mats_ref : e(matrices)
    foreach m of local mats_ref {
        tempname R_`m'
        matrix `R_`m'' = e(`m')
    }
    if `save' & !`rc_ref' {
        quietly copy _s2b_sd.csv _s2b_ref.csv, replace
        quietly copy _s2b_sd.do _s2b_ref.do, replace
        erase _s2b_sd.csv
        erase _s2b_sd.do
    }
    tempfile refdata
    quietly save `refdata'
    restore, preserve
    capture noisily cbinscatter `0'
    local rc_ct = _rc
    local bad
    if `rc_ref' != `rc_ct' local bad rc:`rc_ref'/`rc_ct'
    else if !`rc_ref' {
        quietly ds
        local vars_ct `r(varlist)'
        local new_ct : list vars_ct - vars0
        local new_ref : list vars_ref - vars0
        local new_ref : list new_ref - es_ref
        if "`new_ct'" != "`new_ref'" local bad `bad' variables:[`new_ref']/[`new_ct']
        foreach v of local new_ct {
            rename `v' _ct_`v'
        }
        if e(N) != `N_ref' local bad `bad' e(N):`N_ref'/`=e(N)'
        tempvar es_ct
        generate byte `es_ct' = e(sample)
        local both : list new_ct & new_ref
        quietly merge 1:1 _n using `refdata', keepusing(`es_ref' `both') nogenerate
        quietly count if `es_ct' != `es_ref'
        if r(N) local bad `bad' e(sample):`r(N)'
        foreach v of local both {
            quietly count if _ct_`v' != `v'
            if r(N) local bad `bad' `v':`r(N)'
        }
        foreach m of local mats_ref {
            capture confirm matrix e(`m')
            if _rc {
                local bad `bad' e(`m'):absent
                continue
            }
            tempname C
            matrix `C' = e(`m')
            if rowsof(`C') != rowsof(`R_`m'') | colsof(`C') != colsof(`R_`m'') {
                local bad `bad' e(`m'):size
                continue
            }
            if `"`: rowfullnames `C''"' != `"`: rowfullnames `R_`m'''"' local bad `bad' e(`m'):rownames
            if `"`: colfullnames `C''"' != `"`: colfullnames `R_`m'''"' local bad `bad' e(`m'):colnames
            mata: st_local("diff", strofreal(_s2b_matdiff(st_matrix("`C'"), st_matrix("`R_`m''"))))
            if `diff' > 1e-6 local bad `bad' e(`m'):`diff'
        }
        if `save' {
            capture confirm file _s2b_sd.csv
            if _rc local bad `bad' csv:absent
            else {
                mata: st_local("csv", _s2b_csvdiff("_s2b_ref.csv", "_s2b_sd.csv", 1e-6))
                if "`csv'" != "" local bad `bad' csv:`csv'
                mata: st_local("dof", _s2b_dodiff("_s2b_ref.do", "_s2b_sd.do"))
                if "`dof'" != "" local bad `bad' do:`dof'
            }
        }
    }
    restore
    capture erase _s2b_sd.csv
    capture erase _s2b_sd.do
    capture erase _s2b_sd.dta
    capture erase _s2b_ref.csv
    capture erase _s2b_ref.do
    if "`bad'" == "" & `rc_ref' == 0 test_pass "`name'"
    else if "`bad'" == "" test_pass "`name': same error r(`rc_ref')"
    else test_fail "`name'" "`bad'"
end
capture mata: mata drop _s2b_matdiff()
capture mata: mata drop _s2b_csvdiff()
capture mata: mata drop _s2b_dodiff()
capture mata: mata drop _s2b_lines()
mata:
real scalar _s2b_matdiff(real matrix A, real matrix B)
{
    real scalar i, j, d
    d = 0
    for (i = 1; i <= rows(A); i++) {
        for (j = 1; j <= cols(A); j++) {
            if (A[i, j] >= . & B[i, j] >= .) continue
            if (A[i, j] >= . | B[i, j] >= .) return(.)
            d = max((d, reldif(A[i, j], B[i, j])))
        }
    }
    return(d)
}
string colvector _s2b_lines(string scalar f)
{
    real scalar fh
    string colvector L
    string matrix line
    L = J(0, 1, "")
    fh = fopen(f, "r")
    while ((line = fget(fh)) != J(0, 0, "")) L = L \ line
    fclose(fh)
    return(L)
}
string scalar _s2b_csvdiff(string scalar a, string scalar b, real scalar tol)
{
    string colvector A, B
    string rowvector ra, rb
    real scalar i, j, x, y
    A = _s2b_lines(a)
    B = _s2b_lines(b)
    if (rows(A) != rows(B)) return("rows " + strofreal(rows(A)) + "/" + strofreal(rows(B)))
    if (A[1] != B[1]) return("header")
    for (i = 2; i <= rows(A); i++) {
        ra = tokens(subinstr(A[i], ",", " "))
        rb = tokens(subinstr(B[i], ",", " "))
        if (cols(ra) != cols(rb)) return("row " + strofreal(i))
        for (j = 1; j <= cols(ra); j++) {
            x = strtoreal(ra[j])
            y = strtoreal(rb[j])
            if (x >= . & y >= .) continue
            if (x >= . | y >= . | reldif(x, y) > tol) return("row " + strofreal(i) + " col " + strofreal(j))
        }
    }
    return("")
}
string scalar _s2b_dodiff(string scalar a, string scalar b)
{
    string colvector A, B
    real scalar i
    string scalar pat
    A = _s2b_lines(a)
    B = _s2b_lines(b)
    if (rows(A) != rows(B)) return("lines " + strofreal(rows(A)) + "/" + strofreal(rows(B)))
    pat = "[0-9]*[.]?[0-9]+(e[-+]?[0-9]+)?"
    for (i = 1; i <= rows(A); i++) {
        if (ustrregexra(A[i], pat, "N") != ustrregexra(B[i], pat, "N")) return("line " + strofreal(i))
    }
    return("")
}
end

sysuse auto, clear
gen str5 fs = cond(foreign, "f", "d")
label define rlab 1 "one" 2 "two" 3 "three" 4 "four" 5 "five"
gen rl = rep78
label values rl rlab
* BIN-06: e(sample), e(y#_coefs), e(byvalues), e(rdintervals)
_s2b_s3 BIN-06 lfit: price mpg
foreach o in "linetype(qfit)" "rd(25)" "by(foreign)" "by(foreign) rd(20 30) linetype(qfit)" ///
    "by(fs)" "by(rl)" "controls(weight)" "absorb(rep78)" "linetype(connect)" "linetype(none)" ///
    "discrete" "nquantiles(5)" {
    _s2b_s3 BIN-06 `o': price mpg if price < 14000, `o'
}
_s2b_s3 BIN-06 aweights: price mpg [aw=weight], nquantiles(5)
_s2b_s3 BIN-06 fweights: price mpg [fw=rep78], nquantiles(8) by(foreign)
_s2b_s3 BIN-06 several y: price weight mpg, by(foreign) linetype(qfit) rd(25)
_s2b_s3 BIN-06 several y, by(string): price weight mpg, by(fs) nquantiles(10)
* BIN-07: binscatter's abbreviations n() and m()/ms()
foreach o in "n(10)" "n(5) by(foreign)" "m(Oh)" "ms(Oh D)" "n(10) ms(Oh) by(foreign)" {
    _s2b_s3 BIN-07 `o': price mpg, `o'
}
_s2b_s3 BIN-07 several y n() ms(): price weight mpg, n(8) ms(Oh D)
* BIN-08: binscatter's legacy options, their notes, and both forms given
foreach o in "nbins(10)" "create_xq" "x_q(rep78)" "symbols(Oh)" "unique(x)" "linetype(noline)" ///
    "nbins(10) create_xq symbols(D) linetype(noline)" {
    capture drop q_mpg
    _s2b_s3 BIN-08 `o': price mpg, `o'
}
capture drop q_mpg
_s2b_s3 BIN-08 several y nbins() symbols(): price weight mpg, nbins(6) symbols(Oh D)
* With both the option and its legacy form binscatter shows its message and
* stops with r(0); cbinscatter shows the same message with r(198)
foreach o in "nbins(10) n(5)" "create_xq genxq(z)" "x_q(rep78) xq(rep78)" "symbols(Oh) ms(D)" {
    capture noisily binscatter price mpg, `o'
    local rc_ref = _rc
    capture noisily cbinscatter price mpg, `o'
    local rc_ct = _rc
    capture confirm variable z
    local made = !_rc
    if `rc_ref' == 0 & `rc_ct' == 198 & !`made' test_pass "BIN-08 `o' refused (binscatter r(0), cbinscatter r(198))"
    else test_fail "BIN-08 `o' refused" "binscatter r(`rc_ref') cbinscatter r(`rc_ct')"
}
* BIN-09: nquantiles() above 1000. With fewer distinct x values binscatter
* bins x as discrete, as cbinscatter does; above 100 quantile bins binscatter
* itself stops ("32,740 invalid name", r(198)) where cbinscatter forms them
_s2b_s3 BIN-09 nquantiles(1500), discrete x: price mpg, nquantiles(1500) savedata(_s2b_sd) replace
_s2b_s3 BIN-09 nquantiles(5000) by(), discrete x: price mpg, nquantiles(5000) by(foreign)
preserve
clear
quietly set obs 30000
gen double x = runiform()
gen double y = x + rnormal()
capture binscatter y x, nquantiles(1500)
local rc_ref = _rc
capture cbinscatter y x, nquantiles(1500) nograph
local rc_ct = _rc
local nb = 0
if !`rc_ct' {
    mata: st_local("nb", strofreal(sum(st_matrix("e(bindata)")[., 2] :< .)))
}
if `rc_ref' == 198 & !`rc_ct' & `nb' == 1500 test_pass "BIN-09 nquantiles(1500) quantile bins (binscatter r(198))"
else test_fail "BIN-09 nquantiles(1500) quantile bins" "binscatter r(`rc_ref') cbinscatter r(`rc_ct') bins `nb'"
restore
sysuse auto, clear
gen str5 fs = cond(foreign, "f", "d")
* BIN-04: savedata() writes binscatter's .csv and .do (and the bins in .dta)
foreach o in "" "linetype(qfit) rd(25)" "by(foreign)" "by(fs) discrete" "by(foreign) linetype(qfit) rd(20 30)" ///
    "linetype(connect) ms(Oh)" "linetype(none) colors(red blue)" "mcolors(red) lcolors(green)" ///
    "controls(weight) nquantiles(10)" "xq(rep78)" "n(8) title(T) ytitle(Y)" {
    _s2b_s3 BIN-04 savedata `o': price mpg, savedata(_s2b_sd) replace `o'
}
_s2b_s3 BIN-04 savedata several y: price weight mpg, savedata(_s2b_sd) replace by(foreign) linetype(qfit) rd(25)
_s2b_s3 BIN-04 savedata several y, connect: price weight mpg, savedata(_s2b_sd) replace linetype(connect) n(10)
* Without replace an existing .csv or .do is refused, as in binscatter
capture noisily {
    capture erase _s2b_sd.dta
    capture erase _s2b_sd.do
    quietly copy "validation/stata/s2_bin.do" _s2b_sd.csv, replace
    capture binscatter price mpg, savedata(_s2b_sd)
    local rc_ref = _rc
    capture cbinscatter price mpg, savedata(_s2b_sd)
    local rc_ct = _rc
    capture confirm file _s2b_sd.dta
    local made = !_rc
    erase _s2b_sd.csv
    assert `rc_ref' == 602 & `rc_ct' == 602 & !`made'
    quietly cbinscatter price mpg, savedata(_s2b_sd) nograph
    confirm file _s2b_sd.csv
    confirm file _s2b_sd.do
    confirm file _s2b_sd.dta
    capture binscatter price mpg, savedata(_s2b_sd)
    local rc_ref = _rc
    capture cbinscatter price mpg, savedata(_s2b_sd) nograph
    assert `rc_ref' == 602 & _rc == 602
    * The .do recreates binscatter's dataset and graph
    preserve
    clear
    quietly do _s2b_sd.do
    confirm variable mpg price
    assert _N == 15
    restore
    * a name with an extension is the .dta file alone
    erase _s2b_sd.csv
    erase _s2b_sd.do
    erase _s2b_sd.dta
    quietly cbinscatter price mpg, savedata(_s2b_sd.dta) nograph
    confirm file _s2b_sd.dta
    capture confirm file _s2b_sd.csv
    assert _rc
    erase _s2b_sd.dta
}
if _rc test_fail "BIN-04 savedata files, replace and the .do" "r(`=_rc')"
else test_pass "BIN-04 savedata files, replace and the .do"
capture erase _s2b_sd.csv
capture erase _s2b_sd.do
capture erase _s2b_sd.dta

print_summary s2_bin
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_bin"
    global CTOOLS_CASE_COMPLETE "s2_bin"
    exit
}

di as error "Unknown s2_bin case: `__case'"
exit 198
