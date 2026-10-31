*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! ctools: e(b) and e(V) from a plugin's compact results
* _ctools_ebv <b> <V>, names() vmatrix() [base() stub() kkeep() bmatrix()
*     collinear() cons() kept() colnames()]
*
* names() lists every term of the model in e(b) order and base() flags the
* base/omitted placeholders among them (_ctools_fvsplit's r(base); default:
* none). The other terms are the plugin's columns, in order. Two layouts of
* the plugin's results are understood:
*
*   stub(<s>) kkeep(#)   the plugin dropped its collinear columns: scalars
*                        <s>_collinear_<k> flag each column, <s>_beta_<i> hold
*                        the kept coefficients and vmatrix() their covariance
*                        in its top-left block. With kkeep(0) every column
*                        counts as dropped.
*   bmatrix(<m>)         the plugin kept a (zero) row and column for every
*                        column: bmatrix() and vmatrix() cover them all and
*                        collinear() lists a 0/1 flag per column.
*
* cons(<scalar>) appends _cons with that coefficient; vmatrix() then carries
* the constant in the row and column after the coefficients.
*
* <b> and <V> are created with zeros for placeholders and dropped columns and
* named: a placeholder keeps its name, a dropped column gets the o. prefix.
* kept() names a matrix for the kept coefficients alone (stub layout, not
* created when there are none) and colnames() a caller local for the names.
program define _ctools_ebv
    version 14.1
    syntax namelist(min=2 max=2), NAMES(string) Vmatrix(name) ///
        [BASE(string) STUB(string) KKeep(integer -1) Bmatrix(name) ///
         COLLinear(string) CONS(name) KEPT(name) COLNames(name)]
    gettoken b V : namelist
    if ("`stub'" == "") == ("`bmatrix'" == "") | ("`stub'" != "" & `kkeep' < 0) {
        di as error "_ctools_ebv: specify stub() with kkeep(), or bmatrix()"
        exit 198
    }
    mata: _ctools_ebv_build("`b'", "`V'")
    matrix colnames `b' = `outnames'
    matrix colnames `V' = `outnames'
    matrix rownames `V' = `outnames'
    if "`colnames'" != "" c_local `colnames' `outnames'
end

capture mata: mata drop _ctools_ebv_build()
capture mata: mata drop _ctools_ebv_scalar()
mata:
real scalar _ctools_ebv_scalar(string scalar name)
{
    real matrix value
    value = st_numscalar(name)
    if (rows(value) != 1 | cols(value) != 1) {
        errprintf("ctools: plugin result %s not found\n", name)
        exit(498)
    }
    return(value)
}

void _ctools_ebv_build(string scalar bname, string scalar Vname)
{
    string rowvector names
    string scalar stub
    real rowvector base, sent, coll, kept, pos, bsrc, b
    real matrix Vsrc, V
    real scalar T, Ks, Kk, m, n, i, hascons

    names = tokens(st_local("names"))
    T = cols(names)
    base = (st_local("base") == "" ? J(1, T, 0) : strtoreal(tokens(st_local("base"))))
    if (cols(base) != T | missing(base)) {
        errprintf("ctools: %g terms but %g base flags\n", T, cols(base))
        exit(498)
    }
    Ks = T - sum(base :!= 0)
    sent = (Ks ? selectindex(base :== 0) : J(1, 0, .))

    stub = st_local("stub")
    if (stub == "") {
        /* a row and column for every plugin column */
        coll = (st_local("collinear") == "" ? J(1, Ks, 0) : strtoreal(tokens(st_local("collinear"))))
        if (cols(coll) != Ks | missing(coll)) {
            errprintf("ctools: %g plugin columns but %g collinearity flags\n", Ks, cols(coll))
            exit(498)
        }
        kept = (Ks ? 1..Ks : J(1, 0, .))
        bsrc = st_matrix(st_local("bmatrix"))
    }
    else {
        /* collinear columns dropped by the plugin */
        coll = J(1, Ks, 1)
        if (strtoreal(st_local("kkeep")) > 0) {
            for (i = 1; i <= Ks; i++) {
                coll[i] = (_ctools_ebv_scalar(stub + "_collinear_" + strofreal(i)) != 0)
            }
        }
        Kk = Ks - sum(coll)
        kept = (Kk ? selectindex(coll :== 0) : J(1, 0, .))
        bsrc = J(1, Kk, 0)
        for (i = 1; i <= Kk; i++) {
            bsrc[i] = _ctools_ebv_scalar(stub + "_beta_" + strofreal(i))
        }
        if (st_local("kept") != "" & Kk) st_matrix(st_local("kept"), bsrc)
    }

    hascons = (st_local("cons") != "")
    n = T + hascons
    pos = (cols(kept) ? sent[kept] : J(1, 0, .))
    if (hascons) {
        pos = pos, n
        bsrc = bsrc, _ctools_ebv_scalar(st_local("cons"))
    }
    m = cols(pos)
    b = J(1, n, 0)
    V = J(n, n, 0)
    if (m) {
        Vsrc = st_matrix(st_local("vmatrix"))
        if (cols(bsrc) != m | rows(Vsrc) < m | cols(Vsrc) < m) {
            errprintf("ctools: plugin results do not match the %g estimated coefficients\n", m)
            exit(498)
        }
        b[pos] = bsrc
        V[pos, pos] = Vsrc[|1, 1 \ m, m|]
    }
    st_matrix(bname, b)
    st_matrix(Vname, V)

    for (i = 1; i <= Ks; i++) {
        if (coll[i]) names[sent[i]] = "o." + names[sent[i]]
    }
    if (hascons) names = names, "_cons"
    st_local("outnames", invtokens(names))
}
end
