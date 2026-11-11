* bootstrap_reference.do
* Reference draws for the bootstrap work in docs/development/BOOTSTRAP_PLAN.md.
*
* bsref_draws, reps(#) generate(stub) [strata(varlist) cluster(varlist)]
*
* Creates double variables stub1..stub# holding the frequency weights that #
* successive -bsample- calls with the same strata() and cluster() options
* would draw from the current random-number state, as -bootstrap- draws its
* replicates. The data keep their order. Random numbers are taken from
* Stata's generator in the order -bsample- takes them, so the generator ends
* in the state those calls leave it in.
*
* The emulation follows bsample's documented design: units (observations, or
* clusters numbered in sorted order within strata) are drawn with replacement
* within each stratum by one uniform per unit, and a second uniform per unit
* decides which row of a repeated draw carries its count. Every observation
* is in the sample; strata and cluster variables may not be missing.

capture program drop bsref_draws
program define bsref_draws
    version 14.1
    syntax, Reps(integer) GENerate(name) [STRata(varlist) CLuster(varlist)]
    if `reps' < 1 {
        di as error "reps() must be positive"
        exit 198
    }
    forvalues r = 1/`reps' {
        confirm new variable `generate'`r'
    }
    foreach v of local strata {
        quietly count if missing(`v')
        if r(N) {
            di as error "strata variable `v' has missing values"
            exit 416
        }
    }
    foreach v of local cluster {
        quietly count if missing(`v')
        if r(N) {
            di as error "cluster variable `v' has missing values"
            exit 416
        }
    }
    tempvar order sidx cidx
    quietly generate double `order' = _n
    * bsample sorts by strata and cluster (stable) when either is given and
    * otherwise draws in the current order
    if "`strata'`cluster'" != "" sort `strata' `cluster', stable
    * Number strata and clusters without reordering rows: the order of rows
    * within a stratum decides which row each observation draw refers to.
    if "`strata'" != "" {
        quietly by `strata': generate double `sidx' = (_n == 1)
        quietly replace `sidx' = sum(`sidx')
    }
    else quietly generate double `sidx' = 1
    * Units drawn per replicate: rows, or clusters summed over strata
    local units = _N
    if "`cluster'" != "" {
        quietly by `strata' `cluster': generate double `cidx' = (_n == 1)
        quietly count if `cidx'
        local units = r(N)
        if "`strata'" != "" quietly by `strata': replace `cidx' = sum(`cidx')
        else quietly replace `cidx' = sum(`cidx')
    }
    * bsample takes one uniform per unit for the draws, then one per unit for
    * the tie-breaks, in data order; the Mata routines only consume them
    forvalues r = 1/`reps' {
        quietly generate double `generate'`r' = .
        mata: _bsref_u_draw = runiform(`units', 1)
        mata: _bsref_u_tie = runiform(`units', 1)
        if "`cluster'" == "" {
            mata: st_store(., "`generate'`r'", ///
                _bsref_obs(st_data(., "`sidx'"), _bsref_u_draw, _bsref_u_tie))
        }
        else {
            mata: st_store(., "`generate'`r'", _bsref_clust(st_data(., "`sidx'"), ///
                st_data(., "`cidx'"), _bsref_u_draw, _bsref_u_tie))
        }
    }
    mata: mata drop _bsref_u_draw _bsref_u_tie
    sort `order'
end

mata:
// Observation resampling within strata. s[i] is the 1-based stratum of row
// i, with rows in bsample's order; u_draw and u_tie hold one uniform per row
// each, taken in that order. Returns each row's frequency weight.
real colvector _bsref_obs(real colvector s, real colvector u_draw,
    real colvector u_tie)
{
    real scalar N, i, k, run, prev, last
    real colvector ns, r, w, p, m

    N = rows(s)
    ns = J(max(s), 1, 0)
    for (i = 1; i <= N; i++) ns[s[i]] = ns[s[i]] + 1
    // the first uniform of a row gives its drawn index within the stratum;
    // the second breaks ties among repeated draws
    r = J(N, 1, .)
    for (i = 1; i <= N; i++) r[i] = trunc(u_draw[i] * ns[s[i]] + 1)
    w = u_tie
    p = order((s, r, w), (1, 2, 3))
    m = J(N, 1, 0)
    prev = 0
    for (k = 1; k <= N; k++) {
        run = 1
        if (k > 1) {
            if (s[p[k]] == s[p[k-1]] & r[p[k]] == r[p[k-1]]) run = prev + 1
        }
        prev = run
        last = (k == N)
        if (!last) last = (s[p[k+1]] != s[p[k]] | r[p[k+1]] != r[p[k]])
        // the last row of each run of equal draws carries the run's count
        if (last) m[p[k]] = run
    }
    return(m)
}

// Cluster resampling within strata. s[i] is the stratum and c[i] the 1-based
// cluster within the stratum of row i, with rows sorted by stratum and
// cluster; u_draw and u_tie hold one uniform per cluster each, stratum by
// stratum. Returns each row's frequency weight (that of its cluster).
real colvector _bsref_clust(real colvector s, real colvector c,
    real colvector u_draw, real colvector u_tie)
{
    real scalar N, S, i, j, k, t, run, prev, nc, first, last
    real colvector firstrow, ncl, r, w, rs, p, W, cm, m

    N = rows(s)
    S = max(s)
    firstrow = J(S, 1, .)
    ncl = J(S, 1, 0)
    for (i = 1; i <= N; i++) {
        if (firstrow[s[i]] == .) firstrow[s[i]] = i
        if (c[i] > ncl[s[i]]) ncl[s[i]] = c[i]
    }
    // bsample draws on the first nc rows of each stratum, stratum by
    // stratum: the draws for every stratum first, then the tie-breaks
    r = J(N, 1, .)
    w = J(N, 1, .)
    t = 0
    for (j = 1; j <= S; j++) {
        for (i = firstrow[j]; i < firstrow[j] + ncl[j]; i++) {
            t++
            r[i] = trunc(u_draw[t] * ncl[j] + 1)
            w[i] = u_tie[t]
        }
    }
    m = J(N, 1, 0)
    for (j = 1; j <= S; j++) {
        nc = ncl[j]
        first = firstrow[j]
        rs = r[|first \ first + nc - 1|]
        p = order((rs, w[|first \ first + nc - 1|], (1::nc)), (1, 2, 3))
        rs = rs[p]
        W = J(nc, 1, 0)
        prev = 0
        for (k = 1; k <= nc; k++) {
            run = 1
            if (k > 1) {
                if (rs[k] == rs[k-1]) run = prev + 1
            }
            prev = run
            last = (k == nc)
            if (!last) last = (rs[k+1] != rs[k])
            if (last) W[k] = run
        }
        // bsample maps the counts back through the sort permutation: the
        // cluster numbered q receives the count at position p[q]
        cm = W[p]
        for (i = 1; i <= N; i++) {
            if (s[i] == j) m[i] = cm[c[i]]
        }
    }
    return(m)
}

// Wild cluster bootstrap, restricted (WCR) to H0: b_k = 0, computed
// replicate by replicate: OLS of y on X (unweighted, X holding the
// constant), the restricted fit without column k, y* = fit + resid :* v,
// and the cluster-robust t of b_k without small-sample factors. cl: the
// cluster (1 .. C) of each row, rows sorted by cl; V: C x B weights.
// Returns the observed t, then t*_1 .. t*_B.
real colvector _bsref_wild_t(real matrix X, real colvector y, real colvector cl,
    real scalar k, real matrix V)
{
    return(_bsref_wild_t2(X, y, cl, cl, k, V))
}

// The same with bootstrap clusters bc (1 .. C, rows of V) that may differ
// from the vce() clusters cl
real colvector _bsref_wild_t2(real matrix X, real colvector y, real colvector cl,
    real colvector bc, real scalar k, real matrix V)
{
    real matrix XX, Xr, info
    real colvector b, br, yr, ut, ys, bs, out
    real scalar r

    XX = invsym(cross(X, X))
    b = XX * cross(X, y)
    Xr = X[., select(1..cols(X), (1..cols(X)) :!= k)]
    br = invsym(cross(Xr, Xr)) * cross(Xr, y)
    yr = Xr * br
    ut = y - yr
    info = panelsetup(cl, 1)
    out = J(cols(V) + 1, 1, .)
    out[1] = b[k] / sqrt(_bsref_wild_vkk(X, XX, y - X * b, info, k))
    for (r = 1; r <= cols(V); r++) {
        ys = yr + ut :* V[bc, r]
        bs = XX * cross(X, ys)
        out[r + 1] = bs[k] / sqrt(_bsref_wild_vkk(X, XX, ys - X * bs, info, k))
    }
    return(out)
}

real scalar _bsref_wild_vkk(real matrix X, real matrix XX, real colvector u,
    real matrix info, real scalar k)
{
    real matrix meat
    real colvector s
    real scalar c

    meat = J(cols(X), cols(X), 0)
    for (c = 1; c <= rows(info); c++) {
        s = cross(panelsubmatrix(X, c, info), panelsubmatrix(u, c, info))
        meat = meat + s * s'
    }
    return((XX * meat * XX)[k, k])
}

// boottest's enumeration: all 2^C Rademacher vectors, replicate j (0-based)
// with v_c = +1 when bit C - c of j is set (c = 1 .. C)
real matrix _bsref_wild_enum(real scalar C)
{
    real matrix V
    real scalar j, c

    V = J(C, 2^C, -1)
    for (j = 0; j < 2^C; j++) {
        for (c = 1; c <= C; c++) {
            if (mod(floor(j / 2^(C - c)), 2)) V[c, j + 1] = 1
        }
    }
    return(V)
}

// n values from 0-based position q0 of those _ctools_hdfe_boot packs into
// stub1, stub2, ... (value q in variable floor(q / _N) + 1, observation
// mod(q, _N) + 1), read one value at a time
real colvector _bsref_packed(string scalar stub, real scalar q0, real scalar n)
{
    real colvector x
    real scalar i, q

    x = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        q = q0 + i - 1
        x[i] = _st_data(mod(q, st_nobs()) + 1, st_varindex(stub + strofreal(floor(q / st_nobs()) + 1)))
    }
    return(x)
}
// The WRE bootstrap (Davidson and MacKinnon 2010, as boottest computes it):
// the t statistics of H0: b_k = b0 for the k-class fit of y on Z = [Y2, X1]
// with instruments X = [X1, X2], the observed one first, then one per
// column of V. The restricted fit is LIML; the efficient first stage
// regresses Y2 on X and its residual u; replicate r uses Y2* = X Pi + U2 v,
// y* = Z* b~ + u v with v = V[bc, r]. est: 0 2SLS, 1 LIML, 2 Fuller
// (lambda - kc), 3 k-class (kc); w weights; cl the vce() clusters (data
// sorted by cl; 1 .. N for robust); bc the bootstrap clusters (rows of V).
// Dense and direct: absorbed effects enter as dummies in X1.
real colvector _bsref_wre_t(real colvector y, real matrix Y2, real matrix X1,
    real matrix X2, real colvector w, real colvector cl, real colvector bc,
    real scalar k, real scalar b0, real scalar est, real scalar kc, real matrix V)
{
    real matrix sw, Z, X, XXi, Zr, Yd, X1r, A, Bm, XU, P, Pi, U2, Y2s, Zs, info
    real colvector yy, yr, br, bt, u, out, v, ys
    real rowvector sel, endo, exo
    real scalar K, K2, kap, r

    sw = sqrt(w)
    yy = y :* sw
    Z = (Y2 :* sw, X1 :* sw)
    X = (X1 :* sw, X2 :* sw)
    K2 = cols(Y2)
    K = cols(Z)
    XXi = invsym(cross(X, X))
    sel = select(1..K, (1..K) :!= k)
    endo = select(sel, sel :<= K2)
    exo = select(sel, sel :> K2)
    Zr = Z[., sel]
    yr = yy - b0 * Z[., k]
    kap = 1
    if (cols(endo)) {
        Yd = (yr, Z[., endo])
        X1r = cols(exo) ? Z[., exo] : J(rows(Z), 0, .)
        A = cross(Yd, Yd)
        if (cols(exo)) A = A - cross(Yd, X1r) * invsym(cross(X1r, X1r)) * cross(X1r, Yd)
        Bm = cross(Yd, Yd) - cross(Yd, X) * XXi * cross(X, Yd)
        kap = _bsref_wre_lambda(A, Bm)
    }
    br = _bsref_kclass(Zr, yr, X, XXi, kap)
    bt = J(K, 1, 0)
    bt[sel] = br
    bt[k] = b0
    u = yy - Z * bt
    XU = (X, u)
    P = invsym(cross(XU, XU)) * cross(XU, Z[., 1..K2])
    Pi = P[1..cols(X), .]
    U2 = Z[., 1..K2] - X * Pi
    info = panelsetup(cl, 1)
    out = J(cols(V) + 1, 1, .)
    for (r = 0; r <= cols(V); r++) {
        v = r ? V[bc, r] : J(rows(Z), 1, 1)
        Y2s = X * Pi + U2 :* v
        Zs = (Y2s, Z[., (K2 + 1)..K])
        ys = Zs * bt + u :* v
        out[r + 1] = _bsref_wre_tk(Zs, ys, X, XXi, info, k, b0, est, kc, cols(Z) - K2)
    }
    return(out)
}

// LIML's lambda: the smallest root of |A - lambda Bm| = 0
real scalar _bsref_wre_lambda(real matrix A, real matrix Bm)
{
    real matrix L, Li

    L = cholesky(Bm)
    Li = luinv(L)
    return(min(symeigenvalues(makesymmetric(Li * A * Li'))))
}

real colvector _bsref_kclass(real matrix Z, real colvector y, real matrix X,
    real matrix XXi, real scalar kap)
{
    real matrix XZ
    real colvector Xy

    XZ = cross(X, Z)
    Xy = cross(X, y)
    return(invsym((1 - kap) * cross(Z, Z) + kap * XZ' * XXi * XZ) * ///
        ((1 - kap) * cross(Z, y) + kap * XZ' * XXi * Xy))
}

// The estimator's t statistic for b_k = b0 with the cluster-robust
// sandwich (P_X Z in the middle, as ivreg2)
real scalar _bsref_wre_tk(real matrix Zs, real colvector ys, real matrix X,
    real matrix XXi, real matrix info, real scalar k, real scalar b0,
    real scalar est, real scalar kc, real scalar K1)
{
    real matrix XZ, H, Hi, PZ, meat, Yd, X1, A, Bm
    real colvector b, e, s
    real scalar kap, K2, c

    XZ = cross(X, Zs)
    K2 = cols(Zs) - K1
    kap = 1
    if (est == 3) kap = kc
    else if (est == 1 | est == 2) {
        if (cols(X) == cols(Zs)) kap = 1
        else {
            Yd = (ys, Zs[., 1..K2])
            X1 = Zs[., (K2 + 1)..cols(Zs)]
            A = cross(Yd, Yd)
            if (K1) A = A - cross(Yd, X1) * invsym(cross(X1, X1)) * cross(X1, Yd)
            Bm = cross(Yd, Yd) - cross(Yd, X) * XXi * cross(X, Yd)
            kap = _bsref_wre_lambda(A, Bm)
        }
        if (est == 2) kap = kap - kc
    }
    H = (1 - kap) * cross(Zs, Zs) + kap * XZ' * XXi * XZ
    Hi = invsym(H)
    b = Hi * ((1 - kap) * cross(Zs, ys) + kap * XZ' * XXi * cross(X, ys))
    e = ys - Zs * b
    PZ = X * (XXi * XZ)
    meat = J(cols(Zs), cols(Zs), 0)
    for (c = 1; c <= rows(info); c++) {
        s = cross(panelsubmatrix(PZ, c, info), panelsubmatrix(e, c, info))
        meat = meat + s * s'
    }
    return((b[k] - b0) / sqrt((Hi * meat * Hi)[k, k]))
}
end
