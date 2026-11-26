*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! C-accelerated joinby: all pairwise combinations within groups (ctools)
*! original: Stata's joinby (joinby.ado 2.1.0, December 22, 2022)
*! original: see help cjoinby for details

* Mata helpers: what the plugin cannot reach through Stata's plugin
* interface (variable creation, labels, formats, characteristics, notes and
* strL values). The observations themselves are joined in the plugin.
capture mata: mata drop _cjoinby_meta_collect()
capture mata: mata drop _cjoinby_meta_apply()
capture mata: mata drop _cjoinby_notes_append()
capture mata: mata drop _cjoinby_addvars()
capture mata: mata drop _cjoinby_undo_record()
capture mata: mata drop _cjoinby_undo_restore()
capture mata: mata drop _cjoinby_aliases()
capture mata: mata drop _cjoinby_strl_hold()
capture mata: mata drop _cjoinby_strl_fill()
capture mata: mata drop _cjoinby_strl_update()
capture mata: mata drop _cjoinby_drop_externals()
capture mata: mata drop _cjoinby_wider()
capture mata: mata drop _cjoinby_first()
capture mata: mata drop _cjoinby_types()
capture mata: mata drop _cjoinby_classify()
capture mata: mata drop _cjoinby_promotions()
capture mata: mata drop _cjoinby_plugin_vars()
capture mata: mata drop _cjoinby_positions()
capture mata: mata drop _cjoinby_strw()
mata:
/* The storage type that holds both a master and a using variable of the same
   generic type, as merge chooses it */
string scalar _cjoinby_wider(string scalar mt, string scalar ut)
{
    string rowvector ranked
    real scalar mr, ur, rank, j

    if (substr(mt, 1, 3) == "str") {
        if (mt == "strL" | ut == "strL") return("strL")
        return("str" + strofreal(max((strtoreal(substr(mt, 4, .)), strtoreal(substr(ut, 4, .))))))
    }
    ranked = ("byte", "int", "long", "float", "double")
    mr = ur = 0
    for (j = 1; j <= 5; j++) {
        if (ranked[j] == mt) mr = j
        if (ranked[j] == ut) ur = j
    }
    rank = max((mr, ur))
    /* float cannot hold every long */
    if ((mr == 3 & ur == 4) | (mr == 4 & ur == 3)) rank = 5
    return(ranked[rank])
}

/* The first n names of v as a macro list */
string scalar _cjoinby_first(string rowvector v, real scalar n)
{
    return(n ? invtokens(v[|1 \ n|]) : "")
}

/* The storage types of the variables in memory, in local mastertypes */
void _cjoinby_types()
{
    real scalar j
    string rowvector types

    types = J(1, st_nvar(), "")
    for (j = 1; j <= st_nvar(); j++) types[j] = st_vartype(j)
    st_local("mastertypes", invtokens(types))
}

/* Run in the using data. Sorts its variables, other than the keys, into
   locals of the caller:
     newvars, newtypes   the variables the master does not have
     newstrl, paynew     the strL ones among them, and the others (which the
                         plugin writes)
     shared, sharedtargets  the variables both have, of one generic type, and
                         the storage type that holds both
     sharedstrl          the shared ones held as strL on either side
     payshared           under update, the other shared ones (the plugin
                         applies update to them)
   A shared variable that is a string in one dataset and a number in the
   other, or a frame alias in the master, is in none: merge leaves it alone. */
void _cjoinby_classify(string rowvector keys, string rowvector mastervars,
    string rowvector mastertypes, string rowvector aliasvars, real scalar update)
{
    transmorphic M, K, A
    string rowvector newvars, newtypes, newstrl, paynew, shared, targets, sharedstrl, payshared
    string scalar name, ut, mt, target
    real scalar j, k, nnew, nstrl, npay, nshared, nsstrl, nspay

    M = asarray_create()
    for (j = 1; j <= cols(mastervars); j++) asarray(M, mastervars[j], j)
    K = asarray_create()
    for (j = 1; j <= cols(keys); j++) asarray(K, keys[j], 1)
    A = asarray_create()
    for (j = 1; j <= cols(aliasvars); j++) asarray(A, aliasvars[j], 1)
    k = st_nvar()
    newvars = newtypes = newstrl = paynew = J(1, k, "")
    shared = targets = sharedstrl = payshared = J(1, k, "")
    nnew = nstrl = npay = nshared = nsstrl = nspay = 0
    for (j = 1; j <= k; j++) {
        name = st_varname(j)
        if (asarray_contains(K, name)) continue
        ut = st_vartype(j)
        if (!asarray_contains(M, name)) {
            newvars[++nnew] = name
            newtypes[nnew] = ut
            if (ut == "strL") newstrl[++nstrl] = name
            else paynew[++npay] = name
            continue
        }
        if (asarray_contains(A, name)) continue
        mt = mastertypes[asarray(M, name)]
        if ((substr(mt, 1, 3) == "str") != (substr(ut, 1, 3) == "str")) continue
        target = _cjoinby_wider(mt, ut)
        shared[++nshared] = name
        targets[nshared] = target
        if (target == "strL") sharedstrl[++nsstrl] = name
        else if (update) payshared[++nspay] = name
    }
    st_local("newvars", _cjoinby_first(newvars, nnew))
    st_local("newtypes", _cjoinby_first(newtypes, nnew))
    st_local("newstrl", _cjoinby_first(newstrl, nstrl))
    st_local("paynew", _cjoinby_first(paynew, npay))
    st_local("shared", _cjoinby_first(shared, nshared))
    st_local("sharedtargets", _cjoinby_first(targets, nshared))
    st_local("sharedstrl", _cjoinby_first(sharedstrl, nsstrl))
    st_local("payshared", _cjoinby_first(payshared, nspay))
}

/* The variables of vars whose storage type differs from its target, in
   locals promote (names) and promoteto (types) */
void _cjoinby_promotions(string rowvector vars, string rowvector targets)
{
    string rowvector names, types
    real scalar j, n

    names = types = J(1, cols(vars), "")
    n = 0
    for (j = 1; j <= cols(vars); j++) {
        if (st_vartype(vars[j]) == targets[j]) continue
        names[++n] = vars[j]
        types[n] = targets[j]
    }
    st_local("promote", _cjoinby_first(names, n))
    st_local("promoteto", _cjoinby_first(types, n))
}

/* The master's variables as the plugin sees them now: in local mvars the
   ones it moves, in local mstrl the strL ones, which it cannot write. Frame
   aliases are in neither: their values follow the link variable. */
void _cjoinby_plugin_vars(string rowvector mastervars, string rowvector aliasvars)
{
    transmorphic A
    string rowvector mvars, mstrl
    real scalar j, nm, ns

    A = asarray_create()
    for (j = 1; j <= cols(aliasvars); j++) asarray(A, aliasvars[j], 1)
    mvars = mstrl = J(1, cols(mastervars), "")
    nm = ns = 0
    for (j = 1; j <= cols(mastervars); j++) {
        if (st_vartype(mastervars[j]) == "strL") mstrl[++ns] = mastervars[j]
        else if (!asarray_contains(A, mastervars[j])) mvars[++nm] = mastervars[j]
    }
    st_local("mvars", _cjoinby_first(mvars, nm))
    st_local("mstrl", _cjoinby_first(mstrl, ns))
}

/* The string widths of the plugin varlist, in local __ctools_strw, as
   _ctools_strw sets it (a comma-separated width per variable, 0 for numeric
   and strL variables): the plugin stages fixed-width strings by them. One
   pass here, where _ctools_strw's macro loop grows with the square of the
   number of variables. */
void _cjoinby_strw(string rowvector vars)
{
    string rowvector widths
    string scalar type
    real scalar j

    widths = J(1, cols(vars), "0")
    for (j = 1; j <= cols(vars); j++) {
        type = st_vartype(vars[j])
        if (substr(type, 1, 3) == "str" & type != "strL") widths[j] = substr(type, 4, .)
    }
    st_local("__ctools_strw", invtokens(widths, ","))
}

/* What the plugin's write call reads from locals: the position among mvars
   of each key (__cj_keypos) and of each shared payload variable, followed by
   a 0 for each of the nnew new payload variables (__cj_paymap) */
void _cjoinby_positions(string rowvector mvars, string rowvector keys,
    string rowvector payshared, real scalar nnew)
{
    transmorphic M
    real rowvector keypos, paymap
    real scalar j

    M = asarray_create()
    for (j = 1; j <= cols(mvars); j++) asarray(M, mvars[j], j)
    keypos = J(1, cols(keys), 0)
    for (j = 1; j <= cols(keys); j++) keypos[j] = asarray(M, keys[j])
    paymap = J(1, cols(payshared) + nnew, 0)
    for (j = 1; j <= cols(payshared); j++) paymap[j] = asarray(M, payshared[j])
    st_local("__cj_keypos", invtokens(strofreal(keypos)))
    st_local("__cj_paymap", invtokens(strofreal(paymap)))
}

/* Run in the using data. Keeps in external __cjoinby_M what joinby (through
   merge) passes on from the using data:
     "new"     name, format, variable label and value-label name of each new
               variable
     "common"  name, variable label and value-label name of each key and
               shared variable
     "chars"   owner, name and text of every characteristic (notes included)
               of the new, key and shared variables and of _dta
     "labels"  the value labels to define, with their values ("V:" + name)
               and texts ("T:" + name) */
void _cjoinby_meta_collect(string rowvector newvars, string rowvector common,
    string rowvector labels)
{
    transmorphic A
    string matrix N, K, C
    string rowvector owners
    string colvector names, text, kept
    real colvector values
    real scalar j, i
    pointer() scalar p

    A = asarray_create()
    N = J(cols(newvars), 4, "")
    for (j = 1; j <= cols(newvars); j++) {
        N[j, 1] = newvars[j]
        N[j, 2] = st_varformat(newvars[j])
        N[j, 3] = st_varlabel(newvars[j])
        N[j, 4] = st_varvaluelabel(newvars[j])
    }
    K = J(cols(common), 3, "")
    for (j = 1; j <= cols(common); j++) {
        K[j, 1] = common[j]
        K[j, 2] = st_varlabel(common[j])
        K[j, 3] = st_varvaluelabel(common[j])
    }
    owners = (newvars, common, "_dta")
    C = J(0, 3, "")
    for (j = 1; j <= cols(owners); j++) {
        names = st_dir("char", owners[j], "*")
        for (i = 1; i <= rows(names); i++) {
            C = C \ (owners[j], names[i], st_global(owners[j] + "[" + names[i] + "]"))
        }
    }
    kept = J(0, 1, "")
    for (j = 1; j <= cols(labels); j++) {
        if (!st_vlexists(labels[j])) continue
        st_vlload(labels[j], values, text)
        if (rows(values) == 0) continue
        kept = kept \ labels[j]
        asarray(A, "V:" + labels[j], values)
        asarray(A, "T:" + labels[j], text)
    }
    asarray(A, "new", N)
    asarray(A, "common", K)
    asarray(A, "chars", C)
    asarray(A, "labels", kept)
    rmexternal("__cjoinby_M")
    p = crexternal("__cjoinby_M")
    *p = A
}

/* Append the using notes of owner (a variable or _dta) to the master's, as
   merge does: a note whose text the master's own notes already hold is left
   out. notes are the using note texts, in order. */
void _cjoinby_notes_append(string scalar owner, string colvector notes)
{
    real scalar n0, k, j, h
    string colvector have

    n0 = strtoreal(st_global(owner + "[note0]"))
    if (n0 >= . | n0 < 0 | n0 != trunc(n0)) n0 = 0
    have = J(n0, 1, "")
    for (h = 1; h <= n0; h++) have[h] = st_global(owner + "[note" + strofreal(h) + "]")
    k = n0
    for (j = 1; j <= rows(notes); j++) {
        if (anyof(have, notes[j])) continue
        k++
        st_global(owner + "[note" + strofreal(k) + "]", notes[j])
    }
    if (k != n0) st_global(owner + "[note0]", strofreal(k))
}

/* Run in the master, once the new variables exist. Applies what
   _cjoinby_meta_collect kept: the value labels the master does not define;
   formats, labels and characteristics of the new variables; for the key and
   shared variables and for _dta, the using notes and the characteristics the
   master does not have; and, under update, the using variable label and
   value label of a key or shared variable that has none. */
void _cjoinby_meta_apply(real scalar update)
{
    transmorphic A
    string matrix N, K, C
    string colvector labels, notes
    string scalar owner, name
    real scalar j, i, n, isnew
    transmorphic isnewvar
    pointer() scalar p

    p = findexternal("__cjoinby_M")
    if (p == NULL) _error(3498, "cjoinby: the using metadata are gone")
    A = *p
    labels = asarray(A, "labels")
    for (j = 1; j <= rows(labels); j++) {
        if (st_vlexists(labels[j])) continue
        st_vlmodify(labels[j], asarray(A, "V:" + labels[j]), asarray(A, "T:" + labels[j]))
    }
    N = asarray(A, "new")
    isnewvar = asarray_create()
    for (j = 1; j <= rows(N); j++) {
        asarray(isnewvar, N[j, 1], 1)
        st_varformat(N[j, 1], N[j, 2])
        st_varlabel(N[j, 1], N[j, 3])
        if (N[j, 4] != "") st_varvaluelabel(N[j, 1], N[j, 4])
    }
    C = asarray(A, "chars")
    i = 1
    while (i <= rows(C)) {
        owner = C[i, 1]
        isnew = asarray_contains(isnewvar, owner)
        n = 0
        for (j = i; j <= rows(C); j++) {
            if (C[j, 1] != owner) break
            if (C[j, 2] == "note0") n = strtoreal(C[j, 3])
        }
        if (n >= . | n < 0 | n != trunc(n)) n = 0
        notes = J(n, 1, "")
        for (j = i; j <= rows(C); j++) {
            if (C[j, 1] != owner) break
            name = C[j, 2]
            if (isnew) {
                st_global(owner + "[" + name + "]", C[j, 3])
                continue
            }
            if (regexm(name, "^note[0-9]+$")) {
                if (name != "note0" & strtoreal(substr(name, 5, .)) <= n) {
                    notes[strtoreal(substr(name, 5, .))] = C[j, 3]
                }
                continue
            }
            if (st_global(owner + "[" + name + "]") == "") {
                st_global(owner + "[" + name + "]", C[j, 3])
            }
        }
        if (!isnew & n > 0) _cjoinby_notes_append(owner, notes)
        i = j
    }
    if (update) {
        K = asarray(A, "common")
        for (j = 1; j <= rows(K); j++) {
            if (K[j, 2] != "" & st_varlabel(K[j, 1]) == "") st_varlabel(K[j, 1], K[j, 2])
            if (K[j, 3] != "" & st_varvaluelabel(K[j, 1]) == "") {
                st_varvaluelabel(K[j, 1], K[j, 3])
            }
        }
    }
}

/* Add variables at the end of the dataset. The plugin writes every
   observation of the fixed-width ones, so they are not filled first; strL
   variables start empty. */
void _cjoinby_addvars(string rowvector names, string rowvector types)
{
    real scalar j

    for (j = 1; j <= cols(names); j++) {
        (void) st_addvar(types[j], names[j], types[j] != "strL")
    }
}

/* Record what a join may change in the key and shared variables: storage
   type, format, labels and characteristics (kept in external
   __cjoinby_undo for _cjoinby_undo_restore). */
void _cjoinby_undo_record(string rowvector vars)
{
    transmorphic A
    string matrix V, C
    string colvector names
    real scalar j, i
    pointer() scalar p

    V = J(cols(vars), 5, "")
    C = J(0, 3, "")
    for (j = 1; j <= cols(vars); j++) {
        V[j, 1] = vars[j]
        V[j, 2] = st_vartype(vars[j])
        V[j, 3] = st_varformat(vars[j])
        V[j, 4] = st_varlabel(vars[j])
        V[j, 5] = st_varvaluelabel(vars[j])
        names = st_dir("char", vars[j], "*")
        for (i = 1; i <= rows(names); i++) {
            C = C \ (vars[j], names[i], st_global(vars[j] + "[" + names[i] + "]"))
        }
    }
    A = asarray_create()
    asarray(A, "vars", V)
    asarray(A, "chars", C)
    rmexternal("__cjoinby_undo")
    p = crexternal("__cjoinby_undo")
    *p = A
}

/* Put back what _cjoinby_undo_record kept. The observations are the
   master's again, so every variable fits its original storage type. */
void _cjoinby_undo_restore()
{
    transmorphic A
    string matrix V, C
    string colvector now, old
    real scalar j, i
    pointer() scalar p

    p = findexternal("__cjoinby_undo")
    if (p == NULL) return
    A = *p
    V = asarray(A, "vars")
    C = asarray(A, "chars")
    for (j = 1; j <= rows(V); j++) {
        if (_st_varindex(V[j, 1]) >= .) continue
        if (st_vartype(V[j, 1]) != V[j, 2]) {
            stata("capture noisily recast " + V[j, 2] + " " + V[j, 1])
        }
        st_varformat(V[j, 1], V[j, 3])
        st_varlabel(V[j, 1], V[j, 4])
        if (st_varvaluelabel(V[j, 1]) != V[j, 5]) st_varvaluelabel(V[j, 1], V[j, 5])
        now = st_dir("char", V[j, 1], "*")
        old = select(C[., 2], C[., 1] :== V[j, 1])
        for (i = 1; i <= rows(now); i++) {
            if (!anyof(old, now[i])) st_global(V[j, 1] + "[" + now[i] + "]", "")
        }
    }
    for (i = 1; i <= rows(C); i++) {
        if (_st_varindex(C[i, 1]) >= .) continue
        if (st_global(C[i, 1] + "[" + C[i, 2] + "]") != C[i, 3]) {
            st_global(C[i, 1] + "[" + C[i, 2] + "]", C[i, 3])
        }
    }
    rmexternal("__cjoinby_undo")
}

/* The frame alias variables (Stata 18) of the data in memory, in local
   aliasvars: values of an alias follow its link variable, and the plugin's
   stores into it would have no effect. */
void _cjoinby_aliases()
{
    real scalar j
    string scalar found

    found = ""
    for (j = 1; j <= st_nvar(); j++) {
        if (st_isalias(st_varname(j))) found = found + " " + st_varname(j)
    }
    st_local("aliasvars", strtrim(found))
}

/* Run in the using data: keep the values of the named variables, which the
   plugin cannot write (strL), in external __cjoinby_U. */
void _cjoinby_strl_hold(string rowvector names)
{
    pointer() scalar p

    rmexternal("__cjoinby_U")
    p = crexternal("__cjoinby_U")
    *p = st_sdata(., names)
}

/* Store column col of __cjoinby_U in observations 1..n of strL variable
   name, by the using observation uobs holds; observations without one keep
   "". In blocks of rows, to bound the strings held at once. */
void _cjoinby_strl_fill(string scalar name, real scalar col, string scalar uobs,
    real scalar n)
{
    pointer() scalar p
    real scalar a, b
    real colvector idx, sel
    string colvector values

    p = findexternal("__cjoinby_U")
    if (p == NULL) _error(3498, "cjoinby: the using strL values are gone")
    for (a = 1; a <= n; a = b + 1) {
        b = min((a + 99999, n))
        idx = st_data((a, b), uobs)
        sel = selectindex(idx :< .)
        if (rows(sel) == 0) continue
        values = J(rows(idx), 1, "")
        values[sel] = (*p)[idx[sel], col]
        st_sstore((a, b), name, values)
    }
}

/* update for a shared variable held as strL: in observations 1..n with a
   using observation, an empty value takes the using value (column col of
   __cjoinby_U); with replace, so does any value when the using value is not
   empty. */
void _cjoinby_strl_update(string scalar name, real scalar col, string scalar uobs,
    real scalar replace, real scalar n)
{
    pointer() scalar p
    real scalar a, b
    real colvector idx, sel, take
    string colvector values, incoming

    p = findexternal("__cjoinby_U")
    if (p == NULL) _error(3498, "cjoinby: the using strL values are gone")
    for (a = 1; a <= n; a = b + 1) {
        b = min((a + 99999, n))
        idx = st_data((a, b), uobs)
        sel = selectindex(idx :< .)
        if (rows(sel) == 0) continue
        values = st_sdata((a, b), name)
        incoming = (*p)[idx[sel], col]
        take = values[sel] :== ""
        if (replace) take = take :| (incoming :!= "")
        if (!any(take)) continue
        values[select(sel, take)] = select(incoming, take)
        st_sstore((a, b), name, values)
    }
}

void _cjoinby_drop_externals()
{
    rmexternal("__cjoinby_M")
    rmexternal("__cjoinby_U")
    rmexternal("__cjoinby_undo")
}
end

program define cjoinby, rclass
    version 14.1
    syntax [varlist(default=none)] using/ [, UNMatched(string) noLabel UPDATE REPLACE ///
        _merge(name) Verbose THReads(integer 0) SNAPshot FAST]
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
    * A failed or interrupted join leaves the master as it was without a
    * snapshot of the dataset (Stata 16 or later): the using data are read in
    * a frame of their own; the master only gains variables, value labels,
    * notes and characteristics and observations at the end; _cjoinby_impl
    * records what else it changes (storage types, labels and
    * characteristics of the key and shared variables); and the plugin keeps
    * the master observations it moves, writing them back after a failure.
    * Option snapshot preserves the whole dataset instead (on disk above
    * c(max_preservemem)), as does every join before Stata 16, where the
    * using data can only be read in place of the master. Option fast keeps
    * no copy and undoes nothing.
    * (validation hook: $CTOOLS_cjoinby_noframes runs the path of Stata 14
    * and 15 on a later Stata). Without frames _cjoinby_impl preserves the
    * master itself, around reading the using data in its place; Stata then
    * restores it when _cjoinby_impl ends in an error.
    local noframes = c(stata_version) < 16 | "$CTOOLS_cjoinby_noframes" == "1"
    local mode `snapshot'`fast'
    if `noframes' local mode snapshot
    local undo = "`mode'" == ""
    local presnap = "`mode'" == "snapshot" & !`noframes'
    local caller_frame
    if c(stata_version) >= 16 local caller_frame `c(frame)'
    * Stata keeps e(sample) in a hidden column that joinby moves with each
    * master observation (0 on using-only ones) but the plugin cannot reach.
    * When it is set, hold the estimates with e(sample) copied into a master
    * variable, which the join carries, and unhold them afterwards.
    local esample : e(functions)
    local esample : list posof "sample" in esample
    if `esample' {
        tempname eshold
        tempvar esvar
    }
    * Break is recognized only inside _cjoinby_impl.
    nobreak {
        macro drop CTOOLS_cjoinby_plugin CTOOLS_cjoinby_mvars CTOOLS_cjoinby_esvar
        if `esample' {
            _estimates hold `eshold', copy varname(`esvar')
            global CTOOLS_cjoinby_esvar `esvar'
        }
        if `presnap' preserve
        if `undo' _ctools_undo record CTOOLS_u_cjoinby_
        capture noisily break _cjoinby_impl `macval(0)'
        local rc = _rc
        if "`caller_frame'" != "" capture frame change `caller_frame'
        if `rc' {
            if "$CTOOLS_cjoinby_plugin" != "" {
                capture program ctools_plugin, plugin using("$CTOOLS_cjoinby_plugin")
                if `undo' & "$CTOOLS_cjoinby_mvars" != "" {
                    * Observations first, while the master variables have
                    * the storage types the plugin copied them with
                    if _N < $CTOOLS_u_cjoinby_N quietly set obs $CTOOLS_u_cjoinby_N
                    capture noisily plugin call ctools_plugin $CTOOLS_cjoinby_mvars, "cjoinby rollback"
                    if _rc di as error "cjoinby: the master observations could not be restored"
                }
                capture plugin call ctools_plugin, "cjoinby clear"
            }
            if `undo' {
                * storage types, labels and characteristics go back once the
                * added observations are gone; putting them back marks the
                * data changed again
                local changed0 = 0${CTOOLS_u_cjoinby_changed}
                local sorted0 ${CTOOLS_u_cjoinby_sortedby}
                _ctools_undo rollback CTOOLS_u_cjoinby_
                mata: _cjoinby_undo_restore()
                if "`sorted0'" != "" & "`: sortedby'" != "`sorted0'" sort `sorted0', stable
                mata: st_updata(`changed0')
            }
            if `presnap' restore
            macro drop CTOOLS_cjoinby_plugin CTOOLS_cjoinby_mvars CTOOLS_cjoinby_esvar
            mata: _cjoinby_drop_externals()
            if `esample' _cjoinby_esample `eshold' `esvar'
            capture _ctools_cleanup cjoinby
            exit `rc'
        }
        return add
        * Free what the plugin kept for a rollback
        if "$CTOOLS_cjoinby_plugin" != "" {
            capture program ctools_plugin, plugin using("$CTOOLS_cjoinby_plugin")
            capture plugin call ctools_plugin, "cjoinby clear"
        }
        if `undo' _ctools_undo clear CTOOLS_u_cjoinby_
        if `presnap' restore, not
        macro drop CTOOLS_cjoinby_plugin CTOOLS_cjoinby_mvars CTOOLS_cjoinby_esvar
        mata: _cjoinby_drop_externals()
        if `esample' _cjoinby_esample `eshold' `esvar'
        capture _ctools_cleanup cjoinby
    }
end

* Unhold the estimates cjoinby held, with e(sample) taken from the master
* variable the join carried: 0 on using-only observations, as joinby leaves
* them (_cjoinby_impl left it on the first pairing of each master
* observation only). Without that variable (option fast after a failure) the
* estimates return with an empty e(sample).
program define _cjoinby_esample
    args hold var
    capture confirm variable `var', exact
    if !_rc {
        quietly count if missing(`var')
        if r(N) quietly replace `var' = 0 if missing(`var')
    }
    _estimates unhold `hold'
end

program define _cjoinby_impl, rclass
    version 14.1
    syntax [varlist(default=none)] using/ [, UNMatched(string) noLabel UPDATE REPLACE ///
        _merge(name) Verbose THReads(integer 0) SNAPshot FAST]

    local timing = "`verbose'" != ""
    if `timing' quietly _ctools_clock t_begin

    * =========================================================================
    * Options
    * =========================================================================
    if `threads' < 0 {
        di as error "threads() must be an integer from 1 through 256"
        exit 198
    }
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    _cjoinby_unmatched `"`unmatched'"'
    if "`replace'" != "" & "`update'" == "" {
        di as error "option {bf:replace} may not be specified without {bf:update}"
        exit 198
    }
    local mode `snapshot'`fast'
    if c(stata_version) < 16 | "$CTOOLS_cjoinby_noframes" == "1" local mode snapshot
    local undo = "`mode'" == ""
    local nolabel `label'

    * joinby creates _merge only with unmatched() or _merge()
    * (a default _merge the master already has is r(110); a name given in
    * _merge() is refused by joinby's merge, r(198), after the using data
    * were checked for it, r(110))
    local mergevar `_merge'
    if "`mergevar'" == "" & "`unmatch'" != "none" local mergevar _merge
    local mergeclash 0
    if "`mergevar'" != "" {
        capture confirm new variable `mergevar'
        if _rc {
            if "`_merge'" == "" {
                di as error "variable `mergevar' already defined"
                exit 110
            }
            local mergeclash 1
        }
    }

    * =========================================================================
    * The master, before anything changes
    * =========================================================================
    if c(k) == 0 {
        di as error "no variables defined"
        exit 111
    }
    local N0 = _N
    unab mastervars : _all
    local sortedby : sortedby
    mata: _cjoinby_types()
    local aliasvars
    if c(stata_version) >= 18 mata: _cjoinby_aliases()
    quietly label dir
    local masterlabels `r(names)'

    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')

    * =========================================================================
    * The using data: in a frame of their own (Stata 16 or later), or in
    * place of the preserved master. The plugin stages and indexes them; only
    * their metadata return to the master.
    * =========================================================================
    local useframes = c(stata_version) >= 16 & "$CTOOLS_cjoinby_noframes" != "1"
    if `useframes' {
        local masterframe `c(frame)'
        tempname uframe
        frame create `uframe'
        frame change `uframe'
    }
    else preserve
    capture noisily {
        quietly use `"`macval(using)'"', clear
        if `timing' quietly _ctools_lap t_use t_begin `t_begin'
        if _N > 2147483647 exit 920
        local usingvars
        if c(k) > 0 unab usingvars : _all

        * Key variables: the varlist, or the variables the datasets share
        if "`varlist'" == "" {
            local keys : list mastervars & usingvars
            if "`keys'" == "" {
                di as error "no variables in common to master and using data"
                di as error "thus a join is not possible"
                exit 198
            }
            di as text "join on common variables: " as result "`keys'"
        }
        else local keys : list uniq varlist
        local keytargets
        foreach v of local keys {
            if !`: list v in usingvars' {
                di as error "variable `v' not found"
                exit 111
            }
            local mt : word `: list posof "`v'" in mastervars' of `mastertypes'
            local ut : type `v'
            if `: list v in aliasvars' {
                di as error "variable `v' is alias in master data but `ut' in using data"
                di as error "{p 4 4 2}Key variables (on which observations are matched) may be"
                di as error "type alias in {bf:joinby}, but {bf:cjoinby} matches stored values:"
                di as error "use {bf:frunalias} to recast the key variables in the master data.{p_end}"
                exit 106
            }
            if "`mt'" == "strL" | "`ut'" == "strL" {
                if "`mt'" == "strL" di as error "key variable `v' is strL"
                else di as error "key variable `v' is strL in using data."
                di as error "{p 4 4 2}The key variables (on which observations are matched)"
                di as error "can be str#, but they cannot be strLs.{p_end}"
                exit 106
            }
            if (substr("`mt'", 1, 3) == "str") != (substr("`ut'", 1, 3) == "str") {
                di as error "key variable `v' is `mt' in master but `ut' in using data"
                di as error "{p 4 4 2}Each key variable (on which observations are matched)"
                di as error "must be of the same generic type in the master and using"
                di as error "datasets.  Same generic type means both numeric or both string.{p_end}"
                exit 106
            }
            mata: st_local("target", _cjoinby_wider("`mt'", "`ut'"))
            local keytargets `keytargets' `target'
        }
        if "`mergevar'" != "" {
            if `: list mergevar in usingvars' {
                di as error "variable `mergevar' already defined"
                exit 110
            }
            if `mergeclash' {
                di as error "`mergevar' invalid name"
                exit 198
            }
        }

        * The other using variables. A variable the master also has is
        * shared: the master's values stand, unless update fills them in. A
        * shared variable that is a string in one dataset and a number in
        * the other, or a frame alias in the master, is left alone, as merge
        * leaves it.
        mata: _cjoinby_classify(tokens(st_local("keys")), tokens(st_local("mastervars")), ///
            tokens(st_local("mastertypes")), tokens(st_local("aliasvars")), `="`update'" != ""')

        * Value labels: every one the using file defines, or with nolabel
        * only those of its key variables (joinby's first merge, of the
        * using keys, takes no nolabel)
        quietly label dir
        local ulabels `r(names)'
        if "`nolabel'" != "" {
            local ulabels
            foreach v of local keys {
                local ulabels `ulabels' `: value label `v''
            }
            local ulabels : list uniq ulabels
        }
        local newlabels : list ulabels - masterlabels
        local commonvars `keys' `shared'
        mata: _cjoinby_meta_collect(tokens(st_local("newvars")), ///
            tokens(st_local("commonvars")), tokens(st_local("newlabels")))

        * strL values the plugin cannot write: kept in Mata, found by the
        * using observation number the plugin gives each output row
        local ustrl `newstrl'
        if "`update'" != "" local ustrl `ustrl' `sharedstrl'
        if "`ustrl'" != "" mata: _cjoinby_strl_hold(tokens(st_local("ustrl")))

        local payload `payshared' `paynew'
        local nkeys : word count `keys'
        local npay : word count `payload'
        if `timing' quietly _ctools_lap t_umeta t_begin `t_begin'
        mata: _cjoinby_strw(tokens(st_local("keys") + " " + st_local("payload")))
        * Set before the plugin holds the using data, so that cjoinby
        * releases them when a later step fails
        global CTOOLS_cjoinby_plugin "`__ctools_plugin'"
        plugin call ctools_plugin `keys' `payload', ///
            "cjoinby `threadopt' using `nkeys' `npay' `N0'"
        local using_N = _N
        if `timing' quietly _ctools_lap t_uplugin t_begin `t_begin'
    }
    local rc = _rc
    if `useframes' {
        frame change `masterframe'
        capture frame drop `uframe'
    }
    else if `rc' == 0 restore, preserve
    if `rc' exit `rc'
    if `timing' quietly _ctools_lap t_uframe t_begin `t_begin'

    * =========================================================================
    * Plan the join in the untouched master: an error here leaves nothing
    * to undo.
    * =========================================================================
    * joinby is a sortpreserve program: a sorted master gives a result in
    * the same sort order, with the using-only observations placed by their
    * values of the sort variables (their keys; missing elsewhere). The
    * plugin lays the output out that way when it can read the sort
    * variables; otherwise the joined data are sorted at the end.
    local keepusing = inlist("`unmatch'", "both", "using")
    local keepmaster = inlist("`unmatch'", "both", "master")
    local planextra
    local sortcols
    local nsort 0
    local needsort 0
    if "`sortedby'" != "" {
        local changed : list sortedby & shared
        if "`update'" != "" & "`changed'" != "" local needsort 1
        else if `keepusing' {
            local readable 1
            foreach v of local sortedby {
                if `: list v in keys' continue
                local mt : word `: list posof "`v'" in mastervars' of `mastertypes'
                if "`mt'" == "strL" | `: list v in aliasvars' local readable 0
                local planextra `planextra' `v'
            }
            local planextra : list uniq planextra
            if `readable' {
                local planvars `keys' `planextra'
                foreach v of local sortedby {
                    local sortcols `sortcols' `: list posof "`v'" in planvars'
                }
                local nsort : word count `sortcols'
            }
            else local planextra
        }
    }
    local planvars `keys' `planextra'
    local __cj_sortcols `sortcols'
    mata: _cjoinby_strw(tokens(st_local("planvars")))
    plugin call ctools_plugin `planvars', ///
        "cjoinby `threadopt' plan `nkeys' `=`keepmaster' + 2 * `keepusing'' `nsort'"
    local nout = scalar(_cjoinby_N)
    local n1 = scalar(_cjoinby_N_1)
    local n2 = scalar(_cjoinby_N_2)
    local n3 = scalar(_cjoinby_N_3)
    local stage_from = scalar(_cjoinby_stage_from)
    local placed = scalar(_cjoinby_placed)
    * using-only observations the plugin did not place by the sort variables
    if "`sortedby'" != "" & `n2' > 0 & !`placed' local needsort 1
    if `timing' quietly _ctools_lap t_plan t_begin `t_begin'

    * =========================================================================
    * Change the master: storage types, new variables and their metadata,
    * observations
    * =========================================================================
    if `undo' mata: _cjoinby_undo_record(tokens(st_local("commonvars")))
    * Validation hook: fail as Break would before the first change
    if "$CTOOLS_cjoinby_fault" == "record" exit 1

    * Key and shared variables take the wider of the two storage types
    mata: _cjoinby_promotions(tokens(st_local("commonvars")), ///
        tokens(st_local("keytargets") + " " + st_local("sharedtargets")))
    local i 0
    foreach v of local promote {
        local ++i
        local to : word `i' of `promoteto'
        quietly recast `to' `v'
        if substr("`to'", 1, 3) == "str" & "`to'" != "strL" {
            * merge widens a string's display format with its storage type
            local fmt : format `v'
            if regexm("`fmt'", "^%(-?)([0-9]+)s$") {
                if real(regexs(2)) < real(substr("`to'", 4, .)) {
                    format `v' %`=regexs(1)'`=substr("`to'", 4, .)'s
                }
            }
        }
    }
    * The master's variables the plugin moves, and its strL variables
    mata: _cjoinby_plugin_vars(tokens(st_local("mastervars")), tokens(st_local("aliasvars")))

    * The plugin carries the observation numbers of each output row's
    * master and using observation for the strL values copied afterwards.
    local moved = `stage_from' < `N0'
    local strlupdate
    if "`update'" != "" local strlupdate `sharedstrl'
    * strL master variables whose values must be rebuilt: all of them when
    * observations move; the ones update may change always (into copies,
    * so that the originals stay for a rollback)
    local mcopy
    if `moved' local mcopy `mstrl'
    else local mcopy `strlupdate'
    * e(sample), which cjoinby carries in a master variable, stays with one
    * pairing of each master observation in joinby (expand's copies are new
    * observations): the first one here.
    local esvar $CTOOLS_cjoinby_esvar
    local mobs
    local uobs
    if ("`mcopy'`esvar'" != "") & `moved' tempvar mobs
    if "`ustrl'" != "" tempvar uobs

    local addnames `mergevar' `newvars' `mobs' `uobs'
    local addtypes
    if "`mergevar'" != "" local addtypes byte
    local addtypes `addtypes' `newtypes'
    if "`mobs'" != "" local addtypes `addtypes' double
    if "`uobs'" != "" local addtypes `addtypes' double
    if "`addnames'" != "" {
        mata: _cjoinby_addvars(tokens(st_local("addnames")), tokens(st_local("addtypes")))
    }
    mata: _cjoinby_meta_apply(`="`update'" != ""')
    if `nout' > _N quietly set obs `nout'
    if `timing' quietly _ctools_lap t_schema t_begin `t_begin'
    * Validation hook: fail as Break would once the schema has changed
    if "$CTOOLS_cjoinby_fault" == "schema" exit 1

    * =========================================================================
    * Write the join
    * =========================================================================
    local nmaster : word count `mvars'
    local nnew : word count `paynew'
    mata: _cjoinby_positions(tokens(st_local("mvars")), tokens(st_local("keys")), ///
        tokens(st_local("payshared")), `nnew')
    local flags `update' `replace'
    if "`mergevar'" != "" local flags `flags' merge
    if "`mobs'" != "" local flags `flags' mobs
    if "`uobs'" != "" local flags `flags' uobs
    if `undo' local flags `flags' undo
    if "$CTOOLS_cjoinby_fault" == "store" local flags `flags' fault_store
    * cjoinby's rollback hands the plugin the variables it staged
    global CTOOLS_cjoinby_mvars `mvars'
    local wvars `mvars' `paynew' `mergevar' `mobs' `uobs'
    mata: _cjoinby_strw(tokens(st_local("wvars")))
    plugin call ctools_plugin `wvars', "cjoinby `threadopt' write `nmaster' `nnew' `flags'"
    if `timing' quietly _ctools_lap t_write t_begin `t_begin'
    * Validation hook: an interruption after the plugin has returned
    if "$CTOOLS_cjoinby_fault" == "post" exit 1

    if "`esvar'" != "" & "`mobs'" != "" & `nout' > 1 {
        quietly replace `esvar' = 0 if `mobs' == `mobs'[_n - 1] in 2/`nout'
    }

    * =========================================================================
    * strL values, by the observation numbers the plugin wrote
    * =========================================================================
    local copies
    if `nout' > 0 {
        foreach v of local mcopy {
            tempvar copy
            if `moved' quietly generate strL `copy' = `v'[`mobs'] in 1/`nout'
            else quietly generate strL `copy' = `v' in 1/`nout'
            local copies `copies' `copy'
        }
        local col 0
        foreach v of local ustrl {
            local ++col
            if `: list v in newstrl' {
                mata: _cjoinby_strl_fill("`v'", `col', "`uobs'", `nout')
            }
            else {
                local dest : word `: list posof "`v'" in mcopy' of `copies'
                mata: _cjoinby_strl_update("`dest'", `col', "`uobs'", ///
                    `="`replace'" != ""', `nout')
            }
        }
    }
    * Validation hook: fail as Break would once the strL values are in place
    if "$CTOOLS_cjoinby_fault" == "strl" exit 1

    * =========================================================================
    * The result stands from here on
    * =========================================================================
    nobreak {
        if _N > `nout' {
            if `nout' == 0 quietly drop in 1/l
            else quietly drop in `=`nout' + 1'/l
        }
        * The sort order of the master: the output is in it already unless
        * the plugin could not place the using-only observations or update
        * changed a sort variable. Sorted here, while a failure can still be
        * undone, unless a sort variable is one of the strL variables about
        * to be rebuilt.
        local sortlate : list sortedby & mcopy
        if `needsort' & "`sortlate'" == "" & `nout' > 0 {
            sort `sortedby', stable
            local needsort 0
        }
        * Each copy takes its variable's place: format, label and
        * characteristics (notes included)
        if "`mcopy'" != "" & `nout' > 0 {
            local left `copies'
            foreach v of local mcopy {
                gettoken copy left : left
                format `copy' `: format `v''
                mata: st_varlabel(st_local("copy"), st_varlabel(st_local("v")))
                local chars : char `v'[]
                foreach c of local chars {
                    mata: st_global(st_local("copy") + "[" + st_local("c") + "]", ///
                        st_global(st_local("v") + "[" + st_local("c") + "]"))
                }
                drop `v'
                rename `copy' `v'
            }
        }
        if "`mobs'`uobs'" != "" drop `mobs' `uobs'
        if "`mcopy'" != "" & `nout' > 0 order `mastervars' `mergevar' `newvars'

        if `nout' == 0 di as text "(joinby on `keys' is empty)"
        else _cjoinby_mergelabel "`mergevar'" "`update'"

        * A sort of data already in order only sets the marker again
        * (recast and the strL copies drop it).
        if "`sortedby'" != "" & `nout' > 0 {
            if `needsort' | "`: sortedby'" != "`sortedby'" sort `sortedby', stable
        }
        * joinby marks the data changed; plugin stores do not
        mata: st_updata(1)
    }

    * (without frames: the result stands, not the preserved master)
    if !`useframes' restore, not

    local elapsed = .
    if `timing' {
        quietly _ctools_lap t_post t_begin `t_begin'
        local elapsed = `t_use' + `t_umeta' + `t_uplugin' + `t_uframe' + `t_plan' + ///
            `t_schema' + `t_write' + `t_post'
        di as text ""
        di as text "{hline 55}"
        di as text "cjoinby timing breakdown:"
        di as text "{hline 55}"
        di as text "  Key variables:            " as result "`keys'"
        di as text "  Master / using obs:       " as result `N0' as text " / " as result `using_N'
        di as text "  Key groups in using:      " as result scalar(_cjoinby_groups) ///
            as text cond(scalar(_cjoinby_direct), " (direct lookup)", " (hashed lookup)")
        di as text "  Output obs:               " as result `nout'
        di as text "  {hline 53}"
        di as text "  C plugin:"
        di as text "    Load using data:        " as result %8.4f scalar(_cjoinby_t_using_load) " sec"
        di as text "    Index using keys:       " as result %8.4f scalar(_cjoinby_t_using_index) " sec"
        di as text "    Load master keys:       " as result %8.4f scalar(_cjoinby_t_master_load) " sec"
        di as text "    Match and lay out:      " as result %8.4f scalar(_cjoinby_t_match) " sec"
        di as text "    Stage master rows:      " as result %8.4f scalar(_cjoinby_t_stage) " sec"
        di as text "    Write joined data:      " as result %8.4f scalar(_cjoinby_t_write) " sec"
        local t_c = scalar(_cjoinby_t_using_load) + scalar(_cjoinby_t_using_index) + ///
            scalar(_cjoinby_t_master_load) + scalar(_cjoinby_t_match) + ///
            scalar(_cjoinby_t_stage) + scalar(_cjoinby_t_write)
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f `t_c' " sec"
        di as text "  {hline 53}"
        di as text "  Stata:"
        di as text "    Read using file:        " as result %8.4f `t_use' " sec"
        di as text "    Using metadata:         " as result %8.4f (`t_umeta' + `t_uframe') " sec"
        di as text "    New variables, obs:     " as result %8.4f `t_schema' " sec"
        di as text "    Finish (strL, sort):    " as result %8.4f `t_post' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f ///
            (`t_uplugin' + `t_plan' + `t_write' - `t_c') " sec"
        di as text "  {hline 53}"
        di as text "    Wall clock total:       " as result %8.4f `elapsed' " sec"
        di as text "{hline 55}"
        di as text "  Threads:                  " as result scalar(_cjoinby_threads_max) ///
            as text cond(scalar(_cjoinby_openmp_enabled), "", " (OpenMP disabled)")
        di as text "  Master obs staged:        " as result scalar(_cjoinby_staged)
        di as text "{hline 55}"
    }

    return scalar N = _N
    return scalar N_1 = `n1'
    return scalar N_2 = `n2'
    return scalar N_3 = `n3'
    return scalar time = `elapsed'
    return local using `"`using'"'
    return local keyvars "`keys'"
end

* unmatched(): none (the default), both, master or using, abbreviated as
* joinby allows. Sets local unmatch in the caller.
program define _cjoinby_unmatched
    args value
    if `"`value'"' == "" {
        c_local unmatch none
        exit
    }
    local value = lower(`"`value'"')
    local len = length(`"`value'"')
    local found
    foreach word in both master using none {
        if `"`value'"' == substr("`word'", 1, `len') {
            local found `word'
            continue, break
        }
    }
    if "`found'" == "" {
        di as error `"`value' is invalid value for unmatch()"'
        exit 198
    }
    c_local unmatch `found'
end

* joinby's value label for _merge. Like joinby, it is defined even when no
* _merge variable is kept, replaces a __MERGE that joinby (or cjoinby)
* defined, and leaves alone a __MERGE of the user's.
program define _cjoinby_mergelabel
    args mergevar update
    capture label list __MERGE
    if !_rc {
        local signature : label __MERGE 9999
        if "`signature'" == "signature" label drop __MERGE
        else {
            if "`mergevar'" != "" di as text "value label __MERGE for `mergevar' already exists"
            exit
        }
    }
    if "`update'" == "" {
        label define __MERGE 1 "only in master data" 2 "only in using data" ///
            3 "both in master and using data" 9999 "signature"
    }
    else {
        label define __MERGE 1 "only in master data" 2 "only in using data" ///
            3 "in both, master agrees with using data" ///
            4 "in both, missing in master data updated" ///
            5 "in both, master disagrees with using data" 9999 "signature"
    }
    if "`mergevar'" != "" label values `mergevar' __MERGE
end
